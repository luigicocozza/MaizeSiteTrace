#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import random
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from sklearn.decomposition import PCA
from torch_geometric.nn import GAE, GCNConv
from torch_geometric.utils import negative_sampling


BASE_DIR = Path(".").resolve()

PYG_V3 = BASE_DIR / "data/processed/pyg/maize_ncrna_graph_v3_prjna419326.pt"
GRAPH_NODES_V3 = BASE_DIR / "data/processed/graph_nodes_v3_prjna419326.csv"

OUT_DIR = BASE_DIR / "outputs/v3_gae"
MODEL_DIR = BASE_DIR / "models/graph_autoencoder_v3"
EMB_DIR = BASE_DIR / "data/processed/embeddings"

OUT_DIR.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)
EMB_DIR.mkdir(parents=True, exist_ok=True)

MODEL_PATH = MODEL_DIR / "gae_v3_model.pt"
METRICS_JSON = OUT_DIR / "v3_gae_metrics.json"
SUMMARY_TXT = OUT_DIR / "v3_gae_summary.txt"
HISTORY_CSV = OUT_DIR / "v3_gae_training_history.csv"

EMB_ALL_CSV = EMB_DIR / "gae_v3_node_embeddings.csv"
EMB_EXPERIMENTAL_CSV = EMB_DIR / "gae_v3_experimental_srna_embeddings.csv"
EMB_NCRNA_CSV = EMB_DIR / "gae_v3_rnacentral_ncrna_embeddings.csv"
EMB_PCA_CSV = EMB_DIR / "gae_v3_node_embeddings_pca2.csv"

SEED = 42
HIDDEN_CHANNELS = 64
LATENT_CHANNELS = 16
EPOCHS = 300
LR = 0.01
WEIGHT_DECAY = 1e-5
DROPOUT = 0.20
VAL_RATIO = 0.10
TEST_RATIO = 0.10


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def load_pt(path: Path):
    try:
        return torch.load(path, weights_only=False)
    except TypeError:
        return torch.load(path)


class GCNEncoder(torch.nn.Module):
    def __init__(self, in_channels: int, hidden_channels: int, out_channels: int, dropout: float):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, out_channels)
        self.dropout = dropout

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, p=self.dropout, training=self.training)
        x = self.conv2(x, edge_index)
        return x


def split_edges(edge_index: torch.Tensor, val_ratio: float, test_ratio: float):
    num_edges = edge_index.shape[1]
    perm = torch.randperm(num_edges)

    n_test = int(num_edges * test_ratio)
    n_val = int(num_edges * val_ratio)
    n_train = num_edges - n_val - n_test

    train_idx = perm[:n_train]
    val_idx = perm[n_train:n_train + n_val]
    test_idx = perm[n_train + n_val:]

    return (
        edge_index[:, train_idx],
        edge_index[:, val_idx],
        edge_index[:, test_idx],
    )


def standardize_x(x: torch.Tensor):
    x = x.float()
    x = torch.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

    mean = x.mean(dim=0, keepdim=True)
    std = x.std(dim=0, keepdim=True)
    std = torch.where(std == 0, torch.ones_like(std), std)

    x_scaled = (x - mean) / std
    x_scaled = torch.nan_to_num(x_scaled, nan=0.0, posinf=0.0, neginf=0.0)

    return x_scaled, mean, std


def make_negative_edges(full_edge_index: torch.Tensor, num_nodes: int, num_neg: int):
    return negative_sampling(
        edge_index=full_edge_index,
        num_nodes=num_nodes,
        num_neg_samples=num_neg,
        method="sparse",
    )


def safe_test(model: GAE, z: torch.Tensor, pos_edge_index: torch.Tensor, full_edge_index: torch.Tensor, num_nodes: int):
    neg_edge_index = make_negative_edges(full_edge_index, num_nodes, pos_edge_index.shape[1])
    auc, ap = model.test(z, pos_edge_index, neg_edge_index)
    return float(auc), float(ap)


def main() -> int:
    set_seed(SEED)

    print("TRAINING GAE V3 - PRJNA419326")
    print("=" * 80)
    print(f"Input graph: {PYG_V3}")

    obj = load_pt(PYG_V3)
    data = obj["data"]
    metadata = obj["metadata"]

    num_nodes = int(data.num_nodes)
    num_edges = int(data.edge_index.shape[1])
    in_channels = int(data.x.shape[1])

    print(f"Nodi: {num_nodes}")
    print(f"Archi: {num_edges}")
    print(f"Feature input: {in_channels}")

    x, x_mean, x_std = standardize_x(data.x)
    edge_index = data.edge_index.long()

    train_pos, val_pos, test_pos = split_edges(edge_index, VAL_RATIO, TEST_RATIO)

    print("")
    print("Split archi:")
    print(f"Train edges: {train_pos.shape[1]}")
    print(f"Val edges:   {val_pos.shape[1]}")
    print(f"Test edges:  {test_pos.shape[1]}")
    print("")

    model = GAE(
        GCNEncoder(
            in_channels=in_channels,
            hidden_channels=HIDDEN_CHANNELS,
            out_channels=LATENT_CHANNELS,
            dropout=DROPOUT,
        )
    )

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=LR,
        weight_decay=WEIGHT_DECAY,
    )

    history = []

    for epoch in range(1, EPOCHS + 1):
        model.train()
        optimizer.zero_grad()

        z = model.encode(x, train_pos)
        loss = model.recon_loss(z, train_pos)

        loss.backward()
        optimizer.step()

        row = {
            "epoch": epoch,
            "train_loss": float(loss.item()),
            "val_auc": "",
            "val_ap": "",
        }

        if epoch == 1 or epoch % 10 == 0 or epoch == EPOCHS:
            model.eval()
            with torch.no_grad():
                z_eval = model.encode(x, train_pos)
                val_auc, val_ap = safe_test(model, z_eval, val_pos, edge_index, num_nodes)

            row["val_auc"] = val_auc
            row["val_ap"] = val_ap

            print(
                f"Epoch {epoch:03d} | "
                f"loss={loss.item():.6f} | "
                f"val_auc={val_auc:.4f} | "
                f"val_ap={val_ap:.4f}"
            )

        history.append(row)

    model.eval()
    with torch.no_grad():
        z = model.encode(x, train_pos)
        test_auc, test_ap = safe_test(model, z, test_pos, edge_index, num_nodes)

    print("")
    print("Test finale:")
    print(f"test_auc={test_auc:.4f}")
    print(f"test_ap={test_ap:.4f}")

    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "model_type": "GAE_GCNEncoder_V3",
            "input_graph": str(PYG_V3),
            "in_channels": in_channels,
            "hidden_channels": HIDDEN_CHANNELS,
            "latent_channels": LATENT_CHANNELS,
            "dropout": DROPOUT,
            "x_mean": x_mean,
            "x_std": x_std,
            "train_pos_edge_index": train_pos,
            "val_pos_edge_index": val_pos,
            "test_pos_edge_index": test_pos,
            "seed": SEED,
            "created": datetime.now().isoformat(timespec="seconds"),
        },
        MODEL_PATH,
    )

    with HISTORY_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["epoch", "train_loss", "val_auc", "val_ap"])
        writer.writeheader()
        writer.writerows(history)

    z_np = z.detach().cpu().numpy()

    idx_to_node_id = metadata["mappings"]["idx_to_node_id"]

    nodes_df = pd.read_csv(GRAPH_NODES_V3)
    node_meta = nodes_df.set_index("node_id").to_dict(orient="index")

    rows = []

    for idx in range(num_nodes):
        node_id = idx_to_node_id.get(str(idx), f"unknown:{idx}")
        meta = node_meta.get(node_id, {})

        row = {
            "node_idx": idx,
            "node_id": node_id,
            "node_type": meta.get("node_type", ""),
            "name": meta.get("name", ""),
            "source": meta.get("source", ""),
            "rna_type": meta.get("rna_type", ""),
            "sequence_length": meta.get("sequence_length", ""),
            "expression_status": meta.get("expression_status", ""),
            "is_experimental_srna": str(node_id).startswith("experimental_srna:"),
            "is_rnacentral_ncrna": str(node_id).startswith("ncrna:"),
        }

        for j in range(z_np.shape[1]):
            row[f"z_{j:02d}"] = float(z_np[idx, j])

        rows.append(row)

    emb_df = pd.DataFrame(rows)
    emb_df.to_csv(EMB_ALL_CSV, index=False)

    emb_df[emb_df["is_experimental_srna"] == True].to_csv(EMB_EXPERIMENTAL_CSV, index=False)
    emb_df[emb_df["is_rnacentral_ncrna"] == True].to_csv(EMB_NCRNA_CSV, index=False)

    pca = PCA(n_components=2, random_state=SEED)
    pca2 = pca.fit_transform(z_np)

    pca_df = emb_df[
        [
            "node_idx",
            "node_id",
            "node_type",
            "name",
            "source",
            "rna_type",
            "expression_status",
            "is_experimental_srna",
            "is_rnacentral_ncrna",
        ]
    ].copy()

    pca_df["pca1"] = pca2[:, 0]
    pca_df["pca2"] = pca2[:, 1]
    pca_df.to_csv(EMB_PCA_CSV, index=False)

    metrics = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "input_graph": str(PYG_V3),
        "model_path": str(MODEL_PATH),
        "n_nodes": num_nodes,
        "n_edges": num_edges,
        "n_node_features": in_channels,
        "train_edges": int(train_pos.shape[1]),
        "val_edges": int(val_pos.shape[1]),
        "test_edges": int(test_pos.shape[1]),
        "epochs": EPOCHS,
        "learning_rate": LR,
        "weight_decay": WEIGHT_DECAY,
        "hidden_channels": HIDDEN_CHANNELS,
        "latent_channels": LATENT_CHANNELS,
        "dropout": DROPOUT,
        "final_train_loss": float(history[-1]["train_loss"]),
        "test_auc": test_auc,
        "test_ap": test_ap,
        "embedding_shape": list(z_np.shape),
        "pca_explained_variance_ratio": [float(x) for x in pca.explained_variance_ratio_],
        "outputs": {
            "model": str(MODEL_PATH),
            "metrics_json": str(METRICS_JSON),
            "summary_txt": str(SUMMARY_TXT),
            "history_csv": str(HISTORY_CSV),
            "embeddings_all": str(EMB_ALL_CSV),
            "embeddings_experimental_srna": str(EMB_EXPERIMENTAL_CSV),
            "embeddings_rnacentral_ncrna": str(EMB_NCRNA_CSV),
            "embeddings_pca2": str(EMB_PCA_CSV),
        },
    }

    METRICS_JSON.write_text(json.dumps(metrics, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("GAE V3 TRAINING COMPLETATO")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Grafo: {PYG_V3}")
    lines.append(f"Nodi: {num_nodes}")
    lines.append(f"Archi: {num_edges}")
    lines.append(f"Feature input: {in_channels}")
    lines.append(f"Embedding dimension: {LATENT_CHANNELS}")
    lines.append("")
    lines.append("Split archi:")
    lines.append(f"Train: {train_pos.shape[1]}")
    lines.append(f"Validation: {val_pos.shape[1]}")
    lines.append(f"Test: {test_pos.shape[1]}")
    lines.append("")
    lines.append("Metriche finali:")
    lines.append(f"Final train loss: {history[-1]['train_loss']:.6f}")
    lines.append(f"Test AUC: {test_auc:.6f}")
    lines.append(f"Test AP: {test_ap:.6f}")
    lines.append("")
    lines.append("Output:")
    lines.append(str(MODEL_PATH))
    lines.append(str(METRICS_JSON))
    lines.append(str(HISTORY_CSV))
    lines.append(str(EMB_ALL_CSV))
    lines.append(str(EMB_EXPERIMENTAL_CSV))
    lines.append(str(EMB_NCRNA_CSV))
    lines.append(str(EMB_PCA_CSV))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
