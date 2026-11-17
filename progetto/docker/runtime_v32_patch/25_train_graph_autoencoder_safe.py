#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import random
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from sklearn.decomposition import PCA
from sklearn.metrics import average_precision_score, roc_auc_score

from torch_geometric.nn import GCNConv
from torch_geometric.utils import negative_sampling, remove_self_loops, to_undirected


BASE_DIR = Path(".").resolve()

PYG_PATH = BASE_DIR / "data/processed/pyg/maize_ncrna_graph.pt"
GRAPH_NODES_CSV = BASE_DIR / "data/processed/graph_nodes.csv"

from datetime import datetime as _primetta_dt_safe
from uuid import uuid4 as _primetta_uuid_safe
__primetta_gae_safe_root = (BASE_DIR / 'outputs' / 'gae_gui_safe_runs' / (_primetta_dt_safe.now().strftime('%Y%m%dT%H%M%S%f') + '_' + _primetta_uuid_safe().hex[:8]))
OUT_DIR = __primetta_gae_safe_root / 'results'
PLOTS_DIR = OUT_DIR / "plots"
MODEL_DIR = __primetta_gae_safe_root / 'models'
EMB_DIR = __primetta_gae_safe_root / 'embeddings'

MODEL_PATH = MODEL_DIR / "gae_model.pt"
HISTORY_CSV = OUT_DIR / "gae_training_history.csv"
METRICS_JSON = OUT_DIR / "gae_metrics.json"
SUMMARY_TXT = OUT_DIR / "gae_summary.txt"

NODE_EMB_CSV = EMB_DIR / "gae_node_embeddings.csv"
NCRNA_EMB_CSV = EMB_DIR / "gae_ncrna_embeddings.csv"


class GCNEncoder(nn.Module):
    def __init__(self, in_channels: int, hidden_channels: int, latent_channels: int, dropout: float):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, latent_channels)
        self.dropout = dropout

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor) -> torch.Tensor:
        h = self.conv1(x, edge_index)
        h = F.relu(h)
        h = F.dropout(h, p=self.dropout, training=self.training)
        z = self.conv2(h, edge_index)
        return z


class GraphAutoEncoder(nn.Module):
    def __init__(self, in_channels: int, hidden_channels: int, latent_channels: int, dropout: float):
        super().__init__()
        self.encoder = GCNEncoder(in_channels, hidden_channels, latent_channels, dropout)

    def encode(self, x: torch.Tensor, edge_index: torch.Tensor) -> torch.Tensor:
        return self.encoder(x, edge_index)

    @staticmethod
    def decode(z: torch.Tensor, edge_index: torch.Tensor) -> torch.Tensor:
        src = edge_index[0]
        dst = edge_index[1]
        return (z[src] * z[dst]).sum(dim=1)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def load_pyg_data():
    if not PYG_PATH.exists():
        raise FileNotFoundError(f"Dataset PyG non trovato: {PYG_PATH}")

    try:
        obj = torch.load(PYG_PATH, map_location="cpu", weights_only=False)
    except TypeError:
        obj = torch.load(PYG_PATH, map_location="cpu")

    if isinstance(obj, dict) and "data" in obj:
        return obj["data"], obj.get("metadata", {})

    return obj, {}


def unique_undirected_edges(edge_index: torch.Tensor) -> torch.Tensor:
    edge_index, _ = remove_self_loops(edge_index)
    edge_index = to_undirected(edge_index)

    pairs = set()
    for u, v in edge_index.t().tolist():
        if int(u) == int(v):
            continue
        a, b = sorted((int(u), int(v)))
        pairs.add((a, b))

    pairs = sorted(pairs)
    return torch.tensor(pairs, dtype=torch.long).t().contiguous()


def split_edges(pos_edges: torch.Tensor, seed: int):
    n = pos_edges.size(1)
    idx = list(range(n))
    random.Random(seed).shuffle(idx)

    n_train = int(n * 0.80)
    n_val = int(n * 0.10)

    train_idx = idx[:n_train]
    val_idx = idx[n_train:n_train + n_val]
    test_idx = idx[n_train + n_val:]

    return pos_edges[:, train_idx], pos_edges[:, val_idx], pos_edges[:, test_idx]


def make_message_passing_edges(train_pos: torch.Tensor) -> torch.Tensor:
    return torch.cat([train_pos, train_pos.flip(0)], dim=1).contiguous()


def sample_negatives(all_pos: torch.Tensor, num_nodes: int, n_samples: int) -> torch.Tensor:
    positive_directed = torch.cat([all_pos, all_pos.flip(0)], dim=1).contiguous()
    return negative_sampling(
        edge_index=positive_directed,
        num_nodes=num_nodes,
        num_neg_samples=n_samples,
        method="sparse",
    )


def loss_fn(model: GraphAutoEncoder, z: torch.Tensor, pos_edges: torch.Tensor, neg_edges: torch.Tensor) -> torch.Tensor:
    pos_logits = model.decode(z, pos_edges)
    neg_logits = model.decode(z, neg_edges)

    logits = torch.cat([pos_logits, neg_logits])
    labels = torch.cat([
        torch.ones(pos_logits.numel()),
        torch.zeros(neg_logits.numel()),
    ])

    return F.binary_cross_entropy_with_logits(logits, labels)


@torch.no_grad()
def evaluate(model: GraphAutoEncoder, z: torch.Tensor, pos_edges: torch.Tensor, neg_edges: torch.Tensor) -> dict:
    pos_logits = model.decode(z, pos_edges)
    neg_logits = model.decode(z, neg_edges)

    logits = torch.cat([pos_logits, neg_logits])
    probs = torch.sigmoid(logits).cpu().numpy()

    labels = np.concatenate([
        np.ones(pos_logits.numel(), dtype=np.int64),
        np.zeros(neg_logits.numel(), dtype=np.int64),
    ])

    return {
        "auc": float(roc_auc_score(labels, probs)),
        "average_precision": float(average_precision_score(labels, probs)),
    }


def save_embeddings(z: torch.Tensor) -> None:
    EMB_DIR.mkdir(parents=True, exist_ok=True)

    z_np = z.detach().cpu().numpy()
    emb_cols = [f"gae_emb_{i+1:03d}" for i in range(z_np.shape[1])]

    if GRAPH_NODES_CSV.exists():
        nodes = pd.read_csv(GRAPH_NODES_CSV)
    else:
        nodes = pd.DataFrame({"node_index": range(z_np.shape[0])})

    if len(nodes) != z_np.shape[0]:
        nodes = pd.DataFrame({"node_index": range(z_np.shape[0])})

    emb = pd.DataFrame(z_np, columns=emb_cols)
    out = pd.concat([nodes.reset_index(drop=True), emb], axis=1)
    out.to_csv(NODE_EMB_CSV, index=False)

    if "node_type" in out.columns:
        ncrna = out[out["node_type"].astype(str).str.lower().eq("ncrna")].copy()
    else:
        ncrna = out.copy()

    ncrna.to_csv(NCRNA_EMB_CSV, index=False)


def save_plots(history: list[dict], z: torch.Tensor) -> None:
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    hist = pd.DataFrame(history)

    plt.figure(figsize=(8, 5))
    plt.plot(hist["epoch"], hist["loss"])
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Graph Autoencoder training loss")
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "gae_loss_curve.png", dpi=160)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(hist["epoch"], hist["val_average_precision"])
    plt.xlabel("Epoch")
    plt.ylabel("Validation average precision")
    plt.title("Graph Autoencoder validation AP")
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "gae_validation_ap_curve.png", dpi=160)
    plt.close()

    z_np = z.detach().cpu().numpy()
    coords = PCA(n_components=2, random_state=42).fit_transform(z_np)

    plt.figure(figsize=(7, 6))
    plt.scatter(coords[:, 0], coords[:, 1], s=18)
    plt.xlabel("GAE component 1")
    plt.ylabel("GAE component 2")
    plt.title("GAE node embeddings")
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "gae_node_embeddings_2d.png", dpi=160)
    plt.close()

    if GRAPH_NODES_CSV.exists():
        nodes = pd.read_csv(GRAPH_NODES_CSV)
        if len(nodes) == z_np.shape[0] and "node_type" in nodes.columns:
            mask = nodes["node_type"].astype(str).str.lower().eq("ncrna").values
            if mask.sum() > 0:
                coords_n = coords[mask]
                plt.figure(figsize=(7, 6))
                plt.scatter(coords_n[:, 0], coords_n[:, 1], s=18)
                plt.xlabel("GAE component 1")
                plt.ylabel("GAE component 2")
                plt.title("GAE ncRNA embeddings")
                plt.tight_layout()
                plt.savefig(PLOTS_DIR / "gae_ncrna_embeddings_2d.png", dpi=160)
                plt.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--hidden-dim", type=int, default=64)
    parser.add_argument("--latent-dim", type=int, default=16)
    parser.add_argument("--lr", type=float, default=0.01)
    parser.add_argument("--dropout", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    set_seed(args.seed)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    EMB_DIR.mkdir(parents=True, exist_ok=True)

    data, metadata = load_pyg_data()

    x = data.x.float()
    edge_index = data.edge_index.long()

    num_nodes = x.size(0)
    in_channels = x.size(1)

    all_pos = unique_undirected_edges(edge_index)
    train_pos, val_pos, test_pos = split_edges(all_pos, args.seed)
    train_edge_index = make_message_passing_edges(train_pos)

    val_neg = sample_negatives(all_pos, num_nodes, val_pos.size(1))
    test_neg = sample_negatives(all_pos, num_nodes, test_pos.size(1))

    model = GraphAutoEncoder(
        in_channels=in_channels,
        hidden_channels=args.hidden_dim,
        latent_channels=args.latent_dim,
        dropout=args.dropout,
    )

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=1e-4)

    history = []
    best_state = None
    best_val_ap = -1.0

    for epoch in range(1, args.epochs + 1):
        model.train()
        optimizer.zero_grad()

        z = model.encode(x, train_edge_index)
        train_neg = sample_negatives(all_pos, num_nodes, train_pos.size(1))
        loss = loss_fn(model, z, train_pos, train_neg)

        loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            z_eval = model.encode(x, train_edge_index)
            val_metrics = evaluate(model, z_eval, val_pos, val_neg)

        val_ap = val_metrics["average_precision"]
        if val_ap > best_val_ap:
            best_val_ap = val_ap
            best_state = {
                "epoch": epoch,
                "model_state_dict": model.state_dict(),
                "val_metrics": val_metrics,
            }

        row = {
            "epoch": epoch,
            "loss": float(loss.item()),
            "val_auc": val_metrics["auc"],
            "val_average_precision": val_metrics["average_precision"],
        }
        history.append(row)

        if epoch == 1 or epoch % 25 == 0 or epoch == args.epochs:
            print(
                f"Epoch {epoch:04d} | "
                f"loss={row['loss']:.5f} | "
                f"val_auc={row['val_auc']:.4f} | "
                f"val_ap={row['val_average_precision']:.4f}"
            )

    if best_state is not None:
        model.load_state_dict(best_state["model_state_dict"])

    model.eval()
    with torch.no_grad():
        z_final = model.encode(x, train_edge_index)
        train_neg_final = sample_negatives(all_pos, num_nodes, train_pos.size(1))
        train_metrics = evaluate(model, z_final, train_pos, train_neg_final)
        val_metrics = evaluate(model, z_final, val_pos, val_neg)
        test_metrics = evaluate(model, z_final, test_pos, test_neg)

    pd.DataFrame(history).to_csv(HISTORY_CSV, index=False)
    save_embeddings(z_final)
    save_plots(history, z_final)

    torch.save(
        {
            "created": datetime.now().isoformat(timespec="seconds"),
            "model_state_dict": model.state_dict(),
            "args": vars(args),
            "in_channels": in_channels,
            "train_edge_index": train_edge_index,
        },
        MODEL_PATH,
    )

    metrics = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "module": "Graph Autoencoder experimental neural network",
        "training_type": "unsupervised/self-supervised link reconstruction",
        "important_note": "This is not a supervised drought-response predictor.",
        "input": {
            "pyg_dataset": str(PYG_PATH),
            "nodes": int(num_nodes),
            "node_features": int(in_channels),
            "positive_undirected_edges": int(all_pos.size(1)),
            "train_edges": int(train_pos.size(1)),
            "validation_edges": int(val_pos.size(1)),
            "test_edges": int(test_pos.size(1)),
        },
        "model": {
            "encoder": "GCNConv + GCNConv",
            "decoder": "inner product",
            "hidden_dim": args.hidden_dim,
            "latent_dim": args.latent_dim,
            "epochs": args.epochs,
            "lr": args.lr,
            "dropout": args.dropout,
        },
        "metrics": {
            "train": train_metrics,
            "validation": val_metrics,
            "test": test_metrics,
        },
        "outputs": {
            "model": str(MODEL_PATH),
            "node_embeddings": str(NODE_EMB_CSV),
            "ncrna_embeddings": str(NCRNA_EMB_CSV),
            "history": str(HISTORY_CSV),
            "summary": str(SUMMARY_TXT),
            "plots": str(PLOTS_DIR),
        },
    }

    METRICS_JSON.write_text(json.dumps(metrics, indent=2, ensure_ascii=False), encoding="utf-8")

    summary = [
        "GRAPH AUTOENCODER EXPERIMENTAL MODULE",
        "=" * 80,
        "",
        "Status: completed",
        "Training type: unsupervised/self-supervised link reconstruction",
        "Architecture: GCN encoder + inner-product decoder",
        "",
        f"Nodes: {num_nodes}",
        f"Input node features: {in_channels}",
        f"Positive undirected edges: {all_pos.size(1)}",
        f"Train edges: {train_pos.size(1)}",
        f"Validation edges: {val_pos.size(1)}",
        f"Test edges: {test_pos.size(1)}",
        "",
        f"Hidden dim: {args.hidden_dim}",
        f"Latent dim: {args.latent_dim}",
        f"Epochs: {args.epochs}",
        f"Learning rate: {args.lr}",
        f"Dropout: {args.dropout}",
        "",
        "Metrics:",
        f"Train AUC: {train_metrics['auc']}",
        f"Train AP: {train_metrics['average_precision']}",
        f"Validation AUC: {val_metrics['auc']}",
        f"Validation AP: {val_metrics['average_precision']}",
        f"Test AUC: {test_metrics['auc']}",
        f"Test AP: {test_metrics['average_precision']}",
        "",
        "Outputs:",
        f"Model: {MODEL_PATH}",
        f"Node embeddings: {NODE_EMB_CSV}",
        f"ncRNA embeddings: {NCRNA_EMB_CSV}",
        f"History: {HISTORY_CSV}",
        f"Metrics JSON: {METRICS_JSON}",
        f"Plots directory: {PLOTS_DIR}",
        "",
        "Methodological note:",
        "This is a real neural network module on the graph.",
        "It is exploratory and unsupervised.",
        "It does not prove drought-response biology.",
    ]

    SUMMARY_TXT.write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
