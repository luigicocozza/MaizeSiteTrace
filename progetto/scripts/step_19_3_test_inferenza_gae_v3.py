#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 19.3
TEST REALE DI INFERENZA GRAPH AUTOENCODER V3

Utilizza checkpoint e grafo originali.
Esegue soltanto inferenza CPU, senza training.

ATTENZIONE:
Il grafo PyG viene deserializzato con weights_only=False,
come previsto dal codice storico del progetto.
Eseguire soltanto su file locali ritenuti attendibili.
"""

from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import subprocess
import sys


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

RUNS = ROOT / "target_prediction_v33_vscode" / "runs"

IMAGE = "ai_primetta_v33_gnn_cpu:18.6"

MODEL = (
    ROOT / "models/graph_autoencoder_v3/gae_v3_model.pt"
)

GRAPH = (
    ROOT
    / "data/processed/pyg/"
    "maize_ncrna_graph_v3_prjna419326.pt"
)


RUNNER = r'''
import json
import math
import sys

import torch
import torch.nn.functional as F
from torch_geometric.nn import GAE, GCNConv


MODEL = "/project/models/graph_autoencoder_v3/gae_v3_model.pt"
GRAPH = (
    "/project/data/processed/pyg/"
    "maize_ncrna_graph_v3_prjna419326.pt"
)


class GCNEncoder(torch.nn.Module):
    """
    Architettura coerente con
    51_train_v3_graph_autoencoder.py.
    """

    def __init__(
        self,
        in_channels,
        hidden_channels,
        out_channels,
        dropout
    ):
        super().__init__()

        self.conv1 = GCNConv(
            in_channels,
            hidden_channels
        )

        self.conv2 = GCNConv(
            hidden_channels,
            out_channels
        )

        self.dropout = dropout

    def forward(self, x, edge_index):

        x = self.conv1(x, edge_index)
        x = F.relu(x)

        x = F.dropout(
            x,
            p=self.dropout,
            training=self.training
        )

        x = self.conv2(x, edge_index)
        return x


def main():

    torch.set_num_threads(2)

    print("[1] CARICAMENTO CHECKPOINT")

    # Caricamento ristretto al formato sicuro
    # per tensori e strutture elementari.
    checkpoint = torch.load(
        MODEL,
        map_location="cpu",
        weights_only=True
    )

    if not isinstance(checkpoint, dict):
        raise RuntimeError(
            "Checkpoint non in formato dizionario"
        )

    required = [
        "model_state_dict",
        "in_channels",
        "hidden_channels",
        "latent_channels",
    ]

    for key in required:
        if key not in checkpoint:
            raise RuntimeError(
                "Campo checkpoint assente: " + key
            )

    print(
        "[OK] Tipo:",
        checkpoint.get("model_type", "NON DICHIARATO")
    )

    print(
        "[OK] Parametri:",
        {
            key: checkpoint[key]
            for key in required
            if key != "model_state_dict"
        }
    )

    print("[2] CARICAMENTO GRAFO REALE")

    # File locale storico PyG:
    # deserializzazione completa richiesta
    # per oggetti Data di torch_geometric.
    graph_object = torch.load(
        GRAPH,
        map_location="cpu",
        weights_only=False
    )

    if isinstance(graph_object, dict):
        if "data" not in graph_object:
            raise RuntimeError(
                "Campo data non presente nel grafo"
            )
        data = graph_object["data"]
    else:
        data = graph_object

    x = data.x.float()
    edges = data.edge_index.long()

    if x.ndim != 2:
        raise RuntimeError(
            "Matrice delle feature non valida"
        )

    if edges.ndim != 2 or edges.shape[0] != 2:
        raise RuntimeError(
            "Matrice edge_index non valida"
        )

    n_nodes = int(x.shape[0])
    n_features = int(x.shape[1])
    n_edges = int(edges.shape[1])

    if n_nodes == 0 or n_edges == 0:
        raise RuntimeError(
            "Grafo privo di nodi o archi"
        )

    if edges.min().item() < 0:
        raise RuntimeError(
            "Indice arco negativo"
        )

    if edges.max().item() >= n_nodes:
        raise RuntimeError(
            "Archi riferiti a nodi inesistenti"
        )

    print("[OK] Nodi:", n_nodes)
    print("[OK] Archi:", n_edges)
    print("[OK] Feature:", n_features)

    if n_features != int(checkpoint["in_channels"]):
        raise RuntimeError(
            "Dimensioni feature incompatibili"
        )

    print("[3] NORMALIZZAZIONE FEATURE")

    x = torch.nan_to_num(
        x,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    saved_mean = checkpoint.get("x_mean")
    saved_std = checkpoint.get("x_std")

    if (
        isinstance(saved_mean, torch.Tensor)
        and isinstance(saved_std, torch.Tensor)
    ):
        mean = saved_mean.float()
        std = saved_std.float()

        normalization = "PARAMETRI_CHECKPOINT"

    else:
        # Replica la formula di standardizzazione
        # presente nello script di training V3.
        # Non certifica che gli stessi parametri
        # o lo stesso dataset siano stati usati
        # nella sessione originale.
        mean = x.mean(dim=0, keepdim=True)
        std = x.std(dim=0, keepdim=True)

        normalization = (
            "RICALCOLATA_SUL_GRAFO_CORRENTE"
        )

    std = torch.where(
        std == 0,
        torch.ones_like(std),
        std
    )

    x = (x - mean) / std

    x = torch.nan_to_num(
        x,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    if not torch.isfinite(x).all().item():
        raise RuntimeError(
            "Feature normalizzate non finite"
        )

    print("[OK]", normalization)

    print("[4] RICOSTRUZIONE MODELLO")

    encoder = GCNEncoder(
        in_channels=int(checkpoint["in_channels"]),
        hidden_channels=int(
            checkpoint["hidden_channels"]
        ),
        out_channels=int(
            checkpoint["latent_channels"]
        ),
        dropout=float(
            checkpoint.get("dropout", 0.20)
        )
    )

    model = GAE(encoder)

    model.load_state_dict(
        checkpoint["model_state_dict"],
        strict=True
    )

    model.eval()

    print("[OK] Pesi caricati con strict=True")

    print("[5] INFERENZA CPU SU GRAFO REALE")

    with torch.inference_mode():
        embeddings = model.encode(x, edges)

    expected_shape = (
        n_nodes,
        int(checkpoint["latent_channels"])
    )

    if tuple(embeddings.shape) != expected_shape:
        raise RuntimeError(
            "Dimensioni embedding inattese"
        )

    if not torch.isfinite(embeddings).all().item():
        raise RuntimeError(
            "Embedding NaN o infiniti"
        )

    print(
        "[OK] Dimensioni embedding:",
        tuple(embeddings.shape)
    )

    minimum = float(embeddings.min().item())
    maximum = float(embeddings.max().item())

    average_norm = float(
        torch.linalg.vector_norm(
            embeddings,
            dim=1
        ).mean().item()
    )

    if not all(
        math.isfinite(value)
        for value in (
            minimum,
            maximum,
            average_norm
        )
    ):
        raise RuntimeError(
            "Statistiche embedding non finite"
        )

    report = {
        "status": "PASS",
        "model": "GAE_V3",
        "graph_nodes": n_nodes,
        "graph_edges": n_edges,
        "input_features": n_features,
        "embedding_dimension": (
            int(embeddings.shape[1])
        ),
        "embedding_min": minimum,
        "embedding_max": maximum,
        "mean_embedding_norm": average_norm,
        "normalization": normalization,
        "weights_loaded_strict": True,
        "graph_used_for_inference": (
            "FULL_GRAPH"
        ),
        "original_training_edge_split_reproduced": False,
        "scientific_performance_validated": False,
        "training_performed": False,
        "checkpoint_modified": False,
    }

    print(
        "INFERENCE_RESULT="
        + json.dumps(report)
    )


if __name__ == "__main__":
    main()
'''


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(4 * 1024 * 1024),
            b""
        ):
            digest.update(block)

    return digest.hexdigest()


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 19.3")
    print("INFERENZA REALE GAE V3")
    print("=" * 70, flush=True)

    for path in (MODEL, GRAPH):
        if not path.is_file():
            raise FileNotFoundError(path)

    before = {
        "checkpoint": sha256(MODEL),
        "graph": sha256(GRAPH),
    }

    RUNS.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP19_3_INFERENZA_GAE_V3_{stamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report_path = (
        output / "RAPPORTO_STEP19_3.json"
    )

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "19.3",
        "status": "NOT_STARTED",
        "runtime_image": IMAGE,
        "checkpoint": str(MODEL),
        "graph": str(GRAPH),
        "sha256_before": before,
        "training_performed": False,
    }

    def save_report():
        report_path.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False
            ),
            encoding="utf-8"
        )

    save_report()

    command = [
        "docker", "run",
        "--rm",
        "--pull", "never",
        "--network", "none",
        "--read-only",
        "--tmpfs", "/tmp:rw,nosuid,size=128m",
        "--cap-drop", "ALL",
        "--security-opt", "no-new-privileges",
        "-e", "HOME=/tmp",
        "-e", "OMP_NUM_THREADS=2",
        "-v", f"{ROOT}:/project:ro",
        "-w", "/project",
        "--entrypoint", "python3",
        IMAGE,
        "-B", "-c", RUNNER
    ]

    print("\n[1] AVVIO INFERENZA DOCKER", flush=True)

    try:
        process = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=240,
            check=False
        )

        print(process.stdout, flush=True)

        if process.stderr:
            print(
                process.stderr[-4000:],
                file=sys.stderr,
                flush=True
            )

        if process.returncode != 0:
            raise RuntimeError(
                "Inferenza Docker fallita: "
                + str(process.returncode)
            )

        marker = "INFERENCE_RESULT="

        matches = [
            line[len(marker):]
            for line in process.stdout.splitlines()
            if line.startswith(marker)
        ]

        if len(matches) != 1:
            raise RuntimeError(
                "Risultato di inferenza assente"
            )

        result = json.loads(matches[0])

        if result.get("status") != "PASS":
            raise RuntimeError(
                "Inferenza non superata"
            )

        report["inference"] = result
        report["status"] = (
            "TECHNICAL_INFERENCE_PASS"
        )

    except Exception as error:
        report["status"] = "FAILED"
        report["error"] = str(error)

    finally:
        after = {
            "checkpoint": sha256(MODEL),
            "graph": sha256(GRAPH)
        }

        report["sha256_after"] = after
        report["original_files_unchanged"] = (
            before == after
        )

        if before != after:
            report["status"] = (
                "FAILED_INPUT_HASH_CHANGED"
            )

        save_report()

    print()
    print("=" * 70)
    print("STEP 19.3 - RISULTATO")
    print("=" * 70)

    print("Stato:", report["status"])
    print(
        "File originali invariati:",
        report["original_files_unchanged"]
    )
    print("Report:", report_path)

    print("[OK] Nessun training")
    print("[OK] GUI V3.2 non modificata")
    print("[OK] Motore V3.3 non modificato")

    if report["status"] != "TECHNICAL_INFERENCE_PASS":
        print(
            "\n[ATTENZIONE] Il collaudo GAE V3 "
            "non e stato superato."
        )
        sys.exit(1)

    print(
        "\n[OK] Inferenza tecnica GAE V3 superata"
    )

    print(
        "[NOTA] Le metriche scientifiche "
        "non sono state rivalidate."
    )


if __name__ == "__main__":
    main()