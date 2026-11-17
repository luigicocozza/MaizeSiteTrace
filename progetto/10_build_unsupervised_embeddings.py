#!/usr/bin/env python3
"""
10_build_unsupervised_embeddings.py

Costruisce embedding non supervisionati dei nodi del grafo.

Metodo:
1. Carica il dataset PyTorch Geometric.
2. Standardizza le feature dei nodi.
3. Propaga le feature lungo il grafo tramite media dei vicini.
4. Applica PCA/SVD per ottenere embedding a bassa dimensione.
5. Esporta embedding per tutti i nodi e per i soli ncRNA.

Output:
- data/processed/embeddings/node_embeddings.csv
- data/processed/embeddings/ncrna_embeddings.csv
- outputs/embeddings/unsupervised_embeddings_report.json
- outputs/embeddings/unsupervised_embeddings_summary.txt
- outputs/embeddings/plots/all_nodes_embedding_2d.png
- outputs/embeddings/plots/ncrna_embedding_2d.png
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import torch
import yaml

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except Exception:
    PCA = None
    StandardScaler = None
    SKLEARN_AVAILABLE = False

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except Exception:
    plt = None
    MATPLOTLIB_AVAILABLE = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Costruisce embedding non supervisionati del grafo ncRNA."
    )
    parser.add_argument("--config", default="configs/config_public.yaml")
    parser.add_argument("--pyg-pt", default="data/processed/pyg/maize_ncrna_graph.pt")
    parser.add_argument("--graph-nodes-csv", default="data/processed/graph_nodes.csv")
    parser.add_argument("--ncrna-features-csv", default="data/processed/ncrna_features.csv")
    parser.add_argument("--output-dir", default="data/processed/embeddings")
    parser.add_argument("--report-dir", default="outputs/embeddings")
    parser.add_argument("--embedding-dim", type=int, default=16)
    parser.add_argument("--propagation-steps", type=int, default=2)
    parser.add_argument("--alpha", type=float, default=0.50)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)


def load_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        obj = yaml.safe_load(f)
    return obj if isinstance(obj, dict) else {}


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV non trovato: {path}")
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"CSV vuoto: {path}")
    return rows


def write_csv(rows: List[Dict[str, Any]], path: Path, fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_pyg_dataset(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Dataset PyG non trovato: {path}")

    try:
        obj = torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        obj = torch.load(path, map_location="cpu")

    if not isinstance(obj, dict):
        raise ValueError("Il file .pt non contiene un dizionario.")

    if "data" not in obj or "metadata" not in obj:
        raise ValueError("Il file .pt non contiene 'data' e 'metadata'.")

    return obj


def get_tensor(data: Any, key: str) -> torch.Tensor:
    if hasattr(data, key):
        return getattr(data, key)
    return data[key]


def build_node_lookup_from_metadata(metadata: Dict[str, Any], n_nodes: int) -> Tuple[Dict[int, str], Dict[int, str]]:
    mappings = metadata.get("mappings", {})

    idx_to_node_id_raw = mappings.get("idx_to_node_id", {})
    idx_to_node_type_raw = mappings.get("idx_to_node_type", {})

    node_type_tensor = None

    idx_to_node_id: Dict[int, str] = {}
    idx_to_node_type: Dict[int, str] = {}

    for i in range(n_nodes):
        idx_to_node_id[i] = idx_to_node_id_raw.get(str(i), f"node_{i}")

    # Il node_type viene ricostruito più avanti anche dal tensore node_type.
    for key, value in idx_to_node_type_raw.items():
        try:
            idx_to_node_type[int(key)] = value
        except Exception:
            continue

    return idx_to_node_id, idx_to_node_type


def build_node_metadata_table(
    node_rows: List[Dict[str, str]],
) -> Dict[str, Dict[str, str]]:
    table = {}
    for row in node_rows:
        node_id = row.get("node_id", "")
        if node_id:
            table[node_id] = row
    return table


def standardize_features(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

    if SKLEARN_AVAILABLE and StandardScaler is not None:
        scaler = StandardScaler()
        return scaler.fit_transform(x)

    mean = x.mean(axis=0, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    return (x - mean) / std


def build_normalized_adjacency(edge_index: np.ndarray, n_nodes: int) -> np.ndarray:
    adjacency = np.zeros((n_nodes, n_nodes), dtype=np.float64)

    for src, dst in edge_index.T:
        src = int(src)
        dst = int(dst)

        if src < 0 or dst < 0 or src >= n_nodes or dst >= n_nodes:
            continue

        adjacency[src, dst] = 1.0
        adjacency[dst, src] = 1.0

    # self-loop
    np.fill_diagonal(adjacency, 1.0)

    degree = adjacency.sum(axis=1, keepdims=True)
    degree[degree == 0] = 1.0

    normalized = adjacency / degree
    return normalized


def propagate_features(
    x: np.ndarray,
    edge_index: np.ndarray,
    propagation_steps: int,
    alpha: float,
) -> np.ndarray:
    n_nodes = x.shape[0]
    adjacency = build_normalized_adjacency(edge_index=edge_index, n_nodes=n_nodes)

    propagated = x.copy()

    for _ in range(max(0, propagation_steps)):
        neighbor_mean = adjacency @ propagated
        propagated = alpha * x + (1.0 - alpha) * neighbor_mean

    return propagated


def reduce_dimensionality(
    x: np.ndarray,
    embedding_dim: int,
) -> Tuple[np.ndarray, List[float], str]:
    n_components = min(embedding_dim, x.shape[0], x.shape[1])

    if n_components < 2:
        raise ValueError("Numero componenti troppo basso per creare embedding 2D/ND.")

    if SKLEARN_AVAILABLE and PCA is not None:
        pca = PCA(n_components=n_components, random_state=42)
        z = pca.fit_transform(x)
        explained = [float(v) for v in pca.explained_variance_ratio_]
        method = "PCA_sklearn"
        return z, explained, method

    x_centered = x - x.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(x_centered, full_matrices=False)
    z = u[:, :n_components] * s[:n_components]

    variances = s ** 2
    total = float(variances.sum()) if float(variances.sum()) > 0 else 1.0
    explained = [float(v / total) for v in variances[:n_components]]

    method = "SVD_numpy"
    return z, explained, method


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in {"", "NA", None}:
            return default
        return float(value)
    except Exception:
        return default


def build_ncrna_feature_lookup(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    out = {}
    for row in rows:
        node_id = row.get("node_id", "")
        if node_id:
            out[node_id] = row
    return out


def build_embedding_rows(
    embeddings: np.ndarray,
    y: np.ndarray,
    is_ncrna: np.ndarray,
    idx_to_node_id: Dict[int, str],
    node_metadata: Dict[str, Dict[str, str]],
    ncrna_feature_lookup: Dict[str, Dict[str, str]],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    all_rows: List[Dict[str, Any]] = []
    ncrna_rows: List[Dict[str, Any]] = []

    emb_cols = [f"emb_{i + 1:03d}" for i in range(embeddings.shape[1])]

    for idx in range(embeddings.shape[0]):
        node_id = idx_to_node_id[idx]
        meta = node_metadata.get(node_id, {})
        feature_meta = ncrna_feature_lookup.get(node_id, {})

        row: Dict[str, Any] = {
            "node_index": idx,
            "node_id": node_id,
            "node_type": meta.get("node_type", "unknown"),
            "name": meta.get("name", ""),
            "is_ncrna": int(bool(is_ncrna[idx])),
            "label": int(y[idx]),
            "rna_type": feature_meta.get("rna_type", meta.get("rna_type", "")),
            "sequence_length": feature_meta.get("seq_length_from_fasta", meta.get("sequence_length", "")),
            "gc_frac": feature_meta.get("gc_frac", ""),
        }

        for col, value in zip(emb_cols, embeddings[idx]):
            row[col] = round(float(value), 8)

        all_rows.append(row)

        if bool(is_ncrna[idx]):
            ncrna_rows.append(row)

    return all_rows, ncrna_rows


def save_embedding_plot(
    rows: List[Dict[str, Any]],
    output_path: Path,
    title: str,
    color_by: str = "node_type",
) -> bool:
    if not MATPLOTLIB_AVAILABLE or plt is None:
        return False

    if not rows:
        return False

    if "emb_001" not in rows[0] or "emb_002" not in rows[0]:
        return False

    output_path.parent.mkdir(parents=True, exist_ok=True)

    groups: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        group = str(row.get(color_by, "unknown"))
        groups.setdefault(group, []).append(row)

    plt.figure(figsize=(9, 7))

    for group, group_rows in groups.items():
        xs = [safe_float(r.get("emb_001"), 0.0) for r in group_rows]
        ys = [safe_float(r.get("emb_002"), 0.0) for r in group_rows]
        plt.scatter(xs, ys, s=24, alpha=0.75, label=group)

    plt.title(title)
    plt.xlabel("Embedding 1")
    plt.ylabel("Embedding 2")
    if len(groups) <= 12:
        plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()

    return True


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    pyg_path = Path(args.pyg_pt)
    graph_nodes_path = Path(args.graph_nodes_csv)
    ncrna_features_path = Path(args.ncrna_features_csv)
    output_dir = Path(args.output_dir)
    report_dir = Path(args.report_dir)
    plots_dir = report_dir / "plots"

    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_yaml(config_path)

    logging.info("Caricamento dataset PyG: %s", pyg_path)
    obj = load_pyg_dataset(pyg_path)
    data = obj["data"]
    metadata = obj["metadata"]

    x_tensor = get_tensor(data, "x")
    edge_index_tensor = get_tensor(data, "edge_index")
    y_tensor = get_tensor(data, "y")
    is_ncrna_tensor = get_tensor(data, "is_ncrna_node")

    x = x_tensor.detach().cpu().numpy()
    edge_index = edge_index_tensor.detach().cpu().numpy()
    y = y_tensor.detach().cpu().numpy()
    is_ncrna = is_ncrna_tensor.detach().cpu().numpy().astype(bool)

    n_nodes = x.shape[0]

    logging.info("Caricamento graph nodes: %s", graph_nodes_path)
    node_rows = read_csv(graph_nodes_path)
    node_metadata = build_node_metadata_table(node_rows)

    logging.info("Caricamento ncRNA features: %s", ncrna_features_path)
    ncrna_feature_rows = read_csv(ncrna_features_path)
    ncrna_feature_lookup = build_ncrna_feature_lookup(ncrna_feature_rows)

    idx_to_node_id, _ = build_node_lookup_from_metadata(metadata, n_nodes=n_nodes)

    logging.info("Standardizzazione feature")
    x_scaled = standardize_features(x)

    logging.info("Propagazione feature lungo il grafo")
    x_graph = propagate_features(
        x=x_scaled,
        edge_index=edge_index,
        propagation_steps=args.propagation_steps,
        alpha=args.alpha,
    )

    logging.info("Riduzione dimensionale non supervisionata")
    embeddings, explained_variance, method = reduce_dimensionality(
        x=x_graph,
        embedding_dim=args.embedding_dim,
    )

    all_rows, ncrna_rows = build_embedding_rows(
        embeddings=embeddings,
        y=y,
        is_ncrna=is_ncrna,
        idx_to_node_id=idx_to_node_id,
        node_metadata=node_metadata,
        ncrna_feature_lookup=ncrna_feature_lookup,
    )

    emb_cols = [f"emb_{i + 1:03d}" for i in range(embeddings.shape[1])]
    fieldnames = [
        "node_index",
        "node_id",
        "node_type",
        "name",
        "is_ncrna",
        "label",
        "rna_type",
        "sequence_length",
        "gc_frac",
    ] + emb_cols

    node_embeddings_csv = output_dir / "node_embeddings.csv"
    ncrna_embeddings_csv = output_dir / "ncrna_embeddings.csv"

    write_csv(all_rows, node_embeddings_csv, fieldnames)
    write_csv(ncrna_rows, ncrna_embeddings_csv, fieldnames)

    plots_created = []

    if save_embedding_plot(
        all_rows,
        plots_dir / "all_nodes_embedding_2d.png",
        title="All graph nodes - unsupervised embedding",
        color_by="node_type",
    ):
        plots_created.append("plots/all_nodes_embedding_2d.png")

    if save_embedding_plot(
        ncrna_rows,
        plots_dir / "ncrna_embedding_2d.png",
        title="ncRNA nodes - unsupervised embedding",
        color_by="rna_type",
    ):
        plots_created.append("plots/ncrna_embedding_2d.png")

    node_type_counter = Counter(row["node_type"] for row in all_rows)
    rna_type_counter = Counter(row["rna_type"] for row in ncrna_rows)
    label_counter = Counter(str(row["label"]) for row in all_rows)

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "project": cfg.get("project", {}),
        "inputs": {
            "pyg_pt": str(pyg_path),
            "graph_nodes_csv": str(graph_nodes_path),
            "ncrna_features_csv": str(ncrna_features_path),
        },
        "outputs": {
            "node_embeddings_csv": str(node_embeddings_csv),
            "ncrna_embeddings_csv": str(ncrna_embeddings_csv),
            "report_json": str(report_dir / "unsupervised_embeddings_report.json"),
            "summary_txt": str(report_dir / "unsupervised_embeddings_summary.txt"),
            "plots_dir": str(plots_dir),
        },
        "summary": {
            "n_nodes": len(all_rows),
            "n_ncrna_nodes": len(ncrna_rows),
            "input_feature_dim": int(x.shape[1]),
            "embedding_dim_requested": args.embedding_dim,
            "embedding_dim_created": int(embeddings.shape[1]),
            "graph_edges": int(edge_index.shape[1]),
            "propagation_steps": args.propagation_steps,
            "alpha": args.alpha,
            "dimensionality_reduction_method": method,
            "sklearn_available": SKLEARN_AVAILABLE,
            "matplotlib_available": MATPLOTLIB_AVAILABLE,
            "explained_variance_ratio": [round(float(v), 8) for v in explained_variance],
            "explained_variance_first_2": round(float(sum(explained_variance[:2])), 8),
            "explained_variance_total_created": round(float(sum(explained_variance)), 8),
            "node_type_distribution": dict(node_type_counter),
            "ncrna_rna_type_distribution": dict(rna_type_counter),
            "label_distribution": dict(label_counter),
            "plots_created": plots_created,
            "methodological_note": (
                "Embeddings are unsupervised and do not represent supervised drought-response "
                "prediction. They are suitable for similarity search, clustering, visualization, "
                "candidate prioritization after adding stronger evidence, and initialization of "
                "later graph neural network experiments."
            ),
        },
    }

    report_path = report_dir / "unsupervised_embeddings_report.json"
    with report_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    summary_lines = []
    summary_lines.append("UNSUPERVISED EMBEDDINGS SUMMARY")
    summary_lines.append("=" * 80)
    summary_lines.append(f"Nodes: {len(all_rows)}")
    summary_lines.append(f"ncRNA nodes: {len(ncrna_rows)}")
    summary_lines.append(f"Input feature dim: {x.shape[1]}")
    summary_lines.append(f"Embedding dim created: {embeddings.shape[1]}")
    summary_lines.append(f"Graph edges: {edge_index.shape[1]}")
    summary_lines.append(f"Method: {method}")
    summary_lines.append(f"Propagation steps: {args.propagation_steps}")
    summary_lines.append(f"Alpha: {args.alpha}")
    summary_lines.append(f"Explained variance first 2 dims: {report['summary']['explained_variance_first_2']}")
    summary_lines.append(f"Explained variance total created dims: {report['summary']['explained_variance_total_created']}")
    summary_lines.append("")
    summary_lines.append(f"Node type distribution: {dict(node_type_counter)}")
    summary_lines.append(f"ncRNA RNA type distribution: {dict(rna_type_counter)}")
    summary_lines.append(f"Label distribution: {dict(label_counter)}")
    summary_lines.append("")
    summary_lines.append(f"Plots created: {plots_created}")
    summary_lines.append("")
    summary_lines.append("Methodological note:")
    summary_lines.append(report["summary"]["methodological_note"])

    summary_path = report_dir / "unsupervised_embeddings_summary.txt"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")

    print("\n".join(summary_lines))

    logging.info("Node embeddings: %s", node_embeddings_csv)
    logging.info("ncRNA embeddings: %s", ncrna_embeddings_csv)
    logging.info("Report: %s", report_path)
    logging.info("Summary: %s", summary_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
