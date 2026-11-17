#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import torch


BASE_DIR = Path(".").resolve()

PYG_PATH = BASE_DIR / "data/processed/pyg/maize_ncrna_graph.pt"
GRAPH_NODES = BASE_DIR / "data/processed/graph_nodes.csv"
GRAPH_EDGES = BASE_DIR / "data/processed/graph_edges.csv"
NCRNA_FEATURES = BASE_DIR / "data/processed/ncrna_features.csv"
V3_RANKING = BASE_DIR / "outputs/prjna419326/PRJNA419326_experimental_srna_v3_ranking.csv"

OUT_DIR = BASE_DIR / "outputs/prjna419326"
OUT_TXT = OUT_DIR / "V3_graph_metadata_detailed_inspection_summary.txt"
OUT_JSON = OUT_DIR / "V3_graph_metadata_detailed_inspection_report.json"


def load_pt(path: Path):
    try:
        return torch.load(path, weights_only=False)
    except TypeError:
        return torch.load(path)


def describe_value(v):
    if hasattr(v, "shape"):
        return f"{type(v)} shape={tuple(v.shape)}"
    if isinstance(v, dict):
        return f"{type(v)} keys={list(v.keys())[:20]}"
    if isinstance(v, list):
        return f"{type(v)} len={len(v)} first={v[:5]}"
    return f"{type(v)} value={str(v)[:200]}"


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    obj = load_pt(PYG_PATH)
    data = obj["data"]
    metadata = obj.get("metadata", {})

    nodes = pd.read_csv(GRAPH_NODES)
    edges = pd.read_csv(GRAPH_EDGES)
    features = pd.read_csv(NCRNA_FEATURES)
    v3 = pd.read_csv(V3_RANKING)

    lines = []

    lines.append("ISPEZIONE DETTAGLIATA GRAFO V2 PER COSTRUZIONE V3")
    lines.append("=" * 80)

    lines.append("")
    lines.append("1) Oggetto PyG interno")
    lines.append(f"Tipo data: {type(data)}")
    lines.append(f"Numero nodi: {data.num_nodes}")
    lines.append(f"Numero archi: {data.edge_index.shape[1]}")
    lines.append(f"x shape: {tuple(data.x.shape)}")

    lines.append("")
    lines.append("2) Chiavi PyG e shape")
    for key in data.keys():
        value = getattr(data, key)
        lines.append(f"- {key}: {describe_value(value)}")

    lines.append("")
    lines.append("3) Metadata principali")
    lines.append(f"Metadata keys: {list(metadata.keys())}")

    for mk, mv in metadata.items():
        lines.append("")
        lines.append(f"Metadata[{mk}]: {describe_value(mv)}")

        if isinstance(mv, dict):
            for sk, sv in list(mv.items())[:20]:
                lines.append(f"  - {sk}: {describe_value(sv)}")

    lines.append("")
    lines.append("4) graph_nodes.csv")
    lines.append(f"Shape: {nodes.shape}")
    lines.append("Colonne:")
    lines.append(", ".join(nodes.columns))
    lines.append("")
    lines.append(nodes.head(10).to_string(index=False))

    lines.append("")
    lines.append("5) graph_edges.csv")
    lines.append(f"Shape: {edges.shape}")
    lines.append("Colonne:")
    lines.append(", ".join(edges.columns))
    lines.append("")
    lines.append(edges.head(10).to_string(index=False))

    lines.append("")
    lines.append("6) ncrna_features.csv")
    lines.append(f"Shape: {features.shape}")
    lines.append("Colonne:")
    lines.append(", ".join(features.columns))
    lines.append("")
    lines.append(features.head(5).to_string(index=False))

    lines.append("")
    lines.append("7) ranking small RNA V3")
    lines.append(f"Shape: {v3.shape}")
    lines.append("Colonne:")
    lines.append(", ".join(v3.columns))
    lines.append("")
    lines.append(v3.head(5).to_string(index=False))

    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    report = {
        "pyg_path": str(PYG_PATH),
        "data_type": str(type(data)),
        "num_nodes": int(data.num_nodes),
        "num_edges": int(data.edge_index.shape[1]),
        "x_shape": list(data.x.shape),
        "data_keys": list(data.keys()),
        "metadata_keys": list(metadata.keys()),
        "graph_nodes_shape": list(nodes.shape),
        "graph_nodes_columns": list(nodes.columns),
        "graph_edges_shape": list(edges.shape),
        "graph_edges_columns": list(edges.columns),
        "ncrna_features_shape": list(features.shape),
        "ncrna_features_columns": list(features.columns),
        "v3_ranking_shape": list(v3.shape),
        "v3_ranking_columns": list(v3.columns),
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n".join(lines))
    print("")
    print("Output creati:")
    print(OUT_TXT)
    print(OUT_JSON)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
