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
OUT_TXT = OUT_DIR / "V3_graph_input_inspection_summary.txt"
OUT_JSON = OUT_DIR / "V3_graph_input_inspection_report.json"


def load_pyg(path: Path):
    try:
        return torch.load(path, weights_only=False)
    except TypeError:
        return torch.load(path)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("ISPEZIONE INPUT GRAFO V3")
    print("=" * 80)

    data = load_pyg(PYG_PATH)

    nodes = pd.read_csv(GRAPH_NODES)
    edges = pd.read_csv(GRAPH_EDGES)
    ncrna_features = pd.read_csv(NCRNA_FEATURES)
    v3_ranking = pd.read_csv(V3_RANKING)

    info = {
        "pyg_path": str(PYG_PATH),
        "graph_nodes_csv": str(GRAPH_NODES),
        "graph_edges_csv": str(GRAPH_EDGES),
        "ncrna_features_csv": str(NCRNA_FEATURES),
        "v3_ranking_csv": str(V3_RANKING),
        "pyg_class": str(type(data)),
        "pyg_num_nodes": int(data.num_nodes),
        "pyg_num_edges": int(data.edge_index.shape[1]),
        "pyg_x_shape": list(data.x.shape) if hasattr(data, "x") and data.x is not None else None,
        "graph_nodes_shape": list(nodes.shape),
        "graph_edges_shape": list(edges.shape),
        "ncrna_features_shape": list(ncrna_features.shape),
        "v3_ranking_shape": list(v3_ranking.shape),
        "pyg_keys": list(data.keys()) if hasattr(data, "keys") else [],
        "graph_nodes_columns": list(nodes.columns),
        "graph_edges_columns": list(edges.columns),
        "ncrna_features_columns": list(ncrna_features.columns),
        "v3_ranking_columns": list(v3_ranking.columns),
    }

    lines = []
    lines.append("ISPEZIONE INPUT GRAFO V3")
    lines.append("=" * 80)
    lines.append("")
    lines.append("1) Oggetto PyG esistente")
    lines.append(f"Path: {PYG_PATH}")
    lines.append(f"Classe: {type(data)}")
    lines.append(f"Nodi PyG: {data.num_nodes}")
    lines.append(f"Archi PyG: {data.edge_index.shape[1]}")
    lines.append(f"Feature matrix x: {tuple(data.x.shape) if hasattr(data, 'x') and data.x is not None else None}")
    lines.append(f"Chiavi PyG: {list(data.keys()) if hasattr(data, 'keys') else []}")
    lines.append("")
    lines.append("2) Tabelle grafo V2")
    lines.append(f"graph_nodes.csv: {nodes.shape}")
    lines.append(f"graph_edges.csv: {edges.shape}")
    lines.append(f"ncrna_features.csv: {ncrna_features.shape}")
    lines.append("")
    lines.append("3) Ranking sperimentale V3")
    lines.append(f"V3 ranking: {v3_ranking.shape}")
    lines.append("")
    lines.append("4) Colonne graph_nodes.csv")
    lines.append(", ".join(nodes.columns))
    lines.append("")
    lines.append("5) Colonne graph_edges.csv")
    lines.append(", ".join(edges.columns))
    lines.append("")
    lines.append("6) Colonne ncrna_features.csv")
    lines.append(", ".join(ncrna_features.columns))
    lines.append("")
    lines.append("7) Colonne V3 ranking")
    lines.append(", ".join(v3_ranking.columns))
    lines.append("")
    lines.append("8) Prime righe graph_nodes")
    lines.append(nodes.head(5).to_string(index=False))
    lines.append("")
    lines.append("9) Prime righe graph_edges")
    lines.append(edges.head(5).to_string(index=False))
    lines.append("")
    lines.append("10) Prime righe V3 ranking")
    lines.append(v3_ranking.head(5).to_string(index=False))

    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    OUT_JSON.write_text(json.dumps(info, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n".join(lines))
    print("")
    print("Output creati:")
    print(OUT_TXT)
    print(OUT_JSON)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
