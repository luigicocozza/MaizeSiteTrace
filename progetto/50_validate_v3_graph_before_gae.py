#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import torch


BASE_DIR = Path(".").resolve()

PYG_V3 = BASE_DIR / "data/processed/pyg/maize_ncrna_graph_v3_prjna419326.pt"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

SUMMARY_TXT = OUT_DIR / "V3_graph_validation_before_gae_summary.txt"
REPORT_JSON = OUT_DIR / "V3_graph_validation_before_gae_report.json"


def load_pt(path: Path):
    try:
        return torch.load(path, weights_only=False)
    except TypeError:
        return torch.load(path)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    obj = load_pt(PYG_V3)
    data = obj["data"]
    metadata = obj["metadata"]

    n_nodes = int(data.num_nodes)
    n_edges = int(data.edge_index.shape[1])
    x_shape = tuple(data.x.shape)
    edge_attr_shape = tuple(data.edge_attr.shape)

    checks = {}

    checks["has_x"] = data.x is not None
    checks["has_edge_index"] = data.edge_index is not None
    checks["has_edge_attr"] = data.edge_attr is not None
    checks["x_rows_equal_num_nodes"] = int(data.x.shape[0]) == n_nodes
    checks["edge_attr_rows_equal_edges"] = int(data.edge_attr.shape[0]) == n_edges
    checks["edge_type_rows_equal_edges"] = int(data.edge_type.shape[0]) == n_edges
    checks["evidence_source_rows_equal_edges"] = int(data.evidence_source.shape[0]) == n_edges
    checks["node_type_rows_equal_nodes"] = int(data.node_type.shape[0]) == n_nodes
    checks["y_rows_equal_nodes"] = int(data.y.shape[0]) == n_nodes
    checks["train_mask_rows_equal_nodes"] = int(data.train_mask.shape[0]) == n_nodes
    checks["unlabeled_mask_rows_equal_nodes"] = int(data.unlabeled_mask.shape[0]) == n_nodes
    checks["is_ncrna_rows_equal_nodes"] = int(data.is_ncrna_node.shape[0]) == n_nodes

    min_edge_idx = int(data.edge_index.min().item())
    max_edge_idx = int(data.edge_index.max().item())

    checks["edge_index_min_ge_zero"] = min_edge_idx >= 0
    checks["edge_index_max_lt_num_nodes"] = max_edge_idx < n_nodes

    has_nan_x = bool(torch.isnan(data.x).any().item())
    has_inf_x = bool(torch.isinf(data.x).any().item())

    checks["x_no_nan"] = not has_nan_x
    checks["x_no_inf"] = not has_inf_x

    has_nan_edge_attr = bool(torch.isnan(data.edge_attr).any().item())
    has_inf_edge_attr = bool(torch.isinf(data.edge_attr).any().item())

    checks["edge_attr_no_nan"] = not has_nan_edge_attr
    checks["edge_attr_no_inf"] = not has_inf_edge_attr

    train_count = int(data.train_mask.sum().item())
    unlabeled_count = int(data.unlabeled_mask.sum().item())
    ncrna_count = int(data.is_ncrna_node.sum().item())

    status = "PASS" if all(checks.values()) else "FAIL"

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "pyg_v3": str(PYG_V3),
        "status": status,
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "x_shape": list(x_shape),
        "edge_attr_shape": list(edge_attr_shape),
        "min_edge_idx": min_edge_idx,
        "max_edge_idx": max_edge_idx,
        "train_mask_true": train_count,
        "unlabeled_mask_true": unlabeled_count,
        "is_ncrna_true": ncrna_count,
        "metadata_summary": metadata.get("summary", {}),
        "checks": checks,
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("VALIDAZIONE GRAFO V3 PRIMA DEL RETRAINING GAE")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Status: {status}")
    lines.append("")
    lines.append(f"Nodi: {n_nodes}")
    lines.append(f"Archi: {n_edges}")
    lines.append(f"x shape: {x_shape}")
    lines.append(f"edge_attr shape: {edge_attr_shape}")
    lines.append(f"edge_index min: {min_edge_idx}")
    lines.append(f"edge_index max: {max_edge_idx}")
    lines.append(f"train_mask true: {train_count}")
    lines.append(f"unlabeled_mask true: {unlabeled_count}")
    lines.append(f"is_ncrna_node true: {ncrna_count}")
    lines.append("")
    lines.append("Controlli:")
    for k, v in checks.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("Output:")
    lines.append(str(SUMMARY_TXT))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
