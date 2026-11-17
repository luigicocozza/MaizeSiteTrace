#!/usr/bin/env python3
"""
08_build_pyg_dataset.py

Costruisce il primo dataset PyTorch Geometric a partire da:
- graph_nodes.csv
- graph_edges.csv
- ncrna_features.csv

Output:
- data/processed/pyg/maize_ncrna_graph.pt
- outputs/pyg_dataset_mappings.json
- outputs/pyg_dataset_report.json
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import torch
import yaml

try:
    from torch_geometric.data import Data
    PYG_AVAILABLE = True
except Exception:
    Data = None
    PYG_AVAILABLE = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Costruisce dataset PyTorch Geometric per maize_drought_ncrna_gnn."
    )
    parser.add_argument("--config", default="configs/config_public.yaml")
    parser.add_argument("--graph-nodes-csv", default="data/processed/graph_nodes.csv")
    parser.add_argument("--graph-edges-csv", default="data/processed/graph_edges.csv")
    parser.add_argument("--ncrna-features-csv", default="data/processed/ncrna_features.csv")
    parser.add_argument("--output-pt", default="data/processed/pyg/maize_ncrna_graph.pt")
    parser.add_argument("--mappings-json", default="outputs/pyg_dataset_mappings.json")
    parser.add_argument("--report-json", default="outputs/pyg_dataset_report.json")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


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


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in {"", "NA", None}:
            return default
        return float(value)
    except Exception:
        return default


def safe_int(value: Any, default: int = -1) -> int:
    try:
        if value in {"", "NA", None}:
            return default
        return int(float(str(value)))
    except Exception:
        return default


def index_ncrna_features(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    indexed = {}
    for row in rows:
        node_id = row.get("node_id", "")
        if node_id:
            indexed[node_id] = row
    return indexed


def detect_ncrna_feature_columns(rows: List[Dict[str, str]]) -> List[str]:
    exclude = {
        "rnacentral_id",
        "node_id",
        "rna_type",
        "label",
        "label_confidence",
    }

    first = rows[0]
    candidate_columns = []

    for col in first.keys():
        if col in exclude:
            continue

        # Manteniamo solo colonne numeriche.
        ok_numeric = True
        for row in rows[:20]:
            try:
                value = row.get(col, "")
                if value in {"", "NA", None}:
                    continue
                float(value)
            except Exception:
                ok_numeric = False
                break

        if ok_numeric:
            candidate_columns.append(col)

    return candidate_columns


def build_node_matrices(
    node_rows: List[Dict[str, str]],
    ncrna_feature_rows: List[Dict[str, str]],
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, Dict[str, Any]]:
    node_ids = [row.get("node_id", "") for row in node_rows]
    node_id_to_idx = {node_id: i for i, node_id in enumerate(node_ids)}

    node_types = sorted(set(row.get("node_type", "unknown") or "unknown" for row in node_rows))
    node_type_to_idx = {node_type: i for i, node_type in enumerate(node_types)}

    ncrna_by_node = index_ncrna_features(ncrna_feature_rows)
    ncrna_feature_columns = detect_ncrna_feature_columns(ncrna_feature_rows)

    feature_matrix: List[List[float]] = []
    labels: List[int] = []
    node_type_codes: List[int] = []
    is_ncrna_node: List[int] = []

    for row in node_rows:
        node_id = row.get("node_id", "")
        node_type = row.get("node_type", "unknown") or "unknown"

        node_type_one_hot = [
            1.0 if node_type == nt else 0.0
            for nt in node_types
        ]

        ncrna_features = [0.0 for _ in ncrna_feature_columns]
        label_value = -1
        is_ncrna = 0

        if node_id in ncrna_by_node:
            is_ncrna = 1
            feat_row = ncrna_by_node[node_id]
            ncrna_features = [
                safe_float(feat_row.get(col), 0.0)
                for col in ncrna_feature_columns
            ]
            label_value = safe_int(feat_row.get("label"), -1)

        feature_matrix.append(node_type_one_hot + ncrna_features)
        labels.append(label_value)
        node_type_codes.append(node_type_to_idx[node_type])
        is_ncrna_node.append(is_ncrna)

    x = torch.tensor(feature_matrix, dtype=torch.float32)
    y = torch.tensor(labels, dtype=torch.long)
    node_type_tensor = torch.tensor(node_type_codes, dtype=torch.long)
    is_ncrna_tensor = torch.tensor(is_ncrna_node, dtype=torch.bool)

    train_mask = y >= 0
    unlabeled_mask = y < 0

    mappings = {
        "node_id_to_idx": node_id_to_idx,
        "idx_to_node_id": {str(v): k for k, v in node_id_to_idx.items()},
        "node_type_to_idx": node_type_to_idx,
        "idx_to_node_type": {str(v): k for k, v in node_type_to_idx.items()},
        "node_type_feature_columns": [f"node_type__{nt}" for nt in node_types],
        "ncrna_feature_columns": ncrna_feature_columns,
        "all_feature_columns": [f"node_type__{nt}" for nt in node_types] + ncrna_feature_columns,
    }

    return x, y, train_mask, unlabeled_mask, node_type_tensor, is_ncrna_tensor, mappings


def build_edge_tensors(
    edge_rows: List[Dict[str, str]],
    node_id_to_idx: Dict[str, int],
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, Dict[str, Any], List[str]]:
    edge_types = sorted(set(row.get("edge_type", "unknown") or "unknown" for row in edge_rows))
    evidence_sources = sorted(set(row.get("evidence_source", "unknown") or "unknown" for row in edge_rows))

    edge_type_to_idx = {edge_type: i for i, edge_type in enumerate(edge_types)}
    evidence_source_to_idx = {source: i for i, source in enumerate(evidence_sources)}

    sources = []
    targets = []
    edge_attrs = []
    edge_type_codes = []
    evidence_source_codes = []
    skipped_edges = []

    for edge in edge_rows:
        src = edge.get("source_node_id", "")
        tgt = edge.get("target_node_id", "")

        if src not in node_id_to_idx or tgt not in node_id_to_idx:
            skipped_edges.append(f"{src}->{tgt}")
            continue

        sources.append(node_id_to_idx[src])
        targets.append(node_id_to_idx[tgt])

        weight = safe_float(edge.get("weight"), 1.0)
        confidence = safe_float(edge.get("confidence_score"), 1.0)
        is_direct = safe_float(edge.get("is_direct_biological_evidence"), 0.0)
        is_proxy = safe_float(edge.get("is_proxy_evidence"), 0.0)

        edge_attrs.append([weight, confidence, is_direct, is_proxy])

        edge_type = edge.get("edge_type", "unknown") or "unknown"
        evidence_source = edge.get("evidence_source", "unknown") or "unknown"

        edge_type_codes.append(edge_type_to_idx[edge_type])
        evidence_source_codes.append(evidence_source_to_idx[evidence_source])

    edge_index = torch.tensor([sources, targets], dtype=torch.long)
    edge_attr = torch.tensor(edge_attrs, dtype=torch.float32)
    edge_type_tensor = torch.tensor(edge_type_codes, dtype=torch.long)
    evidence_source_tensor = torch.tensor(evidence_source_codes, dtype=torch.long)

    mappings = {
        "edge_type_to_idx": edge_type_to_idx,
        "idx_to_edge_type": {str(v): k for k, v in edge_type_to_idx.items()},
        "evidence_source_to_idx": evidence_source_to_idx,
        "idx_to_evidence_source": {str(v): k for k, v in evidence_source_to_idx.items()},
        "edge_attr_columns": [
            "weight",
            "confidence_score",
            "is_direct_biological_evidence",
            "is_proxy_evidence",
        ],
    }

    return (
        edge_index,
        edge_attr,
        edge_type_tensor,
        evidence_source_tensor,
        mappings,
        skipped_edges,
    )


def save_dataset(
    output_pt: Path,
    x: torch.Tensor,
    edge_index: torch.Tensor,
    edge_attr: torch.Tensor,
    y: torch.Tensor,
    train_mask: torch.Tensor,
    unlabeled_mask: torch.Tensor,
    node_type_tensor: torch.Tensor,
    is_ncrna_tensor: torch.Tensor,
    edge_type_tensor: torch.Tensor,
    evidence_source_tensor: torch.Tensor,
    metadata: Dict[str, Any],
) -> str:
    output_pt.parent.mkdir(parents=True, exist_ok=True)

    if PYG_AVAILABLE and Data is not None:
        data = Data(
            x=x,
            edge_index=edge_index,
            edge_attr=edge_attr,
            y=y,
            train_mask=train_mask,
            unlabeled_mask=unlabeled_mask,
        )

        data.node_type = node_type_tensor
        data.is_ncrna_node = is_ncrna_tensor
        data.edge_type = edge_type_tensor
        data.evidence_source = evidence_source_tensor

        dataset_type = "torch_geometric.data.Data"
    else:
        data = {
            "x": x,
            "edge_index": edge_index,
            "edge_attr": edge_attr,
            "y": y,
            "train_mask": train_mask,
            "unlabeled_mask": unlabeled_mask,
            "node_type": node_type_tensor,
            "is_ncrna_node": is_ncrna_tensor,
            "edge_type": edge_type_tensor,
            "evidence_source": evidence_source_tensor,
        }

        dataset_type = "plain_tensor_dictionary"

    torch.save(
        {
            "data": data,
            "metadata": metadata,
        },
        output_pt,
    )

    return dataset_type


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    graph_nodes_path = Path(args.graph_nodes_csv)
    graph_edges_path = Path(args.graph_edges_csv)
    ncrna_features_path = Path(args.ncrna_features_csv)
    output_pt = Path(args.output_pt)
    mappings_json = Path(args.mappings_json)
    report_json = Path(args.report_json)

    cfg = load_yaml(config_path)

    logging.info("PyG disponibile: %s", PYG_AVAILABLE)

    logging.info("Caricamento graph nodes: %s", graph_nodes_path)
    node_rows = read_csv(graph_nodes_path)

    logging.info("Caricamento graph edges: %s", graph_edges_path)
    edge_rows = read_csv(graph_edges_path)

    logging.info("Caricamento ncRNA features: %s", ncrna_features_path)
    ncrna_feature_rows = read_csv(ncrna_features_path)

    logging.info("Costruzione tensori nodi")
    (
        x,
        y,
        train_mask,
        unlabeled_mask,
        node_type_tensor,
        is_ncrna_tensor,
        node_mappings,
    ) = build_node_matrices(
        node_rows=node_rows,
        ncrna_feature_rows=ncrna_feature_rows,
    )

    logging.info("Costruzione tensori archi")
    (
        edge_index,
        edge_attr,
        edge_type_tensor,
        evidence_source_tensor,
        edge_mappings,
        skipped_edges,
    ) = build_edge_tensors(
        edge_rows=edge_rows,
        node_id_to_idx=node_mappings["node_id_to_idx"],
    )

    mappings = {}
    mappings.update(node_mappings)
    mappings.update(edge_mappings)

    label_distribution = Counter(str(int(v)) for v in y.tolist())
    node_type_distribution = Counter(row.get("node_type", "unknown") for row in node_rows)
    edge_type_distribution = Counter(row.get("edge_type", "unknown") for row in edge_rows)

    metadata = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "project": cfg.get("project", {}),
        "source_files": {
            "config": str(config_path),
            "graph_nodes_csv": str(graph_nodes_path),
            "graph_edges_csv": str(graph_edges_path),
            "ncrna_features_csv": str(ncrna_features_path),
        },
        "mappings": mappings,
        "summary": {
            "n_nodes": int(x.shape[0]),
            "n_edges": int(edge_index.shape[1]),
            "n_node_features": int(x.shape[1]),
            "n_edge_features": int(edge_attr.shape[1]),
            "n_labeled_nodes": int(train_mask.sum().item()),
            "n_unlabeled_nodes": int(unlabeled_mask.sum().item()),
            "n_ncrna_nodes": int(is_ncrna_tensor.sum().item()),
            "label_distribution": dict(label_distribution),
            "node_type_distribution": dict(node_type_distribution),
            "edge_type_distribution": dict(edge_type_distribution),
            "skipped_edges_count": len(skipped_edges),
            "methodological_note": (
                "All labels are currently unknown/unlabeled. The PyG object is ready "
                "for graph construction, embeddings, unsupervised/semi-supervised "
                "experiments, and later supervised training after stronger labels are added."
            ),
        },
    }

    dataset_type = save_dataset(
        output_pt=output_pt,
        x=x,
        edge_index=edge_index,
        edge_attr=edge_attr,
        y=y,
        train_mask=train_mask,
        unlabeled_mask=unlabeled_mask,
        node_type_tensor=node_type_tensor,
        is_ncrna_tensor=is_ncrna_tensor,
        edge_type_tensor=edge_type_tensor,
        evidence_source_tensor=evidence_source_tensor,
        metadata=metadata,
    )

    metadata["summary"]["dataset_type"] = dataset_type
    metadata["summary"]["pyg_available"] = PYG_AVAILABLE

    mappings_json.parent.mkdir(parents=True, exist_ok=True)
    with mappings_json.open("w", encoding="utf-8") as f:
        json.dump(mappings, f, indent=2, ensure_ascii=False)

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "outputs": {
            "pyg_dataset": str(output_pt),
            "mappings_json": str(mappings_json),
            "report_json": str(report_json),
        },
        "summary": metadata["summary"],
        "skipped_edges": skipped_edges[:50],
    }

    report_json.parent.mkdir(parents=True, exist_ok=True)
    with report_json.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    logging.info("Dataset salvato: %s", output_pt)
    logging.info("Mappings salvati: %s", mappings_json)
    logging.info("Report salvato: %s", report_json)
    logging.info("x shape: %s", tuple(x.shape))
    logging.info("edge_index shape: %s", tuple(edge_index.shape))
    logging.info("edge_attr shape: %s", tuple(edge_attr.shape))
    logging.info("n labeled nodes: %s", int(train_mask.sum().item()))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
