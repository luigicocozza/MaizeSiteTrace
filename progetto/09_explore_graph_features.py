#!/usr/bin/env python3
"""
09_explore_graph_features.py

Analisi esplorativa iniziale del grafo PyTorch Geometric, delle tabelle
del grafo e della matrice feature ncRNA.

Output:
- outputs/eda/graph_eda_report.json
- outputs/eda/graph_eda_summary.txt
- outputs/eda/*.csv
- outputs/eda/plots/*.png
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import torch
import yaml

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except Exception:
    plt = None
    MATPLOTLIB_AVAILABLE = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analisi esplorativa del grafo e delle feature ncRNA."
    )
    parser.add_argument("--config", default="configs/config_public.yaml")
    parser.add_argument("--pyg-pt", default="data/processed/pyg/maize_ncrna_graph.pt")
    parser.add_argument("--graph-nodes-csv", default="data/processed/graph_nodes.csv")
    parser.add_argument("--graph-edges-csv", default="data/processed/graph_edges.csv")
    parser.add_argument("--ncrna-features-csv", default="data/processed/ncrna_features.csv")
    parser.add_argument("--output-dir", default="outputs/eda")
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


def write_csv(rows: List[Dict[str, Any]], path: Path, fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def safe_float(value: Any) -> float | None:
    try:
        if value in {"", "NA", None}:
            return None
        return float(value)
    except Exception:
        return None


def summarize_numeric(values: List[float]) -> Dict[str, Any]:
    clean = [v for v in values if v is not None and not math.isnan(v)]
    if not clean:
        return {
            "count": 0,
            "mean": "NA",
            "std": "NA",
            "min": "NA",
            "q25": "NA",
            "median": "NA",
            "q75": "NA",
            "max": "NA",
        }

    clean_sorted = sorted(clean)

    def quantile(q: float) -> float:
        idx = int(round((len(clean_sorted) - 1) * q))
        return clean_sorted[idx]

    if len(clean) > 1:
        std = statistics.stdev(clean)
    else:
        std = 0.0

    return {
        "count": len(clean),
        "mean": round(statistics.mean(clean), 6),
        "std": round(std, 6),
        "min": round(min(clean), 6),
        "q25": round(quantile(0.25), 6),
        "median": round(statistics.median(clean), 6),
        "q75": round(quantile(0.75), 6),
        "max": round(max(clean), 6),
    }


def load_pyg_object(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Dataset PyG non trovato: {path}")

    try:
        obj = torch.load(path, map_location="cpu", weights_only=False)
    except TypeError:
        obj = torch.load(path, map_location="cpu")

    if not isinstance(obj, dict):
        raise ValueError("Il file .pt non contiene un dizionario.")

    if "data" not in obj or "metadata" not in obj:
        raise ValueError("Il file .pt non contiene le chiavi 'data' e 'metadata'.")

    return obj


def get_tensor_shapes(data: Any) -> Dict[str, Any]:
    if hasattr(data, "x"):
        return {
            "x": tuple(data.x.shape),
            "edge_index": tuple(data.edge_index.shape),
            "edge_attr": tuple(data.edge_attr.shape),
            "y": tuple(data.y.shape),
            "train_mask_true": int(data.train_mask.sum()),
            "unlabeled_mask_true": int(data.unlabeled_mask.sum()),
        }

    return {
        "x": tuple(data["x"].shape),
        "edge_index": tuple(data["edge_index"].shape),
        "edge_attr": tuple(data["edge_attr"].shape),
        "y": tuple(data["y"].shape),
        "train_mask_true": int(data["train_mask"].sum()),
        "unlabeled_mask_true": int(data["unlabeled_mask"].sum()),
    }


def counter_to_rows(counter: Counter, key_name: str, value_name: str = "count") -> List[Dict[str, Any]]:
    total = sum(counter.values())
    rows = []
    for key, value in counter.most_common():
        frac = value / total if total else 0
        rows.append(
            {
                key_name: key,
                value_name: value,
                "fraction": round(frac, 6),
            }
        )
    return rows


def compute_degree_tables(edge_rows: List[Dict[str, str]]) -> Dict[str, Dict[str, int]]:
    degree = defaultdict(lambda: {"degree_total": 0, "degree_out": 0, "degree_in": 0})

    for edge in edge_rows:
        src = edge.get("source_node_id", "")
        tgt = edge.get("target_node_id", "")

        if src:
            degree[src]["degree_total"] += 1
            degree[src]["degree_out"] += 1

        if tgt:
            degree[tgt]["degree_total"] += 1
            degree[tgt]["degree_in"] += 1

    return degree


def build_top_ncrna_by_degree(
    feature_rows: List[Dict[str, str]],
    degree_table: Dict[str, Dict[str, int]],
    top_n: int = 30,
) -> List[Dict[str, Any]]:
    rows = []

    for rec in feature_rows:
        node_id = rec.get("node_id", "")
        deg = degree_table.get(node_id, {"degree_total": 0, "degree_out": 0, "degree_in": 0})

        rows.append(
            {
                "rnacentral_id": rec.get("rnacentral_id", ""),
                "node_id": node_id,
                "rna_type": rec.get("rna_type", ""),
                "sequence_length_metadata": rec.get("sequence_length_metadata", ""),
                "gc_frac": rec.get("gc_frac", ""),
                "label": rec.get("label", ""),
                "degree_total": deg["degree_total"],
                "degree_out": deg["degree_out"],
                "degree_in": deg["degree_in"],
                "mean_confidence": rec.get("mean_confidence", ""),
                "mean_weight": rec.get("mean_weight", ""),
            }
        )

    rows.sort(key=lambda x: int(x["degree_total"]), reverse=True)
    return rows[:top_n]


def numeric_feature_summary(feature_rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    exclude = {"rnacentral_id", "node_id", "rna_type"}
    columns = [col for col in feature_rows[0].keys() if col not in exclude]

    summary_rows = []

    for col in columns:
        values = []
        numeric = True

        for row in feature_rows:
            val = safe_float(row.get(col))
            if val is None and row.get(col) not in {"", "NA", None}:
                numeric = False
                break
            if val is not None:
                values.append(val)

        if not numeric:
            continue

        stats = summarize_numeric(values)
        stats["feature"] = col
        summary_rows.append(stats)

    ordered = []
    for row in summary_rows:
        ordered.append(
            {
                "feature": row["feature"],
                "count": row["count"],
                "mean": row["mean"],
                "std": row["std"],
                "min": row["min"],
                "q25": row["q25"],
                "median": row["median"],
                "q75": row["q75"],
                "max": row["max"],
            }
        )

    return ordered


def rna_type_summary(feature_rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    groups = defaultdict(list)

    for row in feature_rows:
        groups[row.get("rna_type", "unknown")].append(row)

    rows = []

    for rna_type, records in sorted(groups.items()):
        lengths = [safe_float(r.get("seq_length_from_fasta")) for r in records]
        gc_values = [safe_float(r.get("gc_frac")) for r in records]
        degree_values = [safe_float(r.get("degree_total")) for r in records]

        length_stats = summarize_numeric([v for v in lengths if v is not None])
        gc_stats = summarize_numeric([v for v in gc_values if v is not None])
        degree_stats = summarize_numeric([v for v in degree_values if v is not None])

        rows.append(
            {
                "rna_type": rna_type,
                "count": len(records),
                "mean_length": length_stats["mean"],
                "median_length": length_stats["median"],
                "mean_gc_frac": gc_stats["mean"],
                "median_gc_frac": gc_stats["median"],
                "mean_degree": degree_stats["mean"],
                "median_degree": degree_stats["median"],
            }
        )

    return rows


def save_bar_plot(rows: List[Dict[str, Any]], x_key: str, y_key: str, title: str, output_path: Path) -> bool:
    if not MATPLOTLIB_AVAILABLE:
        return False

    if not rows:
        return False

    labels = [str(row[x_key]) for row in rows]
    values = [float(row[y_key]) for row in rows]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 6))
    plt.bar(labels, values)
    plt.title(title)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()

    return True


def save_histogram(values: List[float], title: str, output_path: Path, bins: int = 30) -> bool:
    if not MATPLOTLIB_AVAILABLE:
        return False

    clean = [v for v in values if v is not None and not math.isnan(v)]
    if not clean:
        return False

    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 6))
    plt.hist(clean, bins=bins)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()

    return True


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    pyg_path = Path(args.pyg_pt)
    nodes_path = Path(args.graph_nodes_csv)
    edges_path = Path(args.graph_edges_csv)
    features_path = Path(args.ncrna_features_csv)
    output_dir = Path(args.output_dir)
    plots_dir = output_dir / "plots"

    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_yaml(config_path)

    logging.info("Caricamento dataset PyG: %s", pyg_path)
    pyg_obj = load_pyg_object(pyg_path)
    data = pyg_obj["data"]
    metadata = pyg_obj["metadata"]

    tensor_shapes = get_tensor_shapes(data)

    logging.info("Caricamento graph nodes: %s", nodes_path)
    node_rows = read_csv(nodes_path)

    logging.info("Caricamento graph edges: %s", edges_path)
    edge_rows = read_csv(edges_path)

    logging.info("Caricamento ncRNA features: %s", features_path)
    feature_rows = read_csv(features_path)

    node_type_counter = Counter(row.get("node_type", "unknown") for row in node_rows)
    edge_type_counter = Counter(row.get("edge_type", "unknown") for row in edge_rows)
    evidence_source_counter = Counter(row.get("evidence_source", "unknown") for row in edge_rows)
    rna_type_counter = Counter(row.get("rna_type", "unknown") for row in feature_rows)
    label_counter = Counter(row.get("label", "unknown") for row in feature_rows)

    node_type_rows = counter_to_rows(node_type_counter, "node_type")
    edge_type_rows = counter_to_rows(edge_type_counter, "edge_type")
    evidence_rows = counter_to_rows(evidence_source_counter, "evidence_source")
    rna_type_rows = counter_to_rows(rna_type_counter, "rna_type")
    label_rows = counter_to_rows(label_counter, "label")

    degree_table = compute_degree_tables(edge_rows)
    top_degree_rows = build_top_ncrna_by_degree(feature_rows, degree_table, top_n=30)
    feature_summary_rows = numeric_feature_summary(feature_rows)
    rna_summary_rows = rna_type_summary(feature_rows)

    write_csv(node_type_rows, output_dir / "node_type_distribution.csv", ["node_type", "count", "fraction"])
    write_csv(edge_type_rows, output_dir / "edge_type_distribution.csv", ["edge_type", "count", "fraction"])
    write_csv(evidence_rows, output_dir / "evidence_source_distribution.csv", ["evidence_source", "count", "fraction"])
    write_csv(rna_type_rows, output_dir / "rna_type_distribution.csv", ["rna_type", "count", "fraction"])
    write_csv(label_rows, output_dir / "label_distribution.csv", ["label", "count", "fraction"])

    write_csv(
        top_degree_rows,
        output_dir / "top_ncrna_by_degree.csv",
        [
            "rnacentral_id",
            "node_id",
            "rna_type",
            "sequence_length_metadata",
            "gc_frac",
            "label",
            "degree_total",
            "degree_out",
            "degree_in",
            "mean_confidence",
            "mean_weight",
        ],
    )

    write_csv(
        feature_summary_rows,
        output_dir / "ncrna_feature_numeric_summary.csv",
        ["feature", "count", "mean", "std", "min", "q25", "median", "q75", "max"],
    )

    write_csv(
        rna_summary_rows,
        output_dir / "rna_type_summary.csv",
        [
            "rna_type",
            "count",
            "mean_length",
            "median_length",
            "mean_gc_frac",
            "median_gc_frac",
            "mean_degree",
            "median_degree",
        ],
    )

    plots_created = []

    if save_bar_plot(node_type_rows, "node_type", "count", "Node type distribution", plots_dir / "node_type_distribution.png"):
        plots_created.append("plots/node_type_distribution.png")

    if save_bar_plot(edge_type_rows, "edge_type", "count", "Edge type distribution", plots_dir / "edge_type_distribution.png"):
        plots_created.append("plots/edge_type_distribution.png")

    if save_bar_plot(rna_type_rows, "rna_type", "count", "ncRNA type distribution", plots_dir / "rna_type_distribution.png"):
        plots_created.append("plots/rna_type_distribution.png")

    gc_values = [safe_float(row.get("gc_frac")) for row in feature_rows]
    if save_histogram([v for v in gc_values if v is not None], "GC fraction distribution", plots_dir / "gc_fraction_histogram.png"):
        plots_created.append("plots/gc_fraction_histogram.png")

    length_values = [safe_float(row.get("seq_length_from_fasta")) for row in feature_rows]
    if save_histogram([v for v in length_values if v is not None], "Sequence length distribution", plots_dir / "sequence_length_histogram.png"):
        plots_created.append("plots/sequence_length_histogram.png")

    degree_values = [safe_float(row.get("degree_total")) for row in feature_rows]
    if save_histogram([v for v in degree_values if v is not None], "ncRNA degree distribution", plots_dir / "ncrna_degree_histogram.png"):
        plots_created.append("plots/ncrna_degree_histogram.png")

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "project": cfg.get("project", {}),
        "inputs": {
            "pyg_pt": str(pyg_path),
            "graph_nodes_csv": str(nodes_path),
            "graph_edges_csv": str(edges_path),
            "ncrna_features_csv": str(features_path),
        },
        "outputs": {
            "output_dir": str(output_dir),
            "plots_dir": str(plots_dir),
            "report_json": str(output_dir / "graph_eda_report.json"),
            "summary_txt": str(output_dir / "graph_eda_summary.txt"),
        },
        "summary": {
            "tensor_shapes": tensor_shapes,
            "n_nodes": len(node_rows),
            "n_edges": len(edge_rows),
            "n_ncrna_feature_rows": len(feature_rows),
            "node_type_distribution": dict(node_type_counter),
            "edge_type_distribution": dict(edge_type_counter),
            "evidence_source_distribution": dict(evidence_source_counter),
            "rna_type_distribution": dict(rna_type_counter),
            "label_distribution": dict(label_counter),
            "plots_created": plots_created,
            "matplotlib_available": MATPLOTLIB_AVAILABLE,
            "methodological_note": (
                "Exploratory analysis only. Because labels are still unknown, this step "
                "does not evaluate supervised prediction performance."
            ),
        },
    }

    report_path = output_dir / "graph_eda_report.json"
    with report_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    summary_lines = []
    summary_lines.append("GRAPH EDA SUMMARY")
    summary_lines.append("=" * 80)
    summary_lines.append(f"Nodes: {len(node_rows)}")
    summary_lines.append(f"Edges: {len(edge_rows)}")
    summary_lines.append(f"ncRNA feature rows: {len(feature_rows)}")
    summary_lines.append("")
    summary_lines.append(f"Tensor x shape: {tensor_shapes.get('x')}")
    summary_lines.append(f"Tensor edge_index shape: {tensor_shapes.get('edge_index')}")
    summary_lines.append(f"Tensor edge_attr shape: {tensor_shapes.get('edge_attr')}")
    summary_lines.append(f"Train mask true: {tensor_shapes.get('train_mask_true')}")
    summary_lines.append("")
    summary_lines.append(f"Node type distribution: {dict(node_type_counter)}")
    summary_lines.append(f"Edge type distribution: {dict(edge_type_counter)}")
    summary_lines.append(f"RNA type distribution: {dict(rna_type_counter)}")
    summary_lines.append(f"Label distribution: {dict(label_counter)}")
    summary_lines.append("")
    summary_lines.append(f"Plots created: {plots_created}")
    summary_lines.append("")
    summary_lines.append("Methodological note:")
    summary_lines.append(report["summary"]["methodological_note"])

    summary_path = output_dir / "graph_eda_summary.txt"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")

    print("\n".join(summary_lines))

    logging.info("Report EDA: %s", report_path)
    logging.info("Summary EDA: %s", summary_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
