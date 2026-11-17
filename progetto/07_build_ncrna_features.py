#!/usr/bin/env python3
"""
07_build_ncrna_features.py

Costruisce la prima matrice feature per i nodi ncRNA del progetto
maize_drought_ncrna_gnn.

Input:
- data/processed/maize_sequences.fasta
- data/processed/maize_sequences_metadata.csv
- data/processed/maize_expression_summary.csv
- data/processed/maize_labels.csv
- data/processed/graph_edges.csv

Output:
- data/processed/ncrna_features.csv
- outputs/ncrna_features_report.json
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Costruisce ncrna_features.csv per i nodi ncRNA."
    )
    parser.add_argument("--config", default="configs/config_public.yaml")
    parser.add_argument("--fasta", default="data/processed/maize_sequences.fasta")
    parser.add_argument("--metadata-csv", default="data/processed/maize_sequences_metadata.csv")
    parser.add_argument("--expression-csv", default="data/processed/maize_expression_summary.csv")
    parser.add_argument("--labels-csv", default="data/processed/maize_labels.csv")
    parser.add_argument("--graph-edges-csv", default="data/processed/graph_edges.csv")
    parser.add_argument("--output-csv", default="data/processed/ncrna_features.csv")
    parser.add_argument("--report-json", default="outputs/ncrna_features_report.json")
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


def parse_fasta(path: Path) -> Dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(f"FASTA non trovato: {path}")

    records: Dict[str, str] = {}
    current_id = None
    seq_parts: List[str] = []

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(seq_parts).upper()

                header = line[1:]
                current_id = header.split("|")[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    if current_id is not None:
        records[current_id] = "".join(seq_parts).upper()

    if not records:
        raise ValueError(f"Nessuna sequenza letta dal FASTA: {path}")

    return records


def index_by_id(rows: List[Dict[str, str]], key: str) -> Dict[str, Dict[str, str]]:
    return {row.get(key, ""): row for row in rows if row.get(key, "")}


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in {"", "NA", None}:
            return default
        return float(value)
    except Exception:
        return default


def safe_int(value: Any, default: int = 0) -> int:
    try:
        if value in {"", "NA", None}:
            return default
        return int(float(str(value)))
    except Exception:
        return default


def sequence_features(seq: str) -> Dict[str, Any]:
    seq = seq.upper().replace("T", "U")
    length = len(seq)

    counts = Counter(seq)
    a = counts.get("A", 0)
    c = counts.get("C", 0)
    g = counts.get("G", 0)
    u = counts.get("U", 0)
    n = counts.get("N", 0)

    denom = length if length > 0 else 1
    gc = g + c
    au = a + u

    return {
        "seq_length_from_fasta": length,
        "log10_seq_length": round(math.log10(length + 1), 6),
        "a_count": a,
        "c_count": c,
        "g_count": g,
        "u_count": u,
        "n_count": n,
        "a_frac": round(a / denom, 6),
        "c_frac": round(c / denom, 6),
        "g_frac": round(g / denom, 6),
        "u_frac": round(u / denom, 6),
        "n_frac": round(n / denom, 6),
        "gc_count": gc,
        "au_count": au,
        "gc_frac": round(gc / denom, 6),
        "au_frac": round(au / denom, 6),
    }


def build_graph_feature_index(edge_rows: List[Dict[str, str]]) -> Dict[str, Dict[str, Any]]:
    features: Dict[str, Dict[str, Any]] = defaultdict(lambda: {
        "degree_total": 0,
        "degree_out": 0,
        "degree_in": 0,
        "n_direct_biological_edges": 0,
        "n_proxy_edges": 0,
        "sum_confidence": 0.0,
        "sum_weight": 0.0,
        "n_weighted_edges": 0,
        "edge_type_counts": Counter(),
        "evidence_source_counts": Counter(),
    })

    for edge in edge_rows:
        src = edge.get("source_node_id", "")
        tgt = edge.get("target_node_id", "")
        edge_type = edge.get("edge_type", "")
        evidence_source = edge.get("evidence_source", "")
        confidence = safe_float(edge.get("confidence_score"), 0.0)
        weight = safe_float(edge.get("weight"), 0.0)
        is_direct = safe_int(edge.get("is_direct_biological_evidence"), 0)
        is_proxy = safe_int(edge.get("is_proxy_evidence"), 0)

        endpoints = []

        if src.startswith("ncrna:"):
            rid = src.replace("ncrna:", "", 1)
            endpoints.append((rid, "out"))

        if tgt.startswith("ncrna:"):
            rid = tgt.replace("ncrna:", "", 1)
            endpoints.append((rid, "in"))

        for rid, direction in endpoints:
            f = features[rid]
            f["degree_total"] += 1
            if direction == "out":
                f["degree_out"] += 1
            else:
                f["degree_in"] += 1

            f["n_direct_biological_edges"] += is_direct
            f["n_proxy_edges"] += is_proxy
            f["sum_confidence"] += confidence
            f["sum_weight"] += weight
            f["n_weighted_edges"] += 1
            f["edge_type_counts"][edge_type] += 1
            f["evidence_source_counts"][evidence_source] += 1

    final: Dict[str, Dict[str, Any]] = {}

    for rid, f in features.items():
        n_edges = f["n_weighted_edges"] if f["n_weighted_edges"] > 0 else 1

        final[rid] = {
            "degree_total": f["degree_total"],
            "degree_out": f["degree_out"],
            "degree_in": f["degree_in"],
            "n_direct_biological_edges": f["n_direct_biological_edges"],
            "n_proxy_edges": f["n_proxy_edges"],
            "mean_confidence": round(f["sum_confidence"] / n_edges, 6),
            "mean_weight": round(f["sum_weight"] / n_edges, 6),
            "edge_count_has_rna_type": f["edge_type_counts"].get("has_rna_type", 0),
            "edge_count_listed_in_database": f["edge_type_counts"].get("listed_in_database", 0),
            "edge_count_belongs_to_species": f["edge_type_counts"].get("belongs_to_species", 0),
            "edge_count_has_expression_status": f["edge_type_counts"].get("has_expression_status", 0),
            "evidence_count_rnacentral": f["evidence_source_counts"].get("RNAcentral", 0),
            "evidence_count_expression_summary": f["evidence_source_counts"].get("pipeline_expression_summary", 0),
        }

    return final


def build_feature_rows(
    metadata_rows: List[Dict[str, str]],
    expression_rows: List[Dict[str, str]],
    label_rows: List[Dict[str, str]],
    fasta_sequences: Dict[str, str],
    graph_features: Dict[str, Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], List[str]]:
    expression_by_id = index_by_id(expression_rows, "rnacentral_id")
    label_by_id = index_by_id(label_rows, "rnacentral_id")

    rna_types = sorted(set(row.get("rna_type", "unknown") or "unknown" for row in metadata_rows))
    warnings: List[str] = []
    rows: List[Dict[str, Any]] = []

    for rec in metadata_rows:
        rid = rec.get("rnacentral_id", "")
        if not rid:
            continue

        seq = fasta_sequences.get(rid, "")
        if not seq:
            warnings.append(f"Sequenza mancante per {rid}")

        expr = expression_by_id.get(rid, {})
        label = label_by_id.get(rid, {})
        graph = graph_features.get(rid, {})

        rna_type = rec.get("rna_type", "unknown") or "unknown"

        row: Dict[str, Any] = {
            "rnacentral_id": rid,
            "node_id": f"ncrna:{rid}",
            "rna_type": rna_type,
            "sequence_length_metadata": safe_int(rec.get("sequence_length"), 0),
            "label": safe_int(label.get("label"), -1),
            "label_confidence": safe_float(label.get("label_confidence"), 0.0),
            "has_quantitative_expression": safe_int(expr.get("has_quantitative_expression"), 0),
            "has_drought_expression_evidence": safe_int(expr.get("has_drought_expression_evidence"), 0),
            "expression_status_not_quantified_yet": 1 if expr.get("expression_status") == "not_quantified_yet" else 0,
        }

        row.update(sequence_features(seq))

        for rt in rna_types:
            clean_rt = rt.replace(" ", "_").replace("-", "_")
            row[f"rna_type__{clean_rt}"] = 1 if rna_type == rt else 0

        graph_defaults = {
            "degree_total": 0,
            "degree_out": 0,
            "degree_in": 0,
            "n_direct_biological_edges": 0,
            "n_proxy_edges": 0,
            "mean_confidence": 0.0,
            "mean_weight": 0.0,
            "edge_count_has_rna_type": 0,
            "edge_count_listed_in_database": 0,
            "edge_count_belongs_to_species": 0,
            "edge_count_has_expression_status": 0,
            "evidence_count_rnacentral": 0,
            "evidence_count_expression_summary": 0,
        }

        for key, default in graph_defaults.items():
            row[key] = graph.get(key, default)

        rows.append(row)

    return rows, warnings


def write_features_csv(rows: List[Dict[str, Any]], path: Path) -> None:
    if not rows:
        raise ValueError("Nessuna riga feature da scrivere.")

    path.parent.mkdir(parents=True, exist_ok=True)

    fixed_first = [
        "rnacentral_id",
        "node_id",
        "rna_type",
        "sequence_length_metadata",
        "seq_length_from_fasta",
        "label",
        "label_confidence",
    ]

    all_keys = list(rows[0].keys())
    remaining = [k for k in all_keys if k not in fixed_first]
    fieldnames = fixed_first + remaining

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    fasta_path = Path(args.fasta)
    metadata_path = Path(args.metadata_csv)
    expression_path = Path(args.expression_csv)
    labels_path = Path(args.labels_csv)
    graph_edges_path = Path(args.graph_edges_csv)
    output_csv = Path(args.output_csv)
    report_json = Path(args.report_json)

    cfg = load_yaml(config_path)

    logging.info("Caricamento FASTA: %s", fasta_path)
    fasta_sequences = parse_fasta(fasta_path)

    logging.info("Caricamento metadata: %s", metadata_path)
    metadata_rows = read_csv(metadata_path)

    logging.info("Caricamento expression summary: %s", expression_path)
    expression_rows = read_csv(expression_path)

    logging.info("Caricamento labels: %s", labels_path)
    label_rows = read_csv(labels_path)

    logging.info("Caricamento graph edges: %s", graph_edges_path)
    edge_rows = read_csv(graph_edges_path)

    logging.info("Costruzione feature di grafo per ncRNA")
    graph_features = build_graph_feature_index(edge_rows)

    logging.info("Costruzione matrice feature ncRNA")
    feature_rows, warnings = build_feature_rows(
        metadata_rows=metadata_rows,
        expression_rows=expression_rows,
        label_rows=label_rows,
        fasta_sequences=fasta_sequences,
        graph_features=graph_features,
    )

    write_features_csv(feature_rows, output_csv)

    rna_type_counter = Counter(row["rna_type"] for row in feature_rows)
    label_counter = Counter(str(row["label"]) for row in feature_rows)

    numeric_feature_columns = [
        key for key in feature_rows[0].keys()
        if key not in {"rnacentral_id", "node_id", "rna_type"}
    ]

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "project": cfg.get("project", {}),
        "inputs": {
            "fasta": str(fasta_path),
            "metadata_csv": str(metadata_path),
            "expression_csv": str(expression_path),
            "labels_csv": str(labels_path),
            "graph_edges_csv": str(graph_edges_path),
        },
        "outputs": {
            "ncrna_features_csv": str(output_csv),
            "report_json": str(report_json),
        },
        "summary": {
            "n_feature_rows": len(feature_rows),
            "n_fasta_sequences": len(fasta_sequences),
            "n_numeric_feature_columns": len(numeric_feature_columns),
            "rna_type_distribution": dict(rna_type_counter),
            "label_distribution": dict(label_counter),
            "warnings": warnings,
            "methodological_note": (
                "These are initial node features for ncRNA nodes. Labels are still unknown; "
                "the feature matrix is suitable for graph construction, exploratory analysis, "
                "semi-supervised learning or later integration of stronger labels."
            ),
        },
    }

    report_json.parent.mkdir(parents=True, exist_ok=True)
    with report_json.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    logging.info("Feature rows: %s", len(feature_rows))
    logging.info("Feature CSV: %s", output_csv)
    logging.info("Report: %s", report_json)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
