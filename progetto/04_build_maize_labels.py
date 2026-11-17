#!/usr/bin/env python3
"""
04_build_maize_labels.py

Costruisce maize_labels.csv per il progetto maize_drought_ncrna_gnn.

Nota metodologica:
in questa fase non vengono inventate label positive/negative.
Le label vengono impostate a -1 quando mancano evidenze specifiche
per ncRNA drought-responsive.

Output:
- data/processed/maize_labels.csv
- outputs/maize_labels_report.json
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Costruisce proxy/unknown labels per ncRNA del mais."
    )
    parser.add_argument(
        "--config",
        default="configs/config_public.yaml",
        help="File YAML di configurazione.",
    )
    parser.add_argument(
        "--metadata-csv",
        default="data/processed/maize_sequences_metadata.csv",
        help="Metadata CSV delle sequenze ncRNA.",
    )
    parser.add_argument(
        "--expression-csv",
        default="data/processed/maize_expression_summary.csv",
        help="Expression summary CSV.",
    )
    parser.add_argument(
        "--known-edges-csv",
        default="data/processed/maize_known_edges.csv",
        help="CSV degli archi biologici/contextuali.",
    )
    parser.add_argument(
        "--output-csv",
        default=None,
        help="CSV labels di output. Se omesso usa il percorso del YAML.",
    )
    parser.add_argument(
        "--report-json",
        default="outputs/maize_labels_report.json",
        help="Report JSON finale.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Logging dettagliato.",
    )
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
        raise FileNotFoundError(f"Config YAML non trovato: {path}")
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError("Config YAML non valido.")
    return cfg


def read_csv(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV non trovato: {path}")
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"CSV vuoto: {path}")
    return rows


def normalize_ncrna_id(value: str) -> str:
    value = value or ""
    if value.startswith("ncrna:"):
        return value.replace("ncrna:", "", 1)
    return value


def collect_edges_by_ncrna(edges: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    by_id: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    for edge in edges:
        src = edge.get("source_node_id", "")
        tgt = edge.get("target_node_id", "")

        if src.startswith("ncrna:"):
            by_id[normalize_ncrna_id(src)].append(edge)

        if tgt.startswith("ncrna:"):
            by_id[normalize_ncrna_id(tgt)].append(edge)

    return by_id


def collect_expression_by_id(rows: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    return {
        row.get("rnacentral_id"): row
        for row in rows
        if row.get("rnacentral_id")
    }


def safe_float(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if text in {"", "NA", "nan", "None"}:
        return None
    try:
        return float(text)
    except Exception:
        return None


def build_label_row(
    rec: Dict[str, Any],
    expression_by_id: Dict[str, Dict[str, Any]],
    edges_by_ncrna: Dict[str, List[Dict[str, Any]]],
    cfg: Dict[str, Any],
) -> Dict[str, Any]:
    rid = rec.get("rnacentral_id", "")
    rna_type = rec.get("rna_type", "unknown")
    sequence_length = rec.get("sequence_length", "")

    labeling_cfg = cfg.get("labeling", {})
    positive_label = int(labeling_cfg.get("positive_label", 1))
    negative_label = int(labeling_cfg.get("negative_label", 0))
    unknown_label = int(labeling_cfg.get("unknown_label", -1))
    min_positive = int(labeling_cfg.get("minimum_positive_evidence_count", 2))

    expr = expression_by_id.get(rid, {})
    edges = edges_by_ncrna.get(rid, [])

    evidence_count = 0
    drought_expression_evidence = 0
    direct_biological_evidence_count = 0
    proxy_evidence_count = 0
    database_evidence_count = 0
    annotation_evidence_count = 0

    evidence_notes = []

    for edge in edges:
        edge_type = edge.get("edge_type", "")
        evidence_source = edge.get("evidence_source", "")

        if evidence_source in {"RNAcentral"}:
            database_evidence_count += 1

        if edge_type in {"has_rna_type", "belongs_to_species", "listed_in_database"}:
            annotation_evidence_count += 1

        if str(edge.get("is_direct_biological_evidence", "0")) == "1":
            direct_biological_evidence_count += 1

        if str(edge.get("is_proxy_evidence", "0")) == "1":
            proxy_evidence_count += 1

    has_quant_expr = int(str(expr.get("has_quantitative_expression", "0")) == "1")
    has_drought_expr = int(str(expr.get("has_drought_expression_evidence", "0")) == "1")

    log2fc = safe_float(expr.get("log2_fold_change"))
    padj = safe_float(expr.get("adjusted_pvalue"))

    if has_drought_expr:
        drought_expression_evidence += 1
        evidence_count += 1
        evidence_notes.append("drought_expression_evidence")

    if log2fc is not None and padj is not None:
        thresholds = cfg.get("preprocessing", {}).get("expression_filters", {}).get(
            "differential_expression", {}
        )
        log2fc_thr = float(thresholds.get("log2fc_threshold", 1.0))
        padj_thr = float(thresholds.get("adjusted_pvalue_threshold", 0.05))

        if abs(log2fc) >= log2fc_thr and padj <= padj_thr:
            drought_expression_evidence += 1
            evidence_count += 1
            evidence_notes.append("differential_expression_threshold_passed")

    # Annotazione database e tipo RNA sono evidenze di identità, non evidenze drought.
    # Non vengono usate per rendere positivo un candidato.
    label = unknown_label
    label_status = "unknown_unlabeled"
    label_confidence = 0.0

    if evidence_count >= min_positive:
        label = positive_label
        label_status = "positive_proxy_label"
        label_confidence = min(1.0, evidence_count / 5.0)
    else:
        # Non assegniamo 0 automaticamente: assenza di evidenza non è evidenza di assenza.
        label = unknown_label
        label_status = "unknown_pending_drought_specific_evidence"
        label_confidence = 0.0

    return {
        "rnacentral_id": rid,
        "node_id": f"ncrna:{rid}",
        "rna_type": rna_type,
        "sequence_length": sequence_length,
        "label": label,
        "label_status": label_status,
        "label_confidence": f"{label_confidence:.3f}",
        "positive_evidence_count": evidence_count,
        "minimum_positive_evidence_count": min_positive,
        "drought_expression_evidence_count": drought_expression_evidence,
        "direct_biological_evidence_count": direct_biological_evidence_count,
        "proxy_evidence_count": proxy_evidence_count,
        "database_evidence_count": database_evidence_count,
        "annotation_evidence_count": annotation_evidence_count,
        "has_quantitative_expression": has_quant_expr,
        "has_drought_expression_evidence": has_drought_expr,
        "log2_fold_change": expr.get("log2_fold_change", "NA"),
        "adjusted_pvalue": expr.get("adjusted_pvalue", "NA"),
        "evidence_notes": ";".join(evidence_notes) if evidence_notes else "no_drought_specific_evidence_yet",
        "methodological_note": (
            "Unknown label assigned because no ncRNA-specific drought evidence "
            "is available yet. Database/provenance annotations are not treated as "
            "proof of drought responsiveness."
        ),
    }


def write_labels_csv(rows: List[Dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "rnacentral_id",
        "node_id",
        "rna_type",
        "sequence_length",
        "label",
        "label_status",
        "label_confidence",
        "positive_evidence_count",
        "minimum_positive_evidence_count",
        "drought_expression_evidence_count",
        "direct_biological_evidence_count",
        "proxy_evidence_count",
        "database_evidence_count",
        "annotation_evidence_count",
        "has_quantitative_expression",
        "has_drought_expression_evidence",
        "log2_fold_change",
        "adjusted_pvalue",
        "evidence_notes",
        "methodological_note",
    ]

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_report(report: Dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    metadata_csv = Path(args.metadata_csv)
    expression_csv = Path(args.expression_csv)
    known_edges_csv = Path(args.known_edges_csv)
    report_json = Path(args.report_json)

    cfg = load_yaml(config_path)

    if args.output_csv:
        output_csv = Path(args.output_csv)
    else:
        output_csv = Path(
            cfg.get("pipeline_outputs", {})
            .get("initial_data", {})
            .get("maize_labels", "data/processed/maize_labels.csv")
        )

    logging.info("Caricamento metadata: %s", metadata_csv)
    metadata_rows = read_csv(metadata_csv)

    logging.info("Caricamento expression summary: %s", expression_csv)
    expression_rows = read_csv(expression_csv)

    logging.info("Caricamento known edges: %s", known_edges_csv)
    edge_rows = read_csv(known_edges_csv)

    expression_by_id = collect_expression_by_id(expression_rows)
    edges_by_ncrna = collect_edges_by_ncrna(edge_rows)

    label_rows: List[Dict[str, Any]] = []
    for rec in metadata_rows:
        label_rows.append(
            build_label_row(
                rec=rec,
                expression_by_id=expression_by_id,
                edges_by_ncrna=edges_by_ncrna,
                cfg=cfg,
            )
        )

    write_labels_csv(label_rows, output_csv)

    label_counter = Counter(str(row["label"]) for row in label_rows)
    status_counter = Counter(row["label_status"] for row in label_rows)
    rna_type_counter = Counter(row["rna_type"] for row in label_rows)

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "inputs": {
            "metadata_csv": str(metadata_csv),
            "expression_csv": str(expression_csv),
            "known_edges_csv": str(known_edges_csv),
        },
        "outputs": {
            "maize_labels": str(output_csv),
            "report_json": str(report_json),
        },
        "summary": {
            "n_rows": len(label_rows),
            "label_distribution": dict(label_counter),
            "label_status_distribution": dict(status_counter),
            "rna_type_distribution": dict(rna_type_counter),
            "n_with_quantitative_expression": sum(
                int(row["has_quantitative_expression"]) for row in label_rows
            ),
            "n_with_drought_expression_evidence": sum(
                int(row["has_drought_expression_evidence"]) for row in label_rows
            ),
            "methodological_note": (
                "All labels may remain unknown at this stage. This is expected when "
                "no ncRNA-specific drought expression or curated functional evidence "
                "has been integrated yet."
            ),
        },
    }

    write_report(report, report_json)

    logging.info("File labels creato: %s", output_csv)
    logging.info("Report creato: %s", report_json)
    logging.info("Righe labels: %s", len(label_rows))
    logging.info("Distribuzione label: %s", dict(label_counter))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
