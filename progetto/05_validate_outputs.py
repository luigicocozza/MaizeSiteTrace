#!/usr/bin/env python3
"""
05_validate_outputs.py

Valida formalmente i primi output del progetto maize_drought_ncrna_gnn.

Input principali:
- data/processed/maize_sequences.fasta
- data/processed/maize_sequences_metadata.csv
- data/processed/maize_expression_summary.csv
- data/processed/maize_known_edges.csv
- data/processed/maize_labels.csv

Output:
- outputs/validation_report.json
- outputs/validation_summary.txt
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
from typing import Dict, List, Tuple, Any

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Valida output iniziali del progetto maize drought ncRNA GNN."
    )
    parser.add_argument("--config", default="configs/config_public.yaml")
    parser.add_argument("--fasta", default="data/processed/maize_sequences.fasta")
    parser.add_argument("--metadata-csv", default="data/processed/maize_sequences_metadata.csv")
    parser.add_argument("--expression-csv", default="data/processed/maize_expression_summary.csv")
    parser.add_argument("--edges-csv", default="data/processed/maize_known_edges.csv")
    parser.add_argument("--labels-csv", default="data/processed/maize_labels.csv")
    parser.add_argument("--report-json", default="outputs/validation_report.json")
    parser.add_argument("--summary-txt", default="outputs/validation_summary.txt")
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


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def parse_fasta(path: Path) -> Tuple[List[Dict[str, Any]], List[str]]:
    records = []
    problems = []

    if not path.exists():
        problems.append(f"FASTA non trovato: {path}")
        return records, problems

    current_header = None
    current_seq_parts = []

    with path.open("r", encoding="utf-8") as f:
        for line_number, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    records.append(make_fasta_record(current_header, current_seq_parts))
                current_header = line[1:]
                current_seq_parts = []
            else:
                if current_header is None:
                    problems.append(f"Sequenza senza header alla linea {line_number}")
                else:
                    current_seq_parts.append(line)

    if current_header is not None:
        records.append(make_fasta_record(current_header, current_seq_parts))

    return records, problems


def make_fasta_record(header: str, seq_parts: List[str]) -> Dict[str, Any]:
    sequence = "".join(seq_parts).upper()
    first_field = header.split("|")[0]
    return {
        "header": header,
        "rnacentral_id": first_field,
        "sequence": sequence,
        "length": len(sequence),
    }


def file_status(path: Path) -> Dict[str, Any]:
    return {
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else 0,
    }


def validate_required_columns(
    rows: List[Dict[str, str]],
    required_columns: List[str],
    file_label: str,
    errors: List[str],
) -> None:
    if not rows:
        errors.append(f"{file_label}: file vuoto o non leggibile.")
        return

    present = set(rows[0].keys())
    missing = [col for col in required_columns if col not in present]

    if missing:
        errors.append(f"{file_label}: colonne mancanti: {missing}")


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    fasta_path = Path(args.fasta)
    metadata_path = Path(args.metadata_csv)
    expression_path = Path(args.expression_csv)
    edges_path = Path(args.edges_csv)
    labels_path = Path(args.labels_csv)
    report_json = Path(args.report_json)
    summary_txt = Path(args.summary_txt)

    report_json.parent.mkdir(parents=True, exist_ok=True)
    summary_txt.parent.mkdir(parents=True, exist_ok=True)

    cfg = load_yaml(config_path)

    errors: List[str] = []
    warnings: List[str] = []
    notes: List[str] = []

    files = {
        "config": file_status(config_path),
        "fasta": file_status(fasta_path),
        "metadata_csv": file_status(metadata_path),
        "expression_csv": file_status(expression_path),
        "edges_csv": file_status(edges_path),
        "labels_csv": file_status(labels_path),
    }

    for name, status in files.items():
        if not status["exists"]:
            errors.append(f"File mancante: {name} -> {status['path']}")
        elif status["size_bytes"] == 0:
            errors.append(f"File vuoto: {name} -> {status['path']}")

    fasta_records, fasta_problems = parse_fasta(fasta_path)
    errors.extend(fasta_problems)

    metadata_rows = read_csv_rows(metadata_path)
    expression_rows = read_csv_rows(expression_path)
    edge_rows = read_csv_rows(edges_path)
    label_rows = read_csv_rows(labels_path)

    validate_required_columns(
        metadata_rows,
        ["rnacentral_id", "rna_type", "sequence_length"],
        "metadata_csv",
        errors,
    )
    validate_required_columns(
        expression_rows,
        ["rnacentral_id", "expression_status", "has_quantitative_expression"],
        "expression_csv",
        errors,
    )
    validate_required_columns(
        edge_rows,
        ["source_node_id", "target_node_id", "edge_type", "evidence_source"],
        "edges_csv",
        errors,
    )
    validate_required_columns(
        label_rows,
        ["rnacentral_id", "node_id", "label", "label_status"],
        "labels_csv",
        errors,
    )

    fasta_ids = [r["rnacentral_id"] for r in fasta_records]
    fasta_id_set = set(fasta_ids)
    metadata_ids = [r.get("rnacentral_id", "") for r in metadata_rows]
    metadata_id_set = set(metadata_ids)
    expression_ids = [r.get("rnacentral_id", "") for r in expression_rows]
    expression_id_set = set(expression_ids)
    label_ids = [r.get("rnacentral_id", "") for r in label_rows]
    label_id_set = set(label_ids)

    duplicate_fasta_ids = [item for item, count in Counter(fasta_ids).items() if count > 1]
    duplicate_metadata_ids = [item for item, count in Counter(metadata_ids).items() if count > 1]
    duplicate_expression_ids = [item for item, count in Counter(expression_ids).items() if count > 1]
    duplicate_label_ids = [item for item, count in Counter(label_ids).items() if count > 1]

    if duplicate_fasta_ids:
        errors.append(f"Duplicati nel FASTA: {duplicate_fasta_ids[:10]}")
    if duplicate_metadata_ids:
        errors.append(f"Duplicati nel metadata CSV: {duplicate_metadata_ids[:10]}")
    if duplicate_expression_ids:
        errors.append(f"Duplicati nell'expression CSV: {duplicate_expression_ids[:10]}")
    if duplicate_label_ids:
        errors.append(f"Duplicati nelle labels CSV: {duplicate_label_ids[:10]}")

    allowed_alphabet = set("ACGTUN")
    invalid_sequences = []
    for rec in fasta_records:
        invalid_chars = sorted(set(rec["sequence"]) - allowed_alphabet)
        if invalid_chars:
            invalid_sequences.append(
                {
                    "rnacentral_id": rec["rnacentral_id"],
                    "invalid_chars": invalid_chars,
                }
            )

    if invalid_sequences:
        errors.append(f"Sequenze con caratteri non validi: {invalid_sequences[:10]}")

    if fasta_id_set != metadata_id_set:
        errors.append(
            "Gli ID FASTA non coincidono con gli ID del metadata CSV. "
            f"Solo FASTA: {len(fasta_id_set - metadata_id_set)}; "
            f"solo metadata: {len(metadata_id_set - fasta_id_set)}"
        )

    if metadata_id_set != expression_id_set:
        errors.append(
            "Gli ID metadata non coincidono con gli ID expression summary. "
            f"Solo metadata: {len(metadata_id_set - expression_id_set)}; "
            f"solo expression: {len(expression_id_set - metadata_id_set)}"
        )

    if metadata_id_set != label_id_set:
        errors.append(
            "Gli ID metadata non coincidono con gli ID labels. "
            f"Solo metadata: {len(metadata_id_set - label_id_set)}; "
            f"solo labels: {len(label_id_set - metadata_id_set)}"
        )

    ncrna_edge_ids = set()
    for edge in edge_rows:
        for key in ["source_node_id", "target_node_id"]:
            value = edge.get(key, "")
            if value.startswith("ncrna:"):
                ncrna_edge_ids.add(value.replace("ncrna:", "", 1))

    unknown_edge_ncrna_ids = ncrna_edge_ids - metadata_id_set
    if unknown_edge_ncrna_ids:
        errors.append(
            f"Archi con ncRNA non presenti nel metadata: {list(unknown_edge_ncrna_ids)[:10]}"
        )

    label_distribution = Counter(row.get("label", "") for row in label_rows)
    label_status_distribution = Counter(row.get("label_status", "") for row in label_rows)
    expression_status_distribution = Counter(row.get("expression_status", "") for row in expression_rows)
    edge_type_distribution = Counter(row.get("edge_type", "") for row in edge_rows)

    if label_distribution.get("-1", 0) == len(label_rows):
        warnings.append(
            "Tutte le label sono -1/unknown. Corretto in questa fase, ma non sufficiente "
            "per addestrare un classificatore supervisionato finale."
        )

    if expression_status_distribution.get("not_quantified_yet", 0) == len(expression_rows):
        warnings.append(
            "Tutte le righe expression sono not_quantified_yet. Corretto come placeholder, "
            "ma servirà integrare RNA-seq quantitativo o proxy evidence più forte."
        )

    if len(fasta_records) < 100:
        warnings.append(
            f"Numero sequenze basso per training robusto: {len(fasta_records)}."
        )

    if not errors:
        if warnings:
            validation_status = "pass_with_warnings"
        else:
            validation_status = "pass"
    else:
        validation_status = "failed"

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "project": cfg.get("project", {}),
        "files": files,
        "summary": {
            "validation_status": validation_status,
            "n_errors": len(errors),
            "n_warnings": len(warnings),
            "n_fasta_records": len(fasta_records),
            "n_metadata_rows": len(metadata_rows),
            "n_expression_rows": len(expression_rows),
            "n_edges": len(edge_rows),
            "n_label_rows": len(label_rows),
            "n_unique_fasta_ids": len(fasta_id_set),
            "n_unique_metadata_ids": len(metadata_id_set),
            "n_unique_expression_ids": len(expression_id_set),
            "n_unique_label_ids": len(label_id_set),
            "label_distribution": dict(label_distribution),
            "label_status_distribution": dict(label_status_distribution),
            "expression_status_distribution": dict(expression_status_distribution),
            "edge_type_distribution": dict(edge_type_distribution),
        },
        "errors": errors,
        "warnings": warnings,
        "notes": notes,
    }

    with report_json.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    lines = []
    lines.append("VALIDATION SUMMARY")
    lines.append("=" * 80)
    lines.append(f"Status: {validation_status}")
    lines.append(f"Errors: {len(errors)}")
    lines.append(f"Warnings: {len(warnings)}")
    lines.append("")
    lines.append(f"FASTA records: {len(fasta_records)}")
    lines.append(f"Metadata rows: {len(metadata_rows)}")
    lines.append(f"Expression rows: {len(expression_rows)}")
    lines.append(f"Known edges: {len(edge_rows)}")
    lines.append(f"Label rows: {len(label_rows)}")
    lines.append("")
    lines.append(f"Label distribution: {dict(label_distribution)}")
    lines.append(f"Expression status distribution: {dict(expression_status_distribution)}")
    lines.append(f"Edge type distribution: {dict(edge_type_distribution)}")
    lines.append("")

    if errors:
        lines.append("ERRORS")
        lines.append("-" * 80)
        for err in errors:
            lines.append(f"- {err}")
        lines.append("")

    if warnings:
        lines.append("WARNINGS")
        lines.append("-" * 80)
        for warn in warnings:
            lines.append(f"- {warn}")
        lines.append("")

    with summary_txt.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    logging.info("Report JSON: %s", report_json)
    logging.info("Summary TXT: %s", summary_txt)

    return 0 if validation_status in {"pass", "pass_with_warnings"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
