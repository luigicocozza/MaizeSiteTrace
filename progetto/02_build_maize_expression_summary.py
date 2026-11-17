#!/usr/bin/env python3
"""
02_build_maize_expression_summary.py

Crea un file iniziale di expression summary per ncRNA del mais.

Nota metodologica:
in questa fase non vengono inventati valori RNA-seq. Se non sono presenti
matrici quantitative di espressione, lo script produce un summary conservativo
con valori NA e stato "not_quantified_yet".

Input:
- data/processed/maize_sequences_metadata.csv
- data/raw/ncbi_geo_maize_drought_search.json
- data/raw/europe_pmc_maize_drought_ncrna.json
- configs/config_public.yaml

Output:
- data/processed/maize_expression_summary.csv
- outputs/maize_expression_summary_report.json
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
from typing import Any, Dict, List

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Costruisce maize_expression_summary.csv iniziale e conservativo."
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
        "--geo-json",
        default="data/raw/ncbi_geo_maize_drought_search.json",
        help="JSON NCBI GEO/SRA search.",
    )
    parser.add_argument(
        "--europe-pmc-json",
        default="data/raw/europe_pmc_maize_drought_ncrna.json",
        help="JSON Europe PMC letteratura mais/siccità/ncRNA.",
    )
    parser.add_argument(
        "--output-csv",
        default=None,
        help="CSV di output. Se omesso usa il percorso del YAML.",
    )
    parser.add_argument(
        "--report-json",
        default="outputs/maize_expression_summary_report.json",
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


def load_json_if_exists(path: Path) -> Dict[str, Any]:
    if not path.exists():
        logging.warning("File JSON non trovato: %s", path)
        return {}
    with path.open("r", encoding="utf-8") as f:
        obj = json.load(f)
    if not isinstance(obj, dict):
        logging.warning("Il file JSON non contiene un dizionario: %s", path)
        return {}
    return obj


def read_metadata_csv(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Metadata CSV non trovato: {path}")

    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        raise ValueError(f"Metadata CSV vuoto: {path}")

    required = {"rnacentral_id", "rna_type", "sequence_length"}
    missing = required - set(rows[0].keys())
    if missing:
        raise ValueError(f"Colonne mancanti nel metadata CSV: {missing}")

    return rows


def extract_geo_summary(geo: Dict[str, Any]) -> Dict[str, Any]:
    esearch = geo.get("esearchresult", {}) if isinstance(geo, dict) else {}
    idlist = esearch.get("idlist", []) or []

    return {
        "geo_search_count": int(esearch.get("count", 0) or 0),
        "geo_retmax": int(esearch.get("retmax", 0) or 0),
        "geo_ids": idlist,
        "geo_ids_joined": ";".join(idlist),
    }


def extract_literature_summary(europe: Dict[str, Any]) -> Dict[str, Any]:
    result_list = (
        europe.get("resultList", {}).get("result", [])
        if isinstance(europe, dict)
        else []
    )

    titles = []
    years = []
    for rec in result_list:
        title = rec.get("title")
        year = rec.get("pubYear")
        if title:
            titles.append(title)
        if year:
            years.append(str(year))

    return {
        "literature_hit_count": int(europe.get("hitCount", 0) or 0),
        "downloaded_literature_records": len(result_list),
        "sample_literature_titles": titles[:10],
        "sample_literature_years": years[:10],
    }


def build_expression_rows(
    metadata_rows: List[Dict[str, Any]],
    geo_summary: Dict[str, Any],
    literature_summary: Dict[str, Any],
) -> List[Dict[str, Any]]:
    rows = []

    for rec in metadata_rows:
        rnacentral_id = rec.get("rnacentral_id", "")
        rna_type = rec.get("rna_type", "unknown")
        sequence_length = rec.get("sequence_length", "")

        row = {
            "rnacentral_id": rnacentral_id,
            "rna_type": rna_type,
            "sequence_length": sequence_length,

            # Campi quantitativi: non inventati.
            "mean_expression_control": "NA",
            "mean_expression_drought": "NA",
            "log2_fold_change": "NA",
            "adjusted_pvalue": "NA",

            # Stato metodologico.
            "expression_status": "not_quantified_yet",
            "expression_source": "public_GEO_metadata_only",
            "expression_dataset_id": "NA",
            "condition_control_label": "control",
            "condition_stress_label": "drought",

            # Evidenza contestuale: utile per report, non equivalente a espressione.
            "geo_search_count": geo_summary.get("geo_search_count", 0),
            "geo_retmax": geo_summary.get("geo_retmax", 0),
            "geo_candidate_dataset_ids": geo_summary.get("geo_ids_joined", ""),
            "literature_hit_count": literature_summary.get("literature_hit_count", 0),
            "downloaded_literature_records": literature_summary.get(
                "downloaded_literature_records", 0
            ),

            # Flag conservativi per evitare overclaiming.
            "has_quantitative_expression": 0,
            "has_drought_expression_evidence": 0,
            "expression_evidence_note": (
                "Quantitative expression not yet computed. "
                "This row is a placeholder linked to public dataset discovery."
            ),
        }

        rows.append(row)

    return rows


def write_csv(rows: List[Dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "rnacentral_id",
        "rna_type",
        "sequence_length",
        "mean_expression_control",
        "mean_expression_drought",
        "log2_fold_change",
        "adjusted_pvalue",
        "expression_status",
        "expression_source",
        "expression_dataset_id",
        "condition_control_label",
        "condition_stress_label",
        "geo_search_count",
        "geo_retmax",
        "geo_candidate_dataset_ids",
        "literature_hit_count",
        "downloaded_literature_records",
        "has_quantitative_expression",
        "has_drought_expression_evidence",
        "expression_evidence_note",
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
    geo_json = Path(args.geo_json)
    europe_json = Path(args.europe_pmc_json)
    report_json = Path(args.report_json)

    cfg = load_yaml(config_path)

    if args.output_csv:
        output_csv = Path(args.output_csv)
    else:
        output_csv = Path(
            cfg.get("pipeline_outputs", {})
            .get("initial_data", {})
            .get("maize_expression_summary", "data/processed/maize_expression_summary.csv")
        )

    logging.info("Caricamento metadata sequenze: %s", metadata_csv)
    metadata_rows = read_metadata_csv(metadata_csv)

    logging.info("Caricamento GEO JSON: %s", geo_json)
    geo = load_json_if_exists(geo_json)
    geo_summary = extract_geo_summary(geo)

    logging.info("Caricamento Europe PMC JSON: %s", europe_json)
    europe = load_json_if_exists(europe_json)
    literature_summary = extract_literature_summary(europe)

    logging.info("Costruzione expression summary conservativo")
    expression_rows = build_expression_rows(
        metadata_rows=metadata_rows,
        geo_summary=geo_summary,
        literature_summary=literature_summary,
    )

    write_csv(expression_rows, output_csv)

    rna_type_counter = Counter(row["rna_type"] for row in expression_rows)
    status_counter = Counter(row["expression_status"] for row in expression_rows)

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "inputs": {
            "metadata_csv": str(metadata_csv),
            "geo_json": str(geo_json),
            "europe_pmc_json": str(europe_json),
        },
        "outputs": {
            "expression_summary_csv": str(output_csv),
            "report_json": str(report_json),
        },
        "summary": {
            "n_rows": len(expression_rows),
            "rna_type_distribution": dict(rna_type_counter),
            "expression_status_distribution": dict(status_counter),
            "geo_search_count": geo_summary.get("geo_search_count", 0),
            "geo_retmax": geo_summary.get("geo_retmax", 0),
            "literature_hit_count": literature_summary.get("literature_hit_count", 0),
            "downloaded_literature_records": literature_summary.get(
                "downloaded_literature_records", 0
            ),
            "quantitative_expression_available": False,
            "methodological_note": (
                "No quantitative RNA-seq expression values were inferred. "
                "The file is a conservative placeholder for downstream pipeline compatibility."
            ),
        },
    }

    write_report(report, report_json)

    logging.info("File creato: %s", output_csv)
    logging.info("Report creato: %s", report_json)
    logging.info("Righe create: %s", len(expression_rows))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
