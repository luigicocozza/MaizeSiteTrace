#!/usr/bin/env python3
"""
01_build_maize_sequences.py

Costruisce il file FASTA delle sequenze ncRNA del mais a partire
dal JSON RNAcentral scaricato nello step 00.

Input:
- data/raw/rnacentral_maize_ncrna.json
- configs/config_public.yaml

Output:
- data/processed/maize_sequences.fasta
- data/processed/maize_sequences_metadata.csv
- outputs/maize_sequences_report.json
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Costruisce maize_sequences.fasta da RNAcentral JSON."
    )
    parser.add_argument(
        "--config",
        default="configs/config_public.yaml",
        help="File YAML di configurazione.",
    )
    parser.add_argument(
        "--input-json",
        default="data/raw/rnacentral_maize_ncrna.json",
        help="File JSON RNAcentral scaricato nello step 00.",
    )
    parser.add_argument(
        "--output-fasta",
        default=None,
        help="File FASTA di output. Se omesso usa il percorso del YAML.",
    )
    parser.add_argument(
        "--metadata-csv",
        default="data/processed/maize_sequences_metadata.csv",
        help="File CSV metadata associato alle sequenze.",
    )
    parser.add_argument(
        "--report-json",
        default="outputs/maize_sequences_report.json",
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
        raise ValueError("Il file YAML non contiene un dizionario valido.")

    return cfg


def load_rnacentral_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"File RNAcentral non trovato: {path}")

    with path.open("r", encoding="utf-8") as f:
        obj = json.load(f)

    if not isinstance(obj, dict):
        raise ValueError("Il JSON RNAcentral non è un dizionario.")

    if "results" not in obj:
        raise ValueError("Il JSON RNAcentral non contiene la chiave 'results'.")

    if not isinstance(obj["results"], list):
        raise ValueError("La chiave 'results' del JSON RNAcentral non è una lista.")

    return obj


def clean_sequence(seq: str) -> str:
    seq = seq.strip().upper()
    seq = re.sub(r"\s+", "", seq)
    return seq


def is_valid_sequence(seq: str, allowed_alphabet: str) -> bool:
    allowed = set(allowed_alphabet.upper())
    return bool(seq) and all(base in allowed for base in seq.upper())


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def safe_header_value(value: Any) -> str:
    if value is None:
        return "NA"
    text = str(value)
    text = text.replace(" ", "_")
    text = text.replace("|", "_")
    text = text.replace(";", "_")
    return text


def extract_records(
    rnacentral_json: Dict[str, Any],
    cfg: Dict[str, Any],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    sequence_filters = cfg.get("preprocessing", {}).get("sequence_filters", {})

    min_length = int(sequence_filters.get("min_length", 18))
    max_length = int(sequence_filters.get("max_length", 10000))
    allowed_alphabet = sequence_filters.get("allowed_alphabet", "ACGTUNacgtun")
    remove_ambiguous = bool(sequence_filters.get("remove_ambiguous_sequences", True))

    raw_results = rnacentral_json.get("results", [])

    accepted: List[Dict[str, Any]] = []
    seen_ids = set()
    seen_md5 = set()

    counters = Counter()
    rna_type_counter = Counter()
    reject_reasons = Counter()

    for rec in raw_results:
        counters["raw_records"] += 1

        rnacentral_id = rec.get("rnacentral_id")
        sequence = rec.get("sequence")
        rna_type = rec.get("rna_type")
        length_reported = rec.get("length")
        md5 = rec.get("md5")
        description = rec.get("description")
        url = rec.get("url")
        publications = rec.get("publications")
        xrefs = rec.get("xrefs")
        distinct_databases = rec.get("distinct_databases")
        distinct_organisms = rec.get("distinct_organisms")

        if not rnacentral_id:
            reject_reasons["missing_rnacentral_id"] += 1
            continue

        if rnacentral_id in seen_ids:
            reject_reasons["duplicate_rnacentral_id"] += 1
            continue

        if not sequence:
            reject_reasons["missing_sequence"] += 1
            continue

        seq = clean_sequence(sequence)
        seq_len = len(seq)

        if seq_len < min_length:
            reject_reasons["too_short"] += 1
            continue

        if seq_len > max_length:
            reject_reasons["too_long"] += 1
            continue

        if remove_ambiguous and not is_valid_sequence(seq, allowed_alphabet):
            reject_reasons["invalid_alphabet"] += 1
            continue

        if md5 and md5 in seen_md5:
            reject_reasons["duplicate_md5"] += 1
            continue

        seen_ids.add(rnacentral_id)
        if md5:
            seen_md5.add(md5)

        rna_type_counter[str(rna_type or "unknown")] += 1

        accepted.append(
            {
                "rnacentral_id": rnacentral_id,
                "sequence": seq,
                "sequence_length": seq_len,
                "length_reported": length_reported,
                "rna_type": rna_type or "unknown",
                "md5": md5 or "",
                "description": description or "",
                "url": url or "",
                "publications_count": len(publications) if isinstance(publications, list) else 0,
                "xrefs_count": len(xrefs) if isinstance(xrefs, list) else 0,
                "distinct_databases_count": len(distinct_databases) if isinstance(distinct_databases, list) else 0,
                "distinct_organisms_count": len(distinct_organisms) if isinstance(distinct_organisms, list) else 0,
            }
        )

    report = {
        "raw_records": counters["raw_records"],
        "accepted_records": len(accepted),
        "rejected_records": counters["raw_records"] - len(accepted),
        "reject_reasons": dict(reject_reasons),
        "rna_type_distribution": dict(rna_type_counter),
        "filters": {
            "min_length": min_length,
            "max_length": max_length,
            "allowed_alphabet": allowed_alphabet,
            "remove_ambiguous_sequences": remove_ambiguous,
        },
    }

    return accepted, report


def write_fasta(records: List[Dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as f:
        for rec in records:
            header = (
                f">{safe_header_value(rec['rnacentral_id'])}"
                f"|rna_type={safe_header_value(rec['rna_type'])}"
                f"|length={safe_header_value(rec['sequence_length'])}"
                f"|md5={safe_header_value(rec['md5'])}"
            )
            f.write(header + "\n")
            f.write(wrap_sequence(rec["sequence"]) + "\n")


def write_metadata_csv(records: List[Dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "rnacentral_id",
        "rna_type",
        "sequence_length",
        "length_reported",
        "md5",
        "description",
        "url",
        "publications_count",
        "xrefs_count",
        "distinct_databases_count",
        "distinct_organisms_count",
    ]

    with output_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            writer.writerow({k: rec.get(k, "") for k in fieldnames})


def write_report(report: Dict[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    input_json = Path(args.input_json)
    metadata_csv = Path(args.metadata_csv)
    report_json = Path(args.report_json)

    cfg = load_yaml(config_path)

    if args.output_fasta:
        output_fasta = Path(args.output_fasta)
    else:
        output_fasta = Path(
            cfg.get("pipeline_outputs", {})
            .get("initial_data", {})
            .get("maize_sequences", "data/processed/maize_sequences.fasta")
        )

    logging.info("Caricamento config: %s", config_path)
    logging.info("Caricamento RNAcentral JSON: %s", input_json)

    rnacentral_json = load_rnacentral_json(input_json)

    records, extraction_report = extract_records(rnacentral_json, cfg)

    if not records:
        raise ValueError("Nessuna sequenza valida estratta da RNAcentral.")

    logging.info("Sequenze accettate: %s", len(records))
    logging.info("Scrittura FASTA: %s", output_fasta)
    write_fasta(records, output_fasta)

    logging.info("Scrittura metadata CSV: %s", metadata_csv)
    write_metadata_csv(records, metadata_csv)

    final_report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "input_json": str(input_json),
        "outputs": {
            "fasta": str(output_fasta),
            "metadata_csv": str(metadata_csv),
            "report_json": str(report_json),
        },
        "summary": extraction_report,
    }

    write_report(final_report, report_json)

    logging.info("Report scritto in: %s", report_json)
    logging.info("Completato correttamente.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
