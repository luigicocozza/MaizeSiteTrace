#!/usr/bin/env python3
"""
00d_fetch_more_rnacentral_maize.py

Scarica più pagine RNAcentral per Zea mays e crea un JSON combinato
con più sequenze ncRNA/RNA utilizzabili nella pipeline.

Output:
- data/raw/rnacentral_maize_ncrna_expanded.json
- outputs/rnacentral_expanded_report.json
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import requests
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scarica più pagine RNAcentral per Zea mays."
    )
    parser.add_argument(
        "--config",
        default="configs/config_public.yaml",
        help="File YAML di configurazione.",
    )
    parser.add_argument(
        "--pages",
        type=int,
        default=10,
        help="Numero massimo di pagine RNAcentral da scaricare.",
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=100,
        help="Numero di record per pagina.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="Timeout HTTP in secondi.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.5,
        help="Pausa in secondi tra le richieste.",
    )
    parser.add_argument(
        "--output-json",
        default="data/raw/rnacentral_maize_ncrna_expanded.json",
        help="File JSON combinato di output.",
    )
    parser.add_argument(
        "--report-json",
        default="outputs/rnacentral_expanded_report.json",
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


def load_config(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Config non trovato: {path}")
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError("YAML non valido.")
    return cfg


def fetch_page(
    base_url: str,
    taxon_id: int,
    page: int,
    page_size: int,
    timeout: int,
) -> Dict[str, Any]:
    params = {
        "taxid": taxon_id,
        "format": "json",
        "page": page,
        "page_size": page_size,
    }

    headers = {
        "User-Agent": "maize-drought-ncrna-gnn/0.1.0 (+research; RNAcentral pagination)"
    }

    response = requests.get(
        base_url,
        params=params,
        timeout=timeout,
        headers=headers,
    )
    response.raise_for_status()
    return response.json()


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    output_json = Path(args.output_json)
    report_json = Path(args.report_json)

    output_json.parent.mkdir(parents=True, exist_ok=True)
    report_json.parent.mkdir(parents=True, exist_ok=True)

    cfg = load_config(config_path)

    taxon_id = int(
        cfg.get("scope", {})
        .get("species", {})
        .get("ncbi_taxon_id", 4577)
    )

    base_url = "https://rnacentral.org/api/v1/rna/"

    all_records: List[Dict[str, Any]] = []
    seen_ids = set()
    page_reports = []

    logging.info("Avvio download paginato RNAcentral")
    logging.info("Taxon ID: %s", taxon_id)
    logging.info("Pagine richieste: %s", args.pages)
    logging.info("Page size: %s", args.page_size)

    for page in range(1, args.pages + 1):
        logging.info("Scarico pagina %s/%s", page, args.pages)

        try:
            obj = fetch_page(
                base_url=base_url,
                taxon_id=taxon_id,
                page=page,
                page_size=args.page_size,
                timeout=args.timeout,
            )

            results = obj.get("results", [])
            if not isinstance(results, list):
                results = []

            new_count = 0
            duplicate_count = 0

            for rec in results:
                rid = rec.get("rnacentral_id")
                if not rid:
                    continue
                if rid in seen_ids:
                    duplicate_count += 1
                    continue
                seen_ids.add(rid)
                all_records.append(rec)
                new_count += 1

            page_reports.append(
                {
                    "page": page,
                    "status": "ok",
                    "records_in_page": len(results),
                    "new_records": new_count,
                    "duplicates": duplicate_count,
                    "next": obj.get("next"),
                }
            )

            logging.info(
                "Pagina %s: record=%s, nuovi=%s, duplicati=%s",
                page,
                len(results),
                new_count,
                duplicate_count,
            )

            if not results:
                logging.info("Pagina vuota. Interrompo.")
                break

            if obj.get("next") is None:
                logging.info("RNAcentral non indica pagine successive. Interrompo.")
                break

            time.sleep(args.sleep)

        except Exception as exc:
            logging.error("Errore pagina %s: %s", page, exc)
            page_reports.append(
                {
                    "page": page,
                    "status": "failed",
                    "error": str(exc),
                }
            )
            continue

    combined = {
        "source": "RNAcentral API",
        "taxon_id": taxon_id,
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "results": all_records,
        "page_reports": page_reports,
    }

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(combined, f, indent=2, ensure_ascii=False)

    rna_type_counts = {}
    for rec in all_records:
        rna_type = rec.get("rna_type") or "unknown"
        rna_type_counts[rna_type] = rna_type_counts.get(rna_type, 0) + 1

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "output_json": str(output_json),
        "summary": {
            "pages_requested": args.pages,
            "page_size": args.page_size,
            "pages_processed": len(page_reports),
            "unique_records": len(all_records),
            "rna_type_distribution": rna_type_counts,
        },
        "page_reports": page_reports,
    }

    with report_json.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    logging.info("JSON combinato scritto in: %s", output_json)
    logging.info("Report scritto in: %s", report_json)
    logging.info("Record unici totali: %s", len(all_records))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
