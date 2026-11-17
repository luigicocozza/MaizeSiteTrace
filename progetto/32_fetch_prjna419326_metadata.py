#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import urllib.parse
import urllib.request
from pathlib import Path


BASE_DIR = Path(".").resolve()
PROJECT = "PRJNA419326"

META_DIR = BASE_DIR / "data/public_ena" / PROJECT / "metadata"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

META_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

TSV_PATH = META_DIR / f"{PROJECT}_ena_run_metadata.tsv"
JSON_PATH = META_DIR / f"{PROJECT}_ena_run_metadata_summary.json"
TXT_PATH = OUT_DIR / f"{PROJECT}_metadata_summary.txt"

FIELDS = [
    "study_accession",
    "secondary_study_accession",
    "study_title",
    "sample_accession",
    "secondary_sample_accession",
    "sample_alias",
    "sample_title",
    "sample_description",
    "experiment_accession",
    "experiment_title",
    "run_accession",
    "scientific_name",
    "tax_id",
    "instrument_platform",
    "instrument_model",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
    "submitted_ftp",
    "submitted_md5",
    "submitted_bytes",
    "sra_ftp",
    "sra_md5",
    "sra_bytes",
    "first_public",
    "last_updated",
]


def fetch_ena_metadata() -> str:
    base_url = "https://www.ebi.ac.uk/ena/portal/api/filereport"

    params = {
        "accession": PROJECT,
        "result": "read_run",
        "fields": ",".join(FIELDS),
        "format": "tsv",
        "download": "true",
        "limit": "0",
    }

    url = base_url + "?" + urllib.parse.urlencode(params)

    print("ENA URL:")
    print(url)
    print("")

    with urllib.request.urlopen(url, timeout=120) as response:
        text = response.read().decode("utf-8", errors="replace")

    return text


def parse_tsv(text: str) -> list[dict]:
    lines = [line for line in text.splitlines() if line.strip()]

    if not lines:
        return []

    reader = csv.DictReader(lines, delimiter="\t")
    return list(reader)


def safe_int(value: str) -> int:
    if value is None:
        return 0

    value = str(value).strip()

    if not value:
        return 0

    total = 0

    for part in value.split(";"):
        part = part.strip()
        if not part:
            continue
        try:
            total += int(part)
        except Exception:
            pass

    return total


def unique_values(rows: list[dict], field: str) -> list[str]:
    values = sorted({str(row.get(field, "")).strip() for row in rows if str(row.get(field, "")).strip()})
    return values


def main() -> int:
    text = fetch_ena_metadata()

    TSV_PATH.write_text(text, encoding="utf-8")

    rows = parse_tsv(text)

    total_fastq_bytes = sum(safe_int(row.get("fastq_bytes", "")) for row in rows)
    total_submitted_bytes = sum(safe_int(row.get("submitted_bytes", "")) for row in rows)
    total_sra_bytes = sum(safe_int(row.get("sra_bytes", "")) for row in rows)

    summary = {
        "project": PROJECT,
        "rows": len(rows),
        "tsv_path": str(TSV_PATH),
        "unique_studies": unique_values(rows, "study_accession"),
        "unique_secondary_studies": unique_values(rows, "secondary_study_accession"),
        "scientific_names": unique_values(rows, "scientific_name"),
        "tax_ids": unique_values(rows, "tax_id"),
        "library_strategies": unique_values(rows, "library_strategy"),
        "library_sources": unique_values(rows, "library_source"),
        "library_selections": unique_values(rows, "library_selection"),
        "library_layouts": unique_values(rows, "library_layout"),
        "instrument_platforms": unique_values(rows, "instrument_platform"),
        "instrument_models": unique_values(rows, "instrument_model"),
        "n_samples": len(unique_values(rows, "sample_accession")),
        "n_runs": len(unique_values(rows, "run_accession")),
        "total_fastq_bytes": total_fastq_bytes,
        "total_fastq_gb": round(total_fastq_bytes / 1024 / 1024 / 1024, 3),
        "total_submitted_bytes": total_submitted_bytes,
        "total_submitted_gb": round(total_submitted_bytes / 1024 / 1024 / 1024, 3),
        "total_sra_bytes": total_sra_bytes,
        "total_sra_gb": round(total_sra_bytes / 1024 / 1024 / 1024, 3),
    }

    JSON_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append(f"ENA METADATA SUMMARY - {PROJECT}")
    lines.append("=" * 80)
    lines.append(f"Rows: {summary['rows']}")
    lines.append(f"Runs: {summary['n_runs']}")
    lines.append(f"Samples: {summary['n_samples']}")
    lines.append(f"Scientific names: {summary['scientific_names']}")
    lines.append(f"Tax IDs: {summary['tax_ids']}")
    lines.append(f"Library strategies: {summary['library_strategies']}")
    lines.append(f"Library sources: {summary['library_sources']}")
    lines.append(f"Library selections: {summary['library_selections']}")
    lines.append(f"Library layouts: {summary['library_layouts']}")
    lines.append(f"Instrument platforms: {summary['instrument_platforms']}")
    lines.append(f"Instrument models: {summary['instrument_models']}")
    lines.append(f"Total FASTQ GB: {summary['total_fastq_gb']}")
    lines.append(f"Total submitted GB: {summary['total_submitted_gb']}")
    lines.append(f"Total SRA GB: {summary['total_sra_gb']}")
    lines.append("")
    lines.append("Output files:")
    lines.append(f"- {TSV_PATH}")
    lines.append(f"- {JSON_PATH}")
    lines.append(f"- {TXT_PATH}")
    lines.append("")
    lines.append("Prime righe run:")
    lines.append("-" * 80)

    for row in rows[:10]:
        lines.append(
            f"{row.get('run_accession', '')} | "
            f"{row.get('sample_title', '')} | "
            f"{row.get('library_strategy', '')} | "
            f"{row.get('library_layout', '')} | "
            f"{row.get('fastq_bytes', '')}"
        )

    TXT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
