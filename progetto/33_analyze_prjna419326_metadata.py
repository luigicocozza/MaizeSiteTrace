#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path
from collections import defaultdict


BASE_DIR = Path(".").resolve()
PROJECT = "PRJNA419326"

META_DIR = BASE_DIR / "data/public_ena" / PROJECT / "metadata"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

TSV_PATH = META_DIR / f"{PROJECT}_ena_run_metadata.tsv"

RUNS_ALL_CSV = OUT_DIR / f"{PROJECT}_runs_all.csv"
RUNS_RNASEQ_CSV = OUT_DIR / f"{PROJECT}_runs_rnaseq.csv"
RUNS_MIRNASEQ_CSV = OUT_DIR / f"{PROJECT}_runs_mirnaseq.csv"
DOWNLOAD_PLAN_TSV = OUT_DIR / f"{PROJECT}_download_plan.tsv"
SUMMARY_JSON = OUT_DIR / f"{PROJECT}_metadata_analysis.json"
SUMMARY_TXT = OUT_DIR / f"{PROJECT}_metadata_analysis_summary.txt"


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


def gb(value: int) -> float:
    return round(value / 1024 / 1024 / 1024, 3)


def read_rows() -> list[dict]:
    with TSV_PATH.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def normalize_row(row: dict) -> dict:
    fastq_bytes = safe_int(row.get("fastq_bytes", ""))
    read_count = safe_int(row.get("read_count", ""))
    base_count = safe_int(row.get("base_count", ""))

    fastq_ftp = row.get("fastq_ftp", "")
    fastq_files = [x for x in fastq_ftp.split(";") if x.strip()]

    text_for_condition = " ".join(
        [
            row.get("sample_alias", ""),
            row.get("sample_title", ""),
            row.get("sample_description", ""),
            row.get("experiment_title", ""),
            row.get("study_title", ""),
        ]
    ).lower()

    possible_condition = "unknown"

    drought_terms = [
        "drought",
        "water deficit",
        "dehydration",
        "dry",
        "stress",
        "peg",
    ]

    control_terms = [
        "control",
        "well-watered",
        "well watered",
        "ww",
        "mock",
        "normal",
    ]

    if any(term in text_for_condition for term in drought_terms):
        possible_condition = "drought_or_stress"

    if any(term in text_for_condition for term in control_terms):
        if possible_condition == "drought_or_stress":
            possible_condition = "mixed_or_unclear"
        else:
            possible_condition = "control_like"

    return {
        "run_accession": row.get("run_accession", ""),
        "sample_accession": row.get("sample_accession", ""),
        "secondary_sample_accession": row.get("secondary_sample_accession", ""),
        "sample_alias": row.get("sample_alias", ""),
        "sample_title": row.get("sample_title", ""),
        "sample_description": row.get("sample_description", ""),
        "experiment_accession": row.get("experiment_accession", ""),
        "experiment_title": row.get("experiment_title", ""),
        "study_accession": row.get("study_accession", ""),
        "secondary_study_accession": row.get("secondary_study_accession", ""),
        "study_title": row.get("study_title", ""),
        "scientific_name": row.get("scientific_name", ""),
        "tax_id": row.get("tax_id", ""),
        "library_strategy": row.get("library_strategy", ""),
        "library_source": row.get("library_source", ""),
        "library_selection": row.get("library_selection", ""),
        "library_layout": row.get("library_layout", ""),
        "instrument_model": row.get("instrument_model", ""),
        "read_count": read_count,
        "base_count": base_count,
        "fastq_bytes": fastq_bytes,
        "fastq_gb": gb(fastq_bytes),
        "n_fastq_files": len(fastq_files),
        "fastq_ftp": fastq_ftp,
        "fastq_md5": row.get("fastq_md5", ""),
        "possible_condition_from_text": possible_condition,
        "first_public": row.get("first_public", ""),
        "last_updated": row.get("last_updated", ""),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fields = list(rows[0].keys())

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows_raw = read_rows()
    rows = [normalize_row(row) for row in rows_raw]

    rnaseq = [row for row in rows if row["library_strategy"] == "RNA-Seq"]
    mirnaseq = [row for row in rows if row["library_strategy"] == "miRNA-Seq"]

    write_csv(RUNS_ALL_CSV, rows)
    write_csv(RUNS_RNASEQ_CSV, rnaseq)
    write_csv(RUNS_MIRNASEQ_CSV, mirnaseq)

    by_strategy = defaultdict(lambda: {"runs": 0, "fastq_bytes": 0})
    by_layout = defaultdict(lambda: {"runs": 0, "fastq_bytes": 0})
    by_condition = defaultdict(lambda: {"runs": 0, "fastq_bytes": 0})

    for row in rows:
        s = row["library_strategy"]
        l = row["library_layout"]
        c = row["possible_condition_from_text"]

        by_strategy[s]["runs"] += 1
        by_strategy[s]["fastq_bytes"] += row["fastq_bytes"]

        by_layout[l]["runs"] += 1
        by_layout[l]["fastq_bytes"] += row["fastq_bytes"]

        by_condition[c]["runs"] += 1
        by_condition[c]["fastq_bytes"] += row["fastq_bytes"]

    download_rows = []

    for row in rows:
        priority = "low"

        if row["library_strategy"] == "miRNA-Seq":
            priority = "high"

        if row["library_strategy"] == "RNA-Seq":
            priority = "medium"

        download_rows.append(
            {
                "priority": priority,
                "run_accession": row["run_accession"],
                "library_strategy": row["library_strategy"],
                "library_layout": row["library_layout"],
                "fastq_gb": row["fastq_gb"],
                "n_fastq_files": row["n_fastq_files"],
                "possible_condition_from_text": row["possible_condition_from_text"],
                "sample_accession": row["sample_accession"],
                "sample_alias": row["sample_alias"],
                "experiment_title": row["experiment_title"],
                "fastq_ftp": row["fastq_ftp"],
                "fastq_md5": row["fastq_md5"],
            }
        )

    with DOWNLOAD_PLAN_TSV.open("w", encoding="utf-8", newline="") as f:
        fields = list(download_rows[0].keys())
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(download_rows)

    summary = {
        "project": PROJECT,
        "total_runs": len(rows),
        "rnaseq_runs": len(rnaseq),
        "mirnaseq_runs": len(mirnaseq),
        "total_fastq_gb": gb(sum(row["fastq_bytes"] for row in rows)),
        "rnaseq_fastq_gb": gb(sum(row["fastq_bytes"] for row in rnaseq)),
        "mirnaseq_fastq_gb": gb(sum(row["fastq_bytes"] for row in mirnaseq)),
        "by_strategy": {
            key: {
                "runs": value["runs"],
                "fastq_gb": gb(value["fastq_bytes"]),
            }
            for key, value in by_strategy.items()
        },
        "by_layout": {
            key: {
                "runs": value["runs"],
                "fastq_gb": gb(value["fastq_bytes"]),
            }
            for key, value in by_layout.items()
        },
        "by_condition_from_text": {
            key: {
                "runs": value["runs"],
                "fastq_gb": gb(value["fastq_bytes"]),
            }
            for key, value in by_condition.items()
        },
        "outputs": {
            "runs_all_csv": str(RUNS_ALL_CSV),
            "rnaseq_csv": str(RUNS_RNASEQ_CSV),
            "mirnaseq_csv": str(RUNS_MIRNASEQ_CSV),
            "download_plan_tsv": str(DOWNLOAD_PLAN_TSV),
            "summary_json": str(SUMMARY_JSON),
            "summary_txt": str(SUMMARY_TXT),
        },
    }

    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append(f"ANALISI METADATI - {PROJECT}")
    lines.append("=" * 80)
    lines.append(f"Run totali: {summary['total_runs']}")
    lines.append(f"RNA-Seq runs: {summary['rnaseq_runs']}")
    lines.append(f"miRNA-Seq runs: {summary['mirnaseq_runs']}")
    lines.append(f"FASTQ totali GB: {summary['total_fastq_gb']}")
    lines.append(f"RNA-Seq FASTQ GB: {summary['rnaseq_fastq_gb']}")
    lines.append(f"miRNA-Seq FASTQ GB: {summary['mirnaseq_fastq_gb']}")
    lines.append("")
    lines.append("Distribuzione per strategia:")
    for key, value in summary["by_strategy"].items():
        lines.append(f"- {key}: {value['runs']} run, {value['fastq_gb']} GB")
    lines.append("")
    lines.append("Distribuzione per layout:")
    for key, value in summary["by_layout"].items():
        lines.append(f"- {key}: {value['runs']} run, {value['fastq_gb']} GB")
    lines.append("")
    lines.append("Condizione ricavata automaticamente dal testo:")
    for key, value in summary["by_condition_from_text"].items():
        lines.append(f"- {key}: {value['runs']} run, {value['fastq_gb']} GB")
    lines.append("")
    lines.append("Prime righe RNA-Seq:")
    lines.append("-" * 80)
    for row in rnaseq[:10]:
        lines.append(
            f"{row['run_accession']} | {row['library_layout']} | {row['fastq_gb']} GB | "
            f"{row['possible_condition_from_text']} | {row['sample_alias']} | {row['experiment_title']}"
        )
    lines.append("")
    lines.append("Prime righe miRNA-Seq:")
    lines.append("-" * 80)
    for row in mirnaseq[:10]:
        lines.append(
            f"{row['run_accession']} | {row['library_layout']} | {row['fastq_gb']} GB | "
            f"{row['possible_condition_from_text']} | {row['sample_alias']} | {row['experiment_title']}"
        )
    lines.append("")
    lines.append("Output generati:")
    for path in summary["outputs"].values():
        lines.append(f"- {path}")

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
