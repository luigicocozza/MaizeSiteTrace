#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from collections import defaultdict


BASE_DIR = Path(".").resolve()
PROJECT = "PRJNA419326"

IN_CSV = BASE_DIR / "outputs/prjna419326" / f"{PROJECT}_runs_all.csv"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

DESIGN_CSV = OUT_DIR / f"{PROJECT}_experimental_design_curated.csv"
DESIGN_TSV = OUT_DIR / f"{PROJECT}_experimental_design_curated.tsv"
DOWNLOAD_PLAN_TSV = OUT_DIR / f"{PROJECT}_download_plan_curated.tsv"
SUMMARY_JSON = OUT_DIR / f"{PROJECT}_experimental_design_summary.json"
SUMMARY_TXT = OUT_DIR / f"{PROJECT}_experimental_design_summary.txt"


def read_csv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_alias(alias: str, library_strategy: str) -> dict:
    alias = str(alias).strip()

    result = {
        "genotype": "unknown",
        "condition_code": "unknown",
        "condition": "unknown",
        "sequencing_type": "unknown",
        "replicate": "unknown",
        "parse_status": "warning",
    }

    pattern = r"^(?P<genotype>[A-Za-z0-9]+)_(?P<condition>WW|DS)(?:_(?P<srna>sRNA))?_rep(?P<replicate>[0-9]+)$"
    match = re.match(pattern, alias)

    if not match:
        return result

    genotype = match.group("genotype")
    condition_code = match.group("condition")
    srna = match.group("srna")
    replicate = match.group("replicate")

    condition = {
        "WW": "well_watered_control",
        "DS": "drought_stress",
    }.get(condition_code, "unknown")

    if srna == "sRNA" or library_strategy == "miRNA-Seq":
        sequencing_type = "small_rna_mirna_seq"
    elif library_strategy == "RNA-Seq":
        sequencing_type = "mrna_rnaseq"
    else:
        sequencing_type = "unknown"

    result.update(
        {
            "genotype": genotype,
            "condition_code": condition_code,
            "condition": condition,
            "sequencing_type": sequencing_type,
            "replicate": replicate,
            "parse_status": "ok",
        }
    )

    return result


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = read_csv(IN_CSV)

    curated = []

    for row in rows:
        parsed = parse_alias(row.get("sample_alias", ""), row.get("library_strategy", ""))

        fastq_gb = float(row.get("fastq_gb", 0) or 0)

        priority = "low"

        if parsed["sequencing_type"] == "small_rna_mirna_seq":
            priority = "high_first_download"
        elif parsed["sequencing_type"] == "mrna_rnaseq":
            priority = "medium_second_download"

        curated.append(
            {
                "project": PROJECT,
                "run_accession": row.get("run_accession", ""),
                "sample_accession": row.get("sample_accession", ""),
                "sample_alias": row.get("sample_alias", ""),
                "genotype": parsed["genotype"],
                "condition_code": parsed["condition_code"],
                "condition": parsed["condition"],
                "replicate": parsed["replicate"],
                "sequencing_type": parsed["sequencing_type"],
                "library_strategy": row.get("library_strategy", ""),
                "library_layout": row.get("library_layout", ""),
                "n_fastq_files": row.get("n_fastq_files", ""),
                "fastq_gb": fastq_gb,
                "fastq_ftp": row.get("fastq_ftp", ""),
                "fastq_md5": row.get("fastq_md5", ""),
                "experiment_title": row.get("experiment_title", ""),
                "parse_status": parsed["parse_status"],
                "download_priority": priority,
            }
        )

    curated = sorted(
        curated,
        key=lambda r: (
            r["sequencing_type"],
            r["genotype"],
            r["condition_code"],
            int(r["replicate"]) if str(r["replicate"]).isdigit() else 999,
        ),
    )

    write_csv(DESIGN_CSV, curated)
    write_tsv(DESIGN_TSV, curated)
    write_tsv(DOWNLOAD_PLAN_TSV, curated)

    summary_groups = defaultdict(lambda: {"runs": 0, "fastq_gb": 0.0})

    for row in curated:
        key = (
            row["sequencing_type"],
            row["genotype"],
            row["condition"],
        )

        summary_groups[key]["runs"] += 1
        summary_groups[key]["fastq_gb"] += float(row["fastq_gb"])

    summary_table = []

    for key, value in sorted(summary_groups.items()):
        sequencing_type, genotype, condition = key
        summary_table.append(
            {
                "sequencing_type": sequencing_type,
                "genotype": genotype,
                "condition": condition,
                "runs": value["runs"],
                "fastq_gb": round(value["fastq_gb"], 3),
            }
        )

    mirna_rows = [r for r in curated if r["sequencing_type"] == "small_rna_mirna_seq"]
    rnaseq_rows = [r for r in curated if r["sequencing_type"] == "mrna_rnaseq"]

    summary = {
        "project": PROJECT,
        "total_runs": len(curated),
        "parse_ok": sum(1 for r in curated if r["parse_status"] == "ok"),
        "parse_warning": sum(1 for r in curated if r["parse_status"] != "ok"),
        "mirna_runs": len(mirna_rows),
        "rnaseq_runs": len(rnaseq_rows),
        "mirna_fastq_gb": round(sum(float(r["fastq_gb"]) for r in mirna_rows), 3),
        "rnaseq_fastq_gb": round(sum(float(r["fastq_gb"]) for r in rnaseq_rows), 3),
        "total_fastq_gb": round(sum(float(r["fastq_gb"]) for r in curated), 3),
        "groups": summary_table,
        "recommended_strategy": (
            "Scaricare prima i 12 run miRNA-Seq/small RNA-Seq, circa 14.271 GB, "
            "per integrare direttamente gli ncRNA/miRNA nel ranking V3. "
            "Scaricare i 12 run RNA-Seq, circa 47.128 GB, in una fase successiva "
            "per integrare espressione dei geni target."
        ),
        "outputs": {
            "design_csv": str(DESIGN_CSV),
            "design_tsv": str(DESIGN_TSV),
            "download_plan_tsv": str(DOWNLOAD_PLAN_TSV),
            "summary_json": str(SUMMARY_JSON),
            "summary_txt": str(SUMMARY_TXT),
        },
    }

    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append(f"DESIGN SPERIMENTALE CURATO - {PROJECT}")
    lines.append("=" * 80)
    lines.append(f"Run totali: {summary['total_runs']}")
    lines.append(f"Parse OK: {summary['parse_ok']}")
    lines.append(f"Parse warning: {summary['parse_warning']}")
    lines.append(f"miRNA/small RNA-Seq runs: {summary['mirna_runs']} ({summary['mirna_fastq_gb']} GB)")
    lines.append(f"RNA-Seq runs: {summary['rnaseq_runs']} ({summary['rnaseq_fastq_gb']} GB)")
    lines.append(f"FASTQ totali: {summary['total_fastq_gb']} GB")
    lines.append("")
    lines.append("Gruppi sperimentali:")
    lines.append("-" * 80)

    for group in summary_table:
        lines.append(
            f"{group['sequencing_type']} | {group['genotype']} | {group['condition']} | "
            f"{group['runs']} run | {group['fastq_gb']} GB"
        )

    lines.append("")
    lines.append("Strategia consigliata:")
    lines.append(summary["recommended_strategy"])
    lines.append("")
    lines.append("Prime righe design:")
    lines.append("-" * 80)

    for row in curated[:12]:
        lines.append(
            f"{row['run_accession']} | {row['sample_alias']} | {row['genotype']} | "
            f"{row['condition']} | rep{row['replicate']} | {row['sequencing_type']} | "
            f"{row['fastq_gb']} GB"
        )

    lines.append("")
    lines.append("Output generati:")
    for value in summary["outputs"].values():
        lines.append(f"- {value}")

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
