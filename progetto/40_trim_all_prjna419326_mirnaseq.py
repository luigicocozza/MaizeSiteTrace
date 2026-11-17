#!/usr/bin/env python3
from __future__ import annotations

import csv
import subprocess
import sys
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

PLAN = BASE_DIR / "outputs/prjna419326" / f"{PROJECT}_download_plan_mirnaseq_only.tsv"
RAW_DIR = BASE_DIR / "data/public_ena" / PROJECT / "fastq" / "mirnaseq"
TRIM_DIR = BASE_DIR / "data/public_ena" / PROJECT / "processed" / "mirnaseq_trimmed"
REPORT_DIR = BASE_DIR / "outputs/prjna419326" / "cutadapt"

SUMMARY_CSV = BASE_DIR / "outputs/prjna419326" / f"{PROJECT}_mirnaseq_cutadapt_trimming_summary.csv"
SUMMARY_TXT = BASE_DIR / "outputs/prjna419326" / f"{PROJECT}_mirnaseq_cutadapt_trimming_summary.txt"

ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


def read_plan():
    with PLAN.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def parse_cutadapt_report(path: Path) -> dict:
    metrics = {
        "total_reads_processed": "",
        "reads_with_adapters": "",
        "reads_too_short": "",
        "reads_too_long": "",
        "reads_written": "",
    }

    if not path.exists():
        return metrics

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()

        if s.startswith("Total reads processed:"):
            metrics["total_reads_processed"] = s.split(":", 1)[1].strip()
        elif s.startswith("Reads with adapters:"):
            metrics["reads_with_adapters"] = s.split(":", 1)[1].strip()
        elif s.startswith("Reads that were too short:"):
            metrics["reads_too_short"] = s.split(":", 1)[1].strip()
        elif s.startswith("Reads that were too long:"):
            metrics["reads_too_long"] = s.split(":", 1)[1].strip()
        elif s.startswith("Reads written"):
            metrics["reads_written"] = s.split(":", 1)[1].strip()

    return metrics


def main() -> int:
    TRIM_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    rows = read_plan()
    results = []

    print(f"TRIMMING COMPLETO - {PROJECT} miRNA/small RNA")
    print("=" * 80)
    print(f"Campioni da processare: {len(rows)}")
    print(f"Adapter 3': {ADAPTER}")
    print("")

    for i, row in enumerate(rows, start=1):
        run = row["run_accession"]
        sample = row.get("sample_alias", "")
        genotype = row.get("genotype", "")
        condition = row.get("condition", "")
        replicate = row.get("replicate", "")

        raw = RAW_DIR / f"{run}.fastq.gz"
        trimmed = TRIM_DIR / f"{run}.trimmed.fastq.gz"
        report = REPORT_DIR / f"{run}_cutadapt_report.txt"

        print("=" * 80)
        print(f"[{i}/{len(rows)}] {run} | {sample} | {genotype} | {condition} | rep {replicate}")
        print("=" * 80)

        if not raw.exists():
            print(f"ERRORE: manca input {raw}")
            status = "missing_input"
        else:
            cmd = [
                sys.executable,
                "-m",
                "cutadapt",
                "-a",
                ADAPTER,
                "-q",
                "20",
                "-m",
                "18",
                "-M",
                "30",
                "--trim-n",
                "-o",
                str(trimmed),
                str(raw),
            ]

            print("Eseguo cutadapt...")
            with report.open("w", encoding="utf-8") as log:
                proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)

            status = "ok" if proc.returncode == 0 else f"error_{proc.returncode}"

            print(f"Stato: {status}")
            print(f"Output: {trimmed}")
            print(f"Report: {report}")

        metrics = parse_cutadapt_report(report)

        results.append(
            {
                "run_accession": run,
                "sample_alias": sample,
                "genotype": genotype,
                "condition": condition,
                "replicate": replicate,
                "status": status,
                "raw_fastq": str(raw),
                "trimmed_fastq": str(trimmed),
                "trimmed_size_mb": round(trimmed.stat().st_size / 1024 / 1024, 3) if trimmed.exists() else "",
                "report": str(report),
                **metrics,
            }
        )

    fieldnames = [
        "run_accession",
        "sample_alias",
        "genotype",
        "condition",
        "replicate",
        "status",
        "raw_fastq",
        "trimmed_fastq",
        "trimmed_size_mb",
        "report",
        "total_reads_processed",
        "reads_with_adapters",
        "reads_too_short",
        "reads_too_long",
        "reads_written",
    ]

    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    ok_count = sum(1 for r in results if r["status"] == "ok")

    lines = []
    lines.append(f"TRIMMING COMPLETO - {PROJECT} miRNA/small RNA")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Campioni processati: {len(results)}")
    lines.append(f"Campioni OK: {ok_count}/{len(results)}")
    lines.append(f"Adapter 3': {ADAPTER}")
    lines.append("")
    lines.append("Sintesi:")
    lines.append("run | sample | genotype | condition | replicate | status | trimmed_MB | reads_written")

    for r in results:
        lines.append(
            f"{r['run_accession']} | {r['sample_alias']} | {r['genotype']} | "
            f"{r['condition']} | {r['replicate']} | {r['status']} | "
            f"{r['trimmed_size_mb']} | {r['reads_written']}"
        )

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("\n".join(lines))
    print("")
    print("Output creati:")
    print(SUMMARY_TXT)
    print(SUMMARY_CSV)

    return 0 if ok_count == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
