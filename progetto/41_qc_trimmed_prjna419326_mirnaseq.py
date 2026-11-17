#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
import statistics
from collections import Counter
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

PLAN = BASE_DIR / "outputs/prjna419326" / f"{PROJECT}_download_plan_mirnaseq_only.tsv"
TRIM_DIR = BASE_DIR / "data/public_ena" / PROJECT / "processed" / "mirnaseq_trimmed"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

SUMMARY_TXT = OUT_DIR / f"{PROJECT}_mirnaseq_trimmed_qc_summary.txt"
SUMMARY_CSV = OUT_DIR / f"{PROJECT}_mirnaseq_trimmed_qc_summary.csv"
LENGTH_CSV = OUT_DIR / f"{PROJECT}_mirnaseq_trimmed_length_distribution.csv"
REPORT_JSON = OUT_DIR / f"{PROJECT}_mirnaseq_trimmed_qc_report.json"

MAX_READS_PER_FILE = 100000


def read_plan() -> dict[str, dict]:
    with PLAN.open("r", encoding="utf-8") as f:
        return {row["run_accession"]: row for row in csv.DictReader(f, delimiter="\t")}


def inspect_fastq(path: Path, max_reads: int) -> dict:
    lengths = []
    length_counter = Counter()
    gc_values = []
    q_values = []
    reads = 0

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        while reads < max_reads:
            h = f.readline()
            if not h:
                break

            seq = f.readline().strip().upper()
            plus = f.readline()
            qual = f.readline().strip()

            if not seq or not plus or not qual:
                break

            L = len(seq)
            lengths.append(L)
            length_counter[L] += 1

            gc = seq.count("G") + seq.count("C")
            gc_values.append(100.0 * gc / L if L else 0.0)
            q_values.append(sum(max(0, ord(c) - 33) for c in qual) / len(qual) if qual else 0.0)

            reads += 1

    if not lengths:
        return {
            "sampled_reads": 0,
            "min_len": 0,
            "max_len": 0,
            "mean_len": 0,
            "median_len": 0,
            "pct_18_30": 0,
            "pct_20_24": 0,
            "pct_21_24": 0,
            "mean_gc": 0,
            "mean_phred33": 0,
            "top_lengths": [],
            "length_counter": {},
        }

    return {
        "sampled_reads": reads,
        "min_len": min(lengths),
        "max_len": max(lengths),
        "mean_len": round(statistics.mean(lengths), 3),
        "median_len": round(statistics.median(lengths), 3),
        "pct_18_30": round(100.0 * sum(1 for x in lengths if 18 <= x <= 30) / len(lengths), 3),
        "pct_20_24": round(100.0 * sum(1 for x in lengths if 20 <= x <= 24) / len(lengths), 3),
        "pct_21_24": round(100.0 * sum(1 for x in lengths if 21 <= x <= 24) / len(lengths), 3),
        "mean_gc": round(statistics.mean(gc_values), 3),
        "mean_phred33": round(statistics.mean(q_values), 3),
        "top_lengths": length_counter.most_common(15),
        "length_counter": dict(sorted(length_counter.items())),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    plan = read_plan()
    files = sorted(TRIM_DIR.glob("*.trimmed.fastq.gz"))

    results = []
    length_rows = []

    print(f"QC POST-TRIMMING - {PROJECT} miRNA/small RNA")
    print("=" * 80)
    print(f"FASTQ trimmed trovati: {len(files)}")
    print(f"Reads campionate per file: {MAX_READS_PER_FILE}")
    print("")

    for path in files:
        run = path.name.replace(".trimmed.fastq.gz", "")
        meta = plan.get(run, {})

        print(f"Analizzo {path.name} ...")

        qc = inspect_fastq(path, MAX_READS_PER_FILE)

        row = {
            "run_accession": run,
            "sample_alias": meta.get("sample_alias", ""),
            "genotype": meta.get("genotype", ""),
            "condition": meta.get("condition", ""),
            "replicate": meta.get("replicate", ""),
            "trimmed_fastq": str(path),
            "trimmed_size_mb": round(path.stat().st_size / 1024 / 1024, 3),
            **{k: v for k, v in qc.items() if k != "length_counter"},
        }

        results.append(row)

        for length, count in qc["length_counter"].items():
            length_rows.append(
                {
                    "run_accession": run,
                    "sample_alias": meta.get("sample_alias", ""),
                    "genotype": meta.get("genotype", ""),
                    "condition": meta.get("condition", ""),
                    "replicate": meta.get("replicate", ""),
                    "read_length": length,
                    "count": count,
                }
            )

    with SUMMARY_CSV.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "run_accession",
            "sample_alias",
            "genotype",
            "condition",
            "replicate",
            "trimmed_fastq",
            "trimmed_size_mb",
            "sampled_reads",
            "min_len",
            "max_len",
            "mean_len",
            "median_len",
            "pct_18_30",
            "pct_20_24",
            "pct_21_24",
            "mean_gc",
            "mean_phred33",
            "top_lengths",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    with LENGTH_CSV.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "run_accession",
            "sample_alias",
            "genotype",
            "condition",
            "replicate",
            "read_length",
            "count",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(length_rows)

    REPORT_JSON.write_text(
        json.dumps(
            {
                "project": PROJECT,
                "created": datetime.now().isoformat(timespec="seconds"),
                "trim_dir": str(TRIM_DIR),
                "max_reads_per_file": MAX_READS_PER_FILE,
                "n_trimmed_fastq": len(files),
                "results": results,
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    lines = []
    lines.append(f"QC POST-TRIMMING - {PROJECT} miRNA/small RNA")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"FASTQ trimmed analizzati: {len(results)}")
    lines.append(f"Reads campionate per file: {MAX_READS_PER_FILE}")
    lines.append("")
    lines.append("Sintesi:")
    lines.append(
        "run | sample | genotype | condition | replicate | reads | "
        "mean_len | median_len | pct_18_30 | pct_20_24 | pct_21_24 | mean_Q"
    )

    for r in results:
        lines.append(
            f"{r['run_accession']} | {r['sample_alias']} | {r['genotype']} | "
            f"{r['condition']} | {r['replicate']} | {r['sampled_reads']} | "
            f"{r['mean_len']} | {r['median_len']} | {r['pct_18_30']} | "
            f"{r['pct_20_24']} | {r['pct_21_24']} | {r['mean_phred33']}"
        )

    lines.append("")
    lines.append("Controllo biologico:")
    lines.append("- pct_18_30 deve essere vicino al 100% perché abbiamo filtrato le reads nella finestra small RNA.")
    lines.append("- pct_20_24 e pct_21_24 indicano la quota più direttamente compatibile con miRNA/siRNA canonici.")
    lines.append("- Campioni con percentuali molto diverse dagli altri andranno controllati prima dell'analisi differenziale.")

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("\n".join(lines))
    print("")
    print("Output creati:")
    print(SUMMARY_TXT)
    print(SUMMARY_CSV)
    print(LENGTH_CSV)
    print(REPORT_JSON)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
