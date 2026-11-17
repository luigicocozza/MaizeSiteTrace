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

FASTQ_DIR = BASE_DIR / "data/public_ena" / PROJECT / "fastq" / "mirnaseq"
PLAN_TSV = BASE_DIR / "outputs/prjna419326" / f"{PROJECT}_download_plan_mirnaseq_only.tsv"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

SUMMARY_TXT = OUT_DIR / f"{PROJECT}_mirnaseq_light_fastq_qc_summary.txt"
SUMMARY_CSV = OUT_DIR / f"{PROJECT}_mirnaseq_light_fastq_qc_summary.csv"
LENGTH_CSV = OUT_DIR / f"{PROJECT}_mirnaseq_light_fastq_length_distribution.csv"
REPORT_JSON = OUT_DIR / f"{PROJECT}_mirnaseq_light_fastq_qc_report.json"

MAX_READS_PER_FILE = 50000


def read_plan() -> dict[str, dict]:
    data = {}

    with PLAN_TSV.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            run = row["run_accession"]
            data[run] = row

    return data


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)


def mean_phred33(qual: str) -> float:
    if not qual:
        return 0.0
    return sum(max(0, ord(c) - 33) for c in qual) / len(qual)


def inspect_fastq(path: Path, max_reads: int) -> dict:
    lengths = []
    gc_values = []
    q_values = []
    length_counter = Counter()

    reads = 0

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        while reads < max_reads:
            header = f.readline()
            if not header:
                break

            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()

            if not seq or not plus or not qual:
                break

            L = len(seq)
            lengths.append(L)
            length_counter[L] += 1
            gc_values.append(gc_percent(seq))
            q_values.append(mean_phred33(qual))
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
            "mean_gc": 0,
            "mean_phred33": 0,
            "top_lengths": [],
            "length_counter": {},
        }

    pct_18_30 = 100.0 * sum(1 for x in lengths if 18 <= x <= 30) / len(lengths)
    pct_20_24 = 100.0 * sum(1 for x in lengths if 20 <= x <= 24) / len(lengths)

    return {
        "sampled_reads": reads,
        "min_len": min(lengths),
        "max_len": max(lengths),
        "mean_len": round(statistics.mean(lengths), 3),
        "median_len": round(statistics.median(lengths), 3),
        "pct_18_30": round(pct_18_30, 3),
        "pct_20_24": round(pct_20_24, 3),
        "mean_gc": round(statistics.mean(gc_values), 3),
        "mean_phred33": round(statistics.mean(q_values), 3),
        "top_lengths": length_counter.most_common(10),
        "length_counter": dict(sorted(length_counter.items())),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    plan = read_plan()
    fastqs = sorted(FASTQ_DIR.glob("*.fastq.gz"))

    results = []
    length_rows = []

    print(f"LIGHT QC FASTQ - {PROJECT}")
    print("=" * 80)
    print(f"FASTQ trovati: {len(fastqs)}")
    print(f"Reads campionate per file: {MAX_READS_PER_FILE}")
    print("")

    for path in fastqs:
        run = path.name.replace(".fastq.gz", "")
        meta = plan.get(run, {})

        print(f"Controllo {path.name} ...")

        qc = inspect_fastq(path, MAX_READS_PER_FILE)

        result = {
            "run_accession": run,
            "sample_alias": meta.get("sample_alias", ""),
            "genotype": meta.get("genotype", ""),
            "condition": meta.get("condition", ""),
            "replicate": meta.get("replicate", ""),
            "file": str(path),
            "size_gb": round(path.stat().st_size / 1024 / 1024 / 1024, 3),
            **{k: v for k, v in qc.items() if k not in {"length_counter"}},
        }

        results.append(result)

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
            "file",
            "size_gb",
            "sampled_reads",
            "min_len",
            "max_len",
            "mean_len",
            "median_len",
            "pct_18_30",
            "pct_20_24",
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

    report = {
        "project": PROJECT,
        "created": datetime.now().isoformat(timespec="seconds"),
        "fastq_dir": str(FASTQ_DIR),
        "max_reads_per_file": MAX_READS_PER_FILE,
        "n_fastq": len(fastqs),
        "results": results,
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append(f"LIGHT QC FASTQ - {PROJECT} miRNA/small RNA")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"FASTQ analizzati: {len(results)}")
    lines.append(f"Reads campionate per file: {MAX_READS_PER_FILE}")
    lines.append("")
    lines.append("Sintesi per campione:")
    lines.append(
        "run_accession | sample_alias | genotype | condition | replicate | "
        "reads | mean_len | median_len | pct_18_30 | pct_20_24 | mean_GC | mean_Q"
    )

    for r in results:
        lines.append(
            f"{r['run_accession']} | {r['sample_alias']} | {r['genotype']} | "
            f"{r['condition']} | {r['replicate']} | {r['sampled_reads']} | "
            f"{r['mean_len']} | {r['median_len']} | {r['pct_18_30']} | "
            f"{r['pct_20_24']} | {r['mean_gc']} | {r['mean_phred33']}"
        )

    lines.append("")
    lines.append("Interpretazione:")
    lines.append("- pct_18_30 indica la quota di reads già nella finestra tipica small RNA/miRNA.")
    lines.append("- Se le reads risultano molto più lunghe, è normale nei dati grezzi: servirà trimming degli adattatori.")
    lines.append("- Questo controllo è leggero: non sostituisce FastQC/MultiQC, ma verifica rapidamente struttura e coerenza dei FASTQ.")

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
