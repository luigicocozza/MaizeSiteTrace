#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

PLAN = BASE_DIR / "outputs/prjna419326" / f"{PROJECT}_download_plan_mirnaseq_only.tsv"
TRIM_DIR = BASE_DIR / "data/public_ena" / PROJECT / "processed" / "mirnaseq_trimmed"

OUT_DIR = BASE_DIR / "outputs/prjna419326"
COLLAPSED_DIR = OUT_DIR / "collapsed_srna_counts"

SAMPLE_SUMMARY_CSV = OUT_DIR / f"{PROJECT}_mirnaseq_collapsed_sample_summary.csv"
COUNT_MATRIX_TSV = OUT_DIR / f"{PROJECT}_mirnaseq_srna_count_matrix_filtered.tsv"
RPM_MATRIX_TSV = OUT_DIR / f"{PROJECT}_mirnaseq_srna_rpm_matrix_filtered.tsv"
TOP_FEATURES_TSV = OUT_DIR / f"{PROJECT}_mirnaseq_top_srna_features.tsv"
REPORT_JSON = OUT_DIR / f"{PROJECT}_mirnaseq_srna_count_matrix_report.json"
SUMMARY_TXT = OUT_DIR / f"{PROJECT}_mirnaseq_srna_count_matrix_summary.txt"

MIN_TOTAL_COUNT = 10
MIN_SAMPLES_DETECTED = 2
TOP_N = 5000


def read_plan():
    with PLAN.open("r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    # Ordine biologico già presente nel piano
    return rows


def iter_fastq_sequences(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip().upper()
            plus = f.readline()
            qual = f.readline()

            if not seq:
                break

            yield seq


def collapse_sample(run: str, path: Path, out_path: Path) -> dict:
    counter = Counter()
    total_reads = 0
    accepted_reads = 0

    for seq in iter_fastq_sequences(path):
        total_reads += 1
        L = len(seq)

        if 18 <= L <= 30 and "N" not in seq:
            counter[seq] += 1
            accepted_reads += 1

    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sequence", "length", "count"])

        for seq, count in counter.most_common():
            writer.writerow([seq, len(seq), count])

    return {
        "run_accession": run,
        "trimmed_fastq": str(path),
        "collapsed_tsv_gz": str(out_path),
        "total_trimmed_reads": total_reads,
        "accepted_18_30_noN_reads": accepted_reads,
        "unique_sequences": len(counter),
        "top_sequence": counter.most_common(1)[0][0] if counter else "",
        "top_sequence_count": counter.most_common(1)[0][1] if counter else 0,
    }


def read_collapsed_counts(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield row["sequence"], int(row["count"])


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    COLLAPSED_DIR.mkdir(parents=True, exist_ok=True)

    plan = read_plan()

    print(f"COLLASSO READS TRIMMED E MATRICE CONTEGGI - {PROJECT}")
    print("=" * 80)
    print(f"Campioni nel piano: {len(plan)}")
    print(f"Filtro matrice: total_count >= {MIN_TOTAL_COUNT}, samples_detected >= {MIN_SAMPLES_DETECTED}")
    print("")

    sample_rows = []
    sample_aliases = []
    run_to_collapsed = {}

    # 1. Collasso per campione
    for i, row in enumerate(plan, start=1):
        run = row["run_accession"]
        sample = row["sample_alias"]

        trimmed = TRIM_DIR / f"{run}.trimmed.fastq.gz"
        collapsed = COLLAPSED_DIR / f"{run}_{sample}.collapsed.tsv.gz"

        print("=" * 80)
        print(f"[{i}/{len(plan)}] Collasso {run} | {sample}")
        print("=" * 80)

        if not trimmed.exists():
            raise FileNotFoundError(f"Manca file trimmed: {trimmed}")

        summary = collapse_sample(run, trimmed, collapsed)

        merged = {
            "sample_alias": sample,
            "genotype": row.get("genotype", ""),
            "condition": row.get("condition", ""),
            "replicate": row.get("replicate", ""),
            **summary,
        }

        sample_rows.append(merged)
        sample_aliases.append(sample)
        run_to_collapsed[run] = collapsed

        print(f"Reads trimmed totali: {summary['total_trimmed_reads']}")
        print(f"Reads accettate 18-30 noN: {summary['accepted_18_30_noN_reads']}")
        print(f"Sequenze uniche: {summary['unique_sequences']}")
        print(f"Output: {collapsed}")

    # 2. Riassunto campioni
    with SAMPLE_SUMMARY_CSV.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "run_accession",
            "sample_alias",
            "genotype",
            "condition",
            "replicate",
            "trimmed_fastq",
            "collapsed_tsv_gz",
            "total_trimmed_reads",
            "accepted_18_30_noN_reads",
            "unique_sequences",
            "top_sequence",
            "top_sequence_count",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sample_rows)

    # 3. Primo passaggio globale: total count e presenza
    total_counts = Counter()
    sample_presence = Counter()

    print("")
    print("Calcolo conteggi globali e presenza nei campioni...")

    for row in sample_rows:
        collapsed = Path(row["collapsed_tsv_gz"])
        seen_in_sample = set()

        for seq, count in read_collapsed_counts(collapsed):
            total_counts[seq] += count
            seen_in_sample.add(seq)

        for seq in seen_in_sample:
            sample_presence[seq] += 1

    selected_sequences = [
        seq for seq, total in total_counts.items()
        if total >= MIN_TOTAL_COUNT and sample_presence[seq] >= MIN_SAMPLES_DETECTED
    ]

    selected_sequences.sort(key=lambda s: (-total_counts[s], s))

    print(f"Sequenze totali globali: {len(total_counts)}")
    print(f"Sequenze filtrate per matrice: {len(selected_sequences)}")

    selected_set = set(selected_sequences)

    # 4. Carico solo i conteggi delle sequenze filtrate
    counts_by_sample = {}

    for row in sample_rows:
        sample = row["sample_alias"]
        collapsed = Path(row["collapsed_tsv_gz"])

        d = {}
        for seq, count in read_collapsed_counts(collapsed):
            if seq in selected_set:
                d[seq] = count

        counts_by_sample[sample] = d

    library_sizes = {
        row["sample_alias"]: int(row["accepted_18_30_noN_reads"])
        for row in sample_rows
    }

    # 5. Matrice count
    with COUNT_MATRIX_TSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            ["feature_id", "sequence", "length", "total_count", "samples_detected"]
            + sample_aliases
        )

        for idx, seq in enumerate(selected_sequences, start=1):
            feature_id = f"sRNA_{idx:08d}"
            row_counts = [counts_by_sample[s].get(seq, 0) for s in sample_aliases]

            writer.writerow(
                [feature_id, seq, len(seq), total_counts[seq], sample_presence[seq]]
                + row_counts
            )

    # 6. Matrice RPM
    with RPM_MATRIX_TSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            ["feature_id", "sequence", "length", "total_count", "samples_detected"]
            + sample_aliases
        )

        for idx, seq in enumerate(selected_sequences, start=1):
            feature_id = f"sRNA_{idx:08d}"
            rpm_values = []

            for sample in sample_aliases:
                count = counts_by_sample[sample].get(seq, 0)
                lib = max(library_sizes[sample], 1)
                rpm = 1_000_000 * count / lib
                rpm_values.append(round(rpm, 6))

            writer.writerow(
                [feature_id, seq, len(seq), total_counts[seq], sample_presence[seq]]
                + rpm_values
            )

    # 7. Top feature
    with TOP_FEATURES_TSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            ["feature_id", "sequence", "length", "total_count", "samples_detected"]
            + sample_aliases
        )

        for idx, seq in enumerate(selected_sequences[:TOP_N], start=1):
            feature_id = f"sRNA_{idx:08d}"
            row_counts = [counts_by_sample[s].get(seq, 0) for s in sample_aliases]

            writer.writerow(
                [feature_id, seq, len(seq), total_counts[seq], sample_presence[seq]]
                + row_counts
            )

    report = {
        "project": PROJECT,
        "created": datetime.now().isoformat(timespec="seconds"),
        "min_total_count": MIN_TOTAL_COUNT,
        "min_samples_detected": MIN_SAMPLES_DETECTED,
        "n_samples": len(sample_rows),
        "sample_aliases": sample_aliases,
        "total_unique_sequences_global": len(total_counts),
        "filtered_sequences_in_matrix": len(selected_sequences),
        "sample_summary": sample_rows,
        "outputs": {
            "sample_summary_csv": str(SAMPLE_SUMMARY_CSV),
            "count_matrix_tsv": str(COUNT_MATRIX_TSV),
            "rpm_matrix_tsv": str(RPM_MATRIX_TSV),
            "top_features_tsv": str(TOP_FEATURES_TSV),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append(f"MATRICE CONTEGGI small RNA - {PROJECT}")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Campioni: {len(sample_rows)}")
    lines.append(f"Sequenze uniche globali: {len(total_counts)}")
    lines.append(f"Sequenze filtrate in matrice: {len(selected_sequences)}")
    lines.append(f"Filtro: total_count >= {MIN_TOTAL_COUNT}, samples_detected >= {MIN_SAMPLES_DETECTED}")
    lines.append("")
    lines.append("Campioni:")
    lines.append("run | sample | genotype | condition | replicate | accepted_reads | unique_sequences")

    for r in sample_rows:
        lines.append(
            f"{r['run_accession']} | {r['sample_alias']} | {r['genotype']} | "
            f"{r['condition']} | {r['replicate']} | {r['accepted_18_30_noN_reads']} | "
            f"{r['unique_sequences']}"
        )

    lines.append("")
    lines.append("Output:")
    lines.append(str(SAMPLE_SUMMARY_CSV))
    lines.append(str(COUNT_MATRIX_TSV))
    lines.append(str(RPM_MATRIX_TSV))
    lines.append(str(TOP_FEATURES_TSV))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
