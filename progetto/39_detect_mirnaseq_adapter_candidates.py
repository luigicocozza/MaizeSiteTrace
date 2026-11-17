#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
from collections import Counter
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

FASTQ_DIR = BASE_DIR / "data/public_ena" / PROJECT / "fastq" / "mirnaseq"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

OUT_TXT = OUT_DIR / f"{PROJECT}_mirnaseq_adapter_candidate_report.txt"
OUT_CSV = OUT_DIR / f"{PROJECT}_mirnaseq_adapter_candidate_kmers.csv"

MAX_READS_PER_FILE = 100000
KMER_MIN = 10
KMER_MAX = 24
START_POS_MIN = 18


def read_fastq_sequences(path: Path, max_reads: int):
    n = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        while n < max_reads:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip().upper()
            plus = f.readline()
            qual = f.readline()
            if not seq:
                break
            yield seq
            n += 1


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    global_counter = Counter()
    per_file_rows = []

    fastqs = sorted(FASTQ_DIR.glob("*.fastq.gz"))

    print(f"ADAPTER CANDIDATE DETECTION - {PROJECT}")
    print("=" * 80)
    print(f"FASTQ: {len(fastqs)}")
    print(f"Reads campionate per file: {MAX_READS_PER_FILE}")
    print("")

    for path in fastqs:
        print(f"Analizzo {path.name} ...")
        local_counter = Counter()
        reads = 0

        for seq in read_fastq_sequences(path, MAX_READS_PER_FILE):
            reads += 1

            for k in range(KMER_MIN, KMER_MAX + 1):
                for start in range(START_POS_MIN, max(START_POS_MIN, len(seq) - k + 1)):
                    kmer = seq[start:start + k]
                    if len(kmer) == k and "N" not in kmer:
                        local_counter[kmer] += 1
                        global_counter[kmer] += 1

        for kmer, count in local_counter.most_common(30):
            per_file_rows.append(
                {
                    "file": path.name,
                    "kmer": kmer,
                    "length": len(kmer),
                    "count": count,
                    "sampled_reads": reads,
                    "percent_reads_approx": round(100.0 * count / max(reads, 1), 3),
                }
            )

    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        fieldnames = ["file", "kmer", "length", "count", "sampled_reads", "percent_reads_approx"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(per_file_rows)

    lines = []
    lines.append(f"ADAPTER CANDIDATE REPORT - {PROJECT} miRNA/small RNA")
    lines.append("=" * 80)
    lines.append(f"FASTQ analizzati: {len(fastqs)}")
    lines.append(f"Reads campionate per file: {MAX_READS_PER_FILE}")
    lines.append("")
    lines.append("Top candidati globali:")
    lines.append("rank | kmer | length | count")

    for i, (kmer, count) in enumerate(global_counter.most_common(40), start=1):
        lines.append(f"{i} | {kmer} | {len(kmer)} | {count}")

    lines.append("")
    lines.append("Nota:")
    lines.append("- I candidati più frequenti dopo la posizione 18 possono indicare l'adattatore 3'.")
    lines.append("- Il candidato va confermato prima del trimming definitivo.")
    lines.append("- Dopo il trimming ci aspettiamo molte reads nella finestra 18-30 nt.")

    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("\n".join(lines))
    print("")
    print("Output creati:")
    print(OUT_TXT)
    print(OUT_CSV)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
