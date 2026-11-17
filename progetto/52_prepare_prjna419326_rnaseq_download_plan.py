#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()
OUT_DIR = BASE_DIR / "outputs/prjna419326"

CANDIDATE_INPUTS = [
    OUT_DIR / f"{PROJECT}_download_plan.tsv",
    OUT_DIR / f"{PROJECT}_download_plan_curated.tsv",
    OUT_DIR / f"{PROJECT}_runs_rnaseq.csv",
    OUT_DIR / f"{PROJECT}_runs_all.csv",
    BASE_DIR / "data/public_ena/PRJNA419326/metadata/PRJNA419326_ena_run_metadata.tsv",
]

OUT_TSV = OUT_DIR / f"{PROJECT}_download_plan_rnaseq_only.tsv"
SUMMARY_TXT = OUT_DIR / f"{PROJECT}_download_plan_rnaseq_only_summary.txt"
REPORT_JSON = OUT_DIR / f"{PROJECT}_download_plan_rnaseq_only_report.json"


def read_table(path: Path) -> tuple[list[dict], str]:
    if path.suffix.lower() == ".csv":
        with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
            return list(csv.DictReader(f)), ","
    else:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t")), "\t"


def looks_rnaseq(row: dict) -> bool:
    joined = " ".join(str(v) for v in row.values()).lower()

    if "mirna-seq" in joined:
        return False
    if "mirnaseq" in joined:
        return False
    if "small_rna" in joined:
        return False
    if "small rna" in joined:
        return False
    if "_srna" in joined:
        return False

    return (
        "rna-seq" in joined
        or "rnaseq" in joined
        or "transcriptomic" in joined
        or "mrna" in joined
        or "_rna" in joined
    )


def parse_size_bytes(row: dict) -> int:
    for key in ["fastq_bytes", "size_bytes", "bytes", "read_bytes", "total_bytes"]:
        if key in row and str(row[key]).strip():
            try:
                return int(float(str(row[key]).replace(",", "")))
            except Exception:
                pass

    for key in ["fastq_ftp", "submitted_ftp"]:
        if key in row and str(row[key]).strip():
            # ENA spesso separa più file con ; ma non sempre dà dimensione qui.
            pass

    return 0


def get_run_accession(row: dict) -> str:
    for key in ["run_accession", "run", "Run", "run_id", "accession"]:
        if key in row and str(row[key]).strip():
            return str(row[key]).strip()
    return ""


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    chosen = None
    rows = []
    delimiter = "\t"

    for p in CANDIDATE_INPUTS:
        if not p.exists():
            continue

        tmp, delim = read_table(p)
        if not tmp:
            continue

        rnaseq_rows = [r for r in tmp if looks_rnaseq(r)]

        if rnaseq_rows:
            chosen = p
            rows = rnaseq_rows
            delimiter = delim
            break

    if chosen is None:
        raise SystemExit("ERRORE: non trovo righe RNA-Seq nei file metadata/design disponibili.")

    # dedup per run accession
    dedup = {}
    for r in rows:
        run = get_run_accession(r)
        if run:
            dedup[run] = r

    rows = list(dedup.values()) if dedup else rows

    # ordinamento per run accession se presente
    rows.sort(key=get_run_accession)

    # scrittura tsv conservando tutte le colonne disponibili
    fieldnames = list(rows[0].keys())

    with OUT_TSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    total_bytes = sum(parse_size_bytes(r) for r in rows)
    total_gb = total_bytes / 1e9 if total_bytes else 0

    runs = [get_run_accession(r) for r in rows]

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "input_used": str(chosen),
        "output_tsv": str(OUT_TSV),
        "n_rows": len(rows),
        "runs": runs,
        "estimated_total_bytes": total_bytes,
        "estimated_total_gb": total_gb,
        "columns": fieldnames,
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("PIANO DOWNLOAD RNA-SEQ/mRNA PRJNA419326")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Input usato: {chosen}")
    lines.append(f"Righe RNA-Seq: {len(rows)}")
    lines.append(f"Stima dimensione totale GB: {total_gb:.3f}" if total_gb else "Stima dimensione totale GB: non disponibile nel file")
    lines.append("")
    lines.append("Run RNA-Seq:")
    for r in runs:
        lines.append(f"- {r}")
    lines.append("")
    lines.append("Colonne disponibili:")
    lines.append(", ".join(fieldnames))
    lines.append("")
    lines.append("Output:")
    lines.append(str(OUT_TSV))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
