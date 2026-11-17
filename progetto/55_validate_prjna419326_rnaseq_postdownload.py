#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

PLAN = BASE_DIR / "outputs/prjna419326/PRJNA419326_download_plan_rnaseq_only.tsv"
FASTQ_DIR = BASE_DIR / "data/public_ena/PRJNA419326/rnaseq_fastq"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

SAMPLE_SHEET = OUT_DIR / "PRJNA419326_rnaseq_sample_sheet.tsv"
VALIDATION_CSV = OUT_DIR / "PRJNA419326_rnaseq_postdownload_validation.csv"
SUMMARY_TXT = OUT_DIR / "PRJNA419326_rnaseq_postdownload_validation_summary.txt"
REPORT_JSON = OUT_DIR / "PRJNA419326_rnaseq_postdownload_validation_report.json"


def split_multi(value: str) -> list[str]:
    value = str(value or "").strip()
    if not value:
        return []
    return [x.strip() for x in value.split(";") if x.strip()]


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024 * 8), b""):
            h.update(block)
    return h.hexdigest()


def gzip_test(path: Path) -> bool:
    try:
        completed = subprocess.run(
            ["gzip", "-t", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return completed.returncode == 0
    except Exception:
        return False


def infer_genotype(sample_alias: str) -> str:
    if sample_alias.startswith("B73_"):
        return "B73"
    if sample_alias.startswith("H21_"):
        return "H21"
    return "unknown"


def infer_condition(sample_alias: str) -> str:
    if "_DS_" in sample_alias:
        return "DS"
    if "_WW_" in sample_alias:
        return "WW"
    return "unknown"


def infer_replicate(sample_alias: str) -> str:
    if "_rep" in sample_alias:
        return sample_alias.split("_rep")[-1]
    return ""


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with PLAN.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    validation_rows = []
    sample_rows = []

    for row in rows:
        run = row["run_accession"]
        sample_alias = row.get("sample_alias", "")
        genotype = infer_genotype(sample_alias)
        condition = infer_condition(sample_alias)
        replicate = infer_replicate(sample_alias)

        urls = split_multi(row.get("fastq_ftp", ""))
        md5s = split_multi(row.get("fastq_md5", ""))

        fastq_files = []
        statuses = []

        for idx, url in enumerate(urls):
            fname = Path(url).name
            path = FASTQ_DIR / fname
            expected_md5 = md5s[idx] if idx < len(md5s) else ""

            exists = path.exists()
            size_bytes = path.stat().st_size if exists else 0
            observed_md5 = md5sum(path) if exists else ""
            md5_ok = bool(expected_md5 and observed_md5 == expected_md5)
            gz_ok = gzip_test(path) if exists else False

            status = "OK" if exists and md5_ok and gz_ok else "FAIL"
            statuses.append(status)
            fastq_files.append(str(path))

            validation_rows.append(
                {
                    "run_accession": run,
                    "sample_alias": sample_alias,
                    "genotype": genotype,
                    "condition": condition,
                    "replicate": replicate,
                    "pair": idx + 1,
                    "filename": fname,
                    "path": str(path),
                    "exists": exists,
                    "size_bytes": size_bytes,
                    "expected_md5": expected_md5,
                    "observed_md5": observed_md5,
                    "md5_ok": md5_ok,
                    "gzip_ok": gz_ok,
                    "status": status,
                }
            )

        sample_rows.append(
            {
                "run_accession": run,
                "sample_alias": sample_alias,
                "genotype": genotype,
                "condition": condition,
                "replicate": replicate,
                "library_strategy": row.get("library_strategy", ""),
                "library_layout": row.get("library_layout", ""),
                "fastq_1": fastq_files[0] if len(fastq_files) > 0 else "",
                "fastq_2": fastq_files[1] if len(fastq_files) > 1 else "",
                "sample_status": "OK" if all(s == "OK" for s in statuses) and len(statuses) == 2 else "FAIL",
            }
        )

    with VALIDATION_CSV.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "run_accession",
            "sample_alias",
            "genotype",
            "condition",
            "replicate",
            "pair",
            "filename",
            "path",
            "exists",
            "size_bytes",
            "expected_md5",
            "observed_md5",
            "md5_ok",
            "gzip_ok",
            "status",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(validation_rows)

    with SAMPLE_SHEET.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "run_accession",
            "sample_alias",
            "genotype",
            "condition",
            "replicate",
            "library_strategy",
            "library_layout",
            "fastq_1",
            "fastq_2",
            "sample_status",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sample_rows)

    n_fastq = len(validation_rows)
    n_fastq_ok = sum(1 for r in validation_rows if r["status"] == "OK")
    n_samples = len(sample_rows)
    n_samples_ok = sum(1 for r in sample_rows if r["sample_status"] == "OK")

    group_counts = {}
    for r in sample_rows:
        key = f"{r['genotype']}_{r['condition']}"
        group_counts[key] = group_counts.get(key, 0) + 1

    final_status = "PASS" if n_fastq == n_fastq_ok and n_samples == n_samples_ok else "FAIL"

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "status": final_status,
        "plan": str(PLAN),
        "fastq_dir": str(FASTQ_DIR),
        "n_samples": n_samples,
        "n_samples_ok": n_samples_ok,
        "n_fastq": n_fastq,
        "n_fastq_ok": n_fastq_ok,
        "group_counts": group_counts,
        "outputs": {
            "sample_sheet": str(SAMPLE_SHEET),
            "validation_csv": str(VALIDATION_CSV),
            "summary_txt": str(SUMMARY_TXT),
            "report_json": str(REPORT_JSON),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("VALIDAZIONE POST-DOWNLOAD RNA-SEQ PRJNA419326")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Status: {final_status}")
    lines.append(f"Campioni: {n_samples}")
    lines.append(f"Campioni OK: {n_samples_ok}")
    lines.append(f"FASTQ: {n_fastq}")
    lines.append(f"FASTQ OK: {n_fastq_ok}")
    lines.append("")
    lines.append("Gruppi sperimentali:")
    for k in sorted(group_counts):
        lines.append(f"- {k}: {group_counts[k]}")
    lines.append("")
    lines.append("Output:")
    lines.append(str(SAMPLE_SHEET))
    lines.append(str(VALIDATION_CSV))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0 if final_status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
