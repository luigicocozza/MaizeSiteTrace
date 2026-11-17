#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path
from datetime import datetime


BASE_DIR = Path(".").resolve()
PROJECT = "PRJNA419326"

PLAN_TSV = BASE_DIR / "outputs/prjna419326" / f"{PROJECT}_download_plan_mirnaseq_only.tsv"
FASTQ_DIR = BASE_DIR / "data/public_ena" / PROJECT / "fastq" / "mirnaseq"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

REPORT_JSON = OUT_DIR / f"{PROJECT}_mirnaseq_fastq_validation_report.json"
REPORT_TXT = OUT_DIR / f"{PROJECT}_mirnaseq_fastq_validation_report.txt"


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_ok(path: Path) -> bool:
    proc = subprocess.run(
        ["gzip", "-t", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return proc.returncode == 0


def split_field(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(";") if x.strip()]


def main() -> int:
    rows = []

    with PLAN_TSV.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    results = []

    for row in rows:
        urls = split_field(row.get("fastq_ftp", ""))
        md5s = split_field(row.get("fastq_md5", ""))

        for i, url in enumerate(urls):
            filename = Path(url).name
            path = FASTQ_DIR / filename
            aria2_path = Path(str(path) + ".aria2")
            expected_md5 = md5s[i] if i < len(md5s) else ""

            exists = path.exists()
            aria2_exists = aria2_path.exists()
            size_bytes = path.stat().st_size if exists else 0
            actual_md5 = ""
            gzip_status = False

            if exists:
                print(f"Controllo MD5: {filename}")
                actual_md5 = md5_file(path)

                if actual_md5 == expected_md5:
                    print(f"Controllo gzip: {filename}")
                    gzip_status = gzip_ok(path)

            if not exists:
                status = "missing"
            elif aria2_exists:
                status = "has_aria2_incomplete_or_not_finalized"
            elif actual_md5 == expected_md5 and gzip_status:
                status = "valid_md5_gzip_ok"
            elif actual_md5 == expected_md5 and not gzip_status:
                status = "md5_ok_but_gzip_failed"
            else:
                status = "invalid_md5"

            results.append(
                {
                    "run_accession": row.get("run_accession", ""),
                    "sample_alias": row.get("sample_alias", ""),
                    "filename": filename,
                    "path": str(path),
                    "exists": exists,
                    "aria2_exists": aria2_exists,
                    "size_bytes": size_bytes,
                    "expected_md5": expected_md5,
                    "actual_md5": actual_md5,
                    "gzip_ok": gzip_status,
                    "status": status,
                    "fastq_ftp": url,
                }
            )

    counts = {}
    for r in results:
        counts[r["status"]] = counts.get(r["status"], 0) + 1

    report = {
        "project": PROJECT,
        "created": datetime.now().isoformat(timespec="seconds"),
        "fastq_dir": str(FASTQ_DIR),
        "n_files": len(results),
        "status_counts": counts,
        "results": results,
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append(f"VALIDAZIONE FASTQ miRNA - {PROJECT}")
    lines.append("=" * 80)
    lines.append(f"Cartella: {FASTQ_DIR}")
    lines.append(f"File controllati: {len(results)}")
    lines.append("")
    lines.append("Status counts:")
    for k, v in sorted(counts.items()):
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("Dettaglio:")
    for r in results:
        size_gb = r["size_bytes"] / 1024 / 1024 / 1024
        lines.append(
            f"{r['run_accession']} | {r['sample_alias']} | {r['filename']} | "
            f"{size_gb:.3f} GB | aria2={r['aria2_exists']} | "
            f"gzip={r['gzip_ok']} | {r['status']}"
        )

    REPORT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    print("")
    print(f"Report JSON: {REPORT_JSON}")
    print(f"Report TXT: {REPORT_TXT}")

    invalid = [r for r in results if r["status"] != "valid_md5_gzip_ok"]

    if invalid:
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
