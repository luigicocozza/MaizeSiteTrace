#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()
OUT_DIR = BASE_DIR / "outputs/prjna419326"

PLAN = OUT_DIR / f"{PROJECT}_download_plan_rnaseq_only.tsv"
SUMMARY_TXT = OUT_DIR / f"{PROJECT}_rnaseq_download_plan_detailed_summary.txt"
REPORT_JSON = OUT_DIR / f"{PROJECT}_rnaseq_download_plan_detailed_report.json"
ARIA2_INPUT = OUT_DIR / f"{PROJECT}_rnaseq_aria2_urls.txt"
EXPECTED_MD5 = OUT_DIR / f"{PROJECT}_rnaseq_expected_md5.txt"


def split_multi(value: str) -> list[str]:
    if value is None:
        return []
    value = str(value).strip()
    if not value:
        return []
    return [x.strip() for x in value.split(";") if x.strip()]


def https_from_ftp(url: str) -> str:
    url = str(url).strip()
    if url.startswith("ftp.sra.ebi.ac.uk/"):
        return "https://" + url
    if url.startswith("ftp://ftp.sra.ebi.ac.uk/"):
        return "https://" + url.replace("ftp://", "")
    if url.startswith("https://"):
        return url
    if url.startswith("http://"):
        return url
    return url


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with PLAN.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    if not rows:
        raise SystemExit(f"ERRORE: piano vuoto: {PLAN}")

    total_gb = 0.0
    total_fastq_files = 0
    aria_lines = []
    md5_lines = []

    detailed = []

    for row in rows:
        run = row.get("run_accession", "").strip()
        strategy = row.get("library_strategy", "").strip()
        layout = row.get("library_layout", "").strip()
        sample = row.get("sample_alias", "").strip()
        condition = row.get("possible_condition_from_text", "").strip()
        title = row.get("experiment_title", "").strip()

        try:
            fastq_gb = float(str(row.get("fastq_gb", "0")).replace(",", "."))
        except Exception:
            fastq_gb = 0.0

        urls = [https_from_ftp(x) for x in split_multi(row.get("fastq_ftp", ""))]
        md5s = split_multi(row.get("fastq_md5", ""))

        total_gb += fastq_gb
        total_fastq_files += len(urls)

        for idx, url in enumerate(urls):
            fname = Path(url).name
            aria_lines.append(url)

            md5 = md5s[idx] if idx < len(md5s) else ""
            if md5:
                md5_lines.append(f"{md5}  {fname}")

        detailed.append(
            {
                "run_accession": run,
                "library_strategy": strategy,
                "library_layout": layout,
                "fastq_gb": fastq_gb,
                "n_fastq_urls": len(urls),
                "n_md5": len(md5s),
                "sample_alias": sample,
                "possible_condition_from_text": condition,
                "experiment_title": title,
                "fastq_files": [Path(u).name for u in urls],
            }
        )

    ARIA2_INPUT.write_text("\n".join(aria_lines) + "\n", encoding="utf-8")
    EXPECTED_MD5.write_text("\n".join(md5_lines) + "\n", encoding="utf-8")

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "input_plan": str(PLAN),
        "n_runs": len(rows),
        "estimated_total_fastq_gb": total_gb,
        "total_fastq_files": total_fastq_files,
        "aria2_input": str(ARIA2_INPUT),
        "expected_md5": str(EXPECTED_MD5),
        "runs": detailed,
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("DETTAGLIO PIANO DOWNLOAD RNA-SEQ PRJNA419326")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Input: {PLAN}")
    lines.append(f"Run RNA-Seq: {len(rows)}")
    lines.append(f"File FASTQ totali: {total_fastq_files}")
    lines.append(f"Dimensione totale stimata GB: {total_gb:.3f}")
    lines.append("")
    lines.append("Run dettagliate:")
    lines.append("run | strategy | layout | fastq_gb | n_fastq | n_md5 | sample_alias | condition")
    for r in detailed:
        lines.append(
            f"{r['run_accession']} | {r['library_strategy']} | {r['library_layout']} | "
            f"{r['fastq_gb']:.3f} | {r['n_fastq_urls']} | {r['n_md5']} | "
            f"{r['sample_alias']} | {r['possible_condition_from_text']}"
        )
    lines.append("")
    lines.append("File creati:")
    lines.append(str(ARIA2_INPUT))
    lines.append(str(EXPECTED_MD5))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
