#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


BASE_DIR = Path(".").resolve()
PROJECT = "PRJNA419326"

PLAN_TSV = BASE_DIR / "outputs/prjna419326" / f"{PROJECT}_download_plan_mirnaseq_only.tsv"
FASTQ_DIR = BASE_DIR / "data/public_ena" / PROJECT / "fastq" / "mirnaseq"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

REPORT_JSON = OUT_DIR / f"{PROJECT}_mirnaseq_download_report.json"
REPORT_TXT = OUT_DIR / f"{PROJECT}_mirnaseq_download_report.txt"


def read_plan() -> list[dict]:
    with PLAN_TSV.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def normalize_url(value: str) -> str:
    value = value.strip()

    if not value:
        return ""

    if value.startswith("http://") or value.startswith("https://") or value.startswith("ftp://"):
        return value

    return "https://" + value


def split_field(value: str) -> list[str]:
    if value is None:
        return []

    return [x.strip() for x in str(value).split(";") if x.strip()]


def md5_file(path: Path) -> str:
    h = hashlib.md5()

    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)

    return h.hexdigest()


def downloader_command(url: str, out_path: Path) -> list[str]:
    wget = shutil.which("wget")

    if wget:
        return [
            wget,
            "-c",
            "--tries=10",
            "--timeout=60",
            "--read-timeout=60",
            "-O",
            str(out_path),
            url,
        ]

    curl = shutil.which("curl")

    if curl:
        return [
            curl,
            "-L",
            "-C",
            "-",
            "--retry",
            "10",
            "--connect-timeout",
            "60",
            "-o",
            str(out_path),
            url,
        ]

    raise RuntimeError("Ne' wget ne' curl risultano disponibili nel sistema.")


def download_file(url: str, out_path: Path) -> int:
    cmd = downloader_command(url, out_path)

    print("")
    print("Comando download:")
    print(" ".join(cmd))
    print("")

    proc = subprocess.run(cmd)

    return proc.returncode


def prepare_items(rows: list[dict]) -> list[dict]:
    items = []

    for row in rows:
        urls = split_field(row.get("fastq_ftp", ""))
        md5s = split_field(row.get("fastq_md5", ""))

        for idx, raw_url in enumerate(urls):
            url = normalize_url(raw_url)
            filename = Path(raw_url).name
            expected_md5 = md5s[idx] if idx < len(md5s) else ""

            items.append(
                {
                    "project": PROJECT,
                    "run_accession": row.get("run_accession", ""),
                    "sample_alias": row.get("sample_alias", ""),
                    "genotype": row.get("genotype", ""),
                    "condition": row.get("condition", ""),
                    "replicate": row.get("replicate", ""),
                    "sequencing_type": row.get("sequencing_type", ""),
                    "library_layout": row.get("library_layout", ""),
                    "fastq_gb": row.get("fastq_gb", ""),
                    "url": url,
                    "filename": filename,
                    "expected_md5": expected_md5,
                    "out_path": str(FASTQ_DIR / filename),
                }
            )

    return items


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="Esegue realmente il download.")
    parser.add_argument("--max-files", type=int, default=0, help="Scarica al massimo N file. 0 = tutti.")
    args = parser.parse_args()

    FASTQ_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = read_plan()
    items = prepare_items(rows)

    if args.max_files and args.max_files > 0:
        items = items[: args.max_files]

    results = []

    print(f"DOWNLOAD PLAN - {PROJECT} miRNA/small RNA")
    print("=" * 80)
    print(f"Run nel piano: {len(rows)}")
    print(f"File FASTQ nel piano: {len(items)}")
    print(f"Cartella download: {FASTQ_DIR}")
    print(f"Modalita' execute: {args.execute}")
    print("")

    for item in items:
        out_path = Path(item["out_path"])
        expected_md5 = item["expected_md5"]

        print("-" * 80)
        print(f"Run: {item['run_accession']}")
        print(f"Sample: {item['sample_alias']}")
        print(f"Genotype: {item['genotype']}")
        print(f"Condition: {item['condition']}")
        print(f"Replicate: {item['replicate']}")
        print(f"File: {item['filename']}")
        print(f"URL: {item['url']}")
        print(f"Expected MD5: {expected_md5}")

        status = "planned"
        actual_md5 = ""
        size_bytes = 0

        if out_path.exists():
            size_bytes = out_path.stat().st_size
            print(f"File gia' presente: {out_path} ({size_bytes} bytes)")

            if expected_md5:
                actual_md5 = md5_file(out_path)

                if actual_md5 == expected_md5:
                    status = "already_present_md5_ok"
                    print("MD5 OK. Download saltato.")
                else:
                    status = "already_present_md5_mismatch"
                    print("MD5 NON corrisponde. Il file verra' riscaricato se --execute e' attivo.")
            else:
                status = "already_present_no_md5"
                print("MD5 atteso non disponibile. Download saltato per sicurezza.")

        if status in {"planned", "already_present_md5_mismatch"}:
            if not args.execute:
                print("Dry-run: download non eseguito.")
                status = "dry_run_planned"
            else:
                rc = download_file(item["url"], out_path)

                if rc != 0:
                    status = f"download_failed_rc_{rc}"
                    print(f"ERRORE download, return code: {rc}")
                else:
                    size_bytes = out_path.stat().st_size if out_path.exists() else 0

                    if expected_md5:
                        actual_md5 = md5_file(out_path)

                        if actual_md5 == expected_md5:
                            status = "downloaded_md5_ok"
                            print("Download completato. MD5 OK.")
                        else:
                            status = "downloaded_md5_mismatch"
                            print("Download completato ma MD5 NON corrisponde.")
                    else:
                        status = "downloaded_no_md5"
                        print("Download completato. MD5 atteso non disponibile.")

        results.append(
            {
                **item,
                "status": status,
                "actual_md5": actual_md5,
                "size_bytes": size_bytes,
            }
        )

    summary = {
        "project": PROJECT,
        "created": datetime.now().isoformat(timespec="seconds"),
        "execute": args.execute,
        "fastq_dir": str(FASTQ_DIR),
        "n_items": len(results),
        "status_counts": {},
        "results": results,
    }

    for result in results:
        status = result["status"]
        summary["status_counts"][status] = summary["status_counts"].get(status, 0) + 1

    REPORT_JSON.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append(f"DOWNLOAD REPORT - {PROJECT} miRNA/small RNA")
    lines.append("=" * 80)
    lines.append(f"Execute: {args.execute}")
    lines.append(f"FASTQ dir: {FASTQ_DIR}")
    lines.append(f"File totali: {len(results)}")
    lines.append("")
    lines.append("Status counts:")
    for key, value in sorted(summary["status_counts"].items()):
        lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append("File:")
    for result in results:
        lines.append(
            f"{result['run_accession']} | {result['sample_alias']} | "
            f"{result['filename']} | {result['status']}"
        )

    REPORT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("=" * 80)
    print("REPORT")
    print("=" * 80)
    print("\n".join(lines))
    print("")
    print(f"Report JSON: {REPORT_JSON}")
    print(f"Report TXT: {REPORT_TXT}")

    failed = [
        r
        for r in results
        if "failed" in r["status"] or "mismatch" in r["status"]
    ]

    if failed:
        print("")
        print("ATTENZIONE: alcuni file hanno errore o MD5 non corrispondente.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
