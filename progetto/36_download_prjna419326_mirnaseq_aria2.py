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

REPORT_JSON = OUT_DIR / f"{PROJECT}_mirnaseq_aria2_download_report.json"
REPORT_TXT = OUT_DIR / f"{PROJECT}_mirnaseq_aria2_download_report.txt"


def read_plan() -> list[dict]:
    with PLAN_TSV.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def split_field(value: str) -> list[str]:
    if value is None:
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


def to_ftp_url(value: str) -> str:
    value = value.strip()
    value = value.replace("https://", "")
    value = value.replace("http://", "")
    value = value.replace("ftp://", "")
    return "ftp://" + value


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_test(path: Path) -> bool:
    proc = subprocess.run(["gzip", "-t", str(path)])
    return proc.returncode == 0


def prepare_items(rows: list[dict]) -> list[dict]:
    items = []

    for row in rows:
        urls = split_field(row.get("fastq_ftp", ""))
        md5s = split_field(row.get("fastq_md5", ""))

        for idx, raw_url in enumerate(urls):
            filename = Path(raw_url).name
            expected_md5 = md5s[idx] if idx < len(md5s) else ""

            items.append(
                {
                    "run_accession": row.get("run_accession", ""),
                    "sample_alias": row.get("sample_alias", ""),
                    "genotype": row.get("genotype", ""),
                    "condition": row.get("condition", ""),
                    "replicate": row.get("replicate", ""),
                    "library_layout": row.get("library_layout", ""),
                    "fastq_gb": row.get("fastq_gb", ""),
                    "url": to_ftp_url(raw_url),
                    "filename": filename,
                    "expected_md5": expected_md5,
                    "out_path": str(FASTQ_DIR / filename),
                }
            )

    return items


def aria2_download(url: str, out_path: Path, expected_md5: str) -> int:
    aria2c = shutil.which("aria2c")

    if not aria2c:
        raise RuntimeError("aria2c non trovato. Installa con: sudo apt install -y aria2")

    cmd = [
        aria2c,
        f"--dir={out_path.parent}",
        f"--out={out_path.name}",
        "--max-connection-per-server=1",
        "--split=1",
        "--continue=true",
        "--auto-file-renaming=false",
        "--allow-overwrite=true",
        "--retry-wait=30",
        "--max-tries=20",
        "--timeout=120",
        "--connect-timeout=120",
    ]

    if expected_md5:
        cmd.append(f"--checksum=md5={expected_md5}")

    cmd.append(url)

    print("")
    print("Comando aria2c:")
    print(" ".join(cmd))
    print("")

    proc = subprocess.run(cmd)
    return proc.returncode


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="Esegue il download reale.")
    parser.add_argument("--max-files", type=int, default=0, help="Limita il numero di file processati. 0 = tutti.")
    args = parser.parse_args()

    FASTQ_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = read_plan()
    items = prepare_items(rows)

    if args.max_files and args.max_files > 0:
        items = items[: args.max_files]

    results = []

    print(f"ARIA2 DOWNLOAD PLAN - {PROJECT} miRNA/small RNA")
    print("=" * 80)
    print(f"Run nel piano: {len(rows)}")
    print(f"File FASTQ da processare: {len(items)}")
    print(f"Cartella FASTQ: {FASTQ_DIR}")
    print(f"Execute: {args.execute}")
    print("")

    for item in items:
        out_path = Path(item["out_path"])
        expected_md5 = item["expected_md5"]

        print("=" * 80)
        print(f"Run: {item['run_accession']}")
        print(f"Sample: {item['sample_alias']}")
        print(f"Genotype: {item['genotype']}")
        print(f"Condition: {item['condition']}")
        print(f"Replicate: {item['replicate']}")
        print(f"File: {item['filename']}")
        print(f"Expected MD5: {expected_md5}")
        print(f"URL FTP: {item['url']}")

        status = "planned"
        actual_md5 = ""
        size_bytes = 0
        gzip_ok = False

        if out_path.exists():
            size_bytes = out_path.stat().st_size
            actual_md5 = md5_file(out_path)
            gzip_ok = gzip_test(out_path)

            print(f"File già presente: {out_path}")
            print(f"MD5 locale: {actual_md5}")
            print(f"GZIP OK: {gzip_ok}")

            if expected_md5 and actual_md5 == expected_md5 and gzip_ok:
                status = "already_present_md5_gzip_ok"
                print("File già valido. Download saltato.")
            else:
                status = "existing_invalid"
                print("File esistente non valido.")

                if args.execute:
                    print("Cancello file non valido prima del nuovo download.")
                    out_path.unlink(missing_ok=True)
                    aria2_tmp = Path(str(out_path) + ".aria2")
                    aria2_tmp.unlink(missing_ok=True)

        if status in {"planned", "existing_invalid"}:
            if not args.execute:
                status = "dry_run_planned"
                print("Dry-run: download non eseguito.")
            else:
                rc = aria2_download(item["url"], out_path, expected_md5)

                if rc != 0:
                    status = f"aria2_failed_rc_{rc}"
                    print(f"ERRORE aria2c, return code {rc}")
                else:
                    size_bytes = out_path.stat().st_size if out_path.exists() else 0
                    actual_md5 = md5_file(out_path) if out_path.exists() else ""
                    gzip_ok = gzip_test(out_path) if out_path.exists() else False

                    print(f"MD5 locale dopo download: {actual_md5}")
                    print(f"GZIP OK dopo download: {gzip_ok}")

                    if expected_md5 and actual_md5 == expected_md5 and gzip_ok:
                        status = "downloaded_md5_gzip_ok"
                        print("Download valido.")
                    else:
                        status = "downloaded_invalid"
                        print("Download non valido.")

        results.append(
            {
                **item,
                "status": status,
                "actual_md5": actual_md5,
                "gzip_ok": gzip_ok,
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
    lines.append(f"ARIA2 DOWNLOAD REPORT - {PROJECT} miRNA/small RNA")
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
            f"{result['filename']} | {result['status']} | gzip_ok={result['gzip_ok']}"
        )

    REPORT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("=" * 80)
    print("REPORT FINALE")
    print("=" * 80)
    print("\n".join(lines))
    print("")
    print(f"Report JSON: {REPORT_JSON}")
    print(f"Report TXT: {REPORT_TXT}")

    bad = [
        r
        for r in results
        if not (
            r["status"] in {"already_present_md5_gzip_ok", "downloaded_md5_gzip_ok"}
            and r["gzip_ok"] is True
        )
    ]

    if bad:
        print("")
        print("ATTENZIONE: alcuni file non sono validi.")
        return 1

    print("")
    print("Tutti i file processati sono validi.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
