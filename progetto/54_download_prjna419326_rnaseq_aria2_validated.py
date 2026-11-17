#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

PLAN = BASE_DIR / "outputs/prjna419326/PRJNA419326_download_plan_rnaseq_only.tsv"

FASTQ_DIR = BASE_DIR / "data/public_ena/PRJNA419326/rnaseq_fastq"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

LOG_TXT = OUT_DIR / "PRJNA419326_rnaseq_download_aria2_validated_log.txt"
SUMMARY_TXT = OUT_DIR / "PRJNA419326_rnaseq_download_VALIDATO.txt"
REPORT_JSON = OUT_DIR / "PRJNA419326_rnaseq_download_VALIDATO_report.json"
EXPECTED_MD5 = OUT_DIR / "PRJNA419326_rnaseq_expected_md5.txt"

ARIA2C = shutil.which("aria2c")


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


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024 * 8), b""):
            h.update(block)
    return h.hexdigest()


def append_log(msg: str) -> None:
    print(msg, flush=True)
    with LOG_TXT.open("a", encoding="utf-8") as f:
        f.write(msg + "\n")


def read_plan() -> list[dict]:
    with PLAN.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build_download_items(rows: list[dict]) -> list[dict]:
    items = []

    for row in rows:
        run = row["run_accession"]
        sample_alias = row.get("sample_alias", "")
        layout = row.get("library_layout", "")
        strategy = row.get("library_strategy", "")
        fastq_gb = row.get("fastq_gb", "")

        urls = [https_from_ftp(x) for x in split_multi(row.get("fastq_ftp", ""))]
        md5s = split_multi(row.get("fastq_md5", ""))

        for i, url in enumerate(urls):
            fname = Path(url).name
            expected = md5s[i] if i < len(md5s) else ""

            items.append(
                {
                    "run_accession": run,
                    "sample_alias": sample_alias,
                    "library_strategy": strategy,
                    "library_layout": layout,
                    "fastq_gb_run": fastq_gb,
                    "url": url,
                    "filename": fname,
                    "expected_md5": expected,
                    "local_path": str(FASTQ_DIR / fname),
                    "read_pair_index": i + 1,
                }
            )

    return items


def run_aria2(url: str, out_dir: Path, filename: str) -> int:
    cmd = [
        ARIA2C,
        "--continue=true",
        "--max-connection-per-server=1",
        "--split=1",
        "--min-split-size=100M",
        "--max-tries=0",
        "--retry-wait=20",
        "--timeout=120",
        "--connect-timeout=60",
        "--file-allocation=none",
        "--summary-interval=60",
        "--console-log-level=notice",
        "--dir", str(out_dir),
        "--out", filename,
        url,
    ]

    append_log("Comando aria2c:")
    append_log(" ".join(cmd))

    completed = subprocess.run(cmd)
    return int(completed.returncode)


def main() -> int:
    if ARIA2C is None:
        raise SystemExit("ERRORE: aria2c non trovato. Installa con: sudo apt install aria2")

    FASTQ_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if LOG_TXT.exists():
        LOG_TXT.unlink()

    rows = read_plan()
    items = build_download_items(rows)

    EXPECTED_MD5.write_text(
        "\n".join(
            f"{item['expected_md5']}  {item['filename']}"
            for item in items
            if item["expected_md5"]
        )
        + "\n",
        encoding="utf-8",
    )

    append_log("DOWNLOAD RNA-SEQ PRJNA419326 CON VALIDAZIONE MD5")
    append_log("=" * 80)
    append_log(f"Data inizio: {datetime.now().isoformat(timespec='seconds')}")
    append_log(f"Run RNA-Seq: {len(rows)}")
    append_log(f"FASTQ totali: {len(items)}")
    append_log(f"Cartella output: {FASTQ_DIR}")
    append_log("")

    results = []

    for idx, item in enumerate(items, start=1):
        path = Path(item["local_path"])
        expected = item["expected_md5"]

        append_log("=" * 80)
        append_log(f"[{idx}/{len(items)}] {item['filename']}")
        append_log(f"Run: {item['run_accession']} | sample: {item['sample_alias']} | pair: {item['read_pair_index']}")
        append_log(f"URL: {item['url']}")

        status = "unknown"
        observed_md5 = ""
        aria2_exit = None

        if path.exists() and expected:
            observed_md5 = md5sum(path)
            if observed_md5 == expected:
                status = "already_present_md5_ok"
                append_log(f"Già presente e MD5 OK: {path}")
            else:
                append_log(f"File presente ma MD5 NON OK: {path}")
                append_log(f"expected={expected}")
                append_log(f"observed={observed_md5}")
                append_log("Rilancio aria2c con resume.")
                aria2_exit = run_aria2(item["url"], FASTQ_DIR, item["filename"])
                observed_md5 = md5sum(path) if path.exists() else ""
                status = "downloaded_md5_ok" if observed_md5 == expected else "downloaded_md5_fail"

        else:
            aria2_exit = run_aria2(item["url"], FASTQ_DIR, item["filename"])

            if path.exists() and expected:
                observed_md5 = md5sum(path)
                status = "downloaded_md5_ok" if observed_md5 == expected else "downloaded_md5_fail"
            elif path.exists():
                status = "downloaded_no_md5_available"
            else:
                status = "download_failed_no_file"

        append_log(f"Status: {status}")
        if expected:
            append_log(f"Expected MD5: {expected}")
            append_log(f"Observed MD5: {observed_md5}")

        results.append(
            {
                **item,
                "status": status,
                "aria2_exit": aria2_exit,
                "observed_md5": observed_md5,
                "size_bytes": path.stat().st_size if path.exists() else 0,
            }
        )

        if status.endswith("fail") or status == "download_failed_no_file":
            append_log("")
            append_log("ERRORE: download o MD5 fallito. Interrompo per evitare di proseguire con dati corrotti.")
            break

        append_log("")

    ok = [
        r for r in results
        if r["status"] in {"already_present_md5_ok", "downloaded_md5_ok", "downloaded_no_md5_available"}
    ]

    failed = [
        r for r in results
        if r["status"].endswith("fail") or r["status"] == "download_failed_no_file"
    ]

    total_size = sum(r["size_bytes"] for r in results)
    total_size_gb = total_size / 1e9

    final_status = "PASS" if len(ok) == len(items) and not failed else "FAIL"

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "final_status": final_status,
        "plan": str(PLAN),
        "fastq_dir": str(FASTQ_DIR),
        "n_runs": len(rows),
        "n_expected_fastq": len(items),
        "n_processed": len(results),
        "n_ok": len(ok),
        "n_failed": len(failed),
        "downloaded_or_validated_size_gb": total_size_gb,
        "expected_md5_file": str(EXPECTED_MD5),
        "summary_txt": str(SUMMARY_TXT),
        "log_txt": str(LOG_TXT),
        "results": results,
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("RNA-SEQ PRJNA419326 DOWNLOAD VALIDATO")
    lines.append("=" * 80)
    lines.append(f"Data fine: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Status: {final_status}")
    lines.append(f"Run attese: {len(rows)}")
    lines.append(f"FASTQ attesi: {len(items)}")
    lines.append(f"FASTQ processati: {len(results)}")
    lines.append(f"FASTQ OK: {len(ok)}")
    lines.append(f"FASTQ falliti: {len(failed)}")
    lines.append(f"Dimensione validata GB: {total_size_gb:.3f}")
    lines.append("")
    lines.append("Cartella FASTQ:")
    lines.append(str(FASTQ_DIR))
    lines.append("")
    lines.append("Output:")
    lines.append(str(SUMMARY_TXT))
    lines.append(str(REPORT_JSON))
    lines.append(str(LOG_TXT))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    append_log("\n".join(lines))

    return 0 if final_status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
