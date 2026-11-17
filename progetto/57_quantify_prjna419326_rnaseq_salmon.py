#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import subprocess
from datetime import datetime
from pathlib import Path

import pandas as pd


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

SAMPLE_SHEET = BASE_DIR / "outputs/prjna419326/PRJNA419326_rnaseq_sample_sheet.tsv"

REF_DIR = BASE_DIR / "data/reference/zea_mays_ensembl_plants"
SALMON_INDEX = REF_DIR / "salmon_index_cdna_plus_ncrna"
TX2GENE = REF_DIR / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.tx2gene.tsv"

QUANT_DIR = BASE_DIR / "data/processed/rnaseq_salmon_quant"
MATRIX_DIR = BASE_DIR / "data/processed/rnaseq_matrices"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

QUANT_DIR.mkdir(parents=True, exist_ok=True)
MATRIX_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_TXT = OUT_DIR / "PRJNA419326_rnaseq_salmon_quant_summary.txt"
REPORT_JSON = OUT_DIR / "PRJNA419326_rnaseq_salmon_quant_report.json"
RUN_LOG = OUT_DIR / "PRJNA419326_rnaseq_salmon_quant_run.log"

TRANSCRIPT_TPM = MATRIX_DIR / "PRJNA419326_rnaseq_transcript_TPM_matrix.tsv"
TRANSCRIPT_COUNTS = MATRIX_DIR / "PRJNA419326_rnaseq_transcript_NumReads_matrix.tsv"
GENE_TPM = MATRIX_DIR / "PRJNA419326_rnaseq_gene_TPM_matrix.tsv"
GENE_COUNTS = MATRIX_DIR / "PRJNA419326_rnaseq_gene_NumReads_matrix.tsv"


THREADS = 2


def append_log(msg: str) -> None:
    print(msg, flush=True)
    with RUN_LOG.open("a", encoding="utf-8") as f:
        f.write(msg + "\n")


def run_cmd(cmd: list[str]) -> tuple[int, str, str]:
    completed = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return completed.returncode, completed.stdout, completed.stderr


def read_sample_sheet() -> list[dict]:
    with SAMPLE_SHEET.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def run_salmon_for_sample(sample: dict) -> dict:
    sample_alias = sample["sample_alias"]
    fastq_1 = Path(sample["fastq_1"])
    fastq_2 = Path(sample["fastq_2"])

    out_dir = QUANT_DIR / sample_alias
    quant_sf = out_dir / "quant.sf"
    log_file = out_dir / "salmon_quant.log"

    if quant_sf.exists() and quant_sf.stat().st_size > 0:
        return {
            "sample_alias": sample_alias,
            "run_accession": sample["run_accession"],
            "status": "already_present",
            "quant_sf": str(quant_sf),
            "log_file": str(log_file),
            "returncode": 0,
        }

    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "salmon",
        "quant",
        "-i",
        str(SALMON_INDEX),
        "-l",
        "A",
        "-1",
        str(fastq_1),
        "-2",
        str(fastq_2),
        "-p",
        str(THREADS),
        "--validateMappings",
        "-o",
        str(out_dir),
    ]

    append_log("=" * 80)
    append_log(f"Campione: {sample_alias}")
    append_log("Comando:")
    append_log(" ".join(cmd))

    code, stdout, stderr = run_cmd(cmd)

    log_file.write_text(
        "COMMAND:\n"
        + " ".join(cmd)
        + "\n\nSTDOUT:\n"
        + stdout
        + "\n\nSTDERR:\n"
        + stderr
        + "\n",
        encoding="utf-8",
    )

    status = "OK" if code == 0 and quant_sf.exists() and quant_sf.stat().st_size > 0 else "FAIL"

    append_log(f"Status: {status}")
    append_log(f"Return code: {code}")
    append_log(f"Output: {quant_sf}")

    return {
        "sample_alias": sample_alias,
        "run_accession": sample["run_accession"],
        "status": status,
        "quant_sf": str(quant_sf),
        "log_file": str(log_file),
        "returncode": code,
    }


def build_matrices(samples: list[dict], quant_results: list[dict]) -> dict:
    ok_results = [r for r in quant_results if r["status"] in {"OK", "already_present"}]

    transcript_tpm = {}
    transcript_counts = {}

    for r in ok_results:
        sample_alias = r["sample_alias"]
        q = pd.read_csv(r["quant_sf"], sep="\t")

        transcript_tpm[sample_alias] = q.set_index("Name")["TPM"]
        transcript_counts[sample_alias] = q.set_index("Name")["NumReads"]

    tpm_df = pd.DataFrame(transcript_tpm).fillna(0.0)
    counts_df = pd.DataFrame(transcript_counts).fillna(0.0)

    tpm_df.index.name = "transcript_id"
    counts_df.index.name = "transcript_id"

    tpm_df.to_csv(TRANSCRIPT_TPM, sep="\t")
    counts_df.to_csv(TRANSCRIPT_COUNTS, sep="\t")

    tx2gene = pd.read_csv(TX2GENE, sep="\t")
    tx2gene = tx2gene.drop_duplicates(subset=["transcript_id"])

    tpm_long = tpm_df.reset_index().merge(tx2gene, on="transcript_id", how="left")
    counts_long = counts_df.reset_index().merge(tx2gene, on="transcript_id", how="left")

    tpm_long["gene_id"] = tpm_long["gene_id"].fillna(tpm_long["transcript_id"])
    counts_long["gene_id"] = counts_long["gene_id"].fillna(counts_long["transcript_id"])

    sample_cols = [r["sample_alias"] for r in ok_results]

    gene_tpm = tpm_long.groupby("gene_id")[sample_cols].sum()
    gene_counts = counts_long.groupby("gene_id")[sample_cols].sum()

    gene_tpm.index.name = "gene_id"
    gene_counts.index.name = "gene_id"

    gene_tpm.to_csv(GENE_TPM, sep="\t")
    gene_counts.to_csv(GENE_COUNTS, sep="\t")

    return {
        "n_transcripts": int(tpm_df.shape[0]),
        "n_genes": int(gene_tpm.shape[0]),
        "n_samples": len(ok_results),
        "transcript_tpm": str(TRANSCRIPT_TPM),
        "transcript_counts": str(TRANSCRIPT_COUNTS),
        "gene_tpm": str(GENE_TPM),
        "gene_counts": str(GENE_COUNTS),
    }


def main() -> int:
    if RUN_LOG.exists():
        RUN_LOG.unlink()

    append_log("QUANTIFICAZIONE RNA-SEQ PRJNA419326 CON SALMON")
    append_log("=" * 80)
    append_log(f"Data inizio: {datetime.now().isoformat(timespec='seconds')}")
    append_log(f"Sample sheet: {SAMPLE_SHEET}")
    append_log(f"Salmon index: {SALMON_INDEX}")
    append_log("")

    if not SAMPLE_SHEET.exists():
        raise SystemExit(f"ERRORE: sample sheet non trovato: {SAMPLE_SHEET}")

    if not SALMON_INDEX.exists():
        raise SystemExit(f"ERRORE: indice Salmon non trovato: {SALMON_INDEX}")

    samples = read_sample_sheet()

    samples = [s for s in samples if s.get("sample_status") == "OK"]

    quant_results = []

    for i, sample in enumerate(samples, start=1):
        append_log("")
        append_log(f"[{i}/{len(samples)}] {sample['sample_alias']}")
        result = run_salmon_for_sample(sample)
        quant_results.append(result)

        if result["status"] == "FAIL":
            append_log("ERRORE: quantificazione fallita. Interrompo.")
            break

    n_ok = sum(1 for r in quant_results if r["status"] in {"OK", "already_present"})
    n_fail = sum(1 for r in quant_results if r["status"] == "FAIL")

    final_status = "PASS" if n_ok == len(samples) and n_fail == 0 else "FAIL"

    matrix_stats = {}

    if final_status == "PASS":
        matrix_stats = build_matrices(samples, quant_results)

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "status": final_status,
        "sample_sheet": str(SAMPLE_SHEET),
        "salmon_index": str(SALMON_INDEX),
        "n_samples_expected": len(samples),
        "n_quant_ok": n_ok,
        "n_quant_fail": n_fail,
        "threads": THREADS,
        "quant_dir": str(QUANT_DIR),
        "matrix_stats": matrix_stats,
        "results": quant_results,
        "outputs": {
            "summary_txt": str(SUMMARY_TXT),
            "report_json": str(REPORT_JSON),
            "run_log": str(RUN_LOG),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("QUANTIFICAZIONE RNA-SEQ PRJNA419326 CON SALMON")
    lines.append("=" * 80)
    lines.append(f"Data fine: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Status: {final_status}")
    lines.append(f"Campioni attesi: {len(samples)}")
    lines.append(f"Campioni quantificati OK: {n_ok}")
    lines.append(f"Campioni falliti: {n_fail}")
    lines.append("")
    lines.append("Output quantificazioni:")
    lines.append(str(QUANT_DIR))

    if matrix_stats:
        lines.append("")
        lines.append("Matrici create:")
        lines.append(f"Transcripti: {matrix_stats['n_transcripts']}")
        lines.append(f"Geni: {matrix_stats['n_genes']}")
        lines.append(str(TRANSCRIPT_TPM))
        lines.append(str(TRANSCRIPT_COUNTS))
        lines.append(str(GENE_TPM))
        lines.append(str(GENE_COUNTS))

    lines.append("")
    lines.append("Report:")
    lines.append(str(SUMMARY_TXT))
    lines.append(str(REPORT_JSON))
    lines.append(str(RUN_LOG))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("\n".join(lines))

    return 0 if final_status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
