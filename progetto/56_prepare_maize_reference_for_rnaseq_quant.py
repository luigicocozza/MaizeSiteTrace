#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path


BASE_DIR = Path(".").resolve()

REF_DIR = BASE_DIR / "data/reference/zea_mays_ensembl_plants"
OUT_DIR = BASE_DIR / "outputs/prjna419326"

REF_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_TXT = OUT_DIR / "PRJNA419326_rnaseq_reference_prepare_summary.txt"
REPORT_JSON = OUT_DIR / "PRJNA419326_rnaseq_reference_prepare_report.json"

# Ensembl Plants / current
CDNA_URL = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/fasta/zea_mays/cdna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz"
NCRNA_URL = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/fasta/zea_mays/ncrna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.ncrna.fa.gz"
GTF_URL = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/plants/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.63.gtf.gz"

CDNA_GZ = REF_DIR / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz"
NCRNA_GZ = REF_DIR / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.ncrna.fa.gz"
GTF_GZ = REF_DIR / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.63.gtf.gz"

TRANSCRIPTOME_GZ = REF_DIR / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna_plus_ncrna.fa.gz"
TX2GENE_TSV = REF_DIR / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.tx2gene.tsv"

SALMON_INDEX_DIR = REF_DIR / "salmon_index_cdna_plus_ncrna"


def tool_path(name: str) -> str:
    return shutil.which(name) or ""


def run(cmd: list[str]) -> tuple[int, str, str]:
    completed = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return completed.returncode, completed.stdout, completed.stderr


def download(url: str, out: Path) -> bool:
    if out.exists() and out.stat().st_size > 0:
        return True

    cmd = [
        "wget",
        "-c",
        "-O",
        str(out),
        url,
    ]

    code, stdout, stderr = run(cmd)
    return code == 0 and out.exists() and out.stat().st_size > 0


def parse_gene_from_header(header: str) -> tuple[str, str, str]:
    """
    Header Ensembl tipico:
    >Zm00001eb000010_T001 cdna chromosome:...
    gene:Zm00001eb000010 gene_biotype:protein_coding transcript_biotype:protein_coding
    """
    h = header.replace(">", "").strip()
    transcript_id = h.split()[0]

    gene_id = ""
    gene_biotype = ""

    for token in h.split():
        if token.startswith("gene:"):
            gene_id = token.replace("gene:", "")
        elif token.startswith("gene_biotype:"):
            gene_biotype = token.replace("gene_biotype:", "")

    if not gene_id:
        if "_T" in transcript_id:
            gene_id = transcript_id.split("_T")[0]
        else:
            gene_id = transcript_id

    return transcript_id, gene_id, gene_biotype


def build_combined_transcriptome_and_tx2gene() -> dict:
    headers = []
    n_cdna = 0
    n_ncrna = 0

    with gzip.open(TRANSCRIPTOME_GZ, "wt", encoding="utf-8") as out:
        for source_name, path in [("cdna", CDNA_GZ), ("ncrna", NCRNA_GZ)]:
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                for line in f:
                    out.write(line)
                    if line.startswith(">"):
                        headers.append((source_name, line.strip()))
                        if source_name == "cdna":
                            n_cdna += 1
                        else:
                            n_ncrna += 1

    with TX2GENE_TSV.open("w", encoding="utf-8") as f:
        f.write("transcript_id\tgene_id\tgene_biotype\tsource\n")
        for source_name, header in headers:
            tx, gene, biotype = parse_gene_from_header(header)
            f.write(f"{tx}\t{gene}\t{biotype}\t{source_name}\n")

    return {
        "n_transcripts_total": len(headers),
        "n_cdna": n_cdna,
        "n_ncrna": n_ncrna,
    }


def build_salmon_index() -> tuple[bool, str]:
    salmon = tool_path("salmon")
    if not salmon:
        return False, "salmon_non_trovato"

    if SALMON_INDEX_DIR.exists() and any(SALMON_INDEX_DIR.iterdir()):
        return True, "index_gia_presente"

    cmd = [
        salmon,
        "index",
        "-t",
        str(TRANSCRIPTOME_GZ),
        "-i",
        str(SALMON_INDEX_DIR),
        "-k",
        "31",
    ]

    code, stdout, stderr = run(cmd)

    ok = code == 0 and SALMON_INDEX_DIR.exists() and any(SALMON_INDEX_DIR.iterdir())

    log_path = OUT_DIR / "PRJNA419326_salmon_index_build.log"
    log_path.write_text(
        "COMMAND:\n"
        + " ".join(cmd)
        + "\n\nSTDOUT:\n"
        + stdout
        + "\n\nSTDERR:\n"
        + stderr
        + "\n",
        encoding="utf-8",
    )

    return ok, str(log_path)


def main() -> int:
    tools = {
        "wget": tool_path("wget"),
        "gzip": tool_path("gzip"),
        "salmon": tool_path("salmon"),
    }

    missing_required = [k for k in ["wget", "gzip"] if not tools[k]]

    if missing_required:
        raise SystemExit(f"ERRORE: strumenti mancanti: {missing_required}")

    download_results = {
        "cdna": download(CDNA_URL, CDNA_GZ),
        "ncrna": download(NCRNA_URL, NCRNA_GZ),
        "gtf": download(GTF_URL, GTF_GZ),
    }

    if not all(download_results.values()):
        raise SystemExit(f"ERRORE download reference: {download_results}")

    combined_stats = build_combined_transcriptome_and_tx2gene()

    salmon_index_ok, salmon_index_note = build_salmon_index()

    status = "PASS" if salmon_index_ok else "PARTIAL_PASS_NO_SALMON_INDEX"

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "status": status,
        "tools": tools,
        "reference_dir": str(REF_DIR),
        "download_results": download_results,
        "urls": {
            "cdna": CDNA_URL,
            "ncrna": NCRNA_URL,
            "gtf": GTF_URL,
        },
        "files": {
            "cdna_gz": str(CDNA_GZ),
            "ncrna_gz": str(NCRNA_GZ),
            "gtf_gz": str(GTF_GZ),
            "transcriptome_gz": str(TRANSCRIPTOME_GZ),
            "tx2gene_tsv": str(TX2GENE_TSV),
            "salmon_index_dir": str(SALMON_INDEX_DIR),
        },
        "combined_stats": combined_stats,
        "salmon_index_ok": salmon_index_ok,
        "salmon_index_note": salmon_index_note,
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("PREPARAZIONE REFERENCE MAIS PER QUANTIFICAZIONE RNA-SEQ")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Status: {status}")
    lines.append("")
    lines.append("Strumenti:")
    for k, v in tools.items():
        lines.append(f"- {k}: {v if v else 'NON TROVATO'}")
    lines.append("")
    lines.append("Download reference:")
    for k, v in download_results.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("Transcriptome combinato:")
    lines.append(f"- cDNA transcripts: {combined_stats['n_cdna']}")
    lines.append(f"- ncRNA transcripts: {combined_stats['n_ncrna']}")
    lines.append(f"- totale transcripts: {combined_stats['n_transcripts_total']}")
    lines.append("")
    lines.append(f"Salmon index OK: {salmon_index_ok}")
    lines.append(f"Salmon index note: {salmon_index_note}")
    lines.append("")
    lines.append("Output:")
    lines.append(str(CDNA_GZ))
    lines.append(str(NCRNA_GZ))
    lines.append(str(GTF_GZ))
    lines.append(str(TRANSCRIPTOME_GZ))
    lines.append(str(TX2GENE_TSV))
    lines.append(str(SALMON_INDEX_DIR))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    if status == "PASS":
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
