#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

OUT_DIR = BASE_DIR / "outputs/prjna419326"
PACKAGE_DIR = OUT_DIR / "v3_final_package"
STAGING_DIR = PACKAGE_DIR / "maize_ncrna_gnn_v3_prjna419326_package"
ZIP_PATH = PACKAGE_DIR / "maize_ncrna_gnn_v3_prjna419326_package.zip"
MANIFEST_JSON = PACKAGE_DIR / "maize_ncrna_gnn_v3_prjna419326_manifest.json"
SUMMARY_TXT = PACKAGE_DIR / "maize_ncrna_gnn_v3_prjna419326_package_summary.txt"

SCRIPTS = [
    "32_fetch_prjna419326_metadata.py",
    "33_analyze_prjna419326_metadata.py",
    "34_create_prjna419326_experimental_design.py",
    "35_download_prjna419326_mirnaseq.py",
    "36_download_prjna419326_mirnaseq_aria2.py",
    "38_light_qc_prjna419326_mirnaseq_fastq.py",
    "39_detect_mirnaseq_adapter_candidates.py",
    "40_trim_all_prjna419326_mirnaseq.py",
    "41_qc_trimmed_prjna419326_mirnaseq.py",
    "42_collapse_trimmed_reads_build_count_matrix.py",
    "43_screen_srna_drought_response.py",
    "44_map_srna_top500_to_rnacentral_sequences.py",
    "45_build_prjna419326_experimental_srna_ranking.py",
    "46_generate_prjna419326_v3_methodological_report.py",
    "47_build_v3_prjna419326_final_package.py",
]

OUTPUT_FILES = [
    f"{PROJECT}_mirnaseq_download_VALIDATO.txt",
    f"{PROJECT}_mirnaseq_md5_final_validation.txt",
    f"{PROJECT}_mirnaseq_light_fastq_qc_summary.txt",
    f"{PROJECT}_mirnaseq_adapter_candidate_report.txt",
    f"{PROJECT}_mirnaseq_cutadapt_trimming_summary.txt",
    f"{PROJECT}_mirnaseq_cutadapt_trimming_summary.csv",
    f"{PROJECT}_mirnaseq_trimmed_qc_summary.txt",
    f"{PROJECT}_mirnaseq_trimmed_qc_summary.csv",
    f"{PROJECT}_mirnaseq_trimmed_length_distribution.csv",
    f"{PROJECT}_mirnaseq_srna_count_matrix_summary.txt",
    f"{PROJECT}_mirnaseq_collapsed_sample_summary.csv",
    f"{PROJECT}_mirnaseq_srna_drought_response_summary.txt",
    f"{PROJECT}_mirnaseq_srna_drought_response_report.json",
    f"{PROJECT}_mirnaseq_srna_drought_response_candidates.tsv",
    f"{PROJECT}_mirnaseq_srna_drought_response_top500.tsv",
    f"{PROJECT}_mirnaseq_srna_drought_response_top500.fasta",
    f"{PROJECT}_top500_srna_to_rnacentral_mapping_summary.txt",
    f"{PROJECT}_top500_srna_to_rnacentral_mapping_report.json",
    f"{PROJECT}_top500_srna_to_rnacentral_exact_matches.tsv",
    f"{PROJECT}_top500_srna_unmatched.tsv",
    f"{PROJECT}_experimental_srna_v3_ranking_summary.txt",
    f"{PROJECT}_experimental_srna_v3_ranking.tsv",
    f"{PROJECT}_experimental_srna_v3_ranking.csv",
    f"{PROJECT}_experimental_srna_v3_ranking_top500.fasta",
    f"{PROJECT}_experimental_srna_v3_ranking_report.json",
    f"{PROJECT}_V3_report_metodologico_integrato.md",
    f"{PROJECT}_V3_report_metodologico_integrato.txt",
    f"{PROJECT}_V3_report_metodologico_integrato.json",
]

METADATA_FILES = [
    "data/public_ena/PRJNA419326/metadata/PRJNA419326_ena_run_metadata.tsv",
    "data/public_ena/PRJNA419326/metadata/PRJNA419326_ena_run_metadata_summary.json",
]

DESIGN_OUTPUTS = [
    "outputs/prjna419326/PRJNA419326_runs_all.csv",
    "outputs/prjna419326/PRJNA419326_runs_rnaseq.csv",
    "outputs/prjna419326/PRJNA419326_runs_mirnaseq.csv",
    "outputs/prjna419326/PRJNA419326_download_plan.tsv",
    "outputs/prjna419326/PRJNA419326_experimental_design_curated.csv",
    "outputs/prjna419326/PRJNA419326_experimental_design_curated.tsv",
    "outputs/prjna419326/PRJNA419326_download_plan_curated.tsv",
    "outputs/prjna419326/PRJNA419326_download_plan_mirnaseq_only.tsv",
    "outputs/prjna419326/PRJNA419326_experimental_design_summary.txt",
    "outputs/prjna419326/PRJNA419326_experimental_design_summary.json",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def copy_file(src: Path, dst: Path, copied: list, missing: list):
    if not src.exists():
        missing.append(str(src))
        return

    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

    copied.append(
        {
            "source": str(src),
            "package_path": str(dst.relative_to(STAGING_DIR)),
            "size_bytes": dst.stat().st_size,
            "sha256": sha256(dst),
        }
    )


def main() -> int:
    if STAGING_DIR.exists():
        shutil.rmtree(STAGING_DIR)

    PACKAGE_DIR.mkdir(parents=True, exist_ok=True)
    STAGING_DIR.mkdir(parents=True, exist_ok=True)

    copied = []
    missing = []

    # README
    readme = STAGING_DIR / "README_V3_PRJNA419326.txt"
    readme.write_text(
        "\n".join(
            [
                "Pacchetto finale V3 - PRJNA419326",
                "=" * 80,
                "",
                "Questo pacchetto contiene script, report, ranking e file tabellari derivati",
                "dall'integrazione del dataset small RNA/miRNA PRJNA419326 nel framework",
                "AI/GNN per ncRNA di mais in risposta alla siccità.",
                "",
                "Nota importante:",
                "- I FASTQ grezzi e trimmed NON sono inclusi nel pacchetto perché molto grandi.",
                "- La loro validazione è documentata nei report inclusi.",
                "- Il download è riproducibile tramite i piani e gli script inclusi.",
                "",
                "Output principali:",
                "- report metodologico integrato V3",
                "- ranking sperimentale V3 top 500",
                "- candidati drought-responsive",
                "- mapping RNAcentral",
                "- summary QC/trimming/count matrix",
                "",
                f"Creato: {datetime.now().isoformat(timespec='seconds')}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    copied.append(
        {
            "source": "generated",
            "package_path": str(readme.relative_to(STAGING_DIR)),
            "size_bytes": readme.stat().st_size,
            "sha256": sha256(readme),
        }
    )

    # scripts
    for script in SCRIPTS:
        copy_file(BASE_DIR / script, STAGING_DIR / "scripts" / script, copied, missing)

    # outputs prjna
    for fname in OUTPUT_FILES:
        copy_file(OUT_DIR / fname, STAGING_DIR / "outputs/prjna419326" / fname, copied, missing)

    # metadata and design files
    for rel in METADATA_FILES + DESIGN_OUTPUTS:
        copy_file(BASE_DIR / rel, STAGING_DIR / rel, copied, missing)

    manifest = {
        "package_name": "maize_ncrna_gnn_v3_prjna419326_package",
        "created": datetime.now().isoformat(timespec="seconds"),
        "project_dir": str(BASE_DIR),
        "description": "V3 con integrazione dataset PRJNA419326 small RNA/miRNA siccità mais",
        "raw_fastq_included": False,
        "trimmed_fastq_included": False,
        "reason_large_files_excluded": "FASTQ grezzi/trimmed e file collassati sono grandi; il pacchetto include report, script e output interpretativi.",
        "files_copied": copied,
        "missing_files": missing,
        "status": "pass" if not missing else "warning_missing_files",
    }

    manifest_in_package = STAGING_DIR / "manifest_v3_prjna419326.json"
    manifest_in_package.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    # zip
    if ZIP_PATH.exists():
        ZIP_PATH.unlink()

    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(STAGING_DIR.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(PACKAGE_DIR))

    zip_sha = sha256(ZIP_PATH)

    manifest["zip_path"] = str(ZIP_PATH)
    manifest["zip_size_bytes"] = ZIP_PATH.stat().st_size
    manifest["zip_sha256"] = zip_sha

    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    summary = []
    summary.append("PACCHETTO FINALE V3 - PRJNA419326")
    summary.append("=" * 80)
    summary.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    summary.append(f"Status: {manifest['status']}")
    summary.append(f"File copiati: {len(copied)}")
    summary.append(f"File mancanti: {len(missing)}")
    summary.append(f"ZIP: {ZIP_PATH}")
    summary.append(f"ZIP size bytes: {ZIP_PATH.stat().st_size}")
    summary.append(f"ZIP SHA256: {zip_sha}")
    summary.append("")
    summary.append("Nota:")
    summary.append("I FASTQ grezzi/trimmed non sono inclusi per evitare un pacchetto enorme.")
    summary.append("La validazione e riproducibilità sono documentate nei report e negli script.")
    summary.append("")
    summary.append("Missing files:")
    if missing:
        summary.extend(missing)
    else:
        summary.append("Nessuno")

    SUMMARY_TXT.write_text("\n".join(summary) + "\n", encoding="utf-8")

    print("\n".join(summary))

    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
