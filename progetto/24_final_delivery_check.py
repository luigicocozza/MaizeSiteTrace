#!/usr/bin/env python3
from __future__ import annotations

import json
import zipfile
from datetime import datetime
from pathlib import Path


BASE_DIR = Path(".").resolve()
ZIP_PATH = BASE_DIR / "outputs/final_package/maize_ncrna_gnn_final_package.zip"
OUT_JSON = BASE_DIR / "outputs/final_package/final_delivery_check.json"
OUT_TXT = BASE_DIR / "outputs/final_package/final_delivery_check_summary.txt"


REQUIRED_FILES = [
    "README.md",
    "Dockerfile",
    ".dockerignore",
    "docker-compose.yml",
    "docker/start_services.sh",
    "run_docker.sh",
    "stop_docker.sh",
    "validate_docker.sh",

    "app/streamlit_app.py",
    "api/main.py",
    "16_rag_chat_engine.py",
    "17_validate_end_to_end.py",
    "18_build_final_project_package.py",
    "19_build_software_user_guide_pdf.py",
    "20_build_final_technical_report.py",
    "21_build_scientific_manuscript_draft.py",
    "22_add_references_to_manuscript.py",
    "23_build_submission_materials.py",

    "data/processed/ranking/ncrna_exploratory_ranking.csv",
    "data/processed/similarity/ncrna_similarity_topk.csv",
    "data/processed/ncrna_features.csv",
    "data/processed/embeddings/ncrna_embeddings.csv",
    "data/processed/graph_nodes.csv",
    "data/processed/graph_edges.csv",
    "data/processed/rag/project_knowledge_base.jsonl",

    "outputs/e2e/e2e_validation_summary.txt",
    "outputs/manuals/guida_software_step_by_step.pdf",
    "outputs/reports/relazione_tecnica_finale.pdf",

    "docs/SOFTWARE_USER_GUIDE.md",
    "docs/RELAZIONE_TECNICA_FINALE.md",
    "docs/SCIENTIFIC_MANUSCRIPT_DRAFT.md",
    "docs/SCIENTIFIC_MANUSCRIPT_WITH_REFERENCES.md",
    "docs/SUBMISSION_CHECKLIST.md",
    "docs/HIGHLIGHTS.md",
    "docs/COVER_LETTER_DRAFT.md",

    "outputs/manuscript/scientific_manuscript_draft.pdf",
    "outputs/manuscript/scientific_manuscript_with_references.pdf",
    "outputs/manuscript/scientific_manuscript_metadata.json",
    "outputs/manuscript/references_used.json",
    "outputs/manuscript/submission_materials_metadata.json",
]


def main() -> int:
    if not ZIP_PATH.exists():
        raise FileNotFoundError(f"ZIP non trovato: {ZIP_PATH}")

    with zipfile.ZipFile(ZIP_PATH, "r") as z:
        names = set(z.namelist())

    results = []
    missing = []

    for item in REQUIRED_FILES:
        ok = item in names
        results.append({"file": item, "status": "OK" if ok else "MANCANTE"})
        if not ok:
            missing.append(item)

    status = "PASS" if not missing else "FAIL"

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "zip_path": str(ZIP_PATH),
        "zip_files_count": len(names),
        "required_files_count": len(REQUIRED_FILES),
        "missing_files_count": len(missing),
        "status": status,
        "missing": missing,
        "results": results,
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("FINAL DELIVERY CHECK")
    lines.append("=" * 80)
    lines.append(f"Status: {status}")
    lines.append(f"ZIP: {ZIP_PATH}")
    lines.append(f"File totali nello ZIP: {len(names)}")
    lines.append(f"File richiesti controllati: {len(REQUIRED_FILES)}")
    lines.append(f"File mancanti: {len(missing)}")
    lines.append("")

    if missing:
        lines.append("File mancanti:")
        for item in missing:
            lines.append(f"- {item}")
    else:
        lines.append("Tutti i file richiesti sono presenti.")

    lines.append("")
    lines.append("Dettaglio:")
    for r in results:
        lines.append(f"- {r['file']}: {r['status']}")

    OUT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
