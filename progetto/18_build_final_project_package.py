#!/usr/bin/env python3
from __future__ import annotations

import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

BASE_DIR = Path(".").resolve()
DOCS_DIR = BASE_DIR / "docs"
OUT_DIR = BASE_DIR / "outputs/final_package"
ZIP_PATH = OUT_DIR / "maize_ncrna_gnn_final_package.zip"

EXCLUDE_DIRS = {".venv", "__pycache__", ".git", "backups", "logs"}
EXCLUDE_PREFIXES = {"data/raw", "outputs/final_package"}


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
        return obj if isinstance(obj, dict) else {}
    except Exception:
        return {}


def csv_rows(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        return max(0, len(path.read_text(encoding="utf-8", errors="replace").splitlines()) - 1)
    except Exception:
        return 0


def write_file(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def get_summary() -> dict:
    e2e = read_json(BASE_DIR / "outputs/e2e/e2e_validation_report.json")
    pyg = read_json(BASE_DIR / "outputs/pyg_dataset_report.json")
    ranking = read_json(BASE_DIR / "outputs/ranking/ncrna_exploratory_ranking_report.json")
    similarity = read_json(BASE_DIR / "outputs/similarity/ncrna_similarity_report.json")
    rag = read_json(BASE_DIR / "outputs/rag/rag_knowledge_base_report.json")

    return {
        "project": "maize_drought_ncrna_gnn",
        "validation_status": e2e.get("validation_status", "unknown"),
        "overall_pass": e2e.get("checks", {}).get("overall_pass", False),
        "n_ranked_ncrna": ranking.get("summary", {}).get(
            "n_ranked_ncrna",
            csv_rows(BASE_DIR / "data/processed/ranking/ncrna_exploratory_ranking.csv"),
        ),
        "n_similarity_rows": similarity.get("summary", {}).get(
            "topk_rows",
            csv_rows(BASE_DIR / "data/processed/similarity/ncrna_similarity_topk.csv"),
        ),
        "n_graph_nodes": pyg.get("summary", {}).get(
            "n_nodes",
            csv_rows(BASE_DIR / "data/processed/graph_nodes.csv"),
        ),
        "n_graph_edges": pyg.get("summary", {}).get(
            "n_edges",
            csv_rows(BASE_DIR / "data/processed/graph_edges.csv"),
        ),
        "rag_chunks": rag.get("summary", {}).get("n_chunks", "NA"),
        "timestamp_utc": now(),
    }


def make_docs(summary: dict) -> None:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)

    readme = "\n".join([
        "# Maize ncRNA GNN",
        "",
        "Framework AI/GNN per identificare e prioritizzare ncRNA coinvolti nella risposta alla siccità nel mais.",
        "",
        "## Stato progetto",
        "",
        f"- Validazione end-to-end: **{summary['validation_status']}**",
        f"- Overall pass: **{summary['overall_pass']}**",
        f"- ncRNA nel ranking: **{summary['n_ranked_ncrna']}**",
        f"- Righe similarità: **{summary['n_similarity_rows']}**",
        f"- Nodi grafo: **{summary['n_graph_nodes']}**",
        f"- Archi grafo: **{summary['n_graph_edges']}**",
        f"- Chunk RAG: **{summary['rag_chunks']}**",
        "",
        "## Avvio rapido",
        "",
        "cd ~/progetti/ncrna_rnacentral_ai",
        "source .venv/bin/activate",
        "./run_project.sh",
        "",
        "GUI: http://localhost:8501",
        "API: http://localhost:8000/docs",
        "",
        "## Stop",
        "",
        "./stop_project.sh",
        "",
        "## Nota metodologica",
        "",
        "Ranking, similarità ed embedding sono esplorativi e non supervisionati.",
        "Non dimostrano validazione biologica della risposta alla siccità.",
    ])

    operativo = "\n".join([
        "# README operativo",
        "",
        "## Comandi principali",
        "",
        "cd ~/progetti/ncrna_rnacentral_ai",
        "source .venv/bin/activate",
        "./run_project.sh",
        "python 17_validate_end_to_end.py --verbose",
        "cat outputs/e2e/e2e_validation_summary.txt",
        "./stop_project.sh",
    ])

    struttura = "\n".join([
        "# Struttura progetto",
        "",
        "app/                  GUI Streamlit",
        "api/                  Backend FastAPI",
        "configs/              Configurazioni",
        "data/processed/       Dataset processati",
        "outputs/              Report e output",
        "docs/                 Documentazione",
        "run_project.sh        Avvio completo",
        "stop_project.sh       Stop completo",
        "16_rag_chat_engine.py Motore Chat RAG",
        "17_validate_end_to_end.py Validazione end-to-end",
    ])

    api_ref = "\n".join([
        "# API reference",
        "",
        "| Endpoint | Metodo | Funzione |",
        "|---|---|---|",
        "| /health | GET | Stato API |",
        "| /summary | GET | Sintesi progetto |",
        "| /candidates | GET | Lista candidati |",
        "| /candidate/{rnacentral_id} | GET | Dettaglio candidato |",
        "| /candidate/{rnacentral_id}/similar | GET | Simili del candidato |",
        "| /rag/ask | POST | Chat RAG |",
        "| /files/inventory | GET | Inventario file |",
    ])

    gui_guide = "\n".join([
        "# GUI guide",
        "",
        "| Pagina | Funzione |",
        "|---|---|",
        "| Dashboard | Sintesi generale |",
        "| Ranking candidati | Ranking esplorativo |",
        "| Dettaglio candidato | Schede candidate, feature, simili, export |",
        "| Ricerca avanzata | Filtri avanzati |",
        "| Similarità | Similarità top-k |",
        "| Embedding | Embedding ncRNA |",
        "| Chat RAG | Chat sul progetto |",
        "| Report metodologico | Report completo |",
        "| Download | Export ZIP/CSV/Markdown |",
    ])

    commands = "\n".join([
        "# Comandi",
        "",
        "./run_project.sh",
        "./run_gui.sh",
        "./run_api.sh",
        "python 17_validate_end_to_end.py --verbose",
        "python 18_build_final_project_package.py",
        "./stop_project.sh",
    ])

    requirements = "\n".join([
        "pandas",
        "numpy",
        "PyYAML",
        "requests",
        "scikit-learn",
        "matplotlib",
        "streamlit",
        "fastapi",
        "uvicorn",
        "pydantic",
        "torch",
        "torch-geometric",
        "biopython",
    ])

    write_file(BASE_DIR / "README.md", readme)
    write_file(DOCS_DIR / "README_OPERATIVO.md", operativo)
    write_file(DOCS_DIR / "PROJECT_STRUCTURE.md", struttura)
    write_file(DOCS_DIR / "COMMANDS.md", commands)
    write_file(DOCS_DIR / "API_REFERENCE.md", api_ref)
    write_file(DOCS_DIR / "GUI_GUIDE.md", gui_guide)
    write_file(BASE_DIR / "requirements_runtime.txt", requirements)


def excluded(path: Path) -> bool:
    rel = path.relative_to(BASE_DIR).as_posix()
    parts = set(path.relative_to(BASE_DIR).parts)

    if parts & EXCLUDE_DIRS:
        return True

    for prefix in EXCLUDE_PREFIXES:
        if rel == prefix or rel.startswith(prefix + "/"):
            return True

    if path.suffix in {".pyc", ".pyo"}:
        return True

    return False


def collect_files() -> list[Path]:
    roots = ["app", "api", "configs", "data/processed", "outputs", "docs", "docker"]
    files: list[Path] = []

    for root in roots:
        p = BASE_DIR / root
        if p.exists():
            for item in p.rglob("*"):
                if item.is_file() and not excluded(item):
                    files.append(item)

    for pattern in ["*.py", "*.sh", "README.md", "requirements*.txt", "Dockerfile", ".dockerignore", "docker-compose.yml"]:
        for item in BASE_DIR.glob(pattern):
            if item.is_file() and not excluded(item):
                files.append(item)

    unique = {}
    for f in files:
        unique[f.relative_to(BASE_DIR).as_posix()] = f

    return [unique[k] for k in sorted(unique)]


def build_zip(files: list[Path]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if ZIP_PATH.exists():
        ZIP_PATH.unlink()

    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for f in files:
            z.write(f, arcname=f.relative_to(BASE_DIR).as_posix())


def write_manifest(summary: dict, files: list[Path]) -> dict:
    manifest = {
        "timestamp_utc": now(),
        "summary": summary,
        "zip_path": str(ZIP_PATH),
        "zip_size_bytes": ZIP_PATH.stat().st_size if ZIP_PATH.exists() else 0,
        "n_files": len(files),
        "files": [
            {
                "path": f.relative_to(BASE_DIR).as_posix(),
                "size_bytes": f.stat().st_size,
            }
            for f in files
        ],
        "methodological_note": "Ranking e similarità sono esplorativi e non validano risposta alla siccità.",
    }

    write_file(
        OUT_DIR / "final_project_manifest.json",
        json.dumps(manifest, indent=2, ensure_ascii=False),
    )

    return manifest


def write_summary(manifest: dict) -> None:
    s = manifest["summary"]
    zip_size_mb = manifest["zip_size_bytes"] / 1024 / 1024

    text = "\n".join([
        "FINAL PROJECT PACKAGE SUMMARY",
        "=" * 80,
        f"Project: {s['project']}",
        f"Validation status: {s['validation_status']}",
        f"Overall pass: {s['overall_pass']}",
        "",
        f"Ranked ncRNA: {s['n_ranked_ncrna']}",
        f"Similarity rows: {s['n_similarity_rows']}",
        f"Graph nodes: {s['n_graph_nodes']}",
        f"Graph edges: {s['n_graph_edges']}",
        f"RAG chunks: {s['rag_chunks']}",
        "",
        f"Files in package: {manifest['n_files']}",
        f"ZIP: {manifest['zip_path']}",
        f"ZIP size MB: {zip_size_mb:.3f}",
        "",
        "Main commands:",
        "./run_project.sh",
        "python 17_validate_end_to_end.py --verbose",
        "./stop_project.sh",
        "",
        "Methodological note:",
        manifest["methodological_note"],
    ])

    write_file(OUT_DIR / "final_project_summary.txt", text)
    print(text)


def main() -> int:
    summary = get_summary()
    make_docs(summary)
    files = collect_files()
    build_zip(files)
    manifest = write_manifest(summary, files)
    write_summary(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
