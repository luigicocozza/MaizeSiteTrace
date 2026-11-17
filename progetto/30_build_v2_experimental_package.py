#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import shutil
import zipfile
from datetime import datetime
from pathlib import Path


BASE_DIR = Path(".").resolve()

OUT_DIR = BASE_DIR / "outputs/v2_package"
PACKAGE_STAGING = OUT_DIR / "maize_ncrna_gnn_v2_userdata_gae"
ZIP_PATH = OUT_DIR / "maize_ncrna_gnn_v2_userdata_gae_package.zip"
CHECKSUM_PATH = OUT_DIR / "maize_ncrna_gnn_v2_userdata_gae_package.zip.sha256"
MANIFEST_JSON = OUT_DIR / "v2_package_manifest.json"
README_V2 = PACKAGE_STAGING / "README_V2_USERDATA_GAE.txt"


EXCLUDE_NAMES = {
    ".venv",
    "__pycache__",
    ".git",
    ".pytest_cache",
}

EXCLUDE_SUFFIXES = {
    ".pyc",
    ".pyo",
}

EXCLUDE_DIR_PARTS = {
    "outputs/v2_package",
    "logs",
    "backups",
}


ROOT_FILES = [
    "README.md",
    "requirements.txt",
    "requirements_runtime.txt",
    "requirements_docker_runtime.txt",
    "Dockerfile",
    "docker-compose.yml",
    "run_project.sh",
    "stop_project.sh",
    "run_api.sh",
    "run_gui.sh",
    "run_docker.sh",
    "stop_docker.sh",
    "validate_docker.sh",
]


ROOT_SCRIPT_PREFIXES = [
    "00",
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
    "13",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "23",
    "24",
    "25",
    "26",
    "27",
    "28",
    "29",
    "30",
    "31",
]


DIRS_TO_COPY = [
    "api",
    "app",
    "configs",
    "data/processed",
    "data/user_data",
    "docker",
    "docs",
    "models",
    "outputs/gae",
    "outputs/user_data",
    "outputs/e2e",
    "outputs/reports",
    "outputs/manuscript",
    "outputs/submission",
    "outputs/embeddings",
    "outputs/rag",
    "outputs/v2_guide",
    "src",
]


def should_skip(path: Path) -> bool:
    rel = path.relative_to(BASE_DIR).as_posix()

    if path.name in EXCLUDE_NAMES:
        return True

    if path.suffix in EXCLUDE_SUFFIXES:
        return True

    for part in EXCLUDE_DIR_PARTS:
        if rel == part or rel.startswith(part + "/"):
            return True

    return False


def copy_file(src: Path, dst: Path) -> None:
    if not src.exists() or not src.is_file():
        return

    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def copy_dir(src: Path, dst: Path) -> None:
    if not src.exists() or not src.is_dir():
        return

    for item in src.rglob("*"):
        if should_skip(item):
            continue

        rel = item.relative_to(src)
        target = dst / rel

        if item.is_dir():
            target.mkdir(parents=True, exist_ok=True)
        elif item.is_file():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item, target)


def collect_root_scripts() -> list[Path]:
    scripts = []

    for path in BASE_DIR.glob("*.py"):
        name = path.name

        for prefix in ROOT_SCRIPT_PREFIXES:
            if name.startswith(prefix):
                scripts.append(path)
                break

    return sorted(set(scripts))


def write_readme() -> None:
    text = f"""PACCHETTO V2 SPERIMENTALE - maize ncRNA GNN userdata GAE

Creato: {datetime.now().isoformat(timespec="seconds")}

Questa e' la seconda versione sperimentale del progetto.

Contiene:
1. Pipeline AI/GNN originale.
2. Modulo rete neurale Graph Autoencoder.
3. Pagina GUI per avviare il training GAE.
4. Pagina GUI per caricare dati del ricercatore.
5. Validazione file CSV/TSV/FASTA caricati.
6. Integrazione dati ricercatore con ranking ncRNA.
7. Pagina GUI Ranking dati ricercatore.
8. Knowledge base RAG sui dati caricati.
9. Chat dati ricercatore.
10. Guida operativa V2 in PDF e TXT.

Cartelle principali:
- app/
- api/
- configs/
- data/processed/
- data/user_data/
- models/graph_autoencoder/
- outputs/gae/
- outputs/user_data/

Avvio:
1. Entrare nella cartella del progetto.
2. Attivare ambiente Python.
3. Eseguire:
   ./run_project.sh

Interfaccia:
- GUI: http://localhost:8501
- API: http://localhost:8000/docs

Pagine V2:
- Rete neurale GAE
- Dati ricercatore
- Ranking dati ricercatore
- Chat dati ricercatore

Nota metodologica:
La rete neurale GAE e' sperimentale e non supervisionata.
Il ranking con dati utente e' tecnico-esplorativo.
La presenza di un ncRNA nei dati caricati non dimostra da sola un ruolo biologico nella risposta alla siccita'.
Sono necessari controllo qualita', replicati, analisi differenziale robusta e validazione sperimentale.
"""
    README_V2.write_text(text, encoding="utf-8")


def zip_directory(src_dir: Path, zip_path: Path) -> None:
    if zip_path.exists():
        zip_path.unlink()

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(src_dir.rglob("*")):
            if path.is_file():
                arcname = path.relative_to(src_dir.parent)
                zf.write(path, arcname)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()

    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)

    return h.hexdigest()


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if PACKAGE_STAGING.exists():
        shutil.rmtree(PACKAGE_STAGING)

    PACKAGE_STAGING.mkdir(parents=True, exist_ok=True)

    copied_files = []

    for filename in ROOT_FILES:
        src = BASE_DIR / filename
        dst = PACKAGE_STAGING / filename

        if src.exists():
            copy_file(src, dst)
            copied_files.append(str(dst.relative_to(PACKAGE_STAGING)))

    for script in collect_root_scripts():
        dst = PACKAGE_STAGING / script.name
        copy_file(script, dst)
        copied_files.append(str(dst.relative_to(PACKAGE_STAGING)))

    for dirname in DIRS_TO_COPY:
        src = BASE_DIR / dirname
        dst = PACKAGE_STAGING / dirname

        if src.exists():
            copy_dir(src, dst)

    write_readme()

    zip_directory(PACKAGE_STAGING, ZIP_PATH)

    checksum = sha256_file(ZIP_PATH)
    CHECKSUM_PATH.write_text(f"{checksum}  {ZIP_PATH.name}\n", encoding="utf-8")

    package_files = [p for p in PACKAGE_STAGING.rglob("*") if p.is_file()]

    required = [
        "25_train_graph_autoencoder.py",
        "26_validate_user_dataset.py",
        "27_link_user_data_to_ranking.py",
        "28_build_user_rag_knowledge.py",
        "29_user_rag_chat.py",
        "app/pages/01_Rete_neurale_GAE.py",
        "app/pages/02_Dati_ricercatore.py",
        "app/pages/03_Ranking_dati_ricercatore.py",
        "app/pages/04_Chat_dati_ricercatore.py",
        "models/graph_autoencoder/gae_model.pt",
        "outputs/gae/gae_summary.txt",
        "outputs/user_data",
        "outputs/v2_guide/guida_operativa_v2_userdata_gae.pdf",
        "outputs/v2_guide/guida_operativa_v2_userdata_gae.txt",
        "31_generate_v2_guide_pdf.py",
    ]

    missing = []

    for item in required:
        if not (PACKAGE_STAGING / item).exists():
            missing.append(item)

    manifest = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "package_dir": str(PACKAGE_STAGING),
        "zip_path": str(ZIP_PATH),
        "checksum_path": str(CHECKSUM_PATH),
        "sha256": checksum,
        "files_count": len(package_files),
        "zip_size_bytes": ZIP_PATH.stat().st_size,
        "required_missing": missing,
        "status": "pass" if not missing else "warning",
        "methodological_note": (
            "Pacchetto V2 sperimentale con Graph Autoencoder, dati ricercatore, "
            "ranking integrato e chat RAG utente."
        ),
    }

    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    print("V2 PACKAGE SUMMARY")
    print("=" * 80)
    print(f"Status: {manifest['status']}")
    print(f"Package dir: {PACKAGE_STAGING}")
    print(f"ZIP: {ZIP_PATH}")
    print(f"ZIP size bytes: {manifest['zip_size_bytes']}")
    print(f"Files count: {manifest['files_count']}")
    print(f"SHA256: {checksum}")
    print(f"Checksum file: {CHECKSUM_PATH}")
    print(f"Manifest: {MANIFEST_JSON}")

    if missing:
        print("")
        print("Missing required items:")
        for item in missing:
            print(f"- {item}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
