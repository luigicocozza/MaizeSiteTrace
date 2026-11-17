#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 18.2
AUDIT DOCKER E AMBIENTE ML/GNN

Verifiche:
- disponibilita Docker e Docker Compose
- inventario container e immagini
- servizi Compose del progetto
- dipendenze dichiarate nei requirements
- disponibilita ML/GNN nei container gia attivi
- report finale JSON e Markdown

NON esegue:
- docker build / compose up
- installazione dipendenze
- modifiche ai container
- caricamento checkpoint
- training
"""

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"

MODULES = [
    "numpy",
    "pandas",
    "sklearn",
    "torch",
    "torch_geometric",
    "streamlit",
    "fastapi",
    "transformers",
    "joblib",
    "networkx",
]

COMPOSE_FILES = [
    ROOT / "docker-compose.yml",
    BASE / "compose.yaml",
]

REQUIREMENTS_FILES = [
    ROOT / "requirements.txt",
    ROOT / "requirements_runtime.txt",
    ROOT / "requirements_docker_runtime.txt",
]


def command(args, cwd=None, timeout=30):
    try:
        result = subprocess.run(
            args,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
        return {
            "ok": result.returncode == 0,
            "stdout": result.stdout.strip(),
            "error": result.stderr.strip()[:1500],
            "returncode": result.returncode,
        }
    except Exception as exc:
        return {
            "ok": False,
            "stdout": "",
            "error": str(exc),
            "returncode": None,
        }


def parse_json_lines(raw):
    rows = []
    for line in raw.splitlines():
        if line.strip():
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                raise RuntimeError(
                    "Output Docker JSON non valido"
                )
    return rows


def inspect_running_container(container_id):
    # Legge esclusivamente la disponibilita dei moduli.
    # Non importa torch e non carica file serializzati.
    code = (
        "import importlib.util,json,sys;"
        "names=json.loads(sys.argv[1]);"
        "result={};"
        "\nfor n in names:\n"
        " try:\n"
        "  result[n]=importlib.util.find_spec(n) is not None\n"
        " except (ImportError,ValueError,ModuleNotFoundError):\n"
        "  result[n]=False\n"
        "print(json.dumps({'python':sys.version.split()[0],"
        "'modules':result}))"
    )

    for interpreter in ("python3", "python"):
        result = command(
            [
                "docker", "exec", container_id,
                interpreter, "-c", code,
                json.dumps(MODULES)
            ],
            timeout=40,
        )

        if result["ok"]:
            try:
                return {
                    "ok": True,
                    "interpreter": interpreter,
                    **json.loads(result["stdout"]),
                }
            except json.JSONDecodeError:
                pass

    return {
        "ok": False,
        "reason": (
            "Python non disponibile o non eseguibile "
            "nel container attivo"
        ),
    }


def main():
    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 18.2")
    print("AUDIT DOCKER ML/GNN")
    print("=" * 72)

    print("\n[1] VERIFICA DOCKER")

    version = command(["docker", "--version"])
    compose_version = command(
        ["docker", "compose", "version"]
    )
    daemon = command(
        ["docker", "info", "--format", "{{.ServerVersion}}"]
    )

    print("Docker CLI:", version["stdout"] or version["error"])
    print(
        "Docker Compose:",
        compose_version["stdout"] or compose_version["error"]
    )
    print(
        "Docker Engine:",
        daemon["stdout"] if daemon["ok"]
        else "NON ACCESSIBILE"
    )

    print("\n[2] CONTAINER E IMMAGINI")

    containers = []
    images = []

    if daemon["ok"]:
        ps = command([
            "docker", "ps", "-a",
            "--format", "{{json .}}"
        ])

        img = command([
            "docker", "image", "ls",
            "--format", "{{json .}}"
        ])

        if ps["ok"]:
            containers = parse_json_lines(ps["stdout"])

        if img["ok"]:
            images = parse_json_lines(img["stdout"])

    print("Container:", len(containers))
    print("Immagini:", len(images))

    for item in containers:
        print(
            "CONTAINER:",
            item.get("Names", ""),
            "|",
            item.get("Status", ""),
            "|",
            item.get("Image", ""),
        )

    print("\n[3] SERVIZI COMPOSE")

    compose_results = []

    for path in COMPOSE_FILES:
        if not path.is_file():
            print("[MANCANTE]", path)
            continue

        result = command(
            [
                "docker", "compose",
                "-f", str(path),
                "config", "--services"
            ],
            cwd=path.parent,
        )

        services = (
            result["stdout"].splitlines()
            if result["ok"] else []
        )

        compose_results.append({
            "file": str(path),
            "valid": result["ok"],
            "services": services,
            "error": result["error"],
        })

        print(
            path.relative_to(ROOT),
            "|",
            "SERVIZI:",
            ", ".join(services)
            if services else "NON DETERMINATI"
        )

        if result["error"]:
            print("  Diagnostica:", result["error"][:300])

    print("\n[4] REQUIREMENTS")

    requirements = {}

    for path in REQUIREMENTS_FILES:
        if not path.is_file():
            continue

        packages = []

        for raw in path.read_text(
            encoding="utf-8"
        ).splitlines():

            value = raw.strip()

            if not value or value.startswith("#"):
                continue

            # Non registriamo URL o credenziali:
            # solo nomi di dipendenze dichiarate.
            if value.startswith(("-", "git+", "http")):
                continue

            package = value.split("#", 1)[0].strip()
            package = package.split(";", 1)[0].strip()

            if package:
                packages.append(package[:150])

        requirements[str(path.relative_to(ROOT))] = packages

        print(
            path.relative_to(ROOT),
            "|",
            len(packages),
            "dipendenze dichiarate"
        )

    print("\n[5] VERIFICA PYTHON NEI CONTAINER ATTIVI")

    runtime_checks = []

    for item in containers:
        status = item.get("Status", "").lower()

        if not status.startswith("up "):
            continue

        container_id = item.get("ID")

        if not container_id:
            continue

        result = inspect_running_container(container_id)

        runtime_checks.append({
            "container_name": item.get("Names", ""),
            "container_image": item.get("Image", ""),
            "result": result,
        })

        print()
        print("Container:", item.get("Names", ""))

        if not result["ok"]:
            print("[ATTENZIONE]", result["reason"])
            continue

        available = sum(
            bool(value)
            for value in result["modules"].values()
        )

        print("Python:", result["python"])
        print(
            "Librerie disponibili:",
            f"{available}/{len(MODULES)}"
        )

        for name, present in result["modules"].items():
            print(
                f"  {name}:",
                "OK" if present else "ASSENTE"
            )

    print("\n[6] VALUTAZIONE")

    candidate = []

    for check in runtime_checks:
        result = check["result"]
        modules = result.get("modules", {})

        if all(
            modules.get(name, False)
            for name in (
                "numpy",
                "torch",
                "torch_geometric",
            )
        ):
            candidate.append(check["container_name"])

    if candidate:
        status = "CANDIDATO_RUNTIME_ML_GNN_TROVATO"
        print("[OK] Ambiente ML/GNN candidato:", candidate)
    elif not daemon["ok"]:
        status = "DOCKER_NON_ACCESSIBILE"
        print("[ATTENZIONE] Docker Engine non accessibile")
    else:
        status = "RUNTIME_ML_GNN_NON_ANCORA_VERIFICATO"
        print(
            "[ATTENZIONE] Nessun container attivo "
            "ha confermato tutte le dipendenze ML/GNN."
        )
        print(
            "Le immagini Docker esistenti potrebbero "
            "comunque contenere l'ambiente richiesto."
        )

    print("\n[7] GENERAZIONE REPORT")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS / f"STEP18_2_AUDIT_DOCKER_{timestamp}"
    )
    output.mkdir(parents=True, exist_ok=False)

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "18.2",
        "status": status,
        "docker_cli": version["ok"],
        "docker_compose": compose_version["ok"],
        "docker_engine": daemon["ok"],
        "containers": containers,
        "images": images,
        "compose": compose_results,
        "requirements": requirements,
        "runtime_checks": runtime_checks,
        "candidate_containers": candidate,
        "inference_tested": False,
        "model_accuracy_validated": False,
        "docker_build_performed": False,
        "containers_started": False,
        "training_performed": False,
        "existing_components_modified": False,
    }

    report_path = output / "RAPPORTO_STEP18_2.json"

    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("\n" + "=" * 72)
    print("STEP 18.2 TERMINATO")
    print("=" * 72)
    print("Stato:", status)
    print("Report:", report_path)
    print("[OK] Nessuna installazione")
    print("[OK] Nessun training")
    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")


if __name__ == "__main__":
    main()