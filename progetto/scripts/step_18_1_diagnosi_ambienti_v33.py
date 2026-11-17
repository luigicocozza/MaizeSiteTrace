#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 18.1

Diagnosi non invasiva:
- interprete Python corrente
- ambienti virtuali presenti
- dipendenze disponibili per interprete
- riferimenti a requirements e Docker
- report JSON con valutazione corretta

Non installa pacchetti.
Non carica checkpoint.
Non modifica GUI, modelli o dataset.
"""

import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

RUNS = ROOT / "target_prediction_v33_vscode" / "runs"

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

CANDIDATE_DIRS = [
    ROOT,
    ROOT / "target_prediction_v33_vscode",
    ROOT.parent,
    Path.home(),
]

ENV_NAMES = [
    ".venv",
    "venv",
    "env",
    ".env",
    "venv_ai",
    "venv_gnn",
]

CHECK_CODE = r"""
import importlib.util
import json
import sys

names = json.loads(sys.argv[1])

result = {}
for name in names:
    try:
        result[name] = (
            importlib.util.find_spec(name) is not None
        )
    except (ImportError, ValueError, ModuleNotFoundError):
        result[name] = False

print(json.dumps({
    "executable": sys.executable,
    "version": sys.version.split()[0],
    "modules": result,
}))
"""


def inspect_python(path):

    try:
        result = subprocess.run(
            [
                str(path),
                "-c",
                CHECK_CODE,
                json.dumps(MODULES),
            ],
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        )

        return json.loads(result.stdout)

    except Exception as error:
        return {
            "executable": str(path),
            "error": str(error),
            "modules": {},
        }


def collect_interpreters():

    candidates = [Path(sys.executable)]

    for directory in CANDIDATE_DIRS:
        for env in ENV_NAMES:
            python = directory / env / "bin" / "python"

            if python.is_file():
                candidates.append(python)

    unique = []
    seen = set()

    for candidate in candidates:
        try:
            key = str(candidate.resolve())
        except OSError:
            key = str(candidate)

        if key not in seen:
            seen.add(key)
            unique.append(candidate)

    return unique


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 18.1")
    print("DIAGNOSI AMBIENTI PYTHON")
    print("=" * 70)

    print("\n[1] INTERPRETE CORRENTE")
    print("Python:", sys.executable)
    print("Versione:", sys.version.split()[0])
    print("VIRTUAL_ENV:", os.environ.get("VIRTUAL_ENV", "NON ATTIVO"))

    print("\n[2] AMBIENTI PYTHON")

    environments = []

    for interpreter in collect_interpreters():

        info = inspect_python(interpreter)

        modules = info.get("modules", {})
        found = sum(bool(value) for value in modules.values())

        info["available_modules"] = found
        info["total_modules"] = len(MODULES)

        environments.append(info)

        print()
        print("Python:", info["executable"])

        if "error" in info:
            print("[ERRORE]", info["error"])
            continue

        print(
            "Moduli disponibili:",
            f"{found}/{len(MODULES)}"
        )

        for name in MODULES:
            status = (
                "OK" if modules.get(name) else "ASSENTE"
            )
            print(f"  {name}: {status}")

    print("\n[3] CONFIGURAZIONI DISPONIBILI")

    configs = []

    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue

        if any(
            excluded in path.parts
            for excluded in (
                ".git", ".venv", "venv",
                "__pycache__", "runs", "node_modules"
            )
        ):
            continue

        name = path.name.lower()

        if (
            name.startswith("requirements")
            or name in (
                "pyproject.toml",
                "environment.yml",
                "environment.yaml",
                "pipfile",
            )
            or "dockerfile" in name
            or (
                "compose" in name
                and path.suffix in (".yaml", ".yml")
            )
        ):
            configs.append(str(path.relative_to(ROOT)))

    print("File identificati:", len(configs))

    for filename in sorted(configs)[:40]:
        print(" ", filename)

    print("\n[4] VALUTAZIONE")

    usable = [
        env for env in environments
        if "error" not in env
        and all(
            env["modules"].get(name, False)
            for name in ("numpy", "torch", "torch_geometric")
        )
    ]

    if usable:
        status = "AMBIENTE_ML_GNN_CANDIDATO_TROVATO"
        print("[OK] Ambiente ML/GNN candidato individuato")
    else:
        status = "AMBIENTE_ML_GNN_DA_INDIVIDUARE"
        print(
            "[ATTENZIONE] Nessun interprete verificato "
            "dispone contemporaneamente di numpy, "
            "torch e torch_geometric."
        )
        print(
            "Controllare eventuali ambienti Docker "
            "o ambienti virtuali in percorsi differenti."
        )

    print(
        "\n[NOTA] La presenza dei moduli non "
        "certifica la compatibilita delle versioni "
        "o il corretto funzionamento dell'inferenza."
    )

    print("\n[5] RAPPORTO")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP18_1_DIAGNOSI_AMBIENTI_{timestamp}"
    )
    output.mkdir(parents=True, exist_ok=False)

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "18.1",
        "status": status,
        "current_python": sys.executable,
        "virtual_env": os.environ.get("VIRTUAL_ENV"),
        "environments": environments,
        "configuration_files": sorted(configs),
        "inference_tested": False,
        "docker_tested": False,
        "training_performed": False,
        "existing_components_modified": False,
    }

    report_path = output / "RAPPORTO_STEP18_1.json"

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print("Stato:", status)
    print("Report:", report_path)
    print("[OK] File originali invariati")
    print("[OK] Nessun training")


if __name__ == "__main__":
    main()