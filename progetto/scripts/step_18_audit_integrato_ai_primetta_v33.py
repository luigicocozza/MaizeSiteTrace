#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 18
AUDIT INTEGRATO ML/GNN, MOTORE E GUI

READ-ONLY: nessuna modifica ai componenti esistenti.

Esegue:
- Inventario dei moduli Python
- Analisi statica delle classi/funzioni
- Individuazione di modelli e checkpoint
- Individuazione di configurazioni e Docker
- Verifica sintattica degli script Python
- Verifica delle dipendenze installate
- Individuazione di collegamenti GUI/API/ML/GNN
- Controllo degli output scientifici STEP 17
- Produzione di report CSV/JSON/Markdown

NON:
- importa modelli non fidati
- carica checkpoint pickle/joblib/torch
- esegue training
- avvia Docker
- modifica GUI, motore, dataset o modelli
"""

import ast
import csv
import importlib.util
import json
import platform
import sys

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"

EXCLUDED = {
    ".git", ".venv", "venv",
    "__pycache__", "node_modules",
    "runs", "backups", ".mypy_cache"
}

MODEL_EXTENSIONS = {
    ".pt", ".pth", ".ckpt", ".h5",
    ".keras", ".onnx", ".joblib",
    ".pkl", ".safetensors"
}

CONFIG_EXTENSIONS = {
    ".yaml", ".yml", ".toml",
    ".ini", ".json", ".env"
}

MAX_FILE_MB = 5

MODULES = [
    "numpy", "pandas", "sklearn",
    "torch", "torch_geometric",
    "streamlit", "fastapi",
    "transformers", "joblib",
    "networkx"
]

KEYWORDS = {
    "GUI": [
        "streamlit", "gradio",
        "tkinter", "dash"
    ],
    "ML": [
        "sklearn", "xgboost",
        "catboost", "randomforest",
        "predict(", "fit("
    ],
    "GNN": [
        "torch_geometric", "graphsage",
        "heteroconv", "gcnconv",
        "gatconv", "sageconv"
    ],
    "API": [
        "fastapi", "flask",
        "uvicorn", "api_url"
    ],
    "RAG": [
        "retriever", "embedding",
        "vectorstore", "faiss",
        "chromadb"
    ],
}


def read_small_text(path):
    if path.stat().st_size > MAX_FILE_MB * 1048576:
        return ""

    try:
        return path.read_text(
            encoding="utf-8"
        )
    except (UnicodeError, OSError):
        return ""


def walk_project():
    import os

    for current, dirs, files in os.walk(ROOT):
        dirs[:] = [
            d for d in dirs
            if d not in EXCLUDED
            and not d.startswith(".")
        ]

        for name in files:
            yield Path(current) / name


def relative(path):
    return str(path.relative_to(ROOT))


def detect_components(path, text):
    lower = text.lower()
    return sorted([
        name
        for name, terms in KEYWORDS.items()
        if any(term in lower for term in terms)
    ])


def audit_python(path):
    text = read_small_text(path)

    if not text:
        return {
            "file": relative(path),
            "syntax": "NON_ANALIZZATO",
            "error": "File non leggibile o troppo grande",
            "classes": "",
            "functions": "",
            "components": ""
        }

    try:
        tree = ast.parse(text, filename=str(path))
        syntax = "OK"
        error = ""
        classes = [
            node.name for node in tree.body
            if isinstance(node, ast.ClassDef)
        ]
        functions = [
            node.name for node in tree.body
            if isinstance(
                node, (ast.FunctionDef, ast.AsyncFunctionDef)
            )
        ]
    except SyntaxError as exc:
        syntax = "ERRORE"
        error = f"Riga {exc.lineno}: {exc.msg}"
        classes = []
        functions = []

    return {
        "file": relative(path),
        "syntax": syntax,
        "error": error,
        "classes": ";".join(classes[:40]),
        "functions": ";".join(functions[:50]),
        "components": ";".join(
            detect_components(path, text)
        )
    }


def write_csv(path, rows, columns):
    with path.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns,
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def latest_report(prefix, filename):
    matches = sorted(
        p / filename
        for p in RUNS.glob(prefix + "*")
        if p.is_dir() and (p / filename).is_file()
    )
    return matches[-1] if matches else None


def main():
    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 18")
    print("AUDIT INTEGRATO ML/GNN E MOTORE")
    print("=" * 72)

    if not ROOT.is_dir() or not BASE.is_dir():
        raise SystemExit(
            "[ERRORE] Repository del progetto non trovato"
        )

    print("\n[1] INVENTARIO REPOSITORY")

    python_files = []
    model_files = []
    config_files = []
    other_files = []

    for path in walk_project():
        suffix = path.suffix.lower()
        name = path.name.lower()

        if suffix == ".py":
            python_files.append(path)

        elif suffix in MODEL_EXTENSIONS:
            model_files.append(path)

        elif (
            suffix in CONFIG_EXTENSIONS
            or "dockerfile" in name
            or name.startswith("requirements")
            or "compose" in name
        ):
            config_files.append(path)

        else:
            other_files.append(path)

    print("Script Python:", len(python_files))
    print("Checkpoint/modelli:", len(model_files))
    print("Configurazioni:", len(config_files))

    print("\n[2] VERIFICA SINTATTICA")

    audits = [
        audit_python(path)
        for path in sorted(python_files)
    ]

    syntax_counts = Counter(
        row["syntax"] for row in audits
    )

    print("Sintassi OK:", syntax_counts["OK"])
    print("Errori:", syntax_counts["ERRORE"])
    print(
        "Non analizzati:",
        syntax_counts["NON_ANALIZZATO"]
    )

    for row in audits:
        if row["syntax"] == "ERRORE":
            print(
                "[ERRORE]",
                row["file"],
                row["error"]
            )

    print("\n[3] MODULI ML/GNN E COMPONENTI")

    component_counts = Counter()

    for row in audits:
        for name in row["components"].split(";"):
            if name:
                component_counts[name] += 1

    for name in KEYWORDS:
        print(
            f"{name}:",
            component_counts[name],
            "file con riferimenti"
        )

    print("\n[4] CHECKPOINT DISPONIBILI")

    models = []

    for path in sorted(model_files):
        info = {
            "file": relative(path),
            "size_mb": round(
                path.stat().st_size / 1048576, 4
            ),
            "extension": path.suffix.lower(),
            "loaded": False,
            "verified_performance": False
        }
        models.append(info)

        print(
            info["file"],
            "|",
            info["size_mb"],
            "MB"
        )

    if not models:
        print(
            "[ATTENZIONE] Nessun checkpoint "
            "identificato dalle estensioni note."
        )

    print("\n[5] DIPENDENZE PYTHON")

    dependencies = {}

    for name in MODULES:
        try:
            present = (
                importlib.util.find_spec(name)
                is not None
            )
        except (ImportError, ValueError):
            present = False

        dependencies[name] = present
        print(
            name,
            "DISPONIBILE" if present else "ASSENTE"
        )

    print("\n[6] CONTROLLI RISULTATI GENOMICI")

    expected_reports = {
        "17.10": (
            "STEP17_10_MAPPING_GENOMICO_",
            "RAPPORTO_STEP17_10.json"
        ),
        "17.11": (
            "STEP17_11_ANNOTAZIONE_GTF_",
            "RAPPORTO_STEP17_11.json"
        ),
        "17.13": (
            "STEP17_13_ANNOTAZIONE_TE_",
            "RAPPORTO_STEP17_13.json"
        ),
        "17.14": (
            "STEP17_14_SEQUENZE_RIPETUTE_",
            "RAPPORTO_STEP17_14.json"
        ),
        "17.15": (
            "STEP17_15_COPIE_GENOMICHE_",
            "RAPPORTO_STEP17_15.json"
        ),
        "17.16": (
            "STEP17_16_AUDIT_COPIE_",
            "RAPPORTO_STEP17_16.json"
        )
    }

    scientific = {}

    for step, (prefix, filename) in (
        expected_reports.items()
    ):
        path = latest_report(prefix, filename)

        if path:
            try:
                data = json.loads(
                    path.read_text(encoding="utf-8")
                )
                scientific[step] = {
                    "present": True,
                    "valid_json": True,
                    "path": str(path),
                    "status": data.get(
                        "stato", data.get("status", "")
                    )
                }
                print("[OK]", step)
            except (ValueError, OSError) as error:
                scientific[step] = {
                    "present": True,
                    "valid_json": False,
                    "error": str(error),
                }
                print("[ERRORE]", step)
        else:
            scientific[step] = {
                "present": False,
                "valid_json": False,
            }
            print("[MANCANTE]", step)

    print("\n[7] GENERAZIONE RAPPORTI")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP18_AUDIT_INTEGRATO_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    write_csv(
        output / "AUDIT_SCRIPT_PYTHON.csv",
        audits,
        [
            "file", "syntax", "error",
            "classes", "functions", "components"
        ]
    )

    write_csv(
        output / "INVENTARIO_MODELLI.csv",
        models,
        [
            "file", "size_mb", "extension",
            "loaded", "verified_performance"
        ]
    )

    configs = [
        {
            "file": relative(p),
            "size_bytes": p.stat().st_size,
        }
        for p in sorted(config_files)
    ]

    write_csv(
        output / "INVENTARIO_CONFIGURAZIONI.csv",
        configs,
        ["file", "size_bytes"]
    )

    issues = []

    if syntax_counts["ERRORE"]:
        issues.append(
            "Esistono script Python con errori sintattici"
        )

    if not models:
        issues.append(
            "Checkpoint ML/GNN non identificati automaticamente"
        )

    if not any(
        row["components"]
        for row in audits
    ):
        issues.append(
            "Collegamenti applicativi non individuati "
            "dall'analisi statica"
        )

    if not all(
        row["valid_json"]
        for row in scientific.values()
    ):
        issues.append(
            "Alcuni report genomici mancano "
            "o non sono leggibili"
        )

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "18",
        "status": "STATIC_AUDIT_COMPLETED",
        "timestamp_utc": timestamp,
        "python": sys.version,
        "platform": platform.platform(),
        "python_scripts": len(audits),
        "syntax_counts": dict(syntax_counts),
        "components": dict(component_counts),
        "models_found": len(models),
        "models": models,
        "dependencies": dependencies,
        "genomic_reports": scientific,
        "issues": issues,
        "model_inference_tested": False,
        "model_performance_validated": False,
        "gui_end_to_end_tested": False,
        "limitations": [
            "Audit statico: non esegue inferenza.",
            "I riferimenti ML/GNN nel codice "
            "non garantiscono funzionalita operative.",
            "I checkpoint non vengono caricati.",
            "Le prestazioni predittive devono "
            "essere verificate con test separati.",
            "Le dipendenze vengono controllate "
            "nell'interprete Python corrente.",
            "Le verifiche genomiche confermano "
            "la presenza dei report, non una "
            "validazione biologica indipendente."
        ],
        "gui_modified": False,
        "engine_modified": False,
        "models_modified": False,
        "training_performed": False
    }

    json_path = output / "RAPPORTO_STEP18.json"

    json_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    md = [
        "# AI PRIMETTA V3.3 — Audit STEP 18",
        "",
        f"- Script Python: {len(audits)}",
        f"- Sintassi OK: {syntax_counts['OK']}",
        f"- Errori sintattici: {syntax_counts['ERRORE']}",
        f"- Checkpoint individuati: {len(models)}",
        "",
        "## Componenti individuati",
        ""
    ]

    for name, count in component_counts.items():
        md.append(f"- {name}: {count} file")

    md.extend([
        "",
        "## Punti da verificare",
        ""
    ])

    if issues:
        md.extend(
            f"- {issue}" for issue in issues
        )
    else:
        md.append(
            "- Nessuna criticita rilevata "
            "nell'ambito dell'audit statico."
        )

    md.extend([
        "",
        "## Limitazioni",
        "",
        "L'audit non certifica l'accuratezza "
        "dei modelli, il corretto funzionamento "
        "della GUI o la validita biologica "
        "delle predizioni."
    ])

    (output / "REPORT_STEP18.md").write_text(
        "\n".join(md) + "\n",
        encoding="utf-8"
    )

    print("\n" + "=" * 72)
    print("STEP 18 - AUDIT STATICO COMPLETATO")
    print("=" * 72)
    print("Script Python:", len(audits))
    print("Checkpoint:", len(models))
    print("Criticita preliminari:", len(issues))
    print("Cartella:", output)
    print("Report:", json_path)
    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun training")

    print(
        "\n[NOTA] Il superamento dell'audit "
        "statico non equivale al collaudo "
        "funzionale del sistema."
    )


if __name__ == "__main__":
    main()