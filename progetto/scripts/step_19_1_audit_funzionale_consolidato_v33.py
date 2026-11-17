#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 19.1
Audit consolidato dei percorsi di inferenza e runtime.

READ-ONLY:
- Non modifica GUI, backend o modelli
- Non carica checkpoint .pt/.joblib
- Non installa dipendenze
- Non avvia training
- Non ricostruisce immagini Docker
"""

import ast
import csv
import json
import os
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = (
    Path.home() / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"

CONTAINER = "ai_primetta_target_v33-prediction-1"
GNN_IMAGE = "ai_primetta_v33_gnn_cpu:18.6"

SKIP = {
    ".git", ".venv", "venv", "__pycache__",
    "node_modules", "runs", "delivery",
    "outputs", ".pytest_cache", "backups"
}

KEYWORDS = {
    "TORCH_LOAD": r"\btorch\s*\.\s*load\s*\(",
    "JOBLIB_LOAD": r"\bjoblib\s*\.\s*load\s*\(",
    "PICKLE_LOAD": r"\bpickle\s*\.\s*load\s*\(",
    "STATE_DICT": r"\bload_state_dict\s*\(",
    "PREDICT": r"\.\s*predict(?:_proba)?\s*\(",
    "EVAL": r"\.\s*eval\s*\(",
    "FORWARD": r"\bdef\s+forward\s*\(",
    "GNN": (
        r"torch_geometric|GCNConv|SAGEConv|"
        r"GATConv|HeteroConv|GraphAutoencoder|"
        r"GAE"
    ),
    "ML": (
        r"scikit.learn|sklearn|"
        r"RandomForest|LogisticRegression|"
        r"GradientBoosting|XGBoost"
    ),
    "STREAMLIT": r"\bstreamlit\b|\bst\s*\.",
    "FASTAPI": r"\bFastAPI\s*\(|@(?:app|router)\.",
    "API_CALL": (
        r"\brequests\s*\.\s*(?:get|post)\s*\(|"
        r"\bhttpx\s*\."
    ),
}

MODEL_TERMS = [
    "gae_model.pt",
    "gae_v3_model.pt",
    "PRJNA419326_srna_supervised_best_model",
    "PRJNA419326_srna_supervised_antileakage",
    "maize_ncrna_graph",
]

CONTAINER_CHECK = r'''
import json
import urllib.request

urls = [
    "http://127.0.0.1:8501/_stcore/health",
    "http://127.0.0.1:8000/openapi.json",
]

results = []

for url in urls:
    try:
        with urllib.request.urlopen(
            url, timeout=4
        ) as response:
            results.append({
                "url": url,
                "status": response.status,
                "ok": 200 <= response.status < 400
            })
    except Exception as exc:
        results.append({
            "url": url,
            "ok": False,
            "error": str(exc)[:200]
        })

print("AUDIT_JSON=" + json.dumps(results))
'''


def execute(args, timeout=30):
    try:
        p = subprocess.run(
            args,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False
        )
        return {
            "ok": p.returncode == 0,
            "code": p.returncode,
            "stdout": p.stdout.strip(),
            "stderr": p.stderr.strip()[:500]
        }
    except Exception as exc:
        return {
            "ok": False,
            "code": None,
            "stdout": "",
            "stderr": str(exc)
        }


def scan_python_files():
    results = []
    parse_errors = []

    for current, dirs, files in os.walk(ROOT):
        dirs[:] = [
            d for d in dirs
            if d not in SKIP
            and not d.startswith(".")
        ]

        for filename in files:
            if not filename.endswith(".py"):
                continue

            path = Path(current) / filename

            try:
                source = path.read_text(encoding="utf-8")
                tree = ast.parse(source, filename=str(path))
            except Exception as exc:
                parse_errors.append({
                    "file": str(path.relative_to(ROOT)),
                    "error": str(exc)
                })
                continue

            functions = [
                n.name
                for n in ast.walk(tree)
                if isinstance(
                    n, (ast.FunctionDef, ast.AsyncFunctionDef)
                )
            ]

            for lineno, line in enumerate(
                source.splitlines(), start=1
            ):
                stripped = line.strip()

                if not stripped:
                    continue

                signals = [
                    name
                    for name, pattern in KEYWORDS.items()
                    if re.search(pattern, line, re.I)
                ]

                model_refs = [
                    name for name in MODEL_TERMS
                    if name.lower() in line.lower()
                ]

                if not signals and not model_refs:
                    continue

                results.append({
                    "file": str(path.relative_to(ROOT)),
                    "line": lineno,
                    "signals": ";".join(signals),
                    "model_refs": ";".join(model_refs),
                    "excerpt": stripped[:240]
                })

            results.append({
                "file": str(path.relative_to(ROOT)),
                "line": 0,
                "signals": "FUNCTION_INVENTORY",
                "model_refs": "",
                "excerpt": ";".join(
                    name for name in functions
                    if any(
                        word in name.lower()
                        for word in (
                            "train", "infer", "predict",
                            "load", "evaluate", "model",
                            "graph", "rank"
                        )
                    )
                )[:1000]
            })

    return results, parse_errors


def inspect_container():
    response = execute([
        "docker", "inspect", CONTAINER
    ])

    if not response["ok"]:
        return {
            "exists": False,
            "error": response["stderr"]
        }

    info = json.loads(response["stdout"])[0]

    return {
        "exists": True,
        "running": bool(info["State"]["Running"]),
        "health": (
            info["State"].get("Health") or {}
        ).get("Status", "not_configured"),
        "image": info["Config"].get("Image"),
        "ports": (
            info.get("NetworkSettings", {})
            .get("Ports") or {}
        ),
        "mounts": [
            {
                "source": m.get("Source"),
                "destination": m.get("Destination"),
                "rw": m.get("RW")
            }
            for m in info.get("Mounts", [])
        ]
    }


def check_container_http():
    result = execute([
        "docker", "exec",
        CONTAINER,
        "python3", "-B", "-c",
        CONTAINER_CHECK
    ], timeout=20)

    if not result["ok"]:
        return {
            "executed": False,
            "error": result["stderr"]
        }

    for line in result["stdout"].splitlines():
        if line.startswith("AUDIT_JSON="):
            return {
                "executed": True,
                "endpoints": json.loads(
                    line.split("=", 1)[1]
                )
            }

    return {
        "executed": False,
        "error": "Risultato HTTP non identificato"
    }


def save_csv(path, rows):
    columns = [
        "file", "line", "signals",
        "model_refs", "excerpt"
    ]

    with path.open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 19.1")
    print("AUDIT FUNZIONALE CONSOLIDATO")
    print("=" * 70)

    print("\n[1] SCANSIONE CODICE COMPLETO")

    matches, errors = scan_python_files()
    counts = Counter()

    files = set()

    for row in matches:
        files.add(row["file"])

        for signal in row["signals"].split(";"):
            if signal:
                counts[signal] += 1

    print("File Python esaminati:", len(files))
    print("Errori di parsing:", len(errors))

    for name in KEYWORDS:
        print(f"{name}: {counts[name]}")

    print("\n[2] RIFERIMENTI AI MODELLI")

    relevant = [
        row for row in matches
        if row["model_refs"]
        or any(
            term in row["signals"]
            for term in (
                "TORCH_LOAD", "JOBLIB_LOAD",
                "STATE_DICT", "PREDICT", "FORWARD"
            )
        )
    ]

    print("Riferimenti rilevati:", len(relevant))

    for row in relevant[:45]:
        print(
            f"{row['file']}:{row['line']} "
            f"{row['signals']} "
            f"{row['model_refs']}"
        )

    print("\n[3] RUNTIME DOCKER GNN")

    image = execute([
        "docker", "image", "inspect",
        GNN_IMAGE
    ])

    print(
        "Immagine:",
        GNN_IMAGE,
        "PRESENTE" if image["ok"] else "ASSENTE"
    )

    print("\n[4] CONTAINER PREDIZIONE V3.3")

    container = inspect_container()

    print("Esistente:", container.get("exists"))
    print("In esecuzione:", container.get("running"))
    print("Health Docker:", container.get("health"))
    print("Porte:", container.get("ports", {}))

    print("\n[5] VERIFICA HTTP INTERNA")

    if container.get("running"):
        http = check_container_http()
    else:
        http = {
            "executed": False,
            "error": "Container non attivo"
        }

    for item in http.get("endpoints", []):
        print(
            "[OK]" if item["ok"] else "[NON RAGGIUNGIBILE]",
            item["url"],
            item.get("status", item.get("error"))
        )

    print("\n[6] VALUTAZIONE")

    blockers = []

    if errors:
        blockers.append("Errori di parsing Python")

    if not image["ok"]:
        blockers.append("Immagine GNN assente")

    if not container.get("running"):
        blockers.append("Container V3.3 non attivo")

    if not relevant:
        blockers.append(
            "Nessun percorso di inferenza identificato"
        )

    if not any(
        item.get("ok")
        for item in http.get("endpoints", [])
    ):
        blockers.append(
            "Nessun endpoint HTTP interno raggiungibile"
        )

    if blockers:
        status = "COLLAUDO_FUNZIONALE_DA_COMPLETARE"
    else:
        status = "PREREQUISITI_FUNZIONALI_RILEVATI"

    print("Stato:", status)
    print("Criticita:", len(blockers))

    for issue in blockers:
        print("[ATTENZIONE]", issue)

    print("\n[7] REPORT")

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP19_1_AUDIT_FUNZIONALE_{stamp}"
    )

    output.mkdir(parents=True, exist_ok=False)

    save_csv(
        output / "PERCORSI_INFERENZA_REALI.csv",
        matches
    )

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "19.1",
        "status": status,
        "python_files": len(files),
        "parse_errors": errors,
        "signal_counts": dict(counts),
        "relevant_locations": relevant,
        "runtime_image_present": image["ok"],
        "prediction_container": container,
        "http_checks": http,
        "blockers": blockers,
        "checkpoint_inference_tested": False,
        "predictive_performance_tested": False,
        "gui_interaction_tested": False,
        "training_performed": False,
        "existing_components_modified": False
    }

    report_path = (
        output / "RAPPORTO_STEP19_1.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print("\nSTEP 19.1 TERMINATO")
    print("Report:", report_path)
    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Checkpoint invariati")
    print("[OK] Nessun training")

    print(
        "\n[NOTA] La presenza di riferimenti nel "
        "codice e la risposta HTTP non sostituiscono "
        "il test di inferenza sui modelli reali."
    )


if __name__ == "__main__":
    main()