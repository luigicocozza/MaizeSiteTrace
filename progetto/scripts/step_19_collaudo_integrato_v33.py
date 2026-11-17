#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 19
Collaudo preliminare integrato e non invasivo.

Non carica checkpoint pickle/torch.
Non esegue training.
Non modifica GUI, modelli, codice o container.
Verifica le condizioni necessarie per procedere
all'inferenza reale in sicurezza.
"""

import ast
import hashlib
import json
import subprocess
import urllib.error
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home() / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"

IMAGE = "ai_primetta_v33_gnn_cpu:18.6"

MODEL_DIRS = [
    ROOT / "models",
    ROOT / "data" / "processed" / "pyg",
]

MODEL_EXT = {".pt", ".pth", ".joblib", ".pkl", ".onnx"}

CONTAINERS = [
    "ai_primetta_target_v33-prediction-1",
]

FUNCTION_WORDS = (
    "load", "predict", "infer",
    "evaluate", "forward", "checkpoint"
)


def execute(args, timeout=25):
    try:
        p = subprocess.run(
            args,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False
        )
        return {
            "ok": p.returncode == 0,
            "stdout": p.stdout.strip(),
            "stderr": p.stderr.strip()[:500],
        }
    except Exception as exc:
        return {
            "ok": False,
            "stdout": "",
            "stderr": str(exc),
        }


def checksum(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(4 * 1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def get_container(name):
    response = execute([
        "docker", "inspect", name
    ])

    if not response["ok"]:
        return {
            "name": name,
            "exists": False,
            "running": False,
            "error": response["stderr"],
        }

    info = json.loads(response["stdout"])[0]
    ports = (
        info.get("NetworkSettings", {})
        .get("Ports") or {}
    )

    return {
        "name": name,
        "exists": True,
        "running": info["State"]["Running"],
        "status": info["State"]["Status"],
        "ports": ports,
        "image_id": info["Image"],
    }


def health_check(url):
    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "AI-PRIMETTA-AUDIT"}
        )

        with urllib.request.urlopen(
            request, timeout=5
        ) as response:
            return {
                "url": url,
                "http_status": response.status,
                "ok": 200 <= response.status < 400,
            }

    except Exception as exc:
        return {
            "url": url,
            "ok": False,
            "error": str(exc)[:300],
        }


def find_model_files():
    found = []

    for directory in MODEL_DIRS:
        if not directory.is_dir():
            continue

        for path in sorted(directory.rglob("*")):
            if path.is_file() and path.suffix in MODEL_EXT:
                found.append({
                    "path": str(path.relative_to(ROOT)),
                    "bytes": path.stat().st_size,
                    "sha256": checksum(path),
                    "loaded": False,
                    "inference_tested": False,
                })

    return found


def find_inference_code():
    found = []

    roots = [
        ROOT / "app",
        ROOT / "api",
        ROOT / "src",
        BASE,
    ]

    seen = set()

    for directory in roots:
        if not directory.is_dir():
            continue

        for path in directory.rglob("*.py"):
            if "runs" in path.parts:
                continue

            if path in seen:
                continue

            seen.add(path)

            try:
                source = path.read_text(encoding="utf-8")
                tree = ast.parse(source)
            except (OSError, UnicodeError, SyntaxError):
                continue

            functions = []

            for node in ast.walk(tree):
                if isinstance(
                    node,
                    (ast.FunctionDef, ast.AsyncFunctionDef)
                ):
                    if any(
                        word in node.name.lower()
                        for word in FUNCTION_WORDS
                    ):
                        functions.append(node.name)

            if functions:
                found.append({
                    "file": str(path.relative_to(ROOT)),
                    "functions": sorted(set(functions)),
                })

    return found


def main():
    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 19")
    print("COLLAUDO PRELIMINARE INTEGRATO")
    print("=" * 70)

    RUNS.mkdir(parents=True, exist_ok=True)

    print("\n[1] RUNTIME DOCKER")

    image_check = execute([
        "docker", "image", "inspect", IMAGE
    ])

    if image_check["ok"]:
        image_info = json.loads(image_check["stdout"])[0]
        runtime = {
            "exists": True,
            "tag": IMAGE,
            "image_id": image_info["Id"],
        }
        print("[OK]", IMAGE)
    else:
        runtime = {
            "exists": False,
            "tag": IMAGE,
        }
        print("[ERRORE] Immagine GNN mancante")

    print("\n[2] CHECKPOINT REALI")

    models = find_model_files()

    for model in models:
        print(
            "[OK]", model["path"],
            "| SHA256:", model["sha256"][:16]
        )

    print("File individuati:", len(models))

    print("\n[3] CODICE DI INFERENZA")

    inference = find_inference_code()

    for item in inference:
        print(
            item["file"],
            "| funzioni:",
            ", ".join(item["functions"][:8])
        )

    print("File con funzioni candidate:", len(inference))

    print("\n[4] CONTAINER APPLICATIVI")

    containers = [
        get_container(name) for name in CONTAINERS
    ]

    for item in containers:
        print(
            item["name"],
            "| running:",
            item["running"]
        )

    print("\n[5] HEALTH CHECK HTTP")

    checks = []

    for container in containers:
        if not container["running"]:
            continue

        ports = container.get("ports", {})

        for internal_port, mappings in ports.items():
            if not mappings:
                continue

            port_number = internal_port.split("/")[0]

            if port_number not in {"8501", "8000"}:
                continue

            for mapping in mappings:
                host_port = mapping.get("HostPort")

                if not host_port:
                    continue

                suffix = (
                    "/_stcore/health"
                    if port_number == "8501"
                    else "/openapi.json"
                )

                url = (
                    f"http://127.0.0.1:{host_port}"
                    f"{suffix}"
                )

                check = health_check(url)
                checks.append(check)

                print(
                    "[OK]" if check["ok"] else "[ATTENZIONE]",
                    url,
                )

    print("\n[6] RAPPORTO")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    out = RUNS / (
        f"STEP19_COLLAUDO_INTEGRATO_{timestamp}"
    )

    out.mkdir(parents=True, exist_ok=False)

    blockers = []

    if not runtime["exists"]:
        blockers.append("Runtime GNN mancante")

    if not models:
        blockers.append("Checkpoint non individuati")

    if not inference:
        blockers.append(
            "Codice di inferenza non individuato automaticamente"
        )

    if not all(c["running"] for c in containers):
        blockers.append(
            "Container applicativo non in esecuzione"
        )

    if any(not c["ok"] for c in checks):
        blockers.append(
            "Uno o piu health check HTTP falliti"
        )

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "19",
        "status": (
            "PRELIMINARY_CHECKS_PASS"
            if not blockers
            else "PRELIMINARY_CHECKS_WITH_BLOCKERS"
        ),
        "runtime": runtime,
        "models": models,
        "inference_code_candidates": inference,
        "containers": containers,
        "http_checks": checks,
        "blockers": blockers,
        "checkpoint_loading_tested": False,
        "model_inference_tested": False,
        "model_performance_validated": False,
        "gui_interaction_tested": False,
        "training_performed": False,
        "existing_components_modified": False,
    }

    report_path = out / "RAPPORTO_STEP19.json"

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print()
    print("=" * 70)
    print("STEP 19 - VERIFICA PRELIMINARE TERMINATA")
    print("=" * 70)

    print("Checkpoint:", len(models))
    print("Moduli candidati:", len(inference))
    print("Blocchi rilevati:", len(blockers))
    print("Report:", report_path)

    print("\n[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli invariati")
    print("[OK] Nessun training")

    print(
        "\n[NOTA] Il rapporto preliminare non "
        "certifica ancora l'inferenza reale."
    )


if __name__ == "__main__":
    main()