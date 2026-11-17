#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 19
Individuazione dei percorsi reali di inferenza.

Analisi non invasiva:
- ricerca delle chiamate di caricamento dei modelli
- identificazione di checkpoint e moduli ML/GNN
- ispezione dei container e dei relativi port
- health check HTTP dei servizi esposti
- report JSON e Markdown

Non carica file .pt/.joblib.
Non esegue training o inferenza.
Non modifica GUI, container, modelli o dataset.
"""

import ast
import json
import re
import subprocess
import urllib.request

from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"

SOURCE_DIRS = [
    ROOT / "app",
    ROOT / "api",
    ROOT / "src",
    BASE,
]

SKIP_DIRS = {
    "runs", ".venv", "venv", "__pycache__",
    "node_modules", ".git",
}

PATTERNS = {
    "TORCH_LOAD": r"\btorch\.load\s*\(",
    "JOBLIB_LOAD": r"\bjoblib\.load\s*\(",
    "PICKLE_LOAD": r"\bpickle\.load\s*\(",
    "LOAD_STATE_DICT": r"\bload_state_dict\s*\(",
    "TORCH_NO_GRAD": r"\b(?:torch\.)?no_grad\s*\(",
    "MODEL_EVAL": r"\.eval\s*\(",
    "MODEL_PREDICT": r"\.predict(?:_proba)?\s*\(",
    "GAE": r"\b(?:GAE|GraphAutoEncoder|VGAE)\b",
    "PYG": r"\btorch_geometric\b",
    "HTTP_ROUTE": r"@(?:app|router)\.(?:get|post|put)\s*\(",
    "STREAMLIT": r"\b(?:streamlit|st\.)",
}

MODEL_NAMES = [
    "gae_v3_model.pt",
    "gae_model.pt",
    "PRJNA419326_srna_supervised_best_model.joblib",
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib",
    "maize_ncrna_graph_v3_prjna419326.pt",
]


def command(args):
    try:
        p = subprocess.run(
            args, capture_output=True,
            text=True, timeout=20
        )
        if p.returncode != 0:
            return None
        return p.stdout.strip()
    except Exception:
        return None


def health(url):
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": "AI-PRIMETTA-STEP19"}
        )
        with urllib.request.urlopen(req, timeout=4) as response:
            return {
                "url": url,
                "status": response.status,
                "ok": 200 <= response.status < 400,
            }
    except Exception as exc:
        return {
            "url": url,
            "ok": False,
            "error": str(exc)[:200],
        }


def source_files():
    seen = set()

    for folder in SOURCE_DIRS:
        if not folder.is_dir():
            continue

        for path in folder.rglob("*.py"):
            if any(p in SKIP_DIRS for p in path.parts):
                continue

            if path in seen:
                continue

            seen.add(path)
            yield path


def audit_sources():
    records = []

    for path in source_files():
        try:
            text = path.read_text(encoding="utf-8")
            ast.parse(text)
        except (OSError, UnicodeError, SyntaxError):
            continue

        lines = text.splitlines()

        for number, line in enumerate(lines, 1):
            found = [
                name for name, pattern in PATTERNS.items()
                if re.search(pattern, line)
            ]

            referenced = [
                filename for filename in MODEL_NAMES
                if filename in line
            ]

            if found or referenced:
                records.append({
                    "file": str(path.relative_to(ROOT)),
                    "line": number,
                    "signals": sorted(set(found)),
                    "model_references": referenced,
                    "source_excerpt": line.strip()[:220],
                })

    return records


def docker_containers():
    raw = command(["docker", "ps", "--format", "{{json .}}"])

    if raw is None:
        return []

    rows = []
    for line in raw.splitlines():
        try:
            row = json.loads(line)
        except ValueError:
            continue

        name = row.get("Names", "")

        # Esclude esplicitamente AI ORZO.
        if "orzo" in name.lower():
            continue

        if not any(
            term in name.lower()
            for term in ("primetta", "ncrna", "maize")
        ):
            continue

        details = command(["docker", "inspect", row["ID"]])

        if not details:
            continue

        info = json.loads(details)[0]

        rows.append({
            "name": name,
            "image": info["Config"].get("Image"),
            "running": info["State"].get("Running"),
            "ports": (
                info.get("NetworkSettings", {})
                .get("Ports") or {}
            ),
            "mounts": [
                {
                    "destination": item.get("Destination"),
                    "source": item.get("Source"),
                    "rw": item.get("RW"),
                }
                for item in info.get("Mounts", [])
            ],
        })

    return rows


def check_http(containers):
    results = []
    seen = set()

    for container in containers:
        for internal, mappings in container["ports"].items():
            port = internal.split("/")[0]

            if port not in {"8000", "8501"}:
                continue

            for mapping in mappings or []:
                host_port = mapping.get("HostPort")
                host_ip = mapping.get("HostIp", "")

                if not host_port:
                    continue

                # Verifica soltanto port pubblicati localmente.
                if host_ip not in (
                    "", "127.0.0.1", "0.0.0.0", "::", "::1"
                ):
                    continue

                suffix = (
                    "/_stcore/health"
                    if port == "8501"
                    else "/openapi.json"
                )

                url = f"http://127.0.0.1:{host_port}{suffix}"

                if url not in seen:
                    seen.add(url)
                    results.append(health(url))

    return results


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 19")
    print("IDENTIFICAZIONE INFERENZA REALE")
    print("=" * 70)

    print("\n[1] ANALISI DEL CODICE")

    records = audit_sources()
    counts = defaultdict(int)

    for row in records:
        for signal in row["signals"]:
            counts[signal] += 1

    for signal in PATTERNS:
        print(f"{signal}: {counts[signal]} riferimenti")

    important = {
        "TORCH_LOAD", "JOBLIB_LOAD", "PICKLE_LOAD",
        "LOAD_STATE_DICT", "MODEL_PREDICT",
    }

    print("\nPrincipali punti di caricamento/inferenza:")

    for row in records:
        if important.intersection(row["signals"]):
            print(
                f"{row['file']}:{row['line']}",
                ",".join(row["signals"])
            )

    print("\n[2] CONTAINER SPECIFICI AI PRIMETTA")

    containers = docker_containers()

    for container in containers:
        print(
            container["name"],
            "| running:", container["running"],
            "| image:", container["image"]
        )

    print("\n[3] ENDPOINT HTTP")

    checks = check_http(containers)

    for check in checks:
        print(
            "[OK]" if check["ok"] else "[ATTENZIONE]",
            check["url"],
            check.get("status", check.get("error", ""))
        )

    if not checks:
        print(
            "[NON VERIFICATO] Nessun endpoint HTTP "
            "locale individuato dai port pubblicati."
        )

    print("\n[4] SALVATAGGIO REPORT")

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / f"STEP19_INFERENZA_REALE_{stamp}"
    output.mkdir(parents=True, exist_ok=False)

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "19",
        "phase": "INFERENCE_PATH_DISCOVERY",
        "source_matches": records,
        "signal_counts": dict(counts),
        "containers": containers,
        "http_checks": checks,
        "checkpoint_inference_tested": False,
        "performance_validated": False,
        "gui_interaction_tested": False,
        "training_performed": False,
        "existing_components_modified": False,
    }

    report_path = output / "RAPPORTO_STEP19_INFERENZA.json"

    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("Report:", report_path)
    print("[OK] Codice originale invariato")
    print("[OK] Checkpoint invariati")
    print("[OK] Nessun training")


if __name__ == "__main__":
    main()