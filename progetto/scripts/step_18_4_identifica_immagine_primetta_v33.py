#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 18.4

Identificazione non invasiva delle immagini Docker
del progetto e verifica del contesto di build.

Non esegue build, pull, run, installazioni o training.
Non modifica immagini, container, GUI o modelli.
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

RUNS = ROOT / "target_prediction_v33_vscode" / "runs"
COMPOSE = ROOT / "docker-compose.yml"
SERVICE = "maize-ncrna-gnn"


def execute(args, timeout=30):
    try:
        p = subprocess.run(
            args,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        return {
            "ok": p.returncode == 0,
            "stdout": p.stdout.strip(),
            "error": p.stderr.strip()[:1000],
        }
    except Exception as exc:
        return {
            "ok": False,
            "stdout": "",
            "error": str(exc),
        }


def read_json(args):
    result = execute(args)
    if not result["ok"]:
        return None, result["error"]
    try:
        return json.loads(result["stdout"]), ""
    except json.JSONDecodeError:
        return None, "JSON non valido"


def main():
    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 18.4")
    print("RICERCA IMMAGINE DOCKER ML/GNN")
    print("=" * 70)

    if not COMPOSE.is_file():
        raise RuntimeError("Compose principale mancante")

    config, error = read_json([
        "docker", "compose",
        "-f", str(COMPOSE),
        "config", "--format", "json",
    ])

    if config is None:
        raise RuntimeError(error)

    service = config.get("services", {}).get(SERVICE)

    if service is None:
        raise RuntimeError("Servizio non presente nel Compose")

    project = config.get("name", "")
    build = service.get("build") or {}

    print("\n[1] IDENTITA DEL PROGETTO")
    print("Progetto Compose:", project)
    print("Servizio:", SERVICE)
    print("Immagine esplicita:", service.get("image"))

    if isinstance(build, dict):
        print("Contesto build:", build.get("context"))
        print("Dockerfile:", build.get("dockerfile", "Dockerfile"))
    else:
        print("Configurazione build:", build)

    print("\n[2] INVENTARIO IMMAGINI LOCALI")

    result = execute([
        "docker", "image", "ls",
        "--no-trunc",
        "--format", "{{json .}}",
    ])

    if not result["ok"]:
        raise RuntimeError(result["error"])

    images = []

    for line in result["stdout"].splitlines():
        if line.strip():
            images.append(json.loads(line))

    print("Immagini trovate:", len(images))

    candidates = []

    for item in images:
        name = (
            str(item.get("Repository", ""))
            + ":"
            + str(item.get("Tag", ""))
        )

        low = name.lower()

        if "orzo" in low or "ollama" in low:
            continue

        if not any(
            word in low
            for word in (
                "primetta",
                "ncrna",
                "maize",
                "srna",
            )
        ):
            continue

        candidate = {
            "image": name,
            "id": item.get("ID", ""),
            "size": item.get("Size", ""),
            "reason": "NOME_COMPATIBILE",
        }

        candidates.append(candidate)

        print(
            "[CANDIDATA]",
            candidate["image"],
            "|",
            candidate["size"],
        )

    print("\n[3] IMMAGINI DEL PROGETTO COMPOSE")

    compose_images = execute([
        "docker", "compose",
        "-f", str(COMPOSE),
        "images", "--format", "json",
    ])

    if compose_images["ok"]:
        print(
            compose_images["stdout"]
            or "Nessuna immagine associata ai container del progetto"
        )
    else:
        print(
            "[INFO] Inventario Compose non disponibile:",
            compose_images["error"],
        )

    print("\n[4] VERIFICA ETICHETTE IMMAGINI")

    inspected = []

    for candidate in candidates:
        image_name = candidate["image"]

        details, error = read_json([
            "docker", "image",
            "inspect", image_name,
        ])

        if not details:
            print("[ATTENZIONE]", image_name, error)
            continue

        info = details[0]
        labels = info.get("Config", {}).get("Labels") or {}

        inspected.append({
            **candidate,
            "image_id": info.get("Id", ""),
            "labels": labels,
        })

        print("\nImmagine:", image_name)
        print("ID:", info.get("Id", ""))

        print(
            "Etichette:",
            json.dumps(
                labels,
                ensure_ascii=False,
            )[:600]
        )

    print("\n[5] VALUTAZIONE")

    if inspected:
        status = "IMMAGINI_CANDIDATE_IDENTIFICATE"
        print(
            "[OK] Immagini potenzialmente utilizzabili "
            "individuate: ",
            len(inspected)
        )
        print(
            "[NOTA] Il nome dell'immagine non prova "
            "ancora la compatibilita ML/GNN."
        )
    else:
        status = "NESSUNA_IMMAGINE_NOMINATA_IDENTIFICATA"
        print(
            "[ATTENZIONE] Nessuna immagine con nome "
            "riconducibile ad AI PRIMETTA."
        )
        print(
            "[NOTA] Potrebbero esistere immagini senza tag "
            "o con nomi differenti."
        )

    print("\n[6] REPORT")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP18_4_IMMAGINI_PRIMETTA_{timestamp}"
    )

    output.mkdir(parents=True, exist_ok=False)

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "18.4",
        "status": status,
        "compose_project": project,
        "compose_service": SERVICE,
        "build_configuration": build,
        "total_images": len(images),
        "candidates": inspected,
        "compose_images_result": compose_images,
        "runtime_verified": False,
        "model_inference_tested": False,
        "training_performed": False,
        "existing_components_modified": False,
    }

    report_path = output / "RAPPORTO_STEP18_4.json"

    report_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    print("Stato:", status)
    print("Report:", report_path)
    print("[OK] Nessuna modifica al sistema")


if __name__ == "__main__":
    main()