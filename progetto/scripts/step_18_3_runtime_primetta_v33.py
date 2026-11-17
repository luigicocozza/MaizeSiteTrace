#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 18.3

Verifica runtime ML/GNN specifico del progetto.

- Legge il Compose principale
- Identifica il servizio maize-ncrna-gnn
- Verifica immagini e container esistenti
- Non utilizza immagini di AI_ORZO_ITALIA
- Controlla le dipendenze senza caricare modelli
- Produce un report JSON

Un eventuale test dell'immagine usa un container
temporaneo con rete disabilitata e filesystem
in sola lettura.

Nessun training, build, pull o installazione.
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

MODULES = [
    "numpy",
    "pandas",
    "sklearn",
    "torch",
    "torch_geometric",
    "joblib",
    "networkx",
]

CHECK_CODE = """
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
    "python": sys.version.split()[0],
    "modules": result
}))
"""


def run(args, timeout=40):
    try:
        process = subprocess.run(
            args,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )

        return {
            "ok": process.returncode == 0,
            "stdout": process.stdout.strip(),
            "stderr": process.stderr.strip()[:1000],
        }
    except Exception as exc:
        return {
            "ok": False,
            "stdout": "",
            "stderr": str(exc),
        }


def docker_json(args):
    response = run(args)

    if not response["ok"]:
        return None, response["stderr"]

    try:
        return json.loads(response["stdout"]), None
    except json.JSONDecodeError:
        return None, "Risposta JSON non valida"


def candidate_containers():

    data, error = docker_json([
        "docker", "ps", "-a",
        "--format", "json",
    ])

    # docker ps --format json produces JSONL,
    # not a single JSON object.
    if data is not None:
        return data if isinstance(data, list) else [data]

    response = run([
        "docker", "ps", "-a",
        "--format", "{{json .}}",
    ])

    if not response["ok"]:
        return []

    rows = []

    for line in response["stdout"].splitlines():
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            continue

    return rows


def probe_image(image):

    for executable in ("python3", "python"):

        command = [
            "docker", "run",
            "--rm",
            "--pull", "never",
            "--network", "none",
            "--read-only",
            "--tmpfs", "/tmp:rw,nosuid,noexec,size=64m",
            "--cap-drop", "ALL",
            "--security-opt", "no-new-privileges",
            "--entrypoint", executable,
            image,
            "-B",
            "-c",
            CHECK_CODE,
            json.dumps(MODULES),
        ]

        result = run(command, timeout=90)

        if result["ok"]:
            try:
                data = json.loads(result["stdout"])
                return {
                    "ok": True,
                    "interpreter": executable,
                    **data,
                }
            except json.JSONDecodeError:
                pass

    return {
        "ok": False,
        "reason": (
            "Impossibile verificare Python "
            "nell'immagine Docker locale"
        ),
    }


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 — STEP 18.3")
    print("VERIFICA RUNTIME ML/GNN")
    print("=" * 70)

    if not COMPOSE.is_file():
        raise SystemExit(
            "[ERRORE] Compose principale mancante"
        )

    RUNS.mkdir(parents=True, exist_ok=True)

    print("\n[1] CONFIGURAZIONE COMPOSE")

    config, error = docker_json([
        "docker", "compose",
        "-f", str(COMPOSE),
        "config", "--format", "json",
    ])

    if config is None:
        raise SystemExit(
            "[ERRORE] Compose non leggibile: "
            + str(error)
        )

    services = config.get("services", {})

    print("Servizi:", ", ".join(sorted(services)))

    if SERVICE not in services:
        raise SystemExit(
            "[ERRORE] Servizio ML/GNN non trovato"
        )

    service = services[SERVICE]

    configured_image = service.get("image")
    build_present = bool(service.get("build"))

    print("Servizio:", SERVICE)
    print("Immagine dichiarata:", configured_image)
    print("Configurazione build:", build_present)

    print("\n[2] CONTAINER ASSOCIATI AL SERVIZIO")

    containers = candidate_containers()

    matched = []

    for container in containers:

        container_id = container.get("ID", "")

        if not container_id:
            continue

        inspection, error = docker_json([
            "docker", "inspect",
            container_id,
        ])

        if not inspection:
            continue

        info = inspection[0]

        labels = (
            info.get("Config", {})
            .get("Labels") or {}
        )

        name = info.get("Name", "").lstrip("/")
        service_name = labels.get(
            "com.docker.compose.service", ""
        )

        if (
            service_name != SERVICE
            and SERVICE not in name
        ):
            continue

        image = info.get("Image", "")

        matched.append({
            "name": name,
            "id": container_id,
            "image_id": image,
            "running": info.get(
                "State", {}
            ).get("Running", False),
        })

        print(
            "Container:", name,
            "| Running:",
            matched[-1]["running"]
        )

    print("\n[3] INDIVIDUAZIONE IMMAGINE LOCALE")

    image_candidates = []

    if configured_image:

        result = run([
            "docker", "image", "inspect",
            configured_image,
        ])

        if result["ok"]:
            image_candidates.append(configured_image)

    for item in matched:

        image_id = item["image_id"]

        result = run([
            "docker", "image", "inspect",
            image_id,
        ])

        if result["ok"]:
            image_candidates.append(image_id)

    image_candidates = list(
        dict.fromkeys(image_candidates)
    )

    print(
        "Immagini riconducibili al servizio:",
        len(image_candidates)
    )

    if not image_candidates:
        print(
            "[ATTENZIONE] Nessuna immagine locale "
            "riconducibile con certezza al servizio."
        )

    print("\n[4] VERIFICA DIPENDENZE ML/GNN")

    checks = []

    # Evita la verifica di immagini generiche:
    # si usano solo immagini legate al Compose
    # o ai container del servizio AI PRIMETTA.
    for image in image_candidates:

        print("\nImmagine:", image)

        result = probe_image(image)

        checks.append({
            "image": image,
            "result": result,
        })

        if not result["ok"]:
            print("[ATTENZIONE]", result["reason"])
            continue

        print("Python:", result["python"])

        for name, present in result["modules"].items():
            print(
                f"  {name}:",
                "OK" if present else "ASSENTE"
            )

    print("\n[5] VALUTAZIONE")

    usable = []

    for check in checks:
        result = check["result"]
        modules = result.get("modules", {})

        if all(
            modules.get(name, False)
            for name in (
                "numpy", "torch", "torch_geometric"
            )
        ):
            usable.append(check["image"])

    if usable:
        status = "RUNTIME_PRIMETTA_ML_GNN_CANDIDATO"
        print(
            "[OK] Immagine AI PRIMETTA con "
            "dipendenze ML/GNN disponibili"
        )
    elif image_candidates:
        status = "RUNTIME_PRIMETTA_DA_COMPLETARE"
        print(
            "[ATTENZIONE] Immagine presente, "
            "ma runtime ML/GNN non confermato"
        )
    else:
        status = "IMMAGINE_ML_GNN_NON_IDENTIFICATA"
        print(
            "[ATTENZIONE] Identificare o "
            "ricostruire l'immagine dal Compose"
        )

    print("\n[6] REPORT")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP18_3_RUNTIME_PRIMETTA_{timestamp}"
    )

    output.mkdir(parents=True, exist_ok=False)

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "18.3",
        "status": status,
        "compose": str(COMPOSE),
        "service": SERVICE,
        "configured_image": configured_image,
        "build_present": build_present,
        "matched_containers": matched,
        "image_candidates": image_candidates,
        "checks": checks,
        "usable_images": usable,
        "model_inference_tested": False,
        "model_accuracy_validated": False,
        "training_performed": False,
        "image_build_performed": False,
        "packages_installed": False,
        "persistent_containers_modified": False,
        "gui_modified": False,
    }

    report_path = output / "RAPPORTO_STEP18_3.json"

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print("\nSTEP 18.3 TERMINATO")
    print("Stato:", status)
    print("Report:", report_path)
    print("[OK] Nessun training")
    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Nessuna installazione")


if __name__ == "__main__":
    main()