#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA - STEP 22.6

Risoluzione conflitto porta GUI originale V3.2.

Porta precedente: 8501 (occupata da AI ORZO)
Nuova porta:       8503

Crea un file Docker Compose aggiuntivo.
Non modifica il Compose originale.

Preserva:
- AI ORZO
- GUI AI PRIMETTA V3.3
- Monitor ML
- Modelli e dati scientifici

Comandi:

python3 scripts/step_22_6_avvio_gui_v32_porta_8503.py --start
python3 scripts/step_22_6_avvio_gui_v32_porta_8503.py --status
python3 scripts/step_22_6_avvio_gui_v32_porta_8503.py --stop
"""

import argparse
import copy
import json
import subprocess
import sys
import time
import urllib.request

from pathlib import Path


# ============================================================
# CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

ORIGINAL = ROOT / "docker-compose.yml"

OVERRIDE = (
    ROOT / "docker-compose.primetta-v32-8503.yaml"
)

SERVICE = "maize-ncrna-gnn"

OLD_PORT = 8501
NEW_PORT = 8503

URL = f"http://127.0.0.1:{NEW_PORT}"


# ============================================================
# ESECUZIONE COMANDI
# ============================================================

def execute(command, check=True):

    process = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False
    )

    if process.stdout.strip():
        print(
            process.stdout,
            flush=True
        )

    if process.stderr.strip():
        print(
            process.stderr,
            flush=True
        )

    if check and process.returncode != 0:
        raise RuntimeError(
            "Comando fallito: "
            + " ".join(command)
        )

    return process


# ============================================================
# LETTURA COMPOSE ORIGINALE
# ============================================================

def original_configuration():

    response = execute([
        "docker",
        "compose",
        "-f",
        str(ORIGINAL),
        "config",
        "--format",
        "json"
    ])

    config = json.loads(response.stdout)

    services = config.get(
        "services",
        {}
    )

    if SERVICE not in services:
        raise RuntimeError(
            "Servizio GUI V3.2 non trovato"
        )

    return services[SERVICE]


# ============================================================
# CREAZIONE COMPOSE AGGIUNTIVO
# ============================================================

def prepare_override():

    original_service = original_configuration()

    existing_ports = original_service.get(
        "ports",
        []
    )

    if not existing_ports:
        raise RuntimeError(
            "Nessuna porta pubblicata nel "
            "servizio originale"
        )

    new_ports = copy.deepcopy(
        existing_ports
    )

    changed = 0

    for port in new_ports:

        if str(
            port.get("published", "")
        ) == str(OLD_PORT):

            print(
                "[INFO] Porta precedente:",
                port
            )

            port["published"] = str(
                NEW_PORT
            )

            port["host_ip"] = (
                "127.0.0.1"
            )

            changed += 1

    if changed != 1:
        raise RuntimeError(
            "Attesa una sola porta pubblicata "
            f"{OLD_PORT}; trovate {changed}"
        )

    # !override sostituisce l'elenco porte,
    # senza sommarlo a quello originale.
    #
    # Le eventuali altre porte già presenti
    # vengono preservate.

    content = (
        "# AI PRIMETTA V3.2\n"
        "# Override della sola configurazione porte\n"
        "# Compose originale non modificato\n\n"
        "services:\n"
        f"  {SERVICE}:\n"
        "    ports: !override\n"
    )

    for port in new_ports:

        content += (
            "      - "
            + json.dumps(
                port,
                ensure_ascii=False
            )
            + "\n"
        )

    if OVERRIDE.exists():

        existing = OVERRIDE.read_text(
            encoding="utf-8"
        )

        if existing != content:
            raise RuntimeError(
                "Override già esistente con "
                "contenuto differente. "
                "Nessuna sovrascrittura effettuata."
            )

        print(
            "[OK] Override già presente"
        )

    else:

        OVERRIDE.write_text(
            content,
            encoding="utf-8"
        )

        print(
            "[OK] Override creato:",
            OVERRIDE
        )

    return new_ports


# ============================================================
# COMANDI DOCKER COMPOSE
# ============================================================

def compose_command():

    return [
        "docker",
        "compose",
        "-f",
        str(ORIGINAL),
        "-f",
        str(OVERRIDE)
    ]


# ============================================================
# VERIFICA CONFIGURAZIONE EFFETTIVA
# ============================================================

def verify_configuration():

    response = execute([
        *compose_command(),
        "config",
        "--format",
        "json"
    ])

    config = json.loads(
        response.stdout
    )

    ports = (
        config["services"][SERVICE]
        .get("ports", [])
    )

    correct = [
        port
        for port in ports
        if (
            str(port.get("published")) == str(NEW_PORT)
            and
            str(port.get("host_ip")) == "127.0.0.1"
        )
    ]

    old = [
        port
        for port in ports
        if str(
            port.get("published")
        ) == str(OLD_PORT)
    ]

    if len(correct) != 1 or old:
        raise RuntimeError(
            "La porta 8503 non risulta "
            "configurata correttamente"
        )

    print(
        "[OK] Porta GUI V3.2:",
        NEW_PORT
    )

    print(
        "[OK] Porta originale 8501 "
        "non pubblicata da questo servizio"
    )


# ============================================================
# CONTROLLO HTTP
# ============================================================

def check_http():

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    health = (
        URL + "/_stcore/health"
    )

    for attempt in range(10):

        try:

            with opener.open(
                health,
                timeout=3
            ) as response:

                if response.status == 200:

                    print(
                        "[OK] Streamlit risponde:",
                        health
                    )

                    return True

        except Exception:
            pass

        time.sleep(1)

    print(
        "[ATTENZIONE] Healthcheck GUI "
        "non ancora riuscito"
    )

    return False


# ============================================================
# AVVIO GUI V3.2
# ============================================================

def start():

    prepare_override()

    verify_configuration()

    print()
    print("[1] AVVIO GUI ORIGINALE V3.2")

    execute([
        *compose_command(),
        "up",
        "-d",
        "--no-deps",
        "--no-build",
        "--pull",
        "never",
        SERVICE
    ])

    print()
    print("[2] STATO CONTAINER")

    execute([
        *compose_command(),
        "ps",
        SERVICE
    ])

    print()
    print("[3] HEALTHCHECK")

    success = check_http()

    print()
    print("[4] LOG CONTAINER")

    execute([
        *compose_command(),
        "logs",
        "--tail",
        "60",
        SERVICE
    ], check=False)

    print()

    if success:

        print("=" * 65)

        print(
            "[OK] GUI ORIGINALE V3.2 "
            "RAGGIUNGIBILE"
        )

        print(
            "[URL]",
            URL
        )

        print("=" * 65)

    else:

        print(
            "[DA VERIFICARE] Il container "
            "potrebbe non aver avviato "
            "correttamente Streamlit."
        )

        sys.exit(2)


# ============================================================
# STATO
# ============================================================

def status():

    if not OVERRIDE.is_file():
        raise RuntimeError(
            "Override non presente. "
            "Eseguire prima --start."
        )

    print(
        "[1] STATO GUI V3.2"
    )

    execute([
        *compose_command(),
        "ps",
        SERVICE
    ])

    print()
    print(
        "[2] HEALTHCHECK"
    )

    check_http()

    print()
    print(
        "[3] LOG"
    )

    execute([
        *compose_command(),
        "logs",
        "--tail",
        "40",
        SERVICE
    ], check=False)


# ============================================================
# ARRESTO
# ============================================================

def stop():

    if not OVERRIDE.is_file():
        raise RuntimeError(
            "Override non presente"
        )

    execute([
        *compose_command(),
        "stop",
        SERVICE
    ])

    print(
        "[OK] Fermata solo "
        "la GUI originale AI PRIMETTA V3.2"
    )

    print(
        "[OK] AI ORZO non coinvolto"
    )

    print(
        "[OK] GUI V3.3 e monitor ML "
        "non coinvolti"
    )


# ============================================================
# PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    modes = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    modes.add_argument(
        "--start",
        action="store_true"
    )

    modes.add_argument(
        "--status",
        action="store_true"
    )

    modes.add_argument(
        "--stop",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 65)
    print(
        "AI PRIMETTA - STEP 22.6"
    )
    print(
        "GUI V3.2 SU PORTA 8503"
    )
    print("=" * 65)

    if args.start:
        start()

    elif args.status:
        status()

    elif args.stop:
        stop()


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__,
            str(exc),
            file=sys.stderr
        )

        sys.exit(1)