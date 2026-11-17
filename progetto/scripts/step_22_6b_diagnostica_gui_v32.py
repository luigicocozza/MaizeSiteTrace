#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA - STEP 22.6B
Diagnostica non distruttiva GUI V3.2.

Non avvia training, non arresta/ricrea container,
non modifica GUI, modelli né dati.

Scrive solo un rapporto nella cartella runs.
"""

import ast
import hashlib
import json
import subprocess
import sys
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT / "target_prediction_v33_vscode" / "runs"
)

ORIGINAL = ROOT / "docker-compose.yml"

OVERRIDE = (
    ROOT / "docker-compose.primetta-v32-8503.yaml"
)

SERVICE = "maize-ncrna-gnn"

PAGES = [
    ROOT / "app" / "streamlit_app.py",
    ROOT / "app" / "pages" / "01_Rete_neurale_GAE.py",
    ROOT / "app" / "pages" / "02_Dati_ricercatore.py",
    ROOT / "app" / "pages" / "03_Ranking_dati_ricercatore.py",
    ROOT / "app" / "pages" / "04_Chat_dati_ricercatore.py",
]

ARTIFACTS = [
    ROOT / "models" / "graph_autoencoder_v3"
    / "gae_v3_model.pt",

    ROOT / "data" / "processed" / "pyg"
    / "maize_ncrna_graph_v3_prjna419326.pt",

    ROOT / "models"
    / "srna_supervised_ml_v1_antileakage"
    / "PRJNA419326_srna_supervised_antileakage_best_model.joblib",
]


# ============================================================
# COMANDI DI SOLA LETTURA
# ============================================================

def shell(command, timeout=35):

    try:

        process = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False
        )

        return {
            "ok": process.returncode == 0,
            "code": process.returncode,
            "stdout": process.stdout.strip(),
            "stderr": process.stderr.strip()[-1200:]
        }

    except Exception as exc:

        return {
            "ok": False,
            "error": str(exc)
        }


# ============================================================
# CONTROLLO HTTP
# ============================================================

def http_check(port):

    url = (
        f"http://127.0.0.1:{port}/_stcore/health"
    )

    try:

        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            url,
            timeout=6
        ) as response:

            return {
                "ok": response.status == 200,
                "http_status": response.status,
                "response": response.read(100).decode(
                    "utf-8",
                    errors="replace"
                )
            }

    except Exception as exc:

        return {
            "ok": False,
            "error": str(exc)
        }


# ============================================================
# SHA256
# ============================================================

def sha256(path):

    if not path.is_file():
        return None

    digest = hashlib.sha256()

    with path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(1024 * 1024),
            b""
        ):
            digest.update(block)

    return digest.hexdigest()


# ============================================================
# INFORMAZIONI CONTAINER
# ============================================================

def inspect_container(identifier):

    if not identifier:
        return {"found": False}

    output = shell([
        "docker",
        "inspect",
        identifier
    ])

    if not output["ok"]:

        return {
            "found": False,
            "error": (
                output.get("stderr")
                or output.get("error")
            )
        }

    obj = json.loads(
        output["stdout"]
    )[0]

    return {
        "found": True,

        "name": obj.get(
            "Name", ""
        ).lstrip("/"),

        "running": obj.get(
            "State", {}
        ).get("Running", False),

        "image": obj.get(
            "Config", {}
        ).get("Image"),

        "published_ports": obj.get(
            "NetworkSettings", {}
        ).get("Ports"),

        "mounts": [
            {
                "source": m.get("Source"),
                "destination": m.get("Destination"),
                "read_only": not m.get("RW", True)
            }
            for m in obj.get("Mounts", [])
        ],
    }


# ============================================================
# PROGRAMMA PRINCIPALE
# ============================================================

def main():

    print("=" * 70)
    print("AI PRIMETTA - STEP 22.6B")
    print("DIAGNOSTICA GUI ORIGINALE V3.2")
    print("=" * 70)

    if not ORIGINAL.is_file():

        raise FileNotFoundError(
            ORIGINAL
        )

    compose = [
        "docker",
        "compose",
        "-f",
        str(ORIGINAL)
    ]

    if OVERRIDE.is_file():

        compose += [
            "-f",
            str(OVERRIDE)
        ]

    # --------------------------------------------------------
    # 1. STATO GUI V3.2
    # --------------------------------------------------------

    response = shell([
        *compose,
        "ps",
        "-q",
        SERVICE
    ])

    container_ids = response.get(
        "stdout", ""
    ).splitlines()

    container = (
        inspect_container(container_ids[0])
        if container_ids
        else {
            "found": False,
            "compose_error": response
        }
    )

    print("\n[1] SERVIZIO GUI V3.2")

    print(
        "Container:",
        container.get("name", "NON TROVATO")
    )

    print(
        "Running:",
        container.get("running", False)
    )

    print(
        "Immagine:",
        container.get("image", "N/D")
    )

    # --------------------------------------------------------
    # 2. HEALTHCHECK
    # --------------------------------------------------------

    print("\n[2] SALUTE DELLE DUE GUI")

    endpoints = {
        "v32_8503": http_check(8503),
        "v33_8502": http_check(8502)
    }

    for name, value in endpoints.items():

        print(
            name,
            "HTTP",
            value.get("http_status", "N/D"),
            "OK"
            if value["ok"]
            else value.get("error", "ERRORE")
        )

    # --------------------------------------------------------
    # 3. DIPENDENZE DEL CONTAINER V3.2
    # --------------------------------------------------------

    print(
        "\n[3] LIBRERIE DISPONIBILI "
        "NEL CONTAINER V3.2"
    )

    packages = {
        "status": "NON_VERIFICATO"
    }

    if container.get("running"):

        code = (
            "import importlib.util,json;"
            "names=["
            "'streamlit','torch','torch_geometric',"
            "'sklearn','joblib','networkx','fastapi'"
            "];"
            "print(json.dumps({"
            "k:bool(importlib.util.find_spec(k)) "
            "for k in names"
            "}))"
        )

        probe = shell([
            "docker",
            "exec",
            container_ids[0],
            "python3",
            "-B",
            "-c",
            code
        ])

        if probe["ok"]:

            try:

                packages = {
                    "status": "OK",
                    "modules": json.loads(
                        probe["stdout"]
                    )
                }

                for name, found in packages[
                    "modules"
                ].items():

                    print(
                        name,
                        "PRESENTE" if found else "ASSENTE"
                    )

            except json.JSONDecodeError:

                packages = {
                    "status": "ERRORE",
                    "raw": probe["stdout"][:250]
                }

        else:

            packages = {
                "status": "NON_VERIFICATO",
                "error": (
                    probe.get("stderr")
                    or probe.get("error")
                )
            }

    if packages["status"] != "OK":

        print(
            "Diagnostica pacchetti:",
            packages["status"]
        )

    # --------------------------------------------------------
    # 4. ANALISI DELLE CINQUE SCHERMATE
    # --------------------------------------------------------

    print("\n[4] FILE DELLE CINQUE SCHERMATE")

    pages = {}

    for path in PAGES:

        relative = str(
            path.relative_to(ROOT)
        )

        info = {
            "exists": path.is_file(),
            "sha256": sha256(path)
        }

        if path.is_file():

            try:

                source = path.read_text(
                    encoding="utf-8"
                )

                tree = ast.parse(
                    source,
                    filename=relative
                )

                buttons = []

                for node in ast.walk(tree):

                    if not isinstance(
                        node,
                        ast.Call
                    ):
                        continue

                    name = ast.unparse(
                        node.func
                    )

                    if name not in (
                        "st.button",
                        "st.form_submit_button"
                    ):
                        continue

                    title = "<dinamico>"

                    if (
                        node.args
                        and isinstance(
                            node.args[0],
                            ast.Constant
                        )
                        and isinstance(
                            node.args[0].value,
                            str
                        )
                    ):

                        title = node.args[0].value

                    buttons.append({
                        "line": node.lineno,
                        "label": title
                    })

                info["ui_buttons"] = buttons

            except (
                OSError,
                SyntaxError,
                UnicodeError
            ) as exc:

                info["parse_error"] = str(exc)

        pages[relative] = info

        print(
            relative,
            "PRESENTE" if info["exists"] else "MANCANTE",
            "| pulsanti:",
            len(info.get("ui_buttons", []))
        )

    # --------------------------------------------------------
    # 5. CHECKPOINT
    # --------------------------------------------------------

    print("\n[5] CHECKPOINT E GRAFO GAE V3")

    artifacts = {}

    for path in ARTIFACTS:

        relative = str(
            path.relative_to(ROOT)
        )

        artifacts[relative] = {
            "exists": path.is_file(),
            "bytes": (
                path.stat().st_size
                if path.is_file()
                else None
            )
        }

        print(
            relative,
            "PRESENTE"
            if path.is_file()
            else "MANCANTE"
        )

    # --------------------------------------------------------
    # 6. ALTRI CONTAINER
    # --------------------------------------------------------

    print("\n[6] APPLICAZIONI IN ESECUZIONE")

    other = shell([
        "docker",
        "ps",
        "--format",
        "{{.Names}} | {{.Status}} | {{.Ports}}"
    ])

    print(
        other.get(
            "stdout",
            "NON VERIFICABILE"
        )
    )

    # --------------------------------------------------------
    # 7. RAPPORTO
    # --------------------------------------------------------

    print("\n[7] CREAZIONE RAPPORTO DIAGNOSTICO")

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    output = RUNS / (
        "STEP22_6B_DIAGNOSTICA_V32_" + stamp
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {
        "project": "AI PRIMETTA",
        "step": "22.6B",
        "timestamp_utc": stamp,

        "container_v32": container,
        "http": endpoints,
        "container_packages": packages,

        "ui_source": pages,

        "gae_v3_artifacts": artifacts,

        "other_containers": other.get(
            "stdout", ""
        ),

        "graph_gui": {
            "nodes_reported": 384,
            "edges_reported": 1134
        },

        "graph_gae_v3": {
            "nodes_previous_test": 884,
            "edges_previous_test": 4782
        },

        "training_performed": False,

        "application_files_modified": False,

        "note": (
            "HTTP 200 e pagine visibili non dimostrano "
            "che training GAE, caricamento dati "
            "e Chat RAG funzionino."
        )
    }

    target = output / "RAPPORTO_STEP22_6B.json"

    target.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print(
        "[OK] Rapporto:",
        target
    )

    print(
        "[OK] Nessun training"
    )

    print(
        "[OK] GUI, checkpoint e dati non modificati"
    )


# ============================================================
# AVVIO
# ============================================================

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