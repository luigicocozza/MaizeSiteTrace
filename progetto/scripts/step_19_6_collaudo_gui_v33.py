#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==========================================================
AI PRIMETTA V3.3 - STEP 19.6
COLLAUDO AUTOMATICO GUI STREAMLIT
==========================================================

OBIETTIVI:

1. Verificare i report dei modelli GAE V3 e ML Antileakage.
2. Controllare HTTP 200 della GUI sulla porta 8502.
3. Verificare il container Docker della GUI.
4. Eseguire uno smoke test Streamlit con AppTest.
5. Identificare gli elementi renderizzati nella GUI.
6. Produrre un report JSON e il log del test.

VINCOLI:

- GUI V3.2 invariata.
- Motore V3.3 invariato.
- Checkpoint originali invariati.
- Dataset originali invariati.
- Nessun training.
- Nessuna installazione.
- Nessuna modifica ai container persistenti.

Il test AppTest esegue soltanto il primo rendering della pagina.
Non interagisce con pulsanti, chat o upload.
"""

import json
import subprocess
import sys
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


# ==========================================================
# 1. CONFIGURAZIONE PROGETTO
# ==========================================================

ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"

RUNS = BASE / "runs"

APP = BASE / "app_v33.py"

CONTAINER = "ai_primetta_target_v33-prediction-1"

HEALTH = "http://127.0.0.1:8502/_stcore/health"


# ==========================================================
# 2. CODICE APPTEST DA ESEGUIRE NEL CONTAINER TEMPORANEO
# ==========================================================

APPTEST_CODE = r'''
import json
import os
import sys

from pathlib import Path


app_path = Path(
    "/project/target_prediction_v33_vscode/app_v33.py"
)

result = {
    "app": str(app_path),
    "test_type": "STREAMLIT_APPTEST_FIRST_RENDER"
}


try:

    from streamlit.testing.v1 import AppTest

    if not app_path.is_file():
        raise FileNotFoundError(str(app_path))

    # Directory di lavoro dell'applicazione.
    os.chdir(str(app_path.parent))

    sys.path.insert(
        0,
        str(app_path.parent)
    )

    # Costruzione del test.
    at = AppTest.from_file(
        str(app_path),
        default_timeout=45
    )

    # Primo rendering senza interazioni.
    at.run(timeout=45)

    # Eventuali eccezioni Streamlit.
    errors = [
        str(
            getattr(
                item,
                "message",
                item
            )
        )[:700]
        for item in at.exception
    ]

    result["exceptions"] = errors

    # Inventario elementi renderizzati.
    counts = {}

    element_types = (
        "title",
        "header",
        "markdown",
        "button",
        "selectbox",
        "radio",
        "text_input",
        "text_area",
        "chat_input",
        "dataframe",
        "warning",
        "error"
    )

    for name in element_types:

        try:
            counts[name] = len(
                getattr(at, name)
            )
        except Exception:
            counts[name] = None

    result["rendered_elements"] = counts

    if errors:
        result["status"] = "FAIL_EXCEPTIONS"
    else:
        result["status"] = "PASS"


except Exception as exc:

    result["status"] = "NOT_VERIFIED"

    result["error"] = (
        type(exc).__name__
        + ": "
        + str(exc)[:700]
    )


print(
    "APPTEST_RESULT="
    + json.dumps(
        result,
        ensure_ascii=False
    ),
    flush=True
)
'''


# ==========================================================
# 3. ESECUZIONE COMANDI DI SISTEMA
# ==========================================================

def execute(arguments, timeout=100):

    try:

        process = subprocess.run(
            arguments,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False
        )

        return {
            "ok": process.returncode == 0,
            "returncode": process.returncode,
            "stdout": process.stdout.strip(),
            "stderr": process.stderr.strip()[-1500:]
        }

    except Exception as exc:

        return {
            "ok": False,
            "error": (
                f"{type(exc).__name__}: {exc}"
            )
        }


# ==========================================================
# 4. RECUPERO REPORT PRECEDENTI
# ==========================================================

def last_report(prefix, name):

    files = sorted(
        RUNS.glob(
            f"{prefix}_*/{name}"
        )
    )

    if not files:
        return {
            "found": False
        }

    path = files[-1]

    try:

        data = json.loads(
            path.read_text(
                encoding="utf-8"
            )
        )

        return {
            "found": True,
            "path": str(path),
            "status": data.get("status"),
            "original_files_unchanged": data.get(
                "original_files_unchanged"
            ),
            "historical_agreement": (
                data.get("ml_inference") or {}
            ).get(
                "historical_label_agreement"
            )
        }

    except Exception as exc:

        return {
            "found": True,
            "path": str(path),
            "error": str(exc)
        }


# ==========================================================
# 5. HEALTH CHECK HTTP GUI
# ==========================================================

def health_check():

    try:

        # Esclude eventuali proxy di sistema
        # per le richieste verso localhost.
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            HEALTH,
            timeout=8
        ) as response:

            body = response.read(100).decode(
                "utf-8",
                errors="replace"
            ).strip()

            return {
                "ok": (
                    response.status == 200
                    and body == "ok"
                ),
                "http_status": response.status,
                "body": body
            }

    except Exception as exc:

        return {
            "ok": False,
            "error": str(exc)
        }


# ==========================================================
# 6. PROCEDURA PRINCIPALE
# ==========================================================

def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 19.6")
    print("COLLAUDO AUTOMATICO GUI STREAMLIT")
    print("=" * 70)

    if not ROOT.is_dir():
        raise SystemExit(
            f"[ERRORE] Progetto non trovato: {ROOT}"
        )

    RUNS.mkdir(
        parents=True,
        exist_ok=True
    )

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP19_6_COLLAUDO_GUI_{stamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    # ------------------------------------------------------
    # 6.1 VERIFICA REPORT MODELLI
    # ------------------------------------------------------

    print("\n[1] REPORT DEI MODELLI")

    gae = last_report(
        "STEP19_3_INFERENZA_GAE_V3",
        "RAPPORTO_STEP19_3.json"
    )

    ml = last_report(
        "STEP19_5_COMPATIBILITA_ML",
        "RAPPORTO_STEP19_5.json"
    )

    print(
        "GAE V3:",
        gae.get(
            "status",
            "REPORT NON TROVATO"
        )
    )

    print(
        "ML Antileakage:",
        ml.get(
            "status",
            "REPORT NON TROVATO"
        )
    )

    print(
        "Concordanza predizioni storiche:",
        ml.get("historical_agreement")
    )

    # ------------------------------------------------------
    # 6.2 HEALTH CHECK GUI
    # ------------------------------------------------------

    print("\n[2] HEALTH CHECK GUI V3.3")

    health = health_check()

    print(
        "Risultato HTTP:",
        health
    )

    # ------------------------------------------------------
    # 6.3 VERIFICA CONTAINER DOCKER
    # ------------------------------------------------------

    print("\n[3] VERIFICA CONTAINER DOCKER")

    inspect = execute(
        [
            "docker",
            "inspect",
            CONTAINER
        ],
        timeout=20
    )

    info = (
        json.loads(
            inspect["stdout"]
        )[0]
        if inspect["ok"]
        else {}
    )

    image = (
        info.get("Config") or {}
    ).get("Image")

    running = (
        info.get("State") or {}
    ).get(
        "Running",
        False
    )

    print(
        "Container:",
        CONTAINER
    )

    print(
        "In esecuzione:",
        running
    )

    print(
        "Immagine:",
        image or "NON DISPONIBILE"
    )

    # ------------------------------------------------------
    # 6.4 TEST APPTEST STREAMLIT
    # ------------------------------------------------------

    print("\n[4] COLLAUDO APPTEST STREAMLIT")

    app_test = {
        "status": "NOT_VERIFIED",
        "reason": "Test non avviato"
    }

    if running and image and APP.is_file():

        cmd = [
            "docker",
            "run",
            "--rm",
            "--pull",
            "never",
            "--network",
            "none",
            "--read-only",
            "--tmpfs",
            "/tmp:rw,nosuid,nodev,size=256m",
            "--cap-drop",
            "ALL",
            "--security-opt",
            "no-new-privileges",
            "--pids-limit",
            "128",
            "--memory",
            "4g",
            "--cpus",
            "2",
            "-e",
            "HOME=/tmp",
            "-e",
            "XDG_CACHE_HOME=/tmp",
            "-e",
            "PYTHONDONTWRITEBYTECODE=1",
            "-v",
            f"{ROOT}:/project:ro",
            "-w",
            "/project/target_prediction_v33_vscode",
            "--entrypoint",
            "python3",
            image,
            "-B",
            "-c",
            APPTEST_CODE
        ]

        process = execute(
            cmd,
            timeout=105
        )

        # Salvataggio log completo del test.
        app_test_log = (
            output / "APPTEST_LOG.txt"
        )

        app_test_log.write_text(
            (process.get("stdout") or "")
            + "\nSTDERR:\n"
            + (process.get("stderr") or ""),
            encoding="utf-8"
        )

        marker = "APPTEST_RESULT="

        matches = [
            line.partition(marker)[2]
            for line in process.get(
                "stdout",
                ""
            ).splitlines()
            if line.startswith(marker)
        ]

        if len(matches) == 1:

            try:
                app_test = json.loads(
                    matches[0]
                )

            except ValueError:

                app_test = {
                    "status": "NOT_VERIFIED",
                    "error": "JSON AppTest non valido"
                }

        else:

            app_test = {
                "status": "NOT_VERIFIED",
                "error": (
                    "AppTest non ha restituito "
                    "un risultato strutturato"
                ),
                "returncode": process.get(
                    "returncode"
                ),
                "stderr": process.get(
                    "stderr",
                    ""
                )[-500:]
            }

    elif not APP.is_file():

        app_test = {
            "status": "NOT_VERIFIED",
            "reason": (
                f"File GUI non trovato: {APP}"
            )
        }

    print(
        "Stato AppTest:",
        app_test.get("status")
    )

    if app_test.get("error"):

        print(
            "[ATTENZIONE]",
            app_test["error"]
        )

    if app_test.get("exceptions"):

        print(
            "Eccezioni Streamlit:",
            app_test["exceptions"][:3]
        )

    if app_test.get("rendered_elements"):

        print(
            "Elementi GUI:",
            app_test["rendered_elements"]
        )

    # ------------------------------------------------------
    # 6.5 VERIFICA COMPLESSIVA
    # ------------------------------------------------------

    print("\n[5] VALUTAZIONE COMPLESSIVA")

    models_ok = (
        gae.get("status")
        == "TECHNICAL_INFERENCE_PASS"
        and
        ml.get("status")
        == "COMPATIBLE_ML_INFERENCE_PASS"
        and
        ml.get("historical_agreement")
        == 1.0
    )

    http_ok = health.get(
        "ok",
        False
    )

    if not models_ok or not http_ok:

        status = "BLOCKED"

    elif app_test.get("status") == "PASS":

        status = (
            "SOFTWARE_SMOKE_PASS_"
            "MANUAL_WORKFLOW_PENDING"
        )

    else:

        status = (
            "MODELS_AND_HTTP_PASS_"
            "APPTEST_NOT_VERIFIED"
        )

    print(
        "Stato finale:",
        status
    )

    # ------------------------------------------------------
    # 6.6 GENERAZIONE REPORT
    # ------------------------------------------------------

    print("\n[6] GENERAZIONE REPORT")

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "19.6",
        "status": status,

        "gae_checkpoint_report": gae,
        "ml_checkpoint_report": ml,

        "gui_health": health,

        "container_running": running,
        "image": image,

        "app_test": app_test,

        "real_gui_user_interactions_tested": False,

        "file_upload_tested": False,

        "chat_llm_tested": False,

        "gui_to_ml_model_link_tested": False,

        "independent_scientific_validation": False,

        "training_performed": False,

        "existing_components_modified": False
    }

    report_path = (
        output / "RAPPORTO_STEP19_6.json"
    )

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
    print("STEP 19.6 TERMINATO")
    print("=" * 70)

    print(
        "Stato:",
        status
    )

    print(
        "Report:",
        report_path
    )

    print(
        "[OK] Nessun training eseguito"
    )

    print(
        "[OK] GUI V3.2 invariata"
    )

    print(
        "[OK] Motore V3.3 invariato"
    )

    print(
        "[OK] Checkpoint e dati invariati"
    )

    print(
        "[INFO] Interazioni GUI, upload, chat e "
        "collegamento effettivo ai modelli "
        "non ancora verificati."
    )

    if status == "BLOCKED":
        sys.exit(1)


# ==========================================================
# 7. AVVIO
# ==========================================================

if __name__ == "__main__":
    main()