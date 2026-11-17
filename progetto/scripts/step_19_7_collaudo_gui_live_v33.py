#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 19.7

Collaudo del servizio GUI realmente in esecuzione.

- Controllo container Docker
- Verifica HTTP Streamlit
- Recupero risultati GAE V3 e ML antileakage
- Verifica integrità file locali
- Registrazione facoltativa prova utente
- Report JSON

Non esegue training, installazioni, AppTest
o modifiche alla GUI e ai checkpoint.
"""

import argparse
import hashlib
import json
import subprocess
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

RUNS = ROOT / "target_prediction_v33_vscode/runs"

CONTAINER = "ai_primetta_target_v33-prediction-1"

URL = "http://127.0.0.1:8502"

FILES = {
    "gui": ROOT / "target_prediction_v33_vscode/app_v33.py",
    "gae_v3": ROOT / "models/graph_autoencoder_v3/gae_v3_model.pt",
    "ml_antileakage": ROOT /
        "models/srna_supervised_ml_v1_antileakage/"
        "PRJNA419326_srna_supervised_antileakage_best_model.joblib",
}


def sha256(path):
    if not path.is_file():
        return None

    h = hashlib.sha256()

    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def latest_report(prefix, filename):
    paths = sorted(
        RUNS.glob(f"{prefix}_*/{filename}")
    )

    if not paths:
        return None

    return json.loads(
        paths[-1].read_text(encoding="utf-8")
    )


def http_check(path):
    try:
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            URL + path,
            timeout=10
        ) as r:
            body = r.read(200000)
            return {
                "ok": r.status == 200,
                "status": r.status,
                "content_type": r.headers.get(
                    "Content-Type", ""
                ),
                "bytes_read": len(body),
                "body": (
                    body.decode(
                        "utf-8",
                        errors="replace"
                    ).strip()
                    if path == "/_stcore/health"
                    else None
                )
            }

    except Exception as e:
        return {
            "ok": False,
            "error": str(e)
        }


def inspect_container():
    p = subprocess.run(
        ["docker", "inspect", CONTAINER],
        capture_output=True,
        text=True,
        check=False
    )

    if p.returncode != 0:
        return {
            "running": False,
            "error": p.stderr.strip()[:300]
        }

    info = json.loads(p.stdout)[0]

    ports = (
        info.get("NetworkSettings", {})
        .get("Ports") or {}
    )

    return {
        "running": info["State"]["Running"],
        "status": info["State"]["Status"],
        "image": info["Config"]["Image"],
        "ports": ports
    }


def ask(question):
    print("\n" + question)
    print("s = sì | n = no | x = non verificato")

    while True:
        value = input("Risposta: ").strip().lower()

        if value in ("s", "n", "x"):
            return value


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--manual",
        action="store_true",
        help="Registra anche il collaudo nel browser"
    )

    args = parser.parse_args()

    print("=" * 68)
    print("AI PRIMETTA V3.3 — STEP 19.7")
    print("COLLAUDO GUI LIVE")
    print("=" * 68)

    before = {
        name: sha256(path)
        for name, path in FILES.items()
    }

    print("\n[1] CONTAINER DOCKER")
    container = inspect_container()
    print(json.dumps(container, indent=2))

    print("\n[2] HEALTH CHECK")
    health = http_check("/_stcore/health")
    print(health)

    print("\n[3] PAGINA PRINCIPALE")
    page = http_check("/")
    print({
        "http_status": page.get("status"),
        "content_type": page.get("content_type"),
        "bytes_read": page.get("bytes_read")
    })

    print("\n[4] REPORT ML/GNN")

    gae = latest_report(
        "STEP19_3_INFERENZA_GAE_V3",
        "RAPPORTO_STEP19_3.json"
    )

    ml = latest_report(
        "STEP19_5_COMPATIBILITA_ML",
        "RAPPORTO_STEP19_5.json"
    )

    gae_ok = (
        gae is not None
        and gae.get("status")
        == "TECHNICAL_INFERENCE_PASS"
    )

    ml_ok = (
        ml is not None
        and ml.get("status")
        == "COMPATIBLE_ML_INFERENCE_PASS"
        and (
            ml.get("ml_inference") or {}
        ).get("historical_label_agreement") == 1.0
    )

    print("GAE V3:", "PASS" if gae_ok else "NON VERIFICATO")
    print("ML:", "PASS" if ml_ok else "NON VERIFICATO")

    manual = {}

    if args.manual:
        print("\n[5] PROVA NEL BROWSER")
        print("Apri:", URL)
        print(
            "Verifica usando solo dati già presenti."
        )

        manual["pagina"] = ask(
            "La GUI si apre e mostra i widget senza errori?"
        )

        manual["predizione"] = ask(
            "Hai selezionato una sRNA reale "
            "e ottenuto un risultato applicativo?"
        )

        manual["navigazione"] = ask(
            "Le sezioni previste dalla GUI "
            "si aprono senza errori?"
        )

        manual["risultati"] = ask(
            "Le tabelle e i risultati previsti "
            "sono visualizzabili correttamente?"
        )

    else:
        print(
            "\n[5] Prova manuale non richiesta"
        )

    print("\n[6] INTEGRITÀ FILE")

    after = {
        name: sha256(path)
        for name, path in FILES.items()
    }

    unchanged = before == after

    print("File locali invariati:", unchanged)

    http_ok = (
        container.get("running") is True
        and health.get("ok") is True
        and health.get("body") == "ok"
        and page.get("ok") is True
    )

    if not http_ok or not gae_ok or not ml_ok:
        status = "BLOCKED"

    elif args.manual and all(
        value == "s"
        for value in manual.values()
    ):
        status = (
            "GUI_LIVE_AND_USER_WORKFLOW_PASS_"
            "INTEGRATION_NOT_CERTIFIED"
        )

    else:
        status = (
            "RUNTIME_AND_HTTP_PASS_"
            "MANUAL_WORKFLOW_PENDING"
        )

    if not unchanged:
        status = "FAILED_LOCAL_HASH_CHANGED"

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    out = RUNS / f"STEP19_7_GUI_LIVE_{stamp}"
    out.mkdir(parents=True, exist_ok=False)

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "19.7",
        "status": status,
        "container": container,
        "health": health,
        "page": page,
        "gae_inference_pass": gae_ok,
        "ml_inference_pass": ml_ok,
        "manual_user_checks": manual,
        "sha256_before": before,
        "sha256_after": after,
        "local_files_unchanged": unchanged,
        "gui_to_ml_gnn_direct_link_proven": False,
        "independent_scientific_validation": False,
        "training_performed": False
    }

    report_path = out / "RAPPORTO_STEP19_7.json"

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print("\n" + "=" * 68)
    print("STEP 19.7 TERMINATO")
    print("=" * 68)
    print("Stato:", status)
    print("Report:", report_path)
    print("[OK] Nessun training")
    print("[OK] Nessuna modifica intenzionale ai file")

    if status.startswith(("BLOCKED", "FAILED")):
        raise SystemExit(1)


if __name__ == "__main__":
    main()