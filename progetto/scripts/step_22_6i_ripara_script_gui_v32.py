#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.2, STEP 22.6I.

--audit   verifica script nel repository e nel container
--repair  verifica, monta i file Python mancanti in sola lettura,
          riavvia soltanto la V3.2 e prova la RAG
--status  verifica container ed endpoint
--rollback elimina i nuovi mount dalla configurazione attiva

Non esegue training, non modifica GUI, script scientifici o checkpoint.
"""

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

SERVICE = "maize-ncrna-gnn"

CONTAINER = "maize_ncrna_gnn"

IMAGE = "ai_primetta_v32_full_ai_sklearn18:22.6g"

BASE = [
    ROOT / "docker-compose.yml",
    ROOT / "docker-compose.primetta-v32-8503.yaml",
    ROOT / "docker-compose.primetta-v32-full-ai.yaml",
    ROOT / "docker-compose.primetta-v32-sklearn18.yaml",
]

PATCH = ROOT / "docker-compose.primetta-v32-scripts.yaml"

RUNS = ROOT / "target_prediction_v33_vscode" / "runs"

CRITICAL = {
    "16_rag_chat_engine.py",
    "25_train_graph_autoencoder.py",
}


# ============================================================
# 1. COMANDI DI SISTEMA
# ============================================================

def call(args, timeout=100, ok=True):

    proc = subprocess.run(
        list(map(str, args)),
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )

    if ok and proc.returncode:

        raise RuntimeError(
            "Comando fallito: "
            + " ".join(map(str, args))
            + "\n"
            + proc.stderr[-1500:]
        )

    return proc


# ============================================================
# 2. CONFIGURAZIONE DOCKER COMPOSE
# ============================================================

def compose(with_patch=False):

    paths = BASE + ([PATCH] if with_patch else [])

    for path in paths:

        if not path.is_file():
            raise FileNotFoundError(str(path))

    out = ["docker", "compose"]

    for path in paths:
        out.extend(["-f", str(path)])

    return out


def inspect_config(with_patch=False):

    result = call(
        compose(with_patch)
        + ["config", "--format", "json"]
    )

    obj = json.loads(result.stdout)

    svc = obj["services"][SERVICE]

    ports = svc.get("ports", [])

    if svc.get("image") != IMAGE:

        raise RuntimeError(
            "Immagine Docker diversa da quella collaudata"
        )

    if not any(
        str(p.get("published")) == "8503"
        and str(p.get("target")) == "8501"
        and p.get("host_ip") == "127.0.0.1"
        for p in ports
    ):

        raise RuntimeError(
            "Porta V3.2 127.0.0.1:8503 non configurata"
        )

    if any(
        str(p.get("published")) == "8501"
        for p in ports
    ):

        raise RuntimeError(
            "Conflitto porta AI ORZO"
        )

    return svc


# ============================================================
# 3. IMMAGINE ATTIVA
# ============================================================

def image_active():

    p = call(
        [
            "docker",
            "inspect",
            CONTAINER,
            "--format",
            "{{.Config.Image}}",
        ],
        ok=False,
    )

    return (
        p.stdout.strip()
        if p.returncode == 0
        else None
    )


# ============================================================
# 4. CHECKSUM FILE
# ============================================================

def sha(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for part in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(part)

    return h.hexdigest()


def ui_hashes():

    paths = [
        ROOT / "app" / "streamlit_app.py"
    ]

    paths.extend(
        sorted(
            (ROOT / "app" / "pages").glob("*.py")
        )
    )

    if (
        len(paths) < 5
        or any(not p.is_file() for p in paths)
    ):

        raise RuntimeError(
            "File GUI principali mancanti"
        )

    return {
        str(p.relative_to(ROOT)): sha(p)
        for p in paths
    }


# ============================================================
# 5. INDIVIDUAZIONE SCRIPT ORIGINALI
# ============================================================

def candidates():

    names = set(CRITICAL)

    pages = [
        ROOT / "app" / "streamlit_app.py",

        *sorted(
            (ROOT / "app" / "pages").glob("*.py")
        ),
    ]

    # Trova riferimenti espliciti agli script Python
    # presenti nelle pagine Streamlit.

    for p in pages:

        if p.is_file():

            source = p.read_text(
                encoding="utf-8"
            )

            names.update(
                re.findall(
                    r"(?<![\w/])([\w.-]+\.py)\b",
                    source,
                )
            )

    # Include gli script storici numerati
    # disponibili nella radice del repository.

    names.update(
        p.name
        for p in ROOT.glob("[0-9]*_*.py")
        if p.is_file()
    )

    names = {
        n
        for n in names
        if (
            n.endswith(".py")
            and "/" not in n
            and "\\" not in n
        )
    }

    if len(names) > 180:

        raise RuntimeError(
            "Numero di script inaspettatamente alto"
        )

    mapping = {}

    for name in sorted(names):

        found = [
            p
            for p in (
                ROOT / name,
                ROOT / "scripts" / name,
            )
            if p.is_file()
        ]

        if len(found) == 1:

            mapping[name] = found[0]

        elif len(found) > 1:

            if sha(found[0]) != sha(found[1]):

                raise RuntimeError(
                    "Script duplicato discordante: "
                    + name
                )

            mapping[name] = found[0]

    return names, mapping


# ============================================================
# 6. VERIFICA FILE NEL CONTAINER
# ============================================================

def container_files(names):

    code = (
        "import json,pathlib; "
        "names=json.loads("
        + repr(json.dumps(sorted(names)))
        + "); "
        "print(json.dumps({"
        "n:(pathlib.Path('/app')/n).is_file() "
        "for n in names"
        "}))"
    )

    proc = call([
        "docker",
        "exec",
        CONTAINER,
        "python3",
        "-B",
        "-c",
        code,
    ])

    return json.loads(proc.stdout)


# ============================================================
# 7. AUDIT SCRIPTS GUI
# ============================================================

def audit():

    inspect_config(False)

    if image_active() != IMAGE:

        raise RuntimeError(
            "La V3.2 non utilizza l'immagine "
            "22.6g collaudata"
        )

    names, mapping = candidates()

    in_container = container_files(names)

    missing = {
        name: mapping[name]
        for name in names
        if (
            not in_container.get(name)
            and name in mapping
        )
    }

    unavailable = sorted(
        name
        for name in names
        if (
            not in_container.get(name)
            and name not in mapping
        )
    )

    print("[IMMAGINE]", IMAGE)

    print(
        "[SCRIPT] Verificati:",
        len(names)
    )

    print(
        "[SCRIPT] Mancanti nel container, "
        "ma disponibili sul disco:",
        len(missing)
    )

    for name, source in sorted(
        missing.items()
    ):

        print(
            "  [MONTABILE]",
            name,
            "<=",
            source.relative_to(ROOT)
        )

    for name in unavailable:

        print(
            "  [ASSENTE ANCHE SUL DISCO]",
            name
        )

    return missing, unavailable


# ============================================================
# 8. HEALTHCHECK HTTP
# ============================================================

def health(port, path):

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    url = (
        f"http://127.0.0.1:{port}{path}"
    )

    try:

        with opener.open(
            url,
            timeout=4
        ) as response:

            return response.status == 200

    except Exception:

        return False


# ============================================================
# 9. TEST REALE CHAT RAG
# ============================================================

def rag_test():

    payload = {

        "question": (
            "Quali ncRNA caricati dal ricercatore "
            "sono presenti nel ranking e hanno "
            "FDR inferiore a 0.05?"
        ),

        "top_k": 5,
    }

    request = urllib.request.Request(

        "http://127.0.0.1:8000/rag/ask",

        data=json.dumps(
            payload
        ).encode("utf-8"),

        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        },

        method="POST",
    )

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    try:

        with opener.open(
            request,
            timeout=90
        ) as response:

            body = response.read(
                16000
            ).decode(
                "utf-8",
                errors="replace"
            )

            return {

                "status": response.status,

                "ok": response.status == 200,

                "body": body[:1800],
            }

    except urllib.error.HTTPError as exc:

        return {

            "status": exc.code,

            "ok": False,

            "body": exc.read(
                2500
            ).decode(
                "utf-8",
                "replace"
            ),
        }

    except Exception as exc:

        return {

            "ok": False,

            "error": repr(exc),
        }


# ============================================================
# 10. STATO APPLICAZIONI
# ============================================================

def status():

    print(
        "[CONTAINER]",
        image_active()
    )

    services = [

        (
            "PRIMETTA V3.2",
            8503,
            "/_stcore/health",
        ),

        (
            "API V3.2",
            8000,
            "/health",
        ),

        (
            "PRIMETTA V3.3",
            8502,
            "/_stcore/health",
        ),

        (
            "AI ORZO",
            8501,
            "/_stcore/health",
        ),
    ]

    for label, port, path in services:

        print(
            "[OK]"
            if health(port, path)
            else "[NON DISPONIBILE]",
            label
        )

    process = call(
        [
            "docker",
            "ps",
            "--format",
            (
                "table {{.Names}}\t"
                "{{.Status}}\t{{.Ports}}"
            ),
        ],
        ok=False,
    )

    print(process.stdout)


# ============================================================
# 11. AVVIO CONTAINER V3.2
# ============================================================

def start_container(with_patch):

    command = compose(with_patch) + [

        "up",

        "-d",

        "--no-deps",

        "--no-build",

        "--pull",
        "never",

        "--force-recreate",

        SERVICE,
    ]

    call(
        command,
        timeout=160
    )

    for _ in range(25):

        ready = (

            image_active() == IMAGE

            and health(
                8503,
                "/_stcore/health"
            )

            and health(
                8000,
                "/health"
            )
        )

        if ready:
            return

        time.sleep(2)

    raise RuntimeError(
        "Riavvio V3.2: GUI o API non rispondono"
    )


# ============================================================
# 12. RIPRISTINO
# ============================================================

def rollback():

    inspect_config(False)

    print(
        "[ROLLBACK] Ripristino servizio "
        "senza mount aggiuntivi"
    )

    start_container(False)

    status()


# ============================================================
# 13. RIPARAZIONE SCRIPTS GUI
# ============================================================

def repair():

    before = ui_hashes()

    missing, unavailable = audit()

    # Non ricostruire arbitrariamente il motore
    # Chat RAG se manca dal repository originale.

    if "16_rag_chat_engine.py" in unavailable:

        raise RuntimeError(
            "16_rag_chat_engine.py non e' "
            "disponibile neanche nel repository. "
            "Verificare il file originale."
        )

    if not missing:

        print(
            "[ATTENZIONE] Nessuno script "
            "da montare: il problema RAG "
            "ha un'altra causa"
        )

        print(
            "[RAG]",
            json.dumps(
                rag_test(),
                ensure_ascii=False
            )
        )

        return

    # --------------------------------------------------------
    # A. CONFIGURAZIONE MOUNT IN SOLA LETTURA
    # --------------------------------------------------------

    volumes = [

        {

            "type": "bind",

            "source": str(src),

            "target": "/app/" + name,

            "read_only": True,
        }

        for name, src in sorted(
            missing.items()
        )
    ]

    config = {

        "services": {

            SERVICE: {

                "volumes": volumes
            }
        }
    }

    content = json.dumps(
        config,
        indent=2,
        ensure_ascii=False
    ) + "\n"

    if (
        PATCH.exists()
        and PATCH.read_text(
            encoding="utf-8"
        ) != content
    ):

        raise RuntimeError(
            "Override esistente differente; "
            "nessuna sovrascrittura."
        )

    if not PATCH.exists():

        PATCH.write_text(
            content,
            encoding="utf-8"
        )

        print(
            "[CREATO]",
            PATCH
        )

    # --------------------------------------------------------
    # B. CONTROLLO COMPOSE RISULTANTE
    # --------------------------------------------------------

    svc = inspect_config(True)

    mounted = {
        v.get("target")
        for v in svc.get(
            "volumes", []
        )
    }

    expected = {
        v["target"]
        for v in volumes
    }

    if not expected.issubset(mounted):

        raise RuntimeError(
            "Compose ha scartato "
            "alcuni montaggi"
        )

    if ui_hashes() != before:

        raise RuntimeError(
            "File GUI modificati "
            "prima dell'attivazione"
        )

    # --------------------------------------------------------
    # C. RIAVVIO SOLO GUI V3.2
    # --------------------------------------------------------

    print(
        "[AVVIO] Riavvio della sola "
        "GUI V3.2, senza ricostruire immagini"
    )

    try:

        start_container(True)

        after = container_files(
            missing
        )

        if not all(
            after.get(n)
            for n in missing
        ):

            raise RuntimeError(
                "Uno o piu' script montati "
                "non sono visibili nel container"
            )

        if ui_hashes() != before:

            raise RuntimeError(
                "File GUI originali modificati"
            )

    except Exception:

        print(
            "[ERRORE] Tentativo di rollback"
        )

        try:

            rollback()

        except Exception as exc:

            print(
                "[ERRORE ROLLBACK]",
                exc
            )

        raise

    print(
        "[OK] Script mancanti resi accessibili:",
        len(missing)
    )

    # --------------------------------------------------------
    # D. TEST APPLICATIVO
    # --------------------------------------------------------

    status()

    test = rag_test()

    print(
        "[RAG TEST]",
        json.dumps(
            test,
            ensure_ascii=False,
            indent=2
        )
    )

    # --------------------------------------------------------
    # E. RAPPORTO
    # --------------------------------------------------------

    timestamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP22_6I_GUI_V32_SCRIPTS_"
        + timestamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {

        "step": "22.6I",

        "image": IMAGE,

        "mounted": {
            k: str(v)
            for k, v in missing.items()
        },

        "not_found_on_host": unavailable,

        "rag_http_test": test,

        "gui_source_sha256": before,

        "original_gui_modified": False,

        "training_performed": False,

        "gae_button_not_tested": True,
    }

    path = (
        folder / "RAPPORTO_STEP22_6I.json"
    )

    path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print(
        "[REPORT]",
        path
    )

    if not test.get("ok"):

        print(
            "[DA RISOLVERE] La Chat RAG "
            "richiede ulteriori controlli; "
            "verificare il nuovo errore."
        )


# ============================================================
# 14. AVVIO PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    group = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    for option in (
        "audit",
        "repair",
        "status",
        "rollback",
    ):

        group.add_argument(
            "--" + option,
            action="store_true"
        )

    args = parser.parse_args()

    print("=" * 65)

    print(
        "AI PRIMETTA V3.2 - STEP 22.6I"
    )

    print(
        "RIPRISTINO SCRIPTS GUI / CHAT RAG"
    )

    print("=" * 65)

    if args.audit:

        audit()

    elif args.repair:

        repair()

    elif args.rollback:

        rollback()

    else:

        status()


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