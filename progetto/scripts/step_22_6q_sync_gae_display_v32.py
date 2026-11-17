#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.2 - STEP 22.6Q

--audit
    Controlla la pagina originale senza modificarla.

--install
    Monta una copia della pagina che visualizza
    i risultati dell'ultimo training GAE isolato.

--status
    Mostra la configurazione attiva.

--rollback
    Ripristina la pagina GAE originale.

Non esegue training.
Non modifica i checkpoint scientifici.
Non cambia il layout Streamlit.
"""

import argparse
import ast
import hashlib
import json
import subprocess
import sys
import time
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

PAGE = ROOT / "app/pages/01_Rete_neurale_GAE.py"

PATCH = (
    ROOT
    / "docker/runtime_v32_patch/"
    "01_Rete_neurale_GAE_latest.py"
)

OVERRIDE = (
    ROOT
    / "docker-compose.primetta-v32-gae-display.yaml"
)

REPORTS = ROOT / "target_prediction_v33_vscode/runs"

RUN_ROOT = ROOT / "outputs/gae_gui_safe_runs"

IMAGE = "ai_primetta_v32_full_ai_sklearn18:22.6g"

CONTAINER = "maize_ncrna_gnn"

SERVICE = "maize-ncrna-gnn"

TARGET = "/app/app/pages/01_Rete_neurale_GAE.py"

COMPOSE = [
    ROOT / name
    for name in (
        "docker-compose.yml",
        "docker-compose.primetta-v32-8503.yaml",
        "docker-compose.primetta-v32-full-ai.yaml",
        "docker-compose.primetta-v32-sklearn18.yaml",
        "docker-compose.primetta-v32-scripts.yaml",
        "docker-compose.primetta-v32-rag-persistenza.yaml",
        "docker-compose.primetta-v32-gae-safe.yaml",
    )
]


# ============================================================
# 2. COLLEGAMENTO DINAMICO DEGLI OUTPUT
# ============================================================

# Il codice seguente verrà inserito esclusivamente
# prima della sezione "Stato modulo GAE".
#
# Non modifica il comando di training, i pulsanti,
# la struttura HTML/Streamlit o gli output originali.

DISPLAY_HOOK = '''
# STEP22_6Q_BEGIN: risultati dell'ultimo training isolato

from pathlib import Path as _v32q_Path

_v32q_root = ROOT / "outputs" / "gae_gui_safe_runs"

_v32q_required = (
    "models/gae_model.pt",
    "embeddings/gae_node_embeddings.csv",
    "embeddings/gae_ncrna_embeddings.csv",
    "results/gae_training_history.csv",
    "results/gae_metrics.json",
)

_v32q_finished = sorted(
    (
        p for p in _v32q_root.iterdir()
        if p.is_dir()
        and not p.is_symlink()
        and all(
            (p / n).is_file()
            and (p / n).stat().st_size > 0
            for n in _v32q_required
        )
    ),
    key=lambda p: p.name,
) if _v32q_root.is_dir() else []

_v32q_run = (
    _v32q_finished[-1]
    if _v32q_finished
    else None
)


def _v32q_relocate(value):
    """
    Reindirizza i percorsi legacy GAE
    alla cartella dell'ultimo training completo.
    """

    if isinstance(value, _v32q_Path):

        path = value
        as_string = False

    elif (
        isinstance(value, str)
        and value.startswith("/")
    ):

        path = _v32q_Path(value)
        as_string = True

    elif isinstance(value, list):

        return [
            _v32q_relocate(v)
            for v in value
        ]

    elif isinstance(value, tuple):

        return tuple(
            _v32q_relocate(v)
            for v in value
        )

    elif isinstance(value, dict):

        return {
            k: _v32q_relocate(v)
            for k, v in value.items()
        }

    else:

        return value

    parts = path.parts

    alternatives = (
        (
            ("outputs", "gae"),
            "results"
        ),
        (
            ("models", "graph_autoencoder"),
            "models"
        ),
        (
            ("data", "processed", "embeddings"),
            "embeddings"
        ),
    )

    for prefix, directory in alternatives:

        for start in range(
            0,
            len(parts) - len(prefix) + 1
        ):

            if (
                tuple(
                    parts[start:start + len(prefix)]
                )
                != prefix
            ):
                continue

            suffix = parts[
                start + len(prefix):
            ]

            new_path = (
                _v32q_run / directory
            )

            for component in suffix:

                new_path = (
                    new_path / component
                )

            return (
                str(new_path)
                if as_string
                else new_path
            )

    return value


if _v32q_run is not None:

    _v32q_keywords = (
        "DIR",
        "PATH",
        "FILE",
        "PLOT",
        "SUMMARY",
        "METRIC",
        "MODEL",
        "EMB",
        "HISTORY",
        "CSV",
        "JSON",
        "TXT",
        "OUTPUT",
    )

    for _v32q_key, _v32q_value in list(
        globals().items()
    ):

        if (
            _v32q_key.startswith("_")
            or _v32q_key == "ROOT"
            or not any(
                k in _v32q_key.upper()
                for k in _v32q_keywords
            )
        ):
            continue

        _v32q_new = _v32q_relocate(
            _v32q_value
        )

        if (
            _v32q_new is not _v32q_value
            and _v32q_new != _v32q_value
        ):

            globals()[
                _v32q_key
            ] = _v32q_new

# STEP22_6Q_END
'''


# ============================================================
# 3. COMANDI DI SISTEMA
# ============================================================

def execute(args, timeout=150, check=True):

    proc = subprocess.run(
        list(map(str, args)),
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )

    if check and proc.returncode:

        raise RuntimeError(
            "Comando fallito: "
            + " ".join(map(str, args[:8]))
            + "\n"
            + proc.stderr[-1000:]
        )

    return proc


# ============================================================
# 4. CHECKSUM FILE
# ============================================================

def sha(path):

    h = hashlib.sha256()

    with path.open("rb") as stream:

        for part in iter(
            lambda: stream.read(
                1024 * 1024
            ),
            b"",
        ):

            h.update(part)

    return h.hexdigest()


# ============================================================
# 5. CONTROLLO PULSANTI STREAMLIT
# ============================================================

def page_ast(text):

    tree = ast.parse(text)

    buttons = sorted(
        ast.unparse(node)

        for node in ast.walk(tree)

        if isinstance(node, ast.Call)

        and ast.unparse(node.func) in (
            "st.button",
            "st.form_submit_button",
        )
    )

    return buttons


# ============================================================
# 6. PREPARAZIONE PAGINA CORRETTA
# ============================================================

def patch_text():

    source = PAGE.read_text(
        encoding="utf-8"
    )

    if "STEP22_6Q_BEGIN" in source:

        raise RuntimeError(
            "La pagina originale contiene "
            "gia' la modifica 22.6Q"
        )

    lines = source.splitlines(
        keepends=True
    )

    anchors = [
        (i, line)

        for i, line in enumerate(lines)

        if "Stato modulo GAE" in line
        and "st." in line
    ]

    if len(anchors) != 1:

        raise RuntimeError(
            "Intestazione 'Stato modulo GAE' "
            "non univoca: nessuna modifica"
        )

    index, line = anchors[0]

    if line[:len(line) - len(line.lstrip())]:

        raise RuntimeError(
            "Sezione stato annidata: "
            "necessario controllo del sorgente"
        )

    tail = "".join(
        lines[index:]
    )

    legacy_literals = (
        '"outputs/gae',
        "'outputs/gae",
        '"models/graph_autoencoder',
        "'models/graph_autoencoder",
        '"data/processed/embeddings',
        "'data/processed/embeddings",
    )

    if any(
        term in tail
        for term in legacy_literals
    ):

        raise RuntimeError(
            "Percorsi legacy definiti direttamente "
            "nella sezione Stato GAE: "
            "modifica automatica sospesa"
        )

    result = (
        "".join(lines[:index])
        + DISPLAY_HOOK
        + "\n"
        + "".join(lines[index:])
    )

    compile(
        result,
        str(PATCH),
        "exec"
    )

    if page_ast(result) != page_ast(source):

        raise RuntimeError(
            "Pulsanti Streamlit modificati: "
            "patch rifiutata"
        )

    return (
        result,
        index + 1,
        len(page_ast(source)),
    )


# ============================================================
# 7. CONFIGURAZIONE COMPOSE
# ============================================================

def compose(patched=False):

    files = COMPOSE + (
        [OVERRIDE]
        if patched
        else []
    )

    for path in files:

        if not path.is_file():

            raise FileNotFoundError(
                path
            )

    command = [
        "docker",
        "compose",
    ]

    for path in files:

        command.extend([
            "-f",
            str(path),
        ])

    return command


def configuration(patched=False):

    response = execute(
        compose(patched)
        + [
            "config",
            "--format",
            "json",
        ]
    )

    service = json.loads(
        response.stdout
    )["services"][SERVICE]

    if service.get("image") != IMAGE:

        raise RuntimeError(
            "Immagine Docker differente "
            "dalla versione collaudata"
        )

    ports = service.get(
        "ports", []
    )

    valid_port = any(
        str(port.get("published")) == "8503"
        and str(port.get("target")) == "8501"
        and port.get("host_ip") == "127.0.0.1"
        for port in ports
    )

    if not valid_port:

        raise RuntimeError(
            "Porta GUI 8503 non conforme"
        )

    mounts = {
        m["target"]: m
        for m in service.get(
            "volumes", []
        )
    }

    if patched:

        item = mounts.get(
            TARGET
        )

        if (
            not item
            or Path(
                item["source"]
            ).resolve() != PATCH.resolve()
            or not item.get(
                "read_only"
            )
        ):

            raise RuntimeError(
                "Montaggio pagina GAE non sicuro"
            )

    safe_trainer = mounts.get(
        "/app/25_train_graph_autoencoder.py"
    )

    if (
        not safe_trainer
        or not safe_trainer.get(
            "read_only"
        )
        or "25_train_graph_autoencoder_safe.py"
        not in safe_trainer["source"]
    ):

        raise RuntimeError(
            "Trainer GAE isolato "
            "non configurato"
        )

    return mounts


# ============================================================
# 8. STATO CONTAINER
# ============================================================

def container_info():

    proc = execute(
        [
            "docker",
            "inspect",
            CONTAINER,
        ],
        check=False,
    )

    if proc.returncode:

        return {}

    data = json.loads(
        proc.stdout
    )[0]

    mounts = {
        m["Destination"]: m
        for m in data.get(
            "Mounts", []
        )
    }

    return {
        "running": data[
            "State"
        ]["Running"],

        "image": data[
            "Config"
        ]["Image"],

        "page_source": mounts.get(
            TARGET, {}
        ).get("Source"),

        "trainer_source": mounts.get(
            "/app/25_train_graph_autoencoder.py",
            {}
        ).get("Source"),
    }


# ============================================================
# 9. HEALTHCHECK
# ============================================================

def health(port, path):

    try:

        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            f"http://127.0.0.1:{port}{path}",
            timeout=5,
        ) as response:

            return response.status == 200

    except Exception:

        return False


# ============================================================
# 10. IDENTIFICAZIONE ULTIMO TRAINING
# ============================================================

def last_run():

    if not RUN_ROOT.is_dir():

        return None

    required = [
        "models/gae_model.pt",
        "embeddings/gae_node_embeddings.csv",
        "embeddings/gae_ncrna_embeddings.csv",
        "results/gae_training_history.csv",
        "results/gae_metrics.json",
    ]

    finished = sorted(
        (
            p
            for p in RUN_ROOT.iterdir()

            if p.is_dir()
            and not p.is_symlink()
            and all(
                (p / name).is_file()
                and (p / name).stat().st_size > 0
                for name in required
            )
        ),
        key=lambda p: p.name,
    )

    return (
        finished[-1]
        if finished
        else None
    )


# ============================================================
# 11. CREAZIONE NON DISTRUTTIVA DEI FILE
# ============================================================

def make_file(path, text):

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if path.exists():

        if path.read_text(
            encoding="utf-8"
        ) != text:

            raise RuntimeError(
                "File preesistente differente: "
                + str(path)
            )

    else:

        path.write_text(
            text,
            encoding="utf-8"
        )

        print(
            "[CREATO]",
            path.relative_to(ROOT)
        )


# ============================================================
# 12. RIAVVIO ESCLUSIVO DELLA GUI V3.2
# ============================================================

def start(patched):

    configuration(patched)

    execute(
        compose(patched)
        + [
            "up",
            "-d",
            "--no-deps",
            "--no-build",
            "--pull",
            "never",
            "--force-recreate",
            SERVICE,
        ],
        timeout=170,
    )

    for _ in range(30):

        info = container_info()

        if (
            info.get("running")
            and info.get("image") == IMAGE
            and health(
                8503,
                "/_stcore/health",
            )
            and health(
                8000,
                "/health",
            )
        ):

            return

        time.sleep(2)

    raise RuntimeError(
        "GUI/API V3.2 non raggiungibili"
    )


# ============================================================
# 13. STATO GENERALE
# ============================================================

def status():

    state = container_info()

    print(
        "[IMMAGINE]",
        state.get("image"),
    )

    print(
        "[PAGINA GAE ATTIVA]",
        state.get("page_source")
        or "originale immagine",
    )

    print(
        "[TRAINER SICURO]",
        state.get("trainer_source"),
    )

    latest = last_run()

    print(
        "[ULTIMO TRAINING COMPLETO]",
        (
            latest.relative_to(ROOT)
            if latest
            else "nessuno"
        ),
    )

    if latest:

        metrics = json.loads(
            (
                latest
                / "results/gae_metrics.json"
            ).read_text(
                encoding="utf-8"
            )
        )

        print(
            "[METRICHE NUOVO RUN - CHIAVI]",
            sorted(metrics.keys())[:25],
        )

    services = (
        ("GUI V3.2", 8503, "/_stcore/health"),
        ("API V3.2", 8000, "/health"),
        ("GUI V3.3", 8502, "/_stcore/health"),
        ("AI ORZO", 8501, "/_stcore/health"),
    )

    for name, port, route in services:

        print(
            "[OK]"
            if health(port, route)
            else "[NON VERIFICATO]",
            name,
        )


# ============================================================
# 14. AUDIT PREVENTIVO
# ============================================================

def audit():

    configuration(False)

    patched, line, buttons = patch_text()

    latest = last_run()

    print(
        "[OK] Pagina originale leggibile"
    )

    print(
        "[OK] Pulsanti preservati:",
        buttons,
    )

    print(
        "[INSERIMENTO PRIMA DELLA RIGA]",
        line,
    )

    print(
        "[ULTIMO RUN COMPLETO]",
        (
            latest.relative_to(ROOT)
            if latest
            else "nessuno"
        ),
    )

    print(
        "[OK] Anteprima compilabile:",
        len(patched.splitlines()),
        "righe",
    )

    print(
        "[OK] Solo audit: "
        "nessun file modificato"
    )


# ============================================================
# 15. INSTALLAZIONE COLLEGAMENTO GAE
# ============================================================

def install():

    configuration(False)

    info = container_info()

    if (
        not info.get("running")
        or info.get("image") != IMAGE
    ):

        raise RuntimeError(
            "Container V3.2 non conforme"
        )

    if last_run() is None:

        raise RuntimeError(
            "Nessun training GAE isolato completo"
        )

    patched, line, buttons = patch_text()

    original = sha(PAGE)

    trainer = (
        ROOT / "25_train_graph_autoencoder.py"
    )

    trainer_hash = sha(trainer)

    # --------------------------------------------------------
    # A. COPIA DELLA PAGINA
    # --------------------------------------------------------

    make_file(
        PATCH,
        patched,
    )

    # --------------------------------------------------------
    # B. OVERRIDE DOCKER
    # --------------------------------------------------------

    override = json.dumps(
        {
            "services": {
                SERVICE: {
                    "volumes": [
                        {
                            "type": "bind",
                            "source": str(PATCH),
                            "target": TARGET,
                            "read_only": True,
                        }
                    ]
                }
            }
        },
        indent=2,
    ) + "\n"

    make_file(
        OVERRIDE,
        override,
    )

    configuration(True)

    if (
        sha(PAGE) != original
        or sha(trainer) != trainer_hash
    ):

        raise RuntimeError(
            "Sorgenti originali modificati"
        )

    # --------------------------------------------------------
    # C. ATTIVAZIONE
    # --------------------------------------------------------

    try:

        start(True)

        live = container_info()

        if (
            Path(
                live.get("page_source") or ""
            ).resolve()
            != PATCH.resolve()
        ):

            raise RuntimeError(
                "Pagina GAE montata non conforme"
            )

        code = (
            "import pathlib,hashlib;"
            "print(hashlib.sha256("
            "pathlib.Path("
            "'/app/app/pages/01_Rete_neurale_GAE.py'"
            ").read_bytes()).hexdigest())"
        )

        verification = execute(
            [
                "docker",
                "exec",
                CONTAINER,
                "python3",
                "-B",
                "-c",
                code,
            ]
        )

        if (
            verification.stdout.strip()
            != sha(PATCH)
        ):

            raise RuntimeError(
                "Checksum della pagina Docker "
                "non corrispondente"
            )

        if (
            sha(PAGE) != original
            or sha(trainer) != trainer_hash
        ):

            raise RuntimeError(
                "Sorgenti originali cambiati"
            )

    except Exception:

        print(
            "[RIPRISTINO] Tentativo di "
            "ritorno alla pagina originale"
        )

        try:

            start(False)

        except Exception as exc:

            print(
                "[ATTENZIONE] Rollback:",
                exc,
            )

        raise

    # --------------------------------------------------------
    # D. REPORT
    # --------------------------------------------------------

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = REPORTS / (
        "STEP22_6Q_GAE_DISPLAY_" + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False,
    )

    report_data = {
        "step": "22.6Q",
        "status": "PATCH_MONTATA",
        "original_page_sha256": original,
        "original_trainer_sha256": trainer_hash,
        "patch_sha256": sha(PATCH),
        "buttons_unchanged": buttons,
        "insert_line": line,
        "display_run": str(last_run()),
        "training_performed": False,
        "original_gui_source_modified": False,
        "visual_verification_pending": True,
    }

    report = (
        folder / "RAPPORTO_STEP22_6Q.json"
    )

    report.write_text(
        json.dumps(
            report_data,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print(
        "[OK] Pagina GAE collegata "
        "all'ultimo training completo"
    )

    print(
        "[OK] Nessun training eseguito"
    )

    print(
        "[OK] File originali invariati"
    )

    print(
        "[REPORT]",
        report,
    )

    status()


# ============================================================
# 16. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    options = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    for mode in (
        "audit",
        "install",
        "status",
        "rollback",
    ):

        options.add_argument(
            "--" + mode,
            action="store_true",
        )

    args = parser.parse_args()

    print("=" * 68)

    print(
        "AI PRIMETTA V3.2 - STEP 22.6Q"
    )

    print(
        "COLLEGAMENTO RISULTATI GAE ALLA GUI"
    )

    print("=" * 68)

    if args.audit:

        audit()

    elif args.install:

        install()

    elif args.status:

        status()

    else:

        start(False)

        print(
            "[OK] Pagina originale ripristinata"
        )

        status()


if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__,
            exc,
            file=sys.stderr,
        )

        sys.exit(1)