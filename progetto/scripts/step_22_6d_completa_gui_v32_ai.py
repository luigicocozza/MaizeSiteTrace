#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA - STEP 22.6D

POTENZIAMENTO DEL RUNTIME GUI ORIGINALE V3.2

Funzioni:

--install
    Costruisce una nuova immagine Docker con
    PyTorch, PyTorch Geometric e NetworkX.
    Verifica il funzionamento del runtime AI.
    Attiva la nuova immagine sulla porta 8503.

--status
    Controlla lo stato delle GUI e dell'API.

--rollback
    Ripristina l'immagine originale V3.2,
    mantenendo la porta 8503.

Non modifica GUI, dati, checkpoint,
AI ORZO o monitor ML V3.3.

Non avvia nuovi training.
"""

import argparse
import hashlib
import json
import re
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

BASE_COMPOSE = ROOT / "docker-compose.yml"

PORT_COMPOSE = (
    ROOT / "docker-compose.primetta-v32-8503.yaml"
)

FULL_COMPOSE = (
    ROOT / "docker-compose.primetta-v32-full-ai.yaml"
)

DOCKERFILE = (
    ROOT / "docker" / "Dockerfile.v32_full_ai_cpu"
)

SERVICE = "maize-ncrna-gnn"

NEW_IMAGE = "ai_primetta_v32_full_ai_cpu:22.6d"

RUNS = (
    ROOT / "target_prediction_v33_vscode" / "runs"
)

PAGES = [
    ROOT / "app/streamlit_app.py"
] + sorted(
    (ROOT / "app/pages").glob("*.py")
)

PROTECTED = PAGES + [

    ROOT
    / "models/graph_autoencoder_v3/"
    "gae_v3_model.pt",

    ROOT
    / "data/processed/pyg/"
    "maize_ncrna_graph_v3_prjna419326.pt",

    ROOT
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib",
]


# ============================================================
# 2. COLLAUDO GNN NEL CONTAINER
# ============================================================

SMOKE = r'''

import json
import os

import torch
import torch_geometric
import networkx
import streamlit
import fastapi
import sklearn
import joblib

from torch_geometric.nn import GCNConv
from torch_geometric.data import Data

from pathlib import Path


# ------------------------------------------------------------
# TEST DI UNA GCN REALE
# ------------------------------------------------------------

x = torch.tensor([
    [1., 0.],
    [0., 1.],
    [1., 1.]
])

edge_index = torch.tensor(
    [
        [0, 1, 2, 2],
        [1, 2, 0, 1]
    ],
    dtype=torch.long
)

model = GCNConv(
    2,
    4
)

output = model(
    x,
    edge_index
)

assert tuple(output.shape) == (3, 4)

loss = output.square().mean()

loss.backward()

assert model.lin.weight.grad is not None


# ------------------------------------------------------------
# CONTROLLO RUNTIME CPU
# ------------------------------------------------------------

assert (
    not torch.cuda.is_available()
    or os.environ.get("CUDA_VISIBLE_DEVICES") == "-1"
)


# ------------------------------------------------------------
# CONTROLLO GRAFO E CHECKPOINT DEL PROGETTO
# ------------------------------------------------------------

graph_path = Path(
    "/project/data/processed/pyg/"
    "maize_ncrna_graph_v3_prjna419326.pt"
)

checkpoint_path = Path(
    "/project/models/graph_autoencoder_v3/"
    "gae_v3_model.pt"
)

assert graph_path.is_file()

assert checkpoint_path.is_file()


# ATTENZIONE:
# weights_only=False va usato solo per checkpoint
# attendibili del repository scientifico.

graph = torch.load(
    graph_path,
    map_location="cpu",
    weights_only=False
)

checkpoint = torch.load(
    checkpoint_path,
    map_location="cpu",
    weights_only=False
)


assert hasattr(
    graph,
    "edge_index"
)

assert graph.edge_index.shape[0] == 2

assert isinstance(
    checkpoint,
    dict
)


# ------------------------------------------------------------
# RISULTATO COLLAUDO
# ------------------------------------------------------------

result = {

    "torch": torch.__version__,

    "torch_geometric": torch_geometric.__version__,

    "networkx": networkx.__version__,

    "sklearn": sklearn.__version__,

    "streamlit": streamlit.__version__,

    "graph_nodes": int(
        graph.num_nodes
    ),

    "graph_edges": int(
        graph.edge_index.shape[1]
    ),

    "checkpoint_loaded": True,

    "gcn_forward_backward": True,

    "training_performed": False
}

print(
    "FULL_AI_SMOKE="
    + json.dumps(result)
)

'''


# ============================================================
# 3. ESECUZIONE COMANDI
# ============================================================

def execute(
    args,
    timeout=120,
    capture=False,
    check=True
):

    print(
        "$ " + " ".join(map(str, args)),
        flush=True
    )

    result = subprocess.run(
        args,
        text=True,
        capture_output=capture,
        timeout=timeout,
        check=False
    )

    if check and result.returncode:

        raise RuntimeError(
            "Comando fallito, codice "
            + str(result.returncode)
            + ": "
            + (result.stderr or "")[-1400:]
        )

    return result


# ============================================================
# 4. COMPOSE
# ============================================================

def docker_compose(full=False):

    files = [
        BASE_COMPOSE,
        PORT_COMPOSE
    ]

    if full:

        files.append(
            FULL_COMPOSE
        )

    command = [
        "docker",
        "compose"
    ]

    for path in files:

        command.extend([
            "-f",
            str(path)
        ])

    return command


def read_config(full=False):

    response = execute(
        [
            *docker_compose(full),
            "config",
            "--format",
            "json"
        ],
        capture=True
    )

    return json.loads(
        response.stdout
    )


# ============================================================
# 5. CONTROLLO PORTE
# ============================================================

def ports_ok(config):

    service = config[
        "services"
    ][SERVICE]

    ports = service.get(
        "ports",
        []
    )

    good = [
        p
        for p in ports

        if (
            str(p.get("published")) == "8503"
            and str(p.get("target")) == "8501"
        )
    ]

    conflicting = [
        p
        for p in ports

        if str(
            p.get("published")
        ) == "8501"
    ]

    if len(good) != 1 or conflicting:

        raise RuntimeError(
            "Configurazione porte non sicura: "
            "GUI V3.2 deve usare 8503 -> 8501. "
            "La porta host 8501 appartiene ad AI ORZO."
        )

    return service


# ============================================================
# 6. IMMAGINE ORIGINALE
# ============================================================

def base_image():

    service = ports_ok(
        read_config(False)
    )

    image = service.get(
        "image"
    )

    if not image:

        response = execute(
            [
                "docker",
                "inspect",
                "maize_ncrna_gnn",
                "--format",
                "{{.Config.Image}}"
            ],
            capture=True
        )

        image = response.stdout.strip()

    if (
        not image
        or not re.fullmatch(
            r"[A-Za-z0-9_./:@-]+",
            image
        )
    ):

        raise RuntimeError(
            "Immagine originale non valida"
        )

    execute(
        [
            "docker",
            "image",
            "inspect",
            image
        ],
        capture=True
    )

    return image


# ============================================================
# 7. VERSIONE PYTHON
# ============================================================

def python_version(image):

    response = execute(
        [
            "docker",
            "run",
            "--rm",
            "--pull",
            "never",
            "--network",
            "none",
            "--entrypoint",
            "python3",
            image,
            "-c",
            (
                "import sys,json;"
                "print(json.dumps({"
                "'version':list(sys.version_info[:3])"
                "}))"
            )
        ],
        capture=True
    )

    version = json.loads(
        response.stdout
    )["version"]

    if not (
        (3, 10)
        <= tuple(version[:2])
        <= (3, 14)
    ):

        raise RuntimeError(
            "Python non supportato: "
            + str(version)
        )

    print(
        "[OK] Python:",
        ".".join(map(str, version))
    )

    return version


def original_user(image):

    response = execute(
        [
            "docker",
            "image",
            "inspect",
            image,
            "--format",
            "{{json .Config.User}}"
        ],
        capture=True
    )

    return (
        json.loads(
            response.stdout.strip()
        )
        or "root"
    )


# ============================================================
# 8. INTEGRITÀ FILE ORIGINALI
# ============================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):

            h.update(chunk)

    return h.hexdigest()


def originals():

    missing = [
        str(path)
        for path in PROTECTED
        if not path.is_file()
    ]

    if missing:

        raise RuntimeError(
            "File mancanti: "
            + "; ".join(missing)
        )

    return {
        str(path.relative_to(ROOT)): sha256(path)
        for path in PROTECTED
    }


# ============================================================
# 9. SCRITTURA FILE SENZA SOVRASCRITTURA
# ============================================================

def write_if_new(path, content):

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if path.exists():

        existing = path.read_text(
            encoding="utf-8"
        )

        if existing != content:

            raise RuntimeError(
                "File già presente con contenuto diverso: "
                + str(path)
            )

        print(
            "[OK] Già presente:",
            path.relative_to(ROOT)
        )

    else:

        path.write_text(
            content,
            encoding="utf-8"
        )

        print(
            "[OK] Generato:",
            path.relative_to(ROOT)
        )


# ============================================================
# 10. GENERAZIONE DOCKERFILE E OVERRIDE
# ============================================================

def generate(image, python, user):

    # NetworkX 3.4.2 per Python 3.10;
    # NetworkX 3.6.1 per Python >=3.11.

    nx_version = (
        "3.4.2"
        if tuple(python[:2]) == (3, 10)
        else "3.6.1"
    )

    dockerfile = (

        "# AI PRIMETTA - Runtime CPU ML/GNN\n"

        f"FROM {image}\n"

        "USER root\n"

        "ENV PIP_DISABLE_PIP_VERSION_CHECK=1 "
        "PIP_NO_CACHE_DIR=1 "
        "CUDA_VISIBLE_DEVICES=-1\n"

        "RUN python3 -m pip install "
        "--no-cache-dir "
        "--index-url https://download.pytorch.org/whl/cpu "
        "'torch==2.11.0'\n"

        "RUN python3 -m pip install "
        "--no-cache-dir "
        "'torch-geometric==2.8.1' "
        f"'networkx=={nx_version}'\n"

        f"USER {user}\n"
    )

    override = (

        "# AI PRIMETTA V3.2 - Immagine completa\n"

        "services:\n"

        f"  {SERVICE}:\n"

        f"    image: {NEW_IMAGE}\n"

        "    pull_policy: never\n"
    )

    write_if_new(
        DOCKERFILE,
        dockerfile
    )

    write_if_new(
        FULL_COMPOSE,
        override
    )

    config = read_config(True)

    service = ports_ok(config)

    if service.get("image") != NEW_IMAGE:

        raise RuntimeError(
            "Immagine nuova non applicata nel Compose"
        )

    print(
        "[OK] Porta GUI originale mantenuta: 8503"
    )

    print(
        "[OK] Nessuna modifica alla porta 8501 di AI ORZO"
    )


# ============================================================
# 11. COLLAUDO IMMAGINE AI
# ============================================================

def smoke_image():

    command = [

        "docker",
        "run",
        "--rm",
        "--pull",
        "never",

        "--network",
        "none",

        "--read-only",

        "--tmpfs",
        "/tmp:rw,size=256m",

        "--cpus",
        "2",

        "-e",
        "OMP_NUM_THREADS=2",

        "-e",
        "CUDA_VISIBLE_DEVICES=-1",

        "-v",
        f"{ROOT}:/project:ro",

        "--entrypoint",
        "python3",

        NEW_IMAGE,

        "-B",
        "-c",
        SMOKE
    ]

    response = execute(
        command,
        timeout=240,
        capture=True
    )

    print(
        response.stdout
    )

    marker = "FULL_AI_SMOKE="

    if marker not in response.stdout:

        raise RuntimeError(
            "Smoke test GNN non riuscito"
        )

    payload = response.stdout.split(
        marker
    )[-1].splitlines()[0]

    return json.loads(payload)


# ============================================================
# 12. CONTROLLO HTTP
# ============================================================

def http_ok(port, path):

    url = (
        f"http://127.0.0.1:{port}{path}"
    )

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
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
# 13. IMMAGINE ATTIVA
# ============================================================

def active_image():

    response = execute(
        [
            "docker",
            "inspect",
            "maize_ncrna_gnn",
            "--format",
            "{{.Config.Image}}"
        ],
        capture=True,
        check=False
    )

    if response.returncode == 0:

        return response.stdout.strip()

    return "NON_TROVATO"


# ============================================================
# 14. STATO APPLICAZIONI
# ============================================================

def current_status():

    command = docker_compose(
        FULL_COMPOSE.exists()
    )

    execute(
        [
            *command,
            "ps",
            SERVICE
        ],
        check=False
    )

    print(
        "\nImmagine attiva:",
        active_image()
    )

    services = [
        (
            "GUI PRIMETTA V3.2",
            8503,
            "/_stcore/health"
        ),
        (
            "GUI PRIMETTA V3.3",
            8502,
            "/_stcore/health"
        ),
        (
            "GUI AI ORZO",
            8501,
            "/_stcore/health"
        ),
        (
            "API PRIMETTA V3.2",
            8000,
            "/health"
        )
    ]

    for name, port, path in services:

        success = http_ok(
            port,
            path
        )

        print(
            "[OK]" if success else "[DA VERIFICARE]",
            name,
            "porta",
            port
        )


# ============================================================
# 15. RIPRISTINO CONFIGURAZIONE PRECEDENTE
# ============================================================

def rollback():

    print(
        "[ROLLBACK] Ripristino immagine originale V3.2"
    )

    execute(
        [
            *docker_compose(False),

            "up",
            "-d",

            "--no-deps",
            "--no-build",
            "--pull",
            "never",

            SERVICE
        ],
        timeout=180
    )

    current_status()


# ============================================================
# 16. INSTALLAZIONE RUNTIME COMPLETO
# ============================================================

def install():

    for path in (
        BASE_COMPOSE,
        PORT_COMPOSE
    ):

        if not path.is_file():

            raise RuntimeError(
                "Compose mancante: "
                + str(path)
            )

    # --------------------------------------------------------
    # STATO INIZIALE
    # --------------------------------------------------------

    print(
        "\n[1] STATO PRIMA DELL'INTERVENTO"
    )

    current_status()

    before = originals()

    image = base_image()

    if image == NEW_IMAGE:

        raise RuntimeError(
            "Il Compose originale punta "
            "già alla nuova immagine"
        )

    python = python_version(image)

    user = original_user(image)

    # --------------------------------------------------------
    # GENERAZIONE FILE
    # --------------------------------------------------------

    print(
        "\n[2] GENERAZIONE FILE DOCKER"
    )

    generate(
        image,
        python,
        user
    )

    # --------------------------------------------------------
    # NUOVA IMMAGINE
    # --------------------------------------------------------

    print(
        "\n[3] COSTRUZIONE IMMAGINE AI V3.2"
    )

    execute(
        [
            "docker",
            "build",

            "-f",
            str(DOCKERFILE),

            "-t",
            NEW_IMAGE,

            str(DOCKERFILE.parent)
        ],
        timeout=1800
    )

    # --------------------------------------------------------
    # COLLAUDO GCN
    # --------------------------------------------------------

    print(
        "\n[4] COLLAUDO GNN IN CONTAINER ISOLATO"
    )

    proof = smoke_image()

    if before != originals():

        raise RuntimeError(
            "File originali variati prima dell'attivazione"
        )

    # --------------------------------------------------------
    # ATTIVAZIONE
    # --------------------------------------------------------

    print(
        "\n[5] ATTIVAZIONE NUOVA GUI V3.2"
    )

    try:

        execute(
            [
                *docker_compose(True),

                "up",
                "-d",

                "--no-deps",
                "--no-build",
                "--pull",
                "never",

                SERVICE
            ],
            timeout=180
        )

        healthy = False

        for _ in range(18):

            if (
                http_ok(
                    8503,
                    "/_stcore/health"
                )
                and http_ok(
                    8000,
                    "/health"
                )
            ):

                healthy = True
                break

            time.sleep(1)

        if (
            active_image() != NEW_IMAGE
            or not healthy
        ):

            raise RuntimeError(
                "Nuova GUI/API non operativa"
            )

    except Exception:

        print(
            "[ATTENZIONE] Avvio non valido. "
            "Ripristino automatico immagine originale."
        )

        rollback()

        raise

    # --------------------------------------------------------
    # VERIFICA INTEGRITÀ
    # --------------------------------------------------------

    if before != originals():

        raise RuntimeError(
            "File originali modificati: "
            "necessaria verifica"
        )

    # --------------------------------------------------------
    # RIEPILOGO
    # --------------------------------------------------------

    print(
        "\n[6] STATO COMPLESSIVO"
    )

    current_status()

    # --------------------------------------------------------
    # RAPPORTO DI COLLAUDO
    # --------------------------------------------------------

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP22_6D_GUI_V32_FULL_AI_"
        + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {

        "project": "AI PRIMETTA",

        "step": "22.6D",

        "status": "IMAGE_UPGRADE_AND_SMOKE_PASS",

        "old_image": image,

        "new_image": NEW_IMAGE,

        "test_ml_gnn": proof,

        "gui_v32_http_200": http_ok(
            8503,
            "/_stcore/health"
        ),

        "api_http_200": http_ok(
            8000,
            "/health"
        ),

        "original_sha256": before,

        "training_performed": False,

        "complete_gui_functional_validation": False,

        "note": (
            "Dipendenze installate e test GCN superato. "
            "I pulsanti GAE, upload e Chat RAG "
            "richiedono verifiche funzionali distinte."
        )
    }

    report_path = (
        folder / "REPORT_STEP22_6D.json"
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
    print("[OK] STEP 22.6D INSTALLAZIONE COMPLETATA")
    print("=" * 70)

    print(
        "Nuova immagine:",
        NEW_IMAGE
    )

    print(
        "GUI V3.2:",
        "http://127.0.0.1:8503"
    )

    print(
        "Report:",
        report_path
    )

    print(
        "[OK] Checkpoint originali preservati"
    )

    print(
        "[OK] Nessun training"
    )

    print(
        "[PROSSIMO] Collaudo funzionale "
        "GAE, caricamento dati, ranking e Chat RAG"
    )


# ============================================================
# 17. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    group = parser.add_mutually_exclusive_group(
        required=True
    )

    group.add_argument(
        "--install",
        action="store_true"
    )

    group.add_argument(
        "--status",
        action="store_true"
    )

    group.add_argument(
        "--rollback",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 70)
    print("AI PRIMETTA V3.2 — STEP 22.6D")
    print("POTENZIAMENTO RUNTIME ML/GNN")
    print("=" * 70)

    if args.install:

        install()

    elif args.rollback:

        rollback()

    else:

        current_status()


# ============================================================
# 18. AVVIO
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