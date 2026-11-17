#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.2 STEP 22.6G.

Corregge sklearn senza modificare le GUI.

Esempi:
    python3 scripts/step_22_6g_fix_v32.py --install
    python3 scripts/step_22_6g_fix_v32.py --status
    python3 scripts/step_22_6g_fix_v32.py --rollback
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

OLD = "ai_primetta_v32_full_ai_cpu:22.6d"

NEW = "ai_primetta_v32_full_ai_sklearn18:22.6g"

SERVICE = "maize-ncrna-gnn"

CONTAINER = "maize_ncrna_gnn"

C1 = ROOT / "docker-compose.yml"

C2 = ROOT / "docker-compose.primetta-v32-8503.yaml"

C3 = ROOT / "docker-compose.primetta-v32-full-ai.yaml"

C4 = ROOT / "docker-compose.primetta-v32-sklearn18.yaml"

DOCKERFILE = ROOT / "docker/Dockerfile.v32_sklearn18"

RUNS = ROOT / "target_prediction_v33_vscode/runs"


# ============================================================
# 2. TEST REALE GAE E ML NEL CONTAINER
# ============================================================

WORKER = r'''
import json
import re
import warnings

from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
import joblib
import torch
import torch_geometric
import networkx

from torch import nn
from torch_geometric.nn import GCNConv

from sklearn.exceptions import InconsistentVersionWarning


p = Path("/project")

out = {
    "passed": False,
    "tests": {},
    "details": {},
    "errors": []
}


def check(name, fn):

    try:

        out["details"][name] = fn()

        out["tests"][name] = True

    except Exception as exc:

        out["tests"][name] = False

        out["errors"].append(
            name + ": " + repr(exc)
        )


# ============================================================
# A. VERSIONI LIBRERIE
# ============================================================

def versions():

    assert sklearn.__version__ == "1.8.0"

    assert torch.__version__.startswith("2.11.0")

    assert torch_geometric.__version__ == "2.8.1"

    return {
        "sklearn": sklearn.__version__,
        "torch": torch.__version__,
        "pyg": torch_geometric.__version__,
        "networkx": networkx.__version__
    }


check("versions", versions)


# ============================================================
# B. INFERENZA REALE GRAPH AUTOENCODER V3
# ============================================================

def gae():

    graph_path = (
        p
        / "data/processed/pyg/"
        "maize_ncrna_graph_v3_prjna419326.pt"
    )

    ck_path = (
        p
        / "models/graph_autoencoder_v3/"
        "gae_v3_model.pt"
    )

    # Solo artefatti scientifici locali attendibili.
    payload = torch.load(
        graph_path,
        map_location="cpu",
        weights_only=False
    )

    ck = torch.load(
        ck_path,
        map_location="cpu",
        weights_only=False
    )

    graph = (
        payload["data"]
        if isinstance(payload, dict)
        else payload
    )

    assert graph.x.ndim == 2

    assert graph.edge_index.shape[0] == 2

    x = graph.x.float()

    edges = graph.edge_index.long()

    n_in = int(ck["in_channels"])

    n_hidden = int(ck["hidden_channels"])

    n_out = int(ck["latent_channels"])

    assert x.shape[1] == n_in

    assert n_in == 47

    assert int(graph.num_nodes) == 884

    assert int(edges.shape[1]) == 4782

    state = ck["model_state_dict"]

    assert isinstance(state, dict)

    keys = list(state)

    # Gestione dei prefissi nei checkpoint.
    for prefix in (
        "module.encoder.",
        "model.encoder.",
        "encoder.",
        "module."
    ):

        if all(
            key.startswith(prefix)
            for key in keys
        ):

            state = {
                key[len(prefix):]: value
                for key, value in state.items()
            }

            break

    # --------------------------------------------------------
    # ARCHITETTURA ENCODER GAE V3
    # --------------------------------------------------------

    class Encoder(nn.Module):

        def __init__(self):

            super().__init__()

            self.conv1 = GCNConv(
                n_in,
                n_hidden
            )

            self.conv2 = GCNConv(
                n_hidden,
                n_out
            )

        def forward(self, x, edge_index):

            hidden = self.conv1(
                x,
                edge_index
            ).relu()

            return self.conv2(
                hidden,
                edge_index
            )

    encoder = Encoder()

    encoder.load_state_dict(
        state,
        strict=True
    )

    encoder.eval()

    # --------------------------------------------------------
    # NORMALIZZAZIONE FEATURE
    # --------------------------------------------------------

    mean = torch.as_tensor(
        ck["x_mean"],
        dtype=x.dtype
    ).reshape(-1)

    std = torch.as_tensor(
        ck["x_std"],
        dtype=x.dtype
    ).reshape(-1)

    assert mean.numel() == n_in

    assert std.numel() == n_in

    assert bool(torch.isfinite(mean).all())

    assert bool(torch.isfinite(std).all())

    x = (
        x - mean
    ) / std.clamp_min(1e-8)

    assert bool(torch.isfinite(x).all())

    # --------------------------------------------------------
    # INFERENZA SENZA TRAINING
    # --------------------------------------------------------

    with torch.inference_mode():

        z = encoder(
            x,
            edges
        )

    assert tuple(z.shape) == (
        884,
        16
    )

    assert bool(torch.isfinite(z).all())

    return {
        "nodes": int(z.shape[0]),
        "embedding_dim": int(z.shape[1]),
        "edges": int(edges.shape[1]),
        "checkpoint_loaded_strict": True,
        "training_performed": False
    }


check("gae_real_inference", gae)


# ============================================================
# C. CLASSIFICATORE ML ANTILEAKAGE
# ============================================================

IDS = [
    "sRNA_00000002",
    "sRNA_00000040",
    "sRNA_00000007",
    "sRNA_00000435",
    "sRNA_00000037",
]


def ml():

    model_path = (
        p
        / "models/srna_supervised_ml_v1_antileakage/"
        "PRJNA419326_srna_supervised_"
        "antileakage_best_model.joblib"
    )

    table_path = (
        p
        / "outputs/prjna419326/"
        "PRJNA419326_srna_supervised_"
        "ML_antileakage_predictions.tsv"
    )

    with warnings.catch_warnings():

        warnings.simplefilter(
            "error",
            InconsistentVersionWarning
        )

        model = joblib.load(
            model_path
        )

    features = list(
        getattr(
            model,
            "feature_names_in_",
            []
        )
    )

    assert len(features) == 22

    assert sum(
        name.startswith("z_")
        for name in features
    ) == 16

    frame = pd.read_csv(
        table_path,
        sep="\t",
        dtype={"feature_id": str}
    )

    assert all(
        column in frame
        for column in features
    )

    chosen = []

    # --------------------------------------------------------
    # ASSOCIAZIONE DELLE CINQUE RNA
    # --------------------------------------------------------

    for identifier in IDS:

        pattern = (
            r"(?<![A-Za-z0-9])"
            + re.escape(identifier)
            + r"(?![0-9])"
        )

        subset = frame[
            frame.feature_id
            .astype(str)
            .str.contains(
                pattern,
                regex=True,
                na=False
            )
        ]

        assert len(subset) == 1, (
            identifier,
            len(subset)
        )

        chosen.append(
            subset.iloc[0]
        )

    # --------------------------------------------------------
    # PREPARAZIONE FEATURE
    # --------------------------------------------------------

    X = pd.DataFrame(
        chosen
    ).loc[:, features]

    X = X.apply(
        pd.to_numeric,
        errors="raise"
    )

    assert bool(
        np.isfinite(
            X.to_numpy(dtype=float)
        ).all()
    )

    # --------------------------------------------------------
    # INFERENZA REALE
    # --------------------------------------------------------

    predicted = model.predict(X)

    probabilities = model.predict_proba(X)

    assert probabilities.shape[0] == 5

    assert bool(
        np.isfinite(
            probabilities
        ).all()
    )

    assert bool(
        np.allclose(
            probabilities.sum(axis=1),
            1.0,
            atol=1e-6
        )
    )

    historical = [
        str(int(float(value)))
        for value in pd.DataFrame(
            chosen
        )["antileakage_predicted_label"]
    ]

    actual = [
        str(int(float(value)))
        for value in predicted
    ]

    assert historical == actual

    return {
        "ids": IDS,
        "labels": actual,
        "historical_labels": historical,
        "matches_historical": 5,
        "feature_count": 22,
        "independent_validation": False
    }


check("ml_5_historical_srna", ml)


# ============================================================
# D. RISULTATO COMPLESSIVO
# ============================================================

required = (
    "versions",
    "gae_real_inference",
    "ml_5_historical_srna"
)

out["passed"] = all(
    out["tests"].get(name, False)
    for name in required
)

print(
    "PRIMETTA_226G="
    + json.dumps(
        out,
        ensure_ascii=False
    ),
    flush=True
)
'''


# ============================================================
# 3. ESECUZIONE COMANDI
# ============================================================

def cmd(
    args,
    timeout=180,
    output=True,
    check=True
):

    result = subprocess.run(
        args,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False
    )

    if output:

        if result.stdout.strip():

            print(
                result.stdout.strip(),
                flush=True
            )

        if result.stderr.strip():

            print(
                result.stderr[-2500:].strip(),
                file=sys.stderr,
                flush=True
            )

    if check and result.returncode:

        raise RuntimeError(
            "Comando fallito (exit %s): %s"
            % (
                result.returncode,
                " ".join(args)
            )
        )

    return result


# ============================================================
# 4. COMPOSE DOCKER
# ============================================================

def compose(
    include_new=False,
    include_old=True
):

    paths = [
        C1,
        C2
    ]

    if include_old:

        paths.append(C3)

    if include_new:

        paths.append(C4)

    for path in paths:

        if not path.is_file():

            raise FileNotFoundError(
                str(path)
            )

    args = [
        "docker",
        "compose"
    ]

    for path in paths:

        args.extend([
            "-f",
            str(path)
        ])

    return args


# ============================================================
# 5. VERIFICA CONFIGURAZIONE
# ============================================================

def verify_config(
    args,
    expected_image
):

    response = cmd(
        [
            *args,
            "config",
            "--format",
            "json"
        ],
        output=False
    )

    obj = json.loads(
        response.stdout
    )

    service = obj[
        "services"
    ][SERVICE]

    ports = service.get(
        "ports",
        []
    )

    matches = [
        port
        for port in ports
        if (
            str(port.get("published")) == "8503"
            and int(port.get("target", 0)) == 8501
            and port.get("host_ip") == "127.0.0.1"
        )
    ]

    assert len(matches) == 1

    assert not any(
        str(port.get("published")) == "8501"
        for port in ports
    )

    assert service.get("image") == expected_image

    return service


# ============================================================
# 6. UTENTE DELL'IMMAGINE ORIGINALE
# ============================================================

def image_user(image):

    response = cmd(
        [
            "docker",
            "image",
            "inspect",
            image,
            "--format",
            "{{json .Config.User}}"
        ],
        output=False
    )

    user = json.loads(
        response.stdout.strip()
    ) or "root"

    if not re.fullmatch(
        r"[A-Za-z0-9_.:-]+",
        user
    ):

        raise RuntimeError(
            "Utente Docker non riconosciuto"
        )

    return user


# ============================================================
# 7. CREAZIONE FILE SENZA SOVRASCRITTURA
# ============================================================

def create_if_needed(
    path,
    content
):

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if path.exists():

        if path.read_text(
            encoding="utf-8"
        ) != content:

            raise RuntimeError(
                "Esiste un file differente: "
                + str(path)
            )

        print(
            "[OK] File già esistente:",
            path.relative_to(ROOT)
        )

    else:

        path.write_text(
            content,
            encoding="utf-8"
        )

        print(
            "[OK] Creato:",
            path.relative_to(ROOT)
        )


# ============================================================
# 8. GENERAZIONE DOCKERFILE E COMPOSE
# ============================================================

def ensure_files():

    cmd(
        [
            "docker",
            "image",
            "inspect",
            OLD,
            "--format",
            "{{.Id}}"
        ],
        output=False
    )

    user = image_user(OLD)

    dockerfile = (
        f"FROM {OLD}\n"
        "USER root\n"
        "ENV PIP_DISABLE_PIP_VERSION_CHECK=1 "
        "PIP_NO_CACHE_DIR=1\n"
        "RUN python3 -m pip install "
        "--no-cache-dir --no-deps "
        '"scikit-learn==1.8.0"\n'
        f"USER {user}\n"
    )

    override = (
        "services:\n"
        "  maize-ncrna-gnn:\n"
        f"    image: {NEW}\n"
        "    pull_policy: never\n"
    )

    create_if_needed(
        DOCKERFILE,
        dockerfile
    )

    create_if_needed(
        C4,
        override
    )

    verify_config(
        compose(True, True),
        NEW
    )


# ============================================================
# 9. COLLAUDO NUOVA IMMAGINE
# ============================================================

def image_probe():

    response = cmd(
        [
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
            "--memory",
            "3g",
            "--cpus",
            "2",
            "-e",
            "HOME=/tmp",
            "-e",
            "OMP_NUM_THREADS=2",
            "-e",
            "PYTHONDONTWRITEBYTECODE=1",
            "-v",
            f"{ROOT}:/project:ro",
            "--entrypoint",
            "python3",
            NEW,
            "-B",
            "-c",
            WORKER
        ],
        timeout=240,
        output=False
    )

    marker = "PRIMETTA_226G="

    messages = [
        line.split(
            marker,
            1
        )[1]
        for line in response.stdout.splitlines()
        if line.startswith(marker)
    ]

    if (
        response.returncode != 0
        or len(messages) != 1
    ):

        raise RuntimeError(
            "Test immagine non riuscito: "
            + response.stderr[-1500:]
            + response.stdout[-700:]
        )

    detail = json.loads(
        messages[0]
    )

    print(
        "[TEST] Versioni / GAE / ML:",
        json.dumps(
            detail,
            ensure_ascii=False,
            indent=2
        )
    )

    return detail


# ============================================================
# 10. IMMAGINE ATTIVA
# ============================================================

def image_active():

    response = cmd(
        [
            "docker",
            "inspect",
            CONTAINER,
            "--format",
            "{{.Config.Image}}"
        ],
        output=False,
        check=False
    )

    return (
        response.stdout.strip()
        if response.returncode == 0
        else None
    )


# ============================================================
# 11. HEALTHCHECK
# ============================================================

def http(port, path):

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    try:

        with opener.open(
            f"http://127.0.0.1:{port}{path}",
            timeout=5
        ) as response:

            return response.status == 200

    except Exception:

        return False


# ============================================================
# 12. PROTEZIONE FILE SCIENTIFICI
# ============================================================

def hashes():

    protected = [

        ROOT / "app/streamlit_app.py",

        *sorted(
            (ROOT / "app/pages").glob("*.py")
        ),

        ROOT
        / "models/graph_autoencoder_v3/"
        "gae_v3_model.pt",

        ROOT
        / "data/processed/pyg/"
        "maize_ncrna_graph_v3_prjna419326.pt",

        ROOT
        / "models/srna_supervised_ml_v1_antileakage/"
        "PRJNA419326_srna_supervised_"
        "antileakage_best_model.joblib"
    ]

    result = {}

    for path in protected:

        if not path.is_file():

            raise FileNotFoundError(
                str(path)
            )

        checksum = hashlib.sha256()

        with path.open("rb") as handle:

            for block in iter(
                lambda: handle.read(
                    1024 * 1024
                ),
                b""
            ):

                checksum.update(block)

        result[
            str(path.relative_to(ROOT))
        ] = checksum.hexdigest()

    return result


# ============================================================
# 13. RIPRISTINO IMMAGINE PRECEDENTE
# ============================================================

def rollback():

    print(
        "[RIPRISTINO] Ritorno "
        "all'immagine 22.6d."
    )

    verify_config(
        compose(False, True),
        OLD
    )

    cmd(
        [
            *compose(False, True),
            "up",
            "-d",
            "--no-deps",
            "--no-build",
            "--pull",
            "never",
            "--force-recreate",
            SERVICE
        ]
    )

    print(
        "[RIPRISTINO] Immagine attiva:",
        image_active()
    )


# ============================================================
# 14. STATO APPLICAZIONI
# ============================================================

def status():

    print(
        "[STATO] Immagine V3.2:",
        image_active()
    )

    endpoints = [

        (
            "PRIMETTA V3.2",
            8503,
            "/_stcore/health"
        ),

        (
            "PRIMETTA V3.3",
            8502,
            "/_stcore/health"
        ),

        (
            "AI ORZO",
            8501,
            "/_stcore/health"
        ),

        (
            "API V3.2",
            8000,
            "/health"
        )
    ]

    for name, port, path in endpoints:

        print(
            "[OK]"
            if http(port, path)
            else "[NON VERIFICATO]",
            name,
            port
        )

    cmd(
        [
            "docker",
            "ps",
            "--format",
            "table {{.Names}}\t{{.Status}}\t{{.Ports}}"
        ],
        check=False
    )


# ============================================================
# 15. INSTALLAZIONE E ATTIVAZIONE
# ============================================================

def install():

    print(
        "[1] Verifica configurazione e file originali"
    )

    verify_config(
        compose(False, True),
        OLD
    )

    before = hashes()

    ensure_files()

    # --------------------------------------------------------
    # COSTRUZIONE IMMAGINE
    # --------------------------------------------------------

    print(
        "[2] Costruzione immagine "
        "con scikit-learn 1.8.0"
    )

    cmd(
        [
            "docker",
            "build",
            "-f",
            str(DOCKERFILE),
            "-t",
            NEW,
            str(DOCKERFILE.parent)
        ],
        timeout=1800
    )

    # --------------------------------------------------------
    # COLLAUDO GAE E ML
    # --------------------------------------------------------

    print(
        "[3] Test reale GAE + ML "
        "senza training"
    )

    proof = image_probe()

    if not proof.get("passed"):

        print(
            "[STOP] Immagine NON attivata. "
            "Esaminare il rapporto del test."
        )

        return False

    if before != hashes():

        raise RuntimeError(
            "File originali modificati. "
            "Attivazione bloccata."
        )

    # --------------------------------------------------------
    # ATTIVAZIONE GUI V3.2
    # --------------------------------------------------------

    print(
        "[4] Attivazione del solo "
        "container GUI V3.2"
    )

    try:

        cmd(
            [
                *compose(True, True),
                "up",
                "-d",
                "--no-deps",
                "--no-build",
                "--pull",
                "never",
                "--force-recreate",
                SERVICE
            ],
            timeout=180
        )

        ready = False

        for _ in range(30):

            if (
                image_active() == NEW
                and http(
                    8503,
                    "/_stcore/health"
                )
                and http(
                    8000,
                    "/health"
                )
            ):

                ready = True

                break

            time.sleep(2)

        if not ready:

            raise RuntimeError(
                "GUI V3.2 o API "
                "non raggiungibili."
            )

        # ----------------------------------------------------
        # VERIFICA VERSIONI NEL CONTAINER ATTIVO
        # ----------------------------------------------------

        cmd(
            [
                "docker",
                "exec",
                CONTAINER,
                "python3",
                "-B",
                "-c",
                (
                    "import sklearn,torch,"
                    "torch_geometric,networkx;"
                    "assert sklearn.__version__=='1.8.0';"
                    "print(sklearn.__version__)"
                )
            ]
        )

        if before != hashes():

            raise RuntimeError(
                "File originali modificati "
                "durante l'attivazione."
            )

    except Exception:

        rollback()

        raise

    # --------------------------------------------------------
    # RAPPORTO FINALE
    # --------------------------------------------------------

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    destination = RUNS / (
        "STEP22_6G_V32_FIX_" + stamp
    )

    destination.mkdir(
        parents=True,
        exist_ok=False
    )

    report_path = (
        destination / "RAPPORTO_STEP22_6G.json"
    )

    report_path.write_text(
        json.dumps(
            {
                "status": (
                    "RUNTIME_V32_SKLEARN18_GAE_ML_PASS"
                ),
                "image": NEW,
                "proof": proof,
                "protected_sha256": before,
                "gui_page_functional_tested": False,
                "rag_question_tested": False,
                "training_performed": False
            },
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print()
    print("=" * 68)

    print(
        "[OK] NUOVO RUNTIME GUI V3.2 ATTIVO"
    )

    print(
        "[OK] scikit-learn 1.8.0"
    )

    print(
        "[OK] Inferenza GAE verificata"
    )

    print(
        "[OK] Inferenza ML verificata"
    )

    print(
        "[OK] Nessun training"
    )

    print(
        "[REPORT]",
        report_path
    )

    print("=" * 68)

    status()

    return True


# ============================================================
# 16. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    options = parser.add_mutually_exclusive_group(
        required=True
    )

    for flag in (
        "install",
        "status",
        "rollback"
    ):

        options.add_argument(
            "--" + flag,
            action="store_true"
        )

    args = parser.parse_args()

    print("=" * 68)

    print(
        "AI PRIMETTA V3.2 - STEP 22.6G"
    )

    print(
        "COMPATIBILITA SCIKIT-LEARN + GAE + ML"
    )

    print("=" * 68)

    if args.install:

        if not install():

            sys.exit(2)

    elif args.status:

        status()

    else:

        rollback()

        status()


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        sys.exit(130)

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__,
            str(exc),
            file=sys.stderr
        )

        sys.exit(1)