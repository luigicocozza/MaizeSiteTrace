#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.2 - STEP 22.6E

Diagnosi runtime PyTorch/PyG e attivazione sicura
dell'immagine Docker gia costruita.

Comandi:

python3 scripts/step_22_6e_verifica_attiva_gui_v32.py --install
python3 scripts/step_22_6e_verifica_attiva_gui_v32.py --status
python3 scripts/step_22_6e_verifica_attiva_gui_v32.py --rollback

Non installa pacchetti.
Non esegue training scientifico.
Non modifica GUI, checkpoint e dati.
Non modifica AI ORZO o la GUI V3.3.
"""

import argparse
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

COMPOSE = ROOT / "docker-compose.yml"

PORT_OVERRIDE = (
    ROOT / "docker-compose.primetta-v32-8503.yaml"
)

AI_OVERRIDE = (
    ROOT / "docker-compose.primetta-v32-full-ai.yaml"
)

IMAGE = "ai_primetta_v32_full_ai_cpu:22.6d"

SERVICE = "maize-ncrna-gnn"

CONTAINER = "maize_ncrna_gnn"

RUNS = (
    ROOT / "target_prediction_v33_vscode" / "runs"
)

REPORT_NAME = "RAPPORTO_STEP22_6E.json"


# ============================================================
# 2. COLLAUDO DELLE LIBRERIE E DEGLI ARTEFATTI
# ============================================================

SMOKE = r'''

import json

from pathlib import Path


report = {
    "tests": {},
    "artifacts": {},
    "errors": []
}


# ------------------------------------------------------------
# A. IMPORT LIBRERIE
# ------------------------------------------------------------

try:

    import torch
    import torch_geometric
    import networkx
    import streamlit
    import sklearn
    import joblib

    from torch_geometric.nn import GCNConv

    report["versions"] = {
        "torch": torch.__version__,
        "torch_geometric": torch_geometric.__version__,
        "networkx": networkx.__version__,
        "streamlit": streamlit.__version__,
        "sklearn": sklearn.__version__
    }

    report["tests"]["dependencies_import"] = True

except Exception as exc:

    report["tests"]["dependencies_import"] = False

    report["errors"].append(
        "IMPORT: " + repr(exc)
    )


# ------------------------------------------------------------
# B. PROVA GCN SU TRE NODI SINTETICI
# ------------------------------------------------------------

if report["tests"].get("dependencies_import"):

    try:

        x = torch.tensor(
            [
                [1., 0.],
                [0., 1.],
                [1., 1.]
            ],
            requires_grad=True
        )

        edges = torch.tensor(
            [
                [0, 1, 2, 2],
                [1, 2, 0, 1]
            ],
            dtype=torch.long
        )

        layer = GCNConv(
            2,
            4
        )

        output = layer(
            x,
            edges
        )

        assert tuple(
            output.shape
        ) == (3, 4), "GCN output shape"

        # Verifica del funzionamento autograd
        # su dati sintetici, non scientifici.
        output.square().mean().backward()

        assert (
            layer.lin.weight.grad is not None
        ), "GCN gradient"

        report["tests"][
            "gcn_forward_backward_synthetic"
        ] = True

    except Exception as exc:

        report["tests"][
            "gcn_forward_backward_synthetic"
        ] = False

        report["errors"].append(
            "GCN: " + repr(exc)
        )


    # --------------------------------------------------------
    # C. CARICAMENTO GRAFO E CHECKPOINT ORIGINALI
    # --------------------------------------------------------

    base = Path("/project")

    paths = {

        "graph": (
            base
            / "data/processed/pyg/"
            "maize_ncrna_graph_v3_prjna419326.pt"
        ),

        "checkpoint": (
            base
            / "models/graph_autoencoder_v3/"
            "gae_v3_model.pt"
        )
    }

    for name, path in paths.items():

        detail = {
            "exists": path.is_file(),
            "path": str(path)
        }

        if path.is_file():

            try:

                # Solo artefatti locali attendibili.
                obj = torch.load(
                    path,
                    map_location="cpu",
                    weights_only=False
                )

                detail["load_ok"] = True

                detail["python_type"] = (
                    type(obj).__module__
                    + "."
                    + type(obj).__name__
                )

                if isinstance(obj, dict):

                    detail["keys_preview"] = [
                        str(k)
                        for k in list(obj.keys())[:25]
                    ]

                else:

                    try:

                        detail["keys_preview"] = [
                            str(k)
                            for k in list(obj.keys())[:25]
                        ]

                    except (
                        AttributeError,
                        TypeError
                    ):
                        pass

                # Le informazioni vengono lette solo
                # se esistono realmente nell'oggetto.

                for field in (
                    "x",
                    "edge_index"
                ):

                    value = getattr(
                        obj,
                        field,
                        None
                    )

                    if isinstance(
                        value,
                        torch.Tensor
                    ):

                        detail[field + "_shape"] = list(
                            value.shape
                        )

                try:

                    detail["num_nodes"] = int(
                        obj.num_nodes
                    )

                except (
                    AttributeError,
                    ValueError,
                    TypeError
                ):
                    pass

                del obj

            except Exception as exc:

                detail["load_ok"] = False

                detail["load_error"] = repr(exc)

                report["errors"].append(
                    name.upper()
                    + ": "
                    + repr(exc)
                )

        else:

            detail["load_ok"] = False

            report["errors"].append(
                name.upper() + ": file assente"
            )

        report["artifacts"][name] = detail

        report["tests"][
            name + "_loaded"
        ] = bool(detail["load_ok"])


# ------------------------------------------------------------
# D. ESITO DEL COLLAUDO
# ------------------------------------------------------------

required = (
    "dependencies_import",
    "gcn_forward_backward_synthetic",
    "graph_loaded",
    "checkpoint_loaded"
)

report["passed"] = all(
    report["tests"].get(name, False)
    for name in required
)

print(
    "PRIMETTA_RUNTIME_DIAG="
    + json.dumps(
        report,
        ensure_ascii=False
    ),
    flush=True
)

'''


# ============================================================
# 3. ESECUZIONE COMANDI
# ============================================================

def invoke(
    args,
    timeout=120,
    check=True
):

    print(
        "$",
        " ".join(str(x) for x in args),
        flush=True
    )

    proc = subprocess.run(
        args,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False
    )

    if proc.stdout.strip():

        print(
            proc.stdout.strip(),
            flush=True
        )

    if proc.stderr.strip():

        print(
            proc.stderr.strip()[-2400:],
            file=sys.stderr,
            flush=True
        )

    if check and proc.returncode:

        raise RuntimeError(
            "Comando non riuscito "
            f"(exit {proc.returncode})"
        )

    return proc


# ============================================================
# 4. CONFIGURAZIONE DOCKER COMPOSE
# ============================================================

def compose(full):

    files = [
        COMPOSE,
        PORT_OVERRIDE
    ]

    if full:

        files.append(
            AI_OVERRIDE
        )

    missing = [
        str(p)
        for p in files
        if not p.is_file()
    ]

    if missing:

        raise FileNotFoundError(
            "Compose mancante: "
            + ", ".join(missing)
        )

    args = [
        "docker",
        "compose"
    ]

    for p in files:

        args.extend([
            "-f",
            str(p)
        ])

    return args


# ============================================================
# 5. CONTROLLO SERVIZIO E PORTE
# ============================================================

def effective_service(full):

    result = invoke([
        *compose(full),
        "config",
        "--format",
        "json"
    ])

    data = json.loads(
        result.stdout
    )

    cfg = data[
        "services"
    ][SERVICE]

    ports = cfg.get(
        "ports",
        []
    )

    expected = [
        p
        for p in ports
        if (
            str(p.get("published")) == "8503"
            and int(p.get("target", -1)) == 8501
            and str(p.get("host_ip")) == "127.0.0.1"
        )
    ]

    forbidden = [
        p
        for p in ports
        if str(
            p.get("published")
        ) == "8501"
    ]

    if len(expected) != 1 or forbidden:

        raise RuntimeError(
            "Configurazione porte non sicura. "
            "GUI V3.2 deve utilizzare "
            "127.0.0.1:8503 -> 8501."
        )

    if full and cfg.get("image") != IMAGE:

        raise RuntimeError(
            "Override AI non punta "
            "all'immagine prevista"
        )

    return cfg


# ============================================================
# 6. IMMAGINE ATTUALMENTE ATTIVA
# ============================================================

def active_image():

    result = invoke(
        [
            "docker",
            "inspect",
            CONTAINER,
            "--format",
            "{{.Config.Image}}"
        ],
        check=False
    )

    if result.returncode == 0:

        return result.stdout.strip()

    return None


# ============================================================
# 7. HEALTHCHECK HTTP
# ============================================================

def http_check(port, path):

    url = (
        f"http://127.0.0.1:{port}{path}"
    )

    try:

        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            url,
            timeout=4
        ) as response:

            return response.status == 200

    except Exception:

        return False


# ============================================================
# 8. INTEGRITÀ SHA256
# ============================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b""
        ):

            h.update(block)

    return h.hexdigest()


def originals():

    files = [
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
        "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
    ]

    missing = [
        str(p)
        for p in files
        if not p.is_file()
    ]

    if missing:

        raise RuntimeError(
            "Originali mancanti: "
            + ", ".join(missing)
        )

    return {
        str(p.relative_to(ROOT)): sha256(p)
        for p in files
    }


# ============================================================
# 9. COLLAUDO IMMAGINE AI GIÀ COSTRUITA
# ============================================================

def image_test():

    invoke([
        "docker",
        "image",
        "inspect",
        IMAGE,
        "--format",
        "{{.Id}}"
    ])

    result = invoke(
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

            "--cpus",
            "2",

            "--memory",
            "3g",

            "-e",
            "CUDA_VISIBLE_DEVICES=-1",

            "-e",
            "OMP_NUM_THREADS=2",

            "-v",
            f"{ROOT}:/project:ro",

            "--entrypoint",
            "python3",

            IMAGE,

            "-B",
            "-c",
            SMOKE
        ],
        timeout=240
    )

    marker = "PRIMETTA_RUNTIME_DIAG="

    messages = [
        x[len(marker):]
        for x in result.stdout.splitlines()
        if x.startswith(marker)
    ]

    if len(messages) != 1:

        raise RuntimeError(
            "Test privo di un rapporto "
            "diagnostico valido"
        )

    report = json.loads(
        messages[0]
    )

    print(
        "\n[DIAGNOSTICA] Versioni:",
        report.get("versions", {})
    )

    print(
        "[DIAGNOSTICA] Verifiche:",
        report["tests"]
    )

    for name, detail in report[
        "artifacts"
    ].items():

        print(
            "[DIAGNOSTICA]",
            name,
            detail
        )

    if report["errors"]:

        print(
            "[DIAGNOSTICA] Errori:",
            report["errors"]
        )

    return report


# ============================================================
# 10. ROLLBACK DEL SOLO SERVIZIO V3.2
# ============================================================

def rollback():

    print(
        "\n[ROLLBACK] Ripristino immagine "
        "originale GUI V3.2"
    )

    effective_service(False)

    invoke(
        [
            *compose(False),

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

    print(
        "[ROLLBACK] Immagine attiva:",
        active_image()
    )


# ============================================================
# 11. SALVATAGGIO RAPPORTO
# ============================================================

def report_to_file(payload):

    RUNS.mkdir(
        parents=True,
        exist_ok=True
    )

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP22_6E_RUNTIME_V32_" + stamp
    )

    folder.mkdir(
        exist_ok=False
    )

    path = folder / REPORT_NAME

    path.write_text(
        json.dumps(
            payload,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    print(
        "[REPORT]",
        path
    )


# ============================================================
# 12. STATO GENERALE
# ============================================================

def status():

    print(
        "\n[STATO] Immagine attiva:",
        active_image()
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

        ok = http_check(
            port,
            path
        )

        print(
            "[OK]" if ok else "[NON VERIFICATO]",
            name,
            port
        )

    invoke(
        [
            "docker",
            "ps",
            "--format",
            "table {{.Names}}\t{{.Status}}\t{{.Ports}}"
        ],
        check=False
    )


# ============================================================
# 13. INSTALLAZIONE E ATTIVAZIONE
# ============================================================

def install():

    # --------------------------------------------------------
    # A. VERIFICA CONFIGURAZIONI
    # --------------------------------------------------------

    effective_service(False)

    effective_service(True)

    # --------------------------------------------------------
    # B. SHA256 PRIMA DEI TEST
    # --------------------------------------------------------

    print(
        "\n[1] Controllo SHA256 originali"
    )

    before = originals()

    # --------------------------------------------------------
    # C. COLLAUDO IMMAGINE ESISTENTE
    # --------------------------------------------------------

    print(
        "\n[2] Diagnosi immagine AI "
        "gia costruita"
    )

    diagnostic = image_test()

    if not diagnostic["passed"]:

        print(
            "[STOP] Diagnostica non superata. "
            "GUI V3.2 lasciata invariata."
        )

        report_to_file(
            {
                "step": "22.6E",

                "status": (
                    "DIAGNOSTIC_FAILED_NO_ACTIVATION"
                ),

                "diagnostic": diagnostic,

                "image_active": active_image(),

                "training_performed": False
            }
        )

        sys.exit(2)

    # --------------------------------------------------------
    # D. VERIFICA FILE ORIGINALI
    # --------------------------------------------------------

    if originals() != before:

        raise RuntimeError(
            "Checksum originali cambiati. "
            "Attivazione bloccata."
        )

    # --------------------------------------------------------
    # E. ATTIVAZIONE NUOVA IMMAGINE
    # --------------------------------------------------------

    print(
        "\n[3] Attivazione immagine AI "
        "per il solo servizio V3.2"
    )

    try:

        invoke(
            [
                *compose(True),

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

        healthy = False

        for _ in range(30):

            if (
                active_image() == IMAGE
                and http_check(
                    8503,
                    "/_stcore/health"
                )
                and http_check(
                    8000,
                    "/health"
                )
            ):

                healthy = True
                break

            time.sleep(2)

        if not healthy:

            raise RuntimeError(
                "GUI/API non raggiungibili "
                "o immagine attiva errata"
            )

        # Controllo reale anche nel container
        # utilizzato dall'applicazione.

        invoke(
            [
                "docker",
                "exec",
                CONTAINER,
                "python3",
                "-B",
                "-c",
                (
                    "import torch,"
                    "torch_geometric,networkx;"
                    "print(torch.__version__,"
                    "torch_geometric.__version__,"
                    "networkx.__version__)"
                )
            ],
            timeout=40
        )

        if originals() != before:

            raise RuntimeError(
                "File originali modificati "
                "durante l'attivazione"
            )

    except Exception:

        rollback()

        raise

    # --------------------------------------------------------
    # F. STATO FINALE
    # --------------------------------------------------------

    print(
        "\n[OK] Nuova immagine AI attiva"
    )

    status()

    report_to_file(
        {
            "step": "22.6E",

            "status": "RUNTIME_AI_ACTIVATED",

            "diagnostic": diagnostic,

            "active_image": active_image(),

            "original_sha256": before,

            "training_performed": False,

            "full_functional_gui_validation": False,

            "note": (
                "Librerie e caricamento artefatti "
                "verificati. Rimangono da collaudare "
                "i pulsanti GAE, upload e Chat RAG."
            )
        }
    )


# ============================================================
# 14. AVVIO PRINCIPALE
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

    print("=" * 72)

    print(
        "AI PRIMETTA V3.2 - STEP 22.6E"
    )

    print(
        "CORREZIONE TEST E ATTIVAZIONE"
    )

    print("=" * 72)

    if args.install:

        install()

    elif args.rollback:

        rollback()
        status()

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