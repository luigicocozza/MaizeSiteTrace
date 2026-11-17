#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA - STEP 22.6F

COLLAUDO FUNZIONALE DEL RUNTIME GUI V3.2

Verifica:
- Immagine Docker effettivamente attiva.
- Pulsanti e riferimenti del modulo GAE.
- Librerie PyTorch/PyG/NetworkX.
- Operazione GCN con autograd.
- Struttura del grafo GAE V3.
- Metadati e struttura del checkpoint.
- Compatibilita checkpoint ML/scikit-learn.
- Stato delle due GUI, di AI ORZO e dell'API.
- Endpoint FastAPI disponibili.

NON avvia training scientifici.
NON modifica GUI, modelli o dati originali.
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
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

CONTAINER = "maize_ncrna_gnn"

IMAGE = "ai_primetta_v32_full_ai_cpu:22.6d"

PAGE = (
    ROOT
    / "app/pages/01_Rete_neurale_GAE.py"
)


# ============================================================
# 2. TEST ESEGUITO NEL CONTAINER V3.2
# ============================================================

PROBE = r'''

import json
import warnings

from pathlib import Path


result = {

    "tests": {},

    "details": {},

    "warnings": []
}


# ------------------------------------------------------------
# A. IMPORT LIBRERIE
# ------------------------------------------------------------

try:

    import torch
    import torch_geometric
    import networkx
    import sklearn
    import joblib

    from sklearn.exceptions import (
        InconsistentVersionWarning
    )

    from torch_geometric.nn import GCNConv

    result["details"]["versions"] = {

        "torch": torch.__version__,

        "pyg": torch_geometric.__version__,

        "networkx": networkx.__version__,

        "sklearn": sklearn.__version__
    }

    result["tests"]["imports"] = True

except Exception as exc:

    result["tests"]["imports"] = False

    result["warnings"].append(
        "Import: " + repr(exc)
    )


# ============================================================
# B. OPERAZIONE REALE GCN SU DATI SINTETICI
# ============================================================

if result["tests"]["imports"]:

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

        y = layer(
            x,
            edges
        )

        loss = y.square().mean()

        loss.backward()

        result["tests"]["gcn_operation"] = (
            tuple(y.shape) == (3, 4)
            and layer.lin.weight.grad is not None
        )

    except Exception as exc:

        result["tests"]["gcn_operation"] = False

        result["warnings"].append(
            "GCN: " + repr(exc)
        )

    # ========================================================
    # C. STRUTTURA GRAFO GAE V3
    # ========================================================

    try:

        graph_path = Path(
            "/app/data/processed/pyg/"
            "maize_ncrna_graph_v3_prjna419326.pt"
        )

        # Solo artefatto locale attendibile.
        raw = torch.load(
            graph_path,
            map_location="cpu",
            weights_only=False
        )

        graph = (
            raw["data"]
            if isinstance(raw, dict)
            else raw
        )

        assert hasattr(
            graph,
            "x"
        )

        assert hasattr(
            graph,
            "edge_index"
        )

        assert graph.x.ndim == 2

        assert graph.edge_index.shape[0] == 2

        assert bool(
            torch.isfinite(
                graph.x
            ).all()
        )

        result["details"]["graph"] = {

            "nodes": int(
                graph.num_nodes
            ),

            "edges": int(
                graph.edge_index.shape[1]
            ),

            "features": int(
                graph.x.shape[1]
            ),

            "metadata_type": (
                type(
                    raw.get("metadata")
                ).__name__
                if isinstance(raw, dict)
                else None
            )
        }

        result["tests"]["graph_structure"] = True

    except Exception as exc:

        result["tests"]["graph_structure"] = False

        result["warnings"].append(
            "Grafo: " + repr(exc)
        )

    # ========================================================
    # D. STRUTTURA CHECKPOINT GAE
    # ========================================================

    try:

        checkpoint_path = Path(
            "/app/models/graph_autoencoder_v3/"
            "gae_v3_model.pt"
        )

        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=False
        )

        assert isinstance(
            checkpoint,
            dict
        )

        assert checkpoint.get(
            "model_state_dict"
        )

        detail = {

            "model_type": str(
                checkpoint.get(
                    "model_type"
                )
            ),

            "in_channels": checkpoint.get(
                "in_channels"
            ),

            "hidden_channels": checkpoint.get(
                "hidden_channels"
            ),

            "latent_channels": checkpoint.get(
                "latent_channels"
            ),

            "state_dict_entries": len(
                checkpoint["model_state_dict"]
            )
        }

        result["details"]["checkpoint"] = detail

        result["tests"]["checkpoint_structure"] = True

        if result["tests"].get(
            "graph_structure"
        ):

            result["tests"][
                "feature_dimension_matches"
            ] = (
                int(
                    detail["in_channels"]
                )
                ==
                result["details"]["graph"]["features"]
            )

    except Exception as exc:

        result["tests"]["checkpoint_structure"] = False

        result["warnings"].append(
            "Checkpoint: " + repr(exc)
        )

    # ========================================================
    # E. COMPATIBILITA CLASSIFICATORE ANTILEAKAGE
    # ========================================================

    try:

        model_path = Path(
            "/app/models/"
            "srna_supervised_ml_v1_antileakage/"
            "PRJNA419326_srna_supervised_"
            "antileakage_best_model.joblib"
        )

        with warnings.catch_warnings():

            warnings.simplefilter(
                "error",
                InconsistentVersionWarning
            )

            model = joblib.load(
                model_path
            )

        result["tests"][
            "antileakage_model_load_compatible"
        ] = True

        result["details"][
            "ml_feature_count"
        ] = len(
            getattr(
                model,
                "feature_names_in_",
                []
            )
        )

    except Exception as exc:

        result["tests"][
            "antileakage_model_load_compatible"
        ] = False

        result["warnings"].append(
            "Compatibilita sklearn/joblib: "
            + repr(exc)
        )


# ============================================================
# F. RISULTATO DEL TEST
# ============================================================

print(
    "PRIMETTA_GAE_FUNCTIONAL="
    + json.dumps(
        result,
        ensure_ascii=False,
        default=str
    ),
    flush=True
)

'''


# ============================================================
# 3. FUNZIONI DI SUPPORTO
# ============================================================

def command(args, timeout=120):

    process = subprocess.run(
        args,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False
    )

    return {

        "ok": process.returncode == 0,

        "stdout": process.stdout,

        "stderr": process.stderr[-1600:],

        "returncode": process.returncode
    }


# ============================================================
# 4. CONTROLLO HTTP
# ============================================================

def health(port, path):

    url = (
        f"http://127.0.0.1:{port}{path}"
    )

    try:

        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            url,
            timeout=6
        ) as response:

            return response.status == 200

    except Exception:

        return False


# ============================================================
# 5. CONTROLLO FASTAPI
# ============================================================

def check_api():

    try:

        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            "http://127.0.0.1:8000/openapi.json",
            timeout=8
        ) as response:

            specification = json.load(
                response
            )

        paths = specification.get(
            "paths",
            {}
        )

        return {

            "ok": True,

            "paths": sorted(
                paths
            ),

            "rag_route_declared": (
                "/rag/ask" in paths
            )
        }

    except Exception as exc:

        return {
            "ok": False,
            "error": repr(exc)
        }


# ============================================================
# 6. ANALISI DEL PULSANTE GAE
# ============================================================

def scan_gae_button():

    source = PAGE.read_text(
        encoding="utf-8"
    )

    lines = source.splitlines()

    tree = ast.parse(
        source
    )

    result = {

        "file": str(
            PAGE.relative_to(ROOT)
        ),

        "sha256": hashlib.sha256(
            PAGE.read_bytes()
        ).hexdigest(),

        "buttons": [],

        "training_keywords": []
    }

    for node in ast.walk(tree):

        if not isinstance(
            node,
            ast.Call
        ):
            continue

        try:

            function = ast.unparse(
                node.func
            )

        except Exception:
            continue

        if function in (
            "st.button",
            "st.form_submit_button"
        ):

            label = (
                ast.unparse(
                    node.args[0]
                )
                if node.args
                else "<dinamico>"
            )

            result["buttons"].append({

                "line": node.lineno,

                "label": label
            })

    # Ricerca collegamenti verso backend e training.

    keywords = (
        "subprocess",
        "train_",
        "torch.load",
        "gae_v3",
        "requests.post",
        "api_url"
    )

    for line_number, line in enumerate(
        lines,
        start=1
    ):

        if any(
            keyword in line.lower()
            for keyword in keywords
        ):

            result[
                "training_keywords"
            ].append({

                "line": line_number,

                "text": line.strip()[:180]
            })

    return result


# ============================================================
# 7. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    print("=" * 68)

    print(
        "AI PRIMETTA - STEP 22.6F"
    )

    print(
        "COLLAUDO FUNZIONALE SENZA TRAINING"
    )

    print("=" * 68)

    if not PAGE.is_file():

        raise FileNotFoundError(
            PAGE
        )

    # --------------------------------------------------------
    # 8. IMMAGINE ATTIVA
    # --------------------------------------------------------

    inspect = command([
        "docker",
        "inspect",
        CONTAINER,
        "--format",
        "{{.Config.Image}}"
    ])

    active_image = (
        inspect["stdout"].strip()
        if inspect["ok"]
        else ""
    )

    print(
        "\n[1] IMMAGINE ATTIVA:",
        active_image or "NON DISPONIBILE"
    )

    if active_image != IMAGE:

        raise RuntimeError(
            "La GUI V3.2 non utilizza "
            "l'immagine prevista."
        )

    # --------------------------------------------------------
    # 9. PULSANTI DELLA GUI GAE
    # --------------------------------------------------------

    ui = scan_gae_button()

    print(
        "\n[2] PULSANTI DEL MODULO GAE"
    )

    for button in ui["buttons"]:

        print(
            "Riga",
            button["line"],
            "|",
            button["label"]
        )

    print(
        "Riferimenti potenzialmente "
        "legati al training:",
        len(
            ui["training_keywords"]
        )
    )

    # --------------------------------------------------------
    # 10. PROVA NEL CONTAINER ATTIVO
    # --------------------------------------------------------

    print(
        "\n[3] VERIFICA RUNTIME ML/GNN"
    )

    response = command(
        [
            "docker",
            "exec",
            CONTAINER,

            "python3",
            "-B",
            "-c",
            PROBE
        ],
        timeout=180
    )

    if not response["ok"]:

        raise RuntimeError(
            "Test Docker fallito: "
            + response["stderr"]
        )

    marker = "PRIMETTA_GAE_FUNCTIONAL="

    messages = [
        line[len(marker):]

        for line in response[
            "stdout"
        ].splitlines()

        if line.startswith(marker)
    ]

    if len(messages) != 1:

        raise RuntimeError(
            "Risposta diagnostica non valida: "
            + response["stdout"][-1000:]
        )

    result = json.loads(
        messages[0]
    )

    for name, passed in result[
        "tests"
    ].items():

        print(
            "[OK]" if passed else "[DA RISOLVERE]",
            name
        )

    print(
        "\nDettagli:",
        json.dumps(
            result["details"],
            ensure_ascii=False
        )
    )

    for warning in result[
        "warnings"
    ]:

        print(
            "[NOTA]",
            warning
        )

    # --------------------------------------------------------
    # 11. CONTROLLO GUI E API
    # --------------------------------------------------------

    print(
        "\n[4] CONTROLLO GUI E API"
    )

    endpoints = {

        "PRIMETTA_V32_8503":
            health(
                8503,
                "/_stcore/health"
            ),

        "PRIMETTA_V33_8502":
            health(
                8502,
                "/_stcore/health"
            ),

        "AI_ORZO_8501":
            health(
                8501,
                "/_stcore/health"
            ),

        "API_V32_8000":
            health(
                8000,
                "/health"
            )
    }

    for name, passed in endpoints.items():

        print(
            "[OK]" if passed else "[DA RISOLVERE]",
            name
        )

    api = check_api()

    print(
        "OpenAPI:",
        (
            "DISPONIBILE"
            if api["ok"]
            else "NON DISPONIBILE"
        )
    )

    print(
        "Endpoint RAG dichiarato:",
        api.get(
            "rag_route_declared",
            False
        )
    )

    # --------------------------------------------------------
    # 12. RAPPORTO DI COLLAUDO
    # --------------------------------------------------------

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP22_6F_COLLAUDO_V32_"
        + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {

        "step": "22.6F",

        "timestamp_utc": stamp,

        "image": active_image,

        "gae_ui_static_audit": ui,

        "docker_runtime": result,

        "http": endpoints,

        "api_openapi": api,

        "training_performed": False,

        "gui_modified": False,

        "actual_gae_page_training_button_tested": False,

        "rag_chat_answer_tested": False,

        "note": (
            "Il caricamento del checkpoint "
            "e il test GCN non dimostrano "
            "ancora il percorso del pulsante GAE."
        )
    }

    output = (
        folder
        / "RAPPORTO_STEP22_6F.json"
    )

    output.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    print(
        "\n[REPORT]",
        output
    )

    print(
        "[OK] Nessun training"
    )

    print(
        "[OK] Nessuna modifica alle GUI"
    )

    print(
        "[LIMITE] Da verificare "
        "il pulsante GAE e una risposta RAG"
    )


# ============================================================
# 13. AVVIO
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