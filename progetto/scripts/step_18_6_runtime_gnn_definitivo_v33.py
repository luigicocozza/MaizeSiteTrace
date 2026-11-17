#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 18.6
Runtime ML/GNN CPU separato e verificato.

- Non modifica il Dockerfile originale
- Non modifica la GUI V3.2
- Non modifica il motore V3.3
- Non carica o addestra checkpoint
- Crea una nuova immagine Docker dedicata
- Esegue uno smoke test CPU di PyTorch/PyG
- Salva report e Dockerfile per riproducibilita
"""

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

RUNS = ROOT / "target_prediction_v33_vscode" / "runs"

SOURCE_IMAGE = (
    "ncrna_rnacentral_ai_v3_prjna419326-"
    "maize-ncrna-gnn:latest"
)

NEW_IMAGE = "ai_primetta_v33_gnn_cpu:18.6"

DOCKERFILE = f"""
FROM {SOURCE_IMAGE}

USER root

RUN python3 -m pip install --no-cache-dir \\
    torch==2.11.0 \\
    --index-url https://download.pytorch.org/whl/cpu

RUN python3 -m pip install --no-cache-dir \\
    torch-geometric==2.8.1 \\
    networkx

ENV OMP_NUM_THREADS=2
ENV MKL_NUM_THREADS=2
ENV PYTHONDONTWRITEBYTECODE=1
"""

TEST_CODE = r"""
import importlib
import json
import sys

packages = [
    "numpy",
    "pandas",
    "sklearn",
    "torch",
    "torch_geometric",
    "joblib",
    "networkx",
]

versions = {}

for name in packages:
    module = importlib.import_module(name)
    versions[name] = getattr(
        module, "__version__", "non dichiarata"
    )
    print("[OK]", name, versions[name], flush=True)

import torch
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv, SAGEConv

torch.set_num_threads(2)
torch.manual_seed(42)

x = torch.tensor([
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
], dtype=torch.float32)

edge_index = torch.tensor([
    [0, 1, 2, 0],
    [1, 2, 0, 2],
], dtype=torch.long)

graph = Data(x=x, edge_index=edge_index)

assert graph.num_nodes == 3

for model in (
    GCNConv(2, 4),
    SAGEConv(2, 4),
):
    model.eval()

    with torch.no_grad():
        output = model(
            graph.x,
            graph.edge_index,
        )

    assert output.shape == (3, 4)
    assert torch.isfinite(output).all().item()

    print(
        "[OK]",
        type(model).__name__,
        "output:",
        tuple(output.shape),
        flush=True,
    )

print(
    "RUNTIME_TEST_RESULT=" +
    json.dumps({
        "status": "PASS",
        "python": sys.version.split()[0],
        "torch_cpu": torch.version.cuda is None,
        "packages": versions,
        "gnn_tests": ["GCNConv", "SAGEConv"],
        "inference_checkpoint_tested": False,
        "training_performed": False,
    }),
    flush=True,
)
"""


def run(command, capture=False):
    result = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
        check=False,
    )

    if result.returncode != 0:
        if capture and result.stdout:
            print(result.stdout)
        raise RuntimeError(
            "Comando fallito: " + " ".join(command)
        )

    return result.stdout if capture else ""


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 18.6")
    print("BUILD RUNTIME ML/GNN CPU")
    print("=" * 70, flush=True)

    if not ROOT.is_dir():
        raise RuntimeError("Repository non trovato")

    run([
        "docker",
        "image",
        "inspect",
        SOURCE_IMAGE,
    ], capture=True)

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP18_6_RUNTIME_GNN_{timestamp}"
    )

    output.mkdir(parents=True, exist_ok=False)

    dockerfile = output / "Dockerfile.gnn_cpu"

    dockerfile.write_text(
        DOCKERFILE.strip() + "\n",
        encoding="utf-8",
    )

    print("\n[1] Immagine sorgente:", SOURCE_IMAGE)
    print("[2] Immagine nuova:", NEW_IMAGE)
    print("[3] Dockerfile:", dockerfile)

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "18.6",
        "source_image": SOURCE_IMAGE,
        "new_image": NEW_IMAGE,
        "dockerfile": str(dockerfile),
        "build": "NOT_STARTED",
        "runtime_test": "NOT_STARTED",
        "training_performed": False,
        "original_images_modified": False,
        "gui_modified": False,
        "models_modified": False,
    }

    report_path = output / "RAPPORTO_STEP18_6.json"

    def save():
        report_path.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )

    save()

    try:
        print("\n[4] COSTRUZIONE IMMAGINE", flush=True)

        run([
            "docker", "build",
            "--pull=false",
            "-f", str(dockerfile),
            "-t", NEW_IMAGE,
            str(output),
        ])

        report["build"] = "PASS"
        save()

        print("\n[5] TEST LIBRERIE E GNN", flush=True)

        result = run([
            "docker", "run",
            "--rm",
            "--pull", "never",
            "--network", "none",
            "--read-only",
            "--tmpfs", "/tmp:rw,nosuid,size=128m",
            "--cap-drop", "ALL",
            "--security-opt", "no-new-privileges",
            "--entrypoint", "python3",
            NEW_IMAGE,
            "-B",
            "-c",
            TEST_CODE,
        ], capture=True)

        print(result)

        marker = "RUNTIME_TEST_RESULT="

        matches = [
            line[len(marker):]
            for line in result.splitlines()
            if line.startswith(marker)
        ]

        if len(matches) != 1:
            raise RuntimeError(
                "Risultato strutturato del test assente"
            )

        test = json.loads(matches[0])

        if test.get("status") != "PASS":
            raise RuntimeError("Test ML/GNN non superato")

        if not test.get("torch_cpu"):
            raise RuntimeError(
                "La distribuzione Torch non e CPU-only"
            )

        report["runtime_test"] = "PASS"
        report["test_details"] = test
        save()

    except Exception as error:
        report["error"] = str(error)
        save()
        raise

    print("\n" + "=" * 70)
    print("STEP 18.6 COMPLETATO")
    print("=" * 70)
    print("[OK] Immagine nuova:", NEW_IMAGE)
    print("[OK] Runtime PyTorch CPU")
    print("[OK] PyTorch Geometric")
    print("[OK] Test GCNConv e SAGEConv")
    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Nessun training")
    print("Report:", report_path)

    print(
        "\n[NOTA] I checkpoint reali verranno "
        "verificati nello STEP 19."
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("\n[ERRORE]", exc)
        sys.exit(1)