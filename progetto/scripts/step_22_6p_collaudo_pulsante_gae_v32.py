#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.2 - STEP 22.6P

--before:
    Salva una baseline di integrità prima
    del training GAE dalla GUI.

--after:
    Controlla che il nuovo training abbia
    generato i file nelle cartelle isolate
    e che gli originali siano invariati.

Non avvia training, non modifica GUI o Docker.
"""

import argparse
import hashlib
import json
import subprocess
import sys
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT / "target_prediction_v33_vscode/runs"
)

SAFE_OUTPUTS = (
    ROOT / "outputs/gae_gui_safe_runs"
)

PATCH = (
    ROOT
    / "docker/runtime_v32_patch/"
    "25_train_graph_autoencoder_safe.py"
)

CONTAINER = "maize_ncrna_gnn"

IMAGE = (
    "ai_primetta_v32_full_ai_sklearn18:22.6g"
)

MAX_PROTECTED_BYTES = 500 * 1024 * 1024


# ============================================================
# 1. COMANDI
# ============================================================

def command(args, timeout=40):

    p = subprocess.run(
        [str(x) for x in args],
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout
    )

    if p.returncode:

        raise RuntimeError(
            "Comando fallito: "
            + " ".join(map(str, args[:5]))
            + "\n"
            + p.stderr[-1200:]
        )

    return p.stdout.strip()


# ============================================================
# 2. CHECKSUM SHA256
# ============================================================

def sha256(path):

    result = hashlib.sha256()

    with path.open("rb") as stream:

        for chunk in iter(
            lambda: stream.read(1024 * 1024),
            b""
        ):

            result.update(chunk)

    return result.hexdigest()


# ============================================================
# 3. FILE SCIENTIFICI DA PROTEGGERE
# ============================================================

def protected_paths():

    paths = [

        ROOT / "25_train_graph_autoencoder.py",

        ROOT / "api/main.py",

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
        / "models/"
        "srna_supervised_ml_v1_antileakage/"
        "PRJNA419326_srna_supervised_"
        "antileakage_best_model.joblib",

        ROOT
        / "data/processed/embeddings/"
        "gae_node_embeddings.csv",

        ROOT
        / "data/processed/embeddings/"
        "gae_ncrna_embeddings.csv",
    ]

    # Protezione di tutti i file delle due
    # cartelle del GAE originale.

    for base in (

        ROOT / "models/graph_autoencoder",

        ROOT / "outputs/gae"
    ):

        if not base.is_dir():

            raise RuntimeError(
                "Cartella originale assente: "
                + str(base)
            )

        paths.extend(
            sorted(base.rglob("*"))
        )

    return paths


# ============================================================
# 4. BASELINE FILE ORIGINALI
# ============================================================

def snapshot():

    files = {}

    total = 0

    for path in protected_paths():

        if path.is_symlink():

            raise RuntimeError(
                "Link simbolico inatteso: "
                + str(path)
            )

        if path.is_dir():
            continue

        if not path.is_file():

            raise RuntimeError(
                "File da proteggere assente: "
                + str(path)
            )

        size = path.stat().st_size

        total += size

        if total > MAX_PROTECTED_BYTES:

            raise RuntimeError(
                "File protetti oltre 500 MB: "
                "richiesta verifica manuale"
            )

        files[
            str(path.relative_to(ROOT))
        ] = {

            "sha256": sha256(path),

            "bytes": size
        }

    return files


# ============================================================
# 5. VERIFICA DOCKER
# ============================================================

def docker_state():

    raw = command([
        "docker",
        "inspect",
        CONTAINER
    ])

    doc = json.loads(raw)[0]

    mounts = {

        v["Destination"]: v

        for v in doc.get("Mounts", [])
    }

    if not doc.get(
        "State", {}
    ).get("Running"):

        raise RuntimeError(
            "Container non in esecuzione"
        )

    if doc["Config"]["Image"] != IMAGE:

        raise RuntimeError(
            "Immagine Docker inattesa"
        )

    # Verifica che il trainer sia montato
    # in sola lettura e gli output siano persistenti.

    required = {

        "/app/25_train_graph_autoencoder.py":
            (PATCH, False),

        "/app/outputs":
            (ROOT / "outputs", True),

        "/app/models":
            (ROOT / "models", True),

        "/app/data/processed":
            (ROOT / "data/processed", True),

        "/app/data/user_data":
            (ROOT / "data/user_data", True),
    }

    for target, (host, writable) in required.items():

        item = mounts.get(target)

        if (
            not item

            or Path(
                item["Source"]
            ).resolve() != host.resolve()

            or bool(
                item.get("RW")
            ) != writable
        ):

            raise RuntimeError(
                "Montaggio Docker inatteso: "
                + target
            )

    if not PATCH.is_file():

        raise RuntimeError(
            "Trainer isolato assente"
        )

    # Confronta il file realmente visto dal
    # container con la patch sul disco Ubuntu.

    inside_hash = command([
        "docker",
        "exec",
        CONTAINER,
        "python3",
        "-B",
        "-c",
        (
            "import pathlib,hashlib;"
            "print(hashlib.sha256("
            "pathlib.Path("
            "'/app/25_train_graph_autoencoder.py'"
            ").read_bytes()).hexdigest())"
        )
    ])

    if inside_hash != sha256(PATCH):

        raise RuntimeError(
            "Il trainer del container non "
            "coincide con quello isolato"
        )

    # Verifica GUI e API.

    for port, suffix in (

        (8503, "/_stcore/health"),

        (8000, "/health")
    ):

        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            f"http://127.0.0.1:{port}{suffix}",
            timeout=7
        ) as response:

            if response.status != 200:

                raise RuntimeError(
                    f"Servizio {port} non sano"
                )

    return {

        "image": IMAGE,

        "trainer_sha256": inside_hash,

        "trainer_mounted_read_only": True,

        "uploads_persistent": True
    }


# ============================================================
# 6. TRAINING GAE GIÀ PRESENTI
# ============================================================

def runs_present():

    if not SAFE_OUTPUTS.is_dir():
        return []

    return sorted(

        p.name

        for p in SAFE_OUTPUTS.iterdir()

        if p.is_dir()
        and not p.is_symlink()
    )


# ============================================================
# 7. FASE PRECEDENTE AL TRAINING
# ============================================================

def baseline():

    state = docker_state()

    hashes = snapshot()

    initial = runs_present()

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP22_6P_GAE_BUTTON_" + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    content = {

        "step": "22.6P",

        "created_utc": stamp,

        "docker": state,

        "original_sha256": hashes,

        "runs_before": initial,

        "training_started_by_this_script": False
    }

    path = (
        folder / "BASELINE_STEP22_6P.json"
    )

    path.write_text(

        json.dumps(
            content,
            indent=2,
            ensure_ascii=False
        ),

        encoding="utf-8"
    )

    print(
        "[OK] Trainer isolato montato, "
        "immagine e volumi conformi"
    )

    print(
        "[OK] File originali controllati:",
        len(hashes)
    )

    print(
        "[OK] Esecuzioni GAE già presenti:",
        len(initial)
    )

    print(
        "[BASELINE]",
        path
    )

    print(
        "[PROSSIMA AZIONE] Avviare un training "
        "dalla GUI 8503 e poi eseguire --after"
    )


# ============================================================
# 8. COLLAUDO DEL CHECKPOINT NEL CONTAINER
# ============================================================

CHECKPOINT_PROBE = r'''

import json
import sys
import torch

path = sys.argv[1]

# Caricamento esclusivo del checkpoint locale
# appena prodotto dal progetto attendibile.

checkpoint = torch.load(
    path,
    map_location="cpu",
    weights_only=False
)

assert isinstance(
    checkpoint,
    dict
)

state = checkpoint.get(
    "model_state_dict"
)

assert isinstance(
    state,
    dict
) and state

for key, value in state.items():

    if isinstance(
        value,
        torch.Tensor
    ):

        assert bool(
            torch.isfinite(value).all()
        ), key

print(
    "PRIMETTA_CHECKPOINT="
    + json.dumps({

        "state_dict_entries": len(state),

        "in_channels": checkpoint.get(
            "in_channels"
        ),

        "keys": list(
            checkpoint.keys()
        )[:16]
    })
)

'''


def checkpoint_check(checkpoint):

    relative = checkpoint.relative_to(
        ROOT / "outputs"
    )

    inside = (
        "/app/outputs/"
        + relative.as_posix()
    )

    stdout = command(
        [
            "docker",
            "exec",
            CONTAINER,
            "python3",
            "-B",
            "-c",
            CHECKPOINT_PROBE,
            inside
        ],
        timeout=75
    )

    line = next(
        (
            s for s in stdout.splitlines()

            if s.startswith(
                "PRIMETTA_CHECKPOINT="
            )
        ),
        None
    )

    if not line:

        raise RuntimeError(
            "Risposta checkpoint non valida: "
            + stdout[-800:]
        )

    return json.loads(
        line.split("=", 1)[1]
    )


# ============================================================
# 9. SELEZIONE BASELINE
# ============================================================

def select_baseline():

    matches = sorted(
        RUNS.glob(
            "STEP22_6P_GAE_BUTTON_*/"
            "BASELINE_STEP22_6P.json"
        )
    )

    matches = [

        p for p in matches

        if not (
            p.parent / "VERIFICA_STEP22_6P.json"
        ).exists()
    ]

    if not matches:

        raise RuntimeError(
            "Nessuna baseline disponibile. "
            "Eseguire prima --before"
        )

    return matches[-1]


# ============================================================
# 10. VERIFICA DOPO IL TRAINING
# ============================================================

def verify():

    path = select_baseline()

    base = json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )

    state = docker_state()

    after_hashes = snapshot()

    originals_ok = (
        base["original_sha256"]
        == after_hashes
    )

    differences = sorted(

        set(
            base["original_sha256"]
        )

        ^ set(after_hashes)
    )

    for key in (
        set(base["original_sha256"])
        & set(after_hashes)
    ):

        if (
            base["original_sha256"][key]
            != after_hashes[key]
        ):

            differences.append(key)

    differences.sort()

    new_runs = sorted(

        set(runs_present())

        - set(base["runs_before"])
    )

    print(
        "[ORIGINALI]",
        "PASS" if originals_ok else "FAIL",
        "| Differenze:",
        len(differences)
    )

    for name in differences:

        print(
            "  [MODIFICATO]",
            name
        )

    print(
        "[NUOVI TRAINING]",
        len(new_runs)
    )

    details = []

    # --------------------------------------------------------
    # 11. COLLAUDO FILE GENERATI
    # --------------------------------------------------------

    for name in new_runs:

        folder = SAFE_OUTPUTS / name

        checkpoint = (
            folder / "models/gae_model.pt"
        )

        history = (
            folder
            / "results/gae_training_history.csv"
        )

        node = (
            folder
            / "embeddings/gae_node_embeddings.csv"
        )

        ncrna = (
            folder
            / "embeddings/gae_ncrna_embeddings.csv"
        )

        required = (
            checkpoint,
            history,
            node,
            ncrna
        )

        present = all(

            p.is_file()

            and p.stat().st_size > 0

            for p in required
        )

        result = {

            "run_id": name,

            "required_files_present": present,

            "files": {

                str(p.relative_to(folder)):
                    p.is_file()

                for p in required
            },

            "checkpoint": None,

            "error": None
        }

        if present:

            try:

                result["checkpoint"] = (
                    checkpoint_check(
                        checkpoint
                    )
                )

            except Exception as exc:

                result["error"] = repr(exc)

        else:

            result["error"] = (
                "Checkpoint o CSV "
                "mancanti oppure vuoti"
            )

        result["passed"] = (
            present
            and result["error"] is None
        )

        details.append(result)

        print(
            "[RUN]",
            name,
            (
                "PASS"
                if result["passed"]
                else "FAIL"
            )
        )

        if result["error"]:

            print(
                "[DETTAGLIO]",
                result["error"]
            )

    # --------------------------------------------------------
    # 12. ESITO
    # --------------------------------------------------------

    passed = (

        originals_ok

        and bool(details)

        and all(
            x["passed"]
            for x in details
        )
    )

    report = {

        "step": "22.6P",

        "passed": passed,

        "baseline": str(path),

        "docker": state,

        "original_files_unchanged":
            originals_ok,

        "original_differences":
            differences,

        "new_safe_runs":
            details,

        "training_executed_by_verifier":
            False,

        "ui_original_changed":
            False
    }

    output = (
        path.parent / "VERIFICA_STEP22_6P.json"
    )

    output.write_text(

        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),

        encoding="utf-8"
    )

    print(
        "[REPORT]",
        output
    )

    print(
        "[ESITO]",
        "PASS"
        if passed
        else "DA RISOLVERE"
    )

    if not passed:
        sys.exit(2)


# ============================================================
# 13. PROGRAMMA PRINCIPALE
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

    group.add_argument(
        "--before",
        action="store_true"
    )

    group.add_argument(
        "--after",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 70)

    print(
        "AI PRIMETTA V3.2 | STEP 22.6P "
        "| COLLAUDO PULSANTE GAE"
    )

    print("=" * 70)

    if args.before:

        baseline()

    else:

        verify()


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