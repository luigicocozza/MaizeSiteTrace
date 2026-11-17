#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 19.4

Collaudo del modello supervised ML antileakage reale.

- Verifica il checkpoint
- Usa Docker ML/GNN gia collaudato
- Identifica le feature attese dal modello
- Cerca matrici dati reali compatibili
- Esegue predict e predict_proba quando possibile
- Verifica l'integrita del checkpoint
- Salva il report JSON

Nessun training, nessuna modifica alla GUI o ai dati.

Sicurezza: joblib usa la deserializzazione Python.
Caricare esclusivamente checkpoint locali attendibili.
"""

import hashlib
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

IMAGE = "ai_primetta_v33_gnn_cpu:18.6"

MODEL = (
    ROOT
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)


RUNNER = r'''
import json
import os
import sys

from collections import Counter
from pathlib import Path

import joblib
import numpy as np
import pandas as pd


ROOT = Path("/project")

MODEL = (
    ROOT
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)


def emit(result):
    print(
        "INFERENCE_RESULT=" +
        json.dumps(result, ensure_ascii=False),
        flush=True
    )


def find_feature_names(model):

    names = getattr(
        model, "feature_names_in_", None
    )

    if names is not None:
        return [str(x) for x in names]

    for step in getattr(
        model, "named_steps", {}
    ).values():

        names = getattr(
            step, "feature_names_in_", None
        )

        if names is not None:
            return [str(x) for x in names]

    return []


def candidate_tables():

    skip = {
        "references", "runs", "delivery",
        "final_package", "v2_package",
        ".venv", "venv", ".git",
        "__pycache__", "node_modules",
    }

    for base in (
        ROOT / "data",
        ROOT / "outputs",
    ):

        if not base.is_dir():
            continue

        for current, dirs, files in os.walk(base):

            dirs[:] = [
                x for x in dirs
                if x not in skip
            ]

            for name in files:

                path = Path(current) / name

                if path.suffix.lower() not in (
                    ".csv", ".tsv"
                ):
                    continue

                try:
                    size = path.stat().st_size
                except OSError:
                    continue

                if 0 < size <= 150 * 1024 * 1024:
                    yield path


def read_table(path, columns, nrows=20):

    separator = (
        "\t" if path.suffix.lower() == ".tsv"
        else ","
    )

    return pd.read_csv(
        path,
        sep=separator,
        usecols=columns,
        nrows=nrows,
        low_memory=False
    )


def main():

    print("[1] CARICAMENTO MODELLO REALE", flush=True)

    model = joblib.load(MODEL)

    print(
        "[OK] Tipo modello:",
        type(model).__name__,
        flush=True
    )

    expected = find_feature_names(model)

    if not expected:

        emit({
            "status": "BLOCKED",
            "reason": "FEATURE_SCHEMA_NOT_AVAILABLE",
            "model_type": type(model).__name__,
            "model_loaded": True,
            "inference_tested": False,
            "training_performed": False,
        })
        return

    if len(expected) != len(set(expected)):
        raise RuntimeError(
            "Nomi delle feature duplicati"
        )

    print(
        "[OK] Feature richieste:",
        len(expected),
        flush=True
    )

    print("[2] RICERCA MATRICE DATI COMPATIBILE")

    candidates = []
    inspected = 0

    for path in candidate_tables():

        inspected += 1

        separator = (
            "\t" if path.suffix.lower() == ".tsv"
            else ","
        )

        try:
            header = pd.read_csv(
                path, sep=separator,
                nrows=0
            )

            available = set(header.columns)

        except Exception:
            continue

        common = sum(
            feature in available
            for feature in expected
        )

        if common == len(expected):
            candidates.append(path)

    print("[OK] Tabelle esaminate:", inspected)
    print("[OK] Tabelle compatibili:", len(candidates))

    if not candidates:

        emit({
            "status": "BLOCKED",
            "reason": "NO_COMPATIBLE_REAL_FEATURE_TABLE",
            "model_loaded": True,
            "expected_features": expected,
            "tables_inspected": inspected,
            "inference_tested": False,
            "training_performed": False,
        })
        return

    # Priorita ai file del workflow antileakage.
    candidates.sort(
        key=lambda path: (
            "antileakage" not in str(path).lower(),
            "srna" not in str(path).lower(),
            str(path),
        )
    )

    print("[3] INFERENZA SUI DATI REALI")

    failures = []

    for path in candidates:

        try:
            X = read_table(path, expected)

            # Ordine delle feature: esattamente
            # quello registrato dal modello.
            X = X.loc[:, expected]

            for column in expected:
                X[column] = pd.to_numeric(
                    X[column],
                    errors="coerce"
                )

            array = X.to_numpy(dtype=float)

            if not len(array):
                raise ValueError(
                    "Nessuna riga disponibile"
                )

            valid = (
                ~np.isinf(array).any(axis=1)
                & np.isfinite(array).any(axis=1)
            )

            X = X.loc[valid].copy()

            if X.empty:
                raise ValueError(
                    "Nessuna riga valida"
                )

            predictions = model.predict(X)

            if len(predictions) != len(X):
                raise RuntimeError(
                    "Numero predizioni incoerente"
                )

            if hasattr(model, "predict_proba"):

                probabilities = np.asarray(
                    model.predict_proba(X)
                )

                if not np.isfinite(
                    probabilities
                ).all():
                    raise RuntimeError(
                        "Probabilita non finite"
                    )

                if (
                    (probabilities < -1e-8).any()
                    or (probabilities > 1 + 1e-8).any()
                ):
                    raise RuntimeError(
                        "Probabilita fuori intervallo"
                    )

                probability_shape = list(
                    probabilities.shape
                )

            else:
                probability_shape = None

            distribution = Counter(
                str(value)
                for value in predictions
            )

            result = {
                "status": "PASS",
                "model_type": type(model).__name__,
                "model_loaded": True,
                "feature_count": len(expected),
                "feature_names": expected,
                "source_table": str(
                    path.relative_to(ROOT)
                ),
                "rows_tested": len(X),
                "predictions_distribution": dict(
                    distribution
                ),
                "predict_proba_shape": probability_shape,
                "inference_tested": True,
                "training_performed": False,
                "independent_validation": False,
                "scientific_performance_validated": False,
            }

            print(
                "[OK] Tabella:",
                result["source_table"]
            )

            print(
                "[OK] Righe predette:",
                len(X)
            )

            print(
                "[OK] Classi predette:",
                result["predictions_distribution"]
            )

            emit(result)
            return

        except Exception as exc:
            failures.append({
                "file": str(path.relative_to(ROOT)),
                "error": (
                    type(exc).__name__
                    + ": " + str(exc)
                )[:300]
            })

    emit({
        "status": "BLOCKED",
        "reason": "COMPATIBLE_TABLES_FAILED_INFERENCE",
        "model_loaded": True,
        "feature_count": len(expected),
        "candidate_tables": len(candidates),
        "failures": failures[:10],
        "inference_tested": False,
        "training_performed": False,
    })


if __name__ == "__main__":
    main()
'''


def sha256(path):

    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(4 * 1024 * 1024),
            b""
        ):
            digest.update(chunk)

    return digest.hexdigest()


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 19.4")
    print("COLLAUDO SUPERVISED ML ANTILEAKAGE")
    print("=" * 70)

    if not MODEL.is_file():
        raise FileNotFoundError(MODEL)

    before = sha256(MODEL)

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP19_4_ML_ANTILEAKAGE_{stamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report_path = (
        output / "RAPPORTO_STEP19_4.json"
    )

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "19.4",
        "model": str(MODEL),
        "runtime": IMAGE,
        "status": "NOT_TESTED",
        "sha256_before": before,
        "training_performed": False,
    }

    command = [
        "docker", "run",
        "--rm",
        "--pull", "never",
        "--network", "none",
        "--read-only",
        "--tmpfs", "/tmp:rw,nosuid,nodev,size=128m",
        "--cap-drop", "ALL",
        "--security-opt", "no-new-privileges",
        "--pids-limit", "128",
        "--memory", "2g",
        "--cpus", "2",
        "--user", "65534:65534",
        "-e", "HOME=/tmp",
        "-e", "OMP_NUM_THREADS=2",
        "-v", f"{ROOT}:/project:ro",
        "-w", "/project",
        "--entrypoint", "python3",
        IMAGE,
        "-B", "-c", RUNNER
    ]

    print("\n[1] AVVIO TEST DOCKER")

    try:

        proc = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=180,
            check=False
        )

        if proc.stdout:
            print(proc.stdout)

        if proc.stderr:
            print(proc.stderr[-2500:])

        if proc.returncode != 0:
            raise RuntimeError(
                "Docker exit code: "
                + str(proc.returncode)
            )

        marker = "INFERENCE_RESULT="

        found = [
            line[len(marker):]
            for line in proc.stdout.splitlines()
            if line.startswith(marker)
        ]

        if len(found) != 1:
            raise RuntimeError(
                "Risultato strutturato assente"
            )

        result = json.loads(found[0])

        report["inference"] = result

        report["status"] = (
            "TECHNICAL_INFERENCE_PASS"
            if result["status"] == "PASS"
            else "INFERENCE_BLOCKED"
        )

    except Exception as exc:

        report["status"] = "FAILED"
        report["error"] = str(exc)

    finally:

        after = sha256(MODEL)

        report["sha256_after"] = after

        report["model_unchanged"] = (
            before == after
        )

        if before != after:
            report["status"] = (
                "FAILED_MODEL_HASH_CHANGED"
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
    print("STEP 19.4 - RISULTATO")
    print("=" * 70)

    print("Stato:", report["status"])
    print(
        "Checkpoint invariato:",
        report["model_unchanged"]
    )

    print("Report:", report_path)
    print("[OK] Nessun training")
    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")

    if report["status"] != "TECHNICAL_INFERENCE_PASS":
        sys.exit(1)


if __name__ == "__main__":
    main()