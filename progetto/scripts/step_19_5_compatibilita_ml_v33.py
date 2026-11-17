#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 19.5

Collaudo di compatibilita e riproducibilita
del modello supervised antileakage.

- Crea immagine Docker separata con sklearn 1.8.0
- Non modifica l'immagine GNN originale
- Carica il checkpoint ML esistente
- Esegue inferenza su tutta la tabella originale
- Confronta le predizioni con quelle salvate
- Controlla le probabilita di classificazione
- Verifica la salute della GUI sulla porta 8502
- Controlla SHA256 di modello e tabella
- Genera report JSON

Nessun training.
Nessuna modifica a GUI, motore o checkpoint.

ATTENZIONE: joblib puo eseguire codice durante il
caricamento. Utilizzare soltanto il checkpoint
originale del progetto, ritenuto attendibile.
"""

import hashlib
import json
import subprocess
import sys
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

RUNS = ROOT / "target_prediction_v33_vscode" / "runs"

BASE_IMAGE = "ai_primetta_v33_gnn_cpu:18.6"
NEW_IMAGE = "ai_primetta_v33_ml_compat:19.5"

MODEL = (
    ROOT
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

DATA = (
    ROOT
    / "outputs/prjna419326/"
    "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)

DOCKERFILE = f"""
FROM {BASE_IMAGE}

USER root

RUN python3 -m pip install \\
    --no-cache-dir \\
    --only-binary=:all: \\
    scikit-learn==1.8.0

RUN python3 -c "import sklearn; assert sklearn.__version__ == '1.8.0'"

ENV PYTHONDONTWRITEBYTECODE=1
ENV OMP_NUM_THREADS=2
ENV MKL_NUM_THREADS=2
"""

RUNNER = r'''
import json
import warnings

from collections import Counter
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn

from sklearn.exceptions import InconsistentVersionWarning


MODEL = Path(
    "/project/models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

DATA = Path(
    "/project/outputs/prjna419326/"
    "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)


def normalize_label(value):
    if pd.isna(value):
        return "<MISSING>"

    try:
        number = float(value)

        if np.isfinite(number):
            if number.is_integer():
                return str(int(number))

            return str(number)
    except (TypeError, ValueError):
        pass

    return str(value).strip()


def emit(data):
    print(
        "ML_RESULT="
        + json.dumps(
            data,
            ensure_ascii=False,
        ),
        flush=True
    )


def main():

    print("=" * 65)
    print("AI PRIMETTA V3.3 - ML ANTILEAKAGE")
    print("COMPATIBILITA SKLEARN 1.8.0")
    print("=" * 65)

    print("[1] VERSIONE SCIKIT-LEARN")

    version = sklearn.__version__

    print("Versione:", version)

    if version != "1.8.0":
        raise RuntimeError(
            "Versione scikit-learn non corretta"
        )

    print("[2] CARICAMENTO CHECKPOINT")

    # Qualsiasi InconsistentVersionWarning
    # viene considerato un errore.
    warnings.simplefilter(
        "error",
        InconsistentVersionWarning
    )

    with warnings.catch_warnings():
        warnings.simplefilter(
            "error",
            InconsistentVersionWarning
        )

        model = joblib.load(MODEL)

    print("[OK] Modello:", type(model).__name__)

    feature_names = getattr(
        model,
        "feature_names_in_",
        None
    )

    if feature_names is None:
        raise RuntimeError(
            "Schema feature non disponibile"
        )

    features = [
        str(name)
        for name in feature_names
    ]

    if len(features) != len(set(features)):
        raise RuntimeError(
            "Feature duplicate"
        )

    print("[OK] Feature:", len(features))

    if len(features) != 22:
        raise RuntimeError(
            "Attese 22 feature, "
            f"rilevate {len(features)}"
        )

    if not hasattr(model, "predict_proba"):
        raise RuntimeError(
            "Modello privo di predict_proba"
        )

    print("[3] LETTURA TABELLA ORIGINALE")

    header = pd.read_csv(
        DATA,
        sep="\t",
        nrows=0
    )

    missing = [
        name for name in features
        if name not in header.columns
    ]

    if missing:
        raise RuntimeError(
            "Feature mancanti: "
            + ", ".join(missing)
        )

    saved_label_column = (
        "antileakage_predicted_label"
        if "antileakage_predicted_label" in header.columns
        else None
    )

    print(
        "[OK] Predizioni storiche disponibili:",
        bool(saved_label_column)
    )

    print("[4] INFERENZA SU TUTTE LE RIGHE")

    total = 0
    distribution = Counter()

    compared_labels = 0
    matching_labels = 0
    max_probability_sum_error = 0.0

    classes = [
        normalize_label(value)
        for value in model.classes_
    ]

    for chunk in pd.read_csv(
        DATA,
        sep="\t",
        chunksize=512,
        low_memory=False
    ):

        X = chunk.loc[:, features].copy()

        for column in features:
            X[column] = pd.to_numeric(
                X[column],
                errors="coerce"
            )

        predictions = np.asarray(
            model.predict(X)
        )

        probabilities = np.asarray(
            model.predict_proba(X)
        )

        n = len(X)

        if predictions.shape[0] != n:
            raise RuntimeError(
                "Numero predizioni non coerente"
            )

        if probabilities.shape != (
            n,
            len(classes)
        ):
            raise RuntimeError(
                "Dimensioni probabilita inattese"
            )

        if not np.isfinite(probabilities).all():
            raise RuntimeError(
                "Probabilita NaN o infinite"
            )

        if (
            (probabilities < -1e-8).any()
            or (probabilities > 1 + 1e-8).any()
        ):
            raise RuntimeError(
                "Probabilita fuori [0,1]"
            )

        deviations = np.abs(
            probabilities.sum(axis=1) - 1.0
        )

        if len(deviations):
            max_probability_sum_error = max(
                max_probability_sum_error,
                float(deviations.max())
            )

        if max_probability_sum_error > 1e-6:
            raise RuntimeError(
                "Somma probabilita diversa da 1"
            )

        predicted_labels = [
            normalize_label(value)
            for value in predictions
        ]

        distribution.update(predicted_labels)
        total += n

        if saved_label_column:

            historical_labels = [
                normalize_label(value)
                for value in chunk[
                    saved_label_column
                ]
            ]

            matching_labels += sum(
                predicted == historical
                for predicted, historical
                in zip(
                    predicted_labels,
                    historical_labels
                )
            )

            compared_labels += n

    if total == 0:
        raise RuntimeError(
            "Tabella priva di osservazioni"
        )

    print("[OK] Righe:", total)

    print(
        "[OK] Distribuzione:",
        dict(distribution)
    )

    agreement = None

    if compared_labels:
        agreement = (
            matching_labels / compared_labels
        )

        print(
            "[OK] Concordanza con predizioni storiche:",
            f"{agreement:.6%}"
        )

    print("[5] RISULTATO")

    result = {
        "status": "PASS",
        "sklearn_version": version,
        "model_type": type(model).__name__,
        "feature_count": len(features),
        "features": features,
        "classes": classes,
        "rows_tested": total,
        "prediction_distribution": dict(
            distribution
        ),
        "stored_label_column": saved_label_column,
        "labels_compared": compared_labels,
        "labels_matching": matching_labels,
        "historical_label_agreement": agreement,
        "max_probability_sum_error": (
            max_probability_sum_error
        ),
        "inconsistent_version_warning": False,
        "independent_validation": False,
        "scientific_performance_validated": False,
        "training_performed": False,
    }

    emit(result)


if __name__ == "__main__":
    main()
'''


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(4 * 1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def execute(args, timeout=None, capture=False):
    process = subprocess.run(
        args,
        text=True,
        capture_output=capture,
        timeout=timeout,
        check=False
    )

    if process.returncode != 0:
        details = (
            (process.stdout or "")
            + "\n"
            + (process.stderr or "")
        ) if capture else ""

        raise RuntimeError(
            "Comando fallito "
            f"(exit {process.returncode}): "
            + " ".join(args[:3])
            + "\n"
            + details[-2000:]
        )

    return process


def gui_health():
    url = "http://127.0.0.1:8502/_stcore/health"

    try:
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            url,
            timeout=6
        ) as response:
            body = response.read(100).decode(
                "utf-8",
                errors="replace"
            )

            return {
                "url": url,
                "http_status": response.status,
                "body": body,
                "ok": (
                    response.status == 200
                    and body.strip() == "ok"
                )
            }

    except Exception as exc:
        return {
            "url": url,
            "ok": False,
            "error": str(exc)
        }


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 19.5")
    print("COMPATIBILITA ML E VERIFICA COMPLETA")
    print("=" * 70, flush=True)

    for path in (MODEL, DATA):
        if not path.is_file():
            raise FileNotFoundError(path)

    RUNS.mkdir(
        parents=True,
        exist_ok=True
    )

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP19_5_COMPATIBILITA_ML_{stamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    dockerfile = output / "Dockerfile.ml_compat"

    dockerfile.write_text(
        DOCKERFILE.strip() + "\n",
        encoding="utf-8"
    )

    report_path = output / "RAPPORTO_STEP19_5.json"

    before = {
        "model": sha256(MODEL),
        "data": sha256(DATA)
    }

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "19.5",
        "status": "NOT_STARTED",
        "source_image": BASE_IMAGE,
        "compatibility_image": NEW_IMAGE,
        "original_sklearn_version": "1.8.0",
        "dockerfile": str(dockerfile),
        "sha256_before": before,
        "training_performed": False,
        "existing_components_modified": False,
    }

    def save():
        report_path.write_text(
            json.dumps(
                report,
                indent=2,
                ensure_ascii=False
            ),
            encoding="utf-8"
        )

    save()

    try:

        print(
            "\n[1] VERIFICA IMMAGINE SORGENTE",
            flush=True
        )

        execute([
            "docker", "image", "inspect",
            BASE_IMAGE
        ], capture=True)

        print("[OK]", BASE_IMAGE)

        print(
            "\n[2] BUILD AMBIENTE SKLEARN 1.8.0",
            flush=True
        )

        execute([
            "docker", "build",
            "--pull=false",
            "-f", str(dockerfile),
            "-t", NEW_IMAGE,
            str(output)
        ])

        report["docker_build"] = "PASS"
        save()

        print(
            "\n[3] INFERENZA SU TUTTA LA TABELLA",
            flush=True
        )

        process = execute([
            "docker", "run",
            "--rm",
            "--pull", "never",
            "--network", "none",
            "--read-only",
            "--tmpfs", "/tmp:rw,nosuid,nodev,size=128m",
            "--cap-drop", "ALL",
            "--security-opt", "no-new-privileges",
            "--pids-limit", "128",
            "--memory", "3g",
            "--cpus", "2",
            "--user", "65534:65534",
            "-e", "HOME=/tmp",
            "-e", "OMP_NUM_THREADS=2",
            "-v", f"{ROOT}:/project:ro",
            "-w", "/project",
            "--entrypoint", "python3",
            NEW_IMAGE,
            "-B", "-c", RUNNER
        ], timeout=300, capture=True)

        print(process.stdout, flush=True)

        if process.stderr:
            print(
                process.stderr[-2500:],
                file=sys.stderr
            )

        marker = "ML_RESULT="

        matches = [
            line[len(marker):]
            for line in process.stdout.splitlines()
            if line.startswith(marker)
        ]

        if len(matches) != 1:
            raise RuntimeError(
                "Risultato strutturato non trovato"
            )

        result = json.loads(matches[0])

        if result["status"] != "PASS":
            raise RuntimeError(
                "Test ML non superato"
            )

        report["ml_inference"] = result

        if (
            result["historical_label_agreement"]
            is not None
            and result["historical_label_agreement"] < 1.0
        ):
            report["status"] = (
                "INFERENCE_PASS_HISTORICAL_MISMATCH"
            )
        else:
            report["status"] = (
                "COMPATIBLE_ML_INFERENCE_PASS"
            )

        print(
            "\n[4] VERIFICA GUI V3.3",
            flush=True
        )

        report["gui_health"] = gui_health()

        print(
            "HTTP GUI:",
            report["gui_health"].get(
                "http_status", "NON RAGGIUNGIBILE"
            )
        )

    except Exception as exc:

        report["status"] = "FAILED"
        report["error"] = str(exc)

        print(
            "\n[ERRORE]",
            exc,
            flush=True
        )

    finally:

        after = {
            "model": sha256(MODEL),
            "data": sha256(DATA)
        }

        report["sha256_after"] = after
        report["original_files_unchanged"] = (
            before == after
        )

        if before != after:
            report["status"] = (
                "FAILED_ORIGINAL_FILES_CHANGED"
            )

        save()

    print()
    print("=" * 70)
    print("STEP 19.5 - RISULTATO")
    print("=" * 70)

    print("Stato:", report["status"])

    print(
        "Originali invariati:",
        report["original_files_unchanged"]
    )

    print("Report:", report_path)

    print("[OK] Nessun training")
    print("[OK] GUI V3.2 non modificata")
    print("[OK] Motore V3.3 non modificato")

    if report["status"] != "COMPATIBLE_ML_INFERENCE_PASS":
        sys.exit(1)

    if not report.get(
        "gui_health", {}
    ).get("ok", False):
        print(
            "[ATTENZIONE] Inferenza ML verificata, "
            "GUI da controllare."
        )
        sys.exit(1)

    print(
        "\n[OK] Inferenza ML compatibile superata"
    )

    print(
        "[OK] Servizio GUI raggiungibile"
    )

    print(
        "[NOTA] Questo test non equivale "
        "a validazione scientifica indipendente."
    )


if __name__ == "__main__":
    main()