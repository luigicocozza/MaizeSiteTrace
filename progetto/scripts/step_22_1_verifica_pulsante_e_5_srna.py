#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 22.1

Verifica:
1. Codice eseguito dal pulsante V3.3
2. Chiamate subprocess
3. Inferenza ML sulle 5 small RNA reali

Nessuna modifica a GUI, checkpoint o dati.
Nessun training.
"""

import ast
import hashlib
import json
import re
import subprocess
import sys

from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

BASE = ROOT / "target_prediction_v33_vscode"
GUI = BASE / "app_v33.py"
RUNS = BASE / "runs"

MODEL = (
    ROOT
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

TABLE = (
    ROOT
    / "outputs/prjna419326/"
    "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)

IMAGE = "ai_primetta_v33_ml_compat:19.5"


# ==========================================================
# CODICE INFERENZA DA ESEGUIRE IN DOCKER
# ==========================================================

TEST_IN_CONTAINER = r'''
import json
import re
import warnings

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn

from sklearn.exceptions import InconsistentVersionWarning


root = Path("/project")

model_path = (
    root
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

table_path = (
    root
    / "outputs/prjna419326/"
    "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)

ids = [
    "sRNA_00000002",
    "sRNA_00000040",
    "sRNA_00000007",
    "sRNA_00000435",
    "sRNA_00000037",
]

result = {
    "status": "NOT_VERIFIED",
    "sklearn": sklearn.__version__,
    "records": [],
    "input_source": "historical_supervised_table",
    "gui_invoked_model": False,
    "independent_validation": False,
    "training_performed": False,
}

try:

    if sklearn.__version__ != "1.8.0":
        raise RuntimeError(
            "Runtime diverso da scikit-learn 1.8.0"
        )

    with warnings.catch_warnings():

        warnings.simplefilter(
            "error",
            InconsistentVersionWarning
        )

        # Caricare soltanto il checkpoint locale
        # attendibile del progetto.
        model = joblib.load(model_path)

    features = list(
        getattr(
            model,
            "feature_names_in_",
            []
        )
    )

    if len(features) != 22:
        raise RuntimeError(
            "Schema delle 22 feature assente o diverso"
        )

    frame = pd.read_csv(
        table_path,
        sep="\t"
    )

    if "feature_id" not in frame.columns:
        raise RuntimeError(
            "Colonna feature_id assente"
        )

    missing_features = [
        feature
        for feature in features
        if feature not in frame.columns
    ]

    if missing_features:
        raise RuntimeError(
            "Feature mancanti: "
            + str(missing_features)
        )

    matched = []

    for target in ids:

        pattern = re.compile(
            r"(?<![A-Za-z0-9])"
            + re.escape(target)
            + r"(?![0-9])"
        )

        subset = frame.loc[
            frame["feature_id"]
            .astype(str)
            .map(
                lambda value: bool(
                    pattern.search(value)
                )
            )
        ]

        if len(subset) != 1:
            raise RuntimeError(
                f"{target}: {len(subset)} record trovati, "
                "richiesto 1"
            )

        matched.append(
            subset.iloc[0]
        )

    df = pd.DataFrame(matched)

    X = df.loc[:, features].copy()

    for feature in features:
        X[feature] = pd.to_numeric(
            X[feature],
            errors="coerce"
        )

    predictions = model.predict(X)

    probabilities = np.asarray(
        model.predict_proba(X)
    )

    if probabilities.shape != (
        5,
        len(model.classes_)
    ):
        raise RuntimeError(
            "Dimensioni predict_proba inattese"
        )

    if not np.isfinite(
        probabilities
    ).all():
        raise RuntimeError(
            "Probabilita non finite"
        )

    if not np.allclose(
        probabilities.sum(axis=1),
        1.0,
        atol=1e-6
    ):
        raise RuntimeError(
            "Probabilita non normalizzate"
        )

    classes = [
        str(value)
        for value in model.classes_
    ]

    positive_idx = (
        classes.index("1")
        if "1" in classes
        else None
    )

    for index, target in enumerate(ids):

        label = str(
            predictions[index]
        )

        try:
            label = str(
                int(float(label))
            )
        except (ValueError, TypeError):
            pass

        historical = df.iloc[index].get(
            "antileakage_predicted_label"
        )

        if pd.isna(historical):
            historical_label = None

        else:
            try:
                historical_label = str(
                    int(float(historical))
                )
            except (ValueError, TypeError):
                historical_label = str(
                    historical
                )

        result["records"].append({
            "rna": target,
            "feature_id": str(
                df.iloc[index]["feature_id"]
            ),
            "predicted_class": label,
            "probability_class_1": (
                float(
                    probabilities[
                        index,
                        positive_idx
                    ]
                )
                if positive_idx is not None
                else None
            ),
            "historical_label": historical_label,
            "matches_historical": (
                label == historical_label
                if historical_label is not None
                else None
            ),
        })

    result["status"] = (
        "PASS_5_REAL_SRNA_HISTORICAL_FEATURES"
    )

except Exception as exc:

    result["error"] = (
        type(exc).__name__
        + ": "
        + str(exc)
    )

print(
    "ML_TEST_JSON="
    + json.dumps(
        result,
        ensure_ascii=False
    ),
    flush=True
)
'''


# ==========================================================
# FUNZIONI DI SUPPORTO
# ==========================================================

def digest(path):

    if not path.is_file():
        return None

    h = hashlib.sha256()

    with path.open("rb") as source:

        for chunk in iter(
            lambda: source.read(
                4 * 1024 * 1024
            ),
            b""
        ):
            h.update(chunk)

    return h.hexdigest()


def call(args, timeout=180):

    try:

        process = subprocess.run(
            args,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False
        )

        return {
            "ok": process.returncode == 0,
            "returncode": process.returncode,
            "stdout": process.stdout,
            "stderr": process.stderr[-2000:]
        }

    except Exception as exc:

        return {
            "ok": False,
            "error": str(exc)
        }


def fun_name(node):

    try:
        return ast.unparse(node)

    except Exception:
        return "?"


# ==========================================================
# ANALISI DEL PULSANTE DELLA GUI V3.3
# ==========================================================

def analyse_button():

    content = GUI.read_text(
        encoding="utf-8"
    )

    tree = ast.parse(content)

    lines = content.splitlines()

    result = {
        "button_line": None,
        "button_calls": [],
        "subprocess_calls": [],
        "direct_model_loading_calls": [],
    }

    for node in ast.walk(tree):

        if isinstance(node, ast.If):

            test_calls = [
                item
                for item in ast.walk(node.test)
                if isinstance(
                    item,
                    ast.Call
                )
            ]

            is_button = any(
                "st.button" in fun_name(call.func)
                or "st.form_submit_button"
                in fun_name(call.func)
                for call in test_calls
            )

            if not is_button:
                continue

            text = "\n".join(
                lines[
                    node.lineno - 1:
                    min(
                        node.lineno + 8,
                        len(lines)
                    )
                ]
            )

            if (
                "Calcola geni target candidati"
                not in text
            ):
                continue

            result["button_line"] = (
                node.lineno
            )

            calls = []

            for statement in node.body:

                for item in ast.walk(
                    statement
                ):

                    if isinstance(
                        item,
                        ast.Call
                    ):

                        calls.append({
                            "line": item.lineno,
                            "call": fun_name(
                                item.func
                            )[:130]
                        })

            result["button_calls"] = sorted(
                calls,
                key=lambda item: item["line"]
            )

        if isinstance(node, ast.Call):

            name = fun_name(
                node.func
            )

            if name in (
                "subprocess.run",
                "subprocess.Popen",
                "os.system"
            ):

                start = max(
                    0,
                    node.lineno - 11
                )

                end = min(
                    len(lines),
                    (node.end_lineno or node.lineno) + 9
                )

                excerpt = "\n".join(
                    f"{i+1}: {lines[i]}"
                    for i in range(
                        start,
                        end
                    )
                )

                result["subprocess_calls"].append({
                    "line": node.lineno,
                    "function": name,
                    "excerpt": excerpt
                })

            if any(
                token in name
                for token in (
                    "torch.load",
                    "joblib.load",
                    "load_state_dict",
                    "predict_proba"
                )
            ):

                result[
                    "direct_model_loading_calls"
                ].append({
                    "line": node.lineno,
                    "call": name
                })

    return result


# ==========================================================
# TEST ML DELLE CINQUE SMALL RNA
# ==========================================================

def test_five():

    for path in (
        MODEL,
        TABLE
    ):

        if not path.is_file():

            return {
                "status": "NOT_VERIFIED",
                "error": (
                    f"File mancante: {path}"
                )
            }

    image_check = call(
        [
            "docker",
            "image",
            "inspect",
            IMAGE
        ],
        timeout=15
    )

    if not image_check["ok"]:

        return {
            "status": "NOT_VERIFIED",
            "error": (
                f"Immagine assente: {IMAGE}"
            )
        }

    command = [
        "docker",
        "run",
        "--rm",
        "--pull", "never",
        "--network", "none",
        "--read-only",
        "--tmpfs",
        "/tmp:rw,nosuid,nodev,size=128m",
        "--cap-drop", "ALL",
        "--security-opt",
        "no-new-privileges",
        "--pids-limit", "128",
        "--memory", "3g",
        "--cpus", "2",
        "--user", "65534:65534",
        "-e", "HOME=/tmp",
        "-v", f"{ROOT}:/project:ro",
        "-w", "/project",
        "--entrypoint", "python3",
        IMAGE,
        "-B",
        "-c",
        TEST_IN_CONTAINER
    ]

    response = call(
        command,
        timeout=180
    )

    if not response["ok"]:

        return {
            "status": "NOT_VERIFIED",
            "error": (
                response.get("stderr")
                or response.get("error")
            ),
            "stdout_tail": (
                response.get(
                    "stdout",
                    ""
                )[-600:]
            )
        }

    for line in response[
        "stdout"
    ].splitlines():

        if line.startswith(
            "ML_TEST_JSON="
        ):

            return json.loads(
                line.partition("=")[2]
            )

    return {
        "status": "NOT_VERIFIED",
        "error": (
            "Il container non ha prodotto "
            "ML_TEST_JSON"
        ),
        "stdout_tail": (
            response["stdout"][-600:]
        )
    }


# ==========================================================
# PROGRAMMA PRINCIPALE
# ==========================================================

def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 22.1")
    print("VERIFICA PULSANTE E 5 SMALL RNA")
    print("=" * 70)

    if not GUI.is_file():
        raise SystemExit(
            f"[ERRORE] GUI assente: {GUI}"
        )

    files = (
        GUI,
        MODEL,
        TABLE
    )

    before = {
        str(path.relative_to(ROOT)): digest(path)
        for path in files
    }

    # ------------------------------------------------------
    # 1. PULSANTE
    # ------------------------------------------------------

    print(
        "\n[1] CODICE EFFETTIVO DEL PULSANTE"
    )

    gui = analyse_button()

    print(
        "Riga pulsante:",
        gui["button_line"]
    )

    print(
        "Chiamate subprocess:",
        len(
            gui["subprocess_calls"]
        )
    )

    for entry in gui[
        "subprocess_calls"
    ]:

        print(
            "\n--- subprocess alla riga",
            entry["line"],
            "---"
        )

        print(
            entry["excerpt"]
        )

    print(
        "\nFunzioni chiamate dal blocco pulsante:"
    )

    relevant = [
        item
        for item in gui["button_calls"]
        if item["call"] not in (
            "str",
            "len",
            "bool",
            "Path",
            "print"
        )
    ]

    for entry in relevant[-35:]:

        print(
            f"  {entry['line']}: "
            f"{entry['call']}"
        )

    # ------------------------------------------------------
    # 2. ML REALE
    # ------------------------------------------------------

    print(
        "\n[2] INFERENZA ML DELLE 5 SMALL RNA"
    )

    ml = test_five()

    print(
        "Esito:",
        ml["status"]
    )

    if ml.get("error"):
        print(
            "Motivo:",
            ml["error"]
        )

    for row in ml.get(
        "records",
        []
    ):

        print(
            f"{row['rna']} | "
            f"classe={row['predicted_class']} | "
            f"p1={row['probability_class_1']} | "
            f"storico={row['historical_label']}"
        )

    # ------------------------------------------------------
    # 3. REPORT
    # ------------------------------------------------------

    after = {
        str(path.relative_to(ROOT)): digest(path)
        for path in files
    }

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "22.1",
        "gui_button_static_trace": gui,
        "five_srna_historical_ml_inference": ml,
        "original_files_unchanged": (
            before == after
        ),
        "sha256_before": before,
        "sha256_after": after,
        "gui_button_to_ml_gnn_runtime_proven": False,
        "independent_scientific_validation": False,
        "training_performed": False,
        "note": (
            "La prova su 5 RNA usa feature "
            "gia calcolate nella tabella storica. "
            "Non dimostra una chiamata ML/GNN "
            "dalla GUI e non costituisce "
            "valutazione indipendente."
        )
    }

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    output = RUNS / (
        f"STEP22_1_INTEGRAZIONE_{stamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report_path = (
        output
        / "RAPPORTO_STEP22_1.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    print(
        "\n[3] REPORT"
    )

    print(
        "Originali invariati:",
        before == after
    )

    print(
        "Percorso:",
        report_path
    )

    print(
        "Collegamento GUI -> ML/GNN: "
        "NON ANCORA DIMOSTRATO"
    )

    if before != after:
        sys.exit(
            "[ERRORE] SHA256 dei file "
            "originali variato durante il test"
        )


if __name__ == "__main__":
    main()