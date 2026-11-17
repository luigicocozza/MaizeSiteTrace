#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 22.3 CORRETTO

Unisce le evidenze della GUI V3.3 con
l'inferenza ML sulle feature storiche.

Non modifica:
- GUI V3.2
- GUI V3.3
- checkpoint ML/GNN
- dati scientifici originali
- container Docker persistenti

Non esegue training né nuova inferenza GAE.
"""

import csv
import hashlib
import json
import subprocess
import sys

from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

BASE = ROOT / "target_prediction_v33_vscode"

RUNS = BASE / "runs"

EVIDENCE = RUNS / "STEP20_EVIDENZE_GUI"

IMAGE = "ai_primetta_v33_ml_compat:19.5"

MODEL = (
    ROOT
    / "models/srna_supervised_ml_v1_antileakage"
    / "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

HISTORY = (
    ROOT
    / "outputs/prjna419326"
    / "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)

GUI = BASE / "app_v33.py"

EVIDENCE_FILES = (
    "01_rna_qc.tsv",
    "02_alignment_hypotheses.tsv",
    "03_historical_12mer_hypotheses.tsv",
    "04_gene_target_candidates.tsv",
    "05_rna_candidate_summary.tsv",
    "06_provenance_manifest.json",
)


# ============================================================
# 2. MOTORE ML - CONTAINER ISOLATO
# ============================================================

WORKER = r'''

import json
import re
import sys
import warnings

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn

from sklearn.exceptions import InconsistentVersionWarning


requests = json.loads(sys.argv[1])

root = Path("/project")

model_path = (
    root
    / "models/srna_supervised_ml_v1_antileakage"
    / "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

table_path = (
    root
    / "outputs/prjna419326"
    / "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)


result = {
    "status": "ERROR",
    "sklearn_version": sklearn.__version__,
    "records": [],
    "unresolved": [],
    "feature_names": [],
    "historical_rows": 0,
    "historical_sample_ids": [],
    "runtime_ml_inference": False,
    "runtime_gae_inference": False,
    "gui_button_called_model": False,
}


# ------------------------------------------------------------
# NORMALIZZAZIONE SEQUENZA
# ------------------------------------------------------------

def clean_sequence(value):

    sequence = re.sub(
        r"\s+",
        "",
        str(value or "")
    ).upper().replace("U", "T")

    if sequence in ("", "NAN", "NONE"):
        return ""

    if not re.fullmatch(
        r"[ACGTN]+",
        sequence
    ):
        return ""

    return sequence


# ------------------------------------------------------------
# IDENTIFICATIVI STORICI
# ------------------------------------------------------------

def id_tokens(value):

    text = str(
        value or ""
    ).strip().upper()

    values = {text} if text else set()

    values.update(
        re.findall(
            r"PRJNA\d+_SRNA_RANK_\d+|SRNA_\d+",
            text
        )
    )

    rank = re.search(
        r"(?:PRJNA\d+_SRNA_)?RANK_0*(\d+)",
        text
    )

    if rank:
        values.add(
            "RANK#" + str(
                int(rank.group(1))
            )
        )

    if text.isdigit():
        values.add(
            "RANK#" + str(
                int(text)
            )
        )

    return values


def normalize_label(value):

    try:
        return str(
            int(float(value))
        )

    except (ValueError, TypeError):
        return str(value)


# ============================================================
# 3. CARICAMENTO MODELLO E FEATURE
# ============================================================

try:

    if sklearn.__version__ != "1.8.0":

        raise RuntimeError(
            "Versione scikit-learn non compatibile "
            "con il checkpoint verificato"
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

    n_embeddings = sum(
        str(feature).startswith("z_")
        for feature in features
    )

    if (
        len(features) != 22
        or n_embeddings != 16
    ):

        raise RuntimeError(
            "Schema ML inatteso: "
            "richieste 22 feature, "
            "di cui 16 embedding z_"
        )

    result["feature_names"] = features

    # --------------------------------------------------------
    # LETTURA TABELLA STORICA
    # --------------------------------------------------------

    frame = pd.read_csv(
        table_path,
        sep="\t",
        dtype=str,
        low_memory=False
    ).fillna("")

    required_columns = [
        "feature_id",
        "sequence"
    ] + features

    missing = [
        name
        for name in required_columns
        if name not in frame.columns
    ]

    if missing:

        raise RuntimeError(
            "Colonne storiche mancanti: "
            + str(missing)
        )

    result["historical_rows"] = len(frame)

    result["historical_sample_ids"] = [
        {
            "feature_id": str(
                row.get("feature_id", "")
            ),
            "v3_rank_id": str(
                row.get("v3_rank_id", "")
            ),
        }
        for _, row in frame.head(5).iterrows()
    ]

    # --------------------------------------------------------
    # COSTRUZIONE INDICI STORICI
    # --------------------------------------------------------

    id_columns = [
        name
        for name in (
            "feature_id",
            "v3_rank_id"
        )
        if name in frame.columns
    ]

    historical_ids = []

    for _, row in frame.iterrows():

        values = set()

        for column in id_columns:

            values.update(
                id_tokens(
                    row[column]
                )
            )

        historical_ids.append(
            values
        )

    historical_sequences = [
        clean_sequence(value)
        for value in frame["sequence"]
    ]

    # ========================================================
    # 4. ASSOCIAZIONE RNA GUI -> MODELLO ML
    # ========================================================

    matched = []

    for request in requests:

        rna_id = str(
            request["rna_id"]
        )

        sequence = clean_sequence(
            request.get(
                "sequence",
                ""
            )
        )

        identifiers = set()

        for identifier in request.get(
            "identifiers",
            [rna_id]
        ):

            identifiers.update(
                id_tokens(identifier)
            )

        reason = None

        # ----------------------------------------------------
        # MATCH PER IDENTIFICATIVO
        # ----------------------------------------------------

        by_id = {
            index
            for index, values in enumerate(
                historical_ids
            )
            if identifiers.intersection(
                values
            )
        }

        # ----------------------------------------------------
        # MATCH PER SEQUENZA
        # ----------------------------------------------------

        by_sequence = {
            index
            for index, candidate in enumerate(
                historical_sequences
            )
            if (
                sequence
                and candidate == sequence
            )
        }

        # ----------------------------------------------------
        # CONCORDANZA TRA CRITERI
        # ----------------------------------------------------

        if by_id and by_sequence:

            candidates = (
                by_id & by_sequence
            )

            method = (
                "ID_E_SEQUENZA_ESATTA"
            )

            if not candidates:
                reason = (
                    "CONFLITTO_TRA_IDENTIFICATIVO_"
                    "E_SEQUENZA"
                )

        elif by_id:

            candidates = by_id

            method = (
                "IDENTIFICATIVO_STORICO"
            )

            if sequence:

                candidates = {
                    index
                    for index in candidates
                    if (
                        historical_sequences[index]
                        == sequence
                    )
                }

                if not candidates:
                    reason = (
                        "ID_CON_SEQUENZA_DIFFERENTE"
                    )

        else:

            candidates = by_sequence

            method = (
                "SEQUENZA_ESATTA"
            )

        # ----------------------------------------------------
        # BLOCCO ASSOCIAZIONI AMBIGUE
        # ----------------------------------------------------

        if len(candidates) != 1:

            result["unresolved"].append({

                "rna_id": rna_id,

                "gui_sequence_present": bool(
                    sequence
                ),

                "id_matches": len(
                    by_id
                ),

                "sequence_matches": len(
                    by_sequence
                ),

                "intersection_matches": len(
                    by_id & by_sequence
                ),

                "reason": (
                    reason
                    or (
                        "NESSUNA_CORRISPONDENZA"
                        if not candidates
                        else "CORRISPONDENZA_NON_UNIVOCA"
                    )
                ),

                "potential_feature_ids": [
                    str(
                        frame.iloc[index][
                            "feature_id"
                        ]
                    )
                    for index in sorted(
                        candidates
                    )[:10]
                ]
            })

            continue

        index = next(
            iter(candidates)
        )

        matched.append(
            (
                request,
                index,
                method
            )
        )

    # --------------------------------------------------------
    # VERIFICA UNIVOCITA'
    # --------------------------------------------------------

    matched_indices = [
        index
        for _, index, _ in matched
    ]

    if len(matched_indices) != len(
        set(matched_indices)
    ):

        raise RuntimeError(
            "RNA GUI differenti associate "
            "allo stesso record ML"
        )

    # ========================================================
    # 5. INFERENZA DEL CLASSIFICATORE
    # ========================================================

    if matched:

        X = frame.iloc[
            matched_indices
        ].loc[
            :,
            features
        ].copy()

        for feature in features:

            X[feature] = pd.to_numeric(
                X[feature],
                errors="coerce"
            )

        matrix = X.to_numpy(
            dtype=float
        )

        if not np.isfinite(
            matrix
        ).all():

            raise RuntimeError(
                "Feature ML non numeriche "
                "o non finite"
            )

        predictions = model.predict(
            X
        )

        probabilities = np.asarray(
            model.predict_proba(X)
        )

        classes = [
            normalize_label(value)
            for value in model.classes_
        ]

        if "1" not in classes:

            raise RuntimeError(
                "Classe positiva 1 assente"
            )

        if probabilities.shape != (
            len(matched),
            len(classes)
        ):

            raise RuntimeError(
                "Dimensioni probabilita errate"
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

        positive_index = classes.index(
            "1"
        )

        # ----------------------------------------------------
        # RISULTATI DELLE RNA
        # ----------------------------------------------------

        for position, (
            request,
            index,
            method
        ) in enumerate(matched):

            historical = frame.iloc[
                index
            ]

            label = normalize_label(
                predictions[position]
            )

            historical_label = normalize_label(
                historical.get(
                    "antileakage_predicted_label",
                    ""
                )
            )

            value = historical.get(
                "antileakage_probability_positive",
                ""
            )

            historical_probability = (
                float(value)
                if value != ""
                else None
            )

            probability = float(
                probabilities[
                    position,
                    positive_index
                ]
            )

            result["records"].append({

                "rna_id": request["rna_id"],

                "feature_id": str(
                    historical["feature_id"]
                ),

                "v3_rank_id": str(
                    historical.get(
                        "v3_rank_id",
                        ""
                    )
                ),

                "matching_method": method,

                "ml_predicted_label": label,

                "ml_probability_class1": probability,

                "historical_label": (
                    historical_label
                ),

                "historical_probability_class1": (
                    historical_probability
                ),

                "historical_label_reproduced": (
                    label == historical_label
                ),

                "historical_probability_delta": (
                    abs(
                        probability
                        - historical_probability
                    )
                    if historical_probability is not None
                    else None
                ),
            })

        result["runtime_ml_inference"] = True

    # ========================================================
    # 6. ESITO DEL CONTAINER
    # ========================================================

    if len(result["records"]) == len(
        requests
    ):

        result["status"] = (
            "PASS_FULL_HISTORICAL_GUI_ML_MATCH"
        )

    else:

        result["status"] = (
            "PARTIAL_HISTORICAL_GUI_ML_MATCH"
        )

except Exception as exc:

    result["status"] = "ERROR"

    result["error"] = (
        type(exc).__name__
        + ": "
        + str(exc)
    )


print(
    "PRIMETTA_STEP223_JSON="
    + json.dumps(
        result,
        ensure_ascii=False
    ),
    flush=True
)

'''


# ============================================================
# 7. SHA256
# ============================================================

def checksum(path):

    h = hashlib.sha256()

    with path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(
                4 * 1024 * 1024
            ),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# ============================================================
# 8. LETTURA TSV
# ============================================================

def read_table(path):

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        if not reader.fieldnames:

            raise RuntimeError(
                "Intestazione TSV mancante: "
                + str(path)
            )

        return list(reader)


# ============================================================
# 9. IDENTIFICAZIONE RNA
# ============================================================

def row_id(row):

    columns = (
        "rna_id",
        "srna_id",
        "id",
        "feature_id",
        "v3_rank_id"
    )

    for column in columns:

        value = str(
            row.get(
                column,
                ""
            ) or ""
        ).strip()

        if value:
            return value

    raise RuntimeError(
        "Identificativo RNA assente. "
        "Colonne: "
        + str(list(row))
    )


def row_sequence(row):

    columns = (
        "sequence",
        "srna_sequence",
        "rna_sequence",
        "rna_seq",
        "seq"
    )

    for column in columns:

        value = str(
            row.get(
                column,
                ""
            ) or ""
        ).strip()

        if value:
            return value

    return ""


# ============================================================
# 10. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    print("=" * 70)

    print(
        "AI PRIMETTA V3.3 | STEP 22.3 CORRETTO"
    )

    print(
        "INTEGRAZIONE EVIDENZE GUI + ML"
    )

    print("=" * 70)

    # --------------------------------------------------------
    # FILE PROTETTI
    # --------------------------------------------------------

    protected = [
        GUI,
        MODEL,
        HISTORY
    ] + [
        EVIDENCE / name
        for name in EVIDENCE_FILES
    ]

    missing = [
        str(path)
        for path in protected
        if not path.is_file()
    ]

    if missing:

        raise RuntimeError(
            "File necessari mancanti:\n"
            + "\n".join(missing)
        )

    before = {
        str(path.relative_to(ROOT)): checksum(path)
        for path in protected
    }

    # --------------------------------------------------------
    # INPUT GUI
    # --------------------------------------------------------

    qc = read_table(
        EVIDENCE / "01_rna_qc.tsv"
    )

    summaries = read_table(
        EVIDENCE / "05_rna_candidate_summary.tsv"
    )

    if not qc:

        raise RuntimeError(
            "Nessuna RNA nella tabella QC"
        )

    counts = {
        name: len(
            read_table(
                EVIDENCE / name
            )
        )
        for name in EVIDENCE_FILES
        if name.endswith(".tsv")
    }

    summary_map = {}

    for row in summaries:

        identifier = row_id(row)

        if identifier in summary_map:

            raise RuntimeError(
                "RNA duplicata nel summary: "
                + identifier
            )

        summary_map[identifier] = row

    # --------------------------------------------------------
    # COSTRUZIONE RICHIESTE ML
    # --------------------------------------------------------

    requests = []

    for row in qc:

        identifier = row_id(row)

        if identifier in {
            item["rna_id"]
            for item in requests
        }:

            raise RuntimeError(
                "RNA duplicata nel QC: "
                + identifier
            )

        summary = summary_map.get(
            identifier,
            {}
        )

        sequences = [
            row_sequence(row),
            row_sequence(summary)
        ]

        sequences = [
            value
            for value in sequences
            if value
        ]

        if len(sequences) == 2:

            first = sequences[0].upper().replace(
                "U",
                "T"
            )

            second = sequences[1].upper().replace(
                "U",
                "T"
            )

            if first != second:

                raise RuntimeError(
                    "Sequenze GUI non concordanti: "
                    + identifier
                )

        identifiers = [
            identifier
        ]

        for source in (
            row,
            summary
        ):

            for column in (
                "feature_id",
                "v3_rank_id",
                "srna_id"
            ):

                value = source.get(
                    column
                )

                if value:
                    identifiers.append(
                        str(value)
                    )

        requests.append({

            "rna_id": identifier,

            "sequence": (
                sequences[0]
                if sequences
                else ""
            ),

            "identifiers": identifiers
        })

    # --------------------------------------------------------
    # DIAGNOSTICA INPUT
    # --------------------------------------------------------

    print()

    print(
        "[1] RNA IDENTIFICATE NELLA GUI"
    )

    for request in requests:

        print(
            " ",
            request["rna_id"],
            "| sequenza disponibile:",
            bool(
                request["sequence"]
            )
        )

    print(
        "\nRNA totali:",
        len(requests)
    )

    print(
        "Siti RNA-trascritto:",
        counts["02_alignment_hypotheses.tsv"]
    )

    # --------------------------------------------------------
    # VERIFICA IMMAGINE DOCKER
    # --------------------------------------------------------

    image_check = subprocess.run(
        [
            "docker",
            "image",
            "inspect",
            IMAGE
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=20
    )

    if image_check.returncode != 0:

        raise RuntimeError(
            "Immagine Docker ML non disponibile: "
            + IMAGE
        )

    # ========================================================
    # 11. AVVIO INFERENZA ML ISOLATA
    # ========================================================

    print()

    print(
        "[2] INFERENZA ML - CONTAINER ISOLATO"
    )

    print(
        "Nessun training. "
        "Repository montato in sola lettura."
    )

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
        "/tmp:rw,nosuid,nodev,size=128m",

        "--cap-drop",
        "ALL",

        "--security-opt",
        "no-new-privileges",

        "--pids-limit",
        "128",

        "--memory",
        "3g",

        "--cpus",
        "2",

        "--user",
        "65534:65534",

        "-e",
        "HOME=/tmp",

        "-v",
        f"{ROOT}:/project:ro",

        "-w",
        "/project",

        "--entrypoint",
        "python3",

        IMAGE,

        "-B",
        "-c",
        WORKER,

        json.dumps(
            requests,
            ensure_ascii=False
        )
    ]

    process = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        timeout=180
    )

    if process.returncode != 0:

        raise RuntimeError(
            "Container ML terminato con errore:\n"
            + process.stderr[-1800:]
        )

    # --------------------------------------------------------
    # RISPOSTA DEL MODELLO
    # --------------------------------------------------------

    marker = "PRIMETTA_STEP223_JSON="

    responses = [
        line[len(marker):]
        for line in process.stdout.splitlines()
        if line.startswith(marker)
    ]

    if len(responses) != 1:

        raise RuntimeError(
            "Risposta ML assente o ambigua:\n"
            + process.stdout[-800:]
            + process.stderr[-500:]
        )

    result = json.loads(
        responses[0]
    )

    # --------------------------------------------------------
    # VERIFICA INTEGRITA' ORIGINALI
    # --------------------------------------------------------

    after = {
        str(path.relative_to(ROOT)): checksum(path)
        for path in protected
    }

    if before != after:

        raise RuntimeError(
            "File originali modificati durante "
            "il controllo. Nessun output salvato."
        )

    # ========================================================
    # 12. COSTRUZIONE RISULTATI INTEGRATI
    # ========================================================

    matches = {
        record["rna_id"]: record
        for record in result["records"]
    }

    report_rows = []

    qc_map = {
        row_id(row): row
        for row in qc
    }

    for request in requests:

        identifier = request["rna_id"]

        prediction = matches.get(
            identifier,
            {}
        )

        quality = qc_map[
            identifier
        ]

        summary = summary_map.get(
            identifier,
            {}
        )

        probability = (
            f"{prediction['ml_probability_class1']:.12f}"
            if prediction
            else ""
        )

        report_rows.append({

            "rna_id_gui": identifier,

            "feature_id_ml": prediction.get(
                "feature_id",
                ""
            ),

            "v3_rank_id_ml": prediction.get(
                "v3_rank_id",
                ""
            ),

            "metodo_associazione": prediction.get(
                "matching_method",
                "NON_RISOLTO"
            ),

            "qc_status": quality.get(
                "status",
                ""
            ),

            "siti_target_gui": summary.get(
                "n_new_sites",
                "NON_DISPONIBILE"
            ),

            "classe_ml": prediction.get(
                "ml_predicted_label",
                ""
            ),

            "probabilita_classe_1": probability,

            "riproduce_classe_storica": prediction.get(
                "historical_label_reproduced",
                ""
            ),

            "tipo_inferenza": (
                "ML_SU_FEATURE_STORICHE"
                if prediction
                else "NESSUNA_INFERENZA"
            )
        })

    # ========================================================
    # 13. SALVATAGGIO NUOVI RISULTATI
    # ========================================================

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    output = RUNS / (
        "STEP22_3_INTEGRAZIONE_CORRETTA_"
        + stamp
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    # --------------------------------------------------------
    # TABELLA TSV
    # --------------------------------------------------------

    tsv_path = output / (
        "07_integrazione_gui_ml_storica.tsv"
    )

    with tsv_path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                report_rows[0]
            ),
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(
            report_rows
        )

    # --------------------------------------------------------
    # MANIFESTO JSON
    # --------------------------------------------------------

    report = {

        "project": "AI PRIMETTA V3.3",

        "step": "22.3_CORRETTO",

        "status": result["status"],

        "rna_gui": len(requests),

        "rna_ml_matched": len(
            result["records"]
        ),

        "gui_table_counts": counts,

        "ml_probe": result,

        "original_sha256_before_and_after": before,

        "originals_unchanged": True,

        "training_performed": False,

        "runtime_gae_inference": False,

        "gui_button_calls_ml": False,

        "independent_validation": False,

        "scientific_limits": [

            (
                "Predizioni eseguite sulle "
                "22 feature storiche, inclusi "
                "16 embedding GNN gia calcolati."
            ),

            (
                "Il pulsante della GUI V3.3 "
                "non e collegato automaticamente "
                "al classificatore."
            ),

            (
                "La probabilita della classe 1 "
                "non dimostra causalita biologica "
                "o generalizzazione predittiva."
            ),

            (
                "Zero siti target nelle condizioni "
                "dello STEP20 non contraddice "
                "una classificazione ML positiva."
            )
        ]
    }

    manifest = output / (
        "08_manifest_integrazione_backend.json"
    )

    manifest.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    # ========================================================
    # 14. STAMPA RISULTATI
    # ========================================================

    print()

    print(
        "[3] ESITO:",
        result["status"]
    )

    for row in report_rows:

        print(
            " ",
            row["rna_id_gui"],

            "->",

            row["feature_id_ml"]
            or "NESSUN MATCH",

            "| metodo:",
            row["metodo_associazione"],

            "| classe:",
            row["classe_ml"]
            or "N/D",

            "| siti GUI:",
            row["siti_target_gui"]
        )

    # --------------------------------------------------------
    # DIAGNOSTICA MANCATI MATCH
    # --------------------------------------------------------

    if result.get("unresolved"):

        print()

        print(
            "[DIAGNOSTICA] RNA NON ASSOCIATE:"
        )

        for issue in result["unresolved"]:

            print(
                json.dumps(
                    issue,
                    ensure_ascii=False
                )
            )

    if result.get("error"):

        print(
            "[ERRORE MODELLO]",
            result["error"]
        )

    print()

    print(
        "[4] ORIGINALI INVARIATI: SI"
    )

    print(
        "[5] OUTPUT:",
        output
    )

    print(
        "TABELLA:",
        tsv_path.name
    )

    print(
        "MANIFESTO:",
        manifest.name
    )

    print()

    print(
        "[LIMITE] Integrazione dei risultati, "
        "non collegamento automatico del pulsante."
    )

    if result["status"] != (
        "PASS_FULL_HISTORICAL_GUI_ML_MATCH"
    ):

        raise SystemExit(2)

    print()

    print("=" * 70)

    print(
        "[OK] STEP 22.3 CORRETTO COMPLETATO"
    )

    print("=" * 70)


# ============================================================
# 15. AVVIO
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        print(
            "[INTERROTTO]"
        )

        sys.exit(130)

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__ + ":",
            exc,
            file=sys.stderr
        )

        sys.exit(1)