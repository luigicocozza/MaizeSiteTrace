#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 22.4

Collegamento automatico ESTERNO GUI -> ML.

Modalita:
  --once   elabora una cartella esistente
  --watch  monitora nuovi risultati della GUI

Non modifica le GUI, i checkpoint o i dati.
Non esegue training.
Non esegue una nuova inferenza GAE.
"""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import time

from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT / "target_prediction_v33_vscode" / "runs"
)

IMAGE = "ai_primetta_v33_ml_compat:19.5"

MODEL = (
    ROOT
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

HISTORY = (
    ROOT
    / "outputs/prjna419326/"
    "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)

GUI = (
    ROOT / "target_prediction_v33_vscode/app_v33.py"
)

MARKER = "PRIMETTA_AUTO_ML_JSON="


# ============================================================
# MOTORE DI INFERENZA DOCKER
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


items = json.loads(sys.argv[1])

root = Path("/project")

modelpath = (
    root
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

tablepath = (
    root
    / "outputs/prjna419326/"
    "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)


result = {
    "status": "ERROR",
    "records": [],
    "unresolved": [],
    "ml_runtime_executed": False,
    "gae_runtime_executed": False,
    "historical_features": True,
    "sklearn_version": sklearn.__version__,
}


def seq(value):

    sequence = re.sub(
        r"\s+",
        "",
        str(value or "")
    ).upper().replace("U", "T")

    if re.fullmatch(
        r"[ACGTN]+",
        sequence
    ):
        return sequence

    return ""


def tokens(value):

    text = str(
        value or ""
    ).strip().upper()

    values = {
        text
    } if text else set()

    values.update(
        re.findall(
            r"PRJNA\d+_SRNA_RANK_\d+|SRNA_\d+",
            text
        )
    )

    return values


def label(value):

    try:
        return str(
            int(float(value))
        )

    except (TypeError, ValueError):
        return str(value)


try:

    # --------------------------------------------------------
    # VERSIONE COMPATIBILE
    # --------------------------------------------------------

    if sklearn.__version__ != "1.8.0":

        raise RuntimeError(
            "Runtime sklearn diverso da 1.8.0"
        )

    # --------------------------------------------------------
    # CARICAMENTO MODELLO LOCALE
    # --------------------------------------------------------

    with warnings.catch_warnings():

        warnings.simplefilter(
            "error",
            InconsistentVersionWarning
        )

        model = joblib.load(
            modelpath
        )

    # --------------------------------------------------------
    # FEATURE ML
    # --------------------------------------------------------

    features = list(
        getattr(
            model,
            "feature_names_in_",
            []
        )
    )

    if (
        len(features) != 22
        or sum(
            str(feature).startswith("z_")
            for feature in features
        ) != 16
    ):

        raise RuntimeError(
            "Schema 22 feature / "
            "16 embedding non compatibile"
        )

    frame = pd.read_csv(
        tablepath,
        sep="\t",
        dtype=str,
        low_memory=False
    ).fillna("")

    required = [
        "feature_id",
        "sequence",
        *features
    ]

    missing = [
        name
        for name in required
        if name not in frame.columns
    ]

    if missing:

        raise RuntimeError(
            "Colonne mancanti: "
            + str(missing)
        )

    # --------------------------------------------------------
    # INDICI STORICI
    # --------------------------------------------------------

    idcols = [
        column
        for column in (
            "feature_id",
            "v3_rank_id"
        )
        if column in frame
    ]

    knownids = [
        set().union(
            *(
                tokens(row[column])
                for column in idcols
            )
        )
        for _, row in frame.iterrows()
    ]

    knownseq = [
        seq(value)
        for value in frame["sequence"]
    ]

    selected = []

    # --------------------------------------------------------
    # MATCH RNA GUI -> RNA STORICA
    # --------------------------------------------------------

    for item in items:

        qtokens = set().union(
            *(
                tokens(value)
                for value in item["identifiers"]
            )
        )

        qseq = seq(
            item.get("sequence", "")
        )

        byid = {
            index
            for index, values in enumerate(knownids)
            if qtokens & values
        }

        byseq = {
            index
            for index, value in enumerate(knownseq)
            if qseq and value == qseq
        }

        if byid and qseq:

            candidates = byid & byseq
            method = "ID_E_SEQUENZA_ESATTA"

        elif byid:

            candidates = byid
            method = "IDENTIFICATIVO_STORICO"

        else:

            candidates = byseq
            method = "SEQUENZA_ESATTA"

        if len(candidates) != 1:

            result["unresolved"].append({
                "rna_id": item["rna_id"],
                "id_matches": len(byid),
                "sequence_matches": len(byseq),
                "matches": len(candidates),
                "reason": (
                    "AMBIGUO"
                    if candidates
                    else "NESSUN_MATCH_O_CONFLITTO"
                )
            })

            continue

        selected.append((
            item,
            next(iter(candidates)),
            method
        ))

    # --------------------------------------------------------
    # CONTROLLO CORRISPONDENZA UNIVOCA
    # --------------------------------------------------------

    positions = [
        position
        for _, position, _ in selected
    ]

    if len(positions) != len(
        set(positions)
    ):

        raise RuntimeError(
            "RNA GUI differenti associate "
            "allo stesso record storico"
        )

    # --------------------------------------------------------
    # INFERENZA ML
    # --------------------------------------------------------

    if selected:

        X = frame.iloc[
            positions
        ][features].copy()

        for column in features:

            X[column] = pd.to_numeric(
                X[column],
                errors="coerce"
            )

        if not np.isfinite(
            X.to_numpy(dtype=float)
        ).all():

            raise RuntimeError(
                "Feature non finite"
            )

        predicted = model.predict(X)

        probabilities = np.asarray(
            model.predict_proba(X)
        )

        classes = [
            label(value)
            for value in model.classes_
        ]

        if "1" not in classes:

            raise RuntimeError(
                "Classe positiva 1 assente"
            )

        if probabilities.shape != (
            len(selected),
            len(classes)
        ):

            raise RuntimeError(
                "Dimensioni probabilita non valide"
            )

        if not np.isfinite(
            probabilities
        ).all():

            raise RuntimeError(
                "Probabilita non finite"
            )

        if not np.allclose(
            probabilities.sum(axis=1),
            1,
            atol=1e-6
        ):

            raise RuntimeError(
                "Probabilita non normalizzate"
            )

        positive_index = classes.index("1")

        # ----------------------------------------------------
        # RISULTATI
        # ----------------------------------------------------

        for index, (
            item,
            position,
            method
        ) in enumerate(selected):

            historical = frame.iloc[
                position
            ]

            result["records"].append({
                "rna_id": item["rna_id"],

                "feature_id": str(
                    historical["feature_id"]
                ),

                "match_method": method,

                "class_ml": label(
                    predicted[index]
                ),

                "p_class1": float(
                    probabilities[
                        index,
                        positive_index
                    ]
                ),

                "historical_label": label(
                    historical.get(
                        "antileakage_predicted_label",
                        ""
                    )
                )
            })

        result["ml_runtime_executed"] = True

    # --------------------------------------------------------
    # ESITO
    # --------------------------------------------------------

    if len(selected) == len(items):

        result["status"] = (
            "PASS_FULL_HISTORICAL_MATCH"
        )

    else:

        result["status"] = (
            "PARTIAL_HISTORICAL_MATCH"
        )

except Exception as exc:

    result["status"] = "ERROR"

    result["error"] = (
        type(exc).__name__
        + ": "
        + str(exc)
    )


print(
    "PRIMETTA_AUTO_ML_JSON="
    + json.dumps(
        result,
        ensure_ascii=False
    ),
    flush=True
)

'''


# ============================================================
# CALCOLO SHA256
# ============================================================

def sha(path):

    checksum = hashlib.sha256()

    with path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(
                4 * 1024 * 1024
            ),
            b""
        ):

            checksum.update(block)

    return checksum.hexdigest()


# ============================================================
# LETTURA E IDENTIFICAZIONE TABELLE GUI
# ============================================================

def table_info(path):

    if path.stat().st_size > 100 * 1024 * 1024:

        return None

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        columns = set(
            reader.fieldnames or []
        )

        if not columns:
            return None

        return columns, list(reader)


def find_tables(source):

    candidates = []

    for path in sorted(
        source.rglob("*.tsv")
    ):

        if len(
            path.relative_to(source).parts
        ) > 3:

            continue

        try:

            details = table_info(path)

            if details:

                candidates.append((
                    path,
                    *details
                ))

        except (
            OSError,
            UnicodeError,
            csv.Error
        ):

            continue

    qc = [
        (path, rows)
        for path, columns, rows in candidates
        if (
            "status" in columns
            and "rna_id" in columns
            and (
                "sequence" in columns
                or "srna_sequence" in columns
            )
        )
    ]

    summary = [
        (path, rows)
        for path, columns, rows in candidates
        if (
            "n_new_sites" in columns
            and "rna_id" in columns
        )
    ]

    if len(qc) != 1 or len(summary) != 1:

        return None

    return qc[0], summary[0]


# ============================================================
# ESTRAZIONE SEQUENZA
# ============================================================

def get_seq(row):

    for column in (
        "sequence",
        "srna_sequence",
        "rna_sequence",
        "rna_seq",
        "seq"
    ):

        if row.get(column):

            return str(
                row[column]
            ).strip()

    return ""


# ============================================================
# PREPARAZIONE RICHIESTE ML
# ============================================================

def requests_from(
    qc_rows,
    summary_rows
):

    if (
        not qc_rows
        or len(qc_rows) > 10000
    ):

        raise RuntimeError(
            "Numero RNA nullo o superiore a 10000"
        )

    summaries = {}

    for row in summary_rows:

        identifier = str(
            row.get("rna_id", "")
        ).strip()

        if (
            not identifier
            or identifier in summaries
        ):

            raise RuntimeError(
                "ID RNA assente o duplicato "
                "nel summary"
            )

        summaries[identifier] = row

    payload = []

    identifiers_seen = set()

    for row in qc_rows:

        rna_id = str(
            row.get("rna_id", "")
        ).strip()

        if (
            not rna_id
            or rna_id in identifiers_seen
        ):

            raise RuntimeError(
                "ID RNA assente o duplicato "
                "nel QC"
            )

        identifiers_seen.add(
            rna_id
        )

        summary = summaries.get(
            rna_id,
            {}
        )

        first = get_seq(
            row
        ).upper().replace(
            "U",
            "T"
        )

        second = get_seq(
            summary
        ).upper().replace(
            "U",
            "T"
        )

        if (
            first
            and second
            and re.sub(r"\s+", "", first)
            != re.sub(r"\s+", "", second)
        ):

            raise RuntimeError(
                "Sequenze discordanti per "
                + rna_id
            )

        identifiers = [rna_id]

        for record in (
            row,
            summary
        ):

            identifiers.extend(
                str(record[column])
                for column in (
                    "feature_id",
                    "v3_rank_id",
                    "srna_id"
                )
                if record.get(column)
            )

        payload.append({
            "rna_id": rna_id,

            "sequence": first or second,

            "identifiers": identifiers,

            "qc_status": row.get(
                "status",
                ""
            ),

            "n_new_sites": summary.get(
                "n_new_sites",
                ""
            )
        })

    if set(summaries) - identifiers_seen:

        raise RuntimeError(
            "Summary contiene RNA estranee al QC"
        )

    return payload


# ============================================================
# ESECUZIONE MODELLO IN DOCKER
# ============================================================

def run_model(payload):

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
            [
                {
                    key: row[key]
                    for key in (
                        "rna_id",
                        "sequence",
                        "identifiers"
                    )
                }
                for row in payload
            ],
            ensure_ascii=False
        )
    ]

    process = subprocess.run(
        command,
        text=True,
        capture_output=True,
        timeout=180,
        check=False
    )

    if process.returncode:

        raise RuntimeError(
            "Container ML: "
            + process.stderr[-1500:]
        )

    responses = [
        line[len(MARKER):]
        for line in process.stdout.splitlines()
        if line.startswith(MARKER)
    ]

    if len(responses) != 1:

        raise RuntimeError(
            "Risposta ML non valida: "
            + process.stdout[-800:]
        )

    result = json.loads(
        responses[0]
    )

    if result.get(
        "status"
    ) == "ERROR":

        raise RuntimeError(
            result.get(
                "error",
                "Errore nel modello"
            )
        )

    return result


# ============================================================
# SCRITTURA TSV
# ============================================================

def output_tsv(path, rows):

    columns = [
        "rna_id_gui",
        "qc_status",
        "n_new_sites_gui",
        "feature_id_ml",
        "matching_method",
        "class_ml",
        "probability_class1",
        "historical_label",
        "ml_status"
    ]

    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(rows)


# ============================================================
# ELABORAZIONE SINGOLO RISULTATO GUI
# ============================================================

def process(source):

    found = find_tables(
        source
    )

    if found is None:

        raise RuntimeError(
            "Non trovate una tabella QC "
            "e una summary univoche"
        )

    (qpath, qrows), (
        spath,
        srows
    ) = found

    payload = requests_from(
        qrows,
        srows
    )

    protected = [
        GUI,
        MODEL,
        HISTORY,
        qpath,
        spath
    ]

    if any(
        not path.is_file()
        for path in protected
    ):

        raise RuntimeError(
            "Mancano GUI, checkpoint "
            "o tabella storica"
        )

    before = {
        str(path): sha(path)
        for path in protected
    }

    # Inferenza ML
    ml = run_model(
        payload
    )

    after = {
        str(path): sha(path)
        for path in protected
    }

    if before != after:

        raise RuntimeError(
            "File sorgente modificati "
            "durante il collaudo"
        )

    indexed = {
        record["rna_id"]: record
        for record in ml["records"]
    }

    rows = []

    for item in payload:

        prediction = indexed.get(
            item["rna_id"]
        )

        rows.append({
            "rna_id_gui": item["rna_id"],

            "qc_status": item["qc_status"],

            "n_new_sites_gui": item["n_new_sites"],

            "feature_id_ml": (
                prediction["feature_id"]
                if prediction
                else ""
            ),

            "matching_method": (
                prediction["match_method"]
                if prediction
                else ""
            ),

            "class_ml": (
                prediction["class_ml"]
                if prediction
                else ""
            ),

            "probability_class1": (
                f"{prediction['p_class1']:.12f}"
                if prediction
                else ""
            ),

            "historical_label": (
                prediction["historical_label"]
                if prediction
                else ""
            ),

            "ml_status": (
                "ML_FEATURE_STORICHE"
                if prediction
                else "NON_CLASSIFICABILE"
            )
        })

    # --------------------------------------------------------
    # NUOVA CARTELLA RISULTATI
    # --------------------------------------------------------

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    output = RUNS / (
        "STEP22_4_AUTO_GUI_ML_"
        + stamp
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    tsv_path = (
        output
        / "09_gui_ml_auto.tsv"
    )

    output_tsv(
        tsv_path,
        rows
    )

    # --------------------------------------------------------
    # MANIFESTO
    # --------------------------------------------------------

    manifest = {
        "step": "22.4",

        "status": ml["status"],

        "source": str(source),

        "source_files": {
            "qc": str(qpath),
            "summary": str(spath)
        },

        "input_count": len(payload),

        "matched_count": len(
            ml["records"]
        ),

        "unresolved": ml["unresolved"],

        "model_image": IMAGE,

        "ml_runtime_executed": ml[
            "ml_runtime_executed"
        ],

        "gae_runtime_executed": False,

        "historical_22_features": True,

        "gui_source_modified": False,

        "gui_displays_ml_results": False,

        "independent_validation": False,

        "training_performed": False,

        "sha256_inputs": before,

        "notes": [
            (
                "La GUI produce risultati; "
                "il monitor esterno li elabora "
                "successivamente."
            ),
            (
                "Nessuna modifica grafica "
                "o richiamo ML dal codice Streamlit."
            ),
            (
                "Le RNA senza feature storiche "
                "non vengono classificate."
            ),
            (
                "Il GAE non viene rieseguito."
            )
        ]
    }

    manifest_path = (
        output
        / "10_manifest_gui_ml_auto.json"
    )

    manifest_path.write_text(
        json.dumps(
            manifest,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    # --------------------------------------------------------
    # RIEPILOGO
    # --------------------------------------------------------

    print(
        f"[OK] {ml['status']}: "
        f"{len(ml['records'])}/{len(payload)} RNA"
    )

    print(
        "[OUTPUT]",
        output
    )

    for row in rows:

        print(
            " ",
            row["rna_id_gui"][:70],
            "| ML:",
            row["class_ml"] or "N/D",
            "| siti:",
            row["n_new_sites_gui"]
        )

    return output


# ============================================================
# PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    modes = parser.add_mutually_exclusive_group(
        required=True
    )

    modes.add_argument(
        "--once",
        action="store_true",
        help="Elabora una cartella esistente"
    )

    modes.add_argument(
        "--watch",
        action="store_true",
        help="Sorveglia nuovi risultati GUI"
    )

    parser.add_argument(
        "--source",
        type=Path,
        help="Cartella evidenze per --once"
    )

    parser.add_argument(
        "--interval",
        type=int,
        default=5,
        help="Intervallo di controllo in secondi"
    )

    args = parser.parse_args()

    if args.watch and args.source:

        parser.error(
            "--source vale solo con --once"
        )

    if not 2 <= args.interval <= 300:

        parser.error(
            "--interval deve essere tra 2 e 300"
        )

    print(
        "AI PRIMETTA V3.3 - STEP 22.4"
    )

    print(
        "COLLEGAMENTO AUTOMATICO ESTERNO GUI -> ML"
    )

    # --------------------------------------------------------
    # CHECK FILE ORIGINALI
    # --------------------------------------------------------

    for path in (
        MODEL,
        HISTORY,
        GUI
    ):

        if not path.is_file():

            raise RuntimeError(
                "File essenziale non trovato: "
                + str(path)
            )

    # --------------------------------------------------------
    # CHECK DOCKER
    # --------------------------------------------------------

    check = subprocess.run(
        [
            "docker",
            "image",
            "inspect",
            IMAGE
        ],
        capture_output=True,
        text=True,
        timeout=20,
        check=False
    )

    if check.returncode:

        raise RuntimeError(
            "Immagine Docker non disponibile: "
            + IMAGE
        )

    RUNS.mkdir(
        parents=True,
        exist_ok=True
    )

    # ========================================================
    # MODALITA' --once
    # ========================================================

    if args.once:

        if args.source:

            source = (
                args.source
                if args.source.is_absolute()
                else ROOT / args.source
            ).resolve()

        else:

            candidates = sorted(
                (
                    path
                    for path in RUNS.glob("result_*")
                    if path.is_dir()
                ),
                key=lambda path: path.stat().st_mtime,
                reverse=True
            )

            source = (
                candidates[0]
                if candidates
                else RUNS / "STEP20_EVIDENZE_GUI"
            )

        if (
            not source.is_dir()
            or not source.is_relative_to(
                RUNS.resolve()
            )
        ):

            raise RuntimeError(
                "Cartella sorgente assente "
                "o esterna a runs"
            )

        process(source)

        return

    # ========================================================
    # MODALITA' --watch
    # ========================================================

    seen = {
        str(path.resolve())
        for path in RUNS.glob("result_*")
        if path.is_dir()
    }

    print(
        "[WATCH] Monitoraggio avviato."
    )

    print(
        "[WATCH] Verranno elaborate solo "
        "NUOVE cartelle result_*."
    )

    print(
        "[WATCH] Ctrl+C per interrompere."
    )

    print(
        "[WATCH] Ora esegui una nuova analisi "
        "dalla GUI V3.3 sulla porta 8502."
    )

    while True:

        for path in sorted(
            RUNS.glob("result_*")
        ):

            if not path.is_dir():

                continue

            identifier = str(
                path.resolve()
            )

            if identifier in seen:

                continue

            found = find_tables(
                path
            )

            if found is None:

                continue

            sources = [
                found[0][0],
                found[1][0]
            ]

            last_modified = max(
                source.stat().st_mtime
                for source in sources
            )

            if (
                time.time() - last_modified
                < args.interval
            ):

                continue

            print(
                "[NUOVO RISULTATO]",
                path
            )

            try:

                process(path)

            except Exception as exc:

                print(
                    "[ERRORE]",
                    type(exc).__name__,
                    str(exc),
                    flush=True
                )

            seen.add(
                identifier
            )

        time.sleep(
            args.interval
        )


# ============================================================
# AVVIO
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        print(
            "\n[STOP] Monitoraggio terminato."
        )

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__,
            str(exc),
            file=sys.stderr
        )

        sys.exit(1)