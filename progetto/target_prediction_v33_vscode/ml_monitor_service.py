#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - SERVIZIO MONITOR ML

Sorveglia le nuove elaborazioni della GUI V3.3.

Utilizza:
- Modello supervised antileakage originale
- 22 feature storiche
- 16 embedding GNN storici

Non esegue un nuovo GAE.
Non modifica le GUI.
Non esegue training.
"""

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time

from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn


# ============================================================
# A. PERCORSI NEL CONTAINER DOCKER
# ============================================================

PROJECT = Path("/project")

RUNS = Path("/runs")

OUT = RUNS / "STEP22_5_MONITOR"

MODEL = (
    PROJECT
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

HISTORY = (
    PROJECT
    / "outputs/prjna419326/"
    "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)


# ============================================================
# B. FUNZIONI DI SUPPORTO
# ============================================================

def say(*args):
    print(*args, flush=True)


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for chunk in iter(
            lambda: f.read(4 * 1024 * 1024),
            b""
        ):
            h.update(chunk)

    return h.hexdigest()


# ============================================================
# C. LETTURA DELLE TABELLE
# ============================================================

def read_tsv(path):

    if path.stat().st_size > 100 * 1024 * 1024:

        raise ValueError(
            "TSV GUI oltre 100 MB: " + str(path)
        )

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        if not reader.fieldnames:

            raise ValueError(
                "TSV senza intestazione: " + str(path)
            )

        return list(reader), set(reader.fieldnames)


# ============================================================
# D. RICERCA RISULTATI GUI
# ============================================================

def locate(source):

    candidates = []

    for path in sorted(source.rglob("*.tsv")):

        if len(path.relative_to(source).parts) > 3:
            continue

        try:

            rows, cols = read_tsv(path)

        except (
            OSError,
            ValueError,
            UnicodeError,
            csv.Error
        ):
            continue

        candidates.append(
            (path, rows, cols)
        )

    qc = [
        (p, r)
        for p, r, c in candidates
        if (
            {"rna_id", "status"} <= c
            and (
                "sequence" in c
                or "srna_sequence" in c
            )
        )
    ]

    sm = [
        (p, r)
        for p, r, c in candidates
        if {
            "rna_id",
            "n_new_sites"
        } <= c
    ]

    if len(qc) != 1 or len(sm) != 1:

        raise ValueError(
            "QC/summary non univoci "
            "o non ancora disponibili"
        )

    return qc[0], sm[0]


# ============================================================
# E. NORMALIZZAZIONE SEQUENZE RNA
# ============================================================

def seq(value):

    s = re.sub(
        r"\s+",
        "",
        str(value or "")
    ).upper().replace("U", "T")

    return (
        s
        if re.fullmatch(r"[ACGTN]+", s)
        else ""
    )


# ============================================================
# F. NORMALIZZAZIONE IDENTIFICATIVI
# ============================================================

def ids(value):

    text = str(
        value or ""
    ).strip().upper()

    found = {text} if text else set()

    found.update(
        re.findall(
            r"PRJNA\d+_SRNA_RANK_\d+|SRNA_\d+",
            text
        )
    )

    return found


def label(value):

    try:
        return str(int(float(value)))

    except (TypeError, ValueError):
        return str(value)


# ============================================================
# G. PREPARAZIONE DATI RNA
# ============================================================

def prepare(qrows, srows):

    summary = {}

    for row in srows:

        name = str(
            row.get("rna_id", "")
        ).strip()

        if not name or name in summary:

            raise ValueError(
                "ID vuoto/duplicato nel summary"
            )

        summary[name] = row

    if not qrows or len(qrows) > 10000:

        raise ValueError(
            "Numero RNA nel QC non valido"
        )

    payload = []

    seen = set()

    for row in qrows:

        name = str(
            row.get("rna_id", "")
        ).strip()

        if not name or name in seen:

            raise ValueError(
                "ID vuoto/duplicato nel QC"
            )

        seen.add(name)

        second = summary.get(
            name,
            {}
        )

        s1 = seq(
            row.get("sequence")
            or row.get("srna_sequence")
        )

        s2 = seq(
            second.get("sequence")
            or second.get("srna_sequence")
        )

        if s1 and s2 and s1 != s2:

            raise ValueError(
                "Sequenza QC/summary discordante per "
                + name
            )

        tokens = ids(name)

        for obj in (row, second):

            for col in (
                "feature_id",
                "v3_rank_id",
                "srna_id"
            ):

                tokens.update(
                    ids(obj.get(col))
                )

        payload.append({

            "rna_id": name,

            "sequence": s1 or s2,

            "tokens": tokens,

            "qc_status": row.get(
                "status",
                ""
            ),

            "n_new_sites": second.get(
                "n_new_sites",
                ""
            )
        })

    if set(summary) - seen:

        raise ValueError(
            "Summary contiene RNA assenti dal QC"
        )

    return payload


# ============================================================
# H. CARICAMENTO DEL MODELLO ML
# ============================================================

def load_model():

    if sklearn.__version__ != "1.8.0":

        raise RuntimeError(
            "Versione sklearn inattesa: "
            + sklearn.__version__
        )

    # Caricare esclusivamente il checkpoint
    # attendibile del repository locale.

    model = joblib.load(MODEL)

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
            str(c).startswith("z_")
            for c in features
        ) != 16
    ):

        raise RuntimeError(
            "Schema checkpoint incompatibile "
            "(22 feature, 16 z_)"
        )

    df = pd.read_csv(
        HISTORY,
        sep="\t",
        dtype=str,
        low_memory=False
    ).fillna("")

    required = [
        "feature_id",
        "sequence",
        *features
    ]

    if any(
        c not in df.columns
        for c in required
    ):

        raise RuntimeError(
            "Colonne richieste assenti "
            "dalla tabella ML"
        )

    historical_ids = []

    for _, row in df.iterrows():

        values = (
            ids(row.get("feature_id"))
            | ids(row.get("v3_rank_id"))
        )

        historical_ids.append(values)

    historical_seq = [
        seq(s)
        for s in df["sequence"]
    ]

    say(
        "[MODEL] sklearn",
        sklearn.__version__,
        "| record",
        len(df),
        "| feature",
        len(features)
    )

    return (
        model,
        features,
        df,
        historical_ids,
        historical_seq
    )


# ============================================================
# I. CLASSIFICAZIONE MACHINE LEARNING
# ============================================================

def classify(payload, loaded):

    (
        model,
        features,
        df,
        historical_ids,
        historical_seq
    ) = loaded

    chosen = []

    unresolved = []

    # --------------------------------------------------------
    # ASSOCIAZIONE RNA GUI -> DATI ML
    # --------------------------------------------------------

    for item in payload:

        by_id = {
            i
            for i, values in enumerate(historical_ids)
            if item["tokens"] & values
        }

        by_seq = {
            i
            for i, sequence in enumerate(historical_seq)
            if (
                item["sequence"]
                and sequence == item["sequence"]
            )
        }

        if by_id and item["sequence"]:

            candidates = by_id & by_seq

            method = "ID_E_SEQUENZA_ESATTA"

        elif by_id:

            candidates = by_id

            method = "IDENTIFICATIVO_STORICO"

        else:

            candidates = by_seq

            method = "SEQUENZA_ESATTA"

        if len(candidates) != 1:

            unresolved.append({

                "rna_id": item["rna_id"],

                "id_matches": len(by_id),

                "sequence_matches": len(by_seq),

                "candidate_count": len(candidates),

                "reason": (
                    "AMBIGUO"
                    if candidates
                    else "NESSUN_MATCH_O_CONFLITTO"
                )
            })

            continue

        chosen.append(
            (
                item,
                next(iter(candidates)),
                method
            )
        )

    # --------------------------------------------------------
    # CONTROLLO UNIVOCITA'
    # --------------------------------------------------------

    positions = [
        index
        for _, index, _ in chosen
    ]

    if len(positions) != len(set(positions)):

        raise RuntimeError(
            "RNA GUI diverse puntano "
            "alla stessa riga storica"
        )

    predictions = {}

    # --------------------------------------------------------
    # INFERENZA EFFETTIVA
    # --------------------------------------------------------

    if chosen:

        X = df.iloc[
            positions
        ][features].copy()

        for col in features:

            X[col] = pd.to_numeric(
                X[col],
                errors="coerce"
            )

        if not np.isfinite(
            X.to_numpy(dtype=float)
        ).all():

            raise RuntimeError(
                "Feature storiche non numeriche "
                "o non finite"
            )

        predicted_labels = model.predict(X)

        probabilities = np.asarray(
            model.predict_proba(X)
        )

        classes = [
            label(c)
            for c in model.classes_
        ]

        if (
            "1" not in classes
            or probabilities.shape != (
                len(chosen),
                len(classes)
            )
        ):

            raise RuntimeError(
                "Dimensioni/classi predittive inattese"
            )

        if (
            not np.isfinite(probabilities).all()
            or not np.allclose(
                probabilities.sum(axis=1),
                1.0,
                atol=1e-6
            )
        ):

            raise RuntimeError(
                "Probabilita non valide"
            )

        positive = classes.index("1")

        # ----------------------------------------------------
        # PREDIZIONE PER RNA
        # ----------------------------------------------------

        for row_i, (
            item,
            df_i,
            method
        ) in enumerate(chosen):

            hist = df.iloc[df_i]

            predictions[item["rna_id"]] = {

                "feature_id": str(
                    hist["feature_id"]
                ),

                "match_method": method,

                "ml_class": label(
                    predicted_labels[row_i]
                ),

                "p_class1": float(
                    probabilities[
                        row_i,
                        positive
                    ]
                ),

                "historical_class": label(
                    hist.get(
                        "antileakage_predicted_label",
                        ""
                    )
                )
            }

    # --------------------------------------------------------
    # COSTRUZIONE RISULTATI
    # --------------------------------------------------------

    rows = []

    for item in payload:

        pred = predictions.get(
            item["rna_id"],
            {}
        )

        rows.append({

            "rna_id_gui": item["rna_id"],

            "qc_status": item["qc_status"],

            "n_new_sites_gui": item["n_new_sites"],

            "feature_id_ml": pred.get(
                "feature_id",
                ""
            ),

            "matching_method": pred.get(
                "match_method",
                ""
            ),

            "class_ml": pred.get(
                "ml_class",
                ""
            ),

            "probability_class1": (
                f"{pred['p_class1']:.12f}"
                if pred
                else ""
            ),

            "historical_label": pred.get(
                "historical_class",
                ""
            ),

            "ml_status": (
                "ML_FEATURE_STORICHE"
                if pred
                else "NON_CLASSIFICABILE"
            )
        })

    status = (
        "PASS_FULL_HISTORICAL_MATCH"
        if len(predictions) == len(payload)
        else "PARTIAL_HISTORICAL_MATCH"
    )

    return rows, unresolved, status


# ============================================================
# J. SCRITTURA JSON SICURA
# ============================================================

def atomic_json(destination, obj):

    temp = destination.with_name(
        destination.name + ".tmp"
    )

    temp.write_text(
        json.dumps(
            obj,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    temp.replace(destination)


# ============================================================
# K. ELABORAZIONE DI UNA CARTELLA GUI
# ============================================================

def process(source, loaded, mode):

    source = source.resolve()

    if (
        not source.is_dir()
        or not source.is_relative_to(
            RUNS.resolve()
        )
    ):

        raise ValueError(
            "Sorgente non valida all'interno di /runs"
        )

    (qpath, qc), (spath, summary) = locate(source)

    signatures = {

        "qc": sha(qpath),

        "summary": sha(spath)
    }

    payload = prepare(qc, summary)

    rows, unresolved, status = classify(
        payload,
        loaded
    )

    if signatures != {

        "qc": sha(qpath),

        "summary": sha(spath)
    }:

        raise RuntimeError(
            "Risultati GUI modificati durante "
            "l'inferenza"
        )

    safe = re.sub(
        r"[^A-Za-z0-9_-]+",
        "_",
        source.name
    )[:90]

    if mode == "manual":

        suffix = "_" + datetime.now(
            timezone.utc
        ).strftime(
            "%Y%m%dT%H%M%S%fZ"
        )

    else:

        suffix = ""

    destination = OUT / (
        mode + "_" + safe + suffix
    )

    if destination.exists():

        raise RuntimeError(
            "Output gia esistente: "
            + str(destination)
        )

    temporary = OUT / (
        ".pending_"
        + safe
        + "_"
        + str(os.getpid())
    )

    temporary.mkdir(
        parents=True,
        exist_ok=False
    )

    try:

        # ----------------------------------------------------
        # RISULTATI TSV
        # ----------------------------------------------------

        tsv = temporary / "09_gui_ml_auto.tsv"

        with tsv.open(
            "w",
            encoding="utf-8",
            newline=""
        ) as f:

            writer = csv.DictWriter(
                f,
                fieldnames=list(rows[0]),
                delimiter="\t"
            )

            writer.writeheader()

            writer.writerows(rows)

        # ----------------------------------------------------
        # MANIFESTO SCIENTIFICO
        # ----------------------------------------------------

        atomic_json(
            temporary / "10_manifest_gui_ml_auto.json",
            {

                "project": "AI PRIMETTA V3.3",

                "step": "22.5",

                "timestamp_utc": utc(),

                "status": status,

                "input_source": str(source),

                "input_sha256": signatures,

                "input_count": len(payload),

                "predictions_count": (
                    len(payload) - len(unresolved)
                ),

                "unresolved": unresolved,

                "sklearn_version": sklearn.__version__,

                "ml_runtime_executed": bool(
                    len(payload) - len(unresolved)
                ),

                "gae_runtime_executed": False,

                "uses_historical_22_features": True,

                "gui_source_modified": False,

                "gui_displays_results": False,

                "scientific_independent_validation": False,

                "training_performed": False,

                "note": (
                    "Inferenza ML su feature e "
                    "embedding storici; "
                    "il GAE non viene rieseguito."
                )
            }
        )

        temporary.rename(destination)

    finally:

        if temporary.exists():

            for path in temporary.iterdir():
                path.unlink()

            temporary.rmdir()

    say(
        "[OK]",
        status,
        f"{len(payload)-len(unresolved)}/{len(payload)} RNA",
        "->",
        destination
    )

    return status


# ============================================================
# L. MONITORAGGIO AUTOMATICO
# ============================================================

def watch(loaded, interval):

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    state_file = OUT / "processed_state.json"

    if state_file.exists():

        state = json.loads(
            state_file.read_text(
                encoding="utf-8"
            )
        )

    else:

        # Il primo avvio non rielabora
        # automaticamente i vecchi esperimenti.

        state = {

            "created_utc": utc(),

            "results": {

                path.name: {
                    "status": "PREESISTENTE"
                }

                for path in RUNS.glob("result_*")
                if path.is_dir()
            }
        }

        atomic_json(
            state_file,
            state
        )

    say(
        "[WATCH] Servizio attivo."
    )

    say(
        "[WATCH] Sorveglianza nuove cartelle result_*"
    )

    observed = {}

    while True:

        for source in sorted(
            RUNS.glob("result_*")
        ):

            if (
                not source.is_dir()
                or source.name in state["results"]
            ):
                continue

            try:

                (qp, _), (sp, _) = locate(source)

                fingerprint = [

                    (
                        path.stat().st_size,
                        path.stat().st_mtime_ns
                    )

                    for path in (qp, sp)
                ]

                age = time.time() - max(
                    qp.stat().st_mtime,
                    sp.stat().st_mtime
                )

                if (
                    age < interval
                    or observed.get(
                        source.name
                    ) != fingerprint
                ):

                    observed[source.name] = fingerprint

                    continue

                status = process(
                    source,
                    loaded,
                    "auto"
                )

                state["results"][source.name] = {

                    "status": status,

                    "processed_utc": utc()
                }

                atomic_json(
                    state_file,
                    state
                )

            except (ValueError, OSError):

                # La GUI potrebbe essere
                # ancora impegnata a scrivere
                # i risultati. Riprova.

                continue

            except Exception as exc:

                say(
                    "[ERRORE]",
                    source.name,
                    type(exc).__name__,
                    str(exc)
                )

                state["results"][source.name] = {

                    "status": "ERRORE",

                    "message": str(exc),

                    "processed_utc": utc()
                }

                atomic_json(
                    state_file,
                    state
                )

        time.sleep(interval)


# ============================================================
# M. AVVIO DEL BACKEND
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    mode = parser.add_mutually_exclusive_group(
        required=True
    )

    mode.add_argument(
        "--once",
        metavar="CARTELLA"
    )

    mode.add_argument(
        "--watch",
        action="store_true"
    )

    parser.add_argument(
        "--interval",
        type=int,
        default=5
    )

    args = parser.parse_args()

    if not 2 <= args.interval <= 300:

        parser.error(
            "--interval tra 2 e 300 secondi"
        )

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    loaded = load_model()

    if args.once:

        status = process(
            Path(args.once),
            loaded,
            "manual"
        )

        if status != "PASS_FULL_HISTORICAL_MATCH":
            sys.exit(2)

    else:

        watch(
            loaded,
            args.interval
        )


if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        say(
            "[STOP] Monitor arrestato"
        )

    except Exception as exc:

        say(
            "[ERRORE]",
            type(exc).__name__,
            str(exc)
        )

        sys.exit(1)
