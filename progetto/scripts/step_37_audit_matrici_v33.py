#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 37

Confronto dei risultati FASTQ grezzi
con le matrici small RNA gia' elaborate.

--audit:
    Solo lettura.
    Verifica la provenienza dello STEP36.
    Individua conteggi elaborati e metadati.

Nessun download, training o modifica
a GUI, Docker, dati e modelli.
"""

import argparse
import csv
import hashlib
import json
import sys

from collections import Counter
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

SRNAS = (
    "sRNA_00000037",
    "sRNA_00000435"
)

MATRICES = (

    "outputs/prjna419326/"
    "PRJNA419326_mirnaseq_srna_count_matrix_filtered.tsv",

    "outputs/prjna419326/"
    "PRJNA419326_mirnaseq_srna_rpm_matrix_filtered.tsv",

    "delivery/"
    "maize_ncrna_gnn_v3_prjna419326_docker_gui_ml_consegna_20261011/"
    "outputs/prjna419326/"
    "PRJNA419326_mirnaseq_srna_count_matrix_filtered.tsv"
)


# ============================================================
# 2. CHECKSUM
# ============================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1048576),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# ============================================================
# 3. CONTROLLO FILE
# ============================================================

def safe_file(value):

    p = Path(value)

    if not p.is_absolute():
        p = ROOT / p

    p = p.resolve()

    if (
        not p.is_file()
        or not p.is_relative_to(
            ROOT.resolve()
        )
    ):

        raise ValueError(
            f"File assente o esterno al progetto: {p}"
        )

    return p


# ============================================================
# 4. LETTURA TSV
# ============================================================

def read_tsv(path):

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        columns = reader.fieldnames or []

        if (
            not columns
            or len(columns) != len(set(columns))
        ):

            raise ValueError(
                f"Intestazioni TSV invalide: {path}"
            )

        rows = list(reader)

    if not rows:

        raise ValueError(
            f"Tabella vuota: {path}"
        )

    return columns, rows


# ============================================================
# 5. RECUPERO E VERIFICA STEP36
# ============================================================

def load_step36():

    reports = sorted(
        RUNS.glob(
            "STEP36_SUPPORTO_FASTQ_*/"
            "RAPPORTO_STEP36.json"
        )
    )

    if not reports:

        raise FileNotFoundError(
            "Report STEP36 non trovato"
        )

    report_path = safe_file(
        reports[-1]
    )

    report = json.loads(
        report_path.read_text(
            encoding="utf-8"
        )
    )

    if (
        report.get("step") != 36
        or report.get("sample_count") != 12
    ):

        raise ValueError(
            "Manifest STEP36 inatteso"
        )

    # ========================================================
    # 6. CONTROLLO CHECKSUM
    # ========================================================

    table_path = safe_file(
        report["output_tsv"]
    )

    if (
        sha256(table_path)
        != report["output_tsv_sha256"]
    ):

        raise ValueError(
            "Checksum STEP36 non coerente"
        )

    columns, rows = read_tsv(
        table_path
    )

    required = {

        "sample",
        "srna_id",
        "reads_examined",
        "sense_prefix",
        "sense_anywhere",
        "reverse_prefix",
        "reverse_anywhere",
        "reads_any_orientation"
    }

    if (
        not required.issubset(columns)
        or len(rows) != 24
    ):

        raise ValueError(
            "Attese 24 righe STEP36 "
            "con le colonne previste"
        )

    # ========================================================
    # 7. RICALCOLO DEI TOTALI
    # ========================================================

    totals = {
        sid: Counter()
        for sid in SRNAS
    }

    samples = {}
    seen = set()

    for row in rows:

        sid = row["srna_id"]
        sample = row["sample"]

        if (
            sid not in SRNAS
            or (sample, sid) in seen
        ):

            raise ValueError(
                "Coppia campione/sRNA inattesa "
                "o duplicata"
            )

        seen.add(
            (sample, sid)
        )

        values = {}

        for key in required - {
            "sample",
            "srna_id"
        }:

            n = int(
                row[key]
            )

            if n < 0:

                raise ValueError(
                    "Conteggio negativo"
                )

            values[key] = n

            totals[sid][key] += n

        if not (

            values["sense_prefix"]
            <= values["sense_anywhere"]
            <= values["reads_any_orientation"]

            and

            values["reverse_prefix"]
            <= values["reverse_anywhere"]
            <= values["reads_any_orientation"]
            <= values["reads_examined"]
        ):

            raise ValueError(
                "Conteggi FASTQ incoerenti: "
                + sample
            )

        if (
            sample in samples
            and samples[sample]
            != values["reads_examined"]
        ):

            raise ValueError(
                "Totali differenti "
                "per lo stesso campione"
            )

        samples[sample] = (
            values["reads_examined"]
        )

    if len(samples) != 12:

        raise ValueError(
            "Attesi 12 campioni distinti"
        )

    # ========================================================
    # 8. CONFRONTO CON IL MANIFEST STEP36
    # ========================================================

    for sid, counts in totals.items():

        previous = (
            report["summary"]
            .get(sid, {})
        )

        for key, value in counts.items():

            if previous.get(key) != value:

                raise ValueError(
                    f"Totale non riproducibile: "
                    f"{sid}/{key}"
                )

    return (
        report_path,
        totals,
        samples
    )


# ============================================================
# 9. ESAME DELLE MATRICI ELABORATE
# ============================================================

def inspect_matrix(path):

    columns, rows = read_tsv(
        path
    )

    selected = {
        sid: []
        for sid in SRNAS
    }

    for row in rows:

        for sid in SRNAS:

            if sid in row.values():

                selected[sid].append(
                    row
                )

    print(
        "[MATRICE]",
        path.relative_to(ROOT),
        f"| righe={len(rows)}",
        f"| colonne={len(columns)}"
    )

    print(
        "[INTESTAZIONI]",
        ", ".join(
            columns[:24]
        )
    )

    # ========================================================
    # 10. VALORI DELLE DUE SMALL RNA
    # ========================================================

    for sid, found in selected.items():

        if len(found) != 1:

            print(
                f"[MATRICE_sRNA] {sid}: "
                f"righe_trovate={len(found)}"
            )

            continue

        record = found[0]

        numeric = []

        for column, value in record.items():

            if (
                value is None
                or value.strip() in (
                    "",
                    "NA",
                    "NaN"
                )
                or value.strip() in SRNAS
            ):

                continue

            try:

                number = float(
                    value
                )

                if (
                    number >= 0
                    and number != float("inf")
                ):

                    numeric.append(
                        (column, value)
                    )

            except (
                ValueError,
                TypeError
            ):

                continue

        print(
            f"[MATRICE_sRNA] {sid}: "
            f"colonne_numeriche={len(numeric)}"
        )

        for column, value in numeric[:18]:

            print(
                f"  {column}: {value}"
            )

        if len(numeric) > 18:

            print(
                "  ... altre",
                len(numeric) - 18,
                "colonne"
            )


# ============================================================
# 11. RICERCA DEI METADATI
# ============================================================

def find_metadata():

    keywords = (

        "sample",
        "metadata",
        "design",
        "manifest",
        "adapter",
        "trim",
        "cutadapt",
        "fastp",
        "sra_run",
        "runinfo",
        "biosample"
    )

    folders = (

        ROOT / "outputs/prjna419326",

        ROOT / "data/public_ena/PRJNA419326",

        ROOT / "target_prediction_v33_vscode/references"
    )

    found = set()

    for folder in folders:

        if not folder.is_dir():
            continue

        for path in folder.rglob("*"):

            if not path.is_file():
                continue

            if any(

                word in path.name.lower()

                for word in keywords
            ):

                found.add(
                    str(
                        path.relative_to(ROOT)
                    )
                )

    print(
        "[METADATI_CANDIDATI]",
        len(found),
        "file individuati per nome"
    )

    for path in sorted(found)[:30]:

        print(
            " ",
            path
        )

    if len(found) > 30:

        print(
            " ... altri",
            len(found) - 30,
            "file"
        )


# ============================================================
# 12. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--audit",
        action="store_true",
        required=True
    )

    parser.parse_args()

    print("=" * 70)

    print(
        "AI PRIMETTA V3.3 | STEP37 "
        "- AUDIT MATRICI E CAMPIONI"
    )

    print("=" * 70)

    # ========================================================
    # 13. VERIFICA STEP36
    # ========================================================

    report, totals, samples = (
        load_step36()
    )

    print(
        "[PASS] STEP36 e SHA256 verificati; "
        "12 SRR, 24 righe"
    )

    # ========================================================
    # 14. SINTESI DELLE LETTURE GREZZE
    # ========================================================

    for sid, counts in totals.items():

        n = counts[
            "reads_examined"
        ]

        matches = counts[
            "reads_any_orientation"
        ]

        prefix = counts[
            "sense_prefix"
        ]

        per_million = (
            1000000 * matches / n
        )

        print(

            f"[FASTQ] {sid}: "

            f"reads={n:,}, "

            f"match={matches:,}, "

            f"prefix={prefix:,}, "

            f"raw_match_per_million="
            f"{per_million:.2f}"
        )

    print(
        "[CAMPIONI]",
        ", ".join(
            sorted(samples)
        )
    )

    # ========================================================
    # 15. LETTURA MATRICI GIA' ELABORATE
    # ========================================================

    inspected = set()

    for relative in MATRICES:

        path = (
            ROOT / relative
        )

        if (
            path.is_file()
            and path.resolve() not in inspected
        ):

            inspect_matrix(
                path
            )

            inspected.add(
                path.resolve()
            )

    if not inspected:

        print(
            "[AVVISO] Nessuna matrice "
            "elaborata nota trovata"
        )

    # ========================================================
    # 16. METADATI DISPONIBILI
    # ========================================================

    find_metadata()

    print(
        "[LIMITE] I match nei FASTQ grezzi "
        "non equivalgono ai conteggi "
        "di small RNA mature."
    )

    print(
        "[LIMITE] Le condizioni sperimentali "
        "non sono assegnate senza "
        "metadati verificati."
    )

    print(
        "[OK] Audit in sola lettura; "
        "GUI, Docker, dati e modelli invariati"
    )


# ============================================================
# 17. AVVIO
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            f"[ERRORE] "
            f"{type(exc).__name__}: {exc}",
            file=sys.stderr
        )

        sys.exit(1)