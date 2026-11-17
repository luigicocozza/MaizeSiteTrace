
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 38 CORRETTO

Audit in sola lettura:
- Matrici count e RPM
- Identita' delle due small RNA
- Variabilita' delle repliche
- FASTQ originali e trimmed
- Disegno sperimentale PRJNA419326
- Corrispondenza fra run_accession e colonne RPM

CORREZIONE:
Gli identificativi SRR vengono letti
ESCLUSIVAMENTE dalla colonna run_accession.
Non vengono cercati nei percorsi FTP.

Nessun training, nessun download.
GUI, Docker, dati e modelli invariati.
"""

import argparse
import csv
import math
import re
import statistics
import sys

from collections import Counter
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

DATA = ROOT / "outputs" / "prjna419326"

COUNT_FILE = (
    DATA /
    "PRJNA419326_mirnaseq_srna_count_matrix_filtered.tsv"
)

RPM_FILE = (
    DATA /
    "PRJNA419326_mirnaseq_srna_rpm_matrix_filtered.tsv"
)

DESIGN_FILE = (
    DATA /
    "PRJNA419326_experimental_design_curated.tsv"
)

RAW_FOLDER = (
    ROOT /
    "data/public_ena/PRJNA419326/fastq/mirnaseq"
)

TRIM_FOLDER = (
    ROOT /
    "data/public_ena/PRJNA419326/"
    "processed/mirnaseq_trimmed"
)

EXPECTED = {

    "sRNA_00000037":
        "UUCUUGGAUUUAUGAAAGACGAACA",

    "sRNA_00000435":
        "AGAUGAUACGAGCACCUCUGAGCU"
}

SAMPLES = [

    f"{genotype}_{condition}_sRNA_rep{rep}"

    for genotype in ("B73", "H21")

    for condition in ("DS", "WW")

    for rep in (1, 2, 3)
]

RUN_RE = re.compile(
    r"SRR\d{6,}$"
)

REP_RE = re.compile(
    r"(?:rep(?:licate)?[_\- ]*)?([123])$",
    re.IGNORECASE
)


# ============================================================
# 2. CONTROLLO CONDIZIONI
# ============================================================

def require(ok, message):

    if not ok:
        raise ValueError(message)


# ============================================================
# 3. LETTURA MIRATA DELLE MATRICI
# ============================================================

def table_extract(path):

    require(
        path.is_file(),
        f"Matrice assente: {path}"
    )

    selected = {}

    nrows = 0

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        columns = reader.fieldnames or []

        required = {

            "feature_id",
            "sequence",
            "length",
            "total_count",
            "samples_detected",
            *SAMPLES
        }

        require(
            required.issubset(columns)
            and len(columns) == len(set(columns)),
            f"Intestazioni incomplete o duplicate: {path.name}"
        )

        for row in reader:

            nrows += 1

            sid = (
                row.get("feature_id") or ""
            ).strip()

            if sid not in EXPECTED:
                continue

            require(
                sid not in selected,
                f"Feature duplicata: {sid}"
            )

            sequence = (
                (row.get("sequence") or "")
                .strip()
                .upper()
                .replace("T", "U")
            )

            require(
                sequence == EXPECTED[sid]
                and int(row["length"]) == len(sequence),
                f"Sequenza o lunghezza discordante: {sid}"
            )

            require(
                int(row["samples_detected"]) == 12,
                f"Rilevazione diversa da 12: {sid}"
            )

            values = {}

            for sample in SAMPLES:

                value = float(
                    row[sample]
                )

                require(
                    math.isfinite(value)
                    and value >= 0,
                    f"Valore non valido: {sid}/{sample}"
                )

                values[sample] = value

            selected[sid] = {

                "values":
                    values,

                "total_count":
                    int(row["total_count"])
            }

    require(
        set(selected) == set(EXPECTED),
        "Mancano una o entrambe le small RNA"
    )

    return nrows, selected


# ============================================================
# 4. IDENTIFICAZIONE DEI FASTQ
# ============================================================

def discover_runs(folder, suffix):

    require(
        folder.is_dir(),
        f"Cartella assente: {folder}"
    )

    discovered = {}

    for file in sorted(
        folder.glob(f"SRR*{suffix}")
    ):

        name = file.name

        run = name[:-len(suffix)]

        require(
            RUN_RE.fullmatch(run) is not None,
            f"Nome FASTQ inatteso: {name}"
        )

        require(
            run not in discovered,
            f"Run duplicato nella cartella: {run}"
        )

        discovered[run] = file

    return discovered


# ============================================================
# 5. RICOSTRUZIONE DEL CAMPIONE SCIENTIFICO
# ============================================================

def canonical_sample(row):

    """
    Ricostruisce la colonna della matrice
    a partire dai campi espliciti:

    genotype
    condition_code
    replicate

    Non utilizza mai fastq_ftp.
    """

    genotype = (
        row.get("genotype") or ""
    ).strip().upper()

    condition = (
        row.get("condition_code") or ""
    ).strip().upper()

    replicate = (
        row.get("replicate") or ""
    ).strip()

    match = REP_RE.fullmatch(
        replicate
    )

    if (
        genotype in ("B73", "H21")
        and condition in ("DS", "WW")
        and match
    ):

        return (
            f"{genotype}_{condition}"
            f"_sRNA_rep{match.group(1)}"
        )

    alias = (
        row.get("sample_alias") or ""
    ).strip()

    if alias in SAMPLES:
        return alias

    return None


# ============================================================
# 6. AUDIT DEL DISEGNO SPERIMENTALE
# ============================================================

def audit_design(raw_runs, trimmed_runs):

    require(
        DESIGN_FILE.is_file(),
        f"Design assente: {DESIGN_FILE}"
    )

    with DESIGN_FILE.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        columns = reader.fieldnames or []

        require(
            "run_accession" in columns,
            "Colonna run_accession assente"
        )

        require(
            len(columns) == len(set(columns)),
            "Colonne design duplicate"
        )

        rows = list(reader)

    require(
        bool(rows),
        "Disegno sperimentale vuoto"
    )

    print(
        f"[DESIGN] righe={len(rows)} "
        f"| colonne={', '.join(columns[:18])}"
    )

    # ========================================================
    # 7. CORREZIONE DEL BUG SRR632
    # ========================================================

    index = {}

    for row in rows:

        # SI LEGGE SOLTANTO run_accession.
        # I campi fastq_ftp, fastq_md5 e
        # experiment_title NON vengono
        # utilizzati per ricavare gli SRR.

        run = (
            row.get("run_accession") or ""
        ).strip()

        require(
            RUN_RE.fullmatch(run) is not None,
            f"run_accession non valido: {run!r}"
        )

        require(
            run not in index,
            f"run_accession realmente duplicato: {run}"
        )

        index[run] = row

    print(
        "[PASS] run_accession completi "
        f"e univoci nel design: {len(index)}"
    )

    # ========================================================
    # 8. VERIFICA DEI 12 FASTQ ORIGINALI
    # ========================================================

    missing = sorted(
        set(raw_runs) - set(index)
    )

    require(
        not missing,
        "FASTQ small RNA non presenti "
        f"nel design: {missing}"
    )

    require(
        set(raw_runs) == set(trimmed_runs),
        "I FASTQ raw e trimmed "
        "hanno identificativi SRR differenti"
    )

    print(
        f"[PASS] I {len(raw_runs)} SRR "
        "small RNA corrispondono "
        "ai FASTQ raw e trimmed"
    )

    # ========================================================
    # 9. TIPI DI SEQUENZIAMENTO
    # ========================================================

    categories = Counter(

        (
            row.get("sequencing_type")
            or "non_specificato"
        ).strip()

        for row in rows
    )

    print(
        "[DESIGN] sequencing_type nel file:",
        dict(categories)
    )

    # ========================================================
    # 10. MAPPATURA SRR -> COLONNE DELLA MATRICE
    # ========================================================

    sample_index = {}

    unresolved = 0

    for run in sorted(raw_runs):

        row = index[run]

        sample = canonical_sample(
            row
        )

        print(

            f"[RUN] {run} "

            f"| genotype={row.get('genotype', '')} "

            f"| condition_code="
            f"{row.get('condition_code', '')} "

            f"| replicate="
            f"{row.get('replicate', '')} "

            f"| sequencing_type="
            f"{row.get('sequencing_type', '')} "

            f"| matrice="
            f"{sample or 'DA_VERIFICARE'}"
        )

        if sample is None:

            unresolved += 1

        else:

            require(
                sample not in sample_index,
                "Campione di matrice associato "
                f"a piu' SRR: {sample}"
            )

            sample_index[sample] = run

    # ========================================================
    # 11. VERIFICA COMPLETA DELLE 12 CONDIZIONI
    # ========================================================

    if (
        unresolved == 0
        and set(sample_index) == set(SAMPLES)
    ):

        print(
            "[PASS] Mappatura completa: "
            "12 SRR -> 12 colonne di matrice "
            "B73/H21 x DS/WW x 3"
        )

    else:

        print(

            "[ATTENZIONE] Associazione "
            "non conclusiva: "

            f"{unresolved} run non interpretabili; "

            "colonne mancanti="

            + str(
                sorted(
                    set(SAMPLES)
                    - set(sample_index)
                )
            )
        )

    print(

        "[DESIGN] Altri run del progetto "
        "non appartenenti ai 12 FASTQ "
        f"small RNA: {len(index) - len(raw_runs)}"
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
        required=True,
        action="store_true"
    )

    parser.parse_args()

    print("=" * 70)

    print(
        "AI PRIMETTA V3.3 | STEP 38 "
        "- AUDIT REPLICHE E DESIGN CORRETTO"
    )

    print("=" * 70)

    # ========================================================
    # 13. CONTROLLO MATRICI COUNT E RPM
    # ========================================================

    ncount, counts = table_extract(
        COUNT_FILE
    )

    nrpm, rpm = table_extract(
        RPM_FILE
    )

    require(
        ncount == nrpm,
        "Matrici count/RPM con "
        "numeri di righe differenti"
    )

    print(
        f"[PASS] Matrici count/RPM: "
        f"{ncount:,} feature ciascuna"
    )

    # ========================================================
    # 14. VERIFICA DELLE DUE SMALL RNA
    # ========================================================

    for sid in EXPECTED:

        samples = counts[sid]["values"]

        require(
            all(
                value.is_integer()
                for value in samples.values()
            ),
            f"Conteggi non interi: {sid}"
        )

        require(
            int(sum(samples.values()))
            == counts[sid]["total_count"],
            "total_count non corrisponde "
            f"alla somma: {sid}"
        )

        require(
            rpm[sid]["total_count"]
            == counts[sid]["total_count"],
            "total_count count/RPM "
            f"discordante: {sid}"
        )

        for sample, count in samples.items():

            require(
                (count == 0)
                ==
                (
                    rpm[sid]["values"][sample]
                    == 0
                ),
                "Zero/nonzero count/RPM "
                f"incoerente: {sid}/{sample}"
            )

        print(

            f"[PASS] {sid} "

            f"| sequenza={EXPECTED[sid]} "

            f"| conteggi="
            f"{counts[sid]['total_count']:,} "

            "| campioni=12/12"
        )

        # ====================================================
        # 15. REPLICHE DELLE QUATTRO CONDIZIONI
        # ====================================================

        for genotype in ("B73", "H21"):

            medians = {}

            for condition in ("WW", "DS"):

                arr = [

                    rpm[sid]["values"][

                        f"{genotype}_{condition}"
                        f"_sRNA_rep{i}"
                    ]

                    for i in (1, 2, 3)
                ]

                medians[condition] = (
                    statistics.median(arr)
                )

                if min(arr) > 0:

                    ratio = (
                        f"{max(arr)/min(arr):.1f}"
                    )

                elif max(arr) > 0:

                    ratio = "infinito"

                else:

                    ratio = "n.d."

                print(

                    f"[REPLICHE] {sid} "

                    f"{genotype}_{condition} "

                    "| RPM="

                    + ", ".join(
                        f"{v:.3f}"
                        for v in arr
                    )

                    + f" | mediana="
                    f"{medians[condition]:.3f}"

                    + f" | max/min={ratio}"
                )

                # ============================================
                # 16. SEGNALAZIONE ALTA VARIABILITA'
                # ============================================

                if (
                    min(arr) == 0

                    or

                    (
                        min(arr) > 0
                        and max(arr)/min(arr) >= 10
                    )
                ):

                    print(

                        f"[VARIABILITA] {sid} "

                        f"{genotype}_{condition}: "

                        "rapporto >=10 o zero"
                    )

            # ================================================
            # 17. RAPPORTO DESCRITTIVO DS / WW
            # ================================================

            if medians["WW"] > 0:

                ratio = (
                    medians["DS"]
                    / medians["WW"]
                )

                print(

                    f"[DESCRITTIVO] {sid} "

                    f"{genotype}: "

                    f"mediane DS/WW={ratio:.3f}"
                )

    # ========================================================
    # 18. IDENTIFICAZIONE FASTQ
    # ========================================================

    raw = discover_runs(
        RAW_FOLDER,
        ".fastq.gz"
    )

    trimmed = discover_runs(
        TRIM_FOLDER,
        ".trimmed.fastq.gz"
    )

    print(
        f"[FASTQ] raw={len(raw)} "
        f"| trimmed={len(trimmed)}"
    )

    require(
        len(raw) == len(trimmed) == 12,
        "Attese 12 librerie raw e 12 trimmed"
    )

    # ========================================================
    # 19. CONTROLLO DEI METADATI CORRETTO
    # ========================================================

    audit_design(
        raw,
        trimmed
    )

    # ========================================================
    # 20. CONCLUSIONI
    # ========================================================

    print(

        "[LIMITE] Mediane RPM descrittive; "
        "nessun test differenziale "
        "e nessuna classificazione miRNA/siRNA."
    )

    print(

        "[OK] Audit in sola lettura; "
        "GUI, Docker, dati e modelli invariati"
    )


# ============================================================
# 21. AVVIO
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

        raise SystemExit(1)
