#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 32

Verifica in sola lettura:
- Riferimento genomico B73 NAM-5.0
- Sequenze ncRNA B73
- Annotazione genica GTF
- Annotazione dei trasposoni GFF3
- FASTQ originali PRJNA419326

--audit: verifica campionaria senza modifiche.

Non scarica dati, non esegue training,
non modifica GUI, Docker o modelli.

La verifica iniziale dei gzip non equivale
al controllo di integrita' dei file completi.
"""

import argparse
import gzip

from collections import Counter
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE DEL PROGETTO
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

GENOME = (
    ROOT
    / "target_prediction_v33_vscode"
    / "references"
    / "Zm-B73-REFERENCE-NAM-5.0.fa.gz"
)

NCRNA = (
    ROOT
    / "data/reference/zea_mays_ensembl_plants"
    / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.ncrna.fa.gz"
)

GTF = (
    ROOT
    / "data/reference/zea_mays_ensembl_plants"
    / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.63.gtf.gz"
)

TE_GFF = (
    ROOT
    / "target_prediction_v33_vscode"
    / "references"
    / "Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz"
)

FASTQ_DIR = (
    ROOT
    / "data/public_ena/PRJNA419326/fastq/mirnaseq"
)


# ============================================================
# 2. VERIFICA FILE GZIP
# ============================================================

def check_gzip(path):

    if not path.is_file():

        raise FileNotFoundError(
            f"File non trovato: {path}"
        )

    with path.open("rb") as f:

        magic = f.read(2)

        if magic != b"\x1f\x8b":

            raise ValueError(
                f"Formato gzip non valido: {path}"
            )

    return round(
        path.stat().st_size / (1024 * 1024),
        1
    )


# ============================================================
# 3. LETTURA CAMPIONARIA FASTA
# ============================================================

def fasta_preview(path):

    size = check_gzip(path)

    headers = []

    first_sequence = None

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as f:

        for index, line in enumerate(f):

            if index >= 1500:
                break

            if len(headers) >= 3:
                break

            line = line.strip()

            if line.startswith(">"):

                headers.append(
                    line[1:].split()[0]
                )

            elif (
                line
                and first_sequence is None
                and headers
            ):

                first_sequence = (
                    line.upper()[:70]
                )

    if (
        not headers
        or first_sequence is None
    ):

        raise ValueError(
            f"FASTA non riconosciuto: {path}"
        )

    return (
        size,
        headers,
        first_sequence
    )


# ============================================================
# 4. LETTURA CAMPIONARIA GTF / GFF3
# ============================================================

def annotation_preview(
    path,
    max_lines=25000
):

    size = check_gzip(path)

    feature_types = Counter()

    first_record = None

    first_gene_attributes = None

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as f:

        for index, line in enumerate(f):

            if index >= max_lines:
                break

            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue

            columns = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(columns) != 9:
                continue

            feature = columns[2]

            feature_types[feature] += 1

            if first_record is None:

                first_record = (
                    columns[0],
                    feature,
                    columns[3],
                    columns[4]
                )

            if (
                feature == "gene"
                and first_gene_attributes is None
            ):

                first_gene_attributes = (
                    columns[8][:120]
                )

            if (
                first_gene_attributes is not None
                and feature_types.get(
                    "transcript", 0
                ) > 0
            ):
                break

    if first_record is None:

        raise ValueError(
            "Formato GTF/GFF3 non riconosciuto: "
            + str(path)
        )

    return (
        size,
        first_record,
        feature_types,
        first_gene_attributes
    )


# ============================================================
# 5. VERIFICA PRIMA LETTURA FASTQ
# ============================================================

def fastq_preview(path):

    check_gzip(path)

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace"
    ) as f:

        record = [

            f.readline().rstrip("\n\r")

            for _ in range(4)
        ]

    if (
        any(not item for item in record)
        or not record[0].startswith("@")
        or not record[2].startswith("+")
        or len(record[1]) != len(record[3])
        or not set(
            record[1].upper()
        ) <= set("ACGTUN")
    ):

        raise ValueError(
            f"Prima lettura FASTQ non valida: {path}"
        )

    return len(
        record[1]
    )


# ============================================================
# 6. PROGRAMMA PRINCIPALE
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

    print("=" * 68)

    print(
        "AI PRIMETTA V3.3 | STEP 32"
    )

    print(
        "AUDIT RIFERIMENTI BIOLOGICI B73"
    )

    print("=" * 68)

    # ========================================================
    # 7. VERIFICA GENOMA E ncRNA
    # ========================================================

    for label, path in (

        ("GENOMA_CANDIDATO", GENOME),

        ("NC_RNA_NON_GENOMA", NCRNA)

    ):

        size, headers, first_seq = (
            fasta_preview(path)
        )

        print(
            f"[{label}] "
            f"{path.relative_to(ROOT)}"
        )

        print(
            f"  Dimensione compressa: "
            f"{size} MiB"
        )

        print(
            "  Intestazioni campionate:",
            ", ".join(headers)
        )

        print(
            "  Prime basi:",
            first_seq[:45]
        )

    # ========================================================
    # 8. VERIFICA ANNOTAZIONI
    # ========================================================

    for label, path in (

        ("ANNOTAZIONE_GTF", GTF),

        ("ANNOTAZIONE_TE", TE_GFF)

    ):

        (
            size,
            example,
            features,
            attributes
        ) = annotation_preview(path)

        print(
            f"[{label}] "
            f"{path.relative_to(ROOT)}"
        )

        print(
            f"  Dimensione compressa: "
            f"{size} MiB"
        )

        print(
            f"  Primo record: "
            f"seq={example[0]} "
            f"tipo={example[1]} "
            f"intervallo={example[2]}-{example[3]}"
        )

        print(
            "  Tipi nel campione:",
            dict(
                features.most_common(6)
            )
        )

        if attributes:

            print(
                "  Attributi gene:",
                attributes
            )

    # ========================================================
    # 9. VERIFICA FASTQ PRJNA419326
    # ========================================================

    if not FASTQ_DIR.is_dir():

        raise FileNotFoundError(
            "Cartella FASTQ assente: "
            + str(FASTQ_DIR)
        )

    fastqs = sorted(
        FASTQ_DIR.glob("*.fastq.gz")
    )

    print(
        "[FASTQ] File nella cartella:",
        len(fastqs)
    )

    if not fastqs:

        raise ValueError(
            "Nessun FASTQ trovato"
        )

    lengths = Counter()

    for path in fastqs:

        read_length = fastq_preview(path)

        lengths[
            read_length
        ] += 1

    print(
        "[FASTQ] Prime letture formalmente valide:",
        f"{sum(lengths.values())}/{len(fastqs)}"
    )

    print(
        "[FASTQ] Lunghezze delle prime letture:",
        dict(
            sorted(lengths.items())
        )
    )

    print(
        "[FASTQ] Primi file:",
        ", ".join(
            p.name
            for p in fastqs[:5]
        )
    )

    # ========================================================
    # 10. CONCLUSIONI DELL'AUDIT
    # ========================================================

    print(
        "[LIMITE] Verificate le intestazioni "
        "e le prime letture, non l'integrita' "
        "completa dei file compressi"
    )

    print(
        "[LIMITE] La presenza dei file "
        "non dimostra che le sequenze "
        "siano miRNA o siRNA"
    )

    print(
        "[DECISIONE] Dopo la verifica "
        "del riferimento, individuare "
        "i loci delle 2 small RNA "
        "sul DNA B73 NAM-5.0"
    )

    print(
        "[OK] Audit in sola lettura; "
        "GUI, Docker, modelli e dati invariati"
    )


# ============================================================
# 11. AVVIO
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            f"[ERRORE] "
            f"{type(exc).__name__}: {exc}"
        )

        raise SystemExit(1)