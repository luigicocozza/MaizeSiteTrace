#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 35 — CONTESTO GENOMICO E ncRNA

--audit: analisi in sola lettura
--run: analisi e salvataggio rapporto JSON

Non modifica dati originali, GUI, Docker o modelli.
"""

import argparse
import csv
import gzip
import hashlib
import json
import sys

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


# =====================================================
# 1. PERCORSI
# =====================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

NCRNA = (
    ROOT
    / "data/reference/zea_mays_ensembl_plants"
    / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.ncrna.fa.gz"
)

EXPECTED = {
    "sRNA_00000037": 587,
    "sRNA_00000435": 2
}


# =====================================================
# 2. CONTROLLO FILE
# =====================================================

def checked_file(value):

    p = Path(value)

    if not p.is_absolute():
        p = ROOT / p

    p = p.resolve()

    if (
        not p.is_file()
        or not p.is_relative_to(RUNS.resolve())
    ):

        raise ValueError(
            f"File mancante o fuori da runs: {p}"
        )

    return p


# =====================================================
# 3. SHA256
# =====================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1048576),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# =====================================================
# 4. LETTURA FASTA
# =====================================================

def fasta_records(stream):

    name = None
    parts = []

    for raw in stream:

        line = raw.strip()

        if not line:
            continue

        if line.startswith(">"):

            if name is not None:

                if not parts:
                    raise ValueError(
                        "Record FASTA vuoto"
                    )

                yield (
                    name,
                    "".join(parts)
                    .upper()
                    .replace("T", "U")
                )

            name = line[1:]
            parts = []

        else:

            if name is None:
                raise ValueError(
                    "FASTA senza intestazione"
                )

            parts.append(line)

    if name is not None:

        if not parts:
            raise ValueError(
                "Ultimo record FASTA vuoto"
            )

        yield (
            name,
            "".join(parts)
            .upper()
            .replace("T", "U")
        )


# =====================================================
# 5. COMPLEMENTO INVERSO RNA
# =====================================================

def reverse_complement(seq):

    translation = str.maketrans(
        "ACGU",
        "UGCA"
    )

    return seq.translate(
        translation
    )[::-1]


# =====================================================
# 6. RECUPERO STEP34
# =====================================================

def read_inputs():

    reports = sorted(
        RUNS.glob(
            "STEP34_ANNOTAZIONE_LOCI_*/"
            "RAPPORTO_STEP34.json"
        )
    )

    if not reports:
        raise FileNotFoundError(
            "STEP34 non trovato"
        )

    p34 = checked_file(
        reports[-1]
    )

    r34 = json.loads(
        p34.read_text(
            encoding="utf-8"
        )
    )

    if (
        str(r34.get("step")) != "34"
        or r34.get("total_loci") != 589
    ):

        raise ValueError(
            "Manifest STEP34 non coerente"
        )

    # =================================================
    # 7. TABELLA DEI 589 LOCI
    # =================================================

    p34_loci = checked_file(
        r34["output_tsv"]
    )

    if (
        sha256(p34_loci)
        != r34["output_sha256"]
    ):

        raise ValueError(
            "Checksum STEP34 non corrispondente"
        )

    # =================================================
    # 8. CATENA STEP33
    # =================================================

    p33 = checked_file(
        r34["source_step33"]
    )

    if (
        sha256(p33)
        != r34["source_step33_sha256"]
    ):

        raise ValueError(
            "Checksum STEP33 non corrispondente"
        )

    r33 = json.loads(
        p33.read_text(
            encoding="utf-8"
        )
    )

    # =================================================
    # 9. SEQUENZE ORIGINALI SMALL RNA
    # =================================================

    fasta_path = checked_file(
        r33["srna_fasta"]
    )

    if (
        sha256(fasta_path)
        != r33["srna_fasta_sha256"]
    ):

        raise ValueError(
            "Checksum FASTA small RNA differente"
        )

    with fasta_path.open(
        encoding="utf-8-sig"
    ) as f:

        queries = {

            header.split()[0]: seq

            for header, seq in fasta_records(f)
        }

    if set(queries) != set(EXPECTED):

        raise ValueError(
            "Small RNA originali inattese"
        )

    # =================================================
    # 10. LETTURA DEI LOCI
    # =================================================

    with p34_loci.open(
        encoding="utf-8-sig",
        newline=""
    ) as f:

        rows = list(
            csv.DictReader(
                f,
                delimiter="\t"
            )
        )

    counts = Counter(
        r["srna_id"]
        for r in rows
    )

    if (
        len(rows) != 589
        or counts != Counter(EXPECTED)
    ):

        raise ValueError(
            "Numero di loci inatteso"
        )

    return (
        p34,
        p34_loci,
        queries,
        rows
    )


# =====================================================
# 11. CONFRONTO CON IL RIFERIMENTO ncRNA
# =====================================================

def inspect_ncrna(queries):

    if not NCRNA.is_file():

        raise FileNotFoundError(
            f"Riferimento ncRNA assente: {NCRNA}"
        )

    hits = {

        sid: {

            "sense_records": 0,

            "reverse_complement_records": 0,

            "sense_examples": [],

            "reverse_examples": []

        }

        for sid in queries
    }

    total = 0

    with gzip.open(
        NCRNA,
        "rt",
        encoding="utf-8"
    ) as f:

        for header, sequence in fasta_records(f):

            total += 1

            for sid, query in queries.items():

                result = hits[sid]

                # ================================
                # 12. CORRISPONDENZA DIRETTA
                # ================================

                if query in sequence:

                    result[
                        "sense_records"
                    ] += 1

                    if (
                        len(
                            result["sense_examples"]
                        ) < 4
                    ):

                        result[
                            "sense_examples"
                        ].append(
                            header[:200]
                        )

                # ================================
                # 13. COMPLEMENTO INVERSO
                # ================================

                reverse = (
                    reverse_complement(query)
                )

                if reverse in sequence:

                    result[
                        "reverse_complement_records"
                    ] += 1

                    if (
                        len(
                            result["reverse_examples"]
                        ) < 4
                    ):

                        result[
                            "reverse_examples"
                        ].append(
                            header[:200]
                        )

    if not total:
        raise ValueError(
            "Riferimento ncRNA vuoto"
        )

    return total, hits


# =====================================================
# 14. PROGRAMMA PRINCIPALE
# =====================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    mode.add_argument(
        "--audit",
        action="store_true"
    )

    mode.add_argument(
        "--run",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 70)

    print(
        "AI PRIMETTA V3.3 | STEP 35 "
        "- CONTESTO GENOMICO E ncRNA"
    )

    print("=" * 70)

    # =================================================
    # 15. CONTROLLO DATI STEP33/34
    # =================================================

    (
        report,
        loci_path,
        queries,
        rows
    ) = read_inputs()

    print(
        "[PASS] STEP33/34, "
        "589 loci e checksum verificati"
    )

    # =================================================
    # 16. IDENTIFICAZIONE GENI
    # =================================================

    annotated = [

        r for r in rows

        if r["gene_overlap"] == "yes"
    ]

    print(
        "[GENI] Loci con sovrapposizione "
        f"genica: {len(annotated)}"
    )

    for r in annotated:

        print(

            f"[GENE] {r['srna_id']} "

            f"| {r['contig']}:"
            f"{r['start_1based']}-"
            f"{r['end_1based']} "

            f"({r['strand']}) "

            f"| ID={r['gene_ids']} "

            f"| biotipi={r['gene_biotypes']}"
        )

    # =================================================
    # 17. DISTRIBUZIONE SUL CROMOSOMA 6
    # =================================================

    chr6 = sorted(

        int(
            r["start_1based"]
        )

        for r in rows

        if (
            r["srna_id"] == "sRNA_00000037"
            and r["contig"] == "chr6"
        )
    )

    bins = Counter(

        (position - 1) // 1000000

        for position in chr6
    )

    if chr6:

        print(

            f"[CHR6] corrispondenze="
            f"{len(chr6)} "

            f"| prima={chr6[0]} "

            f"| ultima={chr6[-1]}"
        )

        print(

            "[CHR6] Regioni di 1 Mb "
            "piu popolate: "

            + "; ".join(

                f"{b * 1000000 + 1}:{n}"

                for b, n in bins.most_common(5)
            )
        )

    # =================================================
    # 18. ANALISI FASTA ncRNA
    # =================================================

    count, matches = inspect_ncrna(
        queries
    )

    print(

        "[NCRNA] Trascritti ncRNA "
        f"esaminati: {count}"
    )

    for sid, result in matches.items():

        print(

            f"[NCRNA] {sid} "

            f"| sense="
            f"{result['sense_records']} "

            f"| reverse_complement="
            f"{result['reverse_complement_records']}"
        )

        for example in result[
            "sense_examples"
        ]:

            print(
                "[ESEMPIO SENSE]",
                example
            )

        for example in result[
            "reverse_examples"
        ]:

            print(
                "[ESEMPIO REVERSE]",
                example
            )

    # =================================================
    # 19. LIMITAZIONI
    # =================================================

    print(
        "[LIMITE] Le corrispondenze ncRNA "
        "non dimostrano biogenesi o funzione; "
        "l'assenza di match non esclude "
        "RNA non annotati"
    )

    # =================================================
    # 20. SALVATAGGIO RISULTATI
    # =================================================

    if args.run:

        stamp = datetime.now(
            timezone.utc
        ).strftime(
            "%Y%m%dT%H%M%S%fZ"
        )

        folder = RUNS / (
            "STEP35_CONTESTO_SRNA_"
            + stamp
        )

        folder.mkdir(
            parents=True,
            exist_ok=False
        )

        output = (
            folder
            / "RAPPORTO_STEP35.json"
        )

        payload = {

            "step": 35,

            "status":
                "EXPLORATORY_NOT_CLASSIFIED",

            "source_step34":
                str(report.relative_to(ROOT)),

            "source_loci_sha256":
                sha256(loci_path),

            "ncrna_fasta":
                str(NCRNA.relative_to(ROOT)),

            "ncrna_sha256":
                sha256(NCRNA),

            "ncrna_transcripts":
                count,

            "matches":
                matches,

            "gene_overlap_loci":
                annotated,

            "chr6_srna00000037_positions": {

                "count":
                    len(chr6),

                "min":
                    min(chr6) if chr6 else None,

                "max":
                    max(chr6) if chr6 else None,

                "bins_1Mb":
                    dict(bins)
            },

            "training_performed":
                False,

            "gui_modified":
                False,

            "docker_modified":
                False
        }

        output.write_text(

            json.dumps(
                payload,
                indent=2,
                ensure_ascii=False
            ) + "\n",

            encoding="utf-8"
        )

        print(
            "[REPORT]",
            output
        )

    print(
        "[OK] Sorgenti, GUI, Docker "
        "e modelli invariati"
    )


# =====================================================
# 21. AVVIO
# =====================================================

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