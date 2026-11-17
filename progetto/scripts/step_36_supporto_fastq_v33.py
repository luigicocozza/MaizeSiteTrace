
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP36

Riscontro delle 2 small RNA nelle letture
FASTQ grezze di PRJNA419326.

--audit:
    Prime 10000 letture per libreria.
    Nessun file generato.

--run:
    Tutte le letture delle 12 librerie.
    Genera TSV e JSON.

Nessuna modifica a GUI, Docker,
modelli ML/GNN o dati originali.
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


# ====================================================
# 1. CONFIGURAZIONE
# ====================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

FASTQ_DIR = (
    ROOT
    / "data/public_ena/PRJNA419326/fastq/mirnaseq"
)

EXPECTED = {

    "sRNA_00000037":
        "TTCTTGGATTTATGAAAGACGAACA",

    "sRNA_00000435":
        "AGATGATACGAGCACCTCTGAGCT"
}

PREVIEW_READS = 10000


# ====================================================
# 2. CHECKSUM SHA256
# ====================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1048576),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# ====================================================
# 3. RECUPERO SEQUENZE ORIGINALI
# ====================================================

def locate_source():

    manifests = sorted(
        RUNS.glob(
            "STEP33_LOCI_GENOMICI_*/"
            "RAPPORTO_STEP33.json"
        )
    )

    if not manifests:

        raise FileNotFoundError(
            "RAPPORTO_STEP33.json assente"
        )

    report = manifests[-1].resolve()

    obj = json.loads(
        report.read_text(
            encoding="utf-8"
        )
    )

    if obj.get("step") != 33:

        raise ValueError(
            "Manifest STEP33 non valido"
        )

    queries = Path(
        obj["srna_fasta"]
    )

    if not queries.is_absolute():

        queries = (
            ROOT / queries
        )

    queries = queries.resolve()

    if (
        not queries.is_file()
        or not queries.is_relative_to(
            RUNS.resolve()
        )
    ):

        raise ValueError(
            "FASTA sorgente non valido"
        )

    if (
        sha256(queries)
        != obj["srna_fasta_sha256"]
    ):

        raise ValueError(
            "SHA256 FASTA small RNA differente"
        )

    observed = {}

    name = None
    seq = []

    for line in queries.read_text(
        encoding="utf-8-sig"
    ).splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if name is not None:

                if (
                    name in observed
                    or not seq
                ):

                    raise ValueError(
                        "FASTA duplicato/vuoto"
                    )

                observed[name] = (
                    "".join(seq)
                    .upper()
                    .replace("U", "T")
                )

            name = (
                line[1:].split()[0]
            )

            seq = []

        else:

            if name is None:

                raise ValueError(
                    "FASTA privo di intestazione"
                )

            seq.append(line)

    if name is not None:

        if (
            name in observed
            or not seq
        ):

            raise ValueError(
                "Ultimo record FASTA errato"
            )

        observed[name] = (
            "".join(seq)
            .upper()
            .replace("U", "T")
        )

    if observed != EXPECTED:

        raise ValueError(
            "Sequenze non corrispondenti "
            "ai dati originali STEP35"
        )

    return report, queries


# ====================================================
# 4. COMPLEMENTO INVERSO DEL DNA
# ====================================================

def reverse_complement(seq):

    return seq.translate(
        str.maketrans(
            "ACGT",
            "TGCA"
        )
    )[::-1]


# ====================================================
# 5. LETTURA E CONTEGGIO FASTQ
# ====================================================

def scan_fastq(
    path,
    max_reads=None
):

    counts = Counter()

    queries = {

        sid: (
            sequence,
            reverse_complement(sequence)
        )

        for sid, sequence in EXPECTED.items()
    }

    with gzip.open(
        path,
        "rt",
        encoding="ascii"
    ) as f:

        while (
            max_reads is None
            or counts["reads"] < max_reads
        ):

            header = f.readline()

            if not header:
                break

            seq = (
                f.readline()
                .rstrip("\r\n")
            )

            plus = f.readline()

            quality = (
                f.readline()
                .rstrip("\r\n")
            )

            if (
                not header.startswith("@")
                or not plus.startswith("+")
                or not seq
                or len(seq) != len(quality)
                or not set(seq.upper())
                <= set("ACGTN")
            ):

                raise ValueError(
                    "Record FASTQ non valido: "
                    f"{path.name}, "
                    f"record={counts['reads'] + 1}"
                )

            seq = seq.upper()

            counts["reads"] += 1

            # ========================================
            # 6. MATCH DELLE DUE SMALL RNA
            # ========================================

            for sid, (
                query,
                reverse
            ) in queries.items():

                forward_match = (
                    query in seq
                )

                reverse_match = (
                    reverse in seq
                )

                if forward_match:

                    counts[
                        f"{sid}|sense_anywhere"
                    ] += 1

                    if seq.startswith(query):

                        counts[
                            f"{sid}|sense_prefix"
                        ] += 1

                if reverse_match:

                    counts[
                        f"{sid}|reverse_anywhere"
                    ] += 1

                    if seq.startswith(reverse):

                        counts[
                            f"{sid}|reverse_prefix"
                        ] += 1

                if (
                    forward_match
                    or reverse_match
                ):

                    counts[
                        f"{sid}|read_any_orientation"
                    ] += 1

    return counts


# ====================================================
# 7. PROGRAMMA PRINCIPALE
# ====================================================

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

    print("=" * 68)

    print(
        "AI PRIMETTA V3.3 | STEP 36 "
        "- SUPPORTO NELLE LETTURE FASTQ"
    )

    print("=" * 68)

    # ================================================
    # 8. VERIFICA DATI DI PARTENZA
    # ================================================

    report, queries = locate_source()

    paths = sorted(
        FASTQ_DIR.glob(
            "*.fastq.gz"
        )
    )

    if (
        len(paths) != 12
        or len({
            p.name
            for p in paths
        }) != 12
    ):

        raise ValueError(
            "Attesi 12 FASTQ originali; "
            f"trovati {len(paths)}"
        )

    print(
        "[PASS] STEP33, sequenze originali "
        "e checksum verificati; "
        "12 FASTQ presenti"
    )

    rows = []

    totals = {
        sid: Counter()
        for sid in EXPECTED
    }

    # ================================================
    # 9. ANALISI DELLE 12 LIBRERIE
    # ================================================

    for index, path in enumerate(
        paths,
        1
    ):

        counted = scan_fastq(

            path,

            (
                PREVIEW_READS
                if args.audit
                else None
            )
        )

        n = counted["reads"]

        if n == 0:

            raise ValueError(
                f"FASTQ vuoto: {path.name}"
            )

        for sid in EXPECTED:

            data = {

                "sample":
                    path.name,

                "srna_id":
                    sid,

                "reads_examined":
                    n,

                "sense_prefix":
                    counted[
                        f"{sid}|sense_prefix"
                    ],

                "sense_anywhere":
                    counted[
                        f"{sid}|sense_anywhere"
                    ],

                "reverse_prefix":
                    counted[
                        f"{sid}|reverse_prefix"
                    ],

                "reverse_anywhere":
                    counted[
                        f"{sid}|reverse_anywhere"
                    ],

                "reads_any_orientation":
                    counted[
                        f"{sid}|read_any_orientation"
                    ],

                "reads_any_per_million_raw":
                    round(

                        1000000
                        * counted[
                            f"{sid}|read_any_orientation"
                        ]
                        / n,

                        4
                    ),

                "scan_mode":
                    (
                        "sample_first_10000"
                        if args.audit
                        else "complete"
                    )
            }

            rows.append(data)

            for key in (

                "reads_examined",
                "sense_prefix",
                "sense_anywhere",
                "reverse_prefix",
                "reverse_anywhere",
                "reads_any_orientation"

            ):

                totals[sid][key] += (
                    data[key]
                )

        print(

            f"[FASTQ] {index}/12 "

            f"{path.name} "

            f"| letture={n:,}",

            flush=True
        )

    # ================================================
    # 10. RISULTATO COMPLESSIVO
    # ================================================

    for sid in EXPECTED:

        t = totals[sid]

        print(

            f"[RISULTATO] {sid}: "

            f"letture="
            f"{t['reads_examined']:,} "

            f"| 5prima_sense="
            f"{t['sense_prefix']:,} "

            f"| ovunque_sense="
            f"{t['sense_anywhere']:,} "

            f"| 5prima_reverse="
            f"{t['reverse_prefix']:,} "

            f"| ovunque_reverse="
            f"{t['reverse_anywhere']:,} "

            f"| letture_con_match="
            f"{t['reads_any_orientation']:,}"
        )

    print(

        "[LIMITE] Conteggi su sequenze FASTQ "
        "grezze: non rappresentano "
        "abbondanza normalizzata, prova "
        "di origine genomica o "
        "classificazione biologica."
    )

    # ================================================
    # 11. AUDIT SENZA SCRITTURA
    # ================================================

    if args.audit:

        print(
            "[OK] Audit campionario, "
            "nessun file modificato"
        )

        return

    # ================================================
    # 12. ESPORTAZIONE COMPLETA
    # ================================================

    timestamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    out = RUNS / (
        "STEP36_SUPPORTO_FASTQ_"
        + timestamp
    )

    out.mkdir(
        parents=True,
        exist_ok=False
    )

    tsv = (
        out
        / "SUPPORTO_2_SRNA_12_FASTQ.tsv"
    )

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

    # ================================================
    # 13. MANIFEST SCIENTIFICO
    # ================================================

    manifest = (
        out
        / "RAPPORTO_STEP36.json"
    )

    payload = {

        "step": 36,

        "status":
            "RAW_READ_EXACT_SEQUENCE_SUPPORT_NOT_CLASSIFICATION",

        "source_step33":
            str(
                report.relative_to(ROOT)
            ),

        "source_fasta_sha256":
            sha256(queries),

        "sample_count":
            len(paths),

        "summary": {

            sid: dict(values)

            for sid, values in totals.items()
        },

        "output_tsv":
            str(tsv.relative_to(ROOT)),

        "output_tsv_sha256":
            sha256(tsv),

        "limitations": [

            "Le letture grezze possono "
            "contenere adattatori.",

            "Un match interno alla lettura "
            "non dimostra una small RNA matura.",

            "Non sono state rimosse "
            "duplicazioni PCR.",

            "Non viene inferita la condizione "
            "sperimentale dei campioni.",

            "Non e' dimostrata l'origine "
            "ribosomale o la regolazione "
            "di un gene bersaglio."
        ],

        "gui_modified":
            False,

        "docker_modified":
            False,

        "training_performed":
            False
    }

    manifest.write_text(

        json.dumps(
            payload,
            ensure_ascii=False,
            indent=2
        ) + "\n",

        encoding="utf-8"
    )

    print(
        "[TSV]",
        tsv
    )

    print(
        "[REPORT]",
        manifest
    )

    print(
        "[OK] STEP36 completato; "
        "GUI, Docker e modelli invariati"
    )


# ====================================================
# 14. AVVIO
# ====================================================

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
