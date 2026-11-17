
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 19

Screening diagnostico di complementarita':
- HSP di 20 nt dal 5' della small RNA
- Appaiamenti Watson-Crick e G:U
- Regione seed 2-13: massimo 2 mismatch non compatibili
- Nessun gap
- Punteggio diagnostico, NON psRNATarget

Riutilizza le funzioni di lettura dello STEP 18.

--audit : verifica input e test sintetici
--scan  : ricerca sul trascrittoma B73

Non modifica GUI, Docker, modelli o dati originali.
"""

import argparse
import csv
import heapq
import itertools
import json
import sys

from collections import Counter
from datetime import datetime, timezone

import step_18_verifica_target_v33 as previous


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = previous.ROOT
RUNS = previous.RUNS

HSP = 20

# Posizioni nella sequenza bersaglio orientata 5'->3'.
# Corrispondono alle posizioni 13..2 della sRNA.
SEED_START = 7
SEED_END = 19

MAX_SEED_HARD = 2

MAX_DIAGNOSTIC_SCORE = 5.0

TOP_N = 30

# Basi DNA; T rappresenta U.
ALLOWED_TARGET = {
    "A": ("T",),
    "C": ("G",),
    "G": ("C", "T"),
    "T": ("A", "G"),
}


# ============================================================
# 2. BASI COMPLEMENTARI E G:U
# ============================================================

def options_for_srna(srna):

    if not set(srna) <= set(ALLOWED_TARGET):

        raise ValueError(
            "Sequenza sRNA con basi non A/C/G/T"
        )

    return [
        ALLOWED_TARGET[base]
        for base in reversed(srna[:HSP])
    ]


# ============================================================
# 3. SEED DI RICERCA
# ============================================================

def anchors(options):

    """
    Tre segmenti disgiunti di quattro basi.

    Con massimo due mismatch non compatibili
    nel seed, almeno uno dei tre segmenti
    deve essere completamente compatibile,
    considerando anche i G:U.
    """

    groups = []

    for start in range(
        SEED_START,
        SEED_END,
        4
    ):

        combinations = itertools.product(
            *options[start:start + 4]
        )

        words = {
            "".join(values)
            for values in combinations
        }

        groups.append(
            (
                start,
                sorted(words)
            )
        )

    assert len(groups) == 3

    assert all(
        len(word) == 4
        for _, words in groups
        for word in words
    )

    return groups


# ============================================================
# 4. RICERCA POSIZIONI CANDIDATE
# ============================================================

def candidate_positions(
    transcript,
    anchor_groups
):

    positions = set()

    largest = len(transcript) - HSP

    if largest < 0:
        return positions

    for offset, words in anchor_groups:

        for word in words:

            index = transcript.find(word)

            while index >= 0:

                start = index - offset

                if 0 <= start <= largest:
                    positions.add(start)

                index = transcript.find(
                    word,
                    index + 1
                )

    return positions


# ============================================================
# 5. VALUTAZIONE APPAIAMENTO
# ============================================================

def classify(
    window,
    wc,
    allowed
):

    seed_hard = 0
    hard = 0
    wobble = 0

    for i, base in enumerate(window):

        # Appaiamento Watson-Crick.
        if base == wc[i]:
            continue

        # Appaiamento G:U.
        if base in allowed[i]:

            wobble += 1

        else:

            hard += 1

            if SEED_START <= i < SEED_END:

                seed_hard += 1

                if seed_hard > MAX_SEED_HARD:
                    return None

    score = hard + 0.5 * wobble

    return (
        seed_hard,
        hard,
        wobble,
        score
    )


# ============================================================
# 6. SCANSIONE DI UN TRASCRITTO
# ============================================================

def scan_one_transcript(
    sequence,
    srna
):

    allowed = options_for_srna(srna)

    wc = previous.revcomp(
        srna[:HSP]
    )

    groups = anchors(allowed)

    positions = candidate_positions(
        sequence,
        groups
    )

    for start in sorted(positions):

        window = sequence[
            start:start + HSP
        ]

        if not set(window) <= set("ACGT"):
            continue

        result = classify(
            window,
            wc,
            allowed
        )

        if result is not None:

            yield (
                start,
                window,
                *result
            )


# ============================================================
# 7. TEST SINTETICO
# ============================================================

def self_test():

    srna = "AGCTGATCGTACGTCAGTGAACGT"

    wc = previous.revcomp(
        srna[:HSP]
    )

    allowed = options_for_srna(srna)

    altered = list(wc)

    # Due mismatch reali in due blocchi diversi.
    for pos in (
        SEED_START,
        SEED_START + 4
    ):

        altered[pos] = next(
            base
            for base in "ACGT"
            if base not in allowed[pos]
        )

    # Un appaiamento wobble nel terzo blocco.
    wobble_pos = next(
        pos
        for pos in range(
            SEED_START + 8,
            SEED_END
        )
        if len(allowed[pos]) == 2
    )

    altered[wobble_pos] = next(
        base
        for base in allowed[wobble_pos]
        if base != wc[wobble_pos]
    )

    sequence = (
        "AAAAA"
        + "".join(altered)
        + "CCCCC"
    )

    matches = list(
        scan_one_transcript(
            sequence,
            srna
        )
    )

    assert any(
        hit[0] == 5
        and hit[2:5] == (2, 2, 1)
        for hit in matches
    ), matches

    # Tre mismatch veri, uno per blocco:
    # deve essere escluso dal vincolo seed.
    third = list(wc)

    for pos in (
        SEED_START,
        SEED_START + 4,
        SEED_START + 8
    ):

        third[pos] = next(
            base
            for base in "ACGT"
            if base not in allowed[pos]
        )

    control = (
        "AAAAA"
        + "".join(third)
        + "CCCCC"
    )

    assert not any(
        hit[0] == 5
        for hit in scan_one_transcript(
            control,
            srna
        )
    )

    return (
        "PASS: 2 mismatch + 1 G:U; "
        "esclusione con 3 mismatch nel seed"
    )


# ============================================================
# 8. AUDIT INPUT
# ============================================================

def audit():

    reference = previous.select_reference()

    input_file = previous.select_input()

    mapping_file = previous.find_mapping(
        reference
    )

    queries = previous.load_queries(
        input_file
    )

    print(
        "[INPUT]",
        input_file.relative_to(ROOT)
    )

    print(
        "[REFERENCE]",
        reference.relative_to(ROOT)
    )

    print(
        "[TX2GENE]",
        (
            mapping_file.relative_to(ROOT)
            if mapping_file
            else "assente"
        )
    )

    for ident, srna in queries:

        groups = anchors(
            options_for_srna(srna)
        )

        print(
            f"[sRNA] {ident}"
            f" | {len(srna)} nt"
            f" | HSP={HSP}"
            f" | seed=2-13"
            f" | varianti="
            f"{[len(group[1]) for group in groups]}"
        )

    print(
        "[SELF-TEST]",
        self_test()
    )

    print(
        "[LIMITE] Screening diagnostico: "
        "HSP20, G:U, senza gap"
    )

    return (
        reference,
        input_file,
        mapping_file,
        queries
    )


# ============================================================
# 9. SCANSIONE DEL TRASCRITTOMA
# ============================================================

def scan():

    (
        reference,
        input_file,
        mapping_file,
        queries
    ) = audit()

    mapping = previous.parse_mapping(
        mapping_file
    )

    prepared = []

    for ident, srna in queries:

        allowed = options_for_srna(srna)

        wc = previous.revcomp(
            srna[:HSP]
        )

        prepared.append(
            (
                ident,
                srna,
                wc,
                allowed,
                anchors(allowed)
            )
        )

    stats = {
        ident: Counter()
        for ident, _ in queries
    }

    top = {
        ident: []
        for ident, _ in queries
    }

    seen = set()

    ordinal = 0
    n_transcripts = 0

    for transcript, sequence in (
        previous.fasta_records(reference)
    ):

        if transcript in seen:

            raise RuntimeError(
                "Trascritto duplicato: "
                + transcript
            )

        seen.add(transcript)

        n_transcripts += 1

        for (
            ident,
            srna,
            wc,
            allowed,
            groups
        ) in prepared:

            positions = candidate_positions(
                sequence,
                groups
            )

            for start in positions:

                window = sequence[
                    start:start + HSP
                ]

                if not set(window) <= set("ACGT"):
                    continue

                result = classify(
                    window,
                    wc,
                    allowed
                )

                if result is None:
                    continue

                (
                    seed_hard,
                    hard,
                    wobble,
                    score
                ) = result

                stats[ident][
                    "seed_compatible"
                ] += 1

                if score > MAX_DIAGNOSTIC_SCORE:
                    continue

                stats[ident][
                    "diagnostic_score_le_5"
                ] += 1

                ordinal += 1

                hit = {

                    "srna_id": ident,

                    "transcript_id": transcript,

                    "gene_id": mapping.get(
                        transcript,
                        ""
                    ),

                    "start_1based": start + 1,

                    "end_1based": start + HSP,

                    "srna_length": len(srna),

                    "hsp_nt": HSP,

                    "hard_mismatches": hard,

                    "seed_hard_mismatches": seed_hard,

                    "gu_wobbles": wobble,

                    "diagnostic_score": score,

                    "site_5to3_DNA": window,

                    "perfect_reverse_complement_5to3_DNA":
                        wc
                }

                # Mantiene soltanto i migliori TOP_N
                # senza conservare tutti i siti in RAM.

                rank = (
                    -score,
                    -hard,
                    -wobble,
                    ordinal,
                    hit
                )

                heap = top[ident]

                if len(heap) < TOP_N:

                    heapq.heappush(
                        heap,
                        rank
                    )

                elif rank > heap[0]:

                    heapq.heapreplace(
                        heap,
                        rank
                    )

        if n_transcripts % 10000 == 0:

            print(
                "[SCANSIONE]",
                n_transcripts,
                "trascritti",
                flush=True
            )

    # ========================================================
    # 10. OUTPUT
    # ========================================================

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP19_COMPLEMENTARITA_HSP20_"
        + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    tsv_path = (
        folder
        / "TOP_CANDIDATI_ESPLORATIVI.tsv"
    )

    fields = [
        "srna_id",
        "transcript_id",
        "gene_id",
        "start_1based",
        "end_1based",
        "srna_length",
        "hsp_nt",
        "hard_mismatches",
        "seed_hard_mismatches",
        "gu_wobbles",
        "diagnostic_score",
        "site_5to3_DNA",
        "perfect_reverse_complement_5to3_DNA",
    ]

    with tsv_path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as stream:

        writer = csv.DictWriter(
            stream,
            delimiter="\t",
            fieldnames=fields
        )

        writer.writeheader()

        for ident, _ in queries:

            best = sorted(
                top[ident],
                reverse=True
            )

            for (
                _,
                _,
                _,
                _,
                hit
            ) in best:

                writer.writerow(hit)

    # ========================================================
    # 11. REPORT
    # ========================================================

    results = {}

    for ident, _ in queries:

        results[ident] = {

            "seed_compatible_sites":
                stats[ident][
                    "seed_compatible"
                ],

            "sites_score_le_5":
                stats[ident][
                    "diagnostic_score_le_5"
                ],

            "top_sites_exported":
                len(top[ident])
        }

        print(
            f"[RISULTATO] {ident}"
            f" | seed compatibili="
            f"{results[ident]['seed_compatible_sites']}"
            f" | punteggio <=5="
            f"{results[ident]['sites_score_le_5']}"
            f" | top esportati="
            f"{results[ident]['top_sites_exported']}"
        )

    report = {

        "step": "19",

        "status": "COMPLETED",

        "transcripts_scanned": n_transcripts,

        "input": str(
            input_file.relative_to(ROOT)
        ),

        "reference": str(
            reference.relative_to(ROOT)
        ),

        "mapping": (
            str(mapping_file.relative_to(ROOT))
            if mapping_file
            else None
        ),

        "method": (
            "HSP20, reverse-complement, "
            "G:U, no gaps, "
            "seed positions 2-13"
        ),

        "seed_hard_mismatch_max": 2,

        "diagnostic_score_formula": (
            "1*nonpair + 0.5*G:U"
        ),

        "diagnostic_cutoff":
            MAX_DIAGNOSTIC_SCORE,

        "not_psRNATarget_score": True,

        "results": results,

        "limitations": [

            "Non equivale a psRNATarget",

            "Analizza soltanto i primi 20 nt "
            "di ogni small RNA",

            "Non gestisce gap, energia di legame "
            "o accessibilita' del target",

            "Le corrispondenze non costituiscono "
            "validazione biologica"
        ],

        "gui_modified": False,

        "models_modified": False,

        "training_performed": False
    }

    report_path = (
        folder / "RAPPORTO_STEP19.json"
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
        "[REPORT]",
        report_path
    )

    print(
        "[TSV]",
        tsv_path
    )

    print(
        "[OK] Nessuna GUI modificata"
    )

    print(
        "[OK] Nessun training"
    )


# ============================================================
# 12. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    group = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    group.add_argument(
        "--audit",
        action="store_true"
    )

    group.add_argument(
        "--scan",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 68)

    print(
        "AI PRIMETTA V3.3 | STEP 19"
    )

    print(
        "HSP20, G:U, SEED 2-13"
    )

    print("=" * 68)

    if args.audit:

        audit()

        print(
            "[OK] Audit senza modifiche"
        )

    else:

        scan()


if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__,
            str(exc),
            file=sys.stderr
        )

        sys.exit(1)
