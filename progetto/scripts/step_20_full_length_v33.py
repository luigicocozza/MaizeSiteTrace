
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 20

Rivaluta sulla lunghezza completa tutti i siti HSP20
dello STEP 19, mantenendo il legame con i trascritti B73.

--audit : verifica sequenze, duplicazioni e orientamento
--scan  : analisi completa e report

Nessuna modifica a GUI, Docker, checkpoint o modelli.
"""

import argparse
import csv
import json
import sys

from collections import Counter, defaultdict
from datetime import datetime, timezone

import step_18_verifica_target_v33 as s18
import step_19_diagnostica_appaiamento_v33 as s19


ROOT = s18.ROOT
RUNS = s18.RUNS

MAX_ROWS = 100000
TOP_N = 30
HSP_LIMIT = 5.0


# ==========================================================
# 1. VALUTAZIONE COMPLETA
# ==========================================================

def full_pairing(window, srna):

    if len(window) != len(srna):

        raise ValueError(
            "Lunghezza target diversa dalla sRNA"
        )

    allowed = [
        s19.ALLOWED_TARGET[b]
        for b in reversed(srna)
    ]

    perfect = s18.revcomp(srna)

    hard = 0
    wobble = 0

    bad_positions = []
    wobble_positions = []

    for i, (base, reference) in enumerate(
        zip(window, perfect)
    ):

        srna_position = len(srna) - i

        if base == reference:
            continue

        if base in allowed[i]:

            wobble += 1

            wobble_positions.append(
                srna_position
            )

        else:

            hard += 1

            bad_positions.append(
                srna_position
            )

    return {

        "full_hard": hard,

        "full_gu": wobble,

        "full_score": hard + 0.5 * wobble,

        "hard_srna_positions": ",".join(
            map(
                str,
                sorted(bad_positions)
            )
        ),

        "gu_srna_positions": ",".join(
            map(
                str,
                sorted(wobble_positions)
            )
        )
    }


# ==========================================================
# 2. ESTENSIONE DEL SITO HSP20
# ==========================================================

def extend_site(
    transcript,
    hsp_start,
    srna
):

    tail = len(srna) - s19.HSP

    # Le basi aggiuntive del target si trovano
    # a monte della regione HSP20 nel trascritto.

    start = hsp_start - tail
    stop = hsp_start + s19.HSP

    if start < 0 or stop > len(transcript):
        return None

    sequence = transcript[start:stop]

    if (
        len(sequence) != len(srna)
        or not set(sequence) <= set("ACGT")
    ):
        return None

    assert (
        sequence[tail:]
        == transcript[hsp_start:stop]
    )

    return (
        start,
        stop,
        sequence,
        full_pairing(sequence, srna)
    )


# ==========================================================
# 3. TEST SINTETICI
# ==========================================================

def test_algorithm(queries):

    for ident, srna in queries:

        perfect = s18.revcomp(srna)

        hsp = s18.revcomp(
            srna[:s19.HSP]
        )

        assert perfect.endswith(hsp), ident

        transcript = (
            "AAAAA" + perfect + "CCCC"
        )

        hsp_start = (
            5 + len(srna) - s19.HSP
        )

        anchors = s19.anchors(
            s19.options_for_srna(srna)
        )

        assert hsp_start in (
            s19.candidate_positions(
                transcript,
                anchors
            )
        )

        result = extend_site(
            transcript,
            hsp_start,
            srna
        )

        assert result is not None

        start, stop, sequence, score = result

        assert (
            start,
            stop,
            sequence
        ) == (
            5,
            5 + len(srna),
            perfect
        )

        assert score["full_score"] == 0

        # Il sito non deve essere esteso
        # oltre l'inizio del trascritto.

        assert extend_site(
            hsp + "CCCC",
            0,
            srna
        ) is None

        # Introduce un mismatch nella parte
        # non coperta dallo STEP 19.

        altered = list(perfect)

        options = (
            s19.ALLOWED_TARGET[
                srna[-1]
            ]
        )

        altered[0] = next(
            c for c in "ACGT"
            if c not in options
        )

        mutated = "".join(altered)

        assert mutated.endswith(hsp)

        modified = full_pairing(
            mutated,
            srna
        )

        assert modified["full_hard"] == 1

        assert modified["full_score"] == 1

    return (
        "PASS: estensione a monte, orientamento "
        "e scoring per tutte le sRNA"
    )


# ==========================================================
# 4. AUDIT DELLE CINQUE SEQUENZE
# ==========================================================

def audit():

    if (
        s19.HSP != 20
        or s19.MAX_SEED_HARD != 2
        or s19.MAX_DIAGNOSTIC_SCORE != HSP_LIMIT
    ):

        raise RuntimeError(
            "Parametri STEP19 modificati"
        )

    input_file = s18.select_input()

    reference = s18.select_reference()

    mapping_file = s18.find_mapping(
        reference
    )

    queries = s18.load_queries(
        input_file
    )

    by_hsp = defaultdict(list)

    full_sequences = set()

    print(
        "[INPUT]",
        input_file.relative_to(ROOT)
    )

    print(
        "[RIFERIMENTO]",
        reference.relative_to(ROOT)
    )

    for ident, srna in queries:

        full_sequences.add(srna)

        by_hsp[
            srna[:20]
        ].append(ident)

        print(
            f"[sRNA] {ident}: "
            f"{len(srna)} nt "
            f"| primi20={srna[:20]} "
            f"| coda3'={srna[20:]}"
        )

    for hsp, ids in by_hsp.items():

        if len(ids) > 1:

            print(
                "[HSP IDENTICO]",
                hsp,
                "|",
                ", ".join(ids)
            )

    print(
        "[DIVERSITA']",
        len(full_sequences),
        "sequenze complete distinte;",
        len(by_hsp),
        "HSP da 20 nt distinti"
    )

    print(
        "[SELF-TEST]",
        test_algorithm(queries)
    )

    return (
        input_file,
        reference,
        mapping_file,
        queries
    )


# ==========================================================
# 5. RECUPERO DEL REPORT STEP19
# ==========================================================

def prior_step19():

    reports = sorted(
        RUNS.glob(
            "STEP19_COMPLEMENTARITA_HSP20_*/"
            "RAPPORTO_STEP19.json"
        )
    )

    if not reports:

        raise RuntimeError(
            "Rapporto STEP19 non trovato"
        )

    path = reports[-1]

    document = json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )

    return path, document


# ==========================================================
# 6. SCANSIONE COMPLETA
# ==========================================================

def run_scan():

    (
        input_file,
        reference,
        mapping_file,
        queries
    ) = audit()

    prev_path, prev = prior_step19()

    if (
        prev.get("input")
        != str(input_file.relative_to(ROOT))
        or prev.get("reference")
        != str(reference.relative_to(ROOT))
    ):

        raise RuntimeError(
            "Input/riferimento differente "
            "dal report STEP19"
        )

    mapping = s18.parse_mapping(
        mapping_file
    )

    prepared = []

    for ident, srna in queries:

        allowed = s19.options_for_srna(
            srna
        )

        wc = s18.revcomp(
            srna[:20]
        )

        prepared.append(
            (
                ident,
                srna,
                wc,
                allowed,
                s19.anchors(allowed)
            )
        )

    hsp_counts = Counter()
    boundary_counts = Counter()
    valid_counts = Counter()

    scores = defaultdict(Counter)

    records = {
        ident: []
        for ident, _ in queries
    }

    transcripts = 0

    seen = set()

    # ------------------------------------------------------
    # ANALISI TRASCRITTI B73
    # ------------------------------------------------------

    for tx, seq in s18.fasta_records(
        reference
    ):

        if tx in seen:

            raise RuntimeError(
                f"Trascritto duplicato: {tx}"
            )

        seen.add(tx)
        transcripts += 1

        for (
            ident,
            srna,
            wc,
            allowed,
            anchors
        ) in prepared:

            positions = (
                s19.candidate_positions(
                    seq,
                    anchors
                )
            )

            for hsp_start in positions:

                win = seq[
                    hsp_start:
                    hsp_start + s19.HSP
                ]

                if not set(win) <= set("ACGT"):
                    continue

                result = s19.classify(
                    win,
                    wc,
                    allowed
                )

                if result is None:
                    continue

                (
                    seed_bad,
                    hsp_bad,
                    hsp_gu,
                    hsp_score
                ) = result

                if hsp_score > HSP_LIMIT:
                    continue

                hsp_counts[ident] += 1

                extension = extend_site(
                    seq,
                    hsp_start,
                    srna
                )

                if extension is None:

                    boundary_counts[ident] += 1
                    continue

                (
                    left,
                    right,
                    full_window,
                    values
                ) = extension

                valid_counts[ident] += 1

                full_score = values[
                    "full_score"
                ]

                if full_score <= 5:
                    scores[ident]["le_5"] += 1

                if full_score <= 7:
                    scores[ident]["le_7"] += 1

                if full_score <= 10:
                    scores[ident]["le_10"] += 1

                if len(records[ident]) >= MAX_ROWS:

                    raise RuntimeError(
                        "Troppi siti per "
                        + ident
                    )

                records[ident].append({

                    "srna_id": ident,

                    "srna_length": len(srna),

                    "transcript_id": tx,

                    "gene_id": mapping.get(
                        tx,
                        ""
                    ),

                    "start_1based": left + 1,

                    "end_1based": right,

                    "hsp_start_1based":
                        hsp_start + 1,

                    "hsp_hard": hsp_bad,

                    "hsp_gu": hsp_gu,

                    "hsp_score": hsp_score,

                    "seed_hard": seed_bad,

                    "full_site_DNA_5to3":
                        full_window,

                    "full_perfect_target_DNA_5to3":
                        s18.revcomp(srna),

                    **values
                })

        if transcripts % 10000 == 0:

            print(
                "[SCANSIONE]",
                transcripts,
                "trascritti",
                flush=True
            )

    # ------------------------------------------------------
    # VERIFICA COERENZA CON STEP19
    # ------------------------------------------------------

    print(
        "[REPORT CONFRONTO]",
        prev_path.relative_to(ROOT)
    )

    for ident, _ in queries:

        expected = (
            prev["results"][ident][
                "sites_score_le_5"
            ]
        )

        if hsp_counts[ident] != expected:

            raise RuntimeError(
                f"Incoerenza {ident}: "
                f"STEP19={expected}, "
                f"STEP20={hsp_counts[ident]}"
            )

    # ======================================================
    # 7. ESPORTAZIONE RISULTATI
    # ======================================================

    fields = [

        "srna_id",
        "srna_length",
        "transcript_id",
        "gene_id",
        "start_1based",
        "end_1based",
        "hsp_start_1based",
        "hsp_hard",
        "hsp_gu",
        "hsp_score",
        "seed_hard",
        "full_hard",
        "full_gu",
        "full_score",
        "hard_srna_positions",
        "gu_srna_positions",
        "full_site_DNA_5to3",
        "full_perfect_target_DNA_5to3"
    ]

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP20_COMPLEMENTARITA_FULL_LENGTH_"
        + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    full_tsv = (
        folder
        / "TUTTI_SITI_RIVALUTATI_FULL_LENGTH.tsv"
    )

    top_tsv = (
        folder
        / "TOP30_PER_SRNA_FULL_LENGTH.tsv"
    )

    def sortkey(row):

        return (
            row["full_score"],
            row["full_hard"],
            row["full_gu"],
            row["hsp_score"],
            row["transcript_id"],
            row["start_1based"]
        )

    for path, max_rows in (
        (full_tsv, None),
        (top_tsv, TOP_N)
    ):

        with path.open(
            "w",
            encoding="utf-8",
            newline=""
        ) as stream:

            writer = csv.DictWriter(
                stream,
                fieldnames=fields,
                delimiter="\t"
            )

            writer.writeheader()

            for ident, _ in queries:

                data = sorted(
                    records[ident],
                    key=sortkey
                )

                writer.writerows(
                    data
                    if max_rows is None
                    else data[:max_rows]
                )

    # ======================================================
    # 8. RIEPILOGO PER SMALL RNA
    # ======================================================

    summary = {}

    for ident, srna in queries:

        best = min(
            (
                row["full_score"]
                for row in records[ident]
            ),
            default=None
        )

        genes = {
            row["gene_id"]
            for row in records[ident]
            if row["gene_id"]
        }

        summary[ident] = {

            "srna_length": len(srna),

            "hsp20_score_le_5":
                hsp_counts[ident],

            "full_length_evaluable":
                valid_counts[ident],

            "full_length_not_evaluable_boundary":
                boundary_counts[ident],

            "full_score_le_5":
                scores[ident]["le_5"],

            "full_score_le_7":
                scores[ident]["le_7"],

            "full_score_le_10":
                scores[ident]["le_10"],

            "best_full_score": best,

            "distinct_genes_mapped":
                len(genes)
        }

        print(
            f"[RISULTATO] {ident}"
            f" | HSP={hsp_counts[ident]}"
            f" | completi={valid_counts[ident]}"
            f" | bordo={boundary_counts[ident]}"
            f" | full<=5={scores[ident]['le_5']}"
            f" | full<=7={scores[ident]['le_7']}"
            f" | migliore={best}"
        )

    # ======================================================
    # 9. RAPPORTO SCIENTIFICO
    # ======================================================

    first20_groups = defaultdict(list)

    for ident, srna in queries:

        first20_groups[
            srna[:20]
        ].append(ident)

    report = {

        "step": "20",

        "status": "COMPLETED",

        "source_step19": str(
            prev_path.relative_to(ROOT)
        ),

        "reference": str(
            reference.relative_to(ROOT)
        ),

        "input": str(
            input_file.relative_to(ROOT)
        ),

        "transcripts_processed":
            transcripts,

        "method": (
            "Rivalutazione full-length, "
            "senza gap, dei siti HSP20 "
            "con score <=5; "
            "estensione a monte"
        ),

        "score": (
            "diagnostico = nonpair + 0.5*G:U; "
            "NON psRNATarget expectation"
        ),

        "first20_groups": dict(
            first20_groups
        ),

        "results": summary,

        "full_tsv": str(
            full_tsv.relative_to(ROOT)
        ),

        "top_tsv": str(
            top_tsv.relative_to(ROOT)
        ),

        "limitations": [

            "Non gestisce inserzioni o delezioni",

            "Non valuta accessibilita' del target",

            "Non calcola energia di legame",

            "Non costituisce validazione biologica",

            "Risultati condizionati dal "
            "prefiltro HSP20 dello STEP19"
        ],

        "gui_modified": False,

        "training_performed": False
    }

    json_path = (
        folder / "RAPPORTO_STEP20.json"
    )

    json_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    print(
        "[REPORT]",
        json_path
    )

    print(
        "[TSV TUTTI]",
        full_tsv
    )

    print(
        "[TSV TOP]",
        top_tsv
    )

    print(
        "[OK] GUI, Docker e modelli invariati"
    )


# ==========================================================
# 10. PROGRAMMA PRINCIPALE
# ==========================================================

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
        "AI PRIMETTA V3.3 — STEP 20"
    )

    print(
        "COMPLEMENTARITA' FULL LENGTH"
    )

    print("=" * 68)

    if args.audit:

        audit()

        print(
            "[OK] Audit senza modifiche"
        )

    else:

        run_scan()


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
