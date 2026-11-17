
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 18

Verifica indipendente della predizione target
delle cinque small RNA reali PRJNA419326.

Modalità:
    --audit : controlla input, riferimento e algoritmo
    --scan  : scansione esaustiva per 0–3 mismatch

Non modifica GUI V3.2, GUI V3.3, Docker,
modelli ML/GNN o dati originali.
"""

import argparse
import csv
import gzip
import json
import sys

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


# ==========================================================
# 1. CONFIGURAZIONE
# ==========================================================

ROOT = Path(__file__).resolve().parent.parent

V33 = ROOT / "target_prediction_v33_vscode"

RUNS = V33 / "runs"

EXPECTED_SRNA = (
    "sRNA_00000002",
    "sRNA_00000040",
    "sRNA_00000007",
    "sRNA_00000435",
    "sRNA_00000037",
)

MAX_MM = 3

MAX_ROWS_PER_RNA = 2000

DNA = set("ACGT")

COMPLEMENT = str.maketrans(
    "ACGT",
    "TGCA"
)


# ==========================================================
# 2. LETTURA FASTA
# ==========================================================

def fasta_records(path):

    opener = (
        gzip.open
        if path.suffix == ".gz"
        else open
    )

    with opener(
        path,
        "rt",
        encoding="utf-8",
        errors="strict"
    ) as stream:

        header = None
        chunks = []

        for line in stream:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:

                    yield (
                        header,
                        "".join(chunks)
                        .upper()
                        .replace("U", "T")
                    )

                header = line[1:].split()[0]
                chunks = []

            else:

                if header is None:

                    raise ValueError(
                        "FASTA non valido"
                    )

                chunks.append(line)

        if header is not None:

            yield (
                header,
                "".join(chunks)
                .upper()
                .replace("U", "T")
            )


# ==========================================================
# 3. RICERCA INPUT STEP16
# ==========================================================

def select_input():

    exact = (
        RUNS
        / "STEP16_5_SRNA_REALI_20261021T181110_648766Z"
        / "input_5_srna_reali.fa"
    )

    if exact.is_file():
        return exact

    found = sorted(
        RUNS.glob(
            "STEP16_5_SRNA_REALI_*/"
            "input_5_srna_reali.fa"
        )
    )

    if len(found) != 1:

        raise RuntimeError(
            "FASTA delle cinque small RNA "
            "non individuato univocamente. "
            f"File trovati: {len(found)}"
        )

    return found[0]


# ==========================================================
# 4. RIFERIMENTO TRASCRITTOMICO B73
# ==========================================================

def select_reference():

    folder = V33 / "references"

    matching = sorted(
        folder.glob(
            "*Zm-B73-REFERENCE-NAM-5.0*"
            "cdna*.fa.gz"
        )
    )

    if not matching:

        matching = sorted(
            folder.glob(
                "*cdna*.fa.gz"
            )
        )

    if len(matching) != 1:

        raise RuntimeError(
            "Riferimento cDNA non univoco: "
            f"{len(matching)} file in {folder}"
        )

    return matching[0]


# ==========================================================
# 5. ANNOTAZIONE TRASCRITTO → GENE
# ==========================================================

def find_mapping(reference):

    direct = reference.name.replace(
        ".cdna.all.fa.gz",
        ".tx2gene.tsv"
    )

    paths = [
        reference.parent / direct
    ]

    paths += sorted(
        reference.parent.glob(
            "*.tx2gene.tsv"
        )
    )

    paths = list(
        dict.fromkeys(
            p for p in paths
            if p.is_file()
        )
    )

    return (
        paths[0]
        if len(paths) == 1
        else None
    )


def parse_mapping(path):

    if path is None:
        return {}

    mapping = {}

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        reader = csv.reader(
            stream,
            delimiter="\t"
        )

        for i, row in enumerate(reader):

            if len(row) < 2:
                continue

            tx = row[0].strip()
            gene = row[1].strip()

            if (
                i == 0
                and tx.lower() in (
                    "transcript",
                    "transcript_id",
                    "tx_id"
                )
            ):
                continue

            if tx and gene:

                mapping[tx] = gene

    return mapping


# ==========================================================
# 6. CINQUE SMALL RNA REALI
# ==========================================================

def load_queries(path):

    records = list(
        fasta_records(path)
    )

    if len(records) != 5:

        raise RuntimeError(
            "Numero small RNA inatteso: "
            f"{len(records)}"
        )

    queries = []
    seen = set()

    for header, sequence in records:

        identifiers = [
            ident
            for ident in EXPECTED_SRNA
            if ident in header
        ]

        if len(identifiers) != 1:

            raise RuntimeError(
                "ID non riconosciuto: "
                + header
            )

        ident = identifiers[0]

        if ident in seen:

            raise RuntimeError(
                "ID duplicato: "
                + ident
            )

        seen.add(ident)

        if (
            not 20 <= len(sequence) <= 30
            or not set(sequence) <= DNA
        ):

            raise RuntimeError(
                "Sequenza non valida: "
                + ident
            )

        queries.append(
            (ident, sequence)
        )

    if seen != set(EXPECTED_SRNA):

        raise RuntimeError(
            "Insieme di small RNA "
            "diverso da quello previsto"
        )

    return queries


# ==========================================================
# 7. REVERSE COMPLEMENT
# ==========================================================

def revcomp(sequence):

    return sequence.translate(
        COMPLEMENT
    )[::-1]


# ==========================================================
# 8. PARTIZIONE ESAUSTIVA DEI SEED
# ==========================================================

def blocks(
    sequence,
    max_mm=MAX_MM
):

    count = max_mm + 1
    n = len(sequence)

    return [
        (
            i * n // count,

            sequence[
                i * n // count:
                (i + 1) * n // count
            ]
        )

        for i in range(count)
    ]


# ==========================================================
# 9. ALGORITMO HAMMING ESAUSTIVO
# ==========================================================

def ungapped_hits(
    transcript,
    target_query,
    max_mm=MAX_MM
):

    """
    Ricerca tutti gli allineamenti ungapped
    con distanza Hamming <= max_mm.

    Dividendo la sequenza in max_mm+1
    segmenti disgiunti, almeno uno
    deve essere identico in ogni allineamento
    che soddisfa la soglia.
    """

    length = len(target_query)

    if len(transcript) < length:
        return []

    starts = set()

    for offset, seed in blocks(
        target_query,
        max_mm
    ):

        index = transcript.find(seed)

        while index != -1:

            start = index - offset

            if (
                0 <= start
                <= len(transcript) - length
            ):

                starts.add(start)

            index = transcript.find(
                seed,
                index + 1
            )

    result = []

    for start in sorted(starts):

        segment = transcript[
            start:start + length
        ]

        if not set(segment) <= DNA:
            continue

        mismatches = []

        for j, (a, b) in enumerate(
            zip(
                segment,
                target_query
            ),
            1
        ):

            if a != b:

                mismatches.append(j)

                if len(mismatches) > max_mm:
                    break

        if len(mismatches) <= max_mm:

            result.append(
                (
                    start,
                    len(mismatches),
                    segment,
                    mismatches
                )
            )

    return result


# ==========================================================
# 10. TEST SINTETICI DELL'ALGORITMO
# ==========================================================

def self_test(queries):

    for ident, sequence in queries:

        target = revcomp(sequence)

        partitions = blocks(target)

        assert len(partitions) == 4

        assert "".join(
            part
            for _, part in partitions
        ) == target

        for changes in range(4):

            modified = list(target)

            for offset, _ in partitions[:changes]:

                modified[offset] = next(
                    nucleotide
                    for nucleotide in "ACGT"
                    if nucleotide != modified[offset]
                )

            test_sequence = (
                "CCCCC"
                + "".join(modified)
                + "GGGGG"
            )

            hits = ungapped_hits(
                test_sequence,
                target
            )

            assert any(
                start == 5
                and mm == changes
                for start, mm, _, _ in hits
            ), ident

    return (
        "Test 0/1/2/3 mismatch "
        "superati per tutte le cinque sRNA"
    )


# ==========================================================
# 11. AUDIT PREVENTIVO
# ==========================================================

def audit():

    if not (
        V33 / "app_v33.py"
    ).is_file():

        raise RuntimeError(
            "app_v33.py assente"
        )

    query_file = select_input()

    reference = select_reference()

    queries = load_queries(
        query_file
    )

    mapping_file = find_mapping(
        reference
    )

    print(
        "[REPOSITORY]",
        ROOT
    )

    print(
        "[INPUT 5 sRNA]",
        query_file.relative_to(ROOT)
    )

    for ident, sequence in queries:

        print(
            f"  {ident}: "
            f"{len(sequence)} nt "
            "| sequenza verificata"
        )

    print(
        "[RIFERIMENTO B73]",
        reference.relative_to(ROOT)
    )

    print(
        "[DIMENSIONE COMPRESSA]",
        reference.stat().st_size,
        "byte"
    )

    print(
        "[TX2GENE]",
        (
            mapping_file.relative_to(ROOT)
            if mapping_file
            else "assente/non univoco"
        )
    )

    print(
        "[SELF-TEST]",
        self_test(queries)
    )

    print(
        "[METODO] Reverse-complement, "
        "Hamming 0–3, senza gap"
    )

    print(
        "[LIMITE] G:U, indel ed energia "
        "non valutati"
    )

    return (
        query_file,
        reference,
        mapping_file,
        queries
    )


# ==========================================================
# 12. SCANSIONE COMPLETA
# ==========================================================

def scan():

    (
        input_file,
        reference,
        mapping_file,
        queries
    ) = audit()

    mapping = parse_mapping(
        mapping_file
    )

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    output = RUNS / (
        "STEP18_ESATTO_HAMMING3_"
        + stamp
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    rows_file = (
        output
        / "SITI_COMPLEMENTARI_0_3_MM.tsv"
    )

    stats = {
        ident: Counter()
        for ident, _ in queries
    }

    rows_written = Counter()

    seen_transcripts = set()

    n_transcripts = 0
    n_nt = 0

    columns = (
        "srna_id",
        "srna_length",
        "transcript_id",
        "gene_id",
        "start_1based",
        "end_1based",
        "mismatches",
        "mismatch_positions_1based",
        "srna_5to3_DNA",
        "complement_5to3_DNA",
        "target_site_5to3_DNA"
    )

    prepared = [
        (
            ident,
            sequence,
            revcomp(sequence)
        )

        for ident, sequence in queries
    ]

    with rows_file.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as stream:

        writer = csv.writer(
            stream,
            delimiter="\t"
        )

        writer.writerow(columns)

        for transcript, sequence in fasta_records(
            reference
        ):

            if transcript in seen_transcripts:

                raise RuntimeError(
                    "Trascritto duplicato: "
                    + transcript
                )

            seen_transcripts.add(
                transcript
            )

            n_transcripts += 1

            n_nt += len(sequence)

            for ident, srna, target in prepared:

                hits = ungapped_hits(
                    sequence,
                    target
                )

                for (
                    start,
                    mm,
                    segment,
                    positions
                ) in hits:

                    stats[ident][
                        str(mm)
                    ] += 1

                    if (
                        rows_written[ident]
                        >= MAX_ROWS_PER_RNA
                    ):
                        continue

                    writer.writerow(
                        (
                            ident,
                            len(srna),
                            transcript,
                            mapping.get(
                                transcript,
                                ""
                            ),
                            start + 1,
                            start + len(srna),
                            mm,
                            ",".join(
                                map(
                                    str,
                                    positions
                                )
                            ),
                            srna,
                            target,
                            segment
                        )
                    )

                    rows_written[ident] += 1

            if n_transcripts % 10000 == 0:

                print(
                    "[SCANSIONE]",
                    n_transcripts,
                    "trascritti esaminati",
                    flush=True
                )

    # ======================================================
    # 13. RIEPILOGO RISULTATI
    # ======================================================

    summary = {}

    for ident, sequence in queries:

        counts = {
            str(d): stats[ident][str(d)]
            for d in range(MAX_MM + 1)
        }

        total = sum(
            counts.values()
        )

        summary[ident] = {

            "length_nt": len(sequence),

            "sites_by_mismatches": counts,

            "sites_total_0_3": total,

            "sites_exported": (
                rows_written[ident]
            ),

            "tsv_truncated": (
                total > rows_written[ident]
            )
        }

        print(
            f"[RISULTATO] {ident}: "
            f"{total} siti | "
            f"mismatch={counts}"
        )

    # ======================================================
    # 14. RAPPORTO DI PROVENIENZA
    # ======================================================

    report = {

        "step": "18",

        "method": (
            "Seed esaustivi per distanza "
            "Hamming <=3, reverse complement, "
            "senza gap"
        ),

        "status": "COMPLETED",

        "input_srna": str(
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

        "transcripts_processed": n_transcripts,

        "reference_bases_processed": n_nt,

        "sample": "Zea mays B73 cDNA",

        "query_count": len(queries),

        "result": summary,

        "results_tsv": str(
            rows_file.relative_to(ROOT)
        ),

        "limitations": [

            "Non considera gap e G:U wobble",

            "Non calcola energia o accessibilità",

            "Non dimostra interazioni biologiche",

            "Eventuali zero match riguardano "
            "soltanto gli allineamenti ungapped "
            "con al massimo 3 mismatch"
        ],

        "v32_gui_modified": False,

        "v33_gui_modified": False,

        "training_performed": False
    }

    report_file = (
        output / "RAPPORTO_STEP18.json"
    )

    report_file.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    print(
        "[REPORT]",
        report_file
    )

    print(
        "[RISULTATI]",
        rows_file
    )

    print(
        "[OK] Nessuna GUI modificata"
    )

    print(
        "[OK] Nessun training"
    )


# ==========================================================
# 15. PROGRAMMA PRINCIPALE
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

    print("=" * 69)

    print(
        "AI PRIMETTA V3.3 | STEP 18"
    )

    print(
        "VERIFICA ESAUSTIVA DEI TARGET"
    )

    print("=" * 69)

    if args.audit:

        audit()

        print(
            "[OK] Audit completato "
            "senza modifiche"
        )

    else:

        scan()


if __name__ == "__main__":

    try:

        main()

    except Exception as error:

        print(
            "[ERRORE]",
            type(error).__name__,
            str(error),
            file=sys.stderr
        )

        sys.exit(1)
