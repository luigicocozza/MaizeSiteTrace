#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.14 - ANALISI DEI CONTESTI GENOMICI DEGLI sRNA

Input reali:
- Genoma B73 NAM-5.0 verificato
- 646 loci dello STEP 17.11
- 5 sRNA PRJNA419326 dello STEP 16
- Annotazione TE dello STEP 17.13

Funzioni:
1. Verifica coerenza degli input
2. Estrazione delle sequenze genomiche con +/-100 nt
3. Controllo delle corrispondenze su entrambi gli strand
4. Entropia di Shannon e contenuto GC
5. Ricerca omopolimeri e ripetizioni periodiche
6. Analisi dei contesti identici tra loci
7. Riepilogo per small RNA
8. Esportazione CSV e JSON

Nessuna modifica alla GUI, al motore o ai dataset.
Nessun training.
Nessuna attribuzione funzionale automatica.
"""

from __future__ import annotations

import csv
import gzip
import json
import math

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"

GENOME = (
    BASE / "references"
    / "Zm-B73-REFERENCE-NAM-5.0.fa.gz"
)

FLANK = 100

COMPLEMENT = str.maketrans(
    "ACGTN", "TGCAN"
)


def section(number, title):
    print("\n" + "=" * 70)
    print(f"[{number}] {title}")
    print("=" * 70)


def latest(prefix, required_file):

    folders = sorted(
        (
            p for p in RUNS.glob(prefix + "*")
            if p.is_dir()
            and (p / required_file).is_file()
        ),
        reverse=True
    )

    if not folders:
        raise RuntimeError(
            "Cartella non trovata: " + prefix
        )

    return folders[0]


def read_csv(path):

    with path.open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, columns):

    with path.open(
        "w", encoding="utf-8", newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_rnas(path):

    records = {}
    header = None
    parts = []

    with path.open(
        "r", encoding="utf-8"
    ) as handle:

        for raw in handle:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    records[header] = "".join(parts).upper()

                header = line[1:].split()[0]
                parts = []

            elif header is not None:
                parts.append(line)

    if header is not None:
        records[header] = "".join(parts).upper()

    return {
        key: value.replace("U", "T")
        for key, value in records.items()
    }


def short_name(rna_id):

    for part in rna_id.split("|"):
        if part.startswith("feature="):
            return part.split("=", 1)[1]

    return rna_id.split("|")[0]


def shannon_entropy(seq):

    if not seq:
        return None

    counts = Counter(seq)
    length = len(seq)

    return -sum(
        (count / length) * math.log2(count / length)
        for count in counts.values()
    )


def longest_homopolymer(seq):

    if not seq:
        return 0

    best = 1
    current = 1

    for i in range(1, len(seq)):

        if seq[i] == seq[i - 1]:
            current += 1
            best = max(best, current)
        else:
            current = 1

    return best


def longest_periodic_run(seq, max_period=6):
    """
    Lunghezza massima di una ripetizione periodica
    consecutiva, con periodo da 1 a 6 nt.

    Misura descrittiva, non chiamata formale TE.
    """

    best = 0

    for period in range(1, max_period + 1):

        for start in range(len(seq) - period + 1):

            motif = seq[start:start + period]

            if len(motif) != period:
                continue

            if "N" in motif:
                continue

            end = start + period

            while (
                end + period <= len(seq)
                and seq[end:end + period] == motif
            ):
                end += period

            if end - start >= 2 * period:
                best = max(best, end - start)

    return best


def load_loci(path):

    rows = read_csv(path)

    if len(rows) != 646:
        raise RuntimeError(
            f"Attesi 646 loci: trovati {len(rows)}"
        )

    result = []
    seen = set()

    for i, row in enumerate(rows, start=1):

        start = int(row["start_1based"])
        end = int(row["end_1based"])

        if (
            start < 1
            or end < start
            or end - start + 1
            != int(row["length_nt"])
        ):
            raise RuntimeError(
                "Coordinate non valide: " + str(i)
            )

        if row["strand"] not in ("+", "-"):
            raise RuntimeError(
                "Strand non valido: " + str(i)
            )

        key = (
            row["rna_id"],
            row["chromosome"],
            start,
            end,
            row["strand"],
        )

        if key in seen:
            raise RuntimeError(
                "Locus duplicato: " + str(key)
            )

        seen.add(key)

        result.append({
            "locus_id": i,
            "rna_id": row["rna_id"],
            "srna": short_name(row["rna_id"]),
            "contig": row["chromosome"],
            "start": start,
            "end": end,
            "strand": row["strand"],
            "region": row.get("region", ""),
            "window_start": max(1, start - FLANK),
            "window_end": end + FLANK,
            "parts": [],
            "context": "",
        })

    return result


def extract_contexts(path, loci):
    """
    Una sola lettura del FASTA gzip.

    Estrae soltanto le finestre richieste,
    senza caricare i cromosomi completi in RAM.
    """

    grouped = defaultdict(list)

    for locus in loci:
        grouped[locus["contig"]].append(locus)

    for group in grouped.values():
        group.sort(key=lambda x: x["window_start"])

    contig = None
    position = 0
    targets = []
    first_active = 0
    seen_contigs = set()

    def finish_contig():

        if contig is None:
            return

        for locus in targets:
            locus["context"] = "".join(
                locus["parts"]
            )

        print(
            f"[OK] {contig} | "
            f"finestre={len(targets)}",
            flush=True
        )

    with gzip.open(
        path, "rt", encoding="utf-8"
    ) as handle:

        for raw in handle:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                finish_contig()

                contig = line[1:].split()[0]

                if contig in seen_contigs:
                    raise RuntimeError(
                        "Contig FASTA duplicato: "
                        + contig
                    )

                seen_contigs.add(contig)
                position = 0
                first_active = 0
                targets = grouped.get(contig, [])
                continue

            if contig is None:
                raise RuntimeError(
                    "FASTA senza intestazione"
                )

            sequence = line.upper()
            line_start = position + 1
            line_end = position + len(sequence)

            while (
                first_active < len(targets)
                and targets[first_active]["window_end"]
                < line_start
            ):
                first_active += 1

            j = first_active

            while (
                j < len(targets)
                and targets[j]["window_start"]
                <= line_end
            ):

                locus = targets[j]

                begin = max(
                    line_start,
                    locus["window_start"]
                )

                finish = min(
                    line_end,
                    locus["window_end"]
                )

                if begin <= finish:

                    locus["parts"].append(
                        sequence[
                            begin - line_start:
                            finish - line_start + 1
                        ]
                    )

                j += 1

            position = line_end

    finish_contig()

    missing = set(grouped) - seen_contigs

    if missing:
        raise RuntimeError(
            "Contig non presenti nel genoma: "
            + ", ".join(sorted(missing))
        )


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.14")
    print("ANALISI SEQUENZE GENOMICHE RIPETUTE")
    print("=" * 70)

    section(1, "INDIVIDUAZIONE INPUT")

    step16 = latest(
        "STEP16_5_SRNA_REALI_",
        "input_5_srna_reali.fa"
    )

    step11 = latest(
        "STEP17_11_ANNOTAZIONE_GTF_",
        "LOCI_ANNOTATI_STEP17_11.csv"
    )

    step13 = latest(
        "STEP17_13_ANNOTAZIONE_TE_",
        "LOCI_CON_ANNOTAZIONI_TE.csv"
    )

    rna_file = step16 / "input_5_srna_reali.fa"

    loci_file = (
        step11 / "LOCI_ANNOTATI_STEP17_11.csv"
    )

    te_file = (
        step13 / "LOCI_CON_ANNOTAZIONI_TE.csv"
    )

    if not GENOME.is_file():
        raise FileNotFoundError(GENOME)

    rnas = read_rnas(rna_file)
    loci = load_loci(loci_file)
    te_rows = read_csv(te_file)

    if len(rnas) != 5 or len(te_rows) != 646:
        raise RuntimeError(
            "Numero input non coerente"
        )

    te_lookup = {}

    for row in te_rows:

        key = (
            short_name(row["rna_id"]),
            row["chromosome"],
            int(row["start_1based"]),
            int(row["end_1based"]),
            row["strand"]
        )

        if key in te_lookup:
            raise RuntimeError(
                "Annotazione TE duplicata"
            )

        te_lookup[key] = (
            row["te_overlap"].strip().lower()
            == "true"
        )

    section(2, "ESTRAZIONE CONTESTI GENOMICI")

    extract_contexts(GENOME, loci)

    section(3, "CONTROLLO SEQUENZE E COMPLESSITA")

    results = []

    for locus in loci:

        original = rnas.get(locus["rna_id"])

        if original is None:
            raise RuntimeError(
                "RNA mancante: " + locus["rna_id"]
            )

        context = locus["context"]

        left = locus["start"] - locus["window_start"]
        length = len(original)

        site = context[left:left + length]

        expected = (
            original
            if locus["strand"] == "+"
            else revcomp(original)
        )

        if site != expected:
            raise RuntimeError(
                "Mismatch nel locus originale: "
                f"{locus['contig']}:"
                f"{locus['start']}-"
                f"{locus['end']} "
                f"RNA={locus['srna']}"
            )

        # Orienta il contesto rispetto allo small RNA.
        oriented = (
            context
            if locus["strand"] == "+"
            else revcomp(context)
        )

        valid = [
            base for base in oriented
            if base in "ACGT"
        ]

        valid_sequence = "".join(valid)

        gc = (
            100 * (
                valid_sequence.count("G")
                + valid_sequence.count("C")
            ) / len(valid_sequence)
            if valid_sequence else None
        )

        key = (
            locus["srna"],
            locus["contig"],
            locus["start"],
            locus["end"],
            locus["strand"],
        )

        if key not in te_lookup:
            raise RuntimeError(
                "Annotazione TE assente per un locus"
            )

        results.append({
            "locus_id": locus["locus_id"],
            "srna": locus["srna"],
            "rna_id": locus["rna_id"],
            "chromosome": locus["contig"],
            "start_1based": locus["start"],
            "end_1based": locus["end"],
            "strand": locus["strand"],
            "region": locus["region"],
            "te_overlap": te_lookup[key],
            "sequence_verified": True,
            "context_length": len(oriented),
            "expected_context_length": (
                length + 2 * FLANK
            ),
            "entropy_bits": round(
                shannon_entropy(valid_sequence), 4
            ),
            "gc_percent": round(gc, 3),
            "longest_homopolymer": (
                longest_homopolymer(oriented)
            ),
            "longest_periodic_run_1_6nt": (
                longest_periodic_run(oriented)
            ),
            "context_sequence": oriented,
        })

    if len(results) != 646:
        raise RuntimeError(
            "Sono stati persi loci"
        )

    section(4, "ANALISI CONTESTI IDENTICI")

    context_counts = Counter(
        (row["srna"], row["context_sequence"])
        for row in results
        if row["context_length"]
        == row["expected_context_length"]
    )

    for row in results:

        if (
            row["context_length"]
            == row["expected_context_length"]
        ):
            row["identical_context_copies"] = (
                context_counts[
                    (
                        row["srna"],
                        row["context_sequence"]
                    )
                ]
            )
        else:
            row["identical_context_copies"] = None

    section(5, "RIEPILOGO PER SMALL RNA")

    summary = []

    names = sorted({
        row["srna"] for row in results
    })

    for name in names:

        subset = [
            row for row in results
            if row["srna"] == name
        ]

        complete = [
            row for row in subset
            if row["identical_context_copies"]
            is not None
        ]

        entropy_values = [
            row["entropy_bits"] for row in subset
        ]

        repeated_contexts = sum(
            row["identical_context_copies"] > 1
            for row in complete
        )

        row = {
            "srna": name,
            "loci_totali": len(subset),
            "loci_con_TE": sum(
                x["te_overlap"] for x in subset
            ),
            "contesti_completi": len(complete),
            "contesti_identici_tra_loci": (
                repeated_contexts
            ),
            "contesti_distinti": len({
                x["context_sequence"]
                for x in complete
            }),
            "entropia_media": round(
                sum(entropy_values)
                / len(entropy_values),
                4
            ),
            "entropia_minima": min(entropy_values),
            "omopolimero_massimo": max(
                x["longest_homopolymer"]
                for x in subset
            ),
            "ripetizione_periodica_massima": max(
                x["longest_periodic_run_1_6nt"]
                for x in subset
            ),
        }

        summary.append(row)

        print()
        print("RNA:", name)
        print("Loci:", row["loci_totali"])
        print("Con TE:", row["loci_con_TE"])

        print(
            "Contesti identici tra loci:",
            repeated_contexts
        )

        print(
            "Contesti distinti:",
            row["contesti_distinti"]
        )

        print(
            "Entropia media:",
            row["entropia_media"]
        )

        print(
            "Entropia minima:",
            row["entropia_minima"]
        )

        print(
            "Ripetizione periodica massima:",
            row["ripetizione_periodica_massima"]
        )

    section(6, "SALVATAGGIO")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_14_SEQUENZE_RIPETUTE_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    write_csv(
        output / "CONTESTI_GENOMICI_SRNA.csv",
        results,
        list(results[0].keys())
    )

    write_csv(
        output / "SINTESI_CONTESTI_PER_SRNA.csv",
        summary,
        list(summary[0].keys())
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.14",
        "stato": "ANALISI_COMPLETATA",
        "genoma": str(GENOME),
        "fonte_loci": str(loci_file),
        "fonte_te": str(te_file),
        "finestra_flank_nt": FLANK,
        "loci_verificati": len(results),
        "rna_verificati": len(names),
        "mismatch_sui_loci_originali": 0,
        "summary": summary,
        "limitazioni": [
            "Solo contesti locali di +/-100 nt",
            "La ripetizione esatta del contesto "
            "non equivale a classificazione TE",
            "Entropia e ripetizioni periodiche "
            "sono descrittori di sequenza",
            "I contesti parziali non sono utilizzati "
            "nel confronto di identita",
            "Le sequenze genomiche non dimostrano "
            "il locus effettivo di biogenesi",
            "Nessuna validazione biologica",
        ],
        "gui_modified": False,
        "engine_modified": False,
        "training_performed": False,
    }

    report_path = (
        output / "RAPPORTO_STEP17_14.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print()
    print("=" * 70)
    print("STEP 17.14 TERMINATO")
    print("=" * 70)

    print("Loci verificati:", len(results))
    print("Small RNA:", len(names))
    print("Errori sequenza:", 0)

    print("\nRisultati:", output)
    print("Report:", report_path)

    print("\n[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun training")
    print("[OK] File originali preservati")


if __name__ == "__main__":
    main()