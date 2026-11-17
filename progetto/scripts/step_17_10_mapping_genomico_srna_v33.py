#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
================================================================
AI PRIMETTA V3.3
STEP 17.10 - MAPPING GENOMICO DI 5 SMALL RNA REALI
================================================================

OBIETTIVO
Localizzare i 5 small RNA PRJNA419326 nel genoma
Zea mays B73 NAM-5.0.

FUNZIONALITA
1. Controllo checksum ufficiale del genoma
2. Lettura dei 5 small RNA reali
3. Scansione del genoma completo, inclusi scaffold
4. Ricerca esatta su strand + e -
5. Coordinate genomiche 1-based
6. Identificazione dei loci multipli
7. Conteggi distinti per RNA e cromosoma
8. Salvataggio CSV e JSON
9. Controlli di coerenza finali

APPROCCIO
- Ricerca esatta: zero mismatch
- Nessun gap o bulge
- Nessun pacchetto Python esterno
- Lettura FASTA a blocchi
- Nessuna decompressione permanente del genoma

LIMITI SCIENTIFICI
Il mapping genomico esatto NON dimostra:
- biogenesi dello small RNA nel locus
- targeting post-trascrizionale
- regolazione genica
- meccanismi RdDM
- validazione sperimentale

PROTEZIONI
- GUI V3.2 invariata
- Motore V3.3 invariato
- Modelli ML/GNN invariati
- STEP precedenti invariati
- Nessun training
================================================================
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


# ================================================================
# CONFIGURAZIONE
# ================================================================

ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"

ASSEMBLY = "Zm-B73-REFERENCE-NAM-5.0"

GENOME = (
    BASE / "references"
    / (ASSEMBLY + ".fa.gz")
)

OFFICIAL_MD5 = (
    "8efe28f60e4415df7d7870e448f96c09"
)

EXPECTED_RNA = 5
EXPECTED_CHROMOSOMES = {f"chr{i}" for i in range(1, 11)}

# Dimensione dei blocchi di sequenza in caratteri.
CHUNK_SIZE = 4 * 1024 * 1024

# Salva un numero massimo di righe, ma continua
# a contare anche tutti gli eventuali loci successivi.
MAX_SAVED_LOCI = 100000

DNA = set("ACGT")

COMPLEMENT = str.maketrans(
    "ACGT",
    "TGCA"
)


# ================================================================
# FUNZIONI
# ================================================================

def section(number, title):

    print()
    print("=" * 72)
    print(f"[{number}] {title}")
    print("=" * 72)


def md5_file(path):

    h = hashlib.md5()

    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(4 * 1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def open_text(path):

    if str(path).lower().endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8"
        )

    return path.open(
        "r",
        encoding="utf-8"
    )


def read_fasta(path):

    header = None
    fragments = []

    with open_text(path) as handle:

        for raw in handle:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:

                    yield (
                        header.split()[0],
                        "".join(fragments)
                        .upper()
                        .replace("U", "T")
                    )

                header = line[1:]
                fragments = []

            elif header is not None:
                fragments.append(line)

    if header is not None:

        yield (
            header.split()[0],
            "".join(fragments)
            .upper()
            .replace("U", "T")
        )


def revcomp(sequence):

    return sequence.translate(
        COMPLEMENT
    )[::-1]


def latest_step16():

    matches = sorted(
        [
            p for p in RUNS.glob(
                "STEP16_5_SRNA_REALI_*"
            )
            if p.is_dir()
        ],
        reverse=True
    )

    if not matches:
        raise RuntimeError(
            "Nessuna cartella STEP 16 trovata"
        )

    return matches[0]


def build_patterns(rnas):

    patterns = defaultdict(list)

    for rna_id, sequence in rnas.items():

        patterns[sequence].append(
            (rna_id, "+")
        )

        patterns[revcomp(sequence)].append(
            (rna_id, "-")
        )

    return dict(patterns)


def scan_genome(
    genome_path,
    patterns,
    csv_writer,
    statistics
):
    """
    Ricerca esatta a blocchi.

    La sovrapposizione fra blocchi evita di perdere
    siti attraversanti il confine di lettura.

    Ogni posizione iniziale viene controllata una
    sola volta, anche per pattern ripetuti.
    """

    max_length = max(
        len(pattern) for pattern in patterns
    )

    overlap = max_length - 1

    current_chromosome = None

    buffer = ""
    buffer_start = 0
    chromosome_length = 0

    def search_buffer(sequence, start_offset, final=False):

        if not sequence:
            return 0

        if final:
            limit = len(sequence)
        else:
            # Le ultime max_length-1 posizioni
            # restano nel blocco successivo.
            limit = max(
                0,
                len(sequence) - overlap
            )

        for pattern, assignments in patterns.items():

            position = sequence.find(pattern)

            while position >= 0:

                if position >= limit:
                    break

                genomic_start = (
                    start_offset + position + 1
                )

                genomic_end = (
                    genomic_start + len(pattern) - 1
                )

                for rna_id, strand in assignments:

                    statistics["total_hits"] += 1

                    statistics["by_rna"][rna_id] += 1

                    statistics["by_strand"][strand] += 1

                    statistics["by_chromosome"][
                        current_chromosome
                    ] += 1

                    statistics["rna_chromosome"][
                        (rna_id, current_chromosome, strand)
                    ] += 1

                    if (
                        statistics["saved_hits"]
                        < MAX_SAVED_LOCI
                    ):

                        csv_writer.writerow({
                            "rna_id": rna_id,
                            "chromosome": current_chromosome,
                            "start_1based": genomic_start,
                            "end_1based": genomic_end,
                            "strand": strand,
                            "mismatches": 0,
                            "length_nt": len(pattern),
                            "reference": ASSEMBLY,
                            "mapping_type": "EXACT_GENOMIC",
                            "validation_status":
                                "COMPUTATIONAL_NOT_VALIDATED"
                        })

                        statistics["saved_hits"] += 1

                position = sequence.find(
                    pattern,
                    position + 1
                )

        return limit

    def process_chunk(final=False):

        nonlocal buffer, buffer_start

        limit = search_buffer(
            buffer,
            buffer_start,
            final=final
        )

        buffer = buffer[limit:]
        buffer_start += limit

    def finish_chromosome():

        nonlocal buffer
        nonlocal buffer_start
        nonlocal chromosome_length

        if current_chromosome is None:
            return

        process_chunk(final=True)

        statistics["chromosome_lengths"][
            current_chromosome
        ] = chromosome_length

        statistics["chromosomes_scanned"] += 1

        print(
            f"[OK] {current_chromosome} "
            f"| basi={chromosome_length:,} "
            f"| loci={statistics['by_chromosome'][current_chromosome]:,}",
            flush=True
        )

        buffer = ""
        buffer_start = 0
        chromosome_length = 0

    with gzip.open(
        genome_path,
        "rt",
        encoding="utf-8"
    ) as handle:

        for raw in handle:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                finish_chromosome()

                current_chromosome = (
                    line[1:].split()[0]
                )

                if (
                    current_chromosome
                    in statistics["chromosome_lengths"]
                ):
                    raise RuntimeError(
                        "Identificativo FASTA duplicato: "
                        + current_chromosome
                    )

                continue

            if current_chromosome is None:
                raise RuntimeError(
                    "FASTA privo di intestazione"
                )

            sequence = line.upper()

            chromosome_length += len(sequence)

            buffer += sequence

            if len(buffer) >= CHUNK_SIZE:
                process_chunk(final=False)

    finish_chromosome()


# ================================================================
# MAIN
# ================================================================

def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.10")
    print("MAPPING GENOMICO 5 SMALL RNA REALI")
    print("=" * 72)

    RUNS.mkdir(
        parents=True,
        exist_ok=True
    )

    step16 = latest_step16()

    RNA_FASTA = (
        step16 / "input_5_srna_reali.fa"
    )

    # ------------------------------------------------------------
    # 1. INTEGRITA
    # ------------------------------------------------------------

    section(
        1,
        "VERIFICA INPUT ORIGINALI"
    )

    for path in [GENOME, RNA_FASTA]:

        if not path.is_file():
            raise FileNotFoundError(path)

    local_md5 = md5_file(GENOME)

    print("Genoma:", GENOME)
    print("MD5:", local_md5)

    if local_md5 != OFFICIAL_MD5:

        raise RuntimeError(
            "Checksum del genoma non conforme. "
            "Mapping interrotto."
        )

    print("[OK] Checksum ufficiale conforme")

    # ------------------------------------------------------------
    # 2. SMALL RNA
    # ------------------------------------------------------------

    section(
        2,
        "CARICAMENTO SMALL RNA REALI"
    )

    rnas = dict(
        read_fasta(RNA_FASTA)
    )

    if len(rnas) != EXPECTED_RNA:

        raise RuntimeError(
            f"Attesi {EXPECTED_RNA} sRNA, "
            f"trovati {len(rnas)}"
        )

    for rna_id, sequence in rnas.items():

        if not 18 <= len(sequence) <= 30:
            raise RuntimeError(
                "Lunghezza RNA non valida: "
                + rna_id
            )

        if set(sequence) - DNA:
            raise RuntimeError(
                "Basi ambigue nell'RNA: "
                + rna_id
            )

        print(
            f"[RNA] {rna_id[:65]} "
            f"| {len(sequence)} nt"
        )

    # ------------------------------------------------------------
    # 3. PATTERN
    # ------------------------------------------------------------

    section(
        3,
        "PREPARAZIONE STRAND + E -"
    )

    patterns = build_patterns(rnas)

    print(
        "Sequenze di ricerca distinte:",
        len(patterns)
    )

    print(
        "Orientamenti:",
        "+ / -"
    )

    # ------------------------------------------------------------
    # 4. CARTELLA RISULTATI
    # ------------------------------------------------------------

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_10_MAPPING_GENOMICO_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    csv_file = (
        output / "LOCI_GENOMICI_SRNA.csv"
    )

    statistics = {
        "total_hits": 0,
        "saved_hits": 0,
        "chromosomes_scanned": 0,
        "chromosome_lengths": {},
        "by_rna": Counter(),
        "by_strand": Counter(),
        "by_chromosome": Counter(),
        "rna_chromosome": Counter()
    }

    # ------------------------------------------------------------
    # 5. SCANSIONE
    # ------------------------------------------------------------

    section(
        4,
        "SCANSIONE GENOMICA COMPLETA"
    )

    print(
        "[INFO] Genoma B73 NAM-5.0",
        flush=True
    )

    print(
        "[INFO] Ricerca esatta: 0 mismatch",
        flush=True
    )

    with csv_file.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rna_id",
                "chromosome",
                "start_1based",
                "end_1based",
                "strand",
                "mismatches",
                "length_nt",
                "reference",
                "mapping_type",
                "validation_status"
            ]
        )

        writer.writeheader()

        scan_genome(
            GENOME,
            patterns,
            writer,
            statistics
        )

    # ------------------------------------------------------------
    # 6. CONTROLLI
    # ------------------------------------------------------------

    section(
        5,
        "CONTROLLO COMPLETEZZA"
    )

    chromosome_set = set(
        statistics["chromosome_lengths"]
    )

    missing = sorted(
        EXPECTED_CHROMOSOMES - chromosome_set
    )

    print(
        "Sequenze scansionate:",
        statistics["chromosomes_scanned"]
    )

    print(
        "Cromosomi mancanti:",
        missing
    )

    print(
        "Loci esatti totali:",
        statistics["total_hits"]
    )

    print(
        "Loci salvati:",
        statistics["saved_hits"]
    )

    if missing:
        raise RuntimeError(
            "Mancano cromosomi canonici: "
            + ", ".join(missing)
        )

    truncated = (
        statistics["saved_hits"]
        < statistics["total_hits"]
    )

    if truncated:
        print(
            "[ATTENZIONE] Limite righe CSV raggiunto. "
            "I conteggi rimangono completi."
        )

    # ------------------------------------------------------------
    # 7. RISULTATI
    # ------------------------------------------------------------

    section(
        6,
        "RISULTATI PER SMALL RNA"
    )

    summary_rows = []

    for rna_id, sequence in rnas.items():

        count = statistics[
            "by_rna"
        ][rna_id]

        print()
        print("RNA:", rna_id)
        print("Lunghezza:", len(sequence))
        print("Loci esatti:", count)

        status = (
            "EXACT_MAPPING_FOUND"
            if count > 0
            else "NO_EXACT_MAPPING"
        )

        summary_rows.append({
            "rna_id": rna_id,
            "length_nt": len(sequence),
            "exact_loci": count,
            "mapping_status": status
        })

    # ------------------------------------------------------------
    # 8. CSV RIEPILOGO
    # ------------------------------------------------------------

    summary_csv = (
        output / "SINTESI_MAPPING_PER_RNA.csv"
    )

    with summary_csv.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rna_id",
                "length_nt",
                "exact_loci",
                "mapping_status"
            ]
        )

        writer.writeheader()
        writer.writerows(summary_rows)

    distribution_csv = (
        output / "DISTRIBUZIONE_RNA_CROMOSOMI.csv"
    )

    with distribution_csv.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.writer(handle)

        writer.writerow([
            "rna_id",
            "chromosome",
            "strand",
            "exact_loci"
        ])

        for key, count in sorted(
            statistics[
                "rna_chromosome"
            ].items()
        ):

            writer.writerow([
                key[0],
                key[1],
                key[2],
                count
            ])

    # ------------------------------------------------------------
    # 9. REPORT JSON
    # ------------------------------------------------------------

    section(
        7,
        "GENERAZIONE REPORT"
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.10",
        "stato": "MAPPING_ESATTO_COMPLETATO",
        "dataset": "PRJNA419326",
        "assembly": ASSEMBLY,
        "genoma": str(GENOME),
        "genoma_md5": local_md5,
        "checksum_ufficiale_conforme": True,
        "rna_fasta": str(RNA_FASTA),
        "rna_analizzati": len(rnas),
        "mismatch_massimi": 0,
        "strand": ["+", "-"],
        "sequenze_genomiche_scansionate":
            statistics["chromosomes_scanned"],
        "cromosomi_mancanti": missing,
        "basi_genomiche_totali": sum(
            statistics[
                "chromosome_lengths"
            ].values()
        ),
        "loci_esatti_totali":
            statistics["total_hits"],
        "loci_salvati":
            statistics["saved_hits"],
        "csv_troncato": truncated,
        "conteggi_per_rna": dict(
            statistics["by_rna"]
        ),
        "conteggi_per_strand": dict(
            statistics["by_strand"]
        ),
        "conteggi_per_cromosoma": dict(
            statistics["by_chromosome"]
        ),
        "riepilogo_rna": summary_rows,
        "limitazioni": [
            "Solo corrispondenze esatte",
            "Nessun mismatch",
            "Nessun gap o bulge",
            "Loci multipli non equivalgono a target",
            "Nessuna attribuzione funzionale",
            "Nessuna validazione sperimentale",
            "Nessuna annotazione GTF applicata"
        ],
        "gui_modificata": False,
        "motore_modificato": False,
        "training_eseguito": False
    }

    json_path = (
        output / "RAPPORTO_STEP17_10.json"
    )

    json_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    # ------------------------------------------------------------
    # 10. FINALE
    # ------------------------------------------------------------

    print()
    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.10 TERMINATO")
    print("=" * 72)

    print(
        "Sequenze genomiche:",
        statistics["chromosomes_scanned"]
    )

    print(
        "Loci esatti:",
        statistics["total_hits"]
    )

    print(
        "Loci salvati:",
        statistics["saved_hits"]
    )

    print()
    print("Report:", json_path)
    print("Loci:", csv_file)
    print("Riepilogo:", summary_csv)
    print("Distribuzione:", distribution_csv)

    print()
    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun training")
    print("[OK] Genoma originale preservato")

    print()
    print(
        "[NOTA] Mapping genomico esplorativo: "
        "nessuna funzione biologica attribuita."
    )


if __name__ == "__main__":
    main()