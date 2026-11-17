#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 33 - MAPPATURA GENOMICA ESATTA DELLE SMALL RNA

MODALITA':

--audit
    Verifica provenienza scientifica, sequenze
    small RNA e intestazione del genoma.
    Non genera output nel progetto.

--run
    Cerca corrispondenze esatte su entrambi
    i filamenti dell'intero riferimento B73.

NON MODIFICA:
    GUI V3.2 / V3.3
    Docker
    Modelli ML/GNN
    Dati originali
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE PRINCIPALE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

GENOME = (
    ROOT
    / "target_prediction_v33_vscode"
    / "references"
    / "Zm-B73-REFERENCE-NAM-5.0.fa.gz"
)

EXPECTED = {
    "sRNA_00000037": 25,
    "sRNA_00000435": 24
}

CHUNK_BASES = 2_000_000

DNA_LETTERS = set(
    "ACGTRYSWKMBDHVN"
)


# ============================================================
# 2. FUNZIONE DI CONTROLLO
# ============================================================

def verify(condition, message):

    if not condition:
        raise ValueError(message)


# ============================================================
# 3. PROTEZIONE DEI PERCORSI
# ============================================================

def safe_file(value):

    p = Path(value)

    if not p.is_absolute():
        p = ROOT / p

    p = p.resolve()

    verify(
        p.is_file()
        and p.is_relative_to(
            RUNS.resolve()
        ),
        f"File sorgente non valido: {p}"
    )

    return p


# ============================================================
# 4. CHECKSUM SHA256
# ============================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for chunk in iter(
            lambda: f.read(1 << 20),
            b""
        ):

            h.update(chunk)

    return h.hexdigest()


# ============================================================
# 5. LETTURA DEL FASTA SMALL RNA
# ============================================================

def fasta_small(path):

    result = {}

    current = None
    seq = []

    for line in path.read_text(
        encoding="utf-8-sig"
    ).splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if current is not None:

                verify(
                    bool(seq)
                    and current not in result,
                    "FASTA duplicato o vuoto"
                )

                result[current] = (
                    "".join(seq)
                    .upper()
                    .replace("T", "U")
                )

            current = (
                line[1:].split()[0]
            )

            seq = []

            verify(
                bool(current)
                and current not in result,
                "ID FASTA ripetuto"
            )

        else:

            verify(
                current is not None,
                "Sequenza senza intestazione"
            )

            seq.append(line)

    if current is not None:

        verify(
            bool(seq)
            and current not in result,
            "Ultimo record FASTA non valido"
        )

        result[current] = (
            "".join(seq)
            .upper()
            .replace("T", "U")
        )

    # ========================================================
    # 6. CONTROLLO IDENTITA' DELLE DUE SMALL RNA
    # ========================================================

    lengths = {
        k: len(v)
        for k, v in result.items()
    }

    verify(
        lengths == EXPECTED,
        "Identita' e lunghezza sRNA errate"
    )

    verify(
        all(
            set(v) <= set("ACGU")
            for v in result.values()
        ),
        "Basi sRNA non canoniche"
    )

    return result


# ============================================================
# 7. RECUPERO DATI SCIENTIFICI APPROVATI
# ============================================================

def sources():

    candidates = sorted(
        RUNS.glob(
            "STEP30_REGISTRO_SCIENTIFICO_*/"
            "RAPPORTO_STEP30.json"
        )
    )

    verify(
        bool(candidates),
        "Report STEP30 assente"
    )

    p30 = safe_file(
        candidates[-1]
    )

    r30 = json.loads(
        p30.read_text(
            encoding="utf-8"
        )
    )

    verify(
        str(
            r30.get("step")
        ) == "30"
        and r30.get("candidate_count") == 36,
        "Registro STEP30 non coerente"
    )

    # ========================================================
    # 8. VERIFICA REGISTRO 36 CANDIDATI
    # ========================================================

    registers = [

        safe_file(p)

        for p in r30["outputs"]

        if Path(p).name
        == "REGISTRO_36_CANDIDATI.tsv"
    ]

    verify(
        len(registers) == 1,
        "Registro STEP30 mancante"
    )

    with registers[0].open(
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

    verify(
        len(rows) == 36
        and counts == Counter({

            "sRNA_00000037": 17,
            "sRNA_00000435": 19

        }),
        "Registro STEP30 incompleto"
    )

    # ========================================================
    # 9. RICOSTRUZIONE CATENA DI PROVENIENZA
    # ========================================================

    p29 = safe_file(
        r30["source_reports"]["step29"]
    )

    r29 = json.loads(
        p29.read_text(
            encoding="utf-8"
        )
    )

    p28 = safe_file(
        r29["source_step28"]
    )

    r28 = json.loads(
        p28.read_text(
            encoding="utf-8"
        )
    )

    p24 = safe_file(
        r28["source_step24"]
    )

    r24 = json.loads(
        p24.read_text(
            encoding="utf-8"
        )
    )

    p23 = safe_file(
        r24["source_step23b"]
    )

    r23 = json.loads(
        p23.read_text(
            encoding="utf-8"
        )
    )

    p22 = safe_file(
        r23["source_manifest"]
    )

    r22 = json.loads(
        p22.read_text(
            encoding="utf-8"
        )
    )

    # ========================================================
    # 10. SEQUENZE SMALL RNA ORIGINALI
    # ========================================================

    query_path = safe_file(
        p22.parent
        / "SRNA_STANDARD_MAX25.fa"
    )

    verify(
        sha256(query_path)
        == r22["files_sha256"]["standard"],
        "Checksum STEP22 small RNA errato"
    )

    sequences = fasta_small(
        query_path
    )

    return (
        p30,
        p22,
        query_path,
        sequences
    )


# ============================================================
# 11. COMPLEMENTO INVERSO DEL DNA
# ============================================================

def reverse_complement(dna):

    translation = str.maketrans(
        "ACGT",
        "TGCA"
    )

    return dna.translate(
        translation
    )[::-1]


# ============================================================
# 12. AUDIT INIZIALE DEL GENOMA
# ============================================================

def reference_audit(path):

    verify(
        path.is_file(),
        f"Genoma assente: {path}"
    )

    with path.open("rb") as raw:

        verify(
            raw.read(2) == b"\x1f\x8b",
            "Genoma non compresso in gzip"
        )

    with gzip.open(
        path,
        "rt",
        encoding="ascii"
    ) as f:

        header = f.readline().strip()

        seq = f.readline().strip().upper()

    verify(
        header.startswith(">")
        and bool(seq)
        and set(seq) <= DNA_LETTERS,
        "Genoma: intestazione o prima "
        "sequenza invalida"
    )

    return header[1:].split()[0]


# ============================================================
# 13. MAPPATURA GENOMICA ESATTA
# ============================================================

def scan_genome(
    genome,
    queries,
    writer
):

    """
    Ricerca streaming su FASTA genomico.

    Analizza:
      - sequenza diretta
      - sequenza reverse complement

    Gestisce le corrispondenze sovrapposte
    e i confini tra blocchi di lettura.
    """

    patterns = []

    for sid, rna in sorted(
        queries.items()
    ):

        dna = rna.replace(
            "U",
            "T"
        )

        patterns.extend((

            (
                sid,
                "+",
                dna
            ),

            (
                sid,
                "-",
                reverse_complement(dna)
            )
        ))

    max_len = max(
        len(item[2])
        for item in patterns
    )

    counts = Counter()

    examples = {
        sid: []
        for sid in queries
    }

    contigs = []

    chrom = None

    total = 0

    tail = ""

    chunks = []

    pending_bases = 0

    prev_headers = set()

    n_hits_chrom = 0

    # ========================================================
    # 14. ANALISI DEL SINGOLO BLOCCO
    # ========================================================

    def feed(block):

        nonlocal total
        nonlocal tail
        nonlocal n_hits_chrom

        combined = tail + block

        start0 = (
            total - len(tail)
        )

        matches = []

        for (
            sid,
            strand,
            dna
        ) in patterns:

            start = 0

            while True:

                i = combined.find(
                    dna,
                    start
                )

                if i < 0:
                    break

                # Evita duplicati dovuti ai
                # nucleotidi conservati nel
                # blocco precedente.

                if (
                    start0 + i + len(dna)
                    <= total
                ):

                    start = i + 1
                    continue

                begin = (
                    start0 + i + 1
                )

                end = (
                    begin + len(dna) - 1
                )

                matches.append((
                    begin,
                    end,
                    sid,
                    strand,
                    dna
                ))

                # Incrementando di una sola
                # posizione troviamo anche
                # match sovrapposti.

                start = i + 1

        # ====================================================
        # 15. ESPORTAZIONE DELLE CORRISPONDENZE
        # ====================================================

        for (
            begin,
            end,
            sid,
            strand,
            dna
        ) in sorted(matches):

            writer.writerow([

                sid,

                chrom,

                strand,

                begin,

                end,

                len(dna),

                "exact_full_length",

                "B73_NAM5"
            ])

            counts[
                (sid, strand)
            ] += 1

            n_hits_chrom += 1

            if len(
                examples[sid]
            ) < 6:

                examples[sid].append({

                    "contig":
                        chrom,

                    "strand":
                        strand,

                    "start":
                        begin,

                    "end":
                        end
                })

        # ====================================================
        # 16. AGGIORNAMENTO DEL BLOCCO
        # ====================================================

        total += len(block)

        tail = combined[
            -(max_len - 1):
        ]

    # ========================================================
    # 17. SVUOTAMENTO BUFFER
    # ========================================================

    def flush():

        nonlocal chunks
        nonlocal pending_bases

        if chunks:

            feed(
                "".join(chunks)
            )

            chunks = []

            pending_bases = 0

    # ========================================================
    # 18. CHIUSURA DEL CROMOSOMA
    # ========================================================

    def close_chrom():

        nonlocal chrom
        nonlocal total
        nonlocal tail
        nonlocal n_hits_chrom

        if chrom is None:
            return

        flush()

        verify(
            total > 0,
            f"Contig senza sequenza: {chrom}"
        )

        contigs.append({

            "name":
                chrom,

            "length_bp":
                total,

            "hits":
                n_hits_chrom
        })

        if (
            len(contigs) <= 15
            or n_hits_chrom
        ):

            print(

                f"[CONTIG] {chrom} "

                f"| bp={total} "

                f"| loci={n_hits_chrom}",

                flush=True
            )

        chrom = None

        total = 0

        tail = ""

        n_hits_chrom = 0

    # ========================================================
    # 19. SCANSIONE DEL FASTA COMPRESSO
    # ========================================================

    with gzip.open(
        genome,
        "rt",
        encoding="ascii"
    ) as f:

        for row in f:

            row = row.strip()

            if not row:
                continue

            if row.startswith(">"):

                close_chrom()

                chrom = (
                    row[1:].split()[0]
                )

                verify(
                    bool(chrom)
                    and chrom not in prev_headers,
                    "Identificativo contig "
                    "duplicato o vuoto"
                )

                prev_headers.add(
                    chrom
                )

            else:

                verify(
                    chrom is not None,
                    "Basi prima della prima "
                    "intestazione FASTA"
                )

                block = row.upper()

                verify(
                    set(block) <= DNA_LETTERS,
                    f"Caratteri inattesi "
                    f"nel contig {chrom}"
                )

                chunks.append(
                    block
                )

                pending_bases += len(
                    block
                )

                if (
                    pending_bases
                    >= CHUNK_BASES
                ):

                    flush()

        close_chrom()

    verify(
        bool(contigs),
        "Nessun contig nel riferimento"
    )

    return (
        counts,
        contigs,
        examples
    )


# ============================================================
# 20. PROGRAMMA PRINCIPALE
# ============================================================

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
        "AI PRIMETTA V3.3 | STEP 33 "
        "- LOCI GENOMICI ESATTI"
    )

    print("=" * 68)

    # ========================================================
    # 21. VERIFICA SCIENTIFICA INIZIALE
    # ========================================================

    (
        p30,
        p22,
        small_file,
        queries
    ) = sources()

    header = reference_audit(
        GENOME
    )

    print(
        "[PASS] Registro STEP30 e SHA256 "
        "sequenze STEP22 verificati"
    )

    print(

        "[GENOMA]",

        GENOME.relative_to(ROOT),

        f"| primo contig={header}"
    )

    for sid, seq in sorted(
        queries.items()
    ):

        print(

            f"[sRNA] {sid} "

            f"| {len(seq)} nt "

            f"| DNA={seq.replace('U', 'T')}"
        )

    # ========================================================
    # 22. MODALITA' AUDIT
    # ========================================================

    if args.audit:

        print(
            "[OK] Audit sola lettura. "
            "Genoma completo NON ancora scandito."
        )

        return

    # ========================================================
    # 23. NUOVA CARTELLA OUTPUT
    # ========================================================

    tag = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    outdir = RUNS / (
        f"STEP33_LOCI_GENOMICI_{tag}"
    )

    outdir.mkdir(
        parents=True,
        exist_ok=False
    )

    partial = (
        outdir
        / "LOCI_ESATTI_2_SRNA.tsv.partial"
    )

    output = (
        outdir
        / "LOCI_ESATTI_2_SRNA.tsv"
    )

    # ========================================================
    # 24. ESECUZIONE RICERCA COMPLETA
    # ========================================================

    try:

        with partial.open(
            "w",
            encoding="utf-8",
            newline=""
        ) as stream:

            writer = csv.writer(
                stream,
                delimiter="\t"
            )

            writer.writerow([

                "srna_id",

                "contig",

                "strand",

                "start_1based",

                "end_1based",

                "length_nt",

                "match_type",

                "reference"
            ])

            (
                counts,
                contigs,
                examples
            ) = scan_genome(

                GENOME,

                queries,

                writer
            )

        # Il file viene marcato come definitivo
        # solo dopo il completamento della scansione.

        os.replace(
            partial,
            output
        )

        # ====================================================
        # 25. CHECKSUM DEL GENOMA UTILIZZATO
        # ====================================================

        genome_hash = sha256(
            GENOME
        )

        # ====================================================
        # 26. RIEPILOGO PER SMALL RNA
        # ====================================================

        stats = {}

        for sid in sorted(
            queries
        ):

            plus = counts[
                (sid, "+")
            ]

            minus = counts[
                (sid, "-")
            ]

            stats[sid] = {

                "plus":
                    plus,

                "minus":
                    minus,

                "total_exact_loci":
                    plus + minus,

                "examples":
                    examples[sid],

                "classification":
                    "unclassified_srna"
            }

            print(

                f"[RISULTATO] {sid}: "

                f"loci_plus={plus}, "

                f"loci_minus={minus}, "

                f"totale={plus + minus}"
            )

        # ====================================================
        # 27. RAPPORTO SCIENTIFICO JSON
        # ====================================================

        summary = {

            "step":
                33,

            "status":
                "EXACT_GENOMIC_MAPPING_NOT_BIOLOGICAL_CLASSIFICATION",

            "genome":
                str(
                    GENOME.relative_to(ROOT)
                ),

            "genome_compressed_sha256":
                genome_hash,

            "srna_fasta":
                str(
                    small_file.relative_to(ROOT)
                ),

            "srna_fasta_sha256":
                sha256(small_file),

            "source_step30":
                str(
                    p30.relative_to(ROOT)
                ),

            "source_step30_sha256":
                sha256(p30),

            "contigs_scanned":
                len(contigs),

            "reference_bases_scanned":
                sum(
                    item["length_bp"]
                    for item in contigs
                ),

            "contig_lengths":
                contigs,

            "by_srna":
                stats,

            "output_tsv":
                str(
                    output.relative_to(ROOT)
                ),

            "outputs_sha256": {

                "loci_tsv":
                    sha256(output)
            },

            "limits": [

                "Only exact full-length genomic "
                "matches on both strands",

                "No mismatches, indels, precursor "
                "evidence or RNA processing evidence",

                "Multi-mapping regions may be "
                "repetitive; a match is not "
                "proof of transcription",

                "No matches on B73 does not "
                "exclude sequence origin from "
                "other maize genotypes",

                "No GTF/TE overlap or miRNA/siRNA "
                "biological classification yet"
            ],

            "gui_modified":
                False,

            "docker_modified":
                False,

            "training_performed":
                False
        }

        report = (
            outdir
            / "RAPPORTO_STEP33.json"
        )

        report.write_text(

            json.dumps(
                summary,
                ensure_ascii=False,
                indent=2
            ) + "\n",

            encoding="utf-8"
        )

        # ====================================================
        # 28. OUTPUT FINALE
        # ====================================================

        print(

            "[GENOMA]",

            len(contigs),

            "contig;",

            summary[
                "reference_bases_scanned"
            ],

            "bp analizzati"
        )

        print(
            "[TSV]",
            output
        )

        print(
            "[REPORT]",
            report
        )

        print(
            "[OK] Scansione completata, "
            "nessuna modifica a GUI/Docker/modelli"
        )

    # ========================================================
    # 29. PROTEZIONE IN CASO DI INTERRUZIONE
    # ========================================================

    except Exception:

        warning = (
            outdir
            / "RUN_INCOMPLETA_NON_USARE.txt"
        )

        warning.write_text(

            "Scansione interrotta: "
            "non utilizzare eventuali "
            "file .partial.\n",

            encoding="utf-8"
        )

        raise


# ============================================================
# 30. AVVIO DEL PROGRAMMA
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