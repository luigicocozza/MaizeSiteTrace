#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 31

AUDIT RISORSE PER LA CLASSIFICAZIONE BIOLOGICA
DELLE SMALL RNA.

--audit:
    Verifica STEP30, FASTA originale e
    presenza delle risorse genomiche.

Non modifica GUI, Docker, dati o modelli.
"""

import argparse
import csv
import hashlib
import json
import os
import sys

from collections import Counter
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

EXPECTED = {
    "sRNA_00000037": 25,
    "sRNA_00000435": 24
}

SKIP = {
    ".git",
    ".venv",
    "venv",
    "node_modules",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    "runs"
}

SEQ_SUFFIXES = (
    ".fa",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz"
)


# ============================================================
# 2. CHECKSUM SHA256
# ============================================================

def checksum(path):

    h = hashlib.sha256()

    with path.open("rb") as stream:

        for block in iter(
            lambda: stream.read(1024 * 1024),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# ============================================================
# 3. LETTURA FASTA
# ============================================================

def read_fasta(path):

    records = {}

    current = None
    chunks = []

    for raw in path.read_text(
        encoding="utf-8-sig"
    ).splitlines():

        row = raw.strip()

        if not row:
            continue

        if row.startswith(">"):

            if current is not None:

                if (
                    not chunks
                    or current in records
                ):

                    raise ValueError(
                        "Record FASTA non valido: "
                        + current
                    )

                records[current] = (
                    "".join(chunks)
                    .upper()
                    .replace("T", "U")
                )

            current = row[1:].split()[0]

            chunks = []

            if (
                not current
                or current in records
            ):

                raise ValueError(
                    "ID FASTA duplicato/vuoto"
                )

        else:

            if current is None:

                raise ValueError(
                    "Sequenza senza intestazione FASTA"
                )

            chunks.append(row)

    if current is not None:

        if (
            not chunks
            or current in records
        ):

            raise ValueError(
                "Ultimo record FASTA non valido"
            )

        records[current] = (
            "".join(chunks)
            .upper()
            .replace("T", "U")
        )

    if not records:

        raise ValueError(
            "FASTA vuoto"
        )

    return records


# ============================================================
# 4. VERIFICA REGISTRO STEP30
# ============================================================

def verify_existing_results():

    reports = sorted(
        RUNS.glob(
            "STEP30_REGISTRO_SCIENTIFICO_*/"
            "RAPPORTO_STEP30.json"
        )
    )

    if not reports:

        raise FileNotFoundError(
            "STEP30 mancante"
        )

    report_path = reports[-1]

    r30 = json.loads(
        report_path.read_text(
            encoding="utf-8"
        )
    )

    if (
        r30.get("step") != 30
        or r30.get("candidate_count") != 36
    ):

        raise ValueError(
            "Manifest STEP30 inatteso"
        )

    # ========================================================
    # 5. RECUPERO REGISTRO SCIENTIFICO
    # ========================================================

    outputs = [
        ROOT / name
        for name in r30["outputs"]
    ]

    register = next(
        (
            p
            for p in outputs
            if p.name == "REGISTRO_36_CANDIDATI.tsv"
        ),
        None
    )

    if (
        register is None
        or not register.is_file()
        or not register.resolve().is_relative_to(
            RUNS.resolve()
        )
    ):

        raise ValueError(
            "Registro STEP30 non trovato"
        )

    with register.open(
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        records = list(
            csv.DictReader(
                stream,
                delimiter="\t"
            )
        )

    counts = Counter(
        row["srna_id"]
        for row in records
    )

    if (
        len(records) != 36
        or counts != Counter({
            "sRNA_00000037": 17,
            "sRNA_00000435": 19
        })
    ):

        raise ValueError(
            "Conteggi registro STEP30 incoerenti"
        )

    if any(
        float(
            r["BH_q_74_exploratory"]
        ) < 0.05

        for r in records
    ):

        raise ValueError(
            "Significativita' non coerente "
            "con il rapporto STEP30"
        )

    # ========================================================
    # 6. RICOSTRUZIONE PROVENIENZA STEP22
    # ========================================================

    r29_path = (
        ROOT
        / r30["source_reports"]["step29"]
    )

    if (
        not r29_path.is_file()
        or not r29_path.resolve().is_relative_to(
            RUNS.resolve()
        )
    ):

        raise ValueError(
            "Fonte STEP29 non valida"
        )

    r29 = json.loads(
        r29_path.read_text(
            encoding="utf-8"
        )
    )

    r28_path = (
        ROOT
        / r29["source_step28"]
    )

    r28 = json.loads(
        r28_path.read_text(
            encoding="utf-8"
        )
    )

    r24_path = (
        ROOT
        / r28["source_step24"]
    )

    r24 = json.loads(
        r24_path.read_text(
            encoding="utf-8"
        )
    )

    r23_path = (
        ROOT
        / r24["source_step23b"]
    )

    r23 = json.loads(
        r23_path.read_text(
            encoding="utf-8"
        )
    )

    step22_path = Path(
        r23["source_manifest"]
    )

    if not step22_path.is_absolute():

        step22_path = (
            ROOT
            / step22_path
        )

    if (
        not step22_path.is_file()
        or not step22_path.resolve().is_relative_to(
            RUNS.resolve()
        )
    ):

        raise ValueError(
            "MANIFEST_STEP22.json originale "
            "non trovato"
        )

    manifest = json.loads(
        step22_path.read_text(
            encoding="utf-8"
        )
    )

    # ========================================================
    # 7. VERIFICA SEQUENZE ORIGINALI
    # ========================================================

    srna_path = (
        step22_path.parent
        / "SRNA_STANDARD_MAX25.fa"
    )

    if not srna_path.is_file():

        raise FileNotFoundError(
            str(srna_path)
        )

    if checksum(srna_path) != (
        manifest["files_sha256"]["standard"]
    ):

        raise ValueError(
            "SHA256 small RNA discordante "
            "dal manifest STEP22"
        )

    srnas = read_fasta(
        srna_path
    )

    lengths = {
        key: len(seq)
        for key, seq in srnas.items()
    }

    if lengths != EXPECTED:

        raise ValueError(
            "Identita' o lunghezza delle "
            "small RNA inattese"
        )

    if any(
        set(seq) - set("ACGU")
        for seq in srnas.values()
    ):

        raise ValueError(
            "Small RNA con basi non canoniche"
        )

    return (
        report_path,
        srna_path,
        srnas
    )


# ============================================================
# 8. CLASSIFICAZIONE PRELIMINARE DEI FILE
# ============================================================

def categorize(filename):

    name = filename.lower()

    is_fasta = name.endswith(
        SEQ_SUFFIXES
    )

    is_cdna = any(
        token in name

        for token in (
            "cdna",
            "transcript",
            "mrna",
            "peptide",
            "protein",
            ".pep.",
            "srna",
            "smallrna"
        )
    )

    # cDNA: trascritti, non genoma completo

    if is_fasta and "cdna" in name:

        return (
            "trascrittoma_cDNA_NON_genoma"
        )

    # Potenziali archivi di miRNA note

    if (
        any(
            token in name

            for token in (
                "mirbase",
                "pmiren",
                "mirna"
            )
        )

        and

        (
            is_fasta
            or name.endswith(
                (
                    ".gff",
                    ".gff3",
                    ".bed",
                    ".tsv"
                )
            )
        )
    ):

        return (
            "annotazione_miRNA_da_verificare"
        )

    # Possibile riferimento genomico B73

    if (
        is_fasta
        and not is_cdna
        and any(

            token in name

            for token in (
                "genome",
                "genomic",
                "primary_assembly",
                "toplevel",
                ".dna.",
                "chromosome",
                "reference"
            )
        )
    ):

        return (
            "genoma_DNA_da_verificare"
        )

    # Letture originali di sequenziamento

    if name.endswith(
        (
            ".fastq",
            ".fastq.gz",
            ".fq",
            ".fq.gz"
        )
    ):

        return (
            "letture_FASTQ_da_verificare"
        )

    # Annotazioni genomiche

    if name.endswith(
        (
            ".gff",
            ".gff3",
            ".gtf",
            ".gff.gz",
            ".gff3.gz",
            ".gtf.gz"
        )
    ):

        return (
            "annotazione_genica_da_verificare"
        )

    # Possibili allineamenti o loci

    if name.endswith(
        (
            ".bam",
            ".cram",
            ".sam",
            ".bed",
            ".bed.gz"
        )
    ):

        return (
            "allineamento_o_loci_da_verificare"
        )

    return None


# ============================================================
# 9. INVENTARIO RISORSE PRESENTI
# ============================================================

def inspect_resources():

    categories = (

        "genoma_DNA_da_verificare",

        "trascrittoma_cDNA_NON_genoma",

        "annotazione_genica_da_verificare",

        "annotazione_miRNA_da_verificare",

        "letture_FASTQ_da_verificare",

        "allineamento_o_loci_da_verificare"
    )

    buckets = {
        key: []
        for key in categories
    }

    examined = 0

    for folder, dirs, files in os.walk(
        ROOT
    ):

        # Evita di ripercorrere output storici,
        # ambienti virtuali e dipendenze software.

        dirs[:] = [
            directory

            for directory in dirs

            if directory not in SKIP
        ]

        for filename in files:

            examined += 1

            category = categorize(
                filename
            )

            if category is not None:

                path = (
                    Path(folder)
                    / filename
                )

                buckets[category].append(
                    str(
                        path.relative_to(ROOT)
                    )
                )

    return examined, buckets


# ============================================================
# 10. PROGRAMMA PRINCIPALE
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

    print("=" * 67)

    print(
        "AI PRIMETTA V3.3 | STEP 31"
    )

    print(
        "AUDIT CLASSE BIOLOGICA SMALL RNA"
    )

    print("=" * 67)

    # ========================================================
    # 11. VERIFICA DEI RISULTATI PRECEDENTI
    # ========================================================

    report, fasta_path, srnas = (
        verify_existing_results()
    )

    print(
        "[PASS] STEP30: "
        "36 siti candidati; "
        "nessun BH q<0.05"
    )

    print(
        "[PASS] FASTA STEP22 e SHA256 verificati"
    )

    # ========================================================
    # 12. DESCRIZIONE DELLE DUE SMALL RNA
    # ========================================================

    for sid, sequence in sorted(
        srnas.items()
    ):

        print(

            f"[sRNA] {sid} "

            f"| {len(sequence)} nt "

            f"| 5'={sequence[0]} "

            f"| sequenza={sequence}"
        )

        print(
            "[CLASSE] non determinabile "
            "dalla sola lunghezza "
            "o complementarita'"
        )

    # ========================================================
    # 13. INVENTARIO RISORSE
    # ========================================================

    examined, resources = (
        inspect_resources()
    )

    print(

        "[INVENTARIO] File esaminati "
        f"(esclusi runs e ambienti): {examined}"
    )

    for category, paths in (
        resources.items()
    ):

        print(
            f"[RISORSE] {category}: {len(paths)}"
        )

        for path in sorted(paths)[:8]:

            print(
                "  -",
                path
            )

        if len(paths) > 8:

            print(
                f"  ... altri {len(paths)-8} file"
            )

    # ========================================================
    # 14. RISULTATO FINALE
    # ========================================================

    print(
        "[NOTA] Le risorse sono individuate "
        "in base ai nomi: la loro presenza "
        "non ne dimostra l'adeguatezza"
    )

    print(
        "[DECISIONE] Il prossimo passo "
        "dipende dalla disponibilita' "
        "di genoma B73, GFF3 e dati sRNA-seq"
    )

    print(
        "[OK] Audit in sola lettura; "
        "GUI, Docker e modelli invariati"
    )


# ============================================================
# 15. AVVIO
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            f"[ERRORE] {type(exc).__name__}: {exc}",
            file=sys.stderr
        )

        sys.exit(1)