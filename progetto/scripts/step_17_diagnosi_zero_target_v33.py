#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
=======================================================================
AI PRIMETTA V3.3
STEP 17 - DIAGNOSI ZERO TARGET RNA-GENE SU DATI REALI
=======================================================================

OBIETTIVO
---------
Diagnosticare il risultato dello STEP 16:

- 5 small RNA reali PRJNA419326
- 72539 trascritti analizzati
- 0 siti candidati
- 0 geni target
- massimo 2 mismatch

VERIFICHE
---------
1. Integrita del progetto
2. Identificazione dell'ultimo STEP 16
3. Individuazione dei FASTA reali
4. Verifica trascrittoma B73
5. Verifica mapping trascritto-gene
6. Controllo positivo e negativo sintetico
7. Matching indipendente RNA-trascritto
8. Analisi progressiva delle soglie di mismatch
9. Verifica orientamento diretto/reverse-complement
10. Audit statico del codice del motore
11. Generazione report JSON, CSV e Markdown

SICUREZZA
---------
- GUI V3.2 non modificata
- Modelli ML/GNN non modificati
- Dataset originali non modificati
- Risultati STEP 16 non sovrascritti
- Nessun training
- Nessun dato biologico inventato
- Nessun container modificato

LIMITI
------
Il matching indipendente senza gap non equivale alla
predizione biologica del motore V3.3.

Questo STEP non dimostra l'esistenza di target biologici,
ne' richiama direttamente analyze() nel container Docker.

=======================================================================
"""

from __future__ import annotations

import ast
import csv
import gzip
import hashlib
import json
import os
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


# ================================================================
# 1. CONFIGURAZIONE
# ================================================================

ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

SCRIPTS_DIR = ROOT / "scripts"

RUNS_DIR = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

STEP16_SCRIPT = (
    SCRIPTS_DIR
    / "step_16_test_5_srna_reali_v33.py"
)

STEP17_NAME = "step_17_diagnosi_zero_target_v33.py"

TRANSCRIPTOME_NAME = (
    "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz"
)

# Analisi diagnostica, non predizione biologica.
MISMATCH_THRESHOLDS = [0, 1, 2, 3, 4, 5]

# Il controllo e' intenzionalmente limitato per evitare
# tempi eccessivi nella VM VirtualBox.
MAX_TRANSCRIPTS_TEST = 200

MAX_BASES_PER_TRANSCRIPT = 3000

MAX_EXAMPLES_PER_SRNA = 5

EXPECTED_SRNA = 5

EXPECTED_TRANSCRIPTS = 72539

DNA_BASES = set("ACGT")

VALID_FASTA_BASES = set("ACGTN")

IGNORE_DIRS = {
    ".git",
    ".venv",
    "venv",
    "__pycache__",
    ".pytest_cache",
    "node_modules",
}

COMPLEMENT = str.maketrans(
    "ACGT",
    "TGCA",
)


# ================================================================
# 2. FUNZIONI GENERALI
# ================================================================

def print_section(number, title):
    print()
    print("=" * 72)
    print(f"[{number}] {title}")
    print("=" * 72)


def sha256_file(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def open_text(path):

    if str(path).lower().endswith(".gz"):

        return gzip.open(
            path,
            "rt",
            encoding="utf-8",
            errors="replace",
        )

    return path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    )


def iter_fasta(path):
    """
    Lettura streaming FASTA.

    Converte U in T per confronti fra RNA e DNA.
    """

    header = None
    fragments = []

    with open_text(path) as handle:

        for raw_line in handle:

            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:

                    yield (
                        header,
                        "".join(fragments)
                        .upper()
                        .replace("U", "T"),
                    )

                header = line[1:].strip()
                fragments = []

            elif header is not None:
                fragments.append(line)

    if header is not None:

        yield (
            header,
            "".join(fragments)
            .upper()
            .replace("U", "T"),
        )


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def sequence_entropy(sequence):
    """
    Entropia di Shannon in bit.
    """

    if not sequence:
        return 0.0

    from math import log2

    counts = Counter(sequence)
    length = len(sequence)

    return -sum(
        (count / length)
        * log2(count / length)
        for count in counts.values()
    )


def find_project_files():

    results = []

    for current, dirs, filenames in os.walk(ROOT):

        dirs[:] = [
            directory
            for directory in dirs
            if directory not in IGNORE_DIRS
            and not directory.startswith("STEP17_")
        ]

        for filename in filenames:

            results.append(
                Path(current) / filename
            )

    return results


def is_fasta(path):

    return bool(
        re.search(
            r"\.(fa|fasta|fna)(\.gz)?$",
            path.name,
            re.IGNORECASE,
        )
    )


def safe_relative(path):

    try:
        return str(path.relative_to(ROOT))

    except ValueError:
        return str(path)


# ================================================================
# 3. MATCHING DIAGNOSTICO
# ================================================================

def best_hamming_match(
    query,
    transcript,
    max_mismatches=5,
):
    """
    Cerca il miglior allineamento senza gap.

    Restituisce:
    - numero mismatch
    - posizione iniziale 0-based

    Restituisce None se non trova corrispondenze
    entro la soglia.

    Non gestisce:
    - bulge
    - gap
    - pairing G:U
    - energia libera
    - accessibilita RNA
    - degradoma
    """

    query_length = len(query)

    transcript_length = len(transcript)

    if (
        query_length == 0
        or transcript_length < query_length
    ):
        return None

    best = max_mismatches + 1
    best_start = None

    for start in range(
        transcript_length - query_length + 1
    ):

        mismatch = 0

        for i in range(query_length):

            if query[i] != transcript[start + i]:

                mismatch += 1

                if mismatch >= best:
                    break

        if mismatch < best:

            best = mismatch
            best_start = start

            if best == 0:
                break

    if best <= max_mismatches:

        return {
            "mismatches": best,
            "start_0_based": best_start,
        }

    return None


def positive_negative_controls():

    synthetic_rna = "ATGCCGTAACGTAGCTGACCTAGC"

    direct_transcript = (
        "GGGG"
        + synthetic_rna
        + "CCCC"
    )

    rc = reverse_complement(synthetic_rna)

    reverse_transcript = (
        "GGGG"
        + rc
        + "CCCC"
    )

    negative_transcript = "T" * 120

    direct = best_hamming_match(
        synthetic_rna,
        direct_transcript,
        2,
    )

    reverse = best_hamming_match(
        rc,
        reverse_transcript,
        2,
    )

    negative = best_hamming_match(
        synthetic_rna,
        negative_transcript,
        2,
    )

    passed = (
        direct is not None
        and direct["mismatches"] == 0
        and reverse is not None
        and reverse["mismatches"] == 0
        and negative is None
    )

    return {
        "test_superato": passed,
        "match_diretto": direct,
        "match_reverse_complement": reverse,
        "match_negativo": negative,
        "natura": "CONTROLLO_SINTETICO_SOFTWARE",
        "validazione_biologica": False,
    }


# ================================================================
# 4. INDIVIDUAZIONE STEP 16
# ================================================================

def find_latest_step16():

    if not RUNS_DIR.exists():
        return None

    candidates = sorted(
        [
            p
            for p in RUNS_DIR.iterdir()
            if p.is_dir()
            and p.name.startswith(
                "STEP16_5_SRNA_REALI_"
            )
        ],
        key=lambda p: p.name,
        reverse=True,
    )

    return candidates[0] if candidates else None


# ================================================================
# 5. SELEZIONE FASTA SMALL RNA
# ================================================================

def find_five_srna(files, step16):

    candidates = []

    if step16:

        candidates.extend(
            [
                p
                for p in step16.rglob("*")
                if p.is_file()
                and is_fasta(p)
            ]
        )

    candidates.extend(
        [
            p
            for p in files
            if is_fasta(p)
            and (
                "srna" in p.name.lower()
                or "small_rna" in p.name.lower()
            )
        ]
    )

    seen = set()

    for path in candidates:

        if path in seen:
            continue

        seen.add(path)

        try:

            # Evita FASTA grandi estranei ai 5 sRNA.
            if path.stat().st_size > 5_000_000:
                continue

            records = list(iter_fasta(path))

            if len(records) != EXPECTED_SRNA:
                continue

            if not all(
                18 <= len(seq) <= 35
                for _, seq in records
            ):
                continue

            return path, records

        except (OSError, EOFError, ValueError):
            continue

    return None, []


# ================================================================
# 6. AUDIT CODICE DEL MOTORE
# ================================================================

def inspect_python_file(path):

    result = {
        "file": safe_relative(path),
        "sha256": sha256_file(path),
        "funzioni": [],
        "parametri_analyze": [],
        "parole_chiave": [],
    }

    source = path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    try:
        tree = ast.parse(source)

    except SyntaxError as error:

        result["errore_sintassi"] = str(error)

        return result

    for node in ast.walk(tree):

        if isinstance(
            node,
            (
                ast.FunctionDef,
                ast.AsyncFunctionDef,
            ),
        ):

            result["funzioni"].append(node.name)

            if node.name == "analyze":

                result["parametri_analyze"] = [
                    arg.arg
                    for arg in (
                        node.args.posonlyargs
                        + node.args.args
                        + node.args.kwonlyargs
                    )
                ]

    keywords = [
        "max_mismatches",
        "maximum_sites",
        "entropy_min",
        "reverse_complement",
        "revcomp",
        "rna_fasta",
        "transcripts",
        "mapping_file",
        "docker",
        "analyze",
    ]

    result["parole_chiave"] = [
        keyword
        for keyword in keywords
        if keyword in source
    ]

    return result


# ================================================================
# 7. PROCEDURA PRINCIPALE
# ================================================================

def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3")
    print("STEP 17 - DIAGNOSI ZERO TARGET")
    print("RNA-GENE | DATI REALI PRJNA419326")
    print("=" * 72)

    if not ROOT.exists():

        raise SystemExit(
            "[ERRORE] Cartella progetto non trovata:\n"
            + str(ROOT)
        )

    RUNS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output_dir = (
        RUNS_DIR
        / f"STEP17_DIAGNOSI_ZERO_TARGET_{timestamp}"
    )

    output_dir.mkdir(exist_ok=False)

    report = {
        "progetto": "AI PRIMETTA",
        "versione": "V3.3",
        "step": 17,
        "timestamp_utc": timestamp,
        "stato": "DIAGNOSI_IN_CORSO",
        "trascrittoma": {},
        "mapping": {},
        "small_rna": {},
        "controlli_sintetici": {},
        "matching_reale": [],
        "audit_codice": [],
        "avvertenze": [],
        "limiti": [
            "Matching Hamming senza gap",
            "Analisi limitata a un campione di trascritti",
            "Nessuna validazione biologica",
            "Motore Docker V3.3 non eseguito",
            "GNN non applicata",
        ],
    }

    # ============================================================
    # CONTROLLO 1 - INVENTARIO
    # ============================================================

    print_section(1, "INVENTARIO PROGETTO")

    files = find_project_files()

    print("[INFO] File individuati:", len(files))

    step16 = find_latest_step16()

    report["step16_run"] = (
        str(step16)
        if step16
        else None
    )

    if step16:

        print("[OK] STEP 16 trovato:")
        print(step16)

    else:

        print("[ATTENZIONE] STEP 16 non trovato")

        report["avvertenze"].append(
            "Directory STEP 16 non individuata"
        )

    # ============================================================
    # CONTROLLO 2 - SMALL RNA
    # ============================================================

    print_section(2, "SMALL RNA REALI")

    srna_path, srna_records = find_five_srna(
        files,
        step16,
    )

    if srna_path:

        print("[OK] FASTA small RNA:")
        print(srna_path)

        records_info = []

        for header, sequence in srna_records:

            info = {
                "id": header,
                "length_nt": len(sequence),
                "entropy": round(
                    sequence_entropy(sequence),
                    6,
                ),
                "valid_acgt": (
                    set(sequence) <= DNA_BASES
                ),
            }

            records_info.append(info)

            print(
                f"RNA: {header[:55]} "
                f"| nt={len(sequence)} "
                f"| H={info['entropy']:.4f}"
            )

        report["small_rna"] = {
            "file": str(srna_path),
            "sha256": sha256_file(srna_path),
            "numero": len(srna_records),
            "record": records_info,
        }

    else:

        print(
            "[ATTENZIONE] FASTA contenente esattamente "
            "5 small RNA non trovato"
        )

        report["avvertenze"].append(
            "FASTA dei 5 small RNA non identificato: "
            "il test di matching reale non potra "
            "essere eseguito"
        )

    # ============================================================
    # CONTROLLO 3 - TRASCRITTOMA
    # ============================================================

    print_section(3, "TRASCRITTOMA REALE B73")

    transcriptomes = [
        p
        for p in files
        if p.name == TRANSCRIPTOME_NAME
    ]

    transcript_sample = []

    if transcriptomes:

        path = transcriptomes[0]

        print("[INFO] Lettura:", path)

        number = 0
        min_length = None
        max_length = 0

        invalid_bases = Counter()

        for header, sequence in iter_fasta(path):

            number += 1

            sequence_length = len(sequence)

            min_length = (
                sequence_length
                if min_length is None
                else min(min_length, sequence_length)
            )

            max_length = max(
                max_length,
                sequence_length,
            )

            invalid_bases.update(
                base
                for base in sequence
                if base not in VALID_FASTA_BASES
            )

            if (
                len(transcript_sample)
                < MAX_TRANSCRIPTS_TEST
            ):

                transcript_sample.append(
                    (
                        header,
                        sequence[
                            :MAX_BASES_PER_TRANSCRIPT
                        ],
                    )
                )

        report["trascrittoma"] = {
            "file": str(path),
            "sha256": sha256_file(path),
            "numero_trascritti": number,
            "numero_atteso_step16": EXPECTED_TRANSCRIPTS,
            "numero_conforme_step16": (
                number == EXPECTED_TRANSCRIPTS
            ),
            "lunghezza_min": min_length,
            "lunghezza_max": max_length,
            "basi_non_acgtn": dict(invalid_bases),
            "trascritti_campione": len(
                transcript_sample
            ),
            "basi_massime_campione": (
                MAX_BASES_PER_TRANSCRIPT
            ),
        }

        print("[OK] Trascritti:", number)
        print("[INFO] Lunghezza minima:", min_length)
        print("[INFO] Lunghezza massima:", max_length)

        if number != EXPECTED_TRANSCRIPTS:

            print(
                "[ATTENZIONE] Numero diverso "
                "dal risultato STEP 16"
            )

            report["avvertenze"].append(
                "Numero trascritti differente "
                "dai 72539 dello STEP 16"
            )

    else:

        print(
            "[ATTENZIONE] Trascrittoma originale "
            "non individuato"
        )

        report["avvertenze"].append(
            "Trascrittoma B73 non trovato"
        )

    # ============================================================
    # CONTROLLO 4 - MAPPING
    # ============================================================

    print_section(4, "MAPPING TRASCRITTO-GENE")

    mapping_candidates = [
        p
        for p in files
        if p.name.lower() in {
            "tx2gene.tsv",
            "tx2gene.tsv.gz",
        }
    ]

    if mapping_candidates:

        mapping_path = mapping_candidates[0]

        mapping_rows = 0
        examples = []

        with open_text(mapping_path) as handle:

            for line in handle:

                if not line.strip():
                    continue

                mapping_rows += 1

                if len(examples) < 3:
                    examples.append(
                        line.strip()[:200]
                    )

        report["mapping"] = {
            "file": str(mapping_path),
            "sha256": sha256_file(mapping_path),
            "righe_non_vuote": mapping_rows,
            "esempi": examples,
        }

        print("[OK] Mapping trovato:")
        print(mapping_path)
        print("[INFO] Righe:", mapping_rows)

    else:

        print(
            "[ATTENZIONE] tx2gene.tsv non trovato"
        )

        report["avvertenze"].append(
            "Mapping tx2gene non individuato"
        )

    # ============================================================
    # CONTROLLO 5 - TEST SINTETICI
    # ============================================================

    print_section(
        5,
        "CONTROLLI POSITIVO E NEGATIVO",
    )

    synthetic = positive_negative_controls()

    report["controlli_sintetici"] = synthetic

    print(
        "Match diretto:",
        synthetic["match_diretto"],
    )

    print(
        "Reverse complement:",
        synthetic["match_reverse_complement"],
    )

    print(
        "Controllo negativo:",
        synthetic["match_negativo"],
    )

    if synthetic["test_superato"]:

        print(
            "[OK] Controlli software sintetici superati"
        )

    else:

        print(
            "[ERRORE] Controlli sintetici non superati"
        )

        report["avvertenze"].append(
            "Matching diagnostico indipendente "
            "non supera i controlli sintetici"
        )

    # ============================================================
    # CONTROLLO 6 - MATCHING SU DATI REALI
    # ============================================================

    print_section(
        6,
        "MATCHING INDIPENDENTE SU DATI REALI",
    )

    csv_rows = []

    if srna_records and transcript_sample:

        print(
            "[INFO] Trascritti campionati:",
            len(transcript_sample),
        )

        print(
            "[INFO] Mismatch:",
            MISMATCH_THRESHOLDS,
        )

        for index, (
            srna_header,
            srna_sequence,
        ) in enumerate(
            srna_records,
            start=1,
        ):

            print()
            print(
                f"[RNA {index}/{len(srna_records)}] "
                f"{srna_header[:55]}"
            )

            reverse_sequence = reverse_complement(
                srna_sequence
            )

            counts_direct = Counter()
            counts_reverse = Counter()

            examples = []

            for (
                tx_header,
                tx_sequence,
            ) in transcript_sample:

                direct = best_hamming_match(
                    srna_sequence,
                    tx_sequence,
                    max(MISMATCH_THRESHOLDS),
                )

                reverse = best_hamming_match(
                    reverse_sequence,
                    tx_sequence,
                    max(MISMATCH_THRESHOLDS),
                )

                for threshold in MISMATCH_THRESHOLDS:

                    if (
                        direct is not None
                        and direct["mismatches"]
                        <= threshold
                    ):
                        counts_direct[threshold] += 1

                    if (
                        reverse is not None
                        and reverse["mismatches"]
                        <= threshold
                    ):
                        counts_reverse[threshold] += 1

                if (
                    direct is not None
                    or reverse is not None
                ) and len(examples) < MAX_EXAMPLES_PER_SRNA:

                    examples.append(
                        {
                            "trascritto": tx_header.split()[0],
                            "mismatch_diretto": (
                                direct["mismatches"]
                                if direct
                                else None
                            ),
                            "mismatch_reverse": (
                                reverse["mismatches"]
                                if reverse
                                else None
                            ),
                            "posizione_diretta": (
                                direct["start_0_based"]
                                if direct
                                else None
                            ),
                            "posizione_reverse": (
                                reverse["start_0_based"]
                                if reverse
                                else None
                            ),
                        }
                    )

            record = {
                "srna_id": srna_header,
                "lunghezza_nt": len(srna_sequence),
                "trascritti_testati": len(
                    transcript_sample
                ),
                "hits_diretti": dict(counts_direct),
                "hits_reverse": dict(counts_reverse),
                "esempi": examples,
            }

            report["matching_reale"].append(record)

            for threshold in MISMATCH_THRESHOLDS:

                direct_hits = counts_direct[threshold]
                reverse_hits = counts_reverse[threshold]

                print(
                    f"  mismatch <= {threshold}: "
                    f"diretti={direct_hits} "
                    f"reverse={reverse_hits}"
                )

                csv_rows.append(
                    {
                        "srna": srna_header,
                        "lunghezza_nt": len(
                            srna_sequence
                        ),
                        "trascritti_testati": len(
                            transcript_sample
                        ),
                        "max_mismatch": threshold,
                        "hit_diretti": direct_hits,
                        "hit_reverse": reverse_hits,
                    }
                )

    else:

        print(
            "[ATTENZIONE] Matching reale non eseguito"
        )

        report["avvertenze"].append(
            "Matching reale non eseguito: "
            "mancano small RNA o trascritti"
        )

    # ============================================================
    # CONTROLLO 7 - CODICE MOTORE
    # ============================================================

    print_section(
        7,
        "AUDIT STATICO DEL CODICE MOTORE",
    )

    python_candidates = []

    for path in files:

        if path.suffix != ".py":
            continue

        try:

            if path.stat().st_size > 2_000_000:
                continue

            text = path.read_text(
                encoding="utf-8",
                errors="replace",
            )

            if (
                "max_mismatches" in text
                or "maximum_sites" in text
            ):

                python_candidates.append(path)

        except OSError:
            continue

    # Include esplicitamente lo script STEP 16.
    if (
        STEP16_SCRIPT.exists()
        and STEP16_SCRIPT not in python_candidates
    ):

        python_candidates.append(STEP16_SCRIPT)

    for path in python_candidates:

        audit = inspect_python_file(path)

        report["audit_codice"].append(audit)

        print("[FILE]", safe_relative(path))

        if audit["parametri_analyze"]:

            print(
                "  analyze():",
                audit["parametri_analyze"],
            )

    print(
        "[INFO] File Python controllati:",
        len(python_candidates),
    )

    # ============================================================
    # CONTROLLO 8 - INTERPRETAZIONE
    # ============================================================

    print_section(
        8,
        "INTERPRETAZIONE DIAGNOSTICA",
    )

    recommendations = []

    if not srna_records:

        recommendations.append(
            "Individuare il FASTA effettivamente "
            "utilizzato nello STEP 16."
        )

    if not transcript_sample:

        recommendations.append(
            "Verificare il percorso del trascrittoma "
            "realmente montato nel container Docker."
        )

    if srna_records and transcript_sample:

        hits_2 = sum(
            int(item["hits_reverse"].get(2, 0))
            for item in report["matching_reale"]
        )

        hits_5 = sum(
            int(item["hits_reverse"].get(5, 0))
            for item in report["matching_reale"]
        )

        if hits_2 == 0 and hits_5 > 0:

            recommendations.append(
                "Il campione contiene complementarita "
                "reverse con 3-5 mismatch ma non "
                "con massimo 2 mismatch. "
                "Verificare la severita della soglia."
            )

        elif hits_2 > 0:

            recommendations.append(
                "Il controllo indipendente trova "
                "corrispondenze reverse entro 2 mismatch. "
                "Confrontare orientamento e filtri "
                "del motore V3.3."
            )

        else:

            recommendations.append(
                "Il campione non mostra hit reverse "
                "entro 5 mismatch. Estendere l'audit "
                "prima di modificare il motore."
            )

    recommendations.append(
        "Eseguire successivamente un controllo "
        "positivo direttamente con analyze() "
        "nel container Docker."
    )

    recommendations.append(
        "Non considerare un hit Hamming "
        "come target biologico validato."
    )

    report["raccomandazioni"] = recommendations

    for recommendation in recommendations:

        print("-", recommendation)

    # ============================================================
    # CONTROLLO 9 - OUTPUT
    # ============================================================

    print_section(
        9,
        "GENERAZIONE REPORT",
    )

    report["stato"] = (
        "DIAGNOSI_PARZIALE"
        if (
            not srna_records
            or not transcript_sample
        )
        else "DIAGNOSI_COMPLETATA"
    )

    json_path = (
        output_dir
        / "RAPPORTO_STEP17_ZERO_TARGET.json"
    )

    csv_path = (
        output_dir
        / "MATCHING_SENSIBILITA_STEP17.csv"
    )

    markdown_path = (
        output_dir
        / "RAPPORTO_STEP17_ZERO_TARGET.md"
    )

    json_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    if csv_rows:

        with csv_path.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as handle:

            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "srna",
                    "lunghezza_nt",
                    "trascritti_testati",
                    "max_mismatch",
                    "hit_diretti",
                    "hit_reverse",
                ],
            )

            writer.writeheader()
            writer.writerows(csv_rows)

    md_lines = [
        "# AI PRIMETTA V3.3",
        "",
        "## STEP 17 - Diagnosi ZERO TARGET",
        "",
        f"**Stato:** {report['stato']}",
        "",
        "### Dataset",
        "",
        f"- STEP 16: {step16}",
        f"- FASTA small RNA: {srna_path}",
        (
            "- Trascritti B73: "
            + str(
                report["trascrittoma"].get(
                    "numero_trascritti",
                    "NON DISPONIBILE",
                )
            )
        ),
        (
            "- Mapping: "
            + str(
                report["mapping"].get(
                    "file",
                    "NON DISPONIBILE",
                )
            )
        ),
        "",
        "### Controlli sintetici",
        "",
        (
            "- Superati: "
            + str(synthetic["test_superato"])
        ),
        "",
        "### Matching su dati reali",
        "",
        (
            "Analisi indipendente senza gap, "
            "eseguita su un campione limitato."
        ),
        "",
    ]

    for item in report["matching_reale"]:

        md_lines.append(
            f"#### {item['srna_id']}"
        )

        md_lines.append("")

        md_lines.append(
            "| Mismatch | Diretti | Reverse |"
        )

        md_lines.append(
            "|---:|---:|---:|"
        )

        for threshold in MISMATCH_THRESHOLDS:

            md_lines.append(
                f"| {threshold} | "
                f"{item['hits_diretti'].get(threshold, 0)} | "
                f"{item['hits_reverse'].get(threshold, 0)} |"
            )

        md_lines.append("")

    md_lines.extend(
        [
            "### Raccomandazioni",
            "",
        ]
    )

    md_lines.extend(
        f"- {item}"
        for item in recommendations
    )

    md_lines.extend(
        [
            "",
            "### Avvertenze",
            "",
        ]
    )

    md_lines.extend(
        f"- {item}"
        for item in report["avvertenze"]
    )

    md_lines.extend(
        [
            "",
            "### Limiti scientifici",
            "",
            "- Non sono state validate interazioni RNA-gene.",
            "- Non e' stato eseguito il motore Docker V3.3.",
            "- Non sono state usate matrici WW/DS.",
            "- Nessun modello GNN e' stato applicato.",
            "- I controlli sintetici sono esclusivamente software.",
            "",
        ]
    )

    markdown_path.write_text(
        "\n".join(md_lines),
        encoding="utf-8",
    )

    print("[OK] JSON:", json_path)

    if csv_rows:
        print("[OK] CSV:", csv_path)

    print("[OK] Markdown:", markdown_path)

    # ============================================================
    # RIEPILOGO FINALE
    # ============================================================

    print()
    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17 TERMINATO")
    print("=" * 72)

    print("Stato:", report["stato"])

    print(
        "Small RNA identificati:",
        len(srna_records),
    )

    print(
        "Trascritti B73:",
        report["trascrittoma"].get(
            "numero_trascritti",
            "NON DISPONIBILE",
        ),
    )

    print(
        "Trascritti campione:",
        len(transcript_sample),
    )

    print(
        "Controlli sintetici:",
        "SUPERATI"
        if synthetic["test_superato"]
        else "NON SUPERATI",
    )

    print(
        "Avvertenze:",
        len(report["avvertenze"]),
    )

    print()
    print("[OK] Nessun file originale sovrascritto")
    print("[OK] Nessun training eseguito")
    print("[OK] Nessuna modifica programmata alla GUI")
    print("[OK] Nessuna modifica programmata ai modelli")

    print()
    print("CARTELLA RISULTATI:")
    print(output_dir)

    print()
    print(
        "[NOTA] STEP 17 = diagnosi software, "
        "non validazione biologica."
    )


if __name__ == "__main__":
    main()