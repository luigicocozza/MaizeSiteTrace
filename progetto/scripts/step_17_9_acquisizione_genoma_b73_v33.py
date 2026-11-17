#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===========================================================
AI PRIMETTA V3.3
STEP 17.9 - ACQUISIZIONE GENOMA REALE B73 NAM-5.0
===========================================================

Sorgente:
https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/

Obiettivi:
1. Verificare spazio disco disponibile.
2. Acquisire il FASTA genomico reale.
3. Riutilizzare eventuale download valido.
4. Verificare integrita gzip.
5. Verificare cromosomi e sequenze.
6. Calcolare SHA256.
7. Generare rapporto JSON.

Nessuna modifica a:
- GUI V3.2
- Motore V3.3
- Modelli ML/GNN
- Dati sperimentali originali

NON esegue training.
NON crea target biologici artificiali.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


# ===========================================================
# CONFIGURAZIONE
# ===========================================================

ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"

REFERENCES = BASE / "references"
RUNS = BASE / "runs"

ASSEMBLY = "Zm-B73-REFERENCE-NAM-5.0"

FILENAME = ASSEMBLY + ".fa.gz"

URL = (
    "https://download.maizegdb.org/"
    "Zm-B73-REFERENCE-NAM-5.0/"
    + FILENAME
)

DESTINATION = REFERENCES / FILENAME

PARTIAL = REFERENCES / (FILENAME + ".part")

# Spazio di sicurezza per il FASTA compresso e i report.
MIN_FREE_GB = 1.5

# Dimensione minima conservativa: serve soltanto
# a individuare file chiaramente incompleti.
MIN_COMPRESSED_MB = 300

EXPECTED_CHROMOSOMES = {
    str(i) for i in range(1, 11)
}


# ===========================================================
# UTILITA
# ===========================================================

def section(number, title):
    print()
    print("=" * 70)
    print(f"[{number}] {title}")
    print("=" * 70)


def utc_now():
    return datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")


def sha256_file(path):

    digest = hashlib.sha256()

    with path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(4 * 1024 * 1024),
            b""
        ):
            digest.update(block)

    return digest.hexdigest()


def free_space_gb(path):

    usage = shutil.disk_usage(path)

    return usage.free / (1024 ** 3)


def validate_genome(path):
    """
    Verifica integralmente il flusso gzip.

    Conta le sequenze e i nucleotidi senza
    caricare l'intero genoma in RAM.

    Non sostituisce la verifica di un checksum
    ufficiale di provenienza.
    """

    if not path.is_file():
        raise FileNotFoundError(path)

    size_mb = path.stat().st_size / (1024 ** 2)

    if size_mb < MIN_COMPRESSED_MB:
        raise RuntimeError(
            "FASTA troppo piccolo: "
            f"{size_mb:.2f} MB"
        )

    print("[INFO] Verifica gzip e sequenze...")

    records = 0
    total_bases = 0

    current_header = None
    current_length = 0

    headers = []
    chromosomes = set()

    lengths = []

    invalid_bases = Counter()

    def finish_record():

        nonlocal records, current_length

        if current_header is None:
            return

        records += 1

        lengths.append(current_length)

        identifier = current_header.split()[0]

        if len(headers) < 30:
            headers.append(identifier)

        # Riconosce identificativi cromosomici
        # 1..10, chr1..chr10 e chromosome_1..10.
        match = re.fullmatch(
            r"(?:chr|chromosome_?)?([1-9]|10)",
            identifier,
            flags=re.IGNORECASE
        )

        if match:
            chromosomes.add(match.group(1))

        current_length = 0

    with gzip.open(path, "rt", encoding="utf-8") as f:

        for raw in f:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                finish_record()

                current_header = line[1:].strip()

            else:

                if current_header is None:
                    raise RuntimeError(
                        "FASTA non valido: "
                        "sequenza prima dell'intestazione"
                    )

                sequence = line.upper()

                current_length += len(sequence)
                total_bases += len(sequence)

                invalid_bases.update(
                    base
                    for base in sequence
                    if base not in "ACGTN"
                )

    finish_record()

    missing = sorted(
        EXPECTED_CHROMOSOMES - chromosomes,
        key=int
    )

    if records == 0:
        raise RuntimeError(
            "Nessuna sequenza FASTA trovata"
        )

    if total_bases < 1_000_000_000:
        raise RuntimeError(
            "Quantita di sequenza insufficiente "
            "per il genoma B73 atteso"
        )

    result = {
        "file": str(path),
        "dimensione_mb": round(size_mb, 3),
        "numero_record": records,
        "nucleotidi_totali": total_bases,
        "lunghezza_massima": max(lengths),
        "cromosomi_identificati": sorted(
            chromosomes, key=int
        ),
        "cromosomi_mancanti": missing,
        "prime_intestazioni": headers,
        "basi_non_acgtn": dict(invalid_bases),
        "gzip_letto_integralmente": True,
        "sha256": sha256_file(path),
    }

    return result


# ===========================================================
# DOWNLOAD
# ===========================================================

def download_genome():

    if not shutil.which("curl"):
        raise RuntimeError(
            "curl non installato. "
            "Verificare la disponibilita del comando."
        )

    print("[INFO] Sorgente:", URL)
    print("[INFO] Destinazione:", DESTINATION)

    command = [
        "curl",
        "--fail",
        "--location",
        "--show-error",
        "--progress-bar",
        "--retry", "3",
        "--retry-delay", "5",
        "--connect-timeout", "30",
        "--continue-at", "-",
        "--output", str(PARTIAL),
        URL,
    ]

    print(
        "[INFO] Download avviato. "
        "Il file parziale resta separato "
        "dal FASTA definitivo.",
        flush=True
    )

    subprocess.run(
        command,
        check=True
    )

    print("[OK] Trasferimento terminato")


# ===========================================================
# MAIN
# ===========================================================

def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.9")
    print("ACQUISIZIONE GENOMA REALE B73 NAM-5.0")
    print("=" * 70)

    if not ROOT.is_dir():
        raise SystemExit(
            f"[ERRORE] Progetto non trovato: {ROOT}"
        )

    REFERENCES.mkdir(
        parents=True,
        exist_ok=True
    )

    RUNS.mkdir(
        parents=True,
        exist_ok=True
    )

    # -------------------------------------------------------
    # 1. CONTROLLO SPAZIO
    # -------------------------------------------------------

    section(1, "CONTROLLO SPAZIO DISCO")

    free_gb = free_space_gb(REFERENCES)

    print(
        f"[INFO] Spazio disponibile: {free_gb:.2f} GB"
    )

    if not DESTINATION.exists():

        if free_gb < MIN_FREE_GB:
            raise SystemExit(
                "[ERRORE] Spazio disco insufficiente "
                "per il download in sicurezza."
            )

    # -------------------------------------------------------
    # 2. ACQUISIZIONE
    # -------------------------------------------------------

    section(2, "CONTROLLO FASTA GENOMICO")

    genome_info = None

    source_status = None

    if DESTINATION.is_file():

        print("[INFO] Genoma gia presente")
        print(DESTINATION)

        try:

            genome_info = validate_genome(
                DESTINATION
            )

            source_status = (
                "FILE_ESISTENTE_VERIFICATO"
            )

            print(
                "[OK] FASTA esistente valido"
            )

        except Exception as error:

            # Non eliminare o sovrascrivere il file.
            raise SystemExit(
                "[ERRORE] FASTA esistente non valido: "
                + str(error)
                + "\nFile preservato. "
                "Richiede verifica manuale."
            )

    else:

        if PARTIAL.exists():

            print(
                "[INFO] Download parziale individuato:"
            )

            print(PARTIAL)

        download_genome()

        # ---------------------------------------------------
        # 3. VALIDAZIONE DOWNLOAD
        # ---------------------------------------------------

        section(3, "VALIDAZIONE FASTA SCARICATO")

        genome_info = validate_genome(
            PARTIAL
        )

        print("[OK] Integrita gzip verificata")

        if genome_info["cromosomi_mancanti"]:

            raise SystemExit(
                "[ERRORE] Non identificati i cromosomi "
                + ", ".join(
                    genome_info["cromosomi_mancanti"]
                )
                + ". File parziale preservato."
            )

        # Rename atomico sullo stesso filesystem.
        # Non sovrascrive file precedentemente validi.
        if DESTINATION.exists():
            raise SystemExit(
                "[ERRORE] File destinazione gia esistente"
            )

        PARTIAL.rename(DESTINATION)

        genome_info["file"] = str(
            DESTINATION
        )

        source_status = (
            "DOWNLOAD_COMPLETATO_VALIDATO"
        )

    # -------------------------------------------------------
    # 4. RISULTATI
    # -------------------------------------------------------

    section(4, "RIEPILOGO GENOMA")

    print(
        "Dimensione compressa:",
        genome_info["dimensione_mb"],
        "MB"
    )

    print(
        "Sequenze FASTA:",
        genome_info["numero_record"]
    )

    print(
        "Nucleotidi:",
        genome_info["nucleotidi_totali"]
    )

    print(
        "Cromosomi:",
        genome_info["cromosomi_identificati"]
    )

    print(
        "Cromosomi mancanti:",
        genome_info["cromosomi_mancanti"]
    )

    print(
        "SHA256:",
        genome_info["sha256"]
    )

    print("\nPrime intestazioni:")

    for header in genome_info["prime_intestazioni"][:12]:
        print(" >", header)

    # -------------------------------------------------------
    # 5. REPORT
    # -------------------------------------------------------

    section(5, "SALVATAGGIO REPORT")

    timestamp = utc_now()

    output = (
        RUNS
        / f"STEP17_9_GENOMA_B73_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.9",
        "stato": (
            "GENOMA_ACQUISITO_E_VERIFICATO"
            if not genome_info["cromosomi_mancanti"]
            else "GENOMA_DA_VERIFICARE"
        ),
        "assemblaggio_atteso": ASSEMBLY,
        "url_sorgente": URL,
        "stato_file": source_status,
        "genoma": genome_info,
        "checksum_ufficiale_confrontato": False,
        "gtf_originale_preservato": True,
        "gui_modificata": False,
        "motore_modificato": False,
        "training_eseguito": False,
        "validazione_biologica": False,
        "prossimo_step": (
            "17.10 - Mapping genomico esplorativo "
            "dei cinque sRNA reali"
        ),
        "nota": (
            "Controllo strutturale del genoma. "
            "Il confronto con checksum ufficiale "
            "deve essere effettuato separatamente "
            "prima della validazione definitiva "
            "della provenienza."
        ),
    }

    json_path = (
        output / "RAPPORTO_STEP17_9.json"
    )

    json_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print("[OK] Report:", json_path)

    print()
    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.9 TERMINATO")
    print("=" * 70)

    print("Stato:", report["stato"])

    print("Genoma:", DESTINATION)

    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun training")
    print("[OK] Dataset PRJNA419326 preservato")

    print(
        "\n[NOTA] Prossima verifica: "
        "checksum ufficiale e coordinate "
        "prima del mapping genomico."
    )


if __name__ == "__main__":
    main()