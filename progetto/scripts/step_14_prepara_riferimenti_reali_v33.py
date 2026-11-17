
#!/usr/bin/env python3

"""
AI PRIMETTA V3.3
STEP 14 - PREPARAZIONE RIFERIMENTI BIOLOGICI REALI

Ambiente:
Visual Studio Code - Ubuntu 24.04

Obiettivo:
Verificare e collegare i riferimenti Ensembl Plants
di Zea mays alla pipeline RNA-gene target.

Nessuna modifica alla GUI originale.
Nessun download.
Nessuna duplicazione del trascrittoma.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import sys

from pathlib import Path


# ============================================================
# 1. PERCORSI DEL PROGETTO
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

REF_SOURCE = (
    ROOT /
    "data/reference/zea_mays_ensembl_plants"
)

REF_DEST = (
    ROOT /
    "target_prediction_v33_vscode/references"
)

STEM = "Zea_mays.Zm-B73-REFERENCE-NAM-5.0"

FASTA_NAME = f"{STEM}.cdna.all.fa.gz"

MAPPING_NAME = f"{STEM}.tx2gene.tsv"

INPUTS = {

    FASTA_NAME:
        REF_SOURCE / FASTA_NAME,

    MAPPING_NAME:
        REF_SOURCE / MAPPING_NAME
}

GUI_FILES = [

    "app/streamlit_app.py",

    "app/pages/01_Rete_neurale_GAE.py",

    "app/pages/02_Dati_ricercatore.py",

    "app/pages/03_Ranking_dati_ricercatore.py",

    "app/pages/04_Chat_dati_ricercatore.py"
]


# ============================================================
# 2. HASH SHA256
# ============================================================

def sha256(path: Path) -> str:

    h = hashlib.sha256()

    with path.open("rb") as stream:

        for block in iter(
            lambda: stream.read(1048576),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# ============================================================
# 3. VERIFICA MAPPING TRASCRITTO-GENE
# ============================================================

def read_tx2gene(path: Path):

    transcripts = {}

    genes = set()

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        reader = csv.DictReader(
            stream,
            delimiter="\t"
        )

        fields = set(
            reader.fieldnames or []
        )

        required = {
            "transcript_id",
            "gene_id"
        }

        if not required.issubset(fields):

            raise ValueError(
                "Mapping tx2gene non compatibile. "
                "Sono richieste le colonne "
                "transcript_id e gene_id. "
                f"Colonne trovate: {reader.fieldnames}"
            )

        for number, row in enumerate(
            reader,
            start=2
        ):

            tx = row[
                "transcript_id"
            ].strip()

            gene = row[
                "gene_id"
            ].strip()

            if not tx or not gene:

                raise ValueError(
                    f"Mapping incompleto alla riga {number}"
                )

            if tx in transcripts:

                if transcripts[tx] != gene:

                    raise ValueError(
                        f"Mapping conflittuale per {tx}"
                    )

            transcripts[tx] = gene

            genes.add(gene)

    if not transcripts:

        raise ValueError(
            "Mapping privo di identificativi."
        )

    return (
        set(transcripts),
        genes,
        len(transcripts)
    )


# ============================================================
# 4. VERIFICA IDENTIFICATIVI FASTA
# ============================================================

def sample_fasta_headers(
    path: Path,
    limit: int = 2500
):

    headers = []

    with gzip.open(
        path,
        "rt",
        encoding="utf-8"
    ) as stream:

        for line in stream:

            if line.startswith(">"):

                tokens = line[1:].split()

                identifier = (
                    tokens[0]
                    if tokens
                    else ""
                )

                if not identifier:

                    raise ValueError(
                        "Header FASTA vuoto."
                    )

                headers.append(
                    identifier
                )

                if len(headers) >= limit:
                    break

            elif not headers and line.strip():

                raise ValueError(
                    "FASTA non valido: manca >"
                )

    if not headers:

        raise ValueError(
            "Nessun trascritto trovato."
        )

    return headers


# ============================================================
# 5. VERIFICA PROGETTO E GUI
# ============================================================

def verify_project():

    for filename in GUI_FILES:

        path = ROOT / filename

        if not path.is_file():

            raise FileNotFoundError(
                f"File GUI mancante: {filename}"
            )

    addon = (
        ROOT /
        "target_prediction_v33_vscode"
    )

    required = [

        addon / "app_v33.py",

        addon / "compose.yaml"
    ]

    for path in required:

        if not path.is_file():

            raise FileNotFoundError(
                f"Componente V3.3 mancante: {path}"
            )

    if not REF_DEST.is_dir():

        raise FileNotFoundError(
            f"Cartella references mancante: {REF_DEST}"
        )


# ============================================================
# 6. PROCEDURA PRINCIPALE
# ============================================================

def main():

    print("=" * 65)

    print(
        "AI PRIMETTA V3.3 - STEP 14"
    )

    print(
        "PREPARAZIONE RIFERIMENTI REALI ZEA MAYS"
    )

    print("=" * 65)

    print()
    print("Progetto:", ROOT)

    verify_project()

    gui_before = {

        name: sha256(ROOT / name)

        for name in GUI_FILES
    }

    # --------------------------------------------------------
    # Verifica esistenza dei dati biologici
    # --------------------------------------------------------

    print()
    print("1. CONTROLLO FILE DI RIFERIMENTO")

    for name, path in INPUTS.items():

        if not path.is_file():

            raise FileNotFoundError(
                f"Riferimento mancante: {path}"
            )

        size_mib = (
            path.stat().st_size / 1048576
        )

        print(
            f"[OK] {name}: {size_mib:.1f} MiB"
        )

    # --------------------------------------------------------
    # Verifica mapping
    # --------------------------------------------------------

    print()
    print("2. CONTROLLO IDENTIFICATIVI BIOLOGICI")

    tx_ids, genes, n_tx = read_tx2gene(
        INPUTS[MAPPING_NAME]
    )

    sampled = sample_fasta_headers(
        INPUTS[FASTA_NAME]
    )

    hits = sum(
        tx in tx_ids
        for tx in sampled
    )

    coverage = hits / len(sampled)

    print(
        f"Trascritti nel mapping: {n_tx:,}"
    )

    print(
        f"Geni distinti: {len(genes):,}"
    )

    print(
        f"Trascritti FASTA campionati: {len(sampled):,}"
    )

    print(
        f"Identificativi corrispondenti: {hits:,}"
    )

    print(
        f"Compatibilita campionata: {coverage:.2%}"
    )

    if coverage < 0.95:

        print()
        print(
            "[ERRORE] Compatibilita FASTA/tx2gene "
            "inferiore al 95%."
        )

        print(
            "Nessun riferimento verra collegato."
        )

        return 2

    # --------------------------------------------------------
    # Protezione da sovrascritture
    # --------------------------------------------------------

    print()
    print("3. CONTROLLO DESTINAZIONI")

    for name, source in INPUTS.items():

        destination = REF_DEST / name

        if destination.exists():

            if not os.path.samefile(
                source,
                destination
            ):

                print(
                    "[ERRORE] File gia presente "
                    "ma differente:"
                )

                print(destination)

                return 3

    # --------------------------------------------------------
    # Collegamento dati - hard link
    # --------------------------------------------------------

    print()
    print("4. COLLEGAMENTO RIFERIMENTI REALI")

    newly_linked = []

    try:

        for name, source in INPUTS.items():

            destination = REF_DEST / name

            if destination.exists():

                print(
                    "[OK] Riferimento gia collegato:",
                    name
                )

            else:

                os.link(
                    source,
                    destination
                )

                newly_linked.append(
                    destination
                )

                print(
                    "[OK] Collegamento senza duplicazione:",
                    name
                )

        # ----------------------------------------------------
        # Verifica finale GUI originale
        # ----------------------------------------------------

        print()
        print("5. VERIFICA GUI ORIGINALE")

        gui_after = {

            name: sha256(ROOT / name)

            for name in GUI_FILES
        }

        if gui_before != gui_after:

            raise RuntimeError(
                "Sono cambiati file della GUI originale."
            )

    except Exception:

        for path in newly_linked:

            path.unlink(
                missing_ok=True
            )

        raise

    # --------------------------------------------------------
    # Risultati
    # --------------------------------------------------------

    print()
    print("=" * 65)

    print("[OK] RIFERIMENTI BIOLOGICI COLLEGATI")

    print("[OK] GUI V3.2 INVARIATA")

    print("[OK] CONTROLLI SHA256 SUPERATI")

    print("=" * 65)

    print()
    print("Directory riferimenti:")

    print(REF_DEST)

    print()
    print("Nella GUI V3.3 selezionare:")

    print()
    print("FASTA trascrittoma:")

    print(FASTA_NAME)

    print()
    print("Mapping tx2gene:")

    print(MAPPING_NAME)

    print()
    print("Indirizzo GUI:")

    print("http://localhost:8502")

    print()
    print(
        "NOTA: il controllo di compatibilita "
        "ha utilizzato un campione di intestazioni."
    )

    print(
        "Prima di elaborare tutti i 500 RNA "
        "servono test prestazionali e scientifici."
    )

    return 0


# ============================================================
# 7. AVVIO
# ============================================================

if __name__ == "__main__":

    try:

        sys.exit(
            main()
        )

    except (
        OSError,
        ValueError,
        RuntimeError
    ) as error:

        print(
            "[ERRORE]",
            error,
            file=sys.stderr
        )

        sys.exit(1)
