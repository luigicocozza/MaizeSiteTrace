#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.5 - SCANSIONE REALE CON SOGLIA 3

- Usa i 5 sRNA reali dello STEP 16
- Usa i 72539 trascritti B73
- Usa il mapping reale verificato nello STEP 17.1
- Esegue la funzione scan() originale
- Confronta candidati <=2 e =3 mismatch
- Non modifica GUI, motore, Docker o modelli
- Non effettua validazione biologica
"""

import csv
import hashlib
import json
import runpy
import sys

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"

ENGINE = BASE / "app_v33.py"

RUNS = BASE / "runs"

STEM = "Zea_mays.Zm-B73-REFERENCE-NAM-5.0"

TRANSCRIPTOME = (
    BASE / "references"
    / f"{STEM}.cdna.all.fa.gz"
)

MAPPING = (
    BASE / "references"
    / f"{STEM}.tx2gene.tsv"
)

MAX_MISMATCHES = 3

# Limite imposto alle righe restituite da scan()
# per ciascuno small RNA.
MAXIMUM_SITES = 3000


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b""
        ):
            digest.update(block)

    return digest.hexdigest()


def latest_step16():

    candidates = sorted(
        [
            p for p in RUNS.glob(
                "STEP16_5_SRNA_REALI_*"
            )
            if p.is_dir()
        ],
        reverse=True,
    )

    if not candidates:
        raise RuntimeError(
            "STEP 16 non trovato"
        )

    return candidates[0]


def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.5")
    print("SCANSIONE RNA-GENE SU TRASCRITTOMA REALE")
    print("=" * 72)

    step16 = latest_step16()

    RNA_FASTA = (
        step16 / "input_5_srna_reali.fa"
    )

    inputs = [
        ENGINE,
        RNA_FASTA,
        TRANSCRIPTOME,
        MAPPING,
    ]

    for path in inputs:

        if not path.is_file():
            raise FileNotFoundError(path)

    print("\n[1] INPUT REALI")

    print("Small RNA:", RNA_FASTA)
    print("Trascrittoma:", TRANSCRIPTOME)
    print("Mapping:", MAPPING)

    print("\n[2] CARICAMENTO MOTORE")

    # Esegue il modulo con nome diverso da __main__,
    # come gia' verificato nello STEP 17.4.
    # Non viene richiamata la funzione main().

    previous_path = list(sys.path)

    try:
        sys.path.insert(0, str(BASE))

        namespace = runpy.run_path(
            str(ENGINE),
            run_name="__audit_step17_5__",
        )

    finally:
        sys.path[:] = previous_path

    required = ["scan", "fasta", "table"]

    for name in required:

        if not callable(namespace.get(name)):
            raise RuntimeError(
                f"Funzione mancante: {name}"
            )

    scan = namespace["scan"]
    fasta = namespace["fasta"]
    table = namespace["table"]

    print("[OK] Motore originale caricato")

    print("\n[3] LETTURA DEI CINQUE SMALL RNA")

    rnas = dict(fasta(RNA_FASTA))

    if len(rnas) != 5:
        raise RuntimeError(
            f"Attesi 5 sRNA, trovati {len(rnas)}"
        )

    for rid, sequence in rnas.items():

        print(
            f"{rid[:65]} | "
            f"{len(sequence)} nt"
        )

    print("\n[4] LETTURA MAPPING REALE")

    mapping = {}

    for row in table(MAPPING):

        tx = row["transcript_id"].strip()
        gene = row["gene_id"].strip()

        if not tx or not gene:
            raise RuntimeError(
                "Mapping con ID vuoti"
            )

        if tx in mapping and mapping[tx] != gene:
            raise RuntimeError(
                f"Mapping ambiguo: {tx}"
            )

        mapping[tx] = gene

    print(
        "[OK] Trascritti nel mapping:",
        len(mapping),
    )

    print("\n[5] SCANSIONE REALE")

    print(
        "[INFO] Soglia massima:",
        MAX_MISMATCHES,
    )

    print(
        "[INFO] Limite siti per sRNA:",
        MAXIMUM_SITES,
    )

    print(
        "[INFO] Esecuzione di scan() "
        "sul trascrittoma completo",
        flush=True,
    )

    sites, statistics = scan(
        rnas,
        TRANSCRIPTOME,
        mapping,
        max_mismatches=MAX_MISMATCHES,
        maximum_sites=MAXIMUM_SITES,
    )

    print("[OK] Scansione completata")

    print("\n[6] DISTRIBUZIONE MISMATCH")

    distribution = Counter(
        int(site["mismatches"])
        for site in sites
    )

    for mismatches in range(4):

        print(
            f"Mismatch {mismatches}:",
            distribution[mismatches],
        )

    sites_0_2 = [
        s for s in sites
        if int(s["mismatches"]) <= 2
    ]

    sites_3 = [
        s for s in sites
        if int(s["mismatches"]) == 3
    ]

    print("\n[7] RISULTATI")

    print(
        "Trascritti scansionati:",
        statistics.get("transcripts_scanned"),
    )

    print(
        "Trascritti senza mapping:",
        statistics.get(
            "transcripts_without_tx2gene"
        ),
    )

    print(
        "Siti restituiti totali:",
        len(sites),
    )

    print(
        "Siti con 0-2 mismatch:",
        len(sites_0_2),
    )

    print(
        "Siti con 3 mismatch:",
        len(sites_3),
    )

    genes = set(
        s["gene_id"]
        for s in sites
    )

    print(
        "Geni distinti tra i siti restituiti:",
        len(genes),
    )

    counts_per_rna = Counter(
        s["rna_id"]
        for s in sites
    )

    truncated_possible = [
        rid
        for rid, count in counts_per_rna.items()
        if count >= MAXIMUM_SITES
    ]

    if truncated_possible:

        print(
            "[ATTENZIONE] Possibile limite "
            "di siti raggiunto per:",
            len(truncated_possible),
            "sRNA",
        )

    print("\n[8] SALVATAGGIO")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_5_SCAN_REALE_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False,
    )

    csv_file = output / "SITI_CANDIDATI_REALI.csv"

    fields = [
        "rna_id",
        "transcript_id",
        "gene_id",
        "site_start_1based",
        "site_end_1based",
        "mismatches",
        "identity_fraction",
        "validation_status",
    ]

    with csv_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(sites)

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.5",
        "stato": "SCANSIONE_TERMINATA",
        "motore_sha256": sha256(ENGINE),
        "rna_fasta_sha256": sha256(RNA_FASTA),
        "trascrittoma_sha256": sha256(TRANSCRIPTOME),
        "mapping_sha256": sha256(MAPPING),
        "soglia_mismatch": MAX_MISMATCHES,
        "maximum_sites_per_rna": MAXIMUM_SITES,
        "statistiche_motore": statistics,
        "siti_restituiti": len(sites),
        "siti_0_2_mismatch": len(sites_0_2),
        "siti_3_mismatch": len(sites_3),
        "geni_distinti": len(genes),
        "distribuzione_mismatch": {
            str(i): distribution[i]
            for i in range(4)
        },
        "conteggi_per_rna": dict(counts_per_rna),
        "possibile_troncamento_rna": truncated_possible,
        "limitazioni": [
            "Ricerca basata su tre seed da 8 nt",
            "Sensibilita incompleta a 3 mismatch",
            "Numero di siti restituiti limitato per sRNA",
            "Nessuna validazione biologica indipendente",
            "Espressione WW/DS non valutata",
            "GNN non applicata",
        ],
    }

    json_file = (
        output / "RAPPORTO_STEP17_5.json"
    )

    json_file.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print("[OK] CSV:", csv_file)
    print("[OK] JSON:", json_file)

    print("\n" + "=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.5 TERMINATO")
    print("=" * 72)

    print("Siti 0-2 mismatch:", len(sites_0_2))
    print("Siti 3 mismatch:", len(sites_3))
    print("Geni distinti:", len(genes))

    print("\nRisultati:", output)

    print("[OK] GUI non modificata")
    print("[OK] Motore non modificato")
    print("[OK] Modelli invariati")
    print("[OK] Nessun training")

    print(
        "\n[NOTA] I siti sono ipotesi "
        "computazionali non validate."
    )


if __name__ == "__main__":
    main()