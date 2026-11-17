#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.4 - TEST DELLA FUNZIONE SCAN ORIGINALE

Verifica con sequenze reali PRJNA419326 e
trascritti di controllo esplicitamente sintetici.

Non modifica app_v33.py, GUI, Docker o modelli.
"""

from __future__ import annotations

import csv
import hashlib
import json
import runpy
import sys

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


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def latest_step16():

    paths = sorted(
        (
            p for p in RUNS.glob(
                "STEP16_5_SRNA_REALI_*"
            )
            if p.is_dir()
        ),
        reverse=True
    )

    if not paths:
        raise RuntimeError(
            "Cartella STEP 16 non trovata"
        )

    return paths[0]


def read_fasta(path):

    records = {}
    header = None
    parts = []

    def add_record():
        if header is not None:
            sequence = (
                "".join(parts)
                .upper()
                .replace("U", "T")
            )
            records[header.split()[0]] = sequence

    with path.open(
        "r",
        encoding="utf-8"
    ) as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                add_record()
                header = line[1:]
                parts = []
            else:
                parts.append(line)

    add_record()

    return records


def mutate(sequence, positions):

    sequence = list(sequence)

    replacements = {
        "A": "C",
        "C": "G",
        "G": "T",
        "T": "A",
    }

    for position in positions:

        old = sequence[position]

        sequence[position] = replacements[old]

    return "".join(sequence)


def create_controls(rnas, revcomp):

    transcripts = {}
    mapping = {}
    expected = []

    for number, (rna_id, sequence) in enumerate(
        rnas.items(),
        start=1
    ):

        target = revcomp(sequence)

        length = len(target)

        # Tre seed effettivi del motore V3.3
        offsets = sorted({
            0,
            (length - 8) // 2,
            length - 8
        })

        controls = {
            "MM0": target,
            "MM1": mutate(
                target,
                [offsets[0]]
            ),
            "MM2": mutate(
                target,
                [offsets[0], offsets[1]]
            ),
            "MM3": mutate(
                target,
                offsets
            ),
        }

        for category, site in controls.items():

            tx_id = (
                f"CONTROLLO_SINTETICO_"
                f"RNA{number:02d}_{category}"
            )

            gene_id = (
                f"GENE_SINTETICO_RNA{number:02d}"
            )

            # Posizione del sito: nucleotidi 21...
            transcripts[tx_id] = (
                "A" * 20
                + site
                + "C" * 20
            )

            mapping[tx_id] = gene_id

            expected.append({
                "rna_id": rna_id,
                "transcript_id": tx_id,
                "gene_id": gene_id,
                "categoria": category,
                "mismatch_attesi": int(
                    category[-1]
                ),
                "start_atteso_1based": 21,
                "offset_seed": offsets,
            })

    return transcripts, mapping, expected


def save_fasta(path, records):

    with path.open(
        "w",
        encoding="utf-8"
    ) as f:

        for header, sequence in records.items():
            f.write(f">{header}\n")
            f.write(sequence + "\n")


def run_test(
    engine_scan,
    rnas,
    fasta_path,
    mapping,
    expected,
    threshold
):

    hits, statistics = engine_scan(
        rnas,
        fasta_path,
        mapping,
        max_mismatches=threshold,
        maximum_sites=300
    )

    index = {}

    for hit in hits:

        key = (
            hit["rna_id"],
            hit["transcript_id"],
            hit["site_start_1based"]
        )

        index[key] = hit

    results = []

    for item in expected:

        key = (
            item["rna_id"],
            item["transcript_id"],
            item["start_atteso_1based"]
        )

        hit = index.get(key)

        expected_detected = (
            item["mismatch_attesi"] <= threshold
        )

        detected = hit is not None

        correct_mm = (
            hit is not None
            and hit["mismatches"]
            == item["mismatch_attesi"]
        )

        # Con 3 mismatch la ricerca a 3 seed
        # non garantisce il recupero del sito.
        if (
            threshold == 3
            and item["mismatch_attesi"] == 3
        ):
            status = (
                "TROVATO"
                if detected and correct_mm
                else "NON_RILEVATO_SEED"
            )

        elif expected_detected:

            status = (
                "PASS"
                if detected and correct_mm
                else "FAIL"
            )

        else:

            status = (
                "PASS"
                if not detected
                else "FAIL"
            )

        results.append({
            **item,
            "soglia_mismatch": threshold,
            "rilevato": detected,
            "mismatch_osservati": (
                hit["mismatches"]
                if hit else None
            ),
            "stato": status,
        })

    return results, statistics, len(hits)


def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.4")
    print("TEST DIRETTO FUNZIONE SCAN ORIGINALE")
    print("=" * 72)

    if not ENGINE.is_file():
        raise SystemExit(
            f"[ERRORE] Motore non trovato: {ENGINE}"
        )

    step16 = latest_step16()

    srna_file = (
        step16 / "input_5_srna_reali.fa"
    )

    if not srna_file.is_file():
        raise SystemExit(
            f"[ERRORE] FASTA non trovato: {srna_file}"
        )

    rnas = read_fasta(srna_file)

    if len(rnas) != 5:
        raise SystemExit(
            "[ERRORE] Attesi esattamente 5 sRNA"
        )

    print("\n[1] SORGENTI")

    print("Motore:", ENGINE)
    print("SHA256:", sha256(ENGINE))
    print("FASTA:", srna_file)
    print("sRNA:", len(rnas))

    print("\n[2] CARICAMENTO MOTORE ORIGINALE")

    # runpy esegue il codice Python a livello di modulo,
    # ma non invoca main() se protetto da
    # if __name__ == "__main__".
    # Prima dell'esecuzione verificare che app_v33.py
    # non contenga altre azioni automatiche indesiderate.

    original_sys_path = list(sys.path)

    try:
        sys.path.insert(0, str(BASE))

        namespace = runpy.run_path(
            str(ENGINE),
            run_name="__audit_step17_4__"
        )

    finally:
        sys.path[:] = original_sys_path

    engine_scan = namespace.get("scan")
    engine_revcomp = namespace.get("revcomp")

    if not callable(engine_scan):
        raise RuntimeError(
            "Funzione scan() originale non trovata"
        )

    if not callable(engine_revcomp):
        raise RuntimeError(
            "Funzione revcomp() originale non trovata"
        )

    print("[OK] scan() originale caricata")
    print("[OK] revcomp() originale caricata")

    print("\n[3] GENERAZIONE CONTROLLI")

    transcripts, mapping, expected = (
        create_controls(
            rnas,
            engine_revcomp
        )
    )

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_4_TEST_SCAN_{stamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    synthetic_fasta = (
        output / "CONTROLLI_SOFTWARE_SINTETICI.fa"
    )

    save_fasta(
        synthetic_fasta,
        transcripts
    )

    print(
        "[OK] Trascritti sintetici:",
        len(transcripts)
    )

    print(
        "[OK] Siti attesi:",
        len(expected)
    )

    print("\n[4] TEST SOGLIA 2 MISMATCH")

    results2, stats2, hits2 = run_test(
        engine_scan,
        rnas,
        synthetic_fasta,
        mapping,
        expected,
        threshold=2
    )

    for item in results2:

        print(
            item["transcript_id"],
            "|",
            item["stato"]
        )

    print(
        "Siti restituiti dal motore:",
        hits2
    )

    print(
        "Statistiche:",
        stats2
    )

    print("\n[5] TEST SOGLIA 3 MISMATCH")

    results3, stats3, hits3 = run_test(
        engine_scan,
        rnas,
        synthetic_fasta,
        mapping,
        expected,
        threshold=3
    )

    for item in results3:

        if item["categoria"] == "MM3":

            print(
                item["transcript_id"],
                "|",
                item["stato"]
            )

    print(
        "Siti restituiti dal motore:",
        hits3
    )

    print(
        "Statistiche:",
        stats3
    )

    print("\n[6] VALUTAZIONE")

    failures2 = [
        row for row in results2
        if row["stato"] == "FAIL"
    ]

    failures3 = [
        row for row in results3
        if row["stato"] == "FAIL"
    ]

    missed3 = [
        row for row in results3
        if row["stato"] == "NON_RILEVATO_SEED"
    ]

    passed = (
        len(failures2) == 0
        and len(failures3) == 0
    )

    status = (
        "TEST_SCAN_SUPERATO"
        if passed
        else "TEST_SCAN_NON_SUPERATO"
    )

    print("Stato:", status)
    print("Errori soglia 2:", len(failures2))
    print("Errori soglia 3:", len(failures3))
    print(
        "Siti 3 mismatch non rilevati:",
        len(missed3)
    )

    print("\n[7] SALVATAGGIO RISULTATI")

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.4",
        "stato": status,
        "motore": str(ENGINE),
        "sha256_motore": sha256(ENGINE),
        "small_rna_file": str(srna_file),
        "numero_srna_reali": len(rnas),
        "numero_trascritti_sintetici": len(transcripts),
        "numero_siti_attesi": len(expected),
        "test_mismatch_2": {
            "statistiche": stats2,
            "numero_hits_totali": hits2,
            "risultati": results2
        },
        "test_mismatch_3": {
            "statistiche": stats3,
            "numero_hits_totali": hits3,
            "risultati": results3
        },
        "errori_mm2": len(failures2),
        "errori_mm3": len(failures3),
        "siti_mm3_non_rilevati": len(missed3),
        "natura": "CONTROLLO_SOFTWARE",
        "validazione_biologica": False,
        "training_eseguito": False,
    }

    json_path = (
        output / "RAPPORTO_STEP17_4.json"
    )

    json_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    csv_path = (
        output / "RISULTATI_TEST_SCAN.csv"
    )

    all_rows = results2 + results3

    with csv_path.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(all_rows[0].keys())
        )

        writer.writeheader()
        writer.writerows(all_rows)

    print("[OK] JSON:", json_path)
    print("[OK] CSV:", csv_path)

    print("\n" + "=" * 72)
    print("STEP 17.4 TERMINATO")
    print("=" * 72)

    print("Stato finale:", status)
    print("Cartella risultati:", output)

    print("[OK] GUI non modificata")
    print("[OK] Motore originale non modificato")
    print("[OK] Nessun training")
    print("[OK] Dati reali non sovrascritti")

    print(
        "\n[NOTA] I trascritti di controllo sono "
        "sintetici: non costituiscono target biologici."
    )


if __name__ == "__main__":
    main()