#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.7 - AUDIT RISORSE GENOMICHE

Obiettivi:
- Leggere i risultati dello STEP 17.6
- Esaminare i migliori allineamenti esplorativi
- Individuare risorse FASTA genomiche, GFF3 e GTF
- Verificare la disponibilita di dati per una
  futura analisi delle origini genomiche degli sRNA

Nessuna modifica a GUI, modelli, motore o dataset.
Nessun download.
Nessun training.
"""

from __future__ import annotations

import csv
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"

EXCLUDED = {
    ".git",
    ".venv",
    "venv",
    "__pycache__",
    "node_modules",
}

FASTA_SUFFIXES = (
    ".fa",
    ".fasta",
    ".fna",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
)

ANNOTATION_SUFFIXES = (
    ".gff",
    ".gff3",
    ".gtf",
    ".gff.gz",
    ".gff3.gz",
    ".gtf.gz",
)


def latest_run(prefix):

    folders = sorted(
        [
            p for p in RUNS.glob(prefix + "*")
            if p.is_dir()
        ],
        reverse=True,
    )

    if not folders:
        raise RuntimeError(
            f"Nessun risultato trovato: {prefix}"
        )

    return folders[0]


def load_step17_6():

    folder = latest_run(
        "STEP17_6_DIAGNOSI_SEED_"
    )

    path = folder / "RAPPORTO_STEP17_6.json"

    if not path.is_file():
        raise FileNotFoundError(path)

    data = json.loads(
        path.read_text(encoding="utf-8")
    )

    return data, path


def find_resources():

    resources = []

    for directory, dirs, filenames in os.walk(ROOT):

        dirs[:] = [
            d for d in dirs
            if d not in EXCLUDED
            and not d.startswith("STEP17_")
        ]

        for filename in filenames:

            name = filename.lower()

            if name.endswith(
                FASTA_SUFFIXES
                + ANNOTATION_SUFFIXES
            ):

                path = Path(directory) / filename

                size = path.stat().st_size

                if name.endswith(
                    ANNOTATION_SUFFIXES
                ):
                    category = "ANNOTAZIONE"

                elif any(
                    term in name
                    for term in [
                        "genome",
                        "genomic",
                        "dna.primary",
                        "dna.toplevel",
                        "dna_sm",
                        "dna_rm",
                    ]
                ):
                    category = "GENOMA_POSSIBILE"

                elif any(
                    term in name
                    for term in [
                        "cdna",
                        "transcript",
                        "ncrna",
                        "cds",
                        "pep",
                    ]
                ):
                    category = "TRASCRITTOMA_O_SEQUENZE"

                else:
                    category = "FASTA_DA_CLASSIFICARE"

                resources.append({
                    "tipo": category,
                    "file": str(path),
                    "dimensione_mb": round(
                        size / (1024 ** 2), 3
                    ),
                })

    return sorted(
        resources,
        key=lambda x: (
            x["tipo"],
            x["file"],
        )
    )


def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.7")
    print("AUDIT RISORSE GENOMICHE E CANDIDATI REALI")
    print("=" * 72)

    data, source = load_step17_6()

    print("\n[1] RISULTATI STEP 17.6")
    print("Report:", source)

    results = data["rna_results"]

    print("Small RNA:", len(results))

    summary = []
    examples = []

    for rna_id, info in results.items():

        row = {
            "rna_id": rna_id,
            "lunghezza_nt": info["lunghezza"],
            "occorrenze_seed": info[
                "anchor_occurrences"
            ],
            "posizioni_valide": info[
                "candidate_valid_positions"
            ],
            "minimo_mismatch": info[
                "best_mismatches"
            ],
            "identita_migliore_percento": (
                round(
                    100 * (
                        1 - info["best_mismatches"]
                        / info["lunghezza"]
                    ),
                    3
                )
                if info["best_mismatches"] is not None
                else None
            ),
        }

        summary.append(row)

        print("\nRNA:", rna_id)
        print(
            "Minimo mismatch:",
            row["minimo_mismatch"]
        )
        print(
            "Identita migliore:",
            row["identita_migliore_percento"],
            "%"
        )

        for example in info.get("examples", []):

            examples.append({
                "rna_id": rna_id,
                **example,
            })

    print("\n[2] RISORSE DISPONIBILI")

    resources = find_resources()

    counts = Counter(
        r["tipo"] for r in resources
    )

    for category, number in sorted(counts.items()):
        print(f"{category}: {number}")

    print("\n[3] POSSIBILI GENOMI")

    genome_candidates = [
        resource
        for resource in resources
        if resource["tipo"] == "GENOMA_POSSIBILE"
    ]

    if genome_candidates:

        for resource in genome_candidates:

            print(
                resource["dimensione_mb"],
                "MB |",
                resource["file"]
            )

    else:

        print(
            "[INFO] Nessun FASTA chiaramente "
            "identificabile come genomico "
            "dal nome del file."
        )

    print("\n[4] ANNOTAZIONI")

    annotations = [
        resource
        for resource in resources
        if resource["tipo"] == "ANNOTAZIONE"
    ]

    for resource in annotations:

        print(
            resource["dimensione_mb"],
            "MB |",
            resource["file"]
        )

    if not annotations:
        print("[INFO] Nessuna annotazione individuata")

    print("\n[5] GENERAZIONE RAPPORTI")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_7_AUDIT_GENOMICO_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.7",
        "stato": "AUDIT_COMPLETATO",
        "fonte_step17_6": str(source),
        "risultati_srna": summary,
        "migliori_allineamenti": examples,
        "risorse": resources,
        "possibili_genomi": genome_candidates,
        "annotazioni": annotations,
        "genoma_verificato": False,
        "classificazione_srna_verificata": False,
        "validazione_biologica": False,
        "note": [
            "La classificazione delle risorse "
            "si basa sui nomi dei file.",
            "Le risorse genomiche devono essere "
            "verificate prima di utilizzarle.",
            "Gli allineamenti con 6-10 mismatch "
            "non sono target biologici validati.",
            "Le ipotesi di regolazione genomica "
            "richiedono analisi dedicate.",
        ],
    }

    json_path = output / "RAPPORTO_STEP17_7.json"

    json_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    def save_csv(path, rows):

        if not rows:
            return

        fields = list(rows[0].keys())

        with path.open(
            "w",
            newline="",
            encoding="utf-8"
        ) as handle:

            writer = csv.DictWriter(
                handle,
                fieldnames=fields,
                extrasaction="ignore"
            )

            writer.writeheader()
            writer.writerows(rows)

    save_csv(
        output / "SINTESI_SRNA.csv",
        summary
    )

    save_csv(
        output / "MIGLIORI_ALLINEAMENTI.csv",
        examples
    )

    save_csv(
        output / "INVENTARIO_RISORSE.csv",
        resources
    )

    print("\n" + "=" * 72)
    print("STEP 17.7 - AUDIT COMPLETATO")
    print("=" * 72)

    print("Small RNA:", len(summary))
    print("Allineamenti riportati:", len(examples))
    print("Risorse individuate:", len(resources))
    print("Possibili FASTA genomici:", len(genome_candidates))
    print("File annotazione:", len(annotations))

    print("\nCartella risultati:", output)

    print("\n[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun training")
    print("[OK] Nessun download")
    print("[OK] Nessun dato originale modificato")

    print(
        "\n[NOTA] L'audit non stabilisce "
        "l'origine genomica degli sRNA."
    )


if __name__ == "__main__":
    main()