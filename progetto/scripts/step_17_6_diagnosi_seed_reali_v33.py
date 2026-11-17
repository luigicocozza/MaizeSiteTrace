#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.6 - DIAGNOSI SEED SU TRASCRITTOMA REALE

Scansione non distruttiva dei 5 sRNA PRJNA419326
contro il trascrittoma B73.

- Conta le occorrenze delle ancore 8-mer
- Verifica le posizioni candidate complete
- Calcola mismatch senza gap
- Conserva esempi dei migliori allineamenti
- Non modifica il motore V3.3
- Non modifica la GUI V3.2
- Non esegue training
- Non attribuisce validazione biologica
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json

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

STEM = "Zea_mays.Zm-B73-REFERENCE-NAM-5.0"

TRANSCRIPTOME = (
    BASE / "references"
    / f"{STEM}.cdna.all.fa.gz"
)

MAPPING = (
    BASE / "references"
    / f"{STEM}.tx2gene.tsv"
)

ANCHOR_LENGTH = 8
MAX_DIAGNOSTIC_MISMATCHES = 10
MAX_EXAMPLES_PER_RNA = 10

DNA = set("ACGT")

COMPLEMENT = str.maketrans(
    "ACGT", "TGCA"
)


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(
            path, "rt", encoding="utf-8"
        )
    return path.open(
        "r", encoding="utf-8"
    )


def fasta(path):

    header = None
    parts = []

    with open_text(path) as handle:

        for raw in handle:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    yield (
                        header.split()[0],
                        "".join(parts)
                        .upper()
                        .replace("U", "T")
                    )

                header = line[1:]
                parts = []

            elif header is not None:
                parts.append(line)

    if header is not None:
        yield (
            header.split()[0],
            "".join(parts)
            .upper()
            .replace("U", "T")
        )


def revcomp(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


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


def load_mapping(path):

    mapping = {}

    with open_text(path) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        if not {
            "transcript_id", "gene_id"
        }.issubset(reader.fieldnames or []):

            raise RuntimeError(
                "Intestazione mapping non valida"
            )

        for row in reader:

            tx = row["transcript_id"].strip()
            gene = row["gene_id"].strip()

            if tx in mapping and mapping[tx] != gene:
                raise RuntimeError(
                    f"Mapping ambiguo: {tx}"
                )

            mapping[tx] = gene

    return mapping


def mismatch_count(a, b):

    return sum(
        x != y for x, y in zip(a, b)
    )


def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.6")
    print("DIAGNOSI ANCORE RNA-GENE SU DATI REALI")
    print("=" * 72)

    step16 = latest_step16()

    rna_file = (
        step16 / "input_5_srna_reali.fa"
    )

    for path in [
        rna_file, TRANSCRIPTOME, MAPPING
    ]:
        if not path.is_file():
            raise FileNotFoundError(path)

    rnas = dict(fasta(rna_file))

    if len(rnas) != 5:
        raise RuntimeError(
            f"Attesi 5 small RNA, trovati {len(rnas)}"
        )

    mapping = load_mapping(MAPPING)

    print("\n[1] INPUT")

    print("sRNA:", len(rnas))
    print("Mapping:", len(mapping))
    print("Trascrittoma:", TRANSCRIPTOME)

    # ============================================================
    # INDICE DELLE ANCORE
    # ============================================================

    print("\n[2] COSTRUZIONE INDICE 8-MER")

    seed_index = defaultdict(set)
    reverse_sequences = {}

    per_rna = {}

    for rna_id, sequence in rnas.items():

        target = revcomp(sequence)
        reverse_sequences[rna_id] = target

        length = len(target)

        offsets = sorted({
            0,
            (length - ANCHOR_LENGTH) // 2,
            length - ANCHOR_LENGTH
        })

        anchors = []

        for offset in offsets:

            anchor = target[
                offset:offset + ANCHOR_LENGTH
            ]

            seed_index[anchor].add(
                (rna_id, offset)
            )

            anchors.append({
                "offset": offset,
                "sequenza": anchor,
                "occorrenze": 0
            })

        per_rna[rna_id] = {
            "lunghezza": length,
            "ancore": anchors,
            "anchor_occurrences": 0,
            "candidate_positions": 0,
            "candidate_valid_positions": 0,
            "best_mismatches": None,
            "mismatch_distribution": Counter(),
            "examples": []
        }

        print(
            f"[RNA] {rna_id[:48]} "
            f"| {length} nt "
            f"| ancore={len(offsets)}"
        )

    # ============================================================
    # SCANSIONE REALE
    # ============================================================

    print("\n[3] SCANSIONE TRASCRITTOMA")

    transcript_count = 0
    unmapped = 0
    bases_scanned = 0

    for tx_id, transcript in fasta(TRANSCRIPTOME):

        transcript_count += 1
        bases_scanned += len(transcript)

        if tx_id not in mapping:
            unmapped += 1
            continue

        gene_id = mapping[tx_id]

        visited = set()

        if len(transcript) < ANCHOR_LENGTH:
            continue

        for index in range(
            len(transcript) - ANCHOR_LENGTH + 1
        ):

            anchor = transcript[
                index:index + ANCHOR_LENGTH
            ]

            matches = seed_index.get(anchor)

            if not matches:
                continue

            for rna_id, offset in matches:

                info = per_rna[rna_id]

                info["anchor_occurrences"] += 1

                for entry in info["ancore"]:
                    if (
                        entry["offset"] == offset
                        and entry["sequenza"] == anchor
                    ):
                        entry["occorrenze"] += 1
                        break

                start = index - offset

                target = reverse_sequences[rna_id]
                length = len(target)

                if (
                    start < 0
                    or start + length > len(transcript)
                ):
                    continue

                identifier = (
                    rna_id, start
                )

                if identifier in visited:
                    continue

                visited.add(identifier)

                info["candidate_positions"] += 1

                site = transcript[
                    start:start + length
                ]

                if set(site) - DNA:
                    continue

                info["candidate_valid_positions"] += 1

                mismatches = mismatch_count(
                    target, site
                )

                best = info["best_mismatches"]

                if (
                    best is None
                    or mismatches < best
                ):
                    info["best_mismatches"] = mismatches

                if (
                    mismatches
                    <= MAX_DIAGNOSTIC_MISMATCHES
                ):
                    info["mismatch_distribution"][
                        mismatches
                    ] += 1

                    example = {
                        "transcript_id": tx_id,
                        "gene_id": gene_id,
                        "position_1based": start + 1,
                        "mismatches": mismatches,
                        "target_site": site
                    }

                    info["examples"].append(example)

                    info["examples"].sort(
                        key=lambda row: (
                            row["mismatches"],
                            row["transcript_id"],
                            row["position_1based"]
                        )
                    )

                    del info["examples"][
                        MAX_EXAMPLES_PER_RNA:
                    ]

        if transcript_count % 10000 == 0:

            print(
                f"[INFO] Trascritti processati: "
                f"{transcript_count}",
                flush=True
            )

    # ============================================================
    # REPORT
    # ============================================================

    print("\n[4] RISULTATI PER SMALL RNA")

    rows = []

    for rna_id, info in per_rna.items():

        print("\nRNA:", rna_id)

        print(
            "Lunghezza:",
            info["lunghezza"]
        )

        print(
            "Occorrenze seed:",
            info["anchor_occurrences"]
        )

        print(
            "Posizioni candidate:",
            info["candidate_positions"]
        )

        print(
            "Posizioni valide:",
            info["candidate_valid_positions"]
        )

        print(
            "Miglior numero mismatch:",
            info["best_mismatches"]
        )

        for anchor in info["ancore"]:
            print(
                "  Seed:",
                anchor["sequenza"],
                "| offset:",
                anchor["offset"],
                "| occorrenze:",
                anchor["occorrenze"]
            )

        for mm in range(
            MAX_DIAGNOSTIC_MISMATCHES + 1
        ):
            count = info[
                "mismatch_distribution"
            ][mm]

            if count:
                print(
                    f"  mismatch={mm}: {count}"
                )

        rows.append({
            "rna_id": rna_id,
            "length_nt": info["lunghezza"],
            "seed_occurrences": info[
                "anchor_occurrences"
            ],
            "candidate_positions": info[
                "candidate_positions"
            ],
            "valid_positions": info[
                "candidate_valid_positions"
            ],
            "best_mismatches": info[
                "best_mismatches"
            ],
            "hits_le_2": sum(
                info["mismatch_distribution"][i]
                for i in range(3)
            ),
            "hits_le_3": sum(
                info["mismatch_distribution"][i]
                for i in range(4)
            ),
            "hits_le_5": sum(
                info["mismatch_distribution"][i]
                for i in range(6)
            ),
            "hits_le_10": sum(
                info["mismatch_distribution"][i]
                for i in range(11)
            )
        })

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_6_DIAGNOSI_SEED_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    json_file = output / "RAPPORTO_STEP17_6.json"
    csv_file = output / "DIAGNOSI_SEED_PER_RNA.csv"

    serializable = {}

    for rna_id, info in per_rna.items():
        serializable[rna_id] = {
            **info,
            "mismatch_distribution": dict(
                info["mismatch_distribution"]
            )
        }

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.6",
        "stato": "DIAGNOSI_COMPLETATA",
        "input_sha256": {
            "rna": sha256(rna_file),
            "transcriptome": sha256(TRANSCRIPTOME),
            "mapping": sha256(MAPPING)
        },
        "transcripts_scanned": transcript_count,
        "transcripts_without_mapping": unmapped,
        "bases_scanned": bases_scanned,
        "seed_length": ANCHOR_LENGTH,
        "maximum_diagnostic_mismatches": (
            MAX_DIAGNOSTIC_MISMATCHES
        ),
        "rna_results": serializable,
        "limitations": [
            "Solo allineamenti recuperati da seed 8-mer",
            "Nessun gap o bulge",
            "Nessuno scoring G:U",
            "Nessuna analisi di accessibilita",
            "Nessuna validazione biologica",
            "Non equivale ad una predizione esaustiva"
        ]
    }

    json_file.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    with csv_file.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys())
        )

        writer.writeheader()
        writer.writerows(rows)

    print("\n" + "=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.6 TERMINATO")
    print("=" * 72)

    print("Trascritti scansionati:", transcript_count)
    print("Trascritti senza mapping:", unmapped)
    print("Nucleotidi scansionati:", bases_scanned)

    print("\nReport JSON:", json_file)
    print("Report CSV:", csv_file)

    print("\n[OK] GUI invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun training eseguito")

    print(
        "\n[NOTA] I risultati sono corrispondenze "
        "diagnostiche di sequenza, "
        "non target biologici validati."
    )


if __name__ == "__main__":
    main()