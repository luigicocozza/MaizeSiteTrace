#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 17.15
CARATTERIZZAZIONE DELLE COPIE GENOMICHE RIPETUTE

Input:
    CONTESTI_GENOMICI_SRNA.csv (STEP 17.14)
    RAPPORTO_STEP17_14.json

Output:
    GRUPPI_CONTESTI_IDENTICI.csv
    DISTRIBUZIONE_NUMERO_COPIE.csv
    LOCI_DEI_GRUPPI_RIPETUTI.csv
    SINTESI_PER_SRNA.csv
    RAPPORTO_STEP17_15.json

Nessun download o nuovo mapping.
Nessuna modifica a GUI, motore, ML/GNN o dataset.
"""

from pathlib import Path
from datetime import datetime, timezone
from collections import Counter, defaultdict
import csv
import hashlib
import json


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

RUNS = (
    ROOT / "target_prediction_v33_vscode" / "runs"
)

EXPECTED_LOCI = 646
EXPECTED_RNA = 5


def section(n, title):
    print()
    print("=" * 70)
    print(f"[{n}] {title}")
    print("=" * 70)


def latest_step14():
    paths = sorted(
        p for p in RUNS.glob(
            "STEP17_14_SEQUENZE_RIPETUTE_*"
        )
        if p.is_dir()
        and (p / "CONTESTI_GENOMICI_SRNA.csv").is_file()
        and (p / "RAPPORTO_STEP17_14.json").is_file()
    )

    if not paths:
        raise RuntimeError(
            "Risultati STEP 17.14 non disponibili"
        )

    return paths[-1]


def load_csv(path):
    with path.open(
        "r", encoding="utf-8", newline=""
    ) as f:
        return list(csv.DictReader(f))


def save_csv(path, rows, columns):
    with path.open(
        "w", encoding="utf-8", newline=""
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=columns,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def is_chromosome(name):
    return (
        name.startswith("chr")
        and name[3:] in {
            str(i) for i in range(1, 11)
        }
    )


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.15")
    print("CARATTERIZZAZIONE COPIE GENOMICHE")
    print("=" * 70)

    section(1, "RECUPERO RISULTATI STEP 17.14")

    previous = latest_step14()
    input_csv = (
        previous / "CONTESTI_GENOMICI_SRNA.csv"
    )

    report14 = json.loads(
        (previous / "RAPPORTO_STEP17_14.json")
        .read_text(encoding="utf-8")
    )

    rows = load_csv(input_csv)

    if len(rows) != EXPECTED_LOCI:
        raise RuntimeError(
            "Conteggio loci inatteso: "
            + str(len(rows))
        )

    if report14.get("loci_verificati") != EXPECTED_LOCI:
        raise RuntimeError(
            "Report STEP 17.14 non coerente"
        )

    names = {r["srna"] for r in rows}

    if len(names) != EXPECTED_RNA:
        raise RuntimeError(
            "Numero sRNA inatteso"
        )

    groups = defaultdict(list)
    locus_keys = set()

    for row in rows:

        name = row["srna"]
        sequence = row["context_sequence"].upper()

        start = int(row["start_1based"])
        end = int(row["end_1based"])
        length = int(row["context_length"])
        expected = int(row["expected_context_length"])

        if row["sequence_verified"].lower() != "true":
            raise RuntimeError(
                "Locus non verificato: " + name
            )

        if len(sequence) != length:
            raise RuntimeError(
                "Lunghezza contesto incoerente"
            )

        if row["strand"] not in ("+", "-"):
            raise RuntimeError(
                "Strand non valido"
            )

        if start < 1 or end < start:
            raise RuntimeError(
                "Coordinate non valide"
            )

        key = (
            name,
            row["chromosome"],
            start,
            end,
            row["strand"],
        )

        if key in locus_keys:
            raise RuntimeError(
                "Locus duplicato: " + str(key)
            )

        locus_keys.add(key)

        # Non confrontiamo contesti incompleti
        # con quelli completi.
        if length == expected:
            groups[(name, sequence)].append(row)

    section(2, "ANALISI GRUPPI DI COPIE")

    output_groups = []
    output_loci = []
    copy_histogram = []
    summary = []

    group_number = 0

    for name in sorted(names):

        these_groups = [
            (sequence, members)
            for (rna, sequence), members in groups.items()
            if rna == name
        ]

        these_groups.sort(
            key=lambda item: (
                -len(item[1]),
                item[0],
            )
        )

        local_hist = Counter()
        repeated_loci = 0
        chromosome_loci = 0
        scaffold_loci = 0

        for sequence, members in these_groups:

            group_number += 1
            copies = len(members)

            local_hist[copies] += 1

            if copies > 1:
                repeated_loci += copies

            chromosome_members = [
                r for r in members
                if is_chromosome(r["chromosome"])
            ]

            scaffold_members = [
                r for r in members
                if r["chromosome"].startswith("scaf_")
            ]

            chromosome_loci += len(chromosome_members)
            scaffold_loci += len(scaffold_members)

            contigs = sorted({
                r["chromosome"] for r in members
            })

            strands = sorted({
                r["strand"] for r in members
            })

            digest = hashlib.sha256(
                sequence.encode("ascii")
            ).hexdigest()

            group_id = f"CTX_{group_number:05d}"

            output_groups.append({
                "group_id": group_id,
                "srna": name,
                "copies": copies,
                "chromosome_copies": len(
                    chromosome_members
                ),
                "scaffold_copies": len(
                    scaffold_members
                ),
                "distinct_contigs": len(contigs),
                "contigs": ";".join(contigs),
                "strands": ";".join(strands),
                "context_length_nt": len(sequence),
                "context_sha256": digest,
                "context_sequence": sequence,
            })

            for member in members:

                output_loci.append({
                    "group_id": group_id,
                    "srna": name,
                    "copies_in_group": copies,
                    "chromosome": member["chromosome"],
                    "start_1based": member[
                        "start_1based"
                    ],
                    "end_1based": member[
                        "end_1based"
                    ],
                    "strand": member["strand"],
                    "te_overlap": member["te_overlap"],
                    "region": member["region"],
                    "context_sha256": digest,
                })

        for copies, number_of_groups in sorted(
            local_hist.items()
        ):
            copy_histogram.append({
                "srna": name,
                "copies_per_context": copies,
                "number_of_contexts": number_of_groups,
                "total_loci": copies * number_of_groups,
            })

        total_loci = sum(
            row["srna"] == name for row in rows
        )

        complete_loci = sum(
            len(members)
            for _, members in these_groups
        )

        largest_group = max(
            (len(m) for _, m in these_groups),
            default=0,
        )

        result = {
            "srna": name,
            "loci_totali": total_loci,
            "contesti_completi": complete_loci,
            "contesti_parziali": (
                total_loci - complete_loci
            ),
            "gruppi_distinti": len(these_groups),
            "gruppi_con_almeno_2_copie": sum(
                count for copies, count
                in local_hist.items()
                if copies >= 2
            ),
            "loci_in_contesti_ripetuti": repeated_loci,
            "loci_in_contesti_non_ripetuti": (
                complete_loci - repeated_loci
            ),
            "copie_gruppo_piu_abbondante": largest_group,
            "copie_cromosomiche_complete": chromosome_loci,
            "copie_scaffold_complete": scaffold_loci,
        }

        summary.append(result)

        print()
        print("Small RNA:", name)
        print("Loci:", total_loci)
        print("Contesti distinti:", len(these_groups))
        print(
            "Gruppi ripetuti:",
            result["gruppi_con_almeno_2_copie"]
        )
        print(
            "Loci con contesti ripetuti:",
            repeated_loci
        )
        print(
            "Copie del contesto piu abbondante:",
            largest_group
        )
        print(
            "Contesti incompleti:",
            result["contesti_parziali"]
        )

    section(3, "CONTROLLO DI COERENZA")

    complete_total = sum(
        row["contesti_completi"] for row in summary
    )

    if complete_total != len(output_loci):
        raise RuntimeError(
            "Conteggi dei contesti incoerenti"
        )

    if sum(
        row["loci_totali"] for row in summary
    ) != EXPECTED_LOCI:
        raise RuntimeError(
            "Sono stati persi loci"
        )

    if sum(
        int(row["copies"]) for row in output_groups
    ) != complete_total:
        raise RuntimeError(
            "Copie dei gruppi incoerenti"
        )

    # Verifica incrociata con STEP 17.14
    old_summary = {
        row["srna"]: row
        for row in report14["summary"]
    }

    for row in summary:

        old = old_summary[row["srna"]]

        if (
            row["loci_in_contesti_ripetuti"]
            != int(old["contesti_identici_tra_loci"])
        ):
            raise RuntimeError(
                "Numero ripetizioni non coerente "
                "con STEP 17.14: " + row["srna"]
            )

        if (
            row["gruppi_distinti"]
            != int(old["contesti_distinti"])
        ):
            raise RuntimeError(
                "Numero contesti distinti incoerente: "
                + row["srna"]
            )

    print("[OK] Conteggi conformi allo STEP 17.14")
    print("[OK] Gruppi e loci riconciliati")

    section(4, "SALVATAGGIO RISULTATI")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP17_15_COPIE_GENOMICHE_{timestamp}"
    )

    output.mkdir(parents=True, exist_ok=False)

    save_csv(
        output / "GRUPPI_CONTESTI_IDENTICI.csv",
        output_groups,
        [
            "group_id", "srna", "copies",
            "chromosome_copies", "scaffold_copies",
            "distinct_contigs", "contigs", "strands",
            "context_length_nt", "context_sha256",
            "context_sequence",
        ],
    )

    save_csv(
        output / "DISTRIBUZIONE_NUMERO_COPIE.csv",
        copy_histogram,
        [
            "srna", "copies_per_context",
            "number_of_contexts", "total_loci",
        ],
    )

    save_csv(
        output / "LOCI_DEI_GRUPPI_RIPETUTI.csv",
        [
            row for row in output_loci
            if row["copies_in_group"] > 1
        ],
        [
            "group_id", "srna", "copies_in_group",
            "chromosome", "start_1based", "end_1based",
            "strand", "te_overlap", "region",
            "context_sha256",
        ],
    )

    save_csv(
        output / "SINTESI_PER_SRNA.csv",
        summary,
        list(summary[0]),
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.15",
        "stato": "ANALISI_COMPLETATA",
        "fonte": str(input_csv),
        "loci_totali": EXPECTED_LOCI,
        "gruppi_distinti_completi": len(output_groups),
        "summary": summary,
        "limitazioni": [
            "I gruppi considerano solo identita esatta "
            "delle finestre genomiche complete.",
            "Il confronto e orientato rispetto agli sRNA.",
            "Contesti incompleti esclusi dai gruppi.",
            "La ripetizione non identifica automaticamente "
            "una famiglia di trasposoni.",
            "Non e dimostrata la biogenesi dello small RNA.",
            "Nessuna validazione sperimentale.",
        ],
        "gui_modified": False,
        "engine_modified": False,
        "training_performed": False,
    }

    report_path = output / "RAPPORTO_STEP17_15.json"

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print()
    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.15 TERMINATO")
    print("=" * 70)
    print("Loci esaminati:", EXPECTED_LOCI)
    print("Gruppi distinti:", len(output_groups))
    print("Risultati:", output)
    print("Report:", report_path)
    print()
    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun training")
    print("[OK] Dataset originali preservati")


if __name__ == "__main__":
    main()