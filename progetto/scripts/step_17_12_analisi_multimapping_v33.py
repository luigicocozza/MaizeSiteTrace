#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.12 - MULTIMAPPING E LOCI GENOMICI CONDIVISI

INPUT:
- LOCI_ANNOTATI_STEP17_11.csv
- GENI_SOVRAPPOSTI_STEP17_11.csv
- RAPPORTO_STEP17_11.json

ANALISI:
1. Integrita dei 646 loci
2. Distribuzione per RNA e cromosoma
3. Multimapping su cromosomi e scaffold
4. Loci sovrapposti tra RNA differenti
5. Cluster genomici entro 100 nucleotidi
6. Associazioni con geni annotati
7. Analisi particolare dello sRNA_00000037
8. Produzione CSV e JSON

SICUREZZA:
- Nessun file originale modificato
- GUI V3.2 invariata
- Motore V3.3 invariato
- Nessun training
- Nessuna attribuzione automatica di target
- Nessuna identificazione presunta di trasposoni
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

EXPECTED_LOCI = 646
EXPECTED_SRNA = 5

CLUSTER_DISTANCE_NT = 100


def section(number, title):
    print()
    print("=" * 70)
    print(f"[{number}] {title}")
    print("=" * 70)


def latest_step11():

    folders = sorted(
        [
            p for p in RUNS.glob(
                "STEP17_11_ANNOTAZIONE_GTF_*"
            )
            if p.is_dir()
        ],
        reverse=True
    )

    if not folders:
        raise RuntimeError(
            "STEP 17.11 non trovato"
        )

    return folders[0]


def read_csv(path):

    if not path.is_file():
        raise FileNotFoundError(path)

    with path.open(
        "r",
        newline="",
        encoding="utf-8"
    ) as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):

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


def short_id(rna_id):

    for component in rna_id.split("|"):
        if component.startswith("feature="):
            return component.split("=", 1)[1]

    return rna_id.split("|")[0]


def chromosome_type(chromosome):

    name = chromosome.lower()

    if name.startswith("chr") and name[3:].isdigit():
        return "CROMOSOMA"

    if name.startswith("scaf_"):
        return "SCAFFOLD"

    return "ALTRO"


def normalize_loci(rows):

    result = []
    seen = set()

    for row in rows:

        row = dict(row)

        row["start_1based"] = int(
            row["start_1based"]
        )

        row["end_1based"] = int(
            row["end_1based"]
        )

        row["length_nt"] = int(
            row["length_nt"]
        )

        if (
            row["start_1based"] < 1
            or row["end_1based"] < row["start_1based"]
        ):
            raise RuntimeError(
                "Coordinate genomiche non valide"
            )

        if (
            row["end_1based"]
            - row["start_1based"] + 1
            != row["length_nt"]
        ):
            raise RuntimeError(
                "Lunghezza del locus non coerente"
            )

        if row["strand"] not in ("+", "-"):
            raise RuntimeError(
                "Orientamento genomico non valido"
            )

        key = (
            row["rna_id"],
            row["chromosome"],
            row["start_1based"],
            row["end_1based"],
            row["strand"]
        )

        if key in seen:
            raise RuntimeError(
                "Locus duplicato: " + str(key)
            )

        seen.add(key)

        row["srna_short"] = short_id(
            row["rna_id"]
        )

        row["genomic_type"] = chromosome_type(
            row["chromosome"]
        )

        result.append(row)

    return result


def build_clusters(loci, distance=0):
    """
    Raggruppa intervalli sullo stesso contig.

    distance=0:
        intervalli sovrapposti.

    distance=100:
        intervalli sovrapposti o separati
        da massimo 100 nucleotidi.

    Gli strand sono mantenuti come informazione,
    ma NON separano i cluster.
    """

    by_contig = defaultdict(list)

    for locus in loci:
        by_contig[locus["chromosome"]].append(
            locus
        )

    clusters = []
    cluster_number = 0

    for chromosome in sorted(by_contig):

        rows = sorted(
            by_contig[chromosome],
            key=lambda r: (
                r["start_1based"],
                r["end_1based"],
                r["srna_short"]
            )
        )

        current = []
        current_end = None

        def save_cluster(members):

            nonlocal cluster_number

            if not members:
                return

            cluster_number += 1

            rnas = sorted({
                r["srna_short"] for r in members
            })

            strands = sorted({
                r["strand"] for r in members
            })

            clusters.append({
                "cluster_id": (
                    f"CLUSTER_{cluster_number:06d}"
                ),
                "chromosome": chromosome,
                "start_1based": min(
                    r["start_1based"]
                    for r in members
                ),
                "end_1based": max(
                    r["end_1based"]
                    for r in members
                ),
                "loci_count": len(members),
                "distinct_srna": len(rnas),
                "srna_ids": ";".join(rnas),
                "strands": ";".join(strands),
                "genomic_type": chromosome_type(
                    chromosome
                ),
            })

        for row in rows:

            if not current:

                current = [row]
                current_end = row["end_1based"]
                continue

            gap = (
                row["start_1based"]
                - current_end - 1
            )

            if gap <= distance:

                current.append(row)
                current_end = max(
                    current_end,
                    row["end_1based"]
                )

            else:

                save_cluster(current)
                current = [row]
                current_end = row["end_1based"]

        save_cluster(current)

    return clusters


def analyze_genes(rows, loci):

    keys = {
        (
            r["rna_id"],
            r["chromosome"],
            r["start_1based"],
            r["end_1based"],
            r["strand"]
        )
        for r in loci
    }

    result = []

    for row in rows:

        row = dict(row)

        row["start_1based"] = int(
            row["start_1based"]
        )

        row["end_1based"] = int(
            row["end_1based"]
        )

        key = (
            row["rna_id"],
            row["chromosome"],
            row["start_1based"],
            row["end_1based"],
            row["strand"]
        )

        if key not in keys:
            raise RuntimeError(
                "Associazione genica senza locus: "
                + str(key)
            )

        row["srna_short"] = short_id(
            row["rna_id"]
        )

        result.append(row)

    return result


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.12")
    print("MULTIMAPPING E LOCI CONDIVISI")
    print("=" * 70)

    previous = latest_step11()

    loci_file = (
        previous / "LOCI_ANNOTATI_STEP17_11.csv"
    )

    genes_file = (
        previous / "GENI_SOVRAPPOSTI_STEP17_11.csv"
    )

    report_file = (
        previous / "RAPPORTO_STEP17_11.json"
    )

    # ========================================================
    # 1. INPUT
    # ========================================================

    section(1, "LETTURA INPUT STEP 17.11")

    original_report = json.loads(
        report_file.read_text(
            encoding="utf-8"
        )
    )

    loci = normalize_loci(
        read_csv(loci_file)
    )

    genes = analyze_genes(
        read_csv(genes_file),
        loci
    )

    rna_ids = sorted({
        row["rna_id"] for row in loci
    })

    if len(loci) != EXPECTED_LOCI:
        raise RuntimeError(
            "Numero loci inatteso: "
            + str(len(loci))
        )

    if len(rna_ids) != EXPECTED_SRNA:
        raise RuntimeError(
            "Numero small RNA inatteso"
        )

    if len(loci) != original_report["loci_totali"]:
        raise RuntimeError(
            "Incoerenza con report STEP 17.11"
        )

    print("[OK] Loci:", len(loci))
    print("[OK] Small RNA:", len(rna_ids))
    print("[OK] Associazioni geniche:", len(genes))

    # ========================================================
    # 2. DISTRIBUZIONE
    # ========================================================

    section(2, "DISTRIBUZIONE MULTIMAPPING")

    per_rna = Counter(
        row["srna_short"] for row in loci
    )

    per_chromosome = Counter(
        row["chromosome"] for row in loci
    )

    per_strand = Counter(
        row["strand"] for row in loci
    )

    distribution = defaultdict(Counter)

    for row in loci:
        distribution[row["srna_short"]][
            row["genomic_type"]
        ] += 1

    summary = []

    for name, count in sorted(per_rna.items()):

        chr_count = distribution[name][
            "CROMOSOMA"
        ]

        scaffold_count = distribution[name][
            "SCAFFOLD"
        ]

        summary.append({
            "srna": name,
            "loci_totali": count,
            "loci_cromosomici": chr_count,
            "loci_scaffold": scaffold_count,
            "loci_altri": distribution[name][
                "ALTRO"
            ],
            "percentuale_scaffold": round(
                100 * scaffold_count / count,
                3
            ),
            "multimapping": (
                "SI" if count > 1 else "NO"
            )
        })

        print()
        print("Small RNA:", name)
        print("Loci:", count)
        print("Cromosomi:", chr_count)
        print("Scaffold:", scaffold_count)

    print("\nPrincipali contig per numero di loci:")

    for contig, count in per_chromosome.most_common(12):
        print(f"  {contig}: {count}")

    # ========================================================
    # 3. CLUSTER CON SOVRAPPOSIZIONE REALE
    # ========================================================

    section(3, "LOCI SOVRAPPOSTI")

    overlap_clusters = build_clusters(
        loci,
        distance=0
    )

    shared = [
        row for row in overlap_clusters
        if row["distinct_srna"] > 1
    ]

    print(
        "Cluster con sovrapposizione:",
        len(overlap_clusters)
    )

    print(
        "Cluster condivisi tra RNA distinti:",
        len(shared)
    )

    for row in shared[:20]:

        print(
            row["chromosome"],
            f"{row['start_1based']}-"
            f"{row['end_1based']}",
            "| RNA:",
            row["srna_ids"],
        )

    # ========================================================
    # 4. CLUSTER DI VICINANZA
    # ========================================================

    section(4, "CLUSTER GENOMICI ENTRO 100 NT")

    near_clusters = build_clusters(
        loci,
        distance=CLUSTER_DISTANCE_NT
    )

    near_shared = [
        row for row in near_clusters
        if row["distinct_srna"] > 1
    ]

    print(
        "Cluster di vicinanza:",
        len(near_clusters)
    )

    print(
        "Cluster con RNA distinti:",
        len(near_shared)
    )

    # ========================================================
    # 5. ASSOCIAZIONI GENICHE
    # ========================================================

    section(5, "GENI SOVRAPPOSTI")

    distinct_genes = sorted({
        row["gene_id"] for row in genes
    })

    print(
        "Associazioni locus-gene:",
        len(genes)
    )

    print(
        "Geni distinti:",
        len(distinct_genes)
    )

    for row in genes:

        print(
            row["srna_short"],
            "|",
            row["gene_id"],
            "|",
            row["chromosome"],
            "|",
            f"{row['start_1based']}-"
            f"{row['end_1based']}",
            "|",
            row["region"]
        )

    # ========================================================
    # 6. VALUTAZIONE ESPLORATIVA
    # ========================================================

    section(6, "VALUTAZIONE SCIENTIFICA")

    most_mapped = max(
        per_rna,
        key=per_rna.get
    )

    print(
        "Small RNA con piu loci:",
        most_mapped
    )

    print(
        "Numero loci:",
        per_rna[most_mapped]
    )

    print(
        "Loci condivisi in cluster:",
        len(shared)
    )

    print(
        "Associazioni geniche:",
        len(genes)
    )

    print(
        "\n[NOTA] Numerosi loci genomici "
        "suggeriscono un approfondimento "
        "delle sequenze ripetitive."
    )

    print(
        "[NOTA] Non e stata verificata "
        "l'appartenenza a elementi trasponibili."
    )

    print(
        "[NOTA] Le associazioni locus-gene "
        "non sono predizioni di targeting."
    )

    # ========================================================
    # 7. SALVATAGGIO
    # ========================================================

    section(7, "SALVATAGGIO REPORT")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_12_MULTIMAPPING_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    write_csv(
        output / "SINTESI_MULTIMAPPING.csv",
        summary,
        list(summary[0].keys())
    )

    cluster_fields = [
        "cluster_id",
        "chromosome",
        "start_1based",
        "end_1based",
        "loci_count",
        "distinct_srna",
        "srna_ids",
        "strands",
        "genomic_type",
    ]

    write_csv(
        output / "CLUSTER_SOVRAPPOSTI.csv",
        overlap_clusters,
        cluster_fields
    )

    write_csv(
        output / "CLUSTER_CONDIVISI.csv",
        shared,
        cluster_fields
    )

    write_csv(
        output / "CLUSTER_ENTRO_100NT.csv",
        near_clusters,
        cluster_fields
    )

    write_csv(
        output / "GENI_SOVRAPPOSTI_VERIFICATI.csv",
        genes,
        list(genes[0].keys()) if genes else [
            "rna_id",
            "chromosome",
            "start_1based",
            "end_1based",
            "strand",
            "gene_id",
            "region",
            "is_exonic_overlap",
            "srna_short",
        ]
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.12",
        "stato": "ANALISI_COMPLETATA",
        "input_step17_11": str(previous),
        "loci_totali": len(loci),
        "small_rna": len(rna_ids),
        "conteggi_per_rna": dict(per_rna),
        "conteggi_per_chromosome": dict(
            per_chromosome
        ),
        "conteggi_per_strand": dict(
            per_strand
        ),
        "multimapping_summary": summary,
        "overlap_clusters": len(
            overlap_clusters
        ),
        "shared_clusters": len(shared),
        "near_clusters_100nt": len(
            near_clusters
        ),
        "near_shared_clusters": len(
            near_shared
        ),
        "gene_associations": len(genes),
        "distinct_genes": distinct_genes,
        "most_multimapped_srna": most_mapped,
        "repeat_annotation_verified": False,
        "transposable_elements_verified": False,
        "biological_targets_validated": False,
        "limitations": [
            "Cluster di sovrapposizione non "
            "dimostrano origine comune.",
            "Cluster di vicinanza non "
            "identificano automaticamente loci di biogenesi.",
            "Il multimapping non identifica "
            "automaticamente elementi trasponibili.",
            "La sovrapposizione con geni "
            "non equivale a targeting.",
            "Non sono state consultate "
            "annotazioni RepeatMasker o TE.",
        ],
        "gui_modified": False,
        "engine_modified": False,
        "training_performed": False,
    }

    json_path = (
        output / "RAPPORTO_STEP17_12.json"
    )

    json_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print()
    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.12 TERMINATO")
    print("=" * 70)

    print("Loci:", len(loci))
    print("Small RNA:", len(rna_ids))
    print("Cluster condivisi:", len(shared))
    print("Cluster entro 100 nt:", len(near_clusters))
    print("Associazioni geniche:", len(genes))

    print("\nRisultati:")
    print(output)

    print("\n[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun training")
    print("[OK] Dati originali preservati")


if __name__ == "__main__":
    main()