#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 17.16

Audit delle copie genomiche dominanti:
- Riutilizza i dati dello STEP 17.15
- Verifica gruppi e coordinate
- Identifica i gruppi dominanti
- Analizza cromosomi, scaffold e strand
- Analizza distanza tra loci dello stesso gruppo
- Verifica sequenze identiche e basi ambigue
- Produce report CSV e JSON

Nessun download, mapping o training.
GUI V3.2 e motore V3.3 non modificati.
"""

from pathlib import Path
from datetime import datetime, timezone
from collections import Counter, defaultdict
import csv
import json
import hashlib


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

RUNS = ROOT / "target_prediction_v33_vscode" / "runs"

TARGET_SRNA = "sRNA_00000037"
EXPECTED_TARGET_LOCI = 587
EXPECTED_TOTAL_LOCI = 646


def latest_run(prefix, filename):
    candidates = sorted(
        (
            p for p in RUNS.glob(prefix + "*")
            if p.is_dir() and (p / filename).is_file()
        ),
        reverse=True,
    )
    if not candidates:
        raise RuntimeError("Risultati mancanti: " + prefix)
    return candidates[0]


def read_csv(path):
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def save_csv(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=columns, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def contig_class(name):
    if name.startswith("chr") and name[3:] in {
        str(i) for i in range(1, 11)
    }:
        return "CROMOSOMA"
    if name.startswith("scaf_"):
        return "SCAFFOLD"
    return "ALTRO"


def main():
    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.16")
    print("AUDIT DELLE COPIE GENOMICHE DOMINANTI")
    print("=" * 72)

    step15 = latest_run(
        "STEP17_15_COPIE_GENOMICHE_",
        "GRUPPI_CONTESTI_IDENTICI.csv",
    )

    groups = read_csv(
        step15 / "GRUPPI_CONTESTI_IDENTICI.csv"
    )

    # Lo STEP 17.15 include soltanto i loci dei gruppi
    # ripetuti nel relativo CSV. Recuperiamo tutti i
    # 646 loci dallo STEP 17.14.
    step14 = latest_run(
        "STEP17_14_SEQUENZE_RIPETUTE_",
        "CONTESTI_GENOMICI_SRNA.csv",
    )
    loci = read_csv(
        step14 / "CONTESTI_GENOMICI_SRNA.csv"
    )

    if len(loci) != EXPECTED_TOTAL_LOCI:
        raise RuntimeError("Numero loci STEP 17.14 incoerente")

    print("\n[1] VERIFICA INPUT")
    print("Gruppi distinti:", len(groups))
    print("Loci analizzati:", len(loci))

    by_srna = Counter(r["srna"] for r in loci)

    if by_srna[TARGET_SRNA] != EXPECTED_TARGET_LOCI:
        raise RuntimeError("Numero loci sRNA_00000037 inatteso")

    group_lookup = {}
    group_rows = defaultdict(list)
    locus_keys = set()

    for group in groups:
        key = (group["srna"], group["context_sequence"])
        if key in group_lookup:
            raise RuntimeError("Gruppo duplicato")
        if int(group["copies"]) < 1:
            raise RuntimeError("Conteggio copie non valido")
        if len(group["context_sequence"]) != int(
            group["context_length_nt"]
        ):
            raise RuntimeError("Lunghezza contesto incoerente")
        digest = hashlib.sha256(
            group["context_sequence"].encode("ascii")
        ).hexdigest()
        if digest != group["context_sha256"]:
            raise RuntimeError("SHA256 contesto non conforme")
        group_lookup[key] = group

    for row in loci:
        start = int(row["start_1based"])
        end = int(row["end_1based"])
        if start < 1 or end < start:
            raise RuntimeError("Coordinate non valide")

        unique_key = (
            row["srna"], row["chromosome"],
            start, end, row["strand"],
        )
        if unique_key in locus_keys:
            raise RuntimeError("Locus duplicato")
        locus_keys.add(unique_key)

        if row["context_length"] != row["expected_context_length"]:
            raise RuntimeError(
                "Contesto parziale non previsto in questo audit"
            )

        key = (row["srna"], row["context_sequence"])
        if key not in group_lookup:
            raise RuntimeError("Contesto non presente nei gruppi")
        group_rows[key].append(row)

    for key, group in group_lookup.items():
        if len(group_rows[key]) != int(group["copies"]):
            raise RuntimeError(
                "Conteggio copie non coerente: "
                + group["group_id"]
            )

    print("[OK] Gruppi, checksum e coordinate coerenti")

    print("\n[2] GRUPPI DOMINANTI")

    dominant_summary = []
    distributions = []
    distance_rows = []

    names = sorted(by_srna)

    for name in names:
        candidates = [
            g for g in groups if g["srna"] == name
        ]

        dominant = max(
            candidates, key=lambda g: int(g["copies"])
        )

        key = (name, dominant["context_sequence"])
        members = group_rows[key]

        strand_counts = Counter(
            r["strand"] for r in members
        )
        contig_counts = Counter(
            r["chromosome"] for r in members
        )
        class_counts = Counter(
            contig_class(r["chromosome"])
            for r in members
        )

        seq = dominant["context_sequence"]
        ambiguous = sum(base not in "ACGT" for base in seq)
        gc = (
            100 * (seq.count("G") + seq.count("C")) / len(seq)
        )

        info = {
            "srna": name,
            "group_id": dominant["group_id"],
            "total_loci": by_srna[name],
            "dominant_copies": len(members),
            "dominant_fraction_percent": round(
                100 * len(members) / by_srna[name], 3
            ),
            "context_length_nt": len(seq),
            "context_sha256": dominant["context_sha256"],
            "ambiguous_bases": ambiguous,
            "gc_percent": round(gc, 3),
            "chromosomal_copies": class_counts["CROMOSOMA"],
            "scaffold_copies": class_counts["SCAFFOLD"],
            "other_copies": class_counts["ALTRO"],
            "plus_strand": strand_counts["+"],
            "minus_strand": strand_counts["-"],
            "distinct_contigs": len(contig_counts),
            "top_contig": (
                contig_counts.most_common(1)[0][0]
            ),
            "top_contig_copies": (
                contig_counts.most_common(1)[0][1]
            ),
        }

        dominant_summary.append(info)

        print(
            f"{name}: {len(members)} copie / "
            f"{by_srna[name]} loci "
            f"({info['dominant_fraction_percent']}%)"
        )

        for contig, number in sorted(
            contig_counts.items(),
            key=lambda x: (-x[1], x[0]),
        ):
            distributions.append({
                "srna": name,
                "group_id": dominant["group_id"],
                "contig": contig,
                "contig_type": contig_class(contig),
                "copies": number,
                "plus_strand": sum(
                    r["strand"] == "+"
                    for r in members
                    if r["chromosome"] == contig
                ),
                "minus_strand": sum(
                    r["strand"] == "-"
                    for r in members
                    if r["chromosome"] == contig
                ),
            })

        per_contig = defaultdict(list)
        for row in members:
            per_contig[row["chromosome"]].append(row)

        for contig, contig_loci in per_contig.items():
            ordered = sorted(
                contig_loci,
                key=lambda x: int(x["start_1based"])
            )

            for previous, current in zip(
                ordered, ordered[1:]
            ):
                gap = (
                    int(current["start_1based"])
                    - int(previous["end_1based"]) - 1
                )

                distance_rows.append({
                    "srna": name,
                    "group_id": dominant["group_id"],
                    "contig": contig,
                    "previous_start": previous["start_1based"],
                    "previous_end": previous["end_1based"],
                    "current_start": current["start_1based"],
                    "current_end": current["end_1based"],
                    "gap_nt": gap,
                    "overlap": gap < 0,
                })

    print("\n[3] FOCUS sRNA_00000037")

    target = next(
        r for r in dominant_summary
        if r["srna"] == TARGET_SRNA
    )

    for key, value in target.items():
        print(f"{key}: {value}")

    target_distances = [
        int(r["gap_nt"])
        for r in distance_rows
        if r["srna"] == TARGET_SRNA
    ]

    print(
        "Coppie di loci successivi sullo stesso contig:",
        len(target_distances)
    )

    if target_distances:
        print("Distanza minima (gap):", min(target_distances))
        print("Distanza massima (gap):", max(target_distances))

    print("\n[4] SALVATAGGIO")

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS / f"STEP17_16_AUDIT_COPIE_{stamp}"
    )
    output.mkdir(parents=True, exist_ok=False)

    save_csv(
        output / "SINTESI_GRUPPI_DOMINANTI.csv",
        dominant_summary,
        list(dominant_summary[0]),
    )
    save_csv(
        output / "DISTRIBUZIONE_COPIE_PER_CONTIG.csv",
        distributions,
        list(distributions[0]),
    )
    save_csv(
        output / "DISTANZE_COPIE_SULLO_STESSO_CONTIG.csv",
        distance_rows,
        [
            "srna", "group_id", "contig",
            "previous_start", "previous_end",
            "current_start", "current_end",
            "gap_nt", "overlap",
        ],
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.16",
        "stato": "AUDIT_COMPLETATO",
        "input_step17_14": str(step14),
        "input_step17_15": str(step15),
        "loci_verificati": len(loci),
        "gruppi_verificati": len(groups),
        "sintesi_gruppi_dominanti": dominant_summary,
        "focus_srna_00000037": target,
        "interpretazione": (
            "Le copie identiche dimostrano ripetitivita "
            "di sequenza nel genoma di riferimento, "
            "non la classe biologica dell'elemento "
            "o il locus effettivo di biogenesi."
        ),
        "gui_modified": False,
        "engine_modified": False,
        "training_performed": False,
        "biological_validation": False,
    }

    report_path = output / "RAPPORTO_STEP17_16.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("\nSTEP 17.16 TERMINATO")
    print("Report:", report_path)
    print("Risultati:", output)
    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun training")
    print("[OK] Dataset originali preservati")


if __name__ == "__main__":
    main()