#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.11 - ANNOTAZIONE GENOMICA DEI LOCI SMALL RNA

Input:
- LOCI_GENOMICI_SRNA.csv dello STEP 17.10
- GTF Ensembl Plants B73 NAM-5.0 release 63

Output:
- LOCI_ANNOTATI_STEP17_11.csv
- SINTESI_SRNA_STEP17_11.csv
- GENI_SOVRAPPOSTI_STEP17_11.csv
- RAPPORTO_STEP17_11.json

Caratteristiche:
- Analisi di tutti i loci disponibili
- Gestione chr1/1 e scaffold
- Intersezione con geni ed esoni
- Classificazione delle regioni
- Associazioni multiple preservate
- Nessun download e nessun training
- Nessuna modifica a GUI, motore e dati originali

Importante:
La sovrapposizione con un gene non dimostra
che tale gene sia il bersaglio dello small RNA.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"

GTF = (
    ROOT
    / "data"
    / "reference"
    / "zea_mays_ensembl_plants"
    / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.63.gtf.gz"
)

EXPECTED_ASSEMBLY = "Zm-B73-REFERENCE-NAM-5.0"

EXPECTED_SRNA = 5
EXPECTED_LOCI = 646

LOCUS_FIELDS = [
    "rna_id",
    "chromosome",
    "start_1based",
    "end_1based",
    "strand",
    "mismatches",
    "length_nt",
    "reference",
    "mapping_type",
    "validation_status",
]


# ============================================================
# 2. FUNZIONI GENERALI
# ============================================================

def section(number, title):
    print()
    print("=" * 72)
    print(f"[{number}] {title}")
    print("=" * 72)


def now_utc():
    return datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")


def open_text(path):
    if str(path).endswith(".gz"):
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


def normalize_contig(value):
    """
    Normalizza i cromosomi principali:
    chr1 -> 1
    chromosome_1 -> 1
    1 -> 1

    Mantiene gli scaffold con il proprio nome.
    """

    value = value.strip()

    match = re.fullmatch(
        r"(?:chr|chromosome_?)?([1-9]|10)",
        value,
        flags=re.IGNORECASE,
    )

    if match:
        return match.group(1)

    return value


def contig_type(value):
    normalized = normalize_contig(value)

    if normalized in {
        str(i) for i in range(1, 11)
    }:
        return "CROMOSOMA"

    if normalized.lower().startswith("scaf_"):
        return "SCAFFOLD"

    return "ALTRO_CONTIG"


def parse_attributes(text):
    """
    Analizza attributi GTF del tipo:
    gene_id "ID"; transcript_id "ID";
    """

    result = {}

    for key, value in re.findall(
        r'([A-Za-z_][A-Za-z0-9_]*)\s+"([^"]*)"',
        text,
    ):
        result[key] = value

    return result


def latest_step17_10():

    matches = sorted(
        [
            p for p in RUNS.glob(
                "STEP17_10_MAPPING_GENOMICO_*"
            )
            if p.is_dir()
            and (
                p / "LOCI_GENOMICI_SRNA.csv"
            ).is_file()
            and (
                p / "RAPPORTO_STEP17_10.json"
            ).is_file()
        ],
        reverse=True,
    )

    if not matches:
        raise RuntimeError(
            "Risultati STEP 17.10 non individuati"
        )

    return matches[0]


# ============================================================
# 3. LETTURA LOCI REALI
# ============================================================

def load_loci(path):

    loci = []
    seen = set()

    with path.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:

        reader = csv.DictReader(handle)

        missing = (
            set(LOCUS_FIELDS)
            - set(reader.fieldnames or [])
        )

        if missing:
            raise RuntimeError(
                "Colonne mancanti nel CSV loci: "
                + ", ".join(sorted(missing))
            )

        for line_number, row in enumerate(
            reader,
            start=2,
        ):

            chromosome = row["chromosome"].strip()
            start = int(row["start_1based"])
            end = int(row["end_1based"])
            strand = row["strand"].strip()
            length = int(row["length_nt"])
            mismatches = int(row["mismatches"])

            if not chromosome or start < 1:
                raise RuntimeError(
                    f"Coordinate non valide, riga {line_number}"
                )

            if end - start + 1 != length:
                raise RuntimeError(
                    f"Lunghezza non coerente, riga {line_number}"
                )

            if strand not in ("+", "-"):
                raise RuntimeError(
                    f"Strand non valido, riga {line_number}"
                )

            if mismatches != 0:
                raise RuntimeError(
                    "Sono attesi solo loci esatti"
                )

            if row["reference"] != EXPECTED_ASSEMBLY:
                raise RuntimeError(
                    "Assemblaggio inatteso, riga "
                    + str(line_number)
                )

            key = (
                row["rna_id"],
                chromosome,
                start,
                end,
                strand,
            )

            if key in seen:
                raise RuntimeError(
                    "Locus duplicato nel CSV: "
                    + str(key)
                )

            seen.add(key)

            loci.append({
                **row,
                "start_1based": start,
                "end_1based": end,
                "mismatches": mismatches,
                "length_nt": length,
                "contig_normalized": normalize_contig(
                    chromosome
                ),
                "contig_type": contig_type(chromosome),
            })

    return loci


# ============================================================
# 4. LETTURA ANNOTAZIONE GTF
# ============================================================

def load_gtf(path):

    print("[INFO] Lettura integrale GTF...")

    intervals = defaultdict(
        lambda: {
            "gene": [],
            "exon": [],
        }
    )

    features = Counter()
    contigs = Counter()

    parsed_lines = 0
    skipped_lines = 0

    genes_seen = set()

    with open_text(path) as handle:

        for raw in handle:

            if not raw.strip() or raw.startswith("#"):
                continue

            fields = raw.rstrip("\n").split("\t")

            if len(fields) != 9:
                skipped_lines += 1
                continue

            (
                seqname,
                source,
                feature,
                start_raw,
                end_raw,
                score,
                strand,
                frame,
                attributes_raw,
            ) = fields

            parsed_lines += 1

            features[feature] += 1

            contig = normalize_contig(seqname)
            contigs[contig] += 1

            if feature not in ("gene", "exon"):
                continue

            attrs = parse_attributes(attributes_raw)

            gene_id = attrs.get("gene_id", "")

            if not gene_id:
                skipped_lines += 1
                continue

            start = int(start_raw)
            end = int(end_raw)

            if start > end or start < 1:
                skipped_lines += 1
                continue

            if feature == "gene":
                genes_seen.add(gene_id)

            intervals[contig][feature].append(
                (
                    start,
                    end,
                    gene_id,
                    strand,
                )
            )

    print(
        "[OK] Record GTF letti:",
        parsed_lines
    )

    print(
        "[OK] Geni annotati distinti:",
        len(genes_seen)
    )

    print(
        "[OK] Contig con annotazioni:",
        len(contigs)
    )

    return intervals, {
        "record_letti": parsed_lines,
        "record_scartati": skipped_lines,
        "geni_distinti": len(genes_seen),
        "feature_counts": dict(features),
        "contig_counts": dict(contigs),
    }


# ============================================================
# 5. INDICE INTERVALLI
# ============================================================

def prepare_interval_index(intervals):
    """
    Ordinamento e prefix maximum end.

    Permette di individuare le sovrapposizioni
    senza confrontare ciascun locus con tutti
    gli intervalli dell'annotazione.
    """

    index = {}

    for contig, by_type in intervals.items():

        index[contig] = {}

        for feature, rows in by_type.items():

            ordered = sorted(
                rows,
                key=lambda row: (
                    row[0],
                    row[1],
                )
            )

            starts = []
            max_ends = []

            current_max = 0

            for row in ordered:

                starts.append(row[0])

                current_max = max(
                    current_max,
                    row[1]
                )

                max_ends.append(current_max)

            index[contig][feature] = {
                "rows": ordered,
                "starts": starts,
                "max_ends": max_ends,
            }

    return index


def find_overlaps(
    indexed_feature,
    query_start,
    query_end,
):
    """
    Restituisce gli intervalli che si
    sovrappongono al locus inclusivo 1-based.
    """

    if indexed_feature is None:
        return []

    rows = indexed_feature["rows"]
    starts = indexed_feature["starts"]
    max_ends = indexed_feature["max_ends"]

    if not rows:
        return []

    right = bisect.bisect_right(
        starts,
        query_end
    ) - 1

    result = []

    while right >= 0:

        if max_ends[right] < query_start:
            break

        row = rows[right]

        if row[1] >= query_start:
            result.append(row)

        right -= 1

    return result


# ============================================================
# 6. ANNOTAZIONE DEI LOCI
# ============================================================

def annotate_locus(locus, index):

    contig = locus["contig_normalized"]

    start = locus["start_1based"]
    end = locus["end_1based"]

    contig_index = index.get(contig)

    if contig_index is None:

        return {
            "region": "NO_GTF_ANNOTATION",
            "gene_ids": [],
            "exon_gene_ids": [],
            "gene_strands": [],
            "gene_overlap_count": 0,
            "exon_overlap_count": 0,
            "relative_strands": [],
            "annotation_available_on_contig": False,
        }

    gene_hits = find_overlaps(
        contig_index.get("gene"),
        start,
        end,
    )

    exon_hits = find_overlaps(
        contig_index.get("exon"),
        start,
        end,
    )

    gene_ids = sorted({
        row[2] for row in gene_hits
    })

    exon_gene_ids = sorted({
        row[2] for row in exon_hits
    })

    gene_strands = sorted({
        row[3] for row in gene_hits
    })

    relative_strands = sorted({
        (
            "SAME"
            if row[3] == locus["strand"]
            else "OPPOSITE"
        )
        for row in gene_hits
        if row[3] in ("+", "-")
    })

    if exon_hits:
        region = "EXONIC_OVERLAP"

    elif gene_hits:
        region = "INTRAGENIC_NONEXONIC"

    else:
        region = "INTERGENIC_ANNOTATED_CONTIG"

    return {
        "region": region,
        "gene_ids": gene_ids,
        "exon_gene_ids": exon_gene_ids,
        "gene_strands": gene_strands,
        "gene_overlap_count": len(gene_hits),
        "exon_overlap_count": len(exon_hits),
        "relative_strands": relative_strands,
        "annotation_available_on_contig": True,
    }


# ============================================================
# 7. SCRITTURA CSV
# ============================================================

def write_csv(path, rows, fields):

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# 8. MAIN
# ============================================================

def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.11")
    print("ANNOTAZIONE GENOMICA DEI LOCI SMALL RNA")
    print("=" * 72)

    if not GTF.is_file():
        raise SystemExit(
            f"[ERRORE] GTF mancante: {GTF}"
        )

    previous = latest_step17_10()

    loci_csv = (
        previous / "LOCI_GENOMICI_SRNA.csv"
    )

    previous_report = json.loads(
        (
            previous / "RAPPORTO_STEP17_10.json"
        ).read_text(encoding="utf-8")
    )

    # --------------------------------------------------------
    # 1. INPUT
    # --------------------------------------------------------

    section(1, "VERIFICA INPUT")

    print("Loci:", loci_csv)
    print("GTF:", GTF)

    loci = load_loci(loci_csv)

    unique_rna = {
        locus["rna_id"] for locus in loci
    }

    print("[OK] Loci caricati:", len(loci))
    print("[OK] sRNA distinti:", len(unique_rna))

    expected = previous_report.get(
        "loci_esatti_totali"
    )

    saved = previous_report.get(
        "loci_salvati"
    )

    if previous_report.get("csv_troncato"):
        raise RuntimeError(
            "Il CSV del mapping risulta troncato"
        )

    if len(loci) != expected or len(loci) != saved:
        raise RuntimeError(
            "Numero di loci incoerente "
            "con il report STEP 17.10"
        )

    if len(unique_rna) != EXPECTED_SRNA:
        raise RuntimeError(
            "Numero di sRNA non conforme"
        )

    if len(loci) != EXPECTED_LOCI:
        raise RuntimeError(
            "Attesi 646 loci, trovati "
            + str(len(loci))
        )

    # --------------------------------------------------------
    # 2. GTF
    # --------------------------------------------------------

    section(2, "LETTURA ANNOTAZIONE GTF")

    intervals, gtf_stats = load_gtf(GTF)

    index = prepare_interval_index(intervals)

    # --------------------------------------------------------
    # 3. ANNOTAZIONE
    # --------------------------------------------------------

    section(3, "INTERSEZIONE LOCI - GENI - ESONI")

    annotated = []

    gene_associations = []

    region_counts = Counter()
    contig_counts = Counter()

    per_rna = defaultdict(Counter)

    for locus in loci:

        annotation = annotate_locus(
            locus,
            index
        )

        region = annotation["region"]

        region_counts[region] += 1

        contig_counts[
            locus["contig_type"]
        ] += 1

        per_rna[locus["rna_id"]]["TOTAL"] += 1

        per_rna[locus["rna_id"]][region] += 1

        per_rna[locus["rna_id"]][
            locus["contig_type"]
        ] += 1

        record = {
            **locus,
            **annotation,
        }

        for field in [
            "gene_ids",
            "exon_gene_ids",
            "gene_strands",
            "relative_strands",
        ]:
            record[field] = ";".join(
                record[field]
            )

        annotated.append(record)

        for gene_id in annotation["gene_ids"]:

            gene_associations.append({
                "rna_id": locus["rna_id"],
                "chromosome": locus["chromosome"],
                "start_1based": locus["start_1based"],
                "end_1based": locus["end_1based"],
                "strand": locus["strand"],
                "gene_id": gene_id,
                "region": region,
                "is_exonic_overlap": (
                    gene_id in annotation[
                        "exon_gene_ids"
                    ]
                ),
            })

    print("\nDistribuzione delle regioni:")

    for region, count in sorted(
        region_counts.items()
    ):
        print(f"  {region}: {count}")

    print("\nDistribuzione contig:")

    for contig_type, count in sorted(
        contig_counts.items()
    ):
        print(f"  {contig_type}: {count}")

    # --------------------------------------------------------
    # 4. RIEPILOGO PER RNA
    # --------------------------------------------------------

    section(4, "RIEPILOGO PER SMALL RNA")

    summary = []

    for rna_id in sorted(unique_rna):

        counts = per_rna[rna_id]

        print()
        print("RNA:", rna_id)
        print("Loci:", counts["TOTAL"])

        for region in sorted(region_counts):
            print(
                f"  {region}:",
                counts[region]
            )

        summary.append({
            "rna_id": rna_id,
            "loci_totali": counts["TOTAL"],
            "loci_cromosomici": counts["CROMOSOMA"],
            "loci_scaffold": counts["SCAFFOLD"],
            "exonic_overlap": counts[
                "EXONIC_OVERLAP"
            ],
            "intragenic_nonexonic": counts[
                "INTRAGENIC_NONEXONIC"
            ],
            "intergenic_annotated": counts[
                "INTERGENIC_ANNOTATED_CONTIG"
            ],
            "no_gtf_annotation": counts[
                "NO_GTF_ANNOTATION"
            ],
        })

    # --------------------------------------------------------
    # 5. SALVATAGGIO
    # --------------------------------------------------------

    section(5, "GENERAZIONE RISULTATI")

    timestamp = now_utc()

    output = (
        RUNS
        / f"STEP17_11_ANNOTAZIONE_GTF_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    annotated_csv = (
        output / "LOCI_ANNOTATI_STEP17_11.csv"
    )

    summary_csv = (
        output / "SINTESI_SRNA_STEP17_11.csv"
    )

    gene_csv = (
        output / "GENI_SOVRAPPOSTI_STEP17_11.csv"
    )

    write_csv(
        annotated_csv,
        annotated,
        list(annotated[0].keys()),
    )

    write_csv(
        summary_csv,
        summary,
        list(summary[0].keys()),
    )

    write_csv(
        gene_csv,
        gene_associations,
        [
            "rna_id",
            "chromosome",
            "start_1based",
            "end_1based",
            "strand",
            "gene_id",
            "region",
            "is_exonic_overlap",
        ],
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.11",
        "stato": "ANNOTAZIONE_COMPLETATA",
        "assembly": EXPECTED_ASSEMBLY,
        "input_step17_10": str(loci_csv),
        "gtf": str(GTF),
        "gtf_statistics": gtf_stats,
        "loci_totali": len(loci),
        "small_rna": len(unique_rna),
        "region_counts": dict(region_counts),
        "contig_counts": dict(contig_counts),
        "gene_overlap_associations": len(
            gene_associations
        ),
        "summary_per_rna": summary,
        "limitations": [
            "Le coordinate genomiche indicano "
            "corrispondenze di sequenza.",
            "Il GTF non descrive necessariamente "
            "tutte le regioni ripetitive.",
            "INTRAGENIC_NONEXONIC non equivale "
            "sempre a un introne.",
            "Le regioni senza annotazione GTF "
            "non vengono definite intergeniche.",
            "La sovrapposizione con un gene "
            "non dimostra che sia un target.",
            "Non sono stati eseguiti modelli GNN "
            "o analisi di espressione WW/DS.",
        ],
        "gui_modified": False,
        "engine_modified": False,
        "training_performed": False,
        "biological_validation": False,
    }

    json_path = (
        output / "RAPPORTO_STEP17_11.json"
    )

    json_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8",
    )

    # --------------------------------------------------------
    # 6. CONTROLLO FINALE
    # --------------------------------------------------------

    section(6, "CONTROLLI FINALI")

    if sum(region_counts.values()) != len(loci):
        raise RuntimeError(
            "Conteggio delle annotazioni incoerente"
        )

    if sum(
        row["loci_totali"] for row in summary
    ) != len(loci):
        raise RuntimeError(
            "Conteggio per RNA incoerente"
        )

    print("[OK] Tutti i loci classificati")
    print("[OK] Conteggi coerenti")
    print("[OK] Associazioni geniche preservate")

    print()
    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.11 TERMINATO")
    print("=" * 72)

    print("Loci analizzati:", len(loci))
    print("Small RNA:", len(unique_rna))
    print(
        "Associazioni locus-gene:",
        len(gene_associations)
    )

    print("\nCartella risultati:")
    print(output)

    print("\nCSV annotazione:")
    print(annotated_csv)

    print("\nCSV riepilogo:")
    print(summary_csv)

    print("\nCSV associazioni geniche:")
    print(gene_csv)

    print("\nReport JSON:")
    print(json_path)

    print("\n[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Dati originali preservati")
    print("[OK] Nessun training")

    print(
        "\n[NOTA] Sovrapposizioni genomiche "
        "esplorative, non target biologici validati."
    )


if __name__ == "__main__":
    main()