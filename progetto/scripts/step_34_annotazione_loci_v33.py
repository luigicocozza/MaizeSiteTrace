#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 34

ANNOTAZIONE DEI LOCI GENOMICI DELLE SMALL RNA

MODALITA':

--audit
    Verifica STEP33, checksum, GTF e TE GFF3.
    Non modifica né genera file nel progetto.

--run
    Analizza le annotazioni genomiche e
    determina quali loci si sovrappongono
    a geni, elementi trasponibili o
    altre regioni ripetute.

NON MODIFICA:
    GUI V3.2 e V3.3
    Docker
    Modelli ML/GNN
    Dati scientifici originali
"""

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import re
import sys

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

GTF = (
    ROOT
    / "data/reference/zea_mays_ensembl_plants"
    / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.63.gtf.gz"
)

TE_GFF = (
    ROOT
    / "target_prediction_v33_vscode"
    / "references"
    / "Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz"
)

EXPECTED_LENGTHS = {
    "sRNA_00000037": 25,
    "sRNA_00000435": 24
}


# ============================================================
# 2. ESPRESSIONI REGOLARI
# ============================================================

GENE_ID_RE = re.compile(
    r'\bgene_id\s+"([^"]+)"'
)

GENE_BIOTYPE_RE = re.compile(
    r'\bgene_biotype\s+"([^"]+)"'
)

TE_RE = re.compile(
    r"transposon|transposable_element|helitron"
    r"|(?:^|[_-])(?:line|sine)(?:$|[_-])",
    re.IGNORECASE
)

OTHER_REPEAT_RE = re.compile(
    r"repeat|knob|satellite|"
    r"target_site_duplication|low_complexity",
    re.IGNORECASE
)


# ============================================================
# 3. CONTROLLO CONDIZIONI
# ============================================================

def require(condition, message):

    if not condition:
        raise ValueError(message)


# ============================================================
# 4. CHECKSUM SHA256
# ============================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as stream:

        for block in iter(
            lambda: stream.read(1024 * 1024),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# ============================================================
# 5. CONTROLLO PERCORSI
# ============================================================

def safe_run_file(path):

    p = Path(path)

    if not p.is_absolute():
        p = ROOT / p

    p = p.resolve()

    require(
        p.is_file()
        and p.is_relative_to(RUNS.resolve()),
        f"File fuori da runs o assente: {p}"
    )

    return p


# ============================================================
# 6. NORMALIZZAZIONE NOMI CROMOSOMICI
# ============================================================

def normalize(seqname):

    """
    Uniforma:
        chr1 -> 1
        chr6 -> 6
        chr10 -> 10

    Mantiene invariati:
        scaf_48
        scaf_77
        altri scaffold
    """

    name = seqname.strip()

    if (
        name.lower().startswith("chr")
        and name[3:].isdigit()
    ):

        return str(
            int(name[3:])
        )

    if name.isdigit():

        return str(
            int(name)
        )

    return name


# ============================================================
# 7. CARICAMENTO DELLA MAPPATURA STEP33
# ============================================================

def load_mapping():

    reports = sorted(
        RUNS.glob(
            "STEP33_LOCI_GENOMICI_*/"
            "RAPPORTO_STEP33.json"
        )
    )

    require(
        bool(reports),
        "RAPPORTO_STEP33.json non trovato"
    )

    report_path = safe_run_file(
        reports[-1]
    )

    report = json.loads(
        report_path.read_text(
            encoding="utf-8"
        )
    )

    require(
        str(report.get("step")) == "33",
        "Manifest STEP33 non valido"
    )

    require(
        report.get("status")
        ==
        "EXACT_GENOMIC_MAPPING_NOT_BIOLOGICAL_CLASSIFICATION",
        "Stato STEP33 inatteso"
    )

    # ========================================================
    # 8. VERIFICA DEL TSV ORIGINALE
    # ========================================================

    loci_file = safe_run_file(
        report["output_tsv"]
    )

    require(
        sha256(loci_file)
        == report["outputs_sha256"]["loci_tsv"],
        "Checksum STEP33 loci discordante"
    )

    # ========================================================
    # 9. CONTROLLO CONTIG E LUNGHEZZE
    # ========================================================

    contigs = report.get(
        "contig_lengths",
        []
    )

    reference_contigs = {

        item["name"]:
            int(item["length_bp"])

        for item in contigs
    }

    require(
        len(reference_contigs)
        == report.get("contigs_scanned"),
        "Elenco contig STEP33 incoerente"
    )

    require(
        sum(reference_contigs.values())
        == report.get("reference_bases_scanned"),
        "Dimensioni genoma STEP33 incoerenti"
    )

    # ========================================================
    # 10. LETTURA DEI LOCI GENOMICI
    # ========================================================

    rows = []

    with loci_file.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        reader = csv.DictReader(
            stream,
            delimiter="\t"
        )

        required = {

            "srna_id",
            "contig",
            "strand",
            "start_1based",
            "end_1based",
            "length_nt",
            "match_type",
            "reference"
        }

        require(
            required.issubset(
                reader.fieldnames or []
            ),
            "Colonne loci STEP33 incomplete"
        )

        seen = set()

        for row in reader:

            sid = row["srna_id"]
            contig = row["contig"]
            strand = row["strand"]

            require(
                sid in EXPECTED_LENGTHS
                and contig in reference_contigs
                and strand in ("+", "-"),
                "Locus sconosciuto"
            )

            start = int(
                row["start_1based"]
            )

            end = int(
                row["end_1based"]
            )

            length = int(
                row["length_nt"]
            )

            require(
                1 <= start <= end
                <= reference_contigs[contig]
                and end - start + 1
                == length
                == EXPECTED_LENGTHS[sid],
                "Coordinate locus non valide"
            )

            require(
                row["match_type"]
                == "exact_full_length",
                "Match non esatto nel registro"
            )

            identifier = (
                sid,
                contig,
                strand,
                start,
                end
            )

            require(
                identifier not in seen,
                f"Locus duplicato: {identifier}"
            )

            seen.add(identifier)

            # =================================================
            # 11. STRUTTURA DEL RECORD
            # =================================================

            rows.append({

                "srna_id":
                    sid,

                "contig":
                    contig,

                "canonical_contig":
                    normalize(contig),

                "strand":
                    strand,

                "start":
                    start,

                "end":
                    end,

                "length":
                    length,

                "gene_ids":
                    set(),

                "gene_biotypes":
                    set(),

                "gene_annotations":
                    0,

                "gene_same_strand":
                    0,

                "gene_opposite_strand":
                    0,

                "te_types":
                    Counter(),

                "te_feature_count":
                    0,

                "te_full_cover":
                    0,

                "other_repeat_types":
                    Counter(),

                "other_feature_types":
                    Counter()
            })

    require(
        bool(rows),
        "Nessun locus nel TSV STEP33"
    )

    # ========================================================
    # 12. VERIFICA NUMERO LOCI
    # ========================================================

    by_srna = Counter(
        row["srna_id"]
        for row in rows
    )

    for sid in EXPECTED_LENGTHS:

        require(
            by_srna[sid]
            ==
            report["by_srna"][sid]["total_exact_loci"],
            f"Conteggio STEP33 discordante: {sid}"
        )

    return (
        report_path,
        report,
        loci_file,
        rows
    )


# ============================================================
# 13. AUDIT DELLE ANNOTAZIONI GTF/GFF3
# ============================================================

def inspect_annotation(path, label):

    require(
        path.is_file(),
        f"File assente: {path}"
    )

    with path.open("rb") as stream:

        require(
            stream.read(2) == b"\x1f\x8b",
            f"File non gzip: {path}"
        )

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="strict"
    ) as stream:

        for line in stream:

            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue

            parts = (
                line.rstrip("\r\n")
                .split("\t")
            )

            if (
                len(parts) == 9
                and parts[3].isdigit()
                and parts[4].isdigit()
            ):

                print(

                    f"[RIFERIMENTO] {label}: "

                    f"seqid={parts[0]} "

                    f"-> {normalize(parts[0])}, "

                    f"feature={parts[2]}"
                )

                return

    raise ValueError(
        f"Nessun record valido in {path}"
    )


# ============================================================
# 14. INDICE DEI LOCI
# ============================================================

def locus_index(rows):

    mapping = defaultdict(list)

    for item in rows:

        mapping[
            item["canonical_contig"]
        ].append(item)

    indexed = {}

    for contig, loci in mapping.items():

        loci.sort(
            key=lambda item: item["start"]
        )

        starts = [
            item["start"]
            for item in loci
        ]

        max_length = max(
            item["length"]
            for item in loci
        )

        indexed[contig] = (
            starts,
            loci,
            max_length
        )

    return indexed


# ============================================================
# 15. CLASSIFICAZIONE DEGLI ELEMENTI ANNOTATI
# ============================================================

def label_te_feature(feature):

    if TE_RE.search(feature):

        return "TE"

    if OTHER_REPEAT_RE.search(feature):

        return "other_repeat"

    return "other_feature"


# ============================================================
# 16. ANALISI GTF / GFF3
# ============================================================

def scan_annotations(
    path,
    indexed,
    kind
):

    """
    Lettura streaming dei file compressi.

    Usa bisect per individuare soltanto i loci
    prossimi alla regione annotata, evitando
    di caricare tutti gli elementi GTF/GFF3
    in memoria.
    """

    counters = Counter()

    matched_contigs = set()
    annotation_contigs = set()

    with gzip.open(
        path,
        "rt",
        encoding="utf-8"
    ) as stream:

        for raw in stream:

            if (
                not raw.strip()
                or raw.startswith("#")
            ):
                continue

            columns = (
                raw.rstrip("\r\n")
                .split("\t")
            )

            if len(columns) != 9:

                counters[
                    "malformed_lines"
                ] += 1

                continue

            counters[
                "valid_records"
            ] += 1

            feature = columns[2]

            # =================================================
            # 17. PER GTF UTILIZZIAMO GLI INTERVALLI GENE
            # =================================================

            if (
                kind == "GTF"
                and feature != "gene"
            ):
                continue

            counters[
                "processed_features"
            ] += 1

            if (
                counters["processed_features"]
                % 500000 == 0
            ):

                print(

                    f"[LETTURA {kind}] "

                    f"annotazioni="
                    f"{counters['processed_features']:,} "

                    f"| intersezioni="
                    f"{counters['locus_feature_intersections']:,}",

                    flush=True
                )

            # =================================================
            # 18. CONTROLLO COORDINATE ANNOTAZIONE
            # =================================================

            try:

                start = int(
                    columns[3]
                )

                end = int(
                    columns[4]
                )

            except ValueError:

                counters[
                    "bad_coordinates"
                ] += 1

                continue

            if (
                start < 1
                or end < start
            ):

                counters[
                    "bad_coordinates"
                ] += 1

                continue

            chrom = normalize(
                columns[0]
            )

            annotation_contigs.add(
                chrom
            )

            if chrom not in indexed:
                continue

            # =================================================
            # 19. RICERCA INTERSEZIONI
            # =================================================

            (
                starts,
                loci,
                max_length
            ) = indexed[chrom]

            low = bisect.bisect_left(
                starts,
                max(
                    1,
                    start - max_length + 1
                )
            )

            high = bisect.bisect_right(
                starts,
                end
            )

            if low == high:
                continue

            classification = (
                label_te_feature(feature)
                if kind == "TE"
                else None
            )

            gene_id = (
                GENE_ID_RE.search(columns[8])
                if kind == "GTF"
                else None
            )

            biotype = (
                GENE_BIOTYPE_RE.search(columns[8])
                if kind == "GTF"
                else None
            )

            # =================================================
            # 20. ANNOTAZIONE DEI LOCI SOVRAPPOSTI
            # =================================================

            for locus in loci[low:high]:

                if (
                    locus["end"] < start
                    or locus["start"] > end
                ):
                    continue

                counters[
                    "locus_feature_intersections"
                ] += 1

                matched_contigs.add(
                    chrom
                )

                # =============================================
                # 21. INTERSEZIONE CON UN GENE
                # =============================================

                if kind == "GTF":

                    locus[
                        "gene_annotations"
                    ] += 1

                    if gene_id:

                        locus[
                            "gene_ids"
                        ].add(
                            gene_id.group(1)
                        )

                    if biotype:

                        locus[
                            "gene_biotypes"
                        ].add(
                            biotype.group(1)
                        )

                    if (
                        columns[6]
                        == locus["strand"]
                    ):

                        locus[
                            "gene_same_strand"
                        ] += 1

                    elif columns[6] in ("+", "-"):

                        locus[
                            "gene_opposite_strand"
                        ] += 1

                # =============================================
                # 22. INTERSEZIONE CON TE
                # =============================================

                elif classification == "TE":

                    locus[
                        "te_types"
                    ][feature] += 1

                    locus[
                        "te_feature_count"
                    ] += 1

                    if (
                        start <= locus["start"]
                        and end >= locus["end"]
                    ):

                        locus[
                            "te_full_cover"
                        ] += 1

                # =============================================
                # 23. ALTRE REGIONI RIPETUTE
                # =============================================

                elif classification == "other_repeat":

                    locus[
                        "other_repeat_types"
                    ][feature] += 1

                else:

                    locus[
                        "other_feature_types"
                    ][feature] += 1

    counters[
        "annotation_contigs"
    ] = len(annotation_contigs)

    counters[
        "matched_contigs"
    ] = len(matched_contigs)

    print(

        f"[ANNOTAZIONE {kind}] "

        f"features="
        f"{counters['processed_features']} "

        f"| intersezioni="
        f"{counters['locus_feature_intersections']} "

        f"| contig annotati="
        f"{counters['annotation_contigs']}"
    )

    return counters


# ============================================================
# 24. CLASSIFICAZIONE FINALE DEI LOCI
# ============================================================

def make_results(rows):

    detailed = []

    summary = {}

    chromosomes = {
        str(i)
        for i in range(1, 11)
    }

    ordered = sorted(

        rows,

        key=lambda item: (

            item["srna_id"],
            item["contig"],
            item["start"],
            item["strand"]
        )
    )

    for item in ordered:

        te = bool(
            item["te_types"]
        )

        gene = bool(
            item["gene_annotations"]
        )

        other = bool(
            item["other_repeat_types"]
        )

        # ====================================================
        # 25. CATEGORIE DI ANNOTAZIONE
        # ====================================================

        if te and gene:

            annotation = (
                "TE_and_gene_overlap"
            )

        elif te:

            annotation = (
                "TE_overlap"
            )

        elif gene:

            annotation = (
                "gene_overlap"
            )

        elif other:

            annotation = (
                "other_repeat_annotation"
            )

        elif item["other_feature_types"]:

            annotation = (
                "other_GFF_feature"
            )

        else:

            annotation = (
                "no_overlap_in_loaded_annotations"
            )

        # ====================================================
        # 26. REGISTRO DEL LOCUS
        # ====================================================

        detailed.append({

            "srna_id":
                item["srna_id"],

            "contig":
                item["contig"],

            "strand":
                item["strand"],

            "start_1based":
                item["start"],

            "end_1based":
                item["end"],

            "length_nt":
                item["length"],

            "contig_class":
                (
                    "chromosome"
                    if item["canonical_contig"]
                    in chromosomes
                    else "scaffold_or_other"
                ),

            "gene_overlap":
                "yes" if gene else "no",

            "gene_ids":
                ";".join(
                    sorted(item["gene_ids"])
                ),

            "gene_biotypes":
                ";".join(
                    sorted(item["gene_biotypes"])
                ),

            "gene_features_overlapping":
                item["gene_annotations"],

            "gene_same_strand":
                item["gene_same_strand"],

            "gene_opposite_strand":
                item["gene_opposite_strand"],

            "TE_feature_overlap":
                "yes" if te else "no",

            "TE_feature_annotations_overlapping":
                item["te_feature_count"],

            "TE_feature_annotations_fully_covering_site":
                item["te_full_cover"],

            "TE_types":
                ";".join(
                    sorted(item["te_types"])
                ),

            "other_repeat_types":
                ";".join(
                    sorted(
                        item["other_repeat_types"]
                    )
                ),

            "other_GFF_feature_types":
                ";".join(
                    sorted(
                        item["other_feature_types"]
                    )
                ),

            "annotation_category":
                annotation,

            "biological_status":
                "computational_origin_not_validated"
        })

    # ========================================================
    # 27. STATISTICHE PER SMALL RNA
    # ========================================================

    for sid in EXPECTED_LENGTHS:

        subset = [

            row
            for row in detailed

            if row["srna_id"] == sid
        ]

        categories = Counter(
            row["annotation_category"]
            for row in subset
        )

        summary[sid] = {

            "total_exact_matches":
                len(subset),

            "chromosome_matches":
                sum(
                    row["contig_class"]
                    == "chromosome"
                    for row in subset
                ),

            "scaffold_other_matches":
                sum(
                    row["contig_class"]
                    != "chromosome"
                    for row in subset
                ),

            "with_TE_feature":
                sum(
                    row["TE_feature_overlap"]
                    == "yes"
                    for row in subset
                ),

            "with_gene_feature":
                sum(
                    row["gene_overlap"]
                    == "yes"
                    for row in subset
                ),

            "with_other_repeat_feature":
                sum(
                    bool(row["other_repeat_types"])
                    for row in subset
                ),

            "with_TE_full_covering_feature":
                sum(

                    int(
                        row[
                            "TE_feature_annotations_fully_covering_site"
                        ]
                    ) > 0

                    for row in subset
                ),

            "annotation_categories":
                dict(categories),

            "contig_counts_top10":
                Counter(

                    row["contig"]
                    for row in subset

                ).most_common(10)
        }

    return detailed, summary


# ============================================================
# 28. ESPORTAZIONE TSV
# ============================================================

def write_tsv(path, rows):

    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as stream:

        writer = csv.DictWriter(
            stream,
            fieldnames=list(rows[0]),
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(rows)


# ============================================================
# 29. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    modes = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    modes.add_argument(
        "--audit",
        action="store_true"
    )

    modes.add_argument(
        "--run",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 70)

    print(
        "AI PRIMETTA V3.3 | STEP 34 "
        "- ANNOTAZIONE LOCI GENOMICI"
    )

    print("=" * 70)

    # ========================================================
    # 30. AUDIT STEP33
    # ========================================================

    (
        report_path,
        report,
        loci_tsv,
        rows
    ) = load_mapping()

    print(

        f"[PASS] STEP33 e SHA256: "
        f"{len(rows)} loci; "

        + ", ".join(

            f"{sid}="
            f"{report['by_srna'][sid]['total_exact_loci']}"

            for sid in EXPECTED_LENGTHS
        )
    )

    # ========================================================
    # 31. CONTROLLO GTF E GFF3
    # ========================================================

    inspect_annotation(
        GTF,
        "GTF"
    )

    inspect_annotation(
        TE_GFF,
        "TE GFF3"
    )

    print(
        "[PASS] Normalizzazione cromosomi: "
        "chr1 e 1 corrispondono; "
        "gli scaffold restano invariati"
    )

    # ========================================================
    # 32. MODALITA' AUDIT
    # ========================================================

    if args.audit:

        print(
            "[OK] Audit in sola lettura; "
            "annotazioni complete non ancora esaminate"
        )

        return

    # ========================================================
    # 33. COSTRUZIONE INDICE
    # ========================================================

    indexed = locus_index(
        rows
    )

    # ========================================================
    # 34. INTERSEZIONE CON GENI
    # ========================================================

    stats_gtf = scan_annotations(
        GTF,
        indexed,
        "GTF"
    )

    # ========================================================
    # 35. INTERSEZIONE CON ELEMENTI TRASPONIBILI
    # ========================================================

    stats_te = scan_annotations(
        TE_GFF,
        indexed,
        "TE"
    )

    # ========================================================
    # 36. CLASSIFICAZIONE DEI RISULTATI
    # ========================================================

    detailed, summary = make_results(
        rows
    )

    require(
        len(detailed) == len(rows),
        "Numero loci finale incongruente"
    )

    # ========================================================
    # 37. STAMPA RISULTATI SCIENTIFICI
    # ========================================================

    for sid in EXPECTED_LENGTHS:

        values = summary[sid]

        print(

            f"[RISULTATO] {sid}: "

            f"totale="
            f"{values['total_exact_matches']} "

            f"| TE="
            f"{values['with_TE_feature']} "

            f"| geni="
            f"{values['with_gene_feature']} "

            f"| cromosomi="
            f"{values['chromosome_matches']} "

            f"| scaffold="
            f"{values['scaffold_other_matches']}"
        )

        print(

            f"[CATEGORIE] {sid}: "

            + json.dumps(
                values["annotation_categories"],
                ensure_ascii=False,
                sort_keys=True
            )
        )

    # ========================================================
    # 38. CREAZIONE CARTELLA OUTPUT
    # ========================================================

    timestamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    outdir = RUNS / (
        "STEP34_ANNOTAZIONE_LOCI_"
        + timestamp
    )

    outdir.mkdir(
        parents=True,
        exist_ok=False
    )

    loci_output = (
        outdir
        / "ANNOTAZIONE_589_LOCI.tsv"
    )

    json_output = (
        outdir
        / "RAPPORTO_STEP34.json"
    )

    # ========================================================
    # 39. SCRITTURA DEL REGISTRO
    # ========================================================

    write_tsv(
        loci_output,
        detailed
    )

    # ========================================================
    # 40. RAPPORTO SCIENTIFICO
    # ========================================================

    payload = {

        "step":
            34,

        "status":
            "EXPLORATORY_GENOMIC_ANNOTATION_NOT_BIOLOGICAL_CLASSIFICATION",

        "source_step33":
            str(
                report_path.relative_to(ROOT)
            ),

        "source_step33_sha256":
            sha256(report_path),

        "source_loci_sha256":
            sha256(loci_tsv),

        "source_genome_sha256_from_step33":
            report["genome_compressed_sha256"],

        "gtf_file":
            str(GTF.relative_to(ROOT)),

        "gtf_sha256":
            sha256(GTF),

        "te_gff_file":
            str(TE_GFF.relative_to(ROOT)),

        "te_gff_sha256":
            sha256(TE_GFF),

        "gtf_scan":
            dict(stats_gtf),

        "te_scan":
            dict(stats_te),

        "total_loci":
            len(detailed),

        "by_srna":
            summary,

        "output_tsv":
            str(
                loci_output.relative_to(ROOT)
            ),

        "output_sha256":
            sha256(loci_output),

        "limitations": [

            "Classificazione TE conservativa "
            "basata sui tipi di feature GFF3.",

            "Le altre regioni ripetute sono "
            "classificate separatamente.",

            "Le intersezioni con geni utilizzano "
            "gli intervalli gene del GTF.",

            "Le corrispondenze multiple su scaffold "
            "possono riflettere ripetizioni "
            "o ridondanza nell'assembly.",

            "La localizzazione non dimostra "
            "la biogenesi o l'espressione "
            "della small RNA.",

            "L'assenza di sovrapposizione "
            "nelle annotazioni utilizzate "
            "non prova l'assenza biologica "
            "di geni o trasposoni.",

            "Nessun training, RNAplex "
            "o RNAplfold eseguito."
        ],

        "gui_modified":
            False,

        "docker_modified":
            False,

        "models_modified":
            False
    }

    json_output.write_text(

        json.dumps(
            payload,
            ensure_ascii=False,
            indent=2
        ) + "\n",

        encoding="utf-8"
    )

    # ========================================================
    # 41. OUTPUT FINALE
    # ========================================================

    print(
        "[TSV]",
        loci_output
    )

    print(
        "[REPORT]",
        json_output
    )

    print(
        "[OK] STEP34 completato; "
        "GUI, Docker, modelli e dati sorgenti invariati"
    )


# ============================================================
# 42. AVVIO PROGRAMMA
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(

            f"[ERRORE] "
            f"{type(exc).__name__}: {exc}",

            file=sys.stderr
        )

        sys.exit(1)