#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.8 - VERIFICA RISORSE GENOMICHE B73

Obiettivi:
1. Recuperare inventario STEP 17.7
2. Esaminare contenuto e intestazioni FASTA
3. Distinguere possibili FASTA genomici
4. Verificare annotazione GTF B73
5. Preparare il successivo step genomico

NON modifica:
- GUI V3.2
- motore V3.3
- modelli ML/GNN
- dataset originali

Nessun download, training o dato inventato.
"""

from __future__ import annotations

import csv
import gzip
import json
import re

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

GTF = (
    ROOT
    / "data/reference/zea_mays_ensembl_plants"
    / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.63.gtf.gz"
)

MAX_FASTA_HEADERS = 5
MAX_GTF_RECORDS = 10000


def open_text(path):

    if str(path).lower().endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8",
            errors="replace"
        )

    return path.open(
        "r",
        encoding="utf-8",
        errors="replace"
    )


def latest_step17_7():

    folders = sorted(
        [
            p for p in RUNS.glob(
                "STEP17_7_AUDIT_GENOMICO_*"
            )
            if p.is_dir()
        ],
        reverse=True
    )

    if not folders:
        raise RuntimeError(
            "Risultati STEP 17.7 non trovati"
        )

    path = folders[0] / "RAPPORTO_STEP17_7.json"

    if not path.is_file():
        raise FileNotFoundError(path)

    data = json.loads(
        path.read_text(encoding="utf-8")
    )

    return data, path


def inspect_fasta(path):

    headers = []
    sample_bases = 0
    ambiguous = Counter()

    current_length = 0
    sample_record_lengths = []

    with open_text(path) as handle:

        for line in handle:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if headers:
                    sample_record_lengths.append(
                        current_length
                    )

                headers.append(line[1:][:250])
                current_length = 0

                if len(headers) > MAX_FASTA_HEADERS:
                    headers.pop()
                    break

            else:

                current_length += len(line)

                sample_bases += len(line)

                if sample_bases <= 100000:

                    ambiguous.update(
                        base
                        for base in line.upper()
                        if base not in "ACGTUN"
                    )

        else:

            if headers:
                sample_record_lengths.append(
                    current_length
                )

    name = path.name.lower()

    header_text = " ".join(headers).lower()

    genomic_keywords = [
        "chromosome",
        "dna:chromosome",
        "scaffold",
        "supercontig",
        "genomic",
        "genome",
    ]

    transcript_keywords = [
        "cdna",
        "transcript",
        "mrna",
        "ncrna",
        "cds",
    ]

    genomic_signal = any(
        keyword in name or keyword in header_text
        for keyword in genomic_keywords
    )

    transcript_signal = any(
        keyword in name or keyword in header_text
        for keyword in transcript_keywords
    )

    if genomic_signal and not transcript_signal:
        classification = "POSSIBILE_GENOMICO"

    elif transcript_signal and not genomic_signal:
        classification = "POSSIBILE_TRASCRITTOMICO"

    elif genomic_signal and transcript_signal:
        classification = "AMBIGUO"

    else:
        classification = "DA_VERIFICARE"

    return {
        "file": str(path),
        "size_mb": round(
            path.stat().st_size / 1048576,
            3
        ),
        "classification": classification,
        "first_headers": headers,
        "sample_record_lengths": (
            sample_record_lengths
        ),
        "nonstandard_bases_sample": dict(
            ambiguous.most_common(10)
        ),
        "note": (
            "Classificazione euristica. "
            "Il file non viene validato come "
            "genoma completo."
        ),
    }


def inspect_gtf(path):

    if not path.is_file():

        return {
            "presente": False,
            "file": str(path),
        }

    chromosomes = Counter()
    features = Counter()
    biotypes = Counter()

    total = 0
    examples = []

    with open_text(path) as handle:

        for line in handle:

            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            total += 1

            chromosome = fields[0]
            feature = fields[2]
            attributes = fields[8]

            chromosomes[chromosome] += 1
            features[feature] += 1

            match = re.search(
                r'(?:gene_biotype|gene_type) "([^"]+)"',
                attributes
            )

            if match:
                biotypes[match.group(1)] += 1

            if len(examples) < 3:

                examples.append({
                    "chromosome": chromosome,
                    "feature": feature,
                    "start": fields[3],
                    "end": fields[4],
                    "strand": fields[6],
                    "attributes": attributes[:180],
                })

            if total >= MAX_GTF_RECORDS:
                break

    return {
        "presente": True,
        "file": str(path),
        "records_sampled": total,
        "sample_limit": MAX_GTF_RECORDS,
        "chromosomes_sample": dict(chromosomes),
        "features_sample": dict(features),
        "gene_biotypes_sample": dict(biotypes),
        "examples": examples,
        "analisi_completa": False,
    }


def write_csv(path, rows):

    if not rows:
        return

    with path.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0])
        )

        writer.writeheader()
        writer.writerows(rows)


def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 17.8")
    print("VERIFICA RISORSE GENOMICHE B73")
    print("=" * 72)

    data, source = latest_step17_7()

    resources = [
        r for r in data.get("risorse", [])
        if r["tipo"] != "ANNOTAZIONE"
    ]

    print("\n[1] INVENTARIO STEP 17.7")
    print("Fonte:", source)
    print("FASTA da controllare:", len(resources))

    results = []

    print("\n[2] ISPEZIONE FASTA")

    for index, resource in enumerate(
        resources,
        start=1
    ):

        path = Path(resource["file"])

        print(
            f"\n[{index}/{len(resources)}]",
            path.name
        )

        if not path.is_file():

            results.append({
                "file": str(path),
                "classification": "FILE_ASSENTE"
            })

            print("[ATTENZIONE] File assente")
            continue

        try:

            info = inspect_fasta(path)
            results.append(info)

            print(
                "Dimensione:",
                info["size_mb"],
                "MB"
            )

            print(
                "Classificazione:",
                info["classification"]
            )

            for header in info["first_headers"][:2]:
                print("HEADER:", header[:160])

        except Exception as error:

            results.append({
                "file": str(path),
                "classification": "ERRORE_LETTURA",
                "errore": str(error),
            })

            print("[ERRORE]", error)

    print("\n[3] CONTROLLO GTF")

    gtf_info = inspect_gtf(GTF)

    print("GTF presente:", gtf_info["presente"])

    if gtf_info["presente"]:

        print(
            "Record campionati:",
            gtf_info["records_sampled"]
        )

        print(
            "Cromosomi/contig nel campione:",
            len(gtf_info["chromosomes_sample"])
        )

        print(
            "Feature:",
            gtf_info["features_sample"]
        )

    print("\n[4] VALUTAZIONE DISPONIBILITA")

    genome_candidates = [
        item for item in results
        if item.get("classification")
        == "POSSIBILE_GENOMICO"
    ]

    unclassified = [
        item for item in results
        if item.get("classification")
        in ("DA_VERIFICARE", "AMBIGUO")
    ]

    if genome_candidates:

        print(
            "[INFO] Possibili FASTA genomici:",
            len(genome_candidates)
        )

        for item in genome_candidates:
            print(" -", item["file"])

        recommendation = (
            "Verificare assemblaggio, cromosomi, "
            "dimensioni e provenienza dei FASTA "
            "prima di usarli per il mapping genomico."
        )

    else:

        print(
            "[INFO] Nessun FASTA identificato "
            "come chiaramente genomico."
        )

        recommendation = (
            "Identificare e acquisire, se assente, "
            "il FASTA genomico Zea mays B73 "
            "NAM-5.0 coerente con il GTF."
        )

    print(
        "\nRisorse ancora ambigue:",
        len(unclassified)
    )

    print("\nRaccomandazione:")
    print(recommendation)

    print("\n[5] SALVATAGGIO")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_8_VERIFICA_GENOMA_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.8",
        "stato": "AUDIT_COMPLETATO",
        "fonte": str(source),
        "fasta_ispezionati": len(results),
        "risultati_fasta": results,
        "gtf": gtf_info,
        "genome_candidates": genome_candidates,
        "risorse_ambigue": unclassified,
        "raccomandazione": recommendation,
        "genoma_completo_verificato": False,
        "validazione_biologica": False,
        "motore_modificato": False,
        "gui_modificata": False,
    }

    json_path = (
        output / "RAPPORTO_STEP17_8.json"
    )

    json_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    csv_rows = [
        {
            "file": item["file"],
            "classification": item.get(
                "classification", ""
            ),
            "size_mb": item.get(
                "size_mb", ""
            ),
            "first_header": (
                item.get("first_headers", [""])[0]
                if item.get("first_headers")
                else ""
            )
        }
        for item in results
    ]

    csv_path = (
        output / "CLASSIFICAZIONE_FASTA.csv"
    )

    write_csv(csv_path, csv_rows)

    print("\n" + "=" * 72)
    print("STEP 17.8 COMPLETATO")
    print("=" * 72)

    print(
        "FASTA esaminati:",
        len(results)
    )

    print(
        "Possibili genomici:",
        len(genome_candidates)
    )

    print(
        "Risorse ambigue:",
        len(unclassified)
    )

    print("Report JSON:", json_path)
    print("Report CSV:", csv_path)

    print("\n[OK] GUI invariata")
    print("[OK] Motore invariato")
    print("[OK] Modelli ML/GNN invariati")
    print("[OK] Nessun download")
    print("[OK] Nessun training")
    print("[OK] File originali preservati")

    print(
        "\n[NOTA] Questo step verifica le "
        "risorse disponibili, non identifica "
        "ancora i loci genomici degli sRNA."
    )


if __name__ == "__main__":
    main()