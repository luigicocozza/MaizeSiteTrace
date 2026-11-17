#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 17.13
INTERSEZIONE LOCI SMALL RNA CON ANNOTAZIONI TE

Dati reali:
- PRJNA419326: 5 small RNA
- B73 NAM-5.0: 646 loci genomici esatti
- Annotazioni TE ufficiali MaizeGDB

Non modifica GUI V3.2, motore V3.3,
modelli, dati originali o risultati precedenti.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import urllib.parse
import urllib.request

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home() / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"
REF = BASE / "references"

ASSEMBLY = "Zm-B73-REFERENCE-NAM-5.0"
TE_NAME = ASSEMBLY + ".TE.gff3.gz"

TE_URL = (
    "https://download.maizegdb.org/"
    "Zm-B73-REFERENCE-NAM-5.0/"
    + TE_NAME
)

TE_FILE = REF / TE_NAME
TE_PART = REF / (TE_NAME + ".part")

EXPECTED_LOCI = 646
EXPECTED_RNA = 5


def section(number, name):
    print("\n" + "=" * 70)
    print(f"[{number}] {name}")
    print("=" * 70)


def short_id(name):
    for part in name.split("|"):
        if part.startswith("feature="):
            return part.split("=", 1)[1]
    return name.split("|")[0]


def normalize_contig(value):
    value = value.strip()

    match = re.fullmatch(
        r"(?:chr|chromosome_?)?([1-9]|10)",
        value,
        flags=re.IGNORECASE,
    )

    if match:
        return "chr" + match.group(1)

    return value


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(4 * 1024 * 1024),
            b""
        ):
            digest.update(block)
    return digest.hexdigest()


def latest_step11():
    folders = sorted(
        p for p in RUNS.glob(
            "STEP17_11_ANNOTAZIONE_GTF_*"
        )
        if p.is_dir()
        and (p / "LOCI_ANNOTATI_STEP17_11.csv").is_file()
    )

    if not folders:
        raise RuntimeError(
            "STEP 17.11 non trovato"
        )

    return folders[-1]


def load_loci(path):

    with path.open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        loci = list(csv.DictReader(handle))

    if len(loci) != EXPECTED_LOCI:
        raise RuntimeError(
            f"Attesi 646 loci, trovati {len(loci)}"
        )

    names = {r["rna_id"] for r in loci}
    if len(names) != EXPECTED_RNA:
        raise RuntimeError(
            "Numero small RNA inatteso"
        )

    seen = set()

    for row in loci:

        row["start"] = int(row["start_1based"])
        row["end"] = int(row["end_1based"])
        row["contig"] = normalize_contig(
            row["chromosome"]
        )
        row["srna"] = short_id(row["rna_id"])

        if (
            row["start"] < 1
            or row["end"] < row["start"]
        ):
            raise RuntimeError(
                "Coordinate locus non valide"
            )

        key = (
            row["rna_id"],
            row["contig"],
            row["start"],
            row["end"],
            row["strand"],
        )

        if key in seen:
            raise RuntimeError(
                f"Locus duplicato: {key}"
            )

        seen.add(key)

    return loci


def verify_gzip(path):
    """
    Legge integralmente il GFF compresso.
    Il controllo di decompressione non equivale
    a un checksum ufficiale di provenienza.
    """

    rows = 0
    headers = []

    with gzip.open(
        path, "rt", encoding="utf-8"
    ) as handle:

        for line in handle:

            if line.startswith("#"):
                if len(headers) < 15:
                    headers.append(line.strip())
                continue

            if line.strip():
                rows += 1

    if rows == 0:
        raise RuntimeError(
            "File GFF3 vuoto"
        )

    return rows, headers


def acquire_te():

    REF.mkdir(parents=True, exist_ok=True)

    if TE_FILE.is_file():

        print("[INFO] Annotazione TE gia presente")
        count, headers = verify_gzip(TE_FILE)

        return count, headers, "RIUTILIZZATO"

    print("[INFO] Download ufficiale:", TE_URL)

    # Download in file temporaneo.
    # Se l'interruzione lascia un .part, al prossimo
    # tentativo viene riscaricato integralmente.
    request = urllib.request.Request(
        TE_URL,
        headers={
            "User-Agent": "AI-PRIMETTA-V3.3-Research"
        },
    )

    try:
        with urllib.request.urlopen(
            request, timeout=120
        ) as response:

            if response.status != 200:
                raise RuntimeError(
                    "Risposta HTTP inattesa"
                )

            with TE_PART.open("wb") as output:
                shutil.copyfileobj(
                    response, output,
                    length=1024 * 1024
                )

        if TE_PART.stat().st_size < 1_000_000:
            raise RuntimeError(
                "Download TE sospettosamente piccolo"
            )

        count, headers = verify_gzip(TE_PART)

        if TE_FILE.exists():
            raise RuntimeError(
                "Destinazione TE gia presente: "
                "file non sovrascritto"
            )

        TE_PART.rename(TE_FILE)

        return count, headers, "SCARICATO"

    except Exception:
        # Conserva eventuale file parziale
        # esclusivamente per diagnostica.
        raise


def parse_gff_attrs(raw):

    attrs = {}

    for item in raw.strip().split(";"):

        if not item:
            continue

        if "=" not in item:
            continue

        key, value = item.split("=", 1)

        attrs[key] = urllib.parse.unquote(
            value.strip()
        )

    return attrs


def load_te_target_contigs(path, target_contigs):
    """
    Carica solo le annotazioni situate nei contig
    nei quali sono presenti i 646 loci.
    Evita di memorizzare l'intero GFF in RAM.

    Mantiene tutti i tipi di feature; la
    molteplicità gerarchica viene dichiarata.
    """

    annotations = defaultdict(list)
    feature_counts = Counter()
    all_contigs = set()

    line_count = 0
    invalid = 0
    records_selected = 0

    with gzip.open(
        path, "rt", encoding="utf-8"
    ) as handle:

        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            line_count += 1

            parts = line.rstrip("\n").split("\t")

            if len(parts) != 9:
                invalid += 1
                continue

            contig = normalize_contig(parts[0])
            all_contigs.add(contig)

            feature_counts[parts[2]] += 1

            if contig not in target_contigs:
                continue

            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                invalid += 1
                continue

            if start < 1 or end < start:
                invalid += 1
                continue

            attrs = parse_gff_attrs(parts[8])

            annotations[contig].append({
                "start": start,
                "end": end,
                "source": parts[1],
                "feature": parts[2],
                "strand": parts[6],
                "id": attrs.get("ID", ""),
                "parent": attrs.get("Parent", ""),
                "name": attrs.get("Name", ""),
                "classification": (
                    attrs.get("Classification")
                    or attrs.get("classification")
                    or attrs.get("Class")
                    or attrs.get("class")
                    or attrs.get("family")
                    or attrs.get("Family")
                    or attrs.get("type")
                    or ""
                ),
                "attributes": attrs,
            })

            records_selected += 1

    return annotations, {
        "record_totali": line_count,
        "record_selezionati": records_selected,
        "record_malformati": invalid,
        "contig_annotati": len(all_contigs),
        "contig_names": sorted(all_contigs),
        "feature_types": dict(feature_counts),
    }


def build_index(annotations):

    indexed = {}

    for contig, rows in annotations.items():

        rows.sort(
            key=lambda row: (
                row["start"], row["end"]
            )
        )

        starts = []
        prefix_max_end = []
        current = 0

        for row in rows:

            starts.append(row["start"])
            current = max(current, row["end"])
            prefix_max_end.append(current)

        indexed[contig] = (
            rows, starts, prefix_max_end
        )

    return indexed


def overlap_query(index_entry, start, end):

    if index_entry is None:
        return []

    rows, starts, prefix = index_entry

    position = bisect.bisect_right(
        starts, end
    ) - 1

    hits = []

    while position >= 0:

        if prefix[position] < start:
            break

        row = rows[position]

        if row["end"] >= start:
            hits.append(row)

        position -= 1

    return hits


def write_csv(path, rows, columns):

    with path.open(
        "w", encoding="utf-8", newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.13")
    print("ANNOTAZIONE ELEMENTI TRASPONIBILI B73")
    print("=" * 70)

    section(1, "LETTURA DEI 646 LOCI REALI")

    step11 = latest_step11()

    loci_file = (
        step11 / "LOCI_ANNOTATI_STEP17_11.csv"
    )

    loci = load_loci(loci_file)

    target_contigs = {
        row["contig"] for row in loci
    }

    print("[OK] Loci:", len(loci))
    print("[OK] Contig con loci:", len(target_contigs))

    section(2, "ANNOTAZIONE UFFICIALE MaizeGDB")

    expected_lines, header, source_status = acquire_te()

    print("[OK] File:", TE_FILE)
    print("[OK] Origine:", TE_URL)
    print("[OK] GFF letto:", expected_lines, "record")
    print("[OK] Stato:", source_status)

    section(3, "INDICIZZAZIONE DELLE ANNOTAZIONI TE")

    annotations, gff_stats = load_te_target_contigs(
        TE_FILE, target_contigs
    )

    index = build_index(annotations)

    matched_contigs = (
        target_contigs.intersection(
            gff_stats["contig_names"]
        )
    )

    print(
        "[OK] Record sui contig di interesse:",
        gff_stats["record_selezionati"]
    )

    print(
        "[OK] Contig riconciliati:",
        len(matched_contigs),
        "/", len(target_contigs)
    )

    print(
        "[INFO] Tipi di feature GFF:",
        gff_stats["feature_types"]
    )

    # Evita di interpretare zero sovrapposizioni
    # quando il file usa identificativi incompatibili.
    if not matched_contigs:
        raise RuntimeError(
            "Nessun contig corrispondente tra "
            "loci genomici e annotazione TE. "
            "Verificare nomenclatura e assembly."
        )

    section(4, "INTERSEZIONE LOCI - TE")

    summary = []
    associations = []

    hit_counter = Counter()
    per_rna = defaultdict(Counter)
    family_counter = Counter()
    feature_counter = Counter()

    for locus in loci:

        entry = index.get(locus["contig"])

        overlapping = overlap_query(
            entry,
            locus["start"],
            locus["end"]
        )

        # Due righe GFF diverse possono rappresentare
        # componenti dello stesso elemento TE.
        # Pertanto distinguiamo loci e righe sovrapposte.
        key = (
            locus["rna_id"],
            locus["chromosome"],
            locus["start"],
            locus["end"],
            locus["strand"],
        )

        if overlapping:
            hit_counter["loci_con_overlap"] += 1
            per_rna[locus["srna"]]["con_TE"] += 1
        else:
            hit_counter["loci_senza_overlap"] += 1
            per_rna[locus["srna"]]["senza_TE"] += 1

        if entry is None:
            hit_counter["loci_contig_senza_TE_GFF"] += 1

        per_rna[locus["srna"]]["totali"] += 1

        names = set()
        classifications = set()
        features = set()

        for te in overlapping:

            overlap_nt = (
                min(locus["end"], te["end"])
                - max(locus["start"], te["start"])
                + 1
            )

            if overlap_nt <= 0:
                raise RuntimeError(
                    "Errore nel calcolo intersezione"
                )

            feature_counter[te["feature"]] += 1

            family = (
                te["classification"]
                or "NON_SPECIFICATA"
            )

            family_counter[family] += 1

            features.add(te["feature"])
            classifications.add(family)

            if te["id"]:
                names.add(te["id"])

            associations.append({
                "rna_id": locus["srna"],
                "chromosome": locus["chromosome"],
                "locus_start": locus["start"],
                "locus_end": locus["end"],
                "locus_strand": locus["strand"],
                "te_start": te["start"],
                "te_end": te["end"],
                "te_strand": te["strand"],
                "te_feature": te["feature"],
                "te_id": te["id"],
                "te_parent": te["parent"],
                "te_name": te["name"],
                "te_classification": family,
                "overlap_nt": overlap_nt,
            })

        summary.append({
            "rna_id": locus["srna"],
            "chromosome": locus["chromosome"],
            "start_1based": locus["start"],
            "end_1based": locus["end"],
            "strand": locus["strand"],
            "genomic_region": locus.get("region", ""),
            "te_overlap": bool(overlapping),
            "te_overlap_records": len(overlapping),
            "te_ids": ";".join(sorted(names)),
            "te_classes": ";".join(
                sorted(classifications)
            ),
            "te_feature_types": ";".join(
                sorted(features)
            ),
            "contig_in_te_annotation": (
                locus["contig"] in matched_contigs
            ),
        })

    if len(summary) != EXPECTED_LOCI:
        raise RuntimeError(
            "Sono stati persi loci durante l'analisi"
        )

    section(5, "RISULTATI PER SMALL RNA")

    overview = []

    for name in sorted(per_rna):

        counts = per_rna[name]

        total = counts["totali"]
        positive = counts["con_TE"]
        negative = counts["senza_TE"]

        if positive + negative != total:
            raise RuntimeError(
                "Conteggio TE incoerente: " + name
            )

        percentage = round(
            100 * positive / total, 3
        )

        overview.append({
            "srna": name,
            "loci_totali": total,
            "loci_con_sovrapposizione_TE": positive,
            "loci_senza_sovrapposizione_TE": negative,
            "percentuale_overlap_TE": percentage,
        })

        print()
        print("RNA:", name)
        print("Loci totali:", total)
        print("Loci con TE:", positive)
        print("Loci senza TE:", negative)
        print("Percentuale:", percentage, "%")

    print()
    print(
        "Loci con almeno una feature TE:",
        hit_counter["loci_con_overlap"]
    )

    print(
        "Loci senza feature TE:",
        hit_counter["loci_senza_overlap"]
    )

    print(
        "Associazioni locus-feature TE:",
        len(associations)
    )

    print(
        "Loci su contig privi di record TE:",
        hit_counter["loci_contig_senza_TE_GFF"]
    )

    print("\nFeature TE più frequenti:")

    for feature, count in feature_counter.most_common(10):
        print(" ", feature, ":", count)

    print("\nClassificazioni TE disponibili:")

    for family, count in family_counter.most_common(10):
        print(" ", family, ":", count)

    section(6, "GENERAZIONE REPORT")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_13_ANNOTAZIONE_TE_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    write_csv(
        output / "LOCI_CON_ANNOTAZIONI_TE.csv",
        summary,
        list(summary[0].keys()),
    )

    write_csv(
        output / "ASSOCIAZIONI_LOCI_TE.csv",
        associations,
        [
            "rna_id",
            "chromosome",
            "locus_start",
            "locus_end",
            "locus_strand",
            "te_start",
            "te_end",
            "te_strand",
            "te_feature",
            "te_id",
            "te_parent",
            "te_name",
            "te_classification",
            "overlap_nt",
        ],
    )

    write_csv(
        output / "SINTESI_TE_PER_SRNA.csv",
        overview,
        list(overview[0].keys()),
    )

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "17.13",
        "status": "TE_OVERLAP_COMPLETED",
        "reference": ASSEMBLY,
        "source_step17_11": str(loci_file),
        "te_url": TE_URL,
        "te_file": str(TE_FILE),
        "te_sha256_local": sha256(TE_FILE),
        "te_integrity_gzip": True,
        "te_source_checksum_compared": False,
        "te_download_status": source_status,
        "gtf_contig_matches": len(matched_contigs),
        "target_contigs": len(target_contigs),
        "gff_statistics": {
            "record_totali": gff_stats["record_totali"],
            "record_selezionati": gff_stats["record_selezionati"],
            "record_malformati": gff_stats["record_malformati"],
            "contig_annotati": gff_stats["contig_annotati"],
            "feature_types": gff_stats["feature_types"],
        },
        "loci_totali": len(loci),
        "loci_con_overlap_TE": hit_counter[
            "loci_con_overlap"
        ],
        "loci_senza_overlap_TE": hit_counter[
            "loci_senza_overlap"
        ],
        "associazioni_feature_TE": len(associations),
        "summary_per_srna": overview,
        "te_feature_counts": dict(feature_counter),
        "te_classification_counts": dict(family_counter),
        "limitations": [
            "Una sovrapposizione con TE annotato "
            "non dimostra biogenesi o regolazione.",
            "Alcuni record GFF possono essere "
            "gerarchici o riferirsi a sub-feature.",
            "Le associazioni locus-feature non "
            "equivalgono al numero di TE distinti.",
            "L'assenza di overlap non esclude "
            "regioni ripetitive non annotate.",
            "Il checksum ufficiale del file TE "
            "non e stato confrontato.",
            "I loci senza annotazioni TE nel loro "
            "contig richiedono interpretazione cauta.",
        ],
        "gui_modified": False,
        "engine_modified": False,
        "training_performed": False,
        "biological_validation": False,
    }

    report_path = (
        output / "RAPPORTO_STEP17_13.json"
    )

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
    print("AI PRIMETTA V3.3 - STEP 17.13 TERMINATO")
    print("=" * 70)

    print("Loci analizzati:", len(loci))
    print(
        "Loci con overlap TE:",
        hit_counter["loci_con_overlap"]
    )
    print(
        "Loci senza overlap TE:",
        hit_counter["loci_senza_overlap"]
    )

    print("\nRisultati:", output)
    print("Report:", report_path)

    print("\n[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli invariati")
    print("[OK] Nessun training")
    print("[OK] Dati originali preservati")

    print(
        "\n[NOTA] Le sovrapposizioni con "
        "elementi trasponibili sono "
        "annotazioni genomiche computazionali."
    )


if __name__ == "__main__":
    main()