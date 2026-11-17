#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 25

Analisi delle energie RNAplex e della copertura
dei siti candidati individuati nello STEP22.

MODALITA':
    --audit
        Controlla gli input, analizza i dati
        e stampa le statistiche senza scrivere.

    --run
        Esegue gli stessi controlli e crea
        TSV, JSON e rapporto Markdown.

Non modifica GUI, Docker, modelli o dati originali.
Non riesegue RNAplex.
"""

import argparse
import csv
import hashlib
import json
import math
import random
import statistics
import sys

from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal
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

SRNAS = {
    "sRNA_00000435": 24,
    "sRNA_00000037": 25
}


# ============================================================
# 2. VALIDAZIONE DEI PERCORSI
# ============================================================

def source(path):

    p = Path(path)

    if not p.is_absolute():
        p = ROOT / p

    p = p.resolve()

    if (
        not p.is_file()
        or not p.is_relative_to(RUNS.resolve())
    ):

        raise ValueError(
            f"File assente o fuori runs: {p}"
        )

    return p


# ============================================================
# 3. CALCOLO CHECKSUM SHA256
# ============================================================

def sha(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for part in iter(
            lambda: f.read(1 << 20),
            b""
        ):

            h.update(part)

    return h.hexdigest()


# ============================================================
# 4. LETTURA TSV
# ============================================================

def tsv(path):

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        names = reader.fieldnames or []

        if (
            not names
            or len(set(names)) != len(names)
        ):

            raise ValueError(
                f"Intestazioni TSV non valide: {path}"
            )

        return list(reader)


# ============================================================
# 5. LETTURA FASTA
# ============================================================

def fasta(path):

    records = {}

    name = None
    parts = []

    with path.open(
        encoding="utf-8-sig"
    ) as f:

        for raw in f:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    if (
                        not parts
                        or name in records
                    ):

                        raise ValueError(
                            "Record FASTA duplicato o vuoto"
                        )

                    records[name] = (
                        "".join(parts)
                        .upper()
                        .replace("T", "U")
                    )

                name = line[1:].split()[0]
                parts = []

            elif name is None:

                raise ValueError(
                    "FASTA senza intestazione"
                )

            else:

                parts.append(line)

    if name is not None:

        if (
            not parts
            or name in records
        ):

            raise ValueError(
                "Record FASTA finale duplicato/vuoto"
            )

        records[name] = (
            "".join(parts)
            .upper()
            .replace("T", "U")
        )

    if not records:

        raise ValueError(
            "FASTA vuoto"
        )

    return records


# ============================================================
# 6. RECUPERO E VERIFICA DEI RISULTATI
# ============================================================

def check():

    paths = sorted(
        RUNS.glob(
            "STEP24_CONFRONTO_RNAPLEX_CORRETTO_*/"
            "RAPPORTO_STEP24.json"
        )
    )

    if not paths:

        raise FileNotFoundError(
            "Rapporto STEP24 corretto non trovato"
        )

    report_path = source(
        paths[-1]
    )

    report = json.loads(
        report_path.read_text(
            encoding="utf-8"
        )
    )

    if (
        report.get("step") != "24"
        or report.get("RNAplex_pairs") != 74
    ):

        raise ValueError(
            "Rapporto STEP24 non coerente"
        )

    step23_path = source(
        report["source_step23b"]
    )

    step23 = json.loads(
        step23_path.read_text(
            encoding="utf-8"
        )
    )

    manifest_path = source(
        step23["source_manifest"]
    )

    inputs = {

        "srna": source(
            manifest_path.parent
            / "SRNA_STANDARD_MAX25.fa"
        ),

        "targets": source(
            manifest_path.parent
            / "TRASCRITTI_B73_CANDIDATI.fa"
        ),

        "sites": source(
            manifest_path.parent
            / "SITI_E_RNASEQ_SNAPSHOT.tsv"
        ),

        "raw": source(
            step23_path.parent
            / "OUTPUT_GREZZO_RNAPLEX.txt"
        ),

        "table": source(
            step23_path.parent
            / "INTERAZIONI_RNAPLEX.tsv"
        )
    }

    # ========================================================
    # 7. CONTROLLO INTEGRITA'
    # ========================================================

    for key, path in inputs.items():

        if sha(path) != (
            report["source_sha256"][key]
        ):

            raise ValueError(
                f"Checksum modificato: {key}"
            )

    sequences = fasta(
        inputs["srna"]
    )

    transcripts = fasta(
        inputs["targets"]
    )

    if (
        {
            k: len(v)
            for k, v in sequences.items()
        } != SRNAS
        or len(transcripts) != 37
    ):

        raise ValueError(
            "Numero small RNA/trascritti inatteso"
        )

    # ========================================================
    # 8. CONTROLLO DEI SITI ORIGINALI STEP22
    # ========================================================

    sites = tsv(
        inputs["sites"]
    )

    if len(sites) != 37:

        raise ValueError(
            "Attesi 37 siti STEP22"
        )

    positives = {}

    mapping = {}

    for site in sites:

        sid = site["srna_id"]
        tx = site["transcript_id"]

        if (
            tx not in transcripts
            or tx in mapping
        ):

            raise ValueError(
                "Trascritto STEP22 assente o ripetuto"
            )

        a = int(
            site["start_1based"]
        )

        b = int(
            site["end_1based"]
        )

        if not (
            1 <= a <= b <= len(transcripts[tx])
        ):

            raise ValueError(
                "Coordinate STEP22 non valide"
            )

        observed = transcripts[tx][
            a - 1:b
        ]

        expected = (
            site["full_site_DNA_5to3"]
            .upper()
            .replace("T", "U")
        )

        if observed != expected:

            raise ValueError(
                "Sito STEP22 discordante dal FASTA"
            )

        mapping[tx] = site["gene_id"]

        if sid in sequences:

            positives[
                (sid, tx)
            ] = site

    if (
        len(positives) != 36
        or len(mapping) != 37
    ):

        raise ValueError(
            "Attesi 36 siti delle due sRNA "
            "su 37 trascritti"
        )

    # ========================================================
    # 9. LETTURA DEL CONFRONTO STEP24
    # ========================================================

    result_path = source(
        report["output_tsv"]
    )

    data = tsv(
        result_path
    )

    if len(data) != 74:

        raise ValueError(
            "Attese 74 righe STEP24"
        )

    required = {

        "srna_id",
        "transcript_id",
        "energy_kcal_mol",
        "target_start_1based",
        "target_end_1based",
        "srna_start_1based",
        "srna_end_1based",
        "target_fragment_RNA_5to3",
        "comparison",
        "overlap_nt"
    }

    if not required.issubset(
        set(data[0])
    ):

        raise ValueError(
            "Colonne STEP24 mancanti"
        )

    pairs = {}

    # ========================================================
    # 10. VALIDAZIONE DELLE 74 COPPIE
    # ========================================================

    for row in data:

        sid = row["srna_id"]
        tx = row["transcript_id"]

        key = (sid, tx)

        if (
            sid not in sequences
            or tx not in transcripts
            or key in pairs
        ):

            raise ValueError(
                f"Coppia imprevista/duplicata: {key}"
            )

        t1 = int(
            row["target_start_1based"]
        )

        t2 = int(
            row["target_end_1based"]
        )

        q1 = int(
            row["srna_start_1based"]
        )

        q2 = int(
            row["srna_end_1based"]
        )

        if not (
            1 <= t1 <= t2 <= len(transcripts[tx])
            and
            1 <= q1 <= q2 <= len(sequences[sid])
        ):

            raise ValueError(
                f"Coordinate STEP24 non valide: {key}"
            )

        observed = transcripts[tx][
            t1 - 1:t2
        ]

        expected = (
            row["target_fragment_RNA_5to3"]
            .upper()
            .replace("T", "U")
        )

        if observed != expected:

            raise ValueError(
                "Frammento target non coerente "
                "con il FASTA"
            )

        energy = Decimal(
            row["energy_kcal_mol"]
        )

        if not energy.is_finite():

            raise ValueError(
                "Energia non finita"
            )

        site = positives.get(
            key
        )

        overlap = 0
        coverage = None
        iou = None

        # ====================================================
        # 11. CALCOLO COPERTURA DEL SITO
        # ====================================================

        if site is not None:

            a = int(
                site["start_1based"]
            )

            b = int(
                site["end_1based"]
            )

            overlap = max(
                0,
                min(t2, b)
                - max(t1, a)
                + 1
            )

            length_site = (
                b - a + 1
            )

            length_rnaplex = (
                t2 - t1 + 1
            )

            coverage = (
                overlap / length_site
            )

            union = (
                length_site
                + length_rnaplex
                - overlap
            )

            iou = (
                overlap / union
            )

        if overlap != int(
            row["overlap_nt"]
        ):

            raise ValueError(
                f"Overlap errato: {key}"
            )

        expected_class = (
            "sovrapposizione"
            if site
            else "nessun_sito_STEP22_per_coppia"
        )

        if (
            row["comparison"] != expected_class
            or (
                site is not None
                and overlap == 0
            )
        ):

            raise ValueError(
                f"Classificazione imprevista: {key}"
            )

        pairs[key] = {

            "row": row,
            "energy": energy,
            "overlap": overlap,
            "coverage": coverage,
            "iou": iou,
            "site": site
        }

    expected_pairs = {

        (sid, tx)
        for sid in sequences
        for tx in transcripts
    }

    if set(pairs) != expected_pairs:

        raise ValueError(
            "Mancano coppie RNAplex"
        )

    counts = Counter(
        item["row"]["comparison"]
        for item in pairs.values()
    )

    if counts != {

        "sovrapposizione": 36,
        "nessun_sito_STEP22_per_coppia": 38
    }:

        raise ValueError(
            "Confronto STEP24 diverso dal log"
        )

    return (
        report_path,
        inputs,
        pairs,
        positives
    )


# ============================================================
# 12. MEDIANA
# ============================================================

def median(values):

    return (
        statistics.median(values)
        if values
        else None
    )


# ============================================================
# 13. CONFRONTO APPAIATO DELLE ENERGIE
# ============================================================

def summaries(pairs, positives):

    paired = []

    for (
        sid,
        tx
    ), site in sorted(
        positives.items()
    ):

        other = next(
            item
            for item in SRNAS
            if item != sid
        )

        matched = pairs[
            (sid, tx)
        ]

        control = pairs[
            (other, tx)
        ]

        # Energia più negativa:
        # interazione più favorevole
        # all'interno del modello RNAplex.

        delta = (
            matched["energy"]
            - control["energy"]
        )

        paired.append({

            "srna_id": sid,

            "transcript_id": tx,

            "gene_id": site["gene_id"],

            "matched_energy_kcal_mol":
                str(matched["energy"]),

            "other_srna_id":
                other,

            "other_energy_kcal_mol":
                str(control["energy"]),

            "delta_matched_minus_other_kcal_mol":
                str(delta),

            "matched_energy_more_negative":
                delta < 0,

            "overlap_nt":
                matched["overlap"],

            "site_coverage_fraction":
                f"{matched['coverage']:.6f}",

            "intersection_over_union":
                f"{matched['iou']:.6f}",

            "STEP22_full_score":
                site["full_score"]
        })

    # ========================================================
    # 14. DIFFERENZE ENERGETICHE
    # ========================================================

    deltas = [

        float(
            row[
                "delta_matched_minus_other_kcal_mol"
            ]
        )

        for row in paired
    ]

    negative = sum(
        x < 0
        for x in deltas
    )

    positive = sum(
        x > 0
        for x in deltas
    )

    ties = (
        len(deltas)
        - negative
        - positive
    )

    # ========================================================
    # 15. SIGN TEST ESPLORATIVO
    # ========================================================

    n = negative + positive

    if n:

        tail = sum(
            math.comb(n, k)
            for k in range(
                min(negative, positive) + 1
            )
        )

        p_value = min(
            1.0,
            2 * tail / (2 ** n)
        )

    else:

        p_value = None

    # ========================================================
    # 16. BOOTSTRAP DELLA MEDIANA
    # ========================================================

    rng = random.Random(
        20261023
    )

    samples = sorted(

        median([

            rng.choice(deltas)

            for _ in deltas
        ])

        for _ in range(5000)
    )

    # ========================================================
    # 17. COPERTURA DEI 36 SITI
    # ========================================================

    fractions = [

        pairs[key]["coverage"]

        for key in positives
    ]

    unions = [

        pairs[key]["iou"]

        for key in positives
    ]

    # ========================================================
    # 18. STATISTICHE PER SMALL RNA
    # ========================================================

    per_srna = {}

    for sid in SRNAS:

        subset = [

            row
            for row in paired
            if row["srna_id"] == sid
        ]

        ds = [

            float(
                row[
                    "delta_matched_minus_other_kcal_mol"
                ]
            )

            for row in subset
        ]

        per_srna[sid] = {

            "candidate_transcripts":
                len(subset),

            "favorable_matched_energy":
                sum(x < 0 for x in ds),

            "median_delta_kcal_mol":
                median(ds),

            "median_site_coverage":
                median([

                    float(
                        row[
                            "site_coverage_fraction"
                        ]
                    )

                    for row in subset
                ])
        }

    # ========================================================
    # 19. RAPPORTO STATISTICO
    # ========================================================

    summary = {

        "compared_transcripts":
            len(paired),

        "candidate_more_negative":
            negative,

        "other_more_negative":
            positive,

        "energy_ties":
            ties,

        "median_paired_delta_kcal_mol":
            median(deltas),

        "mean_paired_delta_kcal_mol":
            statistics.mean(deltas),

        "exploratory_two_sided_sign_test_p":
            p_value,

        "exploratory_median_delta_bootstrap_95_percentile_CI": [

            samples[
                int(0.025 * len(samples))
            ],

            samples[
                int(0.975 * len(samples)) - 1
            ]
        ],

        "median_site_coverage":
            median(fractions),

        "median_interval_iou":
            median(unions),

        "sites_covering_at_least_25pct":
            sum(
                x >= 0.25
                for x in fractions
            ),

        "sites_covering_at_least_50pct":
            sum(
                x >= 0.50
                for x in fractions
            ),

        "sites_covering_at_least_75pct":
            sum(
                x >= 0.75
                for x in fractions
            ),

        "sites_covering_at_least_90pct":
            sum(
                x >= 0.90
                for x in fractions
            ),

        "sites_completely_covered":
            sum(
                x == 1
                for x in fractions
            ),

        "by_srna":
            per_srna
    }

    return paired, summary


# ============================================================
# 20. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = parser.add_mutually_exclusive_group(
        required=True
    )

    mode.add_argument(
        "--audit",
        action="store_true"
    )

    mode.add_argument(
        "--run",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 67)

    print(
        "AI PRIMETTA V3.3 | STEP 25"
    )

    print(
        "COPERTURA ED ENERGIE RNAplex"
    )

    print("=" * 67)

    (
        report_path,
        inputs,
        pairs,
        positives
    ) = check()

    rows, stat = summaries(
        pairs,
        positives
    )

    # ========================================================
    # 21. OUTPUT DI CONTROLLO
    # ========================================================

    print(
        "[PASS] STEP24, FASTA, coordinate "
        "e checksum: 74/74"
    )

    print(

        "[COPERTURA] "
        f">=50%="
        f"{stat['sites_covering_at_least_50pct']}/36; "

        f">=75%="
        f"{stat['sites_covering_at_least_75pct']}/36; "

        f"completa="
        f"{stat['sites_completely_covered']}/36"
    )

    print(

        "[ENERGIA ABBINATA] "
        f"più favorevole sito preselezionato="
        f"{stat['candidate_more_negative']}; "

        f"altra sRNA="
        f"{stat['other_more_negative']}; "

        f"pareggi="
        f"{stat['energy_ties']}"
    )

    print(

        "[DELTA MEDIANO] "
        f"{stat['median_paired_delta_kcal_mol']} "
        "kcal/mol (match - altra sRNA)"
    )

    print(

        "[TEST ESPLORATIVO] "
        f"p={stat['exploratory_two_sided_sign_test_p']}"
    )

    for sid, value in stat[
        "by_srna"
    ].items():

        print(

            f"[RISULTATO] {sid}: "

            f"siti={value['candidate_transcripts']}, "

            f"match_energy_favorevole="
            f"{value['favorable_matched_energy']}, "

            f"mediana_delta="
            f"{value['median_delta_kcal_mol']}, "

            f"copertura_mediana="
            f"{value['median_site_coverage']:.3f}"
        )

    # ========================================================
    # 22. MODALITA' AUDIT
    # ========================================================

    if args.audit:

        print(
            "[OK] Audit senza scrittura; "
            "GUI e modelli invariati"
        )

        return

    # ========================================================
    # 23. NUOVA CARTELLA OUTPUT
    # ========================================================

    timestamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    destination = RUNS / (
        "STEP25_ENERGIA_COPERTURA_"
        + timestamp
    )

    destination.mkdir(
        parents=True,
        exist_ok=False
    )

    # ========================================================
    # 24. TSV DEL CONFRONTO APPAIATO
    # ========================================================

    outfile = (
        destination
        / "CONFRONTO_APPAIATO_36_TRASCRITTI.tsv"
    )

    with outfile.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=list(rows[0]),
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(
            rows
        )

    # ========================================================
    # 25. RAPPORTO JSON
    # ========================================================

    report = {

        "step": "25",

        "status":
            "EXPLORATORY_COMPARISON",

        "source_step24":
            str(
                report_path.relative_to(ROOT)
            ),

        "source_sha256": {

            key: sha(path)

            for key, path in inputs.items()
        },

        "statistics":
            stat,

        "paired_tsv":
            str(
                outfile.relative_to(ROOT)
            ),

        "limitations": [

            "Confronto condizionato ai "
            "trascritti preselezionati.",

            "Le energie RNAplex non dimostrano "
            "target biologici.",

            "Sign test e bootstrap sono "
            "esplorativi e non correggono "
            "la selezione dei trascritti.",

            "La mancanza di significativita' "
            "RNA-seq B73/H21 resta invariata.",

            "Accessibilita' e validazione "
            "sperimentale non valutate."
        ],

        "gui_modified":
            False,

        "docker_modified":
            False,

        "training_performed":
            False
    }

    jsonfile = (
        destination
        / "RAPPORTO_STEP25.json"
    )

    jsonfile.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ) + "\n",
        encoding="utf-8"
    )

    # ========================================================
    # 26. RAPPORTO MARKDOWN
    # ========================================================

    mdfile = (
        destination
        / "SINTESI_STEP25.md"
    )

    lines = [

        "# AI PRIMETTA V3.3 — STEP 25",

        "",

        "Analisi esplorativa dei 37 trascritti "
        "preselezionati.",

        "",

        "Non costituisce validazione biologica.",

        "",

        "## Statistiche"
    ]

    for key, value in stat.items():

        if key != "by_srna":

            lines.append(
                f"- {key}: {value}"
            )

    mdfile.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8"
    )

    # ========================================================
    # 27. ESITO FINALE
    # ========================================================

    print(
        "[TSV]",
        outfile
    )

    print(
        "[REPORT]",
        jsonfile
    )

    print(
        "[MARKDOWN]",
        mdfile
    )

    print(
        "[OK] Nessuna GUI, Docker, modello "
        "o file scientifico originale modificato"
    )


# ============================================================
# 28. AVVIO
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