#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 30

Registro consolidato delle evidenze scientifiche.

--audit
    Verifica in sola lettura STEP27, STEP28, STEP29.

--run
    Genera:
    - REGISTRO_36_CANDIDATI.tsv
    - RAPPORTO_SCIENTIFICO_STEP30.md
    - RAPPORTO_STEP30.json

Non esegue RNAplex, RNAplfold o training ML/GNN.
Non modifica GUI, Docker o dati originali.
"""

import argparse
import csv
import hashlib
import json
import math
import statistics
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

EXPECTED = {
    "sRNA_00000037": (17, 25),
    "sRNA_00000435": (19, 24)
}


# ============================================================
# 2. VERIFICA CONDIZIONI
# ============================================================

def require(condition, message):

    if not condition:
        raise ValueError(message)


# ============================================================
# 3. VERIFICA PERCORSI
# ============================================================

def file_in_runs(name):

    p = Path(name)

    if not p.is_absolute():
        p = ROOT / p

    p = p.resolve()

    require(
        p.is_file()
        and p.is_relative_to(RUNS.resolve()),
        f"File mancante o esterno a runs: {p}"
    )

    return p


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
# 5. LETTURA TSV
# ============================================================

def read_tsv(path, needed):

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        reader = csv.DictReader(
            stream,
            delimiter="\t"
        )

        columns = reader.fieldnames or []

        require(
            len(columns) == len(set(columns))
            and set(needed) <= set(columns),
            f"Colonne errate in {path}: "
            f"{set(needed) - set(columns)}"
        )

        records = list(reader)

        require(
            records
            and all(
                None not in r
                for r in records
            ),
            f"TSV vuoto o malformato: {path}"
        )

        return records


# ============================================================
# 6. CONVERSIONE NUMERICA
# ============================================================

def asfloat(value, field):

    try:
        v = float(value)

    except (TypeError, ValueError) as exc:

        raise ValueError(
            f"Valore non numerico {field}: {value!r}"
        ) from exc

    require(
        math.isfinite(v),
        f"Valore non finito {field}"
    )

    return v


# ============================================================
# 7. INDICE UNIVOCO
# ============================================================

def index_unique(records, keyfields):

    found = {}

    for r in records:

        key = tuple(
            r[k]
            for k in keyfields
        )

        require(
            key not in found,
            f"Identificativo duplicato: {key}"
        )

        found[key] = r

    return found


# ============================================================
# 8. RICERCA AUTOMATICA STEP29
# ============================================================

def locate_step29():

    reports = sorted(
        RUNS.glob(
            "STEP29_ACCESSIBILITA_RELATIVA_*/"
            "RAPPORTO_STEP29.json"
        )
    )

    require(
        reports,
        "RAPPORTO_STEP29.json non trovato"
    )

    return file_in_runs(
        reports[-1]
    )


# ============================================================
# 9. CARICAMENTO DEI REPORT
# ============================================================

def load_data():

    p29 = locate_step29()

    r29 = json.loads(
        p29.read_text(
            encoding="utf-8"
        )
    )

    require(
        str(r29.get("step")) == "29"
        and r29.get("sites") == 36
        and r29.get("transcripts") == 37,
        "Manifest STEP29 inatteso"
    )

    # ========================================================
    # 10. PROVENIENZA STEP27 / STEP28
    # ========================================================

    p27 = file_in_runs(
        r29["source_step27"]
    )

    p28 = file_in_runs(
        r29["source_step28"]
    )

    require(
        sha256(p27)
        == r29["source_reports_sha256"]["step27"]
        and
        sha256(p28)
        == r29["source_reports_sha256"]["step28"],
        "Report STEP27/28 modificato dopo STEP29"
    )

    r27 = json.loads(
        p27.read_text(
            encoding="utf-8"
        )
    )

    r28 = json.loads(
        p28.read_text(
            encoding="utf-8"
        )
    )

    require(
        str(r27.get("step")) == "27"
        and str(r28.get("step")) == "28",
        "Step sorgente non corretto"
    )

    # ========================================================
    # 11. RISULTATI STEP29
    # ========================================================

    names29 = {
        Path(x).name: file_in_runs(x)
        for x in r29["outputs"]
    }

    require(
        set(names29) == {
            "ACCESSIBILITA_RELATIVA_36_SITI.tsv",
            "RIEPILOGO_ACCESSIBILITA_PER_GENE.tsv"
        },
        "File STEP29 non corrispondenti"
    )

    p_sites = names29[
        "ACCESSIBILITA_RELATIVA_36_SITI.tsv"
    ]

    p_gene = names29[
        "RIEPILOGO_ACCESSIBILITA_PER_GENE.tsv"
    ]

    # ========================================================
    # 12. RISULTATI STEP28
    # ========================================================

    p28sites = file_in_runs(
        r28["output_tsv"]
    )

    # ========================================================
    # 13. RISULTATI STEP27
    # ========================================================

    p27pairs = [
        file_in_runs(x)

        for x in r27["outputs"]

        if Path(x).name
        == "CONFRONTO_74_PVALUE_QVALUE.tsv"
    ]

    require(
        len(p27pairs) == 1,
        "Tabella 74 coppie STEP27 "
        "assente o duplicata"
    )

    # ========================================================
    # 14. LETTURA DATI STEP29
    # ========================================================

    a29 = read_tsv(
        p_sites,
        [
            "srna_id",
            "transcript_id",
            "gene_id",
            "site_start_1based",
            "site_end_1based",
            "P10_site_mean",
            "percentile_accessibilita_nel_trascritto",
            "RNAplex_MFE_kcal_mol",
            "RNAplex_empirical_p",
            "RNAplex_BH_q_74",
            "biological_status"
        ]
    )

    # ========================================================
    # 15. LETTURA DATI STEP28
    # ========================================================

    a28 = read_tsv(
        p28sites,
        [
            "srna_id",
            "transcript_id",
            "gene_id",
            "target_start_1based",
            "target_end_1based",
            "measurement_complete",
            "mean_P_unpaired_10nt",
            "RNAplex_MFE_kcal_mol"
        ]
    )

    # ========================================================
    # 16. LETTURA DATI STEP27
    # ========================================================

    a27 = read_tsv(
        p27pairs[0],
        [
            "srna_id",
            "transcript_id",
            "gene_id",
            "STEP22_candidate",
            "energy_original_kcal_mol",
            "p_empirical_31",
            "BH_q_74_exploratory"
        ]
    )

    # ========================================================
    # 17. RIEPILOGO PER GENE STEP29
    # ========================================================

    g29 = read_tsv(
        p_gene,
        [
            "srna_id",
            "gene_id",
            "n_trascritti",
            "mediana_percentile_accessibilita"
        ]
    )

    # ========================================================
    # 18. CONTROLLO NUMERO RIGHE
    # ========================================================

    require(
        len(a29) == 36
        and len(a28) == 36
        and len(a27) == 74,
        "Numero righe STEP27/28/29 errato"
    )

    require(
        r28.get("transcripts_processed") == 37
        and r28.get("sites_processed") == 36
        and len(r28.get("profiles", {})) == 37,
        "Profili STEP28 incompleti"
    )

    # ========================================================
    # 19. INDICI DELLE COPPIE
    # ========================================================

    idx29 = index_unique(
        a29,
        ("srna_id", "transcript_id")
    )

    idx28 = index_unique(
        a28,
        ("srna_id", "transcript_id")
    )

    idx27 = index_unique(
        a27,
        ("srna_id", "transcript_id")
    )

    require(
        set(idx29) == set(idx28)
        and set(idx29) <= set(idx27),
        "Coppie STEP27/28/29 discordanti"
    )

    require(
        Counter(
            k[0]
            for k in idx29
        )
        ==
        Counter({
            sid: values[0]
            for sid, values in EXPECTED.items()
        }),
        "Numero candidati per sRNA diverso da 17/19"
    )

    require(
        sum(
            r["STEP22_candidate"].strip().lower()
            == "yes"
            for r in a27
        ) == 36,
        "STEP27: non risultano 36 candidati"
    )

    # ========================================================
    # 20. CONFRONTO SCIENTIFICO PER COPPIA
    # ========================================================

    final = []

    for sid, tx in sorted(idx29):

        x = idx29[(sid, tx)]
        y = idx28[(sid, tx)]
        z = idx27[(sid, tx)]

        require(
            z["STEP22_candidate"].strip().lower()
            == "yes",
            f"Coppia non candidata STEP27: {sid}/{tx}"
        )

        # ====================================================
        # 21. VERIFICA IDENTIFICATIVO GENE
        # ====================================================

        require(
            x["gene_id"]
            == y["gene_id"]
            == z["gene_id"],
            f"Gene discordante: {sid}/{tx}"
        )

        # ====================================================
        # 22. VERIFICA COORDINATE
        # ====================================================

        a = int(
            x["site_start_1based"]
        )

        b = int(
            x["site_end_1based"]
        )

        require(
            a >= 1
            and b >= a
            and b-a+1 == EXPECTED[sid][1]
            and a == int(
                y["target_start_1based"]
            )
            and b == int(
                y["target_end_1based"]
            ),
            f"Coordinate discordanti: {sid}/{tx}"
        )

        require(
            y["measurement_complete"]
            .strip()
            .lower() == "yes",
            f"Misura incompleta: {sid}/{tx}"
        )

        # ====================================================
        # 23. CONVERSIONE DEI VALORI SCIENTIFICI
        # ====================================================

        p10 = asfloat(
            x["P10_site_mean"],
            "P10"
        )

        percentile = asfloat(
            x[
                "percentile_accessibilita_nel_trascritto"
            ],
            "percentile"
        )

        q = asfloat(
            z["BH_q_74_exploratory"],
            "BH"
        )

        p = asfloat(
            z["p_empirical_31"],
            "p empirico"
        )

        energy = asfloat(
            z["energy_original_kcal_mol"],
            "energia"
        )

        require(
            0 <= p10 <= 1
            and 0 <= percentile <= 100
            and 0 <= p <= 1
            and 0 <= q <= 1,
            f"Misure fuori intervallo: {sid}/{tx}"
        )

        # ====================================================
        # 24. CONFRONTO STEP27 / STEP28 / STEP29
        # ====================================================

        require(
            abs(
                p10 - asfloat(
                    y["mean_P_unpaired_10nt"],
                    "P10 STEP28"
                )
            ) < 5e-8

            and

            abs(
                energy - asfloat(
                    x["RNAplex_MFE_kcal_mol"],
                    "energia STEP29"
                )
            ) < 1e-5

            and

            abs(
                energy - asfloat(
                    y["RNAplex_MFE_kcal_mol"],
                    "energia STEP28"
                )
            ) < 1e-5

            and

            abs(
                p - asfloat(
                    x["RNAplex_empirical_p"],
                    "p STEP29"
                )
            ) < 2e-8

            and

            abs(
                q - asfloat(
                    x["RNAplex_BH_q_74"],
                    "q STEP29"
                )
            ) < 2e-8,

            f"Valori STEP27/28/29 discordanti: "
            f"{sid}/{tx}"
        )

        # ====================================================
        # 25. STATO BIOLOGICO
        # ====================================================

        require(
            x["biological_status"]
            == "exploratory_not_validated",
            f"Stato biologico inatteso: {sid}/{tx}"
        )

        # ====================================================
        # 26. RECORD INTEGRATO
        # ====================================================

        final.append({

            "srna_id":
                sid,

            "transcript_id":
                tx,

            "gene_id":
                x["gene_id"],

            "site_start_1based":
                a,

            "site_end_1based":
                b,

            "RNAplex_energy_kcal_mol":
                z["energy_original_kcal_mol"],

            "shuffle_empirical_p":
                z["p_empirical_31"],

            "BH_q_74_exploratory":
                z["BH_q_74_exploratory"],

            "mean_P10":
                x["P10_site_mean"],

            "within_transcript_percentile":
                x[
                    "percentile_accessibilita_nel_trascritto"
                ],

            "high_relative_accessibility_ge75":
                (
                    "yes"
                    if percentile >= 75
                    else "no"
                ),

            "biological_status":
                "ipotesi_computazionale_non_validata"
        })

    # ========================================================
    # 27. VERIFICA RAGGRUPPAMENTO PER GENE
    # ========================================================

    groups = defaultdict(list)

    for r in final:

        groups[
            (r["srna_id"], r["gene_id"])
        ].append(
            float(
                r["within_transcript_percentile"]
            )
        )

    idx_gene = index_unique(
        g29,
        ("srna_id", "gene_id")
    )

    require(
        set(idx_gene) == set(groups),
        "Aggregazione per gene STEP29 discordante"
    )

    for key, values in groups.items():

        source = idx_gene[key]

        require(
            int(source["n_trascritti"])
            == len(values)

            and

            abs(
                asfloat(
                    source[
                        "mediana_percentile_accessibilita"
                    ],
                    "mediana gene"
                )
                - statistics.median(values)
            ) <= 0.002,

            f"Mediana per gene incoerente: {key}"
        )

    # ========================================================
    # 28. RIEPILOGO PER SMALL RNA
    # ========================================================

    stats = {}

    for sid in EXPECTED:

        items = [
            r for r in final
            if r["srna_id"] == sid
        ]

        percentiles = [
            float(
                r["within_transcript_percentile"]
            )
            for r in items
        ]

        stats[sid] = {

            "sites":
                len(items),

            "genes":
                len({
                    r["gene_id"]
                    for r in items
                }),

            "median_percentile":
                round(
                    statistics.median(percentiles),
                    3
                ),

            "sites_percentile_ge75":
                sum(
                    x >= 75
                    for x in percentiles
                ),

            "sites_bh_q_below_005":
                sum(
                    float(
                        r["BH_q_74_exploratory"]
                    ) < 0.05
                    for r in items
                )
        }

        # ====================================================
        # 29. VERIFICA CON RAPPORTO STEP29
        # ====================================================

        prior = (
            r29["summary_by_srna"][sid]
        )

        require(
            prior["sites"]
            == stats[sid]["sites"]

            and

            prior["unique_genes"]
            == stats[sid]["genes"]

            and

            prior["sites_percentile_ge_75"]
            == stats[sid]["sites_percentile_ge75"]

            and

            abs(
                prior[
                    "median_within_transcript_percentile"
                ]
                - stats[sid]["median_percentile"]
            ) <= 0.002,

            f"Riepilogo STEP29 discordante: {sid}"
        )

    # ========================================================
    # 30. CORREZIONE BH SULLE 74 COPPIE
    # ========================================================

    n_bh74 = sum(
        asfloat(
            r["BH_q_74_exploratory"],
            "BH-74"
        ) < 0.05
        for r in a27
    )

    stats[
        "all_74_pairs_BH_significant"
    ] = n_bh74

    return (
        p29,
        p27,
        p28,
        [
            p_sites,
            p_gene,
            p28sites,
            p27pairs[0]
        ],
        final,
        stats
    )


# ============================================================
# 31. SCRITTURA TSV
# ============================================================

def save_tsv(path, records):

    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as stream:

        writer = csv.DictWriter(
            stream,
            fieldnames=list(records[0]),
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(records)


# ============================================================
# 32. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = (
        parser.add_mutually_exclusive_group(
            required=True
        )
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

    print("=" * 70)

    print(
        "AI PRIMETTA V3.3 | STEP 30 "
        "- REGISTRO SCIENTIFICO"
    )

    print("=" * 70)

    # ========================================================
    # 33. AUDIT INTEGRATO
    # ========================================================

    (
        p29,
        p27,
        p28,
        inputs,
        records,
        stats
    ) = load_data()

    print(
        "[PASS] Coerenza STEP27/28/29; "
        "74 coppie, 36 candidati, 37 profili"
    )

    print(
        "[BH-74] Coppie con q<0.05:",
        stats[
            "all_74_pairs_BH_significant"
        ]
    )

    for sid in EXPECTED:

        s = stats[sid]

        print(
            f"[RISULTATO] {sid}: "
            f"siti={s['sites']}, "
            f"geni={s['genes']}, "
            f"mediana_percentile="
            f"{s['median_percentile']}, "
            f"accessibilita_ge75="
            f"{s['sites_percentile_ge75']}, "
            f"coppie_BH_significative="
            f"{s['sites_bh_q_below_005']}"
        )

    print(
        "[STATO] Tutte le 36 relazioni "
        "restano ipotesi computazionali"
    )

    # ========================================================
    # 34. MODALITA AUDIT
    # ========================================================

    if args.audit:

        print(
            "[OK] Audit in sola lettura"
        )

        return

    # ========================================================
    # 35. CARTELLA OUTPUT
    # ========================================================

    tag = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    dest = RUNS / (
        f"STEP30_REGISTRO_SCIENTIFICO_{tag}"
    )

    dest.mkdir(
        parents=True,
        exist_ok=False
    )

    tab = (
        dest
        / "REGISTRO_36_CANDIDATI.tsv"
    )

    rep = (
        dest
        / "RAPPORTO_SCIENTIFICO_STEP30.md"
    )

    manifest = (
        dest
        / "RAPPORTO_STEP30.json"
    )

    # ========================================================
    # 36. SALVATAGGIO REGISTRO
    # ========================================================

    save_tsv(
        tab,
        records
    )

    # ========================================================
    # 37. RAPPORTO SCIENTIFICO MARKDOWN
    # ========================================================

    report_lines = [

        "# AI PRIMETTA V3.3 — STEP 30",

        "",

        "## Chiusura della fase esplorativa",

        "",

        "Registro delle evidenze computazionali "
        "per 36 siti candidati su 37 trascritti B73.",

        "Non dimostra l'interazione biologica "
        "delle small RNA con i geni.",

        "",

        "## Risultati",

        "",

        "| Small RNA | Siti | Geni | "
        "Mediana percentile P10 | "
        "Siti >= 75° | BH q<0,05 |",

        "|---|---:|---:|---:|---:|---:|"
    ]

    for sid in EXPECTED:

        s = stats[sid]

        report_lines.append(

            f"| {sid} "

            f"| {s['sites']} "

            f"| {s['genes']} "

            f"| {s['median_percentile']:.3f} "

            f"| {s['sites_percentile_ge75']} "

            f"| {s['sites_bh_q_below_005']} |"
        )

    report_lines += [

        "",

        "Tra tutte le 74 coppie RNAplex, "
        "quelle con BH q<0,05 sono: "
        f"{stats['all_74_pairs_BH_significant']}.",

        "",

        "## Interpretazione e limitazioni",

        "",

        "- I trascritti erano già stati "
        "preselezionati per complementarità.",

        "- I p empirici derivano da 30 shuffle "
        "per small RNA; la risoluzione "
        "minima è 1/31.",

        "- Le isoforme di uno stesso gene "
        "non equivalgono a repliche "
        "biologiche indipendenti.",

        "- P10 descrive la probabilità "
        "di apertura di segmenti di 10 nt, "
        "non dell'intero sito.",

        "- Il percentile è relativo alle "
        "finestre dello stesso trascritto; "
        "non è un p-value.",

        "- I dati RNA-seq STEP21B non "
        "sono stati ricalcolati in "
        "questo STEP30. Per la loro "
        "interpretazione riferirsi "
        "all'audit originale.",

        "- Nessuna interazione è stata "
        "classificata come target "
        "biologicamente validato.",

        "- Nessun training ML/GNN e "
        "nessuna modifica a GUI e Docker.",

        "",

        "## Prossima decisione scientifica",

        "",

        "Verificare la classe biologica "
        "delle small RNA e la qualità "
        "delle annotazioni prima di "
        "utilizzare questi siti come "
        "possibili archi di una GNN.",

        ""
    ]

    rep.write_text(
        "\n".join(report_lines),
        encoding="utf-8"
    )

    # ========================================================
    # 38. MANIFEST DI PROVENIENZA
    # ========================================================

    payload = {

        "step":
            30,

        "status":
            "EXPLORATORY_REGISTER_NOT_VALIDATED",

        "source_reports": {

            "step27":
                str(p27.relative_to(ROOT)),

            "step28":
                str(p28.relative_to(ROOT)),

            "step29":
                str(p29.relative_to(ROOT))
        },

        "inputs_sha256": {

            str(p.relative_to(ROOT)):
                sha256(p)

            for p in [
                p27,
                p28,
                p29,
                *inputs
            ]
        },

        "summary":
            stats,

        "candidate_count":
            len(records),

        "outputs": [

            str(tab.relative_to(ROOT)),
            str(rep.relative_to(ROOT))
        ],

        "rnaseq_reanalyzed":
            False,

        "rna_interactions_biologically_validated":
            0,

        "training_performed":
            False,

        "gui_modified":
            False,

        "docker_modified":
            False
    }

    manifest.write_text(

        json.dumps(
            payload,
            indent=2,
            ensure_ascii=False
        ) + "\n",

        encoding="utf-8"
    )

    # ========================================================
    # 39. OUTPUT FINALE
    # ========================================================

    print(
        "[TSV]",
        tab
    )

    print(
        "[REPORT]",
        rep
    )

    print(
        "[MANIFEST]",
        manifest
    )

    print(
        "[OK] STEP30 completato; "
        "sorgenti, GUI e Docker invariati"
    )


# ============================================================
# 40. AVVIO
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