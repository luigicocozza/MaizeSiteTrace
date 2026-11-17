
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AI PRIMETTA V3.3 - STEP 21B: correzione delle evidenze RNA-seq.

--audit: verifica fonti, schema e integrita' senza scrivere.
--run: crea una NUOVA cartella e report, senza modificare sorgenti.

I qvalue sono usati per la soglia statistica.
Il campo fdr_response viene conservato come dato grezzo,
senza reinterpretarlo come valore numerico di FDR.
"""

import argparse
import csv
import hashlib
import json
import sys

from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

RNASEQ = (
    ROOT
    / "outputs"
    / "prjna419326"
    / "PRJNA419326_rnaseq_gene_DE_combined_B73_H21.tsv"
)

THRESHOLD = Decimal("0.05")
SCORE_LIMIT = Decimal("5")

PREFIXES = ("B73", "H21")

BASE_REQUIRED = (
    "srna_id",
    "transcript_id",
    "gene_id",
    "start_1based",
    "end_1based",
    "full_score",
    "full_hard",
    "full_gu",
)

RNA_REQUIRED = (
    "gene_id",
) + tuple(
    item
    for p in PREFIXES
    for item in (
        f"{p}_qvalue",
        f"{p}_log2FC_TPM_DS_vs_WW",
    )
)


# ==========================================================
# 1. CHECKSUM
# ==========================================================

def sha256(path):

    digest = hashlib.sha256()

    with path.open("rb") as handle:

        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b""
        ):
            digest.update(chunk)

    return digest.hexdigest()


# ==========================================================
# 2. LETTURA TSV
# ==========================================================

def read_tsv(path):

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        fields = reader.fieldnames or []

        if (
            not fields
            or len(set(fields)) != len(fields)
        ):
            raise ValueError(
                "Intestazioni TSV assenti o duplicate: "
                + str(path)
            )

        rows = list(reader)

    return fields, rows


# ==========================================================
# 3. CONVERSIONE NUMERICA SICURA
# ==========================================================

def decimal_optional(
    value,
    field,
    qvalue=False
):

    string = str(value or "").strip()

    if string.lower() in (
        "",
        "na",
        "nan",
        "null",
        "none",
        "-"
    ):
        return None

    try:

        result = Decimal(
            string.replace(",", ".")
        )

    except InvalidOperation as exc:

        raise ValueError(
            f"Valore non numerico in {field}: "
            f"{string!r}"
        ) from exc

    if not result.is_finite():

        raise ValueError(
            f"Valore non finito in {field}: "
            f"{string!r}"
        )

    if qvalue and not (
        Decimal(0) <= result <= Decimal(1)
    ):

        raise ValueError(
            f"qvalue fuori intervallo [0,1] "
            f"in {field}: {result}"
        )

    return result


# ==========================================================
# 4. RECUPERO DEI FILE ORIGINALI
# ==========================================================

def origin_files():

    reports = sorted(
        RUNS.glob(
            "STEP21_TRIAGE_TARGET_SRNA_*/"
            "RAPPORTO_STEP21.json"
        )
    )

    if not reports:

        raise FileNotFoundError(
            "Nessun RAPPORTO_STEP21.json disponibile"
        )

    path = reports[-1]

    report = json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )

    if report.get("step") != "21":

        raise ValueError(
            "Report STEP21 non valido"
        )

    step20 = (
        ROOT
        / report["source_step20_tsv"]
    )

    if (
        not step20.is_file()
        or not step20.resolve().is_relative_to(
            RUNS.resolve()
        )
    ):

        raise ValueError(
            "TSV STEP20 assente o fuori "
            "dalla cartella runs"
        )

    if sha256(step20) != report[
        "sha256_step20_tsv"
    ]:

        raise ValueError(
            "Il TSV STEP20 e' cambiato "
            "dal precedente STEP21"
        )

    candidate = (
        path.parent
        / "SITI_FULL_SCORE_LE5_ANNOTATI.tsv"
    )

    if not candidate.is_file():

        raise FileNotFoundError(candidate)

    if not RNASEQ.is_file():

        raise FileNotFoundError(RNASEQ)

    return (
        path,
        report,
        step20,
        candidate
    )


# ==========================================================
# 5. PREPARAZIONE E VERIFICA SCIENTIFICA
# ==========================================================

def prepare():

    (
        path,
        previous,
        step20,
        old
    ) = origin_files()

    full_columns, all_sites = (
        read_tsv(step20)
    )

    _, old_rows = read_tsv(old)

    (
        expression_columns,
        expression_rows
    ) = read_tsv(RNASEQ)

    missing = (
        set(BASE_REQUIRED)
        - set(full_columns)
    )

    if missing:

        raise ValueError(
            "Colonne STEP20 mancanti: "
            + repr(sorted(missing))
        )

    missing = (
        set(RNA_REQUIRED)
        - set(expression_columns)
    )

    if missing:

        raise ValueError(
            "Colonne RNA-seq mancanti: "
            + repr(sorted(missing))
        )

    # Ricostruzione dei candidati partendo
    # dallo STEP20, evitando di riutilizzare
    # le annotazioni FDR errate dello STEP21.

    selected = []
    keys = set()

    for row in all_sites:

        score = decimal_optional(
            row["full_score"],
            "full_score"
        )

        hard = int(
            row["full_hard"]
        )

        gu = int(
            row["full_gu"]
        )

        expected = (
            hard
            + Decimal("0.5") * gu
        )

        if (
            score is None
            or score != expected
        ):

            raise ValueError(
                "Score full-length incoerente"
            )

        if score > SCORE_LIMIT:
            continue

        key = tuple(
            row[field]
            for field in (
                "srna_id",
                "transcript_id",
                "start_1based",
                "end_1based"
            )
        )

        if key in keys:

            raise ValueError(
                "Sito duplicato: "
                + repr(key)
            )

        keys.add(key)

        selected.append(row)

    # Verifica che i candidati siano
    # gli stessi dello STEP21.

    old_keys = {

        tuple(
            row[field]
            for field in (
                "srna_id",
                "transcript_id",
                "start_1based",
                "end_1based"
            )
        )

        for row in old_rows
    }

    if (
        len(old_keys) != len(old_rows)
        or keys != old_keys
    ):

        raise ValueError(
            "I siti prioritari STEP20 "
            "e STEP21 non coincidono"
        )

    if len(selected) != int(
        previous["candidate_site_count"]
    ):

        raise ValueError(
            "Numero candidati diverso "
            "dal rapporto STEP21"
        )

    # ======================================================
    # 6. INDICE RNA-SEQ PER GENE
    # ======================================================

    by_gene = {}

    duplicate_ids = set()

    for row in expression_rows:

        gene = (
            row.get("gene_id") or ""
        ).strip()

        if not gene:
            continue

        if gene in by_gene:

            duplicate_ids.add(gene)

        else:

            by_gene[gene] = row

    # Non utilizzare ID ambigui.

    for gene in duplicate_ids:

        del by_gene[gene]

    # ======================================================
    # 7. ASSOCIAZIONE ESATTA
    # ======================================================

    result = []

    matches = 0

    missing_genes = set()

    for row in selected:

        gene = (
            row.get("gene_id") or ""
        ).strip()

        expression = by_gene.get(gene)

        if expression is not None:

            status = "exact_gene_id_found"

            matches += 1

        elif not gene:

            status = "missing_gene_id"

        elif gene in duplicate_ids:

            status = (
                "ambiguous_duplicate_rnaseq_gene"
            )

        else:

            status = "gene_not_in_rnaseq"

            missing_genes.add(gene)

        out = dict(row)

        out[
            "rnaseq_join_status_verified"
        ] = status

        # ================================================
        # 8. EVIDENZE B73 E H21
        # ================================================

        for p in PREFIXES:

            q = (
                decimal_optional(
                    expression[f"{p}_qvalue"],
                    f"{p}_qvalue",
                    True
                )

                if expression
                else None
            )

            fc = (
                decimal_optional(
                    expression[
                        f"{p}_log2FC_TPM_DS_vs_WW"
                    ],
                    f"{p}_log2FC_TPM_DS_vs_WW"
                )

                if expression
                else None
            )

            out[f"{p}_qvalue"] = (
                ""
                if q is None
                else str(q)
            )

            out[
                f"{p}_log2FC_TPM_DS_vs_WW"
            ] = (
                ""
                if fc is None
                else str(fc)
            )

            # qvalue < 0.05, in senso stretto.

            out[
                f"{p}_significant_q_lt_005"
            ] = (

                "not_available"
                if q is None

                else "yes"
                if q < THRESHOLD

                else "no"
            )

            # Direzione dello stress idrico:
            # DS rispetto a WW.

            out[
                f"{p}_direction_DS_vs_WW"
            ] = (

                "not_available"
                if fc is None

                else "up"
                if fc > 0

                else "down"
                if fc < 0

                else "unchanged"
            )

            # Conserva altri dati di origine.
            # fdr_response NON sostituisce qvalue.

            for key in (
                f"{p}_fdr_response",
                f"{p}_expressed_filter",
                f"{p}_mean_TPM_DS",
                f"{p}_mean_TPM_WW"
            ):

                out[f"{key}_raw"] = (

                    (
                        expression.get(key)
                        or ""
                    ).strip()

                    if expression
                    else ""
                )

        result.append(out)

    # ======================================================
    # 9. CONTROLLI INCROCIATI CON STEP21
    # ======================================================

    if matches != int(
        previous["rnaseq_exact_match_sites"]
    ):

        raise ValueError(
            "Numero associazioni RNA-seq "
            "diverso dal rapporto STEP21"
        )

    n_genes = len({

        row["gene_id"]

        for row in result

        if row["gene_id"]
    })

    if n_genes != int(
        previous["unique_annotated_genes"]
    ):

        raise ValueError(
            "Numero geni diverso "
            "dal rapporto STEP21"
        )

    result.sort(
        key=lambda r: (
            Decimal(r["full_score"]),
            r["srna_id"],
            r["transcript_id"],
            int(r["start_1based"])
        )
    )

    return {

        "step21_report": path,
        "step21": previous,

        "step20_file": step20,
        "candidate_file": old,

        "rnaseq_file": RNASEQ,

        "rnaseq_rows": len(
            expression_rows
        ),

        "rnaseq_duplicates": sorted(
            duplicate_ids
        ),

        "unmapped_genes": sorted(
            missing_genes
        ),

        "full_columns": full_columns,

        "rows": result,

        "exact_matches": matches,

        "genes": n_genes
    }


# ==========================================================
# 10. AUDIT
# ==========================================================

def audit(info):

    print(
        "[STEP21]",
        info["step21_report"].relative_to(ROOT)
    )

    print(
        "[STEP20]",
        info["step20_file"].relative_to(ROOT)
    )

    print(
        "[RNA-SEQ]",
        info["rnaseq_file"].relative_to(ROOT)
    )

    print(
        "[INTEGRITA] checksum STEP20 "
        "e lista siti STEP21: PASS"
    )

    print(
        "[SITI]",
        len(info["rows"])
    )

    print(
        "[GENI]",
        info["genes"]
    )

    print(
        "[JOIN RNA-SEQ]",
        info["exact_matches"],
        "/",
        len(info["rows"])
    )

    print(
        "[RNA-SEQ] righe:",
        info["rnaseq_rows"],

        "| gene ID duplicati:",
        len(info["rnaseq_duplicates"])
    )

    print(
        "[GENI NON CORRISPONDENTI]",
        info["unmapped_genes"] or "nessuno"
    )

    print(
        "[SORGENTI FDR]",
        "B73_qvalue, H21_qvalue"
    )

    print(
        "[SORGENTI log2FC]",
        "B73_log2FC_TPM_DS_vs_WW, "
        "H21_log2FC_TPM_DS_vs_WW"
    )

    for p in PREFIXES:

        counts = Counter(

            row[
                f"{p}_significant_q_lt_005"
            ]

            for row in info["rows"]
        )

        print(
            f"[{p} q<0.05] "
            f"si={counts['yes']}, "
            f"no={counts['no']}, "
            f"non disponibile="
            f"{counts['not_available']} "
            "(conteggio SITI)"
        )

    print(
        "[OK] Sorgenti e provenienza verificate; "
        "GUI e Docker invariati"
    )


# ==========================================================
# 11. ESPORTAZIONE CORRETTA
# ==========================================================

def export(info):

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    dest = RUNS / (
        "STEP21B_RNASEQ_CORRETTO_"
        + stamp
    )

    dest.mkdir(
        parents=True,
        exist_ok=False
    )

    # ======================================================
    # 12. DETTAGLIO DEI SITI
    # ======================================================

    extension = [
        "rnaseq_join_status_verified"
    ]

    for p in PREFIXES:

        extension.extend([

            f"{p}_qvalue",

            f"{p}_log2FC_TPM_DS_vs_WW",

            f"{p}_significant_q_lt_005",

            f"{p}_direction_DS_vs_WW",

            f"{p}_fdr_response_raw",

            f"{p}_expressed_filter_raw",

            f"{p}_mean_TPM_DS_raw",

            f"{p}_mean_TPM_WW_raw"
        ])

    columns = (
        info["full_columns"]
        + extension
    )

    detail = (
        dest
        / "SITI_PRIORITARI_RNASEQ_CORRETTI.tsv"
    )

    with detail.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=columns
        )

        writer.writeheader()

        writer.writerows(
            info["rows"]
        )

    # ======================================================
    # 13. RAGGRUPPAMENTO sRNA-GENE
    # ======================================================

    groups = defaultdict(list)

    for row in info["rows"]:

        key = (

            row["srna_id"],

            row["gene_id"]

            or (
                "NO_GENE::"
                + row["transcript_id"]
            )
        )

        groups[key].append(row)

    gene_rows = []

    for (
        srna,
        gene_key
    ), rows in sorted(groups.items()):

        best = min(
            rows,
            key=lambda r: (

                Decimal(
                    r["full_score"]
                ),

                int(
                    r["full_hard"]
                ),

                r["transcript_id"]
            )
        )

        item = {

            "srna_id": srna,

            "gene_or_transcript":
                gene_key,

            "gene_id":
                best["gene_id"],

            "site_count":
                len(rows),

            "distinct_transcripts":
                len({

                    r["transcript_id"]

                    for r in rows
                }),

            "best_full_score":
                best["full_score"],

            "rnaseq_join_status_verified":
                best[
                    "rnaseq_join_status_verified"
                ]
        }

        for p in PREFIXES:

            for key in (

                f"{p}_qvalue",

                f"{p}_log2FC_TPM_DS_vs_WW",

                f"{p}_significant_q_lt_005",

                f"{p}_direction_DS_vs_WW"
            ):

                item[key] = best[key]

        gene_rows.append(item)

    gene_path = (
        dest
        / "COPPIE_SRNA_GENE_RNASEQ_CORRETTE.tsv"
    )

    gene_columns = (

        list(gene_rows[0])

        if gene_rows

        else [
            "srna_id",
            "gene_or_transcript",
            "gene_id",
            "site_count"
        ]
    )

    with gene_path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=gene_columns
        )

        writer.writeheader()

        writer.writerows(gene_rows)

    # ======================================================
    # 14. RAPPORTO SCIENTIFICO
    # ======================================================

    report = {

        "step": "21B",

        "status":
            "DESCRIPTIVE_CORRECTION_ONLY",

        "source_step21": str(
            info["step21_report"].relative_to(ROOT)
        ),

        "source_step20_tsv": str(
            info["step20_file"].relative_to(ROOT)
        ),

        "source_step20_sha256": sha256(
            info["step20_file"]
        ),

        "source_rnaseq": str(
            info["rnaseq_file"].relative_to(ROOT)
        ),

        "source_rnaseq_sha256": sha256(
            info["rnaseq_file"]
        ),

        "candidate_sites":
            len(info["rows"]),

        "unique_genes":
            info["genes"],

        "srna_gene_pairs":
            len(gene_rows),

        "sites_exact_gene_match":
            info["exact_matches"],

        "unmatched_genes":
            info["unmapped_genes"],

        "rnaseq_ambiguous_duplicate_gene_ids":
            len(info["rnaseq_duplicates"]),

        "qvalue_threshold":
            "<0.05 (strict)",

        "statistical_columns": [
            "B73_qvalue",
            "H21_qvalue"
        ],

        "fold_change_columns": [
            "B73_log2FC_TPM_DS_vs_WW",
            "H21_log2FC_TPM_DS_vs_WW"
        ],

        "fdr_response_interpretation":
            "raw source field; "
            "NOT used as numeric qvalue",

        "site_level_q_significant_counts": {

            p: dict(
                Counter(

                    row[
                        f"{p}_significant_q_lt_005"
                    ]

                    for row in info["rows"]
                )
            )

            for p in PREFIXES
        },

        "gene_pair_level_q_significant_counts": {

            p: dict(
                Counter(

                    row[
                        f"{p}_significant_q_lt_005"
                    ]

                    for row in gene_rows
                )
            )

            for p in PREFIXES
        },

        "outputs": [

            str(
                detail.relative_to(ROOT)
            ),

            str(
                gene_path.relative_to(ROOT)
            )
        ],

        "limitations": [

            "Differential expression is not "
            "experimental evidence of targeting.",

            "Biological interaction and "
            "small-RNA class remain unverified.",

            "Sites remain conditional on "
            "STEP19 HSP20 ungapped filtering."
        ],

        "modifies_gui": False,

        "modifies_docker": False,

        "modifies_models": False,

        "training_performed": False
    }

    report_path = (
        dest / "RAPPORTO_STEP21B.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    print(
        "[COPPIE S-RNA/GENE]",
        len(gene_rows)
    )

    print(
        "[REPORT]",
        report_path
    )

    print(
        "[SITI]",
        detail
    )

    print(
        "[GENI]",
        gene_path
    )

    print(
        "[OK] STEP21B concluso, "
        "file originali preservati"
    )


# ==========================================================
# 15. MAIN
# ==========================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    group = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    group.add_argument(
        "--audit",
        action="store_true"
    )

    group.add_argument(
        "--run",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 70)

    print(
        "AI PRIMETTA V3.3 | STEP 21B "
        "- EVIDENZE RNA-SEQ CORRETTE"
    )

    print("=" * 70)

    info = prepare()

    audit(info)

    if args.run:
        export(info)


if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__,
            str(exc),
            file=sys.stderr
        )

        sys.exit(1)
