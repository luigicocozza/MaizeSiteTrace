
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 | STEP 21

Annotazione non distruttiva dei siti STEP20.

--audit:
    Mostra file sorgente e schemi disponibili.

--run:
    Controlla tutti i siti e annota quelli
    con full_score <= 5.

NON costituisce validazione sperimentale
di interazioni sRNA-mRNA.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys

from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT / "target_prediction_v33_vscode/runs"
)

RNA_SEQ = (
    ROOT
    / "outputs/prjna419326/"
    "PRJNA419326_rnaseq_gene_DE_combined_B73_H21.tsv"
)

CUTOFF = Decimal("5")

FDR_THRESHOLD = Decimal("0.05")

REQUIRED = {
    "srna_id",
    "srna_length",
    "transcript_id",
    "gene_id",
    "start_1based",
    "end_1based",
    "hsp_score",
    "full_score",
    "full_hard",
    "full_gu",
    "seed_hard",
    "full_site_DNA_5to3",
    "full_perfect_target_DNA_5to3"
}


# ==========================================================
# 1. CHECKSUM
# ==========================================================

def digest(path: Path) -> str:

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# ==========================================================
# 2. VALIDAZIONE VALORI
# ==========================================================

def number(value, field):

    try:

        v = Decimal(
            str(value)
            .strip()
            .replace(",", ".")
        )

        if not v.is_finite():
            raise ValueError("not finite")

        return v

    except (InvalidOperation, ValueError) as e:

        raise ValueError(
            f"Valore non numerico per "
            f"{field}: {value!r}"
        ) from e


def optional_decimal(value):

    s = str(value or "").strip()

    if s.lower() in {
        "",
        "na",
        "nan",
        "none",
        "null",
        "-"
    }:

        return None

    try:

        value = Decimal(
            s.replace(",", ".")
        )

        return (
            value
            if value.is_finite()
            else None
        )

    except InvalidOperation:

        return None


# ==========================================================
# 3. INDIVIDUAZIONE STEP20
# ==========================================================

def latest_step20():

    files = sorted(
        RUNS.glob(
            "STEP20_COMPLEMENTARITA_FULL_LENGTH_*/"
            "RAPPORTO_STEP20.json"
        )
    )

    if not files:

        raise FileNotFoundError(
            "Nessun report STEP20 disponibile"
        )

    path = files[-1]

    report = json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )

    if (
        report.get("step") != "20"
        or report.get("status") != "COMPLETED"
    ):

        raise RuntimeError(
            "Report STEP20 inatteso"
        )

    tsv = ROOT / report["full_tsv"]

    if (
        not tsv.is_file()
        or not tsv.resolve().is_relative_to(
            RUNS.resolve()
        )
    ):

        raise RuntimeError(
            "File TSV STEP20 mancante "
            "o fuori da RUNS"
        )

    return path, tsv, report


# ==========================================================
# 4. LETTURA INTESTAZIONI TSV
# ==========================================================

def columns(path, delimiter="\t"):

    if not path.is_file():
        return []

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        return (
            csv.DictReader(
                stream,
                delimiter=delimiter
            ).fieldnames
            or []
        )


def normalize_header(s):

    return re.sub(
        r"[^a-z0-9]",
        "",
        s.casefold()
    )


# ==========================================================
# 5. IDENTIFICAZIONE COLONNA GENE
# ==========================================================

def identify_gene_column(headers):

    preferences = (
        (
            "geneid",
            "geneidentifier"
        ),
        ("gene",),
        ("genename",)
    )

    for accepted in preferences:

        matches = [
            h
            for h in headers
            if normalize_header(h) in accepted
        ]

        if matches:

            return (
                matches[0]
                if len(matches) == 1
                else None
            ), matches

    return None, []


# ==========================================================
# 6. CARICAMENTO EVIDENZE RNA-SEQ
# ==========================================================

def reference_data():

    if not RNA_SEQ.is_file():

        return {}, {
            "status": "RNAseq_file_not_found",
            "file": str(
                RNA_SEQ.relative_to(ROOT)
            )
        }

    cols = columns(RNA_SEQ)

    gene_col, matches = (
        identify_gene_column(cols)
    )

    info = {

        "status": (
            "exact_gene_join_available"
            if gene_col
            else "ambiguous_or_missing_gene_id"
        ),

        "file": str(
            RNA_SEQ.relative_to(ROOT)
        ),

        "gene_id_column": gene_col,

        "gene_column_candidates": matches,

        "columns": cols,

        "sha256": digest(RNA_SEQ),

        "duplicate_gene_ids": 0,

        "read_rows": 0
    }

    if gene_col is None:

        return {}, info

    genes = {}

    duplicates = set()

    with RNA_SEQ.open(
        encoding="utf-8-sig",
        newline=""
    ) as f:

        for row in csv.DictReader(
            f,
            delimiter="\t"
        ):

            info["read_rows"] += 1

            gene_id = (
                row.get(gene_col) or ""
            ).strip()

            if not gene_id:
                continue

            if gene_id in genes:

                duplicates.add(gene_id)

            else:

                genes[gene_id] = row

    # Gli identificativi duplicati vengono
    # esclusi per evitare join ambigui.

    for gene in duplicates:

        genes.pop(gene, None)

    info["duplicate_gene_ids"] = (
        len(duplicates)
    )

    info["unique_gene_ids_usable"] = (
        len(genes)
    )

    return genes, info


# ==========================================================
# 7. EVIDENZA B73 E H21
# ==========================================================

def expression_evidence(gene, row):

    data = {}

    if not row:

        for prefix in ("B73", "H21"):

            data[f"{prefix}_fdr"] = ""

            data[f"{prefix}_log2fc"] = ""

            data[
                f"{prefix}_DE_FDR_lt_005"
            ] = "unavailable"

        return data

    norms = {
        normalize_header(k): k
        for k in row.keys()
    }

    for prefix in ("B73", "H21"):

        p = prefix.lower()

        fdr_options = [

            f"{p}fdrresponse",

            f"{p}padjresponse",

            f"{p}fdr",

            f"{p}padj"
        ]

        fc_options = [

            f"{p}log2fcresponse",

            f"{p}log2foldchangeresponse",

            f"{p}log2fc",

            f"{p}log2foldchange"
        ]

        fdr_name = next(
            (
                norms[n]
                for n in fdr_options
                if n in norms
            ),
            None
        )

        fc_name = next(
            (
                norms[n]
                for n in fc_options
                if n in norms
            ),
            None
        )

        fdr = (
            optional_decimal(
                row.get(fdr_name)
            )
            if fdr_name
            else None
        )

        fc = (
            optional_decimal(
                row.get(fc_name)
            )
            if fc_name
            else None
        )

        if (
            fdr is not None
            and not 0 <= fdr <= 1
        ):

            fdr = None

        data[f"{prefix}_fdr"] = (
            str(fdr)
            if fdr is not None
            else ""
        )

        data[f"{prefix}_log2fc"] = (
            str(fc)
            if fc is not None
            else ""
        )

        data[
            f"{prefix}_DE_FDR_lt_005"
        ] = (
            (
                "yes"
                if fdr < FDR_THRESHOLD
                else "no"
            )
            if fdr is not None
            else "unavailable"
        )

    return data


# ==========================================================
# 8. AUDIT PREVENTIVO
# ==========================================================

def audit():

    report_path, tsv, report = (
        latest_step20()
    )

    headers = columns(tsv)

    missing = REQUIRED.difference(
        headers
    )

    if missing:

        raise RuntimeError(
            "Colonne mancanti nel TSV: "
            + ", ".join(
                sorted(missing)
            )
        )

    print(
        "[STEP20 REPORT]",
        report_path.relative_to(ROOT)
    )

    print(
        "[STEP20 TSV]",
        tsv.relative_to(ROOT)
    )

    print(
        "[COLONNE SITI]",
        ", ".join(headers)
    )

    print(
        "[RNA-SEQ]",
        RNA_SEQ.relative_to(ROOT),
        "| presente:",
        RNA_SEQ.is_file()
    )

    expr_headers = columns(
        RNA_SEQ
    )

    print(
        "[COLONNE RNA-SEQ]",
        (
            expr_headers
            if expr_headers
            else "non disponibili"
        )
    )

    gene_col, possible = (
        identify_gene_column(
            expr_headers
        )
    )

    print(
        "[JOIN GENE]",
        (
            gene_col
            if gene_col
            else f"non univoco {possible}"
        )
    )

    print(
        "[SOGLIA] full_score <=",
        CUTOFF,
        "(diagnostica; NON psRNATarget)"
    )

    print(
        "[OK] Audit soltanto: nessuna modifica"
    )

    return report_path, tsv, report


# ==========================================================
# 9. ELABORAZIONE PRINCIPALE
# ==========================================================

def run():

    report_path, tsv, src = audit()

    genes, annotation_info = (
        reference_data()
    )

    all_count = Counter()

    strong_count = Counter()

    seen = set()

    rows = []

    cols = columns(tsv)

    # ------------------------------------------------------
    # LETTURA E CONTROLLO DI TUTTI I SITI STEP20
    # ------------------------------------------------------

    with tsv.open(
        encoding="utf-8-sig",
        newline=""
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for row in reader:

            ident = (
                row.get("srna_id") or ""
            ).strip()

            if ident not in src["results"]:

                raise ValueError(
                    "RNA non presente nel report "
                    "STEP20: " + ident
                )

            tx = (
                row.get("transcript_id") or ""
            ).strip()

            gene = (
                row.get("gene_id") or ""
            ).strip()

            start = int(
                row["start_1based"]
            )

            stop = int(
                row["end_1based"]
            )

            length = int(
                row["srna_length"]
            )

            full = number(
                row["full_score"],
                "full_score"
            )

            hard = int(
                row["full_hard"]
            )

            gu = int(
                row["full_gu"]
            )

            hsp_score = number(
                row["hsp_score"],
                "hsp_score"
            )

            key = (
                ident,
                tx,
                start,
                stop
            )

            if key in seen:

                raise RuntimeError(
                    f"Sito duplicato: {key}"
                )

            seen.add(key)

            if stop - start + 1 != length:

                raise ValueError(
                    "Coordinate/lunghezza "
                    "discordanti: " + str(key)
                )

            expected_score = (
                Decimal(hard)
                + Decimal("0.5") * Decimal(gu)
            )

            if full != expected_score:

                raise ValueError(
                    "Scoring incoerente: "
                    + str(key)
                )

            if hsp_score > CUTOFF:

                raise ValueError(
                    "STEP20 ha incluso HSP >5: "
                    + str(key)
                )

            if (
                len(
                    row["full_site_DNA_5to3"]
                )
                != length
            ):

                raise ValueError(
                    "Sito nucleotidico "
                    "di lunghezza errata"
                )

            all_count[ident] += 1

            # ----------------------------------------------
            # SELEZIONE SOLO DEI SITI PRIORITARI
            # ----------------------------------------------

            if full > CUTOFF:
                continue

            strong_count[ident] += 1

            exp = (
                genes.get(gene)
                if gene
                else None
            )

            if not gene:

                mapping_status = (
                    "gene_id_absent"
                )

            elif not genes:

                mapping_status = (
                    "RNAseq_unavailable"
                )

            elif exp:

                mapping_status = (
                    "exact_gene_id_found"
                )

            else:

                mapping_status = (
                    "exact_gene_id_not_found"
                )

            annotated = dict(row)

            annotated[
                "rnaseq_join_status"
            ] = mapping_status

            annotated.update(
                expression_evidence(
                    gene,
                    exp
                )
            )

            rows.append(annotated)

    # ======================================================
    # 10. COERENZA CON IL REPORT STEP20
    # ======================================================

    for ident, expected in src[
        "results"
    ].items():

        if (
            all_count[ident]
            != int(
                expected[
                    "full_length_evaluable"
                ]
            )
        ):

            raise RuntimeError(
                "Numero siti diverso STEP20 "
                f"per {ident}: "
                f"{all_count[ident]}"
            )

        if (
            strong_count[ident]
            != int(
                expected["full_score_le_5"]
            )
        ):

            raise RuntimeError(
                "Numero <=5 diverso STEP20 "
                f"per {ident}: "
                f"{strong_count[ident]}"
            )

    rows.sort(
        key=lambda r: (

            number(
                r["full_score"],
                "full_score"
            ),

            int(r["full_hard"]),

            int(r["full_gu"]),

            r["srna_id"],

            r["transcript_id"],

            int(r["start_1based"])
        )
    )

    # ======================================================
    # 11. NUOVA CARTELLA OUTPUT
    # ======================================================

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP21_TRIAGE_TARGET_SRNA_"
        + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    # ======================================================
    # 12. TSV DETTAGLIATO
    # ======================================================

    detailed = (
        folder
        / "SITI_FULL_SCORE_LE5_ANNOTATI.tsv"
    )

    fields = (
        cols
        + ["rnaseq_join_status"]
        + [
            f"{prefix}_{field}"

            for prefix in (
                "B73",
                "H21"
            )

            for field in (
                "fdr",
                "log2fc",
                "DE_FDR_lt_005"
            )
        ]
    )

    with detailed.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()

        writer.writerows(rows)

    # ======================================================
    # 13. RAGGRUPPAMENTO PER GENE
    # ======================================================

    grouped = defaultdict(list)

    for row in rows:

        gene_or_tx = (
            row["gene_id"]
            or "NO_GENE::"
            + row["transcript_id"]
        )

        grouped[
            (
                row["srna_id"],
                gene_or_tx
            )
        ].append(row)

    summary_fields = [

        "srna_id",
        "gene_or_transcript",
        "gene_id",
        "site_count",
        "distinct_transcripts",
        "best_full_score",
        "min_full_hard",
        "min_full_gu",
        "rnaseq_join_status",
        "B73_fdr",
        "B73_log2fc",
        "B73_DE_FDR_lt_005",
        "H21_fdr",
        "H21_log2fc",
        "H21_DE_FDR_lt_005"
    ]

    gene_rows = []

    for (
        ident,
        gene_or_tx
    ), data in sorted(
        grouped.items()
    ):

        best = min(
            data,
            key=lambda r: (

                number(
                    r["full_score"],
                    "score"
                ),

                int(r["full_hard"]),

                int(r["full_gu"])
            )
        )

        gene_rows.append({

            "srna_id": ident,

            "gene_or_transcript":
                gene_or_tx,

            "gene_id": best["gene_id"],

            "site_count": len(data),

            "distinct_transcripts": len({
                x["transcript_id"]
                for x in data
            }),

            "best_full_score":
                best["full_score"],

            "min_full_hard": min(
                int(x["full_hard"])
                for x in data
            ),

            "min_full_gu": min(
                int(x["full_gu"])
                for x in data
            ),

            **{
                name: best.get(name, "")

                for name in summary_fields[8:]
            }
        })

    gene_rows.sort(
        key=lambda r: (

            number(
                r["best_full_score"],
                "best_full_score"
            ),

            r["srna_id"],

            r["gene_or_transcript"]
        )
    )

    gene_file = (
        folder
        / "PRIORITA_PER_GENE_O_TRASCRITTO.tsv"
    )

    with gene_file.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=summary_fields,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()

        writer.writerows(
            gene_rows
        )

    # ======================================================
    # 14. REPORT SCIENTIFICO
    # ======================================================

    report = {

        "step": "21",

        "status":
            "DESCRIPTIVE_TRIAGE_ONLY",

        "source_step20_report": str(
            report_path.relative_to(ROOT)
        ),

        "source_step20_tsv": str(
            tsv.relative_to(ROOT)
        ),

        "sha256_step20_tsv": digest(tsv),

        "rnaseq": annotation_info,

        "cutoff": str(CUTOFF),

        "score_is_psRNATarget_expectation":
            False,

        "all_full_length_sites_validated":
            dict(all_count),

        "candidate_sites_by_srna":
            dict(strong_count),

        "candidate_site_count":
            len(rows),

        "unique_srna_gene_or_transcript_pairs":
            len(gene_rows),

        "unique_annotated_genes": len({

            r["gene_id"]
            for r in rows
            if r["gene_id"]
        }),

        "rnaseq_exact_match_sites": sum(

            r["rnaseq_join_status"]
            == "exact_gene_id_found"

            for r in rows
        ),

        "b73_DE_support_sites": sum(

            r["B73_DE_FDR_lt_005"]
            == "yes"

            for r in rows
        ),

        "h21_DE_support_sites": sum(

            r["H21_DE_FDR_lt_005"]
            == "yes"

            for r in rows
        ),

        "outputs": [

            str(
                detailed.relative_to(ROOT)
            ),

            str(
                gene_file.relative_to(ROOT)
            )
        ],

        "limitations": [

            "L'espressione differenziale "
            "RNA-seq e' un contesto, "
            "non prova di targeting",

            "Manca validazione psRNATarget "
            "e degradome/PARE",

            "Selezione condizionata al "
            "prefiltro HSP20 senza gap",

            "Punteggio diagnostico locale, "
            "non calibrato statisticamente"
        ],

        "v32_modified": False,

        "v33_gui_modified": False,

        "docker_modified": False,

        "training_performed": False
    }

    report_file = (
        folder / "RAPPORTO_STEP21.json"
    )

    report_file.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    # ======================================================
    # 15. OUTPUT TERMINALE
    # ======================================================

    for ident in src["results"]:

        print(
            f"[RISULTATO] {ident}: "
            f"full={all_count[ident]}, "
            f"prioritari <=5="
            f"{strong_count[ident]}"
        )

    print(
        "[ANNOTAZIONE] RNA-seq:",
        annotation_info["status"],
        "| corrispondenze esatte siti=",
        report["rnaseq_exact_match_sites"]
    )

    print(
        "[RIEPILOGO]",
        len(rows),
        "siti,",
        len(gene_rows),
        "coppie sRNA/gene o trascritto,",
        report["unique_annotated_genes"],
        "geni annotati distinti"
    )

    print(
        "[REPORT]",
        report_file
    )

    print(
        "[TSV SITI]",
        detailed
    )

    print(
        "[TSV GENI]",
        gene_file
    )

    print(
        "[OK] Nessuna modifica "
        "a GUI, dati, modelli o Docker"
    )


# ==========================================================
# 16. PROGRAMMA PRINCIPALE
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

    print("=" * 65)

    print(
        "AI PRIMETTA V3.3 | STEP 21 - TRIAGE TARGET"
    )

    print("=" * 65)

    if args.audit:

        audit()

    else:

        run()


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
