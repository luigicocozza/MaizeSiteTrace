#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


BASE_DIR = Path(".").resolve()

UPLOAD_ROOT = BASE_DIR / "data/user_data/uploads"
OUT_ROOT = BASE_DIR / "outputs/user_data"
RANKING_CSV = BASE_DIR / "data/processed/ranking/ncrna_exploratory_ranking.csv"


ID_CANDIDATES = [
    "rnacentral_id",
    "ncrna_id",
    "transcript_id",
    "sequence_id",
    "id",
]

LOGFC_CANDIDATES = [
    "log2FC",
    "log2fc",
    "logFC",
    "logfc",
]

FDR_CANDIDATES = [
    "FDR",
    "fdr",
    "padj",
    "adjusted_pvalue",
    "adjusted_p_value",
]

PVALUE_CANDIDATES = [
    "pvalue",
    "p_value",
    "PValue",
    "P.Value",
]

CONDITION_CANDIDATES = [
    "condition",
    "comparison",
    "contrast",
]


def find_latest_upload_dir() -> Path:
    if not UPLOAD_ROOT.exists():
        raise FileNotFoundError(f"Cartella upload non trovata: {UPLOAD_ROOT}")

    dirs = [p for p in UPLOAD_ROOT.iterdir() if p.is_dir()]

    if not dirs:
        raise FileNotFoundError(f"Nessun dataset caricato in: {UPLOAD_ROOT}")

    dirs = sorted(dirs, key=lambda p: p.stat().st_mtime, reverse=True)
    return dirs[0]


def find_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower_map = {str(c).lower(): c for c in df.columns}

    for candidate in candidates:
        key = candidate.lower()
        if key in lower_map:
            return lower_map[key]

    return None


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix == ".csv":
        return pd.read_csv(path)

    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")

    raise ValueError(f"Formato tabella non supportato: {path}")


def normalize_id(value) -> str:
    if pd.isna(value):
        return ""

    return str(value).strip()


def collect_user_tables(dataset_dir: Path) -> list[dict]:
    tables = []

    for path in sorted(dataset_dir.iterdir()):
        if not path.is_file():
            continue

        if path.suffix.lower() not in [".csv", ".tsv"]:
            continue

        try:
            df = read_table(path)
        except Exception as exc:
            tables.append(
                {
                    "file": str(path),
                    "status": "error",
                    "message": str(exc),
                    "df": None,
                    "id_column": None,
                }
            )
            continue

        id_col = find_column(df, ID_CANDIDATES)

        tables.append(
            {
                "file": str(path),
                "status": "ok" if id_col else "warning",
                "message": "ok" if id_col else "nessuna colonna ID compatibile trovata",
                "df": df,
                "id_column": id_col,
            }
        )

    return tables


def build_user_evidence_table(tables: list[dict]) -> pd.DataFrame:
    records = []

    for item in tables:
        df = item.get("df")
        id_col = item.get("id_column")

        if df is None or not id_col:
            continue

        logfc_col = find_column(df, LOGFC_CANDIDATES)
        fdr_col = find_column(df, FDR_CANDIDATES)
        pvalue_col = find_column(df, PVALUE_CANDIDATES)
        condition_col = find_column(df, CONDITION_CANDIDATES)

        for _, row in df.iterrows():
            rid = normalize_id(row.get(id_col))

            if not rid:
                continue

            record = {
                "rnacentral_id": rid,
                "source_file": Path(item["file"]).name,
                "user_id_column": id_col,
                "user_log2FC": row.get(logfc_col) if logfc_col else None,
                "user_FDR": row.get(fdr_col) if fdr_col else None,
                "user_pvalue": row.get(pvalue_col) if pvalue_col else None,
                "user_condition": row.get(condition_col) if condition_col else None,
            }

            records.append(record)

    if not records:
        return pd.DataFrame(
            columns=[
                "rnacentral_id",
                "source_file",
                "user_id_column",
                "user_log2FC",
                "user_FDR",
                "user_pvalue",
                "user_condition",
            ]
        )

    return pd.DataFrame(records)


def numeric_min(series: pd.Series):
    values = pd.to_numeric(series, errors="coerce")
    values = values.dropna()

    if values.empty:
        return None

    return float(values.min())


def numeric_max_abs(series: pd.Series):
    values = pd.to_numeric(series, errors="coerce")
    values = values.dropna()

    if values.empty:
        return None

    idx = values.abs().idxmax()
    return float(values.loc[idx])


def summarize_user_evidence(user_df: pd.DataFrame) -> pd.DataFrame:
    if user_df.empty:
        return pd.DataFrame(
            columns=[
                "rnacentral_id",
                "user_data_present",
                "user_records_count",
                "user_best_log2FC_abs",
                "user_min_FDR",
                "user_min_pvalue",
                "user_conditions",
                "user_source_files",
            ]
        )

    grouped = []

    for rid, sub in user_df.groupby("rnacentral_id"):
        conditions = sorted(
            {
                str(x)
                for x in sub["user_condition"].dropna().tolist()
                if str(x).strip()
            }
        )

        source_files = sorted(
            {
                str(x)
                for x in sub["source_file"].dropna().tolist()
                if str(x).strip()
            }
        )

        grouped.append(
            {
                "rnacentral_id": rid,
                "user_data_present": True,
                "user_records_count": int(len(sub)),
                "user_best_log2FC_abs": numeric_max_abs(sub["user_log2FC"]),
                "user_min_FDR": numeric_min(sub["user_FDR"]),
                "user_min_pvalue": numeric_min(sub["user_pvalue"]),
                "user_conditions": "; ".join(conditions),
                "user_source_files": "; ".join(source_files),
            }
        )

    return pd.DataFrame(grouped)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset-dir",
        default=None,
        help="Cartella dataset utente. Se omessa, usa l'ultimo upload.",
    )
    parser.add_argument(
        "--ranking-csv",
        default=str(RANKING_CSV),
        help="File ranking ncRNA.",
    )
    args = parser.parse_args()

    if args.dataset_dir:
        dataset_dir = Path(args.dataset_dir).resolve()
    else:
        dataset_dir = find_latest_upload_dir().resolve()

    if not dataset_dir.exists():
        raise FileNotFoundError(f"Dataset dir non trovata: {dataset_dir}")

    ranking_path = Path(args.ranking_csv).resolve()

    if not ranking_path.exists():
        raise FileNotFoundError(f"Ranking CSV non trovato: {ranking_path}")

    ranking = pd.read_csv(ranking_path)

    ranking_id_col = find_column(ranking, ID_CANDIDATES)

    if not ranking_id_col:
        raise ValueError(
            "Nessuna colonna ID trovata nel ranking. "
            "Serve una colonna tipo rnacentral_id, ncrna_id, transcript_id o id."
        )

    ranking = ranking.copy()
    ranking["rnacentral_id_normalized"] = ranking[ranking_id_col].map(normalize_id)

    tables = collect_user_tables(dataset_dir)
    user_raw = build_user_evidence_table(tables)
    user_summary = summarize_user_evidence(user_raw)

    if not user_summary.empty:
        user_summary["rnacentral_id_normalized"] = user_summary["rnacentral_id"].map(normalize_id)
    else:
        user_summary["rnacentral_id_normalized"] = []

    merged = ranking.merge(
        user_summary,
        on="rnacentral_id_normalized",
        how="left",
        suffixes=("", "_user"),
    )

    merged["user_data_present"] = merged["user_data_present"].fillna(False)
    merged["user_records_count"] = merged["user_records_count"].fillna(0).astype(int)

    if "rnacentral_id_user" in merged.columns:
        merged = merged.drop(columns=["rnacentral_id_user"])

    out_dir = OUT_ROOT / dataset_dir.name
    out_dir.mkdir(parents=True, exist_ok=True)

    ranking_out = out_dir / "ranking_with_user_evidence.csv"
    user_raw_out = out_dir / "user_evidence_raw_records.csv"
    matched_out = out_dir / "matched_user_ncrna.csv"
    unmatched_out = out_dir / "unmatched_user_ncrna.csv"
    report_json = out_dir / "user_ranking_integration_report.json"
    summary_txt = out_dir / "user_ranking_integration_summary.txt"

    merged.to_csv(ranking_out, index=False)
    user_raw.to_csv(user_raw_out, index=False)

    matched = merged[merged["user_data_present"] == True].copy()
    matched.to_csv(matched_out, index=False)

    ranking_ids = set(ranking["rnacentral_id_normalized"].dropna().astype(str))
    user_ids = set(user_raw["rnacentral_id"].dropna().astype(str)) if not user_raw.empty else set()
    unmatched_ids = sorted(user_ids - ranking_ids)

    unmatched = pd.DataFrame({"rnacentral_id": unmatched_ids})
    unmatched.to_csv(unmatched_out, index=False)

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "dataset_dir": str(dataset_dir),
        "ranking_csv": str(ranking_path),
        "ranking_id_column": ranking_id_col,
        "ranking_rows": int(len(ranking)),
        "user_raw_records": int(len(user_raw)),
        "user_unique_ids": int(len(user_ids)),
        "matched_ranking_rows": int(len(matched)),
        "unmatched_user_ids": int(len(unmatched_ids)),
        "tables": [
            {
                "file": item["file"],
                "status": item["status"],
                "message": item["message"],
                "id_column": item["id_column"],
            }
            for item in tables
        ],
        "outputs": {
            "ranking_with_user_evidence": str(ranking_out),
            "user_evidence_raw_records": str(user_raw_out),
            "matched_user_ncrna": str(matched_out),
            "unmatched_user_ncrna": str(unmatched_out),
            "summary": str(summary_txt),
            "report_json": str(report_json),
        },
        "methodological_note": (
            "Il collegamento e' tecnico-esplorativo. "
            "La presenza di un ncRNA nei dati del ricercatore non dimostra da sola "
            "un ruolo biologico nella risposta alla siccita'."
        ),
    }

    report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = [
        "USER DATASET - RANKING INTEGRATION SUMMARY",
        "=" * 80,
        f"Dataset dir: {dataset_dir}",
        f"Ranking CSV: {ranking_path}",
        f"Ranking ID column: {ranking_id_col}",
        "",
        f"Ranking rows: {len(ranking)}",
        f"User raw records: {len(user_raw)}",
        f"User unique IDs: {len(user_ids)}",
        f"Matched ranking rows: {len(matched)}",
        f"Unmatched user IDs: {len(unmatched_ids)}",
        "",
        "Input tables:",
    ]

    for item in tables:
        lines.append(f"- {item['file']}")
        lines.append(f"  status: {item['status']}")
        lines.append(f"  id_column: {item['id_column']}")
        lines.append(f"  message: {item['message']}")

    lines.extend(
        [
            "",
            "Outputs:",
            f"- {ranking_out}",
            f"- {user_raw_out}",
            f"- {matched_out}",
            f"- {unmatched_out}",
            f"- {report_json}",
            "",
            "Methodological note:",
            report["methodological_note"],
        ]
    )

    summary_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
