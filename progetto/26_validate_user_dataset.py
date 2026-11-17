#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime
from pathlib import Path

import pandas as pd


BASE_DIR = Path(".").resolve()
OUT_ROOT = BASE_DIR / "outputs/user_data"


ID_COLUMNS = [
    "rnacentral_id",
    "ncrna_id",
    "transcript_id",
    "gene_id",
    "id",
    "sequence_id",
]

DE_COLUMNS = [
    "log2FC",
    "log2fc",
    "logFC",
    "pvalue",
    "p_value",
    "padj",
    "FDR",
    "fdr",
    "condition",
    "comparison",
]


def safe_name(name: str) -> str:
    name = name.strip().replace(" ", "_")
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    return name or "user_dataset"


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")

    if suffix == ".csv":
        return pd.read_csv(path)

    return pd.read_csv(path, sep=None, engine="python")


def inspect_fasta(path: Path) -> dict:
    n_records = 0
    lengths = []
    current_len = 0
    current_id = None
    first_ids = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    lengths.append(current_len)

                current_id = line[1:].split()[0]
                if len(first_ids) < 10:
                    first_ids.append(current_id)

                n_records += 1
                current_len = 0
            else:
                current_len += len(line)

    if current_id is not None:
        lengths.append(current_len)

    if lengths:
        min_len = min(lengths)
        max_len = max(lengths)
        mean_len = sum(lengths) / len(lengths)
    else:
        min_len = max_len = mean_len = 0

    return {
        "file": str(path),
        "type": "fasta",
        "n_records": n_records,
        "min_length": min_len,
        "max_length": max_len,
        "mean_length": mean_len,
        "first_ids": first_ids,
    }


def inspect_table(path: Path) -> dict:
    info = {
        "file": str(path),
        "type": "table",
        "status": "ok",
        "messages": [],
    }

    try:
        df = read_table(path)
    except Exception as exc:
        info["status"] = "error"
        info["messages"].append(f"Errore lettura tabella: {exc}")
        return info

    columns = list(df.columns)
    lower_map = {str(c).lower(): c for c in columns}

    found_id_cols = [c for c in ID_COLUMNS if c.lower() in lower_map]
    found_de_cols = [c for c in DE_COLUMNS if c.lower() in lower_map]

    info.update(
        {
            "n_rows": int(len(df)),
            "n_columns": int(len(columns)),
            "columns": columns,
            "id_columns_found": found_id_cols,
            "differential_expression_columns_found": found_de_cols,
            "preview_records": df.head(10).fillna("").to_dict(orient="records"),
        }
    )

    if not found_id_cols:
        info["status"] = "warning"
        info["messages"].append(
            "Nessuna colonna identificativa standard trovata. "
            "Consigliate: rnacentral_id, ncrna_id, transcript_id, gene_id, id."
        )

    if not found_de_cols:
        info["messages"].append(
            "Nessuna colonna tipica di espressione differenziale trovata. "
            "Consigliate: log2FC, pvalue, FDR/padj, condition/comparison."
        )

    return info


def validate_dataset(input_dir: Path, dataset_name: str) -> dict:
    dataset_name = safe_name(dataset_name)
    report_dir = OUT_ROOT / dataset_name
    report_dir.mkdir(parents=True, exist_ok=True)

    files = sorted([p for p in input_dir.iterdir() if p.is_file()])

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "dataset_name": dataset_name,
        "input_dir": str(input_dir),
        "status": "ok",
        "files_total": len(files),
        "files": [],
        "summary": {
            "tables": 0,
            "fasta": 0,
            "other": 0,
            "errors": 0,
            "warnings": 0,
        },
        "interpretation": (
            "Validazione tecnica del dataset utente. "
            "Questo step non esegue ancora analisi RNA-seq completa, "
            "ma prepara l'importazione nel sistema AI/RAG."
        ),
    }

    if not files:
        report["status"] = "error"
        report["summary"]["errors"] += 1
        report["files"].append(
            {
                "type": "none",
                "status": "error",
                "messages": ["Nessun file trovato nella cartella input."],
            }
        )

    for path in files:
        suffix = path.suffix.lower()

        if suffix in [".fa", ".fasta", ".fna"]:
            info = inspect_fasta(path)
            report["summary"]["fasta"] += 1

        elif suffix in [".csv", ".tsv", ".txt"]:
            info = inspect_table(path)
            report["summary"]["tables"] += 1

        else:
            info = {
                "file": str(path),
                "type": "other",
                "status": "warning",
                "messages": ["Formato non analizzato automaticamente."],
            }
            report["summary"]["other"] += 1

        if info.get("status") == "error":
            report["summary"]["errors"] += 1
        elif info.get("status") == "warning":
            report["summary"]["warnings"] += 1

        report["files"].append(info)

    if report["summary"]["errors"] > 0:
        report["status"] = "error"
    elif report["summary"]["warnings"] > 0:
        report["status"] = "warning"

    report_json = report_dir / "user_dataset_validation_report.json"
    summary_txt = report_dir / "user_dataset_validation_summary.txt"

    report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = [
        "USER DATASET VALIDATION SUMMARY",
        "=" * 80,
        f"Dataset: {dataset_name}",
        f"Input dir: {input_dir}",
        f"Status: {report['status']}",
        f"Files total: {report['files_total']}",
        "",
        f"Tables: {report['summary']['tables']}",
        f"FASTA: {report['summary']['fasta']}",
        f"Other: {report['summary']['other']}",
        f"Warnings: {report['summary']['warnings']}",
        f"Errors: {report['summary']['errors']}",
        "",
        "Files:",
    ]

    for item in report["files"]:
        lines.append(f"- {item.get('file', '-')}")
        lines.append(f"  type: {item.get('type')}")
        lines.append(f"  status: {item.get('status', 'ok')}")
        for msg in item.get("messages", []):
            lines.append(f"  note: {msg}")

    lines.extend(
        [
            "",
            "Interpretation:",
            report["interpretation"],
        ]
    )

    summary_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    print("")
    print(f"Report JSON: {report_json}")
    print(f"Summary TXT: {summary_txt}")

    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--dataset-name", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input_dir).resolve()
    if not input_dir.exists():
        raise FileNotFoundError(f"Cartella input non trovata: {input_dir}")

    report = validate_dataset(input_dir, args.dataset_name)

    return 0 if report["status"] in ["ok", "warning"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
