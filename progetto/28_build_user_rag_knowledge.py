#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd


BASE_DIR = Path(".").resolve()
OUT_ROOT = BASE_DIR / "outputs/user_data"


def latest_output_dir() -> Path:
    if not OUT_ROOT.exists():
        raise FileNotFoundError(f"Cartella non trovata: {OUT_ROOT}")

    dirs = [p for p in OUT_ROOT.iterdir() if p.is_dir()]

    if not dirs:
        raise FileNotFoundError("Nessun output user_data trovato.")

    dirs = sorted(dirs, key=lambda p: p.stat().st_mtime, reverse=True)
    return dirs[0]


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def add_chunk(chunks: list[dict], title: str, content: str, source: str, chunk_type: str) -> None:
    content = str(content).strip()

    if not content:
        return

    chunks.append(
        {
            "chunk_id": f"user_chunk_{len(chunks) + 1:05d}",
            "title": title,
            "content": content,
            "source": source,
            "chunk_type": chunk_type,
        }
    )


def main() -> int:
    out_dir = latest_output_dir()

    matched_path = out_dir / "matched_user_ncrna.csv"
    ranking_path = out_dir / "ranking_with_user_evidence.csv"
    raw_path = out_dir / "user_evidence_raw_records.csv"
    summary_path = out_dir / "user_ranking_integration_summary.txt"

    matched = read_csv(matched_path)
    ranking = read_csv(ranking_path)
    raw = read_csv(raw_path)
    summary = read_text(summary_path)

    chunks: list[dict] = []

    add_chunk(
        chunks,
        "Riepilogo integrazione dati ricercatore",
        summary,
        str(summary_path),
        "summary",
    )

    if not matched.empty:
        for _, row in matched.iterrows():
            rid = row.get("rnacentral_id", "")
            rank = row.get("exploratory_rank", "")
            rna_type = row.get("rna_type", "")
            score = row.get("exploratory_rank_score", "")
            logfc = row.get("user_best_log2FC_abs", "")
            fdr = row.get("user_min_FDR", "")
            pvalue = row.get("user_min_pvalue", "")
            condition = row.get("user_conditions", "")
            source_file = row.get("user_source_files", "")
            note = row.get("methodological_note", "")

            content = (
                f"ncRNA {rid}. "
                f"Rank esplorativo: {rank}. "
                f"Tipo RNA: {rna_type}. "
                f"Score ranking: {score}. "
                f"Presente nei dati del ricercatore: si. "
                f"log2FC massimo in valore assoluto: {logfc}. "
                f"FDR minimo: {fdr}. "
                f"p-value minimo: {pvalue}. "
                f"Condizione/confronto: {condition}. "
                f"File sorgente: {source_file}. "
                f"Nota metodologica: {note}. "
                f"Interpretazione: la coincidenza tra ranking e dati utente e' evidenza tecnica, "
                f"non prova biologica definitiva di risposta alla siccita'."
            )

            add_chunk(
                chunks,
                f"Candidato supportato da dati utente: {rid}",
                content,
                str(matched_path),
                "matched_candidate",
            )

    if not raw.empty:
        for _, row in raw.iterrows():
            rid = row.get("rnacentral_id", "")
            logfc = row.get("user_log2FC", "")
            fdr = row.get("user_FDR", "")
            pvalue = row.get("user_pvalue", "")
            condition = row.get("user_condition", "")
            source_file = row.get("source_file", "")

            content = (
                f"Record grezzo utente per ncRNA {rid}. "
                f"log2FC: {logfc}. "
                f"FDR: {fdr}. "
                f"p-value: {pvalue}. "
                f"Condizione: {condition}. "
                f"File sorgente: {source_file}."
            )

            add_chunk(
                chunks,
                f"Record utente: {rid}",
                content,
                str(raw_path),
                "raw_user_record",
            )

    if not ranking.empty:
        n_total = len(ranking)
        n_supported = int(ranking.get("user_data_present", pd.Series(dtype=bool)).fillna(False).sum())

        content = (
            f"Il ranking integrato contiene {n_total} ncRNA. "
            f"Gli ncRNA con evidenza nei dati caricati dal ricercatore sono {n_supported}. "
            f"Questa integrazione serve a priorizzare candidati per successive verifiche biologiche."
        )

        add_chunk(
            chunks,
            "Statistiche ranking integrato con dati utente",
            content,
            str(ranking_path),
            "ranking_statistics",
        )

    kb_jsonl = out_dir / "user_rag_knowledge_base.jsonl"
    kb_json = out_dir / "user_rag_knowledge_base.json"
    report_json = out_dir / "user_rag_knowledge_report.json"
    summary_txt = out_dir / "user_rag_knowledge_summary.txt"

    with kb_jsonl.open("w", encoding="utf-8") as f:
        for chunk in chunks:
            f.write(json.dumps(chunk, ensure_ascii=False) + "\n")

    kb_json.write_text(json.dumps(chunks, indent=2, ensure_ascii=False), encoding="utf-8")

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "output_dir": str(out_dir),
        "chunks": len(chunks),
        "matched_candidates": int(len(matched)),
        "raw_records": int(len(raw)),
        "outputs": {
            "jsonl": str(kb_jsonl),
            "json": str(kb_json),
            "report": str(report_json),
            "summary": str(summary_txt),
        },
        "methodological_note": (
            "La knowledge base RAG utente contiene dati caricati dal ricercatore "
            "e risultati di integrazione con il ranking. Le risposte sono esplorative."
        ),
    }

    report_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = [
        "USER RAG KNOWLEDGE BASE SUMMARY",
        "=" * 80,
        f"Output dir: {out_dir}",
        f"Chunks created: {len(chunks)}",
        f"Matched candidates: {len(matched)}",
        f"Raw user records: {len(raw)}",
        "",
        "Outputs:",
        f"- {kb_jsonl}",
        f"- {kb_json}",
        f"- {report_json}",
        "",
        "Methodological note:",
        report["methodological_note"],
    ]

    summary_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
