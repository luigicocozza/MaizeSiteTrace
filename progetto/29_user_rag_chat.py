#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
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


def tokenize(text: str) -> set[str]:
    tokens = re.findall(r"[A-Za-z0-9_:.+-]+", text.lower())
    return {t for t in tokens if len(t) >= 3}


def load_chunks(path: Path) -> list[dict]:
    chunks = []

    if not path.exists():
        return chunks

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                chunks.append(json.loads(line))

    return chunks


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def score_chunk(question_tokens: set[str], chunk: dict) -> int:
    content = f"{chunk.get('title', '')} {chunk.get('content', '')}"
    chunk_tokens = tokenize(content)
    return len(question_tokens & chunk_tokens)


def extract_fdr_threshold(question: str) -> float | None:
    q = question.lower()

    if "fdr" not in q:
        return None

    patterns = [
        r"fdr\s*(?:<=|<|inferiore a|minore di|sotto)\s*([0-9]+(?:[.,][0-9]+)?)",
        r"([0-9]+(?:[.,][0-9]+)?)\s*(?:di fdr|fdr)",
    ]

    for pattern in patterns:
        m = re.search(pattern, q)
        if m:
            value = m.group(1).replace(",", ".")
            try:
                return float(value)
            except Exception:
                return None

    return 0.05


def deterministic_user_answer(question: str, out_dir: Path) -> str | None:
    q = question.lower()

    matched_path = out_dir / "matched_user_ncrna.csv"
    raw_path = out_dir / "user_evidence_raw_records.csv"
    report_path = out_dir / "user_ranking_integration_report.json"

    matched = read_csv(matched_path)
    raw = read_csv(raw_path)

    if matched.empty and raw.empty:
        return None

    wants_fdr = "fdr" in q
    wants_ranking = "ranking" in q or "rank" in q
    wants_logfc = "log2fc" in q or "logfc" in q
    wants_present = (
        "presenti" in q
        or "presente" in q
        or "caricati" in q
        or "ricercatore" in q
        or "dataset" in q
    )

    if not any([wants_fdr, wants_ranking, wants_logfc, wants_present]):
        return None

    df = matched.copy()

    if df.empty:
        return (
            "Non risultano ncRNA del ranking presenti nei dati caricati dal ricercatore.\n\n"
            "Interpretazione: l'assenza di match puo' dipendere da ID non corrispondenti, "
            "formato diverso degli identificativi o dataset non ancora integrato."
        )

    threshold = extract_fdr_threshold(question)

    if threshold is not None and "user_min_FDR" in df.columns:
        fdr_values = pd.to_numeric(df["user_min_FDR"], errors="coerce")
        df = df[fdr_values <= threshold].copy()

    if df.empty:
        return (
            f"Nessun ncRNA del ranking integrato rispetta il filtro richiesto.\n\n"
            f"Soglia FDR considerata: {threshold}\n\n"
            "Interpretazione: il filtro e' tecnico-esplorativo e dipende dai dati caricati."
        )

    preferred_cols = [
        "exploratory_rank",
        "rnacentral_id",
        "rna_type",
        "exploratory_rank_score",
        "user_best_log2FC_abs",
        "user_min_FDR",
        "user_min_pvalue",
        "user_conditions",
        "user_source_files",
    ]

    existing_cols = [c for c in preferred_cols if c in df.columns]
    df_small = df[existing_cols].copy()

    lines = []
    lines.append("Risposta basata sui dati caricati dal ricercatore e sul ranking integrato.")
    lines.append("")
    lines.append(f"Domanda: {question}")
    lines.append("")

    if threshold is not None:
        lines.append(f"Filtro applicato: FDR <= {threshold}")
        lines.append("")

    lines.append(f"ncRNA trovati: {len(df_small)}")
    lines.append("")

    for _, row in df_small.iterrows():
        rid = row.get("rnacentral_id", "")
        rank = row.get("exploratory_rank", "")
        rna_type = row.get("rna_type", "")
        score = row.get("exploratory_rank_score", "")
        logfc = row.get("user_best_log2FC_abs", "")
        fdr = row.get("user_min_FDR", "")
        pvalue = row.get("user_min_pvalue", "")
        condition = row.get("user_conditions", "")
        source = row.get("user_source_files", "")

        lines.append(f"- {rid}")
        lines.append(f"  rank esplorativo: {rank}")
        lines.append(f"  tipo RNA: {rna_type}")
        lines.append(f"  score ranking: {score}")
        lines.append(f"  log2FC utente: {logfc}")
        lines.append(f"  FDR minimo: {fdr}")
        lines.append(f"  p-value minimo: {pvalue}")
        lines.append(f"  condizione: {condition}")
        lines.append(f"  file sorgente: {source}")

    lines.append("")
    lines.append("Interpretazione corretta:")
    lines.append(
        "Questi risultati indicano una coincidenza tecnica tra ncRNA del ranking "
        "e dati caricati dal ricercatore. Non costituiscono da soli una prova biologica "
        "definitiva di risposta alla siccita'."
    )

    if report_path.exists():
        lines.append("")
        lines.append(f"Report integrazione: {report_path}")

    return "\n".join(lines)


def rag_answer(question: str, out_dir: Path, top_k: int = 5) -> str:
    kb_path = out_dir / "user_rag_knowledge_base.jsonl"

    if not kb_path.exists():
        raise FileNotFoundError(
            f"Knowledge base utente non trovata: {kb_path}. "
            f"Eseguire prima 28_build_user_rag_knowledge.py"
        )

    chunks = load_chunks(kb_path)
    q_tokens = tokenize(question)

    scored = []

    for chunk in chunks:
        score = score_chunk(q_tokens, chunk)
        scored.append((score, chunk))

    scored = sorted(scored, key=lambda x: x[0], reverse=True)
    selected = [chunk for score, chunk in scored if score > 0][:top_k]

    if not selected:
        selected = [chunk for _, chunk in scored[:top_k]]

    context_lines = []

    for chunk in selected:
        context_lines.append(f"- {chunk.get('title')}: {chunk.get('content')}")

    context = "\n".join(context_lines)

    answer = (
        "Risposta basata sui dati caricati dal ricercatore e sul ranking integrato.\n\n"
        f"Domanda: {question}\n\n"
        "Contesto rilevante trovato:\n"
        f"{context}\n\n"
        "Interpretazione corretta:\n"
        "Questi risultati indicano una coincidenza tecnica tra ncRNA del ranking e dati caricati. "
        "Non costituiscono da soli una prova biologica definitiva di risposta alla siccita'."
    )

    return answer


def answer_question(question: str, top_k: int = 5) -> dict:
    out_dir = latest_output_dir()

    deterministic = deterministic_user_answer(question, out_dir)

    if deterministic is not None:
        answer = deterministic
        provider = "deterministic_user_data_filter"
    else:
        answer = rag_answer(question, out_dir, top_k=top_k)
        provider = "simple_user_rag"

    return {
        "question": question,
        "output_dir": str(out_dir),
        "provider": provider,
        "answer": answer,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--question", required=True)
    parser.add_argument("--top-k", type=int, default=5)
    args = parser.parse_args()

    result = answer_question(args.question, top_k=args.top_k)

    print(result["answer"])
    print("")
    print(f"Provider: {result['provider']}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
