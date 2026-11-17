#!/usr/bin/env python3
"""
FastAPI backend per maize_drought_ncrna_gnn.

Avvio:
./run_api.sh

Endpoint:
- GET  /health
- GET  /summary
- GET  /candidates
- GET  /candidate/{rnacentral_id}
- GET  /candidate/{rnacentral_id}/similar
- POST /rag/ask
- GET  /files/inventory
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional

import pandas as pd
from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel


BASE_DIR = Path(__file__).resolve().parents[1]

PATHS = {
    "ranking": BASE_DIR / "data/processed/ranking/ncrna_exploratory_ranking.csv",
    "top50": BASE_DIR / "data/processed/ranking/ncrna_exploratory_top50.csv",
    "similarity": BASE_DIR / "data/processed/similarity/ncrna_similarity_topk.csv",
    "features": BASE_DIR / "data/processed/ncrna_features.csv",
    "embeddings": BASE_DIR / "data/processed/embeddings/ncrna_embeddings.csv",
    "graph_nodes": BASE_DIR / "data/processed/graph_nodes.csv",
    "graph_edges": BASE_DIR / "data/processed/graph_edges.csv",
    "rag_kb": BASE_DIR / "data/processed/rag/project_knowledge_base.jsonl",
    "rag_engine": BASE_DIR / "16_rag_chat_engine.py",
    "method_report": BASE_DIR / "outputs/reports/methodological_report.md",
}


app = FastAPI(
    title="Maize ncRNA GNN API",
    description="Backend API per ranking, similarità, candidato e Chat RAG.",
    version="0.1.0",
)


class RagAskRequest(BaseModel):
    question: str
    provider: str = "extractive"
    top_k: int = 6


def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def df_records(df: pd.DataFrame) -> List[Dict[str, Any]]:
    if df is None or df.empty:
        return []
    safe = df.copy()
    safe = safe.where(pd.notnull(safe), None)
    return safe.to_dict(orient="records")


def load_rag_module():
    engine_path = PATHS["rag_engine"]

    if not engine_path.exists():
        raise RuntimeError("16_rag_chat_engine.py non trovato.")

    spec = importlib.util.spec_from_file_location("rag_chat_engine_api", engine_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Impossibile caricare il modulo RAG.")

    module = importlib.util.module_from_spec(spec)

    # Necessario per dataclass con Python 3.14
    sys.modules[spec.name] = module

    spec.loader.exec_module(module)
    return module


def get_rag_answer(question: str, provider: str = "extractive", top_k: int = 6) -> Dict[str, Any]:
    module = load_rag_module()

    if not PATHS["rag_kb"].exists():
        raise RuntimeError("Knowledge base RAG non trovata.")

    records = module.load_jsonl(PATHS["rag_kb"])
    retriever = module.RagRetriever(records)

    args = SimpleNamespace(
        kb=str(PATHS["rag_kb"]),
        question=question,
        interactive=False,
        provider=provider,
        top_k=top_k,
        ollama_url="http://localhost:11434/api/generate",
        ollama_model="llama3.1:8b",
        openai_model="gpt-4o-mini",
        max_context_chars=9000,
    )

    answer, chunks, used_provider = module.answer_question(question, retriever, args)

    sources = [
        {
            "rank": c.rank,
            "score": c.score,
            "chunk_id": c.chunk_id,
            "title": c.title,
            "source": c.source,
            "source_type": c.source_type,
        }
        for c in chunks
    ]

    return {
        "question": question,
        "provider": used_provider,
        "answer": answer,
        "sources": sources,
    }


@app.get("/health")
def health() -> Dict[str, Any]:
    return {
        "status": "ok",
        "project": "maize_drought_ncrna_gnn",
        "base_dir": str(BASE_DIR),
    }


@app.get("/files/inventory")
def files_inventory() -> Dict[str, Any]:
    rows = []

    for key, path in PATHS.items():
        rows.append(
            {
                "artifact": key,
                "path": str(path.relative_to(BASE_DIR)),
                "exists": path.exists(),
                "size_kb": round(path.stat().st_size / 1024, 3) if path.exists() else 0.0,
            }
        )

    return {
        "n_files": len(rows),
        "n_existing": sum(1 for r in rows if r["exists"]),
        "files": rows,
    }


@app.get("/summary")
def summary() -> Dict[str, Any]:
    ranking = load_csv(PATHS["ranking"])
    similarity = load_csv(PATHS["similarity"])
    features = load_csv(PATHS["features"])
    embeddings = load_csv(PATHS["embeddings"])
    nodes = load_csv(PATHS["graph_nodes"])
    edges = load_csv(PATHS["graph_edges"])

    rna_type_distribution = {}
    if not ranking.empty and "rna_type" in ranking.columns:
        rna_type_distribution = ranking["rna_type"].value_counts().to_dict()

    return {
        "project": "maize_drought_ncrna_gnn",
        "n_ranked_ncrna": len(ranking),
        "n_similarity_rows": len(similarity),
        "n_feature_rows": len(features),
        "n_embedding_rows": len(embeddings),
        "n_graph_nodes": len(nodes),
        "n_graph_edges": len(edges),
        "rna_type_distribution": rna_type_distribution,
        "methodological_note": (
            "Ranking, embedding e similarità sono esplorativi e non supervisionati. "
            "Non dimostrano validazione biologica della risposta alla siccità."
        ),
    }


@app.get("/candidates")
def candidates(
    query: Optional[str] = Query(default=None),
    rna_type: Optional[str] = Query(default=None),
    min_score: Optional[float] = Query(default=None),
    max_rank: Optional[int] = Query(default=None),
    limit: int = Query(default=50, ge=1, le=250),
) -> Dict[str, Any]:
    ranking = load_csv(PATHS["ranking"])

    if ranking.empty:
        return {"count": 0, "results": []}

    df = ranking.copy()

    if query:
        q = query.lower()
        mask = pd.Series(False, index=df.index)
        for col in ["rnacentral_id", "node_id", "rna_type", "top3_similar_ids"]:
            if col in df.columns:
                mask = mask | df[col].astype(str).str.lower().str.contains(q, na=False)
        df = df[mask]

    if rna_type and "rna_type" in df.columns:
        df = df[df["rna_type"].astype(str) == str(rna_type)]

    if min_score is not None and "exploratory_rank_score" in df.columns:
        score = pd.to_numeric(df["exploratory_rank_score"], errors="coerce")
        df = df[score >= min_score]

    if max_rank is not None and "exploratory_rank" in df.columns:
        rank = pd.to_numeric(df["exploratory_rank"], errors="coerce")
        df = df[rank <= max_rank]

    if "exploratory_rank" in df.columns:
        df = df.sort_values("exploratory_rank")

    df = df.head(limit)

    return {
        "count": len(df),
        "results": df_records(df),
        "methodological_note": "Candidati ordinati da ranking esplorativo non supervisionato.",
    }


@app.get("/candidate/{rnacentral_id}")
def candidate_detail(rnacentral_id: str) -> Dict[str, Any]:
    ranking = load_csv(PATHS["ranking"])
    features = load_csv(PATHS["features"])
    embeddings = load_csv(PATHS["embeddings"])

    if ranking.empty or "rnacentral_id" not in ranking.columns:
        raise HTTPException(status_code=404, detail="Ranking non disponibile.")

    row_df = ranking[ranking["rnacentral_id"].astype(str) == str(rnacentral_id)]

    if row_df.empty:
        raise HTTPException(status_code=404, detail=f"Candidato non trovato: {rnacentral_id}")

    row = row_df.iloc[0].to_dict()
    node_id = str(row.get("node_id", f"ncrna:{rnacentral_id}"))

    feat_df = pd.DataFrame()
    if not features.empty and "node_id" in features.columns:
        feat_df = features[features["node_id"].astype(str) == node_id]

    emb_df = pd.DataFrame()
    if not embeddings.empty and "node_id" in embeddings.columns:
        emb_df = embeddings[embeddings["node_id"].astype(str) == node_id]

    return {
        "candidate": row,
        "features": df_records(feat_df),
        "embedding": df_records(emb_df),
        "methodological_note": (
            "Dettaglio candidato esplorativo. Non è prova di risposta alla siccità."
        ),
    }


@app.get("/candidate/{rnacentral_id}/similar")
def candidate_similar(
    rnacentral_id: str,
    limit: int = Query(default=10, ge=1, le=50),
) -> Dict[str, Any]:
    ranking = load_csv(PATHS["ranking"])
    similarity = load_csv(PATHS["similarity"])

    if ranking.empty or similarity.empty:
        raise HTTPException(status_code=404, detail="Ranking o similarità non disponibili.")

    row_df = ranking[ranking["rnacentral_id"].astype(str) == str(rnacentral_id)]

    if row_df.empty:
        raise HTTPException(status_code=404, detail=f"Candidato non trovato: {rnacentral_id}")

    node_id = str(row_df.iloc[0].get("node_id", f"ncrna:{rnacentral_id}"))

    if "query_node_id" not in similarity.columns:
        return {"count": 0, "results": []}

    sim_df = similarity[similarity["query_node_id"].astype(str) == node_id].copy()

    if "query_rank" in sim_df.columns:
        sim_df = sim_df.sort_values("query_rank")

    sim_df = sim_df.head(limit)

    return {
        "rnacentral_id": rnacentral_id,
        "node_id": node_id,
        "count": len(sim_df),
        "results": df_records(sim_df),
        "methodological_note": "Similarità calcolata su embedding non supervisionati.",
    }


@app.post("/rag/ask")
def rag_ask(payload: RagAskRequest) -> Dict[str, Any]:
    if not payload.question.strip():
        raise HTTPException(status_code=400, detail="Domanda vuota.")

    try:
        return get_rag_answer(
            question=payload.question,
            provider=payload.provider,
            top_k=payload.top_k,
        )
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
