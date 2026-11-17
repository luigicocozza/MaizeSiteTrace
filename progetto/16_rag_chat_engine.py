#!/usr/bin/env python3
"""
16_rag_chat_engine.py

Motore Chat LLM/RAG per il progetto maize_drought_ncrna_gnn.

Supporta:
- retrieval su data/processed/rag/project_knowledge_base.jsonl
- provider Ollama locale
- provider OpenAI API se OPENAI_API_KEY è disponibile
- fallback estrattivo se nessun LLM è disponibile
- citazione dei chunk usati

Esempi:

python 16_rag_chat_engine.py --question "Quanti ncRNA sono stati analizzati?"
python 16_rag_chat_engine.py --provider extractive --question "Quali sono i limiti metodologici?"
python 16_rag_chat_engine.py --interactive
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import requests

try:
    from sklearn.feature_extraction.text import TfidfVectorizer
    from sklearn.metrics.pairwise import cosine_similarity
    SKLEARN_AVAILABLE = True
except Exception:
    TfidfVectorizer = None
    cosine_similarity = None
    SKLEARN_AVAILABLE = False


@dataclass
class RetrievedChunk:
    rank: int
    score: float
    chunk_id: str
    title: str
    source: str
    source_type: str
    text: str
    metadata: Dict[str, Any]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Motore Chat LLM/RAG per progetto ncRNA/GNN mais."
    )
    parser.add_argument(
        "--kb",
        default="data/processed/rag/project_knowledge_base.jsonl",
        help="Knowledge base JSONL.",
    )
    parser.add_argument(
        "--question",
        default="",
        help="Domanda singola.",
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Avvia modalità chat interattiva.",
    )
    parser.add_argument(
        "--provider",
        choices=["auto", "ollama", "openai", "extractive"],
        default="auto",
        help="Provider risposta.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=6,
        help="Numero chunk recuperati.",
    )
    parser.add_argument(
        "--ollama-url",
        default="http://localhost:11434/api/generate",
        help="Endpoint Ollama.",
    )
    parser.add_argument(
        "--ollama-model",
        default=os.environ.get("OLLAMA_MODEL", "llama3.1:8b"),
        help="Modello Ollama.",
    )
    parser.add_argument(
        "--openai-model",
        default=os.environ.get("OPENAI_MODEL", "gpt-4o-mini"),
        help="Modello OpenAI API.",
    )
    parser.add_argument(
        "--max-context-chars",
        type=int,
        default=9000,
        help="Massimo numero caratteri contesto passato al LLM.",
    )
    return parser.parse_args()


def load_jsonl(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"Knowledge base non trovata: {path}")

    records = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            records.append(json.loads(line))

    if not records:
        raise ValueError(f"Knowledge base vuota: {path}")

    return records


def normalize_text(text: str) -> str:
    text = text.lower()
    text = re.sub(r"[^a-zàèéìòù0-9_:\-\. ]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def keyword_score(query: str, text: str) -> float:
    q_terms = set(normalize_text(query).split())
    t_terms = normalize_text(text).split()

    if not q_terms or not t_terms:
        return 0.0

    t_counter = {}
    for term in t_terms:
        t_counter[term] = t_counter.get(term, 0) + 1

    score = 0.0
    for term in q_terms:
        score += t_counter.get(term, 0)

    return score / max(1, len(t_terms))


class RagRetriever:
    def __init__(self, records: List[Dict[str, Any]]):
        self.records = records
        self.texts = [
            f"{r.get('title', '')}\n{r.get('text', '')}\n{json.dumps(r.get('metadata', {}), ensure_ascii=False)}"
            for r in records
        ]

        self.vectorizer = None
        self.matrix = None

        if SKLEARN_AVAILABLE:
            self.vectorizer = TfidfVectorizer(
                lowercase=True,
                strip_accents="unicode",
                ngram_range=(1, 2),
                max_features=50000,
            )
            self.matrix = self.vectorizer.fit_transform(self.texts)

    def search(self, query: str, top_k: int = 6) -> List[RetrievedChunk]:
        if self.vectorizer is not None and self.matrix is not None and cosine_similarity is not None:
            q = self.vectorizer.transform([query])
            sims = cosine_similarity(q, self.matrix).ravel()
            order = sims.argsort()[::-1][:top_k]

            chunks = []
            for rank, idx in enumerate(order, start=1):
                rec = self.records[int(idx)]
                chunks.append(
                    RetrievedChunk(
                        rank=rank,
                        score=float(sims[int(idx)]),
                        chunk_id=str(rec.get("chunk_id", "")),
                        title=str(rec.get("title", "")),
                        source=str(rec.get("source", "")),
                        source_type=str(rec.get("source_type", "")),
                        text=str(rec.get("text", "")),
                        metadata=rec.get("metadata", {}) if isinstance(rec.get("metadata", {}), dict) else {},
                    )
                )
            return chunks

        scored = []
        for idx, text in enumerate(self.texts):
            scored.append((idx, keyword_score(query, text)))

        scored.sort(key=lambda x: x[1], reverse=True)
        chunks = []

        for rank, (idx, score) in enumerate(scored[:top_k], start=1):
            rec = self.records[idx]
            chunks.append(
                RetrievedChunk(
                    rank=rank,
                    score=float(score),
                    chunk_id=str(rec.get("chunk_id", "")),
                    title=str(rec.get("title", "")),
                    source=str(rec.get("source", "")),
                    source_type=str(rec.get("source_type", "")),
                    text=str(rec.get("text", "")),
                    metadata=rec.get("metadata", {}) if isinstance(rec.get("metadata", {}), dict) else {},
                )
            )

        return chunks


def build_context(chunks: List[RetrievedChunk], max_chars: int) -> str:
    parts = []

    for c in chunks:
        block = f"""
[CHUNK {c.rank}]
chunk_id: {c.chunk_id}
title: {c.title}
source: {c.source}
source_type: {c.source_type}
score: {c.score:.6f}

{c.text}
""".strip()

        parts.append(block)

    context = "\n\n---\n\n".join(parts)

    if len(context) > max_chars:
        context = context[:max_chars] + "\n\n[CONTESTO TRONCATO]"

    return context


def build_prompt(question: str, context: str) -> str:
    return f"""
Sei l'assistente RAG del progetto maize_drought_ncrna_gnn.

Regole obbligatorie:
1. Rispondi in italiano.
2. Usa solo il contesto fornito.
3. Non dichiarare che un ncRNA è realmente drought-responsive se il contesto non lo dimostra.
4. Distingui sempre tra:
   - risultato computazionale/esplorativo;
   - evidenza biologica validata.
5. Quando parli del ranking, specifica che è esplorativo e non supervisionato.
6. Se l'informazione non è nel contesto, dillo chiaramente.
7. Alla fine indica i chunk usati.

DOMANDA:
{question}

CONTESTO:
{context}

RISPOSTA:
""".strip()


def ollama_available(url: str) -> bool:
    base = url.replace("/api/generate", "")
    try:
        r = requests.get(base + "/api/tags", timeout=2)
        return r.status_code == 200
    except Exception:
        return False


def answer_with_ollama(prompt: str, url: str, model: str) -> str:
    payload = {
        "model": model,
        "prompt": prompt,
        "stream": False,
        "options": {
            "temperature": 0.2,
            "num_ctx": 8192,
        },
    }

    r = requests.post(url, json=payload, timeout=180)
    r.raise_for_status()
    obj = r.json()
    return str(obj.get("response", "")).strip()


def openai_available() -> bool:
    return bool(os.environ.get("OPENAI_API_KEY"))


def answer_with_openai(prompt: str, model: str) -> str:
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError("OPENAI_API_KEY non impostata.")

    url = "https://api.openai.com/v1/chat/completions"

    payload = {
        "model": model,
        "messages": [
            {
                "role": "system",
                "content": (
                    "Sei un assistente RAG scientifico prudente. "
                    "Rispondi solo con il contesto fornito e non fare overclaiming."
                ),
            },
            {
                "role": "user",
                "content": prompt,
            },
        ],
        "temperature": 0.2,
    }

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
    }

    r = requests.post(url, headers=headers, json=payload, timeout=180)
    r.raise_for_status()
    obj = r.json()
    return obj["choices"][0]["message"]["content"].strip()


def extractive_answer(question: str, chunks: List[RetrievedChunk]) -> str:
    lines = []
    lines.append("Risposta estrattiva basata sui chunk recuperati.")
    lines.append("")

    if not chunks:
        lines.append("Non ho recuperato chunk pertinenti dalla knowledge base.")
        return "\n".join(lines)

    q_norm = normalize_text(question)
    q_terms = [t for t in q_norm.split() if len(t) >= 3]

    selected_sentences = []

    for c in chunks:
        text = re.sub(r"\s+", " ", c.text)
        sentences = re.split(r"(?<=[\.\?\!])\s+", text)

        local_hits = []
        for s in sentences:
            s_norm = normalize_text(s)
            hit_count = sum(1 for term in q_terms if term in s_norm)
            if hit_count > 0:
                local_hits.append((hit_count, s))

        local_hits.sort(key=lambda x: x[0], reverse=True)

        for _, sentence in local_hits[:2]:
            selected_sentences.append((c, sentence))

    if not selected_sentences:
        for c in chunks[:3]:
            selected_sentences.append((c, c.text[:500]))

    lines.append("Sintesi:")
    for c, sentence in selected_sentences[:8]:
        lines.append(f"- {sentence.strip()} [chunk: {c.chunk_id}]")

    lines.append("")
    lines.append(
        "Nota metodologica: questa risposta è estrattiva. "
        "Il ranking e la similarità restano risultati esplorativi e non provano risposta alla siccità."
    )

    return "\n".join(lines)


def choose_provider(args: argparse.Namespace) -> str:
    if args.provider != "auto":
        return args.provider

    if ollama_available(args.ollama_url):
        return "ollama"

    if openai_available():
        return "openai"

    return "extractive"




def read_csv_rows_for_structured_answers(path: str) -> List[Dict[str, str]]:
    p = Path(path)
    if not p.exists():
        return []
    try:
        with p.open("r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))
    except Exception:
        return []


def read_json_summary_for_structured_answers(path: str) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        return {}
    try:
        with p.open("r", encoding="utf-8") as f:
            obj = json.load(f)
        if not isinstance(obj, dict):
            return {}
        summary = obj.get("summary", obj)
        return summary if isinstance(summary, dict) else {}
    except Exception:
        return {}


def structured_project_answer(question: str) -> str | None:
    """
    Risposte deterministiche per domande ricorrenti del progetto.
    Servono a migliorare il fallback senza LLM e la futura GUI.
    """
    q = normalize_text(question)

    seq = read_json_summary_for_structured_answers("outputs/maize_sequences_report.json")
    expr = read_json_summary_for_structured_answers("outputs/maize_expression_summary_report.json")
    graph = read_json_summary_for_structured_answers("outputs/graph_tables_report.json")
    pyg = read_json_summary_for_structured_answers("outputs/pyg_dataset_report.json")
    emb = read_json_summary_for_structured_answers("outputs/embeddings/unsupervised_embeddings_report.json")
    sim = read_json_summary_for_structured_answers("outputs/similarity/ncrna_similarity_report.json")
    ranking = read_json_summary_for_structured_answers("outputs/ranking/ncrna_exploratory_ranking_report.json")
    validation = read_json_summary_for_structured_answers("outputs/validation_report.json")
    methodological = read_json_summary_for_structured_answers("outputs/reports/methodological_report.json")

    top_rows = read_csv_rows_for_structured_answers("data/processed/ranking/ncrna_exploratory_top50.csv")

    wants_count_ncrna = (
        ("quanti" in q or "numero" in q or "analizzati" in q)
        and ("ncrna" in q or "rna" in q)
    )

    wants_graph_counts = (
        ("nodi" in q or "archi" in q or "grafo" in q)
        and ("quanti" in q or "numero" in q or "ha" in q)
    )

    wants_top_candidates = (
        ("primi" in q or "top" in q or "candidati" in q or "ranking" in q)
        and ("candidati" in q or "ranking" in q)
    )

    asks_drought_proof = (
        ("drought" in q or "siccita" in q or "siccità" in q)
        and ("dimostra" in q or "prova" in q or "considerarli" in q or "responsive" in q)
    )

    wants_limits = "limiti" in q or "limitazioni" in q or "metodologici" in q
    wants_status = "stato" in q and "progetto" in q

    if asks_drought_proof:
        return (
            "No. Il ranking non dimostra che i candidati siano realmente drought-responsive.\n\n"
            "Nel progetto attuale il ranking è esplorativo, non supervisionato e basato su embedding, "
            "similarità, feature di sequenza e contesto di grafo. Le label sono mantenute come unknown "
            "(-1) e non sono disponibili valori quantitativi RNA-seq ncRNA-specifici, come log2 fold-change "
            "o p-value. Quindi i candidati vanno considerati prioritari per revisione biologica, "
            "integrazione con dati RNA-seq e validazione sperimentale, non come scoperte già validate."
        )

    if wants_count_ncrna:
        accepted = seq.get("accepted_records", seq.get("n_rows", "NA"))
        ranked = ranking.get("n_ranked_ncrna", "NA")
        sim_n = sim.get("n_ncrna", "NA")
        rna_dist = seq.get("rna_type_distribution", {})

        return (
            f"Sono stati analizzati {accepted} ncRNA nella fase di costruzione delle sequenze.\n\n"
            f"Nel ranking esplorativo risultano {ranked} ncRNA classificati, e nella similarità risultano "
            f"{sim_n} ncRNA confrontati.\n\n"
            f"Distribuzione dei tipi RNA: {rna_dist}.\n\n"
            "Limite metodologico principale: non sono presenti label biologiche validate né dati quantitativi "
            "RNA-seq ncRNA-specifici; quindi il sistema non può ancora affermare che un ncRNA sia realmente "
            "coinvolto nella risposta alla siccità."
        )

    if wants_graph_counts:
        n_nodes = graph.get("n_nodes", pyg.get("n_nodes", "NA"))
        n_edges = graph.get("n_edges", pyg.get("n_edges", "NA"))
        node_dist = graph.get("node_type_distribution", {})
        edge_dist = graph.get("edge_type_distribution", {})

        return (
            f"Il grafo contiene {n_nodes} nodi e {n_edges} archi.\n\n"
            f"Distribuzione dei nodi: {node_dist}.\n\n"
            f"Distribuzione degli archi: {edge_dist}.\n\n"
            "Gli archi rappresentano relazioni di provenienza, annotazione e contesto; non dimostrano "
            "automaticamente funzione biologica nella siccità."
        )

    if wants_top_candidates and top_rows:
        n = 5
        if "10" in q:
            n = 10

        lines = []
        lines.append(f"Primi {n} candidati del ranking esplorativo:")
        lines.append("")

        for row in top_rows[:n]:
            lines.append(
                f"{row.get('exploratory_rank')}. {row.get('rnacentral_id')} | "
                f"{row.get('rna_type')} | score={row.get('exploratory_rank_score')} | "
                f"top3={row.get('top3_similar_ids')}"
            )

        lines.append("")
        lines.append(
            "Nota metodologica: questi candidati sono prioritari in senso computazionale/esplorativo. "
            "Non devono essere considerati automaticamente drought-responsive senza evidenze biologiche, "
            "RNA-seq o validazione sperimentale."
        )

        return "\n".join(lines)

    if wants_limits:
        return (
            "Limiti metodologici principali:\n\n"
            "1. Le label sono unknown (-1), quindi non è ancora possibile addestrare o valutare "
            "un classificatore supervisionato di risposta alla siccità.\n"
            "2. I valori di espressione non sono quantitativi: mancano log2 fold-change, p-value "
            "e matrici RNA-seq ncRNA-specifiche.\n"
            "3. Il grafo contiene soprattutto relazioni di provenienza, annotazione, letteratura e contesto; "
            "non contiene ancora interazioni ncRNA-target validate.\n"
            "4. Similarità e ranking sono basati su embedding non supervisionati e non dimostrano funzione biologica.\n"
            "5. I candidati prioritari devono essere verificati con letteratura, dataset RNA-seq e validazione biologica."
        )

    if wants_status:
        validation_status = validation.get("validation_status", methodological.get("validation_status", "NA"))
        n_nodes = pyg.get("n_nodes", graph.get("n_nodes", "NA"))
        n_edges = pyg.get("n_edges", graph.get("n_edges", "NA"))
        emb_dim = emb.get("embedding_dim_created", "NA")
        ranked = ranking.get("n_ranked_ncrna", "NA")
        inventory_existing = methodological.get("inventory_existing", "NA")
        inventory_total = methodological.get("inventory_total", "NA")

        return (
            "Stato attuale del progetto:\n\n"
            f"- Validazione: {validation_status}\n"
            f"- Inventario output: {inventory_existing}/{inventory_total} file presenti\n"
            f"- Grafo PyG: {n_nodes} nodi e {n_edges} archi\n"
            f"- Embedding: {emb_dim} dimensioni\n"
            f"- ncRNA nel ranking: {ranked}\n"
            "- Knowledge base RAG: costruita\n"
            "- Motore chat RAG: funzionante con fallback estrattivo\n\n"
            "Il progetto è nella fase software/RAG/GUI. L'articolo scientifico va preparato solo dopo "
            "completamento e validazione end-to-end della GUI, della chat e dei moduli di esportazione."
        )

    return None



def answer_question(
    question: str,
    retriever: RagRetriever,
    args: argparse.Namespace,
) -> Tuple[str, List[RetrievedChunk], str]:
    chunks = retriever.search(question, top_k=args.top_k)

    structured = structured_project_answer(question)
    if structured is not None:
        return structured, chunks, "structured_rag"

    context = build_context(chunks, max_chars=args.max_context_chars)
    prompt = build_prompt(question, context)
    provider = choose_provider(args)

    try:
        if provider == "ollama":
            answer = answer_with_ollama(prompt, url=args.ollama_url, model=args.ollama_model)
        elif provider == "openai":
            answer = answer_with_openai(prompt, model=args.openai_model)
        elif provider == "extractive":
            answer = extractive_answer(question, chunks)
        else:
            answer = extractive_answer(question, chunks)
            provider = "extractive"
    except Exception as exc:
        provider = "extractive_fallback"
        answer = (
            f"Provider LLM non disponibile o errore: {exc}\n\n"
            + extractive_answer(question, chunks)
        )

    return answer, chunks, provider


def format_sources(chunks: List[RetrievedChunk]) -> str:
    lines = []
    lines.append("")
    lines.append("Fonti/chunk usati:")
    for c in chunks:
        lines.append(
            f"- [{c.rank}] {c.chunk_id} | {c.title} | "
            f"source_type={c.source_type} | score={c.score:.4f}"
        )
    return "\n".join(lines)


def print_answer(question: str, answer: str, chunks: List[RetrievedChunk], provider: str) -> None:
    print("")
    print("=" * 80)
    print("DOMANDA")
    print("=" * 80)
    print(question)
    print("")
    print("=" * 80)
    print(f"RISPOSTA | provider={provider}")
    print("=" * 80)
    print(answer.strip())
    print(format_sources(chunks))
    print("")


def interactive_loop(retriever: RagRetriever, args: argparse.Namespace) -> None:
    print("")
    print("Chat RAG progetto maize_drought_ncrna_gnn")
    print("Scrivi 'exit' o 'quit' per uscire.")
    print("")

    while True:
        try:
            question = input("Domanda > ").strip()
        except EOFError:
            break

        if question.lower() in {"exit", "quit", "q"}:
            break

        if not question:
            continue

        answer, chunks, provider = answer_question(question, retriever, args)
        print_answer(question, answer, chunks, provider)


def main() -> int:
    args = parse_args()

    kb_path = Path(args.kb)
    records = load_jsonl(kb_path)
    retriever = RagRetriever(records)

    if args.interactive:
        interactive_loop(retriever, args)
        return 0

    question = args.question.strip()

    if not question:
        print("Errore: usa --question oppure --interactive")
        return 1

    answer, chunks, provider = answer_question(question, retriever, args)
    print_answer(question, answer, chunks, provider)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
