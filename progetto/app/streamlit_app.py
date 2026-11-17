#!/usr/bin/env python3
"""
GUI Streamlit completa per maize_drought_ncrna_gnn.

Avvio:
./run_gui.sh

Sezioni:
- Dashboard
- Ranking candidati
- Dettaglio candidato
- Similarità
- Embedding
- Chat RAG
- Report metodologico
- Download
"""

from __future__ import annotations

import importlib.util
import io
import json
import zipfile
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Tuple

import pandas as pd
import streamlit as st


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
    "method_report_md": BASE_DIR / "outputs/reports/methodological_report.md",
    "method_summary_txt": BASE_DIR / "outputs/reports/methodological_report_summary.txt",
    "embedding_plot": BASE_DIR / "outputs/embeddings/plots/ncrna_embedding_2d.png",
    "ranking_plot": BASE_DIR / "outputs/ranking/plots/top20_exploratory_ranking.png",
    "similarity_plot": BASE_DIR / "outputs/similarity/plots/similarity_score_histogram.png",
    "rag_engine": BASE_DIR / "16_rag_chat_engine.py",
}


st.set_page_config(
    page_title="Maize ncRNA GNN",
    page_icon="🌽",
    layout="wide",
)


@st.cache_data
def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


@st.cache_data
def load_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


@st.cache_resource
def load_rag_module():
    engine_path = PATHS["rag_engine"]

    if not engine_path.exists():
        return None

    spec = importlib.util.spec_from_file_location("rag_chat_engine", engine_path)
    if spec is None or spec.loader is None:
        return None

    module = importlib.util.module_from_spec(spec)

    # Necessario con Python 3.14/dataclasses:
    # il modulo deve essere registrato in sys.modules prima di exec_module.
    sys.modules[spec.name] = module

    spec.loader.exec_module(module)
    return module


@st.cache_resource
def load_rag_retriever():
    module = load_rag_module()

    if module is None:
        return None, "Motore RAG non trovato."

    kb_path = PATHS["rag_kb"]
    if not kb_path.exists():
        return None, "Knowledge base RAG non trovata."

    try:
        records = module.load_jsonl(kb_path)
        retriever = module.RagRetriever(records)
        return retriever, ""
    except Exception as exc:
        return None, f"Errore caricamento RAG: {exc}"


def build_rag_args(provider: str, top_k: int):
    return SimpleNamespace(
        kb=str(PATHS["rag_kb"]),
        question="",
        interactive=False,
        provider=provider,
        top_k=top_k,
        ollama_url="http://localhost:11434/api/generate",
        ollama_model="llama3.1:8b",
        openai_model="gpt-4o-mini",
        max_context_chars=9000,
    )


def ask_rag(question: str, provider: str = "auto", top_k: int = 6) -> Tuple[str, str, str]:
    module = load_rag_module()
    retriever, error = load_rag_retriever()

    if module is None:
        return "Motore RAG non disponibile: file 16_rag_chat_engine.py non trovato.", "", "error"

    if retriever is None:
        return f"Retriever RAG non disponibile: {error}", "", "error"

    args = build_rag_args(provider=provider, top_k=top_k)

    try:
        answer, chunks, used_provider = module.answer_question(question, retriever, args)
        sources = module.format_sources(chunks)
        return answer, sources, used_provider
    except Exception as exc:
        return f"Errore durante la risposta RAG: {exc}", "", "error"


def methodological_warning() -> None:
    st.warning(
        "Nota metodologica: ranking, similarità ed embedding sono esplorativi e non supervisionati. "
        "Non dimostrano da soli che un ncRNA sia realmente coinvolto nella risposta alla siccità."
    )


def make_arrow_safe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Evita errori/warning PyArrow di Streamlit quando una colonna contiene
    tipi misti, per esempio interi, float e stringhe nella stessa colonna.
    La conversione è solo per la visualizzazione, non modifica i CSV originali.
    """
    if df is None or df.empty:
        return df
    return df.copy().astype(str)


def safe_dataframe(df: pd.DataFrame, **kwargs) -> None:
    st.dataframe(make_arrow_safe(df), **kwargs)


def download_button_for_df(df: pd.DataFrame, label: str, filename: str) -> None:
    if df.empty:
        return

    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=label,
        data=csv_bytes,
        file_name=filename,
        mime="text/csv",
    )


def file_status_table() -> pd.DataFrame:
    rows = []

    for key, path in PATHS.items():
        rows.append(
            {
                "artifact": key,
                "path": str(path.relative_to(BASE_DIR)),
                "exists": path.exists(),
                "size_kb": round(path.stat().st_size / 1024, 2) if path.exists() else 0.0,
            }
        )

    return pd.DataFrame(rows)


def page_dashboard(
    ranking_df: pd.DataFrame,
    similarity_df: pd.DataFrame,
    features_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
) -> None:
    st.title("🌽 Maize ncRNA GNN Dashboard")
    st.subheader("Framework AI/GNN per ncRNA del mais")

    methodological_warning()

    c1, c2, c3, c4 = st.columns(4)

    with c1:
        st.metric("ncRNA nel ranking", len(ranking_df))
    with c2:
        st.metric("Feature ncRNA", len(features_df))
    with c3:
        st.metric("Nodi grafo", len(nodes_df))
    with c4:
        st.metric("Archi grafo", len(edges_df))

    c5, c6, c7, c8 = st.columns(4)

    with c5:
        st.metric("Similarità top-k", len(similarity_df))
    with c6:
        if "rna_type" in ranking_df.columns:
            st.metric("Tipi RNA", ranking_df["rna_type"].nunique())
        else:
            st.metric("Tipi RNA", "NA")
    with c7:
        if "exploratory_rank_score" in ranking_df.columns:
            st.metric("Score massimo", round(float(ranking_df["exploratory_rank_score"].max()), 6))
        else:
            st.metric("Score massimo", "NA")
    with c8:
        if "label" in ranking_df.columns:
            st.metric("Label uniche", ranking_df["label"].nunique())
        else:
            st.metric("Label uniche", "NA")

    st.divider()

    col_a, col_b = st.columns(2)

    with col_a:
        st.subheader("Distribuzione RNA type")
        if not ranking_df.empty and "rna_type" in ranking_df.columns:
            st.bar_chart(ranking_df["rna_type"].value_counts())

    with col_b:
        st.subheader("Distribuzione node type")
        if not nodes_df.empty and "node_type" in nodes_df.columns:
            st.bar_chart(nodes_df["node_type"].value_counts())

    st.divider()

    st.subheader("Grafici prodotti dalla pipeline")

    g1, g2, g3 = st.columns(3)

    with g1:
        if PATHS["embedding_plot"].exists():
            st.image(str(PATHS["embedding_plot"]), caption="Embedding ncRNA 2D")
        else:
            st.info("Plot embedding non trovato.")

    with g2:
        if PATHS["ranking_plot"].exists():
            st.image(str(PATHS["ranking_plot"]), caption="Top 20 ranking")
        else:
            st.info("Plot ranking non trovato.")

    with g3:
        if PATHS["similarity_plot"].exists():
            st.image(str(PATHS["similarity_plot"]), caption="Distribuzione similarità")
        else:
            st.info("Plot similarità non trovato.")

    st.divider()

    st.subheader("Stato file principali")
    safe_dataframe(file_status_table(), width="stretch")


def page_ranking(ranking_df: pd.DataFrame) -> None:
    st.title("🏆 Ranking esplorativo ncRNA")
    methodological_warning()

    if ranking_df.empty:
        st.error("Ranking non disponibile.")
        return

    st.sidebar.header("Filtri ranking")

    rna_types = sorted(ranking_df["rna_type"].dropna().unique().tolist()) if "rna_type" in ranking_df.columns else []
    selected_types = st.sidebar.multiselect("Tipo RNA", rna_types, default=rna_types)

    if "exploratory_rank_score" in ranking_df.columns:
        score_min = float(ranking_df["exploratory_rank_score"].min())
        score_max = float(ranking_df["exploratory_rank_score"].max())
    else:
        score_min, score_max = 0.0, 1.0

    selected_score = st.sidebar.slider(
        "Score minimo",
        min_value=score_min,
        max_value=score_max,
        value=score_min,
        step=0.001,
    )

    top_n = st.sidebar.slider(
        "Righe da mostrare",
        min_value=10,
        max_value=max(10, min(250, len(ranking_df))),
        value=min(50, len(ranking_df)),
        step=10,
    )

    filtered = ranking_df.copy()

    if selected_types and "rna_type" in filtered.columns:
        filtered = filtered[filtered["rna_type"].isin(selected_types)]

    if "exploratory_rank_score" in filtered.columns:
        filtered = filtered[filtered["exploratory_rank_score"] >= selected_score]

    if "exploratory_rank" in filtered.columns:
        filtered = filtered.sort_values("exploratory_rank")

    filtered = filtered.head(top_n)

    safe_dataframe(filtered, width="stretch", height=520)
    download_button_for_df(filtered, "Scarica ranking filtrato", "ranking_filtrato.csv")



def get_candidate_related_tables(
    selected: str,
    row: Dict[str, Any],
    similarity_df: pd.DataFrame,
    features_df: pd.DataFrame,
    embeddings_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    node_id = row.get("node_id", f"ncrna:{selected}")

    if not similarity_df.empty and "query_node_id" in similarity_df.columns:
        sim = similarity_df[similarity_df["query_node_id"].astype(str) == str(node_id)].copy()
    else:
        sim = pd.DataFrame()

    if not features_df.empty and "node_id" in features_df.columns:
        feat = features_df[features_df["node_id"].astype(str) == str(node_id)].copy()
    else:
        feat = pd.DataFrame()

    if not embeddings_df.empty and "node_id" in embeddings_df.columns:
        emb = embeddings_df[embeddings_df["node_id"].astype(str) == str(node_id)].copy()
    else:
        emb = pd.DataFrame()

    return sim, feat, emb


def as_text_value(value: Any, default: str = "NA") -> str:
    if value is None:
        return default
    try:
        if pd.isna(value):
            return default
    except Exception:
        pass
    return str(value)


def build_candidate_explanation(
    row: Dict[str, Any],
    sim: pd.DataFrame,
    feat: pd.DataFrame,
) -> str:
    rid = as_text_value(row.get("rnacentral_id"))
    rna_type = as_text_value(row.get("rna_type"))
    rank = as_text_value(row.get("exploratory_rank"))
    score = as_text_value(row.get("exploratory_rank_score"))
    length = as_text_value(row.get("sequence_length"))
    gc = as_text_value(row.get("gc_frac"))
    label = as_text_value(row.get("label"))
    mean_sim = as_text_value(row.get("mean_topk_similarity"))
    max_sim = as_text_value(row.get("max_topk_similarity"))
    degree = as_text_value(row.get("degree_total"))
    top3 = as_text_value(row.get("top3_similar_ids"))

    has_quant = as_text_value(row.get("has_quantitative_expression", "0"))
    has_drought = as_text_value(row.get("has_drought_expression_evidence", "0"))

    n_sim = len(sim)
    first_similar = []

    if not sim.empty:
        for _, srow in sim.head(5).iterrows():
            first_similar.append(
                f"{srow.get('similar_rnacentral_id', 'NA')} "
                f"({srow.get('similar_rna_type', 'NA')}, "
                f"cos={srow.get('cosine_similarity', 'NA')})"
            )

    similar_text = "; ".join(first_similar) if first_similar else "non disponibile"

    explanation = f"""
### Interpretazione del candidato {rid}

**{rid}** è un candidato ncRNA di tipo **{rna_type}** con rank esplorativo **{rank}** e score **{score}**.

Dal punto di vista computazionale, il candidato presenta:

- lunghezza sequenza: **{length}**
- GC fraction: **{gc}**
- grado nel grafo: **{degree}**
- similarità media top-k: **{mean_sim}**
- similarità massima top-k: **{max_sim}**
- top 3 simili: **{top3}**
- numero di vicini simili disponibili: **{n_sim}**

I primi candidati simili recuperati sono:

**{similar_text}**

### Significato del ranking

Il punteggio deriva da una combinazione di similarità negli embedding, centralità nello spazio embedding, tipo RNA, tipicità GC e qualità della sequenza. Questo rende il candidato interessante per una revisione successiva, ma **non rappresenta una prova biologica diretta**.

### Stato delle evidenze

- Label attuale: **{label}**
- Evidenza quantitativa di espressione disponibile: **{has_quant}**
- Evidenza diretta drought-expression disponibile: **{has_drought}**

### Nota metodologica

Il candidato deve essere considerato **prioritario in senso computazionale/esplorativo**, non automaticamente drought-responsive. Per trasformarlo in un candidato biologicamente supportato servono almeno una o più delle seguenti evidenze:

1. dati RNA-seq quantitativi ncRNA-specifici;
2. associazione con condizioni di stress idrico;
3. letteratura specifica;
4. collegamento con target genici o pathway;
5. validazione sperimentale o dataset indipendente.
"""
    return explanation.strip()


def build_candidate_markdown_report(
    row: Dict[str, Any],
    sim: pd.DataFrame,
    feat: pd.DataFrame,
    emb: pd.DataFrame,
) -> str:
    rid = as_text_value(row.get("rnacentral_id"))
    lines = []

    lines.append(f"# Report candidato ncRNA: {rid}")
    lines.append("")
    lines.append("## 1. Sintesi")
    lines.append("")
    lines.append(build_candidate_explanation(row, sim, feat))
    lines.append("")

    lines.append("## 2. Dati ranking")
    lines.append("")
    for key, value in row.items():
        lines.append(f"- **{key}**: {value}")
    lines.append("")

    lines.append("## 3. Top-k simili")
    lines.append("")
    if sim.empty:
        lines.append("Nessuna tabella di similarità disponibile.")
    else:
        for _, srow in sim.head(20).iterrows():
            lines.append(
                f"- rank {srow.get('query_rank', 'NA')}: "
                f"{srow.get('similar_rnacentral_id', 'NA')} | "
                f"type={srow.get('similar_rna_type', 'NA')} | "
                f"cosine={srow.get('cosine_similarity', 'NA')} | "
                f"same_type={srow.get('same_rna_type', 'NA')}"
            )
    lines.append("")

    lines.append("## 4. Feature")
    lines.append("")
    if feat.empty:
        lines.append("Feature non disponibili.")
    else:
        for key, value in feat.iloc[0].to_dict().items():
            lines.append(f"- **{key}**: {value}")
    lines.append("")

    lines.append("## 5. Embedding")
    lines.append("")
    if emb.empty:
        lines.append("Embedding non disponibile.")
    else:
        emb_row = emb.iloc[0].to_dict()
        emb_keys = [k for k in emb_row.keys() if str(k).startswith("emb_")]
        for key in emb_keys:
            lines.append(f"- **{key}**: {emb_row.get(key)}")
    lines.append("")

    lines.append("## 6. Limite interpretativo")
    lines.append("")
    lines.append(
        "Il report è generato automaticamente dalla GUI. Il ranking è esplorativo e non supervisionato; "
        "non dimostra che il candidato sia realmente coinvolto nella risposta alla siccità."
    )

    return "\n".join(lines) + "\n"


def build_candidate_json_payload(
    row: Dict[str, Any],
    sim: pd.DataFrame,
    feat: pd.DataFrame,
    emb: pd.DataFrame,
) -> Dict[str, Any]:
    return {
        "candidate": {k: as_text_value(v) for k, v in row.items()},
        "topk_similarity": sim.astype(str).to_dict(orient="records") if not sim.empty else [],
        "features": feat.astype(str).to_dict(orient="records") if not feat.empty else [],
        "embedding": emb.astype(str).to_dict(orient="records") if not emb.empty else [],
        "methodological_note": (
            "Exploratory unsupervised candidate. Not validated evidence of drought responsiveness."
        ),
    }


def page_candidate_detail(
    ranking_df: pd.DataFrame,
    similarity_df: pd.DataFrame,
    features_df: pd.DataFrame,
    embeddings_df: pd.DataFrame,
) -> None:
    st.title("🔎 Dettaglio candidato ncRNA")
    methodological_warning()

    if ranking_df.empty:
        st.error("Ranking non disponibile.")
        return

    candidate_ids = ranking_df["rnacentral_id"].astype(str).tolist()

    selected = st.selectbox(
        "Seleziona RNAcentral ID",
        candidate_ids,
        index=0,
    )

    candidate = ranking_df[ranking_df["rnacentral_id"].astype(str) == selected]

    if candidate.empty:
        st.error("Candidato non trovato.")
        return

    row = candidate.iloc[0].to_dict()
    node_id = row.get("node_id", f"ncrna:{selected}")

    sim, feat, emb = get_candidate_related_tables(
        selected=selected,
        row=row,
        similarity_df=similarity_df,
        features_df=features_df,
        embeddings_df=embeddings_df,
    )

    c1, c2, c3, c4 = st.columns(4)

    with c1:
        st.metric("Rank", row.get("exploratory_rank", "NA"))
    with c2:
        st.metric("Tipo RNA", row.get("rna_type", "NA"))
    with c3:
        st.metric("Score", row.get("exploratory_rank_score", "NA"))
    with c4:
        st.metric("Label", row.get("label", "NA"))

    c5, c6, c7, c8 = st.columns(4)

    with c5:
        st.metric("Lunghezza", row.get("sequence_length", "NA"))
    with c6:
        st.metric("GC frac", row.get("gc_frac", "NA"))
    with c7:
        st.metric("Mean top-k sim", row.get("mean_topk_similarity", "NA"))
    with c8:
        st.metric("Degree", row.get("degree_total", "NA"))

    explanation = build_candidate_explanation(row, sim, feat)
    markdown_report = build_candidate_markdown_report(row, sim, feat, emb)
    json_payload = build_candidate_json_payload(row, sim, feat, emb)

    tab0, tab1, tab2, tab3, tab4, tab5 = st.tabs(
        [
            "Sintesi avanzata",
            "Ranking",
            "Top-k simili",
            "Feature",
            "Embedding",
            "Export",
        ]
    )

    with tab0:
        st.markdown(explanation)

        st.divider()
        st.subheader("Spiegazione RAG personalizzata")

        default_question = (
            f"Spiegami il candidato {selected}, i suoi simili e dimmi "
            f"se posso considerarlo drought-responsive."
        )

        question = st.text_input(
            "Domanda sul candidato",
            value=default_question,
            key=f"candidate_question_{selected}",
        )

        if st.button("Interroga RAG sul candidato", key=f"rag_candidate_{selected}"):
            answer, sources, provider = ask_rag(question, provider="auto", top_k=8)
            st.caption(f"Provider usato: {provider}")
            st.markdown(answer)
            with st.expander("Chunk/fonti usati"):
                st.text(sources)

    with tab1:
        st.subheader("Riga ranking")
        safe_dataframe(candidate.T.rename(columns={candidate.index[0]: "value"}), width="stretch")

    with tab2:
        st.subheader("Top-k simili")
        if sim.empty:
            st.info("Similarità non disponibile per questo candidato.")
        else:
            safe_dataframe(sim, width="stretch", height=420)
            download_button_for_df(sim, "Scarica top-k simili CSV", f"{selected}_topk_simili.csv")

    with tab3:
        st.subheader("Feature del candidato")
        if feat.empty:
            st.info("Feature non trovate per questo candidato.")
        else:
            safe_dataframe(feat.T.rename(columns={feat.index[0]: "value"}), width="stretch")

    with tab4:
        st.subheader("Embedding del candidato")
        if emb.empty:
            st.info("Embedding non trovato per questo candidato.")
        else:
            safe_dataframe(emb.T.rename(columns={emb.index[0]: "value"}), width="stretch")

    with tab5:
        st.subheader("Esportazione candidato")

        candidate_csv = candidate.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Scarica riga candidato CSV",
            data=candidate_csv,
            file_name=f"{selected}_candidate_ranking.csv",
            mime="text/csv",
            key=f"download_candidate_csv_{selected}",
        )

        st.download_button(
            "Scarica report candidato Markdown",
            data=markdown_report.encode("utf-8"),
            file_name=f"{selected}_candidate_report.md",
            mime="text/markdown",
            key=f"download_candidate_md_{selected}",
        )

        st.download_button(
            "Scarica pacchetto candidato JSON",
            data=json.dumps(json_payload, indent=2, ensure_ascii=False).encode("utf-8"),
            file_name=f"{selected}_candidate_package.json",
            mime="application/json",
            key=f"download_candidate_json_{selected}",
        )

        if not sim.empty:
            download_button_for_df(
                sim,
                "Scarica tabella simili CSV",
                f"{selected}_similarity_topk.csv",
            )

        st.info(
            "Gli export sono pensati per revisione biologica/manuale. "
            "Non rappresentano validazione sperimentale."
        )





def numeric_series(df: pd.DataFrame, column: str) -> pd.Series:
    if df.empty or column not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[column], errors="coerce")


def page_advanced_search(
    ranking_df: pd.DataFrame,
    similarity_df: pd.DataFrame,
    features_df: pd.DataFrame,
) -> None:
    st.title("🔍 Ricerca avanzata ncRNA")
    methodological_warning()

    if ranking_df.empty:
        st.error("Ranking non disponibile.")
        return

    df = ranking_df.copy()

    st.subheader("Filtri principali")

    col1, col2, col3 = st.columns(3)

    with col1:
        query = st.text_input(
            "Cerca ID, node_id, tipo RNA o simili",
            value="",
            placeholder="es. URS0003C75441 oppure pre_miRNA",
        )

    with col2:
        rna_types = sorted(df["rna_type"].dropna().astype(str).unique().tolist()) if "rna_type" in df.columns else []
        selected_types = st.multiselect(
            "Tipo RNA",
            options=rna_types,
            default=rna_types,
        )

    with col3:
        top_n = st.number_input(
            "Numero massimo risultati",
            min_value=10,
            max_value=max(10, len(df)),
            value=min(100, len(df)),
            step=10,
        )

    st.subheader("Filtri numerici")

    rank_series = numeric_series(df, "exploratory_rank")
    score_series = numeric_series(df, "exploratory_rank_score")
    length_series = numeric_series(df, "sequence_length")
    gc_series = numeric_series(df, "gc_frac")

    c1, c2 = st.columns(2)

    with c1:
        if not rank_series.dropna().empty:
            rank_min = int(rank_series.min())
            rank_max = int(rank_series.max())
            rank_range = st.slider(
                "Intervallo rank",
                min_value=rank_min,
                max_value=rank_max,
                value=(rank_min, rank_max),
                step=1,
            )
        else:
            rank_range = None

        if not length_series.dropna().empty:
            length_min = int(length_series.min())
            length_max = int(length_series.max())
            length_range = st.slider(
                "Intervallo lunghezza",
                min_value=length_min,
                max_value=length_max,
                value=(length_min, length_max),
                step=1,
            )
        else:
            length_range = None

    with c2:
        if not score_series.dropna().empty:
            score_min = float(score_series.min())
            score_max = float(score_series.max())
            score_range = st.slider(
                "Intervallo score esplorativo",
                min_value=score_min,
                max_value=score_max,
                value=(score_min, score_max),
                step=0.001,
            )
        else:
            score_range = None

        if not gc_series.dropna().empty:
            gc_min = float(gc_series.min())
            gc_max = float(gc_series.max())
            gc_range = st.slider(
                "Intervallo GC fraction",
                min_value=gc_min,
                max_value=gc_max,
                value=(gc_min, gc_max),
                step=0.001,
            )
        else:
            gc_range = None

    st.subheader("Filtri evidenza")

    e1, e2, e3 = st.columns(3)

    with e1:
        only_quant = st.checkbox(
            "Solo con evidenza quantitativa",
            value=False,
        )

    with e2:
        only_drought = st.checkbox(
            "Solo con evidenza drought",
            value=False,
        )

    with e3:
        only_unknown_label = st.checkbox(
            "Solo label unknown (-1)",
            value=False,
        )

    filtered = df.copy()

    if query.strip():
        q = query.strip().lower()
        searchable_columns = [
            col for col in [
                "rnacentral_id",
                "node_id",
                "rna_type",
                "top3_similar_ids",
                "methodological_note",
            ]
            if col in filtered.columns
        ]

        mask = pd.Series(False, index=filtered.index)
        for col in searchable_columns:
            mask = mask | filtered[col].astype(str).str.lower().str.contains(q, na=False)
        filtered = filtered[mask]

    if selected_types and "rna_type" in filtered.columns:
        filtered = filtered[filtered["rna_type"].astype(str).isin(selected_types)]

    if rank_range is not None and "exploratory_rank" in filtered.columns:
        ranks = pd.to_numeric(filtered["exploratory_rank"], errors="coerce")
        filtered = filtered[(ranks >= rank_range[0]) & (ranks <= rank_range[1])]

    if score_range is not None and "exploratory_rank_score" in filtered.columns:
        scores = pd.to_numeric(filtered["exploratory_rank_score"], errors="coerce")
        filtered = filtered[(scores >= score_range[0]) & (scores <= score_range[1])]

    if length_range is not None and "sequence_length" in filtered.columns:
        lengths = pd.to_numeric(filtered["sequence_length"], errors="coerce")
        filtered = filtered[(lengths >= length_range[0]) & (lengths <= length_range[1])]

    if gc_range is not None and "gc_frac" in filtered.columns:
        gcs = pd.to_numeric(filtered["gc_frac"], errors="coerce")
        filtered = filtered[(gcs >= gc_range[0]) & (gcs <= gc_range[1])]

    if only_quant and "has_quantitative_expression" in filtered.columns:
        filtered = filtered[filtered["has_quantitative_expression"].astype(str).isin(["1", "True", "true"])]

    if only_drought and "has_drought_expression_evidence" in filtered.columns:
        filtered = filtered[filtered["has_drought_expression_evidence"].astype(str).isin(["1", "True", "true"])]

    if only_unknown_label and "label" in filtered.columns:
        filtered = filtered[filtered["label"].astype(str) == "-1"]

    if "exploratory_rank" in filtered.columns:
        filtered = filtered.sort_values("exploratory_rank")

    filtered = filtered.head(int(top_n))

    st.divider()

    st.subheader("Risultati filtrati")

    m1, m2, m3, m4 = st.columns(4)

    with m1:
        st.metric("Risultati", len(filtered))
    with m2:
        if not filtered.empty and "rna_type" in filtered.columns:
            st.metric("Tipi RNA", filtered["rna_type"].nunique())
        else:
            st.metric("Tipi RNA", "NA")
    with m3:
        if not filtered.empty and "exploratory_rank_score" in filtered.columns:
            st.metric("Score max", round(float(pd.to_numeric(filtered["exploratory_rank_score"], errors="coerce").max()), 6))
        else:
            st.metric("Score max", "NA")
    with m4:
        if not filtered.empty and "exploratory_rank" in filtered.columns:
            st.metric("Miglior rank", int(pd.to_numeric(filtered["exploratory_rank"], errors="coerce").min()))
        else:
            st.metric("Miglior rank", "NA")

    safe_dataframe(filtered, width="stretch", height=520)

    download_button_for_df(
        filtered,
        "Scarica risultati filtrati CSV",
        "ncrna_advanced_search_results.csv",
    )

    st.divider()

    left, right = st.columns(2)

    with left:
        st.subheader("Distribuzione RNA type nei risultati")
        if not filtered.empty and "rna_type" in filtered.columns:
            st.bar_chart(filtered["rna_type"].value_counts())
        else:
            st.info("Nessun risultato da visualizzare.")

    with right:
        st.subheader("Sintesi filtri")
        st.markdown(
            f"""
- Query testuale: `{query if query else 'nessuna'}`
- Tipi RNA selezionati: `{selected_types}`
- Rank range: `{rank_range}`
- Score range: `{score_range}`
- Lunghezza range: `{length_range}`
- GC range: `{gc_range}`
- Solo evidenza quantitativa: `{only_quant}`
- Solo evidenza drought: `{only_drought}`
- Solo label unknown: `{only_unknown_label}`
"""
        )

    st.divider()

    st.subheader("Ricerca nei simili dei risultati")

    if similarity_df.empty or filtered.empty:
        st.info("Similarità non disponibile o nessun risultato filtrato.")
    else:
        selected_nodes = set(filtered["node_id"].astype(str).tolist()) if "node_id" in filtered.columns else set()

        if "query_node_id" in similarity_df.columns:
            sim_filtered = similarity_df[similarity_df["query_node_id"].astype(str).isin(selected_nodes)].copy()
        else:
            sim_filtered = pd.DataFrame()

        st.write(f"Righe similarità collegate ai risultati filtrati: **{len(sim_filtered)}**")
        if not sim_filtered.empty:
            safe_dataframe(sim_filtered, width="stretch", height=360)
            download_button_for_df(
                sim_filtered,
                "Scarica similarità dei risultati filtrati",
                "ncrna_advanced_search_similarity.csv",
            )




def page_similarity(similarity_df: pd.DataFrame) -> None:
    st.title("🧬 Similarità tra ncRNA")
    methodological_warning()

    if similarity_df.empty:
        st.error("Similarità non disponibile.")
        return

    safe_dataframe(similarity_df, width="stretch", height=520)
    download_button_for_df(similarity_df, "Scarica similarità top-k", "ncrna_similarity_topk.csv")

    if "cosine_similarity" in similarity_df.columns:
        st.subheader("Similarità ordinate")
        st.line_chart(similarity_df["cosine_similarity"].sort_values(ascending=False).reset_index(drop=True))

    if PATHS["similarity_plot"].exists():
        st.image(str(PATHS["similarity_plot"]), caption="Istogramma similarità")


def page_embeddings(embeddings_df: pd.DataFrame) -> None:
    st.title("🧭 Embedding ncRNA")
    methodological_warning()

    if embeddings_df.empty:
        st.error("Embedding non disponibili.")
        return

    safe_dataframe(embeddings_df, width="stretch", height=520)
    download_button_for_df(embeddings_df, "Scarica embedding ncRNA", "ncrna_embeddings.csv")

    if {"emb_001", "emb_002"}.issubset(set(embeddings_df.columns)):
        st.subheader("Scatter embedding 2D")
        st.scatter_chart(embeddings_df, x="emb_001", y="emb_002")

    if PATHS["embedding_plot"].exists():
        st.image(str(PATHS["embedding_plot"]), caption="Plot embedding generato")


def page_rag_chat() -> None:
    st.title("💬 Chat LLM/RAG del progetto")
    methodological_warning()

    retriever, error = load_rag_retriever()

    if retriever is None:
        st.error(error)
        return

    st.success("Motore RAG caricato correttamente.")

    with st.sidebar:
        st.header("Impostazioni Chat")
        provider = st.selectbox(
            "Provider",
            ["auto", "extractive", "ollama", "openai"],
            index=0,
        )
        top_k = st.slider("Chunk recuperati", 3, 12, 6, 1)

    examples = [
        "Qual è lo stato del progetto?",
        "Quanti ncRNA sono stati analizzati e qual è il limite metodologico principale?",
        "Quali sono i primi 5 candidati del ranking?",
        "Il ranking dimostra risposta alla siccità?",
        "Quanti nodi e archi ha il grafo?",
        "Quali sono i limiti metodologici?",
    ]

    selected_example = st.selectbox("Domande rapide", [""] + examples)

    if "rag_messages" not in st.session_state:
        st.session_state.rag_messages = []

    if selected_example and st.button("Usa domanda rapida"):
        st.session_state.pending_question = selected_example

    pending = st.session_state.pop("pending_question", "")

    question = st.chat_input("Scrivi una domanda sul progetto, sui candidati, sul ranking o sui limiti...")

    if pending:
        question = pending

    for msg in st.session_state.rag_messages:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])
            if msg.get("sources"):
                with st.expander("Chunk/fonti usati"):
                    st.text(msg["sources"])

    if question:
        st.session_state.rag_messages.append(
            {
                "role": "user",
                "content": question,
                "sources": "",
            }
        )

        with st.chat_message("user"):
            st.markdown(question)

        answer, sources, used_provider = ask_rag(question, provider=provider, top_k=top_k)

        answer_with_provider = f"**Provider:** `{used_provider}`\n\n{answer}"

        st.session_state.rag_messages.append(
            {
                "role": "assistant",
                "content": answer_with_provider,
                "sources": sources,
            }
        )

        with st.chat_message("assistant"):
            st.markdown(answer_with_provider)
            with st.expander("Chunk/fonti usati"):
                st.text(sources)

    if st.button("Pulisci chat"):
        st.session_state.rag_messages = []
        st.rerun()


def page_reports() -> None:
    st.title("📄 Report metodologico")

    summary = load_text(PATHS["method_summary_txt"])
    report = load_text(PATHS["method_report_md"])

    if summary:
        st.subheader("Summary")
        st.text(summary)
    else:
        st.warning("Summary non trovato.")

    if report:
        st.subheader("Report completo")
        st.markdown(report)
    else:
        st.warning("Report metodologico non trovato.")




def df_to_markdown_no_tabulate(df: pd.DataFrame) -> str:
    """
    Converte un DataFrame in tabella Markdown senza usare pandas.to_markdown(),
    così non serve il pacchetto opzionale tabulate.
    """
    if df is None or df.empty:
        return "Nessun dato disponibile."

    safe = df.copy().astype(str)

    def esc(value: str) -> str:
        return (
            str(value)
            .replace("|", "\\|")
            .replace("\n", "<br>")
            .replace("\r", "")
        )

    headers = [esc(c) for c in safe.columns]
    lines = []
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")

    for _, row in safe.iterrows():
        values = [esc(row[col]) for col in safe.columns]
        lines.append("| " + " | ".join(values) + " |")

    return "\n".join(lines)


def build_project_export_inventory() -> pd.DataFrame:
    export_paths = [
        ("ranking_completo", BASE_DIR / "data/processed/ranking/ncrna_exploratory_ranking.csv"),
        ("ranking_top50", BASE_DIR / "data/processed/ranking/ncrna_exploratory_top50.csv"),
        ("similarity_topk", BASE_DIR / "data/processed/similarity/ncrna_similarity_topk.csv"),
        ("similarity_unique_pairs", BASE_DIR / "data/processed/similarity/ncrna_similarity_unique_pairs.csv"),
        ("features_ncrna", BASE_DIR / "data/processed/ncrna_features.csv"),
        ("embeddings_ncrna", BASE_DIR / "data/processed/embeddings/ncrna_embeddings.csv"),
        ("graph_nodes", BASE_DIR / "data/processed/graph_nodes.csv"),
        ("graph_edges", BASE_DIR / "data/processed/graph_edges.csv"),
        ("rag_knowledge_base", BASE_DIR / "data/processed/rag/project_knowledge_base.jsonl"),
        ("methodological_report_md", BASE_DIR / "outputs/reports/methodological_report.md"),
        ("methodological_report_json", BASE_DIR / "outputs/reports/methodological_report.json"),
        ("methodological_summary", BASE_DIR / "outputs/reports/methodological_report_summary.txt"),
        ("rag_report", BASE_DIR / "outputs/rag/rag_knowledge_base_report.json"),
        ("ranking_report", BASE_DIR / "outputs/ranking/ncrna_exploratory_ranking_report.json"),
        ("similarity_report", BASE_DIR / "outputs/similarity/ncrna_similarity_report.json"),
        ("embeddings_report", BASE_DIR / "outputs/embeddings/unsupervised_embeddings_report.json"),
        ("pyg_report", BASE_DIR / "outputs/pyg_dataset_report.json"),
        ("validation_report", BASE_DIR / "outputs/validation_report.json"),
    ]

    rows = []
    for name, p in export_paths:
        rows.append(
            {
                "artifact": name,
                "relative_path": str(p.relative_to(BASE_DIR)),
                "exists": p.exists(),
                "size_kb": round(p.stat().st_size / 1024, 3) if p.exists() else 0.0,
            }
        )

    return pd.DataFrame(rows)


def build_project_export_markdown(
    ranking_df: pd.DataFrame,
    similarity_df: pd.DataFrame,
    features_df: pd.DataFrame,
    embeddings_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
) -> str:
    inventory = build_project_export_inventory()

    lines = []
    lines.append("# Export report progetto maize_drought_ncrna_gnn")
    lines.append("")
    lines.append("## Sintesi computazionale")
    lines.append("")
    lines.append(f"- ncRNA nel ranking: {len(ranking_df)}")
    lines.append(f"- righe feature ncRNA: {len(features_df)}")
    lines.append(f"- righe embedding ncRNA: {len(embeddings_df)}")
    lines.append(f"- righe similarità top-k: {len(similarity_df)}")
    lines.append(f"- nodi grafo: {len(nodes_df)}")
    lines.append(f"- archi grafo: {len(edges_df)}")
    lines.append("")
    lines.append("## Nota metodologica")
    lines.append("")
    lines.append(
        "Ranking, similarità ed embedding sono risultati esplorativi e non supervisionati. "
        "Non costituiscono validazione biologica della risposta alla siccità. "
        "I candidati devono essere revisionati con letteratura, dati RNA-seq e validazione sperimentale."
    )
    lines.append("")
    lines.append("## Top 10 candidati")
    lines.append("")

    if not ranking_df.empty:
        cols = [
            c for c in [
                "exploratory_rank",
                "rnacentral_id",
                "rna_type",
                "exploratory_rank_score",
                "top3_similar_ids",
            ]
            if c in ranking_df.columns
        ]
        lines.append(df_to_markdown_no_tabulate(ranking_df[cols].head(10)))
    else:
        lines.append("Ranking non disponibile.")

    lines.append("")
    lines.append("## Inventario export")
    lines.append("")
    lines.append(df_to_markdown_no_tabulate(inventory))
    lines.append("")

    return "\n".join(lines)


def build_project_zip_bytes(
    ranking_df: pd.DataFrame,
    similarity_df: pd.DataFrame,
    features_df: pd.DataFrame,
    embeddings_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
) -> bytes:
    inventory = build_project_export_inventory()
    report_md = build_project_export_markdown(
        ranking_df=ranking_df,
        similarity_df=similarity_df,
        features_df=features_df,
        embeddings_df=embeddings_df,
        nodes_df=nodes_df,
        edges_df=edges_df,
    )

    buffer = io.BytesIO()

    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr("README_export_report.md", report_md)
        z.writestr("pipeline_export_inventory.csv", inventory.to_csv(index=False))

        for _, row in inventory.iterrows():
            if not bool(row["exists"]):
                continue

            src = BASE_DIR / str(row["relative_path"])
            if not src.exists() or src.is_dir():
                continue

            z.write(src, arcname=str(row["relative_path"]))

    buffer.seek(0)
    return buffer.getvalue()





def page_downloads(
    ranking_df: pd.DataFrame,
    similarity_df: pd.DataFrame,
    features_df: pd.DataFrame,
    embeddings_df: pd.DataFrame,
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
) -> None:
    st.title("⬇️ Download output e pacchetto progetto")

    methodological_warning()

    st.subheader("Export automatico progetto")

    inventory = build_project_export_inventory()
    report_md = build_project_export_markdown(
        ranking_df=ranking_df,
        similarity_df=similarity_df,
        features_df=features_df,
        embeddings_df=embeddings_df,
        nodes_df=nodes_df,
        edges_df=edges_df,
    )

    c1, c2, c3 = st.columns(3)

    with c1:
        st.metric("File export presenti", int(inventory["exists"].sum()))

    with c2:
        st.metric("File export tracciati", len(inventory))

    with c3:
        st.metric("Dimensione totale KB", round(float(inventory["size_kb"].sum()), 2))

    st.download_button(
        "Scarica report export Markdown",
        data=report_md.encode("utf-8"),
        file_name="maize_ncrna_gnn_export_report.md",
        mime="text/markdown",
    )

    st.download_button(
        "Scarica inventario export CSV",
        data=inventory.to_csv(index=False).encode("utf-8"),
        file_name="maize_ncrna_gnn_export_inventory.csv",
        mime="text/csv",
    )

    zip_bytes = build_project_zip_bytes(
        ranking_df=ranking_df,
        similarity_df=similarity_df,
        features_df=features_df,
        embeddings_df=embeddings_df,
        nodes_df=nodes_df,
        edges_df=edges_df,
    )

    st.download_button(
        "Scarica pacchetto ZIP progetto",
        data=zip_bytes,
        file_name="maize_ncrna_gnn_project_export.zip",
        mime="application/zip",
    )

    st.info(
        "Lo ZIP contiene tabelle, report e knowledge base RAG. "
        "Non contiene la cartella .venv né file temporanei."
    )

    st.divider()

    st.subheader("Tabelle principali")

    download_button_for_df(ranking_df, "Ranking esplorativo", "ncrna_exploratory_ranking.csv")
    download_button_for_df(similarity_df, "Similarità top-k", "ncrna_similarity_topk.csv")
    download_button_for_df(features_df, "Feature ncRNA", "ncrna_features.csv")
    download_button_for_df(embeddings_df, "Embedding ncRNA", "ncrna_embeddings.csv")
    download_button_for_df(nodes_df, "Graph nodes", "graph_nodes.csv")
    download_button_for_df(edges_df, "Graph edges", "graph_edges.csv")

    st.divider()

    st.subheader("Inventario file export")
    safe_dataframe(inventory, width="stretch", height=420)

    st.subheader("Inventario file principali GUI")
    safe_dataframe(file_status_table(), width="stretch", height=420)




def main() -> None:
    ranking_df = load_csv(PATHS["ranking"])
    similarity_df = load_csv(PATHS["similarity"])
    features_df = load_csv(PATHS["features"])
    embeddings_df = load_csv(PATHS["embeddings"])
    nodes_df = load_csv(PATHS["graph_nodes"])
    edges_df = load_csv(PATHS["graph_edges"])

    st.sidebar.title("🌽 ncRNA GNN")
    st.sidebar.caption("maize_drought_ncrna_gnn")

    page = st.sidebar.radio(
        "Navigazione",
        [
            "Dashboard",
            "Ranking candidati",
            "Dettaglio candidato",
            "Ricerca avanzata",
            "Similarità",
            "Embedding",
            "Chat RAG",
            "Report metodologico",
            "Download",
        ],
    )

    st.sidebar.divider()
    st.sidebar.caption("Modalità: esplorativa / non supervisionata")

    if page == "Dashboard":
        page_dashboard(ranking_df, similarity_df, features_df, nodes_df, edges_df)
    elif page == "Ranking candidati":
        page_ranking(ranking_df)
    elif page == "Dettaglio candidato":
        page_candidate_detail(ranking_df, similarity_df, features_df, embeddings_df)
    elif page == "Ricerca avanzata":
        page_advanced_search(ranking_df, similarity_df, features_df)
    elif page == "Similarità":
        page_similarity(similarity_df)
    elif page == "Embedding":
        page_embeddings(embeddings_df)
    elif page == "Chat RAG":
        page_rag_chat()
    elif page == "Report metodologico":
        page_reports()
    elif page == "Download":
        page_downloads(ranking_df, similarity_df, features_df, embeddings_df, nodes_df, edges_df)


if __name__ == "__main__":
    main()

# BEGIN V3_PROJECT_REALIZATION_BLOCK
def _render_v3_project_realization_block():
    """
    Blocco informativo aggiunto alla GUI V3.
    """
    try:
        import streamlit as st

        st.sidebar.markdown("---")
        st.sidebar.markdown("### Realizzazione progetto")
        st.sidebar.markdown("**Ing. Luigi Cocozza**")
        st.sidebar.markdown("Esperto in: **AI Machine Learning Engineer**")
        st.sidebar.markdown("_Ingegnere di Machine Learning IA_")
    except Exception:
        pass


_render_v3_project_realization_block()
# END V3_PROJECT_REALIZATION_BLOCK

