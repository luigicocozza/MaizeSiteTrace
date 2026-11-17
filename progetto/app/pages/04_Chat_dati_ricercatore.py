from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import streamlit as st


ROOT = Path(__file__).resolve().parents[2]

CHAT_SCRIPT = ROOT / "29_user_rag_chat.py"
BUILD_KB_SCRIPT = ROOT / "28_build_user_rag_knowledge.py"


st.set_page_config(
    page_title="Chat dati ricercatore",
    layout="wide",
)


st.title("Chat sui dati del ricercatore")

st.warning(
    "Modulo sperimentale. La chat risponde usando i dati caricati dal ricercatore "
    "e il ranking integrato. Le risposte sono esplorative e non dimostrano da sole "
    "un ruolo biologico nella risposta alla siccita'."
)

st.write(
    "Questa pagina permette di interrogare i dati caricati dal ricercatore, "
    "i candidati matched e il ranking integrato."
)

st.divider()

st.subheader("Aggiorna knowledge base RAG utente")

if st.button("Rigenera knowledge base dati ricercatore"):
    if not BUILD_KB_SCRIPT.exists():
        st.error(f"Script non trovato: {BUILD_KB_SCRIPT}")
    else:
        cmd = [
            sys.executable,
            str(BUILD_KB_SCRIPT),
        ]

        st.code(" ".join(cmd), language="bash")

        with st.spinner("Rigenerazione knowledge base utente in corso..."):
            proc = subprocess.run(
                cmd,
                cwd=str(ROOT),
                capture_output=True,
                text=True,
                timeout=600,
            )

        if proc.returncode == 0:
            st.success("Knowledge base utente rigenerata.")
            st.code(proc.stdout[-12000:], language="text")
        else:
            st.error("Errore durante la rigenerazione della knowledge base.")
            st.code(proc.stdout[-8000:], language="text")
            st.code(proc.stderr[-8000:], language="text")


st.divider()

st.subheader("Domanda alla Chat dati ricercatore")

examples = [
    "Quali ncRNA caricati dal ricercatore sono presenti nel ranking e hanno FDR inferiore a 0.05?",
    "Quali ncRNA del ricercatore hanno log2FC positivo?",
    "Quali candidati sono supportati dai dati caricati?",
    "Quali ncRNA hanno il miglior FDR nel dataset del ricercatore?",
    "Riassumi l'integrazione tra ranking e dati del ricercatore.",
]

example_choice = st.selectbox(
    "Esempi di domanda",
    options=[""] + examples,
)

default_question = example_choice if example_choice else examples[0]

question = st.text_area(
    "Scrivi la domanda",
    value=default_question,
    height=120,
)

top_k = st.number_input(
    "Numero massimo di chunk RAG",
    min_value=1,
    max_value=20,
    value=5,
    step=1,
)

if st.button("Interroga Chat dati ricercatore", type="primary"):
    if not CHAT_SCRIPT.exists():
        st.error(f"Script chat non trovato: {CHAT_SCRIPT}")

    elif not question.strip():
        st.error("Inserire una domanda.")

    else:
        cmd = [
            sys.executable,
            str(CHAT_SCRIPT),
            "--question",
            question.strip(),
            "--top-k",
            str(int(top_k)),
        ]

        st.code(" ".join(cmd), language="bash")

        with st.spinner("Risposta in elaborazione..."):
            proc = subprocess.run(
                cmd,
                cwd=str(ROOT),
                capture_output=True,
                text=True,
                timeout=600,
            )

        if proc.returncode == 0:
            st.success("Risposta generata.")
            st.text(proc.stdout)
        else:
            st.error("Errore nella generazione della risposta.")
            st.code(proc.stdout[-8000:], language="text")
            st.code(proc.stderr[-8000:], language="text")


st.divider()

st.info(
    "Interpretazione corretta: la chat usa dati caricati, ranking integrato e knowledge base utente. "
    "Per conclusioni biologiche servono controllo qualita', replicati, analisi differenziale robusta "
    "e validazione sperimentale."
)

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

