from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parents[2]

LINK_SCRIPT = ROOT / "27_link_user_data_to_ranking.py"
OUT_ROOT = ROOT / "outputs/user_data"


st.set_page_config(
    page_title="Ranking dati ricercatore",
    layout="wide",
)


def latest_output_dir() -> Path | None:
    if not OUT_ROOT.exists():
        return None

    dirs = [p for p in OUT_ROOT.iterdir() if p.is_dir()]

    if not dirs:
        return None

    dirs = sorted(dirs, key=lambda p: p.stat().st_mtime, reverse=True)
    return dirs[0]


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def download_button(label: str, path: Path, mime: str) -> None:
    if path.exists():
        st.download_button(
            label=label,
            data=path.read_bytes(),
            file_name=path.name,
            mime=mime,
        )


st.title("Ranking con evidenza del ricercatore")

st.warning(
    "Modulo sperimentale. Il collegamento tra ranking e dati utente e' tecnico-esplorativo. "
    "La presenza di un ncRNA nel dataset del ricercatore non dimostra da sola un ruolo biologico nella siccita'."
)

st.write(
    "Questa pagina mostra quali ncRNA del ranking esplorativo sono presenti anche "
    "nei dati caricati dal ricercatore."
)

st.divider()

st.subheader("Aggiorna collegamento ranking-dati utente")

if st.button("Aggiorna ranking con ultimo dataset caricato", type="primary"):
    if not LINK_SCRIPT.exists():
        st.error(f"Script non trovato: {LINK_SCRIPT}")
    else:
        cmd = [
            sys.executable,
            str(LINK_SCRIPT),
        ]

        st.code(" ".join(cmd), language="bash")

        with st.spinner("Collegamento ranking-dati utente in corso..."):
            proc = subprocess.run(
                cmd,
                cwd=str(ROOT),
                capture_output=True,
                text=True,
                timeout=600,
            )

        if proc.returncode == 0:
            st.success("Collegamento completato.")
            st.code(proc.stdout[-12000:], language="text")
        else:
            st.error("Errore durante il collegamento.")
            st.code(proc.stdout[-8000:], language="text")
            st.code(proc.stderr[-8000:], language="text")


st.divider()

latest_dir = latest_output_dir()

if latest_dir is None:
    st.info("Nessun output utente trovato. Caricare prima un dataset nella pagina Dati ricercatore.")
    st.stop()


summary_path = latest_dir / "user_ranking_integration_summary.txt"
report_path = latest_dir / "user_ranking_integration_report.json"
ranking_path = latest_dir / "ranking_with_user_evidence.csv"
matched_path = latest_dir / "matched_user_ncrna.csv"
unmatched_path = latest_dir / "unmatched_user_ncrna.csv"
raw_path = latest_dir / "user_evidence_raw_records.csv"

summary = read_text(summary_path)
report = read_json(report_path)
ranking_df = read_csv(ranking_path)
matched_df = read_csv(matched_path)
unmatched_df = read_csv(unmatched_path)
raw_df = read_csv(raw_path)

st.subheader("Ultimo dataset integrato")

st.write(f"Cartella output: {latest_dir}")

if report:
    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric("Righe ranking", report.get("ranking_rows", "-"))

    with col2:
        st.metric("Record utente", report.get("user_raw_records", "-"))

    with col3:
        st.metric("Match nel ranking", report.get("matched_ranking_rows", "-"))

    with col4:
        st.metric("ID utente non trovati", report.get("unmatched_user_ids", "-"))


if summary:
    with st.expander("Summary integrazione", expanded=False):
        st.text(summary)


st.divider()

st.subheader("ncRNA del ranking supportati dai dati del ricercatore")

if matched_df.empty:
    st.warning("Nessun ncRNA del ranking e' stato trovato nei dati caricati.")
else:
    df = matched_df.copy()

    search_id = st.text_input("Filtra per rnacentral_id", value="")

    if search_id.strip():
        mask = df["rnacentral_id"].astype(str).str.contains(search_id.strip(), case=False, na=False)
        df = df[mask]

    only_fdr = st.checkbox("Mostra solo FDR <= 0.05", value=False)

    if only_fdr and "user_min_FDR" in df.columns:
        fdr = pd.to_numeric(df["user_min_FDR"], errors="coerce")
        df = df[fdr <= 0.05]

    preferred_cols = [
        "exploratory_rank",
        "rnacentral_id",
        "rna_type",
        "exploratory_rank_score",
        "user_data_present",
        "user_records_count",
        "user_best_log2FC_abs",
        "user_min_FDR",
        "user_min_pvalue",
        "user_conditions",
        "user_source_files",
        "methodological_note",
    ]

    existing_cols = [c for c in preferred_cols if c in df.columns]

    st.dataframe(
        df[existing_cols] if existing_cols else df,
        use_container_width=True,
    )

    st.write(f"Candidati visualizzati: {len(df)}")


st.divider()

st.subheader("Record grezzi caricati dal ricercatore")

if raw_df.empty:
    st.info("Nessun record grezzo disponibile.")
else:
    st.dataframe(raw_df, use_container_width=True)


st.divider()

st.subheader("ID caricati ma non trovati nel ranking")

if unmatched_df.empty:
    st.success("Tutti gli ID utente sono stati trovati nel ranking.")
else:
    st.dataframe(unmatched_df, use_container_width=True)


st.divider()

st.subheader("Download risultati")

col_a, col_b, col_c = st.columns(3)

with col_a:
    download_button("Scarica ranking con evidenza utente", ranking_path, "text/csv")
    download_button("Scarica candidati matched", matched_path, "text/csv")

with col_b:
    download_button("Scarica record grezzi utente", raw_path, "text/csv")
    download_button("Scarica unmatched", unmatched_path, "text/csv")

with col_c:
    download_button("Scarica summary", summary_path, "text/plain")
    download_button("Scarica report JSON", report_path, "application/json")


st.divider()

st.info(
    "Interpretazione: questi risultati indicano coincidenza tecnica tra ranking e dati caricati. "
    "Per parlare di evidenza biologica servono controlli sperimentali, qualita' RNA-seq, replicati, "
    "analisi differenziale robusta e validazione indipendente."
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

