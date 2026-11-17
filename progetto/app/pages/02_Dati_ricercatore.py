from __future__ import annotations

import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parents[2]
UPLOAD_ROOT = ROOT / "data/user_data/uploads"
REPORT_ROOT = ROOT / "outputs/user_data"
VALIDATOR = ROOT / "26_validate_user_dataset.py"


st.set_page_config(
    page_title="Dati ricercatore",
    layout="wide",
)


def safe_name(name: str) -> str:
    name = name.strip().replace(" ", "_")
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    return name or "user_dataset"


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def preview_file(path: Path) -> None:
    suffix = path.suffix.lower()

    if suffix in [".csv", ".tsv", ".txt"]:
        try:
            if suffix == ".tsv":
                df = pd.read_csv(path, sep="\t")
            elif suffix == ".csv":
                df = pd.read_csv(path)
            else:
                df = pd.read_csv(path, sep=None, engine="python")

            st.dataframe(df.head(30), use_container_width=True)

        except Exception as exc:
            st.warning(f"Anteprima tabella non disponibile: {exc}")

    elif suffix in [".fa", ".fasta", ".fna"]:
        try:
            lines = []
            with path.open("r", encoding="utf-8", errors="replace") as f:
                for i, line in enumerate(f):
                    if i >= 40:
                        break
                    lines.append(line.rstrip())

            st.code("\n".join(lines), language="text")

        except Exception as exc:
            st.warning(f"Anteprima FASTA non disponibile: {exc}")


st.title("Dati ricercatore - caricamento dataset")

st.warning(
    "Modulo sperimentale. Questa pagina carica e valida dati gia' processati "
    "del ricercatore. Non esegue ancora una pipeline completa da FASTQ grezzi."
)

st.write("Il ricercatore puo' caricare file FASTA, CSV o TSV gia' processati.")

st.subheader("Formato consigliato per espressione differenziale")

st.code(
    "rnacentral_id,gene_id,condition,log2FC,pvalue,FDR\n"
    "URS0003C75441,Zm00001eb000000,drought_vs_control,1.82,0.0004,0.012\n"
    "URS0003C7543B,Zm00001eb000001,drought_vs_control,-1.10,0.0030,0.041",
    language="text",
)

st.divider()

dataset_name = st.text_input(
    "Nome dataset / esperimento",
    value="esperimento_ricercatore_01",
)

uploaded_files = st.file_uploader(
    "Carica file CSV, TSV, FASTA",
    type=["csv", "tsv", "fa", "fasta", "fna"],
    accept_multiple_files=True,
)

dataset_safe = safe_name(dataset_name)
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

if "last_dataset_dir" not in st.session_state:
    st.session_state["last_dataset_dir"] = None


if st.button("Salva e valida dataset", type="primary"):
    if not uploaded_files:
        st.error("Caricare almeno un file.")

    else:
        UPLOAD_ROOT.mkdir(parents=True, exist_ok=True)

        dataset_dir = UPLOAD_ROOT / f"{dataset_safe}_{timestamp}"
        dataset_dir.mkdir(parents=True, exist_ok=True)

        for uploaded_file in uploaded_files:
            filename = safe_name(uploaded_file.name)
            out_path = dataset_dir / filename
            out_path.write_bytes(uploaded_file.getbuffer())

        st.session_state["last_dataset_dir"] = str(dataset_dir)

        st.success(f"File salvati in: {dataset_dir}")

        if not VALIDATOR.exists():
            st.error(f"Script validazione non trovato: {VALIDATOR}")

        else:
            cmd = [
                sys.executable,
                str(VALIDATOR),
                "--input-dir",
                str(dataset_dir),
                "--dataset-name",
                dataset_dir.name,
            ]

            st.code(" ".join(cmd), language="bash")

            with st.spinner("Validazione dataset in corso..."):
                proc = subprocess.run(
                    cmd,
                    cwd=str(ROOT),
                    capture_output=True,
                    text=True,
                    timeout=600,
                )

            if proc.returncode == 0:
                st.success("Validazione completata.")
                st.code(proc.stdout[-10000:], language="text")

            else:
                st.error("Validazione con errori.")
                st.code(proc.stdout[-8000:], language="text")
                st.code(proc.stderr[-8000:], language="text")


st.divider()

st.subheader("Ultimo dataset caricato")

last_dir = st.session_state.get("last_dataset_dir")

if last_dir:
    dataset_dir = Path(last_dir)

    st.write(f"Cartella dataset: {dataset_dir}")

    files = sorted(dataset_dir.glob("*"))

    if files:
        st.write("File caricati:")

        for path in files:
            st.write(f"- {path.name} - {path.stat().st_size} byte")

        selected = st.selectbox(
            "Anteprima file",
            options=[path.name for path in files],
        )

        selected_path = dataset_dir / selected
        preview_file(selected_path)

    report_dir = REPORT_ROOT / dataset_dir.name
    report_json = report_dir / "user_dataset_validation_report.json"
    summary_txt = report_dir / "user_dataset_validation_summary.txt"

    st.divider()
    st.subheader("Report validazione")

    summary = read_text(summary_txt)
    report = read_json(report_json)

    if summary:
        st.text(summary)

    if report:
        with st.expander("Report JSON completo"):
            st.json(report)

        col1, col2 = st.columns(2)

        with col1:
            if summary_txt.exists():
                st.download_button(
                    "Scarica summary validazione",
                    data=summary_txt.read_bytes(),
                    file_name=summary_txt.name,
                    mime="text/plain",
                )

        with col2:
            if report_json.exists():
                st.download_button(
                    "Scarica report JSON",
                    data=report_json.read_bytes(),
                    file_name=report_json.name,
                    mime="application/json",
                )

else:
    st.info("Nessun dataset caricato in questa sessione.")


st.divider()

st.info(
    "Prossimo step: collegare i dati caricati al ranking, aggiornare il RAG "
    "e permettere alla Chat AI di rispondere anche sui dati del ricercatore."
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

