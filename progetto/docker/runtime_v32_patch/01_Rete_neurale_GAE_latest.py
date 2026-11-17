from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import streamlit as st


ROOT = Path(__file__).resolve().parents[2]

SCRIPT_PATH = ROOT / "25_train_graph_autoencoder.py"
SUMMARY_PATH = ROOT / "outputs/gae/gae_summary.txt"
METRICS_PATH = ROOT / "outputs/gae/gae_metrics.json"
HISTORY_PATH = ROOT / "outputs/gae/gae_training_history.csv"

MODEL_PATH = ROOT / "models/graph_autoencoder/gae_model.pt"

NODE_EMB_PATH = ROOT / "data/processed/embeddings/gae_node_embeddings.csv"
NCRNA_EMB_PATH = ROOT / "data/processed/embeddings/gae_ncrna_embeddings.csv"

PLOTS_DIR = ROOT / "outputs/gae/plots"


st.set_page_config(
    page_title="Rete neurale GAE",
    layout="wide",
)


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def file_download(label: str, path: Path, mime: str = "application/octet-stream") -> None:
    if path.exists():
        st.download_button(
            label=label,
            data=path.read_bytes(),
            file_name=path.name,
            mime=mime,
        )


st.title("Rete neurale sperimentale — Graph Autoencoder")

st.warning(
    "Modulo sperimentale: questa rete neurale apprende embedding dal grafo "
    "ricostruendo collegamenti. Non è un predittore supervisionato della risposta alla siccità."
)

st.markdown(
    """
Questa pagina permette di eseguire una **Graph Autoencoder non supervisionata** sulla rete già costruita.

La rete usa:

- feature dei nodi;
- struttura del grafo;
- archi del grafo;
- encoder GCN;
- decoder a prodotto interno.

Output principali:

- modello `.pt`;
- embedding neurali dei nodi;
- embedding neurali degli ncRNA;
- metriche AUC/AP di ricostruzione collegamenti;
- grafici di training.
"""
)

with st.expander("Parametri training", expanded=True):
    col1, col2, col3 = st.columns(3)

    with col1:
        epochs = st.number_input("Epoche", min_value=10, max_value=2000, value=300, step=10)
        hidden_dim = st.number_input("Hidden dimension", min_value=8, max_value=512, value=64, step=8)

    with col2:
        latent_dim = st.number_input("Latent dimension", min_value=2, max_value=128, value=16, step=2)
        lr = st.number_input("Learning rate", min_value=0.0001, max_value=1.0, value=0.01, step=0.001, format="%.4f")

    with col3:
        dropout = st.number_input("Dropout", min_value=0.0, max_value=0.9, value=0.20, step=0.05)
        seed = st.number_input("Seed", min_value=1, max_value=999999, value=42, step=1)


st.divider()

st.subheader("Avvio rete neurale")

if not SCRIPT_PATH.exists():
    st.error(f"Script non trovato: {SCRIPT_PATH}")
else:
    if st.button("Avvia training Graph Autoencoder", type="primary"):
        cmd = [
            sys.executable,
            str(SCRIPT_PATH),
            "--epochs",
            str(int(epochs)),
            "--hidden-dim",
            str(int(hidden_dim)),
            "--latent-dim",
            str(int(latent_dim)),
            "--lr",
            str(float(lr)),
            "--dropout",
            str(float(dropout)),
            "--seed",
            str(int(seed)),
        ]

        st.info("Training avviato. Attendere il completamento.")
        st.code(" ".join(cmd), language="bash")

        with st.spinner("Training Graph Autoencoder in corso..."):
            proc = subprocess.run(
                cmd,
                cwd=str(ROOT),
                capture_output=True,
                text=True,
                timeout=3600,
            )

        if proc.returncode == 0:
            st.success("Training completato correttamente.")
            st.code(proc.stdout[-12000:], language="text")
        else:
            st.error("Errore durante il training.")
            st.code(proc.stdout[-8000:], language="text")
            st.code(proc.stderr[-8000:], language="text")


st.divider()


# STEP22_6Q_BEGIN: risultati dell'ultimo training isolato

from pathlib import Path as _v32q_Path

_v32q_root = ROOT / "outputs" / "gae_gui_safe_runs"

_v32q_required = (
    "models/gae_model.pt",
    "embeddings/gae_node_embeddings.csv",
    "embeddings/gae_ncrna_embeddings.csv",
    "results/gae_training_history.csv",
    "results/gae_metrics.json",
)

_v32q_finished = sorted(
    (
        p for p in _v32q_root.iterdir()
        if p.is_dir()
        and not p.is_symlink()
        and all(
            (p / n).is_file()
            and (p / n).stat().st_size > 0
            for n in _v32q_required
        )
    ),
    key=lambda p: p.name,
) if _v32q_root.is_dir() else []

_v32q_run = (
    _v32q_finished[-1]
    if _v32q_finished
    else None
)


def _v32q_relocate(value):
    """
    Reindirizza i percorsi legacy GAE
    alla cartella dell'ultimo training completo.
    """

    if isinstance(value, _v32q_Path):

        path = value
        as_string = False

    elif (
        isinstance(value, str)
        and value.startswith("/")
    ):

        path = _v32q_Path(value)
        as_string = True

    elif isinstance(value, list):

        return [
            _v32q_relocate(v)
            for v in value
        ]

    elif isinstance(value, tuple):

        return tuple(
            _v32q_relocate(v)
            for v in value
        )

    elif isinstance(value, dict):

        return {
            k: _v32q_relocate(v)
            for k, v in value.items()
        }

    else:

        return value

    parts = path.parts

    alternatives = (
        (
            ("outputs", "gae"),
            "results"
        ),
        (
            ("models", "graph_autoencoder"),
            "models"
        ),
        (
            ("data", "processed", "embeddings"),
            "embeddings"
        ),
    )

    for prefix, directory in alternatives:

        for start in range(
            0,
            len(parts) - len(prefix) + 1
        ):

            if (
                tuple(
                    parts[start:start + len(prefix)]
                )
                != prefix
            ):
                continue

            suffix = parts[
                start + len(prefix):
            ]

            new_path = (
                _v32q_run / directory
            )

            for component in suffix:

                new_path = (
                    new_path / component
                )

            return (
                str(new_path)
                if as_string
                else new_path
            )

    return value


if _v32q_run is not None:

    _v32q_keywords = (
        "DIR",
        "PATH",
        "FILE",
        "PLOT",
        "SUMMARY",
        "METRIC",
        "MODEL",
        "EMB",
        "HISTORY",
        "CSV",
        "JSON",
        "TXT",
        "OUTPUT",
    )

    for _v32q_key, _v32q_value in list(
        globals().items()
    ):

        if (
            _v32q_key.startswith("_")
            or _v32q_key == "ROOT"
            or not any(
                k in _v32q_key.upper()
                for k in _v32q_keywords
            )
        ):
            continue

        _v32q_new = _v32q_relocate(
            _v32q_value
        )

        if (
            _v32q_new is not _v32q_value
            and _v32q_new != _v32q_value
        ):

            globals()[
                _v32q_key
            ] = _v32q_new

# STEP22_6Q_END

st.subheader("Stato modulo GAE")

summary = read_text(SUMMARY_PATH)
metrics = load_json(METRICS_PATH)

if summary:
    st.success("Risultati GAE disponibili.")
    st.text(summary)
else:
    st.info("Nessun riepilogo GAE disponibile. Avviare il training con il pulsante sopra.")


if metrics:
    st.subheader("Metriche")

    rows = []
    metric_block = metrics.get("metrics", {})

    for split_name in ["train", "validation", "test"]:
        split = metric_block.get(split_name, {})
        rows.append(
            {
                "split": split_name,
                "AUC": split.get("auc"),
                "Average precision": split.get("average_precision"),
            }
        )

    st.dataframe(pd.DataFrame(rows), use_container_width=True)

    with st.expander("JSON metriche completo"):
        st.json(metrics)


st.divider()

st.subheader("Grafici GAE")

plot_files = [
    PLOTS_DIR / "gae_loss_curve.png",
    PLOTS_DIR / "gae_validation_ap_curve.png",
    PLOTS_DIR / "gae_node_embeddings_2d.png",
    PLOTS_DIR / "gae_ncrna_embeddings_2d.png",
]

existing_plots = [p for p in plot_files if p.exists()]

if existing_plots:
    cols = st.columns(2)
    for i, plot in enumerate(existing_plots):
        with cols[i % 2]:
            st.image(str(plot), caption=plot.name, use_container_width=True)
else:
    st.info("Nessun grafico disponibile. Avviare prima il training.")


st.divider()

st.subheader("Download output GAE")

col_a, col_b, col_c = st.columns(3)

with col_a:
    file_download("Scarica summary GAE", SUMMARY_PATH, "text/plain")
    file_download("Scarica metriche JSON", METRICS_PATH, "application/json")
    file_download("Scarica storico training CSV", HISTORY_PATH, "text/csv")

with col_b:
    file_download("Scarica embedding nodi CSV", NODE_EMB_PATH, "text/csv")
    file_download("Scarica embedding ncRNA CSV", NCRNA_EMB_PATH, "text/csv")

with col_c:
    file_download("Scarica modello GAE", MODEL_PATH, "application/octet-stream")


st.divider()

st.info(
    "Interpretazione corretta: la rete neurale produce embedding neurali esplorativi. "
    "Non dimostra che un ncRNA sia coinvolto biologicamente nella risposta alla siccità."
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

