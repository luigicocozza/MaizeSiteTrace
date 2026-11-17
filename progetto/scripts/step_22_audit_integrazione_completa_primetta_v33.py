#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 22

COLLAUDO CONSOLIDATO:
- GUI V3.2
- GUI V3.3
- ML/GNN
- Docker
- Dati scientifici
- Collegamento applicativo

Non modifica GUI, API, checkpoint, dati,
immagini o container persistenti.

Crea soltanto un rapporto nella cartella runs.
"""

import ast
import csv
import hashlib
import json
import re
import subprocess
import sys
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

BASE = ROOT / "target_prediction_v33_vscode"
RUNS = BASE / "runs"
EVID = RUNS / "STEP20_EVIDENZE_GUI"

GUI_V32 = ROOT / "app/streamlit_app.py"
GUI_V33 = BASE / "app_v33.py"

MODEL = (
    ROOT
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

TABLE = (
    ROOT
    / "outputs/prjna419326/"
    "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)

IMAGE = "ai_primetta_v33_ml_compat:19.5"

V33_CONTAINER = "ai_primetta_target_v33-prediction-1"

IDS = [
    "sRNA_00000002",
    "sRNA_00000040",
    "sRNA_00000007",
    "sRNA_00000435",
    "sRNA_00000037",
]

EXCLUDE = {
    "runs",
    "__pycache__",
    ".venv",
    "venv",
    ".git",
    "delivery",
    "node_modules",
    "backups",
}


# ==========================================================
# TEST ML NEL CONTAINER COMPATIBILE
# ==========================================================

MODEL_TEST = r'''
import json
import re
import warnings

from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn

from sklearn.exceptions import InconsistentVersionWarning


root = Path("/project")

model_path = (
    root
    / "models/srna_supervised_ml_v1_antileakage/"
    "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
)

table_path = (
    root
    / "outputs/prjna419326/"
    "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
)

ids = [
    "sRNA_00000002",
    "sRNA_00000040",
    "sRNA_00000007",
    "sRNA_00000435",
    "sRNA_00000037",
]

result = {
    "status": "NON_VERIFICATO",
    "sklearn_version": sklearn.__version__,
    "records": [],
    "independent_validation": False,
    "gui_invoked_this_model": False,
}

try:

    if sklearn.__version__ != "1.8.0":
        raise RuntimeError(
            "Versione scikit-learn diversa da 1.8.0"
        )

    # Solo checkpoint locale attendibile.
    with warnings.catch_warnings():

        warnings.simplefilter(
            "error",
            InconsistentVersionWarning
        )

        model = joblib.load(model_path)

    names = list(
        getattr(model, "feature_names_in_", [])
    )

    if len(names) != 22:
        raise RuntimeError(
            "Schema a 22 feature non rilevato"
        )

    df = pd.read_csv(
        table_path,
        sep="\t"
    )

    missing = [
        name for name in names
        if name not in df.columns
    ]

    if missing:
        raise RuntimeError(
            "Colonne mancanti: " + str(missing)
        )

    idcols = [
        col
        for col in (
            "rna_id",
            "srna_id",
            "short_name",
            "feature",
            "id",
        )
        if col in df.columns
    ]

    if not idcols:
        raise RuntimeError(
            "Nessuna colonna identificativa RNA trovata"
        )

    selected = []

    for target in ids:

        matches = df.iloc[0:0]

        for col in idcols:

            pattern = (
                re.escape(target)
                + r"(?![0-9])"
            )

            candidate = df.loc[
                df[col]
                .astype(str)
                .str.contains(
                    pattern,
                    regex=True,
                    na=False
                )
            ]

            if len(candidate) == 1:
                matches = candidate
                break

        if len(matches) != 1:
            raise RuntimeError(
                f"{target}: trovate {len(matches)} "
                "righe univoche"
            )

        selected.append(
            matches.iloc[0]
        )

    selected_df = pd.DataFrame(selected)

    X = selected_df.loc[:, names].copy()

    for col in names:
        X[col] = pd.to_numeric(
            X[col],
            errors="coerce"
        )

    predictions = model.predict(X)

    probabilities = model.predict_proba(X)

    if not np.isfinite(probabilities).all():
        raise RuntimeError(
            "Probabilita non finite"
        )

    classes = [
        str(c)
        for c in model.classes_
    ]

    positive_idx = (
        classes.index("1")
        if "1" in classes
        else None
    )

    for i, target in enumerate(ids):

        historical = selected_df.iloc[i].get(
            "antileakage_predicted_label"
        )

        value = predictions[i]

        result["records"].append({
            "rna": target,
            "class_predicted": (
                str(int(value))
                if str(value) in ("0", "1")
                else str(value)
            ),
            "probability_class_1": (
                float(
                    probabilities[i, positive_idx]
                )
                if positive_idx is not None
                else None
            ),
            "historical_label": (
                None
                if pd.isna(historical)
                else str(historical)
            ),
        })

    result["status"] = (
        "PASS_5_RNA_DATI_GIA_PRESENTI"
    )

    result["note"] = (
        "Input ricavati dalla tabella ML storica, "
        "non da nuovo caricamento GUI "
        "ne da dataset indipendente"
    )

except Exception as exc:

    result["error"] = (
        type(exc).__name__
        + ": "
        + str(exc)
    )

print(
    "MODEL_PROBE_JSON="
    + json.dumps(
        result,
        ensure_ascii=False
    ),
    flush=True
)
'''


# ==========================================================
# FUNZIONI DI SUPPORTO
# ==========================================================

def call(args, timeout=30):

    try:

        process = subprocess.run(
            args,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False
        )

        return {
            "ok": process.returncode == 0,
            "code": process.returncode,
            "stdout": process.stdout.strip(),
            "stderr": process.stderr.strip()[-1500:]
        }

    except Exception as exc:

        return {
            "ok": False,
            "error": str(exc)
        }


def sha256(path):

    if not path.is_file():
        return None

    digest = hashlib.sha256()

    with path.open("rb") as handle:

        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b""
        ):
            digest.update(chunk)

    return digest.hexdigest()


def last_report(prefix, filename):

    matches = sorted(
        RUNS.glob(
            f"{prefix}_*/{filename}"
        )
    )

    if not matches:
        return {"status": "MANCANTE"}

    try:

        data = json.loads(
            matches[-1].read_text(
                encoding="utf-8"
            )
        )

        return {
            "status": data.get("status"),
            "path": str(matches[-1])
        }

    except Exception as exc:

        return {
            "status": "ERRORE",
            "error": str(exc)
        }


# ==========================================================
# ANALISI CODICE GUI E IMPORT
# ==========================================================

def local_module(name, current):

    paths = [
        current.parent,
        ROOT,
        BASE,
        ROOT / "app",
        ROOT / "api",
        ROOT / "src",
    ]

    for parent in paths:

        relative = Path(
            *name.split(".")
        )

        candidates = (
            parent / f"{relative}.py",
            parent / relative / "__init__.py",
        )

        for path in candidates:

            if (
                path.is_file()
                and path.resolve().is_relative_to(
                    ROOT.resolve()
                )
            ):
                return path.resolve()

    return None


def entry_audit(entry):

    if not entry.is_file():
        return {
            "exists": False,
            "source_files": [],
            "model_calls": [],
            "remote_calls": [],
            "imports": [],
        }

    visited = set()
    queue = [entry.resolve()]

    model_calls = []
    remote_calls = []
    import_strings = []

    while queue and len(visited) < 100:

        path = queue.pop()

        if path in visited:
            continue

        if not path.is_file():
            continue

        if any(
            part in EXCLUDE
            for part in path.parts
        ):
            continue

        visited.add(path)

        try:

            source = path.read_text(
                encoding="utf-8"
            )

            tree = ast.parse(
                source,
                filename=str(path)
            )

        except (
            OSError,
            UnicodeError,
            SyntaxError
        ):
            continue

        lines = source.splitlines()

        relative = str(
            path.relative_to(ROOT)
        )

        modules = []

        for node in ast.walk(tree):

            if isinstance(node, ast.Import):

                modules.extend(
                    alias.name
                    for alias in node.names
                )

            elif isinstance(node, ast.ImportFrom):

                if node.module:

                    modules.append(
                        node.module
                    )

                    modules.extend(
                        node.module + "." + alias.name
                        for alias in node.names
                    )

                elif node.level:

                    modules.extend(
                        alias.name
                        for alias in node.names
                    )

            elif isinstance(node, ast.Call):

                try:
                    function = ast.unparse(
                        node.func
                    )
                except Exception:
                    function = ""

                text = lines[
                    node.lineno - 1
                ].strip()[:220]

                model_signals = (
                    "torch.load",
                    "joblib.load",
                    "pickle.load",
                    "load_state_dict",
                    "predict_proba",
                    ".predict",
                    ".encode",
                    "torch_geometric",
                )

                if any(
                    token in function.lower()
                    for token in model_signals
                ):

                    model_calls.append({
                        "file": relative,
                        "line": node.lineno,
                        "call": function[:100],
                        "source": text,
                    })

                remote_signals = (
                    "requests.get",
                    "requests.post",
                    "httpx.",
                    "urlopen",
                    "subprocess.run",
                    "subprocess.Popen",
                )

                if any(
                    token in function
                    for token in remote_signals
                ):

                    remote_calls.append({
                        "file": relative,
                        "line": node.lineno,
                        "call": function[:100],
                        "source": text,
                    })

                if (
                    "import_module" in function
                    and node.args
                    and isinstance(
                        node.args[0],
                        ast.Constant
                    )
                    and isinstance(
                        node.args[0].value,
                        str
                    )
                ):

                    modules.append(
                        node.args[0].value
                    )

        import_strings.extend(
            {
                "file": relative,
                "module": module
            }
            for module in modules[:60]
        )

        for module in modules:

            resolved = local_module(
                module,
                path
            )

            if (
                resolved is not None
                and resolved not in visited
            ):
                queue.append(resolved)

    return {
        "exists": True,

        "source_files": sorted(
            str(path.relative_to(ROOT))
            for path in visited
        ),

        "model_calls": model_calls,

        "remote_calls": remote_calls,

        "imports": import_strings[:150],

        "source_limit_reached": (
            len(visited) >= 100
        ),
    }


# ==========================================================
# DOCKER E SERVIZI
# ==========================================================

def inspect_container(name):

    response = call([
        "docker",
        "inspect",
        name
    ])

    if not response["ok"]:

        return {
            "found": False,
            "error": (
                response.get("stderr")
                or response.get("error")
            )
        }

    info = json.loads(
        response["stdout"]
    )[0]

    return {
        "found": True,

        "running": bool(
            info.get(
                "State", {}
            ).get("Running")
        ),

        "image": info.get(
            "Config", {}
        ).get("Image"),

        "cmd": info.get(
            "Config", {}
        ).get("Cmd"),

        "working_dir": info.get(
            "Config", {}
        ).get("WorkingDir"),

        "ports": info.get(
            "NetworkSettings", {}
        ).get("Ports"),

        "mounts": [
            {
                "source": item.get("Source"),
                "destination": item.get("Destination"),
                "rw": item.get("RW")
            }
            for item in info.get(
                "Mounts", []
            )
        ]
    }


def health(url):

    try:

        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        )

        with opener.open(
            url,
            timeout=6
        ) as response:

            body = response.read(
                2000
            ).decode(
                "utf-8",
                errors="replace"
            )

            return {
                "ok": response.status == 200,
                "http_code": response.status,
                "body": (
                    body[:50]
                    if url.endswith(
                        "/_stcore/health"
                    )
                    else None
                ),
            }

    except Exception as exc:

        return {
            "ok": False,
            "error": str(exc)
        }


def dependencies_v33():

    code = (
        'import importlib.util,json;'
        'print(json.dumps({'
        'k:bool(importlib.util.find_spec(k)) '
        'for k in ['
        '"torch","torch_geometric","joblib",'
        '"sklearn","streamlit","requests","httpx"'
        ']}))'
    )

    response = call([
        "docker",
        "exec",
        V33_CONTAINER,
        "python3",
        "-B",
        "-c",
        code
    ])

    if not response["ok"]:

        return {
            "status": "NON_VERIFICATO",
            "error": (
                response.get("stderr")
                or response.get("error")
            )
        }

    try:

        return {
            "status": "OK",
            "modules": json.loads(
                response["stdout"]
            )
        }

    except ValueError:

        return {
            "status": "NON_VERIFICATO",
            "detail": response["stdout"][:300]
        }


# ==========================================================
# EVIDENZE GUI
# ==========================================================

def evidence():

    expected = [
        "01_rna_qc.tsv",
        "02_alignment_hypotheses.tsv",
        "03_historical_12mer_hypotheses.tsv",
        "04_gene_target_candidates.tsv",
        "05_rna_candidate_summary.tsv",
        "06_provenance_manifest.json",
    ]

    found = {}

    for name in expected:

        path = EVID / name

        item = {
            "found": path.is_file(),
            "sha256": sha256(path)
        }

        if path.is_file():

            try:

                if path.suffix == ".tsv":

                    with path.open(
                        encoding="utf-8-sig",
                        newline=""
                    ) as handle:

                        rows = list(
                            csv.DictReader(
                                handle,
                                delimiter="\t"
                            )
                        )

                    item["row_count"] = len(rows)

                    if name in (
                        "01_rna_qc.tsv",
                        "05_rna_candidate_summary.tsv"
                    ):

                        item["rna_ids"] = [
                            row.get("rna_id", "")
                            for row in rows
                        ]

                        if name == "01_rna_qc.tsv":

                            item["eligible"] = sum(
                                str(
                                    row.get("status")
                                ).lower() == "eligible"
                                for row in rows
                            )

                        else:

                            item["new_sites"] = sum(
                                int(
                                    row["n_new_sites"]
                                )
                                for row in rows
                            )

                else:

                    obj = json.loads(
                        path.read_text(
                            encoding="utf-8-sig"
                        )
                    )

                    item["json_type"] = (
                        type(obj).__name__
                    )

            except Exception as exc:

                item["error"] = str(exc)

        found[name] = item

    return found


# ==========================================================
# INFERENZA ML SU 5 SMALL RNA REALI
# ==========================================================

def five_srna_model_probe():

    if (
        not MODEL.is_file()
        or not TABLE.is_file()
    ):

        return {
            "status": "NON_VERIFICATO",
            "reason": (
                "Checkpoint o tabella ML assenti"
            )
        }

    if not call([
        "docker",
        "image",
        "inspect",
        IMAGE
    ])["ok"]:

        return {
            "status": "NON_VERIFICATO",
            "reason": (
                "Immagine Docker compatibile assente"
            )
        }

    command = [
        "docker",
        "run",
        "--rm",
        "--pull", "never",
        "--network", "none",
        "--read-only",
        "--tmpfs",
        "/tmp:rw,nosuid,nodev,size=128m",
        "--cap-drop", "ALL",
        "--security-opt",
        "no-new-privileges",
        "--pids-limit", "128",
        "--memory", "3g",
        "--cpus", "2",
        "--user", "65534:65534",
        "-e", "HOME=/tmp",
        "-v", f"{ROOT}:/project:ro",
        "-w", "/project",
        "--entrypoint", "python3",
        IMAGE,
        "-B",
        "-c",
        MODEL_TEST,
    ]

    response = call(
        command,
        timeout=180
    )

    if not response["ok"]:

        return {
            "status": "NON_VERIFICATO",
            "error": (
                response.get("stderr")
                or response.get("error")
            ),
            "stdout_tail": response.get(
                "stdout", ""
            )[-1000:]
        }

    marker = "MODEL_PROBE_JSON="

    for line in response[
        "stdout"
    ].splitlines():

        if line.startswith(marker):

            try:
                return json.loads(
                    line[len(marker):]
                )

            except ValueError:
                break

    return {
        "status": "NON_VERIFICATO",
        "error": (
            "Nessun risultato strutturato"
        )
    }


# ==========================================================
# PROGRAMMA PRINCIPALE
# ==========================================================

def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3 - STEP 22")
    print("COLLAUDO COMPLETO / AUDIT INTEGRAZIONE")
    print("=" * 72)

    if not GUI_V33.is_file():

        raise SystemExit(
            "[ERRORE] app_v33.py non trovato"
        )

    RUNS.mkdir(
        parents=True,
        exist_ok=True
    )

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = RUNS / (
        f"STEP22_AUDIT_INTEGRAZIONE_{stamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report_path = (
        output / "RAPPORTO_STEP22.json"
    )

    # ------------------------------------------------------
    # 1. REPORT PRECEDENTI
    # ------------------------------------------------------

    print("\n[1] REPORT GIÀ COLLAUDATI")

    reports = {

        "gae": last_report(
            "STEP19_3_INFERENZA_GAE_V3",
            "RAPPORTO_STEP19_3.json"
        ),

        "ml": last_report(
            "STEP19_5_COMPATIBILITA_ML",
            "RAPPORTO_STEP19_5.json"
        ),

        "finale": last_report(
            "STEP20_RAPPORTO_FINALE",
            "RAPPORTO_FINALE_STEP20.json"
        ),
    }

    for name, record in reports.items():
        print(
            name.upper(),
            record["status"]
        )

    # ------------------------------------------------------
    # 2. EVIDENZE GUI
    # ------------------------------------------------------

    print("\n[2] EVIDENZE DELLA GUI PREDITTIVA")

    exported = evidence()

    for name, item in exported.items():

        print(
            "[OK]"
            if item["found"]
            and "error" not in item
            else "[DA VERIFICARE]",
            name,
            "righe=",
            item.get("row_count", "-")
        )

    # ------------------------------------------------------
    # 3. ANALISI COLLEGAMENTI NEL CODICE
    # ------------------------------------------------------

    print("\n[3] ANALISI CODICE GUI E MODULI IMPORTATI")

    sources = {
        "v32": entry_audit(GUI_V32),
        "v33": entry_audit(GUI_V33)
    }

    for name, item in sources.items():

        print(
            name,
            "| file:",
            len(item["source_files"]),
            "| chiamate ML/GAE:",
            len(item["model_calls"]),
            "| chiamate esterne:",
            len(item["remote_calls"])
        )

        for callsite in item[
            "model_calls"
        ][:12]:

            print(
                "  ",
                callsite["file"],
                ":",
                callsite["line"],
                callsite["call"]
            )

    # ------------------------------------------------------
    # 4. DOCKER
    # ------------------------------------------------------

    print("\n[4] CONTAINER E HEALTH CHECK")

    container = inspect_container(
        V33_CONTAINER
    )

    print(
        "GUI predizione container running:",
        container.get("running")
    )

    checks = {

        "v33_health": health(
            "http://127.0.0.1:8502/_stcore/health"
        ),

        "v33_page": health(
            "http://127.0.0.1:8502/"
        )
    }

    print(
        "V3.3 health:",
        checks["v33_health"].get(
            "http_code",
            checks["v33_health"].get("error")
        )
    )

    if container.get("running"):

        deps = dependencies_v33()

    else:

        deps = {
            "status": "NON_VERIFICATO"
        }

    print(
        "Librerie nel container V3.3:",
        deps.get(
            "modules",
            deps.get("status")
        )
    )

    # ------------------------------------------------------
    # 5. PAGINE GUI ORIGINALE
    # ------------------------------------------------------

    print("\n[5] PAGINE GUI ORIGINALE V3.2")

    original_pages = []

    for page in sorted(
        (ROOT / "app/pages").glob("*.py")
    ):

        info = entry_audit(page)

        original_pages.append({
            "file": str(
                page.relative_to(ROOT)
            ),

            "model_call_count": len(
                info["model_calls"]
            ),

            "model_calls": info[
                "model_calls"
            ][:20]
        })

    print(
        "Pagine individuate:",
        len(original_pages)
    )

    for item in original_pages:

        print(
            item["file"],
            "| indizi ML/GNN:",
            item["model_call_count"]
        )

    # Controllo del servizio originale,
    # senza avviare container.

    compose_path = ROOT / "docker-compose.yml"

    if compose_path.is_file():

        response = call([
            "docker",
            "compose",
            "-f",
            str(compose_path),
            "ps",
            "-q",
            "maize-ncrna-gnn"
        ])

        ids = (
            response.get(
                "stdout", ""
            ).splitlines()
            if response["ok"]
            else []
        )

    else:
        ids = []

    v32_container = (
        inspect_container(ids[0])
        if ids
        else {
            "found": False,
            "running": False
        }
    )

    print(
        "Container GUI V3.2 running:",
        v32_container.get("running")
    )

    # ------------------------------------------------------
    # 6. INFERENZA ML SULLE 5 sRNA
    # ------------------------------------------------------

    print(
        "\n[6] INFERENZA ML SULLE "
        "5 SMALL RNA REALI"
    )

    probe = five_srna_model_probe()

    print(
        "Esito:",
        probe.get("status")
    )

    if probe.get("error"):
        print(
            "[ATTENZIONE]",
            probe["error"]
        )

    for item in probe.get(
        "records", []
    ):

        probability = item.get(
            "probability_class_1"
        )

        print(
            item["rna"],
            "classe=",
            item["class_predicted"],
            "p(classe1)=",
            (
                round(probability, 5)
                if probability is not None
                else "-"
            )
        )

    # ------------------------------------------------------
    # 7. VALUTAZIONE DEL COLLEGAMENTO
    # ------------------------------------------------------

    print("\n[7] VALUTAZIONE INTEGRAZIONE")

    evidence_calls = sources[
        "v33"
    ]["model_calls"]

    possible_ml = any(

        any(
            signal in item["call"].lower()

            for signal in (
                "torch.load",
                "joblib.load",
                "predict_proba",
                ".predict",
                ".encode",
            )
        )

        for item in evidence_calls
    )

    if possible_ml:

        status = (
            "INDIZI_COLLEGAMENTO_NEL_CODICE_"
            "DA_VERIFICARE"
        )

    else:

        status = (
            "COLLEGAMENTO_DIRETTO_GUI_V33_"
            "ML_GNN_NON_DIMOSTRATO"
        )

    limitations = [
        (
            "Non e dimostrato che il click "
            "della GUI V3.3 richiami i modelli "
            "GAE e ML."
        ),

        (
            "Le chiamate presenti nel codice "
            "non certificano da sole il percorso "
            "di esecuzione della GUI."
        ),

        (
            "La prova ML sulle cinque RNA usa "
            "feature storiche gia disponibili, "
            "non nuovi dati indipendenti."
        ),

        (
            "La GUI originale V3.2 necessita "
            "di una verifica dal vivo qualora "
            "il container non sia attivo."
        ),

        (
            "Non e stata svolta una nuova "
            "validazione biologica indipendente."
        )
    ]

    print("Stato:", status)

    # ------------------------------------------------------
    # 8. RAPPORTO FINALE
    # ------------------------------------------------------

    print("\n[8] GENERAZIONE RAPPORTO")

    report = {

        "project": "AI PRIMETTA V3.3",

        "step": "22",

        "timestamp_utc": stamp,

        "status": status,

        "previous_reports": reports,

        "evidence": exported,

        "gui_source_audit": sources,

        "v32_pages_audit": original_pages,

        "v33_container": container,

        "v32_container": v32_container,

        "v33_http": checks,

        "v33_installed_modules": deps,

        "five_real_srna_ml_probe": probe,

        "gui_button_to_gae_ml_runtime_verified": False,

        "independent_scientific_validation": False,

        "training_performed": False,

        "existing_components_modified": False,

        "limitations": limitations,

        "hash_originals": {

            str(path.relative_to(ROOT)): sha256(path)

            for path in (
                MODEL,
                TABLE,
                GUI_V32,
                GUI_V33,
                ROOT
                / "models/graph_autoencoder_v3/"
                "gae_v3_model.pt",
            )
        }
    }

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print()
    print("=" * 72)
    print("STEP 22 TERMINATO")
    print("=" * 72)

    print("Stato:", status)
    print("Report:", report_path)

    print(
        "[OK] Nessun training"
    )

    print(
        "[OK] GUI originali non modificate"
    )

    print(
        "[OK] Checkpoint invariati"
    )

    print(
        "[OK] Container persistenti non modificati"
    )

    print(
        "\n[NOTA] Il collegamento completo richiede "
        "una prova runtime del percorso "
        "GUI -> modello ML/GNN."
    )


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__,
            str(exc),
            file=sys.stderr
        )

        sys.exit(1)