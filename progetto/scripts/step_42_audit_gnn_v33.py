
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 | STEP 42
AUDIT CONSOLIDATO DELLA GNN ESISTENTE

--audit:
    Esamina file e risultati in sola lettura.

--run:
    Esegue lo stesso audit e salva un rapporto JSON.

NON ESEGUE:
    Training
    Inferenza
    Modifiche alla GUI
    Modifiche Docker
    Modifiche ai checkpoint
"""

import argparse
import csv
import json
import os
import re
import sys

from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

SKIP_DIRS = {
    ".git",
    ".venv",
    "venv",
    "env",
    "__pycache__",
    "node_modules",
    "site-packages",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "fastq",
    "mirnaseq_trimmed",
    "collapsed_srna_counts",
    "reference",
    "references",
    "downloads",
    "raw",
    "tmp",
    "cache"
}

MODEL_EXT = {
    ".pt",
    ".pth",
    ".ckpt",
    ".safetensors",
    ".onnx",
    ".joblib",
    ".pkl"
}

DATA_EXT = {
    ".csv",
    ".tsv",
    ".parquet",
    ".npz",
    ".npy",
    ".json",
    ".jsonl"
}

CODE_EXT = {
    ".py",
    ".r"
}

DOC_EXT = {
    ".json",
    ".yaml",
    ".yml",
    ".md",
    ".txt",
    ".toml"
}

TERMS = re.compile(
    r"gnn|graph|hetero|sage|gat|torch|pyg|"
    r"train|valid|eval|metric|split|model|predict",
    re.I
)

GRAPH_TERMS = re.compile(
    r"gnn|graph|hetero|sage|gat|"
    r"torch_geometric|pyg|message_passing",
    re.I
)

EVAL_TERMS = re.compile(
    r"evaluat|metric|valid|holdout|test[_-]|"
    r"cross.?valid|fold|auc|f1|rmse|accuracy",
    re.I
)

DATANAMES = re.compile(
    r"node|edge|graph|hetero|embedding|split|"
    r"fold|feature|label|metric|evaluat|"
    r"validat|testset",
    re.I
)

SIGNALS = {
    "gnn":
        r"torch_geometric|HeteroConv|HeteroData|"
        r"SAGEConv|GATConv|GCNConv|GraphSAGE|"
        r"MessagePassing|\bGNN\b",

    "group_split":
        r"GroupKFold|GroupShuffleSplit|"
        r"LeaveOneGroupOut|StratifiedGroupKFold|"
        r"group[_ ]?holdout|group[_ ]?split",

    "other_split":
        r"train_test_split|random_split|"
        r"train_mask|test_mask|val_mask|"
        r"train_idx|test_idx|split[_ ]?data",

    "evaluation":
        r"f1_score|roc_auc|average_precision|"
        r"confusion_matrix|mean_squared_error|"
        r"accuracy_score|\bRMSE\b|\bMAE\b|\bAUC\b",

    "feature_evidence":
        r"sample_id|run_accession|genotype|"
        r"condition|edge_index|node_type|edge_type"
}

MAX_FILES = 30000
MAX_READ_BYTES = 750000


# ============================================================
# 2. PERCORSI RELATIVI
# ============================================================

def rel(path):
    return str(path.relative_to(ROOT))


# ============================================================
# 3. INDIVIDUAZIONE RAPPORTI SCIENTIFICI
# ============================================================

def last_report(pattern):

    found = (
        sorted(RUNS.glob(pattern))
        if RUNS.is_dir()
        else []
    )

    return found[-1] if found else None


# ============================================================
# 4. VERIFICA SCIENTIFICA STEP30 / STEP41
# ============================================================

def past_science():

    facts = {
        "step30": None,
        "step41": None,
        "critical_rule":
            "Candidate sRNA e coppie RNAplex esplorative: "
            "non usarle come etichette positive validate."
    }

    p30 = last_report(
        "STEP30_REGISTRO_SCIENTIFICO_*/"
        "RAPPORTO_STEP30.json"
    )

    if p30:

        try:
            obj = json.loads(
                p30.read_text(encoding="utf-8")
            )

            facts["step30"] = {
                "path": rel(p30),
                "step": obj.get("step"),
                "candidate_count": obj.get(
                    "candidate_count"
                ),
                "note":
                    "Lo stato scientifico va confermato "
                    "nel rapporto."
            }

        except (OSError, ValueError) as exc:
            facts["step30"] = {
                "path": rel(p30),
                "error": str(exc)[:160]
            }

    p41 = last_report(
        "STEP41_EDGER_DROUGHT_*/"
        "CONTRASTI_2_SRNA.tsv"
    )

    if p41:

        try:

            with p41.open(
                encoding="utf-8-sig",
                newline=""
            ) as stream:

                rows = list(
                    csv.DictReader(
                        stream,
                        delimiter="\t"
                    )
                )

            needed = {
                "srna_id",
                "contrasto",
                "log2FC",
                "p_value",
                "FDR_BH"
            }

            if (
                len(rows) != 6
                or not needed.issubset(
                    set(rows[0] if rows else [])
                )
            ):
                raise ValueError(
                    "Tabella STEP41 incompleta: "
                    "attese sei righe"
                )

            detail = []

            for row in rows:

                q = float(row["FDR_BH"])

                if not 0 <= q <= 1:
                    raise ValueError(
                        "FDR fuori intervallo"
                    )

                detail.append({
                    "id": row["srna_id"],
                    "contrast": row["contrasto"],
                    "fdr_bh": q,
                    "log2fc": float(
                        row["log2FC"]
                    )
                })

            facts["step41"] = {
                "path": rel(p41),
                "n_comparisons": len(detail),
                "n_significant_fdr_005": sum(
                    item["fdr_bh"] < 0.05
                    for item in detail
                ),
                "details": detail
            }

        except (
            OSError,
            ValueError,
            TypeError,
            KeyError
        ) as exc:

            facts["step41"] = {
                "path": rel(p41),
                "error": str(exc)[:160]
            }

    return facts


# ============================================================
# 5. INVENTARIO DEL PROGETTO
# ============================================================

def walk_project():

    files = []
    truncated = False

    for folder, dirs, filenames in os.walk(
        ROOT,
        followlinks=False
    ):

        parent = Path(folder)

        # I rapporti STEP30 e STEP41 vengono
        # verificati separatamente.
        # Nei runs analizziamo soltanto
        # directory riconducibili a GNN/ML.

        if parent == RUNS:

            dirs[:] = sorted(
                d for d in dirs

                if not d.startswith(
                    "STEP42_AUDIT_GNN_"
                )

                and re.search(
                    r"gnn|graph|hetero|model|train|"
                    r"machine[_-]?learn|[_-]ml[_-]",
                    d,
                    re.I
                )
            )

            continue

        dirs[:] = sorted(
            d for d in dirs

            if d not in SKIP_DIRS

            and not (
                parent / d
            ).is_symlink()
        )

        for name in sorted(filenames):

            path = parent / name

            if (
                path.is_symlink()
                or not path.is_file()
                or path.resolve()
                == Path(__file__).resolve()
            ):
                continue

            suffix = path.suffix.lower()

            if suffix not in (
                MODEL_EXT
                | CODE_EXT
                | DATA_EXT
                | DOC_EXT
            ):
                continue

            files.append(path)

            if len(files) >= MAX_FILES:
                truncated = True
                return files, truncated

    return files, truncated


# ============================================================
# 6. ANALISI DEL CODICE GNN
# ============================================================

def inspect_code(path):

    try:

        if path.stat().st_size > MAX_READ_BYTES:

            return {
                "path": rel(path),
                "note":
                    "File grande, non letto",
                "signals": []
            }

        content = path.read_text(
            encoding="utf-8",
            errors="replace"
        )

    except OSError as exc:

        return {
            "path": rel(path),
            "error": str(exc)[:120],
            "signals": []
        }

    flags = [
        key

        for key, regex in SIGNALS.items()

        if re.search(
            regex,
            content,
            re.I
        )
    ]

    return {
        "path": rel(path),
        "signals": flags,
        "bytes": len(
            content.encode("utf-8")
        )
    }


# ============================================================
# 7. LETTURA LIMITATA DEI MANIFEST JSON
# ============================================================

def inspect_manifest(path):

    item = {
        "path": rel(path)
    }

    try:

        if path.stat().st_size > 2000000:

            item["note"] = (
                "JSON >2 MB, non letto"
            )

            return item

        obj = json.loads(
            path.read_text(
                encoding="utf-8"
            )
        )

        if not isinstance(obj, dict):

            item["note"] = (
                "JSON non oggetto"
            )

            return item

        keys = list(obj)

        item["keys"] = keys[:22]

        indicator = re.compile(
            r"accuracy|f1|auc|precision|recall|"
            r"rmse|mae|r2|valid|test|split|"
            r"seed|checkpoint|train|n_samples",
            re.I
        )

        item["evidence_keys"] = [
            k for k in keys

            if indicator.search(k)
        ][:18]

    except (
        OSError,
        ValueError,
        TypeError
    ) as exc:

        item["error"] = (
            str(exc)[:120]
        )

    return item


# ============================================================
# 8. CLASSIFICAZIONE DELL'INVENTARIO
# ============================================================

def inventory(paths, truncated):

    codes = []
    weights = []
    datasets = []
    manifests = []

    for path in paths:

        suffix = path.suffix.lower()
        relative = rel(path)
        filename = path.name

        # ====================================================
        # 9. CHECKPOINT E MODELLI
        # ====================================================

        if suffix in MODEL_EXT:

            weights.append({
                "path": relative,
                "bytes": path.stat().st_size,
                "potentially_gnn": bool(
                    GRAPH_TERMS.search(relative)
                )
            })

        # ====================================================
        # 10. CODICE PYTHON / R
        # ====================================================

        elif (
            suffix in CODE_EXT
            and TERMS.search(filename)
        ):

            if re.match(
                r"(?:audit_|step_\d+_audit_)",
                filename,
                re.I
            ):
                continue

            entry = inspect_code(path)

            if (
                GRAPH_TERMS.search(relative)
                or "gnn" in entry["signals"]
                or "group_split" in entry["signals"]
                or "other_split" in entry["signals"]
            ):

                codes.append(entry)

        # ====================================================
        # 11. DATASET E FILE STRUTTURATI
        # ====================================================

        elif (
            suffix in DATA_EXT
            and DATANAMES.search(filename)
            and TERMS.search(relative)
        ):

            datasets.append({
                "path": relative,
                "bytes": path.stat().st_size
            })

        # ====================================================
        # 12. METRICHE E MANIFEST
        # ====================================================

        if (
            suffix == ".json"
            and (
                GRAPH_TERMS.search(relative)
                or EVAL_TERMS.search(filename)
            )
            and (
                "manifest" in filename.lower()
                or "metric" in filename.lower()
                or "eval" in filename.lower()
                or "report" in filename.lower()
            )
        ):

            manifests.append(
                inspect_manifest(path)
            )

    # ========================================================
    # 13. PRIORITA' DEL CODICE DA CONTROLLARE
    # ========================================================

    def codescore(item):

        flags = set(
            item.get("signals", [])
        )

        return (
            10 * ("gnn" in flags)
            + 5 * ("group_split" in flags)
            + 3 * ("evaluation" in flags)
            + (
                "train" in item["path"].lower()
            )
        )

    codes.sort(
        key=lambda item: (
            -codescore(item),
            item["path"]
        )
    )

    weights.sort(
        key=lambda item: (
            -int(item["potentially_gnn"]),
            item["path"]
        )
    )

    datasets.sort(
        key=lambda item: item["path"]
    )

    manifests.sort(
        key=lambda item: item["path"]
    )

    gnn_codes = [
        item for item in codes

        if (
            "gnn" in item.get("signals", [])
            or GRAPH_TERMS.search(item["path"])
        )
    ]

    return {
        "scanned_relevant_files": len(paths),

        "truncated_at_limit": truncated,

        "counts": {
            "candidate_code": len(codes),
            "candidate_gnn_code": len(gnn_codes),
            "model_artifacts": len(weights),
            "graph_data_candidates": len(datasets),
            "metadata_reports": len(manifests)
        },

        "code": codes[:30],

        "checkpoints": weights[:30],

        "graph_data": datasets[:30],

        "reports": manifests[:25],

        "code_total": len(codes),

        "checkpoints_total": len(weights),

        "has_group_split_signal": any(
            "group_split" in item.get(
                "signals", []
            )
            for item in codes
        ),

        "has_evaluation_signal": any(
            "evaluation" in item.get(
                "signals", []
            )
            for item in codes
        ),

        "has_gnn_code_signal": (
            len(gnn_codes) > 0
        )
    }


# ============================================================
# 14. STAMPA RAPPORTO
# ============================================================

def display(facts, inv):

    print("=" * 72)

    print(
        "AI PRIMETTA V3.3 | STEP 42 "
        "- AUDIT CONSOLIDATO GNN"
    )

    print("=" * 72)

    step41 = facts["step41"]

    if (
        step41
        and "n_significant_fdr_005" in step41
    ):

        print(
            "[PASS] STEP41:",
            step41["n_comparisons"],
            "confronti delle 2 sRNA;",
            "FDR<0.05 =",
            step41["n_significant_fdr_005"]
        )

    else:

        print(
            "[ATTENZIONE] Tabella STEP41 "
            "non verificabile:",
            step41
        )

    print(
        "[STEP30]",
        "presente"
        if facts["step30"]
        else "non trovato"
    )

    print(
        "[INVENTARIO] File tecnici esaminati:",
        inv["scanned_relevant_files"],
        "(LIMITE RAGGIUNTO)"
        if inv["truncated_at_limit"]
        else ""
    )

    for key, value in inv["counts"].items():

        print(
            f"[CONTEGGIO] {key}={value}"
        )

    # ========================================================
    # 15. FILE PIU' RILEVANTI
    # ========================================================

    sections = (
        ("CODICE", "code", 12),
        ("CHECKPOINT", "checkpoints", 12),
        ("DATI_GRAFO", "graph_data", 8),
        ("REPORT", "reports", 8)
    )

    for label, key, limit in sections:

        items = inv[key]

        print(
            f"[{label}] Prime "
            f"{min(limit, len(items))} evidenze:"
        )

        for item in items[:limit]:

            if key == "code":

                note = ",".join(
                    item.get("signals", [])
                )

            else:

                note = ",".join(
                    item.get("evidence_keys", [])
                )

            print(
                " -",
                item["path"],
                "| " + note if note else ""
            )

    print(
        "[CHECK] Codice di GNN individuato:",
        "SI"
        if inv["has_gnn_code_signal"]
        else "NON CONFERMATO"
    )

    print(
        "[CHECK] Segnali di split per gruppi:",
        "SI"
        if inv["has_group_split_signal"]
        else "NON CONFERMATI"
    )

    print(
        "[CHECK] Segnali di metriche nel codice:",
        "SI"
        if inv["has_evaluation_signal"]
        else "NON CONFERMATI"
    )

    print(
        "[LIMITE] Nomi e segnali nel codice "
        "NON dimostrano training, test indipendenti "
        "o assenza di leakage."
    )

    print(
        "[VINCOLO] Le 2 sRNA e le 36 interazioni "
        "restano candidate esplorative; "
        "nessuna etichetta positiva validata."
    )


# ============================================================
# 16. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    modes = parser.add_mutually_exclusive_group(
        required=True
    )

    modes.add_argument(
        "--audit",
        action="store_true",
        help="Sola lettura"
    )

    modes.add_argument(
        "--run",
        action="store_true",
        help="Salva un rapporto audit JSON"
    )

    args = parser.parse_args()

    if not ROOT.is_dir() or not RUNS.is_dir():

        raise FileNotFoundError(
            "Eseguire all'interno del progetto "
            "AI PRIMETTA con "
            "target_prediction_v33_vscode/runs"
        )

    # ========================================================
    # 17. AUDIT COMPLETO
    # ========================================================

    facts = past_science()

    files, truncated = walk_project()

    inv = inventory(
        files,
        truncated
    )

    display(
        facts,
        inv
    )

    # ========================================================
    # 18. SALVATAGGIO OPZIONALE DEL RAPPORTO
    # ========================================================

    if args.run:

        stamp = datetime.now(
            timezone.utc
        ).strftime(
            "%Y%m%dT%H%M%S%fZ"
        )

        out = RUNS / (
            f"STEP42_AUDIT_GNN_{stamp}"
        )

        out.mkdir(
            parents=False,
            exist_ok=False
        )

        path = (
            out / "RAPPORTO_STEP42.json"
        )

        payload = {
            "step": 42,

            "type":
                "static_inventory_not_model_validation",

            "created_utc": stamp,

            "science": facts,

            "gnn_inventory": inv,

            "limitations": [
                "Nessun modello caricato",
                "Nessun checkpoint eseguito",
                "Nessun test di generalizzazione",
                "Nessun training",
                "Nessuna modifica a GUI o Docker"
            ]
        }

        path.write_text(
            json.dumps(
                payload,
                ensure_ascii=False,
                indent=2
            ) + "\n",
            encoding="utf-8"
        )

        print(
            "[REPORT]",
            path
        )

    print(
        "[OK] Audit completato; sorgenti, GUI, "
        "Docker, modelli e dati originali invariati"
    )


# ============================================================
# 19. AVVIO
# ============================================================

if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print(
            f"[ERRORE] "
            f"{type(exc).__name__}: {exc}",
            file=sys.stderr
        )

        raise SystemExit(1)
