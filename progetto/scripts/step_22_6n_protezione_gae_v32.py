#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.2 - STEP 22.6N

Protezione preventiva del Graph Autoencoder originale.

Analizza i percorsi di output del trainer legacy
e crea una copia di sicurezza delle cartelle
originali prima di qualsiasi nuovo training.

Non:
- esegue training
- modifica GUI
- modifica checkpoint
- riavvia Docker
- modifica gli output scientifici
"""

import ast
import hashlib
import json
import shutil
import sys

from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

TRAINER = ROOT / "25_train_graph_autoencoder.py"

GUI = ROOT / "app/pages/01_Rete_neurale_GAE.py"

RUNS = (
    ROOT / "target_prediction_v33_vscode" / "runs"
)

MAX_BACKUP_BYTES = 400 * 1024 * 1024

OUTPUT_NAMES = (
    "OUT_DIR",
    "MODEL_DIR",
    "MODEL_PATH",
    "PLOTS_DIR",
    "EMB_DIR",
    "HISTORY_CSV",
    "NODE_EMB_CSV",
    "NCRNA_EMB_CSV",
)


def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def resolve_expression(node, values):

    if isinstance(node, ast.Name):

        return values.get(node.id)

    if isinstance(node, ast.Constant):

        if isinstance(node.value, str):
            return node.value

    if isinstance(node, ast.BinOp):

        if isinstance(node.op, ast.Div):

            left = resolve_expression(
                node.left,
                values
            )

            right = resolve_expression(
                node.right,
                values
            )

            if (
                isinstance(left, Path)
                and isinstance(right, str)
            ):
                return left / right

    return None


def inspect_assignments(tree):

    values = {
        "BASE_DIR": ROOT,
        "ROOT": ROOT,
    }

    expressions = {}

    for node in tree.body:

        if not isinstance(node, ast.Assign):
            continue

        for target in node.targets:

            if not isinstance(target, ast.Name):
                continue

            name = target.id

            expressions[name] = ast.unparse(
                node.value
            )

            if name == "BASE_DIR":
                continue

            value = resolve_expression(
                node.value,
                values
            )

            if value is not None:
                values[name] = value

    return values, expressions


def collect_writes(tree):

    names = {
        "save",
        "savefig",
        "to_csv",
        "to_excel",
        "write_text",
        "write_bytes",
        "dump",
        "savez",
        "savez_compressed",
    }

    results = []

    for node in ast.walk(tree):

        if not isinstance(node, ast.Call):
            continue

        function = ast.unparse(node.func)

        short = function.split(".")[-1]

        if short not in names:
            continue

        results.append({
            "line": node.lineno,
            "operation": function,
            "arguments": [
                ast.unparse(arg)[:140]
                for arg in node.args[:3]
            ],
        })

    return sorted(
        results,
        key=lambda x: x["line"]
    )


def inspect_gui():

    content = GUI.read_text(
        encoding="utf-8"
    )

    tree = ast.parse(content)

    commands = []

    for node in ast.walk(tree):

        if isinstance(node, ast.Call):

            text = ast.unparse(node.func)

            if text == "subprocess.run":

                commands.append(
                    node.lineno
                )

    return {
        "file": str(GUI.relative_to(ROOT)),
        "subprocess_lines": commands,
        "sha256": sha256(GUI),
    }


def safe_folder(path):

    if not isinstance(path, Path):
        return False

    resolved = path.resolve()

    # Non consentire backup di tutto il repository
    # o di directory generali troppo ampie.
    forbidden = {
        ROOT.resolve(),
        (ROOT / "outputs").resolve(),
        (ROOT / "models").resolve(),
        (ROOT / "data").resolve(),
    }

    if resolved in forbidden:
        return False

    return (
        resolved.is_relative_to(ROOT.resolve())
        and resolved != RUNS.resolve()
        and not RUNS.resolve().is_relative_to(resolved)
    )


def directory_inventory(folder):

    total = 0
    files = []

    for p in sorted(folder.rglob("*")):

        if p.is_symlink():

            raise RuntimeError(
                "Link simbolico presente: "
                + str(p)
            )

        if not p.is_file():
            continue

        size = p.stat().st_size

        total += size

        files.append({
            "name": str(p.relative_to(folder)),
            "size": size,
            "sha256": sha256(p),
        })

        if total > MAX_BACKUP_BYTES:

            raise RuntimeError(
                "Cartella superiore a 400 MB: "
                + str(folder)
                + ". Backup non avviato."
            )

    return total, files


def main():

    print("=" * 68)
    print("AI PRIMETTA - STEP 22.6N")
    print("PROTEZIONE GAE ORIGINALE")
    print("=" * 68)

    for path in (TRAINER, GUI):

        if not path.is_file():
            raise FileNotFoundError(path)

    source = TRAINER.read_text(
        encoding="utf-8"
    )

    tree = ast.parse(source)

    values, expressions = inspect_assignments(
        tree
    )

    print("\n[1] OUTPUT DEL TRAINER GAE")

    for name in OUTPUT_NAMES:

        if name not in expressions:
            continue

        print(
            name,
            "=",
            expressions[name]
        )

        value = values.get(name)

        if isinstance(value, Path):

            print(
                "  percorso:",
                value
            )

    print("\n[2] OPERAZIONI DI SCRITTURA")

    writes = collect_writes(tree)

    for item in writes:

        print(
            "Riga",
            item["line"],
            "|",
            item["operation"],
            "|",
            item["arguments"]
        )

    print("\n[3] IDENTIFICAZIONE CARTELLE DA PROTEGGERE")

    folders = {}

    for name in ("MODEL_DIR", "OUT_DIR"):

        path = values.get(name)

        if not isinstance(path, Path):

            raise RuntimeError(
                "Percorso non risolvibile: "
                + name
            )

        if not safe_folder(path):

            raise RuntimeError(
                "Percorso non sufficientemente "
                "specifico o non sicuro: "
                + str(path)
            )

        folders[name] = path

        print(
            name,
            path,
            "| esiste:",
            path.is_dir()
        )

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    destination = RUNS / (
        "STEP22_6N_BACKUP_GAE_V32_" + stamp
    )

    # Prima dell'inizio del backup, esaminiamo
    # le cartelle per escludere link e grandi quantità.
    inventories = {}

    for name, folder in folders.items():

        if folder.is_dir():

            size, files = directory_inventory(
                folder
            )

            inventories[name] = {
                "source": str(folder),
                "size_bytes": size,
                "files": files,
            }

        else:

            inventories[name] = {
                "source": str(folder),
                "size_bytes": 0,
                "files": [],
                "not_present": True,
            }

    destination.mkdir(
        parents=True,
        exist_ok=False
    )

    print("\n[4] CREAZIONE COPIE DI SICUREZZA")

    for name, folder in folders.items():

        if not folder.is_dir():

            print(
                "[ASSENTE]",
                name,
                "- nessun file da copiare"
            )

            continue

        target = destination / name

        shutil.copytree(
            folder,
            target,
            symlinks=False
        )

        _, copied = directory_inventory(
            target
        )

        if copied != inventories[name]["files"]:

            raise RuntimeError(
                "Verifica backup fallita: "
                + name
            )

        print(
            "[OK] Backup verificato:",
            name
        )

    report = {
        "step": "22.6N",
        "trainer": str(
            TRAINER.relative_to(ROOT)
        ),
        "trainer_sha256": sha256(TRAINER),
        "gui": inspect_gui(),
        "output_expressions": {
            k: v
            for k, v in expressions.items()
            if k in OUTPUT_NAMES
        },
        "writes": writes,
        "inventories": inventories,
        "training_performed": False,
        "original_files_modified": False,
    }

    report_path = (
        destination / "RAPPORTO_STEP22_6N.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print("\n[REPORT]", report_path)
    print("[OK] Nessun training")
    print("[OK] Originali preservati")
    print("[OK] Backup verificati")
    print("=" * 68)


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