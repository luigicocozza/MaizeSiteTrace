#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.1 - VERIFICA MAPPING RNA-TRASCRITTO-GENE

Diagnosi non distruttiva.
Non modifica GUI, dataset, modelli o motore.
"""

import ast
import gzip
import json
import os
import re
from pathlib import Path

ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"

SCRIPT16 = (
    ROOT / "scripts"
    / "step_16_test_5_srna_reali_v33.py"
)

ENGINE = BASE / "app_v33.py"

TRANSCRIPTOME = (
    BASE / "references"
    / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz"
)

SKIP = {
    ".git", ".venv", "venv",
    "__pycache__", "node_modules",
}


def open_text(path):
    if path.suffix == ".gz":
        return gzip.open(
            path, "rt", encoding="utf-8",
            errors="replace"
        )
    return path.open(
        "r", encoding="utf-8",
        errors="replace"
    )


def transcript_ids(path):
    result = set()

    if not path.exists():
        return result

    with open_text(path) as f:
        for line in f:
            if line.startswith(">"):
                result.add(
                    line[1:].split()[0].split(".")[0]
                )

    return result


def find_mapping_files():
    results = []

    for folder, dirs, files in os.walk(ROOT):
        dirs[:] = [
            d for d in dirs if d not in SKIP
        ]

        for name in files:
            low = name.lower()

            if (
                "tx2gene" in low
                or "transcript_to_gene" in low
                or "transcript2gene" in low
            ):
                path = Path(folder) / name

                if path.is_file() and (
                    low.endswith(".tsv")
                    or low.endswith(".csv")
                    or low.endswith(".txt")
                    or low.endswith(".tsv.gz")
                ):
                    results.append(path)

    return sorted(set(results))


def inspect_mapping(path, known_ids):
    total = 0
    matched = set()
    examples = []
    seen_ids = set()

    try:
        with open_text(path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue

                total += 1

                delimiter = (
                    "\t" if "\t" in line else ","
                )

                columns = [
                    c.strip().strip('"')
                    for c in line.strip().split(delimiter)
                ]

                if len(examples) < 3:
                    examples.append(columns[:3])

                if len(columns) >= 2:
                    tx = columns[0].split(".")[0]
                    seen_ids.add(tx)

                    if tx in known_ids:
                        matched.add(tx)

        return {
            "file": str(path),
            "righe": total,
            "trascritti_distinti_colonna_1": len(seen_ids),
            "coincidenze_con_trascrittoma": len(matched),
            "copertura_trascritti_percento": round(
                100 * len(matched) / len(known_ids), 3
            ) if known_ids else None,
            "esempi": examples,
            "nota": (
                "Controllo esplorativo: si presume "
                "che la prima colonna contenga ID trascritti."
            )
        }

    except Exception as exc:
        return {
            "file": str(path),
            "errore": str(exc)
        }


def inspect_script(path, terms):
    print(f"\nFILE: {path}")

    if not path.exists():
        print("[ATTENZIONE] File non trovato")
        return

    lines = path.read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()

    selected = set()

    for i, line in enumerate(lines):
        if any(term in line.lower() for term in terms):
            for j in range(
                max(0, i - 2),
                min(len(lines), i + 3)
            ):
                selected.add(j)

    for i in sorted(selected):
        print(f"{i + 1:5d} | {lines[i][:180]}")


def inspect_engine():
    print("\nFIRMA DEL MOTORE V3.3")

    if not ENGINE.exists():
        print("[ATTENZIONE] Motore non trovato")
        return

    tree = ast.parse(
        ENGINE.read_text(
            encoding="utf-8",
            errors="replace"
        )
    )

    for node in ast.walk(tree):
        if isinstance(
            node, (ast.FunctionDef, ast.AsyncFunctionDef)
        ) and node.name == "analyze":

            args = (
                node.args.posonlyargs
                + node.args.args
                + node.args.kwonlyargs
            )

            print(
                "analyze(" +
                ", ".join(a.arg for a in args) +
                ")"
            )


def main():
    print("=" * 68)
    print("AI PRIMETTA V3.3 - STEP 17.1")
    print("VERIFICA MAPPING REALE RNA-TRASCRITTO-GENE")
    print("=" * 68)

    print("\n[1] TRASCRITTOMA")
    print("Percorso:", TRANSCRIPTOME)

    ids = transcript_ids(TRANSCRIPTOME)
    print("Trascritti distinti normalizzati:", len(ids))

    print("\n[2] TUTTI I MAPPING INDIVIDUATI")

    mappings = find_mapping_files()
    results = []

    for path in mappings:
        info = inspect_mapping(path, ids)
        results.append(info)

        print("\n--------------------------------")
        print("File:", info["file"])
        print("Righe:", info.get("righe", "N/D"))
        print(
            "Coincidenze:",
            info.get("coincidenze_con_trascrittoma", "N/D")
        )
        print(
            "Copertura:",
            info.get("copertura_trascritti_percento", "N/D")
        )
        print("Esempi:", info.get("esempi", []))

        if "errore" in info:
            print("Errore:", info["errore"])

    if not mappings:
        print("[ATTENZIONE] Nessun mapping individuato")

    print("\n[3] PARAMETRI E PERCORSI DELLO STEP 16")

    inspect_script(
        SCRIPT16,
        [
            "mapping",
            "tx2gene",
            "transcript",
            "docker",
            "analyze",
        ]
    )

    print("\n[4] MOTORE V3.3")
    inspect_engine()

    print("\n[5] SALVATAGGIO REPORT")

    output = BASE / "reports"
    output.mkdir(parents=True, exist_ok=True)

    destination = (
        output / "STEP17_1_AUDIT_MAPPING.json"
    )

    report = {
        "step": "17.1",
        "trascritti_normalizzati": len(ids),
        "mapping_trovati": results,
        "step16_script": str(SCRIPT16),
        "motore": str(ENGINE),
        "mapping_effettivo_step16_verificato": False,
    }

    destination.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print("[OK] Report:", destination)
    print("\n[OK] STEP 17.1 TERMINATO")
    print("[OK] GUI e modelli invariati")
    print(
        "[NOTA] Verificare l'effettivo mapping "
        "usato nello STEP 16 dall'output dei percorsi."
    )


if __name__ == "__main__":
    main()