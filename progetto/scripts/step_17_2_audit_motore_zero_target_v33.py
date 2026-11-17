#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.2 - AUDIT INTERNO MOTORE ZERO TARGET

Analizza:
- Funzione analyze()
- Funzioni di complementarita
- Filtri e soglie
- Orientamento delle sequenze
- Gestione dei mismatch
- Generazione siti candidati
- Risultati del precedente STEP 16

NON modifica GUI, modelli, dati o Docker.
"""

import ast
import json
import re
from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"

ENGINE = BASE / "app_v33.py"

RUNS = BASE / "runs"

TERMS = [
    "mismatch",
    "complement",
    "reverse",
    "revcomp",
    "entropy",
    "candidate",
    "target",
    "site",
    "seed",
    "alignment",
    "align",
    "scan",
    "threshold",
    "maximum_sites",
    "max_mismatches",
    "rna_sequence",
    "transcript_sequence",
    "pair",
    "binding",
]


def find_functions(tree):

    functions = []

    for node in ast.walk(tree):

        if isinstance(
            node,
            (ast.FunctionDef, ast.AsyncFunctionDef)
        ):

            functions.append({
                "name": node.name,
                "line_start": node.lineno,
                "line_end": node.end_lineno,
                "parameters": [
                    a.arg
                    for a in (
                        node.args.posonlyargs
                        + node.args.args
                        + node.args.kwonlyargs
                    )
                ],
            })

    return sorted(
        functions,
        key=lambda x: x["line_start"]
    )


def relevant_lines(lines):

    selected = set()

    for index, line in enumerate(lines):

        low = line.lower()

        if any(term in low for term in TERMS):

            for position in range(
                max(0, index - 3),
                min(len(lines), index + 4)
            ):

                selected.add(position)

    return [
        {
            "line": index + 1,
            "code": lines[index]
        }
        for index in sorted(selected)
    ]


def analyze_function(lines, functions):

    function = next(
        (
            item
            for item in functions
            if item["name"] == "analyze"
        ),
        None
    )

    if function is None:
        return None

    start = function["line_start"] - 1
    end = function["line_end"]

    return {
        "parameters": function["parameters"],
        "start": function["line_start"],
        "end": function["line_end"],
        "source": [
            f"{index + 1:5d} | {lines[index]}"
            for index in range(start, end)
        ]
    }


def last_step16():

    if not RUNS.exists():
        return None

    directories = sorted(
        [
            p for p in RUNS.iterdir()
            if p.is_dir()
            and p.name.startswith(
                "STEP16_5_SRNA_REALI_"
            )
        ],
        reverse=True
    )

    return directories[0] if directories else None


def read_step16_reports(folder):

    results = []

    if folder is None:
        return results

    for path in sorted(folder.rglob("*.json")):

        try:
            if path.stat().st_size > 2_000_000:
                continue

            data = json.loads(
                path.read_text(
                    encoding="utf-8"
                )
            )

            results.append({
                "file": str(path),
                "data": data
            })

        except (
            OSError,
            ValueError,
            UnicodeError
        ):
            continue

    return results


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.2")
    print("AUDIT INTERNO MOTORE ZERO TARGET")
    print("=" * 70)

    if not ENGINE.is_file():
        raise SystemExit(
            f"[ERRORE] Motore non trovato: {ENGINE}"
        )

    lines = ENGINE.read_text(
        encoding="utf-8"
    ).splitlines()

    tree = ast.parse("\n".join(lines))

    functions = find_functions(tree)

    relevant = relevant_lines(lines)

    analysis = analyze_function(
        lines,
        functions
    )

    step16 = last_step16()

    previous_reports = read_step16_reports(
        step16
    )

    print("\n[1] MOTORE")

    print("File:", ENGINE)
    print("Righe:", len(lines))

    print("\n[2] FUNZIONI DISPONIBILI")

    for item in functions:

        print(
            f"{item['name']}() "
            f"| righe "
            f"{item['line_start']}-"
            f"{item['line_end']}"
        )

    print("\n[3] FUNZIONE ANALYZE")

    if analysis:

        print(
            "Parametri:",
            analysis["parameters"]
        )

        print(
            "Righe:",
            analysis["start"],
            "-",
            analysis["end"]
        )

    else:
        print("[ERRORE] analyze() non individuata")

    print("\n[4] RICERCA FILTRI E ALGORITMI")

    print(
        "Righe rilevanti:",
        len(relevant)
    )

    for item in relevant[:150]:

        print(
            f"{item['line']:5d} | "
            f"{item['code'][:150]}"
        )

    if len(relevant) > 150:

        print(
            "[INFO] Output terminale limitato "
            "alle prime 150 righe rilevanti."
        )

    print("\n[5] RISULTATI STEP 16")

    print("Cartella:", step16)

    print(
        "Report JSON trovati:",
        len(previous_reports)
    )

    for item in previous_reports:

        print(
            "[JSON]",
            Path(item["file"]).name
        )

    print("\n[6] SALVATAGGIO AUDIT")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_2_AUDIT_MOTORE_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.2",
        "motore": str(ENGINE),
        "funzioni": functions,
        "analyze": analysis,
        "righe_rilevanti": relevant,
        "step16": str(step16),
        "report_step16": previous_reports,
        "natura": "AUDIT_STATICO",
        "motore_eseguito": False,
        "validazione_biologica": False,
    }

    json_file = output / "AUDIT_MOTORE_V33.json"

    json_file.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    text_file = output / "CODICE_MOTORE_RILEVANTE.txt"

    with text_file.open(
        "w",
        encoding="utf-8"
    ) as f:

        f.write(
            "AI PRIMETTA V3.3 - STEP 17.2\n\n"
        )

        f.write(
            "FUNZIONI DEL MOTORE\n\n"
        )

        for item in functions:

            f.write(
                f"{item['name']} "
                f"[{item['line_start']}-"
                f"{item['line_end']}]\n"
            )

        f.write(
            "\nCODICE RILEVANTE PER I TARGET\n\n"
        )

        for item in relevant:

            f.write(
                f"{item['line']:5d} | "
                f"{item['code']}\n"
            )

        if analysis:

            f.write(
                "\nFUNZIONE ANALYZE COMPLETA\n\n"
            )

            for line in analysis["source"]:
                f.write(line + "\n")

    print("\n[OK] JSON:", json_file)

    print(
        "[OK] Codice diagnostico:",
        text_file
    )

    print("\n" + "=" * 70)
    print("[OK] STEP 17.2 TERMINATO")
    print("=" * 70)

    print("[OK] GUI invariata")
    print("[OK] Modelli invariati")
    print("[OK] Nessun training")
    print("[OK] Nessun dato modificato")
    print("[OK] Nessuna modifica a Docker")

    print("\nRisultati:")
    print(output)


if __name__ == "__main__":
    main()