#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA - STEP 22.6J

Collaudo Chat RAG della GUI originale V3.2.

Modalità:
    --audit : controlla script, dati e collegamenti
    --test  : esegue audit + API RAG + chat ricercatore

Non esegue training.
Non rigenera knowledge base.
Non modifica GUI o modelli scientifici.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import urllib.error
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT / "target_prediction_v33_vscode" / "runs"
)

CONTAINER = "maize_ncrna_gnn"

EXPECTED_IMAGE = (
    "ai_primetta_v32_full_ai_sklearn18:22.6g"
)

QUESTION_FDR = (
    "Quali ncRNA caricati dal ricercatore "
    "sono presenti nel ranking e hanno "
    "FDR inferiore a 0.05?"
)

QUESTION_CONTROL = (
    "Quali sono i primi cinque candidati "
    "del ranking esplorativo?"
)

BACKEND_FILES = [
    "16_rag_chat_engine.py",
    "25_train_graph_autoencoder.py",
    "26_validate_user_dataset.py",
    "27_link_user_data_to_ranking.py",
    "28_build_user_rag_knowledge.py",
    "29_user_rag_chat.py",
]

FDR_NAMES = {
    "fdr",
    "padj",
    "p_adj",
    "qvalue",
    "q_value",
    "qval",
    "adjusted_pvalue",
    "adjusted_p_value",
    "adjusted_pvalue_bh",
    "p_val_adj",
    "adj_p_val",
}


def execute(command, timeout=100):

    try:

        p = subprocess.run(
            command,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )

        return {
            "ok": p.returncode == 0,
            "exit_code": p.returncode,
            "stdout": p.stdout,
            "stderr": p.stderr[-2500:],
        }

    except subprocess.TimeoutExpired:

        return {
            "ok": False,
            "error": "TIMEOUT",
        }

    except Exception as exc:

        return {
            "ok": False,
            "error": repr(exc),
        }


def docker_image():

    result = execute([
        "docker",
        "inspect",
        CONTAINER,
        "--format",
        "{{.Config.Image}}",
    ])

    return (
        result["stdout"].strip()
        if result["ok"]
        else None
    )


def inspect_backend_files():

    result = {}

    for name in BACKEND_FILES:

        host = ROOT / name

        test = execute([
            "docker",
            "exec",
            CONTAINER,
            "test",
            "-f",
            "/app/" + name,
        ])

        result[name] = {
            "host_exists": host.is_file(),
            "container_exists": test["ok"],
        }

        print(
            "[OK]" if test["ok"] else "[ERRORE]",
            name,
        )

    return result


def inspect_gae_trainer():

    path = ROOT / "25_train_graph_autoencoder.py"

    if not path.is_file():

        return {
            "exists": False
        }

    terms = re.compile(
        r"torch\.save|joblib\.dump|"
        r"to_csv\(|savefig\(|"
        r"output_dir|model_dir|"
        r"checkpoint|graph_path|"
        r"model_path|--output|--graph",
        re.IGNORECASE,
    )

    matches = []

    for number, line in enumerate(
        path.read_text(
            encoding="utf-8"
        ).splitlines(),
        start=1,
    ):

        if terms.search(line):

            matches.append({
                "line": number,
                "code": line.strip()[:180],
            })

    return {
        "exists": True,
        "relevant_lines": matches[:35],
        "truncated": len(matches) > 35,
    }


def normalize_field(value):

    return re.sub(
        r"[^a-z0-9]+",
        "_",
        value.lower().strip(),
    ).strip("_")


def scan_fdr_tables():

    roots = [
        ROOT / "outputs",
        ROOT / "data/processed",
        ROOT / "data/uploads",
        ROOT / "data/user",
        ROOT / "uploads",
    ]

    excluded_dirs = {
        "raw",
        "fastq",
        "reads",
        "trimmed",
        "downloads",
        "runs",
        ".git",
        ".venv",
        "__pycache__",
    }

    inspected = 0
    limit = 4000
    findings = []

    for base in roots:

        if not base.is_dir():
            continue

        for directory, dirs, files in os.walk(base):

            dirs[:] = sorted(
                d
                for d in dirs
                if d.lower() not in excluded_dirs
            )

            for filename in sorted(files):

                if inspected >= limit:
                    break

                path = Path(directory) / filename

                if path.suffix.lower() not in (
                    ".csv",
                    ".tsv",
                ):
                    continue

                inspected += 1

                try:

                    with path.open(
                        "r",
                        encoding="utf-8-sig",
                        errors="replace",
                    ) as handle:

                        first = handle.readline(32768)

                    separator = (
                        "\t"
                        if first.count("\t")
                        >= first.count(",")
                        else ","
                    )

                    headers = next(
                        csv.reader(
                            [first],
                            delimiter=separator,
                        )
                    )

                    detected = [
                        h
                        for h in headers
                        if normalize_field(h) in FDR_NAMES
                        or "fdr" in normalize_field(h)
                    ]

                    if detected:

                        findings.append({
                            "file": str(
                                path.relative_to(ROOT)
                            ),
                            "fdr_columns": detected,
                            "columns": headers[:35],
                        })

                except (
                    OSError,
                    UnicodeError,
                    csv.Error,
                    StopIteration,
                ):
                    continue

            if inspected >= limit:
                break

        if inspected >= limit:
            break

    return {
        "tables_inspected": inspected,
        "scan_limit_reached": inspected >= limit,
        "fdr_tables": findings[:30],
        "fdr_tables_total_found": len(findings),
    }


def call_rag(question):

    payload = {
        "question": question,
        "top_k": 5,
    }

    request = urllib.request.Request(
        "http://127.0.0.1:8000/rag/ask",
        data=json.dumps(
            payload
        ).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        method="POST",
    )

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    try:

        with opener.open(
            request,
            timeout=90,
        ) as response:

            data = json.loads(
                response.read(1000000)
            )

            return {
                "ok": response.status == 200,
                "http_status": response.status,
                "provider": data.get("provider"),
                "answer": data.get("answer", ""),
                "sources": data.get("sources", []),
            }

    except urllib.error.HTTPError as exc:

        return {
            "ok": False,
            "http_status": exc.code,
            "error": exc.read(1500).decode(
                "utf-8",
                errors="replace",
            ),
        }

    except Exception as exc:

        return {
            "ok": False,
            "error": repr(exc),
        }


def call_researcher_chat():

    return execute(
        [
            "docker",
            "exec",
            "-w",
            "/app",
            CONTAINER,
            "python3",
            "-B",
            "/app/29_user_rag_chat.py",
            "--question",
            QUESTION_FDR,
            "--top-k",
            "5",
        ],
        timeout=120,
    )


def print_chat_result(label, result):

    print()
    print("=" * 65)
    print(label)
    print("=" * 65)

    print(
        "ESITO:",
        "OK" if result.get("ok") else "ERRORE",
    )

    answer = str(
        result.get("answer", "")
    )

    if answer:

        print(
            "Provider:",
            result.get("provider"),
        )

        print(
            "Risposta:",
            answer[:2200],
        )

        print(
            "Fonti:",
            len(result.get("sources", [])),
        )

        print(
            "Contiene riferimento a FDR/padj:",
            any(
                word in answer.lower()
                for word in (
                    "fdr",
                    "padj",
                    "q-value",
                    "qvalue",
                )
            ),
        )

    elif result.get("error"):

        print(
            "Errore:",
            result["error"],
        )


def save_report(report):

    RUNS.mkdir(
        parents=True,
        exist_ok=True,
    )

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    folder = RUNS / (
        "STEP22_6J_COLLAUDO_CHAT_"
        + stamp
    )

    folder.mkdir(
        exist_ok=False
    )

    path = folder / "RAPPORTO_STEP22_6J.json"

    path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    print()
    print("[REPORT]", path)


def main():

    parser = argparse.ArgumentParser()

    modes = parser.add_mutually_exclusive_group(
        required=True
    )

    modes.add_argument(
        "--audit",
        action="store_true",
    )

    modes.add_argument(
        "--test",
        action="store_true",
    )

    args = parser.parse_args()

    print("=" * 65)
    print("AI PRIMETTA V3.2 - STEP 22.6J")
    print("COLLAUDO CHAT RAG E DATI FDR")
    print("=" * 65)

    image = docker_image()

    print("\n[1] IMMAGINE DOCKER:", image)

    if image != EXPECTED_IMAGE:

        raise RuntimeError(
            "Immagine V3.2 diversa "
            "da quella collaudata"
        )

    print("\n[2] SCRIPT BACKEND")

    scripts = inspect_backend_files()

    print("\n[3] DATI TABELLARI CON FDR")

    tables = scan_fdr_tables()

    print(
        "Tabelle esaminate:",
        tables["tables_inspected"],
    )

    print(
        "Tabelle con colonne FDR:",
        tables["fdr_tables_total_found"],
    )

    for item in tables["fdr_tables"][:15]:

        print(
            "[FDR]",
            item["file"],
            item["fdr_columns"],
        )

    print("\n[4] AUDIT SCRIPT TRAINING GAE")

    trainer = inspect_gae_trainer()

    print(
        "Riferimenti a input/output/checkpoint:",
        len(trainer.get("relevant_lines", [])),
    )

    for line in trainer.get(
        "relevant_lines", []
    )[:18]:

        print(
            f'Riga {line["line"]}:',
            line["code"],
        )

    api_fdr = {
        "not_executed": True
    }

    api_control = {
        "not_executed": True
    }

    gui_chat = {
        "not_executed": True
    }

    if args.test:

        print("\n[5] TEST API - DOMANDA FDR")

        api_fdr = call_rag(
            QUESTION_FDR
        )

        print_chat_result(
            "RISPOSTA API FDR",
            api_fdr,
        )

        print("\n[6] TEST API - DOMANDA CONTROLLO")

        api_control = call_rag(
            QUESTION_CONTROL
        )

        print_chat_result(
            "RISPOSTA API CONTROLLO",
            api_control,
        )

        if (
            api_fdr.get("ok")
            and api_control.get("ok")
        ):

            same = (
                api_fdr.get("answer", "").strip()
                ==
                api_control.get("answer", "").strip()
            )

            print()
            print(
                "[CONFRONTO] Risposte identiche:",
                same,
            )

            if same:

                print(
                    "[ATTENZIONE] L'API produce la "
                    "stessa risposta a due domande "
                    "scientificamente differenti."
                )

        print("\n[7] TEST CHAT DATI RICERCATORE")

        if scripts[
            "29_user_rag_chat.py"
        ]["container_exists"]:

            gui_chat = call_researcher_chat()

            print(
                "Codice uscita:",
                gui_chat.get("exit_code"),
            )

            print(
                "Output:",
                gui_chat.get(
                    "stdout", ""
                )[-3500:],
            )

            if gui_chat.get("stderr"):

                print(
                    "Errori:",
                    gui_chat["stderr"][-1200:],
                )

        else:

            gui_chat = {
                "ok": False,
                "error": (
                    "29_user_rag_chat.py assente"
                ),
            }

    report = {
        "step": "22.6J",
        "image": image,
        "backend_scripts": scripts,
        "fdr_tables": tables,
        "gae_trainer_audit": trainer,
        "api_fdr": api_fdr,
        "api_control": api_control,
        "researcher_chat": gui_chat,
        "training_performed": False,
        "knowledge_base_rebuilt": False,
        "gui_modified": False,
        "scientific_fdr_answer_validated": False,
    }

    save_report(report)

    print()
    print("[OK] Nessun training GAE/ML")
    print("[OK] GUI e checkpoint non modificati")
    print(
        "[NOTA] La correttezza scientifica "
        "richiede confronto con i dati originali."
    )


if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__,
            str(exc),
            file=sys.stderr,
        )

        sys.exit(1)