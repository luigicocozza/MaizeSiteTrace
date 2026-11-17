#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.2 - STEP 22.6K

Individua il codice FastAPI che gestisce /rag/ask,
il collegamento al motore RAG generale e il
backend della Chat dati ricercatore.

Controlla anche la persistenza dei file caricati.

Non modifica GUI, API, Docker, modelli o dati.
Non esegue training o chiamate RAG.
"""

import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

CONTAINER = "maize_ncrna_gnn"

IMAGE = "ai_primetta_v32_full_ai_sklearn18:22.6g"

MERGED = (
    ROOT
    / "outputs/user_data/"
    "esperimento_ricercatore_01_20261009_100916/"
    "ranking_with_user_evidence.csv"
)

PATTERN = re.compile(
    r"/rag/ask|16_rag_chat_engine|"
    r"29_user_rag_chat|"
    r"structured_rag|"
    r"deterministic_user_data_filter",
    re.IGNORECASE,
)

REDACT = re.compile(
    r"password|secret|bearer|"
    r"api_key|access_token",
    re.IGNORECASE,
)


def execute(args, timeout=50):

    p = subprocess.run(
        args,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
    )

    return {
        "ok": p.returncode == 0,
        "stdout": p.stdout,
        "stderr": p.stderr,
    }


def safe_line(line):

    if REDACT.search(line):
        return "[CONTENUTO RISERVATO OSCURATO]"

    return line.rstrip()[:220]


def scan_source(path, before=12, after=28):

    try:
        text = path.read_text(
            encoding="utf-8"
        )

    except (OSError, UnicodeError):
        return []

    lines = text.splitlines()

    hit_indices = [
        i
        for i, line in enumerate(lines)
        if PATTERN.search(line)
    ]

    ranges = []

    for i in hit_indices:

        start = max(
            0,
            i - before
        )

        stop = min(
            len(lines),
            i + after + 1
        )

        ranges.append((start, stop))

    merged = []

    for start, stop in sorted(ranges):

        if merged and start <= merged[-1][1]:

            old_start, old_stop = merged[-1]

            merged[-1] = (
                old_start,
                max(old_stop, stop)
            )

        else:

            merged.append((start, stop))

    result = []

    for start, stop in merged[:5]:

        excerpt = [
            f"{i + 1:04d}: {safe_line(lines[i])}"
            for i in range(start, stop)
        ]

        result.append(
            "\n".join(excerpt)
        )

    return result


def scan_project():

    print("\n[1] ROUTER RAG NEL REPOSITORY")

    roots = [
        ROOT,
        ROOT / "api",
        ROOT / "app",
        ROOT / "src",
    ]

    visited = set()
    found = 0

    for base in roots:

        if not base.exists():
            continue

        for current, dirs, files in os.walk(base):

            dirs[:] = [
                d
                for d in dirs
                if d not in {
                    "outputs",
                    "data",
                    "models",
                    "runs",
                    ".venv",
                    ".git",
                    "__pycache__",
                    "node_modules",
                    "target_prediction_v33_vscode",
                }
            ]

            for name in sorted(files):

                if not name.endswith(".py"):
                    continue

                path = Path(current) / name

                if path in visited:
                    continue

                visited.add(path)

                excerpts = scan_source(path)

                if not excerpts:
                    continue

                found += 1

                print(
                    "\n[FILE]",
                    path.relative_to(ROOT)
                )

                for excerpt in excerpts:
                    print(excerpt)

                if found >= 12:
                    return

    if not found:
        print(
            "[NOTA] Nessun riferimento trovato "
            "nei file Python del repository"
        )


def scan_container():

    print("\n[2] ROUTER RAG NEL CONTAINER")

    script = r'''
import os
import re
from pathlib import Path

pattern = re.compile(
    r"/rag/ask|16_rag_chat_engine|structured_rag",
    re.I
)

private = re.compile(
    r"password|secret|bearer|api_key",
    re.I
)

count = 0

for base, dirs, files in os.walk("/app"):

    dirs[:] = [
        d for d in dirs
        if d not in {
            "outputs", "data", "models", "logs",
            "__pycache__", ".venv", "node_modules"
        }
    ]

    for name in sorted(files):

        if not name.endswith(".py"):
            continue

        path = Path(base) / name

        try:
            lines = path.read_text(
                encoding="utf-8"
            ).splitlines()
        except Exception:
            continue

        hits = [
            i for i, line in enumerate(lines)
            if pattern.search(line)
        ]

        if not hits:
            continue

        print("\n[CONTAINER FILE]", path)

        for i in hits[:5]:

            start = max(0, i - 8)
            stop = min(len(lines), i + 24)

            for j in range(start, stop):

                line = lines[j].strip()

                if private.search(line):
                    line = "[OSCURATO]"

                print(
                    f"{j + 1:04d}: {line[:200]}"
                )

        count += 1

        if count >= 10:
            break

    if count >= 10:
        break
'''

    output = execute([
        "docker",
        "exec",
        CONTAINER,
        "python3",
        "-B",
        "-c",
        script,
    ])

    if output["ok"]:
        print(output["stdout"][:17000])

    else:
        print(
            "[ERRORE]",
            output["stderr"][-1000:]
        )


def inspect_startup():

    print("\n[3] AVVIO API FASTAPI")

    result = execute([
        "docker",
        "top",
        CONTAINER,
        "-eo",
        "pid,args",
    ])

    if result["ok"]:

        print(
            result["stdout"][:3500]
        )

    else:

        print(
            "[NOTA]",
            result["stderr"][-500:]
        )

    for relative in (
        "docker/start_services.sh",
        "docker/start.sh",
    ):

        path = ROOT / relative

        if not path.is_file():
            continue

        print("\n[AVVIO]", relative)

        for i, line in enumerate(
            path.read_text(
                encoding="utf-8"
            ).splitlines()[:100],
            start=1,
        ):

            print(
                f"{i:03d}: {safe_line(line)}"
            )


def verify_user_data():

    print("\n[4] VERIFICA DATI REALI FDR")

    if not MERGED.is_file():

        print(
            "[ERRORE] File integrato non trovato:",
            MERGED
        )
        return

    with MERGED.open(
        encoding="utf-8-sig",
        newline="",
    ) as f:

        reader = csv.DictReader(f)

        columns = reader.fieldnames or []

        print(
            "[COLONNE]",
            columns
        )

        fdr_col = next(
            (
                col for col in columns
                if col.lower() == "user_min_fdr"
            ),
            None,
        )

        if fdr_col is None:

            print(
                "[ERRORE] user_min_FDR assente"
            )
            return

        selected = []

        for row in reader:

            raw = row.get(
                fdr_col, ""
            ).strip()

            if not raw:
                continue

            try:
                fdr = float(raw)

            except ValueError:
                continue

            if not (
                0 <= fdr < 0.05
            ):
                continue

            ids = []

            for value in row.values():

                ids.extend(
                    re.findall(
                        r"URS[0-9A-F]{11}",
                        str(value),
                        flags=re.I,
                    )
                )

            selected.append({
                "rna_ids": sorted(set(ids)),
                "fdr": fdr,
            })

        print(
            "[RISULTATI FDR < 0.05]",
            len(selected)
        )

        for row in selected[:15]:

            print(
                "[RNA]",
                row["rna_ids"],
                "| FDR =",
                row["fdr"]
            )


def inspect_mounts():

    print(
        "\n[5] PERSISTENZA DATI RICERCATORE"
    )

    result = execute([
        "docker",
        "inspect",
        CONTAINER,
    ])

    if not result["ok"]:
        print("[ERRORE] Docker inspect")
        return

    obj = json.loads(
        result["stdout"]
    )[0]

    mounts = obj.get(
        "Mounts", []
    )

    for mount in mounts:

        destination = mount.get(
            "Destination", ""
        )

        print(
            "[MOUNT]",
            destination,
            "RW =" + str(
                mount.get("RW")
            )
        )

    for target in (
        "/app/outputs",
        "/app/data/user_data",
        "/app/data/processed",
    ):

        persisted = any(
            target == m.get("Destination")
            or target.startswith(
                m.get("Destination", "") + "/"
            )
            for m in mounts
        )

        print(
            "[PERSISTENZA]",
            target,
            "SI" if persisted else "DA VERIFICARE"
        )


def main():

    print("=" * 68)
    print("AI PRIMETTA V3.2 — STEP 22.6K")
    print("INDIVIDUAZIONE ROUTER RAG E DATI FDR")
    print("=" * 68)

    image = execute([
        "docker",
        "inspect",
        CONTAINER,
        "--format",
        "{{.Config.Image}}",
    ])

    active = image["stdout"].strip()

    print(
        "[IMMAGINE ATTIVA]",
        active
    )

    if active != IMAGE:

        raise RuntimeError(
            "Immagine Docker non corrispondente "
            "alla V3.2 collaudata"
        )

    scan_project()
    scan_container()
    inspect_startup()
    verify_user_data()
    inspect_mounts()

    print("\n" + "=" * 68)
    print("[OK] Diagnostica conclusa")
    print("[OK] Nessun training")
    print("[OK] Nessun file modificato")
    print("=" * 68)


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