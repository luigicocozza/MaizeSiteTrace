#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA STEP 22.6H

Diagnostica dei collegamenti GUI V3.2
e test RAG opzionale.

Esecuzione:

python3 scripts/step_22_6h_audit_gui_v32.py --audit

python3 scripts/step_22_6h_audit_gui_v32.py --rag-test

Non modifica codice, dati, GUI, modelli o container.
Non esegue training.
"""

import argparse
import ast
import json
import re
import subprocess
import sys
import urllib.error
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

PAGES = {
    "Dashboard": ROOT / "app/streamlit_app.py",

    "GAE": ROOT / "app/pages/01_Rete_neurale_GAE.py",

    "Dati ricercatore":
        ROOT / "app/pages/02_Dati_ricercatore.py",

    "Ranking":
        ROOT / "app/pages/03_Ranking_dati_ricercatore.py",

    "Chat RAG":
        ROOT / "app/pages/04_Chat_dati_ricercatore.py",
}

RUNS = ROOT / "target_prediction_v33_vscode/runs"

API_BASE = "http://127.0.0.1:8000"

CONTAINER = "maize_ncrna_gnn"

EXPECTED_IMAGE = (
    "ai_primetta_v32_full_ai_sklearn18:22.6g"
)

QUESTION = (
    "Quali ncRNA caricati dal ricercatore "
    "sono presenti nel ranking e hanno "
    "FDR inferiore a 0.05?"
)

SENSITIVE = re.compile(
    r"(?i)(password|secret|api.?key|"
    r"authorization|bearer|token\s*=)"
)


# ============================================================
# 2. STATO CONTAINER DOCKER
# ============================================================

def docker_state():

    process = subprocess.run(
        [
            "docker",
            "inspect",
            CONTAINER,
            "--format",
            "{{json .}}"
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=15
    )

    if process.returncode:

        return {
            "ok": False,
            "error": process.stderr.strip()[-300:]
        }

    data = json.loads(process.stdout)

    return {
        "ok": data.get(
            "State", {}
        ).get("Running", False),

        "image": data.get(
            "Config", {}
        ).get("Image"),

        "ports": data.get(
            "NetworkSettings", {}
        ).get("Ports"),

        "mounts": [
            {
                "destination": m.get("Destination"),

                "read_only": not m.get(
                    "RW", True
                )
            }

            for m in data.get("Mounts", [])
        ]
    }


# ============================================================
# 3. ANALISI STATICA DEI PULSANTI GUI
# ============================================================

def extract_ui(name, path):

    if not path.is_file():

        return {
            "exists": False,
            "error": "File assente"
        }

    source = path.read_text(
        encoding="utf-8"
    )

    lines = source.splitlines()

    tree = ast.parse(
        source,
        filename=str(path)
    )

    buttons = []

    references = []

    # --------------------------------------------------------
    # VISITA DELL'ALBERO SINTATTICO
    # --------------------------------------------------------

    def visit(node, parents):

        if isinstance(node, ast.Call):

            try:
                function = ast.unparse(
                    node.func
                )

            except Exception:
                function = ""

            if function in (
                "st.button",
                "st.form_submit_button"
            ):

                label = "<dinamico>"

                if (
                    node.args
                    and isinstance(
                        node.args[0],
                        ast.Constant
                    )
                    and isinstance(
                        node.args[0].value,
                        str
                    )
                ):

                    label = node.args[0].value

                owner = next(
                    (
                        parent

                        for parent in reversed(parents)

                        if isinstance(
                            parent,
                            (
                                ast.If,
                                ast.With,
                                ast.FunctionDef
                            )
                        )
                    ),
                    None
                )

                start = getattr(
                    owner,
                    "lineno",
                    node.lineno
                )

                end = getattr(
                    owner,
                    "end_lineno",
                    node.end_lineno or node.lineno
                )

                limit = (
                    65
                    if name == "GAE"
                    else 22
                )

                excerpt = []

                for i in range(
                    start,
                    min(
                        end,
                        start + limit - 1,
                        len(lines)
                    ) + 1
                ):

                    line = lines[
                        i - 1
                    ].rstrip()

                    if SENSITIVE.search(line):
                        line = "[riga oscurata]"

                    excerpt.append(
                        f"{i:04d}: {line}"
                    )

                buttons.append({

                    "label": label,

                    "line": node.lineno,

                    "handler_lines": [
                        start,
                        end
                    ],

                    "handler_excerpt": excerpt,

                    "excerpt_truncated": (
                        end >= start + limit
                    )
                })

        for child in ast.iter_child_nodes(node):

            visit(
                child,
                parents + [node]
            )

    visit(tree, [])

    # --------------------------------------------------------
    # RICERCA RIFERIMENTI BACKEND
    # --------------------------------------------------------

    keywords = (
        "subprocess",
        "requests.",
        "httpx.",
        "torch.load",
        "train_",
        "fastapi",
        "api_url",
        "rag/ask",
        "pickle.load"
    )

    for i, line in enumerate(
        lines,
        start=1
    ):

        if any(
            word in line.lower()
            for word in keywords
        ):

            references.append({

                "line": i,

                "source": (
                    "[riga oscurata]"
                    if SENSITIVE.search(line)
                    else line.strip()[:180]
                )
            })

    return {
        "exists": True,

        "lines": len(lines),

        "buttons": buttons,

        "backend_references": references[:45],

        "backend_references_truncated": (
            len(references) > 45
        )
    }


# ============================================================
# 4. RICHIESTE HTTP AL BACKEND
# ============================================================

def http_request(
    url,
    body=None,
    timeout=20,
    max_bytes=350000
):

    headers = {
        "Accept": "application/json"
    }

    if body is not None:

        headers[
            "Content-Type"
        ] = "application/json"

    request = urllib.request.Request(
        url,

        headers=headers,

        data=(
            json.dumps(body).encode("utf-8")
            if body is not None
            else None
        ),

        method=(
            "POST"
            if body is not None
            else "GET"
        )
    )

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    try:

        with opener.open(
            request,
            timeout=timeout
        ) as response:

            raw = response.read(
                max_bytes + 1
            )

            if len(raw) > max_bytes:

                return {
                    "ok": False,

                    "status": response.status,

                    "error": "Risposta troppo grande"
                }

            text = raw.decode(
                "utf-8",
                errors="replace"
            )

            try:
                parsed = json.loads(text)

            except json.JSONDecodeError:

                parsed = {
                    "text_excerpt": text[:800]
                }

            return {
                "ok": response.status == 200,

                "status": response.status,

                "response": parsed
            }

    except urllib.error.HTTPError as exc:

        return {
            "ok": False,

            "status": exc.code,

            "error": exc.read(
                1200
            ).decode(
                "utf-8",
                errors="replace"
            )
        }

    except Exception as exc:

        return {
            "ok": False,

            "error": repr(exc)
        }


# ============================================================
# 5. RISOLUZIONE SCHEMA OPENAPI
# ============================================================

def resolve_schema(
    spec,
    schema,
    max_depth=5
):

    current = schema or {}

    for _ in range(max_depth):

        reference = current.get("$ref")

        if not reference:
            break

        prefix = "#/components/schemas/"

        if not reference.startswith(prefix):

            raise ValueError(
                "Riferimento non previsto: "
                + reference
            )

        name = reference.split("/")[-1]

        current = (
            spec.get("components", {})
            .get("schemas", {})
            .get(name, {})
        )

    return current


# ============================================================
# 6. ANALISI SCHEMA CHAT RAG
# ============================================================

def rag_request_schema(spec):

    endpoint = (
        spec.get("paths", {})
        .get("/rag/ask", {})
    )

    method = endpoint.get("post")

    if not method:

        return {
            "declared": bool(endpoint),
            "post": False
        }

    body = (
        method.get("requestBody", {})
        .get("content", {})
        .get("application/json", {})
        .get("schema", {})
    )

    schema = resolve_schema(
        spec,
        body
    )

    return {
        "declared": True,

        "post": True,

        "properties": sorted(
            schema.get(
                "properties", {}
            ).keys()
        ),

        "required": schema.get(
            "required", []
        ),

        "schema": schema
    }


# ============================================================
# 7. TEST EFFETTIVO CHAT RAG
# ============================================================

def execute_rag(
    spec,
    schema_info
):

    if not schema_info.get("post"):

        return {
            "ok": False,

            "skipped": (
                "Endpoint POST /rag/ask "
                "non dichiarato"
            )
        }

    properties = schema_info.get(
        "properties", []
    )

    required = schema_info.get(
        "required", []
    )

    question_key = next(
        (
            field

            for field in (
                "question",
                "query",
                "prompt"
            )

            if field in properties
        ),
        None
    )

    if question_key is None:

        return {
            "ok": False,

            "skipped": (
                "Schema domanda sconosciuto; "
                "nessuna richiesta inviata"
            )
        }

    payload = {
        question_key: QUESTION
    }

    for field in (
        "top_k",
        "max_chunks",
        "k"
    ):

        if field in properties:

            payload[field] = 5

            break

    missing = [
        field
        for field in required

        if field not in payload
    ]

    if missing:

        return {
            "ok": False,

            "skipped": (
                "Campi obbligatori "
                "non riconosciuti: "
                + ", ".join(missing)
            )
        }

    result = http_request(
        API_BASE + "/rag/ask",
        payload,
        timeout=90
    )

    response = result.get(
        "response"
    )

    if isinstance(response, dict):

        response_text = json.dumps(
            response,
            ensure_ascii=False
        )

        result["response"] = {
            "keys": list(
                response.keys()
            ),

            "excerpt": response_text[:1600],

            "truncated": (
                len(response_text) > 1600
            )
        }

    return {
        "payload": payload,
        **result
    }


# ============================================================
# 8. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    mode.add_argument(
        "--audit",
        action="store_true"
    )

    mode.add_argument(
        "--rag-test",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 72)

    print(
        "AI PRIMETTA V3.2 | STEP 22.6H"
    )

    print(
        "AUDIT COLLEGAMENTI GUI"
    )

    print("=" * 72)

    # --------------------------------------------------------
    # 9. STATO DOCKER
    # --------------------------------------------------------

    container = docker_state()

    print(
        "[CONTAINER]",

        (
            "UP"
            if container.get("ok")
            else "NON DISPONIBILE"
        ),

        "| Immagine:",

        container.get(
            "image", "N/D"
        )
    )

    if (
        container.get("image")
        != EXPECTED_IMAGE
    ):

        print(
            "[ATTENZIONE] Immagine differente "
            "da quella dello STEP 22.6G"
        )

    # --------------------------------------------------------
    # 10. ANALISI DELLE PAGINE
    # --------------------------------------------------------

    pages = {}

    for name, path in PAGES.items():

        info = extract_ui(
            name,
            path
        )

        pages[name] = info

        print(
            f"\n[PAGINA] {name}",

            "| presente=",

            info.get("exists"),

            "| pulsanti=",

            len(
                info.get(
                    "buttons", []
                )
            )
        )

        for button in info.get(
            "buttons", []
        ):

            print(
                "Pulsante:",

                button["label"],

                "| riga",

                button["line"],

                "| blocco",

                button["handler_lines"]
            )

            if name in (
                "GAE",
                "Chat RAG",
                "Dati ricercatore",
                "Ranking"
            ):

                print(
                    "\n".join(
                        button["handler_excerpt"]
                    )
                )

                if button[
                    "excerpt_truncated"
                ]:

                    print(
                        "... [blocco abbreviato]"
                    )

        references = [
            f'{item["line"]}:{item["source"]}'

            for item in info.get(
                "backend_references", []
            )[:15]
        ]

        print(
            "Riferimenti backend:",
            references
        )

    # --------------------------------------------------------
    # 11. OPENAPI
    # --------------------------------------------------------

    api = http_request(
        API_BASE + "/openapi.json"
    )

    spec = (
        api.get(
            "response", {}
        )
        if api.get("ok")
        else {}
    )

    schema_info = (
        rag_request_schema(spec)

        if isinstance(spec, dict)

        else {
            "declared": False
        }
    )

    print(
        "\n[API OpenAPI]",

        (
            "OK"
            if api.get("ok")
            else api.get("error")
        )
    )

    paths = (
        sorted(
            spec.get(
                "paths", {}
            )
        )

        if isinstance(spec, dict)

        else []
    )

    print(
        "[ROTTE]",
        paths
    )

    print(
        "[RAG SCHEMA]",

        {
            "post": schema_info.get(
                "post"
            ),

            "properties": schema_info.get(
                "properties"
            ),

            "required": schema_info.get(
                "required"
            )
        }
    )

    # --------------------------------------------------------
    # 12. TEST CHAT RAG OPZIONALE
    # --------------------------------------------------------

    rag = {
        "not_run": True
    }

    if args.rag_test:

        rag = execute_rag(
            spec,
            schema_info
        )

        print(
            "[TEST RAG]",

            {
                key: value

                for key, value in rag.items()

                if key != "response"
            }
        )

        print(
            "[RISPOSTA RAG]",

            rag.get(
                "response",

                rag.get(
                    "error",

                    rag.get("skipped")
                )
            )
        )

    # --------------------------------------------------------
    # 13. RAPPORTO JSON
    # --------------------------------------------------------

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    output = RUNS / (
        "STEP22_6H_AUDIT_GUI_V32_"
        + stamp
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {
        "step": "22.6H",

        "created_utc": stamp,

        "container": container,

        "pages": pages,

        "api_paths": paths,

        "rag_schema": {
            key: value

            for key, value in schema_info.items()

            if key != "schema"
        },

        "rag_test": rag,

        "training_performed": False,

        "existing_files_modified": False,

        "note": (
            "Audit statico dei pulsanti. "
            "Non li esegue. "
            "Test RAG POST solo con --rag-test."
        )
    }

    report_path = (
        output / "RAPPORTO_STEP22_6H.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    print(
        "\n[REPORT]",
        report_path
    )

    print(
        "[OK] Nessun training"
    )

    print(
        "[OK] GUI, modelli e Docker invariati"
    )


# ============================================================
# 14. AVVIO
# ============================================================

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