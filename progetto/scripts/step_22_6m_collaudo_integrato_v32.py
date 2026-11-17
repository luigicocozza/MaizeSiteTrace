#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.2 | STEP 22.6M

Collaudo integrato non distruttivo di:
- API e Chat RAG
- persistenza upload
- collegamenti delle cinque pagine
- sicurezza del training GAE legacy

--audit  Analisi statica senza test di scrittura.
--test   Include interrogazioni API e test temporaneo
         di scrittura-lettura del volume upload.

Non avvia training.
Non modifica pagine, API originale o checkpoint.
Non riavvia Docker.
Non rigenera la knowledge base.
"""

import argparse
import ast
import csv
import hashlib
import json
import subprocess
import sys
import urllib.error
import urllib.request
import uuid

from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

CONTAINER = "maize_ncrna_gnn"

IMAGE = "ai_primetta_v32_full_ai_sklearn18:22.6g"

PERSISTENT = ROOT / "data/user_data"

RUNS = ROOT / "target_prediction_v33_vscode/runs"

DATA_DIR = ROOT / "outputs/user_data"

LEGACY_TRAINER = ROOT / "25_train_graph_autoencoder.py"

LEGACY_MODEL = (
    ROOT / "models/graph_autoencoder/gae_model.pt"
)

V3_MODEL = (
    ROOT / "models/graph_autoencoder_v3/gae_v3_model.pt"
)

PAGES = {
    "Dashboard":
        ROOT / "app/streamlit_app.py",

    "GAE":
        ROOT / "app/pages/01_Rete_neurale_GAE.py",

    "Caricamento":
        ROOT / "app/pages/02_Dati_ricercatore.py",

    "Ranking":
        ROOT / "app/pages/03_Ranking_dati_ricercatore.py",

    "Chat ricercatore":
        ROOT / "app/pages/04_Chat_dati_ricercatore.py",
}

SCRIPTS = [
    "16_rag_chat_engine.py",
    "25_train_graph_autoencoder.py",
    "26_validate_user_dataset.py",
    "27_link_user_data_to_ranking.py",
    "28_build_user_rag_knowledge.py",
    "29_user_rag_chat.py",
]

FDR_QUESTION = (
    "Quali ncRNA caricati dal ricercatore "
    "sono presenti nel ranking e hanno "
    "FDR inferiore a 0.05?"
)

GENERAL_QUESTION = (
    "Quali sono i primi cinque candidati "
    "del ranking esplorativo?"
)


# ============================================================
# 2. ESECUZIONE COMANDI
# ============================================================

def run(args, timeout=30, required=True):

    p = subprocess.run(
        [str(a) for a in args],
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )

    if required and p.returncode:

        raise RuntimeError(
            "Comando fallito: "
            + " ".join(str(x) for x in args[:5])
            + "\n"
            + p.stderr[-800:]
        )

    return p


# ============================================================
# 3. INTEGRITÀ FILE
# ============================================================

def checksum(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):

            h.update(block)

    return h.hexdigest()


def protected_hashes():

    files = (
        list(PAGES.values())
        + [
            ROOT / "api/main.py",
            V3_MODEL,
        ]
    )

    if any(
        not p.is_file()
        for p in files
    ):

        raise FileNotFoundError(
            "File GUI, API o GAE V3 assente"
        )

    return {
        str(p.relative_to(ROOT)): checksum(p)
        for p in files
    }


# ============================================================
# 4. STATO DOCKER E PERSISTENZA
# ============================================================

def docker_state():

    raw = run([
        "docker",
        "inspect",
        CONTAINER,
    ])

    state = json.loads(raw.stdout)[0]

    mounts = state.get("Mounts", [])

    volume = next(
        (
            m for m in mounts
            if m.get("Destination")
            == "/app/data/user_data"
        ),
        None,
    )

    api_mount = next(
        (
            m for m in mounts
            if m.get("Destination")
            == "/app/api/main.py"
        ),
        None,
    )

    source_ok = bool(
        volume
        and Path(
            volume["Source"]
        ).resolve() == PERSISTENT.resolve()
    )

    return {

        "running": bool(
            state["State"]["Running"]
        ),

        "image": state["Config"]["Image"],

        "upload_mount": {

            "exists": volume is not None,

            "read_write": bool(
                volume and volume.get("RW")
            ),

            "source_matches_project": source_ok,
        },

        "patched_api_mount": bool(
            api_mount
            and not api_mount.get("RW")
        ),
    }


# ============================================================
# 5. ANALISI DELLE CINQUE PAGINE
# ============================================================

def static_gui():

    report = {}

    for name, path in PAGES.items():

        source = path.read_text(
            encoding="utf-8"
        )

        tree = ast.parse(
            source,
            filename=str(path),
        )

        buttons = []

        for n in ast.walk(tree):

            if not isinstance(n, ast.Call):
                continue

            fun = ast.unparse(n.func)

            if fun not in (
                "st.button",
                "st.form_submit_button",
            ):
                continue

            label = (
                n.args[0].value

                if (
                    n.args
                    and isinstance(
                        n.args[0],
                        ast.Constant,
                    )
                    and isinstance(
                        n.args[0].value,
                        str,
                    )
                )

                else "<dinamico>"
            )

            buttons.append({

                "label": label,

                "line": n.lineno,
            })

        markers = [

            {
                "line": i,

                "text": s.strip()[:130],
            }

            for i, s in enumerate(
                source.splitlines(),
                1,
            )

            if any(
                k in s
                for k in (

                    "SCRIPT_PATH =",

                    "UPLOAD_ROOT =",

                    "VALIDATOR =",

                    "LINK_SCRIPT =",

                    "BUILD_KB_SCRIPT =",

                    "CHAT_SCRIPT =",

                    "subprocess.run(",
                )
            )
        ]

        report[name] = {

            "syntax_ok": True,

            "buttons": buttons,

            "connections": markers[:14],
        }

    return report


# ============================================================
# 6. VERIFICA SCRIPT NEL CONTAINER
# ============================================================

def scripts_visible():

    probe = (

        "import json,pathlib;ns="
        + repr(SCRIPTS)
        + ";"

        "print(json.dumps({"
        "n:(pathlib.Path('/app')/n).is_file() "
        "for n in ns"
        "}))"
    )

    p = run([
        "docker",
        "exec",
        CONTAINER,
        "python3",
        "-B",
        "-c",
        probe,
    ])

    return json.loads(
        p.stdout.strip()
    )


# ============================================================
# 7. SICUREZZA TRAINING GAE
# ============================================================

def trainer_safety():

    s = LEGACY_TRAINER.read_text(
        encoding="utf-8"
    )

    ast.parse(
        s,
        filename=str(LEGACY_TRAINER),
    )

    lines = s.splitlines()

    refs = [

        {
            "line": i,

            "code": v.strip()[:150],
        }

        for i, v in enumerate(
            lines,
            1,
        )

        if any(
            key in v
            for key in (

                "MODEL_DIR =",

                "MODEL_PATH =",

                "torch.save(",

                "NODE_EMB_CSV =",

                "NCRNA_EMB_CSV =",

                "HISTORY_CSV =",

                "PLOTS_DIR =",
            )
        )
    ]

    return {

        "legacy_trainer_exists": True,

        "legacy_model_exists":
            LEGACY_MODEL.is_file(),

        "v3_model_exists":
            V3_MODEL.is_file(),

        "legacy_and_v3_paths_distinct": (
            LEGACY_MODEL.resolve()
            != V3_MODEL.resolve()
        ),

        "training_button_executed": False,

        "warning": (
            "Non avviare il training GAE "
            "finché gli output legacy "
            "non sono isolati."
        ),

        "references": refs,
    }


# ============================================================
# 8. INDIVIDUAZIONE RANKING INTEGRATO
# ============================================================

def selected_csv():

    paths = list(
        DATA_DIR.glob(
            "*/ranking_with_user_evidence.csv"
        )
    )

    if not paths:

        raise FileNotFoundError(
            "Nessun ranking integrato "
            "in outputs/user_data"
        )

    return max(
        paths,
        key=lambda p: (
            p.stat().st_mtime_ns,
            p.name,
        ),
    )


# ============================================================
# 9. CALCOLO FDR DIRETTO SUL CSV
# ============================================================

def matching_rows(path):

    selected = []

    with path.open(
        encoding="utf-8-sig",
        newline="",
    ) as f:

        reader = csv.DictReader(f)

        required = {
            "rnacentral_id",
            "user_min_FDR",
        }

        if not required.issubset(
            reader.fieldnames or []
        ):

            raise ValueError(
                "Colonne indispensabili assenti"
            )

        for row in reader:

            try:

                fdr = Decimal(
                    (
                        row.get("user_min_FDR")
                        or ""
                    ).strip().replace(",", ".")
                )

            except InvalidOperation:
                continue

            if not (
                fdr.is_finite()
                and Decimal("0") <= fdr
                < Decimal("0.05")
            ):
                continue

            flag = (
                row.get("user_data_present")
                or ""
            ).strip().lower()

            if (
                flag
                and flag not in (
                    "true",
                    "1",
                    "yes",
                    "si",
                )
            ):
                continue

            num = (
                row.get("user_records_count")
                or ""
            ).strip()

            if num and Decimal(num) <= 0:
                continue

            selected.append(
                (
                    row["rnacentral_id"]
                    .strip()
                    .upper(),

                    fdr,
                )
            )

    return sorted(selected)


# ============================================================
# 10. INTERROGAZIONE API RAG
# ============================================================

def ask(question, top_k=50):

    payload = json.dumps({

        "question": question,

        "top_k": top_k,

    }).encode("utf-8")

    req = urllib.request.Request(

        "http://127.0.0.1:8000/rag/ask",

        data=payload,

        headers={
            "Content-Type": "application/json"
        },

        method="POST",
    )

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    try:

        with opener.open(
            req,
            timeout=45,
        ) as result:

            return json.loads(
                result.read(300000)
            )

    except urllib.error.HTTPError as exc:

        raise RuntimeError(
            f"API RAG HTTP {exc.code}: "
            f"{exc.read(500)!r}"
        ) from exc


# ============================================================
# 11. TEST EFFETTIVO VOLUME UPLOAD
# ============================================================

def upload_probe():

    """
    Crea un file temporaneo nel container.
    Verifica la comparsa sull'host.
    Elimina il file di prova.
    """

    filename = (
        ".primetta_step22_6m_"
        + uuid.uuid4().hex
        + ".txt"
    )

    payload = (
        "STEP22_6M_VOLUME_TEST_"
        + uuid.uuid4().hex
    )

    container_code = (

        "from pathlib import Path;"
        "import sys;"

        "d=Path('/app/data/user_data');"

        "p=d/sys.argv[1];"

        "assert not p.exists();"

        "p.write_text("
        "sys.argv[2],encoding='utf-8');"

        "print(str(p))"
    )

    host = PERSISTENT / filename

    try:

        run([
            "docker",
            "exec",
            CONTAINER,
            "python3",
            "-B",
            "-c",
            container_code,
            filename,
            payload,
        ])

        if (
            not host.is_file()
            or host.read_text(
                encoding="utf-8"
            ) != payload
        ):

            raise RuntimeError(
                "Scrittura Docker non visibile "
                "correttamente sull'host"
            )

        return {

            "write_in_container": True,

            "read_on_host": True,

            "cleanup": (
                "eseguita nel blocco finally"
            ),
        }

    finally:

        if host.exists():

            if (
                host.is_symlink()
                or not host.is_file()
            ):

                raise RuntimeError(
                    "File di test inatteso"
                )

            host.unlink()


# ============================================================
# 12. SALVATAGGIO RAPPORTO
# ============================================================

def report_file(payload):

    date = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP22_6M_COLLAUDO_V32_"
        + date
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    outfile = (
        folder / "RAPPORTO_STEP22_6M.json"
    )

    outfile.write_text(

        json.dumps(
            payload,
            indent=2,
            ensure_ascii=False,
            default=str,
        ),

        encoding="utf-8",
    )

    print(
        "\n[REPORT]",
        outfile
    )


# ============================================================
# 13. PROGRAMMA PRINCIPALE
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
        action="store_true",
    )

    mode.add_argument(
        "--test",
        action="store_true",
    )

    args = parser.parse_args()

    print("=" * 70)

    print(
        "AI PRIMETTA V3.2 - STEP 22.6M"
    )

    print(
        "COLLAUDO INTEGRATO SENZA TRAINING"
    )

    print("=" * 70)

    before = protected_hashes()

    docker = docker_state()

    if (
        not docker["running"]
        or docker["image"] != IMAGE
    ):

        raise RuntimeError(
            "Container V3.2 non attivo "
            "con immagine 22.6g"
        )

    # --------------------------------------------------------
    # A. ANALISI GUI
    # --------------------------------------------------------

    pages = static_gui()

    scripts = scripts_visible()

    trainer = trainer_safety()

    # --------------------------------------------------------
    # B. DATI FDR DOCUMENTATI
    # --------------------------------------------------------

    source = selected_csv()

    expected = matching_rows(source)

    print(
        "[DOCKER]",
        docker["image"]
    )

    print(
        "[VOLUME]",
        docker["upload_mount"]
    )

    print(
        "[API MONTATA]",
        docker["patched_api_mount"]
    )

    for title, info in pages.items():

        print(
            "[GUI]",
            title,
            "| pulsanti:",
            [
                x["label"]
                for x in info["buttons"]
            ]
        )

        for connection in info["connections"]:

            print(
                "   ",
                connection["line"],
                connection["text"]
            )

    for name, exists in scripts.items():

        print(
            "[SCRIPT]",
            "OK" if exists else "ASSENTE",
            name
        )

    # --------------------------------------------------------
    # C. VERIFICA SEPARAZIONE GAE
    # --------------------------------------------------------

    for line in trainer["references"]:

        print(
            "[GAE LEGACY]",
            line["line"],
            line["code"]
        )

    print(
        "[GAE] Modello precedente esiste:",
        trainer["legacy_model_exists"]
    )

    print(
        "[GAE] Modello V3 distinto:",
        trainer[
            "legacy_and_v3_paths_distinct"
        ]
    )

    # --------------------------------------------------------
    # D. CSV ORIGINALE
    # --------------------------------------------------------

    print(
        "[FDR] File sorgente:",
        source.relative_to(ROOT)
    )

    print(
        "[FDR] Record documentati "
        "con FDR < 0.05:",
        len(expected)
    )

    for identifier, fdr in expected:

        print(
            "   ",
            identifier,
            "FDR=",
            str(fdr)
        )

    tests = {}

    # --------------------------------------------------------
    # E. COLLAUDI FUNZIONALI
    # --------------------------------------------------------

    if args.test:

        if not (
            docker["upload_mount"][
                "read_write"
            ]

            and docker["upload_mount"][
                "source_matches_project"
            ]
        ):

            raise RuntimeError(
                "Volume upload non persistente "
                "o origine non conforme"
            )

        print(
            "\n[TEST] Scrittura e lettura "
            "del volume persistente"
        )

        tests[
            "upload_persistence"
        ] = upload_probe()

        print(
            "[OK] Scrittura/lettura riuscita; "
            "file temporaneo eliminato"
        )

        print(
            "\n[TEST] API FDR e domanda generica"
        )

        fdr_answer = ask(
            FDR_QUESTION
        )

        generic_answer = ask(
            GENERAL_QUESTION
        )

        observed = sorted(
            (
                str(
                    r["rnacentral_id"]
                ).upper(),

                Decimal(
                    str(r["user_min_FDR"])
                ),
            )

            for r in fdr_answer.get(
                "evidence",
                []
            )
        )

        fdr_pass = (

            fdr_answer.get("provider")
            == "deterministic_user_fdr_v32"

            and fdr_answer.get(
                "total_matches"
            ) == len(expected)

            and observed == expected
        )

        generic_pass = (
            generic_answer.get("provider")
            != "deterministic_user_fdr_v32"
        )

        tests[
            "fdr_api_matching_csv"
        ] = fdr_pass

        tests[
            "generic_rag_separate"
        ] = generic_pass

        tests[
            "fdr_api_provider"
        ] = fdr_answer.get(
            "provider"
        )

        tests[
            "generic_api_provider"
        ] = generic_answer.get(
            "provider"
        )

        print(
            "[FDR API]",
            "PASS" if fdr_pass else "FAIL",
            "|",
            len(observed),
            "RNA |",
            fdr_answer.get("provider")
        )

        print(
            "[RAG GENERICA]",
            "PASS" if generic_pass else "FAIL",
            "|",
            generic_answer.get("provider")
        )

    else:

        print(
            "\n[AUDIT] Nessuna richiesta RAG "
            "o scrittura di file nel volume"
        )

    # --------------------------------------------------------
    # F. INTEGRITÀ FINALE
    # --------------------------------------------------------

    if protected_hashes() != before:

        raise RuntimeError(
            "File GUI, API originale "
            "o modello GAE V3 modificati!"
        )

    output = {

        "step": "22.6M",

        "mode": (
            "test"
            if args.test
            else "audit"
        ),

        "docker": docker,

        "pages": pages,

        "scripts_visible": scripts,

        "gae_legacy_safety": trainer,

        "fdr_source": str(
            source.relative_to(ROOT)
        ),

        "fdr_expected": [

            {
                "id": identifier,
                "fdr": str(fdr)
            }

            for identifier, fdr in expected
        ],

        "tests": tests,

        "original_hashes": before,

        "gae_training_performed": False,

        "gui_buttons_executed": False,

        "knowledge_base_rebuilt": False,

        "gui_modified": False,
    }

    report_file(output)

    print(
        "[OK] Nessun training; "
        "GUI e checkpoint invariati"
    )

    if (
        args.test
        and (
            not tests[
                "fdr_api_matching_csv"
            ]

            or not tests[
                "generic_rag_separate"
            ]
        )
    ):

        sys.exit(2)


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