#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.2 | STEP 22.6L

Ripristina il routing FDR dell'API e rende
persistenti gli upload del ricercatore.

Non modifica:
- GUI originale V3.2
- api/main.py originale
- Checkpoint GAE e ML
- Dataset scientifici
- AI ORZO
- GUI V3.3
- Monitor ML

Comandi:

python3 scripts/step_22_6l_router_rag_persistenza.py --install

python3 scripts/step_22_6l_router_rag_persistenza.py --status

python3 scripts/step_22_6l_router_rag_persistenza.py --rollback
"""

import argparse
import ast
import csv
import hashlib
import json
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

API_ORIGINAL = ROOT / "api/main.py"

PATCH_DIR = ROOT / "docker/runtime_v32_patch"

API_PATCH = PATCH_DIR / "api_main_rag_v32.py"

OVERRIDE = (
    ROOT / "docker-compose.primetta-v32-rag-persistenza.yaml"
)

PERSISTENT = ROOT / "data/user_data"

RUNS = ROOT / "target_prediction_v33_vscode/runs"

SERVICE = "maize-ncrna-gnn"

CONTAINER = "maize_ncrna_gnn"

IMAGE = "ai_primetta_v32_full_ai_sklearn18:22.6g"

COMPOSE_FILES = [
    ROOT / "docker-compose.yml",
    ROOT / "docker-compose.primetta-v32-8503.yaml",
    ROOT / "docker-compose.primetta-v32-full-ai.yaml",
    ROOT / "docker-compose.primetta-v32-sklearn18.yaml",
    ROOT / "docker-compose.primetta-v32-scripts.yaml",
]

FDR_QUESTION = (
    "Quali ncRNA caricati dal ricercatore "
    "sono presenti nel ranking e hanno "
    "FDR inferiore a 0.05?"
)

CONTROL_QUESTION = (
    "Quali sono i primi cinque candidati "
    "del ranking esplorativo?"
)


# ============================================================
# 2. FUNZIONE SCIENTIFICA DA AGGIUNGERE ALLA COPIA API
# ============================================================

HELPER = '''
def _v32_answer_user_fdr(question: str, top_k: int):

    """
    Interrogazione deterministica dei dati del ricercatore.

    Legge il CSV integrato.
    Non inventa FDR.
    Non utilizza lo score ML come FDR.
    Non esegue training.
    """

    import csv
    import math
    import re

    from pathlib import Path

    q = question.casefold()

    if not (
        "fdr" in q
        or "padj" in q
        or "q-value" in q
    ):
        return None

    if not any(
        word in q
        for word in (
            "ricercator",
            "caricat",
            "utente",
            "user"
        )
    ):
        return None

    value_match = re.search(
        r"(?:fdr|padj|q-value)[^0-9]{0,90}"
        r"(0[.,][0-9]+|1(?:[.,]0+)?)",
        q
    )

    threshold = (
        float(
            value_match.group(1).replace(",", ".")
        )
        if value_match
        else 0.05
    )

    if not 0 <= threshold <= 1:

        raise ValueError(
            "Soglia FDR non valida"
        )

    greater = any(
        s in q
        for s in (
            "superiore",
            "maggiore",
            ">="
        )
    )

    inclusive = any(
        s in q
        for s in (
            "<=",
            ">=",
            "uguale",
            "al massimo"
        )
    )

    if greater:

        operator = (
            ">=" if inclusive else ">"
        )

    else:

        operator = (
            "<=" if inclusive else "<"
        )

    # --------------------------------------------------------
    # INDIVIDUAZIONE DATASET INTEGRATO
    # --------------------------------------------------------

    folder = (
        BASE_DIR / "outputs/user_data"
    )

    candidates = list(
        folder.glob(
            "*/ranking_with_user_evidence.csv"
        )
    )

    if not candidates:

        raise RuntimeError(
            "Nessun ranking integrato "
            "del ricercatore disponibile"
        )

    source = max(
        candidates,
        key=lambda p: (
            p.stat().st_mtime_ns,
            p.name
        )
    )

    selected = []

    # --------------------------------------------------------
    # LETTURA CSV
    # --------------------------------------------------------

    with source.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        reader = csv.DictReader(
            stream
        )

        required = {
            "rnacentral_id",
            "exploratory_rank",
            "user_min_FDR"
        }

        if not required.issubset(
            set(reader.fieldnames or [])
        ):

            raise ValueError(
                "Colonne indispensabili "
                "assenti nel ranking integrato"
            )

        # ----------------------------------------------------
        # FILTRO STATISTICO
        # ----------------------------------------------------

        for row in reader:

            raw = (
                row.get(
                    "user_min_FDR"
                ) or ""
            ).strip().replace(",", ".")

            try:

                fdr = float(raw)

            except ValueError:

                continue

            if (
                not math.isfinite(fdr)
                or not 0 <= fdr <= 1
            ):

                continue

            present = str(
                row.get(
                    "user_data_present",
                    ""
                )
            ).strip().lower()

            record_count = str(
                row.get(
                    "user_records_count",
                    ""
                )
            ).strip()

            if (
                present
                and present not in (
                    "true",
                    "1",
                    "yes",
                    "si"
                )
            ):

                continue

            if record_count:

                try:

                    if float(record_count) <= 0:

                        continue

                except ValueError:

                    continue

            # ------------------------------------------------
            # APPLICAZIONE OPERATORE
            # ------------------------------------------------

            if operator == ">=":

                passed = fdr >= threshold

            elif operator == ">":

                passed = fdr > threshold

            elif operator == "<=":

                passed = fdr <= threshold

            else:

                passed = fdr < threshold

            if not passed:

                continue

            # ------------------------------------------------
            # IDENTIFICATIVO RNAcentral
            # ------------------------------------------------

            identifier = str(
                row.get(
                    "rnacentral_id"
                ) or ""
            ).strip().upper()

            if not re.fullmatch(
                r"URS[0-9A-F]{10,12}",
                identifier
            ):

                continue

            selected.append({

                "rnacentral_id": identifier,

                "rank": row.get(
                    "exploratory_rank"
                ),

                "user_min_FDR": fdr,

                "rna_type": row.get(
                    "rna_type"
                ),

                "exploratory_rank_score": row.get(
                    "exploratory_rank_score"
                ),

                "user_source_files": row.get(
                    "user_source_files"
                )
            })

    # --------------------------------------------------------
    # ORDINAMENTO PER RANK
    # --------------------------------------------------------

    selected.sort(
        key=lambda row: (

            int(float(row["rank"]))
            if str(row["rank"]).strip()
            else 10**9,

            row["rnacentral_id"]
        )
    )

    limit = min(
        max(int(top_k), 1),
        50
    )

    display = selected[:limit]

    # --------------------------------------------------------
    # COSTRUZIONE RISPOSTA
    # --------------------------------------------------------

    lines = [

        f"ncRNA del ricercatore presenti "
        f"nel ranking con FDR "
        f"{operator} {threshold:g}: "
        f"{len(selected)} "
        f"(visualizzati {len(display)})."
    ]

    for row in display:

        lines.append(

            f"- {row['rnacentral_id']} | "
            f"rank={row['rank']} | "
            f"FDR={row['user_min_FDR']:.6g} | "
            f"tipo={row['rna_type']}"
        )

    lines.append(

        "I valori FDR provengono dai dati "
        "integrati del ricercatore, non dagli "
        "score del ranking. L'associazione "
        "non equivale a validazione "
        "biologica indipendente."
    )

    # --------------------------------------------------------
    # OUTPUT CON PROVENIENZA SCIENTIFICA
    # --------------------------------------------------------

    return {

        "question": question,

        "provider": "deterministic_user_fdr_v32",

        "answer": "\\n".join(lines),

        "total_matches": len(selected),

        "filters": {

            "user_FDR_operator": operator,

            "threshold": threshold
        },

        "evidence": display,

        "sources": [

            {

                "source": str(
                    source.relative_to(BASE_DIR)
                ),

                "source_type": (
                    "user_ranking_integrated_csv"
                ),

                "rows_returned": len(display)
            }
        ]
    }

'''


# ============================================================
# 3. ROUTER RAG CORRETTO
# ============================================================

ENDPOINT = '''def rag_ask(payload: RagAskRequest) -> Dict[str, Any]:

    if not payload.question.strip():

        raise HTTPException(
            status_code=400,
            detail="Domanda vuota."
        )

    try:

        user_result = _v32_answer_user_fdr(
            payload.question,
            payload.top_k
        )

        if user_result is not None:

            return user_result

        return get_rag_answer(

            question=payload.question,

            provider=payload.provider,

            top_k=payload.top_k
        )

    except Exception as exc:

        raise HTTPException(
            status_code=500,
            detail=str(exc)
        ) from exc

'''


# ============================================================
# 4. ESECUZIONE COMANDI
# ============================================================

def run(
    args,
    timeout=120,
    check=True
):

    proc = subprocess.run(

        list(map(str, args)),

        text=True,

        capture_output=True,

        timeout=timeout,

        check=False
    )

    if check and proc.returncode:

        raise RuntimeError(

            "Comando fallito: "

            + " ".join(map(str, args))

            + "\n"

            + proc.stderr[-1300:]
        )

    return proc


# ============================================================
# 5. GESTIONE DOCKER COMPOSE
# ============================================================

def compose(patched=False):

    files = COMPOSE_FILES + (
        [OVERRIDE] if patched else []
    )

    for path in files:

        if not path.is_file():

            raise FileNotFoundError(
                path
            )

    command = [
        "docker",
        "compose"
    ]

    for path in files:

        command.extend([
            "-f",
            str(path)
        ])

    return command


# ============================================================
# 6. VALIDAZIONE CONFIGURAZIONE
# ============================================================

def configuration(patched=False):

    response = run(
        compose(patched)
        + [
            "config",
            "--format",
            "json"
        ]
    )

    data = json.loads(
        response.stdout
    )

    service = data[
        "services"
    ][SERVICE]

    if service.get("image") != IMAGE:

        raise RuntimeError(
            "Immagine diversa da 22.6g"
        )

    ports = service.get(
        "ports",
        []
    )

    valid = any(

        str(p.get("published")) == "8503"

        and str(p.get("target")) == "8501"

        and p.get("host_ip") == "127.0.0.1"

        for p in ports
    )

    if not valid:

        raise RuntimeError(
            "Porta 8503 non conforme"
        )

    if any(
        str(p.get("published")) == "8501"
        for p in ports
    ):

        raise RuntimeError(
            "Il container non deve "
            "pubblicare la porta 8501"
        )

    return service


# ============================================================
# 7. IMMAGINE ATTIVA
# ============================================================

def active_image():

    response = run(

        [
            "docker",
            "inspect",
            CONTAINER,
            "--format",
            "{{.Config.Image}}"
        ],

        check=False
    )

    return (

        response.stdout.strip()

        if response.returncode == 0

        else None
    )


# ============================================================
# 8. CHECKSUM DEI FILE ORIGINALI
# ============================================================

def digest(path):

    h = hashlib.sha256()

    with path.open("rb") as stream:

        for part in iter(
            lambda: stream.read(
                1024 * 1024
            ),
            b""
        ):

            h.update(part)

    return h.hexdigest()


def guarded_files():

    files = [

        API_ORIGINAL,

        ROOT / "app/streamlit_app.py",

        *sorted(
            (ROOT / "app/pages").glob("*.py")
        ),

        ROOT
        / "models/graph_autoencoder_v3/"
        "gae_v3_model.pt",

        ROOT
        / "data/processed/pyg/"
        "maize_ncrna_graph_v3_prjna419326.pt",

        ROOT
        / "models/srna_supervised_ml_v1_antileakage/"
        "PRJNA419326_srna_supervised_antileakage_best_model.joblib"
    ]

    if any(
        not p.is_file()
        for p in files
    ):

        raise RuntimeError(
            "Mancano file scientifici "
            "o GUI da proteggere"
        )

    return {

        str(p.relative_to(ROOT)): digest(p)

        for p in files
    }


# ============================================================
# 9. GENERAZIONE DELLA COPIA API CORRETTA
# ============================================================

def generate_patch():

    content = API_ORIGINAL.read_text(
        encoding="utf-8"
    )

    tree = ast.parse(
        content
    )

    matches = []

    for node in tree.body:

        if not isinstance(
            node,
            (
                ast.FunctionDef,
                ast.AsyncFunctionDef
            )
        ):

            continue

        if node.name != "rag_ask":

            continue

        for decorator in node.decorator_list:

            if (
                isinstance(
                    decorator,
                    ast.Call
                )
                and decorator.args
                and isinstance(
                    decorator.args[0],
                    ast.Constant
                )
            ):

                if (
                    decorator.args[0].value
                    == "/rag/ask"
                ):

                    matches.append(node)

    if len(matches) != 1:

        raise RuntimeError(
            "Router /rag/ask non univoco"
        )

    node = matches[0]

    old = ast.get_source_segment(
        content,
        node
    ) or ""

    if "get_rag_answer" not in old:

        raise RuntimeError(
            "Il router originale non coincide "
            "con quello precedentemente analizzato"
        )

    # --------------------------------------------------------
    # SOSTITUZIONE SOLO NELLA COPIA
    # --------------------------------------------------------

    lines = content.splitlines(
        keepends=True
    )

    prefix = "".join(
        lines[:node.lineno - 1]
    )

    suffix = "".join(
        lines[node.end_lineno:]
    )

    new_source = (
        prefix
        + ENDPOINT
        + suffix
    )

    # --------------------------------------------------------
    # INSERIMENTO FUNZIONE SCIENTIFICA
    # --------------------------------------------------------

    decorator_line = min(
        d.lineno
        for d in node.decorator_list
    )

    patched_lines = new_source.splitlines(
        keepends=True
    )

    patched_lines.insert(
        decorator_line - 1,
        HELPER
    )

    patched = "".join(
        patched_lines
    )

    compile(
        patched,
        str(API_PATCH),
        "exec"
    )

    return patched


# ============================================================
# 10. SCRITTURA SICURA DEI FILE GENERATI
# ============================================================

def write_equal_or_new(path, text):

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if path.exists():

        if path.read_text(
            encoding="utf-8"
        ) != text:

            raise RuntimeError(
                "File esistente con contenuto "
                "diverso: "
                + str(path)
            )

    else:

        path.write_text(
            text,
            encoding="utf-8"
        )

        print(
            "[CREATO]",
            path.relative_to(ROOT)
        )


# ============================================================
# 11. BACKUP DEGLI UPLOAD ORIGINALI
# ============================================================

def snapshot_uploads(backup):

    backup.mkdir(
        parents=True,
        exist_ok=False
    )

    exists = run(

        [
            "docker",
            "exec",
            CONTAINER,
            "test",
            "-d",
            "/app/data/user_data"
        ],

        check=False
    )

    if exists.returncode:

        print(
            "[UPLOAD] Nessuna cartella "
            "nel container da copiare"
        )

        return

    # --------------------------------------------------------
    # BACKUP DAL CONTAINER AL DISCO
    # --------------------------------------------------------

    run(

        [
            "docker",
            "cp",

            (
                f"{CONTAINER}:"
                "/app/data/user_data/."
            ),

            str(backup)
        ],

        timeout=300
    )

    print(
        "[BACKUP UPLOAD]",
        backup
    )

    PERSISTENT.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # CONSOLIDAMENTO SENZA SOVRASCRITTURA
    # --------------------------------------------------------

    for source in backup.rglob("*"):

        if source.is_symlink():

            raise RuntimeError(
                "Link simbolico non ammesso "
                "nel backup degli upload"
            )

        relative = source.relative_to(
            backup
        )

        target = PERSISTENT / relative

        if source.is_dir():

            target.mkdir(
                parents=True,
                exist_ok=True
            )

        elif source.is_file():

            target.parent.mkdir(
                parents=True,
                exist_ok=True
            )

            if target.exists():

                if (
                    not target.is_file()

                    or digest(target) != digest(source)
                ):

                    raise RuntimeError(
                        "Conflitto nei file upload: "
                        + str(relative)
                    )

            else:

                shutil.copy2(
                    source,
                    target
                )


# ============================================================
# 12. HEALTHCHECK HTTP
# ============================================================

def http(port, path):

    url = (
        f"http://127.0.0.1:{port}{path}"
    )

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    try:

        with opener.open(
            url,
            timeout=5
        ) as response:

            return response.status == 200

    except Exception:

        return False


# ============================================================
# 13. INTERROGAZIONE API
# ============================================================

def query(question):

    request = urllib.request.Request(

        "http://127.0.0.1:8000/rag/ask",

        data=json.dumps({

            "question": question,

            "top_k": 5

        }).encode("utf-8"),

        headers={
            "Content-Type": "application/json"
        },

        method="POST"
    )

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    try:

        with opener.open(
            request,
            timeout=45
        ) as response:

            return json.loads(
                response.read(300000)
            )

    except urllib.error.HTTPError as exc:

        raise RuntimeError(

            f"RAG HTTP {exc.code}: "

            + repr(
                exc.read(1000)
            )

        ) from exc


# ============================================================
# 14. AVVIO DEL SOLO CONTAINER V3.2
# ============================================================

def start(patched):

    configuration(
        patched
    )

    run(

        compose(patched) + [

            "up",
            "-d",

            "--no-deps",

            "--no-build",

            "--pull",
            "never",

            "--force-recreate",

            SERVICE
        ],

        timeout=180
    )

    for _ in range(25):

        ready = (

            active_image() == IMAGE

            and http(
                8503,
                "/_stcore/health"
            )

            and http(
                8000,
                "/health"
            )
        )

        if ready:

            return

        time.sleep(2)

    raise RuntimeError(
        "La GUI V3.2 o la sua API "
        "non rispondono dopo il riavvio"
    )


# ============================================================
# 15. STATO APPLICAZIONI
# ============================================================

def status():

    print(
        "[IMMAGINE V3.2]",
        active_image()
    )

    endpoints = [

        (
            "GUI V3.2",
            8503,
            "/_stcore/health"
        ),

        (
            "API V3.2",
            8000,
            "/health"
        ),

        (
            "GUI V3.3",
            8502,
            "/_stcore/health"
        ),

        (
            "AI ORZO",
            8501,
            "/_stcore/health"
        )
    ]

    for name, port, endpoint in endpoints:

        success = http(
            port,
            endpoint
        )

        print(

            "[OK]" if success else "[ERRORE]",

            name
        )

    # --------------------------------------------------------
    # CONTROLLO MONTAGGI
    # --------------------------------------------------------

    volumes = run(

        [
            "docker",
            "inspect",
            CONTAINER,
            "--format",
            "{{json .Mounts}}"
        ],

        check=False
    )

    if volumes.returncode == 0:

        mounts = json.loads(
            volumes.stdout
        )

        for target in (

            "/app/api/main.py",

            "/app/data/user_data"
        ):

            present = any(

                item.get(
                    "Destination"
                ) == target

                for item in mounts
            )

            print(
                "[MOUNT]",
                target,
                present
            )


# ============================================================
# 16. INSTALLAZIONE E COLLAUDO
# ============================================================

def install():

    # --------------------------------------------------------
    # A. VERIFICA AMBIENTE PRECEDENTE
    # --------------------------------------------------------

    configuration(False)

    if active_image() != IMAGE:

        raise RuntimeError(
            "Immagine V3.2 inattesa: "
            "intervento bloccato"
        )

    before = guarded_files()

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    report_dir = RUNS / (
        "STEP22_6L_RAG_UPLOAD_" + stamp
    )

    report_dir.mkdir(
        parents=True,
        exist_ok=False
    )

    print(
        "\n[1] BACKUP UPLOAD DEL RICERCATORE"
    )

    # --------------------------------------------------------
    # B. BACKUP PRIMA DEI NUOVI VOLUMI
    # --------------------------------------------------------

    backup_dir = (
        report_dir
        / "backup_uploads_pre_mount"
    )

    snapshot_uploads(
        backup_dir
    )

    PERSISTENT.mkdir(
        parents=True,
        exist_ok=True
    )

    print(
        "[OK] Cartella persistente:",
        PERSISTENT
    )

    # --------------------------------------------------------
    # C. PREPARAZIONE API CORRETTA
    # --------------------------------------------------------

    print(
        "\n[2] GENERAZIONE API SCIENTIFICA"
    )

    write_equal_or_new(

        API_PATCH,

        generate_patch()
    )

    # --------------------------------------------------------
    # D. OVERRIDE COMPOSE
    # --------------------------------------------------------

    override = {

        "services": {

            SERVICE: {

                "volumes": [

                    {

                        "type": "bind",

                        "source": str(
                            API_PATCH
                        ),

                        "target": "/app/api/main.py",

                        "read_only": True
                    },

                    {

                        "type": "bind",

                        "source": str(
                            PERSISTENT
                        ),

                        "target": (
                            "/app/data/user_data"
                        )
                    }
                ]
            }
        }
    }

    write_equal_or_new(

        OVERRIDE,

        json.dumps(
            override,
            indent=2
        ) + "\n"
    )

    # --------------------------------------------------------
    # E. VERIFICA VOLUMI COMPOSE
    # --------------------------------------------------------

    cfg = configuration(True)

    mounts = {

        v.get("target"): v

        for v in cfg.get(
            "volumes", []
        )
    }

    expected = {

        "/app/api/main.py": API_PATCH,

        "/app/data/user_data": PERSISTENT
    }

    for target, source in expected.items():

        if target not in mounts:

            raise RuntimeError(
                "Volume Docker mancante: "
                + target
            )

        configured = Path(
            mounts[target]["source"]
        ).resolve()

        if configured != source.resolve():

            raise RuntimeError(
                "Sorgente Docker inattesa: "
                + target
            )

    if guarded_files() != before:

        raise RuntimeError(
            "File originali modificati "
            "prima dell'attivazione"
        )

    # --------------------------------------------------------
    # F. PREFLIGHT API IN CONTAINER ISOLATO
    # --------------------------------------------------------

    print(
        "\n[3] TEST IMPORT API SENZA ATTIVAZIONE"
    )

    preflight = run(

        [
            "docker",
            "run",
            "--rm",

            "--pull",
            "never",

            "--network",
            "none",

            "--entrypoint",
            "python3",

            "-w",
            "/app",

            "-v",
            f"{API_PATCH}:/app/api/main.py:ro",

            "-e",
            "PYTHONDONTWRITEBYTECODE=1",

            IMAGE,

            "-B",
            "-c",

            (
                "import api.main;"
                "assert callable("
                "api.main._v32_answer_user_fdr"
                ");"
                "print('API_PATCH_IMPORT_OK')"
            )
        ],

        timeout=90
    )

    print(
        preflight.stdout.strip()
    )

    # --------------------------------------------------------
    # G. ATTIVAZIONE E VERIFICA
    # --------------------------------------------------------

    print(
        "\n[4] ATTIVAZIONE GUI V3.2"
    )

    try:

        start(True)

        print(
            "[OK] GUI e API riavviate"
        )

        # ----------------------------------------------------
        # TEST RAG SCIENTIFICO
        # ----------------------------------------------------

        print(
            "\n[5] TEST DOMANDA FDR"
        )

        fdr = query(
            FDR_QUESTION
        )

        print(
            "\n[6] TEST DOMANDA GENERICA"
        )

        generic = query(
            CONTROL_QUESTION
        )

        # ----------------------------------------------------
        # CONTROLLO ROUTING
        # ----------------------------------------------------

        if (
            fdr.get("provider")
            != "deterministic_user_fdr_v32"
        ):

            raise RuntimeError(
                "La domanda FDR non utilizza "
                "il router scientifico"
            )

        if (
            generic.get("provider")
            == "deterministic_user_fdr_v32"
        ):

            raise RuntimeError(
                "Una domanda generica è stata "
                "intercettata erroneamente"
            )

        # ----------------------------------------------------
        # CONTROLLO QUATTRO RNA DOCUMENTATI
        # ----------------------------------------------------

        if fdr.get("total_matches") != 4:

            raise RuntimeError(
                "Attese quattro corrispondenze "
                "FDR nel dataset documentato"
            )

        actual = {

            row["rnacentral_id"]

            for row in fdr.get(
                "evidence", []
            )
        }

        expected_ids = {

            "URS0003C75441",

            "URS0003C7543B",

            "URS0003C753CA",

            "URS0003C753C9"
        }

        if actual != expected_ids:

            raise RuntimeError(
                "Gli identificativi ncRNA "
                "non corrispondono ai dati "
                "precedentemente documentati"
            )

        # ----------------------------------------------------
        # VERIFICA MONTAGGIO PERSISTENTE
        # ----------------------------------------------------

        volume_response = run([

            "docker",
            "inspect",
            CONTAINER,
            "--format",
            "{{json .Mounts}}"
        ])

        live_mounts = {

            item.get("Destination"): item

            for item in json.loads(
                volume_response.stdout
            )
        }

        upload_mount = live_mounts.get(
            "/app/data/user_data"
        )

        if (
            not upload_mount
            or not upload_mount.get("RW")
        ):

            raise RuntimeError(
                "Volume upload persistente "
                "non montato in scrittura"
            )

        # ----------------------------------------------------
        # VERIFICA CHECKSUM
        # ----------------------------------------------------

        if guarded_files() != before:

            raise RuntimeError(
                "Integrità dei file originali "
                "non rispettata"
            )

    except Exception:

        print(
            "[ERRORE] Tentativo automatico "
            "di ripristino configurazione precedente"
        )

        try:

            start(False)

        except Exception as rollback_exc:

            print(
                "[ATTENZIONE] Rollback non riuscito:",
                rollback_exc
            )

        raise

    # --------------------------------------------------------
    # H. RISULTATI
    # --------------------------------------------------------

    print(
        "\n[RAG FDR]"
    )

    print(
        fdr["answer"]
    )

    print(
        "\n[RAG GENERICO]"
    )

    print(
        generic.get("provider")
    )

    # --------------------------------------------------------
    # I. RAPPORTO SCIENTIFICO
    # --------------------------------------------------------

    report = {

        "step": "22.6L",

        "status": "PASS",

        "rag_fdr": fdr,

        "rag_generic_provider": generic.get(
            "provider"
        ),

        "backup_uploads": str(
            backup_dir
        ),

        "persistent_uploads": str(
            PERSISTENT
        ),

        "original_sha256": before,

        "training_performed": False,

        "gui_modified": False
    }

    destination = (
        report_dir
        / "RAPPORTO_STEP22_6L.json"
    )

    destination.write_text(

        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),

        encoding="utf-8"
    )

    print()
    print("=" * 65)

    print(
        "[OK] Routing RAG FDR attivato"
    )

    print(
        "[OK] Persistenza upload configurata"
    )

    print(
        "[OK] Nessun training"
    )

    print(
        "[OK] GUI originale preservata"
    )

    print(
        "[REPORT]",
        destination
    )

    print("=" * 65)

    status()


# ============================================================
# 17. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    group = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    for mode in (
        "install",
        "status",
        "rollback"
    ):

        group.add_argument(
            "--" + mode,
            action="store_true"
        )

    args = parser.parse_args()

    print("=" * 65)

    print(
        "AI PRIMETTA V3.2 | STEP 22.6L"
    )

    print(
        "RAG SCIENTIFICO E PERSISTENZA "
        "DATI RICERCATORE"
    )

    print("=" * 65)

    if args.install:

        install()

    elif args.rollback:

        start(False)
        status()

    else:

        status()


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