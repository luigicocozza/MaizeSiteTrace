#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AI PRIMETTA V3.2 - STEP 22.6O

--audit: verifica collegamenti e scritture senza modifiche
--install: salva embedding, isola gli output del trainer e applica override Docker
--status: stato del container
--rollback: rimuove override di isolamento dal servizio attivo

Non avvia il training e non modifica GUI o trainer originale.
"""

import argparse
import ast
import hashlib
import json
import subprocess
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TRAINER = ROOT / "25_train_graph_autoencoder.py"
PAGE = ROOT / "app/pages/01_Rete_neurale_GAE.py"
API = ROOT / "api/main.py"
LEGACY = ROOT / "models/graph_autoencoder/gae_model.pt"
V3 = ROOT / "models/graph_autoencoder_v3/gae_v3_model.pt"
EMB = ROOT / "data/processed/embeddings"
EMBEDDINGS = [
    EMB / "gae_node_embeddings.csv",
    EMB / "gae_ncrna_embeddings.csv"
]
RUNS = ROOT / "target_prediction_v33_vscode/runs"
PATCH = ROOT / "docker/runtime_v32_patch/25_train_graph_autoencoder_safe.py"
OVERRIDE = ROOT / "docker-compose.primetta-v32-gae-safe.yaml"

BASE_FILES = [
    ROOT / name for name in (
        "docker-compose.yml",
        "docker-compose.primetta-v32-8503.yaml",
        "docker-compose.primetta-v32-full-ai.yaml",
        "docker-compose.primetta-v32-sklearn18.yaml",
        "docker-compose.primetta-v32-scripts.yaml",
        "docker-compose.primetta-v32-rag-persistenza.yaml",
    )
]

SERVICE = "maize-ncrna-gnn"
CONTAINER = "maize_ncrna_gnn"
IMAGE = "ai_primetta_v32_full_ai_sklearn18:22.6g"
DEST = "/app/25_train_graph_autoencoder.py"


# ============================================================
# 1. COMANDI
# ============================================================

def cmd(argv, timeout=130, check=True):

    p = subprocess.run(
        [str(x) for x in argv],
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False
    )

    if check and p.returncode:

        raise RuntimeError(
            f"Comando fallito ({p.returncode}): "
            f"{' '.join(map(str, argv[:7]))}\n"
            f"{p.stderr[-1300:]}"
        )

    return p


# ============================================================
# 2. CHECKSUM
# ============================================================

def sha(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def originals():

    paths = [
        TRAINER,
        PAGE,
        API,
        LEGACY,
        V3
    ] + [
        p for p in EMBEDDINGS if p.exists()
    ]

    for path in paths:

        if not path.is_file() or path.is_symlink():

            raise RuntimeError(
                f"File originale non valido: {path}"
            )

    return {
        str(p.relative_to(ROOT)): sha(p)
        for p in paths
    }


# ============================================================
# 3. CONFIGURAZIONE DOCKER COMPOSE
# ============================================================

def compose(patched=False):

    paths = BASE_FILES + (
        [OVERRIDE] if patched else []
    )

    for p in paths:

        if not p.is_file():
            raise FileNotFoundError(p)

    out = ["docker", "compose"]

    for p in paths:
        out += ["-f", str(p)]

    return out


def config(patched=False):

    response = cmd(
        compose(patched) + [
            "config",
            "--format",
            "json"
        ]
    )

    obj = json.loads(response.stdout)
    service = obj["services"][SERVICE]

    if service.get("image") != IMAGE:

        raise RuntimeError(
            "Immagine Docker inattesa"
        )

    ports = service.get("ports", [])

    valid_port = any(
        str(p.get("published")) == "8503"
        and str(p.get("target")) == "8501"
        and p.get("host_ip") == "127.0.0.1"
        for p in ports
    )

    if not valid_port:

        raise RuntimeError(
            "Manca la porta GUI 127.0.0.1:8503"
        )

    if any(
        str(p.get("published")) == "8501"
        for p in ports
    ):

        raise RuntimeError(
            "Conflitto con AI ORZO"
        )

    volumes = {
        v["target"]: v
        for v in service.get("volumes", [])
    }

    output = volumes.get("/app/outputs")

    if (
        not output
        or Path(output["source"]).resolve()
        != (ROOT / "outputs").resolve()
        or output.get("read_only", False)
    ):

        raise RuntimeError(
            "/app/outputs deve essere "
            "un volume persistente RW"
        )

    if patched:

        mount = volumes.get(DEST)

        if (
            not mount
            or Path(mount["source"]).resolve()
            != PATCH.resolve()
            or not mount.get("read_only", False)
        ):

            raise RuntimeError(
                "Mount di protezione trainer "
                "non conforme"
            )

    return service


# ============================================================
# 4. CONTAINER ATTIVO
# ============================================================

def current():

    p = cmd(
        [
            "docker",
            "inspect",
            CONTAINER,
            "--format",
            "{{json .}}"
        ],
        check=False
    )

    if p.returncode:
        return {}

    c = json.loads(p.stdout)

    mounts = {
        m["Destination"]: m
        for m in c.get("Mounts", [])
    }

    return {
        "running": c.get(
            "State", {}
        ).get("Running", False),

        "image": c.get(
            "Config", {}
        ).get("Image"),

        "trainer_mount": mounts.get(
            DEST, {}
        ).get("Source"),

        "upload_persistent": (
            "/app/data/user_data" in mounts
        )
    }


# ============================================================
# 5. ANALISI STATICA DEL TRAINER
# ============================================================

def analyze_trainer():

    source = TRAINER.read_text(
        encoding="utf-8"
    )

    tree = ast.parse(
        source,
        filename=str(TRAINER)
    )

    definitions = {}

    for node in tree.body:

        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
        ):

            definitions[
                node.targets[0].id
            ] = node

    required = (
        "OUT_DIR",
        "MODEL_DIR",
        "EMB_DIR"
    )

    if not all(
        n in definitions
        for n in required + ("BASE_DIR",)
    ):

        raise RuntimeError(
            "Costanti di percorso mancanti "
            "nel trainer"
        )

    if not all(
        definitions[n].lineno
        == definitions[n].end_lineno
        for n in required
    ):

        raise RuntimeError(
            "Definizioni output multilinea: "
            "stop precauzionale"
        )

    first = min(
        definitions[n].lineno
        for n in required
    )

    if definitions["BASE_DIR"].lineno >= first:

        raise RuntimeError(
            "BASE_DIR definita dopo "
            "le directory output"
        )

    # Controlla che ogni destinazione di scrittura
    # dipenda dalle directory che saranno isolate.

    def depends(node, depth=0):

        if depth > 12:
            return False

        if isinstance(node, ast.Name):

            if node.id in required:
                return True

            nxt = definitions.get(node.id)

            return bool(
                nxt and depends(
                    nxt.value,
                    depth + 1
                )
            )

        if (
            isinstance(node, ast.BinOp)
            and isinstance(node.op, ast.Div)
        ):

            return depends(
                node.left,
                depth + 1
            )

        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "joinpath"
        ):

            return depends(
                node.func.value,
                depth + 1
            )

        return False

    writes = []

    for node in ast.walk(tree):

        if not isinstance(node, ast.Call):
            continue

        fn = ast.unparse(node.func)

        target = None

        if (
            fn in (
                "torch.save",
                "joblib.dump",
                "pickle.dump"
            )
            and len(node.args) >= 2
        ):

            target = node.args[1]

        elif (
            fn.split(".")[-1] in (
                "to_csv",
                "to_excel",
                "savefig",
                "to_parquet",
                "to_json"
            )
            and node.args
        ):

            target = node.args[0]

        elif (
            fn.split(".")[-1] in (
                "write_text",
                "write_bytes"
            )
            and isinstance(
                node.func,
                ast.Attribute
            )
        ):

            target = node.func.value

        elif (
            fn in ("open", "Path.open")
            and len(node.args) >= 2
        ):

            mode = node.args[1]

            if (
                isinstance(mode, ast.Constant)
                and isinstance(mode.value, str)
                and any(
                    c in mode.value
                    for c in "wax+"
                )
            ):

                target = node.args[0]

        if target is not None:

            safe = depends(target)

            writes.append({
                "line": node.lineno,
                "function": fn,
                "destination": ast.unparse(target),
                "isolated": safe
            })

    if (
        not writes
        or any(
            not w["isolated"]
            for w in writes
        )
    ):

        raise RuntimeError(
            "Scritture non verificabili: "
            + json.dumps(
                writes,
                ensure_ascii=False
            )
        )

    return (
        source,
        definitions,
        writes,
        first
    )


# ============================================================
# 6. GENERAZIONE DEL TRAINER ISOLATO
# ============================================================

def patch_source():

    source, definitions, writes, first = (
        analyze_trainer()
    )

    # Ogni esecuzione genera una directory nuova.
    # Il codice originale del training viene mantenuto.

    insert = (
        "from datetime import datetime as _primetta_dt_safe\n"
        "from uuid import uuid4 as _primetta_uuid_safe\n"
        "__primetta_gae_safe_root = "
        "(BASE_DIR / 'outputs' / 'gae_gui_safe_runs' / "
        "(_primetta_dt_safe.now().strftime("
        "'%Y%m%dT%H%M%S%f') + '_' + "
        "_primetta_uuid_safe().hex[:8]))\n"
    )

    replacements = {
        definitions["OUT_DIR"].lineno:
            "OUT_DIR = "
            "__primetta_gae_safe_root / 'results'\n",

        definitions["MODEL_DIR"].lineno:
            "MODEL_DIR = "
            "__primetta_gae_safe_root / 'models'\n",

        definitions["EMB_DIR"].lineno:
            "EMB_DIR = "
            "__primetta_gae_safe_root / 'embeddings'\n",
    }

    lines = source.splitlines(
        keepends=True
    )

    out = []

    for n, line in enumerate(
        lines,
        start=1
    ):

        if n == first:
            out.append(insert)

        out.append(
            replacements.get(n, line)
        )

    patched = "".join(out)

    compile(
        patched,
        str(PATCH),
        "exec"
    )

    return patched, writes


# ============================================================
# 7. SCRITTURA NON DISTRUTTIVA
# ============================================================

def write_once(path, content):

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    if path.exists():

        if path.read_text(
            encoding="utf-8"
        ) != content:

            raise RuntimeError(
                "Patch preesistente differente: "
                + str(path)
            )

    else:

        path.write_text(
            content,
            encoding="utf-8"
        )

        print(
            "[CREATO]",
            path.relative_to(ROOT)
        )


# ============================================================
# 8. CONTROLLO HTTP
# ============================================================

def http(port, route):

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    try:

        with opener.open(
            f"http://127.0.0.1:{port}{route}",
            timeout=4
        ) as response:

            return response.status == 200

    except Exception:
        return False


# ============================================================
# 9. RIAVVIO DEL SOLO SERVIZIO V3.2
# ============================================================

def start(patched):

    config(patched)

    cmd(
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
        timeout=160
    )

    for _ in range(25):

        state = current()

        if (
            state.get("image") == IMAGE
            and http(
                8503,
                "/_stcore/health"
            )
            and http(
                8000,
                "/health"
            )
        ):
            return

        time.sleep(2)

    raise RuntimeError(
        "GUI o API V3.2 non raggiungibili "
        "dopo l'avvio"
    )


# ============================================================
# 10. STATO GENERALE
# ============================================================

def status():

    state = current()

    print(
        "[IMMAGINE]",
        state.get("image")
    )

    print(
        "[TRAINER IN CONTAINER]",
        state.get("trainer_mount")
    )

    print(
        "[PERSISTENZA UPLOAD]",
        state.get("upload_persistent")
    )

    services = (
        ("GUI V3.2", 8503, "/_stcore/health"),
        ("API V3.2", 8000, "/health"),
        ("GUI V3.3", 8502, "/_stcore/health"),
        ("AI ORZO", 8501, "/_stcore/health")
    )

    for name, port, path in services:

        print(
            "[OK]"
            if http(port, path)
            else "[NON VERIFICATO]",
            name
        )


# ============================================================
# 11. BACKUP EMBEDDING GAE ORIGINALI
# ============================================================

def backup_embeddings(folder):

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    entries = {}

    for source in EMBEDDINGS:

        if not source.exists():

            entries[source.name] = {
                "exists": False
            }

            continue

        if (
            not source.is_file()
            or source.is_symlink()
        ):

            raise RuntimeError(
                "Embedding non regolare: "
                + str(source)
            )

        dest = folder / source.name

        with (
            source.open("rb") as input_file,
            dest.open("xb") as output_file
        ):

            while True:

                block = input_file.read(
                    1024 * 1024
                )

                if not block:
                    break

                output_file.write(block)

        if sha(source) != sha(dest):

            raise RuntimeError(
                "Backup embedding discordante: "
                + source.name
            )

        entries[source.name] = {
            "exists": True,
            "size": source.stat().st_size,
            "sha256": sha(source),
            "backup": str(dest)
        }

        print(
            "[BACKUP OK]",
            source.name
        )

    return entries


# ============================================================
# 12. AUDIT PREVENTIVO
# ============================================================

def audit():

    config(False)

    print(
        "[CONTAINER]",
        current()
    )

    _, _, writes, _ = analyze_trainer()

    for row in sorted(
        writes,
        key=lambda x: x["line"]
    ):

        print(
            "[SCRITTURA ISOLABILE]",
            row["line"],
            row["destination"]
        )

    print(
        "[EMBEDDING DA PROTEGGERE]"
    )

    for path in EMBEDDINGS:

        print(
            path.relative_to(ROOT),
            "| presente:",
            path.is_file()
        )

    print(
        "[OK] Analisi soltanto: nessun training"
    )


# ============================================================
# 13. INSTALLAZIONE ISOLAMENTO GAE
# ============================================================

def install():

    config(False)

    state = current()

    if (
        not state.get("running")
        or state.get("image") != IMAGE
    ):

        raise RuntimeError(
            "Container V3.2 non conforme"
        )

    before = originals()

    patched, writes = patch_source()

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP22_6O_GAE_GUI_SAFE_" + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    # Prima di qualsiasi intervento, proteggiamo
    # i due embedding non inclusi nello STEP 22.6N.

    saved = backup_embeddings(
        folder / "backup_embedd_originali"
    )

    # Il trainer originale non viene modificato.
    # Creiamo una copia con percorsi output isolati.

    write_once(
        PATCH,
        patched
    )

    override = {
        "services": {
            SERVICE: {
                "volumes": [
                    {
                        "type": "bind",
                        "source": str(PATCH),
                        "target": DEST,
                        "read_only": True
                    }
                ]
            }
        }
    }

    write_once(
        OVERRIDE,
        json.dumps(
            override,
            indent=2,
            ensure_ascii=False
        ) + "\n"
    )

    config(True)

    if originals() != before:

        raise RuntimeError(
            "Hash originali cambiati: "
            "non attivo l'override"
        )

    print(
        "[ATTIVAZIONE] Montaggio trainer protetto "
        "senza eseguire training"
    )

    try:

        start(True)

        active = current()

        if (
            Path(
                active.get("trainer_mount") or ""
            ).resolve()
            != PATCH.resolve()
        ):

            raise RuntimeError(
                "Il trainer attivo non è "
                "quello isolato"
            )

        # Verifica il checksum del file montato
        # leggendo il sorgente senza eseguirlo.

        code = (
            "import hashlib,pathlib;"
            "print(hashlib.sha256("
            "pathlib.Path('"
            + DEST
            + "').read_bytes()).hexdigest())"
        )

        result = cmd([
            "docker",
            "exec",
            CONTAINER,
            "python3",
            "-B",
            "-c",
            code
        ])

        if result.stdout.strip() != sha(PATCH):

            raise RuntimeError(
                "Checksum della patch differente "
                "nel container"
            )

        if originals() != before:

            raise RuntimeError(
                "Un file originale è cambiato"
            )

    except Exception:

        print(
            "[RIPRISTINO] Reimposto "
            "il trainer originale"
        )

        try:
            start(False)

        except Exception as exc:

            print(
                "[ERRORE RIPRISTINO]",
                exc
            )

        raise

    report = {
        "step": "22.6O",
        "status": "TRAINER_ISOLATO_SENZA_TRAINING",

        "trainer_original_sha256": before[
            str(TRAINER.relative_to(ROOT))
        ],

        "trainer_patch_sha256": sha(PATCH),
        "writes": writes,
        "backup_embedd": saved,

        "safe_outputs": (
            "outputs/gae_gui_safe_runs/"
            "<run_id>/{results,models,embeddings}"
        ),

        "training_performed": False,
        "gui_modified": False
    }

    report_path = (
        folder / "RAPPORTO_STEP22_6O.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print(
        "[OK] Trainer originale preservato; "
        "output futuri isolati"
    )

    print(
        "[OK] Nessun training eseguito"
    )

    print(
        "[REPORT]",
        report_path
    )

    status()


# ============================================================
# 14. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    group = parser.add_mutually_exclusive_group(
        required=True
    )

    for mode in (
        "audit",
        "install",
        "status",
        "rollback"
    ):

        group.add_argument(
            "--" + mode,
            action="store_true"
        )

    args = parser.parse_args()

    print("=" * 68)

    print(
        "AI PRIMETTA V3.2 - STEP 22.6O "
        "| ISOLAMENTO TRAINING GAE"
    )

    print("=" * 68)

    if args.audit:

        audit()

    elif args.install:

        install()

    elif args.status:

        status()

    else:

        start(False)

        print(
            "[OK] Ripristinato il trainer "
            "originale; backup mantenuti"
        )

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