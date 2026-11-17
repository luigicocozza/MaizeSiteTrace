#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AI PRIMETTA V3.2 - STEP 22.6S.

Collaudo finale in sola lettura.
Non modifica GUI, Docker, modelli o dati.
Produce soltanto report JSON e Markdown.
"""

import csv
import hashlib
import json
import subprocess
import sys
import urllib.request

from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

CONTAINER = "maize_ncrna_gnn"

IMAGE = "ai_primetta_v32_full_ai_sklearn18:22.6g"

RUNS = ROOT / "target_prediction_v33_vscode/runs"

SAFE_RUNS = ROOT / "outputs/gae_gui_safe_runs"

REQUIRED_MOUNTS = {
    "/app/25_train_graph_autoencoder.py":
        "25_train_graph_autoencoder_safe.py",

    "/app/app/pages/01_Rete_neurale_GAE.py":
        "01_Rete_neurale_GAE_latest.py",

    "/app/api/main.py":
        "api_main_rag_v32.py",

    "/app/data/user_data": None,
}

REPORT = {
    "checks": {},
    "evidence": {},
    "remaining_manual_checks": [
        "Verificare la rigenerazione della knowledge base "
        "utente dopo l'ultimo caricamento.",

        "Confrontare i checksum dei download GAE "
        "e ranking effettuati dal browser."
    ]
}


def check(name, function):

    try:

        value = function()

        REPORT["checks"][name] = {
            "passed": True,
            "result": value
        }

        print("[PASS]", name)

    except Exception as exc:

        REPORT["checks"][name] = {
            "passed": False,
            "error": str(exc)
        }

        print("[FAIL]", name, "-", exc)


def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as file:

        for part in iter(
            lambda: file.read(1024 * 1024),
            b""
        ):

            h.update(part)

    return h.hexdigest()


# ============================================================
# VERIFICA CONTAINER E MONTAGGI
# ============================================================

def inspect_docker():

    proc = subprocess.run(
        ["docker", "inspect", CONTAINER],
        capture_output=True,
        text=True,
        timeout=20
    )

    if proc.returncode:

        raise RuntimeError(proc.stderr[-500:])

    obj = json.loads(proc.stdout)[0]

    assert obj["State"]["Running"]

    assert obj["Config"]["Image"] == IMAGE

    mounts = {
        m["Destination"]: m
        for m in obj["Mounts"]
    }

    for target, suffix in REQUIRED_MOUNTS.items():

        assert target in mounts, (
            "Montaggio assente: " + target
        )

        item = mounts[target]

        source = Path(item["Source"])

        if suffix:

            assert source.name == suffix

            assert not item["RW"]

        else:

            assert source.resolve() == (
                ROOT / "data/user_data"
            ).resolve()

            assert item["RW"]

    assert "/app/outputs" in mounts

    return {
        "image": IMAGE,
        "mounts_verified": list(REQUIRED_MOUNTS)
    }


# ============================================================
# VERIFICA HTTP
# ============================================================

def http_status(port, path):

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    url = (
        f"http://127.0.0.1:{port}{path}"
    )

    with opener.open(
        url,
        timeout=8
    ) as response:

        assert response.status == 200

    return "HTTP 200"


# ============================================================
# VERIFICA ULTIMO DATASET E RANKING
# ============================================================

def newest_dataset():

    folder = ROOT / "data/user_data/uploads"

    files = [
        p for p in folder.glob(
            "*/user_differential_expression.csv"
        )
        if p.is_file()
    ]

    assert files, "Nessun CSV utente trovato"

    latest = max(
        files,
        key=lambda p: (
            p.parent.stat().st_mtime_ns,
            p.parent.name
        )
    )

    integrated = (
        ROOT
        / "outputs/user_data"
        / latest.parent.name
        / "ranking_with_user_evidence.csv"
    )

    assert integrated.is_file(), (
        "Ranking integrato non trovato"
    )

    source = {}

    with latest.open(
        encoding="utf-8-sig",
        newline=""
    ) as file:

        for row in csv.DictReader(file):

            identifier = (
                row["rnacentral_id"]
                .strip()
                .upper()
            )

            assert identifier not in source

            source[identifier] = Decimal(
                row["FDR"].replace(",", ".")
            )

    joined = {}

    with integrated.open(
        encoding="utf-8-sig",
        newline=""
    ) as file:

        for row in csv.DictReader(file):

            identifier = (
                row["rnacentral_id"]
                .strip()
                .upper()
            )

            if identifier not in source:
                continue

            assert identifier not in joined

            fdr = Decimal(
                row["user_min_FDR"].replace(",", ".")
            )

            assert fdr == source[identifier]

            assert str(
                row["user_data_present"]
            ).strip().lower() in (
                "true", "1", "yes", "si"
            )

            joined[identifier] = fdr

    assert len(joined) == len(source)

    selected = {
        identifier: str(fdr)
        for identifier, fdr in joined.items()
        if fdr < Decimal("0.05")
    }

    REPORT["evidence"]["upload"] = str(
        latest.relative_to(ROOT)
    )

    REPORT["evidence"]["ranking"] = str(
        integrated.relative_to(ROOT)
    )

    REPORT["evidence"]["fdr_matches"] = selected

    return {
        "dataset": latest.parent.name,
        "input_records": len(source),
        "matched": len(joined),
        "fdr_below_005": len(selected)
    }


# ============================================================
# VERIFICA API RAG SUL DATASET CORRENTE
# ============================================================

def rag_api():

    payload = {
        "question": (
            "Quali ncRNA caricati dal ricercatore "
            "sono presenti nel ranking e hanno "
            "FDR inferiore a 0.05?"
        ),
        "top_k": 50
    }

    request = urllib.request.Request(
        "http://127.0.0.1:8000/rag/ask",
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json"
        },
        method="POST"
    )

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    with opener.open(
        request,
        timeout=60
    ) as response:

        result = json.load(response)

    assert result.get("provider") == (
        "deterministic_user_fdr_v32"
    )

    expected = REPORT["evidence"]["fdr_matches"]

    observed = {
        str(row["rnacentral_id"]).upper():
            Decimal(str(row["user_min_FDR"]))

        for row in result.get("evidence", [])
    }

    assert result.get("total_matches") == len(expected)

    assert observed == {
        key: Decimal(value)
        for key, value in expected.items()
    }

    sources = json.dumps(
        result.get("sources", []),
        ensure_ascii=False
    )

    dataset = Path(
        REPORT["evidence"]["upload"]
    ).parent.name

    assert dataset in sources, (
        "L'API non indica l'ultimo dataset"
    )

    return {
        "provider": result["provider"],
        "matches": len(observed),
        "dataset": dataset
    }


# ============================================================
# VERIFICA ULTIMO TRAINING GAE
# ============================================================

def gae_run():

    directories = sorted(
        (
            p for p in SAFE_RUNS.iterdir()
            if p.is_dir()
            and not p.is_symlink()
        ),
        key=lambda p: p.name
    )

    assert directories, "Nessun run GAE isolato"

    latest = directories[-1]

    required = [
        "models/gae_model.pt",
        "embeddings/gae_node_embeddings.csv",
        "embeddings/gae_ncrna_embeddings.csv",
        "results/gae_metrics.json",
        "results/gae_training_history.csv",
        "results/plots/gae_loss_curve.png",
        "results/plots/gae_validation_ap_curve.png"
    ]

    for relative in required:

        path = latest / relative

        assert path.is_file(), relative

        assert path.stat().st_size > 0, relative

    metrics = json.loads(
        (
            latest / "results/gae_metrics.json"
        ).read_text(encoding="utf-8")
    )

    REPORT["evidence"]["gae_run"] = str(
        latest.relative_to(ROOT)
    )

    REPORT["evidence"]["gae_model_sha256"] = sha256(
        latest / "models/gae_model.pt"
    )

    REPORT["evidence"]["gae_metrics_sha256"] = sha256(
        latest / "results/gae_metrics.json"
    )

    return {
        "run": latest.name,
        "required_outputs": len(required),
        "metrics_sections": list(metrics)
    }


# ============================================================
# VERIFICA FILE SCIENTIFICI ORIGINALI
# ============================================================

def protect_previous_files():

    files = sorted(
        RUNS.glob(
            "STEP22_6P_GAE_BUTTON_*/"
            "BASELINE_STEP22_6P.json"
        )
    )

    assert files, (
        "Baseline integrità STEP 22.6P mancante"
    )

    baseline = json.loads(
        files[-1].read_text(encoding="utf-8")
    )

    originals = baseline["original_sha256"]

    changed = []

    for relative, metadata in originals.items():

        path = ROOT / relative

        if (
            not path.is_file()
            or sha256(path) != metadata["sha256"]
        ):

            changed.append(relative)

    assert not changed, (
        "File originali modificati: "
        + str(changed)
    )

    return {
        "protected_files_unchanged": len(originals)
    }


# ============================================================
# RAPPORTO FINALE
# ============================================================

def main():

    print("=" * 68)
    print("AI PRIMETTA V3.2 | STEP 22.6S")
    print("CONSOLIDAMENTO E COLLAUDO FINALE")
    print("=" * 68)

    check("Docker e patch", inspect_docker)

    services = [
        ("GUI V3.2", 8503, "/_stcore/health"),
        ("API V3.2", 8000, "/health"),
        ("GUI V3.3", 8502, "/_stcore/health"),
        ("AI ORZO", 8501, "/_stcore/health")
    ]

    for name, port, path in services:

        check(
            name,
            lambda p=port, route=path:
                http_status(p, route)
        )

    check(
        "Ultimo CSV e ranking",
        newest_dataset
    )

    if REPORT["checks"][
        "Ultimo CSV e ranking"
    ]["passed"]:

        check(
            "Risposta FDR API sul nuovo dataset",
            rag_api
        )

    check(
        "Ultimo GAE isolato",
        gae_run
    )

    check(
        "File scientifici originali",
        protect_previous_files
    )

    passed = all(
        item["passed"]
        for item in REPORT["checks"].values()
    )

    REPORT[
        "all_automatic_checks_passed"
    ] = passed

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    folder = RUNS / (
        "STEP22_6S_RELEASE_V32_" + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    json_path = (
        folder / "RAPPORTO_FINAL_V32.json"
    )

    json_path.write_text(
        json.dumps(
            REPORT,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    lines = [
        "# AI PRIMETTA V3.2 — Rapporto finale",
        "",
        "Data UTC: " + stamp,
        "",
        "## Controlli automatici",
        ""
    ]

    for name, result in REPORT["checks"].items():

        status = (
            "PASS"
            if result["passed"]
            else "FAIL"
        )

        detail = result.get(
            "result",
            result.get("error")
        )

        lines.append(
            f"- {status} — {name}: {detail}"
        )

    lines.extend([
        "",
        "## Verifiche manuali residue",
        ""
    ])

    for item in REPORT["remaining_manual_checks"]:

        lines.append("- " + item)

    lines.extend([
        "",
        "La conformità funzionale non equivale "
        "a validazione biologica indipendente."
    ])

    md_path = (
        folder / "RAPPORTO_FINAL_V32.md"
    )

    md_path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8"
    )

    print()
    print("[REPORT JSON]", json_path)
    print("[REPORT MARKDOWN]", md_path)

    print(
        "[ESITO]",
        "CONTROLLI AUTOMATICI PASS"
        if passed
        else "DA RISOLVERE"
    )

    print(
        "[NOTA] Rimangono due verifiche "
        "manuali indicate nel report"
    )

    sys.exit(0 if passed else 2)


if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__,
            exc,
            file=sys.stderr
        )

        sys.exit(1)