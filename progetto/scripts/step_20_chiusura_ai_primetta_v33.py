#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 20
Chiusura tecnico-scientifica documentata.

Legge soltanto evidenze e report precedenti.
Scrive un nuovo rapporto nella cartella runs/.

Non modifica GUI, motori, checkpoint o dataset.
Non esegue training.
"""

import csv
import hashlib
import json
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

RUNS = ROOT / "target_prediction_v33_vscode" / "runs"

EVIDENCE = RUNS / "STEP20_EVIDENZE_GUI"

FILES = (
    "01_rna_qc.tsv",
    "02_alignment_hypotheses.tsv",
    "03_historical_12mer_hypotheses.tsv",
    "04_gene_target_candidates.tsv",
    "05_rna_candidate_summary.tsv",
    "06_provenance_manifest.json",
)


def sha256(path):

    digest = hashlib.sha256()

    with path.open("rb") as stream:
        for block in iter(
            lambda: stream.read(1024 * 1024),
            b""
        ):
            digest.update(block)

    return digest.hexdigest()


def read_tsv(path):

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        if not reader.fieldnames:
            # Un TSV vuoto puo essere un risultato
            # legittimo quando non esistono candidati.
            return {
                "columns": [],
                "row_count": 0,
                "header_absent": True
            }, []

        columns = list(reader.fieldnames)
        rows = list(reader)

        if any(None in row for row in rows):
            raise ValueError(
                "Riga con piu campi dell'intestazione"
            )

        return {
            "columns": columns,
            "row_count": len(rows)
        }, rows


def previous(prefix, filename):

    paths = sorted(
        RUNS.glob(f"{prefix}_*/{filename}")
    )

    if not paths:
        return {
            "found": False,
            "status": "MISSING"
        }

    path = paths[-1]

    try:

        data = json.loads(
            path.read_text(encoding="utf-8")
        )

        return {
            "found": True,
            "path": str(path),
            "status": data.get("status"),
            "data": data
        }

    except (ValueError, OSError) as exc:

        return {
            "found": True,
            "path": str(path),
            "status": "INVALID",
            "error": str(exc)
        }


def http_health():

    url = "http://127.0.0.1:8502/_stcore/health"

    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({})
    )

    try:

        with opener.open(
            url,
            timeout=8
        ) as response:

            content = response.read(100).decode(
                "utf-8",
                errors="replace"
            ).strip()

            return {
                "ok": (
                    response.status == 200
                    and content == "ok"
                ),
                "http_code": response.status,
                "body": content
            }

    except Exception as exc:

        return {
            "ok": False,
            "error": str(exc)
        }


def main():

    print("=" * 72)
    print("AI PRIMETTA V3.3 — STEP 20")
    print("DOSSIER TECNICO-SCIENTIFICO FINALE")
    print("=" * 72)

    if not EVIDENCE.is_dir():
        raise SystemExit(
            "[ERRORE] Cartella evidenze non trovata: "
            + str(EVIDENCE)
        )

    inventory = {}
    tables = {}

    # ==================================================
    # 1. INVENTARIO EVIDENZE
    # ==================================================

    print("\n[1] EVIDENZE ESPORTATE DALLA GUI")

    for name in FILES:

        path = EVIDENCE / name

        if not path.is_file():

            inventory[name] = {
                "present": False
            }

            print("[MANCA]", name)
            continue

        entry = {
            "present": True,
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path)
        }

        try:

            if name.endswith(".tsv"):

                metadata, rows = read_tsv(path)

                entry.update(metadata)
                tables[name] = rows

            else:

                manifest = json.loads(
                    path.read_text(
                        encoding="utf-8-sig"
                    )
                )

                entry["valid_json"] = True

                entry["top_level_keys"] = (
                    sorted(manifest)
                    if isinstance(manifest, dict)
                    else []
                )

        except (
            ValueError,
            OSError,
            UnicodeError,
            csv.Error
        ) as exc:

            entry["parse_error"] = str(exc)

        inventory[name] = entry

        print(
            "[OK]" if "parse_error" not in entry
            else "[ERRORE]",
            name,
            "| righe:",
            entry.get("row_count", "n/d")
        )

    # ==================================================
    # 2. CONTROLLO SCIENTIFICO DELLE CINQUE sRNA
    # ==================================================

    print("\n[2] COERENZA DEI RISULTATI")

    qc = tables.get("01_rna_qc.tsv", [])

    screening = tables.get(
        "05_rna_candidate_summary.tsv",
        []
    )

    qc_ids = [
        str(row.get("rna_id", ""))
        for row in qc
    ]

    summary_ids = [
        str(row.get("rna_id", ""))
        for row in screening
    ]

    qc_eligible = sum(
        str(row.get("status", "")).strip().lower()
        == "eligible"
        for row in qc
    )

    qc_entropy_ok = False

    try:

        qc_entropy_ok = bool(qc) and all(
            float(row["entropy_bits"]) >= 1.40
            for row in qc
        )

    except (ValueError, KeyError, TypeError):
        pass

    sites = None

    try:

        if screening:
            sites = sum(
                int(row["n_new_sites"])
                for row in screening
            )

    except (ValueError, KeyError, TypeError):
        pass

    data_check = {
        "qc_rows": len(qc),
        "qc_eligible": qc_eligible,
        "qc_entropy_threshold_1_40_satisfied":
            qc_entropy_ok,
        "screening_rows": len(screening),
        "same_rna_ids": (
            len(qc_ids) == len(summary_ids) == 5
            and set(qc_ids) == set(summary_ids)
        ),
        "total_new_sites_screening": sites,
        "alignment_rows": inventory.get(
            "02_alignment_hypotheses.tsv", {}
        ).get("row_count"),
        "gene_candidates_rows": inventory.get(
            "04_gene_target_candidates.tsv", {}
        ).get("row_count"),
        "historical_associations_activated": False
    }

    for key, value in data_check.items():
        print(f"{key}: {value}")

    # ==================================================
    # 3. REPORT MODELLI ML/GNN
    # ==================================================

    print("\n[3] REPORT COLLAUDI E SERVIZIO GUI")

    gae = previous(
        "STEP19_3_INFERENZA_GAE_V3",
        "RAPPORTO_STEP19_3.json"
    )

    ml = previous(
        "STEP19_5_COMPATIBILITA_ML",
        "RAPPORTO_STEP19_5.json"
    )

    gui = previous(
        "STEP19_7_GUI_LIVE",
        "RAPPORTO_STEP19_7.json"
    )

    health = http_health()

    gae_ok = (
        gae["status"] == "TECHNICAL_INFERENCE_PASS"
    )

    ml_data = (
        ml.get("data", {})
        .get("ml_inference", {})
    )

    ml_ok = (
        ml["status"] == "COMPATIBLE_ML_INFERENCE_PASS"
        and ml_data.get("rows_tested") == 500
        and ml_data.get(
            "historical_label_agreement"
        ) == 1.0
    )

    print(
        "GAE V3:",
        "PASS" if gae_ok else "NON CONFERMATO"
    )

    print(
        "ML antileakage:",
        "PASS" if ml_ok else "NON CONFERMATO"
    )

    print(
        "GUI HTTP:",
        "PASS" if health["ok"]
        else "NON RAGGIUNGIBILE"
    )

    print(
        "Ultimo report GUI live:",
        gui["status"]
    )

    # ==================================================
    # 4. VALUTAZIONE COMPLESSIVA
    # ==================================================

    print("\n[4] VALUTAZIONE COMPLESSIVA")

    complete = all(
        inventory[name].get("present")
        and not inventory[name].get("parse_error")
        for name in FILES
    )

    scientific_data_consistent = (
        len(qc) == 5
        and qc_eligible == 5
        and qc_entropy_ok
        and len(screening) == 5
        and data_check["same_rna_ids"]
        and sites == 0
        and data_check["alignment_rows"] == 0
        and data_check["gene_candidates_rows"] == 0
    )

    limitations = [
        (
            "Le predizioni ML storiche non equivalgono "
            "a una validazione indipendente."
        ),
        (
            "La GAE V3 e stata eseguita sul grafo completo; "
            "non e stato riprodotto lo split di valutazione "
            "originale."
        ),
        (
            "Il collegamento diretto della GUI di predizione "
            "ai checkpoint GAE/ML non e stato verificato."
        ),
        (
            "Nessun target RNA-trascritto e stato identificato "
            "con i filtri del collaudo; questo non esclude "
            "regolazione biologica."
        ),
        (
            "Le associazioni storiche non sono state attivate "
            "in questa esecuzione."
        ),
        (
            "Non e stata effettuata validazione biologica "
            "sperimentale."
        )
    ]

    if not complete:

        status = "DOSSIER_EVIDENZE_INCOMPLETO"

    elif not scientific_data_consistent:

        status = "DOSSIER_DATI_DA_VERIFICARE"

    elif not (
        gae_ok
        and ml_ok
        and health["ok"]
    ):

        status = "DOSSIER_COLLAUDI_DA_VERIFICARE"

    else:

        status = (
            "CHIUSURA_TECNICA_DOCUMENTATA_"
            "CON_LIMITI_SCIENTIFICI"
        )

    # ==================================================
    # 5. CREAZIONE CARTELLA REPORT FINALE
    # ==================================================

    print("\n[5] GENERAZIONE RAPPORTO FINALE")

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    out = RUNS / (
        f"STEP20_RAPPORTO_FINALE_{stamp}"
    )

    out.mkdir(
        parents=True,
        exist_ok=False
    )

    report = {
        "project": "AI PRIMETTA V3.3",
        "step": "20",
        "status": status,
        "evidence_directory": str(EVIDENCE),
        "evidence_inventory": inventory,
        "data_checks": data_check,
        "gae_report": {
            k: v for k, v in gae.items()
            if k != "data"
        },
        "ml_report": {
            k: v for k, v in ml.items()
            if k != "data"
        },
        "gui_live_report": {
            k: v for k, v in gui.items()
            if k != "data"
        },
        "gui_health": health,
        "limitations": limitations,
        "independent_scientific_validation": False,
        "gui_to_gae_ml_link_tested": False,
        "training_performed": False,
        "existing_components_modified": False
    }

    json_path = out / "RAPPORTO_FINALE_STEP20.json"

    json_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ),
        encoding="utf-8"
    )

    # ==================================================
    # 6. RAPPORTO MARKDOWN LEGGIBILE IN VS CODE
    # ==================================================

    md = [
        "# AI PRIMETTA V3.3 — Rapporto finale STEP 20",
        "",
        f"**Stato:** `{status}`",
        "",
        "## Evidenze GUI",
        "",
        "| File | Presente | Righe | SHA-256 |",
        "|---|---|---:|---|"
    ]

    for name, item in inventory.items():

        present = (
            "Si" if item.get("present")
            else "NO"
        )

        md.append(
            f"| `{name}` | {present} | "
            f"{item.get('row_count', '—')} | "
            f"`{item.get('sha256', '—')}` |"
        )

    md += [
        "",
        "## Verifiche",
        "",
        (
            f"- RNA sottoposte a QC: **{len(qc)}**; "
            f"eleggibili: **{qc_eligible}**."
        ),
        (
            f"- RNA nello screening: **{len(screening)}**; "
            f"nuovi siti: **{sites if sites is not None else 'non determinati'}**."
        ),
        (
            "- Siti RNA-trascritto: "
            f"**{data_check['alignment_rows']}**."
        ),
        (
            "- Geni candidati: "
            f"**{data_check['gene_candidates_rows']}**."
        ),
        (
            "- GAE V3: "
            f"**{'PASS' if gae_ok else 'NON CONFERMATO'}**."
        ),
        (
            "- ML antileakage (500 osservazioni): "
            f"**{'PASS' if ml_ok else 'NON CONFERMATO'}**."
        ),
        (
            "- GUI HTTP 8502: "
            f"**{'PASS' if health['ok'] else 'NON RAGGIUNGIBILE'}**."
        ),
        "",
        "## Limiti scientifici",
        ""
    ]

    md.extend(
        f"- {limitation}"
        for limitation in limitations
    )

    md += [
        "",
        (
            "Nessun training o modifica ai modelli "
            "effettuati da questo script."
        ),
        ""
    ]

    markdown_path = (
        out / "RAPPORTO_FINALE_STEP20.md"
    )

    markdown_path.write_text(
        "\n".join(md),
        encoding="utf-8"
    )

    # ==================================================
    # 7. RISULTATO
    # ==================================================

    print()
    print("=" * 72)
    print("STEP 20 - RISULTATO")
    print("=" * 72)

    print("Stato:", status)
    print("Report JSON:", json_path)
    print("Report Markdown:", markdown_path)

    print("[OK] File originali non modificati")
    print("[OK] Nessun training eseguito")

    if status != (
        "CHIUSURA_TECNICA_DOCUMENTATA_"
        "CON_LIMITI_SCIENTIFICI"
    ):
        print(
            "[ATTENZIONE] Dossier non ancora completo. "
            "Consultare il rapporto."
        )


if __name__ == "__main__":
    main()