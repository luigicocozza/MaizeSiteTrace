#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 23 - psRNATarget V2

PRIMA ANALISI:
- sRNA_00000435: 24 nt
- sRNA_00000037: 25 nt
- 37 trascritti reali B73 NAM-5.0

FUNZIONALITA':
    --audit
        Verifica in sola lettura i file STEP22.

    --run
        Crea una sessione con FASTA verificati
        e apre il server psRNATarget.

    --import-results
        Importa TSV/CSV esportato dal server
        e confronta le predizioni con STEP22.

IMPORTANTE:
La predizione viene eseguita dal server
ufficiale psRNATarget, NON da questo script.

Non modifica:
- GUI V3.2
- GUI V3.3
- Docker
- Modelli ML/GNN
- Dataset scientifici originali
"""

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
import webbrowser

from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE GENERALE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

WEB_URL = (
    "https://www.zhaolab.org/"
    "psRNATarget/analysis?function=3"
)

STANDARD = {
    "sRNA_00000435": 24,
    "sRNA_00000037": 25
}

OUTPUT_PREFIX = (
    "STEP23_PSRNATARGET_STANDARD_"
)


# ============================================================
# 2. TIMESTAMP
# ============================================================

def stamp():

    return datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )


# ============================================================
# 3. CHECKSUM SHA256
# ============================================================

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
# 4. PERCORSI RELATIVI
# ============================================================

def rel(path):

    try:

        return str(
            path.relative_to(ROOT)
        )

    except ValueError:

        return str(path)


# ============================================================
# 5. LETTURA FILE FASTA
# ============================================================

def fasta(path):

    sequences = {}

    name = None
    chunks = []

    with path.open(
        "r",
        encoding="utf-8-sig"
    ) as file:

        for raw in file:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    if not chunks:

                        raise ValueError(
                            f"Sequenza vuota: {name}"
                        )

                    sequences[name] = (
                        "".join(chunks).upper()
                    )

                name = line[1:].split()[0]

                if (
                    not name
                    or name in sequences
                ):

                    raise ValueError(
                        "ID FASTA duplicato o vuoto"
                    )

                chunks = []

            elif name is None:

                raise ValueError(
                    "FASTA senza intestazione"
                )

            else:

                chunks.append(line)

        if name is not None:

            if (
                not chunks
                or name in sequences
            ):

                raise ValueError(
                    "Record FASTA finale invalido"
                )

            sequences[name] = (
                "".join(chunks).upper()
            )

    if not sequences:

        raise ValueError(
            f"FASTA vuoto: {path}"
        )

    return sequences


# ============================================================
# 6. LETTURA TSV
# ============================================================

def tsv(path):

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as file:

        reader = csv.DictReader(
            file,
            delimiter="\t"
        )

        headers = reader.fieldnames or []

        if (
            not headers
            or len(set(headers)) != len(headers)
        ):

            raise ValueError(
                f"Colonne TSV non valide: {path}"
            )

        return headers, list(reader)


# ============================================================
# 7. VERIFICA FILE SCIENTIFICI STEP22
# ============================================================

def check_step22():

    paths = sorted(
        RUNS.glob(
            "STEP22_EXPORT_TARGET_B73_*/"
            "MANIFEST_STEP22.json"
        )
    )

    if not paths:

        raise FileNotFoundError(
            "Manca MANIFEST_STEP22.json"
        )

    manifest_path = paths[-1]

    manifest = json.loads(
        manifest_path.read_text(
            encoding="utf-8"
        )
    )

    if (
        manifest.get("step") != "22"
        or manifest.get("status")
        != "PREPARED_NOT_SUBMITTED"
    ):

        raise ValueError(
            "Manifest STEP22 inatteso"
        )

    root = manifest_path.parent

    files = {

        "standard":
            root / "SRNA_STANDARD_MAX25.fa",

        "targets":
            root / "TRASCRITTI_B73_CANDIDATI.fa",

        "sites":
            root / "SITI_E_RNASEQ_SNAPSHOT.tsv"
    }

    # --------------------------------------------------------
    # CONTROLLO CHECKSUM
    # --------------------------------------------------------

    for kind, path in files.items():

        if (
            not path.is_file()
            or sha256(path)
            != manifest["files_sha256"][kind]
        ):

            raise ValueError(
                "File mancante o checksum diverso: "
                + str(path)
            )

    # --------------------------------------------------------
    # SMALL RNA
    # --------------------------------------------------------

    srnas = fasta(
        files["standard"]
    )

    targets = fasta(
        files["targets"]
    )

    for ident, length in STANDARD.items():

        if (
            ident not in srnas
            or len(srnas[ident]) != length
        ):

            raise ValueError(
                f"sRNA non corretta: {ident}"
            )

    if set(srnas) != set(STANDARD):

        raise ValueError(
            "Il FASTA contiene small RNA non attese"
        )

    if any(
        set(seq) - set("ACGTU")
        for seq in srnas.values()
    ):

        raise ValueError(
            "sRNA con simboli non canonici"
        )

    # --------------------------------------------------------
    # TRASCRITTI B73
    # --------------------------------------------------------

    if (
        len(targets)
        != manifest["unique_transcripts"]
        or len(targets) != 37
    ):

        raise ValueError(
            "Numero trascritti non coerente"
        )

    if any(
        not (50 <= len(seq) <= 5000000)
        or set(seq) - set("ACGTUN")
        for seq in targets.values()
    ):

        raise ValueError(
            "Trascritto non valido per psRNATarget"
        )

    # --------------------------------------------------------
    # SITI PRELIMINARI
    # --------------------------------------------------------

    _, sites = tsv(
        files["sites"]
    )

    if (
        len(sites)
        != manifest["candidate_sites"]
        or len(sites) != 37
    ):

        raise ValueError(
            "Numero siti non coerente"
        )

    counts = Counter()

    by_tx = defaultdict(list)

    for site in sites:

        ident = site["srna_id"]

        tx = site["transcript_id"]

        if tx not in targets:

            raise ValueError(
                f"Trascritto assente: {tx}"
            )

        start = int(
            site["start_1based"]
        )

        end = int(
            site["end_1based"]
        )

        fragment = (
            targets[tx][start - 1:end]
            .replace("U", "T")
        )

        expected = (
            site["full_site_DNA_5to3"]
            .upper()
            .replace("U", "T")
        )

        if (
            start < 1
            or fragment != expected
        ):

            raise ValueError(
                "Coordinate STEP22 non coerenti: "
                f"{ident}, {tx}"
            )

        counts[ident] += 1

        by_tx[
            (ident, tx)
        ].append(site)

    if (
        counts["sRNA_00000435"] != 19
        or counts["sRNA_00000037"] != 17
    ):

        raise ValueError(
            "Numero siti delle due sRNA inatteso"
        )

    print(
        "[PASS] STEP22 / SHA256: "
        "file originali invariati"
    )

    print(
        "[PASS] sRNA:",
        ", ".join(
            f"{x}: {len(y)} nt"
            for x, y in sorted(
                srnas.items()
            )
        )
    )

    print(
        f"[PASS] Trascritti: {len(targets)}; "
        f"siti STEP22: {len(sites)}; "
        "siti delle due sRNA: 36"
    )

    print(
        "[NOTA] psRNATarget analizzerà "
        "le due sRNA contro TUTTI i 37 trascritti."
    )

    return (
        manifest_path,
        files,
        srnas,
        targets,
        sites,
        by_tx
    )


# ============================================================
# 8. PREPARAZIONE PRIMA ANALISI
# ============================================================

def prepare(open_browser):

    (
        manifest_path,
        files,
        srnas,
        targets,
        sites,
        by_tx
    ) = check_step22()

    folder = RUNS / (
        OUTPUT_PREFIX + stamp()
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    names = {

        "standard":
            "01_SMALL_RNA_STANDARD_24_25nt.fa",

        "targets":
            "02_TRASCRITTI_B73_37.fa",

        "sites":
            "03_SITI_STEP22_RIFERIMENTO.tsv"
    }

    exported = {}

    # --------------------------------------------------------
    # COPIA DEI FILE GIÀ VALIDATI
    # --------------------------------------------------------

    for key, filename in names.items():

        destination = (
            folder / filename
        )

        shutil.copyfile(
            files[key],
            destination
        )

        exported[key] = destination

    # --------------------------------------------------------
    # MANIFEST SCIENTIFICO
    # --------------------------------------------------------

    protocol = {

        "study":
            "AI PRIMETTA V3.3 - "
            "psRNATarget prima analisi",

        "status":
            "READY_FOR_MANUAL_WEB_SUBMISSION",

        "source_step22":
            rel(manifest_path),

        "url":
            WEB_URL,

        "scoring_schema":
            "V2 (2017 release)",

        "parameters":
            "Valori predefiniti del server; "
            "registrare quelli effettivamente mostrati",

        "accessibility":
            "Registrare se abilitata",

        "files": {

            key: {

                "path": str(path),

                "sha256": sha256(path)
            }

            for key, path in exported.items()
        },

        "srna_ids":
            sorted(srnas),

        "target_count":
            len(targets),

        "previous_step22_sites":
            len(sites),

        "prediction_performed":
            False,

        "external_job_url":
            None,

        "limitations": [

            "Target preselezionati "
            "mediante screening locale",

            "Non implica verifica biologica",

            "Non costituisce predizione "
            "genome-wide"
        ]
    }

    manifest_file = (
        folder / "MANIFEST_STEP23.json"
    )

    manifest_file.write_text(
        json.dumps(
            protocol,
            ensure_ascii=False,
            indent=2
        ) + "\n",
        encoding="utf-8"
    )

    # --------------------------------------------------------
    # GUIDA PER L'INVIO
    # --------------------------------------------------------

    guide = (

        "# AI PRIMETTA V3.3 - psRNATarget V2\n\n"

        f"Apri: {WEB_URL}\n\n"

        "1. Seleziona Submit small RNAs and targets.\n"

        "2. Upload sRNA file:\n"
        f"{exported['standard']}\n\n"

        "3. Upload target file:\n"
        f"{exported['targets']}\n\n"

        "4. Scoring Schema: V2 (2017 release).\n"

        "5. Mantieni inizialmente gli altri "
        "parametri predefiniti.\n"

        "6. Registra i parametri e "
        "l'opzione accessibility.\n"

        "7. Invia il lavoro dal browser.\n"

        "8. Conserva l'URL del lavoro.\n"

        "9. Scarica i risultati TSV/CSV.\n"

        "10. Importa i risultati con:\n\n"

        "python3 "
        "scripts/step_23_psrnatarget_standard_v33.py "
        "--import-results /percorso/risultati.tsv\n\n"

        "ATTENZIONE: nessuna predizione "
        "viene eseguita localmente.\n"
    )

    guide_file = (
        folder / "GUIDA_INVIO.md"
    )

    guide_file.write_text(
        guide,
        encoding="utf-8"
    )

    # --------------------------------------------------------
    # OUTPUT TERMINALE
    # --------------------------------------------------------

    print(
        "[RUN]",
        folder
    )

    print(
        "[CARICA small RNA]",
        exported["standard"]
    )

    print(
        "[CARICA target]",
        exported["targets"]
    )

    print(
        "[SCHEMA] V2 2017; "
        "parametri predefiniti da registrare"
    )

    # --------------------------------------------------------
    # APERTURA BROWSER
    # --------------------------------------------------------

    if open_browser:

        try:

            ok = webbrowser.open_new_tab(
                WEB_URL
            )

        except Exception:

            ok = False

        if ok:

            print(
                "[BROWSER] Richiesta "
                "apertura inviata"
            )

        else:

            print(
                "[BROWSER] Apri manualmente:",
                WEB_URL
            )

    print(
        "[STATO] "
        "READY_FOR_MANUAL_WEB_SUBMISSION"
    )

    print(
        "[NOTA] Nessun file inviato "
        "al server e nessuna predizione "
        "ancora eseguita"
    )


# ============================================================
# 9. RICONOSCIMENTO COLONNE CSV/TSV
# ============================================================

def column(headers, *aliases):

    def norm(value):

        return re.sub(
            "[^a-z0-9]",
            "",
            str(value).lower()
        )

    lookup = {
        norm(x): x
        for x in headers
    }

    return next(
        (
            lookup[a]
            for a in aliases
            if a in lookup
        ),
        None
    )


# ============================================================
# 10. LETTURA RISULTATI ESTERNI
# ============================================================

def load_external(path):

    if path.suffix.lower() not in {
        ".csv",
        ".tsv",
        ".txt"
    }:

        raise ValueError(
            "Salva i risultati in TSV/CSV/TXT, "
            "non HTML/XLSX"
        )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as file:

        sample = file.read(8192)

        file.seek(0)

        try:

            delim = csv.Sniffer().sniff(
                sample,
                delimiters="\t,;"
            ).delimiter

        except csv.Error:

            delim = (
                "\t"
                if "\t" in sample
                else ","
            )

        reader = csv.DictReader(
            file,
            delimiter=delim
        )

        headers = reader.fieldnames or []

        data = list(reader)

    srna_col = column(
        headers,
        "mirnaacc",
        "srnaacc",
        "srnaid",
        "mirna",
        "srna"
    )

    target_col = column(
        headers,
        "targetacc",
        "targetid",
        "target"
    )

    score_col = column(
        headers,
        "expectation",
        "expect"
    )

    if not (
        srna_col
        and target_col
        and score_col
    ):

        raise ValueError(
            "Colonne psRNATarget non individuate. "
            "Colonne trovate: "
            + repr(headers)
        )

    return (
        data,
        (
            srna_col,
            target_col,
            score_col
        ),
        headers
    )


# ============================================================
# 11. IMPORTAZIONE RISULTATI psRNATarget
# ============================================================

def import_results(
    source,
    run_dir=None,
    job_url=None
):

    runs = sorted(
        p
        for p in RUNS.glob(
            OUTPUT_PREFIX + "*"
        )
        if p.is_dir()
    )

    folder = (
        run_dir
        if run_dir is not None
        else (
            runs[-1]
            if runs
            else None
        )
    )

    if (
        folder is None
        or not (
            folder / "MANIFEST_STEP23.json"
        ).is_file()
    ):

        raise FileNotFoundError(
            "Sessione STEP23 non trovata"
        )

    manifest = json.loads(
        (
            folder / "MANIFEST_STEP23.json"
        ).read_text(
            encoding="utf-8"
        )
    )

    # --------------------------------------------------------
    # VERIFICA FILE DELLA SESSIONE
    # --------------------------------------------------------

    for info in manifest["files"].values():

        path = Path(
            info["path"]
        )

        if (
            not path.is_file()
            or sha256(path) != info["sha256"]
        ):

            raise ValueError(
                "File di sessione alterato: "
                + str(path)
            )

    srna = fasta(
        Path(
            manifest["files"]["standard"]["path"]
        )
    )

    targets = fasta(
        Path(
            manifest["files"]["targets"]["path"]
        )
    )

    _, local_sites = tsv(
        Path(
            manifest["files"]["sites"]["path"]
        )
    )

    indexed = defaultdict(list)

    tx_genes = {}

    for row in local_sites:

        indexed[
            (
                row["srna_id"],
                row["transcript_id"]
            )
        ].append(row)

        tx_genes[
            row["transcript_id"]
        ] = row.get(
            "gene_id",
            ""
        )

    # --------------------------------------------------------
    # LETTURA PREDIZIONI DEL SERVER
    # --------------------------------------------------------

    (
        data,
        (
            srna_col,
            target_col,
            score_col
        ),
        headers
    ) = load_external(source)

    folder_out = folder / (
        "IMPORT_RISULTATI_" + stamp()
    )

    folder_out.mkdir(
        parents=False,
        exist_ok=False
    )

    original = folder_out / (
        "RISULTATI_SERVER_ORIGINALI"
        + source.suffix.lower()
    )

    shutil.copyfile(
        source,
        original
    )

    output = folder_out / (
        "CONFRONTO_PSRNATARGET_STEP22.tsv"
    )

    cols = [

        "srna_id",

        "transcript_id",

        "gene_id",

        "expectation_psRNATarget",

        "target_start",

        "target_end",

        "confronto_step22",

        "migliore_full_score_locale"
    ]

    results = []

    stats = Counter()

    start_col = column(
        headers,
        "targetstart"
    )

    end_col = column(
        headers,
        "targetend"
    )

    # --------------------------------------------------------
    # CONFRONTO TRA SERVER E STEP22
    # --------------------------------------------------------

    for row in data:

        ident = row[
            srna_col
        ].strip()

        tx = row[
            target_col
        ].strip()

        if (
            ident not in srna
            or tx not in targets
        ):

            raise ValueError(
                "Risultati non compatibili "
                "con la sessione: "
                + repr((ident, tx))
            )

        val = (
            row[score_col]
            .strip()
            .replace(",", ".")
        )

        try:

            score = Decimal(val)

            if not score.is_finite():

                raise InvalidOperation

        except (
            InvalidOperation,
            ValueError
        ) as error:

            raise ValueError(
                "Expectation non numerico: "
                + repr(val)
            ) from error

        # ----------------------------------------------------
        # COORDINATE DEL SITO
        # ----------------------------------------------------

        a = (
            int(row[start_col])
            if (
                start_col
                and row[start_col]
                .strip()
                .isdigit()
            )
            else None
        )

        b = (
            int(row[end_col])
            if (
                end_col
                and row[end_col]
                .strip()
                .isdigit()
            )
            else None
        )

        if (
            a is not None
            and b is not None
            and not (
                1 <= a <= b <= len(targets[tx])
            )
        ):

            raise ValueError(
                "Coordinate fuori dal trascritto: "
                + repr((ident, tx, a, b))
            )

        local = indexed.get(
            (ident, tx),
            []
        )

        # ----------------------------------------------------
        # CLASSIFICAZIONE CORRISPONDENZE
        # ----------------------------------------------------

        if not local:

            level = (
                "altro_trascritto_preselezionato"
            )

        elif (
            a is None
            or b is None
        ):

            level = (
                "stesso_trascritto_"
                "coordinate_non_disponibili"
            )

        elif any(
            a == int(
                site["start_1based"]
            )
            and b == int(
                site["end_1based"]
            )
            for site in local
        ):

            level = (
                "stesso_sito_coordinate_esatte"
            )

        elif any(
            max(
                a,
                int(site["start_1based"])
            )
            <= min(
                b,
                int(site["end_1based"])
            )
            for site in local
        ):

            level = (
                "sito_sovrapposto"
            )

        else:

            level = (
                "altro_sito_stesso_trascritto"
            )

        stats[level] += 1

        # ----------------------------------------------------
        # RISULTATO INTEGRATO
        # ----------------------------------------------------

        results.append({

            "srna_id":
                ident,

            "transcript_id":
                tx,

            "gene_id":
                (
                    local[0]["gene_id"]
                    if local
                    else tx_genes.get(tx, "")
                ),

            "expectation_psRNATarget":
                str(score),

            "target_start":
                (
                    ""
                    if a is None
                    else a
                ),

            "target_end":
                (
                    ""
                    if b is None
                    else b
                ),

            "confronto_step22":
                level,

            "migliore_full_score_locale":
                min(
                    (
                        Decimal(
                            site["full_score"]
                        )
                        for site in local
                    ),
                    default=""
                )
        })

    # ========================================================
    # 12. ESPORTAZIONE CONFRONTO
    # ========================================================

    with output.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=cols,
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(
            results
        )

    # ========================================================
    # 13. RAPPORTO FINALE
    # ========================================================

    report = {

        "status":
            "IMPORTED_EXTERNAL_RESULTS",

        "source_job_url":
            job_url,

        "original_sha256":
            sha256(original),

        "result_rows":
            len(data),

        "classification":
            dict(stats),

        "schema":
            "Schema V2 dichiarato nella sessione; "
            "verificare i parametri effettivi "
            "utilizzati sul server",

        "important":
            "psRNATarget Expectation e "
            "punteggio STEP20 "
            "non sono intercambiabili",

        "output":
            str(output)
    }

    report_path = (
        folder_out / "RAPPORTO_CONFRONTO.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ) + "\n",
        encoding="utf-8"
    )

    # ========================================================
    # 14. OUTPUT TERMINALE
    # ========================================================

    print(
        "[RISULTATI IMPORTATI]",
        len(data)
    )

    print(
        "[CONFRONTO]",
        dict(stats)
    )

    print(
        "[TSV]",
        output
    )

    print(
        "[REPORT]",
        report_path
    )


# ============================================================
# 15. PROGRAMMA PRINCIPALE
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

    group.add_argument(
        "--audit",
        action="store_true"
    )

    group.add_argument(
        "--run",
        action="store_true"
    )

    group.add_argument(
        "--import-results",
        type=Path,
        metavar="RISULTATI_TSV_CSV"
    )

    parser.add_argument(
        "--no-browser",
        action="store_true"
    )

    parser.add_argument(
        "--run-dir",
        type=Path,
        help="Sessione STEP23 da associare ai risultati"
    )

    parser.add_argument(
        "--job-url",
        help="URL del lavoro psRNATarget da registrare"
    )

    args = parser.parse_args()

    print("=" * 68)

    print(
        "AI PRIMETTA V3.3 | STEP 23 "
        "- psRNATarget schema V2"
    )

    print("=" * 68)

    if args.audit:

        check_step22()

        print(
            "[OK] Audit in sola lettura superato"
        )

    elif args.run:

        prepare(
            not args.no_browser
        )

    else:

        if not args.import_results.is_file():

            raise FileNotFoundError(
                args.import_results
            )

        import_results(

            args.import_results.resolve(),

            (
                args.run_dir.resolve()
                if args.run_dir
                else None
            ),

            args.job_url
        )


# ============================================================
# 16. AVVIO
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            f"[ERRORE] "
            f"{type(exc).__name__}: {exc}",
            file=sys.stderr
        )

        sys.exit(1)