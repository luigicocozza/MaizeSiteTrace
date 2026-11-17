#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 23B

Alternativa locale al server psRNATarget:
screening delle 2 small RNA reali da 24 e 25 nt
contro 37 trascritti B73 tramite RNAplex.

--audit: verifica software e file STEP22.
--run: analisi RNAplex locale e report.

Non modifica GUI, Docker, checkpoint o dati originali.
RNAplex non e' psRNATarget.
"""

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile

from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. PERCORSI
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

STEP22 = (
    RUNS
    / "STEP22_EXPORT_TARGET_B73_20261023T113609974080Z"
)

MANIFEST = STEP22 / "MANIFEST_STEP22.json"

SRNA_FILE = (
    STEP22 / "SRNA_STANDARD_MAX25.fa"
)

TARGET_FILE = (
    STEP22 / "TRASCRITTI_B73_CANDIDATI.fa"
)

EXPECTED = {
    "sRNA_00000435": 24,
    "sRNA_00000037": 25
}


# ============================================================
# 2. RICONOSCIMENTO RISULTATI RNAplex
# ============================================================

PATTERN = re.compile(
    r"(?P<structure>\S+)\s+"
    r"(?P<from1>\d+),(?P<to1>\d+)\s*:\s*"
    r"(?P<from2>\d+),(?P<to2>\d+)\s*"
    r"\(\s*(?P<energy>[+-]?\d+(?:\.\d+)?)"
)


# ============================================================
# 3. SHA256
# ============================================================

def checksum(path):

    h = hashlib.sha256()

    with path.open("rb") as file:

        for chunk in iter(
            lambda: file.read(1024 * 1024),
            b""
        ):

            h.update(chunk)

    return h.hexdigest()


# ============================================================
# 4. LETTURA FASTA
# ============================================================

def fasta(path):

    records = {}

    name = None
    pieces = []

    text = path.read_text(
        encoding="utf-8-sig"
    )

    for line in text.splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if name is not None:

                records[name] = (
                    "".join(pieces)
                    .upper()
                    .replace("T", "U")
                )

            name = line[1:].split()[0]

            if (
                not name
                or name in records
            ):

                raise ValueError(
                    "Identificativo FASTA "
                    "vuoto oppure duplicato"
                )

            pieces = []

        elif name is None:

            raise ValueError(
                "FASTA senza intestazione"
            )

        else:

            pieces.append(line)

    if name is not None:

        records[name] = (
            "".join(pieces)
            .upper()
            .replace("T", "U")
        )

    if not records:

        raise ValueError(
            "FASTA vuoto"
        )

    for ident, sequence in records.items():

        if (
            not sequence
            or set(sequence) - set("ACGUN")
        ):

            raise ValueError(
                "Sequenza non valida: "
                + ident
            )

    return records


# ============================================================
# 5. AUDIT SCIENTIFICO
# ============================================================

def audit():

    executable = shutil.which(
        "RNAplex"
    )

    if not executable:

        raise RuntimeError(
            "RNAplex non installato. "
            "Eseguire: sudo apt install vienna-rna"
        )

    version = subprocess.run(
        [
            executable,
            "--version"
        ],
        capture_output=True,
        text=True,
        timeout=15,
        check=True
    ).stdout.strip()

    manifest = json.loads(
        MANIFEST.read_text(
            encoding="utf-8"
        )
    )

    for key, path in (
        ("standard", SRNA_FILE),
        ("targets", TARGET_FILE)
    ):

        if (
            checksum(path)
            != manifest["files_sha256"][key]
        ):

            raise ValueError(
                "Checksum diverso da STEP22: "
                + str(path)
            )

    srnas = fasta(
        SRNA_FILE
    )

    targets = fasta(
        TARGET_FILE
    )

    lengths = {
        ident: len(seq)
        for ident, seq in srnas.items()
    }

    if lengths != EXPECTED:

        raise ValueError(
            "Small RNA o lunghezze "
            "non corrispondenti"
        )

    if len(targets) != 37:

        raise ValueError(
            "Numero trascritti inatteso: "
            + str(len(targets))
        )

    if any(
        set(seq) - set("ACGU")
        for seq in srnas.values()
    ):

        raise ValueError(
            "Una small RNA contiene "
            "basi non canoniche"
        )

    print(
        "[VERSIONE]",
        version
    )

    print(
        "[PASS] Manifest STEP22 "
        "e SHA256 verificati"
    )

    print(
        "[PASS] 2 small RNA reali "
        "e 37 trascritti B73"
    )

    print(
        "[METODO] RNAplex locale "
        "senza profili di accessibilita'"
    )

    return (
        executable,
        version,
        srnas,
        targets
    )


# ============================================================
# 6. ANALISI RNAplex
# ============================================================

def analyze():

    (
        executable,
        version,
        srnas,
        targets
    ) = audit()

    timestamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    output = RUNS / (
        "STEP23B_RNAPLEX_STANDARD_"
        + timestamp
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    raw_path = (
        output / "OUTPUT_GREZZO_RNAPLEX.txt"
    )

    table_path = (
        output / "INTERAZIONI_RNAPLEX.tsv"
    )

    results = []
    failed = []

    total = 0

    # ========================================================
    # 7. FILE TEMPORANEI SEPARATI
    # ========================================================

    with tempfile.TemporaryDirectory(
        prefix="primetta_rnaplex_"
    ) as directory:

        with raw_path.open(
            "w",
            encoding="utf-8"
        ) as raw:

            folder = Path(directory)

            qfile = folder / "query.fa"
            tfile = folder / "target.fa"

            # Confronto completo:
            # 2 small RNA x 37 trascritti.

            for sid, query in sorted(
                srnas.items()
            ):

                for tid, target in sorted(
                    targets.items()
                ):

                    total += 1

                    qfile.write_text(
                        f">{sid}\n{query}\n",
                        encoding="utf-8"
                    )

                    tfile.write_text(
                        f">{tid}\n{target}\n",
                        encoding="utf-8"
                    )

                    command = [
                        executable,
                        "-q",
                        str(qfile),
                        "-t",
                        str(tfile)
                    ]

                    # ========================================
                    # 8. ESECUZIONE SOFTWARE SCIENTIFICO
                    # ========================================

                    try:

                        process = subprocess.run(
                            command,
                            capture_output=True,
                            text=True,
                            timeout=90,
                            check=False
                        )

                        block = process.stdout

                        stderr = process.stderr

                        returncode = (
                            process.returncode
                        )

                    except subprocess.TimeoutExpired:

                        block = ""

                        stderr = "TIMEOUT 90s"

                        returncode = -1

                    # ========================================
                    # 9. SALVATAGGIO OUTPUT GREZZO
                    # ========================================

                    raw.write(
                        f"\n=== {sid} :: {tid}"
                        f" | exit={returncode} ===\n"
                    )

                    raw.write(
                        block + "\n"
                    )

                    if stderr:

                        raw.write(
                            "STDERR: "
                            + stderr
                            + "\n"
                        )

                    if returncode != 0:

                        failed.append({

                            "srna_id": sid,

                            "transcript_id": tid,

                            "exit": returncode,

                            "stderr": stderr[:500]
                        })

                        continue

                    # ========================================
                    # 10. ESTRAZIONE ENERGIE
                    # ========================================

                    for match in PATTERN.finditer(
                        block
                    ):

                        obj = (
                            match.groupdict()
                        )

                        results.append({

                            "srna_id": sid,

                            "transcript_id": tid,

                            "interaction_energy_kcal_mol":
                                obj["energy"],

                            "raw_range_1":
                                (
                                    obj["from1"]
                                    + ","
                                    + obj["to1"]
                                ),

                            "raw_range_2":
                                (
                                    obj["from2"]
                                    + ","
                                    + obj["to2"]
                                ),

                            "structure":
                                obj["structure"]
                        })

                    if total % 10 == 0:

                        print(
                            f"[SCANSIONE] "
                            f"{total}/74 coppie",
                            flush=True
                        )

    # ========================================================
    # 11. RISULTATI TABELLARI
    # ========================================================

    columns = [

        "srna_id",

        "transcript_id",

        "interaction_energy_kcal_mol",

        "raw_range_1",

        "raw_range_2",

        "structure"
    ]

    results.sort(
        key=lambda row: (

            float(
                row["interaction_energy_kcal_mol"]
            ),

            row["srna_id"],

            row["transcript_id"]
        )
    )

    with table_path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(
            results
        )

    # ========================================================
    # 12. RAPPORTO SCIENTIFICO
    # ========================================================

    report = {

        "step": "23B",

        "method":
            "RNAplex minimum free energy",

        "viennarna_version":
            version,

        "temperature_C": 37,

        "accessibility":
            "not evaluated",

        "psRNATarget_equivalent":
            False,

        "pairs_analyzed":
            total,

        "result_lines_parsed":
            len(results),

        "pair_failures":
            failed,

        "source_manifest":
            str(MANIFEST),

        "source_sha256": {

            "srna":
                checksum(SRNA_FILE),

            "targets":
                checksum(TARGET_FILE)
        },

        "raw_output":
            str(raw_path),

        "results_tsv":
            str(table_path),

        "limitations": [

            "RNAplex is not psRNATarget",

            "MFE does not prove biological targeting",

            "Target accessibility not evaluated",

            "Coordinates preserved as raw ranges",

            "Biological validation still required"
        ]
    }

    report_path = (
        output / "RAPPORTO_STEP23B.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ) + "\n",
        encoding="utf-8"
    )

    # ========================================================
    # 13. OUTPUT TERMINALE
    # ========================================================

    print(
        "[ESITO]",
        total,
        "coppie;",
        len(results),
        "interazioni lette;",
        len(failed),
        "errori"
    )

    print(
        "[TSV]",
        table_path
    )

    print(
        "[REPORT]",
        report_path
    )

    print(
        "[RAW]",
        raw_path
    )

    if failed:

        raise RuntimeError(
            "Alcune analisi RNAplex "
            "non sono terminate correttamente"
        )

    if not results:

        print(
            "[ATTENZIONE] "
            "Nessuna linea riconosciuta. "
            "Verificare l'output grezzo: "
            "non concludere che manchino target."
        )

    print(
        "[OK] GUI, Docker e modelli invariati"
    )


# ============================================================
# 14. PROGRAMMA PRINCIPALE
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

    args = parser.parse_args()

    print("=" * 69)

    print(
        "AI PRIMETTA V3.3 | STEP 23B"
    )

    print(
        "PREDIZIONE LOCALE RNAplex"
    )

    print("=" * 69)

    if args.audit:

        audit()

        print(
            "[OK] Audit completato "
            "senza modifiche"
        )

    else:

        analyze()


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