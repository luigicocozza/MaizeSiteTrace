
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AI PRIMETTA V3.3 - STEP 22

Prepara i file per una verifica indipendente
con psRNATarget.

--audit:
    Verifica i dati senza generare file.

--prepare:
    Estrae trascritti completi dal B73 e crea
    i file FASTA in una nuova cartella.

Non modifica GUI, Docker, modelli, dati originali.
"""

import argparse
import csv
import gzip
import hashlib
import json
import sys

from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

SRNA_FASTA = (
    RUNS
    / "STEP16_5_SRNA_REALI_20261021T181110_648766Z"
    / "input_5_srna_reali.fa"
)

B73_FASTA = (
    ROOT
    / "target_prediction_v33_vscode"
    / "references"
    / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz"
)

EXPECTED = {
    "sRNA_00000002",
    "sRNA_00000040",
    "sRNA_00000007",
    "sRNA_00000435",
    "sRNA_00000037"
}

REQUIRED = {
    "srna_id",
    "transcript_id",
    "gene_id",
    "start_1based",
    "end_1based",
    "full_site_DNA_5to3",
    "full_score",
    "B73_qvalue",
    "H21_qvalue"
}


# ==========================================================
# 1. CHECKSUM
# ==========================================================

def checksum(path):

    h = hashlib.sha256()

    with path.open("rb") as file:

        for chunk in iter(
            lambda: file.read(1024 * 1024),
            b""
        ):
            h.update(chunk)

    return h.hexdigest()


# ==========================================================
# 2. LETTURA FASTA E FASTA.GZ
# ==========================================================

def fasta(path):

    opener = (
        gzip.open
        if path.suffix == ".gz"
        else open
    )

    with opener(
        path,
        "rt",
        encoding="utf-8"
    ) as file:

        name = None
        parts = []

        for raw in file:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    yield (
                        name,
                        "".join(parts)
                        .upper()
                        .replace("U", "T")
                    )

                name = line[1:].split()[0]
                parts = []

            else:

                if name is None:

                    raise ValueError(
                        "FASTA senza intestazione"
                    )

                parts.append(line)

        if name is not None:

            yield (
                name,
                "".join(parts)
                .upper()
                .replace("U", "T")
            )


# ==========================================================
# 3. INDIVIDUAZIONE STEP21B
# ==========================================================

def latest_step21b():

    candidates = sorted(
        RUNS.glob(
            "STEP21B_RNASEQ_CORRETTO_*/"
            "RAPPORTO_STEP21B.json"
        )
    )

    if not candidates:

        raise FileNotFoundError(
            "Manca il rapporto STEP21B"
        )

    report_path = candidates[-1]

    report = json.loads(
        report_path.read_text(
            encoding="utf-8"
        )
    )

    if report.get("step") != "21B":

        raise ValueError(
            "Rapporto STEP21B inatteso"
        )

    rows_file = next(
        (
            ROOT / item
            for item in report["outputs"]

            if Path(item).name
            == "SITI_PRIORITARI_RNASEQ_CORRETTI.tsv"
        ),
        None
    )

    if (
        rows_file is None
        or not rows_file.is_file()
    ):

        raise FileNotFoundError(
            "TSV STEP21B non trovato"
        )

    if not rows_file.resolve().is_relative_to(
        RUNS.resolve()
    ):

        raise ValueError(
            "TSV esterno alla cartella runs"
        )

    return (
        report_path,
        report,
        rows_file
    )


# ==========================================================
# 4. CARICAMENTO DATI REALI
# ==========================================================

def read_inputs():

    for path in (
        SRNA_FASTA,
        B73_FASTA
    ):

        if not path.is_file():

            raise FileNotFoundError(path)

    (
        report_path,
        report,
        tsv
    ) = latest_step21b()

    rnaseq_path = (
        ROOT / report["source_rnaseq"]
    )

    if (
        checksum(rnaseq_path)
        != report["source_rnaseq_sha256"]
    ):

        raise ValueError(
            "File RNA-seq modificato "
            "dopo STEP21B"
        )

    with tsv.open(
        encoding="utf-8-sig",
        newline=""
    ) as file:

        reader = csv.DictReader(
            file,
            delimiter="\t"
        )

        columns = reader.fieldnames or []

        if not REQUIRED.issubset(columns):

            raise ValueError(
                "Colonne STEP21B mancanti: "
                + str(
                    REQUIRED - set(columns)
                )
            )

        sites = list(reader)

    if len(sites) != report[
        "candidate_sites"
    ]:

        raise ValueError(
            "Numero siti diverso dal report"
        )

    # ======================================================
    # 5. CINQUE SMALL RNA ORIGINALI
    # ======================================================

    ids = {}

    for header, sequence in fasta(
        SRNA_FASTA
    ):

        matches = [
            item
            for item in EXPECTED
            if item in header
        ]

        if (
            len(matches) != 1
            or matches[0] in ids
        ):

            raise ValueError(
                "Identificativo sRNA ambiguo: "
                + header
            )

        if not set(sequence) <= set(
            "ACGT"
        ):

            raise ValueError(
                "Sequenza sRNA non canonica: "
                + header
            )

        ids[matches[0]] = sequence

    if set(ids) != EXPECTED:

        raise ValueError(
            "Insieme delle sRNA diverso "
            "dalle cinque previste"
        )

    # ======================================================
    # 6. COERENZA IDENTIFICATIVI E COORDINATE
    # ======================================================

    keys = set()

    for row in sites:

        sid = row["srna_id"]

        if sid not in ids:

            raise ValueError(
                "sRNA inattesa: " + sid
            )

        start = int(
            row["start_1based"]
        )

        end = int(
            row["end_1based"]
        )

        if (
            start < 1
            or end - start + 1
            != len(ids[sid])
        ):

            raise ValueError(
                "Coordinate incoerenti per "
                + sid
            )

        key = (
            sid,
            row["transcript_id"],
            start,
            end
        )

        if key in keys:

            raise ValueError(
                "Sito ripetuto: "
                + str(key)
            )

        keys.add(key)

    return (
        report_path,
        report,
        tsv,
        columns,
        sites,
        ids
    )


# ==========================================================
# 7. SCRITTURA FASTA
# ==========================================================

def write_fasta(path, pairs):

    with path.open(
        "w",
        encoding="utf-8"
    ) as output:

        for name, sequence in pairs:

            output.write(
                ">" + name + "\n"
            )

            for n in range(
                0,
                len(sequence),
                80
            ):

                output.write(
                    sequence[n:n + 80]
                    + "\n"
                )


# ==========================================================
# 8. AUDIT E PREPARAZIONE
# ==========================================================

def run(prepare):

    (
        report_path,
        report,
        tsv,
        columns,
        sites,
        srnas
    ) = read_inputs()

    counts = Counter(
        row["srna_id"]
        for row in sites
    )

    target_ids = {
        row["transcript_id"]
        for row in sites
    }

    genes = {
        row["gene_id"]
        for row in sites
        if row["gene_id"]
    }

    print(
        "[STEP21B]",
        report_path.relative_to(ROOT)
    )

    print(
        "[SITI]",
        len(sites),
        "| trascritti distinti:",
        len(target_ids),
        "| geni annotati:",
        len(genes)
    )

    print(
        "[JOIN RNA-seq]",
        report["sites_exact_gene_match"],
        "/",
        len(sites)
    )

    standard = []
    custom = []

    for ident, sequence in sorted(
        srnas.items()
    ):

        eligible = len(sequence) <= 25

        print(
            f"[sRNA] {ident}"
            f" | {len(sequence)} nt"
            f" | siti={counts[ident]}"
            f" | schema standard="
            f"{'SI' if eligible else 'NO'}"
        )

        if eligible:

            standard.append(
                (ident, sequence)
            )

        else:

            custom.append(
                (ident, sequence)
            )

    if (
        len(standard) != 2
        or len(custom) != 3
        or len(genes) != report["unique_genes"]
    ):

        raise ValueError(
            "Statistiche inattese "
            "rispetto allo STEP21B"
        )

    if not prepare:

        print(
            "[OK] AUDIT completato: "
            "nessun file scritto"
        )

        return

    # ======================================================
    # 9. ESTRAZIONE TRASCRITTI COMPLETI B73
    # ======================================================

    targets = {}
    found = set()

    for tx, sequence in fasta(
        B73_FASTA
    ):

        if tx not in target_ids:
            continue

        if tx in found:

            raise ValueError(
                "Trascritto duplicato in B73: "
                + tx
            )

        found.add(tx)

        if not (
            50 <= len(sequence) <= 5_000_000
        ):

            raise ValueError(
                "Lunghezza trascritto "
                "non ammessa: " + tx
            )

        if not set(sequence) <= set(
            "ACGTN"
        ):

            raise ValueError(
                "Trascritto con basi "
                "non riconosciute: " + tx
            )

        targets[tx] = sequence

    if found != target_ids:

        raise ValueError(
            "Trascritti assenti da B73: "
            + str(
                sorted(
                    target_ids - found
                )
            )
        )

    # ======================================================
    # 10. VERIFICA COORDINATE ORIGINALI
    # ======================================================

    for row in sites:

        tx = row["transcript_id"]

        a = int(
            row["start_1based"]
        )

        b = int(
            row["end_1based"]
        )

        observed = targets[tx][
            a - 1:b
        ]

        expected = (
            row["full_site_DNA_5to3"]
            .upper()
            .replace("U", "T")
        )

        if observed != expected:

            raise ValueError(
                "Il sito TSV non coincide "
                "col riferimento B73: "
                + row["srna_id"]
                + " / "
                + tx
            )

    # ======================================================
    # 11. VERIFICA SIGNIFICATIVITA' RNA-SEQ
    # ======================================================

    significant_counts = {}

    for prefix in ("B73", "H21"):

        counted = sum(
            row.get(
                f"{prefix}_significant_q_lt_005"
            ) == "yes"

            for row in sites
        )

        reported = (
            report.get(
                "site_level_q_significant_counts",
                {}
            )
            .get(prefix, {})
            .get("yes", 0)
        )

        if counted != reported:

            raise ValueError(
                "Significativita' "
                + prefix
                + " non coerente "
                "col report STEP21B"
            )

        significant_counts[
            f"{prefix}_qvalue_lt_005"
        ] = counted

    # ======================================================
    # 12. NUOVA CARTELLA, SENZA SOVRASCRITTURE
    # ======================================================

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    dest = RUNS / (
        "STEP22_EXPORT_TARGET_B73_"
        + stamp
    )

    dest.mkdir(
        parents=True,
        exist_ok=False
    )

    files = {

        "all_5":
            dest / "SRNA_TUTTE_5.fa",

        "standard":
            dest / "SRNA_STANDARD_MAX25.fa",

        "custom_hold":
            dest / "SRNA_26_28_SOLO_CUSTOM.fa",

        "targets":
            dest / "TRASCRITTI_B73_CANDIDATI.fa",

        "sites":
            dest / "SITI_E_RNASEQ_SNAPSHOT.tsv"
    }

    # ======================================================
    # 13. ESPORTAZIONE SMALL RNA
    # ======================================================

    write_fasta(
        files["all_5"],
        [
            (
                ident,
                seq.replace("T", "U")
            )
            for ident, seq in sorted(
                srnas.items()
            )
        ]
    )

    write_fasta(
        files["standard"],
        [
            (
                ident,
                seq.replace("T", "U")
            )
            for ident, seq in standard
        ]
    )

    write_fasta(
        files["custom_hold"],
        [
            (
                ident,
                seq.replace("T", "U")
            )
            for ident, seq in custom
        ]
    )

    # ======================================================
    # 14. ESPORTAZIONE TRASCRITTI B73
    # ======================================================

    write_fasta(
        files["targets"],
        sorted(
            targets.items()
        )
    )

    # ======================================================
    # 15. SNAPSHOT DATI SCIENTIFICI
    # ======================================================

    with files["sites"].open(
        "w",
        encoding="utf-8",
        newline=""
    ) as output:

        writer = csv.DictWriter(
            output,
            fieldnames=columns,
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(sites)

    # ======================================================
    # 16. MANIFEST DI PROVENIENZA
    # ======================================================

    manifest = {

        "step": "22",

        "status":
            "PREPARED_NOT_SUBMITTED",

        "source_step21b": str(
            report_path.relative_to(ROOT)
        ),

        "source_step21b_sha256":
            checksum(report_path),

        "source_sites_sha256":
            checksum(tsv),

        "reference": str(
            B73_FASTA.relative_to(ROOT)
        ),

        "reference_sha256":
            checksum(B73_FASTA),

        "candidate_sites":
            len(sites),

        "unique_genes":
            len(genes),

        "unique_transcripts":
            len(targets),

        "srna_standard_count":
            len(standard),

        "srna_custom_count":
            len(custom),

        "biological_status":
            "EXPLORATORY_NOT_VALIDATED",

        "statistics":
            significant_counts,

        "files_sha256": {
            name: checksum(path)
            for name, path in files.items()
        },

        "limitations": [

            "Selezione condizionata "
            "da HSP20 e score locale",

            "La DE non dimostra targeting",

            "La classe biologica delle "
            "sRNA non e' verificata",

            "psRNATarget standard non "
            "accetta sRNA oltre 25 nt"
        ]
    }

    manifest_path = (
        dest / "MANIFEST_STEP22.json"
    )

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
            ensure_ascii=False
        ) + "\n",
        encoding="utf-8"
    )

    print(
        "[OK] Le 37 coordinate "
        "coincidono con B73 NAM-5.0"
    )

    for key, path in files.items():

        print(
            f"[FILE {key}] {path}"
        )

    print(
        "[REPORT]",
        manifest_path
    )

    print(
        "[OK] Preparazione conclusa; "
        "nessun invio al server"
    )

    print(
        "[OK] GUI e Docker invariati"
    )


# ==========================================================
# 17. MAIN
# ==========================================================

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
        "--prepare",
        action="store_true"
    )

    args = parser.parse_args()

    print("=" * 68)

    print(
        "AI PRIMETTA V3.3 | STEP 22"
    )

    print(
        "PREPARAZIONE INPUT "
        "DI VERIFICA INDIPENDENTE"
    )

    print("=" * 68)

    run(args.prepare)


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
