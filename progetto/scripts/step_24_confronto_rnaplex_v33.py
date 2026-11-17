#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AI PRIMETTA V3.3 | STEP 24 - confronto RNAplex vs STEP22.

Eseguire con --audit oppure --run.
Il programma NON riesegue RNAplex e NON cambia GUI, Docker o modelli.
"""

import argparse
import csv
import hashlib
import json
import re
import sys

from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "target_prediction_v33_vscode" / "runs"

EXPECTED_RNA = {
    "sRNA_00000435": 24,
    "sRNA_00000037": 25
}

HEADER_RE = re.compile(
    r"^===\s+(\S+)\s+::\s+(\S+)\s+\|\s+exit=(-?\d+)\s+===$",
    re.M
)

HIT_RE = re.compile(
    r"(?P<structure>\S+)\s+"
    r"(?P<a>\d+),(?P<b>\d+)\s*:\s*"
    r"(?P<c>\d+),(?P<d>\d+)\s*"
    r"\(\s*(?P<energy>[+-]?\d+(?:\.\d+)?)"
)


# ============================================================
# 1. CHECKSUM
# ============================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# ============================================================
# 2. PROTEZIONE PERCORSI
# ============================================================

def inside_runs(path):

    p = Path(path)

    if not p.is_absolute():
        p = ROOT / p

    p = p.resolve()

    if (
        not p.is_file()
        or not p.is_relative_to(RUNS.resolve())
    ):

        raise ValueError(
            f"File assente o fuori da runs: {p}"
        )

    return p


# ============================================================
# 3. CONVERSIONE NUMERICA
# ============================================================

def decimal(value):

    try:

        x = Decimal(
            str(value).strip().replace(",", ".")
        )

        if not x.is_finite():

            raise ValueError("non finito")

        return x

    except (InvalidOperation, ValueError) as exc:

        raise ValueError(
            f"Numero non valido: {value!r}"
        ) from exc


# ============================================================
# 4. LETTURA FASTA
# ============================================================

def fasta(path):

    records = {}
    name = None
    chunks = []

    for line in path.read_text(
        encoding="utf-8-sig"
    ).splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if name is not None:

                if not chunks:

                    raise ValueError(
                        f"FASTA vuoto: {name}"
                    )

                records[name] = (
                    "".join(chunks)
                    .upper()
                    .replace("T", "U")
                )

            name = line[1:].split()[0]
            chunks = []

            if name in records or not name:

                raise ValueError(
                    f"ID FASTA ripetuto/vuoto: {name}"
                )

        else:

            if name is None:

                raise ValueError(
                    "FASTA senza intestazione"
                )

            chunks.append(line)

    if name is not None:

        if not chunks or name in records:

            raise ValueError(
                f"Ultimo record FASTA non valido: {name}"
            )

        records[name] = (
            "".join(chunks)
            .upper()
            .replace("T", "U")
        )

    if not records:

        raise ValueError(
            f"FASTA senza sequenze: {path}"
        )

    return records


# ============================================================
# 5. LETTURA TSV
# ============================================================

def read_tsv(path):

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        reader = csv.DictReader(
            stream,
            delimiter="\t"
        )

        names = reader.fieldnames or []

        if (
            not names
            or len(names) != len(set(names))
        ):

            raise ValueError(
                f"Intestazioni TSV non valide: {path}"
            )

        return names, list(reader)


# ============================================================
# 6. VERIFICA INTERVALLI
# ============================================================

def interval(text, n):

    if not re.fullmatch(
        r"\d+,\d+",
        str(text).strip()
    ):

        return None

    a, b = (
        int(v)
        for v in text.strip().split(",")
    )

    return (
        (a, b)
        if 1 <= a <= b <= n
        else None
    )


# ============================================================
# 7. LETTURA OUTPUT GREZZO RNAplex
# ============================================================

def raw_output(path):

    text = path.read_text(
        encoding="utf-8"
    )

    markers = list(
        HEADER_RE.finditer(text)
    )

    if not markers:

        raise ValueError(
            "Output grezzo: nessun blocco riconosciuto"
        )

    found = {}

    for i, marker in enumerate(markers):

        key = (
            marker.group(1),
            marker.group(2)
        )

        if (
            key in found
            or marker.group(3) != "0"
        ):

            raise ValueError(
                f"Blocco RNAplex duplicato o fallito: {key}"
            )

        end = (
            markers[i + 1].start()
            if i + 1 < len(markers)
            else len(text)
        )

        hits = list(
            HIT_RE.finditer(
                text[marker.end():end]
            )
        )

        if len(hits) != 1:

            raise ValueError(
                f"Numero hit grezzi != 1 "
                f"per {key}: {len(hits)}"
            )

        m = hits[0]

        found[key] = (
            f"{m['a']},{m['b']}",
            f"{m['c']},{m['d']}",
            decimal(m["energy"]),
            m["structure"]
        )

    return found


# ============================================================
# 8. CARICAMENTO DATI STEP23B
# ============================================================

def collect():

    reports = sorted(
        RUNS.glob(
            "STEP23B_RNAPLEX_STANDARD_*/"
            "RAPPORTO_STEP23B.json"
        )
    )

    if not reports:

        raise FileNotFoundError(
            "Manca il rapporto dello STEP23B"
        )

    report_path = inside_runs(
        reports[-1]
    )

    report = json.loads(
        report_path.read_text(
            encoding="utf-8"
        )
    )

    if (
        report.get("step") != "23B"
        or report.get("pairs_analyzed") != 74
        or report.get("pair_failures")
    ):

        raise ValueError(
            "STEP23B: le 74 analisi senza errori "
            "non sono confermate"
        )

    # ========================================================
    # 9. MANIFEST E PERCORSI
    # ========================================================

    manifest_path = inside_runs(
        report["source_manifest"]
    )

    manifest = json.loads(
        manifest_path.read_text(
            encoding="utf-8"
        )
    )

    source = manifest_path.parent

    paths = {

        "srna": inside_runs(
            source / "SRNA_STANDARD_MAX25.fa"
        ),

        "targets": inside_runs(
            source / "TRASCRITTI_B73_CANDIDATI.fa"
        ),

        "sites": inside_runs(
            source / "SITI_E_RNASEQ_SNAPSHOT.tsv"
        ),

        "table": inside_runs(
            report_path.parent
            / "INTERAZIONI_RNAPLEX.tsv"
        ),

        "raw": inside_runs(
            report_path.parent
            / "OUTPUT_GREZZO_RNAPLEX.txt"
        )
    }

    # ========================================================
    # 10. INTEGRITA' FILE STEP22
    # ========================================================

    for label, field in [
        ("srna", "standard"),
        ("targets", "targets"),
        ("sites", "sites")
    ]:

        if sha256(paths[label]) != (
            manifest["files_sha256"][field]
        ):

            raise ValueError(
                f"Checksum STEP22 discordante: {label}"
            )

    # ========================================================
    # 11. INTEGRITA' FILE STEP23B
    # ========================================================

    for label, field in [
        ("srna", "srna"),
        ("targets", "targets")
    ]:

        if sha256(paths[label]) != (
            report["source_sha256"][field]
        ):

            raise ValueError(
                f"Checksum STEP23B discordante: {label}"
            )

    # ========================================================
    # 12. CARICAMENTO SEQUENZE REALI
    # ========================================================

    srnas = fasta(
        paths["srna"]
    )

    targets = fasta(
        paths["targets"]
    )

    if (
        {
            k: len(v)
            for k, v in srnas.items()
        } != EXPECTED_RNA
        or len(targets) != 37
    ):

        raise ValueError(
            "Attesi esattamente 2 sRNA "
            "(24/25 nt) e 37 trascritti"
        )

    # ========================================================
    # 13. VALIDAZIONE SITI STEP22
    # ========================================================

    _, step22 = read_tsv(
        paths["sites"]
    )

    if len(step22) != 37:

        raise ValueError(
            "Attesi 37 siti STEP22"
        )

    local = defaultdict(list)
    genes = {}

    for site in step22:

        sid = site["srna_id"]
        tx = site["transcript_id"]

        if tx not in targets:

            raise ValueError(
                f"Trascritto STEP22 assente: {tx}"
            )

        a = int(
            site["start_1based"]
        )

        b = int(
            site["end_1based"]
        )

        expected = (
            site["full_site_DNA_5to3"]
            .upper()
            .replace("T", "U")
        )

        if (
            not 1 <= a <= b <= len(targets[tx])
            or targets[tx][a - 1:b] != expected
        ):

            raise ValueError(
                f"Sito STEP22 non corrisponde "
                f"al FASTA: {sid}, {tx}"
            )

        if (
            tx in genes
            and genes[tx] != site["gene_id"]
        ):

            raise ValueError(
                f"Annotazione gene incoerente per {tx}"
            )

        genes[tx] = site["gene_id"]

        if sid in srnas:

            local[
                (sid, tx)
            ].append(site)

    if sum(
        len(v)
        for v in local.values()
    ) != 36:

        raise ValueError(
            "Siti STEP22 delle due sRNA "
            "diversi da 36"
        )

    # ========================================================
    # 14. VERIFICA DELLE 74 COPPIE RNAplex
    # ========================================================

    required = {
        "srna_id",
        "transcript_id",
        "interaction_energy_kcal_mol",
        "raw_range_1",
        "raw_range_2",
        "structure"
    }

    columns, table = read_tsv(
        paths["table"]
    )

    if not required <= set(columns):

        raise ValueError(
            "Colonne RNAplex TSV assenti: "
            + repr(
                sorted(
                    required - set(columns)
                )
            )
        )

    raw = raw_output(
        paths["raw"]
    )

    expected_pairs = {

        (sid, tx)

        for sid in srnas
        for tx in targets
    }

    if (
        len(table) != 74
        or set(raw) != expected_pairs
    ):

        raise ValueError(
            "Output RNAplex incompleto: "
            "attese 74 coppie"
        )

    # ========================================================
    # 15. IMPOSTAZIONE CORRETTA COORDINATE
    #
    # raw_range_1 = TARGET
    # raw_range_2 = QUERY (small RNA)
    #
    # Confrontiamo anche col RAW e con la lunghezza
    # effettiva dei trascritti e delle small RNA.
    # ========================================================

    seen = set()

    orientation_proved = False

    pair_rows = []

    by_srna = defaultdict(Counter)

    by_type = Counter()

    site_overlap = set()
    site_exact = set()

    for row in table:

        key = (
            row["srna_id"],
            row["transcript_id"]
        )

        if (
            key not in expected_pairs
            or key in seen
        ):

            raise ValueError(
                f"Coppia TSV sconosciuta/duplicata: {key}"
            )

        seen.add(key)

        r1 = row[
            "raw_range_1"
        ].strip()

        r2 = row[
            "raw_range_2"
        ].strip()

        energy = decimal(
            row["interaction_energy_kcal_mol"]
        )

        if (
            r1,
            r2,
            energy,
            row["structure"]
        ) != raw[key]:

            raise ValueError(
                f"TSV e output grezzo discordi "
                f"per {key}"
            )

        query, tx = key

        # PRIMA COPPIA DI COORDINATE: TRASCRITTO
        t = interval(
            r1,
            len(targets[tx])
        )

        # SECONDA COPPIA DI COORDINATE: SMALL RNA
        q = interval(
            r2,
            len(srnas[query])
        )

        if t is None or q is None:

            raise ValueError(
                "Coordinate TARGET:QUERY non valide "
                f"per {key}: {r1} : {r2}"
            )

        # Almeno un intervallo deve dimostrare
        # che la prima coppia non può riferirsi
        # alla small RNA.

        if interval(
            r1,
            len(srnas[query])
        ) is None:

            orientation_proved = True

        t1, t2 = t
        q1, q2 = q

        # ====================================================
        # 16. CONFRONTO CON SITI ORIGINALI STEP22
        # ====================================================

        matches = []

        for site in local.get(key, []):

            a = int(
                site["start_1based"]
            )

            b = int(
                site["end_1based"]
            )

            intersection = max(
                0,
                min(t2, b)
                - max(t1, a)
                + 1
            )

            matches.append(
                (
                    intersection,
                    t1 == a and t2 == b,
                    site
                )
            )

            if intersection:

                site_overlap.add(
                    (query, tx, a, b)
                )

                if t1 == a and t2 == b:

                    site_exact.add(
                        (query, tx, a, b)
                    )

        # ====================================================
        # 17. CLASSIFICAZIONE DELLA COPPIA
        # ====================================================

        if not matches:

            label = (
                "nessun_sito_STEP22_per_coppia"
            )

            nearest = None
            overlap = 0

        else:

            overlap, exact, nearest = max(
                matches,

                key=lambda item: (
                    item[1],
                    item[0],
                    -decimal(
                        item[2]["full_score"]
                    )
                )
            )

            if exact:

                label = "coordinate_esatte"

            elif overlap:

                label = "sovrapposizione"

            else:

                label = (
                    "altro_sito_stesso_trascritto"
                )

        by_type[label] += 1

        by_srna[query][label] += 1

        # ====================================================
        # 18. RECORD FINALE
        # ====================================================

        pair_rows.append({

            "srna_id":
                query,

            "transcript_id":
                tx,

            "gene_id":
                genes.get(tx, ""),

            "energy_kcal_mol":
                str(energy),

            "target_start_1based":
                t1,

            "target_end_1based":
                t2,

            "srna_start_1based":
                q1,

            "srna_end_1based":
                q2,

            "target_fragment_RNA_5to3":
                targets[tx][t1 - 1:t2],

            "RNAplex_structure":
                row["structure"],

            "comparison":
                label,

            "overlap_nt":
                overlap,

            "STEP22_start":
                (
                    nearest["start_1based"]
                    if nearest
                    else ""
                ),

            "STEP22_end":
                (
                    nearest["end_1based"]
                    if nearest
                    else ""
                ),

            "STEP22_full_score":
                (
                    nearest["full_score"]
                    if nearest
                    else ""
                )
        })

    # ========================================================
    # 19. CONTROLLI FINALI
    # ========================================================

    if (
        seen != expected_pairs
        or not orientation_proved
    ):

        raise ValueError(
            "Ordine TARGET:QUERY non verificabile "
            "univocamente: controllare output grezzo"
        )

    if (
        len(site_overlap) > 36
        or len(site_exact) > len(site_overlap)
    ):

        raise ValueError(
            "Conteggio siti originale incoerente"
        )

    return (
        report_path,
        paths,
        pair_rows,
        by_type,
        by_srna,
        site_overlap,
        site_exact
    )


# ============================================================
# 20. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    group = parser.add_mutually_exclusive_group(
        required=True
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

    print("=" * 67)

    print(
        "AI PRIMETTA V3.3 | STEP 24 "
        "- FIX RNAplex TARGET:QUERY"
    )

    print("=" * 67)

    (
        report,
        sources,
        rows,
        counts,
        by_srna,
        overlaps,
        exacts
    ) = collect()

    print(
        "[PASS] 74 coppie RNAplex; "
        "sorgenti FASTA e checksum verificati"
    )

    print(
        "[PASS] Coordinate TARGET:QUERY "
        "verificate con l'output RAW"
    )

    print(
        "[CLASSIFICAZIONE COPPIE]",
        json.dumps(
            dict(counts),
            ensure_ascii=False,
            sort_keys=True
        )
    )

    for sid in sorted(EXPECTED_RNA):

        c = by_srna[sid]

        print(
            f"[RISULTATO] {sid}: "
            f"esatte={c['coordinate_esatte']}, "
            f"sovrapposte={c['sovrapposizione']}, "
            f"altro_sito="
            f"{c['altro_sito_stesso_trascritto']}, "
            f"senza_sito="
            f"{c['nessun_sito_STEP22_per_coppia']}"
        )

    print(
        f"[SITI STEP22] 36 "
        f"| sovrapposti={len(overlaps)} "
        f"| esatti={len(exacts)}"
    )

    if args.audit:

        print(
            "[OK] Audit completato in sola lettura"
        )

        return

    # ========================================================
    # 21. CREAZIONE NUOVA CARTELLA OUTPUT
    # ========================================================

    timeid = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    dest = RUNS / (
        f"STEP24_CONFRONTO_RNAPLEX_CORRETTO_{timeid}"
    )

    dest.mkdir(
        parents=True,
        exist_ok=False
    )

    # ========================================================
    # 22. ESPORTAZIONE TSV
    # ========================================================

    table = (
        dest / "CONFRONTO_74_COPPIE.tsv"
    )

    rows.sort(
        key=lambda r: (
            decimal(r["energy_kcal_mol"]),
            r["srna_id"],
            r["transcript_id"]
        )
    )

    with table.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as out:

        writer = csv.DictWriter(
            out,
            fieldnames=list(rows[0]),
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(rows)

    # ========================================================
    # 23. RAPPORTO SCIENTIFICO JSON
    # ========================================================

    summary = {

        "step": "24",

        "status":
            "comparison_completed_not_biological_validation",

        "source_step23b":
            str(report.relative_to(ROOT)),

        "source_sha256": {

            name: sha256(path)
            for name, path in sources.items()
        },

        "RNAplex_pairs":
            len(rows),

        "original_step22_sites":
            36,

        "original_sites_overlapping":
            len(overlaps),

        "original_sites_exact":
            len(exacts),

        "pair_classification":
            dict(counts),

        "by_srna": {

            k: dict(v)
            for k, v in by_srna.items()
        },

        "RNAplex_interval_order":
            "TARGET:QUERY",

        "output_tsv":
            str(table.relative_to(ROOT)),

        "scientific_limits": [

            "One optimal interaction per pair; "
            "not all possible binding sites",

            "No target accessibility considered",

            "These 37 B73 transcripts "
            "were preselected in STEP22",

            "Coordinate overlap is not "
            "experimental target validation"
        ],

        "gui_modified":
            False,

        "docker_modified":
            False,

        "training_performed":
            False
    }

    # ========================================================
    # 24. SALVATAGGIO REPORT
    # ========================================================

    jpath = (
        dest / "RAPPORTO_STEP24.json"
    )

    jpath.write_text(
        json.dumps(
            summary,
            ensure_ascii=False,
            indent=2
        ) + "\n",
        encoding="utf-8"
    )

    print(
        "[TSV]",
        table
    )

    print(
        "[REPORT]",
        jpath
    )

    print(
        "[OK] File originali e GUI invariati; "
        "nessun RNAplex rieseguito"
    )


# ============================================================
# 25. AVVIO PROGRAMMA
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except Exception as exc:

        print(
            f"[ERRORE] {type(exc).__name__}: {exc}",
            file=sys.stderr
        )

        sys.exit(1)