
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 29

ACCESSIBILITA' RELATIVA PER TRASCRITTO

--audit:
    Verifica gli STEP 27 e 28 e ricalcola
    l'accessibilita' relativa senza scrivere file.

--run:
    Genera TSV per sito, TSV per gene
    e rapporto JSON in una nuova cartella.

--step28-run:
    Cartella STEP28 specifica, facoltativa.

--step27-run:
    Cartella STEP27 specifica, facoltativa.

Non riesegue RNAplfold/RNAplex.
Non modifica GUI, Docker o modelli ML/GNN.
"""

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

LENGTH = 10

SRNAS = {
    "sRNA_00000037": 17,
    "sRNA_00000435": 19
}


# ============================================================
# 2. CONTROLLO CONDIZIONI
# ============================================================

def check(condition, message):

    if not condition:
        raise ValueError(message)


# ============================================================
# 3. CONTROLLO PERCORSI
# ============================================================

def path_in_runs(value):

    p = Path(value)

    if not p.is_absolute():
        p = ROOT / p

    p = p.resolve()

    check(
        p.is_file()
        and p.is_relative_to(RUNS.resolve()),
        f"File assente/fuori da runs: {p}"
    )

    return p


# ============================================================
# 4. SHA256
# ============================================================

def sha(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for chunk in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):

            h.update(chunk)

    return h.hexdigest()


# ============================================================
# 5. LETTURA TSV
# ============================================================

def read_tsv(path, required):

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        fields = reader.fieldnames or []

        check(
            bool(fields)
            and len(fields) == len(set(fields))
            and set(required).issubset(fields),
            f"Colonne TSV inattese: {path}"
        )

        return list(reader)


# ============================================================
# 6. LETTURA FASTA
# ============================================================

def read_fasta(path):

    records = {}

    name = None
    parts = []

    for line in path.read_text(
        encoding="utf-8-sig"
    ).splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if name is not None:

                check(
                    bool(parts)
                    and name not in records,
                    "FASTA duplicato/vuoto"
                )

                records[name] = (
                    "".join(parts)
                    .upper()
                    .replace("T", "U")
                )

            name = line[1:].split()[0]

            parts = []

            check(
                bool(name)
                and name not in records,
                "ID FASTA errato"
            )

        else:

            check(
                name is not None,
                "FASTA senza intestazione"
            )

            parts.append(line)

    if name is not None:

        check(
            bool(parts)
            and name not in records,
            "Ultimo record FASTA errato"
        )

        records[name] = (
            "".join(parts)
            .upper()
            .replace("T", "U")
        )

    check(
        bool(records),
        "FASTA vuoto"
    )

    return records


# ============================================================
# 7. INDIVIDUAZIONE REPORT
# ============================================================

def report_path(folder, pattern, filename):

    if folder is not None:

        return path_in_runs(
            Path(folder) / filename
        )

    hits = sorted(
        RUNS.glob(
            f"{pattern}/{filename}"
        )
    )

    check(
        bool(hits),
        f"Non trovato {filename}"
    )

    return path_in_runs(
        hits[-1]
    )


# ============================================================
# 8. LETTURA PROBABILITA' RNAplfold
# ============================================================

def lunp_p10(path, sequence_length):

    """
    P10[i]:
    probabilita' che il segmento
    [i-9, i] sia interamente non appaiato.

    Coordinate 1-based.
    """

    values = (
        [None]
        * (sequence_length + 1)
    )

    seen = set()

    for line in path.read_text(
        encoding="utf-8-sig"
    ).splitlines():

        cols = line.strip().split()

        if (
            not cols
            or cols[0].startswith("#")
            or not cols[0].isdigit()
        ):
            continue

        position = int(
            cols[0]
        )

        check(
            1 <= position <= sequence_length
            and position not in seen,
            f"Indice _lunp errato: {position}"
        )

        seen.add(position)

        check(
            len(cols) > LENGTH,
            f"Manca colonna P10: {path}"
        )

        token = cols[LENGTH]

        if token.upper() != "NA":

            value = float(
                token
            )

            check(
                math.isfinite(value)
                and 0 <= value <= 1.000001,
                f"P10 fuori range: {token}"
            )

            values[position] = min(
                value,
                1.0
            )

    check(
        len(seen) == sequence_length,
        "Profilo _lunp incompleto: "
        f"{len(seen)}/{sequence_length}"
    )

    return values


# ============================================================
# 9. SOMME CUMULATIVE DEL PROFILO
# ============================================================

def profile_prefix(values):

    sums = [0.0]

    nvalid = [0]

    for value in values[1:]:

        sums.append(
            sums[-1]
            + (
                value
                if value is not None
                else 0.0
            )
        )

        nvalid.append(
            nvalid[-1]
            + (value is not None)
        )

    return sums, nvalid


# ============================================================
# 10. MEDIA P10 PER UN SITO
# ============================================================

def region_mean(prefix, start, end):

    """
    Media P10 di tutte le finestre da 10 nt
    interamente contenute in [start, end].
    """

    sums, valid = prefix

    first = (
        start
        + LENGTH
        - 1
    )

    expected = (
        end
        - first
        + 1
    )

    if (
        expected <= 0
        or (
            valid[end]
            - valid[first - 1]
        ) != expected
    ):

        return None

    return (
        sums[end]
        - sums[first - 1]
    ) / expected


# ============================================================
# 11. RECUPERO E CONTROLLO STEP28
# ============================================================

def load_source(
    path28_dir,
    path27_dir
):

    p28 = report_path(
        path28_dir,
        "STEP28_ACCESSIBILITA_RNAPLFOLD_*",
        "RAPPORTO_STEP28.json"
    )

    r28 = json.loads(
        p28.read_text(
            encoding="utf-8"
        )
    )

    check(
        str(r28.get("step")) == "28"
        and r28.get("sites_processed") == 36
        and r28.get("transcripts_processed") == 37
        and r28.get("unpaired_length") == 10,
        "Manifest STEP28 incompleto"
    )

    # ========================================================
    # 12. RISALITA ALLO STEP24
    # ========================================================

    p24 = path_in_runs(
        r28["source_step24"]
    )

    r24 = json.loads(
        p24.read_text(
            encoding="utf-8"
        )
    )

    check(
        str(r24.get("step")) == "24",
        "STEP24 non coerente"
    )

    p23 = path_in_runs(
        r24["source_step23b"]
    )

    r23 = json.loads(
        p23.read_text(
            encoding="utf-8"
        )
    )

    p22 = path_in_runs(
        r23["source_manifest"]
    )

    # ========================================================
    # 13. FILE SCIENTIFICI ORIGINALI
    # ========================================================

    sources = {

        "srna": path_in_runs(
            p22.parent
            / "SRNA_STANDARD_MAX25.fa"
        ),

        "targets": path_in_runs(
            p22.parent
            / "TRASCRITTI_B73_CANDIDATI.fa"
        ),

        "sites": path_in_runs(
            p22.parent
            / "SITI_E_RNASEQ_SNAPSHOT.tsv"
        ),

        "raw": path_in_runs(
            p23.parent
            / "OUTPUT_GREZZO_RNAPLEX.txt"
        ),

        "table": path_in_runs(
            p23.parent
            / "INTERAZIONI_RNAPLEX.tsv"
        )
    }

    for key, path in sources.items():

        check(
            sha(path)
            == r28["source_sha256"][key],
            f"Checksum STEP28 differente: {key}"
        )

    targets = read_fasta(
        sources["targets"]
    )

    check(
        len(targets) == 37,
        "Numero trascritti B73 errato"
    )

    # ========================================================
    # 14. TABELLA STEP28
    # ========================================================

    sites = read_tsv(

        path_in_runs(
            r28["output_tsv"]
        ),

        [
            "srna_id",
            "transcript_id",
            "gene_id",
            "target_start_1based",
            "target_end_1based",
            "windows10_expected",
            "windows10_available",
            "measurement_complete",
            "mean_P_unpaired_10nt",
            "RNAplex_MFE_kcal_mol"
        ]
    )

    check(
        len(sites) == 36
        and Counter(
            s["srna_id"]
            for s in sites
        ) == Counter(SRNAS),
        "Attesi 36 siti STEP28 (17+19)"
    )

    # ========================================================
    # 15. TABELLA STEP22
    # ========================================================

    step22 = read_tsv(

        sources["sites"],

        [
            "srna_id",
            "transcript_id",
            "gene_id",
            "start_1based",
            "end_1based",
            "full_site_DNA_5to3"
        ]
    )

    step22_index = {

        (
            s["srna_id"],
            s["transcript_id"]
        ): s

        for s in step22
    }

    check(
        len(step22) == 37
        and len(step22_index) == 37,
        "Siti STEP22 assenti o duplicati"
    )

    # ========================================================
    # 16. RECUPERO RISULTATI STEP27
    # ========================================================

    p27 = report_path(
        path27_dir,
        "STEP27_STATISTICA_DINUC_*",
        "RAPPORTO_STEP27.json"
    )

    r27 = json.loads(
        p27.read_text(
            encoding="utf-8"
        )
    )

    check(
        str(r27.get("step")) == "27",
        "Rapporto STEP27 errato"
    )

    p26 = path_in_runs(
        r27["source_step26"]
    )

    r26 = json.loads(
        p26.read_text(
            encoding="utf-8"
        )
    )

    p25 = path_in_runs(
        r26["source_step25"]
    )

    r25 = json.loads(
        p25.read_text(
            encoding="utf-8"
        )
    )

    # ========================================================
    # 17. COERENZA DELLE DUE CATENE SCIENTIFICHE
    # ========================================================

    check(
        path_in_runs(
            r25["source_step24"]
        ) == p24,
        "STEP27 e STEP28 provengono "
        "da sessioni STEP24 diverse"
    )

    for key, path in sources.items():

        check(
            sha(path)
            == r26["source_sha256"][key],
            "Origini STEP26 e STEP28 "
            f"discordanti: {key}"
        )

    # ========================================================
    # 18. RECUPERO TABELLA DEI Q-VALUE
    # ========================================================

    pair_files = [

        path_in_runs(p)

        for p in r27["outputs"]

        if Path(p).name
        == "CONFRONTO_74_PVALUE_QVALUE.tsv"
    ]

    check(
        len(pair_files) == 1,
        "Risultati STEP27 per coppia mancanti"
    )

    pairs = read_tsv(

        pair_files[0],

        [
            "srna_id",
            "transcript_id",
            "gene_id",
            "STEP22_candidate",
            "energy_original_kcal_mol",
            "p_empirical_31",
            "BH_q_74_exploratory"
        ]
    )

    pair_index = {

        (
            r["srna_id"],
            r["transcript_id"]
        ): r

        for r in pairs
    }

    check(
        len(pairs) == 74
        and len(pair_index) == 74,
        "Attese 74 coppie STEP27 distinte"
    )

    # ========================================================
    # 19. VERIFICA DEI 37 PROFILI ORIGINALI
    # ========================================================

    profiles = r28.get(
        "profiles",
        {}
    )

    check(
        set(profiles) == set(targets),
        "Profili STEP28 non corrispondenti "
        "ai trascritti"
    )

    profile_data = {}

    for tx, seq in sorted(
        targets.items()
    ):

        meta = profiles[tx]

        path = path_in_runs(
            meta["raw"]
        )

        check(
            path.parent
            == (
                p28.parent
                / "PROFILI_LUNP_ORIGINALI"
            ),
            "Profilo RNAplfold fuori "
            f"dalla sessione: {tx}"
        )

        check(
            sha(path) == meta["raw_sha256"],
            f"Profilo RNAplfold modificato: {tx}"
        )

        probabilities = lunp_p10(
            path,
            len(seq)
        )

        profile_data[tx] = profile_prefix(
            probabilities
        )

    check(
        all(
            str(s["measurement_complete"]) == "yes"
            for s in sites
        ),
        "Uno o più siti STEP28 sono incompleti"
    )

    return (
        p28,
        p27,
        sources,
        targets,
        sites,
        step22_index,
        pair_index,
        profile_data
    )


# ============================================================
# 20. ANALISI ACCESSIBILITA' RELATIVA
# ============================================================

def evaluate(state):

    (
        p28,
        p27,
        sources,
        targets,
        sites,
        step22,
        pairs,
        profiles
    ) = state

    output = []

    for row in sites:

        sid = row["srna_id"]
        tx = row["transcript_id"]

        key = (
            sid,
            tx
        )

        check(
            key in step22
            and key in pairs,
            f"Coppia candidata mancante: {key}"
        )

        original = step22[key]

        stat = pairs[key]

        start = int(
            row["target_start_1based"]
        )

        end = int(
            row["target_end_1based"]
        )

        length = (
            end - start + 1
        )

        # ====================================================
        # 21. VERIFICA IDENTITA' SITO
        # ====================================================

        check(
            1 <= start <= end <= len(targets[tx])
            and length in (24, 25),
            f"Coordinate STEP28 invalide: {key}"
        )

        check(
            (
                start,
                end
            ) == (

                int(
                    original["start_1based"]
                ),

                int(
                    original["end_1based"]
                )
            )
            and original["gene_id"]
            == row["gene_id"]
            == stat["gene_id"],
            "Coordinate o ID genici "
            f"non corrispondenti: {key}"
        )

        expected_sequence = (
            original["full_site_DNA_5to3"]
            .upper()
            .replace("T", "U")
        )

        check(
            targets[tx][start - 1:end]
            == expected_sequence,
            f"Sito non coincide con FASTA: {key}"
        )

        check(
            stat["STEP22_candidate"] == "yes",
            "STEP27 non considera candidata "
            f"la coppia {key}"
        )

        # ====================================================
        # 22. VERIFICA ENERGIA RNAplex
        # ====================================================

        check(
            abs(
                float(
                    row["RNAplex_MFE_kcal_mol"]
                )
                - float(
                    stat["energy_original_kcal_mol"]
                )
            ) < 1e-6,
            f"Energia STEP27/28 discordante: {key}"
        )

        # ====================================================
        # 23. RICALCOLO P10 DEL SITO
        # ====================================================

        prefix = profiles[tx]

        site_mean = region_mean(
            prefix,
            start,
            end
        )

        check(
            site_mean is not None
            and int(
                row["windows10_expected"]
            ) == length - LENGTH + 1
            and int(
                row["windows10_available"]
            ) == length - LENGTH + 1
            and abs(
                site_mean
                - float(
                    row["mean_P_unpaired_10nt"]
                )
            ) < 5e-8,
            "P10 ricalcolata differente "
            f"da STEP28: {key}"
        )

        # ====================================================
        # 24. BACKGROUND DELLO STESSO TRASCRITTO
        # ====================================================

        background = [

            region_mean(
                prefix,
                a,
                a + length - 1
            )

            for a in range(
                1,
                len(targets[tx]) - length + 2
            )
        ]

        background = [

            value
            for value in background
            if value is not None
        ]

        check(
            bool(background),
            f"Senza finestre di background: {key}"
        )

        # ====================================================
        # 25. PERCENTILE DESCRITTIVO
        # ====================================================

        fraction_le = sum(

            value <= site_mean + 1e-12

            for value in background

        ) / len(background)

        median_background = (
            statistics.median(background)
        )

        percentile = (
            100.0 * fraction_le
        )

        # ====================================================
        # 26. RISULTATO PER SITO
        # ====================================================

        output.append({

            "srna_id":
                sid,

            "transcript_id":
                tx,

            "gene_id":
                row["gene_id"],

            "site_start_1based":
                start,

            "site_end_1based":
                end,

            "site_length_nt":
                length,

            "P10_site_mean":
                f"{site_mean:.10g}",

            "P10_transcript_windows_median":
                f"{median_background:.10g}",

            "P10_site_minus_background_median":
                f"{site_mean - median_background:.10g}",

            "percentile_accessibilita_nel_trascritto":
                f"{percentile:.3f}",

            "number_background_windows":
                len(background),

            "RNAplex_MFE_kcal_mol":
                stat["energy_original_kcal_mol"],

            "RNAplex_empirical_p":
                stat["p_empirical_31"],

            "RNAplex_BH_q_74":
                stat["BH_q_74_exploratory"],

            "RNAseq_significance":
                "non_supportata_da_STEP21B",

            "biological_status":
                "exploratory_not_validated"
        })

    check(
        len(output) == 36,
        "Numero siti finali differente da 36"
    )

    # ========================================================
    # 27. ORDINAMENTO PER PERCENTILE
    # ========================================================

    output.sort(
        key=lambda row: (

            -float(
                row[
                    "percentile_accessibilita_nel_trascritto"
                ]
            ),

            row["srna_id"],

            row["transcript_id"]
        )
    )

    # ========================================================
    # 28. RAGGRUPPAMENTO PER GENE
    # ========================================================

    groups = defaultdict(list)

    for row in output:

        groups[
            (
                row["srna_id"],
                row["gene_id"]
            )
        ].append(row)

    gene_rows = []

    for (
        sid,
        gene
    ), items in sorted(
        groups.items()
    ):

        values = [

            float(
                item[
                    "percentile_accessibilita_nel_trascritto"
                ]
            )

            for item in items
        ]

        gene_rows.append({

            "srna_id":
                sid,

            "gene_id":
                gene,

            "n_trascritti":
                len(items),

            "mediana_percentile_accessibilita":
                f"{statistics.median(values):.3f}",

            "min_percentile":
                f"{min(values):.3f}",

            "max_percentile":
                f"{max(values):.3f}",

            "n_trascritti_percentile_ge_75":
                sum(
                    value >= 75
                    for value in values
                ),

            "biological_status":
                "exploratory_not_validated"
        })

    # ========================================================
    # 29. RIEPILOGO PER SMALL RNA
    # ========================================================

    summary = {}

    for sid in SRNAS:

        items = [

            row
            for row in output
            if row["srna_id"] == sid
        ]

        values = [

            float(
                row[
                    "percentile_accessibilita_nel_trascritto"
                ]
            )

            for row in items
        ]

        summary[sid] = {

            "sites":
                len(items),

            "unique_genes":
                len({
                    row["gene_id"]
                    for row in items
                }),

            "median_within_transcript_percentile":
                round(
                    statistics.median(values),
                    3
                ),

            "sites_percentile_ge_75":
                sum(
                    value >= 75
                    for value in values
                ),

            "sites_percentile_ge_90":
                sum(
                    value >= 90
                    for value in values
                ),

            "median_absolute_P10":
                round(
                    statistics.median(

                        float(
                            row["P10_site_mean"]
                        )

                        for row in items
                    ),
                    8
                )
        }

    return (
        output,
        gene_rows,
        summary
    )


# ============================================================
# 30. SCRITTURA TSV
# ============================================================

def write_tsv(path, rows):

    check(
        bool(rows),
        "Nessuna riga da scrivere"
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0]),
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(rows)


# ============================================================
# 31. PROGRAMMA PRINCIPALE
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

    parser.add_argument(
        "--step28-run",
        type=Path
    )

    parser.add_argument(
        "--step27-run",
        type=Path
    )

    args = parser.parse_args()

    print("=" * 70)

    print(
        "AI PRIMETTA V3.3 | STEP 29 "
        "- ACCESSIBILITA RELATIVA"
    )

    print("=" * 70)

    # ========================================================
    # 32. CARICAMENTO E CALCOLO
    # ========================================================

    state = load_source(
        args.step28_run,
        args.step27_run
    )

    (
        sites,
        genes,
        summary
    ) = evaluate(state)

    print(
        "[PASS] 37 profili originali RNAplfold "
        "e relativi SHA256"
    )

    print(
        "[PASS] 36 siti STEP22/28; "
        "74 coppie STEP27; "
        "confronto P10 ricostruito"
    )

    for sid, stats in sorted(
        summary.items()
    ):

        print(

            f"[RISULTATO] {sid}: "

            f"{json.dumps(stats, ensure_ascii=False)}"
        )

    print(
        "[LIMITE] Percentili descrittivi; "
        "finestre sovrapposte "
        "e target preselezionati"
    )

    # ========================================================
    # 33. AUDIT IN SOLA LETTURA
    # ========================================================

    if args.audit:

        print(
            "[OK] Audit in sola lettura; "
            "GUI e Docker invariati"
        )

        return

    # ========================================================
    # 34. NUOVA CARTELLA RISULTATI
    # ========================================================

    timestamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP29_ACCESSIBILITA_RELATIVA_"
        + timestamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    site_path = (
        folder
        / "ACCESSIBILITA_RELATIVA_36_SITI.tsv"
    )

    gene_path = (
        folder
        / "RIEPILOGO_ACCESSIBILITA_PER_GENE.tsv"
    )

    report_path = (
        folder
        / "RAPPORTO_STEP29.json"
    )

    # ========================================================
    # 35. ESPORTAZIONE RISULTATI
    # ========================================================

    write_tsv(
        site_path,
        sites
    )

    write_tsv(
        gene_path,
        genes
    )

    # ========================================================
    # 36. MANIFEST SCIENTIFICO
    # ========================================================

    report = {

        "step":
            29,

        "status":
            "DESCRIPTIVE_WITHIN_TRANSCRIPT_ACCESSIBILITY",

        "source_step28":
            str(
                state[0].relative_to(ROOT)
            ),

        "source_step27":
            str(
                state[1].relative_to(ROOT)
            ),

        "source_reports_sha256": {

            "step28":
                sha(state[0]),

            "step27":
                sha(state[1])
        },

        "length_unpaired_nt":
            10,

        "sites":
            36,

        "transcripts":
            37,

        "summary_by_srna":
            summary,

        "outputs": [

            str(
                site_path.relative_to(ROOT)
            ),

            str(
                gene_path.relative_to(ROOT)
            )
        ],

        "limitations": [

            "Confronto con finestre della stessa "
            "lunghezza sullo stesso trascritto.",

            "Finestre sovrapposte e dipendenti: "
            "i percentili non sono p-value.",

            "Il background include "
            "la finestra candidata.",

            "P10 e' la media delle probabilita' "
            "di segmenti non appaiati di 10 nt, "
            "non la probabilita' di apertura "
            "dell'intero sito.",

            "Trascritti preselezionati: "
            "nessuna inferenza non distorta "
            "sui target genomici.",

            "Nessun candidato con BH q<0.05 "
            "nello STEP27; nessuna "
            "dimostrazione biologica."
        ],

        "gui_modified":
            False,

        "docker_modified":
            False,

        "training_performed":
            False
    }

    report_path.write_text(

        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ) + "\n",

        encoding="utf-8"
    )

    # ========================================================
    # 37. OUTPUT FINALE
    # ========================================================

    print(
        "[SITI]",
        site_path
    )

    print(
        "[GENI]",
        gene_path
    )

    print(
        "[REPORT]",
        report_path
    )

    print(
        "[OK] Nessun profilo ricalcolato; "
        "GUI, Docker e modelli invariati"
    )


# ============================================================
# 38. AVVIO
# ============================================================

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
