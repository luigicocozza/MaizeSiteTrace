#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 | STEP 26

CONTROLLI NEGATIVI DINUCLEOTIDICI RNAplex

--audit
    Verifica i file originali, genera controlli
    sintetici e riproduce una energia RNAplex.

--run
    Esegue i controlli negativi contro
    gli stessi 37 trascritti B73.

Non modifica:
    GUI V3.2
    GUI V3.3
    Docker
    Modelli ML/GNN
    Dati scientifici originali
"""

import argparse
import csv
import hashlib
import json
import random
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal
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

EXPECTED = {
    "sRNA_00000435": 24,
    "sRNA_00000037": 25
}

HIT = re.compile(
    r"\S+\s+"
    r"\d+,\d+\s*:\s*\d+,\d+\s*"
    r"\(\s*([+-]?\d+(?:\.\d+)?)"
)


# ============================================================
# 2. VERIFICA PERCORSI
# ============================================================

def resolve_file(value):

    path = Path(value)

    if not path.is_absolute():
        path = ROOT / path

    path = path.resolve()

    if (
        not path.is_file()
        or not path.is_relative_to(
            RUNS.resolve()
        )
    ):

        raise RuntimeError(
            "File non disponibile "
            "o fuori dalla cartella runs: "
            + str(path)
        )

    return path


# ============================================================
# 3. SHA256
# ============================================================

def sha(path):

    h = hashlib.sha256()

    with path.open("rb") as file:

        for block in iter(
            lambda: file.read(1 << 20),
            b""
        ):

            h.update(block)

    return h.hexdigest()


# ============================================================
# 4. LETTURA FASTA
# ============================================================

def read_fasta(path):

    records = {}

    name = None
    chunks = []

    for raw in path.read_text(
        encoding="utf-8-sig"
    ).splitlines():

        line = raw.strip()

        if not line:
            continue

        if line.startswith(">"):

            if name is not None:

                if (
                    not chunks
                    or name in records
                ):

                    raise ValueError(
                        "Record FASTA non valido"
                    )

                records[name] = (
                    "".join(chunks)
                    .upper()
                    .replace("T", "U")
                )

            name = line[1:].split()[0]
            chunks = []

        else:

            if name is None:

                raise ValueError(
                    "FASTA senza intestazione"
                )

            chunks.append(line)

    if name is not None:

        if (
            not chunks
            or name in records
        ):

            raise ValueError(
                "Record FASTA finale non valido"
            )

        records[name] = (
            "".join(chunks)
            .upper()
            .replace("T", "U")
        )

    if (
        not records
        or any(
            set(seq) - set("ACGU")
            for seq in records.values()
        )
    ):

        raise ValueError(
            "FASTA vuoto o sequenze non canoniche"
        )

    return records


# ============================================================
# 5. LETTURA TSV
# ============================================================

def read_tsv(path):

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as file:

        reader = csv.DictReader(
            file,
            delimiter="\t"
        )

        columns = reader.fieldnames or []

        if (
            not columns
            or len(set(columns)) != len(columns)
        ):

            raise ValueError(
                "Intestazioni TSV non valide"
            )

        return list(reader)


# ============================================================
# 6. RECUPERO RISULTATI SCIENTIFICI ORIGINALI
# ============================================================

def load_sources():

    reports = sorted(
        RUNS.glob(
            "STEP25_ENERGIA_COPERTURA_*/"
            "RAPPORTO_STEP25.json"
        )
    )

    if not reports:

        raise FileNotFoundError(
            "RAPPORTO_STEP25.json non trovato"
        )

    p25 = resolve_file(
        reports[-1]
    )

    r25 = json.loads(
        p25.read_text(
            encoding="utf-8"
        )
    )

    if r25.get("step") != "25":

        raise ValueError(
            "Rapporto STEP25 non valido"
        )

    p24 = resolve_file(
        r25["source_step24"]
    )

    r24 = json.loads(
        p24.read_text(
            encoding="utf-8"
        )
    )

    if (
        r24.get("step") != "24"
        or r24.get("RNAplex_pairs") != 74
    ):

        raise ValueError(
            "Rapporto STEP24 non coerente"
        )

    p23 = resolve_file(
        r24["source_step23b"]
    )

    r23 = json.loads(
        p23.read_text(
            encoding="utf-8"
        )
    )

    manifest = resolve_file(
        r23["source_manifest"]
    )

    # ========================================================
    # 7. FILE DA VERIFICARE
    # ========================================================

    paths = {

        "srna": resolve_file(
            manifest.parent
            / "SRNA_STANDARD_MAX25.fa"
        ),

        "targets": resolve_file(
            manifest.parent
            / "TRASCRITTI_B73_CANDIDATI.fa"
        ),

        "sites": resolve_file(
            manifest.parent
            / "SITI_E_RNASEQ_SNAPSHOT.tsv"
        ),

        "raw": resolve_file(
            p23.parent
            / "OUTPUT_GREZZO_RNAPLEX.txt"
        ),

        "table": resolve_file(
            p23.parent
            / "INTERAZIONI_RNAPLEX.tsv"
        )
    }

    for key, path in paths.items():

        if sha(path) != (
            r25["source_sha256"][key]
        ):

            raise ValueError(
                "Checksum modificato "
                "rispetto a STEP25: "
                + key
            )

    srnas = read_fasta(
        paths["srna"]
    )

    targets = read_fasta(
        paths["targets"]
    )

    if (
        {
            key: len(seq)
            for key, seq in srnas.items()
        } != EXPECTED
        or len(targets) != 37
    ):

        raise ValueError(
            "Attesi 2 RNA e 37 trascritti"
        )

    # ========================================================
    # 8. IDENTIFICAZIONE DELLE COPPIE CANDIDATE
    # ========================================================

    sites = read_tsv(
        paths["sites"]
    )

    if (
        len(sites) != 37
        or len({
            row["transcript_id"]
            for row in sites
        }) != 37
    ):

        raise ValueError(
            "Siti o trascritti STEP22 non coerenti"
        )

    candidate_pairs = set()

    gene_mapping = {}

    for site in sites:

        sid = site["srna_id"]
        tx = site["transcript_id"]

        a = int(
            site["start_1based"]
        )

        b = int(
            site["end_1based"]
        )

        if tx not in targets:

            raise ValueError(
                "Trascritto non trovato"
            )

        expected_sequence = (
            site["full_site_DNA_5to3"]
            .replace("T", "U")
            .upper()
        )

        observed_sequence = (
            targets[tx][a - 1:b]
        )

        if (
            observed_sequence
            != expected_sequence
        ):

            raise ValueError(
                "Sito non corrispondente al FASTA"
            )

        gene_mapping[tx] = (
            site.get("gene_id", "")
        )

        if sid in srnas:

            candidate_pairs.add(
                (sid, tx)
            )

    if len(candidate_pairs) != 36:

        raise ValueError(
            "Attese 36 coppie candidate"
        )

    # ========================================================
    # 9. ENERGIE ORIGINALI STEP24
    # ========================================================

    original_rows = read_tsv(
        resolve_file(
            r24["output_tsv"]
        )
    )

    original = {}

    for row in original_rows:

        key = (
            row["srna_id"],
            row["transcript_id"]
        )

        if key in original:

            raise ValueError(
                "Coppia RNAplex duplicata"
            )

        value = Decimal(
            row["energy_kcal_mol"]
        )

        if not value.is_finite():

            raise ValueError(
                "Energia non finita"
            )

        original[key] = value

    all_pairs = {

        (sid, tx)
        for sid in srnas
        for tx in targets
    }

    if set(original) != all_pairs:

        raise ValueError(
            "Mancano alcune delle 74 energie"
        )

    executable = shutil.which(
        "RNAplex"
    )

    if executable is None:

        raise FileNotFoundError(
            "RNAplex non disponibile"
        )

    result = subprocess.run(
        [
            executable,
            "--version"
        ],
        capture_output=True,
        text=True,
        check=True,
        timeout=15
    )

    version = result.stdout.strip()

    return (
        p25,
        paths,
        srnas,
        targets,
        candidate_pairs,
        gene_mapping,
        original,
        executable,
        version
    )


# ============================================================
# 10. SHUFFLE CON CONSERVAZIONE DEI DINUCLEOTIDI
# ============================================================

def shuffle_dinuc(sequence, rng):

    """
    Rimescola gli archi di un grafo diretto
    dei dinucleotidi.

    Un percorso euleriano conserva esattamente
    il numero di occorrenze di ogni dinucleotide.
    """

    adjacency = {
        base: []
        for base in set(sequence)
    }

    for a, b in zip(
        sequence,
        sequence[1:]
    ):

        adjacency[a].append(b)

    for edges in adjacency.values():

        rng.shuffle(edges)

    stack = [
        sequence[0]
    ]

    path = []

    while stack:

        base = stack[-1]

        if adjacency[base]:

            stack.append(
                adjacency[base].pop()
            )

        else:

            path.append(
                stack.pop()
            )

    shuffled = "".join(
        reversed(path)
    )

    before = Counter(
        zip(
            sequence,
            sequence[1:]
        )
    )

    after = Counter(
        zip(
            shuffled,
            shuffled[1:]
        )
    )

    if (
        len(shuffled) != len(sequence)
        or before != after
    ):

        raise AssertionError(
            "Dinucleotidi non conservati"
        )

    return shuffled


# ============================================================
# 11. GENERAZIONE DI SEQUENZE NEGATIVE UNICHE
# ============================================================

def make_shuffles(
    srnas,
    count,
    seed
):

    rng = random.Random(
        seed
    )

    result = {}

    for sid, original in sorted(
        srnas.items()
    ):

        generated = set()

        attempts = 0

        while (
            len(generated) < count
            and attempts < 10000
        ):

            attempts += 1

            shuffled = shuffle_dinuc(
                original,
                rng
            )

            if shuffled != original:

                generated.add(
                    shuffled
                )

        if len(generated) != count:

            raise RuntimeError(
                "Numero insufficiente "
                "di controlli unici per "
                + sid
            )

        result[sid] = sorted(
            generated
        )

    return result


# ============================================================
# 12. ESECUZIONE RNAplex
# ============================================================

def run_rnaplex(
    executable,
    sid,
    query,
    tx,
    target,
    qpath,
    tpath
):

    qpath.write_text(
        f">{sid}\n{query}\n",
        encoding="utf-8"
    )

    tpath.write_text(
        f">{tx}\n{target}\n",
        encoding="utf-8"
    )

    proc = subprocess.run(
        [
            executable,
            "-q",
            str(qpath),
            "-t",
            str(tpath)
        ],
        capture_output=True,
        text=True,
        timeout=120,
        check=False
    )

    if proc.returncode != 0:

        raise RuntimeError(
            "RNAplex fallito: "
            + sid
            + "/"
            + tx
            + " "
            + proc.stderr[:350]
        )

    matches = HIT.findall(
        proc.stdout
    )

    if len(matches) != 1:

        raise RuntimeError(
            "Output RNAplex inatteso: "
            + sid
            + "/"
            + tx
            + " "
            + proc.stdout[:350]
        )

    return Decimal(
        matches[0]
    )


# ============================================================
# 13. AUDIT
# ============================================================

def audit(
    data,
    count,
    seed
):

    (
        p25,
        paths,
        srnas,
        targets,
        candidates,
        genes,
        original,
        executable,
        version
    ) = data

    print(
        "[VERSIONE]",
        version
    )

    print(
        "[PASS] SHA256 STEP25, "
        "2 small RNA, 37 trascritti, "
        "74 energie, 36 candidati"
    )

    shuffled = make_shuffles(
        srnas,
        max(3, count),
        seed
    )

    print(
        "[PASS] Shuffle: lunghezza "
        "e dinucleotidi conservati"
    )

    # ========================================================
    # 14. TEST RIPRODUCIBILITA' ENERGIA ORIGINALE
    # ========================================================

    sid = sorted(srnas)[0]
    tx = sorted(targets)[0]

    with tempfile.TemporaryDirectory(
        prefix="primetta26_audit_"
    ) as folder:

        qpath = Path(folder) / "query.fa"
        tpath = Path(folder) / "target.fa"

        computed = run_rnaplex(
            executable,
            sid,
            srnas[sid],
            tx,
            targets[tx],
            qpath,
            tpath
        )

    if computed != original[
        (sid, tx)
    ]:

        raise RuntimeError(
            "Energia originale "
            "non riproducibile: "
            + str(computed)
            + " vs "
            + str(original[(sid, tx)])
        )

    print(
        "[PASS] Energia originale "
        "RNAplex riprodotta:",
        computed,
        "kcal/mol"
    )

    print(
        "[CONTROLLO NEGATIVO]",
        count,
        "shuffle per RNA"
    )

    return shuffled


# ============================================================
# 15. ANALISI DEI CONTROLLI NEGATIVI
# ============================================================

def perform(
    data,
    shuffled,
    seed
):

    (
        p25,
        paths,
        srnas,
        targets,
        candidates,
        genes,
        original,
        executable,
        version
    ) = data

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    folder = RUNS / (
        "STEP26_CONTROLLI_DINUC_RNAPLEX_"
        + stamp
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    # ========================================================
    # 16. SALVATAGGIO SEQUENZE SHUFFLE
    # ========================================================

    shuffle_file = (
        folder
        / "SHUFFLE_DINUCLEOTIDI.tsv"
    )

    with shuffle_file.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as file:

        writer = csv.writer(
            file,
            delimiter="\t"
        )

        writer.writerow([
            "srna_id",
            "shuffle_index",
            "sequence_RNA_5to3"
        ])

        for sid, sequences in sorted(
            shuffled.items()
        ):

            for index, sequence in enumerate(
                sequences,
                1
            ):

                writer.writerow([
                    sid,
                    index,
                    sequence
                ])

    # ========================================================
    # 17. ESECUZIONE DI TUTTI I CONTROLLI
    # ========================================================

    result_file = (
        folder
        / "CONTROLLI_NEGATIVI_ENERGIE.tsv"
    )

    nulls = defaultdict(list)

    total = 0

    with tempfile.TemporaryDirectory(
        prefix="primetta26_"
    ) as temp:

        qpath = Path(temp) / "query.fa"
        tpath = Path(temp) / "target.fa"

        with result_file.open(
            "w",
            encoding="utf-8",
            newline=""
        ) as file:

            writer = csv.writer(
                file,
                delimiter="\t"
            )

            writer.writerow([
                "srna_id",
                "shuffle_index",
                "transcript_id",
                "gene_id",
                "step22_candidate",
                "shuffled_energy_kcal_mol",
                "original_energy_kcal_mol"
            ])

            for sid in sorted(srnas):

                for index, sequence in enumerate(
                    shuffled[sid],
                    1
                ):

                    for tx, target in sorted(
                        targets.items()
                    ):

                        energy = run_rnaplex(
                            executable,
                            sid,
                            sequence,
                            tx,
                            target,
                            qpath,
                            tpath
                        )

                        nulls[
                            (sid, tx)
                        ].append(
                            energy
                        )

                        writer.writerow([
                            sid,
                            index,
                            tx,
                            genes.get(tx, ""),
                            (
                                "yes"
                                if (sid, tx) in candidates
                                else "no"
                            ),
                            str(energy),
                            str(original[(sid, tx)])
                        ])

                        total += 1

                    file.flush()

                    print(
                        f"[CONTROLLO] {sid} "
                        f"shuffle {index}/"
                        f"{len(shuffled[sid])} "
                        f"| confronti={total}",
                        flush=True
                    )

    # ========================================================
    # 18. CONFRONTO ORIGINALI VS CONTROLLI
    # ========================================================

    details = []

    statistics_by_srna = {}

    for sid in sorted(srnas):

        candidate_details = []
        other_details = []

        for tx in sorted(targets):

            null_values = nulls[
                (sid, tx)
            ]

            if len(null_values) != len(
                shuffled[sid]
            ):

                raise AssertionError(
                    "Numero repliche non corretto"
                )

            observed = original[
                (sid, tx)
            ]

            # Energia minore/più negativa:
            # interazione più favorevole.

            null_at_least_as_good = sum(
                x <= observed
                for x in null_values
            )

            empirical_fraction = (
                1 + null_at_least_as_good
            ) / (
                len(null_values) + 1
            )

            original_more_favorable = sum(
                observed < x
                for x in null_values
            )

            entry = {

                "srna_id":
                    sid,

                "transcript_id":
                    tx,

                "gene_id":
                    genes.get(tx, ""),

                "step22_candidate":
                    (
                        "yes"
                        if (sid, tx) in candidates
                        else "no"
                    ),

                "original_energy_kcal_mol":
                    str(observed),

                "null_median_energy_kcal_mol":
                    str(
                        statistics.median(
                            null_values
                        )
                    ),

                "null_min_energy_kcal_mol":
                    str(min(null_values)),

                "null_max_energy_kcal_mol":
                    str(max(null_values)),

                "original_more_negative_than_null_count":
                    original_more_favorable,

                "n_shuffle":
                    len(null_values),

                "empirical_tail_fraction_unadjusted":
                    f"{empirical_fraction:.6f}"
            }

            details.append(
                entry
            )

            if (
                sid,
                tx
            ) in candidates:

                candidate_details.append(
                    entry
                )

            else:

                other_details.append(
                    entry
                )

        # ====================================================
        # 19. STATISTICHE SEPARATE
        # ====================================================

        def describe(items):

            values = [
                float(
                    item[
                        "empirical_tail_fraction_unadjusted"
                    ]
                )
                for item in items
            ]

            better = sum(

                Decimal(
                    item[
                        "original_energy_kcal_mol"
                    ]
                )

                <

                Decimal(
                    item[
                        "null_median_energy_kcal_mol"
                    ]
                )

                for item in items
            )

            return {

                "pairs":
                    len(items),

                "original_better_than_null_median":
                    better,

                "median_unadjusted_tail_fraction":
                    (
                        statistics.median(values)
                        if values
                        else None
                    ),

                "min_unadjusted_tail_fraction":
                    (
                        min(values)
                        if values
                        else None
                    )
            }

        statistics_by_srna[sid] = {

            "candidate_pairs":
                describe(candidate_details),

            "other_pairs":
                describe(other_details)
        }

        print(
            f"[RISULTATO] {sid}: "
            f"{statistics_by_srna[sid]}"
        )

    # ========================================================
    # 20. TSV RIASSUNTIVO
    # ========================================================

    comparison_file = (
        folder
        / "CONFRONTO_ORIGINALI_VS_DINUC.tsv"
    )

    with comparison_file.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as file:

        writer = csv.DictWriter(
            file,
            delimiter="\t",
            fieldnames=list(
                details[0]
            )
        )

        writer.writeheader()

        writer.writerows(
            details
        )

    # ========================================================
    # 21. RAPPORTO SCIENTIFICO
    # ========================================================

    report = {

        "step": 26,

        "status":
            "EXPLORATORY_NEGATIVE_CONTROL",

        "source_step25":
            str(p25.relative_to(ROOT)),

        "source_sha256": {

            key: sha(path)
            for key, path in paths.items()
        },

        "software":
            version,

        "random_seed":
            seed,

        "dinucleotide_preserved":
            True,

        "shuffle_per_srna":
            len(
                next(
                    iter(
                        shuffled.values()
                    )
                )
            ),

        "total_rnaplex_calls":
            total,

        "n_original_target_pairs":
            len(original),

        "candidate_pairs":
            len(candidates),

        "by_srna":
            statistics_by_srna,

        "output_files": [

            str(path.relative_to(ROOT))

            for path in (
                result_file,
                shuffle_file,
                comparison_file
            )
        ],

        "limitations": [

            "Trascritti preselezionati "
            "per complementarita'.",

            "Shuffle con composizione "
            "dinucleotidica conservata.",

            "Frazioni empiriche esplorative "
            "non corrette per test multipli.",

            "Trascritti dello stesso gene "
            "non necessariamente indipendenti.",

            "RNAplex senza accessibilita'.",

            "Nessuna validazione biologica "
            "o sperimentale."
        ],

        "gui_modified":
            False,

        "docker_modified":
            False,

        "training_performed":
            False
    }

    report_file = (
        folder
        / "RAPPORTO_STEP26.json"
    )

    report_file.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ) + "\n",
        encoding="utf-8"
    )

    print(
        "[RISULTATI]",
        comparison_file
    )

    print(
        "[REPORT]",
        report_file
    )

    print(
        "[OK] File sorgente, GUI "
        "e Docker invariati"
    )


# ============================================================
# 22. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = parser.add_mutually_exclusive_group(
        required=True
    )

    mode.add_argument(
        "--audit",
        action="store_true"
    )

    mode.add_argument(
        "--run",
        action="store_true"
    )

    parser.add_argument(
        "--shuffles",
        type=int,
        default=30
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20261023
    )

    args = parser.parse_args()

    if not (
        3 <= args.shuffles <= 100
    ):

        parser.error(
            "--shuffles deve essere "
            "compreso tra 3 e 100"
        )

    print("=" * 69)

    print(
        "AI PRIMETTA V3.3 | STEP 26"
    )

    print(
        "CONTROLLI NEGATIVI DINUCLEOTIDICI"
    )

    print("=" * 69)

    data = load_sources()

    shuffled = audit(
        data,
        args.shuffles,
        args.seed
    )

    if args.audit:

        print(
            "[OK] Audit terminato; "
            "nessun file del progetto generato"
        )

    else:

        perform(
            data,
            shuffled,
            args.seed
        )


# ============================================================
# 23. AVVIO
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