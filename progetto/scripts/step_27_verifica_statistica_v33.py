#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 27 - VERIFICA STATISTICA COMPLETA

FUNZIONALITA':

--audit
    Controlla in sola lettura la provenienza
    dei dati, le 60 sequenze shuffle,
    le 2220 energie e le 74 coppie.

--run
    Esegue gli stessi controlli,
    applica la correzione Benjamini-Hochberg,
    aggrega i risultati per gene e
    genera TSV e JSON in una nuova cartella.

--step26-run PATH
    Facoltativo: permette di scegliere
    una specifica sessione STEP26.

NON ESEGUE:
    - Nuove analisi RNAplex
    - Training ML/GNN
    - Modifiche GUI
    - Modifiche Docker
    - Modifiche ai dati originali
"""

import argparse
import csv
import hashlib
import json
import statistics
import sys

from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
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

SRNA_LENGTHS = {
    "sRNA_00000037": 25,
    "sRNA_00000435": 24
}

N_SHUFFLES = 30


# ============================================================
# 2. VERIFICA PERCORSI
# ============================================================

def safe_file(value):

    p = Path(value)

    if not p.is_absolute():
        p = ROOT / p

    p = p.resolve()

    if (
        not p.is_file()
        or not p.is_relative_to(RUNS.resolve())
    ):

        raise ValueError(
            "File mancante oppure fuori dalla "
            f"cartella runs: {p}"
        )

    return p


# ============================================================
# 3. CHECKSUM SHA256
# ============================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for part in iter(
            lambda: f.read(1048576),
            b""
        ):

            h.update(part)

    return h.hexdigest()


# ============================================================
# 4. CONVERSIONE NUMERICA
# ============================================================

def decimal(value):

    try:

        v = Decimal(
            str(value)
            .strip()
            .replace(",", ".")
        )

        if not v.is_finite():
            raise ValueError("Numero non finito")

        return v

    except (ValueError, InvalidOperation) as exc:

        raise ValueError(
            f"Valore numerico non valido: {value!r}"
        ) from exc


# ============================================================
# 5. LETTURA TSV
# ============================================================

def table(path, required):

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        columns = reader.fieldnames or []

        if (
            len(columns) != len(set(columns))
            or not set(required).issubset(columns)
        ):

            missing = set(required) - set(columns)

            raise ValueError(
                f"Colonne mancanti o duplicate "
                f"in {path}: {missing}"
            )

        rows = list(reader)

    if not rows:
        raise ValueError(f"Tabella vuota: {path}")

    return rows


# ============================================================
# 6. LETTURA FASTA
# ============================================================

def fasta(path):

    result = {}

    label = None
    seq = []

    with path.open(
        "r",
        encoding="utf-8-sig"
    ) as f:

        for raw in f:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                if label is not None:

                    if (
                        not seq
                        or label in result
                    ):

                        raise ValueError(
                            "Record FASTA ripetuto/vuoto: "
                            + label
                        )

                    result[label] = (
                        "".join(seq)
                        .upper()
                        .replace("T", "U")
                    )

                label = line[1:].split()[0]
                seq = []

                if (
                    not label
                    or label in result
                ):

                    raise ValueError(
                        "ID FASTA vuoto/duplicato"
                    )

            else:

                if label is None:

                    raise ValueError(
                        "FASTA privo di intestazione"
                    )

                seq.append(line)

    if label is not None:

        if (
            not seq
            or label in result
        ):

            raise ValueError(
                "Ultimo record FASTA vuoto/duplicato"
            )

        result[label] = (
            "".join(seq)
            .upper()
            .replace("T", "U")
        )

    if not result:
        raise ValueError(f"FASTA vuoto: {path}")

    return result


# ============================================================
# 7. CORREZIONE BENJAMINI-HOCHBERG
# ============================================================

def bh(pvalues):

    """
    Calcola q-value con la procedura
    Benjamini-Hochberg.

    In questo progetto i risultati sono
    esclusivamente diagnostici/esplorativi,
    a causa della preselezione dei trascritti
    e della dipendenza tra osservazioni.
    """

    n = len(pvalues)

    order = sorted(
        range(n),
        key=lambda i: pvalues[i],
        reverse=True
    )

    out = [1.0] * n

    best = 1.0

    for offset, index in enumerate(order):

        rank = n - offset

        best = min(
            best,
            pvalues[index] * n / rank
        )

        out[index] = best

    return out


# ============================================================
# 8. MEDIA E MEDIANA
# ============================================================

def median_decimal(items):

    return statistics.median(items)


def mean_decimal(items):

    return (
        sum(items, Decimal(0))
        / len(items)
    )


# ============================================================
# 9. RECUPERO AUTOMATICO STEP26
# ============================================================

def read_provenance(step26_run=None):

    if step26_run is None:

        manifests = sorted(
            RUNS.glob(
                "STEP26_CONTROLLI_DINUC_RNAPLEX_*/"
                "RAPPORTO_STEP26.json"
            )
        )

        if not manifests:

            raise FileNotFoundError(
                "Nessun RAPPORTO_STEP26.json "
                "disponibile"
            )

        p26 = safe_file(
            manifests[-1]
        )

    else:

        requested = Path(step26_run)

        p26 = safe_file(
            requested / "RAPPORTO_STEP26.json"
        )

    r26 = json.loads(
        p26.read_text(
            encoding="utf-8"
        )
    )

    # ========================================================
    # 10. CONTROLLO MANIFEST STEP26
    # ========================================================

    if (
        r26.get("step") != 26
        or r26.get("status")
        != "EXPLORATORY_NEGATIVE_CONTROL"
        or r26.get("shuffle_per_srna")
        != N_SHUFFLES
        or r26.get("total_rnaplex_calls")
        != 2220
        or r26.get("n_original_target_pairs")
        != 74
        or r26.get("candidate_pairs")
        != 36
    ):

        raise ValueError(
            "Manifest STEP26 inatteso o incompleto"
        )

    # ========================================================
    # 11. RISALITA ALLE SORGENTI
    # ========================================================

    p25 = safe_file(
        r26["source_step25"]
    )

    r25 = json.loads(
        p25.read_text(
            encoding="utf-8"
        )
    )

    p24 = safe_file(
        r25["source_step24"]
    )

    r24 = json.loads(
        p24.read_text(
            encoding="utf-8"
        )
    )

    # Compatibilita' con STEP24: "24" oppure 24.

    if (
        str(r24.get("step")) != "24"
        or r24.get("RNAplex_pairs") != 74
    ):

        raise ValueError(
            "STEP24 non coerente"
        )

    p23 = safe_file(
        r24["source_step23b"]
    )

    r23 = json.loads(
        p23.read_text(
            encoding="utf-8"
        )
    )

    p22 = safe_file(
        r23["source_manifest"]
    )

    # ========================================================
    # 12. FILE SCIENTIFICI ORIGINALI
    # ========================================================

    sources = {

        "srna": safe_file(
            p22.parent
            / "SRNA_STANDARD_MAX25.fa"
        ),

        "targets": safe_file(
            p22.parent
            / "TRASCRITTI_B73_CANDIDATI.fa"
        ),

        "sites": safe_file(
            p22.parent
            / "SITI_E_RNASEQ_SNAPSHOT.tsv"
        ),

        "raw": safe_file(
            p23.parent
            / "OUTPUT_GREZZO_RNAPLEX.txt"
        ),

        "table": safe_file(
            p23.parent
            / "INTERAZIONI_RNAPLEX.tsv"
        )
    }

    # ========================================================
    # 13. CHECKSUM FILE SORGENTI
    # ========================================================

    for name, path in sources.items():

        expected_hash = (
            r26["source_sha256"][name]
        )

        if sha256(path) != expected_hash:

            raise ValueError(
                "Checksum origine STEP26 "
                f"non coerente: {name}"
            )

    # ========================================================
    # 14. FILE RISULTATI STEP26
    # ========================================================

    outputs = {}

    expected_names = {

        "SHUFFLE_DINUCLEOTIDI.tsv",

        "CONTROLLI_NEGATIVI_ENERGIE.tsv",

        "CONFRONTO_ORIGINALI_VS_DINUC.tsv"
    }

    for item in r26["output_files"]:

        p = safe_file(item)

        if (
            p.parent != p26.parent
            or p.name in outputs
        ):

            raise ValueError(
                "File risultati STEP26 duplicato "
                "o appartenente a un'altra sessione"
            )

        outputs[p.name] = p

    if set(outputs) != expected_names:

        raise ValueError(
            "Insieme file STEP26 incompleto"
        )

    result24 = safe_file(
        r24["output_tsv"]
    )

    return (
        p26,
        r26,
        result24,
        sources,
        outputs
    )


# ============================================================
# 15. RICOSTRUZIONE E VALIDAZIONE DEI DATI
# ============================================================

def validate_and_reconstruct(step26_run=None):

    (
        p26,
        r26,
        p24_tsv,
        sources,
        files
    ) = read_provenance(step26_run)

    srnas = fasta(
        sources["srna"]
    )

    targets = fasta(
        sources["targets"]
    )

    # ========================================================
    # 16. CONTROLLO NUMERO SEQUENZE
    # ========================================================

    lengths = {
        k: len(v)
        for k, v in srnas.items()
    }

    if (
        lengths != SRNA_LENGTHS
        or len(targets) != 37
    ):

        raise ValueError(
            "Attese 2 small RNA di 24/25 nt "
            "e 37 trascritti"
        )

    if any(
        set(v) - set("ACGU")
        for v in srnas.values()
    ):

        raise ValueError(
            "Small RNA con basi non canoniche"
        )

    # ========================================================
    # 17. VALIDAZIONE DEI 37 SITI STEP22
    # ========================================================

    sites = table(

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

    if (
        len(sites) != 37
        or len({
            x["transcript_id"]
            for x in sites
        }) != 37
    ):

        raise ValueError(
            "Attesi 37 siti su 37 trascritti distinti"
        )

    candidate_pairs = set()

    transcript_gene = {}

    for row in sites:

        tx = row["transcript_id"]

        sid = row["srna_id"]

        if tx not in targets:

            raise ValueError(
                "Trascritto STEP22 mancante: "
                + tx
            )

        a = int(
            row["start_1based"]
        )

        b = int(
            row["end_1based"]
        )

        fragment = (
            row["full_site_DNA_5to3"]
            .upper()
            .replace("T", "U")
        )

        if (
            not 1 <= a <= b <= len(targets[tx])
            or targets[tx][a - 1:b] != fragment
        ):

            raise ValueError(
                "Coordinate o sequenza STEP22 "
                f"non coerenti: {sid}, {tx}"
            )

        transcript_gene[tx] = (
            row["gene_id"].strip()
        )

        if sid in srnas:

            candidate_pairs.add(
                (sid, tx)
            )

    if len(candidate_pairs) != 36:

        raise ValueError(
            "Attese 36 coppie candidate STEP22"
        )

    # ========================================================
    # 18. ENERGIE ORIGINALI STEP24
    # ========================================================

    originals = {}

    rows24 = table(

        p24_tsv,

        [
            "srna_id",
            "transcript_id",
            "energy_kcal_mol"
        ]
    )

    for row in rows24:

        key = (
            row["srna_id"],
            row["transcript_id"]
        )

        if key in originals:

            raise ValueError(
                "Coppia STEP24 duplicata: "
                + repr(key)
            )

        originals[key] = decimal(
            row["energy_kcal_mol"]
        )

    all_pairs = {

        (sid, tx)

        for sid in srnas

        for tx in targets
    }

    if set(originals) != all_pairs:

        raise ValueError(
            "STEP24: matrice 74 energie incompleta"
        )

    # ========================================================
    # 19. VERIFICA DELLE 60 SEQUENZE SHUFFLE
    # ========================================================

    shuffle_records = table(

        files["SHUFFLE_DINUCLEOTIDI.tsv"],

        [
            "srna_id",
            "shuffle_index",
            "sequence_RNA_5to3"
        ]
    )

    shuffles = {}

    for row in shuffle_records:

        sid = row["srna_id"]

        index = int(
            row["shuffle_index"]
        )

        key = (
            sid,
            index
        )

        sequence = (
            row["sequence_RNA_5to3"]
            .upper()
            .replace("T", "U")
        )

        if (
            sid not in srnas
            or not 1 <= index <= N_SHUFFLES
            or key in shuffles
        ):

            raise ValueError(
                "Shuffle non valido/duplicato: "
                + repr(key)
            )

        if (
            sequence == srnas[sid]
            or len(sequence) != len(srnas[sid])
        ):

            raise ValueError(
                "Shuffle uguale all'originale "
                "o lunghezza errata: "
                + repr(key)
            )

        before = Counter(
            zip(
                srnas[sid],
                srnas[sid][1:]
            )
        )

        after = Counter(
            zip(
                sequence,
                sequence[1:]
            )
        )

        if before != after:

            raise ValueError(
                "Dinucleotidi non conservati: "
                + repr(key)
            )

        shuffles[key] = sequence

    if len(shuffles) != 60:

        raise ValueError(
            "Attese 60 sequenze negative"
        )

    for sid in srnas:

        distinct = {

            shuffles[(sid, i)]

            for i in range(1, 31)
        }

        if len(distinct) != N_SHUFFLES:

            raise ValueError(
                "Shuffle non tutti distinti: "
                + sid
            )

    # ========================================================
    # 20. LETTURA DELLE 2.220 ENERGIE NEGATIVE
    # ========================================================

    negative_rows = table(

        files["CONTROLLI_NEGATIVI_ENERGIE.tsv"],

        [
            "srna_id",
            "shuffle_index",
            "transcript_id",
            "gene_id",
            "step22_candidate",
            "shuffled_energy_kcal_mol",
            "original_energy_kcal_mol"
        ]
    )

    energies = {}

    for row in negative_rows:

        sid = row["srna_id"]

        tx = row["transcript_id"]

        index = int(
            row["shuffle_index"]
        )

        pair = (
            sid,
            tx
        )

        triple = (
            sid,
            index,
            tx
        )

        if (
            pair not in all_pairs
            or (sid, index) not in shuffles
            or triple in energies
        ):

            raise ValueError(
                "Controllo RNAplex duplicato "
                "o imprevisto: "
                + repr(triple)
            )

        # ====================================================
        # 21. CONFRONTO ENERGIA ORIGINALE
        # ====================================================

        original_value = decimal(
            row["original_energy_kcal_mol"]
        )

        if original_value != originals[pair]:

            raise ValueError(
                "Energia originale modificata: "
                + repr(pair)
            )

        # ====================================================
        # 22. CONTROLLO FLAG CANDIDATO
        # ====================================================

        expected_flag = (
            "yes"
            if pair in candidate_pairs
            else "no"
        )

        flag = (
            row["step22_candidate"]
            .strip()
            .lower()
        )

        if flag != expected_flag:

            raise ValueError(
                "Flag candidato discordante: "
                + repr(pair)
            )

        if (
            row["gene_id"].strip()
            != transcript_gene[tx]
        ):

            raise ValueError(
                "Gene discordante: "
                + repr(pair)
            )

        energies[triple] = decimal(
            row["shuffled_energy_kcal_mol"]
        )

    if (
        len(negative_rows) != 2220
        or len(energies) != 2220
    ):

        raise ValueError(
            "Matrice negativa incompleta: "
            "attese 2.220 righe uniche"
        )

    # ========================================================
    # 23. VERIFICA RIEPILOGO ORIGINALE STEP26
    # ========================================================

    summary_rows = table(

        files["CONFRONTO_ORIGINALI_VS_DINUC.tsv"],

        [
            "srna_id",
            "transcript_id",
            "original_energy_kcal_mol",
            "n_shuffle",
            "empirical_tail_fraction_unadjusted"
        ]
    )

    checked = set()

    for row in summary_rows:

        pair = (
            row["srna_id"],
            row["transcript_id"]
        )

        if (
            pair not in all_pairs
            or pair in checked
        ):

            raise ValueError(
                "Riepilogo STEP26 duplicato "
                "o imprevisto: "
                + repr(pair)
            )

        checked.add(pair)

        sid, tx = pair

        shuffled_values = [

            energies[(sid, i, tx)]

            for i in range(1, 31)
        ]

        count = sum(
            x <= originals[pair]
            for x in shuffled_values
        )

        p = (
            1 + count
        ) / 31

        value_p = float(
            row[
                "empirical_tail_fraction_unadjusted"
            ]
        )

        if (
            decimal(
                row["original_energy_kcal_mol"]
            ) != originals[pair]
            or int(row["n_shuffle"]) != 30
            or abs(value_p - p) > 0.0000011
        ):

            raise ValueError(
                "Riepilogo STEP26 "
                "non riproducibile: "
                + repr(pair)
            )

    if checked != all_pairs:

        raise ValueError(
            "Riepilogo STEP26: coppie mancanti"
        )

    print(
        "[PASS] SHA256 sorgenti STEP26, "
        "FASTA, STEP22 e STEP24"
    )

    print(
        "[PASS] 60 shuffle distinti "
        "con dinucleotidi conservati"
    )

    print(
        "[PASS] 2.220 energie, "
        "74 coppie, 36 candidate: "
        "congruenza confermata"
    )

    return (
        p26,
        r26,
        files,
        originals,
        energies,
        candidate_pairs,
        transcript_gene,
        srnas,
        targets
    )


# ============================================================
# 24. ANALISI STATISTICA
# ============================================================

def calculate(
    originals,
    negatives,
    candidate_pairs,
    transcript_gene,
    srnas,
    targets
):

    pairs = sorted(originals)

    pair_records = []

    pvalues = []

    # ========================================================
    # 25. P-VALUES EMPIRICI PER LE 74 COPPIE
    # ========================================================

    for sid, tx in pairs:

        null = [

            negatives[(sid, i, tx)]

            for i in range(1, 31)
        ]

        observed = originals[
            (sid, tx)
        ]

        better_or_equal = sum(
            value <= observed
            for value in null
        )

        p = (
            better_or_equal + 1
        ) / 31

        pvalues.append(p)

        pair_records.append({

            "srna_id":
                sid,

            "transcript_id":
                tx,

            "gene_id":
                transcript_gene[tx],

            "STEP22_candidate":
                (
                    "yes"
                    if (sid, tx) in candidate_pairs
                    else "no"
                ),

            "energy_original_kcal_mol":
                str(observed),

            "null_median_kcal_mol":
                str(
                    median_decimal(null)
                ),

            "null_min_kcal_mol":
                str(min(null)),

            "n_null_as_or_more_favorable":
                better_or_equal,

            "p_empirical_31":
                f"{p:.9f}"
        })

    # ========================================================
    # 26. CORREZIONE BH SU 74 COPPIE
    # ========================================================

    global_bh = bh(
        pvalues
    )

    for row, q in zip(
        pair_records,
        global_bh
    ):

        row[
            "BH_q_74_exploratory"
        ] = f"{q:.9f}"

    # ========================================================
    # 27. BH SU 37 COPPIE PER SMALL RNA
    # ========================================================

    for sid in sorted(srnas):

        indices = [

            i
            for i, row in enumerate(
                pair_records
            )

            if row["srna_id"] == sid
        ]

        subset_p = [

            pvalues[i]
            for i in indices
        ]

        adjusted = bh(
            subset_p
        )

        for i, q in zip(
            indices,
            adjusted
        ):

            pair_records[i][
                "BH_q_37_per_srna_exploratory"
            ] = f"{q:.9f}"

    # ========================================================
    # 28. AGGREGAZIONE PER GENE
    # ========================================================

    gene_records = []

    for sid in sorted(srnas):

        genes = defaultdict(list)

        for tx in sorted(targets):

            gene_id = (
                transcript_gene[tx]
                or f"UNANNOTATED::{tx}"
            )

            genes[gene_id].append(tx)

        for gene, txs in sorted(
            genes.items()
        ):

            obs = mean_decimal([

                originals[(sid, tx)]

                for tx in txs
            ])

            # Si conservano gli stessi 30 shuffle
            # per i trascritti dello stesso gene.

            null = [

                mean_decimal([

                    negatives[(sid, i, tx)]

                    for tx in txs

                ])

                for i in range(1, 31)
            ]

            n = sum(
                value <= obs
                for value in null
            )

            p = (
                1 + n
            ) / 31

            null_median = (
                median_decimal(null)
            )

            gene_records.append({

                "srna_id":
                    sid,

                "gene_id":
                    gene,

                "n_transcripts":
                    len(txs),

                "n_candidate_transcripts":
                    sum(
                        (sid, tx) in candidate_pairs
                        for tx in txs
                    ),

                "original_mean_kcal_mol":
                    str(obs),

                "null_median_of_gene_means_kcal_mol":
                    str(null_median),

                "delta_original_minus_null_median_kcal_mol":
                    str(obs - null_median),

                "n_null_as_or_more_favorable":
                    n,

                "p_empirical_31":
                    f"{p:.9f}"
            })

        # ====================================================
        # 29. BH AL LIVELLO DEI GENI
        # ====================================================

        indices = [

            i
            for i, row in enumerate(
                gene_records
            )

            if row["srna_id"] == sid
        ]

        values = [

            float(
                gene_records[i]["p_empirical_31"]
            )

            for i in indices
        ]

        for i, q in zip(
            indices,
            bh(values)
        ):

            gene_records[i][
                "BH_q_genes_per_srna_exploratory"
            ] = f"{q:.9f}"

    # ========================================================
    # 30. TEST AGGREGATI PER SMALL RNA
    # ========================================================

    group_records = []

    for sid in sorted(srnas):

        for group in (
            "candidate",
            "non_candidate",
            "all"
        ):

            chosen = [

                tx
                for tx in sorted(targets)

                if (
                    group == "all"
                    or (
                        ((sid, tx) in candidate_pairs)
                        == (group == "candidate")
                    )
                )
            ]

            if not chosen:
                continue

            # =================================================
            # 31. RAGGRUPPAMENTO DEI TRASCRITTI PER GENE
            # =================================================

            genes = defaultdict(list)

            for tx in chosen:

                gene_id = (
                    transcript_gene[tx]
                    or f"UNANNOTATED::{tx}"
                )

                genes[gene_id].append(tx)

            gene_sets = list(
                genes.values()
            )

            # =================================================
            # 32. MEDIA PRIMA PER GENE E POI TRA GENI
            # =================================================

            def aggregate(i):

                gene_means = []

                for txs in gene_sets:

                    if i == 0:

                        values = [

                            originals[(sid, tx)]

                            for tx in txs
                        ]

                    else:

                        values = [

                            negatives[(sid, i, tx)]

                            for tx in txs
                        ]

                    gene_means.append(
                        mean_decimal(values)
                    )

                return mean_decimal(
                    gene_means
                )

            observed = aggregate(0)

            null = [

                aggregate(i)

                for i in range(1, 31)
            ]

            count = sum(
                value <= observed
                for value in null
            )

            median_null = (
                median_decimal(null)
            )

            group_records.append({

                "srna_id":
                    sid,

                "group":
                    group,

                "n_transcripts":
                    len(chosen),

                "n_genes":
                    len(genes),

                "original_gene_weighted_energy":
                    str(observed),

                "null_median_gene_weighted_energy":
                    str(median_null),

                "delta_original_minus_null_median":
                    str(
                        observed - median_null
                    ),

                "n_null_as_or_more_favorable":
                    count,

                "p_empirical_31":
                    f"{(1 + count) / 31:.9f}"
            })

    # ========================================================
    # 33. RIEPILOGO STATISTICO
    # ========================================================

    summary = {

        "n_pairs":
            len(pair_records),

        "n_negative_energy_values":
            len(negatives),

        "n_candidates":
            len(candidate_pairs),

        "n_genes":
            len(
                set(
                    transcript_gene.values()
                )
            ),

        "empirical_p_min":
            1 / 31,

        "pair_nominal_p_below_0_05":
            sum(
                p < 0.05
                for p in pvalues
            ),

        "pair_BH_74_q_below_0_05":
            sum(
                q < 0.05
                for q in global_bh
            ),

        "pair_BH_74_q_below_0_10":
            sum(
                q < 0.10
                for q in global_bh
            ),

        "candidate_by_srna": {},

        "group_tests":
            group_records
    }

    # ========================================================
    # 34. STATISTICHE DELLE CANDIDATE
    # ========================================================

    for sid in sorted(srnas):

        records = [

            r
            for r in pair_records

            if (
                r["srna_id"] == sid
                and r["STEP22_candidate"] == "yes"
            )
        ]

        summary[
            "candidate_by_srna"
        ][sid] = {

            "candidate_pairs":
                len(records),

            "nominal_p_below_0_05":
                sum(
                    float(
                        r["p_empirical_31"]
                    ) < 0.05
                    for r in records
                ),

            "BH_74_q_below_0_05":
                sum(
                    float(
                        r["BH_q_74_exploratory"]
                    ) < 0.05
                    for r in records
                ),

            "min_BH_74_q":
                min(
                    float(
                        r["BH_q_74_exploratory"]
                    )
                    for r in records
                )
        }

    return (
        pair_records,
        gene_records,
        group_records,
        summary
    )


# ============================================================
# 35. SCRITTURA TSV
# ============================================================

def save_tsv(path, records):

    if not records:

        raise ValueError(
            "Nessun record per "
            + str(path)
        )

    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=list(records[0])
        )

        writer.writeheader()

        writer.writerows(
            records
        )


# ============================================================
# 36. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    modes = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    modes.add_argument(
        "--audit",
        action="store_true",
        help="Non scrive alcun file"
    )

    modes.add_argument(
        "--run",
        action="store_true",
        help="Esporta tre TSV e un JSON"
    )

    parser.add_argument(
        "--step26-run",
        type=Path,
        default=None,
        help="Cartella STEP26 specifica"
    )

    args = parser.parse_args()

    print("=" * 72)

    print(
        "AI PRIMETTA V3.3 | STEP 27 "
        "- VERIFICA STATISTICA COMPLETA"
    )

    print("=" * 72)

    # ========================================================
    # 37. ACQUISIZIONE E VERIFICA DEI DATI
    # ========================================================

    (
        p26,
        manifest,
        files,
        originals,
        negatives,
        candidates,
        genes,
        srnas,
        targets
    ) = validate_and_reconstruct(
        args.step26_run
    )

    # ========================================================
    # 38. CALCOLO STATISTICO
    # ========================================================

    (
        pair_rows,
        gene_rows,
        group_rows,
        summary
    ) = calculate(

        originals,

        negatives,

        candidates,

        genes,

        srnas,

        targets
    )

    # ========================================================
    # 39. OUTPUT TERMINALE
    # ========================================================

    print(
        "[STEP26]",
        p26
    )

    print(
        "[RISOLUZIONE] "
        "p empirico minimo = 1/31 = 0.032258 "
        "(30 shuffle)"
    )

    print(
        "[BH-74] p<0.05 =",
        summary[
            "pair_nominal_p_below_0_05"
        ],
        "| q_BH<0.05 =",
        summary[
            "pair_BH_74_q_below_0_05"
        ]
    )

    for sid, item in (
        summary[
            "candidate_by_srna"
        ].items()
    ):

        print(
            f"[RISULTATO] {sid}: "
            f"{json.dumps(item, ensure_ascii=False)}"
        )

    for row in group_rows:

        print(

            f"[GENI] {row['srna_id']} "

            f"| {row['group']} "

            f"| geni={row['n_genes']} "

            f"| trascritti={row['n_transcripts']} "

            f"| delta="
            f"{row['delta_original_minus_null_median']} "
            "kcal/mol "

            f"| p={row['p_empirical_31']}"
        )

    print(
        "[AVVERTENZA] BH e p empirici "
        "sono descrittivi: preselezione "
        "e dipendenze tra trascritti."
    )

    # ========================================================
    # 40. AUDIT IN SOLA LETTURA
    # ========================================================

    if args.audit:

        print(
            "[OK] Audit completo in sola lettura; "
            "GUI e Docker invariati"
        )

        return

    # ========================================================
    # 41. CREAZIONE NUOVA CARTELLA
    # ========================================================

    timestamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S%fZ"
    )

    outdir = RUNS / (
        "STEP27_STATISTICA_DINUC_"
        + timestamp
    )

    outdir.mkdir(
        parents=True,
        exist_ok=False
    )

    # ========================================================
    # 42. DEFINIZIONE FILE OUTPUT
    # ========================================================

    out_pairs = (
        outdir
        / "CONFRONTO_74_PVALUE_QVALUE.tsv"
    )

    out_genes = (
        outdir
        / "RIEPILOGO_PER_GENE.tsv"
    )

    out_groups = (
        outdir
        / "TEST_AGGREGATI_PER_SRNA.tsv"
    )

    out_report = (
        outdir
        / "RAPPORTO_STEP27.json"
    )

    # ========================================================
    # 43. ESPORTAZIONE TSV
    # ========================================================

    save_tsv(
        out_pairs,
        pair_rows
    )

    save_tsv(
        out_genes,
        gene_rows
    )

    save_tsv(
        out_groups,
        group_rows
    )

    # ========================================================
    # 44. MANIFEST SCIENTIFICO
    # ========================================================

    report = {

        "step":
            27,

        "status":
            "EXPLORATORY_STATISTICS_NOT_TARGET_VALIDATION",

        "source_step26":
            str(p26.relative_to(ROOT)),

        "source_step26_sha256":
            sha256(p26),

        "source_results_sha256": {

            name: sha256(path)

            for name, path in files.items()
        },

        "summary":
            summary,

        "outputs": [

            str(p.relative_to(ROOT))

            for p in (

                out_pairs,

                out_genes,

                out_groups
            )
        ],

        "limitations": [

            "Trascritti e siti preselezionati "
            "in base alla complementarita'.",

            "Inferenza statistica post-selezione.",

            "Gli stessi 30 shuffle per sRNA "
            "sono condivisi tra 37 trascritti.",

            "I test non rappresentano "
            "osservazioni indipendenti.",

            "BH per coppia e per gene "
            "e' un controllo descrittivo, "
            "non una validazione biologica.",

            "Risoluzione empirica dei p-value "
            "pari a 1/31.",

            "Il test aggregato utilizza "
            "un peso equivalente per gene.",

            "RNAplex non valuta accessibilita' "
            "e non dimostra targeting "
            "sperimentale."
        ],

        "gui_modified":
            False,

        "docker_modified":
            False,

        "training_performed":
            False
    }

    # ========================================================
    # 45. SALVATAGGIO RAPPORTO JSON
    # ========================================================

    out_report.write_text(

        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ) + "\n",

        encoding="utf-8"
    )

    # ========================================================
    # 46. RIEPILOGO FINALE
    # ========================================================

    print(
        "[COPPIE]",
        out_pairs
    )

    print(
        "[GENI]",
        out_genes
    )

    print(
        "[TEST]",
        out_groups
    )

    print(
        "[REPORT]",
        out_report
    )

    print(
        "[OK] STEP27 completato: "
        "nessun modello, GUI o dato originale "
        "modificato"
    )


# ============================================================
# 47. AVVIO
# ============================================================

if __name__ == "__main__":

    try:

        main()

    except Exception as error:

        print(

            f"[ERRORE] "
            f"{type(error).__name__}: {error}",

            file=sys.stderr
        )

        sys.exit(1)