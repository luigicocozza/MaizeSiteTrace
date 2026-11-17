
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 45

SUDDIVISIONE DEGLI ARCHI CON PROTEZIONE
DAL LEAKAGE DIRETTO FRA COPPIE DI NODI.

--audit:
    Controllo del grafo e simulazione split.
    Nessun file generato.

--run:
    Crea un nuovo TSV con train/val/test
    e un rapporto scientifico JSON.

ATTENZIONE:
- Nessun training.
- Nessun caricamento dei modelli.
- Nessuna validazione biologica.
- GUI, Docker e checkpoint invariati.
"""

import argparse
import collections
import csv
import hashlib
import json
import random
import sys

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

NODES = (
    ROOT
    / "data/processed"
    / "graph_nodes_v3_prjna419326.csv"
)

EDGES = (
    ROOT
    / "data/processed"
    / "graph_edges_v3_prjna419326.csv"
)

SEED = 42


# ============================================================
# 2. CONTROLLO CONDIZIONI
# ============================================================

def require(condition, message):

    if not condition:
        raise ValueError(message)


# ============================================================
# 3. CHECKSUM
# ============================================================

def sha256(path):

    digest = hashlib.sha256()

    with path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(1024 * 1024),
            b""
        ):

            digest.update(block)

    return digest.hexdigest()


# ============================================================
# 4. LETTURA CSV
# ============================================================

def read_csv(path, mandatory):

    require(
        path.is_file(),
        f"File assente: {path}"
    )

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(handle)

        fields = reader.fieldnames or []

        require(
            len(fields) == len(set(fields)),
            f"Colonne ripetute: {path.name}"
        )

        require(
            set(mandatory).issubset(fields),
            f"Colonne mancanti: {path.name}"
        )

        rows = list(reader)

    require(
        bool(rows),
        f"File senza dati: {path.name}"
    )

    return fields, rows


# ============================================================
# 5. IDENTIFICATIVO DELLA COPPIA
# ============================================================

def canonical_pair(row):

    return tuple(sorted((
        row["source_node_id"],
        row["target_node_id"]
    )))


# ============================================================
# 6. SPLIT RIPRODUCIBILE DELLE COPPIE
# ============================================================

def construct_split(groups):

    by_type = collections.defaultdict(list)

    for pair, rows in groups.items():

        types = collections.Counter(
            row["edge_type"]
            for row in rows
        )

        primary = sorted(
            types.items(),
            key=lambda item: (
                -item[1],
                item[0]
            )
        )[0][0]

        by_type[primary].append(pair)

    rng = random.Random(SEED)

    assignments = {}

    for edge_type in sorted(by_type):

        pairs = sorted(
            by_type[edge_type]
        )

        rng.shuffle(pairs)

        size = len(pairs)

        # Tipologie con meno di 3 coppie:
        # non possiamo creare tutti e tre
        # gli insiemi senza svuotarli.

        if size < 3:

            for pair in pairs:
                assignments[pair] = "train"

            continue

        n_val = max(
            1,
            round(0.1 * size)
        )

        n_test = max(
            1,
            round(0.1 * size)
        )

        if n_val + n_test >= size:
            n_val, n_test = 1, 1

        for i, pair in enumerate(pairs):

            if i < n_val:
                label = "val"

            elif i < n_val + n_test:
                label = "test"

            else:
                label = "train"

            assignments[pair] = label

    require(
        len(assignments) == len(groups),
        "Assegnazione coppie incompleta"
    )

    return assignments


# ============================================================
# 7. ANALISI DEL GRAFO
# ============================================================

def compute():

    _, node_rows = read_csv(
        NODES,
        ["node_id", "node_type"]
    )

    fields, edge_rows = read_csv(
        EDGES,
        [
            "edge_id",
            "source_node_id",
            "target_node_id",
            "source_node_type",
            "target_node_type",
            "edge_type",
            "is_direct_biological_evidence",
            "is_proxy_evidence"
        ]
    )

    nodes = {
        row["node_id"]: row
        for row in node_rows
    }

    require(
        len(nodes) == len(node_rows),
        "Identificativi dei nodi duplicati"
    )

    seen_ids = set()

    groups = collections.defaultdict(list)

    type_profiles = collections.defaultdict(
        lambda: {
            "source": set(),
            "target": set(),
            "directed": set(),
            "edges": 0
        }
    )

    flags = collections.Counter()

    # ========================================================
    # 8. VERIFICA DEGLI ARCHI
    # ========================================================

    for row in edge_rows:

        edge_id = row["edge_id"]

        require(
            edge_id
            and edge_id not in seen_ids,
            f"ID arco vuoto o duplicato: {edge_id}"
        )

        seen_ids.add(edge_id)

        source = row["source_node_id"]
        target = row["target_node_id"]

        require(
            source in nodes
            and target in nodes,
            f"Nodo mancante: {edge_id}"
        )

        require(
            row["source_node_type"]
            == nodes[source]["node_type"],
            f"Tipo nodo sorgente incoerente: {edge_id}"
        )

        declared_target_type = row["target_node_type"].strip()
        actual_target_type = nodes[target]["node_type"].strip()

        # In V3 il costruttore ha utilizzato "mixed" come
        # segnaposto per cinque specifiche relazioni.
        # Il tipo effettivo viene verificato sul nodo.

        allowed_mixed_targets = {
            "experimental_sequence_similarity": "ncRNA",
            "experimental_drought_response_similarity": "ncRNA",
            "belongs_to_species": "species",
            "has_rna_type": "rna_type",
            "found_candidate_expression_dataset":
                "public_expression_dataset",
        }

        target_type_valid = (
            declared_target_type == actual_target_type
            or (
                declared_target_type == "mixed"
                and allowed_mixed_targets.get(row["edge_type"])
                == actual_target_type
            )
        )

        require(
            target_type_valid,
            (
                f"Tipo nodo destinazione incoerente: {edge_id} | "
                f"arco={row['edge_type']} | "
                f"dichiarato={declared_target_type} | "
                f"reale={actual_target_type}"
            )
        )

        direct = (
            row["is_direct_biological_evidence"]
            .strip()
            .lower()
        )

        proxy = (
            row["is_proxy_evidence"]
            .strip()
            .lower()
        )

        require(
            direct in ("0", "1")
            and proxy in ("0", "1"),
            f"Flag evidenza non validi: {edge_id}"
        )

        flags[(direct, proxy)] += 1

        # Tutti gli archi tra due nodi
        # condividono la stessa coppia,
        # anche se cambiano verso o tipo.

        groups[
            canonical_pair(row)
        ].append(row)

        typ = (
            row["edge_type"],
            row["source_node_type"],
            row["target_node_type"]
        )

        entry = type_profiles[typ]

        entry["source"].add(source)
        entry["target"].add(target)

        entry["directed"].add(
            (source, target)
        )

        entry["edges"] += 1

    # ========================================================
    # 9. ASSEGNAZIONE TRAIN / VAL / TEST
    # ========================================================

    assignments = construct_split(
        groups
    )

    # ========================================================
    # 10. RILEVAMENTO ARCHI RECIPROCI
    # ========================================================

    reverse_pairs = 0

    for rows in groups.values():

        directed = {
            (
                row["source_node_id"],
                row["target_node_id"]
            )
            for row in rows
        }

        if any(
            source != target
            and (target, source) in directed
            for source, target in directed
        ):
            reverse_pairs += 1

    # ========================================================
    # 11. POSSIBILITA' DI CAMPIONAMENTO NEGATIVO
    # ========================================================

    relation_profiles = []

    for (
        name,
        source_type,
        target_type
    ), stats in sorted(type_profiles.items()):

        possible = (
            len(stats["source"])
            * len(stats["target"])
        )

        missing = (
            possible
            - len(stats["directed"])
        )

        relation_profiles.append({

            "edge_type":
                name,

            "source_type":
                source_type,

            "target_type":
                target_type,

            "rows":
                stats["edges"],

            "distinct_source_nodes":
                len(stats["source"]),

            "distinct_target_nodes":
                len(stats["target"]),

            "distinct_directed_pairs":
                len(stats["directed"]),

            "missing_pair_candidates_among_observed_endpoints":
                missing,

            "note":
                "A missing edge is not a "
                "validated biological negative."
        })

    split_count = collections.Counter(

        assignments[
            canonical_pair(row)
        ]

        for row in edge_rows
    )

    split_pairs = collections.Counter(
        assignments.values()
    )

    require(
        sum(split_count.values())
        == len(edge_rows),
        "Conteggio archi incompleto"
    )

    # ========================================================
    # 12. RAPPORTO SCIENTIFICO
    # ========================================================

    result = {

        "step": 45,

        "purpose":
            "STRICT_PAIR_SPLIT_NOT_BIOLOGICAL_VALIDATION",

        "node_count":
            len(nodes),

        "edge_count":
            len(edge_rows),

        "unordered_pairs":
            len(groups),

        "pairs_with_both_directions":
            reverse_pairs,

        "evidence_flags": {

            f"direct_{x}_proxy_{y}": n

            for (x, y), n in sorted(
                flags.items()
            )
        },

        "split_pair_counts":
            dict(split_pairs),

        "split_edge_counts":
            dict(split_count),

        "relations":
            relation_profiles,

        "seed":
            SEED,

        "input_nodes_sha256":
            sha256(NODES),

        "input_edges_sha256":
            sha256(EDGES),

        "warnings": [

            "All splits are grouped by "
            "unordered node pair.",

            "The existing V3 checkpoint "
            "cannot evaluate this holdout: "
            "it was trained on the "
            "original graph.",

            "Metadata and proxy links "
            "are not validated regulatory "
            "interactions.",

            "A missing edge is not a "
            "confirmed biological negative.",

            "This is a transductive "
            "edge split, not validation "
            "on entirely unseen nodes.",

            "No training, inference "
            "or GUI/Docker modifications."
        ]
    }

    return (
        fields,
        edge_rows,
        assignments,
        result
    )


# ============================================================
# 13. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = (
        parser.add_mutually_exclusive_group(
            required=True
        )
    )

    mode.add_argument(
        "--audit",
        action="store_true"
    )

    mode.add_argument(
        "--run",
        action="store_true"
    )

    args = parser.parse_args()

    (
        fields,
        rows,
        assignments,
        report
    ) = compute()

    print("=" * 69)

    print(
        "AI PRIMETTA V3.3 | STEP45"
    )

    print(
        "SPLIT ANTILEAKAGE E TIPOLOGIE ARCHI"
    )

    print("=" * 69)

    print(
        "[GRAFO]",
        f"nodi={report['node_count']}",
        f"archi={report['edge_count']}",
        f"coppie_non_ordinate="
        f"{report['unordered_pairs']}"
    )

    print(
        "[RECIPROCI] Coppie con archi "
        "in entrambi i versi=",
        report["pairs_with_both_directions"]
    )

    print(
        "[SPLIT] Archi=",
        report["split_edge_counts"],
        "| Coppie=",
        report["split_pair_counts"]
    )

    print(
        "[EVIDENZA]",
        report["evidence_flags"]
    )

    # ========================================================
    # 14. RIEPILOGO PER TIPO DI RELAZIONE
    # ========================================================

    for relation in report["relations"]:

        print(

            "[TIPO]",

            relation["edge_type"],

            "("
            + relation["source_type"]
            + " -> "
            + relation["target_type"]
            + ")",

            "archi=",
            relation["rows"],

            "sorgenti=",
            relation["distinct_source_nodes"],

            "destinazioni=",
            relation["distinct_target_nodes"],

            "coppie_mancanti_potenziali=",
            relation[
                "missing_pair_candidates_among_observed_endpoints"
            ]
        )

    print(
        "[VINCOLO] Coppie reciproche e "
        "tipi multipli nello stesso split"
    )

    print(
        "[VINCOLO] Checkpoint V3 esistente "
        "non idoneo a valutare "
        "questo nuovo holdout"
    )

    print(
        "[VINCOLO] Archi assenti non sono "
        "negativi biologici confermati"
    )

    # ========================================================
    # 15. SALVATAGGIO OPZIONALE
    # ========================================================

    if args.run:

        stamp = datetime.now(
            timezone.utc
        ).strftime(
            "%Y%m%dT%H%M%S%fZ"
        )

        outdir = RUNS / (
            "STEP45_SPLIT_STRICT_"
            + stamp
        )

        outdir.mkdir(
            parents=True,
            exist_ok=False
        )

        table_path = (
            outdir
            / "GRAPH_EDGES_V3_SPLIT_GROUPED.tsv"
        )

        with table_path.open(
            "w",
            encoding="utf-8",
            newline=""
        ) as f:

            writer = csv.DictWriter(

                f,

                fieldnames=fields + [
                    "heldout_split",
                    "unordered_pair_key"
                ],

                delimiter="\t"
            )

            writer.writeheader()

            for row in rows:

                pair = canonical_pair(row)

                writer.writerow({

                    **row,

                    "heldout_split":
                        assignments[pair],

                    "unordered_pair_key":
                        "||".join(pair)
                })

        report["output_tsv"] = str(
            table_path.relative_to(ROOT)
        )

        report["output_tsv_sha256"] = (
            sha256(table_path)
        )

        json_path = (
            outdir
            / "RAPPORTO_STEP45.json"
        )

        json_path.write_text(

            json.dumps(
                report,
                indent=2,
                ensure_ascii=False
            ) + "\n",

            encoding="utf-8"
        )

        print(
            "[TSV]",
            table_path
        )

        print(
            "[REPORT]",
            json_path
        )

    print(
        "[OK] Nessun training; "
        "GUI, Docker e modelli invariati"
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

        raise SystemExit(1)
