
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 | STEP 46

GAE sperimentale con holdout dello STEP45.

Obiettivo:
    Ricostruire relazioni computazionali
    di somiglianza tra small RNA.

Non misura:
    Validazione di interazioni regolatorie
    o predizione biologica sperimentale.

--audit: controlli in sola lettura.
--run: addestra un nuovo modello isolato.

Non modifica GUI, Docker, checkpoint precedenti.
"""

import argparse
import collections
import copy
import csv
import hashlib
import itertools
import json
import math
import os
import random
import sys

from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent

RUNS = ROOT / "target_prediction_v33_vscode/runs"

NODES = ROOT / (
    "data/processed/"
    "graph_nodes_v3_prjna419326.csv"
)

EDGES = ROOT / (
    "data/processed/"
    "graph_edges_v3_prjna419326.csv"
)

TARGET_TYPES = {
    "experimental_sequence_similarity",
    "experimental_drought_response_similarity",
}

SEED = 42
EPOCHS = 200
PATIENCE = 35
HIDDEN = 64
LATENT = 16


def fail(msg):
    raise RuntimeError(msg)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(
            lambda: f.read(1048576), b""
        ):
            h.update(block)
    return h.hexdigest()


def read_table(path, sep):
    with path.open(
        encoding="utf-8-sig", newline=""
    ) as f:
        return list(csv.DictReader(f, delimiter=sep))


def pair_key(a, b):
    return tuple(sorted((a, b)))


def prepare():

    candidates = sorted(RUNS.glob(
        "STEP45_SPLIT_STRICT_*/"
        "GRAPH_EDGES_V3_SPLIT_GROUPED.tsv"
    ))

    if not candidates:
        fail("Split STEP45 non trovato")

    split_file = candidates[-1]

    nodes = read_table(NODES, ",")
    source_edges = read_table(EDGES, ",")
    edges = read_table(split_file, "\t")

    if len(nodes) != 884:
        fail("Numero nodi inatteso")

    if len(edges) != 4782:
        fail("Numero archi inatteso")

    ids = [r["node_id"] for r in nodes]

    if len(set(ids)) != len(ids):
        fail("Nodi duplicati")

    node_index = {
        node_id: i
        for i, node_id in enumerate(ids)
    }

    originals = {
        r["edge_id"]: r
        for r in source_edges
    }

    if len(originals) != len(source_edges):
        fail("ID archi originali duplicati")

    if len(originals) != len(edges):
        fail("Inventario archi differente")

    grouped = {}

    splits = {"train", "val", "test"}

    for r in edges:

        eid = r["edge_id"]

        if eid not in originals:
            fail("ID arco estraneo allo STEP45")

        old = originals[eid]

        for col in (
            "source_node_id",
            "target_node_id",
            "edge_type"
        ):
            if r[col] != old[col]:
                fail(
                    "Arco alterato rispetto "
                    "all'originale: " + eid
                )

        a = r["source_node_id"]
        b = r["target_node_id"]

        if a not in node_index or b not in node_index:
            fail("Arco con nodo sconosciuto")

        s = r["heldout_split"]

        if s not in splits:
            fail("Split sconosciuto")

        key = pair_key(a, b)

        if key in grouped and grouped[key] != s:
            fail("LEAKAGE: stessa coppia in split diversi")

        grouped[key] = s

    observed = {
        pair_key(
            r["source_node_id"],
            r["target_node_id"]
        )
        for r in edges
    }

    by_split = {
        "train": set(),
        "val": set(),
        "test": set()
    }

    candidate_nodes = set()

    for r in edges:

        if r["edge_type"] not in TARGET_TYPES:
            continue

        a = r["source_node_id"]
        b = r["target_node_id"]

        if a == b:
            continue

        if (
            next(
                n for n in nodes
                if n["node_id"] == a
            )["node_type"] != "ncRNA"
            or next(
                n for n in nodes
                if n["node_id"] == b
            )["node_type"] != "ncRNA"
        ):
            fail("Relazione di somiglianza non ncRNA-ncRNA")

        candidate_nodes.update((a, b))

        by_split[r["heldout_split"]].add(
            pair_key(a, b)
        )

    for name in ("train", "val", "test"):
        if not by_split[name]:
            fail("Nessuna coppia positiva in " + name)

    # Coppie non osservate nel grafo:
    # negative solo ai fini del benchmark
    # strutturale, NON negative biologiche.

    unseen = [
        (a, b)
        for a, b in itertools.combinations(
            sorted(candidate_nodes), 2
        )
        if (a, b) not in observed
    ]

    rng = random.Random(SEED)
    rng.shuffle(unseen)

    nt = len(by_split["train"])
    nv = len(by_split["val"])
    ne = len(by_split["test"])

    required = 3 * nt + nv + ne

    if len(unseen) < required:
        fail("Coppie non osservate insufficienti")

    negatives = {
        "train_pool": unseen[:3 * nt],
        "val": unseen[3 * nt:3 * nt + nv],
        "test": unseen[
            3 * nt + nv:required
        ]
    }

    # Caratteristiche di base indipendenti
    # dalla matrice degli archi.
    # Non riutilizziamo gli embedding V3.

    node_types = sorted({
        r["node_type"] for r in nodes
    })

    rna_types = sorted({
        r.get("rna_type", "").strip()
        for r in nodes
        if r.get("rna_type", "").strip()
    })

    features = []

    for r in nodes:

        v = [
            float(r["node_type"] == t)
            for t in node_types
        ]

        v.extend(
            float(r.get("rna_type", "").strip() == t)
            for t in rna_types
        )

        raw_length = r.get(
            "sequence_length", ""
        ).strip()

        try:
            length = float(raw_length)
            valid_length = (
                math.isfinite(length)
                and length >= 0
            )
        except ValueError:
            valid_length = False
            length = 0

        v.append(
            math.log1p(length) / 10
            if valid_length else 0.0
        )

        v.append(float(valid_length))

        features.append(v)

    counts = collections.Counter(
        r["heldout_split"] for r in edges
    )

    if set(counts) != splits:
        fail("Suddivisione incompleta")

    return {
        "split_file": split_file,
        "nodes": nodes,
        "edges": edges,
        "node_index": node_index,
        "features": features,
        "positives": by_split,
        "negatives": negatives,
        "candidate_nodes": candidate_nodes,
        "unseen_count": len(unseen),
        "split_counts": counts
    }


def print_audit(data):

    print("=" * 68)
    print("AI PRIMETTA V3.3 | STEP46")
    print("GAE SPERIMENTALE — BENCHMARK STRUTTURALE")
    print("=" * 68)

    print("[NODI]", len(data["nodes"]))
    print("[ARCHI]", len(data["edges"]))

    print(
        "[FEATURE]",
        len(data["features"][0]),
        "caratteristiche intrinseche per nodo"
    )

    for name in ("train", "val", "test"):
        print(
            "[SPLIT]", name,
            "archi=", data["split_counts"][name],
            "coppie_somiglianza=",
            len(data["positives"][name])
        )

    print(
        "[CANDIDATI]",
        len(data["candidate_nodes"]),
        "nodi nelle relazioni di somiglianza"
    )

    print(
        "[NON_OSSERVATE]",
        data["unseen_count"]
    )

    print(
        "[PASS] Nessuna coppia presente "
        "in più di uno split"
    )

    print(
        "[VINCOLO] Negativi = collegamenti "
        "non osservati, non assenze "
        "biologicamente dimostrate"
    )

    print(
        "[VINCOLO] Benchmark transduttivo "
        "di ricostruzione delle somiglianze"
    )


def train(data):

    import numpy as np
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    from sklearn.metrics import (
        average_precision_score,
        roc_auc_score
    )

    from torch_geometric.nn import GCNConv
    from torch_geometric.utils import to_undirected

    torch.manual_seed(SEED)
    np.random.seed(SEED)
    random.seed(SEED)

    torch.set_num_threads(
        min(4, os.cpu_count() or 1)
    )

    index = data["node_index"]

    def tensor_pairs(pairs):
        ordered = sorted(pairs)
        return torch.tensor(
            [
                [index[a], index[b]]
                for a, b in ordered
            ],
            dtype=torch.long
        ).reshape(-1, 2)

    x = torch.tensor(
        data["features"],
        dtype=torch.float32
    )

    # Solamente gli archi del TRAIN
    # partecipano al message passing.

    train_links = [
        (
            index[r["source_node_id"]],
            index[r["target_node_id"]]
        )
        for r in data["edges"]
        if r["heldout_split"] == "train"
    ]

    edge_index = torch.tensor(
        train_links,
        dtype=torch.long
    ).T.contiguous()

    edge_index = to_undirected(
        edge_index,
        num_nodes=len(data["nodes"])
    )

    positives = {
        s: tensor_pairs(data["positives"][s])
        for s in ("train", "val", "test")
    }

    val_negative = tensor_pairs(
        data["negatives"]["val"]
    )

    test_negative = tensor_pairs(
        data["negatives"]["test"]
    )

    negative_pool = list(
        data["negatives"]["train_pool"]
    )

    class Encoder(nn.Module):

        def __init__(self, n_features):
            super().__init__()

            self.conv1 = GCNConv(
                n_features, HIDDEN
            )

            self.conv2 = GCNConv(
                HIDDEN, LATENT
            )

        def forward(self, features, edges):

            h = self.conv1(features, edges)
            h = F.relu(h)

            h = F.dropout(
                h, p=0.2, training=self.training
            )

            return self.conv2(h, edges)

    model = Encoder(x.shape[1])

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=0.01,
        weight_decay=0.00001
    )

    def logits(z, pairs):

        return (
            z[pairs[:, 0]]
            * z[pairs[:, 1]]
        ).sum(dim=1)

    def evaluate(z, pos, neg):

        scores = torch.cat([
            torch.sigmoid(logits(z, pos)),
            torch.sigmoid(logits(z, neg))
        ]).detach().cpu().numpy()

        labels = np.concatenate([
            np.ones(len(pos)),
            np.zeros(len(neg))
        ])

        return {
            "roc_auc": float(
                roc_auc_score(labels, scores)
            ),
            "average_precision": float(
                average_precision_score(
                    labels, scores
                )
            )
        }

    best_ap = -1.0
    best_epoch = 0
    best_state = None
    history = []

    print("[TRAINING] Avvio nuovo GAE isolato")

    for epoch in range(1, EPOCHS + 1):

        rng = random.Random(SEED + epoch)

        sampled_negative = rng.sample(
            negative_pool,
            len(positives["train"])
        )

        neg = tensor_pairs(sampled_negative)

        model.train()
        optimizer.zero_grad()

        z = model(x, edge_index)

        p = logits(z, positives["train"])
        n = logits(z, neg)

        targets = torch.cat([
            torch.ones_like(p),
            torch.zeros_like(n)
        ])

        loss = F.binary_cross_entropy_with_logits(
            torch.cat([p, n]),
            targets
        )

        loss.backward()
        optimizer.step()

        model.eval()

        with torch.no_grad():

            z_val = model(x, edge_index)

            metrics = evaluate(
                z_val,
                positives["val"],
                val_negative
            )

        history.append({
            "epoch": epoch,
            "train_loss": float(loss.item()),
            "val_auc": metrics["roc_auc"],
            "val_ap": metrics["average_precision"]
        })

        if metrics["average_precision"] > best_ap:

            best_ap = metrics["average_precision"]
            best_epoch = epoch

            best_state = copy.deepcopy(
                model.state_dict()
            )

        if epoch == 1 or epoch % 20 == 0:

            print(
                f"[EPOCH {epoch:03d}] "
                f"loss={loss.item():.5f} "
                f"val_auc={metrics['roc_auc']:.4f} "
                f"val_ap={metrics['average_precision']:.4f}"
            )

        if epoch - best_epoch >= PATIENCE:

            print(
                "[EARLY STOP] Epoca:",
                epoch
            )

            break

    if best_state is None:
        fail("Nessuno stato del modello disponibile")

    model.load_state_dict(best_state)
    model.eval()

    # TEST utilizzato una sola volta,
    # dopo la selezione via validation.

    with torch.no_grad():

        z_test = model(x, edge_index)

        test_metrics = evaluate(
            z_test,
            positives["test"],
            test_negative
        )

    print(
        "[BEST] Epoca selezionata:",
        best_epoch
    )

    print(
        "[TEST] AUC:",
        round(test_metrics["roc_auc"], 6)
    )

    print(
        "[TEST] AP:",
        round(
            test_metrics["average_precision"], 6
        )
    )

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    out = RUNS / (
        "STEP46_GAE_SAFE_" + stamp
    )

    out.mkdir(parents=True, exist_ok=False)

    checkpoint = out / "gae_safe_model.pt"

    torch.save({
        "model_state_dict": model.state_dict(),
        "model_type": "GCN_structural_similarity_GAE",
        "input_features": int(x.shape[1]),
        "hidden_channels": HIDDEN,
        "latent_channels": LATENT,
        "best_epoch": best_epoch,
        "seed": SEED
    }, checkpoint)

    history_path = out / "TRAIN_HISTORY.csv"

    with history_path.open(
        "w", newline="", encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "epoch",
                "train_loss",
                "val_auc",
                "val_ap"
            ]
        )

        writer.writeheader()
        writer.writerows(history)

    report = {
        "step": 46,
        "purpose": "structural_similarity_reconstruction",
        "biological_validation": False,
        "seed": SEED,
        "best_epoch": best_epoch,
        "epochs_executed": len(history),
        "nodes": len(data["nodes"]),
        "train_graph_edges_original": (
            data["split_counts"]["train"]
        ),
        "positive_similarity_pairs": {
            s: len(data["positives"][s])
            for s in ("train", "val", "test")
        },
        "negative_sampling": (
            "unobserved_pairs_between_similarity_nodes"
        ),
        "negative_sampling_limitation": (
            "Unobserved does not mean biologically negative."
        ),
        "node_features": (
            "node_type, rna_type, sequence_length"
        ),
        "node_feature_limitation": (
            "Transductive features; no external validation."
        ),
        "selection_metric": "validation_average_precision",
        "test": test_metrics,
        "input_split_sha256": digest(
            data["split_file"]
        ),
        "input_nodes_sha256": digest(NODES),
        "input_edges_sha256": digest(EDGES),
        "historical_checkpoint_modified": False,
        "gui_docker_modified": False
    }

    path = out / "RAPPORTO_STEP46.json"

    path.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ) + "\n",
        encoding="utf-8"
    )

    print("[CHECKPOINT]", checkpoint)
    print("[REPORT]", path)
    print("[OK] Training e test completati")
    print("[OK] Versione originale invariata")


def main():

    parser = argparse.ArgumentParser()

    modes = parser.add_mutually_exclusive_group(
        required=True
    )

    modes.add_argument("--audit", action="store_true")
    modes.add_argument("--run", action="store_true")

    args = parser.parse_args()

    data = prepare()

    print_audit(data)

    if args.run:
        train(data)

    else:
        print("[OK] Audit in sola lettura completato")


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:
        print(
            f"[ERRORE] {type(exc).__name__}: {exc}",
            file=sys.stderr
        )
        sys.exit(1)
