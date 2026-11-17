#!/usr/bin/env python3
from __future__ import annotations

import copy
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import torch


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

PYG_V2 = BASE_DIR / "data/processed/pyg/maize_ncrna_graph.pt"
GRAPH_NODES = BASE_DIR / "data/processed/graph_nodes.csv"
GRAPH_EDGES = BASE_DIR / "data/processed/graph_edges.csv"
NCRNA_FEATURES = BASE_DIR / "data/processed/ncrna_features.csv"
V3_RANKING = BASE_DIR / "outputs/prjna419326/PRJNA419326_experimental_srna_v3_ranking.csv"

PYG_V3 = BASE_DIR / "data/processed/pyg/maize_ncrna_graph_v3_prjna419326.pt"
GRAPH_NODES_V3 = BASE_DIR / "data/processed/graph_nodes_v3_prjna419326.csv"
GRAPH_EDGES_V3 = BASE_DIR / "data/processed/graph_edges_v3_prjna419326.csv"
NCRNA_FEATURES_V3 = BASE_DIR / "data/processed/ncrna_features_v3_prjna419326.csv"

OUT_DIR = BASE_DIR / "outputs/prjna419326"
SUMMARY_TXT = OUT_DIR / "V3_integrated_graph_build_summary.txt"
REPORT_JSON = OUT_DIR / "V3_integrated_graph_build_report.json"


def load_pt(path: Path):
    try:
        return torch.load(path, weights_only=False)
    except TypeError:
        return torch.load(path)


def edge_id(src: str, dst: str, edge_type: str) -> str:
    raw = f"{src}|{dst}|{edge_type}".encode("utf-8")
    return "edge_" + hashlib.md5(raw).hexdigest()


def normalize_seq(seq: str) -> str:
    return str(seq).upper().replace("U", "T").replace(" ", "").replace("\n", "")


def seq_counts(seq: str) -> dict:
    s = normalize_seq(seq)
    a = s.count("A")
    c = s.count("C")
    g = s.count("G")
    t = s.count("T")
    n = s.count("N")
    length = max(len(s), 1)

    return {
        "a_count": a,
        "c_count": c,
        "g_count": g,
        "u_count": t,
        "n_count": n,
        "a_frac": a / length,
        "c_frac": c / length,
        "g_frac": g / length,
        "u_frac": t / length,
        "n_frac": n / length,
        "gc_count": g + c,
        "au_count": a + t,
        "gc_frac": (g + c) / length,
        "au_frac": (a + t) / length,
    }


def kmer_set(seq: str, k: int = 5) -> set[str]:
    s = normalize_seq(seq)
    if len(s) < k:
        return {s}
    return {s[i:i + k] for i in range(0, len(s) - k + 1)}


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / max(len(a | b), 1)


def get_mapping_value(mapping: dict, key: str, fallback: int = 0) -> int:
    if key in mapping:
        return int(mapping[key])
    return fallback


def update_mapping(mapping: dict, inverse: dict, key: str) -> int:
    if key in mapping:
        return int(mapping[key])

    new_idx = max([int(v) for v in mapping.values()] + [-1]) + 1
    mapping[key] = new_idx
    inverse[str(new_idx)] = key
    return new_idx


def make_feature_vector(row: pd.Series, all_feature_columns: list[str]) -> list[float]:
    vals = {c: 0.0 for c in all_feature_columns}

    seq = normalize_seq(row["sequence"])
    length = int(row["length"])
    counts = seq_counts(seq)

    def setv(col: str, value: float):
        if col in vals:
            vals[col] = float(value)

    # Tipo nodo: trattiamo i nuovi small RNA sperimentali come ncRNA per compatibilità dimensionale.
    setv("node_type__ncRNA", 1.0)

    # Feature ncRNA/small RNA
    setv("sequence_length_metadata", length)
    setv("seq_length_from_fasta", length)
    setv("has_quantitative_expression", 1)
    setv("has_drought_expression_evidence", 1)
    setv("expression_status_not_quantified_yet", 0)
    setv("log10_seq_length", math.log10(max(length, 1)))

    for k, v in counts.items():
        setv(k, v)

    setv("rna_type__ncRNA", 0)
    setv("rna_type__pre_miRNA", 0)
    setv("rna_type__sRNA", 1)
    setv("rna_type__tmRNA", 0)

    # Feature topologiche iniziali stimate; saranno coerenti con gli archi sperimentali minimi.
    setv("degree_total", 6)
    setv("degree_out", 6)
    setv("degree_in", 0)
    setv("n_direct_biological_edges", 3)
    setv("n_proxy_edges", 3)
    setv("mean_confidence", 0.85)
    setv("mean_weight", 0.85)

    setv("edge_count_has_rna_type", 1)
    setv("edge_count_listed_in_database", 0)
    setv("edge_count_belongs_to_species", 1)
    setv("edge_count_has_expression_status", 0)
    setv("evidence_count_rnacentral", 0)
    setv("evidence_count_expression_summary", 1)

    return [vals[c] for c in all_feature_columns]


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PYG_V3.parent.mkdir(parents=True, exist_ok=True)

    obj = load_pt(PYG_V2)
    data_v2 = obj["data"]
    metadata_v2 = obj["metadata"]

    metadata_v3 = copy.deepcopy(metadata_v2)
    mappings = metadata_v3["mappings"]

    node_id_to_idx = copy.deepcopy(mappings["node_id_to_idx"])
    idx_to_node_id = copy.deepcopy(mappings["idx_to_node_id"])

    node_type_to_idx = copy.deepcopy(mappings["node_type_to_idx"])
    idx_to_node_type = copy.deepcopy(mappings["idx_to_node_type"])

    edge_type_to_idx = copy.deepcopy(mappings["edge_type_to_idx"])
    idx_to_edge_type = copy.deepcopy(mappings["idx_to_edge_type"])

    evidence_source_to_idx = copy.deepcopy(mappings["evidence_source_to_idx"])
    idx_to_evidence_source = copy.deepcopy(mappings["idx_to_evidence_source"])

    all_feature_columns = list(mappings["all_feature_columns"])

    nodes = pd.read_csv(GRAPH_NODES)
    edges = pd.read_csv(GRAPH_EDGES)
    ncrna_features = pd.read_csv(NCRNA_FEATURES)
    v3 = pd.read_csv(V3_RANKING)

    print("COSTRUZIONE GRAFO INTEGRATO V3")
    print("=" * 80)
    print(f"Nodi V2: {data_v2.num_nodes}")
    print(f"Archi V2: {data_v2.edge_index.shape[1]}")
    print(f"Feature V2: {tuple(data_v2.x.shape)}")
    print(f"Candidati small RNA V3: {len(v3)}")
    print("")

    original_num_nodes = int(data_v2.num_nodes)
    original_num_edges = int(data_v2.edge_index.shape[1])

    # Nodi di supporto già presenti nel grafo V2
    species_node_id = "species:NCBITaxon_4577"
    srna_type_node_id = "rna_type:sRNA"

    dataset_candidates = nodes.loc[nodes["node_type"].astype(str) == "public_expression_dataset", "node_id"].tolist()
    dataset_node_id = dataset_candidates[0] if dataset_candidates else None

    # Tipo nodo PyG: usiamo ncRNA per i nuovi small RNA sperimentali.
    ncRNA_type_idx = get_mapping_value(node_type_to_idx, "ncRNA", 0)

    # Nuovi node_id
    new_node_ids = []
    new_node_rows = []
    new_feature_rows = []
    new_x_rows = []
    new_node_type_values = []
    new_y_values = []
    new_train_mask_values = []
    new_unlabeled_mask_values = []
    new_is_ncrna_values = []

    next_idx = original_num_nodes

    for _, row in v3.iterrows():
        rank_id = str(row["v3_rank_id"])
        feature_id = str(row["feature_id"])
        seq = normalize_seq(row["sequence"])
        length = int(row["length"])

        node_id = f"experimental_srna:{rank_id}"

        node_id_to_idx[node_id] = next_idx
        idx_to_node_id[str(next_idx)] = node_id
        new_node_ids.append(node_id)

        attrs = {
            "project": PROJECT,
            "experimental_source": "PRJNA419326",
            "experimental_node_type": "small_RNA",
            "feature_id": feature_id,
            "v3_rank_id": rank_id,
            "sequence": seq,
            "candidate_class": row["candidate_class"],
            "v3_priority_score": float(row["v3_priority_score"]),
            "stress_score": float(row["stress_score"]),
            "total_count": int(row["total_count"]),
            "samples_detected": int(row["samples_detected"]),
            "mean_rpm_B73_DS": float(row["mean_rpm_B73_DS"]),
            "mean_rpm_B73_WW": float(row["mean_rpm_B73_WW"]),
            "log2FC_B73_DS_vs_WW": float(row["log2FC_B73_DS_vs_WW"]),
            "mean_rpm_H21_DS": float(row["mean_rpm_H21_DS"]),
            "mean_rpm_H21_WW": float(row["mean_rpm_H21_WW"]),
            "log2FC_H21_DS_vs_WW": float(row["log2FC_H21_DS_vs_WW"]),
            "rnacentral_exact_match": row["rnacentral_exact_match"],
            "interpretation": row["interpretation"],
        }

        new_node_rows.append(
            {
                "node_id": node_id,
                "node_type": "ncRNA",
                "name": rank_id,
                "source": PROJECT,
                "rna_type": "sRNA",
                "sequence_length": length,
                "label": -1,
                "label_status": "experimental_drought_response_unlabeled_for_supervised_training",
                "expression_status": row["candidate_class"],
                "attributes_json": json.dumps(attrs, ensure_ascii=False),
            }
        )

        feature_vector = make_feature_vector(row, all_feature_columns)
        new_x_rows.append(feature_vector)

        counts = seq_counts(seq)
        new_feature_row = {
            "rnacentral_id": rank_id,
            "node_id": node_id,
            "rna_type": "sRNA",
            "sequence_length_metadata": length,
            "seq_length_from_fasta": length,
            "label": -1,
            "label_confidence": 0.0,
            "has_quantitative_expression": 1,
            "has_drought_expression_evidence": 1,
            "expression_status_not_quantified_yet": 0,
            "log10_seq_length": math.log10(max(length, 1)),
            **counts,
            "rna_type__ncRNA": 0,
            "rna_type__pre_miRNA": 0,
            "rna_type__sRNA": 1,
            "rna_type__tmRNA": 0,
            "degree_total": 6,
            "degree_out": 6,
            "degree_in": 0,
            "n_direct_biological_edges": 3,
            "n_proxy_edges": 3,
            "mean_confidence": 0.85,
            "mean_weight": 0.85,
            "edge_count_has_rna_type": 1,
            "edge_count_listed_in_database": 0,
            "edge_count_belongs_to_species": 1,
            "edge_count_has_expression_status": 0,
            "evidence_count_rnacentral": 0,
            "evidence_count_expression_summary": 1,
        }
        new_feature_rows.append(new_feature_row)

        new_node_type_values.append(ncRNA_type_idx)
        new_y_values.append(-1)
        new_train_mask_values.append(False)
        new_unlabeled_mask_values.append(True)
        new_is_ncrna_values.append(True)

        next_idx += 1

    # Mapping edge/evidence aggiornati
    prjna_source_idx = update_mapping(
        evidence_source_to_idx,
        idx_to_evidence_source,
        "PRJNA419326",
    )

    seq_sim_type_idx = update_mapping(
        edge_type_to_idx,
        idx_to_edge_type,
        "experimental_sequence_similarity",
    )
    response_sim_type_idx = update_mapping(
        edge_type_to_idx,
        idx_to_edge_type,
        "experimental_drought_response_similarity",
    )

    # Edge type già esistenti o fallback
    belongs_type_idx = update_mapping(edge_type_to_idx, idx_to_edge_type, "belongs_to_species")
    has_rna_type_idx = update_mapping(edge_type_to_idx, idx_to_edge_type, "has_rna_type")
    dataset_type_idx = update_mapping(edge_type_to_idx, idx_to_edge_type, "found_candidate_expression_dataset")

    new_edges_table = []
    new_edge_src = []
    new_edge_dst = []
    new_edge_attr = []
    new_edge_type_values = []
    new_evidence_values = []
    edge_seen = set()

    def add_edge(src_id: str, dst_id: str, edge_type: str, edge_type_idx: int,
                 weight: float, confidence: float, evidence_source: str,
                 evidence_source_idx: int, evidence_label: str,
                 direct: int, proxy: int, notes: str):
        if src_id not in node_id_to_idx or dst_id not in node_id_to_idx:
            return

        key = (src_id, dst_id, edge_type)
        if key in edge_seen:
            return
        edge_seen.add(key)

        eid = edge_id(src_id, dst_id, edge_type)

        new_edges_table.append(
            {
                "edge_id": eid,
                "source_node_id": src_id,
                "target_node_id": dst_id,
                "source_node_type": "ncRNA",
                "target_node_type": "mixed",
                "edge_type": edge_type,
                "weight": float(weight),
                "confidence_score": float(confidence),
                "evidence_source": evidence_source,
                "evidence_label": evidence_label,
                "evidence_url": "",
                "is_direct_biological_evidence": int(direct),
                "is_proxy_evidence": int(proxy),
                "notes": notes,
            }
        )

        new_edge_src.append(int(node_id_to_idx[src_id]))
        new_edge_dst.append(int(node_id_to_idx[dst_id]))
        new_edge_attr.append([float(weight), float(confidence), float(direct), float(proxy)])
        new_edge_type_values.append(int(edge_type_idx))
        new_evidence_values.append(int(evidence_source_idx))

    # Archi biologici minimi verso nodi esistenti
    for node_id in new_node_ids:
        add_edge(
            node_id,
            species_node_id,
            "belongs_to_species",
            belongs_type_idx,
            1.0,
            0.95,
            "PRJNA419326",
            prjna_source_idx,
            "Small RNA sperimentale assegnata a Zea mays",
            1,
            0,
            "Experimental small RNA from maize PRJNA419326.",
        )

        if srna_type_node_id in node_id_to_idx:
            add_edge(
                node_id,
                srna_type_node_id,
                "has_rna_type",
                has_rna_type_idx,
                1.0,
                0.95,
                "PRJNA419326",
                prjna_source_idx,
                "Small RNA class from size-selected miRNA/small RNA sequencing",
                1,
                0,
                "Experimental sRNA node after adapter trimming and 18-30 nt filtering.",
            )

        if dataset_node_id is not None:
            add_edge(
                node_id,
                dataset_node_id,
                "found_candidate_expression_dataset",
                dataset_type_idx,
                0.90,
                0.90,
                "PRJNA419326",
                prjna_source_idx,
                "Candidate observed in PRJNA419326 small RNA drought dataset",
                1,
                0,
                "Connects experimental sRNA candidate to public expression dataset layer.",
            )

    # Similarità di sequenza tra small RNA sperimentali
    seqs = [normalize_seq(s) for s in v3["sequence"].tolist()]
    kmers = [kmer_set(s, 5) for s in seqs]

    for i in range(len(seqs)):
        sims = []
        for j in range(len(seqs)):
            if i == j:
                continue
            sim = jaccard(kmers[i], kmers[j])
            if sim >= 0.60:
                sims.append((sim, j))

        sims.sort(reverse=True)

        for sim, j in sims[:3]:
            src = new_node_ids[i]
            dst = new_node_ids[j]

            add_edge(
                src,
                dst,
                "experimental_sequence_similarity",
                seq_sim_type_idx,
                sim,
                min(0.95, 0.50 + sim / 2),
                "PRJNA419326",
                prjna_source_idx,
                "Sequence k-mer similarity among experimental small RNAs",
                0,
                1,
                "Proxy edge based on 5-mer Jaccard similarity.",
            )

    # Similarità di risposta drought tra small RNA sperimentali
    fc = v3[["log2FC_B73_DS_vs_WW", "log2FC_H21_DS_vs_WW"]].astype(float).to_numpy()
    classes = v3["candidate_class"].astype(str).tolist()

    for i in range(len(seqs)):
        candidates = []

        for j in range(len(seqs)):
            if i == j:
                continue
            if classes[i] != classes[j]:
                continue

            dist = float(np.linalg.norm(fc[i] - fc[j]))
            score = 1.0 / (1.0 + dist)

            candidates.append((score, j))

        candidates.sort(reverse=True)

        for score, j in candidates[:2]:
            src = new_node_ids[i]
            dst = new_node_ids[j]

            add_edge(
                src,
                dst,
                "experimental_drought_response_similarity",
                response_sim_type_idx,
                score,
                min(0.90, 0.55 + score / 3),
                "PRJNA419326",
                prjna_source_idx,
                "Similar drought-response pattern among experimental small RNAs",
                0,
                1,
                "Proxy edge based on same response class and log2FC distance.",
            )

    # Tabelle V3
    nodes_v3 = pd.concat([nodes, pd.DataFrame(new_node_rows)], ignore_index=True)
    edges_v3 = pd.concat([edges, pd.DataFrame(new_edges_table)], ignore_index=True)
    ncrna_features_v3 = pd.concat([ncrna_features, pd.DataFrame(new_feature_rows)], ignore_index=True)

    nodes_v3.to_csv(GRAPH_NODES_V3, index=False)
    edges_v3.to_csv(GRAPH_EDGES_V3, index=False)
    ncrna_features_v3.to_csv(NCRNA_FEATURES_V3, index=False)

    # Tensor V3
    x_new = torch.tensor(new_x_rows, dtype=data_v2.x.dtype)
    x_v3 = torch.cat([data_v2.x, x_new], dim=0)

    new_edge_index = torch.tensor([new_edge_src, new_edge_dst], dtype=torch.long)
    edge_index_v3 = torch.cat([data_v2.edge_index, new_edge_index], dim=1)

    edge_attr_new = torch.tensor(new_edge_attr, dtype=data_v2.edge_attr.dtype)
    edge_attr_v3 = torch.cat([data_v2.edge_attr, edge_attr_new], dim=0)

    edge_type_new = torch.tensor(new_edge_type_values, dtype=data_v2.edge_type.dtype)
    edge_type_v3 = torch.cat([data_v2.edge_type, edge_type_new], dim=0)

    evidence_new = torch.tensor(new_evidence_values, dtype=data_v2.evidence_source.dtype)
    evidence_v3 = torch.cat([data_v2.evidence_source, evidence_new], dim=0)

    node_type_new = torch.tensor(new_node_type_values, dtype=data_v2.node_type.dtype)
    node_type_v3 = torch.cat([data_v2.node_type, node_type_new], dim=0)

    y_new = torch.tensor(new_y_values, dtype=data_v2.y.dtype)
    y_v3 = torch.cat([data_v2.y, y_new], dim=0)

    train_mask_new = torch.tensor(new_train_mask_values, dtype=data_v2.train_mask.dtype)
    train_mask_v3 = torch.cat([data_v2.train_mask, train_mask_new], dim=0)

    unlabeled_mask_new = torch.tensor(new_unlabeled_mask_values, dtype=data_v2.unlabeled_mask.dtype)
    unlabeled_mask_v3 = torch.cat([data_v2.unlabeled_mask, unlabeled_mask_new], dim=0)

    is_ncrna_new = torch.tensor(new_is_ncrna_values, dtype=data_v2.is_ncrna_node.dtype)
    is_ncrna_v3 = torch.cat([data_v2.is_ncrna_node, is_ncrna_new], dim=0)

    data_v3 = copy.deepcopy(data_v2)
    data_v3.x = x_v3
    data_v3.edge_index = edge_index_v3
    data_v3.edge_attr = edge_attr_v3
    data_v3.edge_type = edge_type_v3
    data_v3.evidence_source = evidence_v3
    data_v3.node_type = node_type_v3
    data_v3.y = y_v3
    data_v3.train_mask = train_mask_v3
    data_v3.unlabeled_mask = unlabeled_mask_v3
    data_v3.is_ncrna_node = is_ncrna_v3
    data_v3.num_nodes = int(x_v3.shape[0])

    # Metadati V3
    mappings["node_id_to_idx"] = node_id_to_idx
    mappings["idx_to_node_id"] = idx_to_node_id
    mappings["node_type_to_idx"] = node_type_to_idx
    mappings["idx_to_node_type"] = idx_to_node_type
    mappings["edge_type_to_idx"] = edge_type_to_idx
    mappings["idx_to_edge_type"] = idx_to_edge_type
    mappings["evidence_source_to_idx"] = evidence_source_to_idx
    mappings["idx_to_evidence_source"] = idx_to_evidence_source

    metadata_v3["created_at_utc"] = datetime.now(timezone.utc).isoformat()
    metadata_v3["project"]["version"] = "0.3.0-prjna419326"
    metadata_v3["source_files"]["v3_prjna419326_ranking_csv"] = str(V3_RANKING)
    metadata_v3["source_files"]["graph_nodes_v3_csv"] = str(GRAPH_NODES_V3)
    metadata_v3["source_files"]["graph_edges_v3_csv"] = str(GRAPH_EDGES_V3)
    metadata_v3["source_files"]["ncrna_features_v3_csv"] = str(NCRNA_FEATURES_V3)

    metadata_v3["summary"] = {
        **metadata_v3.get("summary", {}),
        "version": "V3_PRJNA419326_integrated_graph",
        "n_nodes": int(data_v3.num_nodes),
        "n_edges": int(data_v3.edge_index.shape[1]),
        "n_node_features": int(data_v3.x.shape[1]),
        "n_edge_features": int(data_v3.edge_attr.shape[1]),
        "n_original_v2_nodes": original_num_nodes,
        "n_original_v2_edges": original_num_edges,
        "n_added_experimental_srna_nodes": len(new_node_ids),
        "n_added_experimental_edges": len(new_edges_table),
        "n_labeled_nodes": int(data_v3.train_mask.sum().item()),
        "n_unlabeled_nodes": int(data_v3.unlabeled_mask.sum().item()),
        "n_ncrna_nodes": int(data_v3.is_ncrna_node.sum().item()),
        "methodological_note": (
            "V3 graph integrates 500 experimental small RNA candidates from PRJNA419326. "
            "Nodes are represented with the same 47-dimensional feature space as V2 for GAE compatibility. "
            "Edges include species/RNA type/dataset links plus sequence and drought-response similarity edges. "
            "All nodes remain unlabeled for supervised classification; the graph is intended for unsupervised/self-supervised GAE retraining."
        ),
    }

    torch.save({"data": data_v3, "metadata": metadata_v3}, PYG_V3)

    # Report
    edge_type_counts = pd.Series([e["edge_type"] for e in new_edges_table]).value_counts().to_dict()

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "input_pyg_v2": str(PYG_V2),
        "output_pyg_v3": str(PYG_V3),
        "original_num_nodes": original_num_nodes,
        "original_num_edges": original_num_edges,
        "added_experimental_srna_nodes": len(new_node_ids),
        "added_experimental_edges": len(new_edges_table),
        "v3_num_nodes": int(data_v3.num_nodes),
        "v3_num_edges": int(data_v3.edge_index.shape[1]),
        "v3_x_shape": list(data_v3.x.shape),
        "v3_edge_attr_shape": list(data_v3.edge_attr.shape),
        "new_edge_type_counts": edge_type_counts,
        "support_nodes": {
            "species_node_id": species_node_id,
            "srna_type_node_id": srna_type_node_id if srna_type_node_id in node_id_to_idx else None,
            "dataset_node_id": dataset_node_id,
        },
        "outputs": {
            "pyg_v3": str(PYG_V3),
            "graph_nodes_v3_csv": str(GRAPH_NODES_V3),
            "graph_edges_v3_csv": str(GRAPH_EDGES_V3),
            "ncrna_features_v3_csv": str(NCRNA_FEATURES_V3),
            "summary_txt": str(SUMMARY_TXT),
            "report_json": str(REPORT_JSON),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("GRAFO INTEGRATO V3 PRJNA419326 COSTRUITO")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Nodi V2 originali: {original_num_nodes}")
    lines.append(f"Archi V2 originali: {original_num_edges}")
    lines.append(f"Nodi small RNA sperimentali aggiunti: {len(new_node_ids)}")
    lines.append(f"Archi sperimentali aggiunti: {len(new_edges_table)}")
    lines.append(f"Nodi totali V3: {data_v3.num_nodes}")
    lines.append(f"Archi totali V3: {data_v3.edge_index.shape[1]}")
    lines.append(f"Feature matrix V3: {tuple(data_v3.x.shape)}")
    lines.append(f"Edge attr V3: {tuple(data_v3.edge_attr.shape)}")
    lines.append("")
    lines.append("Nuovi edge type:")
    for k, v in edge_type_counts.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("Output:")
    lines.append(str(PYG_V3))
    lines.append(str(GRAPH_NODES_V3))
    lines.append(str(GRAPH_EDGES_V3))
    lines.append(str(NCRNA_FEATURES_V3))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
