#!/usr/bin/env python3
"""
12_build_exploratory_ncrna_ranking.py

Costruisce un ranking esplorativo dei candidati ncRNA usando:
- embedding non supervisionati;
- matrice feature ncRNA;
- risultati di similarità top-k.

Nota metodologica:
il ranking NON rappresenta una prova di risposta alla siccità.
È una priorità esplorativa utile per selezionare candidati da validare
con ulteriori evidenze biologiche, RNA-seq, letteratura o esperimenti.

Input:
- data/processed/embeddings/ncrna_embeddings.csv
- data/processed/ncrna_features.csv
- data/processed/similarity/ncrna_similarity_topk.csv

Output:
- data/processed/ranking/ncrna_exploratory_ranking.csv
- data/processed/ranking/ncrna_exploratory_top50.csv
- outputs/ranking/ncrna_exploratory_ranking_report.json
- outputs/ranking/ncrna_exploratory_ranking_summary.txt
- outputs/ranking/plots/top20_exploratory_ranking.png
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import yaml

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except Exception:
    plt = None
    MATPLOTLIB_AVAILABLE = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Costruisce ranking esplorativo dei candidati ncRNA."
    )
    parser.add_argument("--config", default="configs/config_public.yaml")
    parser.add_argument("--ncrna-embeddings-csv", default="data/processed/embeddings/ncrna_embeddings.csv")
    parser.add_argument("--ncrna-features-csv", default="data/processed/ncrna_features.csv")
    parser.add_argument("--similarity-topk-csv", default="data/processed/similarity/ncrna_similarity_topk.csv")
    parser.add_argument("--output-dir", default="data/processed/ranking")
    parser.add_argument("--report-dir", default="outputs/ranking")
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)


def load_yaml(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        obj = yaml.safe_load(f)
    return obj if isinstance(obj, dict) else {}


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV non trovato: {path}")
    with path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"CSV vuoto: {path}")
    return rows


def write_csv(rows: List[Dict[str, Any]], path: Path, fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value in {"", "NA", None}:
            return default
        return float(value)
    except Exception:
        return default


def safe_int(value: Any, default: int = 0) -> int:
    try:
        if value in {"", "NA", None}:
            return default
        return int(float(str(value)))
    except Exception:
        return default


def detect_embedding_columns(rows: List[Dict[str, str]]) -> List[str]:
    cols = [c for c in rows[0].keys() if c.startswith("emb_")]
    if not cols:
        raise ValueError("Nessuna colonna embedding trovata.")
    return cols


def index_by_node_id(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    return {
        row.get("node_id", ""): row
        for row in rows
        if row.get("node_id", "")
    }


def minmax_scores(values: Dict[str, float]) -> Dict[str, float]:
    if not values:
        return {}

    vals = list(values.values())
    vmin = min(vals)
    vmax = max(vals)

    if math.isclose(vmin, vmax):
        return {k: 1.0 for k in values}

    return {
        k: (v - vmin) / (vmax - vmin)
        for k, v in values.items()
    }


def inverse_distance_scores(values: Dict[str, float]) -> Dict[str, float]:
    raw = {
        k: 1.0 / (1.0 + max(0.0, v))
        for k, v in values.items()
    }
    return minmax_scores(raw)


def rna_type_priority(rna_type: str) -> float:
    """
    Priorità esplorativa, non evidenza biologica definitiva.

    pre_miRNA e ncRNA sono messi più in alto perché più direttamente utili
    per un progetto di regolazione genica/non-coding RNA; tmRNA resta valido
    ma meno prioritario per una prima lista di candidati regolativi.
    """
    priorities = {
        "pre_miRNA": 1.00,
        "ncRNA": 0.85,
        "sRNA": 0.75,
        "tmRNA": 0.35,
    }
    return priorities.get(rna_type, 0.50)


def build_similarity_stats(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, Any]]:
    grouped: Dict[str, List[Dict[str, str]]] = defaultdict(list)

    for row in rows:
        node_id = row.get("query_node_id", "")
        if node_id:
            grouped[node_id].append(row)

    stats: Dict[str, Dict[str, Any]] = {}

    for node_id, group in grouped.items():
        similarities = [safe_float(r.get("cosine_similarity"), 0.0) for r in group]
        same_type = [safe_int(r.get("same_rna_type"), 0) for r in group]
        neighbors = [
            r.get("similar_rnacentral_id", "")
            for r in sorted(group, key=lambda x: safe_int(x.get("query_rank"), 999))
        ]

        similar_types = sorted(set(r.get("similar_rna_type", "") for r in group if r.get("similar_rna_type", "")))

        stats[node_id] = {
            "topk_count": len(group),
            "mean_topk_similarity": float(np.mean(similarities)) if similarities else 0.0,
            "max_topk_similarity": float(np.max(similarities)) if similarities else 0.0,
            "min_topk_similarity": float(np.min(similarities)) if similarities else 0.0,
            "topk_same_type_fraction": float(np.mean(same_type)) if same_type else 0.0,
            "similar_rna_type_diversity": len(similar_types),
            "top3_similar_ids": ";".join(neighbors[:3]),
        }

    return stats


def build_embedding_vectors(rows: List[Dict[str, str]], emb_cols: List[str]) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {}

    for row in rows:
        node_id = row.get("node_id", "")
        if not node_id:
            continue

        vector = np.asarray([safe_float(row.get(col), 0.0) for col in emb_cols], dtype=np.float64)
        out[node_id] = vector

    return out


def build_centroid_distances(vectors: Dict[str, np.ndarray]) -> Dict[str, float]:
    if not vectors:
        return {}

    matrix = np.vstack(list(vectors.values()))
    centroid = matrix.mean(axis=0)

    return {
        node_id: float(np.linalg.norm(vec - centroid))
        for node_id, vec in vectors.items()
    }


def sequence_quality_score(row: Dict[str, str]) -> float:
    n_frac = safe_float(row.get("n_frac"), 0.0)
    has_quant = safe_int(row.get("has_quantitative_expression"), 0)

    # In questa fase has_quantitative_expression è quasi sempre 0.
    # La qualità è soprattutto assenza di N/ambiguità.
    score = 1.0 - min(1.0, max(0.0, n_frac))
    if has_quant:
        score = min(1.0, score + 0.05)

    return score


def gc_typicality_scores(feature_rows: List[Dict[str, str]]) -> Dict[str, float]:
    gc_values = []
    gc_by_node = {}

    for row in feature_rows:
        node_id = row.get("node_id", "")
        gc = safe_float(row.get("gc_frac"), 0.0)
        if node_id:
            gc_by_node[node_id] = gc
            gc_values.append(gc)

    if not gc_values:
        return {}

    median_gc = float(np.median(gc_values))
    max_dev = max(abs(gc - median_gc) for gc in gc_values)

    if math.isclose(max_dev, 0.0):
        return {node_id: 1.0 for node_id in gc_by_node}

    return {
        node_id: max(0.0, 1.0 - abs(gc - median_gc) / max_dev)
        for node_id, gc in gc_by_node.items()
    }


def build_ranking_rows(
    embedding_rows: List[Dict[str, str]],
    feature_rows: List[Dict[str, str]],
    similarity_stats: Dict[str, Dict[str, Any]],
) -> List[Dict[str, Any]]:
    emb_cols = detect_embedding_columns(embedding_rows)
    vectors = build_embedding_vectors(embedding_rows, emb_cols)

    feature_by_node = index_by_node_id(feature_rows)
    centroid_dist = build_centroid_distances(vectors)

    mean_sim_raw = {
        node_id: stats.get("mean_topk_similarity", 0.0)
        for node_id, stats in similarity_stats.items()
    }

    max_sim_raw = {
        node_id: stats.get("max_topk_similarity", 0.0)
        for node_id, stats in similarity_stats.items()
    }

    centrality_score = inverse_distance_scores(centroid_dist)
    mean_sim_score = minmax_scores(mean_sim_raw)
    max_sim_score = minmax_scores(max_sim_raw)
    gc_score = gc_typicality_scores(feature_rows)

    rows: List[Dict[str, Any]] = []

    for emb in embedding_rows:
        node_id = emb.get("node_id", "")
        if not node_id:
            continue

        feat = feature_by_node.get(node_id, {})
        stats = similarity_stats.get(node_id, {})

        rna_type = emb.get("rna_type", feat.get("rna_type", "unknown"))
        rnacentral_id = emb.get("name", node_id.replace("ncrna:", ""))

        mean_topk = float(stats.get("mean_topk_similarity", 0.0))
        max_topk = float(stats.get("max_topk_similarity", 0.0))
        same_type_fraction = float(stats.get("topk_same_type_fraction", 0.0))

        component_mean_similarity = mean_sim_score.get(node_id, 0.0)
        component_max_similarity = max_sim_score.get(node_id, 0.0)
        component_centrality = centrality_score.get(node_id, 0.0)
        component_rna_type_priority = rna_type_priority(rna_type)
        component_gc_typicality = gc_score.get(node_id, 0.0)
        component_sequence_quality = sequence_quality_score(feat)

        exploratory_score = (
            0.30 * component_mean_similarity
            + 0.15 * component_max_similarity
            + 0.20 * component_centrality
            + 0.20 * component_rna_type_priority
            + 0.10 * component_gc_typicality
            + 0.05 * component_sequence_quality
        )

        row = {
            "rnacentral_id": rnacentral_id,
            "node_id": node_id,
            "rna_type": rna_type,
            "sequence_length": emb.get("sequence_length", feat.get("seq_length_from_fasta", "")),
            "gc_frac": emb.get("gc_frac", feat.get("gc_frac", "")),
            "label": emb.get("label", "-1"),

            "exploratory_rank_score": round(float(exploratory_score), 8),

            "component_mean_similarity": round(float(component_mean_similarity), 8),
            "component_max_similarity": round(float(component_max_similarity), 8),
            "component_embedding_centrality": round(float(component_centrality), 8),
            "component_rna_type_priority": round(float(component_rna_type_priority), 8),
            "component_gc_typicality": round(float(component_gc_typicality), 8),
            "component_sequence_quality": round(float(component_sequence_quality), 8),

            "mean_topk_similarity": round(mean_topk, 8),
            "max_topk_similarity": round(max_topk, 8),
            "min_topk_similarity": round(float(stats.get("min_topk_similarity", 0.0)), 8),
            "topk_same_type_fraction": round(same_type_fraction, 8),
            "similar_rna_type_diversity": int(stats.get("similar_rna_type_diversity", 0)),
            "top3_similar_ids": stats.get("top3_similar_ids", ""),

            "degree_total": feat.get("degree_total", ""),
            "mean_confidence": feat.get("mean_confidence", ""),
            "mean_weight": feat.get("mean_weight", ""),
            "has_quantitative_expression": feat.get("has_quantitative_expression", "0"),
            "has_drought_expression_evidence": feat.get("has_drought_expression_evidence", "0"),

            "methodological_note": (
                "Exploratory unsupervised ranking. Not evidence of drought responsiveness."
            ),
        }

        rows.append(row)

    rows.sort(key=lambda x: x["exploratory_rank_score"], reverse=True)

    for i, row in enumerate(rows, start=1):
        row["exploratory_rank"] = i

    return rows


def save_top20_plot(rows: List[Dict[str, Any]], output_path: Path) -> bool:
    if not MATPLOTLIB_AVAILABLE or plt is None:
        return False

    if not rows:
        return False

    top = rows[:20]
    labels = [row["rnacentral_id"] for row in top]
    values = [float(row["exploratory_rank_score"]) for row in top]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(12, 7))
    plt.bar(labels, values)
    plt.title("Top 20 ncRNA exploratory ranking")
    plt.ylabel("Exploratory rank score")
    plt.xticks(rotation=70, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()

    return True


def summarize_scores(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    scores = [float(row["exploratory_rank_score"]) for row in rows]
    if not scores:
        return {}

    return {
        "count": len(scores),
        "mean": round(float(np.mean(scores)), 8),
        "std": round(float(np.std(scores)), 8),
        "min": round(float(np.min(scores)), 8),
        "median": round(float(np.median(scores)), 8),
        "max": round(float(np.max(scores)), 8),
    }


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    embeddings_path = Path(args.ncrna_embeddings_csv)
    features_path = Path(args.ncrna_features_csv)
    similarity_path = Path(args.similarity_topk_csv)

    output_dir = Path(args.output_dir)
    report_dir = Path(args.report_dir)
    plots_dir = report_dir / "plots"

    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_yaml(config_path)

    logging.info("Caricamento embedding ncRNA: %s", embeddings_path)
    embedding_rows = read_csv(embeddings_path)

    logging.info("Caricamento feature ncRNA: %s", features_path)
    feature_rows = read_csv(features_path)

    logging.info("Caricamento top-k similarity: %s", similarity_path)
    similarity_rows = read_csv(similarity_path)

    logging.info("Calcolo statistiche di similarità")
    similarity_stats = build_similarity_stats(similarity_rows)

    logging.info("Costruzione ranking esplorativo")
    ranking_rows = build_ranking_rows(
        embedding_rows=embedding_rows,
        feature_rows=feature_rows,
        similarity_stats=similarity_stats,
    )

    top_rows = ranking_rows[: args.top_n]

    ranking_csv = output_dir / "ncrna_exploratory_ranking.csv"
    top_csv = output_dir / f"ncrna_exploratory_top{args.top_n}.csv"

    fieldnames = [
        "exploratory_rank",
        "rnacentral_id",
        "node_id",
        "rna_type",
        "sequence_length",
        "gc_frac",
        "label",
        "exploratory_rank_score",
        "component_mean_similarity",
        "component_max_similarity",
        "component_embedding_centrality",
        "component_rna_type_priority",
        "component_gc_typicality",
        "component_sequence_quality",
        "mean_topk_similarity",
        "max_topk_similarity",
        "min_topk_similarity",
        "topk_same_type_fraction",
        "similar_rna_type_diversity",
        "top3_similar_ids",
        "degree_total",
        "mean_confidence",
        "mean_weight",
        "has_quantitative_expression",
        "has_drought_expression_evidence",
        "methodological_note",
    ]

    write_csv(ranking_rows, ranking_csv, fieldnames)
    write_csv(top_rows, top_csv, fieldnames)

    plots_created = []
    plot_path = plots_dir / "top20_exploratory_ranking.png"
    if save_top20_plot(ranking_rows, plot_path):
        plots_created.append("plots/top20_exploratory_ranking.png")

    rna_type_counter_all = Counter(row["rna_type"] for row in ranking_rows)
    rna_type_counter_top = Counter(row["rna_type"] for row in top_rows)

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "project": cfg.get("project", {}),
        "inputs": {
            "ncrna_embeddings_csv": str(embeddings_path),
            "ncrna_features_csv": str(features_path),
            "similarity_topk_csv": str(similarity_path),
        },
        "outputs": {
            "ranking_csv": str(ranking_csv),
            "top_csv": str(top_csv),
            "report_json": str(report_dir / "ncrna_exploratory_ranking_report.json"),
            "summary_txt": str(report_dir / "ncrna_exploratory_ranking_summary.txt"),
            "plots_dir": str(plots_dir),
        },
        "summary": {
            "n_ranked_ncrna": len(ranking_rows),
            "top_n": args.top_n,
            "ranking_score_summary": summarize_scores(ranking_rows),
            "rna_type_distribution_all": dict(rna_type_counter_all),
            "rna_type_distribution_top": dict(rna_type_counter_top),
            "top_10_candidates": [
                {
                    "rank": row["exploratory_rank"],
                    "rnacentral_id": row["rnacentral_id"],
                    "rna_type": row["rna_type"],
                    "score": row["exploratory_rank_score"],
                    "top3_similar_ids": row["top3_similar_ids"],
                }
                for row in ranking_rows[:10]
            ],
            "plots_created": plots_created,
            "methodological_note": (
                "Ranking is exploratory and unsupervised. It should be used to select "
                "candidate ncRNAs for literature review, expression integration, or biological "
                "validation. It is not a supervised drought-response prediction."
            ),
        },
    }

    report_path = report_dir / "ncrna_exploratory_ranking_report.json"
    with report_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    summary_lines = []
    summary_lines.append("ncRNA EXPLORATORY RANKING SUMMARY")
    summary_lines.append("=" * 80)
    summary_lines.append(f"Ranked ncRNA: {len(ranking_rows)}")
    summary_lines.append(f"Top N exported: {args.top_n}")
    summary_lines.append("")
    summary_lines.append(f"RNA type distribution all: {dict(rna_type_counter_all)}")
    summary_lines.append(f"RNA type distribution top {args.top_n}: {dict(rna_type_counter_top)}")
    summary_lines.append("")
    summary_lines.append(f"Score summary: {report['summary']['ranking_score_summary']}")
    summary_lines.append("")
    summary_lines.append("Top 10 candidates:")
    for row in ranking_rows[:10]:
        summary_lines.append(
            f"{row['exploratory_rank']}. {row['rnacentral_id']} | "
            f"{row['rna_type']} | score={row['exploratory_rank_score']} | "
            f"top3={row['top3_similar_ids']}"
        )
    summary_lines.append("")
    summary_lines.append(f"Plots created: {plots_created}")
    summary_lines.append("")
    summary_lines.append("Methodological note:")
    summary_lines.append(report["summary"]["methodological_note"])

    summary_path = report_dir / "ncrna_exploratory_ranking_summary.txt"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")

    print("\n".join(summary_lines))

    logging.info("Ranking CSV: %s", ranking_csv)
    logging.info("Top CSV: %s", top_csv)
    logging.info("Report: %s", report_path)
    logging.info("Summary: %s", summary_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
