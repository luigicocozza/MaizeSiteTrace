#!/usr/bin/env python3
"""
11_build_ncrna_similarity_search.py

Costruisce una ricerca di similarità tra ncRNA usando gli embedding
non supervisionati generati nello STEP 16.

Input:
- data/processed/embeddings/ncrna_embeddings.csv

Output:
- data/processed/similarity/ncrna_similarity_topk.csv
- data/processed/similarity/ncrna_similarity_unique_pairs.csv
- outputs/similarity/ncrna_similarity_report.json
- outputs/similarity/ncrna_similarity_summary.txt
- outputs/similarity/plots/similarity_score_histogram.png
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple

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
        description="Ricerca similarità tra ncRNA tramite embedding."
    )
    parser.add_argument("--config", default="configs/config_public.yaml")
    parser.add_argument("--ncrna-embeddings-csv", default="data/processed/embeddings/ncrna_embeddings.csv")
    parser.add_argument("--output-dir", default="data/processed/similarity")
    parser.add_argument("--report-dir", default="outputs/similarity")
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--unique-top-n", type=int, default=200)
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


def detect_embedding_columns(rows: List[Dict[str, str]]) -> List[str]:
    cols = [c for c in rows[0].keys() if c.startswith("emb_")]
    if not cols:
        raise ValueError("Nessuna colonna embedding trovata: attese colonne emb_001, emb_002, ...")
    return cols


def build_embedding_matrix(rows: List[Dict[str, str]], emb_cols: List[str]) -> np.ndarray:
    matrix = []

    for row in rows:
        vector = [safe_float(row.get(col), 0.0) for col in emb_cols]
        matrix.append(vector)

    x = np.asarray(matrix, dtype=np.float64)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

    return x


def cosine_similarity_matrix(x: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    x_norm = x / norms

    sim = x_norm @ x_norm.T
    sim = np.clip(sim, -1.0, 1.0)

    return sim


def build_topk_rows(
    rows: List[Dict[str, str]],
    sim: np.ndarray,
    top_k: int,
) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    n = sim.shape[0]

    for i in range(n):
        order = np.argsort(-sim[i])

        rank = 0
        for j in order:
            if i == j:
                continue

            rank += 1
            if rank > top_k:
                break

            query = rows[i]
            target = rows[j]

            query_rna_type = query.get("rna_type", "")
            target_rna_type = target.get("rna_type", "")

            out.append(
                {
                    "query_rank": rank,
                    "query_rnacentral_id": query.get("name", query.get("node_id", "")),
                    "query_node_id": query.get("node_id", ""),
                    "query_rna_type": query_rna_type,
                    "query_sequence_length": query.get("sequence_length", ""),
                    "query_gc_frac": query.get("gc_frac", ""),
                    "similar_rnacentral_id": target.get("name", target.get("node_id", "")),
                    "similar_node_id": target.get("node_id", ""),
                    "similar_rna_type": target_rna_type,
                    "similar_sequence_length": target.get("sequence_length", ""),
                    "similar_gc_frac": target.get("gc_frac", ""),
                    "cosine_similarity": round(float(sim[i, j]), 8),
                    "same_rna_type": int(query_rna_type == target_rna_type),
                    "length_difference": round(
                        abs(
                            safe_float(query.get("sequence_length"), 0.0)
                            - safe_float(target.get("sequence_length"), 0.0)
                        ),
                        6,
                    ),
                    "gc_difference": round(
                        abs(
                            safe_float(query.get("gc_frac"), 0.0)
                            - safe_float(target.get("gc_frac"), 0.0)
                        ),
                        6,
                    ),
                }
            )

    return out


def build_unique_pair_rows(
    rows: List[Dict[str, str]],
    sim: np.ndarray,
    unique_top_n: int,
) -> List[Dict[str, Any]]:
    pairs: List[Dict[str, Any]] = []
    n = sim.shape[0]

    for i in range(n):
        for j in range(i + 1, n):
            a = rows[i]
            b = rows[j]

            a_rna_type = a.get("rna_type", "")
            b_rna_type = b.get("rna_type", "")

            pairs.append(
                {
                    "rnacentral_id_a": a.get("name", a.get("node_id", "")),
                    "node_id_a": a.get("node_id", ""),
                    "rna_type_a": a_rna_type,
                    "sequence_length_a": a.get("sequence_length", ""),
                    "gc_frac_a": a.get("gc_frac", ""),
                    "rnacentral_id_b": b.get("name", b.get("node_id", "")),
                    "node_id_b": b.get("node_id", ""),
                    "rna_type_b": b_rna_type,
                    "sequence_length_b": b.get("sequence_length", ""),
                    "gc_frac_b": b.get("gc_frac", ""),
                    "cosine_similarity": round(float(sim[i, j]), 8),
                    "same_rna_type": int(a_rna_type == b_rna_type),
                    "length_difference": round(
                        abs(
                            safe_float(a.get("sequence_length"), 0.0)
                            - safe_float(b.get("sequence_length"), 0.0)
                        ),
                        6,
                    ),
                    "gc_difference": round(
                        abs(
                            safe_float(a.get("gc_frac"), 0.0)
                            - safe_float(b.get("gc_frac"), 0.0)
                        ),
                        6,
                    ),
                }
            )

    pairs.sort(key=lambda row: row["cosine_similarity"], reverse=True)
    return pairs[:unique_top_n]


def summarize_scores(scores: List[float]) -> Dict[str, Any]:
    clean = [v for v in scores if v is not None and not math.isnan(v)]

    if not clean:
        return {
            "count": 0,
            "mean": "NA",
            "min": "NA",
            "median": "NA",
            "max": "NA",
        }

    clean_sorted = sorted(clean)

    return {
        "count": len(clean),
        "mean": round(float(np.mean(clean)), 8),
        "std": round(float(np.std(clean)), 8),
        "min": round(float(np.min(clean)), 8),
        "median": round(float(np.median(clean)), 8),
        "max": round(float(np.max(clean)), 8),
        "q25": round(float(np.quantile(clean_sorted, 0.25)), 8),
        "q75": round(float(np.quantile(clean_sorted, 0.75)), 8),
    }


def save_similarity_histogram(scores: List[float], output_path: Path) -> bool:
    if not MATPLOTLIB_AVAILABLE or plt is None:
        return False

    clean = [v for v in scores if v is not None and not math.isnan(v)]
    if not clean:
        return False

    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 6))
    plt.hist(clean, bins=40)
    plt.title("ncRNA embedding cosine similarity distribution")
    plt.xlabel("Cosine similarity")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()

    return True


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    embeddings_csv = Path(args.ncrna_embeddings_csv)
    output_dir = Path(args.output_dir)
    report_dir = Path(args.report_dir)
    plots_dir = report_dir / "plots"

    output_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_yaml(config_path)

    logging.info("Caricamento embedding ncRNA: %s", embeddings_csv)
    rows = read_csv(embeddings_csv)

    emb_cols = detect_embedding_columns(rows)

    logging.info("Costruzione matrice embedding")
    x = build_embedding_matrix(rows, emb_cols)

    logging.info("Calcolo similarità coseno")
    sim = cosine_similarity_matrix(x)

    logging.info("Costruzione top-k similarity")
    topk_rows = build_topk_rows(
        rows=rows,
        sim=sim,
        top_k=args.top_k,
    )

    logging.info("Costruzione coppie uniche più simili")
    unique_pair_rows = build_unique_pair_rows(
        rows=rows,
        sim=sim,
        unique_top_n=args.unique_top_n,
    )

    topk_csv = output_dir / "ncrna_similarity_topk.csv"
    unique_pairs_csv = output_dir / "ncrna_similarity_unique_pairs.csv"

    topk_fieldnames = [
        "query_rank",
        "query_rnacentral_id",
        "query_node_id",
        "query_rna_type",
        "query_sequence_length",
        "query_gc_frac",
        "similar_rnacentral_id",
        "similar_node_id",
        "similar_rna_type",
        "similar_sequence_length",
        "similar_gc_frac",
        "cosine_similarity",
        "same_rna_type",
        "length_difference",
        "gc_difference",
    ]

    unique_fieldnames = [
        "rnacentral_id_a",
        "node_id_a",
        "rna_type_a",
        "sequence_length_a",
        "gc_frac_a",
        "rnacentral_id_b",
        "node_id_b",
        "rna_type_b",
        "sequence_length_b",
        "gc_frac_b",
        "cosine_similarity",
        "same_rna_type",
        "length_difference",
        "gc_difference",
    ]

    write_csv(topk_rows, topk_csv, topk_fieldnames)
    write_csv(unique_pair_rows, unique_pairs_csv, unique_fieldnames)

    all_pair_scores = []
    n = sim.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            all_pair_scores.append(float(sim[i, j]))

    topk_scores = [float(row["cosine_similarity"]) for row in topk_rows]
    unique_scores = [float(row["cosine_similarity"]) for row in unique_pair_rows]

    rna_type_counter = Counter(row.get("rna_type", "unknown") for row in rows)
    same_type_topk = sum(int(row["same_rna_type"]) for row in topk_rows)
    same_type_unique = sum(int(row["same_rna_type"]) for row in unique_pair_rows)

    plots_created = []
    hist_path = plots_dir / "similarity_score_histogram.png"
    if save_similarity_histogram(all_pair_scores, hist_path):
        plots_created.append("plots/similarity_score_histogram.png")

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "project": cfg.get("project", {}),
        "inputs": {
            "ncrna_embeddings_csv": str(embeddings_csv),
        },
        "outputs": {
            "topk_csv": str(topk_csv),
            "unique_pairs_csv": str(unique_pairs_csv),
            "report_json": str(report_dir / "ncrna_similarity_report.json"),
            "summary_txt": str(report_dir / "ncrna_similarity_summary.txt"),
            "plots_dir": str(plots_dir),
        },
        "summary": {
            "n_ncrna": len(rows),
            "embedding_dim": len(emb_cols),
            "top_k": args.top_k,
            "topk_rows": len(topk_rows),
            "unique_top_n": args.unique_top_n,
            "unique_pair_rows": len(unique_pair_rows),
            "all_pair_score_summary": summarize_scores(all_pair_scores),
            "topk_score_summary": summarize_scores(topk_scores),
            "unique_pair_score_summary": summarize_scores(unique_scores),
            "rna_type_distribution": dict(rna_type_counter),
            "topk_same_rna_type_count": same_type_topk,
            "topk_same_rna_type_fraction": round(
                same_type_topk / len(topk_rows), 6
            ) if topk_rows else 0,
            "unique_same_rna_type_count": same_type_unique,
            "unique_same_rna_type_fraction": round(
                same_type_unique / len(unique_pair_rows), 6
            ) if unique_pair_rows else 0,
            "plots_created": plots_created,
            "methodological_note": (
                "Similarity is computed on unsupervised embeddings. High similarity does not "
                "prove common drought function; it identifies candidates with similar sequence, "
                "metadata and graph-context representations."
            ),
        },
    }

    report_path = report_dir / "ncrna_similarity_report.json"
    with report_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    summary_lines = []
    summary_lines.append("ncRNA SIMILARITY SEARCH SUMMARY")
    summary_lines.append("=" * 80)
    summary_lines.append(f"ncRNA records: {len(rows)}")
    summary_lines.append(f"Embedding dim: {len(emb_cols)}")
    summary_lines.append(f"Top-k per ncRNA: {args.top_k}")
    summary_lines.append(f"Top-k rows: {len(topk_rows)}")
    summary_lines.append(f"Unique pair rows: {len(unique_pair_rows)}")
    summary_lines.append("")
    summary_lines.append(f"RNA type distribution: {dict(rna_type_counter)}")
    summary_lines.append(f"Top-k same RNA type fraction: {report['summary']['topk_same_rna_type_fraction']}")
    summary_lines.append(f"Unique same RNA type fraction: {report['summary']['unique_same_rna_type_fraction']}")
    summary_lines.append("")
    summary_lines.append(f"All pair score summary: {report['summary']['all_pair_score_summary']}")
    summary_lines.append(f"Top-k score summary: {report['summary']['topk_score_summary']}")
    summary_lines.append("")
    summary_lines.append(f"Plots created: {plots_created}")
    summary_lines.append("")
    summary_lines.append("Methodological note:")
    summary_lines.append(report["summary"]["methodological_note"])

    summary_path = report_dir / "ncrna_similarity_summary.txt"
    with summary_path.open("w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")

    print("\n".join(summary_lines))

    logging.info("Top-k similarity CSV: %s", topk_csv)
    logging.info("Unique pair CSV: %s", unique_pairs_csv)
    logging.info("Report: %s", report_path)
    logging.info("Summary: %s", summary_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
