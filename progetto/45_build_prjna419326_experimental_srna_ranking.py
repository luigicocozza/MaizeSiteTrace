#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()
OUT_DIR = BASE_DIR / "outputs/prjna419326"

TOP500_TSV = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_top500.tsv"
CANDIDATES_TSV = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_candidates.tsv"
UNMATCHED_TSV = OUT_DIR / f"{PROJECT}_top500_srna_unmatched.tsv"

RANKING_TSV = OUT_DIR / f"{PROJECT}_experimental_srna_v3_ranking.tsv"
RANKING_CSV = OUT_DIR / f"{PROJECT}_experimental_srna_v3_ranking.csv"
RANKING_FASTA = OUT_DIR / f"{PROJECT}_experimental_srna_v3_ranking_top500.fasta"
SUMMARY_TXT = OUT_DIR / f"{PROJECT}_experimental_srna_v3_ranking_summary.txt"
REPORT_JSON = OUT_DIR / f"{PROJECT}_experimental_srna_v3_ranking_report.json"


CLASS_PRIORITY = {
    "conserved_up": 100,
    "conserved_down": 95,
    "opposite_response": 90,
    "B73_up_only": 80,
    "B73_down_only": 78,
    "H21_up_only": 76,
    "H21_down_only": 74,
    "not_candidate": 0,
}


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def compute_v3_priority(row: dict) -> float:
    cls = row["candidate_class"]
    class_score = CLASS_PRIORITY.get(cls, 0)

    stress_score = float(row["stress_score"])
    total_count = int(row["total_count"])
    samples_detected = int(row["samples_detected"])
    length = int(row["length"])

    b73_fc = abs(float(row["log2FC_B73_DS_vs_WW"]))
    h21_fc = abs(float(row["log2FC_H21_DS_vs_WW"]))

    length_bonus = 0
    if 20 <= length <= 24:
        length_bonus = 10
    elif 21 <= length <= 24:
        length_bonus = 8
    elif 18 <= length <= 30:
        length_bonus = 5

    replication_bonus = min(samples_detected, 12) * 1.5

    abundance_bonus = 0
    if total_count >= 100000:
        abundance_bonus = 15
    elif total_count >= 10000:
        abundance_bonus = 10
    elif total_count >= 1000:
        abundance_bonus = 5

    fc_bonus = min(max(b73_fc, h21_fc), 10) * 2

    return round(class_score + length_bonus + replication_bonus + abundance_bonus + fc_bonus + stress_score, 6)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    top = read_tsv(TOP500_TSV)

    print("BUILD V3 EXPERIMENTAL small RNA RANKING")
    print("=" * 80)
    print(f"Input top500: {TOP500_TSV}")
    print(f"Record input: {len(top)}")
    print("")

    ranked = []

    for row in top:
        seq = row["sequence"]
        feature_id = row["feature_id"]

        out = {
            "v3_rank_id": "",
            "feature_id": feature_id,
            "sequence": seq,
            "length": row["length"],
            "candidate_class": row["candidate_class"],
            "v3_priority_score": compute_v3_priority(row),
            "stress_score": row["stress_score"],
            "total_count": row["total_count"],
            "samples_detected": row["samples_detected"],
            "mean_rpm_B73_DS": row["mean_rpm_B73_DS"],
            "mean_rpm_B73_WW": row["mean_rpm_B73_WW"],
            "log2FC_B73_DS_vs_WW": row["log2FC_B73_DS_vs_WW"],
            "mean_rpm_H21_DS": row["mean_rpm_H21_DS"],
            "mean_rpm_H21_WW": row["mean_rpm_H21_WW"],
            "log2FC_H21_DS_vs_WW": row["log2FC_H21_DS_vs_WW"],
            "rnacentral_exact_match": "no",
            "interpretation": "",
        }

        cls = row["candidate_class"]

        if cls == "conserved_up":
            out["interpretation"] = "small RNA up-regolata in risposta a siccità in B73 e H21"
        elif cls == "conserved_down":
            out["interpretation"] = "small RNA down-regolata in risposta a siccità in B73 e H21"
        elif cls == "opposite_response":
            out["interpretation"] = "risposta opposta tra genotipi, potenziale marcatore genotype-specific"
        elif cls.startswith("B73"):
            out["interpretation"] = "risposta alla siccità specifica del genotipo B73"
        elif cls.startswith("H21"):
            out["interpretation"] = "risposta alla siccità specifica del genotipo H21"
        else:
            out["interpretation"] = "candidato preliminare da verificare"

        ranked.append(out)

    ranked.sort(key=lambda r: float(r["v3_priority_score"]), reverse=True)

    for i, row in enumerate(ranked, start=1):
        row["v3_rank_id"] = f"PRJNA419326_sRNA_rank_{i:04d}"

    fieldnames = [
        "v3_rank_id",
        "feature_id",
        "sequence",
        "length",
        "candidate_class",
        "v3_priority_score",
        "stress_score",
        "total_count",
        "samples_detected",
        "mean_rpm_B73_DS",
        "mean_rpm_B73_WW",
        "log2FC_B73_DS_vs_WW",
        "mean_rpm_H21_DS",
        "mean_rpm_H21_WW",
        "log2FC_H21_DS_vs_WW",
        "rnacentral_exact_match",
        "interpretation",
    ]

    with RANKING_TSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(ranked)

    with RANKING_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ranked)

    with RANKING_FASTA.open("w", encoding="utf-8") as f:
        for row in ranked:
            header = (
                f">{row['v3_rank_id']}|feature={row['feature_id']}|"
                f"class={row['candidate_class']}|"
                f"v3_score={row['v3_priority_score']}|"
                f"B73_log2FC={row['log2FC_B73_DS_vs_WW']}|"
                f"H21_log2FC={row['log2FC_H21_DS_vs_WW']}|"
                f"len={row['length']}|total={row['total_count']}"
            )
            f.write(header + "\n")
            f.write(row["sequence"] + "\n")

    class_counts = {}
    for row in ranked:
        cls = row["candidate_class"]
        class_counts[cls] = class_counts.get(cls, 0) + 1

    report = {
        "project": PROJECT,
        "created": datetime.now().isoformat(timespec="seconds"),
        "input_top500": str(TOP500_TSV),
        "n_ranked": len(ranked),
        "class_counts": class_counts,
        "note": "Ranking sperimentale V3 basato su small RNA reali PRJNA419326. Nessun match esatto RNAcentral top500 trovato nella fase precedente.",
        "outputs": {
            "ranking_tsv": str(RANKING_TSV),
            "ranking_csv": str(RANKING_CSV),
            "ranking_fasta": str(RANKING_FASTA),
            "summary_txt": str(SUMMARY_TXT),
            "report_json": str(REPORT_JSON),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("RANKING V3 small RNA sperimentali PRJNA419326")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Candidati ranked: {len(ranked)}")
    lines.append("")
    lines.append("Classi nel ranking:")
    for k in sorted(class_counts):
        lines.append(f"{k}: {class_counts[k]}")
    lines.append("")
    lines.append("Top 20 candidati:")
    lines.append(
        "rank | v3_rank_id | feature_id | class | len | total_count | "
        "B73_log2FC | H21_log2FC | v3_priority_score"
    )

    for i, row in enumerate(ranked[:20], start=1):
        lines.append(
            f"{i} | {row['v3_rank_id']} | {row['feature_id']} | "
            f"{row['candidate_class']} | {row['length']} | {row['total_count']} | "
            f"{row['log2FC_B73_DS_vs_WW']} | {row['log2FC_H21_DS_vs_WW']} | "
            f"{row['v3_priority_score']}"
        )

    lines.append("")
    lines.append("Output:")
    lines.append(str(RANKING_TSV))
    lines.append(str(RANKING_CSV))
    lines.append(str(RANKING_FASTA))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
