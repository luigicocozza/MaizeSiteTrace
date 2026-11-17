#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()
OUT_DIR = BASE_DIR / "outputs/prjna419326"

RPM_MATRIX = OUT_DIR / f"{PROJECT}_mirnaseq_srna_rpm_matrix_filtered.tsv"

ALL_SCORES_TSV = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_scores.tsv"
CANDIDATES_TSV = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_candidates.tsv"
TOP500_TSV = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_top500.tsv"
TOP500_FASTA = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_top500.fasta"
SUMMARY_TXT = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_summary.txt"
REPORT_JSON = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_report.json"

PSEUDO_RPM = 0.1
MIN_MAX_MEAN_RPM = 5.0
MIN_ABS_LOG2FC = 1.0
TOP_N = 500


def mean(values):
    values = list(values)
    return sum(values) / len(values) if values else 0.0


def log2fc(ds_values, ww_values):
    return math.log2((mean(ds_values) + PSEUDO_RPM) / (mean(ww_values) + PSEUDO_RPM))


def direction(v):
    if v >= MIN_ABS_LOG2FC:
        return "up"
    if v <= -MIN_ABS_LOG2FC:
        return "down"
    return "stable"


def candidate_class(b73_fc, h21_fc):
    b73_dir = direction(b73_fc)
    h21_dir = direction(h21_fc)

    if b73_dir == "up" and h21_dir == "up":
        return "conserved_up"
    if b73_dir == "down" and h21_dir == "down":
        return "conserved_down"
    if b73_dir in {"up", "down"} and h21_dir in {"up", "down"} and b73_dir != h21_dir:
        return "opposite_response"
    if b73_dir != "stable" and h21_dir == "stable":
        return f"B73_{b73_dir}_only"
    if h21_dir != "stable" and b73_dir == "stable":
        return f"H21_{h21_dir}_only"
    return "not_candidate"


def stress_score(total_count, max_mean_rpm, b73_fc, h21_fc):
    abs_b73 = abs(b73_fc)
    abs_h21 = abs(h21_fc)

    same_direction_bonus = 0.0
    if b73_fc * h21_fc > 0:
        same_direction_bonus = min(abs_b73, abs_h21)

    return round(
        (max(abs_b73, abs_h21) + 0.5 * same_direction_bonus)
        * math.log10(total_count + 10)
        * math.log10(max_mean_rpm + 10),
        6,
    )


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"SCREENING small RNA drought-response - {PROJECT}")
    print("=" * 80)
    print(f"Input RPM matrix: {RPM_MATRIX}")
    print("")

    if not RPM_MATRIX.exists():
        raise FileNotFoundError(f"Manca matrice RPM: {RPM_MATRIX}")

    class_counts = {}
    total_features = 0
    candidate_count = 0
    top_records = []

    with RPM_MATRIX.open("r", encoding="utf-8") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        sample_cols = reader.fieldnames[5:]

        b73_ds = [c for c in sample_cols if c.startswith("B73_DS")]
        b73_ww = [c for c in sample_cols if c.startswith("B73_WW")]
        h21_ds = [c for c in sample_cols if c.startswith("H21_DS")]
        h21_ww = [c for c in sample_cols if c.startswith("H21_WW")]

        print("Gruppi individuati:")
        print("B73 drought:", b73_ds)
        print("B73 control:", b73_ww)
        print("H21 drought:", h21_ds)
        print("H21 control:", h21_ww)
        print("")

        if not (len(b73_ds) == len(b73_ww) == len(h21_ds) == len(h21_ww) == 3):
            raise RuntimeError("Gruppi non bilanciati: attese 3 repliche per ciascun gruppo.")

        fieldnames = [
            "feature_id",
            "sequence",
            "length",
            "total_count",
            "samples_detected",
            "mean_rpm_B73_DS",
            "mean_rpm_B73_WW",
            "log2FC_B73_DS_vs_WW",
            "mean_rpm_H21_DS",
            "mean_rpm_H21_WW",
            "log2FC_H21_DS_vs_WW",
            "max_mean_rpm",
            "candidate_class",
            "stress_score",
        ]

        with ALL_SCORES_TSV.open("w", newline="", encoding="utf-8") as fall, \
             CANDIDATES_TSV.open("w", newline="", encoding="utf-8") as fcand:

            writer_all = csv.DictWriter(fall, fieldnames=fieldnames, delimiter="\t")
            writer_cand = csv.DictWriter(fcand, fieldnames=fieldnames, delimiter="\t")
            writer_all.writeheader()
            writer_cand.writeheader()

            for row in reader:
                total_features += 1

                feature_id = row["feature_id"]
                seq = row["sequence"]
                length = int(row["length"])
                total_count = int(row["total_count"])
                samples_detected = int(row["samples_detected"])

                vals = {c: float(row[c]) for c in sample_cols}

                b73_ds_vals = [vals[c] for c in b73_ds]
                b73_ww_vals = [vals[c] for c in b73_ww]
                h21_ds_vals = [vals[c] for c in h21_ds]
                h21_ww_vals = [vals[c] for c in h21_ww]

                mean_b73_ds = mean(b73_ds_vals)
                mean_b73_ww = mean(b73_ww_vals)
                mean_h21_ds = mean(h21_ds_vals)
                mean_h21_ww = mean(h21_ww_vals)

                fc_b73 = log2fc(b73_ds_vals, b73_ww_vals)
                fc_h21 = log2fc(h21_ds_vals, h21_ww_vals)

                max_mean = max(mean_b73_ds, mean_b73_ww, mean_h21_ds, mean_h21_ww)
                cls = candidate_class(fc_b73, fc_h21)
                score = stress_score(total_count, max_mean, fc_b73, fc_h21)

                out = {
                    "feature_id": feature_id,
                    "sequence": seq,
                    "length": length,
                    "total_count": total_count,
                    "samples_detected": samples_detected,
                    "mean_rpm_B73_DS": round(mean_b73_ds, 6),
                    "mean_rpm_B73_WW": round(mean_b73_ww, 6),
                    "log2FC_B73_DS_vs_WW": round(fc_b73, 6),
                    "mean_rpm_H21_DS": round(mean_h21_ds, 6),
                    "mean_rpm_H21_WW": round(mean_h21_ww, 6),
                    "log2FC_H21_DS_vs_WW": round(fc_h21, 6),
                    "max_mean_rpm": round(max_mean, 6),
                    "candidate_class": cls,
                    "stress_score": score,
                }

                writer_all.writerow(out)

                is_candidate = (
                    cls != "not_candidate"
                    and max_mean >= MIN_MAX_MEAN_RPM
                    and 18 <= length <= 30
                )

                class_counts[cls] = class_counts.get(cls, 0) + 1

                if is_candidate:
                    candidate_count += 1
                    writer_cand.writerow(out)
                    top_records.append(out)

                if total_features % 100000 == 0:
                    print(f"Processate {total_features} feature...")

    top_records.sort(key=lambda r: float(r["stress_score"]), reverse=True)
    top_records = top_records[:TOP_N]

    with TOP500_TSV.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(top_records[0].keys()) if top_records else []
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        if fieldnames:
            writer.writeheader()
            writer.writerows(top_records)

    with TOP500_FASTA.open("w", encoding="utf-8") as f:
        for r in top_records:
            header = (
                f">{r['feature_id']}|class={r['candidate_class']}|"
                f"score={r['stress_score']}|"
                f"B73_log2FC={r['log2FC_B73_DS_vs_WW']}|"
                f"H21_log2FC={r['log2FC_H21_DS_vs_WW']}|"
                f"len={r['length']}|total={r['total_count']}"
            )
            f.write(header + "\n")
            f.write(r["sequence"] + "\n")

    report = {
        "project": PROJECT,
        "created": datetime.now().isoformat(timespec="seconds"),
        "input_rpm_matrix": str(RPM_MATRIX),
        "thresholds": {
            "pseudo_rpm": PSEUDO_RPM,
            "min_max_mean_rpm": MIN_MAX_MEAN_RPM,
            "min_abs_log2fc": MIN_ABS_LOG2FC,
        },
        "total_features_screened": total_features,
        "candidate_count": candidate_count,
        "class_counts": class_counts,
        "outputs": {
            "all_scores_tsv": str(ALL_SCORES_TSV),
            "candidates_tsv": str(CANDIDATES_TSV),
            "top500_tsv": str(TOP500_TSV),
            "top500_fasta": str(TOP500_FASTA),
            "summary_txt": str(SUMMARY_TXT),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append(f"SCREENING small RNA drought-response - {PROJECT}")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Feature analizzate: {total_features}")
    lines.append(f"Candidati preliminari: {candidate_count}")
    lines.append(f"Soglie: max_mean_RPM >= {MIN_MAX_MEAN_RPM}, abs(log2FC) >= {MIN_ABS_LOG2FC}")
    lines.append("")
    lines.append("Classi:")
    for k in sorted(class_counts):
        lines.append(f"{k}: {class_counts[k]}")
    lines.append("")
    lines.append(f"Top {TOP_N} candidati salvati in FASTA e TSV.")
    lines.append("")
    lines.append("Output:")
    lines.append(str(ALL_SCORES_TSV))
    lines.append(str(CANDIDATES_TSV))
    lines.append(str(TOP500_TSV))
    lines.append(str(TOP500_FASTA))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("")
    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
