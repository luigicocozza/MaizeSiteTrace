#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy import stats
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

SAMPLE_SHEET = BASE_DIR / "outputs/prjna419326/PRJNA419326_rnaseq_sample_sheet.tsv"

GENE_TPM = BASE_DIR / "data/processed/rnaseq_matrices/PRJNA419326_rnaseq_gene_TPM_matrix.tsv"
GENE_COUNTS = BASE_DIR / "data/processed/rnaseq_matrices/PRJNA419326_rnaseq_gene_NumReads_matrix.tsv"

OUT_DIR = BASE_DIR / "outputs/prjna419326"
OUT_DIR.mkdir(parents=True, exist_ok=True)

QC_TSV = OUT_DIR / "PRJNA419326_rnaseq_gene_matrix_qc.tsv"

B73_DE = OUT_DIR / "PRJNA419326_rnaseq_gene_DE_B73_DS_vs_WW.tsv"
H21_DE = OUT_DIR / "PRJNA419326_rnaseq_gene_DE_H21_DS_vs_WW.tsv"
COMBINED_DE = OUT_DIR / "PRJNA419326_rnaseq_gene_DE_combined_B73_H21.tsv"
TOP_CANDIDATES = OUT_DIR / "PRJNA419326_rnaseq_gene_drought_candidates_top1000.tsv"

SUMMARY_TXT = OUT_DIR / "PRJNA419326_rnaseq_gene_DE_screening_summary.txt"
REPORT_JSON = OUT_DIR / "PRJNA419326_rnaseq_gene_DE_screening_report.json"


MIN_TOTAL_COUNTS = 10.0
MIN_TPM = 1.0
MIN_SAMPLES_WITH_TPM = 3
LOG2FC_THRESHOLD = 1.0
PVALUE_THRESHOLD = 0.05
QVALUE_THRESHOLD = 0.10


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    q = np.ones(n, dtype=float)

    valid = np.isfinite(p)
    if valid.sum() == 0:
        return q

    p_valid = p[valid]
    order = np.argsort(p_valid)
    ranked = p_valid[order]

    q_ranked = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_ranked = np.clip(q_ranked, 0, 1)

    q_valid = np.empty_like(p_valid)
    q_valid[order] = q_ranked

    q[valid] = q_valid
    return q


def safe_ttest(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if len(a) < 2 or len(b) < 2:
        return 1.0

    if np.allclose(a, a[0]) and np.allclose(b, b[0]) and np.isclose(a[0], b[0]):
        return 1.0

    if not SCIPY_OK:
        return float("nan")

    try:
        res = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
        p = float(res.pvalue)
        if not math.isfinite(p):
            return 1.0
        return p
    except Exception:
        return 1.0


def classify_response(log2fc: float, pvalue: float, qvalue: float) -> tuple[str, str, str]:
    if log2fc >= LOG2FC_THRESHOLD and pvalue <= PVALUE_THRESHOLD:
        nominal = "up"
    elif log2fc <= -LOG2FC_THRESHOLD and pvalue <= PVALUE_THRESHOLD:
        nominal = "down"
    else:
        nominal = "stable"

    if log2fc >= LOG2FC_THRESHOLD and qvalue <= QVALUE_THRESHOLD:
        fdr = "up"
    elif log2fc <= -LOG2FC_THRESHOLD and qvalue <= QVALUE_THRESHOLD:
        fdr = "down"
    else:
        fdr = "stable"

    # Classe utile per generare label biologiche in modo non troppo restrittivo.
    if log2fc >= LOG2FC_THRESHOLD and pvalue <= 0.10:
        screen = "up"
    elif log2fc <= -LOG2FC_THRESHOLD and pvalue <= 0.10:
        screen = "down"
    else:
        screen = "stable"

    return nominal, fdr, screen


def run_de_for_genotype(genotype: str, meta: pd.DataFrame, tpm: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    ds_samples = meta.loc[(meta["genotype"] == genotype) & (meta["condition"] == "DS"), "sample_alias"].tolist()
    ww_samples = meta.loc[(meta["genotype"] == genotype) & (meta["condition"] == "WW"), "sample_alias"].tolist()

    missing = [s for s in ds_samples + ww_samples if s not in tpm.columns]
    if missing:
        raise SystemExit(f"ERRORE: campioni mancanti nella matrice per {genotype}: {missing}")

    rows = []

    all_samples = ds_samples + ww_samples

    for gene_id in tpm.index:
        tpm_ds = tpm.loc[gene_id, ds_samples].astype(float).to_numpy()
        tpm_ww = tpm.loc[gene_id, ww_samples].astype(float).to_numpy()

        count_ds = counts.loc[gene_id, ds_samples].astype(float).to_numpy()
        count_ww = counts.loc[gene_id, ww_samples].astype(float).to_numpy()

        total_counts = float(np.sum(count_ds) + np.sum(count_ww))
        samples_with_tpm = int(np.sum(tpm.loc[gene_id, all_samples].astype(float).to_numpy() >= MIN_TPM))

        expressed = total_counts >= MIN_TOTAL_COUNTS and samples_with_tpm >= MIN_SAMPLES_WITH_TPM

        mean_tpm_ds = float(np.mean(tpm_ds))
        mean_tpm_ww = float(np.mean(tpm_ww))
        mean_counts_ds = float(np.mean(count_ds))
        mean_counts_ww = float(np.mean(count_ww))

        log2fc_tpm = float(np.log2((mean_tpm_ds + 0.1) / (mean_tpm_ww + 0.1)))
        log2fc_counts = float(np.log2((mean_counts_ds + 1.0) / (mean_counts_ww + 1.0)))

        pvalue = safe_ttest(np.log2(tpm_ds + 1.0), np.log2(tpm_ww + 1.0))

        rows.append(
            {
                "gene_id": gene_id,
                "genotype": genotype,
                "comparison": f"{genotype}_DS_vs_WW",
                "expressed_filter": expressed,
                "total_counts": total_counts,
                "samples_with_TPM_ge_1": samples_with_tpm,
                "mean_TPM_DS": mean_tpm_ds,
                "mean_TPM_WW": mean_tpm_ww,
                "mean_counts_DS": mean_counts_ds,
                "mean_counts_WW": mean_counts_ww,
                "log2FC_TPM_DS_vs_WW": log2fc_tpm,
                "log2FC_counts_DS_vs_WW": log2fc_counts,
                "pvalue_welch_log2TPM": pvalue,
            }
        )

    de = pd.DataFrame(rows)

    # FDR solo sui geni che passano il filtro di espressione
    de["qvalue_BH"] = 1.0
    mask = de["expressed_filter"].astype(bool)
    de.loc[mask, "qvalue_BH"] = bh_fdr(de.loc[mask, "pvalue_welch_log2TPM"].to_numpy())

    classes = de.apply(
        lambda r: classify_response(
            float(r["log2FC_TPM_DS_vs_WW"]),
            float(r["pvalue_welch_log2TPM"]) if pd.notna(r["pvalue_welch_log2TPM"]) else 1.0,
            float(r["qvalue_BH"]) if pd.notna(r["qvalue_BH"]) else 1.0,
        ),
        axis=1,
    )

    de["nominal_response"] = [x[0] for x in classes]
    de["fdr_response"] = [x[1] for x in classes]
    de["screen_response"] = [x[2] for x in classes]

    de["minus_log10_pvalue"] = -np.log10(de["pvalue_welch_log2TPM"].replace(0, 1e-300).fillna(1.0))
    de["de_priority_score"] = (
        de["expressed_filter"].astype(int)
        * de["log2FC_TPM_DS_vs_WW"].abs()
        * de["minus_log10_pvalue"]
        * np.log10(de["total_counts"] + 10.0)
    )

    de = de.sort_values(
        by=["de_priority_score", "expressed_filter", "minus_log10_pvalue", "log2FC_TPM_DS_vs_WW"],
        ascending=[False, False, False, False],
    )

    return de


def combined_class(b73: str, h21: str) -> str:
    if b73 == "up" and h21 == "up":
        return "conserved_up"
    if b73 == "down" and h21 == "down":
        return "conserved_down"
    if b73 == "up" and h21 == "down":
        return "opposite_B73_up_H21_down"
    if b73 == "down" and h21 == "up":
        return "opposite_B73_down_H21_up"
    if b73 == "up" and h21 == "stable":
        return "B73_up_only"
    if b73 == "down" and h21 == "stable":
        return "B73_down_only"
    if b73 == "stable" and h21 == "up":
        return "H21_up_only"
    if b73 == "stable" and h21 == "down":
        return "H21_down_only"
    return "stable_or_not_detected"


def main() -> int:
    meta = pd.read_csv(SAMPLE_SHEET, sep="\t")
    tpm = pd.read_csv(GENE_TPM, sep="\t", index_col=0)
    counts = pd.read_csv(GENE_COUNTS, sep="\t", index_col=0)

    # QC campioni
    qc_rows = []
    for sample in meta["sample_alias"].tolist():
        if sample not in tpm.columns:
            continue
        qc_rows.append(
            {
                "sample_alias": sample,
                "genotype": meta.loc[meta["sample_alias"] == sample, "genotype"].iloc[0],
                "condition": meta.loc[meta["sample_alias"] == sample, "condition"].iloc[0],
                "replicate": meta.loc[meta["sample_alias"] == sample, "replicate"].iloc[0],
                "total_NumReads": float(counts[sample].sum()),
                "total_TPM": float(tpm[sample].sum()),
                "genes_TPM_ge_1": int((tpm[sample] >= 1.0).sum()),
                "genes_counts_ge_10": int((counts[sample] >= 10.0).sum()),
            }
        )

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(QC_TSV, sep="\t", index=False)

    b73 = run_de_for_genotype("B73", meta, tpm, counts)
    h21 = run_de_for_genotype("H21", meta, tpm, counts)

    b73.to_csv(B73_DE, sep="\t", index=False)
    h21.to_csv(H21_DE, sep="\t", index=False)

    b73_small = b73[
        [
            "gene_id",
            "expressed_filter",
            "mean_TPM_DS",
            "mean_TPM_WW",
            "log2FC_TPM_DS_vs_WW",
            "pvalue_welch_log2TPM",
            "qvalue_BH",
            "nominal_response",
            "fdr_response",
            "screen_response",
            "de_priority_score",
        ]
    ].rename(
        columns={
            "expressed_filter": "B73_expressed_filter",
            "mean_TPM_DS": "B73_mean_TPM_DS",
            "mean_TPM_WW": "B73_mean_TPM_WW",
            "log2FC_TPM_DS_vs_WW": "B73_log2FC_TPM_DS_vs_WW",
            "pvalue_welch_log2TPM": "B73_pvalue",
            "qvalue_BH": "B73_qvalue",
            "nominal_response": "B73_nominal_response",
            "fdr_response": "B73_fdr_response",
            "screen_response": "B73_screen_response",
            "de_priority_score": "B73_de_priority_score",
        }
    )

    h21_small = h21[
        [
            "gene_id",
            "expressed_filter",
            "mean_TPM_DS",
            "mean_TPM_WW",
            "log2FC_TPM_DS_vs_WW",
            "pvalue_welch_log2TPM",
            "qvalue_BH",
            "nominal_response",
            "fdr_response",
            "screen_response",
            "de_priority_score",
        ]
    ].rename(
        columns={
            "expressed_filter": "H21_expressed_filter",
            "mean_TPM_DS": "H21_mean_TPM_DS",
            "mean_TPM_WW": "H21_mean_TPM_WW",
            "log2FC_TPM_DS_vs_WW": "H21_log2FC_TPM_DS_vs_WW",
            "pvalue_welch_log2TPM": "H21_pvalue",
            "qvalue_BH": "H21_qvalue",
            "nominal_response": "H21_nominal_response",
            "fdr_response": "H21_fdr_response",
            "screen_response": "H21_screen_response",
            "de_priority_score": "H21_de_priority_score",
        }
    )

    combined = b73_small.merge(h21_small, on="gene_id", how="outer").fillna(
        {
            "B73_screen_response": "stable",
            "H21_screen_response": "stable",
            "B73_nominal_response": "stable",
            "H21_nominal_response": "stable",
            "B73_fdr_response": "stable",
            "H21_fdr_response": "stable",
            "B73_de_priority_score": 0.0,
            "H21_de_priority_score": 0.0,
        }
    )

    combined["combined_screen_class"] = combined.apply(
        lambda r: combined_class(str(r["B73_screen_response"]), str(r["H21_screen_response"])),
        axis=1,
    )

    combined["combined_de_priority_score"] = (
        combined["B73_de_priority_score"].astype(float).fillna(0.0)
        + combined["H21_de_priority_score"].astype(float).fillna(0.0)
    )

    combined = combined.sort_values(
        by=["combined_de_priority_score"],
        ascending=False,
    )

    combined.to_csv(COMBINED_DE, sep="\t", index=False)

    top = combined.head(1000).copy()
    top.to_csv(TOP_CANDIDATES, sep="\t", index=False)

    b73_counts = b73["screen_response"].value_counts().to_dict()
    h21_counts = h21["screen_response"].value_counts().to_dict()
    combined_counts = combined["combined_screen_class"].value_counts().to_dict()

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "scipy_ok": SCIPY_OK,
        "inputs": {
            "sample_sheet": str(SAMPLE_SHEET),
            "gene_tpm": str(GENE_TPM),
            "gene_counts": str(GENE_COUNTS),
        },
        "parameters": {
            "min_total_counts": MIN_TOTAL_COUNTS,
            "min_tpm": MIN_TPM,
            "min_samples_with_tpm": MIN_SAMPLES_WITH_TPM,
            "log2fc_threshold": LOG2FC_THRESHOLD,
            "pvalue_threshold": PVALUE_THRESHOLD,
            "qvalue_threshold": QVALUE_THRESHOLD,
        },
        "matrix_shape": {
            "genes": int(tpm.shape[0]),
            "samples": int(tpm.shape[1]),
        },
        "b73_screen_counts": b73_counts,
        "h21_screen_counts": h21_counts,
        "combined_screen_counts": combined_counts,
        "outputs": {
            "qc_tsv": str(QC_TSV),
            "b73_de": str(B73_DE),
            "h21_de": str(H21_DE),
            "combined_de": str(COMBINED_DE),
            "top_candidates": str(TOP_CANDIDATES),
            "summary_txt": str(SUMMARY_TXT),
            "report_json": str(REPORT_JSON),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("SCREENING DIFFERENZIALE RNA-SEQ GENE-LEVEL PRJNA419326")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"SciPy disponibile: {SCIPY_OK}")
    lines.append("")
    lines.append(f"Geni in matrice: {tpm.shape[0]}")
    lines.append(f"Campioni in matrice: {tpm.shape[1]}")
    lines.append("")
    lines.append("Parametri:")
    lines.append(f"- min_total_counts: {MIN_TOTAL_COUNTS}")
    lines.append(f"- min_TPM: {MIN_TPM}")
    lines.append(f"- min_samples_with_TPM: {MIN_SAMPLES_WITH_TPM}")
    lines.append(f"- soglia |log2FC|: {LOG2FC_THRESHOLD}")
    lines.append(f"- soglia p-value nominale: {PVALUE_THRESHOLD}")
    lines.append(f"- soglia q-value FDR: {QVALUE_THRESHOLD}")
    lines.append("")
    lines.append("Conteggi screen_response B73:")
    for k, v in b73_counts.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("Conteggi screen_response H21:")
    for k, v in h21_counts.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("Classi combinate B73/H21:")
    for k, v in combined_counts.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("Top 20 geni candidati:")
    show_cols = [
        "gene_id",
        "combined_screen_class",
        "B73_log2FC_TPM_DS_vs_WW",
        "B73_pvalue",
        "H21_log2FC_TPM_DS_vs_WW",
        "H21_pvalue",
        "combined_de_priority_score",
    ]
    for _, r in combined.head(20)[show_cols].iterrows():
        lines.append(
            f"- {r['gene_id']} | {r['combined_screen_class']} | "
            f"B73_log2FC={float(r['B73_log2FC_TPM_DS_vs_WW']):.3f}, "
            f"B73_p={float(r['B73_pvalue']):.3g}, "
            f"H21_log2FC={float(r['H21_log2FC_TPM_DS_vs_WW']):.3f}, "
            f"H21_p={float(r['H21_pvalue']):.3g}, "
            f"score={float(r['combined_de_priority_score']):.3f}"
        )
    lines.append("")
    lines.append("Output:")
    lines.append(str(QC_TSV))
    lines.append(str(B73_DE))
    lines.append(str(H21_DE))
    lines.append(str(COMBINED_DE))
    lines.append(str(TOP_CANDIDATES))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
