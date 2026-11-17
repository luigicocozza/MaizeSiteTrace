#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd


BASE_DIR = Path(".").resolve()
OUT_DIR = BASE_DIR / "outputs/prjna419326"
FINAL_DIR = OUT_DIR / "final_v3_report"

FINAL_DIR.mkdir(parents=True, exist_ok=True)

ANTILEAK_PRED = OUT_DIR / "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
FULL_PRED = OUT_DIR / "PRJNA419326_srna_supervised_ML_predictions.tsv"
ANTILEAK_COMP = OUT_DIR / "PRJNA419326_srna_supervised_ML_antileakage_model_comparison.tsv"
FULL_COMP = OUT_DIR / "PRJNA419326_srna_supervised_ML_model_comparison.tsv"

LABELS_REPORT = OUT_DIR / "PRJNA419326_srna_supervised_labels_report.json"
DE_REPORT = OUT_DIR / "PRJNA419326_rnaseq_gene_DE_screening_report.json"
QUANT_REPORT = OUT_DIR / "PRJNA419326_rnaseq_salmon_quant_report.json"
DOWNLOAD_REPORT = OUT_DIR / "PRJNA419326_rnaseq_download_VALIDATO_report.json"
REF_REPORT = OUT_DIR / "PRJNA419326_rnaseq_reference_prepare_report.json"
ANTILEAK_REPORT = OUT_DIR / "PRJNA419326_srna_supervised_ML_antileakage_report.json"
FULL_ML_REPORT = OUT_DIR / "PRJNA419326_srna_supervised_ML_report.json"

GAE_REPORT = BASE_DIR / "outputs/v3_gae/v3_gae_metrics.json"

FINAL_TABLE = FINAL_DIR / "PRJNA419326_V3_FINAL_srna_candidate_priority_table.tsv"
FINAL_TOP100 = FINAL_DIR / "PRJNA419326_V3_FINAL_top100_srna_candidates.tsv"
FINAL_REPORT_MD = FINAL_DIR / "PRJNA419326_V3_FINAL_scientific_report.md"
FINAL_REPORT_TXT = FINAL_DIR / "PRJNA419326_V3_FINAL_scientific_report.txt"
FINAL_REPORT_JSON = FINAL_DIR / "PRJNA419326_V3_FINAL_report_index.json"


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def safe_col(df: pd.DataFrame, col: str, default=None):
    if col in df.columns:
        return df[col]
    return default


def assign_final_tier(row) -> str:
    anti_p = float(row.get("antileakage_probability_positive", 0.0))
    full_p = float(row.get("ml_predicted_probability_positive", 0.0))
    silver_y = int(row.get("ml_label_y", -1))
    label_status = str(row.get("label_status", ""))

    if anti_p >= 0.90 and silver_y == 1:
        return "TIER_1_high_confidence_supported"
    if anti_p >= 0.90 and silver_y == -1:
        return "TIER_2_model_prioritized_unlabeled"
    if anti_p >= 0.75 and silver_y == 1:
        return "TIER_3_supported_moderate_model"
    if anti_p >= 0.75 and silver_y == -1:
        return "TIER_4_unlabeled_moderate_model"
    if silver_y == 1 and full_p >= 0.90:
        return "TIER_5_rule_supported_full_model_only"
    return "TIER_6_low_priority_or_negative"


def main() -> int:
    anti = pd.read_csv(ANTILEAK_PRED, sep="\t")
    full = pd.read_csv(FULL_PRED, sep="\t")

    keep_full = [
        "v3_rank_id",
        "ml_predicted_probability_positive",
        "ml_predicted_label",
        "ml_model_best",
    ]
    full_small = full[[c for c in keep_full if c in full.columns]].copy()

    merged = anti.merge(full_small, on="v3_rank_id", how="left")

    merged["final_tier"] = merged.apply(assign_final_tier, axis=1)

    # Score finale composito: privilegia modello anti-leakage, poi evidenza biologica e priorità iniziale.
    merged["final_priority_score"] = (
        100.0 * merged["antileakage_probability_positive"].astype(float).fillna(0.0)
        + 0.10 * merged["v3_priority_score"].astype(float).fillna(0.0)
        + 2.0 * merged["max_target_support_score"].astype(float).fillna(0.0)
        + 1.0 * merged["total_opposite_support"].astype(float).fillna(0.0)
    )

    preferred_cols = [
        "v3_rank_id",
        "final_tier",
        "final_priority_score",
        "antileakage_probability_positive",
        "antileakage_predicted_label",
        "ml_predicted_probability_positive",
        "ml_predicted_label",
        "ml_label_y",
        "label_status",
        "ml_label_status",
        "candidate_class",
        "v3_priority_score",
        "stress_score",
        "length",
        "sequence",
        "total_count",
        "samples_detected",
        "log2FC_B73_DS_vs_WW",
        "log2FC_H21_DS_vs_WW",
        "n_predicted_target_genes",
        "n_opposite_support_B73",
        "n_opposite_support_H21",
        "total_opposite_support",
        "max_target_support_score",
        "mean_target_support_score",
        "feature_id",
    ]

    final = merged[[c for c in preferred_cols if c in merged.columns]].copy()

    tier_order = {
        "TIER_1_high_confidence_supported": 1,
        "TIER_2_model_prioritized_unlabeled": 2,
        "TIER_3_supported_moderate_model": 3,
        "TIER_4_unlabeled_moderate_model": 4,
        "TIER_5_rule_supported_full_model_only": 5,
        "TIER_6_low_priority_or_negative": 6,
    }

    final["tier_order"] = final["final_tier"].map(tier_order).fillna(99).astype(int)

    final = final.sort_values(
        by=[
            "tier_order",
            "final_priority_score",
            "antileakage_probability_positive",
            "v3_priority_score",
        ],
        ascending=[True, False, False, False],
    ).drop(columns=["tier_order"])

    final.to_csv(FINAL_TABLE, sep="\t", index=False)
    final.head(100).to_csv(FINAL_TOP100, sep="\t", index=False)

    tier_counts = final["final_tier"].value_counts().to_dict()

    anti_comp = pd.read_csv(ANTILEAK_COMP, sep="\t") if ANTILEAK_COMP.exists() else pd.DataFrame()
    full_comp = pd.read_csv(FULL_COMP, sep="\t") if FULL_COMP.exists() else pd.DataFrame()

    labels_report = read_json(LABELS_REPORT)
    de_report = read_json(DE_REPORT)
    quant_report = read_json(QUANT_REPORT)
    download_report = read_json(DOWNLOAD_REPORT)
    ref_report = read_json(REF_REPORT)
    anti_report = read_json(ANTILEAK_REPORT)
    full_ml_report = read_json(FULL_ML_REPORT)
    gae_report = read_json(GAE_REPORT)

    top20 = final.head(20)

    lines = []
    lines.append("# Report finale V3 - AI/GNN small RNA e RNA-Seq PRJNA419326")
    lines.append("")
    lines.append(f"Data generazione: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("## 1. Sintesi operativa")
    lines.append("")
    lines.append("Il progetto V3 integra dati small RNA-Seq e RNA-Seq pubblici PRJNA419326 con un grafo AI/GNN, un Graph Autoencoder e una fase di machine learning supervisionata su label computazionali/silver.")
    lines.append("")
    lines.append("Punti chiave:")
    lines.append("- small RNA sperimentali prioritizzati: 500")
    lines.append(f"- RNA-Seq FASTQ validati: {download_report.get('n_ok', 'NA')} su {download_report.get('n_expected_fastq', 'NA')}")
    lines.append(f"- Geni quantificati: {quant_report.get('matrix_stats', {}).get('n_genes', 'NA')}")
    lines.append(f"- Trascritti quantificati: {quant_report.get('matrix_stats', {}).get('n_transcripts', 'NA')}")
    lines.append(f"- Link smallRNA-target candidati: {labels_report.get('n_target_links', 'NA')}")
    lines.append(f"- Positivi candidati da integrazione biologica: {labels_report.get('label_counts', {}).get('positive_candidate_functional_drought_srna', 'NA')}")
    lines.append(f"- Campioni usati nella ML supervisionata: {anti_report.get('n_ml_samples', 'NA')}")
    lines.append("")
    lines.append("## 2. Nota metodologica")
    lines.append("")
    lines.append("Le label usate per la ML sono computazionali/silver: non sono validazioni sperimentali. La pipeline è però utile per generare una lista prioritaria di candidati da validare biologicamente.")
    lines.append("")
    lines.append("Sono stati prodotti due modelli:")
    lines.append("- modello full: usa anche feature derivate dal target-support; utile operativamente ma soggetto a leakage metodologico;")
    lines.append("- modello anti-leakage: esclude feature direttamente usate per costruire le label; è il risultato più difendibile.")
    lines.append("")
    lines.append("## 3. RNA-Seq")
    lines.append("")
    lines.append(f"- Status download RNA-Seq: {download_report.get('final_status', 'NA')}")
    lines.append(f"- FASTQ validati: {download_report.get('n_ok', 'NA')}")
    lines.append(f"- Dimensione validata GB: {download_report.get('downloaded_or_validated_size_gb', 'NA')}")
    lines.append(f"- Status quantificazione Salmon: {quant_report.get('status', 'NA')}")
    lines.append(f"- Campioni quantificati OK: {quant_report.get('n_quant_ok', 'NA')}")
    lines.append("")
    lines.append("## 4. Screening differenziale gene-level")
    lines.append("")
    lines.append(f"- Geni in matrice: {de_report.get('matrix_shape', {}).get('genes', 'NA')}")
    lines.append(f"- Campioni in matrice: {de_report.get('matrix_shape', {}).get('samples', 'NA')}")
    lines.append("")
    lines.append("Classi combinate B73/H21:")
    for k, v in de_report.get("combined_screen_counts", {}).items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("## 5. Label candidate small RNA")
    lines.append("")
    lines.append(f"- small RNA analizzati: {labels_report.get('n_srna', 'NA')}")
    lines.append(f"- Trascritti scansionati: {labels_report.get('n_transcripts_scanned', 'NA')}")
    lines.append(f"- Trascritti con hit: {labels_report.get('n_transcripts_with_hits', 'NA')}")
    lines.append(f"- Link smallRNA-target: {labels_report.get('n_target_links', 'NA')}")
    lines.append("")
    lines.append("Conteggio label:")
    for k, v in labels_report.get("label_counts", {}).items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("## 6. ML supervisionata anti-leakage")
    lines.append("")
    lines.append(f"- Miglior modello anti-leakage: {anti_report.get('best_model', 'NA')}")
    lines.append(f"- Feature anti-leakage usate: {anti_report.get('n_features_used', 'NA')}")
    lines.append(f"- Positivi y=1: {anti_report.get('n_positive', 'NA')}")
    lines.append(f"- Negativi y=0: {anti_report.get('n_negative', 'NA')}")
    lines.append("")
    lines.append("Metriche anti-leakage:")
    if not anti_comp.empty:
        for _, r in anti_comp.iterrows():
            lines.append(
                f"- {r['model']}: AP={float(r['average_precision']):.3f}, "
                f"AUC={float(r['roc_auc']):.3f}, F1={float(r['f1']):.3f}, "
                f"precision={float(r['precision']):.3f}, recall={float(r['recall']):.3f}"
            )
    lines.append("")
    lines.append("## 7. Distribuzione finale dei candidati")
    lines.append("")
    for k, v in tier_counts.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("## 8. Top 20 candidati finali")
    lines.append("")
    for _, r in top20.iterrows():
        lines.append(
            f"- {r['v3_rank_id']} | {r['final_tier']} | "
            f"p_anti={float(r['antileakage_probability_positive']):.3f} | "
            f"class={r['candidate_class']} | "
            f"priority={float(r['v3_priority_score']):.3f} | "
            f"target_score={float(r['max_target_support_score']):.3f}"
        )
    lines.append("")
    lines.append("## 9. File finali")
    lines.append("")
    lines.append(f"- Tabella finale candidati: `{FINAL_TABLE}`")
    lines.append(f"- Top 100 candidati: `{FINAL_TOP100}`")
    lines.append(f"- Report Markdown: `{FINAL_REPORT_MD}`")
    lines.append(f"- Report TXT: `{FINAL_REPORT_TXT}`")
    lines.append("")

    md_text = "\n".join(lines) + "\n"
    FINAL_REPORT_MD.write_text(md_text, encoding="utf-8")
    FINAL_REPORT_TXT.write_text(md_text.replace("#", ""), encoding="utf-8")

    index = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "final_dir": str(FINAL_DIR),
        "outputs": {
            "final_table": str(FINAL_TABLE),
            "final_top100": str(FINAL_TOP100),
            "final_report_md": str(FINAL_REPORT_MD),
            "final_report_txt": str(FINAL_REPORT_TXT),
        },
        "tier_counts": {str(k): int(v) for k, v in tier_counts.items()},
        "top20_v3_rank_ids": top20["v3_rank_id"].astype(str).tolist(),
    }

    FINAL_REPORT_JSON.write_text(json.dumps(index, indent=2, ensure_ascii=False), encoding="utf-8")

    print("REPORT FINALE V3 CREATO")
    print("=" * 80)
    print(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    print("")
    print("Distribuzione tier:")
    for k, v in tier_counts.items():
        print(f"- {k}: {v}")
    print("")
    print("Top 20 candidati finali:")
    for _, r in top20.iterrows():
        print(
            f"- {r['v3_rank_id']} | {r['final_tier']} | "
            f"p_anti={float(r['antileakage_probability_positive']):.3f} | "
            f"class={r['candidate_class']} | "
            f"priority={float(r['v3_priority_score']):.3f} | "
            f"target_score={float(r['max_target_support_score']):.3f}"
        )
    print("")
    print("Output:")
    print(FINAL_TABLE)
    print(FINAL_TOP100)
    print(FINAL_REPORT_MD)
    print(FINAL_REPORT_TXT)
    print(FINAL_REPORT_JSON)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
