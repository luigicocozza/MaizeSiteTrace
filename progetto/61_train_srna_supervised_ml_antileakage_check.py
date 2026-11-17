#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


BASE_DIR = Path(".").resolve()

INPUT_DATASET = BASE_DIR / "outputs/prjna419326/PRJNA419326_srna_supervised_ML_dataset_silver_labels.tsv"

OUT_DIR = BASE_DIR / "outputs/prjna419326"
MODEL_DIR = BASE_DIR / "models/srna_supervised_ml_v1_antileakage"

OUT_DIR.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)

COMPARISON_TSV = OUT_DIR / "PRJNA419326_srna_supervised_ML_antileakage_model_comparison.tsv"
PREDICTIONS_TSV = OUT_DIR / "PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv"
FEATURES_TSV = OUT_DIR / "PRJNA419326_srna_supervised_ML_antileakage_features_used.tsv"
BEST_MODEL = MODEL_DIR / "PRJNA419326_srna_supervised_antileakage_best_model.joblib"

SUMMARY_TXT = OUT_DIR / "PRJNA419326_srna_supervised_ML_antileakage_summary.txt"
REPORT_JSON = OUT_DIR / "PRJNA419326_srna_supervised_ML_antileakage_report.json"


# Escludiamo le feature usate direttamente o quasi direttamente per generare le silver label.
LEAKAGE_COLUMNS = {
    "v3_priority_score",
    "n_predicted_target_genes",
    "n_opposite_support_B73",
    "n_opposite_support_H21",
    "max_target_support_score",
    "mean_target_support_score",
    "total_opposite_support",
    "candidate_class",
    "label_status",
    "ml_label_status",
    "supervised_label_y",
}


def make_preprocessor(numeric_cols: list[str]) -> ColumnTransformer:
    numeric_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )

    return ColumnTransformer(
        transformers=[
            ("num", numeric_pipe, numeric_cols),
        ],
        remainder="drop",
    )


def evaluate_model(name: str, model: Pipeline, X: pd.DataFrame, y: np.ndarray, cv: StratifiedKFold) -> dict:
    pred = cross_val_predict(model, X, y, cv=cv, method="predict")

    try:
        proba = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
    except Exception:
        proba = pred.astype(float)

    return {
        "model": name,
        "accuracy": accuracy_score(y, pred),
        "precision": precision_score(y, pred, zero_division=0),
        "recall": recall_score(y, pred, zero_division=0),
        "f1": f1_score(y, pred, zero_division=0),
        "roc_auc": roc_auc_score(y, proba),
        "average_precision": average_precision_score(y, proba),
        "n_samples": int(len(y)),
        "n_positive": int((y == 1).sum()),
        "n_negative": int((y == 0).sum()),
    }


def main() -> int:
    df = pd.read_csv(INPUT_DATASET, sep="\t")

    ml_df = df[df["ml_label_y"].isin([0, 1])].copy()
    ml_df["ml_label_y"] = ml_df["ml_label_y"].astype(int)

    # Feature ammesse:
    # - proprietà semplici dello small RNA
    # - log2FC small RNA
    # - embedding GAE z_*
    # Escluse feature target/supporto usate per generare direttamente la label.
    candidate_numeric = [
        "length",
        "stress_score",
        "total_count",
        "samples_detected",
        "log2FC_B73_DS_vs_WW",
        "log2FC_H21_DS_vs_WW",
    ]

    z_cols = [c for c in ml_df.columns if c.startswith("z_")]

    numeric_cols = [
        c for c in candidate_numeric + z_cols
        if c in ml_df.columns and c not in LEAKAGE_COLUMNS
    ]

    if not numeric_cols:
        raise SystemExit("ERRORE: nessuna feature anti-leakage disponibile.")

    X = ml_df[numeric_cols].copy()
    y = ml_df["ml_label_y"].to_numpy()

    min_class = int(min((y == 0).sum(), (y == 1).sum()))
    n_splits = max(2, min(5, min_class))

    cv = StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=42,
    )

    models = {
        "logistic_regression_antileakage": LogisticRegression(
            max_iter=2000,
            class_weight="balanced",
            solver="liblinear",
            random_state=42,
        ),
        "random_forest_antileakage": RandomForestClassifier(
            n_estimators=500,
            min_samples_leaf=5,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
        "gradient_boosting_antileakage": GradientBoostingClassifier(
            random_state=42,
        ),
    }

    rows = []
    fitted = {}

    for name, estimator in models.items():
        pipe = Pipeline(
            steps=[
                ("preprocess", make_preprocessor(numeric_cols)),
                ("model", estimator),
            ]
        )

        metrics = evaluate_model(name, pipe, X, y, cv)
        rows.append(metrics)

        pipe.fit(X, y)
        fitted[name] = pipe

    comparison = pd.DataFrame(rows).sort_values(
        by=["average_precision", "roc_auc", "f1"],
        ascending=False,
    )

    best_name = str(comparison.iloc[0]["model"])
    best_model = fitted[best_name]

    joblib.dump(best_model, BEST_MODEL)

    comparison.to_csv(COMPARISON_TSV, sep="\t", index=False)

    all_df = df.copy()
    all_X = all_df[numeric_cols].copy()

    proba = best_model.predict_proba(all_X)[:, 1]
    pred = best_model.predict(all_X)

    all_df["antileakage_model_best"] = best_name
    all_df["antileakage_probability_positive"] = proba
    all_df["antileakage_predicted_label"] = pred

    all_df = all_df.sort_values(
        by=["antileakage_probability_positive", "stress_score", "total_count"],
        ascending=False,
    )

    all_df.to_csv(PREDICTIONS_TSV, sep="\t", index=False)

    features = pd.DataFrame(
        {
            "feature": numeric_cols,
            "feature_group": [
                "gae_embedding" if c.startswith("z_") else "srna_observed_feature"
                for c in numeric_cols
            ],
        }
    )
    features.to_csv(FEATURES_TSV, sep="\t", index=False)

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "methodological_note": (
            "Anti-leakage model excludes target-support and priority columns directly involved "
            "in creating the silver labels. Metrics are therefore a more conservative diagnostic, "
            "but labels remain computational/silver."
        ),
        "input_dataset": str(INPUT_DATASET),
        "n_ml_samples": int(len(ml_df)),
        "n_positive": int((y == 1).sum()),
        "n_negative": int((y == 0).sum()),
        "n_total_scored": int(len(df)),
        "n_features_used": int(len(numeric_cols)),
        "features_used": numeric_cols,
        "excluded_leakage_columns": sorted(LEAKAGE_COLUMNS),
        "n_splits_cv": int(n_splits),
        "best_model": best_name,
        "outputs": {
            "comparison_tsv": str(COMPARISON_TSV),
            "predictions_tsv": str(PREDICTIONS_TSV),
            "features_tsv": str(FEATURES_TSV),
            "best_model": str(BEST_MODEL),
            "summary_txt": str(SUMMARY_TXT),
            "report_json": str(REPORT_JSON),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("CONTROLLO ML ANTI-LEAKAGE small RNA")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("Nota metodologica:")
    lines.append("- Escluse le feature direttamente usate per generare le silver label.")
    lines.append("- Il risultato resta basato su label computazionali, ma è più difendibile del modello full.")
    lines.append("")
    lines.append("Dataset:")
    lines.append(f"- campioni ML: {len(ml_df)}")
    lines.append(f"- positivi y=1: {int((y == 1).sum())}")
    lines.append(f"- negativi y=0: {int((y == 0).sum())}")
    lines.append(f"- candidati totali predetti: {len(df)}")
    lines.append(f"- feature anti-leakage usate: {len(numeric_cols)}")
    lines.append(f"- CV folds: {n_splits}")
    lines.append("")
    lines.append("Confronto modelli anti-leakage:")
    for _, r in comparison.iterrows():
        lines.append(
            f"- {r['model']} | "
            f"AP={float(r['average_precision']):.3f} | "
            f"AUC={float(r['roc_auc']):.3f} | "
            f"F1={float(r['f1']):.3f} | "
            f"precision={float(r['precision']):.3f} | "
            f"recall={float(r['recall']):.3f}"
        )
    lines.append("")
    lines.append(f"Miglior modello anti-leakage: {best_name}")
    lines.append("")
    lines.append("Feature usate:")
    for c in numeric_cols:
        lines.append(f"- {c}")
    lines.append("")
    lines.append("Top 20 candidati predetti anti-leakage:")
    for _, r in all_df.head(20).iterrows():
        lines.append(
            f"- {r['v3_rank_id']} | class={r['candidate_class']} | "
            f"p_positive={float(r['antileakage_probability_positive']):.3f} | "
            f"silver_y={int(r['ml_label_y'])} | "
            f"priority={float(r['v3_priority_score']):.3f} | "
            f"target_score={float(r['max_target_support_score']):.3f}"
        )
    lines.append("")
    lines.append("Output:")
    lines.append(str(COMPARISON_TSV))
    lines.append(str(PREDICTIONS_TSV))
    lines.append(str(FEATURES_TSV))
    lines.append(str(BEST_MODEL))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
