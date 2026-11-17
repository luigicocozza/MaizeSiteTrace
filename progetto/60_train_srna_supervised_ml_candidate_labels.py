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
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

TRAINING_TABLE = BASE_DIR / "outputs/prjna419326/PRJNA419326_srna_supervised_training_table.tsv"

OUT_DIR = BASE_DIR / "outputs/prjna419326"
MODEL_DIR = BASE_DIR / "models/srna_supervised_ml_v1"

OUT_DIR.mkdir(parents=True, exist_ok=True)
MODEL_DIR.mkdir(parents=True, exist_ok=True)

ML_DATASET = OUT_DIR / "PRJNA419326_srna_supervised_ML_dataset_silver_labels.tsv"
MODEL_COMPARISON = OUT_DIR / "PRJNA419326_srna_supervised_ML_model_comparison.tsv"
PREDICTIONS = OUT_DIR / "PRJNA419326_srna_supervised_ML_predictions.tsv"
FEATURE_IMPORTANCE = OUT_DIR / "PRJNA419326_srna_supervised_ML_feature_importance.tsv"

BEST_MODEL_PATH = MODEL_DIR / "PRJNA419326_srna_supervised_best_model.joblib"
SUMMARY_TXT = OUT_DIR / "PRJNA419326_srna_supervised_ML_summary.txt"
REPORT_JSON = OUT_DIR / "PRJNA419326_srna_supervised_ML_report.json"


POSITIVE_LABEL = "positive_candidate_functional_drought_srna"


def build_silver_labels(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df["total_opposite_support"] = (
        df["n_opposite_support_B73"].astype(float).fillna(0)
        + df["n_opposite_support_H21"].astype(float).fillna(0)
    )

    df["ml_label_y"] = -1
    df["ml_label_status"] = "excluded_uncertain"

    # Positivi: quelli già definiti dal passo biologico integrato.
    pos_mask = df["label_status"].astype(str).eq(POSITIVE_LABEL)
    df.loc[pos_mask, "ml_label_y"] = 1
    df.loc[pos_mask, "ml_label_status"] = "positive_silver_from_integrated_evidence"

    # Negativi silver:
    # tra gli incerti, scegliamo quelli con nessun supporto opposto,
    # basso supporto target e bassa priorità relativa.
    uncertain = df["label_status"].astype(str).eq("uncertain_unlabeled")

    u = df.loc[uncertain].copy()

    if len(u) == 0:
        return df

    priority_q40 = float(u["v3_priority_score"].quantile(0.40))
    target_q40 = float(u["max_target_support_score"].quantile(0.40))

    neg_mask = (
        uncertain
        & (df["total_opposite_support"] == 0)
        & (df["v3_priority_score"].astype(float) <= priority_q40)
        & (df["max_target_support_score"].astype(float) <= target_q40)
    )

    # Se i negativi sono pochi, rilassiamo la soglia.
    if int(neg_mask.sum()) < 40:
        priority_q60 = float(u["v3_priority_score"].quantile(0.60))
        target_q60 = float(u["max_target_support_score"].quantile(0.60))

        neg_mask = (
            uncertain
            & (df["total_opposite_support"] == 0)
            & (df["v3_priority_score"].astype(float) <= priority_q60)
            & (df["max_target_support_score"].astype(float) <= target_q60)
        )

    # Se ancora pochi, prendiamo i peggiori incerti per supporto complessivo.
    if int(neg_mask.sum()) < 40:
        u = u.copy()
        u["negative_score"] = (
            u["v3_priority_score"].rank(pct=True)
            + u["max_target_support_score"].rank(pct=True)
            + u["total_opposite_support"].rank(pct=True)
        )
        n_take = min(100, len(u))
        negative_ids = set(u.sort_values("negative_score").head(n_take)["v3_rank_id"].astype(str))
        neg_mask = uncertain & df["v3_rank_id"].astype(str).isin(negative_ids)

    df.loc[neg_mask, "ml_label_y"] = 0
    df.loc[neg_mask, "ml_label_status"] = "negative_silver_low_support_uncertain"

    return df


def get_feature_columns(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    numeric_candidates = [
        "length",
        "v3_priority_score",
        "stress_score",
        "total_count",
        "samples_detected",
        "log2FC_B73_DS_vs_WW",
        "log2FC_H21_DS_vs_WW",
        "n_predicted_target_genes",
        "n_opposite_support_B73",
        "n_opposite_support_H21",
        "max_target_support_score",
        "mean_target_support_score",
        "total_opposite_support",
    ]

    z_cols = [c for c in df.columns if c.startswith("z_")]

    numeric_cols = [c for c in numeric_candidates + z_cols if c in df.columns]
    categorical_cols = [c for c in ["candidate_class"] if c in df.columns]

    return numeric_cols, categorical_cols


def make_preprocessor(numeric_cols: list[str], categorical_cols: list[str]) -> ColumnTransformer:
    numeric_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )

    categorical_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", OneHotEncoder(handle_unknown="ignore")),
        ]
    )

    return ColumnTransformer(
        transformers=[
            ("num", numeric_pipe, numeric_cols),
            ("cat", categorical_pipe, categorical_cols),
        ],
        remainder="drop",
    )


def evaluate_model(name: str, model: Pipeline, X: pd.DataFrame, y: np.ndarray, cv: StratifiedKFold) -> dict:
    pred = cross_val_predict(model, X, y, cv=cv, method="predict")

    try:
        proba = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
    except Exception:
        proba = pred.astype(float)

    out = {
        "model": name,
        "accuracy": accuracy_score(y, pred),
        "precision": precision_score(y, pred, zero_division=0),
        "recall": recall_score(y, pred, zero_division=0),
        "f1": f1_score(y, pred, zero_division=0),
        "roc_auc": roc_auc_score(y, proba) if len(np.unique(y)) == 2 else np.nan,
        "average_precision": average_precision_score(y, proba) if len(np.unique(y)) == 2 else np.nan,
        "n_samples": int(len(y)),
        "n_positive": int((y == 1).sum()),
        "n_negative": int((y == 0).sum()),
    }

    return out


def main() -> int:
    df = pd.read_csv(TRAINING_TABLE, sep="\t")

    df = build_silver_labels(df)

    ml_df = df[df["ml_label_y"].isin([0, 1])].copy()
    ml_df["ml_label_y"] = ml_df["ml_label_y"].astype(int)

    label_counts = ml_df["ml_label_y"].value_counts().to_dict()
    status_counts = df["ml_label_status"].value_counts().to_dict()

    if len(label_counts) < 2:
        raise SystemExit("ERRORE: non ci sono entrambe le classi 0 e 1 per il training supervisionato.")

    numeric_cols, categorical_cols = get_feature_columns(ml_df)

    feature_cols = numeric_cols + categorical_cols

    X = ml_df[feature_cols].copy()
    y = ml_df["ml_label_y"].to_numpy()

    min_class = int(min((y == 0).sum(), (y == 1).sum()))
    n_splits = max(2, min(5, min_class))

    cv = StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=42,
    )

    models = {
        "logistic_regression": LogisticRegression(
            max_iter=2000,
            class_weight="balanced",
            solver="liblinear",
            random_state=42,
        ),
        "random_forest": RandomForestClassifier(
            n_estimators=500,
            max_depth=None,
            min_samples_leaf=3,
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
        "gradient_boosting": GradientBoostingClassifier(
            random_state=42,
        ),
    }

    rows = []
    fitted_models = {}

    for name, estimator in models.items():
        pipe = Pipeline(
            steps=[
                ("preprocess", make_preprocessor(numeric_cols, categorical_cols)),
                ("model", estimator),
            ]
        )

        metrics = evaluate_model(name, pipe, X, y, cv)
        rows.append(metrics)

        pipe.fit(X, y)
        fitted_models[name] = pipe

    comparison = pd.DataFrame(rows).sort_values(
        by=["average_precision", "roc_auc", "f1"],
        ascending=False,
    )

    comparison.to_csv(MODEL_COMPARISON, sep="\t", index=False)

    best_name = str(comparison.iloc[0]["model"])
    best_model = fitted_models[best_name]

    joblib.dump(best_model, BEST_MODEL_PATH)

    # Fit finale e predizioni su tutti i 500 small RNA
    all_X = df[feature_cols].copy()
    pred_proba = best_model.predict_proba(all_X)[:, 1]
    pred_label = best_model.predict(all_X)

    pred_df = df.copy()
    pred_df["ml_model_best"] = best_name
    pred_df["ml_predicted_probability_positive"] = pred_proba
    pred_df["ml_predicted_label"] = pred_label

    pred_df = pred_df.sort_values(
        by=["ml_predicted_probability_positive", "v3_priority_score", "max_target_support_score"],
        ascending=False,
    )

    pred_df.to_csv(PREDICTIONS, sep="\t", index=False)

    # Feature importance per Random Forest o Gradient Boosting
    importance_written = False
    try:
        model_step = best_model.named_steps["model"]
        preprocess = best_model.named_steps["preprocess"]

        feature_names = preprocess.get_feature_names_out()

        if hasattr(model_step, "feature_importances_"):
            imp = pd.DataFrame(
                {
                    "feature": feature_names,
                    "importance": model_step.feature_importances_,
                }
            ).sort_values("importance", ascending=False)

            imp.to_csv(FEATURE_IMPORTANCE, sep="\t", index=False)
            importance_written = True
    except Exception:
        importance_written = False

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "methodological_note": (
            "Supervised ML uses computational silver labels. Positives come from integrated smallRNA/RNA-Seq evidence; "
            "negatives are weak/silver negatives selected from uncertain low-support candidates. "
            "This is not experimental validation."
        ),
        "inputs": {
            "training_table": str(TRAINING_TABLE),
        },
        "label_counts_all_status": status_counts,
        "ml_label_counts": {str(k): int(v) for k, v in label_counts.items()},
        "n_ml_samples": int(len(ml_df)),
        "n_total_candidates_scored": int(len(df)),
        "numeric_features": numeric_cols,
        "categorical_features": categorical_cols,
        "n_splits_cv": n_splits,
        "best_model": best_name,
        "outputs": {
            "ml_dataset": str(ML_DATASET),
            "model_comparison": str(MODEL_COMPARISON),
            "predictions": str(PREDICTIONS),
            "feature_importance": str(FEATURE_IMPORTANCE) if importance_written else None,
            "best_model_path": str(BEST_MODEL_PATH),
            "summary_txt": str(SUMMARY_TXT),
            "report_json": str(REPORT_JSON),
        },
    }

    df.to_csv(ML_DATASET, sep="\t", index=False)

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    top_predictions = pred_df.head(20)

    lines = []
    lines.append("ML SUPERVISIONATA small RNA - LABEL CANDIDATE / SILVER LABELS")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("Nota metodologica:")
    lines.append("- Training supervisionato su label computazionali/silver, non su validazione sperimentale.")
    lines.append("- Positivi: evidenza integrata small RNA + RNA-Seq.")
    lines.append("- Negativi: candidati incerti con basso supporto biologico relativo.")
    lines.append("")
    lines.append("Conteggio status label complessivo:")
    for k, v in status_counts.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("Dataset ML:")
    lines.append(f"- campioni ML: {len(ml_df)}")
    lines.append(f"- positivi y=1: {int((y == 1).sum())}")
    lines.append(f"- negativi y=0: {int((y == 0).sum())}")
    lines.append(f"- candidati totali predetti: {len(df)}")
    lines.append(f"- feature numeriche: {len(numeric_cols)}")
    lines.append(f"- feature categoriche: {len(categorical_cols)}")
    lines.append(f"- CV folds: {n_splits}")
    lines.append("")
    lines.append("Confronto modelli:")
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
    lines.append(f"Miglior modello: {best_name}")
    lines.append("")
    lines.append("Top 20 candidati predetti dal modello:")
    for _, r in top_predictions.iterrows():
        lines.append(
            f"- {r['v3_rank_id']} | class={r['candidate_class']} | "
            f"p_positive={float(r['ml_predicted_probability_positive']):.3f} | "
            f"ml_label={int(r['ml_label_y'])} | "
            f"original_status={r['label_status']} | "
            f"priority={float(r['v3_priority_score']):.3f} | "
            f"target_score={float(r['max_target_support_score']):.3f}"
        )
    lines.append("")
    lines.append("Output:")
    lines.append(str(ML_DATASET))
    lines.append(str(MODEL_COMPARISON))
    lines.append(str(PREDICTIONS))
    if importance_written:
        lines.append(str(FEATURE_IMPORTANCE))
    lines.append(str(BEST_MODEL_PATH))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
