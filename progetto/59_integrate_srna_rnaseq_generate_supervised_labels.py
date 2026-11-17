#!/usr/bin/env python3
from __future__ import annotations

import gzip
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()

SRNA_RANKING = BASE_DIR / "outputs/prjna419326/PRJNA419326_experimental_srna_v3_ranking.csv"
RNASEQ_DE = BASE_DIR / "outputs/prjna419326/PRJNA419326_rnaseq_gene_DE_combined_B73_H21.tsv"

REF_DIR = BASE_DIR / "data/reference/zea_mays_ensembl_plants"
TRANSCRIPTOME = REF_DIR / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna_plus_ncrna.fa.gz"
TX2GENE = REF_DIR / "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.tx2gene.tsv"

EMB_EXPERIMENTAL = BASE_DIR / "data/processed/embeddings/gae_v3_experimental_srna_embeddings.csv"

OUT_DIR = BASE_DIR / "outputs/prjna419326"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_LINKS = OUT_DIR / "PRJNA419326_srna_rnaseq_predicted_target_links.tsv"
LABELS_TSV = OUT_DIR / "PRJNA419326_srna_supervised_label_candidates.tsv"
TRAINING_TABLE = OUT_DIR / "PRJNA419326_srna_supervised_training_table.tsv"
SUMMARY_TXT = OUT_DIR / "PRJNA419326_srna_supervised_labels_summary.txt"
REPORT_JSON = OUT_DIR / "PRJNA419326_srna_supervised_labels_report.json"


K = 12
MAX_TARGETS_PER_SRNA = 100
MIN_TARGET_SUPPORT_SCORE_FOR_POSITIVE = 2.0
MIN_PRIORITY_SCORE_FOR_POSITIVE = 0.50
MAX_PRIORITY_SCORE_FOR_NEGATIVE = 0.20


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")
    return str(seq).translate(table)[::-1].upper().replace("U", "T")


def norm_seq(seq: str) -> str:
    return str(seq).upper().replace("U", "T").replace(" ", "").replace("\n", "")


def fasta_iter_gz(path: Path):
    name = None
    seq_parts = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_parts).upper().replace("U", "T")
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    if name is not None:
        yield name, "".join(seq_parts).upper().replace("U", "T")


def srna_direction_by_genotype(candidate_class: str, genotype: str) -> str:
    c = str(candidate_class)

    if c == "conserved_up":
        return "up"
    if c == "conserved_down":
        return "down"

    if genotype == "B73":
        if c == "B73_up_only":
            return "up"
        if c == "B73_down_only":
            return "down"
        if c == "opposite_response":
            return "mixed"

    if genotype == "H21":
        if c == "H21_up_only":
            return "up"
        if c == "H21_down_only":
            return "down"
        if c == "opposite_response":
            return "mixed"

    return "stable_or_unknown"


def opposite_support(srna_dir: str, gene_log2fc: float, gene_pvalue: float) -> tuple[bool, float, str]:
    """
    Coerenza biologica semplificata:
    sRNA up  -> target mRNA down
    sRNA down -> target mRNA up
    """
    try:
        log2fc = float(gene_log2fc)
    except Exception:
        log2fc = 0.0

    try:
        p = float(gene_pvalue)
    except Exception:
        p = 1.0

    if srna_dir == "up" and log2fc <= -1.0 and p <= 0.10:
        return True, abs(log2fc) * (-np.log10(max(p, 1e-300))), "srna_up_gene_down"

    if srna_dir == "down" and log2fc >= 1.0 and p <= 0.10:
        return True, abs(log2fc) * (-np.log10(max(p, 1e-300))), "srna_down_gene_up"

    return False, 0.0, "not_opposite_or_not_significant"


def build_query_kmers(srna: pd.DataFrame) -> dict[str, list[dict]]:
    """
    Crea kmers reverse-complement da usare come euristica target.
    Non è una predizione validata tipo psRNATarget: è una prima inferenza locale.
    """
    q = defaultdict(list)

    for _, row in srna.iterrows():
        srna_id = str(row["v3_rank_id"])
        seq = norm_seq(row["sequence"])
        rc = revcomp(seq)

        kmers = set()

        if len(rc) >= K:
            kmers.add(rc[:K])
            kmers.add(rc[-K:])
            mid = max((len(rc) - K) // 2, 0)
            kmers.add(rc[mid:mid + K])

        for kmer in kmers:
            if len(kmer) == K and "N" not in kmer:
                q[kmer].append(
                    {
                        "srna_id": srna_id,
                        "srna_sequence": seq,
                        "srna_rc": rc,
                        "candidate_class": row["candidate_class"],
                        "v3_priority_score": float(row["v3_priority_score"]),
                    }
                )

    return q


def main() -> int:
    srna = pd.read_csv(SRNA_RANKING)
    de = pd.read_csv(RNASEQ_DE, sep="\t")
    tx2gene = pd.read_csv(TX2GENE, sep="\t")

    emb = None
    if EMB_EXPERIMENTAL.exists():
        emb = pd.read_csv(EMB_EXPERIMENTAL)

    tx_to_gene = dict(zip(tx2gene["transcript_id"].astype(str), tx2gene["gene_id"].astype(str)))

    de_by_gene = de.set_index("gene_id").to_dict(orient="index")

    query_kmers = build_query_kmers(srna)

    srna_gene_hits = defaultdict(lambda: {
        "matched_kmers": set(),
        "transcripts": set(),
        "full_rc_match": False,
    })

    n_transcripts_scanned = 0
    n_transcripts_with_hits = 0

    for transcript_id, seq in fasta_iter_gz(TRANSCRIPTOME):
        n_transcripts_scanned += 1
        gene_id = tx_to_gene.get(transcript_id, transcript_id)

        seen_this_transcript = False

        # kmer scan mirato
        if len(seq) >= K:
            for i in range(0, len(seq) - K + 1):
                kmer = seq[i:i + K]
                if kmer in query_kmers:
                    for item in query_kmers[kmer]:
                        key = (item["srna_id"], gene_id)
                        srna_gene_hits[key]["matched_kmers"].add(kmer)
                        srna_gene_hits[key]["transcripts"].add(transcript_id)
                        seen_this_transcript = True

        # full reverse-complement exact match solo per gli sRNA con almeno un kmer già compatibile
        if seen_this_transcript:
            n_transcripts_with_hits += 1

    link_rows = []

    srna_meta = srna.set_index("v3_rank_id").to_dict(orient="index")

    for (srna_id, gene_id), hit in srna_gene_hits.items():
        if gene_id not in de_by_gene:
            continue

        meta = srna_meta.get(srna_id, {})
        candidate_class = str(meta.get("candidate_class", ""))
        srna_seq = norm_seq(meta.get("sequence", ""))
        srna_rc = revcomp(srna_seq)

        matched_kmers = sorted(hit["matched_kmers"])
        transcripts = sorted(hit["transcripts"])

        # full match controllato sui trascritti candidati
        full_match = False
        # Per evitare rilettura pesante del FASTA qui non facciamo full scan.
        # La label usa soprattutto numero di kmer e coerenza espressione.

        b73_dir = srna_direction_by_genotype(candidate_class, "B73")
        h21_dir = srna_direction_by_genotype(candidate_class, "H21")

        de_info = de_by_gene[gene_id]

        b73_ok, b73_score, b73_note = opposite_support(
            b73_dir,
            de_info.get("B73_log2FC_TPM_DS_vs_WW", 0.0),
            de_info.get("B73_pvalue", 1.0),
        )

        h21_ok, h21_score, h21_note = opposite_support(
            h21_dir,
            de_info.get("H21_log2FC_TPM_DS_vs_WW", 0.0),
            de_info.get("H21_pvalue", 1.0),
        )

        kmer_support = len(matched_kmers)

        expression_support_score = float(b73_score + h21_score)
        target_support_score = float(kmer_support + expression_support_score)

        link_rows.append(
            {
                "v3_rank_id": srna_id,
                "gene_id": gene_id,
                "candidate_class": candidate_class,
                "v3_priority_score": float(meta.get("v3_priority_score", 0.0)),
                "srna_sequence": srna_seq,
                "n_matched_12mers": kmer_support,
                "matched_12mers": ";".join(matched_kmers),
                "n_transcripts_hit": len(transcripts),
                "transcripts_hit": ";".join(transcripts[:20]),
                "B73_srna_direction": b73_dir,
                "H21_srna_direction": h21_dir,
                "B73_gene_log2FC": de_info.get("B73_log2FC_TPM_DS_vs_WW", 0.0),
                "B73_gene_pvalue": de_info.get("B73_pvalue", 1.0),
                "H21_gene_log2FC": de_info.get("H21_log2FC_TPM_DS_vs_WW", 0.0),
                "H21_gene_pvalue": de_info.get("H21_pvalue", 1.0),
                "B73_opposite_support": b73_ok,
                "H21_opposite_support": h21_ok,
                "B73_support_note": b73_note,
                "H21_support_note": h21_note,
                "expression_support_score": expression_support_score,
                "target_support_score": target_support_score,
                "gene_combined_class": de_info.get("combined_screen_class", ""),
                "gene_combined_de_priority_score": de_info.get("combined_de_priority_score", 0.0),
            }
        )

    links = pd.DataFrame(link_rows)

    if not links.empty:
        links = links.sort_values(
            by=["target_support_score", "expression_support_score", "n_matched_12mers"],
            ascending=False,
        )

        # limitiamo target per sRNA per evitare tabella enorme
        links = links.groupby("v3_rank_id", group_keys=False).head(MAX_TARGETS_PER_SRNA)

    links.to_csv(TARGET_LINKS, sep="\t", index=False)

    # Summary per sRNA
    label_rows = []

    for _, row in srna.iterrows():
        srna_id = str(row["v3_rank_id"])
        subset = links[links["v3_rank_id"] == srna_id] if not links.empty else pd.DataFrame()

        n_targets = int(len(subset))
        n_opposite_b73 = int(subset["B73_opposite_support"].sum()) if n_targets else 0
        n_opposite_h21 = int(subset["H21_opposite_support"].sum()) if n_targets else 0
        max_target_support = float(subset["target_support_score"].max()) if n_targets else 0.0
        mean_target_support = float(subset["target_support_score"].mean()) if n_targets else 0.0

        priority = float(row["v3_priority_score"])

        positive_rule = (
            priority >= MIN_PRIORITY_SCORE_FOR_POSITIVE
            and max_target_support >= MIN_TARGET_SUPPORT_SCORE_FOR_POSITIVE
            and (n_opposite_b73 + n_opposite_h21) >= 1
        )

        negative_rule = (
            priority <= MAX_PRIORITY_SCORE_FOR_NEGATIVE
            and n_targets == 0
        )

        if positive_rule:
            y = 1
            label_status = "positive_candidate_functional_drought_srna"
        elif negative_rule:
            y = 0
            label_status = "negative_candidate_low_priority_no_target_support"
        else:
            y = -1
            label_status = "uncertain_unlabeled"

        label_rows.append(
            {
                "v3_rank_id": srna_id,
                "feature_id": row["feature_id"],
                "sequence": row["sequence"],
                "length": row["length"],
                "candidate_class": row["candidate_class"],
                "v3_priority_score": priority,
                "stress_score": row["stress_score"],
                "total_count": row["total_count"],
                "samples_detected": row["samples_detected"],
                "log2FC_B73_DS_vs_WW": row["log2FC_B73_DS_vs_WW"],
                "log2FC_H21_DS_vs_WW": row["log2FC_H21_DS_vs_WW"],
                "n_predicted_target_genes": n_targets,
                "n_opposite_support_B73": n_opposite_b73,
                "n_opposite_support_H21": n_opposite_h21,
                "max_target_support_score": max_target_support,
                "mean_target_support_score": mean_target_support,
                "supervised_label_y": y,
                "label_status": label_status,
            }
        )

    labels = pd.DataFrame(label_rows)

    labels = labels.sort_values(
        by=[
            "supervised_label_y",
            "max_target_support_score",
            "v3_priority_score",
            "stress_score",
        ],
        ascending=[False, False, False, False],
    )

    labels.to_csv(LABELS_TSV, sep="\t", index=False)

    training = labels.copy()

    if emb is not None:
        emb = emb.copy()
        emb["v3_rank_id"] = emb["node_id"].astype(str).str.replace("experimental_srna:", "", regex=False)
        z_cols = [c for c in emb.columns if c.startswith("z_")]
        training = training.merge(
            emb[["v3_rank_id"] + z_cols],
            on="v3_rank_id",
            how="left",
        )

    training.to_csv(TRAINING_TABLE, sep="\t", index=False)

    label_counts = labels["label_status"].value_counts().to_dict()
    y_counts = labels["supervised_label_y"].value_counts().to_dict()

    top_positive = labels[labels["supervised_label_y"] == 1].head(20)

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "methodological_note": (
            "Labels are computational candidate labels, not experimentally validated labels. "
            "Target links use local 12-mer reverse-complement evidence plus RNA-Seq opposite-direction support."
        ),
        "inputs": {
            "srna_ranking": str(SRNA_RANKING),
            "rnaseq_de": str(RNASEQ_DE),
            "transcriptome": str(TRANSCRIPTOME),
            "tx2gene": str(TX2GENE),
            "gae_embeddings": str(EMB_EXPERIMENTAL) if EMB_EXPERIMENTAL.exists() else None,
        },
        "parameters": {
            "kmer_size": K,
            "max_targets_per_srna": MAX_TARGETS_PER_SRNA,
            "min_target_support_score_for_positive": MIN_TARGET_SUPPORT_SCORE_FOR_POSITIVE,
            "min_priority_score_for_positive": MIN_PRIORITY_SCORE_FOR_POSITIVE,
            "max_priority_score_for_negative": MAX_PRIORITY_SCORE_FOR_NEGATIVE,
        },
        "n_srna": int(len(srna)),
        "n_transcripts_scanned": int(n_transcripts_scanned),
        "n_transcripts_with_hits": int(n_transcripts_with_hits),
        "n_target_links": int(len(links)),
        "label_counts": label_counts,
        "y_counts": {str(k): int(v) for k, v in y_counts.items()},
        "outputs": {
            "target_links": str(TARGET_LINKS),
            "labels_tsv": str(LABELS_TSV),
            "training_table": str(TRAINING_TABLE),
            "summary_txt": str(SUMMARY_TXT),
            "report_json": str(REPORT_JSON),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("INTEGRAZIONE small RNA + RNA-SEQ PER LABEL SUPERVISIONATE")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")
    lines.append("Nota metodologica:")
    lines.append("- Le label sono computazionali/candidate, non validazione sperimentale.")
    lines.append("- Target evidence: 12-mer reverse-complement + supporto RNA-Seq in direzione opposta.")
    lines.append("")
    lines.append(f"small RNA analizzati: {len(srna)}")
    lines.append(f"Trascritti scansionati: {n_transcripts_scanned}")
    lines.append(f"Trascritti con hit: {n_transcripts_with_hits}")
    lines.append(f"Link smallRNA-target prodotti: {len(links)}")
    lines.append("")
    lines.append("Conteggio label:")
    for k, v in label_counts.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("Conteggio y:")
    for k, v in sorted(y_counts.items(), key=lambda x: str(x[0])):
        lines.append(f"- y={k}: {v}")
    lines.append("")
    lines.append("Top positive candidates:")
    if len(top_positive) == 0:
        lines.append("- Nessun candidato positivo con le soglie attuali.")
    else:
        for _, r in top_positive.iterrows():
            lines.append(
                f"- {r['v3_rank_id']} | class={r['candidate_class']} | "
                f"priority={float(r['v3_priority_score']):.3f} | "
                f"targets={int(r['n_predicted_target_genes'])} | "
                f"opp_B73={int(r['n_opposite_support_B73'])} | "
                f"opp_H21={int(r['n_opposite_support_H21'])} | "
                f"max_target_score={float(r['max_target_support_score']):.3f}"
            )
    lines.append("")
    lines.append("Output:")
    lines.append(str(TARGET_LINKS))
    lines.append(str(LABELS_TSV))
    lines.append(str(TRAINING_TABLE))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
