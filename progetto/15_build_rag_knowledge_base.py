#!/usr/bin/env python3
"""
15_build_rag_knowledge_base.py

Costruisce una knowledge base testuale/strutturata per la chat LLM/RAG
del progetto maize_drought_ncrna_gnn.

Output:
- data/processed/rag/project_knowledge_base.jsonl
- outputs/rag/rag_knowledge_base_report.json
- outputs/rag/rag_knowledge_base_summary.txt
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Costruisce knowledge base RAG per chat LLM del progetto."
    )
    parser.add_argument("--output-jsonl", default="data/processed/rag/project_knowledge_base.jsonl")
    parser.add_argument("--report-json", default="outputs/rag/rag_knowledge_base_report.json")
    parser.add_argument("--summary-txt", default="outputs/rag/rag_knowledge_base_summary.txt")
    parser.add_argument("--chunk-size", type=int, default=1200)
    parser.add_argument("--chunk-overlap", type=int, default=150)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


def read_text(path: Path) -> str:
    if not path.exists():
        logging.warning("File testo mancante: %s", path)
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def read_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        logging.warning("File JSON mancante: %s", path)
        return {}
    try:
        with path.open("r", encoding="utf-8") as f:
            obj = json.load(f)
        return obj if isinstance(obj, dict) else {}
    except Exception as exc:
        logging.warning("Errore JSON %s: %s", path, exc)
        return {}


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        logging.warning("File CSV mancante: %s", path)
        return []
    try:
        with path.open("r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))
    except Exception as exc:
        logging.warning("Errore CSV %s: %s", path, exc)
        return []


def clean_text(text: str) -> str:
    text = text.replace("\r\n", "\n")
    text = re.sub(r"\n{3,}", "\n\n", text)
    text = re.sub(r"[ \t]+", " ", text)
    return text.strip()


def chunk_text(text: str, chunk_size: int, overlap: int) -> List[str]:
    text = clean_text(text)
    if not text:
        return []

    chunks = []
    start = 0

    while start < len(text):
        end = min(start + chunk_size, len(text))
        chunk = text[start:end].strip()

        if chunk:
            chunks.append(chunk)

        if end >= len(text):
            break

        start = max(0, end - overlap)

    return chunks


def make_record(
    chunk_id: str,
    source: str,
    source_type: str,
    title: str,
    text: str,
    metadata: Dict[str, Any] | None = None,
) -> Dict[str, Any]:
    if metadata is None:
        metadata = {}

    return {
        "chunk_id": chunk_id,
        "source": source,
        "source_type": source_type,
        "title": title,
        "text": clean_text(text),
        "metadata": metadata,
    }


def add_text_file_chunks(
    records: List[Dict[str, Any]],
    path: Path,
    source_type: str,
    title: str,
    chunk_size: int,
    chunk_overlap: int,
) -> None:
    text = read_text(path)
    chunks = chunk_text(text, chunk_size=chunk_size, overlap=chunk_overlap)

    for i, chunk in enumerate(chunks, start=1):
        records.append(
            make_record(
                chunk_id=f"{source_type}_{i:04d}",
                source=str(path),
                source_type=source_type,
                title=f"{title} - chunk {i}",
                text=chunk,
                metadata={
                    "path": str(path),
                    "chunk_number": i,
                },
            )
        )


def add_json_summary_record(
    records: List[Dict[str, Any]],
    path: Path,
    source_type: str,
    title: str,
) -> None:
    obj = read_json(path)
    if not obj:
        return

    summary = obj.get("summary", obj)
    text = json.dumps(summary, indent=2, ensure_ascii=False)

    records.append(
        make_record(
            chunk_id=source_type,
            source=str(path),
            source_type=source_type,
            title=title,
            text=text,
            metadata={
                "path": str(path),
                "json_section": "summary",
            },
        )
    )


def add_top_candidates(records: List[Dict[str, Any]], path: Path, limit: int = 50) -> None:
    rows = read_csv(path)
    for row in rows[:limit]:
        rank = row.get("exploratory_rank", "")
        rid = row.get("rnacentral_id", "")
        rna_type = row.get("rna_type", "")
        score = row.get("exploratory_rank_score", "")
        top3 = row.get("top3_similar_ids", "")

        text = (
            f"Candidato ncRNA rank {rank}: {rid}. "
            f"Tipo RNA: {rna_type}. "
            f"Score esplorativo: {score}. "
            f"Top 3 ncRNA simili: {top3}. "
            f"Nota: il ranking è esplorativo e non dimostra risposta alla siccità."
        )

        records.append(
            make_record(
                chunk_id=f"candidate_rank_{rank}",
                source=str(path),
                source_type="top_candidate",
                title=f"Candidato rank {rank} - {rid}",
                text=text,
                metadata=dict(row),
            )
        )


def add_ranking_rows(records: List[Dict[str, Any]], path: Path, limit: int = 250) -> None:
    rows = read_csv(path)

    for row in rows[:limit]:
        rank = row.get("exploratory_rank", "")
        rid = row.get("rnacentral_id", "")
        node_id = row.get("node_id", "")
        rna_type = row.get("rna_type", "")
        length = row.get("sequence_length", "")
        gc = row.get("gc_frac", "")
        score = row.get("exploratory_rank_score", "")
        mean_sim = row.get("mean_topk_similarity", "")
        top3 = row.get("top3_similar_ids", "")

        text = (
            f"ncRNA {rid} ({node_id}) ha rank esplorativo {rank}, "
            f"tipo RNA {rna_type}, lunghezza {length}, GC fraction {gc}, "
            f"score {score}, similarità media top-k {mean_sim}. "
            f"Top3 simili: {top3}. "
            f"Questo non è una prova di drought responsiveness."
        )

        records.append(
            make_record(
                chunk_id=f"ranking_{rank}",
                source=str(path),
                source_type="ranking_row",
                title=f"Ranking ncRNA {rid}",
                text=text,
                metadata=dict(row),
            )
        )


def add_project_status_record(records: List[Dict[str, Any]]) -> None:
    text = """
Stato del progetto:
La pipeline ha costruito sequenze ncRNA, metadata, expression summary conservativo,
archi di conoscenza, grafo, feature ncRNA, dataset PyTorch Geometric, EDA,
embedding non supervisionati, similarità tra ncRNA, ranking esplorativo e report metodologico.
Le label sono unknown (-1). Non esiste ancora prova funzionale di risposta alla siccità.
La chat LLM deve rispondere in modo prudente, distinguendo risultati computazionali
da evidenze biologiche validate.
"""

    records.append(
        make_record(
            chunk_id="project_status",
            source="internal",
            source_type="project_status",
            title="Stato metodologico del progetto",
            text=text,
            metadata={
                "important": True,
                "claim_level": "exploratory_unsupervised",
            },
        )
    )


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    output_jsonl = Path(args.output_jsonl)
    report_json = Path(args.report_json)
    summary_txt = Path(args.summary_txt)

    output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    report_json.parent.mkdir(parents=True, exist_ok=True)
    summary_txt.parent.mkdir(parents=True, exist_ok=True)

    records: List[Dict[str, Any]] = []

    add_project_status_record(records)

    text_sources = [
        (
            Path("outputs/reports/methodological_report.md"),
            "methodological_report",
            "Report metodologico completo",
        ),
        (
            Path("outputs/reports/methodological_report_summary.txt"),
            "methodological_summary",
            "Summary report metodologico",
        ),
        (
            Path("outputs/eda/graph_eda_summary.txt"),
            "eda_summary",
            "Summary EDA grafo e feature",
        ),
        (
            Path("outputs/embeddings/unsupervised_embeddings_summary.txt"),
            "embeddings_summary",
            "Summary embedding non supervisionati",
        ),
        (
            Path("outputs/similarity/ncrna_similarity_summary.txt"),
            "similarity_summary",
            "Summary similarità ncRNA",
        ),
        (
            Path("outputs/ranking/ncrna_exploratory_ranking_summary.txt"),
            "ranking_summary",
            "Summary ranking esplorativo",
        ),
    ]

    for path, source_type, title in text_sources:
        add_text_file_chunks(
            records=records,
            path=path,
            source_type=source_type,
            title=title,
            chunk_size=args.chunk_size,
            chunk_overlap=args.chunk_overlap,
        )

    json_sources = [
        (Path("outputs/maize_sequences_report.json"), "json_sequences", "Report sequenze"),
        (Path("outputs/maize_expression_summary_report.json"), "json_expression", "Report expression summary"),
        (Path("outputs/maize_known_edges_report.json"), "json_known_edges", "Report archi iniziali"),
        (Path("outputs/maize_labels_report.json"), "json_labels", "Report label"),
        (Path("outputs/validation_report.json"), "json_validation", "Report validazione"),
        (Path("outputs/graph_tables_report.json"), "json_graph", "Report grafo"),
        (Path("outputs/ncrna_features_report.json"), "json_features", "Report feature ncRNA"),
        (Path("outputs/pyg_dataset_report.json"), "json_pyg", "Report PyTorch Geometric"),
        (Path("outputs/embeddings/unsupervised_embeddings_report.json"), "json_embeddings", "Report embedding"),
        (Path("outputs/similarity/ncrna_similarity_report.json"), "json_similarity", "Report similarità"),
        (Path("outputs/ranking/ncrna_exploratory_ranking_report.json"), "json_ranking", "Report ranking"),
        (Path("outputs/reports/methodological_report.json"), "json_methodological", "Report metodologico JSON"),
    ]

    for path, source_type, title in json_sources:
        add_json_summary_record(records, path, source_type, title)

    add_top_candidates(
        records,
        Path("data/processed/ranking/ncrna_exploratory_top50.csv"),
        limit=50,
    )

    add_ranking_rows(
        records,
        Path("data/processed/ranking/ncrna_exploratory_ranking.csv"),
        limit=250,
    )

    with output_jsonl.open("w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

    source_counter = Counter(rec["source_type"] for rec in records)

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": {
            "knowledge_base_jsonl": str(output_jsonl),
            "report_json": str(report_json),
            "summary_txt": str(summary_txt),
        },
        "summary": {
            "n_chunks": len(records),
            "source_type_distribution": dict(source_counter),
            "methodological_note": (
                "Knowledge base for RAG chat. It contains project outputs and conservative methodological notes."
            ),
        },
    }

    with report_json.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    lines = []
    lines.append("RAG KNOWLEDGE BASE SUMMARY")
    lines.append("=" * 80)
    lines.append(f"Knowledge base: {output_jsonl}")
    lines.append(f"Chunks: {len(records)}")
    lines.append(f"Source distribution: {dict(source_counter)}")
    lines.append("")
    lines.append("Nota:")
    lines.append(report["summary"]["methodological_note"])

    with summary_txt.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print("\n".join(lines))

    logging.info("Knowledge base scritta: %s", output_jsonl)
    logging.info("Report: %s", report_json)
    logging.info("Summary: %s", summary_txt)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
