#!/usr/bin/env python3
"""
06_build_graph_tables.py

Costruisce le tabelle nodi e archi del grafo biologico/contestuale
per il progetto maize_drought_ncrna_gnn.

Input:
- data/processed/maize_sequences_metadata.csv
- data/processed/maize_expression_summary.csv
- data/processed/maize_labels.csv
- data/processed/maize_known_edges.csv
- configs/config_public.yaml

Output:
- data/processed/graph_nodes.csv
- data/processed/graph_edges.csv
- outputs/graph_tables_report.json
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Costruisce graph_nodes.csv e graph_edges.csv."
    )
    parser.add_argument("--config", default="configs/config_public.yaml")
    parser.add_argument("--metadata-csv", default="data/processed/maize_sequences_metadata.csv")
    parser.add_argument("--expression-csv", default="data/processed/maize_expression_summary.csv")
    parser.add_argument("--labels-csv", default="data/processed/maize_labels.csv")
    parser.add_argument("--known-edges-csv", default="data/processed/maize_known_edges.csv")
    parser.add_argument("--nodes-csv", default="data/processed/graph_nodes.csv")
    parser.add_argument("--edges-csv", default="data/processed/graph_edges.csv")
    parser.add_argument("--report-json", default="outputs/graph_tables_report.json")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


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


def safe_int(value: Any, default: int = 0) -> int:
    try:
        if value in {"", "NA", None}:
            return default
        return int(float(str(value)))
    except Exception:
        return default


def safe_float_text(value: Any, default: str = "0.000") -> str:
    try:
        if value in {"", "NA", None}:
            return default
        return f"{float(value):.3f}"
    except Exception:
        return default


def make_name_from_node_id(node_id: str) -> str:
    if ":" in node_id:
        return node_id.split(":", 1)[1]
    return node_id


def add_node(
    nodes: Dict[str, Dict[str, Any]],
    node_id: str,
    node_type: str,
    name: str = "",
    source: str = "",
    rna_type: str = "",
    sequence_length: Any = "",
    label: Any = "",
    label_status: str = "",
    expression_status: str = "",
    attributes: Dict[str, Any] | None = None,
) -> None:
    if not node_id:
        return

    if attributes is None:
        attributes = {}

    if node_id not in nodes:
        nodes[node_id] = {
            "node_id": node_id,
            "node_type": node_type,
            "name": name or make_name_from_node_id(node_id),
            "source": source,
            "rna_type": rna_type,
            "sequence_length": sequence_length,
            "label": label,
            "label_status": label_status,
            "expression_status": expression_status,
            "attributes_json": json.dumps(attributes, ensure_ascii=False),
        }
    else:
        current = nodes[node_id]

        for key, value in {
            "node_type": node_type,
            "name": name,
            "source": source,
            "rna_type": rna_type,
            "sequence_length": sequence_length,
            "label": label,
            "label_status": label_status,
            "expression_status": expression_status,
        }.items():
            if value not in {"", None} and current.get(key) in {"", None}:
                current[key] = value

        try:
            old_attrs = json.loads(current.get("attributes_json", "{}"))
            old_attrs.update(attributes)
            current["attributes_json"] = json.dumps(old_attrs, ensure_ascii=False)
        except Exception:
            current["attributes_json"] = json.dumps(attributes, ensure_ascii=False)


def build_expression_by_id(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    return {
        row.get("rnacentral_id", ""): row
        for row in rows
        if row.get("rnacentral_id", "")
    }


def build_label_by_id(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    return {
        row.get("rnacentral_id", ""): row
        for row in rows
        if row.get("rnacentral_id", "")
    }


def build_nodes(
    metadata_rows: List[Dict[str, str]],
    expression_rows: List[Dict[str, str]],
    label_rows: List[Dict[str, str]],
    edge_rows: List[Dict[str, str]],
    cfg: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    nodes: Dict[str, Dict[str, Any]] = {}

    expression_by_id = build_expression_by_id(expression_rows)
    label_by_id = build_label_by_id(label_rows)

    species = cfg.get("scope", {}).get("species", {})
    scientific_name = species.get("scientific_name", "Zea mays")
    taxon_id = species.get("ncbi_taxon_id", 4577)

    add_node(
        nodes,
        node_id=f"species:NCBITaxon_{taxon_id}",
        node_type="species",
        name=scientific_name,
        source="config",
        attributes={
            "scientific_name": scientific_name,
            "common_name": species.get("common_name", "maize"),
            "italian_name": species.get("italian_name", "mais"),
            "ncbi_taxon_id": taxon_id,
        },
    )

    add_node(
        nodes,
        node_id="database:RNAcentral",
        node_type="database",
        name="RNAcentral",
        source="RNAcentral",
        attributes={"database_role": "ncRNA sequence source"},
    )

    for rec in metadata_rows:
        rid = rec.get("rnacentral_id", "")
        if not rid:
            continue

        expr = expression_by_id.get(rid, {})
        label = label_by_id.get(rid, {})

        node_id = f"ncrna:{rid}"
        rna_type = rec.get("rna_type", "unknown")
        sequence_length = rec.get("sequence_length", "")

        add_node(
            nodes,
            node_id=node_id,
            node_type="ncRNA",
            name=rid,
            source="RNAcentral",
            rna_type=rna_type,
            sequence_length=sequence_length,
            label=label.get("label", ""),
            label_status=label.get("label_status", ""),
            expression_status=expr.get("expression_status", ""),
            attributes={
                "rnacentral_id": rid,
                "md5": rec.get("md5", ""),
                "description": rec.get("description", ""),
                "url": rec.get("url", ""),
                "publications_count": rec.get("publications_count", ""),
                "xrefs_count": rec.get("xrefs_count", ""),
                "has_quantitative_expression": expr.get("has_quantitative_expression", ""),
                "has_drought_expression_evidence": expr.get("has_drought_expression_evidence", ""),
                "label_confidence": label.get("label_confidence", ""),
                "evidence_notes": label.get("evidence_notes", ""),
            },
        )

        add_node(
            nodes,
            node_id=f"rna_type:{rna_type}",
            node_type="rna_type",
            name=rna_type,
            source="RNAcentral",
            attributes={"category": "RNA type annotation"},
        )

        expr_status = expr.get("expression_status", "")
        if expr_status:
            add_node(
                nodes,
                node_id=f"expression_status:{expr_status}",
                node_type="expression_status",
                name=expr_status,
                source="pipeline_expression_summary",
                attributes={
                    "status_meaning": "quantitative expression status for ncRNA candidate"
                },
            )

    for edge in edge_rows:
        src = edge.get("source_node_id", "")
        tgt = edge.get("target_node_id", "")
        src_type = edge.get("source_node_type", "")
        tgt_type = edge.get("target_node_type", "")
        evidence_source = edge.get("evidence_source", "")

        add_node(
            nodes,
            node_id=src,
            node_type=src_type,
            name=make_name_from_node_id(src),
            source=evidence_source,
            attributes={"created_from_edge_endpoint": True},
        )

        add_node(
            nodes,
            node_id=tgt,
            node_type=tgt_type,
            name=make_name_from_node_id(tgt),
            source=evidence_source,
            attributes={"created_from_edge_endpoint": True},
        )

    return nodes


def normalize_edges(edge_rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    graph_edges = []

    for edge in edge_rows:
        graph_edges.append(
            {
                "edge_id": edge.get("edge_id", ""),
                "source_node_id": edge.get("source_node_id", ""),
                "target_node_id": edge.get("target_node_id", ""),
                "source_node_type": edge.get("source_node_type", ""),
                "target_node_type": edge.get("target_node_type", ""),
                "edge_type": edge.get("edge_type", ""),
                "weight": safe_float_text(edge.get("weight", "1.0"), "1.000"),
                "confidence_score": safe_float_text(edge.get("confidence_score", "1.0"), "1.000"),
                "evidence_source": edge.get("evidence_source", ""),
                "evidence_label": edge.get("evidence_label", ""),
                "evidence_url": edge.get("evidence_url", ""),
                "is_direct_biological_evidence": safe_int(
                    edge.get("is_direct_biological_evidence", 0), 0
                ),
                "is_proxy_evidence": safe_int(edge.get("is_proxy_evidence", 0), 0),
                "notes": edge.get("notes", ""),
            }
        )

    return graph_edges


def validate_graph(
    nodes: Dict[str, Dict[str, Any]],
    graph_edges: List[Dict[str, Any]],
) -> Dict[str, Any]:
    node_ids = set(nodes.keys())

    missing_sources = []
    missing_targets = []
    self_loops = 0

    for edge in graph_edges:
        src = edge.get("source_node_id", "")
        tgt = edge.get("target_node_id", "")

        if src not in node_ids:
            missing_sources.append(src)

        if tgt not in node_ids:
            missing_targets.append(tgt)

        if src == tgt:
            self_loops += 1

    errors = []
    warnings = []

    if missing_sources:
        errors.append(f"Sorgenti mancanti nei nodi: {sorted(set(missing_sources))[:10]}")

    if missing_targets:
        errors.append(f"Target mancanti nei nodi: {sorted(set(missing_targets))[:10]}")

    if self_loops:
        warnings.append(f"Self-loop presenti: {self_loops}")

    status = "pass" if not errors else "failed"

    return {
        "status": status,
        "errors": errors,
        "warnings": warnings,
        "missing_sources_count": len(missing_sources),
        "missing_targets_count": len(missing_targets),
        "self_loops": self_loops,
    }


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    metadata_path = Path(args.metadata_csv)
    expression_path = Path(args.expression_csv)
    labels_path = Path(args.labels_csv)
    known_edges_path = Path(args.known_edges_csv)
    nodes_path = Path(args.nodes_csv)
    edges_path = Path(args.edges_csv)
    report_path = Path(args.report_json)

    cfg = load_yaml(config_path)

    logging.info("Caricamento metadata: %s", metadata_path)
    metadata_rows = read_csv(metadata_path)

    logging.info("Caricamento expression summary: %s", expression_path)
    expression_rows = read_csv(expression_path)

    logging.info("Caricamento labels: %s", labels_path)
    label_rows = read_csv(labels_path)

    logging.info("Caricamento known edges: %s", known_edges_path)
    edge_rows = read_csv(known_edges_path)

    logging.info("Costruzione nodi")
    nodes = build_nodes(
        metadata_rows=metadata_rows,
        expression_rows=expression_rows,
        label_rows=label_rows,
        edge_rows=edge_rows,
        cfg=cfg,
    )

    logging.info("Normalizzazione archi")
    graph_edges = normalize_edges(edge_rows)

    validation = validate_graph(nodes, graph_edges)

    node_rows = list(nodes.values())

    node_fieldnames = [
        "node_id",
        "node_type",
        "name",
        "source",
        "rna_type",
        "sequence_length",
        "label",
        "label_status",
        "expression_status",
        "attributes_json",
    ]

    edge_fieldnames = [
        "edge_id",
        "source_node_id",
        "target_node_id",
        "source_node_type",
        "target_node_type",
        "edge_type",
        "weight",
        "confidence_score",
        "evidence_source",
        "evidence_label",
        "evidence_url",
        "is_direct_biological_evidence",
        "is_proxy_evidence",
        "notes",
    ]

    write_csv(node_rows, nodes_path, node_fieldnames)
    write_csv(graph_edges, edges_path, edge_fieldnames)

    node_type_counter = Counter(row["node_type"] for row in node_rows)
    edge_type_counter = Counter(row["edge_type"] for row in graph_edges)
    evidence_source_counter = Counter(row["evidence_source"] for row in graph_edges)

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "inputs": {
            "metadata_csv": str(metadata_path),
            "expression_csv": str(expression_path),
            "labels_csv": str(labels_path),
            "known_edges_csv": str(known_edges_path),
        },
        "outputs": {
            "graph_nodes_csv": str(nodes_path),
            "graph_edges_csv": str(edges_path),
            "report_json": str(report_path),
        },
        "summary": {
            "n_nodes": len(node_rows),
            "n_edges": len(graph_edges),
            "node_type_distribution": dict(node_type_counter),
            "edge_type_distribution": dict(edge_type_counter),
            "evidence_source_distribution": dict(evidence_source_counter),
            "validation": validation,
            "methodological_note": (
                "Graph tables are based on provenance, annotation, search-context "
                "and GO annotation edges. The graph is ready for feature construction, "
                "but supervised training requires stronger labels or expression evidence."
            ),
        },
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    logging.info("Nodi creati: %s", len(node_rows))
    logging.info("Archi creati: %s", len(graph_edges))
    logging.info("File nodi: %s", nodes_path)
    logging.info("File archi: %s", edges_path)
    logging.info("Report: %s", report_path)
    logging.info("Validazione grafo: %s", validation["status"])

    return 0 if validation["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
