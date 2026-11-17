#!/usr/bin/env python3
"""
03_build_maize_known_edges.py

Costruisce un primo file di archi/relazioni biologiche e contestuali
per il progetto maize_drought_ncrna_gnn.

Nota metodologica:
questo script NON inventa interazioni ncRNA-gene e NON dichiara che tutti
gli ncRNA siano drought-responsive. Costruisce invece archi tracciabili:
- ncRNA -> RNA type
- ncRNA -> RNAcentral
- ncRNA -> specie Zea mays
- ncRNA -> expression status
- drought context -> GEO dataset candidates
- drought context -> literature records
- gene product -> GO term da QuickGO

Output:
- data/processed/maize_known_edges.csv
- outputs/maize_known_edges_report.json
"""

from __future__ import annotations

import argparse
import csv
import hashlib
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
        description="Costruisce maize_known_edges.csv per grafo biologico mais/siccità."
    )
    parser.add_argument(
        "--config",
        default="configs/config_public.yaml",
        help="File YAML di configurazione.",
    )
    parser.add_argument(
        "--metadata-csv",
        default="data/processed/maize_sequences_metadata.csv",
        help="Metadata CSV delle sequenze ncRNA.",
    )
    parser.add_argument(
        "--expression-csv",
        default="data/processed/maize_expression_summary.csv",
        help="Expression summary CSV.",
    )
    parser.add_argument(
        "--geo-json",
        default="data/raw/ncbi_geo_maize_drought_search.json",
        help="JSON NCBI GEO search.",
    )
    parser.add_argument(
        "--europe-pmc-json",
        default="data/raw/europe_pmc_maize_drought_ncrna.json",
        help="JSON Europe PMC.",
    )
    parser.add_argument(
        "--quickgo-json",
        default="data/raw/quickgo_maize_drought_annotations.json",
        help="JSON QuickGO.",
    )
    parser.add_argument(
        "--output-csv",
        default=None,
        help="CSV archi di output. Se omesso usa il percorso del YAML.",
    )
    parser.add_argument(
        "--report-json",
        default="outputs/maize_known_edges_report.json",
        help="Report JSON finale.",
    )
    parser.add_argument(
        "--max-literature-records",
        type=int,
        default=50,
        help="Numero massimo di pubblicazioni Europe PMC da trasformare in archi.",
    )
    parser.add_argument(
        "--max-quickgo-records",
        type=int,
        default=200,
        help="Numero massimo di record QuickGO da trasformare in archi.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Logging dettagliato.",
    )
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
        raise FileNotFoundError(f"Config YAML non trovato: {path}")
    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError("Config YAML non valido.")
    return cfg


def load_json_if_exists(path: Path) -> Dict[str, Any]:
    if not path.exists():
        logging.warning("File JSON non trovato: %s", path)
        return {}
    with path.open("r", encoding="utf-8") as f:
        obj = json.load(f)
    if not isinstance(obj, dict):
        logging.warning("Il file JSON non contiene un dizionario: %s", path)
        return {}
    return obj


def read_csv(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(f"CSV non trovato: {path}")
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        raise ValueError(f"CSV vuoto: {path}")
    return rows


def make_edge_id(
    source_node_id: str,
    target_node_id: str,
    edge_type: str,
    evidence_source: str,
) -> str:
    raw = f"{source_node_id}|{target_node_id}|{edge_type}|{evidence_source}"
    digest = hashlib.md5(raw.encode("utf-8")).hexdigest()
    return f"edge_{digest}"


def add_edge(
    edges: List[Dict[str, Any]],
    seen: set,
    source_node_id: str,
    source_node_type: str,
    target_node_id: str,
    target_node_type: str,
    edge_type: str,
    evidence_source: str,
    evidence_label: str,
    evidence_url: str = "",
    confidence_score: float = 1.0,
    weight: float = 1.0,
    is_direct_biological_evidence: int = 0,
    is_proxy_evidence: int = 0,
    notes: str = "",
) -> None:
    edge_id = make_edge_id(
        source_node_id=source_node_id,
        target_node_id=target_node_id,
        edge_type=edge_type,
        evidence_source=evidence_source,
    )

    if edge_id in seen:
        return

    seen.add(edge_id)

    edges.append(
        {
            "edge_id": edge_id,
            "source_node_id": source_node_id,
            "source_node_type": source_node_type,
            "target_node_id": target_node_id,
            "target_node_type": target_node_type,
            "edge_type": edge_type,
            "evidence_source": evidence_source,
            "evidence_label": evidence_label,
            "evidence_url": evidence_url,
            "confidence_score": f"{confidence_score:.3f}",
            "weight": f"{weight:.3f}",
            "is_direct_biological_evidence": is_direct_biological_evidence,
            "is_proxy_evidence": is_proxy_evidence,
            "notes": notes,
        }
    )


def build_ncrna_edges(
    metadata_rows: List[Dict[str, Any]],
    expression_rows: List[Dict[str, Any]],
    species_name: str,
    taxon_id: int,
) -> List[Dict[str, Any]]:
    edges: List[Dict[str, Any]] = []
    seen = set()

    expression_by_id = {
        row.get("rnacentral_id"): row
        for row in expression_rows
        if row.get("rnacentral_id")
    }

    species_node_id = f"species:NCBITaxon_{taxon_id}"
    database_node_id = "database:RNAcentral"

    for rec in metadata_rows:
        rid = rec.get("rnacentral_id")
        if not rid:
            continue

        ncrna_node_id = f"ncrna:{rid}"
        rna_type = rec.get("rna_type") or "unknown"
        rna_type_node_id = f"rna_type:{rna_type}"

        add_edge(
            edges,
            seen,
            source_node_id=ncrna_node_id,
            source_node_type="ncRNA",
            target_node_id=rna_type_node_id,
            target_node_type="rna_type",
            edge_type="has_rna_type",
            evidence_source="RNAcentral",
            evidence_label="RNA type annotation from RNAcentral record",
            confidence_score=0.95,
            weight=1.0,
            is_direct_biological_evidence=1,
            is_proxy_evidence=0,
            notes="Annotation edge; not a drought-response claim.",
        )

        add_edge(
            edges,
            seen,
            source_node_id=ncrna_node_id,
            source_node_type="ncRNA",
            target_node_id=database_node_id,
            target_node_type="database",
            edge_type="listed_in_database",
            evidence_source="RNAcentral",
            evidence_label="ncRNA sequence retrieved from RNAcentral",
            confidence_score=0.95,
            weight=1.0,
            is_direct_biological_evidence=1,
            is_proxy_evidence=0,
            notes="Source provenance edge.",
        )

        add_edge(
            edges,
            seen,
            source_node_id=ncrna_node_id,
            source_node_type="ncRNA",
            target_node_id=species_node_id,
            target_node_type="species",
            edge_type="belongs_to_species",
            evidence_source="RNAcentral",
            evidence_label=f"Taxon assignment to {species_name}",
            confidence_score=0.95,
            weight=1.0,
            is_direct_biological_evidence=1,
            is_proxy_evidence=0,
            notes="Species assignment edge.",
        )

        expr = expression_by_id.get(rid, {})
        expression_status = expr.get("expression_status", "unknown")
        expression_status_node_id = f"expression_status:{expression_status}"

        add_edge(
            edges,
            seen,
            source_node_id=ncrna_node_id,
            source_node_type="ncRNA",
            target_node_id=expression_status_node_id,
            target_node_type="expression_status",
            edge_type="has_expression_status",
            evidence_source="pipeline_expression_summary",
            evidence_label="Expression status from conservative summary",
            confidence_score=0.50,
            weight=0.5,
            is_direct_biological_evidence=0,
            is_proxy_evidence=1,
            notes="Does not indicate differential expression; quantitative expression not computed yet.",
        )

    return edges


def build_geo_edges(geo: Dict[str, Any]) -> List[Dict[str, Any]]:
    edges: List[Dict[str, Any]] = []
    seen = set()

    context_node_id = "context:maize_drought_public_expression_search"
    esearch = geo.get("esearchresult", {}) if isinstance(geo, dict) else {}
    ids = esearch.get("idlist", []) or []

    for geo_id in ids:
        dataset_node_id = f"geo_dataset:{geo_id}"

        add_edge(
            edges,
            seen,
            source_node_id=context_node_id,
            source_node_type="search_context",
            target_node_id=dataset_node_id,
            target_node_type="public_expression_dataset",
            edge_type="found_candidate_expression_dataset",
            evidence_source="NCBI_GEO",
            evidence_label="Dataset returned by public maize drought RNA-seq/transcriptome search",
            evidence_url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}",
            confidence_score=0.70,
            weight=0.7,
            is_direct_biological_evidence=0,
            is_proxy_evidence=1,
            notes="Search result only; dataset must be curated before expression quantification.",
        )

    return edges


def get_publication_id(rec: Dict[str, Any], index: int) -> str:
    for key in ["pmid", "pmcid", "doi", "id"]:
        value = rec.get(key)
        if value:
            safe = str(value).replace(" ", "_").replace("/", "_")
            return safe
    return f"europepmc_record_{index}"


def build_literature_edges(
    europe: Dict[str, Any],
    max_records: int,
) -> List[Dict[str, Any]]:
    edges: List[Dict[str, Any]] = []
    seen = set()

    context_node_id = "context:maize_drought_ncrna_literature_search"
    records = europe.get("resultList", {}).get("result", []) if isinstance(europe, dict) else []

    for i, rec in enumerate(records[:max_records], start=1):
        pub_id = get_publication_id(rec, i)
        publication_node_id = f"publication:{pub_id}"

        title = rec.get("title", "")
        year = rec.get("pubYear", "")
        journal = rec.get("journalTitle", "")
        doi = rec.get("doi", "")
        pmid = rec.get("pmid", "")

        if doi:
            url = f"https://doi.org/{doi}"
        elif pmid:
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        else:
            url = ""

        label = title[:250] if title else "Europe PMC literature record"

        add_edge(
            edges,
            seen,
            source_node_id=context_node_id,
            source_node_type="search_context",
            target_node_id=publication_node_id,
            target_node_type="publication",
            edge_type="retrieved_literature_record",
            evidence_source="Europe_PMC",
            evidence_label=label,
            evidence_url=url,
            confidence_score=0.70,
            weight=0.7,
            is_direct_biological_evidence=0,
            is_proxy_evidence=1,
            notes=f"Publication search result; year={year}; journal={journal}. Manual validation required.",
        )

    return edges


def build_quickgo_edges(
    quickgo: Dict[str, Any],
    max_records: int,
) -> List[Dict[str, Any]]:
    edges: List[Dict[str, Any]] = []
    seen = set()

    results = quickgo.get("results", []) if isinstance(quickgo, dict) else []

    for rec in results[:max_records]:
        gene_product_id = rec.get("geneProductId")
        go_id = rec.get("goId")
        qualifier = rec.get("qualifier", "")
        evidence = rec.get("goEvidence", "")
        aspect = rec.get("goAspect", "")
        assigned_by = rec.get("assignedBy", "")
        taxon_id = rec.get("taxonId", "")

        if not gene_product_id or not go_id:
            continue

        gene_node_id = f"gene_product:{gene_product_id}"
        go_node_id = f"go_term:{go_id}"

        add_edge(
            edges,
            seen,
            source_node_id=gene_node_id,
            source_node_type="gene_product",
            target_node_id=go_node_id,
            target_node_type="GO_term",
            edge_type="annotated_with_go_term",
            evidence_source="QuickGO",
            evidence_label=f"{qualifier} {go_id}; evidence={evidence}; aspect={aspect}",
            evidence_url=f"https://www.ebi.ac.uk/QuickGO/term/{go_id}",
            confidence_score=0.80 if evidence != "IEA" else 0.55,
            weight=0.80 if evidence != "IEA" else 0.55,
            is_direct_biological_evidence=1,
            is_proxy_evidence=0 if evidence != "IEA" else 1,
            notes=f"QuickGO annotation; assignedBy={assigned_by}; taxonId={taxon_id}.",
        )

    return edges


def write_edges_csv(edges: List[Dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "edge_id",
        "source_node_id",
        "source_node_type",
        "target_node_id",
        "target_node_type",
        "edge_type",
        "evidence_source",
        "evidence_label",
        "evidence_url",
        "confidence_score",
        "weight",
        "is_direct_biological_evidence",
        "is_proxy_evidence",
        "notes",
    ]

    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for edge in edges:
            writer.writerow(edge)


def write_report(report: Dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    metadata_csv = Path(args.metadata_csv)
    expression_csv = Path(args.expression_csv)
    geo_json = Path(args.geo_json)
    europe_json = Path(args.europe_pmc_json)
    quickgo_json = Path(args.quickgo_json)
    report_json = Path(args.report_json)

    cfg = load_yaml(config_path)

    if args.output_csv:
        output_csv = Path(args.output_csv)
    else:
        output_csv = Path(
            cfg.get("pipeline_outputs", {})
            .get("initial_data", {})
            .get("maize_known_edges", "data/processed/maize_known_edges.csv")
        )

    species = cfg.get("scope", {}).get("species", {})
    species_name = species.get("scientific_name", "Zea mays")
    taxon_id = int(species.get("ncbi_taxon_id", 4577))

    logging.info("Caricamento metadata: %s", metadata_csv)
    metadata_rows = read_csv(metadata_csv)

    logging.info("Caricamento expression summary: %s", expression_csv)
    expression_rows = read_csv(expression_csv)

    logging.info("Caricamento GEO JSON: %s", geo_json)
    geo = load_json_if_exists(geo_json)

    logging.info("Caricamento Europe PMC JSON: %s", europe_json)
    europe = load_json_if_exists(europe_json)

    logging.info("Caricamento QuickGO JSON: %s", quickgo_json)
    quickgo = load_json_if_exists(quickgo_json)

    all_edges: List[Dict[str, Any]] = []

    ncrna_edges = build_ncrna_edges(
        metadata_rows=metadata_rows,
        expression_rows=expression_rows,
        species_name=species_name,
        taxon_id=taxon_id,
    )
    all_edges.extend(ncrna_edges)

    geo_edges = build_geo_edges(geo)
    all_edges.extend(geo_edges)

    literature_edges = build_literature_edges(
        europe=europe,
        max_records=args.max_literature_records,
    )
    all_edges.extend(literature_edges)

    quickgo_edges = build_quickgo_edges(
        quickgo=quickgo,
        max_records=args.max_quickgo_records,
    )
    all_edges.extend(quickgo_edges)

    write_edges_csv(all_edges, output_csv)

    edge_type_counter = Counter(edge["edge_type"] for edge in all_edges)
    source_type_counter = Counter(edge["source_node_type"] for edge in all_edges)
    target_type_counter = Counter(edge["target_node_type"] for edge in all_edges)
    evidence_source_counter = Counter(edge["evidence_source"] for edge in all_edges)

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "inputs": {
            "metadata_csv": str(metadata_csv),
            "expression_csv": str(expression_csv),
            "geo_json": str(geo_json),
            "europe_pmc_json": str(europe_json),
            "quickgo_json": str(quickgo_json),
        },
        "outputs": {
            "maize_known_edges": str(output_csv),
            "report_json": str(report_json),
        },
        "summary": {
            "n_edges_total": len(all_edges),
            "n_ncrna_edges": len(ncrna_edges),
            "n_geo_edges": len(geo_edges),
            "n_literature_edges": len(literature_edges),
            "n_quickgo_edges": len(quickgo_edges),
            "edge_type_distribution": dict(edge_type_counter),
            "source_node_type_distribution": dict(source_type_counter),
            "target_node_type_distribution": dict(target_type_counter),
            "evidence_source_distribution": dict(evidence_source_counter),
            "methodological_note": (
                "Edges are provenance, annotation, search-context, and GO annotation edges. "
                "They do not prove that every ncRNA is drought-responsive."
            ),
        },
    }

    write_report(report, report_json)

    logging.info("File archi creato: %s", output_csv)
    logging.info("Report creato: %s", report_json)
    logging.info("Archi totali: %s", len(all_edges))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
