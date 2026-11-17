#!/usr/bin/env python3
"""
13_build_methodological_report.py

Genera un report metodologico automatico del progetto
maize_drought_ncrna_gnn / AI-GNN per ncRNA del mais.

Il report riassume:
- dati pubblici scaricati;
- sequenze ncRNA;
- expression summary;
- grafo;
- label;
- feature;
- dataset PyTorch Geometric;
- EDA;
- embedding non supervisionati;
- similarità tra ncRNA;
- ranking esplorativo.

Nota metodologica:
il report non dichiara una scoperta biologica definitiva.
Documenta lo stato computazionale della pipeline e i limiti attuali.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Genera report metodologico automatico del progetto ncRNA/GNN."
    )
    parser.add_argument("--config", default="configs/config_public.yaml")
    parser.add_argument("--output-dir", default="outputs/reports")
    parser.add_argument("--top-candidates-csv", default="data/processed/ranking/ncrna_exploratory_top50.csv")
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
        logging.warning("Config YAML non trovato: %s", path)
        return {}
    with path.open("r", encoding="utf-8") as f:
        obj = yaml.safe_load(f)
    return obj if isinstance(obj, dict) else {}


def load_json_if_exists(path: Path) -> Dict[str, Any]:
    if not path.exists():
        logging.warning("JSON non trovato: %s", path)
        return {}
    try:
        with path.open("r", encoding="utf-8") as f:
            obj = json.load(f)
        return obj if isinstance(obj, dict) else {}
    except Exception as exc:
        logging.warning("Errore lettura JSON %s: %s", path, exc)
        return {}


def read_csv_if_exists(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        logging.warning("CSV non trovato: %s", path)
        return []
    try:
        with path.open("r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))
    except Exception as exc:
        logging.warning("Errore lettura CSV %s: %s", path, exc)
        return []


def write_csv(rows: List[Dict[str, Any]], path: Path, fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def file_status(path: Path, artifact_name: str) -> Dict[str, Any]:
    return {
        "artifact": artifact_name,
        "path": str(path),
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else 0,
        "size_kb": round(path.stat().st_size / 1024, 3) if path.exists() else 0,
    }


def get_summary(report: Dict[str, Any]) -> Dict[str, Any]:
    summary = report.get("summary", {})
    return summary if isinstance(summary, dict) else {}


def json_block(obj: Any) -> str:
    return "```json\n" + json.dumps(obj, indent=2, ensure_ascii=False) + "\n```"


def compact_json(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=False)


def add_section(lines: List[str], title: str) -> None:
    lines.append("")
    lines.append(f"## {title}")
    lines.append("")


def build_inventory() -> List[Dict[str, Any]]:
    items = [
        ("config_public_yaml", "configs/config_public.yaml"),
        ("rnacentral_expanded_json", "data/raw/rnacentral_maize_ncrna_expanded.json"),
        ("maize_sequences_fasta", "data/processed/maize_sequences.fasta"),
        ("maize_sequences_metadata_csv", "data/processed/maize_sequences_metadata.csv"),
        ("maize_expression_summary_csv", "data/processed/maize_expression_summary.csv"),
        ("maize_known_edges_csv", "data/processed/maize_known_edges.csv"),
        ("maize_labels_csv", "data/processed/maize_labels.csv"),
        ("graph_nodes_csv", "data/processed/graph_nodes.csv"),
        ("graph_edges_csv", "data/processed/graph_edges.csv"),
        ("ncrna_features_csv", "data/processed/ncrna_features.csv"),
        ("pyg_graph_pt", "data/processed/pyg/maize_ncrna_graph.pt"),
        ("node_embeddings_csv", "data/processed/embeddings/node_embeddings.csv"),
        ("ncrna_embeddings_csv", "data/processed/embeddings/ncrna_embeddings.csv"),
        ("similarity_topk_csv", "data/processed/similarity/ncrna_similarity_topk.csv"),
        ("similarity_unique_pairs_csv", "data/processed/similarity/ncrna_similarity_unique_pairs.csv"),
        ("ncrna_exploratory_ranking_csv", "data/processed/ranking/ncrna_exploratory_ranking.csv"),
        ("ncrna_exploratory_top50_csv", "data/processed/ranking/ncrna_exploratory_top50.csv"),
        ("validation_report_json", "outputs/validation_report.json"),
        ("graph_tables_report_json", "outputs/graph_tables_report.json"),
        ("ncrna_features_report_json", "outputs/ncrna_features_report.json"),
        ("pyg_dataset_report_json", "outputs/pyg_dataset_report.json"),
        ("eda_report_json", "outputs/eda/graph_eda_report.json"),
        ("embeddings_report_json", "outputs/embeddings/unsupervised_embeddings_report.json"),
        ("similarity_report_json", "outputs/similarity/ncrna_similarity_report.json"),
        ("ranking_report_json", "outputs/ranking/ncrna_exploratory_ranking_report.json"),
    ]

    return [
        file_status(Path(path), artifact_name=name)
        for name, path in items
    ]


def load_all_reports() -> Dict[str, Dict[str, Any]]:
    report_paths = {
        "downloads": "outputs/downloads_report.json",
        "raw_inspection": "outputs/raw_downloads_inspection.json",
        "json_sources_profile": "outputs/json_sources_profile.json",
        "sequences": "outputs/maize_sequences_report.json",
        "expression": "outputs/maize_expression_summary_report.json",
        "known_edges": "outputs/maize_known_edges_report.json",
        "labels": "outputs/maize_labels_report.json",
        "validation": "outputs/validation_report.json",
        "graph_tables": "outputs/graph_tables_report.json",
        "features": "outputs/ncrna_features_report.json",
        "pyg_dataset": "outputs/pyg_dataset_report.json",
        "eda": "outputs/eda/graph_eda_report.json",
        "embeddings": "outputs/embeddings/unsupervised_embeddings_report.json",
        "similarity": "outputs/similarity/ncrna_similarity_report.json",
        "ranking": "outputs/ranking/ncrna_exploratory_ranking_report.json",
    }

    return {
        name: load_json_if_exists(Path(path))
        for name, path in report_paths.items()
    }


def copy_top_candidates_for_review(
    top_candidates: List[Dict[str, str]],
    output_path: Path,
) -> None:
    if not top_candidates:
        return

    rows = []
    for row in top_candidates:
        new_row = dict(row)
        new_row["review_action"] = (
            "Verificare letteratura, eventuale conservazione, presenza in dataset RNA-seq "
            "e possibile associazione indiretta con risposta alla siccità."
        )
        new_row["interpretation_limit"] = (
            "Ranking esplorativo non supervisionato; non è prova di drought responsiveness."
        )
        rows.append(new_row)

    fieldnames = list(rows[0].keys())
    write_csv(rows, output_path, fieldnames)


def build_markdown_report(
    cfg: Dict[str, Any],
    reports: Dict[str, Dict[str, Any]],
    inventory: List[Dict[str, Any]],
    top_candidates: List[Dict[str, str]],
) -> str:
    project = cfg.get("project", {})
    scope = cfg.get("scope", {})
    species = scope.get("species", {}) if isinstance(scope, dict) else {}

    validation_summary = get_summary(reports.get("validation", {}))
    sequence_summary = get_summary(reports.get("sequences", {}))
    expression_summary = get_summary(reports.get("expression", {}))
    graph_summary = get_summary(reports.get("graph_tables", {}))
    features_summary = get_summary(reports.get("features", {}))
    pyg_summary = get_summary(reports.get("pyg_dataset", {}))
    eda_summary = get_summary(reports.get("eda", {}))
    embeddings_summary = get_summary(reports.get("embeddings", {}))
    similarity_summary = get_summary(reports.get("similarity", {}))
    ranking_summary = get_summary(reports.get("ranking", {}))
    labels_summary = get_summary(reports.get("labels", {}))
    known_edges_summary = get_summary(reports.get("known_edges", {}))

    now = datetime.now(timezone.utc).isoformat()

    lines: List[str] = []

    lines.append("# Report metodologico automatico")
    lines.append("")
    lines.append(f"**Generato il:** {now}")
    lines.append("")
    lines.append("**Titolo operativo del progetto:**")
    lines.append("")
    lines.append("> Framework AI/GNN per identificare e prioritizzare ncRNA coinvolti nella risposta alla siccità nel mais")
    lines.append("")

    lines.append("## 1. Obiettivo del lavoro")
    lines.append("")
    lines.append(
        "L'obiettivo della pipeline è costruire un framework computazionale basato su dati pubblici, "
        "rappresentazione a grafo, feature di sequenza, embedding non supervisionati e ranking esplorativo "
        "per prioritizzare ncRNA di *Zea mays* da sottoporre a successiva validazione biologica."
    )
    lines.append("")
    lines.append(
        "In questa fase il sistema non produce ancora una classificazione supervisionata della risposta alla siccità, "
        "perché non sono disponibili label positive/negative validate o dati quantitativi RNA-seq associati direttamente "
        "ai singoli ncRNA."
    )

    add_section(lines, "2. Configurazione del progetto")

    lines.append("Informazioni principali derivate dal file di configurazione:")
    lines.append("")
    lines.append(json_block({
        "project": project,
        "species": species,
    }))

    add_section(lines, "3. Dati pubblici e sorgenti")

    downloads_summary = get_summary(reports.get("downloads", {}))
    json_profile_summary = get_summary(reports.get("json_sources_profile", {}))

    lines.append("La pipeline usa sorgenti pubbliche e dati scaricabili o interrogabili automaticamente.")
    lines.append("")
    lines.append("Sintesi del download, quando disponibile:")
    lines.append("")
    lines.append(json_block(downloads_summary))
    lines.append("")
    lines.append("Profilazione delle sorgenti JSON, quando disponibile:")
    lines.append("")
    lines.append(json_block(json_profile_summary))

    add_section(lines, "4. Dataset di sequenze ncRNA")

    lines.append(
        "Le sequenze ncRNA sono state ottenute da RNAcentral per *Zea mays* e convertite in FASTA e metadata tabellari."
    )
    lines.append("")
    lines.append(json_block(sequence_summary))

    add_section(lines, "5. Expression summary")

    lines.append(
        "Lo step di expression summary è intenzionalmente conservativo: non inventa valori di espressione, "
        "log2 fold-change o p-value. In assenza di quantificazione ncRNA-specifica, le righe restano marcate come "
        "`not_quantified_yet`."
    )
    lines.append("")
    lines.append(json_block(expression_summary))

    add_section(lines, "6. Archi di conoscenza e grafo")

    lines.append(
        "Il grafo iniziale integra archi di provenienza, annotazione e contesto: tipo RNA, database di origine, specie, "
        "stato di espressione, dataset candidati, record di letteratura e annotazioni GO. Questi archi non dimostrano "
        "automaticamente una funzione nella risposta alla siccità."
    )
    lines.append("")
    lines.append("Sintesi degli archi iniziali:")
    lines.append("")
    lines.append(json_block(known_edges_summary))
    lines.append("")
    lines.append("Sintesi delle tabelle del grafo:")
    lines.append("")
    lines.append(json_block(graph_summary))

    add_section(lines, "7. Label")

    lines.append(
        "Le label sono mantenute come `-1`, cioè unknown/unlabeled. Questo evita di trasformare annotazioni di database "
        "o similarità computazionale in evidenza biologica sperimentale."
    )
    lines.append("")
    lines.append(json_block(labels_summary))

    add_section(lines, "8. Validazione formale degli output")

    lines.append(
        "La validazione controlla presenza dei file, coerenza degli ID RNAcentral, duplicati, caratteri FASTA, "
        "coerenza tra sequenze, expression summary, grafo e label."
    )
    lines.append("")
    lines.append(json_block(validation_summary))

    add_section(lines, "9. Feature dei nodi ncRNA")

    lines.append(
        "La matrice feature dei nodi ncRNA combina proprietà di sequenza, composizione nucleotidica, tipo RNA, "
        "informazioni di grafo e stato di espressione. È adatta per analisi esplorativa, embedding e futuri modelli GNN."
    )
    lines.append("")
    lines.append(json_block(features_summary))

    add_section(lines, "10. Dataset PyTorch Geometric")

    lines.append(
        "Le tabelle di nodi, archi e feature sono state convertite in un oggetto PyTorch Geometric, contenente `x`, "
        "`edge_index`, `edge_attr`, `y`, `train_mask` e `unlabeled_mask`."
    )
    lines.append("")
    lines.append(json_block(pyg_summary))

    add_section(lines, "11. Analisi esplorativa")

    lines.append(
        "L'analisi esplorativa riassume distribuzioni di nodi, archi, tipi RNA, feature numeriche e produce grafici "
        "diagnostici. Non valuta performance supervisionata perché le label sono ancora unknown."
    )
    lines.append("")
    lines.append(json_block(eda_summary))

    add_section(lines, "12. Embedding non supervisionati")

    lines.append(
        "Gli embedding non supervisionati sono stati ottenuti standardizzando le feature, propagandole lungo il grafo "
        "e applicando riduzione dimensionale. Gli embedding servono per similarità, clustering e ranking esplorativo."
    )
    lines.append("")
    lines.append(json_block(embeddings_summary))

    add_section(lines, "13. Similarità tra ncRNA")

    lines.append(
        "La similarità tra ncRNA è calcolata sugli embedding tramite similarità coseno. Similarità elevata indica "
        "somiglianza computazionale, non prova di funzione biologica condivisa."
    )
    lines.append("")
    lines.append(json_block(similarity_summary))

    add_section(lines, "14. Ranking esplorativo dei candidati")

    lines.append(
        "Il ranking esplorativo combina similarità locale, centralità nello spazio embedding, priorità del tipo RNA, "
        "tipicità GC e qualità della sequenza. Il ranking serve a selezionare candidati per revisione biologica."
    )
    lines.append("")
    lines.append(json_block(ranking_summary))

    lines.append("")
    lines.append("### Primi candidati da revisionare")
    lines.append("")

    if top_candidates:
        lines.append("| Rank | RNAcentral ID | Tipo RNA | Score | Top3 simili |")
        lines.append("|---:|---|---|---:|---|")
        for row in top_candidates[:20]:
            lines.append(
                f"| {row.get('exploratory_rank', '')} "
                f"| {row.get('rnacentral_id', '')} "
                f"| {row.get('rna_type', '')} "
                f"| {row.get('exploratory_rank_score', '')} "
                f"| {row.get('top3_similar_ids', '')} |"
            )
    else:
        lines.append("Nessun candidato top disponibile.")

    add_section(lines, "15. Inventario degli output")

    existing = sum(1 for item in inventory if item["exists"])
    missing = sum(1 for item in inventory if not item["exists"])

    lines.append(f"File tracciati: **{len(inventory)}**")
    lines.append("")
    lines.append(f"- Presenti: **{existing}**")
    lines.append(f"- Mancanti: **{missing}**")
    lines.append("")
    lines.append("| Artifact | Path | Exists | Size KB |")
    lines.append("|---|---|---:|---:|")
    for item in inventory:
        lines.append(
            f"| {item['artifact']} | `{item['path']}` | {item['exists']} | {item['size_kb']} |"
        )

    add_section(lines, "16. Limiti metodologici attuali")

    lines.append(
        "1. Le label restano unknown (`-1`), quindi non è ancora possibile addestrare o valutare un classificatore "
        "supervisionato di risposta alla siccità."
    )
    lines.append(
        "2. I valori di espressione non sono ancora quantitativi: non sono disponibili log2 fold-change o p-value "
        "specifici per ciascun ncRNA."
    )
    lines.append(
        "3. Il grafo contiene relazioni di annotazione, provenienza, letteratura e contesto; non contiene ancora "
        "interazioni ncRNA-target validate."
    )
    lines.append(
        "4. La similarità embedding e il ranking esplorativo non dimostrano una funzione biologica comune."
    )
    lines.append(
        "5. I candidati prioritari devono essere verificati con letteratura, dataset RNA-seq, annotazioni funzionali "
        "e possibilmente validazione sperimentale."
    )

    add_section(lines, "17. Prossimi step consigliati")

    lines.append(
        "1. Collegare i candidati prioritari alla letteratura e ai dataset pubblici disponibili."
    )
    lines.append(
        "2. Integrare dati RNA-seq reali, quando disponibili, per ottenere evidenze quantitative."
    )
    lines.append(
        "3. Costruire proxy label più forti solo quando supportate da evidenze tracciabili."
    )
    lines.append(
        "4. Preparare una versione eterogenea del grafo con nodi gene, GO term, dataset, pubblicazioni e stress."
    )
    lines.append(
        "5. Implementare un modello GNN non supervisionato o semi-supervisionato per embedding più robusti."
    )
    lines.append(
        "6. Preparare una bozza di manoscritto focalizzata su pipeline, tracciabilità, limiti e casi candidati."
    )

    add_section(lines, "18. Conclusione")

    lines.append(
        "La pipeline ha prodotto un dataset computazionale coerente per ncRNA del mais, una rappresentazione a grafo, "
        "feature tabellari, dataset PyTorch Geometric, embedding non supervisionati, similarità tra ncRNA e ranking "
        "esplorativo. Lo stato attuale è adatto a una relazione metodologica e a una successiva fase di integrazione "
        "biologica, ma non costituisce ancora una dimostrazione funzionale di risposta alla siccità."
    )

    return "\n".join(lines) + "\n"


def build_plain_text_from_markdown(md: str) -> str:
    text = md.replace("# ", "")
    text = text.replace("## ", "")
    text = text.replace("### ", "")
    text = text.replace("**", "")
    text = text.replace("`", "")
    text = text.replace("> ", "")
    return text


def main() -> int:
    args = parse_args()
    setup_logging(args.verbose)

    config_path = Path(args.config)
    output_dir = Path(args.output_dir)
    top_candidates_path = Path(args.top_candidates_csv)

    output_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_yaml(config_path)

    logging.info("Caricamento report JSON precedenti")
    reports = load_all_reports()

    logging.info("Costruzione inventario output")
    inventory = build_inventory()

    logging.info("Caricamento top candidates: %s", top_candidates_path)
    top_candidates = read_csv_if_exists(top_candidates_path)

    logging.info("Scrittura inventario CSV")
    inventory_csv = output_dir / "pipeline_outputs_inventory.csv"
    write_csv(
        inventory,
        inventory_csv,
        ["artifact", "path", "exists", "size_bytes", "size_kb"],
    )

    logging.info("Scrittura top candidates per review")
    top_review_csv = output_dir / "top_candidates_for_review.csv"
    copy_top_candidates_for_review(top_candidates, top_review_csv)

    logging.info("Costruzione report markdown")
    md = build_markdown_report(
        cfg=cfg,
        reports=reports,
        inventory=inventory,
        top_candidates=top_candidates,
    )

    txt = build_plain_text_from_markdown(md)

    report_md = output_dir / "methodological_report.md"
    report_txt = output_dir / "methodological_report.txt"
    report_json = output_dir / "methodological_report.json"
    summary_txt = output_dir / "methodological_report_summary.txt"

    with report_md.open("w", encoding="utf-8") as f:
        f.write(md)

    with report_txt.open("w", encoding="utf-8") as f:
        f.write(txt)

    validation_summary = get_summary(reports.get("validation", {}))
    ranking_summary = get_summary(reports.get("ranking", {}))
    pyg_summary = get_summary(reports.get("pyg_dataset", {}))
    embeddings_summary = get_summary(reports.get("embeddings", {}))
    similarity_summary = get_summary(reports.get("similarity", {}))

    methodological_json = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_used": str(config_path),
        "outputs": {
            "methodological_report_md": str(report_md),
            "methodological_report_txt": str(report_txt),
            "methodological_report_json": str(report_json),
            "methodological_report_summary_txt": str(summary_txt),
            "pipeline_outputs_inventory_csv": str(inventory_csv),
            "top_candidates_for_review_csv": str(top_review_csv),
        },
        "summary": {
            "inventory_total": len(inventory),
            "inventory_existing": sum(1 for item in inventory if item["exists"]),
            "inventory_missing": sum(1 for item in inventory if not item["exists"]),
            "validation_status": validation_summary.get("validation_status", "unknown"),
            "pyg_n_nodes": pyg_summary.get("n_nodes", "NA"),
            "pyg_n_edges": pyg_summary.get("n_edges", "NA"),
            "embedding_dim_created": embeddings_summary.get("embedding_dim_created", "NA"),
            "similarity_n_ncrna": similarity_summary.get("n_ncrna", "NA"),
            "ranking_n_ranked_ncrna": ranking_summary.get("n_ranked_ncrna", "NA"),
            "ranking_top_n": ranking_summary.get("top_n", "NA"),
            "methodological_note": (
                "Automatic methodological report. It documents the computational state of the "
                "pipeline and does not claim validated drought-responsive ncRNAs."
            ),
        },
    }

    with report_json.open("w", encoding="utf-8") as f:
        json.dump(methodological_json, f, indent=2, ensure_ascii=False)

    summary_lines = []
    summary_lines.append("METHODOLOGICAL REPORT SUMMARY")
    summary_lines.append("=" * 80)
    summary_lines.append(f"Markdown report: {report_md}")
    summary_lines.append(f"Text report: {report_txt}")
    summary_lines.append(f"JSON report: {report_json}")
    summary_lines.append(f"Inventory CSV: {inventory_csv}")
    summary_lines.append(f"Top candidates review CSV: {top_review_csv}")
    summary_lines.append("")
    summary_lines.append(f"Inventory total: {methodological_json['summary']['inventory_total']}")
    summary_lines.append(f"Inventory existing: {methodological_json['summary']['inventory_existing']}")
    summary_lines.append(f"Inventory missing: {methodological_json['summary']['inventory_missing']}")
    summary_lines.append(f"Validation status: {methodological_json['summary']['validation_status']}")
    summary_lines.append(f"PyG nodes: {methodological_json['summary']['pyg_n_nodes']}")
    summary_lines.append(f"PyG edges: {methodological_json['summary']['pyg_n_edges']}")
    summary_lines.append(f"Embedding dim: {methodological_json['summary']['embedding_dim_created']}")
    summary_lines.append(f"Similarity ncRNA: {methodological_json['summary']['similarity_n_ncrna']}")
    summary_lines.append(f"Ranked ncRNA: {methodological_json['summary']['ranking_n_ranked_ncrna']}")
    summary_lines.append("")
    summary_lines.append("Nota:")
    summary_lines.append(methodological_json["summary"]["methodological_note"])

    with summary_txt.open("w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")

    print("\n".join(summary_lines))

    logging.info("Report Markdown: %s", report_md)
    logging.info("Report TXT: %s", report_txt)
    logging.info("Report JSON: %s", report_json)
    logging.info("Summary TXT: %s", summary_txt)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
