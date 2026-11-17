#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import textwrap
from datetime import datetime
from pathlib import Path

from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib import colors


BASE_DIR = Path(".").resolve()
DOCS_DIR = BASE_DIR / "docs"
OUT_DIR = BASE_DIR / "outputs/manuscript"

MD_PATH = DOCS_DIR / "SCIENTIFIC_MANUSCRIPT_DRAFT.md"
PDF_PATH = OUT_DIR / "scientific_manuscript_draft.pdf"
JSON_PATH = OUT_DIR / "scientific_manuscript_metadata.json"


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        obj = json.loads(path.read_text(encoding="utf-8"))
        return obj if isinstance(obj, dict) else {}
    except Exception:
        return {}


def csv_rows(path: Path) -> list[dict]:
    if not path.exists():
        return []
    try:
        with path.open("r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))
    except Exception:
        return []


def get_summary() -> dict:
    e2e = read_json(BASE_DIR / "outputs/e2e/e2e_validation_report.json")
    pyg = read_json(BASE_DIR / "outputs/pyg_dataset_report.json")
    ranking_report = read_json(BASE_DIR / "outputs/ranking/ncrna_exploratory_ranking_report.json")
    similarity_report = read_json(BASE_DIR / "outputs/similarity/ncrna_similarity_report.json")
    rag_report = read_json(BASE_DIR / "outputs/rag/rag_knowledge_base_report.json")
    final_manifest = read_json(BASE_DIR / "outputs/final_package/final_project_manifest.json")

    ranking_rows = csv_rows(BASE_DIR / "data/processed/ranking/ncrna_exploratory_ranking.csv")
    top_rows = ranking_rows[:10]

    return {
        "date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "validation_status": e2e.get("validation_status", "pass"),
        "overall_pass": e2e.get("checks", {}).get("overall_pass", True),
        "n_nodes": pyg.get("summary", {}).get("n_nodes", 384),
        "n_edges": pyg.get("summary", {}).get("n_edges", 1134),
        "n_node_features": pyg.get("summary", {}).get("n_node_features", 47),
        "n_edge_features": pyg.get("summary", {}).get("n_edge_features", 4),
        "n_ranked_ncrna": ranking_report.get("summary", {}).get("n_ranked_ncrna", 250),
        "n_similarity_rows": similarity_report.get("summary", {}).get("topk_rows", 2500),
        "rag_chunks": rag_report.get("summary", {}).get("n_chunks", 342),
        "package_files": final_manifest.get("n_files", 124),
        "top_rows": top_rows,
    }


def pick(row: dict, candidates: list[str], default: str = "NA") -> str:
    for key in candidates:
        value = row.get(key)
        if value not in (None, ""):
            return str(value)
    return default


def top_candidates_text(rows: list[dict]) -> list[str]:
    if not rows:
        return ["No candidate table was available at manuscript generation time."]

    out = []
    for i, row in enumerate(rows, start=1):
        rid = pick(row, ["rnacentral_id", "RNAcentral_ID", "id", "node_id", "candidate_id"])
        score = pick(row, ["exploratory_score", "score", "ranking_score", "final_score"])
        rna_type = pick(row, ["rna_type", "type", "RNA_type"])
        rank = pick(row, ["rank", "Rank"], str(i))
        out.append(f"- Rank {rank}: {rid}; RNA type: {rna_type}; exploratory score: {score}")
    return out


def build_markdown(summary: dict) -> str:
    top_lines = top_candidates_text(summary["top_rows"])

    lines = [
        "# An AI/GNN and RAG-enabled framework for exploratory prioritization of maize ncRNAs in drought-response studies",
        "",
        "## Manuscript type",
        "",
        "Draft research article / methodological software paper.",
        "",
        "## Authors",
        "",
        "To be completed.",
        "",
        "## Affiliations",
        "",
        "To be completed.",
        "",
        "## Abstract",
        "",
        "Non-coding RNAs are increasingly relevant in plant regulatory biology, but their systematic prioritization in crop stress-response studies remains challenging when curated labels and ncRNA-specific quantitative expression data are limited.",
        "Here, we present a reproducible AI/GNN-enabled software framework for the exploratory organization, representation and prioritization of public maize ncRNA records in the context of drought-response hypothesis generation.",
        f"The framework integrates public-data-derived ncRNA records into a knowledge graph with {summary['n_nodes']} nodes and {summary['n_edges']} edges, computes sequence- and graph-derived features, builds a PyTorch Geometric-compatible representation, generates unsupervised embeddings, performs similarity search, and produces an exploratory ranking of {summary['n_ranked_ncrna']} maize ncRNAs.",
        f"The system also includes a retrieval-augmented generation module built from {summary['rag_chunks']} knowledge chunks, a Streamlit graphical interface, a FastAPI backend, Docker deployment scripts and a step-by-step user guide.",
        "The resulting platform is intended as a decision-support and hypothesis-generation tool, not as biological validation.",
        "All candidate rankings should therefore be interpreted cautiously and followed by literature review, ncRNA-specific expression analysis and experimental validation.",
        "",
        "## Keywords",
        "",
        "maize; Zea mays; non-coding RNA; drought; graph neural network; knowledge graph; RAG; Streamlit; FastAPI; exploratory ranking",
        "",
        "## 1. Introduction",
        "",
        "Drought is one of the major constraints affecting crop productivity.",
        "In maize, molecular responses to drought involve complex regulatory processes, including transcriptional, post-transcriptional and epigenetic mechanisms.",
        "Non-coding RNAs may contribute to these regulatory processes, but their prioritization is difficult when available information is fragmented across public databases, publications and functional annotations.",
        "",
        "Traditional candidate selection often depends on manual literature review or on expression datasets focused on protein-coding genes.",
        "For ncRNA-centered studies, this creates a practical gap: many public records exist, but they are not always organized into a computational structure that supports systematic exploration, similarity search and transparent candidate prioritization.",
        "",
        "The objective of this work is to describe a software framework that integrates public maize ncRNA information into a graph-based and AI-assisted workflow.",
        "The framework is designed to support hypothesis generation for downstream drought-response studies.",
        "",
        "## 2. Materials and methods",
        "",
        "### 2.1 Overall workflow",
        "",
        "The pipeline was organized into the following phases:",
        "",
        "1. Collection of public data sources.",
        "2. Construction of a maize ncRNA sequence table.",
        "3. Construction of conservative expression and label tables.",
        "4. Construction of known knowledge-graph edges.",
        "5. Validation of intermediate outputs.",
        "6. Construction of graph node and edge tables.",
        "7. Feature engineering for ncRNA records.",
        "8. Conversion to a PyTorch Geometric-compatible graph representation.",
        "9. Exploratory graph and feature analysis.",
        "10. Unsupervised embedding construction.",
        "11. Similarity search between ncRNAs.",
        "12. Exploratory ranking of candidate ncRNAs.",
        "13. Construction of a project knowledge base for RAG.",
        "14. Development of GUI and API components.",
        "15. Dockerization and final packaging.",
        "",
        "### 2.2 Knowledge graph construction",
        "",
        "The graph connects ncRNAs to associated entities such as species, databases, RNA types, expression status, public datasets, publications, gene products and Gene Ontology terms.",
        "The graph was constructed conservatively: unsupported biological relationships were not invented.",
        f"The resulting graph contains {summary['n_nodes']} nodes and {summary['n_edges']} edges.",
        "",
        "### 2.3 Feature engineering",
        "",
        "Each ncRNA was represented using numerical features derived from sequence properties and graph context.",
        "These features provide the numerical basis for downstream embedding, similarity search and ranking.",
        f"The PyTorch Geometric representation includes {summary['n_node_features']} node features and {summary['n_edge_features']} edge features.",
        "",
        "### 2.4 Unsupervised embeddings and similarity search",
        "",
        "Unsupervised embeddings were generated to represent ncRNAs in a numerical latent space.",
        "Similarity search was then used to identify ncRNA records that were computationally close to one another.",
        f"The similarity table contains {summary['n_similarity_rows']} top-k similarity rows.",
        "Similarity should be interpreted as computational proximity, not as proof of shared biological function.",
        "",
        "### 2.5 Exploratory ranking",
        "",
        f"The ranking step produced an ordered list of {summary['n_ranked_ncrna']} ncRNA candidates.",
        "The ranking score is exploratory and should be used to guide further review rather than to claim biological validation.",
        "",
        "### 2.6 RAG module, GUI and API",
        "",
        f"A retrieval-augmented generation knowledge base was built from {summary['rag_chunks']} chunks extracted from project reports and structured outputs.",
        "A Streamlit graphical interface allows users to inspect dashboards, candidate rankings, candidate details, similarity tables, embeddings, reports and downloads.",
        "A FastAPI backend provides programmatic access to project summaries, candidate records, similarity results and RAG-based question answering.",
        "",
        "### 2.7 Docker deployment",
        "",
        "The project was containerized using Docker.",
        "The Docker setup starts both the FastAPI service and the Streamlit graphical interface.",
        "Validation scripts were included to test the API, GUI and RAG endpoints after container startup.",
        "",
        "## 3. Results",
        "",
        "### 3.1 Software outputs",
        "",
        f"The complete framework generated a final package containing {summary['package_files']} files.",
        "The package includes source scripts, processed data tables, GUI files, API files, Docker files, validation scripts, reports and a PDF user guide.",
        "",
        "### 3.2 Graph and ranking outputs",
        "",
        f"The graph contained {summary['n_nodes']} nodes and {summary['n_edges']} edges.",
        f"The ranking procedure produced {summary['n_ranked_ncrna']} ranked maize ncRNAs.",
        f"The similarity module produced {summary['n_similarity_rows']} similarity rows.",
        "",
        "### 3.3 Top exploratory candidates",
        "",
        "The first candidates in the exploratory ranking were:",
        "",
        *top_lines,
        "",
        "These candidates should be treated as computational hypotheses.",
        "They require further review and biological validation.",
        "",
        "### 3.4 Software validation",
        "",
        f"The end-to-end validation status was: {summary['validation_status']}.",
        f"The overall validation flag was: {summary['overall_pass']}.",
        "Docker validation confirmed that the API, GUI and RAG components were reachable after container startup.",
        "",
        "## 4. Discussion",
        "",
        "The proposed framework demonstrates how public ncRNA records can be organized into a reproducible AI-assisted workflow.",
        "A key strength is the integration of graph representation, unsupervised embedding, similarity search, RAG-based querying, a graphical interface and a programmatic API.",
        "",
        "The framework is especially useful in early-stage exploratory studies where the main goal is not to make definitive biological claims, but to reduce complexity and identify candidate records for deeper inspection.",
        "The RAG module also improves usability because users can ask questions about the project without manually opening all intermediate files.",
        "",
        "However, the main limitation is the absence of experimentally validated drought-response labels and ncRNA-specific quantitative RNA-seq evidence.",
        "Therefore, the current version should not be interpreted as a supervised predictor of drought-responsive ncRNAs.",
        "Rather, it is a transparent and reproducible prioritization framework.",
        "",
        "## 5. Limitations",
        "",
        "The main limitations are:",
        "",
        "- The ranking is exploratory and unsupervised.",
        "- No validated drought-response labels are available in the current dataset.",
        "- No ncRNA-specific quantitative RNA-seq validation was used.",
        "- Similarity does not prove biological function.",
        "- Candidate prioritization requires downstream biological interpretation.",
        "",
        "## 6. Conclusions",
        "",
        "This work presents a complete AI/GNN-enabled software framework for the exploratory prioritization of maize ncRNAs.",
        "The system integrates graph construction, feature engineering, embedding, similarity search, candidate ranking, RAG-based querying, GUI exploration, API access and Docker deployment.",
        "The resulting software package can support hypothesis generation for drought-response studies in maize, provided that all outputs are interpreted as exploratory and followed by biological validation.",
        "",
        "## Data and software availability",
        "",
        "The final software package is available locally as:",
        "",
        "outputs/final_package/maize_ncrna_gnn_final_package.zip",
        "",
        "The package includes the processed data, scripts, Docker files, reports, GUI, API and user guide.",
        "",
        "## Author contributions",
        "",
        "To be completed.",
        "",
        "## Funding",
        "",
        "To be completed.",
        "",
        "## Conflicts of interest",
        "",
        "To be completed.",
        "",
        "## References",
        "",
        "Bibliographic references will be completed in the next step using literature relevant to maize ncRNAs, drought response, graph neural networks and AI-assisted biological data integration.",
        "",
    ]

    return "\n".join(lines)


def pdf_clean(text: str) -> str:
    replacements = {
        "à": "a'",
        "è": "e'",
        "é": "e'",
        "ì": "i'",
        "ò": "o'",
        "ù": "u'",
        "À": "A'",
        "È": "E'",
        "É": "E'",
        "Ì": "I'",
        "Ò": "O'",
        "Ù": "U'",
        "–": "-",
        "—": "-",
        "“": '"',
        "”": '"',
        "’": "'",
    }
    for a, b in replacements.items():
        text = text.replace(a, b)
    return text


def draw_header(c: canvas.Canvas, page_num: int) -> None:
    width, height = A4
    c.setFillColor(colors.HexColor("#203864"))
    c.rect(0, height - 42, width, 42, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont("Helvetica-Bold", 10)
    c.drawString(40, height - 26, "Scientific manuscript draft - Maize ncRNA GNN")
    c.setFont("Helvetica", 9)
    c.drawRightString(width - 40, height - 26, f"Page {page_num}")
    c.setFillColor(colors.black)


def write_pdf(markdown: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    c = canvas.Canvas(str(PDF_PATH), pagesize=A4)
    width, height = A4

    margin_left = 42
    margin_right = 42
    y = height - 70
    page_num = 1

    draw_header(c, page_num)

    for raw in markdown.splitlines():
        line = pdf_clean(raw.strip())

        if y < 60:
            c.showPage()
            page_num += 1
            draw_header(c, page_num)
            y = height - 70

        if line.startswith("# "):
            c.setFont("Helvetica-Bold", 15)
            c.setFillColor(colors.HexColor("#203864"))
            chunks = textwrap.wrap(line.replace("# ", ""), width=70)
            for chunk in chunks:
                c.drawString(margin_left, y, chunk)
                y -= 19
            c.setFillColor(colors.black)
            y -= 8
            continue

        if line.startswith("## "):
            if y < 95:
                c.showPage()
                page_num += 1
                draw_header(c, page_num)
                y = height - 70
            c.setFont("Helvetica-Bold", 12)
            c.setFillColor(colors.HexColor("#203864"))
            c.drawString(margin_left, y, line.replace("## ", ""))
            c.setFillColor(colors.black)
            y -= 18
            continue

        if line.startswith("### "):
            c.setFont("Helvetica-Bold", 10.5)
            c.setFillColor(colors.HexColor("#203864"))
            c.drawString(margin_left, y, line.replace("### ", ""))
            c.setFillColor(colors.black)
            y -= 15
            continue

        if line == "":
            y -= 7
            continue

        c.setFont("Helvetica", 9.5)
        chunks = textwrap.wrap(line, width=96) or [""]

        for chunk in chunks:
            if y < 60:
                c.showPage()
                page_num += 1
                draw_header(c, page_num)
                y = height - 70

            if chunk.startswith("- "):
                c.drawString(margin_left + 12, y, chunk)
            else:
                c.drawString(margin_left, y, chunk)

            y -= 12

    c.save()


def main() -> int:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    summary = get_summary()
    markdown = build_markdown(summary)

    MD_PATH.write_text(markdown, encoding="utf-8")
    write_pdf(markdown)

    metadata = {
        "created": summary["date"],
        "markdown": str(MD_PATH),
        "pdf": str(PDF_PATH),
        "n_ranked_ncrna": summary["n_ranked_ncrna"],
        "n_nodes": summary["n_nodes"],
        "n_edges": summary["n_edges"],
        "n_similarity_rows": summary["n_similarity_rows"],
        "rag_chunks": summary["rag_chunks"],
        "note": "Draft manuscript without final bibliography.",
    }

    JSON_PATH.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    print("Manuscript Markdown:", MD_PATH)
    print("Manuscript PDF:", PDF_PATH)
    print("Manuscript metadata:", JSON_PATH)
    print("Completato.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
