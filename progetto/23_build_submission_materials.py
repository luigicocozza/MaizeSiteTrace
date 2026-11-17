#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

BASE_DIR = Path(".").resolve()
DOCS_DIR = BASE_DIR / "docs"
OUT_DIR = BASE_DIR / "outputs/manuscript"

CHECKLIST_MD = DOCS_DIR / "SUBMISSION_CHECKLIST.md"
COVER_LETTER_MD = DOCS_DIR / "COVER_LETTER_DRAFT.md"
HIGHLIGHTS_MD = DOCS_DIR / "HIGHLIGHTS.md"
SUBMISSION_JSON = OUT_DIR / "submission_materials_metadata.json"


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.strip() + "\n", encoding="utf-8")


def main() -> int:
    DOCS_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    today = datetime.now().strftime("%Y-%m-%d %H:%M")

    checklist = f"""
# Submission checklist

Data generazione: {today}

## File articolo

- [x] Bozza articolo Markdown: docs/SCIENTIFIC_MANUSCRIPT_DRAFT.md
- [x] Bozza articolo PDF: outputs/manuscript/scientific_manuscript_draft.pdf
- [x] Articolo con bibliografia Markdown: docs/SCIENTIFIC_MANUSCRIPT_WITH_REFERENCES.md
- [x] Articolo con bibliografia PDF: outputs/manuscript/scientific_manuscript_with_references.pdf
- [x] Bibliografia JSON: outputs/manuscript/references_used.json

## File software

- [x] Dockerfile
- [x] .dockerignore
- [x] docker/start_services.sh
- [x] run_docker.sh
- [x] stop_docker.sh
- [x] validate_docker.sh

## File documentazione

- [x] Guida software PDF: outputs/manuals/guida_software_step_by_step.pdf
- [x] Relazione tecnica PDF: outputs/reports/relazione_tecnica_finale.pdf
- [x] README.md
- [x] docs/README_OPERATIVO.md
- [x] docs/API_REFERENCE.md
- [x] docs/GUI_GUIDE.md

## File dati principali

- [x] data/processed/ranking/ncrna_exploratory_ranking.csv
- [x] data/processed/similarity/ncrna_similarity_topk.csv
- [x] data/processed/ncrna_features.csv
- [x] data/processed/embeddings/ncrna_embeddings.csv
- [x] data/processed/graph_nodes.csv
- [x] data/processed/graph_edges.csv
- [x] data/processed/rag/project_knowledge_base.jsonl

## Campi da completare prima dell'invio

- [ ] Nome autore/i.
- [ ] Affiliazione/i.
- [ ] Email corresponding author.
- [ ] ORCID, se disponibile.
- [ ] Funding.
- [ ] Conflicts of interest.
- [ ] Data availability statement definitivo.
- [ ] Software availability statement definitivo.
- [ ] Ethical statement, se richiesto.
- [ ] Ringraziamenti.
- [ ] Scelta rivista definitiva.
- [ ] Adattamento formato alla rivista scelta.

## Avvertenza metodologica da mantenere

Il ranking e' esplorativo e non supervisionato.
Non dimostra validazione biologica della risposta alla siccita'.
I candidati richiedono revisione bibliografica, dati RNA-seq e validazione sperimentale.
"""

    highlights = """
# Highlights

- A reproducible AI/GNN-enabled framework was developed for exploratory maize ncRNA prioritization.
- The workflow integrates public ncRNA records, knowledge graph construction, feature engineering and embeddings.
- The software provides ncRNA similarity search and an exploratory ranking of 250 maize ncRNAs.
- A RAG module, Streamlit GUI, FastAPI backend and Docker deployment improve usability and reproducibility.
- Candidate rankings are hypothesis-generating and require downstream biological validation.
"""

    cover_letter = """
# Cover letter draft

Dear Editor,

We are pleased to submit the manuscript entitled:

"An AI/GNN and RAG-enabled framework for exploratory prioritization of maize ncRNAs in drought-response studies"

for consideration as a methodological/software-oriented research article.

The manuscript describes a reproducible software framework designed to organize public maize non-coding RNA records into a graph-based AI workflow. The system integrates knowledge graph construction, feature engineering, PyTorch Geometric-compatible graph representation, unsupervised embeddings, similarity search, exploratory candidate ranking, retrieval-augmented generation, a Streamlit graphical interface, a FastAPI backend and Docker deployment.

The work is intended as a transparent hypothesis-generation platform for maize drought-response studies. Importantly, the manuscript clearly states that the current ranking is exploratory and does not represent biological validation. Candidate ncRNAs must be further evaluated by literature review, ncRNA-specific expression analysis and experimental validation.

The submitted package includes the manuscript draft, software documentation, Docker files, processed data tables, validation outputs and a step-by-step user guide.

We believe this work may be of interest to readers working on crop bioinformatics, AI-assisted biological data integration, non-coding RNA prioritization and reproducible computational pipelines.

Sincerely,

[Author names to be completed]
[Affiliation to be completed]
[Corresponding author email to be completed]
"""

    metadata = {
        "created": today,
        "checklist": str(CHECKLIST_MD),
        "highlights": str(HIGHLIGHTS_MD),
        "cover_letter": str(COVER_LETTER_MD),
        "status": "submission_materials_draft_created",
        "fields_to_complete": [
            "authors",
            "affiliations",
            "corresponding_author_email",
            "orcid",
            "funding",
            "conflicts_of_interest",
            "journal_format",
        ],
    }

    write(CHECKLIST_MD, checklist)
    write(HIGHLIGHTS_MD, highlights)
    write(COVER_LETTER_MD, cover_letter)
    write(SUBMISSION_JSON, json.dumps(metadata, indent=2, ensure_ascii=False))

    print("Checklist:", CHECKLIST_MD)
    print("Highlights:", HIGHLIGHTS_MD)
    print("Cover letter:", COVER_LETTER_MD)
    print("Metadata:", SUBMISSION_JSON)
    print("Completato.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
