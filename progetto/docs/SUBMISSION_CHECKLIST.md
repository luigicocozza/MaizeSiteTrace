# Submission checklist

Data generazione: 2026-10-08 15:07

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
