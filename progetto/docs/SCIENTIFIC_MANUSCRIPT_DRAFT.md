# An AI/GNN and RAG-enabled framework for exploratory prioritization of maize ncRNAs in drought-response studies

## Manuscript type

Draft research article / methodological software paper.

## Authors

To be completed.

## Affiliations

To be completed.

## Abstract

Non-coding RNAs are increasingly relevant in plant regulatory biology, but their systematic prioritization in crop stress-response studies remains challenging when curated labels and ncRNA-specific quantitative expression data are limited.
Here, we present a reproducible AI/GNN-enabled software framework for the exploratory organization, representation and prioritization of public maize ncRNA records in the context of drought-response hypothesis generation.
The framework integrates public-data-derived ncRNA records into a knowledge graph with 384 nodes and 1134 edges, computes sequence- and graph-derived features, builds a PyTorch Geometric-compatible representation, generates unsupervised embeddings, performs similarity search, and produces an exploratory ranking of 250 maize ncRNAs.
The system also includes a retrieval-augmented generation module built from 342 knowledge chunks, a Streamlit graphical interface, a FastAPI backend, Docker deployment scripts and a step-by-step user guide.
The resulting platform is intended as a decision-support and hypothesis-generation tool, not as biological validation.
All candidate rankings should therefore be interpreted cautiously and followed by literature review, ncRNA-specific expression analysis and experimental validation.

## Keywords

maize; Zea mays; non-coding RNA; drought; graph neural network; knowledge graph; RAG; Streamlit; FastAPI; exploratory ranking

## 1. Introduction

Drought is one of the major constraints affecting crop productivity.
In maize, molecular responses to drought involve complex regulatory processes, including transcriptional, post-transcriptional and epigenetic mechanisms.
Non-coding RNAs may contribute to these regulatory processes, but their prioritization is difficult when available information is fragmented across public databases, publications and functional annotations.

Traditional candidate selection often depends on manual literature review or on expression datasets focused on protein-coding genes.
For ncRNA-centered studies, this creates a practical gap: many public records exist, but they are not always organized into a computational structure that supports systematic exploration, similarity search and transparent candidate prioritization.

The objective of this work is to describe a software framework that integrates public maize ncRNA information into a graph-based and AI-assisted workflow.
The framework is designed to support hypothesis generation for downstream drought-response studies.

## 2. Materials and methods

### 2.1 Overall workflow

The pipeline was organized into the following phases:

1. Collection of public data sources.
2. Construction of a maize ncRNA sequence table.
3. Construction of conservative expression and label tables.
4. Construction of known knowledge-graph edges.
5. Validation of intermediate outputs.
6. Construction of graph node and edge tables.
7. Feature engineering for ncRNA records.
8. Conversion to a PyTorch Geometric-compatible graph representation.
9. Exploratory graph and feature analysis.
10. Unsupervised embedding construction.
11. Similarity search between ncRNAs.
12. Exploratory ranking of candidate ncRNAs.
13. Construction of a project knowledge base for RAG.
14. Development of GUI and API components.
15. Dockerization and final packaging.

### 2.2 Knowledge graph construction

The graph connects ncRNAs to associated entities such as species, databases, RNA types, expression status, public datasets, publications, gene products and Gene Ontology terms.
The graph was constructed conservatively: unsupported biological relationships were not invented.
The resulting graph contains 384 nodes and 1134 edges.

### 2.3 Feature engineering

Each ncRNA was represented using numerical features derived from sequence properties and graph context.
These features provide the numerical basis for downstream embedding, similarity search and ranking.
The PyTorch Geometric representation includes 47 node features and 4 edge features.

### 2.4 Unsupervised embeddings and similarity search

Unsupervised embeddings were generated to represent ncRNAs in a numerical latent space.
Similarity search was then used to identify ncRNA records that were computationally close to one another.
The similarity table contains 2500 top-k similarity rows.
Similarity should be interpreted as computational proximity, not as proof of shared biological function.

### 2.5 Exploratory ranking

The ranking step produced an ordered list of 250 ncRNA candidates.
The ranking score is exploratory and should be used to guide further review rather than to claim biological validation.

### 2.6 RAG module, GUI and API

A retrieval-augmented generation knowledge base was built from 342 chunks extracted from project reports and structured outputs.
A Streamlit graphical interface allows users to inspect dashboards, candidate rankings, candidate details, similarity tables, embeddings, reports and downloads.
A FastAPI backend provides programmatic access to project summaries, candidate records, similarity results and RAG-based question answering.

### 2.7 Docker deployment

The project was containerized using Docker.
The Docker setup starts both the FastAPI service and the Streamlit graphical interface.
Validation scripts were included to test the API, GUI and RAG endpoints after container startup.

## 3. Results

### 3.1 Software outputs

The complete framework generated a final package containing 124 files.
The package includes source scripts, processed data tables, GUI files, API files, Docker files, validation scripts, reports and a PDF user guide.

### 3.2 Graph and ranking outputs

The graph contained 384 nodes and 1134 edges.
The ranking procedure produced 250 ranked maize ncRNAs.
The similarity module produced 2500 similarity rows.

### 3.3 Top exploratory candidates

The first candidates in the exploratory ranking were:

- Rank 1: URS0003C75441; RNA type: pre_miRNA; exploratory score: NA
- Rank 2: URS0003C7543B; RNA type: pre_miRNA; exploratory score: NA
- Rank 3: URS0003C753CA; RNA type: pre_miRNA; exploratory score: NA
- Rank 4: URS0003C753C9; RNA type: pre_miRNA; exploratory score: NA
- Rank 5: URS0003C753AE; RNA type: pre_miRNA; exploratory score: NA
- Rank 6: URS0003C7548D; RNA type: pre_miRNA; exploratory score: NA
- Rank 7: URS0003C7545D; RNA type: pre_miRNA; exploratory score: NA
- Rank 8: URS0003C75489; RNA type: pre_miRNA; exploratory score: NA
- Rank 9: URS0003C75406; RNA type: pre_miRNA; exploratory score: NA
- Rank 10: URS0003C7548A; RNA type: pre_miRNA; exploratory score: NA

These candidates should be treated as computational hypotheses.
They require further review and biological validation.

### 3.4 Software validation

The end-to-end validation status was: pass.
The overall validation flag was: True.
Docker validation confirmed that the API, GUI and RAG components were reachable after container startup.

## 4. Discussion

The proposed framework demonstrates how public ncRNA records can be organized into a reproducible AI-assisted workflow.
A key strength is the integration of graph representation, unsupervised embedding, similarity search, RAG-based querying, a graphical interface and a programmatic API.

The framework is especially useful in early-stage exploratory studies where the main goal is not to make definitive biological claims, but to reduce complexity and identify candidate records for deeper inspection.
The RAG module also improves usability because users can ask questions about the project without manually opening all intermediate files.

However, the main limitation is the absence of experimentally validated drought-response labels and ncRNA-specific quantitative RNA-seq evidence.
Therefore, the current version should not be interpreted as a supervised predictor of drought-responsive ncRNAs.
Rather, it is a transparent and reproducible prioritization framework.

## 5. Limitations

The main limitations are:

- The ranking is exploratory and unsupervised.
- No validated drought-response labels are available in the current dataset.
- No ncRNA-specific quantitative RNA-seq validation was used.
- Similarity does not prove biological function.
- Candidate prioritization requires downstream biological interpretation.

## 6. Conclusions

This work presents a complete AI/GNN-enabled software framework for the exploratory prioritization of maize ncRNAs.
The system integrates graph construction, feature engineering, embedding, similarity search, candidate ranking, RAG-based querying, GUI exploration, API access and Docker deployment.
The resulting software package can support hypothesis generation for drought-response studies in maize, provided that all outputs are interpreted as exploratory and followed by biological validation.

## Data and software availability

The final software package is available locally as:

outputs/final_package/maize_ncrna_gnn_final_package.zip

The package includes the processed data, scripts, Docker files, reports, GUI, API and user guide.

## Author contributions

To be completed.

## Funding

To be completed.

## Conflicts of interest

To be completed.

## References

Bibliographic references will be completed in the next step using literature relevant to maize ncRNAs, drought response, graph neural networks and AI-assisted biological data integration.
