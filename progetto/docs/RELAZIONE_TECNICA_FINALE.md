# Relazione tecnica finale

## Titolo del progetto

Framework AI/GNN per identificare e prioritizzare ncRNA coinvolti nella risposta alla siccita' nel mais

Data generazione: 2026-10-08 14:54

## 1. Scopo del software

Il software e' stato sviluppato per costruire una piattaforma completa di analisi computazionale dedicata agli ncRNA del mais.
L'obiettivo e' organizzare dati pubblici, costruire un grafo di conoscenza, calcolare feature, generare embedding, produrre similarita' tra ncRNA e ottenere un ranking esplorativo di candidati.

Il sistema non sostituisce la validazione biologica. Serve a proporre candidati da studiare meglio.

## 2. Componenti principali

- Dataset processati in formato CSV e FASTA.
- Grafo di conoscenza.
- Dataset PyTorch Geometric.
- Feature numeriche degli ncRNA.
- Embedding non supervisionati.
- Ricerca di similarita'.
- Ranking esplorativo.
- Chat RAG basata sui report del progetto.
- GUI Streamlit.
- API FastAPI.
- Docker per avvio semplificato.
- Guida PDF utente step-by-step.

## 3. Numeri principali del progetto

- ncRNA nel ranking: 250
- Righe similarita': 2500
- Nodi grafo: 384
- Archi grafo: 1134
- Feature di nodo: 47
- Feature di arco: 4
- Chunk knowledge base RAG: 342
- File nel pacchetto finale: 121

## 4. Pipeline realizzata

La pipeline e' composta dai seguenti passaggi:

1. Download e raccolta dati pubblici.
2. Ispezione dei file scaricati.
3. Profilazione delle fonti JSON.
4. Costruzione del dataset ncRNA del mais.
5. Costruzione del riepilogo di espressione.
6. Costruzione degli archi di conoscenza.
7. Costruzione delle label conservative.
8. Validazione degli output.
9. Costruzione delle tabelle nodi e archi.
10. Calcolo delle feature degli ncRNA.
11. Costruzione del dataset PyTorch Geometric.
12. Analisi esplorativa.
13. Calcolo embedding non supervisionati.
14. Calcolo similarita' tra ncRNA.
15. Ranking esplorativo.
16. Report metodologico.
17. Knowledge base RAG.
18. GUI Streamlit.
19. API FastAPI.
20. Validazione end-to-end.
21. Dockerizzazione.
22. Guida PDF utente.

## 5. Grafo di conoscenza

Il grafo rappresenta relazioni tra ncRNA, specie, database, tipo RNA, stato di espressione, dataset pubblici, pubblicazioni, termini GO e prodotti genici.
Il grafo contiene nodi e archi costruiti in modo conservativo, evitando di inventare relazioni biologiche non supportate.

## 6. Feature degli ncRNA

Le feature includono proprieta' di sequenza e informazioni derivate dalla struttura del grafo.
Queste feature sono usate per costruire rappresentazioni numeriche utili per embedding, similarita' e ranking.

## 7. Embedding e similarita'

Gli embedding trasformano ogni ncRNA in un vettore numerico.
La similarita' confronta questi vettori per trovare ncRNA computazionalmente vicini.
Una similarita' alta non dimostra funzione biologica comune: indica solo vicinanza nel modello computazionale.

## 8. Ranking esplorativo

Il ranking ordina gli ncRNA secondo un punteggio esplorativo.
Il ranking serve a scegliere candidati da controllare con letteratura, dati sperimentali e analisi biologiche.
Il ranking non dimostra coinvolgimento reale nella risposta alla siccita'.

## 9. Chat RAG

La Chat RAG permette di fare domande sul progetto usando i report e la knowledge base interna.
Esempi di domande utili:

- Quanti ncRNA sono stati analizzati?
- Quali sono i limiti metodologici?
- Quali sono i primi candidati?
- Quanti nodi e archi ha il grafo?

## 10. GUI Streamlit

La GUI consente di usare il software senza scrivere codice.
Le finestre principali sono:

- Dashboard.
- Ranking candidati.
- Dettaglio candidato.
- Ricerca avanzata.
- Similarita'.
- Embedding.
- Chat RAG.
- Report metodologico.
- Download.

## 11. API FastAPI

La API consente di interrogare il progetto da browser o da altri software.
Endpoint principali:

- /health
- /summary
- /candidates
- /candidate/{rnacentral_id}
- /candidate/{rnacentral_id}/similar
- /rag/ask
- /files/inventory

## 12. Docker

Docker consente di avviare il software in modo piu' semplice e riproducibile.
Sono stati creati:

- Dockerfile.
- .dockerignore.
- docker/start_services.sh.
- run_docker.sh.
- stop_docker.sh.
- validate_docker.sh.

Comandi principali:

cd ~/progetti/ncrna_rnacentral_ai
./run_docker.sh
./validate_docker.sh
./stop_docker.sh

## 13. Validazione

Stato validazione end-to-end: pass
Overall pass: True

La validazione controlla:

- presenza file richiesti;
- compilazione script;
- funzionamento API;
- funzionamento GUI;
- correttezza dei numeri principali;
- funzionamento RAG;
- disponibilita' dei candidati e delle similarita'.

## 14. Output principali

- outputs/reports/methodological_report.md
- outputs/reports/relazione_tecnica_finale.pdf
- outputs/manuals/guida_software_step_by_step.pdf
- outputs/final_package/maize_ncrna_gnn_final_package.zip
- data/processed/ranking/ncrna_exploratory_ranking.csv
- data/processed/similarity/ncrna_similarity_topk.csv
- data/processed/ncrna_features.csv
- data/processed/embeddings/ncrna_embeddings.csv

## 15. Limiti metodologici

Il progetto e' esplorativo.
Non sono presenti label biologiche validate per dimostrare direttamente la risposta alla siccita'.
Non sono presenti dati RNA-seq ncRNA-specifici usati per validare quantitativamente i candidati.
Il ranking non deve essere interpretato come prova biologica.

## 16. Uso corretto dei risultati

I risultati devono essere usati per:

- selezionare candidati da studiare;
- confrontare ncRNA simili;
- generare ipotesi biologiche;
- preparare analisi successive;
- supportare una pubblicazione metodologica.

Non devono essere usati per affermare direttamente che un ncRNA sia validato nella risposta alla siccita'.

## 17. Conclusione

Il progetto ha prodotto una piattaforma software completa e funzionante.
Il sistema include dati processati, grafo, ranking, similarita', RAG, GUI, API, Docker, guida utente e pacchetto finale.

Il software e' ora pronto per essere descritto nella relazione tecnica e poi trasformato in articolo scientifico.
