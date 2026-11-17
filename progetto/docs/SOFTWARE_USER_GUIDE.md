# Guida semplice del software Maize ncRNA GNN

Data generazione: 2026-10-08 14:21

Questa guida spiega come avviare e usare il software in modo semplice.

Il software serve a esplorare ncRNA del mais con metodi AI, grafo, embedding, ranking e chat RAG.

ATTENZIONE:
Il software NON dimostra che un ncRNA sia davvero coinvolto nella risposta alla siccita'.
Il ranking e' esplorativo. Serve per scegliere candidati da controllare meglio con letteratura, RNA-seq e validazione biologica.

================================================================================
1. COSA CONTIENE IL SOFTWARE
================================================================================

Il software contiene:

- una GUI, cioe' una finestra web facile da usare;
- una API FastAPI, cioe' un servizio interrogabile anche da altri programmi;
- una chat RAG, cioe' una chat che risponde usando i file e i report del progetto;
- un ranking di 250 ncRNA;
- una tabella di similarita' tra ncRNA;
- feature di sequenza;
- embedding;
- grafo con 384 nodi e 1134 archi;
- report e file scaricabili.

================================================================================
2. AVVIO CON DOCKER
================================================================================

Docker permette di avviare tutto senza configurare ogni volta Python.

Comando da usare:

cd ~/progetti/ncrna_rnacentral_ai
./run_docker.sh

Quando Docker parte, aprire nel browser:

GUI:
http://localhost:8501

API:
http://localhost:8000/docs

Per controllare che funzioni:

./validate_docker.sh

Per fermare Docker:

./stop_docker.sh

================================================================================
3. AVVIO SENZA DOCKER
================================================================================

Usare questi comandi:

cd ~/progetti/ncrna_rnacentral_ai
source .venv/bin/activate
./run_project.sh

Poi aprire:

GUI:
http://localhost:8501

API:
http://localhost:8000/docs

Per fermare:

./stop_project.sh

================================================================================
4. FINESTRA DASHBOARD
================================================================================

La Dashboard e' la prima pagina del software.

Mostra:

- numero di ncRNA nel ranking;
- numero di feature;
- numero di nodi del grafo;
- numero di archi del grafo;
- grafici principali;
- stato dei file.

INPUT:
Non richiede input. Serve solo a vedere lo stato generale.

OUTPUT:
Mostra numeri, grafici e tabelle di controllo.

COME INTERPRETARE:
Se vedi 250 ncRNA, 384 nodi e 1134 archi, il progetto e' caricato correttamente.

================================================================================
5. FINESTRA RANKING CANDIDATI
================================================================================

Questa finestra mostra la lista dei candidati ncRNA ordinati per punteggio esplorativo.

INPUT DISPONIBILI:

- tipo RNA;
- score minimo;
- numero di righe da mostrare.

OUTPUT:

- tabella dei candidati filtrati;
- download del CSV filtrato.

SIGNIFICATO:
Il primo candidato ha rank 1. Un rank alto non significa prova biologica.
Significa solo che il candidato e' interessante per revisione successiva.

================================================================================
6. FINESTRA DETTAGLIO CANDIDATO
================================================================================

Questa finestra permette di guardare un singolo ncRNA.

INPUT:

- RNAcentral ID del candidato.

SCHEDE DISPONIBILI:

1. Sintesi avanzata
   Spiega il candidato in parole semplici.

2. Ranking
   Mostra la riga completa del ranking.

3. Top-k simili
   Mostra gli ncRNA piu' simili.

4. Feature
   Mostra lunghezza, GC, composizione e informazioni del grafo.

5. Embedding
   Mostra le coordinate numeriche usate per similarita' e ranking.

6. Export
   Permette di scaricare:
   - CSV del candidato;
   - report Markdown;
   - pacchetto JSON.

OUTPUT:
Report del candidato, tabelle e file scaricabili.

ATTENZIONE:
Anche se un candidato ha score alto, non e' automaticamente drought-responsive.

================================================================================
7. FINESTRA RICERCA AVANZATA
================================================================================

Questa pagina serve per cercare candidati con filtri piu' precisi.

INPUT:

- testo da cercare;
- tipo RNA;
- intervallo di rank;
- intervallo di score;
- lunghezza sequenza;
- GC fraction;
- evidenza quantitativa;
- evidenza drought;
- label unknown.

OUTPUT:

- tabella filtrata;
- numero di risultati;
- distribuzione RNA type;
- download CSV.

ESEMPIO:
Per cercare il primo candidato, scrivere:

URS0003C75441

================================================================================
8. FINESTRA SIMILARITA'
================================================================================

Questa finestra mostra quali ncRNA sono simili tra loro.

INPUT:
Non e' obbligatorio. Si guarda la tabella.

OUTPUT:

- query ncRNA;
- ncRNA simile;
- tipo RNA;
- cosine similarity;
- differenza di lunghezza;
- differenza GC.

COME INTERPRETARE:
Similarita' alta significa che due ncRNA sono simili nei dati computazionali.
Non significa per forza che abbiano la stessa funzione biologica.

================================================================================
9. FINESTRA EMBEDDING
================================================================================

Questa pagina mostra gli embedding degli ncRNA.

INPUT:
Nessuno.

OUTPUT:

- tabella embedding;
- coordinate emb_001, emb_002, ecc.;
- grafico 2D.

SIGNIFICATO:
Gli embedding sono numeri che riassumono le caratteristiche degli ncRNA.
Servono per similarita', clustering e ranking.

================================================================================
10. FINESTRA CHAT RAG
================================================================================

La Chat RAG permette di fare domande sul progetto.

ESEMPI DI DOMANDE:

- Quanti ncRNA sono stati analizzati?
- Quali sono i primi 5 candidati?
- Il ranking dimostra risposta alla siccita'?
- Quanti nodi e archi ha il grafo?
- Quali sono i limiti metodologici?

INPUT:
Una domanda scritta dall'utente.

OUTPUT:
Una risposta basata sui report e sui file del progetto.

IMPORTANTE:
La chat deve rispondere in modo prudente.
Non deve dire che un candidato e' validato se non ci sono prove biologiche.

================================================================================
11. FINESTRA REPORT METODOLOGICO
================================================================================

Questa finestra mostra il report metodologico automatico.

INPUT:
Nessuno.

OUTPUT:

- summary del progetto;
- spiegazione pipeline;
- limiti;
- candidati;
- note metodologiche.

USO:
Serve per capire come e' stato costruito il progetto.

================================================================================
12. FINESTRA DOWNLOAD
================================================================================

Questa finestra serve per scaricare i risultati.

INPUT:
Click sui pulsanti.

OUTPUT DISPONIBILI:

- report export Markdown;
- inventario CSV;
- ZIP progetto;
- ranking;
- similarita';
- feature;
- embedding;
- nodi grafo;
- archi grafo.

USO:
Serve per salvare i risultati o portarli su un altro computer.

================================================================================
13. API FASTAPI
================================================================================

La API si apre qui:

http://localhost:8000/docs

Endpoint principali:

/health
Controlla se la API funziona.

/summary
Mostra numeri principali del progetto.

/candidates
Mostra la lista dei candidati.

/candidate/URS0003C75441
Mostra il dettaglio di un candidato.

/candidate/URS0003C75441/similar
Mostra gli ncRNA simili.

/rag/ask
Permette di fare domande alla chat RAG.

================================================================================
14. FILE PRINCIPALI DI INPUT
================================================================================

I file principali usati dal software sono:

data/processed/ranking/ncrna_exploratory_ranking.csv
Ranking dei candidati.

data/processed/similarity/ncrna_similarity_topk.csv
Similarita' tra ncRNA.

data/processed/ncrna_features.csv
Feature degli ncRNA.

data/processed/embeddings/ncrna_embeddings.csv
Embedding degli ncRNA.

data/processed/graph_nodes.csv
Nodi del grafo.

data/processed/graph_edges.csv
Archi del grafo.

data/processed/rag/project_knowledge_base.jsonl
Knowledge base usata dalla chat.

================================================================================
15. FILE PRINCIPALI DI OUTPUT
================================================================================

Output principali:

outputs/reports/methodological_report.md
Report metodologico.

outputs/e2e/e2e_validation_summary.txt
Validazione end-to-end.

outputs/final_package/maize_ncrna_gnn_final_package.zip
Pacchetto finale.

outputs/manuals/guida_software_step_by_step.pdf
Questa guida PDF.

================================================================================
16. COME FERMARE IL SOFTWARE
================================================================================

Se usi Docker:

./stop_docker.sh

Se usi avvio locale:

./stop_project.sh

================================================================================
17. COSA NON DEVE ESSERE DETTO
================================================================================

Non dire:

- questo ncRNA e' sicuramente coinvolto nella siccita';
- il ranking dimostra una risposta biologica;
- la similarita' dimostra funzione comune.

Dire invece:

- questo ncRNA e' un candidato esplorativo;
- va controllato con letteratura;
- servono dati RNA-seq;
- serve validazione biologica.

================================================================================
18. CONCLUSIONE
================================================================================

Il software e' utile per esplorare, ordinare e interrogare candidati ncRNA del mais.

Il risultato principale e' una piattaforma funzionante con:

- Docker;
- GUI;
- API;
- Chat RAG;
- ranking;
- similarita';
- export;
- validazione end-to-end.

L'articolo scientifico deve essere scritto solo dopo questa fase software e dopo la guida.
