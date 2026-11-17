# Maize ncRNA GNN

Framework AI/GNN per identificare e prioritizzare ncRNA coinvolti nella risposta alla siccità nel mais.

## Stato progetto

- Validazione end-to-end: **pass**
- Overall pass: **True**
- ncRNA nel ranking: **250**
- Righe similarità: **2500**
- Nodi grafo: **384**
- Archi grafo: **1134**
- Chunk RAG: **342**

## Avvio rapido

cd ~/progetti/ncrna_rnacentral_ai
source .venv/bin/activate
./run_project.sh

GUI: http://localhost:8501
API: http://localhost:8000/docs

## Stop

./stop_project.sh

## Nota metodologica

Ranking, similarità ed embedding sono esplorativi e non supervisionati.
Non dimostrano validazione biologica della risposta alla siccità.
