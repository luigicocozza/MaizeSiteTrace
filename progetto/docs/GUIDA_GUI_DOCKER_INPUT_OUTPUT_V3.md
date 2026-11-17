# Guida GUI, Docker, input e output - Progetto V3 PRJNA419326

## 1. Realizzazione progetto

Ing. Luigi Cocozza

Esperto in: AI Machine Learning Engineer

Ingegnere di Machine Learning IA

## 2. Scopo della GUI

La GUI Streamlit permette di consultare e usare i risultati del progetto V3.

Funzioni principali:

- consultazione rete neurale GAE;
- consultazione grafo AI/GNN;
- caricamento o lettura dati del ricercatore;
- ranking dei candidati small RNA/ncRNA;
- consultazione predizioni ML full e anti-leakage;
- interrogazione dei risultati tramite chat;
- accesso ai report finali.

## 3. Avvio locale senza Docker

Entrare nella cartella del progetto:

    cd ~/progetti/ncrna_rnacentral_ai_v3_prjna419326
    source ~/progetti/ncrna_rnacentral_ai/.venv/bin/activate

Avviare la GUI:

    streamlit run app/streamlit_app.py

Aprire nel browser:

    http://localhost:8501

Se si accede dalla macchina host alla VM VirtualBox, usare l'indirizzo IP della VM oppure configurare il port forwarding della porta 8501.

## 4. Avvio con Docker

Verificare Docker:

    docker --version
    docker compose version

Avvio con Docker Compose:

    docker compose up --build

Avvio in background:

    docker compose up -d --build

Controllare lo stato:

    docker compose ps

Leggere i log:

    docker compose logs -f

Fermare il servizio:

    docker compose down

Se si usa solo Dockerfile:

    docker build -t maize-ncrna-gnn-v3 .
    docker run --rm -it -p 8501:8501 maize-ncrna-gnn-v3


## 5. Struttura della GUI

File principali:

    app/streamlit_app.py
    app/pages/01_Rete_neurale_GAE.py
    app/pages/02_Dati_ricercatore.py
    app/pages/03_Ranking_dati_ricercatore.py
    app/pages/04_Chat_dati_ricercatore.py

La GUI contiene la dicitura:

    Realizzazione progetto
    Ing. Luigi Cocozza
    Esperto in: AI Machine Learning Engineer
    Ingegnere di Machine Learning IA

## 6. Pagina 01 - Rete neurale GAE

Questa pagina serve per consultare:

- grafo AI/GNN;
- Graph Autoencoder;
- embedding dei nodi;
- metriche AUC e Average Precision;
- risultati del modello GAE V3.

Output collegati:

    outputs/v3_gae/
    data/processed/embeddings/
    models/graph_autoencoder_v3/

## 7. Pagina 02 - Dati ricercatore

Questa pagina serve per caricare o consultare dati sperimentali.

Input consigliati:

    CSV
    TSV
    TXT tabellare
    FASTA, se previsto dalla pagina

Colonne minime consigliate:

    id
    sequence

Colonne utili:

    length
    total_count
    samples_detected
    log2FC_B73_DS_vs_WW
    log2FC_H21_DS_vs_WW
    candidate_class
    stress_score
    priority_score
    annotation

## 8. Pagina 03 - Ranking dati ricercatore

Questa pagina serve per consultare:

- ranking dei candidati;
- score di priorità;
- classe del candidato;
- probabilità del modello anti-leakage;
- tier finale;
- target support.

Output principali:

    outputs/prjna419326/final_v3_report/PRJNA419326_V3_FINAL_srna_candidate_priority_table.tsv
    outputs/prjna419326/final_v3_report/PRJNA419326_V3_FINAL_top100_srna_candidates.tsv

## 9. Pagina 04 - Chat dati ricercatore

Questa pagina serve per interrogare i risultati.

Esempi di domande:

    Quali sono i primi 20 small RNA candidati?
    Quali candidati sono TIER_1?
    Che differenza c'è tra modello full e modello anti-leakage?
    Quali small RNA hanno supporto opposto small RNA/mRNA?
    Quali file contengono le predizioni ML?
    Quali geni risultano drought-responsive?


## 10. Input principali già presenti

File principali consultabili dalla GUI o collegati ai risultati:

    outputs/prjna419326/final_v3_report/PRJNA419326_V3_FINAL_srna_candidate_priority_table.tsv
    outputs/prjna419326/final_v3_report/PRJNA419326_V3_FINAL_top100_srna_candidates.tsv
    outputs/prjna419326/final_v3_report/PRJNA419326_V3_FINAL_scientific_report.md
    outputs/prjna419326/PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv
    outputs/prjna419326/PRJNA419326_srna_supervised_ML_predictions.tsv
    outputs/prjna419326/PRJNA419326_srna_rnaseq_predicted_target_links.tsv
    outputs/prjna419326/PRJNA419326_srna_supervised_label_candidates.tsv
    data/processed/rnaseq_matrices/PRJNA419326_rnaseq_gene_TPM_matrix.tsv
    data/processed/rnaseq_matrices/PRJNA419326_rnaseq_gene_NumReads_matrix.tsv
    data/processed/embeddings/gae_v3_experimental_srna_embeddings.csv
    models/srna_supervised_ml_v1_antileakage/PRJNA419326_srna_supervised_antileakage_best_model.joblib

## 11. Output principali

Tabelle finali:

    outputs/prjna419326/final_v3_report/PRJNA419326_V3_FINAL_srna_candidate_priority_table.tsv
    outputs/prjna419326/final_v3_report/PRJNA419326_V3_FINAL_top100_srna_candidates.tsv

Report finale:

    outputs/prjna419326/final_v3_report/PRJNA419326_V3_FINAL_scientific_report.md
    outputs/prjna419326/final_v3_report/PRJNA419326_V3_FINAL_scientific_report.txt

Output machine learning:

    outputs/prjna419326/PRJNA419326_srna_supervised_ML_model_comparison.tsv
    outputs/prjna419326/PRJNA419326_srna_supervised_ML_predictions.tsv
    outputs/prjna419326/PRJNA419326_srna_supervised_ML_antileakage_model_comparison.tsv
    outputs/prjna419326/PRJNA419326_srna_supervised_ML_antileakage_predictions.tsv
    outputs/prjna419326/PRJNA419326_srna_supervised_ML_antileakage_features_used.tsv

Modelli:

    models/graph_autoencoder_v3/
    models/srna_supervised_ml_v1/
    models/srna_supervised_ml_v1_antileakage/

## 12. Interpretazione dei tier

TIER_1_high_confidence_supported:
candidati con alta probabilità anti-leakage e supporto biologico computazionale.

TIER_2_model_prioritized_unlabeled:
candidati non etichettati positivi dalle regole, ma suggeriti dal modello anti-leakage.

TIER_3_supported_moderate_model:
candidati positivi con supporto biologico e probabilità moderata.

TIER_4_unlabeled_moderate_model:
candidati non etichettati ma con probabilità moderata.

TIER_5_rule_supported_full_model_only:
candidati supportati dal modello full, ma meno forti nel controllo anti-leakage.

TIER_6_low_priority_or_negative:
candidati a priorità bassa o negativi silver.

## 13. Nota metodologica

Le label supervisionate sono computazionali/silver labels.

Non rappresentano validazione sperimentale.

Il modello più difendibile è quello anti-leakage, perché esclude le feature usate direttamente nella costruzione delle label.

## 14. File esclusi dal pacchetto ZIP

Per ridurre la dimensione del pacchetto vengono esclusi:

    FASTQ grezzi RNA-Seq
    FASTQ grezzi small RNA-Seq
    reference Ensembl Plants
    indice Salmon
    cartelle Salmon quant per singolo campione
    file temporanei .aria2
    backup temporanei .bak_before_author_v3

## 15. Controlli rapidi

Verificare GUI:

    grep -R "Ing. Luigi Cocozza\|AI Machine Learning Engineer\|Ingegnere di Machine Learning IA" \
      --exclude-dir=.venv \
      --exclude-dir=delivery \
      --exclude-dir=data \
      --exclude-dir=models \
      --exclude-dir=outputs \
      -n app | head -50

Verificare report finali:

    ls -lh outputs/prjna419326/final_v3_report/

Verificare modelli:

    find models -maxdepth 3 -type f | sort

## 16. Procedura finale di consegna

1. Aggiornare la GUI.
2. Generare il report finale V3.
3. Generare questa guida.
4. Creare lo ZIP finale.
5. Conservare SHA256 dello ZIP.
6. Consegnare ZIP + SHA256.

Comando finale:

    python 63_package_v3_final_delivery.py

