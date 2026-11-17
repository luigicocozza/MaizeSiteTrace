#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()
OUT_DIR = BASE_DIR / "outputs/prjna419326"

DOWNLOAD_VALIDATED = OUT_DIR / f"{PROJECT}_mirnaseq_download_VALIDATO.txt"
QC_RAW = OUT_DIR / f"{PROJECT}_mirnaseq_light_fastq_qc_summary.txt"
TRIM_SUMMARY = OUT_DIR / f"{PROJECT}_mirnaseq_cutadapt_trimming_summary.txt"
QC_TRIMMED = OUT_DIR / f"{PROJECT}_mirnaseq_trimmed_qc_summary.txt"
COUNT_MATRIX_SUMMARY = OUT_DIR / f"{PROJECT}_mirnaseq_srna_count_matrix_summary.txt"
DROUGHT_SUMMARY = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_summary.txt"
MAPPING_SUMMARY = OUT_DIR / f"{PROJECT}_top500_srna_to_rnacentral_mapping_summary.txt"
RANKING_SUMMARY = OUT_DIR / f"{PROJECT}_experimental_srna_v3_ranking_summary.txt"
RANKING_CSV = OUT_DIR / f"{PROJECT}_experimental_srna_v3_ranking.csv"

REPORT_TXT = OUT_DIR / f"{PROJECT}_V3_report_metodologico_integrato.txt"
REPORT_MD = OUT_DIR / f"{PROJECT}_V3_report_metodologico_integrato.md"
REPORT_JSON = OUT_DIR / f"{PROJECT}_V3_report_metodologico_integrato.json"


def read_text(path: Path) -> str:
    if not path.exists():
        return f"[MANCANTE] {path}\n"
    return path.read_text(encoding="utf-8", errors="replace")


def read_top_ranking(path: Path, n: int = 20) -> list[dict]:
    if not path.exists():
        return []

    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)[:n]


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    top20 = read_top_ranking(RANKING_CSV, 20)

    sections = []

    sections.append("# Report metodologico integrato V3 - PRJNA419326")
    sections.append("")
    sections.append(f"Data generazione: {datetime.now().isoformat(timespec='seconds')}")
    sections.append("")
    sections.append("## 1. Obiettivo della Versione 3")
    sections.append("")
    sections.append(
        "La Versione 3 estende il framework AI/GNN sviluppato nelle versioni precedenti "
        "integrando un dataset sperimentale pubblico reale di small RNA/miRNA da mais "
        "sottoposto a siccità. L'obiettivo non è sostituire il livello RNAcentral/GNN, "
        "ma aggiungere un livello sperimentale di evidenza drought-responsive."
    )
    sections.append("")
    sections.append("Schema concettuale:")
    sections.append("")
    sections.append("```text")
    sections.append("RNAcentral maize ncRNA layer")
    sections.append("        +")
    sections.append("Graph / GAE AI layer")
    sections.append("        +")
    sections.append("PRJNA419326 experimental small RNA drought-response layer")
    sections.append("        =")
    sections.append("Prioritizzazione V3 integrata dei candidati")
    sections.append("```")
    sections.append("")

    sections.append("## 2. Download e validazione FASTQ")
    sections.append("")
    sections.append("```text")
    sections.append(read_text(DOWNLOAD_VALIDATED).strip())
    sections.append("```")
    sections.append("")

    sections.append("## 3. QC iniziale dei FASTQ grezzi")
    sections.append("")
    sections.append(
        "Il QC iniziale ha mostrato reads grezze lunghe 50 nt. Questo indica dati small RNA "
        "non ancora tagliati, con presenza di adattatore 3'."
    )
    sections.append("")
    sections.append("```text")
    sections.append(read_text(QC_RAW).strip())
    sections.append("```")
    sections.append("")

    sections.append("## 4. Trimming adattatore")
    sections.append("")
    sections.append(
        "L'adattatore 3' dominante è risultato compatibile con la sequenza Illumina small RNA:"
    )
    sections.append("")
    sections.append("```text")
    sections.append("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC")
    sections.append("```")
    sections.append("")
    sections.append("```text")
    sections.append(read_text(TRIM_SUMMARY).strip())
    sections.append("```")
    sections.append("")

    sections.append("## 5. QC post-trimming")
    sections.append("")
    sections.append(
        "Dopo il trimming, i campioni risultano nella finestra 18-30 nt, "
        "coerente con small RNA/miRNA/siRNA. La quota 20-24 nt misura la componente "
        "più direttamente compatibile con small RNA canonici."
    )
    sections.append("")
    sections.append("```text")
    sections.append(read_text(QC_TRIMMED).strip())
    sections.append("```")
    sections.append("")

    sections.append("## 6. Collasso reads e matrice di conteggio")
    sections.append("")
    sections.append(
        "Le reads trimmed identiche sono state collassate per campione. "
        "La matrice finale include sequenze con total_count >= 10 e presenti in almeno 2 campioni."
    )
    sections.append("")
    sections.append("```text")
    sections.append(read_text(COUNT_MATRIX_SUMMARY).strip())
    sections.append("```")
    sections.append("")

    sections.append("## 7. Screening drought-response")
    sections.append("")
    sections.append(
        "Lo screening preliminare usa RPM medi per gruppo e log2FC drought vs controllo. "
        "Non sostituisce un'analisi differenziale formale con FDR, ma serve per prioritizzare "
        "i candidati più promettenti nel framework V3."
    )
    sections.append("")
    sections.append("```text")
    sections.append(read_text(DROUGHT_SUMMARY).strip())
    sections.append("```")
    sections.append("")

    sections.append("## 8. Mapping contro RNAcentral")
    sections.append("")
    sections.append(
        "Il mapping esatto delle top 500 small RNA contro le 250 sequenze RNAcentral non ha "
        "individuato match esatti. Questo risultato non annulla il valore dei candidati: "
        "indica che le small RNA sperimentali PRJNA419326 devono essere trattate come nuovo "
        "livello sperimentale, distinto dal livello RNAcentral già presente."
    )
    sections.append("")
    sections.append("```text")
    sections.append(read_text(MAPPING_SUMMARY).strip())
    sections.append("```")
    sections.append("")

    sections.append("## 9. Ranking sperimentale V3")
    sections.append("")
    sections.append(
        "Il ranking V3 ordina le top 500 small RNA drought-responsive combinando classe biologica, "
        "score di risposta allo stress, abbondanza, presenza nei campioni e lunghezza."
    )
    sections.append("")
    sections.append("```text")
    sections.append(read_text(RANKING_SUMMARY).strip())
    sections.append("```")
    sections.append("")

    sections.append("## 10. Top 20 candidati V3")
    sections.append("")
    sections.append(
        "| rank | v3_rank_id | feature_id | class | length | total_count | B73_log2FC | H21_log2FC | score |"
    )
    sections.append(
        "|---:|---|---|---|---:|---:|---:|---:|---:|"
    )

    for i, row in enumerate(top20, start=1):
        sections.append(
            f"| {i} | {row['v3_rank_id']} | {row['feature_id']} | "
            f"{row['candidate_class']} | {row['length']} | {row['total_count']} | "
            f"{row['log2FC_B73_DS_vs_WW']} | {row['log2FC_H21_DS_vs_WW']} | "
            f"{row['v3_priority_score']} |"
        )

    sections.append("")
    sections.append("## 11. Conclusione tecnica")
    sections.append("")
    sections.append(
        "La Versione 3 introduce un dataset reale small RNA/miRNA da stress idrico nel mais. "
        "Il dataset è stato scaricato, validato tramite MD5, sottoposto a trimming, controllato "
        "per lunghezza post-trimming, collassato in sequenze uniche, trasformato in matrice di "
        "conteggio/RPM e usato per generare un ranking sperimentale di candidati drought-responsive."
    )
    sections.append("")
    sections.append(
        "Il risultato più importante è la creazione di un nuovo livello sperimentale di candidati "
        "small RNA, distinto dal livello RNAcentral, ma integrabile nel framework AI/GNN come "
        "evidenza biologica reale."
    )
    sections.append("")

    text = "\n".join(sections) + "\n"

    REPORT_MD.write_text(text, encoding="utf-8")
    REPORT_TXT.write_text(text.replace("# ", "").replace("## ", ""), encoding="utf-8")

    report = {
        "project": PROJECT,
        "created": datetime.now().isoformat(timespec="seconds"),
        "report_md": str(REPORT_MD),
        "report_txt": str(REPORT_TXT),
        "top20_count": len(top20),
        "inputs": {
            "download_validated": str(DOWNLOAD_VALIDATED),
            "qc_raw": str(QC_RAW),
            "trim_summary": str(TRIM_SUMMARY),
            "qc_trimmed": str(QC_TRIMMED),
            "count_matrix_summary": str(COUNT_MATRIX_SUMMARY),
            "drought_summary": str(DROUGHT_SUMMARY),
            "mapping_summary": str(MAPPING_SUMMARY),
            "ranking_summary": str(RANKING_SUMMARY),
            "ranking_csv": str(RANKING_CSV),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    print("REPORT METODOLOGICO V3 CREATO")
    print("=" * 80)
    print(REPORT_MD)
    print(REPORT_TXT)
    print(REPORT_JSON)
    print("")
    print("Anteprima top 5:")
    for i, row in enumerate(top20[:5], start=1):
        print(
            f"{i}. {row['v3_rank_id']} | {row['feature_id']} | "
            f"{row['candidate_class']} | len={row['length']} | "
            f"B73={row['log2FC_B73_DS_vs_WW']} | H21={row['log2FC_H21_DS_vs_WW']} | "
            f"score={row['v3_priority_score']}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
