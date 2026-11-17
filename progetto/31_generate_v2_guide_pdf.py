#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    PageBreak,
    Table,
    TableStyle,
)


BASE_DIR = Path(".").resolve()
OUT_DIR = BASE_DIR / "outputs/v2_guide"

PDF_PATH = OUT_DIR / "guida_operativa_v2_userdata_gae.pdf"
TXT_PATH = OUT_DIR / "guida_operativa_v2_userdata_gae.txt"

MANIFEST_PATH = BASE_DIR / "outputs/v2_package/v2_package_manifest.json"
GAE_SUMMARY_PATH = BASE_DIR / "outputs/gae/gae_summary.txt"


def latest_user_dir() -> Path | None:
    root = BASE_DIR / "outputs/user_data"

    if not root.exists():
        return None

    dirs = [p for p in root.iterdir() if p.is_dir()]

    if not dirs:
        return None

    return sorted(dirs, key=lambda p: p.stat().st_mtime, reverse=True)[0]


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}

    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def clean(text: str) -> str:
    replacements = {
        "è": "e'",
        "é": "e'",
        "à": "a'",
        "ò": "o'",
        "ù": "u'",
        "ì": "i'",
        "È": "E'",
        "–": "-",
        "—": "-",
        "“": '"',
        "”": '"',
        "’": "'",
    }

    for old, new in replacements.items():
        text = text.replace(old, new)

    return text


def html_escape(text: str) -> str:
    return (
        clean(str(text))
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace("\n", "<br/>")
    )


def paragraph(story, text: str, style) -> None:
    story.append(Paragraph(html_escape(text), style))


def section(story, title: str, body: list[str], h_style, b_style) -> None:
    paragraph(story, title, h_style)
    story.append(Spacer(1, 0.15 * cm))

    for item in body:
        paragraph(story, item, b_style)
        story.append(Spacer(1, 0.08 * cm))

    story.append(Spacer(1, 0.25 * cm))


def simple_table(story, headers: list[str], rows: list[list[str]]) -> None:
    data = [headers] + rows

    table = Table(data, colWidths=[5.0 * cm, 10.5 * cm])

    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                ("FONTSIZE", (0, 0), (-1, -1), 8),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )

    story.append(table)
    story.append(Spacer(1, 0.35 * cm))


def build_txt(manifest: dict, gae_summary: str, user_summary: str, rag_summary: str) -> str:
    lines = []

    lines.append("GUIDA OPERATIVA V2 - maize ncRNA GNN userdata GAE")
    lines.append("=" * 80)
    lines.append(f"Data generazione: {datetime.now().isoformat(timespec='seconds')}")
    lines.append("")

    lines.append("1. SCOPO DELLA V2")
    lines.append("-" * 80)
    lines.append("La V2 estende il progetto originale con rete neurale Graph Autoencoder,")
    lines.append("caricamento dati del ricercatore, ranking integrato e Chat dati ricercatore.")
    lines.append("")

    lines.append("2. FUNZIONI AGGIUNTE")
    lines.append("-" * 80)
    lines.append("1. Modulo rete neurale Graph Autoencoder.")
    lines.append("2. Pulsante GUI per avviare il training GAE.")
    lines.append("3. Upload dati ricercatore.")
    lines.append("4. Validazione CSV, TSV e FASTA.")
    lines.append("5. Collegamento dati ricercatore al ranking.")
    lines.append("6. Pagina Ranking dati ricercatore.")
    lines.append("7. Knowledge base RAG dati ricercatore.")
    lines.append("8. Chat dati ricercatore.")
    lines.append("9. Docker V2 validato.")
    lines.append("")

    lines.append("3. AVVIO LOCALE")
    lines.append("-" * 80)
    lines.append("cd ~/progetti/ncrna_rnacentral_ai_v2_userdata_gae")
    lines.append("source ~/progetti/ncrna_rnacentral_ai/.venv/bin/activate")
    lines.append("./run_project.sh")
    lines.append("")
    lines.append("GUI: http://localhost:8501")
    lines.append("API: http://localhost:8000/docs")
    lines.append("")

    lines.append("4. AVVIO DOCKER")
    lines.append("-" * 80)
    lines.append("cd ~/progetti/ncrna_rnacentral_ai_v2_userdata_gae")
    lines.append("./run_docker.sh")
    lines.append("./validate_docker.sh")
    lines.append("./stop_docker.sh")
    lines.append("")

    lines.append("5. FORMATO DATI RICERCATORE")
    lines.append("-" * 80)
    lines.append("File consigliato: user_differential_expression.csv")
    lines.append("Colonne minime: rnacentral_id, condition, log2FC, pvalue, FDR")
    lines.append("Esempio:")
    lines.append("URS0003C75441,drought_vs_control,1.82,0.0004,0.012")
    lines.append("")

    lines.append("6. STATO PACCHETTO V2")
    lines.append("-" * 80)
    lines.append(f"Status: {manifest.get('status', 'non disponibile')}")
    lines.append(f"Files count: {manifest.get('files_count', 'non disponibile')}")
    lines.append(f"SHA256: {manifest.get('sha256', 'non disponibile')}")
    lines.append("")

    lines.append("7. RIEPILOGO GAE")
    lines.append("-" * 80)
    lines.append(gae_summary[:2500] if gae_summary else "gae_summary.txt non trovato")
    lines.append("")

    lines.append("8. RIEPILOGO INTEGRAZIONE DATI RICERCATORE")
    lines.append("-" * 80)
    lines.append(user_summary[:3000] if user_summary else "user_ranking_integration_summary.txt non trovato")
    lines.append("")

    lines.append("9. RIEPILOGO RAG UTENTE")
    lines.append("-" * 80)
    lines.append(rag_summary[:2000] if rag_summary else "user_rag_knowledge_summary.txt non trovato")
    lines.append("")

    lines.append("10. LIMITI METODOLOGICI")
    lines.append("-" * 80)
    lines.append("La V2 e' sperimentale.")
    lines.append("Il modulo GAE e' non supervisionato.")
    lines.append("Il ranking con dati utente e' tecnico-esplorativo.")
    lines.append("La presenza di un ncRNA nel dataset del ricercatore non dimostra da sola un ruolo biologico nella risposta alla siccita'.")
    lines.append("Servono dataset reali, repliche biologiche, analisi differenziale robusta e validazione sperimentale.")
    lines.append("")

    lines.append("11. INTEGRAZIONE FUTURA PRJNA419326")
    lines.append("-" * 80)
    lines.append("Il progetto PRJNA419326 non e' ancora stato integrato.")
    lines.append("Gli step futuri sono: scaricare metadati ENA, identificare controllo/siccita', scaricare FASTQ o count table,")
    lines.append("eseguire QC, mappare al genoma Zea mays, generare tabella differenziale e caricarla nella V2.")
    lines.append("")

    return "\n".join(lines) + "\n"


def build_pdf(txt: str, manifest: dict) -> None:
    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        "TitleCustom",
        parent=styles["Title"],
        alignment=TA_CENTER,
        fontName="Helvetica-Bold",
        fontSize=18,
        leading=22,
        spaceAfter=18,
    )

    h_style = ParagraphStyle(
        "HeadingCustom",
        parent=styles["Heading1"],
        alignment=TA_LEFT,
        fontName="Helvetica-Bold",
        fontSize=13,
        leading=16,
        spaceBefore=10,
        spaceAfter=8,
    )

    b_style = ParagraphStyle(
        "BodyCustom",
        parent=styles["BodyText"],
        alignment=TA_LEFT,
        fontName="Helvetica",
        fontSize=9.2,
        leading=12,
        spaceAfter=4,
    )

    code_style = ParagraphStyle(
        "CodeCustom",
        parent=styles["Code"],
        fontName="Courier",
        fontSize=7.4,
        leading=9,
        leftIndent=8,
        rightIndent=8,
        spaceAfter=6,
    )

    doc = SimpleDocTemplate(
        str(PDF_PATH),
        pagesize=A4,
        leftMargin=1.6 * cm,
        rightMargin=1.6 * cm,
        topMargin=1.5 * cm,
        bottomMargin=1.5 * cm,
        title="Guida operativa V2 userdata GAE",
    )

    story = []

    paragraph(story, "Guida operativa V2", title_style)
    paragraph(story, "maize ncRNA GNN userdata GAE", title_style)
    paragraph(story, f"Generata il {datetime.now().isoformat(timespec='seconds')}", b_style)

    story.append(Spacer(1, 0.4 * cm))

    section(
        story,
        "1. Scopo della V2",
        [
            "La V2 estende il progetto originale con rete neurale Graph Autoencoder, caricamento dati del ricercatore, ranking integrato e Chat dati ricercatore.",
            "La piattaforma resta sperimentale: i risultati servono a prioritizzare candidati ncRNA, non a dimostrare automaticamente un ruolo biologico nella risposta alla siccita'.",
        ],
        h_style,
        b_style,
    )

    paragraph(story, "2. Componenti principali", h_style)
    simple_table(
        story,
        ["Componente", "Funzione"],
        [
            ["Pipeline AI/GNN", "Costruzione dataset, feature, grafo, embedding e ranking esplorativo."],
            ["Graph Autoencoder", "Modulo neurale non supervisionato per embedding e ricostruzione link."],
            ["Dati ricercatore", "Upload e validazione di CSV, TSV e FASTA."],
            ["Ranking integrato", "Collegamento tra ncRNA del ranking e dati caricati dal ricercatore."],
            ["Chat dati ricercatore", "Interrogazione dei risultati caricati e integrati."],
            ["Docker", "Esecuzione containerizzata di API e GUI."],
        ],
    )

    section(
        story,
        "3. Avvio locale",
        [
            "cd ~/progetti/ncrna_rnacentral_ai_v2_userdata_gae",
            "source ~/progetti/ncrna_rnacentral_ai/.venv/bin/activate",
            "./run_project.sh",
            "GUI: http://localhost:8501",
            "API: http://localhost:8000/docs",
        ],
        h_style,
        code_style,
    )

    section(
        story,
        "4. Avvio Docker",
        [
            "cd ~/progetti/ncrna_rnacentral_ai_v2_userdata_gae",
            "./run_docker.sh",
            "./validate_docker.sh",
            "./stop_docker.sh",
        ],
        h_style,
        code_style,
    )

    paragraph(story, "5. Pagine GUI della V2", h_style)
    simple_table(
        story,
        ["Pagina", "Uso"],
        [
            ["Rete neurale GAE", "Avvia training Graph Autoencoder e visualizza metriche."],
            ["Dati ricercatore", "Carica CSV, TSV e FASTA."],
            ["Ranking dati ricercatore", "Mostra i match tra ranking e dati caricati."],
            ["Chat dati ricercatore", "Interroga i risultati caricati tramite chat."],
        ],
    )

    section(
        story,
        "6. Formato dati ricercatore",
        [
            "File consigliato: user_differential_expression.csv",
            "Colonne minime: rnacentral_id, condition, log2FC, pvalue, FDR",
            "Esempio: URS0003C75441,drought_vs_control,1.82,0.0004,0.012",
            "Colonne aggiuntive consigliate: gene_id, baseMean, tissue, genotype, time_point, replicate.",
        ],
        h_style,
        b_style,
    )

    section(
        story,
        "7. Stato pacchetto V2",
        [
            f"Status: {manifest.get('status', 'non disponibile')}",
            f"Files count: {manifest.get('files_count', 'non disponibile')}",
            f"SHA256: {manifest.get('sha256', 'non disponibile')}",
        ],
        h_style,
        b_style,
    )

    story.append(PageBreak())

    section(
        story,
        "8. Interpretazione dei risultati",
        [
            "Il ranking integrato indica quali ncRNA del ranking esplorativo sono presenti anche nei dati caricati dal ricercatore.",
            "Un FDR basso aumenta la priorita' del candidato, ma non dimostra causalita'.",
            "Un log2FC elevato indica una variazione marcata tra condizioni, ma richiede controllo di qualita' e validazione.",
            "La Chat dati ricercatore restituisce risposte basate sui file caricati e sul ranking integrato.",
        ],
        h_style,
        b_style,
    )

    section(
        story,
        "9. Limiti metodologici",
        [
            "La V2 e' sperimentale.",
            "Il Graph Autoencoder e' non supervisionato.",
            "Il ranking con dati utente e' tecnico-esplorativo.",
            "La presenza di un ncRNA nei dati caricati non dimostra da sola un ruolo biologico nella risposta alla siccita'.",
            "Servono dataset reali, repliche biologiche, analisi differenziale robusta, dataset indipendenti e validazione sperimentale.",
        ],
        h_style,
        b_style,
    )

    section(
        story,
        "10. Integrazione futura PRJNA419326",
        [
            "Il progetto PRJNA419326 non e' ancora stato integrato nella V2.",
            "Gli step futuri sono: scaricare metadati ENA, identificare controllo e siccita', scaricare FASTQ o count table, eseguire QC, mappare al genoma Zea mays, generare tabella differenziale e caricarla nella GUI V2.",
        ],
        h_style,
        b_style,
    )

    doc.build(story)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    manifest = read_json(MANIFEST_PATH)
    gae_summary = read_text(GAE_SUMMARY_PATH)

    user_dir = latest_user_dir()

    user_summary = ""
    rag_summary = ""

    if user_dir:
        user_summary = read_text(user_dir / "user_ranking_integration_summary.txt")
        rag_summary = read_text(user_dir / "user_rag_knowledge_summary.txt")

    txt = build_txt(manifest, gae_summary, user_summary, rag_summary)

    TXT_PATH.write_text(txt, encoding="utf-8")
    build_pdf(txt, manifest)

    print("GUIDA V2 GENERATA")
    print("=" * 80)
    print(f"TXT: {TXT_PATH}")
    print(f"PDF: {PDF_PATH}")
    print(f"PDF size bytes: {PDF_PATH.stat().st_size}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
