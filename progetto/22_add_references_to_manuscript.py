#!/usr/bin/env python3
from __future__ import annotations

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

INPUT_MD = DOCS_DIR / "SCIENTIFIC_MANUSCRIPT_DRAFT.md"
OUTPUT_MD = DOCS_DIR / "SCIENTIFIC_MANUSCRIPT_WITH_REFERENCES.md"
OUTPUT_PDF = OUT_DIR / "scientific_manuscript_with_references.pdf"
OUTPUT_JSON = OUT_DIR / "references_used.json"


REFERENCES = [
    {
        "id": 1,
        "text": "The RNAcentral Consortium. RNAcentral 2021: secondary structure integration, improved sequence search and new member databases. Nucleic Acids Research. 2021;49(D1):D212-D220. doi:10.1093/nar/gkaa921."
    },
    {
        "id": 2,
        "text": "The RNAcentral Consortium. RNAcentral: a comprehensive database of non-coding RNA sequences. Nucleic Acids Research. 2017;45(D1):D128-D134. doi:10.1093/nar/gkw1008."
    },
    {
        "id": 3,
        "text": "Zhang P, Fan Y, Sun X, Chen L, Terzaghi W, Bucher E, Li L, Dai M. A large-scale circular RNA profiling reveals universal molecular mechanisms responsive to drought stress in maize and Arabidopsis. Plant Journal. 2019;98:697-713. doi:10.1111/tpj.14267."
    },
    {
        "id": 4,
        "text": "Pang J, Zhang X, Ma X, Zhao J. Spatio-Temporal Transcriptional Dynamics of Maize Long Non-Coding RNAs Responsive to Drought Stress. Genes. 2019;10:138. doi:10.3390/genes10020138."
    },
    {
        "id": 5,
        "text": "Zhang W, Han Z, Guo Q, Liu Y, Zheng Y, Wu F, Jin W. Identification of maize long non-coding RNAs responsive to drought stress. PLoS ONE. 2014;9:e98958. doi:10.1371/journal.pone.0098958."
    },
    {
        "id": 6,
        "text": "Fey M, Lenssen JE. Fast Graph Representation Learning with PyTorch Geometric. ICLR Workshop on Representation Learning on Graphs and Manifolds. 2019. arXiv:1903.02428."
    },
    {
        "id": 7,
        "text": "Battaglia PW, Hamrick JB, Bapst V, Sanchez-Gonzalez A, Zambaldi V, Malinowski M, Tacchetti A, Raposo D, Santoro A, Faulkner R, et al. Relational inductive biases, deep learning, and graph networks. arXiv. 2018. arXiv:1806.01261."
    },
    {
        "id": 8,
        "text": "Zhou J, Cui G, Hu S, Zhang Z, Yang C, Liu Z, Wang L, Li C, Sun M. Graph neural networks: A review of methods and applications. AI Open. 2020;1:57-81. doi:10.1016/j.aiopen.2021.01.001."
    },
    {
        "id": 9,
        "text": "Lewis P, Perez E, Piktus A, Petroni F, Karpukhin V, Goyal N, Kuttler H, Lewis M, Yih W, Rocktaschel T, Riedel S, Kiela D. Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks. NeurIPS. 2020. arXiv:2005.11401."
    },
    {
        "id": 10,
        "text": "Streamlit Documentation. Streamlit is an open-source Python framework for building dynamic data apps. Accessed during project preparation."
    },
    {
        "id": 11,
        "text": "FastAPI Documentation. FastAPI is a modern, high-performance web framework for building APIs with Python based on standard type hints. Accessed during project preparation."
    },
]


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


def add_in_text_citations(md: str) -> str:
    replacements = {
        "Non-coding RNAs are increasingly relevant in plant regulatory biology,": "Non-coding RNAs are increasingly relevant in plant regulatory biology [3-5],",
        "public maize ncRNA records": "public maize ncRNA records [1,2]",
        "The graph connects ncRNAs to associated entities": "The graph connects ncRNAs to associated entities [7,8]",
        "The PyTorch Geometric representation includes": "The PyTorch Geometric representation includes [6]",
        "A retrieval-augmented generation knowledge base was built": "A retrieval-augmented generation knowledge base was built [9]",
        "A Streamlit graphical interface allows users": "A Streamlit graphical interface allows users [10]",
        "A FastAPI backend provides programmatic access": "A FastAPI backend provides programmatic access [11]",
        "Drought is one of the major constraints affecting crop productivity.": "Drought is one of the major constraints affecting crop productivity, and ncRNA-mediated regulation has been increasingly reported in plant drought responses [3-5].",
        "The project was containerized using Docker.": "The project was containerized using Docker to improve portability and reproducibility.",
    }

    for old, new in replacements.items():
        md = md.replace(old, new)

    return md


def replace_references(md: str) -> str:
    ref_section = ["## References", ""]
    for ref in REFERENCES:
        ref_section.append(f"[{ref['id']}] {ref['text']}")
        ref_section.append("")

    ref_text = "\n".join(ref_section).rstrip() + "\n"

    marker = "## References"
    if marker in md:
        md = md.split(marker)[0].rstrip() + "\n\n" + ref_text
    else:
        md = md.rstrip() + "\n\n" + ref_text

    return md


def draw_header(c: canvas.Canvas, page_num: int) -> None:
    width, height = A4
    c.setFillColor(colors.HexColor("#203864"))
    c.rect(0, height - 42, width, 42, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont("Helvetica-Bold", 10)
    c.drawString(40, height - 26, "Scientific manuscript with references - Maize ncRNA GNN")
    c.setFont("Helvetica", 9)
    c.drawRightString(width - 40, height - 26, f"Page {page_num}")
    c.setFillColor(colors.black)


def write_pdf(markdown: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    c = canvas.Canvas(str(OUTPUT_PDF), pagesize=A4)
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
    if not INPUT_MD.exists():
        raise FileNotFoundError(f"File non trovato: {INPUT_MD}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_DIR.mkdir(parents=True, exist_ok=True)

    md = INPUT_MD.read_text(encoding="utf-8")
    md = add_in_text_citations(md)
    md = replace_references(md)

    OUTPUT_MD.write_text(md, encoding="utf-8")
    write_pdf(md)

    metadata = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "input": str(INPUT_MD),
        "markdown": str(OUTPUT_MD),
        "pdf": str(OUTPUT_PDF),
        "n_references": len(REFERENCES),
        "references": REFERENCES,
    }

    OUTPUT_JSON.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    print("Manuscript with references Markdown:", OUTPUT_MD)
    print("Manuscript with references PDF:", OUTPUT_PDF)
    print("References JSON:", OUTPUT_JSON)
    print("References:", len(REFERENCES))
    print("Completato.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
