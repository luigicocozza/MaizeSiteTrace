#!/usr/bin/env python3
from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path


BASE_DIR = Path(".").resolve()
OUT_DIR = BASE_DIR / "outputs/prjna419326"
OUT_DIR.mkdir(parents=True, exist_ok=True)

REPORT_JSON = OUT_DIR / "V3_gui_project_realization_update_report.json"
IDENTITY_TXT = OUT_DIR / "V3_project_realization_identity.txt"

AUTHOR_LINES = [
    "Realizzazione progetto",
    "Ing. Luigi Cocozza",
    "Esperto in: AI Machine Learning Engineer",
    "Ingegnere di Machine Learning IA",
]

MARKER_BEGIN = "# BEGIN V3_PROJECT_REALIZATION_BLOCK"
MARKER_END = "# END V3_PROJECT_REALIZATION_BLOCK"

EXCLUDE_DIRS = {
    ".git",
    ".venv",
    "__pycache__",
    "delivery",
    "data",
    "models",
    "outputs",
    "reference",
    "rnaseq_fastq",
    "rnaseq_salmon_quant",
}


STREAMLIT_BLOCK = f'''

{MARKER_BEGIN}
def _render_v3_project_realization_block():
    """
    Blocco informativo aggiunto alla GUI V3.
    """
    try:
        import streamlit as st

        st.sidebar.markdown("---")
        st.sidebar.markdown("### Realizzazione progetto")
        st.sidebar.markdown("**Ing. Luigi Cocozza**")
        st.sidebar.markdown("Esperto in: **AI Machine Learning Engineer**")
        st.sidebar.markdown("_Ingegnere di Machine Learning IA_")
    except Exception:
        pass


_render_v3_project_realization_block()
{MARKER_END}
'''


GRADIO_MARKDOWN = '''gr.Markdown("""
---
### Realizzazione progetto

**Ing. Luigi Cocozza**  
Esperto in: **AI Machine Learning Engineer**  
*Ingegnere di Machine Learning IA*
""")
'''


HTML_BLOCK = '''
<!-- BEGIN V3_PROJECT_REALIZATION_BLOCK -->
<footer style="margin-top: 2rem; padding: 1rem; border-top: 1px solid #ddd; font-size: 0.95rem;">
  <strong>Realizzazione progetto</strong><br>
  Ing. Luigi Cocozza<br>
  Esperto in: <strong>AI Machine Learning Engineer</strong><br>
  <em>Ingegnere di Machine Learning IA</em>
</footer>
<!-- END V3_PROJECT_REALIZATION_BLOCK -->
'''


def is_excluded(path: Path) -> bool:
    return any(part in EXCLUDE_DIRS for part in path.parts)


def backup_file(path: Path) -> Path:
    backup = path.with_suffix(path.suffix + ".bak_before_author_v3")
    if not backup.exists():
        shutil.copy2(path, backup)
    return backup


def patch_streamlit_file(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")

    if MARKER_BEGIN in text:
        return {
            "file": str(path),
            "type": "streamlit",
            "status": "already_present",
            "backup": None,
        }

    if "streamlit" not in text and "st." not in text:
        return {
            "file": str(path),
            "type": "streamlit",
            "status": "not_streamlit",
            "backup": None,
        }

    backup = backup_file(path)
    path.write_text(text.rstrip() + STREAMLIT_BLOCK + "\n", encoding="utf-8")

    return {
        "file": str(path),
        "type": "streamlit",
        "status": "patched",
        "backup": str(backup),
    }


def patch_gradio_file(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")

    if MARKER_BEGIN in text:
        return {
            "file": str(path),
            "type": "gradio",
            "status": "already_present",
            "backup": None,
        }

    if "gr.Blocks" not in text and "gradio" not in text:
        return {
            "file": str(path),
            "type": "gradio",
            "status": "not_gradio",
            "backup": None,
        }

    lines = text.splitlines()
    patched_lines = []
    inserted = False

    for line in lines:
        patched_lines.append(line)

        if not inserted and "with gr.Blocks" in line:
            indent = line[: len(line) - len(line.lstrip())] + "    "
            patched_lines.append(indent + f"{MARKER_BEGIN}")
            for md_line in GRADIO_MARKDOWN.splitlines():
                patched_lines.append(indent + md_line)
            patched_lines.append(indent + f"{MARKER_END}")
            inserted = True

    if not inserted:
        return {
            "file": str(path),
            "type": "gradio",
            "status": "gradio_blocks_not_found",
            "backup": None,
        }

    backup = backup_file(path)
    path.write_text("\n".join(patched_lines) + "\n", encoding="utf-8")

    return {
        "file": str(path),
        "type": "gradio",
        "status": "patched",
        "backup": str(backup),
    }


def patch_html_file(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")

    if "BEGIN V3_PROJECT_REALIZATION_BLOCK" in text:
        return {
            "file": str(path),
            "type": "html",
            "status": "already_present",
            "backup": None,
        }

    backup = backup_file(path)

    lower = text.lower()
    if "</body>" in lower:
        idx = lower.rfind("</body>")
        new_text = text[:idx] + HTML_BLOCK + "\n" + text[idx:]
    else:
        new_text = text.rstrip() + HTML_BLOCK + "\n"

    path.write_text(new_text, encoding="utf-8")

    return {
        "file": str(path),
        "type": "html",
        "status": "patched",
        "backup": str(backup),
    }


def find_gui_files() -> list[Path]:
    candidates = []

    for path in BASE_DIR.rglob("*"):
        if not path.is_file():
            continue

        if is_excluded(path.relative_to(BASE_DIR)):
            continue

        name = path.name.lower()
        rel = str(path.relative_to(BASE_DIR)).lower()

        if path.suffix == ".py":
            try:
                text = path.read_text(encoding="utf-8", errors="replace")
            except Exception:
                continue

            is_gui_name = any(k in rel for k in ["gui", "app", "streamlit", "dashboard", "chat", "rag"])
            is_gui_content = any(k in text for k in ["streamlit", "import streamlit", "gradio", "import gradio", "gr.Blocks"])

            if is_gui_name and is_gui_content:
                candidates.append(path)

        elif path.suffix.lower() in {".html", ".htm"}:
            is_gui_name = any(k in rel for k in ["gui", "app", "dashboard", "template", "static"])
            if is_gui_name:
                candidates.append(path)

    return sorted(set(candidates))


def main() -> int:
    IDENTITY_TXT.write_text(
        "\n".join(AUTHOR_LINES) + "\n",
        encoding="utf-8",
    )

    candidates = find_gui_files()

    results = []

    for path in candidates:
        text = path.read_text(encoding="utf-8", errors="replace") if path.suffix == ".py" else ""

        if path.suffix == ".py" and ("streamlit" in text or "st." in text):
            results.append(patch_streamlit_file(path))
        elif path.suffix == ".py" and ("gradio" in text or "gr.Blocks" in text):
            results.append(patch_gradio_file(path))
        elif path.suffix.lower() in {".html", ".htm"}:
            results.append(patch_html_file(path))

    patched = [r for r in results if r["status"] == "patched"]
    already = [r for r in results if r["status"] == "already_present"]

    status = "PASS" if patched or already else "NO_GUI_FILE_PATCHED"

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "status": status,
        "identity": {
            "realizzazione": "Ing. Luigi Cocozza",
            "profilo": "AI Machine Learning Engineer",
            "profilo_it": "Ingegnere di Machine Learning IA",
        },
        "identity_txt": str(IDENTITY_TXT),
        "n_candidates": len(candidates),
        "n_patched": len(patched),
        "n_already_present": len(already),
        "results": results,
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    print("AGGIORNAMENTO GUI - REALIZZAZIONE PROGETTO")
    print("=" * 80)
    print(f"Status: {status}")
    print(f"File GUI candidati: {len(candidates)}")
    print(f"File modificati: {len(patched)}")
    print(f"File già aggiornati: {len(already)}")
    print("")
    print("Dicitura inserita:")
    for line in AUTHOR_LINES:
        print(f"- {line}")
    print("")
    print("Risultati:")
    if not results:
        print("- Nessun file GUI trovato automaticamente.")
    else:
        for r in results:
            print(f"- {r['status']} | {r['type']} | {r['file']}")
            if r.get("backup"):
                print(f"  backup: {r['backup']}")
    print("")
    print("File identità:")
    print(IDENTITY_TXT)
    print("")
    print("Report:")
    print(REPORT_JSON)

    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
