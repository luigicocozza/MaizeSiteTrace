#!/usr/bin/env python3
"""
17_validate_end_to_end.py

Validazione end-to-end del progetto maize_drought_ncrna_gnn.

Controlla:
- file principali;
- compilazione script principali;
- API FastAPI;
- GUI Streamlit raggiungibile;
- endpoint candidati;
- endpoint RAG;
- coerenza dei numeri principali.

Output:
- outputs/e2e/e2e_validation_report.json
- outputs/e2e/e2e_validation_summary.txt
"""

from __future__ import annotations

import argparse
import json
import py_compile
import sys
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Tuple


def parse_args():
    p = argparse.ArgumentParser(description="Validazione end-to-end progetto.")
    p.add_argument("--api-url", default="http://localhost:8000")
    p.add_argument("--gui-url", default="http://localhost:8501")
    p.add_argument("--output-dir", default="outputs/e2e")
    p.add_argument("--verbose", action="store_true")
    return p.parse_args()


def http_get_json(url: str, timeout: int = 10) -> Tuple[bool, Any, str]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            raw = r.read().decode("utf-8", errors="replace")
            return True, json.loads(raw), ""
    except Exception as exc:
        return False, None, str(exc)


def http_post_json(url: str, payload: Dict[str, Any], timeout: int = 20) -> Tuple[bool, Any, str]:
    try:
        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(
            url,
            data=data,
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=timeout) as r:
            raw = r.read().decode("utf-8", errors="replace")
            return True, json.loads(raw), ""
    except Exception as exc:
        return False, None, str(exc)


def http_get_status(url: str, timeout: int = 10) -> Tuple[bool, int | None, str]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            return r.status < 500, r.status, ""
    except Exception as exc:
        return False, None, str(exc)


def check_files() -> List[Dict[str, Any]]:
    paths = [
        "app/streamlit_app.py",
        "api/main.py",
        "16_rag_chat_engine.py",
        "run_gui.sh",
        "run_api.sh",
        "run_project.sh",
        "stop_project.sh",
        "data/processed/ranking/ncrna_exploratory_ranking.csv",
        "data/processed/similarity/ncrna_similarity_topk.csv",
        "data/processed/embeddings/ncrna_embeddings.csv",
        "data/processed/ncrna_features.csv",
        "data/processed/graph_nodes.csv",
        "data/processed/graph_edges.csv",
        "data/processed/rag/project_knowledge_base.jsonl",
        "outputs/reports/methodological_report.md",
        "outputs/reports/methodological_report.json",
    ]

    rows = []

    for item in paths:
        p = Path(item)
        rows.append(
            {
                "path": item,
                "exists": p.exists(),
                "size_bytes": p.stat().st_size if p.exists() else 0,
            }
        )

    return rows


def check_py_compile() -> List[Dict[str, Any]]:
    scripts = [
        "app/streamlit_app.py",
        "api/main.py",
        "16_rag_chat_engine.py",
        "17_validate_end_to_end.py",
    ]

    rows = []

    for script in scripts:
        try:
            py_compile.compile(script, doraise=True)
            rows.append({"script": script, "ok": True, "error": ""})
        except Exception as exc:
            rows.append({"script": script, "ok": False, "error": str(exc)})

    return rows


def validate_api(api_url: str) -> Dict[str, Any]:
    results: Dict[str, Any] = {}

    ok, data, err = http_get_json(f"{api_url}/health")
    results["health"] = {"ok": ok, "data": data, "error": err}

    ok, data, err = http_get_json(f"{api_url}/summary")
    results["summary"] = {"ok": ok, "data": data, "error": err}

    ok, data, err = http_get_json(f"{api_url}/candidates?limit=5")
    results["candidates"] = {"ok": ok, "data": data, "error": err}

    ok, data, err = http_get_json(f"{api_url}/candidate/URS0003C75441")
    results["candidate_detail"] = {"ok": ok, "data": data, "error": err}

    ok, data, err = http_get_json(f"{api_url}/candidate/URS0003C75441/similar?limit=5")
    results["candidate_similar"] = {"ok": ok, "data": data, "error": err}

    ok, data, err = http_post_json(
        f"{api_url}/rag/ask",
        {
            "question": "Quanti ncRNA sono stati analizzati e qual è il limite metodologico principale?",
            "provider": "extractive",
            "top_k": 6,
        },
    )
    results["rag_ask"] = {"ok": ok, "data": data, "error": err}

    return results


def derive_checks(api_results: Dict[str, Any], gui_ok: bool, file_rows: List[Dict[str, Any]], compile_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    summary = api_results.get("summary", {}).get("data") or {}
    candidates = api_results.get("candidates", {}).get("data") or {}
    candidate_detail = api_results.get("candidate_detail", {}).get("data") or {}
    similar = api_results.get("candidate_similar", {}).get("data") or {}
    rag = api_results.get("rag_ask", {}).get("data") or {}

    checks = {
        "all_required_files_exist": all(row["exists"] for row in file_rows),
        "all_scripts_compile": all(row["ok"] for row in compile_rows),
        "api_health_ok": bool(api_results.get("health", {}).get("ok")),
        "api_summary_ok": bool(api_results.get("summary", {}).get("ok")),
        "gui_http_ok": gui_ok,
        "summary_n_ranked_ncrna_250": summary.get("n_ranked_ncrna") == 250,
        "summary_graph_nodes_384": summary.get("n_graph_nodes") == 384,
        "summary_graph_edges_1134": summary.get("n_graph_edges") == 1134,
        "candidates_count_5": candidates.get("count") == 5,
        "candidate_detail_id_ok": (candidate_detail.get("candidate") or {}).get("rnacentral_id") == "URS0003C75441",
        "candidate_similar_count_5": similar.get("count") == 5,
        "rag_provider_ok": rag.get("provider") in {"structured_rag", "extractive", "extractive_fallback"},
        "rag_mentions_250": "250" in str(rag.get("answer", "")),
    }

    checks["overall_pass"] = all(checks.values())
    return checks


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    file_rows = check_files()
    compile_rows = check_py_compile()

    api_results = validate_api(args.api_url)
    gui_ok, gui_status, gui_error = http_get_status(args.gui_url)

    checks = derive_checks(
        api_results=api_results,
        gui_ok=gui_ok,
        file_rows=file_rows,
        compile_rows=compile_rows,
    )

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "api_url": args.api_url,
        "gui_url": args.gui_url,
        "files": file_rows,
        "py_compile": compile_rows,
        "api_results": api_results,
        "gui": {
            "ok": gui_ok,
            "status": gui_status,
            "error": gui_error,
        },
        "checks": checks,
        "validation_status": "pass" if checks["overall_pass"] else "fail",
    }

    report_path = output_dir / "e2e_validation_report.json"
    summary_path = output_dir / "e2e_validation_summary.txt"

    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("END-TO-END VALIDATION SUMMARY")
    lines.append("=" * 80)
    lines.append(f"Validation status: {report['validation_status']}")
    lines.append(f"API URL: {args.api_url}")
    lines.append(f"GUI URL: {args.gui_url}")
    lines.append("")
    lines.append("Checks:")
    for key, value in checks.items():
        lines.append(f"- {key}: {value}")
    lines.append("")
    lines.append(f"Report JSON: {report_path}")
    lines.append(f"Summary TXT: {summary_path}")

    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0 if checks["overall_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
