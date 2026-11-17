#!/usr/bin/env python3

import json
from pathlib import Path

RAW_DIR = Path("data/raw")
OUT_DIR = Path("outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

report = []

for path in sorted(RAW_DIR.iterdir()):
    if not path.is_file():
        continue

    info = {
        "filename": path.name,
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "detected_type": None,
        "json_valid": False,
        "json_top_level_type": None,
        "json_keys": [],
        "preview": "",
        "notes": []
    }

    try:
        text = path.read_text(encoding="utf-8", errors="replace")
        stripped = text.strip()
        info["preview"] = stripped[:300].replace("\n", " ")

        if stripped.startswith("{") or stripped.startswith("["):
            info["detected_type"] = "possible_json"
            try:
                obj = json.loads(stripped)
                info["json_valid"] = True
                info["json_top_level_type"] = type(obj).__name__

                if isinstance(obj, dict):
                    info["json_keys"] = list(obj.keys())[:20]
                elif isinstance(obj, list):
                    info["notes"].append(f"lista_json_con_{len(obj)}_elementi")

            except Exception as e:
                info["notes"].append(f"json_non_valido: {e}")

        elif stripped.lower().startswith("<!doctype html") or stripped.lower().startswith("<html"):
            info["detected_type"] = "html_page"
            info["notes"].append("probabile pagina HTML, non dataset strutturato")

        elif "<html" in stripped[:1000].lower():
            info["detected_type"] = "html_page"
            info["notes"].append("probabile pagina HTML, non dataset strutturato")

        else:
            info["detected_type"] = "text_or_unknown"

    except Exception as e:
        info["detected_type"] = "read_error"
        info["notes"].append(str(e))

    report.append(info)

out_path = OUT_DIR / "raw_downloads_inspection.json"
with out_path.open("w", encoding="utf-8") as f:
    json.dump(report, f, indent=2, ensure_ascii=False)

print("Report scritto in:", out_path)
print()

for item in report:
    print("=" * 80)
    print("FILE:", item["filename"])
    print("DIMENSIONE:", item["size_bytes"])
    print("TIPO RILEVATO:", item["detected_type"])
    print("JSON VALIDO:", item["json_valid"])
    print("TIPO JSON:", item["json_top_level_type"])
    print("CHIAVI JSON:", item["json_keys"])
    print("NOTE:", item["notes"])
    print("PREVIEW:", item["preview"][:200])
