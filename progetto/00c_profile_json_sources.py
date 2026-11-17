#!/usr/bin/env python3

import json
from pathlib import Path
from collections import Counter

RAW_DIR = Path("data/raw")
OUT_DIR = Path("outputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

files = {
    "rnacentral": RAW_DIR / "rnacentral_maize_ncrna.json",
    "europe_pmc": RAW_DIR / "europe_pmc_maize_drought_ncrna.json",
    "ncbi_geo": RAW_DIR / "ncbi_geo_maize_drought_search.json",
    "quickgo": RAW_DIR / "quickgo_maize_drought_annotations.json",
}

report = {}

def load_json(path):
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)

# RNAcentral
rnacentral = load_json(files["rnacentral"])
if rnacentral:
    results = rnacentral.get("results", [])
    key_counter = Counter()
    rna_types = Counter()
    has_sequence = 0
    sample_records = []

    for rec in results:
        key_counter.update(rec.keys())
        if "sequence" in rec and rec.get("sequence"):
            has_sequence += 1
        rna_type = rec.get("rna_type") or rec.get("rna_type_label") or rec.get("type")
        if rna_type:
            rna_types[str(rna_type)] += 1
        if len(sample_records) < 3:
            sample_records.append(rec)

    report["rnacentral"] = {
        "file": str(files["rnacentral"]),
        "top_level_keys": list(rnacentral.keys()),
        "n_results": len(results),
        "next_page": rnacentral.get("next"),
        "previous_page": rnacentral.get("previous"),
        "record_keys": sorted(list(key_counter.keys())),
        "n_records_with_sequence_field": has_sequence,
        "rna_types": dict(rna_types),
        "sample_records": sample_records,
    }

# Europe PMC
europe = load_json(files["europe_pmc"])
if europe:
    result_list = europe.get("resultList", {}).get("result", [])
    report["europe_pmc"] = {
        "file": str(files["europe_pmc"]),
        "top_level_keys": list(europe.keys()),
        "hitCount": europe.get("hitCount"),
        "n_downloaded_results": len(result_list),
        "sample_titles": [
            r.get("title") for r in result_list[:5]
        ],
        "sample_years": [
            r.get("pubYear") for r in result_list[:5]
        ],
    }

# NCBI GEO
geo = load_json(files["ncbi_geo"])
if geo:
    esearch = geo.get("esearchresult", {})
    report["ncbi_geo"] = {
        "file": str(files["ncbi_geo"]),
        "top_level_keys": list(geo.keys()),
        "count": esearch.get("count"),
        "retmax": esearch.get("retmax"),
        "idlist": esearch.get("idlist", [])[:20],
    }

# QuickGO
quickgo = load_json(files["quickgo"])
if quickgo:
    results = quickgo.get("results", [])
    go_aspects = Counter()
    qualifiers = Counter()
    evidence = Counter()
    go_ids = Counter()

    for r in results:
        if r.get("goAspect"):
            go_aspects[r.get("goAspect")] += 1
        if r.get("qualifier"):
            qualifiers[r.get("qualifier")] += 1
        if r.get("goEvidence"):
            evidence[r.get("goEvidence")] += 1
        if r.get("goId"):
            go_ids[r.get("goId")] += 1

    report["quickgo"] = {
        "file": str(files["quickgo"]),
        "top_level_keys": list(quickgo.keys()),
        "numberOfHits": quickgo.get("numberOfHits"),
        "n_downloaded_results": len(results),
        "go_aspects": dict(go_aspects),
        "qualifiers": dict(qualifiers),
        "evidence_codes": dict(evidence),
        "top_go_ids": dict(go_ids.most_common(10)),
        "sample_records": results[:3],
    }

out_path = OUT_DIR / "json_sources_profile.json"
with out_path.open("w", encoding="utf-8") as f:
    json.dump(report, f, indent=2, ensure_ascii=False)

print("Report scritto in:", out_path)
print()

print("=" * 80)
print("RNAcentral")
if "rnacentral" in report:
    r = report["rnacentral"]
    print("N risultati:", r["n_results"])
    print("Next page:", r["next_page"])
    print("Campi record:", r["record_keys"])
    print("Record con campo sequence:", r["n_records_with_sequence_field"])
    print("Tipi RNA:", r["rna_types"])
else:
    print("Non disponibile")

print("=" * 80)
print("Europe PMC")
if "europe_pmc" in report:
    r = report["europe_pmc"]
    print("Hit count:", r["hitCount"])
    print("Risultati scaricati:", r["n_downloaded_results"])
    print("Titoli esempio:")
    for title in r["sample_titles"]:
        print("-", title)
else:
    print("Non disponibile")

print("=" * 80)
print("NCBI GEO")
if "ncbi_geo" in report:
    r = report["ncbi_geo"]
    print("Count:", r["count"])
    print("ID esempio:", r["idlist"])
else:
    print("Non disponibile")

print("=" * 80)
print("QuickGO")
if "quickgo" in report:
    r = report["quickgo"]
    print("NumberOfHits:", r["numberOfHits"])
    print("Risultati scaricati:", r["n_downloaded_results"])
    print("GO aspects:", r["go_aspects"])
    print("Evidence codes:", r["evidence_codes"])
else:
    print("Non disponibile")
