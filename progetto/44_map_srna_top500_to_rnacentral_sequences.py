#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path


PROJECT = "PRJNA419326"
BASE_DIR = Path(".").resolve()
OUT_DIR = BASE_DIR / "outputs/prjna419326"

TOP500_TSV = OUT_DIR / f"{PROJECT}_mirnaseq_srna_drought_response_top500.tsv"

CANDIDATE_FASTA_PATHS = [
    BASE_DIR / "data/rnacentral/maize_sequences.fasta",
    BASE_DIR / "data/public/rnacentral/maize_sequences.fasta",
    BASE_DIR / "outputs/rnacentral/maize_sequences.fasta",
    BASE_DIR / "maize_sequences.fasta",
]

MATCHES_TSV = OUT_DIR / f"{PROJECT}_top500_srna_to_rnacentral_exact_matches.tsv"
UNMATCHED_TSV = OUT_DIR / f"{PROJECT}_top500_srna_unmatched.tsv"
SUMMARY_TXT = OUT_DIR / f"{PROJECT}_top500_srna_to_rnacentral_mapping_summary.txt"
REPORT_JSON = OUT_DIR / f"{PROJECT}_top500_srna_to_rnacentral_mapping_report.json"


def find_fasta() -> Path:
    for p in CANDIDATE_FASTA_PATHS:
        if p.exists():
            return p

    hits = sorted(BASE_DIR.glob("**/maize_sequences.fasta"))
    if hits:
        return hits[0]

    raise FileNotFoundError(
        "Non trovo maize_sequences.fasta nel progetto. "
        "Verifica dove sono state salvate le sequenze RNAcentral."
    )


def parse_fasta(path: Path) -> dict[str, str]:
    records = {}
    current_id = None
    chunks = []

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(chunks).upper().replace("U", "T")

                header = line[1:]
                current_id = header.split()[0].split("|")[0]
                chunks = []
            else:
                chunks.append(line)

    if current_id is not None:
        records[current_id] = "".join(chunks).upper().replace("U", "T")

    return records


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")
    return seq.translate(table)[::-1].upper().replace("U", "T")


def read_top500(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fasta = find_fasta()
    rnacentral = parse_fasta(fasta)
    top = read_top500(TOP500_TSV)

    print(f"MAPPING top500 small RNA su RNAcentral")
    print("=" * 80)
    print(f"Top500 input: {TOP500_TSV}")
    print(f"RNAcentral FASTA: {fasta}")
    print(f"Record RNAcentral: {len(rnacentral)}")
    print(f"Feature small RNA: {len(top)}")
    print("")

    matches = []
    unmatched = []

    for r in top:
        seq = r["sequence"].upper().replace("U", "T")
        rc = revcomp(seq)

        local_matches = []

        for rid, refseq in rnacentral.items():
            pos_fwd = refseq.find(seq)
            if pos_fwd >= 0:
                local_matches.append(
                    {
                        **r,
                        "rnacentral_id": rid,
                        "match_type": "forward_substring",
                        "start_0based": pos_fwd,
                        "end_0based_exclusive": pos_fwd + len(seq),
                        "matched_sequence": seq,
                    }
                )

            pos_rev = refseq.find(rc)
            if pos_rev >= 0:
                local_matches.append(
                    {
                        **r,
                        "rnacentral_id": rid,
                        "match_type": "reverse_complement_substring",
                        "start_0based": pos_rev,
                        "end_0based_exclusive": pos_rev + len(seq),
                        "matched_sequence": rc,
                    }
                )

        if local_matches:
            matches.extend(local_matches)
        else:
            unmatched.append(r)

    match_fieldnames = list(top[0].keys()) + [
        "rnacentral_id",
        "match_type",
        "start_0based",
        "end_0based_exclusive",
        "matched_sequence",
    ]

    with MATCHES_TSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=match_fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(matches)

    with UNMATCHED_TSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(top[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(unmatched)

    matched_features = sorted({m["feature_id"] for m in matches})
    matched_rnacentral = sorted({m["rnacentral_id"] for m in matches})

    class_counts = {}
    for m in matches:
        cls = m.get("candidate_class", "")
        class_counts[cls] = class_counts.get(cls, 0) + 1

    report = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "top500_input": str(TOP500_TSV),
        "rnacentral_fasta": str(fasta),
        "n_rnacentral_records": len(rnacentral),
        "n_top_features": len(top),
        "n_match_rows": len(matches),
        "n_matched_features": len(matched_features),
        "n_unmatched_features": len(unmatched),
        "n_matched_rnacentral_ids": len(matched_rnacentral),
        "matched_rnacentral_ids": matched_rnacentral,
        "class_counts_in_matches": class_counts,
        "outputs": {
            "matches_tsv": str(MATCHES_TSV),
            "unmatched_tsv": str(UNMATCHED_TSV),
            "summary_txt": str(SUMMARY_TXT),
            "report_json": str(REPORT_JSON),
        },
    }

    REPORT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("MAPPING top500 small RNA PRJNA419326 su RNAcentral")
    lines.append("=" * 80)
    lines.append(f"Data: {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"RNAcentral FASTA: {fasta}")
    lines.append(f"Record RNAcentral: {len(rnacentral)}")
    lines.append(f"Feature small RNA top500: {len(top)}")
    lines.append(f"Righe match: {len(matches)}")
    lines.append(f"Feature small RNA con match: {len(matched_features)}")
    lines.append(f"Feature small RNA senza match: {len(unmatched)}")
    lines.append(f"RNAcentral ID coinvolti: {len(matched_rnacentral)}")
    lines.append("")
    lines.append("Classi nei match:")
    for k in sorted(class_counts):
        lines.append(f"{k}: {class_counts[k]}")
    lines.append("")
    lines.append("Primi RNAcentral ID coinvolti:")
    for rid in matched_rnacentral[:30]:
        lines.append(f"- {rid}")
    lines.append("")
    lines.append("Output:")
    lines.append(str(MATCHES_TSV))
    lines.append(str(UNMATCHED_TSV))
    lines.append(str(REPORT_JSON))

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
