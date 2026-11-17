
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 48 — Chiusura scientifica e documentale

Obiettivi:
- Verifica finale dello STEP47 V7
- Consolidamento delle evidenze PARE e siccità
- Verifica esatta delle due small RNA candidate
- Mapping genomico opzionale RefGen_v4 -> NAM-5.0
- Confronto opzionale con miRNA maturi di miRBase
- Audit degli artefatti del progetto
- Rapporto conclusivo e manifest SHA-256

Vincoli:
- Nessun training ML/GNN
- Nessuna modifica alla GUI o a Docker
- Nessuna modifica ai checkpoint
- Nessuna etichetta funzionale inventata
- Nessuna equivalenza genomica dedotta senza fonte

Dipendenze: esclusivamente standard library Python.
"""

import argparse
import csv
import gzip
import hashlib
import json
import re

from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE GENERALE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

RUNS = (
    ROOT
    / "target_prediction_v33_vscode"
    / "runs"
)

CANDIDATES = {
    "sRNA_00000037": "UUCUUGGAUUUAUGAAAGACGAACA",
    "sRNA_00000435": "AGAUGAUACGAGCACCUCUGAGCU",
}

XREF_URL = (
    "https://download.maizegdb.org/"
    "Pan-genes/B73_gene_xref/"
    "B73v4_to_B73v5.tsv"
)

MIRBASE_URL = "https://www.mirbase.org/download/"

V4_RE = re.compile(
    r"(?<![A-Za-z0-9])Zm00001d\d+(?![A-Za-z0-9])",
    re.I,
)

V5_RE = re.compile(
    r"(?<![A-Za-z0-9])Zm00001eb\d+(?![A-Za-z0-9])",
    re.I,
)


# ============================================================
# 2. FUNZIONI GENERALI
# ============================================================

def digest(path):
    """Calcola SHA-256 del file."""

    h = hashlib.sha256()

    with Path(path).open("rb") as f:
        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def read_tsv(path):
    """Legge un TSV conservando i nomi delle colonne."""

    with Path(path).open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as f:

        return list(
            csv.DictReader(f, delimiter="\t")
        )


def write_tsv(path, rows, columns):
    """Scrive un TSV senza modificare gli input."""

    with Path(path).open(
        "w",
        encoding="utf-8",
        newline="",
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def write_json(path, data):
    Path(path).write_text(
        json.dumps(
            data,
            indent=2,
            ensure_ascii=False,
        ) + "\n",
        encoding="utf-8",
    )


def norm_id(value):
    """Normalizza ID per confronti esatti controllati."""

    result = re.sub(
        r"[^a-z0-9]",
        "",
        str(value).lower(),
    )

    return result.removeprefix("zma")


def rna(value):
    """Normalizza una sequenza DNA/RNA in RNA."""

    val = re.sub(
        r"\s",
        "",
        str(value).upper(),
    ).replace("T", "U")

    if val and set(val) <= set("ACGU"):
        return val

    return ""


def assembly(gene):
    """Classificazione per prefisso, non conversione genomica."""

    gene = str(gene)

    if gene.startswith("Zm00001eb"):
        return "NAM5"

    if gene.startswith("Zm00001d"):
        return "RefGen_v4"

    if gene.startswith(("GRMZM2G", "GRMZM5G")):
        return "older_build"

    if gene.lower().startswith("zma-mir"):
        return "MIR_precursor_locus"

    if re.fullmatch(
        r"[A-Z]{2}\d{6}\.\d+_FG\d+",
        gene,
        re.I,
    ):
        return "legacy_FG"

    return "unresolved"


# ============================================================
# 3. CROSSWALK GENOMICO OPZIONALE
# ============================================================

def load_xref(path):
    """
    Legge una tabella esterna di corrispondenza B73v4/B73v5.

    Criterio conservativo:
    una riga viene accettata soltanto se contiene
    esattamente un identificativo v4 e uno v5.

    Le righe con molteplici identificativi non vengono
    trasformate in associazioni cartesiane presunte.
    """

    mapping = defaultdict(set)

    total_rows = 0
    accepted_rows = 0
    ambiguous_rows = 0

    with Path(path).open(
        encoding="utf-8-sig",
        errors="replace",
        newline="",
    ) as f:

        for line in f:

            total_rows += 1

            v4_ids = set(V4_RE.findall(line))
            v5_ids = set(V5_RE.findall(line))

            if not v4_ids or not v5_ids:
                continue

            if len(v4_ids) != 1 or len(v5_ids) != 1:
                ambiguous_rows += 1
                continue

            old = next(iter(v4_ids))
            new = next(iter(v5_ids))

            mapping[old].add(new)
            accepted_rows += 1

    if not mapping:
        raise ValueError(
            "Crosswalk: nessuna coppia v4/v5 univoca "
            "riconosciuta. Verificare lo schema del file."
        )

    details = {
        "rows_read": total_rows,
        "unambiguous_rows_used": accepted_rows,
        "rows_with_multiple_ids_not_assigned": ambiguous_rows,
        "distinct_v4_ids": len(mapping),
        "interpretation": (
            "Corrispondenze di annotazione, "
            "non identità funzionali validate"
        ),
    }

    return mapping, details


# ============================================================
# 4. SEQUENZE MATURE miRBase
# ============================================================

def load_mature(path):
    """
    Legge un FASTA di miRNA maturi.

    Supporta .fa, .fasta e .gz.

    Gli identificativi vengono confrontati senza
    espandere automaticamente famiglie e isoforme.
    """

    path = Path(path)

    opener = (
        gzip.open
        if path.suffix.lower() == ".gz"
        else open
    )

    found = defaultdict(set)

    name = None
    parts = []

    def store():

        if name is None:
            return

        sequence = rna("".join(parts))

        if sequence:
            found[norm_id(name)].add(sequence)

    with opener(
        path,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as f:

        for line in f:

            line = line.strip()

            if line.startswith(">"):

                store()

                name = line[1:].split()[0]
                parts = []

            elif line and not line.startswith(";"):

                parts.append(line)

    store()

    if not found:
        raise ValueError(
            "FASTA miRNA maturi vuoto o non valido"
        )

    return found


# ============================================================
# 5. AUDIT DEGLI ARTEFATTI
# ============================================================

def artifact_audit(root):

    results = []

    paths = {
        "ambiente_dedicato": root / ".venv",

        "script_STEP47": (
            root
            / "scripts"
            / "step_47_fonti_b73_pipeline_unica.py"
        ),

        "requisiti_CPU": (
            root / "requirements-v33-cpu-locked.txt"
        ),
    }

    for name, path in paths.items():

        results.append({
            "componente": name,
            "presente": path.exists(),
            "percorso": str(path),
            "verificato_funzionante": False,
        })

    run46 = sorted(
        RUNS.glob("STEP46_GAE_SAFE_*")
    )

    model = (
        run46[-1] / "gae_safe_model.pt"
        if run46 else None
    )

    results.append({
        "componente": "checkpoint_GAE_STEP46",
        "presente": bool(
            model and model.is_file()
        ),
        "percorso": (
            str(model) if model else "non_trovato"
        ),
        "verificato_funzionante": False,
    })

    results.append({
        "componente": "GUI_Docker",
        "presente": "non_verificato",
        "percorso": "preservati_senza_modifiche",
        "verificato_funzionante": False,
    })

    return results


# ============================================================
# 6. TEST INTERNI
# ============================================================

def self_test():

    assert rna("ACGTu") == "ACGUU"

    assert (
        norm_id("zma-miR156a-5p")
        == norm_id("miR156a-5p")
    )

    assert (
        assembly("Zm00001d015233")
        == "RefGen_v4"
    )

    assert (
        assembly("Zm00001eb000001")
        == "NAM5"
    )

    assert (
        assembly("zma-MIR159c")
        == "MIR_precursor_locus"
    )

    assert (
        assembly("AC206198.3_FG001")
        == "legacy_FG"
    )

    for identifier, sequence in CANDIDATES.items():

        assert rna(sequence) == sequence

        assert 18 <= len(sequence) <= 35

    print(
        "[SELF_TEST_STEP48] PASS: "
        "sequenze, identificativi, assembly"
    )


# ============================================================
# 7. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "AI PRIMETTA V3.3 - Chiusura scientifica "
            "senza training"
        )
    )

    parser.add_argument(
        "--step47-run",
        type=Path,
        help=(
            "Cartella STEP47 specifica. "
            "Default: ultimo run con QC V7."
        ),
    )

    parser.add_argument(
        "--crosswalk",
        type=Path,
        help=(
            "File TSV B73v4/B73v5 di MaizeGDB, opzionale"
        ),
    )

    parser.add_argument(
        "--mature-fasta",
        type=Path,
        help=(
            "FASTA di miRNA maturi, opzionale"
        ),
    )

    parser.add_argument(
        "--self-test",
        action="store_true",
    )

    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    # ========================================================
    # 8. IDENTIFICAZIONE RUN STEP47
    # ========================================================

    if args.step47_run:

        run = args.step47_run.resolve()

    else:

        possible = [
            path
            for path in RUNS.glob(
                "STEP47_FONTI_B73_UNICO_*"
            )
            if (
                path / "QC_STEP47_V7.json"
            ).is_file()
        ]

        if not possible:
            raise SystemExit(
                "[STOP] Nessuno STEP47 con QC V7 trovato"
            )

        run = sorted(possible)[-1]

    filenames = (
        "RAPPORTO_STEP47.json",
        "QC_STEP47_V6.json",
        "QC_STEP47_V7.json",
        "ASSOCIAZIONI_EVIDENZE_DA_VERIFICARE.tsv",
    )

    paths = {
        name: run / name
        for name in filenames
    }

    for name, path in paths.items():

        if not path.is_file():
            raise SystemExit(
                f"[STOP] File STEP47 mancante: {path}"
            )

    # ========================================================
    # 9. CARICAMENTO DATI E QC
    # ========================================================

    data = read_tsv(
        paths[
            "ASSOCIAZIONI_EVIDENZE_DA_VERIFICARE.tsv"
        ]
    )

    qc6 = json.loads(
        paths["QC_STEP47_V6.json"].read_text(
            encoding="utf-8"
        )
    )

    qc7 = json.loads(
        paths["QC_STEP47_V7.json"].read_text(
            encoding="utf-8"
        )
    )

    report47 = json.loads(
        paths["RAPPORTO_STEP47.json"].read_text(
            encoding="utf-8"
        )
    )

    s4 = [
        x for x in data
        if x.get("evidence_type")
        == "PARE_cleavage_reported_in_B73_non_drought"
    ]

    s5 = [
        x for x in data
        if x.get("evidence_type")
        == "PARE_cleavage_efficiency_S5"
    ]

    s7 = [
        x for x in data
        if x.get("source") == "DROUGHT_B73_2020"
    ]

    # I conteggi sono i risultati reali dello STEP47 V7.

    expected = {
        "S4": 1583,
        "S5": 1049,
        "S7": 690,
        "TOTAL": 3322,
    }

    actual = {
        "S4": len(s4),
        "S5": len(s5),
        "S7": len(s7),
        "TOTAL": len(data),
    }

    if actual != expected:

        raise SystemExit(
            "[STOP] Numerosità STEP47 differente "
            f"dalla V7 verificata: {actual}"
        )

    if (
        qc6.get("excluded_total") != 0
        or qc7.get(
            "previously_excluded_after_correction"
        ) != 0
    ):

        raise SystemExit(
            "[STOP] QC STEP47 non concluso: "
            "esclusioni presenti"
        )

    unique_s7 = int(
        qc6.get(
            "s7_unique_sequence_gene_transcript_site_keys",
            -1,
        )
    )

    if unique_s7 != 473:

        raise SystemExit(
            "[STOP] Deduplicazione S7 inattesa: "
            f"{unique_s7}"
        )

    # ========================================================
    # 10. CONFRONTO NORMALIZZATO S4 / S5
    # ========================================================

    def pare_key(item):

        return (
            norm_id(
                item.get("tissue_or_stage", "")
            ),

            norm_id(
                item.get("pare_site_id", "")
            ),

            norm_id(
                item.get("mirna_id", "")
            ),

            norm_id(
                item.get("target_gene_id", "")
            ),
        )

    s4_keys = {
        pare_key(item)
        for item in s4
    }

    overlap_s5 = sum(
        pare_key(item) in s4_keys
        for item in s5
    )

    if overlap_s5 != 941:

        raise SystemExit(
            "[STOP] Corrispondenze S4/S5 inattese: "
            f"{overlap_s5} (attese 941)"
        )

    # ========================================================
    # 11. FILE GENOMICI OPZIONALI
    # ========================================================

    crosswalk = {}

    xref_details = {
        "usato": False,
        "origine": XREF_URL,
    }

    if args.crosswalk:

        xref_path = args.crosswalk.resolve()

        if not xref_path.is_file():
            raise SystemExit(
                "[STOP] Crosswalk non trovato: "
                f"{xref_path}"
            )

        crosswalk, details = load_xref(
            xref_path
        )

        xref_details.update({
            "usato": True,
            "file": str(xref_path),
            "sha256": digest(xref_path),
            **details,
        })

    mature = {}

    mature_details = {
        "usato": False,
        "origine": MIRBASE_URL,
    }

    if args.mature_fasta:

        mature_path = args.mature_fasta.resolve()

        if not mature_path.is_file():
            raise SystemExit(
                "[STOP] FASTA non trovato: "
                f"{mature_path}"
            )

        mature = load_mature(mature_path)

        mature_details.update({
            "usato": True,
            "file": str(mature_path),
            "sha256": digest(mature_path),
            "id_distinti": len(mature),
        })

    # ========================================================
    # 12. CREAZIONE CARTELLA STEP48
    # ========================================================

    stamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP48_CHIUSURA_PRIMETTA_{stamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False,
    )

    # ========================================================
    # 13. ARMONIZZAZIONE IDENTIFICATIVI GENOMICI
    # ========================================================

    target_ids = sorted({
        item.get("target_gene_id", "")
        for item in data
        if item.get("target_gene_id")
    })

    gene_rows = []
    status_counts = Counter()

    for gene in target_ids:

        cls = assembly(gene)

        mapped = (
            sorted(crosswalk.get(gene, set()))
            if cls == "RefGen_v4"
            else []
        )

        if cls == "NAM5":

            status = "already_NAM5"

        elif cls != "RefGen_v4":

            status = "other_assembly_or_entity"

        elif not xref_details["usato"]:

            status = "xref_not_supplied"

        elif not mapped:

            status = "no_unambiguous_xref"

        elif len(mapped) == 1:

            status = (
                "one_documented_xref_"
                "not_functional_validation"
            )

        else:

            status = "multiple_xrefs_ambiguous"

        status_counts[status] += 1

        gene_rows.append({
            "original_id": gene,
            "original_class": cls,
            "status": status,
            "n_xrefs": len(mapped),
            "nam5_xrefs": ";".join(mapped),
            "note": (
                "Associazione tra annotazioni; "
                "non validazione biologica"
            ),
        })

    write_tsv(
        output / "STEP48_CROSSWALK.tsv",
        gene_rows,
        [
            "original_id",
            "original_class",
            "status",
            "n_xrefs",
            "nam5_xrefs",
            "note",
        ],
    )

    # ========================================================
    # 14. VERIFICA IDENTITÀ miRNA MATURI PARE
    # ========================================================

    pare_rows = []

    pare_candidate_hits = Counter()
    pare_status = Counter()
    category_counts = Counter()

    for item in s4:

        mirna_id = item.get(
            "mirna_id",
            "",
        )

        sequences = sorted(
            mature.get(
                norm_id(mirna_id),
                set(),
            )
        )

        if not mature:

            status = "mature_reference_not_supplied"
            resolved = ""

        elif len(sequences) == 0:

            status = "no_exact_mature_id"
            resolved = ""

        elif len(sequences) > 1:

            status = "ambiguous_mature_sequences"
            resolved = ""

        else:

            status = "exact_id_single_mature_sequence"
            resolved = sequences[0]

        hits = [
            candidate_id
            for candidate_id, candidate_seq
            in CANDIDATES.items()
            if resolved
            and resolved == candidate_seq
        ]

        for candidate_id in hits:
            pare_candidate_hits[candidate_id] += 1

        pare_status[status] += 1

        category = str(
            item.get("pare_category", "")
        ).strip()

        category_counts[
            category or "missing"
        ] += 1

        pare_rows.append({
            "sheet": item.get("sheet", ""),
            "row": item.get("row", ""),
            "mirna_id": mirna_id,
            "target_original": item.get(
                "target_gene_id", ""
            ),
            "tissue": item.get(
                "tissue_or_stage", ""
            ),
            "pare_category": category,
            "pare_pvalue_raw": item.get(
                "pare_pvalue_raw", ""
            ),
            "mature_lookup_status": status,
            "mature_sequence_RNA": resolved,
            "candidate_exact_mature_sequence":
                ";".join(hits),
            "interpretation": (
                "Test di identità di sequenza, "
                "non validazione della siccità"
            ),
        })

    write_tsv(
        output / "STEP48_PARE_MATURE_IDENTITIES.tsv",
        pare_rows,
        [
            "sheet",
            "row",
            "mirna_id",
            "target_original",
            "tissue",
            "pare_category",
            "pare_pvalue_raw",
            "mature_lookup_status",
            "mature_sequence_RNA",
            "candidate_exact_mature_sequence",
            "interpretation",
        ],
    )

    # ========================================================
    # 15. CONFRONTO DELLE CANDIDATE CON S7
    # ========================================================

    s7_hits = Counter()

    for item in s7:

        sequence = rna(
            item.get("mirna_sequence", "")
        )

        for candidate_id, candidate_seq in CANDIDATES.items():

            if sequence and sequence == candidate_seq:
                s7_hits[candidate_id] += 1

    db_sequences = {
        sequence
        for sequences in mature.values()
        for sequence in sequences
    }

    candidate_results = []

    for candidate_id, sequence in CANDIDATES.items():

        candidate_results.append({
            "candidate": candidate_id,

            "sequence_RNA": sequence,

            "length_nt": len(sequence),

            "S7_exact_sequence_rows":
                s7_hits[candidate_id],

            "PARE_S4_exact_mature_ID_and_sequence_rows":
                pare_candidate_hits[candidate_id],

            "miRBase_exact_sequence_present": (
                "YES" if sequence in db_sequences else "NO"
            ) if mature else "NOT_CHECKED",

            "PARE_sequence_screen_status": (
                "CHECKED_WITH_MATURE_REFERENCE"
                if mature
                else "NOT_EVALUABLE_NO_MATURE_REFERENCE"
            ),

            "functional_target_validation":
                "NOT_DEMONSTRATED",
        })

    write_tsv(
        output / "STEP48_CANDIDATE_VERIFICATION.tsv",
        candidate_results,
        list(candidate_results[0]),
    )

    # ========================================================
    # 16. AUDIT DEGLI ARTEFATTI
    # ========================================================

    assets = artifact_audit(ROOT)

    write_tsv(
        output / "STEP48_ARTIFACT_AUDIT.tsv",
        assets,
        [
            "componente",
            "presente",
            "percorso",
            "verificato_funzionante",
        ],
    )

    # ========================================================
    # 17. RIEPILOGO SCIENTIFICO JSON
    # ========================================================

    overview = {
        "project": "AI PRIMETTA V3.3",

        "stage": "STEP48_END_OF_ENGINEERING_PHASE",

        "input_step47_run": str(run),

        "step47_input_sha256": {
            name: digest(path)
            for name, path in paths.items()
        },

        "evidence_counts": {
            "PARE_S4": len(s4),
            "PARE_S5": len(s5),
            "S7_rows": len(s7),
            "S7_unique": unique_s7,
            "S5_rows_matching_S4": overlap_s5,
            "excluded": 0,
        },

        "target_classes": dict(
            Counter(
                item.get(
                    "target_gene_build",
                    "unknown",
                )
                for item in data
            )
        ),

        "PARE_S4_category_counts": dict(
            category_counts
        ),

        "PARE_mature_status": dict(
            pare_status
        ),

        "xref": xref_details,

        "mature_sequences": mature_details,

        "crosswalk_status_counts": dict(
            status_counts
        ),

        "candidate_results": candidate_results,

        "model_assessment": {
            "benchmark": (
                "STEP46 transductive GAE "
                "similarity-link reconstruction"
            ),

            "reported_test_auc_from_prior_log":
                0.542509,

            "reported_test_ap_from_prior_log":
                0.554925,

            "metrics_recomputed_in_STEP48":
                False,

            "functional_target_prediction_validated":
                False,

            "external_target_benchmark_performed":
                False,
        },

        "quality_gates": {
            "STEP47_structural_QC": "PASS",

            "STEP48_archival_report": "PASS",

            "GUI_Docker_execution_test":
                "NOT_PERFORMED_IN_THIS_STEP",

            "NAM5_mapping_complete": False,

            "candidate_functional_validation":
                False,

            "new_GNN_training": False,

            "validated_target_prediction_claim":
                "NOT_SUPPORTED",
        },

        "upstream_doi": [
            "10.1111/tpj.15943",
            "10.1534/g3.120.401199",
        ],

        "source_report_present": bool(report47),

        "QC_V7_present": bool(qc7),
    }

    write_json(
        output / "STEP48_FINAL_SUMMARY.json",
        overview,
    )

    # ========================================================
    # 18. RAPPORTO CONCLUSIVO SCIENTIFICO
    # ========================================================

    lines = [
        "# AI PRIMETTA V3.3 — Rapporto conclusivo STEP 48",

        "",

        "**Stato:** conclusa la preparazione documentale "
        "e il controllo strutturale delle evidenze.",

        "",

        "**Validazione biologica dei target:** "
        "NON DIMOSTRATA.",

        "",

        "## 1. Domanda scientifica",

        "",

        "Il progetto analizza small RNA di mais B73/H21 "
        "in condizioni di controllo e siccità, "
        "con riferimento a PRJNA419326.",

        "",

        "Candidate principali: "
        "sRNA_00000037 e sRNA_00000435.",

        "",

        "## 2. Evidenze esterne",

        "",

        f"- PARE S4: {len(s4)} righe.",

        f"- PARE S5: {len(s5)} righe.",

        f"- S5 con chiave in S4: {overlap_s5}.",

        f"- Associazioni siccità S7: {len(s7)} righe.",

        f"- Associazioni S7 distinte: {unique_s7}.",

        "- Record esclusi dopo STEP47 V7: zero.",

        "",

        "Le associazioni S7 non sono prove sperimentali "
        "di taglio dei trascritti.",

        "",

        "## 3. Verifica delle candidate",

        "",
    ]

    for result in candidate_results:

        lines.append(
            f"- **{result['candidate']}**: "
            f"{result['S7_exact_sequence_rows']} match S7; "
            f"{result['PARE_S4_exact_mature_ID_and_sequence_rows']} "
            "match PARE con identificativo e sequenza matura. "
            f"Stato: {result['PARE_sequence_screen_status']}."
        )

    lines.extend([
        "",

        "Nessuno di questi confronti costituisce "
        "validazione funzionale in PRJNA419326.",

        "",

        "## 4. Armonizzazione genomica",

        "",

        "- Crosswalk v4/v5 fornito: "
        + ("sì." if crosswalk else "no."),

        "- FASTA miRNA maturi fornito: "
        + ("sì." if mature else "no."),

        "- Le associazioni genomiche non dimostrano "
        "automaticamente equivalenza funzionale.",

        "- Le famiglie miRNA composite non vengono "
        "espanse in isoforme presunte.",

        "",

        "## 5. Valutazione della GNN",

        "",

        "- STEP46: AUC test riportata 0,542509.",

        "- STEP46: AP test riportata 0,554925.",

        "- Valori riferiti alla ricostruzione di "
        "legami di similarità del grafo.",

        "- Metriche non ricalcolate nello STEP48.",

        "- Non costituiscono prova di capacità "
        "predittiva sui target biologici.",

        "",

        "## 6. Chiusura del progetto",

        "",

        "**Completato:** acquisizione e QC delle "
        "evidenze, pipeline documentale e rapporto.",

        "",

        "**Da verificare separatamente:** "
        "operatività GUI/Docker e validazione "
        "biologica sperimentale.",

        "",

        "**Stato finale:** prototipo scientifico "
        "esplorativo documentato; non validato "
        "come predittore funzionale.",

        "",

        "## 7. Fonti",

        "",

        "- DOI 10.1111/tpj.15943",

        "- DOI 10.1534/g3.120.401199",

        f"- MaizeGDB: {XREF_URL}",

        f"- miRBase: {MIRBASE_URL}",

        "",

        "## 8. Integrità",

        "",

        "Le impronte SHA-256 degli input e degli "
        "output sono riportate nel manifest STEP48.",
    ])

    report_path = (
        output / "STEP48_RAPPORTO_CONCLUSIVO.md"
    )

    report_path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    # ========================================================
    # 19. MANIFEST E HASH SHA-256
    # ========================================================

    produced = sorted(
        path
        for path in output.iterdir()
        if path.is_file()
    )

    script_path = Path(__file__).resolve()

    references = {}

    for ref in (
        args.crosswalk,
        args.mature_fasta,
    ):

        if ref is not None:
            resolved = ref.resolve()
            references[str(resolved)] = digest(resolved)

    manifest = {
        "created_utc": datetime.now(
            timezone.utc
        ).isoformat(),

        "script": str(script_path),

        "script_sha256": digest(script_path),

        "step47_run": str(run),

        "inputs": {
            str(path): digest(path)
            for path in paths.values()
        },

        "additional_references": references,

        "outputs": {
            path.name: digest(path)
            for path in produced
        },

        "immutable_inputs": True,

        "models_retrained": False,
    }

    write_json(
        output / "STEP48_MANIFEST.json",
        manifest,
    )

    # ========================================================
    # 20. OUTPUT TERMINALE
    # ========================================================

    print("=" * 72)
    print("AI PRIMETTA V3.3 | STEP48 - CHIUSURA SCIENTIFICA")
    print("=" * 72)

    print(
        "[STEP47_QC] PASS | "
        "S4=1583 S5=1049 S7=690 | esclusi=0"
    )

    print(
        "[S5_S4_MATCH_NORMALIZZATO]",
        overlap_s5,
    )

    print(
        "[S7_UNIVOCI]",
        unique_s7,
    )

    print(
        "[CROSSWALK_USATO]",
        bool(crosswalk),
    )

    print(
        "[MATURE_FASTA_USATO]",
        bool(mature),
    )

    for result in candidate_results:

        print(
            "[CANDIDATA]",
            result["candidate"],
            "| S7=",
            result["S7_exact_sequence_rows"],
            "| PARE_MATURE=",
            result[
                "PARE_S4_exact_mature_ID_and_sequence_rows"
            ],
        )

    print(
        "[VALIDAZIONE_TARGET_FUNZIONALI] "
        "NON_DIMOSTRATA"
    )

    print(
        "[GUI_DOCKER_CHECKPOINT] "
        "NON_MODIFICATI"
    )

    print(
        "[GUI_DOCKER_COLLAUDO] "
        "NON_ESEGUITO_IN_STEP48"
    )

    print("[OUTPUT]", output)
    print("[RAPPORTO]", report_path)

    print(
        "[MANIFEST]",
        output / "STEP48_MANIFEST.json",
    )

    print(
        "[OK] STEP48 completato senza training"
    )


if __name__ == "__main__":
    main()
