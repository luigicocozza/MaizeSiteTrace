
#!/usr/bin/env python3
# ============================================================
# AI PRIMETTA V3.3
# SCRIPT UNICO - VISUAL STUDIO CODE / UBUNTU 24.04
#
# FUNZIONI:
# 1. INSTALLA il nuovo modulo V3.3
# 2. VERIFICA l'integrita' della GUI V3.2
# 3. TESTA il motore su dati sintetici
# 4. AVVIA il Docker del nuovo modulo
# 5. MOSTRA lo stato dei container
# 6. ARRESTA il solo modulo V3.3
#
# SCIENZA:
# - Analisi sequenze small RNA 18-30 nt
# - Ricerca esplorativa siti RNA-trascritto
# - Mapping biologico trascritto -> gene
# - Confronto espressione DS/WW
# - Ranking geni candidati
# - Consultazione associazioni storiche V3
#
# ATTENZIONE:
# Le predizioni NON sono validate biologicamente.
# La GNN originale viene conservata, non riaddestrata.
# ============================================================

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys

from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE GENERALE
# ============================================================

ADDON = "target_prediction_v33_vscode"

GUI_FILES = [
    "app/streamlit_app.py",
    "app/pages/01_Rete_neurale_GAE.py",
    "app/pages/02_Dati_ricercatore.py",
    "app/pages/03_Ranking_dati_ricercatore.py",
    "app/pages/04_Chat_dati_ricercatore.py",
]

DNA = set("ACGT")

BASES = str.maketrans(
    "ACGT",
    "TGCA"
)

QC_FIELDS = [
    "rna_id",
    "sequence",
    "length_nt",
    "entropy_bits",
    "status",
    "reason",
]

SITE_FIELDS = [
    "rna_id",
    "transcript_id",
    "gene_id",
    "site_start_1based",
    "site_end_1based",
    "mismatches",
    "identity_fraction",
    "opposite_expression_genotypes",
    "validation_status",
]

GENE_FIELDS = [
    "gene_id",
    "n_rna",
    "rna_ids",
    "n_transcripts",
    "n_sites",
    "best_mismatches",
    "opposite_expression_genotypes",
    "scientific_status",
]

SUMMARY_FIELDS = [
    "rna_id",
    "qc",
    "rna_log2fc_ds_ww_by_genotype",
    "n_new_sites",
    "n_historical_hypotheses",
]

HIST_FIELDS = [
    "rna_id",
    "gene_id",
    "historical_12mer_hits",
    "historical_score",
    "scientific_status",
]


# ============================================================
# 2. CONTROLLO CRITTOGRAFICO GUI ORIGINALE
# ============================================================

def file_sha(path: Path) -> str:

    h = hashlib.sha256()

    with path.open("rb") as stream:

        for block in iter(
            lambda: stream.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def project_path(raw: str) -> Path:

    project = Path(raw).expanduser().resolve()

    missing = [
        str(project / filename)
        for filename in GUI_FILES
        if not (project / filename).is_file()
    ]

    if missing:

        raise ValueError(
            "Apri la cartella originale V3.2 in VS Code.\n"
            "Mancano questi file:\n"
            + "\n".join(missing)
        )

    return project


def snapshot(root: Path) -> dict:

    return {
        filename: file_sha(root / filename)
        for filename in GUI_FILES
    }


def child(root: Path) -> Path:

    return root / ADDON


def new_text(path: Path, value: str) -> None:

    if path.exists():
        raise FileExistsError(
            "Non sovrascrivo il file esistente: "
            + str(path)
        )

    path.write_text(
        value,
        encoding="utf-8"
    )


# ============================================================
# 3. INSTALLAZIONE AUTOMATICA MODULO V3.3
# ============================================================

def install(root: Path) -> None:

    destination = child(root)

    if destination.exists():

        raise FileExistsError(
            f"{destination} esiste gia'. "
            "Non eseguo sovrascritture."
        )

    source = Path(__file__).resolve()

    if not source.is_file():
        raise FileNotFoundError(
            "Script sorgente non trovato."
        )

    before = snapshot(root)

    destination.mkdir(parents=True)

    try:

        for folder in (
            "references",
            "runs",
            "demo"
        ):
            (destination / folder).mkdir()

        # Lo script copia se stesso:
        # non occorre creare un secondo file Python.

        shutil.copy2(
            source,
            destination / "app_v33.py"
        )

        # Dockerfile completamente generato

        dockerfile = '''FROM python:3.11-slim

WORKDIR /app

RUN pip install --no-cache-dir streamlit==1.50.0

COPY app_v33.py /app/app_v33.py

EXPOSE 8502

CMD ["streamlit", "run", "/app/app_v33.py", "--server.address=0.0.0.0", "--server.port=8502", "--server.headless=true", "--server.maxUploadSize=250"]
'''

        new_text(
            destination / "Dockerfile",
            dockerfile
        )

        # Docker Compose completamente generato

        compose = '''name: ai_primetta_target_v33

services:

  prediction:

    build: .

    ports:
      - "127.0.0.1:8502:8502"

    volumes:
      - ./runs:/app/runs
      - ./references:/app/references:ro
      - ..:/app/original_project:ro

    restart: unless-stopped
'''

        new_text(
            destination / "compose.yaml",
            compose
        )

        new_text(
            destination / ".dockerignore",
            "runs\nreferences\ndemo\n__pycache__\n"
        )

        # Manifesto di integrita' GUI

        manifest = destination / "GUI_V32_SHA256.json"

        manifest.write_text(
            json.dumps(
                before,
                indent=2
            ),
            encoding="utf-8"
        )

        # ====================================================
        # DATI SINTETICI PER TEST AUTOMATICO
        # ====================================================

        demo = destination / "demo"

        sequence = "ATGCATGCACTGACATGCACAT"

        target = revcomp(sequence)

        new_text(
            demo / "srna.fa",
            f">RNA_TEST_01\n{sequence}\n"
        )

        new_text(
            demo / "transcripts.fa",
            f">TX_TEST_001\nGGG{target}TTT\n"
        )

        new_text(
            demo / "tx2gene.tsv",
            "transcript_id\tgene_id\n"
            "TX_TEST_001\tGENE_TEST_001\n"
        )

        new_text(
            demo / "metadata.csv",
            "sample_id,genotype,condition\n"
            "S1,B73,WW\n"
            "S2,B73,DS\n"
        )

        new_text(
            demo / "rna_expression.csv",
            "feature_id,S1,S2\n"
            "RNA_TEST_01,10,40\n"
        )

        new_text(
            demo / "gene_expression.csv",
            "feature_id,S1,S2\n"
            "GENE_TEST_001,50,10\n"
        )

        # ====================================================
        # DOCUMENTAZIONE
        # ====================================================

        guida = """AI PRIMETTA V3.3

Modulo di predizione esplorativa RNA-gene.

GUI originale V3.2:
http://localhost:8501

Nuova GUI V3.3:
http://localhost:8502

INPUT

1. FASTA small RNA 18-30 nt
2. Trascrittoma FASTA di riferimento
3. Mapping transcript_id -> gene_id
4. Espressione RNA, facoltativa
5. Espressione geni, facoltativa
6. Metadati WW/DS, facoltativi

OUTPUT

01_rna_qc.tsv
02_alignment_hypotheses.tsv
03_historical_12mer_hypotheses.tsv
04_gene_target_candidates.tsv
05_rna_candidate_summary.tsv
06_provenance_manifest.json
LEGGIMI_RISULTATI.md

LIMITI SCIENTIFICI

Il sistema genera ipotesi computazionali.
Non dimostra il targeting biologico.
Non calcola p-value o FDR.
Non addestra una nuova GNN RNA-gene.
Gli ncRNA superiori a 30 nt richiedono un metodo diverso.
"""

        new_text(
            destination / "LEGGIMI.txt",
            guida
        )

    except Exception:

        shutil.rmtree(destination)

        raise

    after = snapshot(root)

    if before != after:

        raise RuntimeError(
            "ATTENZIONE: sono cambiati file della GUI."
        )

    print()
    print("=" * 65)
    print("INSTALLAZIONE AI PRIMETTA V3.3 COMPLETATA")
    print("=" * 65)
    print()
    print("Modulo:", destination)
    print("[OK] GUI V3.2 invariata")
    print("[OK] SHA256 salvati")
    print("[OK] Dockerfile creato")
    print("[OK] Compose creato")
    print("[OK] Motore RNA-target presente")
    print("[OK] Dati di test creati")
    print()
    print("Prossimo passaggio: comando TEST")


# ============================================================
# 4. VERIFICA DELLA GUI
# ============================================================

def verify(root: Path) -> None:

    destination = child(root)

    manifest = destination / "GUI_V32_SHA256.json"

    if not manifest.is_file():

        raise FileNotFoundError(
            "Modulo V3.3 o manifesto SHA256 mancante."
        )

    expected = json.loads(
        manifest.read_text(encoding="utf-8")
    )

    current = snapshot(root)

    different = [
        filename
        for filename in expected
        if expected[filename] != current.get(filename)
    ]

    if different:

        raise RuntimeError(
            "La GUI originale risulta modificata: "
            + ", ".join(different)
        )

    required = [
        "app_v33.py",
        "Dockerfile",
        "compose.yaml",
        "references",
        "runs",
    ]

    for filename in required:

        if not (destination / filename).exists():

            raise RuntimeError(
                "Componente V3.3 mancante: "
                + filename
            )

    print("[OK] GUI V3.2 invariata.")
    print("[OK] Modulo V3.3 presente.")
    print("[OK] Verifica SHA256 completata.")


# ============================================================
# 5. GESTIONE DOCKER
# ============================================================

def run_compose(root: Path, action: str) -> None:

    verify(root)

    if not shutil.which("docker"):

        raise RuntimeError(
            "Docker non risulta installato nella VM."
        )

    compose_file = child(root) / "compose.yaml"

    command = [
        "docker",
        "compose",
        "-f",
        str(compose_file)
    ]

    if action == "avvia":

        command += [
            "up",
            "-d",
            "--build"
        ]

    elif action == "arresta":

        command += [
            "stop"
        ]

    elif action == "stato":

        command += [
            "ps"
        ]

    print()
    print("Eseguo:", " ".join(command))
    print()

    subprocess.run(
        command,
        check=True
    )

    if action == "avvia":

        print()
        print("NUOVA GUI V3.3:")
        print("http://localhost:8502")
        print()
        print("GUI ORIGINALE V3.2:")
        print("http://localhost:8501")


# ============================================================
# 6. LETTURA FASTA, CSV, TSV E GZIP
# ============================================================

def open_text(path: Path):

    if path.name.endswith(".gz"):

        return gzip.open(
            path,
            "rt",
            encoding="utf-8-sig"
        )

    return path.open(
        "r",
        encoding="utf-8-sig"
    )


def fasta(path: Path):

    name = None
    fragments = []
    identifiers = set()

    with open_text(path) as stream:

        for raw_line in stream:

            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    yield (
                        name,
                        "".join(fragments)
                        .upper()
                        .replace("U", "T")
                    )

                tokens = line[1:].split()

                name = tokens[0] if tokens else ""

                if not name or name in identifiers:

                    raise ValueError(
                        "Identificativo FASTA vuoto/duplicato: "
                        + repr(name)
                    )

                identifiers.add(name)

                fragments = []

            elif name is None:

                raise ValueError(
                    "FASTA non valido: manca intestazione >"
                )

            else:

                fragments.append(
                    re.sub(r"\s+", "", line)
                )

    if name is not None:

        yield (
            name,
            "".join(fragments)
            .upper()
            .replace("U", "T")
        )


def revcomp(sequence: str) -> str:

    return sequence.translate(BASES)[::-1]


def sequence_entropy(sequence: str) -> float:

    length = len(sequence)

    if length == 0:
        return 0.0

    entropy = 0.0

    for base in DNA:

        count = sequence.count(base)

        if count:

            probability = count / length

            entropy -= probability * math.log2(
                probability
            )

    return entropy


def table(path: Path) -> list[dict]:

    delimiter = (
        "\t"
        if path.name.endswith((".tsv", ".tsv.gz"))
        else ","
    )

    with open_text(path) as stream:

        reader = csv.DictReader(
            stream,
            delimiter=delimiter
        )

        if not reader.fieldnames:

            raise ValueError(
                "Tabella senza intestazione: "
                + str(path)
            )

        if len(set(reader.fieldnames)) != len(
            reader.fieldnames
        ):

            raise ValueError(
                "Colonne duplicate: "
                + str(path)
            )

        return list(reader)


def require(
    rows: list,
    columns: list[str],
    path: Path
) -> None:

    if not rows:

        raise ValueError(
            "Tabella vuota: "
            + str(path)
        )

    for column in columns:

        if column not in rows[0]:

            raise ValueError(
                f"Colonna {column} mancante in {path}"
            )


def tsv(
    path: Path,
    rows: list[dict],
    columns: list[str]
) -> None:

    with path.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as stream:

        writer = csv.DictWriter(
            stream,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()

        writer.writerows(rows)


# ============================================================
# 7. ESPRESSIONE GENICA - CONFRONTO DS / WW
# ============================================================

def expressions(
    matrix: Path | None,
    metadata_path: Path | None
) -> dict:

    if matrix is None:
        return {}

    if metadata_path is None:

        raise ValueError(
            "La matrice di espressione richiede i metadati."
        )

    metadata = table(metadata_path)

    require(
        metadata,
        ["sample_id", "genotype", "condition"],
        metadata_path
    )

    samples = [
        row["sample_id"]
        for row in metadata
    ]

    if len(samples) != len(set(samples)):

        raise ValueError(
            "Identificativi campioni duplicati."
        )

    for row in metadata:

        if row["condition"] not in ("DS", "WW"):

            raise ValueError(
                "Condizione ammessa: DS oppure WW."
            )

    rows = table(matrix)

    require(
        rows,
        ["feature_id"] + samples,
        matrix
    )

    results = {}
    seen_features = set()

    for row in rows:

        feature = row["feature_id"]

        if not feature or feature in seen_features:

            raise ValueError(
                "Feature duplicata o vuota: "
                + feature
            )

        seen_features.add(feature)

        values = {}

        for sample in samples:

            value = float(row[sample])

            if not math.isfinite(value) or value < 0:

                raise ValueError(
                    "Espressione non valida: "
                    + feature
                )

            values[sample] = value

        genotype_results = {}

        genotypes = sorted({
            row["genotype"]
            for row in metadata
        })

        for genotype in genotypes:

            ds_samples = [
                values[row["sample_id"]]
                for row in metadata
                if row["genotype"] == genotype
                and row["condition"] == "DS"
            ]

            ww_samples = [
                values[row["sample_id"]]
                for row in metadata
                if row["genotype"] == genotype
                and row["condition"] == "WW"
            ]

            if not ds_samples or not ww_samples:
                continue

            mean_ds = (
                sum(ds_samples) / len(ds_samples)
            )

            mean_ww = (
                sum(ww_samples) / len(ww_samples)
            )

            log2fc = math.log2(
                (mean_ds + 1.0) /
                (mean_ww + 1.0)
            )

            genotype_results[genotype] = round(
                log2fc,
                5
            )

        results[feature] = genotype_results

    return results


# ============================================================
# 8. RICERCA SITI TARGET RNA -> TRASCRITTO
# ============================================================

def scan(
    rnas: dict,
    transcriptome: Path,
    mapping: dict,
    max_mismatches: int,
    maximum_sites: int
):

    # Il filtro con 3 ancore da 8 nucleotidi
    # non garantisce sensibilita' completa.

    seed_index = defaultdict(set)

    reverse_sequences = {
        rna_id: revcomp(sequence)
        for rna_id, sequence in rnas.items()
    }

    for rna_id, sequence in reverse_sequences.items():

        offsets = {
            0,
            (len(sequence) - 8) // 2,
            len(sequence) - 8
        }

        for offset in offsets:

            anchor = sequence[
                offset:offset + 8
            ]

            seed_index[anchor].add(
                (rna_id, offset)
            )

    hits_by_rna = defaultdict(list)

    transcripts_scanned = 0
    transcripts_unmapped = 0

    for transcript_id, sequence in fasta(
        transcriptome
    ):

        transcripts_scanned += 1

        gene_id = mapping.get(transcript_id)

        if gene_id is None:

            transcripts_unmapped += 1

            continue

        visited_sites = set()

        for index in range(
            max(0, len(sequence) - 7)
        ):

            anchor = sequence[
                index:index + 8
            ]

            if anchor not in seed_index:
                continue

            for rna_id, offset in seed_index[anchor]:

                start = index - offset

                target_sequence = (
                    reverse_sequences[rna_id]
                )

                target_length = len(
                    target_sequence
                )

                if start < 0:
                    continue

                if start + target_length > len(sequence):
                    continue

                site_identifier = (
                    rna_id,
                    start
                )

                if site_identifier in visited_sites:
                    continue

                visited_sites.add(site_identifier)

                target_site = sequence[
                    start:start + target_length
                ]

                if set(target_site) - DNA:
                    continue

                mismatches = sum(
                    base_a != base_b
                    for base_a, base_b in zip(
                        target_sequence,
                        target_site
                    )
                )

                if mismatches > max_mismatches:
                    continue

                identity = (
                    target_length - mismatches
                ) / target_length

                hits_by_rna[rna_id].append({

                    "rna_id": rna_id,

                    "transcript_id": transcript_id,

                    "gene_id": gene_id,

                    "site_start_1based": start + 1,

                    "site_end_1based":
                        start + target_length,

                    "mismatches": mismatches,

                    "identity_fraction":
                        round(identity, 6),

                    "validation_status":
                        "EXPLORATORY_NOT_VALIDATED"
                })

            # Protezione dalla crescita incontrollata RAM

            total_hits = sum(
                len(rows)
                for rows in hits_by_rna.values()
            )

            if total_hits > 100000:

                for rna_id in hits_by_rna:

                    hits_by_rna[rna_id].sort(
                        key=lambda row: (
                            row["mismatches"],
                            row["transcript_id"],
                            row["site_start_1based"]
                        )
                    )

                    del hits_by_rna[rna_id][
                        maximum_sites:
                    ]

    result = []

    for rna_id in sorted(hits_by_rna):

        rows = sorted(
            hits_by_rna[rna_id],
            key=lambda row: (
                row["mismatches"],
                row["transcript_id"],
                row["site_start_1based"]
            )
        )

        result.extend(
            rows[:maximum_sites]
        )

    statistics = {

        "transcripts_scanned":
            transcripts_scanned,

        "transcripts_without_tx2gene":
            transcripts_unmapped,

        "seed":
            "3x8mer_NON_EXHAUSTIVE",

        "max_sites_per_rna":
            maximum_sites,
    }

    return result, statistics


# ============================================================
# 9. PIPELINE COMPLETA DI ANALISI
# ============================================================

def analyze(
    rna_fasta: Path,
    transcripts: Path | None,
    mapping_file: Path | None,
    rna_expression: Path | None,
    gene_expression: Path | None,
    metadata: Path | None,
    historical_file: Path | None,
    output: Path,
    max_mismatches: int = 2,
    entropy_min: float = 1.4,
    maximum_sites: int = 300
) -> dict:

    if (
        (transcripts is None)
        !=
        (mapping_file is None)
    ):

        raise ValueError(
            "Trascrittoma e tx2gene devono essere "
            "presenti insieme."
        )

    if transcripts is None and historical_file is None:

        raise ValueError(
            "Per ottenere target servono trascrittoma "
            "e mapping oppure l'archivio storico."
        )

    if (
        rna_expression or gene_expression
    ) and metadata is None:

        raise ValueError(
            "Mancano i metadati dei campioni."
        )

    if not 0 <= max_mismatches <= 3:

        raise ValueError(
            "Mismatch ammessi: 0-3."
        )

    if not 0 <= entropy_min <= 2:

        raise ValueError(
            "Entropia fuori intervallo."
        )

    if not 1 <= maximum_sites <= 3000:

        raise ValueError(
            "Numero massimo siti non valido."
        )

    named_inputs = {

        "rna_fasta": rna_fasta,

        "transcriptome": transcripts,

        "tx2gene": mapping_file,

        "rna_expression": rna_expression,

        "gene_expression": gene_expression,

        "metadata": metadata,

        "historical": historical_file,
    }

    inputs = [
        path
        for path in named_inputs.values()
        if path is not None
    ]

    for path in inputs:

        if not path.is_file():

            raise FileNotFoundError(
                "File di input mancante: "
                + str(path)
            )

    # ========================================================
    # QC SEQUENZE RNA
    # ========================================================

    qualified = {}

    quality = []

    for rna_id, sequence in fasta(rna_fasta):

        reasons = []

        if set(sequence) - DNA:

            reasons.append(
                "basi_ambigue"
            )

        if not 18 <= len(sequence) <= 30:

            reasons.append(
                "fuori_18_30_nt"
            )

        entropy = round(
            sequence_entropy(sequence),
            4
        )

        if entropy < entropy_min:

            reasons.append(
                "bassa_complessita"
            )

        status = (
            "eligible"
            if not reasons
            else "excluded"
        )

        quality.append({

            "rna_id": rna_id,

            "sequence": sequence,

            "length_nt": len(sequence),

            "entropy_bits": entropy,

            "status": status,

            "reason": ";".join(reasons)
        })

        if not reasons:

            qualified[rna_id] = sequence

    if not quality:

        raise ValueError(
            "Nessuna sequenza RNA nel FASTA."
        )

    # ========================================================
    # ESPRESSIONE
    # ========================================================

    de_rna = expressions(
        rna_expression,
        metadata
    )

    de_gene = expressions(
        gene_expression,
        metadata
    )

    # ========================================================
    # MAPPING TRASCRITTO -> GENE
    # ========================================================

    sites = []

    scan_statistics = {
        "run": False
    }

    if transcripts is not None:

        mapping_rows = table(
            mapping_file
        )

        require(
            mapping_rows,
            ["transcript_id", "gene_id"],
            mapping_file
        )

        mapping = {}

        for row in mapping_rows:

            transcript_id = row[
                "transcript_id"
            ].strip()

            gene_id = row[
                "gene_id"
            ].strip()

            if not transcript_id or not gene_id:

                raise ValueError(
                    "Mapping tx2gene incompleto."
                )

            if transcript_id in mapping:

                if mapping[transcript_id] != gene_id:

                    raise ValueError(
                        "Mapping tx2gene conflittuale."
                    )

            mapping[transcript_id] = gene_id

        sites, scan_statistics = scan(
            qualified,
            transcripts,
            mapping,
            max_mismatches,
            maximum_sites
        )

        scan_statistics["run"] = True

    # ========================================================
    # ASSOCIAZIONE RNA-GENE CON ESPRESSIONE
    # ========================================================

    for site in sites:

        rna_id = site["rna_id"]

        gene_id = site["gene_id"]

        related_genotypes = (
            set(de_rna.get(rna_id, {}))
            &
            set(de_gene.get(gene_id, {}))
        )

        opposite = []

        for genotype in related_genotypes:

            rna_fc = de_rna[
                rna_id
            ][genotype]

            gene_fc = de_gene[
                gene_id
            ][genotype]

            if (
                abs(rna_fc) >= 1
                and abs(gene_fc) >= 1
                and rna_fc * gene_fc < 0
            ):

                opposite.append(genotype)

        site[
            "opposite_expression_genotypes"
        ] = ";".join(sorted(opposite))

    # ========================================================
    # ARCHIVIO STORICO V3
    # ========================================================

    historical = []

    if historical_file is not None:

        previous = table(
            historical_file
        )

        require(
            previous,
            ["srna_sequence", "gene_id"],
            historical_file
        )

        sequence_to_rna = defaultdict(list)

        for rna_id, sequence in qualified.items():

            sequence_to_rna[sequence].append(
                rna_id
            )

        for row in previous:

            sequence = re.sub(
                r"\s+",
                "",
                row["srna_sequence"]
            ).upper().replace("U", "T")

            for rna_id in sequence_to_rna.get(
                sequence,
                []
            ):

                historical.append({

                    "rna_id": rna_id,

                    "gene_id": row["gene_id"],

                    "historical_12mer_hits":
                        row.get(
                            "n_matched_12mers",
                            ""
                        ),

                    "historical_score":
                        row.get(
                            "target_support_score",
                            ""
                        ),

                    "scientific_status":
                        "HISTORICAL_12MER_NOT_VALIDATED"
                })

    # ========================================================
    # RANKING GENI CANDIDATI
    # ========================================================

    grouped = defaultdict(
        lambda: {

            "rnas": set(),

            "transcripts": set(),

            "sites": 0,

            "best_mm": 999,

            "opposite": set()
        }
    )

    for site in sites:

        gene_id = site["gene_id"]

        info = grouped[gene_id]

        info["rnas"].add(
            site["rna_id"]
        )

        info["transcripts"].add(
            site["transcript_id"]
        )

        info["sites"] += 1

        info["best_mm"] = min(
            info["best_mm"],
            site["mismatches"]
        )

        opposite = site[
            "opposite_expression_genotypes"
        ]

        if opposite:

            info["opposite"].update(
                opposite.split(";")
            )

    genes = []

    for gene_id, info in grouped.items():

        genes.append({

            "gene_id": gene_id,

            "n_rna": len(info["rnas"]),

            "rna_ids":
                ";".join(sorted(info["rnas"])),

            "n_transcripts":
                len(info["transcripts"]),

            "n_sites": info["sites"],

            "best_mismatches":
                info["best_mm"],

            "opposite_expression_genotypes":
                ";".join(
                    sorted(info["opposite"])
                ),

            "scientific_status":
                "EXPLORATORY_NOT_VALIDATED"
        })

    genes.sort(
        key=lambda row: (

            row["best_mismatches"],

            -row["n_rna"],

            -row["n_sites"],

            row["gene_id"]
        )
    )

    # ========================================================
    # RIEPILOGO RNA
    # ========================================================

    new_counts = defaultdict(int)

    historical_counts = defaultdict(int)

    for site in sites:

        new_counts[
            site["rna_id"]
        ] += 1

    for item in historical:

        historical_counts[
            item["rna_id"]
        ] += 1

    summary = []

    for item in quality:

        rna_id = item["rna_id"]

        fc = ";".join(

            f"{genotype}:{value:+.3f}"

            for genotype, value in sorted(
                de_rna.get(
                    rna_id,
                    {}
                ).items()
            )
        )

        summary.append({

            "rna_id": rna_id,

            "qc": item["status"],

            "rna_log2fc_ds_ww_by_genotype": fc,

            "n_new_sites":
                new_counts[rna_id],

            "n_historical_hypotheses":
                historical_counts[rna_id]
        })

    # ========================================================
    # SALVATAGGIO RISULTATI
    # ========================================================

    if output.exists():

        raise FileExistsError(
            "La cartella output esiste gia': "
            + str(output)
        )

    output.mkdir(
        parents=True
    )

    tsv(
        output / "01_rna_qc.tsv",
        quality,
        QC_FIELDS
    )

    tsv(
        output / "02_alignment_hypotheses.tsv",
        sites,
        SITE_FIELDS
    )

    tsv(
        output / "03_historical_12mer_hypotheses.tsv",
        historical,
        HIST_FIELDS
    )

    tsv(
        output / "04_gene_target_candidates.tsv",
        genes,
        GENE_FIELDS
    )

    tsv(
        output / "05_rna_candidate_summary.tsv",
        summary,
        SUMMARY_FIELDS
    )

    # ========================================================
    # TRACCIABILITA SCIENTIFICA
    # ========================================================

    manifest = {

        "project":
            "AI PRIMETTA V3.3",

        "created_utc":
            datetime.now(
                timezone.utc
            ).isoformat(),

        "status":
            "EXPLORATORY_NOT_BIOLOGICALLY_VALIDATED",

        "n_input_rna":
            len(quality),

        "n_eligible_rna":
            len(qualified),

        "n_alignment_sites":
            len(sites),

        "n_genes":
            len(genes),

        "n_historical_hypotheses":
            len(historical),

        "input_sha256": {

            name: file_sha(path)

            for name, path in named_inputs.items()

            if path is not None
        },

        "scan":
            scan_statistics,

        "method":
            "ungapped_reverse_complement_3x8mer",

        "limits": [

            "Non esaustivo",

            "Non calcola energia di legame",

            "Non considera bulge o accessibilita",

            "Nessun p-value o FDR",

            "Nessuna GNN RNA-gene riaddestrata",

            "Identita di sequenza non dimostra targeting",

            "lncRNA >30 nt esclusi"
        ]
    }

    manifest_path = (
        output /
        "06_provenance_manifest.json"
    )

    manifest_path.write_text(
        json.dumps(
            manifest,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    report = f"""# AI PRIMETTA V3.3

RISULTATI COMPUTAZIONALI ESPLORATIVI.

NON VALIDATI BIOLOGICAMENTE.

RNA inseriti: {len(quality)}
RNA idonei: {len(qualified)}
Siti candidati: {len(sites)}
Geni candidati: {len(genes)}
Associazioni storiche: {len(historical)}

Il motore produce ipotesi RNA-gene.

Le associazioni storiche 12-mer sono separate.

Il modello GNN/GAE originale non viene riaddestrato.

Per la validazione occorrono evidenze indipendenti:
psRNATarget, degradoma/PARE, RNA-seq e qPCR.
"""

    (
        output /
        "LEGGIMI_RISULTATI.md"
    ).write_text(
        report,
        encoding="utf-8"
    )

    return manifest


# ============================================================
# 10. TEST AUTOMATICI
# ============================================================

def test_offline(root: Path) -> None:

    verify(root)

    demo = child(root) / "demo"

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%d_%H%M%S_%f"
    )

    output = (
        child(root)
        / "runs"
        / f"TEST_SINTETICO_{stamp}"
    )

    result = analyze(

        rna_fasta=demo / "srna.fa",

        transcripts=demo / "transcripts.fa",

        mapping_file=demo / "tx2gene.tsv",

        rna_expression=demo / "rna_expression.csv",

        gene_expression=demo / "gene_expression.csv",

        metadata=demo / "metadata.csv",

        historical_file=None,

        output=output
    )

    assert result[
        "n_eligible_rna"
    ] == 1

    assert result[
        "n_genes"
    ] == 1

    assert result[
        "n_alignment_sites"
    ] == 1

    genes = table(
        output /
        "04_gene_target_candidates.tsv"
    )

    assert genes[0][
        "gene_id"
    ] == "GENE_TEST_001"

    assert genes[0][
        "opposite_expression_genotypes"
    ] == "B73"

    # Test negativo: senza mapping non bisogna
    # inventare identificativi dei geni.

    try:

        analyze(

            rna_fasta=demo / "srna.fa",

            transcripts=demo / "transcripts.fa",

            mapping_file=None,

            rna_expression=None,

            gene_expression=None,

            metadata=None,

            historical_file=None,

            output=output.with_name(
                output.name + "_NEG"
            )
        )

        raise AssertionError(
            "Il mapping assente e stato accettato."
        )

    except ValueError:

        pass

    verify(root)

    print()
    print("=" * 60)
    print("TEST AI PRIMETTA V3.3")
    print("=" * 60)

    print("[OK] FASTA RNA")

    print("[OK] Reverse complement")

    print("[OK] Mapping RNA-trascritto-gene")

    print("[OK] Espressione DS/WW")

    print("[OK] Ranking gene candidato")

    print("[OK] Test mapping mancante")

    print("[OK] GUI V3.2 invariata")

    print()
    print("Risultati test:", output)
    print()
    print(
        "ATTENZIONE: test sintetico, "
        "non validazione biologica."
    )


# ============================================================
# 11. INTERFACCIA GRAFICA STREAMLIT V3.3
# ============================================================

def run_streamlit() -> None:

    import streamlit as st

    st.set_page_config(

        page_title="AI PRIMETTA V3.3",

        page_icon="🌽",

        layout="wide"
    )

    st.title(
        "AI PRIMETTA V3.3"
    )

    st.subheader(
        "Predizione esplorativa small RNA - geni target"
    )

    st.info(
        "La GUI originale V3.2 non e stata modificata. "
        "Questo modulo e indipendente sulla porta 8502."
    )

    st.warning(
        "I geni individuati sono candidati computazionali. "
        "Non sono target biologicamente validati."
    )

    st.markdown(
        """
        ### Obiettivo biologico

        Identificare possibili geni del mais associati
        a sequenze small RNA, integrando:

        - Complementarita RNA-trascritto.
        - Identificativi trascritto-gene.
        - Espressione dei piccoli RNA.
        - Espressione dei geni.
        - Condizioni controllo WW e siccita DS.
        """
    )

    # ========================================================
    # INPUT BIOLOGICI
    # ========================================================

    st.header(
        "1. Dati sperimentali"
    )

    column_a, column_b = st.columns(2)

    with column_a:

        rna = st.file_uploader(

            "FASTA small RNA / ncRNA",

            type=[
                "fa",
                "fasta",
                "fna",
                "gz"
            ]
        )

        rna_expression = st.file_uploader(

            "Matrice espressione small RNA",

            type=[
                "csv",
                "tsv"
            ]
        )

        gene_expression = st.file_uploader(

            "Matrice espressione genica",

            type=[
                "csv",
                "tsv"
            ]
        )

        metadata = st.file_uploader(

            "Metadati campioni DS/WW",

            type=[
                "csv",
                "tsv"
            ]
        )

    # ========================================================
    # RIFERIMENTI GENOMICI
    # ========================================================

    with column_b:

        reference_dir = Path(
            "/app/references"
        )

        if not reference_dir.is_dir():

            reference_dir = (
                Path(__file__).parent /
                "references"
            )

        available = sorted(
            path.name
            for path in reference_dir.iterdir()
            if path.is_file()
        )

        transcript_name = st.selectbox(

            "FASTA trascrittoma riferimento",

            [""] + available
        )

        mapping_name = st.selectbox(

            "File mapping tx2gene",

            [""] + available
        )

        historical_path = Path(

            "/app/original_project/outputs/"
            "prjna419326/"
            "PRJNA419326_srna_rnaseq_"
            "predicted_target_links.tsv"
        )

        historical_enabled = st.checkbox(

            "Consulta associazioni storiche V3",

            value=False,

            disabled=not historical_path.is_file()
        )

        max_mismatches = st.slider(

            "Mismatch massimi",

            0,
            3,
            2
        )

        entropy_min = st.slider(

            "Entropia minima",

            0.0,
            2.0,
            1.4,
            0.1
        )

        maximum_sites = st.number_input(

            "Siti massimi per RNA",

            min_value=1,

            max_value=3000,

            value=300
        )

        experiment = st.text_input(

            "Nome esperimento",

            "nuovi_campioni_mais"
        )

    # ========================================================
    # GUIDA INPUT
    # ========================================================

    with st.expander(
        "Formati dei file da utilizzare"
    ):

        st.markdown(
            "**Esempio FASTA RNA**"
        )

        st.code(
            ">RNA_001\nAUGCAUGCACUGACAUGCACAU",
            language="text"
        )

        st.markdown(
            "**Esempio mapping tx2gene**"
        )

        st.code(
            "transcript_id\tgene_id\n"
            "Zm_tx_001\tZm_gene_001",
            language="text"
        )

        st.markdown(
            "**Metadati sperimentali**"
        )

        st.code(
            "sample_id,genotype,condition\n"
            "S1,B73,WW\n"
            "S2,B73,DS",
            language="text"
        )

        st.markdown(
            "**Matrice di espressione normalizzata**"
        )

        st.code(
            "feature_id,S1,S2\n"
            "RNA_001,10,30",
            language="text"
        )

        st.caption(
            "L'espressione deve essere gia normalizzata. "
            "Il confronto log2FC e descrittivo, "
            "senza p-value o FDR."
        )

    # ========================================================
    # CALCOLO
    # ========================================================

    st.header(
        "2. Analisi e predizione"
    )

    if st.button(

        "Calcola geni target candidati",

        type="primary"
    ):

        if rna is None:

            st.error(
                "Devi caricare un FASTA RNA."
            )

        elif not (
            (transcript_name and mapping_name)
            or historical_enabled
        ):

            st.error(
                "Serve trascrittoma + mapping "
                "oppure archivio storico."
            )

        elif bool(transcript_name) != bool(
            mapping_name
        ):

            st.error(
                "Trascrittoma e mapping devono "
                "essere selezionati insieme."
            )

        else:

            safe_name = re.sub(
                r"[^a-zA-Z0-9_-]",
                "_",
                experiment
            ).strip("_")[:35]

            if not safe_name:

                safe_name = "esperimento"

            stamp = datetime.now(
                timezone.utc
            ).strftime(
                "%Y%m%d_%H%M%S_%f_UTC"
            )

            run_root = Path(
                "/app/runs"
            )

            if not run_root.is_dir():

                run_root = (
                    Path(__file__).parent /
                    "runs"
                )

            upload_dir = run_root / (
                "input_" + safe_name + "_" + stamp
            )

            output_dir = run_root / (
                "result_" + safe_name + "_" + stamp
            )

            upload_dir.mkdir(
                parents=True,
                exist_ok=False
            )

            def save_upload(
                uploaded,
                basename
            ):

                if uploaded is None:
                    return None

                extension = (
                    ".tsv"
                    if uploaded.name.lower().endswith(
                        ".tsv"
                    )
                    else ".csv"
                )

                if basename == "RNA":

                    extension = (
                        ".fa.gz"
                        if uploaded.name.lower().endswith(
                            ".gz"
                        )
                        else ".fa"
                    )

                path = upload_dir / (
                    basename + extension
                )

                path.write_bytes(
                    uploaded.getvalue()
                )

                return path

            try:

                rna_path = save_upload(
                    rna,
                    "RNA"
                )

                rna_expr_path = save_upload(
                    rna_expression,
                    "rna_expression"
                )

                gene_expr_path = save_upload(
                    gene_expression,
                    "gene_expression"
                )

                metadata_path = save_upload(
                    metadata,
                    "metadata"
                )

                transcript_path = (
                    reference_dir / transcript_name
                    if transcript_name
                    else None
                )

                mapping_path = (
                    reference_dir / mapping_name
                    if mapping_name
                    else None
                )

                previous_path = (
                    historical_path
                    if historical_enabled
                    else None
                )

                with st.spinner(
                    "Elaborazione RNA-gene in corso..."
                ):

                    result = analyze(

                        rna_fasta=rna_path,

                        transcripts=transcript_path,

                        mapping_file=mapping_path,

                        rna_expression=rna_expr_path,

                        gene_expression=gene_expr_path,

                        metadata=metadata_path,

                        historical_file=previous_path,

                        output=output_dir,

                        max_mismatches=max_mismatches,

                        entropy_min=entropy_min,

                        maximum_sites=int(
                            maximum_sites
                        )
                    )

                st.session_state[
                    "v33_output"
                ] = str(output_dir)

                st.success(
                    f"Analisi completata: "
                    f"{result['n_alignment_sites']} siti, "
                    f"{result['n_genes']} geni candidati, "
                    f"{result['n_historical_hypotheses']} "
                    f"associazioni storiche."
                )

            except Exception as error:

                st.error(
                    "Analisi non completata: "
                    + str(error)
                )

    # ========================================================
    # OUTPUT E VISUALIZZAZIONE GUI
    # ========================================================

    result_dir = st.session_state.get(
        "v33_output"
    )

    if result_dir and Path(result_dir).is_dir():

        result_path = Path(
            result_dir
        )

        st.header(
            "3. Risultati biologico-computazionali"
        )

        output_files = [

            (
                "Ranking geni candidati",
                "04_gene_target_candidates.tsv"
            ),

            (
                "Siti RNA-trascritto",
                "02_alignment_hypotheses.tsv"
            ),

            (
                "Associazioni storiche",
                "03_historical_12mer_hypotheses.tsv"
            ),

            (
                "Screening RNA",
                "05_rna_candidate_summary.tsv"
            ),

            (
                "Controlli qualita",
                "01_rna_qc.tsv"
            )
        ]

        for title, filename in output_files:

            with st.expander(
                title,
                expanded=(
                    title == "Ranking geni candidati"
                )
            ):

                file_path = (
                    result_path / filename
                )

                data = table(file_path)

                st.dataframe(
                    data,
                    use_container_width=True
                )

                st.download_button(

                    "Scarica " + filename,

                    file_path.read_bytes(),

                    file_name=filename,

                    mime="text/tab-separated-values",

                    key="download_" + filename
                )

        st.download_button(

            "Scarica manifesto scientifico JSON",

            (
                result_path /
                "06_provenance_manifest.json"
            ).read_bytes(),

            file_name=(
                "06_provenance_manifest.json"
            )
        )

        st.warning(
            "La presenza di un gene in tabella "
            "non costituisce una conferma "
            "sperimentale del targeting."
        )


# ============================================================
# 12. PROGRAMMA PRINCIPALE
# ============================================================

def main() -> int:

    # Streamlit utilizza lo stesso file Python.

    if (
        os.environ.get(
            "STREAMLIT_SERVER_PORT"
        )
        or "streamlit" in sys.modules
    ):

        run_streamlit()

        return 0

    parser = argparse.ArgumentParser(

        description=(
            "AI PRIMETTA V3.3 "
            "- gestione completa Ubuntu"
        )
    )

    parser.add_argument(

        "azione",

        choices=[
            "installa",
            "verifica",
            "test",
            "avvia",
            "arresta",
            "stato"
        ]
    )

    parser.add_argument(

        "--project",

        default=".",

        help=(
            "Radice progetto V3.2 "
            "contenente app/streamlit_app.py"
        )
    )

    args = parser.parse_args()

    root = project_path(
        args.project
    )

    if args.azione == "installa":

        install(root)

    elif args.azione == "verifica":

        verify(root)

    elif args.azione == "test":

        test_offline(root)

    else:

        run_compose(
            root,
            args.azione
        )

    return 0


if __name__ == "__main__":

    try:

        sys.exit(
            main()
        )

    except (
        ValueError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        subprocess.CalledProcessError,
        AssertionError
    ) as error:

        print(
            "[ERRORE]",
            error,
            file=sys.stderr
        )

        sys.exit(1)
