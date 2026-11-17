
#!/usr/bin/env python3
"""
AI PRIMETTA V3.3 - STEP 16 CORRETTO
Test funzionale su 5 small RNA reali PRJNA419326.

- Non modifica la GUI V3.2.
- Non modifica il motore V3.3.
- Non modifica i modelli ML/GNN.
- Usa il trascrittoma reale Zea mays.
- Verifica i siti trovati contro il FASTA originale.
- Verifica la corrispondenza trascritto-gene.
- Genera rapporto JSON e Markdown.

ATTENZIONE:
Test funzionale, NON validazione biologica indipendente.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import subprocess
import sys

from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. PERCORSI DEL PROGETTO
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

ADDON = ROOT / "target_prediction_v33_vscode"

SOURCE = (
    ROOT
    / "outputs/prjna419326"
    / "PRJNA419326_experimental_srna_v3_ranking_top500.fasta"
)

STEM = "Zea_mays.Zm-B73-REFERENCE-NAM-5.0"

REF_FASTA = (
    ADDON / "references" / (STEM + ".cdna.all.fa.gz")
)

REF_MAPPING = (
    ADDON / "references" / (STEM + ".tx2gene.tsv")
)

COMPOSE = ADDON / "compose.yaml"

GUI_MANIFEST = ADDON / "GUI_V32_SHA256.json"

N_RNA = 5
MAX_MISMATCHES = 2
MIN_ENTROPY = 1.4
MAX_SITES_PER_RNA = 300


# ============================================================
# 2. HASH SHA256
# ============================================================

def sha256(path: Path) -> str:

    digest = hashlib.sha256()

    with path.open("rb") as stream:

        for block in iter(
            lambda: stream.read(1024 * 1024),
            b""
        ):
            digest.update(block)

    return digest.hexdigest()


# ============================================================
# 3. LETTURA FASTA
# ============================================================

def read_fasta(path: Path):

    identifier = None
    fragments = []
    seen = set()

    with path.open(
        "r",
        encoding="utf-8-sig"
    ) as stream:

        for lineno, raw in enumerate(stream, 1):

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                if identifier is not None:
                    yield (
                        identifier,
                        "".join(fragments)
                        .upper()
                        .replace("U", "T")
                    )

                fields = line[1:].split()

                identifier = fields[0] if fields else ""

                if not identifier or identifier in seen:
                    raise ValueError(
                        f"ID FASTA vuoto o duplicato "
                        f"alla riga {lineno}"
                    )

                seen.add(identifier)
                fragments = []

            else:

                if identifier is None:
                    raise ValueError(
                        f"FASTA senza intestazione "
                        f"alla riga {lineno}"
                    )

                fragments.append(
                    "".join(line.split())
                )

    if identifier is not None:
        yield (
            identifier,
            "".join(fragments)
            .upper()
            .replace("U", "T")
        )


# ============================================================
# 4. ENTROPIA DELLE SEQUENZE
# ============================================================

def entropy(seq: str) -> float:

    if not seq:
        return 0.0

    return -sum(
        (seq.count(base) / len(seq))
        * math.log2(seq.count(base) / len(seq))

        for base in "ACGT"

        if base in seq
    )


# ============================================================
# 5. SELEZIONE DI CINQUE SMALL RNA REALI
# ============================================================

def select_five(path: Path):

    selected = []
    excluded = []
    observed_sequences = set()

    for rank, (rid, seq) in enumerate(
        read_fasta(path),
        1
    ):

        reasons = []

        if not 18 <= len(seq) <= 30:
            reasons.append(
                "lunghezza_non_18_30nt"
            )

        if set(seq) - set("ACGT"):
            reasons.append(
                "basi_ambigue"
            )

        h = entropy(seq)

        if h < MIN_ENTROPY:
            reasons.append(
                "bassa_complessita"
            )

        if seq in observed_sequences:
            reasons.append(
                "sequenza_duplicata"
            )

        if reasons:

            excluded.append({
                "rank": rank,
                "rna_id": rid,
                "reasons": reasons
            })

            continue

        observed_sequences.add(seq)

        selected.append({
            "rank": rank,
            "rna_id": rid,
            "sequence": seq,
            "length_nt": len(seq),
            "entropy_bits": round(h, 4)
        })

        if len(selected) == N_RNA:
            break

    if len(selected) != N_RNA:
        raise ValueError(
            f"Trovati {len(selected)} small RNA idonei; "
            f"richiesti {N_RNA}"
        )

    return selected, excluded


# ============================================================
# 6. VERIFICA INTEGRITA GUI ORIGINALE
# ============================================================

def verify_gui():

    if not GUI_MANIFEST.is_file():
        raise FileNotFoundError(
            f"Manifest GUI assente: {GUI_MANIFEST}"
        )

    hashes = json.loads(
        GUI_MANIFEST.read_text(
            encoding="utf-8"
        )
    )

    if not isinstance(hashes, dict) or not hashes:
        raise ValueError(
            "Manifest SHA256 GUI non valido"
        )

    changed = []

    for name, expected in hashes.items():

        path = ROOT / name

        if (
            not path.is_file()
            or sha256(path) != expected
        ):
            changed.append(name)

    if changed:
        raise RuntimeError(
            "GUI V3.2 non conforme: "
            + ", ".join(changed)
        )

    print(
        "[OK] GUI V3.2 conforme agli SHA256",
        flush=True
    )


# ============================================================
# 7. ESECUZIONE NEL CONTAINER DOCKER
# ============================================================

def run_docker(script: str):

    command = [
        "docker",
        "compose",
        "-f",
        str(COMPOSE),
        "exec",
        "-T",
        "--user",
        f"{os.getuid()}:{os.getgid()}",
        "prediction",
        "python3",
        "-u",
        "-"
    ]

    subprocess.run(
        command,
        input=script,
        text=True,
        cwd=ROOT,
        check=True,
        timeout=1800
    )


# ============================================================
# 8. PROCEDURA PRINCIPALE
# ============================================================

def main() -> int:

    print("=" * 66)
    print(
        "AI PRIMETTA V3.3 - STEP 16 CORRETTO"
    )
    print(
        "TEST CON 5 SMALL RNA REALI PRJNA419326"
    )
    print("=" * 66)

    required = (
        SOURCE,
        REF_FASTA,
        REF_MAPPING,
        COMPOSE,
        GUI_MANIFEST,
        ADDON / "app_v33.py"
    )

    for file in required:

        if (
            not file.is_file()
            or file.stat().st_size == 0
        ):
            raise FileNotFoundError(
                f"File mancante o vuoto: {file}"
            )

    verify_gui()

    chosen, rejected = select_five(SOURCE)

    print()
    print("SMALL RNA SELEZIONATI:")

    for item in chosen:

        print(
            f"Rank {item['rank']} | "
            f"{item['rna_id']} | "
            f"{item['length_nt']} nt | "
            f"Entropia={item['entropy_bits']}"
        )

    # --------------------------------------------------------
    # Cartella unica per ogni esecuzione
    # --------------------------------------------------------

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%S_%fZ"
    )

    run = (
        ADDON
        / "runs"
        / ("STEP16_5_SRNA_REALI_" + stamp)
    )

    run.mkdir(
        parents=True,
        exist_ok=False
    )

    selected_fasta = (
        run / "input_5_srna_reali.fa"
    )

    with selected_fasta.open(
        "w",
        encoding="utf-8"
    ) as stream:

        for item in chosen:

            stream.write(
                f">{item['rna_id']}\n"
                f"{item['sequence']}\n"
            )

    provenance = {
        "source": str(SOURCE.relative_to(ROOT)),
        "sha256": sha256(SOURCE),
        "selection": (
            "Primi cinque small RNA idonei "
            "e distinti; non campione casuale"
        ),
        "selected": chosen,
        "rejected_before_selection": rejected,
        "status": (
            "TEST_FUNZIONALE_NON_VALIDAZIONE_BIOLOGICA"
        )
    }

    (
        run / "SELEZIONE_E_PROVENIENZA.json"
    ).write_text(
        json.dumps(
            provenance,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    docker_run = (
        "/app/runs/" + run.name
    )

    print()
    print(
        "[INFO] Avvio motore V3.3 su trascrittoma reale",
        flush=True
    )

    print(
        "[INFO] Nessuna modifica a GUI, modelli e riferimenti",
        flush=True
    )

    # ========================================================
    # 9. CODICE ESEGUITO DIRETTAMENTE IN DOCKER
    # ========================================================

    inside = r'''
import inspect
import json
import math
import sys

from collections import defaultdict
from pathlib import Path

sys.path.insert(0, "/app")

import app_v33 as engine

base = Path(__DOCKER_RUN__)

ref = Path("/app/references")

stem = "Zea_mays.Zm-B73-REFERENCE-NAM-5.0"

out = base / "risultati"

rna_path = base / "input_5_srna_reali.fa"

transcriptome = (
    ref / (stem + ".cdna.all.fa.gz")
)

tx2gene_path = (
    ref / (stem + ".tx2gene.tsv")
)

# ============================================================
# CONTROLLO DELLA FIRMA DEL MOTORE
# ============================================================

# Le revisioni del backend V3.3 possono usare
# nomi differenti per gli ultimi parametri.
#
# Gli argomenti posizionali evitano il problema
# dei nomi max_mm / max_mismatches e
# max_sites / maximum_sites.
#
# I primi otto argomenti indicano i file;
# gli ultimi tre indicano i parametri numerici.

args = (
    rna_path,
    transcriptome,
    tx2gene_path,
    None,
    None,
    None,
    None,
    out,
    2,
    1.4,
    300
)

signature = inspect.signature(
    engine.analyze
)

print(
    "[DOCKER] Firma del motore:",
    signature,
    flush=True
)

try:

    signature.bind(*args)

except TypeError as error:

    raise RuntimeError(
        "Firma del motore non compatibile: "
        + str(error)
    ) from error

# ============================================================
# ESECUZIONE DELLA PREDIZIONE
# ============================================================

print(
    "[DOCKER] Avvio analisi esplorativa...",
    flush=True
)

manifest = engine.analyze(*args)

print(
    "[DOCKER] Analisi completata.",
    flush=True
)

# ============================================================
# VERIFICA RISULTATI
# ============================================================

sites = engine.table(
    out / "02_alignment_hypotheses.tsv"
)

genes = engine.table(
    out / "04_gene_target_candidates.tsv"
)

qc = engine.table(
    out / "01_rna_qc.tsv"
)

rnas = dict(
    engine.fasta(rna_path)
)

map_rows = engine.table(
    tx2gene_path
)

tx2gene = {}

for row in map_rows:

    tx = row["transcript_id"].strip()

    gene = row["gene_id"].strip()

    if (
        not tx
        or not gene
        or (
            tx in tx2gene
            and tx2gene[tx] != gene
        )
    ):
        raise AssertionError(
            "Mapping trascritto-gene non valido"
        )

    tx2gene[tx] = gene

# ============================================================
# VERIFICA CONTROLLO QUALITA
# ============================================================

if len(rnas) != 5 or len(qc) != 5:

    raise AssertionError(
        "Numero RNA/QC diverso da cinque"
    )

if any(
    row["status"] != "eligible"
    for row in qc
):

    raise AssertionError(
        "Almeno uno dei cinque RNA non supera il QC"
    )

if (
    int(manifest.get("n_input_rna", -1)) != 5
    or int(manifest.get("n_eligible_rna", -1)) != 5
):

    raise AssertionError(
        "Conteggi RNA non coerenti"
    )

if (
    int(manifest.get("n_alignment_sites", -1))
    != len(sites)
):

    raise AssertionError(
        "Conteggio dei siti non coerente"
    )

if (
    int(manifest.get("n_genes", -1))
    != len(genes)
):

    raise AssertionError(
        "Conteggio dei geni non coerente"
    )

# ============================================================
# VERIFICA MAPPING DEI SITI
# ============================================================

wanted = defaultdict(list)

for site in sites:

    rid = site["rna_id"]

    tx = site["transcript_id"]

    gene = site["gene_id"]

    if (
        rid not in rnas
        or tx2gene.get(tx) != gene
    ):

        raise AssertionError(
            "Associazione RNA-trascritto-gene incoerente"
        )

    wanted[tx].append(site)

# ============================================================
# VERIFICA POSIZIONI E COMPLEMENTARITA
# ============================================================

verified = 0

for tx, sequence in engine.fasta(
    transcriptome
):

    records = wanted.pop(
        tx,
        None
    )

    if not records:
        continue

    for record in records:

        rid = record["rna_id"]

        expected = engine.revcomp(
            rnas[rid]
        )

        start = (
            int(record["site_start_1based"]) - 1
        )

        end = int(
            record["site_end_1based"]
        )

        observed = sequence[start:end]

        mismatches = sum(
            a != b
            for a, b in zip(
                observed,
                expected
            )
        )

        if (
            start < 0
            or end - start != len(expected)
            or len(observed) != len(expected)
        ):

            raise AssertionError(
                "Coordinate del sito non valide"
            )

        if (
            mismatches != int(record["mismatches"])
            or mismatches > 2
        ):

            raise AssertionError(
                "Mismatch non coerenti"
            )

        identity = (
            len(expected) - mismatches
        ) / len(expected)

        if not math.isclose(
            float(record["identity_fraction"]),
            identity,
            abs_tol=1.1e-6
        ):

            raise AssertionError(
                "Identita di sequenza non coerente"
            )

        verified += 1

if wanted or verified != len(sites):

    raise AssertionError(
        "Impossibile verificare tutti i siti"
    )

# ============================================================
# VERIFICA RANKING GENI CANDIDATI
# ============================================================

observed_genes = {
    site["gene_id"]
    for site in sites
}

ranking_genes = {
    gene["gene_id"]
    for gene in genes
}

if (
    len(genes) != len(observed_genes)
    or ranking_genes != observed_genes
):

    raise AssertionError(
        "Ranking dei geni incoerente con i siti"
    )

for gene in genes:

    gid = gene["gene_id"]

    gene_sites = [
        site
        for site in sites
        if site["gene_id"] == gid
    ]

    if (
        int(gene["n_sites"])
        != len(gene_sites)
    ):

        raise AssertionError(
            "Numero siti per gene errato"
        )

    if (
        int(gene["n_rna"])
        != len({
            site["rna_id"]
            for site in gene_sites
        })
    ):

        raise AssertionError(
            "Numero RNA per gene errato"
        )

# ============================================================
# RAPPORTO DI COLLAUDO
# ============================================================

report = {

    "stato":
        "TEST_FUNZIONALE_SUPERATO",

    "numero_srna_reali":
        len(rnas),

    "sRNA_QC_idonei":
        len(qc),

    "trascritti_scansionati":
        manifest.get(
            "scan",
            {}
        ).get("transcripts_scanned"),

    "trascritti_senza_mapping":
        manifest.get(
            "scan",
            {}
        ).get("transcripts_without_tx2gene"),

    "siti_candidati":
        len(sites),

    "siti_verificati_su_Fasta":
        verified,

    "geni_target_candidati_distinti":
        len(genes),

    "soglia_mismatch":
        2,

    "score_biologico_validato":
        False,

    "espressione_WW_DS":
        "NON_VALUTATA: matrici non fornite",

    "modello_GNN_per_target":
        "NON_APPLICATO",

    "nota":
        "Ipotesi RNA-gene esplorative, "
        "non prova sperimentale"
}

(
    base / "RAPPORTO_STEP16.json"
).write_text(
    json.dumps(
        report,
        ensure_ascii=False,
        indent=2
    ),
    encoding="utf-8"
)

(
    base / "RAPPORTO_STEP16.md"
).write_text(

    "# AI PRIMETTA V3.3 - STEP 16\n\n"

    + "**Test funzionale con cinque sRNA reali. "
      "Non validazione biologica.**\n\n"

    + "\n".join(
        f"- **{key}:** {value}"
        for key, value in report.items()
    )

    + "\n",

    encoding="utf-8"
)

print(
    "[DOCKER] CONTROLLI SUPERATI",
    flush=True
)

print(
    json.dumps(
        report,
        ensure_ascii=False,
        indent=2
    ),
    flush=True
)
'''.replace(
        "__DOCKER_RUN__",
        repr(docker_run)
    )

    # ========================================================
    # 10. AVVIO E VERIFICA RAPPORTO
    # ========================================================

    run_docker(inside)

    report_path = (
        run / "RAPPORTO_STEP16.json"
    )

    if not report_path.is_file():
        raise RuntimeError(
            "RAPPORTO_STEP16.json non generato"
        )

    report = json.loads(
        report_path.read_text(
            encoding="utf-8"
        )
    )

    print()
    print("=" * 66)
    print(
        "[OK] STEP 16 COMPLETATO"
    )
    print("=" * 66)

    print(
        "Cartella risultati:",
        run
    )

    for key, value in report.items():

        print(
            f"{key}: {value}"
        )

    verify_gui()

    print()
    print(
        "[NOTA] Risultati esplorativi, "
        "non validazione biologica indipendente."
    )

    return 0


# ============================================================
# 11. ESECUZIONE
# ============================================================

if __name__ == "__main__":

    try:

        sys.exit(main())

    except (
        OSError,
        ValueError,
        RuntimeError,
        AssertionError,
        TypeError,
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired
    ) as error:

        print(
            "[ERRORE]",
            error,
            file=sys.stderr
        )

        sys.exit(1)
