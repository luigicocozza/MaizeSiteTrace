#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 17.9.1 - VERIFICA CHECKSUM UFFICIALE B73 NAM-5.0

- Download del manifesto MD5 ufficiale MaizeGDB
- Estrazione checksum del FASTA genomico
- Calcolo MD5 e SHA256 del FASTA locale
- Confronto dei valori
- Report JSON

GUI, motore, modelli e dati originali invariati.
"""

import hashlib
import json
import re
import urllib.request

from datetime import datetime, timezone
from pathlib import Path


ROOT = (
    Path.home()
    / "progetti"
    / "ncrna_rnacentral_ai_v3_prjna419326"
)

BASE = ROOT / "target_prediction_v33_vscode"

ASSEMBLY = "Zm-B73-REFERENCE-NAM-5.0"

FASTA_NAME = ASSEMBLY + ".fa.gz"

FASTA = (
    BASE / "references" / FASTA_NAME
)

MANIFEST_URL = (
    "https://download.maizegdb.org/"
    "Zm-B73-REFERENCE-NAM-5.0/"
    "Zm-B73-REFERENCE-NAM-5.0.md5"
)

RUNS = BASE / "runs"


def checksum(path, algorithm):

    digest = hashlib.new(algorithm)

    with path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(4 * 1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def read_manifest():

    request = urllib.request.Request(
        MANIFEST_URL,
        headers={
            "User-Agent": "AI-PRIMETTA-Research/3.3"
        }
    )

    with urllib.request.urlopen(
        request,
        timeout=60
    ) as response:

        raw = response.read(1_000_000 + 1)

    if len(raw) > 1_000_000:
        raise RuntimeError(
            "Manifesto eccessivamente grande"
        )

    return raw.decode(
        "utf-8",
        errors="replace"
    )


def extract_md5(text):

    matches = []

    for line in text.splitlines():

        line = line.strip()

        if not line:
            continue

        # Formato md5sum:
        # 32 caratteri esadecimali + nome file.
        match = re.match(
            r"^([0-9a-fA-F]{32})\s+\*?(.+)$",
            line
        )

        if match:

            value = match.group(1).lower()

            filename = match.group(2).strip()

            if Path(filename).name == FASTA_NAME:
                matches.append(value)

            continue

        # Formato alternativo:
        # MD5 (filename) = hash
        match = re.match(
            r"^MD5\s*\((.+)\)\s*=\s*"
            r"([0-9a-fA-F]{32})$",
            line,
            re.IGNORECASE
        )

        if match:

            filename = match.group(1).strip()

            value = match.group(2).lower()

            if Path(filename).name == FASTA_NAME:
                matches.append(value)

    matches = list(set(matches))

    if len(matches) != 1:

        raise RuntimeError(
            "Checksum ufficiale non identificato "
            "univocamente nel manifesto. "
            "Nessuna validazione dichiarata."
        )

    return matches[0]


def main():

    print("=" * 70)
    print("AI PRIMETTA V3.3 - STEP 17.9.1")
    print("VERIFICA MD5 UFFICIALE GENOMA B73")
    print("=" * 70)

    if not FASTA.is_file():
        raise SystemExit(
            f"[ERRORE] FASTA mancante: {FASTA}"
        )

    print("\n[1] FILE LOCALE")

    print("FASTA:", FASTA)

    print(
        "Dimensione MB:",
        round(FASTA.stat().st_size / 1048576, 3)
    )

    print("\n[2] MANIFESTO UFFICIALE")

    print("Fonte:", MANIFEST_URL)

    try:

        manifest = read_manifest()

    except Exception as error:

        raise SystemExit(
            "[ERRORE] Download manifesto fallito: "
            + str(error)
            + "\nIl FASTA locale rimane invariato."
        )

    print("[OK] Manifesto acquisito")

    try:

        official_md5 = extract_md5(manifest)

    except Exception as error:

        print("\n[INFO] Contenuto manifesto:")
        print(manifest[:12000])

        raise SystemExit(
            "[ERRORE] " + str(error)
        )

    print("MD5 ufficiale:", official_md5)

    print("\n[3] CALCOLO CHECKSUM LOCALI")

    local_md5 = checksum(FASTA, "md5")

    local_sha256 = checksum(FASTA, "sha256")

    print("MD5 locale:", local_md5)

    print("SHA256 locale:", local_sha256)

    print("\n[4] CONFRONTO")

    verified = (
        local_md5.lower() == official_md5.lower()
    )

    if verified:

        status = "CHECKSUM_UFFICIALE_CONFORME"

        print(
            "[OK] MD5 locale identico "
            "al checksum ufficiale MaizeGDB"
        )

    else:

        status = "CHECKSUM_NON_CONFORME"

        print(
            "[ERRORE] MD5 differente "
            "dal checksum ufficiale"
        )

    print("\n[5] REPORT")

    timestamp = datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")

    output = (
        RUNS
        / f"STEP17_9_1_CHECKSUM_{timestamp}"
    )

    output.mkdir(
        parents=True,
        exist_ok=False
    )

    manifest_path = (
        output / "MANIFESTO_UFFICIALE.md5"
    )

    manifest_path.write_text(
        manifest,
        encoding="utf-8"
    )

    report = {
        "progetto": "AI PRIMETTA V3.3",
        "step": "17.9.1",
        "stato": status,
        "assembly": ASSEMBLY,
        "file": str(FASTA),
        "manifest_url": MANIFEST_URL,
        "manifest_file": str(manifest_path),
        "md5_ufficiale": official_md5,
        "md5_locale": local_md5,
        "sha256_locale": local_sha256,
        "checksum_conforme": verified,
        "genoma_modificato": False,
        "gui_modificata": False,
        "modelli_modificati": False,
    }

    report_path = (
        output / "RAPPORTO_STEP17_9_1.json"
    )

    report_path.write_text(
        json.dumps(
            report,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )

    print("[OK] Report:", report_path)

    print("\n" + "=" * 70)
    print("STEP 17.9.1 TERMINATO")
    print("=" * 70)

    print("Stato:", status)

    print("[OK] GUI V3.2 invariata")
    print("[OK] Motore V3.3 invariato")
    print("[OK] Modelli invariati")

    if not verified:
        raise SystemExit(
            "[STOP] Checksum non conforme. "
            "Non utilizzare questo FASTA "
            "per il mapping genomico."
        )


if __name__ == "__main__":
    main()