#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3
STEP 21 - ESPORTAZIONE PORTATILE WINDOWS

Eseguire dalla VM Ubuntu con Visual Studio Code.

Funzioni:
- Identifica i due servizi Docker originali
- Identifica le immagini locali
- Esporta le immagini Docker
- Copia il progetto con modelli, dati e riferimenti
- Converte i bind mount Ubuntu in percorsi relativi
- Genera 2 BAT di avvio e 1 BAT di arresto
- Genera i Compose portatili
- Crea e verifica lo ZIP finale

Non modifica:
- GUI
- Motore V3.3
- Checkpoint
- Dataset originali
- Container esistenti
"""

import argparse
import copy
import datetime as dt
import fnmatch
import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import zipfile


# ============================================================
# 1. CONFIGURAZIONE
# ============================================================

ROOT = Path(__file__).resolve().parent.parent

FOLDER = "AI_PRIMETTA_V33_PORTATILE_WINDOWS"

SERVICE_INFO = [
    (
        "v32",
        ROOT / "docker-compose.yml",
        "maize-ncrna-gnn",
        8501,
        18501,
    ),
    (
        "v33",
        ROOT / "target_prediction_v33_vscode" / "compose.yaml",
        "prediction",
        8502,
        18502,
    ),
]

EXCLUDE_DIRS = {
    ".git",
    ".venv",
    "venv",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    "node_modules",
    ".vscode",
}

EXCLUDE_GLOBS = {
    ".env",
    ".env.*",
    "*.pem",
    "*.key",
    "id_rsa*",
    "id_ed25519*",
    "*.p12",
    "*.pfx",
}

SECRET_RE = re.compile(
    r"TOKEN|PASSWORD|PASSWD|SECRET|PRIVATE_KEY|API_KEY|CREDENTIAL",
    re.I,
)


# ============================================================
# 2. COMANDI DOCKER
# ============================================================

def cmd(arguments):

    process = subprocess.run(
        arguments,
        capture_output=True,
        text=True,
        check=False,
    )

    if process.returncode != 0:
        raise RuntimeError(
            "Comando fallito: "
            + " ".join(map(str, arguments))
            + "\n"
            + process.stderr[-1200:]
        )

    return process.stdout


def compose_config(path):

    if not path.is_file():
        raise RuntimeError(
            f"Compose mancante: {path}"
        )

    response = cmd([
        "docker",
        "compose",
        "-f",
        str(path),
        "config",
        "--format",
        "json",
    ])

    return json.loads(response)


def image_present(name):

    process = subprocess.run(
        ["docker", "image", "inspect", name],
        capture_output=True,
        check=False,
    )

    return process.returncode == 0


def resolve_image(config, service, fallback):

    service_cfg = config["services"][service]

    candidates = [
        service_cfg.get("image"),
        f"{config.get('name', ROOT.name)}-{service}:latest",
        fallback,
    ]

    for image in dict.fromkeys(
        name for name in candidates if name
    ):
        if image_present(image):
            return image

    raise RuntimeError(
        f"Nessuna immagine Docker locale per {service}. "
        f"Cercate: {candidates}"
    )


# ============================================================
# 3. CONVERSIONE COMPOSE PER WINDOWS
# ============================================================

def safe_relative(path):

    resolved = Path(path).resolve()

    try:
        return resolved.relative_to(ROOT)
    except ValueError:
        raise RuntimeError(
            "Bind mount esterno al progetto, "
            f"non esportato: {resolved}"
        )


def rewrite_service(
    config,
    name,
    service_name,
    image,
    internal_port,
    windows_port,
):

    service = copy.deepcopy(
        config["services"][service_name]
    )

    for key in (
        "build",
        "develop",
        "container_name",
        "pull_policy",
        "env_file",
        "profiles",
    ):
        service.pop(key, None)

    service["image"] = image

    service["ports"] = [{
        "target": internal_port,
        "published": str(windows_port),
        "host_ip": "127.0.0.1",
        "protocol": "tcp",
    }]

    service["restart"] = "unless-stopped"

    if service.get("depends_on"):
        raise RuntimeError(
            f"Servizi dipendenti non inclusi: "
            f"{service['depends_on']}"
        )

    if any(
        service.get(key)
        for key in ("secrets", "configs", "devices")
    ):
        raise RuntimeError(
            f"Il servizio {service_name} richiede "
            "adattamento manuale di secrets/configs/devices"
        )

    # Controllo credenziali nella configurazione.
    environment = service.get("environment") or {}

    if isinstance(environment, dict):

        for key, value in environment.items():

            if (
                SECRET_RE.search(key)
                and value not in (None, "")
            ):
                raise RuntimeError(
                    "Variabile potenzialmente riservata "
                    f"nel Compose: {key}"
                )

    # Conversione dei volumi Ubuntu.
    for volume in service.get("volumes", []):

        if (
            not isinstance(volume, dict)
            or volume.get("type") != "bind"
        ):
            raise RuntimeError(
                "Volume non portatile individuato: "
                f"{volume}"
            )

        source = volume.get("source")

        if not source:
            raise RuntimeError(
                "Bind mount senza sorgente"
            )

        relative = safe_relative(source)

        if not Path(source).exists():
            raise RuntimeError(
                f"Sorgente inesistente: {source}"
            )

        if not relative.parts:
            portable_source = "./progetto"
        else:
            portable_source = (
                "./progetto/" + relative.as_posix()
            )

        volume["source"] = portable_source
        volume.pop("bind", None)

    # Evita riferimenti assoluti residui alla VM.
    if str(ROOT) in json.dumps(
        service, ensure_ascii=False
    ):
        raise RuntimeError(
            "Il servizio contiene altri percorsi "
            "assoluti Ubuntu da adattare."
        )

    service.pop("networks", None)

    return {
        "name": f"ai_primetta_portabile_{name}",
        "services": {
            service_name: service
        },
    }


# ============================================================
# 4. INVENTARIO FILE DEL PROGETTO
# ============================================================

def safe_file(path):

    return not any(
        fnmatch.fnmatch(path.name, pattern)
        for pattern in EXCLUDE_GLOBS
    )


def files_to_package():

    skipped = []

    for current, directories, filenames in os.walk(
        ROOT,
        followlinks=False,
    ):

        current = Path(current)

        valid_directories = []

        for directory in sorted(directories):

            path = current / directory

            if (
                directory in EXCLUDE_DIRS
                or path.is_symlink()
                or not safe_file(path)
            ):
                skipped.append(
                    str(path.relative_to(ROOT))
                )
            else:
                valid_directories.append(directory)

        directories[:] = valid_directories

        for filename in sorted(filenames):

            path = current / filename

            if (
                path.is_symlink()
                or not path.is_file()
                or not safe_file(path)
            ):
                skipped.append(
                    str(path.relative_to(ROOT))
                )
                continue

            yield path, skipped


# ============================================================
# 5. GENERAZIONE FILE WINDOWS BAT
# ============================================================

def make_bat(mode):

    names = {
        "v32": "AVVIO GUI ORIGINALE V3.2",
        "v33": "AVVIO GUI TARGET PREDICTION V3.3",
        "stop": "ARRESTO ENTRAMBE LE GUI",
    }

    return (
        "@echo off\r\n"
        "setlocal\r\n"
        'cd /d "%~dp0"\r\n'
        f"echo AI PRIMETTA - {names[mode]}\r\n"
        "powershell.exe -NoProfile "
        "-ExecutionPolicy Bypass "
        '-File "%~dp0gestione_portatile.ps1" '
        f"-Azione {mode}\r\n"
        "if errorlevel 1 ( "
        "echo [ERRORE] Controlla i messaggi sopra. "
        "& pause & exit /b 1 )\r\n"
        "endlocal\r\n"
    )


# ============================================================
# 6. SCRIPT POWERSHELL WINDOWS
# ============================================================

POWERSHELL = r'''
param(
    [ValidateSet('v32','v33','stop')]
    [string]$Azione
)

$ErrorActionPreference = 'Stop'

$root = Split-Path -Parent $MyInvocation.MyCommand.Path

Set-Location $root

$cfg = Get-Content `
    -LiteralPath (Join-Path $root 'config_portatile.json') `
    -Raw | ConvertFrom-Json


function Run-Docker([string[]]$ArgsList) {

    & docker @ArgsList

    if ($LASTEXITCODE -ne 0) {
        throw "Docker exit code $LASTEXITCODE"
    }
}


try {

    if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
        throw 'Docker Desktop non trovato.'
    }

    & docker info --format '{{.ServerVersion}}' `
        1>$null 2>$null

    if ($LASTEXITCODE -ne 0) {
        throw 'Docker Engine non disponibile.'
    }

    if ($Azione -ne 'stop') {

        # Importazione immagini, solo se mancanti.
        $missing = $false

        foreach ($img in @($cfg.images)) {

            & docker image inspect "$img" `
                1>$null 2>$null

            if ($LASTEXITCODE -ne 0) {
                $missing = $true
                break
            }
        }

        if ($missing) {

            $archive = Join-Path `
                $root 'images\docker_images.tar.gz'

            if (-not (Test-Path -LiteralPath $archive)) {
                throw "Archivio immagini mancante: $archive"
            }

            Write-Host '[INFO] Importazione immagini Docker...'

            Run-Docker -ArgsList @(
                'load',
                '-i',
                $archive
            )
        }

        $item = $cfg.services.$Azione

        $compose = Join-Path $root $item.compose

        Run-Docker -ArgsList @(
            'compose',
            '-f',
            $compose,
            'up',
            '-d',
            '--no-build',
            '--pull',
            'never',
            $item.service
        )

        $url = "http://127.0.0.1:$($item.windows_port)"

        $health = "$url/_stcore/health"

        $ok = $false

        # Verifica che Streamlit risponda.
        for ($i = 0; $i -lt 40; $i++) {

            try {

                $response = Invoke-WebRequest `
                    -Uri $health `
                    -UseBasicParsing `
                    -TimeoutSec 4

                if (
                    $response.StatusCode -eq 200 `
                    -and $response.Content.Trim() -eq 'ok'
                ) {
                    $ok = $true
                    break
                }

            } catch { }

            Start-Sleep -Seconds 2
        }

        if (-not $ok) {

            & docker compose -f $compose ps

            throw "GUI non raggiungibile: $health"
        }

        Write-Host "[OK] GUI disponibile: $url"

        Start-Process $url

    } else {

        foreach ($which in @('v32', 'v33')) {

            $item = $cfg.services.$which

            $compose = Join-Path $root $item.compose

            Run-Docker -ArgsList @(
                'compose',
                '-f',
                $compose,
                'stop',
                $item.service
            )
        }

        Write-Host '[OK] Entrambe le GUI fermate.'
    }

} catch {

    Write-Host "[ERRORE] $($_.Exception.Message)" `
        -ForegroundColor Red

    exit 1
}
'''


# ============================================================
# 7. ISTRUZIONI GENERATE NEL PACCHETTO
# ============================================================

README = r"""
AI PRIMETTA V3.3 - PACCHETTO PORTATILE WINDOWS
================================================

Il pacchetto non dipende dalla VM Ubuntu originale.

REQUISITI:
- Windows 10/11 x64.
- Docker Desktop installato e operativo.
- Modalita Linux containers / WSL2.
- Spazio su disco sufficiente.

UTILIZZO:

1. Estrai lo ZIP, ad esempio in C:\AI_PRIMETTA
2. Avvia Docker Desktop.
3. Esegui 01_AVVIA_GUI_ORIGINALE_V32.bat
4. Esegui 02_AVVIA_GUI_PREDIZIONE_V33.bat
5. Esegui 03_FERMA_ENTRAMBE_LE_GUI.bat

PAGINE:

GUI V3.2:
http://127.0.0.1:18501

GUI predizione V3.3:
http://127.0.0.1:18502

Lo ZIP contiene immagini Docker esportate dalla VM,
codice, modelli, dati e configurazioni relative.

Non vengono modificati i checkpoint.

ATTENZIONE:
La GUI di predizione V3.3 e i modelli ML/GNN
costituiscono componenti distinti.

La disponibilita delle due pagine non prova
l'integrazione diretta GUI-ML/GNN.

I file .env e le credenziali private sono esclusi.

Controllare il contenuto prima di distribuire
l'archivio a terzi.

Il pacchetto puo occupare diversi gigabyte.

Utilizzare NTFS o exFAT se il file supera 4 GB.
"""


# ============================================================
# 8. ESPORTAZIONE IMMAGINI DOCKER
# ============================================================

def save_images(images, archive):

    print(
        "\n[3] ESPORTAZIONE IMMAGINI DOCKER",
        flush=True,
    )

    for image in images:
        print("Immagine:", image)

    error_file = archive.parent / "docker_save_stderr.txt"

    with error_file.open("wb") as stderr:

        process = subprocess.Popen(
            ["docker", "save", *images],
            stdout=subprocess.PIPE,
            stderr=stderr,
        )

        try:

            with gzip.open(
                archive,
                "wb",
                compresslevel=1,
            ) as compressed:

                while True:

                    block = process.stdout.read(
                        4 * 1024 * 1024
                    )

                    if not block:
                        break

                    compressed.write(block)

        finally:

            if process.stdout:
                process.stdout.close()

        code = process.wait()

    if code != 0:

        error = error_file.read_text(
            errors="replace"
        )

        raise RuntimeError(
            "Esportazione immagini fallita: "
            + error[-1200:]
        )

    print(
        "[OK] Archivio immagini:",
        round(
            archive.stat().st_size / (1024 ** 3),
            2
        ),
        "GiB",
    )


# ============================================================
# 9. PROCEDURA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--gui-only",
        action="store_true",
        help="Esclude il runtime ML/GNN aggiuntivo",
    )

    args = parser.parse_args()

    print("=" * 68)
    print("AI PRIMETTA V3.3")
    print("STEP 21 - CREAZIONE ZIP PORTATILE WINDOWS")
    print("=" * 68)

    if not (
        (ROOT / "models").is_dir()
        and
        (ROOT / "target_prediction_v33_vscode").is_dir()
    ):
        raise RuntimeError(
            "Struttura del progetto AI PRIMETTA non trovata"
        )

    # Docker disponibile?
    cmd([
        "docker",
        "info",
        "--format",
        "{{.ServerVersion}}",
    ])

    manifest = {}

    services = {}

    images = []

    print("\n[1] CONTROLLO DEI DUE SERVIZI")

    for (
        key,
        compose_path,
        service_name,
        internal_port,
        windows_port,
    ) in SERVICE_INFO:

        config = compose_config(compose_path)

        if service_name not in config.get("services", {}):
            raise RuntimeError(
                f"Servizio non trovato: {service_name}"
            )

        if key == "v32":
            fallback = (
                "ncrna_rnacentral_ai_v3_prjna419326-"
                "maize-ncrna-gnn:latest"
            )
        else:
            fallback = (
                "ai_primetta_target_v33-prediction:latest"
            )

        image = resolve_image(
            config,
            service_name,
            fallback,
        )

        manifest[key] = rewrite_service(
            config,
            key,
            service_name,
            image,
            internal_port,
            windows_port,
        )

        services[key] = {
            "compose": f"compose_{key}.yaml",
            "service": service_name,
            "windows_port": windows_port,
            "image": image,
        }

        images.append(image)

        print(
            "[OK]",
            key,
            "| Immagine:",
            image,
            "| Porta Windows:",
            windows_port,
        )

    # Immagine scientifica ML/GNN opzionale.
    if not args.gui_only:

        ml_image = "ai_primetta_v33_ml_compat:19.5"

        if image_present(ml_image):

            images.append(ml_image)

            print(
                "[OK] Runtime ML/GNN incluso:",
                ml_image,
            )

        else:

            print(
                "[ATTENZIONE] Runtime ML/GNN "
                "aggiuntivo non disponibile"
            )

    images = list(dict.fromkeys(images))

    configuration = {
        "images": images,
        "services": services,
    }

    # ========================================================
    # 10. DESTINAZIONE ARCHIVIO
    # ========================================================

    stamp = dt.datetime.now(
        dt.timezone.utc
    ).strftime("%Y%m%dT%H%M%SZ")

    destination = ROOT.parent / (
        f"{FOLDER}_{stamp}.zip"
    )

    print("\n[2] DESTINAZIONE ZIP")
    print(destination)

    # Cartella temporanea esterna al progetto.
    with tempfile.TemporaryDirectory(
        prefix="ai_primetta_portatile_",
        dir=ROOT.parent,
    ) as temporary:

        temp = Path(temporary)

        images_archive = (
            temp / "docker_images.tar.gz"
        )

        save_images(
            images,
            images_archive,
        )

        print("\n[4] CREAZIONE ZIP COMPLETO", flush=True)

        skipped = []

        file_count = 0

        with zipfile.ZipFile(
            destination,
            "w",
            allowZip64=True,
        ) as archive:

            base = FOLDER + "/"

            # Compose Windows.
            for key in ("v32", "v33"):

                archive.writestr(
                    base + f"compose_{key}.yaml",
                    json.dumps(
                        manifest[key],
                        indent=2,
                        ensure_ascii=False,
                    ) + "\n",
                )

            # Configurazione di gestione.
            archive.writestr(
                base + "config_portatile.json",
                json.dumps(
                    configuration,
                    indent=2,
                ) + "\n",
            )

            # Gestore PowerShell.
            archive.writestr(
                base + "gestione_portatile.ps1",
                POWERSHELL,
            )

            # I due BAT di avvio.
            archive.writestr(
                base + "01_AVVIA_GUI_ORIGINALE_V32.bat",
                make_bat("v32"),
            )

            archive.writestr(
                base + "02_AVVIA_GUI_PREDIZIONE_V33.bat",
                make_bat("v33"),
            )

            # BAT di stop.
            archive.writestr(
                base + "03_FERMA_ENTRAMBE_LE_GUI.bat",
                make_bat("stop"),
            )

            archive.writestr(
                base + "ISTRUZIONI.txt",
                README,
            )

            # Docker images compress.
            archive.write(
                images_archive,
                arcname=(
                    base
                    + "images/docker_images.tar.gz"
                ),
                compress_type=zipfile.ZIP_STORED,
            )

            # Copia del progetto.
            for path, skipped_ref in files_to_package():

                skipped = skipped_ref

                relative = (
                    path.relative_to(ROOT)
                    .as_posix()
                )

                archive.write(
                    path,
                    arcname=(
                        base
                        + "progetto/"
                        + relative
                    ),
                    compress_type=zipfile.ZIP_DEFLATED,
                    compresslevel=1,
                )

                file_count += 1

                if file_count % 500 == 0:
                    print(
                        "File copiati:",
                        file_count,
                        flush=True,
                    )

            archive.writestr(
                base + "FILE_ESCLUSI.txt",
                "\n".join(skipped) + "\n",
            )

            export_manifest = {
                "created_utc": stamp,
                "project": ROOT.name,
                "docker_images": images,
                "files_copied": file_count,
                "files_excluded": skipped,
            }

            archive.writestr(
                base + "MANIFEST_EXPORT.json",
                json.dumps(
                    export_manifest,
                    indent=2,
                ),
            )

    # ========================================================
    # 11. VERIFICA INTEGRITA ZIP
    # ========================================================

    print("\n[5] VERIFICA ZIP")

    with zipfile.ZipFile(destination) as archive:

        bad_file = archive.testzip()

        if bad_file:
            raise RuntimeError(
                "ZIP danneggiato: "
                + bad_file
            )

    print()
    print("=" * 68)
    print("STEP 21 COMPLETATO")
    print("=" * 68)

    print("[OK] ZIP portatile creato")
    print("[OK] File progetto:", file_count)
    print("[OK] Immagini Docker:", len(images))
    print("[OK] Integrità archivio verificata")
    print("[OK] GUI originali non modificate")
    print("[OK] Nessun training")

    print("\nPERCORSO ZIP:")
    print(destination)

    print(
        "\n[NOTA] Il pacchetto deve essere "
        "collaudato su Windows con Docker Desktop."
    )


# ============================================================
# 12. AVVIO
# ============================================================

if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print(
            "\n[ERRORE]",
            exc,
            file=sys.stderr,
        )

        sys.exit(1)