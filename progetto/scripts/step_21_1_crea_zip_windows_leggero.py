#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 21.1
CREAZIONE ZIP WINDOWS LEGGERO

Esclude FASTQ e altri file grezzi di sequenziamento.
Preserva GUI, modelli, dati elaborati e riferimenti.

Esecuzione:
    python3 scripts/step_21_1_crea_zip_windows_leggero.py

Creazione ZIP:
    python3 scripts/step_21_1_crea_zip_windows_leggero.py --build

Non modifica il repository originale.
"""

import argparse
import copy
import datetime as dt
import fnmatch
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT.parent
PACK_NAME = "AI_PRIMETTA_V33_WINDOWS_LIGHT"

SERVICES = {
    "v32": {
        "compose": ROOT / "docker-compose.yml",
        "service": "maize-ncrna-gnn",
        "fallback": (
            "ncrna_rnacentral_ai_v3_prjna419326-"
            "maize-ncrna-gnn:latest"
        ),
        "internal": 8501,
        "windows": 18501,
    },
    "v33": {
        "compose": (
            ROOT / "target_prediction_v33_vscode/compose.yaml"
        ),
        "service": "prediction",
        "fallback": "ai_primetta_target_v33-prediction:latest",
        "internal": 8502,
        "windows": 18502,
    },
}

SKIP_DIRS = {
    ".git",
    ".venv",
    "venv",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    "node_modules",
    "delivery",
    "backups",
    "logs",
}

SKIP_RAW_EXT = (
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".bam",
    ".cram",
    ".sam",
    ".sra",
    ".bai",
    ".csi",
    ".bt2",
    ".bt2l",
)

SENSITIVE_NAMES = (
    ".env",
    ".env.*",
    "id_rsa",
    "id_ed25519",
    "*.pem",
    "*.key",
    "*.p12",
    "*.pfx",
)

ZIP_COMPRESSED = (
    ".gz",
    ".zip",
    ".7z",
    ".jpg",
    ".jpeg",
    ".png",
    ".webp",
    ".pt",
    ".pth",
    ".joblib",
    ".parquet",
    ".npz",
    ".tif",
    ".tiff",
)


# ==========================================================
# FILE BAT WINDOWS
# ==========================================================

BAT = r'''@echo off
setlocal
cd /d "%~dp0"

powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0gestione_portatile.ps1" -Azione __ACTION__

if errorlevel 1 (
    echo [ERRORE] Controllare Docker Desktop e i messaggi.
    pause
    exit /b 1
)

endlocal
'''


# ==========================================================
# GESTIONE DOCKER DA WINDOWS
# ==========================================================

POWERSHELL = r'''
param(
    [ValidateSet('v32','v33','stop')]
    [string]$Azione
)

$ErrorActionPreference = 'Stop'

Set-Location (
    Split-Path -Parent $MyInvocation.MyCommand.Path
)

$c = Get-Content -Raw 'config_portatile.json' |
    ConvertFrom-Json

function Dc([string[]]$a) {
    & docker @a
    if ($LASTEXITCODE -ne 0) {
        throw "Comando Docker fallito"
    }
}

try {

    if (-not (
        Get-Command docker -ErrorAction SilentlyContinue
    )) {
        throw 'Docker Desktop non installato.'
    }

    & docker info --format '{{.ServerVersion}}' `
        1>$null 2>$null

    if ($LASTEXITCODE -ne 0) {
        throw 'Avvia Docker Desktop in modalita Linux containers.'
    }

    if ($Azione -eq 'stop') {

        foreach ($k in @('v32','v33')) {

            $s = $c.services.$k

            Dc -a @(
                'compose',
                '-f',
                $s.compose,
                'stop',
                $s.service
            )
        }

        Write-Host '[OK] Entrambe le GUI fermate.'

    } else {

        $mustLoad = $false

        foreach ($img in $c.images) {

            & docker image inspect $img `
                1>$null 2>$null

            if ($LASTEXITCODE -ne 0) {
                $mustLoad = $true
                break
            }
        }

        if ($mustLoad) {

            $archive = 'images\docker_images.tar.gz'

            if (-not (Test-Path $archive)) {
                throw 'Archivio immagini Docker mancante.'
            }

            Write-Host '[INFO] Importazione immagini...'

            Dc -a @(
                'load',
                '-i',
                $archive
            )
        }

        $s = $c.services.$Azione

        Dc -a @(
            'compose',
            '-f',
            $s.compose,
            'up',
            '-d',
            '--no-build',
            '--pull',
            'never',
            $s.service
        )

        $url = "http://127.0.0.1:$($s.windows_port)"

        $ok = $false

        for ($n = 0; $n -lt 45; $n++) {

            try {

                $response = Invoke-WebRequest `
                    -Uri "$url/_stcore/health" `
                    -UseBasicParsing `
                    -TimeoutSec 4

                if (
                    $response.StatusCode -eq 200 `
                    -and $response.Content.Trim() -eq 'ok'
                ) {
                    $ok = $true
                    break
                }

            } catch {}

            Start-Sleep -Seconds 2
        }

        if (-not $ok) {
            throw "GUI non raggiungibile su $url"
        }

        Write-Host "[OK] GUI disponibile: $url"

        Start-Process $url
    }

} catch {

    Write-Host "[ERRORE] $($_.Exception.Message)" `
        -ForegroundColor Red

    exit 1
}
'''


README = """
AI PRIMETTA V3.3
PACCHETTO PORTATILE WINDOWS LEGGERO

REQUISITI:
- Windows 10/11 x64
- Docker Desktop installato
- Linux containers attivi

ISTRUZIONI:

1. Estrarre completamente lo ZIP su disco NTFS.
2. Avviare Docker Desktop.
3. Eseguire 01_AVVIA_GUI_ORIGINALE_V32.bat
4. Eseguire 02_AVVIA_GUI_PREDIZIONE_V33.bat
5. Usare 03_FERMA_ENTRAMBE_LE_GUI.bat per lo stop.

INDIRIZZI:

GUI V3.2:
http://127.0.0.1:18501

GUI V3.3:
http://127.0.0.1:18502

I FASTQ grezzi non sono inclusi.
Non sono state modificate GUI o checkpoint.

Il pacchetto Windows deve essere collaudato sul PC di
destinazione prima di considerarlo operativo.
"""


# ==========================================================
# UTILITA DOCKER
# ==========================================================

def command(args):

    p = subprocess.run(
        args,
        capture_output=True,
        text=True,
        check=False,
    )

    if p.returncode:
        raise RuntimeError(
            " ".join(map(str, args))
            + "\n"
            + p.stderr[-1600:]
        )

    return p.stdout


def select_image(config, key, cfg):

    image = (
        config["services"][cfg["service"]]
        .get("image")
    )

    if not image:
        image = (
            f"{config.get('name', ROOT.name)}-"
            f"{cfg['service']}:latest"
        )

    for candidate in (image, cfg["fallback"]):

        p = subprocess.run(
            ["docker", "image", "inspect", candidate],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        if p.returncode == 0:
            return candidate

    raise RuntimeError(
        f"Immagine Docker mancante: {key}"
    )


# ==========================================================
# COMPOSE PORTATILE
# ==========================================================

def portable_service(cfg, service, image):

    s = copy.deepcopy(service)

    unsupported = (
        "depends_on",
        "volumes_from",
        "secrets",
        "configs",
        "devices",
        "links",
        "extends",
    )

    for key in unsupported:
        if s.get(key):
            raise RuntimeError(
                f"Servizio con dipendenza non portatile: {key}"
            )

    if s.get("env_file"):
        raise RuntimeError(
            "env_file presente: richiede gestione separata"
        )

    env = s.get("environment") or {}

    if isinstance(env, dict):

        for key, value in env.items():

            if any(
                word in key.upper()
                for word in (
                    "PASSWORD",
                    "TOKEN",
                    "SECRET",
                    "API_KEY",
                    "PRIVATE_KEY",
                )
            ) and value not in ("", None):

                raise RuntimeError(
                    f"Variabile riservata nel Compose: {key}"
                )

    for key in (
        "build",
        "develop",
        "container_name",
        "pull_policy",
        "networks",
    ):
        s.pop(key, None)

    s["image"] = image

    s["ports"] = [{
        "target": cfg["internal"],
        "published": str(cfg["windows"]),
        "host_ip": "127.0.0.1",
        "protocol": "tcp",
    }]

    s["restart"] = "unless-stopped"

    volumes = []

    for volume in s.get("volumes", []):

        if (
            not isinstance(volume, dict)
            or volume.get("type") != "bind"
        ):
            raise RuntimeError(
                "Volume Docker non portatile: "
                + str(volume)
            )

        source = Path(volume["source"]).resolve()

        try:
            relative = source.relative_to(ROOT.resolve())
        except ValueError:
            raise RuntimeError(
                f"Mount esterno al progetto: {source}"
            )

        if not source.exists():
            raise RuntimeError(
                f"Sorgente Docker mancante: {source}"
            )

        new_volume = copy.deepcopy(volume)

        new_volume["source"] = (
            "./progetto"
            + (
                "/" + relative.as_posix()
                if relative.parts
                else ""
            )
        )

        new_volume.pop("bind", None)

        volumes.append(new_volume)

    if volumes:
        s["volumes"] = volumes

    if str(ROOT) in json.dumps(s):
        raise RuntimeError(
            "Riferimenti assoluti Ubuntu nel Compose"
        )

    return {
        "name": (
            "ai_primetta_light_"
            + cfg["service"].replace("-", "_")
        ),
        "services": {
            cfg["service"]: s
        },
    }


# ==========================================================
# FILTRO FILE GREZZI E PRIVATI
# ==========================================================

def reason_to_skip(path):

    relative = path.relative_to(ROOT)

    if any(
        part in SKIP_DIRS
        for part in relative.parts
    ):
        return "cartella di sviluppo, backup o delivery"

    lower = path.name.lower()

    if lower.endswith(SKIP_RAW_EXT):
        return "file grezzo di sequenziamento/allineamento"

    if any(
        fnmatch.fnmatch(lower, pattern)
        for pattern in SENSITIVE_NAMES
    ):
        return "file sensibile o credenziale"

    if path.is_symlink():
        return "collegamento simbolico"

    return None


def inventory():

    included = []
    excluded = []

    for current, dirs, names in os.walk(
        ROOT,
        followlinks=False,
    ):

        current = Path(current)

        valid_dirs = []

        for dirname in sorted(dirs):

            path = current / dirname

            reason = reason_to_skip(path)

            if reason:
                excluded.append((
                    path.relative_to(ROOT).as_posix() + "/",
                    reason,
                ))
            else:
                valid_dirs.append(dirname)

        dirs[:] = valid_dirs

        for filename in sorted(names):

            path = current / filename

            if not (
                path.is_file()
                or path.is_symlink()
            ):
                continue

            reason = reason_to_skip(path)

            if reason:
                excluded.append((
                    path.relative_to(ROOT).as_posix(),
                    reason,
                ))
                continue

            included.append((
                path,
                path.stat().st_size,
            ))

    return included, excluded


# ==========================================================
# CHECKSUM SHA256
# ==========================================================

def sha256(path):

    digest = hashlib.sha256()

    with path.open("rb") as handle:

        for block in iter(
            lambda: handle.read(4 * 1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


# ==========================================================
# ESPORTAZIONE IMMAGINI DOCKER
# ==========================================================

def export_images(images, destination):

    process = subprocess.Popen(
        ["docker", "save", *images],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    try:

        with gzip.open(
            destination,
            "wb",
            compresslevel=1,
        ) as output:

            while True:

                chunk = process.stdout.read(
                    4 * 1024 * 1024
                )

                if not chunk:
                    break

                output.write(chunk)

        process.stdout.close()

        errors = process.stderr.read().decode(
            errors="replace"
        )

        code = process.wait()

        if code:
            raise RuntimeError(
                "docker save fallito: "
                + errors[-1000:]
            )

    finally:

        if process.poll() is None:
            process.kill()
            process.wait()


# ==========================================================
# PROCEDURA PRINCIPALE
# ==========================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--build",
        action="store_true",
        help="Crea realmente lo ZIP",
    )

    parser.add_argument(
        "--max-gb",
        type=float,
        default=12,
        help="Dimensione massima dei file: 12 GiB",
    )

    parser.add_argument(
        "--include-ml-runtime",
        action="store_true",
        help="Include immagine ML compatibile",
    )

    args = parser.parse_args()

    if not (
        (ROOT / "app/streamlit_app.py").is_file()
        and
        (
            ROOT
            / "target_prediction_v33_vscode/app_v33.py"
        ).is_file()
    ):
        raise RuntimeError(
            "Repository AI PRIMETTA non individuato"
        )

    print("=" * 70)
    print("AI PRIMETTA V3.3")
    print("STEP 21.1 - PACCHETTO WINDOWS LEGGERO")
    print("=" * 70, flush=True)

    command([
        "docker",
        "info",
        "--format",
        "{{.ServerVersion}}",
    ])

    configs = {}
    images = []
    services = {}

    print("\n[1] VERIFICA SERVIZI")

    for key, cfg in SERVICES.items():

        config = json.loads(
            command([
                "docker",
                "compose",
                "-f",
                str(cfg["compose"]),
                "config",
                "--format",
                "json",
            ])
        )

        if cfg["service"] not in config.get(
            "services", {}
        ):
            raise RuntimeError(
                f"Servizio non trovato: {cfg['service']}"
            )

        image = select_image(
            config, key, cfg
        )

        configs[key] = portable_service(
            {**cfg},
            config["services"][cfg["service"]],
            image,
        )

        services[key] = {
            "compose": f"compose_{key}.yaml",
            "service": cfg["service"],
            "windows_port": cfg["windows"],
        }

        images.append(image)

        print(
            "[OK]",
            key,
            "|",
            image,
        )

    if args.include_ml_runtime:

        image = "ai_primetta_v33_ml_compat:19.5"

        check = subprocess.run(
            ["docker", "image", "inspect", image],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        if check.returncode != 0:
            raise RuntimeError(
                f"Immagine ML mancante: {image}"
            )

        images.append(image)

    images = list(dict.fromkeys(images))

    # ======================================================
    # INVENTARIO PRELIMINARE
    # ======================================================

    print("\n[2] INVENTARIO DATI")

    files, skipped = inventory()

    total = sum(
        size for _, size in files
    )

    total_gb = total / (1024 ** 3)

    print("File da includere:", len(files))
    print("Percorsi esclusi:", len(skipped))

    print(
        "Dimensione dati non compressi:",
        f"{total_gb:.2f} GiB",
    )

    print("\nFile più grandi da includere:")

    for path, size in sorted(
        files,
        key=lambda x: x[1],
        reverse=True,
    )[:12]:

        print(
            f"{size / 1024**3:.2f} GiB |",
            path.relative_to(ROOT),
        )

    if total_gb > args.max_gb:

        raise RuntimeError(
            f"Pacchetto troppo grande: {total_gb:.2f} GiB. "
            f"Limite: {args.max_gb:.2f} GiB."
        )

    if not args.build:

        print("\n[OK] CONTROLLO PRELIMINARE COMPLETATO")
        print("Nessun ZIP creato.")
        print("Per creare lo ZIP aggiungere --build.")
        return

    # ======================================================
    # CREAZIONE ARCHIVIO
    # ======================================================

    stamp = dt.datetime.now(
        dt.timezone.utc
    ).strftime("%Y%m%dT%H%M%SZ")

    destination = (
        OUTDIR
        / f"{PACK_NAME}_{stamp}.zip"
    )

    if destination.exists():
        raise RuntimeError(
            f"Archivio già presente: {destination}"
        )

    manifest = {
        "created_utc": stamp,
        "project": "AI PRIMETTA V3.3",
        "images": images,
        "project_files": len(files),
        "uncompressed_bytes": total,
        "excluded_raw_reads": True,
        "excluded_paths": skipped,
        "windows_tested": False,
    }

    try:

        with tempfile.TemporaryDirectory(
            prefix="primetta_light_",
            dir=OUTDIR,
        ) as temporary:

            tmp = Path(temporary)

            images_archive = (
                tmp / "docker_images.tar.gz"
            )

            print(
                "\n[3] ESPORTAZIONE IMMAGINI DOCKER",
                flush=True,
            )

            export_images(
                images,
                images_archive,
            )

            print(
                "[OK] Immagini Docker:",
                round(
                    images_archive.stat().st_size
                    / (1024 ** 3),
                    2,
                ),
                "GiB",
            )

            print(
                "\n[4] CREAZIONE ZIP",
                flush=True,
            )

            with zipfile.ZipFile(
                destination,
                "w",
                allowZip64=True,
            ) as archive:

                prefix = PACK_NAME + "/"

                for key, config in configs.items():

                    archive.writestr(
                        prefix + f"compose_{key}.yaml",
                        json.dumps(
                            config,
                            indent=2,
                            ensure_ascii=False,
                        ) + "\n",
                    )

                archive.writestr(
                    prefix + "config_portatile.json",
                    json.dumps(
                        {
                            "images": images,
                            "services": services,
                        },
                        indent=2,
                    ),
                )

                archive.writestr(
                    prefix + "gestione_portatile.ps1",
                    POWERSHELL,
                )

                bat_files = [
                    (
                        "01_AVVIA_GUI_ORIGINALE_V32.bat",
                        "v32",
                    ),
                    (
                        "02_AVVIA_GUI_PREDIZIONE_V33.bat",
                        "v33",
                    ),
                    (
                        "03_FERMA_ENTRAMBE_LE_GUI.bat",
                        "stop",
                    ),
                ]

                for filename, action in bat_files:

                    archive.writestr(
                        prefix + filename,
                        BAT.replace(
                            "__ACTION__",
                            action,
                        ),
                    )

                archive.writestr(
                    prefix + "ISTRUZIONI.txt",
                    README,
                )

                archive.writestr(
                    prefix + "MANIFEST_EXPORT.json",
                    json.dumps(
                        manifest,
                        indent=2,
                        ensure_ascii=False,
                    ),
                )

                excluded_text = (
                    "percorso\tmotivo\n"
                    + "\n".join(
                        f"{path}\t{reason}"
                        for path, reason in skipped
                    )
                )

                archive.writestr(
                    prefix + "FILE_ESCLUSI.tsv",
                    excluded_text,
                )

                archive.write(
                    images_archive,
                    prefix + "images/docker_images.tar.gz",
                    compress_type=zipfile.ZIP_STORED,
                )

                for index, (path, size) in enumerate(
                    files,
                    start=1,
                ):

                    relative = (
                        path.relative_to(ROOT)
                        .as_posix()
                    )

                    arcname = (
                        prefix
                        + "progetto/"
                        + relative
                    )

                    if path.name.lower().endswith(
                        ZIP_COMPRESSED
                    ):
                        compression = zipfile.ZIP_STORED
                    else:
                        compression = zipfile.ZIP_DEFLATED

                    archive.write(
                        path,
                        arcname=arcname,
                        compress_type=compression,
                        compresslevel=1,
                    )

                    if index % 300 == 0:
                        print(
                            f"File aggiunti: {index}/{len(files)}",
                            flush=True,
                        )

        # ==================================================
        # VERIFICA ZIP
        # ==================================================

        print(
            "\n[5] VERIFICA INTEGRITÀ ZIP",
            flush=True,
        )

        with zipfile.ZipFile(
            destination
        ) as archive:

            bad_file = archive.testzip()

            if bad_file:
                raise RuntimeError(
                    f"ZIP non valido: {bad_file}"
                )

        print()
        print("=" * 70)
        print("STEP 21.1 COMPLETATO")
        print("=" * 70)

        print("[OK] ZIP creato")
        print("[OK] Integrità ZIP verificata")
        print("[OK] FASTQ grezzi esclusi")
        print("[OK] GUI originali preservate")
        print("[OK] Nessun training")

        print("\nPERCORSO:")
        print(destination)

        print(
            "\nDIMENSIONE:",
            round(
                destination.stat().st_size / (1024 ** 3),
                2,
            ),
            "GiB",
        )

        print("\nSHA256:")
        print(sha256(destination))

        print(
            "\n[NOTA] Verificare entrambi gli avvii "
            "sul nuovo PC Windows."
        )

    except Exception:

        if destination.exists():

            try:
                destination.unlink()

                print(
                    "[INFO] ZIP incompleto eliminato"
                )

            except OSError:
                pass

        raise


# ==========================================================
# AVVIO
# ==========================================================

if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print("\n[INTERROTTO] Procedura annullata")
        sys.exit(130)

    except Exception as exc:
        print(
            "[ERRORE]",
            exc,
            file=sys.stderr,
        )
        sys.exit(1)