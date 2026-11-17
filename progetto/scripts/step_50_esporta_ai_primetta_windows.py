
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 - STEP 50 V2
Esportazione Windows offline con ZIP64.

--audit : solo lettura
--build : crea ZIP con una sola copia del progetto

Riutilizza le immagini Docker dello STEP 50 precedente.
Non riaddestra modelli e non modifica le GUI.
"""

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
import zipfile

from collections import defaultdict
from pathlib import Path


PROJECT = Path(__file__).resolve().parent.parent

OLD = (
    Path.home()
    / "AI_PRIMETTA_WINDOWS_EXPORT"
    / "AI_PRIMETTA_V33_WINDOWS_PORTABLE_20261024T160054Z"
    / "images"
    / "ai_primetta_images.tar"
)

OUTPUT = Path.home() / "AI_PRIMETTA_WINDOWS_EXPORT"

CONTAINERS = {
    "main": "maize_ncrna_gnn",
    "prediction": "ai_primetta_target_v33-prediction-1",
    "monitor": "ai_primetta_v33_ml_monitor-monitor-1",
}

TAGS = {
    "main": "ai-primetta-portable-main:20261024t160054z",
    "prediction": "ai-primetta-portable-prediction:20261024t160054z",
    "monitor": "ai-primetta-portable-monitor:20261024t160054z",
}

PORTS = {
    "main": {8000, 8503},
    "prediction": {8502},
    "monitor": set(),
}

SKIP_DIRS = {
    ".venv",
    ".git",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    ".cache",
    "node_modules",
}

RUNTIME_EXCLUSIONS = [
    ("data", "public_ena"),
    ("data", "raw"),
    ("data", "downloads"),
    ("delivery",),
    ("outputs", "final_package"),
    ("outputs", "v2_package"),
]

SKIP_SUFFIXES = (
    ".pyc",
    ".swp",
    ".tmp",
    ".part",
)

PRIVATE_NAMES = {
    ".env",
    "id_rsa",
    "id_ed25519",
    "credentials.json",
    "secrets.json",
}

PRIVATE_SUFFIXES = {
    ".pem",
    ".p12",
    ".pfx",
    ".key",
}


# ============================================================
# 1. COMANDI DI SISTEMA
# ============================================================

def cmd(*items, check=True):

    process = subprocess.run(
        [str(x) for x in items],
        capture_output=True,
        text=True,
    )

    if check and process.returncode:
        raise RuntimeError(
            "Comando fallito: "
            + " ".join(map(str, items))
            + "\n"
            + process.stderr[-600:]
        )

    return process


# ============================================================
# 2. VERIFICA CONTAINER
# ============================================================

def inspect_containers():

    process = cmd(
        "docker",
        "inspect",
        *CONTAINERS.values(),
    )

    found = {
        c["Name"].lstrip("/"): c
        for c in json.loads(process.stdout)
    }

    configs = {}

    for service, name in CONTAINERS.items():

        if (
            name not in found
            or not found[name]["State"]["Running"]
        ):
            raise RuntimeError(
                "Container assente o fermo: " + name
            )

        configs[service] = found[name]

    return configs


# ============================================================
# 3. VERIFICA ARCHIVIO IMMAGINI ESISTENTE
# ============================================================

def verify_images_and_tar(tar_path):

    if (
        not tar_path.is_file()
        or tar_path.stat().st_size < 1_000_000
    ):
        raise RuntimeError(
            "Archivio immagini mancante o troppo piccolo: "
            + str(tar_path)
        )

    for image in TAGS.values():
        cmd(
            "docker",
            "image",
            "inspect",
            image,
        )

    try:
        with tarfile.open(tar_path, mode="r:") as archive:

            member = archive.getmember(
                "manifest.json"
            )

            stream = archive.extractfile(member)

            if stream is None:
                raise RuntimeError(
                    "Manifest TAR non leggibile"
                )

            manifest = json.load(stream)

            tags = {
                tag
                for entry in manifest
                for tag in entry.get("RepoTags", [])
            }

    except (
        tarfile.TarError,
        KeyError,
        ValueError,
    ) as exc:

        raise RuntimeError(
            "Archivio Docker non valido: " + str(exc)
        ) from exc

    missing = set(TAGS.values()) - tags

    if missing:
        raise RuntimeError(
            "Immagini mancanti nel TAR: "
            + str(sorted(missing))
        )

    return tags


# ============================================================
# 4. COMPOSE PORTATILE
# ============================================================

def compose_from_containers(configs):

    result = {
        "name": "ai_primetta_portable",
        "services": {},
        "networks": {
            "primetta_main": {},
            "primetta_prediction": {},
        },
    }

    audit_sources = {}

    for service, config in configs.items():

        ports = {}

        entries = (
            config.get("NetworkSettings", {})
            .get("Ports") or {}
        )

        for inner, hosts in entries.items():

            if (
                not hosts
                or not inner.endswith("/tcp")
            ):
                continue

            inside = int(
                inner.split("/")[0]
            )

            for host in hosts:

                outside = int(
                    host["HostPort"]
                )

                ports[outside] = (
                    f"127.0.0.1:{outside}:{inside}"
                )

        if set(ports) != PORTS[service]:
            raise RuntimeError(
                f"Porte inattese di {service}: "
                f"{sorted(ports)}"
            )

        host_config = config.get(
            "HostConfig"
        ) or {}

        if (
            host_config.get("Privileged")
            or host_config.get("Devices")
            or host_config.get("NetworkMode") == "host"
        ):
            raise RuntimeError(
                "Configurazione Docker non portabile: "
                + service
            )

        svc = {
            "image": TAGS[service],
            "restart": "unless-stopped",
        }

        if ports:
            svc["ports"] = list(
                ports.values()
            )

        if service == "monitor":
            svc["network_mode"] = "none"
        else:
            svc["networks"] = [
                (
                    "primetta_main"
                    if service == "main"
                    else "primetta_prediction"
                )
            ]

        volumes = []

        for mount in config.get("Mounts", []):

            if mount.get("Type") != "bind":
                raise RuntimeError(
                    "Mount non gestito: "
                    + str(mount)
                )

            source_path = Path(
                mount["Source"]
            ).resolve()

            try:
                relative = source_path.relative_to(
                    PROJECT
                )

            except ValueError as exc:
                raise RuntimeError(
                    "Bind esterno al progetto: "
                    + str(source_path)
                ) from exc

            if not source_path.exists():
                raise RuntimeError(
                    "Origine bind inesistente: "
                    + str(source_path)
                )

            if any(
                p in SKIP_DIRS
                for p in relative.parts
            ):
                raise RuntimeError(
                    "Bind in directory esclusa: "
                    + str(source_path)
                )

            if any(
                relative.parts[:len(x)] == x
                for x in RUNTIME_EXCLUSIONS
            ):
                raise RuntimeError(
                    "Bind richiesto in directory raw esclusa: "
                    + str(source_path)
                )

            source = "./project"

            if relative.parts:
                source += "/" + relative.as_posix()

            volumes.append({
                "type": "bind",
                "source": source,
                "target": mount["Destination"],
                "read_only": not mount.get(
                    "RW", True
                ),
            })

            audit_sources[
                relative.as_posix()
            ] = (
                "dir"
                if source_path.is_dir()
                else "file"
            )

        svc["volumes"] = volumes
        result["services"][service] = svc

    return result, audit_sources


# ============================================================
# 5. ESCLUSIONI
# ============================================================

def is_skipped(relative, profile):

    if any(
        part in SKIP_DIRS
        for part in relative.parts
    ):
        return True

    if profile == "runtime":

        if any(
            relative.parts[:len(x)] == x
            for x in RUNTIME_EXCLUSIONS
        ):
            return True

    if relative.name.lower().endswith(
        SKIP_SUFFIXES
    ):
        return True

    return False


# ============================================================
# 6. SCANSIONE PROGETTO SENZA COPIE
# ============================================================

def scan_project(profile):

    records = []
    directories = set()
    skipped_links = []
    sizes = defaultdict(int)

    total = 0

    for base, dirnames, filenames in os.walk(
        PROJECT,
        topdown=True,
        followlinks=False,
    ):

        folder = Path(base)
        relative_dir = folder.relative_to(PROJECT)

        kept = []

        for dirname in dirnames:

            path = folder / dirname
            relative = path.relative_to(PROJECT)

            if is_skipped(relative, profile):
                continue

            if path.is_symlink():
                skipped_links.append(
                    relative.as_posix()
                )
                continue

            kept.append(dirname)

        dirnames[:] = kept

        if relative_dir.parts:
            directories.add(
                relative_dir.as_posix()
            )

        for filename in filenames:

            path = folder / filename
            relative = path.relative_to(PROJECT)

            if is_skipped(relative, profile):
                continue

            if path.is_symlink():
                skipped_links.append(
                    relative.as_posix()
                )
                continue

            try:
                stat = path.stat()

            except OSError as exc:
                raise RuntimeError(
                    "Impossibile leggere: "
                    + str(path)
                ) from exc

            if not path.is_file():
                continue

            records.append((
                path,
                relative.as_posix(),
                stat.st_size,
            ))

            total += stat.st_size

            first = (
                relative.parts[0]
                if relative.parts
                else "root"
            )

            sizes[first] += stat.st_size

    return (
        records,
        directories,
        skipped_links,
        sizes,
        total,
    )


# ============================================================
# 7. INTEGRITA DEI MOUNT
# ============================================================

def verify_mounted_sources(
    audit_sources,
    records,
    directories,
    profile,
):

    files = {
        relative
        for _, relative, _ in records
    }

    for relative, kind in audit_sources.items():

        if relative == ".":
            continue

        if is_skipped(
            Path(relative), profile
        ):
            raise RuntimeError(
                "Mount necessario escluso: "
                + relative
            )

        if kind == "file":

            if relative not in files:
                raise RuntimeError(
                    "File bind non incluso: "
                    + relative
                )

        if kind == "dir":

            if relative not in directories:
                raise RuntimeError(
                    "Directory bind non inclusa: "
                    + relative
                )


# ============================================================
# 8. HASH SHA256
# ============================================================

def sha256(path):

    digest = hashlib.sha256()

    with Path(path).open("rb") as handle:

        for block in iter(
            lambda: handle.read(
                4 * 1024 * 1024
            ),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


# ============================================================
# 9. FILE BAT WINDOWS
# ============================================================

def bats():

    images = " ".join(
        '"' + tag + '"'
        for tag in TAGS.values()
    )

    common = (
        '@echo off\r\n'
        'setlocal EnableExtensions\r\n'
        'cd /d "%~dp0"\r\n'
        'where docker >nul 2>&1\r\n'
        'if errorlevel 1 (\r\n'
        '  echo [ERRORE] Docker Desktop non trovato.\r\n'
        '  pause\r\n'
        '  exit /b 1\r\n'
        ')\r\n'
    )

    start = common + (
        'docker info >nul 2>&1\r\n'
        'if errorlevel 1 (\r\n'
        '  if exist "%ProgramFiles%\\Docker\\Docker\\Docker Desktop.exe" (\r\n'
        '    start "" "%ProgramFiles%\\Docker\\Docker\\Docker Desktop.exe"\r\n'
        '  )\r\n'
        '  for /l %%I in (1,1,90) do (\r\n'
        '    docker info >nul 2>&1 && goto READY\r\n'
        '    timeout /t 2 /nobreak >nul\r\n'
        '  )\r\n'
        '  echo Docker non disponibile.\r\n'
        '  pause\r\n'
        '  exit /b 1\r\n'
        ')\r\n'
        ':READY\r\n'
        'set "NEEDLOAD="\r\n'
        'for %%G in (' + images + ') do (\r\n'
        '  docker image inspect "%%~G" >nul 2>&1 || set "NEEDLOAD=1"\r\n'
        ')\r\n'
        'if defined NEEDLOAD (\r\n'
        '  echo Importazione immagini Docker offline...\r\n'
        '  docker load -i "%~dp0images\\ai_primetta_images.tar"\r\n'
        '  if errorlevel 1 goto FAIL\r\n'
        ')\r\n'
        'docker compose -p ai_primetta_portable -f "%~dp0compose.portable.json" up -d --pull never\r\n'
        'if errorlevel 1 goto FAIL\r\n'
        'echo Avvio completato.\r\n'
        'start "" "http://127.0.0.1:8503/"\r\n'
        'start "" "http://127.0.0.1:8502/"\r\n'
        'echo API: http://127.0.0.1:8000/docs\r\n'
        'pause\r\n'
        'exit /b 0\r\n'
        ':FAIL\r\n'
        'echo Errore di avvio.\r\n'
        'pause\r\n'
        'exit /b 1\r\n'
    )

    stop = common + (
        'docker compose -p ai_primetta_portable -f "%~dp0compose.portable.json" stop\r\n'
        'if errorlevel 1 (\r\n'
        '  echo Errore arresto.\r\n'
        '  pause\r\n'
        '  exit /b 1\r\n'
        ')\r\n'
        'echo Solo AI PRIMETTA arrestata.\r\n'
        'echo Dati conservati.\r\n'
        'pause\r\n'
    )

    status = common + (
        'docker compose -p ai_primetta_portable -f "%~dp0compose.portable.json" ps\r\n'
        'echo GUI principale: http://127.0.0.1:8503/\r\n'
        'echo GUI predittiva: http://127.0.0.1:8502/\r\n'
        'echo API: http://127.0.0.1:8000/docs\r\n'
        'pause\r\n'
    )

    return {
        "AVVIA_AI_PRIMETTA.bat": start,
        "STOP_AI_PRIMETTA.bat": stop,
        "STATO_AI_PRIMETTA.bat": status,
    }


# ============================================================
# 10. DOCUMENTAZIONE WINDOWS
# ============================================================

def readme(profile):

    return (
        "AI PRIMETTA V3.3 | WINDOWS OFFLINE\n\n"
        "REQUISITI:\n"
        "Windows 10/11 a 64 bit.\n"
        "Docker Desktop con WSL 2 e Linux containers.\n"
        "Disco NTFS e spazio sufficiente.\n\n"

        "PROCEDURA:\n"
        "1. Estrarre lo ZIP su C:\\AI_PRIMETTA_V33.\n"
        "2. Avviare Docker Desktop.\n"
        "3. Eseguire AVVIA_AI_PRIMETTA.bat.\n"
        "4. Aprire le due GUI.\n"
        "5. Arrestare con STOP_AI_PRIMETTA.bat.\n\n"

        "GUI PRINCIPALE: http://127.0.0.1:8503\n"
        "GUI PREDITTIVA: http://127.0.0.1:8502\n"
        "API: http://127.0.0.1:8000/docs\n\n"

        f"PROFILO: {profile}\n\n"

        "Il profilo runtime esclude i dati grezzi ENA,\n"
        "le cache e le copie duplicate delle consegne.\n"
        "Non e un archivio scientifico integrale.\n\n"

        "ATTENZIONE:\n"
        "Il progetto e esplorativo.\n"
        "I target biologici non sono validati.\n"
        "Verificare credenziali, dati riservati e licenze\n"
        "prima di condividere il pacchetto.\n\n"

        "La portabilita deve essere collaudata\n"
        "su un PC Windows differente.\n"
    )


# ============================================================
# 11. GENERAZIONE ZIP DIRETTA
# ============================================================

def write_zip(
    output,
    old_tar,
    compose,
    records,
    directories,
    profile,
    audit_sources,
    total,
):

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if output.exists():
        raise RuntimeError(
            "Archivio gia presente: "
            + str(output)
        )

    created = dt.datetime.now(
        dt.timezone.utc
    ).isoformat()

    tar_hash = sha256(old_tar)

    manifest = {
        "project": "AI PRIMETTA V3.3",
        "profile": profile,
        "created_utc": created,
        "image_tags": TAGS,
        "image_tar_sha256": tar_hash,
        "files_in_project": len(records),
        "uncompressed_project_bytes": total,
        "unique_bind_sources": len(audit_sources),
        "source_containers": CONTAINERS,
        "scientific_status": (
            "EXPLORATORY_NOT_BIOLOGICALLY_VALIDATED"
        ),
    }

    started = time.monotonic()
    written = 0

    with zipfile.ZipFile(
        output,
        "w",
        allowZip64=True,
    ) as archive:

        archive.writestr(
            "compose.portable.json",
            json.dumps(
                compose,
                indent=2
            ) + "\n",
        )

        archive.writestr(
            "MANIFEST_PORTABILITA.json",
            json.dumps(
                manifest,
                indent=2
            ) + "\n",
        )

        archive.writestr(
            "LEGGIMI_WINDOWS.txt",
            readme(profile),
        )

        for name, content in bats().items():

            archive.writestr(
                name,
                content.encode("ascii"),
                compress_type=zipfile.ZIP_STORED,
            )

        archive.write(
            old_tar,
            "images/ai_primetta_images.tar",
            compress_type=zipfile.ZIP_STORED,
        )

        print(
            "[IMMAGINI] Archivio esistente incluso "
            "una sola volta",
            flush=True,
        )

        for directory in sorted(directories):

            archive.writestr(
                "project/"
                + directory.rstrip("/")
                + "/",
                b"",
                compress_type=zipfile.ZIP_STORED,
            )

        for index, (
            source,
            relative,
            length,
        ) in enumerate(records, 1):

            suffix = source.suffix.lower()

            compression = (
                zipfile.ZIP_STORED
                if suffix in {
                    ".gz",
                    ".zip",
                    ".tar",
                    ".pt",
                    ".pth",
                    ".joblib",
                    ".npz",
                    ".parquet",
                    ".jpg",
                    ".png",
                }
                else zipfile.ZIP_DEFLATED
            )

            archive.write(
                source,
                "project/" + relative,
                compress_type=compression,
                compresslevel=(
                    1
                    if compression == zipfile.ZIP_DEFLATED
                    else None
                ),
            )

            written += length

            if index % 1000 == 0:

                print(
                    f"[AVANZAMENTO] "
                    f"{index}/{len(records)} file | "
                    f"{written / 1024**3:.2f}/"
                    f"{total / 1024**3:.2f} GiB",
                    flush=True,
                )

                time.sleep(0.05)

    print("[ZIP_OK]", output)

    print(
        "[ZIP_GIB]",
        round(
            output.stat().st_size / 1024**3,
            3,
        ),
    )

    print(
        "[ZIP_SHA256]",
        sha256(output),
    )

    print(
        "[TEMPO_SECONDI]",
        round(
            time.monotonic() - started,
            1,
        ),
    )

    print(
        "[NOTA] Collaudo obbligatorio su Windows"
    )


# ============================================================
# 12. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = parser.add_mutually_exclusive_group(
        required=True
    )

    mode.add_argument(
        "--audit",
        action="store_true",
    )

    mode.add_argument(
        "--build",
        action="store_true",
    )

    parser.add_argument(
        "--profile",
        choices=[
            "runtime",
            "full",
        ],
        default="runtime",
    )

    parser.add_argument(
        "--max-gb",
        type=float,
        default=12.0,
    )

    parser.add_argument(
        "--allow-private",
        action="store_true",
    )

    parser.add_argument(
        "--image-tar",
        type=Path,
        default=OLD,
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT,
    )

    args = parser.parse_args()

    print("=" * 70)

    print(
        "AI PRIMETTA V3.3 | STEP50 V2 |",
        "AUDIT" if args.audit else "BUILD",
    )

    print("[PROFILO]", args.profile)
    print("[PROJECT]", PROJECT)
    print("=" * 70, flush=True)

    configs = inspect_containers()

    verify_images_and_tar(
        args.image_tar
    )

    compose, sources = compose_from_containers(
        configs
    )

    print(
        "[CONTAINER]",
        len(configs),
    )

    print(
        "[MOUNT_TOTALI]",
        sum(
            len(config["Mounts"])
            for config in configs.values()
        ),
    )

    print(
        "[SORGENTI_BIND_UNICHE]",
        len(sources),
    )

    print(
        "[ARCHIVIO_IMMAGINI_MIB]",
        round(
            args.image_tar.stat().st_size / 1024**2,
            1,
        ),
    )

    (
        records,
        directories,
        links,
        sizes,
        total,
    ) = scan_project(
        args.profile
    )

    verify_mounted_sources(
        sources,
        records,
        directories,
        args.profile,
    )

    print(
        "[FILE_DA_INSERIRE]",
        len(records),
    )

    print(
        "[DIRECTORY]",
        len(directories),
    )

    print(
        "[DIMENSIONE_DATI_GIB]",
        round(
            total / 1024**3,
            3,
        ),
    )

    print(
        "[LINK_SIMBOLICI_ESCLUSI]",
        len(links),
    )

    for name, value in sorted(
        sizes.items(),
        key=lambda item: -item[1],
    )[:12]:

        print(
            f"[AREA] {name}: "
            f"{value / 1024**3:.3f} GiB"
        )

    private = [
        relative
        for _, relative, _ in records
        if (
            Path(relative).name.lower()
            in PRIVATE_NAMES
            or Path(relative).suffix.lower()
            in PRIVATE_SUFFIXES
        )
    ]

    print(
        "[FILE_POTENZIALMENTE_RISERVATI]",
        len(private),
    )

    for relative in private[:12]:
        print(
            "[RISERVATO_CONTROLLARE]",
            relative,
        )

    print(
        "[NOTA] Una sola copia del progetto"
    )

    print(
        "[NOTA] Nessun nuovo Docker commit"
    )

    if args.audit:

        print(
            "[AUDIT_OK] Nessun file modificato"
        )

        return

    if private and not args.allow_private:

        raise RuntimeError(
            "File potenzialmente riservati presenti. "
            "Verificare e usare --allow-private "
            "solo per distribuzione autorizzata."
        )

    if links:

        print(
            "[ATTENZIONE] Symlink esclusi:",
            links[:10],
        )

    limit = args.max_gb * 1024**3

    if total > limit:

        raise RuntimeError(
            f"Dati {total / 1024**3:.2f} GiB "
            f"oltre il limite {args.max_gb} GiB. "
            "Verificare prima di aumentarlo."
        )

    output_dir = args.output_dir.expanduser().resolve()

    if (
        output_dir == PROJECT
        or PROJECT in output_dir.parents
    ):
        raise RuntimeError(
            "Destinazione dentro il progetto "
            "non consentita"
        )

    free = shutil.disk_usage(
        output_dir.parent
    ).free

    required = (
        total
        + args.image_tar.stat().st_size
        + 2 * 1024**3
    )

    if free < required:

        raise RuntimeError(
            "Spazio disco insufficiente"
        )

    stamp = dt.datetime.now(
        dt.timezone.utc
    ).strftime(
        "%Y%m%dT%H%M%SZ"
    )

    output = output_dir / (
        "AI_PRIMETTA_V33_WINDOWS_PORTABLE_V2_"
        + stamp
        + ".zip"
    )

    write_zip(
        output,
        args.image_tar,
        compose,
        records,
        directories,
        args.profile,
        sources,
        total,
    )


if __name__ == "__main__":

    try:
        main()

    except (
        OSError,
        RuntimeError,
        subprocess.SubprocessError,
        zipfile.BadZipFile,
    ) as exc:

        print(
            "[STOP]",
            exc,
            file=sys.stderr,
        )

        sys.exit(1)
