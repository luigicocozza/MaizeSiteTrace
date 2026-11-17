#!/usr/bin/env python3
"""
00_download_public_data.py

Scarica/interroga le fonti pubbliche configurate in configs/config_public.yaml
per il progetto maize_drought_ncrna_gnn.

Output principali:
- data/raw/*
- logs/download_public_data.log
- outputs/downloads_report.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import urlparse

import requests
import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scarica/interroga fonti pubbliche per il progetto ncRNA mais/siccità."
    )
    parser.add_argument(
        "--config",
        default="configs/config_public.yaml",
        help="Percorso del file YAML di configurazione.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Riscarica anche i file già presenti.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="Timeout in secondi per ogni richiesta HTTP.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Abilita logging più dettagliato.",
    )
    parser.add_argument(
        "--max-records",
        type=int,
        default=50,
        help="Numero massimo di record per API che supportano paginazione/retmax.",
    )
    return parser.parse_args()


def load_yaml(config_path: Path) -> Dict[str, Any]:
    if not config_path.exists():
        raise FileNotFoundError(f"File di configurazione non trovato: {config_path}")

    try:
        with config_path.open("r", encoding="utf-8") as f:
            cfg = yaml.safe_load(f)
    except yaml.YAMLError as exc:
        raise ValueError(f"YAML malformato: {exc}") from exc

    if not isinstance(cfg, dict):
        raise ValueError("Il file YAML non contiene un dizionario valido.")

    return cfg


def validate_config(cfg: Dict[str, Any]) -> None:
    required_sections = [
        "project",
        "paths",
        "runtime",
        "scope",
        "tools",
        "logging",
        "sources",
        "reference_databases",
        "training_data_strategy",
        "preprocessing",
        "labeling",
        "model",
        "pipeline_outputs",
        "quality_control",
        "reproducibility",
        "publication_strategy",
    ]

    missing = [section for section in required_sections if section not in cfg]
    if missing:
        raise ValueError(f"Sezioni obbligatorie mancanti nel YAML: {missing}")

    if not isinstance(cfg.get("sources"), list):
        raise ValueError("La sezione 'sources' deve essere una lista.")

    species = cfg.get("scope", {}).get("species", {})
    if species.get("scientific_name") != "Zea mays":
        logging.warning(
            "Attenzione: la specie scientifica non risulta 'Zea mays'. Valore trovato: %s",
            species.get("scientific_name"),
        )


def ensure_directories(cfg: Dict[str, Any]) -> Dict[str, Path]:
    paths_cfg = cfg.get("paths", {})

    directories = {
        "raw_data": Path(paths_cfg.get("raw_data", "data/raw")),
        "interim_data": Path(paths_cfg.get("interim_data", "data/interim")),
        "processed_data": Path(paths_cfg.get("processed_data", "data/processed")),
        "logs": Path(paths_cfg.get("logs", "logs")),
        "outputs": Path(paths_cfg.get("outputs", "outputs")),
        "models": Path(paths_cfg.get("models", "models")),
        "figures": Path(paths_cfg.get("figures", "outputs/figures")),
        "tables": Path(paths_cfg.get("tables", "outputs/tables")),
        "reports": Path(paths_cfg.get("reports", "outputs/reports")),
    }

    for path in directories.values():
        path.mkdir(parents=True, exist_ok=True)

    return directories


def setup_logging(log_dir: Path, verbose: bool) -> None:
    log_file = log_dir / "download_public_data.log"
    level = logging.DEBUG if verbose else logging.INFO

    handlers = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(log_file, mode="a", encoding="utf-8"),
    ]

    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_valid_url(url: str) -> bool:
    parsed = urlparse(url)
    return parsed.scheme in {"http", "https"} and bool(parsed.netloc)


def sanitize_filename(name: str) -> str:
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in name.strip())
    return safe or "downloaded_source"


def build_request(
    source: Dict[str, Any],
    cfg: Dict[str, Any],
    max_records: int,
) -> Tuple[str, Dict[str, Any]]:
    """
    Costruisce URL e parametri.
    Per alcune fonti API aggiunge query biologicamente utili.
    Per fonti web semplici scarica la pagina o metadata grezzi.
    """
    name = source.get("name", "")
    url = source.get("url", "")
    params: Dict[str, Any] = {}

    taxon_id = cfg.get("scope", {}).get("species", {}).get("ncbi_taxon_id", 4577)

    geo_query = (
        cfg.get("reference_databases", {})
        .get("ncbi_geo_sra", {})
        .get("search_query", '("Zea mays" OR maize) AND drought')
    )

    literature_query = (
        cfg.get("reference_databases", {})
        .get("europe_pmc", {})
        .get(
            "search_query",
            '("Zea mays" OR maize) AND drought AND (ncRNA OR lncRNA OR miRNA)',
        )
    )

    lowered = name.lower()

    if lowered == "ncbi_geo":
        params = {
            "db": "gds",
            "term": geo_query,
            "retmode": "json",
            "retmax": max_records,
        }

    elif lowered == "europe_pmc":
        params = {
            "query": literature_query,
            "format": "json",
            "pageSize": min(max_records, 100),
        }

    elif lowered == "quickgo":
        params = {
            "taxonId": taxon_id,
            "limit": min(max_records, 100),
        }

    elif lowered == "rnacentral":
        params = {
            "taxid": taxon_id,
            "page_size": min(max_records, 100),
            "format": "json",
        }

    return url, params


def download_source(
    source: Dict[str, Any],
    cfg: Dict[str, Any],
    output_dir: Path,
    force_download: bool,
    timeout: int,
    max_records: int,
) -> Dict[str, Any]:
    source_name = source.get("name", "unknown_source")
    enabled = bool(source.get("enabled", False))
    output_filename = source.get("output_filename") or f"{sanitize_filename(source_name)}.dat"
    output_path = output_dir / output_filename
    tmp_path = output_path.with_suffix(output_path.suffix + ".part")

    result: Dict[str, Any] = {
        "source_name": source_name,
        "enabled": enabled,
        "category": source.get("category"),
        "access_type": source.get("access_type"),
        "url": source.get("url"),
        "output_file": str(output_path),
        "status": None,
        "http_status_code": None,
        "file_size_bytes": None,
        "sha256": None,
        "error": None,
        "request_params": {},
    }

    if not enabled:
        result["status"] = "disabled"
        return result

    url = source.get("url")
    if not url or not isinstance(url, str):
        result["status"] = "failed"
        result["error"] = "URL mancante o non valida nel YAML."
        return result

    if not is_valid_url(url):
        result["status"] = "failed"
        result["error"] = f"URL non valida: {url}"
        return result

    if output_path.exists() and not force_download:
        result["status"] = "skipped_existing"
        result["file_size_bytes"] = output_path.stat().st_size
        result["sha256"] = sha256_file(output_path)
        return result

    request_url, params = build_request(source, cfg, max_records=max_records)
    result["request_params"] = params

    headers = {
        "User-Agent": "maize-drought-ncrna-gnn/0.1.0 (+research; public data pipeline)"
    }

    logging.info("Download/interrogazione fonte: %s", source_name)
    logging.info("URL: %s", request_url)
    if params:
        logging.info("Parametri: %s", params)

    try:
        with requests.get(
            request_url,
            params=params,
            timeout=timeout,
            headers=headers,
            stream=True,
        ) as response:
            result["http_status_code"] = response.status_code
            response.raise_for_status()

            with tmp_path.open("wb") as f:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)

        shutil.move(str(tmp_path), str(output_path))

        result["status"] = "downloaded"
        result["file_size_bytes"] = output_path.stat().st_size
        result["sha256"] = sha256_file(output_path)

    except requests.exceptions.Timeout:
        result["status"] = "failed"
        result["error"] = f"Timeout dopo {timeout} secondi."
        if tmp_path.exists():
            tmp_path.unlink()

    except requests.exceptions.RequestException as exc:
        result["status"] = "failed"
        result["error"] = f"Errore HTTP/rete: {exc}"
        if tmp_path.exists():
            tmp_path.unlink()

    except Exception as exc:
        result["status"] = "failed"
        result["error"] = f"Errore inatteso: {exc}"
        if tmp_path.exists():
            tmp_path.unlink()

    return result


def save_config_snapshot(config_path: Path, outputs_dir: Path) -> Optional[str]:
    try:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        snapshot_path = outputs_dir / f"config_snapshot_{timestamp}.yaml"
        shutil.copy2(config_path, snapshot_path)
        return str(snapshot_path)
    except Exception as exc:
        logging.warning("Impossibile salvare snapshot del config: %s", exc)
        return None


def write_report(report: Dict[str, Any], report_path: Path) -> None:
    with report_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)


def main() -> int:
    args = parse_args()
    config_path = Path(args.config)

    start_time = time.time()
    timestamp_utc = datetime.now(timezone.utc).isoformat()

    try:
        cfg = load_yaml(config_path)
        validate_config(cfg)
        directories = ensure_directories(cfg)
        setup_logging(directories["logs"], args.verbose)

        logging.info("Avvio download dati pubblici")
        logging.info("Config: %s", config_path)

        config_snapshot = save_config_snapshot(config_path, directories["outputs"])

        source_results: List[Dict[str, Any]] = []
        for source in cfg.get("sources", []):
            result = download_source(
                source=source,
                cfg=cfg,
                output_dir=directories["raw_data"],
                force_download=args.force_download,
                timeout=args.timeout,
                max_records=args.max_records,
            )
            source_results.append(result)

            if result["status"] == "downloaded":
                logging.info(
                    "OK: %s -> %s (%s bytes)",
                    result["source_name"],
                    result["output_file"],
                    result["file_size_bytes"],
                )
            elif result["status"] == "skipped_existing":
                logging.info(
                    "SKIP: %s già presente -> %s",
                    result["source_name"],
                    result["output_file"],
                )
            elif result["status"] == "disabled":
                logging.info("DISABLED: %s", result["source_name"])
            else:
                logging.error(
                    "FAILED: %s | %s",
                    result["source_name"],
                    result.get("error"),
                )

        downloaded = [r for r in source_results if r["status"] == "downloaded"]
        skipped = [r for r in source_results if r["status"] == "skipped_existing"]
        failed = [r for r in source_results if r["status"] == "failed"]
        disabled = [r for r in source_results if r["status"] == "disabled"]

        final_status = "success" if not failed else "completed_with_errors"

        report = {
            "timestamp_utc": timestamp_utc,
            "config_used": str(config_path),
            "config_snapshot": config_snapshot,
            "project": cfg.get("project", {}),
            "species": cfg.get("scope", {}).get("species", {}),
            "biological_focus": cfg.get("scope", {}).get("biological_focus", {}),
            "parameters": {
                "force_download": args.force_download,
                "timeout": args.timeout,
                "verbose": args.verbose,
                "max_records": args.max_records,
            },
            "summary": {
                "final_status": final_status,
                "n_sources_total": len(source_results),
                "n_downloaded": len(downloaded),
                "n_skipped_existing": len(skipped),
                "n_failed": len(failed),
                "n_disabled": len(disabled),
                "elapsed_seconds": round(time.time() - start_time, 2),
            },
            "sources": source_results,
        }

        report_path = directories["outputs"] / "downloads_report.json"
        write_report(report, report_path)

        logging.info("Report salvato in: %s", report_path)
        logging.info("Stato finale: %s", final_status)

        return 0 if final_status == "success" else 2

    except Exception as exc:
        fallback_outputs = Path("outputs")
        fallback_outputs.mkdir(parents=True, exist_ok=True)
        fallback_report = fallback_outputs / "downloads_report.json"

        error_report = {
            "timestamp_utc": timestamp_utc,
            "config_used": str(config_path),
            "summary": {
                "final_status": "fatal_error",
                "elapsed_seconds": round(time.time() - start_time, 2),
            },
            "error": str(exc),
        }

        with fallback_report.open("w", encoding="utf-8") as f:
            json.dump(error_report, f, indent=2, ensure_ascii=False)

        print(f"ERRORE FATALE: {exc}", file=sys.stderr)
        print(f"Report errore salvato in: {fallback_report}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
