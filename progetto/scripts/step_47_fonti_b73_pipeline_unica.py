
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AI PRIMETTA V3.3 — STEP 47
PIPELINE UNICA PER LE EVIDENZE BIOLOGICHE B73

--self-test:
    Verifica interna del codice senza connessioni.

--audit:
    Controllo in sola lettura del progetto.

--run:
    Acquisizione dei dati, estrazione, controllo
    scientifico, confronto e report.

--run --offline:
    Rianalisi delle fonti gia' acquisite.

FONTI:
    He et al. 2022
    DOI: 10.1111/tpj.15943
    PMC9804894
    PARE B73 - Dataset S4/S5

    He et al. 2020
    DOI: 10.1534/g3.120.401199
    PMC7202030
    Siccita' B73 e associazioni small RNA/geni

VINCOLI:
    Nessun training ML/GNN.
    Nessuna modifica alla GUI.
    Nessuna modifica a Docker.
    Nessuna modifica ai checkpoint.
    Nessuna etichetta biologica inventata.
"""

import argparse
import collections
import csv
import hashlib
import io
import json
import os
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

from datetime import datetime, timezone
from pathlib import Path


# ============================================================
# 1. CONFIGURAZIONE DELLE FONTI
# ============================================================

SOURCES = {
    "PARE_B73_2022": {
        "pmcid": "PMC9804894",
        "doi": "10.1111/tpj.15943",
        "context": (
            "B73; diversi tessuti e stadi; "
            "evidenze PARE; non siccita'"
        ),
        "expected": "TPJ-112-268-s002.xlsx",
        "article": (
            "https://pmc.ncbi.nlm.nih.gov/"
            "articles/PMC9804894/"
        ),
    },

    "DROUGHT_B73_2020": {
        "pmcid": "PMC7202030",
        "doi": "10.1534/g3.120.401199",
        "context": (
            "B73; siccita' progressiva; "
            "associazioni sRNA-geni"
        ),
        "figshare": 11918463,
        "article": (
            "https://pmc.ncbi.nlm.nih.gov/"
            "articles/PMC7202030/"
        ),
    },
}

CANDIDATES = {
    "sRNA_00000037": "UUCUUGGAUUUAUGAAAGACGAACA",
    "sRNA_00000435": "AGAUGAUACGAGCACCUCUGAGCU",
}

S3 = "https://pmc-oa-opendata.s3.amazonaws.com/"

ALLOWED = (
    "pmc-oa-opendata.s3.amazonaws.com",
    "api.figshare.com",
    "ndownloader.figshare.com",
    "figshare.com",
    "s3.amazonaws.com",
    "s3-eu-west-1.amazonaws.com",
)

ALLOWED = tuple(ALLOWED) + ("www.ebi.ac.uk",)
MAX_FILE = 150 * 1024 * 1024
MAX_UNPACK = 500 * 1024 * 1024
MAX_ROWS = 300000

XNS = (
    "{http://schemas.openxmlformats.org/"
    "spreadsheetml/2006/main}"
)

RNS = (
    "{http://schemas.openxmlformats.org/"
    "officeDocument/2006/relationships}"
)


# ============================================================
# 2. UTILITA'
# ============================================================

def now():
    return datetime.now(
        timezone.utc
    ).strftime("%Y%m%dT%H%M%S%fZ")


def sha(path):

    digest = hashlib.sha256()

    with open(path, "rb") as stream:

        for block in iter(
            lambda: stream.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def put_json(path, obj):

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    path.write_text(
        json.dumps(
            obj,
            indent=2,
            ensure_ascii=False,
            default=str,
        ) + "\n",
        encoding="utf-8",
    )


# ============================================================
# 3. CONTROLLO DEGLI INDIRIZZI
# ============================================================

def url_ok(url):

    parsed = urllib.parse.urlparse(url)

    host = parsed.hostname

    if parsed.scheme != "https" or not host:
        return False

    return (
        host in ALLOWED
        or host.endswith(".figshare.com")
        or host.endswith(".amazonaws.com")
    )


def request_bytes(url, limit=MAX_FILE):

    if not url_ok(url):
        raise RuntimeError(
            "Host HTTPS non autorizzato: "
            + url[:150]
        )

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": (
                "AI-PRIMETTA-V3.3-"
                "scientific-dataset/1.0"
            )
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=55,
    ) as response:

        if not url_ok(response.geturl()):
            raise RuntimeError(
                "Redirect non autorizzato: "
                + response.geturl()[:150]
            )

        data = response.read(limit + 1)

    if len(data) > limit:
        raise RuntimeError(
            "File troppo grande: " + url[:100]
        )

    return data


def get_json(url):

    return json.loads(
        request_bytes(
            url,
            limit=5 * 1024 * 1024,
        )
    )


def is_supported(path):

    return path.suffix.lower() in (
        ".xlsx",
        ".xls",
        ".csv",
        ".tsv",
        ".txt",
        ".zip",
    )


# ============================================================
# 4. VALIDAZIONE DEI FILE ORIGINALI
# ============================================================

def valid_blob(name, blob):

    suffix = Path(name).suffix.lower()

    if suffix in (".xlsx", ".zip"):

        if not zipfile.is_zipfile(
            io.BytesIO(blob)
        ):
            return False

        with zipfile.ZipFile(
            io.BytesIO(blob)
        ) as archive:

            total = sum(
                item.file_size
                for item in archive.infolist()
            )

            if total > MAX_UNPACK:
                return False

            if any(
                item.flag_bits & 1
                for item in archive.infolist()
            ):
                return False

            if archive.testzip() is not None:
                return False

            if suffix == ".xlsx":

                required = {
                    "[Content_Types].xml",
                    "xl/workbook.xml",
                }

                if not required.issubset(
                    archive.namelist()
                ):
                    return False

        return True

    if suffix == ".xls":

        return blob.startswith(
            bytes.fromhex(
                "D0CF11E0A1B11AE1"
            )
        )

    return bool(blob) and (
        b"\x00" not in blob[:1000]
    )


def obtain(url, destination, md5=None):

    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if destination.exists():

        existing = destination.read_bytes()

        if not valid_blob(
            destination.name,
            existing,
        ):
            raise RuntimeError(
                "File esistente non valido: "
                + str(destination)
            )

        return {
            "status": "cached",
            "path": str(destination),
            "sha256": sha(destination),
            "bytes": destination.stat().st_size,
        }

    blob = request_bytes(url)

    if not valid_blob(
        destination.name,
        blob,
    ):
        raise RuntimeError(
            "Formato scaricato non valido: "
            + destination.name
        )

    if md5:

        calculated = hashlib.md5(
            blob
        ).hexdigest().lower()

        if calculated != str(md5).lower():

            raise RuntimeError(
                "Checksum MD5 non corrispondente: "
                + destination.name
            )

    temporary = destination.with_name(
        destination.name + ".download.tmp"
    )

    try:

        with open(temporary, "xb") as stream:
            stream.write(blob)

        os.replace(
            temporary,
            destination,
        )

    finally:

        if temporary.exists():
            temporary.unlink()

    return {
        "status": "downloaded",
        "path": str(destination),
        "url": url,
        "sha256": sha(destination),
        "bytes": len(blob),
    }


# ============================================================
# 5. ACQUISIZIONE METADATI PMC
# ============================================================

def s3_media(pmcid, offline, root):

    cache = (
        root
        / "PMC"
        / (pmcid + "_metadata.json")
    )

    if cache.exists():

        meta = json.loads(
            cache.read_text(
                encoding="utf-8"
            )
        )

    elif offline:

        raise FileNotFoundError(
            "Metadata PMC non disponibili: "
            + str(cache)
        )

    else:

        meta = None
        errors = []

        for version in (1, 2, 3, 4):

            url = (
                S3
                + "metadata/"
                + pmcid
                + f".{version}.json"
            )

            try:

                item = get_json(url)

                if item.get(
                    "pmcid", pmcid
                ) == pmcid:

                    meta = item
                    break

            except Exception as exc:

                errors.append(
                    f"v{version}: "
                    + str(exc)[:80]
                )

        if meta is None:

            raise RuntimeError(
                "Metadata PMC non disponibili: "
                + "; ".join(errors)
            )

        put_json(cache, meta)

    urls = []

    for entry in meta.get(
        "media_urls", []
    ):

        if isinstance(entry, str):

            urls.append(entry)

        elif isinstance(entry, dict):

            for field in (
                "url",
                "uri",
                "href",
            ):

                value = entry.get(field)

                if isinstance(value, str):

                    urls.append(value)
                    break

    return meta, sorted(set(urls))


# ============================================================
# 6. ACQUISIZIONE METADATI FIGSHARE
# ============================================================


# STEP 47 V2 - acquisizione PARE con nomi supplementi flessibili.
def acquire_pare_b73(
    *,
    external,
    offline,
    inventory,
    input_files,
    errors
):
    pmcid = "PMC9804894"

    dest_dir = external / "PMC" / pmcid
    dest_dir.mkdir(parents=True, exist_ok=True)

    failures = []
    available = []

    # 1. Distribuzione ufficiale PMC su S3.
    try:
        _, media = s3_media(
            pmcid,
            offline=offline,
            root=external
        )

        print("[PMC_MEDIA]", len(media))

        for item in media:
            parsed = urllib.parse.urlparse(item)
            name = Path(parsed.path).name

            print("[PMC_ASSET]", name)

            if not name.lower().endswith(".xlsx"):
                continue

            if item.startswith("s3://pmc-oa-opendata/"):
                url = (
                    S3
                    + item.split(
                        "s3://pmc-oa-opendata/", 1
                    )[1]
                )
            elif item.startswith("https://"):
                url = item
            else:
                failures.append(
                    "URI PMC non supportato: " + item[:160]
                )
                continue

            destination = dest_dir / name

            try:
                if offline and not destination.is_file():
                    continue

                if offline:
                    if not valid_blob(
                        name, destination.read_bytes()
                    ):
                        raise RuntimeError(
                            "XLSX locale non valido"
                        )

                    info = {
                        "status": "cached",
                        "sha256": sha(destination),
                        "bytes": destination.stat().st_size
                    }
                else:
                    query = urllib.parse.parse_qs(
                        urllib.parse.urlparse(url).query
                    )

                    md5 = query.get("md5", [None])[0]

                    info = obtain(
                        url,
                        destination,
                        md5=md5
                    )

                info.update({
                    "source": "PARE_B73_2022",
                    "file": name
                })

                inventory.append(info)
                available.append(destination)

                print(
                    "[ACQUISITO] PARE_B73_2022",
                    name,
                    info["status"]
                )

            except Exception as exc:
                failures.append(
                    f"PMC {name}: {exc}"
                )

    except Exception as exc:
        failures.append(
            "PMC metadata: " + str(exc)
        )

    # 2. Fallback: supplementi ufficiali Europe PMC.
    if not available:
        archive = (
            dest_dir
            / "PMC9804894_supplementary_EuropePMC.zip"
        )

        try:
            if offline:
                if not archive.is_file():
                    raise FileNotFoundError(
                        "ZIP Europe PMC non presente in cache"
                    )
            else:
                obtain(
                    "https://www.ebi.ac.uk/europepmc/"
                    "webservices/rest/PMC9804894/"
                    "supplementaryFiles",
                    archive
                )

            with zipfile.ZipFile(archive) as z:
                entries = [
                    entry
                    for entry in z.infolist()
                    if (
                        not entry.is_dir()
                        and Path(
                            entry.filename
                        ).suffix.lower() == ".xlsx"
                    )
                ]

                print(
                    "[EUROPEPMC_XLSX]",
                    len(entries)
                )

                for entry in entries:
                    name = Path(entry.filename).name

                    if entry.file_size > MAX_FILE:
                        failures.append(
                            "XLSX oltre limite: " + name
                        )
                        continue

                    destination = dest_dir / name
                    data = z.read(entry)

                    if not valid_blob(name, data):
                        failures.append(
                            "XLSX non valido: " + name
                        )
                        continue

                    if destination.exists():
                        if sha(destination) != hashlib.sha256(
                            data
                        ).hexdigest():
                            failures.append(
                                "Collisione supplemento: " + name
                            )
                            continue
                    else:
                        destination.write_bytes(data)

                    inventory.append({
                        "status": "cached_from_EuropePMC_ZIP",
                        "source": "PARE_B73_2022",
                        "file": name,
                        "bytes": len(data),
                        "sha256": sha(destination)
                    })

                    available.append(destination)

                    print(
                        "[ACQUISITO] PARE_B73_2022",
                        name,
                        "EuropePMC"
                    )

        except Exception as exc:
            failures.append(
                "Europe PMC: " + str(exc)
            )

    # 3. Registrazione senza duplicazioni.
    for item in dict.fromkeys(available):
        input_files.append(
            ("PARE_B73_2022", item)
        )

    if not available:
        errors.append(
            "PARE B73 non acquisito: "
            + "; ".join(failures)
        )

    elif failures:
        print(
            "[PARE_NOTE] Alcune vie alternative "
            "non sono riuscite:",
            "; ".join(failures)[:500]
        )

def figshare_media(offline, root):

    cache = (
        root
        / "FIGSHARE"
        / "11918463_metadata.json"
    )

    if cache.exists():

        meta = json.loads(
            cache.read_text(
                encoding="utf-8"
            )
        )

    elif offline:

        raise FileNotFoundError(
            "Metadata Figshare non disponibili: "
            + str(cache)
        )

    else:

        meta = get_json(
            "https://api.figshare.com/"
            "v2/articles/11918463"
        )

        put_json(cache, meta)

    files = [
        item
        for item in meta.get("files", [])
        if isinstance(item, dict)
    ]

    return meta, files


# ============================================================
# 7. LETTURA XLSX SENZA DIPENDENZE ESTERNE
# ============================================================

def extract_xlsx(path):

    with zipfile.ZipFile(path) as archive:

        strings = []

        if (
            "xl/sharedStrings.xml"
            in archive.namelist()
        ):

            with archive.open(
                "xl/sharedStrings.xml"
            ) as stream:

                for _, element in ET.iterparse(
                    stream,
                    events=("end",),
                ):

                    if element.tag == XNS + "si":

                        value = "".join(
                            text.text or ""
                            for text in element.iter(
                                XNS + "t"
                            )
                        )

                        strings.append(value)
                        element.clear()

        workbook = ET.fromstring(
            archive.read(
                "xl/workbook.xml"
            )
        )

        relationships = ET.fromstring(
            archive.read(
                "xl/_rels/workbook.xml.rels"
            )
        )

        relations = {
            element.get("Id"):
                element.get("Target")
            for element in relationships
        }

        paths = set(
            archive.namelist()
        )

        for sheet in workbook.iter(
            XNS + "sheet"
        ):

            name = sheet.get(
                "name", "Sheet"
            )

            target = relations.get(
                sheet.get(RNS + "id")
            )

            if not target:
                continue

            if target.startswith("/"):

                entry = target.lstrip("/")

            else:

                import posixpath

                entry = posixpath.normpath(
                    "xl/" + target
                )

            if entry not in paths:
                continue

            with archive.open(entry) as stream:

                for _, element in ET.iterparse(
                    stream,
                    events=("end",),
                ):

                    if element.tag != XNS + "row":
                        continue

                    row = {}

                    for cell in element.findall(
                        XNS + "c"
                    ):

                        coordinate = cell.get(
                            "r", ""
                        )

                        column = re.match(
                            "[A-Z]+",
                            coordinate,
                        )

                        if not column:
                            continue

                        kind = cell.get("t")

                        if kind == "inlineStr":

                            value = "".join(
                                text.text or ""
                                for text in cell.iter(
                                    XNS + "t"
                                )
                            )

                        else:

                            raw = cell.find(
                                XNS + "v"
                            )

                            value = (
                                raw.text
                                if (
                                    raw is not None
                                    and raw.text
                                )
                                else ""
                            )

                            if kind == "s" and value:

                                value = strings[
                                    int(value)
                                ]

                        row[
                            column.group()
                        ] = str(value).strip()

                    if row:

                        yield (
                            name,
                            int(
                                element.get("r")
                                or 0
                            ),
                            row,
                        )

                    element.clear()


# ============================================================
# 8. LETTURA TABELLE CSV / TSV / XLS
# ============================================================

def xl_col(number):

    result = ""

    while True:

        result = (
            chr(number % 26 + 65)
            + result
        )

        number = number // 26 - 1

        if number < 0:
            return result


def spreadsheet_tables(path):

    extension = path.suffix.lower()

    if extension == ".xlsx":

        yield from extract_xlsx(path)

    elif extension == ".xls":

        try:
            import xlrd

        except ImportError as exc:

            raise RuntimeError(
                "Supplemento XLS legacy. "
                "Installare xlrd nella .venv "
                "e rieseguire lo stesso script."
            ) from exc

        workbook = xlrd.open_workbook(
            str(path),
            on_demand=True,
        )

        for sheet in workbook.sheets():

            for index in range(
                sheet.nrows
            ):

                row = {}

                for column in range(
                    sheet.ncols
                ):

                    value = str(
                        sheet.cell_value(
                            index,
                            column,
                        )
                    ).strip()

                    if value:

                        row[
                            xl_col(column)
                        ] = value

                if row:

                    yield (
                        sheet.name,
                        index + 1,
                        row,
                    )

    elif extension in (
        ".csv",
        ".tsv",
        ".txt",
    ):

        with path.open(
            encoding="utf-8-sig",
            errors="replace",
            newline="",
        ) as stream:

            sample = stream.read(8192)
            stream.seek(0)

            delimiter = (
                "\t"
                if extension == ".tsv"
                else ","
            )

            if extension != ".tsv":

                try:

                    delimiter = (
                        csv.Sniffer()
                        .sniff(
                            sample,
                            delimiters=",\t;",
                        )
                        .delimiter
                    )

                except csv.Error:
                    pass

            reader = csv.reader(
                stream,
                delimiter=delimiter,
            )

            for number, values in enumerate(
                reader,
                1,
            ):

                row = {
                    xl_col(index):
                        str(value).strip()
                    for index, value
                    in enumerate(values)
                    if str(value).strip()
                }

                if row:

                    yield (
                        path.stem,
                        number,
                        row,
                    )

    else:

        raise RuntimeError(
            "Formato non supportato: "
            + extension
        )


# ============================================================
# 9. NORMALIZZAZIONE IDENTIFICATIVI
# ============================================================

def tidy(value):

    return re.sub(
        "[^a-z0-9]+",
        "",
        str(value).lower(),
    )


def seq(value):

    normalized = re.sub(
        r"\s+",
        "",
        str(value)
        .upper()
        .replace("U", "T"),
    )

    if re.fullmatch(
        "[ACGT]{18,35}",
        normalized,
    ):

        return normalized

    return ""


def norm_gene(value):

    text = str(value).strip()

    if re.fullmatch(
        r"zma-MIR[0-9]+[A-Za-z0-9-]*",
        text,
        re.I
    ):
        return text

    if re.fullmatch(
        r"[A-Z]{2}[0-9]{6}\.[0-9]+_FG[0-9]+",
        text,
        re.I
    ):
        return text

    found = re.search(
        r"(Zm00001(?:eb|d)[0-9]+|GRMZM[25]G[0-9]+)",
        text,
        re.I
    )

    return found.group(1) if found else ""



def gene_build(identifier):

    text = str(identifier)

    if text.startswith("Zm00001eb"):
        return "NAM5"

    if text.startswith("Zm00001d"):
        return "RefGen_v4"

    if text.startswith(("GRMZM2G", "GRMZM5G")):
        return "older_build"

    if re.fullmatch(
        r"zma-MIR[0-9]+[A-Za-z0-9-]*",
        text,
        re.I
    ):
        return "miRNA_precursor_locus_not_mRNA"

    if re.fullmatch(
        r"[A-Z]{2}[0-9]{6}\.[0-9]+_FG[0-9]+",
        text,
        re.I
    ):
        return "legacy_FG_accession_unmapped"

    return "unresolved"



# ============================================================
# 10. RICONOSCIMENTO INTESTAZIONI
# ============================================================

def identify_header(rows):

    patterns = {
        "mirna":
            r"(?:^|[^a-z])"
            r"(mirna|micro.?rna|small.?rna|srna)"
            r"(?:$|[^a-z])",

        "gene":
            r"(?:target\s*"
            r"(?:gene|transcript|id|mrna)?"
            r"|gene\s*(?:id|name)"
            r"|transcript\s*(?:id|name)"
            r"|\bdrg\b)",

        "category":
            r"\b(?:pare|cleavage|degradome)?"
            r"\s*categor(?:y|ies)\b",

        "sequence":
            r"\b(?:mirna|srna|small.?rna|mature)?"
            r"\s*sequence\b",

        "tissue":
            r"\b(?:tissue|stage|sample|organ|library)\b",

        "correlation":
            r"\b(?:correlation|pearson|spearman"
            r"|r.?value|anti.?correlation)\b",

        "cleavage":
            r"\b(?:cleavage\s*site|p.?value"
            r"|p.?adjust|t.?plot|pare\s*reads)\b",

        "efficiency":
            r"\b(?:local|global)?"
            r"\s*cleavage\s*efficien(?:cy|cies)\b",
    }

    best = None

    for index, (_, row) in enumerate(
        rows[:30]
    ):

        found = {}

        for column, value in row.items():

            for key, pattern in patterns.items():

                if (
                    key not in found
                    and re.search(
                        pattern,
                        value,
                        re.I,
                    )
                ):

                    found[key] = column

        score = (
            4 * int("mirna" in found)
            + 4 * int("gene" in found)
            + 2 * int("category" in found)
            + len(found)
        )

        if best is None or score > best[0]:

            best = (
                score,
                index,
                found,
            )

    if (
        best
        and best[0] >= 9
        and "mirna" in best[2]
        and "gene" in best[2]
    ):

        return best[1], best[2]

    return None, {}


def section_from_name(name, rows):

    text = (
        name
        + " "
        + " ".join(
            " ".join(row.values())
            for _, row in rows[:12]
        )
    )

    match = re.search(
        r"\b(?:dataset|table)\s*S\s*([1-9])\b",
        text,
        flags=re.I,
    )

    if not match:

        match = re.search(
            r"^\s*S\s*([1-9])\b",
            name,
            flags=re.I,
        )

    return (
        "S" + match.group(1)
        if match
        else None
    )


# ============================================================
# 11. ESTRAZIONE DELLE EVIDENZE
# ============================================================


def parse_collection(files, out_dir):
    # STEP47_SCHEMI_B73_V5
    # PARE B73 S4/S5 e siccità B73 Table S7.

    profiles = []
    matches = []
    evidence = []
    failures = []

    tables_dir = out_dir / "TABLES_EXTRACTED"
    tables_dir.mkdir(parents=True, exist_ok=True)

    candidates = {
        seq(value): ident
        for ident, value in CANDIDATES.items()
    }

    s7_names = {
        "srna-gene pairs",
        "srna-drg pairs",
        "mirna-drg pairs",
    }

    for source, path in files:
        try:
            grouped = collections.defaultdict(list)
            row_count = 0

            for sheet, number, row in spreadsheet_tables(path):
                grouped[sheet].append((number, row))
                row_count += 1

                if row_count > MAX_ROWS:
                    raise RuntimeError(
                        "Limite righe superato: " + path.name
                    )

            for sheet, rows in grouped.items():

                section = section_from_name(sheet, rows)

                is_pare = (
                    source == "PARE_B73_2022"
                    and section in {"S4", "S5"}
                )

                is_s7 = (
                    source == "DROUGHT_B73_2020"
                    and path.name.casefold() == "table s7.xlsx"
                    and sheet.casefold().strip() in s7_names
                )

                header_index = None
                detected = {}

                # Identificazione delle intestazioni reali.
                if is_pare or is_s7:

                    for i, (_, row) in enumerate(rows[:20]):

                        if is_pare and all(
                            tidy(row.get(col, "")) == name
                            for col, name in {
                                "A": "tissue",
                                "B": "siteid",
                                "C": "query",
                                "D": "transcript",
                                "E": "degradomecategory",
                                "F": "degradomepval",
                            }.items()
                        ):
                            header_index = i
                            detected = {
                                "mirna": "C",
                                "gene": "D",
                                "tissue": "A",
                                "siteid": "B",
                                "category": "E",
                                "pvalue": "F",
                            }
                            break

                        if is_s7 and (
                            tidy(row.get("A", "")) == "geneid"
                            and tidy(row.get("B", ""))
                            in {"srna", "mirna"}
                            and "correlationbetween"
                            in tidy(row.get("G", ""))
                        ):
                            header_index = i
                            detected = {
                                "gene": "A",
                                "sequence": "B",
                                "correlation": "G",
                                "interaction_padj": "T",
                            }
                            break

                else:
                    header_index, detected = identify_header(rows)

                # Conservazione delle tabelle originali.
                columns = sorted(
                    {
                        col
                        for _, row in rows
                        for col in row
                    },
                    key=lambda col: (len(col), col),
                )

                filename = (
                    f"{source}_{path.stem}_"
                    f"{tidy(sheet)[:24]}_{len(profiles)}.tsv"
                )

                output_file = tables_dir / filename

                with output_file.open(
                    "w",
                    encoding="utf-8",
                    newline="",
                ) as stream:

                    writer = csv.DictWriter(
                        stream,
                        fieldnames=["excel_row"] + columns,
                        delimiter="\t",
                    )

                    writer.writeheader()

                    for number, row in rows:
                        writer.writerow({
                            "excel_row": number,
                            **row,
                        })

                profiles.append({
                    "source": source,
                    "file": path.name,
                    "sheet": sheet,
                    "rows": len(rows),
                    "section": section,
                    "header_excel_row": (
                        rows[header_index][0]
                        if header_index is not None
                        else None
                    ),
                    "recognized_columns": detected,
                    "preview_rows": [
                        row for _, row in rows[:2]
                    ],
                    "extracted_tsv": str(
                        output_file.relative_to(out_dir)
                    ),
                })

                if not (is_pare or is_s7):
                    continue

                if header_index is None:
                    raise RuntimeError(
                        "Schema non riconosciuto: "
                        + path.name + " / " + sheet
                    )

                # Elaborazione delle sole righe dati.
                for number, row in rows[header_index + 1:]:

                    if is_pare:

                        mirna_id = row.get("C", "").strip()
                        gene_raw = row.get("D", "").strip()

                        # La colonna I di S4 contiene un duplex.
                        # Non equivale alla sequenza matura.
                        rna_seq = ""

                        if not mirna_id or not gene_raw:
                            continue

                        category = row.get("E", "").strip()
                        sample = row.get("A", "").strip()

                        correlation = None
                        direction = "not_applicable"

                        evidence_type = (
                            "PARE_cleavage_reported_in_B73_non_drought"
                            if section == "S4"
                            else "PARE_cleavage_efficiency_S5"
                        )

                        level = (
                            "PARE_category_" + category
                            if category in {"0", "1", "2", "3", "4"}
                            else "PARE_category_unresolved"
                        )

                    else:

                        # S7 colonna B: sequenza nucleotidica.
                        mirna_id = ""
                        gene_raw = row.get("A", "").strip()
                        rna_seq = seq(row.get("B", ""))

                        if not gene_raw or not rna_seq:
                            continue

                        category = ""
                        sample = ""

                        try:
                            correlation = float(
                                row.get("G", "")
                            )

                            if not -1 <= correlation <= 1:
                                correlation = None

                        except (ValueError, TypeError):
                            correlation = None

                        if correlation is None:
                            direction = "unresolved"
                        elif correlation < 0:
                            direction = "negative"
                        else:
                            direction = "non_negative"

                        evidence_type = (
                            "association_only_not_PARE_cleavage"
                        )

                        level = (
                            "drought_"
                            + tidy(sheet)
                            + "_association"
                        )

                    gene_id = norm_gene(gene_raw)

                    if not gene_id:
                        continue

                    candidate = candidates.get(
                        rna_seq, ""
                    )

                    if candidate:
                        matches.append({
                            "source": source,
                            "file": path.name,
                            "sheet": sheet,
                            "row": number,
                            "cell": "B",
                            "candidate": candidate,
                            "exact_sequence": rna_seq,
                            "status":
                                "exact_sequence_not_"
                                "functional_validation",
                        })

                    evidence.append({
                        "source": source,
                        "sheet": sheet,
                        "row": number,

                        "mirna_id": mirna_id,
                        "mirna_sequence": rna_seq,

                        "target_raw": gene_raw,
                        "target_gene_id": gene_id,
                        "target_gene_build":
                            gene_build(gene_id),

                        "pare_category": category,
                        "tissue_or_stage": sample,

                        "pare_site_id": (
                            row.get("B", "")
                            if is_pare else ""
                        ),

                        "pare_pvalue_raw": (
                            row.get("F", "")
                            if is_pare else ""
                        ),

                        "pare_allen_score": (
                            row.get("G", "")
                            if is_pare and section == "S4"
                            else ""
                        ),

                        "s5_ValueP1": (
                            row.get("G", "")
                            if is_pare and section == "S5"
                            else ""
                        ),

                        "s5_MedianFlank2": (
                            row.get("H", "")
                            if is_pare and section == "S5"
                            else ""
                        ),

                        "s5_TotalValue3": (
                            row.get("I", "")
                            if is_pare and section == "S5"
                            else ""
                        ),

                        "psRNATarget_score_raw": (
                            row.get("F", "")
                            if is_s7 else ""
                        ),

                        "interaction_padj_raw": (
                            row.get("T", "")
                            if is_s7 else ""
                        ),

                        "target_position_raw": (
                            row.get("E", "")
                            if is_s7 else ""
                        ),

                        "drought_regulation_raw": (
                            row.get("X", "")
                            if is_s7 else ""
                        ),

                        "correlation_r": correlation,
                        "correlation_direction": direction,

                        "evidence_type": evidence_type,
                        "evidence_level": level,

                        "candidate_exact_match": candidate,

                        "not_validated_in_PRJNA419326": "1",
                    })

        except Exception as exc:

            failures.append({
                "source": source,
                "file": str(path),
                "error": str(exc),
            })

    return profiles, matches, evidence, failures


def write_tsv(path, rows, fields):

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as stream:

        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# 13. CONTROLLO GENOMICO REFGEN / NAM5
# ============================================================

def check_crosswalk(evidence, crosswalk):

    mappings = {}

    if crosswalk:

        with crosswalk.open(
            encoding="utf-8-sig",
            newline="",
        ) as stream:

            reader = csv.DictReader(
                stream,
                delimiter="\t",
            )

            required = {
                "source_gene_id",
                "nam5_gene_id",
            }

            if not required.issubset(
                reader.fieldnames or []
            ):

                raise RuntimeError(
                    "Crosswalk: richieste colonne "
                    "source_gene_id e nam5_gene_id"
                )

            for item in reader:

                source_id = item[
                    "source_gene_id"
                ].strip()

                target_id = item[
                    "nam5_gene_id"
                ].strip()

                if (
                    not source_id
                    or not re.fullmatch(
                        r"Zm00001eb\d+",
                        target_id,
                    )
                ):

                    raise RuntimeError(
                        "Crosswalk con valore non valido"
                    )

                if (
                    source_id in mappings
                    and mappings[source_id]
                    != target_id
                ):

                    raise RuntimeError(
                        "Crosswalk non univoca per "
                        + source_id
                    )

                mappings[
                    source_id
                ] = target_id

    for item in evidence:

        gene = item[
            "target_gene_id"
        ]

        if gene.startswith("Zm00001eb"):

            item[
                "nam5_gene_id"
            ] = gene

            item[
                "nam5_mapping_status"
            ] = (
                "native_NAM5_id_"
                "unverified_presence"
            )

        elif gene in mappings:

            item[
                "nam5_gene_id"
            ] = mappings[gene]

            item[
                "nam5_mapping_status"
            ] = (
                "provided_external_crosswalk_"
                "requires_provenance"
            )

        else:

            item[
                "nam5_gene_id"
            ] = ""

            item[
                "nam5_mapping_status"
            ] = (
                "unmapped_no_assembly_"
                "equivalence_assumed"
            )


# ============================================================
# 14. SELF-TEST
# ============================================================

def self_test():

    import tempfile

    assert seq(
        "AUGCAUGCAUGCAUGCAUGC"
    ) == "ATGCATGCATGCATGCATGC"

    assert seq("AUGC-no") == ""

    assert xl_col(0) == "A"
    assert xl_col(26) == "AA"

    assert (
        gene_build("Zm00001d123456")
        == "RefGen_v4"
    )

    assert (
        gene_build("Zm00001eb123456")
        == "NAM5"
    )

    header = [
        (
            1,
            {
                "A": "miRNA",
                "B": "Target gene",
                "C": "Category",
            },
        ),
        (
            2,
            {
                "A": "zma-miR156",
                "B": "Zm00001d123456",
                "C": "0",
            },
        ),
    ]

    index, detected = identify_header(
        header
    )

    assert index == 0

    assert {
        "mirna",
        "gene",
        "category",
    }.issubset(detected)

    with tempfile.TemporaryDirectory() as temp:

        path = Path(temp) / "test.csv"

        path.write_text(
            "miRNA,target gene,Category\n"
            "zma-miR156,Zm00001d123456,0\n",
            encoding="utf-8",
        )

        assert len(
            list(
                spreadsheet_tables(path)
            )
        ) == 2

    print(
        "[SELF_TEST] PASS: "
        "schemi, sequenze, ID genomici, CSV"
    )


# ============================================================
# 15. PROGRAMMA PRINCIPALE
# ============================================================


# ============================================================
# STEP47_QC_V6_UNICO
# QC delle righe escluse, deduplicazione e assembly mapping.
# Nessun training o modifica dei dati scientifici originali.
# ============================================================

def step47_qc_v6(out):

    report_path = out / "RAPPORTO_STEP47.json"

    evidence_path = (
        out / "ASSOCIAZIONI_EVIDENZE_DA_VERIFICARE.tsv"
    )

    report = json.loads(
        report_path.read_text(encoding="utf-8")
    )

    def read_tsv(path):

        with path.open(
            encoding="utf-8-sig",
            newline=""
        ) as stream:

            return list(
                csv.DictReader(
                    stream,
                    delimiter="\t"
                )
            )

    evidence = read_tsv(evidence_path)

    evidence_ids = {
        (
            r["source"],
            r["sheet"],
            r["row"]
        )
        for r in evidence
    }

    # ========================================================
    # 1. IDENTIFICAZIONE DELLE TABELLE PERTINENTI
    # ========================================================

    relevant = [
        profile
        for profile in report["sheets"]
        if (
            profile["source"] == "PARE_B73_2022"
            and profile["section"] in ("S4", "S5")
        )
        or (
            profile["source"] == "DROUGHT_B73_2020"
            and profile["file"].casefold() == "table s7.xlsx"
            and profile["sheet"].casefold().strip() in {
                "srna-gene pairs",
                "srna-drg pairs",
                "mirna-drg pairs"
            }
        )
    ]

    if not any(
        x["source"] == "PARE_B73_2022"
        and x["section"] == "S4"
        for x in relevant
    ):
        raise RuntimeError(
            "QC V6: Dataset PARE S4 mancante"
        )

    if not any(
        x["source"] == "PARE_B73_2022"
        and x["section"] == "S5"
        for x in relevant
    ):
        raise RuntimeError(
            "QC V6: Dataset PARE S5 mancante"
        )

    if not any(
        x["source"] == "DROUGHT_B73_2020"
        for x in relevant
    ):
        raise RuntimeError(
            "QC V6: Table S7 mancante"
        )

    # ========================================================
    # 2. ANALISI DELLE RIGHE ESCLUSE
    # ========================================================

    excluded = []
    raws = {}
    per_sheet = []

    for profile in relevant:

        source = profile["source"]
        sheet = profile["sheet"]
        filename = profile["file"]
        header = profile.get("header_excel_row")

        if header is None:
            raise RuntimeError(
                "Intestazione non riconosciuta: "
                + filename + " / " + sheet
            )

        raw_rows = read_tsv(
            out / profile["extracted_tsv"]
        )

        n_scanned = 0
        n_kept = 0
        n_excluded = 0

        for row in raw_rows:

            try:
                line = int(
                    row.get("excel_row", "0")
                )
            except ValueError:
                raise RuntimeError(
                    "Numero riga non valido: " + filename
                )

            if line <= int(header):
                continue

            n_scanned += 1

            key = (
                source,
                sheet,
                str(line)
            )

            if key in raws:
                raise RuntimeError(
                    "Collisione chiave riga: " + str(key)
                )

            raws[key] = row

            if key in evidence_ids:
                n_kept += 1
                continue

            n_excluded += 1

            if source == "PARE_B73_2022":

                if (
                    not row.get("C", "").strip()
                    or not row.get("D", "").strip()
                ):
                    reason = "query_o_trascritto_vuoto"

                elif not norm_gene(row["D"]):
                    reason = (
                        "identificativo_trascritto_"
                        "non_riconosciuto"
                    )

                else:
                    reason = (
                        "altra_esclusione_da_rivedere"
                    )

            else:

                if (
                    not row.get("A", "").strip()
                    or not row.get("B", "").strip()
                ):
                    reason = "gene_o_sequenza_vuota"

                elif not seq(row["B"]):
                    reason = (
                        "sequenza_non_canonica_o_"
                        "lunghezza_fuori_18_35"
                    )

                elif not norm_gene(row["A"]):
                    reason = (
                        "identificativo_gene_"
                        "non_riconosciuto"
                    )

                else:
                    reason = (
                        "altra_esclusione_da_rivedere"
                    )

            excluded.append({
                "source": source,
                "file": filename,
                "sheet": sheet,
                "excel_row": line,
                "reason": reason,
                "col_A": row.get("A", "")[:160],
                "col_B": row.get("B", "")[:160],
                "col_C": row.get("C", "")[:160],
                "col_D": row.get("D", "")[:160],
            })

        per_sheet.append({
            "source": source,
            "sheet": sheet,
            "rows_after_header": n_scanned,
            "included": n_kept,
            "excluded": n_excluded,
        })

        assert (
            n_scanned == n_kept + n_excluded
        ), "Bilancio delle righe non coerente"

    # ========================================================
    # 3. DEDUPLICAZIONE PARE S4 E S5
    # ========================================================

    s4 = {}
    s5 = {}

    for r in evidence:

        typ = r.get("evidence_type", "")

        if typ not in (
            "PARE_cleavage_reported_in_B73_non_drought",
            "PARE_cleavage_efficiency_S5"
        ):
            continue

        # La chiave contiene l'identità dell'evento.
        # Non converte gli identificativi genomici.

        key = (
            tidy(r.get('tissue_or_stage', '')),
            tidy(r.get('mirna_id', '')).removeprefix('zma'),
            tidy(r.get('target_raw', '')),
            tidy(r.get('pare_site_id', ''))
        )

        bucket = (
            s4
            if typ == (
                "PARE_cleavage_reported_in_B73_non_drought"
            )
            else s5
        )

        bucket.setdefault(key, []).append(r)

    pare_events = []
    pval_invalid = 0

    for key, items in sorted(s4.items()):

        categories = sorted({
            x.get("pare_category", "")
            for x in items
        })

        pvals = []

        for item in items:

            try:
                value = float(
                    item.get("pare_pvalue_raw", "")
                )

                if not 0 <= value <= 1:
                    raise ValueError()

                pvals.append(value)

            except (ValueError, TypeError):
                pval_invalid += 1

        pare_events.append({
            "tissue": key[0],
            "mirna_id": key[1],
            "transcript": key[2],
            "site_id": key[3],
            "n_rows_s4": len(items),
            "categories_s4": ";".join(categories),
            "min_reported_pval_s4": (
                min(pvals) if pvals else ""
            ),
            "has_s5_same_normalized_key": int(key in s5),
            "interpretation": (
                "PARE_reported_non_drought_"
                "not_training_label"
            ),
        })

    # ========================================================
    # 4. DEDUPLICAZIONE DELLE ASSOCIAZIONI S7
    # ========================================================

    s7_groups = {}

    invalid_corr = 0
    invalid_padj = 0

    for r in evidence:

        if r["source"] != "DROUGHT_B73_2020":
            continue

        row = raws.get((
            r["source"],
            r["sheet"],
            r["row"]
        ))

        if row is None:
            raise RuntimeError(
                "Riga S7 assente dai dati estratti"
            )

        sequence = seq(row.get("B", ""))
        gene = row.get("A", "").strip()
        transcript = row.get("D", "").strip()
        site = row.get("E", "").strip()

        key = (
            sequence,
            gene,
            transcript,
            site
        )

        s7_groups.setdefault(key, []).append(r)

        try:
            correlation = float(
                row.get("G", "")
            )

            if not -1 <= correlation <= 1:
                raise ValueError()

        except (ValueError, TypeError):
            invalid_corr += 1

        try:
            padj = float(
                row.get("T", "")
            )

            if not 0 <= padj <= 1:
                raise ValueError()

        except (ValueError, TypeError):
            invalid_padj += 1

    s7_pairs = []

    for key, items in sorted(
        s7_groups.items()
    ):

        s7_pairs.append({
            "srna_sequence_DNA": key[0],
            "target_gene_original": key[1],
            "target_transcript_original": key[2],
            "target_position": key[3],
            "n_occurrences": len(items),
            "source_sheets": ";".join(
                sorted({
                    x["sheet"]
                    for x in items
                })
            ),
            "assembly_mapping": "NOT_VERIFIED",
            "evidence": (
                "negative_correlation_and_in_silico_"
                "targeting_not_cleavage_validation"
            ),
        })

    # ========================================================
    # 5. CONTROLLO ASSEMBLY GENOMICO
    # ========================================================

    mapping = collections.Counter(
        (
            r.get(
                "target_gene_build",
                "unresolved"
            ),
            r.get(
                "nam5_mapping_status",
                "missing_status"
            )
        )
        for r in evidence
    )

    mapping_rows = [
        {
            "source_build": build,
            "mapping_status": status,
            "n_rows": count,
        }
        for (build, status), count
        in sorted(mapping.items())
    ]

    # ========================================================
    # 6. ESPORTAZIONE DEI RISULTATI QC
    # ========================================================

    write_tsv(
        out / "QC_RIGHE_ESCLUSE_V6.tsv",
        excluded,
        [
            "source", "file", "sheet",
            "excel_row", "reason",
            "col_A", "col_B", "col_C", "col_D"
        ],
    )

    write_tsv(
        out / "PARE_S4_SITI_UNIVOCI_V6.tsv",
        pare_events,
        [
            "tissue", "mirna_id", "transcript",
            "site_id", "n_rows_s4",
            "categories_s4",
            "min_reported_pval_s4",
            "has_s5_same_normalized_key",
            "interpretation"
        ],
    )

    write_tsv(
        out / "S7_ASSOCIAZIONI_UNIVOCHE_V6.tsv",
        s7_pairs,
        [
            "srna_sequence_DNA",
            "target_gene_original",
            "target_transcript_original",
            "target_position",
            "n_occurrences",
            "source_sheets",
            "assembly_mapping",
            "evidence"
        ],
    )

    write_tsv(
        out / "QC_ASSEMBLY_MAPPING_V6.tsv",
        mapping_rows,
        [
            "source_build",
            "mapping_status",
            "n_rows"
        ],
    )

    # ========================================================
    # 7. RAPPORTO SCIENTIFICO
    # ========================================================

    stats = {
        "method": (
            "QC sulle righe originali; nessuna "
            "equivalenza genomica presunta"
        ),

        "per_sheet": per_sheet,

        "excluded_total": len(excluded),

        "excluded_reasons": dict(
            collections.Counter(
                x["reason"]
                for x in excluded
            )
        ),

        "pare_s4_rows": sum(
            len(x)
            for x in s4.values()
        ),

        "pare_s4_unique_site_keys": len(s4),

        "pare_s4_duplicate_occurrences": sum(
            len(x) - 1
            for x in s4.values()
        ),

        "pare_s4_pvalues_missing_or_invalid":
            pval_invalid,

        "pare_s5_rows": sum(
            len(x)
            for x in s5.values()
        ),

        "pare_s5_unique_site_keys": len(s5),

        "pare_s5_rows_with_normalized_s4_key": sum(
            len(x)
            for key, x in s5.items()
            if key in s4
        ),

        "s7_rows": sum(
            len(x)
            for x in s7_groups.values()
        ),

        "s7_unique_sequence_gene_transcript_site_keys":
            len(s7_groups),

        "s7_repeated_rows_over_three_sheets": sum(
            len(x) - 1
            for x in s7_groups.values()
        ),

        "s7_correlation_invalid_or_missing":
            invalid_corr,

        "s7_interaction_padj_invalid_or_missing":
            invalid_padj,

        "assembly_mapping_rows": mapping_rows,

        "candidate_matching_limit": (
            "Sequenze complete disponibili in S7; "
            "S4 richiede associazione degli ID miRNA "
            "a sequenze mature verificate"
        ),

        "cross_assembly_overlap_claimed": False,

        "suitable_as_validated_ml_training_labels": False,

        "original_files_modified": False,
    }

    put_json(
        out / "QC_STEP47_V6.json",
        stats
    )

    report["quality_control_v6"] = stats

    put_json(
        report_path,
        report
    )

    with (
        out / "RAPPORTO_STEP47.md"
    ).open(
        "a",
        encoding="utf-8"
    ) as stream:

        stream.write(
            "\n## Controllo scientifico V6\n\n"
        )

        stream.write(
            f"- Righe escluse: {len(excluded)}\n"
        )

        for item in per_sheet:

            stream.write(
                f"- {item['source']} / "
                f"{item['sheet']}: "
                f"{item['included']} incluse, "
                f"{item['excluded']} escluse\n"
            )

        stream.write(
            f"- S4 siti univoci: {len(s4)}\n"
        )

        stream.write(
            "- S7 associazioni univoche: "
            f"{len(s7_groups)}\n"
        )

        stream.write(
            "- Mapping NAM-5.0 non dedotto "
            "senza evidenza documentata.\n"
        )

        stream.write(
            "- Dati non convertiti automaticamente "
            "in etichette validate per GNN.\n"
        )

    # ========================================================
    # 8. OUTPUT DI CONTROLLO
    # ========================================================

    print(
        "[QC_V6_RIGHE_ESCLUSE]",
        len(excluded)
    )

    print(
        "[QC_V6_MOTIVI]",
        stats["excluded_reasons"]
    )

    print(
        "[QC_V6_S4_SITI_UNIVOCI]",
        len(s4)
    )

    print(
        "[QC_V6_S5_CON_CHIAVE_S4]",
        stats["pare_s5_rows_with_normalized_s4_key"]
    )

    print(
        "[QC_V6_S7_ASSOCIAZIONI_UNIVOCHE]",
        len(s7_groups)
    )

    print(
        "[QC_V6_ASSEMBLY]",
        mapping_rows
    )

    print(
        "[QC_V6_REPORT]",
        out / "QC_STEP47_V6.json"
    )




# STEP47_QC_V7_ENTITY_FIXED
# Classificazione descrittiva, non etichette ML validate.

def step47_qc_v7(out):

    source = (
        out / "ASSOCIAZIONI_EVIDENZE_DA_VERIFICARE.tsv"
    )

    with source.open(
        encoding="utf-8-sig",
        newline=""
    ) as stream:

        rows = list(
            csv.DictReader(
                stream,
                delimiter="\t"
            )
        )

    if not rows:
        raise RuntimeError(
            "V7: tabella evidenze vuota"
        )

    category = {
        "miRNA_precursor_locus_not_mRNA":
            "miRNA_precursor_locus",

        "legacy_FG_accession_unmapped":
            "historical_FG_gene_identifier",

        "RefGen_v4":
            "RefGen_v4_gene_model_biotype_unverified",

        "NAM5":
            "NAM5_gene_model_biotype_unverified",

        "older_build":
            "legacy_gene_model_biotype_unverified",

        "unresolved":
            "unresolved",
    }

    curated = []

    for r in rows:

        build = r.get(
            "target_gene_build",
            "unresolved"
        )

        curated.append({
            "source": r.get("source", ""),
            "sheet": r.get("sheet", ""),
            "excel_row": r.get("row", ""),
            "target_original": r.get("target_raw", ""),
            "target_id": r.get("target_gene_id", ""),
            "assembly_class": build,
            "entity_class": category.get(
                build,
                "unresolved"
            ),
            "nam5_mapping_status": r.get(
                "nam5_mapping_status",
                ""
            ),
            "validated_ml_label": "NO",
        })

    write_tsv(
        out / "QC_ENTITA_TARGET_V7.tsv",
        curated,
        list(curated[0])
    )

    totals = collections.Counter(
        (r["source"], r["entity_class"])
        for r in curated
    )

    counts = [
        {
            "source": src,
            "entity_class": cls,
            "rows": count
        }
        for (src, cls), count
        in sorted(totals.items())
    ]

    write_tsv(
        out / "QC_CLASSI_TARGET_V7.tsv",
        counts,
        ["source", "entity_class", "rows"]
    )

    qc_path = out / "QC_STEP47_V6.json"

    qc_v6 = json.loads(
        qc_path.read_text(encoding="utf-8")
    )

    info = {
        "classification": counts,

        "previously_excluded_after_correction":
            qc_v6.get("excluded_total"),

        "pare_mir_precursor_rows": sum(
            r["entity_class"] == "miRNA_precursor_locus"
            for r in curated
        ),

        "s7_legacy_fg_rows": sum(
            r["source"] == "DROUGHT_B73_2020"
            and r["entity_class"]
                == "historical_FG_gene_identifier"
            for r in curated
        ),

        "s5_with_s4_normalized_key":
            qc_v6.get(
                "pare_s5_rows_with_normalized_s4_key"
            ),

        "nam5_crosswalk_verified": False,
        "validated_ml_training_labels_created": False,

        "note": (
            "Identificativi recuperati: biotipo e assembly "
            "non dedotti dai soli prefissi."
        ),
    }

    put_json(
        out / "QC_STEP47_V7.json",
        info
    )

    main_report = (
        out / "RAPPORTO_STEP47.json"
    )

    report = json.loads(
        main_report.read_text(encoding="utf-8")
    )

    report["target_entity_curation_v7"] = info

    put_json(main_report, report)

    print(
        "[V7_RIGHE_ESCLUSE]",
        info["previously_excluded_after_correction"]
    )

    print(
        "[V7_PRECURSORI_MIR_PARE]",
        info["pare_mir_precursor_rows"]
    )

    print(
        "[V7_ACCESSIONS_FG_S7]",
        info["s7_legacy_fg_rows"]
    )

    print(
        "[V7_S5_MATCH_NORMALIZZATO]",
        info["s5_with_s4_normalized_key"]
    )

    print(
        "[V7_REPORT]",
        out / "QC_STEP47_V7.json"
    )

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    modes = parser.add_mutually_exclusive_group(
        required=True
    )

    modes.add_argument(
        "--self-test",
        action="store_true",
    )

    modes.add_argument(
        "--audit",
        action="store_true",
    )

    modes.add_argument(
        "--run",
        action="store_true",
    )

    parser.add_argument(
        "--offline",
        action="store_true",
        help="Analisi senza accesso alla rete",
    )

    parser.add_argument(
        "--project-root",
        type=Path,
        default=(
            Path(__file__)
            .resolve()
            .parent
            .parent
        ),
    )

    parser.add_argument(
        "--crosswalk",
        type=Path,
        default=None,
        help=(
            "TSV RefGen-NAM5 con "
            "provenienza documentata"
        ),
    )

    args = parser.parse_args()

    if args.self_test:

        self_test()
        return

    root = args.project_root.resolve()

    external = (
        root
        / "data/external/"
        "ai_primetta_step47_b73"
    )

    runs = (
        root
        / "target_prediction_v33_vscode/"
        "runs"
    )

    if (
        not (root / "data/processed").is_dir()
        or not runs.is_dir()
    ):

        raise RuntimeError(
            "Progetto AI PRIMETTA V3.3 "
            "non riconosciuto"
        )

    print("=" * 68)

    print(
        "AI PRIMETTA V3.3 | STEP47 "
        "- PARE B73 + SICCITA B73"
    )

    print("=" * 68)

    print("[PROGETTO]", root)

    for key, value in SOURCES.items():

        print(
            "[FONTE]",
            key,
            "|",
            value["doi"],
            "|",
            value["context"],
        )

    print(
        "[CANDIDATE]",
        ", ".join(CANDIDATES),
    )

    print(
        "[VINCOLI] Nessun training; "
        "GUI, Docker e checkpoint invariati"
    )

    print(
        "[VINCOLI] PARE non dimostra "
        "regolazione specifica sotto siccita"
    )

    # ========================================================
    # 16. MODALITA' AUDIT
    # ========================================================

    if args.audit:

        print(
            "[CACHE]",
            external,
            "presente=",
            external.exists(),
        )

        if external.exists():

            for path in sorted(
                external.rglob("*")
            ):

                if path.is_file():

                    print(
                        "[FILE]",
                        path.relative_to(root),
                        path.stat().st_size,
                        "byte",
                    )

        print(
            "[OK] Audit in sola lettura"
        )

        return

    # ========================================================
    # 17. ACQUISIZIONE DELLE FONTI
    # ========================================================

    errors = []
    inventory = []
    input_files = []


    # Acquisizione PARE B73 tramite PMC o Europe PMC.
    acquire_pare_b73(
        external=external,
        offline=args.offline,
        inventory=inventory,
        input_files=input_files,
        errors=errors
    )
    try:

        metadata, media = figshare_media(
            offline=args.offline,
            root=external,
        )

        print(
            "[FIGSHARE_FILES]",
            len(media),
        )

        for item in media:

            name = Path(
                item.get("name", "")
            ).name

            url = (
                item.get("download_url")
                or (
                    "https://api.figshare.com/"
                    "v2/file/download/"
                    + str(item.get("id"))
                )
            )

            if (
                not name
                or not is_supported(
                    Path(name)
                )
            ):
                continue

            destination = (
                external
                / "FIGSHARE"
                / name
            )

            if (
                args.offline
                and not destination.exists()
            ):

                errors.append(
                    "DROUGHT: file assente offline "
                    + name
                )

                continue

            try:

                if args.offline:

                    info = {
                        "status": "cached",
                        "path": str(destination),
                        "sha256": sha(destination),
                        "bytes":
                            destination.stat().st_size,
                    }

                else:

                    info = obtain(
                        url,
                        destination,
                        md5=(
                            item.get("supplied_md5")
                            or item.get("md5")
                        ),
                    )

                info.update({
                    "source": "DROUGHT_B73_2020",
                    "file": destination.name,
                })

                inventory.append(info)

                if (
                    destination.suffix.lower()
                    == ".zip"
                ):

                    with zipfile.ZipFile(
                        destination
                    ) as archive:

                        for member in archive.infolist():

                            internal_name = Path(
                                member.filename
                            ).name

                            if not internal_name:
                                continue

                            internal_path = Path(
                                internal_name
                            )

                            if not is_supported(
                                internal_path
                            ):
                                continue

                            if (
                                internal_path.suffix.lower()
                                == ".zip"
                            ):
                                continue

                            if member.file_size > MAX_FILE:
                                continue

                            extracted = (
                                external
                                / "FIGSHARE"
                                / "unzipped"
                                / internal_name
                            )

                            if extracted.exists():

                                if (
                                    extracted.stat().st_size
                                    != member.file_size
                                ):

                                    raise RuntimeError(
                                        "Collisione file ZIP: "
                                        + internal_name
                                    )

                            else:

                                extracted.parent.mkdir(
                                    parents=True,
                                    exist_ok=True,
                                )

                                with archive.open(
                                    member
                                ) as source_stream:

                                    with extracted.open(
                                        "xb"
                                    ) as target_stream:

                                        import shutil

                                        shutil.copyfileobj(
                                            source_stream,
                                            target_stream,
                                        )

                            input_files.append((
                                "DROUGHT_B73_2020",
                                extracted,
                            ))

                else:

                    input_files.append((
                        "DROUGHT_B73_2020",
                        destination,
                    ))

                print(
                    "[ACQUISITO]",
                    "DROUGHT_B73_2020",
                    destination.name,
                    info["status"],
                )

            except Exception as exc:

                errors.append(
                    "DROUGHT file "
                    + name
                    + ": "
                    + str(exc)
                )

    except Exception as exc:

        errors.append(
            "DROUGHT: " + str(exc)
        )

    if not input_files:

        raise RuntimeError(
            "Nessun file sorgente analizzabile: "
            + "; ".join(errors)
        )

    # ========================================================
    # 19. CREAZIONE DELLA CARTELLA RISULTATI
    # ========================================================

    out = (
        runs
        / (
            "STEP47_FONTI_B73_UNICO_"
            + now()
        )
    )

    out.mkdir(
        parents=True,
        exist_ok=False,
    )

    # ========================================================
    # 20. ESTRAZIONE E ANALISI
    # ========================================================

    (
        profiles,
        matches,
        evidence,
        failures,
    ) = parse_collection(
        input_files,
        out_dir=out,
    )

    check_crosswalk(
        evidence,
        args.crosswalk,
    )

    write_tsv(
        out / "MATCH_ESATTI_2_SRNA.tsv",
        matches,
        [
            "source",
            "file",
            "sheet",
            "row",
            "cell",
            "candidate",
            "exact_sequence",
            "status",
        ],
    )

    evidence_columns = ['source', 'sheet', 'row', 'mirna_id', 'mirna_sequence', 'target_raw', 'target_gene_id', 'target_gene_build', 'nam5_gene_id', 'nam5_mapping_status', 'pare_category', 'tissue_or_stage', 'pare_site_id', 'pare_pvalue_raw', 'pare_allen_score', 's5_ValueP1', 's5_MedianFlank2', 's5_TotalValue3', 'psRNATarget_score_raw', 'interaction_padj_raw', 'target_position_raw', 'drought_regulation_raw', 'correlation_r', 'correlation_direction', 'evidence_type', 'evidence_level', 'candidate_exact_match', 'not_validated_in_PRJNA419326']

    write_tsv(
        out / "ASSOCIAZIONI_EVIDENZE_DA_VERIFICARE.tsv",
        evidence,
        evidence_columns,
    )

    # ========================================================
    # 21. CLASSIFICAZIONE DELLE EVIDENZE
    # ========================================================

    s4 = [
        item
        for item in evidence
        if item["evidence_type"]
        == "PARE_cleavage_reported_in_B73_non_drought"
    ]

    drought = [
        item
        for item in evidence
        if item["source"]
        == "DROUGHT_B73_2020"
    ]

    pare_levels = collections.Counter(
        item["evidence_level"]
        for item in s4
    )

    correlation_directions = (
        collections.Counter(
            item["correlation_direction"]
            for item in drought
        )
    )

    source_pairs = (
        collections.defaultdict(set)
    )

    for item in evidence:

        if (
            item["mirna_id"]
            and item["target_gene_id"]
        ):

            source_pairs[
                item["source"]
            ].add((
                tidy(item["mirna_id"]),
                item["target_gene_id"],
            ))

    external_overlap = (
        source_pairs["PARE_B73_2022"]
        & source_pairs["DROUGHT_B73_2020"]
    )

    # ========================================================
    # 22. CONFRONTO DELLE DUE SMALL RNA
    # ========================================================

    comparison = []

    for candidate, sequence in CANDIDATES.items():

        events = [
            item
            for item in matches
            if item["candidate"]
            == candidate
        ]

        direct_records = [
            item
            for item in evidence
            if item["candidate_exact_match"]
            == candidate
        ]

        comparison.append({
            "candidate": candidate,
            "sequence_RNA": sequence,
            "external_exact_cell_matches":
                len(events),
            "pare_rows_with_direct_sequence":
                sum(
                    item["source"]
                    == "PARE_B73_2022"
                    for item in direct_records
                ),
            "drought_rows_with_direct_sequence":
                sum(
                    item["source"]
                    == "DROUGHT_B73_2020"
                    for item in direct_records
                ),
            "drought_B73_H21_function_validated":
                "NO",
            "comment":
                "Match esatto non prova "
                "biogenesi, bersaglio o funzione "
                "sotto siccita",
        })

    write_tsv(
        out / "CONFRONTO_DUE_CANDIDATE.tsv",
        comparison,
        [
            "candidate",
            "sequence_RNA",
            "external_exact_cell_matches",
            "pare_rows_with_direct_sequence",
            "drought_rows_with_direct_sequence",
            "drought_B73_H21_function_validated",
            "comment",
        ],
    )

    # ========================================================
    # 23. SOVRAPPOSIZIONE FRA LE DUE FONTI
    # ========================================================

    overlap_rows = [
        {
            "mirna_identifier": mirna,
            "target_gene_same_id": gene,
            "interpretation":
                "same_string_only_review_"
                "genome_assembly_context",
        }
        for mirna, gene
        in sorted(external_overlap)
    ]

    write_tsv(
        out / "OVERLAP_FRA_FONTI_SOLO_ID_IDENTICI.tsv",
        overlap_rows,
        [
            "mirna_identifier",
            "target_gene_same_id",
            "interpretation",
        ],
    )

    status = collections.Counter(
        item["evidence_type"]
        for item in evidence
    )

    # ========================================================
    # 24. CONTROLLI SCIENTIFICI FINALI
    # ========================================================

    issues = []

    if not any(
        source == "PARE_B73_2022"
        for source, _ in input_files
    ):

        issues.append(
            "Manca supplemento PARE B73"
        )

    if not any(
        source == "DROUGHT_B73_2020"
        for source, _ in input_files
    ):

        issues.append(
            "Mancano supplementi siccita B73"
        )

    if not any(
        item["source"] == "PARE_B73_2022"
        and item["section"] == "S4"
        for item in profiles
    ):

        issues.append(
            "Dataset S4 PARE non identificato: "
            "verificare struttura Excel"
        )

    if not any(
        item["source"] == "PARE_B73_2022"
        and item["section"] == "S5"
        for item in profiles
    ):

        issues.append(
            "Dataset S5 non identificato"
        )

    if not any(
        item["source"] == "PARE_B73_2022"
        for item in evidence
    ):

        issues.append(
            "Nessuna coppia PARE estratta "
            "con schema riconosciuto"
        )

    if not any(
        item["source"] == "DROUGHT_B73_2020"
        for item in evidence
    ):

        issues.append(
            "Nessuna coppia siccita estratta "
            "con schema riconosciuto"
        )

    if any(
        item["nam5_mapping_status"]
        == "unmapped_no_assembly_equivalence_assumed"
        for item in evidence
    ):

        issues.append(
            "Mapping RefGen/NAM5 "
            "non dimostrato per alcune coppie"
        )

    if failures:

        issues.append(
            "Alcuni file non analizzati: "
            "vedere parse_failures"
        )

    if errors:

        issues.append(
            "Alcune fonti non acquisite: "
            "vedere acquisition_errors"
        )

    # ========================================================
    # 25. RAPPORTO JSON
    # ========================================================

    report = {
        "step": 47,
        "completed_utc": now(),
        "project": str(root),
        "sources": SOURCES,
        "inventory": inventory,
        "sheets": profiles,
        "n_exact_candidate_matches":
            len(matches),
        "n_potential_evidence_pairs":
            len(evidence),
        "by_evidence_type":
            dict(status),
        "s4_pare_category_counts":
            dict(pare_levels),
        "drought_correlation_directions":
            dict(correlation_directions),
        "same_identifier_overlap_between_sources_not_assembly_verified":
            len(overlap_rows),
        "candidate_comparison":
            comparison,
        "acquisition_errors":
            errors,
        "parse_failures":
            failures,
        "limitations":
            issues,
        "assumptions": [
            "Nessuna associazione computazionale "
            "trasformata in interazione validata",

            "PARE B73 in altri tessuti non "
            "dimostra funzione sotto siccita",

            "Correlazione sRNA-gene "
            "non prova taglio",

            "Nessuna equivalenza genomica "
            "presunta tra versioni B73",

            "Mancanza di match esatto "
            "non dimostra assenza di funzione",

            "Non usare direttamente questi "
            "file come etichette positive "
            "o negative per ML",
        ],
        "scientific_status":
            "PRELIMINARY_NEEDS_MANUAL_"
            "SCHEMA_AND_EVIDENCE_REVIEW",
        "preserved": {
            "GUI": True,
            "Docker": True,
            "historical_checkpoints": True,
        },
    }


    # STEP47_SCHEMI_B73_V5
    # Conteggi separati per tipo di evidenza.

    pare_s4 = [
        item for item in evidence
        if item["evidence_type"] ==
        "PARE_cleavage_reported_in_B73_non_drought"
    ]

    pare_s5 = [
        item for item in evidence
        if item["evidence_type"] ==
        "PARE_cleavage_efficiency_S5"
    ]

    s7_evidence = [
        item for item in evidence
        if item["source"] == "DROUGHT_B73_2020"
    ]

    def pare_site_key(item):
        return (
            tidy(item["tissue_or_stage"]),
            item["pare_site_id"].strip(),
            tidy(item["mirna_id"]).removeprefix("zma"),
            item["target_gene_id"],
        )

    s4_sites = {
        pare_site_key(item)
        for item in pare_s4
    }

    s5_shared = sum(
        pare_site_key(item) in s4_sites
        for item in pare_s5
    )

    report["schema_parsing"] = (
        "validated_column_names_from_original_workbooks"
    )

    report["pare_s4_rows"] = len(pare_s4)
    report["pare_s5_rows"] = len(pare_s5)

    report["pare_s5_rows_matching_s4_site_key"] = (
        s5_shared
    )

    report["drought_s7_association_rows"] = (
        len(s7_evidence)
    )

    report["cross_source_overlap_interpretable"] = False

    report["scientific_status"] = (
        "PRELIMINARY_EXTERNAL_EVIDENCE_NO_TRAINING"
    )

    print("[PARE_S4_RIGHE]", len(pare_s4))
    print("[PARE_S5_RIGHE]", len(pare_s5))
    print("[S5_SITI_MATCH_STRICT_S4]", s5_shared)
    print("[S7_ASSOCIAZIONI_RIGHE]", len(s7_evidence))

    put_json(
        out / "RAPPORTO_STEP47.json",
        report,
    )

    # ========================================================
    # 26. RAPPORTO SCIENTIFICO MARKDOWN
    # ========================================================

    markdown = [
        "# AI PRIMETTA V3.3 — STEP 47",
        "",
        "## Fonti: PARE B73 2022 e siccita B73 2020",
        "",
        f"- File acquisiti: {len(inventory)}",
        f"- Fogli/tabelle estratti: {len(profiles)}",
        (
            "- Righe candidate di interazioni "
            f"estratte: {len(evidence)}"
        ),
        (
            "- Match esatti delle due sRNA: "
            f"{len(matches)}"
        ),
        "",
        "## Categorie PARE S4",
    ]

    markdown.extend(
        f"- {key}: {value}"
        for key, value
        in pare_levels.items()
    )

    markdown.extend([
        "",
        "## Correlazioni",
    ])

    markdown.extend(
        f"- {key}: {value}"
        for key, value
        in correlation_directions.items()
    )

    markdown.extend([
        "",
        (
            "- Overlap per identificativi "
            f"identici: {len(overlap_rows)}"
        ),
        "",
        "## Tipologie di evidenza",
    ])

    markdown.extend(
        f"- {key}: {value}"
        for key, value
        in status.items()
    )

    markdown.extend([
        "",
        "## Limiti e verifiche necessarie",
    ])

    markdown.extend(
        f"- {issue}"
        for issue in issues
    )

    markdown.extend([
        "",
        "Nessun training.",
        "Nessun dato inventato.",
        "Nessun checkpoint modificato.",
        "GUI e Docker invariati.",
    ])

    report_md = (
        out / "RAPPORTO_STEP47.md"
    )

    report_md.write_text(
        "\n".join(markdown) + "\n",
        encoding="utf-8",
    )

    # ========================================================
    # 27. OUTPUT FINALE
    # ========================================================

    print(
        "[FOGLI]",
        len(profiles),
    )

    print(
        "[RIGHE_EVIDENZA]",
        len(evidence),
    )

    print(
        "[MATCH_ESATTI]",
        len(matches),
    )

    for item in profiles:

        print(
            "[FOGLIO]",
            item["source"],
            item["file"],
            item["sheet"],
            "section=",
            item["section"],
            "rows=",
            item["rows"],
            "fields=",
            item["recognized_columns"],
        )

    for issue in issues:

        print(
            "[ATTENZIONE]",
            issue,
        )

    print(
        "[REPORT]",
        report_md,
    )

    print(
        "[JSON]",
        out / "RAPPORTO_STEP47.json",
    )

    print(
        "[OK] Pipeline unica conclusa"
    )

    print(
        "[OK] Nessun training; "
        "GUI, Docker e checkpoint invariati"
    )


# ============================================================
# 28. AVVIO
# ============================================================
    step47_qc_v6(out)
    step47_qc_v7(out)

if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print(
            "[ERRORE]",
            type(exc).__name__ + ":",
            str(exc),
            file=sys.stderr,
        )

        sys.exit(1)
