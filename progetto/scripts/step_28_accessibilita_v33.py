
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""AI PRIMETTA V3.3 | STEP 28 - RNAplfold accessibilita' locale.

--audit: controllo integrita' e test interno RNAplfold, senza output nel progetto.
--run:   calcolo su 37 trascritti, misura descrittiva sui 36 siti delle 2 sRNA.

Nessuna modifica a GUI, Docker, modelli o dati originali.
"""

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / 'target_prediction_v33_vscode' / 'runs'
SRNAS = {'sRNA_00000037': 25, 'sRNA_00000435': 24}
UNPAIRED_LENGTH = 10
WINDOW = 150
SPAN = 100


# ============================================================
# 1. GESTIONE ERRORI
# ============================================================

def fail(message):
    raise ValueError(message)


# ============================================================
# 2. PROTEZIONE DEI PERCORSI
# ============================================================

def safe(path):
    p = Path(path)
    if not p.is_absolute():
        p = ROOT / p
    p = p.resolve()
    if not p.is_file() or not p.is_relative_to(RUNS.resolve()):
        fail(f'File non valido o esterno a runs: {p}')
    return p


# ============================================================
# 3. CHECKSUM SHA256
# ============================================================

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


# ============================================================
# 4. LETTURA FASTA
# ============================================================

def read_fasta(path):
    out = {}
    name, parts = None, []

    for line in path.read_text(
        encoding='utf-8-sig'
    ).splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith('>'):

            if name is not None:
                if not parts or name in out:
                    fail(f'FASTA duplicato/vuoto: {name}')

                out[name] = (
                    ''.join(parts)
                    .upper()
                    .replace('T', 'U')
                )

            name = line[1:].split()[0]

            if not name or name in out:
                fail('Identificativo FASTA invalido')

            parts = []

        else:
            if name is None:
                fail('FASTA senza intestazione')

            parts.append(line)

    if name is not None:
        if not parts or name in out:
            fail('Ultimo record FASTA invalido')

        out[name] = (
            ''.join(parts)
            .upper()
            .replace('T', 'U')
        )

    if not out:
        fail('FASTA senza sequenze')

    return out


# ============================================================
# 5. LETTURA TSV
# ============================================================

def read_tsv(path):
    with path.open(
        encoding='utf-8-sig',
        newline=''
    ) as f:

        reader = csv.DictReader(
            f,
            delimiter='\t'
        )

        fields = reader.fieldnames or []

        if not fields or len(fields) != len(set(fields)):
            fail(f'Intestazioni TSV invalide: {path}')

        return list(reader)


# ============================================================
# 6. RECUPERO E VALIDAZIONE STEP24
# ============================================================

def prepare_inputs():

    found = sorted(
        RUNS.glob(
            'STEP24_CONFRONTO_RNAPLEX_CORRETTO_*/'
            'RAPPORTO_STEP24.json'
        )
    )

    if not found:
        fail('Non trovato RAPPORTO_STEP24.json corretto')

    r24_path = safe(found[-1])

    r24 = json.loads(
        r24_path.read_text(
            encoding='utf-8'
        )
    )

    if (
        str(r24.get('step')) != '24'
        or r24.get('RNAplex_pairs') != 74
    ):
        fail('Rapporto STEP24 non coerente')

    r23_path = safe(
        r24['source_step23b']
    )

    r23 = json.loads(
        r23_path.read_text(
            encoding='utf-8'
        )
    )

    manifest = safe(
        r23['source_manifest']
    )

    sources = {
        'srna': safe(
            manifest.parent / 'SRNA_STANDARD_MAX25.fa'
        ),
        'targets': safe(
            manifest.parent / 'TRASCRITTI_B73_CANDIDATI.fa'
        ),
        'sites': safe(
            manifest.parent / 'SITI_E_RNASEQ_SNAPSHOT.tsv'
        ),
        'raw': safe(
            r23_path.parent / 'OUTPUT_GREZZO_RNAPLEX.txt'
        ),
        'table': safe(
            r23_path.parent / 'INTERAZIONI_RNAPLEX.tsv'
        ),
    }

    # ========================================================
    # 7. CONTROLLI CHECKSUM
    # ========================================================

    for key, file in sources.items():
        if sha(file) != r24['source_sha256'][key]:
            fail(f'Checksum STEP24 non coerente: {key}')

    srna = read_fasta(
        sources['srna']
    )

    targets = read_fasta(
        sources['targets']
    )

    if (
        {k: len(v) for k, v in srna.items()} != SRNAS
        or len(targets) != 37
    ):
        fail('Attesi esattamente 2 sRNA e 37 trascritti')

    if any(
        set(s) - set('ACGU')
        for s in srna.values()
    ):
        fail('sRNA con basi non canoniche')

    if any(
        set(s) - set('ACGUN') or len(s) < 50
        for s in targets.values()
    ):
        fail('Trascritti B73 non validi')

    # ========================================================
    # 8. VALIDAZIONE DEI 37 SITI ORIGINALI
    # ========================================================

    sites = read_tsv(
        sources['sites']
    )

    if (
        len(sites) != 37
        or len({
            s['transcript_id']
            for s in sites
        }) != 37
    ):
        fail('Attesi 37 siti su 37 trascritti distinti')

    candidates = []

    for row in sites:

        sid = row['srna_id']
        tx = row['transcript_id']

        if tx not in targets:
            fail(f'Trascritto STEP22 assente: {tx}')

        a = int(
            row['start_1based']
        )

        b = int(
            row['end_1based']
        )

        seq = (
            row['full_site_DNA_5to3']
            .upper()
            .replace('T', 'U')
        )

        if (
            not (1 <= a <= b <= len(targets[tx]))
            or targets[tx][a-1:b] != seq
        ):
            fail(f'Coordinate/sito difformi: {sid} {tx}')

        if sid in srna:

            if b-a+1 != len(srna[sid]):
                fail(
                    f'Sito con lunghezza inattesa: {sid} {tx}'
                )

            candidates.append(row)

    if (
        len(candidates) != 36
        or Counter(
            r['srna_id']
            for r in candidates
        ) != Counter({
            'sRNA_00000037': 17,
            'sRNA_00000435': 19
        })
    ):
        fail('Attesi 36 siti per le due sRNA')

    # ========================================================
    # 9. VALIDAZIONE RISULTATI RNAplex STEP24
    # ========================================================

    r24_data = read_tsv(
        safe(r24['output_tsv'])
    )

    bypair = {}

    for row in r24_data:

        key = (
            row['srna_id'],
            row['transcript_id']
        )

        if key in bypair:
            fail('Coppia RNAplex ripetuta')

        bypair[key] = row

    expected_pairs = {
        (s, t)
        for s in srna
        for t in targets
    }

    if (
        len(bypair) != 74
        or set(bypair) != expected_pairs
    ):
        fail('Coppie STEP24 incomplete')

    return (
        r24_path,
        sources,
        targets,
        candidates,
        bypair
    )


# ============================================================
# 10. LETTURA DEI PROFILI RNAplfold
# ============================================================

def lunp_probabilities(path, sequence_length):

    """
    p10[i] = probabilita' che i nucleotidi
    [i-9 ... i] siano tutti non appaiati.

    Coordinate i in base 1.
    """

    probs = {}
    all_positions = set()

    for text in path.read_text(
        encoding='utf-8-sig'
    ).splitlines():

        text = text.strip()

        if not text or text.startswith('#'):
            continue

        fields = text.split()

        if not fields[0].isdigit():
            continue

        i = int(fields[0])

        if (
            i < 1
            or i > sequence_length
            or i in all_positions
        ):
            fail(
                f'Posizione _lunp duplicata '
                f'o invalida: {i}'
            )

        all_positions.add(i)

        if len(fields) <= UNPAIRED_LENGTH:
            fail(
                f'Colonna p{UNPAIRED_LENGTH} '
                f'assente in riga {i}'
            )

        token = fields[UNPAIRED_LENGTH]

        if token.upper() != 'NA':

            try:
                value = float(token)

            except ValueError as exc:
                raise ValueError(
                    f'Probabilita non numerica: {token}'
                ) from exc

            if (
                not math.isfinite(value)
                or value < 0
                or value > 1.000001
            ):
                fail(
                    f'Probabilita fuori [0,1] '
                    f'alla riga {i}: {token}'
                )

            probs[i] = min(value, 1.0)

    if len(all_positions) != sequence_length:
        fail(
            f'_lunp incompleto: '
            f'{len(all_positions)}/{sequence_length} '
            'posizioni'
        )

    if not probs:
        fail('Nessuna probabilita p10 disponibile')

    return probs


# ============================================================
# 11. ESECUZIONE RNAplfold
# ============================================================

def launch_plfold(
    executable,
    tx,
    sequence,
    workdir
):

    w = min(
        WINDOW,
        len(sequence)
    )

    span = min(
        SPAN,
        w
    )

    command = [
        executable,
        '-W', str(w),
        '-L', str(span),
        '-u', str(UNPAIRED_LENGTH)
    ]

    result = subprocess.run(
        command,
        input=f'>{tx}\n{sequence}\n',
        text=True,
        cwd=str(workdir),
        capture_output=True,
        timeout=240,
        check=False
    )

    if result.returncode != 0:
        fail(
            f'RNAplfold fallito per {tx}: '
            + result.stderr[:400]
        )

    outputs = list(
        workdir.glob('*_lunp')
    )

    if len(outputs) != 1:
        fail(
            f'RNAplfold ha generato '
            f'{len(outputs)} file _lunp per {tx}; '
            'controllare il programma'
        )

    return (
        outputs[0],
        lunp_probabilities(
            outputs[0],
            len(sequence)
        ),
        command
    )


# ============================================================
# 12. AUDIT
# ============================================================

def audit(executable):

    (
        r24,
        sources,
        targets,
        candidates,
        bypair
    ) = prepare_inputs()

    version_result = subprocess.run(
        [
            executable,
            '--version'
        ],
        capture_output=True,
        text=True,
        timeout=15,
        check=True
    )

    version = (
        version_result.stdout
        or version_result.stderr
    ).strip()

    print(
        '[VERSIONE]',
        version
    )

    print(
        '[PASS] Integrita STEP22/24; '
        '2 sRNA, 37 trascritti, 36 siti: OK'
    )

    # ========================================================
    # 13. TEST CON RNA SINTETICO
    # ========================================================

    test_sequence = (
        'ACGUGCAUUG' * 11
    )[:100]

    with tempfile.TemporaryDirectory(
        prefix='primetta28_test_'
    ) as d:

        output, probs, cmd = launch_plfold(
            executable,
            'test_sintetico',
            test_sequence,
            Path(d)
        )

        if not all(
            i in probs
            for i in range(
                10,
                len(test_sequence)+1
            )
        ):
            fail(
                'Test RNAplfold: probabilita '
                'p10 non complete'
            )

    print(
        '[PASS] Test RNAplfold _lunp: '
        'coordinate e probabilita p10 valide'
    )

    return (
        r24,
        sources,
        targets,
        candidates,
        bypair,
        version
    )


# ============================================================
# 14. ANALISI ACCESSIBILITA'
# ============================================================

def run_analysis(executable, state):

    (
        r24,
        sources,
        targets,
        candidates,
        bypair,
        version
    ) = state

    stamp = datetime.now(
        timezone.utc
    ).strftime(
        '%Y%m%dT%H%M%S%fZ'
    )

    folder = RUNS / (
        f'STEP28_ACCESSIBILITA_RNAPLFOLD_{stamp}'
    )

    folder.mkdir(
        parents=True,
        exist_ok=False
    )

    rawfolder = (
        folder / 'PROFILI_LUNP_ORIGINALI'
    )

    rawfolder.mkdir()

    profiles = {}
    commands = {}

    # ========================================================
    # 15. CALCOLO DEI 37 PROFILI
    # ========================================================

    for ix, (tx, sequence) in enumerate(
        sorted(targets.items()),
        1
    ):

        with tempfile.TemporaryDirectory(
            prefix='primetta28_tx_'
        ) as tmp:

            lunp_file, profile, cmd = launch_plfold(
                executable,
                tx,
                sequence,
                Path(tmp)
            )

            dest = (
                rawfolder / f'{tx}_lunp'
            )

            if dest.exists():
                fail(
                    f'Profilo duplicato: {tx}'
                )

            shutil.copyfile(
                lunp_file,
                dest
            )

            profiles[tx] = profile

            commands[tx] = {
                'command': cmd[1:],
                'raw_sha256': sha(dest),
                'raw': str(
                    dest.relative_to(ROOT)
                )
            }

        print(
            f'[PROFILO] {ix}/37 '
            f'| {tx} '
            f'| lunghezza={len(sequence)}',
            flush=True
        )

    # ========================================================
    # 16. ANALISI DEI 36 SITI
    # ========================================================

    data = []

    for site in candidates:

        sid = site['srna_id']
        tx = site['transcript_id']

        a = int(
            site['start_1based']
        )

        b = int(
            site['end_1based']
        )

        # Ogni finestra di 10 nt deve essere
        # interamente contenuta nel sito originale.

        ends = list(
            range(
                a + UNPAIRED_LENGTH - 1,
                b + 1
            )
        )

        available = [
            profiles[tx][i]
            for i in ends
            if i in profiles[tx]
        ]

        if not available:
            fail(
                'Sito senza finestre accessibili '
                f'misurate: {sid} {tx}'
            )

        m = bypair[
            (sid, tx)
        ]

        data.append({

            'srna_id':
                sid,

            'transcript_id':
                tx,

            'gene_id':
                site['gene_id'],

            'target_start_1based':
                a,

            'target_end_1based':
                b,

            'windows10_expected':
                len(ends),

            'windows10_available':
                len(available),

            'measurement_complete':
                (
                    'yes'
                    if len(available) == len(ends)
                    else 'no'
                ),

            'mean_P_unpaired_10nt':
                f'{statistics.mean(available):.9g}',

            'median_P_unpaired_10nt':
                f'{statistics.median(available):.9g}',

            'min_P_unpaired_10nt':
                f'{min(available):.9g}',

            'max_P_unpaired_10nt':
                f'{max(available):.9g}',

            'RNAplex_MFE_kcal_mol':
                m['energy_kcal_mol'],

            'STEP24_comparison':
                m['comparison'],

            'RNAplfold_temperature_C':
                37,
        })

    # ========================================================
    # 17. ORDINAMENTO DESCRITTIVO
    # ========================================================

    ordered = sorted(
        data,
        key=lambda r: (
            -float(
                r['mean_P_unpaired_10nt']
            ),
            r['srna_id'],
            r['transcript_id']
        )
    )

    output = (
        folder / 'ACCESSIBILITA_36_SITI.tsv'
    )

    with output.open(
        'w',
        encoding='utf-8',
        newline=''
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(ordered[0]),
            delimiter='\t'
        )

        writer.writeheader()
        writer.writerows(ordered)

    # ========================================================
    # 18. STATISTICHE PER SMALL RNA
    # ========================================================

    stats = {}

    for sid in SRNAS:

        items = [
            r
            for r in data
            if (
                r['srna_id'] == sid
                and r['measurement_complete'] == 'yes'
            )
        ]

        stats[sid] = {

            'sites_total':
                sum(
                    r['srna_id'] == sid
                    for r in data
                ),

            'sites_complete':
                len(items),

            'median_mean_p10':
                (
                    statistics.median(
                        float(
                            r['mean_P_unpaired_10nt']
                        )
                        for r in items
                    )
                    if items
                    else None
                ),
        }

        print(
            f'[RISULTATO] {sid}: '
            f'{stats[sid]}'
        )

    # ========================================================
    # 19. RAPPORTO SCIENTIFICO
    # ========================================================

    report = {

        'step': 28,

        'status':
            'EXPLORATORY_LOCAL_ACCESSIBILITY',

        'method':
            'RNAplfold unpaired probabilities, '
            'P10 windows within STEP22 site',

        'viennarna_version':
            version,

        'window_max':
            WINDOW,

        'span_max':
            SPAN,

        'unpaired_length':
            UNPAIRED_LENGTH,

        'source_step24':
            str(r24.relative_to(ROOT)),

        'source_sha256': {
            k: sha(v)
            for k, v in sources.items()
        },

        'transcripts_processed':
            len(profiles),

        'sites_processed':
            len(data),

        'sites_with_complete_p10':
            sum(
                r['measurement_complete'] == 'yes'
                for r in data
            ),

        'by_srna':
            stats,

        'profiles':
            commands,

        'output_tsv':
            str(output.relative_to(ROOT)),

        'limitations': [

            'Probabilita locali di segmenti '
            'da 10 nt non appaiati, non intero '
            'sito da 24-25 nt',

            'Non e una misura sperimentale '
            'di accessibilita in vivo',

            'Siti e trascritti gia preselezionati '
            'per complementarita',

            'Nessuna conferma di targeting '
            'o DE significativa',
        ],

        'gui_modified':
            False,

        'docker_modified':
            False,

        'training_performed':
            False,
    }

    reportfile = (
        folder / 'RAPPORTO_STEP28.json'
    )

    reportfile.write_text(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2
        ) + '\n',
        encoding='utf-8'
    )

    print(
        '[COPERTURA] misure complete:',
        report['sites_with_complete_p10'],
        '/36'
    )

    print(
        '[TSV]',
        output
    )

    print(
        '[REPORT]',
        reportfile
    )

    print(
        '[OK] STEP28 concluso; '
        'GUI, Docker e modelli invariati'
    )


# ============================================================
# 20. PROGRAMMA PRINCIPALE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    mode = parser.add_mutually_exclusive_group(
        required=True
    )

    mode.add_argument(
        '--audit',
        action='store_true'
    )

    mode.add_argument(
        '--run',
        action='store_true'
    )

    args = parser.parse_args()

    print('=' * 70)

    print(
        'AI PRIMETTA V3.3 | STEP 28 '
        '- ACCESSIBILITA RNAplfold'
    )

    print('=' * 70)

    executable = shutil.which(
        'RNAplfold'
    )

    if executable is None:
        fail(
            'RNAplfold assente. Verificare: '
            'which RNAplfold; '
            'sudo apt install vienna-rna'
        )

    state = audit(
        executable
    )

    if args.audit:

        print(
            '[OK] Audit completo, '
            'nessun file scritto nel progetto'
        )

    else:

        run_analysis(
            executable,
            state
        )


# ============================================================
# 21. AVVIO
# ============================================================

if __name__ == '__main__':

    try:
        main()

    except Exception as exc:

        print(
            f'[ERRORE] '
            f'{type(exc).__name__}: {exc}',
            file=sys.stderr
        )

        sys.exit(1)
