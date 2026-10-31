#!/usr/bin/env python3
"""Reproduce the two seven-metric Supplementary Table 6 files."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
PACKAGE = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE / 'code'))
from extract_parenchymal_means import extract, sha256

OUTCOMES = ('EDSS', 'MSPro', 'T25FW', '9HPT-D', '9HPT-ND')
METRICS = ('T1', 'MTR', 'FA', 'MD', 'NDI', 'ODI', 'ISOVF')
EXPECTED_N = {'EDSS': 79, 'MSPro': 78, 'T25FW': 75, '9HPT-D': 78, '9HPT-ND': 78}


def rows(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def quote(path):
    return "'" + str(Path(path).resolve()).replace("'", "''") + "'"


def build_table(model_file, baseline_file, output, rule):
    models = rows(model_file)
    assert len(models) == 15
    original = rows(baseline_file)
    assert len(original) == 20
    assert {r['MRI Metric'] for r in original} == set(METRICS[:4])
    new = []
    for r in models:
        assert r['Status'] == 'COMPLETED' and float(r['Adjusted']) == 0
        assert r['Demographics'] == 'none'
        outcome = r['Outcome']
        n = int(float(r['N']))
        assert n == EXPECTED_N[outcome]
        if outcome in ('EDSS', 'MSPro'):
            positives = int(float(r['PositiveN']))
            assert positives == (31 if outcome == 'EDSS' and rule == 'ge4' else 39)
            n_text = f'{n} ({positives})'
            performance = (f"AUC {float(r['Performance']):.3f} "
                           f"({float(r['CI_lower']):.3f}-{float(r['CI_upper']):.3f})")
        else:
            n_text = str(n)
            performance = f"R2 {float(r['Performance']):.3f}"
        p = float(r['p_raw'])
        new.append({'Outcome': outcome,
                    'MRI Metric': 'ISOVF' if r['Metric'] == 'FWF' else r['Metric'],
                    'Model': 'Whole brain', 'N': n_text,
                    'Apparent performance': performance,
                    'AIC': f"{float(r['AIC']):.2f}",
                    'p (raw)': '< 0.0001' if p < .0001 else f'{p:.4f}'})
    lookup = {(r['Outcome'], r['MRI Metric']): r for r in original + new}
    assert len(lookup) == 35
    combined = [lookup[o, m] for o in OUTCOMES for m in METRICS]
    assert [r for r in combined if r['MRI Metric'] in METRICS[:4]] == original
    with output.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(original[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(combined)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path,
                        help='Aggregate table destination; defaults to results beside this script.')
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    required = ('cohort_file', 'freesurfer_root', 'transformed_maps', 'clinical_file',
                'metrics_dir', 'classical_file', 'baseline_table_gt3',
                'baseline_table_ge4', 'matlab_executable', 'work_dir')
    missing = set(required) - set(cfg)
    if missing:
        raise ValueError(f'Missing configuration keys: {sorted(missing)}')
    repo = PACKAGE.parents[1]
    work = Path(cfg['work_dir']).expanduser().resolve()
    if work.is_relative_to(repo):
        raise ValueError('work_dir must be outside the public repository.')
    if work.exists():
        raise FileExistsError(f'Choose a new work_dir: {work}')
    output = (args.output_dir or PACKAGE / 'results').expanduser().resolve()
    names = {rule: output / f'Supplementary_Table_6_EDSS_{rule}.csv'
             for rule in ('gt3', 'ge4')}
    if any(path.exists() for path in names.values()):
        raise FileExistsError('Choose a new output directory to preserve existing tables.')
    work.mkdir(parents=True, mode=0o700)
    output.mkdir(parents=True, exist_ok=True)
    means = extract(cfg, work)
    matlab_script = PACKAGE / 'code/run_comment14_noddi_whole_brain_models.m'
    manifest = {'entry_script_sha256': sha256(__file__),
                'matlab_script_sha256': sha256(matlab_script), 'tables': {}}
    for rule in ('gt3', 'ge4'):
        model_dir = work / rule / 'models'
        qa_dir = work / rule / 'qa'
        expr = (f"addpath({quote(PACKAGE / 'code')}); "
                f"run_comment14_noddi_whole_brain_models('ClinicalFile',{quote(cfg['clinical_file'])},"
                f"'WholeBrainFile',{quote(means)},'MetricsDir',{quote(cfg['metrics_dir'])},"
                f"'ClassicalFile',{quote(cfg['classical_file'])},'EDSSRule','{rule}',"
                f"'OutputDir',{quote(model_dir)},'QADir',{quote(qa_dir)});")
        print(f'Fitting EDSS {rule}', flush=True)
        with (work / f'matlab_{rule}.log').open('w') as log:
            subprocess.run([cfg['matlab_executable'], '-batch', expr],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        build_table(model_dir / 'whole_brain_qmri_models.csv',
                    cfg[f'baseline_table_{rule}'], names[rule], rule)
        manifest['tables'][rule] = sha256(names[rule])
        print(f'Wrote {names[rule]}', flush=True)
    (work / 'run_manifest_private.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
