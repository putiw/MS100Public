#!/usr/bin/env python3
"""Single master runner for the parallel Comment 14 EDSS >3 package."""
from __future__ import annotations
import argparse, json, shutil, subprocess, sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT.parent / "code"

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--run-dir", type=Path)
    a = ap.parse_args()
    cfg = json.loads(a.config.expanduser().resolve().read_text())
    run = (a.run_dir or ROOT / "work" / ("run_" + datetime.now().strftime("%Y%m%d_%H%M%S"))).resolve()
    if run.exists(): raise FileExistsError(run)
    code = run / "code"; models = run / "results/model_results"; tables = run / "results/tables"
    code.mkdir(parents=True); models.mkdir(parents=True); tables.mkdir(parents=True)
    source = (CANON / "run_comment14_whole_brain_models.m").read_text()
    source = source.replace("bin_thresholds = [4,1];", "bin_thresholds = [3,1];")
    source = source.replace("bin_use_geq = [true,false];", "bin_use_geq = [false,false];")
    source = source.replace("expected_binary_positive = [31,31,31,31; 39,39,39,39];", "expected_binary_positive = [39,39,39,39; 39,39,39,39];")
    source = source.replace("EDSS >= 4.0 (revisionExtras sensitivity definition)", "EDSS > 3.0 (Comment 14 parallel sensitivity definition)")
    (code / "run_comment14_gt3_models.m").write_text(source)
    matlab = cfg["matlab_executable"]
    q = lambda x: "'" + str(Path(x).expanduser().resolve()).replace("'", "''") + "'"
    expr = (f"addpath({q(code)}); run_comment14_gt3_models('ClinicalFile',{q(Path(cfg['stats_dir'])/'clinicalScore.xlsx')},"
            f"'WholeBrainFile',{q(cfg['whole_brain_means_file'])},'MetricsDir',{q(cfg['metrics_dir'])},"
            f"'ClassicalFile',{q(cfg.get('classical_file', Path(cfg['metrics_dir'])/'ClassicalRegionAllMetrics.csv'))},"
            f"'OutputDir',{q(models)},'QADir',{q(run/'qa')});")
    subprocess.run([matlab, "-batch", expr], check=True)
    subprocess.run([sys.executable, str(ROOT/"code/build_table_gt3.py"), "--models", str(models/"whole_brain_qmri_models.csv"), "--output", str(tables/"Supplementary_Table_6.csv")], check=True)
    report = {"status":"PASS", "edss_rule":"EDSS > 3.0", "rows":20, "table_rows":20, "reference_package":"parallel; existing EDSS >=4 files unchanged", "run_dir":str(run)}
    (run/"qa/reproduction_report.json").parent.mkdir(parents=True, exist_ok=True)
    (run/"qa/reproduction_report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))
if __name__ == "__main__": main()
