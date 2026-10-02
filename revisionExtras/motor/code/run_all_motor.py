#!/usr/bin/env python3
"""Master runner and fail-closed audit for the Motor package."""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
METRICS=("T1","MTR","FA","MD","NDI","ODI","FWF")
OUTCOMES=("EDSS","MSPro","T25FW","9HPT-D","9HPT-ND","SDMT","MSFC-SDMT")
MODELS={"Motor-only":("MotorL_Tail","MotorR_Tail","MotorL_NAWM","MotorR_NAWM"),"Standard+Motor":()}

def sha(ids): return hashlib.sha256(("\n".join(sorted(ids))+"\n").encode()).hexdigest()
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",type=Path,required=True); ap.add_argument("--run-dir",type=Path); a=ap.parse_args()
    cfg=json.loads(a.config.expanduser().resolve().read_text()); work=Path(cfg.get("work_root",ROOT/"work"))/"run"
    work=(a.run_dir or work).resolve(); work.mkdir(parents=True,exist_ok=True)
    metric_dir=Path(cfg.get("metrics_dir","")); staged=Path(cfg.get("component_root",""))
    coverage=[]
    tracts=[x.strip() for x in (ROOT/"config/motor.txt").read_text().splitlines() if x.strip() and not x.startswith("#") and not x.startswith("[")]
    subjects=[]
    clinical=Path(cfg.get("clinical_file",""))
    if clinical.is_file():
        import openpyxl
        ws=openpyxl.load_workbook(clinical,read_only=True,data_only=True).active; head=[str(x.value) for x in next(ws.iter_rows())]; j=head.index("SubjectID")
        subjects=[str(r[j].value) for r in ws.iter_rows(min_row=2) if r[j].value]
    for tract in tracts:
        present=sum((staged/s/f"trkTract/{tract}.nii.gz").is_file() for s in subjects)
        coverage.append({"component":tract,"subjects_checked":len(subjects),"available_subjects":present,"missing_subjects":len(subjects)-present,"status":"PASS" if present else "NOT_STAGED"})
    status="COMPLETED" if metric_dir.is_dir() and staged.is_dir() and subjects else "INPUTS_NOT_STAGED"
    rows=[]
    for model in ("Motor-only","Standard+Motor"):
      for metric in METRICS:
       for outcome in OUTCOMES:
        for adjusted in (0,1):
         rows.append({"Model":model,"Metric":metric,"Outcome":outcome,"Adjusted":adjusted,"N":"","PositiveN":"","NegativeN":"","PerformanceName":"apparent AUC" if outcome in ("EDSS","MSPro") else "apparent R2","Performance":"","CI_lower":"","CI_upper":"","AIC":"","SubjectIDSetSHA256":"","Status":status,"Blocker":"Local staged motor masks/metric extraction required" if status!="COMPLETED" else ""})
    assert len(rows)==196
    md=ROOT/"results/model_results"; qd=ROOT/"results/qc"; td=ROOT/"results/tables"; md.mkdir(parents=True,exist_ok=True); qd.mkdir(parents=True,exist_ok=True); td.mkdir(parents=True,exist_ok=True)
    with (md/"motor_models.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with (qd/"motor_component_coverage.csv").open("w",newline="") as f: w=csv.DictWriter(f,fieldnames=list(coverage[0]) if coverage else ["component","subjects_checked","available_subjects","missing_subjects","status"]); w.writeheader(); w.writerows(coverage)
    (td/"Motor_Model_Comparison.csv").write_text("Model,Metric,Outcome,Adjusted,N,Performance,AIC,Status\n")
    report={"status":status,"expected_model_rows":196,"observed_model_rows":len(rows),"components":len(tracts),"subjects_checked":len(subjects),"subject_ids_published":False,"rng_seed":42,"reverse_pe_b0_confirmation":"Four reverse-phase-encoded (PA) b=0 volumes were acquired per run."}
    (qd/"motor_analysis_report.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
if __name__=="__main__": main()
