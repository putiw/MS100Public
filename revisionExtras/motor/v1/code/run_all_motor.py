#!/usr/bin/env python3
"""One user-facing runner: stage, extract, fit, summarize, and verify Motor models."""
from __future__ import annotations
import argparse,csv,hashlib,json,subprocess,sys
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
METRICS=('T1','MTR','FA','MD','NDI','ODI','FWF')
OUTCOMES=('EDSS','MSPro','T25FW','9HPT-D','9HPT-ND','SDMT','MSFC-SDMT')

def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def rows(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))

def write_csv(p,records,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(records)

def aggregate(cfg,extract_dir,model_file,output_dir):
 qc=json.loads((extract_dir/'subject_qc.json').read_text())
 tracts=[x.strip() for x in (ROOT/'config/motor.txt').read_text().splitlines() if x.strip() and not x.startswith(('#','['))]
 coverage=[]
 for tract in tracts:
  available=sum(tract in q['present'] for q in qc)
  empty=sum(tract in q['empty'] for q in qc)
  coverage.append({'Component':tract,'SubjectsChecked':len(qc),'Available':available,'Unavailable':len(qc)-available,'PresentButEmpty':empty})
 model_rows=rows(model_file)
 if len(model_rows)!=196:raise AssertionError(f'Expected 196 model rows, got {len(model_rows)}')
 bykey={(r['Metric'],r['Outcome'],r['Adjusted'],r['Model']):r for r in model_rows}
 if len(bykey)!=196:raise AssertionError('Duplicate model cell')
 comparisons=[]
 for m in METRICS:
  for o in OUTCOMES:
   for adjusted in ('0','1'):
    motor=bykey[(m,o,adjusted,'Motor-only')]; combined=bykey[(m,o,adjusted,'Standard+Motor')]
    if motor['N']!=combined['N'] or motor['SubjectIDSetSHA256']!=combined['SubjectIDSetSHA256']:raise AssertionError(f'Participant mismatch {m} {o} {adjusted}')
    comparisons.append({'Metric':'ISOVF' if m=='FWF' else m,'Outcome':o,'Adjusted':adjusted,'N':motor['N'],
      'PositiveN':motor['PositiveN'],'NegativeN':motor['NegativeN'],
      'PerformanceName':motor['PerformanceName'],'MotorOnly':motor['Performance'],'MotorOnly_CI_lower':motor['CI_lower'],'MotorOnly_CI_upper':motor['CI_upper'],'MotorOnly_AIC':motor['AIC'],
      'StandardPlusMotor':combined['Performance'],'StandardPlusMotor_CI_lower':combined['CI_lower'],'StandardPlusMotor_CI_upper':combined['CI_upper'],'StandardPlusMotor_AIC':combined['AIC'],
      'MotorOnly_Status':motor['Status'],'StandardPlusMotor_Status':combined['Status']})
 if len(comparisons)!=98:raise AssertionError('Expected 98 paired cells')
 available=[q for q in qc if q['status']=='PASS']
 unavailable=[q for q in qc if q['status']=='UNAVAILABLE_NO_NONEMPTY_COMPONENT']
 if len(available)+len(unavailable)!=len(qc) or not available:raise AssertionError('Unexpected Motor mask status')
 mask_pass=all(q['geometry_pass'] and q['split_exact'] and q['combined_voxels']>0 and q['left_voxels']>0 and q['right_voxels']>0 and q['left_voxels']+q['right_voxels']==q['combined_voxels'] and q['left_centroid_x']<q['right_centroid_x'] for q in available)
 if not mask_pass:raise AssertionError('Motor geometry/nonzero/split QC failed')
 missing_map_counts={m:sum(m in q['invalid_metric_maps'] for q in qc) for m in METRICS}
 result_dir=output_dir/'results'; result_dir.mkdir(parents=True)
 model_out=result_dir/'model_results/motor_models.csv';model_out.parent.mkdir();model_out.write_bytes(model_file.read_bytes())
 write_csv(result_dir/'tables/Motor_Model_Comparison.csv',comparisons,list(comparisons[0]))
 write_csv(result_dir/'qc/motor_component_coverage.csv',coverage,list(coverage[0]))
 status_counts=dict(Counter(r['Status'] for r in model_rows))
 report={'schema_version':1,'status':'PASS' if status_counts=={'COMPLETED':196} else 'PARTIAL',
  'model_rows':len(model_rows),'comparison_cells':len(comparisons),'model_status_counts':status_counts,
  'subjects_checked':len(qc),'motor_mask_available_subjects':len(available),
  'motor_mask_unavailable_subjects':len(unavailable),'component_count':len(tracts),
  'mask_geometry_nonzero_and_split_pass':mask_pass,
  'matched_subject_sets_within_all_cells':True,'min_combined_voxels':min(q['combined_voxels'] for q in available),
  'min_left_voxels':min(q['left_voxels'] for q in available),'min_right_voxels':min(q['right_voxels'] for q in available),
  'max_missing_components_per_subject':max(len(q['missing']) for q in qc),
  'missing_metric_map_subject_counts':missing_map_counts,
  'edss_rule':'EDSS > 3','mspro_rule':'MSPro > 1','rng_seed':42,'auc_bootstrap_resamples':2000,
  'noddi_zero_policy':'valid inside fitting support','legacy_missing_policy':'zero-valued legacy summaries excluded from model fitting',
  'reverse_pe_b0_confirmation':'Four reverse-phase-encoded (PA) b=0 volumes were acquired per run.',
  'subject_identifiers_published':False,'paired_delta_auc_computed':False}
 (result_dir/'qc/motor_analysis_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 manifest={'schema_version':1,'motor_definition_sha256':digest(ROOT/'config/motor.txt'),
  'aggregate_input_sha256':{k:digest(Path(p)) for k,p in {'clinical':str(Path(cfg['stats_dir'])/'clinicalScore.xlsx'),
   **{m:str(Path(cfg['metrics_dir'])/('GroupTract'+m+'_All'+('.xlsx' if m in ('T1','MTR','FA','MD') else '.csv'))) for m in METRICS}}.items()},
  'output_sha256':{p.name:digest(p) for p in [model_out,result_dir/'tables/Motor_Model_Comparison.csv',result_dir/'qc/motor_component_coverage.csv']}}
 (result_dir/'qc/input_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
 return report

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);ap.add_argument('--stage-dir',type=Path,required=True);ap.add_argument('--run-dir',type=Path,required=True);ap.add_argument('--reference-dir',type=Path);ap.add_argument('--workers',type=int,default=1);a=ap.parse_args()
 cfg=json.loads(a.config.read_text());scratch=Path(cfg['scratch_root']).resolve();run=a.run_dir.resolve();stage=a.stage_dir.resolve()
 if not run.is_relative_to(scratch) or not stage.is_relative_to(scratch) or run==stage:raise ValueError('Output and stage must be separate directories below local scratch')
 if run.exists():raise FileExistsError(run)
 subprocess.run([sys.executable,str(ROOT/'code/stage_extract_motor.py'),'--config',str(a.config),'--stage-dir',str(stage),'--output-dir',str(run/'extraction'),'--workers',str(a.workers)],check=True)
 model_dir=run/'models';model_dir.mkdir();model_file=model_dir/'motor_models.csv'
 q=lambda x:"'"+str(Path(x).resolve()).replace("'","''")+"'"
 expr=f"addpath({q(ROOT/'code')});run_motor_models({q(cfg['stats_dir'])},{q(cfg['metrics_dir'])},{q(run/'extraction/metrics')},{q(model_file)});"
 with (run/'matlab.log').open('x') as log:
  subprocess.run([cfg['matlab_executable'],'-batch',expr],stdout=log,stderr=subprocess.STDOUT,check=True)
 print('MATLAB Motor fitting completed',flush=True)
 report=aggregate(cfg,run/'extraction',model_file,run)
 if a.reference_dir:
  reference=a.reference_dir.resolve()
  generated=run/'results'
  for name in ('model_results/motor_models.csv','tables/Motor_Model_Comparison.csv','qc/motor_component_coverage.csv','qc/motor_analysis_report.json','qc/input_manifest.json'):
   if digest(generated/name)!=digest(reference/name):raise AssertionError(f'Regenerated result differs: {name}')
  report['matches_packaged_reference']=True
 print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
