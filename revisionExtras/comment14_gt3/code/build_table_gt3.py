#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv
ap=argparse.ArgumentParser(); ap.add_argument('--models',required=True); ap.add_argument('--output',required=True); a=ap.parse_args()
rows=list(csv.DictReader(open(a.models,newline=''))); metrics=['T1','MTR','FA','MD']; outcomes=['EDSS','MSPro','T25FW','9HPT-D','9HPT-ND']
out=[]
for o in outcomes:
  for m in metrics:
    r=next(x for x in rows if x['Outcome']==o and x['Metric']==m); n=int(float(r['N']))
    if o in ('EDSS','MSPro'):
      perf=f"AUC {float(r['Performance']):.3f} ({float(r['CI_lower']):.3f}-{float(r['CI_upper']):.3f})"; nt=f"{n} ({int(float(r['PositiveN']))})"
    else: perf=f"R2 {float(r['Performance']):.3f}"; nt=str(n)
    p=float(r['p_raw']); pt='< 0.0001' if p<.0001 else f'{p:.4f}'
    out.append({'Outcome':o,'MRI Metric':m,'Model':'Whole brain','N':nt,'Apparent performance':perf,'AIC':f"{float(r['AIC']):.2f}",'p (raw)':pt})
with open(a.output,'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
