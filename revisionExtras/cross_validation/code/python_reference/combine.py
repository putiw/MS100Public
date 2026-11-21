"""Combine per-metric outputs of run.py into one aggregate file."""
import os, sys, pandas as pd
out = os.environ.get('MS100_CV_OUT', '.')
ms = ['T1', 'MTR', 'FA', 'MD', 'NDI', 'ODI', 'FWF']
d = pd.concat([pd.read_csv(os.path.join(out, f'cv_results_{m}.csv')) for m in ms])
d.to_csv(os.path.join(out, 'cv_results_all.csv'), index=False)
print('wrote', os.path.join(out, 'cv_results_all.csv'), len(d), 'rows')
