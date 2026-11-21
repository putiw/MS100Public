import sys,time,json,os
from cvlib import *
metric=sys.argv[1]; REPS=int(sys.argv[2]) if len(sys.argv)>2 else 10
d=load(metric); P=preds(metric); pt=d[~d.SubjectID.str.startswith('sub-C')].copy()
for c in('T25FW','x9HPTD','x9HPTND'): pt['log'+c]=np.log(pt[c].where(pt[c]>0))
MODELS=['tract16','classical4','tract4_tail','tract8_tail']
rows=[]
def boot_delta(y,a,b,binary,B=2000,seed=1):
    rng=np.random.default_rng(seed); n=len(y); out=[]
    pos=np.where(y==1)[0]; neg=np.where(y==0)[0]
    for _ in range(B):
        if binary: i=np.r_[rng.choice(pos,len(pos)),rng.choice(neg,len(neg))]; out.append(auc(y[i],a[i])-auc(y[i],b[i]))
        else:
            i=rng.choice(n,n); sst=np.sum((y[i]-y[i].mean())**2)
            out.append((np.sum((y[i]-b[i])**2)-np.sum((y[i]-a[i])**2))/sst)
    return np.percentile(out,[2.5,97.5])
OUTS=[('EDSS>=4','EDSS',4.0,True),('MSPro>=2','MSPro',2.0,True),('T25FW','logT25FW',None,False),('9HPT-D','logx9HPTD',None,False),('9HPT-ND','logx9HPTND',None,False),('SDMT','SDMTcorrect',None,False),('MSFC-SDMT','MSFC',None,False)]
for label,col,thr,binary in OUTS:
    cols=P['tract16']+P['classical4']+P['tract4_tail']
    m=pt[[col]+cols].notna().all(1)&(pt[cols]!=0).all(1)
    s=pt[m]; y=(s[col]>=thr).to_numpy().astype(float) if binary else s[col].to_numpy(float)
    oof={}
    for name in MODELS:
        X=s[P[name]].to_numpy(float); t=time.time()
        lp,lam=apparent(X,y,binary)
        o,base=nested(X,y,binary,reps=REPS)
        if binary:
            app=auc(y,lp); cv=np.mean([auc(y,o[r]) for r in range(REPS)]); sd=np.std([auc(y,o[r]) for r in range(REPS)])
        else:
            app=np.corrcoef(lp,y)[0,1]**2
            r2=[1-np.sum((y-o[r])**2)/np.sum((y-base[r])**2) for r in range(REPS)]; cv=np.mean(r2); sd=np.std(r2)
        oof[name]=o.mean(0)
        rows.append(dict(metric=metric,outcome=label,model=name,n_pred=X.shape[1],N=len(y),pos=int(y.sum()) if binary else '',measure='AUC' if binary else 'R2',apparent=round(app,3),cv=round(cv,3),cv_sd_across_repeats=round(sd,3),optimism=round(app-cv,3)))
        print(rows[-1],f'{time.time()-t:.0f}s',flush=True)
    for a,b in (('tract16','classical4'),('tract4_tail','classical4'),('tract8_tail','classical4')):
        lo,hi=boot_delta(y,oof[a],oof[b],binary)
        pe=(auc(y,oof[a])-auc(y,oof[b])) if binary else (np.sum((y-oof[b])**2)-np.sum((y-oof[a])**2))/np.sum((y-y.mean())**2)
        rows.append(dict(metric=metric,outcome=label,model=f'DELTA {a} - {b}',N=len(y),measure='cvAUC diff' if binary else 'cvR2 diff',cv=round(pe,3),ci_lo=round(lo,3),ci_hi=round(hi,3)))
        print(rows[-1],flush=True)
os.makedirs(os.environ.get('MS100_CV_OUT','.'),exist_ok=True)
pd.DataFrame(rows).to_csv(os.path.join(os.environ.get('MS100_CV_OUT','.'),f'cv_results_{metric}.csv'),index=False)
