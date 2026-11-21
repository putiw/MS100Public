"""Nested cross-validation check of tract-based vs classical ridge models.
Reads participant-level tables from ../derivativesStats; writes aggregate results only."""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
import os
D=os.path.join(os.environ.get('MS100_STATS_DIR','../derivativesStats'),'')
CL=['Association','Cerebellar','Occipitoparietal','ProjectionBrainstem']
ICO=['periventricular','juxtacortical','infratentorial','deepwhitematter']
def load(metric):
    c=pd.read_excel(D+'clinicalScore.xlsx')
    if metric in ('NDI','ODI','FWF'):
        t=pd.read_csv(D+f'GroupTract{metric}_All.csv').drop_duplicates('SubjectID')
        i=pd.read_csv(D+'ClassicalRegionNODDI_All.csv')
    else:
        t=pd.read_excel(D+f'GroupTract{metric}_All.xlsx').drop_duplicates('SubjectID')
        i=pd.read_excel(D+'icometrixLesionStats.xlsx')
    d=c.merge(t,on='SubjectID',how='left').merge(i,on='SubjectID',how='left')
    # MSFC-SDMT (z against MS cohort)
    ms=~d.SubjectID.str.startswith('sub-C')
    arm=(1/d.x9HPTD+1/d.x9HPTND)/2; leg=-d.T25FW; cog=d.SDMTcorrect.astype(float)
    z=lambda v:(v-v[ms].mean())/v[ms].std()
    d['MSFC']=(z(arm)+z(leg)+z(cog))/3
    return d
def preds(metric):
    P16=[f'{g}{h}_{k}' for g in CL for k in('Tail','NAWM') for h in('L','R')]
    P4c=[f'{r}{metric}' for r in ICO]
    P4t=[f'{g}_Tail' for g in CL]
    P4l=[f'{g}_Lesion' for g in CL]
    P8=[f'{g}{h}_Tail' for g in CL for h in('L','R')]
    return dict(tract16=P16,classical4=P4c,tract4_tail=P4t,tract4_lesion=P4l,tract8_tail=P8)
def auc(y,s):
    y=np.asarray(y).astype(bool); r=pd.Series(s).rank().to_numpy()
    n1=y.sum(); n0=(~y).sum()
    return (r[y].sum()-n1*(n1+1)/2)/(n1*n0)
def zs(Xtr,Xte):
    m=Xtr.mean(0); s=Xtr.std(0,ddof=1); s[s==0]=1
    return (Xtr-m)/s,(Xte-m)/s
def fit_logit(X,y,lam,it=60):
    n,p=X.shape; b=np.zeros(p); b0=np.log(y.mean()/(1-y.mean()))
    for _ in range(it):
        eta=np.clip(b0+X@b,-30,30); mu=1/(1+np.exp(-eta)); w=np.maximum(mu*(1-mu),1e-9)
        g=np.r_[np.mean(mu-y), X.T@(mu-y)/n+lam*b]
        Xa=np.c_[np.ones(n),X]; H=(Xa.T*w)@Xa/n; H[1:,1:]+=lam*np.eye(p); H[0,0]+=1e-10
        step=np.linalg.solve(H,g)
        b0-=step[0]; b-=step[1:]
        if np.abs(step).max()<1e-8: break
    return b0,b
def fit_lin(X,y,lam):
    n,p=X.shape; ym=y.mean()
    b=np.linalg.solve(X.T@X/n+lam*np.eye(p),X.T@(y-ym)/n)   # X centred
    return ym,b
def folds(y,k,rng,strat):
    idx=np.arange(len(y)); f=np.empty(len(y),int)
    if strat:
        for cls in (0,1):
            ii=rng.permutation(idx[y==cls]); f[ii]=np.arange(len(ii))%k
    else:
        ii=rng.permutation(idx); f[ii]=np.arange(len(ii))%k
    return f
GRID=np.logspace(-4,4,25)
def inner_lambda(X,y,binary,rng,k=5):
    f=folds(y,k,rng,binary); loss=np.zeros(len(GRID))
    for j in range(k):
        tr=f!=j; te=~tr
        Xtr,Xte=zs(X[tr],X[te])
        for li,lam in enumerate(GRID):
            if binary:
                b0,b=fit_logit(Xtr,y[tr],lam); eta=np.clip(b0+Xte@b,-30,30)
                loss[li]+=np.sum(np.log1p(np.exp(eta))-y[te]*eta)      # deviance
            else:
                b0,b=fit_lin(Xtr,y[tr],lam); loss[li]+=np.sum((y[te]-b0-Xte@b)**2)
    return GRID[np.argmin(loss)]
def nested(X,y,binary,reps=10,k=10,seed=0):
    """returns OOF predictions [reps x n]"""
    rng=np.random.default_rng(seed); out=np.zeros((reps,len(y))); base=np.zeros((reps,len(y)))
    for r in range(reps):
        f=folds(y,k,rng,binary)
        for j in range(k):
            tr=f!=j; te=~tr
            lam=inner_lambda(X[tr],y[tr],binary,rng)
            Xtr,Xte=zs(X[tr],X[te])
            if binary: b0,b=fit_logit(Xtr,y[tr],lam)
            else: b0,b=fit_lin(Xtr,y[tr],lam)
            out[r,te]=b0+Xte@b; base[r,te]=y[tr].mean()
    return out,base
def apparent(X,y,binary,lam=None,seed=0):
    rng=np.random.default_rng(seed)
    if lam is None: lam=inner_lambda(X,y,binary,rng,k=10)
    Xz,_=zs(X,X)
    b0,b=(fit_logit if binary else fit_lin)(Xz,y,lam)
    return b0+Xz@b,lam
