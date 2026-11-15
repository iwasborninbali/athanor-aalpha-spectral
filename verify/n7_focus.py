#!/usr/bin/env python3
"""Exhaustive n=7, alpha=1/4: bitmap maxcut and NumPy eigenvalues."""
import time
import numpy as np
n=7; alpha=.25; pairs=[(u,v) for u in range(n) for v in range(u+1,n)]
cuts=[sum(1<<j for j,(u,v) in enumerate(pairs) if ((s>>u)^(s>>v))&1) for s in range(1<<(n-1))]
t0=time.process_time(); wrong=0; corrected=0; p20=0; eq=0; first_wrong=None; first_eq=None
for mask in range(1<<21):
    A=np.zeros((n,n),dtype=float)
    for j,(u,v) in enumerate(pairs):
        if mask&(1<<j):A[u,v]=A[v,u]=1
    d=A.sum(1);m=mask.bit_count()
    ev=np.linalg.eigvalsh(alpha*np.diag(d)+(1-alpha)*A);small=ev[0];large=ev[-1]
    mc=max((mask & cut).bit_count() for cut in cuts)
    paper=2*alpha*m/n-2*(1-alpha)*mc/n
    fixed=2*m/n-4*(1-alpha)*mc/n
    if small>paper+1e-8:
        wrong+=1
        if first_wrong is None:first_wrong=(mask,round(small,9),round(paper,9),mc)
    if small>fixed+1e-8:corrected+=1
    if min(d)>0:
        rs=alpha*d+(1-alpha)*(A@d)/d
        if large>max(rs)+1e-8 or large<min(rs)-1e-8:p20+=1
        if max(rs)-min(rs)<1e-8 and max(d)!=min(d):
            eq+=1
            if first_eq is None:first_eq=(mask,list(d.astype(int)))
    if mask and mask%262144==0:
        print('progress',mask,'cpu_s',round(time.process_time()-t0,2),flush=True)
print('FINAL n',n,'graphs',1<<21,'alpha',alpha,'P29_original_violations',wrong,
'first',first_wrong,'P29_corrected_violations',corrected,'P20_bound_violations',p20,
'P20_nonregular_equal_row_sums',eq,'first',first_eq,'cpu_s',round(time.process_time()-t0,3),flush=True)
