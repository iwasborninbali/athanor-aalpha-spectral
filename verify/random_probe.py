#!/usr/bin/env python3
"""Independent seeded large-order stress tests of spectral inequalities.
Constructs adjacency matrices and degree data directly; no imports from scan.
"""
import math, random, time
import numpy as np
SEED=20261115; rng=random.Random(SEED)
NS=(7,8,10,12,16,20); TRIALS=1000
ALPHAS=(0,.1,.25,.4,.5,.6,.75,.9,1)
t0=time.process_time(); count=0; tests=0

for n in NS:
    local=0
    for t in range(TRIALS):
        p=(.05+.9*rng.random())
        A=np.zeros((n,n))
        for i in range(n):
            for j in range(i+1,n):
                if rng.random()<p: A[i,j]=A[j,i]=1
        d=A.sum(1); m=d.sum()/2; delta=min(d); Delta=max(d)
        w=A@d; adeig=np.linalg.eigvalsh(A)
        for alpha in ALPHAS:
            M=alpha*np.diag(d)+(1-alpha)*A
            v=np.linalg.eigvalsh(M); lo=v[0]; hi=v[-1]
            assert np.all(v>=alpha*delta+(1-alpha)*adeig-1e-8)
            assert np.all(v<=alpha*Delta+(1-alpha)*adeig+1e-8)
            assert hi+1e-8>=math.sqrt(float(d@d)/n)
            assert hi+1e-8>=2*m/n
            assert hi+1e-8>=.5*(alpha*(Delta+1)+math.sqrt(alpha*alpha*(Delta+1)**2+4*Delta*(1-2*alpha))) if Delta>0 else True
            assert hi**2+1e-8>=min(alpha*d*d+(1-alpha)*w)
            assert hi**2<=max(alpha*d*d+(1-alpha)*w)+1e-8
            if delta>0:
                row=alpha*d+(1-alpha)*w/d
                assert hi>=min(row)-1e-8 and hi<=max(row)+1e-8
            if alpha>.5 and delta>0: assert lo>0
            if alpha>=.5 and n>1: assert v[-2]<=alpha*n-1+1e-8
            if alpha<.5 and n>1: assert v[-2]<=n/2-1+1e-8
            assert abs(sum(v)-2*alpha*m)<1e-8
            assert abs(v@v-(2*(1-alpha)**2*m+alpha**2*(d@d)))<1e-8
            tests+=1
        local+=1;count+=1
    print('order',n,'random_graphs',local,'cumulative_graphs',count,'cpu_s',round(time.process_time()-t0,2),flush=True)
print('PASS random graphs',count,'graph_alpha_pairs',tests,'seed',SEED,'n_orders',NS,'alphas',ALPHAS,'cpu_s',round(time.process_time()-t0,3))
