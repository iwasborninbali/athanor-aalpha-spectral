#!/usr/bin/env python3
"""Independent small-graph verification of arXiv:1607.03015v1.
Only Python stdlib and NumPy are used; graph generation, coloring, maxcut
and extremal candidate construction are implemented here from scratch.
"""
import argparse, collections, itertools, math, sys, time
import numpy as np

ALPHAS = (0.0, 0.25, 0.5, 0.75, 1.0)
EPS = 2e-8

def graph(n, mask):
    a = np.zeros((n, n), dtype=float)
    edges = []
    j = 0
    for u in range(n):
        for v in range(u+1, n):
            if (mask >> j) & 1:
                a[u,v] = a[v,u] = 1.0
                edges.append((u,v))
            j += 1
    return a, edges

def eigs(a, alpha):
    d = a.sum(axis=1)
    return np.linalg.eigvalsh(alpha*np.diag(d)+(1-alpha)*a)[::-1]

def connected(a):
    n=len(a); seen={0}; todo=[0]
    while todo:
        v=todo.pop()
        for w in np.flatnonzero(a[v]):
            w=int(w)
            if w not in seen:
                seen.add(w); todo.append(w)
    return len(seen)==n

def chi(a):
    n=len(a); col=[-1]*n
    order=sorted(range(n),key=lambda v:-sum(a[v]))
    def rec(i,k):
        if i==n: return True
        v=order[i]
        for c in range(k):
            if all(col[w]!=c or not a[v,w] for w in range(n)):
                col[v]=c
                if rec(i+1,k): return True
                col[v]=-1
        return False
    for k in range(1,n+1):
        if rec(0,k): return k
    raise AssertionError()

def omega(a):
    n=len(a); best=1
    for s in range(1,1<<n):
        k=s.bit_count()
        if k<=best: continue
        vv=[u for u in range(n) if (s>>u)&1]
        if all(a[u,v] for u,v in itertools.combinations(vv,2)):
            best=k
    return best

def maxcut(edges,n):
    return max((sum(((s>>u)^(s>>v))&1 for u,v in edges)
                for s in range(1<<(n-1))),default=0)

def multipartite(parts):
    n=sum(parts); a=np.ones((n,n))-np.eye(n); offset=0
    for p in parts:
        a[offset:offset+p,offset:offset+p]=0
        offset+=p
    return a

def turan(n,r):
    return multipartite([n//r+(i<n%r) for i in range(r)])

def split(n,r):
    return multipartite([1]*(r-1)+[n-r+1])

def is_regular(a):
    d=a.sum(1); return bool(np.all(d==d[0]))

def note_bad(fails, name, n, mask, alpha, observed, claimed):
    if name not in fails:
        fails[name] = (n,mask,alpha,round(float(observed),11),round(float(claimed),11))

def scan(nmax):
    start=time.process_time(); counts=collections.Counter(); fails={}
    minima={}; equal=collections.Counter(); witnesses={}
    for n in range(1,nmax+1):
        N=1<<(n*(n-1)//2)
        for mask in range(N):
            a,edges=graph(n,mask)
            d=a.sum(1); m=len(edges); lo=min(d); hi=max(d); reg=is_regular(a)
            noniso=(lo>0); conn=connected(a)
            omega_val=omega(a) if n>=2 else 1
            chrom=chi(a) if n>=2 else 1
            Aeval=eigs(a,0); w=a@d
            max_cut=maxcut(edges,n) if n<=6 else None
            counts['graphs']+=1; counts[f'graphs_n{n}']+=1
            alpha_eigs=[]
            for alpha in ALPHAS:
                eig=eigs(a,alpha); alpha_eigs.append(eig)
                rho=eig[0]; smallest=eig[-1]; N1=n
                def le(name,x,y):
                    counts[name]+=1
                    if x>y+EPS: note_bad(fails,name,n,mask,alpha,x,y)
                def ge(name,x,y): le(name,y,x)
                le('P9_eig_degree',max(eig-d[::-1]),0) if False else None
                for idx, val in enumerate(eig):
                    le('P9_degree_bound',val,sorted(d,reverse=True)[idx])
                    le('P10_upper',val,alpha*hi+(1-alpha)*Aeval[idx])
                    ge('P10_lower',val,alpha*lo+(1-alpha)*Aeval[idx])
                le('P9_radius',rho,hi)
                if alpha>0.5:
                    ge('P6_PSD',smallest,0)
                    if noniso: ge('P6_PD',smallest,1e-12)
                root=.5*(alpha*(hi+1)+math.sqrt(alpha**2*(hi+1)**2+4*hi*(1-2*alpha)))
                ge('P11_bound',rho,root)
                if alpha<=0.5: ge('C12_low',rho,alpha*(hi+1))
                if 0.5<=alpha<1: ge('C12_high',rho,alpha*hi+1-alpha)
                ge('P18_lower',rho,Aeval[0]); le('P18_upper',rho,alpha*hi+(1-alpha)*Aeval[0])
                ge('C19_RMS',rho,math.sqrt(float(d@d)/n))
                ge('C19_mean',rho,2*m/n)
                if noniso:
                    rows=alpha*d+(1-alpha)*w/d
                    le('P20_upper',rho,max(rows)); ge('P20_lower',rho,min(rows))
                    if conn and not reg:
                        if abs(rho-max(rows))<1e-8: equal[f'P20upper@{alpha}']+=1; witnesses.setdefault(f'P20upper@{alpha}',(n,mask))
                        if abs(rho-min(rows))<1e-8: equal[f'P20lower@{alpha}']+=1; witnesses.setdefault(f'P20lower@{alpha}',(n,mask))
                if edges:
                    vals=[z for u,v in edges for z in (alpha*d[u]+(1-alpha)*d[v],alpha*d[v]+(1-alpha)*d[u])]
                    le('C21_upper',rho,max(vals)); ge('C21_lower',rho,min(vals))
                walks=alpha*d**2+(1-alpha)*w
                le('P23_upper',rho*rho,max(walks)); ge('P23_lower',rho*rho,min(walks))
                le('P31a' if alpha>=.5 else 'P31b',eig[1],alpha*n-1 if alpha>=.5 else n/2-1) if n>=2 else None
                ge('P29_original',2*alpha*m/n-2*(1-alpha)*max_cut/n,smallest) if max_cut is not None else None
                # Correct Rayleigh bound independent of spectrum calculation.
                le('P29_corrected',smallest,2*m/n-4*(1-alpha)*max_cut/n) if max_cut is not None else None
                le('P34_trace',abs(sum(eig)-2*alpha*m),1e-8)
                le('P35_trace_sq',abs(eig@eig-(2*(1-alpha)**2*m+alpha**2*(d@d))),1e-8)
                if n>=2:
                    if chrom<=2:
                        t=eigs(turan(n,2),alpha)[0]; s=eigs(split(n,2),alpha)[0]
                        if alpha<.5: le('T25_i',rho,t)
                        if alpha>.5: le('T25_ii',rho,s)
                        if alpha==.5: le('T25_iii',rho,n/2)
                    for r in range(2,min(4,n)+1):
                        t=eigs(turan(n,r),alpha)[0]
                        s=eigs(split(n,r),alpha)[0]
                        cutoff=1-1/r
                        if chrom<=r:
                            le('T26_i' if alpha<cutoff else 'T26_ii' if alpha>cutoff else 'T26_iii',rho,t if alpha<cutoff else s if alpha>cutoff else cutoff*n)
                        if omega_val<=r and alpha<1:
                            le('T27_i' if alpha<cutoff else 'T27_ii' if alpha>cutoff else 'T27_iii',rho,t if alpha<cutoff else s if alpha>cutoff else cutoff*n)
                    if 0<alpha<.5:
                        key=(n,alpha)
                        if key not in minima or smallest<minima[key][0]-1e-9:
                            minima[key]=(float(smallest),mask)
                counts['graph_alpha_pairs']+=1
                if alpha==1 and n>=3 and conn and not reg and sum(d==hi)>=2 and abs(eig[1]-(n-1))<EPS:
                    equal['P31a_noncomplete_endpoint']+=1
            for i in range(1,len(ALPHAS)):
                diff=alpha_eigs[i]-alpha_eigs[i-1]
                ge('P3_monotone',min(diff),0)
                le('P4_Lipschitz',max(diff),(ALPHAS[i]-ALPHAS[i-1])*n)
            if noniso and n>=2 and chrom>=2:
                # P8 (alpha < 1/chromatic) checked only on regular graphs.
                if reg:
                    for alpha,eig in zip(ALPHAS,alpha_eigs):
                        if alpha<1/chrom-EPS: le('P8_neg_min',eig[-1],-1e-10)
            if n<=5 and conn:
                # Corollary 33 and Proposition 32: graph-distance matrix powers.
                distances=np.full((n,n),n+1,dtype=int)
                for u in range(n):
                    distances[u,u]=0; queue=[u]
                    for v in queue:
                        for wv in np.flatnonzero(a[v]):
                            v2=int(wv)
                            if distances[u,v2]>n:
                                distances[u,v2]=distances[u,v]+1; queue.append(v2)
                diam=int(distances.max())
                for alpha in (0,.5,.75):
                    b=alpha*np.diag(d)+(1-alpha)*a
                    cur=np.eye(n)
                    for ell in range(1,diam+1):
                        cur=cur@b
                        for u in range(n):
                            for v in range(u+1,n):
                                if distances[u,v]==ell: ge('P32_distance_equal',cur[u,v],1e-12)
                                if distances[u,v]>ell: le('P32_distance_short',abs(cur[u,v]),1e-12)
                    nu=len(np.unique(np.round(eigs(a,alpha),8)))
                    ge('C33_distinct',nu,diam+1)
        print(f'n={n} labeled_graphs={N} cumulative={counts["graphs"]} cpu_s={time.process_time()-start:.2f}',flush=True)
    print('ALPHAS',ALPHAS, 'nmax',nmax)
    for k,v in sorted(counts.items()): print('COUNT',k,v)
    for k,v in sorted(fails.items()): print('FAIL',k,'n/mask/alpha/observed/claimed',v)
    for k,v in sorted(equal.items()): print('EQUALITY',k,v,'first_n_mask',witnesses.get(k))
    for k,v in sorted(minima.items()): print('MINIMUM_P30',k,'lambda_mask',tuple(round(x,9) if isinstance(x,float) else x for x in v))
    print('CPU_SECONDS',round(time.process_time()-start,3))

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--nmax',type=int,default=6)
    args=p.parse_args(); scan(args.nmax)