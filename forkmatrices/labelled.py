"""Independent stdlib-only labeled graph enumeration, n=2..6.
No NetworkX, no eigenvalues, no atlas and no floating point.
"""
from fractions import Fraction
from collections import Counter


def connected(a):
    n = len(a)
    seen = {0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j in range(n):
            if a[i][j] and j not in seen:
                seen.add(j)
                stack.append(j)
    return len(seen) == n


def rows(a):
    n = len(a)
    dist = [[0 if i == j else (1 if a[i][j] else n+1) for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                dist[i][j] = min(dist[i][j], dist[i][k]+dist[k][j])
    d = [sum(r) for r in a]
    q = [Fraction(d[i]*d[i]+sum(a[i][j]*d[j] for j in range(n)),d[i]) for i in range(n)]
    t = [sum(r) for r in dist]
    second = [sum(dist[i][j]*t[j] for j in range(n)) for i in range(n)]
    dr = [Fraction(second[i],t[i]) for i in range(n)]
    return d,q,t,dr


def alpha(t,dr):
    cand = None
    for i in range(1,len(t)):
        delta=dr[i]-dr[0]
        den=Fraction(t[0]-t[i])+delta
        if not den:
            if delta: return None
            continue
        b=delta/den
        if not 0<=b<1:return None
        if cand is None:cand=b
        elif cand != b:return None
    if cand is None:return None
    return cand if len({cand*t[i]+(1-cand)*dr[i] for i in range(len(t))})==1 else None

for n in range(2,7):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    c=Counter()
    for mask in range(1<<len(edges)):
        a=[[0]*n for _ in range(n)]
        for bit,(i,j) in enumerate(edges):
            if mask>>bit&1:a[i][j]=a[j][i]=1
        if not connected(a):continue
        c['connected']+=1
        d,q,t,dr=rows(a)
        if len(set(d))!=1:
            c['nonregular']+=1
            if len(set(q))==1:c['Q']+=1
        if len(set(t))==1:continue
        c['nontransregular']+=1
        if len(set(dr))==1:c['Dist']+=1
        if len(set(Fraction(t[i])+dr[i] for i in range(n)))==1:c['DistQ']+=1
        b=alpha(t,dr)
        if b is not None:
            c['Dalpha']+=1
            if b==0:c['alpha0']+=1
            if b>=Fraction(1,2):c['alpha_ge_half']+=1
    print('n=%d connected=%d nonregular=%d nontransregular=%d Q=%d Dist=%d DistQ=%d Dalpha=%d alpha0=%d alpha_ge_half=%d' % (n,c['connected'],c['nonregular'],c['nontransregular'],c['Q'],c['Dist'],c['DistQ'],c['Dalpha'],c['alpha0'],c['alpha_ge_half']))
print('method=all labeled undirected edge masks n=2..6; stdlib Fractions + Floyd-Warshall')
