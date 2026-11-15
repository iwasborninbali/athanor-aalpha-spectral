#!/usr/bin/env python3
"""Independent exact witness checker; no NetworkX or NumPy.
Hand-decoded graph6, rational matrix-vector arithmetic, exhaustive path DFS.
"""
from fractions import Fraction as F


def decode_graph6(code):
    # Basic graph6 format for 0<=n<=62, adjacency order (0,1),(0,2),(1,2),...
    b = [ord(c)-63 for c in code]
    assert len(b)>=1 and 0<=b[0]<=62
    n = b[0]
    bits = [((v >> k)&1) for v in b[1:] for k in (5,4,3,2,1,0)]
    assert len(bits) >= n*(n-1)//2
    a = [[0]*n for _ in range(n)]
    t=0
    for j in range(1,n):
        for i in range(j):
            a[i][j] = a[j][i] = bits[t]
            t += 1
    return a


def spectral_action(a, alpha, x):
    return [alpha*sum(a[i])*x[i] + (1-alpha)*sum(a[i][j]*x[j] for j in range(len(a)))
            for i in range(len(a))]


def verifies(a, alpha, x, lam):
    return spectral_action(a,alpha,x) == [lam*v for v in x]


def has_path(a, k):
    n=len(a)
    def rec(v,used):
        if len(used)==k:
            return True
        return any(rec(w,used|{w}) for w in range(n) if a[v][w] and w not in used)
    return any(rec(v,{v}) for v in range(n))


def run():
    # Q5: K2+K1 has two invariant blocks, with eigenvalues 1, 2a-1, 0.
    a=decode_graph6('B_')
    assert a==[[0,1,0],[1,0,0],[0,0,0]]
    for alpha in [F(0),F(1,4),F(1,2),F(3,4),F(1)]:
        assert verifies(a,alpha,[F(1),F(1),F(0)],F(1))
        assert verifies(a,alpha,[F(1),F(-1),F(0)],2*alpha-1)
        assert verifies(a,alpha,[F(0),F(0),F(1)],F(0))
    print('PASS Q5 K2+K1 graph6=B_ exact eigenvalues 1, 2alpha-1, 0; lambda2 kink at 1/2')

    # P22: K4 minus edge, alpha=2/5, degree vector is strictly positive Perron vector.
    a=decode_graph6('C|')
    d=[sum(row) for row in a]
    assert sorted(d)==[2,2,3,3]
    assert verifies(a,F(2,5),d,F(13,5))
    row=[F(2,5)*d[i]+F(3,5)*F(sum(a[i][j]*d[j] for j in range(4)),d[i]) for i in range(4)]
    assert row==[F(13,5)]*4
    # Connected nonnegative irreducible matrix + positive eigenvector => Perron root.
    print('PASS P22 graph6=C| degrees='+str(d)+' alpha=2/5 lambda=13/5; eqs20,21 exact')

    # Infinite P22 split-graph family K_t join I_m, t,m>=2.
    for t,m in [(2,2),(2,3),(3,2),(4,5),(10,11)]:
        n=t+m
        b=[[int(i!=j and (i<t or j<t)) for j in range(n)] for i in range(n)]
        d=[sum(row) for row in b]
        alpha=F(m,2*m+t-1)
        lam=F(n-1)-alpha*(m-1)
        assert sorted(set(d))==[t,n-1]
        assert verifies(b,alpha,d,lam)
        print(f'PASS P22 FAMILY t={t} m={m} alpha={alpha} lambda={lam} nonregular')

    # Odd-n lambda2 lower-bound construction K1 join (two copies of Ks).
    for s in [1,2,3,4,5,10]:
        n=2*s+1
        b=[[0]*n for _ in range(n)]
        for i in range(1,n): b[0][i]=b[i][0]=1
        for lo in (1,s+1):
            for i in range(lo,lo+s):
                for j in range(i+1,lo+s): b[i][j]=b[j][i]=1
        x=[0]+[1]*s+[-1]*s
        for alpha in [F(0),F(1,4),F(2,5),F(49,100),F(3,4),F(9,10)]:
            assert verifies(b,alpha,x,F(s-1)+alpha)
        print('PASS ODD n='+str(n)+' alpha-rational exact eigenvalue (n-3)/2+alpha; lambda2 via block ordering')

    # P24 witness graph6 and noninduced path check: both avoid P5.
    for name, code in [('K4+K3','FKdbG'),('hub_plus_leaf_edge','FsaC_')]:
        g=decode_graph6(code)
        assert len(g)==7 and not has_path(g,5)
        print('PASS P24 P5-free n=7 '+name+' graph6='+code+' edges='+str(sum(map(sum,g))//2))
    print('CHECKER SUCCESS')

if __name__=='__main__':
    run()
