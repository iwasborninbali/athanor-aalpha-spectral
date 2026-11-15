"""Standalone exact witness checker; stdlib only, no spectral packages.
Both matrices and graph distances are reconstructed independently.
"""
from fractions import Fraction as F


def verify(edges,n,alpha,kind):
    adj=[[0]*n for _ in range(n)]
    for i,j in edges:
        adj[i][j]=1
        if kind!='directed':adj[j][i]=1
    if kind=='directed':
        # Strong connectivity: BFS on both digraph and transpose.
        for trans in (False,True):
            seen={0};todo=[0]
            while todo:
                i=todo.pop()
                for j in range(n):
                    if (adj[j][i] if trans else adj[i][j]) and j not in seen:
                        seen.add(j);todo.append(j)
            assert len(seen)==n
        v=[sum(row) for row in adj]
        m=adj
    else:
        dist=[[0 if i==j else (1 if adj[i][j] else n+1) for j in range(n)] for i in range(n)]
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    dist[i][j]=min(dist[i][j],dist[i][k]+dist[k][j])
        assert max(max(row) for row in dist)<n+1
        v=[sum(row) for row in dist]
        m=dist
    assert min(v)>0 and len(set(v))>1, 'nonregularity/transmission-irregularity check'
    a=F(alpha)
    mv=[a*v[i]*v[i]+(1-a)*sum(m[i][j]*v[j] for j in range(n)) for i in range(n)]
    roots=[mv[i]/v[i] for i in range(n)]
    assert len(set(roots))==1,(roots,v)
    # Positive vector and nonnegative irreducible matrix: Perron-Frobenius.
    return roots[0],tuple(v)


def split(t,m):
    n=t+m
    return [(i,j) for i in range(n) for j in range(i+1,n) if i<t or j<t]


def split_alpha(t,m):
    return F(2*m+2*t-2-m*t,2*m*m+2*m*t-2*m+t*t-t)

print('EXACT REPLAY, no numpy/networkx/sympy')
for t,m in [(1,2),(1,5),(2,2),(2,9),(3,2),(3,3),(3,4),(4,2),(4,3),(7,2)]:
    alpha=split_alpha(t,m)
    assert 0<=alpha<F(1,2)
    eig,v=verify(split(t,m),t+m,alpha,'distance')
    print(f'Dalpha split t={t} m={m} alpha={alpha} lambda={eig} transmissions={sorted(set(v))}')

for t,m in [(3,4),(4,3)]:
    eig,v=verify(split(t,m),t+m,F(0),'distance')
    print(f'PURE DISTANCE t={t} m={m}: rho={eig} transmissions={sorted(set(v))}')

# Directed split construction with one low-to-high outgoing arc per low.
def digraph_split(t,m):
    n=t+m
    return ([(i,j) for i in range(t) for j in range(n) if i!=j]
            +[(j,0) for j in range(t,n)])

for t,m in [(2,1),(2,4),(3,2),(5,7)]:
    a=F(m,t+2*m-1)
    eig,v=verify(digraph_split(t,m),t+m,a,'directed')
    print(f'DIRECTED Aalpha t={t} m={m} alpha={a} lambda={eig} outdegrees={sorted(set(v))}')

# Known Q semiregular-bipartite case K_{a,b}, exact degree eigenvector.
for a,b in [(1,2),(2,3),(3,4)]:
    n=a+b
    adj=[[int((i<a)!=(j<a)) for j in range(n)] for i in range(n)]
    d=[sum(r) for r in adj]
    roots=[F(d[i]*d[i]+sum(adj[i][j]*d[j] for j in range(n)),d[i]) for i in range(n)]
    assert len(set(roots))==1 and roots[0]==a+b
    print(f'Q K_{{{a},{b}}}: rho={roots[0]} degrees={sorted(set(d))}')
# Nonuniform positive rational edge weights: weighted semiregular K_2,3.
# This uses strengths (weighted degrees), not ordinary degrees.
wrows=[[F(1,2),F(1),F(3,2)],[F(3,2),F(1),F(1,2)]]
w=[[F(0) for j in range(5)] for i in range(5)]
for i in range(2):
    for j in range(3):
        w[i][j+2]=w[j+2][i]=wrows[i][j]
strength=[sum(row) for row in w]
assert strength==[3,3,2,2,2]
qroots=[(strength[i]*strength[i]+sum(w[i][j]*strength[j] for j in range(5)))/strength[i] for i in range(5)]
assert set(qroots)=={F(5)}
print(f'WEIGHTED A_1/2 nonuniform K_2,3: rho={qroots[0]/2} strengths={strength} weights={wrows}')
print('ALL CHECKS PASSED (exact rational eigenvector certificates)')
