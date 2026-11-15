"""Independent exact polynomial cross-check for the C++ congruence scan.

Reconstructs Atlas orders 5/7 using NetworkX and checks characteristic
polynomial coefficients over Z, NOT the C++ elimination code. The sign
variations of a real-rooted characteristic polynomial equal its count of
strictly positive roots including multiplicities. Zero multiplicity is determined by the trailing zero coefficients.
"""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import networkx as nx
import sympy as sp

ROOT = Path(__file__).resolve().parent
FRACS = [(0,1),(1,20),(1,10),(1,8),(1,6),(1,5),(1,4),(1,3),
         (2,5),(9,20),(12,25),(49,100),(499,1000)]

def graph_line(g):
    return nx.to_graph6_bytes(g,header=False).strip().decode('ascii')

def shifted(g,p,q):
    n = g.number_of_nodes()
    A = sp.zeros(n)
    for u,v in g.edges():
        A[u,v] = A[v,u] = q-p
    for u,d in g.degree():
        A[u,u] = p*(d-1)-q*((n-3)//2)
    return A

def inertia_by_polynomial(M):
    x=sp.Symbol('x')
    coeff=[int(v) for v in M.charpoly(x).all_coeffs()]
    nz=[1 if v>0 else -1 for v in coeff if v!=0]
    pos=sum(a!=b for a,b in zip(nz,nz[1:]))
    zeros=0
    for v in reversed(coeff):
        if v==0: zeros+=1
        else: break
    return pos,zeros,coeff

def check_one(g,p,q):
    pos,z,coeff=inertia_by_polynomial(shifted(g,p,q))
    assert pos<=1, (graph_line(g),p,q,pos)
    return (pos==1 and z>=1) or (pos==0 and z>=2)

def witness(n):
    m=(n-1)//2
    g=nx.disjoint_union(nx.complete_graph(m),nx.complete_graph(m))
    g.add_node(2*m)
    g.add_edges_from((2*m,k) for k in range(2*m))
    return g

def published_polynomial(n):
    m=(n-1)//2
    a,x=sp.symbols('a x')
    g=witness(n)
    M=sp.zeros(n)
    for u,v in g.edges(): M[u,v]=M[v,u]=1-a
    for u,d in g.degree(): M[u,u]=a*d
    P=sp.factor(M.charpoly(x).as_expr())
    Q=(x-(m+1)*a+1)**(2*m-2)*(x-m-a+1)*(
        x*x+(1-m-2*a*m-a)*x+2*a*m-2*m+2*a*m*m)
    assert sp.expand(P-Q)==0
    return P

def main():
    for n, g6 in ((5, 'D`{'), (7, 'FK|ko'), (9, 'HJ]CKN~')):
        g=witness(n)
        assert nx.is_isomorphic(g,nx.from_graph6_bytes(g6.encode('ascii')))
        published_polynomial(n)
        print('published polynomial verified n=',n,'g6=',graph_line(g))
    for n in (5,7):
        graphs=[g for g in nx.graph_atlas_g() if len(g)==n]
        for p,q in ((0,1),(1,4),(2,5),(49,100),(499,1000)):
            ties=sum(check_one(g,p,q) for g in graphs)
            print('independent polynomial n=%d alpha=%d/%d graphs=%d ties=%d' % (n,p,q,len(graphs),ties))
    g6=(ROOT/'graphs9.g6').read_bytes()
    lines=g6.splitlines()
    assert len(lines)==274668
    # The C generator emits canonical numeric integer order, which need not be
    # lexicographic graph6 order. Duplicate graph6 lines are still forbidden.
    assert len(lines)==len(set(lines))
    print('generated9 records=',len(lines),'sha256=',sha256(g6).hexdigest())
    indices=[0,1,13,101,999,4096,12345,42000,99999,150000,200000,274667]
    sampled=[nx.from_graph6_bytes(lines[i]) for i in indices]
    sampled.append(witness(9))
    for p,q in FRACS:
        ties=sum(check_one(g,p,q) for g in sampled)
        print('n9 polynomial sample alpha=%d/%d count=%d ties=%d' % (p,q,len(sampled),ties))
    print('all independent polynomial checks passed')

if __name__=='__main__':
    main()
