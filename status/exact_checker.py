#!/usr/bin/env python3
"""Independent exact rational checks of example graphs and formulas with SymPy.

No use of numpy or NetworkX eigenvalue routines.
"""
import sys
import sympy as sp

print('exact_checker v1')
print('python', sys.version.split()[0], 'sympy', sp.__version__)
a=sp.symbols('a', real=True)

def mat(n, edges, alpha):
    M=sp.zeros(n)
    for u,v in edges:
        M[u,v]=M[v,u]=1
    D=sp.diag(*[sum(M.row(i)) for i in range(n)])
    return alpha*D+(1-alpha)*M

# Q5: K2 disjoint K1, eigenvalues are 1, 2a-1, 0.
M=mat(3,[(0,1)],a)
assert sp.factor(M.charpoly().as_expr()) == sp.factor((M.charpoly().gen-1)*M.charpoly().gen*(M.charpoly().gen-(2*a-1)))
print('Q5 K2+K1 characteristic polynomial:', sp.factor(M.charpoly().as_expr()))

# Odd-order two cliques sharing one vertex. Every noncentral vertex has degree m,
# and the +/- vector supported on the two sides has eigenvalue m-1+a.
for m in (1, 2, 3):
    n=2*m+1
    left=list(range(1,m+1)); right=list(range(m+1,2*m+1))
    edges=[]
    for group in (left+[0],right+[0]):
        for i,u in enumerate(group):
            for v in group[i+1:]: edges.append((u,v))
    B=mat(n,edges,a)
    x=sp.Matrix([0]+[1]*m+[-1]*m)
    ev=m-1+a
    assert sp.simplify(B*x-ev*x)==sp.zeros(n,1)
    print(f'clique_share n={n}: exact eigenpair lambda={ev}; vertices={n}; edges={len(edges)}; VERIFIED')
    if n>3:
        for q in (sp.Rational(0),sp.Rational(1,4),sp.Rational(49,100),sp.Rational(1,2)):
            vals=sorted([float(sp.re(z.evalf(35))) for z,mult in B.subs(a,q).eigenvals().items() for _ in range(mult)])
            assert abs(vals[-2]-float(ev.subs(a,q)))<1e-11, (n,q,vals)
            print(f' n={n} alpha={q} second_eigenvalue={vals[-2]:.12f} formula={float(ev.subs(a,q)):.12f} OK')

# For the family K_{r,s}, its nontrivial 2x2 block eigenvalues:
# (a(r+s) +/- sqrt(a²(r-s)²+4(1-a)²rs))/2.
r,s=sp.symbols('r s', positive=True, integer=True)
B=sp.Matrix([[a*s,(1-a)*s],[(1-a)*r,a*r]])
t=sp.symbols('t')
assert sp.factor(B.charpoly(t).as_expr())==sp.factor(t*t-a*(r+s)*t+(2*a-1)*r*s)
print('Krs quotient charpoly:', sp.factor(B.charpoly(t).as_expr()))
print('END exact_checker')
