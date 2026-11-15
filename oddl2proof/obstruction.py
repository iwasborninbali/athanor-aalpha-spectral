#!/usr/bin/env python3
"""Exact five-vertex counterexample to an invalid lambda_2 interpolation rule."""
import networkx as nx
import sympy as sp

G = nx.from_graph6_bytes(b'Dlc')
assert list(G.edges()) == [(0,1),(0,3),(0,4),(1,2),(2,3),(3,4)]
bounds = [(sp.Rational(0),sp.Rational(68890,100000),'upper'),
          (sp.Rational(1,4),sp.Rational(103214,100000),'lower'),
          (sp.Rational(1,2),sp.Rational(137296,100000),'upper')]
for alpha,bound,way in bounds:
    M = sp.Matrix(5,5,lambda i,j: (alpha*G.degree(i) if i==j else
                                  (1-alpha if G.has_edge(i,j) else 0)))
    p = M.charpoly().as_poly()
    roots = p.intervals(eps=sp.Rational(1,10**9))
    (lo,hi),mult = roots[-2]
    assert mult==1
    assert (hi<bound if way=='upper' else lo>bound)
    print('alpha=%s poly=%s lambda2_isolated_in=[%s,%s] %s_bound=%s PASS' %
          (alpha,sp.factor(p.as_expr()),lo,hi,way,bound))
assert bounds[1][1] > (bounds[0][1]+bounds[2][1])/2
print('STRICT FAILURE: lambda2(A_1/4) > 103214/100000 > (68890/100000+137296/100000)/2 > (lambda2(A_0)+lambda2(A_1/2))/2')
