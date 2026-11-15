#!/usr/bin/env python3
"""Independent exact checker for the n=3,5,7 continuous interval result.

Unlike interval_scan.py, constructs the symbolic determinant directly
using SymPy's polynomial-domain Gaussian elimination, and counts the
positive inertia by rational symmetric Schur congruences, NOT Descartes.
Reuses NetworkX's atlas and SymPy root isolation: not an independent
unlabeled-graph generator or independent root-isolation library.
"""
import sys
from fractions import Fraction
import networkx as nx
import sympy as sp

z = sp.Symbol('z')


def matrix_pencil(G):
    n = len(G)
    m = (n - 1) // 2
    return sp.Matrix(n, n, lambda i, j:
        z * (G.degree(i) - 1) - (m - 1) if i == j
        else (1 - z if G.has_edge(i, j) else 0))


def as_fraction(value):
    v = sp.Rational(value)
    return Fraction(int(v.p), int(v.q))


def congruence_inertia(mat):
    """Exact Sylvester inertia by rational, symmetric congruences.

    1x1 pivots are taken after diagonal permutation. An invertible
    basis shear produces a diagonal pivot from nonzero offdiagonal if
    necessary. Each pivot contributes exactly its sign.
    """
    n = mat.rows
    B = [[as_fraction(mat[i, j]) for j in range(n)] for i in range(n)]
    positive = zero = negative = 0
    while B:
        n = len(B)
        pivot = next((i for i in range(n) if B[i][i]), None)
        if pivot is None:
            pair = next(((i, j) for i in range(n) for j in range(i + 1, n)
                         if B[i][j]), None)
            if pair is None:
                zero += n
                break
            i, j = pair
            # e_i -> e_i + e_j, e_k -> e_k for k != i
            old = [row[:] for row in B]
            for k in range(n):
                B[i][k] = old[i][k] + old[j][k]
                B[k][i] = B[i][k]
            B[i][i] = old[i][i] + 2 * old[i][j] + old[j][j]
            pivot = i
        if pivot:
            B[0], B[pivot] = B[pivot], B[0]
            for row in B:
                row[0], row[pivot] = row[pivot], row[0]
        d = B[0][0]
        assert d
        if d > 0:
            positive += 1
        else:
            negative += 1
        B = [[B[i][j] - B[i][0]*B[0][j]/d for j in range(1, n)]
             for i in range(1, n)]
    return positive, zero, negative


def isolate_roots_direct(poly):
    q = poly
    for factor, endpoint in ((sp.Poly(z,z), sp.Integer(0)),
                             (sp.Poly(2*z-1,z), sp.Rational(1,2))):
        while q.degree() > 0 and q.eval(endpoint) == 0:
            q = q.exquo(factor)
    intervals = q.intervals(eps=sp.Rational(1,2**32))
    return [(lo, hi, mult) for (lo, hi), mult in intervals
            if lo > 0 and hi < sp.Rational(1,2)]


def sample_points(roots):
    samples = []
    previous = sp.Integer(0)
    for lo, hi, mult in roots:
        assert mult == 1
        samples += [(previous + lo)/2]
        if lo != hi:
            samples += [lo, hi]
        previous = hi
    samples += [(previous + sp.Rational(1,2))/2]
    return sorted(set(samples))


def main():
    print('method=direct SymPy polynomial determinant + isolated real roots + rational symmetric congruences')
    print('python=' + sys.version.split()[0], 'networkx=' + nx.__version__, 'sympy=' + sp.__version__)
    atlas = nx.graph_atlas_g()
    for n in (3, 5, 7):
        graph_count = permanent = interior_graph_count = root_count = trials = positive_max = equal_zero = 0
        for G in atlas:
            if len(G) != n:
                continue
            graph_count += 1
            code = nx.to_graph6_bytes(G, header=False).decode().strip()
            M = matrix_pencil(G)
            # This polynomial is constructed without integer point evaluation.
            p = sp.Poly(M.det(method='domain-ge'), z)
            pos, zeros, neg = congruence_inertia(M.subs(z, 0))
            assert pos <= 1, (n, code, pos)
            if (pos == 1 and zeros > 0) or (pos == 0 and zeros >= 2):
                equal_zero += 1
            if p.is_zero:
                permanent += 1
                h = nx.disjoint_union(nx.complete_graph((n-1)//2), nx.complete_graph((n-1)//2))
                h.add_node(n-1)
                h.add_edges_from((n-1, k) for k in range(n-1))
                assert nx.is_isomorphic(G, h), (n, code)
                continue
            roots = isolate_roots_direct(p)
            assert all(mult == 1 for _, _, mult in roots), (n, code, roots)
            root_count += len(roots)
            interior_graph_count += bool(roots)
            # Independent root count (Sturm), including endpoint removal.
            number = p.count_roots(0, sp.Rational(1, 2))
            number -= (p.eval(0)==0) + (p.eval(sp.Rational(1,2))==0)
            assert number == len(roots), (n, code, number, roots)
            for alpha in sample_points(roots):
                pos, zeros, neg = congruence_inertia(M.subs(z,alpha))
                assert zeros == 0 and pos <= 1 and pos+neg == n, (n, code, alpha, pos, zeros, neg)
                positive_max = max(positive_max, pos)
                trials += 1
        assert (graph_count,permanent,equal_zero)=={3:(4,1,3),5:(34,1,6),7:(1044,1,11)}[n]
        print('n=%d graphs=%d alpha0_equal=%d identically_zero=%d root_graphs=%d interior_simple_roots=%d rational_pivot_tests=%d max_positive=%d failures=0' %
              (n, graph_count, equal_zero, permanent, interior_graph_count, root_count, trials, positive_max))
    print('INDEPENDENT CHECK PASSED for every odd-order unlabeled graph n=3,5,7 and all alpha in [0,1/2).')

if __name__=='__main__':
    main()
