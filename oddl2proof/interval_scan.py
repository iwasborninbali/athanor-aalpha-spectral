#!/usr/bin/env python3
"""Exact continuous-alpha proof for odd orders 3,5,7.

Reconstructs det(A_alpha -(m-1+alpha)I) by integer Bareiss
interpolation; isolates every root via exact SymPy polynomial arithmetic;
counts positive shifted eigenvalues by Descartes on a real-rooted
characteristic polynomial in every root-free interval.

Requires free Python 3, NetworkX, SymPy. No numpy, no floating point.
"""
import sys
import networkx as nx
import sympy as sp

ALPHA = sp.Symbol('a')


def bareiss_determinant(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    if not n:
        return 1
    sign = previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            index = next((i for i in range(k + 1, n) if a[i][k]), None)
            if index is None:
                return 0
            a[k], a[index] = a[index], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = a[i][j] * pivot - a[i][k] * a[k][j]
                assert value % previous == 0
                a[i][j] = value // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def shifted_integer_rows(G, x):
    """Integer entries when alpha is integer x (interpolation nodes)."""
    n = len(G)
    m = (n - 1) // 2
    adj = [[int(G.has_edge(i,j)) for j in range(n)] for i in range(n)]
    d = [sum(row) for row in adj]
    return [[x * (d[i] - 1) - (m - 1) if i == j else (1 - x) * adj[i][j]
             for j in range(n)] for i in range(n)]


def det_polynomial(G):
    n = len(G)
    values = [bareiss_determinant(shifted_integer_rows(G, x)) for x in range(n + 1)]
    differences = []
    while values:
        differences.append(values[0])
        values = [values[i + 1] - values[i] for i in range(len(values) - 1)]
    # Newton forward differences: p(a)=sum Delta^k p(0)*binomial(a,k).
    polynomial = 0
    falling = 1
    for k, delta in enumerate(differences):
        if k:
            falling *= ALPHA - (k - 1)
        if delta:
            polynomial += sp.Rational(delta, sp.factorial(k)) * falling
    return sp.Poly(sp.expand(polynomial), ALPHA, domain=sp.QQ)


def shifted_sympy(G, alpha):
    n = len(G)
    m = (n - 1) // 2
    degree = dict(G.degree())
    return sp.Matrix(n, n, lambda i, j: (
        alpha * (degree[i] - 1) - (m - 1) if i == j
        else (1 - alpha if G.has_edge(i, j) else 0)))


def positive_zero_count(matrix):
    """Descartes is exact for a real-rooted polynomial (all roots real)."""
    co = matrix.charpoly().all_coeffs()
    zero = 0
    while co[-1] == 0:
        zero += 1
        co.pop()
    signs = [sp.sign(v) for v in co if v != 0]
    positive = sum(signs[i] != signs[i - 1] for i in range(1, len(signs)))
    return positive, zero


def interior_isolators(poly):
    """Exact disjoint intervals for distinct interior roots and multiplicities."""
    if poly.is_zero:
        return []
    # Remove endpoint factors (endpoints are not interior roots).
    while poly.degree() > 0 and poly.eval(0) == 0:
        poly = poly.exquo(sp.Poly(ALPHA, ALPHA))
    while poly.degree() > 0 and poly.eval(sp.Rational(1, 2)) == 0:
        poly = poly.exquo(sp.Poly(2 * ALPHA - 1, ALPHA))
    intervals = poly.intervals(eps=sp.Rational(1, 2**30))
    inside = []
    for (lo, hi), multiplicity in intervals:
        if lo > 0 and hi < sp.Rational(1, 2):
            inside.append((lo, hi, multiplicity))
        else:
            assert hi <= 0 or lo >= sp.Rational(1, 2), (poly, lo, hi)
    return inside


def rational_samples(roots):
    """Every component of (0,1/2) minus roots gets an exact sample.

    For non-rational roots the exact rational isolating interval has a
    sample on *both* sides of the root (its endpoints), as well as gaps.
    """
    if not roots:
        return [sp.Rational(1, 4)]
    out = []
    lower = sp.Rational(0)
    for lo, hi, _ in roots:
        assert lower < lo and lo <= hi
        out.append((lower + lo) / 2)
        if lo != hi:
            out.extend((lo, hi))
        lower = hi
    out.append((lower + sp.Rational(1, 2)) / 2)
    return sorted(set(out))


def main():
    print('method=integer Bareiss + finite differences + SymPy exact root isolation + real-rooted Descartes')
    print('python=' + sys.version.split()[0], 'networkx=' + nx.__version__, 'sympy=' + sp.__version__)
    atlas = nx.graph_atlas_g()
    for n in (3, 5, 7):
        count = identically = root_graphs = root_total = sample_total = 0
        equal_at_zero = 0
        endpoint_half = 0
        pos_max = 0
        examples = []
        for G in atlas:
            if len(G) != n:
                continue
            count += 1
            code = nx.to_graph6_bytes(G, header=False).decode().strip()
            p = det_polynomial(G)
            at_zero = positive_zero_count(shifted_sympy(G, sp.Integer(0)))
            if at_zero[0] <= 1 and ((at_zero[0] == 1 and at_zero[1] >= 1) or
                                     (at_zero[0] == 0 and at_zero[1] >= 2)):
                equal_at_zero += 1
            assert at_zero[0] <= 1, (n, code, 'alpha0', at_zero)
            if p.is_zero:
                identically += 1
                examples.append('identically=' + code)
                # The known H_m is the unique identically-zero polynomial.
                H = nx.disjoint_union(nx.complete_graph((n - 1) // 2),
                                      nx.complete_graph((n - 1) // 2))
                H.add_node(n - 1)
                H.add_edges_from((n - 1, v) for v in range(n - 1))
                assert nx.is_isomorphic(G, H), (n, code, 'unexpected permanent equality')
                continue
            if p.eval(sp.Rational(1, 2)) == 0:
                endpoint_half += 1
            roots = interior_isolators(p)
            if roots:
                root_graphs += 1
                root_total += len(roots)
                examples.append(code + ':' + str(len(roots)))
            # Every interior root must be simple to rule out tangencies
            # lambda_2(alpha)=threshold at an isolated parameter.
            assert all(mult == 1 for _, _, mult in roots), (n, code, roots)
            for alpha in rational_samples(roots):
                positive, zero = positive_zero_count(shifted_sympy(G, alpha))
                assert zero == 0, (n, code, alpha, 'sample is an algebraic root')
                assert positive <= 1, (n, code, alpha, positive)
                pos_max = max(pos_max, positive)
                sample_total += 1
        assert identically == 1, (n, identically)
        assert equal_at_zero == {3:3, 5:6, 7:11}[n], (n, equal_at_zero)
        print('n=%d graphs=%d alpha0_equal=%d identically_zero=%d interior_root_graphs=%d interior_simple_roots=%d isolated_samples=%d endpoint_half_zeros=%d max_positive=%d failures=0' %
              (n, count, equal_at_zero, identically, root_graphs, root_total, sample_total, endpoint_half, pos_max))
        print('n=%d roots_by_graph=%s' % (n, ','.join(examples)))
    print('THEOREM VERIFIED EXACTLY for n=3,5,7 and ALL real alpha in [0,1/2): lambda2<=m-1+alpha; equality for positive alpha only H_m.')


if __name__ == '__main__':
    main()
