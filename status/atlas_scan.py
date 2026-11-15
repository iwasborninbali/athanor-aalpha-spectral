#!/usr/bin/env python3
"""Exhaustive search of the NetworkX nonisomorphic graph atlas (orders 3,5,7).

Numerical screening only: no claim beyond n<=7 or sampled rational alpha.
"""
import hashlib
import sys
import numpy as np
import networkx as nx

ALPHAS = (0.0, 0.25, 0.49, 0.5, 0.75)
ATLAS = nx.graph_atlas_g()
print('atlas_scan v1')
print('python', sys.version.split()[0], 'numpy', np.__version__, 'networkx', nx.__version__)
print('atlas_count', len(ATLAS))
print('degrees/orders per odd n')
for n in (3, 5, 7):
    graphs = [(i, g) for i, g in enumerate(ATLAS) if len(g) == n]
    print('n', n, 'unlabeled_graphs', len(graphs))
    for a in ALPHAS:
        best2 = float('-inf')
        minlast = float('inf')
        i2 = imin = -1
        for i, g in graphs:
            M = nx.to_numpy_array(g, nodelist=range(n), dtype=float)
            d = M.sum(axis=1)
            B = (1 - a) * M + a * np.diag(d)
            vals = np.linalg.eigvalsh(B)
            if vals[-2] > best2 + 1e-9:
                best2, i2 = float(vals[-2]), i
            if vals[0] < minlast - 1e-9:
                minlast, imin = float(vals[0]), i
        print(f'alpha={a:.2f} max_lambda2={best2:.12f} witness_atlas={i2} '
              f'min_lambda={minlast:.12f} min_atlas={imin}')
    # Show one exact graph providing a lambda2 maximum at alpha=1/4
    best2 = -1e100
    witnesses = []
    for i, g in graphs:
        A = nx.to_numpy_array(g, nodelist=range(n))
        eig = np.linalg.eigvalsh(.75*A+.25*np.diag(A.sum(axis=1)))[-2]
        if eig > best2 + 1e-9:
            best2 = eig
            witnesses=[i]
        elif abs(eig-best2) < 1e-9:
            witnesses.append(i)
    print('alpha=.25 maximizer_ids', witnesses)
    print('example_edges', list(ATLAS[witnesses[0]].edges()))
print('END atlas_scan')
