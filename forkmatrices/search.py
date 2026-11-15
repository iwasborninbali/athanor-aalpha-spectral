"""Exact atlas screening for degree/transmission weighted row-sum equality.

Run: OPENBLAS_NUM_THREADS=1 python forkmatrices/search.py
Requires networkx (graph_atlas_g); no paid services or numerical eigensolver.
"""
from collections import Counter
from fractions import Fraction
import math
import networkx as nx


def candidate_alpha(s, p):
    """Find alpha in [0,1) making alpha*s+(1-alpha)*p constant."""
    a = None
    for i in range(1, len(s)):
        delta = p[i] - p[0]
        denom = Fraction(s[0] - s[i]) + delta
        if not denom:
            if delta:
                return None
            continue
        b = delta / denom
        if not (0 <= b < 1):
            return None
        if a is None:
            a = b
        elif a != b:
            return None
    if a is None:
        # All alpha work only for (weighted-)regular cases; excluded here.
        return None
    return a if len({a * s[i] + (1-a) * p[i] for i in range(len(s))}) == 1 else None


def dist_matrix(g):
    lens = dict(nx.all_pairs_shortest_path_length(g))
    return [[lens[i][j] for j in range(len(g))] for i in range(len(g))]


counts = Counter()
examples = {}
for atlas_id, raw in enumerate(nx.graph_atlas_g()):
    n = raw.number_of_nodes()
    if n < 2 or not nx.is_connected(raw):
        continue
    counts['connected'] += 1
    g = nx.convert_node_labels_to_integers(raw)
    adj = [[int(g.has_edge(i,j)) for j in range(n)] for i in range(n)]
    deg = [sum(row) for row in adj]
    if len(set(deg)) == 1:
        continue
    counts['nonregular'] += 1
    walk = [sum(adj[i][j]*deg[j] for j in range(n)) for i in range(n)]
    q = [Fraction(deg[i]*deg[i]+walk[i],deg[i]) for i in range(n)]
    if len(set(q)) == 1:
        counts['Q_degwitness'] += 1
        examples.setdefault('Q_degwitness',(atlas_id, nx.to_graph6_bytes(g,header=False).decode().strip()))
    # For degree-weighted Sombor matrix, test real-valued row sums with strict tolerance.
    so = [sum(math.sqrt(deg[i]**2+deg[j]**2)*deg[j] for j in range(n) if adj[i][j])/deg[i] for i in range(n)]
    if max(so)-min(so) < 1e-9:
        counts['Sombor_degwitness_numeric'] += 1
    z = dist_matrix(g)
    tr = [sum(row) for row in z]
    if len(set(tr)) == 1:
        continue
    counts['nontransmission_regular'] += 1
    two = [sum(z[i][j]*tr[j] for j in range(n)) for i in range(n)]
    dr = [Fraction(two[i],tr[i]) for i in range(n)]
    if len(set(dr)) == 1:
        counts['distance_trwitness'] += 1
        examples.setdefault('distance_trwitness',(atlas_id,nx.to_graph6_bytes(g,header=False).decode().strip()))
    dqr = [Fraction(tr[i])+dr[i] for i in range(n)]
    if len(set(dqr)) == 1:
        counts['distanceQ_trwitness'] += 1
    a = candidate_alpha(tr, dr)
    if a is not None:
        counts['Dalpha_interior_or_zero'] += 1
        if a == 0:
            counts['Dalpha_alpha0'] += 1
        elif a < Fraction(1,2):
            counts['Dalpha_alpha_below_half'] += 1
        else:
            counts['Dalpha_alpha_at_least_half'] += 1
        examples.setdefault('Dalpha_interior_or_zero',(atlas_id,nx.to_graph6_bytes(g,header=False).decode().strip(),str(a)))

print('method=NetworkX graph_atlas_g, connected unlabeled n=2..7, exact Fraction except Sombor float')
print('networkx_version='+nx.__version__)
for key in ['connected','nonregular','nontransmission_regular','Q_degwitness','distance_trwitness','distanceQ_trwitness','Dalpha_interior_or_zero','Dalpha_alpha0','Dalpha_alpha_below_half','Dalpha_alpha_at_least_half','Sombor_degwitness_numeric']:
    print(f'{key}={counts[key]}')
for key in sorted(examples):
    print(f'example_{key}={examples[key]}')
print('Sombor threshold=1e-9; finite test, not a proof')
