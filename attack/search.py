#!/usr/bin/env python3
"""Reproducible atlas enumeration for Nikiforov arXiv:1607.03015.
Only simple UNLABELLED graphs n<=7; NetworkX graph_atlas_g is the universe.
All numerical comparisons use symmetric eigvalsh and an explicit tolerance.
"""
import sys
from collections import defaultdict
import networkx as nx
import numpy as np

ALPHAS = (0.0, .25, .40, .49, .60, .75)
TOL = 1e-8


def graph6(g):
    return nx.to_graph6_bytes(g, header=False).decode('ascii').strip()


def values(g, alpha):
    a = nx.to_numpy_array(g, dtype=float)
    d = a.sum(axis=1)
    ev = np.linalg.eigvalsh(alpha*np.diag(d)+(1-alpha)*a)
    return ev, d, a


def contains_path(g, order):
    # Noninduced P_order (the forbidden-subgraph convention in Problem D).
    if order > len(g):
        return False
    adj = {v:set(g[v]) for v in g}
    def dfs(v, used):
        if len(used) == order:
            return True
        return any(dfs(w, used | {w}) for w in adj[v]-used)
    return any(dfs(v, {v}) for v in g)


def contains_cycle(g, order):
    # Simple cycle of exactly 'order', not necessarily induced.
    if order > len(g):
        return False
    adj = {v:set(g[v]) for v in g}
    def dfs(start, v, used):
        if len(used) == order:
            return start in adj[v]
        return any(dfs(start, w, used | {w}) for w in adj[v]-used)
    return any(dfs(v, v, {v}) for v in g)


def equality_gaps(g, alpha, ev, d, a):
    # Eqs (20)-(23), Prop. 20 and Corollary 21, arXiv pages 12-14.
    if len(g)==0 or min(d)==0:
        return None
    lamb = ev[-1]
    avg = a@d/d
    row = alpha*d + (1-alpha)*avg
    oriented = [alpha*d[u]+(1-alpha)*d[v]
                for u,v in g.edges() for u,v in ((u,v),(v,u))]
    return (max(row)-lamb, lamb-min(row),
            max(oriented)-lamb, lamb-min(oriented))


def fmt(v):
    return f'{v:.10f}'


def main():
    atlas = nx.graph_atlas_g()
    counts = defaultdict(int)
    for g in atlas: counts[len(g)] += 1
    print('PROVENANCE arxiv:1607.03015v1; atlas=networkx.graph_atlas_g; graph6=no_header')
    print('PYTHON', sys.version.split()[0], 'NETWORKX', nx.__version__, 'NUMPY', np.__version__)
    print('TOL', TOL, 'ALPHAS', ','.join(str(a) for a in ALPHAS))
    print('ATLAS', ' '.join(f'n{n}={counts[n]}' for n in range(8)), 'total', len(atlas))
    assert [counts[n] for n in range(8)] == [1,1,2,4,11,34,156,1044]
    graphs = [g for g in atlas if 2 <= len(g) <= 7]
    evcache = {(graph6(g),alpha): values(g,alpha) for g in graphs for alpha in ALPHAS}

    # Q5: exact analytically verifiable spectral crossing.
    q5 = nx.Graph([(0,1)])
    q5.add_node(2)
    print('Q5_WITNESS n=3 graph6='+graph6(q5)+' edges=0-1 eigenvalues={1,0,2alpha-1}')

    # P7: congruence and normalized adjacency threshold; find sample extrema.
    print('P7_THRESHOLD examples:')
    for name, g in [('K4',nx.complete_graph(4)), ('P4',nx.path_graph(4)),
                    ('C5',nx.cycle_graph(5)), ('K1,4',nx.star_graph(4))]:
        a = nx.to_numpy_array(g)
        d = a.sum(axis=1)
        b = a / np.sqrt(np.outer(d,d))
        mu = np.linalg.eigvalsh(b)[0]
        threshold = -mu/(1-mu)
        x = np.linalg.eigvalsh(threshold*np.diag(d)+(1-threshold)*a)[0]
        print(' ',name,'n',len(g),'graph6',graph6(g),'mu_min',fmt(mu),
              'alpha_psd',fmt(threshold),'lambda_min_at_alpha',fmt(x))

    # P22: only connected nonregular graphs; report any numerical equality.
    print('P22_CONNECTED_NONREGULAR_CHECK')
    for alpha in (.25,.40,.60,.75):
        tested = 0
        minima = [float('inf')]*4
        near=[]
        for g in graphs:
            if not nx.is_connected(g): continue
            key=graph6(g)
            ev,d,a=evcache[(key,alpha)]
            if max(d)==min(d): continue
            tested += 1
            gaps=equality_gaps(g,alpha,ev,d,a)
            for i,x in enumerate(gaps): minima[i]=min(minima[i],x)
            if min(gaps) < TOL: near.append((key, tuple(gaps)))
        print(' alpha',alpha,'graphs',tested,'min_gaps_20_21_22_23',
              '/'.join(fmt(x) for x in minima),'near_equalities',len(near))
        for x in near[:8]: print('  NEAR',x)

    # P30 and implicit odd-n Prop.31 optimization: complete search atlas n<=7.
    print('P30_MIN_LAMBDA_MIN: n alpha best_value balanced_Kab_value gap tie_count example_graph6')
    for n in range(2,8):
        gs = [g for g in graphs if len(g)==n]
        kab=nx.complete_bipartite_graph(n//2,n-n//2)
        for alpha in (.25,.40,.49):
            found = [(evcache[(graph6(g),alpha)][0][0],graph6(g)) for g in gs]
            minval = min(x[0] for x in found)
            wins=sorted(x[1] for x in found if abs(x[0]-minval)<TOL)
            kabval=values(kab,alpha)[0][0]
            print(' ',n,alpha,fmt(minval),fmt(kabval),fmt(minval-kabval),
                  len(wins),wins[0])

    print('ODD_LAMBDA2_MAX: n alpha max_value clique_union_value gap tie_count example_graph6')
    for n in (3,5,7):
        gs=[g for g in graphs if len(g)==n]
        union=nx.disjoint_union(nx.complete_graph(n//2),nx.complete_graph(n-n//2))
        for alpha in (0.0,.25,.40,.49):
            found=[(evcache[(graph6(g),alpha)][0][-2],graph6(g)) for g in gs]
            maxval=max(x[0] for x in found)
            wins=sorted(x[1] for x in found if abs(x[0]-maxval)<TOL)
            ref=values(union,alpha)[0][-2]
            print(' ',n,alpha,fmt(maxval),fmt(ref),fmt(maxval-ref),len(wins),wins[0])

    # P24: all n=7 atlas graphs, F a noninduced path/cycle, alpha grid.
    print('P24_N7_EXTREMA: F alpha admissible max_lambda ties example_graph6 edge_count')
    gs=[g for g in graphs if len(g)==7]
    for kind,k in [('P',4),('P',5),('P',6),('C',3),('C',4),('C',5),('C',6)]:
        admissible=[g for g in gs if not (contains_path(g,k) if kind=='P' else contains_cycle(g,k))]
        for alpha in ALPHAS:
            found=[(evcache[(graph6(g),alpha)][0][-1],graph6(g),g.number_of_edges()) for g in admissible]
            maxval=max(x[0] for x in found)
            wins=sorted(x for x in found if abs(x[0]-maxval)<TOL)
            print(' ',kind+str(k),alpha,len(admissible),fmt(maxval),len(wins),wins[0][1],wins[0][2])

if __name__=='__main__':
    main()
