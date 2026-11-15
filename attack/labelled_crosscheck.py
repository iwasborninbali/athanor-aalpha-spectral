#!/usr/bin/env python3
"""Independent enumeration of *labelled* simple graphs, n=2..6.
Does not use NetworkX/atlas or graph6. Numpy eigvalsh is still a common
numerical dependency, so exact witnesses are checked by checker.py.
"""
import numpy as np


def main():
    print('LABELLED_ENUMERATION n=2..6 all 2^(n choose 2) edge masks')
    for n in range(2,7):
        pairs=[(i,j) for j in range(1,n) for i in range(j)]
        total=1<<len(pairs)
        for alpha in (.25,.40,.49):
            bestmin=float('inf'); bestsecond=float('-inf')
            min_count=0; second_count=0
            for mask in range(total):
                a=np.zeros((n,n),dtype=float)
                # Independent direct labelled-bitmask graph constructor.
                for bit,(i,j) in enumerate(pairs):
                    if mask & (1<<bit): a[i,j]=a[j,i]=1.0
                d=a.sum(axis=1)
                spec=np.linalg.eigvalsh(alpha*np.diag(d)+(1-alpha)*a)
                mn=float(spec[0]); second=float(spec[-2])
                if mn < bestmin-1e-8: bestmin=mn;min_count=1
                elif abs(mn-bestmin)<1e-8: min_count+=1
                if second>bestsecond+1e-8:bestsecond=second;second_count=1
                elif abs(second-bestsecond)<1e-8:second_count+=1
            p=n//2;q=n-p
            # Balanced complete bipartite eigenvalue, Prop.37.
            ref=(alpha*n-(alpha*alpha*n*n+4*p*q*(1-2*alpha))**.5)/2
            # Odd n lower-bound construction H_s.
            lower=(n-3)/2+alpha if n%2 else None
            assert abs(bestmin-ref)<1e-8
            if n%2: assert abs(bestsecond-lower)<1e-8
            print(f'n={n} alpha={alpha:.2f} labelled={total} lambda_min={bestmin:.10f} '
                  f'balanced_Kab={ref:.10f} minimizing_masks={min_count} '
                  f'lambda2_max={bestsecond:.10f} odd_Hs_lower={lower} maximizing_masks={second_count}')
    print('LABELLED CROSSCHECK SUCCESS')

if __name__=='__main__': main()
