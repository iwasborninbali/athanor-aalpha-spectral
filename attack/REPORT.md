# Computational attack on Nikiforov's A_alpha open problems

## Scope

Primary source: V. Nikiforov, "Merging the A- and Q-spectral theories", arXiv:1607.03015v1, https://arxiv.org/pdf/1607.03015; published in Applicable Analysis and Discrete Mathematics 11 (2017), 81-107. Page numbers below are the printed pages of the arXiv version. The primary PDF was read using a public web viewer; a container download failed owing to disabled DNS.

Inventory of explicit open items (faithful paraphrases, not verbatim extracts; use the linked source to check the complete wording):
- Question 5, Section 3.2, p.7: Is the k-th ordered eigenvalue of A_alpha(G) differentiable with respect to alpha?
- Problem 7, Section 3.3, p.7: Find the least alpha for which A_alpha(G) is positive semidefinite.
- Problem 22, Section 4.3, p.14: Classify equality in inequalities (20), (21), (22), and (23). Just before this problem, the paper suggests (without a formally numbered conjecture) that equality in (20),(21) should require regularity except at alpha=0 or 1/2.
- Problem D, Section 5, p.15: Maximize the A_alpha spectral radius of an order-n graph forbidding a specified subgraph F. Historical Problems A, B, and C on that page are background formulations, not new A_alpha conjectures. Theorem 27 solves the complete-graph case of D.
- Problem 24, Section 5, p.15: Study Problem D when the forbidden subgraph is a prescribed path or cycle.
- Problem 30, Section 6.1, p.21: For 0<alpha<1/2, minimize the smallest A_alpha eigenvalue over n-vertex graphs.
- Unnumbered question following Proposition 31(b), Section 6.2, p.22: Find the exact maximum second eigenvalue for odd n and 0<=alpha<1/2.

No item explicitly labeled "Conjecture" was found in the checked arXiv text. The informal regularity expectation in Section 4.3 is considered separately below.

Computations: search.py enumerates every unlabeled graph through seven vertices from NetworkX graph_atlas_g (the order-by-order counts for n=0,...,7 are 1,1,2,4,11,34,156,1044). It uses symmetric numpy.linalg.eigvalsh and a 1e-8 comparison tolerance. labelled_crosscheck.py constructs all labelled graphs from scratch, without NetworkX, for n=2,...,6 by iterating 2^(n choose 2) edge masks. checker.py implements its own graph6 reader and uses fractions.Fraction for independent exact witness verification. All programs were written for this attack. Environment: Python 3.13.5, NetworkX 3.6.1, NumPy 2.3.5. No paid services were used.

Replay from repository root:
    OPENBLAS_NUM_THREADS=1 python3 attack/search.py > attack/search.log
    OPENBLAS_NUM_THREADS=1 python3 attack/labelled_crosscheck.py > attack/labelled_crosscheck.log
    python3 attack/checker.py > attack/checker.log
    sha256sum -c attack/SHA256SUMS

## Findings

### Question 5: exact negative answer, with proof

Take G=K2 disjoint union K1, graph6 B_, edges {(0,1)}. Its A_alpha eigenvalues are exactly 1, 2alpha-1, and 0. Therefore the second largest eigenvalue is max(0,2alpha-1), whose left derivative at alpha=1/2 is 0 while the right derivative is 2. This graph is disconnected; the analogous question restricted to connected graphs has not been settled here. The exact eigenvectors are verified by checker.py.

### Problem 7: exact spectral reformulation, with proof

Delete isolated vertices and set B=D^(-1/2) A D^(-1/2) and mu=lambda_min(B). The identity
    A_alpha = D^(1/2) [alpha I+(1-alpha) B] D^(1/2)
is a congruence, so positive semidefiniteness is equivalent to alpha+(1-alpha)*mu >= 0. Therefore the least alpha is -mu/(1-mu) whenever G has an edge, and 0 for an edgeless graph. Here -1<=mu<0. In particular, a graph with a nontrivial bipartite component has threshold 1/2. This is a complete exact spectral characterization, not a closed combinatorial classification and not claimed to be novel. Numeric checks: K4 1/4; P4 1/2; C5 sqrt(5)/5; K1,4 1/2.

### Problem 22: infinite exact counterexample family to the informal regularity conjecture

Let t,m>=2 and G=K_t join I_m, a connected nonregular split graph. Write h=t+m-1. The t clique vertices have degree h; the m independent vertices have degree t. The two types of row sums of D^(-1) A_alpha D, whose maximum/minimum give bounds (20)/(21), are:
    R_high = alpha*h + (1-alpha)*((t-1)*h + m*t)/h
    R_low  = alpha*t + (1-alpha)*h.
Their difference is
    R_high-R_low = (m-1)*((2*m+t-1)*alpha-m)/h.
At alpha=m/(2m+t-1), strictly between 0 and 1/2, these row sums agree and equal L=h-alpha*(m-1). Equivalently A_alpha*d=L*d for the strictly positive degree vector d. Since G is connected and alpha<1, Perron-Frobenius implies lambda(A_alpha)=L. Both inequalities (20) and (21) are equalities. This proves that the proposed regular-only criterion is false for an infinite family; it does not classify all equality cases and says nothing about equality in (22)/(23).

A four-vertex explicit witness is K4 minus one edge, graph6 C|, edges (0,1),(0,2),(0,3),(1,2),(2,3). Degrees are (3,2,3,2), alpha=2/5, and lambda=13/5 exactly. The standalone checker verifies this and several other t,m values using rational arithmetic.

Atlas search over the 980 connected nonregular graphs with 2<=n<=7 found 11 equality cases for (20)/(21) at alpha=.25, 6 at .40, and none at .60 or .75. It found no nonregular equality in (22),(23) at these four sampled values, by tolerance 1e-8. The paper already proves regularity for connected graphs when 1/2<alpha<1 (Proposition 20). No contradiction with that theorem is asserted.

### Odd order after Proposition 31: exact improved construction, proof and finite maxima

For n=2s+1, s>=1, let H_s=K1 join (K_s disjoint union K_s). The hub has degree 2s and the other vertices degree s. A vector equal to +1 on one clique, -1 on the other, and 0 at the hub has A_alpha eigenvalue b=s-1+alpha. Internal zero-sum vectors within each clique have eigenvalue c=(s+1)alpha-1, with multiplicity 2(s-1). The remaining two eigenvalues come from the hub/equal-peripheral block
    [[2s*alpha,2s*(1-alpha)],[(1-alpha),b]].
For 0<=alpha<1, this block has one eigenvalue above b and one below b, while c<b. Hence, exactly,
    lambda_2(A_alpha(H_s)) = s-1+alpha.
The paper's example K_s disjoint union K_(s+1) only gives s-1. Thus this construction gives a stronger lower bound by alpha for every positive alpha in the paper's range; it does not prove a universal upper bound.

All unlabeled graphs for n=3,5,7 were exhaustively tested at alpha=0,.25,.40,.49. At each positive sample the maximizer was uniquely H_s up to isomorphism. For n=7 and alpha=.40, lambda_2=2.4 with graph6 FK|ko. Independent labelled enumeration confirms n=3 and n=5. Only a sampled grid, not every real alpha, was checked.

### Problem 30: finite exhaustive minima

For each 2<=n<=7 and alpha in {.25,.40,.49}, the unique minimizing unlabeled graph was the balanced complete bipartite graph K_floor(n/2),ceil(n/2). Its eigenvalue agrees with Proposition 37. For n=7 the three minima were -1.7260815058, -0.6880613018, and -0.0685988899. A separate enumeration of all labelled graphs through n=6 reproduces those minima. No general theorem is claimed.

### Problem 24: forbidden paths/cycles at n=7

All 1044 unlabeled 7-vertex graphs were tested against noninduced forbidden P4, P5, P6, C3, C4, C5, C6 at alpha=0,.25,.40,.49,.60,.75. search.log records the admissible count, maximum spectral radius, tie count, and a graph6 maximizer for every one of these 42 cases.

Illustration: for P5-free graphs at alpha=.25, the maximum is 3, shared by four isomorphism classes, one being K4 disjoint union K3 (graph6 FKdbG). At alpha=.40, the unique maximizer is a 7-vertex star with one added edge joining two leaves (graph6 FsaC_), spectral radius 3.2821632626. The separate checker confirms neither graph contains a noninduced P5. This shows a grid-level change in extremizer, not its exact transition point. The C6-free maximizer also changes between alpha=.40 and .49.

## Evidence

Files contain complete executable sources and the logs of the runs, with SHA-256:

| New file | SHA-256 |
| --- | --- |
| attack/search.py | 363c3ec132e876807b1b0fa4bf2b9946f0b386d302f47bbfd3bfe7684cf04135 |
| attack/search.log | 47c31a2f3ab7e973781a5a16dc2a4f2296526aeb9ea86891cd79f98cb034557a |
| attack/checker.py | c3776626348837f9950990456bf72f0fa6ed34c1602ca18a341422727294f532 |
| attack/checker.log | 1049e3d7559dbace9d6120e9f499fdb2e5c551fd3f995b4120c578b34e430e67 |
| attack/labelled_crosscheck.py | 100392e759907f01ad9c1f684b739689265c8f613632122e8bbc9f0b2a041e5b |
| attack/labelled_crosscheck.log | 4e1e7fd235f87338c211d269ef3aabef4290a5872cfb815ca3de2ecf9c290015 |

The SHA256SUMS file is the machine-checkable manifest. The entire repo delta must contain only newly added files inside attack/.

## Not covered

There is no proof of the global optimum for Problems 24/30 or the odd-n lambda_2 maximum, no classification of all equality cases for Problem 22, and no verification of arbitrary forbidden F or higher orders. No exact phase-transition alpha for the order-7 extremizers was certified. The scans sample isolated alpha values rather than checking all alpha in a continuum. The numerical searches use floating-point eigvalsh with tolerance 1e-8 and rely on Atlas completeness; the exact statements are accompanied by written proofs and rational checks. A connected-only Question 5 remains untested. No systematic literature-priority search beyond the cited primary paper was performed; do not construe these computations as evidence of novelty in subsequent literature. The source PDF could be read through the web tool but not downloaded into the CPU container because external DNS was unavailable.

## Next steps

1. Classify all connected equality cases in bounds (20),(21) for 0<alpha<1/2, using the degree-vector eigenfunction condition and the split-graph family.
2. Either prove the putative odd-order bound lambda_2 <= (n-3)/2+alpha or search graph families beyond the seven-vertex atlas for counterexamples.
3. Check Problem 30 using exact polynomial inequalities and n=8,9 targeted non-bipartite constructions.
4. Locate and certify exact crossing values among order-7 P5/P6/C6 extremizers, then test larger orders.
