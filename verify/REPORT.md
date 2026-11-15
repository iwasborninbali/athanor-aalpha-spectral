# Independent verification: Nikiforov, *Merging the A- and Q-spectral theories*

## Scope

- **Source examined:** V. Nikiforov, arXiv:1607.03015v1, 26-page PDF at <https://arxiv.org/pdf/1607.03015>, also compared with <https://arxiv.org/html/1607.03015v1>. Page numbers below are the **printed PDF page numbers**. The user-specified 2017 journal article may have editorial revisions; **the journal PDF was not separately checked**.
- **Object:** finite simple undirected graphs, `A_alpha = alpha*D + (1-alpha)*A`, `0 <= alpha <= 1`; descending real eigenvalues. A statement is **refuted** only with a specified graph, parameter, and a directly checkable inequality/equality. **Confirmed on a range** never means a proof for all graphs.
- **Full labeled-graph sweep:** `enumerate.py` visits all 33,867 graphs on `1 <= n <= 6` (1, 2, 8, 64, 1024, 32768 respectively), at `alpha = 0,1/4,1/2,3/4,1`; 169,335 graph/alpha cases. It implements edge masks, graph connectivity, chromatic number, clique number, exact brute-force maxcut, Turan/split candidates, and graph-distance BFS **without graph libraries**. NumPy's symmetric eigensolver supplies floating-point eigenvalues (tolerance `2e-8`).
- **Larger checks:** `n7_focus.py` scans all 2,097,152 seven-vertex labeled graphs at `alpha=1/4` for Proposition 29, the corrected Rayleigh bound, and Proposition 20; `random_probe.py` checks 6,000 reproducibly generated graphs of orders `7,8,10,12,16,20` at nine alphas (54,000 graph/alpha cases). `witness.py`, independently written using SymPy rational matrices and quotient equations, supplies **exact witnesses and spectrum tests**, including the page-6 table for all `2 <= k <= 100` and Proposition 37 for `1 <= b <= a <= 10` at five rational alphas.
- **Distinctions:** Statements about differentiability and uniqueness can be disproved by a single exact example. A generic extremal problem, an arbitrary-order theorem, or a proof in the paper is not established merely by the finite sweep. No generic "weight distribution" table occurs in this paper; the relevant displayed eigenvector weights are (24) and the weights `z_i=1/(r n_i)` in Theorem 26(iii), checked exactly on two multipartitions.

## Findings

### Refuted, with explicit witnesses

1. **Page 6, the only numerical table, column `K_{1,3k^2+1}`, row `lambda(A_{1/2})`: REFUTED.** The paper prints `(3k^2+1)/2` whereas the correct value for the star of order `3k^2+2` is `(3k^2+2)/2`. For `k=2`, `K_{1,13}` has radius **7**, not **13/2**. The other cells and the identity of the maximizing component at `alpha=0,1/2,1` check for `k=2..100`. `witness.py` uses the bipartite equitable quotient, not the paper's table.
2. **Question 5, p. 7: the universal differentiability assertion is FALSE.** Let `G=K3 disjoint_union K1,3`. The `K3` Perron root equals 2 for all alpha; the star root is `s(alpha)=[4alpha+sqrt(16alpha^2-24alpha+12)]/2`. From `s(alpha) < 2` iff `alpha < 1/2` and `s(alpha)>2` iff `alpha>1/2` (within `[0,1]`), `lambda_1(A_alpha(G))=max(2,s(alpha))`. Its left derivative at 1/2 is **0** and its right derivative **1**. This answers the question negatively for disconnected graphs, **not** the connected-graph variant.
3. **Proposition 8, p. 7, missing `r>=2` hypothesis: REFUTED as written.** An edgeless graph has chromatic number `r=1` and is 0-regular. For `alpha=1/4 < 1/r`, `A_alpha=0` is positive semidefinite. The proof divides by `r-1`, revealing the missing assumption. The nontrivial `r>=2` tests passed.
4. **Proposition 11, p. 8, missing `Delta>=1` and endpoint equality: REFUTED as written.** For `K1` (`Delta=0`) at `alpha=1/4`, the printed lower bound equals `1/4` but the radius is 0. For the connected `K3` at `alpha=1`, the bound equals the actual radius 2 although `K3` is not the star `K1,2`. The **inequality** passed every swept graph with `Delta>=1`.
5. **Corollary 12, p. 8: REFUTED for `Delta=0`.** For `G=K1` and `alpha=1/4`, the asserted `rho>=alpha(Delta+1)=1/4` gives `0>=1/4`. The positive-degree cases passed the exhaustive sweep.
6. **Unnumbered strict minimum-degree claim, p. 8, between Propositions 9 and 10: REFUTED with isolates.** At `alpha=3/4`, `G=K2 disjoint_union K1` has `delta=0` and `lambda_min=0=alpha*delta`, contradicting the printed strict `<` for `alpha<1`. The Rayleigh non-strict `<=` remains valid.
7. **Proposition 16 and Corollary 17, p. 10, endpoint `alpha=1`: REFUTED if their unqualified alpha includes 1.** For `G=K4-e` with adjacent universal vertices 0 and 1, they are automorphism-equivalent. At `alpha=1`, `A_1=diag(3,3,2,2)` and `x=(1,0,0,0)` is a Perron eigenvector but `x_0 != x_1`. **For `alpha<1` the intended Perron--Frobenius argument is unaffected.**
8. **Proposition 18, p. 11, equality clauses: REFUTED at endpoints.** For the connected path `P3`, `lambda(A_0)=lambda(A)=sqrt(2)` (equality in (15)) although no `lambda(A)`-regular component exists. At `alpha=1`, the upper bound (16) becomes `lambda(D)=Delta=2`, again equality although the graph has no `Delta`-regular component. The inequalities (15),(16) themselves passed.
9. **Textual conjectural expectation before Problem 22, p. 14, concerning (20),(21): REFUTED even at interior `alpha=1/4`.** Let `G=K1 join K2,2` with edge list `[(0,1),(0,2),(0,3),(0,4),(1,3),(1,4),(2,3),(2,4)]`. Degrees are `(4,3,3,3,3)` and neighbor-degree sums are `(12,10,10,10,10)`. Thus **every** `alpha*d(v)+(1-alpha)*w(v)/d(v)` equals `13/4` at `alpha=1/4`. The positive degree vector is an eigenvector of `A_alpha` with eigenvalue `13/4`, so the Perron root is `13/4`, with equality in BOTH (20) and (21) although the graph is connected and irregular. This refutes the *informal conjectural sentence*, **not** Proposition 20 (whose uniqueness claim is only for `1/2<alpha<1`).
10. **Theorem 26(ii), pp. 16--17, strict uniqueness at `alpha=1`: REFUTED as written.** The same `K1 join K2,2` is 3-partite of order 5, different from `S_{5,2}=K2 join E3`. Both have maximum degree 4 and hence `lambda(A_1)=4`. The claim becomes true if its high-alpha interval excludes the endpoint 1; Theorem 27(ii) *does* explicitly exclude 1. This is an equality-case/domain failure, not a violation of the extremal numeric bound.
11. **Proposition 29, p. 21: REFUTED throughout `0<=alpha<1` by `K3`.** The printed upper bound is `2alpha*m/n - 2(1-alpha)*maxcut/n`. For `K3`, `m=3`, `n=3`, `maxcut=2` and `lambda_min=3alpha-1`. Printed RHS=`(10alpha-4)/3`. Their difference is `(1-alpha)/3>0`, contradicting the claimed `lambda_min <= RHS`. At `alpha=0` this is `-1 <= -4/3`, false. **Corrected universal Rayleigh bound:**

   `lambda_min(A_alpha(G)) <= 2m/n - 4(1-alpha)*maxcut(G)/n`.

   **Complete proof of correction:** choose a bipartition attaining the maximum cut `c`, and let each coordinate of `x` be `+1` or `-1` according to its class. Then `x^T D x=2m`, `x^T A x=2[(m-c)-c]=2m-4c`, and `x^T x=n`. Thus `x^T A_alpha x/n = [2alpha*m + (1-alpha)(2m-4c)]/n = 2m/n-4(1-alpha)c/n`. Rayleigh's minimum principle proves the corrected bound. The line directly preceding Proposition 29 on p. 21 also has the same missing within-class contribution.
12. **Proposition 31(a), p. 22, equality iff `K_n`: REFUTED even inside the interval.** Let `G=K4-e` (degrees `(3,3,2,2)`) and `alpha=3/4`. Its exact spectrum in ascending order is `2-sqrt(2)/2, 3/2, 2, 2+sqrt(2)/2`. Therefore `lambda_2=2=alpha*n-1` although `G!=K4`. The inequality still passes all checked graphs. At `alpha=1`, the simpler diagonal witness gives `lambda_2=3=alpha*n-1`.
13. **Corollary 33, p. 22, absent restriction `alpha<1`: REFUTED at `alpha=1`.** The path `P3` has diameter `D=2`, but `A_1(P3)=diag(1,2,1)` has only **two** distinct eigenvalues, versus the claimed `D+1=3`. Proposition 32 explicitly assumes `alpha<1`; if that is meant to be inherited by Corollary 33, the counterexample disappears.
14. **Corollary 21, p. 13: NOT WELL-DEFINED for edgeless graphs.** Both (22) and (23) take a maximum/minimum over the empty edge set, despite the phrase "For any graph". For graphs with at least one edge, the sampled inequalities passed.

### Confirmed **on the tested ranges**, with status for each checkable claim

| Source / claim (printed page) | Evidence status and tested extent |
| --- | --- |
| Equations (1), (2), (6)--(11), pp. 3--6; Proposition 1, p. 5 | **Confirmed algebraically** from `A_alpha=alpha D+(1-alpha)A`; (11) also confirmed by `P36` tests on complete regular graphs. No independent general proof of Weyl's Theorem WS is claimed. |
| Proposition 2, p. 5 (component eigenvalue maxima/minima) | **Confirmed analytically** by block diagonal decomposition; page-6 table separately audited (one cell refuted). |
| Propositions 3--4, pp. 6--7 (monotone eigenvalues, `n`-Lipschitz) | **Confirmed numerically** for every `n<=6` graph, consecutive five-alpha grid points; 135,468 graph/interval comparisons of eigenvalue arrays. Convexity/concavity follow directly from max/min of affine Rayleigh quotients. |
| Proposition 6, p. 7 (PSD/PD for alpha>1/2) | **Confirmed** at alpha=3/4 and 1 for every graph `n<=6`; positive-definite checks exclude isolates. |
| Proposition 8, p. 7 | **Refuted** at `r=1` as above; **confirmed** for regular nonempty graph cases at sampled alpha below `1/r`. |
| Propositions 9--10, p. 8 (degree and Weyl bounds) | **Confirmed** with 1,010,065 eigenvalue-index checks for each side of Proposition 10 and for Proposition 9, across all `n<=6` graphs and five alphas. |
| Proposition 11 and Corollary 12, p. 8 | **Refuted at `Delta=0`**, and the P11 equality statement at alpha=1; **confirmed inequality** for `Delta>=1` throughout the full sweep. |
| Proposition 13, p. 9 (PF properties) | **Not exhaustively checked as a standalone assertion**; compatible with all scanned connected cases for alpha<1. Relies on established PF theory, not certified anew here. |
| Proposition 14, p. 9 (join quotient) | **Confirmed exactly** for `K2 join E3` at five rational alphas (`witness.py`), independently using the full symmetric matrix and the paper's 2x2 quotient. |
| Proposition 15, pp. 9--10 (edge moving) | **Not separately checked**; its universal quantification over Perron eigenvectors was not exhaustively sampled. |
| Proposition 16 and Corollary 17, p. 10 | **Refuted if alpha=1 included**, as above; no violation for alpha<1 established. |
| Proposition 18, p. 11 | **Bounds confirmed** on all `n<=6`, five-alpha cases; **endpoint equality clauses refuted**. |
| Corollary 19, p. 12 | **Confirmed inequality** across all `n<=6` graphs and five alphas; equality classifications not exhaustively enumerated. |
| Proposition 20 and Corollary 21, pp. 12--13 | **Bounds confirmed** in the defined domains, including the 141,315 cases with no isolates for P20. **Problem 22's informal stronger conjecture refuted**. C21 edgeless graph omission noted. |
| Proposition 23, p. 14 | **Confirmed**: lower and upper bounds on `lambda(A_alpha)^2`, 169,335 graph/alpha cases each. |
| Theorem 25, pp. 15--16 (bipartite extremal spectral radius) | **Confirmed finite inequalities** for every bipartite labeled graph `n<=6` on five-alpha grid; equality classifications not independently established for all n. |
| Theorem 26, pp. 16--17 (`r`-partite extremal spectral radius) | **Confirmed finite inequalities** for `n<=6`, `2<=r<=min(4,n)` on grid; **uniqueness refuted at alpha=1**. Its displayed threshold eigenvector with class weights `1/(r*n_i)` checked exactly for partitions `(2,1,3)` and `(4,2,1,1)`. |
| Theorem 27, pp. 18--19 (`K_{r+1}`-free extremal spectral radius) | **Confirmed finite inequalities** for all graphs `n<=6`, `2<=r<=min(4,n)` and `alpha<1`; no global theorem proved. |
| Lemma 28, pp. 19--20 (multipartite maximizers) | **Not independently checked as a structural or universal uniqueness theorem**; compatible with scanned bounds. |
| Proposition 29, p. 21 | **Refuted**; corrected inequality **confirmed** on all `n<=6` graph/alpha cases plus the exact general proof above. |
| Proposition 31, p. 22 | **Numeric inequalities (a),(b) confirmed** on all `n<=6` tested graphs; **(a)'s uniqueness refuted**, including alpha=3/4. The even-order `2K_{n/2}` equality in (b) is directly checked by the regular-clique formula. |
| Proposition 32, p. 22 (powers and distance) | **Confirmed** for all connected graphs `n<=5` at alpha 0,1/2,3/4 (22,563 positive-distance and 11,697 short-distance comparisons). |
| Corollary 33, p. 22 | **Confirmed** for connected graphs `n<=5` at alpha 0,1/2,3/4 (2,316 comparisons); **refuted at alpha=1**. |
| Propositions 34--35, p. 23 (trace and squared trace) | **Confirmed** by two calculations per case, 169,335 cases each, with eigenvalue sums checked against independent degree and edge-count expressions. |
| Proposition 36, p. 23 (`K_n` spectrum) | **Confirmed exactly**, `2<=n<=10`, five rational alphas, 45 matrices. |
| Propositions 37--38, pp. 23--24 (`K_{a,b}` and star spectrum) | **Confirmed exactly**, `1<=b<=a<=10`, five rational alphas, 275 spectra, comparing a full SymPy matrix with symbolic closed-form roots and multiplicities. |

### Open problems and finite observations

- **Problem 7, p. 7** (least PSD alpha): **not solved in general**. Exact examples: `K_n` has threshold `1/n`, `K_{a,b}` with both parts nonempty has threshold `1/2`; isolates do not affect PSD of nontrivial components.
- **Problem 22, p. 14** (all equality cases of (20)--(23)): **not classified**; the explicit `alpha=1/4`, `K1 join K2,2` example refutes the stated informal expectation. `enumerate.py` recorded 477 labeled connected irregular graph cases attaining the P20 upper and lower row bounds at alpha=1/4 among `n<=6` (relabelings count separately).
- **Problem 24, p. 15** (forbidden paths/cycles): **not solved or exhaustively tested** for specific forbidden path/cycle choices; no structural classification attempted.
- **Problem 30, p. 21** (least eigenvalue for `0<alpha<1/2`): **not solved in general**. At alpha=1/4, exhaustive minima for orders 2--6 are `-1/2,-0.693000468,-1,-1.216364983,-3/2` respectively, all attained by balanced complete bipartite graphs (graph masks in `exhaustive.log`). This is a finite observation, not a proof of general extremality.
- **Odd-order maximum of `lambda_2` under Proposition 31(b), p. 22**: **not solved**; the stated gap remains unclosed by this study. The tested upper bound itself passed.

## Evidence

Run locally in any Python environment with `numpy` and `sympy` (versions in this container: NumPy 2.3.5, SymPy 1.14.0). Commands and their exact standard-output records:

```sh
OPENBLAS_NUM_THREADS=1 python verify/enumerate.py --nmax 6 > verify/exhaustive.log
OPENBLAS_NUM_THREADS=1 python verify/n7_focus.py > verify/n7_focus.log
OPENBLAS_NUM_THREADS=1 python verify/random_probe.py > verify/random_probe.log
OPENBLAS_NUM_THREADS=1 python verify/witness.py > verify/witness.log
sha256sum verify/*.py verify/*.log verify/REPORT.md
```

The controlled scan seed is `20261115`. Graph masks in `exhaustive.log` are **labeled-graph** masks: bit order `[(0,1),(0,2),...,(0,n-1),(1,2),...,(n-2,n-1)]`. The primary examples also appear as human-readable edge lists and exact rational calculations in `witness.py`. The numerical method computes eigenvalues using LAPACK (`numpy.linalg.eigvalsh`); the independent exact witness checker uses SymPy rational characteristic polynomials, a quotient matrix or a positive eigenvector. `SHA256SUMS` records SHA-256 hashes of **each program, output log, and this report** (not of itself). It is generated after all other files are finalized. Scan logs are not edited.

For a check that does not require software, the central P29 contradiction is `G=K3`, `alpha=0`, `maxcut=2`, `lambda_min=-1`, claimed upper bound `-4/3`; the Proposition 31 contradiction is `G=K4-e`, `alpha=3/4`, exact eigenvalues listed above; and the Problem 22 exception is the graph/degree calculation above. These supply independent, manually replayable witnesses even if the numerical sweeps are not run.

## Not covered

- General graphs of unbounded order, irrational parameter values, derivatives for *every* eigenvalue of every **connected** graph, or a formal exhaustive isomorphism classification.
- A machine-checkable proof of the spectral Turan theorems, Lemma 28, Perron--Frobenius theory, or Weyl's theorem; we report finite tests without calling them proofs.
- The later journal typesetting: only the arXiv v1 PDF was fetched/read, so editorial corrections in the 2017 published article are not ruled out.
- A full solution of Problems 7, 22, 24, 30 or the odd-order question for `lambda_2`; the parameter searches and exact witnesses do not imply one.
- Exhaustion of all graphs beyond 7 vertices. The large-order probe is seeded sampling, not an exhaustive test; small floating-point near-equalities are subject to numeric tolerance.

## Next steps

1. Compare every cited arXiv-v1 discrepancy against the **journal PDF** and any published errata; do not attribute discrepancies to the journal without reading it.
2. Amend Proposition 29's formula/proof; publish the `K3` counterexample and corrected Rayleigh bound. Amend the p. 6 table entry.
3. Correct equality and domain statements: Proposition 31(a) requires substantive revision beyond an endpoint exclusion; Proposition 20's informal equality conjecture must allow `K1 join K2,2` at alpha=1/4; Theorem 26(ii), Proposition 11, Proposition 18, Proposition 16/Corollary 17, Corollary 33 and Proposition 8 need explicit hypotheses or restricted parameter ranges.
4. Analyze the equality equations for Proposition 20 by solving `alpha*d(v)+(1-alpha)*w(v)/d(v)=constant` for integer-degree patterns; investigate systematic families beyond `K1 join K2,2`.
5. Use a canonical unlabeled-graph generator or a compiled bitset implementation to scale structural classifications and odd-order `lambda_2` searches, independently replaying witnesses with exact rational methods.
