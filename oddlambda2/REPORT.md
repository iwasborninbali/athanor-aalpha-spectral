# Odd-order maximal second A-alpha eigenvalue (`oddlambda2`)

**Audit date:** 2026-11-15. Athanor lab (Aleksei Kudriashov). **Evidence status:** finite rational-grid theorem certified by exhaustive exact integer inertia and independently checked characteristic polynomials; general odd-order conjecture **unproved**. All programs and log outputs below are new files in this directory. No paid services.

## Question

For simple undirected graphs of odd order `n=2m+1`, determine `max_G lambda_2(alpha D(G)+(1-alpha)A(G))` for `0 <= alpha < 1/2`. Nikiforov, [arXiv:1607.03015](https://arxiv.org/pdf/1607.03015), immediately after Proposition 31 (printed p.22), proves the general bound `lambda_2 <= n/2-1` and observes that `K_m + K_(m+1)` attains only `m-1` when `n` is odd.

The candidate value is `m-1+alpha`, attained by `H_m = K_1 join (K_m + K_m)` (`+` means disjoint union). This candidate and its full characteristic polynomial **are already in Xue--Liu--Yu--Shu (2020)**, not new discoveries.

## Results with evidence status

**R1 (PUBLISHED, independently rechecked):** J. Xue, R. Liu, G. Yu, J. Shu, *The Multiplicity of A_alpha-Eigenvalues of Graphs*, *Electronic Journal of Linear Algebra* **36** (2020), 645--657, [DOI:10.13001/ela.2020.4999](https://doi.org/10.13001/ela.2020.4999), p.653, **equation (3.7)**, gives for `n=2m+1`, `m>=2`, the exact factorization

```
(x-(m+1)*alpha+1)^(2m-2) * (x-m-alpha+1) *
[x^2+(1-m-2*alpha*m-alpha)*x + 2*alpha*m-2*m+2*alpha*m^2].
```

It follows for **all `0<=alpha<1`** that `lambda_2(A_alpha(H_m))=m-1+alpha`. Indeed the internal zero-sum eigenvalue is `(m+1)*alpha-1`, smaller by `m*(1-alpha)>0`; the remaining two eigenvalues are those of `[[2m*alpha,2m*(1-alpha)],[1-alpha,m-1+alpha]]`, whose characteristic polynomial at `m-1+alpha` is strictly negative. Hence one root exceeds the candidate and one falls below it. The `m=1` case follows directly with no zero-sum block. The Xue et al. paper's focus is multiplicity (including Lemma 3.5), **not a theorem that H_m globally maximizes lambda_2**. `sympy_check.py` independently expands (3.7) for `m=2,3,4`. For `alpha>1/2`, `K_n` gives `lambda_2=n*alpha-1 > m-1+alpha`, so no unrestricted-alpha optimality assertion is made.

**R2 (EXACT FINITE-GRID CERTIFICATE):** For **every unlabeled graph** at `n=5,7,9`, each of the following **13 rational values** in `[0,1/2)`, the maximum equals `(n-3)/2 + alpha`:

```
0, 1/20, 1/10, 1/8, 1/6, 1/5, 1/4, 1/3,
2/5, 9/20, 12/25, 49/100, 499/1000.
```

| order | unlabeled graphs | max at alpha=p/q | equality classes at alpha=0 | equality classes at each of 12 positive samples |
| --- | ---: | --- | ---: | ---: |
| 5 | 34 | `1+p/q` | 6 | 1 (H_2) |
| 7 | 1,044 | `2+p/q` | 11 | 1 (H_3) |
| 9 | 274,668 | `3+p/q` | 17 | 1 (H_4) |

**3,584,698 graph/parameter pairs** exactly certified, **zero upper-bound violations**. `exact_scan.log` includes every order, fraction, violation count, equality count and maximizer graph6; H_4 at n9 is graph6 `HJ]CKN~` (compared to an explicit H_4 via NetworkX isomorphism); there is no rounding or numerical eigenvalue comparison. **Not** a continuous-parameter or general-n proof.

**R3 (EXACT SECOND METHOD, scoped):** `sympy_check.py` independently forms the shifted characteristic polynomial over `Z` and counts positive roots by Descartes sign changes. Symmetry makes all its roots real; therefore the number of sign changes equals the exact number of strictly positive roots, counting multiplicity (since factors for negative roots have no positive real zero and the positive-root bound is attained). It independently checks **all 34+1044** graphs for five rational parameters (`0,1/4,2/5,49/100,499/1000`), the Xue polynomial on H_2,H_3,H_4, and **13 selected n9 graphs** for each of all 13 fractions. It checks the 274,668 generated n9 graph6 strings are distinct and fingerprints the generated set. **It does not independently enumerate all 274,668 n9 isomorphism classes or recompute all their characteristic polynomials.** Its log is `sympy_check.log`.

### Why the exhaustive certificate is mathematically exact

For rational `alpha=p/q`, `0<=p<q/2`, and `n=2m+1`, define the integer symmetric matrix

```
C_G(p,q) = q * (A_alpha(G) - (m-1+p/q)*I).
C_ii = p*(deg(i)-1) - q*(m-1)
C_ij = (q-p) if ij is an edge, 0 otherwise.
```

By the spectral theorem, `lambda_2(A_alpha(G)) <= m-1+p/q` **if and only if** `C_G` has *at most one strictly positive eigenvalue*. Equality means either inertia `(positive=1,zero>=1)` or `(positive=0,zero>=2)` (the latter accommodates a repeated top eigenvalue at the bound). `exact_scan.cpp` evaluates symmetric fraction-free Schur complements with unbounded `boost::multiprecision::cpp_int` and tracks pivot signs via the exact ratio of successive Bareiss pivots. A symmetric simultaneous row/column permutation moves nonzero diagonal pivots to the front. If the residual diagonal vanishes but an off-diagonal `C_ij` does not, the invertible congruence change `e_i -> e_i+e_j` produces a nonzero diagonal pivot. If all residual entries vanish, each remaining direction is in the kernel. Sylvester's inertia law proves the integer signs count the eigenvalue signs, including zero multiplicity. Every division has an exact-remainder assertion; no floating point, interval assumptions, or unbounded overflow risk. `sympy_check.py` checks many matrices by a different exact characteristic-polynomial route.

Completeness of the n9 candidate set: `enumerate9.c` starts from NetworkX's **complete** graph atlas on seven vertices (1044 graph6 lines), extends each graph by all `2^7` possible neighborhoods of an added vertex, canonically deduplicates n8, repeats with all `2^8` neighborhoods to n9, and produces **12,346** n8 and **274,668** n9 distinct codes. Every n8 graph loses some vertex to a n7 graph and every n9 graph loses some vertex to a n8 graph, so all classes are reached by induction. The canonical encoding sorts by the isomorphism-invariant `(degree, neighbor-degree histogram)` signature and searches **every permutation within each identical-signature class** for a minimum full adjacency code. Thus isomorphic input graphs yield the same canonical code and distinct isomorphism types cannot collide. Its generated n9 graph6 file has SHA-256 `7d9f2373e60e51f1170f5c9e93bcd76cea0491245eb95d00685bb98ba775aa9e` (2,197,344 bytes); it is deliberately generated on replay, **not committed**. Reproducing counts alone would not prove generator completeness; the preceding orderly-extension and canonicalization argument is essential.

## Priority check

Searched/read **before** calling any result new: this repository's `verify/REPORT.md`, `attack/REPORT.md`, `status/REPORT.md`, `priority/REPORT.md`; Nikiforov [arXiv:1607.03015](https://arxiv.org/pdf/1607.03015) and 2017 [journal PDF](https://doiserbia.nb.rs/img/doi/1452-8630/2017/1452-86301701081N.pdf); Xue--Liu--Yu--Shu (2020), [primary full text](https://journals.uwyo.edu/index.php/ela/article/download/4999/5171/6923), **eq. (3.7), Lemma 3.5**; Chen--Li--Meng, [*On the second largest A_alpha-eigenvalues of graphs* (2019)](https://doi.org/10.1016/j.laa.2019.06.027), which addresses lower bounds/connected minimizers rather than the all-graph upper extremum; and Liu--Das--Shu, [*On the eigenvalues of A_alpha-matrix of graphs* (2020)](https://doi.org/10.1016/j.disc.2020.111917). Searches included exact quoted title of Xue's paper with `cited by`, `2021`--`2026`, `odd second eigenvalue maximum`, and DOI `10.13001/ela.2020.4999`. Investigated [a 2024 diameter/size spectral-radius paper](https://doi.org/10.1080/03081087.2024.2329197), whose publisher metadata lists the Xue article as *related research* (not independently confirmed to be a formal citation), and [a 2024 line-graph eigenvalue paper](https://doi.org/10.1515/spma-2024-0016) (uses the **opposite** alpha convention in its abstract; no verified odd-order maximum). **No citing paper proving the proposed odd-order global maximum was verified in this bounded search; no claim that none exists.** The already-published H_m factorization and the previously documented n<=7 floating-point scans are **not** claimed as priority discoveries here. The verifiable incremental contribution is n9 exhaustive exact arithmetic plus rational-grid certificates.

## Replay and files

Free software: a C11 compiler, C++17 compiler with Boost multiprecision headers, Python 3, NetworkX and SymPy. From `oddlambda2/` run `sh run.sh`; it reconstructs `atlas5.g6`, `atlas7.g6` and `graphs9.g6` (not committed), compiles `enumerate9.c` and `exact_scan.cpp`, reruns the complete exact scan and independent polynomial tests, diffs every resulting log byte-for-byte, checks the generated graph6 fingerprint, and verifies `SHA256SUMS` for all committed inputs and logs. NetworkX's atlas class counts are part of the replay preconditions. All file hashes are SHA-256; `SHA256SUMS` intentionally does not hash itself.

## Not covered

- **No proof** of `lambda_2 <= (n-3)/2+alpha` for **arbitrary odd n** or for a **continuum** of `alpha`; no counterexample found within the specified finite grid only.
- The finite-grid uniqueness at positive samples does **not** imply uniqueness at all positive alpha, including arbitrarily close to zero.
- Independent n9 polynomial checks cover **13 graphs per alpha**, not all 274,668, and depend on the explicitly checked n9 generator; n9 exhaustivity is justified by the generation argument, not a separate enumerator.
- No classification of all alpha=0 equality graphs, no odd n>9 sweep, and **alpha=1/2 is excluded** from the stated question and exact scan.
- Citing-literature search bounded by accessible free full text and metadata. It cannot establish that the question is unaddressed everywhere.

## Next steps

1. **Highest priority: a structural general-n proof or witness counterexample.** For any graph `G` on `2m+1` vertices, show the exact shifted matrix `C_G(p,q)` has positive inertia at most one uniformly for `0<=p/q<1/2`. Seek a forbidden two-positive-direction principal subspace or an extremal edge-switching argument; the exact scan's rare near-extremizers can guide a star-complement/quotient reduction. Avoid assuming numerical grid monotonicity.
2. Develop a **continuous-alpha finite-order certificate**: isolate every possible inertia transition among all n<=9 graphs by resultant/discriminant/Sturm analysis of shifted determinants (including zero crossings), not interpolation from the sampled grid.
3. Independently regenerate n9 unlabeled classes with a second canonical graph generator (e.g. free `geng`/nauty when available), compare graph6 canonical sets by isomorphism-invariant signatures, and tighten the citation-network audit.
