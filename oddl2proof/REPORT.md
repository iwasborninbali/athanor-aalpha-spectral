# oddl2proof: odd-order second A-alpha eigenvalue

Athanor lab (Aleksei Kudriashov), 2026-11-15. Base: `main` of `iwasborninbali/athanor-aalpha-spectral`, archived v1.0 DOI [10.5281/zenodo.24239159](https://doi.org/10.5281/zenodo.24239159). This work is additive: no alteration of previous certificates. Simple undirected graphs throughout; `lambda_2` means second largest eigenvalue, multiplicities counted. All tools are free.

## Question

For `n=2m+1 >=3`, `A_a(G)=a D(G)+(1-a)A(G)` and `0<=a<1/2`, prove `lambda_2(A_a(G)) <= m-1+a` with equality for `a>0` only when `G=H_m=K_1 join (K_m union K_m)`, classify equality at `a=0`, or identify obstructions. Also study `1/2<=a<1`.

## Results with evidence status

### R1. Arbitrary odd order, a whole interval next to zero (PROVED; prior-art dependent)

**Theorem.** For each `m>=1` there exists `epsilon_m>0` such that for **every** graph of order `2m+1` and **every** real `0<a<epsilon_m`,

`lambda_2(A_a(G)) <= m-1+a`, with equality **if and only if** `G` is isomorphic to `H_m`.

**Complete deduction.** We use as an explicitly credited input the adjacency extremizer classification of Lin--Liu--Yin--You (2026), Theorem 1.4, specialized to all simple graphs; see R3. In particular every **connected** graph satisfying `lambda_2(A(G))=m-1` consists of two disjoint `K_m` blocks `U,V`, and one extra vertex `w` joined to `r>=1` vertices of `U` and `s>=1` vertices of `V`, with `r,s<=m`, and no edges between `U,V`.

For such a graph, `m-1` is a **simple** second adjacency eigenvalue. Indeed, deleting `w` leaves `2K_m` with two eigenvalues `m-1`, so interlacing gives `lambda_2(A(G))>=m-1` and `lambda_3(A(G))<=m-1`; the universal adjacency bound gives `lambda_2<=m-1`. To prove simplicity directly, if `(x,t)` is an eigenvector for `m-1`, where `t=x_w`, the equations at each vertex of a clique imply `m x_i = S_U + t 1_{i~w}` on `U`. Summing over `U` yields `r t=0`, so `t=0`; likewise `s t=0`. Thus `x` is constant on each clique, and the `w`-equation is `r x_U+s x_V=0`. The eigenspace is exactly one-dimensional, spanned by `x=s` on `U`, `x=-r` on `V`, `x_w=0`.

For a differentiable simple symmetric-matrix eigenvalue the derivative at `a=0` is its normalized Rayleigh quotient of `D-A`. Since `x^T(D-A)x=sum_{ij in E}(x_i-x_j)^2`, only edges incident to `w` contribute. Hence

`(d/da)lambda_2(A_a(G))|_{a=0} = [r*s^2 + s*r^2]/[m*(s^2+r^2)] = rs(r+s)/(m(r^2+s^2)) <= 1`.

The two inequalities `r+s<=2m` and `2rs<=r^2+s^2` show equality in the derivative bound only if `r=s=m`, precisely `H_m`. For any other connected adjacency maximizer, the derivative is **strictly less than 1**, so `lambda_2(A_a(G))<m-1+a` for all sufficiently small positive `a` (the interval may depend on `G`). If a graph is not an adjacency maximizer, its strict `a=0` gap persists by continuity. Disconnected graphs satisfy `lambda_2(A_a(G))<=m-1<m-1+a` for all `0<a<1/2` by R4. Finally, there are **finitely many** graphs of a given order; take the minimum of the finitely many positive intervals. The H_m equality `lambda_2=m-1+a` for `a<1` is the already-published Xue--Liu--Yu--Shu construction. This proves the theorem, including uniqueness, without estimating `epsilon_m`.

**Priority:** This is a deduction from the June 2026 arXiv preprint's adjacency equality classification and the published 2020 H_m spectrum. We found no prior statement of this *A_alpha near-zero deduction* in the bounded search, but do not assert first discovery.

### R2. Continuous-alpha exact exhaustive theorem at n=3,5,7 (EXACT CERTIFICATE; two methods)

For **every** unlabeled graph of orders `3,5,7` and **all real** `0<=a<1/2`, the proposed upper bound holds. At `a>0`, `H_m` is the only equality class. At `a=0`, equality classes number `3,6,11` respectively. Thus this is strictly stronger in alpha coverage (but not in graph order) than `oddlambda2/`'s 13-point scan for `n=5,7,9`.

`interval_scan.py` computes `p_G(a)=det(A_a(G)-(m-1+a)I)` with **integer** Bareiss determinants at the integer nodes `0,...,n`, reconstructs the degree-at-most-`n` polynomial by exact forward differences, strips roots at `0` and `1/2`, and isolates every interior algebraic root with exact rational intervals. The output has **no multiple interior roots**. A sample from each root-free interval, and both sides of every nonrational isolating interval, is checked by exact positive inertia via the sign variations of the characteristic polynomial (Descartes is exact since all roots of the symmetric matrix are real). On every sample, positive inertia is at most one; at `a=0`, the same exact test counts equality.

Why finite samples cover the **continuum**: the inertia is constant on each connected component of the complement of the roots of `p_G` (no eigenvalue crosses zero). At a **simple** interior root of `p_G`, `lambda_2` cannot just touch zero while retaining nonpositive values on both sides. Indeed, if the shifted matrix has two or more zero eigenvalues, the determinant has at least a double zero as a function of `a` (its adjugate vanishes); if its kernel is one-dimensional, the simple determinant root forces its sole zero eigenvalue to change sign, which would give two positive eigenvalues on one adjacent side if `lambda_2=0` and `lambda_1>0`. Since `lambda_1>=lambda_2`, a simple root with `lambda_1=lambda_2=0` is already excluded by corank at least two. Thus all isolated interior zeros are **strictly below** the second-eigenvalue threshold. The only graph with `p_G` identically zero is isomorphic to H_m, whose exact spectrum is known. This completes the finite-order continuous-parameter certificate.

`interval_checker.py` independently recomputes the polynomial by a direct SymPy symbolic determinant (not Bareiss/interpolation), independently validates root counts using Sturm, and counts inertia by **rational symmetric congruence and Schur complements** (not Descartes). Outputs agree: `n=7` has 1,044 graph classes, 38 graphs with 40 simple interior roots and 1,151 exact rational positive-inertia tests; 0 failures. See both committed logs and hashes. **Shared dependency limitation:** both use NetworkX's graph atlas for completeness up to seven vertices and SymPy for real-root isolation; the matrix methods differ. The independent `geng` atlas generation in `guide/2026-11-15-independent-checks/` supports the graph counts through order 9 but does not replay this interval isolation.

### R3. All alpha=0 equality graphs at arbitrary odd n (ALREADY PUBLISHED, now located)

Lin, Liu, Yin and You, *Krahn--Szego type inequalities for graphs*, arXiv:2606.11659 (submitted 2026-06-10), Theorem 1.4, gives a general-class adjacency nodal-domain result. Apply it to the class of all graphs, where `max_{|V(H)|=s} rho(A(H))=s-1`, uniquely at `K_s`, and this maximum is strictly increasing. The resulting **complete list**, up to isomorphism, of graphs on `2m+1` vertices with `lambda_2(A)=m-1` is:

1. **Disconnected:** `K_m union F`, where `F` is *any* graph on `m+1` vertices with `rho(A(F))>=m-1`.
2. **Connected:** take `K_m union K_m` and a new vertex `w` joined to an arbitrary **nonempty** subset of each clique, with no edges between the cliques.

Both types indeed attain equality: the first by direct union spectra, the second by Cauchy interlacing of the principal `2K_m` and the Powers/Hong adjacency bound. The necessity is the **published 2026 classification**, not proved independently in this directory. The disconnected case includes `K_3 union C_4` at `n=7`, which the narrower `2K_m`-plus-vertex construction misses. The `n=3,5,7` equality counts `3,6,11` are checked exactly in R2; `oddlambda2/` reports 17 at `n=9`, consistent with the classification.

### R4. Useful all-order reductions (PROVED, elementary plus published bounds)

For any graph `G` on `n=2m+1` vertices and `0<=a<1/2`:

* **Disconnected:** `lambda_2(A_a(G))<=m-1`. To see this, write component sizes `s_1>=s_2>=...`. Each component's second eigenvalue is at most `s_i/2-1<=m-1` by Nikiforov's published Proposition 31(b) since `s_i<=2m`. A component's spectral radius is at most its maximum degree, hence `<=s_i-1`, and the second largest of the component top eigenvalues is at most `s_2-1<=m-1`. These account for the global second eigenvalue. Thus disconnected graphs are strictly nonoptimal for every `a>0`.
* **Maximum degree at most m:** Weyl's inequality and the published adjacency bound `lambda_2(A(G))<=m-1` give `lambda_2(A_a(G)) <= (1-a)(m-1)+a Delta(G)<=m-1+a`. If `Delta(G)<m` and `a>0` the inequality is strict. Only `Delta>=m` can obstruct this route.
* **Regular graphs of any degree:** since `A` and `Q=D+A` commute, `A_a=(1-2a)A+aQ` gives `lambda_2(A_a)=(1-2a)lambda_2(A)+a q_2(Q)`. Using the published bounds `lambda_2(A)<=m-1` and `q_2(Q)<=n-2=2m-1`, the result is `<=m-1+a`. Here the eigenvalue ordering is preserved because the eigenvalue of A_a corresponding to adjacency eigenvalue `x` is `a d+(1-a)x`; the displayed combination is the same. These are reductions, not a proof for all irregular connected high-degree graphs.

### R5. Alpha >= 1/2 (COMPLETE PROOF; ALREADY PUBLISHED in substance)

For `1/2<=a<=1` the **global** maximum is `a n-1`, attained by `K_n`; it exceeds the H_m value for `a>1/2` by `m(2a-1)`. Equality for `a>1/2` holds exactly for graphs with **at least two universal vertices** (Lin--Xue--Shu 2018, Theorem 1.1). At `a=1/2` equality holds if and only if the complement has **a balanced bipartite component or at least two bipartite components** (de Lima--Nikiforov 2013). An isolated vertex counts as a bipartite component.

**Self-contained matrix proof of these characterizations:** let `J` denote the all-ones matrix. Then `A_a(K_n)=(an-1)I+(1-a)J`, and

`B=A_a(K_n)-A_a(G)= (1-a)Q(complement G)+(2a-1)D(complement G)`, which is positive semidefinite when `a>=1/2`. For any `x` perpendicular to the all-ones vector, `x^T A_a(G)x=(an-1)||x||^2-x^T Bx`. Courant--Fischer yields `lambda_2<=an-1`, with equality **exactly when** `ker B` contains a nonzero vector perpendicular to 1 (if equality, the max over this subspace must reach the bound; conversely such a vector is an `an-1` eigenvector and the bound forces `lambda_2=an-1`, because this eigenvector changes sign (it is orthogonal to 1), whereas a unique top eigenvector of the entrywise nonnegative symmetric matrix `A_a(G)` may be taken nonnegative by Perron--Frobenius; hence either an eigenvalue lies above it or the top eigenvalue is repeated). For `a>1/2`, `ker B` consists of the vertices isolated in the complement, i.e. universal vertices in G. A nonzero zero-sum vector in that space exists iff at least two. For `a=1/2`, `B=Q(complement G)/2` and its kernel has one alternating-sign vector for each bipartite connected component. Such a vector is orthogonal to 1 iff its component is balanced; with two or more bipartite components a nontrivial zero-sum combination exists. This also proves the equality characterization directly. This theory is **not claimed as a new result**.

### R6. An obstruction to a tempting proof shortcut (EXACT WITNESS)

The unjustified interpolation inequality

`lambda_2((A+Q/2)/2) <= (lambda_2(A)+lambda_2(Q/2))/2`

is **false**, already for the five-vertex graph `Dlc` (graph6; edges `01,03,04,12,23,34`). `obstruction.py` isolates the second eigenvalues of `A_0,A_{1/4},A_{1/2}` by **exact rational intervals** of their characteristic polynomials and proves

`lambda_2(A_{1/4}) > 1.03214 > (0.68890+1.37296)/2 > (lambda_2(A_0)+lambda_2(A_{1/2}))/2`.

Here `A_{1/4}=(A_0+A_{1/2})/2`. Thus simply convex-combining the sharp adjacency and Q/2 **second-eigenvalue** upper bounds is not a valid argument. This does **not** refute the proposed odd-order extremal bound itself (this graph remains below `m-1+a=1.25`).

## Replay and integrity

From `oddl2proof/`, with free Python 3, NetworkX and SymPy:

```sh
sh run.sh
```

The script runs the **two different exact continuum methods** and the exact five-vertex obstruction, byte-compares the logs to committed copies and checks the SHA-256 checksums. The committed logs include versions, count breakdowns and graph6 identifiers. `SHA256SUMS` covers every other file in this directory (programs, report, run script and logs); it excludes itself. Pushed content is fetched back from the branch via the GitHub connector and checked against its local SHA-256.

## Priority check

Checked **before novelty attribution**: main-branch `README.md`; `verify/REPORT.md`, `attack/REPORT.md`, `status/REPORT.md`, `priority/REPORT.md`, `classify22/REPORT.md`, `problem30odd/REPORT.md`, `oddlambda2/REPORT.md`, `forkmatrices/REPORT.md`, and `guide/2026-11-15-independent-checks/README.md`. In particular `oddlambda2`'s H_m polynomial and finite grids are not claimed as new. External sources searched/read on 2026-11-15:

* Nikiforov, *Merging the A- and Q-spectral theories*, journal (2017), Proposition 31, https://doi.org/10.2298/AADM1701081N; arXiv:1607.03015. The journal claim of uniqueness for `a>1/2` is contradicted by the 2018 published classification.
* Xue--Liu--Yu--Shu, *The Multiplicity of A_alpha-Eigenvalues of Graphs*, ELA 36 (2020), 645-657, equation (3.7), https://doi.org/10.13001/ela.2020.4999. Published H_m factorization.
* Lin--Xue--Shu, *On the A_alpha-spectra of graphs*, LAA 556 (2018), Theorem 1.1, https://doi.org/10.1016/j.laa.2018.07.003; arXiv:1709.00182. Published high-alpha equality classification.
* de Lima--Nikiforov, *On the second largest eigenvalue of the signless Laplacian*, LAA 438 (2013), 1215-1222, Theorem 1, https://doi.org/10.1016/j.laa.2012.07.052; arXiv:1202.0964. Published Q equality classification.
* Lin--Liu--Yin--You, *Krahn--Szego type inequalities for graphs*, submitted **2026-06-10**, arXiv:2606.11659, Theorem 1.4. Its freely viewable full text was checked at https://www.researchgate.net/publication/406910462_Krahn-SzegH_o_type_inequalities_for_graphs , including the definition of `G_1 bullet K_1 bullet G_2` and the exact equality statement. This supersedes claims that the `a=0` classification remains open.
* Chen--Li--Meng (2019), https://doi.org/10.1016/j.laa.2019.06.027 (different minimization direction); Liu--Das--Shu (2020), https://doi.org/10.1016/j.disc.2020.111917; the existing repository's previously documented bounded citation checks.

Targeted web searches included `"second largest" "A_alpha" "odd" "maximum" 2026 2025 2024`, `"A_alpha" "K_1" "K_m" "second largest" maximum`, `Powers odd order second largest adjacency eigenvalue equality`, `Krahn Szego type inequalities adjacency extremizers Theorem 1.4 2026`, and `de Lima Nikiforov second signless Laplacian equality complement balanced bipartite`. A complete citing-paper harvest was **not** performed. The original June 2026 paper is a **preprint**, not a verified peer-reviewed publication. No statement of novelty is categorical.

## Not covered

* **No proof or counterexample for all real `a` in `(0,1/2)` at arbitrary odd `n`**, or for arbitrary `a` at `n>=9`. The positive near-zero interval `epsilon_m` is nonconstructive and may shrink with `m`.
* The exact continuous scan is only `n=3,5,7`; the inherited `n=9` certificate is limited to 13 rational alpha points and the guide's 7-point numerical cross-check. No sweep at `n>=11`.
* The complete `a=0` **necessity** rests on the 2026 preprint Theorem 1.4, not an independent proof in this directory. The first-order perturbation proof is fully written above, given that input.
* NetworkX atlas and SymPy root isolator are shared by the two continuum programs; exact determinant and inertia engines are independent. No separate all-graphs generator was run for this task (the guide already used `geng` for n<=9).
* The five-vertex witness refutes the naive spectral interpolation inequality, **not the odd-order conjecture**; no obstruction to the conjecture itself was found.
* No exhaustive current literature priority census, no paid service.

## Next steps

1. **Highest priority:** adapt the 2026 adjacency nodal-domain classification (especially the two `K_m` nodal blocks) to the diagonal-potential matrix `A_a = A+a(D-A)`. Seek an a-dependent two-domain bound and a rigidity statement forcing the central vertex to meet all `2m` clique vertices. Try small perturbations away from the classified `a=0` equality graphs; their exact slopes `rs(r+s)/[m(r^2+s^2)]` quantify the first-order gap.
2. Derive an **explicit** lower bound on `epsilon_m` from the adjacency spectral gap and the slopes; or a uniform gap to settle all `a<1/2`.
3. Extend continuous-root certificates from n=7 to n=9 using isomorphism decomposition and exact resultants rather than 274,668 full symbolic determinants; separately run a second n9 generator and a bounded high-alpha structural scan.
