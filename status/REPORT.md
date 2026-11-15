# A-alpha spectral questions: state of the art and reproducible ideas

**Research cutoff:** 2026-11-15. **Scope:** finite simple undirected graphs; A_alpha(G) = alpha D(G) + (1-alpha) A(G), alpha in [0,1]. The original source is V. Nikiforov, *Merging the A- and Q-spectral theories*, *Applicable Analysis and Discrete Mathematics* **11** (2017), 81–107, [DOI](https://doi.org/10.2298/AADM1701081N), [arXiv:1607.03015, original PDF](https://arxiv.org/pdf/1607.03015).

## Scope

Read the primary arXiv paper, including the exact locations below. Identified **Question 5; Problems 7, 22, 24, and 30; the general Problem D; and an unnumbered odd-order second-eigenvalue question after Proposition 31**. Problems A/B/C (arXiv p. 15) are older motivating adjacency/signless-Laplacian/extremal questions rather than new A-alpha conjectures. The observation preceding Problem 22 is an *informal conjectural equality claim*. There are no statements formally labeled "Conjecture" in the original paper.

**Source-transcription limitation:** The complete original prose cannot be reproduced verbatim here; the mathematical content below is a faithful transcription/paraphrase, with problem numbers, equation numbers, and page-specific primary-source links allowing verbatim comparison. Mathematical formulas and numbered references are retained. Page numbers in this report refer to **printed pages of the arXiv PDF** (not necessarily publisher pagination).

**Status semantics:** "solved" requires a complete written argument or a cited theorem that covers the stated domain; "partly solved" identifies exactly the covered parameters; "not verified solved/open under this search" means no general solution was located, not a proof of absence. Direct mathematical arguments here are labeled as such, not spuriously assigned to external authors.

## Findings

### Question 5 — differentiability (arXiv p. 7, §3.2; [page](https://arxiv.org/pdf/1607.03015#page=7))

Original question asks whether, for each n >= k >= 1, f(alpha) = lambda_k(A_alpha(G)) is differentiable in alpha. **Status: NO in general; complete elementary counterexample (independently derived here).**

Take G = K_2 disjoint union K_1. The symmetric block for K_2 has eigenvalues 1 and 2alpha-1; the isolated vertex supplies 0. Therefore
`lambda_2(A_alpha(G)) = max(0, 2alpha-1)`.
The left and right derivatives at alpha=1/2 are 0 and 2. This explicitly refutes universal differentiability, and the symbolic characteristic polynomial is checked in `exact_checker.py`. In contrast, for a *connected* G and 0 <= alpha < 1 the Perron root is simple (Perron–Frobenius), hence locally analytic in alpha; this does not imply differentiability of other ordered eigenvalues.

### Problem 7 — positive-semidefinite threshold (arXiv p. 7, §3.3; [page](https://arxiv.org/pdf/1607.03015#page=7))

Original task: determine the minimum alpha for which A_alpha(G) is positive semidefinite. **Status: SOLVED as a mathematical problem; self-contained general formula below, plus 2017 published partial structural results.**

For an edgeless graph alpha_0(G)=0. Otherwise discard isolated vertices, let D_+ be the diagonal of strictly positive degrees, put B=D_+^(-1/2) A_+ D_+^(-1/2), and let mu=lambda_min(B). Then
`alpha_0(G) = -mu/(1-mu)`.

**Complete proof.** On the positive-degree block,
`A_alpha = D_+^(1/2) [alpha I+(1-alpha)B] D_+^(1/2)`.
Since D_+^(1/2) is invertible, Sylvester's law of inertia says the two bracketed matrices are simultaneously positive semidefinite. Their eigenvalues are alpha+(1-alpha)mu_i; the smallest is alpha+(1-alpha)mu. This is >=0 exactly when alpha>=-mu/(1-mu). For a nonempty graph the normalized adjacency has eigenvalue 1 and trace zero, so mu<0; the formula is well-defined. Isolated vertices add zero eigenvalues only. This proves the threshold for arbitrary finite G, disconnected included.

**Published evidence:** V. Nikiforov and O. Rojo, *A note on the positive semidefiniteness of A_alpha(G)*, *Linear Algebra Appl.* **519** (2017), 156–163, [DOI](https://doi.org/10.1016/j.laa.2016.12.042), abstract items (1)–(3): d-regular formula, equivalence of alpha_0=1/2 with having a bipartite component, and r-colorable lower bound alpha_0>=1/r. The general normalized-adjacency formula above is explicitly a self-contained rederivation, **not** a claim that this exact formula appears in that article.

### Problem 22 — equality in inequalities (20)–(23) (arXiv pp. 13–14, §4.3; [page](https://arxiv.org/pdf/1607.03015#page=14))

Original asks for all equality cases in the two average-neighbor-degree bounds (20),(21) and two edge-degree bounds (22),(23) for rho_alpha(G). The neighboring informal conjecture says that for (20),(21), away from alpha=0,1/2, equality should force regularity. **Status: PARTLY SOLVED; important connected cases are known, disconnected cases not fully classified in cited results.**

**Published evidence:** S. Liu, K. Ch. Das, and J. Shu, *On the eigenvalues of A_alpha-matrix of graphs*, *Discrete Mathematics* **343** (2020), 111917, [DOI](https://doi.org/10.1016/j.disc.2020.111917), **§5, Theorems 5.4–5.6 and Remark 5.7**. Theorem 5.4 characterizes equality in original (22),(23) for connected G and 0<=alpha<1: regular graphs, or semiregular bipartite graphs at alpha=1/2. Theorems 5.5 and 5.6 give necessary conditions for **disconnected** G in edge-degree and average-neighbor-degree inequalities respectively; Remark 5.7 shows their converses need not hold. The original Proposition 20 already characterizes connected equality for (20),(21) when alpha in (1/2,1). The complete combined equality table (all alpha, connected/disconnected, isolated-vertex conventions) was **not** verified from the literature.

### Problem D — arbitrary forbidden F (arXiv p. 15, §5; [page](https://arxiv.org/pdf/1607.03015#page=15))

Original generic extremal objective: maximize rho_alpha(G) among n-vertex graphs not containing a fixed graph F as a subgraph. **Status: OPEN as a universal problem; many F-families solved.** In particular the original paper's Theorem 27 handles the complete-graph case, and the path/cycle special case is Problem 24 below. A classification for *every* F is not furnished by these family-specific theorems. No single universal solution was verified.

### Problem 24 — forbidding a path or cycle of prescribed order (arXiv p. 15, §5; [page](https://arxiv.org/pdf/1607.03015#page=15))

Original special case of Problem D: F is P_k or C_k. **Status: PARTLY SOLVED (much of the large-n range).** The problem does not restrict n to be sufficiently large or alpha to be an interior value.

- G.-X. Tian, Y.-X. Chen, S.-Y. Cui, *The extremal alpha-index of graphs with no 4-cycle and 5-cycle*, *Linear Algebra Appl.* (2021), [DOI](https://doi.org/10.1016/j.laa.2021.02.022): C_4/C_5 cases, with hypotheses in that paper (do not infer an all-n/all-alpha solution solely from the abstract).
- M.-Z. Chen, A.-M. Liu, X.-D. Zhang, *On the A_alpha-spectral radius of graphs without linear forests*, *Applied Mathematics and Computation* **450** (2023), 128005, [arXiv:2304.03046](https://arxiv.org/pdf/2304.03046#page=2), [DOI](https://doi.org/10.1016/j.amc.2023.128005), **Theorem 1.1**, three cases: a linear forest containing an even-order path, one containing a path of odd order >3, and a disjoint union of P_3's. This yields every single forbidden P_k of order k>=3 for **fixed 0<alpha<1 and n sufficiently large** (P_2 is trivial separately).
- S. Li and Y. Yu, *On A_alpha spectral extrema of graphs forbidding even cycles*, *Linear Algebra Appl.* **668** (2023), 11–27, [DOI](https://doi.org/10.1016/j.laa.2023.03.017), main theorems (cf. article introduction): C_{2k} for fixed k and alpha, sufficiently large n.
- P. Kuang, F. Liu, S. Sun, Y. Wang, J. Zeng, *Tight Bound for Nikiforov's Spectral Even-Cycle Conjecture*, **arXiv preprint posted 2026-06-05**, [arXiv:2606.06987](https://arxiv.org/abs/2606.06987), **Theorem 1.2**: for each epsilon>0, k>=k_epsilon, 0<=alpha<=1-epsilon, and n>=C_epsilon*k, C_{2k+2}-free graphs maximize at S^+_{n,k} uniquely. This improves the sufficient threshold to linear in k, uniformly away from alpha=1. **Preprint, not independently peer-reviewed in this audit.**

Not covered by these cited results: a theorem for *all* n and alpha for arbitrary P_k/C_k, and a full fixed-parameter C_{2k+1} (odd k>=3) classification. Other papers may exist; no exhaustive negative literature claim is made. Note that alpha=1 is degenerate (rho_1=maximum degree) and requires separate equality analysis.

### Problem 30 — smallest eigenvalue at 0<alpha<1/2 (arXiv p. 21, §6.1; [page](https://arxiv.org/pdf/1607.03015#page=21))

Original optimization: `min_{|V(G)|=n} lambda_min(A_alpha(G))`. **Status: PARTLY SOLVED by complete independent argument for EVEN n; ODD n not verified solved in literature search.**

**Theorem for even n (full proof).** If n=2m and 0<=alpha<1/2, then
`min_{|V(G)|=2m} lambda_min(A_alpha(G)) = -m(1-2alpha)`,
attained by K_{m,m}.

To prove the lower bound without external lemmas, let x be any unit real vector and partition its indices into P={i:x_i>0}, N={i:x_i<0}, and zeros. Contributions from edges with same-sign endpoints to x^T A(G)x are nonnegative. Consequently
`x^T A(G)x >= -2 (sum_{i in P} x_i)(sum_{j in N}|x_j|)`.
By Cauchy–Schwarz and 2uv<=u^2+v^2, this is at least
`-sqrt(|P||N|) (||x_P||^2+||x_N||^2) >= -sqrt(floor(n/2)*ceil(n/2))`.
Thus lambda_min(A(G)) >= -sqrt(floor(n/2)*ceil(n/2)). Next Q(G)=D+A is positive semidefinite since `x^T Qx=sum_{uv in E}(x_u+x_v)^2>=0`, and
`A_alpha = alpha Q+(1-2alpha)A`.
The variational principle gives lambda_min(A_alpha)>=-(1-2alpha)sqrt(floor(n/2)ceil(n/2)). For n=2m, take K_{m,m} and x constant +1 on one side, -1 on the other. Then Qx=0 and Ax=-mx, so A_alpha x= -m(1-2alpha)x; equality holds. QED.

For odd n=2m+1, the bound just proved is in general **not attained**, so it does not solve the original problem. For K_{m,m+1}, the smallest eigenvalue is
`(alpha*n - sqrt(alpha^2*n^2 + 4*(1-2alpha)*m*(m+1)))/2`,
as follows from its 2x2 equitable quotient. This is a candidate for a sharp odd-order bound, **not a proved global minimum**. The atlas scan checks n=3,5,7 only at selected alpha.

**Relevant but NOT a solution to Problem 30:** H. Lin, J. Xue, J. Shu, *On the A_alpha-spectra of graphs*, *Linear Algebra Appl.* **556** (2018), [arXiv:1709.00182](https://arxiv.org/abs/1709.00182), concentrates on alpha>1/2. Do not mistake a theorem about the positive-semidefinite regime for this negative-eigenvalue optimization.

### Unnumbered open question — maximal lambda_2 for odd n (arXiv p. 22, after Proposition 31; [page](https://arxiv.org/pdf/1607.03015#page=22))

Nikiforov obtains `lambda_2(A_alpha(G)) <= n/2-1` for 0<=alpha<1/2, with equality for two equal cliques when n is even; for odd n the disjoint union of unequal cliques reaches only (n-3)/2. **Status: OPEN/UNVERIFIED in general; new finite-n candidate and complete lower-bound witness.** Published subsequent work by Y. Chen, D. Li, J. Meng, *On the second largest A_alpha-eigenvalues of graphs*, *Linear Algebra Appl.* **580** (2019), 343–358, [DOI](https://doi.org/10.1016/j.laa.2019.06.027), concerns a lower bound and connected **minimizers**, not this all-graph maximum.

**Exact lower-bound construction and proof.** For n=2m+1, let H_m consist of two K_{m+1}'s intersecting in exactly one vertex v (n vertices). The degree of v is 2m, the other 2m vertices have degree m. Assign x_v=0, x=+1 on one clique's m noncentral vertices, and x=-1 on the other m. Direct multiplication gives
`A_alpha(H_m)x=(m-1+alpha)x`.
This eigenvalue is exactly lambda_2: variations summing to zero within one noncentral m-set have eigenvalue (m+1)alpha-1, less than m-1+alpha for alpha<1. The two remaining eigenvalues are those of the quotient
`[[2m*alpha, 2m*(1-alpha)], [1-alpha, m-1+alpha]]`;
its positive off-diagonal product implies one root strictly above and one strictly below m-1+alpha when alpha<1. Hence
`lambda_2(A_alpha(H_m)) = m-1+alpha = (n-3)/2+alpha`
for 0<=alpha<1. This closes **alpha of the original 1/2 gap** as a lower-bound improvement at odd n, but **does not prove optimality**.

**Conjectural target:** for odd n>=3 and 0<=alpha<=1/2, `max_G lambda_2(A_alpha(G)) = (n-3)/2+alpha`. At alpha=0 or 1/2 extremizers need not be unique. The exhaustive numerical scan of all unlabeled graphs at n=3,5,7 and alpha=0,.25,.49,.5,.75 matches this expression for alpha<=.5. This is **a computational conjecture, not a result in arbitrary n or continuous alpha**. The published Proposition 31(a) takes over for alpha>1/2.

## Evidence

### Replay instructions

Run with Python 3.13.5, NumPy 2.3.5, NetworkX 3.6.1, SymPy 1.14.0 (the code prints library versions, which can differ on replay). No service calls or external paid APIs are made:

```sh
python status/atlas_scan.py > status/atlas_scan.replay.log
python status/exact_checker.py > status/exact_checker.replay.log
diff -u status/atlas_scan.log status/atlas_scan.replay.log
diff -u status/exact_checker.log status/exact_checker.replay.log
sha256sum status/atlas_scan.py status/atlas_scan.log status/exact_checker.py status/exact_checker.log
sha256sum -c status/SHA256SUMS
```

**Enumeration scope:** NetworkX `graph_atlas_g()`, all 1,253 unlabeled graphs of orders 0–7, filtered to 4+34+1,044=**1,082 graphs of n=3,5,7**. Scan uses `numpy.linalg.eigvalsh`, tolerance 1e-9, sampled alpha {0,1/4,49/100,1/2,3/4}. Log names a maximizer/minimizer via the NetworkX atlas indices and explicitly prints witness edges for lambda_2 at alpha=1/4. **This is numerical screening, not an interval-certified enumeration.**

**Independent method:** `exact_checker.py` does not use NumPy or the NetworkX eigenvalue routines. It constructs adjacency/degree matrices directly in SymPy, verifies the Question 5 characteristic polynomial symbolically, checks the signed vector of H_m by exact polynomial matrix arithmetic, verifies its second eigenvalue for n=5,7 at four rational alpha values via SymPy eigenvalues, and verifies the bipartite equitable-quotient characteristic polynomial. A complete written argument for the H_m lambda_2 result and for even-order Problem 30 appears above. SymPy and NumPy are independent eigenvalue implementations; the atlas *enumeration source* itself is not independently regenerated.

**Observed sampled maxima for lambda_2 (alpha<=1/2):**

| n | graphs | alpha=0 | alpha=1/4 | alpha=49/100 | alpha=1/2 |
|---:|---:|---:|---:|---:|---:|
| 3 | 4 | 0 | 0.25 | 0.49 | 0.5 |
| 5 | 34 | 1 | 1.25 | 1.49 | 1.5 |
| 7 | 1044 | 2 | 2.25 | 2.49 | 2.5 |

**SHA-256** (also in `status/SHA256SUMS`; calculated locally on actual program and log bytes):

```text
a80798fe5acc318c31c723a3582b523273736307b24a95218f848d50c743a299  status/atlas_scan.py
b0de0a9578468e9b2ae2b3aef58c37e48b646531e9b16ffebc12f492e9c01184  status/atlas_scan.log
46a432f846f4f9755c282d4cea40c00815937866717c1343a5f87b54ae3d0cf7  status/exact_checker.py
a0ba2d1e2981f186b1dd77e2c6501fe38c40e6f2ce057646119d4728eecda489  status/exact_checker.log
```

### Search trace and provenance

Primary paper: [arXiv:1607.03015](https://arxiv.org/abs/1607.03015), [journal](https://doi.org/10.2298/AADM1701081N). Later work checked against theorem/section statements at publisher, journal/DOI and arXiv pages linked above (2017 Nikiforov–Rojo; 2018 Lin–Xue–Shu; 2019 Chen–Li–Meng; 2020 Liu–Das–Shu; 2021 Tian–Chen–Cui; 2023 Chen–Liu–Zhang and Li–Yu; 2026 Kuang et al.). Public web searches also targeted arXiv, literature citation trails, MathOverflow, Mathematics Stack Exchange, and OEIS for specific problem numbers/formulas; no direct public MO/MSE/OEIS resolution was verified. Requests to the OpenAlex and Semantic Scholar API endpoints were inaccessible in this browsing environment, so **a systematic cited-by API traversal was not accomplished**. No paid services were used. Absence of a verified hit is not evidence that no publication exists.

## Not covered

1. **Full verbatim reproduction** of the source's longer original problem statements: use the linked original page to compare exact language. This file preserves their mathematical substance and exact locations.
2. **Complete citation-graph audit**, particularly OpenAlex and Semantic Scholar API traversal (endpoints inaccessible); unindexed/2026 literature may be missed.
3. **Proof** of the odd-n lambda_2 conjecture, the odd-n case of Problem 30, the entire disconnected Problem 22 equality classification, or general finite-n forbidden odd-cycle cases.
4. **All real alpha**, all n, and independently generated nonisomorphic graph lists: only n=3,5,7 and five alpha values were enumerated numerically; a second method checks witnesses, not atlas completeness.
5. **Machine-checked formal proof** or interval-arithmetic certification. Proofs for Question 5, Problem 7, even-n Problem 30, and the H_m witness are written in full above; cited literature results were not re-proved from first principles.
6. **Publication priority claims** about the elementary proofs or proposed conjecture; they are presented as independent derivations requiring a fresh novelty audit.

## Next steps

1. **Exact phase diagrams for odd-order lambda_2.** First step: enumerate n<=8 using a separately generated graph list (e.g., nauty), generate characteristic polynomials over Q[alpha], and certify any changes of maximizer by Sturm sequences/isolating intervals; compare hashes against atlas for n<=7. Deliver a proof certificate and independent checker.
2. **Structural upper bound for the odd-order conjecture.** First step: prove or disprove `lambda_2<= (n-3)/2+alpha` for graphs with a cut vertex using interlacing on the two principal blocks; test whether coalesced equal cliques are extremal within that subclass. Record each lemma in a complete written proof, and publish counterexample edge lists if false.
3. **Problem 30 at odd n.** First step: extend the exhaustive search to n=8, then n=9 selectively, at rational alpha near 0 and 1/2. Test K_{floor(n/2),ceil(n/2)} against nonbipartite graph families; use exact rational principal-minor/Sturm checking to eliminate numerical ambiguity.
4. **Finish Problem 22's disconnected cases.** First step: derive a necessary-and-sufficient componentwise criterion in terms of each component's rho_alpha and its local min/max edge/neighbor-degree bound, and use the 2020 §5 constructions to see which structural simplifications fail; then write the complete equality proof.
5. **Formalization and sharp endpoint bookkeeping.** First step: mechanize the normalized-adjacency inertia reduction for Problem 7 and the even-n least-eigenvalue proof for Problem 30 in Lean (or independently check every inequality on paper); track alpha=0,1/2,1 and isolated vertices as separate cases.

Each proposed contribution is conditioned on the lab standard: code + exact output log + SHA-256 of both, an independent validator or separately checkable witness, and complete written proofs for any mathematical claim.
