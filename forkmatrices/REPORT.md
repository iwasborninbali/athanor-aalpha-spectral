# Forking Problem 22: degree- and transmission-scaled Perron bounds

**Athanor lab (Aleksei Kudriashov); investigation date 2026-11-15.** Scope: finite simple undirected connected graphs unless specified; nonnegative irreducible graph matrices; `0 <= alpha < 1` unless stated. All newly supplied files are under `forkmatrices/`. Status labels distinguish a **proved mathematical statement** from **publication priority** and a **finite observation**. No paid services were used.

## Question

For a nonnegative irreducible matrix `M` and a positive graph invariant vector `v`, the scaled matrix `diag(v)^(-1) M diag(v)` has row sums `r_i=(Mv)_i/v_i`. Collatz--Wielandt gives `min_i r_i <= rho(M) <= max_i r_i`. **Equality at either end holds exactly when all `r_i` agree** (proof below). For the adjacency-based bounds of Nikiforov's Problem 22, `v` is the degree vector. For distance-based analogues `v` is the transmission vector. Which assertions that this requires regularity/transmission regularity survive the change of matrix?

Repository baseline read before working: [`verify/REPORT.md`](../verify/REPORT.md), [`attack/REPORT.md`](../attack/REPORT.md), [`status/REPORT.md`](../status/REPORT.md), and [`priority/REPORT.md`](../priority/REPORT.md) on `main`. In particular `attack/` **already proves** the undirected `A_alpha` equality family `K_t join I_m`, `alpha=m/(2m+t-1)`, and `verify/` supplies `K_1 join K_{2,2}` at `alpha=1/4`. They are **not** claimed as discoveries here. The 2020 Liu--Das--Shu equality results and 2017 preprint-versus-journal distinctions are recorded in `status/` and `priority/`.

### General certificate lemma (complete proof)

Let `M>=0` be irreducible and `v>0`. Put `a=min_i(Mv)_i/v_i` and `b=max_i(Mv)_i/v_i`. Then `av<=Mv<=bv`; multiplication by a strictly positive left Perron eigenvector `y` gives `a<=rho(M)<=b`. If equality holds at either endpoint, the corresponding nonnegative vector `Mv-av` or `bv-Mv` has zero inner product with `y>0`, so it must vanish. Conversely if `Mv=lambda*v`, positivity and Perron--Frobenius force `lambda=rho(M)` and both ends coincide. Thus **a positive exact eigenvector is a two-sided certificate**, not merely numerical approximate equality. This applies to symmetric distance and degree matrices and to irreducible directed `A_alpha` for `alpha<1`; not automatically to disconnected or endpoint `alpha=1` matrices.

## Results with evidence status

### 1. Distance matrix and generalized distance `D_alpha`: non-transmission-regular infinite family

**Status: PROVED, independently checked with rational arithmetic; specific-family publication priority NOT VERIFIED.** Let `G_{t,m}=K_t join I_m`, `t>=1`, `m>=2`; it is connected and transmission-irregular. Write `E` for its graph-distance matrix, `tau_i=sum_j E_ij`, and `M_alpha=alpha*diag(tau)+(1-alpha)*E`. Clique vertices have transmission `h=t+m-1`; independent-set vertices have transmission `ell=t+2m-2`. The distance-weighted transmission sums for the two vertex classes are

```
U_H=(t-1)*h+m*ell,
U_L=t*h+2*(m-1)*ell.
```

The scaled row sums are `R_H=alpha*h+(1-alpha)*U_H/h` and `R_L=alpha*ell+(1-alpha)*U_L/ell`. Direct substitution and clearing the positive denominators show `R_H=R_L` exactly when

```
alpha(t,m) = (2*m + 2*t - 2 - m*t)
             / (2*m*m + 2*m*t - 2*m + t*t - t).
```

For each pair with nonnegative numerator, this alpha is in `[0,1/2)`. The vector `tau` has entries `h` and `ell` (distinct because `m>=2`), and `M_alpha*tau=R_H*tau`; the certificate lemma proves **equality in both transmission-weighted min/max bounds** at a non-transmission-regular graph. There are infinitely many interior examples, including every `t=1`, `m>=2`, with `alpha=1/(2m)`; every `t=2`, `m>=2`, with `alpha=1/(m*m+m+1)`; and every `t>=3`, `m=2`, with `alpha=2/(t*t+3*t+4)`.

At `alpha=0`, the numerator vanishes exactly when `(t-2)*(m-2)=2`, giving **two explicit pure-distance examples** in this split family: `K_3 join I_4` has transmissions `6` (three vertices) and `9` (four vertices), with `E*tau=8*tau` and `rho(E)=8`; `K_4 join I_3` has transmissions `6` and `8` with `E*tau=7*tau`, `rho(E)=7`. Both counter the *generic* expectation that distance transmission-scaled row equality forces transmission regularity. **They do not refute a theorem whose hypotheses or formula are different.** The broad *iff equal scaled row sums* criterion was already printed for generalized distance bounds (2022 paper in Priority check).

Independent evidence: `search.py` enumerates the full NetworkX connected, nonisomorphic atlas through `n=7`, solving the alpha equation in exact `Fraction` arithmetic (not sampling alphas). Among 980 non-transmission-regular graphs it finds **41** equality instances for `0<=alpha<1`: **39** at `0<alpha<1/2`, **two** at `alpha=0`, none at `alpha>=1/2`; the two pure-distance cases are the named split graphs. `labelled.py`, written without NetworkX, independently enumerates *all labeled* connected graphs through `n=6`, and obtains respectively `3,22,55,1448` generalized-distance equality cases at orders `n=3,4,5,6` (labeled/relabelled witnesses, not isomorphism classes); no pure-distance examples through six vertices. `checker.py` reconstructs distance matrices via Floyd--Warshall and tests the **exact eigenvector equation**, with ten `(t,m)` samples including both pure-distance witnesses. See `search.log`, `labelled.log`, `checker.log`.

### 2. Generalized distance and distance signless Laplacian: regularity at the high-alpha end

**Status: PROVED for this specific row bound, not a new literature attribution.** For a connected undirected graph with at least two vertices and `1/2<=alpha<1`, equality in either transmission-scaled bound for `M_alpha` forces transmission regularity. The distance signless Laplacian `Q_D=diag(tau)+E=2*M_(1/2)` is therefore covered. Converse: if all transmissions equal, the bounds are equalities.

**Proof.** Let `s=min_i tau_i`, `T=max_i tau_i`, and suppose all scaled row sums equal `R`. Each weighted average `U_i/tau_i` lies between `s` and `T`, since distance weights are nonnegative and sum to `tau_i`. At a vertex of transmission `T`, `R>=alpha*T+(1-alpha)*s`; at one of transmission `s`, `R<=alpha*s+(1-alpha)*T`. If `T>s` and `alpha>1/2` these inequalities contradict each other. At `alpha=1/2` they force every **other** vertex to have transmission `s` when viewed from a `T` vertex, and every other vertex to have transmission `T` when viewed from an `s` vertex, because all pairwise distances between distinct vertices are positive. If `n>=3`, a third vertex must simultaneously have transmissions `s` and `T`, impossible; at `n=2` transmissions are equal automatically. QED.

**Endpoint warning:** at `alpha=1`, `M_1=diag(tau)` is reducible; the **upper** bound `rho(M_1)=max tau` is always attained, even if transmissions vary. Hence the `alpha<1` hypothesis is essential. The atlas search found **zero** irregular `Q_D` equality cases through seven vertices, and the independent labeled search found zero through six; these are checks, not substitutes for the proof.

### 3. Signless Laplacian `Q`: equality need not imply regularity; exact classification

**Status: ALREADY KNOWN phenomenon, complete direct proof and exact checks.** Put `Q=diag(d)+A`, where `d>0` is the degree vector. Its degree-scaled row sums are `r_i=d_i+(sum_{j~i}d_j)/d_i`. In a connected graph they are constant iff the graph is regular **or semiregular bipartite** (the two bipartition parts have constant, possibly different, degrees).

**Proof.** Let `Delta` and `delta` be extreme degrees. At a `Delta` vertex, `r>=Delta+delta`; at a `delta` vertex, `r<=delta+Delta`. If all `r` equal, both are equalities, so **every** neighbor of a maximum-degree vertex has degree `delta` and **every** neighbor of a minimum-degree vertex has degree `Delta`. Connectivity propagates these two degree values to the entire graph, alternating across edges. The resulting graph is semiregular bipartite, unless `Delta=delta`, when it is regular. Conversely an `(a,b)`-semiregular bipartite graph has `Q*d=(a+b)*d`, so its spectral radius is `a+b`. Thus `K_{2,3}` is a complete nonregular witness with `rho(Q)=5`; this is not a new counterexample. The atlas test found nine nonregular isomorphism classes through seven vertices and the labeled independent enumeration found `3,4,15,21` at `n=3..6`. `checker.py` verifies three examples exactly.

### 4. Strongly connected, genuinely directed `A_alpha`: another infinite nonregular family

**Status: PROVED exact family and separately checked samples; literature priority UNRESOLVED.** Define a digraph `H_{t,m}` (`t>=2`, `m>=1`) on high vertices `0..t-1` and low vertices `t..t+m-1`: all distinct high vertices have arcs in both directions; each high vertex has an arc to every low vertex; **each low vertex has exactly one outgoing arc, to high vertex 0**. These digraphs are strongly connected and not symmetric/bidirected. Their outdegrees are `h=t+m-1` on high vertices and `1` on low vertices. Use `A_alpha=alpha*diag(d^+)+(1-alpha)*A` and the positive outdegree vector `d^+`. Its two scaled row sums are

```
R_H=alpha*h + (1-alpha)*((t-1)*h+m)/h,
R_L=alpha + (1-alpha)*h.
```

Their difference is `(h-1)*(alpha-(1-alpha)*m/h)`, so at `alpha=m/(t+2*m-1)` the two rows agree. The lemma proves `rho(A_alpha)=R_L` on a **strongly connected non-outregular digraph** for `0<alpha<1/2`. For `t=2,m=1`, explicit arcs are `(0,1),(0,2),(1,0),(1,2),(2,0)` and `alpha=1/3`, `rho=5/3`, with outdegrees `(2,2,1)`. `checker.py` independently rebuilds the arc sets, checks strong connectivity in forward and transpose directions, and tests four exact rational eigenvector certificates. This is not a counterexample to published `Q`-digraph square-root bounds, which have different formulas and parameters.

### 5. Weighted `A_alpha`: nonuniform rational-weight check

**Status: EXACT FINITE WITNESS; transferred semiregular mechanism, ALREADY KNOWN in substance.** On a weighted `K_{2,3}`, use the positive edge-weight block `[[1/2,1,3/2],[3/2,1,1/2]]` between bipartition parts. Its weighted strengths are `(3,3,2,2,2)`, not all equal. Defining `A_alpha=alpha*diag(s)+(1-alpha)*W`, `A_(1/2)*s=(5/2)*s`, and thus the degree/strength-scaled lower and upper row sums equal `5/2`. This explicitly covers **nonconstant edge weights**, but the mechanism is just weighted semiregular bipartiteness for `Q=2*A_(1/2)`; **no novelty claim**. It is checked exactly in `checker.py`; no weighted graph enumeration or arbitrary real-weight classification was attempted.

### 6. Sombor and general degree-based matrices

**Status: FINITE NEGATIVE SEARCH ONLY; a published regularity-only theorem for this *particular scaled bound* was not verified.** For the ordinary Sombor matrix `S_ij=1_(i~j)*sqrt(d_i*d_i+d_j*d_j)` and positive vector `d`, `search.py` tested `max_i (S*d)_i/d_i - min_i (S*d)_i/d_i` against `1e-9` in **980 connected nonregular atlas graphs through order seven**. **No equality** was found. This is numerical testing, with no general proof, no certified interval arithmetic for radicals, and no basis for stating that arbitrary degree-based graph matrices force regularity.

## Priority check

The following were examined, rather than merely assuming absence of prior work. Citations link to the original texts or bibliographic landing pages. **Specific formula matters:** `min/max` of the *scaled rows* must not be conflated with different two-vertex square-root estimates, even when both are described as Merris-type bounds.

| Source inspected | Statement and relevance | Priority consequence |
| --- | --- | --- |
| V. Nikiforov, *Merging the A- and Q-spectral theories*, 2017, [journal PDF](https://doiserbia.nb.rs/img/doi/1452-8630/2017/1452-86301701081N.pdf), Problem 22 and Proposition 20 | Original average-neighbor-degree `A_alpha` bounds and informal expectation; partial regularity in `alpha>1/2`. | The **undirected adjacency split family is prior Athanor work**, not new here. |
| S. Liu, K. C. Das, J. Shu, *On the eigenvalues of A_alpha-matrix of graphs*, *Discrete Math.* 343 (2020) 111917, [DOI](https://doi.org/10.1016/j.disc.2020.111917), §5, Thms. 5.4--5.6 | Published equality results for original edge bounds and certain graph classes. | Do not conflate edge-bound equality with scaled-row equality; do not claim the original issue is wholly open. |
| *New Bounds for the Generalized Distance Spectral Radius/Energy of Graphs* (2022), [full text](https://onlinelibrary.wiley.com/doi/10.1155/2022/9562730), Lemma 5 and subsequent bounds | States the transmission-weighted scaled row criterion involving `alpha*tau_i+(1-alpha)*sum_j dist(i,j)*tau_j/tau_i`, with equality iff these terms are constant. | **The generic equality criterion is already published.** The explicit `K_t join I_m` parameterization and `alpha=0` split witnesses were not found in this bounded search, **not certified new**. |
| *Bounds for the Generalized Distance Eigenvalues of a Graph*, *Symmetry* 11 (2019) 1529, [full text](https://www.mdpi.com/2073-8994/11/12/1529), Theorems 3--4 | Gives other square-root/degree-transmission bounds with equality iff transmission regular, under the paper's hypotheses. | **Not refuted** by the different `diag(tau)^(-1)M_alpha diag(tau)` row bound; compare the formulas first. |
| Alhevaz et al., *On the spectral radius of D_alpha-matrix of a connected graph*, [2024 PDF](https://acta.sapientia.ro/content/docs/2_Alhevaz.pdf), Thms. 2--3 | Transmission-regular equality in particular non-row two-vertex bounds, including high-alpha range. | No contradiction asserted; the high-alpha regularity theorem above concerns **this** scaled-row inequality. |
| Xi and Wang, *Sharp upper bounds on signless Laplacian spectral radius of strongly connected digraphs*, [DMTCS article](https://www.dmgt.uz.zgora.pl/publish/article.php?doi=1915), corollaries on digraph `Q` | Regular/semiregular-type equality statements for different signless Laplacian bounds. | Does not establish the directed `A_alpha` family is previously unknown, nor does the family disprove the cited bound. |
| *On the Distance Signless Laplacian Spectral Radius ...*, *Mathematics* 8 (2020) 792, [full text](https://www.mdpi.com/2227-7390/8/5/792) | Published transmission and distance-signless spectral bounds; the elementary `2*min(tau)<=rho(Q_D)<=2*max(tau)` is familiar. | Do not present that elementary bound or its regular equality as novel. |
| *Sombor spectra of chain graphs*, 2023, [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10395351/) | Confirms the ordinary Sombor-matrix definition; no verified match to the exact scaled-row equality claim. | **No Sombor priority conclusion** from this limited inspection. |

**Queries/search surfaces:** checked repository `verify/`, `attack/`, `status/`, `priority/`; web searches for `Merris signless Laplacian equality semiregular bipartite`, `distance spectral radius average transmission equality`, `generalized distance matrix transmission weighted row sum equality`, `A_alpha digraph degree bound equality`, `weighted A_alpha Perron row equality`, `Sombor matrix average neighbor degree equality`, and `split graph generalized distance alpha spectral radius`, with arXiv, journal, Crossref/DOI-indexed, publisher, and public indexing results. The search was **not** an exhaustive cited-by sweep or a certification that the specific families are absent from all literature. This priority limitation applies particularly to the distance and digraph constructions.

## Replay and evidence

Run from the repository root using **free local Python**; the first scan needs NetworkX, while the exhaustive labeled scan and exact checker need only the Python standard library:

```sh
OPENBLAS_NUM_THREADS=1 python3 forkmatrices/search.py > forkmatrices/search.replay.log
python3 forkmatrices/labelled.py > forkmatrices/labelled.replay.log
python3 forkmatrices/checker.py > forkmatrices/checker.replay.log
diff -u forkmatrices/search.log forkmatrices/search.replay.log
diff -u forkmatrices/labelled.log forkmatrices/labelled.replay.log
diff -u forkmatrices/checker.log forkmatrices/checker.replay.log
sha256sum -c forkmatrices/SHA256SUMS
```

Methods are **independent at the graph-generation level**: `search.py` uses NetworkX's atlas (nonisomorphic graphs through seven vertices), `labelled.py` constructs every labeled adjacency bitmask from scratch (orders two through six), `checker.py` supplies explicit rational-weight/distance/digraph certificates without NetworkX, NumPy, SymPy or eigenvalue routines. No floating-point spectral-radius solver is needed; the only floating-point test is Sombor radical row sums. **SHA256SUMS** contains SHA-256 of every other submitted file (not itself), generated after finalization. The complete general proofs above provide the separate mathematical check in addition to the executable witnesses.

## Not covered

- No general classification for *all* generalized-distance row-bound equality graphs at `alpha<1/2`; `K_t join I_m` is an infinite subfamily, not completeness. No arbitrary-order exhaustive search.
- No proof that the distance/digraph formulas have never been published. The priority audit searches only a bounded set of public sources; original papers' journal updates, theses, obscure and post-2026 indexing may have been missed.
- No claim that a different inequality with a nonlinear square-root term is false; the searched literature contains such *different* regularity-only statements.
- No classification of weighted `A_alpha` for unrestricted real weights; only one exact genuinely nonuniform example. No digraph enumeration of all strongly connected digraphs, only a proved construction and four exact checks.
- No radical-certified Sombor equality search, general Sombor proof or survey of every degree-based weighting. No nonlinear/transmission edge bounds, disconnected graphs, zero-degree undefined scaling, signed/directed distances, or `alpha=1` irreducibility assumptions.
- The atlas is exhaustive only up to seven vertices and the labeled search up to six; exact alpha solving in these finite scans does not certify continuous parameter behavior beyond the scanned graphs.

## Next steps

1. **Highest priority:** start from the proved distance split equation and search full cited-by/related literature for exactly `D_alpha=alpha*Tr+(1-alpha)*D` *transmission-scaled equality* and `K_t join I_m`. The general iff condition is published, so the potential contribution is the explicit exceptional family and especially the **two `alpha=0` pure-distance witnesses**; establish bibliographic priority before naming a new theorem.
2. Search larger graph classes for *other* pure-distance non-transmission-regular exact witnesses, e.g. joins with several regular nontrivial cells, via equitable distance quotient and a separate exact checker. Investigate a structural classification for split graphs first.
3. For the directed `A_alpha` family, audit published **outdegree** Merris bounds versus signless-digraph `Q` variants and test other strongly connected quotient digraphs; record precisely which inequalities each witness attains.
4. For Sombor/degree-based variants, find and inspect a concrete theorem with the **same scaled-row functional**, then use symbolic radicals or exact algebraic-number comparison, not `1e-9`, to test any would-be counterexample. Expand genuinely nonuniform weighted graph sampling separately.
