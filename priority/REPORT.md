# Independent priority audit: Nikiforov A_alpha spectral claims

Audit date: 2026-11-15. Prepared for Athanor lab (Aleksei Kudriashov).
Target: V. Nikiforov, "Merging the A- and Q-spectral theories," Applicable Analysis and Discrete Mathematics 11 (2017), 81-107, DOI: 10.2298/AADM1701081N; preprint arXiv:1607.03015.

## Scope and priority standard

This report checks whether each *specific* proposition (not just its topic) appears in earlier published literature, and distinguishes an explicit result from a result derivable by combining published facts. "Not found" is a bounded-search outcome, not proof of novelty. The original 2016 arXiv version and the 2017 journal publication differ: they must not be cited interchangeably. Numerical examples and derivations below are independent mathematical checks, not evidence of publication.

VERDICT SUMMARY:
1. NOT FOUND (exact nonregular equality family or named examples); closely related bounds/equality theorems are published.
2. ALREADY PUBLISHED (the exact spectral polynomial and hence the eigenvalue formula); the global-optimality conjecture, stated for all alpha<1, is false for alpha>1/2.
3. PARTLY PUBLISHED (the exact value and unique optimizer follow immediately from explicitly published ingredients; no separately stated A_alpha global-extremum theorem located).
4. PARTLY PUBLISHED (a 2018 paper gives the stronger corrected equality characterization; the proposed correction to Proposition 29 is ALREADY printed in the 2017 journal version; the small nondifferentiability example was not found as an explicit publication).

## 1. Nonregular equality in Nikiforov's (20) and (21)

**Verdict: NOT FOUND as an explicitly published family; nearby work is not equivalent.**

For G = K_t join the empty graph on m vertices, t,m >= 2, set H=t+m-1 (degree of each clique vertex), L=t (degree of each independent-set vertex). Their average neighbor degrees are respectively
M_H = ((t-1)H + m*t)/H and M_L = H.
Nikiforov's Proposition 20 (2017, p.92, equations (20)-(21)) bounds the Perron root between min and max over vertices of alpha*d(v)+(1-alpha)*M(v). These expressions are constant on V(G) exactly when alpha=m/(2m+t-1); indeed
alpha*(H-L) = (1-alpha)*(H-M_H) = (1-alpha)*m*(m-1)/H.
With that alpha the degree vector is a positive A_alpha eigenvector:
A_alpha*d = [alpha*t+(1-alpha)*H]*d.
Thus both (20) and (21) have equality with lambda_1=alpha*t+(1-alpha)*(t+m-1), in spite of G being connected and nonregular and 0<alpha<1/2.
Checks: K_2 join 2K_1 = K_4 minus one edge, alpha=2/5, lambda_1=13/5. For G=K_1 join K_{2,2}, the degree classes are 4 and 3, and alpha=1/4 gives both average-degree expressions and lambda_1 equal to 13/4. The latter is an additional nonregular example, not a member of the t,m>=2 family.

**Prior publication examined:** Nikiforov (2017), Proposition 20 and Problem 22, pp.92-93: the paper proves regularity for connected graphs when 1/2<alpha<1 and specifically speculates about remaining equality cases, but does not give this family. Shuting Liu, Kinkar Chandra Das and Jinlong Shu, "On the eigenvalues of A_alpha-matrix of graphs," Discrete Mathematics 343(8) (2020), article 111917, DOI: 10.1016/j.disc.2020.111917, Section 5, Theorems 5.4-5.6: Theorem 5.4 treats equality in Nikiforov's *edge* inequalities (22)-(23), not the *average-neighbor-degree* inequalities (20)-(21); Theorem 5.6 treats disconnected-graph equality for (20)-(21) only when 1/2<alpha<1. The paper's Lemmas 2.9-2.10 characterize other degree relations, not the present weighted equality below 1/2. These are partial topical precedents, not a publication of this claim.

Sources:
- Nikiforov journal PDF, pp.92-93: https://doiserbia.nb.rs/img/doi/1452-8630/2017/1452-86301701081N.pdf
- Liu, Das and Shu (2020), full theorem text at publisher: https://doi.org/10.1016/j.disc.2020.111917

## 2. Second eigenvalue of K_1 join (K_s union K_s), odd n=2s+1

**Verdict: ALREADY PUBLISHED for the exact spectral formula (not a new eigenvalue observation). The proposed unrestricted global optimality is false.**

Jie Xue, Ruifang Liu, Guanglong Yu and Jinlong Shu, "The Multiplicity of A_alpha-Eigenvalues of Graphs," Electronic Journal of Linear Algebra 36 (2020), 645-657, **p.653, equation (3.7)**, print the full A_alpha characteristic polynomial for precisely K_1 join 2K_s:

P(x) = (x-(s+1)*alpha+1)^(2s-2) * (x-s-alpha+1)
       * [x^2+(1-s-2*alpha*s-alpha)*x+2*alpha*s-2*s+2*alpha*s^2].

Thus x=s-1+alpha is explicitly an eigenvalue for all 0<=alpha<1; comparison with the other roots (or the equitable three-cell quotient) identifies it as lambda_2. Their Lemma 3.5 on the same page directly revisits this polynomial. This is an exact prior, not a conjectural reconstruction. The paper studies s>=2; the s=1 instance follows by direct calculation.

The 2017 Nikiforov journal article (p.103, after Proposition 31) notes the odd-order gap using K_s union K_{s+1}; p.104-105 Proposition 37 also supplies general join-of-regular-graphs spectral machinery.

**Important correction to proposed optimality:** for n=2s+1 and alpha>1/2, Proposition 36 of Nikiforov (2017, p.104) gives lambda_2(A_alpha(K_n))=n*alpha-1, and
(n*alpha-1)-(s-1+alpha)=s*(2*alpha-1)>0.
Therefore K_1 join 2K_s cannot maximize lambda_2 over all n-vertex graphs on the full range 0<=alpha<1. A restricted conjecture for alpha<=1/2 would be a different question; this audit does not assert it is optimal even there.

Sources:
- Xue, Liu, Yu and Shu (2020), ELA p.653, (3.7) and Lemma 3.5: https://journals.uwyo.edu/index.php/ela/article/download/4999/5171/6923
- Nikiforov (2017), pp.103-105: https://doiserbia.nb.rs/img/doi/1452-8630/2017/1452-86301701081N.pdf

## 3. Even-n minimum of lambda_min(A_alpha), 0<alpha<1/2

**Verdict: PARTLY PUBLISHED; value and uniqueness are short corollaries of pre-existing printed theorems, but no explicit stand-alone global A_alpha extremal theorem was located.**

Set Q=D+A, so A_alpha = (1-2*alpha)*A + alpha*Q. Since Q is positive semidefinite and 1-2*alpha>0, for every n-vertex G,
lambda_min(A_alpha(G)) >= (1-2*alpha)*lambda_min(A(G))
                          >= -(1-2*alpha)*sqrt(floor(n/2)*ceil(n/2)).
For even n this equals -(n/2)*(1-2*alpha).
The final adjacency bound and its unique equality graph are printed in Liu, Das and Shu (2020), **Lemma 2.5**, citing an earlier source [6]: equality for the adjacency extremum is exactly K_{ceil(n/2),floor(n/2)}. Nikiforov (2017) already states the adjacency extremum in **p.102**, immediately before Problem 30, and prints the A_alpha eigenvalues of every K_{a,b} in **Proposition 38, p.105**. For K_{n/2,n/2}, the sign-changing bipartition vector is simultaneously a zero eigenvector of Q and a -n/2 eigenvector of A, giving equality. If A_alpha achieves equality, it must achieve equality in the adjacency inequality, hence the balanced complete bipartite graph is the unique extremizer up to isomorphism.

The argument solves Problem 30 for even n, but its core ingredients were already published in 2017 and before. We did not find a paper stating this combined extremal conclusion verbatim; priority must be qualified accordingly.

Sources:
- Nikiforov (2017), Problem 30 p.102 and Proposition 38 p.105: https://doiserbia.nb.rs/img/doi/1452-8630/2017/1452-86301701081N.pdf
- Liu, Das and Shu (2020), Lemma 2.5: https://doi.org/10.1016/j.disc.2020.111917

## 4. Differentiability question; Proposition 29; Proposition 31(a)

**Verdict: PARTLY PUBLISHED, with a major preprint-versus-journal correction.**

(a) The differentiability question is **Question 5 in arXiv:1607.03015**, but **Question 6 in the 2017 journal version, p.86**. For G=K_2 union K_1, the full eigenvalue multiset is {1, 2*alpha-1, 0}. Consequently
lambda_2(A_alpha(G))=max(0,2*alpha-1),
which has left derivative 0 and right derivative 2 at alpha=1/2. The example refutes differentiability of every ordered eigenvalue, but we did not locate its explicit prior publication in the queried sources. It is an elementary eigenvalue-crossing observation, so search non-detection should not be treated as a robust novelty claim.

(b) **Proposition 29 is wrong as printed in the arXiv preprint**, where the bound uses 2*alpha*m/n - 2*(1-alpha)*maxcut(G)/n. For K_3, alpha=0, this incorrectly forces -1<=-4/3. **However, the peer-reviewed 2017 journal article ALREADY PRINTS THE CORRECT VERSION** at **p.102, Proposition 29**:
lambda_min(A_alpha(G)) <= 2*m/n - 4*(1-alpha)*maxcut(G)/n.
The preceding journal derivation uses the {+1,-1} cut vector and is valid. It would be incorrect to present this as a correction to the 2017 published paper: it corrects the preprint only.

(c) Nikiforov (2017), p.102-103, Proposition 31(a), says for alpha>1/2 that equality lambda_2=alpha*n-1 forces G=K_n. This uniqueness statement is false for K_4 minus one edge, which has two degree-3 (universal) vertices and lambda_2=4*alpha-1 for alpha>1/2. **This entire correction was already published, in stronger form**: Huiqiu Lin, Jie Xue and Jinlong Shu, "On the A_alpha-spectra of graphs," Linear Algebra and its Applications 556 (2018), 210-219, DOI 10.1016/j.laa.2018.07.003, **Theorem 1.1, p.211**: for alpha>1/2 and k>=2, lambda_k(A_alpha(G))=alpha*n-1 iff G has k vertices of degree n-1. Take k=2 and G=K_4-e. The theorem was inspected in full text, not only the abstract.

Sources:
- Nikiforov peer-reviewed journal (2017), pp.86, 102-103: https://doiserbia.nb.rs/img/doi/1452-8630/2017/1452-86301701081N.pdf
- Nikiforov original arXiv version: https://arxiv.org/pdf/1607.03015
- Lin, Xue and Shu (2018), Theorem 1.1 (full-text preprint p.2): https://arxiv.org/pdf/1709.00182 ; journal DOI: https://doi.org/10.1016/j.laa.2018.07.003

## Search record and boundaries

Date of search: 2026-11-15. Sources inspected in full text at relevant sections: Nikiforov's 2017 published PDF, especially pp.86, 92-93, 102-105; arXiv:1607.03015 (version-specific discrepancies); Lin, Xue and Shu (2018), Theorem 1.1; Xue, Liu, Yu and Shu (2020), equation (3.7) and Lemma 3.5; Liu, Das and Shu (2020), Lemmas 2.5, 2.9-2.10 and Theorems 5.4-5.6.

Web searches included combinations of: "A_alpha average neighbour degree equality"; "Nikiforov Problem 22 equality"; "Merris A_alpha equality nonregular join"; "K1 join 2Ks second eigenvalue"; "A_alpha least eigenvalue balanced complete bipartite"; "Problem 30 Nikiforov"; "A_alpha eigenvalue differentiability"; "Proposition 29 maxcut"; and "Proposition 31(a) equality". Searched for matching discussion/papers across arXiv, MathOverflow, Google Scholar-indexed pages, Semantic Scholar-indexed pages, OpenAlex-indexed metadata and named outlets Discrete Mathematics, Linear Algebra and its Applications, Linear and Multilinear Algebra, Electronic Journal of Linear Algebra, and Applicable Analysis and Discrete Mathematics.

**Citation coverage limitation:** The public web search interface did not provide a complete harvest of the approximately 300 citing works from Google Scholar, Semantic Scholar or OpenAlex. Their complete citation lists were **not** traversed paper by paper, and an exhaustive negative-priority certificate is therefore not possible. No fresh priority search should claim that all citing papers were screened. Some publisher content was accessible as rendered full text rather than a downloadable PDF; theorem and lemma wording was checked when accessible. A full bibliometric harvest and manual inspection of all potentially relevant citing full texts remains an outstanding task.

## Publication recommendation

Do not advertise claim 2's eigenvalue factor as new: cite ELA (2020), p.653, (3.7), and amend the global-optimality conjecture (it is false for alpha>1/2). For claim 4, cite the corrected 2017 journal Proposition 29 and Lin-Xue-Shu (2018), Theorem 1.1, rather than presenting either as an unpublished correction. Claim 3 can be presented as a short deduction/explicit even-n resolution with the old adjacency extremum and Nikiforov's bipartite spectrum fully credited. Claim 1 is the most plausible original observation among the four, but subject to the uncompleted exhaustive citation sweep.
