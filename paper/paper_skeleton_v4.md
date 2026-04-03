# Projected χ² for multi-distribution neutrino cross-section measurements
## Paper Skeleton (v4 — final)

Target: PRD or JINST, ~12 pages
Standalone methods paper, small author list, NOT a MicroBooNE collaboration paper.

---

## Abstract (~150 words)

When multiple differential cross-section distributions are extracted from a
common event sample with correctly evaluated inter-distribution correlations,
the statistical covariance matrix is generically singular. Exact linear
constraints from event sharing — such as equal distribution totals when all
selected events populate bins in every distribution — produce directions with
zero statistical variance. We show that the number of such degenerate
directions is exactly predicted by N_null = N_bins - rank(A), where A is the
bin-combination matrix recording the distinct patterns of bins occupied by
events across distributions. Systematic uncertainties may partially or fully
restore invertibility depending on their nature, but near-singular directions
produce unphysically inflated chi-squared values in global goodness-of-fit
tests. We present a projection method that restricts the test to the subspace
carrying independent statistical information, yielding a chi-squared with
N_bins - N_null degrees of freedom. The method is connected to established
results on hypothesis testing with singular covariance matrices and validated
with both an analytical toy model and a realistic simulated neutrino-argon
cross-section measurement including D'Agostini unfolding and multi-source
systematic uncertainties.

---

## I. Introduction (~1 page)

### Multi-distribution measurements (paragraph 1)

Modern neutrino cross-section measurements increasingly report multiple
differential distributions from the same event sample. Each distribution
represents a different kinematic projection of the selected events. Cite:

Same-event multi-distribution measurements:
  - MicroBooNE QE-like multi-diff: 2301.03700, 2301.03706
  - MicroBooNE 0p/Np simultaneous: 2402.19281, 2402.19216
  - MicroBooNE mesonless CC: 2403.19574
  - MicroBooNE GKI: 2310.06082
  - T2K final-state kinematics: 1802.05078
  - T2K TKI in CC1pi+: 2102.03346
  - NOvA double-differential: 2109.12220, 2410.05526, 2410.10222

Context (shared systematics, different events — mention briefly):
  - MINERvA multi-target: 2301.02272, 2209.07852
  - T2K O/C simultaneous: 2004.05434
  - T2K nu/nubar combined: 2002.09323

Draw the distinction: multi-target measurements share systematic uncertainties
but not events. Multi-distribution measurements from the same sample share
both. This paper addresses the latter.

### Blockwise unfolding (paragraph 2)

The blockwise unfolding framework [Gardiner 2024] provides a rigorous method
for computing the full cross-distribution covariance matrices in all such
measurements. In this approach, all bins for the same observable are organized
into blocks, each corresponding to a single differential distribution. Each
block is unfolded independently using any standard method (D'Agostini,
Wiener-SVD), while the full covariance matrix encoding inter-block correlations
is constructed and propagated through the extraction. The statistical
covariance between bins in different blocks is determined by the number of
events they share. The result is a single combined covariance matrix covering
all reported bins across all distributions.

### Global goodness-of-fit and the problem (paragraph 3)

A natural use of these combined measurements is computing global
goodness-of-fit metrics, such as a chi-squared statistic, testing a model
against all distributions simultaneously. Compared to evaluating each
distribution independently, a global test that accounts for all
inter-distribution correlations extracts the maximum discriminating power
from the data. Such global metrics are essential for interaction model tuning
[Wilkinson 2016, MINERvA GENIE tune 2019, GENIE CC0pi 2022], where multiple
kinematic distributions provide complementary sensitivity to different aspects
of the underlying physics. However, as we show in this paper, the full
inter-block statistical covariance is generically singular when it correctly
accounts for correlations from shared events. The rank deficiency arises from
exact linear constraints imposed by the event-sharing structure of the binning
scheme: certain linear combinations of bins, such as the difference in total
event counts between two distributions of the same events, have identically
zero variance. While systematic uncertainties may formally restore
invertibility, the resulting near-singular directions produce unphysically
inflated chi-squared values, yielding a global test with the wrong effective
number of degrees of freedom.

### This paper (paragraph 4)

We present a projection method that removes these degenerate directions from
the goodness-of-fit test by restricting it to the subspace carrying
independent statistical information. We derive an exact formula for the number
of degenerate directions in terms of the bin-combination matrix recording
which bins each event populates across distributions. We classify the null
directions into structural constraints (determined by the block topology) and
kinematic constraints (determined by the phase-space occupancy), and show how
different classes of systematic uncertainties interact with each type. The
method is validated using an analytical toy example and a realistic scenario
with simulated neutrino-argon interactions including D'Agostini unfolding
and multi-source systematics. We connect the approach to the established
statistical literature on hypothesis testing with singular covariance matrices.

### Outline (paragraph 5)

Brief roadmap: Sec. II derives the singularity and the combination matrix
formula. Sec. III presents the projection method. Sec. IV validates it
analytically. Sec. V demonstrates it realistically. Sec. VI provides practical
guidance. Sec. VII connects to the statistics literature. Sec. VIII summarizes.

---

## II. Origin of the singularity (~2.5 pages)

### II.A. Setup

Define the measurement scenario precisely:
- B blocks, block b has N_b bins (including overflow). N_bins = sum_b N_b.
- Common selected event sample. Each selected event populates exactly one
  bin in each block it belongs to. (Events may not belong to all blocks if
  blocks correspond to different multiplicity slices.)
- Blockwise unfolding per Gardiner: unfolding matrix U and error propagation
  matrix E are block-diagonal. Full inter-block covariance is constructed in
  reco space and propagated through extraction.
- Final result: vector of N_bins cross-section values and N_bins × N_bins
  total covariance matrix C_total = C_stat + C_syst.

### II.B. Statistical covariance from shared events

Review the shared-event covariance formula (Gardiner Eq. 53):

  Cov(D_X, D_Y) = D_v

where D_v is the number of events populating both bins X and Y. This is exact,
following from Poisson statistics of non-overlapping sub-populations. After
propagation through block-diagonal unfolding (Gardiner Eqs. 22, 29, 54-55),
the statistical covariance in unfolded space inherits the same rank structure.

### II.C. The combination matrix and the null space

**The indicator vector.** For each event e, define phi_e in {0,1}^{N_bins}
with a 1 in each bin the event populates across all blocks and 0 elsewhere.
The statistical covariance can be written:

  C_stat = sum_e w_e^2 phi_e phi_e^T                              (key eq)

where w_e is the weight of event e.

**Rank of C_stat.** The rank of this matrix equals the dimension of
span({phi_e}), which depends only on the distinct indicator vectors present,
not on their weights or multiplicities. Define the bin-combination matrix A
as the matrix whose rows are the unique indicator vectors among {phi_e}.
Then:

  rank(C_stat) = rank(A)                                          (key eq)

and the number of degenerate directions is:

  N_null = N_bins - rank(A)                                       (key eq)

This is computable directly from the MC event list before constructing any
covariance matrix.

**Proof:** rank(C_stat) = dim(span({phi_e})) since w_e > 0 for all events
that contribute (zero-weight events are excluded). The span of {phi_e} equals
the row space of A by construction. Therefore rank(C_stat) = rank(A). QED.

**Null space condition.** A vector w is in the null space of C_stat if and
only if, for every event e:

  sum_{a in bins(e)} w_a = 0

where bins(e) is the set of bins event e populates. This imposes one linear
equation per unique indicator vector, yielding the system A w = 0. The null
space is therefore ker(A).

### II.D. Types of constraints

**Structural constraints.** Constraints that hold for ANY kinematic
distribution, depending only on the block/slice topology. They follow from the
fact that every event of a given type appears in one bin of each relevant block.

Example 1 — Equal block totals: If every event appears in blocks A and B, then
sum(block A) = sum(block B) identically. The vector w = (1,...,1,-1,...,-1) is
in ker(A) regardless of the occupancy pattern.

Example 2 — Equal slice totals: If events are categorized by multiplicity (0p,
1p, Np) and each category has its own sub-blocks for different observables,
each category's total is the same across its sub-blocks.

Counting: B blocks sharing all events, no slicing → N_struct = B - 1.
S multiplicity slices, slice s appearing in B_s blocks → N_struct = sum_s (B_s - 1).

**Kinematic constraints.** Additional null directions arising from empty cells
in the multi-block occupancy pattern. When no event simultaneously populates a
particular combination of bins across blocks, that combination imposes no
constraint on w, enlarging the null space.

For the special case of 2 blocks with M and K bins, the combination matrix
corresponds to the bipartite graph between p_T bins and cosθ bins, with edges
for occupied cells. N_null equals the number of connected components of this
graph. For a fully populated table: 1 component, N_null = 1 (purely
structural). For sparse tables: more components, N_null > 1.

For the general multi-block case, the bipartite graph generalizes and rank(A)
is the appropriate computational tool.

The total number of null directions is:

  N_null = N_struct + N_kin

where N_kin >= 0 depends on the physics and binning granularity. N_struct
provides an analytically computable lower bound.

### II.E. Effect of systematic uncertainties on the constraint structure

Different systematic treatments interact differently with the null directions.

**1. Reweighting systematics (flux, xsec model, reinteraction).**
Per-event reweighting assigns a single weight w_i to each event, contributing
identically to all bins that event populates. Events do not move between bins.
The set of occupied bin-combinations — and therefore the matrix A — is
identical in every universe. All null directions are preserved. Adding
reweighting-based C_syst to C_stat does not change the rank.

This is exact: it holds for any reweighting scheme, including flux
normalization, cross-section model shape variations, and hadronic
reinteraction systematics.

**2. Detector response systematics.**
Variations that modify reconstructed kinematic quantities (momentum scale,
angular resolution) cause events to migrate between reco bins. This can
populate previously empty bin-combinations, lifting kinematic constraints.

However, detector variations do NOT break structural constraints: each event
still passes or fails selection (affecting all blocks equally), and if it
passes, it still populates exactly one bin per block. Block totals and slice
totals remain equal.

Effect: N_null may decrease by up to N_kin (some or all kinematic constraints
lifted), but structural constraints survive.

**3. Factorized cross-section model variations (Gardiner Eq. 18).**
The MicroBooNE-style prescription computes predicted event counts as the
varied response matrix applied to fixed signal counts:

  n_a^u = B_a^u + sum_mu Delta_{a,mu}^u phi_mu^CV

Because the response matrix is block-specific but phi^CV is shared, the
per-block totals are no longer forced to agree across universes. This breaks
structural constraints.

Effect: remaining structural constraints lifted. C_total may become full rank.

**The practical consequence is the same in all cases.** Whether the eigenvalues
of C_total in the null directions of C_stat are exactly zero (reweighting
only), small (detector variations), or merely much smaller than the physical
eigenvalues (Eq. 18 prescription), including these directions in a chi-squared
test via standard inversion produces inflated values with the wrong effective
degrees of freedom.

### II.F. Ignoring inter-block correlations

Briefly discuss the alternative of ignoring correlations entirely:
block-diagonal C_stat, full rank, invertible, but wrong — double-counts
statistical information. Summing per-block chi-squared values gives a test
with too many effective degrees of freedom. Historically common in the
literature; a different kind of error from the inflation described above.

---

## III. The projection method (~2 pages)

### III.A. Constructing the projector

Eigendecompose the MC statistical covariance:

  C_stat = V Lambda V^T

with eigenvalues lambda_1 >= ... >= lambda_{N_bins}. Identify the N_null
eigenvalues that are zero (or numerically negligible). The threshold for
distinguishing zero from nonzero eigenvalues can be set relative to the
largest eigenvalue (e.g., lambda_i / lambda_1 < 10^{-10}); in practice there
is a clear gap of many orders of magnitude.

**Validation step.** Verify that N_null matches the prediction from rank(A).
Classify the null eigenvectors:
  - Structural: uniform within blocks (sum rules). Check against expected
    block-total and slice-total constraints.
  - Kinematic: non-uniform within blocks (empty cells). Check against the
    occupancy pattern.

Construct the orthogonal projector onto the range of C_stat:

  P = sum_{i=1}^{N_bins - N_null} v_i v_i^T                     (key eq)

Equivalently P = I - P_null where P_null projects onto the null space.

### III.B. Projected total covariance and chi-squared

Apply the projector to the total covariance:

  C_tilde = P C_total P                                          (key eq)

This is N_bins × N_bins but has rank N_bins - N_null. It has N_null zero
eigenvalues by construction. Compute its pseudo-inverse via eigendecomposition:

  C_tilde^+ = sum_{i=1}^{N_bins - N_null} (1/tilde_lambda_i) tilde_v_i tilde_v_i^T

The projected chi-squared:

  chi^2_proj = r^T C_tilde^+ r                                   (key eq)

where r is the residual (data - prediction), with predictions smeared by A_C
if a regularization matrix is used. Under the null hypothesis with Gaussian
residuals, this is chi-squared distributed with N_bins - N_null d.o.f.

### III.C. Formal justification

**Theorem** (Rao & Mitra 1971, Khatri 1968, Moore 1977): If x ~ N(0, Sigma)
where Sigma is positive semi-definite with rank q, then x^T Sigma^+ x ~
chi^2(q).

In our case, the residuals r live in the full N_bins-dimensional space with
covariance C_total (which may or may not be full rank). We project onto the
q = N_bins - N_null dimensional subspace (the range of C_stat) where the
independent statistical information lives. The projected residual Pr has
covariance P C_total P.

**Condition for validity:** P C_total P must be positive definite on range(P).
This requires that C_syst does not introduce additional null directions within
the range of C_stat. This is generically satisfied — it would fail only if
the systematic model assigned zero uncertainty to some physically meaningful
direction, which would indicate a pathological systematic treatment.

**Connection to Davidov et al. (2018):** Our method is a specific instance of
the "trimmed inverse" approach in their framework for testing with
singular/nearly-singular covariance matrices. The key difference is that our
trimming threshold is determined by the known structural and kinematic
constraints of the measurement (via the null space of C_stat) rather than by
a numerical criterion on the eigenvalues of C_total.

**Important distinction:** This is NOT the same as simply computing chi^2 with
the full C_total^{-1} (when it exists) and quoting N_bins - N_null d.o.f.
The projection changes WHICH directions contribute to chi^2, not just the
d.o.f. count.

### III.D. Comparison to alternatives

Discuss alternatives briefly, in prose (not a list):

(1) Ignoring inter-block correlations: wrong covariance, double-counts
statistical information.

(2) Naive inversion of C_total: impossible when C_total is exactly singular
(reweighting-only systematics). When C_total is nearly singular, gives
inflated chi^2 dominated by structurally degenerate directions.

(3) SVD truncation on C_total directly: identifies near-singular directions
numerically but lacks a physics criterion for the cutoff. Different SVD
thresholds give different d.o.f. Our method uses the known structure from
C_stat, giving an unambiguous, verifiable answer.

(4) Regularized/winsorized inversion (Davidov et al.): appropriate for the
"nearly singular" case where the true rank is uncertain. In our setting, the
structural constraints make many of the degeneracies exact and known, so
trimming (projecting) is the natural choice.

(5) Summing per-block chi-squared values: each block's chi^2 is valid on its
own, but the sum double-counts shared statistical information, giving the
wrong total d.o.f.

---

## IV. Analytical demonstration (~2 pages)

### IV.A. Setup

**Model:** 2D kinematic space (p_T, cos_theta).
  f(p_T) = p_T^k exp(-k p_T / p_0)
  g(cos_theta) = exp(kappa cos_theta)
  p(p_T, cos_theta) = f * g * (1 + alpha * p_T * cos_theta)
  Parameters: alpha = 0.3, k = 3, p_0 = 0.2, kappa = 2.
  N_expected = 1000 events.

**Binning:**
  Block 1 (dsigma/dp_T): 5 bins.
  Block 2 (dsigma/d(cos_theta)): 5 bins.
  N_bins = 10.

**Simplifications:** No background, perfect efficiency, no bin migration
(identity unfolding). This isolates the shared-event effect.

**Systematics:** Reweighting only (10% normalization + shape variations).

**Key feature:** Because all systematics are reweighting, C_total is EXACTLY
singular (rank 9). The projection method is not just helpful — it is the
ONLY way to compute any chi-squared. This is the strongest demonstration.

### IV.B. Covariance structure

Show the explicit form of the 10 × 10 C_stat:

  C_stat = | diag(n_i)     N_{ij}     |
           | N_{ij}^T      diag(m_j)  |

where N_{ij} is the 5 × 5 contingency table of the 2D histogram. Diagonal
entries: n_i = sum_j N_{ij} (row sums), m_j = sum_i N_{ij} (column sums).

**Combination matrix:** With all 25 cells populated, A is 25 × 10 with each
row having exactly two 1s (one in block 1, one in block 2). rank(A) = 9.

**Null vector:** w = (1,1,1,1,1,-1,-1,-1,-1,-1)/sqrt(10). Verify analytically
that C_stat w = 0.

**Systematic covariance:** C_syst has zero variance in the null direction
because reweighting preserves block totals. Show: w^T C_syst w = 0.
Therefore C_total = C_stat + C_syst is also rank 9.

### IV.C. Eigenvalue analysis

**Figure 1:** Eigenvalue spectrum of C_stat. 9 nonzero eigenvalues, 1 at zero
(or ~10^{-12} numerically). Clear gap of many orders of magnitude. Annotate
the null eigenvector.

### IV.D. Pseudo-experiment validation

Generate 10,000 pseudo-experiments from the true distribution. For each,
compute chi-squared three ways:

(a) **Block-diagonal** (ignoring inter-block correlations): Full inversion of
the block-diagonal part of C_total. Compared to chi^2(10).

(b) **Full C_total, naive inversion:** IMPOSSIBLE — matrix is exactly singular.
State this explicitly as a demonstration.

(c) **Projected chi-squared:** P from C_stat null space, C_tilde = P C_total P,
pseudo-inverse. Compared to chi^2(9).

(d) **Per-block chi-squared sum:** block 1 chi^2(5) + block 2 chi^2(5),
treating blocks independently. Compared to chi^2(10).

**Figure 2 (MONEY PLOT):** Chi-squared distributions for methods (a), (c), (d)
overlaid with chi^2(9) and chi^2(10) reference curves. Method (b) absent
because it's impossible. Only method (c) matches its reference curve.

**Figure 3:** p-value distributions for each computable method. Under H_0,
should be uniform. Only method (c) gives uniform p-values.

**Table 1:** Mean and variance for each method vs theoretical expectations
(mean = n_dof, variance = 2 * n_dof).

### IV.E. Combination matrix verification

Show explicitly:
- Compute rank(A) from the occupancy pattern → predicts N_null = 1.
- Compare to eigenvalue count → matches.
- For the 2-block case, show equivalence to bipartite graph connectivity:
  fully populated → 1 component → N_null = 1.

### IV.F. Sparse extension (kinematic constraints)

Brief extension demonstrating kinematic constraints. Use a modified 2D
distribution that leaves some (p_T, cosθ) cells empty by imposing kinematic
correlations (e.g., high p_T only at forward angles).

- Show occupancy matrix with empty cells marked.
- Compute rank(A): now less than 9, so N_null > 1.
- Verify against eigenvalue spectrum.
- Show bipartite graph with multiple connected components.
- Progressive cell-filling experiment: add events to empty cells one by one,
  show N_null decreases by 1 each time a previously disconnected component
  gets bridged, unchanged otherwise.

**Figure 4:** Occupancy matrices for the fully populated and sparse cases,
with connected components highlighted. Annotate N_null for each.

---

## V. Realistic demonstration (~3 pages)

### V.A. Simulation setup

**Event generation:**
- CV MC: GENIE v3 (specify tune) on argon, BNB-like flux.
  Signal: CC0pi (nu_mu CC, no pions in final state).
- Fake data: NuWro (specify version), same flux and signal definition.
  Provides genuine model discrepancy.
- Sufficient MC statistics for migration matrices and covariance estimation.

**Detector effects:**
- Gaussian smearing on reconstructed p_mu and cos_theta_mu
  (sigma_p/p ~ 5%, sigma_theta ~ 2 degrees).
- Efficiency: sigmoid turn-on in p_mu, mild angular dependence.
  Parametrized with tunable parameters for systematic variations.
- Background: ~10% CC1pi misidentified as CC0pi, from GENIE.

**Unfolding:**
- D'Agostini with 4 iterations, independently within each block.
- Error propagation matrix E computed per block (Eq. 24 of Gardiner for
  the case with A_C, or Eq. 32 without).
- Block-diagonal overall U and E per Gardiner Eqs. 54-55.

### V.B. Example 1: Two blocks (N_null >= 1)

**Binning:**
- Block 1: dsigma/dp_mu (8 bins + overflow)
- Block 2: dsigma/d(cos_theta_mu) (8 bins + overflow)
- N_bins = 18

**Combination matrix prediction:** With fully populated 2D histogram and
overflow bins: expect N_null = 1 (purely structural). If some cells are
empty due to kinematic correlations: N_null > 1.

**Systematic model (3 sources):**
1. Flux normalization (reweighting): 10% overall scale. ~100 universes.
2. Cross-section model (reweighting): Vary distribution parameters
   (or GENIE knobs if available). ~100 universes.
3. Detector response (smearing variations): Vary sigma_p, sigma_theta
   by ~20%. ~100 universes. This is the only source that changes
   reconstructed quantities and can lift kinematic constraints.

C_syst computed via standard universe formula (Gardiner Eq. 17).
C_stat from shared-event formula propagated through unfolding.
C_total = C_stat + C_syst.

**Show the constraint hierarchy (key demonstration):**

Step 1: C_stat alone. Count null eigenvalues. Verify against rank(A).

Step 2: Add flux + xsec reweighting systematics. Show rank unchanged.
Compute w^T C_syst w in each null direction — should be zero.

Step 3: Add detector smearing variations. If any kinematic constraints
are lifted, show which ones and verify against the change in occupancy
pattern (previously empty cells now populated under smearing variations).
Structural constraint (equal block totals) survives.

**Table 2:** Rank of C_stat + cumulative systematic sources. Shows the
hierarchy explicitly.

**Figure 5:** Eigenvalue spectrum of C_stat in unfolded space. Identify
null eigenvalues and verify against combination matrix prediction.

**Figure 6:** Eigenvalue spectra of C_total under cumulative addition of
systematic sources. Show the progression: some eigenvalues lift from zero
to small values, while the physical eigenvalues are much larger.

**Pseudo-experiment validation:**
Generate 10,000 pseudo-experiments by Poisson-fluctuating the GENIE
prediction, unfold each, compute chi-squared with projected method.

**Figure 7:** Chi-squared distributions from pseudo-experiments. Show
projected method matches chi^2(N_bins - N_null). Compare to naive
inversion if C_total is invertible, showing inflation.

**NuWro vs GENIE comparison:**
- Show both distributions with GENIE prediction and NuWro "data" overlaid.
- Report chi-squared values:
  * Naive (if invertible): chi^2 / N_bins, p-value
  * Projected: chi^2 / (N_bins - N_null), p-value
  * Per-block sum: chi^2 / N_bins, p-value
- Discuss: projected chi^2 gives a physically meaningful statement about
  model tension. Naive and per-block values are misleading.

**Figure 8:** NuWro vs GENIE: the two distributions with predictions overlaid.

**Table 3:** Chi-squared values and p-values for all methods, NuWro vs GENIE.

### V.C. Example 2: Four blocks with multiplicity slicing (N_null > 1)

**Binning:**
- Block 1: dsigma/dp_mu for 0p events (6 bins)
- Block 2: dsigma/dp_mu for 1+p events (6 bins)
- Block 3: dsigma/d(cos_theta_mu) for 0p events (6 bins)
- Block 4: dsigma/d(cos_theta_mu) for 1+p events (6 bins)
- N_bins = 24

**Structural constraints:**
- Sum(block 1) = Sum(block 3) [same 0p events] → 1 constraint
- Sum(block 2) = Sum(block 4) [same 1+p events] → 1 constraint
- N_struct = 2

**Kinematic constraints likely** because 0p and 1+p sub-populations have
different kinematic distributions, leading to sparser occupancy in the
per-slice 2D histograms.

**Same systematic model and validation as V.B.** Show:
- Total N_null from eigenvalue spectrum.
- Verify against rank(A) from combination matrix.
- Decomposition into structural + kinematic components.
- Pseudo-experiment chi^2 matches chi^2(N_bins - N_null).
- NuWro vs GENIE with projected chi^2.

**Figure 9:** Eigenvalue spectrum for the 4-block case.

**Figure 10:** Chi-squared distributions for the 4-block case.

**Key point:** N_null is predictable (structural part as lower bound) and
exactly verifiable (rank(A)). The method adapts automatically regardless
of the number of blocks and constraints.

### V.D. Practical observations

Summarize findings across both examples:

- **Eigenvalue gap:** Clear separation (many orders of magnitude) between
  null/near-null and physical eigenvalues in all cases.
- **Combination matrix accuracy:** rank(A) matches eigenvalue-based N_null
  exactly in every case tested.
- **Null eigenvector classification:** Structural eigenvectors are uniform
  within blocks and correspond to expected sum rules. Kinematic eigenvectors
  have non-uniform structure and correspond to empty occupancy cells.
- **Systematic hierarchy confirmed:** Reweighting preserves all constraints.
  Detector variations can lift kinematic ones. Whether structural ones are
  broken depends on the specific systematic treatment (Eq. 18 vs direct
  reweighting).
- **Robustness:** The projection method gives the correct chi^2 distribution
  regardless of whether C_total is exactly singular, nearly singular, or
  technically full rank.

---

## VI. Practical implementation guide (~1 page)

### VI.A. Recipe

Step-by-step for experimentalists:

Step 0 — Predict N_null. Before computing any covariance, evaluate
N_null = N_bins - rank(A) from the MC event list. This provides a prediction
to cross-check against the eigendecomposition. Compute N_struct analytically
from the block topology as a lower bound.

Step 1 — Perform blockwise unfolding per Gardiner. Compute full N_bins ×
N_bins total covariance C_total (stat + syst).

Step 2 — Separately compute C_stat using the shared-event formula, propagated
through block-diagonal unfolding.

Step 3 — Eigendecompose C_stat. Identify N_null zero eigenvalues. Verify
against the prediction from Step 0.

Step 4 — Validate null eigenvectors. Structural ones should be uniform within
blocks (sum rules). Kinematic ones should correspond to empty cells in the
multi-block occupancy.

Step 5 — Construct P from the non-null eigenvectors. Compute C_tilde =
P C_total P and its pseudo-inverse C_tilde^+.

Step 6 — chi^2 = r^T C_tilde^+ r with N_bins - N_null d.o.f.

### VI.B. Predicting N_null

Lower bound from structural constraints:
  B blocks sharing all events, no slicing: N_struct = B - 1.
  S slices, slice s in B_s blocks: N_struct = sum_s (B_s - 1).

Full prediction: N_null = N_bins - rank(A), requires the MC event list.
rank(A) depends on the physics and binning — finer binning and stronger
kinematic correlations between variables lead to more empty cells and larger
N_kin.

### VI.C. What to include in data releases

Recommendations for experiments adopting blockwise unfolding:
- Full C_total (already standard per Gardiner's recommendation).
- The projector P (or equivalently the N_null null eigenvectors of C_stat).
- The value of N_null and N_bins - N_null as the stated d.o.f.
- Alternatively: the projected covariance C_tilde and its pseudo-inverse.
- If A_C is used: it is block-diagonal, does not affect inter-block structure,
  compatible with the projection.

### VI.D. Compatibility notes

- Works with any unfolding method within each block.
- Works with or without regularization matrix A_C.
- Compatible with MicroBooNE-style (analytic propagation) and MINERvA-style
  (universe re-extraction) uncertainty treatments.
- For MINERvA-style: use the projector from the CV C_stat. Do not recompute
  the projector in each universe.
- Safe to apply when N_null = 0: projector reduces to identity, standard
  chi-squared recovered. No cost to applying defensively.

---

## VII. Connection to the statistics literature (~0.5 page)

Frame the method in the language of established results:

- The projector P is the orthogonal projection onto the range (column space)
  of C_stat. Equivalently, P = C_stat C_stat^+.
- The pseudo-inverse C_tilde^+ is the Moore-Penrose generalized inverse.
- The projected chi^2 is the standard test for singular normal models,
  studied by Rao & Mitra (1971), Khatri (1968), and Moore (1977).
- The method is a physically motivated instance of the trimmed inverse
  framework of Davidov, Jelsema & Peddada (2018), with the trimming
  threshold determined by the known structural and kinematic constraints
  rather than by a numerical criterion.
- Related approaches in HEP: SVD cuts for near-singular covariance in
  lattice QCD (Bruno & Sommer 2022); rank-deficient covariance in
  unfolding (Blobel 2013).

Key references:
- Rao & Mitra, "Generalized Inverse of Matrices and Its Applications" (1971)
- Khatri, Sankhyā Ser. A 30, 267-280 (1968)
- Moore, JASA 72, 131-137 (1977)
- Davidov, Jelsema & Peddada, JASA 113, 906-918 (2018)
- Mathai & Provost, "Quadratic Forms in Random Variables" (1992)
- Silvapulle & Sen, "Constrained Statistical Inference" (2005), Ch. 3
- Ben-Israel & Greville, "Generalized Inverses" (2003)
- Bruno & Sommer, Comp. Phys. Comm. 285, 108631 (2023)
- Penrose, Proc. Cambridge Phil. Soc. 51, 406-413 (1955)

---

## VIII. Summary and recommendations (~0.5 page)

- Blockwise unfolding correctly captures inter-distribution correlations from
  shared events. A consequence is that the statistical covariance is
  generically singular.

- The null space has two components: structural constraints (predictable from
  the block topology) and kinematic constraints (from empty cells in the
  multi-dimensional phase-space occupancy).

- The total number of null directions is exactly given by
  N_null = N_bins - rank(A), where A is the bin-combination matrix. This is
  computable from the MC event list alone and verified against the eigenvalue
  spectrum of C_stat.

- Different classes of systematic uncertainties interact differently with the
  constraints: reweighting preserves all, detector variations can lift
  kinematic ones, and the factorized cross-section prescription can lift
  structural ones. In all cases, the eigenvalues in the formerly null
  directions are either zero or small, and including them in a chi-squared
  test produces inflated values.

- The projection method restricts the test to the subspace carrying
  independent statistical information, yielding a chi^2 with
  N_bins - N_null d.o.f.

- The method is validated on an analytical toy (exactly singular C_total)
  and a realistic generator-based example (with unfolding, smearing, and
  multi-source systematics).

- We recommend adoption of this method for reporting global goodness-of-fit
  metrics in all multi-distribution neutrino cross-section measurements.

- Application of these methods to a MicroBooNE multi-distribution
  cross-section measurement will be reported in a forthcoming publication.

---

## Appendices

### Appendix A: Null space structure and the combination matrix

Formal treatment for general block/slice topologies.

- Define the combination matrix A for arbitrary block structures including
  partial overlaps (events in some blocks but not others).
- Prove rank(C_stat) = rank(A).
- Derive the 2-block special case: bipartite graph, N_null = number of
  connected components.
- Prove that reweighting preserves rank(A): if the set of events and their
  bin assignments are unchanged, only the weights change, and positive weights
  cannot reduce the span of {phi_e}.
- Discuss the multi-block generalization (hypergraph connectivity).

### Appendix B: Distribution of the projected chi-squared

Self-contained proof that r^T (P C_total P)^+ r ~ chi^2(N_bins - N_null)
under H_0.

Two cases:
(i) C_total is singular (same null space as C_stat or larger). The residuals
are supported on the range of C_total, and the projected quadratic form
restricts to this subspace.

(ii) C_total is full rank but P projects onto a subspace. The projected
residuals Pr have covariance P C_total P with rank N_bins - N_null. The
quadratic form with the pseudo-inverse tests only the projected components.

State the condition for validity: P C_total P must be positive definite on
range(P). Prove this holds when C_syst is positive semi-definite on range(P),
which is generically true.

### Appendix C: Generator configuration and analysis details

(For reproducibility of the realistic example.)
- GENIE version, tune, target, flux specification.
- NuWro version and configuration.
- Smearing parametrization (functional forms, nominal values, variation ranges).
- Efficiency parametrization.
- Background model.
- Systematic universe specifications (number of universes, parameters varied,
  variation ranges).
- D'Agostini iteration count and convergence criteria.

---

## Figure summary

| #  | Content                                              | Section |
|----|------------------------------------------------------|---------|
| 1  | Eigenvalue spectrum of C_stat (analytical toy)       | IV.C    |
| 2  | Chi-squared distributions, 3 methods (toy) — MONEY   | IV.D    |
| 3  | p-value distributions (toy)                          | IV.D    |
| 4  | Occupancy matrices: full vs sparse (toy extension)   | IV.F    |
| 5  | Eigenvalue spectrum of C_stat (2-block realistic)    | V.B     |
| 6  | Eigenvalue spectra under cumulative syst addition    | V.B     |
| 7  | Chi-squared distributions (2-block realistic)        | V.B     |
| 8  | NuWro vs GENIE distributions (2-block)               | V.B     |
| 9  | Eigenvalue spectrum (4-block realistic)              | V.C     |
| 10 | Chi-squared distributions (4-block realistic)        | V.C     |

## Table summary

| #  | Content                                              | Section |
|----|------------------------------------------------------|---------|
| 1  | Chi-squared mean/variance vs theory (toy)            | IV.D    |
| 2  | Rank of C under cumulative syst addition             | V.B     |
| 3  | Chi-squared and p-values, NuWro vs GENIE             | V.B/C   |
| 4  | N_null: predicted rank(A) vs eigenvalue count        | V.D     |

---

## Estimated page count

| Section                                  | Pages |
|------------------------------------------|-------|
| I. Introduction                          | 1.0   |
| II. Origin of singularity                | 2.5   |
| III. Projection method                   | 2.0   |
| IV. Analytical demonstration             | 2.0   |
| V. Realistic demonstration               | 3.0   |
| VI. Practical implementation             | 1.0   |
| VII. Statistics literature               | 0.5   |
| VIII. Summary                            | 0.5   |
| Appendices A+B+C                         | 2.0   |
| References                               | 1.0   |
| **Total**                                |**~16**|

**To trim to ~12 pages:**
- Condense II.E (systematic hierarchy) into a tighter presentation
- Make the 4-block example (V.C) concise or move to appendix
- Merge VI (practical guide) into the summary
- Appendix C can be minimal (just parameter tables)
