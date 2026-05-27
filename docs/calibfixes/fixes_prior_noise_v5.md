# Fixes for Prior Noise and gDNA Handoff, v5

Date: 2026-05-27

This document supersedes `fixes_prior_noise_v4.md`. v4 removed the worst local hacks: no transcript-level RNA floor, no arbitrary `A_r` pseudocount, and no Normal-approximation strand threshold. v5 keeps those decisions and adds the two missing pieces needed for a robust implementation:

1. Capture exposure must be local and empirical-Bayes, not a rigid global scalar.
2. Strand source confidence must be continuous, not a hypothesis-test cliff.

v5 also states the remaining scope boundary clearly: the immediate calibration repair targets the three strata that already have usable non-strand or strand information. Unstranded capture-on is a separate source-identification problem that needs a later FL-plus-boundary model.

## 1. Core contracts

These invariants are not negotiable.

| contract | implication |
|---|---|
| Latent states are strata | `unexpressed_offtarget`, `unexpressed_capture`, `expressed_capture`, and `expressed_offtarget` describe expression/capture context, not RNA/gDNA source. |
| Source mass is explicit | EM source priors consume `PriorMassDeconvolution.gdna_unspliced_mean` and `.rna_unspliced_mean`, never latent-state source proxies. |
| Exposure is geometry plus fitted enrichment | `A_r` is a sampling-opportunity multiplier, not a per-region abundance ratio. |
| RNA prior is grouped | The adaptive prior can constrain aggregate RNA-vs-gDNA mass, but it must not add transcript-level pseudocounts. |
| Uncertainty must travel | Source reliability, exposure identifiability, and prior strength must be diagnosed and propagated instead of hidden behind thresholds. |

## 2. Scope by RNA-seq stratum

| stratum | v5 status | main evidence channel |
|---|---|---|
| Stranded, capture off | in scope | exact strand reliability plus density/background |
| Stranded, capture on | in scope | exact strand reliability plus empirical-Bayes local capture exposure |
| Unstranded, capture off | in scope for non-regression | density/background and existing FL/global calibration |
| Unstranded, capture on | deferred | requires a dedicated FL-plus-boundary source split |

The fourth stratum is intentionally not solved by this document. In unstranded capture-on data, expressed captured exons can overwhelm density evidence, and strand evidence is absent by definition. Fragment length is useful, but it is not sufficient by itself: boundary-crossing fragments and capture geometry are also important measurements. That work should be a later block, not a hidden assumption in the current repair.

## 3. What changed from v4

| v4 idea | v5 refinement |
|---|---|
| Single panel-level `lambda` for capture exposure | Hierarchical local `lambda_u` with empirical-Bayes shrinkage toward the global panel mean. |
| Exact predictive strand test with boolean reliable/unreliable result | Smooth reliability weight from the exact all-RNA predictive distribution or a marginal-likelihood Bayes factor. |
| ESS scaling as a policy sweep | Keep as a later policy sweep, gated by source reliability and exposure identifiability diagnostics. |
| Isoform collapse diagnostics | Defer implementation, but specify matrix-rank/identifiability diagnostics as the right frame. |
| Unstranded capture-on deferred | Keep deferred, but specify that the future model must combine FL, boundary crossing, and capture geometry. |

## 4. Fix A v5: empirical-Bayes local capture exposure

### Problem

v4's single global enrichment model,

```text
A_r = 1 + c_r * (lambda - 1)
```

is principled but too rigid for real capture data. Capture efficiency varies by probe, target sequence, local composition, accessibility, and alignment context. A single global `lambda` can understate gDNA exposure at strongly captured probes and overstate it at weak probes.

The correct robustness mechanism is not a fixed cap or fixed pseudocount. It is local exposure with empirical-Bayes shrinkage.

### Model unit

Fit enrichment at an exposure unit `u`, not necessarily at each fine region. Recommended units, in order of preference:

1. probe group or merged BED12 probe target;
2. fine region if it has enough source-reliable evidence;
3. panel-wide fallback stratum.

Each fine region `r` maps to one or more exposure units through FL-aware capture opportunity `c_r` in `[0, 1]`.

### Hierarchical model

Let:

```text
e_r = rho_off * L_r
c_r = FL-aware capture opportunity
D_r = source-reliable gDNA mass for region/compartment r
theta_u = lambda_u - 1, theta_u >= 0
A_r = 1 + c_r * theta_u
```

Aggregate source-reliable evidence to unit `u`:

```text
O_u = sum_r w_r * e_r * c_r
X_u = sum_r w_r * max(D_r - e_r, 0)
```

where `w_r` is the continuous source-reliability weight from Section 5. `O_u` is the enrichment opportunity; `X_u` is the source-reliable excess gDNA mass above off-target expectation.

Use a Gamma-Poisson empirical-Bayes model:

```text
X_u | theta_u ~ Poisson(O_u * theta_u)
theta_u ~ Gamma(a_global, b_global)
E[theta_u] = a_global / b_global
```

Posterior mean:

```text
theta_u_post = (a_global + X_u) / (b_global + O_u)
lambda_u_post = 1 + theta_u_post
A_r = 1 + c_r * theta_u_post
```

This is the correct form of shrinkage:

- if local evidence is strong, `theta_u_post` follows the local unit;
- if local evidence is weak, it glides back to the global panel mean;
- if exposure is not identifiable, it falls back without a cliff;
- no arbitrary count constant and no maximum exposure cap are introduced.

### Estimating the global prior

Fit `a_global` and `b_global` from source-reliable units. Use method of moments or marginal likelihood over `{X_u, O_u}`. If too few reliable units exist, fall back to a broad global prior centered on the panel-wide enrichment estimate, and flag `capture_exposure_global_weak`.

The global prior is a learned property of the panel/run. It is not a magic number.

### Boundary evidence

Boundary-crossing fragments are not optional. They are often the best evidence that a molecule is genomic rather than mature RNA, especially near captured exon targets. The local exposure fit should accept evidence from contained and boundary compartments, each with its own effective-length/opportunity term:

```text
O_u = sum_compartments sum_r w_{r,c} * e_{r,c} * c_{r,c}
X_u = sum_compartments sum_r w_{r,c} * max(D_{r,c} - e_{r,c}, 0)
```

This keeps capture exposure tied to physical sampling opportunity rather than to a local observed/read-depth ratio.

### Files

- `src/rigel/calibration/capture_exposure.py` for `CaptureExposureFit`.
- `src/rigel/calibration/_exposure.py` for FL-aware capture opportunity helpers.
- `src/rigel/calibration/calibration_iteration.py` to stop deriving `A_r` from latent-state `mu_gdna`.
- `src/rigel/calibration/prior.py` to consume the new exposure field when computing gDNA effective lengths.
- Tests in a new `tests/test_capture_exposure.py`.

### Tests

1. No capture: all `c_r = 0` gives `A_r = 1` exactly.
2. No identifiable source-reliable excess: `theta_u_post` falls back to the global prior and flags weak local evidence.
3. Strong local enrichment: a unit with large `X_u` and `O_u` overrides the global mean smoothly.
4. Weak local enrichment: a unit with tiny `O_u` remains close to the global mean.
5. Boundary evidence contributes to `X_u` and `O_u` with the correct compartment opportunity.
6. No hidden caps: a synthetic unit can return very large `lambda_u_post` if the learned model and evidence support it.

## 5. Fix B v5: continuous exact strand reliability

### Problem

The exact all-RNA predictive distribution is correct, but a boolean test is still too brittle. A region should not jump from zero to full gDNA support because one read moves it across a threshold. Source reliability should glide from weak to strong evidence.

### Predictive distribution

For each strand-informative region or compartment:

```text
n = k_sense + k_antisense
k_minor = count in the RNA-minor orientation
q = RNA minor-orientation rate
```

Use the exact all-RNA predictive distribution:

```text
K_minor | all RNA ~ BetaBinomial(n, alpha_q, beta_q)
```

where `alpha_q` and `beta_q` come from the strand-model posterior. If that posterior is not yet carried through `StrandSummary`, extend `StrandSummary` to carry the underlying same/opposite counts or direct minor-rate posterior parameters.

### Smooth reliability weight

The reliability weight must be continuous and exact. Two acceptable designs are:

#### Preferred: marginal-likelihood reliability

Compare the all-RNA predictive model against the existing mixed RNA/gDNA strand model, marginalizing over gDNA count `D` with the same prior used by `deconvolve_regions_by_strand`:

```text
p0 = P(k_minor | n, all RNA)
p1 = P(k_minor | n, mixed RNA/gDNA)
BF = p1 / p0
w_strand = BF / (1 + BF)
```

This is smooth, exact at low `n`, and parameter-free apart from the model priors already present in the deconvolution layer.

#### Acceptable first implementation: exact-tail reliability

Use the exact predictive CDF or survival function as a monotone evidence score, but do not threshold it. The score must be tested carefully at low `n`, because raw CDF values are not automatically calibrated as posterior probabilities.

The production target should be the marginal-likelihood reliability above.

### Feeding reliability into source mass

Do not hard-drop a region. Weight its strand-derived source mass:

```text
G_r_strand_weighted = w_strand_r * G_r_strand_posterior
R_r_strand_weighted = unspliced_total_r - G_r_strand_weighted
```

Carry `w_strand_r` and the predictive evidence diagnostics into `PriorMassDeconvolution.flags` or a companion reliability array.

This is the key improvement over v4: source uncertainty remains continuous and visible instead of becoming a threshold artifact.

### Kappa interpretation

The current high-gDNA synthetic suite reports `kappa_d = 1e6`, which means the estimated gDNA strand balance is binomial or tighter than binomial. That may be valid for the current simulator, but it does not stress real overdispersion.

Do not force kappa lower as a fix. Instead:

- expose kappa diagnostics clearly;
- add simulator overdispersion settings later;
- let the exact reliability model carry uncertainty from both RNA minor-rate estimation and gDNA strand-balance estimation.

### Files

- `src/rigel/strand_model.py`: expose same/opposite counts or minor-rate posterior parameters.
- `src/rigel/calibration/strand_summary.py`: carry those parameters.
- `src/rigel/calibration/strand_deconv.py`: implement exact reliability helper.
- `src/rigel/calibration/calibration_iteration.py`: propagate reliability into source mass and diagnostics.
- `tests/test_strand_deconv.py`.

### Tests

1. Low-depth exactness: `n=10`, `q=0.01` uses exact beta-binomial/binomial probabilities, not a Normal approximation.
2. Continuity: adding one minor read changes `w_strand` smoothly, not from 0 to 1.
3. Pure RNA: typical all-RNA observations produce low reliability.
4. Strong mixed source: deep upper-tail observations produce reliability near 1.
5. Small strand-training set: uncertainty in `q` widens the predictive null and lowers reliability.
6. Near-unstranded data: reliability is inactive or low when strand contrast is below the existing numerical floor.

## 6. ESS scaling after reliability, not before

The current `MAX_ESS = 3000` may be too small for very large loci, but linear `0.25 * n_em` can become a 25k-count prior in mega-loci. v5 keeps ESS scaling as a controlled policy sweep, not a first-order fix.

Only evaluate a larger ESS after both are true:

1. source mass has continuous reliability weights;
2. capture exposure has empirical-Bayes local enrichment and identifiability diagnostics.

The preferred final cap should be information-weighted:

```text
ESS_cap = f(n_competing_fragments, source_reliability, exposure_identifiability)
```

where `n_competing_fragments` counts fragments with finite competing RNA and gDNA likelihoods, not just all locus fragments.

Benchmark policies:

| policy | formula | status |
|---|---|---|
| current | `min(locus_unspliced, 3000)` | baseline |
| linear | `min(locus_unspliced, 0.25 * n_em)` | stress-test only |
| information-weighted | learned reliability-weighted cap | target |

Do not tune `0.25` if it regresses RNA MARD. If linear scaling is too strong, that is evidence for the information-weighted cap, not for another constant.

## 7. Deferred: unstranded capture-on source split

Unstranded capture-on remains the only stratum not solved by this plan. It should be a separate development block.

The future model should combine at least three signals:

1. fragment-length likelihood under RNA FL vs captured-gDNA FL;
2. boundary-crossing and intronic compatibility evidence;
3. capture geometry and local empirical-Bayes exposure.

The FL part is important but not sufficient alone. Boundary-crossing fragments are direct evidence for genomic molecules and must be part of the model.

Future sketch:

```text
P(fragment | source = RNA)  uses RNA FL + transcript compatibility
P(fragment | source = gDNA) uses captured-gDNA FL + genomic/boundary compatibility
```

Fit a captured-gDNA FL model from source-reliable off-target introns and boundary fragments that experienced the capture protocol. Then use a regional mixture to produce strand-free `G_r` for `PriorMassDeconvolution`.

This work is deferred because the current priority is to repair the glaring stranded capture and source-handoff failures without mixing in another source channel.

## 8. Deferred: isoform identifiability diagnostics

GENE0008.3 dropping to zero should not be fixed with an RNA floor. It is likely an identifiability issue exposed when aggregate RNA mass shrinks under stronger gDNA priors.

Later, add an equivalence-class matrix diagnostic:

```text
M[u, t] = compatibility or likelihood contribution of EM unit u to transcript t
S = M.T @ W @ M
```

Report:

- effective rank of `S`;
- condition number or smallest singular value;
- near-duplicate transcript columns;
- whether zeroed transcripts have unique evidence;
- `isoform_ambiguous = True` when the matrix is rank-deficient or nearly so.

This is important scientifically, but it is not a calibration-prior repair and should not block v5's source/exposure work.

## 9. Implementation sequence

### PR 1: continuous strand reliability

- Extend `StrandSummary` with minor-rate posterior parameters.
- Add exact beta-binomial predictive reliability in `strand_deconv.py`.
- Weight strand-derived source mass continuously.
- Emit reliability diagnostics.

### PR 2: empirical-Bayes local capture exposure

- Add FL-aware capture opportunity per region/compartment.
- Fit global and local `theta_u = lambda_u - 1` with Gamma-Poisson EB.
- Build `A_r = 1 + c_r * theta_u_post`.
- Emit local/global exposure diagnostics.

### PR 3: ESS policy sweep

- Add diagnostics for source reliability and exposure identifiability at the locus level.
- Benchmark current, linear, and information-weighted ESS policies.
- Merge only a policy that preserves all non-target strata.

### PR 4: simulator overdispersion knobs

- Add gDNA strand-balance overdispersion.
- Add RNA minor-orientation overdispersion.
- Add optional region-level strand-bias heterogeneity.
- Rerun the suite under near-binomial and moderately overdispersed settings.

Deferred PRs:

- unstranded capture-on FL-plus-boundary source split;
- isoform matrix-rank diagnostics.

## 10. Acceptance criteria

### Invariants

| invariant | gate |
|---|---|
| No transcript floor | no code path adds prior-scaled mass to every transcript. |
| No latent source labels | adaptive prior and exposure code never uses latent state names as RNA/gDNA source labels. |
| No arbitrary exposure constants | no fixed `A_R_SOURCE_PRIOR_COUNT` or `A_R_MAX` replacement appears in the exposure model. |
| Smooth source reliability | strand source confidence is continuous and exact at low depth. |
| Local exposure shrinkage | capture exposure can vary locally but shrinks to a learned global panel prior. |

### Existing eight-condition suite

| scenario group | gate |
|---|---|
| no-gDNA scenarios | estimated gDNA remains at current tolerance and no false capture exposure inflation appears. |
| high-gDNA SS=0.99 capture-on | gDNA fraction moves toward 0.5 and mRNA MARD improves. |
| high-gDNA SS=0.99 capture-off | mRNA MARD does not regress beyond 5 percent relative. |
| high-gDNA SS=0.50 capture-off | mRNA MARD does not regress beyond 5 percent relative. |
| high-gDNA SS=0.50 capture-on | no material regression; full repair deferred to FL-plus-boundary source split. |

### Diagnostic gates

- No-gDNA capture-on has low source reliability and does not learn large local exposure from residual source noise.
- High-gDNA stranded capture-on learns elevated local exposure where probe evidence supports it.
- Strand reliability changes smoothly under one-read perturbations in low-depth regions.
- ESS changes, if enabled, do not create new transcript-collapse patterns.

## 11. Validation commands

For documentation and diagnostics:

```bash
conda activate rigel && ruff check scripts/debug/review_prior_noise_v4.py
```

For PR 1:

```bash
conda activate rigel && ruff check src/rigel/strand_model.py src/rigel/calibration tests/test_strand_deconv.py
conda activate rigel && pytest tests/test_strand_deconv.py tests/test_calibration_iteration.py -v
```

For PR 2:

```bash
conda activate rigel && ruff check src/rigel/calibration tests/test_capture_exposure.py
conda activate rigel && pytest tests/test_capture_exposure.py tests/test_per_locus_gdna_mass.py -v
```

After PR 1 and PR 2, rerun the eight-condition hybrid-capture suite before any ESS policy change.
