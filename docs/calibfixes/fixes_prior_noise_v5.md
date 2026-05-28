# Fixes for Prior Noise and gDNA Handoff, v5

Date: 2026-05-27

This document supersedes `fixes_prior_noise_v4.md`. v4 removed the worst local hacks: no transcript-level RNA floor, no arbitrary `A_r` pseudocount, and no Normal-approximation strand threshold. v5 keeps those decisions and adds the two missing pieces needed for a robust implementation:

1. Capture exposure must be local and source-reliable, not a rigid global scalar or a latent-state ratio.
2. Strand source confidence must be continuous, not a hypothesis-test cliff.

v5 also states the remaining scope boundary clearly: the immediate calibration repair targets the three strata that already have usable non-strand or strand information. Unstranded capture-on is a separate source-identification problem that needs a later FL-plus-boundary model.

## 1. Core contracts

These invariants are not negotiable.

| contract | implication |
|---|---|
| Latent states are strata | `unexpressed_offtarget`, `unexpressed_capture`, `expressed_capture`, and `expressed_offtarget` describe expression/capture context, not RNA/gDNA source. |
| Source mass is explicit | EM source priors consume `PriorMassDeconvolution.gdna_unspliced_mean` and `.rna_unspliced_mean`, never latent-state source proxies. |
| Exposure is source-reliable local opportunity | `A_r` is an EM gDNA effective-length multiplier derived from source-reliable gDNA/background ratios; enriched regions have `A_r < 1`. |
| RNA prior is grouped | The adaptive prior can constrain aggregate RNA-vs-gDNA mass, but it must not add transcript-level pseudocounts. |
| Uncertainty must travel | Source reliability, exposure identifiability, and prior strength must be diagnosed and propagated instead of hidden behind thresholds. |

## 2. Scope by RNA-seq stratum

| stratum | v5 status | main evidence channel |
|---|---|---|
| Stranded, capture off | in scope | exact strand reliability plus density/background |
| Stranded, capture on | in scope | exact strand reliability plus one-class local capture exposure |
| Unstranded, capture off | in scope for non-regression | density/background and existing FL/global calibration |
| Unstranded, capture on | deferred | requires a dedicated FL-plus-boundary source split |

The fourth stratum is intentionally not solved by this document. In unstranded capture-on data, expressed captured exons can overwhelm density evidence, and strand evidence is absent by definition. Fragment length is useful, but it is not sufficient by itself: boundary-crossing fragments and inferred capture exposure are also important measurements. That work should be a later block, not a hidden assumption in the current repair.

## 3. What changed from v4

| v4 idea | v5 refinement |
|---|---|
| Single panel-level `lambda` for capture exposure | One-class local `lambda_u` from source-reliable gDNA/background ratios. |
| Exact predictive strand test with boolean reliable/unreliable result | Smooth reliability weight from the exact all-RNA predictive distribution or a marginal-likelihood Bayes factor. |
| ESS scaling as a policy sweep | Keep as a later policy sweep, gated by source reliability and exposure identifiability diagnostics. |
| Isoform collapse diagnostics | Defer implementation, but specify matrix-rank/identifiability diagnostics as the right frame. |
| Unstranded capture-on deferred | Keep deferred, but specify that the future model must combine FL, boundary crossing, and inferred capture exposure. |

## 4. Fix A v5: one-class source-reliable capture exposure

### Problem

v4's single global enrichment model,

```text
lambda > 1 means locally enriched gDNA sampling
```

was directionally right about local enrichment but wrong if wired as `A_r > 1`. The native EM subtracts `log(gdna_eff_len)`, and `prior.py` multiplies baseline gDNA effective length by `A_r`. Therefore enriched gDNA sampling must lower `A_r`, not raise it.

The immediate bug is more basic: current code derives `A_r` from latent-state-weighted `mu_gdna`. Latent capture states are not source labels. PR02 should first replace that with local exposure weights computed from PR01 source-reliable gDNA mass against the off-target background.

### Model unit

PR02 treats the panel as unknown. Fit exposure at a locus-level unit `u`, not at each fine region. Recommended units, in order of preference:

1. MultiLocus if transcripts share an EM/capture neighborhood;
2. isolated annotated gene;
3. no capture unit, which gives neutral exposure.

Fine regions are observations for a unit-level ratio. They are not independent exposure units.

### One-class model

Let:

```text
B_{r,c} = rho_off * L_{r,c}
G_{r,c} = source-reliable gDNA mass for region/compartment r,c
T_{r,c} = target_weight from annotation/locus structure
```

Aggregate to unit `u`:

```text
B_u = sum_{r,c in u} T_{r,c} * B_{r,c}
G_u = sum_{r,c in u} T_{r,c} * G_{r,c}
```

Then compute a one-class local exposure ratio:

```text
lambda_u_raw = G_u / B_u
lambda_u = max(lambda_u_raw, 1)  only when pooled source evidence beats background
```

Build region exposure:

```text
A_r = 1 / lambda_u
```

This works for unknown panels without seeding a captured class:

- non-capture samples have `G_u ~= B_u`, so `A_r ~= 1`;
- captured targets with source-reliable gDNA excess have `G_u > B_u`, so `A_r < 1`;
- RNA-rich regions without source-reliable gDNA evidence stay neutral;
- no arbitrary count constant and no maximum exposure cap are introduced.

### Source-reliable evidence only

`G_{r,c}` must come from PR01 source reliability, not from latent density fallback. Do not use:

- `p_captured`;
- `p_unexpressed_capture`;
- `p_expressed_capture`;
- `sweep.mu_sweep` by itself;
- `mu_gdna` when it was derived from latent-state mixing;
- raw observed count excess in RNA-rich exons.

If source-reliable evidence is absent, return neutral `A_r = 1`. This is why unstranded capture-on source splitting remains deferred.

### Boundary evidence

Boundary-crossing fragments are not optional. They are often the best evidence that a molecule is genomic rather than mature RNA, especially near captured exon targets. The local exposure fit should accept source-reliable evidence from contained and boundary compartments, each with its own effective-length/target-weight term:

```text
B_u = sum_compartments sum_r T_{r,c} * B_{r,c}
G_u = sum_compartments sum_r T_{r,c} * G_{r,c}
```

This keeps capture exposure tied to source-reliable gDNA/background ratios rather than to latent states or local observed read depth.

### Files

- `src/rigel/calibration/capture_exposure.py` for capture unit mapping, one-class exposure fitting, and `CaptureExposureFit`.
- `src/rigel/calibration/calibration_iteration.py` to stop deriving `A_r` from latent-state `mu_gdna`.
- `src/rigel/calibration/prior.py` to consume the new exposure field when computing gDNA effective lengths.
- Tests in a new `tests/test_capture_exposure.py`.

### Tests

1. No unit or no source-reliable evidence gives `A_r = 1` exactly.
2. High latent capture probability without PR01 source evidence does not inflate `A_r`.
3. Strong source-reliable local gDNA excess gives `lambda_u = G_u / B_u` and `A_r = 1 / lambda_u`.
4. Local depleted ratios are floored at neutral exposure, not below `1`.
5. Boundary source evidence contributes to `G_u` and `B_u` with the correct target weight.
6. Fine-region splitting does not change unit exposure after aggregation.

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
2. capture exposure has source-reliable local weights and identifiability diagnostics.

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
3. inferred capture exposure and source-reliable local enrichment.

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

### PR 2: one-class source-reliable capture exposure

- Build locus/gene capture units and target weights.
- Aggregate source-reliable gDNA mass and off-target background by unit.
- Fit one-class exposure ratios `lambda_u = G_u / B_u`, activated only by pooled source evidence over background.
- Build `A_r = 1 / lambda_u` for targetable regions and `A_r = 1` elsewhere.
- Emit exposure diagnostics and `inferred_capture_panel.bed`/TSV from the fitted ratios.

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
| Local exposure weights | capture exposure lowers gDNA effective length only from pooled source-reliable local gDNA/background ratios and stays neutral without source evidence. |

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
- High-gDNA stranded capture-on learns elevated local exposure where source-reliable locus evidence supports it.
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
