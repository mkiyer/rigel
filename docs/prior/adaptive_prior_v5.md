# Adaptive Prior v5 — Entropy-Weighted Dirichlet Counts

Date: 2026-05-26
Status: implementation-ready design
Supersedes: `docs/prior/adaptive_prior_v4.md`, `adaptive_prior_v3.md`,
`adaptive_prior_v2.md`, `adaptive_prior_v1.md`, `prior_redesign_v3.md`, and all
earlier prior-design notes. No backwards compatibility is retained.

## 1. Why a v5

v4 was already simple in spirit: one knob, one quantile, one ESS cap. A close
review showed three remaining brittlenesses:

1. The "one user knob" still hid internal calibration confidence thresholds
   that effectively re-introduce hand tuning under a different name.
2. Translating regional evidence into a per-locus Beta via moment matching is
   unstable at the very boundaries that matter most (true zero gDNA, all-gDNA
   capture). Shape parameters below 1 produce U-shaped Betas, and a fixed
   quantile of a U-shaped Beta flips violently with calibration noise.
3. An explicit leave-one-locus-out global shrinkage loop is heavyweight and
   easy to misindex.

v5 removes all three by treating calibration evidence as discrete pseudo-counts
in a Dirichlet–Multinomial framework, scaling each region's contribution by
the information content of its four-state posterior, and replacing LOO with a
single vectorized "global minus local" subtraction.

The whole construction reduces to:

```text
per region r:
    counts_r   = U_r * (1 - H(P_r) / H_max)  * P_r         # 4-vector
    -> collapse to 2 groups (gDNA, RNA) via state mask

per locus l:
    N_l        = sum_r s_lr * counts_r                     # 2-vector

per sample:
    N_global   = sum_l N_l                                 # 2-vector

per locus prior (Dirichlet posterior, group-collapsed):
    alpha_l    = N_l + (N_global - N_l) * shrink_weight    # 2-vector
    cap        = min(U_l, MAX_ESS)
    if sum(alpha_l) > cap:
        alpha_l *= cap / sum(alpha_l)

    alpha_gdna_add_l = alpha_l[0]
    alpha_rna_add_l  = alpha_l[1]
```

No quantile. No Beta inversion. No LOO loop. No user-facing knobs other than
the structural gate already enforced by native EM.

## 2. User Contract

### 2.1 No Prior Knobs

The grouped EM prior has **zero** user-facing tunable parameters. The
calibration model, the four-state tensor, and the structural gate are the
entire interface.

The following are removed from CLI, YAML, and `EMConfig`:

```text
rna_confidence
rna_lower_confidence
gdna_density_confidence
aggregate_prior_strength
aggregate_prior_edge_count
aggregate_prior_max_count
gdna_prior_logit_bias
```

Setting any of them in a YAML config raises a migration error pointing at this
document. No deprecation aliases. No legacy modes. No `--rna-confidence`.

### 2.2 What Replaces Them

Three things, all internal:

- The four-state calibration posterior `P_r` already produced by
  `calibration_iteration.py`.
- The unspliced prior-owned mass `U_r` already produced by
  `PriorMassDeconvolution.unspliced_total`.
- The native structural flag `has_gdna_candidate_l`.

If a user wants to override the model, they edit calibration code, not run
flags.

## 3. Core Objects

### 3.1 Four-State Posterior

The calibration tensor already provides, per region `r`, a length-4 simplex
vector:

```text
P_r = [p_background, p_gdna_only_capture, p_expressed_capture, p_expressed_offtarget]
sum(P_r) = 1
```

Group masks (existing constants in `latent_states.py`):

```text
GDNA_STATES = {STATE_BACKGROUND, STATE_GDNA_ONLY_CAPTURE}
RNA_STATES  = {STATE_EXPRESSED_CAPTURE, STATE_EXPRESSED_OFFTARGET}
```

The per-region group probability vector is the deterministic projection:

```text
q_r = [ P_r[GDNA_STATES].sum(), P_r[RNA_STATES].sum() ]
```

This collapse is exact under the calibration model: every unspliced-owned
fragment in a region belongs to one of the four states, and the prior groups
gDNA against aggregate RNA.

### 3.2 Information Weight

Shannon entropy of the four-state posterior on each region:

```text
H(P_r)  = - sum_s P_r[s] * log P_r[s]        (with 0 log 0 := 0)
H_max   = log(N_STATES) = log 4
w_r     = 1 - H(P_r) / H_max                 in [0, 1]
```

`w_r` is 1 when calibration is certain about the four-state assignment, 0 when
it is maximally ambiguous. No threshold, no cutoff, no tuning.

### 3.3 Regional Pseudo-Counts

The discrete evidence pseudo-count vector emitted by region `r` is:

```text
n_r = U_r * w_r * q_r                          # length-2 (gDNA, RNA)
```

Properties:

- `sum(n_r) = U_r * w_r <= U_r`, so a region never contributes more than its
  observed unspliced mass.
- Maximally certain regions contribute their full mass.
- Maximally ambiguous regions contribute nothing.
- `U_r == 0` contributes nothing.
- `q_r` is already a normalized two-vector, so `n_r` lives on the correct
  count scale.

## 4. Locus And Sample Aggregation

### 4.1 Geometry Projection

Use the existing region-to-locus overlap shares `s_lr` from
`_allocate_counts_by_geometry()`. For each locus `l`:

```text
N_l = sum_r s_lr * n_r                         # length-2
U_l = sum_r s_lr * U_r                         # scalar
```

The same `s_lr` already allocates means and unspliced totals, so v5 reuses
that allocator unchanged. The geometry-allocated `N_l` is mass-conserving in
the limit of certain regions because `n_r <= U_r * w_r` and `sum_r s_lr <= 1`.

### 4.2 Sample-Wide Aggregate

```text
N_global = sum_l N_l                           # length-2
```

This is computed once per calibration run. It is the entire empirical-Bayes
state.

## 5. Per-Locus Prior

### 5.1 Global Minus Local

Replace LOO with one subtraction:

```text
N_other_l = N_global - N_l                     # length-2 per locus
```

`N_other_l` is the cross-sample background that every other locus contributes,
guaranteed nonnegative because `N_l <= N_global` componentwise.

### 5.2 Local Confidence

The locus's own information weight is the mass-weighted region average:

```text
W_l = sum_r s_lr * U_r * w_r
      / max(sum_r s_lr * U_r, EPS_MASS)
```

`W_l` lies in `[0, 1]`. It is the locus-level analog of `w_r`. A locus made
entirely of confident regions has `W_l = 1`; a locus made of ambiguous
regions has `W_l = 0`.

### 5.3 Shrinkage Weight

The global contribution is reduced when the locus already has confident local
evidence:

```text
shrink_l = 1 - W_l                             in [0, 1]
```

A locus that is locally certain ignores the global pool; a locus that is
locally ambiguous accepts the global pool fully.

### 5.4 Combined Counts

```text
alpha_l = N_l + shrink_l * N_other_l           # length-2
```

This is the Dirichlet posterior pseudocount vector for the
{gDNA, aggregate RNA} simplex at locus `l`.

### 5.5 Cap And Structural Gate

```text
cap_l = min(U_l, MAX_ESS)
if sum(alpha_l) > cap_l:
    alpha_l *= cap_l / sum(alpha_l)

if not has_gdna_candidate_l:
    alpha_l = (0, 0)

if U_l <= EPS_MASS:
    alpha_l = (0, 0)

alpha_gdna_add_l = alpha_l[0]
alpha_rna_add_l  = alpha_l[1]
```

`MAX_ESS` is the single numerical safety cap. It is private to
`adaptive_prior.py`.

That is the entire prior. There is no quantile inversion, no Beta, no
moment matching, no per-method trust constant, no logit bias.

## 6. Why This Works At The Boundaries

All hard cases collapse to simple arithmetic:

- **True zero gDNA, well-resolved regions.** `q_r` concentrates on RNA states,
  `w_r` near 1, `N_l` concentrates on RNA. `alpha_l` is a pure RNA prior. No
  Beta required.

- **True all-gDNA capture region.** `q_r` concentrates on gDNA states, `w_r`
  near 1. `alpha_l` is a pure gDNA prior, capped by `U_l`.

- **Ambiguous unstranded unspliced locus.** `w_r` small, `n_r` small, `N_l`
  small. `shrink_l` close to 1, so the locus shrinks to the sample profile.
  No artificial siphoning, no false certainty.

- **Sample with one locus.** `N_other_l = 0`, `alpha_l = N_l`. The prior is
  exactly the locus's own evidence.

- **Sample with totally homogeneous gDNA share.** `N_global` matches every
  `N_l` in direction. The shrinkage to the global pool reinforces the local
  direction without exceeding `U_l`.

- **Sample with heterogeneous capture.** `N_global` is a mixed-direction
  vector, but `shrink_l` is small for confident loci, so capture-confident
  loci are unaffected; only ambiguous loci absorb the sample profile.

In no case does a calibration-noise fluctuation produce a prior swing of the
kind a Beta quantile would. The mapping from evidence to prior is linear.

## 7. Implementation

### 7.1 New Module

```text
src/rigel/calibration/adaptive_prior.py
```

API:

```python
@dataclass(frozen=True, slots=True)
class AdaptivePriorResult:
    alpha_gdna_add:    np.ndarray   # float64[L]
    alpha_rna_add:     np.ndarray   # float64[L]
    n_local:           np.ndarray   # float64[L, 2]
    n_other:           np.ndarray   # float64[L, 2]
    locus_weight:      np.ndarray   # float64[L] = W_l
    shrink_weight:     np.ndarray   # float64[L] = 1 - W_l
    global_counts:     np.ndarray   # float64[2]
    region_weight:     np.ndarray   # float64[R] = w_r
    flags:             np.ndarray   # uint16[L]


def compute_adaptive_prior(
    *,
    region_arrays: RegionArrays,
    multi_loci: list[MultiLocus],
    p_states: np.ndarray,            # float64[R, 4]
    unspliced_total: np.ndarray,     # float64[R]
    has_gdna_candidate: np.ndarray,  # bool[L]
    *,
    max_ess: float = MAX_ESS,
) -> AdaptivePriorResult:
    ...
```

The function does everything: entropy weights, regional counts, geometry
projection, global aggregation, shrinkage, cap, gating.

Private constants:

```python
MAX_ESS  = 3000.0
EPS_MASS = 1.0e-12
```

There are no other constants.

### 7.2 Removed Code

Delete:

- `EMConfig.aggregate_prior_strength`,
  `EMConfig.aggregate_prior_edge_count`,
  `EMConfig.aggregate_prior_max_count`,
  `EMConfig.gdna_prior_logit_bias` and their validation.
- `CalibrationConfig.rna_lower_confidence`,
  `CalibrationConfig.gdna_density_confidence` and all references to them in
  the orchestrator, strand deconvolution, and density evidence paths.
- `GroupedPriorCounts`, `_compute_grouped_prior_counts()`,
  `prior_budget*`, `prior_gdna_share_*biased`, and `gdna_prior_logit_bias`
  references in `calibration/prior.py`, `pipeline.py`, `estimator.py`,
  `cli.py`, and tests.
- The `rna_lower_confidence` / `gdna_density_confidence` CLI flags and their
  `_ParamSpec` entries.
- `_ParamSpec` registry entries for any removed config field.

Where strand deconvolution and density evidence currently consume a
`rna_lower_confidence` or `gdna_density_confidence` argument, replace it with
a fixed internal constant declared at the top of the module, named to make
its non-tunable nature obvious (for example, `_INTERNAL_RNA_LOWER_CI = 0.95`).
These constants are not exported, not configurable, and not part of the prior
algorithm — they only govern the unrelated RNA lower-bound screen and gDNA
upper diagnostic, which v5 does not consume for the prior.

### 7.3 Reused Code

- `_allocate_counts_by_geometry()` is reused unchanged for `s_lr` allocation
  of both `U_r` and `n_r`.
- `enable_gdna_for_multilocus()` continues to produce `has_gdna_candidate_l`.
- The native EM grouped-prior ABI is unchanged. It still takes
  `(alpha_gdna_add, alpha_rna_add, gdna_eff_len, enable_gdna)` per locus.

### 7.4 Wiring

`assemble_priors()` becomes a thin orchestrator:

```python
def assemble_priors(*, multi_loci, em_data, index, calibration) -> PriorTable:
    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    rc = calibration.region_calibration
    has_gdna_candidate = np.array(
        [enable_gdna_for_multilocus(l, em_data) for l in multi_loci], dtype=bool
    )
    result = compute_adaptive_prior(
        region_arrays=region_arrays,
        multi_loci=multi_loci,
        p_states=np.asarray(rc.p_states, dtype=np.float64),
        unspliced_total=np.asarray(rc.prior_mass.unspliced_total, dtype=np.float64),
        has_gdna_candidate=has_gdna_candidate,
    )
    # gdna_eff_len, exposure weighting, locus-span diagnostics handled here
    return PriorTable(...)
```

No `EMConfig` is passed in. The prior has nothing to configure.

## 8. Output Schema

### 8.1 `loci.feather` columns

Replace all v3/v4 prior columns with this minimal set:

```text
alpha_gdna_add              # final gDNA grouped pseudocount
alpha_rna_add               # final aggregate-RNA grouped pseudocount
prior_locus_weight          # W_l
prior_shrink_weight         # 1 - W_l
prior_n_local_gdna          # N_l[0]
prior_n_local_rna           # N_l[1]
prior_n_other_gdna          # N_other_l[0]
prior_n_other_rna           # N_other_l[1]
prior_ess_final             # sum(alpha_l) after cap and gate
prior_flags                 # bitfield
```

Removed columns (no aliases, no compatibility):

```text
gdna_prior_count, gdna_prior_count_em, rna_expected_count,
prior_unspliced_total, prior_budget_raw, prior_budget,
prior_gdna_share_raw, prior_gdna_share_biased, gdna_prior_density,
prior_gdna_share_mean, prior_gdna_share_var, prior_gdna_share_quantile,
prior_ess_local, prior_ess_global, prior_global_share, prior_global_agreement,
prior_conflict_score
```

Goldens are regenerated in a single pass. There is no compatibility window.

### 8.2 `summary.json`

```json
"prior_policy": {
  "name": "entropy_dirichlet_v5",
  "max_ess": 3000.0,
  "global_counts": {"gdna": 1.2e4, "rna": 1.7e5},
  "region_weight": {"mean": 0.78, "p10": 0.21, "p90": 0.97},
  "locus_weight":  {"mean": 0.81, "p10": 0.28, "p90": 0.99},
  "n_loci_total": 12345,
  "n_loci_with_prior_mass": 11890,
  "n_loci_structural_gated": 410,
  "flag_histogram": {"PRIOR_NO_UNSPLICED_MASS": 213}
}
```

### 8.3 Flags

Single bitfield, kept short on purpose:

```text
0x1 PRIOR_NO_UNSPLICED_MASS     U_l <= EPS_MASS
0x2 PRIOR_STRUCTURAL_GATED      has_gdna_candidate_l == False
0x4 PRIOR_ESS_CAPPED            cap_l < sum(alpha_l) before scaling
```

No "no variance" flag, no "conflict" flag, no "global empty" flag. There is no
variance and no global-empty case in v5.

## 9. Tests

### 9.1 Helper Unit Tests (`tests/test_adaptive_prior.py`)

1. `H(P_r)`: certain region (one-hot) has `w_r = 1`; uniform region has
   `w_r = 0`; intermediate is monotone.
2. `n_r`: `U_r == 0` ⇒ `n_r = 0`; `w_r == 0` ⇒ `n_r = 0`; certain pure-RNA
   region ⇒ `n_r = (0, U_r)`; certain pure-gDNA region ⇒ `n_r = (U_r, 0)`.
3. Geometry: disjoint regions sum to per-locus counts equal to per-region
   counts; one region split across two loci splits its counts proportionally
   to `s_lr`.
4. `N_global - N_l >= 0` componentwise for every locus.
5. Locus with `W_l = 1` and nonzero `N_l` ignores any nonzero `N_other`.
6. Locus with `W_l = 0` and `N_other > 0` adopts the full `N_other` direction.
7. Cap: when `sum(alpha_l) > U_l`, `sum(alpha_l)` is rescaled to exactly
   `U_l` and the share direction is preserved.
8. Structural gate: `has_gdna_candidate_l == False` ⇒ both alphas exactly 0.
9. Single-locus sample: `N_other_l = 0` ⇒ `alpha_l = N_l`.
10. Vectorization: shuffling locus order does not change per-locus outputs.

### 9.2 Integration Tests

- Zero-gDNA whole-transcriptome sentinels recover RNA with `alpha_rna_add > 0`
  and `alpha_gdna_add == 0` even though `has_gdna_candidate_l == True`.
- True high-gDNA sentinels carry `alpha_gdna_add > 0` and recover gDNA in EM.
- Unstranded nRNA-heavy sentinels show low `prior_ess_final` rather than a
  large false `alpha_gdna_add`.
- Capture sentinels with cross-locus disagreement: locally confident loci have
  `prior_shrink_weight` near 0 and resist the global pool.

### 9.3 Acceptance Sweep

Run the four-regime sweep:

```text
stranded   whole-transcriptome
unstranded whole-transcriptome
stranded   capture
unstranded capture
```

Each sweep must pass without any code change, config flag, YAML key, or
threshold adjustment between regimes.

## 10. Migration

No backwards compatibility. The PR sequence is:

1. Add `adaptive_prior.py` with the new function and unit tests.
2. Rewire `assemble_priors()` to call it.
3. Delete the v3/v4 grouped-prior code, `EMConfig.aggregate_prior_*`,
   `EMConfig.gdna_prior_logit_bias`, `rna_lower_confidence`,
   `gdna_density_confidence`, related CLI flags, `_ParamSpec` entries,
   YAML keys, and all associated tests.
4. Regenerate `loci.feather` and `summary.json` goldens.
5. Update docs: this file is the source of truth; remove or archive
   `prior_redesign_v3.md`, `adaptive_prior_v1.md`,
   `adaptive_prior_v2.md`, `adaptive_prior_v3.md`,
   `adaptive_prior_v4.md`, and `one_knob_parameterization.md` under
   `docs/prior/archive/`.

YAML migration error message:

```text
Configuration option '{name}' was removed in adaptive prior v5.
The grouped EM prior is now parameter-free.
See docs/prior/adaptive_prior_v5.md.
```

## 11. Non-Goals

- Splitting annotated mRNA and nRNA into separate prior groups.
- Auto-tuning `MAX_ESS` from data.
- Replacing the EM solver or its grouped-prior interface.
- Re-litigating the four-state calibration tensor itself. v5 consumes it as a
  contract; improvements to the tensor are separate work.
- Operational sensitivity/specificity risk dial: see
  [docs/prior/adaptive_prior_v6.md](adaptive_prior_v6.md). v6 is an additive
  extension that introduces a single user-facing dial without changing any v5
  internals.

## 12. Summary

v5 is the smallest design that satisfies the original goal: a robust,
self-tuning grouped EM prior with no user-facing parameters.

The only algorithmic ingredients are:

- the four-state calibration posterior already produced by the v6 pipeline;
- the unspliced prior-owned mass already produced by `PriorMassDeconvolution`;
- Shannon entropy as the universal "how much do I trust this region" weight;
- a single global pseudocount aggregate and one per-locus subtraction;
- one safety cap `MAX_ESS`.

Everything else — Beta moment matching, quantile inversion, posterior
variance propagation, logit biases, edge counts, strength multipliers, trust
constants, LOO loops — is deleted.
