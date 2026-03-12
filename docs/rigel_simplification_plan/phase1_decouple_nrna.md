# Phase 1: Decouple nRNA from mRNA — Eliminate the η Parameterization

## Summary

Remove the hierarchical M-step that links nRNA components to their "child"
mRNA transcripts via a Beta(α, β) prior on the nascent fraction η. Make nRNA a
regular independent EM component that competes with mRNA through the
likelihood, just like any other component. Fix the nRNA component gating that
kills nRNA at low strand specificity.

**Root causes addressed**: RC1, RC2, RC3 (see
[complex_locus_ablation_studies_v1.md](../complex_locus_ablation_studies_v1.md))

**Parameters removed**: 9 of 27 EMConfig fields → 18 remaining

---

## Rationale

The current system links each nRNA span to its overlapping mRNA transcripts
through a "nascent fraction" η = nRNA / (nRNA + mRNA). This η is:

1. **Estimated pre-EM** by a 3-tier empirical Bayes hierarchy (global →
   locus-strand → per-nRNA) using density+strand hybrid estimators and MoM κ
2. **Converted** to Beta(α, β) prior parameters: α = η × κ_nrna, β = (1−η) × κ_nrna
3. **Enforced during EM** via a hierarchical M-step that redistributes mass
   between each nRNA and its sharing mRNAs according to the MAP Beta posterior

This is unnecessary because:

- **nRNA spans already compete with mRNA** for the same fragments through the
  EM likelihood. Intronic reads are unambiguously nRNA. Exonic reads compete
  between mRNA and nRNA based on effective lengths. The EM handles this
  naturally.
- **η varies enormously across genes**: from ~0 (stable mRNA, slow
  transcription) to ~1 (rapid transcription, unstable mRNA). A shared prior
  is biologically inappropriate.
- **The EB shrinkage actively hurts**: at moderate nRNA levels (200 frags),
  the 3-tier hierarchy over-shrinks, causing -14% to -83% nRNA
  underestimation. Heavy nRNA (1000 frags) recovers better simply because it
  has enough evidence to overcome the prior bias.
- **The Beta MAP prevents η from reaching 0 or 1**: creating +864% mRNA false
  positives when only nRNA is expressed (pattern 15).
- **The initialization gating kills nRNA**: at SS ≤ 0.6, `compute_nrna_init()`
  zeros all nRNA, and `build_locus_em_data()` permanently kills the
  components. 10K nRNA fragments get dumped into gDNA.

---

## Changes by File

### 1. `src/rigel/native/em_solver.cpp` — Remove Hierarchical M-Step

#### 1a. `hierarchical_map_em_step()` (~L540–753)

**Current**: Two branches — `n_nrna == 0` (plain MAP-EM) and `n_nrna > 0`
(hierarchical M-step with Beta posterior on η).

**Change**: Delete the `n_nrna > 0` branch entirely (lines ~L660–707). The
function always takes the plain MAP-EM path. Since `n_nrna` is no longer
needed, remove these parameters from the function signature:

- `n_transcripts`
- `n_nrna`
- `t_to_nrna` (const int32_t*)
- `nrna_to_t_off` (const int32_t*)
- `nrna_to_t_idx` (const int32_t*)
- `nrna_frac_alpha` (const double*)
- `nrna_frac_beta` (const double*)

Rename the function to `map_em_step()` (the "hierarchical" qualifier no
longer applies).

**Preserve**: The strand symmetry penalty for gDNA (this is independent of
the nRNA hierarchy and remains in the plain MAP-EM path).

#### 1b. `run_squarem()` (~L876–1170)

**Change**: Remove the nRNA-related parameters from the function signature.
Update all calls to `hierarchical_map_em_step()` (→ `map_em_step()`) to drop
the nRNA args (at ~L1014, L1023, L1063, L1149).

#### 1c. `LocusSubProblem` struct (~L1347–1386)

**Change**: Remove these fields:

- `nrna_frac_alpha` (std::vector<double>)
- `nrna_frac_beta` (std::vector<double>)

**Keep** (needed for posterior assignment and output mapping):

- `local_t_to_local_nrna`
- `nrna_to_t_off`, `nrna_to_t_idx`
- `local_to_global_nrna`

#### 1d. `extract_locus_sub_problem()` (~L1388–1700)

**Change**:
- Remove `all_nrna_frac_alpha` and `all_nrna_frac_beta` from the function
  signature
- Delete the α/β copy block (~L1688–1693)
- **Fix nRNA gating** (see Section 2 below): modify the `nrna_init == 0`
  gating logic

#### 1e. `run_locus_em_native()` (~L1180–1250)

**Change**: Remove from the function signature:
- `nrna_frac_alpha_arr`
- `nrna_frac_beta_arr`

Remove the `hierarchical` flag logic (~L1222–1223). Always pass
`n_transcripts=0`, `n_nrna=0`, nullptrs for the nRNA arrays to `run_squarem`.

Alternatively, since all calls now go to the plain branch, simply stop passing
nRNA hierarchy args to `run_squarem()`.

#### 1f. `batch_locus_em()` (~L1855–1920)

**Change**: Remove from the function signature:
- `nrna_frac_alpha_arr`
- `nrna_frac_beta_arr`

Remove the raw pointer extraction (~L1940–1941). Update the
`extract_locus_sub_problem()` call to drop α/β args. Update `run_squarem()`
calls to drop nRNA args.

#### 1g. Nanobind Bindings (~L2520–2600)

**Change**:
- Remove `nb::arg("nrna_frac_alpha")` and `nb::arg("nrna_frac_beta")` from
  both `run_locus_em_native` and `batch_locus_em` bindings
- Remove `NRNA_FRAC_CLAMP_EPS` constant export

#### 1h. Constant

**Change**: Remove `NRNA_FRAC_CLAMP_EPS` constant (~L43).

---

### 2. `src/rigel/native/em_solver.cpp` — Fix nRNA Component Gating

#### Current gating in `extract_locus_sub_problem()` (~L1667–1675)

```cpp
// Zero nRNA prior when nrna_init is zero for that nRNA
for (int n = 0; n < n_nrna; ++n) {
    int32_t gnrna = unique_nrna_global[n];
    if (all_nrna_init[gnrna] == 0.0) {
        sub.prior[n_t + n] = 0.0;
    }
}
```

**Problem**: When `compute_nrna_init()` returns zeros (SS ≤ 0.6, or no
intronic span), `nrna_init == 0` → prior set to 0 → component permanently
dead → all nRNA mass goes to gDNA or mRNA.

**Change**: Replace the hard gate with a soft floor. When `nrna_init == 0`
but the nRNA has intronic span (i.e., it's a multi-exon gene nRNA), keep the
prior at `EM_PRIOR_EPSILON` rather than zeroing it. The EM can still shrink
the component to near-zero if warranted by the data.

Specifically:
- **Keep zeroing** nRNA components where ALL sharing transcripts are
  single-exon (no intronic span to distinguish nRNA from gDNA — this is
  inherent structural ambiguity, and the prior zeroing is correct)
- **Stop zeroing** nRNA components just because `nrna_init == 0` from the
  strand-estimate failing at low SS. These components have intronic span
  and can receive intronic reads.

Implementation: remove the `nrna_init == 0.0` gating block. The single-exon
gating (which already runs earlier) handles the structurally impossible cases.

#### Corresponding change in `build_locus_em_data()` (locus.py ~L411–413)

```python
# Zero nRNA prior when nrna_init is zero for that nRNA.
nrna_init_local = estimator.nrna_init[unique_global_nrna]
prior[n_t:n_t + n_nrna][nrna_init_local == 0.0] = 0.0
```

**Change**: Remove these 3 lines. The single-exon gating above already
handles truly impossible nRNAs.

---

### 3. `src/rigel/locus.py` — Simplify nRNA Init

#### `compute_nrna_init()` (~L455–511)

**Current**: Computes strand-corrected nRNA init from intronic evidence.
Returns zeros when `2·SS − 1 ≤ 0.2`.

**Change**: Remove the `STRAND_DENOM_MIN` hard threshold. At low SS, use
total intronic coverage (sense + antisense) as the init signal rather than
the strand-corrected estimate (which is undefined at SS ≈ 0.5).

```python
def compute_nrna_init(
    intronic_sense, intronic_antisense,
    nrna_spans, nrna_max_exonic,
    strand_models,
):
    strand_spec = strand_models.strand_specificity
    denom = 2.0 * strand_spec - 1.0

    if denom <= STRAND_DENOM_MIN:
        # Low SS: use total intronic coverage as init
        # (can't separate nRNA from gDNA by strand, but intronic reads
        # are still evidence of nRNA or gDNA — let the EM decide)
        nrna_init = np.maximum(0.0, intronic_sense + intronic_antisense)
    else:
        raw = (intronic_sense - intronic_antisense) / denom
        nrna_init = np.maximum(0.0, raw)

    # Zero nRNAs with no intronic span (single-exon gene nRNAs)
    intronic_span = np.maximum(nrna_spans - nrna_max_exonic, 0.0)
    nrna_init[intronic_span <= 0] = 0.0

    return nrna_init
```

**Note**: `nrna_init` is still used to provide the EM with a rough scale
for nRNA components. It's passed to `batch_locus_em()` and used for
gating in `extract_locus_sub_problem()`. With the gating fix above, a
nonzero `nrna_init` keeps the component alive without constraining its
final value.

#### `build_locus_em_data()` (~L395–446)

**Changes**:
1. Remove the `nrna_init == 0` prior zeroing (3 lines at ~L411–413) — as
   noted in Section 2
2. Remove the `nrna_frac_alpha/beta` fields from the `LocusEMInput`
   construction (~L445–446)

**Keep**: The single-exon nRNA zeroing remains (L395–410) — this is
structurally motivated.

---

### 4. `src/rigel/scored_fragments.py` — Clean LocusEMInput

#### `LocusEMInput` dataclass (~L130–170)

**Change**: Remove these fields:
- `nrna_frac_alpha: np.ndarray`  (~L169)
- `nrna_frac_beta: np.ndarray`  (~L170)

**Keep**: All other nRNA-related fields (`local_to_global_nrna`,
`local_t_to_local_nrna`, `nrna_to_t_offsets`, `nrna_to_t_indices`,
`nrna_init`) — these are needed for posterior assignment and output mapping.

---

### 5. `src/rigel/priors.py` — Remove nRNA Fraction Priors

#### Functions to delete entirely:
- `_compute_hybrid_nrna_frac_vec()` (~L52–120)
- `_aggregate_nrna_frac_by_group()` (~L199–247)
- `compute_nrna_frac_priors()` (~L253–449)

#### Functions to keep:
- `compute_global_gdna_density()` (~L123–147) — used for gDNA priors
- `estimate_kappa()` (~L153–198) — used by `compute_eb_gdna_priors()`

#### Constants to remove:
- `STRAND_DENOM_EPS = 0.01` (~L14) — only used by deleted functions
- `DEFAULT_MEAN_FRAG = 200.0` (~L15) — only used by deleted functions

---

### 6. `src/rigel/estimator.py` — Remove nRNA Fraction State

#### `AbundanceEstimator.__init__()`

**Remove**:
- `self.nrna_frac_alpha` (~L149–150)
- `self.nrna_frac_beta` (~L151–152)
- `self.transcript_exonic_sense` (~L196–197)
- `self.transcript_exonic_antisense` (~L198–199)

**Keep**: `self.nrna_init` (~L146) — still used for nRNA initialization
signal.

#### `run_batch_locus_em()` or similar EM entry point

**Change**: Remove `self.nrna_frac_alpha` and `self.nrna_frac_beta` from
the argument list passed to the C++ `batch_locus_em()` call (~L387–388).

#### Re-exports

**Remove**: `compute_nrna_frac_priors` from imports/exports (~L45).

---

### 7. `src/rigel/config.py` — Remove 9 EMConfig Fields

**Remove**:
```python
nrna_frac_kappa_global: float | None = None      # L68
nrna_frac_kappa_locus: float | None = None        # L73
nrna_frac_kappa_nrna: float = 5.0                 # L77
nrna_frac_mom_min_evidence_global: float = 50.0   # L84
nrna_frac_mom_min_evidence_locus: float = 20.0    # L87
nrna_frac_kappa_min: float = 2.0                  # L90
nrna_frac_kappa_max: float = 200.0                # L92
nrna_frac_kappa_fallback: float = 5.0             # L94
nrna_frac_kappa_min_obs: int = 20                 # L96
```

**Handle gDNA shared params**: The `compute_eb_gdna_priors()` call in
pipeline.py currently borrows `kappa_min`, `kappa_max`, `kappa_fallback`, and
`kappa_min_obs` from the nRNA config fields. These 4 must be duplicated to
dedicated gDNA fields:

```python
# Add to EMConfig (these already have gDNA-specific counterparts for the
# other knobs — this completes the set):
gdna_kappa_min: float = 2.0
gdna_kappa_max: float = 200.0
gdna_kappa_fallback: float = 5.0
gdna_kappa_min_obs: int = 20
```

**Net change**: Remove 9 nRNA fields, add 4 gDNA fields → net -5 fields
(but the 4 gDNA fields replace shared params, so they're not new concepts).

---

### 8. `src/rigel/pipeline.py` — Remove nRNA Prior Computation

#### `_compute_priors()` or equivalent function (~L455–530)

**Remove**: The entire `compute_nrna_frac_priors()` call block (~L504–520).

**Remove**: The locus_id transcript assignment (~L479) — the comment says
"needed by nrna_frac prior cascade". Verify it's not needed by anything else
before removing. (The locus_id array may still be needed for output grouping —
check before removing.)

**Update**: The `compute_eb_gdna_priors()` call to use the new dedicated gDNA
config fields instead of `nrna_frac_kappa_min/max/fallback/min_obs`:

```python
gdna_inits = compute_eb_gdna_priors(
    loci, em_data, estimator, index, strand_models,
    intergenic_density=intergenic_density,
    kappa_ref=em_config.gdna_kappa_ref,
    kappa_locus=em_config.gdna_kappa_locus,
    mom_min_evidence_ref=em_config.gdna_mom_min_evidence_ref,
    mom_min_evidence_locus=em_config.gdna_mom_min_evidence_locus,
    kappa_min=em_config.gdna_kappa_min,         # was nrna_frac_kappa_min
    kappa_max=em_config.gdna_kappa_max,         # was nrna_frac_kappa_max
    kappa_fallback=em_config.gdna_kappa_fallback,  # was nrna_frac_kappa_fallback
    kappa_min_obs=em_config.gdna_kappa_min_obs,    # was nrna_frac_kappa_min_obs
)
```

**Remove**: Import of `compute_nrna_frac_priors` (~L43).

---

### 9. `src/rigel/scan.py` — Remove Exonic Accumulator Args

#### `fused_score_buffer()` call (~L131–140)

**Remove**: `estimator.transcript_exonic_sense` and
`estimator.transcript_exonic_antisense` from the argument list.

**Note**: This requires a corresponding change to the C++ `fused_score_buffer`
function to remove these output parameters. The C++ function currently
accumulates exonic sense/antisense counts per-transcript during scoring.
These accumulators need to be removed from the C++ signature and internal
logic.

---

### 10. `src/rigel/cli.py` — Remove CLI Arguments

#### `_ParamSpec` entries (~L415–423)

**Remove** all nRNA-fraction-related `_ParamSpec` entries (9 entries).

#### Argparse definitions (~L782–824)

**Remove**: `--nrna-frac-kappa-global`, `--nrna-frac-kappa-locus`,
`--nrna-frac-kappa-nrna`, `--nrna-frac-kappa-min`,
`--nrna-frac-kappa-max`, `--nrna-frac-kappa-fallback`,
`--nrna-frac-kappa-min-obs`, and any `--nrna-frac-mom-*` args.

#### Add gDNA MoM args (if not already present)

**Add**: `--gdna-kappa-min`, `--gdna-kappa-max`, `--gdna-kappa-fallback`,
`--gdna-kappa-min-obs` (4 CLI args) mapping to the new EMConfig fields.

---

### 11. `src/rigel/native/em_solver.cpp` — Update C++ `fused_score_buffer`

The `fused_score_buffer` function in C++ currently takes
`transcript_exonic_sense` and `transcript_exonic_antisense` as output arrays
and accumulates per-transcript exonic strand counts during the scoring pass.
These are only consumed by `compute_nrna_frac_priors()`.

**Change**: Remove the two output array parameters and the accumulation logic
from `fused_score_buffer()`. Update the nanobind binding accordingly.

---

### 12. `scripts/profiler.py` — Remove nRNA Prior Profiling

**Remove**: Import of `compute_nrna_frac_priors` (~L81) and its call (~L626).

---

## Test Changes

### Tests to remove:

- `tests/test_estimator.py`: `TestNrnaFracPriors` class (~L955–1060) — tests
  the 3-tier EB shrinkage that is being removed.

### Tests to update:

- `tests/test_em_impl.py`:
  - `_make_linked_locus()` helper (~L565–673): Remove `nrna_frac_alpha`,
    `nrna_frac_beta` parameters. For tests that exercised the hierarchical
    M-step behavior, convert them to plain EM tests or replace them with
    tests that verify nRNA competes naturally with mRNA.
  - `TestLinkedEmBasic` class: Tests like `test_empty_locus_nrna_frac_near_zero`
    and `test_theta_is_simplex` need updating — remove nrna_frac assertions,
    verify the plain EM still produces valid output.

- `tests/test_cli.py`:
  - Remove all `nrna_frac_kappa_*` test cases (~L76–188). Any tests that
    check YAML config loading with nrna_frac keys should verify those keys
    are now ignored (or warning-logged as unknown).
  - Add tests for the new `gdna_kappa_min/max/fallback/min_obs` CLI args.

- `tests/test_gdna.py`:
  - `compute_nrna_init` tests (~L474–495): Update to verify the new behavior
    at low SS (total intronic coverage fallback instead of zeros).

### Tests to add:

- **test_em_impl.py**: Add a test verifying that nRNA components compete
  independently with mRNA through the plain MAP-EM. Create a scenario with
  intronic-only reads (should go to nRNA) and exonic reads (should split
  between mRNA and nRNA based on likelihood).

- **test_gdna.py** or **test_em_impl.py**: Add tests verifying nRNA
  components are not killed at SS=0.5 — they should remain active with
  nonzero theta when intronic reads are present.

### Full test suite:

Run `pytest tests/` after all changes to catch integration regressions.

---

## Validation

1. **Re-run the 1,024-run complex locus ablation**:
   ```bash
   conda run -n rigel python scripts/synthetic_sim_sweep.py \
     -c scripts/complex_locus_config.yaml \
     -o /Users/mkiyer/Downloads/rigel_runs/complex_locus_phase1 \
     -v 2>&1 | tee phase1_sweep.log
   ```

2. **Run analysis**:
   ```bash
   conda run -n rigel python scripts/analyze_complex_locus.py \
     /Users/mkiyer/Downloads/rigel_runs/complex_locus_phase1/sweep_results.tsv \
     > /Users/mkiyer/Downloads/rigel_runs/complex_locus_phase1/analysis.txt
   ```

3. **Compare key metrics against baseline**:

   | Metric | Baseline | Phase 1 Target |
   |--------|----------|----------------|
   | Moderate nRNA accuracy (SS≥0.75, gdna=0) | -14% to -18% | < ±5% |
   | nRNA accuracy at SS=0.5 | -83% to -90% | Meaningful recovery |
   | Pattern 15 mRNA FP (SS=0.9, gdna=0.5) | +864% | Substantially reduced |
   | gDNA overestimate with moderate nRNA | +60-85% | Substantially reduced |
   | SS=0.5 gDNA overestimate | +136-815% | Substantially reduced |
   | Pure mRNA accuracy (no nRNA, no gDNA) | 0.0% | Unchanged |

4. **Run full test suite**: `pytest tests/ -q`

---

## Implementation Order

Recommended sequence to minimize broken intermediate states:

1. **Config** (config.py): Add 4 gDNA fields, remove 9 nRNA fields
2. **CLI** (cli.py): Update args to match new config
3. **Pipeline** (pipeline.py): Remove `compute_nrna_frac_priors()` call,
   update `compute_eb_gdna_priors()` call to use new config fields
4. **Priors** (priors.py): Delete 3 functions + 2 constants
5. **Locus** (locus.py): Fix `compute_nrna_init()`, remove gating in
   `build_locus_em_data()`, remove α/β from `LocusEMInput` construction
6. **Scored fragments** (scored_fragments.py): Remove α/β from dataclass
7. **Estimator** (estimator.py): Remove α/β state, exonic accumulators,
   re-export
8. **C++ EM solver** (em_solver.cpp): Remove hierarchical M-step branch,
   clean signatures, fix nRNA gating, update bindings
9. **C++ scoring** (em_solver.cpp or scoring.cpp): Remove exonic accumulator
   args from `fused_score_buffer()`
10. **Scan** (scan.py): Remove exonic accumulator args
11. **Tests**: Update/remove/add as described
12. **Scripts** (profiler.py): Remove nRNA prior profiling

---

## Risk Assessment

**Low risk**:
- Removing the EB shrinkage and hierarchical M-step is predominantly deletion
- The plain MAP-EM path (`n_nrna == 0` branch) is already well-tested and
  handles the strand symmetry penalty correctly
- nRNA components continue to participate in the EM with proper likelihoods

**Medium risk**:
- The nRNA gating fix at low SS could allow spurious nRNA components to
  absorb gDNA mass at SS=0.5 (intronic gDNA reads could be attributed to
  nRNA). The strand symmetry penalty should mitigate this, but we need to
  verify with the ablation.
- Tests that relied on the hierarchical M-step behavior need careful
  conversion — watch for accidental regressions in test logic.

**Mitigation**:
- The 1,024-run ablation sweep provides comprehensive coverage of the
  accuracy impact across SS, gDNA, and nRNA levels
- If SS=0.5 gDNA accuracy degrades after the gating fix, Phase 3 (scale
  strand penalty by SS) should resolve it
