# v4.1 Implementation Dry-Run

**Status**: dry-run document, no code changed.
**Date**: 2026-05-19
**Companion to**: `gdna_exposure_plan_v4.1.md`.

Purpose: walk every v4.1 step against the current codebase, list the
exact files and call-sites that change, and surface problems that the
abstract plan glossed over. Problems are listed with a recommendation;
implementation should not start until they are settled.

---

## Code map (what already exists)

| Concept | Location | Notes |
|---|---|---|
| Per-region exposure build | `src/rigel/calibration/_regional_exposure.py` `RegionalGdnaExposure.build` | Currently per-class `ρ_ref`, per-class `signal`. Must change. |
| Exposure lookup on a span | `RegionalGdnaExposure.weighted_length_on_ref` | Already step-function-aware; piecewise integrates A over a per-ref `[start, end)`. Reused as-is. |
| Per-position log A lookup | `RegionalGdnaExposure.log_weights_for_positions` | Used only by `_apply_unit_gdna_weights` (numerator path). Becomes dead code in v4.1 production; safe to keep for diagnostics. |
| gDNA overlap effective length (unweighted) | `src/rigel/calibration/_exposure.py` `gdna_eff_len_for_loci` | Truth source for §9.1.2 parity. Has a single-locus fast path. |
| gDNA overlap effective length (weighted) | `_exposure.py` `weighted_gdna_eff_len_for_loci` | Already uses the midpoint convention; already merges per-ref expanded windows. v4.1's `weighted_eff_len_overlap` is a rename + light refactor. |
| Per-transcript unweighted effective length | `src/rigel/frag_length_model.py` `FragmentLengthModel.compute_all_transcript_eff_lens` | Vectorized closed-form `(L+1)(CDF[k]−p0) − CMOM[k]`. Truth source for §9.1.3 parity. |
| Per-transcript exon block table | `TranscriptIndex._t_exon_intervals` (`dict[int, np.ndarray]`) | Already cached at index load. Includes synthetic nRNA single-block entries. Footprint input for containment. |
| Transcript ref id | `TranscriptIndex.t_to_ref_arr` | int32, populated; reused. |
| Transcript strand | `index.t_df["strand"]` (cat code) | Need helper to materialize as `int8 ±1`. |
| Synthetic mask | `index.t_df["is_synthetic"]` | Already used by estimator. |
| Effective-length plumbing into EM | `TranscriptGeometry.effective_lengths` → `AbundanceEstimator._t_eff_len` (single array) | **Critical**: same array is read by native EM (estimator line 311) *and* by TPM / reports (lines 465, 555, 598, 669, 690). v4.1 split is mandatory and not local. |
| Per-locus gDNA effective length | `PriorTable.gdna_eff_len` (already exists as `L_g^A` under regional mode) and `PriorTable.gdna_eff_len_unweighted` (`L_g`) | Both already populated by `assemble_priors`. Naming is misleading: in uniform mode they are equal, in regional mode `gdna_eff_len = L_g^A`. v4.1 should rename for clarity. |
| Per-locus gDNA prior | `PriorTable.gdna_prior_count` (`η_g`, unrescaled) and `PriorTable.gdna_prior_count_regional` (a *different* diagnostic, not the rescaled prior) | The current `gdna_prior_count_regional` is **not** the v4.1 rescaling; it's a regional Poisson-rate diagnostic via `_expected_gdna_count_regional_diag`. v4.1 must add a new column (the actual EM-consumed prior) and clarify or remove the existing diagnostic. |
| Numerator application | `pipeline.py` `_apply_unit_gdna_weights` at line 85, called at line 1073 | Single call site. Easy delete. |
| `RegionalWeightApplicationStats` | `_regional_exposure.py` lines 60–82, summarized via `calibration.with_regional_weighting_stats` | Need a sunset plan. |
| Pipeline orchestrator | `pipeline.py` `quant_from_buffer` and `_setup_geometry_and_estimator` (line 531) | `effective_lengths` is computed once here from `rna_fl.compute_all_transcript_eff_lens(exonic_lengths)` and fed to `TranscriptGeometry`. v4.1 must compute *two* arrays here. |

---

## Problems found during dry-run

### P1 — `AbundanceEstimator` has one effective-length array shared by EM and TPM

`AbundanceEstimator._t_eff_len` is the only effective-length state on
the estimator. It is consumed by:

- native EM, [estimator.py L311](src/rigel/estimator.py#L311);
- `quant.feather` TPM / `effective_length` column, [L465](src/rigel/estimator.py#L465), [L492](src/rigel/estimator.py#L492);
- gene-aggregated `effective_length`, [L555](src/rigel/estimator.py#L555), [L598](src/rigel/estimator.py#L598);
- nRNA `effective_length` / TPM, [L669](src/rigel/estimator.py#L669), [L690](src/rigel/estimator.py#L690), [L696](src/rigel/estimator.py#L696).

There is also a public `effective_lengths` property at
[L221](src/rigel/estimator.py#L221) that returns the EM array.

v4.1 §7.3 says split into `effective_lengths_em` and
`effective_lengths_output`. This is correct but not trivial: every
read of `_t_eff_len` above must be audited. Tests likely pin behavior
of the public `effective_lengths` property.

**Recommendation**:

- Add two private arrays: `_t_eff_len_em` and `_t_eff_len_output`.
- Keep the public `effective_lengths` property pointing at
  `_t_eff_len_output` (so external consumers continue to get the
  user-visible length). Add `effective_lengths_em` property for the EM
  path.
- Inside the EM-batch call, read `_t_eff_len_em`. Everywhere else,
  read `_t_eff_len_output`.
- In uniform / `off` mode, both are the same object.
- Extend `TranscriptGeometry` with an optional second field
  `effective_lengths_em`. Default `None` → use `effective_lengths`
  for EM too (uniform-mode parity).

### P2 — `gdna_eff_len` and `gdna_prior_count` names already mean weighted things

`PriorTable.gdna_eff_len` already stores `L_g^A` under regional mode
(populated by `weighted_gdna_eff_len_for_loci`). It is the *EM-consumed*
length today. `PriorTable.gdna_eff_len_unweighted` stores `L_g`. There
is no separate "output" gDNA effective length used anywhere right now
because reports use `quant`/`nrna_quant`/`gene_quant` from transcript
state, not the per-locus gDNA length. The `loci.feather` `gdna_eff_len`
column is the only public surface and it currently shows `L_g^A`.

For `gdna_prior_count`: today this is `η_g` (unrescaled). EM consumes
it directly. `gdna_prior_count_regional` is a *diagnostic* unrelated to
v4.1's prior rescaling (it is a Poisson-rate estimate, not
`η_g · L_g^A / L_g`).

**Recommendation**:

- Reuse `PriorTable.gdna_eff_len` to mean "EM-consumed length"
  (continues to be `L_g^A` in regional mode). No rename — current
  behavior already matches v4.1 intent.
- Add `PriorTable.gdna_prior_count_em` (new) for the rescaled prior
  `η_g · L_g^A / L_g`; have native EM consume this. Keep
  `gdna_prior_count = η_g` (unrescaled) as the published "what
  calibration thought" diagnostic. Migrate `loci.feather` to expose
  both columns.
- *Remove* `gdna_prior_count_regional` and the diagnostic emitter
  `_expected_gdna_count_regional_diag` — it confuses readers and is
  not what v4.1 means by "rescaled prior".
- For `loci.feather`: publish `gdna_eff_len` (EM-consumed weighted),
  `gdna_eff_len_unweighted` (raw `L_g`), `gdna_prior_count` (raw
  `η_g`), `gdna_prior_count_em` (the value EM actually saw).

This deviates from v4.1 §7.2 (which proposed
`gdna_em_exposure_weight`, `gdna_em_effective_length`,
`gdna_prior_count_unweighted`). The substance is the same; the column
names should follow the existing convention to avoid breaking
downstream parsers.

### P3 — "One global signal" auto-uniform loses a real degradation guard

Current code computes a `signal_per_class` and uses it to attenuate
`log_weight` per class via `log_weight[mask] = max(signal * raw, FLOOR)`.
This protects against noisy classes (e.g. EXON-INTRON boundary signal
is weak but INTRON signal is strong → attenuate the noisy channel
toward identity). v4.1 §2.1 §2.2 says "one global signal" and "no
per-class normalization".

The arguments in v4.1 §2.2 are about the *reference quantile* `ρ_ref`,
which I agree must be global. But they do **not** automatically imply
that the *signal attenuator* should also be global. The signal is a
noise-floor diagnostic, not a calibration. Killing the per-class
attenuator means a class with `E_r` so small that `ρ̂_r` is essentially
Poisson noise still injects high-variance `A_r` into EM denominators.

**Recommendation**:

- Use **one global `ρ_ref`** (v4.1 §2.1 step 2): mandatory; this is
  the actual fix.
- Retain a **per-class signal attenuator** that multiplies the
  per-region `log A_r` toward zero when that class lacks evidence to
  trust its `ρ̂`. Document this explicitly: this is a noise control,
  not a calibration step; classes with no spread above the null get
  attenuated to identity in their own regions and contribute nothing
  to EM denominators.
- Alternative if user disagrees: gate the entire field with one
  global signal (`signal` = max of class signals, or a pooled
  Cochran-style test). Simpler, but loses the ability to mask one
  noisy class while trusting others. I prefer per-class attenuation
  with global `ρ_ref`.

This is a planning question to settle before §10 step 1 lands.

### P4 — Performance: `weighted_length_on_ref` Python loop is the bottleneck

`weighted_length_on_ref` does a per-ref linear scan from the first
hit region (`while cursor < end_int and i < ref_starts.size`). For
gDNA loci this is fine (`O(k loci × k_regions_per_locus)`). For
**every transcript × every fragment length × every exon block** under
containment, total calls are roughly:

- 200k transcripts × ~40 effective FL bins × ~10 exons ≈ **80M calls**.

Each call is at minimum a `searchsorted` + a small Python loop. At ~1
µs per call this is ~80 s of overhead per quant run — already over
the §9.3 G9 +10% gate on a small sample.

**Recommendation**:

- Add a one-shot per-ref prefix-sum array of A over region boundaries
  at exposure-build time:
  ```
  exposure._csum_by_ref: dict[int, np.ndarray]   # shape (R_ref+1,)
  exposure._cstart_by_ref: dict[int, np.ndarray] # region starts incl. sentinel
  ```
  with `_csum_by_ref[ref][i] = Σ_{j<i} A_j · (end_j − start_j)`.
  Then `weighted_length_on_ref(ref, lo, hi)` is two binary searches
  and a few arithmetic ops — `O(log R_ref)` constant-time.
- Make `weighted_length_on_ref` use this fast path when arrays are
  built; fall back to the current implementation otherwise (covers
  the uniform / empty-ref cases without changes elsewhere).
- Vectorize containment: instead of one `weighted_length_on_ref` call
  per `(ℓ, block)`, build numpy arrays `g_lo[ℓ, block]` and
  `g_hi[ℓ, block]` per transcript and call a vectorized variant
  `weighted_length_on_ref_batch(ref, g_lo_arr, g_hi_arr)` that does
  one searchsorted on the whole batch.
- Truncate `ℓ` to `pmf ≥ 1e-6` (already conventional).
- Commit to the G9 gate. If the vectorized Python path still misses,
  port `weighted_eff_len_contained` to C++; native ABI unchanged.

### P5 — `_apply_unit_gdna_weights` deletion and stats sunset

The function is only called from one place
([pipeline.py L1073](src/rigel/pipeline.py#L1073)). Removal is
mechanical. `RegionalWeightApplicationStats` is summarized via
`calibration.with_regional_weighting_stats` and surfaces in
`summary.json["regional_exposure"]["application"]`.

**Recommendation**:

- Delete the call site and the function in one PR.
- Stop emitting `summary.json["regional_exposure"]["application"]`.
- Keep the `RegionalWeightApplicationStats` class for one release
  with a deprecation comment, in case external consumers parse it;
  it never gets populated.
- Drop the `genomic_midpoint` propagation later (separate PR) if
  no other consumer remains; today it's only used by
  `_apply_unit_gdna_weights`. Keep it as a fragment-level diagnostic
  exposed in summary, but stop computing midpoints in the hot path
  if profiling shows benefit.

### P6 — Strand projection for nRNA spans is a no-op

For N=1 footprints, the containment §3.2.1 projection on a single
block reduces to:

- strand +1: `g_lo = g_start + (t_lo − 0)`, `g_hi = g_start + (t_hi − 0)`
- strand −1: `g_lo = g_end − t_hi`, `g_hi = g_end − t_lo`

Since `(g_end − g_start)` equals the block length and the integral of
A over `[g_lo, g_hi)` is symmetric in (lo, hi), both strand choices
yield the same integrand boundaries (modulo orientation). So strand
**doesn't matter for single-block containment**. Document this so the
nRNA-span caller doesn't waste cycles materializing strand. mRNA
multi-exon strand still matters because exon order in transcript
coordinates differs.

### P7 — `enable_gdna[i] == 0` edge case for prior rescaling

`L_g^A / L_g` requires `L_g > 0`. `assemble_priors` already floors
`gdna_eff_len_for_loci` at `min_value=1.0`, but the v4.1 prior
rescaling expression should still guard:

```python
ratio = L_g_A / max(L_g, MIN_EFF_LEN)
eta_g_em = eta_g * ratio
```

For loci with `enable_gdna == 0`, EM disables the gDNA component
entirely; the rescaled prior is unused regardless. No behavior change.

### P8 — `weighted_gdna_eff_len_for_loci` ↔ `weighted_eff_len_overlap` refactor must preserve fast paths

`gdna_eff_len_for_loci` has a single-interval fast path
([_exposure.py L95–105](src/rigel/calibration/_exposure.py#L95)) that
avoids per-`ℓ` iteration. The weighted version does not currently
have an equivalent fast path. The refactor should keep these paths
intact and verify §9.1.2 against representative shapes.

### P9 — Constant `A_0` test depends on full-coverage assumption

§9.1.4's "`L_k^A = A_0 · L_k`" only holds when every byte of the
component's footprint lies inside a region in `RegionArrays`. If a
transcript exon spills past the last annotated region on a contig
(the "spill" path in `weighted_length_on_ref`), the tail is integrated
as identity (`A=1`), not as `A_0`. This is a property of the
"contig has regions but query overruns" branch, not a bug.

**Recommendation**:

- Construct the §9.1.4 fixture so all component footprints sit
  strictly inside the regional partition.
- Add a separate test that verifies the spill-past-last-region
  behavior is documented (returns identity).

---

## Step-by-step implementation sequence

This is what I would do once P3 (per-class signal) is settled and any
naming concerns in P2 are answered.

### Step 1 — Single global `ρ_ref` and per-class signal attenuation
*(blocked on P3 decision)*

File: `src/rigel/calibration/_regional_exposure.py`

- Replace `rho_ref = max(rho_ref_per_class.values())` with one global
  weighted-Q95 over all regions using `E_r` as weights.
- Decide per-class signal attenuator (P3). Default to **keep** the
  per-class attenuator multiplying `log A_r` toward zero.
- `per_class` summary keeps a per-class diagnostic
  `rho_ref_class` for backward compatibility of `summary.json` but
  none of those values are used to compute `A_r`.
- Add `rho_ref_global` (new) field to the dataclass and to the
  summary dict.

New tests:
- `tests/calibration/test_regional_exposure_global_ref.py`:
  the §2.3 sentinel + cross-class equivalence.

### Step 2 — Integration engine and overlap rule

File: `src/rigel/calibration/_exposure.py`

- Introduce `integrate_exposure_over_midpoints(windows, exposure)` as
  a thin wrapper around `weighted_length_on_ref`.
- Refactor `weighted_gdna_eff_len_for_loci` so its per-`ℓ` merge step
  emits midpoint-window tuples and calls the engine. External
  signature unchanged.
- Add public alias `weighted_eff_len_overlap` that takes a
  `footprint` list directly (decoupled from `Locus`/`MultiLocus`
  types); `weighted_gdna_eff_len_for_loci` becomes a thin wrapper
  that translates `MultiLocus.loci` into footprint tuples and
  forwards.

New tests:
- `tests/test_exposure_overlap_engine.py`: §9.1.2 parity vs
  `gdna_eff_len_for_loci` under uniform `A`; §9.1.8 merge stress.

### Step 3 — Containment rule (Python, with prefix-sum exposure)
*(performance gate per P4)*

File: `src/rigel/calibration/_exposure.py`

- Add `weighted_eff_len_contained(footprint, fl, exposure, *, strand,
  min_value)`. Single-pass loop over positive `ℓ`, per-block midpoint
  range as a numpy slice, batched `weighted_length_on_ref_batch`
  call.
- Add `weighted_length_on_ref_batch(ref, g_lo_arr, g_hi_arr)` on
  `RegionalGdnaExposure`. Uses per-ref prefix sums.
- Add lazy build of per-ref prefix sums in
  `RegionalGdnaExposure._build_prefix_sums()`; build once per quant.

New tests:
- `tests/test_exposure_contained_engine.py`: §9.1.3 parity vs
  `compute_all_transcript_eff_lens` under uniform A; §9.1.6
  negative-strand exon-order test; §9.1.4 constant-`A_0` length
  law with full-coverage fixture.

### Step 4 — Per-transcript and per-nRNA-span weighted lengths

File: new `src/rigel/calibration/_exposure_lengths.py` (or extend
`_exposure.py`).

- `weighted_transcript_eff_lens(index, rna_fl, exposure)` →
  `np.ndarray[n_t]`. Iterates `index._t_exon_intervals` and per-row
  `t_to_ref_arr` / strand; calls `weighted_eff_len_contained` per
  transcript with the batched fast path.
- For synthetic rows, `weighted_eff_len_contained` with the N=1
  footprint (strand irrelevant per P6). No new `nrna_span_df`
  needed — synthetic rows already have single-block exon entries.
- Returns `effective_lengths_em` array shape `(n_t,)` in raw `t_idx`
  order.

This *drops* the v4.1 §10 step 5 ("nRNA span bookkeeping") because
the existing `_t_exon_intervals` already keys synthetic spans by
`t_idx`; no separate `nrna_span_df` is needed in the index.

New tests:
- `tests/test_exposure_transcript_lengths.py`: uniform parity vs
  `rna_fl.compute_all_transcript_eff_lens`; constant-`A_0` law.

### Step 5 — Estimator EM/output array split

Files: `src/rigel/config.py`, `src/rigel/estimator.py`,
`src/rigel/pipeline.py::_setup_geometry_and_estimator`.

- `TranscriptGeometry.effective_lengths_em: np.ndarray | None = None`.
- `AbundanceEstimator.__init__` stores `_t_eff_len_em` (defaults to
  `_t_eff_len_output` when geometry omits it). Public
  `effective_lengths` property continues to return
  `_t_eff_len_output` (no breaking change for external consumers).
  New `effective_lengths_em` property returns `_t_eff_len_em`.
- `run_batch_locus_em_partitioned` reads `_t_eff_len_em` for the EM
  arg. All TPM / report code already reads via the
  `_t_eff_len`/`effective_lengths` accessor — repoint to
  `_t_eff_len_output`.
- `_setup_geometry_and_estimator` keeps building unweighted
  `effective_lengths` and now also builds `effective_lengths_em` via
  Step 4 when `regional_exposure.mode == "regional"`.

New tests:
- `tests/test_estimator_em_output_split.py`: uniform parity is
  bit-exact; regional differs only in EM denominator; public
  `effective_length` column stays unweighted.

### Step 6 — gDNA prior rescaling

File: `src/rigel/calibration/locus_prior.py::assemble_priors`.

- After `gdna_eff_len_arr[idx] = weighted_gdna_eff_len_for_loci(...)`
  is computed, add:
  ```python
  ratio = gdna_eff_len_arr[idx] / max(gdna_eff_len_unweighted_arr[idx], MIN_EFF_LEN)
  gdna_prior_count_em_arr[idx] = eta_g * ratio
  ```
- Add `PriorTable.gdna_prior_count_em` (np.ndarray, float64).
- Drop `gdna_prior_count_regional` (P2). Drop the diagnostic
  emitter `_expected_gdna_count_regional_diag`.
- Wire `gdna_prior_count_em` through `_run_locus_em_partitioned` to
  native EM (replacing `gdna_prior_count` arg).

New tests:
- §9.1.5 prior density invariant; ratio sanity check.

### Step 7 — Delete numerator path

Files: `src/rigel/pipeline.py`.

- Remove the `if calibration.regional_exposure is not None:` block at
  L1071–1078 that calls `_apply_unit_gdna_weights`.
- Remove `_apply_unit_gdna_weights` function at L85.
- Remove `calibration.with_regional_weighting_stats(...)` call site.
- Stop emitting `summary.json["regional_exposure"]["application"]`.
- Keep `RegionalWeightApplicationStats` dataclass for one release
  with deprecation comment, returning zeros if any caller persists.
- Keep `ScoredFragments.genomic_midpoint` for now (fragment-level
  QC). Mark for cleanup once nothing reads it.

No new tests; existing pipeline tests must keep passing.

### Step 8 — Diagnostic columns

Files: estimator-emitting code paths for `quant.feather`,
`nrna_quant.feather`, `loci.feather`.

- `quant.feather`: add `em_exposure_weight` (= `L_t^A / L_t`) and
  `em_effective_length` (= `L_t^A`).
- `nrna_quant.feather`: same two columns.
- `loci.feather`: rename per P2 — emit `gdna_eff_len`,
  `gdna_eff_len_unweighted`, `gdna_prior_count`, `gdna_prior_count_em`,
  and `gdna_em_exposure_weight` (= ratio).

New tests:
- Column-presence and dtype guards.

### Step 9 — Synthetic two-region locus

File: `tests/test_regional_exposure_locus.py` (new).

§9.2 implementation: a hand-built `RegionArrays` + minimal
`MultiLocus` and `ScoredFragments`. Validate that v4.1 routing
recovers truth within tolerance; uniform-override is bit-exact to
pre-v4.

### Step 10 — VCaP gates and 24-condition regression

Run existing harness. No code change. Document numbers in a
follow-up PR description.

---

## Estimated surface area

- New files: 3
  (`_exposure_lengths.py`, plus 3 test modules).
- Files touched: ~6
  (`_regional_exposure.py`, `_exposure.py`, `locus_prior.py`,
   `pipeline.py`, `estimator.py`, `config.py`).
- Lines deleted: ~150 (numerator path + per-class signal attenuator
  if P3 picks the "single global signal" route).
- Lines added: ~400.
- Native code touched: 0.
- New CLI flags: 0.

---

## Questions for the user before implementation begins

1. **P2 naming.** Are you OK keeping `PriorTable.gdna_eff_len` as
   "EM-consumed (weighted)" and `gdna_eff_len_unweighted` as raw,
   and adding `gdna_prior_count_em` alongside the unrescaled
   `gdna_prior_count`? This deviates from v4.1 §7.2's column names
   (`gdna_em_effective_length`, `gdna_em_exposure_weight`,
   `gdna_prior_count_unweighted`) but matches the existing
   convention.
2. **P3 signal attenuation.** Single global signal (simpler) or
   per-class signal attenuator with single global `ρ_ref` (safer
   against a noisy class poisoning the field)? My recommendation
   is the latter.
3. **§10 step 5 nRNA span bookkeeping.** I can drop it since
   `_t_exon_intervals` already covers synthetic spans. Confirm.
4. **`genomic_midpoint`.** Keep in `ScoredFragments` for now or
   plan removal in a follow-up?
5. **Native containment port.** Authorize a fall-back to C++ if the
   vectorized Python path misses the G9 +10% gate, or hold the line
   on Python-only?

Once these are answered I can land Step 1 in the next turn.
