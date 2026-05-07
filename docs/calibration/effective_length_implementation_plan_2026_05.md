# Effective-length consistency: implementation plan

**Status:** plan ready for implementation
**Date:** 2026-05-07
**Companion to:** [effective_length_consistency_2026_05.md](effective_length_consistency_2026_05.md) (audit + theory)
**Decision input:** synthesizes user feedback that (a) overhang is intentional and must not be hard-clipped, (b) `log h(ℓ_f)` can be candidate-specific in rigel, (c) gDNA must be Phase 2 with its own derivation, (d) posterior assignment must use the same `log_eff_len` as SQUAREM.

---

## 0. The contract we are committing to

```text
RNA candidate likelihood (per-fragment, per-candidate):
    log_lik = log_strand
            + log h_R(f_len_candidate)        ← candidate-specific FL term
            + overhang_penalty                 ← keeps annotation tolerance
            + mismatch_penalty

RNA component normalization (per-iteration, inside EM):
    log_weight[t] = log(theta[t] + eps) − log(L̃_t)

TPM (post-EM):
    TPM_t = (count_t / L̃_t) · 1e6 / Σ_t' (count_{t'} / L̃_{t'})

Invariants:
    (I1) The denominator inside the EM weight is the same L̃_t used by TPM.
    (I2) The denominator inside posterior assignment is the same L̃_t used by EM.
    (I3) No raw per-fragment start-count denominator (max(L − ℓ_f + 1, 1))
         remains anywhere in the RNA path.
    (I4) Overhang is allowed; it is penalized by overhang_penalty, not by
         hard-clipping.
```

`L̃_t` is the FL-marginal effective length already implemented in
[`FragmentLengthModel.compute_all_transcript_eff_lens`](src/rigel/frag_length_model.py).

---

## 1. Phase plan (sequence + landing criteria)

| Phase | Title | Touches | Lands when |
|------:|-------|---------|------------|
| **0** | Contract + tests | docs + new failing tests | Tests merged red, contract pinned in this doc |
| **1** | mRNA partitioned EM uses `L̃_t` | `pipeline.py`, `estimator.py`, `em_solver.cpp` (partitioned path), `_em_impl` ABI | All new tests green; goldens regenerated; synthetic-sim `gdna_none` FP count drops |
| **2** | Posterior-assignment consistency | `em_solver.cpp::assign_posteriors` + ABI | Assignment uses same `log_eff_len` as SQUAREM; new consistency test passes |
| **3** | Legacy non-partitioned EM aligned | `em_solver.cpp` legacy path + tests | Legacy path uses `effective_lens` only; `apply_bias_correction_uniform` deleted |
| **4** | Diagnostics + overhang audit | `em_solver.cpp` per-locus stats, `pipeline.py` summary | New stats land; we have data to decide whether overhang is ever pathological |
| **5** | Cleanup | remove `tx_starts/tx_ends` from partitioned payload if unused; quarantine `bias.py` | No live callers remain; tests still green |
| **6** | (later) gDNA derivation | new design note → implementation | Geometry first-principles note approved before code changes |

Each phase compiles and tests pass on its own. Phases 1–3 together are the model-changing patch; 4–5 are hygiene; 6 is a separate project.

---

## 2. Phase 0 — Contract + tests (red first)

**Goal:** Lock the contract and write the failing tests that the implementation will satisfy.

**Files added/edited:**

* This doc (contract above).
* New tests in `tests/test_em_eff_length.py`:

```python
def test_em_responsibility_follows_Ltilde_ratio():
    """Two transcripts with identical per-candidate evidence and unequal
    L̃ should split EM mass in ratio L̃_long : L̃_short (long gets more)."""

def test_assignment_uses_same_log_eff_len_as_em():
    """assign_posteriors output must equal one extra EM E-step with the
    same log_eff_len."""

def test_overhang_candidate_is_kept_and_penalized():
    """A fragment with f_len > t_length must remain in the candidate set
    and receive overhang_penalty, not -inf."""

def test_no_per_fragment_length_correction_in_log_lik():
    """After scoring + bias correction, log_lik for two candidates that
    differ only in their candidate transcript's L should NOT include a
    log(L − ℓ_f + 1) term. Verified by calling the scorer in isolation
    and reading log_lik."""

def test_short_unexpressed_tx_no_phantom_count_in_gdna_none():
    """Synthetic-sim regression: top short FP transcript from
    gdna_none_ss_0.99 condition has count < 100 (was 2,749 pre-fix)."""
```

The last test is a true regression test — we do not need to run the full sim; mock a small locus that reproduces the pathology (one short tx + one long tx + ambiguous fragments).

* New test in `tests/test_calibration_geometry.py`:

```python
def test_exon_intron_straddler_does_not_enter_intron_body_density():
    """A fragment that overlaps both an exon and an intron region
    increments mask 0b011 (mixed), enters EXON-INTRON boundary flux,
    but is NOT counted toward INTRON body density (which uses the pure
    MASK_INTRON column)."""
```

This pins the user's correct observation about the calibration accumulator.

**Acceptance:** all tests merged in `failing` state; CI marker `xfail(reason="phase 1 fix pending")`.

---

## 3. Phase 1 — mRNA partitioned EM uses `L̃_t`

**Goal:** The production path (partitioned batch EM) normalizes by `L̃_t`. Per-fragment formula-A correction is removed from this path.

### 3.1 Python side

`src/rigel/pipeline.py` already computes `effective_lengths = rna_fl.compute_all_transcript_eff_lens(...)` on line 386 and stores it on the `TranscriptGeometry` → `AbundanceEstimator._t_eff_len`. **No change here.**

`src/rigel/estimator.py`: in the call site that invokes `_batch_locus_em_partitioned` (search for the partitioned-path entry), replace whatever transcript-length array we pass today with `self._t_eff_len`. Rename argument:

```python
# was: t_lengths = self._exonic_lengths
t_eff_lens = self._t_eff_len
_em_impl.batch_locus_em_partitioned(
    ...,
    t_eff_lens=t_eff_lens,   # renamed from t_lengths or profile_lengths
    ...
)
```

(Verify exact call site and parameter name — the audit found `profile_lengths` in `apply_bias_correction_uniform`; that argument disappears.)

### 3.2 C++ side — `em_solver.cpp::run_em_with_partition`

```cpp
// REMOVE this block (was the bias correction call):
// apply_bias_correction_uniform(
//     sub.log_liks.data(), sub.t_indices.data(),
//     sub.tx_starts.data(), sub.tx_ends.data(),
//     sub.bias_profiles.data(), sub.n_t, n_candidates);

// REPLACE the hard-coded zeros:
// std::vector<double> log_eff_len(nc, 0.0);

// WITH a per-component slice driven by the new ABI input t_eff_lens:
std::vector<double> log_eff_len(nc);
for (int i = 0; i < sub.n_t; ++i) {
    int32_t gt = sub.local_to_global_t[i];
    log_eff_len[i] = std::log(std::max(t_eff_lens[gt], 1.0));
}
log_eff_len[sub.gdna_idx] = 0.0;   // PHASE 1 SCOPE: gDNA stays as today
```

Pass `t_eff_lens` (new `f64_1d` argument) into `batch_locus_em_partitioned`. Update `nb::arg` registrations and Python call site.

### 3.3 ABI change

`_em_impl.batch_locus_em_partitioned` gains a new keyword argument `t_eff_lens: np.ndarray[float64]`. Order in the C++ signature: place it adjacent to the existing transcript-array arguments. Bump no version since rigel is a single-source build.

### 3.4 Build + verify

```bash
conda activate rigel && pip install --no-build-isolation -e .
pytest tests/ -v                    # expect Phase 0 tests now green
pytest tests/ --update-golden       # regenerate golden outputs
```

Run the synthetic sim:

```bash
rm -rf /Users/mkiyer/Downloads/rigel_runs/sim_synthetic/gdna_*/{rigel_out,annotated.bam}
python scripts/sim/run_rigel_analysis.py
```

**Acceptance:**

* All Phase 0 tests green.
* Goldens regenerated (intentional, document in the PR).
* Synthetic-sim `gdna_none_ss_0.99` n_FP drops from 44 to single digits.
* Spearman / Pearson(log) unchanged or improved on `gdna_none` conditions.

---

## 4. Phase 2 — Posterior-assignment consistency

**Goal:** Honor invariant **(I2)** — the final per-fragment assignment uses the same `log_eff_len` as the EM iterations.

### 4.1 The bug

[`assign_posteriors`](src/rigel/native/em_solver.cpp) (line ~1373 in the audit) explicitly comments:

```cpp
// Effective lengths are all 1.0, so log_eff_len = 0.
std::vector<double> log_weights(nc);
for (int c = 0; c < nc; ++c) {
    log_weights[c] = std::log(theta[c] + EM_LOG_EPSILON);
}
```

After Phase 1, `log_eff_len` is no longer 1.0/0.0 for transcripts. Without this fix, EM converges with `L̃`-aware responsibilities, but final assignment falls back to no-normalization — silently giving short transcripts a bonus at the assignment step.

### 4.2 The fix

```cpp
// Accept log_eff_len as input from caller.
// Replace the loop above with:
for (int c = 0; c < nc; ++c) {
    log_weights[c] = std::log(theta[c] + EM_LOG_EPSILON) - log_eff_len[c];
}
```

Pass `log_eff_len.data()` from `run_em_with_partition` into `assign_posteriors`. No new ABI surface — the data already exists in scope.

**Acceptance:**

* `test_assignment_uses_same_log_eff_len_as_em` passes.
* Synthetic-sim numbers do not change much vs. Phase 1 (assignment is largely deterministic for unambiguous fragments) but a controlled multimapper test should show measurable shift.

---

## 5. Phase 3 — Legacy non-partitioned EM aligned

**Goal:** The non-partitioned `_em_impl` path (used by `tests/test_em_impl.py` and a handful of legacy callers) follows the same contract.

### 5.1 Today

[`em_solver.cpp`](src/rigel/native/em_solver.cpp) lines 1140–1163: the non-partitioned path accepts `effective_lens` (already!) but **also** calls `apply_bias_correction_uniform`. So the legacy path applies length correction *twice* in inconsistent ways. Phase 1 leaves this path untouched.

### 5.2 The fix

* Remove the `apply_bias_correction_uniform` call from the legacy path.
* Verify `effective_lens` callers pass FL-marginal values (likely already true in tests).
* Update any test fixtures that pre-baked formula-A corrections into expected outputs.
* Delete `apply_bias_correction_uniform` (no remaining callers) or keep it as a clearly-marked test helper if convenient.

**Acceptance:**

* `pytest tests/test_em_impl.py` green.
* `apply_bias_correction_uniform` either deleted or has zero production callers (audit with `grep`).

---

## 6. Phase 4 — Diagnostics + overhang audit

**Goal:** Generate the data we need to decide whether overhang is ever pathological. The user is right that we should not hard-gate without evidence.

### 6.1 New per-locus stats

In `em_solver.cpp` E-step, when `emit_locus_stats=True`, accumulate:

```cpp
// Per-locus
int64_t n_cand_flen_gt_tlen = 0;
int64_t n_cand_flen_gt_tlen_overhang_zero = 0;
double  posterior_mass_flen_gt_tlen = 0.0;

// Per top-K transcripts
struct {
    int32_t t_idx;
    double  posterior_mass_flen_gt_tlen;
} top5_overhang_offenders[5];
```

Surface in `summary.json` under a new `overhang_audit` block.

### 6.2 Sentinel test

```python
def test_overhang_audit_flags_zero_overhang_anomalies():
    """If any candidate has f_len > t_len AND overhang == 0 AND
    posterior_mass > threshold, flag in summary.json."""
```

This is a "tell us when there's a bug" canary. If it never fires on real data, hard-clipping is unnecessary.

**Acceptance:** Diagnostics land. We collect data on real samples. Hard-clip decision deferred until we see whether the canary ever fires.

---

## 7. Phase 5 — Cleanup

After Phases 1–4 are green and the contract is locked:

1. Audit `tx_starts`, `tx_ends` in the partitioned EM payload. After removing `apply_bias_correction_uniform` they should have zero remaining callers in the EM path. Delete from CSR builder, scorer output, and `_em_impl` ABI. (Coverage-weight already lives in `coverage_wts`.)
2. [`src/rigel/bias.py`](src/rigel/bias.py): `BiasProfile` has zero live callers. Either delete the file or leave a single-line stub that re-exports the canonical `effective_length_contained` as a future hook for non-uniform biases. Recommend delete; reintroduce when we actually have a non-uniform bias model.
3. Add the small public wrapper to `src/rigel/frag_length_model.py`:

```python
def effective_length_contained(
    lengths: np.ndarray,
    fl_model: FragmentLengthModel,
    *,
    min_value: float = 1.0,
) -> np.ndarray:
    """FL-marginal effective length: Σ_ℓ fl_model.pmf(ℓ) · max(L − ℓ + 1, 0).

    Canonical primitive shared by transcript TPM (mRNA) and gDNA Poisson
    densities. Thin wrapper over `compute_all_transcript_eff_lens` so the
    math is named consistently across call sites.
    """
    return fl_model.compute_all_transcript_eff_lens(lengths, min_value=min_value)
```

Then change `src/rigel/calibration/density_global.py::l_eff_contained` to one-line forward to it. (Pure-rename; no behavior change.)

**Acceptance:** Diff is mostly deletions. Tests still green. Code surface area shrinks.

---

## 8. Phase 6 — gDNA, derived from first principles

**Goal:** Decide and implement the right gDNA effective-length geometry. **Not** in this patch.

### 8.1 Why defer

* The gDNA per-hit denominator (`t_span + flank − ℓ_f + 1`) and the per-region calibration denominator (`Σ_ℓ h_G(ℓ)·max(L_region − ℓ + 1, 0)`) measure different geometric events: the former is *anchored at one transcript span*, the latter is *contained in a region*. Whether either is correct, and which the EM should use, requires a derivation we have not done.
* The user correctly points out the right gDNA denominator should probably be MultiLocus-level over the union of intervals from `locus.py`, not an anchor-transcript span.
* Bundling this with the mRNA fix doubles the number of moving parts and risks burying the mRNA win under a gDNA regression.

### 8.2 What Phase 6 should produce

Before any code change, a focused note `docs/calibration/gdna_effective_length_2026_06.md` that:

1. Defines the geometric event the gDNA likelihood is integrating over (contained / overlapping / spanning a window).
2. Derives the FL-marginal effective length for that event.
3. Decides locus-level vs. multi-locus union geometry.
4. Decides what (if any) `flank` term remains.
5. Specifies the calibration denominator change (if any) needed to keep `(C)` and the new `(D′)` in agreement.

Then implementation mirrors Phase 1: `log_eff_len[gdna_idx]` is set per-locus from the new primitive; per-hit `−log(e_h)` is removed from `scoring.cpp`.

**Acceptance for the note (not the code):** Reviewed and signed off before any C++ changes.

---

## 9. Risk + rollback

* **Phase 1 failure mode:** goldens shift in unexpected ways → revert ABI patch (single PR), re-investigate. Scope is small enough that bisection between Phase 0 (red tests, no code change) and Phase 1 is trivial.
* **Phase 2 failure mode:** assignment changes break a downstream count check → fix is mechanical (keep the new normalization, rebaseline the count check).
* **Phase 4 canary fires loudly:** treat as discovery, not regression — analyze the offending transcripts before any clip decision.

---

## 10. Test inventory (new)

| Test | Phase | Purpose |
|------|------:|---------|
| `test_em_responsibility_follows_Ltilde_ratio` | 0/1 | Core invariant (I1) |
| `test_assignment_uses_same_log_eff_len_as_em` | 0/2 | Invariant (I2) |
| `test_overhang_candidate_is_kept_and_penalized` | 0/4 | Invariant (I4) |
| `test_no_per_fragment_length_correction_in_log_lik` | 0/1 | Invariant (I3) |
| `test_short_unexpressed_tx_no_phantom_count_in_gdna_none` | 0/1 | Synthetic-sim regression |
| `test_exon_intron_straddler_does_not_enter_intron_body_density` | 0 | Pin calibration accumulator semantics |
| `test_overhang_audit_flags_zero_overhang_anomalies` | 4 | Canary for hard-clip decision |
| `test_legacy_em_uses_effective_lens_only` | 3 | Legacy-path consistency |

---

## 11. Why this plan is better than the previous draft

* **No premature hard-clipping.** Overhang remains a model property, penalized by `overhang_penalty`. Diagnostics land first; clipping decision waits for evidence.
* **gDNA is its own project.** No flank guesswork mixed into the mRNA fix.
* **Posterior assignment is treated as part of the contract**, not an afterthought. Without Phase 2, EM and assignment optimize different objectives.
* **Notation is honest:** `log h_R(f_len_candidate)` is candidate-specific in rigel because `f_len` is a transcript-space projection that varies across candidate isoforms. The shared piece is the *normalization rule*, not the FL term.
* **Calibration accumulator is not over-claimed:** the audit-corrected statement is "per-region counts use the OR'd mask; downstream densities filter to pure-mask columns." The pinned test makes this concrete.
* **The `BiasProfile` abstraction is not promoted prematurely.** A small typed wrapper near `frag_length_model.py` is enough until we actually have a non-uniform bias model.
* **Each phase is independently mergeable.** Smallest possible diffs, smallest blast radius if something goes wrong.

---

## 12. Open questions to resolve in Phase 0

1. **Exact ABI surface.** Confirm the current keyword name in `batch_locus_em_partitioned` (likely `t_lengths` or `profile_lengths`); decide whether to replace in place (preferred) or add a new arg with deprecation.
2. **`L̃_t` for synthetic nRNA spans.** They have a transcript length too, but is that the right "L" for the FL-marginal integral, or do we want the genomic span? Open in Phase 0; verify the existing nRNA tests still pass after Phase 1.
3. **Floor value for `L̃_t`.** Default `min_value=1.0` matches salmon; do we want a configurable per-tool floor (e.g., `mean(h_R)/4`) for the TPM step on very short transcripts? Defer to Phase 4 with diagnostic data.
4. **Coverage weight `cov_wt`.** Independent of length correction, but it currently rides in `mc.coverage_wt` and gets multiplied into `log_lik` indirectly. Confirm this stays untouched in Phase 1.

These are tractable questions for the Phase 0 design freeze, not blockers for the overall direction.
