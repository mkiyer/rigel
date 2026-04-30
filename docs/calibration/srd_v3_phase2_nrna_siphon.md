# SRD v3 Phase 2 — nRNA siphon: root cause and partial fix

**Date:** 2025  
**Status:** Phase 2 partial fix landed; residual ~50–80% siphon remaining
in degenerate-likelihood loci pending a C++ allocator change (P3).

## Summary

The "nRNA siphon" is the dominant production failure mode in rigel
quantification (per `docs/calibration/srd_v3_phase2_oracle_findings.md`):
intronic gDNA fragments — which the model should classify as gDNA — are
absorbed by the synthetic nRNA component, causing nRNA-pipeline counts
to balloon and (for stranded data) drag mRNA totals down by up to ~37%.

This document reports a focused diagnostic effort that:

1. Built a minimal synthetic harness reproducing the siphon at single
   transcript / two-component scale
   ([scripts/debug/nrna_siphon_minimal.py](../../scripts/debug/nrna_siphon_minimal.py)).
2. Isolated **two** independent root causes, neither of which is the
   harmonic-mean effective-length issue we initially suspected.
3. Implemented a partial fix in
   [src/rigel/locus.py](../../src/rigel/locus.py) that resolves the dominant
   contributor (per-locus prior dilution) without breaking the 914-test
   scenario suite.
4. Identified one remaining structural issue (C++ side
   `α_RNA` allocator) that limits how aggressively the prior fix can be
   tuned.

## Diagnostic chain (E1–E7)

The harness reproduces a **single‐transcript locus** (one mRNA + one
synthetic nRNA) under controlled conditions:

| Exp | Setup | Outcome |
|-----|-------|---------|
| E1 | Pure mRNA, stranded, no gDNA | EM correctly recovers all to mRNA |
| E2 | Per-fragment effective Δ(nRNA − gDNA) for SS=0.99 | sense `+0.703`, anti `−3.892` |
| E3 | Pure intronic-pool gDNA, stranded SS=0.99 | All sense fragments → nRNA (38.9% siphon) |
| E4 | E3 + isolate effective-length contribution | gDNA per-hit `e_h = span+gdna_flank−frag+1` ≈ nRNA `e_n = span−frag+1`, **near-degenerate likelihood** |
| E5 | Mixed locus, **unstranded**, varying `α_gDNA` directly | `α_gDNA = 0.5 → 98% siphon`; `α_gDNA = 50 → 0.5% siphon` |
| E6 | Mixed locus, **stranded SS=0.99**, intronic-only `c_base` | Whole-locus `c_base=5`: 38.9% siphon → intronic-only `c_base=N_intronic`: 0.0% siphon |
| E7 | Compare **per-locus prior derivation** strategies | See table below |

### Two independent root causes

**(1) Strand-asymmetric nRNA pull** (stranded data only).  Per-fragment
likelihood for a *sense* intronic fragment slightly favours the nRNA
candidate over the gDNA candidate (E2: `Δ = +0.703 nats`).  This is
small, but additive over thousands of fragments it is sufficient to
shift the EM equilibrium toward nRNA whenever the prior is weak.

**(2) Per-locus prior dilution** (dominant; both stranded and
unstranded).  The SRD-v1 prior derivation
(`compute_locus_priors_from_partitions`) computed

```
α_gDNA = c_base · mean_f γ_f   (over ALL units in the locus)
```

where γ_f is the per-fragment posterior gDNA share.  Spliced and
exonic-only units have `gdna_log_lik = -inf` ⇒ `γ_f = 0`.  In a
typical locus dominated by exonic mRNA reads (10²–10⁵ such units), the
mean is dragged toward zero, so `α_gDNA` ends up ~0.1–0.5 — far below
the `α_gDNA ≈ 10–50` needed to overcome degenerate intronic likelihoods
and prevent siphoning.

E5 directly demonstrates: holding the per-fragment scoring and the
likelihood data fixed and only varying the *injected* `α_gDNA` is
enough to swing the siphon from 98% (at `α_gDNA = 0.5`) down to 0.5%
(at `α_gDNA = 50`).

## Fix (Phase 2)

`compute_locus_priors_from_partitions` now restricts the γ_f average to
**eligible units only** (those with finite `gdna_log_lik`, i.e. the
unspliced/intronic units the calibrator's π_pool actually pertains
to):

```python
eligible = np.isfinite(gdna_log_liks)
gamma_e  = sigmoid(log_prior_g + gdna_log_liks[eligible] − max_rna[eligible])
gamma_l  = float(gamma_e.mean())
α_gDNA[ℓ] = c_base · gamma_l
α_RNA[ℓ]  = c_base · (1 − gamma_l)
```

with `c_base` raised from **5.0 → 10.0** (the largest value that does
not break any scenario test in `tests/scenarios*/`).

### E7 results (synthetic harness, four scenarios)

| Scenario | Old siphon | Fix `c_base=5` | Fix `c_base=10` | Fix `c_base=100` |
|----------|-----------:|---------------:|----------------:|-----------------:|
| 9k mRNA + 1k gDNA, π=1.0 | 97.5% | 75.5% | **51.0%** | 0.1% |
| 99k mRNA + 1k gDNA, π=1.0 | 99.6% | 75.5% | **51.0%** | 0.1% |
| 9k mRNA + 1k gDNA, π=0.5 | 98.7% | 87.9% | **76.1%** | 8.1% |
| 9k + 1k, π=0.45 (calib underestimate of π=1.0) | 98.9% | 89.1% | **78.6%** | 10.3% |

The bolded `c_base=10` column is what production now ships.  The
`c_base=100` column shows that ~all of the siphon *can* in principle
be removed by the prior, but doing so triggers a **secondary
regression** in the C++ `α_RNA` allocator (see next section).

## Why we cannot just push `c_base` higher (the residual issue)

The C++ EM solver
([`em_solver.cpp` lines 820–855](../../src/rigel/native/em_solver.cpp))
distributes `α_RNA` across the locus's RNA components (mRNA + nRNA)
*proportional to per-component coverage*:

```cpp
prior_out[i] = baseline + α_RNA · coverage_totals[i] / total_rna_coverage;
```

In a pure-mRNA locus that happens to contain intronic-overlapping
fragments (i.e. unspliced reads inside annotated introns), the
synthetic nRNA component picks up some "coverage" from these reads —
which are genuine intronic-overlap candidates for nRNA — even though
the locus contains zero true nascent RNA.  When `α_RNA` is small
(c_base=10, ≈ 0–10 per locus), the spurious nRNA prior is negligible
and the EM correctly drives nRNA to ~0 by likelihood.  When `α_RNA` is
large (c_base=100, 0–100 per locus), the spurious nRNA prior becomes
material and tests like
`tests/scenarios/test_wide_intron.py::test_baseline_no_false_nrna` fail
with t1 dropping from 500 → 488 (nRNA absorbing 12 fragments).

**Why this matters:** the prior fix alone cannot drive the synthetic
siphon below ~50% without breaking pure-mRNA tests.  The remaining
50–80% siphon under degenerate likelihood requires the C++ allocator
to stop routing `α_RNA` into the synthetic nRNA component.

### Recommended P3 follow-up (C++ change)

In `compute_unified_ovr_prior` (`em_solver.cpp:820+`), distinguish
"true RNA components" (annotated mRNA isoforms) from "synthetic nRNA
component" and route `α_RNA` only to the former.  The synthetic nRNA
component should receive prior strength only from independent
calibration of the nRNA share (which is currently not estimated
per-locus — it relies on a global `pi_nrna`-like quantity that does
not exist in SRD v3).

This is non-trivial because:
- The OVR EM treats all components symmetrically by construction.
- The synthetic nRNA component is logically a *catch-all* for unspliced
  reads that fall inside annotated introns; its prior must reflect
  the locus's expected nascent-RNA fraction, which the calibrator does
  not directly estimate.
- One option: set `α_nRNA[ℓ] = 0` (let likelihood entirely determine
  it).  This would eliminate the spurious nRNA pull but may cause
  collapse-to-mRNA in genuine nRNA-rich loci.
- A safer option: estimate a per-locus nRNA prior from the spliced
  vs unspliced ratio of the eligible reads (an analogue of the
  calibrator's `π_pool` but for nRNA-vs-mRNA in the unspliced pool).

## Verification

- All 914 tests in `tests/` (excluding the pre-existing
  `tests/test_calibration.py` failure) pass with the fix.
- Synthetic harness E7 confirms ~2× reduction in siphon at
  `c_base=10`.
- Oracle benchmark re-run pending (P5).

## Files touched

- `src/rigel/locus.py` — `compute_locus_priors_from_partitions`
  rewritten with eligible-only mean; `_C_BASE_DEFAULT = 10.0`.
- `scripts/debug/nrna_siphon_minimal.py` — diagnostic harness with
  experiments E1–E7.

## Open follow-ups

| Priority | Task |
|----------|------|
| P3 | C++ `α_RNA` allocator: stop routing prior into synthetic nRNA component (see above). |
| P3 | Re-run oracle and full HULKRNA benchmarks; compare siphon magnitude before/after fix. |
| P4 | Calibrator: investigate whether `π_pool` underestimates the true gDNA share (E7 row 4 suggests it may — in the dna20m oracle data, calibrated `π_pool ≈ 0.45` vs ground-truth ≈ 1.0 for the intronic pool). |
| P4 | gDNA harmonic-mean effective length — analysis showed this contributes <0.05 nats per fragment vs the prior-dilution mechanism; defer until P3 lands. |
