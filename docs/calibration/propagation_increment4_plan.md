# Increment 4 — the count `var~mean` coupling precision + wiring to production

**Status:** execution plan, 2026-06-13. Companion to
[`propagation_message_passing.md`](propagation_message_passing.md) (the theory + roadmap) and
[`propagation_simplex_plan.md`](propagation_simplex_plan.md) (rev-2). Increment 3 built
`propagate_simplex` (forward–backward gDNA-density propagation) with **two flagged placeholders**:

1. the RTS **process variance `Q`** (per-hop coupling variance) defaulted to the *median local
   observation variance* — a global scalar;
2. `propagate_simplex` is **not wired into `calibrate`** (standalone, behind no flag).

Increment 4 (a) replaces `Q` with the **non-parametric `var~mean`** model we already built, making the
coupling precision — and hence the propagation decay-with-distance — data-driven; (b) wires
`propagate_simplex` into `calibrate.use_propagation`; (c) validates on the flagship + the full gDNA suite
with **illustrations of the fit**; (d) flips the default and retires the count-cascade if it wins.

---

## 1. Reuse what exists — the LOESS `var~mean` is already implemented

In [`density_model.py`](../../src/rigel/calibration/density_model.py):

- **`_loess(x, y, xq, frac, robust_iters=2)`** (lines 93–130) — robust Cleveland local-linear LOESS
  (tricube kernel, bisquare reweight). Handles the on/off-target chasm and heavy-tailed variance.
- **`_count_fraction_variance(...)`** (lines 133–186) — fits `log σ²_density ~ log μ_ρ` on the **paired
  boundary disagreements**: each 2-anchor node gives a density-variance point `¼(d_L−d_R)²` at mean
  `½(d_L+d_R)` (lines 165–167); the LOESS curve is read at each node's own density (line 176); returns the
  fraction-space `σ_g²` (Poisson-floored, Bernoulli-capped).
- **`node_gdna_density(...)`** (lines 188–313) already computes per-node `density`, `d_left`, `d_right`,
  `n_anchor`, and `count_gdna_frac` — everything the curve needs.
- Params: `_LOESS_SPAN=0.4`, `_LOESS_MIN_FIT=5`, `_BISQUARE_C=6.0` (lines 88–90; the standard LOWESS
  values, not tuned). Tests: [`test_count_posterior.py`](../../tests/calibration/test_count_posterior.py)
  (line fit, robustness, the on/off chasm, Poisson floor, Bernoulli cap). Validation script:
  [`scripts/debug/diag_variance_mean.py`](../../scripts/debug/diag_variance_mean.py) (scatter + binned
  trend, log-log, confirms `var∝mean²`).

**Key insight (what the paired disagreement measures).** `d_L` and `d_R` are two estimates of the *same
node's* gDNA density from its two boundary crossings. Their disagreement `¼(d_L−d_R)²` is the
density's **spatial variation across one node's span** (plus crossing-count Poisson noise) — which is
exactly the **per-hop process variance** the RTS smoother needs (how much `ρ_g` changes between adjacent
same-class nodes). So the curve we built for the count *observation* variance is also the right estimator
for the propagation *process* variance. One model, two uses.

---

## 2. The state-space model and where `var~mean` maps

The RTS smoother (`simplex_propagate._rts_smooth`) is a scalar Kalman/RTS on `ρ_g` along each
`(reference, enrichment-class)` chain:

```
state:        ρ_g,i = ρ_g,{i−1} + w_i ,   w_i ~ N(0, Q_i)     (Q_i = per-hop process variance)
observation:  y_i   = ρ_g,i     + v_i ,   v_i ~ N(0, R_i)     (R_i = seed measurement variance, r_i = 1/R_i)
```

- **`R_i` (observation noise)** — already correct from increment 3: a seed's measurement variance
  (intergenic Poisson `U/L²`; single-strand `(U/L)²·var_g` from the strand posterior). Keep.
- **`Q_i` (process noise)** — the increment-4 change: from the `var~mean` curve, `Q_i = σ²_density(μ_i)`
  evaluated at node `i`'s density `μ_i` (the count module's local density, available pre-propagation — no
  circularity). Small `Q` where density is locally stable (a confident seed carries far); large `Q` where
  the local density is volatile (decay fast, trust only the nearest seed). This makes precision-decay
  **derived**, not the global placeholder.

`Q` from the curve is a slight **over-estimate** of the pure process variance (it folds in crossing-count
measurement noise), which errs toward conservative/local propagation — safe; tightening it is future work
(issue 8.2).

---

## 3. Concrete changes

### 3.1 Expose the density-space variance curve (DRY refactor)
In `density_model.py`, factor the LOESS fit out of `_count_fraction_variance` into a reusable:

```python
def density_variance_curve(d_left, d_right, density, *, frac=_LOESS_SPAN) -> np.ndarray:
    """Per-node σ²_density via the robust log-log LOESS on paired boundary disagreements
    (the fit half of _count_fraction_variance). Returns NaN where the curve is undefined
    (<_LOESS_MIN_FIT two-anchor points or density≤0)."""
    # steps 1–2 of _count_fraction_variance: two_anchor mask, μ_d, raw_var, log-log _loess at density
```
`_count_fraction_variance` then calls it and divides by `density²` for `v_rel` (behaviour byte-identical;
pin with the existing `test_count_posterior.py`). This is the single source of `σ²_density`.

### 3.2 Per-node `Q` in the smoother
Generalize `simplex_propagate._rts_smooth(y, r, process_var)` to accept `process_var` as **either a scalar
or a per-node array** (`Q[i]` is the noise on the hop *into* node `i`). The predict step uses `p + Q[i]`.
Trivial change; keeps the scalar path (and the increment-3 tests) working.

### 3.3 Derive `Q` in `propagate_simplex`
- Compute `σ²_density` once via `density_variance_curve(d_left, d_right, density)` (from a
  `node_gdna_density(..., need_count_variance=False)` call — we need `density`, `d_left`, `d_right`, which
  it already computes; expose them or recompute the three cheaply).
- Per chain, set `Q[i] = σ²_density(μ_i)`; where the curve is NaN (too few anchors), fall back to the
  increment-3 median-observation-variance scalar (so a chain with no `var~mean` data still propagates).
- **Optional distance scaling (flagged, off by default):** scale `Q[i]` by the genomic gap to the previous
  same-class node over a typical node span. Start unscaled (per-hop); evaluate in §5.

### 3.4 Optional: count-observable boundaries as extra observations
Increment 3's seeds are intergenic + single-strand regions. A **count-observable boundary side** (pure
gDNA crossing, `boundary_count_observable`) is also a direct `ρ_g` observation. Adding these as RTS
observations (with their crossing-Poisson `R`) enriches the seed set under capture (where intergenic is
depleted). Include if §5 shows the flagship needs it; otherwise defer.

### 3.5 Wire into `calibrate`
In `calibrate.py`, the `if config.use_propagation:` branch currently calls the count-cascade
`propagation.propagate`. Switch it to `simplex_propagate.propagate_simplex`, passing the gDNA region /
boundary-side effective lengths already in scope (`reg_el`, `bnd_el`) and `rna_sense_frac`, `od_g`, `od_r`,
`gdna_prior_count` (add the last to `CalibrationConfig` with the placeholder default; see
[no-magic-numbers]). Keep `propagation.py` importable but unused (delete in increment 6 after the suite
validates). The output `(regions, left, right)` interface is unchanged, so `derive`/`priors`/`result`
need no change.

---

## 4. Illustrations (view the fit on the data)

A new script **`scripts/debug/plot_simplex_var_mean.py`** (reuse the matplotlib patterns in
`scripts/debug/diag_variance_mean.py` and `diag_count_variance_intuition.py`; the palette/style helpers in
`scripts/benchmarking/report.py`). Run on a flagship condition; oracle truth via `parse_origin` exactly as
`scripts/debug/propagation_flagship_dissect.py` does. Four panels saved to PNG:

- **A — the `var~mean` fit.** Scatter of paired-anchor `(μ_density, σ²_density)` in log-log, **coloured by
  enrichment class** (on-target exon vs off-target), with the LOESS curve overlaid. Shows the fit quality
  and the on/off chasm — the picture you want to eyeball as we test span/robustness.
- **B — the derived `Q` field.** Per-node `Q = σ²_density(μ)` vs genomic order (or vs density), to see the
  process-noise the propagation uses (where it decays fast vs carries far).
- **C — propagation accuracy.** Smoothed `ρ_g` vs **oracle** gDNA density per region (identity line),
  coloured by strand class — does the density field match truth, especially at AMBIG nodes.
- **D — `f_g` calibration.** Predicted region `f_g` vs oracle gDNA fraction (identity line), coloured by
  strand class; single-exon vs multi-exon marker. The bottom-line accuracy plot.

CLI: `python scripts/debug/plot_simplex_var_mean.py [condition] [out.png]`, defaulting to the flagship
(`gdna_gdna300_ss_0.99_nrna_none_capture_on`) and a path under the suite. Keep it a standalone diagnostic
(matplotlib is not a package dep).

---

## 5. Validation protocol

1. **Unit.** Extend `test_simplex_propagate.py`: a per-node-`Q` RTS test (array path matches scalar when
   `Q` constant; a high-`Q` node decays faster than a low-`Q` node). A `density_variance_curve` test
   (recovers `var∝mean²` on synthetic two-anchor data; NaN below `_LOESS_MIN_FIT`). Keep
   `test_count_posterior.py` green (the DRY refactor is behaviour-preserving).
2. **Flagship dissection.** Add a `--simplex` path to
   `scripts/debug/propagation_flagship_dissect.py` (or a sibling) calling `propagate_simplex`; confirm
   region **4030271 → ~0.94 without any floor** and the gross leak drops **below the count-cascade's −4%**.
   Generate the §4 illustration on this condition.
3. **Full gDNA suite.** Run the `calibration-benchmark` skill with `--use-propagation` on the
   20-condition suite (`~/Downloads/rigel_runs/gdna_benchmark_5mb`). Compare **net `gdna→rna` leak**,
   **pool gDNA fraction**, and **`gdna_none` false-positives** against (i) the current default and (ii) the
   count-cascade. Localize any regression with `net_flow_per_locus.tsv`.
4. **No-regression.** Single-strand / simple / no-gDNA conditions must match the current behaviour
   (the chain reduces to the strand result); mass conserved by the simplex.

---

## 6. Acceptance criteria (flip the default if all hold)

- Flagship: 4030271 ≈ oracle (~0.94) with no floor; gross/net leak < count-cascade and < current default.
- Suite: net `gdna→rna` leak ≤ current default across conditions; **`gdna_none` FP not worse**; stranded
  conditions' FP not worse; pool fraction within accepted deviation.
- Single-strand conditions unchanged (no-regression).
- All unit + golden tests green.

If met: flip `CalibrationConfig.use_propagation` default to `True` (or make `propagate_simplex` the
calibrate path unconditionally), then increment 6 retires `propagation.py` (count-cascade + floor),
`run_fill` (→ the worklist), and the `strand_deconv` per-node blend.

---

## 7. Sequencing (small, reviewable commits)

1. DRY refactor `density_variance_curve` + test (behaviour-preserving). 
2. Per-node `Q` in `_rts_smooth` + test.
3. `Q` from the curve in `propagate_simplex` + test.
4. The illustration script.
5. Flagship dissection (+ `--simplex`); read the plots; tune `_LOESS_SPAN` only if the fit visibly
   demands it (CV, not by hand).
6. Wire into `calibrate` (`gdna_prior_count` config field); full suite; compare.
7. (If criteria met) flip the default.

---

## 7b. VALIDATION FINDING (2026-06-13) — issue A confirmed, blocks the wiring

Steps 1–3 are implemented + tested (32 green). The illustration
(`scripts/debug/plot_simplex_var_mean.py`) on the flagship
(`gdna_gdna300_ss_0.99_nrna_none_capture_on`) surfaced a real problem **before** wiring into `calibrate`:

| AMBIG class | n | oracle gDNA frac | count clue (`density_model`) | propagate_simplex |
|---|---|---|---|---|
| **EXON** | 87 | **0.946** | **0.246** | (carries the clue) |
| INTRON | 54 | 1.000 | 0.852 | — |
| all AMBIG | 141 | 0.967 | 0.478 | **0.381** |

**Root cause = issue A (the make-or-break), confirmed.** AMBIG **exon** regions are *not* count-observable,
so their gDNA density is imputed from boundary crossings — which cross the **capture enrichment
discontinuity** into low-density off-target introns, dragging the imputed on-target density to ~¼ of
truth. The simplex and the propagation are **not** at fault (they faithfully carry the density they are
given; the machinery's own guarantees — order-independence, no over-subtraction, AMBIG-reads-RNA-when-no-
gDNA — all hold). The defect is upstream, in the **count gDNA density under capture**: there is no reliable
*enriched-exon* gDNA-density source. The strand-derived density at expressed single-strand exons is too
low/noisy (the small gDNA fraction, strand-measured), and the count imputation smears in off-target
introns.

**This is the decision point flagged as issue A.2 / 8.4.** The fix needs a capture-aware enriched-exon gDNA
density — options: (a) **capture-normalized density** (divide a per-node capture factor out, so the chain
couples in off-target-equivalent units and the off-target seeds *predict* on-target density); (b) restrict
the exon-class chain to true exon-class anchors with correct precision (the count-observable **exon-edge
boundary** pure-gDNA crossings, §3.4 — not optional under capture); (c) port the count mean-bias
debiasing (`count_mean_bias_design.md`, Phase 4-mean) that targets exactly this exon gDNA-count
under-call. **Pause for direction before wiring** — this is a joint modelling decision, not a silent patch.

## 8. Open issues / risks

1. **Circularity (issue F).** `Q` uses the count module's density `μ_i`, available *before* propagation —
   so one pass is non-circular. The outer loop (refit `var~mean` from the *propagated* deconvolution,
   re-propagate) is the sound EM-style extension if one pass under-fits; out of scope for increment 4.
2. **`Q` over-estimate.** The paired disagreement includes crossing-count Poisson noise, so `Q` slightly
   over-states the true process variance → conservative (local) propagation. Acceptable; a
   measurement-noise subtraction (`Q = max(0, σ²_density − Poisson)`) is a future tightening.
3. **Distance scaling.** Per-hop `Q` ignores the genomic gap between same-class neighbours (an exon and the
   next exon skip an intron). Start unscaled; §5-2's plots show whether long skips need distance-scaled `Q`.
4. **Enrichment handling.** Increment 3's class-chains are kept (parameter-free). The
   capture-normalized single-chain alternative (divide out a per-node capture factor) is deferred unless
   the class-chains leave too few seeds per chain (watch in §5-3).
5. **Seed sparsity under capture.** Intergenic gDNA is depleted under capture; if the off-target chain is
   seed-starved, add the count-observable boundaries (§3.4).
6. **`gdna_prior_count` derivation.** Still the ~0.5-count placeholder (plan §5 of the simplex doc); its
   derived form (the global gDNA expectation) remains future work, tracked separately.
7. **Perf.** Per-node-`Q` RTS is still O(n) per chain; the `var~mean` LOESS is O(R²) (as today, gated). Fine
   for the suite; the Gaussian-message / vectorization work is phase 7.
8. **No-magic-numbers.** New constant introduced: `CalibrationConfig.gdna_prior_count` (placeholder,
   documented). `_LOESS_SPAN` reused as-is (CV-selectable). No other new constants — `Q` and `R` are
   derived. Flag the `gdna_prior_count` default for the user before wiring (§3.5).
