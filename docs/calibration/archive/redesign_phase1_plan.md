# Redesign Phase 1 — strand module as a standalone likelihood emitter (dry-run plan)

**Status:** dry-run implementation plan for review, 2026-06-11. Phase 1 of
`sequential_calibration_redesign.md`. Built **alongside** the live blend — the live calibration path is
**unchanged and golden stays bit-identical**; nothing is wired into `calibrate.py` this phase.

## Goal

A pure, standalone strand deconvolver: for every node (region-contained view + each boundary side),
emit the gDNA fraction read at the FP-quantile, the resulting unspliced gDNA/RNA **mass** split, and the
**likelihood information** `I` that the count module will use as precision. No count input, no blend.

## What it consumes / produces

```
strand_deconvolve(substrate, region_arrays, *, rna_sense_frac, gdna_strand_overdispersion,
                  rna_strand_overdispersion, deconv_quantile, n_grid)
  -> (contained: StrandSplit, left: StrandSplit, right: StrandSplit)

@dataclass StrandSplit:        # all float64[R]
    gdna_frac   # g_q in [0,1] where strand-informative; NaN where not (AMBIG / no sense)
    gdna_mass   # g_q · mass_unspliced        (NaN where not split)
    rna_mass    # (1−g_q) · mass_unspliced    (UNSPLICED RNA only; spliced is the count module's)
    info        # I = N·(2·rna_sense_frac−1)²  ; N = unspliced count (pos+neg); 0 where not informative
```

## Decisions (every one, grounded in the current code)

**D1 — `info` is the *likelihood* Fisher information `I = N·(2·ss−1)²`, not the BB posterior variance.**
This is the whole point of the redesign (`sequential_calibration_redesign.md` §1.1): at `ss=½` the BB
*posterior* variance is the finite *prior* variance (0.125) — using it would give a silent node weight.
`I = N·(2·ss−1)²` → 0 at `ss=½` (and is set to 0 at AMBIG). `N` is the unspliced **count**
(`n_unspliced_pos + n_unspliced_neg`), not the fractional mass. (`strand_posterior_gdna_frac` keeps
returning its posterior `variance` for the *live* blend path; `strand_deconvolve` ignores it and computes
`I` itself.)

**D2 — the split applies to the unspliced *mass*; `I` comes from the *count*.** `gdna_mass = g_q ·
mass_unspliced`. For the contained view `mass_unspliced == count` (the accumulator stores
`region_contained` as both, `substrate.py:103`); for boundary sides `mass_unspliced` is the fractional
crossing mass and the count is the flux `n_unspliced_pos+neg` — so the two genuinely differ on sides and
must be kept separate. The **spliced** mass is *not* touched here — it is deterministic RNA the count
module adds in the 3-term.

**D3 — the FP-quantile reads the strand posterior CDF at `q`.** Generalise the existing batched
`_grid_posterior_median(post, grid)` → `_grid_posterior_quantile(post, grid, q=0.5)` (replace the `0.5`
constants with `q`; identical at `q=0.5`). `strand_posterior_gdna_frac` gains a `deconv_quantile=0.5`
param and returns `g_q` at that quantile. **This is the clean "strand-deconvolution FP-rate" knob**
(benefit #1) — a genuine posterior quantile of the *strand* split, not the live blend's Gaussian
`center + Φ⁻¹(q)·σ` shift. Default `q=0.5` ⇒ the strand median.

**D4 — `g_q` is computed only where the node is strand-informative; else `NaN`, `I=0`.** A node is
informative iff it is **strand-observable** *and* `(2ss−1)²>0`. Reuse the existing observability logic:
contained regions are observable iff `strand_class ∈ {POS, NEG}` (`deconv_regions`); boundary sides via
`_compute_side`'s `strand_observable` (which already encodes the `TS_NONE` intergenic-wildcard rule). At
non-informative nodes the strand cannot split → emit `gdna_frac = NaN`, `I = 0`. **The count module
(Phase 2) splits these from the density field** — the strand module does not guess. (This keeps mass
conservation a Phase-2 responsibility for `I=0` nodes; Phase-1 conserves it where `I>0`:
`gdna_mass + rna_mass = mass_unspliced`.)

**D5 — orientation: extract the side helper, inline the contained one.** For sides, factor the
sense/antisense/strand_observable computation out of `_compute_side` into `_side_strand_orientation(view,
same, ts_self, ts_other)` (golden-identical refactor; `_compute_side` then calls it and appends its
count-fraction). For the contained view the orientation is the two-liner already in `deconv_regions`
(`sense = where(ts==NEG, neg, pos)`, `observable = (ts==POS)|(ts==NEG)`) — inline it in
`strand_deconvolve`. Avoids both duplication and passing `_compute_side` dummy `eff`/`region_density`.

**D6 — `strand_deconvolve` is built ALONGSIDE; the live blend is untouched.** It is a new function +
dataclass in `strand_deconv.py`; `deconv_regions`/`deconv_sides`/`_deconv_per_node` and `calibrate.py`
are unchanged. The only edits to live code are the golden-identical generalisations (D3: median→quantile;
D5: extract the side helper). **Golden must stay bit-identical** — verified by the full suite.

**D7 — `I` uses the symmetric-point Fisher info `N·(2ss−1)²`.** Exact for a balanced node (`p=½`); for an
unbalanced node the exact Fisher info is `N·(½−ss)²/(p̂(1−p̂))`, larger by `0.25/(p̂(1−p̂)) ≥ 1`. So
`N·(2ss−1)²` is a **conservative lower bound** on the information (slightly under-weights confident nodes)
— clean, closed-form, magic-number-free, and `= N·w` (the discriminability). The exact per-node form is a
flagged refinement (I1).

## Issues / risks uncovered in the dry-run

- **I1 — `info` approximation (D7).** `N·(2ss−1)²` ignores the `p(1−p)` curvature factor. Likely fine
  (the field aggregates many nodes; the bound is conservative), but flag for the Phase-2/3 validation: if
  confident unbalanced nodes are under-weighted in the field, switch to the exact per-node Fisher info
  (available from the same grid). **Do not fix pre-emptively** — measure.
- **I2 — quantile semantics change (D3).** The new knob reads the *strand* posterior; the live knob shifts
  the *blend* center. They are not the same number for `q≠½`. This is intended and harmless in Phase 1
  (not wired in), but Phase 3 must define exactly how `q` (on the strand likelihood) interacts with the
  field in the combine — noted there, not here.
- **I3 — `NaN` at `I=0` (D4).** `gdna_mass`/`rna_mass` are `NaN` at non-informative nodes. Phase 2's
  count module must treat `I=0` (or `isnan(gdna_frac)`) as "impute from the field," and own mass
  conservation for those nodes. Phase 1 tests assert the `I>0` nodes conserve mass and the `I=0` nodes are
  flagged, nothing more.
- **I4 — `strand_posterior_gdna_frac` at AMBIG.** It must only be *called* on strand-observable nodes
  (it needs a defined sense). `strand_deconvolve` masks to `informative` before calling it (mirrors
  `_deconv_per_node`'s `use_strand` gate) — never feed it an AMBIG node.
- **I5 — `(2ss−1)²` is global, `N` is per-node.** `ss = rna_sense_frac` is the one library scalar; `I`'s
  per-node variation is entirely `N`. So a deep AMBIG region still gets `I=0` (observability gate), and a
  shallow confident region gets small `I` (few fragments) — both correct.

## Files & changes

- `calibration/strand_deconv.py`:
  - `_grid_posterior_median` → `_grid_posterior_quantile(post, grid, q=0.5)` (generalise; live callers
    pass nothing ⇒ `0.5`).
  - `strand_posterior_gdna_frac(..., deconv_quantile=0.5)` → return `g_q` at the quantile (median at
    default). Live `_deconv_per_node` unaffected (calls without `deconv_quantile`).
  - extract `_side_strand_orientation` from `_compute_side`; `_compute_side` calls it (golden-identical).
  - **new:** `StrandSplit` dataclass + `strand_deconvolve(...)`.
  - export `StrandSplit`, `strand_deconvolve` in `__all__`.
- **No** change to `calibrate.py`, `density_model.py`, `splice_junction.py`.

## Test plan (`tests/calibration/test_strand_deconvolve.py`, new)

1. **Confident split.** Strong `ss` (e.g. 0.99), balanced strand-observable region ⇒ `g_q ≈ 1`,
   `I = N·(2·0.99−1)² > 0`, `gdna_mass ≈ mass_unspliced`, `rna_mass ≈ 0`.
2. **RNA node.** Strand-observable, sense-skewed to `ss` ⇒ `g_q ≈ 0`.
3. **AMBIG ⇒ no split.** `strand_class = AMBIG` ⇒ `I = 0`, `gdna_frac` NaN.
4. **Unstranded ⇒ no info.** `ss = 0.5` ⇒ `I = 0` for every node (the `(2ss−1)²` factor).
5. **Mass conservation** where `I>0`: `gdna_mass + rna_mass == mass_unspliced`.
6. **Quantile monotonicity.** `g_q` non-decreasing in `q`; `q=0.5` reproduces the median.
7. **`I` formula.** `info == N·(2ss−1)²` on informative nodes, `N = pos+neg`.
8. Plus: `_grid_posterior_quantile(·,·,0.5)` equals the old `_grid_posterior_median` (regression).

## Golden / regression guarantee

The only live-code edits (median→quantile generalisation with `q=0.5` default; side-orientation
extraction) are golden-identical by construction. **Exit criterion:** full suite green + golden
bit-identical + the 8 new `strand_deconvolve` unit tests pass; `strand_deconvolve` is callable and
correct in isolation but not yet consumed by `calibrate`.
