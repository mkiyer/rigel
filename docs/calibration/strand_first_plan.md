# Strand-first calibration — the converged implementation plan

**Status:** the plan to implement, 2026-06-11. Converged after the full round-trip. This supersedes the
rev-1/rev-2 Phase-2 sketches and the precision-weighted-field heuristic. It is a **targeted refinement of
the existing decoupled blend**, not a rewrite — Phase 1 (`strand_deconvolve`) is already built and emits
the one quantity the carry-forward needs.

## The architecture (full circle, refined)

Per region, the gDNA fraction comes from a **strand-first** estimate, supplemented by the **count model**
only as the strand's confidence falls:

```
g_region = w · g_strand  +  (1 − w) · g_count          # the combine
```

- `g_strand` — the region's own strand deconvolution (capture-invariant: a *fraction*, so probe enrichment
  cancels). The primary signal.
- `g_count` — the count model's estimate: the gDNA density imputed from the **strand-deconvolved
  boundaries** (the carry-over / boundary sweep), built on **strand-cleaned** counts (the 3-term). The
  supplement, used where the strand is weak or silent (AMBIG, unstranded, thin).
- `w` — **the carried strand uncertainty** (below). It routes: `w→1` ⇒ strand governs; `w→0` ⇒ count
  governs.

## The solved problem — carrying the strand uncertainty forward

The strand deconvolution emits, per node, the gDNA fraction **and its information**
`I = N · (2κ − 1)²` — the Fisher information per fragment `(2κ−1)²` times the number of unspliced
observations `N` (this is exactly what Phase-1 `strand_deconvolve` already returns as `info`). The count
model converts that information into the combine weight by a **saturating** function:

```
w  =  I / (I + I₀)            I = N·(2κ−1)² ,   I₀ = a small strand-information scale (default 1)
```

This is the whole carry-forward, and it gives the three behaviours by construction:

1. **Unstranded (κ = ½):** `I = 0` ⇒ `w = 0` everywhere ⇒ the strand is a **no-op**; the count model does
   all the work (boundary imputation on raw counts — the floor). ✓
2. **Excellent stranding (κ → 1, decent `N`):** `I` is large ⇒ `w → 1` ⇒ the strand governs each region.
   The count model only **supplements** the regions the strand couldn't see — AMBIG regions, whose own
   `I=0`, take `g_count` from the **carry-over of the strand-deconvolved single-strand boundaries**. And
   the strand's **3-term cleaning** (splitting unspliced into DNA vs RNA) feeds those boundaries, so the
   imputation is clean. ✓
3. **Weakly stranded (κ near ½, or thin `N`):** `0 < w < 1` ⇒ a sensible blend. Not claimed optimal — and
   it doesn't need to be: real libraries are bimodal (≫90%, often ≫99%, *or* essentially unstranded), so
   the in-between is rare. A saturating `w` that behaves sensibly at the extremes is enough. ✓

**Why `I/(I+I₀)` and not `(2κ−1)²`.** Your objection is right: `(2κ−1)²` ignores `N` and over-penalises
good data — `0.99 → 0.96`, dinging an excellent strand model. `I/(I+I₀)` uses the *total* information
`N·(2κ−1)²`, so at κ=0.99 with any real depth (`N` a few dozen) `I` is large and `w → ~1` — it stops
fighting a strand model that is clearly good. It also naturally **distrusts thin regions** (low `N` ⇒ low
`I` ⇒ lean on the boundary imputation), which is exactly right. `I₀` is one interpretable scale ("the
strand information at which we half-trust the strand," ~1 effective discriminating fragment) — a sensible
default to **validate on the suites, not theorise**.

## What changes from the existing code (targeted)

The live calibrator already has this shape: `strand_deconv._deconv_per_node` does `w·g_strand +
(1−w)·g_count`, and `node_gdna_density` builds `g_count`. Two refinements:

1. **The weight.** Replace `w = (2κ−1)²` with `w = I/(I+I₀)`, `I = N·(2κ−1)²`. Per-node (uses `N`),
   saturating. `strand_deconvolve` (Phase 1) already produces `I` (`info`).
2. **The count field.** `node_gdna_density` builds `g_count` from **strand-cleaned** counts (the 3-term:
   subtract the strand-identified RNA from each node before imputing density), and imputes AMBIG /
   unstranded regions by the **carry-over of strand-deconvolved boundaries** (the boundary sweep). Today
   it uses raw counts; cleaning them is the AMBIG fix.

`derive`, `priors`, the EM, the FL/overdispersion fits, and Phase 1 are intact.

## Phases (push-forward — plan-then-implement each)

- **Phase 1 — DONE** (`79289e1`): `strand_deconvolve` emits per node `(g_strand at the FP-quantile,
  cleaned gDNA/RNA split, info `I = N·(2κ−1)²`)`. Golden-identical, alongside.
- **Phase 2 — the count field on strand-cleaned counts + the carry-over.** `node_gdna_density` consumes
  the strand split + `I`; imputes each region's gDNA from its strand-deconvolved boundaries; carry-over
  for AMBIG. **Validate on the toy AMBIG scenario** (`count_gdna_frac[AMBIG]` 0.120 → ~0 at gDNA=0+nascent
  in stranded data; smooth degradation toward the floor as κ→½). Built alongside; golden-identical.
- **Phase 3 — the combine + the weight, wired into `calibrate`.** Replace `_deconv_per_node`'s weight with
  `w = I/(I+I₀)`; route `g = w·g_strand + (1−w)·g_count`; retire the old blend. Regenerate golden.
  **Validate the three behaviours** (unstranded → count-only; ss0.99 → strand-dominated, `w≈1`;
  weakly-stranded → sensible blend) and the suites (no regression; AMBIG/nascent improved).
- **Phase 4 — lock + cleanup.** Promote the toy scenario to a test/golden; full 3-pool suite validation;
  remove dead blend code; update theory docs.

## Sensible defaults (validate, don't theorise)

- `I₀` — the strand-information half-trust scale. Default **1**. The bimodal-data argument says the exact
  value barely matters (almost everything is `w≈0` or `w≈1`); tune on the suites if a weakly-stranded
  condition is sensitive. A sharper saturation (`Iⁿ/(Iⁿ+I₀ⁿ)`, `n>1`) is available if validation wants it.
- The combine stays a **blend**, not a hard route — but because `w` saturates and the data are bimodal,
  it behaves like a route in practice, while degrading smoothly through the rare middle. This also keeps
  the capture-biased count from ever dominating a confident strand (high κ ⇒ `w≈1`), so the
  bias-laundering worry is confined to the rare weakly-stranded regime, which we accept won't be perfect.

## One-paragraph summary

Strand-deconvolve each region (capture-invariant fraction); carry its information `I = N·(2κ−1)²` forward
as a saturating weight `w = I/(I+I₀)`; let the count model — built on the strand-cleaned 3-term counts and
the carry-over of strand-deconvolved boundaries — supply `g_count` for the regions the strand can't see;
combine `g = w·g_strand + (1−w)·g_count`. Unstranded ⇒ `w=0`, count-only. Excellent ⇒ `w≈1`, strand-led
with count filling the AMBIG gaps. Weakly stranded ⇒ a sensible blend we don't over-engineer. Phase 1 is
built; Phase 2 (clean count field + carry-over) is next.
