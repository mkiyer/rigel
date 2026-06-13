# Redesign Phase 3 — the per-node gradient combine `g = w·g_strand + (1−w)·g_count`

## FINAL OUTCOME (as-built, shipped 2026-06-12) — the combine is KEPT, upgraded

> **The dry-run below proposed *retiring* the combine. A deep forensic investigation of the flagship
> scenario reversed that.** Retiring the combine silently dropped the `w·g_strand` term for every
> **exon** (a non-count-observable region's `count_gdna_frac` is *always* the boundary imputation, never
> its own count — so "cleaned count field == combine" held only for count-observable regions, which the
> verification mistakenly generalized). That dropped the capture-invariant own-strand signal → the
> ss0.99+capture gDNA→RNA regression. The fix is to **keep the combine**, as a **per-node gradient**:
>
> ```
> g = w·g_strand + (1−w)·g_count        w = I/(I+I₀),  I = N·(2κ−1)²,  I₀ = 10
> ```
>
> - `g_strand` — the node's own strand deconvolution (direction, capture-invariant). `g_count` — the
>   count imputation on **RAW contained + CLEANED boundaries** + the 3-term splice (magnitude). Independent
>   signals ⇒ no double-count (the original reason to retire was: cleaning the *contained* AND combining
>   double-applies `w` — so we clean only the *boundaries*; the contained stays raw).
> - `w = I/(I+I₀)` is the strand-trust **gradient** (not a cliff): the count module is used to the degree
>   the strand is uninformative. **I₀ = 10** (not the placeholder 1): the strand estimate's `SE(g)=1/√I`,
>   so it is only meaningful at `I≫1`; `I₀=1` half-trusts garbage (SE≈1) — which over-trusted a near-½
>   strand at high N (capture) → a false-gDNA knife-edge. `I₀=10` sets half-trust at `SE≈0.3`.
> - **Validated** (`quick_3to1_5mb`): flagship ss0.99+capture `−274k→−261k` (better than baseline);
>   knife-edge ss0.50+capture `+15.7k→+19`; all other conditions ≤0.1%; toy AMBIG fix holds
>   (gDNA=0+nascent → est 0.024); 1080 tests pass.
>
> The rest of this document is the superseded "retire the combine" exploration, kept for the record.

---

**Status:** dry-run for review, 2026-06-11 (SUPERSEDED — see FINAL OUTCOME above). This phase wires the
strand-cleaned count field (Phase 2) into `calibrate` and is the **first phase that changes golden
output**. Its headline decision — pressed for and **verified** below — is that the explicit strand/count
*combine* (`deconv_regions` / `deconv_sides` / `_deconv_per_node`'s `w·g_strand + (1−w)·g_count`) is
**redundant** once the count field is built on cleaned counts, and should be **retired**. The per-node
weight `w = I/(I+I₀)` becomes the **single** blend mechanism, living only in `cleaned_gdna_count`.

## The finding: the combine is redundant (verified, not asserted)

Today `_deconv_per_node` produces, for an observable region,

```
gdna_mass = center · mass ,  center = w·g_strand + (1−w)·g_count ,  g_count = clip(contained_count / mass)
          = w·g_strand·mass + (1−w)·contained_count            (g_count = contained_count/mass)
```

The Phase-2 cleaned field, with **no** combine, produces

```
gdna_mass = count_gdna_frac_cleaned · mass = clip(cleaned_count / mass) · mass
          = (w·g_strand + (1−w))·contained_count = w·g_strand·contained_count + (1−w)·contained_count
```

The two differ only in the `g_strand` term: `g_strand·mass` vs `g_strand·contained_count`. They are
**identical iff `contained_count == contained_mass`** — which is **empirically exact**: on the toy
substrate all 22 observable regions have `count/mass = 1.0000` (the contained accumulator deposits mass
1 per fully-contained fragment). So:

- The cleaned count field **already is** `w·g_strand + (1−w)·1` for observable regions — the combine's
  `center` with `g_count=1`. Running the combine *on top* would compute `w·g_strand + (1−w)·[w·g_strand +
  (1−w)]` — applying `w·g_strand` **twice**. A double-count, sharper-than-intended toward the strand.
- For **AMBIG / non-observable** regions the combine was already count-routed (`g_strand` undefined →
  `center = g_count`); feeding the **cleaned** boundary imputation makes that the clean answer (the Phase-2
  AMBIG fix). No combine needed.
- For **unstranded** (`κ=½`) `w=0` → cleaning is a no-op → the raw count floor, both ways.

**So the combine adds nothing the cleaning doesn't, and double-counts where it does act.** Retire it. The
"blend" is the cleaning weight `w = I/(I+I₀)`, applied once, per node.

## The architecture (one mechanism, one weight)

```
strand_deconvolve(substrate, κ̂)                      → (contained, left, right) StrandSplits (g_strand, I)
clean each view: cleaned_x = cleaned_gdna_count(split_x, raw_count_x, I₀)   # the single blend, weight w=I/(I+I₀)
node_gdna_density(cleaned)                            → density field + count_gdna_frac (REGION answer)
region masses:  gdna = count_gdna_frac · mass ;  rna = (1−count_gdna_frac)·mass + mass_spliced
side masses (QC): gdna = (cleaned_side/raw_side)·side_mass ;  rna = remainder        # cleaning fraction · mass
derive(region, left, right)                          → gdna_density_global (QC), gdna_geom_len (geometric)
```

`_deconv_per_node`, `deconv_regions`, `deconv_sides`, `NodeDeconv`, and the scalar `w=(2κ−1)²` are
**deleted**. The strand enters only through `cleaned_gdna_count`.

### The three behaviours fall out by construction (no combine)

1. **Unstranded (`κ=½`):** `I=0 ⇒ w=0` everywhere ⇒ cleaning is a no-op ⇒ `node_gdna_density` runs on raw
   counts — the documented count-only floor. ✓
2. **Excellent stranding (`κ→1`, decent `N`):** `w→1` ⇒ observable regions read `g_strand`; AMBIG regions
   inherit the **cleaned** boundary imputation; the count module fills exactly the gaps the strand can't
   see. ✓
3. **Weakly stranded:** `0<w<1` ⇒ a smooth blend that degrades to the floor (validated cliff-free in
   Phase 2's κ-sweep). ✓

## The sides are low-stakes (QC only)

`derive` uses the side masses **only** for `gdna_density_global` (a QC scalar); the per-locus prior flows
from the **region** masses + the geometric `gdna_geom_len`. So the side treatment is not prior-critical.
For a single clean mechanism, the sides are cleaned the same way (side gDNA mass = cleaning fraction ·
side mass). The only behavioural difference vs the retired `deconv_sides` is the `(1−w)` floor (1 here vs
the side's borrowed density), which matters only for non-count-observable sides and only feeds the QC
scalar. Acceptable; called out for the record.

## Key decisions to confirm

- **D1 — Retire the combine** (regions + sides + `_deconv_per_node`). The cleaning is the single blend.
  *Verified redundant* (above). **This is the decision to confirm.**
- **D2 — `w = I/(I+I₀)` is the only weight**, per-node, inside `cleaned_gdna_count`. The scalar `(2κ−1)²`
  is gone. (Agreed.)
- **D3 — `I₀` becomes a live `CalibrationConfig` scalar** (`gdna_strand_info_scale`, default 1). (Agreed.)
- **D4 — Order:** `strand_deconvolve → clean → node_gdna_density(cleaned) → splice-junction upgrade →
  region masses`. The cleaning is upstream of the density/imputation, so AMBIG inherits clean; the splice
  upgrade stays downstream (see D5).

## Open decisions (flagged — these need a call)

- **O1 — Splice-junction upgrade composition.** `region_splice_gdna_frac` currently replaces the
  absolute-density fraction with the boundary gDNA-fraction for eligible exon regions (a *capture*
  correction, orthogonal to nascent cleaning). Proposed: keep it, applied to the **cleaned** region
  fraction (clean → density → splice-upgrade). Risk: it reads boundary crossing-spliced as a mature
  reference and may interact with the cleaning. **Plan: keep as-is on the cleaned fraction, and check on
  the benchmark suite that it still helps; reconsider in 3b if it conflicts.**
- **O2 — FP-rate quantile (`gdna_deconv_quantile`).** The strand-side quantile is retained (it flows
  through `split.gdna_frac = g_q`, set in `strand_posterior_gdna_frac`). The **count-side** widening
  (`center + Φ⁻¹(q)·σ_count`, today in `_deconv_per_node`) is lost when the combine is retired. Proposed:
  re-home it into `node_gdna_density` (which already computes `count_gdna_frac_var`) so the knob keeps
  working on count-routed nodes. **At the default `q=0.5` this is a no-op** (`Φ⁻¹(½)=0`), so 3a can ship
  without it and re-home in 3b if a non-½ quantile is needed. (The count σ is anti-calibrated under
  capture — it only ever *widens* — so dropping it is conservative, not a correctness loss.)

## Files & changes

- `config.py`: add `gdna_strand_info_scale: float = 1.0` to `CalibrationConfig` (+ `__post_init__` `>0`).
- `calibrate.py`: replace the `deconv_regions`/`deconv_sides` block with: `strand_deconvolve` → clean the
  three views → `node_gdna_density(gdna_counts=cleaned)` → region masses from `count_gdna_frac` → side
  masses from the cleaning fraction → `derive`. Keep `region_splice_gdna_frac` (O1).
- `strand_deconv.py`: **3b** — delete `_deconv_per_node`, `deconv_regions`, `deconv_sides`, `NodeDeconv`,
  `_compute_side`/`_SideQuantities`/`_left_right_neighbors` if now unused. `strand_deconvolve` +
  `cleaned_gdna_count` remain.
- `density_model.py`: **3b/O2** — optionally apply the FP-quantile to `count_gdna_frac` (default no-op).
- Golden: **regenerate** (`--update-golden`) — the first intentional output change.

## Staging (one output change, then dead-code removal)

- **Phase 3a — the architecture change (one golden regen).** Wire the cleaned flow into `calibrate`,
  retire the combine *path* (stop calling `deconv_regions`/`deconv_sides`), region + side masses from the
  cleaning. Regenerate golden. **Validate:** the three behaviours; the toy AMBIG/κ-sweep end-to-end (now
  through `calibrate`, not just the count field); the **gDNA benchmark suite** (`calibration-benchmark` —
  net leak/siphon vs HEAD, expecting AMBIG/nascent better, no regression elsewhere).
- **Phase 3b — cleanup (golden-identical).** Delete the now-dead combine code; re-home the FP-quantile
  (O2) if wanted; resolve O1; update `calibrate.py` docstring + `docs/calibration/calibration_theory.md`
  + CLAUDE.md (the blend → the single-mechanism cleaning).

## Test / validation plan

- Unit: `calibrate` region masses == `node_gdna_density(cleaned).count_gdna_frac · mass` (the retirement
  is exact); side masses == cleaning fraction · side mass; unstranded ⇒ raw floor (combine-free).
- Golden: regenerated; diffs explainable (AMBIG/exon gDNA fractions drop where stranded + nascent; pure
  count-only conditions unchanged).
- **End-to-end toy:** the Phase-2 numbers now hold through full `calibrate` (`phase3_ambig_scenario.py`):
  AMBIG `est` collapses at gDNA=0+nascent, stranded; degrades to the floor as κ→½.
- **Benchmark suite:** `net_flow_per_condition.tsv` vs HEAD — gDNA→RNA leak/siphon, per pool; the
  acceptance gate.

## Exit

Combine retired; the cleaning (`w=I/(I+I₀)`) is the only blend; golden regenerated with explainable
diffs; the three behaviours validated end-to-end; the benchmark suite shows the AMBIG/nascent improvement
with no pool regressing. 3b leaves the tree free of the dead combine code.

## Implementation corrections (as-built, 3a)

Three things the dry-run got wrong or deferred, corrected while implementing:

1. **The sides feed the prior — `deconv_sides` is KEPT, retired for *regions* only.** The plan assumed
   the side masses were QC-only (`gdna_density_global`). In fact `priors._transport_boundary_flux` builds
   the per-locus gDNA prior from `contained + left + right`, so the side masses are prior-critical. The
   plan's symmetric retirement (floor-1 side cleaning) injected false gDNA at AMBIG boundaries → the
   transport carried it into the region → `antisense_contained` gdna_em **4→105** (a regression). Fix:
   keep `deconv_sides` — a side's count fraction is the raw own/borrowed crossing density (not a
   cleaned-count field), so its blend `w·g_strand_side + (1−w)·g_count_side` is **not** redundant — fed
   the *cleaned* region density. The combine is retired for **regions only**; `_deconv_per_node` /
   `deconv_sides` stay (3b may unify the side weight to `I/(I+I₀)`). Result: antisense 4→**0** (truth 0).

2. **O1 (splice composition) resolved in 3a via the 3-term form, not deferred.** The toy scenario showed
   it conflicts: `region_splice_gdna_frac` used the **2-term** form (crossing-unspliced lumped as gDNA),
   so a stranded exon with nascent read its crossing as gDNA (single-strand exon `est 0.228` at
   gDNA=0+nascent) — and with the combine retired there was no strand step left to correct it. Fix: feed
   it the strand-**cleaned** gDNA crossing counts (`cleaned_left/right`) — the **3-term** form it already
   documented as "a later layer"; nascent moves to the RNA side. Result: single-strand exon `0.228→0.006`,
   AMBIG `est 0.004` at gDNA=0+nascent. This is the intended 3-term/strand-resolved splice layer, enabled
   by Phase 3's cleaning.

3. **Golden determinism: pinned the EM to 1 thread in the golden harness.** The locus EM's OpenMP
   reduction order is non-deterministic, amplified by the iterative solver to ~1e-8 relative on the
   largest scenario (`higher_frag_count`) — past the bit-exact tolerance. Pre-existing; surfaced when the
   regen landed in a different convergence basin. Fix: `_pipeline_config` pins `EMConfig(n_threads=1)` so
   the regression test is deterministic (it tests the algorithm's output, not the parallel reduction).
