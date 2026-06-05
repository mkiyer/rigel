# Boundary-Flux Transport — Implementation Plan

Length-bias-free re-attribution of boundary-crossing gDNA mass, to fix the capture
gDNA→RNA leak. Builds on the theory note (`phase6_exposure_efflen_theory.md`) and the
third-party review (Approach 2). Python-only; no C++ accumulator change.

## 1. Problem & concept

Under hybrid capture, the fractional accumulator splits a boundary-crossing fragment
across regions by **overlap length** (the uniform-sampling assumption). At an exon→intron
seam this deposits "smear" mass onto the *unexposed* intron (confirmed: intron contained=0
but boundary-side mass ≠ 0). Over the intron's long span that smear inflates the gDNA
effective-support length (IPR), dilutes the gDNA per-position density, and lets gDNA leak
into RNA.

**Concept:** a boundary fragment should be attributed to each side in proportion to that
side's *sampling probability* — `exposure × geometric capacity` — not overlap alone. The
exposed (probed) side, which actually generated the fragment, owns it; the unexposed side
gets ≈0. Exposure is the length-bias-free density (`mass/length`), and capacity is the
directional boundary effective length `𝓔(L)`.

## 2. Theory

### 2.1 Directional boundary effective length 𝓔(L)
A fragment of length ℓ can originate in a region and cross its boundary only if it starts
within ℓ of the seam. Integrating the FL survival function `F(x)=P(ℓ≥x)`:

```
𝓔(L) = ∫₀ᴸ F(x) dx = E[min(ℓ, L)]
```

- `L ≫ μ_FL` → `𝓔 → μ_FL` (a long intron contributes only an FL-window at the seam).
- `L → 0`   → `𝓔 ≈ L` (linear).

**This is exactly our existing `effective_length.boundary_side_eff_length`** (verified:
`E[min(ℓ,L)] = Σ_{x<L} F(x)`). No new geometry needed.

### 2.2 Transport rule
For an internal boundary between regions L (left) and R (right) with pooled gDNA mass
`M = right_g[L] + left_g[R]`:

```
P(origin = L) = ω_L·𝓔(L_L) / (ω_L·𝓔(L_L) + ω_R·𝓔(L_R))      (R analogously)
ΔM_L = M · P(origin=L)     ΔM_R = M · P(origin=R)
```

Per region: `g_region = contained_g + Σ_boundaries (transported share)`. Exposure
`ω_r = (g_region/L_r)/ρ₀` is **density**, which is length-bias-free
(`density = coverage/FL`, independent of L). Seed ω from the naive (pre-transport)
density; iterate (§4).

### 2.3 gDNA effective length (IPR) from transported mass
```
gdna_eff_len[locus] = (Σ_r φ·g_r)² / Σ_r φ·(g_r²/L_r)        L_r = region physical length
```
With the intron's smear transported away, its `g_r→0` and it drops out → the support
contracts to the exon span. The gDNA prior `alpha_gdna = Σ_r φ·g_r` (count) is unchanged by
transport (it only redistributes within the locus; the total is conserved).

## 3. The over-correction (measured)

The `ω·𝓔` rule is not exact for **sub-fragment regions**, because a region shorter than a
fragment is *encompassed* by pass-through fragments (started upstream, transiting it). Those
cross *both* its boundaries, so a high-ω small region is credited at both seams with more
than its overlap — partly undoing the accumulator's encompass-split.

| small exon | ω_small/ω_large | locus IPR vs true | direction |
|---|---|---|---|
| 50 bp  | 4.2× | 1221 / 1550 (−21%) | under → gDNA siphons |
| 150 bp | 4.3× | 1069 / 1650 (−35%) | under → gDNA siphons |
| 500 bp | 2.9× | 1524 / 2000 (−24%) | under → gDNA siphons |

Key properties:
- **Bounded**, not unbounded: `𝓔(L)` caps the small region's weight, so the over-correction
  saturates at ~3–4×, not FL/L→∞.
- **Always in the safe direction**: IPR *under* true → higher gDNA density → gDNA *wins*
  ambiguous reads (the ≥1.0, FP-averse behavior we want). It costs a little small-exon RNA
  sensitivity, never gDNA false-negatives.
- **A crossover at the 1 bp limit**: `𝓔(1)≈1` collapses the tiny region's weight, so it
  *under*-grabs instead. A pure `ω·𝓔` rule is therefore approximate at both extremes; exact
  behavior needs the per-position field (§7).

**Decision: accept the bounded, safe over-correction for v1** (it meets the ≥1.0 bar
conservatively for the common case) and schedule the field as the exact upgrade. Optional
v1.5 mitigation: cap each region's grab at its origination capacity `ρ₀·ω·𝓔(L)` and route
the excess (pass-throughs) onward — reduces over-siphon at the cost of complexity.

## 4. Convergence — no tolerance knob

Measured per-iteration total mass moved (50 bp exon, extreme capture):
`iter2=313.7, iter3=15.7, iter4=0.86, iter5=0.05, iter6=0.003 → 0`.

It crosses **<1 fragment-equivalent (sub-count) by iteration 4**. Stop rule = the natural
quantum **"total moved < 1 fragment"** (no observable change), backstopped at 8 iterations.
This is a physical unit, not a tuned tolerance, so it honors the no-magic-numbers rule.

## 5. Implementation plan (Python-only)

1. **`effective_length.py`** — `boundary_side_eff_length` already computes `𝓔(L)`; add a
   one-line docstring noting its dual role as the directional boundary effective length.
2. **`result.py`** — add `gdna_side_len: np.ndarray` (`𝓔(L)` per region; pure geometry) to
   `CalibrationResult`; validate finite ≥0 in `__post_init__`.
3. **`calibrate.py`** — pass `gdna_side_len = bside_eff` (already computed).
4. **`priors.py`** —
   - Add `_transport_boundary_flux(contained_g, left_g, right_g, length, E_cap, ref_id)`:
     fully **vectorized** (boundary scatter-add, no Python loop over regions), iterating
     until `Σ|Δg| < 1.0` (cap 8). Conserves total mass (untransported ref-edge sides kept).
   - `assemble_priors`: run the transport on the per-region gDNA → `g_tr`; project `g_tr`
     and `g_tr²/L` to loci; `alpha_gdna = w·Σφ·g_tr`, `gdna_eff_len = max(1, (Σφ·g_tr)²/Σφ·(g_tr²/L))`.
     Multimap-blind loci have `g≈0` → IPR falls back to the geometric span (no over-attraction).
5. **Tests** —
   - `test_result_schema.py`: add `gdna_side_len` to the valid kwargs.
   - `test_priors.py`: rewrite for the transport (intron smear → exon; intron g→0).
   - New `test_boundary_flux.py`: unit-test `_transport_boundary_flux` — (a) uniform ω
     reduces to capacity-weighted overlap; (b) one exposed + one unexposed side → exposed
     wins; (c) mass conservation; (d) convergence < 8 iters.

## 6. Validation plan

- **Realistic 1 kb-exon capture scenario** (`dissect_capture_gdna_leak`), binding sweep
  1→1000: gDNA recovery → ~1.0 (≥1.0 acceptable); RNA over-recovery → ~1.0.
- **Small+large-exon scenario** (`proto_boundary_flux_transport`): intron exposure → 0;
  support → exon span (within the documented −20–35% safe band).
- **gDNA-level sweep** (0.05 / 0.10 / 0.45): leak small at realistic contamination.
- **Multimap regressions** (4 tests) + full calibration unit suite: must stay green (the
  transport is a no-op where g≈0).
- **Goldens**: regenerate + inspect.

## 7. Limitations & the per-position field (v2)

The transport is a per-region, adjacent-2-region approximation. It cannot model
pass-throughs spanning >2 regions exactly (hence the sub-fragment over-correction and the
1 bp crossover). The exact, general solution is the **per-position capture-bias field**
`ω(x)` (Approach 1): deconvolve gDNA coverage at fine (10–50 bp) resolution, set the
component effective length to the field IPR `(∫ω)²/∫ω²`, and attribute every fragment ∝
`ω(x)` over its footprint — no regions, no ownership, no length bias at any scale. It costs
a per-locus coverage array (modest) and a deconvolution step. Recommended as v2 once the
transport (v1) is shipped and proven to resolve the leak for typical exons.
```
