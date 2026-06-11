# FL-consistency of the splice-junction partition — diagnostic findings

**Status:** diagnostic complete, 2026-06-11. Probe: `scripts/debug/diag_fl_consistency.py`
(`/tmp/fl_consistency.png`). Resolves the deferred `count_mean_bias_design.md` §9.1 open question
("understand before fixing"). **Verdict: real, large, present in the benchmark — fix recommended.**

## The bug (precise)

`splice_junction.region_splice_gdna_frac` sets an eligible exon region's gDNA *count* fraction to the
boundary fraction `f_b` **directly**. But `f_b` is a **molecular density** fraction: the boundary
divides crossing counts by the **crossing** eff-length `fl_mean`, which cancels to `d_g/(d_g+d_r)`. It
is then applied to the region's contained **count** `M_region`. The correct count fraction of the
contained mass is

```
g_true = d_g·E^g_region(L) / ( d_g·E^g_region(L) + d_r·E^r_region(L) )
       = (f_b·E^g_region) / ( f_b·E^g_region + (1−f_b)·E^r_region )      # exact identity
```

with `E^{g,r}_region(L) = region_eff_length(L, fl_pmf_{g,r})`. So `frac = f_b` is correct **iff**
`E^g_region(L) = E^r_region(L)`, i.e. iff the gDNA and RNA fragment-length distributions are identical.
The fix is exact and **parameter-free** — it reuses `region_eff_length` with the gDNA and RNA pmfs,
both already computed in `calibrate`.

## Why it was missed, and why it is NOT hypothetical

The first cut used the bare form because the toy scenario recovered truth — that scenario had **matched
gDNA/RNA FL** (`r=1`, zero bias; confirmed by the control row). But the **benchmark suites set a real
gap**: RNA `frag_mean=250`, gDNA `frag_mean=350` (`gdna_benchmark_5mb`, `quick_3to1_5mb`,
`nascent_benchmark_5mb`); and `gdna_shortfl_5mb` sets gDNA `frag_mean=100`. So the bias is live in every
suite, and extreme in shortfl.

## Quantified bias (exact benchmark FL params, no sampling noise)

Signed bias `g_bare − g_corr` of the gDNA *fraction* at `f_b=0.5` (`+` = bare over-calls gDNA → RNA
siphon; `−` = bare under-calls gDNA → FP RNA leak):

| exon L (bp) | control (matched) | **benchmark std** gDNA 350 vs RNA 250 | gdna_shortfl gDNA 100 |
|---|---|---|---|
| 150 | 0.000 | +0.036 | **−0.492** |
| 200 | 0.000 | **+0.159** | **−0.460** |
| 250 | 0.000 | **+0.233** | **−0.383** |
| 300 | 0.000 | **+0.246** | −0.287 |
| 400 | 0.000 | +0.188 | −0.166 |
| 600 | 0.000 | +0.085 | −0.088 |
| 1000 | 0.000 | +0.036 | −0.045 |

The bias peaks at exon-scale lengths (200–400 bp, squarely in the real exon-length range) and decays to
<0.05 only beyond ~600 bp. **In the standard benchmark the splice upgrade over-calls gDNA by 16–25
percentage points at typical exons**; in shortfl the bare form under-calls by up to ~0.5. (Rows below
~120 bp carry near-zero contained mass — both eff-lengths → 0 — so the partition is moot there.)

## Two structural insights that shape the fix

1. **The absolute-density path is ALREADY FL-consistent; only the splice path is not.** `density_model`
   computes `count_gdna_frac = clip(density·E^g_region / M_region) = d_g·E^g_region/M_region = g_true`
   — it already multiplies the gDNA density by the gDNA region eff-length. So the fix simply makes the
   splice fraction match the consistency the density estimate already has (the same arithmetic).
2. **The Phase-3 density floor protects only ONE direction.** With the agreed `g = max(g_density_floor,
   g_splice)`: in the gDNA-shorter (shortfl) case the bare `g_splice` *under*-calls and the floor (the
   correct, higher density estimate) wins → **protected**; in the gDNA-longer (benchmark std) case the
   bare `g_splice` *over*-calls and `max` picks it over the correct lower floor → **not protected**. So
   the floor does not substitute for the fix in the production-representative (gDNA-longer) direction.

## Where it propagates (mediation — for the validation, not the diagnostic)

The biased fraction is the calibration **prior** (`count_gdna_frac` → `deconv_regions`). Its effect on
the output is mediated by the blend weight `w=(2κ−1)²`: on strand-observable exon regions in a
**stranded** library the count fraction gets only `~(1−w)≈4%` weight (washed out), but on **count-routed
(unstranded / AMBIG)** nodes it gets full weight. So the FL-consistency bias bites hardest in unstranded
data — exactly where there is no strand signal to correct it. Quantifying the output-level effect is a
pipeline measurement (validation), not part of this analytical probe.

## Recommendation

**Implement the fix** — it is exact, parameter-free, low-risk, and corrects a 16–25-point fraction bias
at typical-exon lengths that is present in the benchmark today:

- In `splice_junction` (`boundary_gdna_fraction` or `region_splice_gdna_frac`), convert the boundary
  density fraction to the region count fraction via the region eff-length ratio:
  `g = (f_b·E^g_region)/(f_b·E^g_region + (1−f_b)·E^r_region)`, using `region_eff_length(L, gdna_pmf)`
  and `region_eff_length(L, rna_pmf)` (both already available in `calibrate`).
- It composes cleanly with the Phase-3 3-term (same arithmetic; the 3-term changes the numerator's
  gDNA/RNA split, the eff-length ratio converts density→count). Do it **as part of, or just before,**
  the Phase-3 splice work.
- **Validation:** `gdna_shortfl_5mb` is the ready-made stress test (gDNA 100 vs RNA 250 → up to −0.5
  bias) — it should show a large gDNA→RNA leak at short exons today that the fix removes; confirm the
  standard suites improve (or stay flat where the blend washes it out) with no regression. Report the
  net-flow metric restricted to short multi-exon regions, where the effect concentrates.

## Open follow-ups

- **Stage (b) realized scenario:** construct short-exon transcripts with known gDNA+mature mixtures at a
  controlled FL gap and confirm the corrected partition recovers the planted gDNA count where the bare
  one does not (extends `test_splice_junction_realized.py`). Do this alongside the fix as its test.
- The probe assumes uniform within-region coverage and that the region's contained-RNA FL equals the
  spliced reference FL (true for nascent + mature-within-exon, same library prep). Capture can distort
  within-region coverage; the realized scenario + suite validation cover that.
