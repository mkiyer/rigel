# v4.3 VCaP Regression — Root Cause Analysis

**Date**: 2026-05-19
**Run**: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/v4_3_with_mm`
**Diagnostic**: `scripts/debug/diagnose_v4_3_regression.py`

## TL;DR

The regression is caused by **three interacting bugs** in the v4.3
regional-exposure field, all in
`src/rigel/calibration/_regional_exposure.py::RegionalGdnaExposure.build`.
They conspire to produce an `A(x)` field that suppresses nRNA and
mRNA denominators **harder than** gDNA denominators, which inverts
the EM rate ordering and pushes true-gDNA fragments into nRNA and
mRNA components.

The data is **VCaP shotgun gDNA + RNA, no hybrid capture**, so the
correct `A(x)` is uniform = 1 and auto-uniform should fire. It
doesn't, because the spread test sees artifactual spread from the
bugs below.

## Smoking gun (from quant outputs)

Median EM exposure weight per component:

| Component | Median weight | Effective L_em shrinkage |
|---|---:|---:|
| gDNA loci  | 0.084 | 12× |
| nRNA spans | 0.030 | 33× |
| mRNA exons | 0.021 | 50× |

EM rate per component is `count / L_em`. Smaller `L_em` ⇒ component
looks more enriched per bp ⇒ wins ambiguous fragments. nRNA wins
2.5× more aggressively than gDNA; mRNA wins 4× more aggressively
than gDNA. Confusion matrix matches: `gDNA → nRNA` 13.05%,
`gDNA → mRNA` 7.39%.

## Root cause 1 — `kappa = 1`, but `E ~ 5000–40000`

```text
rho_hat_r = (Y_r + kappa * rho_global) / (E_r + kappa)
```

`estimate_kappa` is returning `kappa_global = 1.0`. Typical
per-region exposure `E_r` is `~5000` bp (introns), `~40000` bp
(intergenic), `~565` (exon-boundary). With `kappa = 1`, the prior
contribution is **negligible** for every contained-fragment region.
For an empty region (`Y_r = 0`):

```text
rho_hat ≈ kappa * rho_global / E    ⇒    ~4e-8 (intron) to ~4e-7
A_r    ≈ rho_hat / rho_ref          ⇒    ~1e-4 (floor)
```

The shrinkage we expected — *"no evidence for this region, so fall
back to the global mean"* — never happens. Empty regions collapse to
zero `A`. Counterfactual with `kappa = mean(E)`:

```text
rho_hat(Y=0) ≈ kappa * rho_global / 2  ≈  rho_global/2
A ≈ 4.4 → clipped to 1
```

Empty regions correctly fall back to uniform.

**Why kappa is so small**: the method-of-moments estimator picks a
kappa that fits the *overall* `(Y, E)` variance. With the per-class
scale mismatch (RC2 below), inter-region variance is dominated by
inter-class differences, not by within-class Poisson noise. MoM
attributes the variance to "real" `rho_r` differences and shrinks
the prior precision toward zero. Garbage in (mixed scales), garbage
out (negligible shrinkage).

## Root cause 2 — per-class `ρ` scales differ by 4 orders of magnitude

Diagnostic output:

```text
class         n_reg      Y_tot      E_tot   rho_avg    rho_q50    rho_q95  q95/ref
INTERGENIC    33120   1.04e+04   1.31e+09  7.93e-06   1.45e-08   2.54e-07     0.00
INTRON       296039   3.01e+05   1.54e+09  1.95e-04   1.18e-07   3.80e-04     1.78
EXON-INTRON  329057   5.41e+06   1.81e+08  2.99e-02   3.08e-06   1.93e-01   904.14
```

The EXON channel has **150× higher mean `ρ`** than INTRON and
**~4000× higher** than INTERGENIC. Q95/`rho_ref` = 904 for EXON,
meaning the EXON channel alone defines `rho_ref` (every other class
sits below it).

**Why**: Pass-0 fragment categorization is exon-biased. Of 31.5 M
gDNA-eligible fragments:

```text
n_exon_only           18,713,690    (59%)
n_exon_intron         11,259,348    (36%)
n_intron_only            919,752    ( 3%)
n_intergenic_only        282,274    ( 1%)
n_exon_intergenic        275,922    ( 1%)
n_intron_intergenic        2,226    ( <1%)
```

95% of fragments classify as exon-touching. Intronic and intergenic
*evidence* channels see only 1–3% of the gDNA mass, but their
*opportunity* `E` is calibrated for the full contained-fragment
density. Result: `Y/E` is artificially low in INTRON/INTERGENIC and
artificially concentrated in EXON.

The EXON channel `E = sides × b_cross` is a **different physical
unit**: bp-equivalent boundary-crossing exposure. The mathematical
units cancel (rho is "fragments per opportunity bp" in both cases),
but `b_cross` is a per-boundary scalar whose absolute calibration
isn't necessarily comparable to contained `L_eff`. Empirically the
EXON-channel `ρ` is ~150× the INTRON-channel `ρ`. v4.3 §1 anticipated
this risk and said "if it happens, fix `Y/E`, don't hide behind
class normalization." This is the "if it happens" moment.

## Root cause 3 — auto-uniform spread test fails on shotgun gDNA

```text
observed_log_spread = 5.82
null_log_spread     = 0.54
```

The spread is "real" but **not** biological — it comes from:

1. Most regions collapsing to `rho_hat ≈ kappa·rho_global/E` (RC1)
   → bottom of the distribution at `~1e-8`.
2. EXON channel pushed to `~0.2` (RC2)
   → top of the distribution at `~0.2`.

That's `log(0.2 / 1e-8) ≈ 17` nats of apparent spread, of which
none is biological exposure variation. The spread test fires "use
regional," but for VCaP shotgun gDNA the truth is uniform.

## Fix plan

Three small targeted changes, all in
`_regional_exposure.py::build`. None require new abstractions.

### Fix A — kappa floor

```python
kappa_data = float(estimate_kappa(Y[valid], E[valid], rho_global).value)
E_mean = float(E[valid].mean())
kappa = max(kappa_data, E_mean)         # NEW
```

Empty regions now shrink to `~rho_global / 2` ⇒ `A ≈ 0.5–1`.
Cost: regions with `Y > 0` and `E ≪ E_mean` get slightly more
shrunk than today. Acceptable; these are the regions whose `ρ̂` is
least trustworthy.

### Fix B — drop EXON from the global field

Compute `rho_global`, `rho_ref`, and the spread test from
`is_intergenic | is_intron` only. EXON regions still get a per-region
`ρ̂_r` (using the same kappa) and contribute to the `A` field
locally, but their values will mostly clamp to `A = 1` (above
`rho_ref`) and don't pull the scale up.

```python
field_mask = (is_intergenic | is_intron) & valid    # NEW
rho_global = Y[field_mask].sum() / E[field_mask].sum()
kappa_data = estimate_kappa(Y[field_mask], E[field_mask], rho_global)
# shrinkage and rho_ref computed on field_mask, not valid
```

Justification: INTRON + INTERGENIC are physically the same channel
(contained-fragment evidence). They are directly comparable. EXON
is a boundary-crossing channel whose absolute calibration isn't
guaranteed to be on the same scale.

### Fix C — spread test weighted by evidence

```python
# OLD: weights = E[valid]
# NEW: weights = Y[field_mask] + eps      (or Y + kappa)
```

Empty regions should not contribute to the observed spread. The
spread should measure *what we have evidence for*, weighted by *how
much evidence we have*.

### Combined effect on VCaP

With A+B+C:

- `rho_global` from INTRON + INTERGENIC only → ~2 × 10⁻⁴ (dominated
  by the few intron regions with evidence).
- Empty regions shrink to `~rho_global / 2`, well within an order
  of magnitude of `rho_ref`.
- Observed spread (Y-weighted) collapses to ~1 nat — under the null
  spread → auto-uniform fires.
- All components get `A = 1` ⇒ `L_em = L_raw` ⇒ EM behaves like
  `--regional-exposure off`.

For real hybrid-capture data, intron-rich nRNA spans correctly see
lower-A territory (because true gDNA is enriched only over the
capture target, not over their introns), while gDNA loci correctly
see higher-A territory.

## Recommended next steps

1. **Implement A + B + C** in a single PR. ~30 lines of code.
2. **Rerun VCaP** and confirm the confusion matrix returns to (or
   beats) the v3 uniform baseline.
3. **Rerun the synthetic 24-condition harness** to confirm no
   regression in uniform regimes.
4. **Find or simulate one capture dataset** (exome or panel) for the
   non-trivial regional case. The current regression evidence is
   entirely on shotgun data, where the correct answer is "do
   nothing." We have no positive test that regional mode helps when
   it *should* help.
5. **Consider gating regional mode behind an explicit capture flag**
   (e.g. `--capture-bed targets.bed`) until we have a positive
   capture benchmark. `auto` mode without capture data has so far
   only hurt.

## What we should NOT do

- Reintroduce numerator weighting. RC1 and RC2 are about the
  *field itself* being miscalibrated; weighting fragments by a bad
  field is worse than weighting denominators by it.
- Reintroduce per-class `rho_ref`. The problem is the EXON channel
  being on the wrong scale, not the existence of classes. Hiding
  the mismatch behind within-class normalization (v4.1) was a
  symptom-patch.
- Tune `LOG_A_FLOOR` or `LOG_RHO_CLIP_NATS`. These are guards, not
  calibration knobs.

## Open questions for the user

1. **Authorize implementation of A + B + C** in this turn? It is
   tightly scoped to `_regional_exposure.py::build` and the
   `weighted_quantile` / spread-test paths.
2. **Gate `auto` on a capture target file?** If we have no
   evidence that regional mode helps shotgun data, the safer
   default for `auto` may be "off unless capture target provided."
3. **Do we have any capture dataset** (exome, ribo-depleted with
   pull-down, etc.) to validate that regional mode is useful when
   the field is real?
