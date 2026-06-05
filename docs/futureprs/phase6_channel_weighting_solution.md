# Phase 6 — Channel-Weighting Solution (count × strand fusion)

The joint decode infers the per-node gDNA fraction π_g by combining two **orthogonal**
evidence channels — count and strand. The Bayesian-correct fusion is the **product** of the
two likelihoods, with each channel's influence equal to its **precision** (its reliability).
The bug: the count channel's precision was mis-set, so it overrode the strand channel exactly
where the strand channel was right (capture-enriched exons). This is the fix.

## 1. The two channels and their precision

**Strand channel** — a Beta-Binomial likelihood on (sense, antisense). Its discriminating power
between π_g=0 and π_g=1 scales with `|κ_rna − 0.5| · √N`: the separation of the sense rates
(`κ_rna` for RNA vs `0.5` for gDNA) times the Binomial sharpness. So it is **already
self-calibrating** — automatically *flat* (uninformative) for unstranded data (κ_rna=0.5) and
*sharp* for strand-specific data. **No change needed.**

**Count channel** — a Beta prior on π_g with mean `π_count = ρ·L_eff / M` (the expected gDNA
fraction) and concentration (effective sample size) `κ_c`. The bug was **κ_c = the raw read
count M**. That is the concentration you'd get only if you had *observed M labelled fragments*
and counted which were gDNA. You haven't: the calibrator only *models* the gDNA rate ρ. Most of
the M reads at an exon are RNA and carry **no information** about the gDNA density. Setting
κ_c = M therefore claims far more certainty than the channel has, and a sharp prior on the
(imputed, capture-biased) mean overrides the correct strand likelihood.

## 2. The fix — κ_c = the expected gDNA count μ_g (principled, self-gating)

The count channel knows the gDNA count to **Poisson precision**: the expected gDNA is
`μ_g = ρ·L_eff`, known to ±√μ_g, i.e. an **effective sample size of μ_g**. So:

```
κ_c = μ_g = ρ · L_eff            (the expected gDNA count, NOT the total read count)
```

This is **self-gating**, with no branch and no constant:
- **count-decodable region** (intergenic/intronic): ρ is measured from the region's own count,
  ρ = M/L_eff, so `κ_c = ρ·L_eff = M` — full confidence, unchanged from today (correct there).
- **imputed exon**: ρ is the small swept density, so `κ_c = ρ·L_eff ≪ M` — low confidence → the
  strand channel drives. Under capture, ρ is biased low (depleted neighbours) but now *weak*, so
  the strand clue (which correctly reads the antisense excess) wins.

It is principled (the channel's confidence = the gDNA evidence it actually has, to Poisson
precision), data-driven (μ_g from the density × geometric length), needs **no magic number**, and
it **simplifies** the code: `κ_c = density·L_eff` replaces the separately-tracked discrete-count
field, which can be removed. The earlier "gated" variant is identical (count where decodable,
μ_g elsewhere) — μ_g *is* the gate, as one expression.

## 3. Evidence

Decode level (stranded capture, captured exons): gDNA fraction detected **π_g 0.026 → 0.287**
(matching the strand-implied ~0.29); total detected gDNA ×15.

End-to-end, stranded capture, across contamination (the metric that matters):

| gDNA fraction | gDNA recovery (count → μ_g) | RNA over-recovery (count → μ_g) |
|---|---|---|
| 0.05 (realistic) | 0.62 → **0.73** | 1.019 → **1.013** |
| 0.10 (realistic) | 0.58 → **0.68** | 1.042 → **1.032** |
| 0.45 (extreme stress) | 0.44 → **0.55** | 1.248 → **1.202** |

No regression elsewhere: stranded no-capture 0.85 → 0.86; unstranded no-capture improved
(RNA-over 1.097 → 1.088); 80 calibration-unit + multimap-regression tests pass.

**At realistic gDNA (5–10%) the leakage is small (RNA over 1.3–3%)** and μ_g consistently
improves it. The strand channel now drives where it's reliable, the count channel where it's
reliable — a rational, reliability-weighted fusion.

## 4. Honest residual

At **extreme** 45% gDNA capture, μ_g improves but does not fully eliminate leakage (recovery
0.55; RNA over 1.20). The decode now detects the gDNA, but the per-locus EM does not fully
convert the improved prior into gDNA assignment (the silent-probe leak is roughly unchanged),
pointing to a **separate EM-side competition** factor at extreme contamination. This is far
beyond realistic gDNA levels and is out of scope for the channel-weighting fix; it is a
candidate for separate work if extreme-contamination capture must be perfect.

## 5. The complete Phase 6 calibration fix

Two principled, independent changes make the channel combination rational across all library
types:

1. **Option A** — gDNA effective length is the **geometric** span (ω-free); exposure lives in
   the mass, giving gDNA its true *local* density at enriched exons. (Fixes the multimap
   eff-len collapse; correct capture semantics.)
2. **μ_g count weight** — the count channel's strength is the **expected gDNA count**, so it
   defers to the strand channel where RNA dominates. (Fixes the count-prior override; rational
   reliability-weighted fusion.)

Together: correct on strand-specific no-capture and strand-specific capture (the required
types), no regression on unstranded, multimap regressions fixed — with no magic numbers and a
net *simplification* of the decode.
