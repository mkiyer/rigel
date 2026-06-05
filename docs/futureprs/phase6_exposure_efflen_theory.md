# Calibration theory note — exposure, effective length, and the gDNA component

The canonical framework the calibration is built on, after the convergence of two views
("exposure as an effective-length factor" and "effective-support length"). It is the basis for
fixing the gDNA→RNA leak under hybrid capture.

## 1. Exposure and effective length are the same kind of quantity

Both answer: **"how probable is it to observe fragment mass at this position?"**

- **Effective length** is that probability *integrated over a region*, under the assumption of
  **uniform, unbiased** molecular sampling. A region of physical length `L` has effective length
  `L` (FL-corrected); a random fragment is observed there in proportion to `L`.
- **Exposure** `ω` is the empirical factor by which the *real* sampling deviates from uniform.
  Hybrid capture makes sampling non-uniform (probes enrich some positions, deplete others);
  exposure measures that non-uniformity from the data.

So the **exposure-corrected effective length** of a region is

```
L_eff_observed(r) = L_r / ω_r        ω>1 (enriched) → shorter; ω<1 (depleted) → longer; ω=0 → ∞
```

The EM uses effective length as a *density normalizer* (a component's per-position rate is
`abundance / eff_len`), so enrichment (high ω) ⇒ shorter effective length ⇒ higher density ⇒ more
competitive — which is correct.

## 2. What exposure IS, and how to measure it

Exposure is the **capture sampling bias** `p_capture(x)` — a property of genomic position,
applying to *every* molecule there (gDNA and RNA alike). To measure a sampling bias you need a
ruler of *known-uniform* abundance. **gDNA is that ruler**: it tiles the genome at a uniform
density `ρ₀`, so

```
ω_r = (gDNA density in r) / ρ₀ = (g_r / L_r) / ρ₀
```

This is why calibration deconvolves gDNA first. Total fragment density is *not* a valid ruler —
at an expressed exon it conflates capture bias with RNA expression.

## 3. The gDNA component effective length = the IPR of the exposure field

The gDNA component spans the whole locus, but under capture its mass is concentrated where
exposure is high. The single effective length the EM needs is the one that makes the component's
per-position rate equal the gDNA density *where its reads actually are*. That is the **inverse
participation ratio** of the exposure field:

```
eff_len_gDNA = (Σ_r ω_r L_r)² / Σ_r ω_r² L_r          [exposure form]
             = (Σ_r g_r)²     / Σ_r g_r²/L_r           [mass form, g_r = ρ₀ ω_r L_r — identical]
```

Properties:
- **Zero-exposure regions drop out** (ω=0 contributes 0 to both sums) — no ∞, no special-casing.
  This is the answer to "how do we aggregate L/ω across regions": *not* a sum (the ∞ from
  unexposed regions would swamp it) — the IPR.
- **Uniform contamination** (ω≡1): reduces to the geometric span `Σ L_r`.
- **Capture** (ω concentrated on exons): contracts to the exon support.
- Worked example: 98 kb at ω=0 + 2 kb at ω=50 → `(50·2000)²/(50²·2000) = 2000 bp`. Abundance
  `ρ₀·Σω L = 0.1·100000 = 10000` (the real gDNA count); density `10000/2000 = 5` = the true
  enriched exon density. Fair competition.

The component-level abundance (`G = ρ₀ Σω_r L_r`, the gDNA count) and the IPR effective length are
parameterization-independent: "exposure on the effective length" and "effective-support of the
mass" are the *same* object once aggregated, because `G ≡ ρ₀·span` by definition of `ρ₀`.

## 4. Boundary-crossing fragments — the open problem

The IPR is only as good as the exposure field `ω_r` it is built from, and that field is currently
corrupted at region boundaries.

**The bug (confirmed empirically).** A fragment spanning an exon→intron junction is split across
the two regions in proportion to **overlap length** (the fractional accumulator). Under uniform
sampling that is correct (origin is uniform over the footprint). Under capture it is wrong: the
fragment was captured *because it overlaps the exon probe* — its presence says nothing about
intron exposure. Proportional splitting deposits "smear" mass onto the zero-exposure intron:

```
extreme capture, intron region: contained = 0   but boundary-side mass = 140 (pure smear)
```

That smear gives the intron false exposure → inflates the gDNA IPR (the smear adds to `Σg` but,
over a long region, almost nothing to `Σg²/L`) → over-states the support → dilutes the gDNA
density → gDNA loses fair fights → **leaks into RNA**.

**The principle.** A boundary fragment's mass should be attributed to each side in proportion to
that side's *sampling probability*, i.e. **overlap × exposure**, not overlap alone:

```
share(r) ∝ overlap_r · ω_r          uniform: → overlap (current);  capture: exposed side wins
```

The exposed (probed) side owns the captured fragment; the unexposed side gets ≈0.

**The chicken-and-egg.** Exposure-weighted attribution needs exposures, which need the
attribution. Two resolutions:
1. **Contained-anchored (non-iterative).** Estimate ω from *contained* fragments only
   (unambiguous), then attribute boundary mass by those exposures. Works wherever regions are
   large enough to hold contained fragments (the common case for the gDNA support).
2. **Joint inference / per-position field (general).** For **tiny regions** (shorter than a
   fragment, so *all* their mass is boundary mass — common for small exons) contained estimation
   is impossible; exposure must be inferred jointly with the attribution (an EM over boundary
   fates) or as a per-position capture-bias field deconvolved from the gDNA coverage. This is the
   genuinely hard case and is deferred until the contained-anchored fix is proven.

## 5. Summary of the model

- Exposure `ω_r` = capture bias = gDNA density / ρ₀, measured from the gDNA ruler.
- gDNA component effective length = IPR of the exposure field = `(Σω L)²/Σω²L`.
- gDNA prior (count) = `G = Σ g_r`.
- Per-position gDNA rate = `G / eff_len` = the mass-weighted local density → fair EM competition.
- Boundary mass must be attributed by overlap × exposure (open work item §4).
