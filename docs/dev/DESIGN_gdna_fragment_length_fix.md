# Fixing the gDNA fragment-length model — options, with what is already REFUTED (2026-08-30)

Diagnosis: `DIAGNOSIS_gdna_fragment_length.md`. Three independent defects, and a fix for one is not a fix
for another:
**A** the four "pure gDNA" pools are contaminated by unannotated/nascent RNA (capture-OFF; `bias ≈
RNA_share × length gap`, up to −121 bp); **B** the opportunity divisor is blind to the probe panel
(capture-ON, +25 bp even at equal lengths); **C** the EB shrinkage anchor is the global MIXTURE at a magic
ESS of 1000 (capture-ON, once the gDNA pools collapse).

## 0. ⛔ ALREADY REFUTED BY MEASUREMENT — do not re-propose (`scratchpad/fl/unmix_feasibility.py`)

**Subtracting a scaled copy of the spliced pool.** The obvious unmixing —
`g(L) = pool(L) − (1−α)·N·r(L)` with `r` the certified-pure spliced pool — FAILS, even handed the ORACLE
`α`, because the template is wrong and the error is amplified by `1/α`:

| template | TV(contaminant, template) OFF | TV ON | equal-length g05 OFF: true 215.9, pooled 216.7 → |
|---|---|---|---|
| raw spliced pool | 0.148–0.187 | 0.15–0.23 | **144.9** (−71.0) — it BREAKS a correct estimate |
| de-tilted spliced pmf | **0.053–0.091** | 0.18–0.27 | 220.8 (+4.9) — tolerable |

⭐ De-tilting the template (dividing by the sj crossing probability) is what makes it nearly right OFF
capture: on the fl-gap arm it recovers 248.0 / 261.6 / 261.8 at g05/g25/g50 against truths 262.0 / 261.5 /
261.8, where the shipped model reads 140.9 / 200.1 / 233.6.
⛔ **But no template is valid UNDER capture** — the probes reshape the RNA length distribution too, the sj
de-tilt is equally probe-blind, and subtraction then makes things WORSE (equal-length g05 ON: pooled 264.1
against a truth of 261.4, unmixed **285.7**). ⛔ And at low gDNA the `1/α` amplification dominates: at
α = 0.18 a 5 % template error is still a −14 bp residual.

## Option A — SPLIT THE LENGTH BANKS BY STRAND ⭐ RECOMMENDED

The exact analogue of the away-half, and it needs **no template and no mixing weight**. For a pool whose
objects have a known gene strand, with `S`/`A` the sense/antisense length histograms:

    S(L) = ½·g(L)·N_g + κ·r(L)·N_r
    A(L) = ½·g(L)·N_g + (1−κ)·r(L)·N_r
    ⇒  r(L)·N_r = (A − S)/(1 − 2κ)        g(L)·N_g = 2S − 2κ·(A − S)/(1 − 2κ)

⭐ **At the real libraries' κ ≈ 0.002 this is `g ≈ 2S`: the gDNA length distribution is simply twice the
sense-strand histogram.** Conditioning is `1/(1−2κ) ≈ 1` — no amplification anywhere, unlike the `1/α` of
every template method. κ is already fitted, and its own estimator was made contamination-robust on
2026-08-30.

⭐⭐ **THREE OF THE FOUR POOLS ARE ORIENTABLE**, including both crossing pools that dominate under capture:
pool 1 (intronic — the gene's strand), pool 2 (intron|exon), pool 3 (intergenic|exon — the exon's gene
supplies the strand). ⛔ Only pool 0 (intergenic contained) has no strand to orient by, exactly as
intergenic seeds are excluded from the strand fit. Since gDNA's length distribution is a LIBRARY-WIDE
property of one fragmentation process, estimate `g(L)` from the orientable pools and use pool 0 for its
COUNT and as a consistency check, never for its shape.

**Cost:** `pool_lengths` becomes `[pool, strand, L]` — a C++ accumulator change, a payload schema bump and
a rebuild of every scan cache. Bounded and mechanical; the deposit rule itself does not change.
**Blind spot:** an unstranded library (κ = ½ ⇒ `1 − 2κ = 0`) carries no orientation and the estimator is
undefined. That is honest and is the same blind spot the strand fit has. Fall back to Option B off capture,
and to reporting the contamination under it.

## Option B — DE-TILTED TEMPLATE SUBTRACTION, capture-OFF only (the fallback)

As measured above: works off capture (essentially exact at g25/g50), needs **no schema change**, and works
on unstranded libraries where Option A cannot. ⛔ Must be gated OFF under capture, where it is harmful.
`α` is available without circularity from the strand deconvolution of the same objects' *counts*
(`substrate.region_contained.count` is already strand-split), so it needs no new machinery either.

## Option C — the two capture-ON defects, independent of contamination

**C1 · a probe-aware opportunity divisor.** `gdna_opportunity_from_index` is computed from the index alone,
so it removes only ~6 bp of the crossing pools' ~30 bp length selection under capture. `capture_eff_length`
already models the probe panel; the divisor should use it. ⭐ This is what the `fl_anchor_gap.py` `G-gdna`
control has been reporting as "impossible, and therefore a bug" (+6.0 % on all six capture-ON rows) since
2026-08-17.

**C2 · stop shrinking a component toward a mixture.** Replace the `POOL_EB_PRIOR_ESS = 1000` pull toward
`global_pmf` with either (i) within-component smoothing at a derived bandwidth, or ⭐ (ii) **the same
reconciliation the strand fit now uses**: pools 1/2/3 all estimate the SAME `g(L)` under different
opportunities, so a thin pool should borrow from the better-measured ones — its own components, never the
mixture. (ii) removes the constant rather than re-tuning it.

## Option D — ⛔ NOT RECOMMENDED: iterate against the composition solve

Weighting each object's length contribution by its solved `f_g` is circular: the fl model sets the
effective lengths the solve uses. `fl.py`'s purity assertion exists precisely to avoid this. If it is ever
revisited it must be as an explicit EM with a convergence proof, not a refit loop.

## Sequencing, and the one thing that must be measured FIRST

⛔ **The price of A is UNKNOWN and cannot be ranked yet.** `length_ceiling.py` prices a perfect gDNA pmf at
−5.90 %, but "all of it capture-ON, capture-OFF −0.0000" — that zero is the equal-length panel design
making contamination invisible, so the priced 5.90 % belongs to B and C. Contamination's own price has
never been measured on a substrate that can see it.

1. **Regenerate the fl-gap arms on the current sparse-nascent model** (they still carry the retired uniform
   one) and price A there with `length_ceiling.py`. Owner decision; nothing can be ranked without it.
2. **C1 and C2** — independent of A, cheap, and they own the already-priced 5.90 %.
3. **A** if step 1 justifies the schema change, with **B** as the unstranded/no-rebuild fallback.

⚠ Whatever is built, judge it on a panel where the two components' fragment lengths DIFFER. The equal-length
forcing function is deliberate and correct for the EM (`TRAPS: a-length-gap-bypasses-calibration`) — but it
makes this entire class of defect unobservable, and that is why a 95 %-contaminated pool has been shipping.
