# Robust strand-cleaning — precision-weighted concept (DEFERRED)

**Status:** deferred concept capture. 2026-06-10. Records the precision-weighted robust strand-clean
idea (prototyped, then reverted to keep a clean baseline) so it can be reconsidered when strand
deconvolution work resumes. Consolidates and supersedes the earlier `strand_cleaning_robustness_design.md`
(options analysis) and `strand_clean_global_target_design.md` (the g₀=global proposal). **Not active
work** — the current focus is the count channel (`count_channel_capture_design.md`).

## 1. The estimator and its fragility

The calibrator separates a node's unspliced mass into gDNA vs RNA by strand via a closed-form linear
unmix of the observed sense fraction `s = sense/total` between gDNA (sense rate ½, unstranded) and
RNA (sense rate κ = `rna_sense_frac`):

```
ĝ = (s − κ) / (½ − κ)            (method-of-moments gDNA fraction)
```

Unbiased, but its variance `∝ 1/[N(½−κ)²]` **explodes as κ→½** (no strand information) or N→0
(sparse). SS=½ is not an edge case — it is the maximum-uncertainty limit of a continuum, and the bare
ratio returns a finite, wildly amplified number exactly where it should return "I don't know." The
shipped code clips to `[0,1]` and special-cases `|½−κ| < 1e-6 → 1`, a discontinuous cliff that still
amplifies noise at κ=0.499.

## 2. The precision-weighted fix (the concept worth keeping)

Treat the clean as a *measurement* of g with precision equal to its own Fisher information
`τ = 4N(½−κ)²`, and shrink toward a fallback `g₀` with prior precision `τ₀`:

```
g_robust = [ τ·ĝ + τ₀·g₀ ] / [ τ + τ₀ ]
```

Substituting ĝ and τ, **the 1/(½−κ) cancels** — no division by the vanishing denominator:

```
g_robust = [ 4N(½−κ)(s−κ) + τ₀·g₀ ] / [ 4N(½−κ)² + τ₀ ]
```

Properties: smooth and well-defined at κ=½ (→ g₀); → ĝ when strand is strong (`τ ≫ τ₀`); no clip
cliff, no special case. `τ₀` is in units of strand Fisher information — `τ₀=1` ≈ "one fully-stranded
fragment," so the prior only bites when the node's own strand info `N(2κ−1)²` is comparably tiny. This
is the Gaussian/conjugate posterior mean for g under a strand-count likelihood + a prior at g₀; the
exact Beta-Binomial version is what the joint deconv already does for the *final* estimate, so this is
the cheap pre-deconv approximation.

The prototype implemented exactly this as `strand_clean_gdna_frac(sense, total, κ, *, tau0)`:
`data = 4(½−κ)(sense − κ·total)`, `tau = 4·total·(½−κ)²`, `g = clip((data + τ₀)/(tau + τ₀), 0, 1)`;
`τ₀=0` reproduces the legacy unshrunk clip exactly.

## 3. Why it was reverted — the target g₀ is the hard part, and it must be LOCAL

The mechanism is sound (smooth shrinkage, no high-SS regression), but the **fallback target g₀** is
where it lives or dies, and two choices were tried and rejected:

- **g₀ = 1** ("keep raw count"): behaviorally **inert** — on the 20-scenario benchmark it reproduced
  the old cliff's numbers, because "default to all-gDNA at κ=½" *is* what the cliff did. No change.
- **g₀ = library-global gDNA fraction:** **fails under hybrid capture.** Measured: the clean gDNA
  density varies ~8× across off-target observable regions and ~48× on/off-target under capture
  (`scripts/debug/diag_density_locality.py`). A global g₀ is meaningless when density is that local —
  shrinking an exon toward the library mean is wrong in both directions (over- and under-call).

**The lesson for whoever resumes this:** any strand-clean fallback target must be *local* (per-node /
locally-imputed), never global. Capture forbids a global density anywhere in the calibration. This is
the same constraint that reshaped the count channel (`count_channel_capture_design.md` §2.1: locality
in count space, means from boundary-anchored imputation).

## 4. The deeper open question (also deferred)

A standing architectural concern: the current pipeline is **strand pre-clean → joint count×strand
deconv**, and the pre-clean's strand use plus the deconv's strand likelihood may amount to using
strand *twice* — i.e. the "joint" model may really be strand-then-count in series. If strand
deconvolution is rebuilt, the cleaner target is a single robust strand step feeding the count step
(no double-use), with this precision-weighted shrinkage as the robust strand estimator and a
**local** g₀. Resolve the double-use and the local-target question together before reinvesting.

## 5. Status

Deferred. The precision-weighted estimator and the "fallback must be local" lesson are the keepers;
the g₀=1 and g₀=global trials are dead ends, recorded so they are not re-attempted.
