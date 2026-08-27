# THE MESSAGE-POLICY CAMPAIGN — what was tried, measured, and REFUTED (closed 2026-08-27)

    ⚠ **A DEV DOC. It is a RECORD, not a plan.** The code every entry below describes was
    DELETED on 2026-08-27; git carries it. This file exists for one purpose: so a refuted
    experiment is not run again. Nothing here is an instruction to build anything.

A campaign to replace the shipped `RelayPolicy` with a derived message policy ran for several
sessions and did not reach the bar. Its architecture survived — `messages/foundation.py`, the
ratified spec, and `MessagePolicy`, the runner on top of it. Its mechanisms did not.

## The bar, and why the campaign missed it

⛔ **The bar is NOT to beat `SilentPolicy`.** Silence is the floor; on strand-specific data a
sighted exon's own solve is excellent and a message can mostly only disturb it. The goal is to
perform on UNSTRANDED data while doing minimal harm on stranded data. The campaign repeatedly
optimised a pooled total instead, which hides a sign flip between the two halves.

⭐ **The one measurement that survived everything**: propagation is net-harmful wherever the local
solve HAS its own evidence, and its value concentrates where that solve is BLIND. That is the whole
reason the two halves are judged against different bars.

## REFUTED — do not re-run these without new evidence

| mechanism | what happened |
|---|---|
| **A composition-transporting policy** (`CurrencyPolicy`) | Best zero controls ever measured, but lost every in-scope contaminated stratum to silence. Deleted. |
| **The gDNA-continuity rule** (an unsupplied source's gDNA level crosses unscaled) | Built THREE ways — a static per-slot licence (a value RATCHET, gdna densities to 3.9e+32, from breaking the knob's telescoping cancellation), a running-state licence (killed the ratchet, still lost capture-ON), and a fuse-based pure-gDNA re-anchor lattice (too weak — a fuse negotiates where the relay's mass rescale overwrites). Halved one row, lost others. It is only safe beside a scan-time mass rescale. |
| **The premise's exon-end scoping** | The dispersion decomposition proves intron-end hops carry no COMPOSITION cost, yet scoping the charge to exon-end hops REGRESSED the panel: freeing intron chains before a measured LEVEL charge exists releases un-priced level drift. Restoring the pooled charge recovered g98-ss0.50-ON 4.37M → 2.44M. |
| **A class-keyed method-of-moments fit on the observed log-ratio** (as the runtime law for the transport variance) | Tracks truth at intron and plain classes, REFUTED at sj classes: the route-summed flux cancels the visible step exactly where the true error is largest (0.08 observed vs 4.4 true). |
| **A totals-form pair fit** (for the same variance) | Refuted by construction — the knob consumes the totals, so transported totals agree by the (1−w) algebra and carry almost no information. |
| **`FanOutPolicy`** | Measured dominated once the certified-flux anchor gave its destinations own evidence. Deleted 2026-08-24. |

## Derivations worth keeping (kept so they are not re-derived)

* **The conservation identity is a COUNT identity**: `Σ_c ρ_c·E_c = M` — each component's density
  weighted by its OWN opportunity, summing to the slot's observed unspliced count. Summing raw
  densities against the reciprocal-opportunity total is a DIFFERENT identity that is exact at
  boundaries and WRONG at regions, where `inv_abundance` reads `ρ·P(w ≤ ℓ)` — measured on the
  ladder, the median exon's P is 0.452 and the 5th percentile 0.044, so that total was up to 23×
  too low at half of all exon slots.
* **Two variance kinds.** A reframe's cost multiplies every lane of a message identically (a LEVEL
  statement); each component also carries its own (a COMPOSITION statement). Spending the shared
  part as per-component variance converts a common-mode level error into a composition error — the
  measured pathology where a near-zero gDNA claim eats an unstranded slot's unexplained RNA mass.
  A conservation solve carrying both kinds has the multiplicative rescale, the additive
  precision-weighted allocation, and "the unexplained mass is RNA" as its three LIMITS.
* **The transport dispersion, decomposed against certified truth** (four g50 corners, noise
  subtracted): intron↔boundary hops are FREE; all structural error lives on exon↔boundary hops and
  is predominantly COMMON-MODE, exploding under capture. Two sources are identified and both are
  derivable rather than fittable — the truncation frame term `log(1/P(w ≤ ℓ))` (pure geometry) and
  the capture step (common-mode because probes bind gDNA and RNA alike).

## Method lessons (the durable ones already have named homes in `TRAPS.md`)

* ⭐ **Attribute before iterating.** Three laws landed together, the panel regressed, and only an
  attribution factorial — each law removed alone, all 16 conditions — named the carrier. Without
  it the session would have reverted the wrong mechanism.
* ⭐ **A cancelling defect pair reads as success.** Fixing the conservation identity made low-gDNA
  WORSE, because the truncated total had been suppressing claims at short exons and hiding phantom
  gDNA. A correct fix can look like a regression when it removes half of a cancelling pair.
* ⛔ **Layering fixes on fixes is how a policy becomes unmaintainable.** The campaign ended with
  three conservation operators behind flags, two of which were provably limits of the third.
