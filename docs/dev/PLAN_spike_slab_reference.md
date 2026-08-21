# PLAN — the spike-and-slab reference: ONE new function feeding machinery that already exists

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** The derivation is `EQUATIONS.md`
    §9d, the ruling `DESIGN.md` §0c.3, the refutations `ROADMAP.md` §4.2/§4.3 and §0. Execute the
    stages below and move what settles; delete this file when the build lands or is refused.

    Written 2026-08-20 (session 4), at the owner's request: design before build — simple, elegant,
    readable, minimal regression risk. No `src/` change was made for this plan.

## 0. THE DESIGN IN FIVE LINES

1. ψ already owns everything hard: `_location_term` is a per-slot reference location with **bounded
   range, untouched tails, one-pseudo-observation strength, exact neutral collapse, and
   median(f) = m in closed form** — all gated. `CompositionPriors` already slices and regrids it
   safely. The refuted arms went through this same door, so the machinery is proven end to end.
2. The spike-and-slab is therefore **not a new prior engine**. It is: (a) ONE pure builder returning
   three per-slot scalars `(m_spike, m_slab, pi)`; (b) `_location_term` generalised by one
   `logaddexp` to a two-component mixture; (c) one config value.
3. **The slab collapses to a point in closed form**: the median of a log-uniform slab is the
   GEOMETRIC MEAN of its two measured endpoints — no grid, no integration, no constant.
4. **The capture-OFF limit holds by construction, not by branch**: the slab's upper endpoint is
   scaled by the library's own measured enrichment ratio `R` (off-target vs in-gene anchor,
   0.98 without probes / 113 with), so at `R ≈ 1` the slab degenerates onto the spike and the
   mixture IS the validated capture-OFF per-object location.
5. Ships behind `composition_reference: "structural" | "spike_slab"`, default `"structural"` —
   byte-identical until the owner flips it, priced as an arm first (`TRAPS: panel-before-src`).

## 1. WHAT ALREADY EXISTS — read before proposing anything larger

* `simplex_logodds._location_term(lam, location)` — the reference MEAN per slot, entering as
  `−log[(1−m)·f + m·(1−f)]`: range bounded by `log(max(m,1−m)/min(m,1−m))`, tails `e^(−|λ|/2)`
  untouched for every `m` (the module's own L-invariance acceptance), proper for `m ∈ (0,1)` with a
  derived clamp, exact zero at `m = ½`, `median(f) = m` in closed form. Strength is ONE
  pseudo-observation and is NOT re-opened (`EQUATIONS.md` §9c.1).
* `CompositionPriors.location` — a per-slot SCALAR, deliberately un-regriddable; `select()`/
  `regrid()` carry every member together so a hand-threaded second array cannot drift.
* `sweep.solve_chain:408` — the ONE construction site, already branching on `structural_reference`.
* Everything the builder needs is in scope there: `chain` (adjacency), `statics` (`mrna_active`,
  flags), `geometry` (`unspliced_count`, `eff_gdna`, `eff_rna`, the sj banks).

⭐⭐ **Why the mixture escapes the single-mode refutation (`ROADMAP.md` §4.2 rows, §0's 5–8×).**
`pooled`/`local` fed ONE location through this door, and at an enriched slot that location was the
depleted background — a ~6-nat claim on the wrong side, which at an evidence-poor slot IS the
answer. The mixture keeps BOTH hypotheses alive as two bounded bumps on the λ grid; the slot's own
evidence picks between them **at solve time** instead of the prior choosing a priori. ⛔ That is a
HYPOTHESIS until priced — the staged plan below prices it before anything defaults on.

## 2. THE MATH — every quantity measured, no constant introduced

Per slot `i` with unspliced mass `M_i` and gDNA opportunity `E_i` (both already in `geometry`),
restricted to `mrna_active` slots (¬`mrna_active` keeps the shipped structural 0.75 — scope
unchanged):

| symbol | meaning | source |
|---|---|---|
| `rho_0` | the un-enriched gDNA density | pooled `Σcount / ΣE_g` over ¬`mrna_active` slots — the ratio of sums, never a mean of ratios |
| `R` | the library's measured enrichment span | pooled in-gene anchor density (`exon\|intron` boundaries: ¬`mrna_active` slots with a live sj flank inside a gene) ÷ `rho_0`. Measured 0.98 without probes, 113–114 with (`DESIGN.md` §6b.1) — a DETECTOR that needs no threshold |
| `rho_lo,i` | the slab's floor | the slot's adjacent `exon\|intron` flank density where one exists (enrichment is monotone in probe proximity), else `rho_0` |
| `rho_hi,i` | the slab's ceiling | `min( M_i / E_i , rho_lo,i · R )` — an object cannot hold more gDNA than the mass it holds, and cannot be more enriched than the library's own measured span |
| `pi` | the unprobed weight | ½ at pass-0 — the reference's own Jeffreys convention (`DESIGN.md` §0c.3), nothing new |

Then, in composition coordinates (`f = rho·E/M`, clipped by the existing derived clamp):

    m_spike,i  =  rho_0      · E_i / M_i
    m_slab,i   =  sqrt(rho_lo,i · rho_hi,i) · E_i / M_i        ← the log-uniform slab's own MEDIAN
    term_i(λ)  =  logaddexp( log(pi)   + T(λ; m_spike,i),
                             log(1−pi) + T(λ; m_slab,i) )      ← T = the existing bounded form

Properties inherited for free: each component is the gated bounded form, so the mixture's range is
bounded by the larger component and the tails stay `e^(−|λ|/2)`; the AMBIG float32 safety argument
carries over unchanged. The exact-zero rule generalises: both `m = ½` ⇒ return exact 0 (the same
knife-edge protection, same reason).

**The three limits, by construction (each becomes a unit gate):**

* **capture-OFF** — `R ≈ 1` ⇒ `rho_hi → rho_lo → rho_0` ⇒ `m_slab → m_spike` ⇒ the mixture IS the
  single validated capture-OFF per-object location. No branch, no flag.
* **`g00`** — the anchors are empty ⇒ `rho_0 = 0`, `rho_lo = 0` ⇒ both components at the clamp ⇒
  `m → 0`: the vertex is reached because a measurement put the atom there (this exact plumbing
  already measured "zero controls to EXACTLY 0" in the per-object study).
* **`g98` capture-ON** — `rho_hi`'s `M/E` cap alone reads `f ≤ 1`; with the measured `R` the slab's
  median sits at `sqrt(rho_lo · min(M/E, rho_lo·R))·E/M`, within the truth's decade where the truth
  is near the vertex.

⚠ **Deliberately NOT in the first rung**: the §9d.3 flux-tightened upper endpoint (`− S` term). It
tightens `rho_hi` where an adjacent junction speaks, is exactly derived, and is the natural rung 2 —
but it adds a second population question (which flank's flux) to rung 1, and the owner's discipline
is one concept at a time. The `M/E` cap is the structural degenerate case and is enough to build on.

## 3. THE INTEGRATION DIFF — enumerated, ~120 lines total

| # | file | change |
|---|---|---|
| 1 | `simplex_logodds.py` | `spike_slab_reference_location(chain, statics, geometry) -> (m_spike, m_slab)` beside its sibling `structural_reference_location` — same home, same duck-typed style, no new module, no new import (layer 3 stays layer 3) |
| 2 | `simplex_logodds.py` | `CompositionPriors` gains `location_slab: np.ndarray | None = None` (per-slot scalar, same immunity argument); `select()` slices it, `regrid()` passes it through — the two one-line spots the class exists to keep together |
| 3 | `simplex_logodds.py` | `_location_term` grows the mixture branch: `location_slab is None` ⇒ today's path **bit-identical**; else the `logaddexp` above with `pi = ½` inline (the Jeffreys convention, not a parameter) |
| 4 | `sweep.py:408` | the construction site branches on the config value: `"structural"` ⇒ exactly today's call; `"spike_slab"` ⇒ the new builder for `mrna_active` slots + structural for the rest |
| 5 | `config.py` | `composition_reference: str = "structural"` beside `structural_reference` (⚠ reconcile: the new value subsumes the old flag's meaning at ¬`mrna_active` slots — rung 1 keeps BOTH fields with `structural_reference` untouched, and the merge is a later cleanup the owner rules on) |

**What is deliberately NOT built:** no grid-shaped per-slot prior arrays (the location stays a
scalar pair — the regrid-immunity argument is why), no new solver path, no message-policy change,
no landscape coupling (the second solve's landscape-vs-slab question is `docs` §0c.0d's circularity
discipline and stays open), no new module (the 38-module flat-pile lesson).

## 4. REGRESSION RISK, stated as mechanisms

* Default `"structural"` ⇒ the builder never runs and `_location_term`'s `None`-slab path is the
  shipped code — **byte-identity is a gate, not a hope** (rung 2 below).
* The mixture term's range is bounded by its larger component ⇒ the float32 AMBIG-cube argument
  (`_location_term`'s own docstring) carries unchanged.
* The knife-edge rule (exact 0 at neutral) is preserved by the same mechanism that already guards it.
* `CompositionPriors` extension follows the exact pattern the class was designed to make safe.
* The one genuinely new numeric surface is the builder — pure, array-in/array-out, unit-gated
  against hand-computed fixtures before it is ever wired.

## 5. THE STAGED BUILD — one concept per rung, a gate before each

1. **Rung 1 — the builder, unwired.** `spike_slab_reference_location` + unit gates written FIRST
   and verified failing: the three limits (`R≈1` collapse to the single location, exactly; empty
   anchors ⇒ clamp; the `M/E` cap), the geometric-mean median against brute-force quantile of the
   log-uniform, `R` on a two-density fixture, flank fallback, and the `mrna_active` scope asserted.
   Then break it three ways (swap lo/hi; arithmetic-for-geometric mean; drop the `R` cap) and watch
   each gate fire.
2. **Rung 2 — the term + wiring, default OFF.** The `CompositionPriors`/`_location_term`/`sweep`/
   `config` changes; gates: `location_slab=None` bit-identity on a real cached condition (the
   shipped tree reproduced exactly), the mixture's exact-zero at double-neutral, bounded range,
   and `TRAPS: an-ablation-that-never-ran` — the flag demonstrably fires when flipped.
3. **Rung 3 — price it as an ARM on all 16** (`TRAPS: panel-before-src` for the DEFAULT, not the
   code): per stratum, never pooled; the deferred row its own row; BOTH zero controls; the SHUFFLE
   control (permute `m_slab` across slots — must hurt); the evidence split on every row
   (`TRAPS: an-imputation-must-cost-something-every-hop`'s discipline applied to a prior); and the
   three refuted single-mode arms re-run beside it as comparison arms. Under a corrected reference
   the recorded expectation is that every policy improves at the controls — read the three-policy
   contrast again on top (`ROADMAP.md` §1 rank 1a), since the reference was compensating.
4. **Rung 4 — the owner rules on the default**, with the rung-3 table and the `g05` blocker's fate
   (the family's recorded killer: it must NOT regress, and now has its own rungs on the test
   chromosome for the fast loop).

## 6. OPEN DECISIONS, each with a recommendation

* **`pi` refinement** — the owner's total-abundance-landscape projection (a slot's TOTAL abundance
  is informative about ε for both components at once). Recommend: NOT in rung 1; it is the natural
  rung after the flux endpoint, and the owner said of the landscape idea "not essential for the
  first iteration".
* **The flux-tightened `rho_hi`** (§9d.3): rung 2 of the endpoint, after the panel prices rung 1.
* **`structural_reference` merge**: recommend deferring; two flags for one release cycle is
  cheaper than a premature unification.
* **Naming**: `composition_reference` / `spike_slab_reference_location` — no new vocabulary; "spike"
  and "slab" are already `DESIGN.md` §0c.3's ruled words.
