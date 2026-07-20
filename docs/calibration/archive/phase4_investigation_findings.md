# Phase 4 investigation — findings (P4a), and the refined plan (2026-07-17)

**Scope (owner-gated):** *investigate only, no production wiring changes.* Re-measure the refit / flagship
picture now that the DNA-background floor + pinned component (P2/P3) are wired into production, localize the
residual, and bring findings back before any structural change. Branch `calib-ambig-init-wip`.

**Method.** The refit harness (`scripts/debug/refit_loop_study.py`) and the two dissection tools
(`flagship_dissect.py`, `prior_crush_probe.py`) all fit the prior **without** the background — so every number
in `refit_loop_study_findings.md` predates P2/P3. First move was to thread the production background into all
three (each gained a `--no-floor` toggle) so they measure the *current* production regime. Substrate: the cached
`ambig_dense_10mb/_selfsolve_cache` conditions. `OMP_NUM_THREADS=1`.

---

## Finding 1 — the floor wiring is clean (no code-churn bug)

Floor-on vs floor-off, pass-0 mwae, across 8 conditions, reproduces the shipped P2/P3 A/B exactly: flagship
−0.003, stranded capON −0.004, stranded capOFF +0.001, **zero-gDNA and low-DNA-verystrong exactly 0.0000
(dormant)**. No surprise, no regression from the wiring. The owner's "there have to still be bugs from code
churn" suspicion does **not** land on the floor path — it behaves as derived.

## Finding 2 — the disagreement instrument: the flagship solution is 4–27× OVER-smoothed

The ground-truth-free disagreement lens (between-node squared adjacent log-density gap, per component, reported
against the oracle's value on the same edges) is the headline:

| condition | mwae (pass-0→final) | `disG_sol` (pass-0→final) | `disG_orc` |
|---|---|---|---|
| flagship present capON | 0.567 → 0.598 | 2.0 → **1.6** | **8.9** |
| flagship none capON | 0.625 → 0.661 | 2.1 → 1.7 | 6.8 |
| gdna100 flagship | 0.362 → 0.377 | 2.0 → 1.5 | **54.5** |
| gdna5 unstr verystrong | 0.031 → 0.035 | 22.8 → 10.0 | **278.7** |
| gdna1 unstr verystrong | 0.011 → 0.007 | 35.5 → 12.0 | 333.3 |
| stranded capON | 0.070 → 0.060 | 6.4 → 6.6 | 6.9 |
| stranded capOFF | 0.012 → 0.011 | 0.05 → 0.06 | 0.06 |

Two things:

1. **The flagship gDNA field is far too flat.** `disG_sol ≈ 2` vs oracle `≈ 9` (and 27× / 28× under for
   gdna100 / verystrong). The real enrichment crossings — gDNA stepping up from a depleted intron/intergenic
   region into a capture-enriched exon — are **missing from the solution**; they have leaked into RNA. Stranded
   data, by contrast, matches the oracle floor (`disG_sol ≈ disG_orc`), because strand resolves the crossings.
2. **The refit drives `disG_sol` the WRONG way — toward zero, away from the oracle.** Every flagship row loses
   disagreement across refit passes (2.0→1.6, 22.8→10.0, 35.5→12.0) and mwae drifts up. The floor does not
   arrest it.

**This flips the optimization objective.** "Minimize between-node disagreement" — taken literally — pushes the
flagship *deeper* into the degenerate over-smoothed basin (the refit does exactly that). The correct target is
**match the oracle floor** (`disG_sol → disG_orc`), and the flagship needs its gDNA disagreement to go **UP** by
4–27×. Truth-free deployability: on real data we lose `disG_orc`, but the **total-density** disagreement is a
truth-free upper reference for how much real structure exists — a solved gDNA field far flatter than the total
density, with strand uninformative, is the diagnostic.

## Finding 3 — localization: enriched, truly-gDNA, strand-blind exons crushed to ~0

`flagship_dissect.py` (floor on) buckets the pass-0 error by structural class. The residual is concentrated in
`single` (69% of mass) and `AMBIG` (19%) nodes with a **signed** error of −0.58 / −0.67 (systematic *under*-call
of gDNA), and of these, ~40–49% of the mass is **truly gDNA** (oracle `fo > 0.8`). The top error nodes are
unambiguous: e.g. node 1909 is 98% gDNA by the oracle (G=32651/32861 per strand, R=498/537), strand is exactly
uninformative (`str=0.49`), and the solve crushes it to `f_g = 0.00`.

The floor cannot touch these: it activates only where `ρ_g < ρ_bg`, i.e. `f_g < ρ_bg/ρ_tot`. For node 1909,
`ρ_bg ≈ 0.05` (capON depleted off-target) and `ρ_tot ≈ 34` (enriched) ⇒ the floor bites only below
`f_g ≈ 0.0015` — far below the action. **This is the honest limit, measured: the off-target floor structurally
cannot lift an enriched on-target gDNA node.** The received gDNA messages on these nodes are gagged to precision
≈ 0.05 (σ²_transfer, correctly, sees the enrichment crossing) and cannot lift them either.

## Finding 4 — the crusher is the RNA reference arm, not the gDNA prior

`prior_crush_probe.py` decomposes ψ on the worst nodes (floor on). For node 1909:

* **gDNA prior arm** (total-density NPMLE): *favors* high f_g by **+0.60** nats (`f_g=0.9` vs `0.1`).
* **RNA reference arm** `½·log(1−f_g)` (the improper Jeffreys reference, the only thing bounding the `f_g→1`
  vertex today): *penalizes* high f_g by **−1.10** nats.
* Net ψ(0.98) − ψ(0.1) = **−1.33 nats** ⇒ argmax at `f_g = 0.002`. **The RNA arm wins and manufactures RNA.**

So the enriched-node crush is the **asymmetric-ψ bug still live** — but the active crusher on enriched nodes is
the **RNA arm's fixed barrier**, not the gDNA arm's clamped tail (the gDNA arm here actually points the right
way). The gDNA arm is simply too weak (±0.6) to overcome the −1.10 barrier.

## Finding 5 — why V4 (fixed gDNA barrier) and V2 (symmetric arms) both fail — hard control

The same probe evaluates candidate ψ forms on the enriched node (gdna300) **and the identical region in
zero-gDNA** (gdna_none, where node 1909 is now 100% RNA, `ρ_tot ≈ 0.55`):

| variant | gdna300 node 1909 (truth f_g≈0.98) | gdna_none node 1909 (truth f_g≈0) |
|---|---|---|
| V0 current (clamp + Jeffreys) | mean **0.20** (crushed) | mean **0.06** ✓ |
| V2 symmetric total-NPMLE | mean 0.50 (ambiguous) | mean **0.50 / argmax 1.0** ✗ FALSE gDNA |
| V4 restore `½log f_g` barrier | mean 0.60 (fixed) | mean **0.38** ✗ FALSE gDNA |

* **V2 is symmetric under `f_g ↔ 1−f_g`** (both arms the same NPMLE) ⇒ it returns mean 0.5 for *every* unstranded
  node regardless of gDNA content — a catastrophic false positive when the truth is RNA. Symmetric arms just
  move the failure from "crush everything" to "50/50 everything."
* **V4's fixed barrier manufactures gDNA** wherever total density is non-trivial, including zero-gDNA (already
  rejected in the 24-scenario A/B; re-confirmed here at the node level).

**The decisive observation:** the gDNA prior arm **already discriminates correctly at the population level** — it
*favors* high f_g by +0.60 on the enriched flagship node but *disfavors* by −0.76 on the zero-gDNA node (the
total-density NPMLE encodes "enriched ⇒ maybe gDNA, depleted ⇒ not"). It is merely too weak to beat the RNA
barrier. So the fix is **neither a fixed barrier nor symmetric arms** — it is a **data-gated strengthening of the
gDNA arm on genuinely-enriched-gDNA nodes**, which vanishes in zero-gDNA. That is exactly the derivation's
**on-target abundance anchor**.

## Finding 6 — the message ceiling (σ²_transfer is the other half, but regime-signed)

Turning σ²_transfer fully OFF (un-gagging the crossing messages), floor on: flagship mwae **0.567 → 0.439** and
`disG_sol 2.0 → 5.75` (toward oracle 8.9) — proving the missing crossings **are** carried by the messages
σ²_transfer gags. But σ²_T OFF wrecks the stranded regime (0.06 → 0.22–0.45, warm-continue period-2 cycle). So
the message lever cannot be a global switch; it must be **selective — keyed to strand-uninformativeness**. That
is the deferred P4d (it touches the reviewer-signed-off F1) and is complementary to the abundance anchor: even
fully un-gagged, the messages carry the *depleted-neighbor* level across the enrichment step, so they still need
the on-target *level* from somewhere.

---

## Refined Phase 4 plan

The investigation confirms the derivation's target and sharpens the mechanism. The flagship residual is enriched,
truly-gDNA, strand-blind exons that ψ crushes to ~0 because the RNA reference arm (−1.10) beats the too-weak
gDNA arm (+0.60); the off-target floor cannot reach them and the corrective messages are gagged. Both hard
controls (V2, V4) rule out any fixed/symmetric barrier. The one lever that is **data-gated** (so it cannot break
zero-gDNA) is:

**P4c — the on-target enriched-gDNA abundance anchor `ρ_on`.** A one-sided lower anchor at the *enriched* gDNA
level (the analogue of the `ρ_bg` floor, but at the on-target level), which strengthens the gDNA arm on
genuinely-enriched nodes enough to overcome the RNA barrier, and which → 0 in zero-gDNA (no false positive). It
supplies the *level* the gagged messages cannot carry across the enrichment step.

**The open problem is the estimator: how to measure `ρ_on` non-circularly**, since on-target gDNA is entangled
with abundant on-target RNA and strand cannot separate them in the flagship. This must be *derived before any
code* (owner directive: derive, don't band-aid; no magic numbers). Candidate directions to work through:

1. **Structurally-gDNA-only anchors.** The `gonly` nodes (neither strand free ⇒ no RNA can flow ⇒ structurally
   gDNA) solve *perfectly* today (0.000 error). They are on-target gDNA observations that could seed `ρ_on` —
   the question is whether enough of them are capture-enriched to represent the on-target level (they are only
   ~4% of mass), and how to aggregate (a scalar? an enrichment-vs-total-density relation?).
2. **The gDNA-vs-total enrichment relation.** gDNA (unlike RNA) is present off-target too; the off-target `ρ_bg`
   plus the on-target anchors may define an enrichment curve `ρ_g(total density)` that is readable at the
   population level without per-node circularity.
3. **Iterative, floor-anchored (the derivation's framing).** Estimate `ρ_on` from the solved on-target gDNA each
   pass — but only from a *trustworthy* subset (NOT belief-confidence, which Finding 5 of the refit study showed
   is confidently-wrong in the flagship), with cold-restart, gated on `disG_sol → disG_orc` (Finding 2). This is
   where the refit belongs — as a global-scalar estimate anchored by the floor, **not** the whole-landscape
   refit that Findings 2/§refit proved harmful.

**Complementary, deferred:** P4d (selective σ²_transfer un-gag, Finding 6) supplies the spatial structure once
the level is right; it touches F1 and gets its own derivation + reviewer loop. A data-anchored `logP_r` (Finding
4) is the same idea from the RNA side, but Finding 5 shows it must be anchored on *RNA-specific* evidence
(mature/spliced), **not** a shared total-NPMLE — otherwise it collapses to the symmetric-0.5 failure.

**Recommendation:** proceed to **derive the `ρ_on` estimator** (candidates 1–3 above) as a reviewer-ready
document for owner sign-off before writing code — mirroring how the `ρ_bg` background reference was derived and
reviewed. Also fold the truth-free disagreement reference (`disG_sol` vs total-density `disG`) into the
instrument so P4c has a real-data gate.

---

## Addendum — production-faithful re-interrogation (after owner distrust of the harness)

The harness (`refit_loop_study._setup` + pass-0 sweep) was checked against the **real** `calibrate()` pass-0
(`calib_refit_iters=0`, `_debug` hook): **`max|Δf_g| = 0.0`, byte-for-byte identical** — the reconstruction had
not drifted. But two real gaps surfaced: (a) **production ships `calib_refit_iters=1`**, and that refit moves
`f_g` by up to **0.996** (mean 0.026) on the flagship — my node dissections only ever looked at pass-0; (b) the
`scripts/debug` tree genuinely mixes current and retired tools. A single production-driven tool now supersedes
the pass-0 reconstructions: **`scripts/debug/flagship_interrogate.py`** drives the real `calibrate()` (pass-0
*and* shipped), traces each node from the production `_debug` capture, and interrogates chain NEIGHBOURS.

**Population evidence (shipped belief, not anecdote):**
* **35% of all live mass** is enriched, truly-gDNA nodes crushed to `f_g < 0.2` (`n=364`). This *is* the flagship
  error (mwae 0.59).
* **98% of that crushed mass has NO enriched-gDNA neighbour** — they are islands. No spatial message can lift
  them, at any σ²_transfer setting, because no neighbour carries the enriched gDNA level. (Neighbours are the
  exon-intron boundary nodes, which sit at the depleted level `log10 ρ ≈ −1.3` or are themselves crushed.)
* The perfectly-solved `gonly` nodes (structural gDNA-only, 4% mass) are **all at the depleted level**
  (`log10 ρtot ∈ [−1.30, −1.28, −1.26]`, 0% enriched). **There is no clean on-target gDNA anchor.** The
  "aggregate the `gonly` anchors" estimator (candidate 1 of the plan above) is **dead**.
* The enriched gDNA regions carry **exactly zero spliced reads** (`mat_spl` median 0.0).

**The identifiability picture this forces.** The 35%-of-mass error is enriched, **unspliced**, **50/50-strand**,
**spliced-free** exons that are truly gDNA. In an unstranded library these are **fundamentally ambiguous between
gDNA and nascent RNA**, because every calibration signal is silent on that split: strand is 50/50 (unstranded);
spliced=0 rules out *mature* RNA but nascent RNA is also unspliced; the count carries zero info; and the spatial
messages are islands. The only remaining discriminators are the **fragment-length model** (gDNA vs RNA FL) and a
**global gDNA-abundance prior**. Today the RNA reference arm `½·log(1−f_g)` resolves the ambiguity toward RNA
**by fiat** (a fixed, data-independent `f_g→1` barrier) — that is the crush. This is also exactly why V4/V2
failed the zero-gDNA control: they resolve the same ambiguity toward gDNA by fiat. **The correct resolver is the
library's data-driven global gDNA abundance** (gDNA-rich ⇒ these ambiguous unspliced nodes are probably gDNA;
gDNA-poor ⇒ probably nascent), which is the derivation's thesis — but it needs (a) a non-circular estimate of
that global abundance (NOT from clean anchors — there are none in the fully-unstranded corner; likely from the
FL model + off-target `ρ_bg` + the strand-identifiable fraction where it exists) and (b) a fix to the RNA arm so
an abundance-informed gDNA prior can actually win instead of being overridden by the fixed barrier.

**Two decisive controls close the identifiability question:**

* **Strand recovers them completely.** The *same* nodes under stranded `ss_0.99` (same reference, same chain):
  the crushed-enriched-gDNA population collapses from `n=364 / 35%` of mass to `n=14 / 0.1%`, shipped mwae
  0.59 → 0.058. Strand *is* the missing information — nascent RNA is strand-specific, gDNA is not, so a stranded
  (or mostly-stranded) library resolves the gDNA/nascent confound directly, with no new machinery. The fully
  unstranded `ss_0.50` corner deletes exactly that signal.
* **Fragment length is a strong, strand-independent discriminator — but not wired into calibration.** The gDNA
  and RNA FL distributions are nearly disjoint (gDNA ~111 bp sd 32, RNA ~200 bp sd 5; TV distance **0.886**).
  A per-node fragment-length likelihood would separate gDNA from nascent at these nodes *without* strand. But
  the accumulator deposits only **6 region-type-stratified aggregate** FL pools (`fl_pool_mass`:
  {intergenic/intronic/exonic} × {contained/boundary}) used to *train* the gDNA/RNA FL models — **there is no
  per-node FL histogram**, and the count-zero-info deconvolution uses only strand + messages + prior, never a
  per-node FL likelihood. (Caveat: the sim's RNA FL sd=5 is unrealistically tight; real RNA insert sizes spread
  more, so 0.886 is likely sim-optimistic.)

**Where this leaves Phase 4 (for the owner to steer).** The `ss_0.50`-exact flagship enriched-gDNA crush is, with
the *current* channels, near the identifiability frontier: strand is deleted, spliced=0, count is zero-info, the
messages are islands, and per-node FL is not deposited. Three real directions, not one:

1. **Target the mostly-stranded + AMBIG regime** (the "always test AMBIG with ample single-strand nodes" rule).
   The stranded majority both solves its own nodes *and* anchors a global gDNA-abundance estimate that resolves
   the AMBIG subset. `ss_0.50`-exact stays a stress corner, not the ship bar. Cheapest; uses signal we have.
2. **Add a per-node fragment-length channel** to the calibration deconvolution — a fourth composition source
   (strand + FL + messages + prior). Strong and strand-orthogonal, but requires an accumulator schema change
   (deposit per-node/finer FL histograms, byte-for-byte reference update) + a solver change. Biggest lift;
   most fundamental; independent of strandedness.
3. **Global gDNA-abundance prior only** (the derivation's original P4c) — but with no clean anchor it must lean
   on `ρ_bg` + the FL-trained gDNA fraction + a uniformity assumption, and it is a weak *population* prior, not a
   per-node resolver; it cannot install the sharp enriched crossings the truth has (Finding 2).

## Reproduce

```
OMP_NUM_THREADS=1 python scripts/debug/refit_loop_study.py --compact --passes 5 \
    --conditions <flagship,stranded,zero,lowDNA...>            # floor ON (production)
OMP_NUM_THREADS=1 python scripts/debug/refit_loop_study.py --compact --no-floor ...   # pre-P2/P3
OMP_NUM_THREADS=1 python scripts/debug/refit_loop_study.py --no-tv --conditions <flagship,stranded>  # msg ceiling
OMP_NUM_THREADS=1 python scripts/debug/flagship_dissect.py --condition <flagship> --top 12
OMP_NUM_THREADS=1 python scripts/debug/prior_crush_probe.py --condition <flagship|none> --nodes 1909,2337,1523
```
