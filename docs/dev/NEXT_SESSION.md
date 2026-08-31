# NEXT SESSION — the state after the 2026-08-31 fragment-length deconvolution

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand, not what to build — the ranked
    list is `ROADMAP.md` §1. MOVE anything that settles into the permanent docs and DELETE this file.

## What the last three sessions settled

1. **The gDNA strand overdispersion — SHIPPED.** Holds regardless of the annotation: the away-half
   moment, influence weighting, the two components reconciling. `EQUATIONS.md` §6a/§6b/§6c,
   `DESIGN.md` §3.3a. The premise it replaced is `TRAPS: purity-is-a-property-of-the-annotation`.
2. **The gDNA fragment-length model — PRICED**, and the pricing re-ranked it: not a calibration defect
   (a perfect pmf moves the 0.8.0 metric less than the benchmark resolves, in an inconsistent direction
   between the two fl-gap arms), but worth **8–32 % of the library gDNA bias** through the EM.
3. ⭐⭐⭐ **The gDNA fragment-length model — SHIPPED (2026-08-31).** `ROADMAP.md` §0 carries the state;
   `src/rigel/calibration/gdna_density.py` and `fl.py` carry the mechanism; the derivation and the
   refuted routes are `PLAN_gdna_fl_estimator.md` §3/§4/§7.1 and its §9.

## What the fl fix actually is, in three lines

The two CONTAINED pools are the same two shapes at different mixing weights, so a **contrast** cancels the
contaminant with **no template**: `g = [(1−a₁)f₀ − (1−a₀)f₁]/(a₀−a₁)`. The weights need one scalar, the
gDNA density, read off the **LOW side** of the per-object density where a contaminant that only ADDS
cannot reach — the root of a one-sided Poisson moment with a closed form (De Moivre). **No constant
anywhere**, and no accumulator or schema change.

## ✅ SHIPPED — the owner's call, 2026-08-31, with the equal-length cost accepted

A pre-registered criterion said "the estimator must be inert where nothing is wrong", and on the ladder's
equal-length capture-OFF rows it is not, quite: three degrade at **1.25-1.7x the reseed floor** (`g98`
inside it, and the `g00` zero control **exactly unchanged** because it declines). ⭐ That cost is the
expected variance and not a defect — with equal component lengths there is no length bias to remove, so
only the contrast's ~2x variance remains. **The owner accepted the trade**, on these grounds: the ladder
TOTAL favours it by ~9x (both capture-ON rows improve well outside the floor — 8.2 % of the standing bias
at `g50`, **51x the floor**), and wherever the two components' lengths actually differ, as on real cfRNA,
the gain is 20-250x larger.

## ⭐ Why it survives hybrid capture — SOLVED, and worth understanding before touching it

The shared-contaminant premise is FALSE under capture (TV between the two pools' contaminants is 0.95,
against 0.06-0.14 off capture), and the estimator is good there anyway. **Not luck, and no longer
unexplained**: under capture the intergenic pool is *depleted, not impure* (measured purity 0.48 -> 0.955
at `g05`), so `a_0 = rho*E_0/n_0` exceeds 1 and CLIPS. At `a_0 = 1` the contrast collapses algebraically
to `g = f_0` — the intronic pool's coefficient is exactly zero, so **the term the premise was needed for is
removed, not approximated**. Verified to 3e-17.

⛔ **So the safety rests on `a_0` clipping**, i.e. on intergenic really being near-pure under capture. A
probe panel that put RNA back into intergenic space would break it silently. That is the thing to watch,
and it is the reason the two refuted detectors are not missed: there is nothing to detect while the clip
holds.

⚠ Three smaller risks: the one-sided rate is biased **high** by 0.2-6.7 % and its scaling with the
contaminated fraction is unestablished; the decline rule does not fire at `g00` on the test chromosome
(harmless, measured — the derived repair is an SE carrying `rho`'s own relative error rather than dropping
it as common-mode); and **nothing has been run on real data** (`real-data-is-a-test-input`).

## The live thread

⭐⭐ **The message policy — the standing charter, still at rung 0.** `MessagePolicy` remains byte-identical
to `SilentPolicy` on all 30 test-chromosome conditions. The bar is unchanged: **win on unstranded, minimal
harm on stranded, never pooled.** Before building any mechanism, answer rank 1a's question with no solver:
**is the information a blind unstranded or AMBIG slot needs actually present in its neighbours?** The
anchored twin block was designed for exactly that.
⚠ Re-baseline first — the policy numbers in `ROADMAP.md` predate both the strand-estimator change and this
one.

## Cautions that will otherwise cost a day

* ⛔⛔ **COMPARE AGAINST WHAT SHIPS, NOT AGAINST YOUR PROTOTYPE'S BASELINE.** This session called capture-ON
  a blocker for two revisions of a plan on the strength of a prototype scored against the two CONTAINED
  pools; what ships is the FOUR-pool sum, and against that capture-ON is a large WIN. A wrong baseline
  produced a confident, wrong, blocking verdict that survived external review.
* ⛔ **The benchmark resolves a LARGE effect.** A 1e−5 nudge to a nuisance parameter moves 1/30 rows by
  more than 0.5 %. Do not believe a single-row policy difference below ~2 %.
* ⛔ **ALWAYS RUN BOTH fl-GAP SIGN ARMS AND THE EQUAL-LENGTH CONTROL.** The same change has been a win on
  `rna_long` and a loss on `gdna_long` at both the calibration and the transcript level.
* ⛔ **`np.diff(payload.region_bounds)` is WRONG** — bounds are concatenated per reference, so it
  manufactures a phantom region at every junction. Use `gdna_density.region_lengths_from_partition`. ⭐ No
  shipped code had this bug (audited); a prototype did, and it made an estimator look 6× worse than it is.
* ⛔ **`docs/dev/` is a sandbox and nothing may cite into it** — a permanent doc that does is a suite
  FAILURE (`tests/test_docs_boundary.py`), which caught exactly that in a `ROADMAP.md` edit.
