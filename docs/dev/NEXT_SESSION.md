# NEXT SESSION — the state after the 2026-08-31 fragment-length work

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand, not what to build — the ranked
    list is `ROADMAP.md` §1. MOVE anything that settles into the permanent docs and DELETE this file.

## Where the tool is

Two multi-session threads have closed. The gDNA **strand overdispersion** now holds regardless of the
annotation (`EQUATIONS.md` §6a/§6b/§6c, `DESIGN.md` §3.3a). The gDNA **fragment-length model** is
deconvolved, its two estimands are separated, routed and coupled, and the state is `ROADMAP.md` §0's
fragment-length row — read that first, it is the single home for the numbers.

Six commits this session: `a4b6a134` (the contained contrast), `7314632d` (a closure later reversed),
`7b3c70ba` (the reversal), `34bae6b4` (the frame diagnosis), `496c1b81` (routing), `0bbabc51` (the
coupling that removed the last three thresholds).

## ⭐⭐⭐ START HERE: the RNA length law — `PLAN_rna_length_law.md`

The one remaining fragment-length defect, diagnosed to its mechanism and **worth an attributable −18,208
transcripts**. That plan doc has the derivation, the refuted route, the fix, the sequencing and the open
issues. In three lines:

* `rna_pmf` is **−4 bp short** because the sj de-tilt divides by the probability a fragment *crosses* a
  junction while the pool is selected on the splice being **OBSERVED**, and observability falls with
  length. Plus a smaller **estimand mismatch**: it is a MATURE law used for ALL RNA.
* ⛔ **No estimator built on spliced fragments escapes this** — the sj bank's model-free identity is
  *equally* biased (−4.6, −5.3), measured. Do not retry it.
* ⭐ **The escape is already computed and discarded**: the two-pool contrast solves for the contaminant
  law `r_hat`, which has no observability selection and lands **inside 0.2 bp** of truth.

⛔ **Step 0 of that plan is not optional**: re-measure the three-way consumer split, because the
scorer/drain/geometry numbers predate the routing and coupling that landed after them.

## The other live thread

⭐⭐ **The message policy — the standing charter, still at rung 0.** `MessagePolicy` remains byte-identical
to `SilentPolicy` on all 30 test-chromosome conditions. The bar is unchanged: **win on unstranded, minimal
harm on stranded, never pooled.** Before building any mechanism, answer `ROADMAP.md` §1 rank 1a's question
with no solver: **is the information a blind unstranded or AMBIG slot actually present in its neighbours?**
The anchored twin block was designed for exactly that. ⚠ Re-baseline first — the policy numbers in
`ROADMAP.md` predate both this session and the strand-estimator change.

## ⛔ Standing risks from this session's work

* **The gDNA contrast survives capture by a degeneracy, not by its premise.** Under capture the
  shared-contaminant assumption is FALSE (TV 0.95), and it is safe only because the intergenic pool is
  *depleted, not impure*, so `a_0` clips to 1 and the algebra collapses to `g = f_0`. **A probe panel that
  put RNA back into intergenic space would break it silently.** `ROADMAP.md` §0.
* **Geometry is knife-edge sensitive to the pmf tail** — ±188k / −15.5k transcripts from small tail
  changes, because the short-exon contained opportunity falls 66 % at `ell = 100`. Any change to an
  `effective_lengths_em` input needs the reseed floor beside it. That is effective-length-ruler territory
  and now has a mapped lever.
* **Nothing has been run on real data** (`real-data-is-a-test-input`). The cfRNA libraries are where a
  genuinely different contamination profile lives, and no fragment-length claim here has met one.

## Cautions that will otherwise cost a day

* ⛔⛔ **SCORE A CHANGE ON THE METRIC ITS CONSUMER READS.** This session closed the gDNA shape question on
  `gdna_frac_est` — which `em_fl_ceiling.py` prints under a "THE PRODUCT" banner — and was WRONG: the fl
  model's consumer is the EM's per-fragment assignment, and transcript accuracy said the opposite. The
  banner is right for a composition fix and wrong for a shape fix.
* ⛔⛔ **COMPARE AGAINST WHAT SHIPS, NOT YOUR PROTOTYPE'S BASELINE.** Capture-ON was called a blocker for
  two revisions of a plan because a prototype was scored against the two CONTAINED pools while production
  uses the FOUR-pool sum. A wrong baseline produced a confident, wrong verdict that survived external
  review.
* ⛔ **ONE NAME OVER TWO POPULATIONS is this codebase's recurring defect.** Three instances this session:
  the gDNA pmf (chemistry vs census), the "pure" pools (annotation vs genome), the RNA pmf
  (spliced-observed vs all-RNA). None was a modelling or numerical error. ⭐ **The audit that finds them**:
  ask each consumer what population it is conditioning on, and check the producer estimates that one.
* ⛔ **ALWAYS RUN BOTH fl-GAP SIGN ARMS AND THE EQUAL-LENGTH CONTROL.** The same change has been a win on
  `rna_long` and a loss on `gdna_long` at both the calibration and the transcript level.
* ⛔ **`np.diff(payload.region_bounds)` is WRONG** — bounds are concatenated per reference, so it
  manufactures a phantom region at every junction. Use `gdna_density.region_lengths_from_partition`.
  ⭐ No shipped code had this bug (audited); a prototype did, and it made an estimator look 6× worse.
* ⛔ **`docs/dev/` is a sandbox and nothing may cite into it** — a permanent doc that does is a suite
  FAILURE (`tests/test_docs_boundary.py`), which caught exactly that twice this session.
