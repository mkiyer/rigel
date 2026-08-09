# Fragment-length discrepancy — theory, audit, and what to do

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Opened 2026-08-08, from the owner's question: *"In real data the DNA fragment lengths can be much,
    much different from the RNA fragment lengths. We don't want to bias our calibration because of large
    differences in the fragment length distribution. (1) Theoretically, what is the right way to handle
    it? (2) What is our code doing now? (3) What should it do?"*

    Produced by an 11-agent workflow (4 derivations, 5 code audits, 1 synthesis, 1 adversarial review).
    Full agent output: the run transcript under `subagents/workflows/wf_3ee91d9d-b9a/journal.jsonl`.
    ⛔ **The reviewer's verdict was `sound-with-corrections`, and the corrections change the answer.**
    What follows is the CORRECTED position, with the four claims the owner's session verified by hand
    marked ✅.

---

## 0. ⭐⭐⭐ THE ONE THING TO READ

**The obvious fix is net-negative, and shipping it now would make the tool worse.**

The proposed repair is to replace the single pooled per-line share with a per-component share derived
analytically from the fitted pmfs. On the flgap panels that is a large win. ⛔ **On the ladder it is
~3.5× WORSE than the defect it removes** — shipped pooled-share error on the per-line g:r ratio
**+1.17 %**, the analytic replacement fed the *fitted* pmfs **−4.08 %** — because the fitted pmfs carry
a phantom, fixed-sign "gDNA is longer" artifact of **+5 to +18 bp** from the pool censoring in §2.

⭐ **So the direction is right and the SEQUENCE is wrong: fixing the pmfs is a PREREQUISITE, not a
successor.** Today's share is *empirical* (`mass / count` straight off the accumulator) and is therefore
immune to every pmf defect below. Replacing it with a model imports those defects into the split.

---

## 1. Theory — MOVED

⭐ The derivations that were here are now permanent, because the code depends on them:

* **`EQUATIONS.md` §3b** — the conserved mass's opportunity `A_mass(w;a,b) = [min(w−1,a)+min(w−1,b)]/2`,
  and **the pooling theorem**: one share for two components conserves the total exactly while biasing
  the ratio by `share_r/share_g`, so the error is purely compositional and no conservation gate can see it.
* **`EQUATIONS.md` §3c** — `count/inv_length_sum + 1 = E[w]`, the model-free local mean length, its
  natural-placement precondition, and why the tool may not be made gap-robust by shrinking the gap.

⚠ **The one theory claim that did NOT survive review**: the mechanism was described as "structurally
asymmetric". It is **symmetric** — `share_r/share_g` is bounded both ways and the long/short split is
1.07×. The observed ~4× flgap asymmetry is therefore NOT the pooled share.

## 2. What the code does today — the defects that matter

| # | site | defect | verified |
|---|---|---|---|
| 1 | `priors.py` / `substrate.py` | ONE pooled share rescales both components | measured |
| 2 | `fl.py` pool censoring | the fitted `mu_g − mu_r` carries a phantom **+4.7 to +14.9 bp** gap, fixed sign, **8/8** flgap conditions | ✅ **measured 2026-08-08** — and it SPLITS in two, see below |
| 2a | ~~`junction_opportunity`~~ **the DRAIN** | ⛔ **MOSTLY NOT A DIVISOR DEFECT.** Undrained the RNA pool reads −5.0 bp off capture with per-bin `fit/true` 1.049 / 0.907 across w = 200. **Drained it is −0.5 bp and 1.006 / 0.995.** The held fragments are systematically the LONG ones, so an undrained RNA pool is missing its own tail — the second pass already handles it | ✅ **re-measured drained 2026-08-08** |
| 2b | `gdna_opportunity` | ⭐⭐ **THE pmf DEFECT, and now the whole of it.** The four-pool de-tilt is **EXACT off capture** (+0.1 / +0.0 bp) and over-corrects **only under capture**: **+13.6** bp at a 330 bp gDNA mean, **+3.5** at 120, with drained per-bin `fit/true` running **1.22 … 4.18** in the tail. ⛔ **Untouched by the drain (Δ ≤ 0.1 bp)** | ✅ this is `EQUATIONS.md` §4.4 — a placement model, not a divisor |
| 3 | `priors.py` eff-len clamp | ⚠ **the "~3×" here is STALE.** Phase 0 measured the clamp at **1.03–1.25×**, with `nodes only` at 0.70–0.84× — so the edge term is partly compensating a real deficit, not purely inflating | ✅ superseded |
| 4 | `fl.py:86` | ⛔ **NOT a live defect on this substrate.** `prior_ess` 1000 → 0 moves every fitted mean by **≤ 0.1 bp**: pool totals are 0.70–4.75 M, so the shrink weight is ≤ 0.14 %. ⚠ Untested on sparse real data, where it is a different claim | ✅ ablated |
| 5 | every cached payload | `drain: null` while production always drains; the bound on `mu_r` is **[−3.96 %, −0.60 %]**, comparable to the phantom gap itself | agent |
| 6 | `prior_vs_oracle.py:118-125` | `edge_mass_per_crossing` is **excluded from `OVERRIDE_FIELDS`**, so the O arm inherits the pooled-share bias, `P−O` cannot see it, and the defect is mis-attributed into `O−F` | ✅ (reasoning) |

⭐ **A positive datum the synthesis got wrong and the reviewer corrected:**
`edge_unspliced_inv_length_sum` is **LIVE** (`second_pass.py:461` via `pipeline.py:893`) — and being the
one length-free channel, `second_pass`'s rho term is **the only provably gap-robust density estimator in
the tree.** ✅

⚠ **And the panels do not test what we say they test.** `flgap_long` is RNA 206/98 vs gDNA 330/120;
`flgap_short` is RNA 206/98 vs gDNA 120/50 — the **standard deviation differs by 1.2× and 2.0×** as well
as the mean. ✅ `share_c` is a censored functional and so is strongly variance-sensitive *even at equal
means*. Every analysis here, and the "±40 % gap" shorthand, is **mean-only** and therefore incomplete.

✅ **Debt CLEARED 2026-08-08** (`_density_times_span` deleted, the docstring re-based, plus two more it
turned up). The record is `accumulator_prior_plan.md` §5.5 and it is not copied here.

---

## 3. What to do — ⛔ MOVED

⭐ **The ordered plan now lives in `accumulator_prior_plan.md` §3**, re-ranked 2026-08-08 after the
assembler was measured to be a minor term. Two rows of the table that used to be here were retired by
measurement and the reasons are worth keeping:

* ~~**C** "fix the pool de-tilt"~~ — ⛔ **it was ALREADY divide-by-a-probability** (`fl.py:288-294`,
  `EQUATIONS.md` §4.1). Following that row would have had someone rewrite working code. What survives
  of the pmf defect is **2b alone**, and it is a capture-PLACEMENT problem.
* ~~**E** "drain one condition"~~ — ✅ **done, and it was not a footnote**: worth **+4.5 bp on `μ_r`**
  and ≤ 0.1 bp on `μ_g`, removing ~90 % of what had been blamed on 2a.

⚠ **Still open and unresolved:** the predicted pooled-share contribution under capture (−2.49 % on
flgap_short) **exceeds the owner's measured total** for that panel (`rel` 0.004–0.008). Either the
uniform-placement re-weighting is wrong under capture, or the measured figure is capture-OFF only. ⛔
**Resolve before quoting either number again** — a component larger than the total means something is
cancelling (`TRAPS: a-cancelling-defect-pair`).

✅ **Resolved:** the mate-gap junction question — a junction inside an unsequenced gap is HELD in the
side buffer, not filed as contiguous. `accumulator_prior_plan.md` §5.4.

---

## 4. The strategy — MOVED

⭐ The sequencing, the "does the accumulator change again?" answer and the ranked next steps now live in
`accumulator_prior_plan.md` §3–§4, which is the resume point. This file is the **theory and audit**
reference behind it; it should not grow a plan again.
