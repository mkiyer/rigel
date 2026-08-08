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
| 2 | `fl.py` pool censoring | gDNA pools are geometry-limited, the RNA pool splice-limited — **opposite** censorings, so the fitted `mu_g − mu_r` carries a phantom **+5 to +18 bp** gap that is not in the library | agent |
| 3 | `priors.py:427,436,465` | ⭐⭐ **`gdna_eff_len` is clamped by an INCIDENCE-support sum, not the genomic span** — `span = node_eff + Σ edge_eff`, and each line adds `~mu_g − 1` to a ~159 bp node, so ~**3×** inflation. The EM divides the gDNA weight by this | ✅ |
| 4 | `fl.py:86` | `POOL_EB_PRIOR_ESS = 1000.0` — a constant sitting directly on the composition determinant | agent |
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

⚠ **Debt this session created**, to clear in the stage-C commit: `_density_times_span` is now **dead**
(no call sites ✅), and `assemble_priors`' docstring still teaches the retired
`rho_c = Σm/ΣS; prior = rho·span_bp` rule at `priors.py:293-297` ✅.

---

## 3. What to do — resequenced

⛔ **NOT rank 1.** The analytic per-component share must wait for the pmfs.

| | step | why here |
|---|---|---|
| **A** | ⭐ **Make the O arm an actual oracle on this axis** — add the per-line share to `OVERRIDE_FIELDS`, computed from the **origin-split payloads** (`parts['gdna'].edge_unspliced_mass / …count`), which ARE `share_g` on the real partition under real, non-uniform, post-capture placement | one tuple entry; no new simulation, no solver run, oracle cache stays valid. ⛔ Compute it **empirically from the split**, never from `truth_fragment_lengths.tsv` through `crossing_eff_length` — that would make the O arm a MODEL arm and defeat its purpose |
| **B** | Price defect 3 (`gdna_eff_len`'s incidence cap). A ~3× inflation of one of the three `LocusPriors` fields dwarfs a 0.3–2.5 % share bias | may outrank everything else |
| **C** | Fix the pmf estimator: pool de-tilt by **membership probability** (`TRAPS: divide-by-a-probability`), and the shared EB anchor | **prerequisite** for any analytic share |
| **D** | Only then, the per-component share — and evaluate the **hybrid** (empirical scale `M(e)`, analytic ratio `r`) against the fully-empirical route from step A | sequence set by §0 |
| **E** | Drain one condition and re-fit the FL models | until this is a number, every FL measurement here carries a caveat its own size |

⛔ **Gate everything on the flgap PAIR, both capture arms — never the ladder**, which scores `+0.029 %`
on this mechanism and is structurally blind to it.

⚠ **Open and unresolved:** the predicted pooled-share contribution under capture (−2.49 % on
flgap_short) **exceeds the owner's measured total** for that panel (`rel` 0.004–0.008). Either the
uniform-placement re-weighting is wrong under capture, or the measured figure is capture-OFF only. ⛔
**Resolve before quoting either number again** — a component larger than the total means something is
cancelling (`TRAPS: a-cancelling-defect-pair`).

⚠ **Unranked and unverified:** a junction inside an unsequenced mate gap may be filed as a *contiguous*
unspliced crossing, tilting the unspliced RNA population long exactly beside junctions. If true it is a
length-driven composition bias upstream of all of the above. **Needs checking in `build_fragment`.**

---

## 4. The strategy — MOVED

⭐ The sequencing, the "does the accumulator change again?" answer and the ranked next steps now live in
`accumulator_prior_plan.md` §3–§4, which is the resume point. This file is the **theory and audit**
reference behind it; it should not grow a plan again.
