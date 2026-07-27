# P1d and P1e — what they are, why both exist, and when to remove them

**Status 2026-07-27: BOTH RETAINED, deliberately.** They are measurably redundant and both are stand-ins for
one missing structural fact. The decision is to **leave them in place until TSS/TES enter the region map
(P1g), then re-evaluate both together and expect to simplify.** This document is the brief for that
re-evaluation.

Evidence behind everything here: `ABLATION_CAMPAIGN.md` (the toy 2×2 factorial, §1–§7; the real-data firing
diagnostic, §8). Reproduce with `scratchpad/ablate_campaign.py`, `scratchpad/ablate_report.py`,
`scratchpad/p1de_firing.py`.

---

## 1. The one failure they both price

A boundary hands the exon next door its RNA density — the **graft**. What the boundary measured is the RNA
crossing *that seam*; what the exon needs is the RNA *in the exon*. Those differ whenever molecules reach
the exon by its other flank, or start or end inside it. So the graft knows an **inequality**

    ρ_R(exon)  ≥  ρ_ν(B) + ρ_μ(B)

and uses it as an **equality**. That is the whole failure. Both terms exist because of it, and they attack
it from opposite ends:

| | **P1d — `ω_graft`** | **P1e — the conservation surprise** |
|---|---|---|
| direction | the **UNDER**-claim (the seam speaks for less than the exon holds) | the **OVER**-claim (the delivered components assert more fragments than the node observed) |
| observable? | **No.** An under-claim is indistinguishable from "this exon has more gDNA" — count-zero-information exactly. It has no local signature. | **Yes.** `δ = log(M/S)` against `M`, the node's own count. A hard observable, no prior. |
| therefore estimated from | a **population**: an exon's two flanking seams as independent replicates, method of moments, pooled to one library scalar | **nothing** — it is computed per message from that message's own budget |
| where | `enrichment_frame.graft_premise_logvar`, fitted in `bp_solver._seam_pair`, applied on graft edges | inline in `bp_solver._transport`, scoped to `δ < 0` |
| code | ~90 lines (`_seam_pair`, `_flank_dom`, both-seam frame lifting) + a re-fit 3× per sweep | ~30 lines, no fit, no extra pass |

The complementarity is real and is why both were written. **P1e cannot see the under-claim** (nothing the
unsupplied components could be would reveal it — every density is non-negative, so a *shortfall* in the
budget is equally explained by the node's own density being low). **P1d cannot see the over-claim** (it is a
population constant; it does not look at the destination at all).

## 2. Why they are nevertheless redundant

**Two diagnoses, one treatment.** Both end as a variance added to the same message's precision on the same
edges, and precision damping composes: `p → p/(1+p·v₁)/(1+p'·v₂)`. Whichever fires first absorbs part of
what the other would have done — regardless of whether their *reasons* are disjoint.

And there is a structural coupling on top of that: **P1e's noise subtraction is built from precisions that
P1d has already damped.** `b̂²_cons = max(0, δ² − αᵀΣα − 1/n_dst)`, and `αᵀΣα` is assembled from
`_vc = 1/_p3` where `_p3 = [tpg, tpp, tpn]` are *post*-P1d. So a message P1d has widened has a larger noise
floor and truncates P1e to zero more often. **P1d structurally suppresses P1e.**

Measured (`ABLATION_CAMPAIGN.md` §2, §8.3):

* toy 2×2 factorial: `interaction = both − A − B + base` is **negative on all four measurements**
  (−0.0008 to −0.0022) — removing both costs less than the sum of removing each;
* the mechanism, by undoing P1d's damping on the same captured messages: P1e's firing and removal rise by
  0 % (LBX0190), +1.5 % (toy) and **3.3× on vcap — 143× on graft edges alone**, with non-graft edges
  unchanged to the digit. Exactly the predicted signature, and confined to grafts.

## 3. What each is measured to be worth — and why the toy number is not the real number

| | toy (10 Mb synthetic) | real cfRNA |
|---|---|---|
| **P1e** share of total precision it removes | **31.1 %** | **0.10 %** (LBX0190) / 1.18 % (vcap) |
| P1e firing rate (valid & supplied msgs) | 16.5 % | 0.09 % / 1.13 % |
| P1e bias share `E[δ]²/E[δ²]` on firing msgs | 0.438 | **0.692** |
| **P1d** fitted `ω` (pos / neg) | ~0.30 | **0.25 / 0.40** vs **2.80 / 3.41** — a **10× swing between two real samples** |
| P1d per-pair truncation (`d² ≤ noise`) | — | **95.5 %** (sparse) / 73.6 % |
| P1d pairs with `d == 0` exactly | — | **92.1 %** / 37.5 % |
| P1d effective n behind `E[d²]` (IPR) | — | **1.5 %** / 4.9 % of pairs |

Two things follow, and they are the reason this document exists:

1. **P1e is nearly inert in production.** The DL noise floor `1/n_dst` is large on a sparse library, so
   `δ²` rarely clears it — 3.6–5.1 % of `δ < 0` messages, against 36.6 % on the toy. Its toy-measured
   accuracy contribution (+0.0022 / +0.0033 mwae to remove) **will not transfer**. It also fires on a
   *different population*: the sign of `δ` reverses between toy and real on both edge types, so the scope
   gate selects grafts on the toy (62 % of damping) and **peels on real data (83.5 %)**.
2. **P1d's `ω` is not a library constant.** On the sparse sample 92 % of seam pairs report `d = 0` exactly
   and 96 % truncate below the noise floor, so the scalar is carried by ~1.5 % of pairs with the top 1 %
   supplying 70–76 % of `Σd²`. It is an extreme-tail estimate that multiplies the damping on **every**
   graft edge.

## 4. ⭐ The root cause both are standing in for

`ω_graft` was measured to split **≥30×** on whether the boundary carries a **transcript TERMINUS**
(`ω̂` 1.7–1.9) or is a **junction-only seam** (0.04–0.06) — with 20.8 % of edges carrying 71.7 % of the
error. That is not a nuisance; it is the *entire* effect. And it is exactly the distinction the region map
**cannot represent**, because it has no TSS/TES: the solver cannot tell a splice junction from a transcript
end.

So a single pooled scalar is being asked to stand for a bimodal quantity whose two modes differ by 30×. It
necessarily **over-charges the junction-only majority and under-charges the terminus minority** — and
under-charging is the expensive direction, because it leaves messages over-confident.

At a terminus, RNA does **not** flow through the seam at all. The graft's inequality is not slightly loose
there; its premise is simply false. That is a *structural* fact, and the honest fix is structural.

## 5. The roadmap

### Now — keep both. Do not delete on the current evidence.

* Removing P1d costs +11 % exon-single MSE on the confident quartile (its home population).
* Removing P1e costs +0.0022 / +0.0033 mwae on the toy, and the real-data measurement says that number is
  not trustworthy in either direction — it is measuring a term doing 26–300× more work than it will do in
  production.
* Removing both is a genuine trade, not a free win: +0.0004 on the fit substrate at production refit, but
  **5 better / 15 worse** across the 32 conditions.
* Neither deletion is justified while the structural gap they cover is still open.

### ▶ P1g — TSS/TES in the region map. **This is the fix, and this measurement raises its priority.**

Once a boundary knows whether it carries a transcript terminus:

1. **Re-derive `ω_graft` per structural class** — the same method-of-moments equation, one scalar per class.
   The partial-pooling block in `graft_premise_logvar` is the plug-in point. The ≥30× split becomes
   representable, so the estimator stops averaging two populations and the pooled form's whole failure mode
   disappears.
2. **Expect the terminus class to want a graft GATE, not a variance.** At a true terminus the premise is
   false rather than loose; a wide variance is the wrong instrument for "this message should not exist".
   If that lands, P1d's variance shrinks to the junction-only class, where `ω̂ ≈ 0.05` — small enough that
   it may not earn its ~90 lines and per-sweep fit at all.
3. **Then re-test P1e against the corrected graft.** P1e's `δ < 0` is largely *caused* by the graft
   over-claiming; fixing the graft structurally should shrink `δ`'s systematic part, which is precisely the
   "when the bias strata are diagnosed, this term must SHRINK" condition already recorded in
   `variance_ledger.md` §6 and CLAUDE.md. **If P1e does not shrink after P1g, the model has not improved.**
4. **Re-run the 2×2 factorial.** The interaction is the diagnostic: if P1d and P1e stop being redundant
   after P1g, each is pricing something real and separate. If they remain redundant, delete the one with
   the fitted parameter — P1d — because a term with nothing to fit cannot be mis-fitted on new data.

### The order of preference, if only one survives

**Keep P1e, retire P1d**, on three grounds that do not depend on the toy:

* **no fitted parameter** — `δ` is computed per node against that node's own observed mass, so there is
  nothing to mis-fit on a new library. P1d's scalar swings 10× between two samples of the same assay.
* **a third of the code**, no fit, no extra pass.
* **failure direction** — P1e only ever widens a message, so it fails toward under-confidence: safe, and
  correctable by later evidence, which is the standing requirement on pass-0. P1d's pooled scalar
  necessarily under-charges the terminus class, which leaves messages over-confident.

⚠ With the caveat that must travel with that decision: **P1e is nearly inert on real data**, so retiring
P1d is not "keeping the effective one". It means the graft's premise error is priced by almost nothing
until P1g lands. That is an argument for doing P1g *first*, not for keeping a scalar that does not
reproduce across samples.

## 6. What NOT to do

* **Do not delete either on a small suite delta.** Both are near-neutral suite-wide; the toy cannot
  distinguish "inert" from "mis-measured", and §3 shows it gets P1e wrong by 26–300×.
* **Do not rebuild P1d's per-edge form.** `d²` from one pair is a χ²₁ (CV = √2) — mostly noise — and it
  makes a message's precision depend on a non-adjacent node's counts, a real BP violation. Measured worse
  on every axis; the partial-pooling James–Stein variant was also implemented and deleted.
* **Do not rebuild P1e's rank-1 attribution.** MC-refuted: it over-damps λ 5× on a pure scale error.
* **Do not widen P1e's scope past `δ < 0`.** The unscoped half was measured inert (0.0841 vs 0.0841), and
  the shortfall direction does not attribute — `_pin_v`'s partial-claim semantics fill unsupplied
  components from the node's own density, so a shortfall is equally the node's own fault.
* **Do not treat `ω`'s fitted value as a measured constant of the library.** It is not. §3.
