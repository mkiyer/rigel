# The ablation campaign — which solver terms are still load-bearing, and how they interact

**Measured 2026-07-27 at `e1a041e6`.** One term off at a time plus the 2×2 factorial cell, on all 32
conditions at **both** refits, scored on accuracy (suite-wide and on the fit substrate) *and* honest
precision (`z2` on a node set **held fixed** to the baseline's confident quartile — never a per-arm
self-quartile).

Harness: `scratchpad/ablate_campaign.py` (runs the arms, caches the per-condition oracle so all arms share
one load) + `scratchpad/ablate_report.py`. The driver's `base` arm reproduces the benchmark's
`head0_r0` pooled mwae to **0.084524 vs 0.084524** — exact — which is what makes the rest trustworthy.

---

## 1. The table

`Δ` is against `base`; `z2` columns are the held-fixed populations. Lower mwae is better; `z2` **closer to
1** is better (it is MSE ÷ declared variance, so a large value means over-confident).

### refit = 0

| arm | mwae suite | Δ | mwae substrate | Δ | z2 ALL | exon single | exon AMBIG | bnd single | bnd AMBIG | per-condition |
|---|---|---|---|---|---|---|---|---|---|---|
| base | 0.0845 | — | 0.0852 | — | 15.00 | 1.18 | 450.5 | 20.28 | 111.3 | — |
| noP1e | 0.0878 | +0.0033 | 0.0898 | +0.0046 | 14.19 | 1.19 | 180.9 | 20.29 | **4.49** | 4 better / 7 worse / 21 flat |
| noP1d | 0.0841 | **−0.0004** | 0.0847 | **−0.0005** | 14.60 | 1.33 | 188.2 | 20.58 | 97.7 | 8 better / 5 worse / 19 flat |
| noP1e_noP1d | 0.0862 | +0.0017 | 0.0870 | +0.0018 | 13.82 | 1.35 | 114.3 | 20.61 | **4.44** | 3 better / 13 worse / 16 flat |
| noM5c | 0.0864 | +0.0019 | 0.0867 | +0.0015 | 13.72 | 1.39 | 89.7 | 22.20 | 53.0 | 4 better / 4 worse / 24 flat |
| **nofuse** | 0.0846 | **+0.0000** | 0.0853 | **+0.0001** | 14.50 | 1.18 | 296.1 | 20.10 | 107.9 | **0 better / 0 worse / 32 flat** |

### refit = 1 (production)

| arm | mwae suite | Δ | mwae substrate | Δ | z2 ALL | exon single | exon AMBIG | bnd single | bnd AMBIG | per-condition |
|---|---|---|---|---|---|---|---|---|---|---|
| base | 0.0680 | — | 0.0684 | — | 12.70 | 1.19 | 318.5 | 20.50 | 36.09 | — |
| noP1e | 0.0701 | +0.0022 | 0.0710 | +0.0026 | 12.63 | 1.20 | 200.5 | 20.53 | **5.34** | 5 better / 7 worse / 20 flat |
| noP1d | 0.0678 | **−0.0002** | 0.0682 | **−0.0002** | 12.18 | 1.34 | 115.0 | 20.96 | 33.66 | 5 better / 6 worse / 21 flat |
| noP1e_noP1d | 0.0692 | +0.0012 | 0.0688 | +0.0004 | 12.06 | 1.36 | 90.4 | 20.98 | **5.22** | 5 better / 15 worse / 12 flat |
| noM5c | 0.0709 | +0.0030 | 0.0713 | +0.0029 | 11.59 | 1.43 | 66.9 | 22.17 | 9.31 | 1 better / 5 worse / 26 flat |
| **nofuse** | 0.0680 | **+0.0000** | 0.0684 | **+0.0000** | 12.63 | 1.19 | 255.5 | 20.32 | 34.43 | **0 better / 0 worse / 32 flat** |

## 2. ⭐ THE INTERACTION — P1e and P1d are REDUNDANT, on all four measurements

`interaction = both − A − B + base`; negative means removing both costs **less** than the sum of removing
each, i.e. they were partly paying for the same failure.

| | cost(P1e off) | cost(P1d off) | cost(both off) | interaction |
|---|---|---|---|---|
| refit 0, suite | +0.0033 | −0.0004 | +0.0017 | **−0.00121** |
| refit 0, fit substrate | +0.0046 | −0.0005 | +0.0018 | **−0.00223** |
| refit 1, suite | +0.0022 | −0.0002 | +0.0012 | **−0.00078** |
| refit 1, fit substrate | +0.0026 | −0.0002 | +0.0004 | **−0.00197** |

Consistently negative, and largest on the fit substrate. **The owner's suspicion is confirmed**: the two
banner-flagged debts are double-paying. Both price a failure of the graft/peel edge — P1d the graft's
premise (a seam not speaking for its exon), P1e the resulting over-claim against the node's observed mass —
and P1e was landed *after* P1d without re-measuring it.

At production refit the combined removal costs **+0.0004 on the fit substrate** — essentially free — while
`z2` goes ALL 12.70 → 12.06, exon AMBIG 318 → 90 (3.5×), boundary AMBIG 36 → 5.2 (**7×**). But it is
**5 better / 15 worse** per condition, so it is a genuine trade, not a free win.

## 3. ⚠ What `z2` is actually measuring here — decomposed, because the ratio misleads

At AMBIG nodes the **error barely moves**; essentially all of the `z2` movement is the declared variance:

| population (refit 1) | arm | mass-wtd MSE | mass-wtd var | z2 |
|---|---|---|---|---|
| exon AMBIG | base | 3.081e−2 | 9.67e−5 | 318.5 |
| | noP1e_noP1d | 3.119e−2 | 3.45e−4 | 90.4 |
| boundary AMBIG | base | 1.517e−3 | 4.20e−5 | 36.1 |
| | noP1e | 1.517e−3 | 2.84e−4 | 5.34 |
| exon single | base | 3.110e−4 | 2.614e−4 | 1.19 |
| | noP1d | **3.460e−4** | 2.574e−4 | 1.34 |

So: **the AMBIG story is entirely about honesty of the stated precision, not about getting closer to the
truth.** The base solver declares a variance 3–7× too narrow there and the debts are what narrow it.

And the exon-single row is the counter-evidence on P1d: its MSE rises **11 %** (3.11e−4 → 3.46e−4) when P1d
is removed, at flat variance. P1d *is* doing real accuracy work on exactly the population it was fitted on.
**P1d is therefore not a clean delete** — it is a trade of exon-single accuracy against exon-AMBIG
over-confidence.

## 4. Verdicts — ⚠ SUPERSEDED FOR P1d/P1e BY §8

These were reached from the 10 Mb synthetic alone. **§8's real-data firing diagnostic overturns the
P1d and P1e rows** — read it before acting on this table. The level-fuse and M5 rows still stand.

| term | verdict |
|---|---|
| **the vestigial LEVEL FUSE** | ⭐ **DELETE.** Accuracy **0 better / 0 worse / 32 flat at both refits**, `z2` better on every population, and it removes ~15 lines that fuse a single estimator. See §5. |
| **P1d `ω_graft`** | **KEEP FOR NOW, re-derive at P1g.** Near-neutral suite-wide (−0.0004/−0.0002) but load-bearing on exon-single (its home population, +11 % MSE without it). The standing debt note already says it must be re-derived per structural class once TSS/TES land; that is still the right plan. |
| **P1e** | **KEEP, but it is the term to attack.** Genuinely load-bearing for accuracy (+0.0033/+0.0022 to remove) yet it is what makes boundary-AMBIG 7× over-confident. Exactly the "prices a bias as a variance" debt, now measured. |
| **M5's composition half** | **LOAD-BEARING.** Removing it costs +0.0019/+0.0030 — the largest accuracy loss in the campaign. Note it is verified **identically zero at every AMBIG node** (the freeze reference is `f_g = 1` ⇒ `f_g(1−f_g) = 0`), so all of that value is earned on single-strand nodes. |

## 5. The one clean deletion — the level fuse

`_peel_share` still runs three-estimator fusing machinery after two of the three estimators (`_far` and the
destination-own plug-in) were deleted. What survives:

* `_pt = _pm` — a bare alias;
* `_nu = (_pm·_nu_ms)/_pm` — an identity (the scalar twin annotates it as one);
* `_v_nu = ζ(2, ρ̂²/V̂)` — which **re-derives a log-variance `residual_level` already returned and
  `_peel_share` discards** (the `_` in `_nu_m, _, _vl_m = residual_level(...)`).

Algebraically `ρ̂²/V̂ = 1/v_log` exactly, so the surviving code computes `ζ(2, 1/v_log)` where `v_log` was
already in hand. That round trip is the identity only as `v_log → 0`; on thin seams it inflates, with a
ceiling of `ζ(2,1) = π²/6 = 1.64×`.

The `nofuse` arm replaces it with the discarded value. Measured: **accuracy identical (0/0/32 flat, both
refits)**, `z2` ALL 15.00 → 14.50 and 12.70 → 12.63, exon AMBIG 450 → 296 and 318 → 256, boundary single
20.28 → 20.10. Free, more honest, and less code.

## 6. Not yet run

`nofuse × noP1e` and `nofuse × noP1d` (does the fuse fix change the debts' verdicts?); M7's DL deflation
(τ-stream alone, then all sites); the variance-freeze reference; the duplicate λ-emission gate at the
combine; the Schur AMBIG τ gate. Ranked reasoning in the load-bearing inventory.

## 7. Method notes for whoever runs the next round

* The driver monkeypatches module-level terms (`graft_premise_logvar`, `composition_logvar`) and takes a
  hand-applied source patch for inline ones (P1e, the fuse). Worktrees do **not** work: the editable
  install's meta-path finder wins over `PYTHONPATH`, so a worktree silently runs the main repo's code.
* Always confirm `base` reproduces the benchmark before reading any arm.
* Freeze the `z2` node set on the BASELINE arm and let each arm bring its own error and variance. And
  **decompose the ratio** — §3 shows the AMBIG conclusion inverts if you read `z2` alone.

---

# 8. ⭐⭐ THE REAL-DATA FIRING DIAGNOSTIC — and it inverts the toy

**Measured 2026-07-27** on the cached real cfRNA samples (human index v7, ~1.5 M nodes, ~3 M messages),
where neither term had ever been characterised. `scratchpad/p1de_firing.py`; the P1e replay reproduces the
live damped precisions **exactly** (max rel err 0.000e+00), so these are the solver's own numbers.

## 8.1 P1e does 26–300× less work on real data than on the toy

| | toy (10 Mb) | cfRNA LBX0190 | cfRNA vcap |
|---|---|---|---|
| firing rate (of valid & supplied messages) | **16.5 %** | **0.09 %** | 1.13 % |
| **share of total precision it removes** | **31.1 %** | **0.10 %** | 1.18 % |
| of `δ < 0`, share clearing the noise floor | 36.6 % | 5.1 % | 3.6 % |
| median `den = αᵀΣα + 1/n_dst` | 0.080 | 3.16 | — |

The mechanism is not subtle: `b̂² = max(0, δ² − αᵀΣα − 1/n_dst)`, and on a **sparse** library `1/n_dst` is
large, so the DL noise floor swallows `δ²` almost everywhere. The toy is dense by construction.

**Consequence: P1e's measured toy benefit (+0.0022 / +0.0033 mwae to remove) will not transfer.** Every A/B
in the ledger scored with P1e active was scored on a solver whose P1e does 26–300× more work than it will in
production.

And where it does fire it is **more** bias than on the toy — bias share `E[δ]²/E[δ²]` = **0.692** (toy
0.438). The standing debt is worse on real data, not better.

⚠ It also fires on a **different population**. Unconditionally, the sign of `δ` REVERSES between toy and
real on both edge types — graft `E[δ]` −0.90 → **+0.27**, peel +0.78 → **−0.85** — so the `δ < 0` scope
gate selects grafts on the toy (62 % of the damping) and **peels on real data (83.5 %)**. The scope
argument itself is theoretical (non-negativity ⇒ only the over-claim attributes) and still holds; but its
empirical validation was done on a population it does not target in reality.

## 8.2 P1d's fitted scalar is not a library constant — it is an extreme-tail estimate

| | toy | LBX0190 | vcap |
|---|---|---|---|
| fitted `ω` (pos / neg) | ~0.30 | **0.25 / 0.40** | **2.80 / 3.41** |
| graft edges with a live seam pair | ~48 % | 12–17 % | 28–34 % |
| per-pair truncation (`d² ≤ noise`, contributes 0) | — | **95.5–96.2 %** | 73.6–73.7 % |
| pairs with `d == 0` exactly | — | **92.1–93.9 %** | 37.5–42.4 % |
| top-decile share of `Σd²` | 78–92 % | **100.0 %** | 98.2 % |
| effective n behind `E[d²]` (IPR) | — | **767 of 51,406 = 1.5 %** | 4.9–6.1 % |

**`ω` swings 10× between two real samples of the same assay.** On the sparse sample 92–94 % of seam pairs
report `d = 0` exactly and 96 % truncate to zero, so the library-wide scalar is fitted from ~1.5 % of the
pairs, with 70–76 % of `Σd²` coming from the top **1 %**. That is not "one library-level number" — it is a
handful of outlier pairs, and it multiplies the damping on **every** graft edge.

## 8.3 The redundancy mechanism — CONFIRMED, and wildly sample-dependent

Undoing P1d's damping on the same captured messages (P1e's `αᵀΣα` is built from post-P1d precisions, so
P1d structurally suppresses P1e):

| | firing | precision removed |
|---|---|---|
| toy | 946 → 959 | 6.19e4 → 6.28e4 (+1.5 %) |
| LBX0190 | 272 → 272 | 119.3 → 119.3 (**no effect**) |
| **vcap** | 12,334 → **14,014** | 1.52e4 → **4.95e4 (3.3×)** |
| vcap, graft edges only | 408 → **2,088** | 241 → **34,510 (143×)** |

Confirmed: P1d suppresses P1e, entirely on graft edges (non-graft is unchanged to the digit in all three).
But the magnitude ranges from **zero to 143×** across samples — so the toy's measured interaction
(−0.0008 to −0.0022) is not a transferable number either.

## 8.4 Verdict: retire **P1d**, keep P1e — with an honest caveat

On the owner's stated criteria:

* **theoretically sound** — P1e's trigger is a contradiction of a *hard observable* (`M`, the node's own
  count) with **no fitted parameter**, so there is nothing to mis-fit on new data. P1d's estimator is sound
  in *form* (MoM on a measured independent replicate) but its *population* is now measured to be a 1.5–6 %
  tail with 10× cross-sample instability.
* **cheapest / simplest** — P1e is ~30 lines inline, no fit, no extra pass. P1d is ~90 lines
  (`_seam_pair`, `_flank_dom`, both-seam frame lifting, `graft_premise_logvar`) plus a re-fit 3× per sweep.
* **failure direction** — P1e only ever widens, so it fails toward under-confidence (safe, correctable).
  P1d's pooled scalar on a bimodal population necessarily under-charges the high class (termini: 20.8 % of
  edges carrying 71.7 % of the error), and under-charging leaves messages over-confident — the expensive
  direction, and the one the whole project is fighting.

⚠ **The caveat that must travel with this decision:** removing P1d is not "keeping the effective term".
P1e is nearly inert on real data (0.1–1.2 % of damping mass). Retiring P1d means accepting that the graft's
premise error is **priced by essentially nothing** until **P1g** puts TSS/TES in the region map and `ω`
can be re-derived per structural class — which is where the ≥30× split becomes representable and the
quantity becomes principled rather than a stand-in. That is the real fix, and this measurement raises its
priority rather than lowering it.

⚠ Counter-evidence, kept in view: on the toy, removing P1d raises exon-single MSE 11 % on the confident
quartile — it *is* doing work where it was fitted. The case above is that the work is done by blanket
over-damping of the majority class, at a scalar that does not survive contact with a second real sample.
