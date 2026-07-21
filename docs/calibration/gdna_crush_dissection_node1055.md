# The node-1055 f_g "crush" — mechanistic dissection (it is **not** a solver bug)

**Status:** diagnosis complete, hard evidence. Reproduce with `scripts/debug/crush_dissection.py`
(cached `ambig_dense_10mb`, `OMP_NUM_THREADS=1`, no full pipeline).

**Reviewer orientation.** An unstranded, capture-enriched exon (node 1055) with a *true* gDNA fraction
`f_g = 0.902` is deconvolved by the Phase-2 refit to `f_g = 0.0001`. The question we answer here, exactly as
posed: *the node's density sits at the enriched level; if the population prior had a mode there, why does the
node collapse to 0 instead of snapping to that mode?* The answer, proven end-to-end below, is that **the
solver is correct** — hand it a prior with an enriched mode and it lands the node at `f_g = 0.863`. The
collapse is entirely that the **fitted** prior carries no mass at the enriched level, for an *architectural*
reason (circularity), not a coding one. No amount of solver debugging fixes it; the training signal does.

---

## 1. The node and the crush (numbers)

`gdna_gdna100_ss_0.50_nrna_none_capture_on`, node 1055 — a captured single-strand exon in an **unstranded**
library:

| quantity | value |
|---|---|
| mass `M`, eff-length `E` | 24 235 , 2 011 |
| total density `ρ_tot = M/E` | 12.05 (log10 **+1.08**) |
| ORACLE gDNA `G`, RNA `R` | 21 860 , 2 375 |
| **true `f_g`** | **0.902** (true ρ_g log10 **+1.04**) |
| pass-0 `f_g` (no hyperprior) | 0.375 (var 3.56) |
| **refit `f_g` (with hyperprior)** | **0.0001** ← the crush |

Note the refit makes it **worse than pass-0** (0.375 → 0.0001). That detail is diagnostic — see §4c.

---

## 2. The exact code path (5 hops)

The composition (gDNA) arm of the per-node objective ψ is the fitted density projected by a **δ-pin**:
`f_g` enters only through `ρ_g = f_g·M/E`, so a prior over gDNA *density* becomes a prior over `f_g` at this
node's `ρ_tot`.

**Hop 1 — Fit** the population density `P(log ρ_g)` on the pass-0 deconvolved gDNA
(`calibration/calibrate.py:108`, `_fit_gdna_hyperprior`):

```python
g_hat = np.asarray(belief.f_g, dtype=np.float64) * mass_global      # <-- trained on pass-0's OWN f_g
...
return DensityNPMLE.fit(g_hat[sel], eff_global[sel], var_g=var_g, background=background, ...)
```

**Hop 2 — Project (the δ-pin)** onto the f_g grid (`calibration/npmle.py:308`, `DensityNPMLE.logprior`):

```python
log_me = np.log(mass) - np.log(eff)              # (n,) = log(M/E)
log_rho_g = np.log(fg)[None, :] + log_me[:, None]  # (n,K)  ρ_g = f_g·M/E
return np.interp(log_rho_g.ravel(), self.log_rho, self.logP,
                 left=self.logP[0], right=self.logP[-1]).reshape(log_rho_g.shape)
```

**Hop 3 — Build** `global_logprior` (m,K) for the sweep (`calibration/bp_solver.py:327`):

```python
gdna_prior.logprior(solve_grid, mass_global, eff_global) if gdna_prior is not None else None
```

**Hop 4 — Apply.** The prior **REPLACES the gDNA Jeffreys reference** in ψ
(`calibration/simplex_logodds.py:155` `_gdna_arm`, added at `:268`):

```python
def _gdna_arm(lam, global_logprior):
    if global_logprior is not None:
        return np.asarray(global_logprior, np.float64)   # the fitted prior REPLACES ...
    return _JEFFREYS_REF * _log_fg(lam)[None, :]          # ... this ½·log f_g reference
...
psi = psi + _gdna_arm(lam, global_logprior) + _rna_arm(lam)
```

**Hop 5 — Read out** `f_g` as the posterior median (`calibration/simplex_logodds.py:357`, `:127`):

```python
f_g = _posterior_median_fg(post, fg)     # fg at the grid point where the CDF first reaches 0.5
```

Each hop was a suspect. §3 clears all five.

---

## 3. The decisive experiment — swap **only** the prior shape

Hold node 1055's real data fixed; change nothing but the prior's shape; watch `f_g`.
(`scripts/debug/crush_dissection.py`, Parts B & C.)

**Part B — isolated δ-pin** (flat/unstranded strand ⇒ ψ = `logP_g(node) + RNA_ref`):

| prior shape | median | mean | argmax | logP_g at f_g = [.001 .01 .1 .375 .5 .9] |
|---|---|---|---|---|
| **(a) REAL fitted prior** | **0.000** | 0.019 | 0.000 | [ −3.4 −9.0 −5.8 −8.0 −7.9 **−7.3** ] |
| (b) synth DEPLETED-only @ −2.0 | 0.001 | 0.001 | 0.001 | [ −3.5 −29 −69 −69 −69 −69 ] |
| **(c) synth ENRICHED @ truth +1.04** | **0.863** | 0.828 | 0.757 | [ −69 −69 −24 −6.6 −4.8 **−3.4** ] |
| (d) synth BIMODAL (−2.0 & +1.04) | 0.695 | 0.522 | 0.001 | [ −4.2 −30 −24 −7.3 −5.5 −4.1 ] |

**Part C — full real refit, prior replaced by (c):** node 1055 → **f_g = 0.8633** (vs 0.0001 with the real
prior). End-to-end, through the real sweep, messages and all.

**Reading the table.** Under (a) the real prior, `logP_g` at the truth (`f_g=0.9` ⇒ ρ_g=+1.03) is **−7.3**,
far below its value at `f_g=0.001` (−3.4) — the prior *prefers* the depleted answer, so the node goes there.
Flip the mode to the enriched level (c) and the *same* node, *same* solver, *same* median readout lands at
**0.863**. The machinery is correct. (Row (d) also shows the median is well-behaved — 0.695 — when enriched
mass exists; the readout is not the problem either.)

**Conclusion: there is no bug in hops 2–5.** The crush is hop-1 — the fitted prior has ~zero mass at the
enriched level (`logP_g ≈ −7` there), so the δ-pin correctly reports "the population says this density is
depleted."

---

## 4. Root cause — why the fitted prior has no enriched mode

The real prior's modes sit at log10 ρ_g ∈ **[−2.96, −1.83]** (all depleted); its grid *reaches* +1.80 but
carries negligible mass above ~0. Three compounding reasons, none of them a solver defect:

### 4a. Circularity — the prior is trained on pass-0's own under-call
Hop 1 fits on `g_hat = belief.f_g · M` — **pass-0's deconvolved gDNA**. Pass-0 under-calls every enriched
*unstranded* node (node 1055 pass-0 = 0.375, not 0.902), because unstranded + capture is locally
non-identifiable (no strand tilt to reveal gDNA). So node 1055's own training point lands at
ρ_g = 0.375·12.05 ⇒ log10 **+0.65**, not +1.04. **The prior cannot carry mass at a density pass-0 never
produced.** At best it reproduces the under-call; it can never reach the truth. This is the core loop: a
prior built from the answer cannot correct the answer.

### 4b. Minority suppression — the enriched nodes are outvoted
Even the +0.65 training points are a tiny minority (a handful of captured exons) against thousands of
depleted intron/exon nodes. The default fit is the **competitive EM NPMLE** (`gdna_prior_additive=False`),
whose mixture-EM (`npmle.py:113` `_em_weights`) actively competes the minority mode toward zero weight — so
the fitted modes collapse to −2…−3 and `logP_g(+0.65) = −8` (deep tail, per the table row (a)). An *additive*
KDE (each node one un-competed unit kernel) preserves more of the +0.65 mass, but see §5 — it still cannot
exceed +0.65.

### 4c. The Jeffreys replacement — why the refit is *worse* than pass-0
Hop 4: when a fitted prior is present, `_gdna_arm` **returns it instead of** the `½·log f_g` Jeffreys
reference. That reference is a soft `f_g → 0` barrier (½·log f_g → −∞ as f_g → 0) that kept pass-0 off the
vertex (pass-0 = 0.375). The refit deletes that barrier *and* substitutes a term (the depleted prior) that
actively pulls toward f_g → 0. Nothing else in ψ bounds `f_g → 0` (the RNA arm only bounds `f_g → 1`). So a
depleted-only prior is **strictly worse than no prior at all** at an enriched node — exactly the 0.375 →
0.0001 regression.

---

## 5. Why "go back to the KDE" does not fix node 1055

The released v0.7.1 KDE and this branch's NPMLE are both density estimators of the **same under-called pass-0
gDNA**. The additive KDE mitigates §4b (it does not competitively silence the minority) — a real, worthwhile
improvement — but it cannot escape §4a: its enriched training points are still at +0.65, so the best it can
do at node 1055 is pin ≈ pass-0's 0.37, not the true 0.9. The release "worked" everywhere *except*
`ss0.50 + capture_on`, which is precisely this unstranded-capture collapse (see the `ambig_dense` benchmark
note: "gDNA collapses ONLY unstranded+capture"). The KDE is a better *representation*; it is not a fix for
node 1055. Reverting buys robustness elsewhere, not this pole.

---

## 6. The only non-circular escape — condition the prior on enrichment

To place mass at +1.04 for node 1055, the prior needs a signal that is **not** the under-called composition.
There is exactly one: **enrichment** — the node's *total* density (log10 +1.08), which is directly observed
(count-based) and **not** under-called. The fix is to make the prior a distribution over the **intrinsic**
gDNA rate `d_g = ρ_g / e(x)` (enrichment `e(x)` from the Role-A total-density projection `DensityNPMLE.project`
→ `mu_proj`, already computed and already passed into the refit). Then at an enriched node ρ_g = d_g·e(x) is
re-centered *up* by the enrichment, and the δ-pin no longer crushes. This is not "another layer" — it
**replaces** the enrichment-blind marginal with the correctly-conditioned one, and it is the only construction
that breaks the §4a circularity. (Design: `gdna_hyperprior_clean_slate.md` §9.)

A weaker stopgap — restore the `½·log f_g` Jeffreys barrier *underneath* the fitted prior (so the prior can
add enriched mass but can never pull *below* the reference) — would at least stop §4c (no worse-than-pass-0),
turning the confident 0.0001 back into an honest 0.375. That is a mitigation, not a cure.

---

## 7. Reviewer's checklist — where a real bug *would* hide, and what to send

The dissection cleared hops 2–5, but if you want to independently confirm, these are the exact surfaces:

1. **The δ-pin units** (`npmle.py:308-340`). `log_rho` and `logP` are **natural**-log; `logprior` uses
   `np.log` (natural). A decade/ln10 mismatch here would mis-place every node. (Checked: consistent — the
   enriched-prior test pins correctly, which it could not if units were off.)
2. **The grid top / right-clamp** (`npmle.py:216-234`, `interp(..., right=logP[-1])`). If the grid failed to
   span `ρ_g = f_g·M/E`, the enriched f_g would clamp to a tail value. (Checked: grid reaches +1.80 > the
   node's +1.08 accessible max; not clamped.)
3. **The Jeffreys replacement** (`simplex_logodds.py:155-163`). Confirm the prior is meant to *replace* the
   reference; if so, §4c says a depleted prior removes the `f_g→0` barrier — arguably the reference should be
   retained as a floor. **This is the one place a design decision is genuinely debatable.**
4. **The training mask** (`calibrate.py:126-145`). Which nodes feed the fit, and at what `f_g` — this is where
   the circularity (§4a) lives.
5. **The readout** (`simplex_logodds.py:127`). Median vs mean — dissected separately
   (`readout_median_vs_mean_ab.py`); it is a rounding effect on an already-collapsed posterior, not the
   cause.

**To send a reviewer:** this document + `scripts/debug/crush_dissection.py` (self-contained, prints Parts
A/B/C above) + the five code snippets in §2. The one-line summary: *the solver correctly reports what the
prior says; the prior says "depleted" because it was trained on an under-call it cannot see past.*

---

## Appendix — reproduce

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
OMP_NUM_THREADS=1 python scripts/debug/crush_dissection.py
```
Part A prints the node facts + the real prior's modes; Part B the isolated δ-pin table; Part C the full-refit
confirmation (`ENRICHED-prior refit = 0.8633`).
