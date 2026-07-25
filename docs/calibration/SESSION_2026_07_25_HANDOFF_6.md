# Session handoff — the message-variance model is DONE; Phase 2 (the gDNA hyperprior) is the blocker

**Read `docs/calibration/ROADMAP.md` first, then this. This is the LIVE handoff** — it supersedes
`SESSION_2026_07_24_HANDOFF_5.md` (whose §4 task is now complete; keep it and `..._HANDOFF_4.md` for the arc).
Do NOT read `docs/calibration/archive/`. Date: 2026-07-25.

---

> ## ⭐ SUPERSEDING UPDATE (2026-07-25, later the same day) — read this before §0
>
> A suite-wide dissection of where pass-0's error actually lives changed the ordering below. **The next task is
> NOT Phase 2 — it is `density_composition_reconciliation.md`** (owner-directed). Summary of what changed:
> * **92.2 %** of suite error mass is on `τ_own = 0` nodes that cannot self-solve, and messages *halve* their
>   error (0.333 → 0.164) — they are doing their job. The tractable bug set is 7,462 full-rank nodes the
>   messages made **worse** (707 k err-mass, self 0.021 → solved 0.090; 81.7 % exons, 92.4 % stranded × capON).
> * Root cause: **the composition regime is scale-blind.** Relayed claims assert more reads of a component than
>   the node sequenced in total — `Σ_c ρ_c·E_c / M > 1` on **52–71 %** of nodes (p99 31–288×, max 519×).
>   `_pin_v` enforces that bound and its own docstring says it belongs "at EVERY node"; it runs only in the
>   combine. Worth ≈ **2.2 M error mass (16.5 % of suite)**, a MODE defect, fixable **without** the hyperprior.
> * Landed since: the **λ-emission gate** (`b81e926e`) — a source with only one component has no composition
>   claim to make; 0 scenarios worse, refit=1 aggregate **0.0819**, unstranded × capON **0.1681**.
> * Tried and **reverted**: the graft-frame fix (measured mature must not be reframed — premise proven by
>   oracle, 82× over-claim at node 1909, aggregate 0.0964→0.0789). It regresses zero-gDNA libraries
>   (solved `f_g` 0.32 vs oracle 0.0000 over 72 % of mass) because the over-claim was **load-bearing**. Land it
>   only together with the density anchor. Full record: `density_composition_reconciliation.md` §5.
> * Two adjacent modelling gaps recorded there: **no TSS/TES in the region/boundary map** (§4.1, deferred,
>   reaches to the index/accumulator) and **the boundary is a slope, not a cliff — three enrichment ratios, not
>   one** (§4.2, probe placement; oracle-measured 0.97 → 0.87 → 0.54 as capture strengthens).
>
> §3 below (Phase 2) remains correct and is the task **after** that.

## 0. TL;DR — what to do next (one paragraph)

The message-variance model is **finished and winning**: the cross-cliff cost is now split into the two objects
it actually consists of — the reframe's SCALE noise (M5 `Var(log r)`) and the COMPOSITION-mismatch bias
(`b̂²`, DerSimonian–Laird against the node's own self-solve, no tuned constant) — and that replaced the
`(log r)²` proxy. Aggregate pass-0-vs-oracle is the best on record at both refit settings
(0.1267→**0.0969** refit=0, 0.1234→**0.0828** refit=1) with no per-condition regression against the arm it
replaces. **The next task is Phase 2: the gDNA hyperprior**, which is now the single binding constraint. It
blocks two things at once — the ship candidate, and the AMBIG fix. **Do them in this order: (1) fix the
refit's unstranded-capON regression (§3.2), (2) re-apply the measured 6-line hyperprior→`v_own` step (§3.1),
which then fixes AMBIG almost for free, (3) re-solve → ship candidate.** §3 has the exact experiment, its
numbers, and why it must be sequenced that way.

## 1. The owner's governing principle (unchanged, do not violate)

**Pass-0 must be WEAK and correctable.** An over-confident message that PINS a node wrong is worse than a weak
one that is slightly off, because high precision is nearly irreversible in the refit. Prefer the weaker
(more damped, under-confident) option when unsure. **No magic numbers** — pause and discuss before any new
constant or heuristic. Counts are Poisson (ω=0); the synthetic suite is Poisson by construction, so it cannot
validate anything overdispersion-dependent. `OMP_NUM_THREADS=1` for A/B determinism.

> The DL term gives this principle an exact form worth remembering, because it is the invariant to preserve:
> `p_eff = 1/max(v_msg, G² − v_own)` ⇒ **a message out-weighs a node's own belief only if it agrees to within
> `√2·σ_own`** — at every node, every composition, and independently of how deep the source's counts are. No
> amount of neighbour confidence can buy a pin. Any future change that breaks this inequality is a regression
> even if the A/B improves.

## 2. What is DONE + committed (branch `calib-ambig-init-wip`)

* `1a3e0a89` — retire the dead NPMLE σ²_transfer plumbing; struct_lock resolution + interior-anchor test.
* `f7c02c8e` — **the DerSimonian–Laird composition-mismatch cliff term `b̂²`** (replaces `(log r)²`).
* `813f18ba`/`e41c4e58`/`d78da565`/`44f1ecc6` — the arc before it (see HANDOFF_4/5).

### 2.1 The law, as it now stands (`message_variance_derivation.md` M1–M7)

```
    p_message = 1 / ( Var(log f_c^src) + 1/n_src  +  σ²_transfer  +  b̂² )
```

* **`σ²_transfer = Var(log r)`** (M5, `enrichment_frame.transfer_logvar`) — the reframe's own scale noise. 0 on
  the matched graft (`r` is common-mode there and cancels in the composition), load-bearing on peel/partial.
* **`b̂² = max(0, G² − v_msg − v_own)`** (M7, `enrichment_frame.mismatch_gap`/`mismatch_deflate`) — the
  composition-mismatch bias, `G = mo^msg − mo^own` measured against the node's `node_init` self-solve,
  `v_own` from the `τ_λ` foundation via `own_composition_logvar`. Applied in `_transport` on the POST-`_pin_v`
  densities (the pin removes the message's scale claim, leaving a pure disagreement about the split), to all
  three streams; the composition stream gets a single **λ-axis** gap (`G_λ = G_g − G_R`, `v_own,λ = 1/τ_own`)
  because it is one DOF — a per-component gap there is a category error that under-damps by 1–4×.

Settled facts (do NOT relitigate):
* The delivered mode error splits EXACTLY into share-mismatch + scale-noise, to machine precision (M7a/M7b).
* DL recovers the true `b²` to 0.1–0.4% for real mismatches; at `b≈0` it over-damps by
  `0.484·(v_msg+v_own)` — the safe direction, and harmless (a message that agrees moves the fused mode nowhere).
* **τ_own = 0 ⇒ b̂² = 0, bit-identically.** Every AMBIG node and ALL of unstranded data. This is BY DESIGN —
  it is what preserves the M5 unstranded/capture win — and it is why AMBIG is Phase 2's job, not pass-0's.
* A one-sided absence (a pure-gDNA seam asserting `f_c = 0`) is the `b̂² → ∞` limit ⇒ precision 0, carried as a
  mask so the numerical zero-test `_EPS` never sets the surviving precision. (Measured identical to 4 decimals
  against the `_EPS`-floored variant, so the `_EPS`-sensitivity concern raised in review is immaterial —
  but the mask form is the one with no constant in it.)
* **No prior of any kind now enters message precision.** The NPMLE-projection σ²_transfer is retired; the
  enrichment NPMLE is still fit in `calibrate` but only for the QC report + toy injection.

Verification: 4-agent workflow `wf_7d1708d4` (independent τ-form re-derivation, adversarial audit, code audit,
adjudicator) — all three converge on the landed design. Its per-agent reports are worth reading before
changing this term; the scripts are `scratchpad/dlwf_*.py`.

Gates: `pytest tests/calibration tests/native` = **392 pass, 2 xfail, 2 xpass**; `ruff check src/ tests/
scripts/` clean; `scripts/debug/message_variance_mc.py` = 0 failures over M1–M7. Goldens regenerated.

## 3. ▶ THE NEXT TASK — Phase 2, the gDNA hyperprior

### 3.1 The Phase-2 step is already written and measured — it is BLOCKED, not unknown

HANDOFF_5 §7.1 ("feed the hyperprior into DL's `v_own` on AMBIG nodes") was implemented and A/B'd this
session. It is ~6 lines in `bp_solver._unified_solve`, right where `v_own_g/v_own_r/v_own_lam` are built:

```python
    # the hyperprior's own λ-curvature, by the SAME law the intron factory uses — no new law, no constant.
    # `build_node_init` already folds the prior into the self-solve's MODE but not into τ_lam, so on a refit
    # the own belief was judged at a precision that ignored half of what formed it.
    tau_prior = density_factor_precision(global_lp, _lam_lo)      # None at pass-0 ⇒ identically inert
    tau_dl = tau_own if tau_prior is None else tau_own + np.asarray(tau_prior, np.float64)
    v_own_g, v_own_r = own_composition_logvar(_ni.f_g, tau_dl, _struct)
    v_own_lam = np.where(_struct, 0.0, np.where(tau_dl > _EPS, 1.0 / np.maximum(tau_dl, _EPS), np.inf))
```

Add it to the DL null ONLY — never to the SEND precision (`prec_*`): a node knowing more about itself is
grounds for trusting a contradicting message LESS, not for shouting louder. So it can only ever damp.

**Measured (refit=1, `pass0_oracle_bench.py`, arm `p2_r1`):**

| condition | HEAD | + Phase 2.1 | |
|---|---|---|---|
| stranded ss_0.99 | 0.0376 | **0.0333** | ✅ |
| verystrong | 0.1292 | **0.1196** | ✅ |
| capture OFF | 0.0354 | **0.0168** | ✅ |
| **unstranded × capON** | 0.1702 | **0.2177** | ⛔ |
| aggregate | 0.0828 | 0.0849 | ⬆ |

It does exactly what it was predicted to do on every arm the hyperprior is good on, and it regresses
unstranded-capON — because **where the hyperprior is wrong, DL now enforces it** against the very messages
that were correcting it. That is not a flaw in the coupling; it is the known refit defect (ROADMAP §4)
becoming load-bearing. **Reverted, deliberately.** Re-apply it as step 2, after §3.2.

### 3.2 So the real task is: fix the refit's unstranded-capON regression

The refit machinery exists (`calibrate.py`: `for it in range(calib_refit_iters): gdna_hyperprior =
_fit_gdna_hyperprior(...); belief = _sweep(gdna_hyperprior)`); the hyperprior enters ψ's gDNA arm as
`gdna_prior=`. The regression is visible without Phase 2.1 at all: on unstranded × capture-ON the refit is
worse than pass-0 (0.1148 → compare the refit arms in `/tmp/pass0_oracle_bench.tsv`). Dissect it with the
standing loop, which now has a sharper tool set:

* `scratchpad/dl_dissect.py [COND] --refit 1` — splits the error mass by DL-protection state (`τ_own>0` vs
  `τ_own=0`), strand DOF, and node class, and prints self-solve-vs-solved so you can see **how much the
  messages add**, then lists the worst nodes with full message provenance (`--only protected|unprotected`).
* `_capture["_dl"]` — the per-message DL gaps `G_g/G_p/G_n/G_λ`, the contradicted mask, and τ pre/post.
* `scripts/debug/gdna_hyperprior_eval.py`, `hyperprior_fit_options.py`, `refit_loop_study.py` — the existing
  hyperprior tooling.

The hypothesis to test first: on unstranded capture-ON the pass-0 solve that TRAINS the hyperprior is itself
capture-biased, so the fitted P(ρ) puts its mass in the wrong place and the refit then propagates that. Check
whether the hyperprior's modes track the oracle's gDNA-density distribution per capture regime before assuming
the refit loop is at fault.

### 3.3 The AMBIG residual — what is actually left, and where to see it

~50% of the residual stranded error mass sits on `τ_own = 0` AMBIG nodes, which DL leaves undamped by design.
**The minimal reproduction is a 3-node toy, not the 10 Mb suite:**
`tests/calibration/test_bp_solver.py::test_gdna_sweep_factor1_ambig_recovery` (xfail) — on a uniform ρ=0.5
chain, the AMBIG node between two intergenic anchors reads **0.3914** while the anchors are exact to 1e-9.

**It is not a mode defect and not a bug — it is the designed prior-free weakness.** Depth sweep (uniform
composition, counts scaled): 21.7% low at ρ=0.5, 14.2% at ρ=10, 5.6% at ρ=50, 2.1% at ρ=500, 0.8% at ρ=5000.
The transported mode is CORRECT and converges; what is missing is WEIGHT. An AMBIG node has `τ_own = 0`, so it
contributes no evidence of its own; it gets only its neighbours' messages at their honest `1/n` count
precision, against ψ's uninformative Jeffreys reference, which deliberately holds it off the `f_g=1` vertex
until the data earn it. **Do NOT attack this with more damping or by hunting a mode bug** (an earlier draft of
this handoff said "the error is in the message MODE" — that was wrong, and the depth sweep is the disproof).
The fix is a trained prior replacing the uninformative reference: §3.1 + §3.2.

## 4. Open items, ranked (each is real, each has numbers)

1. **The refit's unstranded-capON regression** — §3.2. The blocker.
2. **`f_cur` at unsolvable nodes** (`bp_solver.py`, in the `_RHO_ITERS` loop; a ⚠ comment marks it). ψ's output
   at a node with no free RNA strand is discarded by the write-back gate but still feeds the NEXT iteration's
   reframe — and the solver returns 0 for a node it never solved, so **every gDNA anchor is re-framed as if it
   were 100% RNA** (`ρ_tot` off by `E_g/E_r`, up to 1.8×, on every edge incident to an anchor). The one-line
   fix is `np.where(solvable, ..., f_g)`. It was measured (arm `dl4`) and the CORRECT arithmetic scores
   slightly WORSE: aggregate 0.0969→0.0972 (refit=0), 0.0828→0.0832 (refit=1), unstranded-capON 0.1702→0.1740.
   **Something downstream is compensating for it.** Investigate the pair together; do not flip it blind.
3. **The peel's zero-truncation** (`message_variance_derivation.md` §4 guard, still open). A fully-consumed
   peel emits `t_p = 0` — "no RNA continues past here" — at a live precision, which DL then has to kill as a
   contradiction. The derivation says `ρ_ν < 0` should be a PRIOR TRUNCATION, not an emission at zero. Fixing
   the emission is cleaner than having the cliff term clean up after it.
4. **The `u`-weighted peel and the M1 near-wall guard** (`message_variance_derivation.md` §4) — derived, still
   not wired.

## 5. INVARIANTS + gotchas (must preserve)

* The `√2·σ_own` pin-safety inequality (§1). Precision is discrete counts over a discrete length; enrichment
  scales the MODE, never the PRECISION.
* `N` enters only as power (`τ_λ` Fisher, or `1/n` sampling), never a composition vote.
* Variances in log-odds `Var(λ)` / `Var(log f_c)`, NEVER on simplex fractions.
* The composition is ONE DOF (single-λ message); the tilt θ is a SEPARATE DOF (its own message, AMBIG only).
* `v_own` for DL is composition-ONLY (no `1/n_dst`): `M_dst` cancels from both sides of the gap, so adding it
  would inflate the null and make messages STRONGER — the forbidden direction.
* Do NOT add DL to `_relay`: without a per-hop `_pin_v` the gap carries the uncancelled reframe residual, i.e.
  it would charge scale drift as composition mismatch — verbatim the `(log r)²` defect.
* Do NOT delete `rna_*_frac_var`, `strand_likelihood.py`, `lam_var`, `var_gdna` (all LIVE/referenced).
* Cached suite: `~/Downloads/rigel_runs/ambig_dense_10mb` (`_selfsolve_cache`, ~1 s/scenario).

## 6. ▶ KICKOFF PROMPT (copy-paste)

> We are developing Rigel calibration. The message-variance model is COMPLETE — derived, MC-validated,
> independently verified, implemented, and A/B-won (aggregate pass-0-vs-oracle 0.1267→0.0969 refit=0,
> 0.1234→0.0828 refit=1; the cross-cliff cost is now split into the M5 scale term `Var(log r)` and the
> DerSimonian–Laird composition-mismatch `b̂²`). Read `docs/calibration/ROADMAP.md`, then
> `docs/calibration/SESSION_2026_07_25_HANDOFF_6.md` (START HERE), then `message_variance_derivation.md`. Do
> NOT read `docs/calibration/archive/`.
>
> **The task is Phase 2 — the gDNA hyperprior, which is now the single binding constraint.** In order:
> (1) **fix the refit's unstranded-capON regression** (HANDOFF_6 §3.2) — dissect with the standing loop
> (`scratchpad/dl_dissect.py --refit 1`, `_capture["_dl"]`); first test whether the pass-0 solve that TRAINS
> the hyperprior is itself capture-biased, before blaming the refit loop. (2) Then re-apply the 6-line
> hyperprior→DL-`v_own` step (§3.1, already written and measured: it improves stranded 0.0376→0.0333,
> verystrong 0.1292→0.1196, capture-off 0.0354→0.0168 and ONLY regresses unstranded-capON, because where the
> hyperprior is wrong DL enforces it) — that also fixes the AMBIG residual, whose minimal reproduction is the
> 3-node factor-1 toy (`test_gdna_sweep_factor1_ambig_recovery`, xfail, reads 0.3914 vs 0.5). (3) Then the
> re-solve and a ship candidate. GOVERNING PRINCIPLE: pass-0 must be WEAK and correctable — preserve the exact
> invariant `p_eff = 1/max(v_msg, G²−v_own)`, i.e. a message out-weighs a node's own belief only if it agrees
> within √2·σ_own. No magic numbers. Counts are Poisson. `OMP_NUM_THREADS=1`. Derivation/audit workflows have
> worked well — use them. Regenerate goldens LAST.
