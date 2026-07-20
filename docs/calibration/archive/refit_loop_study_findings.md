# Refit-loop exploration — findings (2026-07-17)

**Status:** research exploration, not a ship decision. Tool: `scripts/debug/refit_loop_study.py` (mirrors
`calibrate.py`'s warm-continuing refit loop; instruments every pass). Substrate: the 7-condition
`ambig_dense_10mb/_selfsolve_cache` (gdna300; unstranded `ss_0.50` = the flagship corner, stranded `ss_0.99`;
capture on/off; nascent present/none). All numbers are per-region mass-weighted `|f_g − oracle|` (**mwae**) plus
the ground-truth-free **disagreement** instrument.

## The instruments (why more than one number)

Per refit pass we record, deliberately NON-reducible to a single scalar:

* **mwae** — oracle error (synthetic truth only).
* **disagreement vs the oracle FLOOR** — the optimization objective posed by the owner: the between-node
  squared adjacent log-density gap, per component (gDNA, RNA), reported *alongside the oracle's value on the
  same edges*. The truth has a NON-ZERO disagreement (real enrichment crossings + expression steps); the target
  is to MATCH the oracle floor, not reach zero. `disG_sol ≪ disG_orc` ⇒ over-smoothing (the degenerate
  collapse); `≫` ⇒ under-solved. This is the ground-truth-free lens (on real data we lose the oracle floor but
  keep the solved value + the total-density reference).
* **Δf_g** — mass-weighted per-pass change: convergence / oscillation.
* **P_g mode count + locations** — landscape drift.

## Findings

### 1. σ²_transfer is the DOMINANT variable, and its sign is regime-dependent
Turning the projection σ²_transfer on/off moves mwae far more than any refit or bandwidth choice, and in
OPPOSITE directions:

| condition | σ²_T ON (pass-0) | σ²_T OFF (pass-0) |
|---|---|---|
| flagship unstranded+capON (+nas) | 0.570 | **0.437** |
| flagship unstranded+capON (none) | 0.630 | **0.480** |
| stranded capON | **0.040** | 0.300 |
| stranded capOFF (+nas) | **0.011** | 0.040 |

σ²_transfer HELPS every stranded / off-capture regime (0.30→0.04 stranded capON) and HURTS the flagship
(0.437→0.570). This reproduces the earlier recorded A/B (flagship 0.48→0.63). **F1 is not wrong** — it correctly
gags messages that cross a real enrichment boundary. The flagship's problem is that those gagged messages (the
long-range gDNA propagation from confident intergenic anchors into enriched regions) are its ONLY signal:
strand carries zero information under an unstranded library, so correctly gagging the unreliable crossing
messages leaves the flagship information-starved. Disagreement lens confirms the mechanism: σ²_T ON collapses
the flagship `disG_sol` to 2.0 (oracle 8.9 — gDNA field far too flat, crossings missing); σ²_T OFF recovers
`disG_sol ≈ 6.7`. The missing spikes reappear as RNA (`disR_sol` 13.1 > oracle 10.7) — enriched gDNA leaking to
RNA, the flagship leak, seen in the disagreement decomposition.

### 2. The refit's sign is DOWNSTREAM of σ²_transfer
* With σ²_T **OFF** the refit HELPS (reproduces the historical `calibrate.py:237` note "0.151→0.118"): flagship
  0.437→0.422, stranded capON 0.300→0.189, stranded capOFF 0.040→0.033.
* With σ²_T **ON** (current production) the refit is net-HARMFUL on 6/7 conditions (flagship 0.570→0.603).
  Cause: a refit re-fits `P_g` on the solved belief, which SHARPENS it (deconvolved gDNA is more concentrated
  than total density); because `node_sweep` uses the one `gdna_prior` for BOTH the logprior AND the σ²_transfer
  projection, the refit sharpens the gagging too (modes proliferate 4→7). The two sharpenings compound. **Wiring
  the projection silently inverted the refit's benefit.**
* **Implication:** pass-0 (weak total-density prior, no refit) is the best mwae in 6/7 conditions today. The
  historical benefit of `calib_refit_iters=1` no longer holds with σ²_transfer ON.

### 3. Mechanism decomposition (cold-restart isolates it)
* The refit ACCURACY DRIFT = the prior-sharpening feedback (self-training). Present identically under cold
  restart (solve fresh from init each pass) — so it is NOT the warm-continue.
* The OSCILLATION = the warm-continue. Stranded capON warm-continues into a period-2 limit cycle
  (mwae 0.0384↔0.0411 forever); **cold-restart eliminates it** (0.0396→0.0394→converges) at no cost elsewhere.
  → cold-restart is strictly better; it removes the reviewer's oscillation risk.

### 4. Bandwidth is a minor knob
mwae is bandwidth-insensitive for the flagship base error (0.57 across h∈[0.10,0.25]). Bandwidth matters
modestly: smaller (0.10) better OFF-capture (a wide kernel blurs the single depleted mode); moderate (~0.25)
best for stranded capON (`disG_sol` 6.82 vs oracle 6.83). No single h is best everywhere → an eventual adaptive
h (per-regime resolution) is plausible, but it is a second-order lever, not the flagship fix.

### 5. Belief-confidence ≠ correctness in the flagship
The owner's "hold out ambiguous nodes" idea (refit `P_g` on the low-`Var(log f_g)` half only) HELPS where
confidence tracks truth — stranded capON 0.040→0.028, clean convergence — but barely moves the flagship
(0.603→0.599), because an unstranded library makes the confident nodes CONFIDENTLY WRONG (confidently flat).
Hold-out by belief-variance is sound only where the belief is calibrated.

### 6. Not a bug (one refuted, one doc fix)
* `belief.var_gdna` is `Var(log f_g)` (`simplex_logodds.py:321`), which is exactly the `τ²=Var(log g)` the
  refit fitter wants (g=f_g·M) — NO units bug. The `NodeBelief.var_gdna` docstring saying "Var(f_c)" is
  MISLEADING (should read "Var(log f_c)") — a doc cleanup, not a functional bug.

## What this says about the refit as an optimization
The owner's framing — iterate to minimize between-node disagreement, constrained by fixed total density + fixed
spliced + strand — is the right objective, and the disagreement-vs-oracle instrument operationalises it. But the
study shows the constraints are INSUFFICIENT in the flagship: with strand uninformative, the disagreement
minimiser has a degenerate freedom to put the enrichment crossings in RNA instead of gDNA, and both σ²_transfer
(by gagging the corrective messages) and the naive refit (by sharpening a biased prior) push it further into
that degenerate basin. The refit converges (Δf_g→0) — but to an over-smoothed fixed point, not the truth.

**The refit is not the lever for the flagship.** The lever is giving the strand-blind flagship a NEW source of
gDNA-vs-RNA information that does not route through the gagged messages. Candidate directions (for discussion,
not yet built):
1. **Adaptive prior strength** — the pass-0 prior is deliberately weak (n_eff≈0.15) so strand dominates; but
   where strand is uninformative the prior is the ONLY signal and should be STRONGER. Prior strength ∝
   (1 − strand information) is principled and count-free.
2. **Enrichment-correlation prior** — encode "high total density (enriched) ⇒ high f_g (gDNA)" directly (the
   physical capture invariant), rather than only the marginal density landscape. The current logprior half-does
   this; making the enrichment→gDNA coupling explicit could label enriched regions without the messages.
3. **Decouple the two jobs of the prior** — the refit should perhaps update the logprior (deconvolution) WITHOUT
   sharpening the σ²_transfer gate, or keep σ²_transfer from the belief-free pass-0 landscape. Currently they
   are tied to one `gdna_prior` and sharpen together (§2).
4. **σ²_transfer belief-gating** — the memory's long-stated intent: gate the crossing-gag on belief QUALITY, so
   a strand-blind node does not have its one informative message gagged as confidently as a strand-anchored one.

## Reproduce
```
OMP_NUM_THREADS=1 python scripts/debug/refit_loop_study.py --passes 5 --bandwidths 0.15
OMP_NUM_THREADS=1 python scripts/debug/refit_loop_study.py --bandwidths 0.10,0.15,0.25,0.40 \
    --conditions gdna_gdna300_ss_0.50_nrna_present_capture_on          # bandwidth sweep, flagship
OMP_NUM_THREADS=1 python scripts/debug/refit_loop_study.py --no-tv     # σ²_transfer OFF (old regime)
OMP_NUM_THREADS=1 python scripts/debug/refit_loop_study.py --cold --conf-frac 0.5   # oscillation + hold-out
```
