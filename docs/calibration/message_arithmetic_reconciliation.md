# Message-passing arithmetic: MODE + PRECISION reconciliation and cleanup

**Status:** derived; evidence gathered; **cleanup not yet executed** — awaiting owner review.
**Scope:** pass-0 message propagation (`bp_solver._scan`). TSS/TES seams out of scope (noted, deferred).
**Date:** 2026-07-22.

---

## Part A — MODE: one rule, applied everywhere

### A.1 The unified imputation

Every message imputes the destination's composition from the source's per-component **densities**:

```
    M_c^dst  =  ρ_c^src  ·  E_c^dst  ·  g_c(src→dst)                     (imputed mass of component c)
    mode_c   =  log( M_c^dst / Σ_c′ M_c′^dst )  +  Δ_c
```

* `ρ_c^src` — the source's density of component *c* (gDNA, RNA₊, RNA₋).
* `E_c^dst` — the destination's **per-component** effective length. This is the whole of the frame change:
  gDNA-FL ≠ RNA-FL, and contained (region) ≠ crossing (boundary).
* `g_c` — the **structural crossing gate**: 1 iff component *c* physically reaches the dst (gDNA always;
  RNA strand *s* only where *s* is continuous on BOTH endpoints — `free_pos`/`free_neg`).
* `Δ_c` — a **component-set reconciliation** term, nonzero only where src and dst carry different components.

**Enrichment cancels identically.** `ρ_c^src` carries `e(src)`; so does the normalizer `Σ M_c′`. The ratio is
therefore a pure composition transfer, and `e(dst)` never enters at all. This is why the shift is correct across
an arbitrary cliff.

**The density mode is structurally wrong across any enrichment change.** `_mode_density(ρ_c, E_c, md)` divides
by the destination's **OBSERVED** total `md/E`, which carries `e(dst)`. That factor does *not* cancel — it is
precisely the ~600× under-imputation measured at `B→exon` (imputed f_g = 0.0016 vs oracle 0.747). **The density
mode should be retired, not tuned.**

### A.2 Per-edge-type table (the standardization target)

| src → dst | src components | dst components | Δ | current | target |
|---|---|---|---|---|---|
| intron → IE/EI boundary | g, ν | g, ν | 0 | shift ✅ | shift |
| boundary → intron | g, ν | g, ν | 0 | shift ✅ | shift |
| exon → IE/EI boundary | g, ν, **μ** | g, ν | **+c_b** | shift+c_b ✅ | shift+c_b |
| boundary → exon | g, ν | g, ν, **μ** | **−c_b** | shift−c_b ✅ | shift−c_b |
| intergenic → TSS/TES boundary | **g only** | g, ν | ? | density ⛔ | **needs derivation** |
| boundary → intergenic | g, ν | **g only** | ? | shift | **needs derivation** |
| exon → exon–exon seam (AMBIG) | g, ν, μ | ? | ? | density ⛔ | **needs derivation** |

`c_b = log(1 + S_B/D_B)` — the mature-dilution term, enrichment-invariant, zero constants
(`exon_boundary_mature_dilution_plan.md`).

**Reading:** the four splice-junction edges are done. Everything still on the density path is a *seam* case
where the component sets differ by RNA **presence** rather than by mature. The structural gate `g_c` already
encodes "RNA does not cross here" — the open question is whether the gate alone is the correct Δ for a seam, or
whether a seam needs its own term. That is the remaining mode derivation.

---

## Part B — PRECISION: term inventory and the obsolete proxy

Current message precision (`_pred_precision`, all channels):

```
    pr = 1 / (  Var(log f_c^src)  +  Var(c_b)  +  σ²_transfer  +  1/N_src  )
```

| term | what it represents | verdict |
|---|---|---|
| `Var(log f_c^src)` | the SOURCE's own composition uncertainty (reference-free, from τ). Epistemic. | ✅ correct |
| `1/N_src` | Poisson counting noise on the source density. Aleatoric, irreducible. | ⚠️ **uses MASS, not COUNT** |
| `Var(c_b)` | counting noise of the mature correction | ✅ correct (now integer counts) |
| `σ²_transfer = var_proj[dst] + (μ_proj[dst] − μ_proj[src])²` | enrichment-cliff damping | ❌ **obsolete proxy** |

### B.1 σ²_transfer is a proxy for a mode that is now correct

σ²_transfer was introduced to **dampen messages crossing enrichment cliffs so they didn't corrupt the solve** —
which was the right call *when the mode was the density mode*, because that mode is genuinely wrong by `e(dst)`
across a cliff. It is an *empirical stand-in for model error*.

Once the mode is enrichment-invariant, the cliff height carries **no information about the reliability of the
transferred composition**. Continuing to damp by `(μ_dst − μ_src)²` does not merely waste information — it
cancels the fix.

**Measured, after the mode corrections landed:**

| edge | mode | σ²_transfer | **prec_g** |
|---|---|---|---|
| `intgc(1)→B(2)` (carries ρ_bg) | density | 23.5 | **0.04** |
| `B(4)→intron(5)` | **shift** | 29.3 | **0.00** |
| `intron(5)→B(6)` | **shift** | 29.9 | **0.00** |
| `B(8)→intron(9)` | **shift** | 29.0 | **0.03** |

σ²_transfer (23–30) dominates the denominator; every cliff-crossing message is throttled to ~0.03 precision.
This is why the interior introns moved only 0.49 → 0.549 instead of toward 0.87: **the modes are right and
carry no weight.**

Note also that the old total-density **disagreement** estimator was already removed from production
(`bp_solver.py:25` — "relocated to `scripts/debug/`"), so σ²_transfer is now the *only* damping term. There is
no second mechanism holding the line if it is removed — which is exactly why it must be replaced deliberately,
not deleted casually.

### B.2 What should replace it

The honest residual is **model error of the transfer rule**, not cliff height: the nascent density ν is only
approximately continuous across a junction, ρ_bg is only approximately uniform (CNV), the FL model is
estimated. That term is real but should be (a) **small**, and (b) **independent of the enrichment cliff**.

Proposal: replace the cliff proxy with a residual measured **after** the corrected transfer — i.e. the observed
disagreement between a node's own belief and what its neighbour's corrected message predicts. That is the
disagreement-variance idea, but fit on the *corrected* transfer instead of on raw total density. **It must be
measured before it is trusted** — the empirical question is whether any damping is warranted at all once the
mode is right.

### B.3 Two further precision defects found in the audit

1. **MASS used where COUNT is required.** `_pred_precision(sm, …)` uses `sm` = facing unspliced **mass** as the
   Poisson count. The codebase states the rule explicitly (`node_geometry.py:81-87`): *"MASS is the correct
   numerator for a DENSITY; the COUNT is what a Poisson VARIANCE needs (`Var(log ρ)=1/n`, not `1/mass` — mass
   sums fractional per-fragment shares; Kish `n_eff=(Σw)²/Σw² ≥ mass`)."* This affects **every** message, not
   just the new one. `statics.u_pos+u_neg` is the integer flux and is available. (Already fixed inside
   `Var(c_b)`; **not** yet fixed in the main path.)
2. **σ²_transfer damps the spliced MEASUREMENT channel** (`pr += SPs/(1+SPs·s2t)`, lines 827/836). A junction
   read count is a direct observation; attenuating it by an enrichment cliff has no justification even under
   the old rationale. Should be independent of σ²_transfer.

---

## Part C — Cleanup plan (proposed, staged, each stage independently gated)

1. **Retire the density mode** for the four splice-junction edge types (done) and delete `rho_g_cross`, the
   unweighted pre-scan geo-mean, now that both exon↔boundary directions carry real precision-bearing messages
   and the α⊗β integrator (already a precision-weighted geometric mean) reconciles the two flanks.
2. **Fix MASS→COUNT** in `_pred_precision` call sites (B.3.1). Expect a uniform, modest precision change.
3. **Decouple the spliced MEASUREMENT channel from σ²_transfer** (B.3.2).
4. **Neutralize the σ²_transfer cliff term** (B.1) — behind a switch, A/B'd, since this is the single largest
   behavioural lever and there is no second damping mechanism in production.
5. **Measure the residual disagreement** after 1–4 and decide whether a replacement damping term is needed at
   all (B.2). Only then remove the switch.
6. **Derive the seam cases** (A.2 rows 5–7: TSS/TES, exon–exon AMBIG) and finish the mode standardization.

Gates for each stage: `msg_audit` direction, the mature-dilution identity test, `gdna_none` phantom guard,
calibration suite, then benchmark A/B before goldens.

**Ordering rationale:** 2 and 3 are unambiguous correctness fixes and should land first. 4 is the big lever and
must be A/B'd on the real suite, not just the toy — σ²_transfer currently protects strand-specific data (where
strand is more reliable than any message), and removing it there could regress even while it helps unstranded.
That asymmetry is the main risk in this whole reconciliation and should be measured explicitly, per-condition,
on stranded **and** unstranded arms.

---

## Part D — MEASURED OUTCOMES (2026-07-22): R1 lands, R2 is blocked on E

Bisected against the `gdna_none` phantom guard (9 zero-gDNA scenarios; truth f_g=0, so all reported mass is
fabricated). Baseline = HEAD solver.

| variant | region | boundary | total phantom | vs baseline |
|---|---|---|---|---|
| baseline (HEAD) | 3,354,039 | 387,498 | 3,741,537 | — |
| mode fixes only (intron shift + exon `±c_b`) | 3,325,827 | 440,693 | 3,766,520 | **+0.7 % ✅** |
| **mode fixes + R1 (MASS→COUNT)** | 3,326,037 | 440,706 | **3,766,743** | **+0.7 % ✅** |
| mode fixes + R1 + R2 | 5,040,371 | 545,917 | 5,586,288 | **+49 % ❌** |

**R1 LANDED.** Integer counts are byte-neutral on the guard (3,766,743 vs 3,766,520) and theoretically correct
(`Var(log ρ)=1/n`). `n_unspl_left/right` now plumbed through `NodeGeometry`, mirroring `spliced_n_*`.

**R2 REVERTED — and the reason is a real dependency, not a tuning failure.** Undamping the measurement is
correct in principle (a junction count is not an imputation, so an enrichment cliff must not attenuate it), but
`pr += S` attaches the MEASUREMENT's confidence to the **PREDICTION's mode**. On exon→boundary edges the mature
absorption (`rho_pos = SPs/esp − absorb_p`) can drive that mode to a clamped ~zero; undamped, this laundered a
weak "no RNA" into a **confident** one → f_g→1 → +49 % phantom gDNA. σ²_transfer was silently holding that
unsound merge together.

> **A measurement's confidence may not attach to a prediction's mode.** R2 therefore requires item **E
> (prediction⊕measurement MERGE)** first — the measurement must carry its OWN mode (the spliced-derived RNA
> density), after which its precision can be undamped safely.

**Revised ordering: E precedes R2.** This also re-frames item E from "a precision refinement" to a
**prerequisite for correct message arithmetic**.

**Method note.** This was found only because the phantom guard is a *delta* against baseline, not an absolute.
Every behavioural stage from here must be A/B'd against HEAD the same way — an absolute number is unreadable.

### D.1 R3 — `rho_g_cross` retired (LANDED, inert)

After the mode fixes, `rho_g_cross` was reachable **only** on the intergenic→boundary (TSS/TES seam) edge —
the intron↔boundary and exon↔boundary paths bypass it. A/B measured it **byte-identical** on the `gdna_none`
phantom guard (3,766,743 both ways), the grounded full-transcript toy (all node beliefs identical), and the
calibration / native / golden suites. Its remaining destinations are structurally pinned, so it had no effect.
Removed.

The seam edge now keeps the source's own densities on the density mode, pending its own derivation (R6).
**Standing conclusion:** the unweighted pre-scan geo-mean was only ever a stand-in for the missing
exon→boundary imputation; the α⊗β integrator (a precision-weighted geometric mean) is the correct place to
reconcile a boundary's two flanks, and it now receives two real, precision-bearing messages to do it with.
