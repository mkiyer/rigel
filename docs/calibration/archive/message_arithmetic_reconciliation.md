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

---

## Part E — HOLISTIC AUDIT (2026-07-22): the node inventory, PEEL vs GRAFT, and what blocks item E

Owner's framing, which reorganizes everything below: **pass-0 solves the UNSPLICED fragment pool only.** The
simplex `(f₊, f₋, f_g)` partitions *unspliced* fragments. Spliced fragments are guaranteed RNA, are directly
measured, and are **not** part of the solve. Every message question reduces to: *which components does the
destination's unspliced pool contain, and how does the source supply them?*

### E.1 The node inventory

| node | UNSPLICED pool | spliced channel |
|---|---|---|
| intergenic region | g | — |
| intron region | g, ν (nascent) | — |
| **exon region** | g, ν, **μ** (within-exon mature — contained, crosses no junction) | — (`spliced_* ≡ 0` for every REGION node) |
| **splice boundary** | g, ν (mature splices away ⇒ the crossing is mature-FREE) | **S_B** (measured mature) |

Two facts do all the work: an exon's unspliced pool **contains** mature, a junction crossing **does not**, and
only the boundary owns a spliced channel.

### E.2 One rule, three operations

A message expresses the source's densities **in the destination's inventory**:

1. **TRANSPORT** — `M_c^dst = ρ_c^src · E_c^dst · g_c` (per-component eff-length frame change + the structural
   gate). Enrichment cancels. This is the shift, and it is the shared spine of every edge.
2. **PEEL** — drop a component the destination does not have.
3. **GRAFT** — add a component the destination has and the source can supply from a *different* channel.

The two splice-junction directions are **not mirror images**, because the endpoints know different things:

| direction | operation | mechanism | why |
|---|---|---|---|
| exon → boundary | **PEEL** μ | additive **shift** `+c_b`, `c_b = log1p(S_B/D_B)` | the exon **cannot decompose its own RNA** — a region has no spliced channel. Only the **dst** measures the mature share (its own S_B/D_B). A shift in the dst's frame ⇒ no eff-length extrapolation. |
| boundary → exon | **GRAFT** μ | ordinary **density** term `ρ_μ = S_B/E_spl`, transported by rule 1 | the dst has μ and the **src measures it**. Nothing special: one more component through the shared transport. |

> **PEEL uses a shift; GRAFT uses a density. Never both on one edge.**

This asymmetry is not a special case to be tidied away — it is forced by *who holds the measurement*. It is
also what makes item **E** well-posed: only a **separated** ρ_μ / ρ_ν can carry provenance-distinct precision.
A shift has nowhere to put two provenances.

### E.3 LANDED — the mature was restored TWICE on boundary→exon

The code grafted the mature into the shift normalizer (`+SPs/_esp` inside `_den`) **and** applied `−c_b`.
Measured on the grounded full-transcript toy (capture ON, unstranded, oracle f_g = 0.754):

| edge | lump-only | c_b-only | **BOTH (old code)** |
|---|---|---|---|
| B(6)→exon(7) | 0.5020 | 0.5748 | **0.4615** |
| B(8)→exon(7) | 0.5050 | 0.5807 | **0.4629** |
| B(4)→exon(3) | 0.4983 | 0.5847 | **0.4502** |

Strictly below either single accounting ⇒ a **systematic gDNA under-call at every exon**. `−c_b` removed
(the graft is the mechanism the owner's model names, and the one E needs).

**Gates:** `test_mature_measurement_disagreement_silenced` **red at HEAD → green** (it was measuring exactly
this bias). Calibration + native suites 306 pass. `gdna_none` phantom guard **3,766,743 → 3,821,731 (+1.46 %)**
— concentrated in `ss_0.50_nrna_present_capture_off` (774,030 → 815,881). The guard is **one-sided**: in a
zero-gDNA library any correction that *raises* f_g must score worse, and this one raises f_g by construction.
The guard cannot distinguish "correctly raised" from "phantom"; the adjudicating evidence is the toy oracle
(0.46 → 0.50 against 0.754) and the now-green disagreement test.

### E.4 ⛔ WHAT BLOCKS ITEM E — `Var(log f) = ∞` on a composition-vacuous source

Measured on the same toy (unstranded), boundary→exon edges:

| edge | τ_src | `v_logfg` | **prec_g** | **prec_p** |
|---|---|---|---|---|
| B(6)→exon(7) | 0 | **inf** | **0.000** | 2.074 |
| B(8)→exon(7) | 0 | **inf** | **0.000** | 1.785 |
| B(4)→exon(3) | ~1e-10 | **inf** | 0.039 | 2.032 |

On unstranded data a boundary has no intrinsic composition evidence (I_strand ≡ 0 at κ=½), so `τ = 0` and
`Var(log f_c^src) = ∞`. Therefore:

* **the gDNA message into an exon carries ZERO precision** — every mode correction in this document is
  multiplied by zero on exactly the data we are trying to fix;
* the **only** surviving channel is `pr += S`, which is why the RNA factor is loud (≈2.0) and the gDNA factor
  silent. `pr += S` is not merely wrong — it is currently **load-bearing**.

**Implementing item E exactly as derived would make this worse, not better.** The share-weighted rule gives
`Var = w_μ²·(1/n_μ) + w_ν²·∞ = ∞` whenever any nascent is imputed ⇒ `pr = 0` ⇒ the measurement is silenced
too, and the exon hears *nothing*. That is the τ-gag bug class returning through the front door.

**The diagnosis:** `Var(log f) = ∞` is an **artifact of the log transform**, not a statement about the physics.
`f` is a **fraction on [0,1]**; its variance is bounded (`Var(f) ≤ ¼` for any distribution on the unit
interval). "No evidence" means *uniform*, not *infinite*. The infinity is manufactured by writing a bounded
quantity's uncertainty in log space and then letting τ→0.

**Consequently the true statement a vacuous source makes is one-sided:** with `ρ_ν = f_r·(sm/E_r)` unknown over
`f_r ∈ [0,1]`, the total `ρ_r = ρ_μ + ρ_ν` lies in `[ρ_μ, ρ_μ + sm/E_r]` — a **bounded interval**, so a finite
variance exists, and the measured mature is a genuine **lower bound** on the destination's RNA (§5 of the merge
derivation) rather than a point claim. A Gaussian factor pinned at that lower bound would over-claim and inflate
f_g — which is precisely the +49 % R2 signature.

**Open decision (owner):** how a composition-vacuous source bounds its imputed part. Three candidates, none
requiring a tuned constant:

| option | `Var(f_r)` when τ=0 | note |
|---|---|---|
| **(a) population prior** | the global gDNA prior's variance at that node | architecturally right — pass-0's *third* information source is the population baseline, and a node with no strand evidence still has it. Uses an already-fit quantity. **Recommended.** |
| (b) max-entropy uniform | 1/12 | a theorem, not a tuned value, but asserts a specific distribution |
| (c) Popoviciu bound | 1/4 | the most conservative finite bound; weakest messages |

Any of the three makes `v_log*` finite everywhere, retires `_pred_precision`'s `return 0` branch, un-silences
the gDNA channel on unstranded data, and lets item E land **exactly as derived** with no special case.

### E.5 Two further defects found (not yet fixed)

1. **`exon → boundary` emits a log-fraction ABOVE 1** — measured `exp(mode_g) = 1.107`. The exon's RNA is
   suppressed to zero in the normalizer (§14 change-1), so the shift already says "pure gDNA, f_g = 1", and
   `+c_b` then pushes it *past* 1. The peel presumes the normalizer still **contains** the RNA it is peeling
   mature out of. Suppression and peel are mutually inconsistent; exactly one may stand.
2. **The two λ-factors into an exon are incoherent** — gDNA rides the shift, RNA-total rides the *density* mode
   (`use_shift` is False whenever either endpoint is an exon). Measured sum of the emitted fractions:
   **0.753, 0.768, 0.767** — they should sum to 1. The fold is currently reconciling a self-contradiction.
   Putting both factors on the shift with the shared normalizer `_den` makes them **one coherent statement**
   (`f_g + f_r ≡ 1` by construction) instead of two competing claims.

### E.6 The consolidation this implies

Every message precision in `_scan` is one primitive — a density built from **provenance-tagged parts**:

```python
    ρ_c = Σ_k ρ_k ,   Var(log ρ_c) = Σ_k (ρ_k/ρ_c)²·v_k ,   then  +σ²_xfer  once per edge
      imputed-from-solve : v = 1/n_unspliced + Var(log f_c^src)
      measured-spliced   : v = 1/n_spliced                     (NO composition term — this IS "spliced is pure RNA")
      peel-shift         : additive on the mode; contributes Var(c_b)
```

`_pred_precision(count, v_log, s2t)` is exactly the **single-part** case of this. So the general form *subsumes*
it — and `pr += S`, the `n_eff` honest clamp, and the mass-vs-count inconsistency in the measurement channel
(`SPs` is MASS used as a count; `SPN` is the integer flux) all disappear rather than being patched.

### E.7 The gDNA FACTORY — owner's locus model, and the seam mode it implies

Owner's framing (2026-07-22), which supersedes the "vacuous factor" options in E.4:

> A **locus** is a genomic span of overlapping transcripts, bookended on both sides by intergenic regions.
> An intergenic region holds a **measured** amount of pure DNA — guaranteed true, no solve, invariant — so it is
> a **factory** of observed gDNA that should propagate into the locus at counting precision, with no penalty.
> It is also the **sink**: no RNA crosses it, and a message that reaches it dies. Stopping messages there is
> also what prevents averaging gDNA across a whole chromosome.

**The factory already exists and is being severed one hop out.** `struct_lock` (composition certainty) is
`locked & is_region` — so an intergenic REGION is certain, but the **seam boundary** between it and the first
exon is not. Measured on the grounded toy: `lock_src = 0` on **every** boundary, hence `Var(log f_g) = ∞` and
**`prec_g = 0.000` on every seam→region edge**. The measured background never enters the locus.

Certainty is a **structural** property, not a node-kind property: `locked ⇔ admits no RNA strand ⇒ every
unspliced fragment is gDNA`. That argument does not care whether the node is a region or a seam.

**But certainty alone is not enough — it needs the right MODE.** A seam is exactly where the source is *blind*:
it carries none of the destination's RNA. Emitting a composition SHIFT there asserts "the dst is pure gDNA"
(measured: **f_g = 1.0000** into an exon whose oracle is 0.754). This gives the missing predicate for A.2's
undecided rows:

> **Gate-equality decides the mode.** The SHIFT is valid iff the source can supply **every component the
> destination admits**: `gate_unequal = (fp[dst] ∧ ¬cont_p) ∨ (fn[dst] ∧ ¬cont_n)`. A blind source may claim
> only a gDNA **density**, converted to a fraction by the destination's own observed total — which is exactly
> `_mode_density`. **The density mode is therefore not retired; it is DEMOTED** to its correct role: the
> one-sided gDNA **floor** at a gate-unequal seam. R6 is answered rather than completed.

And the floor must be **one-sided**: measured 0.395 at the TSS seam vs a destination already at 0.510 — as a
symmetric Gaussian that bound would *pull a better estimate down*.

### E.8 PROTOTYPED AND REVERTED — what the factory still lacks

Implemented all three (structural `struct_lock` + gate-equality + one-sided floor) and measured:

* **It works as designed on the toy.** `lock_src` 0 → 1, `Var(log f_g)` inf → 0 at every seam, seam edges
  correctly fall to the density mode, and the floors (0.395, 0.415) are correctly **suppressed as non-binding**
  against destinations already at 0.510.
* **The phantom guard IMPROVES:** 3,821,731 → **3,814,453**.
* **It breaks the phantom RED LINE on the prior-free path.** With no enrichment prior `σ²_transfer = 0`, so a
  terminal G1 lock asserts its floor at raw counting precision (n = 47–124) and pins a vacuous exon at
  **f_g = 0.9616** (`test_tau_gag_fix_deconvolution_prediction_stays_gated` requires 0.2 < f_g < 0.8).

**The missing piece is not the certainty — it is a sound statement of the FLOOR.** A source's gDNA *density*
bounds a destination's gDNA *fraction* only where enrichment is comparable. Across a cliff `_mode_density`
returns "fractions" of **186, 294, 324** (measured, boundary→intron). A floor that can exceed 1 is not a floor.

> **This reframes R4.** σ²_transfer is *not* merely a legacy proxy for a mode that is now correct — it is the
> mechanism that currently makes a density-based floor safe across an enrichment change. Neutralizing it would
> remove the only thing keeping the factory honest. **R4 must not be attempted before the floor is derived.**

Reverted to keep the tree shippable; both sites carry the derivation as comments. **Next step is a derivation,
not an implementation:** what does a pure-gDNA source at density ρ_g, observed under enrichment e(src), validly
assert about a destination under enrichment e(dst)? That single question gates the factory, the seam modes, and
R4 — and it is the same "compare gDNA densities across differing capture" question the DNA-prior projection
work (`dna_prior_projection_resume.md`) is already circling.
