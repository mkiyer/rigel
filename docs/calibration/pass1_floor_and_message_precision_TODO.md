# Pass-1 global floor + message-precision — open problems (resume here)

**Status:** deferred behind the solver-robustness (nondeterminism) fix — that must land first, because until the
solve is robust no A/B on these is trustworthy. Branch `calib-gdna-accuracy`, 2026-07-05. Companion:
`robust_redescending_prior_design.md` (the shape study), `robust_floor_message_shape_results.md` (memory),
`calibrate_cross_process_nondeterminism.md` (the solver-brittleness memory).

Both problems below share ONE root theme: **the calibration mishandles LOW-PRECISION nodes** (unstranded /
AMBIG, whose strand gives little information). The floor over-pulls them and the messages over-perturb them —
even though their own value may be accurate, just imprecise. The fixes must respect "imprecise ≠ wrong."

---

## Problem 1 — the pass-1 global floor crushes enriched nodes (WRONG)

**The bug.** The pass-1 global prior is a two-sided quadratic in `r = log f_g − log(ρ_target·E/M)`:
`glob = −½·N·r²`, so its pull `−N·r` is **linear in the distance and unbounded** — it pulls *hardest* on the
nodes *farthest* from the floor, which are exactly the capture-enriched exons. On node 1085 (enriched AMBIG,
r≈7) it drags f_g down with force 7. **Crushing an enriched node toward a depleted floor is flat-out wrong.**

**Why it matters even though pass-2 masks it.** Pass-2's KDE currently dominates the final value (stranded is
inert to the pass-1 floor shape). BUT the pass-2 KDE is itself too strong (Problem 3, separate doc) and is
masking pass-1. When the KDE is corrected, **pass-1 values become critical** — and a pass-1 that crushes
enriched nodes to the floor will surface as a large leak. So the floor must be fixed on its own merits.

**The hard constraint — no-gDNA must still work.** A pure weakening or a naive redescending floor breaks the
no-gDNA regime: when the background ρ_floor≈0, every node is "far" from a ~0 floor, so a redescending
(pull-less-when-far) floor RELEASES the balanced-RNA nodes → they float up to false gDNA (measured: a
both-passes redescending floor gave ~106k false-positive gDNA vs baseline's 32). The spring's unbounded
pull-to-zero is exactly what suppresses no-gDNA. So the floor faces a genuine tension:
- **no-gDNA / flat-background nodes**: need a STRONG downward pull (suppress false gDNA).
- **capture-enriched nodes**: need RELEASE (don't drag them to the floor).

**The discriminator** cannot be the distance-from-floor alone (both enriched-gDNA and no-gDNA-flat nodes are
"far"). It must involve the node's OWN evidence: an enriched node has a high, well-supported density; a no-gDNA
node's apparent density is not gДНК. Candidate levers (to derive properly next):
1. **Weaker floor** (smaller pseudocount) — the easy option; reduces the crush but also the no-gDNA
   suppression; likely insufficient alone.
2. **Redescending, count-/precision-aware** (the "gravity" principle) — pull weakens with distance *relative to
   the node's own noise*, so a confident enriched node (small own-variance) releases while a
   sparse/uncertain node stays regularized. This is the right direction; the challenge is preserving no-gDNA
   (where ρ_floor≈0 makes "distance" huge for everyone). May need the target/scale linked to an OBSERVABLE
   density (ρ_floor) rather than an absolute, and a floor that is strong toward the *background* but not
   *below* a node's own supported density.
3. **A floor on gDNA DENSITY, not fraction** — physically gDNA ≥ background everywhere; the prior should assert
   `ρ_g ≥ ρ_floor` (a lower bound that props up, never a cap that crushes). This does not by itself fix the
   no-gDNA over-call (that's an over-call, above the floor) — that suppression must come from elsewhere
   (the strand tilt / messages / a separate no-gDNA mechanism), which is the coupling to Problem 2.

**Empirical anchors (deterministic in-process A/B, `compare_floor_shapes.py`, pass-1-only shape):**
- baseline (spring): capon-unstranded leak −194k, exon peak 0.28.
- cauchy pass-1: leak −129k, peak 0.39 — modest, safe, no over-call, no-gDNA preserved (+32), stranded inert.
- off / max pass-1: recover peak more (0.63–0.69) but flip to over-calling (+157k / +214k).
- Floor and message shapes STACK — combining over-corrects; pick one axis at a time.

**Requirement for the fix:** elegant, simple, principled, magic-number-free, and **no substantial
regressions** (no-gDNA must stay ~perfect; stranded must not regress). Derive it, don't tune it.

---

## Problem 2 — message disagreement-aware precision: `vloc[dst]` is CORRECT; do not remove it

**Current form** (per gDNA/±RNA message; `dst` = receiver):
```
resid    = mo − lfg_loc[dst]                      # surprise vs dst's message-free local belief
base_var = vbg[src] + pois                         # source belief var + source sampling var
s2_edge  = max( resid² − (base_var + vloc[dst]), 0 )   # excess (process) variance, MoM
pr       = 1 / (base_var + s2_edge)                # message precision
```
`vloc[dst]` (the destination's own belief variance) credits the surprise: **if the destination is itself
uncertain, it should trust an incoming message more.** This is CORRECT and must stay. Dropping it
("max-without-vloc") tightens the gate and improves the capon-unstranded leak (−161k → −34k) but at the cost
of a correct term and a stranded regression (−2.8k → −7.5k) — a leak-number win, not a real-data win. Rejected.

**The real finding.** In the unstranded setting node precision is LOW (the strand carries little info), yet the
node's VALUE may be accurate. So the danger is the opposite of over-trusting messages: **message passing must
be gentle enough not to perturb an already-accurate-but-imprecise node.** The lever is either (a) messages are
too strong, or (b) nodes are under-confident. The right fix improves calibration ACCURACY /
CONFIDENCE — e.g. give a genuinely-accurate unstranded node more deserved precision, or make the message
combine gentler in a principled (not term-deleting) way — rather than removing `vloc`.

**Do next:** characterize, on unstranded nodes, whether the message combine MOVES nodes away from their
accurate local value (perturbation) vs toward it (correction). If it perturbs, the fix is a gentler/／more
precision-aware combine that respects the count-zero-info principle (a low-count node's messages should not
override a many-count node's balanced-strand gDNA evidence). Tie to Problem 1 (both over-handle low-precision
nodes).

---

## Sequencing
1. **Solver robustness / nondeterminism FIRST** (`calibrate_cross_process_nondeterminism.md`) — the solve must
   not amplify ε; until then no A/B here is trustworthy.
2. **Pass-2 KDE** (too strong; masks pass-1) — the real peak-smoothing lever.
3. **This doc**: derive the principled floor (Problem 1) + the precision-aware message combine (Problem 2) on
   the now-robust solver, validate on the 16-scenario net-flow benchmark with no substantial regressions.
