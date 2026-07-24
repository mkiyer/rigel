


R1: plumb n_unspl_left/right into NodeGeometry + switch _pred_precision to integer counts

R2: decouple the spliced MEASUREMENT channel from σ²_transfer

R3: retire rho_g_cross geo-mean hack

R4: neutralize σ²_transfer cliff term behind a switch + per-condition A/B (stranded AND unstranded)

R5: measure residual disagreement; R6: seam derivations (TSS/TES, exon-exon AMBIG) → then retire density mode

A': solve-gate reframed as honest precision-0 state (NOT the refuted f_g=1 lock); E: precision merge; L: re-measure nascent≫mature after mode fixes
The consolidated roadmap is now pass0_roadmap.md §8 (supersedes the old §7, which predated this thread). Complete open set:

8.2 — Message-arithmetic reconciliation (active)

R1 MASS→COUNT — audited, next. Only the per-face unspliced count is missing from geometry.
R2 decouple spliced measurement from σ²_transfer.
R3 retire rho_g_cross.
R4 neutralize the σ²_transfer cliff term — the big lever, behind a switch, per-condition A/B.
R5 measure residual disagreement; decide if any damping survives.
R6 seam derivations (TSS/TES, exon–exon AMBIG) — completing R6 is what lets the density mode be retired entirely.
8.3 — Carried forward

A′ solve-gate reframed — the honest precision-0 state, distinct from the refuted f_g=1 lock, plus the hyperprior handoff contract.
E precision prediction⊕measurement merge.
L nascent≫mature — flagged for re-measurement, since the mode fixes may have substantially addressed it; I don't want it carried as open on stale evidence.
Depletion bias — bounded, one-sided, contained; not fixable by precision (it's a mode bias); cross-junction detector is the real fix.
8.4 — Standing gates on every behavioural stage, with goldens regenerated last, only after the per-condition A/B review.

One framing I'd record: the message-arithmetic thread wasn't a digression from §7 — it found that the mode was structurally wrong across the capture cliff, which is the root cause of the unstranded+capture-ON failure §7 was trying to route around. Your "unstranded + capture-ON has never worked" and σ²_transfer's history are the same fact: the damping was suppressing messages because they were wrong, and strand-specific data was strong enough to mask it.

Nothing is committed yet — the working tree has the two mode fixes, the toolkit, injected_priors, and three new docs, with 5 golden tests deliberately red pending the A/B.

Ready to start R1 (plumb the per-face unspliced count into NodeGeometry, switch _pred_precision to integer counts) followed by R2 — both unambiguous correctness fixes, cheap to verify, no behavioural judgement calls.






Remaining tasks — the path to ship pass-0
In priority order (each: derive → plan → execute; A/B-gated; gdna_none guard is hard):

⟵ SHIP-BLOCKER — item L: the nascent ≫ mature over-call. The next task. This is a confident-wrong pattern (over-calling gDNA where the mass is nascent RNA) — precisely what the ship criterion forbids. It's the unstranded nascent-vs-gDNA identifiability floor. Two acceptable outcomes: (a) the intron→exon nascent relay resolves it (a real fix), or (b) the over-call is driven to low precision (honest-weak, deferred to the hyperprior). Either way it must stop being confidently wrong.

The solve-gate (item A) — the honesty refinement. Decide withhold-vs-solve for the coin-flip nodes; needs the correlation/precision metric + counting the fitted prior as an identification source. (Reframed: no boundary inversion — that was the artifact; the open part is the region withhold decision.)

Correctness items (I). Discretization-frame consistency (+1 discrete vs continuous eff-lengths on short flanks); the AMBIG two-root stance; confirm the identifiability-wall deferral path fires. Smaller hardening.

SHIP pass-0 → then the hyperprior (§5). The real under-call lever, fed an honest pass-0. (Lower-priority polish — D seam-anchor, E precision-merge — can follow.)

Work order for the next session
Focus: item L (the ship-blocker). Following the mantra — (1) design & derive, (2) plan, (3) execute:

Extend the grid (quintuple_grid.py) into a full sweep + turn it into a regression test (the solver must "just work" across gDNA × nascent × mature). Instrument the boundary f_g and the intron f_g too, not just the middle exon.
Diagnose the over-call mechanism on the quintuple: is the intron correctly resolving as RNA (low f_g) but failing to relay that to the exon? Or is the exon's over-call at high precision (confident-wrong) vs low (honest-weak)?
Derive the fix: the intron→exon nascent relay (the telephone) — the introns are transcribed and well-defined; they should tell the exon "your unspliced mass is nascent, not gDNA." Or confirm the over-call is honestly weak and hand it to the hyperprior.
Verify across the grid + the gdna_none guard + the panel (no boundary-mode regression) → then move to the solve-gate + correctness items → ship.