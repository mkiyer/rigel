<!-- title: Calibration — state of the work, findings, and resume plan (bookmark) -->
# Calibration solver — state of the work + resume plan (2026-06-27)

A reference snapshot for resuming the **calibration-accuracy** effort: where the solver stands, what the
node-by-node + bridge dissections proved, the bugs found, and the prioritized next steps. Companion to
`forward_backward_plan.md` (the FB design) + `CALIBRATION_ARCHITECTURE.md` (the count-zero-info principle).

> **Framing (2026-06-27, user):** the goal is to **improve CALIBRATION accuracy** (the per-node `f_g`) and close
> the gap — *not* to band-aid the downstream EM/bridge. The §11 gDNA-eff-len finding is a real ~200K lever but it
> is a **bridge/EM** fix (`priors.assemble_priors`), not calibration; pursue it separately, do not let it stand in
> for closing the calibration gap.

## 0. TL;DR

- **Production = the committed Jacobi solver, `main`/branch @ 346f1531.** This session's exploration (FB rewrite +
  the ê/mode fixes + dedup + all harnesses) is **preserved on branch `calib-fb-dissection-jun27`**, NOT merged.
- **FB is correct tree-BP and fast** (~10× Jacobi) but **net worse on the benchmark** (+11% leak) — *only* because
  Jacobi is accidentally regularized by under-convergence (FB == Jacobi@∞). FB is not the win.
- **The node-by-node dissection overturned "FB over-smooths":** propagation is only 4–19% of the error and is
  net-*correcting*. The error lives in the **LOCAL belief** — the strand likelihood (empty at κ≈½) + the **global
  ê(z) prior** (biased/noisy, worst on AMBIG exons).
- **Two genuine bugs found (both on the explore branch, validated-correct, NOT on main):** (1) the **ê(z) eff-len
  bug** — `fit_enrichment_transfer` projected ρ_mature with the *full* RNA eff-len instead of the half-triangle
  `eff_spl` (the same spliced-fix inconsistency, missed in this path); fixing it cut the zero-gDNA ê phantom ~85%
  but **exposed the capture cliff** (flagship total flipped worse) → it needs the robust-message co-fix to be
  net-positive, so it was NOT pushed standalone. (2) the **unbounded message mode** (clipped to [0,1]).
- **§11 NEW — the calibration→EM bridge audit:** `assemble_priors` tabulates boundary/region **mass correctly**
  (conserved), the **prior gDNA split is accurate** (vs oracle, mean|Δ|≈0.05), but the **gDNA component eff-len is
  ~2× too large** → gDNA structurally under-competes → the **~200K stranded-flagship leak**. Halving it collapses
  the leak 176K→1K. This is a **bridge fix, not calibration** (see Framing above).
- **Decision:** keep Jacobi in production; the calibration-accuracy lever is the **local belief (ê(z) + the cliff
  robust message)**; the bridge eff-len is a separate, confirmed, high-value EM fix.

## 1. The journey

1. **GS (single loop, var~mean refit inside the sweep)** — never converged (a measured limit cycle).
2. **Nested loop** (outer var~mean fixed point / inner belief solve) — fixes convergence. Committed.
3. **Jacobi inner** (batched, both-neighbour) — `main`/branch @ 346f1531. ~3.3× faster, benchmark net leak −9%
   vs single-loop. **This is the shipped production solver.**
4. **FB inner** (forward α + backward β, true tree-BP, single pass) — working tree, **uncommitted**. ~10× faster
   than Jacobi but **+11% worse on the benchmark** → triggered the dissection.

## 2. What FB is

`bp_solver.node_sweep`, per OUTER iteration (the var~mean fixed point): (A) one batched message-free LOCAL
lattice solve → per-component (fraction, precision); (B) a FORWARD scan L→R accumulating the left-context
message α + the forward belief (`local ⊗ α`, *not* the reverse — true tree BP, so a thin node relays the
upstream α); (C) a BACKWARD scan R→L → β; (D) combine `α⊗β` and one batched FINAL solve. Two batched lattice
solves + two O(m) moment-space scans. Single-pass (exact on the chain given the local beliefs); the outer loop
iterates the var~mean. `max_passes`/`convergence_delta` are retained but unused.

Timings (flagship, 2443 nodes, max_outer=1/2/4): **FB 2.0/2.6/3.7 s**, Jacobi 10.3/19.6/38.0 s, GS 32.7/64.5/127.8 s.

## 3. The problem

Benchmark (`quick_3to1_5mb`, 16 conditions, vs the committed Jacobi `.jacobi`):

| regime | Jacobi Σ\|leak\| | FB Σ\|leak\| | Δ |
|---|---|---|---|
| **total** | 1,928,732 | 2,143,061 | **+11%** |
| flagship (CAP, gDNA) | 1,138,381 | 1,386,575 | +22% |
| zero-gDNA | 399,770 | 460,631 | +15% |
| capture-off, gDNA | 390,581 | 295,855 | **−24%** |

The pattern (only capture-off improves) first *looked* like FB over-smoothing the capture density cliff.

## 4. The toy (over-smoothing hypothesis — then overturned)

`scripts/debug/fb_cliff_toy.py` — a Gaussian-chain micro-model. It showed **FB == Jacobi@∞** at the cliff
(exon dragged to the intron level); Jacobi only looked better in the benchmark because at its pass budget it was
**under-converged** (accidental regularization). It also derived a principled smooth fix — a **Student-t robust
message** (`w=(ν+1)/(ν+r²)`, `r²=(ρ_src−θ̂)²/v_msg`) — proven to recover the cliff *with* the enrichment prior and
be harmless on a uniform chain. **BUT** the real-data dissection (§5) then showed propagation is *not* the
dominant error — so the robust weight is deferred (it fixes the 4–19%, not the lever).

## 5. The node-by-node dissection (the core finding)

`scripts/debug/fb_node_dissect.py` — per chain REGION node, compares **true** (oracle gDNA fraction from the
split BAMs) vs **strand-only** local vs **local** (strand+global) vs **final** (after messages), decomposing the
gDNA-mass error into strand / global / propagation. Run on three scenarios (all post the mode-fix):

| scenario | strand | global ê(z) | **propagation** | total |
|---|---|---|---|---|
| flagship, unstranded (ss0.50) | 49% (−143K) | 46% (+230K) | **4% (−54K, correcting)** | +32K |
| flagship, stranded (ss0.99) | 46% (−95K) | 51% (+87K) | **3% (−4K)** | −12K |
| zero-gDNA + nascent (ss0.50) | 47% (+442K) | 34% (−170K) | **19% (−181K, correcting)** | +91K |

**Conclusions:**
- **Propagation is never the villain** — 4–19% of \|err\|, always net-*correcting* (it pulls the local's
  over-call back down; in zero-gDNA it pulls the phantom down −181K against the strand's +442K).
- **The strand is empty where it matters** — at κ≈½ (unstranded) the strand-only belief is a constant ~0.68 for
  *every* single-strand exon regardless of truth (count-zero-info: BB can't separate gDNA/RNA at κ=½); even
  stranded it under-calls (RNA's sense-skew masks gDNA's balance). It is not "wrong," it is *uninformative*.
- **The global ê(z) is the sole discriminator and it's biased + noisy — worst on AMBIG exons.** Proven case:
  zero-gDNA **node 2167** (AMBIG exon, true 0): strand 0.30 → **ê(z) pushes it to 0.80 (pure phantom)** →
  propagation drags it back to 0.23. The global *injects* the error; propagation *cleans up*.

**Precise offending nodes to track:** AMBIG exons where ê(z) mis-estimates (2167, 463, 2121, 935, 2179); the
13–15 single-strand exons isolated from any gDNA message (19, 1479, …) stuck on local.

## 6. The boundary-emission finding (user question: which boundaries emit no message, why, correct?)

`fb_node_dissect.py` boundary tally (flagship): **exon↔intergenic (TSS/TES) boundaries emit NO gDNA message —
0 of 202 — despite 345K of gDNA crossing mass.** exon-intron 931/940 and exon-exon 78/78 do emit. Cause: the
gDNA-emission gate is `solvable = (free_pos|free_neg)` — a **strand-continuity** test; at a gene boundary the
strand is discontinuous (exon↔intergenic) → not solvable → suppressed.

**Verdict:** the gate uses the *RNA* criterion (`free_s`) for the *gDNA* message — wrong in principle (gDNA is
strand-agnostic). The design intent (an intergenic *region* sink at init {0,0,1} must not broadcast "all gDNA")
is right, but it wrongly catches TSS/TES *boundaries* that hold real gDNA crossing data, and discards it. It
"works" under capture only by accident (the TSS is a depletion cliff there). Fixing it safely needs the robust
weight (else it relays the cliff) **and** solving the boundary's gDNA (currently a G1 sink) — deferred.

## 7. Bugs / cleanups

| item | status | effect |
|---|---|---|
| **Unbounded message mode** (gDNA→24,183, 46% >1; RNA→479,859, negative) | **FIXED** (clip to [0,1] in `node_sweep` `_scan`) | clean + prerequisite for weights; ~neutral on the benchmark (flagship flat, zero-gDNA +2%) |
| Marginal readout recomputed 6× | **FIXED** (`_node_marginals`, value-identical) | minor speed |
| Message precision can hit ~1.2M (honest `(M/E)²` Jacobian) | deferred | harmless unless the mode is a cliff → robust weight's job |
| gDNA-emission strand-gate (TSS/TES suppressed) | deferred | needs the robust weight + solving the boundary |

## 8. Git topology (as of 2026-06-27 cleanup)

- **`calib-spliced-junction-strand` (main working branch), HEAD 346f1531 = the committed Jacobi solver.** This is
  the validated production base. On top of 346f1531 this cleanup adds only: this doc, `fb_bridge_dissect.py`,
  `fb_prior_tabulation_audit.py` (both base-agnostic — work on the Jacobi base). **No calibration code changed.**
- **`calib-fb-dissection-jun27` (preserved snapshot, pushed).** Everything from the session: the FB inner solve
  (replaces Jacobi in `bp_solver.node_sweep`), the mode clip, the `_capture` hook, `_node_marginals` (dedup), the
  **ê(z) eff-len fix** (`fit_enrichment_transfer` → `eff_spl_left/right`), and the FB-internal harnesses
  (`fb_node_dissect`, `fb_ehat_dissect`, `fb_message_audit`, `fb_cliff_toy`, the gs/outer convergence probes).
  The FB goldens there are stale (regen if FB ever lands).

**To resume FB / the ê fix:** `git checkout calib-fb-dissection-jun27` (or cherry-pick `fit_enrichment_transfer`).
The ê fix is value-validated-correct but net-mixed standalone (see §0 / §11) — re-land it *with* the robust message.

## 9. Resume plan — closing the CALIBRATION gap (priority order)

The lever is the **local belief** (per the §5 dissection), not propagation. In rough ROI order:

1. **The capture cliff + the global ê(z) prior** — the dominant local error (32–46%). The ê eff-len bug is found
   (explore branch); the remaining work is the **robust-message weight** (`fb_cliff_toy.py` has the derivation +
   proof: `w=(ν+1)/(ν+r²)`) so ê can be de-biased *without* the cliff over-smoothing. This is the calibration-
   accuracy lever for both unstranded (the cliff) and zero-gDNA (the ê phantom).
2. **The AMBIG-exon capture leak** (the long-standing PROVEN #1 error — see memory
   `calibration_ambig_exon_capture_leak.md`): balanced-strand AMBIG exons hold the bulk of the under-called gDNA
   under capture; over-confident RNA imputation floods f_g→0. The robust message + the ê de-bias both bear on it.
3. **The gDNA-emission strand-gate** (§6) — TSS/TES boundaries discard real gDNA crossing mass; safe to open only
   with the robust weight + a solved boundary gDNA.
4. **FB commit decision** — only after the local belief is fixed; FB's value is masked until then.

**Separate track (bridge/EM, NOT calibration accuracy):** the §11 gDNA eff-len under-contraction — a confirmed
~200K lever in `priors.assemble_priors`. High value, but it is a downstream fix; keep it off the calibration-
accuracy critical path per the Framing note.

## 11. The calibration→EM bridge audit (NEW, 2026-06-27) — the ~200K stranded leak is the gDNA eff-len

Triggered by: the STRANDED flagship (`gdna_gdna300_ss_0.99_nrna_none_capture_on`) calibration prior is *good*
(contained error −16K) yet the post-EM benchmark leaks ~176K gDNA→RNA. Dissected the bridge
(`calibrate → build_multi_loci → assemble_priors → EM`), focused on `assemble_priors` + the IPR eff-len.

**Findings (`fb_bridge_dissect.py`, `fb_prior_tabulation_audit.py`):**

1. **Mass conservation PASSES** — `Σ gdna_region (2,916,453) = Σ contained (1,647,246) + Σ internal seam
   (1,269,207)` exactly. Boundary (seam) nodes are counted once, no double-count, no drop.
2. **The prior gDNA split is ACCURATE** — per-locus `prior_frac` vs `oracle_frac`: mean|Δ|=0.051, corr +0.78;
   EM-vs-oracle (0.045) is even tighter. `assemble_priors` is *not* mis-assembling the split.
3. **The EM systematically calls 6.7% LESS gDNA than the (accurate) prior** (signed `em_frac − prior_frac` =
   −0.067, every top-leak locus) → the leak is **downstream of a correct prior**.
4. **Root cause = the gDNA component eff-len is ~2× too large.** Decisive ablation: `gdna_eff_len ×0.5` →
   **leak 176,672 → 1,138 (−99%)**, `gdna_observed` 2,772,948 → 2,948,482 (oracle 2,949,620). gDNA was at a
   structural density disadvantage; halving the eff-len removes it.
5. **Mechanism (the boundary/IPR tabulation):** the gDNA eff-len is an IPR over contained-exon nodes
   (density 8.13) **+ boundary seam nodes (density 3.14 — half-probed crossings, 29% of the exon gDNA mass)**.
   The half-density seams **dilute** the IPR. Plus a gDNA↔RNA **asymmetry**: an RNA transcript's eff-len is the
   FL-marginal of its *contiguous spliced* length (no internal boundaries), while gDNA's IPR carries the genomic
   boundary penalty → gDNA's eff-len is systematically larger → it under-competes against RNA.

**Fix direction (bridge, not calibration):** make the gDNA eff-len reflect the exon-interior competition density
(where ~71% of its mass is), comparable to the RNA spliced eff-len — e.g. a contained-only IPR, with the seam
gDNA still counted in `gdna_prior_count` but not diluting the density. **Open design point:** how to treat the
seam (boundary-crossing) mass — it is real gDNA but competes at the boundaries, not the exon interior. No magic
constant (the ×0.5 was only the diagnostic probe). Harnesses: `fb_bridge_dissect.py [cond]`,
`fb_prior_tabulation_audit.py [cond]` (both run on the Jacobi base).

## 10. Instrumentation (all reproducible, OMP_NUM_THREADS=1)

**On main (Jacobi-base, base-agnostic):**
- `fb_bridge_dissect.py [cond]` — per-locus prior/oracle/EM gDNA-fraction comparison; localizes leak to bridge
  (prior wrong) vs downstream (EM/eff-len). The §11 tool.
- `fb_prior_tabulation_audit.py [cond]` — `assemble_priors` boundary/region mass + density tabulation
  (conservation, contained-vs-seam density, vs oracle). The §11 tool.

**On `calib-fb-dissection-jun27` (FB-internal — need the FB `_capture` hook):**
- `fb_node_dissect.py [cond]` — per-node strand/global/propagation decomposition + boundary emission table.
- `fb_ehat_dissect.py` — the ê(z) teacher/apply dissection (z, ê(z), apply mask, vs oracle).
- `fb_message_audit.py` — message defect scan (mode range, prec, NaN, emission gate).
- `fb_cliff_toy.py` — the standalone cliff micro-model + the robust-message proof.
- `node_sweep(..., _capture=dict)` — the inert hook recording per-node local/α/β/final/var~mean.

**Benchmark:** `scripts/sim/evaluate_suite.py --sim-base ~/Downloads/rigel_runs/quick_3to1_5mb` +
`scripts/debug/precision_benchmark_report.py <suite> <baseline-suffix>`.
