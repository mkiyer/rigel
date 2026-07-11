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

- **Forward-backward IS the production solver, shipped on `calib-spliced-junction-strand`.** FB replaces the
  per-node Jacobi inner sweep in `bp_solver.node_sweep` (forward α + backward β, true tree-BP, single inner pass;
  the outer loop iterates the var~mean fixed point). ~10× faster than Jacobi. Goldens regenerated; suite green.
  Shipped *with* the ê(z) eff-len fix + the mode clip + the `_node_marginals` dedup. The dead Jacobi message
  helpers (`_message`/`_message_vec`) are removed (FB inlines the message in `_scan`).
- **The node-by-node dissection (the core finding):** message propagation is only 4–19% of the error and is
  net-*correcting*. The error lives in the **LOCAL belief** — the strand likelihood (empty at κ≈½) + the **global
  ê(z) prior** (biased/noisy, worst on AMBIG exons). So the next lever is the local belief, not the solver.
- **ê(z) eff-len bug FIXED (shipped):** `fit_enrichment_transfer` projected ρ_mature with the *full* RNA eff-len
  instead of the half-triangle `eff_spl` (the spliced-fix da777bc7 corrected every other path and missed this one,
  so 346f1531 shipped it live). Fixed → cuts the zero-gDNA ê phantom ~85%; it also **un-masks the capture cliff**
  (a fragile cancellation where the buggy over-call hid the cliff under-call) → the cliff is now the visible next
  lever (the robust message), not a hidden number.
- **Benchmark honesty:** on `quick_3to1_5mb`, FB's Σ|leak| was ~+11% vs the committed Jacobi — but Jacobi's lower
  number is *accidental under-convergence regularization* (FB == Jacobi@∞, proven in `fb_cliff_toy.py`); the gap is
  the capture cliff, the real next target. (Re-benchmark FB+ê for the current number — §3 predates the ê fix.)
- **§11 — the calibration→EM bridge audit (separate track):** `assemble_priors` mass is conserved + the prior
  split is accurate, but the **gDNA component eff-len is ~2× too large** → ~200K stranded-flagship leak (×0.5 →
  176K→1K). A confirmed high-value **bridge/EM** fix, off the calibration-accuracy critical path (see Framing).
- **Next lever:** the **capture cliff / robust message** on the FB substrate (the local-belief fix for unstranded +
  zero-gDNA), then the AMBIG-exon capture leak. The bridge eff-len is a separate confirmed EM win.

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

## 8. Git topology (as of 2026-06-27)

- **`calib-spliced-junction-strand` (production branch) = FB SHIPPED.** `bp_solver.node_sweep` is forward-backward;
  the ê(z) eff-len fix, the mode clip, and `_node_marginals` (dedup) are in; dead Jacobi `_message`/`_message_vec`
  removed; goldens regenerated; suite green. Harnesses on main: `fb_bridge_dissect`, `fb_prior_tabulation_audit`
  (base-agnostic) + `fb_node_dissect`, `fb_ehat_dissect`, `fb_message_audit`, `fb_cliff_toy`, the gs/outer probes.
- **`calib-fb-dissection-jun27`** — the session-snapshot branch (now redundant, FB is on main); safe to delete.
- **Vestigial:** `node_sweep`'s `max_passes`/`convergence_delta` params are unused by FB (single inner pass; the
  outer loop uses `max_outer`/`outer_convergence_delta`) — retained-for-API, candidate for a follow-up cleanup.

## 9. Resume plan — closing the CALIBRATION gap (priority order)

FB is the substrate; the lever is the **local belief** (per the §5 dissection), not the solver. In ROI order:

1. **The capture cliff + the global ê(z) prior** (the dominant local error, 32–46%) — the ê eff-len bug is now
   fixed; the remaining work is the **robust-message weight** (`fb_cliff_toy.py` has the derivation + proof:
   `w=(ν+1)/(ν+r²)`) so the (now-visible) cliff is not over-smoothed. The lever for both unstranded (the cliff) and
   zero-gDNA (the ê phantom).
2. **The AMBIG-exon capture leak** (long-standing PROVEN #1 error, memory `calibration_ambig_exon_capture_leak`):
   balanced-strand AMBIG exons hold the bulk of the under-called gDNA; over-confident RNA imputation floods f_g→0.
3. **The gDNA-emission strand-gate** (§6) — TSS/TES boundaries discard real gDNA crossing mass; open it only with
   the robust weight + a solved boundary gDNA.
4. **Re-benchmark FB+ê** for the current Σ|leak| (the §3 table predates the ê fix); confirm the cliff is the gap.
5. **Cleanup follow-up:** retire the vestigial `max_passes`/`convergence_delta` (§8); delete the redundant
   `calib-fb-dissection-jun27` branch.

**Separate track (bridge/EM, NOT calibration accuracy):** the §11 gDNA eff-len under-contraction — a confirmed
~200K lever in `priors.assemble_priors`. High value, downstream; keep it off the calibration-accuracy critical path.

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
