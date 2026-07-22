# Next-session prompt — fix the splice-junction message-propagation bug

Paste the block below to start the next session. It is self-contained; the referenced docs carry the detail.

---

You are continuing calibration work on Rigel, branch `calib-ambig-init-wip`. Your task this session is to fix a
**diagnosed, owner-confirmed message-propagation bug** in the calibration solver. Do the reading first, then
walk the worked example to re-confirm the fix on a live node, then implement, then A/B.

## Standing constraints (do not violate)
- **NO new magic numbers / constants / heuristics without pausing to discuss first.** This is a hard rule.
- **No clamps, cliffs, or binary cutoffs.** Prefer smooth, honest Bayesian behaviour. A clamp that manufactures
  confidence is exactly the bug you are fixing.
- **Nascent RNA cannot be inferred from abundance.** There is no biological correlation between expression level
  and nascent fraction. We can **only measure** (spliced = mature at a boundary; unspliced crossing = nascent +
  gDNA). Any "more RNA ⇒ more nascent" reasoning is wrong.
- Build/test/lint **only** inside the activated `rigel` conda env
  (`source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`). Always `OMP_NUM_THREADS=1`.
- **Do not commit** — the owner drives commits. Golden regeneration is deferred to ship-time.
- Develop on the cached toy/benchmark scenarios, never on big suites cold. Always test AMBIG with **ample
  single-strand nodes present** — never AMBIG in isolation.

## Read these first (in order)
1. `docs/calibration/message_absorption_fix.md` — **the work order.** The bug, the worked node example (real
   oracle counts), the owner-confirmed correct model, and the ordered implementation plan (Steps 0–3).
2. `docs/calibration/CALIBRATION_MASTER.md` §2 (count-zero-info), §3 (the message-content rule), §8 (the gap
   table — the two ⛔ rows are your task), §9 (Phase 1b).
3. `docs/calibration/CALIBRATION_ARCHITECTURE.md` §0–§1 — the invariant.

## The bug in one sentence
At a splice-junction boundary, `bp_solver._scan` subtracts the **destination boundary's spliced count** from a
**neighbour's unspliced-RNA density imputation** (`rho_pos -= SPd[i]/ESPd[i]`); these are disjoint pools at
different scales, so it goes negative, and `max(rho, 1/erd)` launders that into a **confident "no RNA"** message
→ the node snaps to f_g≈0.97 and cascades. This is the seed of the confident-false-positive population (420/666
in gdna5) that also corrupts the Phase-2 gDNA-prior fit.

## The correct model (owner-confirmed — implement this)
- A **region** reports only the densities it *has* `(ρ_RNA+, ρ_RNA−, ρ_gDNA)` — **no** spliced subtraction.
- A **boundary** owns splicing: it partitions its own directly-measured spliced (mature) vs unspliced (nascent +
  gDNA) on its **own** effective length; the subtraction lives on the boundary's **outgoing** side (nascent
  continues; spliced does not).
- A cancelled / non-positive imputation carries **no information** ⇒ emit **no message** (precision 0), never a
  confident value.

## Do this, in order
1. **Re-confirm on a live node before touching production code.** Using
   `scripts/debug/confident_fp_trace.py --condition gdna_gdna5_ss_0.50_nrna_present_capture_on --nodes 2930,2929,2931`
   and the oracle pools (`_scan_and_truth(...)['boundary_pools'/'region_pools']`), walk node 2930 (INTRON 2929 —
   boundary 2930 — EXON 2931) through the corrected flow and show it lands at **f_g ≈ 0.007** (nascent ≈ 955,
   gDNA ≈ 7 of 962 unspliced). Nail the **eff-length basis** here (`message_absorption_fix.md` §5) — the incoming
   RNA density and the boundary's spliced count must be reconciled to the same support so the subtraction is
   exact on this node AND degenerate-safe (no spliced ⇒ no change; pure intron ⇒ all-nascent). This is geometry,
   not a knob — if you find yourself wanting a constant, STOP and discuss.
2. **Step 0 — honest clamp.** In `bp_solver._scan`, when the imputed RNA density is ≤ 0 (or ≤ `1/erd`), emit no
   RNA factor (precision 0), do not clamp-and-send. Land + A/B this alone first.
3. **Step 1 — remove the cross-node spliced subtraction** (`absorb_p`/`absorb_n`, gated by the existing
   `_RNA_ABSORB` toggle — make removal permanent once verified). Re-examine the `+ SPs[lsrc]/esp` source-spliced
   addition at the same time.
4. **Step 2 — move the spliced partition to the boundary's own side** in consistent eff-length units (§4–§5 of
   the fix doc). Verify the acceptance test.
5. **A/B every step** on the flagship (unstranded+capture) **plus** zero-gDNA **plus** a stranded control — never
   a single aggregate. Tools: `confident_fp_trace.py`, `npmle_smoothness_diag.py`, `hyperprior_fit_options.py`.
   Record the confident-FP mass fraction before/after. Run `pytest tests/calibration/ -q` (expect 1 pre-existing
   known-red: `test_gdna_sweep_zero_gdna_pin_and_monotone`).
6. When it lands, update `CALIBRATION_MASTER.md` §8 (flip the two ⛔ rows) and note the A/B result; then Phase 2
   (the deconvolved-gDNA hyperprior, Role B) is unblocked.

## Key facts / paths
- Node 2930 oracle: unspliced 962 (7 gDNA + 955 nascent) + 7502 mature-spliced ⇒ f_g=0.007. EXON 2931 unspliced
  119021 (493 gDNA + 73639 mature-unspl + 44889 nascent). INTRON 2929: 78 all-nascent.
- `bp_solver._RNA_ABSORB` (module toggle, default True) = the absorption A/B switch already in place.
- Caches: `~/Downloads/rigel_runs/ambig_dense_10mb/_selfsolve_cache`; suite `ambig_dense_10mb` (24 conditions).
- The `_scan` RNA block is ~L430–455 of `src/rigel/calibration/bp_solver.py`.

Work meticulously and show the owner the worked node numbers before and after each change. Ask before adding any
constant.
