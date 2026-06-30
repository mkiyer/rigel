# Calibration — next-session handoff (tasks ahead)

**Date:** 2026-06-30. **Branch for WIP:** `calib-disagreement-precision`. **Read first:** the project memory
`calib_mature_message_and_toy_prod_driver.md` and `disagreement_precision_shipped_flagship_ehat_bound.md`.

---

## 0. What is on `main` (shipped) vs. WIP

**Shipped to `origin/main`:** the **disagreement-aware message precision** change — retired the var~mean
message-reliability layer; message precision is now per-edge `σ²_edge = max(resid² − (base_var + var_loc[dst]),0)`
against the dst's message-free local anchor; the outer var~mean loop collapsed to one forward-backward pass; the
global prior uses the anchored seed σ²_g. Validated on the 16-condition benchmark (mean MARD 12.23→12.07,
total |gdna_net_surplus| −5.9%; unstranded-capture leak improved most; `nrna_rnd`-unstranded regressed — an
accepted soft spot). Design: `docs/calibration/dispersion_aware_message_precision.md`.

**WIP on the branch (NOT shipped):** the **spliced→mature MEASUREMENT message** in `bp_solver._scan` (see §3).
It works directionally but has two open follow-ups (§4) and its goldens/benchmark are not done.

---

## 1. The design we are building toward (the model)

Under hybrid capture the gDNA prior is a **mixture** (depleted off-target + enriched on-target populations;
intrinsic CNV would add more). Frame this as **soft marginalization over a NONPARAMETRIC empirical prior** built
from the *solved single-strand nodes' gDNA densities* — NOT hard clustering (which is brittle: a misclassified
node snaps to the wrong mode). The per-node grid solve is already a correct product-of-experts; the prior's
strength self-scales by the node's own evidence via the `(2κ−1)²` tilt-Fisher information (→0 at a balanced
ridge AND at κ→½ AND off-capture), so it fills exactly the likelihood's null space with no detector and no
mode-counting. Off-capture the empirical density is unimodal; capture → bimodal; CNV → multimodal — same code.

**Priority #1 is a pristine single-strand solver** (structure + strand + message-propagation + an *extremely
weak* gDNA floor for stability only). The empirical prior is bootstrapped from it, so it must be correct —
*stranded AND unstranded*. Then: train the empirical prior from the single-strand solutions → solve AMBIG nodes.

**Message propagation IS the generalized "spliced residual."** The missing piece was the **precision** model.
Two message KINDS:
- **IMPUTATION** (gDNA, nascent) — genomic/regime smoothness → keep the **disagreement-aware** precision.
- **MEASUREMENT** (mature, from the splice junction) — a gapped read *is* mature RNA, directly observed, and
  mature coverage is uniform along the transcript (same molecule as the exon body) → **count precision**, fused
  with (not bypassing) the strand likelihood, so a confident strand wins and a strand-blind node **snaps** to it.

---

## 2. ⚠️ THE NASCENT-RNA PRINCIPLE (hard constraint — do not violate)

**We must NOT hallucinate nascent RNA from a discrepancy between a boundary and its neighbouring exon.** The
intron over-call we currently see (a pure-gDNA intron solved at f_g 0.54–0.70 instead of 1.0) is exactly this
hallucination — the exon's mature RNA leaking into the intron through the message system.

**Definitive evidence of nascent RNA must come from the INTRON ITSELF — specifically a SINGLE-STRANDED intron —
via one of two signals:**
1. **Density excess:** the intron's fragment density is *much higher* than the intergenic gDNA density (the true
   depleted floor). The excess over the floor is nascent.
2. **Strand:** the strand model shows RNA (sense bias) inside the intron.

Rules that follow:
- Intronic fragments *support* nascent RNA; boundary/exon messages do **not** create it.
- If message passing produces **very scant** hallucinated nascent, that is acceptable; **gross** hallucination
  is not.
- The "floor" for these comparisons is the **truly depleted** intergenic density (a *minimum*-class baseline),
  **not** the population mean. (This also fixes the earlier ê/global "mean-not-minimum" issue.)

This principle drives follow-up B (§4) and the eventual nascent model.

---

## 3. The mature-message WIP (what was implemented this session)

In `src/rigel/calibration/bp_solver.py`, `node_sweep._scan`, blocks `emit_p`/`emit_n`:
- The source RNA message density is now **nascent (`fbp·sm/er`) + mature (`SPs/ESPs`)**, where `SPs` is the
  one-sided spliced on the source face and `ESPs = eff_spl`. Because the spliced geometry is **one-sided** (zero
  on the non-exon face), the mature term routes to the **exon flank only** — never into an intron (automatic).
- When mature is present (`n_mat>0`): **count precision** `1/(vbp + 1/(n_nasc+n_mat))`, **no** disagreement
  silencing (it is a measurement). Pure nascent: keep the disagreement-aware path.
- Re-added the `ESP = (eff_spl_left, eff_spl_right)` geometry tuple in `node_sweep`.

**Validated on the production-faithful toy** (`scripts/debug/toy_prod.py`, TA+ = `(1000,2000),(5000,10000)`):
unstranded-capture exon leak `−0.088 → −0.022`; unstranded over-call `+0.10 → −0.04`; **stranded unchanged**
(no regression); 188 calibration unit tests pass.

---

## 4. Concrete tasks, in order

**(A) Fix the mature projection geometry (`eff_spl` / O2) — PREREQUISITE.** The mature density uses the spliced
half-triangle `eff_spl`, which is ~2× too large for a big exon, so `ρ_mature` is under-counted and the big exon
over-corrects (TA+ R3: truth 0.32 → solved 0.46, `err +0.134`). The measurement is only as good as its eff-len.
Fix the mature projection in `effective_length.py` / `node_densities` / the `_scan` mature term so `ρ_mature`
projects correctly onto the exon's contained RNA eff-len (target: R3 err → ~0, R1 stays ~0).

**(B) Nascent only from intron self-evidence (THE PRINCIPLE, §2).** The intron over-call is the exon's mature
leaking via the *nascent* message (`n_nasc = fbp·sm` carries the whole RNA fraction). First step: impute
**nascent = total RNA − mature**, so a pure-mature exon imputes ~0 nascent into the intron. Then enforce the
principle: an intron's nascent must be supported by the intron's own **density excess over the intergenic floor**
or its **strand sense-bias** — messages must not manufacture it. Target: TA+ intron (nascent=0) → f_g ≈ 1.0.

**(C) Validate + ship.** Re-run the TA+ knob *sweep* (gDNA density × κ gradient 0.5→0.99 × capture-strength
gradient off→mild→extreme) on `toy_prod.py`; confirm exon→truth across the grid, intron→1.0, stranded stable.
Then the 16-condition benchmark (skill `calibration-benchmark`), regen goldens (`pytest tests/ --update-golden`),
confirm `pytest tests/` green. Only ship if the benchmark is net-positive and the `nrna_rnd` regression doesn't
worsen.

**(D) Then the empirical prior + AMBIG (the bigger arc).** Build the nonparametric empirical gDNA-density prior
from solved single-strand nodes (one robust bandwidth knob — provide a bandwidth estimator + a plotting
framework to inspect it; decide the estimator on real data later). Then solve AMBIG nodes by soft marginalization
over it. Grow the toys: partial-overlap and fully-encompassing opposite-strand AMBIG, microexons (regions < a
fragment → zero contained fragments → message-only), multi-gene neighbourhoods.

---

## 5. Tooling & caveats

- **`scripts/debug/toy_prod.py`** — the production-faithful toy driver (full simulator → oracle BAM → real
  `calibrate()` → per-region oracle-truth vs solved). `run(name, genes, kappa=, n_rna=, gdna_fraction=,
  capture=, capture_strength=, nascent=, instrument=)`. Develop against THIS (it *is* production); `instrument=True`
  dumps per-node anchor + messages.
- **DO NOT trust the old hand-rolled harness** (`toy_harness.py`/`toy_validate.py`): it hand-assembled
  `node_sweep` (od=0, clamped κ) and was proven **unfaithful** (flipped stranded RNA→gDNA). Always validate
  against real `calibrate()`.
- **The earlier flagship node-dives** (`flagship_node_dive.py`, the 0.533 anchor / 77-23 split, the predictor
  saturation) used that same hand-assembly → treat their exact numbers as **directional only**. The toy driver
  + the 16-condition benchmark (real `run_pipeline`) are the trusted substrate.
- Benchmark suites with BAMs: `~/Downloads/rigel_runs/quick_3to1_5mb` (16 conditions, used here);
  `~/Downloads/rigel_runs/gdna_benchmark_5mb` (20 conditions, BAMs were cleaned — regenerate via
  `scripts/sim/simulate_suite.py` if needed).
- Every build/test/run inside the activated `rigel` conda env.
