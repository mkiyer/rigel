# Session handoff — the DNA-background reference (Phases 0–3 shipped, Phase 4 next)

*Written 2026-07-17 to let a fresh session pick up mid-stream. Branch: `calib-ambig-init-wip`.*

---

## 0. First actions (do these before anything else)

1. **Read the memory** (loaded via `MEMORY.md`), in this order:
   - `background_reference_derivation.md` — **the current thread** (derivation + implementation progress P0–P3).
   - `refit_loop_study_findings.md` — the flagship dissection + refit study (why we're doing this).
   - `sigma2_transfer_derivation_state.md` — the σ²_transfer F1 work (the prior thread, reviewer-signed-off).
   - `calibration_architecture_count_zero_info.md`, `calibration_roadmap_state.md` — the invariants.
2. **Read the docs** `docs/calibration/background_reference_derivation.md` (the theory, reviewer-confirmed) and
   `docs/calibration/background_reference_implementation_plan.md` (the phased plan — P4/P5 remain).
3. **Activate the env** — everything runs inside the `rigel` conda env:
   `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`. Always `OMP_NUM_THREADS=1`
   for deterministic calibration debug runs.

---

## 1. The problem

Rigel deconvolves each genomic node's unspliced fragment mass into `(f_pos, f_neg, f_g)` — sense-RNA /
antisense-RNA / gDNA — via a region↔boundary belief-propagation sweep (`calibration/bp_solver.py`). The
governing invariant (`CALIBRATION_ARCHITECTURE.md`): **a fragment count carries ZERO intrinsic gDNA-vs-RNA
information**; composition comes only from the strand likelihood, cross-node messages, and the prior.

**The north-star diagnostic is the UNSTRANDED + CAPTURE corner** (`ss_0.50 … capture_on`). It is not a common
real library, but it is where every weakness shows: strand carries no information, so the prior + messages must
do everything. "Make the flagship work and everything else follows" is a proven approach here. **Never test the
flagship in isolation** — the AMBIG stress suite must include ample single-stranded nodes.

### What this session found (the chain that led here)
- The σ²_transfer message-damping (F1) is correct and helps stranded data hugely, but **hurts the flagship**
  (it gags the messages that are the flagship's only signal).
- Dissecting the flagship (`scripts/debug/flagship_dissect.py`): the "extremely weak" gDNA prior actually
  **doubles** the error (strand-only 0.33 → +prior 0.62); the messages *help* 34:1. Root cause: an **asymmetric
  ψ** — the gDNA NPMLE prior replaced the Jeffreys `½log f_g` arm but **clamped its left tail flat**, so it lost
  the `f_g→0` barrier while the RNA arm kept its `f_g→1` barrier. Strand-blind gDNA nodes get crushed to `f_g≈0`.
- The naive fix (V4: restore a symmetric Jeffreys barrier) was **rejected** by the full 24-scenario A/B — it
  fixes the flagship but manufactures 40% false gDNA in zero-gDNA unstranded libraries. **No fixed prior works**
  for both corners; the gDNA *level* is a measured data quantity, not a belief.

---

## 2. The solution (derivation — settled, two external reviewers confirmed)

`docs/calibration/background_reference_derivation.md`. The DNA level is a **genome-wide scalar**, not a
per-region quantity: a region of effective length `E` resolves DNA only above `~1/E`, so a faint background
(the strong-capture case) is resolvable ONLY by pooling regions into one aggregate support. Three parts:

1. **`ρ_bg` — the aggregate background** = `Σg/ΣE` over the pooled **intergenic** regions (signature 0; pure DNA,
   no mature/nascent RNA), with `σ_bg = 1/√Σg`. Precise even at true zero (`Σg=0 ⇒ ρ_bg=0` confidently).
2. **The one-sided floor** — enters the `f_g` solve as `−½·(max(0, log ρ_bg − log ρ_g)/σ_floor)²`. Data-gated
   (dormant when `ρ_bg→0`, so a DNA-free library is never floored), `σ_bg`-softened, and **NEVER a denominator
   or multiplicative scale** (the strong-capture safeguard — off-target `ρ_bg` sinks toward 0 under strong
   capture even when DNA is present, so it must be a soft one-sided *lower bound* that recedes gracefully).
3. **The pinned component (formula 2)** — the NPMLE mixture gains one component pinned at `log ρ_bg`, so the
   sub-resolution mass is placed correctly instead of a Kiefer–Wolfowitz atom at 0 (no double-counting).
   Reviewer-confirmed that `bg + NPMLE` (formula 1, the released tool) double-counts; `bg + (NPMLE − bg)`
   (formula 2) is correct.

**The honest limit (important):** the off-target floor CANNOT lift the enriched *on-target* gDNA nodes (their
`ρ_bg/ρ_tot` is tiny) — that DNA lives on-target, entangled with RNA, so it needs an **iterative global-abundance
estimate** anchored by the floor. **That is Phase 4** and is where the flagship actually improves.

**Region-set caveat (owner directive):** the sim's `nrna_present = 0.5×mature` is unrealistically abundant, so
introns are contaminated → we pool **intergenic-only** for sim development. Real nascent is *sparse*, so on real
data pool **intergenic + intron** (`include_introns=True`, huge aggregate-span/resolution gain) with an optional
`robust_trim_mad` fence for the rare contaminated intron. This is an option to explore **on real data, not the
sim**.

---

## 3. What is implemented and WIRED (Phases 0–3, production)

All in the `rigel` calibration package. The background floor + pinned component are **live in production by
default** (`background_floor=True`).

- **`src/rigel/calibration/background_reference.py`** (NEW) — `measure_background(substrate, region_arrays,
  region_eff_len, *, include_introns=False, robust_trim_mad=None) -> BackgroundReference(log_rho_bg, sigma_bg,
  n_counts, eff_total)`.
- **`src/rigel/calibration/gdna_rate_prior.py`** — `GdnaRatePrior` gained `log_rho_bg`/`sigma_bg` fields (default
  `-inf`/`inf` ⇒ dormant), `fit(background=BackgroundReference|None)` (adds the pinned EM column + stores the
  floor fields), and `logprior()` adds the one-sided floor with **per-node softening**
  `σ_floor² = σ_bg² + 1/(ρ_bg·E)` (folds in each node's own gDNA Poisson spread — without it the floor was a
  near-hard wall under capOFF that over-called nascent nodes).
- **`src/rigel/calibration/calibrate.py`** — measures the background once and passes it to every
  `GdnaRatePrior.fit` (pass-0 + refits).
- **`src/rigel/config.py`** — `CalibrationConfig.background_floor` (True), `background_include_introns` (False),
  `background_robust_trim_mad` (None, validated >0).
- **P0 cleanup:** deleted `message_precision.py` + 3 broken KDE debug scripts; relocated the retired σ²_imp
  helpers to `scripts/debug/_disagreement_variance.py`; fixed stale comments.

**Tests (all green except one known-red):**
`tests/calibration/test_background_reference.py` (8), `tests/calibration/test_gdna_rate_prior.py` (16, incl.
floor one-sidedness / dormancy / per-node-softening / pinned-component). Full calibration suite: **229 pass**.
The **only** non-golden failure is the *pre-existing* `test_bp_solver.py::test_gdna_sweep_zero_gdna_pin_and_monotone`
(fails at 0.6833; confirmed inert/independent this session). The **22 golden failures are pre-existing/stale**
(from the uncommitted σ²_transfer wiring) — **golden regen is deferred to ship-time** (owner: "after a complete,
correct, bug-free, validated implementation").

---

## 4. Current results (generated this session)

- **Full 36-scenario A/B** (`scripts/debug/test_floor_ab.py`, floor-off vs floor-on): gdna_none **+0.0000 at
  every capture** (dormant, no false positives); off-target/DNA-heavy conditions improve (up to −0.022); tiny
  residual on nrna-present capOFF (+0.008 worst); flagship gently improved. Group means: gNONE +0.0000, g1
  +0.0007, g5 +0.0027, g100 −0.0008, g300 −0.0042. **Safe and roughly neutral standalone** — the payoff is P4.
- **Prior-vs-oracle plots** (`scripts/debug/background_prior_plot.py` → `background_prior_plot.png`): the combined
  prior (NPMLE + pinned) vs oracle gDNA density across 6 conditions. Shows: the anchor sits on the depleted peak;
  the pass-0 prior is total-based and over-represents the high-ρ tail (RNA contamination the refit deconvolves);
  the pinned component visibly fixes the Kiefer–Wolfowitz atom (gdna1 verystrong); gdna_none floor dormant.
- **Background validation** (`scripts/debug/background_reference_probe.py`): `ρ_bg` per condition + the
  resolution argument (100% of intergenic regions have count 0 for gdna_none, yet the aggregate resolves it).

---

## 5. NEXT: Phase 4 — the iterative global-abundance lift (the flagship fix)

The off-target floor doesn't lift enriched on-target gDNA; the *global DNA abundance* must be estimated
iteratively (the refit), **anchored by the floor**. Build it carefully — the refit-loop study
(`refit_loop_study_findings.md`) found the naive refit is **net-harmful and can oscillate/drift**:
- The **drift** = prior-sharpening self-training (present cold or warm).
- The **oscillation** = the warm-continue; **cold-restart eliminates it**.
- The **degeneracy**: minimizing between-node disagreement has a trivial over-smoothed minimum; the constraints
  (fixed total, fixed spliced, strand, and now the floor) must pin the correct fixed point.

Phase 4 = each pass re-estimates the global DNA fraction from the *solved on-target* DNA, lifts the enriched
nodes, with the floor as the fixed lower bound they can't be crushed through; use cold-restart +
confidence-weighting; gate on the disagreement-vs-oracle-floor instrument (built into `refit_loop_study.py`):
**the solved between-node disagreement should approach the ORACLE's, not zero** (overshooting below = degenerate
over-smoothing). Then **Phase 5**: golden regen (`pytest tests/ --update-golden`) + rewrite
`CALIBRATION_ARCHITECTURE.md` / `calibration_prior_production_reference.md` to the shipped design.

Parallel future (refit-gated): write **`logP_r`** (the symmetric RNA arm) to retire the lopsided Jeffreys RNA
reference and close the ψ asymmetry at both ends.

---

## 6. Where everything lives

**Docs (`docs/calibration/`):** `background_reference_derivation.md` (theory + strong-capture safeguard + the
formula-2 §8.5 + the validation §5); `background_reference_implementation_plan.md` (phases, code-health ledger);
`flagship_prior_asymmetry_diagnosis.md` (the root-cause dissection + V4 rejection); `transfer_variance_formal_derivation.md`
(σ²_transfer F1 §10 reviewer sign-off); `refit_loop_study_findings.md`; `CLEANUP_LOG.md` (deferred tidy-ups);
`CALIBRATION_ARCHITECTURE.md` (authoritative count-zero-info).

**Production code (`src/rigel/calibration/`):** `background_reference.py`, `gdna_rate_prior.py`, `calibrate.py`,
`bp_solver.py`, `simplex_logodds.py` (the ψ / `_gdna_arm` / `_rna_arm`), `density_model.py`
(`count_observable_masks`), `signature.py` (`BIT_*`, region types), `node_geometry.py`, `substrate.py`,
`effective_length.py`. Config: `src/rigel/config.py` (`CalibrationConfig`).

**Tests (`tests/calibration/`):** `test_background_reference.py`, `test_gdna_rate_prior.py`, `test_bp_solver.py`.
Golden tests: `tests/test_golden_output.py` (+ `tests/golden/`) — stale, regen in P5.

**Debug tools (`scripts/debug/`):** `test_floor_ab.py` (36-scenario floor A/B), `background_prior_plot.py`
(prior vs oracle), `background_reference_probe.py` (ρ_bg validation), `refit_loop_study.py` (the refit harness +
disagreement instrument — **P4's main tool**), `flagship_dissect.py` + `prior_crush_probe.py` (the dissection),
`test_v4_symmetric.py` (the rejected V4, kept for reference), `npmle_movie_frame.py` (deconv-vs-oracle frames),
`selfsolve_diag.py` (`_scan_and_truth`, the cached-scenario loader every tool uses).

**Scenarios / caches (`~/Downloads/rigel_runs/ambig_dense_10mb/`):** 36 conditions
(gdna_none/100/300 × ss × nrna × capture = 24, **plus** the 12 low-DNA/capture-strength added this session:
gdna1/gdna5 × {off,on,verystrong}, unstranded+nascent). New sim config:
`scripts/sim/configs/lowdna_capture_strength.yaml`. Self-solve cache: `_selfsolve_cache/*.pkl` (all 36 cached;
tools read `sim_oracle.bam`, no re-scan needed). Reference/index at `.../rigel_index` + `.../reference`.

---

## 7. Standing directives & gotchas

- **NO new magic numbers / heuristics without discussing first.** Calibration is capped (~25 modules, few
  constants). `ρ_bg`/`σ_bg` = pooled Poisson MLE; floor = `½z²`; MAD→σ `1.4826` is mathematical. The only tunable
  knob is `bandwidth` (a resolution knob) and the new `robust_trim_mad` fence (explicit, real-data-only).
- **Develop on the cached scenarios**, `OMP_NUM_THREADS=1` (cross-process nondeterminism otherwise).
- **The flagship (unstranded+capture) is the diagnostic**; production accuracy on it is NOT yet the ship target.
- **Golden regen is deferred to P5** — until then the 22 golden fails and the one 0.6833 red are expected; judge
  regressions by the 36-scenario A/B, not the goldens.
- **A concurrent editor has touched `bp_solver.py`/`calibrate.py`/`gdna_rate_prior.py` mid-session** — pin md5
  hashes before/after any A/B (post-P3: `bp_solver 871cdf42…` at P0; re-check current).
- **Sim caveat:** intergenic-only for the sim (abundant nascent); `+intron` is the real-data path.
- **`git status` is dirty** (uncommitted σ²_transfer wiring + this session's work + modified goldens). Nothing
  committed this session — the owner drives commits.
- Two external statistical reviewers have signed off on both the σ²_transfer F1 assembly and the background
  reference formulation — the theory is on solid ground; the remaining work is engineering (P4/P5).
