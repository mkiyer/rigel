# Implementation plan — restore the additive KDE representation for the gDNA composition prior (Role B)

**Status:** design + phased implementation plan (2026-07-19). Supersedes the A1/A2/A3/A4 framing
in `CALIBRATION_MASTER.md` §9 with a precise, code-grounded plan. Companion:
`kde_vs_npmle_enriched_mode` (why the shipped KDE keeps the enriched mode), `gdna_background_floor_derivation.md`
(the floor location ρ_floor), `background_reference_derivation.md` (ρ_bg).

> ## ⚠ UPDATE (2026-07-19, later): read [`CALIBRATION_STATUS.md`](CALIBRATION_STATUS.md) FIRST.
> The **pass-0 background anchor (Phase A½ below) has been REVERTED** — pass-0 is prior-free again. The **root
> cause was traced upstream**: pass-0 does not subtract mature RNA at single-strand exons (the un-landed "nascent
> factory," `bp_solver.py:442-446`); the anchor was a downstream patch. **Next work = the nascent factory**, not
> the anchor. The additive-KDE (Role B) content below is still valid but its value hinges on a clean pass-0.
> The bookmark and Phase A½ below are retained for provenance but are superseded by `CALIBRATION_STATUS.md`.
>
> ## 🔖 RESUME BOOKMARK (2026-07-19) — superseded, see above
> **Where we are:** **Phase A0 + A1 are LANDED** (`config.gdna_prior_additive`, default OFF, byte-identical
> verified; 234 tests pass). The additive KDE **works for recovery**: on `gdna5_capON` the hyperprior's enriched
> mode goes **0.176 (EM, starved) → 0.874 (additive)**, 35-node recovery **0.178 → 0.357** (oracle 0.429),
> stranded `gdna300·ss0.99` no-regression (0.087→0.085). **BUT the specificity gate FAILS** (OI-1 confirmed):
> `gdna_none` RNA-node median f_g **0.054 → 0.258**. Root: the additive KDE faithfully represents pass-0's gDNA
> **over-call** (gdna_none pass-0 median ≈ 0.49) that the EM's heavy `n_regions` background used to crush.
>
> **DIGRESSION IN PROGRESS:** diagnosing the **PASS-0 (first-pass) zero-DNA false positives** — find the
> highest-FP zero-DNA scenarios (one nascent-present, one nascent-absent), rank nodes by FP contribution, and
> trace the worst nodes from **initialization → raw counts → neighbour counts → message construction → where the
> node lands** (pass-0 only). **Do NOT dig into the second pass / NPMLE refit** — that is mid-development. Goal:
> first-pass FP mitigation (fix the over-call at the source) so the additive KDE has nothing to amplify.
>
> **The pending decision (after the diagnosis):** OI-1 — the RNA-parsimony reference strength ½ vs 1 — **and/or**
> the message-absorption correctness fix, whichever the diagnosis shows is the root of the pass-0 over-call.
>
> **RESUME AFTER THE DIGRESSION:** land the pass-0 FP fix → re-run the A2 gates (recovery / FP / stranded) →
> if FP now passes, flip `gdna_prior_additive` default ON → A4 (data-driven bandwidth) → Phase B (real data).
> The additive-KDE code is DONE and safe; only the pass-0 FP + the OI-1 call stand between here and A2 default-on.

---

## Phase A½ — Pass-0 background DNA (ADOPTED as production; tuning open)

**Decision (owner, 2026-07-19):** the initial (pass-0) solve WILL include a background-DNA distribution from
intergenic regions — no longer reference-only. Landed experimentally as `_BackgroundAnchor` + config
`pass0_bg_anchor`/`pass0_bg_anchor_prec` (a weak two-sided Gaussian on `log ρ_g` at `ρ_bg`, mode
`log(ρ_bg/ρ_tot)`); default OFF, byte-identical. **To harden to production:** clean up the class/flags, name it
as a first-class calibration stage, and settle the strength (below).

**OI-1 resolved (workflow `wf_a6526c7d` + reviewer doc `rna_parsimony_reviewer_question.md`):** the `½` Jeffreys
reference is CORRECT and unchanged; `−log(1−f_g)` was an improper f_g→1 *promoter*, not a parsimony — a red
herring. The real fix is this restored pass-0 ρ_bg anchor (the shipped v0.7.1 applied it both passes; the
refactor dropped it). Grounding pass-0 also breaks the refit circularity that was contaminating the additive KDE.

**32-scenario benchmark (fixed p=0.5, `scratchpad/bench36.py`, mass-wt region |Δoracle|):** anchor IMPROVES
low/zero-gDNA (0:−0.005, 1:−0.008 — the FP fix), REGRESSES high-gDNA+capture (100:+0.026; big ones
gdna100·present·ss0.99·on +0.165, gdna100·vstrong +0.052, gdna300·on +0.01–0.02, gdna5·vstrong +0.017). Nascent
present regresses more (+0.009) than absent (+0.001). Net +0.006 (slightly negative at fixed strength).

**Refined tuning understanding (supersedes "strength ∝ background"):** the anchor is safe where the gDNA
distribution is **depleted-dominated** and dangerous where **bimodal** (high bg + capture = real enriched gDNA).
The additive KDE refit does NOT fully catch the high-bg+capture crush (regressions are post-refit). So the
adaptivity variable is **enrichment presence / bimodality**, not raw background amount — a stronger anchor at
high bg would regress MORE. Bimodality is visible in the Role-A total-density landscape (read for *caution*,
never composition). **OPEN:** the strength/adaptivity law; and a node-level diagnosis of the big regressions
(esp. the stranded gdna100·ss0.99 +0.165, and the nascent-present skew — possible nascent-modeling gap).

---

## 0. The one-line reframe (from careful study)

The owner's levers **A1 (constrain the background), A2 (additive not competitive), and A3 (drop precision
weighting)** are **not three separable phases — they are one coherent representation**: a **fixed-bandwidth,
occupancy-weighted KDE on the deconvolved gDNA**, plus a **weak separate floor**. A4 (data-driven bandwidth) is
a refinement *on top* of that representation.

**Why A2 alone is insufficient (the subtle finding).** The current NPMLE uses each node's **belief width τ =
√Var(log f_g)** as its kernel width (`_cell_loglik`). The capture-enriched nodes are exactly the low-count,
unstranded, *imprecise* ones (M≈10 ⇒ large τ), so even in an **additive** sum their τ-broadened bumps render
**short** — the enriched mode stays discounted (measured earlier: rendered height 0.19 vs the occupancy ceiling
0.40). Occupancy equals height **only** with a **common (fixed) bandwidth** (A3). So A2 and A3 are codependent,
exactly as A1 and A2 are. All three = "adopt the KDE representation for Role B."

This is the restore the v0.7.1 KDE achieved; we keep the WIP improvements (the derived floor ρ_floor, the
deconvolved-gDNA substrate + 2-pass structure, the σ²_transfer NPMLE for Role A).

---

## 1. The precise derivation (the Role-B density)

**Training set** `T` = REGION nodes that are **single-strand** (`free_pos ^ free_neg`), **live**
(`eff>0, mass>0`), **non-AMBIG**, **EXCLUDING intergenic** (`signature != 0`). Intergenic nodes are (a)
structurally locked to f_g=1 anyway (G1 sinks — they never consume the composition prior) and (b) the
depleted-flood source; they are represented by the floor, not by occupancy. `gonly` non-intergenic structural
seams stay in `T`.

Each node `i∈T` has a **deconvolved-gDNA point estimate** `âᵢ = log(f_g,i · Mᵢ / Eᵢ) = log ρ̂_g,i` (`f_g` from
the pass being refit) and **occupancy 1**.

**The additive KDE density** (collapse to cells `c` of near-equal `âᵢ` for speed; `w_c` = cell node count):

```
π_g(a) = (1/Z) · [  Σ_c  w_c · φ_h(a − â_c)   +   w_floor · φ_{h_f}(a − a_floor)  ]
```

- `φ_h` = Gaussian kernel of bandwidth `h` (decades·ln10). **Fixed / common across nodes (A3)** — NOT per-node
  τ. Initially `h = npmle_bandwidth` (0.15); A4 makes it data-driven (Silverman + median-noise floor).
- `w_c` = cell occupancy; `Σ_c w_c = |T|`. **Occupancy weight (A2), never precision, never mass.**
- `a_floor = background.log_rho_floor` (the derived floor location; `background_reference`).
- **`w_floor = 1` — ONE pseudo-observation (A1)**, not `n_regions`. Weak by construction when `|T| ≫ 1`. This is
  the KDE's `_GLOBAL_STAB_PREC=1` cap, expressed in the additive framework: the whole background counts as one
  pseudo-node. *(This is the one new scalar; it is principled — "one pseudo-observation" — but flag for the
  owner per no-magic-numbers, §4 OI-2.)*
- `h_f` = floor width = `max(sigma_bg, h)` (`sigma_bg = 1/√Σg`; `= h` when `Σg=0` ⇒ `sigma_bg=∞`).
- `Z` normalizes over the grid.

`logP = log π_g` on the grid; the grid spans `[min(a_floor, q0.5(â)) − 3h, max(â) + 2h]` (extend `lo` to include
the floor — already done for the aggregate cell). **`logprior(f_g;M,E)` is UNCHANGED** — the same δ-pin
`np.interp` at `a = log(f_g·M/E)`.

**Contrast with the current fit:** current = per-node Poisson-lognormal `_cell_loglik` (τ-width) → **EM
weights** (`_em_weights`, competitive) + an aggregate background **cell weighted `n_regions`**. New = **fixed-h
kernels at point estimates, occupancy weights, NO EM, weak 1-pseudo-obs floor.**

**Why this keeps the enriched mode** (the three mechanisms, now all satisfied): every node deposits its own
kernel (no EM to compete it away — A2); all kernels share width `h` so occupancy = height (enriched no longer
τ-discounted — A3); the background is one weak pseudo-node, not an `n_regions` tower (A1).

**FP safety (unchanged, verified against the code):** the grid top `= max(â)+2h` is the max *deconvolved gDNA*
rate, so a node's ray `ρ_g=f_g·M/E` above it clamps to the low tail `logP[-1]` (implausibly-high gDNA
disfavored); the additive sum has mass only where deconvolved-gDNA nodes actually are; the floor at `ρ_floor`
is weak. For `gdna_none` (ρ_floor→ the resolution wall, very low): expressed nodes' rays climb *above* the tiny
floor into the clamped low tail ⇒ f_g→0 preserved.

---

## 2. Code touch points (enumerated, exhaustive)

**Change (Role B only):**

| # | file:symbol | change |
|---|---|---|
| T1 | `npmle.py` `DensityNPMLE.fit` | add an **additive KDE build path** (flag `additive=False`). When set: skip `_cell_loglik`/`_em_weights`/aggregate-cell; build `dens_grid` as `Σ_c w_c·φ_h(grid−â_c) + w_floor·φ_{h_f}(grid−a_floor)`; keep the interior-floor + normalization. Populate `weights` with the per-grid occupancy mass (so `project`'s fallback stays well-defined — §4 OI-9), `logP`, `log_rho`, `bandwidth`, `n_cells`. |
| T2 | `calibrate.py:106` `_fit_gdna_hyperprior` | call `DensityNPMLE.fit(..., additive=True)`; **exclude intergenic from `sel`** unconditionally (not only under HYBRID) since they are now the floor; pass `background` for `a_floor`/`sigma_bg`/`w_floor` (NOT as the heavy aggregate cell). |
| T3 | `config.py` `CalibrationConfig` | add a mode flag (`gdna_prior_additive: bool` or a `str` method selector) + surface `w_floor`/bandwidth only if we decide they need config (§4 OI-2, OI-10). Default preserves current behavior until the A/B passes. |

**MUST NOT change (Role A / σ²_transfer, signed-off):**

- `calibrate.py:250` `enrichment_prior = DensityNPMLE.fit(mass_global, eff_global, …)` — **stays the EM NPMLE**
  (`additive=False`). It is fit on **total** density and only ever feeds `project` → σ²_transfer.
- `bp_solver.py:324` `proj_prior.project(...)` — untouched. **Latent coupling (OI-9):** `bp_solver.py:322`
  `proj_prior = enrichment_prior if enrichment_prior is not None else gdna_prior` — if `enrichment_prior` were
  `None`, `project` would run on the additive `gdna_prior`. In production `enrichment_prior` is always set, but
  T1 populates a valid `weights` so the fallback is safe regardless.

**No change needed (verified readers):**

- `bp_solver.py:312` `gdna_prior.logprior(...)` — the additive `logP` is read by the same δ-pin interp.
- `simplex_logodds` `_gdna_arm` (ψ assembly) — consumes the `(m,K)` term as-is.
- `diagnostics.py:44-76` `CalibrationDiagnostics.from_prior` — reads `log_rho`/`logP`/`bandwidth`/`n_cells`,
  all populated by T1; the QC P(ρ) plot will simply now show the KDE density (with the enriched mode).
- `calibrate.py:262,296` — `.n_cells` logging (populated).

**Tests:**

- `tests/calibration/test_npmle.py` — the EM-path + aggregate-cell tests (esp.
  `test_aggregate_cell_zero_counts_anchors_the_derived_floor`) still cover **Role A**; add a **new additive-path
  suite** (occupancy=height; enriched mode survives a depleted-majority; weak-floor placement; `gdna_none`→f_g≈0).
- Golden tests — regen deferred to finalization (already pre-existing red baseline).

---

## 3. Phased execution

**Phase A0 — the additive fit path (isolated, default-OFF). ✅ LANDED 2026-07-19.** T1 done: `_kde_density` +
`DensityNPMLE.fit(additive=False)` in `npmle.py`; `additive=False` byte-identical (14 EM tests + Role A
unchanged); 5 new `test_additive_*` tests pin the properties (232 pass, 1 pre-existing red). Smoke fixture (400
depleted + 100 enriched-imprecise minority): EM starves the minority to height 0.000; additive keeps it at 0.254
= the exact 100/400 occupancy ceiling; a 100k-region weak floor does not crush it. T3 (config flag) deferred to
A1 wiring. No pipeline behavior change yet.

**Phase A1 — wire Role B + render-only validation.** T2 (additive=True, exclude intergenic, weak floor). Fit on
`gdna5_capON` pass-0 and **render-only**: confirm the enriched mode is present at ~occupancy height (the −2.01
mode rises from ~0.19 toward ~0.40) and the depleted-expressed mode is comparable (not a tower). **Gate:** the
fitted `π_g` shows a recoverable enriched mode; `gdna_none`'s π_g collapses to the low floor.

**Phase A2 — gated full-solve A/B (the real test).** Flag ON; full `calibrate()` on the three gates:
`gdna5_capON` recovery (35 nodes → toward 0.63), `gdna_none` FP (median f_g < ~0.02), stranded no-regression.
Iterate `w_floor` (should stay 1) and confirm the δ-pin projection now finds the enriched mode (the KDE proved a
height-read suffices once the density is right). **Gate:** all three; net-flow non-regression on a broader sweep.

**Phase A3 — resolve the open issues that the A/B exposes** (§4), especially the RNA-parsimony reference
strength (OI-1) — **discuss before changing.**

**Phase A4 — data-driven bandwidth.** Replace fixed `h` on the **gDNA prior only** with Silverman +
median-noise-floor (v0.7.1's rule). Secondary; the representation (A0–A2) dominates. **Final selector waits for
real data.**

**Phase B — real-data validation.** The CI sim under-represents the intergenic flood (216 vs 868 here; real =
orders more). Validate `w_floor`, the exclude-intergenic choice, and the bandwidth selector on LBX0190 / MO_3005.

---

## 4. Open issues raised by the study (scrutiny output)

**OI-1 (SHARP — needs discussion). The RNA-parsimony reference strength (½ vs 1).** The shipped KDE's floor was
FP-safe *because* `_kde_logprior` carried `−log(1−f_g)` (coefficient **1**) — without it, "park gDNA at the
depleted mode, dump the rest into free RNA" costs nothing (585k→2k FP). The current ψ uses the **derived
Beta(½,½) reference `½·log(1−f_g)`** (coefficient **½**) via `_rna_arm`. Restoring the enriched mode + a floor
at ρ_bg re-creates an f_g→1 pull for background-level nodes that the **weaker ½ reference may not hold**. But ½
is the *derived-correct* reference (`reference_prior_derivation.md`), and 1 is the KDE's heuristic. **This is a
genuine reference-strength question (NOT the rejected fitted-`logP_r`)** — do not silently change ½→1; measure
the FP gate first (A2), and if it fails, discuss the coefficient explicitly.

**OI-2 (needs discussion). `w_floor = 1` pseudo-observation.** Principled (the KDE's cap), and weak by
construction, but it *is* a new scalar. Confirm "one pseudo-node" is the intended strength, or derive from the
node population.

**OI-3. The floor's f_g→1 pull for background-level nodes.** A node with total ≈ ρ_bg has its floor mass at
f_g≈1; the floor (weak) + RNA reference (OI-1) + strand + messages decide. For true intergenic this is correct
(and usually structurally locked); for a low-RNA node at background level it must be pulled back. Guarded by the
`gdna_none` gate and the OI-1 balance.

**OI-4. Over-smoothing / mode merge.** Fixed `h` = KDE-style marginal density (convolved with noise), mildly
over-smoothed vs the deconvolved latent — acceptable for *detection* (mode present), but monitor that the
enriched (−2) and expressed-depleted (−3) modes stay resolved (1 decade apart; `h=0.15` resolves them). The A4
bandwidth must respect an `h ≤ Δ/2` ceiling.

**OI-5. Dropping τ ignores enriched-node imprecision.** Fixed `h` treats each point estimate as equally sharp,
so the prior *claims* the enriched mode is well-determined when the nodes are low-count. Mitigated by the prior's
overall weakness (n_eff small) and the FP gate; do not let A4's noise-floor re-introduce per-node width.

**OI-6. Determinism.** The additive path removes the EM ⇒ *more* deterministic (good). Independent of the
known C++ parallel-reduction nondeterminism.

**OI-7. Naming.** `DensityNPMLE(additive=True)` is a KDE, not an NPMLE. Keep the flag for the first
implementation (shared interface, minimal risk); revisit a rename/refactor once it lands.

**OI-8. `_collapse` binning.** For the KDE, bin by the density point estimate `â=log(ĝ/E)` (and drop the τ bin),
or place one kernel per existing `(ĝ,E)` cell at its `â_c`. Reuse `_collapse`; note the τ column is unused in
this path.

**OI-9. `project` fallback coupling** — covered by T1 populating a valid `weights` (see §2).

**OI-10. Config surface / magic numbers.** New: the mode flag (structural, fine), `w_floor` (OI-2), the A4
bandwidth selector (real-data). Keep the surface minimal; discuss any tuned constant.

**OI-11. Stratification — REJECTED (owner directive 2026-07-19).** An earlier design called for stratifying the
gDNA fit by signature class. This is **officially rejected and will not be pursued.** The additive KDE with
excluded-intergenic + weak floor carries the enriched/expressed structure directly; there is no stratification.

---

## 5. Acceptance gates & isolation (restated)

- **Isolation:** the additive path applies to the **gDNA composition prior (Role B) ONLY**; `enrichment_prior`
  (Role A → σ²_transfer, signed-off) stays the EM NPMLE, byte-identical.
- **Gates (sim `ambig_dense_10mb`, per A2):** `gdna5_capON` 35-node recovery → toward 0.63; `gdna_none` FP →
  median f_g < ~0.02; stranded (`ss0.99` where available) + resolvable no-regression; broader net-flow sweep.
- **No new magic numbers without discussion** (OI-1, OI-2 flagged).
- **Develop on the sim; finalize `w_floor` + bandwidth on real data (Phase B).**
