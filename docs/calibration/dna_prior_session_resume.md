# DNA Prior (Role B) — Session Resume Guide

Hand-off for resuming the gDNA-hyperprior work after a context clear. Read §0 → §1 → §7, then dive where needed.

---

## 0. Starter prompt (paste at the top of the new session)

> We are building Rigel's **gDNA hyperprior (Role B)** — the DNA-density prior that gently anchors the
> calibration solver at ambiguous nodes (it is the **third line of defense**, after strand and message-
> propagation; it does NOT replace the solver). Read `docs/calibration/dna_prior_session_resume.md` (this file),
> then `gdna_hyperprior_plan.md` (the projection = a **sampling-likelihood**, "what DNA level was this observation
> drawn from?") and `enriched_mode_sensitivity_hypotheses.md`. The work splits into **Goal 1** — fit the DNA-
> density **landscape** to match the oracle across scenarios — and **Goal 2 (the real endpoint)** — the
> **projection** of a node's observed density onto that landscape → an anchor `μ*`. We have a **validated
> projection eval** and an **elegant current best**: the *unified* landscape (`scripts/debug/gdna_landscape_recipe.py`,
> 1 constant, 0 special cases) + an **asymmetric-upward projection** (`scripts/debug/gdna_projection.py`) that
> recovers under-called enriched nodes (`enr_recovery` −0.05 → +0.25). All exploration runs as pure numpy on a
> cached substrate (no calibrate re-runs) via `scripts/debug/gdna_explore_lib.py`. **Constraints (non-negotiable):
> elegance/simplicity > raw accuracy (we prototype on TOY sims; real data may break assumptions — simpler is more
> robust); NO magic tunables without deriving them from a reference scale; pause & reevaluate after each workflow
> round; workflows are OK for exploration.** The immediate open choices are in §6 below.

---

## 1. TL;DR — where we are

- **Goal 1 (landscape) — solved & elegant.** `scripts/debug/gdna_landscape_recipe.py`: one uniform reliability-
  mass weight `w = S0/(v·g/(g+1) + S0)`, zero-native additive Poisson. **1 constant (`S0=(0.15·ln10)²`), 0
  special cases.** Robust cross-suite: mean EMD-to-oracle **0.267** on BOTH ambig (32 scen) and quick (16 scen).
- **Goal 2 (projection) — validated & the endpoint.** Built `L.project` + `L.enriched_sensitivity`. The
  projection *mechanism* is correct (self-consistency ±0.06 on the oracle landscape). But the unified landscape is
  **enriched-blind** (the reliability weight silences the high-variance unstranded enriched nodes → the enriched
  mode has ~no mass → `enr_recovery ≈ 0`).
- **Sensitivity — recovered, elegantly.** Key insight: pass-0 **UNDER-calls** enriched gDNA (directional), and the
  truth is **bracketed** — observed `log10(g/E)` is a *lower* bound, total `log10(mass/E)` an *upper* bound. The
  elegant fix lives in the **projection**: an **asymmetric-upward** read-out (`scripts/debug/gdna_projection.py`,
  `project_asym`, `hup=0.70, hdn=0.02, cap=1.0`) that trusts landscape mass ABOVE the observation. On the *clean*
  landscape it gives `enr_recovery` **−0.05 → +0.25** (ambig +0.27 / quick +0.23), `enr_abs_err` 0.226 → ~0.20 —
  **matching a 15-constant workflow synthesis at 3 projection constants / 0 landscape changes.**

**Bottom line:** the sensitivity problem is solved with a clean landscape + one asymmetric projection. Remaining:
derive the projection's `hup/cap` (last magic), the deferred **specificity** (fabrication) work, and wiring into
the real solver.

---

## 2. The arc (how we got here — the reframes that matter)

1. **The crush.** Node 1055 (unstranded, capture-enriched exon, true f_g 0.902) was crushed by the Phase-2 refit
   to 0.0001. Proven **not a solver bug** (`gdna_crush_dissection_node1055.md`): swap in a prior with an enriched
   mode and it lands 0.863. Cause = the **δ-pin** projection (evaluate the prior along `ρ_g=f_g·M/E` and take the
   dominant mode) against an **enrichment-blind** prior → the tall depleted mode wins.
2. **Reframe to a gentle anchor + sampling-likelihood projection.** The prior is line-3, a gentle nudge. The
   projection should ask "which DNA level was this observation *drawn from*?" (proximity/sampling likelihood — a
   far mode, however tall, has ~0 responsibility), then anchor the solver toward `μ*`. (`gdna_hyperprior_plan.md`.)
3. **Rejected dead-ends** (do not revive): modelling an RNA rate prior `π_r` (inverts the difficulty); the
   enrichment-conditioning log-shift as a *composition* fix; the Stage-0 `logaddexp` floor. All were "make the
   prior *decide* f_g" — wrong task.
4. **Goal 1 landscape exploration** (workflows). Found: additive **Poisson** (zero-native) ≫ point-estimate KDE;
   **zero-count structural** nodes = the zero-DNA representation; a **confidence filter** removes over-call
   pollution. Then a **simplification round** collapsed it to the unified 1-rule recipe. Cross-suite validated
   (0.267 on both suites — no overfit).
5. **Pivot to Goal 2 (projection).** Owner: "the projection is the actual diagnostic." Built + validated the eval.
   Diagnosed the enriched-blindness (enriched nodes land right at ~+0.5 but carry `w≈0.16`).
6. **Sensitivity workflow → the directional-bracket insight → the elegant asymmetric projection.** (This session's
   end state.)

---

## 3. The elegant current best (the two pieces)

**Landscape** — `scripts/debug/gdna_landscape_recipe.py` (`recipe(s) -> density on L.GRID`):
substrate = one uniform membership (single|gonly nodes with mass, OR the zero-count structural anchor);
weight = `w = S0/(v·(g/(g+1)) + S0)` (reliability-as-mass; `g=0` anchor → `w=1` emerges, high-var ambiguous →
`w→0`); fit = zero-native additive Poisson. **Robust (EMD 0.267 both suites) but enriched-blind.**

**Projection** — `scripts/debug/gdna_projection.py` (`project_asym(P, d_obs, hup=0.70, hdn=0.02, cap_up=1.0)`):
the observation is a *lower* bracket, so trust mass ABOVE it (`hup` wide), refuse mass BELOW (`hdn→0`), bounded to
`cap_up` decades; mean read-out. **Recovers the enriched nodes.** `python scripts/debug/gdna_projection.py`
prints `ambig enr_recovery +0.266 err 0.199 / quick +0.226 err 0.204`.

---

## 4. Infrastructure — how to run, the lib, the metrics

**Everything is pure numpy on a cached substrate — NO calibrate re-runs.** Pass-0 is run once per scenario and
pickled. Env: `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`, always
`OMP_NUM_THREADS=1`.

- **`scripts/debug/gdna_explore_lib.py`** — THE shared lib. `import gdna_explore_lib as L`:
  - `L.load_scenarios('ambig'|'quick'|path)` → list of scenario dicts; each has per-NODE arrays (region AND
    boundary): `g_hat` (deconvolved DNA count = f0·mass), `eff`, `G` (ORACLE true DNA count), `var` (=Var(log f_g)
    pass-0 belief variance), `mass`, `f0`, `ntype` (0 intergenic,1 intron,2 exon,3 boundary), `is_region`, `fp`/`fn`.
    `s['group']=(cap, dna, ss, nrna)`.
  - `L.masks(s)` → region/boundary/single(fp^fn)/gonly(~fp&~fn)/ambig(fp&fn)/struct_zero(zero-count intergenic/
    intron)/base(region,has-mass,non-ambiguous)/haveE/havemass/ntype. `L.vpercentile(s, sel, p)`.
  - Fits (→ density on `L.GRID`, len 260, log10 ρ_g −5..2.5): `L.fit_kde`, `L.fit_poisson`, `L.fit_precwidth`,
    `L.fit_poisln`.
  - `L.oracle_landscape(s)`, `L.emd(fit,oracle)`, `L.enriched_l1(fit,oracle)`.
  - `L.project(P, d_obs, h_proj, readout='mean'|'mode'|'median')` — the *default* symmetric projection.
  - `L.enriched_sensitivity(recipe, scen, h_proj, readout)` and `L.enriched_sensitivity_suites(...)` →
    `enr_recovery` (want POSITIVE), `enr_abs_err` (lower better), `fabrication` (zero-DNA canary, want low).
  - `L.evaluate(recipe, scen)` / `L.evaluate_suites(recipe)` → landscape EMD (`mean_emd`, `worst_max_emd`).
  - **NOTE:** the *asymmetric* projection is NOT in the lib's `enriched_sensitivity_suites` (that uses `L.project`);
    use `scripts/debug/gdna_projection.py::score_enriched(suite, recipe, proj)` for it.
- **`scripts/debug/gdna_cache_build.py --suite <path> --out <pkl>`** — (re)build a substrate cache. Caches live in
  the scratchpad dir: `gdna_substrate_cache.pkl` (ambig), `gdna_substrate_cache_quick.pkl` (quick). ~2–4 min/suite
  (slow the first time — it scans BAMs; run in background). Suites under `~/Downloads/rigel_runs/`:
  `ambig_dense_10mb` (32, has AMBIG + full DNA 0..300), `quick_3to1_5mb` (16, gdna300/none, "rnd" nascent). Larger
  batteries not yet cached: `gdna_benchmark_5mb` (20, incl gdna1000), `cfrna` (REAL data), `complex_benchmark_25mb`.

**Metrics summary:** Goal-1 = `evaluate_suites` (EMD-to-oracle, lower better; unified = 0.267/0.539). Goal-2 =
`score_enriched`/`enriched_sensitivity` (`enr_recovery` want +, currently +0.25 with the asymmetric projection;
`enr_abs_err` ~0.20; `fabrication` the zero-DNA specificity canary, currently −0.4..−0.6, deferred).

---

## 5. Artifacts index

**Docs (`docs/calibration/`):** `gdna_hyperprior_plan.md` (★ authoritative — the sampling-likelihood projection +
gentle anchor), `enriched_mode_sensitivity_hypotheses.md` (★ the hypothesis families + benchmarking plan),
`gdna_crush_dissection_node1055.md` (the crush proof + 5-hop δ-pin code path), `gdna_hyperprior_from_scratch.md`
(the "gravity" principle), `gdna_hyperprior_clean_slate.md` / `npmle_crush_dissection.md` (earlier strategic).

**Scripts (`scripts/debug/`):** `gdna_explore_lib.py` (the lib), `gdna_cache_build.py` (cache),
`gdna_landscape_recipe.py` (the unified landscape), `gdna_projection.py` (the asymmetric projection + scorer),
`gdna_landscape_sweep.py` / `gdna_landscape_explore.py` (Goal-1 sweep engines), `gdna_hyperprior_eval.py` (the
END-TO-END real-refit harness), `gdna_validate_recipe.py` (cross-suite), `crush_dissection.py` /
`proximity_projection_test.py` / `pass0_fit_vs_oracle.py` (diagnostics).

**Figures (`docs/figures/`):** `landscape_synth.png`, `landscape_winner.png`, `bandwidth_enriched.png`,
`pass0_npmle_vs_oracle.png`, plus `<stage>_fit_vs_oracle.png`.

**Workflow scripts** (re-runnable/resumable): under `…/workflows/scripts/gdna-*.js` (landscape-explore,
landscape-simplify, enriched-sensitivity). Their full results are in the task `.output` files / journal.jsonl.

---

## 6. Open next steps (immediate decisions — pick one, per pause-&-reevaluate)

> **⚑ CURRENT PRIORITY (owner's pivot): fix pass-0, not the projection.** We discovered the enriched
> under-call is a **pass-0 solver problem** on unstranded data (single-strand enriched nodes under-call gDNA by
> −0.13…−0.30 dec; stranded nodes recover it, −0.01). The asymmetric projection (§3) is a **Band-Aid** for a
> symptom; recovering gDNA at the source cascades to a better landscape AND sensitivity. **Full details + the
> dissection plan + root-cause hypotheses (H-a…H-d) are in `enrichment_sensitivity_worklog.md` §6.** The
> enrichment-sensitivity track below is **paused at the elegant core** (§3) — resume it after the solver work.

1. **Derive the projection's `hup`/`cap` per-node** — the last magic. `hdn→0` is principled (obs is a floor);
   `hup`/`cap` are the *de-bias magnitude* and should scale with the node's deconvolution/belief uncertainty (an
   uncertain node can be under-called more → wider `hup`). Goal: 0 tunables, self-calibrating.
2. **Specificity round** — the fabrication canary drifted up (zero-DNA μ* −0.83 → −0.4..−0.6); the accepted
   sensitivity-first cost. Fix by **evidence-gating** (boundaries self-solve via spliced/intron; spliced fraction
   distinguishes gDNA- vs RNA-enrichment) — the same signals that build the enriched mode should gate out
   fabrication. (Owner suspects this buys specificity back "for free.")
3. **Wire into the solver (Goal-2 end-to-end)** — replace the δ-pin with: fit the unified landscape (on the solved
   nodes), project each node's observed DNA density via `project_asym` → `(μ*, precision)`, feed as the gentle
   **anchor** in ψ (a Gaussian pull toward `μ*`), read out, both solvers. A/B via `gdna_hyperprior_eval.py`. See
   `gdna_hyperprior_plan.md` §3 for the intended wiring (note it predates the asymmetric-projection result).
4. **Bigger battery** — cache `gdna_benchmark_5mb` / `cfrna` / `complex_25mb` and confirm the recipe+projection
   generalize to real data before committing.

---

## 7. Owner principles / constraints (MUST preserve)

- **Elegance/simplicity > raw correctness.** We prototype on TOY sims; real data will violate assumptions, so
  fewer components / fewer derived constants / no special cases win even at a small accuracy cost. (Two workflow
  synthesis winners this session were rejected/distilled for being 15-constant "monstrosities.")
- **No magic numbers/tunables** without deriving them from a reference scale (e.g. `S0` = the bandwidth; the
  variance ceiling = the count-zero-info max-entropy variance). "Adaptive" machinery that reduces to a fixed value
  is over-engineering — prove it earns its complexity.
- **The projection is the ultimate endpoint**, not landscape EMD. Optimize `enr_recovery`/projection, validate on
  the oracle landscape first.
- **Pause & reevaluate after each workflow round.** Workflows are OK for exploration; be creative.
- **Boundaries are powerful** (short but numerous; self-solve via spliced RNA + intron DNA; partially enriched) —
  the key enriched-mode detectors and the evidence path to specificity. The prior is for **all** node types, not
  just exons.
- **Sensitivity vs specificity:** sensitivity now (accept fabrication, watch the canary); specificity later via
  evidence-gating. Always test with **ample single-strand nodes** and across the **full scenario spectrum**.
- Develop on the toy suites; the oracle (`_oracle_per_node`) is the ground truth. Owner drives commits (WIP branch
  `calib-ambig-init-wip`; do not push to main).

---

## 8. Production-code state (⚠ abandoned experimental gates to clean up before wiring)

The δ-pin era left **default-OFF, byte-identical** experimental gates in production (from approaches now REJECTED —
do not build on them; remove them when wiring the new projection):
- `src/rigel/config.py`: `gdna_prior_enrichment_condition` (default False).
- `src/rigel/calibration/calibrate.py`: `_fit_gdna_hyperprior` enrichment-condition eff-scaling + a
  `RIGEL_ADD_GONLY` env gate; `import os`.
- `src/rigel/calibration/bp_solver.py`: `node_sweep` `enrichment_condition` param + logprior eff-scaling.
- `src/rigel/calibration/simplex_logodds.py`: `_gdna_arm` `_STAGE0_FLOOR` (env `RIGEL_STAGE0_FLOOR`, default off)
  logaddexp floor; `import os`.
Tests pass with defaults (byte-identical). Nothing committed this session. The NEW approach (unified landscape +
`project_asym` anchor) **replaces** the δ-pin entirely, so these gates should be reverted/removed during wiring.
