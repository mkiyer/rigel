# Mixture-prior redesign — meticulous implementation plan

**Goal (user, 2026-07-05):** eliminate the calibration cross-run **nondeterminism** by deleting the bistable
`MonotoneVarMean` σ²_g P-spline, replacing the conflated global prior with the clean **mixture** architecture
(`gdna_prior_clean_slate_architecture.md`). **Bar:** come *close to current 16-scenario performance with no
substantial regression* (small regressions OK). We are NOT required to fix the KDE's pre-existing accuracy
bugs (the 86k boundary leak — `KDE_boundary_prior_review.md`; the capon-unstranded exon peak-smoothing) — those
are in the baseline too; we must just **not worsen them** while removing the spline.

## 0. The target architecture (what we're building)

Today (per `PHASE2_gdna_mixture_fit_design.md §8f`): pass-2 prior = `_global_logprior(...) + _kde_logprior(...)`
where `_global_logprior` is CONFLATED — it carries (a) the depleted-floor override for intergenic+intron
(`ρ_floor` + `s2_floor`, both DETERMINISTIC), (b) the genome-wide baseline for exons/boundaries
(`ρ_global` + the BISTABLE σ²_g spline), and (c) a `_GLOBAL_STAB_PREC` cap. The σ²_g spline is the *only*
nondeterminism source; its inflated-under-capture value is an artifact crudely doing the KDE's
enriched-avoidance job.

Target: **one mixture prior on `log ρ_g`** = a DETERMINISTIC background floor (the `_floor_estimate` term, kept
verbatim) + the **seeded KDE** (its always-present background mode via the existing `floor_log_rho`/`floor_weight`
seed, plus its learned enriched mode). Delete `ρ_global` + the σ²_g spline (`_gdna_seed_estimate`'s spline output,
`variance_model.py`). The KDE's enriched mode — not a wide σ²_g — is what spares enriched nodes; the seeded
background mode + `_floor_estimate` term — not `ρ_global` — is the downward anchor for gDNA-poor AMBIG nodes.

**The crux hypothesis to validate:** seeding the KDE's background mode restores the exon AMBIG downward anchor
that `ρ_global` provided (§8f says the *un-seeded* KDE lost it; a direct test earlier — dropping the genome-wide
baseline with an un-seeded KDE — exploded no-gDNA exon FP 32→110k). **Everything hinges on the seed replacing
that anchor.** If it does, the redesign is deterministic + non-regressive.

## 1. Phase 0 — make the KDE fit itself deterministic (PREREQUISITE)

The KDE must be ε-robust or the redesign inherits new nondeterminism. Audit `gdna_density_prior.py`:
- **Bandwidth `silverman`** (`_silverman_bandwidth`): weighted std + weighted IQR via `np.interp` — continuous;
  `np.argsort` on distinct values is stable. Deterministic. ✓ (default; keep.)
- **`h_floor = _weighted_median(std, w)`** (the bandwidth floor): a `searchsorted` selection — DISCRETE, the one
  suspect. Replace with the **interpolated** weighted median (continuous 0.5-weight crossing; same primitive we
  designed earlier). Small, local change.
- **`lscv`** bandwidth (not the default): a discrete `argmax` over candidates — leave unused; if ever enabled,
  it needs the same treatment. Note it.
- `_find_modes` is diagnostic-only (not in the solve path) — ignore.

**Test P0:** controlled-ε perturbation (`isolate_amplification.py`-style) on a run that uses ONLY the KDE prior
(temporarily force `_global_logprior→0`, KDE on): amplification must be ≤1× and cross-process bit-identical.
This isolates KDE determinism before wiring.

## 2. Phase 1 — seed the KDE background mode

- `calibrate.py`: compute `ρ_floor` (already available from the calibrate flow / `_floor_estimate`) and pass
  `GdnaDensityPrior.fit(train_sub, bandwidth=..., floor_log_rho=log(ρ_floor), floor_weight=w)`. The seed
  mechanism already exists (`fit` §2.4).
- **Choosing `w` (no magic number):** the seed is a virtual depleted observation; weight it by the DEPLETED
  teacher mass the substrate already carries — e.g. `w = Σ weight over background (intergenic+intron) teachers`,
  or a small fixed number of pseudo-observations (1–a few). Derive it from the substrate, don't tune a constant.
  Start with `w = n_background_teachers`-scaled; make it a `CalibrationConfig` field with a principled default.
- **Test P1:** plot/inspect the fitted `P(log ρ_g)` (the existing plotting framework) on gdna300-capon (should
  be bimodal: seeded background mode at `log ρ_floor` + learned enriched mode) and no-gDNA (unimodal at the
  seeded background). Confirm the background mode sits at `log ρ_floor` and is always present.

## 3. Phase 2 — swap the conflated global for [deterministic floor + seeded KDE]

In `bp_solver.node_sweep`, replace the `global_lp` construction:
- **Keep** the deterministic depleted-floor term for `floor_mask` nodes: build it directly from `_floor_estimate`
  (`ρ_floor`, `s2_floor`, `var_mean_floor`) — a plain Gaussian-on-`log f_g` at `target = log(ρ_floor·E/M)`,
  precision `1/(var_mean_floor + s2_floor)` capped at `_GLOBAL_STAB_PREC`. This is exactly the floor half of the
  current `_global_logprior`, extracted and kept (fully deterministic — verified ε-stable).
- **Drop** the genome-wide `ρ_global` + σ²_g baseline for non-floor nodes. Their anchor now comes from the
  seeded KDE (background mode) in pass-2.
- **Pass 2:** `global_lp = floor_term + _kde_logprior(seeded KDE)`.
- **Pass 1:** `global_lp = floor_term` only (no KDE yet). This gives gDNA-poor AMBIG nodes the floor's downward
  anchor while the KDE trains on the (strand-informative) single-strand + structural teachers. (The floor term
  is deterministic, so pass-1 is deterministic.)
- **Test P2 (THE crux):** in-process, gdna300-capon (should read gDNA, not crushed), no-gDNA capon exon FP
  (must stay ≈baseline ~32, NOT explode — this is the seed-restores-anchor test), stranded (unchanged). If
  no-gDNA exon FP explodes, the seed weight is too low OR the floor term must also apply to exons — iterate
  (§7 fallbacks).

## 4. Phase 3 — delete the spline machinery

Once P2 passes: delete `_gdna_seed_estimate`'s σ²_g output (return only `ρ_global` if still needed by the floor,
else drop `_gdna_seed_estimate` entirely — check: after Phase 2 `ρ_global` is unused → delete it too),
`_fit_seed_varmean`, the `from .variance_model import MonotoneVarMean` import, `src/rigel/calibration/
variance_model.py` (432 lines), and `tests/calibration/test_variance_model.py`. Simplify `_global_logprior` to
the floor-only term (or inline it). `grep` for `MonotoneVarMean`/`variance_model`/`rho_global`/`sigma2_g` to
confirm no live references.

## 5. Phase 4 — validate (the gate)

1. **Determinism:** `isolate_amplification.py 1e-10` → amplification ≤1× (pass-1 AND pass-2); cross-process
   multithreaded ≥3 runs bit-identical on unstranded/stranded/no-gDNA.
2. **No-regression:** the 16-scenario net-flow + abundance benchmark vs the committed baseline
   (`AB.baseline.*`). Bar = close to baseline, small regressions OK. Watch the known-fragile cells:
   no-gDNA-capon (FP), capoff-unstranded (siphon), gdna300-capon-unstranded (leak), nrna_rnd variants.
3. **Tests:** `pytest tests/calibration/ -q`; then `pytest tests/ --update-golden` (goldens shift — the σ²_g
   term is gone); then `pytest tests/` green. Confirm the accumulator byte-for-byte test is untouched (no C++).
4. **Ruff** clean.

## 6. Files touched (scope)
`bp_solver.py` (extract floor term, drop ρ_global/σ²_g, pass-1/2 wiring), `calibrate.py` (seed the KDE +
possibly a `floor_weight`/pass-1 config), `gdna_density_prior.py` (deterministic bandwidth floor), `config.py`
(a principled seed-weight field), DELETE `variance_model.py` + `test_variance_model.py`. ~4 edited + 2 deleted.

## 7. Known risks + the TRUE safety net (corrected)

**Why there is no trivial "ρ_floor-for-exons" fallback:** the current exon baseline is `ρ_global` (pooled,
*includes* enriched → higher target → little downward pull) with σ²_g≈31 (WEAK precision → barely drags). That
weakness is exactly what stops it dragging enriched exons. Swapping in `ρ_floor` (low) + `s2_floor`≈0 (STRONG)
would pull exons hard toward the depleted floor → **drag enriched = the M2 regression**. So a deterministic
*scalar/floor* baseline for exons regresses; the enriched-non-drag MUST come from either (i) the KDE's enriched
mode shielding them (the redesign), or (ii) reproducing the weak σ²_g≈31 value deterministically.

Ranked risks:
1. **The seed doesn't restore the exon AMBIG anchor** (§8f; no-gDNA exon FP re-explodes 32→~110k as when the
   genome baseline was dropped with an *un-seeded* KDE). *Mitigations, in order:* raise the principled
   `floor_weight`; verify the seeded background mode sits low enough (at `log ρ_floor`) with real Gaussian tails
   (the spring that pulls far no-gDNA nodes down); if still failing, this is not fixable by tuning → go to the
   safety net.
2. **The KDE's enriched mode doesn't shield enriched exons well enough** → capon-unstranded leak worsens beyond
   "small." This is the deferred KDE-quality problem surfacing. If it exceeds the "small regression" bar → safety net.
3. **KDE bandwidth-floor nondeterminism** → Phase 0 interpolated-median fix; if residual, drop `h_floor`.
4. **Boundary leak (86k) worsens** when the genome baseline is dropped (boundaries crush to the depleted KDE
   mode — `KDE_boundary_prior_review.md`). *Mitigation:* keep the deterministic floor term applied to boundary
   nodes; don't attempt the boundary-scale fix (baseline bug, out of scope).
5. **floor_weight sensitivity** → sweep a small principled range; pick by the P2 crux test + benchmark.

**THE SAFETY NET (guaranteed determinism + ~zero regression): the deterministic-λ spline.** If the seeded
mixture cannot hold the exon anchor / shield within the "small regression" bar, fall back to *keeping* the
`MonotoneVarMean` σ²_g fit but making its ONE bistable step — the GCV-λ `argmin` in `_select_lambda` — 
deterministic: replace the discrete `argmin` with a **continuous, GCV/BIC-evidence-weighted average of the fits
over the λ grid** (no `argmin` ⇒ ε-robust; reproduces the σ²_g≈31 value the baseline relies on ⇒ ≈zero
regression by construction). This keeps the ~430-line module (the opposite of the cleanup goal) but is the
guaranteed way to ship determinism now. It is the fallback ONLY — we implement toward the clean seeded mixture.

## 8. Realistic effort
- Phase 0–2 (determinism audit + seed + wiring): ~1 focused session; low conceptual risk (extract the floor
  term, wire the existing seed, fix the bandwidth-floor).
- The crux iteration (does the seeded mixture hold the exon anchor AND shield enriched within the bar?): ~1–2
  sessions; **this is the real uncertainty** — it is the deferred KDE-quality problem, now unavoidable.
- Safety net (deterministic-λ spline): ~half a session if the mixture path stalls; guaranteed to unblock.

**Bottom line:** best case = the clean seeded mixture (deterministic, spline deleted, ≈non-regressive). If the
KDE-quality iteration can't hold the bar, the deterministic-λ spline is the guaranteed determinism ship (keeps
the module). We implement toward the clean mixture and only fall back if the benchmark forces it. Either way the
nondeterminism is fixed.

## 9. Critical-review response + refined direction (2026-07-05, code-verified)

An external review raised five points; verified against the code:

1. **Pass-1 anchoring void / substrate pollution — VALID, empirically confirmed.** M2 §7 already *tested*
   dropping the pass-2 genome baseline: no-gDNA capon exon FP explodes **32 → 110,684**. The genome-wide baseline
   is load-bearing for exon FP suppression; the "seed restores the anchor" crux is betting against a measured
   explosion. → Do NOT drop the genome baseline without proving the seed replaces it in-process first.
2. **Boundary-leak ambush — VALID.** The wide σ²_g≈31 is the cushion that lets a boundary's strand overrule the
   KDE crush-to-depleted; deleting it + a sharp KDE worsens the 86k leak (violates the don't-worsen gate).
   Mitigated by point 5.
3. **floor_weight paradox — VALID.** Scaling the seed weight to the (huge) background teacher mass drowns the
   fragile enriched mode → M2 regression; a small fixed weight is a magic number. → If the genome baseline stays
   (below), the seed is unnecessary and this paradox is avoided entirely.
4. **Flat-tail contradiction — WRONG (verified).** The solve uses `logpdf_kernel` (real Gaussian kernel-sum, real
   quadratic tails), NOT the clamped `logpdf` interpolation (`bp_solver.py:743,758`; measured: clamped → ~585k FP
   vs real → ~2k). The "spring" already exists and is load-bearing. No action (rec #2 already satisfied).
5. **Missing per-node precision convolution — VALID, highest-value refinement.** The KDE is applied as a raw
   population density with a single `h_pop` (`bp_solver.py:758`); the substrate already carries
   `log_rho_std = sqrt(Var(log f_g)+1/(gcount+1))` (`gdna_density_prior.py:48`) but uses it only to floor the
   population bandwidth. Apply it per-node: `h_node = sqrt(h_pop² + σ_node²)` (Direction B, `KDE_boundary_prior_
   review.md`). A low-count boundary/exon → large σ_node → softened KDE → its strand overrules → protects
   boundaries AND respects count-zero-info. Cheap, path-independent.

**Refined direction — DECOUPLE determinism from conflation-dissolution.** The review + M2 §7 converge: fully
dropping the genome baseline (the clean-mixture crux) is the risky part. The determinism fix does not require it.

- **Determinism NOW (safe, ≈zero-regression):** ship **M2** (deterministic log-linear σ²_g(μ), closed-form WLS,
  `M2_loglinear_sigma2g_design.md`) as the genome-wide baseline. Reproduces σ²_g≈31 at ρ_global by construction,
  deletes `variance_model.py` (gets past the spline), keeps the load-bearing exon/boundary anchor (no void, no
  ambush). A *better* safety net than the deterministic-λ spline (simpler; deletes the module).
- **KDE quality (safe, path-independent):** the **precision-convolved KDE** (point 5) + the deterministic
  bandwidth (Phase 0). Improves boundaries/exons on any baseline.
- **Conflation-dissolution (the clean mixture, deferred):** dropping the genome baseline + the seeded background
  mode is the elegant endpoint but carries the confirmed 110k risk + the floor_weight paradox + the boundary
  ambush. It becomes its own focused, benchmark-gated effort — not coupled to the release-blocking determinism fix.

The clean-slate mixture (`gdna_prior_clean_slate_architecture.md`) remains the north star; this refinement just
sequences the safe determinism ship ahead of the risky conflation-dissolution.
