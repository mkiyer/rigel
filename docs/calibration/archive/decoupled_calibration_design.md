# Decoupled calibration — design & implementation guide

**Status:** authoritative design for the next calibration architecture. 2026-06-10. Supersedes the
joint count×strand deconvolution. Detailed enough to implement from directly (a dry run). Companion:
`count_channel_capture_design.md` (the count-channel work, now scoped to the count module);
`strand_clean_robust_deferred.md` (the strand estimator lineage); the joint approach is archived in
`archive/joint_deconvolution.md` (to be written during teardown) for possible future resurrection.

---

## 1. Why decouple — the evidence

The joint per-node posterior `posterior(g) ∝ count_prior(g) × strand_likelihood(g)` multiplies two
**unequal** estimators, and the benchmark proved the multiplication is the problem. Removing the
count overdispersion (which had been crushing the count concentration to ~0) made the flagship leak
*triple* (gdna1000 cap_on ss0.99: 3.74% → 10.95%): a biased count prior, once un-crushed, fights and
drags down an excellent strand result. The old overdispersion "worked" only as an accidental crutch
that nullified the count so strand could win.

The root cause is a **statistical asymmetry**:

- The **strand** estimator is *unbiased* (gDNA is symmetric at ½, RNA is at κ; the split is a clean
  contrast), merely *noisy* at low N or low strand-specificity (SS).
- The **count** estimator is *biased* under hybrid capture (the exon gDNA density is imputed from
  depleted off-target neighbours → ~2× under-estimate).

Mixing a biased estimator into an unbiased one **re-introduces bias**. No precision-titration of the
joint fixes this in general — the only robust move is to **not mix them**: let strand decide wherever
it has signal, and use count *only* where strand structurally has none. Across the regimes:

| regime | strand | count | decision |
|---|---|---|---|
| **unstranded** (κ≈½) | no information (no-op) | the only signal | **count** |
| **strong SS** (the common case; real libraries are usually ≥99%, sometimes ≥90%) | excellent, unbiased | mediocre, capture-biased | **strand**; count only for strand-blind (`AMBIG`) nodes |
| **weak SS** (rare) | partial | partial | the *only* synergy case — deferred (see §8) |

Decoupling lets us make each model excellent independently, eliminates the strand **double-use**
(strand currently cleans the count density *and* enters the likelihood), and removes the
un-titratable joint weighting. The two modules act on **disjoint sets of nodes**, so they cannot
conflict.

---

## 2. The decoupled architecture

Two deconvolution modules with a **structural handoff**, not a numeric blend:

```
substrate  → κ = rna_sense_frac + identifiability   (strand_balance)
           → gDNA/RNA strand Beta-Binomial overdispersions  (gdna_strand)   [strand module input]
           → ρ: gDNA density, strand-free, + local imputation (density_model) [count module input]

  per node, route by TWO booleans (no magic threshold):
    strand_identifiable (GLOBAL: is κ usable? = strand_balance identifiability check)
    strand_observable   (PER-NODE: is the node's sense defined?)

    route = strand_identifiable AND strand_observable
      ├─ TRUE  → STRAND module:  BB posterior over g  → gDNA/RNA split   (strand_deconv)
      └─ FALSE → COUNT  module:  g = clip(ρ·eff_len / mass)  → split      (density_model → calibrate)

           → derive (gdna_density_global, geom length)   (unchanged)
           → assemble_priors: boundary-flux transport + per-locus prior   (unchanged)
```

The gate is two already-computed booleans. **Unstranded library ⇒ `strand_identifiable=False` ⇒
every node routes to count** (strand is a true no-op, skipped entirely). **Stranded library ⇒**
strand handles every strand-observable node (exons, introns, gene-edge boundaries via the
intergenic-wildcard rule); count handles only `AMBIG`/no-defined-sense nodes.

Mass conservation, the per-node (contained/left/right) structure, derive, and the boundary-flux
transport are **unchanged** — only *how each node's gDNA fraction is computed* changes.

---

## 3. The strand module (`strand_deconv.py`)

**Input per strand-observable node:** `sense`, `antisense`, `κ`, `gdna_strand_overdispersion`,
`rna_strand_overdispersion`. **Output:** gDNA fraction `g_strand` + gDNA/RNA mass.

It is the **Beta-Binomial posterior over g** with a **weak prior** — the existing `strand_loglik`
on a grid, with the count prior term **removed**:

```
p(g)        = g·½ + (1−g)·κ                      # node sense rate as a function of gDNA fraction
log_post(g) = log Beta(g; ½, ½)  +  strand_loglik(grid, sense, antisense, κ, od_gdna, od_rna)
g_strand    = posterior median on the [ε, 1−ε] grid       # (variance σ²_strand → the FP-rate quantile)
gDNA = g_strand·M ;  RNA = (1−g_strand)·M + spliced
```

- **Weak prior** = Jeffreys `Beta(½,½)` (the standard uninformative prior; not a tunable). It only
  matters at low N; at the SS where we *enter* the strand module (κ identifiable) the likelihood
  dominates. This is the **robust Bayesian strand deconvolution** — the exact version of the
  `strand_clean_robust_deferred.md` estimator (its MLE is the same linear unmix `(ŝ−κ)/(½−κ)`, but
  the posterior is bounded on `[0,1]` and overdispersion-widened, so it has **no clip bias** — the
  very bias that caused the nascent regression in the Phase-0 experiment).
- **No-op property:** if ever entered at κ≈½ the likelihood is flat → posterior = prior → wide,
  uninformative. But the global gate normally prevents that (unstranded ⇒ routed to count), so the
  no-op is enforced structurally, not relied on numerically.
- **Unbiased even when noisy:** a sparse strand-observable node yields a wide-but-unbiased posterior;
  aggregated across nodes into the per-locus prior, the noise averages out. We therefore keep strand
  here even at low N rather than handing to the *biased* count (see §5, Q2).

This is `joint_deconv._joint_per_node` with the `count_*` terms deleted and the prior switched from
`Beta(count·gf+…)` to `Beta(½,½)`.

---

## 4. The count module (`density_model.py` + `calibrate`)

**Used for:** `AMBIG`/no-defined-sense nodes, and *all* nodes when the library is unstranded.
**Input per node:** node mass `M`, effective length, local gDNA density `ρ_local`.
**Output:** `g_count = clip(ρ_local · eff_len / M, 0, 1)` → gDNA/RNA mass.

- **ρ (global gDNA density) is strand-free** — read from **intergenic** count-observable regions
  (`TS_NONE`: no transcript ⇒ pure gDNA, no nascent), `ρ = Σ gDNA_count(intergenic) / Σ eff_len`.
  No strand-cleaning anywhere in the count module (this is what removes the double-use).
- **Local imputation (capture):** the existing boundary-anchored bidirectional sweep
  (`density_model`), with the **strand-cleaning stripped out** (`strand_clean_gdna_frac` deleted):
  a non-observable node takes the gDNA density of its adjacent observable boundary sides, swept
  inward through runs.
- **No Beta prior, no Jeffreys/own-mean floor** — `g_count` is a direct density ratio, so the whole
  `_MIN_CONCENTRATION`/`_JEFFREYS` question **dissolves**. The gdna_none false positive is fixed
  *structurally*: a pure-RNA library has `ρ≈0` (intergenic carries ~no gDNA) ⇒ `g_count≈0` ⇒ no
  manufactured gDNA, with no floor logic needed.

**Known limitation (the count module's job to improve):** under capture the imputed `ρ_local`
under-estimates on-target exon density (the ~2× bias). This is now *isolated* to where the count
module is the sole signal — **unstranded + capture** — and is exactly the target of
`count_channel_capture_design.md` (the point-5 unspliced-fraction projection: `g ≈ f·C` using the
node's own enriched count, and anchor-strength precision). Nascent RNA contaminating the count
density is likewise a count-module concern only at `AMBIG` nodes / unstranded libraries (where strand
can't separate it anyway); in a stranded library the strand module removes nascent at every
strand-observable node.

---

## 5. The two questions, answered

**Q1 — graceful no-op when unstranded.** Use the Beta-Binomial we already have as the strand module's
posterior (§3); its variance → prior variance as κ→½, and the **global identifiability gate** routes
unstranded libraries entirely to count, so strand is skipped, not merely flattened. We do **not** port
the precision-weighted *shrinkage* from `strand_clean_robust_deferred.md` — that was the Gaussian
*approximation* to this exact posterior, built for the closed-form pre-clean we are deleting. Keep the
lesson (confidence → 0 at κ=½), use the exact BB.

**Q2 — count must not ruin good strand (and vice versa).** They act on **disjoint node sets**, so a
strand-observable node is never touched by count — no titration to get wrong. The asymmetry settles
the edge cases: *noisy* strand (low N, still stranded) is unbiased and aggregates correctly, so we
keep it over the *biased* count; *bad* strand means *unidentifiable* strand (κ≈½), which is precisely
the gate's handoff to count. We never blend a biased estimate into an unbiased one.

---

## 6. Teardown / cleanup (file-by-file)

The goal is an elegant tree with **no coupled residue**. Current state: the Phase-0/1 experiment is
already applied (overdispersion removed, own-mean floor added) — both are **superseded** here, so the
teardown starts from the committed-plus-Phase-0/1 tree and removes the coupling itself.

- **`joint_deconv.py` → `strand_deconv.py`** (rename; role changed). Strip `_joint_per_node` to the
  strand-only posterior (§3): remove `count_gdna_frac_in`, `count_concentration_in`, the
  `a_c/b_c` count-prior construction, and `_MIN_CONCENTRATION`; the prior becomes `Beta(½,½)`. Keep
  the grid, `strand_loglik`, the median, and the FP-aversion shift (since replaced by the
  `gdna_deconv_quantile` knob — see `phase2_design.md`). `deconv_regions` /
  `deconv_sides` keep their sense/observability machinery but now feed the strand posterior only.
  `_compute_side` loses its count-prior-density computation (the `own`/`density`/`count_gdna_frac`
  block) — it now only produces `sense/antisense/strand_observable` for the strand module.
  `boundary_side_seeds` stays (gdna_strand fit).
- **`density_model.py`** — becomes the count module's density estimator. Delete
  `strand_clean_gdna_frac` and the `_UNSTRANDED_TOL` constant; `node_gdna_density` no longer
  strand-cleans (`clean_count` → raw `n_unspliced`), and `ρ_global` is read from **intergenic** regions
  only (not all count-observable). `count_gdna_frac` becomes the count module's `g_count`
  (= `clip(ρ·eff/M)`), `count_support` is no longer needed as a *concentration* (the count module has
  no Beta) — drop it unless reused.
- **`calibrate.py`** — add the routing: compute `strand_identifiable` (from `strand_balance`) and
  per-node `strand_observable`; route each node to the strand posterior or the count ratio; assemble
  the per-node gDNA/RNA mass arrays from the routed results. Remove the joint wiring.
- **`config.py`** — confirm no coupled knobs remain (`count_overdispersion_*`,
  `strand_clean_prior_weight` already gone). The FP-aversion knob (now `gdna_deconv_quantile`),
  `n_grid`, strand-overdispersion priors stay.
- **`strand_likelihood.py`, `gdna_strand.py`, `strand_balance.py`, `derive.py`, `priors.py`** —
  unchanged (the BB likelihood, the overdispersion fits, κ, derive, and the boundary-flux transport
  are all reused as-is).
- **Tests** — `test_joint_deconv.py` → `test_strand_deconv.py` (assert the strand posterior in
  isolation; drop the count-prior-flat fixture). Add a `test_routing` (unstranded → all count;
  stranded → observable nodes strand, AMBIG count). `test_density_model.py` — drop strand-clean
  assertions, assert raw intergenic ρ + imputation. Re-run/regenerate golden after the benchmark
  confirms the change.
- **Docs** — fold the decoupled flow into `calibration_theory.md` §4 (replace the joint §4.4/§4.5
  with the strand-module/count-module/handoff); update `ARCHITECTURE.md`, `CLAUDE.md`,
  `CALIBRATION_TODO.md`. Magic-number audit (§9) recorded in the theory doc.

---

## 7. Archiving the joint approach

The joint count×strand posterior is a legitimate design we may resurrect (if §8's synergy proves
real). Preserve it without leaving dead code in the tree:

1. Write `docs/calibration/archive/joint_deconvolution.md` — the current joint architecture in enough
   detail to rebuild it: the per-node `posterior(g) ∝ Beta(count_gdna_frac, count_support) ×
   BB_strand`, the count-prior mean/concentration, the overdispersion-limiting that down-weighted the
   count, *why* it was retired (the bias asymmetry of §1), and the exact files/functions it lived in
   (`joint_deconv._joint_per_node`, the count-prior terms). Link the git SHA of the last joint commit.
2. The code itself is removed (recoverable from git history via the recorded SHA) — the archive doc is
   the human-readable resurrection guide. No `if coupled:` branches left in the live path.

---

## 8. Future re-coupling (final, likely-unneeded step)

After both modules are individually optimized, *consider* synergy in the two regimes where neither
channel alone is strong: **low-count nodes** and **weak-SS libraries** (both rare). Approaches, from
simple to current:

- **Convex weight (simplest):** `g = w·g_count + (1−w)·g_strand`, `w` a coupling weight. Easy, but
  picking `w` is the unsolved problem — it must depend on each channel's reliability.
- **The retired joint (more sophisticated):** `posterior ∝ count_prior × strand_likelihood` already
  *is* a principled precision weighting — the issue was never the form, it was that the count prior's
  precision was dishonest (overstated) and its mean biased. If the count module (§4 + point-5 work)
  becomes *honest* (accurate mean, calibrated precision), the joint would weight correctly on its own.
- **Recommended order:** decouple → optimize each module → only then revisit coupling, by giving the
  *honest* count estimate a precision and combining with the strand posterior's precision
  (inverse-variance), gated to fire only at low strand information. Treat as optional; do not build
  speculatively.

---

## 9. Magic-number audit

- **Routing gate:** two booleans, both already computed (κ-identifiability, node strand-observability).
  No threshold introduced.
- **Strand module prior:** `Beta(½,½)` — the standard Jeffreys uninformative prior, not a tunable.
- **Count module:** `g = clip(ρ·eff/M)` — a pure ratio, no prior, no floor; the `_JEFFREYS` /
  `_MIN_CONCENTRATION` / `count_overdispersion_*` / `strand_clean_prior_weight` constants are all
  **deleted**.
- Remaining knobs are pre-existing and principled: `n_grid` (grid resolution), `gdna_deconv_quantile`
  (the FP-rate quantile, default ½ — the uncertainty-aware successor to the old `gdna_strand_llr_bias`
  tilt), the strand-overdispersion priors. Net: the decoupling **removes** magic numbers rather than
  adding them.

---

## 10. Validation plan

Run the 20-scenario benchmark (skill `calibration-benchmark`) + the nascent scenario tests
(`tests/scenarios/test_nrna_double_counting.py`, which the `nrna_none` benchmark cannot see).
Acceptance:

- **stranded + capture recovers** to the old strand-deference level or better (gdna1000 cap_on ss0.99
  ≤ ~3.7%) — pure strand, no biased count fighting it.
- **gdna_none false positives stay ~0** (count module ρ≈0 ⇒ g_count≈0; structural, no floor).
- **nascent scenarios** (`g0 nrna_30 ss0.90`) **recover** to passing — the BB posterior has no clip
  bias (the Phase-0 regression cause).
- **unstranded + capture** becomes a *pure count-module* number — expected to remain the hard case
  until the point-5 imputation work lands; track it as the count module's optimization target, not a
  decoupling regression.
- capture-off conditions unchanged.

Implement behind the benchmark gate; regenerate golden only after inspecting the diff and confirming
the wins above.
