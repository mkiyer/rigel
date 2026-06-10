# Phase 0 + Phase 1 — implementation plan (dry run)

**Status:** implementation plan / dry run. 2026-06-10. Concrete, file-by-file change set for two
count-channel fixes, written to verify the result is clean before editing code. Builds on
`count_channel_capture_design.md`.

- **Phase 0 — tear down the count overdispersion model.** It mis-specifies capture (between-region
  mean heterogeneity booked as dispersion → α explodes → `N_eff` crushed → count channel
  annihilated). Revert to the honest Poisson precision: the count-prior concentration **is** the
  gDNA count behind the node (`count_support`), with no overdispersion limiting. Keep the local
  boundary-density imputation and the boundary-flux transport exactly as they are (audited accurate
  post-transport on stranded+capture: 0.926 vs 0.906 truth).
- **Phase 1 — floor the count prior on the node's OWN mean, not Jeffreys ½.** At zero count
  concentration the prior currently collapses to Beta(½,½) → an unstranded node drifts to a 50/50
  split and manufactures false gDNA. Floor instead at one pseudo-observation centred on the node's
  own `count_gdna_frac`, so a zero-evidence node defers to the count clue's belief.

This change set is **net subtractive**: it deletes a module, a test file, four config fields, and a
constant; it adds one constant and a handful of lines. No new abstractions.

---

## Baseline note — revert the inert strand-clean trial first

The working tree carries an uncommitted, **behaviorally inert** strand-clean shrinkage trial
(`strand_clean_prior_weight` / `tau0`, CALIBRATION_TODO Issue #1) across `density_model.py`,
`joint_deconv.py`, `calibrate.py`, `config.py`. Strand-cleaning is deferred, the trial changes no
benchmark numbers, and it is tangled into the same lines Phase 0 edits. Revert it to the committed
baseline so we start clean:

- `git checkout -- src/rigel/calibration/density_model.py` — its only uncommitted change *is* the
  trial; this restores the committed closed-form `strand_clean_gdna_frac(sense, total, rna_sense_frac)`
  (with `_UNSTRANDED_TOL`) and needs no Phase 0 edits (the imputation stays as-is).
- For `joint_deconv.py`, `calibrate.py`, `config.py` the trial revert is folded into the Phase 0
  edits below (those files change anyway).

(If you'd rather keep the trial parked, say so — but it is dead weight under the deferral.)

---

## Phase 0 — file-by-file

### `src/rigel/calibration/count_dispersion.py` → **delete the module**
`effective_count`, `fit_gdna_count_overdispersion`, `_nb_mom`, `_shrink`, `CountDispersionModel`.
Confirmed only consumed by `calibrate.py` and `joint_deconv.py`, both edited below — no other users.

### `src/rigel/calibration/calibrate.py`
**Imports** — drop the dead ones:
```python
# remove:
from .count_dispersion import effective_count, fit_gdna_count_overdispersion
from .signature import TS_AMBIG          # only used by the deleted count block
# and trim joint_deconv import (boundary_side_seeds no longer used here — gdna_strand keeps it):
from .joint_deconv import boundary_side_seeds, deconv_regions, deconv_sides
#  →
from .joint_deconv import deconv_regions, deconv_sides
```

**`node_gdna_density` call** — drop the trial kwarg:
```python
    node_density = node_gdna_density(
        substrate, region_arrays, region_eff_len, fl_mean, rna_sense_frac=rna_sense_frac,
    )   # (no strand_clean_prior_weight)
```

**The entire count-overdispersion block** (current lines ~141–184: `ts`, `contained_seed`,
`contained_count`, `contained_len`, the `boundary_side_seeds` call, `crossing_count`, `crossing_len`,
`fit_gdna_count_overdispersion`, `region_alpha`, `support`, `count_evidence = effective_count(...)`,
the `logger.debug(... count overdispersion ...)`) → **replace with**:
```python
    # Count-prior concentration = the honest gDNA count behind each node's density estimate
    # (count_support): the contained gDNA count for a count-observable region, the anchoring crossing
    # count for an imputed region (so the concentration already tracks anchor strength), 0 for a
    # no-anchor region (→ the joint deconv's own-mean floor governs). No overdispersion limiting:
    # the count is trusted at its Poisson precision. Capture heterogeneity is handled by the local
    # imputation + boundary-flux transport, NOT by shrinking the count away (the single global α
    # mistook on/off-target mean variance for dispersion and annihilated the count channel — see
    # count_channel_capture_design.md).
    count_evidence = np.asarray(node_density.count_support, dtype=np.float64)
```

**`deconv_sides` call** — drop `alpha_crossing` and the trial kwarg:
```python
    left, right = deconv_sides(
        substrate, region_arrays, node_density, boundary_eff_len,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        gdna_strand_llr_bias=config.gdna_strand_llr_bias,
        n_grid=config.n_grid,
    )
```

### `src/rigel/calibration/joint_deconv.py`
**Import** — remove `from .count_dispersion import effective_count`.

**`deconv_sides`** — remove the `alpha_crossing=0.0` and `strand_clean_prior_weight=0.0` params; in the
inner `_deconv`, the concentration is the raw crossing count:
```python
        count_concentration = sq.n_side   # raw crossing count (Poisson precision; no overdispersion limit)
```
and the `_compute_side` call drops `strand_clean_prior_weight=...`.

**`_compute_side`** — revert the trial: remove the `*, strand_clean_prior_weight=0.0` param and call
`strand_clean_gdna_frac(sense, n_side, rna_sense_frac)` (no `tau0`).

Docstring touch-ups in this module that reference "overdispersion-limited `N/(1+α·N)`" /
`count_dispersion` (module + `deconv_regions`/`deconv_sides` docstrings) → reword to "the count
behind the node (`count_support` / crossing count)".

### `src/rigel/config.py`
Remove four fields + their validation: `strand_clean_prior_weight` (trial), `count_overdispersion_prior`,
`count_overdispersion_prior_weight`, and the three corresponding `__post_init__` checks
(lines ~339–354). Nothing else references them.

### `tests/calibration/test_count_dispersion.py` → **delete**
All assertions target the deleted module. Scan `test_calibrate*.py` for any `config.count_overdispersion_*`
or `count_disp` references and remove (none expected outside this file).

---

## Phase 1 — `src/rigel/calibration/joint_deconv.py` (`_joint_per_node`)

Replace the module constant:
```python
_JEFFREYS = 0.5  # principled prior floor (Beta(½,½) at zero count evidence), not a tunable
#  →
_MIN_CONCENTRATION = 1.0  # one pseudo-observation; floors the count-prior concentration so a
#                           zero-evidence node keeps its OWN count mean, not Jeffreys ½.
```

Replace the prior construction (current lines 96–99):
```python
        count_gdna_frac = min(max(float(count_gdna_frac_in[i]), 0.0), 1.0)
        count_concentration = max(float(count_concentration_in[i]), 0.0)
        a_c = count_concentration * count_gdna_frac + _JEFFREYS
        b_c = count_concentration * (1.0 - count_gdna_frac) + _JEFFREYS
#  →
        count_gdna_frac = min(max(float(count_gdna_frac_in[i]), 0.0), 1.0)
        # Floor the concentration at one pseudo-observation centred on the node's own count mean.
        # Where the count carries no concentration (e.g. an exon with ~no gDNA expected) the prior
        # still asserts that mean, so an unstranded node defers to the count clue's belief instead of
        # drifting to a 50/50 split — the Beta(½,½) drift manufactured false gDNA (Phase 1).
        count_concentration = max(float(count_concentration_in[i]), _MIN_CONCENTRATION)
        a_c = count_concentration * count_gdna_frac
        b_c = count_concentration * (1.0 - count_gdna_frac)
```

**Why this is well-defined on the grid.** `a_c`,`b_c` ≥ 0 and at least one is ≥ 1 (when `gf=0`,
`b_c=conc≥1`; when `gf=1`, `a_c=conc≥1`). The grid is bounded `[1e-4, 1−1e-4]`, so
`(a_c−1)·log_grid + (b_c−1)·log1p(−grid)` is always finite — no improper-prior blowup. Behaviour:
`gf→0 ⇒` posterior mass → 0 gDNA (correct for a no-gDNA exon); `gf→1 ⇒` → all gDNA; intermediate `gf`
with `conc=1 ⇒` a weak prior centred at `gf` that the strand likelihood refines (or, when flat,
governs). A high-`conc` node is unchanged (the `+0.5` was always negligible there).

Update the `_joint_per_node` / module docstrings that say "Jeffreys-floored ⇒ Beta(½,½) ⇒ defers to
strand" to describe the own-mean floor.

---

## Docs to update (no stale overdispersion references left)

- `CLAUDE.md` — the calibration package map lists `count_dispersion.py` and the acyclic-pipeline
  bullet enumerates "count overdispersion"; remove both.
- `docs/calibration/calibration_theory.md` §4 pipeline list + §4.4 (count overdispersion) — delete §4.4,
  renumber; the concentration is now "the gDNA count behind the node."
- `docs/ARCHITECTURE.md`, `docs/README.md` — drop `count_dispersion` from the module tables/index.
- `docs/calibration/CALIBRATION_TODO.md` — close Issue #2 (count overdispersion) as reverted; note the
  count-channel direction now lives in `count_channel_capture_design.md`.
- Mark `docs/calibration/count_overdispersion_integration_plan.md` superseded/obsolete (it documents the
  reverted design) — recommend deletion to avoid a stale design doc lingering.

---

## Validation gate

1. Build (`pip install --no-build-isolation -e .` not needed — pure Python) + `ruff check src/ tests/`.
2. `pytest tests/calibration -v` and full `pytest tests/` (expect the deleted-module test gone; golden
   may shift — regenerate with `--update-golden` only after inspecting the diff).
3. Skill `calibration-benchmark` on all 20 scenarios. Acceptance:
   - **unstranded+capture leaks drop** (gdna100/400/1000 ss0.50: 13.6/17.4/18.7%) — the count channel
     is no longer crushed where it's the only signal.
   - **gdna_none unstranded+capture FP → ~0** (Phase 1: the prior no longer drifts to ½ on no-gDNA exons).
   - **stranded + capture-off conditions unchanged** (no regression; the overdispersion was ~inert there
     and the floor `+0.5` was negligible at high concentration).
4. Re-run `scripts/debug/diag_imputation_truth.py` on gdna1000 ss0.50 cap_on: the post-transport
   exon gDNA fraction should rise from 0.686 toward ~0.907 (the count channel restored).

## Elegance checklist (what the tree looks like after)

- One fewer module (`count_dispersion.py`), one fewer test file, four fewer config knobs, one constant
  swapped (`_JEFFREYS` → `_MIN_CONCENTRATION`).
- `calibrate.py` loses ~35 lines (the whole fit block) for 2 lines; no `α`, no contained/crossing
  split, no `boundary_side_seeds` round-trip for the count fit.
- The count-prior concentration has **one** definition — `count_support` — used identically for regions
  and (as the crossing count) for boundary sides. No overdispersion vocabulary remains in the count path.
- `boundary_side_seeds` survives (gdna_strand.py), `effective_count`/`fit_*` do not.
- No vestigial flags, no dead branches, no "legacy/trial" comments.

## Out of scope (named, not silently dropped)

The unstranded+capture residual will shrink but not vanish — closing it needs honest on-target density
(the point-5 unspliced-fraction projection and anchor-strength imputation variance), tracked in
`count_channel_capture_design.md` §3.2–3.3. Phase 0/1 deliberately do **not** attempt it; they remove
the catastrophe and the false positive, and restore the count channel as the foundation those build on.
