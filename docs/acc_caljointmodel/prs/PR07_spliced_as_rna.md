# PR 7 — Spliced fragments are RNA: remove the `eps_s` splice-artifact model

A focused calibration fix from the zero-gDNA deep dive
([calibration_zero_gdna_diagnosis.md](../calibration_zero_gdna_diagnosis.md)).
Lands before the validation PR.

## 0. Problem

The `eps_s` ("gDNA splice-artifact rate") machinery is both **conceptually
wrong** and **degenerate**, and it is the taproot of the zero-gDNA over-call
(single_exon: t1 expected 339, observed 41; 298 RNA reads → gDNA):

- **Trained from the wrong place.** `update_eps_s` estimates `eps_s` from the
  *contained* spliced counts — but a spliced fragment jumps exons across an
  intron, so by definition it spans **multiple regions** → it is a *boundary*
  observation, never *contained*. Contained spliced ≈ 0 (only unannotated novel
  splicing). So `eps_s` is trained from a structurally-empty category.
- **Degenerate default.** With no contained spliced reads, `(1 + Σπ_g·n_s)/(1 +
  Σn_s) = 1/1 = 1` ⇒ `eps_s = 1`.
- **Applied to the populated place.** The E-step applies `eps_s` to *every*
  view's spliced mass, including the boundary views where the real junction
  reads live. With `eps_s = 1`, **all** annotated junction reads are converted
  to phantom gDNA — which inflates ρ₀, fabricates exposure heterogeneity, and
  ultimately makes the gDNA component a fragment-attractor in the locus EM.

## 1. Decision (locked) & the resulting model

**Remove the `eps_s` model entirely — do not just set it to 0, delete the
infrastructure.** Splice artifacts (gDNA fragments that misalign as spliced) are
handled by the **`alignable` splice-junction blacklist** (index
`splice_blacklist.feather`, applied at scan by `filter_blacklisted_sjs`), which
removes artifact-junction reads *upstream* so they are never counted as spliced.
The calibration does **not** model artifacts. This is intentionally simple; the
`eps_s` infrastructure is recoverable from git history if a non-blacklist
artifact model is ever justified.

New per-region soft-allocation (per view):

```
m_g = π_g · mass_unspliced                       # gDNA: the unspliced split only
m_d = (1 − π_g) · mass_unspliced + mass_spliced   # RNA: unspliced RNA split + ALL spliced
```

**Spliced fragments are deterministic RNA** — in every view, including the
contained view. Contained spliced (the rare unannotated novel-splice fragment)
**stays RNA**, exactly as before, just without the `eps_s` gDNA carve-out.

## 2. Implementation plan

| File | Change |
|---|---|
| `calibration/estep.py` | `estep_view`: drop the `eps_s` parameter; `m_g = m_g_unspl`; `m_d = m_d_unspl + view.mass_spliced`. |
| `calibration/mstep.py` | Remove `update_eps_s` and `_EPS_S_PRIOR`. |
| `calibration/calibrate.py` | Remove `eps_s` from `initial_hyperparameters` (now returns `(exposure_dispersion, rho_d_bb)`), the outer loop, the `shared` dict, the result construction, and the iteration log. Remove `_EPS_S_INIT`. |
| `calibration/result.py` | Remove the `CalibrationResult.eps_s` field + its `(0,1)` validation. |
| `pipeline.py` | Drop `eps_s` from the calibration summary log. |
| `scripts/debug/*` | Drop `eps_s` references. |
| tests | `test_mstep` (delete the `update_eps_s` tests), `test_result_schema` (drop `eps_s` from kwargs + the bad-param sweep), `test_estep` (drop the `eps_s=` kwarg), `test_calibrate`, `test_priors` (the `_calibration` helper). |
| docs | `docs/caljointmodel/03_inference.md` §5.5 (ε_s) and §1 generative model: mark the splice-artifact rate **deferred/removed**, with a pointer to the blacklist as the artifact mechanism. |

## 3. What is preserved

- The **spliced channel** (`mass_spliced`, `n_spliced`, the motif-oriented sense
  counts) — still needed for RNA mass and for the strand model
  (`fit_strand_balance`, which never used `eps_s`).
- The **blacklist** path — unchanged; it remains the (only) artifact mechanism.
- Contained spliced → RNA (unchanged behaviour, minus the carve-out).

## 4. Constants (Q6)

**Net removal.** `_EPS_S_PRIOR` and `_EPS_S_INIT` are deleted; **no new
constants**; no new tunable. `eps_s = 0` is implicit (spliced is fully RNA), and
even that symbol disappears.

## 5. Validation

- `scripts/debug/trace_zero_gdna.py` (and forcing eps_s→0, already tested):
  total gDNA → ~0, ρ₀ → ~0, exposure ω → ~1, t1 recovers to ~339.
- Full scenario suite: the zero/low-gDNA failures (single_exon, non_overlapping,
  overlapping_antisense, …) should clear; watch for any gDNA-bearing regression.
- Golden re-baseline (the calibration output changes) + `ruff`.

## 6. Decisions

**Settled (your call):** remove `eps_s` infrastructure entirely; spliced → RNA;
contained spliced → RNA; blacklist is the artifact mechanism; re-addable later.

**Open (minor):** doc 03 §5.5 disposition — delete the ε_s subsection, or keep it
as a "deferred (use the blacklist)" note? Recommend the latter (preserves the
rationale + the re-add path).

## Rollback

Pure removal. Reverting restores the `eps_s` model. Independent of PR 8 (strand).
