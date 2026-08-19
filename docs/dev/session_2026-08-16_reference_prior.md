# HANDOFF — ψ's composition reference gained a per-slot MEAN, and one blocker stops the default flip

    ⚠ **A DEV DOC. Nothing may cite it.** ⛔ It is a HANDOFF, not the state: the ranked list is
    `ROADMAP.md` §1, the state is `ROADMAP.md` §0, and everything settled has already MOVED —
    `EQUATIONS.md` §9c (the derivation), `DESIGN.md` §6b (the rulings), `TRAPS.md` (three named rules),
    `ROADMAP.md` §4.2 (five refused mechanisms). The live design is `docs/dev/rna_arm_design.md`.

    ⭐ Session of 2026-08-16, continuing 2026-08-15. Branch `fragment-length-gold-standard`.
    ⛔ **NOTHING IS COMMITTED.** The tree is dirty; `git status` is the truth.

## ⭐⭐⭐ START HERE — the one-paragraph state

ψ's composition reference is `Beta(a, b)` on the λ grid, and its MEAN was never chosen: `_JEFFREYS_REF = ½`
asserts the library is half gDNA. The derivation (`EQUATIONS.md` §9c) shows the shipped reference is the
degenerate case of a three-term prior whose third term carries the mean. **That term is now in `src/` and
is INERT** — `CompositionPriors.location` plus `simplex_logodds._location_term`, `None` ⇒ not written ⇒
bit-identical. The panel arm that drives it (`vertex_ceiling.py --arm ref_loc=struct`) measures
**0.381 / 0.659 / 0.363 / 0.800** per stratum with the zero control at **1.000** and `noop` byte-identical
on 16/16. ⛔ **One blocker stops the default being flipped**, and it is diagnosed to a named repair (§4).

## 1. WHAT IS IN `src/` (inert, gated, bit-identical)

* `CompositionPriors.location` — a per-slot **scalar** `(m,)`, not an `(m, K)` grid array. `select()`
  slices it with the arms; `regrid()` passes it through unchanged, because a scalar has no lattice to be
  carried between and `regrid-in-the-right-variable` cannot apply to it.
* `simplex_logodds._location_term(lam, location)` → `−log[(1−m)·f_g + m·(1−f_g)]`, formed by `logaddexp`
  over `_log_fg`/`_log1m_fg`. Added beside the two arms at **both** ψ call sites — the single-strand solve
  and the AMBIG one, inside its float32 cast.
* `tests/calibration/test_reference_location.py` — **36 gates**, each with its own perturbation.
* ⭐ Bit-identity proven on real data: term present with `location = None` against `_location_term`
  removed entirely = **19/19 output arrays identical**; forcing a nonzero term moves **10/19**,
  `max|Δ| = 48.47`.
* ⛔ **Nothing sets it.** Behaviour is unchanged and the suite is green.

## 2. WHAT IS IN THE HARNESS (measured, not shipped)

`scripts/design/vertex_ceiling.py --arm ref_loc={noop,struct,struct_grid,struct_soft,pooled,local}`
computes `m_i` per slot and injects it at ψ's ONE `CompositionPriors` construction site.
`scripts/design/object_composition.py` (new, 25/25 self-test) scores the PRIOR alone, with no solver.

**The panel, final Σ|Δ| in FRAGMENTS per stratum, ratio to a `base` re-recorded in session:**

| stratum | `noop` | ⭐ `struct` | `struct_grid` | `struct_soft` | `pooled` |
|---|---|---|---|---|---|
| stranded × capture OFF ⭐IN SCOPE | **1.000** | **0.381** | 0.384 | 0.609 | 0.321 |
| stranded × capture ON ⭐IN SCOPE | **1.000** | **0.659** | 0.660 | 1.045 | ⛔ 5.323 |
| unstranded × capture OFF ⭐IN SCOPE | **1.000** | **0.363** | 0.366 | 0.580 | 0.297 |
| unstranded × capture ON ⛔DEFERRED | **1.000** | 0.800 | 0.800 | 1.000 | 0.805 |
| ⛔ `g00` ZERO-gDNA CONTROL | **1.000** | 1.000 | 1.000 | 1.000 | 0.001 |

⛔ `pooled`/`local` are REFUSED: **5.3× worse** on stranded × capture-ON, because the on-target gDNA
density is under-read 2.6–3.6× pre-solve. ⭐ `struct` is uniformly non-negative.

## 3. WHAT `struct` ASSERTS — and it is exactly the owner's four classes

`m ≈ 1` wherever `¬mrna_active` — no annotated MATURE transcript is continuous across the position.
Measured, that mask is **exactly**: intergenic REGION (1,312) + gene-edge BOUNDARY (2,620) + intron REGION
(9,805) + one-flank-exonic BOUNDARY (19,610) = **33,347 of 70,176 slots**, no remainder. True `f_g` =
**1.0000** on all four on every condition including capture-ON, asserting 1 costs **0 fragments**, and all
four are **EMPTY at `g00`** — which is why the control reads 1.000.

⚠ **It is NOT "no annotated transcription there".** The intron flank IS transcribed, as nascent. That is
why the nascent level is the load-bearing quantity, and why deconvolving it out of the introns (the
owner's density-model route) is what turns the claim from assumed into measured.

⭐ **The strength needs no constant.** The term is worth `log(1/eps)` nats, so capping at the highest grid
point (`eps = sigma(−L)`) makes it exactly `L` — already `sweep_logodds_window`. Override budget
`L/log(2κ)` = **14.6 fragments** on stranded data. ⛔ On UNSTRANDED data nothing overrides it at any depth
(0.9998 at 10,000 fragments), because `κ = ½` leaves no channel to prove otherwise with.

## 4. ⛔⛔⛔ THE BLOCKER, DIAGNOSED — and the repair is named

With the reference forced on, a hand-built pure-gDNA intron reads **0.7661** against **0.9006** shipped
(`test_sweep.py::test_mature_no_nascent_hallucination_in_introns`, bound 0.85). **The wrong direction.**

| | location OFF | location ON |
|---|---|---|
| `fg_loc` — the LOCAL self-solve | 0.9858 | ⭐ **0.9998** — correct |
| `_tau0_lam` — own λ-channel precision | 0.1358 | ⛔ **0.0000** |
| final `f_g` (truth 1.0) | 0.9006 | ⛔ **0.7661** |

⭐ The local solve is RIGHT. What inverts it is `tau_lam` collapsing to zero, and `tau_lam > eps` is
`has_own_composition_evidence` — *"this slot has no own composition evidence"*. **The strongest possible
prior makes a slot read as evidence-FREE.**

⭐⭐ The strand term vanishing at the vertex is **correct physics**: rank-1 in
`p = ½+(κ−½)(1−f_g)sinθ` ⇒ on λ, `c·a² ∝ f_g²(1−f_g)²`, which at `f_g = 0.9998` is ~4e-8 of its value at
0.9858. ⛔ **And a FLOOR is already refused** — `has_own_composition_evidence`'s docstring records one at
`1/(2L)²` derived, implemented and refuted.

⭐⭐⭐ **THE REPAIR IS AN ASYMMETRY, NOT A RULE.** `tau_lam` already accepts a λ-factor's own curvature —
`tau_fac = density_factor_precision(intron_prior, lam_grid)`. The intron factory is a λ-factor and
contributes its curvature; the location term is also a λ-factor and contributes none. **One more caller of
a function that exists.** No floor, no constant.

⚠ `tau_lam` also feeds the relay's licence and `solvability_audit`'s CERTAIN class, so raising it at
vertex objects moves what counts as "confidently wrong". Measure with `solvability_audit.py` and report.

## 5. THE ORDER — and it is NOT gate-first

⛔ The falsification gate cannot come first: nothing sets `location` in production, so a gate against
inert code is vacuous.

| # | what | gate |
|---|---|---|
| **1** | `CalibrationConfig.structural_reference: bool = False`, a builder beside `_location_term` taking `(statics, logodds_window)`, threaded at `sweep.py`'s ONE construction site | OFF ⇒ **bit-identical**; suite green BOTH ways |
| **2** | the 1b gate, **verified failing**: a pure-gDNA intron at `m → 1` keeps nonzero `tau_lam` and its final `f_g` RISES | fails before, passes after |
| **3** | the repair: `tau_lam += density_factor_precision(location_factor, lam_grid)` | the gate flips; `location = None` stays bit-identical |
| **4** | `tests/calibration` with the flag ON as green as with it OFF | baseline **1,185 passed / 7 xfailed / 0 errors** |
| **5** | the mask gate: `¬mrna_active` ≡ the four classes | raises on a fifth shape |
| **6** | the panel through `calibrate`, not the harness | must reproduce 0.381 / 0.659 / 0.363 / 0.800 / 1.000 |
| **7** | thermometer `quant_accuracy.py --arm base`; flip the default; MOVE `rna_arm_design.md` §51–58 out | no end-to-end regression |

⚠ **Deliberately after:** the intron resolution wall (makes `eps` a measurement, not a bound), the
capture-ON on-target density, and the relay. ⛔ The standing xfail
`test_an_rna_free_interior_exon_reads_its_truth` will **not** go green — it is an EXON, so `mature_here`
is true and `struct` leaves it at ½ by design.

## 6. TRAPS THAT COST THIS SESSION

* `a-mean-hits-the-mass-weighted-centre-by-luck` — hit **four** times before it was named.
* `a-clamp-at-the-closed-end-escapes-the-window` — an `L`-sweep over the interior never sees the clamp.
* `the-deconvolution-is-as-good-as-the-density-it-is-handed` — price the density before believing the residual.
* `shard-an-arm-sweep-by-condition` — the instruments WRITE their caches; three arms on one condition
  truncated a `payload.npz`.
* ⚠ **Twice a gate failed and the code was right** — my threshold was wrong both times, and
  re-specifying against the measured structure is what surfaced the real findings.

## 7. STATE

Suite **0 failed / 3,458 passed / 0 skipped / 10 xfail**. Lint clean. ⛔ Nothing committed.
Panel: the 16-condition ladder, unchanged. ⚠ One `_main` oracle-cache entry was deleted and will rebuild
on demand (`gdna_g00_ss_0.99_nrna_none_capture_on`).
