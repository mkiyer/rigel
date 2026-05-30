# PR 3 — Strand-balance model (D2)

**Parent plan:** [`../00_implementation_plan.md`](../00_implementation_plan.md) §4 D2, §7 (PR 3).
**Type:** Python-only (no index rebuild, no C++).
**Build required:** no.
**Status:** **Implemented on branch `pr3-strand-balance`** (after PR 2.5
landed). All decisions resolved. Because PR 2.5 made the substrate's spliced
channels motif-oriented (`n_spliced_sense`/`n_spliced_antisense`), PR 3 needed
**no fold and no substrate extension** — it fits `ρ_r_bb` directly. See §IV
execution notes; Parts I–III below are the pre-PR-2.5 plan (T1/T2 obsoleted).

PR 3 builds the **RNA strand-balance model**: the per-library RNA sense mean
`κ_rna` and its Beta-Binomial overdispersion `ρ_r_bb`. `κ_rna` comes from the
live `StrandModel`; `ρ_r_bb` is fit from **integer sense/antisense counts on
spliced fragments at splice junctions**. These feed the PR 4 E-step's strand
log-Bayes-factor (gDNA-BB vs RNA-BB). No E-step/deconvolution lands here.

## §0 — Revision summary (2026-05-30 review)

Your review materially changed the design. The key insight and its
consequences:

### §0.1 Spliced fragments carry a self-defined strand (revises §I.3)

A spliced fragment's strand is fixed by the **asymmetric splice motif** (GT/AG),
read as the **XS** (or minimap2 **ts**) BAM tag — an orthogonal, library- and
region-independent strand call. So a spliced fragment has a definite strand
**even inside a strand-ambiguous region or boundary**, and **annotated** spliced
fragments (we restrict to annotated junctions to suppress artifacts) can train
the RNA strand model across AMBIG territory. Refinement of D7: only **unspliced**
fragments in AMBIG regions are strand-unrecoverable (→ density/sweep, D7);
**spliced** fragments are recoverable everywhere via the motif.

### §0.2 Resolved decisions

- **III.1 — `κ_rna` source = live `StrandModel`** (`p_r1_sense`, the global RNA
  sense mean over all annotated spliced uniques — highly accurate). PR 3 does
  **not** re-fit the mean.
- **III.2 — fit `ρ_r_bb`** (overdispersion) from spliced junction counts, using
  boundaries as independent splice-junction observations. May empirically
  reduce to a Binomial (`ρ_r_bb≈0`); confirm on real datasets (PR 7).
- **III.3 — terminology** (see §0.4): `pos`/`neg` for genome strand;
  `sense`/`antisense` for transcript-relative. Drop the `k_plus`/`plus` names.
- **III.5 — fit method**: undecided; needs simulation/real-data study. PR 3
  ships a documented default (MoM) and the method is revisited with data.
- **III.6 — min-observations guard**: needs **both** a total spliced `n` and a
  count of **distinct** junction observations (overdispersion needs ≥2 to be
  defined, and more to be meaningful). Still a heuristic — exact thresholds to
  agree (§III.6).
- **III.8 — wire into `calibrate()`**: yes.

### §0.3 Open architecture fork — spliced strand orientation ★

To fit `ρ_r_bb` we need, per junction, **(n_sense, n_antisense)** = spliced
reads aligning sense vs antisense to that junction's RNA strand. The accumulator
records spliced reads by **read genome strand** and knows the **motif strand**
at deposit. Two ways to get the orientation:

- **(A) Orient at deposit (C++).** The scanner computes
  `sense = (exon_strand == sj_strand)` per spliced fragment and deposits into
  sense/antisense channels. Self-contained per fragment; works in AMBIG without
  any annotation lookup; matches `StrandModel`'s definition exactly. Changes the
  spliced channel semantics (spliced = motif-relative sense/antisense; unspliced
  stays genome strand). *Recommended.*
- **(B) Orient downstream (Python).** Keep depositing spliced by read genome
  strand; in PR 3 orient each boundary's spliced counts by the junction's
  **annotated** RNA strand (needs a boundary↔annotated-SJ strand map; ambiguous
  where annotated junctions exist on both strands at one coordinate).

**This is the one decision I need before coding.** It determines whether the
accumulator's spliced channel means "read strand" or "motif-relative sense."

### §0.4 Consequent scope + sequencing

Your III.4 + III.9 + the fork reach into the **C++ accumulator**, which the
master plan had reserved for PR 6. I propose splitting the work:

- **New PR (insert as PR 2.5 / "accumulator strand & spliced-flux correction", C++):**
  1. **Per-side spliced flux (III.4):** `Boundary.flux[4]` → `flux_left[4]` +
     `flux_right[4]`; deposit flux on the side matching the mass (unspliced
     contiguous ⇒ both sides; spliced intron-skip ⇒ one side each — no false
     intron flux). Payload + `AccumulatorPayload` + substrate updated.
  2. **Spliced strand orientation** per §0.3 (if fork **A**).
  3. **Keep spliced fragments with a clear motif** even where read strand is
     ambiguous (don't drop them at deposit).
  4. **Substrate rework (III.9):** expose **raw `pos`/`neg` genome-strand**
     counts (no pre-oriented `k_plus`); orientation moves downstream (E-step:
     region strand + library mode; strand model: motif). Rebuild required.
- **PR 3 (Python, this doc):** the strand-balance model on the corrected
  substrate — `κ_rna` from `StrandModel`, `ρ_r_bb` from spliced junction counts.

Below, Parts I–III are the PR 3 (Python) plan; they assume the accumulator PR
landed and resolve once the §0.3 fork is decided.

---

---

# Part I — Theory & design concept

## I.1 What the strand channel does, and why we need a fitted RNA model

The E-step (PR 4) discriminates gDNA from RNA in a region's unspliced count by
adding three log-Bayes-factors; the **strand** one compares two hypotheses for
the observed sense count `k_+` among `n_u` unspliced reads:

- **gDNA:** `K_+ ~ BetaBinom(n_u, κ_d = 0.5, ρ_d_bb)` — unstranded, mean 0.5.
- **RNA:** `K_+ ~ BetaBinom(n_u, κ_rna, ρ_r_bb)` — stranded, mean `κ_rna`.

For this to discriminate, we need the RNA strand parameters `κ_rna` (mean) and
`ρ_r_bb` (overdispersion). PR 3 produces them. `ρ_d_bb` (gDNA) is the M-step's
job (PR 5); `κ_d = 0.5` is fixed (biology).

## I.2 Why spliced fragments, integer counts, Beta-Binomial

- **Spliced ⇒ deterministic RNA.** A spliced fragment is RNA (gDNA cannot
  traverse a junction except at the tiny artifact rate `ε_s`). So spliced
  sense/antisense counts are a *clean* readout of the RNA strand distribution,
  with no gDNA/RNA mixing to deconvolve — the ideal data to fit `κ_rna`/`ρ_r_bb`.
- **Integer counts, not fractions (D2).** `9/1` and `90/10` both give fraction
  0.9 but carry very different statistical power. The Beta-Binomial likelihood
  on integer `(k_sense, n)` keeps that power; a fraction throws it away.
- **Beta-Binomial, not Binomial.** Real libraries are over-dispersed around
  their strand mean (protocol/technical variation). The BB's `ρ_r_bb ∈ (0,1)`
  absorbs that. This overdispersion is also the **robustness mechanism** — it
  prevents a few imbalanced regions from distorting the model (the principled
  alternative to outlier-clipping, per the Q2 thread).

## I.3 Which regions/sides contribute — NONE ≠ AMBIG (D7)

The fold from absolute (`strand+`, `strand−`) channels to transcript-relative
(`sense`, `antisense`) is keyed on the region's `ts_class`:

| Strand class | Spliced contribution | Why |
|---|---|---|
| **POS / NEG** | **Yes** — sense = `strand+`/`strand−` (flip for NEG) | Single, defined sense. |
| **NONE** (intergenic) | **No** — no annotated RNA ⇒ ~no spliced reads | Naturally absent from the spliced pool; nothing to contribute. |
| **AMBIG** (both strands) | **No — excluded (D7)** | No valid sense split exists; an arbitrary fold would manufacture phantom strand signal. |

The recovered `_fold_pos_neg_by_transcript_strand` already encodes exactly this:
`eligible = (ts_class == TS_POS) | (ts_class == TS_NEG)`. (Note this differs
from the *unspliced* `k_plus` in the substrate, where NONE is kept as a neutral
arbitrary split for the E-step — see PR 2. The RNA spliced fit uses the
stricter POS/NEG-only eligibility.)

**Observation sources (D2):**
1. **Contained spliced** of each POS/NEG region: `(k_sense, k_antisense)` from
   spliced channels `ch2`/`ch3` oriented by region strand.
2. **Boundary sides** whose adjacent region is POS/NEG: `(k_sense, k_antisense)`
   from the boundary's **integer flux** spliced channels, oriented to that
   side's adjacent region strand. AMBIG-adjacent sides skipped. (The shared
   flux is oriented separately for each side — see §III.4.)

Pool all integer pairs → one Beta-Binomial fit.

## I.4 The fit (generalize the recovered symmetric estimator)

`strand_balance.py @ fc96902` is a **symmetric** (mean fixed 0.5) Beta-Binomial
method-of-moments estimator for gDNA. PR 3 generalizes the fixed mean to the
RNA mean and reports the overdispersion as `ρ_r_bb`:

- Pool `n_pos = Σ k_sense`, `n = Σ (k_sense + k_antisense)` over eligible
  observations.
- **Mean** `κ_rna` (see §III.1 for the source decision).
- **Overdispersion** `ρ_r_bb` by method-of-moments around the mean: compare the
  observed residual variance `Σ (k_sense − κ_rna·nᵢ)²` to the binomial
  expectation `Σ κ_rna(1−κ_rna)·nᵢ` and the max (fully-correlated) variance
  `Σ κ_rna(1−κ_rna)·nᵢ²`; solve for the intra-class correlation `ρ_r_bb`. The
  recovered `_log_beta_binom_pmf` provides an exact BB likelihood for a 1-D
  refinement / unit tests.
- Count-based uncertainty: the estimate carries `n` (total spliced reads) so
  downstream knows its statistical weight.

`ρ_r_bb` relates to the BB concentration `s` by `ρ = 1/(s+1)`; the recovered
code parametrizes by `s` (its `kappa`) — PR 3 reports `ρ_r_bb ∈ (0,1)` to match
the model/`CalibrationResult` schema.

## I.5 Relationship to the live `StrandModel`

The scan already trains a library `StrandModel` (`p_r1_sense` = P(read strand ==
transcript strand), with a Beta posterior and a read-1-sense *mode*) from
spliced unique mappers. `p_r1_sense` **is** an RNA sense-mean estimate, from a
large clean sample. PR 3's per-region spliced counts measure the same quantity
but sparser and post-resolution. Whether `κ_rna` should *be* `p_r1_sense` or be
re-fit from the substrate is the central open decision (§III.1).

---

# Part II — Implementation plan (pending §III)

> Shapes below assume the §III decisions resolve to the **recommended** option;
> they will be finalized after critique.

## II.1 Target surface

```
src/rigel/calibration/
  strand_balance.py   # NEW: StrandBalance dataclass + fit_strand_balance(...)
  substrate.py        # EDIT: add k_plus_spliced to SubstrateView (§III.3)
  result.py           # EDIT?: add kappa_rna field (§III.7)
  calibrate.py        # EDIT?: compute strand balance, store κ_rna/ρ_r_bb (§III.8)
  __init__.py         # re-export StrandBalance
```

Calibration package: 11 modules (was 10). Within the ≤25 budget.

## II.2 — T1: spliced sense in the substrate (if §III.3 = "extend")

Add `k_plus_spliced: int64[R]` to `SubstrateView` (sense among **spliced**,
oriented by `ts_class` with the same `sense_is_pos` rule as `k_plus`). Antisense
spliced = `n_spliced − k_plus_spliced`. Applies to all three views
(contained / left / right). Additive; PR 2 tests unaffected.

## II.3 — T2: observation pooling (`strand_balance.py`)

`fit_strand_balance(substrate, region_arrays, strand_model, config) -> StrandBalance`:

1. **Eligibility mask** `decodable = (ts_class == TS_POS) | (ts_class == TS_NEG)`.
2. **Contained**: for decodable regions with `n_spliced > 0`, collect
   `(k_plus_spliced, n_spliced)`.
3. **Boundary sides**: for decodable regions, collect the left-view and
   right-view `(k_plus_spliced, n_spliced)` from the flux spliced channels
   (skip where the adjacent region is AMBIG — i.e. only decodable regions).
4. Concatenate into integer arrays `(k_sense[], n[])`.

## II.4 — T3: the fit + fallback

- `κ_rna` per §III.1; `ρ_r_bb` by MoM (§I.4), clamped to `(BB_FLOOR, 1−BB_FLOOR)`.
- **Fallbacks** (§III.6): too few spliced observations ⇒ `κ_rna` falls back to
  the `StrandModel` mean (or 0.5 if the library is unstranded/untrained) and
  `ρ_r_bb` to a documented default; the result records `fallback_used`.
- Return `StrandBalance(kappa_rna, rho_r_bb, n_observations, n_sense, n_total,
  fallback_used, fallback_reason)`.

## II.5 — T4: wiring (per §III.8)

If we wire now: `calibrate()` calls `fit_strand_balance(...)` and stores the
real `κ_rna`/`ρ_r_bb` in `CalibrationResult` (replacing the placeholder inits).
The E-step that *uses* them is still PR 4.

## II.6 — T5: tests (`tests/calibration/`)

- `test_strand_balance.py` — fold orientation (POS/NEG flip; NONE/AMBIG
  excluded) on a synthetic substrate; perfectly-stranded data ⇒ `κ_rna→1`,
  small `ρ_r_bb`; over-dispersed data ⇒ larger `ρ_r_bb`; integer power
  (`9/1` vs `90/10` give different uncertainty); empty-pool fallback.
- `test_substrate.py` — extend for `k_plus_spliced` (POS/NEG/NONE/AMBIG).
- `test_calibrate_placeholder.py` / `test_result_schema.py` — updated if the
  schema/wiring decisions add `κ_rna`.
- BB-PMF unit test recovered from `test_strand_deconv.py @ fc96902`.

## II.7 Acceptance gate

- New tests pass; PR 1/PR 2 suites green; `ruff` clean.
- `calibrate()` still returns a schema-valid result; full suite failure mode
  unchanged (post-calibrate `NotImplementedError`).
- Magic-number budget: every constant in §III.6 explicitly agreed.

---

# Part III — Open decisions (please critique)

### III.1 — `κ_rna` source: live `StrandModel` mean, or re-fit from the substrate? ★ (biggest)

`κ_rna` (RNA sense **mean**) can come from:

- **(A) `StrandModel.p_r1_sense`** — the library mean already fit from spliced
  uniques at scan time. Matches **docs 01 §4.3 / 03 §5.3** ("κ_rna pre-computed
  by `StrandModel`; only `ρ_r_bb` fitted"). Higher-quality mean (more data,
  unambiguous uniques); PR 3 then fits only the **overdispersion** `ρ_r_bb`.
- **(B) Re-fit from substrate spliced counts** — a literal reading of **D2**
  ("generalize the fixed-0.5 mean to a fitted RNA mean"), using `StrandModel`
  only for the read-1-sense *direction*. Fully substrate-derived; redundant
  with `p_r1_sense` and sparser.
- **(C) Hybrid** — mean from `StrandModel`, but override per-region if a region
  has strong local spliced evidence (adds complexity; probably v2).

**Recommendation: (A).** It reconciles D2 with docs 01/03, uses the best mean,
and leaves PR 3 fitting the one genuinely new parameter (`ρ_r_bb`). D2's "fitted
mean" is then satisfied by *using a non-0.5 mean in the BB*, not by re-deriving
it. **Your call** — this sets PR 3's whole shape.

### III.2 — Does PR 5 re-fit `ρ_r_bb`, or is PR 3's spliced estimate final?

doc 03 §5.3 fits `ρ_r_bb` in the M-step from **unspliced soft-allocated** RNA
counts; D2 fits it in PR 3 from **spliced** counts. Two different data sources.
**Recommendation:** PR 3's spliced `ρ_r_bb` is **final** (spliced = clean RNA,
no soft allocation needed); PR 5's M-step fits only the gDNA `ρ_d_bb`. Confirm,
so PR 5's scope is settled.

### III.3 — Spliced sense: extend the substrate, or compute in PR 3?

The substrate (PR 2) exposes `n_spliced` (total) but not spliced sense.
- **(A) Extend `SubstrateView`** with `k_plus_spliced` (recommended — the
  substrate is the payload→sufficient-stats bridge; keeps all channel decoding
  in one place; additive).
- **(B) Compute in PR 3** from the raw payload (keeps PR 2 frozen; duplicates
  channel-decoding logic).

**Recommendation: (A).** Confirm.

### III.4 — Boundary shared-flux used once per side (confirm D2 intent)

A boundary's integer `flux` is shared by its two neighbours; D2 attributes it to
**each** side oriented to that side's region strand, so the same flux integer
enters the pool **twice** (once per adjacent region). This is the locked D2
"per side" model, but it means a boundary read counts as a strand observation
for *both* neighbours. Confirm this is intended (vs. counting each boundary
once). Low-stakes for the *mean*; mild for the overdispersion's effective `n`.

### III.5 — Fit method: method-of-moments vs exact MLE

- **(A) MoM** (recovered, closed-form; generalize the symmetric estimator) —
  recommended for PR 3: simple, fast, no optimizer, well-tested.
- **(B) 1-D Newton MLE** on the BB log-likelihood (doc 03 style) — more precise
  at small `n`, adds an optimizer + bounds.

**Recommendation: (A)** for PR 3 (MoM gives `ρ_r_bb` with count weight); revisit
in PR 5 only if benchmarks show MoM bias. Confirm.

### III.6 — Magic numbers / fallbacks (Q6 — need explicit sign-off)

The recovered code carries several constants; PR 3 needs a small, agreed set:

| Constant | Recovered value | Proposed for PR 3 | Status |
|---|---|---|---|
| `ρ_r_bb` floor/ceil | `kappa∈[1e-3,1e6]` | `BB_FLOOR = 1e-6` (doc 03 §8) | derived; OK? |
| `κ_rna` fallback (no data) | `0.5` | `StrandModel.p_r1_sense`, else `0.5` (doc 03 §8 `_K_RNA_FALLBACK`) | OK? |
| `ρ_r_bb` default (no data) | `kappa=1e6` ⇒ `ρ→0` | a documented default (e.g. the M-step init `0.01`)? | **discuss** |
| Min observations to fit | `MIN_REGIONS_FOR_STRAND_MOM = 2` | a min total spliced `n`? | **discuss — heuristic** |

I will **not** introduce any of these without your agreement. In particular the
"min observations" guard is a heuristic; I'd prefer to express robustness via
the count-weighted likelihood (a thin pool → wide uncertainty → naturally
defers to the fallback) rather than a hard region-count cliff. Your guidance?

### III.7 — Does `CalibrationResult` gain a `κ_rna` field?

doc 04 §5 currently excludes `κ_rna` ("pre-computed by `StrandModel`, not in
schema"). If PR 3 fixes/fits it, recording it aids provenance and PR 4. Options:
**(A)** add `kappa_rna: float` to `CalibrationResult`; **(B)** keep it in the
`StrandBalance` object passed PR3→PR4 only. **Recommendation: (A)** (one scalar,
useful in the §10 log line). Confirm.

### III.8 — Wire into `calibrate()` now, or keep standalone until PR 4?

**(A)** PR 3 wires `fit_strand_balance` into `calibrate()` and stores real
`κ_rna`/`ρ_r_bb` in the result (the placeholder still infers no gDNA, but now
reports the true strand model). **(B)** Keep `fit_strand_balance` standalone +
tested; PR 4 wires it. **Recommendation: (A)** — exercises the real inputs end
to end and de-risks PR 4. Confirm.

### III.9 — Channel strand semantics / library mode (confirm, not a fork)

The fold assumes the accumulator's `strand+`/`strand−` channel is the fragment's
**transcript-relative** strand (so for a POS region, `strand+` = sense). I need
to confirm the scanner's convention and whether the `StrandModel` read-1-sense
*mode* must flip the mapping for reverse-stranded protocols. If the scanner
already emits transcript-relative strand, the fold is purely `ts_class`-based
(as recovered) and no flip is needed. **I'll verify against `bam_scanner.cpp` /
`StrandModel` before coding** and report; flagging so the answer is explicit.

---

## Rollback

Revert the PR. `calibrate()` reverts to placeholder `ρ_r_bb`; `kappa_rna` is
dropped from `CalibrationResult`. No on-disk artifacts change.

## IV. Execution notes (implemented 2026-05-31)

PR 2.5 simplified PR 3 substantially: the substrate already exposes
**motif-oriented** `n_spliced_sense` / `n_spliced_antisense` per view, so the
recovered `_fold_pos_neg_by_transcript_strand` + the `SubstrateView`
extension (old §I.4, §II.2) were **not needed**. Spliced observations are
valid in AMBIG regions (motif-defined), so there is **no `ts_class`
filtering** in the pool.

**What landed.**
- `calibration/strand_balance.py`: `StrandBalance` dataclass +
  `fit_strand_balance(substrate, strand_model)`. `κ_rna` = clamped
  `strand_model.p_r1_sense` (the live `StrandModels.exonic_spliced` MLE).
  `ρ_r_bb` = method-of-moments overdispersion for a **known mean** (the
  recovered symmetric estimator generalized off 0.5), pooled over all three
  views' spliced `(k_sense, n)` observations where `n > 0`. Degenerate cases
  (κ at a bound, all `n_i = 1`, under-dispersed) floor to `_BB_FLOOR`.
- `CalibrationResult` gains `kappa_rna` (III.7), validated `0 < κ_rna < 1`.
- `calibrate()` wires the fit (replacing the placeholder `ρ_r_bb`); the
  per-locus log line gains `kappa_rna`. gDNA deconvolution stays placeholder
  (PR 4). `initial_hyperparameters` no longer returns `ρ_r_bb`.

**Decision resolutions.** III.1 `κ_rna` = `StrandModel` (clamped); III.2 fit
`ρ_r_bb` from spliced (boundary + contained) observations; III.3 terminology
(no `k_plus`; raw `pos`/`neg` + `sense`/`antisense` from PR 2.5); III.5 ship
MoM (revisit empirically PR 7); III.6 guard = ≥2 observations (a mathematical
floor for variance, not a tunable) + count-weighting; III.7 add `kappa_rna`;
III.8 wired into `calibrate()`.

**Magic-number flag (Q6).** Three spec-derived constants in
`strand_balance.py`: `_BB_FLOOR = 1e-6` (doc 03 §8 — keeps the BB α/β finite;
clamps κ and ρ), `_RHO_R_BB_FALLBACK = 0.01` (doc 03 §7 init — mild non-zero
default when spliced data is too thin to estimate dispersion, so the strand
channel is not treated as a sharp discriminator), and
`_MIN_OBS_FOR_OVERDISPERSION = 2` (a variance-existence floor, not a tunable).
No fitted/heuristic cliff.

### Acceptance gate — results

- `tests/calibration/` → **54 passed** (incl. 7 new `test_strand_balance`).
- `ruff check` + `ruff format --check` clean.
- Full suite: **248 `NotImplementedError`** (unchanged mode/count, the
  `quant_from_buffer` PR 6 boundary), **+8 passing** (674); no new modes.
- Magic-number budget: 3 spec-derived constants (flagged above), 0 cliffs.
