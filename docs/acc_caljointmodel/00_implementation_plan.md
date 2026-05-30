# Accumulator + Calibrator — Consolidated Implementation Plan

**Date:** 2026-05-29 (rev. 3 — decisions resolved + recovered-code map + PR series)
**Status:** Ready to begin **PR 0**. Decisions D1–D6 resolved (§4).

This is the single authoritative roadmap. It supersedes the *phase plans*
of [`../accumulator/01_implementation.md`](../accumulator/01_implementation.md)
and [`../caljointmodel/05_implementation_plan.md`](../caljointmodel/05_implementation_plan.md),
and consolidates the fractional-accumulator rewrite + calibration-v6
("jointmodel") into one plan.

The work is organized as a **series of PRs** (§7). Each PR has its own
detailed implementation doc under [`prs/`](prs/) that we critique and
implement individually. This master plan owns the cross-PR contract:
state, decisions, the recovered-code map, and the PR sequence.

> **Rev. 3 deltas.** (1) Your D1–D6 answers are now locked into §4 — most
> importantly the boundary model (left/right are *asymmetric, separately
> deconvolved* observations that **share one integer flux**; the doc's
> "½ half-split" is **abandoned**) and the strand-balance model
> (beta-binomial on **integer** sense/antisense counts from regions *and*
> boundary sides). (2) A new **Recovered-Code Map** (§5): the hardest
> pieces are recoverable from pre-burn git (`fc96902` + dangling blobs) —
> we harvest the mechanically-sound math and leave the cliff-ridden
> heuristics dead. (3) Phases → a 9-PR series (§7).

---

## 1. Authoritative cross-references

This file owns the **roadmap**; the design docs own the **model**.

| Topic | Authoritative doc | Reconciliation note |
|---|---|---|
| Accumulator structures, deposit algorithm, mass conservation | [`../accumulator/00_design.md`](../accumulator/00_design.md) | §7.3 example still lists stale `*_fl_hist`; shipped schema has no FL. Fix Phase 8. |
| Goals G2/G3/G4, substrate, principles | [`../caljointmodel/00_overview.md`](../caljointmodel/00_overview.md) | §5 keep/delete file list obsolete. |
| Generative model (NB count, BB strand, deterministic spliced RNA, Gamma exposure) | [`../caljointmodel/01_generative_model.md`](../caljointmodel/01_generative_model.md) | Current — **except** §7 boundary "½ half-split" is superseded by the side-attribution model (§4 D1). |
| Failure audit (why v5 burned) | [`../caljointmodel/02_failure_audit.md`](../caljointmodel/02_failure_audit.md) | Historical; §4/§5 burn lists already executed (with residue, §2.3). |
| Inference (E-step, M-step, closed forms, Newton) | [`../caljointmodel/03_inference.md`](../caljointmodel/03_inference.md) | Current — **except** §4 boundary "½ half-split" superseded (§4 D1); §3.2 `kappa_rna` source clarified (§4 D2). |
| Public API + pseudocount formulas | [`../caljointmodel/04_interface_contract.md`](../caljointmodel/04_interface_contract.md) | §2 substrate framing obsolete (→ `AccumulatorPayload`); §6 consumer signature references types that must be (re)built (§5). |
| Validation plan | [`../caljointmodel/06_validation_plan.md`](../caljointmodel/06_validation_plan.md) | Current. |

Conflicts are resolved by **updating the spec doc** (Phase 8), not by
silent deviation. Known conflicts are tracked in §3/§4.

---

## 2. Ground-truth state of the tree

> Anchored at `4124276 phase b`; line numbers drift — re-grep before editing.

### 2.1 DONE — accumulator substrate (keep, do not touch)

- Native [`accumulator.{h,cpp}`](../../src/rigel/native/calibration/accumulator.h):
  `Region` (uint32[4] contained), `Boundary` (float32[4] mass_left/right
  + uint32[4] flux), `Accumulator`, `AccumulatorSet` (per-ref, thread-safe
  merge). Channel `ch = (spliced?2:0) + (strand_pos?0:1)`.
- Scanner (Phase B): `set_regions`, fractional deposit of non-multimap
  fragments, post-scan merge, `result["calibration"]` payload. No dead
  v5 native code remains.
- [`scan_payload.py`](../../src/rigel/scan_payload.py) `AccumulatorPayload`;
  [`_accumulator.py`](../../src/rigel/_accumulator.py) façade;
  [`calibration/regions.py`](../../src/rigel/calibration/regions.py)
  region partition (built at index time → `regions.feather`; wired into
  the scan by `pipeline.py::_wire_calibration_regions`).
- Green: `tests/native/test_accumulator_spec.py` (20),
  `test_accumulator_payload.py`, `test_scanner_accumulator_integration.py`.

### 2.2 DONE — Phase-A burn skeleton (keep)

- [`calibration/__init__.py`](../../src/rigel/calibration/__init__.py):
  re-exports `CalibrationConfig` (from `config.py`), empty placeholder
  `CalibrationResult`, `calibrate()` → `NotImplementedError`.
- [`config.py::CalibrationConfig`](../../src/rigel/config.py#L212): the 4
  new knobs only; no v5 fields.
- [`calibration/strand_summary.py`](../../src/rigel/calibration/strand_summary.py)
  `StrandSummary` (library-global). CLI v5 flags removed.

### 2.3 RESIDUE — incomplete burn (PR 0 deletes all of this)

Reachable only after the `calibrate()` abort, so pure dead weight:

| Residue | Anchor | Action |
|---|---|---|
| v5 locus-prior EM driver | `pipeline.py::_run_locus_em_partitioned` (~546–956) | DELETE |
| 11 `prior_*` params + plumbing | `pipeline.py::_build_locus_meta` (~588–677) | DELETE |
| 11 `prior_*` output columns | `estimator.py::get_loci_df` (~844–854, 914–924, 953–963) | DELETE |
| dead `_calibration_strand_summary` | `pipeline.py` (~190–194) | DELETE |
| stale docstrings (`fl_models`, `assemble_priors`, `CalibrationScanPayload`) | `scoring.py:140,143`; `locus.py:12`; `pipeline.py:13–20,85,238` | REWRITE |
| stale doc ref `04_outputs.md` (nonexistent) | `config.py:217`, `calibration/__init__.py:11` | REWRITE → `04_interface_contract.md` |
| 3 legacy CLI tests | `tests/test_cli.py:312–337` | DELETE |
| legacy golden col `prior_ess_final` | `tests/test_golden_output.py:197` | REWRITE |
| old-schema `summary.json` consumers | `sim/locus_sweep.py:893–926`, `sim/analysis.py` | STUB (rewrite in Phase 8) |
| old-schema calibration test | `tests/test_sim_analysis.py:42–67` | STUB/DELETE |
| loci goldens with `prior_*` cols | `tests/golden/*_loci_df.{feather,tsv}` (21×2) | REGEN (Phase 8) |

### 2.4 MISSING — remaining work

Calibrator (E/M/exposure/outer loop); `CalibrationSubstrate` adapter; real
`CalibrationResult`; strand-balance model; per-region strand annotation;
exposure-weighted per-transcript lengths + locus-prior consumer; non-stub
`quant_from_buffer`. Most of this is **recoverable** — see §5.

### 2.5 Test reality

`pytest -x` dies at the first end-to-end test with the `calibrate()`
`NotImplementedError`. Every `run_pipeline` test is red; the 3 `test_cli.py`
tests are red; substrate tests are green. End-to-end green returns at PR 6.

---

## 3. Reality reconciliation (design docs vs code)

1. **Substrate is `AccumulatorPayload`** (4-channel region/boundary
   matrices), not doc 04 §2's `CalibrationAggregates`. Sufficient stats
   are channel reductions (§4 D1/D2).
2. **No per-region `kappa_rna`** in the live global `StrandModel` → §4 D2
   builds it from integer substrate counts.
3. **No `φ_{t,r}` allocator / `RegionArrays` in the live tree** — the
   codebase `LocusPartition` is an unrelated CSR scatter. But the
   allocator + effective-length math **exist in pre-burn git** (§5) and
   are recovered in PR 6.
4. **EM prior API is per-locus aggregate scalars**
   ([estimator.py:308](../../src/rigel/estimator.py#L308)
   `gdna_prior_count`/`rna_prior_count`), doc 04 §6 is per-transcript →
   bridge via recovered `assemble_multilocus_prior` (§4 D3, §5).
5. **Boundary "½ half-split" (doc 01 §7, doc 03 §4) is abandoned** in
   favor of the side-attribution model (§4 D1).
6. **Strand-ambiguous regions need a non-strand path.** Docs 01 §4.3 and
   03 §3.2 evaluate a per-region strand BB log-BF using `κ_rna[r]`, implicitly
   assuming every region has a defined sense. **AMBIG regions (both strands)
   have none** — the strand channel is unusable there. §4 D7 specifies the
   replacement (exclude from strand; impute by density-sweep + global
   fallback). Docs 01/03 get an AMBIG caveat now (one paragraph each) and a
   full §-rewrite in Phase 8.

---

## 4. Resolved decisions (D1–D6)

### D1 — Boundaries: separate, asymmetric, flux-shared (the "half-split" is gone)

**Why boundaries are first-class.** Boundaries — especially exon↔intron
boundaries — are where hybrid-capture enrichment is measured: gDNA splashes
out of captured exons into flanking sequence as boundary-crossing fragments.
So boundaries are deconvolved as their **own** sufficient-statistic set,
never merged into regions.

**Definitions (locked).** A boundary at index `b` sits between two regions:
- **left side** ↔ region `b-1` (the region to its left);
- **right side** ↔ region `b` (the region to its right).

The accumulator already stores, per boundary per channel:
- `mass_left` = fractional block-side mass contributed by blocks in the
  **left** region (`b-1`); `mass_right` = same for the **right** region (`b`).
  These two are **disjoint, asymmetric observations** and are
  **deconvolved separately** (an exon side and an intron side have very
  different gDNA/RNA composition).
- `flux` = **one shared integer** count of fragment-events crossing the
  boundary. Integer counts drive Bayesian shrinkage / statistical power;
  fractional mass drives density. **We use both.**

**What replaces the "half-split."** The design docs aggregated boundary
gDNA mass into a region as `½(M_g_L + M_g_R)`. That is **wrong** — the
accumulator already partitions mass by side. A region `r` instead receives
the *natural* side-attributed boundary contributions:

```
region r  ←  boundary[r].mass_right   (r is the right region of boundary r)
          +  boundary[r+1].mass_left  (r is the left region of boundary r+1)
```

(Doc 01 §7 / doc 03 §4 updated accordingly in Phase 8.)

**Likelihood coupling (the precise model).** Per boundary side:
- **Count / strand** use the **shared integer `flux`** (per channel),
  oriented to that side's region strand → NB count LLR + BB strand LLR
  → a per-side gDNA mixing proportion `π_g`. (Statistical power.)
- **Mass / density**: that side's fractional `mass_*` is the magnitude
  deconvolved by `π_g` into gDNA/RNA mass attributed to the adjacent
  region. (Density.)

**Recovered support.** The deleted [`boundaries.py @ fc96902`](#) already
implements this exact left/right model: `BoundaryTable` with
`left_region_*`/`right_region_*` counts, per-side boundary effective
length, and `left_region_index()`/`right_region_index()`. We adapt it to
read `AccumulatorPayload` instead of the old ledger (§5).

> **Remaining sub-decision (pin at PR 4 kickoff, not blocking now):** the
> exact form when a side's adjacent region is strand-ambiguous or
> non-genic (no defined RNA strand). Proposed default: such a side
> contributes to the count/density channels but its strand LLR falls back
> to neutral (`kappa_d = 0.5`), so it cannot manufacture RNA. Confirm at PR 4.

### D2 — Strand-balance model: beta-binomial on integer counts from regions **and** boundary sides

`kappa_rna` is the RNA strand mean; `rho_r_bb` is its beta-binomial
overdispersion. **Both are fit from integer (n_sense, n_antisense)
counts**, never fractions (0.9 from 9/1 and from 90/10 carry very
different statistical power; the fraction discards it).

**Observation extraction (locked).**
- **Regions:** every strand-**unambiguous** region with >0 contained
  spliced fragments contributes one `(n_spliced_sense, n_spliced_antisense)`
  pair, oriented relative to the region strand (and the library mode).
  From the substrate: spliced channels are `ch2` (spliced+`strand_pos`)
  and `ch3` (spliced+`strand_neg`); orient by region strand (flip for `−`).
- **Boundaries:** every boundary side whose adjacent region is
  strand-unambiguous contributes a `(n_spliced_sense, n_spliced_antisense)`
  pair from the boundary's **integer flux** spliced channels, oriented to
  *that side's* region strand. A side adjacent to an ambiguous region is
  skipped. (The left and right sides orient the *same* shared flux to
  *different* region strands — hence "per side.")

Pool all region + boundary-side integer pairs → fit a beta-binomial:
mean `kappa_rna` and overdispersion `rho_r_bb`, with count-based
uncertainty for a proper likelihood.

**New finding / task.** This needs **per-region strand annotation**
(`+ / − / ambiguous`), which the *current* minimal `regions.py` dropped.
It is recoverable: `regions.py @ fc96902` (signatures) +
`_arrays.py::transcript_strand_class` (`ts_class`). PR 1 restores region
strand/signature annotation; PR 3 builds the strand-balance model.

**Recovered support.** `strand_deconv.py @ fc96902`:
`_fold_pos_neg_by_transcript_strand` (the exact integer
sense/antisense constructor), `build_compartment_strand_counts`
(per region + boundary-side triples), `_log_beta_binom_pmf` (vectorized
integer BB log-PMF). `strand_balance.py @ fc96902` (symmetric BB
Method-of-Moments) is the skeleton — generalize its fixed-0.5 mean to a
fitted RNA mean. The live `StrandModel` supplies the library mode. **Avoid**
the `strand_deconv.py` reliability/screen half (log-BF cliffs, flag bits).

### D3 — Per-transcript exposure-weighted lengths → per-locus prior (recover, don't reinvent)

The bridge from calibration outputs to EM priors is **largely recoverable**
pre-burn code. The concept:

- Region effective length `region_eff_len = L − FL + 1`; a boundary is a
  point with `boundary_eff_len = FL` (spliced → RNA FL; unspliced → a
  gDNA/RNA mixture: gDNA fraction with gDNA FL, RNA fraction with RNA FL).
- Each region and boundary has an **exposure factor** (G4). Multiply
  effective length × exposure factor.
- **gDNA (per locus):** sweep regions/boundaries left→right, summing
  `exposure × eff_len`.
- **Per transcript:** collect the regions and the exon-compatible
  boundaries overlapping the transcript; sweep and sum `exposure × eff_len`
  → each transcript's exposure-weighted effective length → its Dirichlet
  pseudocounts. Aggregate per-transcript α to **per-locus**
  `gdna_prior_count` / `rna_prior_count` for the existing EM API.

**Recovered support (the crown jewels).**
- [`_exposure.py @ fc96902`](#) — `component_bp_weighted_exposure` (the
  per-transcript overlap sweep that bp-weights the region exposure array →
  exposure factor), `contained_exposure_clipped` (`L−FL+1`),
  `gdna_eff_len_for_loci`, `fractional_boundary_side_exposure`
  (boundary `eff_len = FL` via FL CDF). REUSE-AS-IS; depends only on the
  **live** `frag_length_model.py`.
- `prior.py::compute_component_exposure_table @ fc96902` — per-transcript
  `eff_len × exposure_factor`. REUSE.
- `adaptive_prior.py::_project_to_loci @ fc96902` — clean region→locus
  overlap-share allocator (locus-level φ). EXTRACT only.
- `locus_prior.py::assemble_multilocus_prior` (M-series blob) — clean
  per-locus → per-MultiLocus aggregation (`n_rna = max(0, n_obs − n_gdna)`).
  REUSE-AS-IS shape. **Avoid** its density/fusion terms.

### D4 — Pre-burn baseline & archived code are in-repo

Pre-burn snapshots live in git history; FL + calibration algorithms are
recoverable. We harvest via `git show <sha>:<path>` (cheat-sheet in §5).
The Phase-7 synthetic baseline is regenerated from a pre-burn SHA
(candidate: `fc96902 new calib system implemented`) into `scratch/preburn/`
**before** PR 0 lands.

### D5 — Locus-prior consumer: `assemble_priors` in `calibration/priors.py` ✓ (agreed)

### D6 — Phase 0 is its own PR ✓ (agreed) — PR 0 below.

### D7 — Strand-ambiguous regions: density-sweep imputation, not strand deconvolution (2026-05-29)

**The distinction (locked).** A region's 4-bit signature gives a *strand class*
that the calibrator must route on. NONE and AMBIG both lack a single transcript
strand but are **not** interchangeable:

- **NONE** (signature `0x0`): no transcript. gDNA-only expectation, unstranded.
  An arbitrary sense assignment is **safe** (neutral evidence). Included in
  both the density and strand models.
- **POS / NEG** (only `+` or only `−` strand bits): single-strand. Sense is
  well-defined → full count + strand + spliced deconvolution.
- **AMBIG** (both `+` and `−` strand bits — overlapping opposite-strand
  transcripts): every read is simultaneously *sense* for one transcript and
  *antisense* for the other, so **no valid sense split exists** and the strand
  channel is **unusable**. An arbitrary convention here would be **wrong** (it
  would manufacture a phantom strand signal). This is a genuine class, not an
  edge case — paralogs and nested/antisense genes produce it.

**Consequence — how AMBIG regions are deconvolved.** They are **excluded** from
(a) the D2 strand-balance fit (already specified) and (b) the E-step strand
log-BF. Their strand-agnostic density (`n_unspliced` / `n_spliced` / mass)
stays valid and is deconvolved by:

1. **Count / density channel** (NB) wherever exposure `ω_r` is identifiable.
2. **Boundary-sweep imputation** — propagate gDNA-density evidence from
   neighbouring strand-**decodable** regions/boundaries into AMBIG territory: a
   sequential **left→right and right→left sweep** from known/valid anchors,
   transferring evidence across each boundary with a reliability-weighted decay
   (recovered `boundary_sweep.py`, §5).
3. **Global fallback** — regions the sweep cannot reach fall back to the
   library-global gDNA density (a Gamma posterior; recovered
   `background_model.py` *mechanism*, §5).

**Where it lands:** **PR 4** (the region E-step routes by strand class; the
sweep + global fallback impute AMBIG and very-low-evidence regions). The work
is sizeable enough that PR 4 may split into PR 4a (decodable-region E-step +
exposure) and PR 4b (ambiguous sweep imputation) — decide at PR 4 kickoff.

**Magic-number caution (Q6).** The burnt sweep/background code concentrates
cliffs: `confidence_floor = 0.6` (`boundary_sweep.py`), `top_t_fraction = 0.01`
(the clip-top-X% seed in `background_model.py`), `damping = 0.5`,
`mu_floor = 0.05`. We recover the **mechanism** and **scrub the cliffs** — the
boundary transfer weight becomes evidence/overdispersion-based, the global
density a proper Gamma posterior with **no percentile clip**. Every surviving
constant is discussed at PR 4 kickoff before coding. **This also resolves the
deferred Q2 "clip-top-X%" thread** (PR01 §7): the principled answer to outlier
robustness is sweep imputation + a global Gamma posterior, not a percentile
cliff.

---

## 5. Recovered-Code Map (harvest, don't reinvent)

**Source commit: `fc96902` ("new calib system implemented")** for v6-draft
files; a few mature pieces are M-series **dangling blobs** (read via
`git cat-file -p <blob>`). All confirmed present.

| Need | Recover from | Verdict |
|---|---|---|
| Boundary↔region side model + per-side eff-len (D1) | `boundaries.py @ fc96902` (326 LoC) | REUSE (adapt to `AccumulatorPayload`) |
| Region/payload CSR containers + region strand class (D2) | `_arrays.py @ fc96902` (`RegionArrays.from_region_df`, `PayloadArrays`, `transcript_strand_class`) | REUSE (adapt `from_payload` to `AccumulatorPayload`) |
| Region partition signatures / strand annotation (D2) | `regions.py @ fc96902` (679 LoC; `classify_boundary_kind`, 4-bit signatures) | REUSE-WITH-CLEANUP |
| Integer sense/antisense constructor + BB log-PMF (D2) | `strand_deconv.py @ fc96902` → `_fold_pos_neg_by_transcript_strand`, `build_compartment_strand_counts`, `_log_beta_binom_pmf` | HARVEST functions; AVOID reliability/screen half |
| Beta-binomial MoM skeleton (D2) | `strand_balance.py @ fc96902` (196 LoC, symmetric mean-0.5) | REUSE-WITH-CLEANUP (generalize mean for RNA) |
| Per-transcript exposure-weighted eff-len (D3) ★ | `_exposure.py @ fc96902` (359 LoC) | REUSE-AS-IS |
| Per-component exposure table (D3) | `prior.py::compute_component_exposure_table @ fc96902` | REUSE (this function only) |
| Region→locus overlap allocator (D3) | `adaptive_prior.py::_project_to_loci @ fc96902` (~lines 390–478) | EXTRACT only |
| Per-locus → per-MultiLocus prior aggregation (D3) | `locus_prior.py::assemble_multilocus_prior` (blob `f7dc103b`) | REUSE-AS-IS shape |
| NB count overdispersion MoM (if needed) | `_kappa.py` (blob `3e7c406a`, 101 LoC) | REUSE-AS-IS |
| Closed-form Gamma-Poisson exposure (G4) | — (`exposure.py` is log-normal EB, NOT this) | WRITE FRESH per doc 03 §4 |
| `region_eff_len = L−FL+1`, gDNA/RNA FL | **LIVE** `frag_length_model.py::compute_all_transcript_eff_lens`, `FragmentLengthModels` | CALL DIRECTLY |
| RNA strand beta posterior + library mode | **LIVE** `strand_model.py` (`n_same/n_opposite/minor_rate_*`) | CALL DIRECTLY |
| Boundary-sweep imputation for AMBIG regions (D7) | `boundary_sweep.py @ fc96902` (`_forward_scan`/`_reverse_scan` L↔R propagation, `run_boundary_sweep`) | RECOVER MECHANISM; scrub `confidence_floor=0.6` + transfer pseudocounts → evidence/overdispersion-based weight |
| Global gDNA-density fallback (D7) | `background_model.py @ fc96902` (`BackgroundModel` Gamma posterior `rho_off`) | RECOVER MECHANISM; **drop** `top_t_fraction=0.01` clip, `damping=0.5`, `mu_floor=0.05` → proper Gamma posterior, no percentile cliff |
| Excess→contained-gDNA predictor (D7) | `boundary_model.py @ fc96902` (`predict_contained_gdna_from_excess`) | EXTRACT (the closed-form predictor only) |
| Test scaffolds | `test_exposure.py`, `test_strand_deconv.py` (count/PMF tests), `test_strand_model.py`, `test_per_locus_gdna_mass.py`, `test_boundary_model.py`, `test_boundary_sweep.py` @ fc96902 | REUSE selectively |

**AVOID (the burned cliffs — do not resurrect):** `integration.py`
(fusion), `calibration_iteration.py` (two-state E/M density staging),
`density_model.py`, `adaptive_prior.py` (except `_project_to_loci`),
`latent_states.py`, `exposure.py` (log-normal EB shrinkage),
`strand_deconv.py` reliability/screen half. These concentrate ~390 tunables.

> **Reclassified for D7 (2026-05-29):** `boundary_sweep.py`,
> `background_model.py`, `boundary_model.py`, and `density_observation.py`
> were previously on the pure-AVOID list, but the **sweep + global-fallback
> mechanism** they implement is exactly what D7 needs for strand-ambiguous
> regions. We **recover the mechanism and scrub the cliffs** (see the D7 rows
> above). The *cliffs* in those files remain forbidden; only the structural
> propagation/posterior math is harvested.

**Adapter strategy.** The recovered exposure/boundary/prior code is built
on `RegionArrays`/`PayloadArrays`/`BoundaryTable`. Rather than rewrite it,
PR 1 builds a thin **`AccumulatorPayload` → RegionArrays/PayloadArrays/
BoundaryTable** adapter, so the recovered math runs largely unchanged.

**Harvest cheat-sheet.**
```bash
git show fc96902:src/rigel/calibration/_exposure.py          # φ overlap allocator + eff-len  ★
git show fc96902:src/rigel/calibration/boundaries.py         # per-boundary side model        ★
git show fc96902:src/rigel/calibration/_arrays.py            # RegionArrays/PayloadArrays
git show fc96902:src/rigel/calibration/strand_balance.py     # BB MoM skeleton
git show fc96902:src/rigel/calibration/strand_deconv.py      # fold + BB log-PMF (harvest fns)
git show fc96902:src/rigel/calibration/prior.py              # compute_component_exposure_table
git show fc96902:src/rigel/calibration/adaptive_prior.py     # _project_to_loci (extract)
git cat-file -p f7dc103b85f18a961c6d6f1d4d1f9f0afd7ebd79     # mature locus_prior.py (structure)
git cat-file -p 3e7c406ab430735d95b28421a9b178b80f0d4d98     # _kappa.py NB MoM
```

---

## 6. Cross-cutting conventions

- **One PR per item in §7**; each ends at its acceptance gate. Next PR
  starts only after the previous gate is green.
- **Build between C++-touching PRs** (`pip install --no-build-isolation -e .`).
  C++-touching PRs: **PR 2.5** (accumulator strand & spliced-flux correction)
  and possibly **PR 6** (only if D3's per-locus aggregation proves
  insufficient). All other PRs are Python-only.
- **No deletion outside the current PR's file list.**
- **No heredocs in the terminal** (workspace convention); multi-line
  diagnostics in `scripts/debug/`.
- **Goal-directed:** every calibrator line traces to G2/G3/G4.
- **Sufficient statistics, not fragments**; per-fragment iteration in the
  calibrator is a bug.
- **FL is not a calibration channel** — it stays per-fragment in the EM
  scorer (`scoring.py`, `frag_length_model.py`). (FL *is* used for
  effective lengths in the prior consumer — that's not a calibration
  channel.)
- **≤ 8 numeric literals** in the calibrator core, each annotated with the
  spec it derives from (doc 03 §8). Recovered code must be scrubbed of
  cliff constants on the way in.
- **Pause-and-discuss rule (hard).** Before introducing *any* new magic
  number, heuristic, or tunable parameter, STOP and discuss with the user
  — do not forge ahead. The legacy system was burned for ~91 magic numbers
  / ~25 cliffs; the goal is to not repeat that. Prefer fitted
  overdispersion (NB `φ`, Beta-Binomial `ρ`) over outlier-clipping cliffs
  for robustness. Keep the module count low (≤ 25).
- **Region definition (locked).** Every genome position carries a 4-bit
  signature {exon_pos, exon_neg, intron_pos, intron_neg}; a region is a
  maximal interval of constant signature, and a region and its neighbours
  always differ in signature (adjacent equal-signature segments are
  merged). Boundaries sit at signature transitions.
- **Strand class — NONE ≠ AMBIG (locked, §4 D7).** Derived from a region's
  signature: **NONE** (signature 0, intergenic) has no transcript — unstranded
  gDNA, so an arbitrary sense is **safe**; it stays in both density and strand
  models. **POS/NEG** (one strand) have a well-defined sense. **AMBIG** (both
  `+` and `−` bits — overlapping opposite-strand transcripts) has **no valid
  sense**: it is excluded from strand deconvolution and recovered by density +
  boundary-sweep imputation + global fallback. Never conflate NONE (no strand)
  with AMBIG (both strands) — an arbitrary sense for AMBIG manufactures phantom
  strand signal.
- **Harvest before writing.** For any new component with a §5 entry, start
  from the recovered source; justify deviations.

---

## 7. The PR series

| PR | Title | Touches | Build | Doc |
|---|---|---|---|---|
| **0** | Burn the residue *(Step 1)* | Python | no | [`prs/PR00_burn_residue.md`](prs/PR00_burn_residue.md) |
| **1** | Reorganize + substrate adapters *(Step 2)* | Python | no | `prs/PR01_reorganize.md` (next) |
| **2** | Calibrator scaffold (types, validation, placeholder) | Python | no | `prs/PR02_scaffold.md` |
| **2.5** | Accumulator strand & spliced-flux correction (D2/D7) | **C++**+Py | **yes** | `prs/PR025_accumulator_strand_flux.md` |
| **3** | Strand-balance model (D2) | Python | no | `prs/PR03_strand_balance.md` |
| **4** | Calibrator core: E-step (G2/G3) + exposure (G4) | Python | no | `prs/PR04_estep_exposure.md` |
| **5** | M-step + outer loop (working calibrator) | Python | no | `prs/PR05_mstep_outer.md` |
| **6** | Integrate: exposure-weighted lengths + locus prior (D3) + `quant_from_buffer` | Python (±C++) | maybe | `prs/PR06_integrate.md` |
| **7** | Validate (paralog, hybrid-capture, sweeps, armis2) | Python | no | `prs/PR07_validate.md` |
| **8** | Final cleanup (goldens, ruff, docs, postmortem) | mechanical | no | `prs/PR08_cleanup.md` |

PR docs are written **just-in-time** (one ahead) so each reflects the
state the previous PR actually produced. PR 0's doc is written now.

### PR-by-PR intent (detail lives in each `prs/` doc)

- **PR 0 — Burn the residue.** Delete every §2.3 remnant; fix stale
  docstrings/doc-refs; capture the pre-burn baseline (D4). No behavior
  change; pipeline still aborts at `calibrate()` but nothing dead hangs
  off it. See [`prs/PR00_burn_residue.md`](prs/PR00_burn_residue.md).
- **PR 1 — Reorganize + substrate adapters.** Target package layout
  (`substrate.py`, `result.py`, `calibrate.py`, `priors.py` skeletons);
  recover the `AccumulatorPayload → RegionArrays/PayloadArrays/BoundaryTable`
  adapter (from `_arrays.py`/`boundaries.py`); restore **region
  strand/signature annotation** (D2 prerequisite, from `regions.py @ fc96902`);
  wire `calibrate(substrate, strand_model, config)` with real args (still
  raising); index-alignment guard (calibration array order == `region_df`).
- **PR 2 — Calibrator scaffold.** Real `CalibrationConfig` (reconcile
  `phi_floor` 1e-9 vs doc 03 1e-6), `CalibrationResult` (doc 04 §5 +
  `__post_init__`), `CalibrationSubstrate` (3 views: contained / left /
  right; channel reductions per §4 D1/D2; boundary→region projection),
  exceptions, placeholder `calibrate` returning zero-mass/unit-exposure.
- **PR 2.5 — Accumulator strand & spliced-flux correction (C++).** Rename
  `exon_strand` → `align_strand`; orient spliced fragments at deposit
  (`sense = align_strand == sj_strand`) so the spliced channels are
  motif-relative sense/antisense (valid in AMBIG, D7); make boundary flux
  **per-side** (`flux_left`/`flux_right`) so intron-skipping spliced reads do
  not create false exon-intron flux (D2); rework the substrate to expose raw
  genome-strand counts (orient downstream). Build + payload-schema change; no
  on-disk index change. See [`prs/PR025_accumulator_strand_flux.md`](prs/PR025_accumulator_strand_flux.md).
- **PR 3 — Strand-balance model.** Beta-binomial on integer
  (n_sense, n_antisense) from regions + boundary sides (§4 D2); harvest
  `_fold_pos_neg_by_transcript_strand` + `_log_beta_binom_pmf` + the
  `strand_balance.py` MoM skeleton (generalized mean). Produces
  `kappa_rna`, `rho_r_bb` with count-based uncertainty. **Excludes AMBIG
  regions/sides** (no valid sense, §4 D7); NONE included as neutral.
- **PR 4 — Calibrator core.** E-step deconvolution (doc 03 §3) **routed by
  strand class (§4 D7)**: strand-decodable regions (POS/NEG/NONE) use NB count
  LLR + BB strand LLR (PR 3's model) → `π_g` → soft masses (G2/G3) on the 3
  substrate views; closed-form Gamma-Poisson exposure (G4, doc 03 §4, fresh);
  boundary side-attribution per §4 D1. **AMBIG regions skip the strand LLR**
  and are imputed by the recovered **boundary sweep** (L↔R density propagation
  from decodable anchors) + **global-density fallback** (recover
  `boundary_sweep.py`/`background_model.py`/`boundary_model.py` *mechanism*,
  scrub cliffs — §5, §4 D7). May split into **PR 4a** (decodable E-step +
  exposure) and **PR 4b** (ambiguous sweep imputation) — decide at kickoff.
  Tests: paralog/hybrid-capture sanity, AMBIG-overlap-locus imputation, mass
  conservation, vectorized==scalar.
- **PR 5 — M-step + outer loop.** `rho_0`, `eps_s` (closed forms);
  `phi`, `rho_d_bb`, `rho_r_bb` (1-D Newton + moment warm start);
  `pi_g_prior` update; outer loop with the mass-change monotonicity
  sentinel (`CalibrationConvergenceError`). Synthetic E2E recovery.
  Write `scripts/debug/dump_calibration_state.py`.
- **PR 6 — Integrate.** Recover `_exposure.py` + `compute_component_exposure_table`
  + `_project_to_loci` + `assemble_multilocus_prior` (§5); build the
  region↔transcript overlap sweep; implement `assemble_priors`
  (doc 04 §6); aggregate per-transcript α → per-locus prior counts (D3);
  rewrite `quant_from_buffer` to wire payload → calibrate → priors → EM;
  end-to-end smoke green.
- **PR 7 — Validate.** Paralog (t1≈t2≈500, natural pass), new
  hybrid-capture scenario, synthetic sweeps vs the PR-0 pre-burn baseline,
  armis2 corners, numerical-stability stress. Watch the FL-free paralog
  risk flag (doc 02 §2.5). Write `validation_report.md`.
- **PR 8 — Final cleanup.** Regenerate goldens; ruff; magic-number audit
  (≤ 8); rewrite sim calibration reports to the v6 schema; reconcile the
  spec docs (§1/§3/§4 conflicts); update `CLAUDE.md` +
  `.github/copilot-instructions.md`; postmortem; archive retention.

---

## 8. Risk register

| PR | Risk | Mitigation |
|---|---|---|
| 0 | A "dead" block has a hidden live caller | `git grep` each symbol; `pytest --collect-only` after each deletion batch |
| 0 | Pre-burn baseline window already lost | Snapshot from `fc96902` via `git worktree` before PR 0 lands (D4) |
| 1 | Recovered adapter (`AccumulatorPayload`→RegionArrays/Boundary) mis-maps indices across ref seams | Index-alignment guard + boundary-projection property test |
| 1 | Region strand annotation regressed when minimal `regions.py` was created | Recover signature/strand logic from `regions.py @ fc96902`; test against synthetic +/−/overlap geometry |
| 3 | Recovered strand code drags in cliff constants | Harvest only the named functions; magic-number scrub on entry |
| 4 | FL-free three-leg paralog rescue weaker than FL-bearing v5 | PR 7 paralog regression must verify; fallbacks in doc 02 §2.5 |
| 4 | Boundary side strand undefined for non-genic side (D1 sub-decision) | Default neutral `kappa_d=0.5` strand LLR; confirm at PR 4 kickoff |
| 4 | AMBIG sweep imputation drags back the burnt cliffs (`confidence_floor`, top-X% clip, `damping`, `mu_floor`) | Recover mechanism only; scrub every cliff; each surviving constant discussed at PR 4 kickoff (§4 D7) |
| 4 | Strand-class routing missed → AMBIG silently oriented (phantom strand signal, the v5-style paralog pathology) | `ts_class` carries the 4-way distinction from PR 1/2; substrate test asserts AMBIG flagged; PR 3/4 key strand use off `ts_class` |
| 5 | Mass-change increases between iterations (EM violation) | `CalibrationConvergenceError`; never silently continue |
| 6 | Per-transcript→per-locus aggregation loses resolution paralogs need | If paralog fails, escalate to per-transcript EM API (C++ change) |
| 7 | armis2 regression, no clear cause | `scripts/debug/dump_calibration_state.py` (PR 5) |
| 8 | Golden regen masks a real regression | Regenerate only after PR 7 sign-off; diff for unexpected non-`prior_*` changes |

---

## 9. Sign-off

Decisions D1–D6 are resolved (§4). **Next action: critique
[`prs/PR00_burn_residue.md`](prs/PR00_burn_residue.md), then implement
PR 0.** Each PR is a standalone, reviewable unit; its doc is written one
ahead and critiqued before implementation.
