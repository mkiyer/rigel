# NEXT SESSION — the AbundanceLandscape: head-to-head vs the NPMLE, the bandwidth question, and the transfer-variance audit

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** It is a handoff: read it, do the
    work, MOVE what settles into the permanent docs, and DELETE this file in the same session
    (precedent: the previous `NEXT_SESSION.md` was deleted per its own instruction).

## Where the tree is (2026-08-21, all committed on `fragment-length-gold-standard`)

⭐ **The measured-prior thread's rungs ①–③ are DONE**; the standing baseline is
**0 failed / 3,661 passed / 9 xfail / 3,670 collected** (this file is the `+1` over 3,660 — one
jargon-gate case). Landed this session, each gated and perturbed:

* **Rung ① — `rigel.calibration.total_abundance`** (layer 3): the composition-free per-slot TOTAL —
  START/END banks over `ℓ`, side-selected by the derived wall rule; validated 32/32 with a FIELD-FREE
  decisive arm (`total_abundance_audit.py`). `region_counts_and_exposure` is the pooled-consumer entry
  point. Rulings: `DESIGN.md` §3.1a-i/ii; derivation `EQUATIONS.md` §2.3b.
* **Rung ② — two background consumers swapped** behind `CalibrationConfig.background_abundance`
  (default `"contained"` ⇒ bit-identical): estimator-level a tie off capture, 1.8–4.3× repair under it;
  end-to-end ~1 % (`TRAPS: a-better-estimator-inside-a-weak-consumer-moves-nothing`).
  `background_reference.py` was CONVERGE-AND-DELETED (a duplicate pool with no live consumer).
  The consumer census (240 sites, the four FORMS) is `CENSUS_total_abundance_consumers.md` here in the
  sandbox.
* **The fl-gap SIDE panel** (`flgap_rna_long` / `flgap_rna_short`, ±170 bp realised, g50 × ss × capture)
  is built, cached, certified (COMPOSITION 8/8; FIELD withheld capture-OFF — a real ±1.2 % intron bias
  in the contained-opportunity model whose SIGN FLIPS with gDNA length; its own ticket, do not fold in).
  Verdict: pooled background rates are gap-INSENSITIVE (`TRAPS:
  a-pooled-rate-cannot-see-a-short-object-factor`); the composition mechanism lives at EXONS.
* **Rung ③ — `rigel.calibration.abundance_landscape`** (layer 5): `fit_abundance_landscape` reuses
  `landscape.fit_landscape` VERBATIM on the wall-exact totals (`var = 0` — a direct measurement has no
  deconvolution ambiguity) + a THRESHOLD-FREE mode census (basins at density minima; masses carry every
  verdict; depleted = the basin CONTAINING the pooled intergenic anchor rate; tolerance = the mode's own
  width). Integrated at calibrate INIT behind `CalibrationConfig.abundance_landscape` (default False ⇒
  bit-identical; True REFUSES without `mature_walls`/`boundary_reach`; injectable via
  `InjectedCalibrationPriors.abundance_landscape`). Gates:
  `tests/calibration/test_abundance_landscape.py` (16), instrument
  `scripts/design/abundance_landscape_census.py` (13/13, `--json` dumps census+curve+rug).
* **The 40-condition census**: every contaminated ladder row anchor-consistent (gaps 0.003–0.063 nats
  off capture — `rho_0` 0.0511 vs anchor 0.051 at `g50 ss0.99`; 0.25–0.65 under capture); enriched
  basin CLAIMED by exons (w̄ₑ 0.35 vs intron 0.002 at g50-ON); capture-OFF non-exon w ≡ 0.0000.
  ⛔ **Capture-OFF totals are BIMODAL (expression, not enrichment)** — the plan's "capture-OFF
  unimodality" gate was corrected to the NON-EXON invariant, and **rung 4 must pair `(rho_0, w)` with
  the enrichment detector** (anchor ratio 0.98 off / 113–114 on). The atlas over all 40 fits is
  published: `https://claude.ai/code/artifact/921a3642-35b2-4505-8303-efb991a5b457`.

## THE WORK — the owner's three questions, in their order

### ① HEAD-TO-HEAD: the NPMLE landscape vs the AbundanceLandscape. Which is more accurate?

Plot the two PER CONDITION, side by side, and score both against a TRUTH:

* The NPMLE arm: `DensityNPMLE.fit(mass_global, eff_global, bandwidth=config.npmle_bandwidth)` —
  exactly `calibrate.py`'s call (`region_gdna_geometry` supplies the pair; ⛔ note it is the forbidden
  `mass/eff_gdna` form, total mass over the gDNA opportunity). The AbundanceLandscape arm: the shipped
  builder off the same cached payload.
* ⭐ **THE TRUTH IS DERIVABLE FROM THE ORACLE PARTITIONS AND NEITHER ESTIMATOR'S MACHINERY**: per
  region, the true total density = `(S_gdna + S_mrna + S_nrna) / ℓ` from the origin payloads' own
  START banks (untruncated, pmf-free — the same truth arm `total_abundance_audit.py` uses). The true
  DISTRIBUTION is that value's exposure-weighted histogram; the true `rho_0` is the gdna partition's
  pooled intergenic rate; the true enrichment split is per-class.
* ⛔ **"More accurate" needs the right yardstick, and EMD is NOT it alone** — `landscape.py`'s own
  header records that EMD is MONOTONE IN SMOOTHING at every reference and must not re-select anything.
  Score mode LOCATIONS and basin MASSES against the truth's (the N5 shape discipline), plus the anchor
  gap each fit produces, per stratum, both zero controls, never pooled.
* ⚠ The two fits answer on different substrates (npmle: every slot with support, regions+boundaries
  pooled; landscape: wall-exact REGIONs). Say so on the page rather than pretending one axis.

### ② THE SPIKINESS / BANDWIDTH QUESTION — how do we know the bandwidth is optimal?

The owner reads the atlas as "spiky, more modes than expected" (7–16 local maxima per condition).
What is known, so it is not re-derived:

* The spikes are SUB-STRUCTURE within the two basins that matter; the census is threshold-free BY
  DESIGN (masses carry verdicts; a wiggle owns ~no mass). The far-left wall mode is the zero-count
  population's own kernels — expected, never selected.
* The kernel width machinery is `landscape.knn_widths` — READ ITS DOCSTRING FIRST: the per-region
  Poisson widths collapse 41× across the density range, so population resolution comes from
  k-th-nearest-neighbour spacing (`k = √n`), scaled by `_KNN_SCALE = 0.5`, floored at one grid step.
  `_KNN_SCALE` was SHAPE-selected against a truth-validated reference (below 0.5 the landscape combs;
  above it the modes merge); ⚠ it has only ever been validated on gDNA-SHAPED data — the standing
  caveat every census page repeats. `_N_GRID = 260` and `_WIDTH_BINS = 12` are computational budgets.
* The npmle's `npmle_bandwidth = 0.15` decades is a FIXED width — a different philosophy (one
  constant), and part of what ① compares.
* ⭐ Principled checks that need no new constant: (a) the truth comparison from ① — does the true
  total-density distribution HAVE that sub-structure? (capture probes are per-exon, so real
  multi-modality within the enriched basin is plausible, not necessarily noise); (b) SPLIT-HALF
  stability — fit on disjoint halves of the regions and count reproducible modes (the documented
  mode-count-vs-sample-size slope target is ≈ −1.0); (c) a `_KNN_SCALE` sweep priced on (a)+(b), per
  stratum. ⛔ `TRAPS: no-magic-numbers` — if a different scale wins, it lands with its selection
  evidence, not as a taste adjustment.

### ③ THE TRANSFER-VARIANCE AUDIT — NPMLE's old role in pass-0 message propagation

⛔ **Start from what is RECORDED, then verify it rather than assume it**: `calibrate.py`'s enrichment-
NPMLE comment says its σ²_transfer role is *"RETIRED (that was a density-uniformity proxy, invalid
under capture, and identically 0 in pass-0); the solver now derives σ²_transfer itself"*, and this
session verified `DensityNPMLE.project` has ZERO callers in `src/` (only `test_npmle.py`;
`messages/variance.py:498` cites the proxy formula in a docstring only). The audit:

1. TRACE where transfer/premise variance actually comes from today, per policy —
   `messages/variance.py` (`composition_logvar`, `residual_level`, `splice_in_frame_logvar`),
   `currency.py`'s `premise_logvar`/`hop_logvar` fit, `relay.py`'s uses. Write the mechanism down.
2. Confirm the retirement (no live npmle path into any σ²) — or find the residue.
3. THEN ask the owner's question: could the AbundanceLandscape supply a BETTER transfer variance?
   It has the machinery the npmle's `project` had, done right: the per-region posterior over the grid
   (`w` already integrates it over the enriched basin — the posterior VARIANCE of `log ρ` given
   `(count, exposure)` is one line on the same grid), fitted on wall-exact totals instead of the
   forbidden pair. ⚠ But it is TOTAL density — the same composition-vacuity that retired the npmle's
   role applies; be explicit about what a total's variance can and cannot say about a gDNA hop.
4. ⛔ **Policy work is FROZEN pending the measured-prior re-contrast** — the audit and a design note
   are this session's deliverable; a src change into `messages/` is NOT, unless the owner rules
   otherwise mid-session.

## Standing next steps beyond the three questions (ranked in `ROADMAP.md` §1 rank 3)

* Rung 4 — the two-atom mixture through the WIDENED location door (`(m_lo, m_hi, w)`; design in
  `PLAN_measured_prior.md` §3e), paired with the enrichment detector; then rung 5's pricing.
* The NPMLE retirement (plan §3d) — ① may produce exactly the evidence that ships it; THE gate is a
  bit-identical deliverable.
* Outstanding: `calibration_vs_oracle.py --background-abundance` on the two fl-gap arms; the
  contained-opportunity intron bias ticket; `/tmp/rigel_ladder_ceiling` (73 GB) still on disk.

## Session mechanics

* FIRST: `python scripts/design/preflight.py` and `python -m pytest tests/ -q`
  (baseline above; ANY failure is a regression; re-derive the collected count, never adjust).
* `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`;
  `OMP_NUM_THREADS=1`; `python -m pytest`, never bare pytest; ⛔ export `RIGEL_SCRATCH` in every
  instrument invocation.
* Falsification first, verified failing; perturb the fixed code; surgical reverts.
* The owner drives commits. The atlas artifact updates by republishing the same file path with its URL.
