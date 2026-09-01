# PLAN v4 — THE MEASURED PRIOR, END TO END: what the accumulator carries, what calibration consumes, and the recipe for both implementation sessions

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** Derivations: `EQUATIONS.md` §2
    (support-conditional cancellation), §9d (the mixture). Rulings: `DESIGN.md` §0c.3, §0c.0d.
    Refutations: `ISSUES.md`'s CLOSED / REFUSED record. ⛔ v1–v3 (this file's history as
    `PLAN_spike_slab_reference.md`) are superseded and survive in git; their defects are recorded in
    v3's header one commit back. **v4 adds the owner's start/span/end taxonomy (2026-08-21), the
    what-do-we-need-to-know inventory, and the two implementation recipes.** This is the handoff for
    BOTH implementation sessions; delete sections as they land, per the MOVE rule.

---

## PART A — THE BIG PICTURE (finalised)

### A.0 The strategy (unchanged from v3, one paragraph)

ψ inherently needs a proper prior (`TRAPS: no-prior-means-haldane`); the correct response is to FIT
it from composition-free observables before any solve runs, not to conjure a reference. TOTAL
abundance is such an observable. A Poisson-kernel landscape fitted on measured totals BEFORE pass-0
(no circularity: its inputs are counts and lengths) is bimodal under capture and delivers, per
slot: the enrichment state `w_i` (the slot's own total projected onto the modes), the span `R` (the
mode ratio — never the 2.6–3.6×-under-reading anchors), and `rho_0` (the depleted mode, gated for
consistency against the pooled ¬`mrna_active` anchors). The pass-0 gDNA prior is then the two-atom
mixture `(1−w_i)·Spike(rho_0) + w_i·Spike(rho_0·R)`, capped by the slot's own total, delivered as a
density on the λ grid through the existing `_location_term` door — nothing picked; the posterior
picks; silent evidence stays wide and reports itself near-undetermined, which is
`solvability_audit`'s doctrine. Refits sharpen from the slots pass-0 earned.

### A.1 The owner's taxonomy — four relations, and it is COMPLETE

Given a REGION `r` and an accepted fragment PATH, exactly one of: no covered base in `r` (no
deposit); its FIRST covered base is in `r` (**START**); its LAST covered base is in `r` (**END**);
or it covers every base of `r` with neither endpoint inside (**SPAN**, strict). A CONTAINED
fragment is START ∧ END in the same region (increments both; never SPAN). Completeness holds
because region bounds sit at every exon endpoint, so an annotated splice can jump a region entirely
(case 1) but can never partially cover one without an endpoint in it.

Expectations under a uniform field (per component, `ℓ` the region length, `d` the template distance
past the relevant end):

    E[S_r] = rho·ℓ          exact for every w — EXCEPT within d < w−1 of the template's DOWNSTREAM end
    E[E_r] = rho·ℓ          exact for every w — EXCEPT within d < w−1 of the template's UPSTREAM end
    E[V_r] = rho·E_f[(w−ℓ−1)₊]        a pmf functional per component — NOT model-free, by design
    E[C_r] = rho·E_f[(ℓ−w+1)₊]        the shipped contained bank; V and C are COMPLEMENTARY
                                       truncations (region_density_derivation T3: together they
                                       estimate rho·(1−f(ℓ+1)))

⭐ **Correction to the taxonomy's one motivating worry**: S is NOT blind at small INTERNAL regions —
a 30 bp mid-transcript exon takes starts at the full rate `rho·30` for every fragment length (the
fragment continues through the sj). The blindness of S/E is at TEMPLATE ENDS only, and it is
opposite-ended, which is exactly why the pair closes it: **use the side whose wall does not bind;
precision-weight where both are exact**. Double-walled regions (template shorter than `2·w_max`)
remain honestly not-model-free and are excluded from fits, never papered over.

⭐ **The differential `S_r − E_r` is a new measured signal and its meaning is derivable**: away from
walls its expectation is 0 for every component; near a template end it is the net fragment-end flux
— zero for gDNA everywhere, and at an RNA terminus its magnitude is that terminus's RNA start-rate.
Recorded as a diagnostic (a terminus-level measurement candidate); NO consumer is wired in this
plan, and none may be added without its own falsification.

### A.2 "What do we need to know?" — the inventory, decided

| need | consumer (today / this plan) | carrier | verdict |
|---|---|---|---|
| **TOTAL abundance** `T_i`, composition-free | the pre-pass-0 total landscape → `(rho_0, R, w_i)` → the measured prior; density-model background comparisons; any surviving enrichment ratio | BOUNDARY: shipped exact banks (`inv_length_sum`, `sj_inv_length_sum`, + certified spliced at `count/(mu_r−1)`, INCIDENCE by design). REGION: **S and E over ℓ, side-selected by the wall rule** | ⭐ the one genuinely NEW need; S shipped, E new |
| **Incidence counts, per strand** | the strand model (Beta-Binomial needs integers per strand); ψ's mixture numerator | `region_contained_count[,2]`, `boundary_unspliced_count[,2]`, spliced, sj | KEEP unchanged — ⛔ `contained` is not redundant with S/E: it is the UNSPLICED-CONTAINED mixture ψ deconvolves (`M_i`), a different population from the spliced-inclusive totals |
| **Conserved mass** | `mass_per_crossing` → the prior assembler's conserved fragment COUNT (the 0.179→0.0027 assembler repair; capture over-call 1.149→1.019); sj artifact detection | `boundary_unspliced_mass`, `boundary_spliced_mass`, `sj_mass` | ⭐ KEEP — the owner's question answered: mass and density interconvert only THROUGH an opportunity model, and not needing that model at the assembler was the measured win. Mass answers "how many fragments", density "how many per bp"; the EM eats fragments |
| **Length pools / deposited_lengths** | fl pmf estimation; `w_max` for the wall rule (read from the support end, never chosen) | shipped | KEEP |
| **`region_contained_inv_opportunity_sum`** | today: only `RegionGeometry.inv_abundance`'s REGION face (consumer: the frozen currency policy) + instrument moments + spec gates | shipped | ⚠ CANDIDATE-RETIRE after the policy re-contrast; do NOT delete in this schema change (one thing per change; the spec gates it) |
| **The length-composition channel** (small/terminal exons: `(S,E,C,V)` + the two pmfs solve the split; `det ∝ E_g(ℓ)−E_r(ℓ)` self-mutes at equal pmfs) | none — ⛔ **owner scope ruling pending** (collides with the 0.8.0 length-channel retirement; ancestor refused at the zero control) | would consume C and V | PARKED; gated on the ruling + the fl-gap side panel + the `g00` gate |

---

## PART B — ✅ EXECUTED 2026-08-21, code AND re-scan; both suites re-certified 16/16 + 16/16

Landed: `region_start_count` widened to genome-strand columns; `region_end_count` and
`region_span_count` added (spec first — six gates verified failing, then green — then C++, gated by
the auto-extending parity gate, FIRED once by a deliberate C++ perturbation); export, payload
(`BANK_AXES`), substrate, and the oracle's sum-to-full bank list carry all three; the F-arm's locus
projection strand-sums. `DESIGN.md` §3.1 documents the banks and their ledger.

⭐⭐ **THE RE-SCAN'S THREE INCIDENTS, recorded so the next schema change does not rediscover them:**

1. **A WIDENED bank is a CHANGED bank, not an added one.** `rescan_panels` correctly REFUSED the
   first pass (shape `(n,) → (n,2)` fails byte-identity); the designed remedy is
   `--expect-changed region_start_count`, which keeps the gate's teeth on every other bank. ⚠ The
   refused first pass had already written some payloads before the cross-condition check fired, so
   the second pass reported four `g00 _main`s as "expected-changed did NOT move" — a re-run over
   already-converted payloads, not a data problem; the CERTIFIER is the validity authority, not the
   rescan report.
2. **The g00 HOLD-OUT has two arms.** `panel.py cache` (via `pass0_vs_oracle`'s prewarm) holds
   zero-gDNA conditions out ENTIRELY: after a from-scratch oracle rebuild the four `g00` dirs simply
   do not exist, and after a partition deletion their `rna_pos`/`rna_neg` are never rebuilt. The
   remedy: `rescan_panels --conditions <the four g00s>` (it does NOT hold g00 out) plus a direct
   `pass0_vs_oracle.ensure_rna_strand_cache(...)` call per g00 condition. Worth wiring into
   `panel.py cache` properly — an owner-visible hygiene item.
3. **`RIGEL_SCRATCH` must be EXPORTED or the instruments overflow `/tmp`.** The rebuild left ~22 GB
   of per-condition prewarm splits in `/tmp/rigel_pass0_oracle` before the variable was exported
   (regenerable; flagged to the owner, not deleted). ⚠ Separately discovered: a PRIOR session's
   instrument self-test left **73 GB** in `/tmp/rigel_ladder_ceiling/_self_test` — a named defect
   (a self-test must clean its work-dir) awaiting an owner-approved cleanup.

## PART B (the original recipe, for the record)

Landed: `region_start_count` widened to genome-strand columns; `region_end_count` and
`region_span_count` added (spec first — six gates verified failing, then green 66/66 — then C++,
gated by the auto-extending parity gate, which was FIRED once by a deliberate C++ perturbation);
export, payload (`BANK_AXES`), substrate, and the oracle's sum-to-full bank list all carry the three;
the F-arm's locus projection strand-sums. `DESIGN.md` §3.1 documents the banks and their ledger.

## PART B (the original recipe, for the record)

**One schema change, one re-scan (~1 h), three banks — bundled so no later ruling needs another
re-scan.**

### B.1 The banks

| bank | shape | rule | note |
|---|---|---|---|
| `region_start_count` | uint32[n_regions] → **uint32[n_regions, 2]** | `+1` at `region_of_pos(first_base)`, column = align strand | WIDENED to genome-strand columns (the counts convention: "counts carry both columns"); the ledger invariant becomes `sum over both columns == deposited` |
| **`region_end_count`** (NEW) | uint32[n_regions, 2] | `+1` at `region_of_pos(last_base)` | `last_base` is already computed at `accumulator.cpp:603`; today `region_of_pos(last_base)` is evaluated only inside the unspliced guard at `:719` — the new bank needs the lookup on SPLICED paths too |
| **`region_span_count`** (NEW) | uint32[n_regions, 2] | per SEGMENT: `+1` on every region strictly interior to the segment's crossed range, excluding the path's own start and end regions | strict cover: every base covered, neither endpoint inside. ⚠ its OWN consumer is conditional (A.2's parked channel) — it is bundled because the marginal cost inside this re-scan is ~zero and bundling converts the pending ruling into a calibration-only change |

### B.2 The invariants (new QC gates, written FIRST, each perturbed)

    sum(region_start_count) == sum(region_end_count) == qc.deposited      ← the ledger, now twice over
    per region:  region_contained_count[r] <= min(S_r, E_r)   (a contained fragment is both)
    per region:  V_r == 0 wherever ℓ_r >= w_max − 1            (spanning needs w >= ℓ+2)
    global:      sum(V) == sum over fragments of (#regions strictly spanned)   ← from the spec side

### B.3 The plumbing checklist — the integer path is HAND-plumbed; trace it against `region_start_count`, never against the float bank

* `accumulator.h`: struct-free arrays beside `region_start_count_` (init `:266`-area, merge
  `:869`-area, accessor beside `:480`); update the invariant comment at `:476`/`:499`.
* `accumulator.cpp`: deposit sites (§B.1); the spliced-path `last_base` lookup.
* `tests/native/_accumulator_reference.py`: the executable specification, SAME change, byte-identity
  gate; its `region_start_count` line and prose updated (it is the file the C++ is gated on).
* export: `bam_scanner.cpp:2224`-area + prop `:2910`-area.
* payload: `scan_payload.py` — `ADDITIVE_AXES` (`:317` pattern), dataclass fields beside `:514`,
  ctor `:766`-area, docstring `:124`/`:311`; ⚠ `n_regions` is read from `region_start_count.shape[0]`
  (`:619`) — still valid after widening.
* substrate: `substrate.py:149`/`:200` pattern — expose as `region_start_count[,2]`,
  `region_end_count[,2]`, `region_span_count[,2]` (int64 on the way out, per convention).
* `deposit_digest` moves BY DESIGN (`scan_cache.py:220-298`; gates `tests/test_scan_cache.py`) — the
  key now sees deposit changes; every cache is correctly refused.
* dtypes: integer banks ⇒ bit-identical merges, NO `FLOAT_BANK_DEPOSITS` entry in
  `rescan_panels.py:78-85`.

### B.4 The rebuild procedure (from §8 of the model-free plan, verified sites)

1. `pip install --no-build-isolation -e ".[dev]"`; suite; **fire the new QC gates by perturbation**
   (e.g. deposit end at `first_region` — the ledger must catch it).
2. `rescan_panels.py --jobs 8` — a pure ADDITION lands in `only_new` with `differing` empty (its
   self-test already gates exactly that); ⚠ the WIDENED start bank appears as removed+added — the
   same delta required on every condition. ⚠ Its `PAYLOADS` is 4; `rna_pos`/`rna_neg` need one
   `panel.py cache` pass after (`--index` override for the test chromosome:
   `~/Downloads/rigel_runs/test_reference/idx`).
3. Re-certify (`calibration_oracle.py`, both suites), `preflight.py`, suite, AND the instruments
   (`TRAPS: a-green-suite-hid-five-dead-instruments`).
4. Re-derive the collected count; account every move.
5. ⚠ Estimated ~1 h total; re-time with `scan_profile.py` before quoting.

### B.5 Explicitly NOT in this change

No deletions (`region_contained_inv_opportunity_sum` stays, gated); no float banks; no boundary-axis
changes; no consumer wiring (that is session 2); no `mass` changes.

---

## PART C — THE CALIBRATION RECIPE (implementation session 2) — REORDERED 2026-08-21 per the owner: VALIDATE FIRST, REPLACE SECOND, BUILD THIRD

⭐⭐ **The owner's sequencing, adopted as the rung order**: ① validate the new `total_abundance`
measure against certified truth until it is beyond doubt; ② replace the existing abundance readings
with it where a consumer exists; ③ only then build the landscape and the prior on top.

⭐ **Expectation-setting for ② on the CURRENT panels, stated before it is measured
(`TRAPS: a-single-level-panel-cannot-see-a-constant`'s discipline):** both suites hold equal gDNA/RNA
fragment-length distributions BY DESIGN, so the COMPOSITION-DEPENDENCE half of the abundance fix is
a null here — do not expect it to move the 0.8.0 metric, and do not read "no change" as "no defect".
What IS visible on equal lengths: the TRUNCATION half (`w` against `ℓ` — the 11.6× at 98 bp exons,
the 58.6 % zero-bank exon hops) and the WALL half. ⚠ And the §4.3 record stands as a warning about
②: a better level channel inside the wrong consumer measured WORSE (the accidental mute), so every
replacement is priced per consumer, never assumed. **Fully exerting the composition half needs the
two-condition fl-gap SIDE PANEL** (config-only: `gdna_ladder.yaml` carries independent
`simulation.frag_*` and `gdna.frag_*` blocks, deliberately identical today; ⛔ a SIDE panel, never a
ladder rung — equal lengths are the ladder's forcing function). Build it when ② needs its verdict.

1. ✅ **Rung 1 — DONE AND VALIDATED 2026-08-21. Its content has MOVED and is deleted from here:** the
   ruling to `DESIGN.md` §3.1a-ii (the four settled points, and why this module is NOT §4.3's refused
   drop-in), the derivation to `EQUATIONS.md` §2.3b (`A_start(w|d)`, the exactness condition, the
   field-free pair), the two lessons to `TRAPS.md`
   (`two-estimators-of-one-rate-weight-the-field-differently`,
   `state-the-population-rule-do-not-inherit-it-from-a-table`), and the numbers to `ISSUES: measured-prior-rung-4`
   rank 3. Built: `rigel.calibration.total_abundance` (layer 3, no consumer) +
   `splice_graph.build_mature_wall_distances`; gated by `tests/calibration/test_total_abundance.py`
   (23 cases, 16 perturbations all fired) and `scripts/design/total_abundance_audit.py` (15/15).
2. ⏳ **Rung 2 — REPLACE PER CONSUMER: ENUMERATED, TWO CONSUMERS LANDED AND PRICED, the rest ranked.**
   ⛔ **A prior session called this rung BLOCKED by a refusal scoped to ONE consumer, and the owner corrected it. The
   error is worth keeping**: §4.3 is scoped to ONE consumer (the currency channel's `enrichment_ratio`,
   substituted UNMASKED and without side selection, judged on that policy's deliverable), and
   "readers of `inv_abundance`" is not "places total abundance is used" — one member versus 240.
   * The census that established this counted **240 sites** (8 SWAP / 20 NEEDS-MEASUREMENT in `src/`);
     its enumeration was deleted 2026-09-01 as a stale-on-arrival artefact — ⛔ **re-derive it with
     `module_census.py` rather than quoting it**. What survives is the RULE: only a TOTAL-wanting site
     (or a gDNA level at a structurally pure-gDNA object, where total ≡ gDNA) is a candidate at all.
   * LANDED behind `CalibrationConfig.background_abundance` (`"contained"` default ⇒ bit-identical):
     `fit_intron_background` (→ ψ) and `measure_background` (→ the refit floor). The PAIR changes, not
     the estimator.
   * PRICED at the estimator level on all 32: a tie off capture, a 1.8–4.3× repair under it, and a
     ~2 pp nascent-contamination cost on RNA-bearing pools. `ISSUES: measured-prior-rung-4` carries the table.
   * ⭐⭐⭐ NEXT, and it is rung ③'s work rather than rung ②'s: `capture_eff_length.py:194` +
     `priors.py:353` — the expression that produces the 0.0951. A composition-free total can supply
     that reference, but only through the LANDSCAPE's enriched mode; a pooled-anchor drop-in is right
     at `g00`/capture-OFF and wrong under capture, which is rank 2's already-recorded shape.

3. ⭐⭐⭐ **Rung 3 — THE `AbundanceLandscape` — DESIGNED 2026-08-21 (the four plans below), building now.**

   ### 3a. The fit — `fit_abundance_landscape`, and every design decision with its reason

   New module `src/rigel/calibration/abundance_landscape.py`, **layer 5** (imports `landscape` SIDEWAYS,
   `total_abundance` layer 3 DOWN, `signature` layer 0). ⭐ **It reuses `fit_landscape` verbatim** —
   that estimator is deliberately component-agnostic (its own docstring: "nothing here knows which
   component it is fitting"), and every hard-won decision in it transfers: zero-native Poisson kernels,
   knn population resolution, the Laplace one-pseudo-region floor, the derived grid.

   The call, and why each argument is what it is:
   * `count, exposure, model_free = total_abundance.region_counts_and_exposure(...)` — **REGIONs
     only.** Boundaries are excluded for the SAME reason the gDNA hyperprior excludes them (owner,
     2026-07-27): they cross rather than contain, and a landscape describes the field over the genome,
     which regions tile and zero-width boundaries do not.
   * `count = counts`, `eff = exposure` — the VALIDATED pair (rung 1, 32/32), Poisson-native.
   * `mass = counts` — `_grid`'s top is `max log10(mass/eff)`; for a TOTAL the observed density IS the
     ceiling (nothing sits above what was measured), so mass ≡ count and the grid derives itself.
   * `var = zeros` — ⭐ **a DIRECT measurement has no deconvolution ambiguity**, so `_reliability`
     returns 1 everywhere and the weights are honestly flat. This is not a shortcut: the gDNA
     hyperprior's `Var(log f_g)` weight exists because ITS training data came out of a solve; ours
     did not.
   * `anchor = (count == 0)` among the selected — semantic parity with the hyperprior's zero-count
     anchor (weight already 1; the flag is documentation).
   * SELECTION: `model_free & (exposure > 0)`. ⛔ Double-walled regions are excluded by the mask
     itself, which is what rung 1 built the mask FOR.
   * ⚠ `_KNN_SCALE`/`_S0` are inherited and **have only ever been validated on gDNA-shaped data** —
     `landscape.py`'s own warning. Every result this landscape reports must carry that caveat until
     rung 5 prices it.

   ### 3b. The MODE CENSUS — no threshold anywhere, and the precedent it follows

   `CalibrationDiagnostics.from_prior` already finds modes as "local maxima of the fitted log-density
   curve, tallest first" — the precedent. The census extends it with BASINS and refuses thresholds:
   * Every interior local maximum is a mode; the grid is PARTITIONED into basins at the interior
     minima between adjacent maxima; each mode carries `(log_rho, basin_mass, width, lo, hi)` with
     width = the mass-weighted std WITHIN the basin (derived, no constant).
   * **DEPLETED** = the basin CONTAINING the pooled anchor rate
     `log(Σcount/Σexposure over intergenic model-free regions)` — the same composition-free pool
     `fit_intron_background` uses; fallback (no anchors, a toy) = the largest-mass basin.
   * **ENRICHED** = the largest-mass basin strictly ABOVE the depleted one; `None` ⇒ unimodal.
   * `rho_0 = exp(depleted.log_rho)`; `span_R = exp(enriched − depleted)` (1.0 when unimodal).
   * ⭐ **WHY NO SIGNIFICANCE THRESHOLD IS NEEDED**: the knn smoothing already suppresses combing (its
     whole job), and a phantom maximum above the bulk carries ~zero basin mass, so it yields
     `w_i ≈ 0` and a near-1 depleted share — HARMLESS, and the capture-OFF gate MEASURES that instead
     of a constant asserting it. Masses carry the verdict continuously.
   * **`w_i`** (per model-free REGION): the slot's own Poisson kernel times the fitted density,
     normalised on the grid, integrated over the ENRICHED basin — an honest posterior responsibility,
     `≡ 0` when unimodal, `NaN` where not model-free.
   * **ANCHOR CONSISTENCY** (the rung's named gate): `|depleted.log_rho − anchor_log_rho|` must sit
     within the depleted mode's own fitted width — two independent estimators of one level, and the
     tolerance is the density's own statement of its resolution rather than a chosen number. The gap
     is reported in nats either way.

   Gates (FIRST, verified failing, ≥3 perturbations): synthetic bimodal (two Poisson populations,
   decisive separation → 2 modes, rho_0/R within grid resolution, `w_i` separates the populations,
   bracketed); synthetic unimodal (enriched None, span_R exactly 1, w ≡ 0); anchor-consistency flip
   (pool the anchors 3 decades off → the flag must flip); the mask respected (a not-model-free region
   must contribute nothing); basins PARTITION (Σ basin_mass = 1); determinism; `None` under 2 training
   regions.

   ⛔⛔ **THE FIRST REAL FIT CORRECTED ONE GATE'S WORDING, 2026-08-21 — record it before it misleads.**
   The rung's spec said "capture-OFF unimodality on a real cached condition". The real capture-OFF
   total is BIMODAL (R ≈ 96 at `g50 ss0.99`): expressed exons sit decades above the gDNA background in
   any TOTAL, capture or no capture — the T-conflates-expression confound arriving exactly where §A.0
   predicted. ⭐ **The invariant that actually holds, exactly, is the NON-EXON one: capture-OFF
   `w[intergenic] = w[intron] = 0.0000` (both, to four decimals) while the anchor gap reads 0.003
   nats (`rho_0` 0.0511 vs the known 0.051).** Under capture the intron class lifts to only 0.0019
   while exons take 0.35 (exposure-weighted). So the census DOES cleanly separate the classes on both
   settings — but "unimodal off capture" was the wrong words for it. ⭐⭐ **Consequence for rung 4,
   settled now**: the landscape ALONE cannot distinguish expression-bimodality from
   capture-bimodality, so the two-atom gDNA prior must consume the landscape's `(rho_0, w_i)` TOGETHER
   with the already-measured enrichment detector (`DESIGN.md` §6b.1's anchor ratio: 0.98 without
   probes, 113–114 with) — off capture the mixture must degenerate to the spike per §9d.4 limit ①
   REGARDLESS of what the total's upper modes are doing, because those modes are RNA. ⚠ Also seen and
   harmless: a tiny mode AT the grid's resolution wall (log_rho ≈ −16, mass ≤ 0.008) — the zero-count
   population's own kernels; the anchor-containing rule never selects it.

   ### 3c. Integration at calibration INIT — this session, default byte-identical

   In `calibrate()`, immediately after the wall-mask block that `background_abundance="measured_total"`
   already builds: fit the `AbundanceLandscape` when `CalibrationConfig.abundance_landscape = True`
   (default **False** ⇒ nothing runs ⇒ bit-identical). It REFUSES without the wall inputs exactly as
   `measured_total` does — no silent fallback. Its outputs this session: `_debug["abundance_landscape"]`
   and a new `InjectedCalibrationPriors.abundance_landscape` field (default None, so every toy path is
   untouched). ⛔ NOTHING reads it in the solve this session — the fit is priced as a QC/injection
   object first, exactly as the npmle was, so the flag flip below is an A/B and not a leap.

   ### 3d. ✅ THE NPMLE RETIREMENT — EXECUTED 2026-08-21. The recipe is deleted; its record MOVED.

   ⛔ The step-by-step recipe that stood here has been REMOVED rather than stamped, per the MOVE rule.
   Where it went: the outcome and the measured justification to `ISSUES: measured-prior-rung-4`; the
   lesson to `TRAPS: measure-a-default-flip-before-you-write-it`; the ruling on which landscape outputs
   a consumer may read to `DESIGN.md` §3.1a-iii.

   ⭐ **THE GATE PASSED**: `rename_identity.py --check` reads `✅ BIT-IDENTICAL — quant digest and every
   array's content unchanged`, which is what had to be true, since nothing the npmle fed was ever read
   by the solve.

   ⛔⛔ **THE ONE PLACE THIS PLAN WAS WRONG, kept because it is the transferable part**: step 2 said
   "flip the default". Flipping `abundance_landscape` to `True` against the flag's REFUSAL broke **65
   callers** (24 failed + 41 errors) — unit and toy fixtures have no wall arrays and never wanted a QC
   panel. The refusal's reason was *"rather than fitting on unmasked totals"*, and the alternative was
   never to fit unmasked but NOT TO FIT, so the missing-walls case now skips with a WARNING and yields
   `None`. `background_abundance` keeps its refusal (it feeds ψ) and a gate asserts the asymmetry.


   ### 3e. The pass-0 REFERENCE — rung 4's design, so it is not re-derived

   The two-atom mixture (`EQUATIONS.md` §9d.4) through the LOCATION DOOR:
   * per slot: `m_lo = clip(rho_0·E_g,i/M_i)`, `m_hi = clip(rho_0·R·E_g,i/M_i, ≤ the slot's own §9d.3
     cap)`, `w_i` from the landscape (REGION slots: own responsibility; BOUNDARY slots: to be ruled —
     the flanks' `w` or their own exact banks);
   * the door widens: `CompositionPriors.location` accepts `(m, 3)` — `(m_lo, m_hi, w)` — and
     `_location_term` becomes `−log[(1−w)·((1−m_lo)·f + m_lo·(1−f)) + w·((1−m_hi)·f + m_hi·(1−f))]`,
     which degenerates EXACTLY to the shipped scalar at `w ∈ {0,1}` or `m_lo = m_hi`, and to zero at
     `m = ½` — the same neutral-row property the scalar term has;
   * strength is §9c.1's and does NOT change (one pseudo-observation; the location never enters τ_λ);
   * config: `composition_reference: str = "structural"` ABSORBS the `structural_reference` bool
     (`"none" | "structural" | "measured"`), default bit-identical;
   * `g00`: anchors → `rho_0 = 0` → the §9c.1 closed-end clamp carries (the coordinate note in §9d.4's
     limit ②);
   * gates: properness, matched tails, exact-zero degeneracy, regrid immunity, bit-identity under the
     default, flag-fires-when-flipped — then rung 5's pricing.

4. **Rung 4 — the two-atom mixture** — design above (3e); build after the landscape census is read on
   real conditions.
5. **Rung 5 — price on all 16**: per stratum, both zero controls, shuffle control (permute `w_i`),
   the evidence split, the refuted single-mode arms beside it, the truth-hot-unprobed CONFOUND read
   (T conflates enrichment with expression; the failure direction is PERMISSIVE — measure it), and
   the refit-1 gDNA-landscape mode census. Then the owner rules on the default, and the
   three-policy re-contrast runs on the honest reference.

**Consumers map (so no bank is dark):** S, E → rung 1 (and their strand columns give the strand
model a spliced-inclusive OPTION later, no re-scan). V → parked channel only (dark until the
ruling; its QC invariants are its interim consumer). `S − E` → diagnostic census only. Everything
else unchanged.

---

## PART D — OPEN RULINGS AND SEQUENCING

* **Owner ruling 1 — ✅ DECIDED 2026-08-21: the small-exon length-composition channel is
  DEFERRED for Part C.** It stays parked with the 0.8.0 retirement; `region_span_count` sits dark
  behind its ledger invariants until a post-0.8.0 revisit. No Part C rung may consume `V` or the
  contained bank as a composition channel.
* **Owner ruling 2 (settled by bundling):** `region_end_count` — RE-OPENED §4.3's condition is met
  (the total landscape is the consumer); included in Part B.
* **Session sequence:** this session = this document. Session 4 = Part B verbatim. Session 5 =
  Part C verbatim. Each starts with `preflight.py` + the suite, ends with the counts accounted.
* The panel's one-`eps`-per-probe simplification and real-panel variation: post-0.8.0.
