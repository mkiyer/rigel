# PLAN v4 — THE MEASURED PRIOR, END TO END: what the accumulator carries, what calibration consumes, and the recipe for both implementation sessions

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** Derivations: `EQUATIONS.md` §2
    (support-conditional cancellation), §9d (the mixture). Rulings: `DESIGN.md` §0c.3, §0c.0d.
    Refutations: `ROADMAP.md` §4.2/§4.3. ⛔ v1–v3 (this file's history as
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

## PART C — THE CALIBRATION RECIPE (implementation session 2; rungs unchanged from v3, inputs now exact)

1. **Rung 1 — `total_abundance` assembly + wall mask.** Per slot: BOUNDARY = shipped banks (+
   certified spliced, INCIDENCE `count/(mu_r−1)`); REGION = side-selected `S/ℓ`, `E/ℓ`
   (precision-weighted where both exact — they are two counts of the same rate). Wall rule DERIVED:
   a side is exact iff its template distance `≥ w_max − 1`, `w_max` read from `deposited_lengths`,
   distances from the reach machinery at the component MINIMUM. Double-walled slots flagged
   not-model-free. Unit gates: hand-built fixtures; the certified-divisor rule asserted; the mask
   against the RANK-3 enumeration; perturb three ways.
2. **Rung 2 — `fit_total_landscape`** (reuse `DensityLandscape` on wall-exact totals, pre-pass-0):
   mode census `(rho_0, R, widths)` + per-slot `w_i`. Gates: synthetic bimodal/unimodal; the
   `rho_0`-vs-anchor CONSISTENCY gate (two independent estimators must agree); capture-OFF
   unimodality on a real cached condition.
3. **Rung 3 — the two-atom mixture** through `CompositionPriors`/`_location_term` (per-slot scalars
   `m_lo, m_hi, w`; `None` ⇒ bit-identical), behind `composition_reference` default `"structural"`.
   Gates: properness, matched tails, exact-zero degeneracy, regrid immunity, bit-identity under the
   default, flag-fires-when-flipped.
4. **Rung 4 — price on all 16**: per stratum, both zero controls, shuffle control (permute `w_i`),
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

* **Owner ruling 1 (before session 2's rung 4 defaults, not before session 1):** the small-exon
  length-composition channel — lift the 0.8.0 retirement for the support-based REGION-only form, or
  keep parked. Requires the fl-gap side panel either way.
* **Owner ruling 2 (settled by bundling):** `region_end_count` — RE-OPENED §4.3's condition is met
  (the total landscape is the consumer); included in Part B.
* **Session sequence:** this session = this document. Session 4 = Part B verbatim. Session 5 =
  Part C verbatim. Each starts with `preflight.py` + the suite, ends with the counts accounted.
* The panel's one-`eps`-per-probe simplification and real-panel variation: post-0.8.0.
