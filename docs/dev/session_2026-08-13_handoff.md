# Session handoff — 2026-08-13

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    ⭐ Deliberately SHORT. The design, the findings and the numbers live in
    `docs/dev/region_boundary_sj_design.md`; this file is only *where you are and what to do next*.

---

## 0. WHERE THE STATE LIVES — read these, not this

| question | answer |
|---|---|
| where is the tool? | `ROADMAP.md` §0 |
| what do I do next? | ⭐ **§2 below**, then `ROADMAP.md` §2 item 0 |
| the per-transcript prior — design, findings, numbers | ⭐⭐ `docs/dev/region_boundary_sj_design.md` — **read it in full** |
| how does the RNA prior work today? | `EQUATIONS.md` §9b **and §9c** (the prior is ALREADY an additive per-component pseudocount; the weights are the design) |
| what must I not repeat? | `TRAPS.md` — two rules landed this session |

**Suite: 3,354 passed / 0 skipped / 10 xfail / 0 failed.** ⚠ Derived from 3,298, not adjusted: `+24`
`test_grouped_prior_update.py`, `+14`→`+8` `test_transcript_path.py`, `+6` `test_warm_start.py`, `+7`
composition gates in `test_result_schema.py`, and the file-parametrised meta-tests (`test_docs_boundary`,
`test_no_jargon_labels` gain one case per new file; `test_layering` one per calibration module;
`test_scripts_index` per script). `ruff` clean. C++ built. **Nothing committed — the owner drives commits.**

---

## 1. ⛔⛔⛔ THE THREE ERRORS THIS SESSION MADE. READ BEFORE TOUCHING AN ARM.

**1. A parameter accepted and silently dropped, which INVERTED a result.**
`pipeline._run_locus_em_partitioned` took `rna_prior_weight` and never passed it on — deleted while
removing an unrelated caller. Every allocation, however extreme, produced byte-identical output, and the
draft conclusion was *"the prior is too weak to matter."* ⭐ **A fire counter could not see it** (it
counted nonzero weights, not weights the solver read). What caught it was `--arm oracle_alloc_flip`: a
maximally WRONG allocation, which *must* move the answer, and did not. **Both defences must stay** — the
flip arm, and the counter that now watches the ESTIMATOR receive the array.

**2. A refuted mechanism re-derived from first principles and shipped as new.**
The first weighting function used the Jeffreys posterior MEAN as a density LOCATION, citing
`TRAPS: a-zero-count-is-a-measurement`. That mechanism is `ROADMAP.md` §4.1 graveyard row one, refused at
**+7,269 %** on the `g00` control. The trap names the DEFECT, not which quantity to repair — its shipped
repair put the half on the PRECISION. Lesson filed as `TRAPS: a-trap-names-the-defect-not-the-repair`.
⛔ **Grep the graveyard before adopting anything a trap seems to motivate.**

**3. Deviating from the owner's spec without saying so.** The first allocation introduced exclusivity
(never in the spec, hard-zeroed 38.7 % of transcripts), dropped the harmonic/geometric soft-min (the
whole noise-smoothing mechanism), and added a per-object `+½` (a magic number). It was 1.42× WORSE and
was deleted. ⭐ The rewrite is measured, not asserted, and the owner's own design won.

⚠ **Also: numbers measured through a broken lane must be re-measured, not carried.** The first
`oracle_alloc` figures were taken before fix (1) and were invalid; §0 of the design doc has the valid ones.

---

## 1b. ⭐⭐ WHAT THE 2026-08-13 (later) SESSION DID — read this before §2, which it partly supersedes

1. ✅ **The conserved sj mass is published** — `CalibrationResult.junction_conserved_mass`, a
   PROPERTY not a field. `ROADMAP.md` §2 item **0c0** has the reasoning and the numbers. §2.1's warning
   below is discharged: it was ~40 lines and **zero** fixture churn, because the thing that made the
   earlier attempt a mess was adding a stored FIELD, and a stored field is exactly what
   `prior_vs_oracle.OVERRIDE_FIELDS` makes wrong.
2. ⛔⛔ **STAGE 6 WAS BUILT AND IS REFUSED.** The soft-min-along-the-path family, 16 arms, both strata.
   `ROADMAP.md` §4 is the row, `TRAPS: an-upper-bound-is-not-an-estimate` the mechanism,
   `region_boundary_sj_design.md` §6a what survives. ⭐ **Do not re-derive it and do not read §2.1 below
   as still open** — it is done and the answer is no.
3. ⭐ **The next candidate is a SPARSITY mechanism, not a better mean**, and the reason is measured: 0.0 %
   of expressed transcripts are ever zeroed while 3,644 of 4,839 silent ones are not.
4. ⚠ **`sj_mass[2]` (item 0c) is still NOT started** — §2.5 below stands unchanged.

---

## 2. ⭐⭐⭐ WHAT TO DO NEXT

⚠ **§2.1 IS CLOSED — see §1b.** The rest stands.

**THE FOUNDATION IS PROVEN. THE WEIGHTING FUNCTION IS THE WORK.**

With true relative abundances as weights and the warm start zeroed, `g00 ss0.99 capture_off`:
transcript Σ|err| **4,275,988 → 696,900 (0.16×)**, gene **445,944 → 57,154 (0.13×)**, gene fp_mass
**23,529 → 0**. A FLIPPED allocation is **1.26×** worse. So good weights become good answers.

⛔ **Three limits travel with that number and must be quoted with it** — it is a CAPABILITY proof and not
a ceiling; true weights also hand over the true SUPPORT for free (a zero weight is exactly absorbing, and
4,579 of 8,750 transcripts are silent at `g00`); and it is ONE condition on a stratum the tool already
handles.

In order:

1. ⭐⭐⭐ **Design the weighting function (stage 6).** Inputs available today: the transcript PATH
   (`splice_graph.build_transcript_path`), per-object mass and count, and the three-way composition
   `(f_g, f_pos, f_neg)` on `CalibrationResult`. ⛔ Before designing, re-read the design doc §0's limits
   and error (2) above.

   ⭐⭐ **THE SPLICE-JUNCTION MASS IS AVAILABLE TODAY AND STAGE 6 IS NOT BLOCKED ON `sj_mass[2]`.** The
   conserved per-sj mass is exactly `mass_rna_junction * junction_mass_per_crossing` — verified on
   `g00 ss0.99 capture_off` to **9.1e-13**. Item 0c adds the per-STRAND split, whose consumer is artifact
   detection, not the weight.

   ⛔⛔ **BUT DO NOT READ `mass_rna_junction` AS A MASS — IT IS AN INCIDENCE COUNT DESPITE ITS NAME.** A
   fragment deposits `+1` on EVERY sj it uses, so measured on that condition there are **2.0719
   incidences per unit of conserved mass** (5,668,526 against 2,735,958.8), and the over-count is worst
   for multi-sj fragments. A weight that reads the count is wrong in proportion to how spliced a
   transcript is — which is exactly the axis a per-transcript prior varies over.
   ⚠ Publishing the conserved mass as its own field would remove the trap and is ~6 lines plus fixture
   updates in five files; it was started at the end of the last session, made a mess, and was reverted.
   Do it deliberately at the START of a session, not as a wrap-up nicety.
2. **Re-run stage 5 on the blind stratum** (unstranded × capture-ON) and one mid-gDNA condition before
   generalising the 84 %. `--arm oracle_alloc` is the arm; `--arm oracle_alloc_flip` is its falsification.
3. **Price what a REACHABLE weighting function could earn** — the ceiling, with the controls stage 5
   omits: a support-only ablation (inject the oracle only where truth is nonzero) to separate
   deconvolution from sparsity, and a gene-aggregated oracle to price the within-gene split.
4. **ψ's non-closure** — `ROADMAP.md` §2 item 0b. Owner ruled it a separate branch, after the prior.
5. **`sj_mass[2]`** — `ROADMAP.md` §2 item **0c**. The owner reversed the one-mass ruling because its
   premise changed (artifact detection is now a named consumer). Moves the payload schema digest ⇒ every
   scan cache rebuilds, ~6 min; bundle with the two dead banks §2.7 wants removed. ⛔ ONE schema change
   at a time or `rescan_panels.py` cannot attribute the delta.
5b. **The rest of the graph refactor** — design doc §5: the rename, the two sj CSRs, the boundary-event
   bitmask.
6. **The SJ output**, once `SJStrandTable`'s 1.576× disagreement with the accumulator is reconciled.
7. **The index duplicate map** — `ROADMAP.md` §2 item **0d**. An ALIAS map so a dropped duplicate is
   findable and the index is freestanding. ⚠ Needs an index rebuild but **NOT** a panel re-scan — it
   changes none of the four cache keys. Do not wait for item 0c.

---

## 3. THE ARMS THAT EXIST NOW (`scripts/design/quant_accuracy.py`)

| arm | what it varies | note |
|---|---|---|
| `oracle_alloc` | true weights + warm start ZEROED | ⭐ stage 5, the owner's spec, the best arm |
| `oracle_alloc_seed` | true weights, shipped seed kept | isolates the allocation alone |
| `oracle_alloc_flip` | a maximally WRONG allocation | ⛔ **the harness falsification — keep it** |
| `warm_uniform` | an equal seed instead of coverage | ⚠ measured a NO-OP (0.976–1.024×): the shipped seed is exonerated as a cause of the tool's error |

⚠ No `quant_accuracy` arm delta below the `base_reseed` floor (~6,700 fragments / 1,317 fields) is
attributable.

---

## 4. NEW INSTRUMENTS AND SURFACES

* `scripts/design/transcript_truth.py` — true per-transcript counts split by splicedness. 10 M records in
  **41 s**, read names only. ⛔ Splicedness from the read name's interval in SPLICED TRANSCRIPT
  coordinates, **never the CIGAR** (1,233,452 of 10 M truly-spliced fragments carry no `N` — the sj
  is in the unsequenced inner gap). Seven named gates, all passing.
* `splice_graph.build_transcript_path` — the transcript → ordered `(kind, id)` walk, transcription order.
* `EMConfig.warm_start` = `coverage` (shipped) / `prior` / `uniform`.
* `CalibrationResult` — `{gdna,rna_pos,rna_neg}_frac_{region,boundary}`. ⛔ Published as solved and NOT
  renormalised; the schema asserts `[0,1]` and deliberately not closure.
* `_apply_grouped_prior_update_test` — the C++ prior update reachable from Python; 24 gates, 12
  perturbations (11 caught, 1 proven inert), 0 holes.

**Two production bugs fixed:** the RNA prior was discarded at any locus with no gDNA candidate (measured
1 locus of 1,269 — rare, but the fix matters once weights are informative); and `chain_edge_deconv`
published `0` for both RNA strands, so the crossing axis's composition summed to `f_g` alone.
