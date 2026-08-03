# Fragment length — the audit, and the cleanup

    Status: ⚠ **HISTORY, as of 2026-08-03** — C0, C1, C2 and C2.6 have all landed and the
      fragment-length work is DONE: the anchor's error against truth is +0.00 % mean / +0.02 % sd.
      ⭐ **Its "C3 is next" is CORRECT AGAIN** and C3 is now `TODO.md` rank 0 — `LEDGER.md` B4
      measured why: calibration fits from the POOLS, and the RNA pool is what C3 corrects.
      §§1–3 are the history of how three definitions of fragment length came to be live at once.

    Status: C0, C1, C2, C2.6 LANDED (2026-07-31 / 2026-08-01). C3 is next; C4 is independent.
            The audit itself is complete and is kept as the record of WHY.
    Trigger: P2's A/B (`LEDGER.md`) — the length likelihood is the first consumer to use the FL models
    as a DISCRIMINANT rather than as a divisor, and it exposed that they disagree with each other.

## §0 ⭐ ROADMAP — WHERE WE ARE AND WHAT IS NEXT

**This file is the critical path.** It was reached from `SOLVER_OBSERVABLES_PLAN.md`, which is where the
calibration-side work lives; that plan's P2 is **blocked on this file's C2 and C3**.

    WHERE THE WORK LIVES
      branch  fragment-length-gold-standard        (⚠ NOT merged to main)
      commit  d045d820  "the prior is a fragment count, and L becomes the gold standard"
              = P1 + P2(off) + C0 + C1, one commit -- LEDGER.md has the per-step narrative
      main    is one commit BEHIND. Fast-forward when you want it:
              git checkout main && git merge --ff-only fragment-length-gold-standard

    THE STANDING CHECK -- run this before and after every step
      source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
      pip install --no-build-isolation -e ".[dev]"      # after ANY src/rigel/native/ change
      python -m pytest tests/ -q                        # ⛔ `python -m`, NOT bare `pytest` -- see below
      ruff check src/ tests/ scripts/                   # ⚠ never `ruff format scripts/`

    ⛔ THE INVOCATION MATTERS. Bare `pytest` does not put the repo root on sys.path, so
    `tests/calibration/test_fl.py`'s `import tests.native._accumulator_reference` raises and the
    suite reads 1808 / 23. That 23rd failure is an artefact of how it was run, not a regression --
    and under "a 23rd failure is a regression" it reads as one. Measured 2026-08-01.

⚠ **1835 / 21 is the baseline as of C2.6** (was 1824 / 21 at C2, 1809 / 22 before it). All 21 are
`test_golden_output`, moved **three** times now -- P1's units fix, C2's FL models and C2.6's `L`.
**A 22nd failure is a regression.** Do not regenerate the goldens to make the count look better -- they
move again at C3.

⭐ **The count went DOWN by one: `TODO.md` §7's `test_nrna_double_counting[g20_n0_s100]` now passes**, at
**0 counts against a limit of 25** (it leaked ~30). ⚠ Do not close §7 on this -- a negative control is
one-sided (trap 19) and this was not the change's target.

| | step | status | gate |
|---|---|---|---|
| **C0** | Prove the accumulator's `L` before promoting it | ✅ **DONE 2026-07-31** | 333,684 exhaustive configs vs an independent set oracle; 7 perturbations, 6 caught, 7th a proven no-op |
| **C1** | Give the accumulator an unconditional length histogram (`deposited_lengths`) | ✅ **DONE 2026-07-31** | `Σ = Σ node_start_count = qc.deposited`, enforced in the accumulator AND at the payload door; parity gate auto-covers it; 6 perturbations |
| **C2** | switch every consumer to the accumulator, delete the scanner histogram, re-point `rigel report` | ✅ **DONE 2026-08-01** | all six sub-gates met; `payload_schema_digest` never moved; report keys identical and now sum to the library. `LEDGER.md` |
| **C2.6** | search EVERY fragment's unsequenced gaps for introns, not only unspliced ones. Spec: **`docs/SPEC_GAP_INTRONS.md`** | ✅ **DONE 2026-08-01** — `LEDGER.md`'s C2.6 entry | ⭐ **G-sd MET: +27.0 % → +1.98 %.** G-tail 85 % met (0.00909 → **0.00137** above the true 713 bp ceiling; `dropped_too_long` **280,558 → 38,309**) and ⛔ **the residual is LOCALISED, not excused** — it is D3, measured at **98.5 %** of what is left. G-gdna control **bit-identical** |
| **C2.7** | ⭐ **D3 — a mate gap holding MORE THAN ONE annotated intron keeps only the first cut** | ⛔ **the only known mechanism left in the tail.** Measured, not hypothesised: emitting every intron in the gap takes the residual to **0.00002** and `dropped_too_long` to **389** | needs the per-gap unanimity test to compare intron **sets** rather than one intron — that is the real work, and M3 did not do it |
| **C3** | §8.1(b): divide each pool by its own opportunity. ⭐ **The formula is DERIVED and PROVEN** — `docs/JUNCTION_OPPORTUNITY.md` | ⭐ **UNBLOCKED** — the anchor is clean. ⛔ **But re-run `JUNCTION_OPPORTUNITY.md` §3.2's θ control first**: its +59.8 % sd and +0.3 % mean were measured against the contaminated anchor AND against a pool that has since moved to −9.58 % / −22.46 % | ⚠ **G3 SPLITS**: G3a (mean) targets the measured **+0.3 %** achievable with true θ; G3b (sd, ceiling) was **C2.6's gate** and is met |
| **C4** | Gate the length discriminant on both pools having data | independent — any time | the `strand_evidence` analogue |
| **C5** | Delete `ScanCache.fl_rna_counts`; verify D7 | ✅ **DONE 2026-08-01, inside C2** | grep clean; D7 asserted end to end and perturbation-proven (`tests/test_d7_transcript_eff_lengths.py`) |
| **→** | **re-run `scripts/design/length_likelihood_ab.py`**, then decide `CalibrationConfig.length_likelihood` | after C3 | `SOLVER_OBSERVABLES_PLAN.md` §6.4 |
| **→** | **regenerate `tests/golden/` ONCE** | after C3 | ⚠ 21 goldens are stale from P1 and will move again at C2/C3 — regenerate **once**, at the end, and ⛔ regenerate twice and diff (`TODO.md` §7) |

✅ **THE DECISION THAT BLOCKED C2 IS MADE** (owner, 2026-08-01) and it was **neither** option offered:
the five per-fragment splice-type counts are **scanner QC** and stay in the scanner. §4's "C2 SCOPE"
block records the ruling and why both offered options were wrong. C2.0 landed on it.

### What is already true and must not be re-litigated

* ⭐ `L` is proven and **is** the tool's one definition of fragment length (C0).
* ⭐ The unconditional anchor **exists**, in the accumulator's own frame (C1) — it closed the support-ceiling
  mismatch entirely and **55 %** of the sd gap. Nothing downstream reads it yet; that is C2.
* ✅ **`build_fl_models` anchors on `deposited_lengths`, and the scanner histogram no longer exists** (C2).
  Its only argument is the payload, so a mixed-frame call is unrepresentable rather than discouraged.
* ⛔ **The residual was TWO defects, not one — measured 2026-08-01, `docs/JUNCTION_OPPORTUNITY.md`.**
  The **+7.7 % mean** IS the junction-opportunity tilt and C3 removes it almost exactly (+8.0 % → +0.3 %
  against truth, corrected with the true θ). The **+32 % sd was NOT**: correcting perfectly closes only
  **~19 %** of it, because the rest was an **uncut intron** — present in the *anchor* as well, and so in
  every consumer including the EM's transcript effective lengths. ✅ **C2.6 fixed it**: the anchor's sd
  against truth is **+1.98 %**, from +26.97 %.
* ⛔ **BUT C2.6's D1 PUT A NEW BIAS INTO THE POOL, and it is an owner decision.** Removing mixed
  fragments from `RNA_SPLICED` removes exactly the ones whose mates sit far apart, so the pool the
  fragment-length model is FITTED FROM is now length-selected **short**: **−9.58 % mean / −22.46 % sd**
  against truth, where cutting the intron and *keeping* the fragment reads **+0.67 % / +2.40 %**. The
  purity argument for D1 still holds; what is now known is its price. `LEDGER.md` C2.6.

---

⭐ **THE ONE-LINE ANSWER TO "WHY ARE THERE TWO HISTOGRAMS?"**

> The accumulator bins **only structurally pure pools**, because purity is the whole point of §8 — so it
> has **no unconditional histogram**. The empirical-Bayes shrinkage needs one as its anchor. So the anchor
> was taken from the scanner, which measures fragment length by **two different rules, neither of which is
> the accumulator's `L`**, over **a different population**.

That is not a design. It is a gap that was filled with the nearest available array.

---

## §1 WHAT EXISTED BEFORE C2 — ⛔ HISTORY, NOT CURRENT STATE

⛔ **Everything in §1–§3 describes the tree BEFORE C2 (2026-08-01) and is kept as the record of what was
wrong and why.** Definitions **A** and **B** no longer exist; `FragmentLengthModels` no longer exists;
`build_fl_models` takes the payload and nothing else. Read `LEDGER.md`'s C2 entry for what replaced them.

### 1.1 Three definitions of "fragment length", all live

| | site | definition | population it is recorded over |
|---|---|---|---|
| **A** | `bam_scanner.cpp:1704` | `frag.genomic_footprint()` — the **genomic** span | intergenic · unspliced · **single-block only** · unique mapper |
| **B** | `bam_scanner.cpp:1810` | `result.get_unique_frag_length_mrna()` — **transcript-space** length | exonic · unique mapper · **unanimous across every candidate transcript**, else discarded to `n_frag_length_ambiguous` |
| **C** | `accumulator.cpp` / `_accumulator_reference.deposit` | **`L` = Σ path segment lengths** (span − introns, mate gap **included**) | every deposited fragment, binned into the 5 pure pools |

**A and B are summed into one array** and called `global_model` — the scanner's "unconditional"
histogram. ⛔ It is not unconditional and it is not one quantity: it is a genomic length for one subset
glued to a transcript-space length for a disjoint subset, with a third subset silently dropped.

⚠ The comment at site A says it feeds "the full library FL distribution". It does not.

### 1.2 The two containers

```
FragmentLengthModels  (scanner, frag_length_model.py)   global_model + 6 category_models
FLModels              (calibration, calibration/fl.py)  global_pmf, rna_pmf, gdna_pmf   <- EB-shrunk
```

`build_fl_models` EB-shrinks the **accumulator-frame** pools (`rna_fl_mass`, `gdna_fl_mass`) toward the
**scanner-frame** anchor (`global_model.counts`) with a single Dirichlet `POOL_EB_PRIOR_ESS`.

### 1.3 Who consumes what

| consumer | reads | note |
|---|---|---|
| the **scorer** (`pipeline.py:726-727`) | `FLModels.rna_pmf`, `gdna_pmf` | via `FragmentLengthModel.from_pmf` |
| **transcript effective lengths for the EM** (`pipeline.py:438`) | `rna_fl.compute_all_transcript_eff_lens` | ⚠ the same shrunk pmf |
| **calibration divisors** (`node_geometry`) | `gdna_fl_pmf`, `rna_fl_pmf` | `contained_eff_length` / `crossing_eff_length` |
| **the length likelihood** (P2, off) | both pmfs | ⭐ the first consumer using them as a **discriminant** |
| the QC report (`cli.py:370`) | scanner `category_models` | display only |
| ⛔ **nobody** | `ScanCache.fl_rna_counts` | written, read, and never used |

---

## §2 THE MEASURED CONSEQUENCE

On a **zero-gDNA** pilot condition, every fragment is RNA — so `global` and `rna_fl_pmf` describe **one
population**. They do not agree:

| | global (the EB anchor) | RNA_SPLICED pool | gap |
|---|---|---|---|
| mean | 210.1 | 234.5 | **+11.6 %** |
| sd | 85.5 | 146.2 | **+71.1 %** |
| support | [50, **713**] | [50, **1000**] | ⛔ different ceilings |
| mass > 500 bp | 0.10 % | **5.53 %** | **53×** |

Two mechanisms, and they are separable:

* the **junction-opportunity tilt** — longer fragments cross more junctions, so they are over-represented
  in `RNA_SPLICED`. The `rna/global` ratio rises with length (log-ratio vs length, corr **+0.70**), and the
  measured +11.6 % / +71 % sits against `S5_DESIGN_LOG.md` §3.6's independently-predicted **+14 % / +50 %**;
* the **frame mismatch** — a 53× excess of >500 bp molecules and a support ceiling 287 bp apart is not a
  smooth opportunity effect. Definition **B** is transcript-space and unanimity-gated; **C** is genomic
  path length over everything that deposited.

---

## §3 THE DEFECTS

| | defect | evidence |
|---|---|---|
| **D1** | **Three definitions of one quantity**, and two of them are summed into a single array | §1.1 |
| **D2** | **EB shrinkage mixes frames** — accumulator-frame pools shrunk toward a scanner-frame anchor. `CARRY_FORWARD.md` §3 trap 27, and this is the *shipped* instance of it | `fl.build_fl_models` |
| **D3** | **`global_model` is not unconditional**, so it is not a library FL distribution — yet it is the EB anchor *and* the QC number reported as the library's fragment length | §1.1, and the comment at `bam_scanner.cpp:1699` claims otherwise |
| **D4** | **The pools are opportunity-tilted and uncorrected.** `ACCUMULATOR_DESIGN.md` §8.1(b) names the fix and it has been **"not yet decided"** since S5.b | §2 |
| **D5** | **`ScanCache.fl_rna_counts` is dead** — written from the scanner's `SPLICED_ANNOT` category, read back, and never consumed (`calibration_inputs` uses the accumulator pool instead). A field whose name says it is live | `grep fl_rna` |
| **D6** | **The scanner silently drops fragments** whose candidate transcripts disagree on length (`n_frag_length_ambiguous`), a population the accumulator handles by a different rule (`path_ambiguous`). Two subsystems, two exclusion policies, one quantity | `resolve_context.h:137-148` |
| **D7** | **The shrunk pmf is reused for the EM's transcript effective lengths**, so a frame error in the anchor propagates into every transcript's effective length — not only into calibration | `pipeline.py:438` |

⚠ **D1–D3 are why P2 fails.** D4 is why it would still be biased after they are fixed. D5–D7 are the
blast radius.

---

## §4 THE CLEANUP

The principle, and everything below follows from it:

> **ONE definition of fragment length — the accumulator's `L` — measured in ONE place, over ONE stated
> population, with every conditioning made explicit and corrected for its own opportunity.**

`L` is the right definition and it is already the specification: `ACCUMULATOR_DESIGN.md` §3.1, executable
in `tests/native/_accumulator_reference.py`. Nothing else needs inventing.

### C0 — ✅ **DONE 2026-07-31: prove `L` before promoting it.** `tests/native/test_fragment_length_proof.py`

⭐ **Owner precondition, and it was right.** ``L`` cannot become the tool's one definition on the strength
of six hand-picked cases written by the same author as the code. It now has an independent **set-arithmetic
oracle** (no sorting, no merging, no ``searchsorted``) against which **333,684 exhaustive configurations**
— the complete space over a 12 bp reference including a 1 bp node — plus 4,000 randomised at realistic
scale all pass, for ``L`` **and** all four deposit populations. `ACCUMULATOR_DESIGN.md` §3.1's requirement
that "whatever counts toward ``L`` must also count as coverage for crossing" is now enforced from one set.
**7 perturbations, 6 caught**; the 7th is provably a QC-counter-only no-op. See `LEDGER.md`.

### C1 — ✅ **DONE 2026-07-31.** Give the accumulator an unconditional length histogram ⭐ **the keystone**

⭐ **Landed as `deposited_lengths`, and it separated §2's two mechanisms.** Against the RNA pool on a
zero-gDNA library the mean gap closed **11.6 % → 7.7 %**, the sd gap **71.1 % → 32.0 %**, and the support
ceiling mismatch (713 vs 1000) entirely. The residual is the junction-opportunity tilt — C3's target, now
against a same-frame reference. 6 perturbations; the 6th found a missing test rather than a bug. Every
existing scan cache was correctly refused and rebuilt (56 s). See `LEDGER.md`.

⛔ **Nothing downstream has moved**: `build_fl_models` still reads the scanner anchor. **C2 is the switch.**

**What.** One more row in `pool_lengths`: every fragment that deposits, binned by its own `L`, with no
purity condition. The 5 pure pools stay exactly as they are.

**Why it is first.** It is the *only* thing that removes the reason the scanner histogram exists. With it,
the EB anchor and the pools are the same measurement of the same quantity over a stated population, and
D1/D2/D3 close together.

**Cost.** One `uint32` row of `max_len` — ~4 KB. ⚠ **Reopens the S3 byte-identity gate**: the reference,
the C++ and the payload schema all move together, and `scan_cache`'s `payload_schema_digest` invalidates
every cache (by design — that digest exists for exactly this).

⚠ **`_pool()` must not change.** It returns `None` for a mixture *on purpose* — an impure pool is worse
than a missing one (§8). The new row is a separate, additional tally, not a sixth pool.

### C2 — Delete the scanner's FL histogram, or demote it to a QC diagnostic

Once C1 lands, `global_model` has no consumer that needs it. Two options, and this is an owner call:

| | option | consequence |
|---|---|---|
| **a** | **Delete** sites A and B and `FragmentLengthModels` entirely | "converge and delete". ⚠ loses the per-`SpliceType` QC breakdown in `cli.py:370` |
| **b** | **Keep as QC only**, renamed to say what it is (`scanner_length_diagnostic`), never an input to a model | keeps the breakdown; ⚠ leaves two histograms in the tree, which is the thing being fixed |

⭐ **Recommendation: (a).** The category breakdown can be rebuilt from the accumulator's own pools if it
is wanted. Two histograms is the defect; keeping one "for QC" is how it comes back.

⚠ Whichever is chosen, **site A and site B must stop being summed into one array** immediately — that is
the indefensible part, independent of everything else.

### ⭐ C2 SCOPE — measured, scoped, and ONE OWNER DECISION OPEN (2026-07-31)

**Owner ruling:** one production code path; the accumulator becomes the gold standard for fragment-length
distributions and serves `rigel report` as well.

#### ✅ G6, measured before C2 destroys it — and it confirms the audit's ordering

| | `gdna100 ss0.50 capture_off` |
|---|---|
| scanner FL observations recorded (`n_frag_length_unambiguous`) | 9,541,061 |
| scanner FL observations **dropped** (`n_frag_length_ambiguous`) | **456,613 = 4.6 %** of its own candidates |
| accumulator deposited (C1's population) | 9,634,502 |

⭐ **The scanner's anchor covered 99.0 % of the accumulator's population.** So §2's gap is almost entirely
**D1, the definition difference**, and only marginally **D6, the population difference** — the audit's
ordering (D1 first) holds, and C2's rationale is confirmed rather than assumed. ⚠ These counters are
deleted by C2 and this is the only record of them.

#### ✅ THE DECISION IS MADE, and it was NEITHER of the two options offered — owner, 2026-08-01

Both options assumed the counts had to come from, or be shaped by, the **accumulator**. The owner
rejected the premise:

> *"There are some QC counts that the scanner must be responsible for. They truly live in the scanner
> … If the QC counts are generated in one place and have no algorithmic use, there's no need to pass
> them. … I don't see a reason to propagate artifacts into the accumulator, it's the scanner's job to
> identify and filter these out."*

⭐ **THE PRINCIPLE THE AUDIT WAS MISSING.** §4's "ONE definition, ONE place, ONE population" governs
**model inputs**. The five splice counts are not a model input — they are a **census of the scanner's
own classification decisions**, with no consumer but the report. They became entangled with the
histogram only because they were *derived from it* (`category_models[stype].n_observations`), which is
an implementation accident. So C2.0 is not "move them to the accumulator"; it is **sever them from the
histogram and leave them where they are generated.**

⛔ **And `option (a)` would have been wrong on its own terms**: `SPLICE_ARTIFACT` fragments never reach
the accumulator (`bam_scanner.cpp`'s deposit adapter returns early — their span derives from an
alignment the blacklist rejected), so the stated gate "the five counts sum to `qc.deposited`" **cannot
hold**. That was found by pricing the option, not by building it.

⚠ **Under EVERY option the five VALUES move**, including (b). Today they count only the fragments that
also yielded a **length observation** — the transcript-space unanimity gate plus the single-block rule
on the intergenic path, a population never stated anywhere. The census counts every fragment the
scanner offers the accumulator. G6 measured that difference at **4.6 %**. The keys are preserved; the
population becomes *stated*, which is the point.

#### The C2 work list, in order

| | step | files |
|---|---|---|
| **C2.0** | ✅ **DONE 2026-08-01** — the splice census as **scanner QC**. One `std::array<int64_t, NUM_SPLICE_TYPES>` in `BamScanStats` + `n_deposit_not_offered`. ⭐ **The accumulator is untouched**: no deposit-signature change, no payload field, `payload_schema_digest` **stays `b7d29676c58b2c65`**, no cache rebuild, S3 not reopened | `bam_scanner.cpp`, `constants.h`, `splice.py`, `stats.py`, `pipeline.py` |
| **C2.1** | `build_fl_models(global_counts=...)` reads **`payload.deposited_lengths`** | `calibration/fl.py`, `pipeline.py:869`, `scan_cache.py:410` |
| **C2.2** | Delete the scanner's FL histogram: **sites A and B**, `FragLenObservations`, `frag_length_observations`, `_replay_fraglen_observations`, `FragmentLengthModels` (the plural container), `n_frag_length_{un,}ambiguous` | `bam_scanner.cpp`, `pipeline.py`, `frag_length_model.py`, `stats.py` |
| **C2.3** | Re-point the report: FL histograms from `deposited_lengths` + the 5 pools; splice breakdown from C2.0's counters | `cli.py:355-390` |
| **C2.4** | Delete `ScanCache.fl_global_counts` **and** `fl_rna_counts` (D5) | `scan_cache.py` |
| **C2.5** | Verify **D7**: `compute_all_transcript_eff_lens` reads the corrected RNA pmf. **Check, do not assume** — that path reaches every transcript's effective length in the EM | `pipeline.py:438` |

#### The gate for each C2 sub-step — so a half-finished C2 is visible, not silent

| | gate, written FIRST and verified failing |
|---|---|
| **C2.0** | ✅ **MET.** The gate as originally written was unsatisfiable (artifacts never deposit). The replacement is an exact identity spanning **both** subsystems, in C1's G2 form: `Σ census − census[SPLICE_ARTIFACT] == qc.deposited + Σ qc.dropped_* + n_deposit_not_offered`. ⭐ Plus an **independent** derivation of the artifact count — scan one BAM against the same index with and without a blacklist; `qc.deposited`, from a subsystem that has never heard of an artifact, must fall by exactly the artifact census. **6 perturbations, all caught**; two needed new fixtures built for them (see `LEDGER.md`) |
| **C2.1** | ⭐ **the anchor really moved** — assert `build_fl_models`' `global_pmf` equals `payload.deposited_lengths` normalised, NOT `frag_length_models.global_model.counts`. Then re-measure §2's table: the mean/sd gap to the RNA pool must read **+7.7 % / +32 %**, the C1 numbers, not the old +11.6 % / +71 % |
| **C2.2** | ⛔ **grep is the gate**: `grep -rn "FragmentLengthModels\|frag_length_observations\|_replay_fraglen\|n_frag_length_" src/` returns **nothing**. A partial delete that still compiles is the failure mode |
| **C2.3** | `rigel report` on a pilot condition produces the **same splice-breakdown keys** as before (option (a)); the FL histogram section is sourced from `deposited_lengths` + the 5 pools. ⚠ Diff the JSON against a report generated at `d045d820` |
| **C2.4** | `grep -rn "fl_global_counts\|fl_rna_counts" src/` returns nothing; a cache written before C2 is refused (the schema digest moves again) and the 8 pilot caches rebuild |
| **C2.5** | a test that `compute_all_transcript_eff_lens` is fed the post-C2 RNA pmf. ⚠ This is the one with the largest blast radius and the least visibility — it reaches every transcript's effective length in the EM |

⚠ **Expect the goldens to move at C2 and again at C3.** That is not a regression; it is the FL models
changing, which is the point. **Do not regenerate until after C3.**

⚠ **`FragmentLengthModel` (singular) STAYS.** It is the scoring/eff-length model built by `from_pmf`;
only `FragmentLengthModels` (the plural scanner container) is deleted. Conflating the two would delete the
scorer.

⚠ **Tests that will move**: `test_frag_length_model`, `test_gdna_frag_length`, `test_summary_report`,
`test_pipeline_routing`, `test_scan_cache`, `test_estimator` — 13 files hold the 48 touch points.

⚠ **Expect the goldens to move again** (the FL models change ⇒ scoring and eff-lengths change). They are
already pending regeneration from P1/P2; regenerate **once**, after C3.

### C3 — Land `ACCUMULATOR_DESIGN.md` §8.1(b): divide each pool by its own opportunity

`placements(w)` for a crossing pool, `(ell − w + 1)+` for a contained pool, **the transcript-level count**
for the junction pool. §8.1(b) already states the formulas and flags the junction one as the hard case.

⚠ **Do this AFTER C1**, not before: until the anchor and the pools share a frame there is no way to tell a
successful de-tilt from a frame artefact — the two effects are confounded in exactly the numbers of §2.

### C4 — Gate the length discriminant on both pools having data

The `strand_evidence` analogue. Its derived noise floor is `¼(1/N_rna + ω_r) + ¼(1/N_gdna + ω_g)`, so
`N_gdna = 0 ⇒ disc = 0` — the strand channel refuses to speak when there is no gDNA to calibrate against.
A discriminant between two fitted distributions must be gated by the sampling uncertainty of their
*separation*. ⚠ This is a derivation, not a flag, and it is independent of C1–C3.

### C5 — Delete `ScanCache.fl_rna_counts` (D5), and re-point D7

The dead field goes. And once C1+C3 land, `compute_all_transcript_eff_lens` should read the corrected RNA
pmf — check that it does rather than assuming it, because that path reaches every transcript's effective
length in the EM, not just calibration.

### Order

```
C1  the accumulator's unconditional histogram      <- keystone; reopens S3
C2  delete / demote the scanner histogram          <- possible only after C1
C3  §8.1(b) de-tilt                                <- measurable only after C1
C4  the discriminant gate                          <- independent, do any time
C5  delete the dead field; re-check the EM path    <- cleanup
then: re-run scripts/design/length_likelihood_ab.py
```

---

## §5 THE GATES

⭐ **The falsification test this area has never had**, and it needs no tuning:

> **On a zero-gDNA library every fragment is RNA.** So after C1 + C3, the de-tilted `rna_fl_pmf` must
> equal the unconditional histogram — and the length likelihood must go **exactly inert**
> (`length_likelihood`'s null is already proven byte-identical end to end).

Four such conditions are on disk. Pass/fail.

| gate | how |
|---|---|
| **G1** | C1's new row is byte-identical between `_accumulator_reference.py` and the C++, and bit-identical at 1/2/4/8 workers — the standing S3 gate |
| **G2** | `Σ unconditional_length_histogram == qc.deposited` — the same externally-checkable form as `Σ node_start_count`, and it fires if any conditioning creeps back in |
| **G3** | On a zero-gDNA condition, `unconditional_pmf` vs de-tilted `rna_fl_pmf`: mean and sd agree within sampling error, and the support ceilings match |
| **G4** | The length likelihood is exactly inert on all four zero-gDNA conditions |
| **G5** | The gdna100 conditions' `f_gdna` does not degrade against P1's baseline — ⛔ score both arms, trap 19 |
| **G6** | ⚠ **Re-record the scanner-drop rate first**: `n_frag_length_unambiguous` / `n_frag_length_ambiguous` are already emitted and have never been reported. If the drop is large, D6 is a bigger population difference than D1's definition difference, and the ordering above changes |

⛔ **Do not damp the length channel to make §2's numbers agree.** It is reporting a real disagreement
between two length models. Damping hides an upstream defect behind a tuned constant —
`CARRY_FORWARD.md` §3 trap 12, recorded three times over.
