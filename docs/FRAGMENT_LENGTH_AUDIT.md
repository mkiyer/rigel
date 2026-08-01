# Fragment length — the audit, and the cleanup

    Status: AUDIT COMPLETE, 2026-07-31. Nothing proposed here has landed.
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
      pytest tests/ -q                                  # expect 1809 / 22
      ruff check src/ tests/ scripts/                   # ⚠ never `ruff format scripts/`

⚠ **1809 / 22 is the baseline, and 22 is the CORRECT number right now**: 21 are `test_golden_output`
moving numerically after P1's units fix, plus `TODO.md` §7's owner call. **A 23rd failure is a
regression.** Do not regenerate the goldens to make the count look better -- they move again at C2/C3.

| | step | status | gate |
|---|---|---|---|
| **C0** | Prove the accumulator's `L` before promoting it | ✅ **DONE 2026-07-31** | 333,684 exhaustive configs vs an independent set oracle; 7 perturbations, 6 caught, 7th a proven no-op |
| **C1** | Give the accumulator an unconditional length histogram (`deposited_lengths`) | ✅ **DONE 2026-07-31** | `Σ = Σ node_start_count = qc.deposited`, enforced in the accumulator AND at the payload door; parity gate auto-covers it; 6 perturbations |
| **C2** | ⭐ **NEXT** — switch every consumer to the accumulator, delete the scanner histogram, re-point `rigel report` | ⛔ **needs one owner decision** (below) | the report's output unchanged under option (a); suite green apart from the known goldens |
| **C3** | §8.1(b): divide each pool by its own opportunity | blocked on C2 | ⭐ the zero-gDNA inertness test (§5) |
| **C4** | Gate the length discriminant on both pools having data | independent — any time | the `strand_evidence` analogue |
| **C5** | Delete `ScanCache.fl_rna_counts`; verify D7 | with C2 | — |
| **→** | **re-run `scripts/design/length_likelihood_ab.py`**, then decide `CalibrationConfig.length_likelihood` | after C3 | `SOLVER_OBSERVABLES_PLAN.md` §6.4 |
| **→** | **regenerate `tests/golden/` ONCE** | after C3 | ⚠ 21 goldens are stale from P1 and will move again at C2/C3 — regenerate **once**, at the end, and ⛔ regenerate twice and diff (`TODO.md` §7) |

⛔ **THE ONE THING BLOCKING C2 IS A DECISION, NOT WORK** — `rigel report`'s five per-fragment splice-type
counts have no accumulator equivalent. Options **(a)** and **(b)** are priced in §4's "C2 SCOPE" block;
the recommendation is **(a)**. Nothing else in C2 is uncertain.

### What is already true and must not be re-litigated

* ⭐ `L` is proven and **is** the tool's one definition of fragment length (C0).
* ⭐ The unconditional anchor **exists**, in the accumulator's own frame (C1) — it closed the support-ceiling
  mismatch entirely and **55 %** of the sd gap. Nothing downstream reads it yet; that is C2.
* ⛔ `build_fl_models` **still** anchors on the scanner histogram. Every FL number in the tool today is
  still the old one.
* The residual **+7.7 % mean / +32 % sd** gap on a zero-gDNA library is the junction-opportunity tilt, and
  it is C3's target — now measurable against a same-frame reference.

---

⭐ **THE ONE-LINE ANSWER TO "WHY ARE THERE TWO HISTOGRAMS?"**

> The accumulator bins **only structurally pure pools**, because purity is the whole point of §8 — so it
> has **no unconditional histogram**. The empirical-Bayes shrinkage needs one as its anchor. So the anchor
> was taken from the scanner, which measures fragment length by **two different rules, neither of which is
> the accumulator's `L`**, over **a different population**.

That is not a design. It is a gap that was filled with the nearest available array.

---

## §1 WHAT EXISTS TODAY

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

#### ⛔ THE OPEN DECISION: `rigel report`'s SPLICE BREAKDOWN has no accumulator equivalent

`cli.py:376-384` reports **five per-FRAGMENT splice-type counts** off `flm.category_models`, which C2
deletes. They are a different classification from anything the accumulator keeps:

| `rigel report` needs (per fragment) | the accumulator has today |
|---|---|
| `unspliced` | ⛔ nothing |
| `spliced_annotated` | ⛔ nothing |
| `spliced_unannotated` | ⚠ `unannotated_introns` — counts **INTRONS, not fragments** |
| `spliced_implicit` | ✅ `sj_implicit_fragments` |
| `splice_artifact` | ⚠ `contradictory_sj_strand` — related, not equal |

⚠ **The five pure pools are NOT this breakdown either** — they are a structural gDNA/RNA classification
(`DNA_INTERGENIC` … `RNA_SPLICED`), deliberately conditioned and deliberately incomplete.

**Two ways forward, and it is the owner's call:**

| | option | cost |
|---|---|---|
| **a** ⭐ | the accumulator grows **five per-fragment splice-type counters** in `DepositCounters`, mirroring `SpliceType`. The report reads them | small and additive (~15 lines each side), ⚠ reopens the S3 byte-identity gate and moves `payload_schema_digest` again. Keeps the report's output **unchanged**, which is the honest test of "converge, don't lose" |
| **b** | the report's splice section changes shape to what the accumulator natively knows (the 5 pools + the existing QC counters) | no accumulator change; ⚠ **changes `rigel report`'s output** and loses `unspliced` / `spliced_annotated` entirely |

⭐ **Recommendation: (a).** The owner's stated goal is one code path *without losing the QC*, and (b)
loses two of five categories. Doing it inside C2 also means the schema digest moves **once**, not twice.

#### The C2 work list, in order

| | step | files |
|---|---|---|
| **C2.0** | **(a) above** — five splice-type counters in the accumulator, reference first, then C++, then payload. Byte-identity gate + a QC-sum invariant | `_accumulator_reference.py`, `accumulator.{h,cpp}`, `bam_scanner.cpp`, `scan_payload.py` |
| **C2.1** | `build_fl_models(global_counts=...)` reads **`payload.deposited_lengths`** | `calibration/fl.py`, `pipeline.py:869`, `scan_cache.py:410` |
| **C2.2** | Delete the scanner's FL histogram: **sites A and B**, `FragLenObservations`, `frag_length_observations`, `_replay_fraglen_observations`, `FragmentLengthModels` (the plural container), `n_frag_length_{un,}ambiguous` | `bam_scanner.cpp`, `pipeline.py`, `frag_length_model.py`, `stats.py` |
| **C2.3** | Re-point the report: FL histograms from `deposited_lengths` + the 5 pools; splice breakdown from C2.0's counters | `cli.py:355-390` |
| **C2.4** | Delete `ScanCache.fl_global_counts` **and** `fl_rna_counts` (D5) | `scan_cache.py` |
| **C2.5** | Verify **D7**: `compute_all_transcript_eff_lens` reads the corrected RNA pmf. **Check, do not assume** — that path reaches every transcript's effective length in the EM | `pipeline.py:438` |

#### The gate for each C2 sub-step — so a half-finished C2 is visible, not silent

| | gate, written FIRST and verified failing |
|---|---|
| **C2.0** | byte-identity to the reference on the new counters (the parity gate picks them up automatically off `dataclasses.fields(Tally)`); a QC-sum invariant in the same form as C1's: the five splice-type counts sum to `qc.deposited`. ⚠ Perturb: drop one category and watch the sum fire |
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
