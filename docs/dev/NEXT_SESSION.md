# NEXT SESSION — start here (2026-09-19, end of day: restore nascent RNA's share of the prior, re-measure, then the EM's gDNA split)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are `docs/ISSUES.md`,
the rulings and the record are `docs/DESIGN.md` (§0b carries the scope and its 2026-09-19 amendment), what
"done" means is `docs/SUCCESS.md`, the lessons are `docs/TRAPS.md` cited by name, and the release procedure is
`docs/MANUAL.md` / `docs/PUBLISHING.md`.

## Where the tool is

Everything through `61e82758` is landed and pushed; the tree is clean, the suite is 3,425 passed / 5 xfail /
3,430 collected and `preflight.py --full` is green. The ladder was REBUILT today under corrected capture physics
(a fragment binds through one contiguous part of a probe; the `gdna_split_penalty` is gone), and on it the
deliverable reads, under fractional assignment, transcript-level Σ|Δ| as a share of the true RNA at
`g00` / `g05` / `g50` / `g98`:

| stratum | numbers |
|---|---|
| unstranded × OFF | 3.1 / 3.4 / 4.2 / 35.0 % |
| stranded × OFF | 3.3 / 2.9 / 3.9 / 24.8 % |
| stranded × ON | 6.8 / 4.5 / 7.8 / 109.8 % |
| unstranded × ON (deferred) | 8.3 / 12.4 / 20.8 / 739.5 % |

The reseed floor is 0–82 fragments. A perfect prior recovers nothing in scope; the simulator's own capture
lengths take stranded ON to 2.6 / 3.4 / 8.6 %; true per-transcript weights halve every stratum.

## ① THE FIRST JOB — restore nascent RNA's share of the RNA prior (`ISSUES: nascent-gets-no-rna-prior`)

The owner's ruling (2026-09-19): the present rule is a HACK and nascent fairness is restored, so that the
per-transcript prior can be described as distributing the RNA pseudocounts uniformly over the RNA components —
none singled out for zero.

**What it is today.** `native/em_solver.cpp:apply_grouped_prior_update` hands the locus's RNA pseudocount only to
components the annotation asserts exist; a synthetic nascent entity is excluded (Dirichlet `alpha = 0`). The flag
arrives as `t_is_synthetic` from `estimator.run_batch_locus_em_partitioned`; `assemble_priors` supplies only the
per-locus total. The derivation is `EQUATIONS.md` §9b and §9b.1.

**The one open choice — read §9b.1 before deciding.** The shipped allocation is `a_i = P · raw[i] / Σ_eligible raw`
— in proportion to the EM's own current belief — and at `raw[i] = 0` it is ABSORBING, which is what stops a zombie
entity being revived by prior mass alone. Admitting the entities at that same weight keeps the absorbing state and
is the smaller change; an equal share per component removes it, and §9b.1 names the activation threshold to design
against (the VBEM fixed point, ~0.16–0.47 alpha units, not the exponential cutoff 0.0014). Decide it explicitly
and write the reason down.

**The invariant that may not move.** The prior is redistributed strictly WITHIN the RNA pool: `Σ out over the RNA
components = rna_count + rna_prior` per M-step, so the library gDNA fraction cannot move by this rule. A gate in
`tests/test_estimator.py` holds it — keep it, and watch it fail under a perturbation that leaks the prior across
the gDNA boundary.

**What moves with the change.** `EQUATIONS.md` §9b and §9b.1; the gates in `tests/test_estimator.py` that pin the
synthetic branch (including the bit-identity of a locus with no synthetic component); the strict xfail
`tests/scenarios/test_antisense_intronic.py::test_nrna_multiexon_t2_low_ss`, which this rule owns (72 fragments on
the annotated antisense `t2` against a limit of 50) — if the restoration closes it, replace it with a test that
asserts the new rule's promise, and if it does not, record the number in the entry. If the `is_synthetic` plumbing
ends up unused, delete it rather than leave it dead (`TRAPS: converge-and-delete`). It is a native change, so
rebuild (`pip install --no-build-isolation -e ".[dev]"`).

**Then RE-MEASURE, before anything else is built on it.** On the rebuilt ladder: `panel.py score` (it passes
`--set em.assignment_mode=fractional` itself) for `base base_reseed oracle oracle_ruler`, plus
`calibration_vs_oracle.py`, `zero_controls.py` and `policy_benchmark.py --panel ladder`. Put the numbers in
`ROADMAP.md`'s deliverable claim beside today's, so the effect of the restoration is visible as its own step.

## ② THEN the EM's gDNA split (`ISSUES: em-overturns-the-calibrated-gdna-split`)

The dominant in-scope residual and the whole of `g98`: under capture the table reads 0.4793 against 0.50 at
`g50 ss.99 ON` while calibration reads +0.9 %; at capture-OFF the same EM over-calls gDNA by taking unspliced RNA
(a nascent-stress reading, so that half is sized at the realistic share first,
`ISSUES: nascent-stress-sensitivity`). Neither the prior, the true ruler nor the gDNA component's length moves it.

The leading hypothesis, to test first: the EM gives each locus ONE gDNA rate spread uniformly along it, while
under capture gDNA's density is an order of magnitude higher at probed exons — where the RNA also sits. Inside a
probed exon the model then under-predicts gDNA and the surplus goes to RNA. It fits the direction, the
concentration in heavily probed isoform-rich genes, and its immunity to the prior and the ruler. Calibration
already publishes per-piece efficiencies, so the test is whether a position-dependent gDNA weight closes the
split. Start with a per-locus attribution on `g50 ss.99 ON` (the EM's gDNA against the certified truth, ranked by
mass), then the alternatives: the prior's weight against the likelihood, and gDNA's strand handling at ss 0.99.
The per-fragment instrument is `confusion.py` beside the dissection data (below).

Then `ROADMAP.md` items ③ (the capture ruler, which waits on an owner decision) and ④ (the per-transcript lane).

## The five xfails are proven defects, each deferred to its thread

`ISSUES: two-sided-exon-row`; `ISSUES: nascent-gets-no-rna-prior`;
`ISSUES: the-lower-bound-noise-ratchet`; `ISSUES: nested-antisense-leak-under-the-sane-ruler` (two rungs).
Closing one means repairing the thing or asserting the invariant structurally, never widening a bound.

## What gates the release (`ROADMAP.md`'s last item)

`PUBLISHING.md` is two commands and a wait; the STATE is what gates it. The deliverable measured and not
regressed per stratum; the zero controls at 0.000 and 1.000 (`zero_controls.py`, and the `g00` rung); the
suite at its standing count with an empty failure set; `preflight.py --full` green; the standing risks re-read
(`ISSUES: capture-degeneracy-standing-risk`, `ISSUES: flgap-panels-stale-nascent-model`); and the manual true
of what ships.

## Standing rulings carried (unchanged)

* Float64 for the whole of ψ; ONE solver; ONE λ lattice; no θ lattice; the tilt's hypothesis space is
  {pure +, pure −, mixed}. The strand channel's liveness is a protocol decision on the spliced 2×2.
* The landscape trains only where a solve locates a slot; the ruler reads its located enriched mode and admits
  failure where there is no gDNA to read.
* Real data is a test input, never a design input. A few high-quality instruments, kept current. The source
  cites no doc. One production path: once native code is validated the Python it replaces is deleted.
* Unstranded × capture-ON stays DEFERRED — reported on every benchmark, never a development target, never
  ranked on a pooled total. The fragment-length composition channel stays retired.
* The owner drives commits. The refit count and CI's on-demand trigger are the owner's.

## Where everything is

* The synced scratchpad: `~/Downloads/rigel_runs/prototypes/2026-09-17_port_ii/` — `commits/committed/` (every
  snapshot, all landed), `s16/`–`s18/` (the port and the block), `s19/`–`s21/` (the performance campaign, its
  attributions and the taper study that the parked thread resumes from).
* Captures and reports: `perf/sweeps_VCaP_step19` (the sweep replay's capture, bit-identical on this tree);
  `perf/plan_final_2026-09-19/` is the deep run's current before-and-after.
* THE REBUILT LADDER's results: `~/Downloads/rigel_runs/suite/ladder/arms/` (the four `quant_accuracy` arms plus
  `oracle_alloc_seed`, and `calibration_vs_oracle.json`) with every stage's log in
  `~/Downloads/rigel_runs/logs/ladder_*.log` and the chain that produced them in `ladder_rebuild_2026-09-19.sh`.
* THE BASELINE on the RETIRED ladder: `~/Downloads/rigel_runs/arms/2026-09-19_e2e_baseline/` (sampled assignment, the per-stratum
  `qa_report.txt`, `decomposition.txt`, `alloc.txt`). THE STRANDED × ON DISSECTION:
  `~/Downloads/rigel_runs/arms/2026-09-19_stranded_on/` — `all_scenarios.txt` (the fractional panel),
  `stranded_on_arms.txt` (every ruler and prior arm), `tables/` (per-transcript tables per arm), `testchr/` (the
  benign-vs-junction control), `confusion_*` (per-fragment truth against assignment at `g50`), and the scratch
  runner `dissect_run.py` / `dissect_analyze.py` / `confusion.py` it came from.
* The identity references: `~/Downloads/rigel_runs/arms/review_identity_*.json`, BIT-IDENTICAL on this tree.
  They are the gate for any change that must not move a number.
