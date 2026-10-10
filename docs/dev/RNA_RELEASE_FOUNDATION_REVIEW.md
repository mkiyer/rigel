# Release-foundation review: the count candidate, its cleanups, and the shortest path to 0.8.0

*2026-10-09. Independent review of the isolated count worktree and its frozen receipts. Sandbox
only: nothing here is authoritative. No source, test or golden was modified; no pipeline, suite
or benchmark was re-run. Findings rest on the patches, the archived receipts and frozen
arrays, read directly.*

## Verdict

**Integrate the count foundation now; it is not yet a release candidate.**

The three count-model changes (intron factory removed, intron strand claims restored, own strand
claims as the conditional-binomial observation law) and the two mechanical cleanups are correct
as code, genuinely simpler, and the cleanups are verified as numerical no-ops on everything they
were checked on. The weak-strand golden failure is a real, in-scope mechanism, but it is a
consequence of two standing rulings once the second prior is gone, not a defect of the
candidate's code, and every repair available is either a recorded refusal or the unfinished
population work. It does not require a solver change before integration. It does require the
owner to accept it knowingly (decision 1 below).

Three things stand between the integrated foundation and 0.8.0, in order of weight:

1. **The detector-free capture reader.** The owner's release clarification of 2026-10-08 keeps
   it a requirement; the shipped reader is a library-level gate (`landscape.located_enriched_mode`
   returning `None`), and every prototype so far has failed an in-scope stratum or the zero-DNA
   control. This is the critical path and the only item that needs design, not engineering.
2. **The landing itself**: eight goldens re-recorded deliberately, the factory deleted from the
   permanent docs under the move rule, the suite re-derived, and one identity run the cleanups
   never had (a real library).
3. **A validation hole**: all four real libraries are stranded (fitted sense fraction 0.0001 to
   0.003), so no real library exercises the path the weak-strand case tests, and no real-library
   identity receipt exists for the cleaned code.

The proposed ordering (foundation, then reader, then bounded validation) is right with two
amendments: commit the foundation before the reader is built, so the reader is measured against
a committed tree rather than a worktree reached through import hooks, and run the two cheap
validation items above before the reader rather than after it.

## What this review verified itself, and what it took from receipts

Verified by direct inspection:

- The worktree `/private/tmp/rigel-count-clean-20261008` is byte-identical to the frozen
  `after/` snapshot in `.cache/rigel_runs/2026-10-08_release_foundation/` (source and tests),
  and the frozen package's `_solve_impl` hash matches `package.json` and `verification.json`.
- The complete candidate diff against the owner's working tree (which includes the owner's
  uncommitted component-opportunity map repairs): twelve source files, 132 lines added and 498
  removed; tests 277 added and 516 removed; no script changed; the goldens are the owner's.
- The two cleanup patches, the native claim builder (`transfer_kernel.h` `strand_profile`,
  `strand_mode`, `claims`), the suite receipt (2,756 collected, 2,748 passed, 8 failed, 0 skipped,
  every failure a golden), the golden-disposition numbers, the stress-count receipts, and the
  frozen per-slot message rows of the extreme fixture (read below).
- The real-library receipts' fitted sense fractions: VCaP 0.00013, lbx0190 0.0024, mo3021 0.0023,
  lbx0588 0.0030. All four are strongly stranded antisense protocols.

Taken from receipts without re-running: byte identity of all 29 calibration fields and both
belief arrays on seven cached conditions (test zero, test stranded ON, test unstranded ON,
ladder g98 OFF, RNA-long stranded ON, RNA-long unstranded ON, RNA-short stranded OFF); the
end-to-end identity on the test chromosome; the compiled mutation coverage (five claim defects,
three interface restorations, all caught); full preflight; ruff; the 54 paired panel conditions
and four real-library pairs.

Not verified by anyone: real-library identity of the cleaned package against the validated
count candidate (`rename_identity.py --bam` was not run for either cleanup); that the
weak-strand mechanism below is the one acting on the ladder's unstranded OFF rows (inferred
from the fixture and the class errors, not traced there).

## The weak-strand case, read from the frozen rows

The fixture (`combo_extreme`): 1,000 fragments, read-name truth 54 mRNA / 246 nascent RNA / 700
gDNA; 36 spliced observations; fitted sense fraction 0.447, strand discriminability 0, so the
strand channel is dead and the library is effectively unstranded capture-OFF (in scope); 13
regions (3 intergenic, 4 intronic, 6 exonic) and 12 boundaries; nascent share 25 %, above the
ladder's own stress level. The per-class receipts:

| `combo_extreme` | current (factory) | candidate | candidate, original landscape held | candidate, EDGE row withdrawn |
|---|---:|---:|---:|---:|
| gDNA incl. intergenic (truth 700) | 677.5 | 765.2 | 714.7 | not scored |
| intron gDNA error, % of observed | 10.2 | 61.4 | 27.7 | 6.6 |
| exon gDNA error | 11.1 | 45.5 | 17.0 | 7.2 |
| boundary gDNA error | 12.8 | 16.3 | 5.6 | 6.8 |
| transcript / gene error | 12.1 / 8.0 | 12.4 / 8.3 | — | — |

Sources: `2026-10-08_count_readiness/current_stress_counts.json`, `2026-10-08_edge_ablation/stress_*.json`.

**The mechanism, from `2026-10-08_stress_evidence/combo_extreme_{removed,kept}_fit0.npz`** (the
inputs to the first landscape fit, before any population feedback):

1. With the strand channel dead, no object has own composition evidence (`tau_lam = 0`), so no
   own claim is built. The only composition information entering a gene is the gene-edge floor:
   the intergenic|exon edge's structurally pure gDNA count, delivered as "at least this gDNA
   density, nothing above" and then forwarded and transported through the gene's faces.
2. In the factory-removed arm, every delivered composition row at every genic slot of both
   genes is one-sided: it sits at its maximum at the all-gDNA end of the grid (the row's value at
   the last cell is 0.00 relative to its maximum at all 18 genic slots) and falls by 240 to 395
   nats toward the all-RNA end. Its lower shoulder (−2 nats) sits at a gDNA share of 0.10 to
   0.35. In the factory arm the same rows are two-sided in the first gene (right end −47 to −5 nats) with modes at
   0.10 to 0.27, because the intron's factory row was two-sided and the faces
   forwarded that.
3. ψ's reference is symmetric (the Jeffreys arms, median ½, kept by the refuted-location
   ruling in `DESIGN.md` §6b.1). The posterior at a genic slot is therefore the reference
   truncated below by the floor: a Beta(½,½) truncated at 0.17 has median 0.71; the frozen
   gDNA shares read 0.48 to 0.84. At the worst intron (slot 4: 38 observed, 5 true gDNA) the
   truth is 0.13 and the estimate 0.80. The edge densities themselves are unremarkable (0.068 to
   0.090 against intergenic 0.050 to 0.060): the error is not an upward fluctuation at the edge,
   it is the floor's one-sidedness meeting a symmetric reference where nothing else speaks.
4. The posterior variance under the truncated reference is 0.07 to 0.42 nat², under the
   one-nat² location floor, so all 13 regions are admitted to landscape training (`admitted`
   all true, `located` 12 of 13) at reliability weights 0.33 to 0.82. Ten genic objects at a
   gDNA share near 0.7 outweigh three intergenic anchors: the fitted population's mode moves
   from the measured background 0.057 to 0.30 and the refits return it (the frozen-landscape
   column above is what the local loss alone costs).

So the candidate's answer at an unstranded, nascent-rich object with no spliced witness is the
reference's centre under a floor, and the population then learns that centre. This is exactly
the configuration `DESIGN.md` §7.1 rule 4 names as unlocated ("a one-sided delivered row, and
its median is the reference measure's under its bound"); the implemented predicate reads
posterior variance, which a floor with a high shoulder satisfies. The factory hid this
everywhere by supplying a two-sided intron row that was a prior wearing a likelihood's clothes.

**Relevance.** The failure needs three things at once: a dead strand channel, objects whose true
composition is far from the floor (nascent-rich introns and exons without certified flux), and
few structurally pure anchors beside them. Real plasma libraries in the owner's archive are all
stranded, and the nascent scope ruling makes nascent-rich objects rare by assumption; so for the
owner's target data this is a stress reading. For an unstranded total-RNA library with heavy
intronic signal it is the typical case, and none is in the validation set. Low depth matters
through the anchor count, not the per-object count: the 1,000-fragment golden panel already
reads this regime, and 2 of 21 toys move materially (the two combined-stress cases), 6 by less
than one gDNA fragment, 13 not at all. On the deep panels the in-scope cost is small and visible
in calibration before it is in the transcript table:

| unstranded × capture-OFF, in scope | current | candidate |
|---|---:|---:|
| ladder, contaminated rows: transcripts / genes % | 2.29 / 0.31 | 2.30 / 0.31 |
| ladder g98 ss0.50 OFF, the worst row | 16.87 / 5.65 | 17.51 / 5.96 |
| ladder g98 ss0.50 OFF, region / boundary gDNA error % | 0.79 / 5.71 | 2.46 / 9.12 |
| test chromosome, contaminated rows | 8.66 / 1.33 | 8.81 / 1.34 |
| test g50 ss0.50 OFF, region / boundary gDNA error % | 1.30 / 2.42 | 1.38 / 3.86 |
| RNA-short gap, g50 | 3.79 / 0.27 | 3.81 / 0.27 |
| RNA-long gap, g50 | 2.16 / 0.22 | 2.14 / 0.23 |

The EDGE-withdrawal ablation (`ISSUES: the-edge-density-floor-under-capture`) is the cleanest
reading of what the floor buys and costs now that the factory is gone: withdrawing it repairs
the stress fixture almost completely (table above) and does not lose any measured in-scope row
(ladder g98 OFF 2.47 / 9.12 → 2.25 / 8.83; stranded rows and the zero control unchanged), while
the deferred unstranded capture-ON rows collapse (test g50 ss0.50 ON 3.69 / 2.54 → 43.5 / 61.9).
The floor is the deferred stratum's only witness of enrichment on unstranded data. That is the
tradeoff the owner is being asked to hold, and it is correctly not resolved by a constant.

## Findings, ranked by release impact

### 1. The capture reader is the critical path, and it has no passing candidate

**Evidence.** `DESIGN.md` §0b (release priority, 2026-10-08): "the existing detector-free capture
and fragment-length requirements remain." The shipped reader reads the landscape's located
enriched mode and returns `None` where none exists (`landscape.located_enriched_mode`,
`capture_efficiency.capture_efficiencies`): a library-level detector, which the 2026-10-05
spectrum ruling forbids. Three prototypes failed on the test chromosome (median ruler,
geometric ruler, geometric with count uncertainty: `ISSUES:
calibration-detects-capture-on-a-capture-off-library`, the ruler screen), and the point-count
reader turns 18.6 inferred gDNA fragments at a zero-truth boundary into a capture weight of
23,831 (`RNA_SHORT_READER_REVIEW.md`). The count candidate does not change this: every
panel and real-library receipt in the plan used the old reader.

**Remedy.** Treat the reader as the one owner design checkpoint and scope it to a consumer, not a
model: the existing posterior-mean reader keeps its structure, the `None` path and the clip go,
and the term that reads the inferred count as a Poisson observation is replaced by the object's
own-observation likelihood over its gDNA density (the certified primitive in
`2026-10-07_density_checkpoint/density_evidence.py`, 113 gates). Its decisive test is already
written: the zero-DNA controls on all three probe layouts, the four capture-OFF gap rows, and the
three in-scope strata at roughly the current numbers, read per stratum with the deferred rows
reported. Reject it if any capture-OFF row's transcript error moves by more than the current
tree's own A/B noise, or if any zero control's false gDNA rises. Do not let it carry a new prior,
a threshold or a class branch.

### 2. The weak-strand failure is an accepted-tradeoff decision, not a pre-integration fix

**Evidence.** The mechanism section above; the frozen-landscape and EDGE-withdrawal controls;
the per-stratum in-scope table. Every available correction is either a recorded refusal
(restore the factory; blanket EDGE withdrawal; "every reading of a bound that reaches the
delivered rows", refused in `ISSUES: the-landscape-training-population-arms`) or the
uncertainty-preserving population work, which is measured and trades class errors
(`FROZEN STRESS-INPUT COMPARISON` under `ISSUES: the-gdna-prior-enters-psi-twice`).

**Remedy.** Integrate without a solver change; re-record the eight goldens as an intended
change with the movements already recorded in `ISSUES` (the three preconditions of the
golden-update rule are met: diff read, magnitudes recorded, truth-scored instruments checked);
keep `combo_extreme` and `combo_moderate` in every future validation as the named weak-strand
stress pair. The smallest theoretically justified correction, if the owner wants one, is in the
owner-decision section (decision 2), because it re-opens a refusal.

### 3. Two cheap validation items are missing and should precede the reader

**Evidence.** No `rename_identity.py --bam` receipt exists for the direct-assembly, plumbing or
interface cleanups; the identity instrument's documented use for a refactor is `--check` after
each stage and `--bam` on a real library. The real-library validation of the count candidate is
sound but says nothing about the weak-strand path: every library is stranded.

**Remedy.** (a) Freeze the frozen count candidate's package as the reference and run
`rename_identity.py --check --bam` on lbx0190 (16 s per pipeline) with the cleaned package;
expect bit identity; any difference stops the landing. (b) If the owner's archive holds an
unstranded library, run the one paired pipeline on it and read the gDNA fraction and the
intronic-pool change beside the stranded four; if none exists, record the gap in `ISSUES` as a
known limit of 0.8.0's validation rather than inferring from the stranded libraries.

### 4. The eight golden failures: six negligible, one small stranded cost, one mechanism

| fixture | protocol / content | gDNA pool move (fragments) | max transcript move | reading |
|---|---|---:|---:|---|
| `antisense_overlap_ss90` | stranded, no gDNA | −0.0003 | 0.0002 | intended, nil |
| `antisense_contained_ss90` | stranded, no gDNA | +0.99 | 0.73 | intended; the fixture already reads 102 false gDNA fragments of 1,000 in both arms (the tiny-toy limit recorded in `DESIGN.md` §7.1) |
| `gdna_light` | stranded, gDNA only | +0.03 | 0.02 | intended, nil |
| `gdna_heavy` | stranded, gDNA only | +0.41 | 0.31 | intended; transcript error 13.4 → 14.6 % on a 35-fragment mRNA pool |
| `nrna_moderate_ss90` | stranded, nascent | +0.75 | 0.01 | intended, nil |
| `nrna_heavy_ss90` | stranded, nascent | +0.30 | 0.003 | intended, nil |
| `combo_moderate` | stranded (discriminability 0.67), gDNA + nascent | +4.94 (truth 687; 690.6 → 695.5) | 1.45 | small cost: transcripts 2.06 → 3.54 %, genes 2.06 → 1.63 %; both native changes act here |
| `combo_extreme` | dead strand channel, gDNA + nascent | +87.7 (truth 700; 677.5 → 765.2) | 0.32 | the mechanism above |

The goldens are output-regression gates with a 1e-6 tolerance, not truth gates, so re-recording
them is the right disposition once the tradeoff is accepted; the truth numbers already live in
`ISSUES`. The movement must be recorded in the commit that updates them, and the new values must
come from the integrated main tree, not the worktree.

### 5. The permanent docs still describe the factory

**Evidence.** `DESIGN.md` carries 14 factory mentions (§6b.9's foundation, §7.1 rule 4's "an empty
intron's factory row", the code-layout table), `EQUATIONS.md` three (§3b's "the intron's
composition is prior-free (the intron factory)", §9c's "the intron factory's rows", §9d's
`T = 6` term count, which is now five), and `calibrate.py`'s module docstring in the worktree
already says the opposite. The worktree's docs are older than main's; main's must win.

**Remedy.** Under the move rule, in the landing commit: delete the factory from DESIGN and
EQUATIONS, retitle nothing, and add one line to `ISSUES: the-gdna-prior-enters-psi-twice`
saying the factory half is gone and the reference-half tilt remains the open half. Update
`CLAUDE.md`'s layering table only if `density_deconv` is deleted (finding 6) and re-derive the
baseline count line. `tests/test_docs_boundary.py` and the constants gate must stay green.

### 6. Remaining dead state after the cleanup (low impact; each a provable no-op)

- `psi_kernel.h` `SlotInputs.lam_prior` and `solve_kernel.cpp`'s `lam_logprior` argument are
  always null in production now; only three gates inject synthetic rows through them
  (`test_sweep.py`, `test_vertex_reference.py`). Either delete the parameter and port those
  gates to `row_logprior`, or document it as the gates' injection seam in the binding docstring.
- `density_deconv.py` (`fit_gdna_background`, `fit_intron_background`, `GdnaBackground`) has no
  production caller; its docstring says "the caller supplies a pure-DNA pool" and there is no
  caller. It is kept alive by its remaining tests and one of the three removal gates. Decide at
  the reader checkpoint: if the reader's off-target floor consumes this measurement, it earns
  its place; otherwise delete it with its tests and drop it from `_layers.py`.
- `_Strand.gdna_strand_overdispersion` and `rna_strand_overdispersion` are 0.0 by ruling
  (`DESIGN.md` §3.3a) and are threaded through `solve_chain`, `solve_blocks` and the ψ
  bindings, and reach only the ψ grid, where they are zero by ruling. The adaptive
  od is deferred past 0.8.0 (`ISSUES: strand-overdispersion-one-shared-value`). Removing both
  is the same shape of cleanup as the one just landed and should be proven the same way.
- `solve_blocks` receives one number twice: `kappa` (the solver's) and `policy_kappa` (the
  policy's), always equal in production. `has_strand` is the only bit the policy adds.
- `transfer_kernel.h` `claims()` tests `is_exon || is_bnd || intron` after `single`; `single`
  already excludes intergenic, so the class list implies a class rule that does not exist.
  Reduce to `single && has_own`; the existing missing-information gate covers it.
- The `messages/__init__.py` docstring says builders read no posterior belief. The rows do not;
  the `has_own` mask is derived from the self-solve, which the incoming belief seeds. True
  except at an exact vertex; say so in one clause.

### 7. One falsification the observation-only claim still lacks

**Evidence.** The belief-invariance gate was replaced by an interface-rejection gate (the binding
refuses `belief_fg`, `od_g`, `od_r`). That proves the arguments are gone, not that the delivered
rows are belief-free through `solve_blocks`, which still receives the beliefs for the count
solver.

**Remedy.** One gate: run `solve_chain` under a diagnostics capture twice with `belief_fg`
perturbed across its range and assert `lam_rows` (the delivered rows) identical at pass zero.
The compiled `restore_gaussian` mutant already exists and must make it fail.

### 8. The A/B machinery and the sandbox are becoming the maintenance burden

**Evidence.** Thirty run directories under `.cache/rigel_runs/` dated 2026-10-08; candidates
selected by `sitecustomize.py` import hooks, process-local constructor replacements and
environment variables (`RIGEL_DENSITY_SITE`); 38 files in `docs/dev/`, several of which are
checkpoints of checkpoints; a "Goal service" paragraph in `RNA_SHORT_FIX_PLAN.md` that is
tool state, not a plan.

**Remedy.** After the foundation commits, every further A/B is worktree-build against main per
the documented recipe, and the hook harnesses are retired with the directories they served.
Collapse the research checkpoints into the plan plus one ledger (the settled numbers are
already in `ISSUES`), and delete the Goal-service paragraph now.

### 9. The deferred-stratum regression is reported, not fixed, and should stay that way for now

**Evidence.** RNA-long g50 unstranded capture-ON: 6.55 / 0.94 → 43.22 / 11.26 % (0.7.1: 28.99 /
5.84); population feedback owns most of it (41.7 / 32.9 → 11.7 / 12.4 with the original
landscape held). The ladder's deferred rows improve (g98 ss0.50 ON 137 → 88 %). The owner
authorized investigating this for robustness while keeping it deferred.

**Remedy.** No action before the reader. The investigation already has its instrument (the
fixed-landscape control) and its first result; it belongs after the reader because the reader
changes what the population feeds.

## The six questions, briefly

1. **Is the cleanup correct, simpler, verified?** Yes on all three for the two mechanical
   cleanups: the patches remove exactly the state they claim (`_Strand.model`, the policy's
   dispersion pair, the native `belief` and `src` pointers, the factory kernels and `tau_fac`),
   the receipts show byte identity on seven conditions and end to end, and the mutation
   coverage is real (compiled defects, assertion failures). Remaining state is finding 6; none
   of it changes a number.
2. **The eight goldens:** finding 4. Six intended and negligible, one small stranded cost, one
   in-scope mechanism that is unresolved evidence about the model, not a defect in the code.
3. **Does the weak-strand failure require correction before integration?** No. Mechanism and
   relevance above; finding 2.
4. **If correction is wanted, the smallest justified change:** decision 2 below, with its
   falsification and rejection conditions; it is an admission predicate, carries no constant,
   and re-opens a refusal, which is why it is the owner's.
5. **What blocks readiness:** the reader (finding 1), the landing (findings 4, 5, 3a), the
   unstranded validation gap (3b). Documented research, not blocking: density evidence, evidence
   admission, equal weights, route profiling, both-strand consumers, EDGE withdrawal.
6. **Ordering:** appropriate with the two amendments in the verdict. The capture-and-opportunity
   question the plan raises (can gDNA enrichment transfer to RNA across a length gap) is a
   limitation of the frame, not a regression of the candidate: both gap panels' stranded ON
   rows move by at most 0.03 points under the candidate. The bounded contrast that answers it
   is the existing `ruler_vs_truth.py --scale` on those four rows, probed and unprobed classes
   apart, current against candidate, from cache; the result that would require a new model is
   a within-gene spread that grows with the gap on the candidate and not on the current tree.

## The smallest next implementation checkpoint: land the foundation

Steps, all engineering, in order:

1. Take the worktree's `src/` and `tests/` onto the owner's tree (the worktree already contains
   the owner's uncommitted map repairs, so this is the worktree minus its stale `docs/`);
   rebuild with `pip install --no-build-isolation -e ".[dev]"`.
2. Identity: `rename_identity.py --check` against the frozen count candidate on the test
   chromosome, and `--bam` on lbx0190 (finding 3a). Bit identity or stop.
3. Goldens: `--update-golden` once, the eight movements named in the commit (finding 4).
4. Docs: the factory out of DESIGN and EQUATIONS; the one-line note in ISSUES; the plan's
   Goal-service paragraph deleted (findings 5, 8).
5. The belief-invariance gate (finding 7) and the `claims()` simplification (finding 6), each
   proven a no-op by step 2's instrument.
6. Re-derive the suite, ruff, `preflight.py --full`; re-derive the baseline line in `CLAUDE.md`.
7. Commit (owner).

Acceptance criteria: zero failures on the re-derived suite with the count derived, not adjusted
(2,756 collected before step 5's additions); bit identity at step 2 on both substrates; full
preflight green; the docs gates green; the ladder report regenerated with the deferred rows
present; the per-stratum numbers in `RNA_FULL_PANEL_VALIDATION.md` reproduced by the committed
tree on one condition per stratum before the receipts are retired.

Then the reader checkpoint (finding 1), which begins with an owner discussion, not code.

## Owner decisions

1. **Accept the weak-strand tradeoff and authorize the golden update.** The in-scope cost is
   bounded (tables above); the stress fixture stays in validation by name; the deferred
   regression stays reported. This is the only decision the landing waits on.
2. **Whether to re-open one refusal, after integration, not before.** The one change that is
   small, constant-free and derivable from the standing rules is to make admission honour
   `DESIGN.md` §7.1 rule 1 for transported floors: a delivered composition row that is flat at
   the all-gDNA end of the grid is a bound, and a slot whose only composition information is a
   bound does not train the landscape. Falsification: three fixture gates (an edge floor alone
   admits nothing; a two-sided source admits; a flat row admits nothing), the stress fixture's
   first landscape trained on its anchors alone (prediction: the gDNA pool moves from 765 toward
   at most 715, the frozen-landscape control's value), and the seven cached conditions at pass
   zero and full. Reject it if any in-scope stratum's region or boundary gDNA error worsens
   against the candidate, if any zero control's false gDNA rises, or if an unstranded OFF
   contaminated row loses against the silent policy. It will lose the deferred stratum, as the
   EDGE ablation did, because those one-sided rows are that stratum's only witness; and a
   reading of "bound" that reached delivered rows was refused with numbers on the factory tree.
   The refusal's measurement predates factory removal, which is the only reason to ask.
3. **The reader's scope.** Confirm that the reader checkpoint is a consumer change (finding 1's
   remedy) and that its acceptance bar is the current tree's in-scope numbers plus the zero
   controls, not an improvement.
4. **The fate of `density_deconv`.** Delete unless the reader consumes the intergenic background
   (finding 6).
5. **Retire the hook harnesses and the run directories after the commit** (finding 8).

## Receipts read

`.cache/rigel_runs/2026-10-08_release_foundation/` (`cleanup.patch`, `package.json`,
`verification.json`, `suite.json`, `claim_mutations.json`, `interface_mutations.json`,
`preflight.log`, `after/`); `2026-10-08_count_cleanup/` (`cleanup.patch`, `retired_tests.json`,
`verification.json`); `2026-10-08_count_readiness/` (`candidate.json`, `golden_disposition.json`,
`stress_truth.json`, `current_stress_counts.json`, `verification.json`);
`2026-10-08_stress_evidence/combo_extreme_{removed,kept}_fit{0,2}.npz`;
`2026-10-08_evidence_attribution/trace.json`; `2026-10-08_edge_ablation/` (`stress_*.json`,
`verification.json`); `2026-10-07_full_panel_counts/diagnosis.json`;
`2026-10-07_release_followup/real/` (`manifest.json`, `comparison.json`, per-library JSON).
Source read in the worktree: `calibrate.py`, `sweep.py`, `messages/`, `density_deconv.py`,
`landscape.py` (`_reliability`), `capture_efficiency.py`, `native/transfer_kernel.h`,
`native/solve_kernel.cpp`, `native/psi_kernel.h`, the two new gate files.
