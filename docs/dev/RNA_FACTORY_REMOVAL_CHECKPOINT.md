# Intron-constraint removal and observational messages

*2026-10-07. Development checkpoint following the owner's authorization to remove
the intron constraint and continue investigating the single-landscape count model.
No production implementation or release readiness is claimed.*

## What was implemented

Three count-changing contrasts were implemented and measured separately. The starting
point is the owner's uncommitted count-repair tree, not bare `a7b03103`.

| Contrast | Exact change | What it isolates |
|---|---|---|
| Remove the intron factory | A scoped Python replacement for `_IntronFactory` retains the measured intergenic background but supplies `rows=None` to every sweep | Removal of the background constraint and its transported claims |
| Restore intron strand messages | In a detached native worktree, extend `claims()` to publish the existing single-strand row from an intron when no factory row exists | The observational channel omitted by removal alone |
| Make strand messages observational | Replace `Chain::strand_profile()`'s belief-frozen Gaussian with the exact conditional binomial likelihood of its two observed columns | Remove the incoming composition belief from the strand-message likelihood |

The latter two changes affect one C++ header. The wiring contrast is one condition;
the observation contrast is a short row evaluator. No new biological constant,
capture detector, expression threshold or RNA-abundance sharing rule was added.

The factory-enabled control in the wiring-only build reproduces all 30 current-tree
test-chromosome count arrays bit-for-bit. The exact-row build changes exon and boundary
messages as well as the newly enabled intron messages; it therefore has its own A/B.

These are **intermediate count experiments**. The final psi solve still uses its
Gaussian observation approximation. Its DNA reference term and landscape assembly,
the point-count landscape fit, and the existing capture reader are unchanged. A strand
row over an expected mixture fraction is not by itself a DNA-density likelihood.

## The density-evidence reference

`local_face.py` extends the preceding own-observation reference to one licensed,
single-strand intron–boundary pair. It uses the existing assumption that both objects
measure the same unspliced RNA/DNA mixture, expressed with each component's opportunity.

The source contributes its observed strand split conditional on its observed total.
That total sets precision; the source's absolute abundance is not imposed on the target.
The target keeps its two Poisson observations and the existing total-RNA reference.
The source factor is evaluated **inside** integration over the target's unknown RNA
amount. Its arguments contain no count posterior, landscape, background or distant RNA
density. The permanent derivation is the conditional-strand evidence discussion in
[EQUATIONS.md](../EQUATIONS.md).

This permits a common capture multiplier to differ between the two objects: multiplying
both component rates at the source changes neither its strand probability nor the factor.
It does **not** prove that DNA and RNA actually have the same capture multiplier under
different fragment-length laws. The law-frame/opportunity issue remains open.

The reference has deliberately narrow scope. It does not certify splice-in/out factors,
both-strand neighbor messages, RNA level bounds, or a chain that repeats a fragment's
observation at several boundaries. At an unstranded source its factor is constant in
composition, as it must be. It is not the refused own-strand-only release reader.

## What the measurements establish

The settled numerical findings, including the six g50 conditions read separately, live
in [ISSUES.md](../ISSUES.md#the-gdna-prior-enters-psi-twice) under the factory-removal screen.
Every one of the 30 conditions, including the six zero-DNA controls, is tabulated in the
scratch receipt `calibration_tables.md`. Each cell is region / boundary **DNA-count**
error; these are not transcript / gene accuracy results.

Removal improves captured boundary estimates but costs information off capture.
Restoring actual intron strand messages recovers much of the stranded boundary loss.
The exact binomial messages remove a verified inference dependency; they do not make
the whole solver prior-free. Unstranded OFF remains worse than the frozen current tree.

A fixed-landscape control separates the two paths: most of that OFF loss persists with
the original landscape held fixed. It cannot be attributed merely to retraining the
landscape. This does not prove that a better explicit single-prior model cannot recover
it; it rules out assuming the kernel-training change alone will necessarily do so.

At this checkpoint no gap/ladder end-to-end or real-library campaign had started.
The owner subsequently accepted the modest test-chromosome OFF tradeoff and authorized
larger validation. That supersedes the earlier instruction to finish the entire evidence
model before advancing the candidate to those measurements. The continuation is recorded
in [RNA_FULL_PANEL_VALIDATION.md](RNA_FULL_PANEL_VALIDATION.md). The earlier reader's
comparisons remain in [RNA_SHORT_READER_REVIEW.md](RNA_SHORT_READER_REVIEW.md); that reader
and this count candidate are different experiments.

## Falsification and deliberate breakage

The native prototype passes 37 checks: factory removal and restoration, every sweep,
intron source wiring, independent binomial probabilities, strand reversal and rebuilt
claims under changed incoming beliefs. The original intron builder fails the three
source-wiring cases. The Gaussian message fails the exact-likelihood and belief-independence
checks before replacement.

The local density reference passes 18 checks against an independent closed-form
polynomial integral, including zero observations, an unstranded source, source absence,
strand reversal, component-rate units and a common source capture multiplier. Freezing
the RNA amount before reading the source fails the independent integral test.

Eight deliberate mutations/reversions fail their intended gates: restore the factory,
leave the patch installed, omit intron strand messages, restore the belief-frozen Gaussian,
freeze RNA outside the integral, use DNA opportunities for RNA, duplicate the source factor,
and change the RNA reference. The duplicate-factor test guards this local calculation;
it does not certify fragment independence across the entire graph.

The main-tree source/tests/goldens diff and installed calibration binary were verified
unchanged. The native worktree is `/private/tmp/rigel-density-20261007`; wheel builds were
unpacked into isolated sites rather than installed into the shared environment. The
drivers verify both native modules' import paths. Calibration threads stayed at their
defaults; the screens read scan caches and ran no EM.

## Execution from here

Continue the authorized count/evidence work, with the factory absent in the experimental
branch. Validate the frozen candidate on the full panels now; use measured failures to
order further model work. The remaining local factors matter particularly where strand
discrimination is unavailable:

- Certify the gene-edge DNA lower bound and a single junction's RNA constraint from their
  actual observations. Preserve their one-sided meanings. Distinguish a profiled bound
  from a marginalized likelihood; do not label either an observed DNA count.
- Derive one consistent splice-face factor before transporting it in both directions.
  The existing forward/reverse uncertainty rules are not certified by the arithmetic
  opportunity-map tests. State any required capture-transfer assumption explicitly.
- Track shared observations before combining factors. A factor cannot be multiplied twice
  merely because it arrived along two routes.
- Only then construct the profile-based landscape contrast and pair a successful partner
  with prior-once. Retain the existing training selection/weights for the first kernel
  contrast. Changing admission rules or introducing a new local biological assumption
  remains an owner discussion.

The factory-removal authorization is already resolved. Reader integration retains its
separate checkpoint. The aim remains a small general model, not acceptance of these
partial experiments as a release fix.

## Reproduction and exact artifacts

All artifacts are under `.cache/rigel_runs/2026-10-07_factory_removal/`:

| Artifact | Contents |
|---|---|
| `baseline.patch`, `baseline_head.txt` | Starting uncommitted tree and HEAD |
| `factory_ablation.py` | Scoped Python removal |
| `native_wiring.patch`, `native_observation.patch` | Separate native edits |
| `wheel/`, `site/`, `wheel_exact/`, `site_exact/`, build logs | Isolated binaries for the two native contrasts |
| `run_native.py` | Import-path enforcement, including the editable-finder check |
| `screen.py`, `test*_screen.json`, per-condition NPZs | Full and pass-zero calibration results, truth checks and binary hashes |
| `prior_feedback.py/json` | Fixed-landscape attribution and bit-identical controls |
| `local_face.py`, `test_local_face.py` | Local density reference and independent polynomial integral |
| `test_factory_ablation.py`, `test_intron_claim.py`, `test_observed_strand_claim.py` | Native/scoped-removal gates |
| `mutation_results.json`, `evidence_mutations.json`, logs | All deliberate failures |
| `summarize.py`, `calibration_tables.md`, `verification.json` | Full condition table and isolation verification |

The `current` arm name inside a raw native receipt means **that binary with the factory
enabled**. The frozen-current column in the rendered comparison always comes from the
original main-tree binary. Binary hashes distinguish every build.
