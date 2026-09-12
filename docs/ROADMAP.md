# ROADMAP — the short ranked view

**What this file is.** The one-line-per-claim state of the tool and the ordered next steps — nothing
else. Three rules keep it short: the substance of every item lives in `ISSUES.md` (the open entries plus
the append-only CLOSED / REFUSED record); the changelog is git, so this file records no history; and no
figure lives here — a claim names the instrument that re-derives it (owner, 2026-08-22). How performance
is judged is `SUCCESS.md`; rulings are `DESIGN.md`; lessons are `TRAPS.md`, cited by name.

## The 0.8.0 frame

The version on disk is `pyproject.toml`'s; the target is 0.8.0, a calibration release, and the metric is
the calibration result scored against oracle calibration (`calibration_vs_oracle.py`,
`solvability_audit.py`, `prior_vs_oracle.py`) — the transcript number is a thermometer, never the
ranking. Three strata are in scope (unstranded × capture-OFF, stranded × capture-OFF, stranded ×
capture-ON); unstranded × capture-ON is deferred — reported on every benchmark, never a development
target, and never ranked on a pooled total (`TRAPS: never-pool-the-strata`). The fragment-length
composition channel is retired until after 0.8.0. The full ruling, including why the ladder gives gDNA
and RNA equal fragment lengths, is `DESIGN.md` §0b.

## Where the tool is — one line per claim; run the named instrument for a current number

- **Library gDNA fraction**: accurate on the three in-scope strata, structurally blind on the deferred
  one (at κ = ½ no channel reaches an AMBIG slot; the θ-independent-channel search is closed) —
  `solvability_audit.py`, `policy_benchmark.py --by-class`.
- **Transcript assignment**: a large share of RNA fragments is misassigned even under a perfect prior —
  calibration and assignment are two problems in two files; in scope a perfect prior no longer improves
  the transcript number, and the `g00` rows carry the largest transcript error of any stratum under
  both arms (the ruler, below) — `quant_accuracy.py` (the thermometer).
- **Stage A (the accumulator)**: done; the fragment ledger closes exactly — `calibration_oracle.py`.
- **Fragment lengths**: closed, both halves — gDNA by the two-pool contrast (`calibration/fl.py`,
  `gdna_density.py`; gates `test_fl.py`, `test_gdna_density.py`), RNA sound as shipped
  (`ISSUES: the-rna-length-law-fix`, CLOSED). Watch: `ISSUES: capture-degeneracy-standing-risk`.
- **gDNA strand overdispersion**: robust to the annotation (`EQUATIONS.md` §6a–§6c, `DESIGN.md` §3.3a);
  on real data read `clamped_at_ceiling` and `effective_seeds`, never the bare value.
- **The message layer**: `transfer` ships on the two-phase backbone with the level lanes
  (`DESIGN.md` §6b.12–§6b.14); `silent` is the measured floor. The bar — win on unstranded, minimal harm
  on stranded, never pooled — is `policy_benchmark.py --panel ladder`; the zero rows are solved by the
  prior under both policies, so the "beats silence" count is read on the contaminated rows, where every
  row favours `transfer`; `calibration_walk.py` says the messages still carry the stranded capture-ON
  rows and are essential on the deferred stratum.
- **The gDNA landscape prior**: done for 0.8.0 (`DESIGN.md` §7.1); the zero controls are solved on the
  metric (`calibration_vs_oracle.py`; `landscape_training_census.py` re-derives the population) and the
  in-scope per-object composition error did not move.
- **ψ**: the composition closes structurally on every published object (`test_vertex_reference.py`);
  the reference location is deleted (`DESIGN.md` §6b.1); the λ-bracket widening is built and ships off
  — `ISSUES: psi-lambda-bracket-unshipped`.
- **The prior assembler**: with perfect masses its own error is negligible — `prior_vs_oracle.py`,
  `mass_prior_ab.py`.
- **The largest number on the metric page is the ruler, not the composition**: at `g00` the
  effective-length shrinkage fabricates a reference from the residual false-positive fragments and
  contracts every transcript (`ISSUES: g00-shrinkage-upstream-repair` — the fix is the detector); the
  never-passed per-transcript prior lane (`ISSUES: per-transcript-prior-lane`) is the other.
- **Calibration's performance**: the one unfinished component — the sweeps dominate a deep run, on a
  single core, while the locus EM beside them is a rounding error. The decomposition is built: a
  terminal receives nothing, the sweep solves the chain a locus block at a time, the block size moves no
  number (`DESIGN.md` §6b.15); what remains is the C/C++ port of the block solve, where the parallelism
  goes — `profiling/profiler.py`, `profiling/sweep_replay.py --block-slots`.
- **Panels**: the sparse-nascent 16-condition ladder and the 30-condition test chromosome, both cached
  and certified — `panel.py status`; the fl-gap side panels carry a different nascent model —
  `ISSUES: flgap-panels-stale-nascent-model`. The ladder's nascent level is a development stress
  (`DESIGN.md` §0b).
- **Oracle FIELD certification**: every ladder row is stamped, but the uniformity gate is vacuous on
  capture-ON and zero-gDNA rows — read the stamp with its vacuity flag — `calibration_oracle.py`.
- **Attribution floor**: the deliverable is not reproducible by default; no `quant_accuracy` delta
  below the reseed floor is attributable — re-derive `--arm base_reseed` in the same session.
- **Reading rules**: rank per stratum; quote `mwae_all` / Σ|err| and the shipped column, never `solv%`
  or pass-0 (`TRAPS: the-intermediate-is-not-the-deliverable`).

## Next — the recommended order

The method is the dissection loop: run the panel → worst in-scope scenario → rank its objects by error
mass (`worst_objects.py`, `calibration_walk.py`) → find the mechanism → gated fix → add the offending
transcripts to the test chromosome → re-run → repeat. The facts this ranking leans on, each named with
its instrument: the ruler's factor at `g00` is the largest in-scope number on the metric page
(`calibration_vs_oracle.py`); a perfect prior is worth nothing in scope end to end
(`quant_accuracy.py`); the vertex ceiling is small on every stranded row and larger on the unstranded
capture-OFF rows, rising with gDNA (`vertex_ceiling.py`); by class the in-scope residual sits on the
intron's own solve (unstranded OFF) and on exon|exon boundaries and walled exons (stranded ON)
(`policy_benchmark.py --by-class`).

1. **Calibration's performance — `ISSUES: performance-memory-bounded-solve`** (owner, 2026-09-11: the
   active thread). The locus decomposition is landed and gated (`DESIGN.md` §6b.15); the agreed order
   from here, each step judged by `profiling/profiler.py --compare` on back-to-back pairs and proven a
   no-op with `profiling/sweep_replay.py` (`--block-slots` for the chunk-exactness of the whole sweep) and
   `design/rename_identity.py --bam` against the `locus_identity_*` references:
   ⓪ re-measure the deep library end to end, `main` against the landed tree — the baseline the rest is
   judged against; ① cache the refit-invariant half of `prepare` across the refit sweeps (the face
   rules and lane faces read counts and geometry, only the own claims read the belief), per block, per
   grid; ② the policy's rules from closures to typed tables — the C/C++ data layout, written in Python
   first (done: `messages.transfer.Faces`); ③ the port of `sweep._solve_block`, the passes and `transfer_rows` first, then `prepare`, then
   ψ, then threads over blocks — with a DERIVED tolerance gate in place of bit-identity, since a language
   port cannot be bit-identical; ④ the intron-factory rows built per block, the last genome-wide arrays;
   ⑤ the scan and the second pass, the stages that scale with depth and the floor once the sweeps are
   compiled. The accuracy frame is unchanged, and no step may move a number.
2. **The ruler at zero gDNA — `ISSUES: g00-shrinkage-upstream-repair`.** A gDNA-free library is the
   modal real case, the composition there is now right, and the effective length the EM divides by is
   still a fraction of the truth because the reference-density detector accepts any few slots with
   positive mass. Derive what "this library has an enriched gDNA mode" is evidence of (a boolean),
   prototype it outside `src/` in `_global_reference_density`'s caller, judge on
   `calibration_vs_oracle.py`'s ruler table per stratum with both zero controls, then
   `quant_accuracy.py`. Settle first whether the instrument's "exactly 1.000 off capture" contract is
   stale, since both P and O read below it there.
3. **The rest of the pre-EM setup** — `priors.py` / `result.py` / `derive.py` against
   `prior_vs_oracle.py` (re-run it first) and the ruler column: `ISSUES: prior-fidelity-vs-deliverable`,
   `ISSUES: eb-shrinkage-magic-ess`, `ISSUES: capture-blind-gdna-divisor`,
   `ISSUES: per-transcript-prior-lane`, `ISSUES: u-ruler-arm`.
4. **The intron's own solve on unstranded capture-OFF** — the intron class carries the largest share of
   the in-scope error there (`policy_benchmark.py --by-class`): the factory profile's resolution against
   the intergenic background (`density_deconv`); dissect with `worst_objects.py`.
5. **The vertex atom** — priced by `vertex_ceiling.py` on silent genes and nascent-free introns; a
   mechanism for it is the prior's reference (`ISSUES: reference-prior-refuted-at-concept-level`
   constrains the form) or the intron's own solve, not a message.
6. **The message policy, only where a row is above the bar**: one prototype arm at a time through
   `policy_prototype.py --module`, halves apart, pass zero beside the pipeline:
   `ISSUES: flux-price-witness-units`, `ISSUES: two-sided-exon-row`, `ISSUES: flux-floor-dispersion`,
   `ISSUES: ambig-node-as-a-gdna-source`, `ISSUES: message-layer-open-cases`.

Then, in standing order: `ISSUES: scan-thread-split-starves-the-workers` ·
`ISSUES: refit-vs-message-arbitration` (re-read under the E-step: the walk now says the prior does the
unstranded rows and the messages the stranded capture-ON ones).

**Later / parked** (each has its entry): `expand-the-gdna-spectrum` · `psi-lambda-bracket-unshipped` ·
`transfer-variance-premise` · `nascent-stress-sensitivity` · `f32-strand-tilt-at-half` · `hygiene-ledger` ·
`oracle-effective-length-diagnostic` · `flgap-panels-stale-nascent-model` · `rename-the-drain` ·
`drain-contaminates-certified-rna` (the ceiling refused the in-solve correction; two recorded
follow-ups) · `the-cancelling-pair` (refused twice) · `crossing-pool-contrast` (blocked) ·
`parked-capture-pilot-sign` · `pure-rna-mirror-asymmetry` · `capture-degeneracy-standing-risk`.

## Deliberately not next

The length composition channel (retired until after 0.8.0) · anything whose only target is the deferred
stratum · every mechanism in `ISSUES.md`'s CLOSED / REFUSED section — read it before proposing anything,
because each entry is a build that was measured and turned down, with the number that killed it.
