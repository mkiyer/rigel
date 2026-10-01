# ROADMAP — the short ranked view

**What this file is.** The one-line-per-claim state of the tool and the ordered next steps — nothing
else. Three rules keep it short: the substance of every item lives in `ISSUES.md` (the open entries plus
the append-only CLOSED / REFUSED record); the changelog is git, so this file records no history; and no
figure lives here — a claim names the instrument that re-derives it (owner, 2026-08-22). How performance
is judged is `SUCCESS.md`; rulings are `DESIGN.md`; lessons are `TRAPS.md`, cited by name.

## The frame (owner, 2026-09-22, 2026-09-23, 2026-09-26 and 2026-09-28)

The version on disk is `pyproject.toml`'s; the target is 0.8.0, A RELEASE OF THE TOOL (`DESIGN.md` §0b). Two numbers
are primary and answer different questions: the transcript table against per-transcript truth is what the release
ships on (`quant_accuracy.py`, per stratum, under fractional assignment), and the calibration result against an
oracle calibration is what ranks a calibration mechanism (`calibration_vs_oracle.py`,
`prior_vs_oracle.py`, with `ruler_vs_truth.py` beside every capture-ON arm, since the oracle calibration cannot see the
capture-contracted length move). An A/B pair runs with the scan pinned (`--set scan.total_threads=1`), so it is exactly
reproducible, and an effect is judged by its size, with genes and pools read beside the transcript table
(`TRAPS: the-deliverable-is-not-reproducible-by-default`). Three strata are in scope; unstranded × capture-ON is
deferred and never ranked on a pooled total (`TRAPS: never-pool-the-strata`). The fragment-length composition channel
is retired until after 0.8.0.

The cleanup is closed; what it left — dead, duplicated and narrated code among it — is `ISSUES: hygiene-ledger`, with
the correctness defects in `ISSUES: latent-defects`. The capture-contracted length of the locus gDNA component, the
synthetic spans and the annotated transcripts is kept as the one shared rule (owner, 2026-09-23, `DESIGN.md` §7.2);
robustness over synthetic accuracy is `DESIGN.md` §0b. On the ladder the in-scope work is at diminishing returns: the
junction price is accepted for now, the RNA floor at pure-gDNA objects is closed for 0.8.0, the pseudocount's odds are
fixed (`DESIGN.md` §3.1c), and which alignments gDNA may explain is ruled (`DESIGN.md` §3.1d). What remains shows on
real libraries and on the test chromosome's planted panels. The work runs on two tracks, each step derived on one page
and A/B'd against what ships, one issue closed at a time: on the cluster, **splicing artifacts on real data**, the
owner's highest priority (2026-09-26) — a real-data problem the simulated panels cannot show, handled across the
aligner, the blacklist's builder and Rigel, with Rigel's own defence the most important; locally, **strand
overdispersion and the gDNA fragment-length law**, the owner's (2026-09-28), then **the low-depth defects**, where
sparse libraries lose their gDNA.

## Where the tool is — one line per claim; run the named instrument for a current number

- **Stage A (the accumulator)**: done; the fragment ledger closes exactly — `calibration_oracle.py`; one open deposit
  defect, the leading intron (`ISSUES: latent-defects`).
- **Library gDNA fraction**: calibration's is accurate on the ladder's three in-scope strata at full depth and
  structurally blind on the deferred one — `policy_benchmark.py --by-class`; it falls with
  depth on sparse libraries (the landscape line below) and over-calls at low gDNA under capture
  (`ISSUES: capture-on-overcalls-gdna-at-low-gdna`); the EM's pseudocounts are neutral in their odds, and at `g98`
  the transcript table leans toward RNA through calibration's RNA floor and the capture likelihood —
  `quant_accuracy.py` (the pools per row); the RNA floor is closed for 0.8.0 until a fundamentally different
  algorithm (owner, 2026-09-26), `ISSUES: rna-prior-floor-at-pure-gdna-loci`.
- **The deliverable, end to end**: measured per stratum on the rebuilt ladder under the shared length; stranded ×
  capture-ON is the worst in-scope stratum and its capture-aware lengths own most of it — the junction price's
  structure, closed for now at diminishing returns (owner, 2026-09-26) — then the capture likelihood's lean toward RNA (the
  synthetic spans over-called at `g98`) and calibration's RNA floor at `g98` — `quant_accuracy.py --arm base` beside
  `--arm oracle_ruler` and `--arm oracle`, `ISSUES: the-capture-length-owns-stranded-capture-on`.
- **Fragment lengths**: the scorer's two length tables are censuses that average capture over different placements,
  and the capture ruler counts capture a second time on RNA's; one uncaptured frame is the largest lever measured on
  stranded × capture ON, but it loses on real length gaps until the ruler prices RNA's capture from gDNA's
  length-resolved mass, which Rigel does not collect — parked past 0.8.0 (owner, 2026-09-30),
  `ISSUES: the-scorer-reads-a-census-length-law`; the boundary inversion's second wave stays open; the RNA law trains on
  spliced fragments that carry splice artifacts (`ISSUES: splicing-artifacts`) — `calibration/fl.py`,
  `gdna_density.py`, `calibration_vs_oracle.py`; watch `ISSUES: capture-degeneracy-standing-risk`.
- **Strand model**: od = 0 by policy (`DESIGN.md` §3.3a) and κ from the genuine junctions (§3.3b), both landed
  2026-10-01 after their VCaP A/Bs (LBX0588's κ 0.064 → 0.0030). The `summary.json` diagnostics and the pruning design
  are next — Tier 0, `ISSUES: strand-overdispersion-one-shared-value`; the prototype harness and the VCaP truth scorer are in
  `~/Downloads/rigel_runs/prototypes/2026-09-30_robust_od/`.
- **The message layer**: `transfer` ships on the two-phase backbone (`DESIGN.md` §6b.12–§6b.14); `silent` is the
  floor; the bar — win on unstranded, minimal harm on stranded, never pooled — is `policy_benchmark.py --panel
  ladder`.
- **The gDNA landscape prior**: right on the ladder at full depth (`DESIGN.md` §7.1); on sparse libraries its refits
  collapse to the grid floor, so calibration's gDNA falls with depth and the EM follows it — Tier 1,
  `ISSUES: the-gdna-landscape-collapses-at-low-depth`, bounded by its zero controls
  (`ISSUES: gdna-landscape-trains-on-false-positives`), with the background's dispersion
  (`ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool`) and the shallow objects' plug-in
  (`ISSUES: strand-plug-in-bias-on-sparse-libraries`) beside it — `calibration_vs_oracle.py --suite` on the test
  chromosome's depth family.
- **ψ**: the composition closes structurally on every published object (`test_vertex_reference.py`); the tilt is
  integrated on derived nodes (`EQUATIONS.md` §9e–§9f); the λ bracket follows the landscape prior's demand; where the
  protocol calls the strand channel dead ψ still reads κ̂ — `policy_benchmark.py --by-class` (unstranded rows),
  `ISSUES: psi-reads-kappa-where-the-strand-channel-is-dead`.
- **The prior assembler**: with perfect masses its one measured error is the pooled crossing share — `prior_vs_oracle.py`
  (`O − S`, `ISSUES: the-pooled-q-in-the-gdna-count`); a perfect `LocusPriors` is worth little in scope —
  `quant_accuracy.py --arm oracle`; the per-transcript lane the EM never receives is worth far more —
  `quant_accuracy.py --arm oracle_alloc_seed`, `ISSUES: per-transcript-prior-lane`.
- **The capture-contracted length**: one shared rule for every EM component — each object's conserved share at
  that object's own capture efficiency, read against the landscape's located enriched mode, a junction priced
  from its neighbours by conservation of bases (`DESIGN.md` §7.2, `EQUATIONS.md` §11); the classes sit near one
  scale and the junction price's within-gene error is its structure, which no neighbouring gDNA object sees — closed
  for now (owner, 2026-09-26), `ruler_vs_truth.py --scale` (the class means and the within-gene spread),
  `ISSUES: the-junction-price-is-noisy-within-a-gene`; what the gDNA witness cannot see of a
  transcript-designed panel is declared — `ISSUES: ruler-witness-geometry-on-transcript-panels`.
- **Splice-artifact detection**: the blacklist errs both ways on real data — artifacts pass as certified RNA and
  rejected genuine reads become gDNA-eligible — and an index whose manifest records no alignable store runs with
  detection off and no error — the cluster track, `ISSUES: splicing-artifacts`,
  `ISSUES: an-unrecorded-splice-blacklist-is-dropped-silently`; `summary.json`'s `sj_blacklist_loaded`.
- **Performance**: the sweep is one native call, bit-identical at every thread count; the review's first two
  numeric no-ops are unparked for Tier 4 and the rest of the thread stays parked —
  `ISSUES: performance-memory-bounded-solve`, `profiling/profiler.py`, `profiling/sweep_replay.py`; a `_debug`
  capture has no memory bound, so one real-genome job runs at a time (`ISSUES: debug-capture-memory-is-unbounded`).
- **Panels**: the 16-condition ladder is rebuilt under the corrected capture physics, cached and certified, one
  realization per condition; the test chromosome is cached and certified; the fl-gap side panels are to be re-simulated and the
  junction-probed twin retired — `panel.py status`, `ISSUES: flgap-panels-stale-nascent-model`,
  `ISSUES: expand-the-gdna-spectrum`.
- **Reading rules**: rank per stratum, never pooled (`TRAPS: never-pool-the-strata`); quote `mwae` / Σ|err| over every
  object with mass (`calibration_vs_oracle.py`), never an intermediate (`TRAPS: the-intermediate-is-not-the-deliverable`).

## Next — the order

Ranked for the release, most critical first (owner, 2026-09-30): the strand overdispersion is finished before
anything else, the cluster track included; then what can ship wrong on real data, or cannot change once released,
comes before the in-scope accuracy work, which is at diminishing returns on the ladder. `ISSUES.md`'s
OPEN section follows this order within each priority. The cluster track does not compete for local A/B windows; the
local track takes one A/B window per mechanism, and a numeric no-op, proven with `rename_identity.py --check`, lands
between any two windows.

**Release-critical — before 0.8.0.**
1. **Tier 0, strand overdispersion — first, before anything else (owner, 2026-09-30)**
   (`ISSUES: strand-overdispersion-one-shared-value`, the 2026-09-30 rulings), one step per window: od = 0 and κ from the
   genuine junctions (both landed 2026-10-01); the `summary.json` diagnostics; the pruning design.
   After it: re-read the strand-input drift in `ISSUES: the-gdna-landscape-collapses-at-low-depth` and
   `ISSUES: strand-likelihood-over-confident-beyond-od`, and refresh the ladder report, the issue-list page and the
   full_lowg explanation page, which lacks the od-arm line of
   `ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool`.
2. **Splicing artifacts on real data** — the cluster track, after the strand overdispersion
   (`ISSUES: splicing-artifacts`), in phases: the substrate — the local-only VCaP mix and cfRNA copies confirmed on the
   cluster and uploaded where missing, a working index built after the index format bump, truth mixes, and the splice
   evidence regenerated from scratch by a new run, a separate task; both errors measured on today's tree; the mechanism
   census; the h-weighted reading; training sets weighted by each fragment's probability of being genuine, which
   unblocks strand-overdispersion step 4; aligner robustness and the catalogue rebuild. Locally, independent of the
   phases: the reject-rule falsification tests and the three-fragment-types A/B.
3. **The index format bump, before the cluster's working-index build** — so the cluster rebuilds once
   (`ISSUES: the-format-changes-to-batch-before-release` (b)); its `summary.json` half lands with step 1 of the
   strand overdispersion, whose one od field it carries.
4. **A dropped blacklist is never silent** — the production index's manifest checked before any cluster quant, and the
   one-line warning (`ISSUES: an-unrecorded-splice-blacklist-is-dropped-silently`).
5. **Tier 1, the low-depth defects** — the largest real-data effect found, after Tier 0 and A/B'd apart from it. In
   order: the failed refit alone (`ISSUES: the-gdna-landscape-collapses-at-low-depth`); the gDNA-only object row in
   `calibration_vs_oracle.py` (`ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool`); the unrun real-data
   arm; the landscape fit, derived and scored by class, without regressing its zero controls
   (`ISSUES: gdna-landscape-trains-on-false-positives`); the background dispersion, its own A/B; then read
   `ISSUES: capture-on-overcalls-gdna-at-low-gdna` and take up `ISSUES: strand-plug-in-bias-on-sparse-libraries`.
6. **User-facing correctness**, one commit each with a falsification test verified failing: the paths a user reaches in
   `ISSUES: latent-defects` first (the input parsers, the zero alphas, the finalizer deadlock), then its number-moving
   defects, each in its own window · `ISSUES: sj-strand-tag-chosen-from-the-first-reads` · a first read of
   `ISSUES: calibration-detects-capture-on-a-capture-off-library` on every capture-OFF panel, since a false reference
   contracts every transcript's length.
7. **The release gates**: the release-gating docs and the flaky reorder gate of `ISSUES: hygiene-ledger` · the two
   unparked numeric no-ops of `ISSUES: performance-memory-bounded-solve` · then the release — `docs/PUBLISHING.md` is
   the procedure; what gates it is the state: the deliverable measured per stratum, the zero rows clean, the suite at
   its standing count, `preflight.py --full` green, CI run once by hand, a real-data smoke run, the standing risks
   re-read, and the manual and the changelog true of what ships. The two residuals at `g98` — the capture likelihood's
   lean toward RNA and the pseudocount's strength — are small against the in-scope strata and are not held for the
   release.

**After the release-critical work — Tier 2, the in-scope mechanisms that move a primary number**, each its own
derivation and A/B: `ISSUES: nascent-stress-sensitivity`, a cheap re-measure of Tier 0's verdicts at the realistic
nascent level · `ISSUES: psi-reads-kappa-where-the-strand-channel-is-dead` · the fl second wave
(`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value`,
`ISSUES: the-fl-boundary-inversion-has-underived-pieces`, `ISSUES: capture-blind-gdna-divisor`,
`ISSUES: eb-shrinkage-magic-ess`) · `ISSUES: the-capture-reference-is-read-at-a-grid-point`, which does not land alone
· overdispersion steps 4–6, research, with or after the splice training sets
(`ISSUES: intron-seeds-near-probes-are-capture-enriched`, `ISSUES: strand-likelihood-over-confident-beyond-od`) ·
`ISSUES: the-pseudocount-strength-is-not-derived` · `ISSUES: per-transcript-prior-lane`, research with no candidate
mechanism yet · the message layer on unstranded capture-OFF (`ISSUES: message-layer-open-cases`,
`ISSUES: refit-vs-message-arbitration`, `ISSUES: the-intron-own-solve-on-unstranded-capture-off`) ·
`ISSUES: the-efficiency-posterior-floor-on-empty-pieces`.

**Later.**
- **Tier 3, real data only**, one real-genome job at a time. On the cluster after the splice phases:
  `ISSUES: scoring-penalties-are-underived-constants` with the h-weighted reading; the multimapper pair
  (`ISSUES: multimapper-intergenic-alignments`, `ISSUES: multimapper-blind-support`) and
  `ISSUES: unwitnessed-loci-and-multimappers-at-the-em` on the aligned ladder, regenerated with a larger simulation
  battery built for the cluster (owner, 2026-09-26). Locally: `ISSUES: vcap-dna-only-objects-lose-gdna-at-full-depth`
  · `ISSUES: em-gdna-exceeds-calibration-on-the-vcap-transcriptome-half`.
- **Tier 4, batched between A/B windows**; anything that moves a number takes a window of its own: the rest of
  `ISSUES: performance-memory-bounded-solve`, proven with `rename_identity.py --bam` · the rest of
  `ISSUES: hygiene-ledger`, the comment and doc sweep and the test gaps · `ISSUES: instrument-ledger` ·
  `ISSUES: debug-capture-memory-is-unbounded`, a src change · the panels (`ISSUES: flgap-panels-stale-nascent-model`,
  `ISSUES: expand-the-gdna-spectrum` with the junction-probed twin's retirement).

**Parked and deferred — Tier 5, each with its entry.** Deferred past 0.8.0: the deferred stratum, with
`ISSUES: a-pure-gdna-library-reads-as-nascent-rna` kept as an open challenge ·
`ISSUES: the-scorer-reads-a-census-length-law` with `ISSUES: the-pooled-q-in-the-gdna-count` (owner, 2026-09-30) ·
`ISSUES: unannotated-transcription-is-booked-as-gdna`, a `DESIGN.md` ruling after 0.8.0 ·
`ISSUES: ruler-witness-geometry-on-transcript-panels` · `ISSUES: overlapping-synthetic-shadows`, an owner decision on
the index. Parked: `ISSUES: yield-variance-beside-the-count` · `ISSUES: capture-premise-untested-on-cdna` ·
`ISSUES: the-atom-at-an-unwitnessed-both-strand-slot` · `ISSUES: two-sided-exon-row` ·
`ISSUES: the-lower-bound-noise-ratchet` · `ISSUES: flux-floor-dispersion` ·
`ISSUES: splice-out-premise-bias-uncorrected` · `ISSUES: the-tilt-census-as-an-instrument` ·
`ISSUES: transfer-variance-premise` · `ISSUES: drain-contaminates-certified-rna` · `ISSUES: crossing-pool-contrast` ·
`ISSUES: capture-degeneracy-standing-risk` · `ISSUES: pure-rna-mirror-asymmetry` ·
`ISSUES: binary-cuts-on-continuous-quantities`.

## Deliberately not next

The length composition channel (retired until after 0.8.0) · one uncaptured frame for the length tables (parked with its
entry) · a capture efficiency the EM re-reads as it runs (deferred by the owner, `DESIGN.md` §7.2) · anything whose only
target is the deferred stratum · every mechanism in `ISSUES.md`'s CLOSED / REFUSED section — read it before proposing
anything, because each entry is a build that was measured and turned down, with the number that killed it.
