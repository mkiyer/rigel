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
- **Fragment lengths**: the realized gDNA law is fed the spliced RNA census unnormalised, so every exon-flanking
  boundary reads as pure gDNA — Tier 0, `ISSUES: the-realized-gdna-length-law-reads-rna-counts`, the boundary
  inversion's second wave after it; the RNA law trains on spliced fragments that carry splice artifacts
  (`ISSUES: splicing-artifacts`) — `calibration/fl.py`, `gdna_density.py`, `calibration_vs_oracle.py`; watch
  `ISSUES: capture-degeneracy-standing-risk`.
- **gDNA strand overdispersion**: the reconcile keeps gDNA's and RNA's values apart where the owner rules one shared
  value, and on a planted panel the shipped gDNA value sits near zero; `EQUATIONS.md` §6b–§6c and `DESIGN.md` §3.3a
  describe a fit that never shipped, §6a's no-bias claim is false, and `clamped_at_ceiling` and `effective_seeds` come
  from the discarded fit, so read neither until step 1 keys them on the shipped value — Tier 0,
  `ISSUES: strand-overdispersion-one-shared-value`, `calibration_vs_oracle.py` on the `odg05` panel.
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
  realization per condition; the test chromosome is cached and certified; the contaminated-seed panel is designed,
  not built (`ISSUES: strand-overdispersion-one-shared-value`); the fl-gap side panels are to be re-simulated and the
  junction-probed twin retired — `panel.py status`, `ISSUES: flgap-panels-stale-nascent-model`,
  `ISSUES: expand-the-gdna-spectrum`.
- **Reading rules**: rank per stratum, never pooled (`TRAPS: never-pool-the-strata`); quote `mwae` / Σ|err| over every
  object with mass (`calibration_vs_oracle.py`), never an intermediate (`TRAPS: the-intermediate-is-not-the-deliverable`).

## Next — the order

Two tracks run side by side, and `ISSUES.md`'s OPEN section follows this order within each priority. The local
track takes one A/B window per mechanism; a numeric no-op, proven with `rename_identity.py --check`, lands between
any two windows. The cluster track does not compete for local windows.

**The cluster track — now.**
- **Splicing artifacts, in phases** (`ISSUES: splicing-artifacts`): the substrate — the local-only VCaP mix and cfRNA
  copies confirmed on the cluster and uploaded where missing, a working index built after the index format bump
  (`ISSUES: the-format-changes-to-batch-before-release`, so the cluster rebuilds once), truth mixes, and the splice
  evidence regenerated from scratch by a new run, a separate task; both errors measured on today's tree; the mechanism
  census; the h-weighted reading; training sets weighted by each fragment's probability of being genuine, which
  unblocks strand-overdispersion step 4; aligner robustness and the catalogue rebuild.
- **The production index's manifest**, checked before any cluster quant
  (`ISSUES: an-unrecorded-splice-blacklist-is-dropped-silently`).

**The local track.**
1. **Tier 0, the owner's work — now**: `ISSUES: strand-overdispersion-one-shared-value` (steps 0–3, then the field
   collapse) and `ISSUES: the-realized-gdna-length-law-reads-rna-counts` (steps 1–4), one step per window: fl step 1,
   a numeric no-op, first; overdispersion step 0 — the contaminated-seed panel and the simulator's matched RNA
   overdispersion, outside the tree; step 1 under its hold rule; steps 2 and 3; fl steps 2, 3 and 4, step 2 never in
   an overdispersion window; the field collapse. After it: re-read the strand-input drift in
   `ISSUES: the-gdna-landscape-collapses-at-low-depth` and `ISSUES: strand-likelihood-over-confident-beyond-od`, and
   refresh the ladder report, the issue-list page and the full_lowg explanation page, which lacks the od-arm line of
   `ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool`.
2. **Tier 1, the low-depth defects — next**: in scope and the largest real-data effect found; after Tier 0 and A/B'd
   apart from it. In order: the failed refit alone (`ISSUES: the-gdna-landscape-collapses-at-low-depth`); the
   gDNA-only object row in `calibration_vs_oracle.py`
   (`ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool`); the unrun real-data arm; the landscape fit,
   derived and scored by class, without regressing its zero controls
   (`ISSUES: gdna-landscape-trains-on-false-positives`); the background dispersion, its own A/B; then read
   `ISSUES: capture-on-overcalls-gdna-at-low-gdna` and take up `ISSUES: strand-plug-in-bias-on-sparse-libraries`.
3. **Tier 2, the other in-scope mechanisms that move a primary number — later**, each its own derivation and A/B:
   `ISSUES: nascent-stress-sensitivity`, a cheap re-measure of Tier 0's verdicts at the realistic nascent level ·
   `ISSUES: psi-reads-kappa-where-the-strand-channel-is-dead` · `ISSUES: the-scorer-reads-a-census-length-law` with
   `ISSUES: the-pooled-q-in-the-gdna-count` · the fl second wave
   (`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value`,
   `ISSUES: the-fl-boundary-inversion-has-underived-pieces`, `ISSUES: capture-blind-gdna-divisor`,
   `ISSUES: eb-shrinkage-magic-ess`) · overdispersion steps 4–6, research, with or after the splice training sets
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
  battery built for the cluster (owner, 2026-09-26). Locally, independent of the phases: the reject-rule falsification
  tests and the three-fragment-types A/B in `ISSUES: splicing-artifacts` ·
  `ISSUES: sj-strand-tag-chosen-from-the-first-reads` · `ISSUES: vcap-dna-only-objects-lose-gdna-at-full-depth` ·
  `ISSUES: em-gdna-exceeds-calibration-on-the-vcap-transcriptome-half`.
- **Tier 4, batched between A/B windows**; anything that moves a number takes a window of its own:
  `ISSUES: latent-defects`, one commit each with a falsification test · `ISSUES: performance-memory-bounded-solve`,
  the numeric no-ops proven with `rename_identity.py --bam` · `ISSUES: hygiene-ledger`, the comment and doc sweep and
  the test gaps, the flaky reorder fixture first · `ISSUES: instrument-ledger` ·
  `ISSUES: debug-capture-memory-is-unbounded`, a src change ·
  `ISSUES: the-format-changes-to-batch-before-release`, before 0.8.0 while `summary.json` schema 3 is unreleased ·
  the panels (`ISSUES: flgap-panels-stale-nascent-model`, `ISSUES: expand-the-gdna-spectrum` with the junction-probed
  twin's retirement) · the release docs (`ISSUES: hygiene-ledger`).
- **The release** — `docs/PUBLISHING.md` is the procedure; what gates it is the state: the deliverable measured per
  stratum, the zero rows clean, the suite at its standing count, `preflight.py --full` green, CI run once by hand, a
  real-data smoke run, the standing risks re-read, and the manual and the changelog true of what ships. The two
  residuals at `g98` — the capture likelihood's lean toward RNA and the pseudocount's strength — are small against
  the in-scope strata and are not held for the release.

**Parked and deferred — Tier 5, each with its entry.** Deferred past 0.8.0: the deferred stratum, with
`ISSUES: a-pure-gdna-library-reads-as-nascent-rna` kept as an open challenge ·
`ISSUES: unannotated-transcription-is-booked-as-gdna`, a `DESIGN.md` ruling after 0.8.0 ·
`ISSUES: ruler-witness-geometry-on-transcript-panels` · `ISSUES: overlapping-synthetic-shadows`, an owner decision on
the index. Parked: `ISSUES: yield-variance-beside-the-count` · `ISSUES: capture-premise-untested-on-cdna` ·
`ISSUES: the-atom-at-an-unwitnessed-both-strand-slot` · the message layer's xfails (`ISSUES: two-sided-exon-row`,
`ISSUES: the-lower-bound-noise-ratchet`) · `ISSUES: flux-floor-dispersion` ·
`ISSUES: splice-out-premise-bias-uncorrected` · `ISSUES: the-tilt-census-as-an-instrument` ·
`ISSUES: transfer-variance-premise` · `ISSUES: drain-contaminates-certified-rna` · `ISSUES: crossing-pool-contrast` ·
`ISSUES: capture-degeneracy-standing-risk` · `ISSUES: pure-rna-mirror-asymmetry` ·
`ISSUES: binary-cuts-on-continuous-quantities`.

## Deliberately not next

The length composition channel (retired until after 0.8.0) · a capture efficiency the EM re-reads as it runs
(deferred by the owner, `DESIGN.md` §7.2) · anything whose only target is the deferred stratum · every mechanism
in `ISSUES.md`'s CLOSED / REFUSED section — read it before proposing anything, because each entry is a build that
was measured and turned down, with the number that killed it.
