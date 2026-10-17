# ROADMAP — the short ranked view

**What this file is.** The one-line-per-claim state of the tool and the ordered next steps — nothing
else. Three rules keep it short: the substance of every item lives in `ISSUES.md` (the open entries plus
the append-only CLOSED / REFUSED record); the changelog is git, so this file records no history; and no
figure lives here — a claim names the instrument that re-derives it (owner, 2026-08-22). How performance
is judged is `SUCCESS.md`; rulings are `DESIGN.md`; lessons are `TRAPS.md`, cited by name.

## The frame (owner, 2026-09-22, 2026-09-23, 2026-09-26, 2026-09-28 and 2026-10-14)

The version on disk is `pyproject.toml`'s; the target is 0.8.0, A RELEASE OF THE TOOL (`DESIGN.md` §0b). Two numbers
are primary and answer different questions: the transcript table against per-transcript truth is what the release
ships on (`quant_accuracy.py`, per stratum, under fractional assignment), and the calibration result against an
oracle calibration is what ranks a calibration mechanism (`calibration_vs_oracle.py`,
`prior_vs_oracle.py`, with `ruler_vs_truth.py` beside every capture-ON arm, since the oracle calibration cannot see the
capture-contracted length move). An A/B pair runs with the scan pinned (`--set scan.total_threads=1`), so it is exactly
reproducible, and an effect is judged by its size, with genes and pools read beside the transcript table
(`TRAPS: the-deliverable-is-not-reproducible-by-default`). Three strata are in scope; unstranded × capture-ON is
deferred and never ranked on a pooled total (`TRAPS: never-pool-the-strata`). Fragment-length robustness is
in scope; a separate fragment-length composition channel has not been selected.

The cleanup is closed; what it left — dead, duplicated and narrated code among it — is `ISSUES: hygiene-ledger`, with
the correctness defects in `ISSUES: latent-defects`. The capture-contracted length of the locus gDNA component, the
synthetic spans and the annotated transcripts is kept as the one shared rule (owner, 2026-09-23, `DESIGN.md` §7.2);
robustness over synthetic accuracy is `DESIGN.md` §0b. On the ladder the in-scope work is at diminishing returns: the
junction price is accepted for now, the RNA floor at pure-gDNA objects is closed for 0.8.0, the pseudocount's odds are
fixed (`DESIGN.md` §3.1c), and which alignments gDNA may explain is ruled (`DESIGN.md` §3.1d). What remains shows on
real libraries and on the test chromosome's planted panels. The 0.8.0 cut (owner, 2026-10-14) put what a user hits and
what cannot change once released before in-scope accuracy work; its release-critical items are landed or measured, and
the release now follows the owner's sequence below. Every change is derived on one page and A/B'd against what ships,
one issue closed at a time; **splicing artifacts on real data** stay the cluster's track — a real-data problem the
simulated panels cannot show, handled across the aligner, the blacklist's builder and Rigel.

## Where the tool is — one line per claim; run the named instrument for a current number

- **Count foundation**: the factory-free observation-based message implementation is integrated; conservative
  allocation in ambiguous weak-strand cases is owner-accepted, with bounds admission deferred —
  `rename_identity.py`, `calibration_vs_oracle.py`, `quant_accuracy.py`; `ISSUES: the-gdna-prior-enters-psi-twice`.
- **The capture reader**: mode-free (owner, 2026-10-13, after two external reviews) — every object's own gDNA density
  against the bulk, no detector (`DESIGN.md` §7.2); on gDNA-poor uncaptured libraries it reads calibration's false gDNA
  as capture — `ruler_vs_truth.py`, `quant_accuracy.py`, `ISSUES: the-reader-reads-false-gdna-as-capture`.
- **Stage A (the accumulator)**: done; the fragment ledger closes exactly — `tests/test_scan_order_independence.py`; one open deposit
  defect, the leading intron (`ISSUES: latent-defects`).
- **Library gDNA fraction**: calibration's is accurate on the ladder's three in-scope strata at full depth and
  structurally blind on the deferred one — `calibration_vs_oracle.py` (its pools); it falls with
  depth on sparse libraries (the landscape line below) and over-calls at low gDNA under capture
  (`ISSUES: capture-on-overcalls-gdna-at-low-gdna`); the EM's pseudocounts are neutral in their odds, and at `g98`
  the transcript table leans toward RNA through calibration's RNA floor and the capture likelihood —
  `quant_accuracy.py` (the pools per row); the RNA floor is closed for 0.8.0 until a fundamentally different
  algorithm (owner, 2026-09-26), `ISSUES: rna-prior-floor-at-pure-gdna-loci`.
- **The deliverable, end to end**: measured per stratum on the rebuilt ladder under the shared length (the report card
  of release step 3 re-derives it); stranded × capture-ON has been the worst in-scope stratum and its capture-aware
  lengths own most of it — the junction price's
  structure, closed for now at diminishing returns (owner, 2026-09-26) — then the capture likelihood's lean toward RNA (the
  synthetic spans over-called at `g98`) and calibration's RNA floor at `g98` — `quant_accuracy.py --arm base` beside
  `--arm oracle_ruler` and `--arm oracle`, `ISSUES: the-capture-length-owns-stranded-capture-on`.
- **Fragment lengths**: the scorer's two length tables are censuses that average capture over different placements,
  and the capture ruler counts capture a second time on RNA's; one uncaptured frame is the largest lever measured on
  stranded × capture ON, but it loses on real length gaps until the ruler prices RNA's capture from gDNA's
  length-resolved mass, which Rigel does not collect — 0.8.1 (owner, 2026-10-14),
  `ISSUES: the-scorer-reads-a-census-length-law`; the boundary inversion's second wave stays open; the RNA law trains on
  spliced fragments that carry splice artifacts (`ISSUES: splicing-artifacts`) — `calibration/fl.py`,
  `gdna_density.py`, `calibration_vs_oracle.py`; watch `ISSUES: capture-degeneracy-standing-risk`.
- **Fragment-length gaps** (re-recorded 2026-10-03 on the re-simulated suite arms, pinned and fractional): the test
  RNA-short OFF regression is repaired by the component-opportunity maps and maintained by the count foundation;
  captured RNA-short and deferred unstranded RNA-long remain material robustness concerns —
  `quant_accuracy.py` on both gap panels beside the equal-length ladder,
  `ISSUES: the-deferred-rna-long-regression-after-factory-removal` and
  `ISSUES: the-scorer-reads-a-census-length-law`.
- **Strand model**: od = 0 by policy (`DESIGN.md` §3.3a) and κ from the genuine junctions (§3.3b), both landed
  2026-10-01 after their VCaP A/Bs; `summary.json` reports the junction fit behind κ (landed 2026-10-14) —
  `ISSUES: strand-overdispersion-one-shared-value`; the splice-artifact pruning the fit points at is the strand prong
  of `ISSUES: splicing-artifacts`; the prototype harness and the VCaP truth scorer are in
  `$RIGEL_RUNS/prototypes/2026-09-30_robust_od/`.
- **The message layer**: `transfer` ships on the two-phase backbone (`DESIGN.md` §6b.12–§6b.14); `silent` is the
  floor; the bar — win on unstranded, minimal harm on stranded, never pooled — is `policy_benchmark.py --panel
  ladder`; on the count foundation it loses to silent on the unstranded × capture-OFF contaminated rows
  (`ISSUES: the-transfer-policy-loses-to-silence-on-unstranded-capture-off`).
- **The gDNA landscape prior**: right on the ladder at full depth (`DESIGN.md` §7.1); on sparse libraries its refits
  collapse to the grid floor, so calibration's gDNA falls with depth and the EM follows it —
  `ISSUES: the-gdna-landscape-collapses-at-low-depth`, with the shallow
  objects' plug-in (`ISSUES: strand-plug-in-bias-on-sparse-libraries`) beside it — `calibration_vs_oracle.py --suite`
  on the test chromosome's depth family.
- **ψ**: the composition closes structurally on every published object (`test_vertex_reference.py`); the tilt is
  integrated on derived nodes (`EQUATIONS.md` §9e–§9f); the λ bracket follows the landscape prior's demand; where the
  protocol calls the strand channel dead ψ still reads κ̂ — `policy_benchmark.py --by-class` (unstranded rows),
  `ISSUES: psi-reads-kappa-where-the-strand-channel-is-dead`.
- **The prior assembler**: with perfect masses its one measured error is the pooled crossing share — `prior_vs_oracle.py`
  (`O − S`, `ISSUES: the-pooled-q-in-the-gdna-count`); a perfect `LocusPriors` is worth little in scope —
  `quant_accuracy.py --arm oracle`; the per-transcript lane the EM never receives is worth far more —
  `quant_accuracy.py --arm oracle_alloc_seed`, `ISSUES: per-transcript-prior-lane`.
- **The capture-contracted length**: one shared rule for every EM component — each object's conserved share at
  that object's own capture efficiency, its gDNA density read against the bulk under the landscape as a two-group
  prior (owner, 2026-10-13), a junction priced from its neighbours by conservation of bases and bounded by the
  saturation level or the best object it touches, since capture saturates (`DESIGN.md` §7.2, `EQUATIONS.md` §11);
  the classes sit near one
  scale and the junction price's within-gene error is its structure, which no neighbouring gDNA object sees — closed
  for now (owner, 2026-09-26), `ruler_vs_truth.py --scale` (the class means and the within-gene spread),
  `ISSUES: the-junction-price-is-noisy-within-a-gene`; what the gDNA witness cannot see of a
  transcript-designed panel is declared — `ISSUES: ruler-witness-geometry-on-transcript-panels`.
- **Splice-artifact detection**: the blacklist errs both ways on real data — artifacts pass as certified RNA and
  rejected genuine reads become gDNA-eligible — and an index whose manifest records no alignable store runs with
  detection off, which its load now says — the cluster track, `ISSUES: splicing-artifacts`,
  `ISSUES: an-unrecorded-splice-blacklist-is-dropped-silently`; `summary.json`'s `sj_blacklist_loaded`.
- **Performance**: the sweep is one native call, bit-identical at every thread count; the locus EM hands loci to its
  threads one at a time; the two 0.8.0 numeric no-ops landed and the rest of the thread is later —
  `ISSUES: performance-memory-bounded-solve`, `profiling/profiler.py`, `profiling/sweep_replay.py`; a `_debug`
  capture has no memory bound, so one real-genome job runs at a time (`ISSUES: debug-capture-memory-is-unbounded`).
- **Panels**: the 16-condition ladder is rebuilt under the corrected capture physics, cached and certified, one
  realization per condition; the test chromosome is cached and certified; the fl-gap side panels and the junction-probed
  twin were re-simulated on the current simulator, cached and certified on 2026-10-03, with the ladder's nascent block
  and no FASTQs — `panel.py status`, `docs/TESTING.md`, `ISSUES: expand-the-gdna-spectrum`.
- **Reading rules**: rank per stratum, never pooled (`TRAPS: never-pool-the-strata`); quote `mwae` / Σ|err| over every
  object with mass (`calibration_vs_oracle.py`), never an intermediate (`TRAPS: the-intermediate-is-not-the-deliverable`).

## Next — the release sequence (owner, 2026-10-14)

The 0.8.0 cut's release-critical items were settled on 2026-10-14: the user-path fixes, the two numeric no-ops and the
format batch landed; the multimapper prior share was kept as shipped and the ψ double-prior fix was not shipped, each on
its A/B, their repairs in 0.8.1 (`ISSUES: the-em-prior-share-excludes-drawn-multimappers`,
`ISSUES: the-gdna-prior-enters-psi-twice`). The release now follows the owner's sequence, in order:
1. **Docs cleanup**: every permanent doc true of the tree. Done.
2. **An independent code review.** Done; what it found is in ISSUES, release-relevant first: `ISSUES: the-release-builds`,
   `ISSUES: the-em-depends-on-its-thread-count`; its smaller items are in `ISSUES: latent-defects` and
   `ISSUES: hygiene-ledger`.
3. **The full benchmark run**, the report card regenerated.
4. **ISSUES cleanup and prioritization**: the OPEN section ordered by priority. Done.
5. **The survey**: implement or defer, entry by entry. Done (owner, 2026-10-14): the release review's fixes are approved;
   the splice-artifact work is handled before any release, once the release candidate is clean
   (`ISSUES: splicing-artifacts`, which now holds the strand pruning); the EM's thread dependence waited for a measured
   prototype, and landed on it (`ISSUES: the-em-depends-on-its-thread-count`); everything else is deferred as tiered.
6. **The release issues fixed**, simple ones first and as many closed as possible, each with a falsification test
   verified failing; a numeric change gets its own A/B window, and a numeric no-op, proven with
   `rename_identity.py --check`, lands between any two windows. Then the involved ones: the EM's reduction, landed
   2026-10-14; then the splice-artifact mechanism (`ISSUES: splicing-artifacts`), its design externally reviewed
   2026-10-14. The derivation, the mechanism fixtures, a prototype and the VCaP A/B come next. The release waits
   for its validation on independent artifact-containing evidence on the cluster (owner, 2026-10-14).
7. **The release**: `docs/PUBLISHING.md` is the procedure; what gates it is the state — the deliverable measured per
   stratum, the zero rows clean, the suite at its standing count, `preflight.py --full` green, CI run once by hand, a
   real-data smoke run, the standing risks re-read, and the manual and the changelog true of what ships.

**0.8.1 — first after the release**, each its own derivation and A/B: `ISSUES: the-em-stops-before-it-converges` with
`ISSUES: the-pseudocount-strength-is-not-derived` and `ISSUES: the-em-prior-share-excludes-drawn-multimappers` ·
`ISSUES: the-reader-reads-false-gdna-as-capture` with `ISSUES: the-gdna-prior-enters-psi-twice` ·
`ISSUES: the-transfer-policy-loses-to-silence-on-unstranded-capture-off` with
`ISSUES: the-intron-own-solve-on-unstranded-capture-off` · `ISSUES: the-scorer-reads-a-census-length-law` with
`ISSUES: the-pooled-q-in-the-gdna-count` · `ISSUES: psi-reads-kappa-where-the-strand-channel-is-dead` ·
`ISSUES: the-gdna-landscape-collapses-at-low-depth` with `ISSUES: strand-plug-in-bias-on-sparse-libraries` ·
`ISSUES: the-sparse-capture-losses` with `ISSUES: the-efficiency-posterior-floor-on-empty-pieces` ·
`ISSUES: sj-strand-tag-chosen-from-the-first-reads` · `ISSUES: an-anchor-only-landscape-reads-as-capture` ·
`ISSUES: the-transfer-price-arithmetic` · `ISSUES: latent-defects` · `ISSUES: instrument-ledger` ·
`ISSUES: input-mistakes-pass-silently`.

**Later**, each with its entry: `ISSUES: the-edge-density-floor-under-capture` ·
`ISSUES: capture-on-overcalls-gdna-at-low-gdna` · `ISSUES: nascent-stress-sensitivity` ·
`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value` ·
`ISSUES: the-fl-boundary-inversion-has-underived-pieces` · `ISSUES: capture-blind-gdna-divisor` ·
`ISSUES: eb-shrinkage-magic-ess` · `ISSUES: the-capture-reference-is-read-at-a-grid-point` ·
`ISSUES: the-bulk-can-be-calibrations-nascent-rna` · `ISSUES: intron-seeds-near-probes-are-capture-enriched` ·
`ISSUES: per-transcript-prior-lane` · `ISSUES: message-layer-open-cases` ·
`ISSUES: scoring-penalties-are-underived-constants` · `ISSUES: multimapper-intergenic-alignments` ·
`ISSUES: unwitnessed-loci-and-multimappers-at-the-em` · `ISSUES: vcap-dna-only-objects-lose-gdna-at-full-depth` ·
`ISSUES: debug-capture-memory-is-unbounded` · `ISSUES: expand-the-gdna-spectrum` ·
`ISSUES: unannotated-transcription-is-booked-as-gdna` · `ISSUES: ruler-witness-geometry-on-transcript-panels` ·
`ISSUES: yield-variance-beside-the-count` · `ISSUES: drain-contaminates-certified-rna` ·
`ISSUES: the-deferred-rna-long-regression-after-factory-removal` · `ISSUES: overlapping-synthetic-shadows` ·
`ISSUES: performance-memory-bounded-solve` · `ISSUES: hygiene-ledger`.

**Parked — Tier 5**, each with its entry: `ISSUES: a-pure-gdna-library-reads-as-nascent-rna` ·
`ISSUES: capture-premise-untested-on-cdna` · `ISSUES: the-atom-at-an-unwitnessed-both-strand-slot` ·
`ISSUES: two-sided-exon-row` · `ISSUES: the-lower-bound-noise-ratchet` · `ISSUES: transfer-variance-premise` ·
`ISSUES: capture-degeneracy-standing-risk` · `ISSUES: binary-cuts-on-continuous-quantities`.

## Deliberately not next

A fragment-length composition channel (in scope, no mechanism selected) · one uncaptured frame for the length tables alone (it
loses on real length gaps; its entry names the next step) · a capture efficiency the EM re-reads as it runs (deferred
by the owner, `DESIGN.md` §7.2) · anything whose only target is the deferred stratum · every mechanism in
`ISSUES.md`'s CLOSED / REFUSED section — read it before proposing anything, because each entry is a build that was
measured and turned down, with the number that killed it.