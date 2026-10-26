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
what cannot change once released before in-scope accuracy work, and its release-critical items are landed or measured;
since 2026-10-18 the release follows the owner's fix order, in the order the tool runs (below). Every change is derived
on one page and A/B'd against what ships, one issue closed at a time; **splicing artifacts on real data** stay the cluster's track — a real-data problem the
simulated panels cannot show, handled across the aligner, the blacklist's builder and Rigel.

## Where the tool is — one line per claim; run the named instrument for a current number

- **Count foundation**: the factory-free observation-based message implementation is integrated; conservative
  allocation in ambiguous weak-strand cases is owner-accepted, with bounds admission deferred —
  `rename_identity.py`, `calibration_vs_oracle.py`, `quant_accuracy.py`; `ISSUES: the-gdna-prior-enters-psi-twice`.
- **The capture reader**: mode-free (owner, 2026-10-13, after two external reviews) — every object's own gDNA density
  against the bulk, no detector (`DESIGN.md` §7.2); on gDNA-poor uncaptured libraries it reads calibration's false gDNA
  as capture — `ruler_vs_truth.py`, `quant_accuracy.py`, `ISSUES: the-reader-reads-false-gdna-as-capture`.
- **The scan**: every usable hit is one molecule of both mates with a fragment length, paired by the search for every
  aligner (HI not read; `DESIGN.md` §1.1); each read pair has one outcome and each record one status, counted in
  `summary.json`'s `read_pairs` and checked against the annotation pass (`DESIGN.md` §1.2) — `tests/test_hit_verdict.py`,
  `tests/test_accounting.py`; a multimapper's intergenic hits never reach the EM
  (`ISSUES: multimapper-intergenic-alignments`).
- **Stage A (the accumulator)**: done; the fragment ledger closes exactly — `tests/test_scan_order_independence.py`; one
  open deposit defect, the leading intron (`ISSUES: latent-defects`); the second pass's four corrections are open
  (`ISSUES: the-second-pass-corrections`).
- **Library gDNA fraction**: calibration's matches the oracle calibration's on the ladder's three in-scope strata at full
  depth and is structurally blind on the deferred one — `calibration_vs_oracle.py` (its pools); end to end, with the
  length law read over the observable genome (2026-10-22) it reads 0.96–1.00 of truth on every in-scope row of the fair
  panel, the STAR rung and both fl-gap panels, where it under-called by up to 28 % and booked the rest as nascent RNA,
  and 0.2526 on the VCaP mix against 0.2518 — `quant_accuracy.py` (the library rows),
  `ISSUES: region-lengths-count-unobservable-bases`; it falls with
  depth on sparse libraries (the landscape line below) and over-calls at low gDNA under capture
  (`ISSUES: capture-on-overcalls-gdna-at-low-gdna`); the EM's pseudocounts are neutral in their odds, and at `g98`
  the transcript table leans toward RNA through calibration's RNA floor and the capture likelihood —
  `quant_accuracy.py` (the pools per row); the RNA floor is closed for 0.8.0 until a fundamentally different
  algorithm (owner, 2026-09-26), `ISSUES: rna-prior-floor-at-pure-gdna-loci`.
- **The deliverable, end to end**: measured per stratum on the rebuilt ladder under the shared length (the report card
  re-derives it, from `quant_accuracy.py`'s arms), on the simulator's alignments and on STAR alignments of the same
  reads (the ladder's STAR rung, `panel.py align`, `quant_accuracy.py --bam`); stranded × capture-ON has been the
  worst in-scope stratum and its capture-aware lengths own most of it — the junction price's structure, closed for now
  at diminishing returns (owner, 2026-09-26) — then the capture likelihood's lean toward RNA (the synthetic spans
  over-called at `g98`) and calibration's RNA floor at `g98` — `quant_accuracy.py --arm base` beside `--arm
  oracle_ruler` and `--arm oracle`, `ISSUES: the-capture-length-owns-stranded-capture-on`.
- **Fragment lengths**: the EM's two length tables are one shared law unless the data show a difference capture
  cannot make (`DESIGN.md` §3.6b, landed 2026-10-26): every gap row kept whole and every full-depth equal-chemistry row
  one law, stranded × capture-ON genes −4.6 %, the deferred unstranded × capture-ON −6.8 % with its pools at the
  previous scorer's, measured pinned and fractional against the previous scorer —
  `$RIGEL_RUNS/prototypes/2026-10-25_shared_law/` (`prod_arm.py`, `score.py`), `quant_accuracy.py` per stratum; the
  coupled gDNA law was deleted 2026-10-26 with its gDNA-long regression recorded, and the frame correction that
  recovers it is the strand frame, `ISSUES: the-strand-frame-under-the-shared-law`; the boundary inversion's second
  wave stays open; the RNA law trains on spliced fragments that carry splice artifacts (`ISSUES: splicing-artifacts`)
  — `calibration/fl.py`,
  `calibration/shared_length_law.py`, `gdna_density.py`, `calibration_vs_oracle.py`; watch
  `ISSUES: capture-degeneracy-standing-risk`.
- **Fragment-length gaps** (re-recorded 2026-10-26 under the shared length law, pinned and fractional): the shared law
  keeps every gap row's difference whole (suite RNA-long and RNA-short, test gDNA-long and RNA-long, capture on); the
  stranded gap rows keep less than a strand frame would (`ISSUES: the-strand-frame-under-the-shared-law`), and the RNA-long
  rows at g05 credit capture with as much of the gap as capture could make (+8 % genes) — `quant_accuracy.py` on both gap
  panels beside the equal-length ladder, `ISSUES: the-deferred-rna-long-regression-after-factory-removal`.
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
- **Panels**: the 16-condition ladder and both fl-gap side panels were wiped and re-simulated on 2026-10-21 on the fair
  simulator (no fragment over an `N` run, errors at 0.002 per base with NM, rounded lengths), cached and certified, one
  realization per condition (the gDNA capture-OFF rows certify composition only: their field gate reads the region
  lengths' unknown bases), and the ladder aligned by rustar for its STAR rung; the test chromosome and the
  junction-probed twin, which hold no `N`, are cached and certified as before — `panel.py status`, `docs/TESTING.md`,
  `ISSUES: expand-the-gdna-spectrum`.
- **Reading rules**: rank per stratum, never pooled (`TRAPS: never-pool-the-strata`); quote `mwae` / Σ|err| over every
  object with mass (`calibration_vs_oracle.py`), never an intermediate (`TRAPS: the-intermediate-is-not-the-deliverable`).

## Next — in the order the tool runs (owner, 2026-10-18 and 2026-10-21), the fair panel's defect first (owner, 2026-10-22)

The tool is fixed in the order it runs, one small change at a time, because a defect upstream confounds every fix
downstream: the capture reader has amplified every scan change (`ISSUES: the-reader-reads-false-gdna-as-capture`). Each
change gets a failing test on a hand-built case, the fix alone, each gate broken on purpose and seen to fire, the STAR
panel at the recommended parameters against truth per stratum (pinned and fractional, and while the reader stands with
every capture efficiency at 1 beside the shipped run), the oracle ladder for regressions, the suite and
`preflight.py --full`. THE GOAL (owner, 2026-10-21): gDNA that aligns as spliced is the main error on real data; the
recommended STAR parameters report the competing alignments and the multimapper handling is to weigh them, so that
splice artifacts are handled downstream (`ISSUES: splicing-artifacts`).

**First: the opportunity the fair panel exposed** (owner, 2026-10-22). The panels re-simulated on 2026-10-21 draw no
read from an `N` run, and on them Rigel under-calls gDNA by 6–21 % on the in-scope strata and books the rest as nascent
RNA. Every later step is measured on them, and the multimapper repair's gDNA term and the splice-artifact work both
weigh gDNA against RNA, so this is taken first; which bases are observable is a property of the reference, so its fix
starts at the index.
1. `ISSUES: region-lengths-count-unobservable-bases`: dissected 2026-10-22 — the term is the gDNA length law, whose
   one-sided rate reads unobservable bases as the cleanest gDNA. The owner's direction (2026-10-22): an observable
   genome from mappability, a mask first; the mask is approved: with an `alignable` store its zero-mappability set at
   the store's longest read length, N included; without one, the N runs; the existing GRCh38 store serves for now.
   Built 2026-10-22, one step at a time: the index records the mask, the second pass and the length law read it, and
   every in-scope row reads its gDNA at 0.96–1.00 of truth; the costs are traced in the entry. Next, one at a time:
   calibration's geometry (the oracle's gdna-field gate must pass after it; the census law's boundary reads went with
   the coupled law, 2026-10-26), the RNA lengths, and the report. Calibration's own count error is the larger half
   of the gDNA count's (`ISSUES: the-pooled-q-in-the-gdna-count`, the day's split): unstranded capture-OFF −4.0 %
   net on the ladder, g98 P/O 0.945.

**The index.**
2. `ISSUES: an-unrecorded-splice-blacklist-is-dropped-silently`: the cluster's manifest checked before any cluster
   quant, the owner's.

**The scan and its buffer.** The scan's repair landed between 2026-10-18 and 2026-10-21
(`ISSUES: the-scanner-builds-fragments-from-one-mate`, closed); what it leaves:
3. `ISSUES: multimapper-intergenic-alignments`: a multimapper's intergenic hits buffered and given a gDNA term, so its
   valid intergenic and genic alignments compete in the EM; release-blocking.
4. `ISSUES: splicing-artifacts`: the derivation, then the scan's infrastructure (junction ids on buffered fragments, the
   contiguous reading of a spliced fragment), the fixtures, a prototype and the VCaP A/B; the release waits for its
   validation on the cluster.
5. **Every panel re-cached once** with `--force`, and the identity references re-frozen, when the scan is final: the
   caches calibration is measured on date from 2026-10-16/17.

**The accumulator.**
6. `ISSUES: the-second-pass-corrections`: the length law, the weights, the mismatches and the motif rule, one arm each.

**Calibration**, its opportunity before its capture reader.
7. `ISSUES: the-gdna-length-contrast-is-ill-conditioned-on-aligned-input`.
8. `ISSUES: the-reader-reads-false-gdna-as-capture`: a chronic problem, approached only once the stages before it are
   resolved (owner, 2026-10-21); the planned fix reads each object's gDNA through calibration's own posterior for it.
9. `ISSUES: the-annotated-path-rule-costs-the-stranded-strata`, read again on the repaired reader.

**After calibration, before the EM**: the scores, the priors and the lengths the EM is handed.
10. In this order (owner, 2026-10-26): `ISSUES: the-pooled-q-in-the-gdna-count` — re-scoped the same day: the share
    is the small half (≤ 0.3 % of gDNA), calibration's own count error the large one and of the other sign, so the
    share lands with step 1's calibration geometry; then
    `ISSUES: the-strand-frame-under-the-shared-law`, the frame correction for stranded capture libraries, which must
    recover the coupling deletion's recorded gDNA-long regression; then
    `ISSUES: the-capture-likelihood-leans-toward-rna`, read beside it. The shared length law landed 2026-10-26
    (`DESIGN.md` §3.6b), closing `ISSUES: the-scorer-reads-a-census-length-law`, and the coupled gDNA law was deleted
    the same day (`ISSUES: the-coupled-gdna-law`).
11. `ISSUES: merge-keeps-least-overhang`.
12. `ISSUES: the-gdna-term-is-a-mean-over-places`, with `ISSUES: mismatch-gate-nan-in-multimapper-gdna` and the scorer
    half of `ISSUES: multimapper-intergenic-alignments`.
13. `ISSUES: the-em-prior-share-excludes-drawn-multimappers`.
14. `ISSUES: pruning-uses-three-references`.

**The release.**
15. `ISSUES: the-release-builds`, then `docs/PUBLISHING.md`'s procedure: the deliverable measured per stratum, the zero
    rows clean, the suite at its standing count, `preflight.py --full` green, CI run once by hand, a real-data smoke run,
    the standing risks re-read, and the manual and the changelog true of what ships.

**0.8.1 — first after the release**, each its own derivation and A/B, in the tool's order. The scan:
`ISSUES: clips-are-walked-along-the-genome-only` · `ISSUES: input-mistakes-pass-silently`. Calibration:
`ISSUES: the-gdna-prior-enters-psi-twice` · `ISSUES: psi-reads-kappa-where-the-strand-channel-is-dead` ·
`ISSUES: the-gdna-landscape-collapses-at-low-depth` with `ISSUES: strand-plug-in-bias-on-sparse-libraries` ·
`ISSUES: an-anchor-only-landscape-reads-as-capture` · `ISSUES: the-sparse-capture-losses` with
`ISSUES: the-efficiency-posterior-floor-on-empty-pieces` · `ISSUES: the-transfer-price-arithmetic` ·
`ISSUES: the-transfer-policy-loses-to-silence-on-unstranded-capture-off` with
`ISSUES: the-intron-own-solve-on-unstranded-capture-off`. Before the EM: `ISSUES: the-pseudocount-strength-is-not-derived`. The EM:
`ISSUES: the-em-stops-before-it-converges`. Across the stages: `ISSUES: latent-defects` · `ISSUES: instrument-ledger`.

**Later**, each with its entry, in the tool's order. The index: `ISSUES: overlapping-synthetic-shadows`. The scan:
`ISSUES: fragments-that-fit-one-gene-but-no-transcript` · `ISSUES: the-tail-rule-waits-for-truth`. The accumulator:
`ISSUES: drain-contaminates-certified-rna` · `ISSUES: draw-clamp-can-pick-a-zero` ·
`ISSUES: second-pass-geometry-in-python`. Calibration: `ISSUES: strand-overdispersion-one-shared-value` ·
`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value` ·
`ISSUES: the-fl-boundary-inversion-has-underived-pieces` · `ISSUES: capture-blind-gdna-divisor` ·
`ISSUES: eb-shrinkage-magic-ess` · `ISSUES: intron-seeds-near-probes-are-capture-enriched` ·
`ISSUES: the-edge-density-floor-under-capture` · `ISSUES: capture-on-overcalls-gdna-at-low-gdna` ·
`ISSUES: the-capture-reference-is-read-at-a-grid-point` · `ISSUES: the-bulk-can-be-calibrations-nascent-rna` ·
`ISSUES: vcap-dna-only-objects-lose-gdna-at-full-depth` · `ISSUES: unannotated-transcription-is-booked-as-gdna` ·
`ISSUES: the-deferred-rna-long-regression-after-factory-removal` · `ISSUES: two-own-composition-thresholds` ·
`ISSUES: message-layer-open-cases`. Before the EM: `ISSUES: per-transcript-prior-lane` ·
`ISSUES: unwitnessed-loci-and-multimappers-at-the-em` · `ISSUES: scoring-penalties-are-underived-constants` ·
`ISSUES: coverage-weight-plateau-for-short-transcripts` · `ISSUES: ruler-witness-geometry-on-transcript-panels` ·
`ISSUES: yield-variance-beside-the-count`. After the EM: `ISSUES: loci-table-recomputes-mrna`. The panels:
`ISSUES: nascent-stress-sensitivity` · `ISSUES: expand-the-gdna-spectrum` ·
`ISSUES: capture-probes-on-dropped-duplicates` · `ISSUES: simulated-errors-redraw-the-base`. The instruments:
`ISSUES: debug-capture-memory-is-unbounded`. Across the stages: `ISSUES: performance-memory-bounded-solve` ·
`ISSUES: hygiene-ledger`.

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