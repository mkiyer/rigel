# Session handoff — 2026-10-22

For the next session. Sandbox: nothing here is authoritative; the state is `ISSUES.md` and `ROADMAP.md`, the rulings
`DESIGN.md`. Read `CLAUDE.md` first, then this, then `ROADMAP.md`'s "Next" and the entries it names.

## Where things stand

- **The observable genome is built through its fourth step, uncommitted** (2026-10-22, after the owner's "Go"):
  `OBSERVABLE_OPPORTUNITY_DESIGN.md` §9 has every number. The index records the unobservable mask
  (`unobservable.feather`, format 11, `rigel.observable`: with a store its zero set at the longest bin, N included;
  without one, the FASTA's N runs); `region_arrays` cuts regions at it; the second pass and the length law read it.
  gDNA now reads 0.96–1.00 of truth on every in-scope row of the ladder, the STAR rung and both fl-gap panels (0.73–1.03
  before), VCaP 0.2526 against 0.2518 (0.2574 before, store-built index); the costs are traced in `ISSUES:
  region-lengths-count-unobservable-bases`. Four code patches and one docs patch, in order on `acc029bb`, are in
  `$RIGEL_RUNS/arms/2026-10-22_observable/patches/` with every receipt beside them. Nothing is committed; the owner
  drives commits. Deletions waiting for that commit: `gdna_opportunity_from_index` and `region_lengths_from_partition`
  (no production caller; their tests and `tests/calibration/_transfer_harness.py` still read them), and the `.v10` index
  directories kept beside the rebuilt format-11 ones (`suite`, `test_reference`, `refs`).
- **The benchmark data was pruned** (owner, 2026-10-22): 100 → 70 GB. No panel keeps FASTQs; `panel.py align` rebuilds a
  pair from the simulator's BAM byte for byte (60a48f4d). The deletions are listed under "Where things are".
- **The panels are fair now, and they show a large error the old ones hid** (2026-10-21 evening). The simulator no
  longer draws a fragment over an `N` run (18.8 % of the suite genome, where the old panels drew about a fifth of their
  gDNA off capture as all-`N` reads), its errors change the base (0.002 per base on the panels, NM on every oracle
  record), lengths are rounded, and probes on folded duplicates bind (ba9501c4; `ISSUES:
  the-simulator-draws-over-unknown-sequence` and its two companions, closed). The ladder and both fl-gap panels were
  wiped and re-simulated, cached and certified, aligned by rustar, scored, and the card republished. On them Rigel
  under-calls gDNA by 6–18 % off capture and 7–21 % on stranded capture-ON and books the rest as nascent RNA;
  transcript error on the three in-scope strata reads 3.23 / 2.34 / 5.49 % (2.36 / 2.02 / 3.88 % on the old panel). An
  A/B toggling only the mask proves the mask causes it, the perfect-prior arm loses as much, and calibration's library
  fraction matches the oracle's, so the loss lies past calibration's count: `ISSUES:
  region-lengths-count-unobservable-bases` (now, release-blocking, stage calibration) has the numbers. Not dissected.
  The STAR rung barely moved (its aligner never placed the all-`N` reads, so it read the defect all along); its g00
  unstranded capture-OFF control rose 6.6 → 29.5 % with its pools unchanged, also undissected.
- **The oracle certification**: every gDNA capture-OFF row's field gate now fails, by exactly its reference's known share,
  so those rows certify composition only. That is the region-length defect read by the gate, not a broken cache.
- **The cluster's STAR panel BAMs were deleted** with the re-simulation (their reads no longer exist); scan changes are
  decided on the rustar rung's BAMs (`TESTING.md` §0d) until the cluster aligns the new reads. The STAR-panel harness
  in `$RIGEL_RUNS/arms/2026-10-21_accounting/` (`run_cond.py`) still reads the deleted BAMs: point it at
  `<condition>/rustar_recommended.bam` before using it.
- **The scan's repair is landed and pushed** (the fix order's Phase A, items 1–15, 2026-10-18 → 2026-10-21): every
  usable hit is one molecule of both mates with a fragment length, paired by the search for every aligner (HI not
  read); chimeric and off-index alignments withheld; one length rule; each read pair one outcome and each record one
  status, counted (`summary.json` schema 5, `read_pairs`) and checked against the annotation pass. The record is
  `ISSUES: the-scanner-builds-fragments-from-one-mate` (closed) and `DESIGN.md` §1.1–§1.2; the commits are
  `git log 07eba0fc^..f763270c`. What the scan still owes: `ISSUES: multimapper-intergenic-alignments`, the scan half
  of `ISSUES: splicing-artifacts`, and (0.8.1) `ISSUES: clips-are-walked-along-the-genome-only`.
- **The suite's standing count is 3,414** (CLAUDE.md; 39 gates added by the observable genome's four steps); `preflight.py --full` green.
- **ISSUES.md was re-prioritized** (owner, 2026-10-21): every open entry names its `stage:`, and the OPEN section runs by
  tier, then in the order the tool runs. Closed in the audit: the one-mate entry and `output-names-are-not-consistent`
  (landed 1fa1c143 on 2026-10-17, still open by mistake). New: `the-second-pass-corrections`,
  `the-gdna-term-is-a-mean-over-places`, `clips-are-walked-along-the-genome-only`. Promoted to `now`, because the
  owner's fix order of 2026-10-18 puts them on the release path: `merge-keeps-least-overhang`,
  `pruning-uses-three-references`, `mismatch-gate-nan-in-multimapper-gdna`,
  `the-em-prior-share-excludes-drawn-multimappers` (the 0.8.0 cut of 2026-10-14 had kept its shipped denominator).
- **The report card was re-run on the fair panel** (ba9501c4, version 12 of
  https://claude.ai/artifact/Cek2wmKtitgbDyfM5gNqyj; receipts in `$RIGEL_RUNS/benchmarks/2026-10-21_resim/`, the VCaP
  run reused from `$RIGEL_RUNS/benchmarks/2026-10-21_card/`, Rigel's quantification unchanged since). Its 0.7.1 column is
  a ballpark from the old panel (owner: stop scoring 0.7.1 for now). Earlier the same day, on f763270c and the old
  panel: the oracle ladder moved by at most ±0.54 % per
  condition since 2026-10-14 (transcript error −0.13 / −0.13 / +0.09 % on the three in-scope strata): the simulator's
  alignments do not exercise the scan's repair, which the STAR panel measures. The card's calibration rows read the
  panels' scan caches, built 2026-10-16/17, and are unchanged to rounding. On the real VCaP mix the gDNA fraction
  reads 0.257 against its truth of 0.2518 (0.252 on 2026-10-14): the scan now withholds 214 k of its 18.7 M read pairs
  (unannotated splice 115 k, chimeric 67 k, unpaired 32 k), which changes what the fraction is a share of; not
  dissected. The issues page is https://claude.ai/artifact/Ah9RoZaM1sHAHoTXWRseFC.
- **The docs were refreshed**: the CHANGELOG (schema 5, `read_pairs`, the annotated BAM's statuses, the multimap option),
  the MANUAL's recommended STAR settings (`--sjdbScore 0 --outFilterMultimapScoreRange 2`) and its stale `NH`
  wording, `TESTING.md` §0d (the STAR panel), `SUCCESS.md`'s step 5 (a scan change). 94 settled sandbox notes were
  deleted (`git show f763270c:docs/dev/<name>`); the untracked review evidence was moved to
  `$RIGEL_RUNS/archive/docs_dev_untracked_2026-10-21/`.
- **Worktrees**: ten stale ones under the project's parent directory are the owner's to remove (each already on main;
  their uncommitted prototype edits were saved in `$RIGEL_RUNS/archive/worktree_patches_2026-10-21/`).
- ⚠ **Those two archive folders are gone**: found empty that evening (`$RIGEL_RUNS/archive/` held only the old panels'
  receipts written then), and Spotlight finds neither name. The patches named in step 2 of "Next" may be lost.

## The frame (owner)

- **Fix in the order the tool runs** (2026-10-18, restated 2026-10-21): index → scan and buffer → accumulator →
  calibration → after calibration, before the EM → EM → after the EM. A defect upstream confounds every fix downstream;
  the capture reader has amplified every scan change.
- **The goal** (2026-10-21): gDNA that aligns as spliced is the main error on real data. The recommended STAR settings
  report the competing alignments, and the multimapper handling is to weigh them, so that splice artifacts are handled
  downstream.
- **The capture reader** is a chronic problem, approached only once the stages before it are resolved.
- **Each change**: audit and design first for anything nontrivial (a design note in this directory, every claim about
  the code citing a re-read line, every number measured, the owner's approval); then a failing test on a hand-built
  case, verified failing on main; the fix alone; each gate broken on purpose and seen to fire; the STAR panel at the
  recommended settings against truth per stratum, pinned and fractional, with every capture efficiency at 1 beside
  the shipped run; the oracle ladder for regressions; `rename_identity.py` for a no-op; the suite and
  `preflight.py --full`; then commit, fast-forward main and push, as for items 1–15. Stage explicit paths only; the
  repository is public.

## Next, in order (ROADMAP.md's "Next")

1. `ISSUES: region-lengths-count-unobservable-bases` (owner, 2026-10-22: first). Steps 1–4 built; next, one at a
   time, as the design's §6 and §9 have them: 4b (the census law's boundary odds and exon on-target shape over reaches;
   16 boundaries on the suite, so likely below the A/B's noise — decide with the owner whether it is measured or only
   gated), 5 (calibration's geometry: the slot divisors, the landscape's anchors, the capture reader's supports, the
   priors' `S_r`; the oracle's gdna-field gate must pass after it), 6 (the transcripts' RNA lengths) and 7 (the
   report). Two findings to dissect beside it: the nascent → annotated leak under capture that step 4 exposed (oracle
   g50 stranded ON, about 50 k fragments; nearest `ISSUES: the-scorer-reads-a-census-length-law`), and the STAR rung's
   own store (cluster).
2. (owner) the cluster's index manifest, `ISSUES: an-unrecorded-splice-blacklist-is-dropped-silently`.
3. `ISSUES: multimapper-intergenic-alignments`, the scan half: buffer a multimapper's hits that overlap no transcript
   (the scan's intergenic path appends a unique mapper's alone). On its own it may be a numeric no-op for the EM, whose
   scorer skips a candidate-less hit, so the scorer's gDNA term for such a hit is its own step; how the group's term
   combines its places is `ISSUES: the-gdna-term-is-a-mean-over-places`, which needs each place's capture efficiency
   (the sum without it cost RNA multimappers 33.6 k on 2026-10-15). Prior art, on a pre-repair base: the patches named
   `rigel-mm-intergenic.*` and `rigel-mm-sum.*` were saved under `$RIGEL_RUNS/archive/worktree_patches_2026-10-21/`,
   which was found empty (see above); their measurements stay in `$RIGEL_RUNS/prototypes/2026-10-15_mm_intergenic/`
   and `2026-10-15_mm_sum/`.
4. `ISSUES: splicing-artifacts`: the derivation first, on paper. The design trail is here (`SPLICE_ARTIFACTS_*`,
   `GDNA_SPLICE_ARTIFACTS_*`) and the owner's own notes (`my_notes.md`: the RNA reading S and the gDNA reading C both
   propagate and compete in the EM).
5. Every panel re-cached once with `--force`, and the identity references re-frozen, when the scan is final. Open for
   the owner: whether that waits for step 4's scan infrastructure or comes after step 3.
6. The accumulator: `ISSUES: the-second-pass-corrections`, one arm each; and measure the reach of `ISSUES: latent-defects`'
   leading intron first.
7. Calibration, opportunity before the reader: the gDNA length law, then the reader (its planned fix: each object's
   gDNA read through ψ's `Var(log f_g)`), then the annotated-path rule's cost read again.
8. Before the EM: the merge, the gDNA sum with its NaN, the prior's denominator, pruning. Then the release.

## Open for the owner

- **Re-align on the cluster?** The deleted five-condition STAR panel could be re-aligned from the new reads at the
  production and recommended settings. The panels keep no FASTQs (2026-10-22): rebuild a pair from the condition's
  `sim_oracle.bam` with `samtools fastq -N`, byte for byte the simulator's.
- The tier promotions above, and the re-cache's timing (step 5).
- **Commit the observable genome's four steps?** One patch per step (with the regressions each traced: two isoform
  cases at g05, the capture-ON nascent → annotated leak step 4 exposed, the reader's g00 amplification on the STAR rung).
- **Which index the card's VCaP row reads.** The card reads `refs/rigel_index` (N-only mask, no blacklist); the owner's
  store gives `refs/rigel_index_alignable` (built 2026-10-22: the store's mask and its splice blacklist, 5.83 M
  junctions). On the store-built index VCaP reads 0.2574 → 0.2526 against 0.2518 with step 4; the blacklist is a change
  of its own, so a switch is measured alone.
- **Step 4b**: the census law's boundary reads touch 16 boundaries on the suite; measure it, or gate it only?
- Since this handoff was written the ladder has a STAR rung (all 16 conditions aligned by rustar with the recommended
  settings, `panel.py align`, scored by `quant_accuracy.py --bam`; `TESTING.md` §0d), on the report card. What stays
  out of tree is the capture-efficiencies-at-1 arm (`ISSUES: instrument-ledger` (j)).
- The MANUAL's new recommended-STAR paragraph and README line: the wording is new.
- `SUCCESS.md` says `panel.py score --jobs 2` (7–8.5 GB per condition), the ladder-report skill `--jobs 8`; the card ran
  at 8 without trouble on this machine.

## Where things are

| what | where |
|---|---|
| the STAR panel | the rung: `$RIGEL_RUNS/suite/ladder/<condition>/rustar_recommended.bam` (`TESTING.md` §0d); the cluster panel's text receipts in `$RIGEL_RUNS/star_panel/` (BAMs deleted) |
| the STAR-panel harness | `$RIGEL_RUNS/arms/2026-10-21_accounting/run_cond.py` (`NO_CAPTURE=1`: every efficiency at 1; `BASE=1`: the prototype site `RIGEL_PROTO_SITE` names), `tables.py`, `use_proto.py`, `identity_proto.py` |
| frozen identity references | on the fair panel: `$RIGEL_RUNS/arms/2026-10-21_resim/rename_ref_{oracle_g00,star_g50u}.json` (the STAR one from the rung's g50 κ0.50 capture-OFF BAM); the real libraries': `$RIGEL_RUNS/arms/2026-10-21_k_outer_end/rename_ref_{cfRNA-A,cfRNA-C}.json` |
| the accounting receipts | `$RIGEL_RUNS/arms/2026-10-21_accounting/` (`summary.md`, `stepN/measure.log`, `breaks_stepN/`) |
| the report card | built by the ladder-report skill from the arms; this run's receipts `$RIGEL_RUNS/benchmarks/2026-10-21_resim/` (`run.sh`, a log per stage, `card/cvo_*.json`) |
| the old panels | `$RIGEL_RUNS/archive/suite_before_resim_2026-10-21/` (text only: arms, truth tables, logs, cache manifests) |
| the mask A/B | `$RIGEL_RUNS/prototypes/2026-10-21_mask_ab/` (the unmasked arm's config, runner, BAMs and scores) |
| deleted 2026-10-22 (owner) | every panel's FASTQs (rebuilt on demand from the simulator's BAM); the local minimap2 panel, aligned from the old reads (its receipts stay in `$RIGEL_RUNS/prototypes/2026-10-15_splice_a/`); the closed campaigns' per-object dumps in `$RIGEL_RUNS/arms/` |
| leftover panel folders | `$RIGEL_RUNS/archive/leftover_panels_2026-10-22/` (`ladder_nrna_lo`, the test chromosome's `step3_*`) |

## The prompt this session ran (2026-10-22)

The next session's prompt follows the owner's review of `OBSERVABLE_OPPORTUNITY_DESIGN.md`.

```
Read CLAUDE.md, then docs/dev/SESSION_HANDOFF.md, then docs/ROADMAP.md's "Next" section and
ISSUES: region-lengths-count-unobservable-bases in full. We fix Rigel in the order the tool runs (index, scan,
accumulator, calibration, before the EM, the EM, after it), one small change at a time. I moved this defect first
(2026-10-22): the fair panel shows it is the largest error we know of, it distorts every later measurement, and
which bases are observable is a property of the reference.

Start with the dissection, as for items 13–15. Re-read and cite every place a region's length becomes an
opportunity, from the index through calibration, the priors and the EM. Then find the term that loses the gDNA on
g50 κ0.99 capture OFF and ON with the mask A/B pair the handoff names (mask on against mask off, the same seeds),
by node class and by locus, pinned and fractional, against truth: calibration's library fraction matches the
oracle's and the perfect-prior arm loses as much, so look past calibration's count. Then write a design note in
docs/dev for my review: what "observable" means (N runs first; a gap two mates span; repeats an aligner cannot
place), where it is recorded and where it is read, the steps, the failing tests and the measurement each step will
be judged on. Write no production code before I approve it.
```
