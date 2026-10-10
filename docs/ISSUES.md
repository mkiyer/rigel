# ISSUES — the issue log

This file is the issue log: one entry per open problem, question, decision or risk, and an append-only record
of what was measured and turned down. An OPEN entry is a `### kebab-name` heading with a `priority:` line (now
/ next / later / parked), the question in a sentence or two, the numbers a ranking turns on, and the
instrument that re-derives them. The OPEN section runs now, then next, later and parked, and within a tier it
follows `ROADMAP.md`'s order. A CLOSED / REFUSED entry keeps its name, its verdict, the killing number and
the date, so a refused mechanism is not rebuilt. Cite an entry as `ISSUES: <name>`; names are the only
identifiers (`tests/test_docs_boundary.py`). What does not belong here: the ranked view (`ROADMAP.md`),
rulings and derivations (`DESIGN.md`, `EQUATIONS.md`), lessons (`TRAPS.md`), and any record of what was done —
the changelog is git.

---

## OPEN

### the-calibration-count-is-blind-to-multimappers
`priority: now — the owner's intended fix, after the native reader lands (2026-10-10); with the reader in place it is the one red gate in the suite · kind: defect · 2026-10-10`
THE SYMPTOM. `tests/scenarios_aligned/test_multimap_counting.py::TestParalogMultimapping::test_gdna_sweep[gdna_100]`
(two sequence-identical 500 bp single-exon paralogs, gDNA abundance 100, aligned reads): the shipped reader calls
154 / 0 against a truth of 75 / 71 (the collapse the test asserts); the capture reader calls 123 / 79, a total
of 202 against 146. Dissected: both paralog exons carry a calibration gDNA count of exactly 0 on a 155 bp
contained support beside intergenic gDNA at 0.23–0.31 fragments/bp, because every fragment inside them maps
twice and `bam_scanner.cpp` deposits into the calibration accumulator only for `is_unique_mapper` — a
multimapper goes to the EM's buffer alone. The capture reader reads the zero honestly (`exp(−ρ·Eg)` at the
background is `e^{−36}`), the weight is 870× below the flanking intergenic regions', `assemble_priors` contracts the locus's gDNA
yield to its boundaries', the EM reads the zero yield as "cannot emit", and the gDNA at the paralogs is called
RNA. The 123 / 79 split is the boundaries' crossing-count noise amplified by the vanished region weights, not a
tie-break, so the test's instruction to delete its collapse branch does not apply. The pre-existing half is the
shipped 154 against 146: the locus prior's gDNA count at a multimapper-only region is 0, so the EM has no gDNA
prior there; the located-mode reader amplified it the same way on every captured library and its detector hid
it on capture-OFF ones.

THE FOOTPRINT. No benchmark saw it: every panel is an oracle BAM (`NH:i:1` throughout) and the production
default is `include_multimap = True`. On the production index's 1,043,881 regions, per region the uniquely
mapped first-mate starts against the multimapper alignment starts (every alignment counted):

| library | unique fragments | multimapper alignments | regions with multimappers and no unique fragment (of them exon regions) | regions where multimappers outnumber unique fragments — share of all multimapper alignments inside them |
|---|---|---|---|---|
| LBX0190 (plasma) | 138,403 | 104,140 | 5,488 (762) | 6,504 — 95 % |
| MO_3021 (plasma) | 772,744 | 583,407 | 23,443 (7,018) | 35,781 — 94 % |
| VCaP mix (deep) | 18,417,903 | 544,065 | 10,099 (5,305) | 13,792 — 72 % |

Multimappers are up to 43 % of a plasma library's alignments and they cluster where calibration counts almost
nothing, which is exactly where the reader reads depletion and the locus prior reads no gDNA.

THE FIX (owner, 2026-10-10, a known issue): a change to the accumulator phase. The accumulator today deposits
uniquely aligned fragments for calibration and discards multimapping ones; it already buffers the fragments
compatible with several fragment lengths (several isoforms) and assigns them in a second pass. Multimapping
fragments are to be handled the same way — buffered, then assigned by the accumulator's second pass to the
probabilistically best placement by the abundance of uniquely aligned fragments there — which fixes
calibration and the capture reader at once; implementable without undue trouble, after the native reader is
complete. The accumulator's own ruling already said as much — "a LATER phase: side buffer, then deterministic
largest-remainder apportionment, integral always" — and the code implements it only for gap-hypothesis
fragments on one placement. Derived: a multimapper enters the side
buffer with one hypothesis per alignment (reference, start, end, introns), the second pass scores each
placement with its existing score (`ρ(h)` from pass one's tally at the placement's objects, the length law,
the strand term), draws one and the drain deposits it; calibration then counts the paralog exons at their
density, the reader reads them and the locus prior holds their gDNA. It is a format change to
`DeferredRecords` (per-hypothesis coordinates), a scanner change (defer instead of skip) and a drain change,
each gated by `tests/native/_accumulator_reference.py`, A/B'd on the aligned scenarios (the only panels with
multimappers) and on VCaP against its truth. The cheaper alternative — a per-region multimapper tally in the
scanner and a unique-deposit opportunity `Eg_r = S_r · U_r / (U_r + A_r / NH_r)` under uniform genomic
sampling — is a geometry estimated from data and needs its own ruling. Neither is a reader change: nothing in
the reader can tell a withheld fragment from an absent one.

### the-capture-weights-are-noisy-on-a-capture-off-library
`priority: next — the price of no detector, to be reduced only by a mechanism A/B'd on its own · kind: cost · 2026-10-10`
Every object's capture weight is its own gDNA density's posterior mode relative to the typical read object's, on every
library (`DESIGN.md` §7.2). On a capture-OFF field the true weights are all one level and the published ones
are its Poisson noise around that level: measured (`quant_accuracy.py`, pinned,
fractional, shipped → reader, transcripts % / genes %) test chromosome stranded × OFF 6.88 / 1.00 → 8.74 /
1.01, unstranded × OFF 7.99 / 1.13 → 8.41 / 1.19, ladder stranded × OFF 6.35 / 1.63 → 7.44 / 1.67; on the 21
golden scenarios (10–12 kb, a handful of objects, every one capture-OFF, scored against the simulator's truth)
the summed transcript error 298 → 324, worse on 7 and better on 5, the worst `combo_extreme` 6.7 → 22.1 of a
54-fragment truth. The noise is the estimator's, not the kernel's (the kernel is its specification to 1e-9),
and it is largest where objects are shallow. A reduction must be a mechanism of its own, A/B'd alone on both
capture-OFF strata and the capture-ON ones, never a detector, never a reference (`ISSUES: the-located-mode-capture-reader`).


Ordered by priority. An entry says what is open and the number a ranking turns on; what was done is git,
what was ruled is `DESIGN.md`.

### strand-overdispersion-one-shared-value
`priority: now — first (owner, 2026-09-30): od = 0 and κ LANDED 2026-10-01; next the summary.json diagnostics, then the pruning design · kind: design · 2026-09-27`
The strand model's two inputs, κ and the od, were both contaminated on real libraries:
- splice artifacts (misaligned gDNA, ½ in a stranded library) inflate κ and the RNA-side spread;
- antisense RNA inflates the gDNA side.

od = 0 by policy landed 2026-10-01 (`DESIGN.md` §3.3a, `EQUATIONS.md` §6). Its A/B against what shipped, with every
step proven to run its arm:
- **VCaP mix, against read-name truth:** gDNA pool −6.40 → −6.11 %, transcripts −0.15 %, genes −1.1 %.
- **Zero gDNA (VCaP RNA-only):** calibration false gDNA +6.1 %, EM +0.5 %.
- **Zero RNA (VCaP DNA-only):** unchanged.
- **LBX0588's gDNA total across subsamples:** 3.45 → 0.74 % at 10 %, 4.84 → 2.09 % at 1 %.
- **Ladder:** flat.
- **Test panel:** stranded × OFF transcripts −3.8 %.
- **odg05's g98 rows (planted gDNA od 0.05):** genes +1.1k / +3.1k and gDNA pool −10k, the one cost.

What the replaced mechanisms measured is in CLOSED: `ISSUES: the-strand-overdispersion-reconcile`,
`ISSUES: the-away-half-gdna-overdispersion`, `ISSUES: the-joint-strand-overdispersion-fit`.

OWNER RULINGS (2026-09-30):
- **κ** comes from a three-class binomial fit over the per-junction strand table: genuine (κ), splice artifact (½),
  reversed (1−κ).
  - The class weights are concave for fixed κ. κ is the shipped Beta(1,1) posterior mean over the reads the
    mixture's maximum credits to RNA, so it IS the shipped `(n_same + 1)/(n_obs + 2)` wherever the maximum puts every
    junction in the genuine class (owner, 2026-10-01).
  - The strand-live gate stays on the pooled counts. Removing artifacts only moves κ away from ½, so the fit cannot
    switch a live channel off.
- **An empty spliced census** is unstranded (κ = ½, od = 0), not an error.
- **Splice-artifact pruning** weighs a fragment's unspliced-gDNA reading against its spliced one by the artifact
  posterior, never a cut.
- **If an od is ever needed**, the owner's first choice is an RNA-only value from the most deeply sequenced junctions,
  applied as the shared value. It is developed on REAL data across several libraries, never one (an exception to
  `TRAPS: real-data-is-a-test-input`, scoped to the od, because the simulator plants no RNA od and no splice
  artifacts).

κ from the genuine junctions LANDED 2026-10-01, after od = 0 (the retired RNA od fit read κ, and at the clean κ LBX0588's
moment went 0.07 → 0.37). Its A/B is in `DESIGN.md` §3.3b and its fit in `EQUATIONS.md` §5.4.

THE REMAINING LANDING, one window each:
1. **The `summary.json` diagnostics** (`ISSUES: the-format-changes-to-batch-before-release` (a)):
   - pooled and clean κ;
   - artifact share by junction and by fragment (LBX0588: 22 % against 12 %);
   - reversed share, strand-live, and the od's label.
2. **The pruning design** with `ISSUES: splicing-artifacts`: junction ids in the fragment buffer (it carries none
   today), deposits deferred until the junction decision, and the weight derived on one page. Strand catches only
   artifacts with a wrong-strand read (an n-read artifact shows none with probability 2⁻ⁿ), and sees nothing on
   unstranded libraries.

DEFERRED PAST 0.8.0:
- an adaptive od whose evidence no single object can swing;
- the native strand variance's `(n·f)²`, which should be `n(n−1)·f²` (`transfer_rows.h` `strand_variance`; it matters
  only at od > 0).

Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-30_robust_od/` — the harnesses (`harness/` od, `harness_w2/` κ),
the per-object censuses, `runs/w1/VERDICT.md`, and the VCaP truth scorer `score_vcap.py`. The reviewer brief is the
sandbox's strand-model review note.

### splicing-artifacts
`priority: now — the cluster track, after the strand overdispersion (owner, 2026-09-30); the owner's highest real-data priority (2026-09-26), a bigger problem than previously thought (owner, 2026-09-24) · kind: defect · 2026-09-24`
A gDNA read that the aligner splices across an annotated junction reads as certified RNA unless the index's
blacklist — junctions the same aligner wrote on simulated genomic reads, kept at two or more, with the longest
anchor seen — rejects it. Today a fully rejected fragment is labelled artifact: the EM treats it as unspliced
(`DESIGN.md` §3.1d) but its blocks stay split at the rejected junction, so its footprint still spans the intron
and over-states gDNA's length when the short anchor is the outer block; calibration holds it out entirely. The
owner's first repair to consider: edit the alignment — delete the rejected junction and keep only the largest
aligned block — so the fragment is unspliced in both stages at its true length. Beyond that: the blacklist is a
yes/no cut (count ≥ 2, anchor ≤ the longest seen) on a continuous quantity, the aligner's rate of writing that
junction on genomic sequence at that anchor length. Measurable only on aligned data: the aligned ladder is drafted
for the cluster (`~/Downloads/rigel_runs/aligned_ladder/`), and on the production store 23.3 % of the ladder's
annotated junctions are blacklisted. Also here, from
`ISSUES: calibration-and-the-em-disagree-on-what-can-be-gdna` (closed): three fragment types calibration still puts on
its unspliced bank, gDNA-eligible, that the EM treats as RNA only — an unannotated sequenced junction; two annotated
junctions on opposite motif strands; mates disagreeing on an acceptor, which calibration merges into an unannotated
intron and the resolver calls annotated — and one length rule the stages apply differently: calibration drops an
unspliced or artifact reading past the maximum fragment length, the EM keeps gDNA while its likelihood competes.
Those three fragment types are a local A/B, independent of the cluster work.
- BOTH WAYS. On the VCaP RNA half the production blacklist rejected 598,010 records where a pure-gDNA library's rate
  predicts about 20,000 (uniquely mapped single-junction records 429,208 of 12.55 M, 3.4 %; multi-junction 154,604 of
  1.88 M, 8.2 %). A rejected genuine read becomes an unspliced, gDNA-eligible fragment whose footprint spans the
  intron, so the cost lands on the gDNA and synthetic-span pools of exactly the listed genes. A repair is judged on
  both errors together.
- THE LARGEST ARTIFACT CLASS IS DECIDED BY THE REFERENCE, through h, the reference mismatches between the spliced and
  unspliced placements. On the DNA library at annotated junctions with a 1–10 bp short side, 96 % of artifacts have
  h = 0 (5,092 of 5,297) and 81 % sit on the acceptor side (4,272), where STAR's sjdbScore breaks the tie toward
  splicing; the blacklist catches them only where alignable sampled the junction (4,643 rejected, 449 h = 0 records
  escaped). On the RNA library 40,557 rejected records are h = 0 and 162,213 need h ≥ 1. The hypothesis to derive on
  one page, the multi-junction case included: an ambiguous junction ADDS an unspliced reading weighted about
  (ε/3)^h, ε from the library, h from the reference (at index time for annotated junctions, at scan time for novel
  ones); the owner's "edit the alignment" is its gDNA limit at h = 0. Rigel reads NM only — never MD, the sequence,
  the qualities or the reference — so this needs new scanner and index infrastructure.
- THE LONG-ANCHOR ARTIFACTS ARE NOT h = 0 (short side ≥ 21 bp: 4,840 DNA records, 60 % uniquely mapped; 98 % of the RNA
  library's rejections, 235,406 of 240,520, uniquely mapped). The hypothesis is remote origin, from retrocopies and
  paralogs: alignable records an artifact from every alignment, secondary and supplementary included, so a parent
  gene's junction is listed with long anchors on both sides and, under the OR rule, rejects every single-junction read
  up to about 2m aligned bases. Measuring it needs an unspliced end-to-end realignment, alignable's per-record origins
  and a junction mappability table; the EM side needs `ISSUES: multimapper-intergenic-alignments`.
- THE BUILDER applies `min_count` (2, underived) per (chrom, intron, strand, read_length) BEFORE aggregating over read
  lengths, so a junction seen once at each of several read lengths never enters, and the shipped feather keeps no
  count, so no threshold can be decided at run time. Fixed with the catalogue rebuild.
- UNTESTED RULES, to pin with local falsification tests before any rule changes: the `<=`, the OR of the anchors, the
  inert stored 0, cumulative anchors across earlier introns, and the mixed fragment (one junction rejected, one kept,
  deposited as spliced with the rejected N as a gap); every fixture's anchors (500–10,000 bp) are longer than any
  read. Undecided: the reject runs per record before the mates are joined, so a junction rejected on one mate and
  kept on the other survives, and any surviving junction beats ARTIFACT, so one kept junction shields the rejected
  ones.
- THE TRAINING SETS. The spliced strand table that sets κ and supplies the od's RNA junction seeds holds artifacts
  whose sense fraction is ½ (the VCaP exome-DNA half: 3,223 spliced observations on 2,255 junctions, κ 0.5008, feeding
  the RNA od with information 5,701), and κ pools per-junction rates that vary with depth
  (`ISSUES: strand-overdispersion-one-shared-value`), so on real libraries the RNA strand MEAN is itself off:
  LBX0588's pooled κ is 0.064 against 0.0030 from its genuine junctions; κ now comes from the genuine junctions
  (2026-10-01, `DESIGN.md` §3.3b), and the per-junction artifact posterior it fits is this entry's input. The repair trains the strand model and
  the RNA length law weighted by each fragment's probability of being genuine.
- THE COUNTERS ARE IN THREE UNITS: `sj_blacklisted` counts junction × record over every record, multimappers
  included; the `splice.*` census counts uniquely mapped fragments; 0.7.1's `summary.json` splice counts are the
  length models' observations. So the 344-library cohort's 0.7.1 summaries do not compare with today's census.
- THE WORK runs on the cluster in phases: the substrate (per-record splice rows, a working index checksummed against
  production and built after the index format bump, `ISSUES: the-format-changes-to-batch-before-release` (b),
  `--notemp` re-runs, more truth mixes); both errors measured on today's tree; the mechanism census; the h-weighted
  reading (`ISSUES: scoring-penalties-are-underived-constants`); the training sets; aligner-settings robustness and the
  catalogue rebuild. The index case's splice evidence on the cluster's `/scratch` is gone (owner, 2026-09-28): it is
  regenerated from scratch by a new run, a separate task. The VCaP mix (its BAM, 0.7.1 `annotated.bam` and
  `summary.json`) and LBX0190 / LBX0588 / MO_3021 (the same, with `SJ.out.tab.gz`) may exist only in
  `~/Downloads/rigel_runs/cfrna/`: confirm cluster copies before the substrate phase and upload any that are missing.

### the-format-changes-to-batch-before-release
`priority: now — release-critical: the index bump before the cluster's working-index build, the summary.json schema while schema 3 is unreleased (owner, 2026-09-28) · kind: decision · 2026-09-28`
Changes to an output or on-disk format, batched so each format changes once. RULED (owner, 2026-09-28): the whole
batch, before 0.8.0; AMENDED the same day: the two density feathers are removed with the total-density landscape, not
renamed (`ISSUES: the-total-density-landscape`).
(a) `summary.json` (schema 3, unreleased):
- one od field instead of two, labelled binomial by policy (od = 0 since 2026-10-01), beside the clean κ and the
  artifact shares (`ISSUES: strand-overdispersion-one-shared-value`, the diagnostics step);
- `gdna_fraction` without the spliced intergenic fragments: `cli.py` adds `stats.n_intergenic`, which includes them,
  and a spliced fragment is certified RNA (Axiom 0). The count is so far inferred (27 at g001 OFF, 29 at g01 OFF, from
  the three-exon shadows), not read from `n_intergenic_spliced`; `quant_accuracy.py` copies it into the release
  report's gDNA pool;
- the depleted gDNA density, read off the gDNA landscape by `landscape.split_basins`, which `located_enriched_mode`
  calls and discards: it restores the report's on-target / off-target fold (dividing by `gdna_density_global` instead
  is biased, 25× against about 1,000× from the modes); it touches about six test construction sites and two
  instruments;
- the report's capture note, which says "no calibration" for a 0.7.1 or schema-2 summary that has a calibration block.
(b) The index (`INDEX_FORMAT_VERSION` 8): drop the written-never-read columns (`transcripts.feather`'s `abundance`,
`nrna_abundance` and `n_exons`; `sj.feather`'s `interval_type`), and add the duplicate map as an alias map
`dropped_t_id → kept_t_id`. Every index is rebuilt, the cluster's included, so land it before the cluster's
working-index build and the cluster rebuilds once. `ISSUES: overlapping-synthetic-shadows`' index merge is a separate
decision.

### an-unrecorded-splice-blacklist-is-dropped-silently
`priority: now — the cluster track: check the manifest before the next cluster quant (owner, 2026-09-28) · kind: risk · 2026-09-28`
`TranscriptIndex.load` applies `splice_blacklist.feather` only when `manifest.json` records `sources.alignable_zarr`;
otherwise it sets `sj_blacklist_size = 0` with no log line. The cluster's production index
(`hulkrna/refs/human/rigel_index/`) holds a 5.8 M-row feather: if it was placed by hand, the next quant there runs
with splice-artifact detection off and no error — format 8 (eeba4b0a) predates the sources block (8c673324), with no
version bump. The only trace is `sj_blacklist_loaded: false` in `summary.json`. Before any cluster quant, check the
manifest or rebuild with `--alignable-zarr`. The one-line warning for a feather without a recorded store is in
`ISSUES: latent-defects`.
Instrument: the manifest; `summary.json`'s `sj_blacklist_loaded`.

### calibration-detects-capture-on-a-capture-off-library
`priority: now — release-critical: the uncommitted count repair is under review for law-frame and uncertainty defects; the detector-free reader landed 2026-10-10 (`DESIGN.md` §7.2) and this row's number under it is unmeasured · kind: defect · 2026-09-29; mechanism 2026-10-01`
THE SYMPTOM (re-recorded 2026-10-03 on the re-simulated panel, pinned and fractional; the first record, 2026-10-01, was
on the deleted old-physics panel and read 42.1 / 3.8 %, genes 4.2 / 0.3 %, 102k). On the genome-scale fl-gap arm
`flgap_rna_short` (realised RNA 78 bp, gDNA 250 bp, 100 bp reads, `g50`), the stranded × capture-OFF library reads
transcripts 42.3 % against 3.8 % unstranded, genes 3.3 % against 0.3 %, and 80k true mRNA fragments land on nascent RNA
(false-positive mass 350k against 23k unstranded), while its calibration metric is unremarkable (region mwae 0.011,
library gDNA 5.01 M against a truth of 5.00 M). That aggregate metric hid the per-object count errors
established below; the reference magnifies their effect in the EM.

THE OBSERVED FAILURE PATH (2026-10-01; its upstream cause is corrected by the 2026-10-05 map experiments below,
`~/Downloads/rigel_runs/prototypes/2026-10-01_shortrna/VERDICT.md`):
1. In a stranded library, calibration's strand deconvolution credits a short exon (about 100 bp) with a fraction of a
   gDNA fragment: its RNA's few wrong-strand reads.
2. That exon has almost no gDNA opportunity. A 250 bp fragment can rarely be contained in it: the slots that make the
   false mode have a median of 0.16 bp of admissible starts, and 91 % have less than one position (325 of the 358 counted
   slots within 3× of the reference at the last refit; none has ten).
3. So its gDNA density is about 10 per bp, 200× the library's 0.05.
4. These slots pass every training rule of `DESIGN.md` §7.1. Rule 4 asks whether the COMPOSITION is located
   (`Var(log f_g) ≤ 1 nat²`), never whether the DENSITY is.
5. They form a located enriched mode: 10.07 per bp, 393 members (the deleted panel read 9.65 per bp, 370).
6. §7.2's ruler reads that mode as the capture reference, so on a capture-OFF library the regions' capture efficiencies
   read 3.7e-6 to 1, median 0.0054 (10th–90th percentile 0.0048–0.020); every other capture-OFF row of both gap arms
   reads no reference at all.
7. Every EM length contracts unevenly, by the efficiencies of the objects each transcript spans, so isoforms of one gene
   are contracted differently and the most contracted absorb the shared fragments; the gDNA component contracts with
   them (on the deleted panel: transcripts by 0.0017–0.10, a 20-fold spread within a gene, gDNA about 200×).

Why the failure has its shape:
- Only stranded: unstranded, the exons' compositions are not located.
- Only at depth: at 10 % no mode forms and both libraries read 9.3 %.
- Only RNA shorter than gDNA: on the RNA-long arm gDNA is the short component.

ISOLATED, one variable each, on the failing condition (2026-10-01, on the deleted panel; the shipped row is re-recorded
on the new panel at 42.3 % and the no-reference arm at 3.53 %, see THE REFERENCE-SOURCE A/B below):
- shipped 42.1 %;
- the EM's strand term at ½ 42.8 %;
- a uniform warm start 42.1 %;
- MAP 42.1 %;
- calibration and the EM unstranded 3.6 %;
- no capture reference and nothing else 3.5 %.

Not od = 0 or κ: the tree before both (35e254ef) gives the shipped row to the fragment.

THE PROTOTYPE (outside the tree): a counted slot trains the gDNA landscape only where its gDNA opportunity admits at
least one contained fragment position (`eff ≥ 1`). Below one position the slot is shorter than gDNA's fragments and has
no measurable gDNA density (`TRAPS: density-below-one-fragment-length`); the anchors are untouched. The true
references keep their members (ladder g05 ON 4,214 of 4,528 near-reference slots, g50 ON 4,474 of 4,864, test g05 ON
352 of 352), while the false mode keeps 29 of 376.

The A/B against the same session's baselines (2026-10-01; the gap-arm row on the deleted panel, not re-recorded since
the guard is refused):
- `flgap_rna_short` stranded × OFF: transcripts 42.08 → 3.54 %, genes 4.20 → 0.23 %, mRNA error 102,416 → 920, nascent
  95,350 → 3,144.
- The ladder: flat in every stratum. Its zero controls fall slightly: false gDNA 126 / 335 / 127 / 183 → 120 / 334 /
  120 / 176. The test panel's are unchanged. That meets `ISSUES: gdna-landscape-trains-on-false-positives`' constraint.
- odg05 and the four other fl panels: transcripts within ±0.15 pt.
- The cost: the test panel's stranded × ON gDNA pool, +10.7 % (6.1k). It is almost all one row, g98 ss0.70 ON (−19.7k →
  −25.4k of 1.15M); that row's reference is unchanged and its RNA is 2 % of the library. The stratum's transcripts and
  mRNA improve.

THE GUARD IS A BAND-AID (owner, 2026-10-02). gDNA's length is a distribution, so a slot with one or two positions is
admitted while its density is still strand noise over a small opportunity. With the guard on, calibration still
over-attributes gDNA 7× below one position and 3× at one to three; truth 244 against 1,118.

The cause is wider than training, but changing the coordinate alone is not a fix: for fixed total and opportunity,
`Var(log rho_g) = Var(log f_g)`. The measured defect is transporting COUNT odds as if the component opportunities
were equal. On the ladder the mean `|log(E_g/E_r)|` is 0.004; on the gap arms it reaches 6 nats. The component
opportunities must enter each message map, as the experiments below establish. The separate capture × length
half remains in `ISSUES: the-scorer-reads-a-census-length-law`; these results do not establish a new training rule
or a spike-and-slab background likelihood as necessary.

The guard's one cost, g98 ss0.70 ON, is not the length mechanism. The guard drops one slot. Three both-stranded objects
whose composition is the median of a multimodal, prior-decided posterior flip modes on that one kernel: 0.30 → 0.035 per
bp against a truth of 0.31. A random located slot dropped instead moves calibration by ±30 fragments.

THE REFERENCE-SOURCE A/B (2026-10-03, outside the tree, pinned and fractional, every panel, every stratum; the
prototype's `VERDICT.md`). Part 2's evidence-fitted gDNA population (one shared RNA amount per region, the gDNA-only and
one-RNA-strand classes) read through the shipped mode reader as the capture reference, nothing else moved: it repairs
this row (transcripts 42.30 → 3.53 %, genes 3.28 → 0.23 %, mRNA error 79.7k → 0.5k) by reading no enriched mode on a
capture-OFF library, is bit-identical to shipped on every other capture-OFF row and every g00 row, moves stranded ×
capture-ON a little and both ways (ladder 6.14 → 6.00 %, RNA-long 6.17 → 5.77 %, test chromosome 13.53 → 14.58 %; its
reference 12–19 % below the shipped one), and cannot serve an UNSTRANDED captured library: with no strand channel the
fit has no evidence class beyond the RNA-free regions, declares no capture where capture is real, and the EM's gDNA pool
drains (ladder deferred row 12.2 → 52.3 %, test 18.1 → 44.1 %). The count-only class (Part 2's step two) grows false
modes instead (0.0014/bp on a captured unstranded row, 0.058/bp on a capture-OFF one). Not landed. The open ruling: an
evidence-fit reference where the library is stranded and the shipped landscape otherwise, or another statement of where
each source is trusted (owner).
Instrument: `~/Downloads/rigel_runs/prototypes/2026-10-03_fl_arms/` (`harness/sitecustomize.py`, `RIGEL_ARM=evidence_ref`;
`VERDICT.md`); the 2026-10-01 prototype's `quant_sub.py`, `training_census.py` (who trains the landscape and at what
opportunity) and `proto/compare.py`.
THE PER-OBJECT COUNTS, MEASURED (2026-10-05; the owner's spectrum ruling forbids the `None` repair). Calibration's
gDNA count is grossly wrong (more than 3-fold and 5 Poisson sd off the oracle's `slot_truth.npz`) at hundreds to
thousands of objects on both gap arms, while the library totals are right: RNA-short stranded × OFF 173 regions
over-called (median 7.3×; the short exons above) and 156 under-called (−16k fragments, intronic clusters booked as
RNA); RNA-long unstranded × OFF 136 regions and 1,431 boundaries over-called (+31k / +55k; whole RNA-rich
neighbourhoods booked as gDNA, the exon and both its boundaries together, so coherent). The shipped ruler's `None`
hides every one of them on an uncaptured library. A ruler with no detector inherits them: a prototype with no
reference (an opportunity-weighted population read by its posterior median, a floor at the off-target level, local
coherence, relative weights) reads RNA-short stranded × OFF at 3.55 % (genes 0.23) and beats 0.7.1 on every test
chromosome stratum under three probe layouts, and takes RNA-long unstranded × OFF from 2.13 to 11.71 % (0.7.1: 5.29).
So the fix that generalizes is calibration's: each object's gDNA in the density frame, so that an object gDNA cannot sit
in cannot hold gDNA and one RNA cannot sit in cannot hide it, with this table as its gate.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-10-05_spectrum_ruler/` (`dump.py` + `analyze_rs.py` read any
condition's per-object counts against `slot_truth.npz` in seconds; `harness.py`, `arms.py` arm `pwfx_own_med`; README).

INDEPENDENT CENSUS AND CHANNEL REPLAY (2026-10-05, `a7b03103`; fresh calibrations of validated cached scans,
default calibration threads, 3.4–22.5 s per row). Every calibration/oracle per-object total agrees within
`7.3e-12`. The reporting predicate is explicit: `abs(k-g) > 5 sqrt(max(g,1))`, and
`k/max(g,1) > 3` or `< 1/3`. This is a diagnostic cutoff, not a calibrated uncertainty test for deconvolution.

| capture-OFF row | regions over: objects / excess | regions under: objects / deficit | boundaries over: objects / excess |
|---|---|---|---|
| RNA-short · stranded | 173 / 2,351 | 156 / 16,115 | 1 / 714 |
| RNA-short · unstranded | 34 / 941 | 205 / 23,010 | 12 / 120 |
| RNA-long · stranded | 17 / 433 | 0 / 0 | 664 / 7,750 |
| RNA-long · unstranded | 136 / 30,594 | 0 / 0 | 1,431 / 55,105 |

The original census reproduces exactly. Its "short-exon" label is not a complete classification: the
173 flagged regions have median opportunity 1.387, but include larger and both-strand objects; one region
alone contributes 911 of the 2,351 excess fragments. Inspect signatures/opportunities and zero-truth
objects separately rather than requiring this entire mixed set to disappear through one mechanism.

Holding the final landscape, strand model, opportunities and lattice FIXED, a local solve without messages
restores the two principal groups:

| original flagged group | oracle gDNA | full solve | local evidence + same factory/landscape |
|---|---:|---:|---:|
| RNA-short stranded · 156 under-called regions | 16,878 | 763 | 16,575 |
| RNA-long unstranded · 136 over-called regions | 3,146 | 33,740 | 4,821 |
| RNA-long unstranded · 1,431 over-called boundaries | 3,897 | 59,002 | 6,108 |

Six representative RNA-short introns and seven RNA-long regions/boundaries were replayed in complete
locus blocks with the global reductions retained. Block outputs reproduce production bit-identically;
reconstructing the final solve from the captured received tables agrees within `1e-14` in fraction.
Removing the responsible **composition** row restores the local answer; removing either RNA level lane
or the gDNA level lane changes none of these examples. RNA-short slot 9604: truth 407, full 0.282,
without left composition 394.104, without right composition 394.458. RNA-long slot 5071: truth 8, full
1,625.601, without left composition 4.217. RNA-long slot 33731 has harmful composition on BOTH sides:
one-sided removals alone retain a large over-call, so an ablation must include joint removal.

This refutes the RNA-level-lane explanation for the inspected introns and does not implicate the
duplicated prior as the necessary cause of these groups. Four independent counterexamples then establish
the composition maps' opportunity defects, under **Density is the frame-invariant currency; a fraction
is not** (`EQUATIONS.md`). Each uses analytic contained/crossing opportunities for a known density mixture,
both length-gap directions and an equal-length control, with and without certified splice flux:

1. The reverse licensed splice map substituted gDNA opportunity for RNA's and `S/Eg` for the actual
   route rate, although the forward map already used the correct quantities. Correcting the reverse
   map removes all 156 RNA-short stranded under-calls, but exposes 105 under-calls on RNA-long stranded.
2. The intron–boundary face forwarded count odds unchanged. Shared densities instead require
   `lambda_dst-lambda_src = log(Eg_dst/Eg_src)-log(Er_dst/Er_src)`. Correcting both directions reduces
   RNA-long unstranded boundary over-calls from 1,155 to 591, but is insufficient alone.
3. Both alternative-splice-site flank maps also substituted gDNA opportunity for RNA's. Their
   disagreement-width center compared count odds instead of the map's prediction. Correcting these
   removes the remaining RNA-short under-calls AND the false enriched landscape mode, with no training
   change. RNA-long unstranded region excess falls from 33,619 to 1,444 fragments.
4. The outside flank at an exon–exon terminus made the same substitution. Correcting it removes
   the remaining 26 RNA-long stranded under-calls and reduces RNA-long unstranded boundary excess
   from 10,142 to 117 fragments. The inside level rule is unchanged.

Each cumulative arm was freshly calibrated on all four OFF rows before proceeding; no transcript A/B
was spent on the intermediate arms that still failed the count screen. Final prototype census:

| capture-OFF row | regions over: objects / excess | regions under: objects / deficit | boundaries over: objects / excess |
|---|---|---|---|
| RNA-short · stranded | 20 / 198 | 0 / 0 | 0 / 0 |
| RNA-short · unstranded | 34 / 942 | 0 / 0 | 13 / 141 |
| RNA-long · stranded | 17 / 196 | 0 / 0 | 3 / 28 |
| RNA-long · unstranded | 42 / 1,444 | 0 / 0 | 10 / 117 |

There are no flagged boundary under-calls. Residual zero-truth synthetic objects remain; this is not
a claim that every inferred count is exact. The capture-ON screen improves absolute gDNA-count error
in every object class of both stranded rows; the deferred unstranded RNA-long row is mixed. The
four map families' falsification gates fail before each edit (4, 4, 18 and 8 cases respectively),
and all 48 pass afterward. Restoring the wrong opportunity, rate or identity in the prepared native
tables makes the corresponding gates fail again. No landscape weight, prior, capture detector,
training guard or new tunable was added. After the full panel gate, these four map repairs were
copied from the isolated native worktree into the working tree on 2026-10-06, with their derivation
and tests; the owner still drives the commit. The ruler is unchanged.

The first pinned/fractional transcript A/B changes those maps **and the alternative-splice discrepancy
center**, with calibration and EM threads at their defaults and four condition shards. The latter also
changes widths at equal component opportunities and was not A/B'd separately. The test chromosome beats 0.7.1 in all twelve strata
(the three strand specificities × capture labels, with zero-gDNA controls separate); it still loses
against the current tree on several stranded rows. The worst per-condition change there is g98
ss0.70 ON: transcripts 46.92 → 69.50 %, genes 19.51 → 42.44 %, against release 55.98 / 39.60 %.
That condition must not disappear behind its stratum average. The complete RNA-short A/B is:

| transcripts / genes (%) | 0.7.1 | current `a7b03103` | component-opportunity maps |
|---|---:|---:|---:|
| stranded × OFF | 6.31 / 0.44 | 42.30 / 3.28 | 3.53 / 0.23 |
| stranded × ON | 20.33 / 3.52 | 19.97 / 1.53 | 13.33 / 1.40 |
| unstranded × OFF | 6.00 / 0.61 | 3.78 / 0.27 | 3.79 / 0.27 |
| unstranded × ON | 29.76 / 6.19 | 19.41 / 2.34 | 17.57 / 2.23 |

The captured stranded improvement weakens the supposed ~20 % capture-by-length floor: this map-only
repair already reaches 13.33 %. The opposite gap passes too:

| RNA-long: transcripts / genes (%) | 0.7.1 | current `a7b03103` | component-opportunity maps |
|---|---:|---:|---:|
| stranded × OFF | 5.23 / 0.46 | 2.26 / 0.19 | 2.27 / 0.19 |
| stranded × ON | 7.38 / 1.62 | 4.48 / 0.68 | 4.23 / 0.66 |
| unstranded × OFF | 5.29 / 0.59 | 2.13 / 0.22 | 2.16 / 0.22 |
| unstranded × ON | 28.99 / 5.84 | 6.86 / 0.95 | 6.55 / 0.94 |

The g98 ss0.70 ON test-chromosome regression is a refitted-landscape effect, not an erroneous final
map at those objects. On slots 2962 / 3010 / 2986, oracle gDNA is 3,137 / 2,904 / 3,034; estimates
move 2,949 / 2,438 / 2,989 → 341 / 319 / 2,241. Crossing the two frozen final-sweep inputs and
their two landscapes, under BOTH native implementations, reproduces those numbers exactly according
to the landscape alone. The reference is unchanged at 3.196705. The earlier claim that these exons
lack local evidence is withdrawn: the archived strand-only replay gives 3,106 / 2,822 at the first
two slots, near truth 3,137 / 2,904; strand plus landscape gives about 30 / 30. The landscape's mass
in the bands below 0.01, 0.01–0.3 and above 0.3/bp agrees to three decimals (0.696 / 0.027 / 0.277),
but its relative log density at 0.3/bp deepens from −8.59 to −9.92. Local evidence exists and is
overridden. ψ uses an interpolated histogram median: low density between modes can make it highly
sensitive, but is not by itself proof of a mathematical discontinuity. No perturbation sweep has
yet established the size of this sensitivity for this repair. This residual remains visible,
rather than being called a uniform count improvement. The completed ladder and alternative probe layouts are below. Every
transcript stratum beats 0.7.1. One secondary release deficit already present on the current tree
remains: the ladder's deferred unstranded capture-ON zero control has gene error 3.218 % against
0.7.1's 3.183 %. The map repair does not cure it. These are Phase 2 measurements, not a completed
detector-free release A/B.

**Review of the count-map repair, 2026-10-06.** The algebra is conditional on the opportunities it
receives. `pipeline.run_pipeline` passes the uniform-frame gDNA law and the capture-selected spliced RNA
law into calibration; the common geometric operator does not make those laws commensurate. The
frozen equal-chemistry g98 ss0.70 ON input has boundary `log(Eg/Er) = −0.04425` and a maximum region
value of `1.19842`, identical in both builds. Calibration already used these inputs in other lanes;
the repair makes additional composition faces read them. All six contaminated ladder ON rows move
gDNA downward and the nRNA pool upward; all six OFF rows move oppositely. The signed association is
verified, but the law mismatch's causal share is not isolated. A calibration-only law override also
changes RNA level lanes and junction rates, so its result alone cannot locate the cause in the maps.
The length-frame problem's home is `ISSUES: the-scorer-reads-a-census-length-law`.

The alternative-splice discrepancy-center correction is a second mechanism: even at equal
opportunities `f_b=1/2, U=S=100` predicts odds `1/3`, against `1/2` from the old center. It requires
its own A/B. The retained width terms also need review. For independent Poisson junction counts,
the delta-method route-rate log variance is `Σ J/A² / (Σ J/A)²`, not generally `1/Σ J`; it omits
covariances when junction observations share fragments. At zero certified rate, TRANSPORT still
blurs with `trigamma(1/2)` while SPLICE_OUT's rate marginal does not. These are not validated as
one common uncertain factor. Required gates include unequal junction opportunities, mirrored
route-side builders, shallow-count uncertainty, and the zero-opportunity switch from composition
to the gDNA level lane. Existing prepared-table perturbations do not cover those builder decisions.

Small transcript-only OFF differences remain measured differences, with attribution unresolved.
At test g00 ss0.99 OFF, transcripts move 6.99834 → 7.70772%, genes stay 0.83498352%, and the RNA pool
changes by less than 0.01 fragment. This is consistent with isoform-allocation amplification
(`ISSUES: benchmark-noise-floors-unmeasured`), but unchanged totals do not prove that a particular
isoform loss is noise. Neither dismiss the loss nor rank the calibration mechanism on it alone.

Equal-length ladder: transcripts / genes (%), zero-gDNA controls (`g00`) read separately.

| Stratum | 0.7.1 | current `a7b03103` | component-opportunity maps |
|---|---:|---:|---:|
| stranded × OFF | 3.95 / 0.44 | 2.02 / 0.25 | 2.02 / 0.25 |
| stranded × OFF · g00 | 3.70 / 0.47 | 1.93 / 0.13 | 1.93 / 0.14 |
| stranded × ON | 9.91 / 2.52 | 3.85 / 1.38 | 3.88 / 1.38 |
| stranded × ON · g00 | 20.22 / 3.47 | 6.46 / 2.36 | 6.47 / 2.36 |
| unstranded × OFF | 10.29 / 1.49 | 2.29 / 0.31 | 2.29 / 0.31 |
| unstranded × OFF · g00 | 20.09 / 2.53 | 1.60 / 0.14 | 1.60 / 0.14 |
| unstranded × ON | 33.21 / 14.25 | 10.42 / 4.33 | 12.40 / 4.47 |
| unstranded × ON · g00 | 10.22 / 3.18 | 7.27 / 3.22 | 7.27 / 3.22 |

Test chromosome, standard probes: transcripts / genes (%), zero-gDNA controls (`g00`) read separately.

| Stratum | 0.7.1 | current `a7b03103` | component-opportunity maps |
|---|---:|---:|---:|
| ss 0.70 × OFF | 9.44 / 1.54 | 8.61 / 1.25 | 8.62 / 1.25 |
| ss 0.70 × OFF · g00 | 10.06 / 1.47 | 6.61 / 0.90 | 6.59 / 0.90 |
| ss 0.70 × ON | 10.28 / 3.32 | 6.91 / 1.51 | 7.69 / 1.69 |
| ss 0.70 × ON · g00 | 13.51 / 6.20 | 8.27 / 2.16 | 8.27 / 2.16 |
| stranded × OFF | 8.20 / 1.20 | 7.39 / 1.08 | 7.67 / 1.08 |
| stranded × OFF · g00 | 11.35 / 1.26 | 7.00 / 0.83 | 7.71 / 0.83 |
| stranded × ON | 8.58 / 2.45 | 7.19 / 0.96 | 7.58 / 0.96 |
| stranded × ON · g00 | 14.22 / 5.83 | 7.31 / 0.67 | 7.35 / 0.67 |
| unstranded × OFF | 9.25 / 1.76 | 8.78 / 1.33 | 8.66 / 1.33 |
| unstranded × OFF · g00 | 9.38 / 1.49 | 5.71 / 1.01 | 5.71 / 1.01 |
| unstranded × ON | 18.94 / 10.31 | 16.07 / 5.22 | 15.97 / 5.21 |
| unstranded × ON · g00 | 17.93 / 7.39 | 9.48 / 2.44 | 9.23 / 2.44 |

Test chromosome, junction-spanning probes: transcripts / genes (%), zero-gDNA controls (`g00`) read separately.

| Stratum | 0.7.1 | current `a7b03103` | component-opportunity maps |
|---|---:|---:|---:|
| ss 0.70 × OFF | 9.44 / 1.54 | 8.61 / 1.25 | 8.62 / 1.25 |
| ss 0.70 × OFF · g00 | 10.06 / 1.47 | 6.61 / 0.90 | 6.59 / 0.90 |
| ss 0.70 × ON | 40.88 / 1.70 | 32.41 / 1.59 | 32.44 / 1.59 |
| ss 0.70 × ON · g00 | 39.17 / 0.77 | 39.02 / 0.40 | 39.05 / 0.40 |
| stranded × OFF | 8.20 / 1.20 | 7.39 / 1.08 | 7.67 / 1.08 |
| stranded × OFF · g00 | 11.35 / 1.26 | 7.00 / 0.83 | 7.71 / 0.83 |
| stranded × ON | 40.01 / 1.46 | 32.16 / 1.15 | 32.06 / 1.14 |
| stranded × ON · g00 | 38.99 / 0.41 | 38.75 / 0.09 | 38.75 / 0.09 |
| unstranded × OFF | 9.25 / 1.76 | 8.78 / 1.33 | 8.66 / 1.33 |
| unstranded × OFF · g00 | 9.38 / 1.49 | 5.71 / 1.01 | 5.71 / 1.01 |
| unstranded × ON | 42.33 / 3.28 | 30.26 / 2.23 | 29.20 / 2.23 |
| unstranded × ON · g00 | 40.58 / 1.35 | 40.40 / 0.96 | 40.40 / 0.96 |

Test chromosome, sparse probes: transcripts / genes (%), zero-gDNA controls (`g00`) read separately.

| Stratum | 0.7.1 | current `a7b03103` | component-opportunity maps |
|---|---:|---:|---:|
| ss 0.70 × OFF | 9.44 / 1.54 | 8.61 / 1.25 | 8.62 / 1.25 |
| ss 0.70 × OFF · g00 | 10.06 / 1.47 | 6.61 / 0.90 | 6.59 / 0.90 |
| ss 0.70 × ON | 55.26 / 29.42 | 27.47 / 4.95 | 28.17 / 4.94 |
| ss 0.70 × ON · g00 | 75.66 / 64.73 | 36.58 / 2.37 | 36.70 / 2.37 |
| stranded × OFF | 8.20 / 1.20 | 7.39 / 1.08 | 7.67 / 1.08 |
| stranded × OFF · g00 | 11.35 / 1.26 | 7.00 / 0.83 | 7.71 / 0.83 |
| stranded × ON | 74.76 / 63.06 | 24.37 / 3.15 | 24.44 / 3.14 |
| stranded × ON · g00 | 45.39 / 24.44 | 35.37 / 0.72 | 35.37 / 0.72 |
| unstranded × OFF | 9.25 / 1.76 | 8.78 / 1.33 | 8.66 / 1.33 |
| unstranded × OFF · g00 | 9.38 / 1.49 | 5.71 / 1.01 | 5.71 / 1.01 |
| unstranded × ON | 95.00 / 66.90 | 40.11 / 13.50 | 39.91 / 13.32 |
| unstranded × ON · g00 | 76.62 / 64.96 | 37.72 / 2.65 | 37.72 / 2.65 |

Goldens were inspected before their isolated-worktree update: six scenarios change,
maximum transcript-count change 0.175 fragments; the extreme mixed toy's BAM-certified locus gDNA is
227 and the estimate moves 192.631 → 204.506 (closer), while the moderate toy's truth is 189 and the
estimate moves 189.138 → 192.609 (farther). This is an accuracy change, not a numerical no-op.
The re-derived suite is 2,699 passed, zero skipped (2,651 baseline + 48 new cases); ruff and
`preflight.py --full` pass. An earlier suite run lacking the conda tools on PATH was discarded;
the final run uses the complete environment and arms every child with the isolated native build.
The landed source was then rebuilt in the conda environment and independently passed all 2,699
tests (zero skipped), ruff and full preflight on 2026-10-06. Its one quadrature warning is unchanged.

The broader claim that no detector-free ruler can tolerate coherent count errors remains a hypothesis;
these measurements justify correcting the counts first, not an impossibility theorem about all rulers.

RULER SCREEN AFTER THE MAP REPAIR (2026-10-06; prototypes only). The earlier posterior-median
readout fails the spectrum ruling: on two population atoms at density 1 and 100 with equal posterior
mass, perturbing one count by ±1e-3, ±1e-6 or ±1e-9 changes the selected density 100-fold.
The geometric posterior mean is continuous through those ties. Its three continuity gates pass,
and deliberately restoring the median fails all three again. This is continuity of the readout
conditional on the population, not a proof about the entire fitted pipeline.

Both readouts were tested on the repaired counts, retaining the previous prototype's population
fit, off-target floor, spatial treatment and junction rule. They fail the unstranded capture-OFF
test-chromosome stratum. A third arm changes ONLY the ruler's working count likelihood: an inferred
count `k` with deconvolution log-fraction variance `V` has delta-method variance `d=k²V`, so the
working variance `mu+d` gives quasi-score `(k-mu)/(mu+d)`. This is a moment approximation, not an
independent evidence likelihood. It restores OFF but fails ON:

| test-chromosome stratum; transcripts / genes (%) | 0.7.1 | current `a7b03103` | map repair | median ruler | geometric ruler | geometric + count uncertainty |
|---|---:|---:|---:|---:|---:|---:|
| unstranded × OFF | 9.25 / 1.76 | 8.78 / 1.33 | 8.66 / 1.33 | 9.41 / 1.41 | 9.41 / 1.40 | 8.81 / 1.33 |
| ss 0.70 × ON | 10.28 / 3.32 | 6.91 / 1.51 | 7.69 / 1.69 | 9.34 / 2.19 | 8.13 / 2.15 | 18.73 / 3.91 |
| stranded × ON | 8.58 / 2.45 | 7.19 / 0.96 | 7.58 / 0.96 | 7.08 / 1.16 | 7.11 / 1.10 | 8.18 / 1.95 |
| unstranded × ON (deferred) | 18.94 / 10.31 | 16.07 / 5.22 | 15.97 / 5.21 | 15.78 / 5.48 | 15.87 / 5.52 | 37.04 / 8.29 |

Each arm ran all 30 test conditions, including every zero control. Failed arms were stopped before
the gap arms, ladder or real libraries. The uncertainty arm's quasi-score identity tests fail in
three cases under Poisson, pass all five cases with the new formula, and fail the same three when
the uncertainty term is deliberately removed. Those algebraic checks do not validate the moment
approximation: the captured-row A/B falsifies its adequacy here.

The OFF failure is not explained solely by the intergenic floor. On g25 / g50 unstranded OFF,
substitute oracle counts ONLY into the ruler, keeping calibration and the EM priors fixed:

| counts read by the geometric ruler; transcripts / genes (%) | g25 | g50 |
|---|---:|---:|
| native inferred counts | 9.26 / 1.16 | 11.29 / 1.83 |
| oracle intergenic / gene-edge counts only | 9.13 / 1.16 | 10.72 / 1.84 |
| oracle counts at the other objects only | 8.11 / 1.16 | 9.86 / 1.51 |
| oracle counts at every object | 8.08 / 1.16 | 9.77 / 1.53 |

At g50, the `gB4_capcluster_ba` exon (slot 2074) has truth 55 gDNA fragments and estimate 95.803;
its reported `Var(log f_g)=0.341` implies a deconvolution sd about 56, but a Poisson read of that
estimate uses sd about 10. The ruler prices it 1.878× background and magnifies isoform error.
Slot 3006 has truth 96, estimate 179.464 and variance 0.140; it is priced 1.613× background.
Removing all messages with the final landscape fixed gives 50.085 / 113.531. These residuals
are not a new proved opportunity-map error. Nor does their posterior variance constitute a
prior-free likelihood: the fitted landscape and neighbouring messages helped produce it.

The count-map repair therefore does not complete the proposed count-as-Poisson ruler. The next
design decision concerns the evidence interface and a single use of the density prior, under
`ISSUES: the-gdna-prior-enters-psi-twice`; no new partner or ruler has been landed. A continuous
per-object evidence model remains to be derived and approved, without an unconditional class split,
a capture detector or a strand-only reference. Real-library release validation waits for a passing
candidate. The capture reference / clip / `None` mechanism is still present in the working tree.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-10-06_rna_short_counts/` (`census.py`, `replay.py`, `messages.py`,
fresh NPZ/JSON outputs and serialized native inputs). No full-genome debug capture was needed.

EVIDENCE CHECKPOINT (2026-10-07; no production change). Five fresh cached calibrations of the
count-repair tree reproduce the three earlier test-chromosome count/opportunity/belief arrays
bit-for-bit. Object totals equal the certified oracle on all five conditions. The two gap
calibrations take 2.48 / 7.63 seconds; test calibrations take 0.27–0.30 seconds, excluding imports.
The calibration binary is `_solve_impl`, SHA256
`d78ef4c95f7ba1dcd3d48ad7ffd41532233aa7718a3bb7238370de9fc8a508f0`.
Earlier receipts hashing `_em_impl` alone do not identify the calibration binary.

A full count × opportunity substitution holds the background, object count, library ceiling
and CI readout fixed. Counts are certified realized origins; opportunities are simulator
pre-capture geometry, independently checked by enumerating fragment placements. These are
diagnostic reader weights, not a new end-to-end candidate:

| condition | maximum weight, inferred counts / fitted opportunity | true counts / fitted opportunity | inferred counts / true opportunity | true counts / true opportunity |
|---|---:|---:|---:|---:|
| test g00 ss0.99 ON | 23,831.03 | 1.0000 | 22,903.37 | 1.0000 |
| RNA-short g50 stranded OFF | 1.3347 | 1.0177 | 1.3348 | 1.0170 |
| RNA-long g50 unstranded OFF | 6.6908 | 1.0131 | 6.6616 | 1.0124 |

At test g25 ss0.70 ON, the 206 no-RNA boundaries with expected enrichment above 20 have exact
counts in every count arm. Their median log-weight error is −0.05038 with fitted opportunity,
−0.02212 with true opportunity, independently of the count substitution. Counts explain the
large zero-control/OFF outliers; opportunity error explains part of this ON bias. Neither
result justifies a universal exon correction. Full tables include every eligible object by
class and admitted strand, with no selection on the inferred count. Region and boundary
incidences are not pooled into a conserved fragment total.

For the worst zero-DNA boundary, the actual strand columns are 52 / 4,878, with fitted
positive-column RNA probability 0.0095571. Marginalizing the unknown RNA amount with the
existing RNA reference leaves zero DNA possible: expected DNA 18.5657 versus zero changes
log evidence by only +0.08448. The CI reader's Poisson likelihood on the inferred 18.5657
fragments assigns zero DNA exactly zero support. A falsification test fails on that reader
and passes on the exact observation-model reference. This establishes lost uncertainty,
not that a positive count posterior median is itself a bug. The boundary is probe-enriched
in the simulator; all-zero oracle DNA counts cannot recover its capture from DNA observations
alone. Do not interpret the neutral oracle-count arm as perfect expected-yield recovery.

The reference covers pure DNA, either RNA strand and both strands with the existing
total-RNA/tilt measure. Its 113 tests pass against independent integration; six deliberate
defects each fail their intended gates. It is an own-observation primitive only, not a
strand-only release reader or a certified local-message likelihood.
Instrument: `.cache/rigel_runs/2026-10-07_density_checkpoint/` (`freeze_inputs.py`,
`frozen_receipt.json`, `factorial.py/json`, `false_boundary_evidence.json`,
`density_evidence.py`, `test_density_evidence.py`, `mutations.json`, `test_oracle_geometry.py`).
The first gap comparison was stopped when a scratch-loop repeatedly decompressed its NPZ
arrays; materializing them once completed the two full comparisons in 49 / 44 seconds.

### the-gdna-prior-enters-psi-twice
`priority: next — the fragment-length campaign's ψ half: a correct fix measured, unlandable alone; its partner is an open ruling (owner) · kind: defect · 2026-10-02`
ψ holds one place for the gDNA rate prior — the reference's gDNA half `½·log f` without a fitted landscape, the landscape
with one — and the kernel adds both (`psi_kernel.h`, the arm), tilting gDNA by `+½` per nat of `log ρ` wherever a
landscape is fitted. THE FIX (Part 1, prototyped in C++ in the worktree `~/proj/rigel-part1`, the patch in
`~/Downloads/rigel_runs/prototypes/2026-10-02_part1/`): the landscape REPLACES the gDNA half, read on the unclipped
`log σ(λ)`, continued below the grid at the half's own slope, read only where the slot's opportunity is inside the
grid's domain. Gate 1 (a Jeffreys-shaped landscape reproduces the reference-only ψ) passes to 1e-15 where shipped fails.
MEASURED (2026-10-02/03, pinned and fractional, the ladder): the zero controls improve — false gDNA at g00 366 / 266 /
358 / 329 → 13 / 9 / 12 / 8 fragments in calibration, the EM's g00 gDNA pool −87 to −96 % off capture — and every
capture-ON stratum loses, end to end as in calibration: stranded × ON transcripts 6.14 → 8.98 % (genes 1.12 → 1.68),
the deferred stratum 12.2 → 75.0 %; off capture flat. The extra tilt was holding enriched exons up against a landscape
that one population fit serves for every region class (under capture about 75 % of exon slots are enriched against
about 0 % of the others), so removing it alone drops them to the depleted mode. The partner tried — the shipped
landscape fitted per class (exon / non-exon) for ψ — is `ISSUES: an-unconditional-class-conditional-landscape`
(refused); paired with Part 1 it cancels both failures (g00 stranded × OFF false gDNA +15,043 → +16; stranded × ON
8.98 → 6.34 %) and ends flat-to-slightly-worse in scope against shipped (6.14 → 6.34 %; test chromosome stranded × OFF
7.39 → 7.70 %), better only on the deferred stratum (12.2 → 11.5 %). OPEN: a continuous per-object partner
that preserves locally supported enrichment while the density prior enters once. The owner's 2026-10-05
spectrum ruling supersedes the earlier suggestion of a library-level conditional split or capture-dependent
witness selection. No partner is chosen; an unconditional class split remains refused. The corrected count-map
experiment and failed point-count ruler screens are recorded at
`ISSUES: calibration-detects-capture-on-a-capture-off-library`.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-10-03_fl_arms/` (`harness/sitecustomize.py`, `RIGEL_BUILD=part1`,
`RIGEL_ARM=classwise`; `VERDICT.md`); `2026-10-02_part1/gates.py`, `classdiag.py`.

COUNT FOUNDATION INTEGRATED (2026-10-09). The factory-free count implementation and
observed conditional-binomial message claims now reside in main's working tree. The
owner accepts conservative DNA allocation in ambiguous cases and the eight reviewed
golden changes; bounds admission is deferred. This removes the separate intron prior,
not the fitted-landscape/reference duplication described above. The count solver's
Gaussian approximation, population admission and existing capture reader remain.
The implementation contract is in DESIGN; **Observed strand evidence** in EQUATIONS
owns the observation law. No new prior, detector, cutoff or tuned constant is added.

The final cleanup deletes the unused background fitter, old Gaussian message helper,
redundant source-class condition, duplicate composition-row input and native protocol
argument, and fixed-zero assembly dispersion state. Low-level dispersion inputs retain
actual nonzero solver-test consumers. Across 15 source files this is 138 additions and
643 deletions relative to the owner's tree before integration. Seven fresh cached pairs
preserve every calibration field and both belief arrays. Full test-chromosome and
LBX0190 pipelines preserve all payload/calibration arrays and the quantification digest
against the original validated candidate, through every cleanup and the actual main build.
The real identity reference has 254,319 transcript rows; no whole-genome debug dump is used.

The full-sweep belief-invariance gate includes endpoint beliefs and delivered level/composition
rows and masks; a compiled Gaussian restoration fails it. Five compiled likelihood defects
cover the 34 numerical own-claim cases; three actual binding restorations trip all three
interface cases. Both factory-removal gates fail on the original implementation. The
integrated suite re-derives **2,752 passed, zero failed/skipped/xfail**, with the eight
accepted goldens regenerated from main and exactly matching the reviewed candidate's
tables/scalars. No tolerance changes. Ruff and full preflight pass. Fresh g50 test-chromosome
runs in all four strata reproduce every candidate metric and pool exactly. The inherited
in-scope and deferred costs below are preserved, not erased by the golden updates.

Instrument: `.cache/rigel_runs/2026-10-09_foundation_landing/` (source manifest, owner snapshot,
integration patch, mutation/JUnit receipts, seven cached pairs, full pipeline identities,
golden review, normal-import stratum runs). The incomplete-environment first main-suite
invocation is retained separately; the complete conda-environment run has no skips. Main
source and goldens are integrated but uncommitted. The detector-free reader and its
uncertainty/opportunity validation remain release work; no broader reader scope is selected
by this engineering landing. All available real libraries remain strongly stranded.

LOCAL CLAIM AUDIT (2026-10-07). A native-builder replay of a real intron/boundary/exon triple
on test g50 ss0.99 ON separates the prior-derived and observational claims. With the strand
channel disabled, the factory contributes a 7.96-nat-range composition row to the boundary,
and 7.95 nats after transport into the exon. Removing only the factory removes that claim;
one certified RNA flux row remains. The observation counts, geometry, strand policy and
library coordinate origins are identical across this ablation. The claim is the intergenic
background's predictive prior, not a new independent observation of the intron's DNA.

The own exon/boundary strand rows also depend on the incoming belief through the frozen
variance. Consequently, subtracting the landscape from a final posterior cannot yield a
prior-free observation likelihood. The current factory rows must be kept distinct from raw
evidence when developing profile-based landscape training and a single-prior assembly.
No replacement assembly has been landed; deleting the factory blindly would remove an
existing source of unstranded composition information.

Changing the RNA coordinate origin affects sampled level-row representations. Boundary and
intron relative-height differences decrease with grid refinement in this triple; the sharp
exon row is not monotonically resolved. Large log differences in negligible tails are not
evidence of a large count error. The `split_live` bit changes the lane's other-column inputs,
but does not change the delivered composition in this particular triple. These limited checks
do not establish global locality or an end-to-end regression from either dependency.
Instrument: `.cache/rigel_runs/2026-10-07_density_checkpoint/`
(`audit_messages.py`, `message_audit.json`, `message_local_inputs.npz`).

FACTORY-REMOVAL SCREEN (2026-10-07; owner-authorized prototype). The Python factory input is
removed from every sweep, retaining the measured background. Two separate native worktree
contrasts then enable a single-strand intron's own observed strand row when its factory is
absent, and replace the message builder's belief-frozen Gaussian with the conditional
binomial observation likelihood. The final psi Gaussian, point-count landscape, additional
DNA reference tilt and capture reader remain unchanged: this is not the finished single-prior
model. The wiring-only build's factory-enabled control reproduces all 30 current count arrays
bit-for-bit. Full and pass-zero screens are recorded separately.

All 30 test-chromosome conditions were scored against certified slot truth, by object class
and admitted RNA strand. The table shows g50 conditions individually. Cells are **region /
boundary absolute DNA-count error (%)**, each divided by that axis's observed incidences;
these are not transcript/gene errors, and the two axes are not a conserved pooled count.

| test condition | current count-repair tree | factory removed | add intron strand messages | exact observed strand messages |
|---|---:|---:|---:|---:|
| ss0.50 OFF | 1.297 / 2.416 | 1.379 / 3.860 | 1.379 / 3.860 | 1.379 / 3.860 |
| ss0.50 ON | 4.950 / 9.956 | 3.688 / 2.543 | 3.688 / 2.543 | 3.688 / 2.543 |
| ss0.70 OFF | 1.221 / 2.045 | 1.304 / 3.118 | 1.326 / 2.528 | 1.323 / 2.521 |
| ss0.70 ON | 2.233 / 4.654 | 2.116 / 1.813 | 2.119 / 1.830 | 2.109 / 1.825 |
| ss0.99 OFF | 1.042 / 1.873 | 1.111 / 2.788 | 1.114 / 1.934 | 1.097 / 1.904 |
| ss0.99 ON | 1.370 / 1.660 | 1.352 / 1.153 | 1.342 / 1.151 | 1.337 / 1.134 |

Removal improves boundary counts on every contaminated capture-ON test condition. OFF
losses are primarily additional intron/boundary under-calling; exon estimates often improve.
The separate intron-strand wiring contrast recovers much of the stranded boundary loss.
The six zero-DNA controls are reported separately; these changes do not repair the old
reader's treatment of an uncertain positive count estimate as observed DNA.

Holding the original fitted landscape fixed leaves most of the OFF loss: g50 ss0.50
boundary error is 3.782% with the frozen landscape versus 3.860% when refitted, against
2.416% with the factory. At ss0.99 the corresponding values are 2.744 / 2.788 / 1.873%.
Factory-enabled fixed-landscape controls reproduce their original counts bit-for-bit.
Thus the lost local constraint, not merely changed landscape training, needs accounting
for in the single-prior replacement. No claim is made that profile training alone will
recover it, or that the intron background assumption is generally false.

The native checks pass 37 cases; an independent intron–boundary density reference passes
18 against a closed-form polynomial integral. Rebuilt binomial strand claims no longer
depend on the incoming belief. The local source factor stays inside RNA nuisance integration
and becomes composition-constant without strand information. This certifies one licensed
face, not the remaining splice/level messages or graph-wide observation independence.
Eight deliberate reversions/mutations fail the intended gates. No main-tree source, existing
tests, goldens or installed native binary changed, and no release accuracy claim is made.
Instrument: `.cache/rigel_runs/2026-10-07_factory_removal/` (`factory_ablation.py`, separate
native patches and wheel sites, `screen.py`, `test*_screen.json`, `calibration_tables.md`,
`prior_feedback.json`, `local_face.py`, mutation logs and `verification.json`). Native
worktree: `/private/tmp/rigel-density-20261007`.

COUNT-CANDIDATE READINESS (2026-10-08). A fresh run of the unchanged production suite
re-derives **2,733 passed**, no failures or skips. The frozen three-change count candidate
passes 2,721 and fails 12: four transfer-face references still encode the belief-frozen
Gaussian, and eight golden scenarios change. Scratch copies of the face references use
SciPy's conditional-binomial law without changing their maps, widths or assertions: all
15 face gates pass with the candidate, and the current implementation fails the four
changed gates. Full preflight passes for both builds, including all five instrument
self-tests; source/test/script lint passes. These results are a compatibility audit, not
a passing production landing: no test or golden was changed.

PORTABLE COUNT GATES (2026-10-08 follow-up). The three frozen count changes now have
38 portable specifications prepared under `tests/calibration/` in the isolated worktree:
35 source/wiring checks in `test_observed_strand_claims.py` and three whole-calibration
checks in `test_intron_prior_removal.py`. The former builds an exon/boundary/intron chain
directly; the latter uses the repository's small prior toy. Neither test file reads a
saved scan or the external prototype directory. The original implementation fails
**28/38**; the frozen count candidate passes **38/38**. Two inherited invariance checks
now explicitly require claims to be present, so an empty result cannot pass them.

Five compiled mutations of the actual frozen claim-builder header exercise all 35
native gates: restore the belief-frozen Gaussian, omit intron observations, lose all
claims, invent unlicensed claims, and ignore the RNA strand. The unchanged compiled
copy matches the frozen native's own-claim tables bit-for-bit. Three Python construction
mutations exercise all three removal gates: restore the factory, restore it only after
the initial solve, and omit the measured background. Failure coverage is **38/38** with
no collection errors or skipped cases. The factory candidate is still selected through
a scoped constructor replacement; this verifies portable specifications, not a direct
production integration. Instrument: `.cache/rigel_runs/2026-10-08_portable_count_gates/`,
JUnit receipts, compiled-header mutations, `run_factory.py`, and `verification.json`.
Main source/tests/goldens and the installed native remain unchanged. No new panel,
count-accuracy or performance result is claimed; the production suite's earlier twelve
candidate failures and the weak-strand concern below remain open.

DIRECT COUNT ASSEMBLY (2026-10-08). The isolated Python implementation now deletes
`FactoryRows`, `_IntronFactory`, `_Solve.factory` and their assembly use instead of
replacing a constructor at runtime. Auditing consumers showed that the old prototype's
retained background fit had no reader after its prior rows were removed. An earlier
portable gate unnecessarily required that unused computation. It is corrected: count
inference must not evaluate it, while an explicitly requested background measurement
must still depend only on intergenic observations. The measurement function and its
arithmetic remain available and unchanged.

Before this deletion the revised specifications fail **1/3** on the scoped-removal
prototype; direct assembly passes all **38** count gates. Restoring the count prior,
restoring only the unused fit, and changing the explicit estimator's pool each fail the
relevant gate, covering all three revised Python tests. The earlier 35 native gates
and compiled mutation receipts are unchanged. Seven fresh cached-condition pairs
compare **all 29 CalibrationResult fields plus two belief arrays**: every value is
bit-identical to the validated frozen count candidate. The native binary is the same
`b096f74e…`; this is a simplification, not a new accuracy or speed claim. Lower-level
factory-input plumbing and production reference/golden migration remain unfinished.
Instrument: `.cache/rigel_runs/2026-10-08_direct_counts/`, package and replay manifests,
JUnit reports and `mutations.json`. Full preflight passes all five instrument self-tests
with the direct package verified in the parent and all five children. Changed source and
portable tests pass Ruff; the documentation gate passes 11/11. This does not replace
the pending full-suite/reference migration on the eventual integration. Main
source/tests/goldens and the installed native are preserved.

COUNT-PLUMBING CLEANUP (2026-10-08). A separate count-only worktree now removes the
remaining intron-specific Python/native argument, factory row construction and storage,
factor-precision and negative-binomial helpers, and the unused `tau_fac` diagnostic.
The generic composition-factor kernel interface and explicit intergenic background
estimator remain. Observed intron claims enter directly; the exon-edge approximation,
landscape fit, count readout and existing capture reader are unchanged. No newer
density-message experiment is bundled into this package.

Seventeen tests solely exercising the retired factory are removed. Four old Gaussian
face references are migrated to the independent binomial reference; the intron-face
fixture isolates intron-originating claims so unrelated boundary claims do not change
that test's subject. No assertion tolerance is widened. The full isolated suite
re-derives **2,754 collected: 2,746 passed, eight failed, zero skipped**. The eight
failures are exactly the previously reviewed golden case set; no golden is updated.
The production tree retains its standing 2,733-test baseline and all owner changes.

Seven fresh paired cached calibrations again reproduce all 29 result fields and two
belief arrays bit-for-bit. The 35 portable native claim gates pass on the cleaned
binary; five actual compiled defects cover every gate, with an unchanged compiled
control reproducing the binary's claim tables exactly. Ruff passes across source,
tests and scripts. Full preflight passes all five instrument self-tests, with package
identity checked in the parent and all five children. This completes the mechanical
count cleanup and reference migration in the worktree; the golden disposition,
weak-strand stress concern and population/capture work still prevent a release-ready
claim. There is no new accuracy or runtime improvement to report from this no-op.
Instrument: `.cache/rigel_runs/2026-10-08_count_cleanup/` (`cleanup.patch`, complete
package manifest, retired-test list, full JUnit report, paired replay arrays,
`claim_mutations.json`, preflight import receipts and `verification.json`).

OBSERVATION-ONLY INTERFACE CLEANUP (2026-10-08). The isolated count candidate's
message builders no longer accept the unused incoming DNA belief or the two unused
strand dispersions. The Python policy carries only its fitted RNA sense fraction;
the tuple adapter, native policy-dispersion arguments and obsolete fixture fields
are removed. The count solver retains its separate beliefs and dispersion inputs.
There is no likelihood, opportunity, propagation, admission or readout change.

The three interface specifications first fail on the old binding. The new native
passes all 37 observed-claim cases; five compiled likelihood defects exercise the
34 numerical/wiring cases, and three actual native-binding restorations each trip
their corresponding interface case. The full suite re-derives **2,756 collected:
2,748 passed, the same eight golden failures, zero skipped**. Two additional cases
come from replacing the single belief-invariance check with three rejected-input
checks. No numerical assertion, tolerance or golden expectation is relaxed.

Seven fresh cached pairs preserve all 29 result fields and both belief arrays;
the maintained identity instrument also preserves the complete test-chromosome
transcript digest and every payload/calibration array, with fractional assignment.
Ruff and full preflight pass. This is a maintenance improvement, not a new accuracy
or speed result. The weak-strand/golden holds below remain unresolved, and the new
capture reader remains separate. Main source/tests/goldens and installed native are
preserved. Instrument: `.cache/rigel_runs/2026-10-08_release_foundation/` (before/after
sources, patch, package manifest, failed-first and mutation receipts, full suite,
paired replays, identity receipt and verification).

WEAK-STRAND FLOOR DIAGNOSIS (2026-10-09). A fresh array-only check of the frozen
`combo_extreme_removed_fit0.npz` confirms that all 18 nonflat delivered composition
rows are maximal at the all-DNA endpoint. Their DNA fractions span 0.4845–0.8354
and their posterior log-fraction variances 0.0670–0.4234, below the existing
one-nat-squared location limit. Slot 4 reads 0.79694 DNA share, against 5/38 true
DNA fragments. The structurally pure background is 0.05743 per opportunity, while
the fitted population mode is 0.30140. These are inputs to the first fit, not a
new fit or pipeline run. A symmetric Beta(1/2,1/2) reference truncated below 0.17
has median 0.70616: the illustration explains why a lower bound can acquire a
high posterior location without measuring that location. The soft delivered rows
are not claimed to be identical to this hard-truncation example.

Together with the recorded frozen-landscape and edge-source contrasts, this
supports the mechanism: under a dead strand channel, one-sided local evidence
leaves location to the reference, and posterior-based training can then propagate
that location through the population. It is a real model limitation under the
retained rules, not a map-arithmetic defect. The exact same causal attribution
has not been established on the ladder or an unstranded real library. The four
validated real libraries have fitted RNA sense fractions 0.00013–0.00300; all are
strongly antisense-stranded. Their successful runs provide no real-data coverage
of the unstranded failure mode. The cleanups also still need a real-library
identity comparison across the complete cleanup chain.

The integration recommendation is to accept this measured tradeoff without a
new solver change, while retaining both combined-stress fixtures and the separate
deferred-stratum reports. The owner accepted this tradeoff and authorized the reviewed
golden updates on 2026-10-09; see DESIGN's **Conservative allocation under uncertainty**.
The bound-only admission experiment is deferred until foundation integration is complete.
The refusal in `ISSUES: the-landscape-training-population-arms` remains in force;
an endpoint maximum alone would not establish a mathematical bound independently
of numerical support. Instrument: frozen stress inputs and
`.cache/rigel_runs/2026-10-08_release_foundation/review_assimilation.json` (input hash,
array-only check and illustrative reference calculation); real protocol values
are in the existing `2026-10-07_release_followup/real/*_candidate.json` receipts.

Every changed golden column was saved and read before any update. The largest absolute
change in an individual transcript is 1.447 fragments (`combo_moderate`); its total
transcript/gene absolute redistribution is 2.319/1.505% of the current counts. These are
output changes, not errors. `combo_extreme` has a much larger hidden pool movement:
gDNA including intergenic rises from 677.506 to 765.228, against read-name truth 700;
unspliced RNA falls from 265.106 to 176.688, against truth 246. Its transcript/gene
truth error changes only from 12.115/7.988% to 12.444/8.343%. The 1,000-fragment,
capture-OFF test has 36 spliced observations, fitted strand fraction 0.44737 and zero
strand discriminability. It is a weak-strand count check within the release's OFF scope,
not the deferred captured stratum and not a claim about whole-library RNA accuracy.

A separate current/observational-native × factory-kept/removed contrast establishes
that **factory removal causes this stress-case change**: the native builds produce
bit-identical count arrays at either factory setting because the strand channel is
inactive. The origin partitions reconstruct the production tally with no ambiguous
origin assignments. Keeping each original fitted landscape while removing the factory
separates the population feedback from the lost local constraint:

| `combo_extreme`, capture-OFF | Current | Candidate | Factory removed, original landscapes fixed |
|---|---:|---:|---:|
| gDNA fragments, including intergenic; truth 700 | 677.506 | 765.228 | 714.666 |
| Intron gDNA-count error (%) | 10.209 | 61.382 | 27.749 |
| Exon gDNA-count error (%) | 11.075 | 45.549 | 16.967 |
| Boundary gDNA-count error (%) | 12.752 | 16.283 | 5.588 |

Class errors divide by each class's observed unspliced incidences. Introns contain
8 true gDNA incidences out of 55; the current/candidate calls are 13.615/41.760.
Region and boundary incidences are not pooled as fragments. The fixed-landscape
intervention is diagnostic and uses a population the factory helped train; it is
not a proposed release method. Both local information loss and population feedback
remain relevant. This extends the known OFF limitation to a weak-evidence stress
case; it does not establish that the intron constraint is generally valid or license
a detector, a class correction, or a new strand threshold. Preserve this case in the
population/prior partner's in-scope validation. Do not erase the concern by updating
its golden or by reading only transcript changes. The previously accepted small
test-chromosome OFF tradeoff and deferred-ON scope remain as ruled.

Instrument: `.cache/rigel_runs/2026-10-08_count_readiness/` (`current.json`,
`candidate.json` and complete JUnit failure lists; isolated `face_tests/`; both
preflight receipts; `golden_outputs/`, `golden_review.json`; `stress_counts.py`,
the per-native JSONs and per-object NPZs). Production source/tests/goldens and the
installed native binary are unchanged. Reader integration is still outstanding.

GOLDEN DISPOSITION AUDIT (2026-10-08 follow-up). A read-only check of all 21 saved
current/candidate scenarios verifies identical table schemas, dtypes, object identities,
text fields and scan accounting. Independent arithmetic reproduces every saved
transcript/gene truth score and DNA/RNA pool total. The eight failing cases have these
absolute output movements, which are not accuracy improvements or acceptance bars:

| Fixture | Maximum transcript movement (fragments) | DNA pool movement (fragments) |
|---|---:|---:|
| `antisense_contained_ss90` | 0.728056 | +0.989977 |
| `antisense_overlap_ss90` | 0.000172 | −0.000305 |
| `combo_extreme` | 0.321837 | +87.722080 |
| `combo_moderate` | 1.447309 | +4.940199 |
| `gdna_heavy` | 0.312752 | +0.407343 |
| `gdna_light` | 0.016933 | +0.027814 |
| `nrna_heavy_ss90` | 0.003195 | +0.298159 |
| `nrna_moderate_ss90` | 0.011188 | +0.752400 |

The comparison uses saved fresh current outputs, rather than rounded or tolerance-matched
expected tables. Six cases move the DNA pool by less than one fragment; this does not
erase their pre-existing errors. In particular, contained antisense has no true DNA but
already reads 102.455 DNA fragments in the current output, versus 103.445 in the candidate.
The two combined stress cases remain separate count/reader concerns. No golden is updated,
no tolerance is widened and no case is declared accepted by this audit. The 0.7.1 release
was not rerun on these fixtures; the paired release-panel comparisons remain distinct.
Instrument: `.cache/rigel_runs/2026-10-08_count_readiness/golden_disposition.py` and
`golden_disposition.json` (input hashes, reconstructed scores and explicit hold dispositions).
This adds no calibration, EM or simulation run and confirms unchanged main source/tests.

WEAK-STRAND INPUT TRACE (2026-10-08). A read-only extension of the stress instrument
captures all three landscape refits for factory-kept and factory-removed controls in
`combo_extreme` and `combo_moderate`. Scan threads are one, EM is fractional with seed
zero, and calibration/EM thread counts retain defaults. The six final control count
arrays, including fixed-original-landscape controls, reproduce the earlier isolated
source-blur receipts bit-for-bit. Typed extraction also reconstructs the existing
delivered composition rows bit-for-bit. No training or message rule changes.

Admission does not explain either fixture's loss. All 13 regions of the extreme case
(3 intergenic, 4 intronic, 6 exonic) and all 11 of the moderate case already train at
every refit, with or without the factory. Removing only the variance cutoff therefore
adds no object in any of these 12 snapshots. In the extreme case the large error is
already present before the first population fit:

| Factory-removed extreme case, input to refit | Intronic DNA estimate; truth 8 | Exonic DNA estimate; truth 7 | Fitted population's modal density |
|---|---:|---:|---:|
| First, from the prior-free solve | 43.516 | 43.003 | 0.30140 |
| Second | 42.939 | 40.406 | 0.30140 |
| Third | 42.297 | 39.111 | 0.29171 |

The measured intergenic background is 0.057431. Keeping the factory instead gives a
modal density of 0.058894 at all three fits. This is consistent with the previously
measured population feedback; it is not evidence that later refits created the initial
local error. The moderate case retains a modal density near its intergenic background
under either factory setting and remains a separate control.

The existing native density evaluator was then read on frozen inputs, without fitting
a new population or changing counts. At each selected mixed region it compares the
point estimate's density with the measured background. For extreme-case chain slot 4
(region 2, an intron) before the third refit: observed incidences 38, true DNA 5,
inferred DNA 30.0549, `Var(log f_g)=0.15274`, retained reliability weight 0.49972.
The corresponding density is 0.32141. Its log-support advantage over background is:

| Input read at those same two densities | Log-support advantage (nats) |
|---|---:|
| Inferred DNA count treated as a Poisson observation | 27.0739 |
| Own observed strand columns, RNA amount profiled out | 0.4231 |
| Own observed strand columns, RNA amount integrated out | 1.2271 |
| Same observation integral with the existing delivered factors | 3.2375 |

The first row is the raw point-count likelihood before population smoothing or any
previous-prior weighting; it is not the full old landscape kernel. The other rows use
the already-declared RNA reference and measured strand probability. This comparison
demonstrates lost allocation uncertainty, but does not certify the inherited factors
as an independent joint observation likelihood. Their preference for excess DNA
persists. None of this establishes that profile training alone repairs the count loss,
and the background is a measured comparison point rather than exact object-rate truth.

The estimator/input comparison below retains these fixed objects, weights and grids.
Admission and the subsequently authorized reliability-weight contrast are separate;
do not restore an intron assumption or introduce a strand threshold to fit this case.
Instrument: `.cache/rigel_runs/2026-10-08_stress_evidence/` (`freeze.py`, `freeze.json`,
per-refit arrays, paired count arrays, `analyze.py`, `analysis.json`). The evidence
evaluation took 11.75 seconds, with no new population fit or count solve. These are
diagnostics, not new transcript/gene accuracy measurements.

FACTOR ATTRIBUTION (2026-10-08 follow-up). Across 64 nonempty single-RNA-strand
object/refit records from the two stress fixtures, the separate DNA/RNA level factors
are absent: their delivered support is in the composition channel. This does not mean
that it contains no DNA measurements; the EDGE builder projects boundary DNA evidence
into that channel. The own-read profile was checked against an independent maximization
of the two-Poisson likelihood. All 320 reference checks pass, including the exactly
unstranded plateau. Profiling versus integration changes the treatment of unknown RNA;
the difference is not an identified code defect or authorization to change its reference.

Two small control calibrations reproduce all six archived pre-refit count, opportunity
and message tables, and both final count arrays, exactly. They stop before EM. With face
roles, parameters, own observations and noncomposition factors fixed, one-at-a-time
EDGE-source withdrawals trace the extreme fixture's composition pressure to the two
ends of the same gene. For its slot 4 before the third refit, log support is 3.2375 with
both sources, 2.5544 without the left source, 2.0763 without the right source and 1.2271
without either. These are separate contrasts, not additive contributions. Removing
sources in the other gene leaves this object's curve identical. In the moderate fixture,
removing every EDGE source leaves most intronic support intact; slot 16 changes only
from 1.1969 to 1.1723. Its observed strand sources remain available.

Thus useful local message propagation and overstated population-training certainty are
distinct issues. The trace does not justify own-observation-only training or blanket
edge withdrawal; the latter's measured information loss is recorded under
`ISSUES: the-edge-density-floor-under-capture`. Preserve neighbour evidence and its
uncertainty. No population replacement, new capture prior or release accuracy result
follows from this attribution. The background remains a measured comparison point,
not exact per-object density truth. Instrument:
`.cache/rigel_runs/2026-10-08_evidence_attribution/` (`attribute.py`, `attribution.json`,
`trace.py`, `trace.json`, six frozen message contexts and `verification.json`).

FROZEN STRESS-INPUT COMPARISON (2026-10-08). The next authorized estimator/input
sequence is complete on both stress cases. It changes only the final landscape; the
first two landscapes, training objects, weights, density grid and uniform pseudo-region
strength are held fixed. Re-rendering the original individual curves recovers the mean
landscape within 1.4e-17 maximum probability difference. The pipeline control and a
direct calibration replay both reproduce the frozen count candidate bit-for-bit.
Every incoming refit belief and opportunity array is also identical. A fresh origin
partition validates the observed and true counts against the saved arrays.

The three successive changes are: replace averaging by the already-tested likelihood
estimator on the same old curves; replace only structurally known-DNA/empty input
curves; then replace only the remaining mixed-object curves with the typed observation
readout. The same reliability weights and smoothing strength apply throughout. All
six numerical fits have a global objective-gap certificate below 1e-6; the largest
measured gap is 2.91e-7. No equal-weight experiment or message-range change is included.

Cells are **DNA incidence-count error (%)**, divided by each row's observed incidences.
The first column is the factory-free count candidate, not the current production tree.

| Stress case · class | Frozen count candidate | Estimator only | + known-DNA inputs | + mixed evidence |
|---|---:|---:|---:|---:|
| Extreme · regions | 10.9194 | 9.6365 | 11.4826 | 9.1890 |
| Extreme · introns | 61.3818 | 55.5842 | 56.8745 | 37.6064 |
| Extreme · exons | 45.5492 | 39.0704 | 54.0145 | 49.5285 |
| Extreme · boundaries | 16.2825 | 7.8281 | 10.5412 | 7.9161 |
| Moderate · regions | 4.6833 | 5.1633 | 5.1633 | 3.9723 |
| Moderate · introns | 29.2016 | 38.2647 | 38.2647 | 20.4325 |
| Moderate · exons | 14.1795 | 13.6094 | 13.6094 | 13.4717 |
| Moderate · boundaries | 4.6917 | 4.2946 | 4.2946 | 5.9787 |

Intergenic counts remain exact in every arm. Mixed evidence improves introns in both
fixtures relative to both controls, but the extreme exons and moderate boundaries
worsen against the frozen candidate. This is not a sufficient population partner:
do not advance it to release benchmarks or claim that preserving uncertainty alone
repairs local count inference. It is a frozen final-refit diagnosis, not a self-consistent
calibration or a prior-once test. No candidate EM or new transcript/gene benchmark ran.
The two fixtures are stress checks, not sources of class weights, cutoffs or biological
assumptions; their nascent share is not a design driver.

Instrument: `.cache/rigel_runs/2026-10-08_stress_evidence/` (`fit_stress.py`,
`*_input_curves.npz`, `*_fits.json`, fitted arrays, `replay_stress.py`, `*_replay.json`
and per-arm count arrays). Each numerical fit took under one second; all evidence curves
for a case took about two seconds. These are execution receipts, not a tool speed-up.
The unresolved inherited-factor, population/prior and representation questions remain.

FULL-PANEL CONTINUATION (2026-10-07). After the owner accepted the modest test-chromosome
OFF tradeoff, the frozen candidate (factory absent, intron strand sources restored,
observational strand messages) was scored on all 16 ladder and both four-condition gap
panels against a fresh current-tree baseline. The same truth-array order and observed
totals are checked for every condition; object classes and admitted strands remain separate.
Stranded count errors remain close to the current tree in both gap directions. Unstranded
losses grow on the larger substrate, and have two distinguishable mechanisms.

Cells below are region / boundary **DNA-count errors (%)**, not transcript/gene errors.
The fixed-landscape diagnostic keeps the factory-enabled candidate binary's final landscape
unchanged. Its factory-enabled controls reproduce the unfrozen counts bit-for-bit; on both
these unstranded conditions those counts also equal the current-tree results.

| condition | current | candidate | candidate, original landscape held fixed |
|---|---:|---:|---:|
| RNA-long g50 ss0.50 ON | 10.009 / 9.582 | 41.676 / 32.913 | 11.665 / 12.353 |
| ladder g98 ss0.50 OFF | 0.788 / 5.711 | 2.465 / 9.124 | 2.173 / 8.714 |

On RNA-long ON, changed population fitting accounts for most of the additional error.
At the first refit, admitted exonic training regions fall from 3,292 to 289, out of 8,384
structurally eligible regions; the number with a composition channel falls from 5,339 to
2,742. Removing the factory changes both composition availability and posterior-precision
admission. Holding the old population fixed substantially repairs the result, but neither
this trace nor that intervention isolates admission, weighting and kernel construction from
one another. It does not yet prove that replacing kernels alone will recover the lost mode.
The large final error sits mainly in exons, where gDNA is under-called; it is not a simple
remaining intron-count offset.

On the ladder's g98 OFF condition, holding the population fixed repairs much less. Most
of that loss is direct: the removed intron constraint supplied local composition information.
A low absolute library-scale error can be large beside the residual RNA in a 98%-DNA sample.
These results require both a population-training check and an accounting of local unstranded
information. No class multiplier, uncertainty cutoff or capture gate follows from them.
The original small test-chromosome OFF tradeoff is not itself a reason to reject the candidate;
the full-panel mechanisms must be evaluated separately.

Paired end-to-end follow-up preserves the stranded gap results and the original RNA-short
repair, but the RNA-long unstranded ON condition fails the released 0.7.1 transcript/gene
bar. The ladder g98 OFF count loss mostly enters unspliced RNA, with a much smaller mature
transcript change. Neither the calibration metric nor the transcript table substitutes for
the other. The owner subsequently confirms that unstranded capture-ON is deferred, not a
0.8.0 requirement, and authorizes deeper investigation for robustness. That loss is not by
itself a release veto. The owner also authorizes a separate evidence-curve admission
prototype; composition requirements, structural exclusions and weights stay fixed in its
first contrast. Production integration still needs the remaining in-scope validation and
landing checks.

Instrument: `.cache/rigel_runs/2026-10-07_full_panel_counts/` (`screen.py`, matched
`*_screen.json`, per-condition arrays and `calibration_tables.md`; `diagnose.py` and
`diagnosis.json` retain each refit's training counts and the fixed-landscape controls).

VALIDATION COMPLETION (2026-10-07). The same frozen count candidate now has 54 paired
end-to-end conditions: test chromosome 30, ladder 16 and four in each gap direction.
Every in-scope stratum summary beats 0.7.1 for transcript and gene error. Individual
high-DNA test conditions remain worse: g98 ss0.70 ON reads 71.94 / 45.54% versus
55.98 / 39.60% on 0.7.1 and 69.50 / 42.44% on the current tree. Do not turn a
stratum summary into an every-condition claim. The deferred RNA-long unstranded ON
loss also remains visible at 43.22 / 11.26% versus current 6.55 / 0.94%.

Four serial whole-genome pairs also completed, without calibration debug dumps.
The candidate/current DNA fractions (%) are LBX0190 8.50/8.34, LBX0588 92.65/91.14,
MO3021 15.87/15.46 and VCaP 23.95/23.44. VCaP moves toward its read-name mixture
of 25.18%. Transcript absolute changes, normalized by the current transcript total,
are 0.13%, 14.18%, 0.14% and 4.03%; gene changes are 0.07%, 6.03%, 0.06% and 0.28%.
These are output changes, not accuracy scores. Half of LBX0588's transcript change
lies in ten genes. The old capture reader is active in both arms, so changes in its
count-derived effective lengths can redistribute isoforms. No real-data tuning follows.
Source-landing checks and a detector-free reader remain outstanding.
Instrument: `.cache/rigel_runs/2026-10-07_release_followup/` (`panel_comparison.json`,
`report_all.py`, paired manifests and hook receipts; `real/manifest.json`,
`real/comparison.json`, `real/change_diagnosis.json` and `report_real.py`).

ADMISSION AND ESTIMATOR CONTROLS (2026-10-07; owner-authorized). Seven cached conditions
separate deletion of the variance cut from changing the population estimator. The native
candidate and factory removal are fixed. The first-refit counts, variances, opportunities,
composition flags and grid are identical before the admission intervention. All controls
retain the structural exclusions, composition requirement, zero-count anchors, posterior
reliability weights and consumer-domain grid. Cells are **region / boundary DNA-count
absolute error (%)**, with separate observed-incidence denominators, not transcript/gene errors.

| condition | current tree | factory-free candidate | candidate, variance cutoff removed |
|---|---:|---:|---:|
| RNA-long g50 ss0.50 ON | 10.0085 / 9.5817 | 41.6758 / 32.9132 | 40.5559 / 33.6206 |
| ladder g98 ss0.50 OFF | 0.7878 / 5.7109 | 2.4646 / 9.1235 | 2.4827 / 9.1821 |
| RNA-short g50 ss0.99 OFF | 0.7854 / 2.9059 | 0.9246 / 2.9482 | 0.9258 / 2.9478 |
| RNA-long g50 ss0.99 ON | 1.4974 / 3.4148 | 1.4657 / 3.3480 | 1.4659 / 3.3479 |
| test g00 ss0.99 ON | 0.0029 / 0.0196 | 0.0035 / 0.0239 | 0.0037 / 0.0234 |
| test g50 ss0.99 ON | 1.3702 / 1.6599 | 1.3371 / 1.1336 | 1.3343 / 1.1336 |
| test g50 ss0.50 ON | 4.9502 / 9.9556 | 3.6884 / 2.5428 | 3.6977 / 2.5379 |

Removing the cutoff admits 2,742 rather than 289 exons on the first RNA-long ON refit,
yet does not repair the regression. At that refit, the newly expanded exonic population's
fitted DNA totals are 783,235 against 1,313,802 true incidences. Thus admission alone is not
the cure; preserving a poor point estimate is not equivalent to preserving its uncertainty.

A second contrast freezes each refit's selected objects, weights, grid and individually
rendered curves. Weighted recombination reproduces the original landscape to rounding;
replaying all three frozen mean landscapes reproduces final count arrays bit-for-bit.
Replacing only averaging with a mixture-likelihood fit of those same curves gives:

| condition | original frozen mean | likelihood fit of the same old curves |
|---|---:|---:|
| RNA-long g50 ss0.50 ON | 41.6758 / 32.9132 | 20.7739 / 19.3768 |
| test g00 ss0.99 ON | 0.0035 / 0.0239 | 0.1568 / 0.2559 |

The captured-row improvement is concentrated in exons (44.2483 → 21.6501%); intron error
worsens (10.8513 → 18.5429%). The zero control's false DNA rises from 32.66 to 1,455.90
region incidences and from 23.33 to 249.69 boundary incidences; these axes are not additive.
**The estimator-only arm is not a landing candidate.** Its inputs still contain inferred
counts, population smoothing and the old prior-weighted low-count curves; they are not
certified observation likelihoods. These measurements justify replacing the inputs, not
declaring that the proposed honest-profile model succeeds or fails. No new end-to-end
accuracy claim follows from these calibration-only controls.

The small-grid reference also exposed a numerical dilution defect: with a sufficiently
large weight on an exactly flat likelihood, a population-normalized stopping tolerance
could accept the starting distribution. Analytically cancelling constant objective terms
fixes this without an information cutoff. Ten likelihood gates pass, including a raw
two-Poisson mixture checked against an independent optimum and preservation of an
exclusively supported 1-in-10,001 minority. Six deliberate defects are caught. Admission
and frozen-rendering checks add twelve gates and ten caught defects; every gate is hit.
On the zero control, multiplicative updates did not reach the requested tolerance within
100,000 iterations. An independent constrained optimizer certifies the same objective to
a gap below 1e-6; its generic termination flag alone was insufficient. Numerical convergence
and support remain explicit implementation requirements, not biological tuning parameters.
Instrument: `.cache/rigel_runs/2026-10-07_profile_admission/` (`screen.py`, per-refit NPZ
arrays, `current.json`, `candidate.json`, `cutoff_removed.json`, `frozen_curves.py`,
`fit_frozen.py`, `replay_frozen.py`, `zero_optimizer.json`, mutation logs and
`verification.json`); `.cache/rigel_runs/2026-10-07_profile_fit/` holds the likelihood reference.

OBSERVATION-LIMIT INPUT CONTRAST (2026-10-07). On the same frozen training objects,
weights, grid and likelihood objective, replace only structurally pure-DNA and empty
objects' rendered curves with their own Poisson observation likelihoods. For an empty
object, integrating its own RNA observation with the existing reference leaves
`exp(-rho Eg)` times a density-independent constant. Other mixed-object curves are
unchanged, so this is a partial input diagnosis, not a finished evidence model.
Fresh replays of the old-curve controls reproduce their previous count arrays bit-for-bit.
Cells remain **region / boundary DNA-count error (%)**, not transcript/gene error.

| condition | likelihood fit on old curves | replace certified observation limits |
|---|---:|---:|
| test g00 ss0.99 ON | 0.1568 / 0.2559 | 0.0023 / 0.0201 |
| RNA-long g50 ss0.50 ON | 20.7739 / 19.3768 | 20.9468 / 19.2018 |

The zero control's false region DNA falls from 1,455.90 to 21.13 and false boundary DNA
from 249.69 to 19.58. Of the remaining region count, 21 incidences are structurally
assigned intergenic DNA despite RNA oracle truth; this contrast does not repair annotation
misclassification. The frozen factory-free baseline, before changing the estimator, was
32.66 / 23.33 false incidences. On RNA-long ON the regional total hides an intron loss:
18.5429 → 24.1715% intron error; exons remain 21.6501 → 21.6612%. Do not promote this
partial candidate on the zero-control improvement alone.

The zero-control decomposition retains all other inputs and the same objective:

| known/empty-row input | false region DNA | false boundary DNA |
|---|---:|---:|
| old population blur and previous-prior weighting | 1,455.90 | 249.69 |
| remove only previous-prior weighting | 2,817.97 | 183.81 |
| remove only population blur | 2,943.64 | 217.53 |
| raw observation likelihood: remove both | 21.13 | 19.58 |

This is an interaction; neither isolated deletion explains the recovery. A smoothed or
prior-weighted distribution is not an observation likelihood merely because it is a curve.
All population optima used here have finite objective and a global objective-gap certificate
below 1e-6. A falsification caught an invalid-zero-likelihood certificate in the diagnostic
optimizer; its guard is fixed. A tighter two-point request was explicitly refused rather
than silently accepted; the independent optimum check uses the experiment's stated 1e-6
numerical tolerance. This is not a panel-selected biological threshold.

The splice-input audit now has a raw conditional observation factor for a single licensed
face, keeping individual route counts and opportunities. Its derivation and the unstranded
identification limit are in `EQUATIONS.md`, **Conditional observations at one splice face**.
At this checkpoint an assembly policy for the local RNA route shares was not selected;
the subsequently authorized profiling reference is recorded below. Graph-wide duplicate
observations remain unverified. Forty scratch observation/optimizer checks pass and thirteen
deliberate defects are caught, with every gate exercised. No production source was changed.
Instrument: `.cache/rigel_runs/2026-10-07_observation_inputs/` (`limits.py`, `fit_inputs.py`,
`replay.py`, per-refit arrays, separate ablation receipts, raw splice reference and mutation logs).

ROUTE-PROFILE ROBUSTNESS STOP (2026-10-07). The owner authorizes exploring route profiling
and requires a general, simple model, not optimization to a simulated dataset. A one-junction
reference profiles the route share by a concave scalar solve, with no new prior or tunable.
Independent maximization in the original route coordinate, analytic limits, component units,
strand reversal and opportunity contrasts pass. The initial fixed-share placeholder fails
22 of 33 arithmetic gates. The implemented reference plus counterexample characterization
has 40 passing scratch checks, all exercised by nine deliberate defects. These certify the
arithmetic and the failure diagnosis, not suitability as general evidence.

The predeclared generality condition **fails before any panel run**. Hold the target DNA
rate at 20 and RNA at 80, give both components and the junction unit opportunity, strand
probability 0.75, and half the RNA each route. Unspliced expected columns are exactly
40/20. Change only junction capture; all tabulated counts are integral expected means,
with no sampling, fitted truth input or simulation. The factor profiles routing at each
candidate density with the RNA amount held fixed:

| junction capture multiplier | observed junction count | source factor's implied DNA rate | true DNA rate | log-likelihood advantage over truth |
|---|---:|---:|---:|---:|
| 0.1 | 4 | 36.3636 | 20 | 0.1254 |
| 1 | 40 | 20 | 20 | 0 |
| 10 | 400 | 3.6364 | 20 | 11.1040 |

The implied rate attains the saturated multinomial likelihood; the incorrect opportunity
at the true rate cannot. Tenfold depth multiplies each nonzero log-likelihood advantage
by ten while preserving the rate bias. Supplying the true capture-weighted junction
opportunity restores the true-rate optimum. This is an oracle control diagnosing the
opportunity assumption, not an implementable capture estimator or an end-to-end error.
The algebraic equivalence for multiple routes is derived in `EQUATIONS.md`,
**Capture-opportunity nonidentification at a splice face**.

Independent exhaustive fragment placement confirms that relative capture is not an
impossible perturbation. At a donor separating a first exon from an intron, a probe in
the second exon enriches spliced crossings while missing continuous crossings. A probe
in the intron does the reverse. Twenty-four small placement/length/binding contrasts include
both equal lengths and gaps in each direction; no single scalar correction survives them.
This prices a structural possibility, not its prevalence in real data. The inherited
shared-capture premise also exists in the point-rate maps; profiling exposes, not creates,
this issue. Raw evidence curves and the prior-once goal are not rejected by this result.

Stop before RNA integration, multiple-route implementation, frozen population replay or
native integration. Do not select a correction factor, admission cutoff or probe-layout
branch. The next proposed checkpoint is to establish what relative capture/opportunity
information is identified by the existing fragment observations before treating splice
ratios as transported density evidence. Until then this is a diagnostic reference, not a
landing candidate. Existing factory-removal validation and deferred-stratum scope stand.
Instrument: `.cache/rigel_runs/2026-10-07_route_profile/` (`PROTOCOL.md`, `route_profile.py`,
`capture_mismatch.py`, the intentionally failing generality receipt, tests and mutations).

LOCAL-FOOTPRINT CLARIFICATION (owner, 2026-10-07). The stop above applies to treating an
isolated face's exact shared-capture formula as general evidence. Extending it to a rejection
of useful local capture inference was too broad: the isolated comparison did not use the
other enriched objects reached by the same probe. The owner's intended model is approximate
local co-enrichment at region/boundary resolution (`DESIGN.md`, **Capture has a local footprint**),
not arbitrary independent capture averages or an identical multiplier across all neighbors.

A follow-up overlap calculation uses one junction-spanning probe with a 50-base part in
each of two 300-base exons, separated by a 500-base intron. The spliced molecule sees one
contiguous probe; continuous molecules can bind either part. Enumerate contained, ordinary
crossing and junction-crossing fragments separately with their own component lengths.
The established diagnostic law is `1 + binding * best contiguous overlap`. At binding=1,
the expected capture multipliers (relative to uniform sampling) are:

| RNA / DNA length | RNA exon-contained | RNA ordinary boundary | RNA junction | DNA exon-contained | DNA ordinary boundary |
|---|---:|---:|---:|---:|---:|
| 78 / 78 | 6.717 | 35.091 | 69.182 | 6.717 | 35.091 |
| 78 / 250 | 6.717 | 35.091 | 69.182 | 26.000 | 46.080 |
| 250 / 78 | 26.000 | 46.080 | 91.161 | 6.717 | 35.091 |

The opposite flank has the same expectations by symmetry. All five channels are enriched
in all nine length/binding contrasts (binding 0.1, 1 and 10); a far boundary at the other
end of the first exon remains unenriched. A closed-form sum independently checks the
contained expectations. This establishes local co-enrichment and unequal object averages
in the same physical example. It does not fit a coupling strength or guarantee that every
adjacent object is enriched. No calibrated density or end-to-end accuracy is measured.
Receipt: `.cache/rigel_runs/2026-10-07_local_footprint/junction_probe.json`.

EXISTING-MESSAGE EVIDENCE INTERFACE (2026-10-07). The factory-free exact-strand build's
delivered composition rows can be read inside the exact own-observation RNA integral,
without a second propagation system. The Python diagnostic implements this for a single
admitted RNA strand, with pure-DNA and zero-opportunity limits. It preserves the existing
transport maps, discrepancy treatment and directionality. It is not yet a production
consumer or a certification of all delivered factors as independent observations.

Fresh native replays on test g00 ss0.99 ON, g25 ss0.70 ON and g50 ss0.50 ON compare the
fitted DNA prior with no DNA prior on the same solve grid. Every delivered composition
row and every present RNA cube row is bit-identical (63, 101 and 27 delivered cube slots).
The largest changes in final DNA share are 0.803, 0.903 and 0.604, respectively, so the
prior intervention is effective. Re-reading the original prior with diagnostics preserves
all count-belief arrays bit-for-bit. This establishes prior independence for these rebuilt
rows on a fixed grid, not invariance to a changed library coordinate or shared technical fit.

At the zero control's largest false-count single-strand exon and boundary, the old
point-count Poisson input overstates the evidence against the measured background:

| object | fitted DNA incidences | point-count log support at background, relative to its maximum | observation/message log support at background, relative to the sampled maximum |
|---|---:|---:|---:|
| exon 762 | 5.760 | -49.988 | -5.748 |
| boundary 256 | 22.812 | -212.428 | -6.428 |

The background density here is 3.65012e-6/bp. The new figures are likelihood contrasts,
not capture weights or new count estimates; the sampled rate is not an oracle realized
count. Own observations alone penalize zero DNA by only about 0.24 nats at each object.
The delivered-row trace reconstructs the native rows bit-for-bit and identifies the
additional pressure: exon 762 receives a strong boundary strand claim (27 opposite-column
reads out of 1,447), and boundary 256 receives its left exon's composition and its right
exon's DNA lower bound. These claims derive from observed strand imbalances, not a leaked
DNA prior. This finding alone does not establish a faulty message or justify suppressing
all positive DNA evidence in an RNA-only realization. Shared-fragment dependence and model
mismatch remain possible contributors, and the evidence must retain finite uncertainty.

The point-count control fails 13 of the original 21 checks. The completed adapter has
24 passing scratch checks without integration warnings; eleven deliberate defects/reversions
plus the control exercise every gate. One transported intron/boundary factor agrees with
the independent two-object integral to the stated interpolation tolerance. General splice
and level rows remain approximate. No capture prior, new biological constant, propagation
rule, production source change or end-to-end reader improvement is claimed.
Instrument: `.cache/rigel_runs/2026-10-07_message_evidence/` (`message_evidence.py`,
`test_message_evidence.py`, falsification/mutation receipts, `audit_delivery.py`,
`delivery_audit.json`, `trace_zero.py`, `zero_trace.json`).

NATIVE EVIDENCE CHECKPOINT (2026-10-07). The same one-RNA-strand integral now has an
isolated C++ evaluator with analytic physical limits, numerical tail control and no
count-path call site. Forty-two reference/interface checks pass; twelve compiled defects
and reversions exercise every gate. Two numerical defects were found and corrected:
a cancelled Simpson error estimate, and loss of significance from subtracting far-tail
Poisson log likelihoods. Neither fix changes the statistical model or adds a rate cap.
Replaying the frozen training snapshots on test zero and RNA-long unstranded ON preserves
all selected objects, beliefs, grids, weights and final count arrays bit-for-bit.
The next contrast changes only the remaining mixed-object curves of the frozen final
population refit. Both-strand cube integration, composite-evidence dependence, coordinate
invariance and production throughput remain unresolved; no new reader accuracy is claimed.
Instrument: `.cache/rigel_runs/2026-10-07_native_evidence/` (isolated source and package,
`test.log`, `mutations.json`, delivery audit and numerical timing receipts) and
`.cache/rigel_runs/2026-10-07_mixed_inputs/` (frozen inputs and count-identity receipts).

MIXED-INPUT AND COORDINATE FOLLOW-UP (2026-10-08). Replacing the remaining mixed
curves of the frozen final population refit, while holding objects, weights, grid and
objective fixed, changes RNA-long unstranded-ON region/boundary DNA-count error from
20.9468/19.2018% to 17.1196/19.1979%. Exon error improves from 21.6612% to 17.3622%,
but intron error worsens from 24.1715% to 30.7955%. This is an attribution control,
not a self-consistent new calibration or a landing candidate. The zero-control final
selection has no mixed rows; its identical result cannot certify the future reader.
The constrained population optimizer needed multiplicative initialization to reach the
unchanged 1e-6 objective-gap requirement; the certified final bound is 6.21e-7.

A held-level counterexample demonstrates a remaining interface loss. `profile_of_level`
substitutes the target's observed total when converting absolute DNA density to composition.
Feeding the resulting row to an integral over expected RNA amount reads a different density
factor. With target counts 4/4, unit opportunities, strand fidelity 1/2 and a fixed
one-sided DNA level from source count 20, the native converted row prefers density 8.2;
reading the same held factor in its original density coordinate prefers 13.8. Both modes
are on a 0.1-spaced diagnostic grid. This does not establish identical capture between
neighbors or certify the original factor as an independent likelihood. It establishes
that the final fused composition row is not a lossless density-evidence interface.
Preserve the already-existing composition and absolute-level channels before the substitution;
no new propagation or spatial coupling follows. The algebra is in `EQUATIONS.md`,
**Absolute level factors retain their coordinates**.
Instrument: `.cache/rigel_runs/2026-10-07_mixed_inputs/` (frozen curves, optimization
certificate and matched final-refit replays); `.cache/rigel_runs/2026-10-07_native_evidence/`
(`check_level_coordinates.py`, `level_coordinate_check.json`).

SEPARATE-CHANNEL IMPLEMENTATION (2026-10-08). The isolated density evaluator now reads
the existing DNA level at candidate DNA density and the existing RNA level inside the
RNA-amount integral. Composition remains a ratio factor. No new message pass, spatial
model, prior or count-path call site was added. Existing diagnostics already expose the
received channels; rebuilding their final assembly is bit-identical. The control fails
11/19 coordinate gates and the empty extractor fails 5/6 wiring gates. The corrected
primitive plus earlier tests and wiring tests pass 67/67; eighteen native/adapter mutants
and five extraction mutants exercise all 67. Native sources and binary are restored after
the campaign. Main production source and installed native code are unchanged.

The frozen RNA-long unstranded-ON final population selection contains 1,159 mixed objects;
85 have a separate RNA level, none has a separate DNA level. The other 1,074 mixed
composition rows are exactly unchanged. Recomputed controls match the archived curves
exactly. Corrected curves plus controls take 43 seconds with two numerical workers;
the unchanged population objective meets its 1e-6 certificate with a 7.58e-7 bound.
The matched count replays take about 15 seconds each. DNA incidence-count errors (%):

| class | previous mixed curves | separate factor coordinates |
|---|---:|---:|
| Regions | 17.1196 | 16.7259 |
| Boundaries | 19.1979 | 18.8416 |
| Introns | 30.7955 | 30.9444 |
| Exons | 17.3622 | 16.9364 |

This repairs the interface but does not resolve the deferred regression or validate a
production population partner. The intron loss persists. The zero-control final selection
again contains no mixed objects. Directly reading its worst exon leaves background log
support unchanged at -5.7478 nats; its worst boundary changes only from -6.4282 to -6.4290
nats, relative to each curve's sampled maximum. Thus this coordinate correction does not
explain or remove their false-capture risk. No capture weights were inferred in this test.

Three test-condition audits retain bit-identical count outputs and prior-independent
received rows/cubes on fixed grids. That does not establish invariance when whole-library
coordinate origins rebuild the tables, nor independence from `split_live`, nor independence
of factors that share observations. In particular, earlier EDGE/LEVEL face builders can
already encode absolute evidence through observed totals; correcting final delivery alone
does not make those composition factors exact. The next tests must distinguish those
inherited approximations from the authorized admission contrast, keeping each intervention
separate. Both-strand evidence and production throughput remain uncompleted. The numerical
primitive is a supported prototype; its population consumer is not a landing candidate.
Instrument: `.cache/rigel_runs/2026-10-08_channel_evidence/` (source snapshots and isolated
package, `falsification.log`, `tests.log`, `mutations.json`, `delivery_mutations.json`,
`*_freeze*.json`, `long_unstr_on_evaluation2.json`, optimizer certificate, matched count
replays and `channel_audit.json`).

BOTH-STRAND DELIVERY (2026-10-08). The diagnostic adapter can now retain the existing
two RNA-level factors at an object admitting both RNA strands. It reads the native cube
table already emitted by final message assembly; it does not rebuild propagation or
apply the single-strand side exclusion. That table contains each strand's intersection
of held levels and its own lower-side bound. Composition and DNA factors reuse the
previous typed adapter. Slot identity, presence bits and each table's density coordinate
are preserved. No native function, count call site, prior or population policy changes.

The former single-strand-only interface fails all 16 adapter specifications; the completed
adapter passes all 16. Fifteen actual adapter mutations each fail, collectively exercising
every specification. The first mutation attempt exposed a weak fixture: held levels
entirely masked the own lower-side bounds. The strengthened fixture makes both matter;
dropping either own bound now fails. A separate row-permutation check prevents silently
reading another object's profile. These are delivery checks, not nuisance-integration checks.

Three frozen test-chromosome input contexts each contain 198 both-strand objects. Native
prepare/pass/assembly and the diagnostic reads take under one second together, excluding
process startup. Presence counts are:

| condition | DNA factor | RNA+ factor | RNA− factor | both RNA factors | neither RNA factor |
|---|---:|---:|---:|---:|---:|
| g50 ss0.99 ON | 83 | 63 | 75 | 18 | 78 |
| g50 ss0.50 ON | 0 | 24 | 4 | 1 | 171 |
| g00 ss0.99 ON | 48 | 24 | 52 | 13 | 135 |

Every present RNA profile, axis and origin equals its native cube entry. Prepared inputs,
received messages and cube storage are unchanged after extraction. A repeated final
assembly reproduces the count-facing rows and all present cube entries bit-for-bit.
These checks do not run calibration or the EM; they establish neither final-count identity
under integration nor an accuracy improvement. None of these contexts supplies a composition
factor at a both-strand object, so that combination is covered by the explicit fixture.
The both-strand RNA-amount/tilt integration with delivered factors remains unimplemented;
finite support, shared observations and inherited face approximations remain open.
Instrument: `.cache/rigel_runs/2026-10-08_both_delivery/` (`both_delivery.py`,
`before_complete.xml`, `candidate.xml`, `mutations_initial.json`, `mutations.json`,
`census.py`, `census.json`, `verification.json`).

BOTH-STRAND DENSITY REFERENCE (2026-10-08). The extracted RNA profiles now have a
standalone consumer using the existing total-RNA measure and three tilt hypotheses.
At fixed tilt it shifts both profiles into log total RNA amount, merges their knots
and reuses the certified native one-amount integral. It integrates over tilt outside
that call. A separate SciPy reference reverses the integration order and evaluates
the original Poisson columns. No native source or production count call site changes.

The initial own-observation control fails 8 of 18 factor checks. The generic reference
then passes 26 observation, factor, normalization, coordinate, limit, missing-profile,
narrow-profile, refinement and input-preservation checks. A subsequent rule audit
finds that finite profile values alone omit existing hypothesis support: the native
count solver uses a delivered RNA profile to exclude the opposite pure-strand atom.
That audit fails 4 of 5 additional checks before the diagnostic consumer passes the
presence bits explicitly. This is preservation of **The AMBIG tilt's hypothesis space
is {pure +, pure −, mixed} — the tilt atom**, not a new strand-presence inference rule.

The completed reference passes 31/31 checks in 26.22 seconds. Twelve actual defects
collectively exercise every check; two further numerical defects are also detected.
Dropping the profile-derived tilt breakpoints fails the narrower of the two deliberately
peaked-profile cases. Dropping the pure-strand hypotheses fails neutral-profile and
unwitnessed-support checks. All source mutations are restored. An initial mutation that
changed inputs at every quadrature call failed through loss of numerical convergence;
the final input-mutation check instead changes a profile once and fails the actual
input-preservation assertion. The weaker initial receipt is retained separately.

These checks establish a reference for the declared composite-factor calculation, not
an independent joint likelihood for the inherited messages. Both-strand own lower-side
bounds can reuse target observations. The existing witness rule remains explicit, and
changing it requires its own model decision. Finite input support is still unresolved.
Whole-genome cost is unmeasured: the scalar reference partitions tilt using pairs of
profile knots and must not be called a production throughput result. No both-strand
objects have been added to population training, no capture weights have been read and
no count or transcript accuracy benchmark was run for this reference.
Instrument: `.cache/rigel_runs/2026-10-08_both_evidence/` (`both_evidence.py`,
`reference.py`, `before.xml`, `witness_before.xml`, `candidate.xml`, `mutations.json`,
`numerical_mutations.json`, `generic_only/`, executed snapshots and `verification.json`).

BOTH-STRAND NUMERICAL COST (2026-10-08). A table census and bounded scalar probes reject
a direct port of the outer reference loop as the next production step. Representative
two-profile objects have 992–1,090 tilt breakpoints from the pairwise knot construction.
The stranded test object (slot 2780, 728/1,407 observed columns) exceeds a 20-second wall
budget for one density, after 4,523 conditional evaluations. This is a measured cost
barrier, not a new count-error or transcript-error result.

An isolated exact compactor deletes only interior knots within constant stretches,
retaining endpoints, all changing segments, absolute heights and profile presence.
The stranded object's neutral composition table shrinks from 172 knots to two and
its positive RNA profile from 172 to 37. Four of six representation checks fail before
the change; all six plus the prior 31 numerical checks pass afterward. Three actual
mutations exercise all six new checks. Interpolation at every original knot and interval
midpoint is identical on all 210 saved object/condition inputs; presence is preserved.
The mathematical reason is **Constant interpolation segments need only their ends** in
`EQUATIONS.md`. No nonconstant curvature or small nonzero change is discarded.

Serial scalar probes on the same saved inputs, without calibration or EM:

| test condition / slot | original wall seconds | compact wall seconds | original / compact conditional calls | log-evidence difference |
|---|---:|---:|---:|---:|
| g50 ss0.99 ON / 2780 | >20, stopped | >20, stopped | 4,523 / 9,392 before stopping | unmeasured |
| g50 ss0.50 ON / 3120 | 7.5394 | 2.9236 | 11,847 / 7,449 | 2.6e-14 |

The latter object has only three reads: the cost is not confined to deep objects.
The compact stranded probe spends 19.8345 of 20.0026 wall seconds inside the existing
native inner function; the compact low-count probe spends 2.8167 of 2.9236 there.
The machine reports high load, so these are bounded diagnostic costs, not release
throughput or a general speed-up claim. CPU and wall times agree closely in the probes
that record both. Main source, tests, goldens and both native binaries remain unchanged.

The compact stranded probe also reports a SciPy integration warning. A separate census
finds 75, 47 and 57 tilt intervals with no representable midpoint in the compact stranded,
unstranded and zero-control representatives, respectively. Their summed widths are about
1e-14 radians. Such intervals cannot be refined by subdivision; this observation does
not yet attribute every warning to that cause. The reference currently requests relative
accuracy separately for every interval, including negligible mass. The subsequent shared
budget experiment below addresses that numerical schedule. Do not drop small intervals,
merge nearby knots with an arbitrary tolerance, select a fixed tilt mesh or port the
expensive loop unchanged. This is numerical work on the same integral, independent of
the range proposal and the separate landscape-weight contrast.
Instrument: `.cache/rigel_runs/2026-10-08_both_cost/` (`census.json`, saved object tables,
`time_one.py`, scalar timing receipts, `compact.py`, `compact_before.xml`,
`compact_candidate.xml`, `mutations.json`, `table_audit.json`, `verification.json`).

BOTH-STRAND SHARED ERROR BUDGET (2026-10-08). This first shared-budget variant is
superseded by the cancellation check below; its receipts remain the comparison control.
The bound in **Tilt variation can be
bounded without fitting RNA** separates observed-column variation from the finite
profiles' slopes and ranges. It needs no inferred RNA amount, new prior or biological
constant. Across all 179 saved intervals without a representable midpoint, native
conditional log likelihoods change by at most 3.8e-12 and stay within the bound plus
the nominal inner integration allowance. This supports an outer-scheduling diagnosis;
it does not attribute every earlier integration warning to the same cause.

A single global SciPy call is rejected by a narrow Gaussian peak: it misses 15.73% of
the mass while reporting about 1.1e-14 relative error. The retained prototype instead
uses adaptive Simpson integration with the existing profile-derived breakpoints and
one error budget over all intervals. Endpoint samples expose the missed boundary
layer. At floating-point resolution, an interval's mass and derived uncertainty remain
in the sum; the prototype rejects an unresolved peak when the total cannot meet tolerance.
The Simpson estimate is still an error estimate, not a global convergence theorem.

The zero-variation bound control fails 10/14 initial checks. Further small-probability
and unstranded-identity falsifications expose rounding defects before repair. The final
25 new checks and the prior 37 pass together (62/62); fourteen actual mutations
collectively exercise all 25 new gates and are restored. These include a tiny interval
that dominates the integral, a negligible but uncertain interval, normalization, strand
mirror, nonmonotone profiles and loss of precision during error summation.

Fresh serial scalar pairs hold the saved inputs, isolated native binary and requested
1e-8 tolerance fixed. Both arms use the exact compactor:

| test condition / slot | separate-budget wall seconds | shared-budget wall seconds | conditional calls, separate / shared | absolute log-evidence difference |
|---|---:|---:|---:|---:|
| g50 ss0.99 ON / 2780 | >20, stopped | 4.9452 | 10,786 before stop / 2,504 | unmeasured |
| g50 ss0.50 ON / 3120 | 2.8500 | 0.7582 | 7,449 / 1,708 | 8.2e-11 |
| g00 ss0.99 ON / 2780 | 11.5175 | 3.4319 | 11,677 / 3,104 | 1.1e-10 |

Refining the tolerance from 1e-8 to 1e-10 at zero, half the observed-total rate and
twice that rate gives nine paired density points across those three saved objects.
All pass the nominal error comparison; the largest log-evidence change is 5.4e-10.
The three screens take 66.2, 14.0 and 60.8 seconds. This is a sparse curve/refinement
check, not a complete production-grid audit or independent truth for these real-valued
inherited factors.

These remain single-object cost probes, not calibration or whole-genome speed claims.
Native inner work still dominates. Seconds per density value remain too expensive for
the production consumer; reduce unnecessary conditional evaluations under the same
error target before choosing a native outer implementation. Keep factor support,
hypothesis weights, witness bits and all statistical inputs fixed. No count solve,
population fit, capture-weight calculation or release benchmark has changed.
Instrument: `.cache/rigel_runs/2026-10-08_tilt_budget/` (`bounds.py`, `tilt_integral.py`,
`budget_evidence.py`, failing controls, `reference.xml`, `mutation_coverage.json`,
refinement and verification receipts); fresh paired timings are the `budget_pair_*`
receipts in `.cache/rigel_runs/2026-10-08_both_cost/`.

BOTH-STRAND REFINEMENT CANCELLATION (2026-10-08). An audit finds that the native inner
integral already checks two successive Simpson differences, while the new outer
integral checked only one. A positive polynomial with an exactly integrated mass
exposes the consequence: the first outer implementation returns half the mass with a
zero error estimate. Coordinate rescaling and log normalization preserve the defect.
All three falsifications fail before the outer integral reuses the inner solver's
second-refinement check; all three then pass. Two actual removals of the check each
exercise every new gate. Together with the prior numerical and representation gates,
the strengthened implementation passes 65/65. The derivation and limits are **One
embedded difference can cancel** in EQUATIONS. This does not prove arbitrary-function
convergence and does not change the likelihood model.

Fresh serial saved-object pairs isolate this check, with the same native binary,
factors and 1e-8 requested tolerance:

| test condition / slot | one-check wall seconds | two-check wall seconds | conditional calls, one / two |
|---|---:|---:|---:|
| g50 ss0.99 ON / 2780 | 4.9443 | 8.1638 | 2,504 / 4,088 |
| g50 ss0.50 ON / 3120 | 0.7984 | 1.1659 | 1,708 / 2,752 |
| g00 ss0.99 ON / 2780 | 3.6477 | 6.3977 | 3,104 / 5,202 |

Nine sparse density-curve points agree with the frozen 1e-10 reference values within
6.0e-11 log likelihood. That agreement and the existing independent-integration gates
are numerical evidence, not new count or transcript accuracy results. Retain the extra
check when optimizing cost; the faster first variant is not the active candidate.
Instrument: `.cache/rigel_runs/2026-10-08_outer_refinement/` (`test_aliasing.py`,
`guarded_integral.py`, `run_guarded.py`, `before.xml`, `reference.xml`,
`mutation_coverage.json`, paired screen receipts and verification).

INNER-EVALUATION COST CENSUS (2026-10-08). A read-only table audit finds occasional
rounding-induced slopes above the sum of the original profiles' maximum slopes. Maximum
ratios are 1.58, 4.01 and 1.22 on the stranded, unstranded and zero-control representatives;
median ratios are below one. This measures a representation effect, not its runtime cost
or permission to discard near-coincident knots.

A separate standalone C++ diagnostic records each inner integrand call's RNA coordinate
and normalization anchor. On twelve conditional inputs drawn from the saved tables,
32.1–49.7% of calls repeat an already evaluated pair. Instrumented and uninstrumented
standalone results are bit-identical on all twelve inputs. The recursive refinement
computes the next subdivision to check error, then recomputes those same values when it
descends. The follow-up below isolates reuse from tolerance allocation and statistical
approximation. No speed-up is claimed from this census alone. Native library binaries
and the worktree's original header remained unchanged during this diagnostic.
Instrument: `.cache/rigel_runs/2026-10-08_inner_cost/` (`tables.json`, `prepare.py`,
`input_metadata.json`, `inputs.txt`, `plain.txt`, `counted.txt`, `counts.json`);
standalone C++ source and binaries are under `.cache/density_cost/` in the isolated
native worktree.

INNER-SAMPLE REUSE (2026-10-08). The isolated native integral now passes the quarter-point
values already evaluated by its error check into recursive children. Both error checks,
the subdivision decisions, budgets, tail bound and summation expressions are unchanged.
Twelve falsifications require fewer calls, bit-identical results, the same unique sampled
coordinates and a byte-identical refinement-decision trace. All twelve fail before this
change and pass afterwards. Two compiled mutations (restore redundant evaluation; count
one child's mass twice) each fail all twelve gates. Calls fall by 21.4–32.9% on these
conditional inputs.

The independent numerical and delivery checks pass 138/138 against the rebuilt library.
Two frozen reference families must run in separate processes: both import a module named
`reference`, so the initial combined invocation reused the wrong oracle and failed eight
signature checks. That receipt is retained; neither frozen test family was edited.
Fresh serial scalar pairs use the guarded outer integral and the same factors/tolerance:

| test condition / slot | previous wall seconds | reuse wall seconds | conditional calls, unchanged |
|---|---:|---:|---:|
| g50 ss0.99 ON / 2780 | 8.2243 | 5.6205 | 4,088 |
| g50 ss0.50 ON / 3120 | 1.1099 | 0.8396 | 2,752 |
| g00 ss0.99 ON / 2780 | 5.8905 | 4.0981 | 5,202 |

The three output values are bit-identical. These scalar timings are 24–32% lower, but
remain unsuitable for production throughput. Fresh cached calibration pairs on seven
conditions also preserve every region and boundary DNA count exactly. This confirms
count-path isolation; it is not a new calibration or transcript accuracy result.
The frozen control library remains `1ddc6e63`; the reuse library is `05171a89` in its
separate diagnostic site. Main source, tests, goldens and installed native library are
unchanged. Measure remaining inner work against each interval's contribution before
changing the allocation of numerical error; do not infer a need for a biological cutoff.
Instrument: `.cache/rigel_runs/2026-10-08_inner_reuse/` (`prepare.py`, `test_reuse.py`,
`before.xml`, `after.xml`, `mutation_coverage.json`, `single.xml`, `both.xml`,
`timing_summary.json`, `counts_control.json`, `counts_candidate.json`, verification);
standalone C++ and its traces are under `.cache/density_reuse/` in the isolated worktree.

INNER TOTAL-ERROR BUDGET (2026-10-08). Interval instrumentation preserves all twelve
conditional values and call counts. On the eight stranded and zero-control inputs,
72.9–88.3% of evaluations in the final integration pass serve intervals whose combined
mass is below the requested total numerical error. Separate relative budgets therefore
explain much of the remaining work. This is a numerical cost finding, not permission to
discard low-density evidence or objects.

The first replacement distributed absolute error by interval width. Its thirteen analytic
cost checks passed, but an existing independent narrow-profile test exhausted floating-point
subdivision. That variant is rejected and its source, native binary and failed receipt are
retained. The corrected scheduler keeps all intervals in a heap and refines the largest
estimated error until the summed error meets the total budget. It re-sums before acceptance
and retains both Simpson checks, sample reuse, all knots and stationary maxima, and the
analytic tail bound. The derivation is **A total integral has a total numerical budget**
in EQUATIONS. There is no new biological assumption or numerical tolerance.

The thirteen cost gates fail against the frozen relative-budget control and pass against
the correction. Two isolated components of the narrow-profile failure also pass. Three
compiled challenges cover all fifteen gates: restore relative budgets (13 failures),
restore width budgets (2), ignore estimated error (15). The rebuilt library passes
141/141 targeted numerical and delivery checks, including the independent both-strand
reference. Seven fresh cached calibrations preserve every region and boundary DNA count
against the frozen reuse control from the same session.

Fresh serial scalar pairs isolate the allocation change, keeping the guarded outer reader,
factors and 1e-8 requested tolerance fixed:

| test condition / slot | separate relative budgets, seconds | total budget, seconds | conditional calls, unchanged |
|---|---:|---:|---:|
| g50 ss0.99 ON / 2780 | 5.6327 | 0.5561 | 4,088 |
| g50 ss0.50 ON / 3120 | 0.8400 | 0.4491 | 2,752 |
| g00 ss0.99 ON / 2780 | 4.1419 | 0.5751 | 5,202 |

These are scalar numerical timings, not whole-library throughput or new calibration
accuracy. Tightening tolerance from 1e-8 to 1e-10 on nine saved density points changes
log likelihood by at most 3.8e-11; the tighter values differ from the prior frozen tight
reference by at most 3.7e-12. This is a sparse curve check, not a complete production grid.
The active isolated library is `636debae`; `05171a89` remains its frozen control.
Production source, tests, goldens and installed native library are unchanged.
Instrument: `.cache/rigel_runs/2026-10-08_inner_budget_cost/` (read-only interval census)
and `.cache/rigel_runs/2026-10-08_inner_budget/` (analytic inputs, conditional probes,
rejected width allocation, compiled mutations, independent checks, timings, count identity
and refinement receipts).

COMPLETE DENSITY-CURVE COST (2026-10-08). The saved stranded ON object 2780 completes
all 260 points of the existing frozen landscape grid. With the shared inner budget,
the guarded outer reader takes 143.44 seconds at 1e-8 and 434.32 seconds at 1e-10.
The maximum log-likelihood difference is 2.92e-10, at density 7.376. This validates
the numerical calculation across this complete fixed grid, not the adequacy of
that grid's support or the source factors. The nominal curve requires 1,062,258
conditional calls, so the consumer is still unsuitable for production throughput.

A cheap piecewise-linear compression screen at log-height tolerances from zero to
1e-8 finds no further knot reduction on the six representative RNA profiles beyond
the existing constant-segment compaction. This was an uncertified cost screen; no
approximate compression was implemented and no tolerance was selected from panel error.

The next prototype applies the existing uniform tilt envelopes to ordinary intervals
before requesting their full Simpson stencil. It retains every interval's mass and
uncertainty and falls back to the existing two-refinement check. Six analytic
constant/background/peak cost gates fail before implementation and pass afterward.
Its first scalar screen exposes a separate accounting defect: stranded ON exceeds
20 seconds because subtraction of large initial error bounds leaves a false residual.
Four of six additional analytic loose-bound cases reproduce the excess work. Direct
sums over live intervals replace the incremental bookkeeping, and all twelve focused
checks then pass. Four actual mutations cover all twelve gates. The independent
both-strand, bound, cancellation and new checks pass 71/71.

Fresh serial scalar pairs use the same native library and factors:

| test condition / slot | guarded outer, seconds | envelope-first outer, seconds | conditional calls, guarded / envelope-first |
|---|---:|---:|---:|
| g50 ss0.99 ON / 2780 | 0.5581 | 0.3306 | 4,088 / 2,162 |
| g50 ss0.50 ON / 3120 | 0.4491 | 0.3456 | 2,752 / 1,918 |
| g00 ss0.99 ON / 2780 | 0.5732 | 0.4800 | 5,202 / 3,876 |

The corrected outer reader also completes the same 260-point curve in 81.81 seconds;
every value agrees with the frozen tight reference within 2.92e-10 log likelihood.
That full-curve timing is an absolute feasibility measurement, not a fresh paired
whole-curve speed-up claim. The initial stalled implementation and its receipts are
retained. The derivation is **Envelope-first quadrature retains mass** in EQUATIONS.
No native, count, opportunity, prior or capture rule changes in this contrast.
The retained diagnostic runner is now `run_lazy.py`, with native `636debae`;
the guarded runner remains the comparison control. These timings do not support
whole-genome use. Complete source-factor and support certification before another consumer
rewrite; numerical agreement does not certify the inherited evidence model.
Instrument: `.cache/rigel_runs/2026-10-08_complete_curves/` (full grids and compression
screen), `.cache/rigel_runs/2026-10-08_lazy_tilt/` (evaluator, falsifications, mutation
coverage, stalled variant, independent checks, timings and verification).

BLURRED-LIKELIHOOD TAILS (2026-10-08). The source-factor audit confirms a separate
arithmetic defect. `blur_row` converts globally normalized log values to probabilities,
convolves, then floors each result at `1e-300`. Distinct likelihood tails consequently
collapse near -690. A direct evaluation of the same finite discrete Gaussian disagrees
by 6,208 nats on a quadratic log row with curvature 100. Certified Poisson-source tests
also fail: counts 100 and 500 lose their low-density tail. Ten of seventeen focused
tests fail before the correction. These are numerical counterexamples, not transcript
accuracy measurements or proof of the remaining RNA-short error's cause.

The isolated native correction retains the ordinary convolution for sums in the normal
floating-point range and uses a local log-sum-exp elsewhere, with log kernel weights.
No blur width, finite radius, edge padding, lattice, message schedule or statistical
rule changes. The derivation is **A convolution preserves relative likelihood in its
tails** in EQUATIONS. All seventeen gates and the seven earlier source-footprint gates
pass. Four compiled defects collectively trigger every new gate: restore the floor,
omit the blur, invent nonzero width at zero variance, or drop tail kernel weights.

A fresh seven-condition calibration A/B reads the existing synthetic scan caches and
scores against per-object oracle truth. The control reproduces its frozen count arrays
exactly. Six candidate conditions are also bit-identical. RNA-long ss0.99 ON changes
at most 0.000106 of a region count and 0.001535 of a boundary count; total absolute
movement is 0.000124 and 0.002575 respectively. Region error is
1.4657218062% -> 1.4657218041%, boundary error 3.3480121193% -> 3.3480121346%.
The zero control and RNA-short stranded OFF are among the identical cases. Retain this
as a numerical repair; it does not solve the outstanding count or capture errors and
does not justify repeating the full release panels. Main source and the installed
library are preserved. The existing message tests retain the identical four old-Gaussian
reference failures in both arms (47 passed / 4 failed each); they were not rewritten.
Ruff, the eleven documentation checks and diff checks pass.
Instrument: `.cache/rigel_runs/2026-10-08_blur_tails/`.

BLUR-TAIL INTEGRATION (2026-10-09). The arithmetic correction is applied independently
after the certified-source footprint repair. It retains ordinary convolution in the
normal floating-point range and evaluates smaller sums in log space. The original
spacing and finite Gaussian operator are unchanged. An additional gate preserves a
window with exactly zero likelihood: its log value remains negative infinity, rather
than becoming a NaN or a finite floor. Eleven of eighteen cases fail against the
integrated control; all pass after correction. Five actual compiled defects collectively
exercise all eighteen gates. The derivation is **A convolution preserves relative
likelihood in its tails** in EQUATIONS.

Six of seven fresh calibration conditions are bit-identical in all 29 result fields
and both belief arrays, including RNA-short stranded OFF and the zero-DNA control.
RNA-long stranded ON changes at most 0.000106 region DNA fragments and 0.001535 boundary
DNA fragments. Every class's calibration-error change is below 0.000000072 percentage
points; nonfinite uncertainty values and their positions are preserved.

The changed condition was then run end to end in ordinary isolated packages, with
scan=1, fractional assignment, seed zero and default calibration/EM threads. Values
are transcript / gene error (%); the archived release's truth hashes are rechecked.

| Condition | 0.7.1 | Before arithmetic repair | Arithmetic repair |
|---|---:|---:|---:|
| RNA-long, g50 ss0.99 ON | 7.3766 / 1.6179 | 4.2591 / 0.6831 | 4.2522 / 0.6832 |

The RNA pool gains 7.50 fragments, the DNA pool gains 0.0623 and the unspliced-RNA pool
loses 7.56. These small changes are not a reason to select or tune the arithmetic.
The serial LBX0190 comparison is bit-identical across the complete 254,319-row
transcript table, saved region/boundary DNA and RNA counts, pools and calibration
scalars. The finite-support problem and the detector-free reader remain open.

The installed source and native binary match the isolated candidate. The final main
suite passes 2,814/2,814 with no skips or xfails; lint, format checks and full preflight
pass. All 147 golden files retain their starting hashes. No commit or push was made.

Receipts: `.cache/rigel_runs/2026-10-09_blur_tails/` (ordinary frozen packages,
`before.xml`, `candidate.xml`, `mutations.json`, paired `screen` outputs,
`quant_comparison.json`, `real/comparison.json`, `main_suite.xml`, `installed.json`,
`verification.json`).

EDGE-PRODUCER COORDINATE AUDIT (2026-10-08). The earlier held-level counterexample also
occurs before delivery: the native EDGE builder converts a pure-DNA boundary's Poisson
lower bound through the recipient's observed total and publishes it as composition.
Reading this row at an expected DNA/RNA ratio evaluates the source at the wrong density.
The derivation is **An edge's count projection is not its intensity coordinate** in
EQUATIONS. This indicts that input to the new intensity reader; it does not establish
a corresponding defect in the existing count solver.

With source count 20, unit opportunities and target columns 4/4, the current producer's
factor at fixed expected DNA 4 varies from -7.0678 to -32.1147 nats as expected RNA
varies from 0.25 to 16. The unchanged source-density factor stays at -16.1888. The
complete target density curve has sampled mode 8.2 through the composition projection
and 13.8 through the source density, agreeing with the earlier generic counterexample.
No capture assumption, source width or population prior changes in this comparison.

Frozen native-message censuses identify 249 active direct EDGE messages in test g50
ss0.99 ON, 243 in g50 ss0.50 ON, and zero in g00 ss0.99 ON. Each active source is
structurally pure DNA; none of these receiving sides already carries a DNA-level
message. Thus final-channel separation cannot recover these discarded coordinates,
and this particular defect cannot explain the measured zero-control false capture.

A small Python attribution adapter removes only a directly received EDGE contribution
from composition and evaluates the same one-sided source factor at candidate DNA
density. Other composition terms, DNA/RNA factors, presence bits and side exclusions
are retained. No new message pass or count call site is added. Ten of thirteen gates
fail against the unchanged interface; all thirteen pass with the adapter. Four actual
reversions/mutations exercise all thirteen, including source duplication and invented
evidence where no direct edge is delivered.

On preselected median-index, single-RNA-strand recipients, the existing 260-point grids
change as follows. These are input-curve attribution measurements, not oracle accuracy:

| test condition / slot | old / corrected sampled density mode | maximum normalized-height change |
|---|---:|---:|
| g50 ss0.99 ON / 1282 | 0.003758 / 0.004892 | 0.34190 |
| g50 ss0.50 ON / 1278 | 1.803805 / 1.803805 | 0.03970 |

The original count-facing delivery remains unchanged. A further face ablation, with
all other prebuilt rules and lane licences fixed, changes only the direct recipient in
the stranded example. In the unstranded example it removes eight composition messages:
the direct edge at slot 1278 and seven forwarded/transported copies through slot 1271.
Consequently the direct-recipient adapter is not a general repair and must not be
installed as one. The subsequent capture-order check also rejects carrying this source
as an absolute density through the existing composition routes: that would change their
meaning and reopen **The paradigm** in DESIGN. The transfer assumption itself now needs
review (`ISSUES: the-edge-density-floor-under-capture`). Do not start a population refit
on only the directly corrected inputs.
Instrument: `.cache/rigel_runs/2026-10-08_edge_source/` (`audit.py`, `edge_sources.py`,
falsification/mutation receipts, `replay.py` and `trace.py`). The first audit's offset-only
comparison failure is retained; normalized native delivery is reconstructed exactly.

MESSAGE-LOCALITY AUDIT (2026-10-08). Rebuilding the existing messages establishes two
different dependencies on disconnected exon expression, with local observations,
opportunities, topology, strand fit and incoming beliefs fixed. These checks concern
message construction; a population landscape is still deliberately shared information.

First, `TransferPolicy.library` makes `split_live` require both the measured live strand
protocol and a counted single-strand exon somewhere in the library. Adding one remote
exon with zero RNA opportunity toggles this bit while leaving both coordinate origins
unchanged. It switches each RNA hop between total-column and strand-difference discrepancy
pricing. A three-object interface counterexample changes the normalized DNA-evidence
curve by up to 0.500 in relative height. A separate census of annotation-derived triples
finds one affected triple among 36 eligible in test g50 ss0.99 ON, and none among 13 in
its zero control. At chain slots 3052–3054, the RNA+ bound delivered to the both-strand
exon changes by 0.250 nats, or 0.01565 of normalized height. These are isolated local
interventions, not whole-panel count errors; the intact libraries already have counted
single-strand exons and do not toggle the bit.

The scratch repair removes only that extra expression condition: an available strand
model and the existing protocol verdict determine witness availability. It preserves
the two hop formulas, all observations, faces, coordinate reductions and the declared
strand decision. Four of seven falsification tests fail first on the existing code;
the repair passes all seven, and four actual implementation mutations exercise every
gate. The annotation-derived bound and hand-built evidence become bit-identical under
the remote-exon intervention. This restores the stated protocol rule; it does not certify
the existing discrepancy approximation, especially at objects admitting both RNA strands.

Seven cached calibration A/Bs cover RNA-long unstranded/stranded ON, ladder g98 unstranded
OFF, RNA-short stranded OFF, and test zero/stranded/unstranded ON. All 35 library-input
comparisons are identical, as are the final region/boundary DNA-count arrays. These are
numerical no-op checks on those conditions, not new transcript/gene accuracy measurements.
These were prototype receipts; production integration is recorded next.

LOCAL-WITNESS INTEGRATION (2026-10-09). After count-foundation integration, the ordinary
installed package reproduces the disconnected-expression dependency. The repair removes
only the redundant counted-exon condition, implementing **Local strand-witness availability**
in DESIGN and EQUATIONS. A fresh portable gate covers both protocol orientations, unavailable
protocol/model guards, received RNA rows, and final single-/both-strand delivery. Six of ten
cases fail before the change; all ten pass afterward. Five actual Python implementation
mutations collectively exercise all ten gates. The 35 existing transfer-policy and RNA-lane
checks also pass. Native source and binaries are unchanged.

Fresh ordinary-package calibration pairs on the same seven conditions preserve all 29 result
fields and both belief arrays exactly. Origin truth, slot identity and observed counts are
checked afresh; class-specific count errors are unchanged. This is a sparse-input robustness
repair, not a new transcript/gene accuracy improvement. The disconnected-coordinate audit also
reproduces the separate zero-origin and finite-support losses on the integrated foundation;
this one-condition removal does not repair those. No capture prior or reader is selected.

The serial LBX0190 end-to-end pair is bit-identical across every payload/calibration array and
the 254,319-row transcript digest. The first full suite identifies one additional golden move:
`antisense_overlap_ss90` has 1,000 true RNA fragments and no DNA. Its false DNA decreases from
2.145629 to 2.063161 fragments; the largest transcript change is 0.046568 fragments. Transcript
and gene error both move from 0.214563% to 0.206316%. Every numeric field is reviewed, with
schemas, dtypes and text unchanged. Only that fixture's seven golden files are regenerated;
no tolerance changes. The completed suite derives **2,762 passed, zero failed/skipped/xfail**;
lint, formatting and full preflight pass. The failed pre-update receipt is retained.

Receipts: `.cache/rigel_runs/2026-10-09_reader_scope/` (`audit.py`, `audit.json`, the failed-first
gate, frozen ordinary packages, `mutations.json`, seven paired `screen` outputs and truth scores,
the real identity pair, `golden_review.json`, `main_suite.xml` and `verification.json`).

Second, the RNA coordinate can be zero despite positive local RNA-source evidence.
The exon fallback reads opportunity rather than positive count: an empty single-strand
exon can suppress the positive both-strand-exon fallback. With no exon counts at all,
counted introns can still publish strand claims, but `rna_lane` drops every source when
the reference is zero. These are executable missing-source counterexamples, replacing
the earlier unmeasured fallback item in the latent-defect inventory.

A separate Python prototype now repairs only the zero-reference path. It forms a positive
numerical coordinate from available RNA-admitting count/opportunity pairs and certified
flux pairs `(count, count/route_rate)`. The original positive reference, DNA coordinate,
strand-witness bit and native source rules are preserved. This does not fit a library RNA
abundance or feed summed counts to a new likelihood. The control fails 12 of 17 new gates;
the repair passes all 17, including flux-only empty pieces on both strands at either
junction end. Six deliberate defects exercise every gate. Both fresh control and repaired
builds pass the existing 17 RNA-lane tests. Seven saved panel contexts produce identical
library inputs with no array mutation; there is no new end-to-end performance claim.
The original two zero-coordinate locality specifications now pass. The separate finite-
support specification remains failing. Main source, goldens and installed native code
remain unchanged. Instrument: `.cache/rigel_runs/2026-10-08_rna_coordinate/`
(`positive_coordinate.py`, `before.xml`, `candidate.xml`, `mutations.json`,
`existing_control.xml`, `existing_candidate.xml`, `existing_locality.xml`, `inputs.json`).

ZERO-COORDINATE INTEGRATION (2026-10-09). The correction above is applied independently
on top of the local-witness repair. The existing positive RNA origin is returned unchanged;
the zero path uses the already-derived positive count/exposure scale. The DNA origin, strand
protocol verdict, native builders, grid, hop prices and source admission remain unchanged.
This implements **Changing a coordinate preserves physical support** in EQUATIONS, without
the proposed independent-range representation. No capture prior, cutoff or population model
is introduced.

The portable source/unit gates extend the archived cases to unstranded and absent-strand-model
libraries at both junction ends on both RNA strands. Twenty of 25 cases fail against the
integrated foundation; all 25 pass after correction. Six actual implementation defects
collectively exercise every gate. A full frozen-package suite passes 2,787 tests, including
the added gates, with no golden movement. Seven fresh cached calibration pairs preserve all
29 result fields and both belief arrays; slot and observation identity are checked against
the oracle before scoring regions, boundaries, introns and exons separately. All those
count errors are unchanged. This restores missing-source behavior in the explicit limit
cases, not an observed accuracy gain on the panels. The finite-support failure remains open.

The integrated whole-library LBX0190 pipeline is bit-identical: 19 payload arrays, 20
calibration arrays and all 254,319 transcript rows. The first main-suite run exposed the
already-recorded scheduler-dependent reorder gate; its deterministic replacement is recorded
under **hygiene-ledger**. The final main suite passes 2,788/2,788 with no skips or xfails;
lint, format checks and full preflight pass. Every golden file is unchanged from the start
of this repair. The installed native binary is unchanged, and the integrated Python source
matches the frozen candidate used for the paired measurements.

Receipts: `.cache/rigel_runs/2026-10-09_rna_coordinate/` (ordinary frozen packages,
`before_expanded.xml`, `candidate.xml`, `mutations.json`, `worktree_suite.xml`, paired
`screen` outputs and class-specific truth scores, `real_reference.json`, `real_candidate.log`,
`main_suite.xml`, `verification.json`).

Finite support is a distinct failure. With the default half-window 10, changing only the
remote RNA coordinate by `exp(14)` moves the local source below the level grid and removes
its delivered RNA factor. Halving the step retains the failure. Expanding the half-window
to 24 in a diagnostic control restores the coordinate contrast to below 4.5e-13 nats at
both spacings. Neither value is a proposed production constant. The algebra belongs to
`EQUATIONS.md`, **Changing a coordinate preserves physical support**. A mode at the lower
grid edge is not alone a truncation diagnosis: a zero-component mode may belong there.
The seven real-input censuses do not establish a current panel loss from this mechanism.

A further isolated source-table error is established with the mode still inside the
grid. A certified count of one at exposure 1,000, reference one and its existing empty-
recipient hop variance gives 0.08679 maximum relative-height error under endpoint-padded
blurring at half-window 10. Direct evaluation of the same finite convolution agrees at
half-windows 20 and 40. These are diagnostic windows, not selected constants. The
source has a known Poisson likelihood beyond the table; repeating its endpoint is not
that likelihood. The operator identity belongs to `EQUATIONS.md`, **A known source
supplies the convolution footprint**.

The isolated native repair evaluates the known source over the existing blur's halo
and crops back to the original grid. It preserves the width, lower-side operation,
zero-width path and all inference grids. Two of seven new gates fail on the control;
all seven pass after repair, and four compiled implementation defects exercise every
gate. The existing 17 RNA-lane gates pass. A fresh seven-condition calibration pair
reproduces the frozen count candidate exactly in its control. The largest changed
per-object DNA count is 0.05911 fragments; region/boundary calibration error changes
are below 0.000011 percentage points. No transcript/gene improvement is claimed.
The original finite-support locality specification still fails after this repair and
the zero-coordinate fallback. Main source and the installed native library are unchanged.
Instrument: `.cache/rigel_runs/2026-10-08_level_support/` (`audit.json`,
`test_source_blur.py`, before/repaired JUnit receipts, `mutations.json`, paired screens
and `existing_locality.xml`).

SOURCE-FOOTPRINT INTEGRATION (2026-10-09). The isolated certified-flux repair is applied
on top of the integrated count foundation and locality corrections. The source is
evaluated over the existing Gaussian footprint before cropping to the retained table;
the grid, source admission, one-sided readout and hop variance remain unchanged. This
implements **A known source supplies the convolution footprint** in EQUATIONS. A fresh
check found a rounding defect in the archived prototype: subtracting extended coordinates
to reconstruct spacing can add a kernel cell. Passing the original spacing through the
same blur corrects it. Two source checks fail against the control, the rounding check
fails against the first prototype, and five actual compiled defects cover all eight gates.

Seven fresh calibration pairs preserve slot/observation identity against truth. The
largest changed DNA count is 0.059104 fragments; every class's error changes by less
than 0.000012 percentage points. These are calibration errors, not transcript errors.
The full candidate suite exposes one golden move: `extreme_abundance_ratio` shifts
0.00000056635 fragments from DNA to RNA. All numeric columns, schemas and text were
reviewed; transcript error is 0.455115648 → 0.455115660%, gene error is
0.001595802 → 0.001595745%, and true DNA is zero. No tolerance change is selected.

The ordinary test-chromosome pipelines use scan=1, fractional assignment, seed zero
and default calibration/EM threads. Each cell is transcript / gene error (%), at g50:

| Stratum | 0.7.1 archived | Integrated foundation | Source-footprint repair |
|---|---:|---:|---:|
| Unstranded OFF | 10.23 / 2.14 | 10.33 / 1.54 | 10.50 / 1.54 |
| Stranded OFF | 10.27 / 1.43 | 9.28 / 1.16 | 9.28 / 1.16 |
| Stranded ON | 12.97 / 4.07 | 10.17 / 1.20 | 10.33 / 1.20 |
| Unstranded ON, deferred | 27.38 / 18.98 | 12.45 / 4.42 | 12.97 / 4.42 |

The stranded-ON downstream movement is concentrated in three ambiguous isoform-cluster
genes (2,457.52 of 2,526.32 changed transcript fragments). Under the existing uniform
start, the contrast reverses: 9.8356 → 9.7458% transcript error, with genes 1.26221 →
1.26216%. Both starts have bit-identical calibration to their respective cached arm.
This supports downstream sensitivity; it does not establish equal objectives or a flat
likelihood, nor justify selecting a different EM start. The shipped start remains under
**the-em-answer-depends-on-where-it-starts**. No transcript-accuracy gain is claimed.

The serial LBX0190 pair has matching inputs and settings. Across 254,319 transcripts,
absolute count change totals 0.000001334 fragments (1.0544e-9%); the largest changed
calibration DNA count is 1.71e-11 fragments. This is output stability, not real-data
accuracy. The pair took 15.57 / 15.64 seconds and does not establish a speed-up. The
finite-support and low-probability blur-arithmetic defects remain separate open work.

The integrated extension is byte-identical to the isolated candidate. The eight portable
gates pass against it, the final main suite passes 2,796/2,796 with no skips or xfails,
and lint, format checks and full preflight pass. Only the seven files for the reviewed
golden fixture changed; their values match the pre-update review exactly. No commit or
push was made.

Receipts: `.cache/rigel_runs/2026-10-09_flux_footprint/` (ordinary frozen packages,
`before.xml`, `spacing_before.xml`, `candidate.xml`, `mutations.json`, paired `screen`
and `quant` outputs, `golden_review.json`, `em_comparison.json`, `real/comparison.json`,
`main_suite.xml`, `golden_landing.json`, `installed.json`, `verification.json`).

SOURCE-RANGE LIMIT (2026-10-08). A subsequent Python diagnostic derives a shared range
from the existing RNA-fraction limits, Poisson source tails and each certified source's
actual blur footprint, at the configured spacing. No capture yield, fitted biological
threshold or simulator label enters this calculation. The unchanged range first fails
all 13 coordinate/unit/range specifications. Raw-source coverage passes those 13 but
fails 17 of 25 flux/end/strand/no-source checks; including the source blur footprint and
preserving an inert range passes all 38. These passes certify the source examples only.
On seven frozen condition inputs, this shared range increases the grid-point count by
2.07–4.22 times. That is an allocation-size estimate, not measured runtime.

Repeated propagation rejects this diagnostic as a general range repair. An unstranded
chain carries a certified RNA source through balanced low-count recipients, using the
unchanged source, hop and propagation rules. At the same 0.2 spacing, the source-derived
half-window 27 is compared with diagnostic half-windows 54, 108 and 216. All comparisons
read the same physical log-density interval, -20 to 0. The errors below are maximum
absolute differences in normalized row height, not transcript/gene or count errors.

| Counted receiving nodes | Half-window 27 versus 216 | Half-window 108 versus 216 |
|---:|---:|---:|
| 10 | 2.94350e-8 | 9.66338e-13 |
| 42 | 0.00105305 | 5.43565e-13 |
| 122 | 0.0362012 | 3.64153e-13 |
| 242 | 0.132593 | 4.72796e-8 |

The two longer chains fail the existing one-percent row target; the wider pair supports
the reference's numerical convergence. These synthetic operator checks use no capture
scenario or expected yield. No new fixed window is selected. No mutation campaign,
calibration screen or release benchmark was run for this rejected range proposal.
The dependency identity belongs to `EQUATIONS.md`, **Repeated blurs require the propagated
footprint**. Whether level tables should have numerical extents independent of the
composition/count interval is an owner representation decision, not yet an approved rule.
Instrument: `.cache/rigel_runs/2026-10-08_shared_support/` (`support.py`, `before.xml`,
`source_only.xml`, `flux_before.xml`, `source_covered.xml`, `range_census.json`,
`propagation.json`, `propagation_curves.npz`, `propagation_red.xml`).

POST-INTEGRATION RANGE RECHECK (2026-10-09). All three retained range specifications
still fail on the current ordinary production import after the witness, zero-coordinate,
source-footprint and blur-tail repairs. The disconnected-expression contrast keeps local
observations, opportunities, topology and the fitted strand protocol fixed. Its RNA
contribution spans 33.7588 nats before the coordinate moves and exactly zero afterward;
the native RNA conversion therefore drops that contribution. A DNA contribution remains,
so this is not a claim that the entire delivered message disappears. Halving the spacing
does not recover the RNA contribution; a wider diagnostic interval does. That establishes
the loss, not complete coordinate invariance of the wider control.

The four repeated-propagation comparisons reproduce the table above to its displayed
precision, including convergence of the widest pair. The source-range diagnostic is
byte-identical to the earlier `support.py`; only the harness interfaces were adapted to
the integrated code. These are operator checks, not accuracy measurements. No source,
production test, golden or native binary changed, and no independent density range was
implemented. Instrument: `.cache/rigel_runs/2026-10-09_range_recheck/` (`audit.json`,
`curves.npz`, source/native hashes and `range_red.xml`: three intentional failures).

INDEPENDENT-RANGE PROTOTYPE (2026-10-09; owner-authorized experiment). A Python reference
plans the needed intervals backward through the existing density operations and evaluates
them forward, preserving each source's normalization and limiting values. A separate C++
evaluator in the isolated worktree implements the same numerical operations. The count
lattice, source admissions, face permissions and hop formulas stay fixed. No production
count consumer is wired to it. The derivation is **A finite core certifies level
normalization** in EQUATIONS; authorization is **Density-range prototype authorization**
in DESIGN.

Both evaluators pass 70 specifications covering the earlier range failures, both component
opportunity-gap directions, both RNA strands and junction ends, empty relays, source
absence, conflicting sources, zero-density limits and even count-grid sizes. Seventeen
actual Python defects collectively exercise every specification; eight compiled evaluator
defects are also caught. The seven frozen panel contexts require 0.846–2.582 times the
cells in the currently present received-density rows, while retaining the count grid.
This is an allocation census, not total pipeline memory. Three cached test-chromosome
contexts agree between Python and C++ to about 1.1e-13 in normalized height. Their first
paired native evaluations take 0.003–0.187 seconds; these are numerical evaluator timings,
not pipeline or calibration timings. The reference adapter's repeated whole-chain total
reduction was replaced by one cached array; all seven census records are otherwise exact.

COORDINATE-PHASE LIMIT (2026-10-09). The independent-range prototype does not yet meet
arbitrary-coordinate invariance. A subsequent check shifts the RNA reference through one
0.2-wide cell while retaining the same local observations, opportunities, strand protocol
and count lattice. For three objects containing 100, 1,000 or 100,000 observations each,
the maximum normalized-message height changes are 0.02753, 0.05317 and 0.05370. All exceed
the existing approximately one-percent row target. Whole-cell shifts agree to rounding.
Moving the density knots along with the coordinate restores agreement; see
**Coordinate translation preserves physical knots** in EQUATIONS. These are message-shape
differences; transcript/gene and calibration consequences have not been measured.

Three new specifications therefore remain red for the range-only prototype outside the production suite. Do not call
the 70 passing range tests complete locality or density-interface certification. Do not
port this evaluator directly into calibration or launch release panels for it. A local,
unit-covariant numerical origin per connected part of each existing level lane was
subsequently authorized for prototyping; see **Local density-unit prototype authorization**
in DESIGN. It adds no biological edges or prior. It is not integrated.
Instrument: `.cache/rigel_runs/2026-10-09_density_support/` (Python reference and replay,
failed-first, Python/native/mutation receipts, census, timings, `phase.json` and
`phase_red.xml`). The standalone C++ evaluator is in the existing isolated count worktree
under `.cache/density_support/`.

LOCAL-UNIT PROTOTYPE (2026-10-09; owner-authorized experiment). The Python replay now
uses `sum(n)/sum(E)` from already admitted sources within each connected component of
each permitted density lane. Counts for a composition source are observed totals, not
estimated RNA or gDNA amounts; certified sources use their own count and implied
exposure. The ratio supplies only a numerical unit. Source rules, hop prices, graph,
count grid and density spacing remain fixed. The derivation is **Local units fix
coordinates, not quadrature error** in EQUATIONS.

All three failed-first coordinate specifications pass. The full isolated unit contract
has 40 passing cases, including opportunity-unit covariance, both RNA strands and junction
ends, source absence, source/storage/query order and disconnected-component packing.
The same 40 checks pass through the existing isolated C++ evaluator; unit selection
itself is Python, not a native policy implementation. Twelve actual Python defects
collectively exercise every new passing check, then are restored. No native source,
production source, production test or golden changes in this experiment.

SPACING SENSITIVITY (2026-10-09). Coordinate invariance does not certify the unchanged
density spacing. With the count lattice fixed, successively refining only the level
evaluation gives the following results on the same three-object controls. The finest
control uses 1/64 of the original density step; this is an instrument, not a proposed
production setting. The last two refinements differ far less than the existing
approximately one-percent row target.

| Observations per object | Coarse vs finest, normalized message-height difference | Last two refinements, height difference | Coarse minus finest, gDNA fraction (percentage points) |
|---|---:|---:|---:|
| 100 | 0.0163663 | 0.00000472 | +0.079412 |
| 1,000 | 0.0348723 | 0.00001712 | +0.019175 |
| 100,000 | 0.0333289 | 0.00031931 | approximately 0 |

The final column is a one-object sensitivity calculation through the existing native
count solver: own observations and every other term are fixed, with the usual symmetric
reference prior and no fitted landscape. It is not a calibration accuracy measurement,
does not exercise landscape feedback, and does not validate the future density reader.
Three new convergence checks remain intentionally red. They prevent reporting the
local-unit covariance result as complete numerical certification. A native unit-selector
port and larger benchmarks were withheld after this cheap screen. Do not select a fixed
refinement from these fixtures. First price the residual approximation in density and
count readouts on cached cases before proposing more numerical machinery or a changed
error contract. Production integration still requires its separate checkpoint.

Instrument: `.cache/rigel_runs/2026-10-09_density_units/` (`before.xml`, `reference.xml`,
`native_evaluator.xml`, mutation receipts, `convergence.json`, `convergence_red.xml`).

SPACING READOUT AUDIT (2026-10-09; owner requested practical effects before more numerical
complexity). Seven fresh cached calibrations supply the fitted landscape, observations,
opportunities, strand protocol and initial count beliefs. All are held fixed, as are the
count lattice, composition messages, admissions, faces, hop prices and propagation. Only
the density spacing changes, from the existing step to one quarter and one eighth. Both
arms use independent support and local units: this isolates spacing, not the full
representation change against production. There is no landscape refit or EM run.

The three test-chromosome conditions cover every counted RNA-admitting object. Each larger
condition selects evenly spaced depth ranks plus the largest-count objects within each
class/strand-admission combination, without truth-based selection. All required message
ancestors are retained. Larger-panel coverage is explicitly partial. The allocation shift
below is `100 * sum(M * abs(delta f_g)) / sum(observed unspliced counts)` over the checked
objects; region/boundary incidences are not unique fragments. Conditions are never pooled.

| Condition | Objects checked / eligible | Count coverage | Absolute allocation shift (%) | Object share shift, 99th percentile / maximum (percentage points) |
|---|---:|---:|---:|---:|
| test g50 ss0.99 ON | 1,976 / 1,976 | 100% | 0.02339 | 0.5766 / 3.1080 |
| test g00 ss0.99 ON | 870 / 870 | 100% | 0.00001 | 0.0003 / 0.0025 |
| test g50 ss0.50 ON, deferred | 1,977 / 1,977 | 100% | 0.00763 | 0.0305 / 0.2400 |
| RNA-short g50 ss0.99 OFF | 566 / 54,075 | 20.4% | 0.00024 | 0.0416 / 0.1282 |
| RNA-long g50 ss0.99 ON | 566 / 38,730 | 23.3% | 0.14501 | 0.4389 / 12.3654 |
| ladder g98 ss0.50 OFF | 566 / 51,494 | 15.2% | 0.00032 | 0.0331 / 0.1817 |
| RNA-long g50 ss0.50 ON, deferred | 566 / 38,841 | 23.4% | 0.00197 | 0.0089 / 0.1652 |

Two both-strand exons account for 93.5% of the sampled RNA-long stranded allocation shift.
Their source counts and own strand likelihoods stay fixed. The finer reference is checked
again at one sixteenth of the density step, on those objects alone; the last two count
readouts differ by 0.0400 and 0.0581 percentage points. No fixed refinement is selected.

| RNA-long object | Oracle gDNA fraction (%) | Production fraction (%) | Prototype, existing step (%) | Prototype, one-sixteenth step (%) |
|---|---:|---:|---:|---:|
| exon slot 39465 | 13.7293 | 13.5046 | 0.0387 | 11.2325 |
| exon slot 44511 | 13.8294 | 14.0196 | 0.0373 | 12.4608 |

Replacing only the delivered negative-RNA row at the first object or positive-RNA row at
the second reproduces essentially the entire coarse-to-fine change. Replacing the gDNA
rows alone does not. Thus the RNA-level numerical representation carries the effect;
this does not yet identify which source/operation in its ancestry causes it. Unlike the
test-chromosome outlier below, these two also have large posterior-mean changes. They
cannot be explained away solely as an unstable median. The current production counts
are close to truth here: this is a failure of the coarse prototype, not a newly measured
production regression. The sample does not estimate its prevalence in other data.

At test-chromosome exon 2962, in contrast, less than one percent total-variation distance
separates the posteriors, but the median moves 3.108 points between two modes. The mean
moves 0.378 points. This is a readout sensitivity under the same fitted prior, not grounds
to change the count model. A count-derived geometric-mean density can be more sensitive
than the allocation: the RNA-long maximum changes about 71-fold, and a zero-DNA test
object changes 1.625-fold while its allocated DNA changes by only 0.0013 incidence counts.
These moments include the existing count prior; they are not the proposed prior-free
density evidence or a tested capture readout. Small count effects alone do not certify
the unfinished reader.

Instrument validation first exposed an adapter defect: omitting the existing pure-strand
witness exclusions broke native-count identity. After correction, the independent
legacy-table replay agrees with production in all seven conditions to at most 2.5e-14
in fraction. Actually omitting that rule again makes every condition disagree. The entire
spacing audit was rerun; earlier provisional values are superseded. The separate direct
native replay is exact. The corrected seven-condition audit takes about four minutes;
the focused finer reference takes 44 seconds and the lane attribution 12 seconds.

This closes the proposed cheap count-readout measurement, not numerical certification or
release validation. Most checked allocations have small sensitivity, but uniformly
negligible effects are falsified by the two captured RNA-long exons. No refinement
parameter, adaptive integrator, production landing, golden update or release A/B follows
from this result. Keep the count foundation stable; retain the outliers for a bounded
RNA-message diagnosis and eventual consumer comparison. Detector-free capture, its prior
and opportunity validation remain separate release work. Instrument:
`.cache/rigel_runs/2026-10-09_spacing_effect/` (`freeze.json`, `assembly_control.json`,
`batch.json`, `report.json`, `outlier_reference.json`, `outlier_attribution.json`,
`verification.json` and the reproducible scripts).

NARROW-SOURCE INTERPOLATION (2026-10-09; follow-up authorized by the owner). Tracing only
the two implicated RNA rows requires 23 numerical ancestors. The sensitive local sources
are certified splice counts, 3,607 at RNA-long exon 39465 and 3,598 at exon 44511. Their
declared source-blur variances are 0.00089929 and 0.00034835 in log density, against table
spacing 0.199797. The functions are much narrower than that spacing. At the coarse knots,
the normalized heights of the two own lower-side curves agree with the fine control to
about 3.8e-10 and 5.1e-16. That apparently excellent nodal agreement misses the error where
the consumer actually evaluates the functions, between stored points.

A decisive contrast takes the *fine final RNA curves* and samples them back onto just
the coarse knots. All other factors stay at the coarse control. The gDNA fractions are:

| RNA-long exon slot | Original coarse curve (%) | Fine RNA curve (%) | Fine curve sampled onto coarse knots (%) |
|---|---:|---:|---:|
| 39465 | 0.03869 | 11.27251 | 0.03913 |
| 44511 | 0.03728 | 12.40258 | 0.04072 |

This isolates loss between the stored knots as sufficient to reproduce essentially all
of the count effect. It does not prove that every source convolution or upstream table
is exact. The own likelihoods, opportunities, counts, priors, admissions and discrepancy
prices are unchanged. The ordinary coarse replay reproduces the previous count receipt
to 1e-14. The trace takes about eight seconds and the resampling contrast about four;
no new calibration or EM is run.

Three independent Poisson-source specifications at 100, 1,000 and 100,000 observations
also fail the existing one-percent normalized-height target *between* knots, with errors
0.22449, 0.55797 and 0.60685. They use no panel truth or probe layout and remain red outside
the production suite. The derivation is **A narrow likelihood can disappear between
accurate knots** in EQUATIONS. A finer fixed default is not selected from these examples.

The proposed small correction preserves the known source curve through its existing blur
and lower-side operation at the consumer, instead of quantizing that direct measurement
onto a coarse intermediate row. It still needs a complete reference, native design and
falsification/mutation coverage; none is claimed implemented here. Keep the validated
count foundation frozen and make this source-precision issue an explicit input contract
of the capture reader. A broadly adaptive second message engine is not justified by this
diagnosis. The owner subsequently authorized the isolated local-prior prototype
recorded below; production selection remains separate. Instrument:
`.cache/rigel_runs/2026-10-09_rna_spacing/` (`trace.py`, `trace.json`, saved ancestor curves,
`interpolation.py`, `interpolation.json`, `interpolation_red.xml`, `verification.json`).

PROPER LOCAL CAPTURE PRIOR (2026-10-09; owner-authorized prototype). The scalar reader
consumes density evidence and a fixed normalized background, with an unbounded
`1/C²` enrichment prior and the continuous counterarm. Its derivation and limits
are **Proper local capture reference and continuous correction** in EQUATIONS.
No count solve, background fit, landscape, junction rule or production source is
changed. The proper tail and equal odds are separate contrasts; every contrast
uses the same likelihood. The old bounded reader first fails three substantive
gates: disconnected RNA changes its log weight by 0.19423, enlarging annotation
changes it by 2.06808, and an external ceiling truncates a strongly supported
capture factor. The candidate passes the mathematical limits and invariances;
deliberate mutations of its tail, counterarm, background, zero-count treatment,
units, ceiling and posterior probability are detected.

The own-evidence screen covers all 139 archived captured objects with 4–30 total
reads in the specified test conditions, plus the previously implicated zero-DNA
boundary. The zero-DNA group has 52 objects including that boundary; its maximum
weight is 1.04789 and the boundary reads 1.00069. This is an own-evidence diagnostic,
not a selected strand-only reader. In weakly stranded capture the own curves often
cannot identify enrichment: the 53 selected g25 ss0.70 objects have median weight
1.168, despite much higher expected capture. Suppressing ambiguity alone would
lose real capture; licensed neighbour evidence remains necessary.

Using current typed message inputs on four prespecified objects gives the following
weights. The two controls are **not shipped readers**: they apply the old bounded
slab or the proper tail to the same honest curve, both at the archived `1/n` odds.

| Condition / slot | Bounded slab, old odds | Proper tail, old odds | Proper tail, equal odds |
|---|---:|---:|---:|
| test g00 ss0.99 ON / boundary 513 | 1.134 | 1.00008 | 1.140 |
| test g50 ss0.99 ON / exon 1000 | 603.83 | 279.94 | 443.87 |
| test g50 ss0.99 ON / exon 1036 | 699.23 | 476.40 | 478.04 |
| test g50 ss0.50 ON / exon 2958, deferred diagnostic | 811.61 | 810.05 | 810.05 |

The captured stranded exons have 9 and 14 total reads and expected capture 787.33
from the simulator's expected yields, not realized DNA counts. The g00 object
tests amplification in the absence of DNA; its library physically has capture,
so weight one is a neutral-evidence limit, not its known chemical capture factor.
The residual 1.14 correction with messages is not claimed harmless end to end.
Background distributions are frozen archive inputs. The deferred object uses the
stranded condition's fixed background as a cost/reference control; this is not
an unstranded background fit or accuracy score. Message tables retain the known
support/interpolation approximations. Tightening outer tolerance from 1e-6 to
1e-8 changes log weight by at most 3.2e-14 on the three stranded objects; it does
not certify those input tables or total integration error.

Exact Poisson risk calculations also show finite-count costs. At background
expectation 0.01 and true capture 100, only one DNA read is expected and mean
reported weight is 27.05; at background expectation 1 and true capture 100,
mean weight is 98.50. An uncaptured pure-DNA object with expectation 0.01 has
about 0.995% probability of a weight above two. These are consequences of the
specified prior/readout and sampling, not a new fitted threshold or a guarantee
of zero false corrections. No parameter was adjusted to these cases.

The complete-consumer cost gate fails. Single-strand cases take about 0.5–1.3
seconds each and thousands of density evaluations in this diagnostic Python/
native oracle. Neither both-strand example completes within its 60-second budget:
test stranded slot 2780 makes 149 calls and deferred unstranded slot 3120 makes
164. The latter has only three reads. These are bounded cost probes, not a
back-to-back speed comparison; the two bounded jobs overlapped briefly. Do not
port the nested evaluator or expand to whole panels at this cost. The next
implementation question is economical evaluation of the same three integrals,
while retaining both RNA strands, local evidence, raw count uncertainty and
certified source precision. No statistical rejection or release accuracy claim
follows from a timeout. Instrument: `.cache/rigel_runs/2026-10-09_capture_prior/`
(`local_prior.py`, exact reference tests, `mutations.json`, `screen.json`, current
typed message snapshots, complete-consumer and tolerance receipts).

SHARED CAPTURE INTEGRATION (2026-10-09; owner-authorized cheaper calculation).
Replacing the scalar outer traversals with shared vector quadrature is rejected
for performance. The raw observations, typed local factors, nuisance reference,
prior and correction score stay fixed. Inner tolerance remains `1e-8`, outer
tolerance `1e-6`. A fresh pre-implementation baseline on the three-read both-strand
object times out at 60 seconds after 172 evidence calls; this is the cost failure
the experiment attempts to repair, not a statistical falsification.

The first global-vector version reports convergence but misses 16.6% of a narrow
analytic curve's slab mass. Two prewritten equivalence checks catch it. Keeping
the original per-interval relative checks repairs those failures; the only
absolute tolerance is the smallest normal floating-point value, allowing fully
underflowed pieces to terminate. This does not introduce a positive density floor.
Reported integration error is an estimate, not certification of unobserved peaks
or the inherited input tables.

Four serial back-to-back pairs agree within `7e-12` in log weight, but all take
more time and more expensive likelihood evaluations. The old scalar routine
already memoizes common nodes; sharing the vector does not remove its inner
RNA/strand integrations, and the new backend adds subdivision work.

| Frozen object | Scalar / shared wall seconds | Scalar / shared evidence calls |
|---|---:|---:|
| test zero-DNA boundary 513 | 1.215 / 3.169 | 9,453 / 23,617 |
| test stranded captured exon 1000 | 0.526 / 1.290 | 4,774 / 11,775 |
| test stranded captured exon 1036 | 0.561 / 1.336 | 4,940 / 11,711 |
| test deferred unstranded exon 2958 | 1.002 / 2.365 | 4,955 / 11,726 |

The stranded both-strand slot 2780 times out at 60 seconds in both arms, each
after 148 evidence calls. That already-launched pair finishes under its bound;
the remaining unstranded pair is cancelled once the completed pairs reject the
implementation. The earlier three-read failed-first receipt remains separate.
These are Python/native-oracle timings on frozen inputs, not whole-tool speeds.
Input and native hashes match across every executed pair. The unstranded control
still uses a frozen stranded background and is not a background-estimation test.

Retain the reference and mathematical checks, archive this performance-negative
implementation, and seek the owner's requested external review before another
numerical framework or native port. No statistical model has been rejected by
this one implementation's cost, and no capture-weight, transcript/gene or
production accuracy improvement is claimed. The final prototype passes 28 checks,
catches seven actual mutated implementations and passes the 11-test docs boundary
gate. Hash checks preserve 104 production source files and the installed native;
the standing 2,814-test suite is not rerun for this source-preserving experiment.
Instrument:
`.cache/rigel_runs/2026-10-09_shared_capture/` (prewritten tests, initial failure
receipt, scalar/vector receipts, mutation copies and source-identity verification).

The shipped finite-support specification remains failing outside the production tests;
the range-only prototype passes those original cases but retains the phase failures.
The local-unit extension fixes those phase failures and retains the spacing limitations
measured above.
Do not claim complete locality certification from the repaired witness and zero-coordinate
gates. The next representation test must preserve physical support and demonstrate
numerical convergence, not introduce a positive biological rate floor or a capture threshold.
Do not re-run whole-genome validation for the no-op witness input; price a count-changing
coordinate repair first on these small controls and cached conditions.
Instrument: `.cache/rigel_runs/2026-10-08_message_locality/` (`audit.json`,
`actual_faces.json`, frozen contexts, `local_witness.py`, before/candidate JUnit receipts,
`witness_mutations.json`, paired panel censuses and `coordinate_unrepaired.xml`).

EVIDENCE-CURVE ADMISSION (2026-10-08). The authorized cutoff-removal contrast is now
measured with the separate-channel density curves. Structural exclusions, composition
requirement, zero-count anchors, existing weights, domain and earlier refits remain fixed.
Previously admitted curves and weights are bit-identical. The final-refit training set
adds 555 objects on test zero (12 introns, 543 exons; total reliability weight 114.4624)
and 2,590 on RNA-long unstranded ON (267 introns, 2,323 exons; weight 295.2835).
The corresponding old populations contain 699 and 3,989 objects with weights 699 and
3,393.3259. All admitted new rows have positive weights and nonconstant curves.

The old admission mask fails four of eight new gates; the corrected selector passes
all eight and six deliberate defects exercise each gate. Ten earlier population-reference
gates remain green. Curves require 85.5 seconds for test zero and 478.5 seconds for RNA-long
with four workers; pilot curves reproduce exactly. Both fits meet the unchanged 1e-6
objective-gap requirement (3.54e-7 and 5.24e-7). These are frozen final-refit comparisons,
not self-consistent calibrations or release candidates. Capture-reference metadata is held
fixed so this test does not silently alter the existing reader's reference selection.

DNA incidence-count errors (%) on RNA-long unstranded ON:

| class | variance cut retained | cut removed | cut removed, control smoothing strength |
|---|---:|---:|---:|
| Regions | 16.7259 | 19.6015 | 19.5770 |
| Boundaries | 18.8416 | 20.6413 | 20.5961 |
| Introns | 30.9444 | 25.6726 | 25.6625 |
| Exons | 16.9364 | 20.1756 | 20.1498 |

The intron improvement comes at the cost of exons and boundaries. Do not promote this
model or tune the cutoff to hide that tradeoff. Under the changed unsmoothed population,
the old introns gain 20.4336 weighted log-likelihood units while old exons lose 5.6796;
new introns gain 1.7990, whereas the many new exons gain only 0.03494. This identifies
which inputs support the movement; it does not prove whether weighting, local-factor
misspecification or population uncertainty is the underlying limit.

The test zero control remains stable. Its inferred region/boundary incidences are
21.12983/19.58098 with the cut and 21.12978/18.78632 without it. The 21 region incidences
fixed as DNA by annotation remain a known model mismatch. A crossed attribution shows
why the boundary change should not be credited to new evidence:

| population shape / pseudo-region strength | region incidences | boundary incidences |
|---|---:|---:|
| Old / old | 21.12983 | 19.58098 |
| New / new | 21.12978 | 18.78632 |
| New / old | 21.12983 | 19.58009 |
| Old / new | 21.12979 | 18.78731 |

The existing one-uniform-pseudo-region mix is diluted when total training weight grows;
this accounts for almost all the zero-control change. Neutrality to flat evidence in the
likelihood objective alone does not guarantee neutrality after that postprocessing
(`EQUATIONS.md`, **Flat evidence and post-fit population mixing**). No smoothing rule is
changed here. Raw likelihoods still inherit the earlier face-model limitations, and the
count readout still has the duplicated density prior.

The subsequent owner-approved equal-weight contrast is recorded below. No class multipliers,
intermediate weight powers or new capture prior are selected. The full count candidate
remains frozen; this deferred-stratum failure is not a new 0.8.0 release veto.
Instrument: `.cache/rigel_runs/2026-10-08_evidence_admission/` (`falsification.log`,
`tests.log`, `mutations.json`, frozen inputs, curve/pilot identity receipts, optimizer
certificates, matched count replays, crossed smoothing controls and influence reports).

EQUAL OBJECT WEIGHTS (2026-10-08). With explicit owner approval, the completed
wide-admission curves were fitted with one unit per object instead of the retained
posterior-derived weights. Objects, curves, grid, objective, 1e-6 numerical objective-gap
requirement, earlier refits, reference metadata and count solver remain fixed. The
post-fit uniform mixture uses the control's weight sum in both arms: 813.4624 on test zero
and 3,688.6094 on RNA-long. It is not diluted by the candidate's 1,254 and 6,579 objects.
All three refits receive bit-identical observations, opportunities and incoming beliefs
between the paired replays. Only the final fitted population differs.

The retained-weight implementation first fails two of four specification gates. Equal
weighting passes all four. Three deliberate source defects cover every gate: restore old
weights, change smoothing strength, or overwrite the frozen control weights. The two-density
optima are checked against closed-form likelihood derivatives. Fits take 1.2–1.7 seconds
each, pass the objective-gap certificate, and the four cached count replays total about
17 seconds. There are no new scans, evidence evaluations or EM runs.

RNA-long g50 unstranded capture-ON, DNA incidence-count errors (%):

| class | fresh retained-weight control | equal object weights |
|---|---:|---:|
| Regions | 18.8154 | 22.0372 |
| Boundaries | 19.7058 | 24.7044 |
| Introns | 23.1580 | 17.1435 |
| Exons | 19.4123 | 23.0449 |

The intron improvement again costs exon and boundary accuracy. On test g00 stranded
capture-ON, region incidences are 21.1297827 versus 21.1297759 and boundary incidences
18.7863310 versus 18.7829682. This negligible movement neither creates a new zero-control
failure nor resolves false capture. Twenty-one region incidences remain structurally
fixed as DNA despite their simulated RNA origin. Equal weights are not promoted to
larger panels or the production population; no intermediate weights are tuned.

A numerical robustness concern also emerges in the retained-weight control. A fresh
SLSQP fit from a uniform start and the earlier multiplicatively initialized fit differ
in objective by only 2.88e-9 per unit weight; both satisfy the same global-gap certificate
(4.96e-7 and 5.24e-7). Their population L1 distance is 0.08278. Replaying those populations
moves intron error from 25.6726 to 23.1580, exon error from 20.1756 to 19.4123 and boundary
error from 20.6413 to 19.7058. The equal-weight tradeoff has the same direction against
either control. A small objective gap does not establish a stable population or stable
downstream median (`EQUATIONS.md`, **Population objective and readout stability**).
This is a limitation of the prototype, not evidence to select the initializer with the
better truth score. Numerical/readout convergence and the inherited local evidence must
be resolved before this population can serve as the prior-once partner. No new prior,
smoothing rule or count readout is selected by this contrast.
Instrument: `.cache/rigel_runs/2026-10-08_equal_weights/` (`before.xml`, `after.xml`,
`mutations.json`, frozen-input hashes, `fits.json`, matched replay arrays/receipts and
`optimizer_sensitivity.json`). Production source, tests, goldens and installed native
remain unchanged.

POPULATION FIT STOPPING (2026-10-08 follow-up). The retained-weight discrepancy above
is largely explained by the loose numerical stopping requirement. More than 99.9% of
the population L1 difference is a transfer between adjacent grid nodes at DNA densities
6.8261e-6 and 7.4138e-6. Holding all other masses fixed, an independent one-dimensional
likelihood derivative locates the conditional mass transfer; its objective gains are
only 6.39e-9 and 9.16e-9 from the two starting fits.

Refining the same constrained likelihood from both starts, with no changed observations,
weights, grid, smoothing or earlier refit, reaches objective-gap bounds 8.05e-11 and
8.18e-11. Their population L1 difference falls from 0.08278 to 3.54e-5. Fits take 0.58
and 11.61 seconds. Matched final count replays give:

| RNA-long unstranded ON class | refined archived start | refined fresh start |
|---|---:|---:|
| Introns, DNA-count error (%) | 25.005960 | 25.005627 |
| Exons, DNA-count error (%) | 20.068190 | 20.068091 |
| Boundaries, DNA-count error (%) | 20.452101 | 20.452057 |

The maximum individual count difference is 1.2644 fragments, with total absolute
difference 11.2780 incidences; this is close class-score agreement, not bit identity.
The observation does not establish distinct exact optima or justify a new biological
prior. It establishes that the original 1e-6 objective-gap criterion alone was insufficient
for this readout. A production numerical tolerance is not selected from this one condition.

The equal-weight fit was also challenged with the tighter target. SLSQP stops on a
line-search condition at a 1.03e-7 gap, satisfying the earlier 1e-6 criterion but **not**
the requested 1e-10 target. Its population L1 change is 0.000506; a diagnostic count
replay gives intron/exon/boundary errors 17.142746/23.044106/24.704295%. This preserves
the equal-weight tradeoff against either refined control. It is not a certified tighter
equal-weight replacement. Do not add an optimizer or change the population model merely
to rescue that failed arm. All three replays verify their incoming refit arrays exactly.
Instrument: `.cache/rigel_runs/2026-10-08_population_stability/` (`refine.py`,
`refinement.json`, saved populations, matched replay arrays and `comparison.json`).

### the-edge-density-floor-under-capture
`priority: now — source-model checkpoint before the density consumer · kind: design · 2026-10-08`

The boundary's structurally pure-DNA count identifies its own capture-weighted density.
It does not guarantee a lower bound on the adjacent exon's average DNA density.
**Rule 5 is a level, and it is one-sided** in DESIGN assumes that ordering. The owner's
subsequent clarification, **The edge floor is an approximation**, authorizes an isolated
withdrawal experiment; it does not select a production replacement. **The paradigm** also refuses
absolute DNA-level transport across capture locales, so relabelling the existing EDGE
source as a level and forwarding it is not a coordinate-only repair.

Exact placement enumeration supplies a counterexample with a probe wholly inside an
exon, no unannotated RNA, equal RNA/DNA fragment lengths and no sampling noise. The exon
is `[0,2000)`, the probe `[0,100)`, fragment length 200, and each fragment has weight
`1 + 10 * overlap`, the simulator's current law. The 199 starts crossing the left edge
have mean weight **752.2563**; the 1,801 starts contained in the exon have mean weight
**29.0400**. The edge's DNA density is **25.9042 times** the exon's average. A common
normalization of DNA yields cancels from that ratio. No probe trails outside the gene.

This is not a dependence on linear binding. For overlap thresholds `t=1,...,100`, the
fractions meeting the threshold are `(200-t)/199` at the boundary and `(101-t)/1801`
inside the exon. The first exceeds the second at every threshold. Thus every
nondecreasing binding function with a positive increment over overlaps 0 through 100
has the same density ordering in this geometry. The comparison holds on reflecting the
whole footprint to the other exon end. It is a counterexample to a guarantee, not an
estimate of the defect's prevalence or its contribution to a release error.

An 81-case enumeration varies DNA and RNA lengths independently over 75, 200 and 500,
probe position over the two ends and centre, and binding over 0, 0.1 and 10. Eighteen
cases violate the assumed ordering. A separate verifier checks **465,156** fragment
weights through `CaptureSampler.fragment_weight` and actual genomic-to-transcript
probe projection, including mirrored footprints; its expected counts agree exactly.
Probe coordinates are diagnostic inputs only, never inputs to calibration.

An independently checked conditional-binomial alternative correctly accounts for the
recipient total's sampling noise under the same ordering. It matches an independent
Poisson optimization to **2.05e-12** over 108 comparisons, but retains the false premise.
With equal expected DNA and RNA counts in the example, both factors prefer an all-DNA
recipient over its true DNA fraction of one half: **97,984.48** log units for the old
factor and **90,786.42** for the conditional alternative, at the stated unit source
rates. These are factors in a noiseless construction, not complete posterior or panel
results. More sequencing strengthens this wrong factor rather than repairing it.

Instrument: `.cache/rigel_runs/2026-10-08_edge_contract/`, `calculate.py` and
`verify_geometry.py`, with JSON and log receipts. The calculation deliberately retains
exit status 1 for the failed generality condition; successful independent verification
does not turn it into a passing repair. No candidate, calibration solve or benchmark
was run. The old EDGE constructor and all production sources are unchanged.

The next proposed contrast is an isolated withdrawal of this cross-object floor while
retaining the boundary's own DNA evidence and all other messages. It is an ablation to
price information loss and propagation effects, not a selected production replacement.
The prototype is authorized; its accuracy benefit and information loss remain unmeasured.
No upper bound, detector, probe input, new coupling
constant or absolute-density propagation is selected by this finding.

ISOLATED WITHDRAWAL (2026-10-08). With owner authorization, the native worktree omits
only the published EDGE row. The EDGE face role remains, blocking accidental activation
of a gDNA level fallback; the boundary's own Poisson density measurement, every other
source, face map, lane licence and propagation pass remain. The pre-change binary fails
**18/24** specifications and the ablation passes **24/24**. Five actual native mutations
cover all 24: restore the floor, leave an onward copy, erase the source's own measurement,
enable a fallback lane, and remove unrelated claims. Neither this nor the following
census changes the density reader, population admission, weighting, prior multiplicity
or capture reference.

Fresh cached-input pairs were scored separately at pass zero and after the default
refits. Cells below are **region / boundary gDNA count error (%)**: summed absolute
per-object error against `slot_truth.npz`, divided by observed counts on that axis.
These are calibration census measurements, not transcript/gene error; no new end-to-end
benchmark was run. The control is the frozen count candidate with the already-certified
numerical message corrections, not 0.7.1 or an integrated detector-free reader.

| Condition | Pass zero: control | Pass zero: EDGE absent | Full: control | Full: EDGE absent |
|---|---:|---:|---:|---:|
| RNA-long g50 ss0.50 ON (deferred) | 33.887 / 37.894 | 34.722 / 37.945 | 41.676 / 32.913 | 61.662 / 45.946 |
| Ladder g98 ss0.50 OFF | 19.788 / 38.524 | 22.013 / 40.970 | 2.465 / 9.124 | 2.248 / 8.829 |
| RNA-short g50 ss0.99 OFF | 1.504 / 4.281 | 1.511 / 4.293 | 0.925 / 2.948 | 0.926 / 2.959 |
| RNA-long g50 ss0.99 ON | 1.762 / 3.904 | 1.763 / 3.904 | 1.466 / 3.348 | 1.466 / 3.348 |
| Test g00 ss0.99 ON | 0.470 / 0.433 | 0.470 / 0.433 | 0.004 / 0.024 | 0.004 / 0.024 |
| Test g50 ss0.99 ON | 1.682 / 1.246 | 1.654 / 1.247 | 1.337 / 1.134 | 1.313 / 1.135 |
| Test g50 ss0.50 ON (deferred) | 21.800 / 14.932 | 21.795 / 30.578 | 3.688 / 2.543 | 43.467 / 61.946 |

The zero control is exactly unchanged on both axes, not merely equal at the displayed
precision. In the test unstranded captured case, the first-pass region error barely
moves, while boundary and intron errors worsen. The later population refits accompany
a much larger region loss. This establishes sensitivity to the lost neighbour evidence;
it does not by itself isolate admission, weights or prior feedback as the cause.

The 1,000-fragment weak-strand OFF fixture from **COUNT-CANDIDATE READINESS** improves
substantially. Its observed counts were matched exactly to the previously validated
origin partition. Intronic error changes **64.57 → 46.25%** at pass zero and
**61.38 → 6.64%** after refits; exon error changes **52.18 → 23.67%** and
**45.55 → 7.23%**, respectively. Boundary error changes **14.52 → 11.87%** and
**16.28 → 6.84%**. With intron constraints already absent in both arms, this identifies
an important contribution from EDGE and its downstream effects in that fixture. It is
not evidence for restoring the intron constraint or for a strand/capture gate.

Conclusion: withdrawing the approximation trades substantial useful imputation for
less bias in some cases. It is not selected as a general replacement and does not solve
the zero-gDNA capture problem. The deferred-stratum loss is a robustness result, not
a new release veto. The next release work must retain this tradeoff without tuning a
cutoff to these panels. No full-panel or real-library run is justified for the blanket
ablation now. Instrument: `.cache/rigel_runs/2026-10-08_edge_ablation/`, frozen control,
native package, assertion-only mutation failures, paired census and weak-strand receipts.

### the-scorer-reads-a-census-length-law
`priority: now — released from parking (owner, 2026-10-02): the largest lever measured on stranded × capture ON; the open design is the capture × length half of the fragment-length review · kind: defect · 2026-09-24`
The E-step scores an unspliced fragment's length with gDNA's library census (`gdna_realized_pmf`) against RNA's spliced
census (`rna_pmf`), and the transcripts' effective lengths and the capture ruler read the same spliced census. Under
capture each census is its origin's law lengthened at that origin's own placements (the ladder's no-EM census: gDNA
216.7 → 240.9 bp, RNA 212 → 228.7), so on equal chemistry the two tables carry a false difference whose sign depends on
where each origin sits (ladder gDNA 12–18 bp longer, the test chromosome 11 bp shorter), and the ruler applies capture a
second time to an already-captured RNA law. Matched within exons the difference vanishes (219.00 against 218.98 bp).
Nil off capture. The scorer reading the uncaptured gDNA law alone tips the other way
(`ISSUES: the-scorer-reads-the-uniform-gdna-law`, refused).
CALIBRATION HAS A RELATED FRAME MISMATCH (verified 2026-10-06): `pipeline.py` supplies `gdna_pmf`
(uniform frame) and `rna_pmf` (junction-de-tilted but capture-selected) to `calibrate`.
`region_geometry.py` integrates them separately for contained and crossing opportunities and uses
the RNA law for junction rates. The component-opportunity count repair is algebraically correct for
its supplied inputs; it does not establish a common pre-capture frame. Overriding calibration's RNA
law with gDNA's on equal-chemistry simulations is a diagnostic, not an RNA-law estimator. A true-gap
control requires the separate simulator laws. The repair's observed pool shifts and unresolved
causal attribution are recorded in `ISSUES: calibration-detects-capture-on-a-capture-off-library`.

**Calibration-only law-frame diagnostic, 2026-10-06.** Eight fresh cached calibrations on the
uncommitted count-repair build: native input laws versus `rna_fl_pmf := gdna_fl_pmf`, only at the
calibration entry point, on four equal-chemistry conditions. No scan or EM was run and no scorer
law changed. The paired runs share payloads, drain choices, configuration and native build; the
input arrays stayed unchanged. All object totals agree exactly with `slot_truth.npz`, and the
matched-law arm has exactly equal gDNA/RNA geometric opportunities. The cache check initially
rejected the default `auto` strand tag; both arms then used the cached scans' explicit `XS`, with
all tally settings validated. Calibration threads remained at their default. Total wall time was
16.93 s, individual calibrations 0.17–1.05 s. This is a diagnostic intervention, not an RNA estimator.

The object metric below is `100 Σ|estimated gDNA − oracle gDNA| / Σ unspliced count`, separately
for regions and boundaries. Pool under-call is true minus estimated conserved gDNA fragments;
it is not an incidence sum. These are **calibration** measurements, not transcript/gene errors.

| Condition | Native laws: regions / boundaries (%) | Matched laws: regions / boundaries (%) | gDNA pool under-call: native → matched (fragments) |
|---|---:|---:|---:|
| Test chromosome g98 ss0.70 OFF | 0.78548 / 2.67339 | 0.78508 / 2.67228 | 1,197 → 1,206 |
| Test chromosome g98 ss0.70 ON | 2.67666 / 4.83476 | 2.67927 / 4.80878 | 24,651 → 24,707 |
| Ladder g50 ss0.50 OFF | 1.11596 / 2.98516 | 1.11800 / 2.99574 | 8,978 → 9,303 |
| Ladder g50 ss0.50 ON | 7.73254 / 12.99809 | 8.37262 / 13.27189 | 682,837 → 723,625 |

The net captured calibration loss is **not rescued** by matching the laws: the ladder's gDNA
estimate falls a further 40,788 fragments, with region and boundary errors both worse. On the
test chromosome, the two implicated exons remain at 336 / 316 against truth 3,137 / 2,904
(native-law estimates 341 / 319). Their original native-law counts reproduce the archived repair
exactly. The input-frame mismatch is real, but the proposed substitution is not a fix. It also
changes junction rates, RNA level lanes and the refitted landscape, so this result neither isolates
the new maps' contribution nor proves that a consistent length treatment is unnecessary. Do not
promote this arm to transcript benchmarks. Next separate frozen-state effects from population
refitting, beside the independent discrepancy-center contrast.

Receipt: `.cache/rigel_runs/2026-10-06_law_frame/` in the repository workspace (`run.py`,
`metadata.json`, `results.jsonl`, per-condition arrays and logs; mirrored from
`/private/tmp/rigel-law-frame-20261006/`). Native SHA-256:
`a0d204eaff734b5e3cde63ff8cf6538c2f9eecdff7f5e25622cf3b20178870c3`.
THE PRIZE (phase 0, 2026-09-30: oracle arms outside the tree, pinned and fractional, calibration untouched). Both tables
and RNA's normalisers in one uncaptured frame — the uncoupled uniform gDNA law for both, right here only because the
panel's chemistry is equal — on the ladder's stranded × capture-ON rows: transcripts −4.5 %, genes −14.5 %, gDNA pool
(EM + intergenic) −14.2 %, annotated −54.3 %, nascent −18.3 %; the simulator's drawn laws read the same (−5.0 / −14.5 %).
By row, transcripts / genes: g05 −3.0 / −9.1 %, g50 −5.8 / −16.0 %, g98 −12.5 / −20.5 %. The per-row gate fails 2 of 3
rows (g05 gDNA 79 → 627 fragments against a drain floor of 16; g98 nascent +7.6 %). The test chromosome's
equal-chemistry panels agree (stranded × ON genes −5 to −20 %). Off capture the ladder's pools lose 2–6 %, with the
exact laws too, and at g00 the uniform law is RNA's own unspliced census (244 bp; transcripts +0.7 %). Equal tables
alone give a third of the transcript gain and 60 % of the gene gain; feeding the ruler the uncaptured law carries most
of the rest (an isolated arm). The plain effective length is not isolable this way: an implicit splice gives each RNA
candidate its own length, so a shared table does not cancel there (3.2 % of buffered fragments on a gap row).
MEASURED ON THE SUITE GAP ARMS (re-recorded 2026-10-03 on the re-simulated panels, pinned and fractional).
`flgap_rna_short` capture ON reads transcripts 19.8 % stranded and 18.8 % unstranded, against 6.2 % and 8.2 % on
`flgap_rna_long`. `ruler_vs_truth.py`, annotated mRNA, median log
error of the shipped ruler and of the ruler fed the oracle gDNA: the short arm's partially probed transcripts (≤ ½ of
their bases probed) read +0.43 / +0.51 shipped and +0.41 / +0.39 at the oracle (ss 0.50 / 0.99), its unprobed ones +2.8 /
+3.0 (oracle +1.5); the long arm's read −0.24 / −0.15 (oracle −0.19 / −0.33) and +0.55 / +0.80; the fully probed class
reads 0.000 on both arms. The sign follows the gap: RNA's capture priced at gDNA's efficiency, below.
WHY IT IS PARKED. On real length gaps under capture the same frame loses even with the simulator's laws: the test
chromosome's stranded × ON genes +92 % (RNA 250 ± 150 against gDNA 100 ± 50) and +10.5 % (the reverse), while capture
OFF holds. The ruler prices RNA's capture at gDNA's per-object efficiency, so it sets `C_R / C_G = 1` where the true ratio
is `E_{Q_Gj}[r] / E_{P_Gj}[r]` (`r = f_R / f_G`; `Q` the captured and `P` the reference gDNA law at object j), exactly 1
only on equal chemistry; today's census tables partly mask it, a cancelling pair (the census frame against the ruler's
transfer). The oracle arm is EB-mixed and keeps calibration's gDNA denominator, so the transfer is the leading
candidate, not an isolated cause.
A COMPLETE DESIGN needs, in observed laws only:
- the relative law `log r(l)` measured where capture cancels: within objects, from strand on one-RNA-strand objects
  (`logit q_jb = α_j − log r_b`), stranded libraries only;
- RNA's own law on its own support, from spliced evidence with its selection specified (junction crossing — the
  de-tilt targets the realized `f·T` by design, `sj_opportunity.py` —, capture exposure, and splice geometry: +12.6 bp
  against +20 bp for contained exon fragments on ladder g50), never built from gDNA's law
  (`ISSUES: a-length-table-built-from-the-other-origins-law`), and not read from gDNA where there is none;
- the missing measurement: gDNA's length-resolved mass on RNA's placement objects, keeping their placement
  restrictions; scalar efficiencies, pooled banks and efficiency strata cannot reconstruct it;
- one mechanism supplying each component's table and its integral together, replacing the realized law and the coupling
  (about 300 lines in `fl.py`), measured together with what calibration reads.
THE NEXT STEP when taken up, fixed before any run (an external review, 2026-09-30): arms B_match (RNA lengths from the
simulator's capture-weighted yield at O's scorer laws, the gDNA locus denominator replaced on the same raw scale after
`assemble_priors`) and B_transfer1 (the same, RNA priced at gDNA's efficiency), against F, the floor, O and a true
no-length control, on both test-chromosome gap panels, the suite gap panels, the ladder and every test-chromosome
capture-ON condition, under the per-row gate. If B_match fails, the route stays parked. Land
`ISSUES: the-pooled-q-in-the-gdna-count` with the repair.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-30_frame_phase0/`: `driver.py` (the arms), `summarize.py` (the
read-out, gDNA as EM + intergenic), `runs/`, and `docs/` (the problem brief, both external reviews, the phase-0 plan and
report with its corrections). The census measurements: `~/Downloads/rigel_runs/prototypes/2026-09-30_fl_fix/`.

### the-gdna-landscape-collapses-at-low-depth
`priority: next — Tier 1, after the strand-overdispersion landing and A/B'd apart from it · kind: defect · 2026-09-28`
Calibration's gDNA estimate falls with sequencing depth. LBX0588 (stranded, capture ON) reads 0.11 / 0.45 / 0.83 gDNA
per deposited fragment at 1 % / 10 % / full depth. The EM follows it (gDNA fraction 0.16 / 0.62 / 0.91), so at 1 %
the RNA pool is about 9× too large. The VCaP mix, which has truth by read name, reads P/O 0.125 / 0.52 / 0.92.

The mechanism is the refits' gDNA landscape.
- A slot trains as located only when Var(log f_g) ≤ 1 nat², which in practice means it holds a gDNA fragment.
- At 1 % depth an object at VCaP's bulk density needs about 1.26 Mb to expect one fragment. So 0.23 % of the ~323k
  kernels locate, and those few are selected above the bulk.
- The rest are zero-count anchors, flat down to the grid floor. The refit loop converges to one mode at the floor
  (log ρ ≈ −17) in all three real libraries.
- Under that prior, sparse unspliced fragments go to RNA.

How much each part owns:
- The refits own 80 % / 69 % of LBX0588's drift (ln P_full/P_d, from 1 % / 10 %) and 86 % / 63 % of VCaP's.
- Pass 0 drifts less. It nets an under-call on no-RNA objects against an over-call on no-gDNA objects.
- The capture reference owns none of it, because it is computed after the solve, but it goes missing: its members
  fall to 0 at 1 % for VCaP, LBX0588 and MO_3021 (0 at every depth for MO_3021), which sets every efficiency to 1. On
  the test chromosome at capture ON the members are g001 0 / 0 / 57, g01 0 / 73 / 329, g05 20 / 313 / 352 (d100 /
  d10 / full), and a reference that does form at low depth reads above proportional (×1.97 at g05 d100, ×1.48 at g01
  d10). What that costs the EM is unmeasured, and `calibration_vs_oracle.py` cannot see it. The reference is a yes/no
  switch; its continuous form, an efficiency whose uncertainty grows as members fall, is a direction to price.
- MO_3021 barely drifts (0.82 of full at 1 %): 65–80 % of its gDNA sits on intergenic regions and gene edges, which
  read their count at any depth.
- The level survives in aggregate. VCaP's intergenic regions at 1 % hold 0.0097 of full depth's true gDNA.
- The drift is not the od's: calibration's gDNA share falls with depth under every od arm. LBX0588 per deposited
  fragment at 1 % → 10 % → full, beside shipped's above: binomial 0.113 → 0.471 → 0.835, the 0.2 value 0.098 → 0.145 →
  0.679. The VCaP transcriptome half drifts too (gDNA fraction 0.50–0.63 % at full depth against 0.34–0.41 % at 10 %
  and 1 %, every arm); MO_3021 (0.075 → 0.077 → 0.091) and the VCaP DNA half (0.045 → 0.044 → 0.054) barely move. The
  0.2 value deepens LBX0588's collapse at 10 % (0.145 against 0.450), so the od landing and this repair are A/B'd
  apart, in that order.
- Strand inputs, an untested cause: two inputs fitted before pass 0 move with depth. MO_3021's gDNA od reads
  0.171 ± 0.042 at 1 % against 0.0084 at full depth; LBX0588's strand specificity rests on about 150 spliced reads
  (0.913 / 0.919 / 0.934 at 1 % against 0.937), and its pass 0 follows the same order (0.5465 / 0.5555 / 0.5641).
  MO_3021's pass 0 rises as depth falls (0.199 / 0.164 / 0.139), the opposite of the others. On the test chromosome,
  planting the true strand model moves the error by at most 6.3 points. Re-check now that od = 0 has landed and once
  κ's fix lands (`ISSUES: strand-overdispersion-one-shared-value`).
- The EM does not always follow calibration. On VCaP it recovers much of what calibration loses (EM recall 0.614 /
  0.747 / 0.853 against calibration P/O 0.125 / 0.518 / 0.922); LBX0588's EM follows calibration down (0.164 / 0.620 /
  0.906). Why is untested: the libraries differ on three axes at once — strand specificity 0.9998 against 0.91–0.94,
  spliced share 51 % against 1.3 %, fragment length (VCaP DNA 307 bp and RNA 244 bp against LBX0588's 104 bp).

On the test chromosome's stranded depth family:
- Under a depth-invariant true landscape, the landscape owns 69–94 % of the capture-ON low-gDNA drift at d100 and
  d10.
- The capture-OFF drift is per-object.
- Injecting the full-depth landscape, shifted by log f, cures g01 and g05 at d100 capture ON (−81.1 → −0.8 and
  −52.6 → −8.3 points of O). It only halves g001, does nothing at capture OFF, and worsens g25 and g50 at d100.
- True values on every structural region at depth d still fail at g01 and g001 d100. So fitting count/E kernels on
  sparse counts is itself part of the defect.
- The same fragile regime at full depth: at g001 capture ON, removing the one 19-fragment `test_blank` slot from the
  training set moves the annotated chromosome's net gDNA error +3 → −630 (Σ|Δ| 983.1 → 1,090.9, n_train 754 → 722).
  The depleted mode rests on very few located kernels at low gDNA.
- Under a correct, depth-invariant landscape the capture-ON residual is misplacement, not loss: g01 d100 no-gDNA
  +62.95, mixed −45.91, no-RNA −23.12 (net −6.08); g001 d100 +83.34 / −55.79 / −24.81 (net +2.75); g001 full +20.63,
  cancelled by the landscape term −16.23. A repair is scored BY CLASS (no-gDNA, mixed, no-RNA), never on the net.
- The collapsed landscape beats the true one at g25 d100 ON (shipped −6.8 against −10.9), and the shifted full-depth
  landscape moves g25 −6.8 → −10.6 and g50 −4.2 → −8.0: two compensating errors, and a repair must not regress these
  rows.
- Not separated: whether the landscape or the density fit removes pass 0's mirror floor (RNA-only objects 1,894.5
  against a predicted 1,806.9, and 65.3 after the refits). No arm isolates the message layer: shipped − silent is
  +33.65 (g05 d100 ON), +79.39 (g01 d10 ON), +34.88 (g25 d100 ON).
- The depth family is ss 0.99 only, so the drift on in-scope unstranded × capture OFF is unmeasured; TESTING §0c still
  describes these configs as the ruler's substrate only. RULED (owner, 2026-09-28): add ss 0.50 rows at d100 for g01
  and g05 only.
- An oracle landscape trained on the library's own true counts leaks the answer (it gave the landscape 29–63 % of the
  capture-OFF drift, against about 0 under the depth-invariant true landscape): use the depth-invariant one.

First, alone: a failed refit overwrites the last landscape with None after the belief was solved under it, losing the
reference, the efficiencies and the debug record. It fires in sparse fits.

Direction: at its maximum, the exact Poisson mixture keeps the prior's implied total equal to the pooled count
(Σ E_i·E[ρ | c_i] = Σ c_i). The shipped fit loses that, and a fit at depth f should equal the full-depth fit shifted by
log f wherever the data resolve it.
- The location floor and the anchor E-step that brought the zero controls down
  (`ISSUES: gdna-landscape-trains-on-false-positives`) are what collapse here; a repair must not regress them.
- `estep_all` and `row_kernel` stay refused (`ISSUES: the-landscape-training-population-arms`).

Not this entry: VCaP's 7.8 % full-depth deficit, already present in pass 0
(`ISSUES: vcap-dna-only-objects-lose-gdna-at-full-depth`).

Still unmeasured:
- The real-data arm: the full-depth VCaP landscape, shifted by log f, on the 1 % and 10 % caches; check that it fires
  (`TRAPS: could-the-arm-have-fired`). The VCaP 10 % seeds put the mode at −14.70 / −14.56 / −17.05 with P/O 0.514 /
  0.526 / 0.516 against a bulk at −11.74, so P is insensitive to where the collapsed mode sits. Run it through
  `v2/scripts/cal_nodebug.py` under `one_at_a_time.sh` (4.4–4.8 GiB peak): one real-genome job, never `_debug`.
- What the drift costs the transcript table: `quant_accuracy.py` and `prior_vs_oracle.py` never ran on the depth
  family, nor `prior_vs_oracle.py` on `full_lowg`. The real-library EM numbers are sampled draws (seed 0), and the
  read-name scorer needs hard labels, against the fractional-benchmark rule. If full depth is roughly right, LBX0588's
  RNA pool is 8.9× too large at 1 % and 4.0× at 10 %; VCaP's DNA sent to annotated transcripts is 11.8 / 5.9 / 4.4 %
  of the DNA.
- Full-depth cfRNA has no truth: MO_3021's refits remove 34 % of pass 0's gDNA (0.139 → 0.091 per deposited fragment;
  exons 21,552 → 4,156; exon|exon 11,363 → 367). Only a truth-labelled mix can settle it (the truth mixes of
  `ISSUES: splicing-artifacts`); a √n depth extrapolation fails on VCaP (8.2 / 2.6 against a true 7.0).

Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/`: `pages/explanation_depth.md`; the depth-drift
chain `rerun/scripts/depth_harness.py`; the arms `shifted_landscape.py`; the tables `rerun/scripts/tables.py` (test
chromosome) and `v2/scripts/tab_v2.py` (real libraries); the substrates `test_reference/scenarios_depth_d100`, `_d10`
and full depth. RULED (owner, 2026-09-28): the harness stays outside the tree until a repair is A/B'd; the one tree
addition is the gDNA-only row of `ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool`.

### gdna-landscape-trains-on-false-positives
`priority: next — Tier 1, a constraint on the landscape repair: its zero controls must not regress · kind: question · 2026-09-02`
What is still open in the landscape estimator, each with its number: (a) under capture a short dark region's
mass is split over both modes of the previous fit (deferred 1.01–1.02×, stranded ON 1.004×); (b) the zero-RNA
controls move ±4 % (`g05 ss.99 ON` +4.4 %, `g98 ss.50 OFF` −1.4 %); (c) the Beta(½,½) reference still decides
a blind slot. CLOSED here 2026-09-14, (d): a slot whose solve is wider than one nat² in `log f_g` no longer
trains (`DESIGN.md` §7.1 rule 4) — the vertex-solved exons that trained 3,771 false fragments at the first
refit on `g00 ss.99 OFF` are out, and the four ladder zero controls read 282 / 194 / 265 / 172. Refused here:
excluding κ-dead exons (`g50 ss.50 ON` 2,691 → 56,422), AMBIG in the final fit (worse 25/32), and the two
other readings of the floor (`DESIGN.md` §7.1 rule 4). Its instrument, `landscape_training_census.py`, was retired 2026-09-14 (in git).

### the-background-dispersion-assumes-a-pure-intergenic-pool
`priority: next — Tier 1, its own A/B after the landscape fit · kind: defect · 2026-09-28`
`fit_intron_background` (`density_deconv.py`) fits the gDNA background's mean and its dispersion α by moments on the
intergenic regions. It assumes they hold only gDNA, and its own comment says unannotated transcription breaks that
(`TRAPS: purity-is-a-property-of-the-annotation`).

On the test chromosome's `full_lowg` (stranded, capture OFF):
- The unannotated transcription on `test_blank` raises the fitted mean 12.7×. It drops α from ∞ (no extra spread) to
  0.142 at g001 and 0.405 at g01.
- The lower α widens the prior on intron gDNA. Changed alone, it moves the annotated chromosome's gDNA by −197 / −687
  fragments; the shipped figures are −189 / −643.
- It owns 27 % / 38 % of the annotated chromosome's gDNA under-call (51 of 192 and 235 of 624 fragments).

Real intergenic regions hold strand-coherent RNA too, about 10–11 % on VCaP, so the route is live on real data. Its
size there is unmeasured.

Read it with the gDNA-only objects apart:
- On `full_lowg` capture OFF, `test_blank` itself carries 8,844 of the 9,194 Σ|Δ| at g001: the designed control,
  locked to gDNA by structure, not a defect. `calibration_vs_oracle.py` reads 2.31× = (9,825.1 + 20,009.7) / (1,170 +
  11,764), driven by it; without it the reading is 0.92× (0.923× shipped, 0.950× with a gDNA-only background). Every
  od arm reads the same 2.31×; the 0.2 value lifts the capture-ON P/O 1.095 → 1.307.
- The annotated remainder still carries the shadow, through the background's α, so a clean reading needs the gDNA-only
  background fit beside it.
- `calibration_vs_oracle.py` does not report structurally gDNA-only objects on their own row. RULED (owner,
  2026-09-28): add the row, the only tree addition for Tier 1. Until it lands, the debug loop must not pick
  `full_lowg` capture OFF as the worst in-scope case on this number.

A repair should read the half of the statistic the contaminant cannot reach, as its own A/B after the landscape's.
Whether it needs a new constant is for the derivation.
- fl's one-sided off-target rate has the same "intergenic = gDNA" contamination in a different estimator
  (`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value`).
- Under capture, intron regions next to probes are not off target either
  (`ISSUES: intron-seeds-near-probes-are-capture-enriched`).

Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/review_alt/bg_split.py`. It fits the background
on the gDNA fragments alone, beside the shipped fit. See also `pages/explanation_lowg.md`.

### capture-on-overcalls-gdna-at-low-gdna
`priority: next — Tier 1, read after the landscape repair · kind: defect · 2026-09-28`
On the test chromosome's `full_lowg` with the `test_blank` control excluded, capture ON is about twice as wrong as
capture OFF, and over-calls: annotated-chromosome Σ|Δ| 983.1 ON against 349.5 OFF at g001 and 2,846.3 against 1,506.7
at g01; the net is positive ON and negative OFF (`blank_arm.txt` +1,089 / −832; a second table +1,193 / −816). It
looks fine overall because capture removes the unprobed shadow. The mechanism is the strand resolution limit: capture
concentrates the low gDNA onto probed exons holding about 1,500 RNA fragments each (at g001, exon true gDNA 70 → 717
under capture), where the summed 0.42·√(n/I) is 1,135 against 1,107 true gDNA, and it errs both ways (g01 ON
single-strand exons +406). The located part is the AMBIG exons at g01 ON — 84 slots holding 524 true gDNA fragments,
read as 1,517 — at +993 shipped, +939 under the oracle landscape and background, +1,151 silent, +2,870 at pass 0. At
g001 the silent policy reads AMBIG +8,398 (capture OFF, P/O 15.58) and +5,029 (ON) against shipped −3 / +51: only
transfer's final solve removes it. The capture-OFF AMBIG over-call is +9 to +22. Untested candidate:
`ISSUES: the-atom-at-an-unwitnessed-both-strand-slot`. Whether any in-scope mechanism beats the resolution limit is
open.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/`: `results/blank_arm.txt`,
`review_alt/arm_classes.txt`, `results/floor_test.txt`.

### strand-plug-in-bias-on-sparse-libraries
`priority: next — Tier 1: the shallow-object loss below, not the plug-in alone, owns the capture-OFF depth drift; invisible on the ladder, no instrument in the tree · kind: defect · 2026-09-21`
Calibration's per-object gDNA fraction is a plug-in that conditions on the object's own fragments, this one
included; the exact statement is the leave-one-out posterior, and the plug-in understates gDNA by −21 / −22 /
−19 / −11 / −2 % at 1 / 2 / 5 / 10 / 50 fragments. The pool-level bias is (live objects) / (unspliced incidence):
0.2–0.4 % on the ladder, where it is guaranteed invisible, and MEASURED 2026-09-21 on real data 0.216 / 0.120 /
0.308 on the cfRNA libraries LBX0190 / LBX0588 / MO_3021 (72–91 % of their live objects hold five or fewer
fragments and carry 15–44 % of the incidence) and 0.055 on the deep VCaP library. Any per-fragment use of the
fraction must take the leave-one-out posterior and must not apply the fragment's own strand term a second time.
The plug-in is not the whole cause: shallow gDNA-only objects lose far more gDNA than it or the closed zero-RNA floor
predicts, even under the oracle landscape and background. On `full_lowg` g001 capture OFF 214 objects (O 304) lose
193 shipped, 155 with a gDNA-only background and 149 under the oracle landscape and background, against bounds of 107
(strand only) and 52–70 (the closed floor); at g01 OFF 556 objects (O 2,158) lose 502 / 338 / 301 against a bound of
429. The ladder's unstranded rows, where no strand plug-in exists, show the same staircase: shallow objects understate
gDNA by 15–40 % (real / predicted 1.04–1.40 at 1–2 fragments, 1.00 at depth). The depth family's [1,3) bin loses
49–99 % at d100 against 4–45 % at full depth; VCaP at 1 % reads its pass-0 introns at 169 against 1,915. The candidate
is the median read-out of a skewed shallow posterior; splitting it from the plug-in needs the kernel's posterior. This
loss owns the capture-OFF depth drift of `ISSUES: the-gdna-landscape-collapses-at-low-depth`.

### latent-defects
`priority: next — release-critical: the paths a user reaches first, each its own commit with a falsification test verified failing; the number-moving ones each in its own window · kind: defect · 2026-09-28`
Correctness defects from the 2026-09-27 review and the cleanup. A test's expectations may move; no panel number moves
unless marked.
MOVES A NUMBER (each alone):
- `index.py`'s `s_tx.nrna_n_contributors = …` should be `+=`: a single-exon transcript covering two merged spans keeps
  only the last count (the nascent table's `n_contributing_transcripts`).
- `Scenario.build_oracle`'s fixed-total branch never reads `gdna_fraction`, so six fixtures simulate NO gDNA, among
  them `test_scan_order_independence`, whose docstring says it puts gDNA in introns. Raise on the combination, then
  repair; the expectations move.
- `lift_choices` keys records on (ref, start, end, align_strand, sj_strand) while the native key adds the introns and
  the hypotheses, so tied records that differ in introns pool: oracle truth can replay the wrong record, or `drain`
  can raise. Frequency unmeasured; only oracle truth moves.
- The simulator's fragment-length sampler: `.astype(int)` truncates (mean frag_mean − 0.5 while `fl_pmf` is evaluated
  at integers), and an empty [frag_min, frag_max] window loops forever. Land with the next re-simulation.
- Antisense at p_sense = ½: the deterministic path books a mismatched fragment as sense, the EM as antisense. Goes with
  the scoring unification (`ISSUES: hygiene-ledger` (h)).
- The leading intron: a first intron starting at the fragment start makes `build_segments` emit no leading segment
  and `left_sj` mislabel every junction end. Fix the executable specification first; it re-caches.
LATENT (no current path reaches it):
- `run_parallel`: a throw from fn(0) leaves `task_` racing, and a worker throw terminates; the ψ and EM callers do not
  catch.
- nanobind: `nb::cast<A>(h).data()` in `em_solver.cpp` and `solve_kernel.cpp` converts into a temporary on a dtype or
  order mismatch and dangles (GIL released); outputs without `.noconvert()` write into a discarded copy; the scatter
  bindings downcast float64 silently; `StreamingScorer` has no `nb::keep_alive`.
- `--overhang-alpha 0` / `--mismatch-alpha 0` give 0 × −inf = NaN: require alpha in (0, 1] or carry the gate as a bool.
- The scan cache: `deposit_digest`'s fixture offers `hypotheses=()`, blind to gap arbitration and deferral;
  `_payload_from_parts` skips `from_dict`'s conservation checks; the two `.npz` files and the manifest are written in
  place, with no temp-and-rename. One validating `AccumulatorPayload` constructor would serve the scan, the cache and
  the drain.
- GTF parsing: a final attribute without `;` is dropped; a missing `transcript_id` escapes `warn-skip`; an exon on
  another chromosome or strand is merged silently.
- Simulator config: a capture entry without `enabled:` and without probes yields an empty `CaptureConfig`; gene and
  transcript keys are unchecked (a typo runs on the default of 100); an empty section fails with a TypeError;
  `sim_command` creates the output directory before validating.
- `splice_graph._ref_slices`: unsorted input drops exons silently.
- An arm that injects `rna_sense_frac` without `n_rna_obs` silently kills the strand channel.
- `fast_exp` on an all −inf row: `static_cast<int64_t>(NaN)` is undefined behaviour.
- `connected_components_native` returns a 2-tuple early and a 5-tuple normally.
- `FragmentAccumulator::finalize` drops `t_offsets_`'s leading 0.
- `_resolve_core`'s non-empty-exons precondition is stated only in its doc comment (all four callers guard it): an
  assert.
- `genomic_span` is 0 for a transcript at position 0 (a truthiness test).
- `set_sj` narrows to int32 unchecked, justified only by a census.
- The buffer's finalizer, run by garbage collection on a thread holding a logging lock, deadlocks joining the writer.
- The report's meta line puts the sample and index names into `innerHTML` unescaped. (Its capture note is in
  `ISSUES: the-format-changes-to-batch-before-release`.)
- An unknown quant YAML key raises without listing the known keys or hinting at an older build, so a 0.7.1
  `config.yaml` is refused unhelpfully; MANUAL's "rerun the exact analysis" needs a caveat.
- A splice blacklist present without a recorded store is dropped with no log line: a one-line warning
  (`ISSUES: an-unrecorded-splice-blacklist-is-dropped-silently`).
- `transfer_rows.h`'s `level_bound_row` `std::max(v, 1e-12)` never binds on a reachable input: delete it together
  with an assert v > 0 in the binding (today v = 0 reads 1e-12; without the floor it would be NaN).

### sj-strand-tag-chosen-from-the-first-reads
`priority: next — release-critical, user-facing: a library whose first spliced reads lack the tag loses its strand, and the resolved mode is recorded nowhere · kind: defect · 2026-09-28`
The scanner picks the junction-strand tag mode from the first 1,000 spliced reads (`detect_sj_strand_tag_native`) and
returns "none" if none carries XS or ts, so a library whose first reads lack the tag loses its strand. An unknown
`--sj-strand-tag` spec silently becomes XS_TS (`parse_sj_tag_spec`'s default), and the resolved mode is never logged
or written to `summary.json` or `config.yaml`, so a run cannot be reproduced from its outputs. The repair reads the tag
per record; it moves numbers only on libraries whose first 1,000 spliced reads lack it.
Instrument: a fixture with the tag absent from the first reads.

### hygiene-ledger
`priority: next — release-critical for its release-gating docs (g) and the flaky reorder gate (e); the rest Tier 4, batched between A/B windows, a numeric no-op landing between any two · kind: hygiene · 2026-08-31; the review of 2026-09-22`
What the reviews left, each its own commit; content-only changes keep the collected count unchanged.
(a) ROTTEN BUT LIVE, moving an instrument's numbers when repaired: `quant_accuracy`'s oracle arms are undrained
(documented there).
(b) DEAD, DEFERRED: `region_span_count` (the retired length channel's) and `region_end_count` (the retired
total-density landscape's), tallied per fragment by the accumulator and carried through the payload and the caches,
read by nothing. Deleting them changes the payload schema and re-caches both panels; step 1b re-cached every panel
without taking them, so they ride the leading-intron fix's re-cache (`ISSUES: latent-defects`; owner, 2026-09-28).
Until then `accumulator.h`'s "each population stores only the channels something READS" is false. The native
`FragmentAccumulator` also still reserves, fills and exports `sj_strand` and `merge_criteria` per fragment, which the
buffer drops at `from_raw`.
(c) PRODUCTION-DEAD, KEPT AS TEST OR INSTRUMENT SURFACE: the resolver's `intron_bp`, the tested half of its overlap
profile; `rigel.sim.benchmark` and `Scenario.build_oracle`, test tooling inside the package;
`region_init.has_own_composition_evidence`, which tests use as a check; `priors._project_regions_to_loci`, kept for
`prior_vs_oracle.py` (`assemble_priors` calls `_region_locus_shares` directly); `fl`'s `GdnaContrast` and
`GdnaRealized`, computed and discarded, settled 2026-09-30 as test surface and taken out of `__all__`;
the vendored `_cgranges_impl` extension, still compiled, optimised, LTO'd and shipped for a test-only `query()`, while
`index.py`'s docstrings present it as a quantification facility.
(d) CLAIMS NOT RE-DERIVED on the current tree, left standing: that most in-scope error sits at the simplex
vertices (`simplex_logodds`, a relay-era measurement); `sweep`'s refused deferral of UNIDENTIFIED slots to the
prior (priced 2026-07, with no refusal entry here); `region_geometry`'s "no per-region spliced floor" A/B
(relay-era); and `fl`'s crossing pools called "gDNA by structure" because mature RNA never crosses an
exon|intron boundary, while RNA that has not spliced there does.
(e) GATES THAT CHECK LESS THAN THEIR NAME:
- CLOSED (2026-10-09), scheduler-dependent reorder gate: another full-suite failure reproduced the
  previously recorded race dependence. Three 16-thread scans legitimately returned serial order;
  no calibration call occurs in that test. Its replacement deliberately reverses complete scanner
  chunks through the real scorer, locus builder and EM, requiring a nonidentity permutation of the
  same fragment set and identical output in both sampled and fractional modes. Existing thread-count
  comparisons stay. Removing the actual fragment-ID sort in a frozen package breaks both new cases:
  sampled counts move by two fragments and fractional counts by up to 5.91e-12. No scanner, locus,
  production model, tolerance or golden changes. Receipt:
  `.cache/rigel_runs/2026-10-09_rna_coordinate/` (`main_suite_race_failure.xml`,
  `test_scan_order_replacement.py`, `permutation.xml`, `permutation_mutation.json`).
- `test_sweep.py::test_gdna_sweep_zero_gdna_pin_and_monotone` asserts `f_g < ½` on slot 3, the AMBIG|intron−
  boundary, while the AMBIG region it means ends at 0.949 under the silent policy, and nothing in it checks
  monotonicity; the mature-exon chain's tests give the same `f_g` with the junction spliced or not, so none of them
  exercises the junction reads.
- `test_conserved_mass.py::test_the_mass_is_the_PER_BASE_attribution` allows `count` whole fragments where its
  docstring claims `count` half-ulps; `test_region_geometry.py` still builds float64 banks as uint64 fixed point in
  its fixtures.
- `density_model.count_observable_masks` has had no production consumer since od = 0 landed (2026-10-01): it selected
  the retired gDNA od fit's seeds. Delete it as a numeric no-op, rewriting `test_gdna_density.py`'s mask block and the
  selector in `test_accumulator_span_unbiased.py`.
- The report test's synthesized v3 summary omits `sj_blacklist_size` and `sj_blacklist_loaded`, so it renders
  "detection off", a shape no real run produces; no test asserts the splice-artifact note.
- The buffer's finalizer test asserts nothing about "still writes what it holds", and its writer-error test calls
  `cleanup()` before consuming, the reverse of production.
- The scan-worker exception gate runs one worker, so deleting `output_queue.abort()` stays green.
- The annotated-BAM ZS check is set membership only, so a label swap passes.
- The fragment-length proof's L check is still guarded by `if contained is not None`; two spec tests carry leftovers.
- The docs-citation gate never checks that it scanned anything (an empty glob or a `.hpp` suffix passes); a URL
  containing `/docs/` fails with no remedy; the TRAPS heading regex is written three times.
- The log-odds window's single source is ungated: restoring a default stays green, and about 23 literal 10.0s sit in
  the tests.
- The deleted float64 pass-through gate was never retargeted at the live boundary masses; `test_substrate.py` still
  says `sj_mass` arrives per strand "for artifact detection".
- `test_a_blacklist_built_from_a_store_loads_into_the_resolver` asserts a size that Python sets, and stays green with
  the native map stubbed.
- The step-1c instrument fixes have no gates: `policy_benchmark`'s apart grouping, `calibration_vs_oracle`'s APART
  strata and its `--json --jobs` merge, `prior_vs_oracle`'s strata, the `find_spec` guard.
(f) SMALL: the `--no-mappability` flag names a store whose mappability the index no longer reads (it opts out of
the splice blacklist); three unread scalars in `solve_kernel.cpp`'s positional aggregates; `pgzip` is an undeclared
optional import in the simulator; eight test sites import `rigel._resolve_impl` for names `rigel.native` re-exports.
- `_solve_impl` builds with `-g` only. Step 1b tried `-march=native` and IPO and reverted them: bit-identical on arm64,
  but on x86-64 they turn on FMA contraction in the ψ, transfer and solve kernels (clang: 0 FMA instructions at -O3,
  159 with -march=haswell) and move numbers (`TRAPS: a-kernel-edit-moves-ulps-under-fma-contraction`). RULED (owner,
  2026-09-28): they stay off; correct the two CMakeLists comments that are false for it. The 3bb35bdd message's
  "both proven bit-identical" is wrong for `_solve_impl` and stands corrected here.
- The numpy bit-matching machinery: the native solve still carries numpy's pairwise-sum order (`PW_BLOCKSIZE = 128`),
  the `lgamma` binding and named temporaries, with no Python twin left. Retiring it moves results by ulps and needs the
  goldens and the sweep replay re-recorded. RULED (owner, 2026-09-28): retire it after the strand-overdispersion and
  fl landings, in one golden re-record window.
Kept as coverage GAPS, not dead code: the five CLI command bodies, the silent policy through `calibrate`, the
simulator's sharded writers and its whole-genome grid, and the zarr splice blacklist.
(g) COMMENT AND DOC RESIDUE, one content-only sweep:
- Native: `bam_scanner.cpp`'s references to the deleted Python Fragment, "Fractional accumulator deposit", an
  orphaned reference, a lost rationale, the intron filter justified twice, `std::mutex` without `<mutex>`;
  `resolve_context.h` and `constants.h`'s history narration, "for compatibility" notes, a misplaced banner, four AF
  bits unused in C++; `scoring.cpp`'s `-Wreorder-ctor` warning and its gDNA per-hit comment omitting `n_cand > 0`;
  `fast_exp.h` listing three "call-site properties", one of them the function's own.
- Calibration: the sweep docstring's "counted in the kernel"; the test named for owned slots; "bound for the gates",
  although production reads `transfer_rows`; `region_geometry.py`'s "five populations" (there are four), the Axiom-0
  tell "nascent-RNA-active" (also in `signature.py`) and mature/phantom history; `result.py` citing the deleted
  `simplex_logodds._compose`; `calibrate.py`'s "two supports", a docstring that lost its conclusion, and the test
  helper `_zero_rna_opportunity` zeroing a field already 0; `strand_model.py`'s spliced model's qualification stated
  two ways, an undefined ε_CI, `posterior_95ci` not marked QC-only; the four texts stating "a quarter of the RNA side"
  as a general property (a deleted measurement on a rebuilt ladder).
- The EM, locus and CLI: `estimator.py`'s "every unambiguous fragment is deterministic" and a 102-character line
  hard-coding (N_t, 10); `locus.py`'s false (ref_id, start) sort claim (it sorts by categorical code) and its pointer
  to `calibration.priors` for the EM pseudocounts; `quant_command` omitting the `--annotated-bam` second pass;
  `report/substrate.py` ("omits the panels") and `report/__init__.py` (the extra is needed only for the charts).
- Tests: `test_estimator.py`'s "no EM at 0" (SQUAREM runs at least one step); `test_resolution.py`'s "(gap=0)";
  `test_buffer.py`'s "None for intergenic"; over-long test docstrings; "Boundary cases" headings from the
  edge→boundary rename.
- The simulator: `whole_genome.py`'s "…and filtering" banner; `splice_motif.py`'s over-long lines.
- Docs: DESIGN names the solver's `mrna_active`, which no longer exists, cites a test by its pre-rename name, and its
  compiler-flags row omits `-fno-finite-math-only` (as does `.github/copilot-instructions.md`); EQUATIONS §3.5g, DESIGN
  §3.5's "four ways", `accumulator.h` and `test_conserved_mass.py` still describe the deleted
  `region_contained_inv_opportunity_sum`, and one measurement table can no longer be reproduced; the executable
  specification and its C++ twin word one mass rule differently, and two payload comments point at inexact text
  (UNCERTAIN: check when touched); `test_chr.yaml` says the shadows are "≈ 5 % of RNA fragments" while the render
  realizes 0.76 %; the headers and TESTING list the YAML's extra keys without `blank:`; TESTING says `rigel sim`
  accepts `index:`, which the CLI now rejects; TESTING has no section for the VCaP read-name truth mix (A00839 exome
  DNA, HWI-D00127 transcriptome; its halves differ in fragment length, 307 / 244 bp, the confound the ladder forbids;
  read-name truth labels the RNA library's own gDNA as RNA; calibration can be scored only as slot-mass overlap).
- RELEASE-GATING: the MANUAL, README, report and cli texts still describe the deleted capture census or the old
  density-file semantics ("Capture on-target enrichment"; MANUAL's density row contradicting its correct neighbour);
  the MANUAL does not say that `sj_blacklisted` counts junction × record while the splice census counts fragments.
(h) SIMPLIFICATIONS, each proven a numeric no-op:
- `solve_one_block` and `transfer_prepare` each wire the chain, so the gates exercise the copy: one shared builder.
- Scoring is written twice, for multimappers and singles; unifying it also removes the p = ½ antisense disagreement
  (`ISSUES: latent-defects`).
- SQUAREM's VBEM and MAP branches (about 150 lines) differ only in the weight, the floor and the normalisation.
- `BamAnnotationWriter::write` re-implements `parse_bam_record`, and multimap detection is written twice.
- `TranscriptIndex.load` round-trips through frozensets: build the CSR directly.
- `intervals.feather`'s exons are read three ways: one walk would serve all three.
- Two functions named `contained_opportunity` sit in layer 2, and one ramp is written twice.
- `FragmentScorer` copies 13 fields into native: return the native scorer and store alphas.
- `_apply_scan_stats` reads 33 hand-listed keys with `.get(key, 0)`: iterate the dataclass and index strictly.
- The simulator's truth is spread over nine `ground_truth_*` methods and three file passes.
- Five instruments import `tests/calibration/_oracle.py` through `sys.path`; four scripts redefine `_shared.RUNS`.
(i) CONSTANTS — what the move into `config.CONSTANTS` left (2026-10-01; the layout is `DESIGN.md` §6):
- The native kernels' own constants are named and commented in their translation units (`em_solver.cpp`'s
  `SQUAREM_BUDGET_DIVISOR`, `transfer_rows.h`'s `EPS` and `TINY`, `solve_kernel.cpp`'s `PW_BLOCKSIZE`, the scanner's
  batch sizes), not gathered in `native/constants.h`, and the C++ inline literals are unaudited. Only the values Python
  shares are single-sourced (exported, never restated).
- `rigel sim`'s CLI restates `Scenario`'s defaults (genome length 5000, seed 42, 1000 fragments); reading them from
  the signature would import the simulator on every CLI start.
- `FragmentLengthModel.to_dict` rounds inline (2 and 6 decimals), beside the CLI's `SUMMARY_DECIMALS`.

### performance-memory-bounded-solve
`priority: next — a release gate for the two unparked numeric no-ops (owner, 2026-09-28); the rest Tier 4, the thread parked (owner, 2026-09-19: the method is the focus) · kind: build · 2026-08-17`
A deep run must be fast enough to iterate on, and memory-bounded. The port (the sweep as one native call) and the
work outside calibration are done (`DESIGN.md` §6b.15). THE BASELINE (VCaP, 18,568,456 fragments, `--threads 8`,
interleaved pairs, 2026-09-19): 125 s — calibrate 39.5 (four sweeps 24.0, landscape fits 4.1, its Python 10.3),
quant 33.7 (locus EM 14.4, capture effective lengths 5.2, scoring 5.7), the scan 28.4, the second pass 13.0, a
second fragment-length fit 2.8, the index load 6.4; peak RSS 10.3 GB, in quant.
OPEN, ranked (a cProfile share ranks, `TRAPS: a-profile-share-is-a-ranking`): ① the capture effective lengths
(quant's 5.2 s) — the short-template taper table that held most of it is deleted with the per-base rule, and the
one shared rule (`capture_eff_length.transcript_objects`, a closed form per object) is re-timed on the deep library
before it is ranked; ② the second pass's
factor combiner, 355,180 calls on 2–4-element arrays, hoistable to segment operations over the hypothesis CSR
(exactness to check); ③ the landscape fits, ~5 s, unexamined; ④ the sweeps' 24.0 s, where the blur's loop
interchange is priced at −1.5 s for a summation-order change and is the owner's call; ⑤ past 100 M fragments.
REFUSED as targets: the scan (at its floor for eight threads) and the index load (mostly native already).
⑥ THE REVIEW'S NUMERIC NO-OPS, each proven with `rename_identity.py --bam` from back-to-back pairs; UNPARKED (owner,
2026-09-28) for the first two, the rest batched:
- `get_loci_df` builds `locus_ids == lid` per locus, O(loci × transcripts): 0.48 ms per locus at 450k transcripts,
  about 19 s at 40k loci, in the CLI writer; `np.bincount` fixes it.
- `iter_chunks_consuming` holds the previous chunk while loading the next, and `_SpillWriter._run` keeps its last
  chunk while blocked: with chunks of up to 1 M fragments against a 2 GiB budget, a `del` fixes both.
- `build_region_partition_arrays` runs 3+ times per run, `build_sj_arrays` twice.
- `--annotated-bam` calls `add` per fragment where a batch API exists.
- The EM clears `local_map` up to the largest global index per locus, and allocates per E-step (heap rows for
  k > 512; the task list rebuilt).
- The scan path copies `gap_introns` per hit, one malloc each.
- The capture sampler's memo evicts across the mRNA and gDNA spaces, and its overlap lookup is O(T) per probe.
- The solve builds whole cubes for slots it then discards, and `pool_of` rebuilds the pool when the task count
  changes.

Instrument: `profiling/profiler.py`, `profiling/sweep_replay.py`.

### nascent-stress-sensitivity
`priority: later — Tier 2, first: a cheap re-measure of Tier 0's verdicts at the realistic nascent level · kind: question · 2026-08-22`
Does any in-scope verdict depend on the nascent stress level? The ladder runs `on_fraction 0.50`; realistic is ~0.10
(`DESIGN.md` §0b). MEASURED 2026-09-20 for the siphon: `g50 ss.99 ON` re-simulated at 0.10
(`~/Downloads/rigel_runs/suite/ladder_nrna_lo`) reads a siphon of +522,205 against +541,216 at 0.50, so that verdict
does not depend on the stress level — the repair is worth the full amount in the expected case. Re-simulate the worst
in-scope scenario at the realistic level and check whether any rank moves; a verdict that holds only at stress is a
robustness finding. The verdict to re-measure: at the realistic level (`g50 ss.99 ON`, 2026-09-21) the synthetic
nascent pool read 4.2× its truth — baseline 5.50 (nascent 4.2×), pseudocount-corrected 4.98 (6.3×), oracle 1.85 (2.2×)
— before the one shared length rule and the count-form pseudocount. `sim/panel.py`, `quant_accuracy.py`,
`policy_benchmark.py`.

### psi-reads-kappa-where-the-strand-channel-is-dead
`priority: later — Tier 2, the first mechanism: in scope, unstranded × capture OFF · kind: defect · 2026-09-28`
ψ reads the estimated κ̂ where the protocol's Bayes factor calls the strand channel dead: the discriminability gates
only `strand_evidence`, and `solve_kernel.cpp` builds ψ's grid with κ̂ unconditionally. At od = 0, κ̂'s sampling error
gives the strand term a slope n(½ − κ̂), so real overdispersion reads as composition along it, at a weight of about
n/m for m spliced reads (SD(κ̂) ≈ 0.08 at about 40 spliced reads). On odg05 ss 0.50 g98 (true od 0.05, κ̂ = 0.515),
assuming od 0 costs +1,178 in scope at capture OFF and +29,106 on the deferred capture-ON row, and nothing against a
true 0 (test panel, κ̂ 0.485). Not dissected. The candidate — ψ reads κ = ½ where the channel is dead — is its own
A/B; od itself is 0 since 2026-10-01 (`ISSUES: strand-overdispersion-one-shared-value`).
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/01_fallback.md` (its last open question);
`policy_benchmark.py --by-class`, unstranded rows.

### the-fl-boundary-inversion-reads-missing-evidence-as-a-value
`priority: later — Tier 2, the fl second wave · kind: defect · 2026-09-28`
Missing evidence must drop out, never be replaced by an invented value; after the RNA-counts fix
(`ISSUES: the-realized-gdna-length-law-reads-rna-counts`, closed 2026-09-30), five places in `fl.py` still break this.
(a) UNMEASURABLE INTRON REGIONS. A region too short to hold a fragment has e_r ≈ 0, so its RNA density reads 0 rather
than unknown, a_b → 1, and nascent exon–intron–exon fragments are priced as capture. At g05 OFF after the fix the
signed w(ε − 1) reads +0.66 against a rectified +0.70 (bias); at g50 −0.006 against +0.072 (rectification only).
Exon|intron pairs with ρ_adj = 0 hold 10 % of boundary counts at a count-weighted ε ≈ 20–21.6; pairs with ε > 5 border
introns of median 94 bp (31–185). The likely cause of the +4 bp g05 ON overshoot. OWNER: carry the variance (the
weight goes to 0) or borrow from the gene's measurable introns? One derivation could also serve
`ISSUES: the-efficiency-posterior-floor-on-empty-pieces`.
(b) ZERO gDNA. The one-sided rate declines only at total_n = 0, so at g00 ρ_off is a geometric floor (2.1e-6, against
5.3e-3 at g05): ε reads about 2,300·n_b wherever ρ_adj = 0, and the realized law turns 99 % "capture excess" (+28.9 bp
over uniform after the fix, +19.1 before: the fix makes it worse). `if den[0] <= 0 or den[1] <= 0: break` drops the
whole boundary stratum at g00 (den = (185k, 0)), since simulated RNA never crosses an exon|intergenic edge. The
one-class limit would admit a2·n2 RNA crossings against m_C = 80; derive it once for both den and n. OWNER: how should
ε fade? It is a floor, not noise, so a variance term may not cure it.
(c) THE CONTAINED PAIR still declines on an empty pool in the four g00 rows, falling back to the capture-selected
four-pool census (243–244 bp under capture; at g00 ON the realized law reads 263 bp against 237, RNA 227). Left open
by `ISSUES: the-gdna-length-law-falls-back-at-identical-purities`.
(d) INTERGENIC RNA. The off-target rate is fitted on intergenic and intron counts, so RNA spread over them lifts it:
toy fixed point 202 bp against 217 (a3 0.63 against a true 0.37); noise-free off-capture steps −10.0 / −15.2 bp at
ig_rna 0.5 / 2. The fl counterpart of `ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool`, in another
estimator (`ISSUES: capture-degeneracy-standing-risk`). Ladder RNA never reaches intergenic space, so its real-data
size is unmeasured.
(e) AN EMPTY CROSSING POOL is structural. `den` counts every crossing but n2 and n3 only single crossers, so den > 0
with an empty pool is reachable (n3 = 0 at g00), and `_normalized`'s uniform 0–1000 law then enters r̂ (n2 = 0) or
becomes g_B and pushes μ_g toward 500 (n2 = n3 = 0). No pool means no boundary stratum; one empty pool skips the
inversion. Gates, both failing today: den > 0 and n2 = 0 gives g_B = the normalised f3; n2 = n3 = 0 leaves μ_g
unmoved and m_B = 0. It moves numbers: on the capture-ON ×30 toy the on-target share goes 0.042 → 0.029.
Instrument: `06_fl.md` and `scratch/fl_loop2.py` under `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/`.

### the-fl-boundary-inversion-has-underived-pieces
`priority: later — Tier 2, the fl second wave, with the missing-evidence entry · kind: derivation · 2026-09-28`
Five pieces of `fl.py`'s boundary inversion are borrowed, not derived.
(a) The one-sided excess clip (cnt − ρ_off·e_g)₊ makes a pure-gDNA region's expected share read below 1 at every
finite depth: E[a_b] = 0.914 / 0.779 / 0.853 / 0.910 / 0.965 at λ = 0.1 / 0.5 / 1 / 10 / 100; real introns average
0.910 (LBX0588), 0.940 (MO_3021), 0.978 (LBX0588 at 10 %). a2, a3 and m_B carry the bias, and the flatness premise
a3 = 1 never holds on a real library; in the toy it moves the realized mean by at most 0.03 bp. Gate: a Poisson
pure-gDNA fixture at λ ≈ 0.5–1 bounding a2 against the exact E[a_b].
(b) The boundary pair's standard error, se_b = √(a2²/n2 + a3²/n3), is borrowed from the contained pair and was never
derived for the incidence-weighted a2 = Σa_b·n_b / Σn_b; it sets the inversion's resolution weight.
(c) The boundary's RNA side is priced at the spliced census mean, not at its own RNA length law; under a capture that
prefers some lengths E_g[s]/E_r[s] need not cancel. Unmeasured, and distinct from
`ISSUES: the-scorer-reads-a-census-length-law`.
(d) UNCONFIRMED: the inversion mixes by incidence shares (n2, n3) while the de-tilted pools mix by start density. A
first draft read the fixed point's |Φ′| falling 0.026 → 0.005 on switching, and no results file exists. Measure first.
(e) THE REFRESH'S 0.25 bp STOPPING TEST. One unconditional refresh replaces it, as straight-line code with f2, f3, n2,
n3 and den out of the loop. The test fires 98 / 95 / 75 % of the time at 260 / 2.6k / 26k crossings; the secant gives
q ≈ 0.03, so one refresh lands within 0.011 bp of the fixed point (g05 ss.99 ON: 244.569 / 244.965 / 244.954 at zero
refreshes / one / the fixed point), 30× below the 0.32 bp replicate spread. Unpriced: off capture the refresh swaps
the contained pool's mean for the noisier crossing pool's.
Instrument: `06_fl.md` under `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/`.

### capture-blind-gdna-divisor
`priority: later — Tier 2, the fl second wave · kind: defect · 2026-08-31`
`gdna_opportunity_from_index` is computed from the index alone, so under capture it removes ~6 bp of a ~30 bp
length selection — the gDNA control moved +6.0 % on all six capture-ON rows (gDNA has no introns to miss), and with
`ISSUES: eb-shrinkage-magic-ess` it owns the −5.90 % capture-ON length ceiling. `capture_eff_length` already models
the panel; it also blocks `ISSUES: crossing-pool-contrast`. ⛔ Not priced on the deliverable: no `quant_accuracy`
arm reaches the locus gDNA length (the oracle overrides only calibration's count arrays;
`ISSUES: end-to-end-error-unattributed` records why no length arm prices it). A price needs the simulator's per-locus gDNA opportunity
as the truth.

### eb-shrinkage-magic-ess
`priority: later — Tier 2, the fl second wave: replaced on the re-simulated fl-gap panels · kind: defect · 2026-08-31`
`CONSTANTS.fragment_length.pool_prior_ess` (1000) shrinks the gDNA pmf toward `global_pmf` (mostly RNA whenever gDNA is a minority) at a
magic ESS: inert on the ladder (0.01 bp), dominant on an fl-gap arm at `g05` capture-ON (`ship−pool` −23.7 of −31.7
bp, measured on an fl-gap panel since deleted; not re-recorded, the instrument being retired). Replacement: weight the pools by their own precision, which is not yet derived (the strand od's precision reconcile was
refused for mismatched inputs, `ISSUES: the-strand-overdispersion-reconcile`). Its instrument, `fl_pool_purity.py`,
was retired 2026-09-14 (in git). The ladder's deliverable cannot see it; only the fl-gap arm can. The realized-law fix
(`ISSUES: the-realized-gdna-length-law-reads-rna-counts`) passes the EB-smoothed RNA pmf (owner, 2026-09-28), so the
ESS also reaches the boundary composition, shifting μ_r at small N_s by ess/(N_s + ess)·(μ_anchor − μ_RNA). Plain
normalisation is no alternative: at N_s = 0 it invents a uniform law (a2 0.109 against a truth of 0.143). In the toy
a2 moves only 0.1472 → 0.1482 as N_s goes 0 → 1e6, and on the ladder EB and plain agree to 0.002 bp; N_s ≈ 1e2–1e4 is
unmeasured. The replacement is measured on the re-simulated fl-gap panels.

### intron-seeds-near-probes-are-capture-enriched
`priority: later — Tier 2, the od research tail; answered before any purity weight or background repair that treats intron regions as off target · kind: question · 2026-09-28`
Probes overhang exon edges and bind the bases there, so the intron regions next to probes are capture-enriched and
"every region seed is off target under capture" is false. On odg05 g05 capture ON truly pure intron seeds read
n/(λ_off·E) at a median 1.89, 90th percentile 8.2, 99th 44 (1.00 / 1.23 / 1.84 at capture OFF): a pooled 1.9× λ_off.
Two places treat introns as off target: `density_deconv`'s background (introns at the intergenic depletion); fl's
one-sided rate (intergenic and intron regions together). Proposal: measure enrichment against distance to the
nearest probe, then treat near-probe introns as enriched or price the loss.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/04_evidence.md`.

### strand-likelihood-over-confident-beyond-od
`priority: later — Tier 2, the od research tail; re-read after strand-overdispersion step 1 · kind: question · 2026-09-28`
A nonzero RNA overdispersion beats the simulator's planted RNA od of 0 on several strata, and the 0.2 value lowers
false gDNA on the zero controls, so the RNA strand likelihood looks over-confident even where its od is truly 0. Its
variance is frozen per slot at the reference composition (`psi_kernel.h`'s slot cube, `transfer_rows.h`): a
misspecification in the strand term, outside the od estimator. MEASURED (the od harness, 2026-09-27): on odg05
stranded ON the joint fit beats the planted pair on calibration (172,329 against 179,181) and transcripts (−28,829
against −11,248); stranded OFF the 0.2 value beats it on calibration (50,971 against 52,004); a planted RNA od of 0
costs +57 % at g98 ss 0.70 ON (184,431 against 117,753); on the ladder's g00 the 0.2 value cuts false gDNA 1,313 →
1,185 (−9.8 %). No mechanism identified; re-read after the joint fit lands.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-27_strand_od/`.

### the-pseudocount-strength-is-not-derived
`priority: later — Tier 2, the exact VB update is a ready A/B; the odds are fixed (DESIGN.md §3.1c) · kind: open question · 2026-09-24`
The EM's two pseudocounts carry `P_g + P_R = U`, one pseudo-fragment per unit with a gDNA candidate — the total
the replaced rule carried, kept so the count form changed only the odds. Nothing derives it. The rule that
counts each fragment once is refused as buildable (`ISSUES: a-count-once-density-prior-for-the-strength`): under
capture calibration holds no estimate from outside the locus. Under MAP the strength cancels at the neutral
point; under the shipped VBEM it does not, and the grouped update lets it move the split WITHIN RNA — a
6-fragment isoform reads 4.0 / 5.0 / 5.8 at `P_R` = 0 / S / 10S. The exact VB update for the same prior holds it
at 4.0 at every strength and needs no constant; it is its own A/B
(`~/Downloads/rigel_runs/prototypes/2026-09-24_spliced_rule/toy_vb_nested.py`). Why a prior is needed at all:
with none the EM's gDNA error is −12.1k / −119.7k / −731.8k at `g05` / `g50` / `g98 ss.99 ON` and −0.4k off
capture, and −35.5k remains at `g98 ss.99 ON` under an exact calibration — a lean of the capture likelihood
toward RNA that no strength may be chosen to cover. The owner's framing (2026-09-24): calibration gives the
locus's gDNA FRACTION, and a multiplier converts it to pseudocounts; that multiplier depends on the number of gDNA
candidates and their likelihoods. Today it is the count of units whose gDNA candidate survives pruning (`U`), and
the pruning (`DESIGN.md` §3.1d) is what stops a unit with a vanishing gDNA reading from counting.

### per-transcript-prior-lane
`priority: later — Tier 2, research: the largest in-scope lever, no candidate mechanism yet; the owner's problem after the junction price (2026-09-22, 2026-09-23) · kind: build · 2026-08-31`
`rna_prior_weight` is PLUMBED end to end — `pipeline.py` → `estimator.py` → `em_solver.cpp` — and ⛔ NOTHING IN
`src/` FILLS IT, so the solver takes the fallback `w_i = raw[i]`, which echoes the EM's own belief and cannot
contradict it. ⚠ The lane is ONE static array, so filling it reallocates the WHOLE RNA pseudocount (2,793,710
fragments at `g50 ss.99 ON`, as large as the unspliced RNA itself), and a coverage-derived weight is a far worse
isoform allocator than the EM's own likelihood (5.21 → 53.37 %); a weight that keeps `raw[i]` where the data cannot
speak needs `raw[i]`, which lives in the kernel.
THE CEILING (`quant_accuracy.py --arm oracle_alloc_seed`, the true weights mature and nascent, the committed tree
`c52c9b93`, fractional, transcript Σ|Δ| as a share of the true RNA at `g00` / `g05` / `g50` / `g98`): stranded OFF
2.01 / 1.59 / 2.49 / 15.70 → 0.43 / 0.43 / 0.70 / 6.38 %, unstranded OFF 1.73 / 1.93 / 2.55 / 18.35 → 0.42 / 0.45 /
0.81 / 8.35 %, stranded ON 6.53 / 3.49 / 5.50 / 31.03 → 2.91 / 1.08 / 2.07 / 18.00 %; the deferred stratum 7.32 /
11.00 / 10.11 / 90.63 → 3.41 / 2.67 / 4.03 / 221.15 %. The largest lever on every in-scope stratum. A capability
proof, never headroom: it hands over the true support, and a zero weight is absorbing.
WHERE ITS GAIN SITS (2026-09-20, the four `g05`/`g50` `ss.99` rows): not in support — an exclusive-junction support
rule recovers ~0 % of it (`ISSUES: exclusive-evidence-support-rule`, refused) and the TRUE support only 3–28 % —
but in the PROPORTIONS among expressed isoforms that share every object: 49–61 % of the gain is on isoforms with no
exclusive sj, region or boundary (94 % of them mosaics), 65–76 % in genes with 11+ isoforms. Two weightings are
refused (`ISSUES: refused-transcript-weights`, `ISSUES: refused-soft-min-path-weighting`). The next candidate moves
proportions among mosaic isoforms or is a sparsity argument on silent mosaics. On the one shared length rule the
synthetic pool's error at `g50 ss.99 ON` is the pseudocount's
(`ISSUES: the-pseudocount-prior-is-biased-toward-gdna`), not a siphon's residual
(`ISSUES: nascent-siphons-gdna-under-capture`).
RULED (owner, 2026-09-28): the lane is kept, and with it `warm_start='prior'`, plumbed alongside it, whose only
producer is `quant_accuracy.py`'s oracle-allocation arm; neither is dead code for a sweep to delete.
Instrument: `quant_accuracy.py`, per stratum.

### message-layer-open-cases
`priority: later — Tier 2, the message layer on unstranded capture OFF · kind: question · 2026-09-09`
Four residual cases, none a hole (`DESIGN.md` §6b.4–§6b.14): (a) an exon with both faces speaking — the
arrivals are summed and the price where both carry the same intron's claim is unmeasured; (b) the intron
factory on a region with both exon and intron bits, under capture; (c) the chain of termini — an empty
outside piece (median 12 bp) whose far face is another terminus, half the ladder's terminus-boundary error,
every upper side refused (`ISSUES: the-edge-upper-side`); (d) substrate `nest` (the prior already serves it,
0.666 vs 0.630) and the antisense's nascent variant (`docs/TESTING.md` §0a; `div` was built 2026-09-14).
`policy_benchmark.py --by-class`.

### refit-vs-message-arbitration
`priority: later — Tier 2, the message layer on unstranded capture OFF · kind: design · 2026-08`
At the unstranded × capture-OFF exon cell the refitted gDNA prior and the message both impute one slot with
nothing arbitrating them; the message is the accurate voice there and the refit displaces it. Re-read under the
E-step (with an instrument since retired): the prior does the unstranded rows and the messages the stranded
capture-ON ones. `policy_benchmark.py --by-class`, `calibration_vs_oracle.py`. Belongs with `ISSUES: gdna-landscape-trains-on-false-positives`.

### the-intron-own-solve-on-unstranded-capture-off
`priority: later — Tier 2, the message layer on unstranded capture OFF (the plan of 2026-09-28; parked 2026-09-05); re-measure before ranking · kind: question · 2026-09-28`
On unstranded data off capture an intron region's own solve has no strand channel, so its gDNA fraction is whatever
the landscape prior and the message layer deliver. Measured 2026-09-05 at 43–45 % of the in-scope off-capture error
under every policy; that predates the location floor, the E-step and the one shared rule, so re-measure before
ranking. Related: `ISSUES: refit-vs-message-arbitration`, `ISSUES: two-sided-exon-row`.
Instrument: `policy_benchmark.py --panel ladder --by-class`, the unstranded capture-OFF rows, intron class.

### the-efficiency-posterior-floor-on-empty-pieces
`priority: later — Tier 2: the unprobed class's scale at low gDNA, its EM cost unmeasured · kind: defect · 2026-09-23`
An object's capture efficiency is the posterior mean of its clipped gDNA density under the landscape, from its own
count (`capture_efficiency`). On an object whose support holds far less than one expected gDNA fragment the count
cannot tell a depleted object from a partly captured one, so the posterior reads above the truth. It is what lifts
the unprobed multi-exon transcripts at `g05 ss.99 ON` to L / Y 2.84 (sd of log 0.69, 706 transcripts) against the
synthetic spans' 1.148 and the gDNA footprints' 1.135 in the same probed class (1.56 against 1.005 / 1.013 before
the one shared rule; at `g50` 1.214 against 1.104 / 1.096). MEASURED 2026-09-23: it is the transcripts' PIECES, not
the junction sum — with the pieces at the truth the class reads 1.05, with the junctions at the truth it is
unchanged. The unprobed transcripts' exon pieces are short (median contained support 4.2 starts, ~0.001 expected
gDNA fragments each at `g05`) and read 2.57× their truth from the simulator's noise-free counts and 4.87× from
calibration's (support-weighted), where the spans' long pieces (support ~97) read 1.07× / 1.05×; at `g50` the same
exon pieces read 1.16× / 1.18×. The short-intron fallback of `ISSUES: the-junction-price-is-noisy-within-a-gene` is
the same posterior on the same kind of object. What it costs through the EM is unmeasured: the class competes with
the gDNA component in unprobed genes, where every fragment is off target. A floor or a threshold on support is a
constant; the open question is what an object holding under one expected fragment should read — the question of
`ISSUES: the-gdna-landscape-collapses-at-low-depth` on capture-ON transcripts, and of fl's unmeasurable regions
(`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value`). `ruler_vs_truth.py
--scale` (the unprobed rows); the census in `~/Downloads/rigel_runs/prototypes/2026-09-23_junction_price/`.

### scoring-penalties-are-underived-constants
`priority: later — Tier 3, with the splice work's h-weighted reading on the cluster · kind: question · 2026-09-28`
Two fragment-score penalties are 0.1 and underived (`scoring.py`): log(`mismatch_alpha`) per NM mismatch, and
log(`overhang_alpha`), RNA only. Both also set where the pruning lattice cuts. The splice-artifact work adds more to
derive or put to the owner: the minimum count of 2, the anchor floors 4 / 6 / 11, `alignSJDBoverhangMin` 8, the 5-bp
mismatch window. The ε of the proposed h-weighted reading (`ISSUES: splicing-artifacts`) must come from the library,
never from `mismatch_alpha`. A zero alpha gives NaN today (`ISSUES: latent-defects`).
Instrument: none yet; the splice truth mixes.

### multimapper-intergenic-alignments
`priority: later — Tier 3, on the aligned ladder after the splice phases; a separate feature after the end-to-end work (owner, 2026-09-24) · kind: defect · 2026-09-24`
A multimapper whose alignments all lie outside genes is intergenic and counted as gDNA — correct. One with an
alignment in a gene is solved by the EM, and the owner's rule is that it is compatible with gDNA if ANY alignment
is unspliced. The scanner never buffers a multimapper's alignment that overlaps no transcript
(`bam_scanner.cpp:1808-1845`: only unique mappers are appended), so a gDNA read from an unannotated processed
pseudogene whose other alignment splices across the parent gene reaches the EM with its spliced alignment alone
and is certified RNA on the parent. The owner expects multimappers to matter most. The repair buffers those
alignments, lets the scorer take a gDNA term from an alignment with no transcript candidate
(`scoring.cpp`'s `n_cand > 0`) while an all-intergenic group still drops silently, keeps the fragment ledger
exact, and decides which locus's gDNA component owns a gDNA reading at an intergenic position, outside that
locus's gDNA length. Unmeasurable on the ladder (`NH` 1); the aligned ladder
(`~/Downloads/rigel_runs/aligned_ladder/`) and `tests/scenarios_aligned/test_multimap_counting.py` are where it is
seen. Taken with `ISSUES: multimapper-blind-support`. Also open for multimappers: the group's gDNA term is the MEAN
over its eligible alignments (`scoring.cpp`, pinned by `tests/test_pipeline_routing.py`), where uniform gDNA
placement derives the sum.

### multimapper-blind-support
`priority: later — Tier 3, a deferred session with ISSUES: multimapper-intergenic-alignments (owner, 2026-09-24) · kind: defect · 2026-09-16`
The accumulator drops every fragment with `NH > 1`, so an object whose gDNA fragments are multimappers holds no
count the calibration can see, and its capture efficiency — read from its own gDNA count against the fully captured
level — reads it at the DEPLETED level: the transcript is contracted as if unprobed. Measured read-only on the two
captured real libraries on the per-base rule (the session's `multimapper_check.py`: per transcript with ≥ 20
fragments over its pieces, the share of multimapping reads among them from the BAM's `NH` tag, then the length
factor's quantiles per share bin): LBX0588 (reference 2.977e-01/bp, 11,579 kernels) median factor 0.076 / 0.003 /
0.013 / 0.001 at a share below 1 % / 1–10 % / 10–50 % / ≥ 50 % (n 36,446 / 6,161 / 1,429 / 928), the share below
1/10 rising 55 → 90 %; VCaP (8.593e-02/bp, 24,496 kernels) 0.231 / 0.149 / 0.028 / 0.002 (n 280,227 / 35,494 /
2,237 / 1,103), the share below 1/10 36 → 95 %. A hundredfold gap between the unique-mapping and the repeat-rich
transcripts, monotone in the share, on both libraries; the two libraries without a reference read every transcript
at 1 and are uninformative. The confound is real and unmeasured here — a repeat-rich transcript may also be
unprobed by design — so the number is an upper bound on the blindness, not its size; the floor `C/(C+1)` the repair
deleted had hidden it behind a +3.4 nat bias on every unprobed transcript
(`ISSUES: ruler-multimapper-floor-caps-the-correction`, CLOSED). The repair is in the support, not the posterior:
an object's support (a region's contained, a boundary's crossing) should count only the starts a fragment could be
uniquely placed at — a mappability of the support from the index, so that a wholly repeated object has no support,
no evidence, and reads the population's clipped mean rather than the depleted level. Instrument:
`ruler_vs_truth.py` cannot see it (the simulator's reads are unique); the truth on a real library is the
multimapper share against the factor as above, and the repair's gate is that the four bins read alike.
RULED (owner, 2026-09-24): "add multimapping deposit to the accumulator and let calibration include multimapping fragments, which is correct behavior and something we have intended to add".

### unwitnessed-loci-and-multimappers-at-the-em
`priority: later — Tier 3, on the aligned ladder with the multimapper pair; the ladder cannot show either case, and no instrument is in the tree (the 2026-09-21 scratch census) · kind: defect · 2026-09-21`
MEASURED 2026-09-21: every real library carries MultiLoci with a gDNA candidate and no witnessed gDNA (`P_g = 0`:
257 / 237 / 288 / 194 loci holding 374 / 723 / 941 / 2,182 gDNA-bearing units on LBX0190 / LBX0588 / MO_3021 /
VCaP; `L_g = 0` nowhere), and a multimapper share of the gDNA-bearing units of 13.3 / 4.0 / 13.2 / 1.0 %,
concentrated in the sparse loci (per-locus median 0, q90 1–5 % over loci with ≥ 20 such units). A kernel that
pins gDNA at calibration's count sends every unspliced unit of an unwitnessed locus to RNA and under-calls
gDNA by the multimapper share one for one; the ladder has neither case (0 of 1,213 loci). A per-fragment prior
prices each multimapper placement at its own objects, which is the cleaner rule.

### vcap-dna-only-objects-lose-gdna-at-full-depth
`priority: later — Tier 3, local real data with truth by read name · kind: defect · 2026-09-28`
At full depth the VCaP mix reads P/O 0.922, and the deficit sits where no truth label is needed: P − O = −365,700,
+17,267 on structural objects and −382,967 elsewhere — −345,473 on no-RNA objects (−15.0 % of their O, 7.35 points of
the library), −58,503 on mixed, +21,009 on no-gDNA. It is already in pass 0 (P/O 0.930). By object size: [1,3)
−26.5 %, [10,30) −12.2 %, [30,100) −11.4 %, [100,1k) −5.0 %; the c/√n floor at c = 0.39 predicts 5.5 % at n = 50 and
2.3 % at n = 300, so objects of 30–1,000 fragments lose about twice that. Mechanism not isolated;
`ISSUES: the-gdna-landscape-collapses-at-low-depth` excludes it.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/pages/explanation_depth.md`; the VCaP-mix scorer
under `v2/scripts/`.

### em-gdna-exceeds-calibration-on-the-vcap-transcriptome-half
`priority: later — Tier 3, local real data with no truth · kind: question · 2026-09-28`
On the VCaP transcriptome half (no added gDNA) the EM assigns more gDNA than calibration called, under every od arm,
and about 2.5× the gDNA prior it was handed. At full depth: calibration 86,650 (od 0), 81,612 (shipped), 69,004 (the
joint fit or the 0.2 value); the EM 96,362–97,856 (0.695–0.705 %); the EM's gDNA prior sum 33.1k–39.8k. Calibration's
`library_gdna_fragments` includes intergenic regions and the EM's total excludes `intergenic_total` = 9,521, so like
for like the EM plus intergenic exceeds calibration by about 24–53 %. Across arms the EM shrinks calibration's 20 %
od-driven swing to 1.5 % (transcripts move 1.8–2.8k of 13.4M). The half is not certified gDNA-free, so neither number
can be scored.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-27_strand_od/results/robust_no-gdna/`.

### instrument-ledger
`priority: later — Tier 4, batched between A/B windows · kind: instrument · 2026-09-28`
What the instruments do wrong or leave out; each fix is its own commit, with its `--self-test` or a gate.
(a) `slot_truth.npz` counts boundary gDNA per crossing, in incidence units, so under capture its sums exceed the
library's true gDNA (slot-frame O 1,281 / 12,485 at g001 / g01 ON against 1,170 / 11,700 fragments; P/O 1.015 in the
slot frame against 1.055 in fragments), while `policy_benchmark.py` calls its Σ|Δ| "in fragments". Document the frame.
(b) No instrument output records the code revision, the `--set` overrides or the prototype arm
(`TRAPS: an-ablation-that-never-ran`).
(c) `sweep_replay.py` pickles slots dataclasses; on Python 3.12 an old capture loads by position, shifted silently, so
every capture taken before 1747b282 must be retaken.
(d) `build_test_reference.py` puts `src` first on `sys.path` unconditionally; fix it before the contaminated panel's
in-tree stage.
(e) The ladder-report skill hard-codes the four ladder rungs (StopIteration on g25 or g001).
(f) `panel.py --index` reaches the cache and oracle stages but not `simulate`.
(g) `calibration_vs_oracle.py` has no arm that carries the true capture efficiencies, so no calibration A/B sees the
capture-contracted length move; `quant_accuracy.py --arm oracle_ruler` is the model.

### debug-capture-memory-is-unbounded
`priority: later — Tier 4: the fix is in src; the one-real-genome-job rule stands meanwhile · kind: instrument · 2026-09-28`
`calibrate(_debug=…)` builds a `SweepCapture` for every sweep with no size or region bound; on the whole human genome
it reached 25 GB in 20 s and nearly crashed the machine (2026-09-28). `ruler_vs_truth.py` always passes `_debug` and
does not refuse a real index, and no capture is region-restricted or streaming.
The per-session memory guards miss Python from other environments, processes under 1 GB and growth between 3-s polls.
For scale: a normal quant peaks at about 3 GB, and `rename_identity.py --check` on one ladder condition holds 8.9–10.7
GB. RULED (owner, 2026-09-28): `_debug` refuses an index whose manifest does not mark a panel, and the capture is
bounded to named regions. The one-real-genome-job rule stands meanwhile.
Instrument: none in the tree; the session's guards (`mem_guard_env.sh`, `mem_guard_fp.sh`).

### expand-the-gdna-spectrum
`priority: later — Tier 4, the panels · kind: decision · 2026-08`
Fill the gDNA spectrum (1, 5, 10, 25 % up past 90) without multiplying benchmarks: a level is justified by a measured
transition and crosses a reduced set of the other axes until an interaction is shown; each condition costs a simulate,
two caches and a certification (`sim/panel.py`). The junction-probed test-chromosome twin, stale since the
capture-physics change, was re-simulated on the current simulator, cached and certified on 2026-10-03 (owner; the
fragment-length campaign reopened), which supersedes the 2026-09-28 ruling to retire it. The depth family's unstranded
rows are ruled in `ISSUES: the-gdna-landscape-collapses-at-low-depth`.

### a-pure-gdna-library-reads-as-nascent-rna
`priority: parked — Tier 5, the deferred stratum; kept as an open challenge (owner, 2026-09-28) · kind: defect · 2026-09-28`
Rigel reads a library of pure genomic DNA as 93 % RNA. The case is the exome-DNA half of the VCaP mix, 4,671,916
fragments, separated by read name. It is unstranded (strand specificity 0.5008) and capture-enriched (calibration's
gDNA reference density is 180× its genome-wide density).

Where its fragments went:
- synthetic nascent RNA, over 53,962 spans: 3,781,970 (81 %);
- annotated multi-exon mRNA: 452,338 (9.7 %), only 4,213 of which come from spliced fragments;
- annotated single-exon transcripts: 109,003 (2.3 %);
- gDNA: 313,624 (6.7 %).

The library has 3,225 fragments spliced at an annotated junction (0.07 %), about what alignment artifacts alone
would give.

This is the deferred unstranded × capture-ON stratum at its extreme. The gDNA fraction cancels from the strand mean,
and capture depletes the intergenic space that would otherwise measure the gDNA density, so a gene span's unspliced
coverage fits nascent RNA as well as gDNA.

What could still tell them apart is the SPLICED FRACTION. Mature RNA at a multi-exon transcript produces
junction-crossing fragments at a rate its geometry and the fragment-length law set (layer 2's `sj_opportunity`).
Here the multi-exon mRNA assignment carries under 1 % spliced support where real RNA carries tens of per cent. And by
the nascent scope ruling, nascent RNA is sparse beside mature RNA, which it cannot outweigh 7:1 as it does here.
Whether either becomes a likelihood term, and how a pure-gDNA library should be recognised, is open.

Its plausible in-mix face: inside the VCaP mix the EM sends the true DNA fragments to RNA at 1 % / 10 % / full depth
(sampled draws, seed 0, not fractional) as synthetic spans 0.254 / 0.185 / 0.099, other synthetic RNA 0.014 / 0.008 /
0.004 and mRNA 0.118 / 0.059 / 0.044: 39 % / 25 % / 14.7 % in all.

Measured with the strand-overdispersion harness, `rigel quant` under every arm (the arms change nothing here):
`~/Downloads/rigel_runs/prototypes/2026-09-27_strand_od/results/robust_no-rna/`.

### the-pooled-q-in-the-gdna-count
`priority: parked — Tier 5, lands with ISSUES: the-scorer-reads-a-census-length-law, parked with it (owner, 2026-09-30) · kind: defect · 2026-09-24`
Calibration's per-locus gDNA count (`priors.assemble_priors`) converts each boundary's gDNA mass to fragments by
`boundary_mass_per_crossing`, pooled over gDNA and RNA, so where RNA dominates a boundary gDNA is over-stated even
under a perfect calibration (`prior_vs_oracle.py` O − S: +33.7k at `g50 ss.99 ON`, +2.9k at g98 ON, +1.9k off
capture). The EM passes a g50 count change through at only 0.37–0.50, so its size in the result is about 12k: with
gDNA's own share the g50 gDNA pool moves −18.2k → −30.4k, spans 142.0k → 152.1k (truth 150.4k), transcripts
unchanged. It cancels part of the capture likelihood's lean toward RNA today, so it lands with those repairs. The
gDNA LENGTH already uses gDNA's own share (`ISSUES: the-pooled-q-in-the-gdna-length`, replaced for the length only).
`prior_vs_oracle.py` (P, O, S) is its one instrument and retires when it lands (owner, 2026-09-28).

### unannotated-transcription-is-booked-as-gdna
`priority: parked — Tier 5: a DESIGN ruling after 0.8.0, recorded now (owner, 2026-09-28) · kind: design · 2026-09-28`
By Axiom 0 an unannotated slot admits only gDNA, so real unannotated transcription is booked as gDNA twice — by
calibration (structural objects give P = count) and by the output's intergenic rule — and it biases any pooled
structural rate a landscape repair might read. DESIGN §0b says the tool does not model it and calls that safe for a
pooled level; `ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool` covers one route only. MEASURED:
`test_blank` books 8,844 of 8,844 RNA fragments as gDNA at g001 capture OFF and 8,889 of 8,889 at g01, its eight
shadow transcripts all of the false-negative mass (12.9 % / 11.5 % of transcript Σ|Δ|); VCaP's structural P − O at
full depth is +17,267 — intergenic +8,654 (12.4 %), gene edges +8,613 (14.1 %) — at a median minor-strand share of
0.000 (74.5 % of 208 objects at ≤ 5 %); on the test chromosome's g00 the false gDNA is 100 % shadow transcription.
Admitting RNA at unannotated slots on stranded data needs a ruling and has costs: the closed floor would point toward
RNA on truly pure intergenic objects, and that RNA would have no transcript to go to. RULED (owner, 2026-09-28): for
0.8.0 only the output's spliced intergenic count changes (`ISSUES: the-format-changes-to-batch-before-release`).
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/results/blank_arm.txt`,
`pages/explanation_depth.md`.

### ruler-witness-geometry-on-transcript-panels
`priority: parked — Tier 5, deferred past 0.8.0 (owner, 2026-09-20): Rigel stays panel-agnostic and takes no panel input · kind: defect · 2026-09-15`
Only a transcript holding the junction a probe spans binds that probe whole, so the extra capture of junction-spanning
fragments is ISOFORM-SPECIFIC, and the EM splits a gene's shared fragments by the ratio of its isoforms' capture-aware
lengths. gDNA holds the probe's parts apart and binds the better one, and at zero gDNA there is no witness at all; the
shipped junction price reads it from gDNA at the junction's low and high boundaries by conservation of bases, with no
panel input — 0.959 of the simulator's own junction capture over 395 junctions of 60 probed transcripts in aggregate,
but per junction too noisily (`ISSUES: the-junction-price-is-noisy-within-a-gene`). MEASURED on the rebuilt ladder
(the simulator binds every fragment through one contiguous part of a probe since 2026-09-19; earlier numbers measured
an asymmetric half-match the owner ruled unphysical), stranded × capture-ON transcript Σ|Δ| as a share of the true RNA
at `g00` / `g05` / `g50` / `g98`, fractional: the per-base rule 6.53 / 3.49 / 5.50 / 31.03 %, the simulator's own
lengths (`quant_accuracy.py --arm oracle_ruler`) 1.34 / 1.71 / 2.90 / 24.66 %. The error sits in highly expressed
multi-isoform genes whose TOTALS are right, and a ruler is judged by the WITHIN-GENE spread of its error
(`TRAPS: judge-a-ruler-by-its-within-gene-spread`). The junction-probed test-chromosome twin was re-simulated on the
current physics on 2026-10-03. THE ONE SHARED RULE'S RESIDUAL (2026-09-23, 60 probed transcripts against the
simulator): exact (1.000) on the objects where a transcript's fragments are gDNA's — contained pieces and cuts away
from junctions, exon edges and ends, 41 % of a probed transcript's yield — and its junction cuts, 40 % of a probed
transcript's yield, are what the junction price covers. Its cuts within a fragment of an exon edge are captured 1.10×
more than gDNA prices them (19 % of the yield), and the simulator binds a fragment through its best single probe part,
which captures a little less than the sum when both exons carry separate probes: that is the ~2 % the classes stay
apart with the true gDNA on every object (`ISSUES: the-gdna-component-length-rule-differs-from-the-transcripts`).
THE DECISION IT WAITS ON: the probe design is the one observable that sees isoform-specific capture directly, and
Rigel reads no panel (`DESIGN.md` §7.2). Real panels often span junctions and the design file is usually
unavailable (owner, 2026-09-19), so the candidates are data-derived: the junction price from gDNA (shipped, its
precision open); a capture field fitted from the coverage shape around each probe and the gDNA footprint; a
per-kit capture profile learned across a cohort; and a sparsity prior on isoform support as the safety net. The
spliced reads at each junction are refused (`ISSUES: the-spliced-read-junction-price`). `quant_accuracy.py --arm
oracle_ruler`, `ruler_vs_truth.py`.

### overlapping-synthetic-shadows
`priority: parked — Tier 5, an owner decision on the index; not a 0.8.0 number · kind: design · 2026-09-20`
`create_nrna_transcripts` clusters TSS/TES within a tolerance per strand and manufactures one synthetic
nascent entity per merged span, so a gene with several isoform clusters carries several near-coincident
shadows: 6,366 of the ladder index's 6,919 overlap a same-strand twin, 5,879 at ≥ 90 %, degrees up to 240.
The twins are indistinguishable to the EM, which puts a family's mass on an arbitrary holder (the live
entity is the holder in 78 % of live families off capture, 42 % on), so every PER-ENTITY nascent number is
a coin flip and the family is the unit that can be scored. They do not cause the capture-ON siphon —
collapsing every family to one component in the shipped EM reads 843k against 692k nascent at `g50 ss.99
ON` — and they cost 7 % of live nascent off capture (live families 942,567 against 1,013,400 at
`g50 ss.99 OFF`). The pool rows and the transcript table do not see them (synthetic rows are dropped);
per-entity diagnostics do. The decision is whether the index should merge them
(one span per gene per strand) or the scoring should. `quant_accuracy.py`, the per-family scoring in
`~/Downloads/rigel_runs/prototypes/2026-09-20_siphon_mechanism/`.

### yield-variance-beside-the-count
`priority: parked — Tier 5, with the per-transcript prior lane; the release ships the contraction as it stands (owner, 2026-09-17) · kind: build · 2026-09-17`
The capture-contracted yield is a posterior mean and carries no variance, so a count on a small yield reads as a
large abundance with the same error bar as any other. Measured on VCaP on the per-base rule by drawing every piece
efficiency from its posterior on the landscape grid and re-running the EM (eight draws, the EM's seed fixed; the
session's `yield_draws.py`, scored by `draws_census.py`): the counts themselves are stable — the dominant isoform
under the draws' mean equals the point estimate's in 98.5 % of multi-isoform genes with ≥ 20 fragments and is
identical in every draw for 90.6 % — while the yield's uncertainty is a CV of 0.006 / 0.061 / 0.317 at the 10th /
50th / 90th percentile of transcripts with ≥ 20 fragments, above the Poisson CV (median 0.100) for 36.4 % of them.
So the variance need not be propagated through the EM; it is a per-transcript number to publish beside the count:
the yield's posterior sd from its objects' posterior variances under the landscape, so that a user's abundance
carries `CV² = 1/k + Var[eff_t]/eff_t²`. On the one shared rule the form runs over the objects — a piece at its
contained share, a cut at its conserved share, a junction's price a signed sum of four objects' posteriors
(`EQUATIONS.md` §11); the objects' posteriors are what `capture_efficiency` already computes, and the second moment
is one more `w @ clipped²`. Derive → gate against the draws' spread on VCaP → `src/`: a result field per object and
an `em_effective_length_sd` column. Its consumer is the allocation of the pooled RNA prior across a locus's
transcripts (`ISSUES: per-transcript-prior-lane`): today the transcripts duel for the ambiguous fragments with no
per-transcript prior, and a variance per transcript is what an allocation better than equal shares would read
(owner, 2026-09-17). Instrument: the draws' spread is the truth the analytic form is gated against.

### capture-premise-untested-on-cdna
`priority: parked — Tier 5, watch: no library in hand can test it · kind: risk · 2026-09-17`
The contracted length reads the panel from gDNA and applies the same efficiencies to cDNA: the premise that
off-target cDNA is depleted as off-target gDNA is, which the simulator satisfies by construction and which
cross-hybridisation of cDNA to paralog probes or nonspecific binding could make milder by an order of magnitude.
Its exposure is not the scale (TPM is on the plain length by ruling and counts see only ratios within a locus) but
the isoform split: on VCaP 4,677 of 12,039 multi-isoform genes with ≥ 20 fragments change dominant isoform between
the contracted and the plain yield on the per-base rule (`DESIGN.md` §7.2, the yield's two consumers), two thirds
of the isoform-level fragment mass with them, stable under the yield's posterior — so the flips are the model's
answer, right if the premise holds. The one experiment that settles it: a sample sequenced both ways, panel and
whole transcriptome — the count ratio on unprobed transcripts is the cDNA depletion directly, and the capture
efficiencies give the gDNA depletion on the same objects. The lab has no such pair (owner, 2026-09-17); the
session's `lever_census.py` is the read-only smoke test to re-run on any new captured library, and `count_unambig`
beside `count` is what tells a user which isoform assignments rest on shared fragments alone.

### the-atom-at-an-unwitnessed-both-strand-slot
`priority: parked — Tier 5, accepted as a limit of the information (owner, 2026-09-14) · kind: known limit · 2026-09-14`
The tilt atom (`DESIGN.md` §6b.15.13, `EQUATIONS.md` §9f) admits "all of this slot's RNA is on strand s" unless
a held RNA level on the other strand rules it out. On stranded data at a both-strand slot whose minor strand
has no witness — no junction certifies it and no single-strand piece of its own exists anywhere, the case of
a single-exon gene wholly inside another gene's exon — the strand split cannot tell "pure s with gDNA" from
"both strands without": the pure hypothesis fits the split with no parameter and the atom pushes toward it.
Measured on the golden toy `antisense_contained` (1,000 fragments): the exon∩exon slot's truth is 0 gDNA;
prior-free the continuum reads 0.244 and the atom 0.483; with the toy's landscape of 8 training regions
0.005 and 0.242 — the antisense transcript's count 81 → 0 and 177.6 false gDNA fragments in a gDNA-free locus.
The mono shared-exon stress (no junction) doubles its false gDNA at 2–20 % minor. On the ladder, where the
landscape trains on ~15k regions, the same slots cost a few fragments each (the (0.2, 0.5] band, `g05 ss.99 ON`
4,942 → 7,005, `g00 ss.99 OFF` 35 → 73) against the strand-pure band's −11k / −13k. THE STANCE: this is not
enough information, and it is accepted as such — no presence witness (a locus with RNA elsewhere does not
imply this slot is expressed; a per-transcript active set would add bookkeeping and no measurement, since the
failing structure has no node where the transcript can be measured alone). On a real library the landscape
prior is the deciding voice; its location floor (`DESIGN.md` §7.1 rule 4, 2026-09-14) is the lever that was
pulled, and the stranded zero controls read 265 / 172 with it — while the 1,000-fragment golden, whose prior
now fits from about four anchors, reads 200.8 false gDNA (from 177.6): the toy's limit. On unstranded data the
atom is inert
(κ = ½ makes the three hypotheses equal) and such a slot is the prior's entirely. Watch it on the census
(`ISSUES: the-tilt-census-as-an-instrument`) and the mono stress; a real library that shows it would be a
nested single-exon antisense gene on a shallow stranded run.

### two-sided-exon-row
`priority: parked — Tier 5, deferred by the owner (2026-09-13), with the enrichment witness · kind: problem · 2026-09-04`
On unstranded data an exon's held row is the intron's composition through the face map, whose upper side is
the map's plateau above its ceiling — a lower bound on gDNA — so at pass zero an unstranded licensed exon
reads ~9× its true gDNA (+2,000 % at `g25 ss.50 OFF`) and forwarding it compounds the bias. A toy gate
(retired 2026-09-28 with its harness) read a factor of 88 (exon |Δf_g| 0.762 beside a pure-gDNA intron, 0.0086
beside a nascent-bearing one); on the ladder the class is 4–6 % of the unstranded error with
`transfer` at parity there, on the test chromosome 20–32 % and transfer worse (2,647 → 4,056 at
`g50 ss.50 OFF`). Every cap, two-sided row or wall is refused where junction probes enrich the flux more than
the crossing (`ISSUES: the-certified-flux-row-as-a-level`, `ISSUES: two-sided-exon-row-forms`,
`ISSUES: the-discrepancy-priced-cap`, `ISSUES: the-two-sided-level-lane`,
`ISSUES: the-wall-above-the-face-map-ceiling`). What would license a two-sided level is whether this library's
gDNA is enriched, witnessed on unstranded data by silent genes' exons against the intergenic density — the
landscape prior's job. The first step when taken up is that witness's derivation; judge at `R exon
(licensed)` and the walled classes, halves apart. `policy_benchmark.py --by-class`.

### the-lower-bound-noise-ratchet
`priority: parked — Tier 5, with the enrichment witness · kind: defect · 2026-09-05`
A level from an RNA-rich node's own strand profile has a mode that is noise around zero; its lower side bounds
its neighbours and the tightest noisy neighbour wins (ladder `g05 ss.99 OFF` 44,714 → 45,076; test chromosome
`g00 ss.99 ON` 64 → 216, 187 of it at `capcluster_ab`'s inner termini). A two-sided own-profile level keeps
fewer fragments overall, so lower-only stays. The gDNA lane's EDGE level does it too (2026-09-14, on a toy
of the encompassing locus): a shallow single-strand flank (404 fragments, truth 0.530, local solve 0.546)
reads 0.596 under the intergenic neighbour's Poisson level — a lower bound at that neighbour's sampled
density, 0.298/bp against the flank's realised 0.27, a 1.6σ excursion the hop's price blurs but does not
move. The cure is the enrichment witness `ISSUES: two-sided-exon-row` waits for. Instrument: the
ladder's `g00` rows (`calibration_vs_oracle.py`).

### flux-floor-dispersion
`priority: parked — Tier 5, with the transport-dispersion decomposition · kind: question · 2026-09-08`
The certified flux at a junction is a lower-sided estimate of the exon's strand RNA level priced by the pair
(`DESIGN.md` §6b.13); the route rate scatters beyond counting (median −3 %, 5–9 % at depth; 0–40 % over on
nine block readings) and at a pure-RNA exon the price cannot see it, so a lucky over-read is a sharp floor a
few points too high. Its decomposition instrument, `transport_dispersion.py`, was retired 2026-09-14 (in git).

### splice-out-premise-bias-uncorrected
`priority: parked — Tier 5, the owner's call · kind: decision · 2026-09-02`
The splice-out message assumes spliced and unspliced fragments at one face share capture affinity; measured,
the premise fails as a bias under capture (`log a` ≈ 0 off, +0.28 under exon probes, +0.78 under junction
probes), so only a subtracted level — a cross-locale fudge, the owner's call — would correct it. No instrument
yet.

### the-tilt-census-as-an-instrument
`priority: parked — Tier 5 · kind: build · 2026-09-13`
"Where does the strand tilt matter, and how does the tool do there" is answered only by a scratch script: per
stratum, the AMBIG slots by the RNA on each strand (bands of the minor strand's share), their depth, the gDNA
error, the TILT read-out's error (`|Δτ|·R/2`, RNA fragments on the wrong strand — measured by nothing else) and the
predicted θ peak width. `policy_benchmark.py --by-class` ranks by node class and not by strand split or tilt
error, so the tilt-error column is the new part; a census of what it overlaps comes before it becomes an
instrument (the instrument ruling: few instruments, kept current).

### transfer-variance-premise
`priority: parked — Tier 5 · kind: question · 2026-08`
Does a hop's transfer variance price a ratio built on a handful of counts? The policy prices every hop by both
witnesses' counting plus the pair's disagreement (``hop_price``, `native/transfer_rows.h`); whether that is right where a
pair agrees by coincidence is the open half (the landscape is no substitute: ~10× over-stated). `EQUATIONS.md`
§3.5h.

### drain-contaminates-certified-rna
`priority: parked — Tier 5 (owner, 2026-09-01) · kind: defect · 2026-08-31`
The second pass deposits some true-gDNA fragments into the certified-RNA banks: 233 records at
`g50 ss.99 OFF`, 1,482 at `g98 ss.99 ON` (1.9 % of that in-scope channel). The leak is exactly posterior
sampling and a no-leak counterfactual moves the 0.8.0 metric −0.60 %/−0.05 % at worst, so the harm is the
certainty claim and no in-solve correction is licensed (`ISSUES: drain-provenance-split`). To build: `DrainQC`
records `Σ q_null`; the middle-bin posterior bias repaired as its own A/B. `calibration_vs_oracle.py`.

### crossing-pool-contrast
`priority: parked — Tier 5, blocked on ISSUES: capture-blind-gdna-divisor · kind: question · 2026-08-31`
A second gDNA length contrast on the crossing pools: with oracle weights it beats the contained one under
capture (TV 0.076 vs 0.136; 0.078 vs 0.182) and starves off capture (pool 3 at 29–630 fragments). Blocked: the
weight estimator does not transfer under capture (`ISSUES: capture-blind-gdna-divisor`), and no shadow
transcript overlaps a gene edge, so pool 3 reads exactly 1.0000 pure
(`TRAPS: purity-is-a-property-of-the-annotation`).

### capture-degeneracy-standing-risk
`priority: parked — Tier 5, watch · kind: risk · 2026-08-31`
The gDNA two-pool contrast survives capture by a degeneracy: the shared-contaminant assumption is false under
capture (TV 0.95 vs 0.06–0.14 off) and it is safe only because the intergenic pool is depleted-not-impure, so
`a_0` clips to 1 and the algebra collapses to `g = f_0`. A probe panel that put RNA into intergenic space would
break it silently; `_deconvolved_gdna_counts` carries the derivation. No panel can fire it today. fl's one-sided
off-target rate absorbs evenly spread intergenic RNA, the same break in another estimator
(`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value`).

### pure-rna-mirror-asymmetry
`priority: parked — Tier 5 · kind: defect · 2026-08`
Two exact per-fragment mirrors of a pure-RNA library deconvolve differently in `count_gdna_region` by a few
percent, neither boundary-only nor monotone in strandedness. An R1-sense library is simulable; no instrument
yet.

### binary-cuts-on-continuous-quantities
`priority: parked — Tier 5, tracking only (owner, 2026-09-28: one entry; the inventory stays outside the tree) · kind: tracking · 2026-09-28`
The owner's rule is no yes/no gates on continuous quantities; the 2026-09-24 inventory ranks 55 such cuts by harm.
Homed elsewhere: whole counts and the fixed seed; the closed ranks; the od ceiling, deleted with the od fits
(`ISSUES: strand-overdispersion-one-shared-value`); the splice and multimapper cuts. Untracked, with measured harm:
- the reference ρ_ref is the landscape's argmax cell: one cell moves 91–100 % of transcripts' EM length by more than
  1 %;
- `calib_refit_iters = 3` stops unconverged: an L1 change of 85–34k fragments at a fourth refit;
- kernels with count < 1 get no location (the low-depth regime of `ISSUES: the-gdna-landscape-collapses-at-low-depth`);
- pruning on a k·ln10 lattice;
- the zero-density veto, which hits 27–43 % of held records;
- one unit merges two loci;
- `TAIL_DECAY`;
- `max_frag_length = 1000`, in three roles.
Each is raised when its mechanism is next touched, never as a sweep.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-24_cut_inventory/INVENTORY.md`.

---

## CLOSED / REFUSED — do not rebuild these; append-only

Every entry keeps its stamped measurement exactly as recorded: a graveyard row without its number is an
invitation to rebuild. A row measured on "all 36 conditions" or quoting `g01`/`g10`/`g25`/`g75`/`g90` predates
the ladder retired 2026-08-13 — the verdict stands as a record, and re-opening one means re-running it on the
current panel. Where a mechanism's only target was unstranded × capture-ON the row is moot as a 0.8.0
candidate on top of being refused; the `g00` zero-control column is never moot.

### flgap-panels-stale-nascent-model
CLOSED 2026-10-03 by the re-simulation. The two fl-gap side panels had not been regenerated in the sparse-nascent rebuild
and carried the retired uniform `fragment_share` nascent model (and, with the junction-probed test-chromosome twin, the
capture physics before 600343bd), so a claim spanning the ladder and a side panel varied two things. RULED (owner,
2026-09-28): re-simulate, then delete the mode. DONE: the three panels were deleted and re-simulated on the current
simulator with the ladder's `nrna:` block verbatim and `emit_fastq: false` (`flgap_rna_long` 2.6 GB, `flgap_rna_short`
2.9 GB, both 4/4 certified COMPOSITION + FIELD; `scenarios_probes_junction` 30 conditions, the same two COMPOSITION-only
OFF rows as the main test panel), their realised lengths measured off `truth_fragment_lengths.tsv` (RNA-long gDNA 78.58 /
RNA 247.63 bp, RNA-short 249.59 / 78.43 bp, capture OFF) and every baseline re-recorded pinned
(`ISSUES: calibration-detects-capture-on-a-capture-off-library`, `ISSUES: the-scorer-reads-a-census-length-law`,
`docs/TESTING.md`). The simulator's `fragment_share` mode, its `shares` field, its two tests and the inert `nrna:` blocks
of the ten test-chromosome configs are deleted; `additive_ratio` and `sparse` remain. Every number measured on the old
panels is void (the fl-gap rows of `ISSUES: eb-shrinkage-magic-ess` say so).

### an-unconditional-class-conditional-landscape
REFUSED 2026-10-03 (outside the tree, pinned and fractional, four arms: shipped, Part 1, the split, both). ψ's gDNA
landscape fitted per class with the SHIPPED estimator — exon regions and boundaries touching an exon on the exon-class
training set, every other slot on the rest, two sweeps on one lattice merged by class (exact: the messages never read
the prior), the capture reference and the efficiencies read off the shipped solve so only ψ's prior moved. It finds the
prize Part 1's verdict named — test chromosome stranded × ON object error −15 % (pool 16.5k → 8.4k), the suite's deferred
rows 18.8 → 14.6 and 8.2 → 6.7 % — and fails in scope: test chromosome capture-OFF object error ×1.06 stranded, ×1.73
at ss 0.70, ×3.72 unstranded (genes 1.33 → 2.87 %), where a class's training set is thin and the exon class has, by
construction, no zero-count anchor; the ladder's g00 stranded × OFF zero control reads 15,043 false gDNA fragments
against 127 (taken from nascent RNA); the RNA-short false-reference row 42.3 → 45.5 %. Off capture the two classes
share one true density, so the split only adds noise there; it helps exactly where the class has its own evidence.
Also measured: a class-wise belief fed to the pooled landscape fit makes it bimodal and the shipped reader invents a
capture reference on a capture-OFF unstranded row (0.165/bp from 533 members). A split CONDITIONAL on evidence that the
classes differ is a different mechanism and is not refused here (`ISSUES: the-gdna-prior-enters-psi-twice`).
Instrument: `~/Downloads/rigel_runs/prototypes/2026-10-03_fl_arms/` (`RIGEL_ARM=classwise`; `VERDICT.md`, `q4/`).

### a-length-table-built-from-the-other-origins-law
REFUSED 2026-09-30 at phase 0 (oracle arms outside the tree, pinned, fractional, 198 libraries). Two constructions of one
origin's scorer table from the other's law times the true length ratio `r = f_R / f_G`: RNA := normalize(g · r), an
external review's proposal, and its mirror gDNA := normalize(rna_pmf / r). Each multiplies a law with little
independent evidence in its tail by a large ratio, so the tail it makes is the smoothing floor's. On the test
chromosome's stranded × capture-ON rows (gDNA as EM + intergenic): RNA from gDNA on the RNA-long panel (RNA 250 ± 150
against gDNA 100 ± 50) reads transcripts ×6.4, genes ×70, gDNA pool ×92, its RNA table's mean 493 bp; gDNA from RNA on
the gDNA-long panel reads transcripts ×2.1, genes ×8.7, gDNA pool ×79. It does not show that two independently
estimated laws cannot hold different supports: a design keeps each origin's own evidence on the union of their supports
(`ISSUES: the-scorer-reads-a-census-length-law`).
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-30_frame_phase0/` (arms R_r and N_r).

### the-realized-gdna-length-law-reads-rna-counts
CLOSED 2026-09-30 by the fix, in the working tree for the owner's commit (`calibration/fl.py`). `_realized_gdna_counts`
normalises the RNA law it is handed, and `build_fl_models` hands it `FLModels.rna_pmf`, the EB law the scorer reads;
the boundary's expected count E_b = ρ_off(μ_g − 1) + ρ_adj(μ_r − 1) retires the `max(μ_r − 1, 1e-9)` and
`max(rna.sum(), 1e-30)` guards. Gates in `test_fl.py`, each failing before: the census is the same at any scale of the
RNA law; off capture, RNA crossing a boundary is not captured gDNA (on-target share 0, the realized law equal to the
uniform); the census reads `FLModels.rna_pmf`. The N_s = 0 fixture was re-derived.
MEASURED on the ladder's no-EM census, realized mean in bp, before → after (truth): g05 ss.99 OFF 219.95 → 216.99
(216.90), realized − uniform +0.40 (ss.50: +0.52); g50 OFF +0.04 / +0.05, its rectification floor; g05 / g50 / g98
ss.99 ON 236.22 → 244.96 / 234.52 → 239.11 / 234.87 → 238.58 (240.97 / 240.84 / 240.87). The entry's three
falsifiers hold. `calibration_vs_oracle.py` moves region Σ|Δ| by +10 / −15 / −2 on the in-scope strata (2026-09-28's
table, whose prototype the landed tree reproduces bit for bit).
THE COST, end to end (`quant_accuracy.py`, pinned and fractional, 198 libraries; the landed tree reproduces the A/B'd
arm on 96 of 96 rows checked): ladder stranded × OFF transcripts −701, genes −93; unstranded × OFF −593 / +42;
stranded × ON +15,085 (+1.75 %) / +6,427 (+4.20 %), pools gDNA +5.6k, nascent +4.6k, annotated +11.3k. The
test chromosome's stranded capture-ON rows about +2 % genes (base panel +2,068, the 200 bp control +2,154); the fl-gap
panels keep their length channel. Goldens moved in two scenarios: `combo_extreme` gDNA 209.5 → 192.6, nascent 259.9 →
277.1; `combo_moderate` −8.1 / +8.5; four others by ≤ 0.004. The bug was cancelling part of the census-vs-census frame
gap, which now owns the stranded × ON cost (`ISSUES: the-scorer-reads-a-census-length-law`). Its steps 3 and 4 moved
to `ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value` (e) and
`ISSUES: the-fl-boundary-inversion-has-underived-pieces` (e).

### the-scorer-reads-the-uniform-gdna-law
REFUSED 2026-09-30 at the A/B (pinned, fractional, 198 libraries). The arm: the realized law deleted with its
coupling, so the scorer and geometry read one gDNA law, the contrast's. It gives gDNA one law in its table and its
normaliser and deletes about 300 lines, but the scorer's RNA table is a captured census, so the uncaptured gDNA law
tips long unspliced gDNA into nascent RNA: ladder stranded × ON transcripts +22,790 (+2.7 %), genes −678, gDNA pool
+81.5k, nascent +69.8k (g50 ss.99 ON nascent +55.8k on a true 150k, +37 %); the test chromosome's stranded capture-ON
rows about +10 % genes (+9,988, +9,645, +9,864); fl-gap long ON genes +3,699 (+8 %). The scorer's half alone (the
shipped `gdna_pmf` read by the scorer, geometry unchanged): ladder stranded × ON −2,535 / −5,363, gDNA +63.7k, nascent
+53.1k. The geometry's half is a second cancelling pair, the coupling's pull on `gdna_pmf` against the capture ruler:
+15.3k transcripts at g05 ss.99 ON (predicted +16.2k). Not a candidate until RNA's scorer table is in the same frame
(`ISSUES: the-scorer-reads-a-census-length-law`).

### the-gdna-length-log-ratio-table
REFUSED 2026-09-29 at stage 0, a prototype outside the tree: six of nine pre-registered falsifiers fired. The design:
the EM's gDNA table as its RNA table times the gDNA–RNA length log-ratio δ, measured inside single-strand objects by
distinct-pair strand unmixing and a two-pool contrast, combined by inverse variance, faded, and held as a Gaussian
log-ratio outside RNA's range. What held: the strand source read the labelled within-object laws to about 2 bp, and
where δ = 0 the design is the no-length arm, worth −13.6k transcripts and −8.6k genes on the ladder's stranded × ON
rows. What killed it:
- THE TABLE. `rna_pmf·e^δ` places gDNA mass beyond RNA's support only through the RNA law's EB floor, which the held
  edge value multiplies. On the test chromosome's gDNA-long panel (gDNA 250 ± 150 against RNA 100 ± 50) stranded × OFF
  genes read shipped 27.5k, no length 32.2k, design 46.1k; the table's mean 433 bp against a true 261, with 64 % of its
  mass beyond 436 bp (true 7 %). Oracle moments still read perr .215 against shipped's .012, and the labelled per-length
  ratio with no family .177: the construction fails, not the estimator.
- POOLING. Pooling normalised laws across objects before dividing bends a real ratio's shape: an exact two-object
  example reads 0.532 nats wrong (an external review, reproduced). No equal-law gate can see it.
- THE CONTRAST AS A δ SOURCE. On equal laws it reads small systematic gaps that its covariance calls significant (the
  test chromosome's base panel: −1.6 to −27 bp on 14 of the 15 rows with λ > 0; unstranded × ON g50 +7,766 genes
  against no length), and under capture it is biased toward zero (−107 against −162 bp on fl-gap long ON).
- THE FAMILY. Two moments per law cost 1.5–3.5× perr against the labelled per-length ratio.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-29_fl_tilt_stage0/` (`report.md`, `REPORT_FOR_REVIEW.md`).

### the-total-density-landscape
CLOSED 2026-09-28 (owner): retired as QC-only, with `total_abundance.py`, `diagnostics.py`, its two wall inputs
(`build_mature_wall_distances`, `build_contiguous_boundary_reach_arrays`), the two feathers it wrote
(`gdna_density_kde.feather`, `gdna_density_regions.feather`) and what only it read (`signature.mrna_active_strands`, the
substrate's start/end/span copies, `fit_landscape`'s `knn_scale`, `split_basins`' anchor rate). No solve, report or
instrument read it: all 29 `CalibrationResult` fields were bit-identical with and without its wall inputs on six
test-chromosome conditions, and `rename_identity.py --check` read bit-identical on the three references frozen at
509248f3. Its mode census (`split_basins`, `located_enriched_mode`, the ruler's reference) moved to `landscape.py`; the
START bank's model-free identity and its exactness condition stay in `EQUATIONS.md` §2.3b. Its rulings as they stood
(2026-08-21):
- THE WALL RULE: a START (END) side is exact iff the template continues `w_max − 1` bases past the region's
  genomic-high (low) bound, with `w_max` read from the support end of `deposited_lengths`, never a quantile; the
  distance is the component minimum over `T(slot)`; a double-walled region is not model-free. Coverage, START-mass
  weighted on the ladder: 94.7 % at capture-OFF, 84.3 % under capture.
- WHAT A CONSUMER COULD READ (a grid sweep over 16 conditions, `CONSTANTS.landscape.grid_points` over a 16× range): `rho_0` (it moved
  8–25 %) and the anchor verdict (12/12 contaminated rows) were consumable; `span_R` was not (58 → 77 → 95.6 → 94.7
  → 1.9 on `g50 ss0.99` OFF as the grid refined); the mode count never
  (`TRAPS: a-mode-count-is-not-a-well-posed-quantity`). A fit on `mass / eff_gdna` carries the divisor's per-region
  spread (offset IQR 0.12 nats off capture, 1.66 under it), which no bandwidth removes.

### benchmark-noise-floors-unmeasured
CLOSED 2026-09-28 (owner): not a project. A last-bit change moves the genes by at most 15 fragments and the gDNA pool
by 113, while the per-transcript Σ|Δ| amplifies it through the EM's start-dependent isoform split (thousands of
fragments: 7,968 on a 4-thread rerun at `g25 ss.50 ON`; `ISSUES: the-em-answer-depends-on-where-it-starts`).
THE RULE: an A/B pair runs with the scan pinned (`--set scan.total_threads=1`), so it is exactly reproducible, and an
effect is judged by its size, with genes and pools read beside the transcript table.

### calibration-vs-oracle-cannot-see-a-ruler-move
CLOSED 2026-09-28: the output that read 0 by construction is deleted — `calibration_vs_oracle.py`'s ruler section and
`prior_vs_oracle.py`'s `gdna_eff_len` score, since the `O` arm keeps `P`'s efficiencies, reference density and gDNA
region lengths, everything the length reads. The length's truth is `ruler_vs_truth.py`; the true-efficiency `O` arm is
`ISSUES: instrument-ledger` (g).

### rna-prior-floor-at-pure-gdna-loci
CLOSED 2026-09-26 (owner) for 0.8.0, to be reopened only with a fundamentally different algorithm. The prize is real —
a perfect calibration prior takes `g98 ss.99 ON` from 25.26 to 19.99 % of the transcripts' true count — but neither
mechanism built for it improves an in-scope stratum without costing another: a no-RNA state read out beside the solve
(`ISSUES: a-no-rna-state-read-beside-the-solve`, refused as a design) and end states inside the solve
(`ISSUES: no-rna-and-no-gdna-end-states-in-psi`, refused at the A/B). What a new approach has to answer, from both:
the weight of "no RNA" at an object must come from evidence whose precision grows with depth, and on unstranded data
an object's strand term carries none (`DESIGN.md` §6b.1); it must not pull objects holding 1–50 % RNA to "no RNA",
where both overshot; and what it changes reaches every contracted length through the capture efficiencies, so it is
judged on the capture-ON rows as well as the pools. The entry as it stood
(`priority: NEXT — the owner's residual pass (2026-09-24): the largest owner of g98's gDNA error · kind: defect · 2026-09-20`):
DISSECTED AGAIN 2026-09-24 on the count form, against per-fragment truth (four lenses, a synthesis and two refuters;
`~/Downloads/rigel_runs/prototypes/2026-09-24_g98_dissection/`). The read-out is the posterior MEDIAN of λ
(`psi_kernel.h`, `DESIGN.md` §6c), not a mean: under the Jeffreys reference with no atom at zero RNA, an object whose
true RNA is below its resolution `1/√(n·I)` reads about 0.42·√(n/I) fragments of RNA, always toward RNA (measured
0.21–0.28·√n stranded, 0.32–0.37·√n unstranded; the stranded/unstranded ratio 0.721 against the Fisher-information
ratio 0.714). It is per OBJECT, and 78–85 % of calibration's count error sits on objects with no RNA, mostly INSIDE
expressed loci; pure-gDNA loci hold only −3.9k / −5.5k / −9.7k. The EM passes the count through at g98 (per-locus
slope 0.90–1.00), and off capture it lands in synthetic spans at introns. Size (the truth restored on no-RNA objects
alone): gDNA −36.7k → −4.7k at `g98 ss.99 OFF`, −61.6k → −12.3k at `g98 ss.50 OFF`; about 47k at `g98 ss.99 ON` once
`ISSUES: the-gdna-length-law-falls-back-at-identical-purities` is separated out. Its transcript cost is small (1.3k
of 1.8k at ss.99 OFF, none at ss.50 OFF). The repair: an atom at zero RNA whose weight is the population's own
no-RNA share, fitted by marginal likelihood across objects (no constant), and a CONTINUOUS read-out — a median of a
posterior carrying an atom is a yes/no cut — which amends the §6c median ruling: the owner's decision. A test bar
must be stated per `n`: at an atom weight of 0.75 the floor falls about 74 % at n = 10 and 89 % at n = 1,000.
The history below measured the pre-count-form RNA pseudocount; its "posterior MEAN" is corrected above.
`rna_prior_count` over-states by +64.5 % at `g98 ss.99 ON` (180,806 against 109,915) and +32.4 % at
`g98 ss.99 OFF`, and it is a diffuse positive floor, not a few loci: 1,089 of 1,149 loci read high and 55k
of the 73k excess sits in loci that are over 99 % gDNA where the true RNA is 544 fragments over 867 loci.
DISSECTED per object (2026-09-20, `g98 ss.99 ON`): on the 14,280 regions holding gDNA and NO RNA the
calibration places 19,538 RNA fragments, on 97 % of them, and the floor is per OBJECT and grows roughly as
the square root of the object's gDNA — a median 0.14 fragments on regions of 1–10 gDNA (11.6 % of their
gDNA), 1.8 at 100–1,000 (1.2 %), 26 at 10,000+ (0.19 %) — the signature of a posterior MEAN of a composition
whose likelihood sits at the boundary `f_r = 0` with a width `~1/√n` (`TRAPS: zero-target-guards-are-one-sided`).
Most of the mass is on BOUNDARIES: 231,410 RNA incidences against 63,972 true, +85,112 on boundaries with
no RNA at all (1.6 % of their gDNA incidences) and +82,326 on mixed ones; regions carry +24,291. ⛔ NOT a
constant fraction and NOT a cheap fix: it is the composition solve's estimator at a pure-gDNA object — a
posterior mean cannot read zero — so the repair is an atom at `f_r = 0` in ψ's hypothesis space (the tilt
already carries {pure +, pure −, mixed}; the composition does not) or a different estimand there, a solver
design item, and it lives at the stress rung (at `g50` the same floor is +4,309 on 519 near-pure loci,
0.15 % of the RNA, invisible in the pool). `prior_vs_oracle.py`, `calibration_vs_oracle.py`; the confidently-wrong
class's z-band table retired with `solvability_audit.py` (2026-09-28). It is why the per-transcript allocation made `g98`
capture-ON worse before the gDNA opportunity was corrected (`ISSUES: per-transcript-prior-lane`).

### a-no-rna-state-read-beside-the-solve
REFUSED 2026-09-26 by the owner, as a design: a second solver on top of the first. The fix belongs in the solve, which
should be able to say "no RNA" itself, and it must leave the tool simpler. Kept as the measurement of the prize of
`ISSUES: rna-prior-floor-at-pure-gdna-loci` (closed). The mechanism
(`~/Downloads/rigel_runs/prototypes/2026-09-24_zero_rna_atom/`, measured through the EM in
`~/Downloads/rigel_runs/prototypes/2026-09-26_zero_atom_em/`): a no-RNA state per object, its weight fitted per
annotation class by marginal likelihood, read as `P0 + (1 − P0)·median`, computed on the pipeline's own final sweep
and handed ONLY to the EM's prior gDNA count; calibration, its landscape and training, and every length as shipped.
All 16 rows, fractional, each against a fresh shipped run and a replicate that gives the row's own floor (0–59
transcript fragments). Transcripts % shipped → with it (perfect calibration prior): stranded OFF 1.988 → 1.979
(1.970), stranded ON 6.203 → 6.150 (6.149), unstranded OFF 2.010 → 2.017 (1.987); genes 0.207 → 0.202, 1.594 → 1.547,
0.238 → 0.249; the deferred stratum 10.28 → 9.80 (9.26). Per row: `g98 ss.99 ON` 25.26 → 20.14 % (19.99), genes
15.03 → 11.06; `g98 ss.99 OFF` 15.14 → 14.72 (14.25); `g50 ss.99 ON` 6.13 → 6.06 (6.07); `g50 ss.99 OFF` −1,262
fragments; every `g00` row's prior unchanged (the fitted weight is 0 there) and its table within 52 transcript
fragments, run-to-run noise that one replicate under-reads; `g05 ss.99 ON` +760 transcript fragments (genes −172);
`g98 ss.50 OFF` WORSE, 16.88 → 17.62 %, genes 5.66 → 6.95 %. The price: it only ever raises gDNA and it overshoots on
objects holding 1–50 % RNA (in fragments, `g50 ss.99 ON` 1–10 % band |error| 40.7k → 50.6k, 10–50 % 38.9k → 44.5k;
`g98 ss.99 ON` signed +4.5k → +26.1k and +1.1k → +3.8k), which lands on the synthetic spans (est ÷ true
`g98 ss.99 OFF` 1.88 → 0.36, `g98 ss.50 OFF` 2.52 → 0.27, `g50 ss.99 ON` 0.91 → 0.77). On unstranded data the strand
term cancels and the weight is read from the density terms alone, which is where it loses (capture OFF; the deferred
unstranded capture-ON stratum gains). Unmeasured: realistic nascent (on_fraction 0.10) at `g98`.

### no-rna-and-no-gdna-end-states-in-psi
REFUSED 2026-09-26 at the A/B; the owner closed the thread for 0.8.0. The owner's fix at the source: ψ's gDNA-fraction
axis given its two end states — "no RNA" (`f_g = 1`) and "no gDNA" (`f_g = 0`) — beside its continuum, each at the
continuum's reference mass, as the strand axis carries "all RNA on +" and "all RNA on −" (`DESIGN.md` §6b.15.13); in
the final solve of every sweep only; read out as each state's answer weighted by its posterior, the §6c median inside
the continuum. A C++ prototype in a worktree, switched off bit-identical to shipped on all 30 calibration arrays
through three refits at `g05 ss.50 OFF` and `g98 ss.99 ON`
(`~/Downloads/rigel_runs/prototypes/2026-09-26_presence_states/`: `DERIVATION.md`, `README.md`, and
`presence_states.patch`, the diff against `062bc9ea`).
Calibration against the oracle (region + boundary |gDNA error|, `g05`–`g98`): stranded OFF −14.8 %, stranded ON
−10.1 %, unstranded OFF −12.1 %, but `g05 ss.99 ON` +6.3 % (80.7k → 85.8k), `g50 ss.99 ON` +1.4 % and `g50 ss.50 OFF`
+1.2 % (the deferred `g05 ss.50 ON` +94.7 %); the `g00` controls stranded 353 → 80 (OFF) and 329 → 44 (ON), unstranded
366 → 552 (OFF) and 266 → 299 (the deferred ON); the "no RNA" state alone takes those two unstranded `g00` rows to
74,442 and 70,536.
Through the EM (fractional; transcripts / genes as a % of truth per stratum, shipped → this, the perfect calibration
prior in brackets): stranded OFF 1.988 → 1.975 (1.970) / 0.207 → 0.202; stranded ON 6.203 → 7.026 (6.149) / 1.594 →
1.778; unstranded OFF 2.010 → 2.013 (1.987) / 0.238 → 0.245; the deferred stratum 10.28 → 18.88. Per row:
`g98 ss.99 ON` 25.26 → 21.73 % (19.99), genes 15.03 → 11.53; `g05 ss.99 ON` 5.56 → 7.78 %, genes 0.59 → 1.16 (+204k
transcript fragments against a replicate floor of 25); the deferred `g05 ss.50 ON` 11.03 → 34.19 %.
Why it fails. (1) The transcript table's harm runs through the capture efficiencies, and so through every contracted
length: with the efficiencies and the capture reference held at the switched-off values, `g05 ss.99 ON`'s +204k
transcript fragments fall to −13 (floor 25), the deferred `g05 ss.50 ON`'s +2.13M to −230, and `g98 ss.99 ON` reads
21.42 %; the pools do not all return — at the deferred `g05 ss.50 ON` the counts alone move gDNA into the synthetic
spans (gDNA 391k → 251k, truth 500k; spans 380k → 520k, truth 285k). At `g05 ss.99 ON` 1,635 objects' efficiencies
move by more than 0.05 (up to 0.75), and objects holding 10–50 % RNA (median count 33) take P(no RNA) > ½ 31 % of the
time. (2) Three states at equal weight give each end even prior odds against the whole continuum, above the no-RNA
share the record fitted by marginal likelihood (`ISSUES: a-no-rna-state-read-beside-the-solve`) on every class at
`g00` and `g05` and at `g50` OFF, and on all but the intron classes at `g50 ss.99 ON`: 0 at `g00`, ≤ 0.024 at `g05`
OFF, 0–0.497 at `g05 ss.99 ON`, ≤ 0.43 at `g50` OFF. (3) On unstranded data the strand term is ½ in every state, so
the states are decided by terms whose precision does not grow with depth — the landscape's density and one-sided level
rows — and the fixed weight is the answer at every depth, the pattern `DESIGN.md` §6b.1 refuses. (4) A one-sided row
read at its held end favours the end state it does not bound, and the "no gDNA" state reads the landscape's density at
its floor cell as the mass of a point (`TRAPS: a-ratio-cannot-carry-zero`). (5) The published `f_g` mixes the states
while its `Var(log f_g)` stays the continuum's, so the landscape trains on a pair that no longer matches. Not built:
the weight learned per class in the refit loop — the fitted shares are near 0 where the states help least and reach
~0.5 on the intron boundaries at `g05 ss.99 ON`, where the harm is.

### the-gdna-length-law-falls-back-at-identical-purities
CLOSED 2026-09-26 by the fix, in the working tree for the owner's commit (`calibration/fl.py`,
`_deconvolved_gdna_counts`): at a purity separation of exactly 0 the contained pair answers with its own de-tilted
mixture — the limit its resolution-weighted fade approaches, and what the boundary pair already returned at its own
tie — instead of declining to the capture-selected four-pool census. MEASURED 2026-09-26 (no EM,
`~/Downloads/rigel_runs/prototypes/2026-09-26_gdna_law_tie/`): only `g98 ss.99 ON` ties on the ladder (both purities
1.000); there the uniform-frame law reads 245.1 → 217.2 bp (truth 216.7) and the realized law 254.3 → 234.9 (truth
240.9, where every other capture-ON row reads 234.5–236.5 against ~241: the census's own offset,
`ISSUES: the-scorer-reads-a-census-length-law`); the other 15 rows' laws are bit-identical, so nothing else moves.
Through the EM (fractional, `g98 ss.99 ON`), transcripts / genes / gDNA est − true / synthetic spans est (truth 6.0k)
/ annotated est − true: 26.18 → 25.26 % / 17.60 → 15.03 % / −281.0k → −246.5k / 79.7k → 53.7k / +26.9k → +18.2k; under
the oracle prior 21.34 → 20.04 % / 12.69 → 9.42 %; `--arm base_reseed` equals `base`. The landed pipeline reproduces
the prototype to 0.24 of 49,013 transcript fragments. The 2026-09-24 reading that transcripts worsened (25.95 →
26.46 %) does not reproduce on this tree. Gates (`tests/calibration/test_fl.py`), verified failing on the shipped
code: the tie at purity 1, the uniform law continuous through the tie end to end, and a tie below purity 1;
restoring the decline fires all three, deleting the tie branch fires the two `errstate` guards (bin `L = 0` is
empty in both pools, so `0/0` makes the inverted law's sum NaN and the old arithmetic reached the mixture by
accident). Left as it was: the contained pair still declines on an empty pool — the `g00` rows, which then read the
four-pool census (243–244 bp under capture) — the same kind of cut, a separate mechanism. The entry as it was opened:
`calibration/fl.py`'s `build_fl_models` estimates gDNA's fragment-length law by contrasting two pools of different
gDNA purity; when the separation is EXACTLY 0.0 ("purities identical", as at `g98 ss.99 ON`, both pools at gDNA
share 1.0) it declines and falls back to the capture-selected four-pool census: `gdna_pmf` reads 245.1 bp against a
true 216.7. Where the contrast is applied the estimator is right (216.7–217.0 at `g50 ss.99 ON` and `g98 ss.50 ON`,
separations 0.021 and −0.0018). A yes/no fallback on a continuous quantity. At `g98 ss.99 ON` this one root owns the
scorer's gDNA law reading +13 bp long (254.3 against 240.9), the gDNA component's offset from the RNA lengths'
scale (gDNA/spans 1.014, gDNA/isoforms 1.041, which fall to g50's levels with the true law), and 13.6k of
calibration's count deficit. Fixed alone: gDNA −100.5k → −67.7k, spans 80.0k → 54.1k, the oracle arm −35.5k →
−21.7k; but transcripts do not improve (25.95 → 26.46 %): correcting the law the contracted lengths use exposes an
opposite error in the length rule, so the repair is judged with its effect on the lengths. Class-mean L/Y moving
toward one scale is NOT a valid gate for it: that moved toward one scale while the EM got worse (−4.8k to −5.9k
gDNA, +6.2k to +9.8k transcripts). Scripts: `analysis/refute-mechanism/fl_root.py`, `scale_arm.py` in `~/Downloads/rigel_runs/prototypes/2026-09-24_g98_dissection/`.

### the-capture-length-owns-stranded-capture-on
CLOSED 2026-09-26 (owner), for now: the tool is at diminishing returns, a change must be simple and improve it, and
the shipped length stays. The frame as it was opened:
revisited with a new argument and a measurement · kind: defect · 2026-09-26`
Stranded × capture-ON is the worst in-scope stratum: transcript Σ|Δ| 6.21 % of truth against 1.99 % stranded OFF
and 2.01 % unstranded OFF (the ladder, fractional, arms of 2026-09-25). Handing the EM the simulator's own
capture-aware lengths for every transcript and synthetic span (`quant_accuracy.py --arm oracle_ruler`; the gDNA
component's length and the priors stay as shipped) takes it to 2.02 %, while a perfect calibration prior
(`--arm oracle`) leaves it at 6.16 %: the length the EM divides each component's reads by under capture
(`DESIGN.md` §7.2, `EQUATIONS.md` §11) owns about two-thirds of the stratum's transcript error. Per condition,
transcripts (genes) as a share of truth, shipped → true lengths: `g00` 6.47 → 1.42 % (2.36 → 0.19 %), `g05`
5.56 → 1.72 % (0.59 → 0.52 %), `g50` 6.13 → 2.84 % (1.43 → 1.37 %), `g98` 26.18 → 25.90 % (17.60 → 18.65 %). At
`g05` and `g50` the length decides the isoform split inside genes; at `g00`, with no gDNA to witness capture, it is
nearly the whole error and moves the split between annotated transcripts and synthetic spans as well; at `g98` it
is not the owner. The capture-OFF rows do not move, as they must not, and the deferred stratum moves the same way
(`g05 ss.50 ON` 11.03 → 1.46 %). A capability proof, never headroom: it hands over the true lengths.
Where it sits, as far as measured: at `g50 ss.99 ON`, 161k of the 168k within-gene error that lengths on their
yield remove sits in genes with a junction-probed isoform (`ISSUES: the-junction-price-is-noisy-within-a-gene`,
parked); the short exon pieces of unprobed transcripts read above their capture truth
(`ISSUES: the-efficiency-posterior-floor-on-empty-pieces`); the reference needs gDNA (`DESIGN.md` §7.2). Not
measured: the transcript cost of each class's length error alone. The per-class length error against the truth is
`ruler_vs_truth.py --scale`; the truth is the simulator's `CaptureSampler.partition_array`; the deliverable is
`quant_accuracy.py --arm base` beside `--arm oracle_ruler`, per stratum, under
`--set em.assignment_mode=fractional`.
WHAT THE CAMPAIGN MEASURED (2026-09-26, `~/Downloads/rigel_runs/prototypes/2026-09-26_capture_length/README.md`;
fractional, stranded × capture ON, transcripts / genes as a % of the true annotated RNA at `g00` / `g05` / `g50`):
* The frame reproduced exactly (base 6.47 / 5.56 / 6.13 / 26.18 with `g98`, `--arm oracle_ruler` 1.42 / 1.72 /
  2.84 / 25.90); under fractional assignment `--arm base_reseed` equals `base` to 0.01 point.
* `--arm oracle_ruler` is not a clean ceiling: it keeps the gDNA component's shipped length, so at `g50` the
  synthetic spans read 3.25× their truth and the gDNA pool −459k. With the gDNA component at its own true yield on
  the same scale the ceiling reads 1.44 / 1.45 / 2.08 (genes 0.19 / 0.23 / 0.58) with clean pools.
* Where it sits: at `g05` / `g50` about 80 % of what the true lengths remove is the isoform split inside probed
  multi-exon genes (within-gene error in genes with a junction-probed isoform 4.61 → 0.96 / 4.43 → 1.19 points);
  at `g00`, with no reference and nothing contracted, half of it is the spans against mRNA (gene-level 2.36 → 0.19).
  Probed genes without a junction probe, unprobed genes and the classes' mean offsets cost almost nothing, and the
  shipped spans already sit on the gDNA component's scale.
* It is the junction price's STRUCTURE, not its noise: with every other length true, the shipped sum on the
  simulator's noise-free gDNA weights alone reads 5.11 / 5.12 / 7.10 — the whole shipped error. On the simulator's
  own capture landscape (4,264 probed multi-exon transcripts) fragments across a junction carry 46 % of a probed
  transcript's yield, lie wholly in exons and are captured nearly uniformly (weight median 1,051 of a full probe's
  1,201; 1,092 where a probe spans the junction, 906 where none does; ~1.19× the reference), while a gDNA fragment
  across a junction's boundary lies half in the intron: noise-free the sum is off by a log sd of 0.64 per junction.
* The isoform split needs within-gene length ratios near 1 %: with every other length true, one true price per gene
  leaves 2.57 / 2.45 / 3.23, one price for every junction 3.23 / 2.98 / 4.05, two prices by whether a probe spans
  the junction 2.65 / 2.38 / 3.41.
* Swapping one class's true lengths into an otherwise true set overstates that class's cost — a gene's span and its
  isoforms are priced from the same noisy objects and their errors cancel (each class alone cost about the whole
  base error at `g50`) — so a class is also attributed from the shipped side, with every gene's own scale kept.
Refused on the way: `ISSUES: a-junction-price-from-its-genes-own-rna`,
`ISSUES: a-junction-price-at-the-best-object-its-fragments-reach`. What stays open beside it:
`ISSUES: the-efficiency-posterior-floor-on-empty-pieces`, `ISSUES: ruler-witness-geometry-on-transcript-panels`.

### the-junction-price-is-noisy-within-a-gene
CLOSED 2026-09-26 (owner), for now, with `ISSUES: the-capture-length-owns-stranded-capture-on`: the sum stays. That
campaign found the within-gene error to be the rule's structure — a fragment across a junction lies wholly in exons
and the gDNA beside it half in the intron — and refused the two prices tried against it. The record as it stood
(PARKED 2026-09-23; the pseudocount's odds fixed 2026-09-24):
SIZED 2026-09-24 (the g98 dissection's g50 contrast, `~/Downloads/rigel_runs/prototypes/2026-09-24_g98_dissection/`): with every length on its
yield at its locus's gDNA scale, the within-gene transcript error at `g50 ss.99 ON` falls 231.6k → 64.0k
(transcripts 6.14 → 1.74 %), and 161k of the 168k sits in genes with a junction-probed isoform; at `g98 ss.99 ON`
it is about 3.5k. Still parked: nothing derivable moves the noise-free ceiling.
The one shared length rule (`DESIGN.md` §7.2, `EQUATIONS.md` §11) prices a junction by conservation of bases,
`c_junction = c_lo + c_hi − ½(c_intron,lo + c_intron,hi)` (`capture_eff_length._cut_efficiencies`). It is right on
average — it is what puts the transcripts on the gDNA component's scale — but it widens the within-gene spread of
the transcripts' lengths, and the EM splits isoforms by that spread. The within-gene sd of log(L / Y) about each
gene's mean (annotated multi-exon transcripts with ≥ 20 fragments, no EM) beside the class mean L / Y ×10⁻³,
`g05` / `g50 ss.99 ON`:

| the junction priced by | within-gene sd | class mean |
|---|---|---|
| the per-base rule, before the port (`ISSUES: the-per-base-rule-for-every-component`) | 0.052 / 0.040 | 0.944 / 0.957 |
| its adjacent pieces (`ISSUES: the-junction-price-at-its-neighbours-mean`) | 0.057 / 0.041 | 0.936 / 0.965 |
| the sum, SHIPPED | 0.072 / 0.057 | 1.076 / 1.112 |

On the shipped tree (`ruler_vs_truth.py --scale`, width step 1) the gDNA component and the synthetic spans read
1.108 / 1.118 at `g05` and 1.131 / 1.141 at `g50`, 3.9 % / 2.6 % from one scale with the isoforms. Giving back the
scale is not free: the synthetic pool's elasticity to its own length is −21 at `g50 ON` (2026-09-21).
MEASURED 2026-09-23 — THE SPREAD IS THE RULE'S, NOT THE COUNTS'. The same rule fed four count sources, within-gene
sd (the shipped row is `--scale`'s, reproduced exactly):

| every region's and boundary's gDNA count | `g05` | `g50` |
|---|---|---|
| the simulator's expected count, as a plug-in — no noise, no prior | 0.049 | 0.053 |
| the same counts through the posterior | 0.054 | 0.053 |
| the true realized counts (`slot_truth.npz`) | 0.059 | 0.054 |
| calibration's, SHIPPED | 0.072 | 0.057 |

Swapping one term at a time (on the transcripts' own routed shares, where the shipped rule reads 0.074 / 0.058): the
junction at the transcript's own truth 0.048 / 0.030; all four junction terms noise-free 0.065 / 0.056 — the ceiling
of any price that only removes noise; the intron terms alone noise-free 0.075 / 0.060, no gain. Per junction the
noise-free rule is 0.27 / 0.31 rms from the truth on prices near 1, and the noise adds 0.23 / 0.13 in quadrature.
The rule's own error is the capture physics, which no gDNA object sees: the simulator binds a fragment through its best
single probe part, so a junction a probe spans (one part in the transcript, two in gDNA) is priced exactly by the
sum, while a junction with separate probes on its two sides binds one of them and the sum prices both — and the two
read the same to every region and boundary (`ISSUES: ruler-witness-geometry-on-transcript-panels`). Where both sides
carry capture the sum over-prices by 18 % (truth / sum 0.847); where one side dominates it under-prices by 25 % (each
boundary is clipped at 1, a junction reaches 1.3); the per-junction median error is 19 % whatever the exon lengths.
Second, the geometry: exons beside junctions are short (median 128 bp, 93 % below the longest fragment), so a
fragment across a junction holds bases of up to four exons where the gDNA crossings at its boundaries hold exon ends
and intron; under additive capture the rule is exact where both exons exceed every fragment (median error 3.6 %)
and 15 % off elsewhere. The count noise is Poisson on ~37 / ~340 expected crossings at a fully captured boundary,
plus calibration's bias at mid-density boundaries (+0.10 in the 0.3–0.7 band at `g05`); the posterior neither adds
nor removes it (posterior and plug-in on the same counts agree to 0.005 rms in every band). The unprobed class's
inflation at `g05` is not the junction's (`ISSUES: the-efficiency-posterior-floor-on-empty-pieces`). Refused on the
way: `ISSUES: a-junction-price-read-from-a-fitted-capture-map`, `ISSUES: the-junction-sum-on-unclipped-posteriors`.
The census, the attribution and the derivation: `~/Downloads/rigel_runs/prototypes/2026-09-23_junction_price/`.
THROUGH THE EM (`panel.py score` on the shipped tree, 16 ladder conditions, fractional assignment, stranded ×
capture ON): transcript and gene Σ|Δ| as a % of the true annotated RNA, and the pools as estimate − truth — gDNA
(intergenic included) / synthetic spans est/true / annotated. The floor is |base − base_reseed|; the third column
adds the population-corrected RNA pseudocount of `ISSUES: the-pseudocount-prior-is-biased-toward-gdna` (not
shipped; the prototype's run of 2026-09-23, which the shipped column reproduces within 0.03 points without it):

| | before the port | the shipped rule | + the corrected pseudocount |
|---|---|---|---|
| `g05` transcripts (floor 0.009) | 3.46 | 5.59 | 5.55 |
| `g50` transcripts (floor 0.003) | 5.50 | 6.18 | 6.20 |
| `g98` transcripts (floor 0.008) | 31.05 | 26.24 | 28.00 |
| `g05` / `g50` / `g98` genes | 0.44 / 1.91 / 16.87 | 0.64 / 2.12 / 14.96 | 0.62 / 1.45 / 18.54 |
| `g05` pools | +22,538 / 0.97 / −15,185 | +76,880 / 0.76 / −8,503 | +3,187 / 0.96 / +7,001 |
| `g50` pools | +35,572 / 1.22 / −68,501 | +172,932 / 0.28 / −64,052 | −25,216 / 0.98 / +28,209 |
| `g98` pools | −122,373 / 22.41 / −5,839 | −54,229 / 8.64 / +8,481 | −117,852 / 15.97 / +28,218 |

At `g05` the transcript error rises 2.13 points against 0.20 for the genes on the shipped rule (2.09 / 0.18 with
the corrected pseudocount): it is isoform allocation. Capture OFF the transcripts and genes move within ±0.03
points except at `g98` (unstranded 18.36 → 18.23, stranded 15.70 → 15.53); the deferred stratum's `g98`
transcripts read 90.63 → 115.07 (116.50 with the corrected pseudocount). The pools improve under the corrected
pseudocount while the isoforms regress either way (`TRAPS: judge-a-ruler-by-its-within-gene-spread`). The ceiling
it is priced against: the simulator's own lengths for every hypothesis with the corrected pseudocount read 1.89 %
and the synthetic pool 0.96× at `g50 ss.99 ON` (2026-09-21).

STILL OPEN FROM THE PORT'S REVIEW (2026-09-23, `g05 ss.99 ON` unless stated; recorded, not fixed):
* THE SHORT-INTRON FALLBACK keys on `S > 0` exactly (`S = gdna_region_eff_len`, the intron piece's gDNA contained
  support), so a piece with 0 < S < 1 (less than one expected start) holds essentially no information yet reads its
  own posterior (≈ the population's mean, 0.244, sd 0.058). Beside junctions 3,057 / 2,724 intron pieces (low /
  high side) have S = 0 and 2,759 / 2,819 have 0 < S < 1; treating S < 1 as no count moves 4,007 of 45,609
  junction prices (1,284 by more than 0.2, at most 0.79). A real library's shorter fragment tail shrinks the S = 0
  set (7,275 → 4,075 regions with a tail to 20 bp). No threshold is proposed: a threshold is a constant. MEASURED
  2026-09-23: beside 15–20 % of scored junctions the intron piece holds under 10 expected starts and the posterior
  moves the price 0.07–0.13 from its noise-free value, yet the intron terms made exactly noise-free leave the
  within-gene sd where it was (above) — a per-junction defect with no isoform cost on the ladder.
* THE CLIP AT 0: 137 transcripts carry more than half their share on junctions priced exactly 0; their contraction
  factor has median 0.0033 (minimum 0.00069), 30–130× below the adjacent-piece price's.
* ABOVE fl AND ABOVE THE PARENT: 2,593 of 8,750 annotated transcripts have `eff_em > fl` (at most 1.80×), and 23 of
  the 7,709 with a parent synthetic span read longer than it (median 1.010, at most 1.150) — possible with junction
  probes, no longer excluded by construction.
* gDNA'S REACH AT A REFERENCE END: gDNA's conserved share takes `UNBOUNDED_REACH` at every boundary, including one
  within a fragment of a reference end, where the accumulator clips (a boundary 10 bp from a reference start: true
  share 10.0, shipped 129.5) — gDNA's reach ruling (`region_geometry`: gDNA takes no reach argument), consequential
  only on short contigs.

THE CONSTRAINTS: the owner's rulings in `DESIGN.md` §7.2 (capture is local; junctions never pooled,
`ISSUES: pooling-junctions`; capture stays outside the EM) and §0b (robustness); no panel input (2026-09-20).
Refused on the way: `ISSUES: the-spliced-read-junction-price`, `ISSUES: a-junction-price-clipped-at-one`,
`ISSUES: a-junction-price-read-from-a-fitted-capture-map`, `ISSUES: the-junction-sum-on-unclipped-posteriors`. A
candidate is read on both columns of the first table at once — the within-gene sd toward the adjacent pieces', the
class mean on the gDNA component's — and only then through the EM; a candidate that only removes noise cannot pass
0.065 / 0.056. What would move it is an observable of the probe physics. `ruler_vs_truth.py --scale` (its
within-gene sd line and the junction-probed rows); `quant_accuracy.py` per stratum under
`--set em.assignment_mode=fractional`.

### a-junction-price-from-its-genes-own-rna
REFUSED 2026-09-26 by the owner, before it was built: one junction price per gene, read from the gene's own spliced
against its unspliced mature RNA. Its noise-free ceiling (the true per-gene price, every other length true) read
2.57 / 2.45 / 3.23 % transcripts at `g00` / `g05` / `g50 ss.99 ON` against the sum's 5.11 / 5.12 / 7.10. Why not: a
gene's isoforms start, end, include and skip a junction independently, so reading capture along a transcript needs
its abundance, which is the EM's output; capture efficiency stays outside the EM (`DESIGN.md` §7.2), and fitting the
two at once is not a narrow change. It would also relax the locality and no-pooling rulings. Do not re-propose a
junction price read from the RNA.

### a-junction-price-at-the-best-object-its-fragments-reach
REFUSED 2026-09-26 at the A/B. The owner's narrow candidate: a fragment across a junction binds its best single
probe part, so the junction takes the largest published efficiency among the objects its fragments touch in the
transcript's own coordinates — the junction's two boundaries, the exon pieces beside it and, where an exon is
shorter than a fragment, the next exons' boundaries — averaged over its placements with the deposit rule's weights;
no constant (`junction_reach/` in `~/Downloads/rigel_runs/prototypes/2026-09-26_capture_length/`). No EM, class-mean
L / Y ×10⁻³ and within-gene sd, shipped → this: `g05 ss.99 ON` isoforms 1.076 → 0.953 (spans / isoforms 1.039 →
1.173), within-gene 0.072 → 0.081; `g50` 1.112 → 0.989 (1.026 → 1.153), 0.057 → 0.049; the unprobed isoforms 2.84 →
59.7 / 1.21 → 21.8. Through the EM (fractional, stranded ON), transcripts / genes / spans est ÷ true: `g05` 5.56 →
4.16 / 0.59 → 1.13 / 0.98 → 0.99; `g50` 6.13 → 5.32 / 1.43 → 2.35 / 0.91 → 0.69; `g98` 26.18 → 27.29 / 17.60 → 19.48 /
13.30 → 12.88; `g00` and capture OFF unchanged (no reference). The isoform split improves while genes, pools and `g98`
regress, through the scale: every published efficiency is clipped at 1 and a fragment across a junction is captured
~1.19× the reference, so the best neighbour reads low and the isoforms fall 11–17 % below the spans and the gDNA
component. Refused on the way, noise-free: each edge read as an exon level (`2·c_edge − c_intron`) under the max
overshoots (class mean 1.20–1.31, within-gene sd 0.085–0.114 against the sum's 0.060); the boundary crossings split
by which side holds the fragment's midpoint price a junction no better than the sum (per-junction log sd 0.58
against 0.64), so no accumulator change is worth it.

### spike-in-references-read-as-depleted
CLOSED 2026-09-26 (owner), for now. A spike-in (ERCC) reference holds no gDNA template, so under capture its objects
read the landscape's depleted level and a probed spike-in's published `em_effective_length` is ~800× short
(ERCC-00131: 0.7 bp at `g05 ss.99 ON`). Its COUNT is right — Σ|Δ| 14.8 / 16.9 fragments of 24,747 / 12,982 at `g05` /
`g50 ss.99 ON`, the same with the true lengths — because the spike-in and its reference's gDNA component are
contracted alike. "No gDNA on the reference" cannot recognise one without a threshold: calibration places gDNA on
the ERCC contigs, up to 5,007 fragments (1,366 on one contig) at `g50 ss.50 ON` and 11.9 on one contig at `g05 ss.99
ON` (`ercc/` in `~/Downloads/rigel_runs/prototypes/2026-09-26_capture_length/`). If reopened: a reference covered end
to end by one single-exon transcript is recognisable from the index alone, with no constant; the gain is the
published yield only.

### the-em-answer-depends-on-where-it-starts
RULED 2026-09-26 (owner): accepted as a known limit and documented (`DESIGN.md` §3.1e; `MANUAL.md`, the FAQ on a
component that dies). The shipped start stays the coverage-weighted one, the most accurate measured.
What is left after the SQUAREM repair (`ISSUES: the-squarem-clamp-decided-which-components-live`, CLOSED): VBEM's
own start-dependence. With the clamp gone MAP reaches one answer from any start (204 fragments apart at
`g00 ss.99 ON`, all in loci still at the iteration cap), while VBEM still ends 8,096 fragments apart there (136
loci, 132 of them winner forks — a candidate dead in one answer and holding fragments in the other). Its E-step
weight ψ(α), with no per-component prior, penalises a small component by about −1/α, so components sharing fragments
race and the start picks the winner; the two answers' likelihoods differ (the uniform start's is higher in 121 of 136
loci), so it is not a flat valley. Candidate repairs, each measured or excluded below: a warm start that is itself
start-free (the MAP optimum, which is unique), or a per-component prior. Instrument: `em_start_lab.py` in
`~/Downloads/rigel_runs/prototypes/2026-09-25_em_start/` (every EM setting re-solved from one pre-EM state in one
process, which repeats itself to 0.0 fragments; `analyze.py` classifies each moved locus). The minus-strand
coverage-weight fix, which seeds only the start, landed on its own (`ISSUES:
the-minus-strand-coverage-weight-sat-one-fragment-length-off`, CLOSED); what it moved is this entry speaking.
MEASURED 2026-09-25 on `g00 ss.99 ON` (the landed solver's prototype): the 136 loci VBEM's two starts disagree on carry
EQUAL truth error — 385,299 (coverage start) against 385,088 (uniform), MAP's single answer 385,376; the coverage start
is closer in 50 loci, the uniform in 64, 22 tie — so what is left is reproducibility, not accuracy: the forks pick
among near-equivalent explanations of genuinely ambiguous fragments.
TWO START-FREE ARMS MEASURED 2026-09-25 (`~/Downloads/rigel_runs/prototypes/2026-09-25_vbem_start/`; the ladder
against the landed build, all 16, fractional, `ladder/compare_vbem.txt`): (1) VBEM STARTED FROM THE MAP OPTIMUM
(prototype `~/proj/rigel-em`, `RIGEL_PROTO_VBEM_START=map`: MAP first, then VBEM from the counts one E-step at θ_MAP
implies) — nearly start-free (the two starts 294 fragments apart at `g00 ss.99 ON` against 8,096; truth error there
936.6k–936.9k against 937.8k–938.0k), in-scope transcript Σ|Δ| −4.7 / +0.2 / −3.8 % (stranded OFF / stranded ON /
unstranded OFF), but genes +0.6 / +1.0 / +4.0 % and the gDNA pool's Σ|error| +3.3 / +2.1 / +3.2 %: it inherits MAP's
signature, more mass on annotated transcripts and less on the synthetic spans and on gDNA (`g50 ss.99 ON` −13.4k →
−15.5k), and it costs 30–50 % more pipeline time on the heavy conditions (measured beside the uniform arm under the
same load). (2) VBEM FROM THE UNIFORM START (a config value, `em.warm_start=uniform`): worse everywhere, transcripts
+2.0 / +0.1 / +0.5 %, genes +1.6 / +0.4 / +1.8 %. So which basin VBEM starts in is not neutral — the coverage start's
favour genes and gDNA, the MAP optimum's favour the transcripts off capture — and neither arm is a clean win. SIZED on an
unstranded row (`g05 ss.50 OFF`, the same lab): the landed VBEM's two starts end 115,159 fragments apart (134 loci, 58
by more than 100), and there the coverage start is the better one — truth error 1,371.7k against 1,376.9k — though the
uniform start is closer in more of the forked loci (76 against 44). There VBEM from the MAP optimum is neither
start-free nor better: its two starts end 19,404 fragments apart, and its truth error is 1,413.8k / 1,410.6k against
the landed 1,371.7k (+3 %) — it moves mass off the synthetic spans and gDNA onto the annotated transcripts, which is
why the transcript table alone read better. The coverage start has been the most accurate start measured on every
condition. Not
tried: choosing among VBEM's optima by VBEM's own objective (the evidence lower bound), which needs its derivation
under the grouped prior. A per-component reference prior is excluded by `EQUATIONS.md` §9b.1 (any lift above the
~0.16–0.47 activation threshold revives every shadow entity; the equal share was refused 2026-09-19).

### whole-counts-optimised-assignment
REFUSED 2026-09-26 (owner): the draw is kept. An assignment chosen for the most reads on their true origin under the
same counts wins reads and loses where they came from. The owner's framing had been: the draw a WARM START, the
repair a FEASIBILITY step (count error zero, `EQUATIONS.md` 13.3–13.4), and what is left an OPTIMISATION over
feasible assignments, wanted with a proof of what it reaches and a stated cost. The exact optimum was derived, built
and priced (`EQUATIONS.md` §13.5: prices, then shortest augmenting paths; exact by construction, certified by its own
dual). Prototype in
`~/Downloads/rigel_runs/prototypes/2026-09-26_whole_counts/` (`opt.c`, `exact.py`, `fidelity.py`): 9,000 random
loci equal brute-force enumeration, and dropping or corrupting the price update, or moving the wrong fragment, fails
that gate. MEASURED on two full libraries, the counts identical under every method (`g05 ss.99`, capture OFF / ON,
9.2 M fragments each):

| | reads on their true origin | reads misplaced by class | time |
|---|---|---|---|
| the fractional posterior (the model's own floor) | 66.53 / 72.41 % (expected) | 8.13 / 4.86 % | — |
| the draw (ships) | 66.53 / 72.43 % | 8.45 / 5.01 % | 0.4 / 0.2 s |
| the optimum, most reads expected correct (`Σ p`) | 68.09 / 74.07 % | 24.20 / 18.38 % | 130 / 7 s |
| the optimum, most probable (`Σ log p`) | 67.44 / 73.14 % | 16.09 / 11.29 % | 106 / 5 s |

"Misplaced by class" is `½ Σ |assigned − true|` over (candidate set, component) cells: whether each component's reads
come from where its reads really came from. An optimum is a vertex of the transportation polytope, so it hands groups
of look-alike fragments to one component instead of sharing them out — its gain and its cost are the same act. The
draw is the posterior sampled with the counts held, within 0.15–0.32 points of the model's own floor. The 130 s is
one locus (737 k fragments over 521 components, a median 47 effective candidates each: ten Gauss–Seidel sweeps still
left 334 k excess, then 110 s of shortest paths); a two-scale solve (the candidate-set classes first, then the
fragments) is designed and not built. REFUSED on the way: the auction algorithm with ε-scaling, which re-bid all
9.2 M fragments at each of six ε stages while near-identical fragments raised prices by ε per bid; it was stopped
after 25 minutes without finishing its first library. Greedy shortcuts do not substitute either: a greedy quota fill
strands 52–95 k fragments per library, and taking the largest re-weighted candidate instead of drawing strands
23–28 k. The draw's shortfall per read is the price of calibration, and an annotated BAM is read for where each
transcript's reads are.

### whole-counts-by-rounding-then-assigning
FIXED 2026-09-25 (`DESIGN.md` §3.1f; `em_solver.cpp`'s `assign_posteriors`; gates `tests/test_estimator.py`'s
`TestWholeCounts`, 14 cases). Whole-count mode gave each fragment to one transcript, drawn from its posterior after every
candidate under 1 % of THAT fragment was zeroed, to keep noise isoforms from collecting whole counts — but on one
fragment a real minority and noise look alike, so the floor took a minor isoform, a gene's unspliced RNA or low gDNA off
every fragment it shared (off capture +24 to +35 % transcript and 2.3–3.4× gene error against fractional). Now COUNT
FIRST: every component's fractional count, the locus's gDNA included, is rounded within its EM locus by largest
remainder; each fragment is drawn from its own posterior re-weighted by (count still owed ÷ posterior mass still to
come); a fragment the draw leaves with every candidate full is repaired onto the counts by a breadth-first chain of
moves, each fragment only to another of its own candidates. The whole counts ARE the rounded fractional counts,
whatever the seed, whenever that rounding is reachable; where it is not (Hall's condition fails — two components one
fragment can reach, both rounded up), a fallback of the same augmenting paths settles every count within one of its
fractional count, which is always reachable. The derivation — termination, the success condition, the cost, and the
optimality it gives up — is `EQUATIONS.md` §13. The fallback was added after a search found the exact-target repair
alone leaving a count two from its expectation; on 30,000 draws over 6,000 small loci built around an unreachable
core (exact targets unreachable in 9,310), no count left its floor or ceiling and no locus total moved.
`assignment_mode="map"` and `assignment_min_posterior` are deleted (owner). The gates hold the counts
to an independent largest-remainder reference on four loci (one ordered so the draw strands and the repair must work;
one where every remainder is under a half, so rounding each count alone would lose a fragment), a minority under 1 % on
every fragment keeping its count, a transcript expecting under half a fragment getting none, seed-independent counts
with seed-dependent fragments, every fragment on one of its own candidates, and — for identical fragments — the first
and the last drawn with the same odds (the draw is an urn). Each of six perturbations fires its own gates: the old
floor, no repair, rounding each count alone, no re-weighting, no randomness, and a lost per-read output. LADDER, the
shipped whole counts against fractional (all 16, `~/Downloads/rigel_runs/prototypes/2026-09-26_whole_counts/ladder/`):
transcript Σ|Δ| −0.02 / +0.04 / −0.01 % (stranded OFF / stranded ON / unstranded OFF), genes −0.04 / −0.00 / −0.13 %,
every condition within ±0.07 % on transcripts and ±0.32 % on genes, the gDNA pool within 51 fragments; transcripts
reported for zero-truth isoforms 3,707–5,849 → 123–537. Fractional output is bit-identical (the goldens regenerate
unchanged). The design was measured first (`~/Downloads/rigel_runs/prototypes/2026-09-26_whole_counts/`, every EM fragment's posterior dumped
and joined to its true origin): the draw alone — before the repair — landed transcripts within 0.3 % of fractional;
a greedy quota fill stranded 52–95k fragments per library, and the exact optimum, priced in
`ISSUES: whole-counts-optimised-assignment`, trades region fidelity for reads assigned correctly.
The earlier record: Whole-count mode (`assignment_mode="sample"`) gives each fragment to one transcript, drawn from
its posterior after every candidate under 1 % (`assignment_min_posterior`) is zeroed, to keep noise isoforms from
collecting whole counts. On one fragment a real minority and noise look the same — a small share — so the floor takes
a minor isoform, a gene's unspliced RNA or low gDNA off every fragment it shares with a dominant candidate: off
capture the synthetic pool loses 16–29k fragments, transcript error rises 15–23 % and gene error 2.0–2.7× (fractional
mode is untouched; `~/Downloads/rigel_runs/prototypes/2026-09-24_cut_inventory/`). A cutoff relative to the best
candidate (the owner's first idea) spares a fragment with many near-equal candidates but not this. The design (owner:
"I like it in theory"): COUNT FIRST — each component's EM expected count, rounded within its locus so the locus total
is exact (a transcript expecting 0.3 fragments gets 0, one expecting 300 gets 300) — THEN ASSIGN each transcript that
many fragments, the ones it most plausibly produced (a transport problem over the unit posteriors, per locus). Every
fragment still goes to exactly one transcript, with no seed. The rounding half is measured
(`~/Downloads/rigel_runs/prototypes/2026-09-25_floor_study/`): transcripts reported for zero-truth isoforms 4,186 →
196 at `g05 ss.99 OFF` and 4,704 → 403 at `g50 ss.99 ON`, with transcript and gene error unchanged. The assignment
half needs a design and a C++ prototype, then an A/B against today's floor and a relative cutoff. False-positive MASS
(90–96 % of it in a few dozen to a few hundred transcripts that each take a small share of many fragments) is the EM's
own and no assignment rule separates it from a real minority. THE ASSIGNMENT HALF, PROTOTYPED 2026-09-25
(`~/Downloads/rigel_runs/prototypes/2026-09-26_whole_counts/`, `summary.txt`): every EM fragment's posterior row
dumped from the solver (it re-sums to `count_em` to 7e-8) and joined to its true origin by read name, on `g05` and
`g50` of the three in-scope strata. Quotas are each component's fractional count (transcripts and the locus's gDNA)
rounded within the locus. The COUNT-EXACT DRAW — one pass, each fragment drawn from its posterior re-weighted by
(quota left ÷ posterior mass left) per candidate, a fragment whose candidates are all full keeping its best — lands
transcripts within +0.19 to +0.26 % of fractional off capture and −0.03 to 0.00 % on it, genes +2.0 to +4.1 % off
capture and −0.3 % on it (the 916–2,634 fragments per library it could not place), reports 190–400 zero-truth
transcripts against fractional's 4,186–4,704, gives each read a draw from its own posterior (reads correct = the
posterior's own expectation, 65.2–78.4 %), and runs in 0.3–0.5 s per library; visiting each locus's fragments in a
random order changes nothing (the simulator's BAM is grouped by origin). Against it: today's floor +24.5 to +35.5 %
transcripts and 2.3–3.4× genes off capture (+1.6 to +3.2 % / +6 to +45 % on); the relative cutoff +16 to +23 % /
1.9–2.9×; a plain draw, with or without dust removed, +5.6 to +7.6 % / +13 to +16 %; a greedy quota fill strands
52–95k fragments per library (their candidates already full) and costs +40 to +53 %; taking the largest re-weighted
candidate instead of drawing strands 23–28k and costs +7.6 to +13.9 %; the argmax assigns 72–85 % of reads correctly
but multiplies transcript error 5–20×. The optimum under the same counts is priced in
`ISSUES: whole-counts-optimised-assignment`.

### the-scan-fraction-banks-are-not-reproducible
REFUSED 2026-09-25 (owner): bit-identical results are not a goal for now; the tally keeps its float sums and the
scan's thread noise is accepted. The facts it was refused with: the scan's workers take batches from one queue as
they finish and each sums its own float64 fraction banks, so the last bits change from run to run at any scan thread
count above one; at `g50 ss.99 OFF` (two runs, EM on one thread, fractional) the six banks differ at ~1e-15, 64 loci's
gDNA prior at ~3e-16, and the EM's forks carry that to 317 transcripts by up to 0.88 fragments and nascent parents by
up to 42; with the scan on one thread all 214 stage arrays repeat bit for bit
(`~/Downloads/rigel_runs/prototypes/2026-09-26_repro/FINDINGS.md`). The repair priced: exact sums (every deposit a
whole number of 2^-94 units in an unsigned 128-bit integer, rounded once on export) made the whole pipeline repeat at
the default thread budget (213 of 214 arrays; the last is the buffer's row order), moved the answer once by about one
run-to-run spread, and cost ~50 MB per scan worker at human scale plus an exact-summing specification
(`exact_sum.patch` beside the findings). `--threads 1` gives bit-identical output (`MANUAL.md`).

### the-minus-strand-coverage-weight-sat-one-fragment-length-off
FIXED 2026-09-25 (`scoring.cpp`'s `coverage_weight`; gate `tests/test_pipeline_routing.py`: every fragment's weight
against the trapezoid coverage model integrated exactly — 29 placements on one exon, a short transcript and a spliced
one, both strands — beside a strand-mirror test and a multimapper test over 14 offsets each). The EM's warm start
weights each candidate by where the fragment sits on the transcript: the trapezoid `min(x, w, L − x)`,
`w = min(f, L/2)`, over the fragment's bases. On a minus-strand transcript the old rule flipped the genomic start
into transcript space and read the fragment forward from it, but there the genomic start is the fragment's 3′ end,
so the window sat one fragment length off: `[L − s, L − s + f)` in place of `[L − s − f, L − s)`. The rule now
measures the start as the resolver does (a start in an intron back from the next exon, `ISSUES:
a-start-in-an-intron-was-measured-from-the-wrong-exon`), reverses the window on the minus strand, and cuts it at both
transcript ends, so an overhang weighs nothing. Each perturbation fails the gate: the old flip 50 cases, an intron
start snapped to the previous exon 4, the overhang shifted 6, no cut 20; the flip alone, without the other two, 8.
It seeds only the start, so it reaches the answer only through VBEM's start-dependence (`ISSUES:
the-em-answer-depends-on-where-it-starts`). LADDER against the SQUAREM-landed arm (all 16, fractional;
`~/Downloads/rigel_runs/prototypes/2026-09-26_coverage_landed/ab_table.txt`): in-scope transcript Σ|Δ| 390.0k →
391.6k (+0.41 %, stranded OFF), 1,489.8k → 1,487.8k (−0.13 %, stranded ON), 395.5k → 396.1k (+0.16 %, unstranded
OFF); genes +0.32 / −0.02 / +0.27 %; the gDNA pool's Σ|error| within 142 fragments. Rows move by −730 to +1,159
transcript fragments (`g50 ss.99 OFF` +1.24 %), above the run-to-run spread (same-build reseed pairs differ by 0–345,
median 7), so the moves are real: the fix changes only the start, and VBEM's answer depends on its start. (Corrected
2026-09-25: first recorded as inside the floor, read against a reseed arm from an older build — a cross-build
difference, not a floor.) Goldens moved by at most 2.3e-8 relative.

### the-squarem-clamp-decided-which-components-live
FIXED 2026-09-25 (`DESIGN.md` §3.1e; `em_solver.cpp`'s `backtracked_squarem_step`, VBEM and MAP; gate
`tests/test_em_start_independence.py`: three real loci the clamp forked by 7–72 fragments between the coverage and the
uniform start, every case failing on the old solver, and each mode's cases alone failing when its backtrack is
removed). SQUAREM's extrapolation overshot a shrinking component below the floor and CLAMPED it there, and a component
at the floor takes no responsibility and no share of the evidence-proportional prior, so it never came back: which of
the components sharing fragments survived depended on the warm start. `MANUAL.md`'s FAQ said a component with genuine
read support recovers — true only for one holding deterministic fragments of its own. It surfaced through the
minus-strand coverage-weight fix, which changes only the warm start yet moved the ladder ±1–2 %, and at
`em.iterations=10000`, `em.convergence_delta=1e-9` the two trees still differed by −1.8k to +5.8k transcript
fragments. Now the step is halved toward the plain double step until no component that step keeps alive is carried
below the floor. MEASURED on `g00 ss.99 ON`, every setting re-solved from one pre-EM state in one process
(`~/Downloads/rigel_runs/prototypes/2026-09-25_em_start/`): the start moved 41,343 fragments under VBEM and 22,447
under MAP, every moved locus clamped and none a flat valley; with backtracking MAP moves 204 (loci at the iteration
cap only) and VBEM 8,096 (its own, `ISSUES: the-em-answer-depends-on-where-it-starts`); from one start the repaired
MAP beats the clamp's likelihood in 116 of 120 loci; false gDNA there (truth 0) VBEM 4,287 → 177, MAP 1,187 → 158.
LADDER, the landed build against the code it replaces (all 16, fractional; `landed/landed_vs_committed.txt`):
in-scope transcript Σ|Δ| 407.7k → 390.0k (−4.3 %, stranded OFF), 1,496.1k → 1,489.8k (−0.4 %, stranded ON),
409.6k → 395.5k (−3.4 %, unstranded OFF); genes −6.1 / −0.3 / −7.6 %; the gDNA pool's Σ|error| −2.3 / −5.8 / −2.3 %,
and every `g00` row's false gDNA ~4.5k → ~150. It unmasks an under-call at `g05`, which the clamp's false gDNA had
partly offset: stranded OFF −1.3k → −4.1k, unstranded OFF −4.0k → −6.6k. MAP, with or without the repair, is worse
than VBEM on the gDNA pool (`g50 ss.99 OFF` −15.2k against −9.5k), so VBEM stays. Goldens moved by at most 3.3e-8
relative (the plain step is now copied at step 1, not recomputed). COST, MEASURED 2026-09-25 on VCaP at 8 threads
(`~/Downloads/rigel_runs/prototypes/2026-09-26_squarem_timing/`): counted, so independent of machine load, +8.8 % SQUAREM
iterations and +13.7 % E-step work, loci at the 333-cycle cap 103 → 123 (the largest locus, 12,890 components, now
among them: 299 → 333); the halvings themselves are ~1e8 component checks, negligible. Three interleaved wall-clock
pairs, taken beside another workload, put the locus-EM stage at +12 % on average (−3 to +27 % per pair, ~+3 s of
~23 s) and the whole run unchanged within ±7 %.

### a-start-in-an-intron-was-measured-from-the-wrong-exon
FIXED 2026-09-25 in the resolver (`resolve_context.h`'s `tx_frag_length`, the one length definition of `DESIGN.md`
§3.4; gate `tests/test_transcript_space_fl.py`'s `TestLengthRuleByEnumeration`, a base-by-base count over every
two-mate fragment with its ends at or beside an exon boundary). The EM's per-candidate fragment length projected both
ends into transcript space and measured an end lying in an intron forward from the PREVIOUS exon's end — right for the
fragment's end, wrong for its start, which landed one intron length downstream: the length read |true − intron| (a
200-bp fragment starting 1 bp before a 2.8-kb intron's end read 2,610), and a length of exactly 0 was stored as
missing (−1), which the scorer reads as log-likelihood 0, the best fit there is. An end whose last base was an
intron's last base read as ending at the next exon's start and dropped the whole intron. Two hand-written tests
asserted the defective values (398 and 2,006; now 3,397 and 5,005), and neither boundary error is caught by any
hand-written test: perturbing either one fails only the enumeration. Before the fix
(`~/Downloads/rigel_runs/prototypes/2026-09-24_cut_inventory/`), 1.35–10.7 % of multi-block candidate entries per
condition had a wrong length, median error 1–1.2 kb, and 1.9–2.6k candidates read −1. MEASURED on all 16 conditions
(fractional, the worktree build against the shipped `base` arm;
`~/Downloads/rigel_runs/prototypes/2026-09-25_start_endpoint/ab_table.txt`): transcript Σ|Δ| summed per in-scope
stratum −560 (stranded OFF), −4,277 (stranded ON), −887 (unstranded OFF); the largest row `g50 ss.99 ON` 301.5k →
298.0k; two rows worse beyond the run-to-run spread (`TRAPS: the-deliverable-is-not-reproducible-by-default`),
`g05 ss.99 OFF` +496 (+0.41 %) and `g00 ss.50 OFF` +168; genes within ±1.2 %; the gDNA pool within ±1k everywhere.
The inventory's estimate that the defect pushed ~1.1k gDNA fragments per `g98` capture-ON condition toward RNA did
not survive: the pool moved by 22 there.

### calibration-and-the-em-disagree-on-what-can-be-gdna
FIXED 2026-09-24 in the EM (`DESIGN.md` §3.1d; `scoring.cpp`'s `gdna_can_explain` and `gdna_competes`; gates in
`tests/test_pipeline_routing.py`); the leftovers moved to `ISSUES: splicing-artifacts` and
`ISSUES: multimapper-intergenic-alignments`. Calibration and the EM decided by different rules whether a fragment
can be gDNA: calibration offered every fragment its unbroken reading and drew among the survivors, while the EM
denied gDNA to every label but "unspliced" and removed it from a multimapper's every hit when one hit was spliced.
Now an unspliced, implicit or artifact alignment may be gDNA, a multimapper may if any alignment may, and gDNA is
pruned exactly as an RNA candidate is. The implicit splice needs no new likelihood: gDNA's term is the unspliced term
at the footprint and the fragment length decides (derived three ways; `DERIVATION.md` in
`~/Downloads/rigel_runs/prototypes/2026-09-24_spliced_rule/`). MEASURED on a one-gene toy on the shipped solver:
denying implicit fragments gDNA put the error on the synthetic span and the annotated transcripts under every
prior — +6 % at a 60-bp intron and gDNA 50 %, +39–62 % at 90 %, +13 % at 150 bp and 90 %. On the ladder, all 16
conditions (fractional, VBEM, on the count form; `~/Downloads/rigel_runs/prototypes/2026-09-24_implicit_gate/ab/`):
gDNA error at stranded ON `g50` −28.0k → −14.1k, `g98` −112.8k → −100.5k, `g05` +1.3k → +2.6k, off capture within
±1.0k; transcript error `g98 ss.99 ON` 54.2k → 50.3k, elsewhere within ±1.7k; gene error lower in 12 of 16 rows.
The gDNA gained at `g50 ON` came from the synthetic spans (1.00× → 0.91×), not the annotated excess (+28.3k →
+27.6k). A first form that limited implicit splices to the maximum fragment length admitted 18k–55k units per
condition whose gDNA reading is about −708 nats, each counting in the prior's strength, and raised `g00 ss.99
OFF`'s false gDNA to +7.0k against +4.4k under the pruning.

### the-pseudocount-prior-is-biased-toward-gdna
FIXED 2026-09-24 by the count form (`DESIGN.md` §3.1c, `EQUATIONS.md` §9b.3, `pipeline.em_pseudocounts`; gate
`tests/test_em_pseudocounts.py`). The EM's two pseudocounts were calibration's gDNA and UNSPLICED RNA counts, but
θ is each component's share of every fragment in the locus — the deterministic spliced fragments enter every
M-step — so the unspliced split stated as the whole pool's leaned toward gDNA in proportion to the spliced
share: its centre read 63.9 % against 50.0 % at `g50 ss.99 ON`, and it over-called an exact toy by 70–944 of 900
as the spliced share ran 11–60 %. Now `P_g = U·min(G_c/N_c, 1)`, `P_R = U − P_g`, the share taken over the
`N_c` fragments calibration counts from (on the ladder `N_c = N`). MEASURED on all 16 conditions
(fractional, VBEM with a MAP pair, the odds changed and the strength held): gDNA error at `g05` / `g50`
unstranded OFF +59.2k → −3.1k and +225.7k → −9.7k, stranded OFF +24.6k → −1.2k and +128.5k → −8.7k, stranded ON
+76.9k → +1.3k and +172.9k → −29.4k (the spans 0.28× → 1.01×); gene error at `g50` 18.3k → 14.4k, 13.6k →
11.4k, 102.6k → 70.2k; the transcript table within 2 %; `g00` unmoved. `g98` worse in every stratum (−18.7k →
−62.8k, −3.1k → −37.8k, −54.2k → −118.3k; genes ON 29.0k → 36.2k): the lean was offsetting
`ISSUES: rna-prior-floor-at-pure-gdna-loci` and the capture likelihood's lean toward RNA, the owner's next
targets. The deferred unstranded × capture-ON stratum moves the other way (`g50` −125.6k → −589.5k). MAP moves
the same way as VBEM, 2–8k further toward RNA. The per-locus A/B, census and tables:
`~/Downloads/rigel_runs/prototypes/2026-09-23_pseudocount/` (`DERIVATION.md`, `ab16/table16.py`).

### a-count-once-density-prior-for-the-strength
REFUSED 2026-09-24 as a buildable rule for the pseudocount's strength; the principle stands. A prior restating
calibration's full answer counts the locus's fragments twice, because the E-step reads them again; the rule that
counts them once puts a Gamma prior on the locus's gDNA density at calibration's estimate from OUTSIDE the locus,
and its exact M-step is `out_g = (raw_g + a − 1) / (1 + a/μ)` with RNA untouched (matched on the shipped MAP solver
to 0.01 fragment). Killing numbers: under capture, where the prior carries the answer, no outside estimate exists —
each object's capture efficiency is read from its own count, so the "outside" centre equals calibration's own gDNA
count to 0.1–1.3 % on the test chromosome, and the locus's flanks are off target, 75–440× below its gDNA. Off
capture it carries three choices that move it at first order: how over-dispersion correlates within a locus (the
test chromosome's `g50 ss.99 OFF` per-locus error 2,375 or 3,349 against the shipped strength's 2,461), the
moment estimate's `α = ∞` boundary (three of four ladder conditions), and MAP against mean against VB (30 % apart
at `a` = 3.3 on the unstranded ridge). At the ladder's own reading it is worse than the shipped strength on all
three capture-OFF test conditions (3,349 / 4,754 / 1,542 against 2,461 / 4,262 / 1,189), and on the shipped
grouped VBEM it erases the split within RNA (a 6-fragment isoform 4.75 → 7.37). Scripts: `toy_gamma_form.py`,
`refute_kappa_*.py` in `~/Downloads/rigel_runs/prototypes/2026-09-24_spliced_rule/`.

### a-junction-price-read-from-a-fitted-capture-map
REFUSED 2026-09-23 at the no-EM gate. The general form of conservation of bases: capture additive over bases, every
region a flat per-base capture read from gDNA's own objects (non-negative least squares over every region's
contained and every boundary's crossing efficiency, each row weighted by its support, per reference), and each
junction priced by routing the transcript's own placements through that map — the shipped sum is its limit where
the pieces beside a junction exceed every fragment. Within-gene sd at `g05` / `g50 ss.99 ON`, on the transcripts'
own routed shares where the sum reads 0.050 / 0.054 on the simulator's noise-free counts and 0.074 / 0.058 on
calibration's: the map 0.060 / 0.064 and 0.089 / 0.073, the unprobed class inflated further (1.52 / 1.50
noise-free against 1.00 / 1.01). With the TRUE per-region probe coverage under additive capture it fixes the
geometry (junction rms 0.58 → 0.19, 1.013 on average); against the simulator's capture it prices junctions 2.1×
high, because probes designed on different isoforms overlap, and overlapping probes add under additive capture but
not where a fragment binds its best single part (`ISSUES: the-junction-price-is-noisy-within-a-gene`).

### the-junction-sum-on-unclipped-posteriors
REFUSED 2026-09-23 (`ruler_vs_truth.py --module`, the prototype's `junction_rulers.py`). The sum by conservation of
bases is linear in density, so it was taken on each of its four terms' unclipped posterior `E[ρ / ρ_ref]` in place
of the published `E[min(ρ / ρ_ref, 1)]`, every other object unchanged. The class medians move toward 0 (partial
mRNA +0.035 → +0.024 at `g05`, +0.019 → −0.020 at `g50`) and the within-gene sd widens, 0.074 → 0.080 / 0.060 →
0.066 (on the simulator's noise-free counts 0.050 → 0.052 / 0.054 → 0.060). The junction keeps the published
clipped efficiencies.

### the-gdna-component-length-rule-differs-from-the-transcripts
RESOLVED 2026-09-23 by the one shared rule (`DESIGN.md` §7.2, `EQUATIONS.md` §11): every EM component — the locus
gDNA component, every synthetic span, every annotated transcript — has as its capture-contracted length its
conserved share of every object its fragments deposit on times that object's capture efficiency. The defect:
transcripts and spans took the per-base rule (their own bases at their pieces' efficiencies with the fragment-end
taper) while the gDNA component took region and boundary objects with the crossing supports converted by the
pooled `q`, and the two rules disagreed by 14–16 % on capture. Killing numbers (`ruler_vs_truth.py --scale`, L / Y
×10⁻³ against the simulator's exact yields, gDNA component / spans / isoforms, and the largest gap): `g50 ss.99 ON`
1.108 / 0.933 / 0.957 — 18.7 % → 1.126 / 1.136 / 1.111 — 2.3 %; `g05 ss.99 ON` 1.107 / 0.932 / 0.944 — 18.7 % →
1.087 / 1.105 / 1.075 — 2.8 %, in the prototype (calibration's counts on the count frame as clipped plug-ins). The
shipped tree, with the published posterior efficiencies (width step 1), reads 1.131 / 1.141 / 1.112 — 2.6 % and
1.108 / 1.118 / 1.076 — 3.9 % at the two rows. Within each locus, at the two rows, in the prototype:
gDNA/isoforms 1.149 / 1.153 → 0.983 / 0.986, gDNA/spans 1.202 / 1.208 → 0.984 / 0.981, spans/isoforms 0.978 /
0.977 → 1.014 / 1.014. The true gDNA on every object reads the same ~2 %, the rule's own residual
(`ISSUES: ruler-witness-geometry-on-transcript-panels`). The port reproduces the prototype A/B'd to 9e-15 on
transcripts except the 113 at the fl floor (shorter than every fragment), and on the gDNA component except the loci
sharing a region with another locus (34 of 1,216 at `g50 ON`, 3.7 % of the gDNA prior; `P_g`-weighted totals
0.999998), where the length now takes the count's own overlap weights. What it cost through the EM is two open
entries: the pseudocount's bias no longer cancelled (`ISSUES: the-pseudocount-prior-is-biased-toward-gdna`) and
the isoforms split by the junction price's noise (`ISSUES: the-junction-price-is-noisy-within-a-gene`). Refused on
the way: `ISSUES: synthetic-spans-on-the-object-rule`, `ISSUES: the-per-base-rule-for-every-component`,
`ISSUES: the-pooled-q-in-the-gdna-length`.

### pooling-junctions
REFUSED 2026-09-23 by the owner, before it was built: a posterior pulling each junction's price toward its
neighbours' price times a ratio learned across all junctions, proposed to cut the sum's noise
(`ISSUES: the-junction-price-is-noisy-within-a-gene`). Some junctions are probed across the junction itself and
some are not, so a ratio pooled over junctions is theoretically invalid: a gain on the ladder would be chance, and
another probe design could read much worse (`DESIGN.md` §7.2). Do not re-propose a junction price that pools
across junctions.

### the-junction-sum-over-prices-separately-probed-exons
CLOSED by landing 2026-10-05 (owner: the interim cap now; `DESIGN.md` §7.2, `EQUATIONS.md` §11; gate
`test_a_junction_is_captured_no_more_than_a_fully_captured_piece`, verified failing on the uncapped sum and fired by
three perturbations). THE DEFECT: the test chromosome's stranded × capture-ON transcripts read 13.89 % on the tree
against 8.58 % on the tagged 0.7.1 (genes 4.01 against 2.45 %; ss 0.70 × ON 15.09 against 10.28 %), the only in-scope
rows where 0.7.1 led besides the RNA-short one. The junction price `c_lo + c_hi − ½(c_intron,lo + c_intron,hi)`
assumed capture adds over a fragment's bases; a fragment binds its best single probe part, so on exons tiled to their
edges each boundary read 0.835 of the exon level, the sum priced probed junctions at a median 1.65, the multi-exon
isoforms sat 7 % long against the spans and the gDNA component (spans/isoforms 0.934 against 0.989 in 0.7.1), mature
RNA moved to the spans (nascent +110k against a stratum truth of 21k) and isoforms mis-split by junction share. THE
KILLING NUMBERS (in process, pinned like the CLI, only `_cut_efficiencies` changed), stranded × ON transcripts / genes,
0.7.1 → uncapped → capped: test chromosome 8.58 / 2.45 → 13.89 / 4.01 → 7.19 / 0.96 (ss 0.70 × ON 10.28 / 3.32 →
15.09 / 5.75 → 6.91 / 1.51); junction-probed panel 40.01 / 1.46 → 31.38 / 0.85 → 32.12 / 1.15; one probe centred per
exon 74.76 / 63.06 → 24.43 / 3.19 → 24.39 / 3.15; ladder 9.91 / 2.52 → 6.14 / 1.12 → 3.85 / 1.38; RNA-short 20.33 /
3.52 → 19.82 / 1.52 → 19.97 / 1.53; RNA-long 7.38 / 1.62 → 6.17 / 0.62 → 4.48 / 0.68. Capture-OFF and `g00` rows
bit-identical; 0.7.1's own lengths in today's EM reproduce 0.7.1 (8.91 %), so nothing else moved. THE COST: genes on
transcript-coordinate panels (ladder 1.12 → 1.38 %, junction-probed 0.85 → 1.15 %), where a probe part spans the
junction, a junction fragment is captured about 1.19× a contained one and the cap under-prices it. The topology (two
probes or one spanning the junction) is not identifiable from gDNA where both sides are captured, and is not modelled
(owner, 2026-10-05). Exon-level prices (0.7.1's imputation, the larger or the mean of the two exons) collapse where no
junction fragment reaches the probe (72.8–75.0 % on the centred panel) and are refused. Instrument:
`~/Downloads/rigel_runs/prototypes/2026-10-05_test_strON_regression/` (`arms_e2e.py`, `ruler_variants.py`,
`scale_compare.py`, `decompose.py`).

### a-junction-price-clipped-at-one
REVERSED 2026-10-05 and landed (owner; `ISSUES: the-junction-sum-over-prices-separately-probed-exons`): the refusal below read the ruler alone and was never priced end to end. REFUSED 2026-09-23. The sum by conservation of bases is floored at 0 and not clipped at 1 — each boundary is an
unclipped half-exon level — and at `g50 ss.99 ON` 18,308 of 45,609 junctions price above 1 (max 1.999) and 4,440
below 0 (`g05 ss.99 ON`: 17,167 above, 7,409 below). A clip at 1 shortened 25 % of transcripts by a median 13 %
(Σ 0.965 of the validated rule): a second clip on pieces already clipped inside their own posteriors.

### the-spliced-read-junction-price
REFUSED 2026-09-23: the owner's proposal to impute a junction's capture from the transcripts' own reads left and
right of it (derived with the correction ÷ `f_spliced`, not ×) — exact in principle, and starved. The pieces
beside junctions have a median length of 107 bp; 58 % of single-isoform junctions hold no unspliced read within a
fragment; and "exonic on the strand" neighbourhoods reach retained-intron-only regions and inflate the price
(per-transcript L up to ×60).

### the-crossing-apportionment
DELETED 2026-09-23 with the one shared rule (`capture_efficiency.capture_efficiencies`, no `shares` argument):
a region's efficiency had read its own contained count plus the crossings at every boundary within a fragment's
reach, apportioned to the pieces its fragments cover (`ISSUES: ruler-multimapper-floor-caps-the-correction`). Each
object now reads only its own count — a region its gDNA contained count on its contained support, a boundary its
gDNA crossing count on its crossing support — and the length prices every crossing at the boundary that holds it,
so a region never borrows the crossings around it. The apportionment had read a tiny piece's efficiency off its
edge crossings, which the rule now prices at those boundaries directly; a piece too short to contain a fragment
has no share, and its efficiency (the population's) multiplies nothing. A deletion by structure, not by a number.

### the-pooled-q-in-the-gdna-length
REPLACED 2026-09-22/23 by gDNA's own conserved share at each boundary (`CalibrationResult.gdna_boundary_conserved_len`,
`EQUATIONS.md` §11). The gDNA component's length had converted a boundary's crossing support by
`q = boundary_mass_per_crossing`, a mass per crossing pooled over gDNA and RNA and so RNA's where RNA dominates a
boundary: against the oracle's gDNA-only payload it overstated gDNA's mass by 1.4 % overall at `g50 ss.99 ON`,
3.4 % at exon|exon boundaries, 8.4 % where RNA is 50–90 % of the crossings, 5 % at `g05 ss.99 ON`, and it read
transcripts 9 % high off capture. Storing the unspliced boundary mass per strand column would give gDNA's own mass
exactly on stranded data (a 2×2 solve per boundary; unstranded data cannot split it) — not built: the length
needs neither. The two pseudocounts still convert a boundary's mass by `q`.

### the-junction-price-at-its-neighbours-mean
REFUSED 2026-09-22 as the scale. Against the simulator's own junction capture (60 probed transcripts, 395
junctions) the junction's two boundaries' mean prices it at 1/1.68, the adjacent pieces 1/1.73, the nearest
measured object along the transcript 1/1.74 and the transcript's own measured mean 1/1.24 (a transcript-level
anchor, barred by the owner's locality ruling, `DESIGN.md` §7.2), where the sum by conservation of bases reads
0.959. On the class scale, with the junction at its two boundaries or its adjacent pieces the gDNA component and
the spans agree to 1 % (0.991 / 0.993 with the true gDNA, 1.137 / 1.147 with calibration's) while the transcripts
sit 10–20 % below (0.910 / 0.970 with calibration's efficiencies). KEPT AS DATA: priced at its adjacent pieces a
junction is as precise within a gene as the per-base rule (within-gene sd 0.057 / 0.041 at `g05` / `g50 ss.99
ON`) and brings the class gap back (`ISSUES: the-junction-price-is-noisy-within-a-gene`).

### the-per-base-rule-for-every-component
REFUSED 2026-09-22 by the owner's ruling (`DESIGN.md` §7.2): a length over each component's own bases at its
pieces' capture efficiencies with the fragment-end taper omits the boundaries that own the crossing mass — any
fragment that crosses a boundary is deposited into it, and a piece shorter than a fragment is seen only through its
boundaries. Its no-EM reading was one scale for gDNA and isoforms (0.956 / 0.957 ×10⁻³ at `g50 ss.99 ON`) with
the spans 2.5 % low (0.933), and only because its clipped per-piece efficiency under-priced gDNA's short probed
pieces by about as much as it under-priced the transcripts' junctions (with the locus gDNA's taper on the outside
of its footprint, 3.0 / 2.0 / 0.7 % from one scale at `g50 ss.99 ON` / `g05 ss.99 ON` / `g50 ON` at `on_fraction`
0.10, 2026-09-21). Within a probed class its spans read 0.950 / 0.949 / 0.964 of the isoforms (`partial ≤ ½` /
`partial > ½` / `probed ≥ 0.9`) and the locus footprints 0.978 / 0.917 / 0.521, and no region-constant per-piece
efficiency closed it (the sampler's own covering truth 1.271 / 1.155 / 1.124, its contained truth 1.150 / 1.007 /
1.059, the probed fraction 1.122 / 1.022 / 0.979; the clipped estimator was the closest; 2026-09-21).
Through the EM with the corrected pseudocount it read 5.49 → 5.14 % transcript error at `g50 ss.99 ON` with the
synthetic pool 1.22× → 1.65×. Re-measured 2026-09-23 beside the one shared rule (fractional, stranded × capture
ON, with the corrected pseudocount on both): transcripts `g05` / `g50` / `g98` 3.47 / 5.15 / 30.43 % against 5.55 /
6.20 / 28.00 % — refused by structure (owner's ruling), not by this number (`DESIGN.md` §0b). Its form for the gDNA
component alone — the locus's bases at the regions' efficiencies — was refused 2026-09-20: it drops the boundary
objects the count keeps; test chromosome `g50 ss.99 ON` gene-level Σ|Δ| 25,633 → 38,174 against 23,967 with the
object form, the gDNA pool +5.5 % against −0.1 %. The taper and the per-base frame are deleted.

### background-abundance-pair-unruled
DELETED 2026-09-22 (the cleanup; owner: "I don't know what measured_total is for"). `CalibrationConfig.background_abundance`
chose which (counts, exposure) pair the pooled gDNA background took — the shipped `"contained"` pair or the START/END banks
over the region's own length. An alternative nobody ruled on is a tunable, and it went with its refusal path, its
`counts_exposure` parameter and its three tests; the contained pair is the only one. `rename_identity.py --check`
bit-identical.

### the-fixed-gdna-level
REFUSED 2026-09-21 as the shipped kernel (33 converged arms on the ladder, one EM thread, fractional). Holding the
gDNA component's count at calibration's locus count every iteration and dropping the RNA pseudocount repairs the
capture-OFF pools exactly as the corrected pseudocount does and halves the `g00` false positives (+3,913 → +1,902),
but at `g98 ss.99 ON` it collapses on every length rule — 46.1 % (object) and 45.4 % (per-base) transcript error
against 28.6 / 30.4 % for the corrected pseudocount, the synthetic pool at 128× / 61× — and 36.6 % even with the
simulator's own lengths for every hypothesis (the corrected pseudocount 23.5 %): at 98 % gDNA a 1 % error in the
level is half the RNA, and the pseudocount's damping is what absorbs it. On real data it would also need a rule
for the 194–288 unwitnessed loci per library and a `(1 − M)` multimapper term (`ISSUES:
unwitnessed-loci-and-multimappers-at-the-em`). Not rebuilt without both.

### the-unclipped-capture-efficiency
REFUSED 2026-09-21. Publishing `E[ρ / ρ_ref]` unclipped in place of `E[min(ρ / ρ_ref, 1)]` — so that pieces above
the reference keep their ratios — moves the per-base family from 3.0 % to 13.2 % from one scale at `g50 ss.99 ON`
(spans/isoforms 0.897, gDNA/isoforms 0.883): it inflates the isoforms' partial exons (1.35) more than the spans'
(1.21), because against the sampler's covering truth the unclipped estimator over-separates the exon classes
(fully probed 1.45, `½ < p < 0.9` 1.31, `p ≤ ½` 1.28) where the clip under-separates them (1.01 / 1.12 / 1.22). The
clip stays.
NOTE 2026-09-23: measured on the per-base rule, since refused (`ISSUES: the-per-base-rule-for-every-component`);
the one shared rule keeps the clip per object, and the unclipped form is re-priced on it before this is quoted
against it.

RE-EXAMINED 2026-09-24 and not re-opened: a reading that captured entities' contracted lengths follow only part
of the true capture range ("elasticity" 0.64 multi-exon at g50) came from weighting transcripts by the EM's own
gDNA exchange, which depends on the length ratio itself; with weights from outside the EM it reads 0.91–1.06
(`analysis/refute-mechanism/elasticity_exog.py` in `~/Downloads/rigel_runs/prototypes/2026-09-24_g98_dissection/`).
### nascent-siphons-gdna-under-capture
REPAIRED 2026-09-20 (`c52c9b93`). THE MECHANISM (per fragment against the read names' truth): the locus gDNA
component's opportunity counted a crossing start at EVERY boundary its fragment crossed, while its pseudocount
counted the fragment once — `assemble_priors` converted a boundary's incidence count by the accumulator's `q` and
left its incidence SUPPORT unconverted (`EQUATIONS.md` §11). Under capture the crossing support is 64 % of a
locus's opportunity (15 % off), so gDNA was priced at `ρ q̄`, handed its sense fragments at the probed exons to the
RNA hypotheses, and the unpinned synthetic shadows took them (95 % of their gDNA was exon-contained sense
fragments; the leak 17.1 % of a locus's gDNA at `q̄` 0.45–0.55, 0.5 % where crossings cross one boundary). One
E-step from the TRUE counts drifted +24 % per step on capture and was a fixed point off it. THE REPAIR converts the
crossing support by the same `q` (gates in `tests/calibration/test_priors.py`, watched to fire): `g50 ss.99 ON`
nascent +541,216 → +32,908, gDNA −626,550 → −56,319, transcript Σ|Δ| 5.21 → 5.50 % (the isoform ruler's own
over-statement unmasked, `ISSUES: ruler-witness-geometry-on-transcript-panels`); the deferred stratum's `g50`
siphon +1,567,550 → +569,674; `calibration_vs_oracle.py` and `policy_benchmark.py` unmoved. RULED OUT, with the
killing number: the ruler's witness geometry (shipped annotated/synthetic contraction ratio 0.0773 against the
simulator's 0.0738); the shadow-exclusive intronic pool (2.76 % of contained gDNA under capture); the per-locus
prior (structurally a gDNA:RNA lever, `--arm oracle` 541,216 → 541,762); nascent RNA itself (+522,205 at
`on_fraction` 0.10). A length knob on the shadows (doubling `L_n` removes 97 % and destroys the true nascent
signal) is a falsification probe, never a repair; the per-gene gDNA opportunity is refused
(`ISSUES: per-gene-gdna-opportunity`). The residual rides with `ISSUES: per-transcript-prior-lane`.
NOTE 2026-09-23: the length's `q` conversion is replaced by gDNA's own conserved share, which holds each start once
as this repair required (`ISSUES: the-pooled-q-in-the-gdna-length`).

### exclusive-evidence-support-rule
REFUSED 2026-09-20: an isoform supported iff its EXCLUSIVE spliced evidence carries reads (the per-transcript
prior's step ①), fed as a static weight through `rna_prior_weight` (the control reproducing base to 0.00 points),
recovers ~0 % of the allocation ceiling — transcript Σ|Δ| base → rule → true support → ceiling: `g50 OFF` 2.49 →
2.49 → 2.23 → 0.70, `g50 ON` 5.49 → 5.47 → 4.54 → 2.07, `g05 OFF` 1.59 → 1.59 → 1.55 → 0.43, `g05 ON` 3.49 → 3.47 →
2.83 → 1.08 (`ss 0.99`, fractional, `71037f6c`). The rule separates the testable set perfectly, but the EM's
likelihood already switches those isoforms off (122 false fragments of 9,950 at `g50 OFF`). Do not re-propose an
exclusive-evidence rule on junctions, regions or boundaries: the likelihood already uses exclusive evidence.

### a-per-object-gdna-rate-in-the-e-step
REFUSED 2026-09-20 (the calibrated-likelihood campaign's step 2, re-verified by the second adversarial review on a
corrected harness). Replacing the locus-pooled gDNA rate `count_g / length_g` in the E-step by calibration's
per-object gDNA density at the fragment's object gave no gain off capture and was worse on capture. One E-step from
the truth, the synthetic pool as a multiple of its truth: off capture every form reads 1.00–1.02; at `g50 ss.99 ON`
the locus-pooled fixed rate 1.14× (1.16× with the certified locus count), the per-object gDNA density 4.21×
(`local_gdna`), 1.18× even with the CERTIFIED per-object density (`local_true`), and 2.72× with a per-object RNA
rate beside it (`local_both`, 87 % of it a harness hard zero; the strike stands without it); at `g05 ss.99 ON`
1.04 / 1.25 / 0.95 / 1.18 in the same order.
A per-object RATE re-enters the gDNA-versus-span contest in calibration's density units against the ruler's; the
per-fragment gDNA PROBABILITY of `ISSUES: the-pseudocount-prior-is-biased-toward-gdna` is dimensionless and removes
that contest, and has not been measured.

### synthetic-spans-on-the-object-rule
REFUSED 2026-09-20. Pricing the synthetic spans by the gDNA component's object rule while the isoforms keep the
ruler puts the spans 18.7 % above the isoforms against the simulator's yield (`g50 ss.99 ON`) and the converged EM
collapses their pool to 0.55× its truth; a mixed state is the knife edge, never a repair.
NOTE 2026-09-23: the spans now share one rule over regions and boundaries with the isoforms and the gDNA
component (`ISSUES: the-gdna-component-length-rule-differs-from-the-transcripts`); what stays refused is a mixed state.

### per-gene-gdna-opportunity
REFUSED 2026-09-20, before it was built, on the measurement its premise fails. The candidate (2026-09-19)
rested on a toy in which every fragment of the locus sits inside the shadow's footprint, where the
shadow's growth factor is `L_g/L_n`; in a locus whose gDNA is uniform the footprint holds its share and
the factor is the density ratio, so pooling a locus's genes under one opportunity destabilises nothing by
itself. Measured on `g50 ss.99`: off capture `L_g/L_fam` reaches 5 with families UNDER-calling 7 %; on
capture families over-claim 3.5× where their own footprint's gDNA density per unit of the longest twin's
opportunity is exactly the locus's (D 0.8–1.2), and single-twin families leak 4.3× against 3.7× for
families of ten. A per-gene split of the shipped opportunity would inherit the twice-counted crossing per
gene (`ISSUES: nascent-siphons-gdna-under-capture`). Do not rebuild it.

### em-overturns-the-calibrated-gdna-split
CLOSED 2026-09-19, in two halves and by two different things. THE CAPTURE-OFF HALF CLOSED BY LANDING
`nascent-gets-no-rna-prior`: the EM's gDNA over-call there was the RNA prior's eligibility rule and not
the EM's arbitration — `g50 ss.50 OFF` reads a library gDNA fraction of 0.5127 against 0.50 where it
read 0.574, `g05 ss.50 OFF` 0.0538 where it read 0.104, and the nascent shortfall that mirrored it
closed with it (`g00 ss.50 OFF` nascent 1,901,090 → 2,018,540 against 2,024,341 true). THE CAPTURE-ON
HALF IS SUPERSEDED by `ISSUES: nascent-siphons-gdna-under-capture`, which names it correctly: the
missing gDNA does not go to the isoforms of probed genes, it goes to the SYNTHETIC NASCENT entities,
one fragment for one. The entry's own leading hypothesis — one gDNA rate spread uniformly along a locus
while capture concentrates gDNA at probed exons — survives intact and is carried there as a candidate;
what did not survive is the claim that the mass lands on annotated isoforms, which the pool rows refute
(annotated Δ is −2,752 / −6,560 / +20,767 at `g05` / `g50` / `g98` `ss.99 ON` against a nascent Δ of
+69,268 / +541,216 / +590,406). Do not re-open this name; the open thread is the siphon's.

### nascent-gets-no-rna-prior
CLOSED by landing 2026-09-19 (owner: "it's a hack; restore nascent RNA fairness"). The EM's RNA pseudocount
now goes to EVERY RNA component in proportion to the evidence it already carries, with none singled out for
zero; the eligibility test, the per-component `component_is_synthetic` flag, the `t_is_synthetic` lane and the
`index` parameter it was the only reader of are deleted (`EQUATIONS.md` §9b–§9b.2, `DESIGN.md` §0b).
THE WEIGHT — the one open choice — is `w_i = raw[i]`, the shipped weights, admitting every component at them
rather than an equal share (owner, 2026-09-19). Three reasons, each with its number: it keeps §9b.1's ABSORBING
STATE at zero evidence for free, since the state is a property of the weights and not of the eligibility test;
an equal share would be a strong informative prior rather than a neutral one, because `rna_prior_count` is a
conserved FRAGMENT COUNT (tens to thousands on an expressed locus), so `P/n` clears the ~0.16–0.47-alpha-unit
activation threshold by one to two orders of magnitude at every component; and it leaves
`component_rna_prior_weight` free for `ISSUES: per-transcript-prior-lane`, whose whole point is to fill that
lane with a MEASURED weight. THE PRICE, stated: the prior no longer helps a shadow entity decay — the rate
falls from `kappa/(1 + P/R)` per M-step to `kappa = w_N/w_T`, still strictly below 1 for free.
MEASURED. The defect was far larger than the entry recorded: on a tied two-component locus with one component
synthetic, a prior of 500 drove it from 100.0 fragments to 2.79e-298 and handed all 200 to the other. The
strict xfail this entry owned reads 70 → **0** fragments of leak onto the unexpressed antisense `t2`
(`test_nrna_multiexon_t2_low_ss`; the "72" in this entry and in `NEXT_SESSION.md` was stale). The gDNA:RNA
split does not move: the per-M-step identity `Σ out over RNA = rna_count + rna_prior` is gated unconditionally
in `tests/native/test_grouped_prior_update.py` and watched to fail under a perturbation that leaks the prior
across the gDNA boundary. A locus with no synthetic component is BIT-IDENTICAL — the operation order is
preserved for exactly that reason. Gates in `tests/test_estimator.py` (the prior moves no component's share of
the RNA pool, exact under MAP and within the derived digamma residual under VBEM; a shadow span still loses to
the transcript it shadows; a supported component keeps its mass; a zero-evidence component is not revived),
each watched to fail under its perturbation. ⚠ It also closed `nested-antisense-leak-under-the-sane-ruler`.

### nested-antisense-leak-under-the-sane-ruler
CLOSED by landing 2026-09-19, as a consequence of `nascent-gets-no-rna-prior` rather than by its own thread.
Both strict xfail rungs pass: the leak onto the unexpressed antisense `t2` falls 124 → **14** at SS 0.65 and
24 → **2** at SS 0.9, against the tight bounds of 20 and 5 the tests were written with. The entry blamed "the
EM's assignment at a nested transcript nothing witnesses" and ranked its lever as the per-transcript prior
lane; the diagnosis was one layer too deep. The leaked fragments were the HOST's nascent entity's, denied
their share of the locus's RNA pseudocount and landing on the one annotated transcript that could also explain
them — so the allocation rule was the mechanism, and the ruler was never implicated. ⚠ The residue at SS 0.65
is real: 14 of 2,000 fragments still reach a transcript whose truth is 0, at the strand specificity where the
separating channel is weakest. It is inside the bound and no longer an xfail; if it grows, it is the
unwitnessed-nested-transcript effect and `ISSUES: the-atom-at-an-unwitnessed-both-strand-slot` is its twin.
`tests/scenarios/test_antisense_intronic.py`, `quant_accuracy.py`.

### oracle-cache-key-hashes-a-thread-count
CLOSED by landing 2026-09-19 (owner): the scan settings' part of a cache's key is DERIVED at read time from the
settings the manifest records, and leaves out the two thread counts (`scan_cache._WORK_ONLY_FIELDS`), which divide
the work and not the tally (`test_scan_order_independence.py`, one worker to eight). No digest string is stored
any more, so the key's definition can move without stranding a cache: every oracle cache written 2026-08-22 at
`bgzf_threads` 4 loads again under the default config, with no re-scan. The five other resource settings
(`log_every`, `fragments_per_chunk`, `read_name_batch_size`, `buffer_size_bytes`, `spill_dir`) stay in the key:
their invariance is not gated. Gates in `test_scan_cache.py` (a thread count is not the key; a tally setting is; the
key is derived and no stored string decides), each watched to fail under its perturbation.

### oracle-ruler-arm-cannot-reach-the-ruler
CLOSED by landing 2026-09-19 (owner): `quant_accuracy.py --arm oracle_ruler` now hands the EM `fl × factor`, the
factor being the simulator's own capture-aware length (`ruler_vs_truth.load_truth`), anchored on the fully probed
class, cached per capture label beside the oracle caches; `oracle_ruler_noop` builds the same and hands back the
shipped lengths. The arm it replaced swapped calibration's count arrays, which the ruler stopped reading at
`c44fc306`, and refused every condition. Gates in `test_quant_accuracy.py` (the lengths reach the EM's published
`em_effective_length` and move the split, the noop reproduces `base`, the guard, the anchor), each watched to fail
under its perturbation — the anchor's first fixture could not tell the two pools apart and was repaired.

### end-to-end-error-unattributed
CLOSED by measurement 2026-09-19: the residual is the EM's and the assignment's. `quant_accuracy.py` on the
16-condition ladder at `cffca248`, seed 20260807, arms `base`, `base_reseed`, `noop`, `oracle`, `oracle_gdna`,
`oracle_rna`, `oracle_efflen` and `oracle_alloc_seed` (`~/Downloads/rigel_runs/arms/2026-09-19_e2e_baseline/`,
`decomposition.txt`). Per stratum, the `g05`–`g98` rows summed, in the order unstranded OFF / stranded OFF /
stranded ON and then the deferred one: transcript-level Σ|Δ| 520,117 / 453,135 / 1,833,492 / 4,487,120 fragments,
4.4 / 3.9 / 12.8 / 31.3 % of the true annotated RNA, against floors |base − base_reseed| of 3,478 / 5,802 / 3,532 /
5,079; gene level 246,441 / 210,165 / 601,452 / 2,640,220 (2.1 / 1.8 / 4.2 / 18.4 %), floors 558 / 282 / 2,228 /
3,605. The gDNA-free `g00` rows read 3.4 / 3.6 / 7.5 / 9.1 % at transcript level and `g05` 3.7 / 3.2 / 10.5 /
16.6 %, so most in-scope error exists with no gDNA at all. A PERFECT PRIOR (`oracle`) recovers 26,401 / 12,837 /
75,687 fragments at transcript level (5.1 / 2.8 / 4.1 %) and 24,913 / 11,497 / 74,687 at gene level (10.1 / 5.5 /
12.4 %); `oracle_rna` alone carries 21,697 / 14,074 / 73,855 of it, `oracle_gdna` 1,805 / 2,433 / 3,523, and
`oracle_efflen` −50 / +22 / −1,104 — ⛔ WITHDRAWN 2026-09-21: that arm overrides count arrays only and the
length reads none, so it cannot fire and those deltas are noise (`ISSUES: capture-blind-gdna-divisor`). By rung it is `g98`'s: at `g05` the move is below the floor
in all three strata (−1,694 against 2,431; +2,013 against 3,722; −1,373 against 1,942), at `g50` it is 2.0 / 2.8 /
1.2 %, at `g98` 37 / 25 / 35 %. Deferred: 40.9 % transcript, 70.4 % gene. So 95–97 % of the in-scope transcript
error and 88–95 % of the gene-level error survive a perfect prior, and what survives splits in two — the
per-transcript allocation (`ISSUES: per-transcript-prior-lane`: truth as the weights removes 49 / 46 / 61 %) and
the capture-OFF gDNA over-call (`ISSUES: em-overturns-the-calibrated-gdna-split`), which neither arm
moves. THE OLD CLAIM — in scope a perfect prior no longer improved the transcript number (2026-09-10:
`oracle`/`base` 1.019 / 1.026 / 0.978 transcript, 1.087 / 1.169 / 0.905 gene) — no longer holds as worded: today
0.949 / 0.972 / 0.959 and 0.899 / 0.945 / 0.876, each above its floor except at `g05`; which of the commits
between the two readings dissolved it is not attributed. Its substance holds and is now measured: the prior is
worth little in scope and the EM owns the rest. Two instruments were broken on the way in: `oracle_ruler` could
not fire (`ISSUES: oracle-ruler-arm-cannot-reach-the-ruler`), so the effective-length shrinkage sits inside no
ceiling here, and the oracle arms ran through a key-only wrapper (`ISSUES: oracle-cache-key-hashes-a-thread-count`).
The panel is not bit-reproducible at a pinned seed: `noop` matched `base` to ≤ 9 fragments per condition and two
direct `base` runs of one condition differed by 67, all below the floor.

### prior-fidelity-vs-deliverable
CLOSED by measurement 2026-09-19 (`ISSUES: end-to-end-error-unattributed`): the anti-correlation is gone on this
tree. A perfect prior now improves the deliverable in every in-scope stratum at both axes — `oracle`/`base` 0.949 /
0.972 / 0.959 at transcript level and 0.899 / 0.945 / 0.876 at gene level on unstranded OFF / stranded OFF /
stranded ON — where on 2026-09-10 it made both capture-OFF strata worse (1.019 / 1.026 transcript, 1.087 / 1.169
gene). The one in-scope row where it still reads worse, `g05 ss.99 OFF` (+2,013), is below its floor (3,722). The
leading answer (the messages destroying a nearly right self-solve at the worst slots) was never confirmed on a
second stratum and is no longer needed. `quant_accuracy.py`.

### scan-thread-split-starves-the-workers
CLOSED by landing 2026-09-19: the budget is split by the measured RATIO — one BGZF decompression thread keeps
about eight scan workers fed — instead of reserving a fixed four, so `bgzf_threads` defaults to deriving
`total // 8` and an explicit `--scan-bgzf-threads` still overrides. The old rule picked the worst measured cell
at every budget. Scan seconds on the 18.6M-fragment library, two interleaved rounds each, by (budget,
decompression threads): 4 — (0) 32.0, (1) 41.2; 8 — (1) 25.8, (2) 26.2, (0) 28.3, (4) 34.5; 16 — (2) 17.6,
(1) 19.6, (4) 20.1; the 2026-09-11 table it reproduces is (3,1) 113.5, (2,2) 59.5, (1,3) 42.4, (0,4) 33.7 at
total 4 and (4,4) 34.4, (2,6) 26.4, (1,7) 24.8 at total 8. BIT-IDENTICAL on all three identity references
including the real-library one that runs the scan, which is structural rather than lucky: every accumulator
bank is a sum of integers, and `test_scan_order_independence.py` holds the tally identical at 1, 2, 4 and 8
workers. `profiling/profiler.py --scan-only`.

### ruler-multimapper-floor-caps-the-correction
CLOSED by landing 2026-09-16 (`DESIGN.md` §7.2, `EQUATIONS.md` §11): the expectation ruler on the per-base
length — a transcript's own bases at their pieces' capture efficiencies with the fragment-length end taper, no
boundary or junction object in the length; a piece's efficiency the posterior mean of its clipped gDNA density
under the fitted landscape from its own contained count and the crossings at every boundary within a fragment's
reach, each crossing apportioned to the pieces its fragments cover by geometry times their own-count densities,
one pass — in `capture_efficiency` (computed by `calibrate`, published on the result), `capture_eff_length` and
`assemble_priors`; the multimapper floor `w = C/(C+1)`, the splice-junction objects and the flank imputation
deleted. The defect: the floor pulled every factor toward 1 by `1/(C+1)`, a +3.4 to +3.7 nat bias on the
unprobed class at every gDNA level; a piece shorter than a fragment had no contained support and read 0 then
the floor; forty junction objects imputed from unsupported flanks swamped four hundred bases of measurement on
the ladder. Killing numbers (`ruler_vs_truth.py`, the unprobed class's median log error; the test chromosome's
`g05 ss.99 / g25 ss.99 / g50 ss.50 / g50 ss.99 ON` rows): shipped +3.41 / +3.32 / +3.74 / +3.76 → landed +0.20 /
−0.03 / −0.15 / −0.03, the probed class 97–98 % → 99 % within ±0.1; the tiny-exon block's `capmixed` −0.63 →
+0.07 and `mixed` +5.9 → +0.2 (the unprobed `tiny` +7.0 → +0.7 / +1.9 at 5 % gDNA, where its edges hold 0.04
expected fragments and the posterior is the prior's valley); the depth ladder's probed class 75 / 93 / 97 % →
86 / 99 / 99 % at a tenth of the depth and 63 / 62 % → 91 / 91 % at a hundredth; the ladder's `g05 ss.99 ON`
unprobed class +4.40 → +0.45. Refused with numbers, each open question tried each way: the solve's own
posterior (+1.25 / +4.96 nat unprobed), the clip outside the expectation (indistinguishable), the own count
alone (+4.5), the plug-in (+3.36), the apportionment iterated to convergence (identical on the test
chromosome, +0.44 against +0.46 on the ladder), a joint update of neighbours (never converges). The
falsification tests were verified failing on the shipped ruler through its own fixture: an unprobed exon with
no gDNA fragment read the floor at +3.19 nat over the depleted level, a transcript of 40 bp exons read 0 and
then the floor. The multimapper blindness the floor stood for is `ISSUES: multimapper-blind-support`.
NOTE 2026-09-23: the per-base length and the crossing apportionment this landing built on are replaced
(`ISSUES: the-per-base-rule-for-every-component`, `ISSUES: the-crossing-apportionment`); the floor's deletion and
the posterior-mean efficiency stand. The own count alone, refused here at +4.5 nat because a piece shorter than a
fragment had bases to price and no count, now ships: under the one shared rule such a piece has no share, so its
efficiency multiplies nothing (`ISSUES: the-crossing-apportionment`).

### the-ruler-reference-on-sparse-real-libraries
CLOSED by landing 2026-09-16 (`DESIGN.md` §7.2, `EQUATIONS.md` §11): a basin's members are the kernels with a
location (count ≥ 1, `DensityLandscape.located`), the enriched candidate is the basin above the depleted one
with the most located members, and it is a mode iff more than `k = √n_located` members resolve it at the
population's k (the median width to the k-th nearest member ≤ 1 nat); the regime is on the result
(`gdna_reference_members`). The defect: `located_enriched_mode` tested membership on every kernel centre, and
a zero-count anchor's centre is its resolution wall `1/E` — a human index trains ~250–325 k anchors whose walls
span every decade, so on a sparse library a basin above the bulk was "located" by walls: LBX0190 (3,434 gDNA
fragments on regions, 1,118 located kernels) chose 10^-1.57/bp from 16,931 members, 16,918 of them anchors and
10 located, while its probed level (599 fragments in 3 regions near 10^-0.5) is no population; MO_3021
(65,874 / 15,088) chose 10^-1.35/bp from 20,053 anchors and 16 located members. Killing numbers, landed: both
panels' 46 rows and the depth ladder's 26 capture-ON rows unchanged to the reference (the enriched modes hold
257–3,606 located members); LBX0190 and MO_3021 read `None` (no basin above the bulk holds more than 11 / 16
located kernels against k = 33 / 123); LBX0588 (33,943 located, 11,579 members at 10^-0.53) and the VCaP
library (74,968, 24,496 members at 10^-1.07) unchanged; two capture-OFF rows of the depth ladder at a tenth of
the depth (`g00`, `g001`) that had read a reference from one located kernel among 37–43 walls read `None`,
restoring the capture-OFF contract there. The choice rule (largest rendered mass, most located members, largest
located weight, highest located basin) agrees on every row measured once the members are located, so the
most-located-members rule ships by derivation, not by number. Falsification test verified failing on the
shipped reader (a fitted basin of 600 short anchors' walls around twelve located kernels read a located mode of
612 members) and the perturbations fired: walls admitted as members, the k rule dropped, the members read at
their own √n_members. The three `review_identity_*` references re-frozen (LBX0190's reference is `None` now,
by design). The measured operating boundary stands: the reference is `None` below about 120 gDNA fragments on
the test chromosome, and about √n_train located probed pieces at one fragment or more are needed for a mode.

### capture-on-strand-pure-ambig-undercall
CLOSED by landing 2026-09-14 (`DESIGN.md` §6b.15.13, `EQUATIONS.md` §9f): the tilt atom — the AMBIG tilt's
hypothesis space {pure +, pure −, mixed} at equal reference weight, two columns at `τ = ±1` in ψ's cube beside
the continuum (`−log π` on its trapezoid weights), a held RNA level on a strand ruling the other strand's atom
out. The mechanism: the continuum's median sat below the strand cap at a strand-pure slot (prior-free, a truth
of 0.50 read 0.31–0.37; 0.46–0.50 now). Killing numbers, the L3 tree → landed: ladder stranded ON
470,862 → 427,046 (`g50 ss.99 ON` 217,636 → 199,409, `g98 ss.99 ON` 176,468 → 149,603), every stratum
better, 15 of 16 contaminated rows better, `g05 ss.99 ON` +1.7 %; stranded zero controls 497 → 550 / 224 → 231,
unstranded within a fragment; test chromosome −0.1 / +0.2 / −0.06 / +0.09 %; the strand-pure census band
−28 % / −38 % on `g50` / `g98 ss.99 ON`, its tilt error 40–80 % lower on every stranded row; the spliced stress
identical on every both-strand row; the encompassing region between TA+'s exons 0.366 → 0.511 against 0.544.
The source reproduces the prototype exactly on the test chromosome and to ≤ 0.007 % on the ladder. Gates in
`test_vertex_reference` (four perturbations fired). Not implemented: the structural witness (within 2 % of the
delivered-level witness). The cost where no witness exists is `the-atom-at-an-unwitnessed-both-strand-slot`.

### deadband-gates-a-gdna-free-library
CLOSED by landing 2026-09-14 (`DESIGN.md` §6b.15.12, `EQUATIONS.md` §5.2b): the strand channel is live iff the
protocol preserves strand — the Bayes factor on the spliced 2×2 of a free κ (the fit's own Beta(1, 1)) against
κ = ½ exactly, closed form, no constant (`region_init.strand_discriminability(κ̂, N)`); `disc = 4(κ̂−½)²` where
live; `n_gdna_obs` deleted throughout. Killing numbers: the replaced form (an unbiased estimate of `(κ−½)²`
floored at zero) is positive on 32 % of unstranded libraries — `g98 ss.50 OFF` shipped LIVE at z = 1.22, and
`g00 ss.50 OFF` was dead only by the gDNA term, without which 499 → 21,484; all four ladder `g00` rows had
`N_gdna = 0` and the channel dead at κ = 0.0099, so every g00-donor toy of the θ thread had run without a strand
channel. Landed: the ladder identical everywhere but unstranded OFF −0.22 % (`g98 ss.50 OFF` 123,657 →
122,981); the test chromosome identical on 30 rows; two contaminated rows bit-identical; the stranded zero
controls 405 → 497 / 194 → 224, located at the refit rung (the landscape's vertex bias,
`gdna-landscape-trains-on-false-positives` (d)); the encompassing exon∩exon slots on a gDNA-free donor
0.054 / 0.046 → 0.002 / 0.005; 17 goldens regenerated (≤ 2e-3 relative, `antisense_contained`'s false gDNA
78.7 → 5.6). REFUSED with it: an RNA level read from a slot's belief — the prior's answer re-emitted as data,
the relay (`TRAPS: one-hop-lifted-out-is-still-the-relay`). Gates in `test_region_init` (five perturbations
fired) and `test_encompassing_locus` on a gDNA-free donor (retired 2026-09-28).

### measured-prior-rung-4
SUPERSEDED 2026-09-14: ψ's composition reference fitted from composition-free observables (`rho_0`, the
enrichment responsibility, the detector as a boolean) was the measured-prior thread's fourth rung
(2026-08-26). The gDNA landscape hyperprior (`DESIGN.md` §7.1, trained on the previous solve's deconvolved
gDNA, the zero controls at a few hundred fragments) is the population prior ψ carries, and ψ has no reference
location (`DESIGN.md` §6b.1); the total-density landscape it would have consumed was QC-only, and is retired
(`ISSUES: the-total-density-landscape`). The requirements it listed (span both lattice ends, exact at `g00`, one
pseudo-fragment) are met or moot. Re-open only with a measured gap the landscape prior cannot close.

### reference-prior-refuted-at-concept-level
RECORD (2026-08-24; the ruling is `DESIGN.md` §6b.1): ψ's reference tilt was refuted at the concept level — on
λ the data's information `I ∝ N_eff·disc·[f_g(1−f_g)]²` is zero at κ = ½ while a tilt holds fixed nats there,
overturned at N = 3 where the strand channel is alive and never at κ = ½ (0.7471 from N = 10 to 10⁶); 82–95 %
of unstranded pass-0 error sat on slots it decided. Refused: another constant (0.75 optimal on none of 16),
`σ(L)` (3.3×/10.8× worse at the zero controls), re-weighting (up to 2.8× worse), an information-weighted tilt,
the arcsine coordinate (`ISSUES: arcsine-magnitude-coordinate`). The reference location is deleted.

### ambig-node-as-a-gdna-source
REFUSED as built 2026-09-08: a both-stranded node's own counts plus its held RNA levels, read as a gDNA level
and emitted lower-sided — `g05 ss.70 OFF` +2.6 % on all three test panels, +4…+38 % on the weak-κ zero
controls, ladder `g05 ss.50 OFF` 1.745× / `g05 ss.99 ON` 1.465× (at low gDNA the level's mode is noise and
travels as a floor); its target, the walled host exon of a `span` locus, gains ~100 fragments. Re-open only
with a gate on the emitted level's own width.

### psi-lambda-bracket-unshipped
CLOSED 2026-09-14, landed with the one-lattice ruling (W5, `DESIGN.md` §6b.15.10): the λ bracket follows the
landscape prior's derived demand (`DensityLandscape.required_logodds_window`, read by `calibrate._sweep`) and
the point count follows the bracket at the fixed step; nothing ships off.

### the-cancelling-pair
REFUSED twice (2026-08-26 the last, with the measured intron reference as load): `struct_lock` rescoped to
`g1_locked ∧ REGION` and the `intergenic|exon` boundary claiming its RNA-contaminated crossing mass as gDNA —
neither half prices alone (`TRAPS: a-cancelling-defect-pair`), worse on two of three in-scope strata with the
wins confined to `g00`. The one-sided certified-RNA bound is the only mechanism the zero control endorsed on
every row (−81.9 %, 8/8) and is panel-negative alone. Revive only with messages on at `g05 ss0.50 capture_on`;
the confusion matrix is re-derivable from `build_structural_claims` (deleted unread 2026-09-22; in git)
against `slot_truth.npz`.

### parked-capture-pilot-sign
MOOT 2026-09-14: two capture-ON pilot rows disagreed about the sign of every length correction (2026-08-13);
both panels were deleted and the correction lives in the length composition channel, retired until after
0.8.0. If the channel returns, find which row lied rather than averaging.

### strand-marginal-volume-factor
REFUSED 2026-09-13, both forms, with their numbers (the derivation and the arms: the sandbox's θ note §11). The
finding stands as a property, not a defect: the exact θ-marginal of the strand term at an interior tilt is
∝ σ_τ(λ) ∝ 1/(1 − f_g), an Occam factor of order √n toward gDNA at a balanced both-strand slot — under a
proper prior on the tilt it is intrinsic (any normalised prior gives it; "uniform in p" is uniform in τ and
keeps it), and it is what the strand channel genuinely says: balanced strands are explained by gDNA with no
parameter. Two ways to remove it were priced on the ladder (Σ|Δ gDNA| both axes, g00 excluded; reference
249,703 / 471,202 / 308,529 / 3,622,383 on stranded OFF / stranded ON / unstranded OFF / deferred):
* the tilt's Jeffreys volume, `+ log(1 − f_g)` on the AMBIG cube (the AMBIG reference Beta(½, 3/2)): 267,373 /
  761,711 / 331,064 / 5,263,361 — +7 / +62 / +7 / +45 %; the g00 rows 2–4 % better; the test chromosome +1–3 %
  on every stratum. Where: the gDNA-rich AMBIG slots under capture (g98 ss.99 ON strand-pure band 34,056 →
  120,845, g50 38,823 → 81,151) — where `a(λ)` is small the strand term does not constrain the tilt and the
  term is a prior shift toward RNA on slots that are mostly gDNA;
* the PROFILE form, the strand term entering the λ posterior through `max_τ L` (the cap, flat inside it), the
  tilt integrated only for the rows and the read-out: 253,117 / 687,783 / 308,557 / 3,653,341 — +1.4 / +46 /
  0 / +0.9 %; g00 marginally better; the test chromosome +0.2 / +0.8 / 0 / 0 %, the same g98 capture-ON bands.
On the mono shared-exon toy (no junction, nothing guards the slot) both forms do what they claim — a balanced
500k exon's own solve 0.997 → 0.49 — and on the spliced toy neither moves a number: the factor bites only
where nothing else informs the slot, and that slot's honest answer under a flat strand marginal is the
reference median under the cap, not zero. What survives: the width factor at gDNA-rich AMBIG slots is
evidence the panels reward, so the measure stays the arcsine and the marginal stays exact. The saturating
variant (normalising by the truncated volume) was refused on paper: it inverts the strand-pure tail. The
exposure — a deep, balanced, junction-free AMBIG exon in a gDNA-free library reads mostly gDNA — is recorded
here with the mono toy as its instrument (`deep_stress.py --mono`), and is the RNA level lanes' business (a
junction on either strand pins it), not ψ's.

### theta-quadrature-at-zero-gdna
CLOSED by landing 2026-09-13 (`DESIGN.md` §6b.15.11; the derivation `EQUATIONS.md` §9e): the θ nodes follow the
strand term's peak (`simplex_logodds._tilt_window`), the node count derived (24 = 2T/π + 1, T = −log ε₆₄). The
recorded mechanism was wrong: the K_t 30 failure was ONE slot (`g00 ss.99 ON`, slot 37345, n 25,242, 28 % of
its RNA on the minor strand, a 0.006 rad peak) — the fixed lattice's sum is a comb across λ on any deep
interior-tilt slot, and 60 held the control by landing a node on that one peak; the strand-purity quartic was
never it, and the trapezoid endpoint weights alone (REFUSED, ≤ 0.7 fragments) removed nothing because the
resolution was the whole term. Killing numbers: the marginal's λ-shape error 90–130 nats at 500k under the
lattice at 60 (0.1–0.4 at 2,400 nodes), 2·10⁻⁶ under the window at 24; the shared-exon stress at 500k, lattice
→ window: `g50` balanced 102,076 → 19 false fragments, 20 %-minor 11,448 → 11, `g00` 10,994 → 632 and
9,907 → 25, tilt error down 6–70×; the ladder a numeric near-no-op (stranded OFF / unstranded unchanged,
stranded ON +0.4 %, `g00` rows identical, the both-strand tilt error at 24 nodes what the lattice reached at
120), the test chromosome within ±3 fragments. Gates: `test_vertex_reference` — the marginal against adaptive
quadrature at 500 and 500k fragments, the derived count converged, a delivered row evaluated at the nodes to
the bit, κ = ½ the whole domain; five perturbations fired. The second step the same day removed the last θ
lattice: the lanes deliver a row's ingredients (`CubeRow`) and `sweep_n_tilt` is deleted — no tilt count
exists in the tool. What it exposed is `ISSUES: strand-marginal-volume-factor`.

### f32-strand-tilt-at-half
CLOSED by landing 2026-09-12: the AMBIG cube is float64 like the rest of ψ — the float32 cube was a memory
choice the tiling made moot, and one solver (`simplex_logodds._solve_logodds`, a single-strand slot the
``K_t = 1`` case) has one precision. At κ = ½ the strand mean is ½ identically in float64 and `w_pos` reads ½.

### landscape-trains-on-real-substrate
CLOSED 2026-09-10, superseded: the no-evidence share of the prior's training mass is now zero by construction
(`RegionBelief.has_composition`, gated); `landscape_training_census.py` reports the population per evidence class.

### the-landscape-training-population-arms
A/B'd on the test chromosome (30) and the ladder (16), halves apart, both zero controls on every row
(2026-09-10): two landed, six refused. Do not rebuild a discount that reaches the delivered rows, nor an
E-step on counted kernels. Whole-library |gDNA − truth|, ratio to the previous default; `noop` byte-identical.

| arm | what | ladder zero controls (ss.50 OFF / ON / ss.99 OFF / ON) | in scope | deferred | verdict |
|---|---|---|---|---|---|
| `no_echo` | a slot with no evidence, or a bound only, does not train | 0.742 / 0.753 / 0.962 / 0.966 | within 0.1 % on every row | 0.98–1.01 | ⭐ LANDED (`DESIGN.md` §7.1) |
| `no_echo_1sided` | also a one-sided delivered composition row | 0.512 / 0.660 / 0.748 / 0.910 | ≤ 1 % either way | ⛔ 1.90 / 1.28 / 1.20 | REFUSED: the probed exons' one-sided rows are the enriched mode's witness |
| `own_var` | the own-evidence variance (`tau_lam` + the row's precision) as `_reliability`'s `v` | 0.737 / 0.742 / 0.956 / 0.959 | ⛔ `g05 ss.50 OFF` 1.058; test chromosome `g05 ss.70 ON` 1.24, `g50 ss.70 ON` 1.12 | ⛔ 2.20 / 3.27 / 3.57 | REFUSED |
| `sweep1_no_echo` | `no_echo` at the first fit only | 0.890 / 0.892 / 0.994 / 0.995 | 1.000 | 1.00 | REFUSED: weaker than the rule at every refit |
| `shuffle` | `no_echo`'s excluded count on random non-anchor slots | 0.914 / 0.789 / 0.978 / 0.968 | 1.000–1.002 | 1.01–1.07 | the control: wins the zero rows (every non-anchor slot is false there) and nothing in scope |
| `refits6` | `calib_refit_iters` 3 → 6, no patch | 0.822 / 0.805 / 0.987 / 0.991 | 1.000; test chromosome unstranded ON 0.85–0.91 | 0.96–0.98 | NOT TAKEN: 0.725 with the rule against 0.742 without, at twice the refit cost (`ISSUES: performance-memory-bounded-solve`) |
| `row_kernel` (unit weight / reliability weight) | a delivered-row slot's own λ-likelihood, mapped to the density grid, as its kernel | test chromosome only: `g00 ss.50 ON` 437× / 0.95 | `g25 ss.50 OFF` 5.40 / 1.57; `g25 ss.50 ON` 2.10 / 0.95 | `g50 ss.50 ON` 1.62 / 1.00 | REFUTED: the sum of per-slot likelihoods is the flat-prior posterior, not a density; a bound spreads uniform mass under itself |
| `oracle_values` (weights as trained / var 0) | the CEILING: the population trained at certified gDNA | 0.686 / — / 0.956 / — (var 0: 0.749 / — / 1.113 / —) | `g05 ss.50 OFF` 1.016 / 1.061 | `g50 ss.50 ON` 0.964 / 0.753 | a price, not a mechanism: the values fix recovers most of it |
| `estep_zero` (ratios to the landed population fix) | the E-step on the kernels with no location (count < 1): kernel × the previous refit's landscape | **0.006 / 0.002 / 0.038 / 0.018** | unstranded OFF 0.952 / 0.973 / 0.988; stranded OFF 0.986 / 0.994 / 0.997; stranded ON 0.980 / 1.004 / 1.004 | 0.960 / 1.020 / 1.011 | ⭐ LANDED (`DESIGN.md` §7.1); the test chromosome wins or ties all 30 rows |
| `estep_all` | the E-step on every kernel | 0.006 / 0.002 / 0.038 / 0.018 | `g05 ss.99 OFF` 1.000, `g50 ss.99 ON` 1.018, `g98 ss.99 ON` 1.047; test chromosome `g25 ss.50 OFF` 1.091 | 0.956 / 1.045 / 1.012 | REFUSED: a counted minority competed away by the bulk, the δ-pin EM predecessor's failure |
| `estep_mirror` | the control: the previous landscape reversed on the grid | test chromosome capture-ON zero rows 1.48 / 12.7 / 8.9 | neutral on capture-OFF rows (the mirrored bulk lands in the depleted zone by the grid's asymmetry) | — | the control discriminates where it can: the placement is what acts |

### data-derived-reference-location
REFUSED 2026-08-24: the intron-reference estimator's selector widened to `solvable_exon` wins all 8
capture-OFF rungs and all four `g00` rows exactly and loses all 6 contaminated capture-ON rungs. Priced on all
16, deleted.

### flux-source-skipped-at-an-empty-exon-piece
CLOSED by landing 2026-09-09: `lanes.rna_lanes` builds a junction's flux level at an empty exon piece too,
priced by `hop_price` on the piece's zero count; ladder through the pipeline every in-scope row within 0.5 %
(worst `g98 ss.99 ON` 1.0048×), stranded zero controls 0.977×/0.958×. Only the ladder can judge it.

### the-empty-flux-source-at-the-junctions-counting-alone
REFUSED 2026-09-09: the sharper price for the empty-piece source (the junction's counting alone) wins the zero
controls more (0.944×/0.970×) and costs `g98 ss.99 ON` 1.0179× (+3,144 at the AMBIG exon|exon boundaries) and
`g50 ss.99 ON` 1.0091× — a sharp RNA floor reaching a node whose price reads the column count
(`ISSUES: flux-price-witness-units`).

### flux-witness-in-strand-units
REFUSED 2026-09-09: the column split's asymmetry as the flux price's witness reads worse at pass zero on the
test chromosome (`g05 ss.50 OFF` +12.6 %, `g05 ss.99 ON` +20 %) — it reads zero at an equal-abundance overlap
exon.

### relay-od-r-discontinuity
CLOSED with the deleted relay policy, 2026-09-09: `g98 ss0.50 capture-OFF` read 217,531 at `od_r ≤ 1e−7` and
212,581 at 1e−5, a threshold not a response (`TRAPS: a-constant-parked-a-value-off-a-knife-edge`). Standing
half: a 1e−5 nudge moved 1/30 rows more than 0.5 % (worst 1.65 %), so a single-row policy difference below
~2 % needs the noise floor re-recorded in the same session.

### levels-always-travel-for-the-gdna-lane
REFUSED 2026-09-08 (three forms): the gDNA lane on every directed face wins every capture-ON ladder row
(unstranded × ON −15 to −17 %) and loses `g05 ss.99 OFF` +7.1 % (45,076 → 48,295) and `g05 ss.50 OFF` +7.8 %
(50,948 → 54,926) — a one-sided floor at a node with no channel of its own is a tilt; with phase 2's ceiling
worse (`g25 ss.50 OFF` 1.59×, the weak-κ zero control 1.8–3.4×). Re-opened only by
`ISSUES: two-sided-exon-row`.

### the-edge-upper-side
REFUSED 2026-09-04, three upper sides on the intergenic|exon edge's level: undampened (`g98 ss.99 ON` +75 %,
`g98 ss.70 ON` +360 %); dampened both sides (sparse panel `g98 ss.99 ON` 36,645 vs 6,981); dampened above only
(no in-scope win, stranded 0/8, deferred 1.065–1.146×). No local witness prices the interior's enrichment.

### the-abundance-discrepancy-map
LANDED 2026-09-02, REPLACED by the level rule 2026-09-04: a fitted step between hypotheses is a pooled
premise; removing it improved `g50 ss.50 OFF` 12,328 → 11,044 and `g25 ss.50 ON` 28,414 → 19,854. Also
refused: a Gaussian summary of the held profile (12,328 → 13,911); the form built from what the boundary holds
(11,503). A level is made from the sender's measurement, never from what it holds; removing the dampening
costs 47 % on `g50 ss.99 ON`.

### the-certified-flux-row-as-a-level
REFUTED by probe placement 2026-09-04: a two-sided flux-implied composition row at certified faces wins the
exon-probed test chromosome (`g50 ss.50 OFF` 11,242 → 8,876) and fails where the probes move — junction-probed
`g50 ss.99 ON` 6,502 → 26,497 (spliced fragments enriched 20–35× more than contained ones), ladder zero
controls 3.3× (`g00 ss.50 OFF` 299,380 → 995,880) and 6.9×. A hard s ≥ 1 truncation refused alone
(5,712 vs 1,803).

### two-sided-exon-row-forms
REFUSED 2026-09-03 on the test chromosome: a two-sided Poisson row fixes capture-OFF (`g50 ss.50 OFF`
11,242 → 8,945) and is catastrophic capture-ON (`g50 ss.99 ON` 6,367 → 15,558); an abundance-bounded row still
fails `g50`/`g05` ON (12,333 / 3,090 against 6,367 / 2,293); a flux cap claims 19 % gDNA at the zero control
(16,819 → 121,263).

### the-discrepancy-priced-cap
REFUSED 2026-09-04: a cap on the face map's plateau priced by the pair's discrepancies fixes the first pass
(−20 to −24 %) and loses the junction-probed panel's weakly stranded rows (`g25 ss.70 ON` 1.122×).

### the-two-sided-level-lane
REFUSED 2026-09-05: with every node reached, a two-sided level lane reads −26 % at pass zero on
`g50 ss.50 OFF` and +33 % on `g50 ss.99 ON` — an exon complex's minimum total density bounds its gDNA from
above only off capture.

### the-wall-above-the-face-map-ceiling
REFUSED 2026-09-09: a wall above the face map's ceiling priced by the junction–exon pair closes the toy
harness gate (|Δf_g| 0.123/0.126 against 0.848) and wins the ladder's target rows at pass zero
(`g05 ss.50 OFF` 0.921×, `g50 ss.50 OFF` 0.847×), and through the pipeline loses the stranded half 0/6
(`g50 ss.99 ON` 1.142×), the deferred rows 1.15–1.43×, the junction and sparse panels' capture-ON stranded
rows 1.4–2.5×.

### the-pooled-hop-step
REFUSED by the owner 2026-09-03: a per-hop-kind pooled disagreement premise. The per-pair rule alone stands
within 0.25 % of it on every ladder row and ahead on four of six stranded rows (`g05 ON` −1.34 vs −1.38 %,
`g98 ON` −2.70 vs −2.49 %); where pooling wins (junction panel `g05` +3.5…+5.6 % vs +8…+9 %) is a systematic
offset the owner declines to extrapolate from other pairs.

### the-flux-factor-hop-premise
REFUSED 2026-09-02: fitting `1/delta` on the flux at the alternative splice site's E hop. Alone it fits
2.07 ± 0.26 on `g50 ss.99 ON` and harms (boundaries 231 → 280; sparse panel 394 → 562); jointly with the step
the two are degenerate at n ≈ 14 (sparse k = 2.59 ± 0.13, χ² = 47.7 for 14 pairs; boundaries 382 → 526).

### arcsine-magnitude-coordinate
REFUSED 2026-09-01: `f_g = sin²φ` in place of λ = logit(f_g). The premise is false (a posterior median cannot
reach a vertex in any coordinate; logit is finer there — top f-gap 1.83e-5 vs 1.37e-3 at K = 60), the leverage
is bounded at ~2.5 % of the vertex-band defect, and the ladder A/B loses every in-scope unstranded × OFF row
(1.107/1.071/1.014), winning only stranded × OFF (1.008/0.939/0.898). Commit `48be3d0a` holds the prototype.
Vertex information re-priced under the shipped policy (2026-09-09, `vertex_ceiling.py`): ≤ 1 % on stranded
in-scope rows, 2–7 % unstranded OFF, 25–35 % deferred.

### refused-transcript-weights
REFUSED (all 36 conditions, seed pinned): a soft-min per-transcript allocation over exclusive objects with a
per-object Jeffreys half — worse on every stratum and the zero control, transcript Σ|err| 57.5 M → 81.6 M
(1.42×; length-proportional 2.10×), the gDNA:RNA split +0.2 % so the allocation alone was priced. Exclusivity
hard-zeroed 38.7 % of transcripts; a density from < 200 bp amplified variance up to 6,534×; the +½ revived the
silent half (false-positive mass 18.6 M → 41.6 M). Deleted; the lane kept
(`ISSUES: per-transcript-prior-lane`).

### refused-soft-min-path-weighting
REFUSED: the owner's path theorem, 12 arms on `g00 ss0.99 capture_off` and 3 on `g50 ss0.50 capture_on` —
worse than base at transcript level on every rung (1.317–1.604× at `g00`, 1.262–1.331× blind), false-positive
mass 1.76–2.20× worse; the gene-level 0.395–0.527× at `g00` collapses to 1.006–1.041× blind. Mechanism
`TRAPS: an-upper-bound-is-not-an-estimate`: 3,644 of 4,839 silent transcripts (75.3 %) share an object with an
expressed one. Kept: the dial is monotone, 0.0 % of expressed transcripts zeroed. The instrument was retired
with the refusal (2026-09-10); git carries it.

### the-rna-length-law-fix
REFUTED 2026-08-31: `rna_pmf` is sound as shipped. The −4 bp sj-observability diagnosis was measured on the
undrained pool; on the drained payload the residual vs mature truth is −0.24/−0.11 bp off capture, −1.06 bp
under it, and dividing by `pi·o(w)` overcorrects +17 bp. The −17,754 transcript win from
`calibrate(rna_fl_pmf)` at `g05 ss.99 ON` is compensation (shifting the shipped shape +10 bp wins −171,342,
far past truth) on the EM's short-exon isoform misassignment; on the 0.8.0 metric the true pmf moves
±0.3–2.9 %, mixed sign. `r_hat` is accurate off capture (−0.04/+0.19 bp), broken under it (−11/−20 bp). Never
rank a calibration input on the thermometer.

### the-general-verdict-table
Mixed verdicts of 2026-08, each with its instrument. "The RNA fragment-length model" is the accumulator's fl
geometry; the length-channel retirement is of a calibration composition channel.

| | closed by | verdict |
|---|---|---|
| **the gDNA scale rule** · **the mass rescale** · **TSS/TES as the population licence** | landed 2026-08-04 | ✅ `EQUATIONS.md` §3.5/§3.5b/§3.5c (the gates named here retired with the relay policy, 2026-09-09). ⚠ The ceiling says the mass rescale cost the panel **nothing** (+0.0002 to delete it outright); it landed on the derivation and on being free |
| **face (I) of the `intron\|exon` BOUNDARY** | re-solve ceiling + panel arm | ⛔ **DO NOT BUILD.** The derivation (`EQUATIONS.md` §3.6) is re-verified and is not what failed: handing both BOUNDARIES the ORACLE truth and re-solving is worth **−0.000** off capture, and the ladder prototype is **negative** (mwae 0.0413 → 0.0426, confidently-wrong +10.7 %). TRAPS: panel-before-src |
| **a LEVEL transfer from the intron** | toy + panel | ⛔ **REFUTED**, +0.207 on capture-ON × unstranded — capture inverts which side is well-counted (TRAPS: capture-inverts-the-counted-side) |
| **the RNA fragment-length model** | a per-pmf ceiling arm (`length_ceiling.py`, retired 2026-09-10) | ⛔ **−0.02 %** at pass-0, **+0.21 % (worse)** over all objects. Root cause exact (`pi(w)` scores sj *crossing*, the pool requires the splice to be *seen*). ⭐ Its value is the BOUND: the whole fragment-length-model cluster costs ≤0.43 % of the shipped solve. TRAPS: price-the-halves-separately |
| **TRAPS: pure-and-length-censored's κ residue, as an ACCURACY fix** | κ injected at exactly ½, all 36 conditions | ⛔ **−0.2 %** unstranded, worse on the shipped solve. ⭐ But the *general* defect — a boolean licence flipped by a small residue — is **the-capture-level-residual**, and the destruction control taught TRAPS: honesty-metrics-reward-ignorance |
| **a nascent-bearing ladder condition** | toy, 36 conditions × 7 rungs | ⚠ **−5 %**, and the wrong way on one stratum. It was a toy harness arm (`--nrna 60`, retired 2026-09-28); it no longer justifies re-simulating the panel |
| **the gDNA prior's BIMODAL CAPACITY, and "give the prior more signal"** | a read of `gdna_landscape.py` + the production refit on real conditions | ⛔ **BOTH BRANCHES CLOSED.** The prior already renders the landscape correctly — **2.98 decades** of mode separation at `g75 ss0.99 capture_ON`, 30× more enriched mass ON than OFF, a single pile at the wall at `g00`. And a prior fitted from ORACLE truth is the same prior (0.04 dec). Not capacity, not signal, not location. ⭐ Why an evidence-free object cannot reach the vertex at all — and why that is the value of missing information rather than headroom — is `EQUATIONS.md` §9a |
| **the Jeffreys MEAN density location** | `--arm eta`, the `g00` zero control | ⛔ **REFUTED at +96,299 %.** It cannot say ZERO (`region_init.rho_g` is an exact 0 at 60,544/70,176 slots — the statement earning the −98 % at `g00`), and the TRAPS: a-ratio-cannot-carry-zero benefit it was credited with belongs to that fix, not to this arm. ⭐ If revisited the derived form is the Gamma **MODE** `max(a−½,0)/E`, which is exactly 0 at a zero count |
| **a threshold anywhere in the licence family** | TRAPS: a-threshold-on-a-fitted-residue implemented and refuted one | ⛔ τ is continuous across the region, so any floor is a tuned constant (TRAPS: a-threshold-on-a-fitted-residue, TRAPS: a-licence-with-no-floor, TRAPS: a-multiplication-gated-by-a-trace — refused three times) |
| **simulator captures a pre-mRNA through every probe its genomic blocks span** | NASCENT SCOPE RULING, 2026-08-22 | ✅ **WON'T-FIX** — nascent × capture fidelity is out of scope (`DESIGN.md` §0b); the sparse rebuild shrank the residual it explains (capture depletes nascent ~an order of magnitude), and no verdict depends on it. Re-open only if a real library forces it |

### the-reference-mean-family
Attempts to give ψ's reference a mean, plus two out of shipping it, each measured on the panel with both zero
controls and refused (2026-08-15/16); the surviving form is `DESIGN.md` §6b.1. Outside the table: giving `τ_λ`
the location term's curvature — bit-identical on the deliverable on all 32 rows, the 3,227× fall in `τ_λ` at a
pinned slot being ~98 % the `[f(1−f)]²` Jacobian (`TRAPS: a-priors-curvature-is-not-the-datas-information`); a
per-object one-pseudo-fragment floor `m_i = E[g]_i/(E[g]_i+1)` — worse on every stratum (0.609 / 1.045 / 0.580
/ 1.000 against 0.381 / 0.659 / 0.363 / 0.800), replaced by the lattice's top point `σ(L)` (`EQUATIONS.md`
§9c.1). A caution: `f_g ≤ 1 − S/M` from certified RNA is false (the truth violated it by 302) —
`boundary_spliced` is a separate bank, and the identity is `ρ_r·E_r = unspliced_RNA + S`.

| # | mechanism | why it was refused |
|---|---|---|
| 1 | **a fitted RNA density `logP_r`**, the mirror of the gDNA landscape | ⛔ the only non-circular form (fit from the solver's own belief) reads **0.988 / 0.997 / 1.037** — nothing, then worse: feeding ψ a density fitted from ψ's own belief tells it what it already believes. ⭐ And the ORACLE version's gain did not survive a shuffle — a shape that is wrong on purpose BEAT the true one at `g98` (0.786 vs 0.854), so the attribution was never established (`TRAPS: attribution-must-survive-a-shuffle`) |
| 2 | **a library-wide Beta mean, `a = f_lib`** | ⛔ `g05` regresses **1.43×**; `f_lib` is calibration's own output so the loop has positive feedback with both vertices attracting; and moving `a`/`b` sets the TAILS as well as the location, so `b = 0.03` leaves **57 %** of the prior outside `L = 10` |
| 3 | **the OBJECT-weighted mean instead of `f_lib`** | ⛔ **the two split by STRAND and neither wins everywhere**: object-weighted 0.584 / 0.452 on the two stranded strata but **5.570×** on unstranded × capture-OFF. ⚠ The sweep that motivated it was `ss_0.99` on three of its four rows — `TRAPS: never-pool-the-strata`, met on a sweep rather than a panel |
| 4 | **a stratified ASSERTION** — pure-gDNA strata claim `f_g = 1`, reweighted by stratum size | ⛔ reads 1.000 at every condition **including `g00`**: an assertion cannot see a library with no gDNA in it. ⭐ What replaced it is a per-object DENSITY, which needs no reweighting at all — the strata select the training set, not the answer |
| 5 | **a pooled RNA density from sj flux** | ⛔ RNA spans six decades with no genomic autocorrelation, so a pooled flux is not a population parameter (owner). It scored well only by sitting on the mass-weighted centre — `TRAPS: a-mean-hits-the-mass-weighted-centre-by-luck`. ⭐ Replaced by RNA-as-residual, which predicts no RNA at all |

### the-truncation-free-region-bank
REFUSED 2026-08-20: `region_start_count / ell` behind `CalibrationConfig.region_abundance_bank` improves the
pass-0 exon solve (no-evidence mass coverage 45.8 % → 66.1 % at `g98 ss0.99 ON`) and regresses the deliverable
— four zero controls 2.18×, deferred 1.84×, stranded × ON 1.20×, the only win 0.843× on stranded × OFF; under
the shipped bank 58.6 % of live exon hops carry a REGION bank of exactly 0, an accidental mute. Survives: the
truncation algebra (`EQUATIONS.md` §2) and the fill gate; the implementation is one commit before `a2b81b34`.
2026-09-27: the fill gate was deleted with the field it filled, `RegionGeometry.inv_abundance`.

### the-message-policy-campaign
Six mechanisms built, measured and refuted; closed 2026-08-27, the code deleted the same day. The bar it
missed is the one the transfer policy was later judged by — not to beat `SilentPolicy`, but to win on
unstranded data at minimal harm on stranded data, halves never pooled. The one measurement that survived:
propagation is net-harmful wherever the local solve has its own evidence, and its value concentrates where
that solve is blind.

| mechanism | what happened |
|---|---|
| **A composition-transporting policy** (`CurrencyPolicy`) | Best zero controls ever measured, but lost every in-scope contaminated stratum to silence. Deleted. |
| **The gDNA-continuity rule** (an unsupplied source's gDNA level crosses unscaled) | Built THREE ways — a static per-slot licence (a value RATCHET, gDNA densities to 3.9e+32, from breaking the knob's telescoping cancellation), a running-state licence (killed the ratchet, still lost capture-ON), and a fuse-based pure-gDNA re-anchor lattice (too weak — a fuse negotiates where the relay's mass rescale overwrites). Halved one row, lost others. Only safe beside a scan-time mass rescale. ⭐ The transfer policy's LEVEL LANE is this rule rebuilt as a one-sided profile with a priced hop (`DESIGN.md` §6b.12) — a different mechanism, judged separately. |
| **The premise's exon-end scoping** | The dispersion decomposition proves intron-end hops carry no COMPOSITION cost, yet scoping the charge to exon-end hops REGRESSED the panel: freeing intron chains before a measured LEVEL charge exists releases un-priced level drift. Restoring the pooled charge recovered `g98 ss.50 ON` 4.37 M → 2.44 M. |
| **A class-keyed method-of-moments fit on the observed log-ratio** (as the runtime law for the transport variance) | Tracks truth at intron and plain classes, REFUTED at sj classes: the route-summed flux cancels the visible step exactly where the true error is largest (0.08 observed vs 4.4 true). |
| **A totals-form pair fit** (for the same variance) | Refuted by construction — the knob consumes the totals, so transported totals agree by the (1−w) algebra and carry almost no information. |
| **`FanOutPolicy`** | Measured dominated once the certified-flux anchor gave its destinations own evidence. Deleted 2026-08-24. |

### the-doubt-graveyard
Eleven mechanisms priced and refused (2026-08-07). `g00` is the owner-required zero-gDNA control: its truth is
exactly 0, so every fragment there is a false positive with nothing to cancel it. The pattern is the result:
every one was a rule for resolving doubt, and at `g00` the doubt must resolve to no gDNA. The only candidate
the control has ever endorsed is the one-sided certified-RNA bound (−81.9 %, 8/8), panel-negative alone
(`ISSUES: the-cancelling-pair`); `zc_struct_lock_g1` is the one row still live, as half of that pair.

| candidate | what it did | `g00` | its target | why it died |
|---|---|---|---|---|
| `zc_jeffreys_mean` | `ρ_g = ½/E_g` at zero mass | ⛔ +7,269 % | −13.9 % | moves the mode UP |
| `zc_logmean` | `ρ_g = e^{ψ₀(½)}/E_g` | ⛔ +6,264 % | −11.3 % | moves the mode UP |
| `zc_anchor_mute` | no `prec_g` at empty locked slots | ⛔ +5,554 % | −7.7 % | kills the zero-gDNA win |
| `zc_struct_lock_g1` | scope `struct_lock` to `g1_locked ∧ REGION` | ⛔ +3,207 % | −1.2 % | ⭐ the MIS-SCOPED mask is load-bearing |
| `zc_reference_var` | `Var(f_g) = ⅛` where `τ = 0` | ✅ +0.0 % | −0.3 % | ⭐ passes the control and is INERT |
| `zc_discrepancy` | `+½ log D` shift, `(log D)²/12` | ⛔ +982 % | panel +4.5 % | moves the mode UP |
| `zc_disc_var` | the variance alone, mode untouched | ⛔ +255 % | panel +0.9 % | damping cannot bite |
| `zc_ref_prior` | own belief = ψ's reference, `τ + 1/π²` | ⛔ +3,792 % | −14.9 % | moves the mode UP |
| `zc_ref_prior_damp` | the two above, PAIRED | ⛔ +3,809 % | −15.5 % | ditto |
| the `eta` rebuild | a clean frame-free re-derivation | ⛔ unbounded | +85–103 % | see `DESIGN.md` §6.1 |
| the mean-location as a structural floor | the same idea in the LEVEL channel | ⛔ +96,299 % | — | it cannot say ZERO |
| `struct_lock = g1_locked ∧ REGION` **re-priced** | the standing strict xfail's own fix, on top of the SPLICE IN precision repair | ⭐ **0.71×** | in-scope **+2.5 / +2.1 / +0.5 %**, deferred +17 % | ⛔ **2026-08-18**, and the sign on `g00` FLIPPED from the 2026-08-11 row above once the relay was fixed underneath it — but the verdict did not: the four `nrna_none` zero controls go **1.8–15× WORSE**. The empty exons' "gDNA = 0 @ 0.2026" is LOAD-BEARING at AMBIG `exon\|exon` boundaries in an RNA-only library (an RNA+ level claim with no RNA− claim drifts ψ to 0.38 without it). ⭐ Still half of a pair |
| the mass rescale refuses a ZERO-MASS source (`pinM`) | `may_share_composition` additionally requires `M[src] > 0` | ⭐⭐ **0.31×** (134 k — better than the pre-licence relay's 154 k) | unstr × OFF 0.999, str × OFF 1.004, ⛔ **str × ON 1.032, six of six worse** | ⛔ **2026-08-18.** Under capture the empty slots between probe-covered stretches are the CONDUITS a relayed composition travels through (`TRAPS: the-divergence-was-a-barrier`), so the rescale at those hops is load-bearing. ⭐ The sharp predicate is NEITHER this nor the row above: refuse a source's OWN zero-count artefact, KEEP a relayed composition passing through an empty slot |

Four more `zc_*` arms exist as decomposition reverts used to attribute the 39 % win, not as proposals:
`zc_own_count`, `zc_live_count`, `zc_total_n` (inert) and `zc_transfer`, which reproduces the pre-fix tree.

### drain-provenance-split
REFUSED by arithmetic 2026-08-31: splitting the certified channel by drain provenance evicts 134,850 correct
drained records to remove 75 contaminants at `g50 ss.99 OFF`, resurrecting the −4 bp spliced-pool bias the
drain exists to repair (`ISSUES: drain-contaminates-certified-rna`).

### the-second-lambda-grid-and-its-regrid
CLOSED by landing 2026-09-13 (`DESIGN.md` §6b.15.10, W5 the grid study): one λ lattice for every consumer,
`sweep_logodds_step` 0.2; the fine single-strand grid, `_regrid_global` and `_scaled_grid` deleted. Priced on the
way and REFUSED, each with its number: a CUBIC regrid on λ (`g00` 892×, deferred 1.15 — the spline overshoots the
landscape prior's floor wall); a SMOOTH-RECONSTRUCTION read-out (a cubic or a monotone cubic of the log-posterior
on a 16× sub-grid: exact on Gaussians, 0.8 / 0.5 of a step on a balanced bimodal posterior against the histogram
quantile's 0.17, a wall overshoot on the cubic); a λ-axis LINEAR regrid (0.994–1.000 on the ladder: the axis was
never the harm, the interpolation was); a step DERIVED from the sharpest posterior (200–930 points where the metric
plateaus by 138, because an unresolved heavy slot costs ≤ n·f(1−f)·dλ/4 fragments and a fragment-budget rule needs
a tolerance); K_t below 60 (`ISSUES: theta-quadrature-at-zero-gdna`).

### rename-the-drain
RULED 2026-09-13 (owner): the term stays. Prepared and not pursued — the census counts (141 src / 146 scripts /
298 tests / 50 docs sites) and the candidates (`settle`, `decide`; `resolve` and `assign` collide) are in the
plan's §7 should the ruling ever be revisited.

### rename-row-and-face
RULED 2026-09-13 (owner): both terms stay. `row` (a slot's max-normalised log-profile over the solve grid;
1,236 src sites, two senses) and `face` (one directed side of a boundary, the `(destination, side)` pair a
rule is keyed by; 303 src sites) keep their names; the census and the candidates are in the plan's §7.

### g00-shrinkage-upstream-repair
CLOSED by landing 2026-09-14 (`DESIGN.md` §7.2, `EQUATIONS.md` §11): the ruler's reference is the located
enriched mode of the fitted gDNA landscape, published as `CalibrationResult.gdna_reference_density` and read
by `capture_eff_length` and `priors`; the kernel-density detector and its two constants are deleted. On both
panels (`calibration_vs_oracle.py` ③): the zero controls' factor 0.154 / 0.141 → 1.000 with nothing moved
(51,436 / 5,108 transcripts had moved); both in-scope capture-OFF strata P = O = 1.000 with nothing moved
(from P 0.957 / 0.970 and O 0.923 / 0.926 on the ladder, P 0.946 / 0.949 on the test chromosome — the
"exactly 1.000 off capture" contract line was not stale, the per-object clip of Poisson noise was the
defect); stranded capture-ON P/O 1.013 / 1.015 against 1.011 / 1.011, the reference within 4 % of the
truth's mass-weighted median on every capture-ON row measured; the deferred stratum 0.941 / 1.080 against
0.990 / 1.000 (one reference for both arms now, where two detectors had agreed by luck). The composition
was fixed first and the factor did not follow: 178–189 false fragments on 35,135 regions still made a
reference. The solve is untouched (`policy_benchmark.py` identical on both panels).

### u-ruler-arm
CLOSED 2026-09-14 with `g00-shrinkage-upstream-repair`: the "~2× loss on the two capture-OFF in-scope strata"
was the per-object clip, and with the reference a property of the solve the U ruler reads 1.000 by
construction (no reference) or `ρ̄/ρ_ref` against P's reference (0.06 on capture-ON, a number about
nothing). Its question — what a noise-free uniform field leaves — is answered structurally: nothing, the O
arm at capture-OFF reads 1.000 with no fitting. The arm is deleted from `calibration_vs_oracle.py`.

### flux-price-witness-units
CLOSED by landing 2026-09-14 (`DESIGN.md` §6b.13, `EQUATIONS.md` §12): the exon's witness of the junction→exon
price is its column count on the PROTOCOL'S SHARE of its RNA opportunity, `(c_s, κ_read · a_r)` — a column is
that much opportunity for the strand's RNA to be counted on it — so the junction's whole-strand rate and the
exon's density are in one unit and the pair's agreement is priced as counting alone at every κ (gated on a
hand-built chain at κ 0.99 / 0.7 / ½). The record: the golden `strand_ss65_multi_iso`'s gDNA-free exon reads
0.008 gDNA against 0.235. Both panels: the ladder in scope within 0.01 % on the metric and 0.01 % on the benchmark (stranded ON −0.01 %), the
deferred stratum +0.3 % / +0.07 %, the four zero rows identical (366 / 258 / 353 / 233); the test chromosome in scope within 0.2 %
(`policy_benchmark.py` −0.17 / 0 / −0.02 %), the deferred stratum −10 % on the benchmark and −8 % on the
metric, the zero controls identical. Three other forms priced and refused (`EQUATIONS.md` §12): the total
unspliced count (right unstranded, wrong at a both-stranded exon on stranded data — the ladder's
`g00 ss.99 ON` 233 → 435 on two AMBIG exons), the total as a one-sided bound (stranded capture-ON −8–10 %),
and the split's strand count at its own precision (the golden's ceiling 0.323, the refused
`flux-witness-in-strand-units` again).

### the-strand-overdispersion-reconcile
REPLACED 2026-10-01 by od = 0 (`DESIGN.md` §3.3a). `reconcile_overdispersions` blended a gDNA and an RNA od by
"information", and each input paired a value from one estimator with a precision from another:
- gDNA's precision came from a weighted fit that was computed, then discarded;
- RNA's was its null information, which credited it up to ~700,000× too much evidence (7.5×10⁷ against ~10²).

So RNA's value won: on `odg05` the shipped gDNA od was ≤ 0.0007 against a planted 0.05, and within one 1 % MO_3021 draw
the two values sat 0.19 apart. Against it, od = 0:
- LBX0588's gDNA total across 10 % subsamples, 3.45 % → 0.74 %;
- on the VCaP mix against read-name truth: gDNA pool −6.40 → −6.11 %, transcripts −0.15 %, genes −1.1 %.

### the-away-half-gdna-overdispersion
REFUSED 2026-10-01, with its influence weights, the `Beta(2,2)` 0.2 clip and the 0.2 no-evidence fallback. The pooled
moment over the genic seeds on the far side of ½ from their gene's RNA:
- **Antisense RNA reaches that side on every real library:** VCaP RNA-only 0.47 over 118,622 seeds, against no gDNA;
  MO_3021 0.133 against its junctions' 0.02.
- **"Unbiased under any RNA content" was false:** seeds mixing gDNA with sense RNA bias it low (planted 0.03 reads
  0.0264 ± 0.0021 against 0.0303 on pure seeds).
- **Its seeds cannot separate an od from composition at all:** at n ≤ 3 the two are one law (`EQUATIONS.md` §6a), and
  that is 73 % of LBX0588's seeds.

### the-joint-strand-overdispersion-fit
SUPERSEDED 2026-10-01, never landed. One influence-weighted estimating equation over the gDNA seeds and the junctions.
- **For it, on `odg05`:** calibration −10 % / −19 % and transcripts −3 % / −7 % (stranded OFF / ON).
- **Against it:**
  - it sat at the 0.2 clip on 13 of 14 real low-gDNA runs, whose gDNA seeds were antisense RNA;
  - it had two roots at low gDNA, with no rule between them;
  - its junction side carries splice artifacts (LBX0588's RNA moment is 0.07 at the pooled κ and 0.37 at the clean one).
- **Not built for it:** the contaminated-seed panel and the simulator's RNA-od arm.
- **Its designs and refusals:** `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/` and
  `ISSUES: the-overdispersion-design-refusals`.

### the-overdispersion-design-refusals
REFUSED 2026-09-28, in the derivations for `ISSUES: strand-overdispersion-one-shared-value`. Do not rebuild any of
these without the number that killed it.
ROOT AND WITNESS RULES:
- The two-witness random-effects rule (each seed set's smallest root, combined by a random-effects mean): the sign of
  a near-zero moment decides the junction witness, a hidden yes/no gate; a lone witness passes straight through (0.87
  once the gDNA-evidence fix empties VCaP's gDNA set); it costs efficiency when one shared ρ is clean (SD 0.0016 →
  0.0058 on LBX0588's depths at ρ = 0.05); it reads ≤ 0.003 on odg05 ss 0.99, where od 0 cost +19.8 %.
- min(x_g, x_r): 0.507 on VCaP at full depth, where both witnesses are inflated; the binomial value on odg05.
- The largest root: with the junctions alone it returns 0.2 on the 6 of 23 real runs whose least root is 0, and 0.2 is
  the most harmful value measured where the truth is 0.
- Maximum Q: Q is not a quasi-log-likelihood, and ΔQ on the two-root rows (−2.06, −2.60) sits inside its ±36 spread
  across draws.
- The walk from 0: its first step is the refused pair-count moment (one deep seed gives 0.83 against a fit of 0.017),
  and it overshoots in 4 of 1,200 low-count draws (0.2 where the least root is 0.031).
- The prototype's bisection is not a rule: it takes the upper root only because the unstable root (0.0063) lies below
  its midpoint.
- A posterior-mean od needs Q = ∫U to be a log-likelihood, which the data-dependent a_s rules out.
- Midpoint values for the fallback: ρ_eff at R = 0.2 is ×2.56, and at R = 1 (a uniform prior on the ICC) ×3.71,
  against ×1.07 for 0 on the minimal ψ.
- The sign-of-Ψ empty-evidence predicate moved 22 of 61 joint and 16 of 61 gDNA-only fits by 2–8 ulps, where the
  shipped-predicate form moves none.
PURITY WEIGHTS:
- A weight from the seed's own split ((2·min(A,T)/n)²): −66 % at ρ = 0.05; sample-splitting each seed, −28 %.
- The capture-blind own-count purity: truly pure intron seeds at odg05 g05 ON read π̂² 0.10 against 0.98.
- Transport with a Jeffreys ½: the ½ was an underived prior, and it decided every real-library row.
fl:
- Bisection to fl's fixed point: it finds a root, and on an empty pool its bracket is [0, 1000]; one refresh lands
  within 0.011 bp of the fixed point.
Sources: `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/` (`02_roots.md` §3, `03_bound.md` §3,
`04_evidence.md` §2–§4, `01_fallback.md` §2 and §5, `06_fl.md` §3).

### the-capture-reference-is-read-at-a-grid-point
CLOSED by deletion 2026-10-10: the reference density no longer exists (`DESIGN.md` §7.2); the reader's mode
is the lattice argmax refined by the parabola through its neighbours, per object, so there is no grid-point
reference to jitter. The open half of the record — a higher reference gave fewer transcript errors in 12 of 12
swaps, a compensating error downstream — is the junction price's structure
(`ISSUES: the-junction-price-is-noisy-within-a-gene`) and stays there.

### the-located-mode-capture-reader
CLOSED by replacement 2026-10-10 (`DESIGN.md` §7.2, `EQUATIONS.md` §11): the reader that read the fully captured
level off the landscape's located enriched mode (a basin above the depleted one with more than `√n` located
members) and switched every correction off when it found none (`None`) was an on/off capture decision, which
the owner refused ("an on-off switch for capture detection is not going to work", 2026-10-10). Do not rebuild
it, and do not rebuild its detector-free replacements, all measured end to end, pinned and fractional, against
the same shipped reader on the test chromosome (transcripts % / genes %):
- the mass-weighted median of the solved gDNA density as the reference: within 10 % of the located mode on
  every synthetic row that had one, 18.5× the located mode's answer on LBX0190; but on the zero-gDNA capture-OFF
  rows calibration's 366 false fragments (`ISSUES: the-gdna-prior-enters-psi-twice`) become the contrast and
  `g00 ss.99 OFF` reads 5.67 → 12.08; at `g05 OFF` +2 from Poisson noise; and where the probed class holds under
  half the gDNA mass (`g50`, single probe) the gDNA component is priced through the same efficiencies and genes
  read 1.48 → 16.82;
- 0.7.1's kernel-density peak (bandwidth 0.4, prominence 0.05, mass-weighted): the same rows, the same way,
  plus an outlier reference where one exon holds 7.7 % of a library's gDNA mass;
- the located-only and gDNA-majority variants of the median: no change at `g00` — the false gDNA IS located,
  and at `g00` exactly one object is gDNA-majority, and one object defines a median;
- every cheap readout of the solve's inferred count under the landscape (mode, mean, median; 3–6 s a genome):
  equal to the honest likelihood's mode on every panel except the zero-gDNA rows, where `g00 ss.99 ON` reads
  7.35 → 14.65 through the same false gDNA.
The honest likelihood with the RNA amount integrated out is immune to the false gDNA because it reads the
strand columns, not the inferred count; its per-object error is at objects with no read under a median or mean
readout (−0.14 to −0.79 nats, the breadth of the smoothed prior) and under 0.04 nats under the mode wherever an
object has a read. Its native kernel reproduces the NumPy prototype to 1.2e-9 on 70,176 ladder slots.
