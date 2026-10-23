# The observable genome — the dissection and the design

2026-10-22 · code base `acc029bb` · revision 2, for the owner's review. Revision 1 dissected the fair panel's loss and
proposed counting N runs out of every opportunity. The owner's direction since (2026-10-22): mappability underlies every
reading of "observable" — reads align neither to N nor to repeats an aligner cannot place — so incorporate it through
`alignable`, simply, as a mask for 0.8.0; the universally unobservable regions first, what depends on the aligner and the
read length later. This revision is that design. Sandbox: the measured dissection lives in `ISSUES:
region-lengths-count-unobservable-bases`; nothing here is settled, and no production code is written.

## 0. In brief

- **The observable genome** is the set of bases a read pair can be observed at: placed by the aligner where it came
  from. The index records its complement as a MASK of unobservable intervals, and every opportunity in the tool counts
  only the start positions whose fragment avoids the mask. Nothing else changes: not the scan, not a deposit, not a count.
- **For 0.8.0 the mask is the universally unobservable set** (owner, 2026-10-22): with an `alignable` store, as real
  genomes are built, every base no tiled read pair covering it was placed back at, at the store's longest read length
  (`alignable`'s own `unmappable_regions(threshold=0.0)`, so no threshold is chosen) — N included, since no read maps to
  N; without a store, the reference's N runs, the one unobservability the tool can know with no aligner's word. On
  GRCh38 the store's set is 164.5 Mb: 151.06 Mb of N and 13.4 Mb of repeats, mostly intergenic.
- **It fixes the measured bug.** The fair panel's index has no store, so its mask is its 18.3 Mb of N, which is exactly
  what the prototype counted out: the length law recovers, and with it the gDNA (0.898 → 1.000 off capture, 0.788 → 0.992
  on) and the transcript table (2.82 → 2.63 %, 6.09 → 4.28 %). It also puts every other density in the tool —
  calibration's, the priors', the EM's lengths — on observable bases.
- **What it leaves, measured**: bases an aligner places only part of the time (the rung's pericentromeric DNA, about
  half mappable), the read-length dependence (≤ 7 Mb on GRCh38), and a second-pass defect the rung exposed that is not
  observability at all (§7).

## 1. The dissection (revision 1, condensed)

The A/B pair, g50 κ0.99 capture-off and capture-on, mask on against off with the same seeds, pinned and fractional, each
EM unit against its read name (receipts: `$RIGEL_RUNS/prototypes/2026-10-22_region_length/`).

- **The term is the gDNA length law the scorer reads.** `gdna_density.one_sided_rate` (`gdna_density.py:119-172`) takes
  its exposure over each region's whole length (`fl.py:379-381`) and trusts the low side, because contamination can only
  add (`gdna_density.py:9-13`). Regions that are almost all N — chr21's first, 5,011,798 bp of which 1,798 are known,
  holding 87 fragments; chr22's 10.7 Mb one at 2.1 % known — read as the cleanest gDNA, and the root falls to 0.0015
  against a true 0.063 on observable bases (0.025 of pooled off capture, 0.17 on). The purities `aₚ = min(ρEₚ/nₚ, 1)`
  (`fl.py:388-392`) fall to 0.037 and 0.018; the contrast (`fl.py:421`) divides by their 0.0195 separation and amplifies
  each length bin's noise about fifty-fold, which the resolution weight cannot see (its noise is the purities' scalar
  error, `fl.py:399-405`). The realized law the scorer reads (`pipeline.py:992-993`, `scoring.cpp:882`) is off by KL 0.53
  nats per fragment off capture and 0.70 on (0.00005 and 0.0053 unmasked), and pruned 5,443 true gDNA fragments' gDNA
  candidate outright (`scoring.cpp:238-239`).
- **Calibration does not lose it.** Its count in loci is right (2.957 M against a true 2.951 M off capture), its per-object
  error moves under a point in every node class, and the EM loses the same share of every locus (median 0.84, IQR
  0.81–0.88), with or without nascent RNA. The perfect-prior arm never reaches the law (`_oracle.py:369-404` replaces
  masses only).
- **The fix, measured**: the law's region inputs over non-N bases alone recover the EM end to end, as injecting the
  simulator's law does.

| condition | arm | gDNA | nascent | transcript Σ\|Δ\| | gene Σ\|Δ\| | realized law, KL |
|---|---|---:|---:|---:|---:|---:|
| capture-off | mask off | 1.000 | 1.00× | 2.64 % | 0.28 % | 0.00005 |
| capture-off | mask on, shipped | 0.898 | 1.51× | 2.82 % | 0.50 % | 0.532 |
| capture-off | mask on, observable law | 1.000 | 1.00× | 2.63 % | 0.30 % | 0.00022 |
| capture-on | mask off | 0.992 | 0.67× | 4.38 % | 2.06 % | 0.0053 |
| capture-on | mask on, shipped | 0.788 | 7.67× | 6.09 % | 3.11 % | 0.701 |
| capture-on | mask on, observable law | 0.992 | 0.66× | 4.28 % | 2.08 % | 0.0066 |

- **Its reach**: every gDNA condition of the ladder and both fl-gap panels (the realized law shipped at KL 0.12–0.91,
  0.0002–0.024 with non-N lengths); the STAR rung (g50 κ0.99: gDNA 0.795 → 0.981 off capture and 0.743 → 0.970 on,
  transcript 13.42 → 9.25 % on), most of `ISSUES: the-gdna-length-contrast-is-ill-conditioned-on-aligned-input`; and
  the real VCaP mix only lightly (its DNA is exome-captured: the rate reads 0.32 of pooled; the realized law against its
  DNA library's own lengths, KL 0.100 → 0.084).

## 2. The observable genome

### 2.1 What it is

An opportunity counts the start positions at which a fragment of length `w` could sit (`EQUATIONS.md` §2;
`effective_length.py:72-90`). It must count only the starts at which the fragment, had it been there, would be in the
BAM at that place: a start nobody can observe is not a place where a component is absent but no place at all, and a
divisor that counts it reads every density low by its share — the fair panel's bug, and a standing one on real data.

Whether a start is observable depends on the reference (N), the aligner and its settings (its multimapper limit above
all), the read length and the fragment length. `alignable` measures exactly this per base: it tiles the genome with
read pairs, aligns them, and counts, for every base, the tiled pairs whose span covers it (`n_frag`) and those the
aligner placed back where they came from (`n_frag_ontarget`), per read-length bin. A pair is on target when both mates
have any reported hit at their true place (`alignable/scorer.py:318-331`), so a multimapper whose hits include its true
place counts — which is Rigel's notion, since Rigel keeps multimappers and draws among their places (`DESIGN.md` §4.2).
Its BINARY mappability `n_frag_ontarget / n_frag` is the share of a base's covering pairs that are observed; the
FRACTIONAL `n_frag_ontarget / nh_ontarget` is about one over the copy number and answers a different question (how
many places share a pair), not observability. The old index read the fractional one as a per-region exposure (April
2026, `git show e8e7ad7f:src/rigel/mappability.py`); it was lost in the calibration rewrite, never refused by a
measurement.

### 2.2 Three layers, in the order to take them

1. **Universal (0.8.0).** The reference's N runs, and the bases of ZERO binary mappability at the store's longest
   read-length bin. A longer read maps more, so the zero sets nearly nest: on GRCh38, 415 kb of the 125 bp set (0.25 %)
   is not zero at 100 bp and 366 kb not at 50 bp, so the longest bin's set is unobservable at every shorter read length
   but for that sliver. `alignable`'s own rule (`AlignableStore.unmappable_regions(..., threshold=0.0)`: no on-target
   pair, which includes no pair at all), so there is no threshold to choose.
2. **The library's read length (later).** The zero set at the library's own bin is larger for shorter reads.
3. **Partial mappability (later).** A base observed some of the time: weights in place of a mask.

### 2.3 What 0.8.0 builds, and the rule for a fragment

The mask is the store's zero set when the index has a store, and the FASTA's N runs when it has none (owner,
2026-10-22: N is part of mappability, not a mask of its own). A fragment is observable iff its span avoids the mask — the rule
the simulator already draws by for N (`sampler.py:243-264`), and the one `alignable`'s coverage counts imply: a base is
zero only when no tiled pair whose span covers it was placed, so a pair covering it is not observed. A region's
opportunity is then the sum over its observable pieces, `Σ_pieces E_f[(ℓₚ − w + 1)₊]`, and gDNA's reach at a boundary
is the distance to the nearest masked base or reference end, where it is `UNBOUNDED_REACH` today
(`region_geometry.py:192-199`). The observable bases per region, per reference and per genome are the pieces' lengths
summed, recorded beside the region and reference lengths.

`alignable`'s zero set already treats "a gap two mates span" physically: a short N run that tiled pairs span with both
mates on known sequence keeps on-target coverage (67 kb of GRCh38's 151.12 Mb of N stays unmasked at 125 bp), so with a
store the tool counts those spanning pairs as observable, as a real aligner places them. Without a store the N runs mask
whole templates, the simulator's rule (`sampler.py:243-264`), so the panels stay fair.

### 2.4 Measured on GRCh38 (the local STAR store, `alignable_grch38_star.zarr.zip`)

Built with STAR 2.7.11b at bins 50, 75, 100 and 125 bp, fragments 250 ± 50 bp, a 1 % error rate, five pairs per
position, and the arriba STAR settings with a 50-hit limit — the ladder's STAR rung's own settings less the two the
MANUAL recommends (`sjdbScore 0`, `outFilterMultimapScoreRange 2`). Over the 3.10 Gb genome (`grch38_mask.py`):

| bin | masked | N in it | beyond N | intervals | median interval | in introns | in exons | partly mappable (0 < m ≤ 0.75) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 50 bp | 171.5 Mb | 151.07 Mb | 20.4 Mb | 28.8 k | 485 bp | 1.61 Mb | 104 kb | 3.3 % of bases |
| 100 bp | 166.5 Mb | 151.06 Mb | 15.5 Mb | 31.2 k | 345 bp | 1.37 Mb | 94 kb | 1.7 % |
| 125 bp | 164.5 Mb | 151.06 Mb | 13.4 Mb | 30.4 k | 296 bp | 1.32 Mb | 91 kb | 1.5 % |

The mask is N plus 13–20 Mb of repeats (163 Mb of it intergenic at 125 bp, 12 % of intergenic sequence); read length
moves it by 2 Mb between 100 and 125 bp and 7 Mb down to 50. Most bases sit at binary mappability above 0.75 (the
store's 1 % error rate keeps many just under 1); 1.5–3.3 % are partly mappable, which a mask leaves at full weight.
**Real reads agree with it**: of the VCaP mix's 37.1 M primary records, 636 put aligned bases inside the 125 bp mask,
582 of them the DNA library's 150 bp reads (longer than the mask's bin) and 626 of them multimappers. On VCaP the
repeat part changes the length law by nothing measurable beside the N part (KL 0.0844 → 0.0849); on GRCh38 the N runs
do the work, and the repeats matter where they sit inside loci.

## 3. Where it is recorded — the index

`rigel index` reads the FASTA's lengths from its `.fai` and never its sequence (`index.py:219-233`, called at `1064`);
with `--alignable-zarr` it reads the store's splice blacklist alone (`index.py:1070-1077`; the MANUAL: "The store's
per-base mappability is not read"). Proposed:

- **Without a store, the N runs** from the FASTA's sequence, by one function shared with the simulator, whose
  `unknown_runs` (`sampler.py:64-76`) moves down to a module both import, so the tool and the panel read one rule.
- **With a store, its zero set at its longest bin** (N included), read through `alignable`'s API in windows (one bin of chr1 read
  whole is 4 GB) and stitched at window edges; about five minutes for GRCh38 (four bins took that long here). Where
  the store and the FASTA share a reference their lengths must agree, or the build is refused (another assembly); a
  reference the store lacks takes its N runs, named in a warning, as the blacklist's unmatched references are
  (`index.py:160-191`).
- **Written**: `unobservable.feather` (`ref, start, end`, merged and sorted; 31 k rows on GRCh38); an
  `observable_length` column on the regions and on the reference lengths; in the manifest, the mask's sources (the
  FASTA; the store with its aligner, version, bin, fragment model and error rate) and its totals.
  `INDEX_FORMAT_VERSION` (`index.py:85`) goes to 11, and the load refuses an index without the mask
  (`index.py:1287-1294`), with no fallback.
- **`--no-mappability`** builds the mask from the N runs alone: synthetic genomes, and the ladder's oracle rung, whose
  simulator places every fragment by truth, repeats included. The flag's name is then true again (`ISSUES:
  hygiene-ledger` (f)).
- **Nothing the scan reads changes**: `partition_hash` and `graph_hash` cover the partition and the junctions, not the
  mask (`index.py:926-990`), so every scan cache stands.

## 4. Where it is read

**One object, at load**: the observable geometry, derived from the mask and the region partition into the frame the
scanner and calibration share (`region_arrays.RegionArrays`, `region_arrays.py:69-129`): each region's observable
pieces and length, each boundary's distance to the nearest masked base either side, each reference's pieces. Every
site below reads it instead of a raw length; on an index with an empty mask every one of them returns today's number to
the bit, so the test chromosome and every synthetic index are fixed points.

| stage | site | what reads a length today | step |
|---|---|---|---|
| frame | `splice_graph.py:814-860`; `region_arrays.py:129` | the bounds; `region_size_bp = end − start` | 2 |
| opportunity | `effective_length.py:72-90` | contained support over the region's length | 2 |
| opportunity | `effective_length.py:93-151`, `region_geometry.py:192-199` | a boundary's crossing support, gDNA at `UNBOUNDED_REACH` | 2 |
| opportunity | `effective_length.py:211-239` | a cut's conserved shares from the pieces' lengths | 2 |
| opportunity | `gdna_density.py:175-193`, `196-216` | contained support; region lengths per reference | 2 |
| opportunity | `gdna_opportunity.py:199-245` (`218`, `230`, `236-244`) | the four pools' `A(w)` and `T(w)` | 2 |
| accumulator | `second_pass.py:465-472` | a placement's density, count over its region's support | 3 |
| length law | `pipeline.py:487-515` → `fl.py:379-381`, `388-392`, `407-415` | the rate's exposure, the purities, the contained pools' de-tilts | 4 |
| length law | `fl.py:516-526`, `531-534`, `573-574` | the realized law's rate, each region's RNA excess, the boundary odds | 4 |
| calibration | `region_geometry.py:190-210` | every slot's divisor, `eff_gdna` and `eff_rna` | 5 |
| calibration | `region_geometry.py:264-284`; `landscape.py:109-135`, `319`; `calibrate.py:128-219` | the landscape's grid and kernels; the zero-count anchor (`174-179`), which an all-N region is | 5 |
| calibration | `calibrate.py:555`, `579-589`; `capture_efficiency.py:131-196` | the published divisors; each object's support in the capture reader | 5 |
| calibration | `calibrate.py:464-474`, `598` | gDNA's conserved share at each boundary | 5 |
| before the EM | `priors.py:331-337` | the locus gDNA length `Σ S_r c_r + Σ M_e c_e` | 5 |
| before the EM | `pipeline.py:660-672`; `capture_eff_length.py:96-170`, `217-255` | every transcript's fl-marginal and EM lengths; a synthetic span's is its genomic span | 6 |
| QC | `derive.py`; `cli.py:515` | `gdna_density_global`, read by no stage | 7 |
| instrument | `calibration_oracle.py`, the gdna-field gate | expects `ρ·E_g` on the tool's own divisor | 5 |

`priors.py:105-113` shares a region among loci by base overlap (a share, not an opportunity), and `sj_opportunity.py`
reads transcripts' exon lengths; neither needs the mask for 0.8.0 (exons hold 91 kb of it on GRCh38). The EM reads the
law and the lengths (`em_solver.cpp:2048-2058`, `scoring.cpp:882`), so it changes through steps 4–6 and needs no change of
its own.

## 5. Why the counts may stay as they are

The mask changes denominators only. That is safe because, by its construction, no read pair is placed back inside it at
the store's read lengths: on the panels no fragment is drawn over N, and on VCaP 636 of 37.1 M records put aligned bases
there (§2.4), almost all multimappers or reads longer than the mask's bin. Withholding a hit inside the mask would be a
scan change with nothing measured to gain; it is not proposed.

## 6. The steps, in the order the tool runs

Each prototyped outside the tree and A/B'd against what ships before it enters `src/`; a falsification test on a
hand-built case, verified failing on main; the fix alone; each gate broken on purpose and seen to fire; then the panels
per stratum, pinned and fractional, against truth, with every capture efficiency at 1 beside the shipped run while the
reader stands, the oracle ladder for regressions, the suite and `preflight.py --full`.

**Step 1 — the index records the mask.**
- *1a, a no-op first*: `tests/_index_builder.py` writes an all-N FASTA for every structural test (`_index_builder.py:76-83`),
  which a mask would turn into zero opportunity everywhere; it writes known sequence instead, proven byte-identical on
  every gate and golden while nothing reads the bases. The simulator's `unknown_runs` moves to the shared module.
- *1b*: the mask, the observable lengths, the manifest, format 11.
- *Tests* (`tests/test_index_observable.py`): from a hand-written FASTA — a run at 0, a single N, a lower-case run, a run
  ending at a reference's end, a reference with none, an IUPAC `R` that is not unknown — and a store of known counts
  written by `alignable`'s own writer (zero bases inside a region, across a window edge, a short N gap pairs span), the
  index records exactly what a per-base enumeration finds — the N runs without a store, the store's zero set with one,
  and N runs only on a reference the store lacks — through a reload; the observable lengths are the enumeration's; a
  store whose references disagree in length with the FASTA is refused; a format-10 index is refused.
- *Breaks*: drop the lower case; end a run one short; read the store at a shorter bin; lose a run at a window edge.
- *Measured*: nothing moves — `rename_identity.py --check` against the frozen references; the suite at its count plus
  the new gates; `preflight.py --full`.

**Step 2 — the observable geometry and its opportunity functions, read by nothing yet.**
- *Tests*: the brute-force enumerations in `tests/calibration/test_effective_length.py` and `test_gdna_opportunity.py`,
  extended to masked references (a run inside a region, at its edge, filling it, two closer than a fragment, beside a
  boundary): every start of every length is placed, refused if its span touches the mask, and classified by the deposit
  rule; each closed form — contained support, gDNA's crossing support with its reaches, the cut shares, the pools'
  `A(w)` and `T(w)` — equals the counts exactly. On main the shipped functions, handed a masked region, return its
  whole-length count and fail.
- *Breaks*: keep the masked bases; leave a boundary's reach unbounded; take `T(w)` over the reference length.
- *Measured*: a numeric no-op.

**Step 3 — the accumulator: the second pass reads a placement's observable density** (`second_pass.py:465-472`).
- *Test*: two places of one multimapper at the same density, one in a region holding a masked run; the shipped draw
  weighs that one low and fails.
- *Measured*: the STAR rung, where the multimappers are; the oracle ladder for regressions. Expected small.

**Step 4 — the length law reads observable opportunity: the measured fix.** `library_fl_models` hands `fl.py` the
observable geometry instead of raw region lengths, so every read of a length there reads one object.
- *Tests*: a hand-built payload of intergenic and intronic regions at a known gDNA density, RNA added to the introns,
  one region all masked and one half masked: the one-sided rate recovers the density on observable bases within its own
  error (`Var(log ρ) ≈ 1/Σg`, `EQUATIONS.md` §7.1); with a pure intergenic pool `a₀` clips and the contrast reduces to
  that pool's de-tilted histogram (`fl.py:346-356`), moved toward the pools' mixture only by the resolution weight
  (`fl.py:426`). On main the rate collapses and both fail.
- *Breaks*: each read of a length reverted alone.
- *Measured, and the bars*: the length-law table on every panel condition (the rate over pooled, the purities, the law
  against the simulator's), shipped against fixed; the A/B pair at the prototype's numbers (§1's table); the whole
  ladder, both fl-gap panels and the STAR rung, every in-scope stratum no worse than shipped on transcripts, genes and
  pools, and expected back near the panel before the re-simulation (2.36 / 2.02 / 3.88 % transcript) less what its other
  changes moved (+0.26 % off capture and −0.10 % on at g50 κ0.99); `calibration_vs_oracle.py`, because calibration's
  divisors read the gDNA law (`region_geometry.py:209`); the VCaP mix on an index built with the GRCh38 store (its gDNA
  fraction against 0.2518, its law against the DNA library's lengths: 0.100 → 0.085 by the prototype).

**Step 5 — calibration's geometry reads it.** The slot divisors (`region_geometry.py:190-210`), and through them the
solve, the landscape (an all-masked region has no support and stops training as an anchor), the capture reader's
supports, the published lengths and the priors' `S_r`; gDNA's conserved shares at boundaries (`calibrate.py:464-474`).
- *Tests*: a region holding a masked run takes its observable support as its divisor; an all-masked region is no anchor;
  a boundary beside a masked run takes the bounded crossing support.
- *The panel's own falsification check*: the oracle's gdna-field gate fails today on every gDNA capture-off row by
  exactly its reference's N share (chr21 0.85 and 1.15–1.16, chr22 0.71 and 1.27–1.29) and passes on the unmasked arm
  (0.99–1.00); it must pass on the masked panel after this step, and fail again with any of the step's reads reverted.
- *Measured*: `calibration_vs_oracle.py` per stratum; `quant_accuracy.py` on the ladder, both fl-gap panels and the STAR
  rung, efficiencies at 1 beside shipped; the reader's depleted objects off capture back to the unmasked arm's (117
  regions and 5 boundaries below 0.9 now, 63 and 0 unmasked).

**Step 6 — the RNA side before the EM**: every transcript's length over its template's observable pieces
(`pipeline.py:660-672`, `capture_eff_length.py:96-170`); the simulator masks RNA templates at their N too
(`sampler.py:219-228`).
- *Test*: a transcript and a synthetic span whose templates hold a masked run take the enumeration of their observable
  starts.
- *Measured*: `quant_accuracy.py`'s nascent pool and transcript table, `ruler_vs_truth.py`; on the ladder one span holds
  N (chr22:18,428,561-18,486,998, 50 kb of 58 kb, whose locus reads 0.936 of its gDNA off capture with the law fixed).

**Step 7 — the report**: `summary.json` carries the observable bases per reference and genome and the mask's sources;
`gdna_density_global` reads fragments per observable base; the MANUAL's mappability section says what the store now
gives. *Measured*: the summary schema's gates; no number moves but the density.

**The STAR rung's mask.** The oracle rung keeps the N-only suite index. The rung needs its aligner's mask for chr21 and
chr22 alone — not GRCh38's, which counts copies elsewhere in the genome (its chr21/chr22 zero set beyond N is 49 kb at
100 bp). Proposed: a store for the suite genome built on the cluster with STAR 2.7.11b, the rung's parameter file
(`scripts/sim/configs/star_recommended.params.txt`), 100 bp reads and the ladder's fragment model, and a second suite
index built with it for the rung. rustar cannot stand in: `alignable` drives the `STAR` command with options rustar
does not take.

## 7. What the mask leaves — measured on the STAR rung

The rung holds the same simulated reads as the oracle BAM, so comparing each region's count across the two shows where
the aligned path loses fragments: intergenic contained counts read 0.967 of the oracle's (−67.8 k), introns 0.999. Two
causes, neither a zero-mappability repeat:

1. **Partly mappable DNA.** chr21:10.6–13.0 Mb and chr22:12.6–15.3 Mb (pericentromeric) lose a quarter of their gDNA
   unplaced — 34.9 k of 141.6 k and 34.4 k of 146.6 k pairs, together 69 k of the rung's 71.6 k unplaced gDNA. In the
   GRCh38 store their binary mappability is 0.55 and 0.52, zero on 5–12 % of bases: a mask catches the zeros, not the
   half. That is layer 3, weights in place of a mask, for later.
2. **The second pass empties duplicated copies.** At chr21:7,092,716-7,129,102 the aligner reported 2,192 of 2,197 gDNA
   pairs with a hit at home among 2–9 places, and the drained tally holds 2 there; an intron at
   chr21:34,459,938-34,481,525 goes from 1,254 to 0 the same way. Both are observable (binary mappability 0.999); their
   fractional mappability, 0.21 and 0.50, says several copies share them. A region all of whose fragments are
   multimappers has no unique pass-one traffic, and the draw weighs a place by that traffic (`second_pass.py:400-472`),
   so it is drawn rarely or never: the region form of the "all-held approximation" `ISSUES: the-second-pass-corrections`
   records for junctions. About 107 k contained deposits move between regions this way off capture at g50, and an
   emptied copy reads to the one-sided rate exactly as an N run did. It belongs to that entry, measured there; the mask
   must not paper over it.

## 8. For the owner

**Decided (owner, 2026-10-22).** (1) The universal mask for 0.8.0, with N taken through mappability: the store's zero set
when there is a store, the N runs only when there is none. (2) The existing GRCh38 store serves for now (bins 50–125 bp,
STAR); a new store, and `alignable` driving rustar as well as STAR, come later, on the cluster. (3) The ladder keeps its
store-less suite index: its mask is its N runs, which is the measured fix. The GRCh38 store covers all 94 of the suite
genome's references at the same lengths, so it could be used there, but on chr21 and chr22 it masks only 29 kb beyond
N at 125 bp, its mappability was computed against the whole genome rather than the rung's two chromosomes, and it would
bring its splice blacklist into every ladder scan; the rung's own store comes with the later `alignable` work.

The questions as they were put:

1. **0.8.0's mask**: the universal layer — N, and the store's zero set at its longest bin. Recommended.
2. **The store.** The local GRCh38 store's longest bin is 125 bp and its STAR settings lack the MANUAL's two; is it the
   0.8.0 reference store, or is it rebuilt with the recommended settings and a 150 bp bin? Expected, not measured: the
   zero set barely moves, since the 50-hit limit decides it; the splice blacklist the index already reads from the store
   would move with a rebuild.
3. **The STAR rung's store**: a suite-genome store built on the cluster, and a rung index of its own (§6).
4. **Later, each its own design**: the library's bin (layer 2), partial mappability as weights (layer 3), and a check
   that a BAM's aligner settings match the store's (STAR writes its command line in the BAM header).
5. **The second pass's emptied copies** (§7.2) are recorded, as measured, in `ISSUES: the-second-pass-corrections`;
   whether they join that entry's arms for 0.8.0 is yours.

## 9. Built and measured (2026-10-22, after the owner's "Go")

Steps 1–4 are in the working tree, uncommitted: one patch per step under `$RIGEL_RUNS/arms/2026-10-22_observable/patches/`,
applying in order on `acc029bb`, with the receipts beside them (suite logs, identity logs, preflight logs, the
calibration rows, the VCaP runs). Every arm below is pinned and fractional, scored against truth; the harness is
`$RIGEL_RUNS/prototypes/2026-10-22_region_length/` (`run_arm.py --sp-observable | --fl-observable-src | --eff1`,
`sweep_score.py`).

| step | what | suite | verified |
|---|---|---|---|
| 1 | the index records the mask (`unobservable.feather`, format 11, `rigel.observable`); the test builder writes known sequence | 3,393 (+18) | the rule against per-base enumeration; 11 breaks fire; bit-identical on oracle g00, STAR g50u, cfRNA-A and -C |
| 2 | the observable geometry (`region_arrays`) and the opportunity twins (`effective_length`, `gdna_opportunity`), read by nothing | 3,408 (+15) | against the reference accumulator and enumeration on masked references; 15 breaks fire; bit-identical ×4 |
| 3 | the second pass's placement density over the region's observable support | 3,410 (+2) | the masked-copy test fails on main (0.502 = 151/301) and passes; 2 breaks fire; oracle g00 bit-identical |
| 4 | the length law over observable pieces, and the pools' `A(w)`, `T(w)` over the observable genome | 3,414 (+4) | 4 gates fail on main and pass; 6 breaks fire, each read reverted alone; src bit-identical to the hook |

**Decisions taken in the build, for review.** `score_held_fragments` takes the frame (`regions=`) in place of the bare
region types, which it derives (they are inert in the scorer, whose accumulators read only `length_under` and
`sj_edge_ids`). `fl.py` takes the frame too, so whole lengths, pieces and later the reaches come from one object; its
contained supports are `gdna_density.contained_opportunity` summed over pieces — not step 2's
`observable_contained_eff_length`, which renormalises its pmf and would break the empty-mask fixed point.
`gdna_opportunity_from_index` and `region_lengths_from_partition` are no longer called by production and are kept until the
tree is committed (deletions after a commit). The census law's exon on-target shape and its boundary odds still read
whole lengths — step 4b below.

**Step 3 is not small** (the design expected it to be): on the STAR rung about 1.03 M of 2.09 M contained-case
placements sit in a region holding N off capture, 63 k on; the self-check found 0 changed values outside a masked
region.

| STAR rung, g50 | transcripts | genes | gDNA | nascent |
|---|---|---|---|---|
| unstranded off | 5.39 → 5.29 % | 2.05 → 1.96 % | 0.695 → 0.746 | 2.44 → 2.19× |
| stranded off | 5.88 → 6.12 % | 2.42 → 2.35 % | 0.795 → 0.845 | 1.94 → 1.69× |
| stranded on | 13.42 → 12.75 % | 7.63 → 7.38 % | 0.743 → 0.752 | 9.22 → 8.86× |

Stranded off's transcript loss is one BCR isoform flip (12,332 of the 16,310 that got worse; BCR's gene total moves
477). g00 stranded off: +30 fragments. The oracle BAM holds no multimapper, so step 3 cannot move the oracle ladder.

**Step 4, the sweep** (the ladder on both BAM kinds, both fl-gap panels; 80 pairs, `step4/sweep_scores.tsv`). gDNA of
truth, the rows with gDNA, shipped → step 4:

| stratum | oracle and fl-gap | STAR rung |
|---|---|---|
| unstranded × OFF | 0.817–1.026 → 0.968–1.002 | 0.725–0.853 → 0.959–0.980 |
| stranded × OFF | 0.889–1.017 → 0.995–1.001 | 0.835–0.917 → 0.977–0.983 |
| stranded × ON | 0.788–1.005 → 0.992–1.003 | 0.752–0.877 → 0.957–0.988 |
| unstranded × ON (deferred) | 0.580–1.003 → 0.803–1.001 | 0.488–0.616 → 0.739–0.812 |

Transcript error at g50, unstranded off / stranded off / stranded on: oracle 2.73 → 2.45 / 2.82 → 2.62 / 6.09 → 4.19 %,
STAR rung 5.29 → 4.95 / 6.12 → 5.91 / 12.75 → 9.08 %; at g98 the oracle reads 28.9 → 17.1 / 26.9 → 15.1 / 102.6 →
31.6 % and the rung 43.3 → 27.9 / 36.5 → 32.7 / 128.7 → 87.5 %. Every fl-gap row improves (stranded on 5.35 → 5.15 %
RNA-long, 14.70 → 13.50 % RNA-short). Calibration's library fraction moves under 0.001 on every in-scope row and its
per-object errors a few percent either way (`calibration_vs_oracle.py`, `cvo_step3.json` against `cvo_step4.json`).
On the VCaP mix (a real-data test input), on the index built with the GRCh38 store, the gDNA fraction reads 0.2574 →
0.2526 against 0.2518 by read name, and the realized gDNA law reads KL 0.101 → 0.086 nats against the DNA library's own
lengths, its mean 293.6 → 306.8 bp against 307.4 (the uniform law 0.166 → 0.168, as the prototype found).

**What step 4 costs, each traced.**
- Oracle g05 stranded off, transcripts 1.57 → 1.62 %: one gene's isoforms (ADARB1, +3,567); genes improve (−364).
- Oracle g05 stranded on, transcripts 3.13 → 3.16 %, genes 0.63 → 0.65 %: diffuse; nascent RNA reads 1.00× → 0.88× and
  gDNA 0.927 → 0.995, so some nascent RNA now reads as gDNA under capture (it persists with every efficiency at 1).
- Oracle g50 stranded on, nascent RNA 7.7× → 0.68×, annotated 1.012 → 1.017: the gDNA nascent RNA absorbed hid a
  nascent → annotated leak of about 50 k fragments under capture. Undissected; `ISSUES:
  the-scorer-reads-a-census-length-law` is the nearest frame.
- STAR rung g00 unstranded off, transcripts 29.47 → 31.16 %: isoform churn (481 k worse, 347 k better) on the capture
  reader's false-gDNA regime; with every efficiency at 1 the same pair reads 3.02 → 3.02 % (`ISSUES:
  the-reader-reads-false-gdna-as-capture`).
- STAR rung g98 stranded off, annotated 1.036 → 0.874 (0.856 with efficiencies at 1), while its transcripts improve 36.5 →
  32.7 % and gDNA 0.835 → 0.985.

**Next.** 4b: the census law's boundary odds and exon on-target shape over reaches (16 of the suite's 22,230
exon-flanking boundaries have a masked base within 1 kb, and no exon holds one; on GRCh38's store mask 1,455 boundaries
have one within 1 kb), A/B'd alone. Then step 5 (calibration's geometry), 6 (the RNA lengths) and 7 (the report), as §6
has them.
