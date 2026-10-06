# A capture ruler with no reference: the RNA-short fix, prototyped (2026-10-05)

*Sandbox document (`docs/dev/`): provisional, not authoritative, cited by nothing outside the sandbox. Written for the
owner and for reviewers. Every number is pinned (one thread everywhere, fractional assignment, EM seed 0), in
process, one scan and one calibration per condition with every arm's EM on the same replayed buffer
(`harness.py --verify` reproduces the CLI to the fragment). The prototype, its harness and every table:
`~/Downloads/rigel_runs/prototypes/2026-10-05_spectrum_ruler/` (`arms.py`, arm `pwfx_own_med`; README there). Nothing
in `src/` implements it.*

## 1. The question and the answer

The RNA-short gap arm's stranded × capture-OFF library reads transcripts at 42.3 % on the tree against 6.31 % on
0.7.1. The owner's ruling of 2026-10-05 forbids the obvious repair: the ruler's `None` ("no enriched mode, so contract
nothing") is a binary capture detector, and capture is a spectrum.

The prototype ruler reads that row at **3.55 % (genes 0.23 %)**, the figure of a ruler that contracts nothing at all
(3.53 %), with no detector: its weights come out equal because the data say so. It beats 0.7.1 on every stratum of the
test chromosome under all three probe layouts, and on the sparse plasma library it finds the capture the tree's mode
reader misses. **It is not a release candidate**: on the mirror arm (gDNA short, RNA
long) its unstranded × capture-OFF row regresses from 2.13 % to 11.71 % against 0.7.1's 5.29 %. The cause is the same
one, met from the other side: calibration's per-object gDNA counts are wrong under a fragment-length gap, in both
directions, and coherently enough that no rule in the ruler can tell them from capture (§2). The shipped ruler survives
both rows only because its `None` ignores every per-object count on an uncaptured library. **The robust fix is in
calibration**: per-object gDNA counts that are right under a length gap. The ruler built here is what makes the tool
spectrum-native once they are (§6).

## 2. Why every ruler that trusts per-object counts fails these rows

Calibration's per-object gDNA counts on the gap arms are grossly wrong (off by more than 3-fold and 5 Poisson sd) at
hundreds to thousands of objects, while the library totals are right. Against the oracle's per-object truth
(`slot_truth.npz`):

| capture-OFF row | regions over-called | regions under-called | boundaries over-called |
|---|---|---|---|
| RNA-short, stranded | 173 (median 7.3× the truth, +2.4k fragments) | 156 (−16k) | 1 |
| RNA-long, unstranded | 136 (7.7×, +31k) | 0 | 1,431 (6.8×, +55k) |

On the RNA-short row the over-calls are short exons that cannot hold a gDNA fragment (under 1 bp of opportunity; their
strand noise divided by that opportunity), and the under-calls are intronic clusters whose gDNA was booked as RNA. On
the RNA-long unstranded row whole RNA-rich neighbourhoods are booked as gDNA, the exon and both its boundaries together
(one exon: 2,071 gDNA fragments of its 3,037, truth 47): coherent, so it looks exactly like captured gDNA.

This is the count-frame deconvolution under a length gap (`FRAGMENT_LENGTH_POSTMORTEM.md` §2). The shipped ruler turns
the first group into a false capture reference; its `None` path hides all three on uncaptured libraries, which is why
"no reference" reads 3.53 %. Every per-object estimator that takes these counts at face value inherits them: relative
posterior means of the landscape read the row at 10–34 %, and a fitted population with a Poisson likelihood at 33.8 %.
A noise model from the solver's own variance (`Var(log f_g)`) does not rescue it: that variance is sharpened by the
prior at the false exons and far too wide where the prior and the messages pinned an accurate count (exon|exon
boundaries between probed exons at ss 0.70: count 39 against a truth of 39, variance up to 5 nats²), and that arm
collapses captured rows (RNA-short stranded × ON 56 %).

## 3. The mechanism, in five steps, with no constant chosen

1. **Each object's level is its posterior median under a population of gDNA densities.** The population is the
   nonparametric maximum-likelihood mixture on the landscape's own grid, fitted to every region and boundary with a
   Poisson likelihood of the solver's count on its gDNA opportunity, and **each object votes by its opportunity**: it
   is the distribution of density per base the genome offers gDNA, so an exon that cannot hold a fragment carries no
   weight in it. It is deconvolved: where every object is consistent with one level it puts its mass on one atom,
   which a kernel sum cannot. The median ignores a level a handful of objects hold; the mean follows it, and so does
   the geometric mean for a level near zero (unmappable regions). Plasma's small probed population is an atom like any
   other: the share it holds is estimated, never tested.
2. **A floor at the off-target density, because capture only enriches.** The level is the intergenic regions' pooled
   count over their pooled opportunity: they admit no RNA, so there is nothing to deconvolve.
3. **Coherence, because capture is spatially coherent at the scale of a fragment and a calibration error at one object
   is not.** An EXONIC region takes the opportunity-weighted median of its own level and its two boundaries' — a short
   exon is read through the gDNA that crosses it, a long one keeps its own. A non-exonic region keeps its own level,
   held at or below its boundaries' larger one where its own opportunity is smaller than theirs. Every boundary lies
   between its two flanks' levels. The exon/non-exon split carries one assumption: probes sit on exons, so an intron's
   contained fragments are never captured. Without it, short introns between probed exons take the probed level
   (1,414 of them on the RNA-short captured row).
4. **The weights are relative.** The EM's weight per component is `log θ − log L`, so a factor common to every length in
   a locus cancels: no reference level, no clip, no `None`.
5. **A junction is the sum by conservation of bases, capped at the larger of its two pieces' levels**: the owner's
   "capture saturates" (2026-10-05), read locally instead of against a global 1.

What it replaces: `landscape.located_enriched_mode` and its census, the reference, the clip at the reference, the
`None` path, and the landscape posterior inside `capture_efficiency`. What it keeps: calibration's composition, the
landscape in ψ, the one shared conserved-share rule, the EM.

## 4. What each step was measured to do (offline, against the oracle's truth)

| step removed | what breaks |
|---|---|
| the opportunity vote | false atoms from the short exons; RNA-short OFF 33.8 % end to end |
| the median (mean instead) | the zero-gDNA rows: one intergenic region with 8,748 unannotated-RNA fragments makes a far atom (25 %) |
| the floor | 1,545 under-called boundaries and 324 regions stay low; the coherence step then spreads them |
| coherence | probed short exons read the uncaptured level (RNA-short stranded × ON 49.8 %); one false boundary at 7,000× on the zero-gDNA captured row |
| the exon/non-exon split | 1,414 short introns at ~230× on the RNA-short captured row |

## 5. The ledger, transcripts / genes (%), per stratum

The tree is `main` with the junction cap at 1 (landed 2026-10-05, uncommitted). 0.7.1 from the receipted Phase 0 sweep.

| panel · stratum | 0.7.1 | tree (cap at 1) | prototype |
|---|---|---|---|
| RNA-short · stranded × OFF (the defect) | 6.31 / 0.44 | 42.30 / 3.28 | **3.55 / 0.23** |
| RNA-short · stranded × ON | 20.33 / 3.52 | 19.97 / 1.53 | 23.74 / 1.69 |
| RNA-short · unstranded × OFF | 6.00 / 0.61 | 3.78 / 0.27 | 3.79 / 0.27 |
| RNA-short · unstranded × ON (deferred) | 29.76 / 6.19 | 18.84 / 2.15 | 21.08 / 2.11 |
| RNA-long · stranded × OFF | 5.23 / 0.46 | 2.26 / 0.19 | 2.59 / 0.20 |
| RNA-long · stranded × ON | 7.38 / 1.62 | 4.48 / 0.68 | 7.09 / 0.74 |
| RNA-long · unstranded × OFF | 5.29 / 0.59 | 2.13 / 0.22 | **11.71 / 0.91** |
| RNA-long · unstranded × ON (deferred) | 28.99 / 5.84 | 6.86 / 0.95 | 8.73 / 1.05 |
| test chromosome · stranded × ON | 8.58 / 2.45 | 7.19 / 0.96 | 7.11 / 1.18 |
| test chromosome · ss 0.70 × ON | 10.28 / 3.32 | 6.91 / 1.51 | 7.79 / 1.98 |
| test chromosome · stranded × OFF | 8.20 | 7.39 / 1.08 | 7.52 / 1.08 |
| test chromosome · unstranded × OFF | 9.25 | 8.78 / 1.33 | 9.07 / 1.41 |
| test chromosome, centred probes · stranded × ON | 74.76 / 63.06 | 24.39 / 3.15 | 22.36 / 2.85 |
| test chromosome, junction probes · stranded × ON | 40.01 / 1.46 | 32.12 / 1.15 | 37.36 / 0.82 |

The real libraries (the CLI's gDNA fraction; only VCaP has a truth, 0.2518 by read name):

| library | 0.7.1 | tree (cap at 1) | prototype |
|---|---|---|---|
| LBX0190 (sparse plasma; the tree reads no reference) | 0.130 | 0.083 | 0.088 |
| MO_3021 (sparse plasma; no reference) | 0.225 | 0.155 | 0.157 |
| LBX0588 | 0.966 | 0.915 | 0.912 |
| VCaP mix | 0.263 | 0.236 | 0.236 |

On LBX0190 the prototype's population finds enriched levels near 106× and 13× the off-target density, at about 5 % of
the opportunity each, where the tree's mode reader finds none: a sparse panel gets its capture with no detector. The
ladder was still running when this was written (`ladder_pwfx.*` in the prototype directory).

Against the simulator's own yields (`ruler_vs_truth.py`): on the RNA-short captured row the prototype is far closer on
unprobed transcripts (median log error −0.23 against +3.08) and noisier on fully probed ones (10–90 %: −0.34 to +0.25
against −0.22 to +0.11); on the junction-probed layout it puts every unprobed class on one scale (1.00 against 0.95–1.02)
and spreads the partly probed isoforms 4× more within a gene (sd 0.219 against 0.054). That spread is the boundary
clamp of step 3: a probe at an exon's end makes the boundary there legitimately more captured than the exon's average,
and the clamp cuts it down. Without the clamp one false boundary at 7,000× survives on the zero-gDNA captured row.

## 6. What this means, and what is open

1. **The finding.** A ruler with no detector is exactly as good as calibration's per-object gDNA counts, and under a
   length gap those counts are wrong in both directions at hundreds to thousands of objects. The shipped ruler's `None`
   is not only a detector; it is what has been hiding those errors on every uncaptured library. Removing it, as the
   spectrum ruling requires, exposes them: RNA-long unstranded × OFF 2.13 → 11.71 %.
2. **The fix that generalizes** is calibration's: each object's gDNA in the density frame, so an object gDNA cannot
   sit in cannot hold gDNA, and an object RNA cannot sit in cannot hide it. That is the design already open in
   `ISSUES: calibration-detects-capture-on-a-capture-off-library` (the opportunity-aware likelihood introns have, given to
   every object), now with its measurement: the per-object tables of §2, read against `slot_truth.npz` in seconds from a
   dump (`dump.py`, `analyze_rs.py`). It is a change to ψ's prior in the native kernel, prototyped in a worktree.
3. **The ruler built here** is what the tool needs once the counts are right: an opportunity-weighted population read
   by its posterior median (one level when nothing is captured, any share of probed objects when something is), a floor
   at the off-target level, coherence for pieces gDNA cannot reach, relative weights, no reference, no clip, no `None`.
   Open in it: the boundary clamp (§5), the exon assumption (probes sit on exons; without it short introns between probed
   exons read captured), an accelerated population fit (99 % of the medians are final at 200 EM iterations), and the
   rulings it reverses (`DESIGN.md` §7.2's reference and clip by the spectrum ruling; the 2026-09-23 "each object reads
   its own count" for short exonic pieces, whose premise — a piece too short for gDNA holds no share — fails when RNA is
   shorter than gDNA).
4. **What neither fixes**: RNA-short stranded × ON is the capture × length coupling. A 78 bp RNA fragment overlaps less
   of a probe than a 250 bp gDNA fragment, so no gDNA-read level is RNA's; that is `CAPTURE_CLEAN_SLATE.md`'s
   per-placement kernel.
