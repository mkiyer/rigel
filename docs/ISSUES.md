# ISSUES — the issue log

This file is the issue log: one entry per open problem, question, decision or risk, and an append-only record
of what was measured and turned down. An OPEN entry is a `### kebab-name` heading with a `priority:` line (now
/ next / later / parked), the question in a sentence or two, the numbers a ranking turns on, and the
instrument that re-derives them. A CLOSED / REFUSED entry keeps its name, its verdict, the killing number and
the date, so a refused mechanism is not rebuilt. Cite an entry as `ISSUES: <name>`; names are the only
identifiers (`tests/test_no_jargon_labels.py`). What does not belong here: the ranked view (`ROADMAP.md`),
rulings and derivations (`DESIGN.md`, `EQUATIONS.md`), lessons (`TRAPS.md`), and any record of what was done —
the changelog is git.

---

## OPEN

Ordered by priority. An entry says what is open and the number a ranking turns on; what was done is git,
what was ruled is `DESIGN.md`.

### the-junction-price-is-noisy-within-a-gene
`priority: PARKED — the sum is accepted for now (owner, 2026-09-23); the pseudocount's odds are fixed since (2026-09-24, ISSUES: the-pseudocount-prior-is-biased-toward-gdna, CLOSED) · kind: defect · 2026-09-23`
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
within-gene sd line and the junction-probed rows); `quant_accuracy.py` per stratum above `--arm base_reseed` under
`--set em.assignment_mode=fractional`.

### splicing-artifacts
`priority: after the end-to-end work (owner, 2026-09-24): "splicing artifacts are a bigger problem than previously thought" · kind: defect · 2026-09-24`
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

### multimapper-intergenic-alignments
`priority: a separate feature, after the end-to-end work (owner, 2026-09-24) · kind: defect · 2026-09-24`
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

### the-em-answer-depends-on-where-it-starts
`priority: NEXT — the owner's investigation (2026-09-25); it holds the minus-strand coverage-weight fix · kind: defect · 2026-09-25`
What is left after the SQUAREM repair (`ISSUES: the-squarem-clamp-decided-which-components-live`, CLOSED): VBEM's
own start-dependence. With the clamp gone MAP reaches one answer from any start (204 fragments apart at
`g00 ss.99 ON`, all in loci still at the iteration cap), while VBEM still ends 8,096 fragments apart there (136
loci, 132 of them winner forks — a candidate dead in one answer and holding fragments in the other). Its E-step
weight ψ(α), with no per-component prior, penalises a small component by about −1/α, so components sharing fragments
race and the start picks the winner; the two answers' likelihoods differ (the uniform start's is higher in 121 of 136
loci), so it is not a flat valley. Candidate repairs, unmeasured: a warm start that is itself start-free (the MAP
optimum, which is unique), or a per-component prior. Instrument: `em_start_lab.py` in
`~/Downloads/rigel_runs/prototypes/2026-09-25_em_start/` (every EM setting re-solved from one pre-EM state in one
process, which repeats itself to 0.0 fragments; `analyze.py` classifies each moved locus). The minus-strand
coverage-weight fix — the warm start's trapezoid coverage weight was read from the flipped start in the wrong
direction, one fragment length off, on minus-strand transcripts (`~/proj/rigel-covwt`, `scoring.cpp`'s
`coverage_weight`, 20 tests) stays held (`TRAPS: a-cancelling-defect-pair`): on the old solver it moved the ladder ±1–2 %
(`g00 ss.99 ON` +7.2k transcript fragments, `g50 ss.50 OFF` genes +12.7 %); RE-PRICED on the repaired solver
(2026-09-25, `~/Downloads/rigel_runs/prototypes/2026-09-25_coverage_on_squarem/ab_table.txt`) it moves in-scope
transcript Σ|Δ| +3.3k (+0.8 %, stranded OFF), −0.8k (−0.05 %, stranded ON), +2.3k (+0.6 %, unstranded OFF), genes and the
gDNA pool within ±0.4 % bar `g05 ss.50 OFF` genes +1.3 % — still the residual start-dependence speaking, since the weight
only seeds the start. It lands with the start-free repair, where the seed cannot move the answer.
MEASURED 2026-09-25 on `g00 ss.99 ON` (the landed solver's prototype): the 136 loci VBEM's two starts disagree on carry
EQUAL truth error — 385,299 (coverage start) against 385,088 (uniform), MAP's single answer 385,376; the coverage start
is closer in 50 loci, the uniform in 64, 22 tie — so what is left is reproducibility, not accuracy: the forks pick
among near-equivalent explanations of genuinely ambiguous fragments.

### the-gdna-length-law-falls-back-at-identical-purities
`priority: NEXT — a defect found by the g98 dissection (2026-09-24) · kind: defect · 2026-09-24`
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

### the-scorer-reads-a-census-length-law
`priority: next, after the two above · kind: defect · 2026-09-24`
The E-step scores an unspliced fragment's length with gDNA's library census (`gdna_realized_pmf`,
`pipeline.py`) against RNA's spliced census. Under capture the length selection depends on where a fragment sits,
and every hypothesis at one footprint shares that footprint's capture — so the ratio wants the two origins' laws in
ONE frame, and two censuses that average capture over different positions tilt short exon-contained gDNA fragments
toward RNA even when both are estimated exactly. Confirmed free of the count and length confounds (g50 ss.99 ON, exact
count, one-scale lengths): the same law for both origins gives annotated −2.7k, spans 156.1k (truth 150.4k),
transcripts 1.70 %; the true census gives +12.7k, 139.2k, 1.79 %; shipped +5.2k, 158.6k, 1.75 %. At `g98 ss.99 ON`
(oracle arm) the scorer's law is worth 22.2k gDNA and 9.8k transcripts, 18.4k / 4.6k of it the estimator
(`ISSUES: the-gdna-length-law-falls-back-at-identical-purities`). Nil off capture. The same law for both is a
diagnostic, not a candidate: real libraries with different gDNA and RNA chemistry need the channel. The repair is a
derivation of the per-fragment length term at a shared footprint.

### the-pooled-q-in-the-gdna-count
`priority: later, with the capture repairs · kind: defect · 2026-09-24`
Calibration's per-locus gDNA count (`priors.assemble_priors`) converts each boundary's gDNA mass to fragments by
`boundary_mass_per_crossing`, pooled over gDNA and RNA, so where RNA dominates a boundary gDNA is over-stated even
under a perfect calibration (`prior_vs_oracle.py` O − S: +33.7k at `g50 ss.99 ON`, +2.9k at g98 ON, +1.9k off
capture). The EM passes a g50 count change through at only 0.37–0.50, so its size in the result is about 12k: with
gDNA's own share the g50 gDNA pool moves −18.2k → −30.4k, spans 142.0k → 152.1k (truth 150.4k), transcripts
unchanged. It cancels part of the capture likelihood's lean toward RNA today, so it lands with those repairs. The
gDNA LENGTH already uses gDNA's own share (`ISSUES: the-pooled-q-in-the-gdna-length`, replaced for the length only).

### the-pseudocount-strength-is-not-derived
`priority: next — after the implicit-splice gate; the odds are fixed (DESIGN.md §3.1c) · kind: open question · 2026-09-24`
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

### the-efficiency-posterior-floor-on-empty-pieces
`priority: MEDIUM — the unprobed class's scale at low gDNA; its EM cost unmeasured · kind: defect · 2026-09-23`
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
constant; the open question is what an object holding under one expected fragment should read. `ruler_vs_truth.py
--scale` (the unprobed rows); the census in `~/Downloads/rigel_runs/prototypes/2026-09-23_junction_price/`.

### strand-plug-in-bias-on-sparse-libraries
`priority: MEDIUM — first-order on real cfRNA, invisible on the ladder; no instrument in the tree (a scratch census of live objects over unspliced incidence, 2026-09-21)`
Calibration's per-object gDNA fraction is a plug-in that conditions on the object's own fragments, this one
included; the exact statement is the leave-one-out posterior, and the plug-in understates gDNA by −21 / −22 /
−19 / −11 / −2 % at 1 / 2 / 5 / 10 / 50 fragments. The pool-level bias is (live objects) / (unspliced incidence):
0.2–0.4 % on the ladder, where it is guaranteed invisible, and MEASURED 2026-09-21 on real data 0.216 / 0.120 /
0.308 on the cfRNA libraries LBX0190 / LBX0588 / MO_3021 (72–91 % of their live objects hold five or fewer
fragments and carry 15–44 % of the incidence) and 0.055 on the deep VCaP library. Any per-fragment use of the
fraction must take the leave-one-out posterior and must not apply the fragment's own strand term a second time.

### unwitnessed-loci-and-multimappers-at-the-em
`priority: MEDIUM — the two real-data cases the ladder cannot show; no instrument in the tree (the same scratch census)`
MEASURED 2026-09-21: every real library carries MultiLoci with a gDNA candidate and no witnessed gDNA (`P_g = 0`:
257 / 237 / 288 / 194 loci holding 374 / 723 / 941 / 2,182 gDNA-bearing units on LBX0190 / LBX0588 / MO_3021 /
VCaP; `L_g = 0` nowhere), and a multimapper share of the gDNA-bearing units of 13.3 / 4.0 / 13.2 / 1.0 %,
concentrated in the sparse loci (per-locus median 0, q90 1–5 % over loci with ≥ 20 such units). A kernel that
pins gDNA at calibration's count sends every unspliced unit of an unwitnessed locus to RNA and under-calls
gDNA by the multimapper share one for one; the ladder has neither case (0 of 1,213 loci). A per-fragment prior
prices each multimapper placement at its own objects, which is the cleaner rule.

### ruler-witness-geometry-on-transcript-panels
`priority: DEFERRED post-0.8.0 by owner ruling 2026-09-20 — Rigel stays panel-agnostic and takes no panel input; the junction price's precision is split out as NOW (2026-09-23) · kind: defect · 2026-09-15`
Only a transcript holding the junction a probe spans binds that probe whole, so the extra capture of
junction-spanning fragments is ISOFORM-SPECIFIC, and the EM splits a gene's shared fragments by the ratio of its
isoforms' capture-aware lengths. gDNA holds the probe's parts apart and binds the better one, and at zero gDNA
there is no witness at all; the shipped junction price reads it from gDNA at the junction's low and high
boundaries by conservation of bases, with no panel input — 0.959 of the simulator's own junction capture over 395
junctions of 60 probed transcripts in aggregate, but per junction too noisily
(`ISSUES: the-junction-price-is-noisy-within-a-gene`). MEASURED on the rebuilt ladder (the simulator binds every
fragment through one contiguous part of a probe since 2026-09-19; earlier numbers measured an asymmetric
half-match the owner ruled unphysical), stranded × capture-ON transcript Σ|Δ| as a share of the true RNA at
`g00` / `g05` / `g50` / `g98`, fractional: the per-base rule 6.53 / 3.49 / 5.50 / 31.03 %, the simulator's own
lengths (`quant_accuracy.py --arm oracle_ruler`) 1.34 / 1.71 / 2.90 / 24.66 %. The error sits in highly expressed
multi-isoform genes whose TOTALS are right, and a ruler is judged by the WITHIN-GENE spread of its error
(`TRAPS: judge-a-ruler-by-its-within-gene-spread`). The junction-probed test-chromosome twin predates the physics
change and is stale.
THE ONE SHARED RULE'S RESIDUAL (2026-09-23, 60 probed transcripts against the simulator): exact (1.000) on the
objects where a transcript's fragments are gDNA's — contained pieces and cuts away from junctions, exon edges and
ends, 41 % of a probed transcript's yield — and its junction cuts, 40 % of a probed transcript's yield, are what the
junction price covers. Its cuts within a fragment of an exon edge are captured 1.10× more than gDNA prices them
(19 % of the yield), and the simulator binds a fragment through its best single probe part, which captures a
little less than the sum when both exons carry separate probes: that is the ~2 % the classes stay apart with the
true gDNA on every object (`ISSUES: the-gdna-component-length-rule-differs-from-the-transcripts`).
THE DECISION IT WAITS ON: the probe design is the one observable that sees isoform-specific capture directly, and
Rigel reads no panel (`DESIGN.md` §7.2). Real panels often span junctions and the design file is usually
unavailable (owner, 2026-09-19), so the candidates are data-derived: the junction price from gDNA (shipped, its
precision open); a capture field fitted from the coverage shape around each probe and the gDNA footprint; a
per-kit capture profile learned across a cohort; and a sparsity prior on isoform support as the safety net. The
spliced reads at each junction are refused (`ISSUES: the-spliced-read-junction-price`). `quant_accuracy.py --arm
oracle_ruler`, `ruler_vs_truth.py`.

### per-transcript-prior-lane
`priority: HIGH — the per-transcript prior, the owner's problem after the capture-contracted length's junction price (2026-09-22, 2026-09-23) · kind: build · 2026-08-31`
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
`quant_accuracy.py`, per stratum above `--arm base_reseed`.

### overlapping-synthetic-shadows
`priority: an owner decision on the index; not a 0.8.0 number · kind: design · 2026-09-20`
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

### rna-prior-floor-at-pure-gdna-loci
`priority: NEXT — the owner's residual pass (2026-09-24): the largest owner of g98's gDNA error · kind: defect · 2026-09-20`
DISSECTED AGAIN 2026-09-24 on the count form, against per-fragment truth (four lenses, a synthesis and two refuters;
`~/Downloads/rigel_runs/prototypes/2026-09-24_g98_dissection/`). The read-out is the posterior MEDIAN of λ (`simplex_logodds.posterior_median_fg`,
`psi_kernel.h`, `DESIGN.md` §6c), not a mean: under the Jeffreys reference with no atom at zero RNA, an object whose
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
0.15 % of the RNA, invisible in the pool). `prior_vs_oracle.py`, `calibration_vs_oracle.py`,
`solvability_audit.py` (the confidently-wrong class). It is why the per-transcript allocation made `g98`
capture-ON worse before the gDNA opportunity was corrected (`ISSUES: per-transcript-prior-lane`).

### nascent-stress-sensitivity
`priority: next — it sizes `ISSUES: em-overturns-the-calibrated-gdna-split` · kind: question · 2026-08-22`
Does any in-scope verdict depend on the nascent stress level? The ladder runs `on_fraction 0.50`; realistic is
~0.10 (`DESIGN.md` §0b). MEASURED 2026-09-20 for the siphon: `g50 ss.99 ON` re-simulated at 0.10
(`~/Downloads/rigel_runs/suite/ladder_nrna_lo`) reads a siphon of +522,205 against +541,216 at 0.50, so that
verdict does not depend on the stress level — the repair is worth the full amount in the expected case. Re-simulate the worst in-scope scenario at the realistic level and check whether any
rank moves; a verdict that holds only at stress is a robustness finding. The first verdict it must size is the
capture-OFF gDNA over-call, measured only at stress (`g05 ss.50 OFF` 0.104 against 0.05). `sim/panel.py`,
`quant_accuracy.py`, `policy_benchmark.py`.

### performance-memory-bounded-solve
`priority: PARKED 2026-09-19 (owner: the method is the focus; machine work resumes on its own) · kind: build · 2026-08-17`
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
`profiling/profiler.py`, `profiling/sweep_replay.py`.

### gdna-landscape-trains-on-false-positives
`priority: later · kind: question · 2026-09-02; the population rule and the E-step landed 2026-09-10, the location floor 2026-09-14 (`DESIGN.md` §7.1)`
What is still open in the landscape estimator, each with its number: (a) under capture a short dark region's
mass is split over both modes of the previous fit (deferred 1.01–1.02×, stranded ON 1.004×); (b) the zero-RNA
controls move ±4 % (`g05 ss.99 ON` +4.4 %, `g98 ss.50 OFF` −1.4 %); (c) the Beta(½,½) reference still decides
a blind slot. CLOSED here 2026-09-14, (d): a slot whose solve is wider than one nat² in `log f_g` no longer
trains (`DESIGN.md` §7.1 rule 4) — the vertex-solved exons that trained 3,771 false fragments at the first
refit on `g00 ss.99 OFF` are out, and the four ladder zero controls read 282 / 194 / 265 / 172. Refused here:
excluding κ-dead exons (`g50 ss.50 ON` 2,691 → 56,422), AMBIG in the final fit (worse 25/32), and the two
other readings of the floor (`DESIGN.md` §7.1 rule 4). Its instrument, `landscape_training_census.py`, was retired 2026-09-14 (in git).

### multimapper-blind-support
`priority: a deferred session (owner, 2026-09-24: "add multimapping deposit to the accumulator and let calibration include multimapping fragments, which is correct behavior and something we have intended to add") · kind: defect · 2026-09-16`
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

### yield-variance-beside-the-count
`priority: later, with the per-transcript prior lane (owner, 2026-09-17: the release ships the contraction as it stands; performance next) · kind: build · 2026-09-17`
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
`priority: watch (no library in hand can test it) · kind: risk · 2026-09-17`
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

### message-layer-open-cases
`priority: next · kind: question · 2026-09-09`
Four residual cases, none a hole (`DESIGN.md` §6b.4–§6b.14): (a) an exon with both faces speaking — the
arrivals are summed and the price where both carry the same intron's claim is unmeasured; (b) the intron
factory on a region with both exon and intron bits, under capture; (c) the chain of termini — an empty
outside piece (median 12 bp) whose far face is another terminus, half the ladder's terminus-boundary error,
every upper side refused (`ISSUES: the-edge-upper-side`); (d) substrate `nest` (the prior already serves it,
0.666 vs 0.630) and the antisense's nascent variant (`docs/TESTING.md` §0a; `div` was built 2026-09-14).
`policy_benchmark.py --by-class`.

### capture-blind-gdna-divisor
`priority: the pre-EM prior chain, after the EM's residual · kind: defect · 2026-08-31`
`gdna_opportunity_from_index` is computed from the index alone, so under capture it removes ~6 bp of a ~30 bp
length selection — the gDNA control moved +6.0 % on all six capture-ON rows (gDNA has no introns to miss), and with
`ISSUES: eb-shrinkage-magic-ess` it owns the −5.90 % capture-ON length ceiling. `capture_eff_length` already models
the panel; it also blocks `ISSUES: crossing-pool-contrast`. ⛔ NOT PRICED ON THE DELIVERABLE: the reading once
quoted here (`quant_accuracy.py --arm oracle_efflen`, 2026-09-19, −50 / +22 / −1,104 fragments) is WITHDRAWN
(2026-09-21) — that arm cannot fire. The oracle overrides only the six count arrays
(`prior_vs_oracle.OVERRIDE_FIELDS`, `OracleTruth.override_masses`) while the locus gDNA length reads its conserved
shares and the published efficiencies and never a count, so `oracle_efflen` is `base` under another name and its
deltas were run-to-run noise (`TRAPS: could-the-arm-have-fired`). The arm is retired or given the simulator's
per-locus gDNA opportunity as a real truth before any number is quoted for it again.

### eb-shrinkage-magic-ess
`priority: the pre-EM prior chain, after the EM's residual · kind: defect · 2026-08-31`
`POOL_EB_PRIOR_ESS = 1000.0` shrinks the gDNA pmf toward `global_pmf` (mostly RNA whenever gDNA is a minority)
at a magic ESS: inert on the ladder (0.01 bp), dominant on the fl-gap arm at `g05` capture-ON (`ship−pool`
−23.7 of −31.7 bp). Replacement: reconcile the pools by their precision (`EQUATIONS.md` §6c). Its
instrument, `fl_pool_purity.py`, was retired 2026-09-14 (in git). The ladder's deliverable cannot see it
(`ISSUES: capture-blind-gdna-divisor` has the `oracle_efflen` number); only the fl-gap arm can.

### refit-vs-message-arbitration
`priority: next · kind: design · 2026-08`
At the unstranded × capture-OFF exon cell the refitted gDNA prior and the message both impute one slot with
nothing arbitrating them; the message is the accurate voice there and the refit displaces it. Re-read under the
E-step (with an instrument since retired): the prior does the unstranded rows and the messages the stranded
capture-ON ones. `policy_benchmark.py --by-class`, `calibration_vs_oracle.py`. Belongs with `ISSUES: gdna-landscape-trains-on-false-positives`.

### the-atom-at-an-unwitnessed-both-strand-slot
`priority: later — accepted as a limit of the information (owner, 2026-09-14) · kind: known limit · 2026-09-14`
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
`priority: later (deferred by the owner 2026-09-13) · kind: problem · 2026-09-04`
On unstranded data an exon's held row is the intron's composition through the face map, whose upper side is
the map's plateau above its ceiling — a lower bound on gDNA — so at pass zero an unstranded licensed exon
reads ~9× its true gDNA (+2,000 % at `g25 ss.50 OFF`) and forwarding it compounds the bias. The toy harness
gate reads a factor of 88 (exon |Δf_g| 0.762 beside a pure-gDNA intron, 0.0086 beside a nascent-bearing one)
and is a strict xfail citing this entry; on the ladder the class is 4–6 % of the unstranded error with
`transfer` at parity there, on the test chromosome 20–32 % and transfer worse (2,647 → 4,056 at
`g50 ss.50 OFF`). Every cap, two-sided row or wall is refused where junction probes enrich the flux more than
the crossing (`ISSUES: the-certified-flux-row-as-a-level`, `ISSUES: two-sided-exon-row-forms`,
`ISSUES: the-discrepancy-priced-cap`, `ISSUES: the-two-sided-level-lane`,
`ISSUES: the-wall-above-the-face-map-ceiling`). What would license a two-sided level is whether this library's
gDNA is enriched, witnessed on unstranded data by silent genes' exons against the intergenic density — the
landscape prior's job. The first step when taken up is that witness's derivation; judge at `R exon
(licensed)` and the walled classes, halves apart. `policy_benchmark.py --by-class`.

### the-lower-bound-noise-ratchet
`priority: later (with the enrichment witness) · kind: defect · 2026-09-05`
A level from an RNA-rich node's own strand profile has a mode that is noise around zero; its lower side bounds
its neighbours and the tightest noisy neighbour wins (ladder `g05 ss.99 OFF` 44,714 → 45,076; test chromosome
`g00 ss.99 ON` 64 → 216, 187 of it at `capcluster_ab`'s inner termini). A two-sided own-profile level keeps
fewer fragments overall, so lower-only stays. The gDNA lane's EDGE level does it too (2026-09-14,
`test_encompassing_locus.py`, its xfail): a shallow single-strand flank (404 fragments, truth 0.530, local
solve 0.546) reads 0.596 under the intergenic neighbour's Poisson level — a lower bound at that neighbour's
sampled density, 0.298/bp against the flank's realised 0.27, a 1.6σ excursion the hop's price blurs but does
not move. The cure is the enrichment witness `ISSUES: two-sided-exon-row` waits for. The ladder's `g00` rows
(`calibration_vs_oracle.py`) and the `test_encompassing_locus.py` xfail.

### flux-floor-dispersion
`priority: later (with the transport-dispersion decomposition) · kind: question · 2026-09-08`
The certified flux at a junction is a lower-sided estimate of the exon's strand RNA level priced by the pair
(`DESIGN.md` §6b.13); the route rate scatters beyond counting (median −3 %, 5–9 % at depth; 0–40 % over on
nine block readings) and at a pure-RNA exon the price cannot see it, so a lucky over-read is a sharp floor a
few points too high. Its decomposition instrument, `transport_dispersion.py`, was retired 2026-09-14 (in git).

### splice-out-premise-bias-uncorrected
`priority: later · kind: decision · 2026-09-02`
The splice-out message assumes spliced and unspliced fragments at one face share capture affinity; measured,
the premise fails as a bias under capture (`log a` ≈ 0 off, +0.28 under exon probes, +0.78 under junction
probes), so only a subtracted level — a cross-locale fudge, the owner's call — would correct it. No instrument
yet.

### the-tilt-census-as-an-instrument
`priority: later · kind: build · 2026-09-13`
"Where does the strand tilt matter, and how does the tool do there" is answered only by a scratch script: per
stratum, the AMBIG slots by the RNA on each strand (bands of the minor strand's share), their depth, the gDNA
error, the TILT read-out's error (`|Δτ|·R/2`, RNA fragments on the wrong strand — measured by nothing else) and the
predicted θ peak width. `policy_benchmark.py --by-class` ranks by node class and not by strand split or tilt
error, so the tilt-error column is the new part; a census of what it overlaps comes before it becomes an
instrument (the instrument ruling: few instruments, kept current).

### transfer-variance-premise
`priority: later · kind: question · 2026-08`
Does a hop's transfer variance price a ratio built on a handful of counts? The policy prices every hop by both
witnesses' counting plus the pair's disagreement (``hop_price``, `native/transfer_rows.h`); whether that is right where a
pair agrees by coincidence is the open half (the landscape is no substitute: ~10× over-stated). `EQUATIONS.md`
§3.5h.

### expand-the-gdna-spectrum
`priority: later · kind: decision · 2026-08`
Fill the gDNA spectrum (1, 5, 10, 25 % up past 90) without multiplying benchmarks: a level is justified by a
measured transition and crosses a reduced set of the other axes until an interaction is shown; each condition
costs a simulate, two caches and a certification (`sim/panel.py`). See
`ISSUES: flgap-panels-stale-nascent-model`.

### flgap-panels-stale-nascent-model
`priority: later · kind: decision · 2026-08-22`
The two fl-gap side panels were not regenerated in the sparse-nascent rebuild and carry the retired uniform
nascent model, so a claim spanning the ladder and a side panel varies two things. Re-simulating is the owner's
call.

### hygiene-ledger
`priority: later · kind: hygiene · 2026-08-31; the review of 2026-09-22`
What the review left, each its own commit:
(a) ROTTEN BUT LIVE, each moving a toy's or an instrument's numbers when repaired: the toy harness's `harvest`
(`tests/calibration/_toy_harness.py`) calibrates its donor undrained and without the two-pool contrast, so every
toy inherits both; `quant_accuracy`'s oracle arms are undrained (documented there).
(b) DEAD, DEFERRED: `region_span_count`, tallied per fragment by the accumulator and carried through the payload,
the substrate and the caches, read by nothing (the retired length channel's) — deleting it changes the payload
schema and re-caches both panels, so it goes with the next change that re-caches anyway (owner, 2026-09-22).
(c) PRODUCTION-DEAD, KEPT AS TEST SURFACE: `FragmentBuffer.append` / `finalize` / `__iter__` and
`BufferedFragment` (the scanner hands chunks over; only `test_buffer.py` appends); the resolver's `intron_bp`, the
tested half of its overlap profile; `rigel.sim.benchmark` and `Scenario.build_oracle`, test tooling inside the
package; `priors.contended_boundaries`, `region_init.has_own_composition_evidence`, `effective_length.fl_mean`
and `Strand.from_is_reverse` / `opposite`, which tests use as checks or oracles.
(d) CLAIMS NOT RE-DERIVED on the current tree, left standing: that most in-scope error sits at the simplex
vertices (`simplex_logodds`, a relay-era measurement); `sweep`'s refused deferral of UNIDENTIFIED slots to the
prior (priced 2026-07, with no refusal entry here); `region_geometry`'s "no per-region spliced floor" A/B
(relay-era); and `fl`'s crossing pools called "gDNA by structure" because mature RNA never crosses an
exon|intron boundary, while RNA that has not spliced there does.
(e) GATES THAT CHECK LESS THAN THEIR NAME: `test_sweep.py::test_gdna_sweep_zero_gdna_pin_and_monotone` asserts
`f_g < ½` on slot 3, the AMBIG|intron− boundary, while the AMBIG region it means ends at 0.949 under the silent
policy, and nothing in it checks monotonicity; the mature-exon chain's tests give the same `f_g` with the junction
spliced or not, so none of them exercises the junction reads;
`test_conserved_mass.py::test_the_mass_is_the_PER_BASE_attribution` allows `count` whole fragments where its docstring
claims `count` half-ulps; `test_gdna_strand_fit.py` and `test_region_geometry.py` still build float64 banks as uint64
fixed point in their fixtures.
(f) SMALL: the `--no-mappability` flag names a store whose mappability the index no longer reads (it opts out of
the splice blacklist); `_solve_impl` builds without the other native modules' `-march=native` / IPO flags
(`CMakeLists.txt`), to confirm deliberate; three unread scalars in `solve_kernel.cpp`'s positional aggregates;
the index's duplicate map as an alias map `dropped_t_id → kept_t_id` (an index rebuild, no panel re-scan).
Kept as coverage GAPS, not dead code: the five CLI command bodies, the silent policy through `calibrate`, the
simulator's sharded writers and its whole-genome grid, and the zarr splice blacklist.

FOUND 2026-09-24: `scripts/design/solvability_audit.py`'s `main()` calls `_oracle_arms.truth_f_gdna`, deleted in
34145493, and raises `AttributeError`; `preflight.py`'s import check cannot see it (the g98 dissection ran it
through a shim).

### drain-contaminates-certified-rna
`priority: later (parked by the owner, 2026-09-01) · kind: defect · 2026-08-31`
The second pass deposits some true-gDNA fragments into the certified-RNA banks: 233 records at
`g50 ss.99 OFF`, 1,482 at `g98 ss.99 ON` (1.9 % of that in-scope channel). The leak is exactly posterior
sampling and a no-leak counterfactual moves the 0.8.0 metric −0.60 %/−0.05 % at worst, so the harm is the
certainty claim and no in-solve correction is licensed (`ISSUES: drain-provenance-split`). To build: `DrainQC`
records `Σ q_null`; the middle-bin posterior bias repaired as its own A/B. `calibration_vs_oracle.py`.

### crossing-pool-contrast
`priority: parked (blocked) · kind: question · 2026-08-31`
A second gDNA length contrast on the crossing pools: with oracle weights it beats the contained one under
capture (TV 0.076 vs 0.136; 0.078 vs 0.182) and starves off capture (pool 3 at 29–630 fragments). Blocked: the
weight estimator does not transfer under capture (`ISSUES: capture-blind-gdna-divisor`), and no shadow
transcript overlaps a gene edge, so pool 3 reads exactly 1.0000 pure
(`TRAPS: purity-is-a-property-of-the-annotation`).

### capture-degeneracy-standing-risk
`priority: parked (watch) · kind: risk · 2026-08-31`
The gDNA two-pool contrast survives capture by a degeneracy: the shared-contaminant assumption is false under
capture (TV 0.95 vs 0.06–0.14 off) and it is safe only because the intergenic pool is depleted-not-impure, so
`a_0` clips to 1 and the algebra collapses to `g = f_0`. A probe panel that put RNA into intergenic space would
break it silently; `_deconvolved_gdna_counts` carries the derivation. No panel can fire it today.

### pure-rna-mirror-asymmetry
`priority: parked · kind: defect · 2026-08`
Two exact per-fragment mirrors of a pure-RNA library deconvolve differently in `count_gdna_region` by a few
percent, neither boundary-only nor monotone in strandedness. An R1-sense library is simulable; no instrument
yet.

---

## CLOSED / REFUSED — do not rebuild these; append-only

Every entry keeps its stamped measurement exactly as recorded: a graveyard row without its number is an
invitation to rebuild. A row measured on "all 36 conditions" or quoting `g01`/`g10`/`g25`/`g75`/`g90` predates
the ladder retired 2026-08-13 — the verdict stands as a record, and re-opening one means re-running it on the
current panel. Where a mechanism's only target was unstranded × capture-ON the row is moot as a 0.8.0
candidate on top of being refused; the `g00` zero-control column is never moot.

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
relative (the plain step is now copied at step 1, not recomputed).

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

### a-junction-price-clipped-at-one
REFUSED 2026-09-23. The sum by conservation of bases is floored at 0 and not clipped at 1 — each boundary is an
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
fired) and `test_encompassing_locus` on a gDNA-free donor.

### measured-prior-rung-4
SUPERSEDED 2026-09-14: ψ's composition reference fitted from composition-free observables (`rho_0`, the
enrichment responsibility, the detector as a boolean) was the measured-prior thread's fourth rung
(2026-08-26). The gDNA landscape hyperprior (`DESIGN.md` §7.1, trained on the previous solve's deconvolved
gDNA, the zero controls at a few hundred fragments) is the population prior ψ carries, and ψ has no reference
location (`DESIGN.md` §6b.1); the total-density landscape it would have consumed is QC-only. The requirements
it listed (span both lattice ends, exact at `g00`, one pseudo-fragment) are met or moot. Re-open only with a
measured gap the landscape prior cannot close.

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
rows 1.4–2.5×. The gate stays a strict xfail.

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
| **a nascent-bearing ladder condition** | toy, 36 conditions × 7 rungs | ⚠ **−5 %**, and the wrong way on one stratum. Keep it as a harness arm (`--nrna 60`); it no longer justifies re-simulating the panel |
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
