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

### performance-memory-bounded-solve
`priority: now · kind: build · 2026-08-17; the port and the block in one native call landed 2026-09-17/18 (`DESIGN.md` §6b.15)`
A deep run must be fast enough to iterate on, and memory-bounded. THE PORT IS DONE and its record is
`DESIGN.md` §6b.15.1–§6b.15.5: the chain is solved a locus block at a time, the block size moves no number,
and the whole sweep is ONE native call over a pool of threads, bit-identical at every thread count. THE
BALANCE THEREFORE CHANGED, and this entry is ranked against the new one, not the old.

THE BASELINE (VCaP, 18,568,456 fragments, `--threads 8`, the tree at `e921869c`, two interleaved pairs of
2026-09-18, `perf/cache_deleted_2026-09-18/`): the run 148.4 s and 145.1 s, peak RSS 10.6 GB and 10.3 GB.
By stage, from the faster pair: calibrate 40.6 (the four sweeps 24.7, the landscape fits 4.2, its own Python
10.6), quant 34.2 (the locus EM 14.3, the capture effective lengths 6.0, scoring 5.9, the priors 3.7, the
partition 3.0), the scan 33.0, the second pass 21.6 (scoring the held fragments 11.8, a fragment-length fit
8.3, the drain 1.3), a second fragment-length fit 8.0, the index load 6.5. Calibration is 28 % of the run and
the sweeps inside it 17 %: the reducible work is now OUTSIDE calibration, and nearly all of it is Python
doing per-object work that numpy or the C++ beside it already does per array. ⛔ Read a saving from the
STAGE row of an interleaved pair, never from the wall: a few seconds is inside the wall's own drift.

THE ATTRIBUTION, one `profiler.py --cprofile` run on the same library (2026-09-18; it inflates the wall to
172 s, so its numbers are SHARES and never timings), by call count per run: 2,708,183 scalar
`np.searchsorted` calls in the second pass's scoring, from `_exact_region_bound` (1,853,986) and `_sj_id`
(926,993), both of which say in their own docstrings that they mirror `Accumulator::sj_edge_id` and
`Accumulator::exact_region_bound`; 1,429,451 `np.mean` calls on lists of at most two, one per exon per
fragment-length fit; 808 `poisson_lower_mean` calls, four `one_sided_rate` fits at 202 bisection steps where
about 60 close a float64 bracket; 1,982 short-length table builds in `effective_length.interval_sums`; and
one region-to-locus overlap computed twice although `_region_locus_shares` says "computed exactly once".

WHAT IS OPEN, ranked, each its own commit with its own gate: ① ~~the exact micro-wins~~ LANDED 2026-09-18,
bit-identical: the bisection ends when its bracket does (the magic 200 gone), the Poisson identity's log-gamma
is a table at its integer argument, the region-to-locus overlap is traversed once per assembly — the rate fit
261 → 46 ms on 80,000 objects at the same rate to every digit, four fits a run, and one 1.7 s traversal gone; ② the scan's
thread split (≈ 8 s, MEASURED, `ISSUES: scan-thread-split-starves-the-workers`; ⛔ a different worker count
changes who deposits, so the real-library identity reference decides whether it is free or priced); ③ the
second pass's boundary lookups — bind the two C++ lookups BATCHED, prove the id spaces agree, restructure
the loop around one pre-pass per reference, delete the Python mirrors (≈ 7 s, and the one-path duplicate goes
with it); ④ the two fragment-length fits — the adjacent-pair loop and the per-exon filter vectorise BIT-EXACTLY
(`np.bincount` accumulates in input order, and an exon has at most two flanking boundaries, so a grouped mean
equals `np.mean` on the list), while the surviving accumulation is the one item that may move a number and is
isolated and priced for that reason (≈ 11 s of 16.3); ⑤ memory — the sweep's arena is 1.19 GB at eight threads
and `CalibrationConfig.sweep_block_slots` scales it linearly while moving no number (`block_slots = 1000`
replays the first sweep at 3.01 s against 2.82 s), while the run's own peak sits in quant, whose 2.2 GB of
scored candidates is measured before anything is proposed; ⑥ what is left after that — `_cum_short`'s 1,982
table builds, the locus EM's 14.3 s, the index load's 6.5 s, the scanner's own throughput, and the whole
problem past 100 M fragments.

⛔ THE GATE FOR EVERY ITEM BUT ⑤: the sweep replay covers the sweep and nothing else, so what sees these is
`rename_identity.py --check` on the three frozen references — and only the `--bam` one runs the scan and the
second pass. Not to do: micro-optimise inside the kernels; bake the λ lattice into the port
(`sweep_logodds_step` is a parameter); trade a calibration number for speed anywhere except where an item
says it is priced. `profiling/profiler.py`, `profiling/sweep_replay.py`.

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

### ruler-witness-geometry-on-transcript-panels
`priority: later (measure on real panels first) · kind: limit · 2026-09-15, amended 2026-09-16`
gDNA is captured in genomic coordinates and a transcript's RNA in its own, and the two differ within a fragment
length of every splice junction a probe spans and at every exon shorter than a fragment. The ladder's panel is
designed in transcript coordinates, so its probes span junctions, and the simulator captures gDNA at a split
probe at `gdna_split_penalty` 0.2 of the cDNA's weight: on the ladder's `g05 ss.99 ON` row the probed
transcripts scatter −0.5 to +0.3 nat (35–36 % within ±0.1, under the shipped ruler and the expectation ruler on
the per-base length alike) even when fed the CERTIFIED TRUE gDNA counts (20 %), where the test chromosome's
benign panel, which spans no junction, reads 99 %. The other direction is the tiny-exon block: a probe centred
on a 40 bp exon binds a gDNA fragment over its full 125 bp while the simulator's non-stacking rule binds a spliced
fragment over one exon's 40 bp, so the probed `captiny` transcripts read +1.03 nat against the sampler's truth
with the mechanism reading their edge crossings exactly (`ruler_vs_truth.py --condition
gdna_g05_ss_0.99_nrna_file_capture_on --out`, the `captiny` rows). The gDNA witness reads the panel's capture
of gDNA, not of cDNA, by a factor the panel's design sets, and no estimator on gDNA alone can see it; the
annotation knows where the junctions and the tiny exons are, so a correction from the probe design is in
principle observable. The ladder row's unprobed transcripts hold no gDNA fragment at all and still read
+0.45 nat under the expectation ruler (+2.0 without the floor on the object set): the calibration assigns gDNA
to exons that have none, the stranded capture-ON composition residual the standing numbers already carry, which
the ruler inherits and cannot repair.

### multimapper-blind-support
`priority: next after the port (the first ruler question on real libraries) · kind: defect · 2026-09-16`
The accumulator drops every fragment with `NH > 1`, so a piece whose gDNA fragments are multimappers holds no
count the calibration can see, and the expectation ruler — which reads a piece's efficiency from its own count and
its crossings against the fully captured level — reads it at the DEPLETED level: the transcript is contracted as if
unprobed. Measured read-only on the two captured real libraries after the ruler's repair (the session's
`multimapper_check.py`: per transcript with ≥ 20 fragments over its pieces, the share of multimapping reads among
them from the BAM's `NH` tag, then the ruler factor's quantiles per share bin): LBX0588 (reference 2.977e-01/bp, 11,579
kernels) median factor 0.076 / 0.003 / 0.013 / 0.001 at a share below 1 % / 1–10 % / 10–50 % / ≥ 50 % (n 36,446 /
6,161 / 1,429 / 928), the share below 1/10 rising 55 → 90 %; VCaP (8.593e-02/bp, 24,496 kernels) 0.231 / 0.149 / 0.028
/ 0.002 (n 280,227 / 35,494 / 2,237 / 1,103), the share below 1/10 36 → 95 %. A hundredfold gap between the
unique-mapping and the repeat-rich transcripts, monotone in the share, on both libraries; the two libraries without
a reference read every transcript at 1 and are uninformative. The confound is real and unmeasured here — a repeat-rich
transcript may also be unprobed by design — so the number is an upper bound on the blindness, not its size; the
floor `C/(C+1)` the repair deleted had hidden it behind a +3.4 nat bias on every unprobed transcript
(`ISSUES: ruler-multimapper-floor-caps-the-correction`, CLOSED). The repair is in the OPPORTUNITY, not the posterior:
a piece's contained and crossing supports should count only the starts a fragment could be uniquely placed at — a
mappability of the support from the index, so that a wholly repeated piece has `S = 0`, no evidence, and reads the
population's clipped mean rather than the depleted level. Instrument: `ruler_vs_truth.py` cannot see it (the
simulator's reads are unique); the truth on a real library is the multimapper share against the factor as above,
and the repair's gate is that the four bins read alike.

### yield-variance-beside-the-count
`priority: later, with the per-transcript prior lane (owner, 2026-09-17: the release ships the contraction as it stands; performance next) · kind: build · 2026-09-17`
The capture-contracted yield is a posterior mean and carries no variance, so a count on a small yield reads as a
large abundance with the same error bar as any other. Measured on VCaP by drawing every piece efficiency from its
posterior on the landscape grid and re-running the EM (eight draws, the EM's seed fixed; the session's
`yield_draws.py`, scored by `draws_census.py`): the counts themselves are stable — the dominant isoform under the
draws' mean equals the point estimate's in 98.5 % of multi-isoform genes with ≥ 20 fragments and is identical in
every draw for 90.6 % — while the yield's uncertainty is a CV of 0.006 / 0.061 / 0.317 at the 10th / 50th / 90th
percentile of transcripts with ≥ 20 fragments, above the Poisson CV (median 0.100) for 36.4 % of them. So the
variance need not be propagated through the EM; it is a per-transcript number to publish beside the count: the
yield's posterior sd from the pieces' posterior variances, `Var[eff_t] = Σ_p (ℓ_p^τ)² Var[c̃_p]` under the
landscape (the pieces' posteriors are what `capture_efficiency` already computes; the second moment is one more
`w @ clipped²`), so that a user's abundance carries `CV² = 1/k + Var[eff_t]/eff_t²`. Derive → gate against the
draws' spread on VCaP → `src/`: a result field per object and an `em_effective_length_sd` column. Its consumer is
the allocation of the pooled RNA prior across a locus's transcripts (`ISSUES: per-transcript-prior-lane`): today
the transcripts duel for the ambiguous fragments with no per-transcript prior, and a variance per transcript is
what an allocation better than equal shares would read (owner, 2026-09-17). Instrument: the
draws' spread is the truth the analytic form is gated against.

### capture-premise-untested-on-cdna
`priority: watch (no library in hand can test it) · kind: risk · 2026-09-17`
The ruler reads the panel from gDNA and applies the same efficiency to cDNA: the premise that off-target cDNA is
depleted as off-target gDNA is, which the simulator satisfies by construction and which cross-hybridisation of
cDNA to paralog probes or nonspecific binding could make milder by an order of magnitude. Its exposure is not the
scale (TPM is on the plain length by ruling and counts see only ratios within a locus) but the isoform split:
on VCaP 4,677 of 12,039 multi-isoform genes with ≥ 20 fragments change dominant isoform between the contracted and
the plain yield (`DESIGN.md` §7.2, the yield's two consumers), two thirds of the isoform-level fragment mass with
them, stable under the yield's posterior — so the flips are the model's answer, right if the premise holds. The
one experiment that settles it: a sample sequenced both ways, panel and whole transcriptome — the count ratio on
unprobed transcripts is the cDNA depletion directly, and the ruler gives the gDNA depletion on the same pieces. The
lab has no such pair (owner, 2026-09-17); the session's `lever_census.py` is the read-only smoke test to re-run on
any new captured library, and `count_unambig` beside `count` is what tells a user which isoform assignments rest on
shared fragments alone.

### nested-antisense-leak-under-the-sane-ruler
`priority: later (EM-side, with the per-transcript prior lane) · kind: defect · 2026-09-14`
With the EM's ruler honest — a gDNA-free library contracts nothing (`DESIGN.md` §7.2) — the negative control
of `tests/scenarios/test_antisense_intronic.py` (a single-exon antisense `t2` inside the host's intron,
truth 0, host nascent RNA at 50) receives the strand-flipped intronic nascent fragments: 24 of 2,000 at
SS 0.9 and 124 at SS 0.65 (the bounds were 5 and 20; both parametrisations are strict xfails). The bounds
had held only because the retired kernel-density reference, fabricated from 1.1 false gDNA fragments,
contracted the host mRNA 4.7× and its nascent entity 3.9× while leaving `t2` at full length, so the
nascent entity's rate per base was inflated fourfold and nothing reached `t2`. The defect is the EM's
assignment at a nested transcript nothing witnesses — the EM-side twin of
`ISSUES: the-atom-at-an-unwitnessed-both-strand-slot` — and its lever is the per-transcript prior lane
(`ISSUES: per-transcript-prior-lane`), not calibration. `quant_accuracy.py`.

### message-layer-open-cases
`priority: next · kind: question · 2026-09-09`
Four residual cases, none a hole (`DESIGN.md` §6b.4–§6b.14): (a) an exon with both faces speaking — the
arrivals are summed and the price where both carry the same intron's claim is unmeasured; (b) the intron
factory on a region with both exon and intron bits, under capture; (c) the chain of termini — an empty
outside piece (median 12 bp) whose far face is another terminus, half the ladder's terminus-boundary error,
every upper side refused (`ISSUES: the-edge-upper-side`); (d) substrate `nest` (the prior already serves it,
0.666 vs 0.630) and the antisense's nascent variant (`docs/TESTING.md` §0a; `div` was built 2026-09-14).
`policy_benchmark.py --by-class`.

### per-transcript-prior-lane
`priority: next · kind: build · 2026-08-31`
`rna_prior_weight` is built end to end but `pipeline.py` omits it, so the shipped EM carries no per-transcript
information; a perfect per-transcript prior roughly halves in-scope gene-level error. Two weightings are
refused (`ISSUES: refused-transcript-weights`, `ISSUES: refused-soft-min-path-weighting`): the support problem
is the whole problem, so the next candidate is a sparsity mechanism, targeting expressed multi-exon
transcripts with median exon ≤ 150 bp. `quant_accuracy.py`.

### scan-thread-split-starves-the-workers
`priority: next · kind: decision · 2026-09-11`
`BamScanConfig.resolved_scan_threads` gives BGZF decompression `min(4, total − 1)` threads and the scan
workers what is left, so a 2–4 thread budget runs ONE worker and decompression is not the bottleneck. Scan
seconds on the 18.6M-fragment library by (bgzf, workers): total 4 — (3,1) 113.5, (2,2) 59.5, (1,3) 42.4,
(0,4) 33.7; total 8 — (4,4) 34.4, (2,6) 26.4, (1,7) 24.8; total 16 — (4,12) 19.1, (1,15) 19.5, (2,14) 17.4.
Any new split rule is a tunable and `--scan-bgzf-threads` is a user-facing flag, so the rule is the owner's.
`profiling/profiler.py --scan-only`.

### capture-blind-gdna-divisor
`priority: next · kind: defect · 2026-08-31`
`gdna_opportunity_from_index` is computed from the index alone, so under capture it removes ~6 bp of a ~30 bp
length selection — the gDNA control moved +6.0 % on all six capture-ON rows (gDNA has no introns to miss),
and with `ISSUES: eb-shrinkage-magic-ess` it owns the −5.90 % capture-ON length ceiling. `capture_eff_length`
already models the panel; it also blocks `ISSUES: crossing-pool-contrast`.

### eb-shrinkage-magic-ess
`priority: next · kind: defect · 2026-08-31`
`POOL_EB_PRIOR_ESS = 1000.0` shrinks the gDNA pmf toward `global_pmf` (mostly RNA whenever gDNA is a minority)
at a magic ESS: inert on the ladder (0.01 bp), dominant on the fl-gap arm at `g05` capture-ON (`ship−pool`
−23.7 of −31.7 bp). Replacement: reconcile the pools by their precision (`EQUATIONS.md` §6c). Its
instrument, `fl_pool_purity.py`, was retired 2026-09-14 (in git).

### refit-vs-message-arbitration
`priority: next · kind: design · 2026-08`
At the unstranded × capture-OFF exon cell the refitted gDNA prior and the message both impute one slot with
nothing arbitrating them; the message is the accurate voice there and the refit displaces it. Re-read under the
E-step: `calibration_walk.py` now says the prior does the unstranded rows and the messages the stranded
capture-ON ones. Belongs with `ISSUES: gdna-landscape-trains-on-false-positives`.

### prior-fidelity-vs-deliverable
`priority: next · kind: question · 2026-08`
Why is prior fidelity anti-correlated with deliverable quality? Leading answer: at the worst slots the
self-solve with the fitted prior is nearly right and the messages destroy it (measured at a retired rung;
confirm on a second stratum). The ruler is out of the way (`ISSUES: g00-shrinkage-upstream-repair`,
CLOSED). `prior_vs_oracle.py`.

### antisense-prior-assembly-casualty
`priority: the prior-assembly session · kind: decision · 2026-08-18 (named 2026-09-13)`
`tests/scenarios/test_antisense_intronic.py::test_nrna_multiexon_t2_low_ss` is a strict xfail: `assemble_priors`
pins synthetic nascent RNA at Dirichlet alpha = 0 (`EQUATIONS.md` §9b), so recovered RNA lands on the annotated
antisense t2 — 72 today against the test's limit of 50. The owner plans to change the alpha = 0 rule in the
post-calibration prior-assembly session; the xfail closes there, with a test that asserts the new rule's promise.

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
not move. The cure is the enrichment witness `ISSUES: two-sided-exon-row` waits for. `zero_controls.py`.

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

### background-abundance-pair-unruled
`priority: later · kind: decision · 2026-09-13`
`CalibrationConfig.background_abundance` chooses which (counts, exposure) pair the pooled gDNA background
takes: `"contained"` (the count over the gDNA contained effective length — unbiased, the fragment-length pmf
in the divisor) or `"measured_total"` (the START/END banks over the region's own length — pmf-free, refuses
without the wall inputs). The two agree off capture and part under it, where the contained divisor
under-reads the true gDNA rate several-fold while the pmf-free pair over-reads on pools carrying nascent RNA.
A design decision, not a tunable: rule which pair ships
(`calibration_vs_oracle.py --set calibration.background_abundance=measured_total` prices the swap) and the
field goes with the ruling.

### the-tilt-census-as-an-instrument
`priority: later · kind: build · 2026-09-13`
"Where does the strand tilt matter, and how does the tool do there" is answered by a scratchpad script
(`tilt_census.py`): per stratum, the AMBIG slots by the RNA on each strand (bands of the minor strand's
share), their depth, the gDNA error, the TILT read-out's error (`|Δτ|·R/2`, RNA fragments on the wrong
strand — measured by nothing else) and the predicted θ peak width. Before it becomes a `scripts/design/`
instrument, census what it would replace (`worst_objects.py` and `policy_benchmark.py --by-class` rank by
class, neither by strand split or tilt error): the tilt error column is new, the rest overlaps. The
shared-exon toy (`deep_stress.py`, spliced and mono) is the natural rung to add to the toy ladder.

### transfer-variance-premise
`priority: later · kind: question · 2026-08`
Does a hop's transfer variance price a ratio built on a handful of counts? The policy prices every hop by both
witnesses' counting plus the pair's disagreement (``hop_price``, `native/transfer_rows.h`); whether that is right where a
pair agrees by coincidence is the open half (the landscape is no substitute: ~10× over-stated). `EQUATIONS.md`
§3.5h.

### nascent-stress-sensitivity
`priority: later · kind: question · 2026-08-22`
Does any in-scope verdict depend on the nascent stress level? The ladder runs `on_fraction 0.50`; realistic is
~0.10 (`DESIGN.md` §0b). Re-simulate the worst in-scope scenario at the realistic level and check whether any
rank moves; a verdict that holds only at stress is a robustness finding. `sim/panel.py`,
`policy_benchmark.py`.

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
`priority: later · kind: hygiene · 2026-08-31; the review census 2026-09-14`
Open, each its own commit: the index's duplicate map as an alias map `dropped_t_id → kept_t_id` (an index
rebuild, no panel re-scan, checked with `rename_identity.py`; `reach` is covered by no other hash). The
2026-09-14 review census (91 % of `src/rigel` executed by the suite and every `--self-test`; C++ is
coverage-blind, so its functions were counted by name) cut what was dead and recorded the rest:
(a) ROTTEN BUT LIVE, each moving an instrument's or a toy's numbers when repaired: `toy_harness.harvest`
calibrates its donor undrained and without the two-pool contrast, so every toy inherits both;
`pass0_vs_oracle`'s C_input arms hand the post-capture truth law to geometry, and its C_info tables price the
retired length channel; `quant_accuracy`'s oracle arms are undrained (documented there);
`pipeline._DEFAULT_MEAN_FRAG = 200.0` is a magic fallback for an empty RNA length pool.
(b) VACUOUS GATE CLAUSES, never reached by their substrate: `test_transfer_faces.py`'s "nowhere when the flags
clear" and `_expected_level_rows`' `level_bound_row` branch; `test_transfer_policy.py`'s LEVEL-face invariants
(no LEVEL face on its toy); `test_transfer_rna_lanes.py`'s terminus clause and `own_level is None` clause.
(c) DUPLICATES across the kept instruments: the closure checks, `stratum` / `is_zero_gdna`, `_SCOPE`, `PANELS`,
`TYPE_NAME`, and the ladder's two paths hard-coded in five instruments.
(d) CLAIMS NOT RE-DERIVED on the current tree, left standing: that most in-scope error sits at the simplex
vertices (`simplex_logodds`, a relay-era measurement); `sweep`'s refused deferral of UNIDENTIFIED slots to the
prior (priced 2026-07, with no refusal entry here); `region_geometry`'s "no per-region spliced floor" A/B
(relay-era); and `fl`'s crossing pools called "gDNA by structure" because mature RNA never crosses an
exon|intron boundary, while RNA that has not spliced there does.
Kept as coverage GAPS, not dead code: the five CLI command bodies, the silent policy through `calibrate`, the
simulator's sharded writers and its whole-genome grid, and the zarr splice blacklist.

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
unprobed class +4.40 → +0.45. Refused with numbers, each open question tried each way (`DESIGN.md` §7.2): the
solve's own posterior (+1.25 / +4.96 nat unprobed), the clip outside the expectation (indistinguishable), the
own count alone (+4.5), the plug-in (+3.36), the apportionment iterated to convergence (identical on the test
chromosome, +0.44 against +0.46 on the ladder), a joint update of neighbours (never converges). The
falsification tests were verified failing on the shipped ruler through its own fixture: an unprobed exon with
no gDNA fragment read the floor at +3.19 nat over the depleted level, a transcript of 40 bp exons read 0 and
then the floor. The multimapper blindness the floor stood for is `ISSUES: multimapper-blind-support`.

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
the confusion matrix is re-derivable from `build_structural_claims` against `slot_truth.npz`.

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
