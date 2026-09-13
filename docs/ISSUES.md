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

### splice-out-premise-bias-uncorrected
`priority: later · kind: decision · 2026-09-02`
The splice-out message assumes spliced and unspliced fragments at one face share capture affinity; measured,
the premise fails as a bias under capture (`log a` ≈ 0 off, +0.28 under exon probes, +0.78 under junction
probes), so only a subtracted level — a cross-locale fudge, the owner's call — would correct it. No instrument
yet.

### gdna-landscape-trains-on-false-positives
`priority: next · kind: question · 2026-09-02; the defect closed 2026-09-10 (`DESIGN.md` §7.1)`
The population rule and the E-step landed (`ISSUES: the-landscape-training-population-arms`; the zero controls
~150k → a few hundred). Open: (a) under capture a short dark region's mass is split over both modes of the
previous fit (deferred 1.01–1.02×, stranded ON 1.004×); (b) the zero-RNA controls move ±4 % (`g05 ss.99 ON`
+4.4 %, `g98 ss.50 OFF` −1.4 %); (c) the Beta(½,½) reference still decides a blind slot. Refused here:
excluding κ-dead exons (`g50 ss.50 ON` 2,691 → 56,422), AMBIG in the final fit (worse 25/32).
`landscape_training_census.py`.

### drain-contaminates-certified-rna
`priority: later (parked by the owner, 2026-09-01) · kind: defect · 2026-08-31`
The second pass deposits some true-gDNA fragments into the certified-RNA banks: 233 records at
`g50 ss.99 OFF`, 1,482 at `g98 ss.99 ON` (1.9 % of that in-scope channel). The leak is exactly posterior
sampling (realized = Σ P(spliced) within 0.0–1.4σ) and a no-leak counterfactual moves the 0.8.0 metric
−0.60 %/−0.05 % at worst, so the harm is the certainty claim and no in-solve correction is licensed
(`ISSUES: drain-provenance-split`). To build: `DrainQC` records `Σ q_null`; the middle-bin posterior bias
repaired as its own A/B. `calibration_vs_oracle.py`.

### measured-prior-rung-4
`priority: now · kind: build · 2026-08-26`
ψ's composition reference fitted from composition-free observables (`DESIGN.md` §3.1a-i, `EQUATIONS.md` §2.3b),
consuming `rho_0`, the per-class enrichment responsibility `w` and the enrichment detector as a boolean —
never `span_R` (`TRAPS: a-mode-count-is-not-a-well-posed-quantity`). Requirements: span both lattice ends;
`rho_bg` on-rate for the destination's class (the intergenic pool is the unprobed rate); exact at `g00`; one
pseudo-fragment; priced with the messages against a shuffle. Not a better tilt
(`ISSUES: reference-prior-refuted-at-concept-level`, `ISSUES: data-derived-reference-location`).
The census that surveyed the field per condition was retired 2026-09-13; the fit is
`calibration.abundance_landscape` (gated by `tests/calibration/test_abundance_landscape.py`) and its measured
facts are in `DESIGN.md` §3.1a-iii.

### reference-prior-refuted-at-concept-level
`priority: now · kind: design-constraint · 2026-08-24`
ψ's reference tilt is refuted at the concept level: on λ = logit(f_g) the data's information
`I ∝ N_eff·disc·[f_g(1−f_g)]²` is zero at κ = ½ while the tilt holds fixed nats there — overturned at N = 3
where the strand channel is alive, never at κ = ½ (0.7471 from N = 10 to 10⁶); 82–95 % of unstranded pass-0
error sits on slots it decides. Refused: another constant (0.75 optimal on none of 16), `σ(L)` (3.3×/10.8×
worse at the zero controls), re-weighting (up to 2.8× worse), an information-weighted tilt, the arcsine
coordinate (`ISSUES: arcsine-magnitude-coordinate`). Candidate: extend the density channel
(`density_lambda_factor`).

### the-lower-bound-noise-ratchet
`priority: later (with the enrichment witness) · kind: defect · 2026-09-05`
A level from an RNA-rich node's own strand profile has a mode that is noise around zero; its lower side bounds
its neighbours and the tightest noisy neighbour wins. Ladder `g05 ss.99 OFF` 44,714 → 45,076; test chromosome
`g00 ss.99 ON` 64 → 216, 187 of it at `capcluster_ab`'s inner termini with no gDNA. A two-sided own-profile
level keeps fewer fragments overall (−0.9 %/+1.5 % vs −4.7 %/−0.5 %), so lower-only stays; the cure is the
enrichment witness `ISSUES: two-sided-exon-row` waits for. `zero_controls.py`.

### flux-floor-dispersion
`priority: later (with the transport-dispersion decomposition) · kind: question · 2026-09-08`
The certified flux at a junction is a lower-sided estimate of the exon's strand RNA level priced by the pair
(`DESIGN.md` §6b.13); the route rate scatters beyond counting (median −3 %, 5–9 % at depth; 0–40 % over on
nine block readings) and at a pure-RNA exon the price cannot see it, so a lucky over-read is a sharp floor a
few points too high. `transport_dispersion.py`.

### flux-price-witness-units
`priority: next · kind: problem · 2026-09-09`
The flux price compares the junction's route rate (whole-strand units) with the exon's column count; a column
holds `(1 − κ)` of the strand's RNA, so every flux floor pays `log(1 − κ)²` nats² that is no disagreement
(0.14 at κ = 0.31, 0.48 unstranded). Record: the golden `strand_ss65_multi_iso`'s nested exon (no gDNA) reads
0.152 gDNA from a ceiling of `f_g ≤ 0.38` where the flux says ≤ 0. Needed: the strand's RNA count as witness
at single-strand exons and a bounded one at both-stranded exons (`ISSUES: flux-witness-in-strand-units` is the
naive form). `policy_prototype.py`.

### ambig-node-as-a-gdna-source
`priority: after phase 2 · kind: decision, measured once · 2026-09-08`
A both-stranded node's own counts plus its held RNA levels, read as a gDNA level and emitted lower-sided —
refused as built: `g05 ss.70 OFF` +2.6 % on all three test panels, +4…+38 % on the weak-κ zero controls,
ladder `g05 ss.50 OFF` 1.745× / `g05 ss.99 ON` 1.465× (at low gDNA the level's mode is noise and travels as a
floor); its target, the walled host exon of a `span` locus, gains ~100 fragments. Re-open with a gate on the
emitted level's own width once phase 2's ceiling is in place; `policy_prototype.py`.

### message-layer-open-cases
`priority: next · kind: question · 2026-09-09`
Four residual cases, none a hole (`DESIGN.md` §6b.4–§6b.14): (a) an exon with both faces speaking — the
arrivals are summed (`_fuse`) and the price where both carry the same intron's claim is unmeasured; (b) the
intron factory on a region with both exon and intron bits, under capture; (c) the chain of termini — an empty
outside piece (median 12 bp) whose far face is another terminus, half the ladder's terminus-boundary error,
every upper side refused (`ISSUES: the-edge-upper-side`); (d) substrate `nest` (the prior already serves it,
0.666 vs 0.630), `div`, the antisense's nascent variant (`docs/TESTING.md` §0a).
`policy_benchmark.py --by-class`.

### two-sided-exon-row
`priority: later (deferred by the owner 2026-09-13 to the calibration-accuracy thread) · kind: problem · 2026-09-04`
Today (the one-lattice tree): the toy harness gate reads a factor of 88 — exon |Δf_g| 0.762 beside a pure-gDNA
intron, 0.0086 beside a nascent-bearing one (the lattice fixed the wet arm and left the mechanism); on the
ladder the `R exon (licensed intron face)` class is 4–6 % of the unstranded error with `transfer` at parity or
better than `silent` there, while on the test chromosome it is 20–32 % and transfer is worse (2,647 → 4,056 at
`g50 ss.50 OFF`). Not the pre-port thread's; the first step when taken up is the witness's derivation.
On unstranded data an exon's held row is the intron's composition through the face map, whose upper side is
the map's plateau above its ceiling — a lower bound on gDNA, so at pass zero an unstranded licensed exon reads
~9× its true gDNA (+2,000 % at `g25 ss.50 OFF`) and forwarding it compounds the bias (+10 % on the in-scope
unstranded row; `policy_prototype.py --by-class`). The plateau is honest: every cap, two-sided row or wall is
refused where junction probes enrich the flux more than the crossing
(`ISSUES: the-certified-flux-row-as-a-level`, `ISSUES: two-sided-exon-row-forms`,
`ISSUES: the-discrepancy-priced-cap`, `ISSUES: the-two-sided-level-lane`,
`ISSUES: the-wall-above-the-face-map-ceiling`). What would license a two-sided level is whether this library's
gDNA is enriched, witnessed on unstranded data by silent genes' exons against the intergenic density — ruled
the landscape prior's job (2026-09-06, `ISSUES: gdna-landscape-trains-on-false-positives`). Judge at
`R exon (licensed)` and the walled classes, halves apart, pass zero beside the pipeline; the toy harness gate
`test_the_harness_REPRODUCES_the_intron_composition_dependence` is a strict xfail citing this entry.

### per-transcript-prior-lane
`priority: next · kind: build · 2026-08-31`
`rna_prior_weight` is built end to end but `pipeline.py` omits it, so the shipped EM carries no per-transcript
information; a perfect per-transcript prior roughly halves in-scope gene-level error. Two weightings are
refused (`ISSUES: refused-transcript-weights`, `ISSUES: refused-soft-min-path-weighting`): the support problem
is the whole problem, so the next candidate is a sparsity mechanism, targeting expressed multi-exon
transcripts with median exon ≤ 150 bp (`ISSUES: the-rna-length-law-fix`). `quant_accuracy.py`.

### background-abundance-pair-unruled
`priority: later · kind: decision · 2026-09-13 (from the tunables census, W6)`
`CalibrationConfig.background_abundance` chooses which (counts, exposure) pair the pooled gDNA background
takes: `"contained"` (the count over the gDNA contained effective length — unbiased, the fragment-length pmf
in the divisor) or `"measured_total"` (the START/END banks over the region's own length — pmf-free, refuses
without the wall inputs). The two agree off capture and part under it, where the contained divisor
under-reads the true gDNA rate several-fold while the pmf-free pair over-reads on pools carrying nascent
RNA. That is a design decision, not a tunable: rule which pair ships (`total_abundance_audit.py` scores
them; `calibration_vs_oracle.py --set calibration.background_abundance=measured_total` prices the swap on
the metric) and the field goes with the ruling. Kept through W6 for that reason alone.

### antisense-prior-assembly-casualty
`priority: the prior-assembly session · kind: decision · 2026-08-18 (named 2026-09-13)`
`tests/scenarios/test_antisense_intronic.py::test_nrna_multiexon_t2_low_ss` is a strict xfail: `assemble_priors`
pins synthetic nascent RNA at Dirichlet alpha = 0 (`EQUATIONS.md` §9b), so recovered RNA lands on the annotated
antisense t2 — 72 today against the test's limit of 50 (80 under the relay), while both pool totals are better
with messages on. The owner plans to change the alpha = 0 rule in the post-calibration prior-assembly session;
the xfail is the executable record of exactly that pending change and closes there, with a test that asserts
the new rule's promise. Not the pre-port thread's.

### theta-quadrature-at-zero-gdna
`priority: later · kind: defect · 2026-09-13`
ψ marginalises the AMBIG cube over θ by a plain sum on a uniform θ lattice (`simplex_logodds._tilt_grid`). At
strand purity the strand likelihood is flat in θ to first order, so a pure-RNA AMBIG slot's θ-posterior is a
quartic peak at ±π/2 of width ∝ n^−¼, and where the lattice does not resolve it the λ-marginal is biased toward
gDNA. On the ladder's `g00 ss.99 capture ON` row, where nothing cancels a false positive, `sweep_n_tilt` 30
reads 9,488 fragments against 262 at 60 and 15 reads worse, with the message layer off as well (ψ's own
quadrature, not the lanes' cube rows); 120 and 240 equal 60. A Chebyshev mesh clustered at the ends with
midpoint weights is REFUSED (30 nodes 9.9×, and it breaks a `g00` row the uniform mesh holds: the coarse middle
hurts too). The per-slot bias is small (≤ 4 fragments at 50k, prior-free, the variance frozen at purity); the
refits amplify it at zero gDNA. Latent at 60 on a deep real library's heavy pure-RNA AMBIG exons (a 0.03 rad
peak at 500k fragments against a 0.053 rad step). The fix is a θ quadrature whose accuracy does not depend on
the node count, not a larger K_t. Derived 2026-09-13 (`docs/dev/THETA_QUADRATURE.md`): at fixed λ the strand
term is an exact Gaussian in τ whose peak narrows as √f_g beyond the strand-pure boundary, so the bias is the
λ-dependence of the lattice's resolution error. REFUSED with its number: the trapezoid endpoint weights (the
lattice's first-order term) — ≤ 0.7 fragments on any `g00` row, the K_t 30 failure 9,821 → 9,821.
`calibration_vs_oracle.py --set calibration.sweep_n_tilt=30`, the `g00` rows.

### performance-memory-bounded-solve
`priority: now · kind: build · 2026-08-17 (mandatory before 0.8.0), re-framed 2026-09-11, the decomposition landed 2026-09-11`
Calibration is the tool's one unfinished component: on 18.6M fragments its four sweeps are 706 s of an
836 s run and hold the 32 GB peak, on ONE core, while the locus EM beside it takes 8 s; the cost is set by
the index's 2.09M chain slots, not by depth. THE DECOMPOSITION IS THE LOCUS and it is BUILT
(`DESIGN.md` §6b.15): a terminal — a region admitting no RNA strand — receives nothing, structurally, so
the sweep solves the chain a locus block at a time (`sweep.solve_chain`, `region_chain.locus_blocks`),
the only cross-block information being the policy's library reduction; ψ's read-out is chunk-exact, so
the block size (`CalibrationConfig.sweep_block_slots`) moves no number and is a working-set knob. What is
LEFT, and why the serial floor of about 800× is not yet cashed: the block solve is Python, and the two
passes and the policy's `prepare` — most of a sweep — hold the GIL, so threads cannot speed them
(measured 0.83–0.94× at 8 threads; processes 6.2×, deferred). The owner's decision (2026-09-11): the
parallel executor waits for the C/C++ port of `_solve_block`, which lands on this structure; until then the
sweeps run serially and one block at a time. Also still genome-wide: the memoised intron-factory rows
(`(n_slots, K)`, 1–2 GB per grid), which a block could build for itself. `profiling/profiler.py`,
`profiling/sweep_replay.py --block-slots`.
THE AGREED ORDER (owner, 2026-09-11; `ROADMAP.md` rank 1 carries it without numbers): ⓪ re-measure the
deep library end to end, `main` against the landed tree, two back-to-back pairs at 8 threads — DONE
(`~/Downloads/rigel_runs/perf/ab_locus_2026-09-11/`, `profiler.py --compare`): wall 917 → 892 s and
892 → 875 s (0.97, 0.98), PEAK 33.2 → 14.9 GB and 32.6 → 15.0 GB, the four sweeps 769 → 738 s and
747 → 727 s (0.96, 0.97), untouched stages at 1.00. The peak is no longer the sweep's (11.5 GB while it
runs) but `build_region_geometry`'s transient (14.5 GB) and the pre-sweep `init_beliefs` solve — the
next memory target, after the sweeps' time; the final ψ lost ~7 s to smaller tiles inside 5,000-slot
blocks (`_block_rows` inside a block), a note for ① and ③; ① DONE, and larger than scoped: the WHOLE message layer — `prepare`, both passes and the
policy's solve — is refit-invariant given the grid (`DESIGN.md` §6b.15), so the refit sweeps are served
their messages from a content-keyed `message_cache.MessageCache` and pay only their two ψ solves. Measured on
the deep library, two back-to-back pairs at 8 threads, cache off → on
(`~/Downloads/rigel_runs/perf/ab_memo_2026-09-11/`): the refit grid is stable (`n_grid` 138, L 23.18 for
all three refits), refit 1 misses its 426 blocks and refits 2–3 hit all 426 — 38 s each instead of
176 s; wall 795 → 515 s and 799 → 515 s (0.65), `calibrate` 704 → 424 s, `prepare`/passes/policy solve
0.49/0.48/0.45, untouched stages 1.00; the cache holds 2.68 GB (its cube rows as the float32 the AMBIG
solve casts them to — 4.1 GB as float64), so the peak rose 15.0 → 17.8 GB. Whether a memory-constrained
run should be able to switch it off is a tunable for the owner to rule on; ② DONE — `messages.faces.Faces`:
the rules as ``(n, 2)`` typed tables over (destination, side) with five kinds and a row store, the lanes'
faces as bits, bit-identical (`DESIGN.md` §6b.15); ③ the port
of `_solve_block` — passes and `transfer_rows`, `prepare`, ψ, then threads over blocks — behind a derived
tolerance gate (promote the tolerant replay comparator into `sweep_replay.py`); ④ the factory rows per
block; ⑤ the scan (`ISSUES: scan-thread-split-starves-the-workers`) and the second pass. Two things not
to do: micro-optimise the Python passes (a silent-hop early exit halves them and the port deletes it),
and bake the λ lattice into the port — its step (`sweep_logodds_step`) is a parameter (`DESIGN.md` §6b.15).
MEMORY (plan step D / worklist W4) DONE 2026-09-12, measured first: the peak was not the sweeps' but two
Python transients — `crossing_eff_length`'s `(404k sj × fragment lengths)` matrix chain, ~9 GB that the RSS
never gave back (macOS keeps freed arenas), raising the whole run's plateau, and `fit_landscape`'s
`(training regions × grid)` kernel matrices, +3 GB at the true peak inside the refits. Three landings, each
a numeric no-op on the metric (worst relative move of any stratum's abs_err 3.4e-15): the crossing divisor
in closed form over the pmf's cumulative sums (`O(objects)`, 8 ms and 34 MB at the human sj count against
1.5 s and 9 GB; the matrix form survives as the brute force its gate compares with); the landscape's
kernels built and summed a row tile at a time (`_render`, 200 MB at a million regions); the factory rows
as `calibrate.FactoryRows`, built per block as the sweep asks and never held for the chain (④ above,
2–5 GB freed). One back-to-back pair on the deep library: PEAK 19.2 → 11.4 GB, wall 505 → 498 s,
`region geometry` 1.44 → 0.20 s, `landscape fit` 6.3 → 4.3 s, untouched stages 1.00. What remains of the
11.4 GB is the pre-calibration floor (~5.6 GB: the index and the payload) plus the message cache's ~4 GB
(plan step E, the owner's switch) plus the sweep's own working set.

### scan-thread-split-starves-the-workers
`priority: next · kind: decision · 2026-09-11`
`BamScanConfig.resolved_scan_threads` gives BGZF decompression `min(4, total - 1)` threads and the scan
workers what is left, so a 2-4 thread budget runs ONE worker and decompression is not the bottleneck.
Scan seconds on the 18.6M-fragment library by (bgzf, workers): total 4 — (3,1) 113.5, (2,2) 59.5, (1,3)
42.4, (0,4) 33.7; total 8 — (4,4) 34.4, (2,6) 26.4, (1,7) 24.8; total 16 — (4,12) 19.1, (1,15) 19.5,
(2,14) 17.4. Any new split rule is a tunable and `--scan-bgzf-threads` is a user-facing flag, so the rule
is the owner's to set. `profiling/profiler.py --scan-only`.

### u-ruler-arm
`priority: next · kind: measurement · 2026-08`
Price the `U` ruler (the oracle's gDNA total at uniform density) end to end: a perfect-composition ruler is a
~2× loss on the two capture-OFF in-scope strata, where the correct factor is exactly 1.000 and `U` reads it
with no fitting. Capture-OFF only; read `ruler_n_moved`, never the aggregate. `calibration_vs_oracle.py`
carries the column.

### g00-shrinkage-upstream-repair
`priority: now · kind: defect · re-priced 2026-09-10`
At the zero-gDNA control the effective-length shrinkage contracts most transcripts where the factor should be
1.000: `calibration_vs_oracle.py` reads 0.150 against the oracle's 1.000 on the `g00` rows (`Σ|Δ len|` 996M
bp), unchanged after the composition repair (828 invented fragments of 18M), because
`capture_eff_length._global_reference_density` detects a reference from any five slots with positive mass. The
repair is the detector — does the library carry an enriched gDNA mode at all, a boolean
(`TRAPS: a-total-density-ratio`); `priors.py` shares the function. Decide first whether the "exactly 1.000 off
capture" contract line is stale (the in-scope capture-OFF strata read P 0.954/0.967 against O 0.923/0.926).

### capture-blind-gdna-divisor
`priority: next · kind: defect · 2026-08-31`
`gdna_opportunity_from_index` is computed from the index alone, so under capture it removes ~6 bp of a ~30 bp
length selection — the gDNA control of the retired `fl_anchor_gap.py` moved +6.0 % on all six capture-ON rows (gDNA has no
introns to miss), and with `ISSUES: eb-shrinkage-magic-ess` it owns the −5.90 % capture-ON length ceiling.
`capture_eff_length` already models the panel; it also blocks `ISSUES: crossing-pool-contrast`.

### eb-shrinkage-magic-ess
`priority: next · kind: defect · 2026-08-31`
`POOL_EB_PRIOR_ESS = 1000.0` shrinks the gDNA pmf toward `global_pmf` (mostly RNA whenever gDNA is a minority)
at a magic ESS: inert on the ladder (0.01 bp), dominant on the fl-gap arm at `g05` capture-ON (`ship−pool`
−23.7 of −31.7 bp). Replacement: reconcile the pools by their precision (`EQUATIONS.md` §6c).
`fl_pool_purity.py`.

### refit-vs-message-arbitration
`priority: next · kind: design · 2026-08`
At the unstranded × capture-OFF exon cell the refitted gDNA prior and the message both impute one slot with
nothing arbitrating them; the message is the accurate voice there and the refit displaces it. Belongs with
`ISSUES: measured-prior-rung-4`. `calibration_walk.py` (the refit rung against the message rung).

### prior-fidelity-vs-deliverable
`priority: next · kind: question · 2026-08`
Why is prior fidelity anti-correlated with deliverable quality? Leading answer: at the worst slots the
self-solve with the fitted prior is nearly right and the messages destroy it (measured at a retired rung;
confirm on a second stratum). Exclude the ruler first (`ISSUES: u-ruler-arm`,
`ISSUES: g00-shrinkage-upstream-repair`). `prior_vs_oracle.py`.

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

### psi-lambda-bracket-unshipped
`priority: later · kind: decision · 2026-08`
ψ's λ bracket was too narrow to express its own prior; `DensityLandscape.required_logodds_window` derives it
with no constant. Built, gated, priced (nearly every in-scope condition improves); ships OFF pending memory at
genome scale (a multiple on the lattice's point count) and the thermometer. Re-derive as a `policy_prototype.py --module`
arm.

### transfer-variance-premise
`priority: later · kind: question · 2026-08`
Does a hop's transfer variance price a ratio built on a handful of counts? The policy prices every hop by both
witnesses' counting plus the pair's disagreement (`transfer_rows.hop_price`); whether that is right where a
pair agrees by coincidence is the open half (the landscape is no substitute: ~10× over-stated). `EQUATIONS.md`
§3.5h.

### nascent-stress-sensitivity
`priority: later · kind: question · 2026-08-22`
Does any in-scope verdict depend on the nascent stress level? The ladder runs `on_fraction 0.50`; realistic is
~0.10 (`DESIGN.md` §0b). Re-simulate the worst in-scope scenario at the realistic level and check whether any
rank moves; a verdict that holds only at stress is a robustness finding. `sim/panel.py`,
`policy_benchmark.py`.

### hygiene-ledger
`priority: later · kind: hygiene · 2026-08-31`
Each its own commit, none moving the 0.8.0 metric:
- the wave-3 frame migration, RULED 2026-09-13 (owner): of the six bank-readers still on pass one, four are
  RETIRED rather than migrated — their questions closed, superseded or parked, none loaded by a test, none in
  any procedure — and two migrate to the drained frame the truth is certified in (`fl_pool_purity` and
  `transport_dispersion`, both migrated 2026-09-13). Retired, each its own commit, all in git: `structural_claims_audit` (stage 0 came
  out confusion-matrix clean; only a parked, twice-refused issue cited it; its inputs are frame-invariant);
  `gdna_pool_census` (Stage A is closed and gated by the accumulator's spec test; it overlapped `fl_pool_purity`);
  `calibration_truth_ab` (its one number is in every `calibration_vs_oracle.py` row's `pools`, in the drained
  frame with both zero controls; the drain it priced shipped; its length ceiling is `em_fl_ceiling.py`'s,
  through the EM); `abundance_landscape_census` (the fit it surveys is QC-only — nothing in the solve reads
  it — its per-condition survey is recorded, and the rung that cited it is not on the roadmap).
- the index's duplicate map as an alias map `dropped_t_id → kept_t_id` (an index rebuild, no panel re-scan —
  verify with `rescan_panels.py`; `reach` is covered by no other hash);

### oracle-effective-length-diagnostic
`priority: later · kind: measurement · 2026-08`
Started, not finished; it re-ranks the two ruler issues. Needs a stashed pre-closure arm or it measures
nothing — one arm is not an A/B. No instrument yet.

### the-cancelling-pair
`priority: parked (refused twice) · kind: design · 2026-08-26`
`struct_lock` rescoped to `g1_locked ∧ REGION` and the `intergenic|exon` boundary claiming its
RNA-contaminated crossing mass as gDNA: neither half prices alone (`TRAPS: a-cancelling-defect-pair`); five
xfails go green iff the pair lands. Re-priced 2026-08-26 with the measured intron reference as load: still
refused (worse on two of three in-scope strata, wins confined to `g00`). The one-sided certified-RNA bound is
the only mechanism the zero control has endorsed on every row (−81.9 %, 8/8) and is panel-negative alone,
because the level channel it was covering for is structurally disconnected (bipartite chain; only pure-gDNA
REGIONs originate a level). Revive only with messages on at `g05 ss0.50 capture_on`. The structural-claims
audit that scored it was retired 2026-09-13; the confusion matrix is re-derivable from `build_structural_claims`
against `calibration_oracle.py`'s `slot_truth.npz`.

### crossing-pool-contrast
`priority: parked (blocked) · kind: question · 2026-08-31`
A second gDNA length contrast on the crossing pools: with oracle weights it beats the contained one under
capture (TV 0.076 vs 0.136; 0.078 vs 0.182) and starves off capture (pool 3 at 29–630 fragments). Blocked: the
weight estimator does not transfer under capture (`ISSUES: capture-blind-gdna-divisor`), and no shadow
transcript overlaps a gene edge, so pool 3 reads exactly 1.0000 pure
(`TRAPS: purity-is-a-property-of-the-annotation`). The two crossing enrichments are equal (1.002/1.034), so
only the common level is missing.

### capture-degeneracy-standing-risk
`priority: parked (watch) · kind: risk · 2026-08-31`
The gDNA two-pool contrast survives capture by a degeneracy: the shared-contaminant assumption is false under
capture (TV 0.95 vs 0.06–0.14 off) and it is safe only because the intergenic pool is depleted-not-impure, so
`a_0` clips to 1 and the algebra collapses to `g = f_0` (to 3e-17). A probe panel that put RNA into intergenic
space would break it silently; `_deconvolved_gdna_counts` carries the derivation. No panel can fire it today.

### pure-rna-mirror-asymmetry
`priority: parked · kind: defect · 2026-08`
Two exact per-fragment mirrors of a pure-RNA library deconvolve differently in `count_gdna_region` by a few
percent, neither boundary-only nor monotone in strandedness. An R1-sense library is simulable; no instrument
yet.

### parked-capture-pilot-sign
`priority: parked · kind: question · 2026-08-13`
Two capture-ON pilot rows disagreed about the sign of every length correction; both panels were deleted and
the correction lives in the retired length channel. If revisited, find which row is lying rather than
averaging; the fl-gap panels are not a drop-in (`ISSUES: flgap-panels-stale-nascent-model`).

---

## CLOSED / REFUSED — do not rebuild these; append-only

### f32-strand-tilt-at-half
CLOSED by landing 2026-09-12: the AMBIG cube is float64 like the rest of ψ — the float32 cube was a memory
choice the tiling made moot, and one solver (`simplex_logodds._solve_logodds`, a single-strand slot the
``K_t = 1`` case) has one precision. At κ = ½ the strand mean is ½ identically in float64 and `w_pos` reads ½.

Every entry keeps its stamped measurement exactly as recorded: a graveyard row without its number is an
invitation to rebuild. A row measured on "all 36 conditions" or quoting `g01`/`g10`/`g25`/`g75`/`g90` predates
the ladder retired 2026-08-13 — the verdict stands as a record, and re-opening one means re-running it on the
current panel. Where a mechanism's only target was unstranded × capture-ON the row is moot as a 0.8.0
candidate on top of being refused; the `g00` zero-control column is never moot.

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
CLOSED by landing 2026-09-13 (`DESIGN.md` §6b.15, W5 the grid study): one λ lattice for every consumer,
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
