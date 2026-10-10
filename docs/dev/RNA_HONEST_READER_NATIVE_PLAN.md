# The native honest capture reader: contract, gates and stages

*2026-10-10. The owner's go: build the native honest reader. Sandbox plan; nothing here is authoritative.
The derivation is `EQUATIONS`' proper local capture reference and the review's fixed-lattice contract
(`RNA_CAPTURE_EXTERNAL_REVIEW.md` §6); the A/B that selected the design is that review's benchmark; the
reason no cheaper reader ships is `RNA_CAPTURE_REFERENCE_INVESTIGATION.md`. This note is the PLAN and the
PROTOTYPE stage's record. Source enters only after the gates below pass.*

## 1. What the reader is

For every slot of the chain (every region and every boundary), the posterior mode of the gDNA density
under the fitted landscape, from the slot's own observations and its delivered local factors:

```
log L(x) = log ∫ Pois(u; ρEg/2 + q·r) · Pois(v; ρEg/2 + (1−q)·r) · C(log(ρEg/r)) · H_pos(r·f_pos/Er)
                 · H_neg(r·f_neg/Er) · r^(−1/2) dr  + D(log ρ)            with ρ = e^x
w(slot) = exp(x*),  x* = argmax_x [ log L(x) + logP_landscape(x) ]      (parabola-refined)
```

`u, v` the slot's unspliced strand columns, `Eg, Er` its gDNA and RNA opportunities, `q` the protocol's
column probability for its RNA strand, `C` the delivered composition row on the λ lattice, `D` the held DNA
level at `ρ`, `H_pos/H_neg` the held RNA levels at their own amounts, all with held ends. A both-strand slot
integrates the strand share over the arcsine continuum plus the two pure atoms admitted by its witness
bits (`−log 3` normalisation); a slot with no RNA opportunity is the pure-DNA limit with the composition at
its all-gDNA end. No reference density, no background, no located-mode test, no clip, no `None`.

The weights are relative (a common factor is invisible to the EM, TPM uses the plain length). The
consumer is the shipped conservation sum with the junction capped at the largest weight among the objects a
junction fragment touches (the two pieces and the two boundaries).

## 2. The lattices, verbatim from the prototype (`reader_proto.py`, archived)

- Density `x`: the landscape's own grid (260 nodes, `CONSTANTS.landscape.grid_points`), restricted to the
  nodes where `logP > max logP − T`, `T = −log ε₆₄` (ψ's window rule): outside it the prior is numerically
  zero and no posterior mode can sit there.
- RNA amount `y = log r`: the union of a 0.2-nat lattice (`CalibrationConfig.sweep_logodds_step`) on
  `[y_hi − 20, y_hi]`, `y_hi = log(n + 12√n + 12)`, with a window of `⌈√(2T)⌉` nodes on either side at step
  `h = min(0.2, 1/√(n+1))` around the own term's mode `y*` (the right-hand zero of `r·d/dr log f`, by
  bisection) and the same window around the full integrand's maximum on the first lattice; trapezoid
  weights from the node spacings, every node clipped to `[y_lo − 5, y_hi + 1]`.
- Strand share for a both-strand slot: 24 midpoint nodes in `φ ∈ (0, π/2)` with `s = cos²φ` and the
  arcsine measure `2dφ/π`, plus a window at step `sd_s = 1/(√(n+1)·max(|2q−1|, 1/√(n+1)))` around the
  observed share `s* = (u/n − (1−q))/(2q−1)` when `|2q−1|√(n+1) > 1`, plus a second window around the
  density nodes' maxima clustered at `sd_s`; the endpoints `0, π/2` included; trapezoid in `φ`.
- Readout: the lattice argmax of `log L + logP`, refined by the parabola through its two neighbours.

Every constant above is a derivation already in the tree (ψ's `K_t = 24`, the 0.2 lattice, `T`); none is
new. The step-halving check of the review (every 40th slot and six deep both-strand slots at half the steps)
is the convergence receipt and stays in the gate suite, not in production.

## 3. Gates (each written and failing first; each deliberate defect must fire a named gate)

| gate | what it checks | reference |
|---|---|---|
| `frozen-objects` | the six frozen objects of `2026-10-09_capture_prior` reproduce the archived log weights to 1e-3 and the NumPy prototype's values to 1e-6 (same lattices) | `product_rule/results.json`, `results_map/*.npz` |
| `analytic-limits` | the analytic limits: flat evidence reads the prior's mode; pure DNA with 100 reads reads the Poisson mode; a slot with `Er = 0` reads the pure-DNA formula | derivation |
| `prototype-identity` | every slot of `benign g50 ss0.99 ON` and `g05 ss0.99 OFF` (3,518 each) matches the NumPy prototype's `map` readout to 1e-6, both-strand slots included; the ladder's `g50 ss0.99 ON` (70,176 slots) likewise | `results_map/*.npz` |
| `thread-identity` | thread count and block packing do not change a single value (bit identity across `n_threads ∈ {1, 4, 8}`) | `rename_identity.py --check` pattern |
| `solve-untouched` | the count outputs of the solve are untouched (the reader runs after the solve and must not perturb it): seven cached conditions bit-identical | `rename_identity.py` |
| `deliberate-defects` | the deliberate defects fire: wrong tail Jacobian in `y`, omitted `−log 3`, a dropped witness atom, the DNA level read at the wrong origin, the composition row read at `log(ρEg)` instead of `log(ρEg/r)`, the window rule off by one, `q` and `1−q` swapped, the parabola refinement disabled (`prototype-identity` tolerance) | `frozen-objects`, `analytic-limits`, `prototype-identity` |
| `speed` | SPEED: the reader stage on LBX0190 (2,087,476 slots) in at most 60 s on 8 cores, target 20 s; reported beside the sweep's own time | this plan §6 |
| `end-to-end-identity` | end to end, pinned and fractional: the kernel's weights through the harness reproduce the NumPy arm `land_map_cap4` on the benchmark's rows (g05–g25 test chromosome, single probe, strength 100, junction-probed, the ladder's six rows) to the pair's own variation | the review's §1.4 table |
| `real-libraries` | real libraries, one process at a time: the four cfRNA BAMs, gDNA fraction and pools against the NumPy arm's numbers where they exist (LBX0190: 0.0997 with the median readout) and the shipped reader's | the review's §1.3 |

## 4. Stages

1. **Stage 1 (this session): the standalone kernel.** A pybind function in the worktree's native module
   taking the per-slot factors as packed arrays (the same inputs the NumPy prototype consumes, produced by
   the transfer harness) and the landscape, returning `x*` per slot. Gates `frozen-objects`, `analytic-limits`, `prototype-identity`, `thread-identity`, `deliberate-defects`; `speed` measured on the
   kernel alone (the Python table extraction is outside the clock and is removed by stage 2).
2. **Stage 2: streaming inside the block solve.** The kernel moves into `solve_blocks` after the last
   refit, reading the block's received tables in place (composition rows, held DNA and RNA levels, cube
   rows) with no `(objects × density)` table anywhere; threaded as the sweep is. `solve-untouched` and `speed` proper.
3. **Stage 3: the consumers and the schema.** `transcript_capture_eff_lengths` and `assemble_priors` read
   `w` where they read `c̃` today; the touched-objects cap replaces the clip in `_cut_efficiencies`; the
   result carries `gdna_capture_weight_region/boundary` and no reference fields. `end-to-end-identity`, `real-libraries`.
4. **Stage 4: landing.** The deletions of the review's §4 (the located-mode census, `capture_efficiency.py`,
   the reference fields and validators, the CLI and report fields, the clip and the `None` return, the
   `np.minimum(eff_len, span)` clamp, their tests), the goldens re-derived and read, DESIGN §7.2 and
   EQUATIONS §11 moved, ISSUES closed by the move rule, the suite re-derived, ruff, `preflight.py --full`.
   The commit is the owner's.

## 5. Isolation

The worktree `~/proj/rigel-honest-reader` (detached at `55f5843d`), built with `pip wheel
--no-build-isolation --no-deps` into `~/Downloads/rigel_runs/prototypes/2026-10-10_honest_reader/`, unpacked
into its `site/`, imported through `use_proto.py` which drops the editable finder; `assert_proto()` proves
which build answered. The main install and its suite are untouched until stage 3.

## 6. The speed budget, derived

Under the lattices of §2 the evaluation count per slot is `N_x · N_y · N_s`: `N_x ≈ 100–140` (the landscape's
support), `N_y ≈ 20–30`, `N_s = 1` for a single-strand slot and `24–60` for a both-strand one. On the
ladder's genome-scale slot population scaled to the plasma library (2.09 M slots, 14 % both-strand) that is
about `3 × 10⁹` evaluations; at 20–70 ns each (one `exp`, two `log`s, two uniform-lattice interpolations) the
kernel is 60–200 s on one thread, 8–25 s on 8. Both-strand slots are over 80 % of it. The gate `speed` is set at
60 s on 8 cores with a 20 s target; if the full-grid evaluation misses it, the derived reduction is a
coarse-to-fine search over `x` (the 0.2-nat lattice plus the own term's analytic DNA mode as a candidate,
then one window per candidate basin), proven a numeric no-op on the mode before it is kept.

## 6b. Stage 4's scope, measured

The worktree's full suite under the stage-3 build: 24 failed, 2,790 passed (the main tree's baseline is
2,814 passed). The failure set, derived not eyeballed:

| gates fired | what they say |
|---|---|
| 20 of `tests/test_golden_output.py` | every golden scenario's output moves, capture-OFF ones included (every golden is capture-OFF): the honest reader reads each object's own gDNA density where the census read exactly 1. The magnitudes and the truth scoring are below |
| `test_constants.py::test_no_bare_numeric_constant_outside_config` | the zero-read search's stride (8) was a bare module constant in `calibrate.py`; it is now `CONSTANTS.calibration.reader_search_stride`, documented and range-checked |
| `test_landscape_training_population.py::…publishes_the_last_landscapes_located_enriched_mode` | the result's reference field no longer carries the located mode; the field and its test go with the census, replaced by the reader's own publication test |
| `test_quant_accuracy.py::…the_ruler_arm_hands_the_EM_the_lengths_it_names…` | the test assumed "on a capture-OFF toy the shipped lengths are the plain ones"; under the honest reader they are not (weights 0.74–0.92 on the toy), so the arm's expectation is rewritten against the plain `effective_length` column |
| `test_multimap_counting.py::TestParalogMultimapping::test_gdna_sweep[gdna_100]` | **a real defect, dissected below**; the test stays red and is the gate for its fix |

### The goldens, read (147 files; `tests/golden/` regenerated in the worktree and diffed against HEAD)

Every scenario's `em_effective_length` moves (the reader's weights are never exactly 1), by up to 82 % on one
transcript (`combo_extreme`); `gdna_eff_len_em` (the locus gDNA yield) contracts by up to 73 %; counts move by
at most 0.3 % on 13 scenarios and by 1–28 % on the rest (`combo_extreme` 28 %, `combo_moderate` 5 %,
`antisense_contained_ss90` 3 %, `gdna_heavy` 1.3 %, `gdna_light` 1.1 %). Scored against the simulator's own
truth (`sim.benchmark.run_benchmark` on every golden scenario, both builds, `prototypes/2026-10-10_honest_reader/golden_truth_*.json`):

| scenario (truth total) | Σ\|observed − truth\| main | honest reader |
|---|---|---|
| combo_extreme (54) | 6.72 | 22.05 |
| combo_moderate (142) | 5.03 | 9.42 |
| gdna_light (294) | 4.51 | 6.45 |
| antisense_contained_ss90 (1000) | 103.45 | 107.59 |
| multi_isoform_5tx (1000) | 38.58 | 43.48 |
| multi_isoform_3tx_skewed (1000) | 27.29 | 24.78 |
| multi_isoform_2tx_equal (1000) | 26.23 | 24.71 |
| all 21 scenarios | 298.07 | 323.79 |

Worse on 7, better on 5, unchanged on 9: on a capture-OFF library the honest reader's per-object weights are
Poisson noise around a common level, scaled by their maximum, and on a 10–12 kb scenario with a handful of
objects that noise is 10–25 % of the weight; it is the same cost the review's A.4 measured at genome scale
(benign stranded × OFF 6.88 → 8.74 %, ladder stranded × OFF 6.35 → 7.44 %), the price of reading every object
with no detector. It is not a defect of the kernel (the kernel is the prototype to 1e-9); it is the estimator.

### The paralog gate: the calibration count is blind to multimappers, and the reader reads the blindness as depletion

`TestParalogMultimapping.test_gdna_sweep[gdna_100]`: two sequence-identical 500 bp single-exon paralogs at gDNA
abundance 100, aligned (not oracle) reads. Dissected on both builds (`paralog_dissect.py`, logs beside it):

| | main (shipped reader) | honest reader |
|---|---|---|
| t1 / t2 counts (truth 75 / 71) | 154 / 0 | 123 / 79 |
| the two paralog exons' calibration gDNA count on a 155 bp contained support | 0 / 0 | 0 / 0 |
| their published capture weight | 1 / 1 | 4.8e-4 / 4.8e-4 |
| the flanking intergenic density | 0.23–0.31 fragments/bp | same |

Every fragment inside a paralog exon multimaps (NH = 2), and `bam_scanner.cpp` deposits into the calibration
accumulator only when `is_unique_mapper` (`deposit_to_accumulator` is skipped for a multimapper; the fragment
goes to the EM's buffer alone). So calibration sees 0 fragments on 155 bp beside 0.23/bp intergenic gDNA. The
census reader on a capture-OFF library said 1 regardless (the detector); the honest reader reads the zero
honestly: `exp(−ρ·Eg)` with `Eg = 155` at the background density is `e^{−36}`, the posterior mode is the
landscape's floor, the weight 1/2000. `assemble_priors` then contracts the locus's gDNA yield to its boundaries'
share, the EM reads a zero yield as "cannot emit", and the gDNA fragments at the paralogs are called RNA:
the total 154 → 202 against a truth of 146. The 123/79 "tie-break" is the boundaries' crossing-count noise
(26/17 against 21/17 fragments) amplified because the regions' weights are gone — not a legitimate break, so the
test's "delete this branch" instruction does not apply.

The pre-existing half: main's 154 against 146 is the same blindness one step earlier — the locus prior's gDNA
count at a multimapper-only region is 0, so the EM has no gDNA prior there. The census reader had the same
amplification on every CAPTURED library (a zero count against the reference reads depleted); the detector hid
it on capture-OFF ones. No benchmark number in the review saw any of this: every panel is an oracle BAM
(`NH:i:1` on every record sampled from `scenarios/` and the ladder), and the production default is
`BamScanConfig.include_multimap = True`.

**The footprint on real libraries** (`mm_footprint.py`: per region, uniquely mapped first-mate starts against
multimapper alignment starts, every alignment counted; the production index's 1,043,881 regions):

| library | unique fragments | multimapper alignments | regions with multimappers and NO unique fragment (exon regions) | regions where multimappers outnumber unique fragments, and the share of all multimapper alignments they hold |
|---|---|---|---|---|
| LBX0190 (plasma) | 138,403 | 104,140 | 5,488 (762) | 6,504 — 95 % |
| MO_3021 (plasma) | 772,744 | 583,407 | 23,443 (7,018) | 35,781 — 94 % |
| VCaP mix (deep) | 18,417,903 | 544,065 | 10,099 (5,305) | 13,792 — 72 % |

Multimappers are a large share of a plasma library (43 % of LBX0190's alignments) and they cluster in regions
where calibration counts almost nothing, which is exactly where the honest reader reads depletion and the locus
prior reads no gDNA. The production-path VCaP A/B (main against the worktree, pinned, truth 0.2518) prices it
end to end: [pending, running].

**The fix, derived (not built).** The accumulator ruling already says what multimappers are: "a LATER phase —
side buffer, then deterministic largest-remainder apportionment, integral always" — and the code never
implemented it for NH > 1 (the side buffer holds only gap-hypothesis fragments on ONE placement). The complete
fix is that ruling: a multimapper enters the side buffer with one hypothesis per alignment (ref, start, end,
introns), the second pass scores each placement with its existing score (`ρ(h)` from pass one's tally at the
placement's objects, the length law, the strand term), draws one, and the drain deposits it — then calibration
counts the paralog exons at their density, the reader reads them, and the locus prior holds their gDNA. It is a
format change to `DeferredRecords` (per-hypothesis coordinates), a scanner change (defer instead of skip) and a
drain change, each gated by the accumulator's executable specification, A/B'd on the aligned scenarios (the only
panels with multimappers) and on VCaP against its truth. The cheaper alternative — a per-region multimapper tally
in the scanner and a unique-deposit opportunity `Eg_r = S_r · U_r/(U_r + A_r/NH_r)` derived under uniform
genomic sampling — is a new geometry estimated from data and needs its own ruling. Neither is a reader change;
nothing in the reader can tell a withheld fragment from an absent one.

## 7. Status

**2026-10-10, stage 1.** The standalone kernel (`honest_reader.cpp`, module `_honest_reader_impl`, built in the
worktree in 30 s) reproduces the NumPy prototype's evaluator: on `benign g50 ss0.99 ON` every one of the
3,518 slots' modes agrees to 1.9e-10, both-strand slots included (`prototype-identity` passed), and the result is bit-identical
at 1 and 8 threads (`thread-identity` passed). One rule had to be made explicit: the prototype centres the own-term window
at the mode under the slot's `q` once per slot, not under each tilt node's assignment; the kernel carries both
rules and the replica uses the prototype's.

**A lattice rule learned.** Dropping the 0.2-nat coarse lattice at deep slots (keeping only the two windows)
moved 8 % of modes by more than 0.01 and one by 0.1, and the curves say why: at a deep slot whose delivered
composition row is a bound, the row is a cliff of hundreds of nats, the full integrand's peak sits at that
cliff several nats from the own term's mode, and a window around the mode alone misses it (log L low by up to
1,164 nats at the affected densities). The coarse lattice at the rows' own resolution is the detection net
that finds such peaks; it stays. Speed comes from elsewhere.

**The cost, measured.** The replica evaluates about 100 ns per lattice node in this first C++ (three
logarithms, an exponential and two interpolations per node): 1.5 ms per single-strand slot and about 160 ms
per both-strand slot (the tilt multiplies the whole y × x work by 24 or more). Extrapolated to LBX0190's
2.1 M slots with 14 % both-strand that is of order an hour and a half on 8 threads: far from `speed`. The program,
each step gated as a numeric no-op on the mode (1e-6 on every slot of the fixtures) before it is kept:

1. a coarse-to-fine search over the density grid (every fourth landscape node plus the own term's analytic
   DNA peak as candidates, then one window per local maximum of the coarse posterior), instead of all 260
   nodes: about 6×;
2. the s-independent shortcut for both-strand slots whose share integrand is constant (no strand reads and
   no RNA level: one evaluation instead of 24), the common case at genome scale;
3. per-node cost: the x-independent terms of the coarse lattice computed once per slot, O(1) interpolation
   on the uniform coarse nodes, the fast exponential of `fast_exp.h`: about 5×;
4. the detection scan of the far-peak net on the s-independent factors only (once per density node, not per
   tilt node) for both-strand slots: about 4×.

Target after 1–4: a both-strand slot near 1 ms, a single-strand slot near 20 µs, LBX0190's reader stage near
a minute on 8 threads, with the shortcut of step 2 deciding whether the 20 s target is reached.

**Progress, same day.** `prototype-identity` holds at genome scale: the ladder's `g50 ss0.99 ON` fixture (70,176 slots, 14 %
both-strand, the deepest 35,184 reads) matches the prototype to 1.2e-9 on every slot; the unoptimised
replica took 198 s on 8 threads there, which extrapolates to 98 min for LBX0190, the honest baseline.
Then, each kept only after the gate on every slot of the test fixture:

| step | mode against the full-grid replica | test chromosome, 8 threads |
|---|---|---|
| the replica (`prototype-identity` exact) | — | 5.5 s |
| 1. coarse-to-fine density search, stride 4 / 8, seeded with the slot's own DNA witnesses | max 8.5e-7 (the parabola's rounding), no slot above 1e-6 | 2.2 s / 1.6 s |
| 2. the share-independent shortcut (no strand reads, no RNA level) | exact by construction | (within the above) |
| every lattice node evaluated once: the RNA integral merges its three sorted lists, the tilt's second pass keeps the first pass's columns | exact (`prototype-identity` unchanged at 1.9e-10) | replica 2.7 s; search 8: 0.5 s |
| the derived tilt count (one node per standard deviation of the share, never above 24) | max 2.7e-3 at 12 shallow both-strand slots: a converged lattice change, not a no-op; kept behind its own flag and gate | 0.4 s |
| the share-independent caches on the coarse nodes (exp once per slot, the composition term once per density node) | [pending] | |

**The density search, refuted for slots with reads and kept for slots without.** On the ladder fixture the
coarse-to-fine search over the density grid missed the full-grid mode at 3 slots per 70,000 at stride 8, and
still at 77 slots at stride 2 (the λ lattice's own resolution), by up to 1.3 nats: at a slot with a few reads
on a few bases of opportunity the posterior is a plateau with several maxima within hundredths of a nat, and
no coarse sampling can reproduce which of them the full grid picks. The full grid is the contract for every
slot with a read. A zero-read slot's log likelihood is the log of a convolution of its factors with the
RNA-amount kernel, smooth at the kernel's width, so the search is exact there once its windows are a stride
wide at the prior's own peaks as well: with that rule the test fixture is bit-identical to the full grid and
the ladder fixture's differences are confined to the tilt skip below (no zero-read slot among them); those
slots are 43 % of the ladder fixture and 99 % of LBX0190.

**The tilt, measured and bounded.** Per slot class on the ladder fixture (CPU time, the replica lattices with
the exact reuse and caches):

| slot class | single-strand | both-strand |
|---|---|---|
| zero reads | 38 µs (the search) | (the share-independent shortcut) |
| 1–9 reads | 0.41 ms | 19 ms |
| 10–99 reads | 0.75 ms | 58 ms |
| 100–999 reads | 1.1 ms | 100 ms |
| ≥ 1,000 reads | 1.2 ms | 123 ms |

The both-strand cost is the 24 base share nodes, each a full RNA integral, twice. Two derived reductions,
each exact only where the share enters through the strand columns alone (a delivered RNA level is a function
of the share too and can lift a far node): skipping base nodes beyond `√(2T)·sd_s` of the observed share
(the window rule), and taking the node count from the share's own standard deviation (one node per `sd_s`,
never above 24). With RNA levels present both are unsafe (the unrestricted versions moved modes by up to 6.8
nats on the ladder) and the slot keeps the fixed count. Restricted to level-free slots, neither is a no-op
either: the skip moves 2 of 70,176 modes by up to 4.8e-3 (deep both-strand slots, where the DNA share of the
columns widens the share's scale beyond `sd_s`), the derived count 9 by up to the same. Both are converged
lattice changes at the level the review's step-halving check already accepts, and both stay behind flags,
off by default: the plasma libraries do not need them, and a deep library is where they would earn their
convergence gate.

**Stage 2 landed in the worktree, same day.** The evaluator moved into a header shared by the standalone
module and the block solve; `solve_blocks` takes an optional `reader` request (the mode array and the lattice
rules), and `solve_one_block` runs the reader over its owned slots after the final ψ, rebuilding each slot's
factors from the arena exactly as the harness did for the prototype: the composition row as the normalised
sum of the sides' compositions, the held DNA level as the sides' present gDNA profiles intersected at the gDNA
lane's origin, a single-strand slot's held RNA level from the level and flux rows of the sides that hold no
composition (admitted when its converted row is not flat), a both-strand slot's cube rows at their own
origin. Python runs it on the last refit only (`calibrate._solve`) and publishes the modes in the debug bundle
(stage 2 leaves the shipped reader in place; stage 3 switches the consumers). On `benign g50 ss0.99 ON` the
in-kernel modes equal the standalone kernel's to 1e-16 and the prototype's to 1.9e-10 on every slot, the
solve's own counts and the shipped efficiencies are bit-identical to the main install's (`solve-untouched`
passed), and the whole calibration with the reader takes 7.3 s on 8 threads. The other three fixtures agree
with the prototype to 5.3e-11, 2.2e-12 and (the ladder's 70,176 slots) 1.2e-9, every slot, the solve untouched
on each; the ladder condition's whole calibration with the reader takes 72 s on 8 threads, of which the
reader is about 60 s, as the per-class costs predicted. (A fixture file written by an early build of the
standalone kernel, before its own-term window rule was aligned, differs from the in-kernel modes by up to
1.8e-2 at both-strand slots; the prototype is the reference and the in-kernel reader matches it.)

**The plasma library, in place.** LBX0190 (2,087,476 slots): the whole calibration with the in-kernel reader
on its last refit takes 43 s on 8 threads, of which the reader is about 34 s (the main install's calibration
alone is about 9 s); every slot is read, 0.47 % of objects sit above 2 nats (the NumPy prototype read 0.49 %),
the largest at +16.9, where the shipped census locates no reference at all. Against the prototype's 4.8 hours
that is a 500× reduction, inside the `speed` gate's 60 s with the 20 s target in reach of the derived
reductions still behind flags. MO_3021 (194,341 slots with reads): 105 s for the whole calibration on 8
threads, the reader about 85 s, every slot read, 4.1 % of objects above 2 nats (its gDNA is concentrated after
all, which no mass-based reference could see: §2 of the reference investigation read it at 1.0× the median),
the largest at +10.2; the census again locates nothing. The `speed` gate is passed on the sparser plasma
library and missed by 25 s on the deeper one; the shallow both-strand slots that carry MO_3021's cost are
exactly what the two flagged reductions address.

**Stage 3 built in the worktree, same day.** `calibrate` publishes the reader's weights as the result's
efficiencies, scaled by their maximum so the consumers' unit stays 1 (the EM is invariant to the scale), with
the unit recorded as the reference density and the count of read slots as the members; a library whose last
refit fits no landscape carries no reader and every efficiency is 1, as before; the located-mode census is no
longer called. The junction's conservation sum is capped at the largest weight among the objects a junction
fragment touches instead of at 1. The production path of this build, through the unchanged consumers and the
EM, against the benchmark's harness arm `land_map_cap4` on the same conditions (`end-to-end-identity`):

| condition | this build's production path | the benchmark's `land_map_cap4` | the shipped reader |
|---|---|---|---|
| g05 ss0.99 capture-off | 9.46 / 0.90 | 9.49 / 0.90 | 6.42 / 0.89 |
| g05 ss0.99 captured | 7.14 / 0.36 | 6.91 / 0.36 | 4.93 / 0.35 |
| g25 ss0.99 capture-off | 8.81 / 0.98 | 8.67 / 0.98 | 7.27 / 0.98 |
| g25 ss0.99 captured | 7.21 / 0.77 | 7.30 / 0.77 | 6.84 / 0.77 |
| g50 ss0.99 capture-off | 9.30 / 1.22 | (not in the benchmark's run) | 9.28 / 1.16 |
| g50 ss0.99 captured | 10.26 / 1.28 | (not in the benchmark's run) | 10.33 / 1.20 |

Genes identical, transcripts within 0.03–0.23 points: the pair's own variation (the benchmark's control arm
moved up to 0.45 on one of these rows between runs), from weights that differ at 1e-10 through the EM's
near-tie sensitivity. The gate is passed: the production path IS the benchmark's recommended arm.

**The extrapolation to the real libraries, from their measured slot populations** (`inputs/real__*.npz`):
LBX0190 has 2,087,476 slots of which 98.9 % carry no read and 22,473 do (20,339 with fewer than ten);
MO_3021 90.7 % none, 194,341 with; VCaP 60.9 % none, 816,798 with (375k between ten and a hundred). At the
costs above (the zero-read search at 56 µs with its final windows) with 14 % of read-slots both-strand:
LBX0190 about 210 CPU-s, **26 s on 8 threads**; MO_3021 about 800 CPU-s, 100 s; VCaP about 5,400 CPU-s,
11 min. The plasma libraries, the product's case, are inside
the target; a deep library is governed by its shallow both-strand slots, which the derived tilt count is
meant for.

### Stage 4 (2026-10-10, evening): landed in the worktree, one red gate by ruling

Done in `~/proj/rigel-honest-reader` (uncommitted; the main tree untouched):
- `CONSTANTS.calibration.reader_search_stride` (documented, range-checked, `BROKEN` entry).
- Deleted: the located-mode census (`landscape.py`: `LandscapeMode`, `_census`, `split_basins`, `LocatedMode`,
  `located_enriched_mode`, the kernel record `centre`/`located`), `capture_efficiency.py` and its test, the
  `rate_log_floor` constant, `gdna_reference_density` / `gdna_reference_members` and their validators, the CLI
  summary keys, the report's enrichment panel (model, html, js) and its tests, the `None` return in
  `transcript_capture_eff_lengths`, the `np.minimum(eff_len, span)` clamp, the `[0, 1]` ceiling on the weights,
  the six mode tests, the publication test, the two reference-schema tests; the layer table's entry.
- The weights' unit is the MEDIAN read object (`calibrate.capture_weights`), not the largest: on every real
  library the largest of ~2.09 M noisy objects was an outlier and every published weight sat below 0.01
  (`prod_*_site10.json`: `w<1e-2` 2,077,402 of 2,087,500 on VCaP). The EM is invariant to the unit —
  `unit_invariance.py`, `combo_extreme` with every weight × 1e-3, × 0.37, × 1e-6: counts and TPM to 1e-15
  relative, gDNA and nascent totals identical; only `em_effective_length` scales.
- New gates: `tests/native/_honest_reader_reference.py` (the executable specification, the prototype's
  likelihood verbatim) + `tests/native/test_honest_reader.py` (the curve at every slot of a 48-slot fixture
  covering every branch to 1e-9, the modes, the zero-read search against the full grid, thread identity, the
  flat-evidence and pure-DNA limits); `tests/calibration/test_capture_reader.py` (no refit ⇒ 1 everywhere; the
  reader on the last sweep only; the published weights are the modes relative to the median; a forced set of
  modes is published exactly; the split onto the two axes); the junction's touched cap
  (`test_a_junction_is_captured_no_more_than_the_most_captured_object_it_touches`, fires at the old cap of 1);
  the ruler-arm test reads the plain `effective_length`.
- Docs (move rule): `DESIGN.md` §7.2 rewritten (the ruling, the receipts, the two open prices), `EQUATIONS.md`
  §11 (the likelihood, why the mode, the median unit, the touched cap), `ISSUES.md` (OPEN:
  `the-calibration-count-is-blind-to-multimappers` with the owner's 2026-10-10 ruling — buffer multimappers
  and assign them in the accumulator's second pass, after the reader — and
  `the-capture-weights-are-noisy-on-a-capture-off-library`; CLOSED: `the-capture-reference-is-read-at-a-grid-point`,
  `the-located-mode-capture-reader` with every refused reference's killing numbers; the detector entry's
  priority line), `SUCCESS.md`, `ROADMAP.md`, `TESTING.md`, `MANUAL.md` (the summary block, the capture
  paragraph); the docs gate passes. `ruff check` and `ruff format` clean.
- The instruments: `ruler_vs_truth.py`'s oracle arm reads the truth densities relative to the largest truth
  density (the reader's estimand with no estimator noise), `calibration_vs_oracle.py` drops the two kwargs.

Production path, pinned (scan threads 1, fractional), main → worktree (`prod_real.py`, `prod_*_{main,site10}.json`):

| library | wall s | gDNA fraction | mRNA | nRNA | gDNA (EM) |
|---|---|---|---|---|---|
| VCaP mix (truth 0.2518) | 178 → 454 | 0.2395 → 0.2401 | 13,336,410 → 13,334,253 | 755,354 → 746,792 | 4,355,162 → 4,365,881 |
| LBX0190 | 9 → 22 | 0.0850 → 0.0976 | 126,512 → 123,408 | 7,699 → 8,959 | 10,109 → 11,954 |
| MO_3021 | 21 → 66 | 0.1587 → 0.1631 | 632,561 → 624,526 | 69,895 → 74,247 | 79,444 → 83,128 |

The reader is 4.6 min of VCaP's 7.6 (the last refit's sweep: 12:11:12 → 12:16:03), 11 s of LBX0190's 22, 43 s
of MO_3021's 66; peak RSS unchanged (10.7 GB VCaP, 6.2 GB plasma against 4.0 under main). The deep-library
cost is the open performance item (the plan's flagged reductions, each a converged lattice change with its
gate). The plasma libraries' gDNA fraction rises 3–15 % relative with no truth to judge it; VCaP's moves
0.0006 toward its truth.

The suite under the stage-4 build (`worktree_suite12.log`, the goldens regenerated under the median unit):
2,805 collected, **1 failed / 2,804 passed** — the one failure is the paralog gate
(`the-calibration-count-is-blind-to-multimappers`), red by the owner's ruling until the accumulator's second
pass assigns multimappers. The goldens under the two units (`golden_site10/` against `tests/golden/`) agree on
every EM column to 1.2e-13 relative (counts and TPM to 4e-15); only `em_effective_length` and `gdna_eff_len_em`
differ, by the unit (up to 1.82× on these small scenarios). Against HEAD the goldens' movement is the one
recorded above (counts up to 28 % on `combo_extreme`, at most 0.3 % on 13 scenarios), the truth-scored
instruments read beside it (the golden truth table, `quant_accuracy.py` in the review's A.4).

Not run in the worktree: `preflight.py --full` (the instruments import the installed package, not a prototype
site); it belongs to the landing into the main tree, with `pip install -e`, the suite and the commit — the
owner's.
