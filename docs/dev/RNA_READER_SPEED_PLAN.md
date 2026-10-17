# The capture reader's deep-library speed — design (2026-10-11, sandbox, not authoritative)

**Status: LANDED on main at fa872c60 (owner, 2026-10-11): the per-node scope, no branch-and-bound; the reader on VCaP 4,643 → 1,106 CPU-s (the pipeline 477 → 253 s), the 60 s gate at the boundary, the 30 s aim a quadrature ruling after 0.8.0 (§9). This note is the campaign's receipt; the rulings and numbers live in DESIGN §7.2, EQUATIONS §11 and ISSUES.** The reader landed 2026-10-10 (`DESIGN.md` §7.2) at 11 s of cfRNA-A's 22 s
pipeline and 43 s of cfRNA-C's 66 s, and at 4.6 min on top of the VCaP mix's 3 min — the deep-library cost the
owner's performance stance refuses ("the entire tool runs in ~3 minutes on 8 cores"). This note derives where
the time goes, two reductions that are EXACT under the window rule the kernel already uses, a third that is
exact on one side, the gates each must pass, and the stages. Nothing here changes a number: every reduction is a
numeric no-op on the mode to the gate's tolerance, or it is not kept.

## 1. Where the time goes — the cost model, from the kernel

`native/honest_reader.h`, `Evaluator::log_L`: for every density node `x` of the landscape's grid (`nx = L.G`,
the full 260-node grid — no restriction is applied in the block solve), a both-strand slot with reads
integrates the share over the arcsine continuum: `K_t = 24` midpoint tilt nodes plus a `2⌈√(2T)⌉+1 = 17`-node
window around the observed share when the columns resolve it (`located`), plus the endpoints, then a SECOND
tilt pass adding windows around the per-`x` maxima of the first (columns reused). Each tilt node calls
`inner`, the `y`-integral over the RNA amount: the coarse lattice — `⌈(20 + 0.2)/0.2⌉ = 101` nodes over
`[y_hi − 20, y_hi]` — plus a 17-node window around the own-term mode `y*` and a 17-node window around the
peak found on those, every node once (`merge3`). Per `(x, tilt, y)` node `log_f` costs one `exp` (cached for
the coarse nodes, `r_coarse`), up to two `log`s (the Poisson terms), the composition interpolation (cached per
`x` across tilt nodes, `comp_coarse`) and the held RNA levels' interpolations (`level_terms`, two
`interp_held` per node when both are delivered).

    per both-strand read slot ≈ nx · N_tilt · N_y · c  =  260 · ~35 · ~135 · ~30 ns  ≈  37 ms

**Measured (stage 0, 2026-10-11, VCaP through the pinned production path, the kernel instrumented per slot;
`prototypes/2026-10-11_reader_speed/vcap_profile.npz`, read by `scratchpad/reader_profile.py`):** 2,087,476
slots, 4,643 CPU-seconds summed over the solve's threads (the 4.6 min of wall).

| class | slots | CPU s | share | mean / slot |
|---|---|---|---|---|
| pure DNA or not admitted (closed form) | 269,858 | 0.9 | 0.0 % | 3 µs |
| zero-read, admitted (the search) | 1,042,252 | 117 | 2.5 % | 112 µs |
| single-strand with reads | 671,648 | 590 | 12.7 % | 0.88 ms |
| both-strand with reads | 103,718 | 3,935 | **84.7 %** | 38 ms |
| … no level | 60,073 | 1,315 | 28.3 % | 22 ms |
| … one level | 41,161 | 2,396 | 51.6 % | **58 ms** |
| … both levels | 2,484 | 224 | 4.8 % | 90 ms |
| both-strand, n in [1, 10) | 51,546 | 1,387 | 29.9 % | 27 ms |
| both-strand, n in [10, 100) | 49,606 | 2,280 | 49.1 % | 46 ms |
| both-strand, n ≥ 100 | 2,566 | 269 | 5.8 % | 105 ms |

The model's shape holds (a no-level both-strand slot with `n < 24` has `26 × 101 × 260 = 683 k` node evaluations
at the measured ~32 ns each = 22 ms), with two corrections the estimate missed:

1. **A delivered RNA level nearly triples the node cost** (22 → 58 ms): `level_terms` evaluated the level's
   argument shift `log(origin) + log(Er) − log(f)` — three `log` calls — at EVERY node, and the composition
   row's interpolation was redone for every tilt node at the same density node (the loop ran tilt-outer, so the
   per-density cache missed every call). Half of VCaP's both-strand cost is this.
2. **Half the both-strand cost is slots with fewer than ten reads**, whose integrand is broad in every
   coordinate: the reductions that prune nodes far from a peak have little to prune there (§3).

**Where the time must come from.** The first pass of the share lattice (24 midpoints + endpoints + the located
window, at every density node) is ~75 % of a both-strand slot; under the landed specification it cannot be
pruned at all (§4), so the pooled second pass is the decision this note puts to the owner.

## 2. The window rule, which every reduction below reuses

The kernel already drops what cannot matter by ONE rule: a node whose log-integrand lies more than
`T = −log ε₆₄ ≈ 36.04` nats below the maximum contributes less than `ε₆₄` of the integral, below the
resolution of the float that holds it (`window_sd = √(2T)` is the window's half-width in standard deviations:
a Gaussian's log-density at `√(2T)` sd is `−T`). The reductions are the same statement applied to nodes the
kernel currently evaluates and then finds negligible. Each is exact to `T`: the dropped mass is at most
`N_nodes · e^{−T}` relative to a single retained node's, i.e. ~`10^{−13}` of the integral, two hundred times
below the gate (`prototype-identity`, the mode to `1e-6`).

## 3. What was built, and the two reductions that were dropped on the measurement

**Stage 1 — the restructure (bit-identical).** The density-node loop is outermost: the composition row on the
coarse lattice, the own term's mode `y*` and the window around it (its nodes, their `r` and composition values)
are set once per density node and read by every strand assignment; the level shifts' logs are hoisted out of the
node loop; each assignment's level terms on the coarse lattice (density-independent) are evaluated once per
slot. Every expression is the specification's, term by term in the same order. **Gate:** modes AND posterior
maxima bit-identical to the landed kernel on all four fixtures (3,518 / 3,518 / 3,518 / 70,176 slots), full
grid and production stride; the specification gates pass. **Speed:** 2.0–2.7× on the fixtures (concurrent
runs, indicative only).

**Reductions A and B, dropped.** Both were estimated before the measurement and are refuted by the integrand's
shape at the slots that carry the cost:

- A (prune the `y`-lattice by an envelope, exact to `T`): at every density node in the gDNA regime
  (`a2 ≳ r`) the Poisson pair saturates as `r → 0` and the integrand's low side is the measure `e^{y/2}` alone —
  20 nats of lattice is only 10 nats of fall, far inside `T = 36`. Nothing is dropped where the mode is. In the
  RNA regime the low tail falls as `n·log(r*/a2)`, which is inside `T` for the small-`n` slots that are half the
  cost. A prunes a few window nodes on deep slots and nothing on the rest.
- B (skip a tilt column by its bound, exact to `T`): the share posterior's width `sd_s = 1/(√(n+1)·|2q−1|)` is
  0.18 at `n = 50` and the `T`-window is ±8.5 `sd_s` — the whole of `[0, 1]` below `n ≈ 1,000`. B skips nothing
  on 99 % of the slots.

C (prune the density grid) survives in a stronger form: for the MODE the margin need not be `T`. A node whose
posterior upper bound lies below the best posterior already found cannot be the argmax, so the margin is only
the bound's rounding (`CONSTANTS.calibration.reader_prune_margin = 1e-6` nats). That is stage 2.

## 4. Stage 2 — the branch and bound over the density grid (exact to the margin)

**The bounds** (`bound_cheap`, `bound_tight` in the kernel; `tests/native/test_reader_envelopes.py` holds each
against the specification with a deliberate defect shown to fire). The integrand at a density node is
`Pois(u; a2 + p) · Pois(v; a2 + m) · C · H_pos · H_neg · r^{1/2}` with `p + m = r` the two strands' RNA amounts,
and any node's trapezoid sum over `y` is at most the lattice's span (26 nats) times the integrand's supremum;
the arcsine mixture and the three hypotheses are averages, so one bound covers them all:

- cheap (30 ns): `Σ_col log Pois(col; max(col, a2)) + ½(y_hi + 1) + max C + max H_pos + max H_neg + log 26 + D(x)`
  — each Poisson factor at most its maximum over rates at or above the gDNA share.
- tight (~1 µs): the own term `−dna − (p+m) + u log(a2+p) + v log(a2+m) + ½ log(p+m)` is JOINTLY CONCAVE in
  `(p, m)`, and the admitted shares form a convex cone, so its maximum is the interior stationary point —
  `r` solves `2r² + (4a2 − 1 − 2n) r − 2a2 = 0`, `p = 2ur/(2r−1) − a2`, `m = 2vr/(2r−1) − a2` — when that point
  lies in the cone, else the larger of the two rays' maxima (each a concave line search, `rna_mode_y`). Then
  `+ max C + max H + log 26 + D(x)`.

**The walk** (`mode_by_bound`): the cheap bound plus the prior orders the grid; walking it downwards, a node is
evaluated only while its cheap bound, then its tight bound, could still beat the best posterior found within
the margin (the sorted walk stops at the first node that cannot); the argmax's neighbours are evaluated for the
parabola. The pruning is sound for any slot; it is EXACT (the full grid's mode) whenever the argmax is unique by
more than the margin.

**Why it needs the second pass's scope changed.** The landed specification pools the second pass's windows over
the whole density grid: node `x`'s share lattice contains the windows around EVERY node's best share. So one
node's `log L` depends on every other node's first pass, and a node that will not be the mode must still be
evaluated to know which window it contributes — and in the gDNA regime its best share is a near-tie decided
by rounding, so even a `T`-exact shortcut there can change the union and move every other node by the lattice
difference. Under that specification only the second pass and the atoms (~25 % of a slot) can be pruned.

The variants, both built behind `Rules.second_pass` and held to the specification (which now takes the same
switch) by the curve and mode gates: **1 — each density node's own window**, at integer multiples of `sd_s`
around its own best share; **2 — no second pass** (the base lattice plus the located window). Under either,
`log L(x)` is a function of `x` alone and the branch and bound is exact to the margin. The number this changes
is measured on the fixtures (§8) and is the owner's ruling.

## 5. The gates, as built

| gate | holds |
|---|---|
| `bound-is-a-bound` | `tests/native/test_reader_envelopes.py`: the cheap bound ≥ the specification's `log L` at every node of 50+ random read slots, and FAILS with the Poisson rate held at the gDNA share; the concave program's value is the brute-force maximum on a fine (amount, share) grid (never below, never 0.05 nats above) and FAILS with the cone's rays left out; the tight bound ≥ `log L` at every node, within 6 nats of the maximum on some slot; the native bound equals the Python one to 1e-9 |
| `prototype-identity` | `test_honest_reader.py`: the curve and the mode against the specification at every slot, under each second-pass scope (the fixture covers every branch) |
| `prune-identity` | `test_honest_reader.py`: with the branch and bound, the mode to 1e-7 and the posterior maximum to 1e-9 against the full grid at every read slot, under scopes 1 and 2 |
| `thread-identity` | bit-identical at 1 and 4 threads, with the pruning on |
| `fixture-identity` | stage 1 bit-identical to the landed kernel on the four fixtures (modes and maxima); stage 2's scope-0 path likewise |
| `end-to-end-identity` | cfRNA-A, cfRNA-C, VCaP through the pinned production path: the transcript table within the pair's own variation, the gDNA fraction to four decimals |
| `speed` | the reader stage on VCaP: at most 60 s, aim 30 s (from 276 s); timings from back-to-back pairs on an idle machine |

## 6. What is NOT proposed

- No change to the lattices' definitions, the window rule's `T`, `K_t`, the coarse step or the readout: every
  evaluated node is evaluated exactly as today.
- No coarse-to-fine search in `x` for read slots: refused (plateau ties, `RNA_HONEST_READER_NATIVE_PLAN.md`
  §7). The branch and bound is a bound, not a search: a pruned node is PROVED below the mode.
- No change to the derived tilt count or `skip_far_tilt` flags (both off, to be deleted at landing).
- No approximation of the levels, the composition rows or the logs (a vectorised `log` would be a further ~1.5×
  on every class and is noted as an option, not proposed).

## 7. Stages

0. **Measured** (2026-10-11): §1.
1. **Built, gated** (bit-identical): the restructure. `prototypes/2026-10-11_reader_speed/site1`.
2. **Built, gated** (exact to the margin): the bounds and the branch and bound, behind `reader_prune`; the second
   pass's scope behind `reader_second_pass` (0 pooled = the landed specification, 1 per node, 2 none), the
   specification taking the same switch. `site2`. The constants are PROTOTYPE knobs: the ruling replaces them.
3. **Measuring**: on the fixtures, the scope variants' change against the pooled specification (every slot) and
   the pruning's identity and speed; on VCaP, cfRNA-A and cfRNA-C through the pinned production path, the reader
   stage's time for (scope 0, no pruning) — the restructure alone — against (scope 1, pruning) and (scope 2,
   pruning), and the end-to-end identity.
4. **Land** (after the ruling): the chosen scope hard-wired, the knobs and the flags `derived_tilt` and
   `skip_far_tilt` deleted; `DESIGN.md` §7.2 gains the branch and bound and the scope as a paragraph,
   `EQUATIONS.md` §11 the two bounds and the concave program; the suite, `preflight.py --full`. The commit is
   the owner's.

## 8. The measurements (filled as they come)

Stage 1 on the fixtures (8 threads, concurrent with other work — indicative): 3.8 → 1.5 s, 1.7 → 0.7 s,
2.2 → 0.8 s, 71.6 → 31.6 s (the ladder's 70,176 slots); bit-identical throughout.

**The scopes against the pooled specification, every slot of the four fixtures, production stride
(`scratchpad/compare_scopes.py`; only both-strand read slots can move, and no other slot did):**

| fixture (both-strand read slots) | scope 1: max Δmode, slots > 1e-6 / > 1e-3 | scope 2: max Δmode, > 1e-6 / > 1e-3 / > 1e-2 |
|---|---|---|
| benign g05 OFF (188) | 5.3e-11, 0 / 0 | 2.3e-2, 78 / 16 / 2 |
| benign g05 ON (88) | 2.1e-5, 1 / 0 | 8.9e-3, 42 / 8 / 0 |
| benign g50 ON (123) | 8.5e-7, 0 / 0 | 6.1e-2, 52 / 5 / 1 |
| ladder g50 ON (5,746) | 3.3e-5, 320 / 0 | 1.8e-1, 2,537 / 46 / 9 |

Scope 1 (each node's own window) is a change below 3.3e-5 nats at every slot — the lattice difference of a
window that a far-away node contributed; scope 2 (no second pass) moves modes by up to 0.18 nats and is
refused. **The pruning is bit-identical to the unpruned mode AND posterior maximum at every slot of every
fixture, under both scopes** (the margin never decided anything).

**The real libraries, pinned production path, back to back on the idle machine (`vcap_ab.sh`; the reader's
CPU-seconds summed over the solve's ~17 threads; A = the restructure alone under the pooled scope, B = scope 1
with the branch and bound):**

| library | pipeline wall A → B | reader CPU-s: landed → A → B | both-strand read slots (ms each, B) | single-strand | zero-read |
|---|---|---|---|---|---|
| VCaP (2.09 M slots) | 311.8 → 294.7 s (landed 477) | 4,643 → 2,047 → 1,777 | 1,490 → 1,253 (103,718 × 12.1 ms) | 475 → 444 | 81 → 79 |
| cfRNA-A | 20.3 → 20.1 s (landed 22) | — → 87 → 84 | 19 → 16 (2,753 × 5.9 ms) | 4 → 7 | 63 → 60 |
| cfRNA-C | 57.6 → 57.8 s (landed 66) | — → 275 → 287 | 142 → 135 (21,147 × 6.4 ms) | 38 → 61 | 95 → 90 |

End to end, B against A: VCaP's gDNA fraction 0.252753 in both, the pools within 0.1 fragments, the weights'
quantiles identical to four figures, the reader's modes within 9.6e-6 nats (17 slots above 1e-6; scope 1's
lattice difference); the transcript table moves 1,561 fragments in |Δ| on 31 rows, all within-gene
redistributions among near-identical isoforms (the largest gene moves 759 among three isoforms with a net of
−1.7; three of eleven genes have a net above one fragment) — the known identical-isoform tie collapse under a
1e-5 prior difference, not an accuracy change. cfRNA-A and cfRNA-C: modes within 1e-10, tables within 2e-3.

**Reading.** The restructure is the whole gain so far (2.27× on VCaP's reader; the plasma readers 11 → ~5 s
and 43 → ~16 s of wall). The branch and bound removes only ~16 % of VCaP's both-strand cost and MAKES the
single-strand slots slower on the plasma libraries (the tight tier's two line searches cost more than a shallow
slot's own evaluation) — the tight tier is now applied to both-strand slots only. Why it prunes so little is
the next measurement: each slot's evaluated-node count and the bound's gap at the argmax (prototype
counters), against the landscape's own breadth — whether the bound's looseness (`log 26` for the span, the
composition's and the levels' maxima) or the prior's breadth on the flat side sets the survivors.

**Why the branch and bound prunes so little — measured (prototype counters, VCaP, scope 1 with pruning;
`scratchpad/prune_stats.py`).** The landscape has 260 nodes at 0.080 nats over 20.8 nats, a `logP` range of
13.7 nats and four local maxima, and it is BROAD: 17 % of its nodes lie within 0.5 nats of its maximum, 37 %
within 1 nat, 75 % within 2 nats. The tight bound's gap at the argmax is 3.2 nats at the median (2.7 on
no-level slots, 3.9 with one level, 5.9 with both; 2.5 on single-strand slots). So a both-strand slot evaluates
176 of 260 nodes at the median (p90 190) — the survivors are the prior's top-3-nats region, which is three
quarters of the grid. Even a bound tight to 1 nat (a Riemann upper sum over the lattice's cells with each cell's
supremum, ~1 ms per slot) would leave 37 % of the grid, i.e. at most ~1.4× net on the both-strand class. The
pruning is capped by the landscape's breadth, which is a property of the data, not of the bound.

**The per-node cost — measured.** A lattice node costs ~22 ns in the restructured kernel (a no-level both-strand
slot: 63 µs per density node over ~28 columns × 101 nodes). Of that, libm's two logs and one exp are at most
~8 ns (micro-benchmark on this machine, 101-element chunks: `log` 3.2 ns, `exp` 2.0 ns per element; Accelerate's
`vvlog` 1.5 ns, a derived-coefficient series `log` 1.5 ns, `exp` 1.5 ns). A series `log`/`exp` written as
vectorisable loops (`native/vmath.h`, within 2 ulp of libm, its `exp` within a derived bound; the gate
`tests/native/test_vmath.py`) changes every fixture mode by at most 1.3e-11 and buys 1.05–1.11× on the
fixtures — the real-library number is below. The other ~14 ns per node are loads, the trapezoid weights'
logs on windowed slots (one `log` per node per column, avoidable by weighting the exponentials instead of
adding log-weights), the level and composition interpolations on window nodes, and the per-column passes.

**The derived tilt count (`derived_tilt`, level-free both-strand slots: `min(24, ⌈(π/2)/sd_s⌉)` base nodes),
measured for the record:** modes move by up to 5.0e-3 (benign g05 ON, 4 slots above 1e-3) and 2.5e-3 on the
ladder (1,108 of 5,746 both-strand slots above 1e-6, 7 above 1e-3); no speed on the fixtures (few such slots
there). It is a quadrature change two orders above scope 1's and is not proposed; it stays the one lever with
large headroom on the small-read slots that are half of VCaP's cost, and a ruling if the aim is to be met.

**The series `log`/`exp` on VCaP (site5 against site2, scope 1 with pruning, back to back, idle):** the reader
1,777 → 1,440 CPU-s (1.23×): both-strand 1,253 → 1,091, single-strand 444 → 274 (1.62×: its nodes were
transcendental-bound), zero-read 79 → 74; the pipeline 294.7 → 275.8 s; modes within 7.9e-11, the gDNA fraction
identical, the transcript table's changes again within-gene tie flips (24 genes, the largest 1,744 fragments
redistributed with a net of −0.16). Kept.

**The micro-stage (site6; exact to rounding, the interpolations bit-identical):** the trapezoid weights applied
to the exponentials instead of added as logs (one `log` per node per column gone on windowed slots); the
composition row's and the windows' held interpolations searched from the previous node's interval instead of
by binary search at every node (the arguments are an ordered run, so the interval moves by at most one knot:
the same interval, the same value bit for bit); the level rows likewise. **Measured (site6 against site5, VCaP, back to back, idle):** the reader 1,440 → 1,106
CPU-s (1.30×): both-strand 1,091 → 857 (8.3 ms per slot), single-strand 274 → 189, zero-read 74 → 59; the
pipeline 275.8 → 252.6 s; modes within 4.0e-12; the gDNA fraction identical; the transcript changes within-gene
tie flips (25 genes, the largest 1,267 fragments redistributed with a net of +0.01, no gene's net above 2.7).
Fixtures against the libm build: modes within 1.3e-11, 1.55–1.80×.

**The state of the reader on VCaP, the whole campaign** (CPU-seconds summed over the solve's threads; the
reader's wall is 52 s by the pipeline's difference from the landed run and 66 s by the CPU ratio to the landed
276 s — the two bracket it):

| | landed | stage 1 (bit-identical) | + scope 1 and the branch and bound | + the series log/exp | + the micro-stage |
|---|---|---|---|---|---|
| reader CPU-s | 4,643 | 2,047 | 1,777 | 1,440 | **1,106** |
| both-strand read (103,718 slots) | 3,935 (38 ms) | 1,490 | 1,253 | 1,091 | 857 (8.3 ms) |
| single-strand read (671,648) | 590 | 475 | 444 | 274 | 189 |
| zero-read (1,042,252) | 117 | 81 | 79 | 74 | 59 |
| pipeline wall (pinned scan) | 477 s | 312 s | 295 s | 276 s | **253 s** |
| modes against the landed | — | identical | ≤ 3.3e-5 (scope 1) | ≤ 8e-11 further | ≤ 4e-12 further |

The fixtures against the landed kernel: 7–14× (the ladder's 70,176 slots 59.9 → 8.7 s). cfRNA-A and cfRNA-C were
measured through stage 2 (87 → 84 and 275 → 287 CPU-s, their readers ~5 and ~17 s of wall; the single-strand
slowdown there is what restricted the tight tier to both-strand slots) and not re-timed since.

**Where the remaining 1,106 CPU-s sit.** 77 % in the both-strand read slots (8.3 ms each: ~150 of 260 density
nodes evaluated, ~28–60 share columns at each, 101–135 lattice nodes per column, ~9 ns per lattice node now),
17 % in the single-strand read slots, 5 % in the zero-read search. The exact levers are spent: the pruning is
capped by the landscape's breadth, the transcendentals are at the series' cost, the interpolations and weights
are O(1) per node. What is left is the lattice itself — the share lattice's K_t = 24 base nodes at every density
node of a slot whose share posterior is as wide as the whole interval (`derived_tilt`, measured at up to 5e-3 on
the modes), and the density grid's 260 nodes — each a quadrature ruling, not a no-op.

Speed on the fixtures (8 threads, the machine otherwise idle), scope 1 with pruning against the landed kernel:
benign 4.4× / 4.4× / 8.3×, the ladder 59.9 → 13.5 s (4.4×); the pruning alone is 1.3–2.2× over the unpruned
scope-1 run (the fixtures are mostly zero-read and single-strand slots, which it does not touch). The real
number is the VCaP stage below.

## 8b. What lands, what is prototype

To land (after the ruling on the scope): the restructured evaluator, `native/vmath.h` with its gate, the scope
hard-wired, the branch and bound with `reader_prune_margin` (the only new constant; `reader_second_pass` and
`reader_prune` deleted as knobs), the flags `derived_tilt` and `skip_far_tilt` deleted, the specification
`_honest_reader_reference.py` with the scope, the gates `test_reader_envelopes.py`, `test_vmath.py` and the
parametrised `test_honest_reader.py`; DESIGN §7.2 a paragraph, EQUATIONS §11 the bounds and the concave program.
PROTOTYPE, never lands (marked in the source): the per-slot profiling in `solve_kernel.cpp` and `calibrate.py`
(`RIGEL_READER_PROFILE`, the counters `n_eval`/`gap`), the env knobs `RIGEL_READER_SECOND_PASS` and
`RIGEL_READER_PRUNE`, the prototype dir's scripts.

## 9. The owner's decisions

1. **The second pass's scope.** The pooled scope (the landed specification) makes every density node's lattice
   depend on every other node's first pass, which blocks the branch and bound on the first pass (75 % of a
   both-strand slot) and is also a fragility: a rounding-decided tie at one irrelevant node changes every node's
   lattice. A per-node scope (1) or none (2) is the natural statement of "refine where this node's peak is",
   makes `log L(x)` a function of `x` alone, and changes the mode by the lattice difference — §8 has the measured
   size. Taking it is a number change and the owner's ruling.
2. **The target.** The reader on VCaP is now 52–66 s of wall (from 276): the 60 s gate at its boundary, the
   30 s aim not reachable by exact means — every remaining lever is a quadrature ruling (the share lattice's
   node count against the share posterior's width, the density grid's resolution), each a number change of
   order 1e-3 to be measured per stratum before it is taken. Whether 50–65 s on a 20-million-fragment
   library is acceptable for 0.8.0, or one of those rulings is wanted, is the owner's call.
3. **Scope 1 without the pruning** is also a consistent landing (the branch and bound is worth 13 % on VCaP and
   nothing on the plasma libraries; it adds `mode_by_bound`, the two bounds and a constant): the owner may
   prefer the simpler tree. The restructure, the series math and the micro-stage are bit-identical or
   rounding-level and carry no ruling.
