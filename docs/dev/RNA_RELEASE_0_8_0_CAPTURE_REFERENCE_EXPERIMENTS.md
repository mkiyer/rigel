# The capture-reference review's sign-off experiments — results (2026-10-13)

**What this is.** The independent review
([`RNA_RELEASE_0_8_0_CAPTURE_REFERENCE_INDEPENDENT_REVIEW.md`](RNA_RELEASE_0_8_0_CAPTURE_REFERENCE_INDEPENDENT_REVIEW.md))
kept the integrated two-sided reader as the 0.8.0 candidate and refused the one-sided patch. It held sign-off on two
experiments:

- the **oracle bulk**, with a normalization control (its experiment 1);
- **stress causality** (its experiment 2).

Both are done; results are below. One finding from the normalization control reaches beyond the reader (§1). Nothing
is committed.

**Protocol.**

- Scan pinned, fractional assignment, `OMP_NUM_THREADS=1`, EM seed 0.
- One condition per run where an arm needs a per-condition input.
- Every table states the EM settings it used: **shipped** (convergence tolerance 1e-6, at most 1,000 iterations) or
  **converged** (tolerance 1e-12, at most 100,000 iterations).

**Arms.** Each arm is a frozen prototype site built from the release build: `site_obulk`, with switches
`RIGEL_BULK_SCALE` and `RIGEL_ORACLE_BULK`.

| Arm | What changes |
|---|---|
| **Rescale arm** | Multiplies every efficiency, the saturation level and the zero-support fallback by 37 |
| **Oracle-location arm** | Moves only the null spike's location to the oracle density. Its likelihood there is in closed form, and saturation and the fallback are rebased on it. Basin masks, prior masses, the alternative, absolute cap locations and calibration are unchanged |

**The oracle** is the simulator's *expected* uncaptured gDNA density, not a realized count. The simulator places a
fixed number of gDNA fragments in proportion to each placement's capture weight, so:

- an unprobed placement expects `N_gDNA · off_target_weight / Z`;
- `Z` is the genome-wide capture-weighted gDNA yield (`CaptureSampler.partition_array` over every reference and width,
  under the panel's gDNA length law).

It agrees with the realized Σ true ÷ Σ support within 0.95–1.07 wherever that is well measured (`expected_bulk.py`,
`review/expected_bulk.json`).

---

## 1. The normalization control, and what it found

Under the **shipped** EM settings, the rescale arm moved results by up to 0.84 points per condition. Every EM input scaled exactly:
all 918 transcript lengths and all 218 locus gDNA lengths went ×37, and the gDNA counts were identical. So nothing in
the pipeline consumes an absolute scale.

The movement is the **EM exhausting its iteration budget** in nearly unidentifiable loci, chiefly antisense-overlap
and capture-cluster blocks:

- a rescale of **1 + 1e-12** moves 3,172 of 581,199 fragments, about as much as ×37 (3,475);
- the largest single move is 928 fragments in one transcript;
- a rerun at scale 1 is identical to the bit;
- under the **converged** settings, the same perturbations move under 1e-5 fragments, and the rescale arm's worst relative
  difference in any result field is **3.8e-10**.

**The reviewer's follow-up locates it more precisely.** The reviewer replayed one captured locus (`B4_capcluster`, 14
transcripts, 27,064 units) with one EM thread (`.cache/capture_reference_followup/summary.json`):

- every run at the shipped budget exhausts its 333 accelerated cycles;
- tightening the tolerance without raising the budget is bit-identical;
- only the budget matters: with 100,000 iterations, a 1e-12 rescale moves the locus's summed component counts by 3.2
  at tolerance 1e-6 and by 1.3e-6 at 1e-12.

In exact arithmetic a common length factor cancels in every E-step. A finite-iteration answer should therefore also be
scale-invariant. That it is not points at the accelerated solver's numerics (SQUAREM extrapolates from differences of
successive updates), which the plain-VBEM control below is designed to test.

My runs above left the EM's own thread count at its automatic default. `OMP_NUM_THREADS=1` does not set the native
pool, so attribution runs must pin `em.n_threads=1`.

**Consequences:**

1. **Observed numerical sensitivity, not a noise floor.** Under the shipped settings, a rounding-level change to the
   EM's inputs moved a condition's transcript error by up to 0.84 points. That is an observation, not a statistical
   floor or a bound: it cannot justify dismissing every smaller difference or trusting every larger one. Attribution
   should be read beside sufficiently converged runs.
2. **Early stopping is a strong implicit regularizer.** Converging fully raises transcript error substantially, for
   both readers:

   | Condition | Release, shipped | Release, converged | Foundation, shipped | Foundation, converged |
   |---|---:|---:|---:|---:|
   | test g05 κ0.99 ON | 5.17 | 10.83 | 4.95 | 10.78 |
   | test g50 κ0.99 ON | 10.13 | 14.16 | 10.33 | 14.06 |
   | lowg g01 κ0.99 ON | 7.72 | 12.30 | 6.57 | 11.54 |
   | lowg g001 κ0.99 OFF | 7.56 | 11.04 | 6.09 | 10.34 |
   | d10 g25 κ0.99 ON | 14.37 | 18.19 | 14.24 | 18.27 |

   The shipped accuracy depends on where the EM is stopped, not only on its objective. Full convergence is not the
   only response: if acceleration causes the instability, a numerical repair may keep the useful early stop while
   restoring scale invariance. This is outside the reader, but it bears on every benchmark number in the project.

## 2. Experiment 1: the oracle bulk (converged EM)

| Condition | Uncaptured gDNA fragments | Census bulk ÷ expected | Foundation | Release | Rescale ×37 | Oracle null | Oracle − release |
|---|---:|---:|---:|---:|---:|---:|---:|
| d10 g001 ON | 4 | 2.44 | 15.68 † | 19.36 | 19.36 | 19.36 | 0.00 |
| lowg g001 ON | 27 | 0.26 | 22.50 | 23.82 | 23.82 | 25.54 | +1.73 |
| d10 g01 ON | 39 | 0.26 | 31.50 | 32.72 | 32.72 | 35.46 | +2.74 |
| d10 g001 OFF | 116 | 0.06 | 17.29 | 17.26 | 17.26 | 17.28 | +0.01 |
| d10 g05 ON | 140 | 0.05 | 17.29 | 17.62 | 17.62 | 18.08 | +0.46 |
| lowg g01 ON | 267 | 0.98 | 11.54 | 12.30 | 12.30 | 12.29 | 0.00 |
| d10 g25 ON | 601 | 1.09 | 18.27 | 18.19 | 18.19 | 18.14 | −0.06 |
| d10 g50 ON | 1,255 | 1.09 | 21.98 | 22.06 | 22.06 | 22.06 | 0.00 |
| test g05 ON | 1,311 | 1.12 | 10.78 | 10.83 | 10.83 | 10.82 | −0.01 |
| lowg g001 OFF | 1,134 | 1.05 | 10.34 | 11.04 | 11.04 | 11.02 | −0.02 |
| test g50 ON | 12,729 | 1.02 | 14.06 | 14.16 | 14.16 | 14.16 | 0.00 |

All conditions are κ 0.99. † The foundation finds no mode in this library (117 gDNA fragments in all), so it applies no
capture correction.

**Reading, against the review's own criterion.** "Recovery prices the spike-location mechanism. Non-recovery refutes
its sufficiency."

- **The null's location does not explain the sparse losses.** Placing it at the true uncaptured density makes the
  sparse libraries worse, not better.
- **The review's next step is therefore not indicated.** That step was the likelihood-fitted null location,
  conditional on this experiment.
- **The bulk estimate does collapse where uncaptured gDNA is sparse** (§5.2a of the packet). Measured end to end,
  that collapse does not cost transcript accuracy in these libraries; correcting it costs some.

**The sparse losses that remain under converged EM** (+1.2 to +3.7 points against the foundation) are largest where
the foundation applies no correction at all (d10 g001). That points to per-object reading noise when the evidence is
scant: the foundation's clip at one level, or no correction, regularizes it. The review's experiment 4, a
length-weighted weak-evidence analysis, is the test that would localize it.

## 3. Experiment 2: stress causality (converged EM)

**Design.**

- The `combo_extreme` design is repeated *K* times along the genome at constant depth: two genes, no capture, heavy gDNA
  and nascent RNA, κ 0.65.
- Each design runs with several read-simulation seeds.
- The third arm feeds the release reader the **true** gDNA deposits: every region's contained and every boundary's
  crossing gDNA fragments, from the oracle read names, plus a landscape refit on them. The rest of the pipeline is
  unchanged.
- Script: `stress_causality.py`.

| *K* | Seed | Foundation | Release | Release + true deposits |
|---:|---:|---:|---:|---:|
| 1 | 42 (the golden) | 12.44 | 36.77 | 13.41 |
| 1 | 1 – 7 (seven seeds) | 10.1 – 34.3 | equal, or within ±0.5 | 10.1 – 34.5 |
| 5 | 42 | 19.06 | 28.28 | 19.06 |
| 5 | 3 | 17.11 | 26.71 | 17.50 |
| 5 | 1, 2 | 17.79, 34.94 | 17.40, 32.27 | 17.92, 17.54 |
| 25 | 42, 1 | 22.38, 16.08 | 22.11, 16.46 | 16.50, 15.42 |

| *K* | Seeds | Mean: foundation | Mean: release | Mean: release + true deposits |
|---:|---:|---:|---:|---:|
| 1 | 8 | 16.90 | 19.97 | 16.21 |
| 5 | 4 | 22.23 | 26.17 | 18.01 |
| 25 | 2 | 19.23 | 19.29 | 15.96 |

**Reading.**

- **The origin is upstream.** Given the true deposits, the reader matches or beats the foundation at every size. The
  failure is calibration booking nascent RNA as gDNA at weak strandedness, which a small population then elevates to
  "the bulk".
- **The deposit oracle largely removes the regression, and the largest tested examples are near parity.** The large
  regression occurs in 1 of 8 seeds at one unit, 2 of 4 at five, and 0 of 2 at twenty-five. Those counts do not
  establish a falling failure probability, and two clean seeds do not establish safety.

## 4. Where this leaves the candidate

The review's two sign-off experiments are answered.

- **Experiment 1:** the bulk's placement is not what costs accuracy in sparse libraries, so no replacement baseline
  estimator is indicated for 0.8.0.
- **Experiment 2:** the stress regression has a located, upstream origin with falling prevalence. Per the review, it
  now needs either a validated local correction, which none of the measured candidates is, or an **explicit
  acceptance of this measured limitation**.

The early-stopping finding (§1) is separate from the reader and needs its own decision.

---

## 5. The reviewer's follow-up: the solver comparison and the convergence control (2026-10-14)

The follow-up kept the integrated reader, accepted the stress limitation, and asked for three things before release
freeze:

1. a focused solver comparison on the captured locus: ordinary against accelerated VBEM, from the same start, prior and
   update equations, with the common-scale controls, budget exhaustion recorded apart from convergence, and a
   plain-update residual;
2. the reader's sparse-probe benefits checked under the convergence control;
3. two benchmark views kept.

**The build.** Every run in this section pins `em.n_threads=1`. A prototype native build from the release head
(`site_solverdiag`) adds an ordinary-update path, switched on by `RIGEL_EM_PLAIN=1`. That path uses the same start,
grouped prior, update equations, floor and stopping rule, with no extrapolation, and its budget counts updates. The
build also reports two things per locus:

- whether the stopping rule fired, or the budget ran out;
- the **plain-update residual**: the normalized L1 change one more ordinary update would make at the answer.

With the switch off, the build is bit-identical to the release.

### 5.1 The captured locus

`B4_capcluster`: ten near-identical isoforms, each with about 1,850 true fragments, replayed from the reviewer's
fixture (`solver_replay.py`).

| Solver | Budget | Tolerance | Steps | Stopping rule fired | Plain residual | Moved by 1 + 1e-12 | Moved by ×37 |
|---|---:|---:|---:|---|---:|---:|---:|
| Accelerated (shipped) | 1,000 | 1e-6 | 333 cycles | no | 8.2e-5 | 3,239 | 3,657 |
| Accelerated | 1,000 | 1e-12 | 333 cycles | no | 8.2e-5 | 3,239 | 3,657 |
| Accelerated | 100,000 | 1e-12 | 7,019–17,018 cycles | yes | 2.9e-13 | 1.3e-6 | 8.3e-7 |
| Ordinary | 1,000 | 1e-6 | 1,000 | no | 1.7e-4 | 5.9e-10 | 1.6e-8 |
| Ordinary | 1,000 | 1e-12 | 1,000 | no | 1.7e-4 | 5.9e-10 | 1.6e-8 |
| Ordinary | 100,000 | 1e-6 | 85,540 | yes | 1.0e-6 | 2.9e-10 | 1.6e-8 |
| Ordinary | 100,000 | 1e-12 | 100,000 | no | 2.9e-11 | 7.1e-11 | 2.0e-8 |

"Moved" is the L1 change, in fragments, of the locus's component counts (gDNA included) from the scale-1 answer of
the same solver, budget and tolerance.

- **The acceleration is the source.** Ordinary updates move under 2e-8 fragments under either rescale, whether or not
  they have converged. At the shipped budget, the accelerated solver moves 3,239. Given enough budget, both solvers
  reach the same answer.
- **Why.** The accelerated step length is −(r·v)/(v·v). Its only upper bound is that every live component stays above
  the floor. Along an almost flat direction, v is a near-zero second difference, so the step is a ratio of
  rounding-sized numbers. The halving search for a feasible step is also a discrete decision, which a rounding change
  can flip. The answer at the budget is a snapshot of that path.

**The answer along the path, against truth** (`locus_trajectory.py`). Σ|estimate − truth| over the ten isoforms and the
locus's single-exon transcript:

| Answer | Σ\|estimate − truth\| | gDNA | The ten isoforms (truth 1,830–1,893 each) |
|---|---:|---:|---|
| Ordinary, 1 update | 5,475 | 4,943 | 1,670–1,678 |
| Ordinary, 10 | 1,579 | 5,415 | 1,820–1,860 |
| Ordinary, 30 | **867** | 5,553 | 1,818–1,921 |
| Ordinary, 100 | 1,041 | 5,613 | 1,756–2,063 |
| Ordinary, 300 | 2,365 | 5,627 | 1,593–2,423 |
| Ordinary, 1,000 | 6,208 | 5,627 | 1,176–3,362 |
| Ordinary, 3,000 | 12,060 | 5,627 | 453–4,652 |
| Ordinary, 10,000 | 19,046 | 5,628 | five at zero |
| Ordinary, 100,000 | 26,198 | 5,628 | seven at zero |
| Accelerated, shipped (333 cycles) | 15,089 | 5,627 | 74–5,199 |
| Accelerated, converged | 26,198 | 5,628 | seven at zero |
| MAP mode, accelerated, 100,000 (not converged) | 19,081 | 5,627 | four at zero |

- **There are two time scales.** The part the data determine settles fast: the gDNA count is fixed within about 300
  updates. The split among the ten isoforms then keeps drifting. At convergence, three isoforms carry the cluster and
  seven carry nothing.
- **Earlier is better here.** The best answer comes after about 30 updates. The shipped accelerated answer sits about
  where ordinary updates would be after several thousand.
- **The drift belongs to the estimator, not the accelerator.** MAP-mode EM drifts the same way. Along directions the
  data barely determine, the converged estimate is very noisy, and stopping early shrinks it toward the even split the
  EM starts from.

### 5.2 The test panel

Thirty conditions (`em_screen.py`, `em_screen_report.py`). "Moved" is the common-scale control: every capture
efficiency and the saturation level multiplied by 1 + 1e-12.

| Stratum | Conditions | Accelerated (shipped): tx % (gene %) | Ordinary: tx % (gene %) | Δ tx | Accelerated moved (fragments) | Ordinary moved | Loci at budget (accelerated / ordinary) | EM time, ordinary ÷ accelerated |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| κ 0.99 OFF | 4 | 7.62 (1.08) | 6.06 (1.08) | −1.56 | 2,425 | 4.2e-8 | 27 / 69 of 876 | 1.93× |
| κ 0.99 OFF g00 | 1 | 7.49 (0.84) | 5.22 (0.84) | −2.27 | 1,201 | 1.8e-8 | 10 / 23 of 200 | 1.97× |
| κ 0.99 ON | 4 | 7.43 (0.96) | 5.25 (0.98) | −2.17 | 4,526 | 3.0e-8 | 18 / 60 of 863 | 1.94× |
| κ 0.99 ON g00 | 1 | 7.35 (0.67) | 3.87 (0.67) | −3.48 | 3,429 | 1.0e-8 | 5 / 10 of 189 | 2.22× |
| κ 0.70 OFF | 4 | 8.14 (1.25) | 6.25 (1.25) | −1.89 | 3,573 | 2.6e-8 | 26 / 78 of 876 | 2.10× |
| κ 0.70 OFF g00 | 1 | 6.58 (0.90) | 4.59 (0.90) | −1.98 | 1,749 | 2.2e-8 | 8 / 21 of 200 | 2.18× |
| κ 0.70 ON | 4 | 8.01 (1.74) | 6.10 (1.75) | −1.91 | 10,862 | 4.9e-8 | 24 / 66 of 866 | 2.11× |
| κ 0.70 ON g00 | 1 | 8.31 (2.16) | 4.88 (2.16) | −3.43 | 1,689 | 2.6e-8 | 7 / 13 of 189 | 1.91× |
| κ 0.50 OFF | 4 | 8.82 (1.34) | 6.64 (1.35) | −2.18 | 2,017 | 4.3e-8 | 28 / 83 of 876 | 2.17× |
| κ 0.50 OFF g00 | 1 | 5.60 (1.01) | 4.78 (1.01) | −0.82 | 520 | 2.2e-8 | 11 / 24 of 200 | 2.23× |
| κ 0.50 ON (deferred) | 4 | 12.42 (3.27) | 9.41 (3.27) | −3.01 | 7,844 | 5.9e-8 | 28 / 76 of 866 | 2.12× |
| κ 0.50 ON g00 (deferred) | 1 | 9.23 (2.44) | 5.56 (2.44) | −3.67 | 1,025 | 1.4e-8 | 7 / 13 of 190 | 2.03× |
| **All** | 30 | **8.34** (1.53) | **6.06** (1.53) | −2.28 | 40,860 | 3.6e-7 | 199 / 536 of 6,391 | 2.07× |

Ordinary updates win on 29 of 30 conditions and tie on the other. Gene error is unchanged, so the whole effect is
isoform allocation within genes.

**The start does not matter.** From the uniform start, the accelerated solver scores 8.32 % and ordinary updates 6.06 %
(8.34 % and 6.06 % from the shipped coverage start). After one update, either start sits at the proportional split.

**The budget is the regularizer.** Ordinary updates, same thirty conditions:

| Ordinary updates, budget | Transcript % | Gene % | Loci at budget (of 6,391) | EM seconds |
|---:|---:|---:|---:|---:|
| 300 | 4.90 | 1.56 | 1,255 | 283 |
| 1,000 | 6.06 | 1.53 | 536 | 429 |
| 3,000 | 8.12 | 1.53 | 215 | 1,021 |
| 10,000 (5 conditions only) | 9.62 (accelerated on the same five: 7.46) | 1.48 | 26 of 979 | 501 |
| Accelerated, shipped | 8.34 | 1.53 | 199 | 207 |

The shipped accelerated solver travels about as far as 3,000 ordinary updates. At that matched distance, the two
solvers score nearly the same (8.12 % vs 8.34 %), and ordinary updates take five times the EM time.

### 5.3 The deciding panels

Accelerated against ordinary updates, both at the shipped budget, one EM thread (`em_screen.py`):

| Panel | Stratum | Conditions | Accelerated (shipped): tx % (gene %) | Ordinary: tx % (gene %) | Δ tx | Loci at budget (accelerated / ordinary) | EM time, ordinary ÷ accelerated |
|---|---|---:|---:|---:|---:|---:|---:|
| Ladder | κ 0.99 ON | 3 | 3.87 (1.38) | 4.39 (1.43) | +0.52 | 61 / 133 of 3,553 | 1.39× |
| Ladder | κ 0.99 OFF | 3 | 2.02 (0.25) | 2.29 (0.25) | +0.27 | 233 / 392 of 3,915 | 1.10× |
| Ladder | κ 0.50 OFF | 3 | 2.36 (0.31) | 2.59 (0.31) | +0.23 | 329 / 565 of 3,915 | 1.11× |
| Ladder | κ 0.50 ON (deferred) | 3 | 7.68 (3.78) | 8.04 (3.84) | +0.36 | 192 / 474 of 3,532 | 1.45× |
| Ladder | **all, g00 included** | 16 | 4.32 (1.56) | 4.62 (1.58) | +0.30 | 1026 / 1964 of 18,391 | 1.15× |
| RNA-short | κ 0.99 ON | 1 | 13.50 (1.34) | 13.99 (1.40) | +0.49 | 19 / 52 of 1,231 | 1.38× |
| RNA-short | κ 0.99 OFF | 1 | 3.55 (0.23) | 3.87 (0.21) | +0.33 | 99 / 142 of 1,293 | 1.14× |
| RNA-short | κ 0.50 OFF | 1 | 3.81 (0.27) | 4.01 (0.26) | +0.19 | 96 / 148 of 1,300 | 0.94× |
| RNA-short | κ 0.50 ON (deferred) | 1 | 19.67 (2.30) | 20.13 (2.43) | +0.46 | 28 / 76 of 1,221 | 1.85× |
| RNA-short | **all** | 4 | 10.73 (1.10) | 11.10 (1.15) | +0.38 | 242 / 418 of 5,045 | 1.23× |
| RNA-long | κ 0.99 ON | 1 | 5.17 (0.62) | 5.36 (0.62) | +0.19 | 18 / 40 of 1,302 | 1.36× |
| RNA-long | κ 0.99 OFF | 1 | 2.26 (0.19) | 2.41 (0.19) | +0.15 | 75 / 122 of 1,355 | 1.14× |
| RNA-long | κ 0.50 OFF | 1 | 2.17 (0.23) | 2.41 (0.23) | +0.24 | 71 / 130 of 1,355 | 0.99× |
| RNA-long | κ 0.50 ON (deferred) | 1 | 43.30 (11.50) | 43.31 (11.49) | +0.01 | 22 / 46 of 1,296 | 2.15× |
| RNA-long | **all** | 4 | 14.32 (3.42) | 14.46 (3.42) | +0.14 | 186 / 338 of 5,308 | 1.25× |

Each fl-gap stratum is one condition. The accelerated arm reproduces the release's shipped numbers (ladder κ 0.99 ON
3.87, RNA-short 13.50, RNA-long 5.17).

**Converged, on the deciding panels.** Tolerance 1e-12, at most 100,000 iterations, and every locus converged; two EM
threads, since a converged answer does not depend on the thread count. Transcripts % (genes %):

| Condition | Shipped: foundation → release | Converged: foundation → release |
|---|---:|---:|
| RNA-short g50 κ 0.99 ON | 13.32 (1.42) → 13.50 (1.34) | 13.32 (1.41) → 13.52 (1.33) |
| RNA-long g50 κ 0.99 ON | 4.25 (0.68) → 5.17 (0.62) | 4.12 (0.68) → 5.05 (0.62) |
| Ladder g05 κ 0.99 ON, release | 3.05 (0.58) | 3.00 (0.57) |
| Ladder g50 κ 0.99 ON, release | 4.46 (2.13) | 4.58 (2.13) |

The two capture-OFF ladder conditions (g05 κ 0.99 and κ 0.50) were still converging after 45 minutes and are not
in the table.
- **On realistic truth the stopping point barely matters.** Converging moves the deciding conditions by −0.13 to +0.12
  points, where it moved the test chromosome's by +3.5 to +5.8. It costs 2–3× the EM time.
- **The reader's two in-scope costs are the reader's.** They are the same in both views: RNA-short +0.18 shipped and
  +0.20 converged, RNA-long +0.92 and +0.93.

### 5.4 The reader's benefits under each solver view

Foundation against release, capture ON at κ 0.99 and 0.70, contaminated conditions. The accelerated and ordinary views
use one EM thread (`view3_compare.py`). The converged view is tolerance 1e-12 with at most 100,000 iterations. Each
cell is release minus foundation, in transcript points; "probed" means the probed genes' own transcripts.

| Panel | Accelerated (shipped): panel / probed | Ordinary updates: panel / probed | Converged: panel / probed |
|---|---:|---:|---:|
| Single probe | −0.73 / −33.1 | −0.69 / −32.8 | −0.70 / −33.1 |
| Sparse | −1.40 / −1.45 | −3.10 / −3.17 | −0.59 / −0.62 |
| Junction | −1.45 / −1.40 | −3.41 / −3.40 | −0.42 / −0.36 |

The absolute levels, foundation → release, transcripts % (genes %):

| Panel | Accelerated (shipped) | Ordinary updates | Converged |
|---|---:|---:|---:|
| Single probe | 8.63 (1.78) → 7.90 (1.08) | 6.78 (1.78) → 6.08 (1.09) | 12.39 (1.78) → 11.69 (1.08) |
| Sparse | 25.13 (3.74) → 23.73 (3.23) | 22.59 (3.73) → 19.49 (3.21) | 26.25 (3.76) → 25.66 (3.25) |
| Junction | 32.29 (1.35) → 30.85 (1.45) | 28.03 (1.36) → 24.62 (1.45) | 32.89 (1.35) → 32.47 (1.45) |

- **The benefits survive every view, with the same sign.** The gene-level gains are identical in all three views,
  because they are calibration effects: the reader decides how much of each gene is gDNA, and the EM's stopping point
  does not touch gene totals.
- **The single-probe gain does not depend on the view.** The probed gene's error is 35.4 % under the foundation and
  2.3–2.6 % under the release in every view.
- **The transcript-level gains scale with the stopping point.** On the sparse and junction panels, they grow under
  ordinary updates (−3.1 and −3.4) and shrink under convergence (−0.6 and −0.4).

### 5.5 Reading

**The reader.** Its benefits survive every solver view, with the same sign. The gene-level gains are identical in all
three views. The transcript-level gains on the sparse and junction panels range from −0.4 to −3.4 points with the
EM's stopping point and never reverse. Nothing here reopens the reader.

**The solver.**

- **The rounding instability belongs to the accelerator.** Ordinary updates are stable at every budget measured.
- **The shipped accuracy is set by how far the EM travels along directions the data barely determine.** The test
  chromosome's designed blocks hold near-identical isoforms at even truth, and there less travel scores better. On
  the ladder and the fl-gap panels, whose truth is skewed like real annotation, ordinary updates at the same budget
  score slightly worse: on 15 of 16 ladder conditions (−0.09 to +0.54 points) and on all eight fl-gap conditions
  (+0.01 to +0.49). An early stop's benefit is a property of the substrate, which is why the ladder decides.
- **A like-for-like repair buys stability only.** Ordinary updates run to the accelerator's distance (about 3,000)
  score the same on the test panel (8.12 % vs 8.34 %), at five times the EM time.
- **Convergence is nearly neutral where it counts.** On the deciding panels the converged answer is within
  −0.13 to +0.12 points of the shipped one. The test chromosome's large convergence penalty comes from its designed
  degenerate blocks.
- **The multimapper port's identical-paralog gate shows the same mechanism.** Two sequence-identical templates at
  gDNA 20, with equal true abundances, split 513 / 39 under the shipped EM: neither even nor a vertex, the "third
  behaviour" that gate exists to name. Converged, they split 552 / 0, a vertex. The total is the same either way.

**What it asks for.** This is a decision for the owner.

1. **0.8.0: keep the shipped solver,** and record the limitation with its magnitude. No stable alternative is better
   on the deciding panels: ordinary updates at the same budget are worse, ordinary updates matched to the
   accelerator's distance are five times slower, and convergence is nearly neutral there at 2–3× the EM time.
2. **0.8.1: regularize the poorly determined directions explicitly.** The drift is the estimator's variance there, so
   the principled remedy is explicit regularization with a stable solver run to convergence, decided on the ladder
   in both views. That keeps the reviewer's order: the numerical repair was investigated first (it is stable, not
   better), and statistical regularization comes next.

**The follow-up's evidence**, beside the tables below: `solver_replay.py`, `locus_trajectory.py`, `em_screen.py`,
`em_screen_report.py`, `view3_compare.py`, `converged_compare.py`, `review/solver_replay.json`, `review/em_screen/`,
the prototype builds `site_solverdiag` and `site_solverdiag_fnd`.

**Evidence** is under `$RIGEL_RUNS/prototypes/2026-10-13_modefree/`:

| What | Where |
|---|---|
| Expected oracle | `expected_bulk.py`, `review/expected_bulk.json` |
| Experiment 1 runs and results | `oracle_runs.py`; `out/obulk_*` (shipped EM) and `out/tight_*` (converged EM) |
| Experiment 2 runs and results | `stress_causality.py`, `review/stress_causality_{foundation,release}.jsonl` |
| Prototype site | `site_obulk` |
