# Ledger

**One row per step, appended as it lands, never retroactively.**

⭐ **This file holds work whose conclusions are still LIVE inputs.** On 2026-08-03 the index, suite, S5
and S6 entries (I1 through S5.g-2) moved to `docs/WIP_ARCHIVE.md`, verbatim — they are settled and
nothing reads their numbers any more. What remains is the solver-observables arc (whose P2 is still
gated OFF, waiting on C3) and the fragment-length / second-pass arc through **B4**, which is where the
current priority comes from. Plan: `docs/accumulator/IMPLEMENTATION_PLAN.md` §0.
Design: `docs/accumulator/DESIGN.md`. Deferred work: `TODO.md`.

The point of this file is **attribution**: a delta is only attributable if it is recorded against a
baseline taken from the same tree in the same session.

⛔ **The old benchmark baseline is VOID.** `r0 0.079005 / r3 0.046675` was the 32-condition
`ambig_dense_10mb` suite, deleted 2026-07-30. Do not quote it, compare against it, or try to reproduce it.

⭐ **THE FIRST BASELINE EXISTS as of S5.f (2026-07-30)** — eight chr22 pilot conditions, in the S5.f entry
below, bit-identical on re-run. ⚠ It carries A7's known **11.0 %** genome-wide gDNA over-call by
deliberate ruling, and it says nothing about dispersion (the simulator is Poisson by construction).

## The index — every entry, oldest first

⚠ **Split 2026-07-30 because this file passed 1,700 lines.** Entries up to and including the deletion
moved **verbatim** to `docs/WIP_ARCHIVE.md` — relocated, never rewritten, so the append-only property that
makes attribution meaningful survives. Everything from `I1` on is below in full.

| entry | what | where |
|---|---|---|
| S1 | index: reach semantics + the junction CSR | `docs/WIP_ARCHIVE.md` |
| S2 · S2.1 · S2.2 · S2.3 | the Python reference, the spec matrix, the doc corrections, the naming | `docs/WIP_ARCHIVE.md` |
| S3 | the C++ wired, byte-identical to the specification | `docs/WIP_ARCHIVE.md` |
| S4 | the payload, typed against the specification's vocabulary | `docs/WIP_ARCHIVE.md` |
| Index rebuild | the artifacts carry the S1 reach | `docs/WIP_ARCHIVE.md` |
| ⛔ Deletion | every benchmark and every index deleted; the baseline voided | `docs/WIP_ARCHIVE.md` |
| **I1** | an index can be rebuilt from its own manifest | below |
| **I2** | the human index rebuilt from scratch — and it is the SAME index | below |
| **I3 + B0** | the suite backbone, and the gate that says what it can judge | below |
| **B1 · B1a · B1b** | the capture simulator made fast, exactly | below |
| **S5.0** | ⭐ S5 STOPPED and turned into a derivation; the observable question answered | below |
| **S5.a** | ⭐ `length_sum` added to every population; `density` renamed `inv_length_sum` | below |
| **S5.b** | `fl.py` re-keyed to the five pure pools; the scan cache unblocked | below |
| **S5.c** | `effective_length.py` reduced to the ONE placements formula | below |
| **S5.d** | substrate collapsed to ONE type; the chain re-keyed to nodes/edges | below |
| **S5.e** | ⭐ A7 ruled; the 18 per-face arrays became 5; the faces dissolve | below |
| **S5.f** | ⭐⭐ **`calibrate()` RUNS — the FIRST BASELINE**; the third axis exported; the pooling dissolves | below |
| **S5.f-addendum** | the κ mirror is CONSISTENT — the inference is correct, only the label is wrong | below |
| **S6** | ⭐ goldens regenerated; two S5.f defects exposed; the EM's sampling mode measured | below |
| **S5.g-1** | the per-contiguous-edge RNA reach built + gated; ⛔ one open question gates the divisor | below |
| **S5.g-2** | ⛔ **the A7 taper A/B: it moves the gDNA fraction by ≤0.0002.** Fact 6's 11.0 % was bp-weighted geometry | below |
| **B2** | the pilot suite, and B0's verdict on it | below |
| **B3** | the scan cache — scan once, calibrate many times | below |
| **C2.6** | ⭐ **gap introns cut on EVERY fragment; the anchor's impossible tail is 85 % gone and the residual is LOCALISED to D3.** ⚠ D1 measured to cost more than it buys | below |

---


## P0 — the composition-evidence census: HOW MUCH MASS REACHES THE SOLVER BLIND (2026-07-31)

    Plan: `docs/calibration/SOLVER_OBSERVABLES_PLAN.md` §4   ·   Tool: `scripts/design/composition_evidence_census.py`
    Prediction recorded BEFORE the run: the ss0.50 conditions carry materially more no-evidence mass
    than the ss0.99 ones, concentrated on AMBIG slots.

⭐ **CONFIRMED, and the finding is larger than the prediction.** P0 exists to falsify P2's premise before
P2 is built. It did not falsify it.

**What is measured.** The mass-weighted share of unspliced fragment mass on chain slots with
`tau_lam <= 1e-9` and **no structural lock** — i.e. slots whose own composition belief carries
`own_precision = 0`, so their gDNA/RNA split is decided entirely by neighbour messages and the
population prior. ⚠ Structurally locked (intergenic pure-gDNA) nodes are EXCLUDED: they are composition
*certain*, not uninformed, and lumping the two would report a correct answer as a failure.

⚠ **`tau_lam = 0` is not "the strand likelihood was ignored".** The local ψ solve still sets the
*location* `f_g` from strand; what is zero is the *precision*, so the node emits nothing and its own
answer carries no weight in the fusion. The census measures blindness of the message system, not of the
per-node read-out.

### The result — 8 chr22 pilot conditions, no production code changed

| condition | κ | **no-ev** | node | edge | AMBIG | 1-str | lock |
|---|---|---|---|---|---|---|---|
| gdna100 ss0.50 capture_off | 0.5000 | **49.4 %** | 28.5 | 20.9 | 12.4 | 36.7 | 28.7 % |
| gdna100 ss0.50 capture_on | 0.4990 | **98.2 %** | 58.7 | 39.5 | 29.0 | 68.1 | 1.0 % |
| gdna100 ss0.99 capture_off | 0.0101 | **12.7 %** | 6.2 | 6.5 | 12.4 | 0.0 | 28.7 % |
| gdna100 ss0.99 capture_on | 0.0101 | **30.0 %** | 16.3 | 13.7 | 29.0 | 0.0 | 1.0 % |
| **none** × 4 (both κ, both capture) | — | **100.0 %** | 61.6–67.2 | 32.8–38.4 | 26.6–40.1 | 59.9–73.4 | 0.0 % |

Unstranded 86.9 % against stranded 60.7 % overall — but the pooled ratio understates it, because the
`none` arm is saturated at 100 % on both. **On the gdna100 arm, where the comparison is clean, it is
49.4 % vs 12.7 % and 98.2 % vs 30.0 % — 3.3× and 3.9×.**

### ⭐ Four findings, three of them not predicted

1. **The strand channel's on/off is exact, and visible in one column.** `no-evidence | single-strand` is
   **0.0 %** on both κ = 0.0101 conditions and **63.6 % / 98.9 %** on the κ = ½ ones. Every single-strand
   slot in a stranded library has own evidence; most have none in an unstranded one. This is
   `docs/SESSION_HANDOFF.md` §2's `I(f_g) ∝ (2κ−1)²` read off production output for the first time.
2. ⭐ **AMBIG mass is blind in EVERY condition, stranded or not**: `no-evidence | AMBIG` is **93.3 % –
   100.0 %** on all eight. That is the Schur gate (`node_init`, approach E) doing exactly what it was
   designed to do — a both-strand node's tilt is a free nuisance, so strand cancels out of `f_g` — and it
   means **13.3 – 40.1 % of library mass has never had own composition evidence in any library Rigel has
   ever run.** ⚠ This is not an unstranded-data problem. It is every both-strand locus, always.
3. ⭐ **A zero-gDNA library is 100 % blind at any κ.** `strand_evidence`'s derived noise floor is
   `σ²_d = ¼(1/N_rna + ω_r) + ¼(1/N_gdna + ω_g)`, so `N_gdna = 0 ⇒ σ²_d → ∞ ⇒ disc = 0`. Documented and
   deliberate, but the consequence had not been quantified: on the `none` arm **nothing** gives any node
   own evidence, and `none ss0.50 capture_off` nevertheless reports `f_gdna = 0.0793` against a truth of
   exactly 0 — 443,277 phantom fragments, produced entirely by messages and the prior.
4. ⭐ **Capture roughly doubles the blindness**, by removing the anchor rather than by adding noise:
   structurally-locked mass collapses **28.7 % → 1.0 %** because intergenic nodes are depleted, and
   no-evidence goes 49.4 % → 98.2 % (ss0.50) and 12.7 % → 30.0 % (ss0.99). The pure-gDNA anchor the
   prior-free pass leans on is a *capture-off* asset.

### The falsification — the census fires, and lands on the natural experiment

⛔ A census that does not move when the thing it claims to measure is switched off measures nothing
(`falsification_needs_perturbation`). Re-run `gdna100 ss0.99 capture_off` with **κ = 0.5 injected and
nothing else changed** (`--inject-kappa 0.5`, via `dataclasses.replace` on the condition's own fitted
`InjectedCalibrationPriors`):

| | no-evidence | no-ev \| single-strand |
|---|---|---|
| as fitted, κ = 0.0101 | 12.7 % | 0.0 % |
| **κ → 0.5 injected** | **49.4 %** | **63.6 %** |
| `gdna100 ss0.50 capture_off`, the natural experiment | **49.4 %** | **63.6 %** |

⭐ **The injected arm reproduces the natural experiment to the last digit** — same payload, one variable,
landing on the independently-simulated condition's value. The census tracks the strand channel and
nothing else.

### What this decides

* **P2 is justified**: the length likelihood is the only proposed source that is independent of κ, of
  `N_gdna`, and of the single-strand/AMBIG gate. Finding 2 is the strongest case for it — the Schur
  argument that silences strand on AMBIG nodes does not apply to a channel that has no dependence on the
  tilt `θ` at all.
* ⚠ **P2's gate must not be scored on the `none` arm alone.** All four zero-gDNA conditions are saturated
  at 100 % blind, so any change scores "better" there (trap 19, one-sidedness). The clean comparison is
  the gdna100 arm, and `gdna100 ss0.50 capture_on` (98.2 % blind, `f_gdna` 0.3754 against ~0.50 — the
  row S5.f recorded as "unexplained") is the sharpest single target in the suite.
* **`docs/calibration/SOLVER_OBSERVABLES_PLAN.md` §4 updated** with the result. P1 (the units fix) is unaffected by this
  measurement and remains next.

---

## P1 — the EM prior is a FRAGMENT COUNT, not an object-incidence sum (2026-07-31)

    Plan: `docs/calibration/SOLVER_OBSERVABLES_PLAN.md` §5   ·   Unit gates: `tests/calibration/test_prior_units.py`
    End-to-end gate: `scripts/design/prior_units_check.py`

⭐ **A UNITS BUG, WITH A PROOF — not a modelling change.** ``assemble_priors`` handed the EM two additive
pseudocounts that are added straight to its own fragment counts (``G = n_gdna + a_g``,
``em_solver.cpp:apply_grouped_prior_update``), but built them by **summing per-object masses**. A fragment
deposits on ``max(K, 1)`` objects, so for a partition of spacing ``s``::

    incidences(w) = max( 1 , (w-1)/s )

Counts are conserved exactly where every node is longer than every fragment and become a
**length-weighted** count where they are not — and the weight is the fragment's own length, so it does
not cancel between two components with different mean lengths.

**The fix, in one line of model:** density is intensive, so pool it as a ratio of sums and integrate it
over the span. ``rho_c = Σ share·mass_c / Σ share·A_c``, then ``prior_c = rho_c · span_bp`` — the **same
genomic span for both components**, so the ratio carries no length tilt. ``A_r`` did not exist:
``CalibrationResult`` gained ``rna_node_eff_len`` / ``rna_edge_eff_len``, **projected** off
``NodeGeometry.eff_rna`` (never recomputed — trap 27) exactly as their gDNA twins are.

### The gates, written first and verified failing

⭐ The failure pattern was **self-confirming**, which is why it was worth writing them deterministically
rather than end to end:

| gate | before | after |
|---|---|---|
| partition invariance, 1200 bp re-tiled | ✅ passes at 1200 and 400 bp nodes, ⛔ **fails at 100 bp** | passes at all four tilings |
| …the RNA prior across that re-tiling | **50.05 → 109.45, a 2.19× swing** with the library unchanged | invariant |
| prior == true fragment count | 34.53 against 36.0 | exact |
| g:r flat in ``mu_g/mu_r`` (swept 0.25×–2×, both directions) | ⛔ fails everywhere **except exactly ``mu_g = mu_r``** | flat |

Invariance holding at 1200/400 bp and breaking at 100 bp is precisely the ``max(1, (w-1)/s)`` crossover;
the ratio test passing at its own null (equal means ⇒ no tilt) identifies the mechanism as the length
tilt and nothing else. gDNA is *partition-invariant in that fixture* because its 50 bp fragments fit
inside every node — the swing is entirely on the component whose fragments do not.

### ⛔ Perturbation P1e found a real defect, and the test was the weak one

Five perturbations. Four were caught (RNA÷gDNA opportunity, dropped ``span``, reverted raw sums, mean of
ratios). **P1e — flooring the divisor to ``1e-9`` instead of testing ``support > 0`` — left all 15 tests
green**, because the only zero-support fixture also had zero mass.

⚠ **``mass > 0`` with ``support == 0`` is an ordinary configuration.** ``contained_eff_length`` is exactly
0 wherever a node is shorter than that component's shortest fragment; measured on the chr22 pilot index
against its own measured pure pools, that is **21.7 % of nodes for RNA and 18.7 % for gDNA**
(`docs/SESSION_HANDOFF.md` §1 fact 9 records 12.4 % genome-wide). The solver can still put mass there — ``f_g``
is an inference, not a fact. The floor would turn one stray fragment on a 40 bp node into a density of
~1e9 (trap 23, the mechanism that once seeded false gDNA into neighbouring exons).

⭐ **And behind the weak test was a defect in the fix itself**: leaving that mass in the pooled numerator
while it contributes nothing to the denominator inflates ``rho`` with no exposure to pay for it. **Both
sides of a pooled rate, or neither** — ``_mass_where_there_is_opportunity`` now drops it from both, the
same contract ``node_geometry._rate`` already keeps on the solver side. New test
``test_mass_on_a_zero_opportunity_object_is_dropped_from_BOTH_sides``; P1e now fails it.

### ⛔ T3 caught a hard violation: the old prior EXCEEDED the whole library

`prior_units_check.py` computes both arms from ONE ``CalibrationResult`` — no re-scan, no second code
path, no simulation noise:

| condition | fragments | A: raw sum | **B: P1** | A/frag | **B/frag** | A f_gdna | **B f_gdna** |
|---|---|---|---|---|---|---|---|
| gdna100 ss0.50 capture_off | 9,634,502 | 6,386,844 | 5,660,367 | 0.663 | 0.588 | 0.4390 | 0.4616 |
| gdna100 ss0.50 capture_on | 9,775,761 | 10,466,755 | 9,190,222 | ⛔ **1.071** | **0.940** | 0.4085 | 0.4285 |
| gdna100 ss0.99 capture_off | 9,633,129 | 6,382,577 | 5,650,154 | 0.663 | 0.587 | 0.3863 | 0.4080 |
| gdna100 ss0.99 capture_on | 9,774,382 | 10,457,211 | 9,268,515 | ⛔ **1.070** | **0.948** | 0.5360 | 0.5556 |
| none ss0.50 capture_off | 4,636,413 | 3,799,206 | 3,247,412 | 0.819 | 0.700 | 0.1167 | 0.1196 |
| none ss0.50 capture_on | 4,791,669 | 4,459,761 | 3,784,949 | 0.931 | 0.790 | 0.0126 | 0.0130 |
| none ss0.99 capture_off | 4,636,640 | 3,797,474 | 3,239,915 | 0.819 | 0.699 | 0.0044 | 0.0045 |
| none ss0.99 capture_on | 4,792,251 | 4,460,570 | 3,784,966 | 0.931 | 0.790 | 0.0020 | 0.0021 |

⭐ **On both `capture_on` gdna100 conditions the pre-P1 prior totalled MORE pseudocounts than the library
has accepted fragments** (1.071× and 1.070×) — for a quantity that arbitrates only the *unspliced*
subset. That is structurally impossible and nothing was checking it. Every arm-B row is now below 1.

⚠ **The f_gdna column is the LOCUS-PROJECTED prior ratio, not `LEDGER`'s library-wide scalar** (intergenic
nodes are dropped by the projection, and they are gDNA-rich — hence 0.4390 here against the S5.f table's
0.4998 on the same condition). The two are not comparable and must not be quoted against each other.

### ⚠ The honest scoring, including where it gets worse

The correction moves gDNA **up** on every condition, as §5.4 predicted before the run. On the zero-gDNA
arm — truth exactly 0, and the one arm whose truth survives the projection — that is **slightly away from
truth**: 0.1167 → 0.1196 on `none ss0.50 capture_off`, and +0.0001 to +0.0004 on the other three.

⛔ **This is trap 19 in reverse and must not be read as a regression.** On a library with no gDNA, any
bias that under-calls gDNA scores better; the raw sum's under-call was accidentally flattering exactly
that arm. P1 is a units error with a proof, not a knob tuned to a scenario, and the residual 0.1196 on a
zero-gDNA library is **P2's target** — it is the unstranded condition P0 measured at **100 % blind**.

### Suite

**1757 passed / 22 failed**, against HEAD's **1758 / 1** re-measured in the same session by `git stash`
(trap 17). Accounting is exact: `1758 − 21 golden + 20 new = 1757`. The 21 are `test_golden_output.py`
moving numerically, which is expected and deliberate — ⚠ **the goldens regenerate ONCE, after P2.** The
1 pre-existing failure (`test_nrna_double_counting`, `TODO.md` §7) is byte-identical.

`tests/calibration/` **569 passed** (549 baseline + 16 new + 4 new schema parametrizations, the latter
because `rna_*_eff_len` were added to the shape/sign/dtype guard lists rather than left ungated).

### ⚠ Unrelated finding: commit 49e9b456 is incomplete

`S6 + S5.g-1` committed `splice_graph.build_contiguous_edge_reach_arrays` and `test_edge_reach.py` but
**left `node_geometry.build_node_geometry`'s `edge_rna_reach` parameter uncommitted**. At HEAD the
builder therefore has no consumer — the A7 switch exists only in the working tree. Not caused by this
work and no data was lost; flagged because a builder whose only caller is uncommitted reads as dead code
to the next reader.

---

## P2 — the fragment-length likelihood: BUILT, MEASURED, and ⛔ **GATED OFF — it is blocked upstream** (2026-07-31)

    Plan: `docs/calibration/SOLVER_OBSERVABLES_PLAN.md` §6   ·   Module: `rigel.calibration.length_likelihood`
    Gates: `tests/calibration/test_length_likelihood.py` (43)   ·   A/B: `scripts/design/length_likelihood_ab.py`
    Switch: `CalibrationConfig.length_likelihood`, **default False** — the False arm is byte-identical to P1

⭐ **THE MECHANISM WORKS, EXACTLY AS DESIGNED.** P0's headline was that 13.3–40.1 % of library mass
reaches the solver with no own composition evidence in every condition, and **100 %** on an unstranded or
zero-gDNA one. Turning this channel on:

| no-evidence mass share | OFF | **ON** |
|---|---|---|
| gdna100 ss0.50 capture_off | 49.4 % | **0.0 %** |
| gdna100 ss0.50 capture_on | 98.2 % | **3.9 %** |
| gdna100 ss0.99 capture_on | 30.0 % | **1.1 %** |
| **none** × 4 | 100.0 % | **0.0 %** |

The blind mass is gone. The Schur argument that silences strand on AMBIG nodes does not apply to a term
with no ``θ`` dependence, and the measurement confirms it.

⛔ **AND THE ANSWER IT GIVES IS WRONG — 15× worse on the only arm with an unambiguous truth.**

| condition | truth | OFF | **ON** |
|---|---|---|---|
| none ss0.99 capture_off | **0.0000** | 0.0030 | **0.1269** |
| none ss0.99 capture_on | **0.0000** | 0.0016 | **0.2664** |
| none ss0.50 capture_off | **0.0000** | 0.0793 | **0.4689** |
| none ss0.50 capture_on | **0.0000** | 0.0101 | **0.6156** |
| gdna100 ss0.50 capture_off | ~0.52 | 0.4998 | 0.6064 (overshoots) |
| gdna100 ss0.99 capture_on | ~0.52 | 0.4897 | 0.4264 (worse) |

Zero-gDNA arm mean ``|f_gdna|``: **0.0235 → 0.3695**.

### ⭐ THE CAUSE IS UPSTREAM, AND IT IS NOT IN THE LIKELIHOOD

⚠ **CORRECTION, 2026-07-31, same session — the first diagnosis below was half wrong and the owner caught
it.** "The pmf is FABRICATED" is not what happens. ``build_fl_models`` EB-shrinks each pool toward the
global histogram, and at ``n_gdna = 0`` it shrinks **all the way**: verified, ``gdna_pmf == global_pmf``
byte-identically on both zero-gDNA conditions. That is exactly the intended behaviour, and the owner's
follow-on is also right — **identical pmfs for the two components make this channel exactly inert**, which
is proven byte-identically end to end (`test_identical_components_are_EXACTLY_inert_not_nearly_inert`).

⭐ **So the real defect is narrower and sharper: on a zero-gDNA library ``global`` and ``rna_fl_pmf``
describe the SAME population — every fragment is RNA — and they disagree.** Measured on
`none ss0.99 capture_off`: global **210.1 ± 85.5** against the RNA pool **234.5 ± 146.2**, i.e. **+11.6 %
in the mean and +71.1 % in the sd**, against `docs/calibration/S5_DESIGN_LOG.md` §3.6's independently-predicted +14 % / +50 %
from junction opportunity. The length likelihood then reads that disagreement as composition, correctly —
there is no defect in the likelihood at all. ⛔ **The channel is a comparator, and it is being handed two
rulers of different length.** The original text follows.

**The gDNA fragment-length pmf reads a fallback on a zero-gDNA library.** Measured:

| condition | gDNA pool n | gdna_fl_pmf mean ± sd | rna_fl_pmf mean ± sd |
|---|---|---|---|
| gdna100 ss0.50 capture_off | 4,526,536 | 194.4 ± 97.4 | 234.9 ± 146.8 |
| **none** ss0.50 capture_off | **0** | **210.1 ± 85.5** | 234.7 ± 146.7 |
| **none** ss0.50 capture_on | **0** | **211.0 ± 85.9** | 242.3 ± 148.0 |

⛔ **Zero observations, and `build_fl_models` still returns a confident pmf** — a fallback that happens to
differ from the RNA pmf, so the likelihood dutifully discriminates against a component that does not
exist. **Every other consumer uses `gdna_fl_pmf` as a DIVISOR, where a wrong pmf is a scale error. This is
the first consumer that uses it as a DISCRIMINANT, where a fabricated pmf manufactures composition out of
nothing.**

⭐ **`strand_evidence` already guards exactly this case and the length channel does not.** Its derived
noise floor is `σ²_d = ¼(1/N_rna + ω_r) + ¼(1/N_gdna + ω_g)`, so `N_gdna = 0 ⇒ disc = 0` — the strand
channel refuses to speak when there is no gDNA to calibrate against. The length channel needs the same
principle, and it is a **derivation**, not a flag: the discriminant is the separation between two fitted
pmfs, and it must be gated by the sampling uncertainty of that separation.

⚠ **And a second upstream defect, which the gdna100 overshoot points at.** `docs/accumulator/DESIGN.md` §8.1(b)
records that **every pool histogram is length-biased by its own opportunity** and must be divided by it
before use; `docs/calibration/S5_DESIGN_LOG.md` §3.6 measured the effect on this very suite (the RNA pool reads 234.9
against a configured 206.1 — **+14 % in the mean, +50 % in the sd** — and the gDNA intergenic pool reads
195.4 against 156.5, *unexplained*). That correction is listed under §4 **"Not yet decided"**. So the two
pmfs are already tilted, by different and partly unknown opportunities, and `length_likelihood` tilts them
a second time by ``A(w)``. ⛔ **§8.1(b) is no longer a tidy-up. It BLOCKS P2.**

### ⭐ Three defects found by the falsification discipline, two of them in my own work

| | found by | what |
|---|---|---|
| **P2d** | perturbation | deleting the ``−½ log det`` term left all 41 tests green. Its grid-variation is ``O(1)`` while the quadratic's is ``O(N)``, so it is negligible at depth and **decisive at the median node**: the peak moves **0.32 at N=1**, 0.05 at N=5, 0.004 at N=50. ⚠ It also documents the model's own limit — the Gaussian is asymptotic in ``N`` and at ``N=1`` is not a trustworthy *location*, because the true single-draw mixture likelihood is bimodal |
| **the equal-pmf null** | A/B | forcing ``rna_fl_pmf = gdna_fl_pmf`` must be byte-identical and was not (0.5008 vs 0.4924) |
| **⭐ the grid-width leak** | chasing that null | rows that are flat *to 1e-11* rather than to 0 pass `density_factor_precision`'s ``ptp > 1e-12`` gate; the near-uniform posterior then returns ``tau = 1/Var(uniform over λ) = 0.029016`` — **the grid's own width sold as composition evidence**, which is exactly what that function's docstring says must never happen. Measured: **689 slots, max tau 0.02902**. Fixed by gating on the MOMENTS structurally, so the null is now *exactly* zero, and pinned by an **exact-equality** assertion (``assert_allclose`` passes at 1e-11 and would miss it) |

⚠ **The grid-width leak is latent in `density_factor_precision` itself, not only in this caller.** Any
factor that is nearly-but-not-exactly flat gets 0.029 of free precision from the grid. The intron factory
is safe today only because its inactive rows are written as exact zeros. A threshold-free fix exists —
report ``max(0, 1/Var_post − 1/Var_uniform)``, the precision the factor ADDS over the grid's own width —
but it changes a shipped, default-on feature and needs its own A/B. **Recorded, not smuggled in.**

### What landed, and what did not

**Landed** (behind `length_likelihood=False`, byte-identical when off): `calibration/length_likelihood.py`
(the tilted moments in closed form — ``O(n_nodes)`` cumulative sums, no ``n_nodes × max_len`` array, which
would be 8 GB at human scale); the two length channels gathered onto `NodeGeometry`; a `length_loglik`
argument threaded through `simplex_logodds`'s 1-D and 2-D solves, `node_init` (as **ungated** source 4) and
`bp_solver`'s local *and* final solve; `_build_length_loglik` in `calibrate`.

**Verified**: U1 against **exact enumeration over integer start positions** (a literal loop, not the
module's own closed form — trap 1) over 4 pmfs × 6 node lengths, both frames; ``q12 ≡ 1`` at a node as a
structural identity of the ``1/L`` deposit weight; ``moments.eff`` byte-identical to
`contained_eff_length`/`crossing_eff_length` (trap 27); 5 perturbations, all now caught.

**Suite**: **1800 passed / 22 failed** — the same 21 golden moves from P1 plus the 1 pre-existing
`TODO.md` §7 failure; nothing new. `tests/calibration/` 654.

### ⛔ THE VERDICT, AND IT IS NOT MINE TO MAKE

P2 is **implemented and inert**. It must not be switched on until the fragment-length pools are trustworthy
as a *discriminant*:

1. **`build_fl_models` must not fabricate a pmf from an empty pool.** A zero-observation pool has no
   distribution, and the honest output is "no gDNA length model", which must make the length channel inert
   — the same statement `strand_evidence`'s `1/N_gdna` already makes for strand.
2. **`docs/accumulator/DESIGN.md` §8.1(b)'s per-pool opportunity correction must land** — divide each pooled
   histogram by its own opportunity before normalising. Listed as undecided since S5.b; it is now the
   blocker.
3. Only then re-run this A/B. The before-picture (P0's blindness census) and the after-picture (this
   entry's no-evidence collapse to 0.0 %) are both recorded, so the re-run is a two-command comparison.

⚠ **Do not "fix" this by damping the channel.** The channel is doing what it was built to do; it is being
fed a length model that does not describe the library. Damping would hide the upstream defect behind a
tuned constant, which is the failure mode `docs/SESSION_HANDOFF.md` §3 trap 12 records three times over.

---

## C0 — ⭐ the accumulator's `L` is PROVEN, and it is fit to be the gold standard (2026-07-31)

    Audit: `docs/accumulator/FRAGMENT_LENGTH_AUDIT.md`   ·   Gate: `tests/native/test_fragment_length_proof.py`
    Owner precondition, 2026-07-31: *"the new fragment length computation is the newest implementation
    and I'm not sure how rigorously it has been tested in all cases; we need to prove that very
    carefully if it's going to become the gold standard."*

⛔ **The precondition was right, and the existing coverage did not meet it.** `test_accumulator_spec`
pinned **six hand-picked malformed intron lists**, chosen by the same author as the code they check —
the coverage pattern that finds what the author thought of and nothing else.

### The oracle is a different ALGORITHM, not a different spelling

`docs/SESSION_HANDOFF.md` §3 trap 1. The oracle is **integer set arithmetic**::

    covered = set(range(start, end))  −  ∪ set(range(a, b)) for each intron

No sorting, no merging, no cursor walk, no ``searchsorted``. Every malformed case the production path
handles by explicit logic — overlapping, nested, abutting, duplicated, zero-length and out-of-range
introns — the oracle handles **by construction**, because set subtraction is idempotent and order-free.

⭐ **And all four deposit populations fall out of that one set**, which is what makes this a proof of the
**geometry** and not only of ``L``:

| | oracle |
|---|---|
| ``L`` | ``len(covered)`` |
| line ``p`` crossed | ``p−1 ∈ covered and p ∈ covered`` — `docs/accumulator/DESIGN.md` §2's definition, verbatim |
| node spanned | ``cuts[i]−1 … cuts[i+1]`` all covered |
| node contained | no junction used, and ``min``/``max`` of ``covered`` fall in one node |

⚠ §3.1 requires that *"whatever counts toward ``L`` must also count as coverage for crossing"*. **Nothing
tested the two against each other before.** They now come from one set, so they cannot disagree.

### Coverage

| | |
|---|---|
| **exhaustive**, 0 introns | 78 configurations |
| **exhaustive**, 1 intron | 7,098 |
| **exhaustive**, 2 introns | **326,508** |
| randomised, ≤4 introns, realistic coordinates (9 kb ref, 3 kb spans, 900 bp introns), fixed seed | 4,000 |

**333,684 exhaustive configurations** — the complete space over a 12 bp reference whose 4 nodes include a
**1 bp** node — plus a randomised sweep at a scale exhaustive enumeration cannot reach (a 300 bp molecule
spanning a long intron, the case §3.2's length limit turns on). **All pass.**

### ⭐ The proof was proven to FIRE — 7 perturbations

| | perturbation | |
|---|---|---|
| **L1'** | ``L`` from ``span − Σ(RAW intron lengths)`` — the formula §3.3 says goes **negative** on a wide overlap | ✅ caught |
| **L2'** | introns not clipped to the fragment | ✅ caught |
| **L4** | fragment not clipped to the reference | ✅ caught |
| **L5** | crossing boundary ``searchsorted`` side flipped right→left | ✅ caught |
| **L6** | spanning loop credits one node too many | ✅ caught |
| **L7** | containment keyed on the fragment EXTENT, not its first/last COVERED base | ✅ caught |
| **L3** | ABUTTING introns no longer merge | ⚠ **not caught** |

⚠ **L3 is a correct non-failure.** Merging abutting introns changes only the ``introns_absorbed`` QC
counter — set subtraction removes ``[10,20) ∪ [20,30)`` and ``[10,30)`` identically, so ``L`` and all four
populations are untouched. The counter is pinned separately by ``test_accumulator_spec``.

⭐ **And a first attempt at L1 was ALSO a no-op**, which is worth recording because it is the design's own
claim demonstrated: replacing ``Σ segments`` with ``span − Σ introns`` **after** ``_normalise_introns``
changes nothing, because normalisation is precisely what makes the two formulas agree. The bug §3.3 warns
about only exists on the RAW list — and that version is caught.

⚠ **An error the oracle caught in ME first**: my first oracle reported crossed lines by *cut* index where
the deposit writes ``edge_base + line − 1``. Writing the oracle independently forced the offset to be
re-derived rather than copied, which is the entire point of trap 1.

### ⛔ Scope, stated rather than implied

The fixture carries **no annotated junctions**, so every fragment is unspliced and the **spliced routing**
(``edge_spliced`` vs ``edge_unspliced``, junction credit, the containment block) is not exercised here;
``test_accumulator_spec`` covers it. That is the right scope for C1: the unconditional histogram bins by
``L``, and ``L`` does not depend on which population the fragment is routed to.

⚠ Proving the Python reference proves the C++ too — `test_accumulator_native_parity` gates the C++ on
byte-identity to it, so a defect here would have been reproduced faithfully in both.

### Verdict

✅ **`L` is fit to be the tool's one definition of fragment length. C1 may proceed.** Suite **1805 / 22**
(the 21 golden moves + the pre-existing `TODO.md` §7 failure; +5 new proof tests).

---

## C1 — ⭐ the accumulator gets an unconditional length histogram, and the frame mismatch is GONE (2026-07-31)

    Audit: `docs/accumulator/FRAGMENT_LENGTH_AUDIT.md` §4 C1   ·   Gates: `tests/native/test_fragment_length_proof.py`,
    `tests/test_accumulator_payload.py`   ·   Precondition: C0 (`L` proven) — done first, by owner ruling

**What.** One `uint32[max_length + 1]` row, `deposited_lengths`, incremented on the same line as
`node_start_count` and the `DEPOSITED` counter. The five pure pools are untouched — `_pool()` still
returns `None` for a mixture, because an impure pool is worse than a missing one (`docs/accumulator/DESIGN.md`
§8). ⭐ **This is a separate tally, not a sixth pool**, and that distinction is the whole point.

**Why.** The accumulator binned only *conditioned* pools, so it had no unconditional histogram; the
empirical-Bayes shrinkage needs one; so the anchor was taken from the **scanner**, which measures length
by two other rules over another population. C1 removes the reason that borrowing existed.

⭐ **It is "unconditional GIVEN DEPOSIT", and the name says so.** It excludes what the accumulator
rejects — too long, ambiguous path, strand-undefined, empty — each of which is counted in `qc`. That is
**exactly** the population the pools are drawn from, which is what makes it the right anchor rather than
merely a convenient one; an anchor over a *wider* population would re-create the frame mismatch in a new
place.

### ⭐ THE MEASUREMENT: C1 separated the two mechanisms the audit could only assert were separable

On the zero-gDNA conditions every fragment is RNA, so the anchor and the RNA pool describe **one
population** and any gap between them is bias:

| | mean | sd |
|---|---|---|
| **the old scanner anchor** | 210.1 | 85.5 |
| ⭐ **C1, the accumulator's own** | **218.0** | **111.2** |
| the RNA_SPLICED pool | 234.7 | 146.7 |

| gap to the RNA pool | before C1 (scanner anchor) | **after C1** | closed |
|---|---|---|---|
| mean | **+11.6 %** | **+7.7 %** | **34 %** |
| sd | **+71.1 %** | **+32.0 %** | **55 %** |
| support ceiling | 713 vs 1000 | **1000 vs 1000** | ⭐ **entirely** |

⭐ **The frame component is gone and the opportunity component is now isolated and measurable.** §2 of the
audit named two mechanisms and could only argue they were separable; C1 separated them. The residual
**+7.7 % / +32 %** is the junction-opportunity tilt, and it is C3's (§8.1(b)'s) target — now against a
same-frame reference instead of a confounded one.

⚠ The old anchor's sd was too *small* (85.5) precisely because definition **B** is transcript-space and
unanimity-gated and its support stopped at 713 bp. An anchor narrower than the thing it anchors is the
worst direction for an EB prior: it shrinks the pools toward a distribution that cannot produce them.

> ⛔ **CORRECTION, added 2026-08-01 by C2.6 — THIS ENTRY READ THE SUPPORT CEILING BACKWARDS.**
> The "713 vs 1000 → 1000 vs 1000, ⭐ entirely" row above is **wrong, and in the reassuring direction**.
> **713 bp is the library's TRUE maximum**, read from `truth_fragment_lengths.tsv`; **1000 is
> `max_frag_length`, the clamp**. Definition **B** took its length from the *transcript*, so it could not
> produce an uncut intron — on the ceiling it was right and `L` was wrong, by an intron that was never
> sequenced. ⚠ The paragraph immediately above ("the old anchor's sd was too small... its support
> stopped at 713 bp") is the same error stated a second way. ⭐ The residual "+7.7 % / +32 %" was
> therefore **not one mechanism**: the mean is the junction tilt (C3's) and the sd was the uncut intron,
> now fixed — **the anchor's sd against truth is +1.98 %, from +26.97 %.** See the C2.6 entry.
> ⚠ None of this undoes C1 or C2, whose other reasons all still hold.

### Gates

| | |
|---|---|
| **G1 byte-identity** | the parity gate reads its field list off `dataclasses.fields(Tally)`, so `deposited_lengths` joined it **automatically**, dtype and shape included. ⭐ Proven to fire: three reference perturbations (stop binning · bin one short · **bin at the SPAN instead of `L`**) all caught |
| **G2 the invariant** | `Σ deposited_lengths == Σ node_start_count == qc.deposited`, incremented on one line so the three cannot drift by construction. Same externally-checkable form as §10.2's start count and a **different statement**: that one says every fragment was located in space, this one that every fragment was binned by length |
| **at the payload door** | `from_scan_result` refuses a histogram that does not sum to `qc.deposited`, or is the wrong length. ⚠ It must live at the boundary and not only in the accumulator's tests, because the payload is what a **cached** scan is rebuilt from |
| **superset** | every pooled fragment is also binned here at the same `L` — the property that makes it a usable anchor. Strict superset in general: an exonic contained fragment and a multi-line crossing enter no pool at all |
| **rejected ≠ binned** | too-long, ambiguous-path and empty fragments bin nothing |

⭐ **6 perturbations, and the 6th found a hole in the tests rather than the code**: deleting the payload's
door-check left everything green, because nothing asserted the check *existed*.
`test_a_deposited_lengths_HISTOGRAM_THAT_DOES_NOT_BIN_EVERY_FRAGMENT_IS_REJECTED` closes it and now fails
without the check.

### Blast radius, all as designed

* `payload_schema_digest` **66b41ea0b645209d → b7d29676c58b2c65**, so every existing scan cache is refused
  at load. ⭐ That is exactly what that key exists for (`scan_cache` docstring: S5.a's `length_sum` is the
  precedent). The 8 pilot caches were rebuilt — **56.0 s**, and all 8 satisfy G2.
* Cache serialisation is driven by `dataclasses.fields(AccumulatorPayload)`, so the new row persists with
  **no change to `scan_cache`**.
* ⚠ **A pre-existing defect fixed in passing**: `Tally` declared `edge_spliced_length_sum` and
  `sj_length_sum` **twice each**. Harmless (the dataclass de-duplicates to 18 fields) but it is a wart in
  the file that IS the executable specification, and it was committed.

**Suite 1809 / 22** — the 21 golden moves + the pre-existing `TODO.md` §7 failure; +4 C1 tests. Native
96/96 including parity and worker determinism.

⛔ **`build_fl_models` still reads the scanner anchor.** C1 makes the correct anchor *exist*; **C2** is
what switches the consumer over and deletes the scanner histogram. Nothing downstream has moved yet, by
design — this step is additive and every number in the tool is unchanged.

---

## C2 — ⭐ ONE definition of fragment length, and the scanner's histogram is DELETED (2026-08-01)

    Audit: `docs/accumulator/FRAGMENT_LENGTH_AUDIT.md` §4 (C2.0–C2.5)
    Gates: `tests/test_one_fragment_length_definition.py` (the standing grep), `tests/test_d7_transcript_eff_lengths.py`,
           `TestSpliceCensus` + `TestFragmentLengthAnchor` in `tests/test_scanner_accumulator_integration.py`
    Measurement: `scripts/design/fl_anchor_gap.py` (new)   ·   Report A/B: this tree vs a worktree at `d045d820`

⭐ **The tool now has ONE definition of fragment length — the accumulator's `L` — measured in ONE place,
over ONE stated population.** D1, D2, D3 and D5 are closed; D6 is closed by deletion (the population that
was silently dropped no longer exists as a concept); D7 is **verified rather than assumed**, which the
audit demanded in those words.

### ⛔ THE OWNER'S RULING REJECTED THE PREMISE OF BOTH OPTIONS, AND IT WAS RIGHT

C2 was blocked on one decision: `rigel report`'s five per-fragment splice-type counts had no accumulator
equivalent. Option **(a)** grew five counters *in the accumulator*; option **(b)** reshaped the report and
lost two of five categories. The recommendation was (a). The owner ruled **neither**:

> *"There are some QC counts that the scanner must be responsible for. They truly live in the scanner …
> If the QC counts are generated in one place and have no algorithmic use, there's no need to pass them.
> … I don't see a reason to propagate artifacts into the accumulator, it's the scanner's job to identify
> and filter these out."*

⭐ **The principle the audit was missing.** §4's "ONE definition, ONE place, ONE population" governs
**model inputs**. The five splice counts are not a model input — they are a **census of the scanner's own
classification decisions**, with no consumer but the report. They were entangled with the histogram only
because they were *derived from it* (`category_models[stype].n_observations`), an implementation accident.
So C2.0 is not "move them to the accumulator"; it is **sever them from the histogram and leave them where
they are generated**.

⛔ **And option (a) would have been wrong on its own terms.** `SPLICE_ARTIFACT` fragments never reach the
accumulator — the deposit adapter returns early, because a blacklisted junction's span derives from an
alignment the scanner has already refused to believe — so (a)'s stated gate, *"the five counts sum to
`qc.deposited`"*, **cannot hold**. That was found by pricing the option, not by building it.

⭐ **What the ruling bought, measured**: `payload_schema_digest` stays **`b7d29676c58b2c65`**. No accumulator
change, no reopened S3 byte-identity gate, and **all 8 pilot scan caches remained valid** — against option
(a), which would have moved the digest and rebuilt every one.

### C2.0 — the splice census, and an identity that closes the books across two subsystems

One `std::array<int64_t, NUM_SPLICE_TYPES>` in `BamScanStats`, incremented at **one** site (the top of the
deposit adapter, before the artifact hold-out), plus **one** counter `n_deposit_not_offered`.

⚠ **An array, not five named fields**, so `census[st]++` is one statement, the merge one loop, the export
one loop over the *existing* `splice_type_label` table — and a category added to `SpliceType` cannot slip
past any of the three. ⭐ **No name table anywhere**: `splice_type_label`'s strings are exactly the
`SpliceType` member names lower-cased, so `rigel.splice.census_field` derives the same key on the Python
side. The correspondence is **tested**, not asserted — which matters because `_apply_scan_stats` copies with
`dict.get(key, 0)`, so drift reads as *zero*, not as an error.

⭐ **The gate as written in §4 was unsatisfiable; the replacement is stronger.** In C1's G2 form:

    Σ census − census[SPLICE_ARTIFACT] == qc.deposited + Σ qc.dropped_* + n_deposit_not_offered

⭐ **And the artifact count is derived a SECOND way, by a subsystem that has never heard of an artifact.**
Scan one BAM against the same index with and without a splice blacklist: blacklisting relabels fragments
`SPLICED_ANNOT → SPLICE_ARTIFACT`, so `qc.deposited` — which knows nothing of splice types or blacklists —
must fall by **exactly** the artifact census. Measured 200 → 40 against a census of 160. `docs/SESSION_HANDOFF.md`
§3 trap 1: a check that re-derives the number by the same route checks nothing.

**6 perturbations, all caught — but TWO needed fixtures built for them, and that is the finding:**

| | perturbation | |
|---|---|---|
| **PC1** | one category never censused | ✅ caught by 3 tests |
| **PC2** | census moved AFTER the artifact hold-out | ✅ caught — ⚠ **only by the blacklist fixture** |
| **PC3** | artifacts counted as `not_offered` as well | ✅ caught |
| **PC4** | one C++ export key renamed (the `.get(key, 0)` silent-zero) | ✅ caught by 3 tests |
| **PC5** | the multi-reference hold-out made silent again | ⛔ **NOT caught at first** |
| **PC6** | census narrowed to only what deposits | ✅ caught by 3 tests |

⛔ **PC5 is the C1-perturbation-6 pattern again**: `n_deposit_not_offered` was zero on both sides of the
identity in every fixture, so deleting it left the suite green. Unlike C0's L3 that is **not** a proven
no-op — it is a hole. Closed by `multi_reference_bam`, a hand-written two-contig BAM whose mates land on
different references and which reaches the adapter *by the intergenic path* — precisely the route the
shipped defect took. PC5 now fires.

### C2.1 — the anchor moved, and the fix is STRUCTURAL

`build_fl_models(payload)` — **the payload is the only argument**. All three histograms are read off one
object in one frame, so the mixed-frame call is **unrepresentable**, not merely discouraged. The EB kernel
over three free histograms survives as `_fl_models_from_histograms`, which production never calls.

⚠ **A value gate alone would have been vacuous.** On the plain oracle fixture the scanner's histogram and
`deposited_lengths` are **byte-identical** — a perfect BAM with no ambiguity makes definitions A/B and C
agree — so the test would have passed whichever anchor was wired in. A byte-identical result is no
evidence (`benchmark_ab_methodology_cautions`). The blacklist fixture separates them (200 vs 40).

⭐ **§2's table re-measured on all 8 pilot conditions — the gate's stated numbers, to the decimal:**

| `gdna_none ss0.50 capture_off` | anchor | RNA pool | gap |
|---|---|---|---|
| mean | **218.0** | 234.7 | **+7.7 %** |
| sd | **111.2** | 146.7 | **+32.0 %** |
| support ceiling | 1000 | 1000 | ⭐ matched |

C1's numbers exactly, against the pre-C1 **+11.6 % / +71.1 %**. The frame component is gone from the
*shipped* code, not merely available to it.

⚠ **NEW, and C3 should expect it: the tilt is capture-dependent.** `capture_on` reads **+11.0 % / +43.4 %**
against `capture_off`'s +7.7 % / +32.0 %. §8.1(b) is aimed at the junction-opportunity tilt; there is a
second effect riding on it.

⚠ **Three test harnesses had already drifted from production** — `_oracle.py`, `test_ambig_scenario`,
`test_accumulator_span_unbiased` still built FL models from the scanner's `category_models`, which
production stopped using at S5.d. They were calibrating against a model the tool does not ship. Converged.

### C2.2 / C2.3 — the deletion, and the report

Deleted: sites **A** and **B**, `FragLenObservations`, `frag_length_observations`,
`_replay_fraglen_observations`, `FragmentLengthModels` (**plural**), `n_frag_length_{un,}ambiguous`, and
`get_unique_frag_length_mrna` — **definition B itself**. `scan_and_buffer` returns a 4-tuple; 14 files
followed.

⚠ **`FragmentLengthModel` (SINGULAR) stays** and has its own survivor test. The two names differ by one
character and the singular is the scorer; a search-and-delete that caught it would have removed the
fragment-length term from scoring and reported it as numerical drift.

⭐ **The gate is a standing search over `src/`, not a one-off command** — `tests/test_one_fragment_length_definition.py`.
A partial delete that still compiles is the failure mode. ⚠ It found **three genuine stale references**
that a build could not: a nanobind docstring still advertising `'frag_length_observations'` as a returned
key, a comment in `resolve_context.h`, and a `config.py` docstring. It exempts an explicit
`DELETED by C2` tombstone and nothing else — a comment that merely *mentions* a deleted symbol is how the
next reader concludes it still exists.

**The report A/B against a worktree at `d045d820`, same scenario, same seed, `assignment_mode="map"`:**

| splice key | `d045d820` | **C2** |
|---|---|---|
| `unspliced` | 1482 | **1547** |
| `spliced_annotated` | 412 | 412 |
| `spliced_unannotated` | 0 | 0 |
| `spliced_implicit` | **8** | **41** |
| `splice_artifact` | 0 | 0 |
| **total** | **1902** | ⭐ **2000** |

⭐ **The key set is IDENTICAL — nothing was lost — and the old numbers were missing 4.9 % of the library.**
2000 is exactly the simulated fragment count. That reproduces G6's 4.6 % scanner-drop rate on a completely
different fixture, and it is the first time the report's splice section has summed to the library.
⚠ **`spliced_implicit` was under-counted 5×**: an implicit splice is exactly the case whose transcript-space
length most often fails the unanimity gate, so it was the worst-hit class.

`fragment_lengths.feather`'s categories change by design: the per-`SpliceType` histograms give way to
`global` + `gdna`/`rna` + the **five pure pools**. ⭐ That surfaces `splash_fl_mass` — the two crossing
pools, the only **on-target** gDNA population — for the first time; its own docstring asked to be reported
("makes that comparison an output instead of an assumption") and nothing was doing it.

### C2.4 / C2.5 — the dead field, and D7 CHECKED

`ScanCache.fl_global_counts`, `fl_rna_counts` (**D5** — written, read back, consumed by nothing) and
`fl_max_size` (a duplicate of `payload.max_length`) are gone, with `fl.npz` itself. ⚠ **Caches written
before C2 still load** — an extra file on disk is not a key — and all 8 pilot caches were verified loading
after the change.

⭐ **D7 is verified, not assumed** — the audit's own instruction. The shipped per-transcript effective
lengths are asserted **exactly equal** to those rebuilt from the payload alone, end to end through
`run_pipeline`. Two supporting gates: the quantity is *sensitive* to the pmf (shifting it 50 bp shortens
every effective length — an equality against a constant array would prove nothing), and perturbation
**PD1** (feeding `global_pmf` instead of `rna_pmf`) is caught.

⚠ **A finding, pinned rather than fixed**: `pipeline`'s guard `if rna_fl.n_observations > 0` **cannot
fire**. `from_pmf` sets `_total_weight = 1.0`, so `n_observations` is exactly 1 whatever went in. The
reachable behaviour is correct — an empty RNA pool EB-shrinks to the unconditional anchor, which beats a
hard-coded 200 bp mean — but the guard *looks* like a data-presence fallback, and that appearance is what
would let a future empty-pool bug hide behind it.

### ⚠ Unrelated: `scripts/sim/fl_estimation_stress.py` was ALREADY DEAD and is deleted

It called `calibrate(index=…, scan_trained=…, fl_prior_ess=…, pool_quality_good=…)` — **not one of those
is a parameter of the current `calibrate`**. Dead since well before this branch; the 2026-07-30 sweep
missed it. Recoverable: `git checkout d045d820 -- scripts/sim/fl_estimation_stress.py`.

### Suite

**1824 passed / 21 failed**, from a re-measured baseline of **1809 / 22**. Accounting is exact:
`1809 + 7 (C2.0) + 2 (C2.1 net) + 10 (C2.2 gate) + 2 (C2.3) + 3 (C2.5) − 9 (deleted container tests) = 1824`.

⭐ **The failure count went DOWN by one: `TODO.md` §7's `test_nrna_double_counting[g20_n0_s100]` now
PASSES.** The silent negative control reads **0 counts against a limit of 25**, where it leaked ~30 before
— not a marginal pass, and stable across three re-runs. ⚠ **Do not close §7 on this.** A negative control
is one-sided (trap 19) and this was not the change's target; the honest statement is that correcting the
fragment-length models removed the leak on this condition, and §7 should be re-read against the other
modes before it is retired.

⚠ **The 21 remaining are `test_golden_output` and they have now moved TWICE** — P1's units fix and C2's FL
models. ⛔ Still **do not regenerate**: they move again at C3. Regenerate **once**, after C3, **twice, and
diff**.

### ⚠ The standing check's INVOCATION matters, and it cost a false alarm

⛔ Use **`python -m pytest`**, not bare `pytest`. Bare `pytest` does not put the repo root on `sys.path`,
so `tests/calibration/test_fl.py`'s `import tests.native._accumulator_reference` raises and the suite reads
**1808 / 23**. That 23rd failure is an artefact of the invocation, not a regression — and under the rule
"a 23rd failure is a regression" it reads as one. `CLAUDE.md`'s command should say `python -m pytest`.

---

## C2.6 — an intron in the mate gap is cut on EVERY fragment, not only unspliced ones (2026-08-01)

    Spec: `docs/accumulator/SPEC_GAP_INTRONS.md`  ·  Cause and evidence: `docs/accumulator/JUNCTION_OPPORTUNITY.md` §4
    Owner ruling, 2026-08-01: *"we should be searching for gap introns within every fragment."*
    Baseline: e2f41cd0, re-measured in this session — **1824 passed / 21 failed**, all 21 goldens.

C0 proved the accumulator's `L` correct **given its inputs**. This entry is about an input that was
incomplete. `resolve_context.h` ran implicit-splice detection only on fragments the resolver had already
classified `SPLICE_UNSPLICED`, so a fragment carrying an observed CIGAR-N splice never had its
**unsequenced mate gap** examined, and an annotated intron sitting in that gap stayed inside `L`.

⚠ `UNSPLICED` never meant "one aligned block" — an unspliced paired-end fragment already has two blocks
and a mate gap, and that case always worked. The missed population is **spliced fragments that also have
a gap intron**, necessarily long because they span two or more introns. That is exactly the tail.

### The change, in three places

| | | |
|---|---|---|
| **§2.1** | `collect_implicit_splice_introns` takes the fragment's **observed** introns and drops every gap that matches one by **exact `(start, end)` equality** | ⛔ the trap: the gap finder walks consecutive aligned blocks, so a CIGAR-N intron is a "hole" too |
| **§2.2** | detection is **unconditional**; the `SPLICE_IMPLICIT` **promotion** stays unspliced-only | `splice_type` feeds scoring, the buffer, strand training and the report's census. This work is about `L`, not classification |
| **§2.3** | the deposit adapter **unions** the two intron lists (sorted, then de-duplicated), `sj_implicit = !implicit_introns.empty()`, `path_ambiguous = implicit_ambiguous` | the two lists stopped being mutually exclusive |

⭐ **The exact-equality rule is the whole of §1.** `transcript_has_implicit_intron_in_gap` accepts any
annotated intron inside the gap **± K** (K = 3 by default), so an overlap-based filter — or none — lets a
**different** nearby intron answer for one the CIGAR already stated. The two then normalise into one
wider interval and `L` comes out too **SHORT**. Measured under X2 on the hand-written fixture: the
near-match fragment reads **500 bp where the molecule is 502**.

### ⭐ THE HEADLINE — `gdna_none ss0.99 capture_off`, scored on `truth_fragment_lengths.tsv`

Truth: mean **217.13**, sd **87.41**, ceiling **713 bp**. Every target is read from that file or is a
control; **nothing here is tuned**.

| the anchor (`deposited_lengths`) | before | after | |
|---|---|---|---|
| **sd vs truth** | **+26.97 %** | ⭐ **+1.98 %** | **G-sd — 92.7 % of the excess removed** |
| mean vs truth | +0.38 % | −1.88 % | ⚠ see the selection effect below |
| mass **> 713 bp** (truth: **0**) | 0.00909 | ⭐ **0.00137** | **G-tail — 85 % removed, ⛔ not 0** |
| mass ≥ 700 bp | 0.00965 | 0.00143 | |
| `qc.dropped_too_long` | **280,558** (5.71 %) | ⭐ **38,309** (0.80 %) | 86 % collapse |
| `qc.deposited` | 4,636,640 | 4,770,233 | |
| `qc.dropped_ambiguous_path` | 82,802 | 191,458 | +108,656 — the new deferrals |
| `qc.sj_implicit_fragments` | 46,700 | 323,654 | 6.8× — the mixed population, seen for the first time |
| `qc.introns_absorbed` | 0 | **0** | the two lists really are disjoint |

### ⛔ G-gdna — THE CONTROL, AND IT DID NOT MOVE ONE DIGIT

`DNA_INTERGENIC` is pure gDNA and gDNA **has no introns to miss**, so a change there would mean the fix
reached fragments with no introns — impossible, therefore a bug. On `gdna100 ss0.99 capture_off`:

| | before | after |
|---|---|---|
| mean vs truth gdna | −0.5976 % | **−0.5976 %** |
| sd vs truth gdna | −0.3424 % | **−0.3424 %** |
| ≥600 / ≥700 bp | 0.00023 / 0.00001 | **0.00023 / 0.00001** |
| ceiling · n | 785 · 4,527,474 | **785 · 4,527,474** |

⭐ Bit-identical on every statistic **including the fragment count**. The control is clean.

### ⛔ G-tail DID NOT REACH 0 — and the residual is MEASURED, not excused

The gate says the mass above the true ceiling must be 0. It is **0.00137**. Per the spec's §6 the
remaining mass was measured rather than closed with a constant, by an experiment (**M3**) that emits
**every** annotated intron inside a gap instead of the first:

| | shipped | M3 (all introns in the gap) |
|---|---|---|
| mass > 713 bp | 0.00137 | ⭐ **0.00002** |
| `dropped_too_long` | 38,309 | ⭐ **389** |
| anchor sd vs truth | +1.98 % | −1.28 % |

⭐ **D3 is 98.5 % of the residual, and it is now a measurement rather than a hypothesis.**
`transcript_has_implicit_intron_in_gap` returns the **first** matching intron and stops, so a mate gap
spanning two annotated introns keeps only one cut. ⚠ M3 is a MEASUREMENT ONLY — it emits multiple
introns without extending the per-gap unanimity test to compare intron *sets*, which is the real work.
`TODO.md`-grade, and it is now the only known mechanism left in the tail.

### ⛔ D1 WAS MEASURED AND IT COSTS MORE THAN IT BUYS — the owner's call to revisit

D1 (recommended in the spec, implemented as ruled) removes a mixed fragment from `RNA_SPLICED`, because a
length partly inferred from the annotation is a product of the model that pool is used to fit. The spec
required the two effects to be reported **separately**; they were, by an experiment (**M1**) that cuts the
intron but leaves the fragment in the pool.

| `RNA_SPLICED` vs truth | before | ⭐ M1: intron cut, fragment KEPT | shipped: intron cut, fragment REMOVED |
|---|---|---|---|
| **mean** | +8.00 % | ⭐ **+0.67 %** | ⛔ **−9.58 %** |
| **sd** | +67.35 % | ⭐ **+2.40 %** | ⛔ **−22.46 %** |
| n | 1,475,626 | 1,609,219 | 1,332,265 |
| mass > 713 bp | 0.02812 | 0.00364 | 0.00000 |

⭐ **Cutting the intron very nearly fixes the pool (+8.00 → +0.67 mean, +67.35 → +2.40 sd). Removing the
mixed fragments then breaks it again in the opposite direction**, because the fragments removed are
exactly the ones whose mates sit far apart — a **length-selection bias**, and the pool is what the
fragment-length model is fitted from. ⚠ The purity argument for D1 is unchanged and still real; what is
new is its price, and it is 10 % of the pool's mean. Reversing it is one line
(`path.sj_implicit = implicit_only`), which is literally M1. **Owner decision.**

### ⚠ The anchor's mean went +0.38 % → −1.88 %, and about half of that is SELECTION

`dropped_ambiguous_path` rose by 108,656 (2.3 % of deposits) because a gap whose candidates disagree now
defers the fragment even when its *other* splice was observed. Those fragments are long by construction,
so the surviving anchor is biased short. Experiment **M2** lets them deposit again:

| anchor | shipped | M2 (ambiguous still deposits) |
|---|---|---|
| mean vs truth | −1.88 % | **−0.83 %** |
| sd vs truth | +1.98 % | +4.44 % |
| mass > 713 bp | 0.00137 | 0.00183 |

⭐ So ~1.05 pp of the −1.88 % is the deferral, not a measurement error — and deferring is the right call
(`L` genuinely is undetermined), so this is a cost recorded, not a defect.

### Gates

| | | |
|---|---|---|
| **U1** | a spliced fragment's mate-gap intron is found, and `L` excludes **both** cuts | ✅ resolver-level and end-to-end (`deposited_lengths` reads exactly `{300, 500, 502}`) |
| **U2** | the observed CIGAR-N intron is **not** re-derived | ✅ |
| **U3** | a near match within ±K is **not** substituted, and `L` does not shrink | ✅ |
| **U4** | `splice_type` does not move; the census is unchanged | ✅ 3 `SPLICED_ANNOT`, 1 `SPLICED_UNANNOT`, **0** implicit |
| **U5** | a MIXED fragment whose candidates disagree is rejected **and counted** | ✅ |
| **D1** | only the fully-observed fragment stays in `RNA_SPLICED` | ✅ two-sided |
| **S3** | byte-identity to `_accumulator_reference.py`, bit-identical at 1/2/4/8 workers | ✅ automatic, unchanged |

⚠ **`_accumulator_reference.py` did NOT move, and that is correct.** The specification describes the
**accumulator**, whose contract is `FragmentPath` in and a tally out. Which introns a fragment *has* is
decided upstream, in the resolver and the deposit adapter. `sj_implicit`'s meaning — "this `L` depends on
an intron that was never sequenced" — is unchanged; what changed is which fragments satisfy it. S3 never
reopened, and `payload_schema_digest` never moved.

### The perturbations

| | perturbation | result |
|---|---|---|
| **X1** | restore the `splice_type == SPLICE_UNSPLICED` gate | ⭐ **6 gates fail** |
| **X2** | drop the §2.1 observed-gap filter | ⭐ **7 fail**, including U2 and U3. `near` reads **500 bp against 502** — the "L too SHORT" defect, exactly as §1 predicts. ⚠ Also caught by the worker-determinism module's own non-vacuity guard: every spliced fragment re-derives its own observed intron, so **`pool_lengths` sums to zero** |
| **X3** | `sj_implicit` back on the splice class | ⭐ **D1 alone fails** — precisely the predicted gate, nothing else |
| **X4** | `path_ambiguous` gated on implicit-only | ⭐ **U5 fails** (+2 collateral) |
| **X5** | leave the union **unsorted** | ⛔ **NOTHING FAILS — the spec's prediction was wrong** |

### ⛔ X5: the spec predicted S3 byte-identity would break. It does not, and here is why

`normalise_introns` **sorts the intron list itself**, in the C++ accumulator and in
`_accumulator_reference.py` alike. So an unsorted union cannot change `L`, the segments, the junction ids
or any tally — the adapter's sort protects exactly one quantity, `qc.introns_absorbed`, by keeping the
adjacent-pair de-duplication a de-duplication.

⚠ And the case where that matters is close to unreachable: an implied intron can only duplicate an
observed one under K-slack, and it must then be **non-adjacent** in the concatenation, which needs a
fragment with ≥2 observed introns *and* a coincident K-slack match. The sort is kept because it costs
nothing on lists of ≤4 and restores the invariant the de-duplication was written against — but it is
**defensive, not load-bearing**, and this entry records that rather than claiming a gate it does not have.
`qc.introns_absorbed == 0` is now asserted on the fixture, which pins §2.1's disjointness claim instead.

### ⚠ Instrumentation added, and why it is not test-only scaffolding

`ResolvedFragment` now exposes `implicit_introns` and `implicit_ambiguous`. Without it **U2 is
unobservable**: the adapter de-duplicates and the accumulator normalises, so an intron re-derived from an
observed splice is merged away and "never emitted" is indistinguishable from "emitted and absorbed". §1
is precisely about the case where they differ, so the emission itself has to be visible. ⛔ It is
**copied, not moved**, out of `RawResolveResult` — `bam_scanner.cpp` builds the `ResolvedFragment` before
calling `deposit_to_accumulator(frag, cr)`, so a move would hand the adapter an empty list and silently
stop cutting every gap intron.

### ⚠ Cost: +8.9 % scan time

`gdna100 ss0.99 capture_off`, 10,000,000 fragments, `OMP_NUM_THREADS=1`, best of three:
**845.4 ns/fragment against 775.6** with detection gated back to unspliced-only. Detection is now a
binary search per candidate per gap on every fragment rather than on the unspliced ones.

### Blast radius, as predicted

`payload_schema_digest` **did not move** (no new field), so the 8 pilot caches stayed loadable — ⛔ but
their contents were wrong and all 8 were rebuilt (45.6 s). The goldens move a **third** time. ⛔ Still
regenerate **once**, after C3, **twice, and diff**.

### Suite

**1835 passed / 21 failed**, from a re-measured baseline of **1824 / 21** in the same session.
Accounting is exact: `1824 + 4 (U1–U4, resolver level) + 7 (end-to-end module) = 1835`. All 21 failures
are `test_golden_output`, unchanged in identity from the C2 entry.

### Files touched

`src/rigel/native/resolve_context.h` · `src/rigel/native/bam_scanner.cpp` · `src/rigel/native/resolve.cpp`
· `scripts/design/fl_anchor_gap.py` (the truth panel: G-tail, G-sd, G-gdna, D1, D3) ·
`tests/test_implicit_splice.py` · `tests/native/test_gap_introns_are_cut.py` (new)

### ⚠ A correction to the record, carried from `docs/accumulator/JUNCTION_OPPORTUNITY.md` §4.4

**C1 read the support ceiling backwards.** It recorded the anchor's ceiling moving from the scanner's
**713** to the accumulator's **1000** as "the support ceiling mismatch closed entirely". ⛔ **713 was the
library's true maximum**, to the base pair; 1000 is `max_frag_length`, the clamp. The deleted definition
**B** took its length from the *transcript*, so it could not produce an uncut intron — on this one
question it was right and `L` was wrong. ⚠ That does not undo C2, whose reasons all still hold.

---

## S1 + S2 — the accumulator ARBITRATES, and the side buffer exists (2026-08-01)

    Plan: `docs/accumulator/PLAN_TWO_PASS.md` §4 (S1, S2)   ·   Spec: `docs/accumulator/SPEC_GAP_PATHS.md` §0–§4
    Baseline: 750cc8ee — **44 failed / 1824 passed**, of which 21 are the stale goldens and
    **23 were a half-finished interface change**. After: **21 failed / 1863 passed**, all 21 goldens.
    ⛔ NOT judged by a calibration A/B, by design — see "what this deliberately does not measure".

`750cc8ee` moved the SPECIFICATION (`tests/native/_accumulator_reference.py`) to the hypothesis
interface and left the C++ behind it. This entry finishes that migration and lands the side buffer.

> A fragment arrives at the accumulator with its **hypothesis set**, not with one path and a flag. Exactly
> one survivor **deposits**; two or more and `L` is not one number, so the fragment is **held WHOLE** for
> the second pass. The empty hypothesis is the genomic one and needs no flag.

### What landed

| | | |
|---|---|---|
| **the C++ caught up** | `DeferredFragments` gains `ref`, is int64 throughout, and gains a canonical `sort` | the bank is the ONE tally whose order is observable; everything else is a sum of integers |
| **`Accumulator` knows which reference it IS** | `ref_id` is a required ctor argument, stamped into every held record | the drain replays a record through `deposit` onto **that** reference's cut axis |
| **the binding exposes both new banks** | `Accumulator.deferred` and `.gap_resolution`, through the same two exporters `build_result` uses | one set of key strings on this side of the ABI, so the parity surface and the payload cannot disagree |
| **the payload carries them** | `AccumulatorPayload.deferred` (a validated CSR) and `.gap_resolution` (a typed census) | |
| **S2: the cache round-trips them** | nested banks split by sub-field type: arrays to the `.npz`, counters to the manifest | ⛔ `dataclasses.asdict` on the bank yields ndarrays, and the manifest is written with `json.dumps(..., default=str)` — which would stringify each array to a **truncated repr**, silently |

### ⭐ THE OWNER'S KNOWN DEFECT IS FIXED, AND IT IS NOW GATED

`payload_schema_digest()` hashed `AccumulatorPayload`'s top-level field names **alone**, so a change
inside a nested bank was invisible to it — a renamed `ScanQC` field let a stale cache be accepted by the
key and then fail deep in the loader with a bare `TypeError`, exactly the failure the digest exists to
prevent. It recurses one level now. ⛔ S1 raised the stakes rather than lowering them: `DeferredFragments`
puts **thirteen** array names inside one payload field and every one is an `.npz` key.

⚠ **The fix had no gate until a perturbation said so** — see X8 below.

### ⭐ FINDING 1 — `GapResolution.RESOLVED_UNSPLICED` COULD NOT BE ENTERED, and it is deleted

The census had a class documented as *"the genomic hypothesis survived alone because every spliced path
was longer than `max_fragment_length`"*. **That is impossible.** A spliced hypothesis CUTS bases the
genomic one keeps, so `L_spliced <= L_genomic` always; the one arbitration filter is
`L <= max_fragment_length`; therefore if the genomic path survives the filter **every** spliced path does
too, and the survivor set can never be `{genomic}` while a spliced path was offered — which is the
condition for being in the census at all.

⭐ **Found by writing the non-vacuity gate**, not by review: the parity battery asserts every census class
is reached, and this one could not be. Checked before deleting — 200,000 random hypothesis sets, the
genomic path was the longest in **every one**.

⭐ **The ordering is now pinned as the REASON rather than the class as a consequence**
(`test_the_GENOMIC_hypothesis_is_ALWAYS_the_LONGEST`, 20,000 randomised sets). A future filter that broke
the ordering fails there instead of quietly needing a deleted class back.

### ⛔ FINDING 2 — THE REFERENCE STAMP WAS UNTESTABLE, AND TWO FIXTURES HAD TO MOVE

Two perturbations **passed the entire suite of 1860 tests**:

| | perturbation | why nothing caught it |
|---|---|---|
| **X2** | `deferred_.append(fragment, 0, …)` — a constant instead of the accumulator's own id | every accumulator-level fixture was reference 0 |
| **X3** | `AccumulatorSet` gives every reference the id `0` | every scan fixture was single-contig, or deferred only on reference 0 |

⛔ A constant is indistinguishable from a correct value until two references both hold something, and the
consequence is not cosmetic: the second pass replays each held record onto `ref`'s cut axis, so a constant
stamp drains one chromosome's coordinates onto another's partition — **the defect the predecessor deposit
adapter actually had**. Two fixtures changed rather than one assertion added:

* `test_accumulator_native_parity.py` now sits on **reference 3**, with three leading references that
  contribute no cuts (legal, and it exercises the per-reference offset arithmetic too);
* `test_gap_introns_are_cut.py` is now **two contigs**, with the disagreeing-isoform gene on chr2.

Both perturbations fire after the change. ⚠ The extent check is the half that matters most: `[66100,67000)`
is a plausible interval on chr1 as well, so "stamped 1" and "inside reference 1's cut range" are different
statements and only the second fails against the real partition.

### Behaviour that MOVED, deliberately, and where it is recorded

| | change | why |
|---|---|---|
| **D1 is deleted** | the RNA pool is keyed on **determinacy**, not provenance | a fragment enters a pool only when exactly ONE hypothesis survived, so its `L` is not in doubt however it was arrived at. C2.6 measured the alternative: **+0.67 % / +2.40 %** against truth on determinacy, **−9.58 % / −22.46 %** on provenance. ⭐ A purity filter on a length pool is a length filter |
| **D3 / C2.7 is solved** | every annotated intron in a gap is cut, and candidates group by their **whole-fragment path** | two transcripts differing only in their SECOND intron are now two hypotheses. ⛔ Its measurement cannot be re-read until the drain — see below |
| ⚠ **a SPLICED_UNANNOT fragment with a gap intron now DEFERS** | `docs/accumulator/SPEC_GAP_PATHS.md` §2: ∅ is available whenever no **annotated** junction was sequenced | the `near` fixture's junction is 2 bp off the annotation, so the molecule is not certified RNA and the genomic explanation survives. Its `L` gate moved from "deposits at 502" to "the enumeration produced the intron set that yields 502 and never the near-miss that yields 500" — asserted on **coordinates**, which is stronger |
| ⚠ **a transcript the read CONTRADICTS is no longer a candidate** | `docs/accumulator/SPEC_GAP_PATHS.md` §3, concern C1 | `test_h_any_intron_satisfies`'s old geometry ran a block contiguously across the transcript's FIRST intron, so the read has bases that transcript splices out. Cutting its *other* intron on that transcript's authority is a path no molecule has. The case keeps its meaning on a geometry that does not contradict, and the contradiction is now its own gate (`test_h2_…`) |
| ⚠ **the `max_fragment_length` prefilter changes CLASSIFICATION** | `docs/accumulator/SPEC_GAP_PATHS.md` §8, concern C3 | the same annotation defers or resolves depending only on how far apart the exons sit — a 300 bp intron holds the fragment, an 1800 bp one rules the genomic path out by length and deposits it. ⭐ Now gated **two-sided on a real scan** rather than assumed |

### Gates

| | | |
|---|---|---|
| **byte-identity** | C++ ↔ the specification, over every `Tally` field including the deferred CSR and the census, read off `dataclasses.fields(Tally)` | ✅ 54 named cases + 10,000 random fragments with random hypothesis sets |
| **non-vacuity** | the battery reaches **every** census class, and defers > 100 of the 10,000 | ✅ — this is what found FINDING 1 |
| **conservation** | `deposited + deferred + dropped_* == offered`, and the bank holds the fragments the counter claims | ✅ at the reference level, at the payload door, and on a real scan against the scanner's **independent** splice census |
| **determinism** | bit-identical at 1/2/4/8 workers, deferred bank included | ✅ at the accumulator (2/4/8 shards) and end-to-end through the scan |
| **S2 round trip** | the side buffer survives the cache byte-identically and comes back **typed** | ✅, and the fixture is asserted non-empty first |

### The perturbations

| | perturbation | result |
|---|---|---|
| **X1** | the canonical sort does nothing | ⭐ **5 fail** — parity and worker determinism |
| **X2** | the `ref` stamp is a constant | ⛔ **passed 1860 tests**; after the fixture moved, **7 fail** |
| **X3** | `AccumulatorSet` stamps every reference 0 | ⛔ **passed 1860 tests**; after the fixture moved, **1 fails** — the one that names it |
| **X4** | the census swaps `rna_or_gdna` and `which_introns` | ⭐ **6 fail** |
| **X5** | the fragment is COUNTED as deferred but not HELD | ⭐ **12 fail** — including the payload's door check |
| **X6** | the length filter keeps only the first hypothesis when the filter empties the set | ⭐ **1 fails** — precisely the case written for it |
| **X7** | `lex_compare` treats a PREFIX as equal | ⭐ **3 fail** |
| **X8** | `payload_schema_digest` stops recursing | ⛔ **not caught** — the fix had no gate; one added, then **1 fails** |
| **X9** | the cache writes a nested bank's arrays into the manifest JSON | ⭐ **2 fail** |

### ⛔ WHAT THIS DELIBERATELY DOES NOT MEASURE

**No calibration A/B, and that is the plan's instruction** (`docs/accumulator/PLAN_TWO_PASS.md` §4). Between S1 and S3 the
tally is deliberately **thinner**: the ambiguous mass is held, not yet deposited. And the held population
is the **long** one — a longer gap admits more hypotheses — so the surviving anchor is biased short.

⛔ **Every junction-opportunity number in the docs predates the drain and must be re-measured after S3,
not before.** So must D3's residual, which is a statistic of the deposited set.

### Files touched

`src/rigel/native/calibration/accumulator.{h,cpp}` · `src/rigel/native/bam_scanner.cpp` ·
`src/rigel/scan_payload.py` · `src/rigel/scan_cache.py` · `tests/native/_accumulator_reference.py` ·
`tests/native/test_accumulator_native_parity.py` · `tests/native/test_gap_hypothesis_arbitration.py` ·
`tests/native/test_gap_introns_are_cut.py` · `tests/native/test_fragment_length_proof.py` ·
`tests/native/test_implicit_splice_deposit.py` · `tests/native/test_accumulator_worker_determinism.py` ·
`tests/test_accumulator_payload.py` · `tests/test_scan_cache.py` · `tests/test_implicit_splice.py` ·
`tests/calibration/_synthetic.py` · `tests/calibration/test_gdna_strand_integration.py` ·
`scripts/design/fl_anchor_gap.py` · `scripts/design/implicit_splice_census.py`

### Blast radius

⛔ **`payload_schema_digest` MOVED** — it recurses now, and the payload gained two fields. Every scan cache
is invalidated **by design**; the 8 pilot caches rebuild (~46 s). The goldens move a **fourth** time.
⛔ Still regenerate **once**, at the very end, **twice, and diff**.

---

## S2.1 — the pipeline is REPRODUCIBLE, and it was never the seed (2026-08-01)

    Gate: `tests/test_scan_order_independence.py`   ·   Retires: `TODO.md` rank 3
    Owner question: *"the RNG seed issue seems minor — is there an easy solution?"*
    Baseline before: 21 failed / 1863 passed. After: **21 failed / 1868 passed**, all 21 goldens.

`TODO.md` rank 3 recorded that `EMConfig.seed` *"does not make `assignment_mode="sample"` reproducible"*
and filed it under the seed. ⛔ **The seed was never the problem**, and the recorded diagnosis sent the
investigation to the wrong subsystem.

### What was measured, before anything was changed

| | |
|---|---|
| a seeded `sample` run, twice at the default thread count | ⛔ counts differ by up to **22** |
| the `rng_seed` handed to the C++ EM on those two runs | ⭐ **bit-identical** — `9008465316566228445` |
| the EM's own thread count (`EMConfig.n_threads` 1 / 2 / 4 / auto) | ⚠ **irrelevant** — fails at all of them |
| scan `total_threads` = 1, 2, 4 | ✅ reproducible |
| scan `total_threads` = 8, 16 (**the default** = `os.cpu_count()`) | ⛔ **not** reproducible |
| the fragment buffer, two 16-thread scans | ⛔ **different row order** (1 chunk → identical; 6 chunks → not) |
| the `frag_id → fragment` mapping, same two scans | ⭐ **byte-identical**, ids contiguous `0..n-1` in BAM order |

⭐ **So: the scanner streams finalized chunks in worker-COMPLETION order, the fragment buffer inherits
that, and the EM's sampled assignment consumes its per-locus RNG in buffer-row order.** Permuting the
units permutes which fragment gets which draw, and units in different equivalence classes have different
posteriors, so the assignment genuinely moves. ⚠ The 1- and 4-thread passes were an artifact of scale —
one chunk, nothing to reorder. On a real library, which emits many chunks at any thread count above one,
the shipped default was **not reproducible at all**.

### The fix — an ORDERING, not a seed, in one place

```python
unit_indices=_in_bam_order(comp_u_indices[lo:hi], em_data.frag_ids)
```

`build_multi_loci` orders each locus's units by `frag_id`. Every per-locus array is scattered by that
index list (`locus_partition.partition_and_free`), so **all of them inherit the canonical order and no
consumer has to know**. ⚠ `unit_indices` is used nowhere else except for its length, so nothing else moves.

### ⛔ WHY NOT SEED EACH DRAW FROM THE FRAGMENT — the plausible wrong fix, and it is caught

Hashing the fragment's own content would also be order-independent, and would be **wrong**: identical
fragments hash alike and so draw alike, collapsing a 60/40 posterior to **100/0** for every group of
duplicates. That is the owner's own D-D ruling (`docs/accumulator/PLAN_TWO_PASS.md` §5.3) arrived at from the other side.

⭐ **And it fails even on its own terms.** A content key **ties exactly on the duplicates**, so the tie
falls back to buffer order — measured (perturbation Y2): ordering by `gdna_log_liks` instead of `frag_id`
leaves **both** `sample` and `fractional` thread-count dependent. `frag_id` is an identity, so it never
ties; ordering by identity and keeping one stream per locus is what preserves the multinomial spread.

### ⭐ THE CODE GOT SMALLER — 46 lines deleted

`em_solver.cpp`'s `build_equiv_classes` sorted the ROWS within each equivalence class on a
log-likelihood fingerprint, ~45 lines, explicitly *"because multi-threaded BAM scanning produces
fragments in non-deterministic order"*. Fixed at the source, that reason is gone: units now arrive
canonical, so `unit_list` is already canonical. **Deleted, and proven a no-op** — the golden
`effective_length` reads `501.82302350785244` before and after, to the last digit.

⚠ The **class** sort is KEPT. It guards a different mechanism — `unordered_map` iteration order, which is
a property of the hash table, not of the scan — and ⛔ it has **no gate**: perturbation Y3 removes it and
nothing fails, because a single build's hash order is stable. Pre-existing; recorded, not fixed.

### ⭐ A SECOND DEFECT FELL OUT, and it was never filed

`fractional` mode was thread-count dependent too — it scatters float posteriors into shared accumulators,
so a permutation reorders the summation and the answer drifts by ULPs. Nobody had noticed because nobody
had compared *across* thread counts. The same ordering fixes it, and the gate covers all three modes.

### Cost

One `argsort` per locus over an int64 key. Measured: **0.15 s at 1 M fragments, 0.91 s at 5 M** — against
a ~66 s calibration. ⭐ An earlier draft densified the ids into a global rank first; it is exactly
equivalent (rank is monotone in `frag_id`) and cost **5×** more — 731 ms of that 906 ms — so it was
deleted too.

### Gates

| | | |
|---|---|---|
| **the contract** | one BAM, one seed, one answer — byte-identical across scan thread counts {1, 4, 16} × 2 repeats, in **all three** assignment modes | ✅ |
| **non-vacuity** | the buffer's row order really does differ between those thread counts | ✅ retried 3×, because it is a race |
| **the assumption** | no two EM units share a `frag_id` — otherwise the sort ties and falls back to buffer order | ✅ (2,118 unique ids from 3,000 fragments; a subset, so gapped but unique) |

⚠ **Two repeats per thread count, not one.** The permutation is random: a single pair can agree by luck,
and did — two 16-thread runs disagreed while 1 vs 4 vs 16 matched on the same tree.

### The perturbations

| | perturbation | result |
|---|---|---|
| **Y1** | revert to buffer-row order | ⭐ **fires** |
| **Y2** | order by CONTENT (`gdna_log_liks`) instead of identity | ⭐ **fires on both `sample` and `fractional`** — the wrong fix is caught |
| **Y3** | delete the equiv-CLASS sort as well | ⛔ **not caught** — pre-existing gap, recorded above |

⚠ Y1 fired via `fractional` rather than `sample` on that run. That is the honest shape of this gate:
`sample` is the mode that matters and the flaky detector; `fractional` is the reliable one, because ULP
drift is deterministic given a permutation. Together they hold.

### What this changes for the owner

⭐ **`assignment_mode="sample"` is now a legitimate A/B arm.** `docs/SESSION_HANDOFF.md` and `TODO.md` §7 both say
to A/B under `map`/`fractional` because the default wobbles by ~0.5 %. With a fixed seed all three modes
now reproduce exactly. ⚠ The three remain **different estimators** — 5441 / 6002 / 6277 on one scenario —
so an A/B must still hold the mode fixed across both arms.

⭐ **And it answers D-D ahead of S3.** The second pass needs no new seed machinery: S1 already gives the
deferred queue a canonical order, so `(global_seed, index in that queue)` is well defined for exactly the
reason established here. The rule is *order by identity, then one stream* — not *hash the content*.

---

## P0 — there was never a strand sign bug (2026-08-02)

    Spec: `docs/accumulator/SPEC_SECOND_PASS.md` §3.3, P0   ·   Gate: `tests/calibration/test_strand_sense_convention.py`
    Retires: `TODO.md` rank 6   ·   Corrects: S5.f, S5.f-addendum, `docs/SESSION_HANDOFF.md` §1 fact 17 / §0 C4
    Suite: **21 failed / 1873 passed**, all 21 goldens.

`TODO.md` rank 6 recorded that the fitted κ is `1 − truth` and that *"only the exported scalar is
mis-labelled"*. ⛔ **Both statements are wrong, and the item had been filed twice.** The owner approved
fixing the sign; the audit found there is no sign to fix.

### The collision: two quantities, both called "strand specificity"

| | |
|---|---|
| `ReadSimConfig.strand_specificity` | *"probability an RNA fragment preserves correct read orientation … an R1↔R2 swap with probability 1 − ss"* — protocol **fidelity**, direction-agnostic |
| `StrandModel.p_r1_sense` | `P(align_strand == the junction's strand)` — **directional**. Its own docstring already said it: *"High (≈0.95) for R1-sense libraries (KAPA). Low (≈0.05) for R1-antisense (Illumina TruSeq dUTP)."* |

⭐ **For an R1-antisense protocol these are complements**, so comparing one against the other reads as a
sign error and is not one. The simulator emits R1-antisense — explicit in the gDNA path,
`r2_seqs, r1_seqs = _batch_extract_reads(...)`, the first extracted read becoming **R2** — which is the
most common real protocol. A dUTP library at 99 % fidelity genuinely has a sense fraction near 0.01.

### ⭐ THE MEASUREMENT — the tool already exposes the matching quantity

`StrandModel.strand_specificity = max(p_r1_sense, p_r1_antisense)` is direction-agnostic, like the
simulator's knob. Two genes, one per strand, 4,000 zero-gDNA fragments:

| simulated `ss` | `p_r1_sense` | `strand_specificity` | `rna_sense_frac` |
|---|---|---|---|
| 1.00 | 0.0000 | ⭐ **1.0000** | 0.0008 |
| 0.99 | 0.0156 | **0.9844** | 0.0164 |
| 0.75 | 0.2299 | **0.7701** | 0.2303 |
| 0.50 | 0.4980 | **0.5020** | 0.4980 |

⭐ It recovers the simulated parameter directly. The comparison the record should always have made.

⚠ **Opposite-strand genes are load-bearing in that fixture.** On a single-strand locus a convention error
that swapped the comparison's operands would move both together and be invisible.

### ⭐ And the 166× measurement now has an explanation rather than a rationalisation

S5.f-addendum forced κ to the nominal 0.99 and a zero-gDNA library read `f_gdna = 0.4992` against the
fitted value's 0.0030, and concluded *"the mirror cancels"*. ⛔ It does not cancel — **there is no mirror**.
`0.0101` is the right answer for an R1-antisense library, and 0.99 is a different quantity substituted for
it. The numbers stand; the explanation changes, and it is the simpler one.

### ⭐ What this unblocks

`docs/accumulator/SPEC_SECOND_PASS.md` §3.3's strand term, immediately and with no fix. `rna_sense_frac` is exactly the
`P(align_strand agrees | RNA)` the second pass needs to score an unspliced fragment's competing strand
hypotheses.

⚠ **But mind the direction.** On an R1-antisense library the value is ≈ 0.01, so the RNA hypothesis whose
implied strand *disagrees* with `align_strand` is the LIKELY one. ⛔ A scorer written as "agreement ⇒
multiply by `rna_sense_frac`" is exactly backwards on every real cfRNA library. Recorded in the spec.

### ⚠ A real gap found while closing a phantom one

**The simulator can only emit R1-ANTISENSE libraries.** `strand_specificity` is a swap probability about a
*fixed* orientation, never a choice of orientation, so **no simulated condition exercises the R1-sense
branch** — and real R1-sense libraries (KAPA-style) exist. Filed as `TODO.md` rank 6b and marked in place
by a **strict xfail** that deletes itself when the simulator gains the switch. ⚠ Low urgency: the branch is
a `max()` and a comparison, and real cfRNA is dUTP.

### Gates

| | | |
|---|---|---|
| **the recovery** | `StrandModel.strand_specificity` recovers the simulated knob at 1.00 / 0.75 / 0.50 | ✅ within 0.03 — ~4σ of the binomial standard error at 4,000 fragments, loose enough never to flake and far tighter than the 0.5-scale error a real flip would give |
| **the direction** | a perfectly stranded simulated library has `p_r1_sense ≈ 0`, and `read1_sense` reports the protocol as antisense | ✅ |
| **one convention** | `rna_sense_frac` is the Beta posterior mean of `p_r1_sense` and does not diverge from it | ✅ — a divergence would mean a second strand convention had appeared between the model and the balance fit |

⛔ **No code changed.** The deliverable is the gate and the corrected record: the absence of that gate is
what let one non-defect be filed twice, and it is the only thing that stops a third time.

---

## P1 — strand in the gap enumeration: two behaviours gated, one defect fixed (2026-08-02)

    Spec: `docs/accumulator/SPEC_SECOND_PASS.md` §3.3 and D-5   ·   Gate: `tests/native/test_gap_hypothesis_strand.py`
    Suite: **21 failed / 1877 passed / 1 xfailed**, all 21 goldens.

The owner's 2026-08-02 ruling: *"any fragment with a splice junction, immediately, we can constrain the
gap [hypotheses] to that strand … fragments without a splice junction are unspliced and could be either
strand."* The audit found the first half **already implemented and tested nowhere**; this entry gates both
halves and fixes the one case that was wrong.

### ⭐ Gated: an observed junction pins the strand

`tP` (+) and `tM` (−) cover the same sequenced blocks and share the observed junction's coordinates, so
neither annotation nor overlap can separate them — only the sequenced **motif** can. Measured through a
real scan:

| observed XS motif | hypotheses | outcome |
|---|---|---|
| `+` | tP's 400 bp gap intron alone | deposits at **L = 450** |
| `−` | tM's 200 bp gap intron alone | deposits at **L = 650** |

⭐ The deposited `L` says which transcript was believed, with no need to read the hypothesis set. ⛔ And it
fires on the count too: without the pin both would survive and the fragment would be **deferred instead of
deposited**.

### ✅ Gated: an unspliced fragment offers both strands

Same locus, no motif: three hypotheses — `(1600,2000)` on `+`, `(1700,1900)` on `−`, and ∅ — with the two
spliced ones carrying **opposite** implied strands. ⚠ This is the case that makes the second pass's strand
term necessary rather than decorative: drop it and these two are separated by length alone.

### ⛔ Fixed: D-5 — one path claimed by both strands

`add_hypothesis` grouped paths by intron **coordinates only**, so two opposite-strand transcripts implying
the same intron merged into one hypothesis carrying the **first-seen** strand — an answer that flips when
two GTF lines are swapped. Grouping is right (it *is* one path); the strand is now set to **AMBIGUOUS**
when supporters disagree, which is what that value already means everywhere else and which `deposit`
already refuses to credit a junction on. Idempotent, three lines.

⚠ Unreachable on human data — **0 of 404,168** junction coordinates are annotated on both strands, and the
index warns that it is biologically impossible. Fixed because the alternative is order-dependent.

### ⚠ THE FIXTURE TRAP, recorded because it cost an hour

Driving `FragmentResolver` directly leaves `t_strand_arr_` empty, so **every hypothesis's implied strand
silently reads NONE** and a strand gate written that way passes against nothing. ⛔ All of these gates go
through the scan.

### The perturbations

| | perturbation | result |
|---|---|---|
| **Z1** | the merged hypothesis keeps the first supporter's strand | ⭐ **fires** — D-5 |
| **Z2** | the implied strand is dropped; every hypothesis reports NONE | ⭐ **fires on two gates** |
| **Z3** | `!certified_rna` removed from the ∅ condition | ⛔ **not caught** — and the reason is worth having |

### ⛔ Z3 DOES NOT FIRE, AND THAT IS NOT A WEAK FIXTURE

Measured: the `open` fragment's candidates are `t_inds = [0, 1, 4, 5]`, where **4 and 5 are the nRNA
shadows** (`RIGEL_NRNA_chr1_{1,2}_1000_2200`, `is_synthetic`). A shadow is single-exon, so it implies
nothing in the gap, so `any_candidate_implies_nothing` is already true and ∅ is emitted without the clause.

⭐ **That is `docs/accumulator/SPEC_GAP_PATHS.md` §2's own ruling arriving by a different route** — *"the nascent shadow IS
the ∅ hypothesis"* — and the two mechanisms agree rather than conflict:

* **spliced** fragment: the observed CIGAR-N intron falls inside the shadow's single exon, so the shadow
  cannot explain the read and drops out of `t_inds`. ∅ is correctly absent — the certified-RNA `mixed`
  fragment in `test_gap_introns_are_cut.py` has exactly one hypothesis.
* **unspliced** fragment: the shadow survives and supplies ∅, which is what §2's table requires.

⚠ So `!certified_rna` is **redundant, not wrong**, and it is KEPT: it states the rule directly rather than
depending on the shadow mechanism continuing to exist, and a single-exon gene has no separate shadow row
at all (`is_nrna` on a non-synthetic row means mature ≡ nascent). ⛔ Do not delete it on the strength of
this perturbation.

### Files touched

`src/rigel/native/resolve_context.h` (the D-5 fix) · `tests/native/test_gap_hypothesis_strand.py` (new)

---

## P2 — the scorer is BUILT and UNGATED (2026-08-02) ⛔ PARTIAL, READ THIS BEFORE USING IT

    Spec: `docs/accumulator/SPEC_SECOND_PASS.md` §3 and P2   ·   Suite: **21 failed / 1877 passed / 1 xfailed**
    ⛔ **`src/rigel/second_pass.py` HAS NO TESTS.** It is inert — nothing imports it — but it is
    production code in `src/` that no gate covers, which is exactly the state where the next reader
    assumes it works. It does not yet have the right to be believed.

### What landed

| | |
|---|---|
| `Accumulator.length_under` | ⭐ `L` under ONE hypothesis without depositing. **Exposed, not reimplemented**: a Python scorer computing its own `L` would be a second definition of the quantity C0/C2 unified, and the one the drain would then disagree with. Verified: spliced 101, genomic 800, clipping correct |
| `MarshalledFragment` | the binding's hypothesis marshalling factored out, shared by `deposit` and `length_under`, so the two cannot drift about how a hypothesis set crosses the ABI |
| `src/rigel/second_pass.py` | the scorer: all three terms of §3, each kept separately on `HypothesisTerms` so a regression is attributable to one rather than to "the score moved" |

⚠ **D-1 is a parameter (`min` vs `geometric`) BECAUSE it must be measured.** It is not a tunable to be
left in: the docstring says it is expected to be deleted once the measurement is in. ⛔ If it is still
there in a month, that is the failure mode `docs/SESSION_HANDOFF.md` warns about.

### ⭐ THE SMOKE RUN ALREADY SHOWS D-3, AND IT IS NOT SUBTLE

Four-fragment fixture, both held records scored **uniform**:

```
record 0: [62100,63000)   L=502  rho=2.00e-03  f=0.000   -> undecided
record 1: [66100,67000)   L=500  rho=0.000     f=0.500   -> undecided
```

On record 1 **both junctions have zero flux** — the only fragment that would have used them is the held
one. That is exactly the systematic hole the spec predicted: *the held fragments are excluded from the
tally they are scored against*. On a toy it is 100 %. ⛔ **The real number is unknown and it is the fact
that decides D-3.**

### ⛔ WHAT IS DELIBERATELY NOT DONE, AND WHY

**The P2 gate is not written.** It says *"on a hand-built locus where the truth is known, the correct
hypothesis takes the larger share"* — and the available fixture has no depth, so every score is uniform
and the gate would pass or fail for reasons that have nothing to do with the scorer. ⚠ Writing it against
this fixture would have produced a green gate over an undecided scorer, which is worse than no gate.

Three things remain, and the order matters:

1. ⭐ **Measure D-3 on the pilot first** — pure measurement, data already on disk
   (`~/Downloads/rigel_runs/suite/pilot`, 8 conditions with `truth_fragment_lengths.tsv`). It needs no new
   fixture, and its answer may change what the gate should even assert. ⛔ If the zero-flux share is large
   the score needs a fallback, and **that is an owner decision — no pseudocount will be invented**.
2. **A discriminating fixture, then the gate.** Deep enough that rho and f actually separate.
3. **Measure D-1**, then delete the parameter.

---

## P2.1 — ⭐ D-3 MEASURED, and the ∅ term is measuring an artefact (2026-08-02)

    Decision: `docs/accumulator/SPEC_SECOND_PASS.md` §8 D-3   ·   Phase: its P2, step 1 of 3
    Tool: `scripts/design/held_flux_census.py` (new)   ·   Data: the 8 pilot conditions, caches
    rebuilt (44.5 s) because S1 moved `payload_schema_digest`.
    Suite: **21 failed / 1877 passed / 1 xfailed**, all 21 goldens — unchanged. ⛔ **No `src/` change.**

D-3 asked whether a junction with zero flux is *impossible* or *merely unobserved*, and warned that the
held fragments are excluded from the tally they are scored against. ⭐ **Both halves now have numbers, and
they point opposite ways from what the smoke run suggested.**

### ⭐ FINDING 1 — the score vector almost never empties; the zeros DECIDE instead

| | zero-gDNA | gdna100 |
|---|---|---|
| records where every hypothesis scored 0 (`n_undecided`) | 0.063 – 0.297 % | 0.122 – 1.010 % |
| ⭐ records where a zero **eliminates** a hypothesis while another lives | **67.3 – 79.8 %** | 55.7 – 68.9 % |

⛔ **The failure D-3 names is a non-event and the one it does not name is everywhere.** The four-fragment
smoke run read 100 % undecided; at scale it is ≤ 1 %. A pseudocount would move `decisive`, not
`undecided` — so it would not be repairing the failure it was proposed for.

### ⛔ FINDING 2 — the ∅ term's evidence set is WRONG, and on gdna100 it is the ENTIRE zero rate

`second_pass._lines_inside` asks for the contiguous lines **strictly between** an intron's endpoints. The
deposit rule (`_accumulator_reference.py:818-824` — a line is crossed iff strictly inside a *segment*)
says the lines that distinguish ∅ from a path splicing `[a, b)` are those at cuts **`a <= c <= b`**,
endpoints included. Brute-forced against the reference, not argued:

```
    intron (200,400):  ∅ credits [0,1,2,3,4]   spliced credits [0,4]     ⇒ distinguishing [1,2,3]
                       scorer asks for [2]                                  ⛔ misses [1,3]
    intron (200,300):  ∅ credits [0,1,2,3,4]   spliced credits [0,3,4]   ⇒ distinguishing [1,2]
                       scorer asks for []                                   ⛔ EMPTY — rho is 0 structurally
```

Both endpoints are cuts whenever the intron is annotated, so the two guaranteed discriminating lines are
exactly the two dropped, and the set is **empty whenever the intron spans one node**. Measured, ∅'s
zero-rho rate under the shipped rule against the rule derived from the deposit:

| condition | ∅ shipped | ⭐ ∅ derived | artefact |
|---|---|---|---|
| gdna100, capture off | 0.4338 / 0.4344 | ⭐ **0.0000** | **+0.434** |
| gdna100, capture on | 0.5276 / 0.5535 | 0.0486 / 0.0728 | +0.479 / +0.481 |
| zero-gDNA, capture off | 0.8219 / 0.7931 | 0.7217 / 0.6897 | +0.100 / +0.103 |
| zero-gDNA, capture on | 0.8798 / 0.8711 | 0.7399 / 0.7378 | +0.140 / +0.133 |

⭐ **On gdna100-capture-off the derived rate is EXACTLY 0.0000 over ~90,000 ∅ hypotheses**, against 43 %
shipped. The whole ∅ zero rate on a gDNA-bearing library is manufactured by the evidence set.

⭐ **And once it is removed the number tracks composition, which is the falsification working.** 0.00
where gDNA exists and ~0.72 where it does not: on a zero-gDNA library nothing crosses an intron
contiguously, so ∅ *is* impossible and a zero density is the correct evidence saying so. ⛔ The residual
0.72 is **signal, not a hole to be filled**.

### ⭐ FINDING 3 — D-3's real population is the spliced arm, and 69–89 % of its zeros are CORRECT

| | |
|---|---|
| spliced hypotheses with `rho == 0` | **17.9 – 22.8 %** (capture on is the high end) |
| of those, cause `unannotated` (`jid < 0`) | ⭐ **0, in all 8 conditions** — every implied intron resolves to a slot |
| of those, cause `annotated_empty` | **100.0 %** — the slot exists and carries no flux |

Scored against the simulator's `truth_abundances.tsv` (`observed_mrna_fragments`, so no target is
chosen): of the 817–880 annotated-empty slots a held hypothesis claims, **695–745 have no expressed
supporter at all**. Claims landing on a slot whose supporter truth says was expressed: **11.2 – 31.4 %**,
and that is an **upper bound** — the criterion is the max over a hypothesis's supporters, transcript-level
rather than junction-level.

### ⛔ FINDING 4 — the self-exclusion mechanism D-3 predicted is NOT what is happening

The spec's reasoning was that a junction used *only* by deferred fragments reads zero. If that dominated,
the heavily-claimed empty slots would be the expressed ones. They are the opposite:

| | slots | claims | median claims/slot | max |
|---|---|---|---|---|
| supporter expressed | 122 | 13,627 | 20 | 925 |
| ⭐ all supporters silent | 695 | 73,794 | **24–29** | **1,509** |

⭐ **The single heaviest slot has 1,509 held claimants and truth says ZERO molecules used it.** A high
claim count is evidence that the annotation offers the path often, not that anything took it. ⛔ **That is
the measured argument against a pseudocount**: it would hand those 1,509 fragments a live route to a
junction the library never used.

### The perturbations

| | perturbation | result |
|---|---|---|
| **W1** | the census stops calling a zero-flux slot zero | ⭐ **fires** — the self-check names 66,359 hypotheses. The decomposition is verified against `score_held_fragments`'s own density array, so a census that drifts from the code aborts |
| **W2** | the derived ∅ rule made identical to the shipped one | ⭐ **fires** — the derived column collapses to exactly the shipped 0.4344, so the arm measures the rule difference and nothing else |

### ⛔ WHAT IS DELIBERATELY NOT DONE

**`_lines_inside` is NOT fixed here.** It is a `src/` change to an ungated module, and P2's own order puts
the discriminating fixture and the gate next — a fix landed before the gate is a fix nothing can hold.
⚠ It is also not a free-standing bug: `∅ derived` on zero-gDNA is 0.72, so correcting the evidence set
does **not** by itself make ∅ selectable, and the two must be judged together.

⛔ **D-3 remains the owner's call.** The measurement says a pseudocount is not indicated: the empty-vector
case is ≤ 1 %, 69–89 % of spliced zeros are truthfully zero, and the heavy zeros are the correct ones.

### Files touched

`scripts/design/held_flux_census.py` (new)

---

## P2.2 — ⭐ THE P2 GATE EXISTS, D-6 is fixed, and D-1 is CLOSED AND DELETED (2026-08-02)

    Spec: `docs/accumulator/SPEC_SECOND_PASS.md` §3 and its P2 (step 2 and step 3)
    Gate: `tests/test_second_pass_scoring.py` (new, 10 tests)   ·   Owner rulings, 2026-08-02, in §0
    Suite: **21 failed / 1887 passed / 1 xfailed** — was 21 / 1877 / 1. All 21 are still goldens.
    ⭐ **`src/rigel/second_pass.py` is no longer untested production code.** That was P2's whole debt.

### §0 The owner's two rulings this entry executes

| | |
|---|---|
| **D-3** | ⭐ **No fallback.** The measurement (P2.1) does not indicate one, and a hard zero stays hard |
| **D-6** | **Fixture and gate FIRST**, then land the fix and watch the gate move — not the reverse |

### ⭐ THE FIXTURE — seven loci, and every one of them is a MIRROR of something

`docs/WIP.md` P2 recorded why the gate was not written earlier: the smoke fixture had no depth, so every
score came out uniform and the gate would have been *green over a scorer that decided nothing*. The
answer is not a deeper fixture but a **paired** one.

| arm | locus | what it asserts |
|---|---|---|
| 1 | the **wide** junction is deep | the deep junction takes the larger share |
| 2 | ⭐ the mirror — the **narrow** one is deep | **the winner FLIPS** |
| 3 | the gap is deeply crossed contiguously | ∅ wins. ⛔ **this arm is D-6** |
| 4 | only the gap's **donor** end is crossed | ∅ LOSES — it needs evidence at both ends |
| 6 | ⭐ the mirror — only the **acceptor** end | same |
| 5 | two junctions at **equal** depth, different implied length | the length term decides |
| 7 | two **opposite-strand** hypotheses, equal width and depth | ⭐ scored twice on ONE payload at `rna_sense_frac` 0.99 and 0.01 — **the winner flips** |

⭐ **NO THRESHOLD APPEARS ANYWHERE IN THIS MODULE**, and the mirrors are why. Arms 1 and 2 share a
geometry, so every implied `L` — and therefore `f` and `s` — is *identical* between them; a test asserts
that rather than trusting it. A flipped winner then isolates `rho` **exactly**, where
`score(A) > 3 * score(B)` would have been a constant chosen after the fact. Arm 7 applies the same trick
to strand without a second fixture, by rescoring one payload at two library sense fractions.

Two further separations of concern, both of which a failure taught:

* ⭐ **Ballast sets the length pmf; the loci set `rho`.** `build_fl_models` does not smooth the global
  anchor, so `global_pmf[L]` is exactly 0 unless some deposited fragment had that length — a zero there
  would kill ∅ for a reason that has nothing to do with density. A gate checks the support exists.
* ⛔ **Arm 5's two lengths (300, 500) are reachable by NO locus.** It first used 200 and 500, and adding
  arm 7 later moved mass into the 200 bin and **broke it**. The fix was not to re-tune counts — it was to
  put both bins out of the loci's reach so the 3:1 ratio is structural.

### ⛔ D-6 FIXED — and the gate was verified failing first

`_lines_inside` → **`_distinguishing_lines`**, asking for cuts `a <= c <= b` instead of `a < c < b`.
Two lines changed. Written against intent, arm 3 failed exactly as predicted —
`genomic=0.0000 spliced=1.0000` — and passed after the fix. Re-measured on the pilot:

| ∅ zero-rho | before | ⭐ after |
|---|---|---|
| gdna100, capture off | 0.4344 | **0.0000** |
| zero-gDNA, capture off | 0.7931 | 0.6897 |

⚠ The zero-gDNA residual is **correct and stays**: with no gDNA nothing crosses an intron contiguously,
so ∅ genuinely has no support and the density says so.

### ✅ D-1 CLOSED, AND THE PARAMETER IS DELETED

Measured over the 8 conditions, `min` against `geometric` on one payload scored twice:

| | |
|---|---|
| records where the **winner** differs | **0.47 – 0.59 %** |
| records where a share moves > 0.10 | 1.9 – 3.4 % (max Δ 0.83) |
| ⭐ the **zero masks** | **bit-identical in all 8** — both return 0 on any zero input, so D-1 never interacted with D-3 |

⭐ **No accuracy argument separates them at that scale, and P2 has no truth criterion that could** — so
the derivation decides. `min` is the bottleneck (a molecule on the path was present at every object on
it); `geometric` had only a robustness intuition. `_aggregate(values, how)` is now `_bottleneck(values)`,
the keyword is gone from `score_held_fragments`, and the A/B script was **deleted with it** rather than
left behind to keep a dead flag alive. ⚠ Shares do move on 2–3 % of records; if P4's drained-tail gate
ever misses against truth, this is the documented place to look.

### The perturbations — 9 of 11 fire

| | perturbation | result |
|---|---|---|
| **Y1** | D-6 reverted to the strictly-between rule | ⭐ fires — arm 3 |
| **Y2** | only the **donor** line kept (half of D-6) | ⛔ **passed arms 1–3**; arm 4 was written for it, then **fires** |
| **Y2b** | only the **acceptor** line kept | ⛔ **passed arms 1–5**; arm 6 is its mirror, then **fires** |
| **Y3** | `rho` ignored | ⭐ fires on 4 |
| **Y4** | the **length** term ignored | ⛔ **passed arms 1–4**; arm 5 was written for it, then **fires** |
| **Y5** | the **strand** term ignored | ⛔ **passed arms 1–6**; arm 7 was written for it, then **fires** |
| **Y5b** | the strand term **INVERTED** | ⭐ fires — this is exactly `docs/WIP.md` P0's warning, now gated |
| **Y6** | a junction reads `count` instead of `inv_length_sum` | ⭐ fires — §1.1's rule, gated |
| **Y7** | §3.2's pmfs swapped (spliced→anchor, ∅→`rna_pmf`) | ⛔ **does not fire** — see below |
| **Y8** | `min` → `geometric` | ⛔ does not fire — ⭐ that is D-1's answer, not a hole |
| **Y9** | ∅ scored over its whole path, not the contested union | ⛔ **does not fire** — see below |

⭐ **FOUR arms exist only because a perturbation passed.** Arms 4, 5, 6 and 7 were each written after a
green battery said the gate was blind to something. That is the discipline working as intended, and it
took the module from 6 tests to 10.

### ⛔ THE TWO THAT DO NOT FIRE, recorded rather than papered over

* **Y7 — the choice of pmf (§3.2) is UNGATED.** The choice only bites where `rna_pmf` and `global_pmf`
  differ materially at the competing lengths *and* `rho` is near-tied. Arm 5 ties `rho`, but both of its
  hypotheses are spliced, so a swap moves them together. Closing it needs a `rho`-tied **∅-versus-spliced**
  locus, which means tuning two different object types to an equal density — a tuned fixture, which is
  what this module is built to avoid. ⚠ Left open deliberately; §3.2 is a live decision.
* **Y9 — scoring ∅ over the contested union rather than its whole path is UNGATED.** Same reason: at
  every locus here the rest of ∅'s path carries the same evidence as the contested part.

### Files touched

`src/rigel/second_pass.py` · `tests/test_second_pass_scoring.py` (new) ·
`scripts/design/held_flux_census.py` (follows the rename; keeps the pre-D-6 rule locally so P2.1's
`artefact` column stays reproducible)

---

## P3 — the DRAIN: one tally path, and it runs off a cached payload (2026-08-02)

    Spec: `docs/accumulator/SPEC_SECOND_PASS.md` §5 (the draw), §6 (the drain), D-2   ·   Phase: its P3
    Gates: `tests/native/test_accumulator_drain.py` (new, 16) · `tests/test_scan_cache.py` (+3)
    Suite: **21 failed / 1906 passed / 1 xfailed** — was 21 / 1887 / 1. All 21 are still goldens.
    Owner rulings, 2026-08-02: D-2 = the pure-function shape; §5.2 = one stream in queue order.

### ⭐ D-2 DECIDED AGAINST THE SPEC'S OWN RECOMMENDATION, and the reason is the scan cache

§6.3 recommended **(a)** — keep the scanner alive, re-enter the live `AccumulatorSet`. Three findings
moved it:

| | |
|---|---|
| ⛔ **(a) forces the drain INSIDE the scan** | `scan_and_buffer` builds the payload from `build_result()` and returns; the `AccumulatorSet` dies with the local |
| ⛔ **which breaks the scan cache** | `build_scan_cache.py` caches the payload and calibration reloads in 0.35 s instead of 61 s. A drain inside the scan bakes one draw into the cache, so every P5 ("would a second iteration move anything?") and P6 ("re-measure everything") becomes a full rescan of 8 conditions |
| ⭐ **(b)'s stated cost does not exist** | §6.3 objected that (b) "needs an `AccumulatorSet` binding with a payload export". It does not: **every tally channel is already a `def_prop_ro` on the bound per-reference `Accumulator`** — that is the parity surface — and `second_pass.py` already built accumulators from a payload for `length_under` |

**The shape that landed**: `drain(payload, choices, …) -> AccumulatorPayload`, pure, with **zero new C++**.
One `Accumulator` per reference rebuilt from the payload's own cut axis, the chosen hypothesis replayed
through the same bound `deposit`, and the delta added to pass one's arrays.

⭐ Three properties fall out rather than being engineered: the drain runs off a **cached** scan; the delta
is a **separate object**, so its contribution to every channel is a subtraction (§6.3's stated reason to
prefer (b)); and it **never sees a thread**, so reproducibility is structural.

### ⭐ §6.1's claim is CHECKABLE, so it is checked

> *"There is no second deposit implementation … byte-identity with the specification is preserved for free."*

`test_draining_a_choice_equals_depositing_it_directly` asserts exactly that: a fragment drained with
hypothesis `h` gives a byte-identically equal tally, over every `Tally` field, to the same fragment offered
`hypotheses=(h,)` in the first place. ⚠ Not a tautology — one route goes through arbitration with three
hypotheses, a canonical sort and a replay; the other deposits once.

### ⛔ TWO FINDINGS THAT WERE NOT ON THE LIST

**1. The census would have depended on the RNG.** `_record_gap_resolution` sends a size-one *spliced* set
to `RESOLVED_SPLICED`, and returns early on an all-unspliced one. So a naive drain grows
`gap_resolved_spliced` by however many draws happened to pick a spliced path, while chosen-∅ fragments
vanish from the census entirely. ⭐ **And there is no `gap_resolved_unspliced` class to put them in** —
S1 deleted it after checking over 200,000 random hypothesis sets that pass-one arbitration can never
produce one, because the genomic path is always the longest. **The drain, however, *chooses*.** So the
drain snapshots the census and restores it, and the information it could not have expressed lives on the
drain's own axis as `chose_genomic` / `chose_spliced`. ⭐ S1's deletion is confirmed correct by a second
route, and it is *why* the drain needs its own block rather than an extra census class.

**2. §6.2's bookkeeping conflicts with a payload door check, and the door wins.** §6.2 says
`deferred_undetermined_gap` is *kept* at pass one's count. But the payload refuses a bank that disagrees
with its counter — *"the counter and the fragments it counts must be the same population"* — and the drain
empties the bank. ⭐ Resolved the other way round: **after the drain the payload describes the FINAL
state** (bank empty, held counters 0, so both door invariants stay absolute and unconditional) and pass
one's numbers are preserved in `DrainQC.offered` and `DrainQC.census_before`. The identity §6.2 wanted is
intact — `deposited + dropped_* == offered` — it just lives entirely inside the drain block, where it is
self-contained. ⚠ The spec's letter changed; its statement did not.

### ⛔ AND TWO LATENT CACHE HOLES THE NEW FIELD EXPOSED

| | |
|---|---|
| **the schema digest recursed exactly ONE level** | enough for `ScanQC` / `GapCensus` / `DeferredFragments`, and not enough the moment `DrainQC` nested a `GapCensus` **inside itself** — X8's invisibility, one level down. `_schema_names` is fully recursive now |
| ⛔ **`drain` is `DrainQC \| None`, and `is_dataclass` is False for a union** | so a plain check dropped the **whole** bank from the key, not merely its nesting. Measured before the fix: `drain__` subfields were `[]` |

⭐ **And the cache now REFUSES a drained payload outright**, rather than teaching the serialiser to nest
two deep. The cache holds a *scan*; caching a drained payload would bake one draw into it and would also
push `census_before` through `json.dumps(default=str)` as a stringified repr — X9's defect one level down.
One check removes the whole question.

### The perturbations — 18 of 19 fire

| | perturbation | result |
|---|---|---|
| **Q1** | the census is not restored after the drain | ⭐ fires |
| **Q2 / Q3** | the held counter not zeroed / the bank not consumed | ⭐ 4 fail each |
| **Q4** | the drain always deposits the FIRST hypothesis | ⭐ 4 fail |
| **Q5** | the three `deferred_*` subclasses not zeroed | ⭐ fires |
| **Q6** | the draw drops the per-run offset | ⭐ fires |
| **Q7** | the draw drops the clamp into the run | ⛔ **did not fire** — see below |
| **Q8** | `ADDITIVE_AXES` forgets `pool_lengths` | ⭐ fires |
| **R1** | `set_junctions` never called | ⭐ fires — the invisible failure the C++ refuses on the scan path |
| **R2 / R3** | local junction offsets not rebased / acceptor cuts not localised | ⭐ 4 fail / fires |
| **R4** | the delta placed at offset 0 for every reference | ⭐ fires |
| **R5** | the library-wide axes skipped | ⭐ fires |
| **R6 / R7 / R8** | `with_drain` keeps the held counter / leaves the bank / forgets the deposits | ⭐ fires each |
| **S1 / S2 / S3** | digest stops looking through `Optional` / recurses one level / the cache accepts a drained payload | ⭐ fires each |

⛔ **Q7 did not fire, and the fix was to the GATE.** The clamp guards a draw that overshoots its run's
cumulative total, which with normalised scores needs a uniform within one ULP — random draws will never
produce it, so the original gate passed with or without the guard. Rewritten to **construct** the case:
scores summing to 0.5 per run make half of all draws overshoot. ⭐ The property gated is now the
function's contract — a choice is always inside its own record's run, whatever the caller's normalisation —
rather than one arithmetic accident.

⚠ **Two-contig fixtures throughout.** S1 found that a constant `ref` stamp passed all 1860 tests because
every fixture was single-contig; reference 0's junction slot base is 0, so every slicing error there is
invisible. R2/R3/R4 only fire because chr2's base is not.

### ⭐ MEASURED ON THE PILOT AT FULL SCALE — 8 conditions, all invariants holding

| | ∅ chosen | of held |
|---|---|---|
| zero-gDNA, capture off / on | 1,757 – 1,959 | **1.0 – 1.7 %** |
| gdna100, capture off / on | 2,135 – **4,804** | 1.2 – 4.4 % |

⭐ **The composition signal comes through**: ∅ is chosen 2.7× more often on a gDNA-bearing library than on
a zero-gDNA one at capture-on (4,804 vs 1,759), which is what the density term is supposed to do.
⚠ **But ~1 % of ∅ choices on a zero-gDNA library are false by construction** — there is no gDNA there, so
every one of them is wrong. Recorded, not fixed: P4's tail gate and P6's re-measurement are where that is
judged, and it is far too small to read as a defect in the score. `dropped_too_long` is **0–2** per
condition, so essentially every held fragment now deposits. Cost: **4–6 s** per condition.

### Blast radius

⛔ **`payload_schema_digest` MOVED again** (`cebaa14c91085c02` → `cc86144468ef210d`): the payload gained
`drain` and the digest recurses fully. Every scan cache is invalidated **by design**; the 8 pilot caches
rebuild in **61 s**. ⚠ The goldens do **not** move — nothing is wired into the pipeline yet, which is P4.

### Files touched

`src/rigel/second_pass.py` · `src/rigel/scan_payload.py` (`DrainQC`, `ADDITIVE_AXES`, `with_drain`) ·
`src/rigel/scan_cache.py` (full recursion, `Optional`, the refusal) ·
`tests/native/_accumulator_reference.py` (`Tally.deferred_canonical`, `Accumulator.drain`) ·
`tests/native/test_accumulator_drain.py` (new) · `tests/test_scan_cache.py`

---

## P4 — the second pass is WIRED IN, and ⭐ THE ANCHOR IS NOW EXACT AGAINST TRUTH (2026-08-02)

    Spec: `docs/accumulator/SPEC_SECOND_PASS.md` §2 (where the drain sits), its P4
    Gates: `tests/test_second_pass_pipeline.py` (new, 5) · `scripts/design/fl_anchor_gap.py --drain`
    Suite: **21 failed / 1911 passed / 1 xfailed** — was 21 / 1906 / 1. All 21 are still goldens.
    Caches rebuilt (61 s) before any number here was quoted.

`_drain_side_buffer` runs in `run_pipeline` between `scan_and_buffer` and `build_fl_models`, so calibration
sees the complete tally and runs **once** — §2's structural claim, now the shipped path. The seed is
`PipelineConfig.second_pass_seed`, ⭐ **deliberately not `em.seed`**: they are two independent RNG consumers
and sharing one field would mean an EM A/B silently re-drew every held fragment, moving the tally the two
arms were being compared on.

### ⭐⭐ THE HEADLINE — the anchor's bias against truth closes to ZERO

On the zero-gDNA conditions every fragment is RNA, so anchor-vs-truth is a pure bias measurement and there
is nothing to argue about:

| `deposited_lengths` vs truth | before the drain | ⭐ after |
|---|---|---|
| mean, capture off | −1.61 % | **+0.00 %** |
| mean, capture on | −0.97 % / −0.98 % | **+0.00 % / −0.00 %** |
| sd, capture off | −1.48 % / −1.47 % | **+0.02 %** |
| sd, capture on | −0.92 % / −0.93 % | **+0.01 %** |

⭐ **This is what the whole two-pass structure was built to do.** The held population is the LONG one — a
longer gap admits more hypotheses — so holding it out biased the anchor short; returning it makes the
anchor *exact* against the simulator's own record of what it wrote. `docs/accumulator/PLAN_TWO_PASS.md` §2.4 said B's target
could not be read until the buffer drained. It has drained, and the anchor needs no correction at all.

### ⭐ THE gDNA CONTROL IS EXACT — all four pools move by ZERO, and there is a derivation

| | |
|---|---|
| `DNA_INTERGENIC`, `DNA_INTRONIC`, `DNA_INTRON_EXON`, `DNA_INTERGENIC_EXON` | ⭐ **delta 0**, all 8 conditions |
| `RNA_SPLICED` | +103,183 … +171,449 — **100 % of the delta** |

Not a tolerance, and not luck. A held fragment's gap contains ≥ 1 annotated intron whose endpoints are
cuts, so a chosen ∅ crosses **both** those lines — a multi-line crossing, which `docs/accumulator/DESIGN.md` §8
deliberately gives **no pool** because it is a gDNA/RNA mixture — and a chosen spliced path used an
annotated junction, so it is `RNA_SPLICED`. A drained fragment therefore *cannot* enter a pool that is
supposed to be pure, whichever hypothesis won.

⚠ **This is NOT C2.6's gDNA control and must not be quoted as it.** There the fix could not reach a
fragment with no introns, so any movement was a bug. Here the drain reaches fragments it is *supposed* to:
on gdna100 a real gDNA fragment whose mate gap spans an annotated intron is genuinely ambiguous and was
genuinely held. What is gated is that such a fragment never lands in a pure pool.

### ⛔ THE TAIL GATE DOES NOT HOLD AS WRITTEN — 0–3 fragments, and the cause is fully diagnosed

§9's P4 asks for *"no fragment above the library's true longest molecule"*. Measured as an absolute count,
not a rounded fraction:

| | before | after | added |
|---|---|---|---|
| zero-gDNA (true ceiling 713 / 729) | **0** | 1 – 3 | **1 – 3** |
| gdna100 ss 0.99 (785) | 0 | **0** | **0** |
| gdna100 ss 0.50 (843) | 0 | 1 – 2 | 1 – 2 |

⭐ **100 % of them come from an UNDECIDED record** — 3 of 3 and 6 of 6 on the two conditions traced
fragment by fragment, and **zero** from a record the evidence decided. An undecided record's score vector
is uniform (the owner's D-3 ruling: no fallback), so ∅ can win with an `L` longer than any real molecule
and only `max_fragment_length` stops it. ⚠ That filter *is* doing its share: of 3 over-ceiling choices on
one condition the `L = 1194` one was rejected `TOO_LONG`, and of 6 on another the 2047 / 3274 / 5875 ones
were — which is exactly why the counts above are smaller than the number of bad choices.

⭐ **And the undecided population collapsed when D-6 landed**: 154 → **10** records on
`gdna_none_ss_0.99_capture_off` (0.089 % → 0.006 %), a **15×** reduction. ∅ having a real density is what
did it. So the residual is 2 fragments in ~5.0 M deposited — ~4 × 10⁻⁷ — and it traces entirely to a ruled-on
design choice rather than to the drain.

⛔ **Owner call, not a fix to invent.** The clean non-magic option is to exclude, from an *undecided* record's
uniform draw, any hypothesis whose `L` has no support in the anchor — which is the same statement
`max_fragment_length` already makes ("not a molecule this library contains") rather than a new constant. But
it is a rule change, so it is recorded, not taken.

### ⚠ AND ONE SUBSTANTIVE QUESTION P4 SURFACED — the RNA pool grows, and it is S4's number

A drained fragment whose chosen intron resolves to an annotated junction enters `RNA_SPLICED`, even though
its splice was **implied and never sequenced**. Measured against the unconditional mRNA truth:

| `RNA_SPLICED` vs truth | before | after |
|---|---|---|
| mean | +2.36 … +4.55 % | **+6.23 … +8.12 %** |
| sd | −6.40 … −4.16 % | −3.22 … −1.57 % |

⭐ **This is junction opportunity re-measured on the drained tally, which is exactly what S4 was waiting
for.** `docs/accumulator/PLAN_TWO_PASS.md` §2.4's table (+2.4 % raw, −3.2 % with the true θ) is **superseded**: the hold-out
was masking the bias, and the correction has a bigger gap to close than it appeared.

⚠ **The pooling itself is defensible and was already ruled on.** Excluding drained fragments would re-create
exactly the failure D1 was deleted for — they are the LONG ones, and *a purity filter on a length pool is a
length filter* (C2.6 measured provenance-keying at −9.58 % / −22.46 %). The circularity is real but
**one-step**: pass one's RNA pmf informs a choice whose result then feeds calibration's RNA pmf, and §7.1
forbids closing that loop. ⭐ **Whether one step is enough is precisely what P5 measures**, and this table
is its input.

### The perturbations — all 6 fire

| | perturbation | result |
|---|---|---|
| **T1** | the drain is not wired in at all | ⭐ fires |
| **T2** | calibration's FL models built from PASS ONE while the drain still runs | ⭐ fires — the gate proves *which* payload calibration used, not merely that a drain happened |
| **T3** | the drain reads `em.seed` instead of its own | ⭐ fires |
| **T4** | the seed never reaches the draw | ⭐ fires |
| **T5** | a drained ∅ is credited to the `DNA_INTRONIC` pool | ⭐ fires — the pure-pool control |
| **T6** | the empty-side-buffer path drains anyway | ⭐ fires |

⚠ **T2's first version did not fire and was rewritten.** It "moved the drain later" by draining *into*
`build_fl_models`, which is semantically the correct code — a perturbation has to change behaviour, not
just shape.

### ⚠ One test had to move, and it is the wiring proving itself

`test_d7_transcript_eff_lengths` rebuilt the payload from an independent scan and compared it to
`run_pipeline`'s. Production now drains, so the two diverged by **~2.7 bp on every transcript**. Its claim —
*"every eff length comes from the payload's RNA pmf"* — is unchanged; "the payload" now means the drained
one. It calls `pipeline._drain_side_buffer` rather than repeating its three steps, so "the same route
production takes" cannot quietly stop being true.

### ⛔ Not done, deliberately

**No CLI flag.** `second_pass_seed` defaults to 0, so `rigel quant` is reproducible with no new surface;
`--seed` reaching only the EM is the decoupling T3 gates. **The goldens still are not regenerated** — 21,
unchanged in count, and S5 is the very end.

### Files touched

`src/rigel/pipeline.py` (`_drain_side_buffer`, the wiring) · `src/rigel/config.py`
(`PipelineConfig.second_pass_seed`) · `scripts/design/fl_anchor_gap.py` (`--drain`: the before/after tail
and gDNA panels) · `tests/test_second_pass_pipeline.py` (new) · `tests/test_d7_transcript_eff_lengths.py`

---

## P4.1 — ⛔ THE ANNIHILATION BUG: a zero in one factor destroyed the other two (2026-08-03)

    Spec: `docs/accumulator/SPEC_SECOND_PASS.md` §3   ·   Found by: scoring the pilot against the simulator's
    PER-FRAGMENT truth, which the oracle BAM's read name carries verbatim
    (`ENST00000458178.2:1281-1331:f:281` — source transcript, then the molecule's own coordinates).
    Gate: `tests/native/test_accumulator_drain.py` (+6)   ·   Suite: **21 failed / 1917 passed / 1 xfailed**

The owner asked whether the RNA length pool growing to +6–8 % was a bug. ⭐ **It is not** — but looking for
it found a different one, and fixing that closed P4's tail residual exactly.

### ⭐ FIRST, THE ACCURACY, because it frames everything

171,534 held fragments on `gdna_none_ss_0.99_capture_off` with an unambiguous true length:

| | |
|---|---|
| exactly the right length | ⭐ **90.5 %** |
| mean error | **+0.1 bp** |
| the true answer was **among the offered candidates** | ⭐ **100.0 %** |
| by candidate count | 93.9 % (2), 87.9 % (3), 84.3 % (4), 79.3 % (5) |

⭐ **Pass one never misses the right explanation** — every failure is a *selection* failure, never a missing
candidate. That is a strong statement about the enumeration, and it is measured rather than assumed.

### ⭐ THE RNA POOL IS NOT A BUG — the fragments really are longer

| | |
|---|---|
| pass-1 `RNA_SPLICED` pool mean | 222.3 bp |
| **true** mean of the held fragments | **315.1 bp** |
| the length we **assigned** them | 315.2 bp |
| our error | ⭐ **+0.13 bp** |

The held fragments are genuinely **+92.9 bp** longer than the pool they join and we measure them to a tenth
of a base pair. The pool moved because its **composition** changed, and mechanistically it must: a fragment
is held when its gap admits several introns, which needs a wide gap; and long fragments are independently
more likely to span a junction, which is what the pool is selected on. ⛔ **So the +6–8 % is junction
opportunity at full size, not error** — S4's target, finally measurable.

### ⛔ AND "IS IT MIXING UP DNA AND RNA?" — no, and the framing needed correcting

On a zero-gDNA condition there is no DNA at all, so a chosen genomic path does **not** mean "this is DNA".
It means *"this RNA molecule was not spliced in this window"* — retained intron, unspliced region — and
measured, **64.3 % of those picks are correct**. The second pass picks a **path**; gDNA-vs-RNA is decided
downstream by calibration. ⚠ The 35.7 % that are wrong are *all* too long (mean +33.9, 0.0 % too short).

### ⛔ THE BUG — an all-zero factor annihilated the informative ones

`score = rho x f x s`, normalised within the candidate set. If **one factor was zero for every candidate**
the product was zero everywhere, the normalisation could not run, and the record fell back to a uniform coin
toss — **discarding the other two factors, which usually did decide.**

Measured on the 10 records that fell back to uniform:

| what decides | records it can decide | correct |
|---|---|---|
| traffic (`rho`) alone | 4 | 50 % |
| ⭐ **length (`f`) alone** | **8** | ⭐ **100 %** |

A worked case — record 118, truth 471 bp:

```
    L= 471   rho=0   f=0.0001    <- CORRECT; the length term knew
    L= 965   rho=0   f=0.0000    <- what the coin picked: longer than any molecule in the library
```

⭐ **The fix introduces no constant.** The score is normalised, so a factor taking the *same* value for
every candidate cancels and cannot affect the answer — `test_a_factor_that_is_CONSTANT_and_positive_already_cancels`
checks that premise. Zero is the one value where the arithmetic loses it: instead of cancelling it destroys
the product. So an all-zero factor is treated as what it is — **uninformative** — and dropped.

⛔ **The PARTIAL-zero case is untouched, and that is the point.** A factor zero for *some* candidates is
highly informative — the zero says "no evidence for this path" — and stays decisive. The owner's D-3 ruling
is left exactly as it was; nothing is floored, softened or smoothed.

### ⭐ Result — and it closes Decision 1 without a cutoff

| | before | ⭐ after |
|---|---|---|
| records reduced to a coin toss | 10 | **2** |
| of the 6 traced by hand, now correct | 0 (chance) | ⭐ **5 of 6**, several at ≈ 1.000 |
| fragments above the library's TRUE ceiling, zero-gDNA | 8 | **0** |
| … gdna100 | 3 | **0** |

⭐ **§9's P4 tail gate now holds LITERALLY AND EXACTLY on all 8 conditions**, and the anchor stays exact
against truth (+0.00 % mean on zero-gDNA).

### ⭐ THE ONE SURVIVOR SHOWED WHAT DECISION 1 ACTUALLY IS

Record 155262 stayed a coin toss: traffic favoured L = 739 and 1024, the length term favoured the true
L = 352, and each zeroed the other's preference — an **irreducible contradiction**, correctly identified.
But the coin then picked 739, longer than any molecule in the library.

So the owner-approved rule landed in its natural form: ⭐ **a coin toss may not offer a length the library
does not contain.** `f = 0` already means *"no fragment of this length was observed anywhere"* — the same
statement `max_fragment_length` makes, read off the **measured distribution** instead of a round number — so
the fallback draws only among candidates the length term has not ruled out, and narrows nothing when it
rules out everything. ⛔ **No empirical-max cutoff was needed and no constant was added**; the owner's
prediction that a cutoff "will never make a difference in real data" is now also true of simulated data,
because the ceiling is not the mechanism.

### The remaining two failure classes, both real and neither an arithmetic bug

| | |
|---|---|
| ⛔ **the correct junction had ZERO traffic while a wrong path had plenty** | record 121: truth is one 20 kb intron whose junction was never observed in the confident set (`rho = 0`); a wrong three-intron path had `rho = 1.6` and won with confidence **1.0000**. ⚠ So a hard zero does not merely make the right answer unselectable — it makes a wrong one **certain**. That is D-3, now measurable against truth |
| ⚠ **the length model was short-biased against the truth** | record 25886: truth is an unspliced 559 bp molecule, scored implausible because the pass-1 anchor excludes exactly the long fragments being held. §3.2 predicted this self-defeating case in advance |

### The perturbations

Six new gates, each written before the fix and verified failing: the all-zero factor, the partial zero
staying decisive, a constant positive factor already cancelling, the irreducible contradiction, all three
factors absent, and the restricted coin toss. ⚠ **`test_CONTRADICTORY_factors_are_still_undecided` had to be
retargeted**: its original construction used the *length* term as one of the two contradicting factors, and
the narrowing then resolved it — so an irreducible contradiction must be built from traffic and strand,
with the length term allowing both. Constructing it is what proves the narrowing has a boundary.

### Files touched

`src/rigel/second_pass.py` (`combine_factors`, the one place the three factors meet) ·
`tests/native/test_accumulator_drain.py`

---

## P4.2 — ⭐ THE COMBINATION RULE, corrected twice, and the FL subpopulations MEASURED (2026-08-03)

    Owner challenge: *"A fragment length of 739 should be exceedingly unlikely under the first-pass RNA FL
    distribution … how could traffic possibly overcome that? Aren't we using likelihoods?"*
    Suite: **21 failed / 1918 passed / 1 xfailed**, all 21 goldens. Gates: +3 in the drain module.

### ⛔ THE OWNER WAS RIGHT, AND P4.1'S FIX WAS HALF A FIX

Checked first: the RNA pmf's support ends at **713 bp**, exactly the library's true longest molecule, so
`f(739)` is **exactly 0** and the length term does annihilate that candidate — `score = 0.0000`. ✅ The
concern about traffic overcoming an impossible length does not apply to the score itself.

⛔ **But it applied to the FALLBACK.** P4.1 restricted the coin toss to lengths the library contains and then
flipped a **fair** coin among them. On record 155262 that left two survivors:

```
    L= 542   f=9.92e-06    score=0.5000   <== PICKED
    L= 352   f=1.42e-03    score=0.5000       TRUTH
```

⭐ **A 143-fold difference in likelihood, discarded.** Narrowing a draw and *weighting* it are different
things, and only the second is using the likelihood. That is the same mistake as the annihilation bug, one
level down, and it was mine.

### ⭐ THE RULE THAT REPLACES BOTH PATCHES — one sentence, no special cases

> Apply the factors **in order of how much evidence stands behind them**. Skip any factor that is flat zero
> across the candidates **still standing**; otherwise multiply by it and drop the candidates it zeroes.

The order is `f` → `s` → `rho`, and it is not arbitrary: an `f` of 0 rests on the whole library's length
distribution (millions of fragments), an `s` of 0 on its entire spliced population, a `rho` of 0 on however
many fragments happened to touch **one object** — often none. ⭐ *"Uninformative"* is now judged among the
survivors rather than globally, which is exactly what weights 352 over 542 by its true 143:1.

⭐⭐ **AND IT ELIMINATES THE "IRREDUCIBLE CONTRADICTION" CASE ENTIRELY.** Each factor either narrows a
non-empty set or is skipped, so the product can never collapse to zero everywhere. A weaker factor can no
longer veto what a stronger one left standing, and two factors can no longer annihilate each other.

| | before P4.1 | after P4.1 | ⭐ after P4.2 |
|---|---|---|---|
| records with tied leaders | 10 | 2 | ⭐ **0** |
| of the 6 traced by hand, correct | 0 (chance) | 5 of 6 | ⭐ **6 of 6** |
| fragments above the library's TRUE ceiling | 8 | 0 | **0** |
| exact length overall | 90.5 % | 90.5 % | 90.5 % (+0.12 bp) |

⚠ Records 121 and 25886 remain wrong, and neither is arithmetic: one has zero traffic on the *correct*
junction, the other a truth the short-biased pass-1 anchor calls implausible.

### ⭐⭐ THE FL SUBPOPULATIONS, MEASURED AGAINST TRUTH — and the answer is GEOMETRIC

Owner hypothesis: *"the unambiguous gap-intron containing fragments give us the true RNA FL
distribution."* Tested on 5,000,000 zero-gDNA fragments (all RNA, no composition confound), each category
scored against the same ground truth (mean **217.1**, sd **87.4**, max **713**, mass ≥ 500 bp **0.1404 %**):

| subpopulation | n | share | mean | Δmean | Δsd | max | ≥500 bp |
|---|---|---|---|---|---|---|---|
| gapless (`L` is a fact) | 2,244,818 | 44.9 % | 139.4 | **−35.80 %** | −53.65 % | **200** | **0.0000 %** |
| gap, RESOLVED | 2,582,789 | 51.7 % | 278.2 | **+28.11 %** | −32.48 % | 713 | 0.2266 % |
| gap, HELD | 172,393 | 3.4 % | 315.0 | +45.08 % | −27.96 % | 661 | 0.6781 % |
| observed splice | 1,750,715 | 35.0 % | 225.4 | ⭐ **+3.80 %** | −4.19 % | 676 | 0.1493 % |

⛔ **NEITHER HALF IS AN UNBIASED ESTIMATOR, AND THE REASON IS THE READ LENGTH, NOT STATISTICS.** With
2 × 100 bp reads, a fragment is gapless **iff it is ≤ 200 bp** — note the `max` of exactly 200 and the
**zero** mass above 500. So `gapless` is hard-censored above 200 and `gap` is censored below it: they are
complementary halves split at twice the read length, and each is badly biased on its own.

⭐ **And the owner's proposal is already what the tool uses.** gapless + gap-resolved = every fragment that
is *not* held = the pass-1 anchor exactly: **n = 4,827,607 (96.6 %), mean −1.61 %, sd −1.47 %.** Those are
P4's pre-drain numbers to the digit. So the proposal is sound and already implemented; its residual bias
**is** the hold-out, and the drain is what closes it — to **+0.00 %**.

⚠ **`observed splice` is the best single subpopulation** (+3.80 % / −4.19 %) and it is the only one that
spans the full range, because observing a CIGAR `N` is roughly length-independent while having a mate gap is
purely a length threshold. ⭐ Its +3.8 % **is** junction opportunity, measured directly for the first time.

### ⭐ WHAT THIS SETTLES ABOUT P5 (iteration)

Nothing is left for a second iteration to improve **on the anchor**: after one drain it is exact against
truth (+0.00 % mean, +0.02 % sd). The `RNA_SPLICED` pool stays ~+6 % long, but this table shows that is a
**selection** effect — the annotated-junction-using population genuinely is longer than the library — and
iteration cannot fix a selection effect. Only S4's opportunity correction can. ⛔ So the +6 % is not
circularity feeding on itself, and P5's answer looks like *"no iteration; go to S4"*.

⚠ **Read length and library width govern how much of this matters**, which is the owner's own caveat: a
tight, short FL distribution leaves more fragments gapless and fewer held, so the anchor is less biased to
begin with and the drain matters less. The conclusion is robust in the safe direction.

### ⚠ THE STRAND TERM IS THE REMAINING MODELLING GAP — the owner's second point, and it is right

Today `∅` gets a flat `s = 1.0` while a spliced candidate gets `rna_sense_frac` or its complement. ⛔ Those
are not likelihoods of the same observable, and the flat 1.0 systematically favours `∅`. The owner's
proposal is the correct shape: on a 99 %-stranded library a fragment aligned **antisense** to the gene is
very unlikely to be RNA, so `∅` (which may be gDNA) should be *favoured*; aligned **sense**, it is
uninformative between the two.

⭐ And it needs no calibration output, by exact analogy with §3.2's treatment of length: `∅`'s component is
unknown, so marginalise over the library's own composition — and `StrandModels.exonic` is that marginal,
measured in pass 1 (it is diluted toward 0.5 by gDNA on a contaminated library and sits at ≈ 0.01 on a
zero-gDNA one). ⚠ Not implemented: it changes the third factor on every score and is the owner's call.

### Files touched

`src/rigel/second_pass.py` (`combine_factors`) · `tests/native/test_accumulator_drain.py` ·
`tests/test_second_pass_scoring.py` (arm 7 ties by construction at a neutral strand fraction, so exactly
one tied record is expected — asserted with its reason rather than as zero)

---

## B4 — ⛔ THE DELIVERABLE, SCORED: the second pass IMPROVED the anchor and made COMPOSITION WORSE (2026-08-03)

    Tool: `scripts/design/calibration_truth_ab.py` (new)   ·   Replaces the S5.f composition baseline
    ⭐ **The measurement the whole fragment-length track existed to justify**, and it does not say what
    was hoped. One thing varied: the same cached scan, calibrated with the side buffer drained or not.

### The numbers, against the simulator's own origin counts

| condition | truth | undrained | err | drained | err | move |
|---|---|---|---|---|---|---|
| gdna100 ss0.50 capture off | 0.5000 | 0.5050 | +0.0050 | 0.4890 | −0.0110 | +0.0060 |
| gdna100 ss0.50 capture **on** | 0.5000 | ⛔ 0.3774 | **−0.1226** | 0.3704 | −0.1296 | +0.0070 |
| gdna100 ss0.99 capture off | 0.5000 | 0.4772 | −0.0228 | 0.4625 | −0.0375 | +0.0148 |
| gdna100 ss0.99 capture on | 0.5000 | 0.4981 | −0.0019 | 0.4894 | −0.0106 | +0.0087 |
| zero-gDNA (4 conditions) | 0.0000 | 0.0018 – 0.0854 | | ~unchanged | | ≤ 0.0006 |

⛔ **Mean |error| on the four contaminated conditions: 0.0381 → 0.0472, +23.9 %.** Every one moves the same
way — *down*, toward less gDNA — by 0.006 to 0.015.

### ⭐ THE MECHANISM, DECOMPOSED — and 55 % of it is the drain being RIGHT

| condition | deposited | held | held/dep | f predicted, MASS ONLY | f observed | residual |
|---|---|---|---|---|---|---|
| gdna100 ss0.50 capoff | 9,825,133 | 172,541 | 1.756 % | 0.4963 | 0.4890 | **−0.0073** |
| gdna100 ss0.50 capon | 9,879,931 | 108,359 | 1.097 % | 0.3733 | 0.3704 | −0.0029 |
| gdna100 ss0.99 capoff | 9,824,905 | 172,771 | 1.759 % | 0.4690 | 0.4625 | **−0.0065** |
| gdna100 ss0.99 capon | 9,880,420 | 108,044 | 1.094 % | 0.4927 | 0.4894 | −0.0033 |

⭐ **The MASS half is correct and unavoidable.** Truth says **99.8 %** of held fragments are RNA (306 gDNA
of 172,079, measured per-fragment on gdna100), so depositing them genuinely lowers the gDNA *fraction*.
⛔ And because the pre-existing estimate was already **too LOW** on 3 of 4 conditions, adding real RNA
necessarily moves it further from truth. **The drain is not creating this error; it is exposing one.**

⚠ **The RESIDUAL half is the fragment-length models changing** — −0.003 to −0.007, consistently negative.
`docs/accumulator/DESIGN.md` §7.2 puts a 10 % length-model error at 0.010–0.026 of composition; the RNA pool
moved ~4 pp long (P4.2: **+2.4 % → +6.2 %** against truth), which predicts ~0.004–0.010. Consistent.

### ⭐⭐ WHAT THIS SETTLES ABOUT THE ROADMAP — calibration consumes the POOLS, not the anchor

`calibrate` takes `gdna_fl_pmf` and `rna_fl_pmf`, which are the two **pure pools** EB-shrunk toward the
anchor. The anchor is now exact (+0.00 % / +0.02 %, P4). ⛔ **The RNA pool is not, and the second pass makes
it worse**, because the fragments it deposits are long, junction-using ones entering a pool selected on
"used an annotated junction" — which P4.2 measured as genuinely **+3.8 %** longer than the library.

⭐ So the fragment-length track fixed the quantity calibration *shrinks toward* and left the quantity it
*fits from* uncorrected. **C3 / S4 (junction opportunity) is therefore not polish — it is the step that
converts the second pass's gain into a composition gain**, and until it lands the second pass is feeding a
better-measured tally into a biased length model.

### ⚠ Two other things this baseline exposes, neither caused by this work

| | |
|---|---|
| ⛔ `gdna100 ss0.50 capture_on` reads **0.3704 against 0.5** | The worst row in the suite and the one S5.f already recorded as unexplained (0.3754 then). Unmoved by everything since. ⭐ It is the sharpest single target the suite has, and `length_likelihood_ab.py`'s docstring already nominates it |
| ⚠ `gdna_none ss0.50 capture_off` reads **0.0854 against 0** | An **unstranded** library: the strand channel carries no information at all, so composition is inferred from length alone and over-calls gDNA by 8.5 pp. The stranded zero-gDNA rows read 0.0018–0.0032 |

⚠ **One caveat on the comparison, stated because it flatters the undrained arm.** Both arms are scored
against the *library's* 0.5, but the undrained tally genuinely contains ~1.8 % less RNA, so its own tally
truth is ≈ 0.5088. Against that the undrained error is −0.0316 rather than −0.0228. The drained arm is the
one whose tally represents the library, so scoring the deliverable against 0.5 is right — but the raw
before/after gap is smaller than the table suggests.

### Files touched

`scripts/design/calibration_truth_ab.py` (new)
