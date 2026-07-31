# S5 — the design log

**S5 was "rewire the consumers". It is not.** The accumulator was rebuilt to make new things possible,
and what those things are had never been derived — so `IMPLEMENTATION_PLAN.md` §5's one-line S5 gate
("calibration runs end to end") would have shipped a plumbing job that locked in a schema nobody had
checked. Owner ruling, 2026-07-30: **derive first, then implement, in phases.**

This file is the running log of that derivation. Two things live here:

| § | what |
|---|---|
| **§1** | ⭐ **The accumulator changes S5 needs** — the upstream list, which reopens S3/S4. A1 `Σ L` · **A4 the node weight** |
| **§2** | The evolving S5 implementation plan, superseding `IMPLEMENTATION_PLAN.md` §4/§5's S5 row |
| §3 | The derivation: what was measured, how, and what is still open |
| §4 | Decisions taken, with their reasons — so they are not re-litigated |

Read `ACCUMULATOR_DESIGN.md` §6 first; this file is its consequence.

---

## §0 STATUS

| | |
|---|---|
| stage | **S5 complete through S5.f, and S6 has landed.** S5.g is the next modelling step |
| tooling | `scripts/design/observable_efficiency.py` — the harness; re-run it, do not quote this file |
| suite | ⭐ **1752 passed / 1 failed** (S5.f: −234 failures; S6: the goldens regenerated). The single failure is an owner modelling call, not a defect — `TODO.md` §7 |
| ⭐ **THE PIVOT IS PASSED** | **`calibrate()` RUNS end to end** on all 8 chr22 pilot conditions, and **the FIRST BASELINE is recorded in `LEDGER.md`'s S5.f entry** — bit-identical on re-run. Every deferred derivation now has something to be A/B'd against |
| ⛔ blocker | **none. S5.g (A7's taper) is next** — the first change that gets a real A/B. See §2's road map |

### ⚠ Three things the first baseline says, and they are not all comfortable

1. **It carries A7's known bias by ruling** — an 11.0 % genome-wide gDNA over-call, +0.36 in the last node
   before a polyA site. That is what S5.g removes, and measuring the removal is why A7 was deferred.
2. ⛔ **The fitted κ is `1 − truth`.** A library simulated at `strand_specificity = 0.99` calibrates to
   **κ = 0.0101**. This answers `CARRY_FORWARD.md` §0 **C4** against ground truth for the first time: the
   near-zero sense fraction on all four real cfRNA libraries is a **convention flip, not biology**.
   Pre-existing — S5.f does not touch `fit_strand_balance` — and it needs its own step.
3. ⚠ **Calibration is bit-identical run to run; the EM samples from the posterior BY DESIGN.** ⭐ Run an
   end-to-end A/B under `assignment_mode="map"` or `"fractional"` — spread drops ~0.5 % → ~1e-10 — and
   hold the mode fixed across both arms, since the three modes are different estimators.

---

## §1 ⭐ ACCUMULATOR CHANGES S5 NEEDS

Each entry says what, why, and what it costs. ⚠ Anything here reopens S3 (the C++), S4 (the payload)
and their byte-identity gates against `tests/native/_accumulator_reference.py`.

### A1 — ✅ **RULED IN for nodes (owner, 2026-07-30): store `Σ L`.** ⛔ The EDGE case is still open — see A5.

**What.** A third per-object, per-strand accumulator: the summed molecule length `Σ L` over the
fragments that deposit on that object. Same populations as today (node contained, node spanning,
contiguous-edge unspliced, contiguous-edge spliced, junction).

**Why — measured, `scripts/design/observable_efficiency.py`.** Fragments arrive as a Poisson process,
so the full length histogram at an object is the sufficient statistic and its Fisher information about
the gDNA share has a closed form. That makes "how much does this storage choice throw away?" exactly
computable as `efficiency = Var(full histogram) / Var(stored set)`.

Swept over gDNA × RNA fragment lengths **independently**, means 50–300 bp × cv 0.15–1.00 (756 ordered
pairs per frame), in three shape families:

| frame | `count` only | **`count, Σ1/L`** (ships) | `count, ΣL` | **`count, Σ1/L, ΣL`** |
|---|---|---|---|---|
| node 151 bp (the median node) | 0.596 / **0.000** | 0.832 / **0.078** | 0.904 / **0.043** | **0.953 / 0.188** |
| node 1000 bp | 0.000 / 0.000 | 0.197 / **0.001** | 0.680 / 0.000 | **0.722 / 0.064** |
| node 3000 bp | 0.000 / 0.000 | 0.182 / **0.000** | 0.646 / 0.000 | **0.692 / 0.047** |
| contiguous edge | 0.000 / 0.000 | 0.324 / **0.000** | 0.535 / 0.000 | **0.686 / 0.078** |

*(median / minimum over the grid, gamma family. The ranking is unchanged under lognormal and normal.)*

⭐ **The shipped pair has an EXACT blind spot, and it is structural.** At an edge the count row is
`(μ_g − 1, μ_r − 1)` and the density row is `(1, 1)`, so the 2×2 determinant is `μ_g − μ_r`. **When the
two components have equal mean length the pair carries literally zero information about the split, at
any depth** — while the full histogram still separates them on spread alone. Two populations with the
same mean length and different breadth is an ordinary configuration, not a corner. Adding `Σ L` removes
the blind spot: the minimum over the whole grid goes `0.000 → 0.078`.

⭐ **And `Σ1/L` collapses on long nodes** — 0.18 median at a 3000 bp node against 0.65 for `Σ L`. Long
intergenic nodes are exactly where the pure gDNA background is measured, so this is not a corner either.

**Why a third moment can work at all.** Given the count, the fully efficient one-dimensional statistic
is `Σ h*(L)` with `h*(w) = (f_g(w) − f_r(w)) / (φ f_g(w) + (1−φ) f_r(w))`. That weight depends on the
fitted length models, and **the pools are built in the same pass, so `h*` is not available at deposit
time** — the weight must be fixed before any library is seen. But a stored triple spans `{1, 1/w, w}`,
and the coefficients of a linear combination *can* be chosen after the models are fit. The efficiency
numbers above are exactly how well the optimal score is approximated in that fixed basis.

**Cost.** `uint64[n, 2]` per population. At human scale with the current five populations:
`(1,043,881×2 + 1,043,595×2 + 404,168) × 2 × 8 B` ≈ **80 MB per worker**, against the 110 MB the S3
design already budgets and an 8.6 GB measured peak RSS. Deposit cost: one add per object touched, no
division (unlike the density quantum), so cheaper per object than what is already there.

⚠ **Open before this is ruled.** `Σ L` is an integer and needs no fixed point, so determinism is free —
but the headroom must be checked: at `L ≤ 2000` and 10⁸ fragments the worst case is 2e11, comfortable
in uint64, but the bound should be stated in the header the way `DENSITY_SCALE`'s is.

### A5 — ✅ **RULED IN (owner, 2026-07-30): edges store `Σ L` too.**

The R-b ruling names the **node** accumulator. The edge was excluded earlier on the grounds that
`1/(L−1)` already gives everything needed there — **and that is correct for the LEVEL**, by the identity
in `NODE_DENSITY_DERIVATION.md` §4.2. It is not correct for the SPLIT.

| contiguous edge | median | **min** |
|---|---|---|
| `count, Σ1/(L−1)` — ships | 0.324 | **0.000** |
| `count, Σ1/(L−1), Σ L` | **0.686** | **0.078** |

⭐ The edge is where the equal-means blind spot is **exact**: the count row is `(μ_g−1, μ_r−1)`, the
density row is `(1,1)`, so the determinant is `μ_g − μ_r` and the pair carries *literally zero*
information about the split whenever the two components share a mean length — at any depth, in every
shape family, at every composition. Nodes never have a blind spot this clean, because their count
coefficient `E[(ℓ−w+1)₊]` depends on the whole distribution below `ℓ`, not only the mean.

✅ **Ruled: edges included.** +70 MB/worker across all populations, +6 % of the measured 8.6 GB peak.

⚠ The junction-edge frame has not been measured separately (its reach geometry differs), but a junction
is a 0-bp jump and structurally an edge, so the same determinant argument applies to it unchanged.

### A4 — ⭐ The node density weight is the wrong argument. **`docs/NODE_DENSITY_DERIVATION.md`**

**What.** The deposit passes `L` to `density_quantum` at a node. The derivation shows the weight must be
the population's own **opportunity**: `(ℓ − L + 1)` when contained, `(L − ℓ − 1)` when spanning. The
edge's `L − 1` is already that rule; the node's `L` is not.

**Why.** `h = 1/A` is the *unique* weight making `E[Σh] ∝ ρ` with a length-distribution-independent
constant (T1, verified: 0 violations against 266 for `1/L`). It reproduces both known limits exactly —
at `ℓ = 0` the node rule becomes the edge rule identically, at `ℓ ≫ FL` it becomes `count/ℓ` — and
contained + spanning together give `ρ·(1 − f(ℓ+1))` to machine precision. It applies the effective-length
correction **per fragment, with no fitted model**: at a 1000 bp node with 200 bp fragments `count/ℓ` reads
0.80× truth and `Σ1/A` reads 1.000000×.

⛔ **But it does NOT help the composition split, and that is a theorem, not a measurement.** Model-free
means both components share a coefficient, so the channel's row is `(K, K)` and its determinant against
anything is zero. Measured: `Σ1/A` alone scores 0.113 efficiency at a 151 bp node and 0.000 at 1000 bp.
**The level and the split are two jobs; no single number does both.**

⚠ **So `Σ1/A` ADDS to the schema, it does not replace `Σ1/L`.** Swapping loses information at short
nodes (0.832 → 0.684 at 151 bp); adding gains most at long ones (0.692 → 0.832 at 3000 bp).

**Cost.** No new machinery — `density_quantum(placements)` already takes a placement count and the node
length is in the cut array the deposit already searches. ✅ **The headroom was re-priced and is a
non-issue**: `A = 1` needs a fragment to exactly fill a node, and even at 10⁸ such fragments the sum
reaches 2.3 % of uint64. A second, smaller scale for this channel is unnecessary and would put two scales
in one schema. Storage is +70 MB/worker per channel across all populations.

⚠ **Recommendation, after pricing: R-b, not R-c.** `Σ1/A` buys the model-free level but almost no
information (0.953 → 0.960 median at a 151 bp node), and its architectural claim is weaker than it first
appeared — the *split* still needs the fitted FL models, so it removes the length model from one of two
paths, not from calibration. Its one real gain is at long nodes (0.692 → 0.832), where the gDNA
background is anchored. `NODE_DENSITY_DERIVATION.md` §7.1–7.2.

### A6 — ⚠ **OPEN: is `spanning` informative GIVEN the adjacent edges?** (owner, 2026-07-30)

Not a blocker for S5.a, but it questions whether a population survives at all.

**The doubt.** A fragment that spans node `i` crosses *both* of `i`'s lines and is counted on both
edges. So what does the node's spanning counter add? And the `ℓ+1` boundary is perplexing: at `w = ℓ`
the fragment is contained, at `w = ℓ+1` it vanishes from the node entirely and is visible only as edge
flux, and it stays invisible to the node however much longer it gets — until `w ≥ ℓ+2`, when it
reappears as spanning.

**What is settled.** Spanning is *not* a re-measurement of the edges. An edge's expected count is
`ρ·(μ−1)`, **independent of all geometry** — every edge has the same expectation, so edge counts carry
no node-length information at all. Spanning is `ρ·E[(w−ℓ−1)₊]`, which traces the FL **survival
function** as `ℓ` varies. Different functional, genuinely different information.
⚠ Note `ACCUMULATOR_DESIGN.md` §5.1 *withdrew* the original justification (pricing adjacent-edge
correlation) and rests spanning on this length-bracketing argument instead. That argument holds.

⛔ **What is NOT settled, and it is the owner's actual question.** At a short node
`spanning(ℓ) ≈ edge_count − ρℓ`, so once the solver has the two flanking edges, spanning may add
little. The measured 0.000 → 0.758 at a 25 bp node compared *contained-only* against
*contained + spanning* — it never gave the alternative the flanking edges. **The comparison that
answers this has not been run.**

⚠ It needs harness work, not just another cell: spanning is a **subset** of edge-crossing, so the two
populations are not disjoint, and `observable_efficiency.var_set` zeroes cross-population covariance —
which is exactly wrong for this pair. The overlap must be modelled before the number means anything.

**Why it is not blocking:** `length_sum` is added per population; if spanning is later dropped, its
channel goes with it. The architectural question is real but downstream.

### A7 — ✅ **RULED (owner, 2026-07-30): reach enters at JUNCTION edges in S5.e; contiguous edges stay UNBOUNDED until after S5.f**

A junction edge is certainly reach-governed — only a spliced molecule uses it. But an **unspliced**
crossing is a gDNA/RNA *mixture*: the RNA part is constrained by its transcript's remaining length, the
gDNA part is not (its template is the chromosome). So reach enters the unspliced divisor **per
component**, not as a property of the edge — which is what makes it complicated, and is why
`taper_g = 1` while `taper_r` comes from the annotation (`ACCUMULATOR_DESIGN.md` §7).

⭐ **The decision reduces to ONE call site**, and that is what made it rulable. gDNA at *any* edge takes
`UNBOUNDED_REACH` on both sides (`taper_g = 1`, i.e. `mu_g − 1`), and the two contained frames take no
reach argument at all. So the only open question was what `(reach_lo, reach_hi)` the RNA component gets
at a **contiguous** edge.

| edge kind | component | reach in S5.e |
|---|---|---|
| contiguous | gDNA | `UNBOUNDED_REACH` — settled by physics, not by this ruling |
| **contiguous** | **RNA** | ⭐ **`UNBOUNDED_REACH`** — the ruling |
| junction | RNA (the only component) | ⭐ **the real per-strand EXONIC reach**, from `edges_df` |

**Why the contiguous edge stays unbounded.** `mu_r − 1` is what production does today, so S5.e varies
exactly one thing: the faces dissolve. A7 then lands as its **own step after S5.f**, which is the only
ordering in which it gets a real A/B — against the first baseline S5.f produces. Doing it inside S5.e
would land a 213-line dissolution *and* a first-ever behaviour change together, with no baseline to
attribute either against.

⚠ **The price is known and must not be forgotten**: `CARRY_FORWARD.md` §1 fact 6 — ignoring the taper
over-calls gDNA by **11.0 %** genome-wide, contiguous seams are worse than junctions (0.750 vs 0.886),
and the gDNA fraction is off by **+0.36** in the last node before a polyA site. **The first baseline
therefore carries a known bias, and the S5.f ledger entry must say so.**

**Why the junction edge takes its real reach now.** A junction edge is a **brand-new population** — the
predecessor had no junction divisor at all (spliced mass went to boundary faces on the reach-free
half-triangle, deleted in S5.c), so there is no number to regress. It also means the reach plumbing is
**exercised rather than dead**: one code path (`crossing_eff_length`), and "does this edge taper?" is
answered by which array is passed, not by a flag with one live branch. Leaving it unbounded would ship a
divisor known to be wrong by up to **4×** at a first exon (`CARRY_FORWARD.md` §2: 199.0 at R=550 against
50.0 at R=50).

### A2 — Two alternatives were tested and are NOT recommended.

* **A coarse length histogram instead of moments.** The sufficient statistic, pure integer counts (no
  fixed point at all), and it needs no deposit weight chosen at scan time — the exact `1/A` can be
  applied post hoc per node. At **equal bytes** it wins at long nodes and edges and loses at short ones,
  because global bins cannot resolve a 25 bp node. Does not dominate; a far larger change.
  `NODE_DENSITY_DERIVATION.md` §7.5.
* **float64 instead of uint64 fixed point.** ⚠ The design's stated justification compares against
  **float32** and is not sound as written — float64 is the same 8 bytes and would be ~10⁵× less
  nondeterministic than the measured figure the decision rests on. The argument that survives is that
  **byte-identity is the S3 gate**, and a tolerance is what hid this project's factor-of-2 bug for
  months. That is a testing argument, not a numerical one. `NODE_DENSITY_DERIVATION.md` §7.4.

`node_spanning_*` (§1 A3 below) needs no accumulator change — it is already stored. The junction edges
are already stored. No other new observable has been shown to earn its bytes.

### A3 — ⭐ `node_spanning_*` is already stored and consumed by NOBODY. Wire it. **No upstream change.**

Not an accumulator change; recorded here because it is the largest single win found and it is free.
Efficiency at a node, contained-only against contained + spanning (gamma):

| node | contained only | **+ spanning** |
|---|---|---|
| 25 bp | 0.000 | **0.758** |
| 100 bp | 0.001 | **0.852** |
| 151 bp | 0.014 | **0.832** |
| 1000 bp | 0.190 | 0.197 *(spanning is empty here, as expected)* |

⭐ **At every node shorter than one fragment, essentially ALL the information is in the spanning
population** — and 56.7 % of human nodes are shorter than a 200 bp fragment. `ACCUMULATOR_DESIGN.md`
§6 asserts this ("short nodes are measured by what spans them"); this is the first time it has been
quantified as an information share.

---

## §2 THE EVOLVING S5

⚠ **Supersedes `IMPLEMENTATION_PLAN.md` §5's S5 row and §4's R1–R4 ranking.** Those describe a
rewiring; what is actually needed is a rewiring *plus* two new information sources *plus* possibly an
upstream schema change. Phases, in dependency order:

| | phase | gate | status |
|---|---|---|---|
| **S5.0** | **The derivation** (§3) — settle A1 | efficiency measured over the FL grid, MC-validated, opportunity formulas enumerated | ✅ **A1 measured** over 756 pool pairs × 3 shape families × 4 compositions. Awaiting the owner's ruling |
| **S5.a** | ✅ **DONE** (`LEDGER.md`) — `Σ L` added to the reference, the C++ and the payload. Node weight unchanged (`1/L`), so the deposit change is purely ADDITIVE — one new array per population, no quantum change | byte-identity to the reference on LBX0190 + MO_3021; bit-identity at 1/2/4/8 workers; S3/S4 gates re-run; deposit cost re-recorded | ✅ byte-identical on 60 k real fragments / 44 refs; 5 reference + 2 C++ perturbations caught; suite failures unchanged at 266; deposit **339.7 ns/frag** |
| **S5.a2** | ⚠ **How `Σ L` ENTERS THE SOLVE** — a separate derivation, and it must not be folded into S5.f | the length-composition likelihood derived, falsification-tested, A/B'd against the first baseline | not started |
| **S5.b** | ✅ **DONE** (`LEDGER.md`) — `fl.py` re-keyed to the five pure pools | ✅ strict xfail XPASSed and was removed; 4 perturbations caught; the §8 pool table reproduced on real cfRNA | done |
| **S5.c** | ✅ **DONE** (`LEDGER.md`) — `effective_length.py` → the one placements formula | ✅ 89 tests, every formula enumerated; reproduces `CARRY_FORWARD.md` §2's taper table | done |
| **S5.d** | ✅ **DONE** (`LEDGER.md`) — substrate → ONE type; the chain re-keyed to `N E N … E N` | ✅ 24 tests, 7 perturbations caught incl. the old `k+1` terminal shape | done |
| **S5.e-1** | `NodeGeometry` + `build_node_geometry` rewritten; the 18 per-face arrays dissolve (old R1b, R2) | ✅ 23 tests written first and verified failing; **15 perturbations, and P15 found a real hole** — a `2*i` slot-layout assumption that every single-reference fixture was blind to | ✅ **done** |
| **S5.e-2** | `bp_solver`'s `(left, right)` tuples collapsed through its six consumers; `node_init` + `density_model` re-keyed | ✅ the solver runs end to end; the factor-1 bedrock invariant holds on the new geometry. ⚠ the scalar/vector twins are NOT merged (a measured 15.7×/op) | ✅ **done** |
| **S5.e-3** | the last per-face test fixtures ported; every transitional shim deleted | ✅ `test_bp_solver.py` 19 failures → **0**; the port exposed two artefacts the old shape hid (a hand-placed mature flux, and a terminal-slot G1 lock that had its own ten-line apology comment) | ✅ **done** |
| **S5.f** | ✅ **DONE** (`LEDGER.md`) — `calibrate` + `CalibrationResult` + `priors`/`capture_eff_length`/`pipeline`, and the junction axis exported | ✅ **calibration runs end to end on all 8 chr22 pilot conditions and the numbers ARE the first baseline**, bit-identical on re-run; `tests/calibration/` green (543/0); 9 perturbations, **P6 found a real hole** (a dead `~either_ambig` guard reading as the rule) | ✅ **done** |
| **S5.g** | ⭐ **A7 proper** — the contiguous-edge RNA reach turned on (§1 A7) | falsification test first; an A/B against S5.f's baseline, which is the whole reason it is sequenced here. ⚠ **A/B on the CALIBRATION numbers, not end to end**: calibration is bit-identical run to run but a downstream transcript count is not (§0) | ⭐ **NEXT** |

⚠ **The 266 suite failures cannot be a per-step gate.** ~200 are end-to-end scenario and golden tests
that will move *numerically*, not merely start running. Each phase is gated on the unit tests written
for it; the suite goes green at S5.f and the goldens are regenerated **once**, at S6.

### ⭐ THE ROAD TO A PRODUCTION CALIBRATION

⛔ **Everything measurable is downstream of one gate.** `calibrate()` does not run, so there is no
number: no baseline, no A/B, no benchmark verdict, no scan-cache toy seed. **S5.f is therefore the
pivot** — every deferred derivation below is deliberately sequenced *after* it, because each one needs a
baseline to be judged against and would otherwise land unfalsifiable.

**Phase 1 — get a number (no new modelling; plumbing only).**

| | step | gate | note |
|---|---|---|---|
| 1 | ✅ **S5.e-rest** — the 10 sweep-behaviour tests ported; `_pending_s5e` deleted | ✅ the file is green |
| 2 | ✅ **S5.f DONE** — `calibrate()` + `CalibrationResult` + `priors` + `derive` + `capture_eff_length` + `pipeline` + the junction axis | ✅ **the FIRST BASELINE, in `LEDGER.md`'s S5.f entry** | ✅ also retired `region_arrays`'s `k+1` boundary↔region mapping, the last surviving old axis, and `strand_deconv`'s two-seeds-per-boundary |

#### ⭐ S5.f, scoped against the actual call sites (2026-07-30)

Every remaining calibration failure is S5.f's, and they are **28 tests over 8 files**: `test_calibrate`
(9), `test_gdna_strand_integration` (4), `test_spliced_boundary_onesidedness` (3),
`test_region_index_alignment` (3), `test_accumulator_span_unbiased` (3), `test_substrate_conservation`
(2), `test_oracle` (2), `test_ambig_scenario` (2).

**The schema change, and it is the whole step.** `CalibrationResult` carries **six** per-region mass
arrays — `mass_{gdna,rna}_{contained,left,right}` — and `priors.assemble_priors` immediately re-pools
two of them as `mass_gdna_right[r] + mass_gdna_left[r+1]`, with `capture_eff_length._pooled_seam_arrays`
doing the identical thing. ⛔ **That split-then-re-pool is a no-op with a history**: `CARRY_FORWARD.md`
§3 trap 2 records the same sum-then-halve pattern hiding an exact factor of 2 for months. Owner ruling
(§4): the left/right pair becomes **per-edge arrays**, and the pooling disappears rather than being
re-derived.

| field | fate |
|---|---|
| `mass_{gdna,rna}_contained` | keep, on the NODE axis — they already are per node |
| `mass_{gdna,rna}_{left,right}` | ⭐ **→ `mass_{gdna,rna}_edge`**, per contiguous edge. `bp_solver.chain_edge_deconv` already returns exactly this |
| `mass_rna_spliced` | keep; it is what `assemble_priors` WITHHOLDS from `rna_prior_count` (a spliced fragment has no gDNA candidate in the EM, so counting it would double it) |
| `gdna_region_eff_len` | → `effective_length.contained_eff_length`, a direct successor, same quantity |
| ⛔ `gdna_boundary_len` | **no successor.** It is `boundary_side_eff_length` = `E[min(ℓ,L)]/2`, deleted in S5.c because it divided a per-FACE mass. Its replacement is the per-edge `crossing_eff_length`, one number, and the `½·(bl[r] + bl[r+1])` averaging in `priors` and `capture_eff_length` goes with it |

**Consumers to rewire, in dependency order:** `result.py` (the schema) → `calibrate.py` (the body; its
dead code is already the specification of what to rewire) → `priors.py` + `capture_eff_length.py` (the
pooling) → `derive.py`, `track.py` (contained only — likely untouched) → `pipeline.py:897` (a QC sum
over the three per-region nodes → two axes).

⚠ **`region_arrays.boundary_region_indices` / `region_boundary_indices` are the last surviving `k+1`
axis in the tree** and must go with this step, or the old shape outlives the schema that needed it.

⚠ **Test fixtures that hand-build a `CalibrationResult` will all move**: `test_result_schema`,
`test_priors`, `test_capture_eff_length`, `tests/calibration/_oracle.py`. Several encode the D6 ½
convention in their comments (`gdna_boundary_len` "IS the halved per-side density length") — that
convention dies with the field, and those comments must not be carried across.

⛔ **DO NOT LAND THE SCHEMA ALONE.** It was written, measured and **reverted** on 2026-07-30: with
`calibrate()` still raising, the new fields are consumed by nobody (speculative code) while every
hand-built fixture fails on the old kwargs — **+35 failures, 256 → 291, none of them informative**.
The schema and its consumers are ONE step. The target shape, so the next session types it rather than
re-derives it:

```python
@dataclass(frozen=True, slots=True)
class CalibrationResult:
    mass_gdna_node:        np.ndarray   # float64[n_nodes]   <- chain_node_deconv
    mass_rna_node:         np.ndarray
    mass_gdna_edge:        np.ndarray   # float64[n_edges]   <- chain_edge_deconv
    mass_rna_edge:         np.ndarray
    mass_rna_spliced_edge: np.ndarray   # the edge_spliced part of mass_rna_edge, WITHHELD by the prior
    gdna_node_eff_len:     np.ndarray   # contained_eff_length on the gDNA pmf
    gdna_edge_eff_len:     np.ndarray   # crossing_eff_length  on the gDNA pmf -- one number per line
    gdna_density_global: float
    rna_sense_frac: float
    gdna_strand_overdispersion: float
    rna_strand_overdispersion: float
    n_nodes: int
    n_edges: int
    config: CalibrationConfig
```

⚠ `mass_rna_spliced` has **no node twin**, and that is structural rather than an omission:
`node_contained` is credited only when the fragment used no junction, so a node's contained population
cannot hold a spliced molecule.

⚠ **Four more modules take the OLD substrate views and were not in the first scope** — found by grep
while attempting the step: `gdna_strand.fit_gdna_strand_from_substrate`,
`density_deconv.fit_intron_background`, `background_reference.measure_background`, and
`calibrate._build_intron_prior`. All four read `substrate.contained.n_unspliced*` /
`substrate.left`/`right`. Budget for them.

| 3 | **S6** — delete the dead paths; regenerate `tests/golden/` **once** | suite green | ~200 of the 266 remaining failures are end-to-end/golden tests that will move NUMERICALLY, not merely start running |

⚠ **Do not fold any modelling change into 1–3.** The whole point of the ordering is that step 2's output
is the reference every later change is measured against, so it must be a pure rewiring of what already
exists — known biases included, and named.

**Phase 2 — the deferred derivations, now A/B-able against that baseline.**

| | step | what it buys, measured | blocked on |
|---|---|---|---|
| 4 | **S5.g — A7**: the contiguous-edge RNA reach taper | removes an **11.0 %** genome-wide gDNA over-call and **+0.36** in the last node before a polyA site (`CARRY_FORWARD.md` §1 fact 6) | a baseline. The plumbing is one array; the ruling is already recorded |
| 5 | **S5.a2 — how `length_sum` enters the solve** | it is stored on every population and consumed by NOBODY. It is the channel that removes the equal-means blind spot: efficiency min `0.000 → 0.078` at an edge, `0.078 → 0.188` at a 151 bp node | its own derivation; must not be folded into S5.f |
| 6 | **A6 then A3 — `node_spanning`** | ⭐ the largest single win found and it needs no upstream change: at every node shorter than one fragment essentially ALL the information is in the spanning population (0.000 → 0.758 at 25 bp), and **56.7 %** of human nodes are shorter than one 200 bp fragment | ⚠ **A6 first**: spanning is a SUBSET of edge-crossing, so `observable_efficiency.var_set`'s zero cross-population covariance is exactly wrong for the pair. The overlap must be modelled before the number means anything |
| 7 | **Does `spliced_count` enter the LEVEL?** | it is a contiguous crossing, so in principle it belongs in ν and in ρ_tot; today it enters only the strand solve, which is a faithful port of the predecessor, not a derived choice | S5.a2's frame |

**Phase 3 — productionise.**

| | step | gate |
|---|---|---|
| 8 | close the benchmark suite's two known gaps: **(c)** replicate pairs so overdispersion is estimable at all, **(f)** conditions in the **1–10 %** gDNA band real libraries live in | `suite_resolves.py` green on all 8 requirements — ⛔ run it before quoting any suite number |
| 9 | the real-cfRNA gate at full scale + `scan_profile.py` regression | the S3 byte-identity gate re-run; ns/fragment re-recorded |
| 10 | the **magic-number audit** | `POOL_EB_PRIOR_ESS = 1000` shifts the fitted pure-gDNA mean **79.4 → 100.3 (+26 %)** against a 4,467-fragment pool; the gDNA strand overdispersion **saturates its 0.2 ceiling on 2 of 4 real libraries** at 66–600 σ, so it is bias and not noise. Both are owner calls, not tuning |

⚠ **The two biggest risks to this plan, both already measured.** (a) The simulator is **Poisson by
construction** (`sim/wgs_engine.py:473`), so nothing dispersion-dependent validates on the suite until
step 8 lands — and the first fix considered provably cannot work (`TODO.md`). (b) The suite's own pure
pools separate by **1.20×** against LBX0190's 2.5×, and that separation is the determinant of every 2×2
in the deconvolution — so a suite result is not a substitute for the real-cfRNA gate at step 9.

### Structural deltas the rewiring must absorb (unchanged by the derivation)

| | old | new |
|---|---|---|
| axes | regions `R`, boundaries `B = R + n_refs` | nodes `N`, contiguous edges `E = N − n_refs`, junction edges `J` |
| chain | `B R B R … R B` — data-free terminal boundaries at both ends | `N E N E … E N` — starts and ends with a **node**; no terminal slot exists |
| channels | 4 columns `[unspl+, unspl−, spl_sense, spl_antisense]` — two conventions in one array | 2 columns = genome strand, without exception. Sense derived from a junction's own strand |
| numerator | fractional per-face `mass` + integer `flux` | integer `count` + fixed-point density (+ `Σ L`, pending A1) |
| faces | a boundary's two sides have *different* divisors | an edge is a 0-bp line with one set of numbers — every face concept dissolves |

---

## §3 THE DERIVATION

### 3.1 Method

`scripts/design/observable_efficiency.py`. Fragments are a Poisson process, so at one object the count
of length-`w` fragments in population `p` is `Poisson(opp_p(w)·(ρ_g f_g(w) + ρ_r f_r(w)))`, independent
across `w`. Hence:

* the **full length histogram is sufficient**, and its Fisher information is closed-form and exact;
* any stored summary is a projection of it, so `efficiency = Var_full / Var_set ∈ [0,1]` is exactly the
  fraction of available information the storage choice keeps.

Grid: gDNA and RNA fragment lengths **independently**, means `{50, 75, 100, 150, 200, 250, 300}` bp ×
cv `{0.15, 0.35, 0.60, 1.00}` — 756 ordered pairs per frame, **both directions**. ⭐ Owner ruling: there
is no rule that RNA is longer than gDNA; it happens to hold for cfRNA and assuming it is how a tool
overfits to one library type. Three shape families (gamma, lognormal, normal), because shape turned out
to matter — see 3.3.

### 3.2 What is verified, and what is approximated

| | |
|---|---|
| opportunity formulas `(ℓ−w+1)₊`, `(w−ℓ−1)₊`, `(w−1)` | **exact enumeration** of integer start positions, deterministic, no tolerance |
| `Var_full` | **exact** (Poisson closed form) |
| `Var_set` | Gaussian / known-covariance — **approximate**, and it flatters a heavy-tailed channel |
| the harness itself | **monotonicity self-check**: a superset of channels can never carry less information |

⭐ **The monotonicity check earned its place immediately.** The first run reported the triple at 0.000
where the shipped pair scored 0.832 — impossible, and a covariance-conditioning failure (channel scales
differ by ~1e10, so a plain condition-number guard rejected informative sets for being badly scaled).
Without the check the write-up would have argued from it. **Falsification discipline applied to a
derivation harness, not only to production code.**

### 3.3 ⚠ Two results that were wrong before they were right

1. **An assumed fragment-length shape reversed the answer.** On gamma pools `(count, ΣL)` beat
   `(count, Σ1/L)` at every node length by 1.5–2.3×. Carrying the *same moments* on a measured
   histogram from the suite scan caches, the advantage vanished. Shape is a variable in this
   derivation, never a background assumption — hence the three families.
2. **The Monte Carlo check was flattering itself.** A first version clipped `φ̂` to `[0,1]` and solved
   by unweighted least squares; the clip truncates the spread and the wrong solver is less efficient
   than the Fisher calculation assumes, so the simulation came out **2× better** than the theory it
   existed to test — precisely on the heavy-tailed channel under suspicion. Corrected to GLS at the
   true covariance with no clipping.

### 3.4 What the corrected Monte Carlo says

Robust (IQR) scales match the Gaussian prediction within ~10 % everywhere, so **the grid ranking is
sound as a bulk statement**. Raw standard deviations do not: `count + ΣL` at an edge gives a realised
sd of 187.5 against a predicted 0.375, with a robust scale of 0.374.

⭐ That is not a defect of `Σ L` — `count + Σ1/L` shows the same thing (realised 8.77, robust 0.404) —
it is the ratio estimator's tail at realistic per-object depth, and it restates
`ACCUMULATOR_DESIGN.md` §6 in the node frame: **there is no per-object answer.** At the median 151 bp
node the entire evidence is ~0.6 contained + ~4.1 spanning fragments. The information numbers above are
therefore about what each object contributes to a **pooled** estimate — which is the right question,
since information is additive — and never about solving one object alone.

### 3.5 The composition sweep — ✅ closed

`h*` depends on the composition `φ_g`, so the fixed basis `{1, 1/w, w}` could have approximated it well
only near `φ = ½`. Swept at `φ_g ∈ {0.05, 0.20, 0.50, 0.80}` (median / min over the grid, gamma):

| frame | set | φ=0.05 | φ=0.20 | φ=0.50 | φ=0.80 |
|---|---|---|---|---|---|
| node 151 bp | ships `count, Σ1/L` | 0.750 / 0.033 | 0.819 / 0.057 | 0.832 / 0.078 | 0.819 / 0.057 |
| node 151 bp | **`count, Σ1/L, ΣL`** | **0.913 / 0.102** | **0.946 / 0.133** | **0.953 / 0.188** | **0.946 / 0.133** |
| node 3000 bp | ships | 0.146 / 0.000 | 0.156 / 0.000 | 0.182 / 0.000 | 0.156 / 0.000 |
| node 3000 bp | **triple** | **0.602 / 0.018** | **0.645 / 0.028** | **0.692 / 0.047** | **0.645 / 0.028** |
| contiguous edge | ships | 0.217 / **0.000** | 0.308 / **0.000** | 0.324 / **0.000** | 0.308 / **0.000** |
| contiguous edge | **triple** | **0.689 / 0.033** | **0.684 / 0.052** | **0.686 / 0.078** | **0.684 / 0.052** |

⭐ **The ranking is unchanged at every composition**, and the shipped pair's zero minimum at an edge is
present at all of them. Low gDNA is the *worst* case for the shipped pair and the triple's margin is
widest there — 0.146 → 0.602 at a long node, a 4.1× gain, in exactly the 1–10 % band real libraries
live in (`TODO.md` §2(f)).

⚠ `φ = 0.20` and `φ = 0.80` are identical by construction, not by coincidence: the grid contains both
orderings of every pool pair, so swapping the two components is the same operation as `φ → 1−φ`.

### 3.6 ⛔ Open

1. **Junction-edge and node-spanning frames are not yet in the grid** as separate rows with their own
   reach geometry.
2. **`CARRY_FORWARD.md` §1 fact 5 does not reproduce** and should be corrected or retired. It claims
   `(count, Σ1/L)` beats `(count, ΣL)` "everywhere, 1.6× more precise, 4× better conditioned, never
   worse". It was measured at gDNA 50 / RNA 200 — a 4× mean separation, deep in the region where the
   `μ_g − μ_r` determinant is large. Over the full grid it is false at every node ≥ 250 bp and at the
   edge, and it is exactly reversed in the equal-means family.
3. **The suite's pure pools disagree with the suite's own configuration.** Configured at gDNA
   156.5 ± 124.6 / RNA 206.1 ± 98.3 (`BENCHMARK_SUITE.md` §2); `pool_lengths` in the scan caches
   measure **195.4 ± 97.7 / 234.9 ± 146.8**. The RNA shift is explained and is `ACCUMULATOR_DESIGN.md`
   §8.1(b)'s junction opportunity bias — longer fragments cross more junctions — here quantified for
   the first time at **+14 % in the mean, +50 % in the sd**. The gDNA intergenic shift 156.5 → 195.4 is
   **not** explained by containment bias (intergenic nodes are ~11.5 kb, worth 0.2 %) and has no
   account. ⚠ It also means the suite's pure-pool separation is **1.20×**, against the 2.5× design §8
   quotes for LBX0190 — and the separation is the determinant of every 2×2 in §6.

---

## §4 DECISIONS TAKEN

| | ruling | why |
|---|---|---|
| **2026-07-30** | ⭐ **A7: junction edges take their real per-strand EXONIC reach in S5.e; contiguous edges pass `UNBOUNDED_REACH` for BOTH components, and A7 proper lands after S5.f** | The decision reduces to one call site (contiguous-edge RNA) — gDNA is `taper_g = 1` everywhere and the contained frames take no reach. Keeping it unbounded leaves S5.e varying exactly ONE thing and is the only ordering where A7 gets a real A/B, against the first baseline S5.f produces. Junctions are a NEW population with no predecessor divisor, so wiring their reach regresses nothing and keeps the one code path exercised instead of dead. ⚠ The first baseline therefore carries the known 11.0 % genome-wide gDNA over-call. Owner |
| **2026-07-30** | ⭐ **The three channels are `count` / `inv_length_sum` / `length_sum`** = `Σ1`, `Σ1/placements`, `ΣL`. ⚠ This RENAMES the shipped `density` | `Σ1/L` at a node is not a density — `ACCUMULATOR_DESIGN.md` §6 says so itself ("a better-conditioned second moment and nothing more"); it is an exact density only at an edge. Naming it `density` puts one word on two concepts, which is `CARRY_FORWARD.md` §3 trap 27. The three names are three sums and are honest at every object. Owner |
| **2026-07-30** | ⭐ **EDGES store `Σ L` as well as nodes** — every population carries the same triple | The edge is where the equal-means blind spot is EXACT (determinant `μ_g − μ_r`, efficiency min 0.000). Excluding it would leave the one hard failure unfixed. Owner |
| **2026-07-30** | ⭐ **NODES accumulate `(count, Σ1/L, Σ L)` — R-b.** The node density weight stays `1/L`; `Σ1/A` is NOT stored | `Σ L` fixes the equal-means blind spot, the one hard failure (min at a 151 bp node 0.078 → 0.188). `Σ1/A` buys a model-free *level* but almost no information (0.953 → 0.960) and its architectural case is partial — the split still needs the FL models. `Σ L` is a plain integer sum: no fixed point, no scale, no overflow scheme. +70 MB/worker, +6 % of peak RSS. Owner |
| **2026-07-30** | **Derive before implementing S5; phase the implementation** | The accumulator was rebuilt to enable new observables. Wiring the consumers without knowing which observables to use would freeze a schema nobody checked. Owner |
| **2026-07-30** | **The FL grid is swept in BOTH directions, gDNA and RNA independent** | There is no rule that RNA is longer than gDNA — it holds for cfRNA and fails elsewhere. cfRNA is one extreme of the spectrum (low input, sparse, zero-inflated); pure cell-line RNA-seq is a different beast. Owner |
| **2026-07-30** | `CalibrationResult`'s per-region `mass_*_left/right` → **per-edge arrays** | `priors.py` pools the two halves straight back together, so splitting and re-pooling is a no-op — and `CARRY_FORWARD.md` §3 trap 2 records that this exact sum-then-halve pattern hid a factor of 2 for months. Owner |
| **2026-07-30** | `build_boundary_flags_array` → a plain per-contiguous-edge array | The `k+1` padded-terminal axis exists only to match the old boundary slots. `E = N − n_refs` needs no padding, and the off-by-one commentary in `splice_graph.py` goes with it |

### Not yet decided

* **A1** — whether `Σ L` is stored. ⭐ The measurement is complete and one-sided.
* **A4** — whether the node weight becomes the reciprocal opportunity, and whether `Σ1/A` is stored
  alongside `Σ1/L` rather than replacing it. `NODE_DENSITY_DERIVATION.md` §7.2 gives three costed
  options; ⭐ **R-b (`count, Σ1/L, ΣL`) is the recommendation** — it fixes the one hard failure, adds a
  single channel that is a plain integer sum with no fixed point at all, and costs +6 % of peak RSS.
* Whether the junction edges become factors on their endpoint nodes in S5 or later
  (`CLAUDE.md`: factors, never message channels — the graph is not a polytree).
* Whether `ACCUMULATOR_DESIGN.md` §8.1(b)'s per-pool opportunity correction lands with S5.b or after.
