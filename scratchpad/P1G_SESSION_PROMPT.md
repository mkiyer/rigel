We are continuing Rigel calibration — **P1g: putting TSS/TES into the region/boundary map.** Read, in order:

- `docs/calibration/ROADMAP.md` — the single entry point. Item 8 in §6 is P1g.
- `docs/calibration/P1G_SCOPE.md` — **the brief for this session.** Read it in full; it is short.
- `docs/calibration/gdna_reframe_terminus.md` — the measurement P1g acts on (the node-level trace).
- `docs/calibration/pin_derivation.md` §12 — how this defect was found, and the two repair arms that failed.
- `docs/calibration/SESSION_2026_07_27_HANDOFF_17.md` — its "NEXT SESSION" block only.

Do not read `docs/calibration/archive/` for guidance.

Setup: `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`, always
`OMP_NUM_THREADS=1`, repo `/Users/mkiyer/proj/rigel`. ⚠ `/tmp/rigel_selfsolve` is shared and
non-namespaced — two sessions doing scenario work corrupt each other.

## Where things stand

**Nothing is committed.** Uncommitted in `src/`: `bp_solver.py`, `calibrate.py`, `simplex_logodds.py`,
`config.py` modified; `calibration/gdna_landscape.py` new. Suite mass-weighted mwae over the 32-condition
`ambig_dense_10mb` battery: **refit=0 0.078786 / refit=1 0.052470**. Gates green: ruff clean on
`src/ tests/ scripts/`, **1214 tests pass, 0 failures — the goldens ARE regenerated at this tree.**

⚠ **Every stored A/B baseline is stale the moment the tree moves.** Re-record from the current tree in the
same session at both refits before comparing anything. **If HEAD-vs-baseline is not 32/32, the BASELINE is
broken, not your change** (this cost 40 minutes once).

## THE MOTIVATION — why P1g is now the highest-value item

The region map records a 4-bit signature `{intron±, exon±}` and a boundary row is
`[boundary_id, ref_name, position]`. So the solver **cannot tell a splice junction from a transcript
terminus** — physically opposite seams, treated identically. This was already the structural debt behind
`ω_graft` (which splits ≥30× on that bit). It is now measured to do a second, larger kind of damage.

The gDNA density a message delivers into an exon is **5.4× too large** at capture-OFF, and the
decomposition is an identity (residual 4.4e-16): the sources are nearly right (+0.110 dec), gDNA really is
uniform (−0.054), and **the reframe is 96 % of the error (+1.508 of +1.564)**. Split by the structural bit:

| source boundary | n | median `log10 r` | median delivered gDNA error |
|---|---|---|---|
| **TERMINUS** | 426 | +0.836 | **7.0× too big** |
| **junction-only** | 732 | +0.021 | **1.0× — exact** |

33 % of edges carry **66–68 % of the error mass**. Reproduced on `gdna300` and on the stranded twin.

**The closed form, and it predicts the measurement to 3 %.** At a terminus no RNA crosses, so
`ρ_tot(src) = ρ_g(src)` and, with `e` the true capture step, `r = e/f_g(dst)` and
`tg = ρ_g(src)·r = ρ_tot(dst)`. The delivered gDNA is too large by **exactly `1/f_g(dst)`** — the
reciprocal of the destination's own composition. Verified: `tg == ρ_tot(dst)` to 1e-9 on 63–66 % of
terminus edges; node 2651 predicted 185× / measured 190×; node 3087 predicted 125× / measured 128×.

## THE GOAL

Land the structural bit and use it. Recommended order (P1G_SCOPE §10):

1. **Restore the v6 boundary columns + version bump.** ⭐ This is a *restoration*: index format v6 already
   carried `is_tss`/`is_tes`/`is_splice_junction`/`genomic_sj_strand`; v7 deleted them when their consumer
   was retired. Recover the builder verbatim from
   `git show 1863ef57^:src/rigel/calibration/regions.py`. **Gate: bit-identical 32/32** (nobody reads it yet).
2. **Thread it to `node_geometry`.** ⚠ This is the real work — **the calibration package never reads
   `boundary_df` today.** Prefer a small `BoundaryAnnotation` passed alongside `region_arrays` over
   widening `RegionArrays` (which touches every call site incl. the debug scripts). **Gate: bit-identical 32/32.**
3. ⭐ **Measure the bit on the REAL cfRNA index before any solver work.** The synthetic annotation is
   nested truncations plus exon skips — no alternative TSS/TES *inside* exons (which would not even produce
   a region boundary). **This is the step that could invalidate the plan, and it is cheap.**
4. **C1 — the reframe at terminus seams.** The prize, with the pre-registered A/B of P1G_SCOPE §8.
5. **C2 — re-derive `ω_graft` per structural class** (same estimating equation, one scalar per class; the
   partial-pooling block in `graft_premise_logvar` is the plug-in point).
6. C3 (`accept_l/accept_r` from structure rather than observation), C4 (P1e's scope), then **goldens LAST**.

**No C++ change is needed** — the accumulator keys on `boundary_id` and its `(boundary_positions,
ref_pos_offsets, region_types)` ABI is untouched by extra `boundary_df` columns. But **every index must be
rebuilt** (the run suites, the cfRNA index, the test fixtures); budget for that, it is the slow part.

## ⛔ Settled — do not re-litigate

- ⛔ **Do not derive the bit from the signature instead of restoring it.** Measured: 83.3 % agreement, and
  the **73 false negatives are exactly the positions that are BOTH a terminus and a junction** (verified as
  a subset) — where one transcript ends at a point another spans, the *union* signature is structurally
  blind, on precisely the mixed seams that matter most. The 211 false positives are 100 % nRNA-span edges.
- ⛔ **Do not "correct" `r` by `f_g(dst)`.** The closed form says that is exactly the missing factor, and it
  is the destination's own belief — that rebuilds the `_pin_v` bug this tree just fixed. **The fix must be
  structural (do not form that ratio), not corrective.** This is the trap.
- ⛔ **`r_g ≡ 1` is a SIZING arm, not a candidate.** It takes unstranded × capOFF × gDNA
  0.0495 → 0.0146 (r0) and 0.0313 → 0.0077 (r1), 6/6 better, capture-OFF overall −0.0178/−0.0124 (9 b/0 w)
  — while costing unstranded × capON **+0.1728/+0.2069**. Any operator must reduce to `r` where a capture
  step is genuinely crossed.
- ⛔ **Not a variance.** The delivered mode is wrong by up to 190×; damping cannot move a mode toward truth.
  That is the P1e lesson (`variance_ledger.md` §6) and it applies verbatim.
- ⛔ **The P-2 residual is closed** — not fixable at the pin; two repair arms measured and failed
  (`pin_derivation.md` §12). Do not re-open it; P1g is what it pointed at.
- ⛔ The share transfer (`pin_derivation.md` (★)) was implemented and refuted; the refutation is inline in
  `bp_solver.py`.

⭐ **The precedent to build on:** `bp_solver`'s λ-emission gate already says *"a pure-gDNA source's
authority is a DENSITY LEVEL, and that travels on the measurement stream (`tmg`)"*. That is correct — **and
the measurement stream is still multiplied by `r`.** The three-stream separation is right in design; the
reframe being applied uniformly across all three is the defect.

## How to run

```bash
export P0_BENCH_OUT=/tmp/pass0_bench.tsv
OMP_NUM_THREADS=1 P0_REFIT=0 python scripts/debug/pass0_oracle_bench.py --arm base_r0
OMP_NUM_THREADS=1 P0_REFIT=1 python scripts/debug/pass0_oracle_bench.py --arm base_r1
# ...make the change, then re-run with --arm new_r0 / new_r1
OMP_NUM_THREADS=1 python scratchpad/p2res_report.py $P0_BENCH_OUT base_r0 new_r0 --per-cond
```

Tools built for this defect (all reproduce the numbers above):
`scratchpad/p2r_g_trace.py` (per-node forensic trace: transcript structure, counts, densities, messages,
and the origin decomposition), `scratchpad/p2r_h_terminus.py` (the terminus/junction split),
`scratchpad/p2r_i_plot.py` (the figure `docs/calibration/figures/gdna_reframe_terminus.png`),
`scripts/debug/pass0_node_dissect.py` (exact-replay ψ channel ablation — ⚠ check its per-node replay
fidelity before trusting an attribution).

## Methodology (standing)

No magic numbers — pause and discuss before any new constant. Vary one thing per A/B arm. **Pre-register
predictions including a falsification test before measuring** (for C1 the falsification test is that
junction-only edges must be UNMOVED — they are already exact, so movement there means the bit is applied
too widely). The synthetic suite is Poisson by construction, so nothing overdispersion-dependent can be
validated on it. Goldens LAST. The owner drives commits. I like seeing plots — visualize fits when
investigating.

## Also live, but NOT this session's task

- ⛔⛔ **The strand overdispersion saturates its `Beta(2,2)` ceiling on every real library** (`od_r` = 0.2000
  on 4/4 cfRNA samples, `od_g` on 2/4; the synthetic suite fits 0.0008–0.0017). The only intrinsic
  gDNA/RNA information source is applied over-confident on all real data. Owner call — the `a ≥ 2` rule has
  a stated justification and data now contradicts it.
- **N2's verdict reversed** (plan §"N2 RE-RUN AFTER W4"): the iterative two-fit architecture beats what
  ships on every axis (enriched recovery 0.572 → 0.757, fabrication 14.5× → 6.7×). Owner decision, since it
  ships a second prior fit.
- `gdna_benchmark_5mb` has still never been run for this work.
