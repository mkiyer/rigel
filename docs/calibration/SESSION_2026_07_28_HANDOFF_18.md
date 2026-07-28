# SESSION 2026-07-28 — HANDOFF 18 (LIVE): what remains before calibration can ship

**Supersedes `SESSION_2026_07_27_HANDOFF_17.md`.** Everything in `src/` is **committed**; the branch is
`calib-ambig-init-wip`, 19 commits ahead of `main`, nothing merged.

---

## 0. STATE — measured at HEAD, this session

| | |
|---|---|
| suite mwae (32-condition `ambig_dense_10mb`) | **refit=0 `0.079005` · refit=3 `0.046675`** |
| tests | **1213 pass**, ruff clean on `src/ tests/ scripts/` |
| goldens | current |
| `od_r` (RNA strand overdispersion), real cfRNA | **0.0086 / 0.0137 / 0.0074 / 0.0134** — off the ceiling on 4/4 ✅ |
| ⛔ `od_g` (gDNA strand overdispersion), real cfRNA | **0.2000 / 0.0031 / 0.2000 / 0.0923** — still **saturating on 2/4** |

⚠ **Re-record every A/B baseline from the current tree in the same session.** Every stored number goes
stale the moment the tree moves; if HEAD-vs-baseline is not 32/32 the *baseline* is broken.

### Per-stratum, at HEAD

| stratum | refit=0 | refit=3 |
|---|---|---|
| ALL 32 | 0.079005 | **0.046675** |
| stranded | 0.029066 | 0.023684 |
| unstranded | 0.114439 | 0.062989 |
| capture ON | 0.103860 | 0.068466 |
| capture OFF | 0.033535 | 0.018132 |
| **VSTRONG** | **0.179467** | 0.079067 |
| **unstranded × capON** | **0.151302** | **0.093479** |
| `gdna_none` (FP guard) | 0.069253 | 0.010766 |

## 1. WHAT LANDED (this session, in order)

| commit | what |
|---|---|
| `3fd8fc59` | **W4** — the `GdnaLandscape` hyperprior, the additive ψ gDNA arm, and the `_pin_v` BP fix |
| `daa32a13` | **N2** — the prior BOOTSTRAP: `calib_refit_iters` 1 → **3**, exposed as `--calib-refit-iters` |
| `c99c399f` | strand-overdispersion **constants cleanup** + the **information-currency** shrinkage |
| `24b09cb9` | **the per-junction SJ strand table** — `od_r` now fitted from one population with κ |
| + docs | `pin_derivation.md` §12, `gdna_reframe_terminus.md`, `P1G_SCOPE.md`, `strand_overdispersion_design.md`, `sj_strand_table_design.md` |

## 2. ⭐ THE REMAINING WORK, RANKED

### R1 — ⛔ `od_g` still saturates on 2/4 real libraries. **The gDNA seed channel is contaminated.**

The RNA half is fixed; **the gDNA half is not.** `strand_overdispersion_design.md` §2 D1: the seed's
`gdna_weight` is `count_gdna_frac`, which for a count-observable region is read from the region's own
unspliced counts and is therefore **≈ 1.0 by construction** — it is not a measurement of gDNA-ness, it is a
**re-encoding of "the annotation called this intergenic"**. Measured top seed on `vcap`: `N = 1523`,
`sense = 5` (sense fraction **0.003**) at `gdna_weight = 1.000`, carrying **12 % of the pooled numerator**.
That is a transcript.

⚠ **A robust estimator is NOT the fix, and this is measured** (`strand_overdispersion_design.md` §3b):
`od_mom = 0.9313` on LBX0190 sits **66–600 σ above the ceiling** — it is **bias, not power**, so no
shrinkage or trimming rule reaches it. ⛔ And a ceiling-based seed screen **rejects zero seeds** at
`a ∈ {2,3}` on all four libraries.

**The fix is upstream**: the gDNA weight must be a measurement. That is the accumulator/index rework's
territory (path-level FL contrast identifies gDNA vs RNA prior-free and strand-free).
**⇒ R1 is BLOCKED on upstream, and the ceiling must stay until it lands** — it is currently the only thing
between the strand likelihood and an `od_g` of 0.93.

### R2 — ⛔ The gDNA REFRAME at transcript termini. **The largest measured pass-0 defect.**

`gdna_reframe_terminus.md`. The decomposition is an identity (residual 4.4e-16): the sources are nearly
right (+0.110 dec), gDNA really is uniform (−0.054), and **the reframe is 96 % of the error (+1.508 of
+1.564)**. At a splice junction `r ≈ 1` and the gDNA transfer is **exact (1.0×)**; at a transcript terminus
no RNA crosses, the boundary face is pure gDNA, and `r` becomes the exon's whole RNA-to-gDNA ratio
(**7.0×**). **33 % of edges carry 66–68 % of the error mass.** Closed form, verified to 1e-9 on 63–66 % of
terminus edges: the delivered gDNA collapses to `ρ_tot(dst)`, i.e. too large by exactly `1/f_g(dst)`.

**Sized:** `r_g ≡ 1` takes unstranded × capOFF × gDNA **0.0495 → 0.0146** at r0 (6/6 better) and capture-OFF
overall −0.0178 — while costing unstranded × capON **+0.1728**, so it bounds the prize and is not a
candidate. **⇒ This is P1g** (`P1G_SCOPE.md`), and per memory it is **subsumed by the accumulator redesign**.

### R3 — the N2 iterative two-fit architecture. **Measured to win; an owner decision, not a task.**

Prior #2 fitted on the *re-solved* AMBIG nodes is non-circular by construction and beats what ships on
every axis: enriched recovery **0.572 → 0.757**, EMD-vs-non-AMBIG 0.2930 → 0.2544, **fabrication
14.5× → 6.7×**. Cost: one extra prior fit + sweep per library. **Not implemented — say yes or no.**

### R4 — ⚠ the `z2` trust metric has a UNITS MISMATCH, and it has been steering the work

`var_gdna` exceeds 0.25 on **33 %** of nodes — impossible for the variance of a fraction — so it is a
**log-space** variance, while `z2` compares it against a **linear** squared error. **Every `z2` number in
the ROADMAP and the handoffs is on a mixed scale**, so its absolute calibration point ("1.0 = honest") is
not meaningful; only the direction of change is. ✅ The production reliability weight is *consistent*
(log-space `var` against `ref = 1/g + S0`); only the diagnostic is broken. **Cheap to fix, and it is the
metric the pass-0 work list is ordered by.**

### R5 — validation gaps

* **`gdna_benchmark_5mb` has never been run** for any of this work — the only held-out synthetic suite.
* **The Gate-0 injection harness is not built.** The synthetic suite has `od = 0` **by construction**, so it
  can only reward estimators that return zero — it validates **one side** of a two-sided question. Inject
  Beta-Binomial dispersion at each library's *real* seed-size distribution and pre-register both
  (i) recover a non-zero `od_true` and (ii) do not manufacture dispersion at `od_true = 0`.
* **Real data has no accuracy oracle.** Four cfRNA payloads run; the only quasi-truth is LBX0190's
  known ~15 % gDNA, against which the prior moves `f_g` 0.236 → 0.124.

### R6 — carried reserves, not blockers

* **`_S0`** holds the entire remaining enriched census (recovery 0.69 → 0.84 at capture-ON, 0.29 → 1.24 at
  VSTRONG with flat weights). ⚠ It cannot serve both ends — the census wants it large, the FP guard small.
* **P1d / P1e** remain deliberately retained and measurably redundant; re-evaluate after P1g.
* **AMBIG exons** are 26 % of error; a trained prior is the designed fix, so largely downstream.

## 3. WHERE THE ERROR IS NOW (refit=3, by node class — the ordering that should drive work)

Exons are the whole problem: **exon single-strand ≈ 59 %** and **exon AMBIG ≈ 26 %** of all error, boundary
≈ 11 %, **intron ≈ 3.6 %** (the intron factory works). The worst strata are **unstranded × capON (0.0935)**
and **VSTRONG (0.0791)**.

## 4. CAN IT SHIP?

**`src/` is fully committed and green.** What stands between here and `main`:

1. ⛔ **`od_g` saturating on 2/4 real libraries** (R1) — the strand likelihood is our strongest evidence and
   it is being applied clipped on half the real samples. Blocked upstream, but it must be *decided*
   (ship with the ceiling as a guard + QC canary, or wait).
2. ⛔ **R2 is a known, sized, unfixed mode error of up to 190×** on a named population.
3. ⚠ **No held-out validation** (R5) — `gdna_benchmark_5mb` unrun, no injection harness.
4. **R3 is an unanswered owner decision.**

**Recommendation: land the branch on `main` as WIP-complete** — it is green, self-consistent, documented,
and every remaining item is either blocked upstream or an explicit decision. Do **not** describe calibration
as production-ready: `od_g`, R2 and R5 are all open, and the docs say so.

## 5. ⛔ DO NOT RE-LITIGATE (measured and settled this session)

* The **P-2 residual** is closed as not-fixable-at-the-pin; two repair arms measured and failed.
* The **share transfer** (`pin_derivation.md` (★)) is implemented, measured and refuted; the refutation is
  inline in `bp_solver.py`.
* **`Beta(3,3)` as a seed screen** — void: the "measurably better" table applied LBX0190's multiplicity to a
  vcap seed. At `a ∈ {2,3}` every multiple-testing frame rejects **zero** seeds on all four libraries.
* **`od₀ = od_max/2`** — rejected; it breaks a tested bug fix and collapses strand log-evidence 305 → 114 nats.
* **"More dispersion = weaker"** is FALSE for `od`: it is the gDNA component's *specificity at its own
  mean*, so inflating it degrades **both** directions.
* **Refitting the RNA mean from its own seeds** — refuted.
* **Crediting all K junctions**, and the **1/K** split (provably biased, 4–12 σ).
* The **100 %-gDNA initialization default** is inert; the emitting boundaries are 96.4 % intron↔exon.
