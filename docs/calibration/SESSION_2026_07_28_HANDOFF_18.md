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

> **✅ RE-VERIFIED 2026-07-28 (independent 4-agent audit, 0/3 refuters on each line). Two strengthenings.**
>
> 1. **Re-measured at HEAD on all four cached cfRNA payloads**, not quoted from commit messages:
>    `od_g` = **0.200000 / 0.003149 / 0.200000 / 0.092339**, `od_r` = 0.008576 / 0.013683 / 0.007362 /
>    0.013412 (LBX0190 / LBX0588 / MO_3021 / vcap). Saturating on **2/4**, unchanged.
> 2. **The weight is not "≈ 1.0" — it is IDENTICALLY 1.0**, measured on the *real production seeds* across
>    both the region and boundary channels on all four libraries: `frac(w < 1−1e-9) = 0`, `min−1 ≈ −2.2e-16`.
>    So `mean_s = ½·w + κ·(1−w)` is always exactly `½`: **the entire gDNA/RNA mixture in that fit is dead
>    code on real data.** The eff-length cancels (`density_model.py:132` divides by `region_eff_len`,
>    `:162-164` multiplies it back) and `substrate.py:101-102` sets contained mass ≡ contained count.
> 3. **Provenance closed**: `git log -L :_region_seeds:` tops out at `bef315ec`/`978ad9ed`, both far older
>    than this session. `24b09cb9` touched the gDNA half **only inside docstrings**. ⚠ If you remember
>    "we fixed the overdispersion saturation", that memory is the **`od_r`** fix.

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

### R3 — ⛔ **IMPLEMENTED, MEASURED, AND REFUTED (2026-07-28). Ships default OFF.**

**⚠ First, this item conflated two different results.** §1's table labels `daa32a13` "N2", and that commit
shipped the **prior BOOTSTRAP** (`calib_refit_iters` 1 → 3) — settled, converged, monotone. What was
*unimplemented* is the other N2 result: **admitting AMBIG into the final hyperprior fit's substrate**. Those
are separate questions (`gdna_hyperprior_production_plan.md:993` says so explicitly) and only the second was
ever open.

**⚠ And the sales pitch below was scored against a retired baseline.** "fabrication 14.5× → 6.7×" compares
against the plan's **pass-0 / refit=1** row. HEAD runs refit=3, so its baseline is the plan's own
"#2 re-solved, AMBIG out" row — and read against *that*, admitting AMBIG wins recovery (0.562 → 0.757) and
EMD (0.2945 → 0.2544) but **LOSES the fabrication guard 4.7× → 6.7×**. The 14.5× → 4.7× was the bootstrap's
own gain. The cost claim was also wrong: three fits and three sweeps already ship, so this is **not** an
extra fit — it is a substrate predicate.

**Implemented** as `CalibrationConfig.calib_admit_ambig_final` (default **False**), gating
`_fit_gdna_hyperprior(admit_ambig=...)` on the **final** fit only. That shape —
`pass-0 → fit #1 (out) → re-solve → … → final fit (IN) → FINAL solve` — is non-circular *at any depth*,
because the admitting fit trains on AMBIG estimates produced by a prior that never saw AMBIG. (Admitting on
a non-final fit would break exactly that; the naive "iterations ≥ 2" predicate is circular from fit #3 on.)

**Measured, four paired arms on all 32 conditions** (each condition scanned once, all arms on the same
payload; pre-registration in `scratchpad/R3_PREREG.md`, written before any number was seen):

| stratum | A `r3,out` (HEAD) | B `r2,IN` | C `r3,IN` | D `r2,out` | **B−D** (flag) | **C−A** (flag) |
|---|---|---|---|---|---|---|
| **ALL 32** | **0.046675** | 0.050253 | 0.048393 | 0.047846 | **+0.002407** | **+0.001718** |
| unstranded | 0.062989 | 0.068834 | 0.065685 | 0.064891 | +0.003943 | +0.002696 |
| capture ON | 0.070245 | 0.075581 | 0.072495 | 0.071988 | +0.003593 | +0.002250 |
| **VSTRONG** | 0.079067 | 0.096601 | 0.084634 | 0.083102 | **+0.013499** | **+0.005568** |
| unstr × capON | 0.089575 | 0.097766 | 0.092843 | 0.092254 | +0.005512 | +0.003268 |
| **`gdna_none`** (FP guard) | **0.010766** | 0.012717 | 0.011065 | 0.011304 | **+0.001413** | **+0.000299** |

Per-condition: **26/32 worse** at depth 2, **25/32 worse** at depth 3. `n_train` 1179 → 1395 (+18.3 %,
matching the plan's +19 %), so the flag is genuinely active — it is not a null result.

**Against the pre-registered predictions:**
* **P1 (the flag is non-inert on the synthetic) — CONFIRMED.** +18.3 % training population, deltas to 0.0247.
* **P2 (the `gdna_none` FP guard REGRESSES) — CONFIRMED**, in both comparisons. This was the deciding test,
  and it agrees with the plan's own fabrication column (4.7× → 6.7×) measured on a different metric.
* **P3 (AMBIG-heavy strata improve) — FALSIFIED, and inverted.** VSTRONG and unstranded × capON are the
  *worst*-hit strata, not the best. Whatever the extra AMBIG training nodes contribute, it is not accuracy
  where the AMBIG mass lives.

**⇒ VERDICT: do not ship.** The pre-registered decision rule was "ship only if ALL-32 improves **and**
`gdna_none` does not regress"; ALL-32 gets worse and the guard regresses. The recovery/EMD gains the plan
measured are real but **do not convert into better deconvolution**, and they are bought with a worse
false-positive guard on two independent metrics. Default stays `False`.

**Side result:** `A − D = −0.001170` on ALL-32 independently re-confirms **depth 3 > depth 2**, i.e. the
shipped `calib_refit_iters = 3`.

⚠ Real-data transfer was doubtful anyway: 71–97 % of real training nodes are zero-mass anchors weighted 1.0
each, so live AMBIG would grow the weighted training population by only **+0.4–3.0 %** on real cfRNA against
**+18.3 %** here.

### R4 — ✅ **FIXED 2026-07-28.** The `z2` units mismatch was real, and correcting it INVERTS the work list

`var_gdna` exceeds 0.25 on **33 %** of nodes — impossible for the variance of a fraction — so it is a
**log-space** variance, while `z2` compares it against a **linear** squared error. **Every `z2` number in
the ROADMAP and in HANDOFF_1..18 is on a mixed scale**, so its absolute calibration point ("1.0 = honest") is
not meaningful; only the direction of change is. ✅ The production reliability weight is *consistent*
(log-space `var` against `ref = 1/g + S0`); only the diagnostic is broken.

**Confirmed at the source, not inferred from the 0.25 bound.** `simplex_logodds._solve_nodes_logodds:394-396`
computes `Lg = _log_fg(lam)` then takes grid moments of `Lg` — so `var_gdna` **is** `Var(log f_g)`, exactly.
The `> 0.25` argument was right but weaker than reading the line.

**Root cause: two docstrings asserted the wrong units** and are what generated this item. Both fixed:
`strand_deconv.NodeDeconv.gdna_frac_var` said `Var(f_g)`; `node_geometry.NodeBelief` said `Var(f_c)`.
`pass0_error_table.py` was itself self-contradictory — line 100 said `Var(log f_g)`, line 113 said `Var(f_g)`.

**Production re-verified consistent on BOTH consumers**, including the one that genuinely needs linear units:
`gdna_landscape._reliability` documents and uses `v = Var(log f_g)`; and `composition_logvar`, which needs the
LINEAR `Var(f_g)`, is **not** fed `belief.var_gdna` — `bp_solver.node_sweep` derives `_var_fg =
(f_g(1−f_g))²/τ_λ` fresh, which is exactly `f_g²·Var(log f_g)` given `Var(log f_g) = (1−f_g)²/τ`. The
delta-method step was already done correctly there.

**The fix** (`pass0_error_table._lin_var`): convert the denominator to linear via the exact lognormal moment
`Var(f_g) = f_g²·(e^v − 1)`, capped at `f_g(1−f_g)` (the greatest variance any `[0,1]` variable with mean
`f_g` can have — the same bound `bp_solver` applies). No tuned constant; the numerator is left alone because a
log-space numerator needs `log(oracle)` and the oracle `f_g` is **exactly 0 on 23.5 %** of scored nodes.

#### ⭐ The consequence: the z2 ordering was EXACTLY BACKWARDS

Measured on all 32 conditions at refit=0. The per-class conversion factor spans **14.8×–287.9×**, so this
re-orders rather than rescales. Suite total **0.046 → 1.007** — pass-0 is *honestly calibrated*, not the
20×-conservative it appeared.

| node class | share of ERR | `z2` OLD | `z2` **NEW** | rank |
|---|---|---|---|---|
| intron AMBIG | 0.4 % | 0.010 | **2.874** | #6 → **#1** |
| intron single | 2.0 % | 0.011 | **2.145** | #5 → **#2** |
| boundary AMBIG | 3.9 % | 0.104 | 1.543 | #1 → #3 |
| exon AMBIG | 26.0 % | 0.069 | 1.140 | #2 → #4 |
| boundary single | 6.0 % | 0.059 | 0.993 | #3 → #5 |
| exon single | 61.8 % | 0.037 | 0.918 | #4 → #6 |

⭐ **The classes carrying 88 % of the error mass are honest** (0.92 / 1.14) — their error is a MODE defect,
not a precision over-credit, so precision work cannot reach it. **The only over-confident classes are the
INTRONS** (2.1 / 2.9), which §3 below calls solved because their *mwae* is lowest. They are accurate and
**overstate their certainty 2–3×**, which points at the intron factory's curvature precision
(`density_deconv.density_factor_precision`) being over-credited. That is a new, small, well-localised lead —
and it is the opposite of where the old metric pointed.

⚠ `pass0_error_diagnostic.py` splits on `var < 1.0` as a confident/under-determined threshold. That cut was
chosen when `var` was believed linear (where `< 1.0` would select *everything*, since a fraction's variance
cannot exceed ¼). It still functions as a data split — ~76 % of nodes sit below it — but it is not the
"confident" set it reads as. `pass0_error_table.py` uses quartiles instead and is unaffected.

⚠ Also fixed: `pass0_error_table.py` hard-coded the **shared** `/tmp/rigel_selfsolve` work dir; it is now
`--work`, defaulting to the old path.

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
