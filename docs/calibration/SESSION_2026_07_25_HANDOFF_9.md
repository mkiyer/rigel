# Session handoff — the AMBIG study (⛔ §2–§4 CORRECTED by §0b — read §0b FIRST)

**This is the LIVE handoff. Read `docs/calibration/ROADMAP.md` first, then this.** It supersedes
`SESSION_2026_07_25_HANDOFF_8.md` (whose §6 "next study: AMBIG" is this file). Do NOT read
`docs/calibration/archive/`. Date: 2026-07-25. Branch `calib-ambig-init-wip`, HEAD `b0955493` (M8).

---

## 0b. ⛔ CORRECTION (owner, same day) — §2–§4 BELOW OVERSTATE THE CASE. READ THIS FIRST.

The owner rejected the estimator §2–§4 are built on, and the MC arbiter (`scratchpad/ambig_8_mc.py`, laws
M9a–M9e, all passing) confirms it. **Three claims below are wrong or overstated:**

1. **`d̂ = (p_obs−½)/(κ−½)` is degenerate and is not an admissible estimator.** `κ` is FITTED — posterior
   `Beta(n_same+1, n_opp+1)`, so `Var(κ) = κ(1−κ)/(n_obs+3)`, already computed as
   `StrandBalance.rna_strand_overdispersion` (currently labelled "QC-only, NOT fed into the deconv"). On
   unstranded data `|κ̂−½|` is **smaller than its own posterior sd** (measured `|κ−½|/sd` = **0.0–1.3**;
   `|κ̂−½|` = 0.00004 vs sd 0.00088) and it sits in a **denominator**. The observed excess is not zero either
   — pure gDNA splits 506/494, not 500/500 (med `|e|` = 0.0018–0.0047). Result: med `|d̂|` up to **118**,
   **80–99 % of mass with `|d̂| > 1`**, and `B` clipped to 0 asserting "100 % RNA" against an oracle `f_g` of
   0.272. That is noise amplification, not a bound. (`scratchpad/ambig_7_kappa_degeneracy.py`)

2. **The SHIPPED solver was never degenerate — only my estimator was.** M9c: at κ≈½ the plug-in-κ
   *likelihood* is already flat (max/min = 1.000) because `A ≪ σ_e` regardless; marginalizing κ gives 1.001.
   So the shipped 2-D solve marginalizes the tilt correctly and never divides. §4's "catastrophic on
   unstranded" is a property of `cons`/`B`, **not** of the solver, and §2's framing of the unstranded failure
   ("d̂ is noise") mis-attributes it: the cause is the division by a quantity smaller than its own error.

3. **The residual AMBIG bias is the TILT PRIOR, not discarded information.** M9d: the *correct Bayesian
   marginal* under the solver's own uniform-in-θ tilt reference is itself biased **−0.10 to −0.15 when
   τ = 1** — which is exactly the 0.09–0.14 under-estimate §1/§2 attribute to the solver failing to extract
   the bound. It is not failing to extract anything. `B` scores well in §3/§4 **only because τ ≈ ±1 in this
   suite** (the minority strand at an `exon+intron(x-strand)` region is the unbaited intron) — a POPULATION
   fact, precisely what pass-0 may not assume. **So §3's "85–91 % of the AMBIG prize" is an offline
   substitution with an inadmissible estimator, not an achievable prior-free target.**

**What SURVIVES.** The census (§1) — AMBIG is 31.6 % of suite error mass, `τ_own = 0` on every AMBIG node,
messages roughly double AMBIG error, and the λ message is the dominant message defect (bit-exact ablation
0.1444 → 0.0986). The algebra `f_g ≤ 1 − |d|` is true. §5's two reverted A/B arms stand as measured.

**THE CORRECTED FORMULATION (derived + MC-validated, M9a/M9b).** Never form `d̂`. The strand likelihood
depends on `f_g` and `κ` **only** through the achievable strand-excess amplitude

```
    A(f_g, κ) = (1 − f_g)·|κ − ½|          L(f_g) = ∫∫ BetaBin(u₊; n, ½ + A·τ, ω)·π(τ)·π(κ) dτ dκ
```

(the sign of `κ−½` is absorbed by τ's symmetry — M9a, exact to 0.00 %). No division anywhere, so there is no
degeneracy: as `κ → ½`, `A → 0` for EVERY `f_g` and the factor flattens — **"no information" is recovered as a
LIMIT, not asserted by a guard**, which also subsumes the existing `disc` deadband. In the sharp limit
(`σ_e ≪ A`, κ known) it has the closed form `p(f_g) ∝ 1/√((1−f_g)² − d²)` on `[0, 1−|d|]` — an *integrable
spike* at the bound, not a hard wall — with median `1 − |d|·cosh(arccosh(1/|d|)/2)` (M9b, rel 0.02–0.08 %).

**THE REAL LEVER, therefore, is the TILT MESSAGE, not a bound.** What converts the constraint into an
estimate is knowing τ ≈ ±1, and the prior-free source for that is the neighbours' per-strand RNA imputation —
which the solver already carries as `theta_imp` and which the ablation shows is **nearly inert**
(0.1444 → 0.1417 when removed). That is the next thing to investigate, and it replaces §6's plan.

---

## 0. TL;DR

**AMBIG is not an information problem after all — on stranded data.** The received view (ROADMAP §5) is that
AMBIG nodes have *no* composition evidence prior-free (`τ_own = 0`, from the rank-1 Schur complement) and must
wait for the hyperprior. That is an **interior** argument, and it is incomplete. The strand mean is exactly
`p = ½ + (κ−½)·d` with `d = f₊ − f₋`, so the tilt `d` is **directly observed**, and the simplex gives a hard
structural bound `f_g ≤ B = 1 − |d̂|`. Zero Fisher information is not zero knowledge when the profile
likelihood has **bounded support**.

Measured: under capture, `B` is an almost **unbiased** estimate of `f_g` on AMBIG nodes (mass-weighted
0.7894 vs oracle 0.7896), while the solver sits 0.09–0.14 **below** it. Snapping AMBIG to `B` would cut suite
mwae by 34–66 % on stranded capture-ON and captures **85–91 % of the entire AMBIG prize**.

`B` alone is not the answer (it asserts the minority strand is silent — a population claim, and it is
catastrophic on unstranded data). But the exact identity `f_g = 1 − |d| − 2·min(f₊, f₋)` closes it
**prior-free**, because the minority strand's RNA is already imputed by the solver's own per-strand messages.
That estimator is validated offline in §4: it **halves AMBIG error on every stranded condition** and beats the
bare bound exactly where the minority strand really carries RNA.

**Nothing was landed this session.** Two candidate fixes were A/B'd and both were ~neutral (§5); the real fix
is a ψ-factor and needs its precision derived + MC-validated first (§6). HEAD remains M8:
**0.0900 (refit=0) / 0.0700 (refit=1)**, gates green.

## 1. The map, re-measured on HEAD (M8)

`scripts/debug/suite_dissect.py --out /tmp/suite_nodes_m8.npz` — 74,494 nodes, total error mass **12.56 M**
(was 12.92 M pre-M8), suite mwae 0.0900.

| stratum | nodes | err-mass | share | self | solved | Δ msg |
|---|---|---|---|---|---|---|
| solvable & single-strand | 51,719 | 8.60 M | **68.4 %** | 0.2080 | 0.0904 | −0.1176 |
| solvable & AMBIG | 10,971 | 3.96 M | **31.6 %** | 0.2686 | 0.1838 | −0.0849 |

⚠ **Scope correction to HANDOFF_7 §3.** "AMBIG ≈ 50 % of error mass" is the *stranded-conditions* census; the
**suite-wide** number is **31.6 %**. Both are right in their scope — quote the scope.

Per-class, stranded capture-ON (`scratchpad/ambig_1_census.py`, raw signature bits × DOF):

| class | DOF | share of err | self | solved | oracle f_g | τ_own |
|---|---|---|---|---|---|---|
| exon | single | 18.3 % | 0.0049 | 0.0112 | 0.746 | 29.07 |
| **exon** | **AMBIG** | **34.8 %** | 0.0958 | **0.1920** | 0.787 | **0.000** |
| **exon+intron(x-strand)** | **AMBIG** | **33.8 %** | 0.0909 | **0.1169** | 0.775 | **0.000** |
| BOUNDARY | AMBIG | 6.2 % | 0.0681 | **0.1713** | 0.858 | **0.000** |

Two facts: **`τ_own = 0` on every AMBIG node**, and **messages roughly double AMBIG error** (self 0.090 →
solved 0.144–0.192) while they *halve* it on single-strand.

## 2. ⭐ THE CENTRAL RESULT — an AMBIG node is NOT opinion-free

gDNA is strand-symmetric, so it drops out of the strand mean entirely:

```
    p = ½ + (κ−½)·d ,   d = f₊ − f₋        ⇒   d̂ = (p_obs − ½)/(κ − ½)      DIRECTLY OBSERVED
    |d| ≤ f₊ + f₋ = 1 − f_g               ⇒   f_g ≤ B = 1 − |d̂|            a HARD structural bound
```

`f_g` dropping out of the mean is *precisely why* the Schur complement is 0 — and it is the same fact that
makes `d` observable. The level set `d = (1−f_g)·τ` terminates where `|τ| = 1`, so the profile likelihood in
`f_g` has **finite support** `[0, B]`. τ_own = 0 says "λ could be anything"; the truth is "λ ≤ λ_B", with λ_B
sharply measured.

**How good is the bound?** (`scratchpad/ambig_1_census.py`, `ambig_3_bound_violation.py`)

| condition | mean slack `B − f_g` | frac slack < 0.05 | `|B − oracle|` |
|---|---|---|---|
| gdna300 ss0.99 nrna_none capON | **0.0081** | **97.0 %** | 0.0102 |
| gdna300 ss0.99 present capON | **0.0142** | **96.5 %** | 0.0138 |
| gdna300 ss0.99 present capOFF | 0.0923 | 63.9 % | 0.0878 |
| gdna100 ss0.50 present capON (κ=½) | 0.5517 | 18.2 % | 0.4041 |

The slack is *exactly* twice the minority strand's RNA (`B − f_g = 2·min(f₊,f₋)`), which is why it is tiny
under capture (the minority strand at an `exon+intron(x-strand)` region is the **intron** strand, unbaited) and
real without capture. **At κ = ½ the bound is worthless** — `d̂` is pure noise — and that is not a special case
to code around: `Var(d̂) ∝ 1/(κ−½)²` diverges there on its own.

## 3. ⭐ THE CEILING — what AMBIG is worth (`scratchpad/ambig_5_ceiling.py`)

Suite mwae with the AMBIG answer replaced offline:

| condition | shipped | AMBIG := **B** | AMBIG := oracle | AMBIG share of err |
|---|---|---|---|---|
| gdna300 ss0.99 none capON | 0.0354 | **0.0120** | 0.0099 | 72.1 % |
| gdna300 ss0.99 present capON | 0.0682 | **0.0378** | 0.0349 | 48.8 % |
| gdna100 ss0.99 present capON | 0.0887 | **0.0558** | 0.0515 | 41.9 % |
| gdna300 ss0.99 present capOFF | 0.0148 | 0.0123 | 0.0058 | 61.1 % |
| gdna100 ss0.50 present capON | 0.1654 | 0.2031 ✗ | 0.1160 | 29.9 % |

On stranded capture-ON the bound captures **85–91 %** of the whole AMBIG prize. This is the largest single
identified, correctable, **prior-free** error source left in pass-0.

## 4. ⭐ THE DESIGN, VALIDATED OFFLINE (`scratchpad/ambig_6_design.py`)

`B` as an *estimate* asserts `min(f₊,f₋) = 0` — a population claim, and it blows up on unstranded. The exact
identity closes it with a quantity the solver **already computes** (the per-strand RNA imputation `mo_p`/`mo_n`):

```
    f_g  =  1 − |d|  −  2·min(f₊, f₋)          ⇒    f̂_g = clip( 1 − |d̂| − 2·min(f₊^imp, f₋^imp), 0, 1 )
```

| condition | shipped | **cons** | B | oracle | AMBIG mwae: ship → cons (B) |
|---|---|---|---|---|---|
| gdna300 ss0.99 none capON | 0.0401 | **0.0235** | 0.0120 | 0.0099 | 0.1444 → **0.0652** (0.0103) |
| gdna300 ss0.99 present capON | 0.0713 | **0.0516** | 0.0378 | 0.0349 | 0.1758 → **0.0803** (0.0139) |
| gdna100 ss0.99 none capON | 0.0342 | **0.0207** | 0.0105 | 0.0074 | 0.1292 → **0.0639** (0.0149) |
| gdna300 ss0.99 present cap**OFF** | 0.0148 | **0.0115** | 0.0123 | 0.0058 | 0.1202 → **0.0762** (0.0863) |
| gdna_none ss0.99 present capON | 0.0238 | **0.0078** | 0.0092 | 0.0036 | 0.1027 → **0.0217** (0.0287) |
| gdna100 ss0.50 present capON | 0.1654 | 0.2031 ✗ | 0.2031 ✗ | 0.1160 | 0.2429 → 0.4281 ✗ |

It **halves AMBIG error on every stranded condition**, and on capOFF / `gdna_none` — exactly where the
minority strand genuinely carries RNA — it **beats the bare bound** (0.0762 vs 0.0863; 0.0217 vs 0.0287). So
the identity is doing real work, not just re-expressing the bound. The imputed minority is rough
(`minorIm` 0.029 vs oracle 0.004 on nrna_none), which is exactly why `cons` lands between `B` and the oracle.

## 5. What was A/B'd this session — both ~NEUTRAL, both instructive

| arm | what | r0 | r1 |
|---|---|---|---|
| `m8` (HEAD) | — | 0.0900 | 0.0700 |
| `dofgate` | gate the λ composition message across a single↔AMBIG DOF mismatch | 0.0902 (11/16/5) | 0.0699 (12/9/11) |
| `ambnd` | finite `v_own_lam` for AMBIG from the bound's log-odds variance | 0.0894 (4/2/26) | 0.0699 (3/2/27) |

Both **improve every stranded ss_0.99 condition and pay for it on unstranded capture-ON** — the same split
§2 predicts. `ambnd` is inert on 26/32 conditions by construction (`Var(λ_B) = ∞` at κ=½, recovered as a
LIMIT, no gate). Neither was landed: **~0.0006 aggregate is not worth the complexity, and §3 says why** — the
AMBIG prize is in the **MODE**, not in damping the messages. Both were reverted; the tree is exactly M8.

`ambnd`'s conceptual flaw, for the record: `Var(λ_B)` is the variance of the *upper endpoint*, used as the
variance of *λ*. The DL gap is measured against the node's own mode (the self-solve), so the two are different
quantities. The bound is **one-sided** information and a symmetric DL cannot consume it. Also measured:
messages violate the bound on 42–46 % of AMBIG nodes, but **those nodes have LOWER error than the compliant
ones** (0.0711 vs 0.1623) — so one-sided clipping is not the fix either.

## 6. ▶ (SUPERSEDED BY §0b — the estimator below is inadmissible; kept for the record)

Not as a post-hoc override. The work, in this project's order:

1. **Derive the precision.** The factor is a Gaussian on λ at mode `logit(f̂_g)`. Its variance is
   `Var(d̂)/(∂f_g/∂d)²` ⊕ `4·Var(min(f₊^imp, f₋^imp))`, with `Var(d̂) = p(1−p)/(N_eff·(κ−½)²)` (Poisson strand
   sampling ⊕ the `ω_r` overdispersion, `N_eff` as in `strand_evidence`). **It must diverge at κ = ½** — it
   does, from `(κ−½)²`, so the unstranded arm is protected as a LIMIT with no gate and no constant. Check it
   against the existing `disc` deadband, which encodes the same fact.
2. **MC-validate** it in `scripts/debug/message_variance_mc.py` (call it M9), the way M7/M8 were.
3. **Implement** in `simplex_logodds._solve_ambig_logodds` (the λ axis, alongside `lam_logprior`) — NOT in
   `node_init`, because it is a likelihood factor, not an init source. Then decide separately whether it also
   licenses a nonzero `τ_own` for AMBIG (which would switch DL on there — §5 says do that on its own evidence).
4. **A/B per-condition at refit=0 AND refit=1.** Watch `unstranded × capON` (must be bit-identical if the
   κ=½ limit is right), `verystrong`, and `gdna_none`.

**Expected prize** (§3/§4, offline): stranded capture-ON suite mwae 0.0342–0.0902 → 0.0207–0.0706, i.e. roughly
**−0.02 aggregate**, the largest single item left.

⚠ **Guard the owner's directive** (memory `ambig_needs_robust_training_subset`): any AMBIG stress test must
include ample single-strand nodes. And ⚠ **`min(f₊,f₋)` is a boundary statistic** — its imputed value is
positively biased (§4), which pushes `f̂_g` DOWN, i.e. toward the current (already too low) answer. That is the
conservative direction, but do not let the factor's precision hide it.

## 7. ⛔ DO NOT RE-RUN — this session's additions

| item | verdict |
|---|---|
| **DOF-mismatch gate on the λ message** | 0.0902 / 0.0699, 11 better / 16 worse (r0). Real stranded gain, cancelled by the unstranded loss. |
| **`v_own_lam` for AMBIG from `Var(λ_B)`** | 0.0894 / 0.0699, +0.0006 aggregate. Correct in kind (τ_own = 0 IS false) but the wrong object (endpoint variance ≠ λ variance) and far too small to justify. |
| **One-sided clipping of bound-violating messages** | Refuted by measurement: violating nodes have *lower* error (0.0711 vs 0.1623). |
| **`B` as the AMBIG estimator** | Right on stranded capture-ON (85–91 % of the prize) but catastrophic on unstranded (0.2031 vs 0.1654) and it asserts a silent minority strand. Use the §4 identity, never `B` alone. |

Everything in HANDOFF_7 §7 and HANDOFF_8 §7 remains in force.

## 8. Tools

`scratchpad/ambig_{1..6}_*.py`, all re-runnable, ~1 s/scenario against the cached substrate.
`ambig_2_mechanism.py` is a **bit-exact ψ replay for AMBIG nodes too** (`max|Δ| = 0`), with `-lam / -gdna /
-rna / -tilt / nomsg` channel ablation — it is how the λ message was identified as the dominant AMBIG
message defect (0.1444 → 0.0986 when removed).

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
python scratchpad/ambig_5_ceiling.py      # the prize
python scratchpad/ambig_6_design.py       # the validated estimator
```
