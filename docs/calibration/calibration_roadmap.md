# Calibration development roadmap — state, derivations, and phases

**Status:** living roadmap. Branch `calib-ambig-init-wip`. 1090/1090 tests green, ruff clean.
**Purpose:** a new session should be able to pick up EXACTLY here. This is the single entry point.
**Written:** 2026-07-15.

> **Read first:** [`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) (count-zero-information —
> authoritative) and [`prior_ramp_and_bp_roadmap.md`](prior_ramp_and_bp_roadmap.md) (the `0.5·λ` ramp, the
> L-sweep, the strand Fisher information). [`calibration_selfsolve_status.md`](calibration_selfsolve_status.md)
> is the earlier initiative — **§1.5's unified-solver plan is still the agreed target**; its §1.2/§1.3
> "boundaries are gagged" thesis is **refuted** (below).

---

## 0. Where we are in one paragraph

`ψ = strand + logP_g + logP_r` is the correct per-node objective and is what the code now implements — except
**`logP_r` is not written**, and that omission is not neutral: it silently selects the **Haldane** reference
(`p(ρ_r) ∝ 1/ρ_r`), which piles infinite prior mass at `ρ_r→0`, i.e. `f_g→1`. Combined with the gDNA arm's
absent reference, ψ is currently **Beta(0,0)** on the composition — improper at both vertices, a **vertex
amplifier**. Everything downstream — the message-precision saga, the `_DEFAULT_L` sensitivity, the flagship
collapse — is measured through that. **Phase 1 is to fix the reference. Nothing else is trustworthy until it
lands.**

---

## 1. The four things that are true (measured, not argued)

### 1.1 "Prior-free" does not exist — you only choose a reference

ψ is integrated on a **uniform-λ grid** (`λ = logit(f_g)`). Omitting a component's prior term does not give
neutrality; the grid supplies one. Writing `p_c(ρ_c) ∝ ρ_c^(c−1)` gives `ψ_λ = strand + c·(log f_g + log(1−f_g))`,
which in composition space is exactly **Beta(c, c)**:

| c | per-rate reference | composition | status |
|---|---|---|---|
| **0** | Haldane `1/ρ` | **Beta(0,0)** — improper at BOTH vertices | ⚠ **SHIPPED TODAY** |
| **½** | **Jeffreys `ρ^(−½)`** | **Beta(½,½)**, ∫ = π | ✅ **DERIVED — the answer** |
| 1 | flat in linear ρ | Beta(1,1) = uniform in f | proper, not derived |

Beta(0,0) hides while ψ is symmetric (the median sits at ½) and detonates the moment anything one-sided
arrives — because a message Gaussian on `log f_g` **saturates** as `f_g→1` (`log f_g → 0`), so it bounds `f_g`
only from *below*.

### 1.2 The derivation of `c = ½` (two independent routes, same answer)

* **Route A — Jeffreys for a Poisson rate.** `g_c ~ Poisson(ρ_c·E_c)` ⇒ `I(ρ_c) = 1/ρ_c` ⇒ `p(ρ_c) ∝ ρ_c^(−½)`.
  `logP` is stored as a density in **log**-rate, so `P(log ρ) = p(ρ)·ρ = ρ^(+½)` ⇒ contribution `+½·log ρ_c`
  `= +½·log f_c` + const.
* **Route B — Jeffreys for the composition.** A binomial proportion's Jeffreys prior is `Beta(½,½)`, which on
  the uniform-λ grid is `+½·log f_g + ½·log(1−f_g)`. **Identical.**

**The rule:** each component contributes its **fitted `logP_c` if we have one, else `+½·log f_c`**. Never
nothing. No Jacobian — the log-rate conversions cancel `log σ'(λ)` exactly, once per component.

**Verified** on node 1479 (ψ = its gDNA message alone; κ=0.5 ⇒ strand flat), posterior median vs grid width L:

| reference | L=4 | L=6 | L=10 | L=14 | L=20 | spread |
|---|---|---|---|---|---|---|
| c=0 Haldane (**shipped**) | 0.474 | 0.627 | **0.893** | 0.982 | 0.999 | **0.525** |
| **c=½ Jeffreys (derived)** | 0.445 | 0.461 | 0.457 | 0.441 | 0.416 | **0.045** ✅ |
| c=1 flat-linear | 0.441 | 0.437 | 0.422 | 0.406 | 0.383 | 0.058 |

The shipped arm reproduces the **observed 0.903** from the message alone. Mechanism, exactly: the message tail
flattens at `−½·p·m² = −0.619` nats ⇒ a plateau at 53.9 % of peak over *infinite* λ ⇒ plateau mass grows
linearly in L ⇒ the median walks up the grid. **Only c=½ is both derived and L-stable.**

### 1.3 Every σ²_imp experiment ever run is invalid

Under the shipped Haldane reference, `f_g` is **non-monotone in message precision** and peaks at ≈0.35:

| prec | 0.0 | 0.1 | **0.351** | 0.702 | 1.5 | 3.0 | 10.0 |
|---|---|---|---|---|---|---|---|
| **shipped** | 0.458 | 0.882 | **0.922** ← peak | 0.893 | 0.679 | 0.354 | 0.246 |
| c=½ | 0.458 | 0.488 | 0.487 | 0.457 | 0.387 | 0.313 | 0.254 |

`σ²_imp = 2.847` caps message precision at `1/2.847 = 0.351` — **dead centre of the peak.** So "×0.5 σ²_imp is
a free win" was walking down the far side of a plateau. **Every σ²_imp A/B — the ×0.5 finding, the "honest
count precision" refutation, the un-gag — was reading plateau position, not reliability.** All retracted.
Re-run them after Phase 1. (c=½ reduces but does not eliminate the non-monotonicity — residual ~15× smaller.)

### 1.4 THE FLAGSHIP, explained — capture destroys the prior's ability to adjudicate

Oracle rate-distribution widths (5–95 %), `quick_3to1_5mb`, nascent-absent:

| | **gDNA span** | **RNA span** |
|---|---|---|
| capture **OFF** | **0.49 decades** | 4.16 dec |
| capture **ON** | **3.54 decades** | 4.78 dec |

**Capture-off works because gDNA is 8× tighter than RNA**: `logP_g` pins the gDNA level and the remainder is
RNA. **Capture-on fails because enrichment blows P(ρ_g) out to 3.54 decades — as broad as RNA.** Neither prior
can adjudicate. ⇒ **`logP_r` alone will NOT fix capture.** Stratifying `logP_g` by signature class is what
collapses 3.54 decades back toward 0.49 *within* a stratum. ![priors](figures/priorfree_logPg_logPr.png)

---

## 2. The failing nodes (the unit of evidence — never use an aggregate)

`gdna_gdna300_ss_0.50_nrna_none_capture_on`, prior-free. Reproduce:
`scratchpad/trace_worst.py`, `scratchpad/census.py`.

Per-class census: **INTERGENIC (locked) = 0.0 % of error** (f_g = 1.000 exactly — trivial, needs nothing).
**EXON = 75.7 %** of error mass. The error *moves with the regime*: capture-off it is **INTRON = 79.6 %**.

| node | class | mass | ρ | true f_g | local | **final** | gDNA msg | RNA msg |
|---|---|---|---|---|---|---|---|---|
| 1479 | EXON | 70,662 | 40.9 | 0.382 | 0.510 | **0.903** | 0.265 | *silent* |
| 19 | EXON | 58,924 | 33.1 | 0.472 | 0.490 | **0.957** | 0.365 | *silent* |
| 2167 | AMBIG | 49,962 | 30.5 | 0.538 | 0.542 | **0.073** | 0.041 | implies 0.569 / 0.996 |
| 1623 | EXON | 51,320 | 33.4 | 0.458 | 0.510 | **0.090** | 0.073 | implies 0.326 |
| 1303 | EXON | 52,817 | 34.0 | 0.478 | 0.529 | **0.137** | 0.100 | implies 0.401 |
| 1299 | EX+IN | 46,270 | 35.8 | 0.450 | 0.490 | **0.067** | 0.046 | present |

Every one: **enriched exon (ρ 30–41) flanked by depleted boundaries (ρ 4–11)**, local precision **≈0.098**
(κ=0.5 ⇒ `I(f_g) ≡ 0` — the node knows *nothing*), messages carrying **88–96 %**. Nodes whose RNA message is
silent **float to 0.90–0.96**; nodes with one land near the gDNA message. That is §1.1, node by node.

**A hard limit, measured:** each exon's TRUE `ρ_g` (15.3–16.4) **exceeds its neighbours' TOTAL density**
(4.3–10.6) by **1.5–3.6×**. Under per-probe capture the exon interior is enriched and the boundary straddles
into the depleted intron. **No message from that boundary can ever give the exon the right answer** — so the
"boundary self-solve feeds the exon" thesis does **not** rescue the capture cell, and louder messages make it
worse. *(Caveat: this is a property of where the simulator puts probes; do not over-fit to it.)*

---

## 3. The phases

### Phase 1 — the reference. ⭠ **NEXT. Nothing else is trustworthy first.**
Implement the §1.2 rule: per component, `logP_c` if fitted else `+½·log f_c`. One derived constant (`½`,
Jeffreys), **zero tuned constants** — and it *retires* `_DEFAULT_L` as an accuracy knob.
**Files:** `simplex_logodds._local_loglik_logodds` (1-D) and `_solve_ambig_logodds` (2-D).
**Predicted on the six nodes:** 1479 `0.903 → ~0.45`; 19 `0.957 → ~0.47`; the L-sweep goes flat (spread
0.525 → ~0.045). **Falsified if** the L-sweep is not flat.
**⚠ OPEN — the 2-DOF reference is NOT settled.** Dirichlet(½,½,½) in `(λ,τ)` derives to
`½·log f_g + 2·log(1−f_g) − ½·(log f_pos + log f_neg)` `= ½·log f_g + log(1−f_g) − ½·log((1−τ²)/4)`
— note the **τ term REWARDS a strand-dominant tilt**; the naive "+½ per component" gets that sign backwards.
Implementing it left `test_gdna_sweep_zero_gdna_pin_and_monotone` at `f_g=0.564` vs its `<0.50` bound
(was 0.44), **undiagnosed** — the attempt was reverted to keep the tree green. **Diagnose before shipping.**

### Phase 2 — the prior-free solve as the substrate.
With the reference correct, run `node_sweep(gdna_prior=None)`. Already a first-class path
(`bp_solver.py:279`) and measured to **tie production** (0.0915 vs 0.0883) *under the broken reference* —
re-measure. This is what replaces the vacuous `f_g=1` substrate (which fits **P(D)** — the total-density
distribution — and is composition-vacuous by the project's own axiom).

### Phase 3 — fit `logP_g` **and** `logP_r` on it, stratified.
"Stratified" = **the DNA/RNA split** (both priors) **and** by signature class where §1.4 demands it.
**⚠ The naive `logP_r` is BLOCKED:** `fit(mass_global, eff_rna_global)` — the exact pass-0 mirror, zero new
code — **loses 14/16 conditions**, because pass-0 `r_hat = M` is the *same* vacuous P(D). Fit it on Phase 2's
solved belief, not on total density. Watch `gdna_rate_prior.py`'s `np.interp(left=logP[0])` clamp: inert today
(verified bit-identical) but becomes **load-bearing** the moment `logP_r` exists, because `ρ_r→0` projects far
below the fitted grid.

### Phase 4 — re-solve and iterate.
Prior-free → fit → solve → refit. Under the vacuous prior the refit *degraded* monotonically
(0.0840→0.1056 over iters 1→5); re-measure once Phases 1–3 land.

### Phase 5 — interdependence (the pie).
`node_geometry._type_belief` sets `var_g`/`var_p`/`var_n` **independently** from ONE 1-DOF solve, and
`bp_solver._scan` combines `fbg`/`fbp`/`fbn` as **three independent geometric means**; nothing enforces
`fbg+fbp+fbn=1`. On a 1-DOF node "var_g=0 but var_p=∞" is **self-contradictory** and the code can represent it.
The **final grid solve IS constrained** (`f_active = 1−f_g`) — the incoherence is in the **relay**. Fix: carry
ONE number (+ its precision) on the free axis for 1-DOF, TWO for 2-DOF; derive the rest.
(`calibration_bp_theory.md` §8 predicted exactly this.)

### Phase 6 — message precision. **LAST.**
Every prior verdict is void (§1.3). Re-derive only after 1–5. Candidate:
`calibration_selfsolve_status.md` §1.3's `Var(log ρ_g) = Var(log f_g)_src + 1/N_src + φ`.

---

## 4. Settled — do not redo

| | verdict |
|---|---|
| `ψ = strand + logP_g + logP_r`, no Jacobian | ✅ **correct** — the `(ρ_g,ρ_r)→(f_g,M)` Jacobian is `M/(E_g·E_r)`, composition-independent, so the two `−log ρ_c` conversions cancel `log σ'(λ)` exactly. **Do not "restore the Jacobian."** |
| the `0.5·λ` improper ramp | ✅ removed (Beta + Haldane + Jacobian summed to it). **But removing all three left Beta(0,0)** — the actual bug was only ever the Haldane. |
| `boundary_side_eff_length` = `E[min(ℓ,R)]/2` | ✅ fixed + test-guarded. Accuracy-neutral (messages were gagged). |
| `emit_locked` | ✅ permanent. quick all 0.1167→0.0980; zero-gDNA guard passes. |
| `_zero_spl` | ✅ deleted — was **inert** (the solve's `gdna_mass`/`rna_mass` outputs were dead). |
| spliced mass in ψ | ❌ **do not add.** At a junction mature RNA *splices*, so unspliced crossing = gDNA + nascent — genuinely disjoint from spliced. Spliced is directly-observed pure RNA, held out, already in the message via `node_densities`. |
| "boundaries are gagged by σ²_imp" (`selfsolve_status` §1.2/§1.3) | ❌ **refuted** — they were **muted** (`emit_locked`), and §1.3's un-gag regresses. |
| the `−ε·(1−f_g)` nascent≈0 tie-break | ❌ dead — "not robust at real κ≈0.499" (`node_prior_design.md` §7). |
| `all` as a metric | ❌ **gDNA-biased** (suite true mass-wt f_g = 0.554). Use zero-gDNA + per-class. |
| strand information | `I(f_g) ∝ n(½−κ)²` — **exactly 0 at κ=½**, any count, any overdispersion. Structural, not fixable. |

## 5. Tooling

Caches (24 + 16 conds, ~9 s/cond): `~/Downloads/rigel_runs/{ambig_dense_10mb,quick_3to1_5mb}/_calib_cache`.
Scratch (`.../scratchpad/`): `trace_worst.py` (the six nodes → their messages), `census.py` (per-node-class
error), `plot_priors.py` (§1.4), `prior_ab.py` (in-process reference A/B via `_STRAND_PRIOR`/`logprior`
monkeypatch — `ψ_prior = a·log f_g + (a−h)·log(1−f_g)`; today = `(0,0)`, derived = `(0.5,0.5)`).
