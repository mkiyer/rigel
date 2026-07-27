# Session handoff — pass-0 is CONVERGED and CLEAN. Next: Phase 2, the gDNA hyperprior.

**This is the LIVE handoff. START HERE.** Date: 2026-07-26. Branch `calib-ambig-init-wip`, HEAD `241a87dd`.
Supersedes `SESSION_2026_07_26_HANDOFF_14.md` (still the P1d/P4 record).

Suite **0.0841 (refit=0) / 0.0679 (refit=1)**. All gates green.

---

## 0. Where the project is, in five lines

* **Pass-0 is done as a prior-free solver.** Every item on `PASS0_FINISH_PLAN.md` is landed, refuted or
  deliberately deferred. The solver is **BP-legal** and `bp_solver.py` reads **no environment variables** —
  there is exactly one production path.
* **The code has had a full cleanup pass**: −892 lines, 12 commits, every one verified **32/32
  byte-identical** on the benchmark at both refits.
* **Two DEBTS ship knowingly** — `ω_graft` and P1e. Both are banner-documented in six places each. Read §6.
* **The next task is Phase 2**, and its first step is **P2a**, ~6 lines that were written and measured once
  and then blocked by two premises that are **both now false**.
* **⚠ P2a's original snippet will not compile as written** — it references `_lam_lo`, which the cleanup
  deleted. One-line fix, given in §2.2. Do not lose 20 minutes to this.

## 1. Read these, in this order

1. **`ROADMAP.md`** — the entry point. §4 (why Phase 2 is unblocked) and §6 (the ordered path).
2. **THIS FILE** — the Phase-2 plan, the run book, and the debts.
3. **`variance_ledger.md`** — every variance term, the no-double-count proof, and §6 (P1e's debt).
4. **`PASS0_FINISH_PLAN.md`** — P0–P6 with live status, plus the refuted/withdrawn items.
5. `CALIBRATION_ARCHITECTURE.md` — count-zero-information, the three information sources. Still authoritative.
6. `SESSION_2026_07_26_HANDOFF_14.md` — P1d, P4, and this session's do-not-re-run table.
7. `SESSION_2026_07_25_HANDOFF_6.md` §3 — the ORIGINAL P2a measurement. ⚠ Its premises have changed; read
   §2.1 below before trusting its conclusion.

**Never read `docs/calibration/archive/`** for guidance. Code docstrings now cite it with an explicit
`archive/` path — that is provenance (why the code is shaped as it is), not instruction.

`docs/calibration/boundary_work_notes.md` and `docs/accumulator/dag_design.md` are the **owner's own files**.
Do not overwrite them.

---

## 2. ▶ P2a — THE NEXT TASK

### 2.1 Why it is unblocked (both original objections are now false)

`SESSION_2026_07_25_HANDOFF_6.md` §3 implemented and measured a ~6-line step that feeds the hyperprior's own
λ-curvature into the DerSimonian–Laird null `v_own`. It improved stranded (0.0376 → 0.0333), verystrong
(0.1292 → 0.1196) and capture-OFF (0.0354 → 0.0168), and **regressed unstranded × capON 0.1702 → 0.2177**.
It was reverted, and Phase 2 was declared blocked on two grounds. **Both have since evaporated:**

**(a) "the refit regresses unstranded capture-ON."** Re-measured at HEAD (mass-weighted, 32 conditions):

| stratum | refit 0 | refit 1 | Δ | better/worse |
|---|---|---|---|---|
| ALL 32 | 0.0841 | **0.0679** | −0.0162 | — |
| **unstranded × capON** | 0.1523 | **0.1275** | **−0.0248** | 4 / 4 |
| capture OFF | 0.0475 | **0.0319** | −0.0156 | 14 / 0 |
| gdna_none | 0.0941 | **0.0029** | −0.0912 | 9 / 0 |

The refit improves **every** stratum. It already did so before P1d landed, so the intervening work fixed it
and the note outlived its own resolution.

**(b) "AMBIG's over-confidence poisons the fit."** It cannot: `_fit_gdna_hyperprior` already selects
**REGION nodes that are single-strand or structural-gDNA — no AMBIG, no boundaries** ("non-circular",
`calibrate.py`). That substrate is at `z2|Q1` **2.26**; the 42.04 sits entirely in the excluded half.

### 2.2 The step, and the trap

Insert in `bp_solver._unified_solve`, where `v_own_g` / `v_own_r` / `v_own_lam` are built (~line 460):

```python
    # the hyperprior's own λ-curvature, by the SAME law the intron factory uses — no new law, no constant.
    # `build_node_init` folds the prior into the self-solve's MODE but not into τ_λ, so on a refit the own
    # belief was being judged at a precision that ignored half of what formed it.
    tau_prior = density_factor_precision(global_lp, _lam_lo)   # None at pass-0 ⇒ identically inert
    tau_dl = tau_own if tau_prior is None else tau_own + np.asarray(tau_prior, np.float64)
    v_own_g, v_own_r = own_composition_logvar(_ni.f_g, tau_dl, _struct)
    v_own_lam = np.where(_struct, 0.0, np.where(tau_dl > _EPS, 1.0 / np.maximum(tau_dl, _EPS), np.inf))
```

> ### ⚠ TRAP — `_lam_lo` NO LONGER EXISTS
> The 2026-07-26 cleanup collapsed `_lam_lo, _fg_lo = _logodds_grid(...)` to `_, solve_grid = ...` because
> `_lam_lo` was then unused (`bp_solver.py:192`). P2a needs it back. **The fix is one line:**
> ```python
> _lam_lo, solve_grid = _logodds_grid(int(n_grid), float(logodds_window))
> ```
> `_logodds_grid` returns `(lam, fg)` — `lam` is the λ lattice `density_factor_precision` wants as its
> second argument (`density_deconv.density_factor_precision(lam_logprior, lam_grid)`).

**The invariant that makes it safe** (do not violate it): add `tau_dl` to the **DL null only** — never to the
SEND precision (`prec_*`). A node knowing more about itself is grounds for trusting a *contradicting* message
**less**, not for shouting louder. So the step can only ever **damp**. It is also **identically inert at
pass-0** (`global_lp` is None ⇒ `tau_prior` is None ⇒ `tau_dl is tau_own`), which means **refit=0 must come
back 32/32 bit-identical** — that is your first check, and if it does not, the wiring is wrong.

### 2.3 How to judge it

Gate on **refit=1**, and on **both** axes (see §4):
* suite mwae and the named strata, per-condition;
* **accuracy on the fit substrate** — the number the prior is actually built from;
* `z2` on a **held-fixed node set** and confidently-wrong mass.

Expect stranded / verystrong / capture-OFF to improve. **Watch unstranded × capON specifically** — that is
where it regressed originally. If it still regresses, that is a *result*, not a failure: it localises the
remaining problem, and §3 is then the work.

---

## 3. THE HYPERPRIOR ROADMAP (P2a → ship candidate)

### Step 0 — orient (do this first, ~15 min)

Read `calibrate.py::_fit_gdna_hyperprior` and `npmle.py::DensityNPMLE`. The whole Phase-2 surface is:

| thing | where | note |
|---|---|---|
| the fit | `calibrate._fit_gdna_hyperprior` | selects the substrate, calls `DensityNPMLE.fit` |
| the refit loop | `calibrate.calibrate`, ~line 410 | `for it in range(config.calib_refit_iters)` |
| **the reset** | same loop | `belief = _init_belief()` **before** `_sweep(prior)` — Phase 2 already re-solves FROM SCRATCH |
| the prior in ψ | `bp_solver`, `global_lp` → `_local_solve` | the composition arm |
| config | `CalibrationConfig` | `calib_refit_iters=1`, `npmle_bandwidth=0.15`, `gdna_prior_additive=False`, `background_floor=True` |

**The substrate, exactly as the code defines it** (`sel = live & isr & (single | gonly)`, then
`sel &= ~intergenic` because `background_floor=True` always fires): **single-strand non-intergenic REGION
nodes — 26,016 nodes, 64.5 % of node mass.** The `| gonly` clause is **dead** in the default config.

### Step 1 — P2a (§2). The cheapest high-value experiment on the list.

### Step 2 — ⭐ P-sub, the SUBSTRATE RESEARCH PHASE (owner-directed)

> *"We just need a dedicated research phase for the hyperprior to see what substrate should be used."*

This replaces the withdrawn "solvable gate" (see `PASS0_FINISH_PLAN.md` §P-solv — do **not** build a binary
flag). It is an empirical sweep. **Score every candidate against the ORACLE-fitted prior**, not against mwae:

```python
# the scoring harness: fit the prior from (a) the arm's beliefs and (b) the ORACLE f_g on the same nodes,
# render both P(log ρ_g) on a common grid, and report total-variation distance.
```

Anchors already measured, so you know the scale:
* today's achievable prior is already **~83 %** of the way from "a different library's prior" to perfect
  (TV 0.1375 vs a between-condition TV of 0.830);
* it is **33 % too BROAD** (sd 0.812 vs 0.610 decades) and has 1.9× too many cells — **smeared, not
  displaced** (median off by 0.09 decades);
* **precision-weighting is load-bearing**: zeroing the declared widths costs **0.0234 TV**, ~10× the entire
  P1d+P4 accuracy regression, and it helps most exactly where the prior is worst (gDNA-free libraries,
  +0.21 TV);
* the prior is badly wrong on the **`none_*` (gDNA-free) family** — TV up to 0.508, median log₁₀ρ_g −2.30
  against an oracle −4.17 — and on `verystrong` (0.273). Everything `gdna100`/`gdna300` × `ss_0.99` is ≤0.09.

Candidates to sweep: all nodes / single-strand only (today) / include AMBIG / include boundaries / exclude
`τ_own = 0` nodes / precision-weighted vs flat / `gdna_prior_additive` True vs False / `npmle_bandwidth`.

⚠ **The known counter-argument, and it is measured**: excluding nodes can *hurt*. `solve_gate_design.md`
records a stronger version (skip unidentified nodes, defer to the prior) as derived, implemented and
**empirically refuted** — +0.010 (r0) / +0.025 (r1). The prior resolves an imperfectly-solved node better
than a deferred `f_g = 1`. And a node reporting a genuine **enriched mode** may be exactly the one an
exclusion rule drops.

### Step 3 — the re-solve, and AMBIG

AMBIG is `z2|Q1` **19.8** after P1e (was 183 → 92.1 → 64.4 → 19.8). It is **excluded from the fit**, so it is
the re-solve's problem, not the fit's — and it is the population the hyperprior exists to rescue
(`strand_likelihood_constrains_tilt_not_fg`: the strand likelihood constrains only the tilt, never `f_g`).
Judge Phase 2 partly on what it does here.

### Step 4 — ship candidate

Only after: the refit improves every stratum, the substrate is chosen on evidence, and the two DEBTS in §6
are either retired or explicitly accepted with the owner.

---

## 4. RUN BOOK — everything you need

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1          # REQUIRED for A/B determinism
cd /Users/mkiyer/proj/rigel
```

### 4.1 Gates — all must stay green

```bash
python -m pytest tests/calibration tests/native -q   # fast loop (~6 s)
python -m pytest tests/ -q                           # 1193 pass, 1 xfail, 3 xpass
ruff check src/ tests/ scripts/                      # clean
python scripts/debug/message_variance_mc.py          # 0 failures over M1–M11
python -m pytest tests/ --update-golden -q           # goldens — LAST, after the solver is final
```

### 4.2 The A/B — 32 conditions, cached, ~1 s/scenario

```bash
P0_REFIT=0 python scripts/debug/pass0_oracle_bench.py --arm NAME_r0
P0_REFIT=1 python scripts/debug/pass0_oracle_bench.py --arm NAME_r1
python scratchpad/ab.py NEW BASE [extra arms...]     # strata + per-condition, both refits
```

Rows accumulate in `/tmp/pass0_oracle_bench.tsv`. **⚠ `/tmp` is volatile. If it has been cleared, the
baselines are gone — re-establish HEAD as `clean0_r0`/`clean0_r1` before comparing anything.**
`scratchpad/ab.py` reports the strata the owner asks for: `stranded ss_0.99`, `unstranded × capON`,
`verystrong`, `gdna_none`, `low gDNA`, `nrna_none`, capture OFF/ON.

Suite: `~/Downloads/rigel_runs/ambig_dense_10mb` (32 conditions, cached in `_selfsolve_cache`).

### 4.3 The TRUST view

```bash
python scripts/debug/pass0_error_table.py --refit 0   # writes /tmp/pass0_state.npz
```

Columns: `CWRONG` (confidently-wrong reads), `errQ1conf`, **`z2`**, `%nodeQ1`.

> ### ⚠⚠ ALWAYS COMPARE `z2` ON A HELD-FIXED NODE SET
> A self-quartile or per-arm fixed-threshold view **re-selects the population** and hides the very effect you
> are measuring — a change that widens variances simply drops its own worst nodes out of the "confident"
> set. **This error was made three times in one session.** The correct form: freeze the node set to the
> BASELINE arm's confident quartile, then let each arm bring its own error and its own variance.
> `scratchpad/v_p4_paired_z2.py` does it correctly; copy that, not `pass0_error_table.py`'s own view.

⚠ `z2` is also **dimensionally mixed** (a linear squared error over a log variance — the house convention).
Trust the *relative* movement between arms on a fixed set; do not read an absolute `z2` as a calibration
ratio.

### 4.4 The metric that decides Phase 2

Track **both**, always (owner: *"if we push towards a combination of ACCURACY with honest precision, we will
prevail"*):

* **ACCURACY on the fit substrate** — `REGION & single-strand & live`. Suite-wide mwae is a *proxy*, not a
  substitute: the substrate is a subset and the two can move in opposite directions.
* **Honest precision** — it has three live consumers: pass-0's own fused mode, the **precision-weighted**
  prior fit (`var_gdna` weights the NPMLE), and a possible iterative re-solve.

### 4.5 Diagnostics

**There are no ablation flags.** All six were removed once their A/Bs settled; `bp_solver.py` reads no
environment variables. To ablate something, reintroduce the branch locally — **do not re-add a flag**.

Inert `_capture` buckets (populated only when `_debug` is passed): `_pin` (pre-pin per-message state), `_dl`
(per-message DL gaps), `_lvl` (level provenance), `_glv` (P1d's fit), `_uni` / `_uni_static`.

`scratchpad/` (untracked, **keep it**): `ab.py` · `p1_worked.py` (**one node in fragments and base pairs —
start here when something looks wrong**) · `v_p4_paired_z2.py` (the held-fixed z2) · `glv_explain.py` /
`glv_tail.py` (P1d's fit) · `p1e_*.py` (P1e's derivation) · `verify_clean.sh` (the behaviour-identity gate).

### 4.6 Behaviour-identity gate for any refactor

```bash
bash scratchpad/verify_clean.sh ARMNAME   # ruff + pytest + 32/32 vs clean0, both refits
```

---

## 5. Invariants — preserve these

* the **`√2·σ_own`** DL pin-safety inequality (`mismatch_deflate`);
* `Σ_c ρ_c·E_c = M` with `_pin_v`'s **partial-claim semantics** — a component the message does not *supply*
  is filled from the node's OWN density and does not move ("supplied" is about PRECISION, never the value);
* the composition is ONE dof (λ); the tilt θ is SEPARATE; **each is fused by its own precision**;
* `N` enters only as power — M11's `k` is literally an effective fragment count;
* **no magic numbers** — structural presence tests, derived quantities, exact limits, or parameters *fitted
  from the data* by method of moments (κ, both strand overdispersions, and `ω_graft` are the precedents);
* **pass-0 must be WEAK and CORRECTABLE** — prefer the under-confident option when unsure;
* **before adding any variance term**, name the physical event it prices and show which existing term would
  otherwise have absorbed it (`variance_ledger.md` §2), then update the ledger;
* **`_relay` and `_transport` are deliberate twins — DO NOT MERGE** (measured 15.7×; see the note at both
  definitions). A change to one must be mirrored in the other.

## 6. ⚠⚠ THE TWO DEBTS — do not let these be forgotten

Both are banner-documented in `CLAUDE.md`, `ROADMAP.md`, `PASS0_FINISH_PLAN.md`, `variance_ledger.md` and the
code. Summarised here because Phase 2 inherits them.

**`ω_graft` (P1d)** — compensates for a **failure in the structural representation**: the region/boundary map
has no TSS/TES, so the solver cannot tell a splice junction from a transcript terminus. One library-wide
fitted scalar stands in for a quantity measured to split **≥30×** on that missing bit (1.7–1.9 at termini vs
0.04–0.06 at junction-only seams; 20.8 % of edges carry 71.7 % of the error). Expected to be **fragile on
real data** (~200 exons, 30–50 % SE, a Poisson-by-construction suite). **Must be re-derived per structural
class when TSS/TES enter the map (P1g).**

**P1e (the conservation surprise)** — **prices a BIAS as a VARIANCE** on most of its firing mass (bias share
53–77 % on graft × one-component, 98.9–99.2 % at intergenic destinations), and **90–100 % of that damping
lands on solvable nodes**, so it is not inert. A variance cannot move a mode toward truth and it never
shrinks. Landed anyway (the only change to improve accuracy AND honest precision together), scoped to the
licensed `δ < 0`. **When the bias strata are diagnosed, this term must SHRINK. If it does not, the model has
not improved.**

## 7. ⛔ DO NOT RE-RUN

This session's table is `SESSION_2026_07_26_HANDOFF_14.md` §3, which also **indexes** the earlier sessions'
tables (they live in their own now-`⛔ SUPERSEDED`-stamped handoffs). Highlights that bear on Phase 2:

| item | verdict |
|---|---|
| the **solvable gate** | ⛔ withdrawn as a concept — pass-0 now solves every node, just inaccurately. `solve_gate_design.md` records a stronger version as derived, implemented and **empirically refuted** (+0.010 r0 / +0.025 r1) |
| **holding out** the nodes not used to fit the prior | measured **unnecessary** — the refit already resets, and holding the excluded half out moves the fit substrate by **±0.0006** |
| the `_relay`/`_transport` **merge** | measured **15.7× slower**; rejected |
| a **belief-free constant frame** for the reframe | +0.0003 suite / +0.0010 fit-substrate, neutral-to-worse on held-fixed z2. ⚠ A sub-agent reported it as a large trust win; that did **not** reproduce |
| P1d **per-edge** (instead of the pooled scalar) | worse on every axis; `d²` from one pair is a χ²₁ (CV = √2) |
| the P1e **rank-1** attribution | MC-refuted — over-damps λ 5× on a pure scale error |

## 8. Vocabulary (owner's — use it)

RNA is **one species**. "Mature" and "nascent" are bookkeeping; only **SPLICED vs UNSPLICED** is observable.
A boundary can be an exon↔exon boundary that is *also* a splice junction, with RNA contiguous across it while
other RNA splices in or out. Both at once.
