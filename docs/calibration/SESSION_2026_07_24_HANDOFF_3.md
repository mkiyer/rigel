> ## ⛔ SUPERSEDED — THIS IS NOT THE LIVE HANDOFF
> The live handoff is **`SESSION_2026_07_26_HANDOFF_14.md`**; the entry point is **`ROADMAP.md`**.
> This file is kept as HISTORY — what was tried, measured and refuted. Its numbers describe the
> code as it was on its own date and many are now superseded. **Do not act on it without checking
> the live handoff first.** Its own DO-NOT-RE-RUN findings (if it has any) remain HERE and still
> stand — HANDOFF_14 §3 indexes which files carry them.

# Session handoff — the MESSAGE variance model (next session begins implementation)

*(this file's own original header — superseded)* — it supersedes
`SESSION_2026_07_24_HANDOFF.md` (session 1) and `SESSION_2026_07_24_HANDOFF_2.md` (session 2), which are now
historical. Date: 2026-07-24. Do NOT read `docs/calibration/archive/`.

---

## 1. Status in one line

The variance **FOUNDATION is settled and landed** (approach E: a node's local composition precision is the
single Schur-marginal scalar `τ_λ`; the strand-gate bug fix is committed + A/B-validated). The **NEXT task is the
MESSAGE variance model** — build the message-propagation precision on the `τ_λ` foundation, absorbing the
deferred composition/sampling (`1/n`) separation. Pass-0 is NOT ready to ship until this lands and wins the A/B.

## 2. Where we are — the full arc (all committed on branch `calib-ambig-init-wip`)

1. **Converged to ONE solver** (`b224bcd8`). Deleted the legacy `_scan` path + all its flags
   (`RIGEL_B1B/N4A/N4B/E2`, the `_UNIFIED` gate) and helpers; the composition/enrichment-ratio unified solver is
   the sole `node_sweep`. `bp_solver.py` 1871 → 726 lines. Extracted the per-node self-solve into
   **`node_init.py`**. Behavior-preserving (byte-identical to the pre-refactor unified path across all 32
   `ambig_dense_10mb` scenarios). Goldens regenerated to the unified default.
2. **Extracted `density_deconv.py`** (`05b70516`). The generic density deconvolution (counts + gDNA prior →
   gDNA/RNA, NB precision); the **intron factory is now the special case** (`fit_intron_background`, the gDNA
   prior = the intergenic node distribution). Byte-identical.
3. **Variance FOUNDATION settled** (`c6df8c50`, `999b753c`). The honest local composition precision is a
   **single Schur-marginal scalar `τ_λ`** (approach E) — NOT a diagonal `(τ_λ,τ_θ)`, because the strand
   Beta-Binomial is **rank-1**. Derived by a 5-approach derivation workflow, numerically re-verified
   (`scratchpad/verify_foundation.py`, all claims PASS), and independently critiqued — **both analyses converge
   on approach E + Option B**. Phase 1 (the strand-gate bug fix) is landed + validated.

## 3. The full roadmap (ROADMAP.md §6)

```
   converge to one solver  →  variance FOUNDATION (τ_λ)  →  MESSAGE variance model  →  gDNA hyperprior refit  →  re-solve  →  ship
        [DONE]                     [DONE]                    [NEXT — this session]        [after]              [after]
```

## 4. The settled FOUNDATION model (what "done" means — do not relitigate)

The strand likelihood depends on `(λ=logit f_g, θ=arcsin τ)` ONLY through `p = ½+(κ−½)(1−f_g)sinθ`, so its
Fisher information is rank-1. The honest gDNA-level precision is the Schur-marginal scalar:
```
    τ_λ = τ_density  +  c·a² · 1[single-strand]      a = ∂p/∂λ = −(κ−½)sinθ·f_g(1−f_g),  c = N_eff/(p(1−p))
```
* **single-strand (1-DOF):** tilt structurally locked ⇒ `τ_λ = τ_density + c·a²` (strand pins f_g);
* **AMBIG (2-DOF):** tilt free ⇒ strand cancels out of f_g (Schur ⇒ 0) ⇒ `τ_λ = τ_density` (density-only).
Verified: Schur = `I_density` to 1e-13; `(λ,t=f₊−f₋)` block-diagonal exactly; the AMBIG-strand-only λ-precision
is N-invariant. Phase 1 gated the strand λ-term to single-strand nodes (`node_init.build_node_init`); the fix
is a bounded phantom removal that **improves the stranded arm** (A/B below) with no regression.

## 5. Exact code state

* Commits: `999b753c`, `c6df8c50`, `05b70516`, `b224bcd8` (+ base `9d14bace`) on `calib-ambig-init-wip`.
  Working tree also has uncommitted doc drafts (`variance_foundation_{proposal,plan,critique}.md`, this handoff,
  `variance_model_concepts.md`) — commit or keep as you prefer.
* Gates: `pytest tests/calibration tests/native` = **371 pass, 3 xfail, 1 xpass**; `ruff check src/ tests/`
  clean; full suite green; goldens on the unified default.
* Benchmark (pass-0 vs oracle, unified default, `scripts/debug/pass0_oracle_bench.py`, `OMP_NUM_THREADS=1`):
  refit=0 aggregate mwae **0.1491**; refit=1 (ship path) **0.0889** (stranded arm clean, ~0.043 capON). The
  unified solver LOSES the A/B vs legacy-with-factory until the message variance model lands — **expected, WIP,
  nothing ships.**

## 6. THE NEXT TASK — the message variance model (ordered implementation)

**Context.** The unified solver transports COMPOSITIONS via enrichment ratios (the *mode* is correct — exon
message f_g 0.682 vs oracle 0.677). But the message-propagation PRECISION is still the OLD density-uniformity
proxy (`σ²_transfer` = the NPMLE-projection cliff-height, `var_proj[dst]+(μ_proj[dst]−μ_proj[src])²`), which is
WRONG for a composition transport under hybrid capture. Replace it. Substrate: `variance_model_handoff.md` (the
prior derivation attempt — its per-component framework stands; §7 lists carry-forward-vs-discard). Now built on
the `τ_λ` foundation.

### Step 1 — DERIVE (math + Monte-Carlo, BEFORE any solver code)
Use `scripts/debug/message_precision_mc.py` (and `verify_foundation.py` as the MC template). Derive + MC-validate:
* **(1a) The per-component density message variance** `Var(log ρ_c) = Σ_k (ρ_k/ρ_c)²·v_k + σ²_transfer`
  (`variance_model_handoff.md` §3). The GRAFT (boundary→exon) is a **SUM** (`ρ_R = ρ_ν + ρ_μ`; share-weighted,
  the item-E rule); the PEEL (exon→boundary) is a **DIFFERENCE** (`ρ_ν = ρ_R/r − ρ_μ`; u-weighted, weights ≥ 1).
  **The `1/n` sampling enters HERE** — the source's per-component Poisson noise `v_k = 1/n_c + Var(log f_c^src)`
  — which is why the composition/sampling separation belongs to this task, not the foundation. **MC-validate
  including the pure-gDNA anchor limit** (where the old ratio-`k` form was singular).
* **(1b) The transfer variance** `σ²_transfer = Var(log r)` (`variance_model_handoff.md` §4), `r` = the
  enrichment ratio. **Direction-dependent:** ~0 on the GRAFT (enrichment cancels, the T1 result), load-bearing
  on the PEEL (carried ~92 % of the peel's variance in the earlier MC). Covers the **relay AND the combine**.
  Retire the density-uniformity proxy. `enrichment_frame.composition_logvar` derives `Var(log ρ_tot)`.
* Carry forward (handoff §7): the per-component framework; T1 (enrichment cancels on the graft); the item-E
  share-weighting; the peel's difference-variance + u-weighting; the MC harness; the **Poisson (ω=0)** decision.
  Discard: the ratio-`k` parameterization; the current `σ²_transfer` proxy.

### Step 2 — IMPLEMENT (pure, tested arithmetic; the closure shrinks)
Put the derived laws in the pure layer (extend `enrichment_frame.py`; rename it `message_arithmetic.py` once it
owns all message math), mirroring `node_init.py`/`density_deconv.py`. Then wire into `bp_solver._unified_solve`.
This is where the DEFERRED FOUNDATION plumbing lands:
* **Split `own_precision`'s two roles** (`bp_solver.py` uses one `pg_own` for both today, lines ~377 transport /
  ~418 fusion): the **fusion weight** `w_c = 1/Var(log f_c)` (composition only) vs the **transport seed**
  `1/(Var(log f_c) + 1/n)`. Land byte-identically first (verbatim `own_composition_logvar + own_precision`, NOT
  the ULP-different `1/(1/τ+1/n)`), THEN flip the fusion to `w_c`.
* **`struct_lock` HARD OVERRIDE in `_fuse`** — a composition-certain node ADOPTS its own belief (skip soft
  fusion); NEVER an `∞` weight (`1/v_log = 1/0 = ∞ ⇒ (∞·a+p·b)/∞ = NaN` cascading through interior anchors).
  REQUIRED before the fusion flip. Emit precision stays `n`. Add an interior-anchor **no-nan** test.

### Step 3 — VALIDATE + LOOP (the standing methodology)
Full `ambig_dense_10mb` suite → **worst scenario by error MASS** → dissect its worst nodes
(`scripts/debug/pass0_node_dissect.py`, the exact-replay ψ ablation) → root cause → fix → repeat. Compare
**pass-0 vs oracle**, per-condition, at **refit=0 AND refit=1** (`P0_REFIT=1`). **PASS iff** aggregate + EVERY
`{stranded,unstranded}×{capON,capOFF}` condition does not regress past the noise floor `ε` — call out
**unstranded-capON** explicitly (the known landmine, memory `pass0_intron_factory_precision_resurrection`). "More
correct" NEVER overrides a measured per-condition regression. Target: unified ≥ legacy-with-factory (0.0949), no
stranded regression. Regenerate goldens last; then hand to the gDNA-hyperprior refit (ROADMAP §4).

## 7. INVARIANTS — must preserve (plan v4 §2)

| principle | rule |
|---|---|
| count-zero-info | `N` enters ONLY as power (`τ_λ` Fisher, or `1/n` sampling), never a composition vote. |
| parameter space | variances in **log-odds** `λ` / `Var(log f_c)`, **NEVER** on simplex fractions. |
| local precision | ONE Schur scalar `τ_λ` per node; a diagonal `(τ_λ,τ_θ)` is prohibited (rank-1 strand). |
| composition ⟂ sampling | fusion weight = `τ_λ`; transport = `τ_λ ⊕ 1/n`. |
| anchor fusion | `struct_lock` = HARD OVERRIDE, never an `∞` weight → nan. |

## 8. Corrected census — do NOT delete these as "dead" (the draft plan was wrong)

* `belief.var_gdna` — **LIVE** (weights `_fit_gdna_hyperprior`, calibrate.py:174). Keep.
* `rna_pos/neg_frac_var` — **LIVE** (`node_init._rna` liveness gate; a numeric no-op today but a real reference).
  Re-express the gate as `free_s & n>0 & rho>_EPS` before ever deleting the field.
* `strand_likelihood.py` — production-dead but the **test oracle** for `_mixture_strand_loglik`. Keep.
* `lam_var`/`theta_var` — dead in production, but **`message_precision_mc.py` references `lam_var`** (this task's
  own harness). Keep.
* `var_pos`/`var_neg` — no ship reader, but asserted in ~5 tests. Re-point, don't silently drop.
* `struct_lock = locked & is_region` — region-only by design (a G1 seam must NOT be struct_lock — phantom emitter).

## 9. Methodology + gotchas

* **Derive → MC-validate → implement → per-condition A/B → loop.** MC-validate EVERY law before solver code
  (this session's foundation win came from exactly this; `verify_foundation.py` is the template).
* **No magic numbers** — pause and discuss before any new constant. Keep the module/constant count small.
* **Counts are Poisson** (ω=0). The synthetic suite is Poisson by construction (memory
  `synthetic_suite_is_poisson_omega_zero`) — it CANNOT validate an overdispersion term.
* **The A/B is deterministic only with `OMP_NUM_THREADS=1`** (memory `calibrate_cross_process_nondeterminism`).
* **Cached suite:** `~/Downloads/rigel_runs/ambig_dense_10mb` (`_selfsolve_cache`), ~1 s/scenario.
* **Multi-agent methodology worked well** this session — a derivation-workflow of candidate approaches +
  adversarial critique + an independent verification pass. Reuse it for the message-variance derivation.

## 10. Key files + tools

* Model: `node_init.py` (`build_node_init`, `strand_evidence`, `own_composition_logvar`, `own_precision`),
  `density_deconv.py`, `simplex_logodds.py` (the ψ mode engine), `enrichment_frame.py` (the pure transport math).
* Solver: `bp_solver.py` — `node_sweep` → `_unified_solve` (`_relay`, `_transport`, `_fuse`, `_damp`, `_pin_v`,
  `node_total_density`).
* Docs: `variance_foundation_proposal.md` (the derivation), `variance_foundation_plan.md` v4 (invariants +
  deferred spec), `variance_foundation_critique.md`, `variance_model_handoff.md` (§3-4 the message-model
  substrate), `unified_solver_design.md`, `CALIBRATION_ARCHITECTURE.md`.
* Tools (`scripts/debug/`): `pass0_oracle_bench.py` (the A/B, `--arm`, `P0_REFIT` env, `--report`),
  `pass0_node_dissect.py` (ψ ablation), `pass0_error_concentration.py`, `message_precision_mc.py` (the MC harness),
  `verify_foundation.py` (the foundation verification, an MC template).

## 11. ▶ KICKOFF PROMPT (copy-paste)

> We are developing Rigel calibration — pass-0. The variance FOUNDATION is settled and landed (a node's local
> composition precision is the single Schur-marginal scalar `τ_λ`; the strand-gate bug fix is committed). Read
> `docs/calibration/ROADMAP.md`, then `docs/calibration/SESSION_2026_07_24_HANDOFF_3.md`,
> `docs/calibration/variance_foundation_proposal.md`, and `docs/calibration/variance_model_handoff.md`. Do not
> read `docs/calibration/archive/`.
>
> **The task is the MESSAGE variance model** — build the message-propagation precision on the `τ_λ` foundation
> and retire the density-uniformity `σ²_transfer` proxy. Per handoff §6: **(Step 1) DERIVE + Monte-Carlo-validate
> BEFORE any solver code** — the per-component density message variance `Var(log ρ_c)` (graft = SUM /
> share-weighted; peel = DIFFERENCE / u-weighted; the `1/n` sampling enters here), including the pure-gDNA anchor
> limit, and the direction-dependent transfer variance `σ²_transfer = Var(log r)` (~0 on the graft, load-bearing
> on the peel). **(Step 2) IMPLEMENT** as pure tested arithmetic (extend `enrichment_frame.py` →
> `message_arithmetic.py`), folding in the deferred composition/sampling split (move `1/n` out of the fusion
> weight into the transport seed) and the **`struct_lock` HARD OVERRIDE** (a composition-certain node must never
> get an `∞` fusion weight → nan cascade). **(Step 3) LOOP**: full ambig_dense_10mb suite → worst scenario by
> error mass → dissect (`pass0_node_dissect.py`) → fix, with the per-condition A/B (refit=0 AND refit=1,
> unstranded-capON called out) as the gate. Preserve the invariants (handoff §7) and the corrected census
> (handoff §8 — do not delete `rna_*_frac_var`, `strand_likelihood.py`, or `lam_var`). No magic numbers —
> pause and discuss before any new constant. Counts are Poisson. Consider a derivation workflow + independent
> verification for the message-variance derivation, as we did for the foundation.
