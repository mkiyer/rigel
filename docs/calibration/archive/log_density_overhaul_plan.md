# Log-density overhaul (Phase 3 #2) — implementation plan

**Status:** design complete, ready to implement. Built from `log_density_1d_solver_design.md` +
an adversarial design review (workflow `wf_8f6e0421`, 10 agents, 48 issues) cross-checked against a
direct read of every touched subsystem. Goal: migrate the calibration MODEL — messages, the global,
and σ²_bio — into log-density space, fixing the §3.3 "fit-in-log / solve-in-linear" inconsistency.
The user's directive: a complete overhaul, **refit in log for all**, simpler and more elegant.

---

## 1. Scope — four coupled changes, all in log-density (the logodds path only)

1. **Messages** (`bp_solver.node_sweep._scan`): send a log-density Gaussian
   `−½·prec_log·(log ρ_c_dst(λ) − μ_log)²`, `μ_log = log ρ_src`. **The fraction conversion and the
   `[0,1]` mode clip hack disappear** (the source's imputed density is the destination's prior mean
   directly — identity-mean).
2. **Precision** composes in log space: `prec_log = 1/(Var(log ρ_c)_src + σ²_bio_log + pois_log)`.
   The `(M/E)²` density Jacobian on the precision **drops** (in log space `d(log ρ_c)/dλ = O(1)`).
3. **σ²_bio** refit directly in log space: `MonotoneVarMean` on `raw = (Δ log ρ)²` vs the linear
   density level, `offset = 1/(ρ_src·E_src + 1) + 1/(ρ_dst·E_dst + 1)` (the inverse-count log-Poisson
   floor). The spline/GCV/IRLS machinery is reused unchanged — only `_edge_varmean`'s inputs change.
4. **Global** (`_global_logprior`): a Gaussian on `log ρ_g`, `−½·N_log·(log ρ_g(λ) − μ_log_g)²`,
   replacing the Beta-on-`f_g` `_binom_pseudo`. `μ_log_g` = a **log-space** enrichment blend;
   `N_log = 1/(var_mean_log + σ²_g_log)`, M-independent.

The bridge throughout: `log ρ_c_dst(λ) = log d_dst + log f_c(λ)`, `d_dst = M_dst/E_c_dst` the observed
total density (count-zero-info: `d` observed, only the split `λ` is free).

---

## 2. Resolved design decisions

**D1 — BRANCH `_scan` / `_global_logprior` / `_edge_varmean` on `use_logodds`; do NOT retire the lattice.**
The message-build, the global, and the σ²_bio fit are *shared* across both backends. Migrating them
unconditionally silently corrupts the lattice path (147 goldens, the 2 fraction-message unit tests, the
zero-gDNA-pin test). Branching keeps the lattice **byte-identical** as the A/B baseline through Phase 4
(when it is deleted → the elegant logodds-only end state). The emission gate (`emit_g/p/n`, structural)
and the FB α/β scan / `_comb` plumbing stay **shared, above the branch**. Branching is a localized
addition; retiring now is a wide, premature deletion that loses attribution during two coupled scale
changes. *(Confirmed by 4 independent stress agents.)*

**D2 — The precision state is `Var(log ρ_c)` per component, NOT `Var(λ)`.** The gDNA message lives on
`log ρ_g`, so its precision must be the variance of `log ρ_g`, moment-matched on the grid:
`Var(log ρ_g) = post@Lg² − (post@Lg)²`, `Lg = log_expit(λ)` (= `log f_g`); analogously
`Var(log ρ_pos/neg)` for the RNA channels. `Var(λ) = Var(logit f_g)` differs by the `log(1−f_g)` term
(`Var(log ρ_g) ≈ (1−f_g)²·Var(λ)`) — mixing them is a silent units bug. **Three separate per-node
scalars**, each in its channel's currency. *(This is NOT the pathological fraction-var vertex-collapse:
`Var(log ρ_g)→small` at `f_g→1` correctly means a confident pure-gDNA node — it SHOULD emit a confident
gDNA message. The fraction-var collapsed at BOTH vertices regardless of evidence; `Var(log ρ_c)` does
not.)*

**D2a — Cap the sent precision at the predictive ceiling** `prec_log ≤ 1/(σ²_bio_log + pois_log)`. A
node whose true posterior extends past the `±L` window piles mass at the edge → `Var(log ρ_c)`
under-reports → over-confident message. The cap (a node can never be more certain than the noise model
allows — `message_state_separation` S3, in log space) bounds this without widening `L` blindly.

**D3 — ONE reference Jacobian; the log-density empirical priors are added bare.** Phase-1's
`log(f_g(1−f_g))` Jacobian converts the uniform-λ Riemann sum into the uniform-`f` reference measure for
the strand LL + Jeffreys + spliced (the composition reference, f-space). The migrated message/global
Gaussians are **likelihood-like factors** added to ψ — they get NO additional Jacobian. Keep exactly one
`log(f_g(1−f_g))` (1-D) / `log(f_g(1−f_g)²)` (2-D) factor; do not re-apply it per log-density term.

**D4 — `var_mean_log = 1/(1+G)` (the DELTA METHOD — the design's own 1/count principle); refit-in-log
only the SPREAD curves.** "Refit in log for all" applies to the curves that HAVE samples — σ²_bio, σ²_g,
σ²_level (between-node spreads, fit on edge points). `var_mean` is a single pooled rate variance (one
number from `G`, `E_tot`), nothing to refit. Its log-space analogue is the delta method
`Var(log ρ) ≈ Var(ρ)/ρ² = var_mean/ρ_global² = [(1+G)/E_tot²]/[(1+G)/E_tot]² = 1/(1+G)` — i.e. the
log-variance of a count is `1/count`, **exactly the same formula as the per-node `pois_log = 1/(ρE+1)`**,
applied to the global's pooled effective count `1+G`. No special function (trigamma was the exact-Gamma
form — rejected as opaque; delta is transparent, consistent, and removes a concept). At `G=0`:
`1/(1+0)=1` → `N_log = 1/(1+σ²_g_log) ≈ 0.5`, the **gentle M-independent zero-gDNA tiebreaker** the
no-phantom property depends on (`test_gdna_sweep_zero_gdna_pin`), if anything slightly firmer toward
`f_g→0` (helps the phantom — verified at Gate 3).

**D5 — Consistent density floor at the existing pseudocount; no new constant.** Every `ρ` is floored to
`ρ_floor = _MSG_PSEUDOCOUNT/E = 1/E` BEFORE any `log`, used IDENTICALLY in `μ_log`, `pois_log`, and the
σ²_bio offset (so they agree). `pois_log = 1/(ρ·E + 1)` (finite at `ρ→0`, `=1` floor). This is the
"a density from a finite count can never be exactly zero" principle already encoded in `_MSG_PSEUDOCOUNT`
— it stays at ≤8 constants. **Truly-empty nodes (`mass ≤ ε`) are short-circuited before ψ is built**
(they are unsolvable — `solvable=(fp|fn)&mass>0` — and would otherwise inject a `−∞·(+10.4)` NaN row).

**D6 — `scipy.special.log_expit` for `log σ(λ)` and `log_expit(−λ)` for `log(1−σ(λ))`; cap `L ≤ 36`.**
The naive `log(clip(σ(λ),ε))` is off by 15 decades in the depleted tail (`1−σ(40)` underflows to exactly
0). `log_expit` stays exact to `L=36`. Replace the Phase-1 Jacobian's clip too.

**D7 — Enrichment blend in LOG space:** `μ_log_g = w·log(max(ê(z),ρ_floor)) + (1−w)·log ρ_global`
(Jensen-clean; `log(ê)` is the spline directly since `ê` is already an `exp(spline(log z))` log-log fit).
`σ²_level` moves to a log² response alongside `σ²_g`.

**D8 — Spliced single-owner = the standalone floor; the message carries NASCENT only.** *(Refined from
the review's "fold into the message" — that breaks the no-neighbour case: a spliced-bearing node with no
neighbour message would lose its floor entirely.)* The principled split: **nascent RNA imputes across
edges** (the message), **spliced/MATURE is a node-LOCAL observation** (its own junctions) that does not
cross. So: the RNA message mode is `log(nascent_src) = log(f_pos_src·M_src/E_rna)` — **drop the
source-spliced `SPs/esp` from the message ρ and drop the `mo=(ρ·erd−SPd)/md` dst-spliced subtraction**;
keep the standalone clipped-Gaussian `short_p/short_n` as the node's own spliced self-floor. No
double-count (message = neighbour's nascent; floor = own spliced — different quantities), and the floor
binds with or without a neighbour. This is *simpler* than the current code (which mixed source+dst
spliced into the message) — the messy `SPs`-in / `SPd`-out cancellation is gone.

---

## 3. Pitfalls → mitigations (deduplicated from the 48 issues)

| # | severity | pitfall | mitigation |
|---|---|---|---|
| P1 | blocker | `log(0)` / `1/(ρE)→∞` at empty faces, `f_c=0` sources, AMBIG `τ=±1`, zero-gDNA | D5: floor `ρ` at `1/E`; `pois=1/(ρE+1)`; short-circuit empty nodes |
| P2 | blocker | `_scan`/global/σ²_bio shared → migrating breaks the lattice | D1: branch on `use_logodds` |
| P3 | blocker | precision units: `Var(λ)` vs `Var(log ρ_c)` | D2: per-component `Var(log ρ_c)`, capped (D2a) |
| P4 | blocker | Jacobian double-count over the log-density terms | D3: one reference Jacobian only |
| P5 | blocker | `var_mean` "refit in log" has no samples → re-introduces phantom | D4: delta-method `1/(1+G)` (= 1/count) |
| P6 | major | `ρ_src ≫ d_dst` → message mean off-grid → monotone edge pull ("new slam") | physically real (more gDNA than capacity → `f_g→1`); bounded by grid + honest precision. **Verify, don't pre-clip**; contingency = humble precision on conflict |
| P7 | major | log-space σ²_bio de-offset `z=raw−off` often negative → fit collapses to 0 → ∞ precision | floor offset at `1/(ρE+1)` (O(1) like `raw_log`); floor σ²_bio_log at the 1-pseudo-obs log-variance |
| P8 | major | AMBIG `(λ,τ)` cube: RNA messages now τ-dependent → can't factor out → memory blows | node-chunk the AMBIG subset (already flagged; now mandatory) |
| P9 | major | `inf`/`0` init sentinels: unsolved AMBIG must send `prec=0`, not flat-grid var | explicit `isfinite` check in the send; don't compute grid-var for unsolved |
| P10 | major | `L`-edge truncation now corrupts sent PRECISION, not just readout | D2a cap + re-verify `L` in the `(K,L)` sweep (cannot be deferred) |
| P11 | minor | enrichment Jensen bias | D7: blend in log |
| P12 | minor | spliced double-credit | D8: single owner |

---

## 4. Implementation order (each step builds; integrated verification at the end per "one-shot refit-for-all")

All work is on the **logodds branch**; the lattice path is untouched (byte-identical).

1. **Numerical infra** (`simplex_logodds.py`): `_log_fg(λ)=log_expit(λ)`, `_log1m_fg(λ)=log_expit(−λ)`;
   a `_floor_density(ρ, E)` helper; cap `L≤36` in config validation. Replace the Phase-1 Jacobian clip
   with `log_expit`.
2. **Precision-state readout** (`_solve_nodes_logodds` / `_solve_ambig_logodds`): compute + return
   `Var(log ρ_g)`, `Var(log ρ_pos)`, `Var(log ρ_neg)` (moment-matched on the grid), each capped at the
   predictive ceiling. These replace the fraction `var_*` **for the logodds send only** (NodeBelief
   carries them; the lattice path keeps fraction var — the field is internal to `node_sweep`).
3. **σ²_bio log refit** (`_edge_varmean` + `fit_gdna_varmean`/`fit_rna_varmean`, behind a `log_space`
   flag threaded from `node_sweep`): `means=0.5(log ρ_dst+log ρ_src)`'s **linear** level for the
   predictor, `raw=(Δlog ρ)²`, `offset=1/(ρ_dst E_dst+1)+1/(ρ_src E_src+1)`. `MonotoneVarMean` reused.
4. **Messages** (`_scan` logodds branch): `μ_log=log(floor(ρ_src))`, `prec_log` per D2/D2a/D5; consumed
   as `−½·prec_log·(log ρ_c_dst(λ)−μ_log)²` in the logodds solver's message term (D8: spliced folded into
   `log ρ_s`, standalone floor removed). Lattice `_scan` branch unchanged.
5. **Global** (`_global_logprior` logodds branch): Gaussian-on-`log ρ_g`, `μ_log_g` per D7,
   `N_log=1/(trigamma(1+G)+σ²_g_log)` per D4. Lattice `_binom_pseudo` branch unchanged.
6. **AMBIG chunking** (P8): node-batch the AMBIG subset in `_solve_ambig_logodds`.
7. **Wire + Jacobian audit** (D3): one reference Jacobian; migrated terms bare.

---

## 5. Verification protocol

**Gate 0 — four non-negotiable invariants (assert BEFORE any accuracy metric):** (a) per-node mass
conservation `mass_gdna+mass_rna == node total`; (b) `f_g∈[0,1]`; (c) `CalibrationResult.__post_init__`
intrinsic invariants; (d) count-zero-info (the count enters only as strand precision + the scale `d`).
The readout stays ONE scalar `f_g=σ(λ_median)`, `rna_mass=(1−f_g)·M+M_spliced` — never reconstruct mass
from independent `log ρ_g`/`log ρ_rna` means (would break conservation).

**Gate 1 — lattice byte-identity:** default `calibration_solver="lattice"` → 147 goldens + 188 unit
tests byte-identical (D1 guarantees it). New PARALLEL logodds unit tests for the log-density message +
global (clip-hack removal does not vertex-slam; the gentle zero-gDNA pin holds).

**Gate 2 — the THREE-way benchmark** (`scripts/sim/evaluate_suite.py`, `net_flow_per_condition`):
(1) lattice `@0ab636e1`, (2) logodds Phase-3#1 (fraction messages), (3) logodds Phase-3#2 (log-density).
Tabulate `net_gdna_to_rna` per condition side-by-side. **Do NOT trust net-leak alone** — the lattice's
phantom-MASKING is the trap.

**Gate 3 — zero-gDNA is the PRIMARY gate** (not a footnote): the zero-gDNA per-region gDNA mass must
strictly **decrease** vs both the `@0ab636e1` baseline AND the Phase-3#1 logodds +73% number. A
log-density global that over-pins `f_g→0` can fake a good net-leak while being wrong — instrument the
pin strength (`N_log` at `G=0` ≈ 0.5–0.6, gentle) explicitly.

**Gate 4 — `(K,L)` sweep vs ground truth** on the dissect cache (P10: `L` now sets sent precision).
Cheap intermediate A/B captures on the cache for attribution (the two coupled scale changes — Jacobian
removal + linear→log σ²_bio — can't be attributed from one suite run).

---

## 6. Open questions to confirm before coding

1. ~~`var_mean` log-space form~~ **RESOLVED: delta method `1/(1+G)`** (D4) — the design's own 1/count
   principle, no special function (trigamma rejected as opaque). User-confirmed 2026-06-28.
2. **W_regime (Student-t robust message weight) is OUT of scope.** RESOLVED: out (user-confirmed). Still
   measure the off-cap-unstranded ss50 cliff at all three A/B points — the log-density message's
   *multiplicative* "close" may soften it for free.
3. **P6 (off-grid monotone message): verify-don't-pre-clip.** RESOLVED: proceed without a clip
   (user-confirmed), monitor; if a vertex-slam shows from `ρ_src ≫ d_dst`, fix via precision (humble on
   conflict), not a mode clip.
