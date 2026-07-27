# SESSION 2026-07-27 — HANDOFF 16 (LIVE)

Supersedes `SESSION_2026_07_26_HANDOFF_15.md`, which is now **⛔ SUPERSEDED** except for its §5 invariants
and §7 do-not-re-run table (both still current). This session was an **optimization / cleanup / audit**
session, not a modelling one. **Next task: the gDNA hyperprior (Phase 2).**

---

## 0. ⚠⚠ READ THIS FIRST — every stored baseline is STALE

One commit this session **changed the solver's output**: `08cfa0e9` (the peel level's log-variance was
capped at π²/6). Everything else was verified bit-identical, but that one moves numbers.

**Consequences, and they have already cost this project time once:**

* the `clean0_*` and `head0_*` arms in `/tmp/pass0_oracle_bench.tsv` no longer reproduce at HEAD;
* `head1_r0` / `head1_r1` were re-recorded at `08cfa0e9` and are the valid baseline — but **`/tmp` is
  volatile.** If it has been cleared, re-record from a `git stash` of HEAD, in the same session, at BOTH
  refits, before comparing anything.
* **If a HEAD-vs-baseline run is not 32/32, the baseline is what is broken, not your change.** This session
  lost ~40 minutes to a stale `clean0` reporting 0/32 for a change that was in fact bit-identical —
  including on the `mass` column, which comes from the accumulator payload and no solver change can reach.

**Current numbers at HEAD** (32 conditions, mass-weighted, pooled over all live nodes):

| | suite mwae | fit substrate |
|---|---|---|
| refit 0 | **0.084566** | 0.085304 |
| refit 1 (production) | **0.067973** | 0.068407 |

## 1. What this session did

**Performance — genome scale 266 s → 67 s (~4×), all bit-identical.** `619465c8`, `6d483004`, `3292adbe`.
The whole story, including why the 10 Mb synthetic ranks the hotspots *backwards*, is in
`PERFORMANCE.md` (rewritten). Headline: profile on `cfrna:LBX0190`, never on the toy — the relay is 81 % of
a real calibration and 34 % of the toy's, and the toy's own #1 hotspot (the NPMLE EM) is 0.7 % at scale.

**Tests — 1202 pass, 0 xfail, 0 xpass, 0 warnings.** `36356527`. Four stale markers retired; one of them
documented an OPEN ITEM that is now **resolved** (|Δf_g| 0.0612 → exactly 0.000000).

**Dead code and stale docs.** `e1a041e6`, bit-identical. Removed `node_sweep`'s never-referenced
`boundary_substrate` parameter (15 call sites), `_RHO_ITERS`, an unreachable `_rna_arm` branch, and four
actively-misleading doc claims — including a docstring citing a function that does not exist and a module
docstring saying the variance model "is being REDONE" when it is COMPLETE.

**The over-engineering audit.** `0ce50a12`, `75b4279b`, `4a55661b`. Full results in
`ABLATION_CAMPAIGN.md`; the P1d/P1e brief is `P1D_P1E_DEBTS.md`.

**One real bug found and fixed.** `08cfa0e9` — see §0 and §2.

## 2. The finding that matters most for Phase 2

`_peel_share` was delivering the peel level's log-variance through a lossy round trip, and
`residual_level`'s `k ≥ 1` floor — correct there, where `k` is a fragment COUNT — became a hard **ceiling of
`ψ'(1) = π²/6`** on it, because the re-derived `k` is a reciprocal variance. A level M11 declared 10 nats
uncertain was delivered as 1.6 (6× over-confident), on exactly the seams where it is least determined.

**Why Phase 2 should care:** the hyperprior is fitted *precision-weighted* by `var_gdna`. Every declared
variance the solver produces is an input to that fit, so a systematically over-confident stratum biases it.
This fix widens exactly the weakest seams. **The hyperprior has never been fitted on a solver without that
ceiling** — treat any pre-`08cfa0e9` hyperprior measurement as suspect.

## 3. ▶ NEXT: the gDNA hyperprior (Phase 2)

Entry points, in order: `ROADMAP.md` → `gdna_hyperprior_resurrection.md` → `npmle.py`'s module docstring
(the two roles) → `calibrate._fit_gdna_hyperprior`.

What is already established and should not be re-derived:

* Phase 2 is **not blocked** — the refit regression and the AMBIG-poisoning concern were both stale
  (substrate z2 = 2.26; AMBIG is excluded from the fit by construction).
* Role B (the composition arm of ψ on the re-solve) is **resurrected and no longer blind** in enriched
  mode (mass 0 → 0.70, exact on stranded). One directional bias survives: **unstranded enriched under-call
  0.15–0.22 decades**. A *symmetric* projection cannot fix a bias — an oracle landscape buys only
  +0.03…+0.08.
* `HANDOFF_15` §2 proposed **P2a**, the 6-line λ-curvature step, as the cheapest high-value experiment.
  That is still the suggested first move.

Two cautions carried forward:

* the caches in `scratchpad/` are the good ones; **two label bugs in `gdna_cache_build.py` make archived
  per-capture / per-DNA numbers suspect** (that script is currently modified in the working tree — the
  owner's own work, left untouched).
* the intron factory has an **open flag**: "refit TRIPLES unstranded-capON error". It is also the only
  source of `τ_λ > 0` on unstranded data, which makes it the gate on M7, on M5's composition half, and on
  the λ-message all at once.

## 4. Deferred, with the reasoning already done

* **M7's DL mismatch deflation has never been measured on real data**, and it is the next audit target.
  Its own docstring concedes it is inert wherever `τ_own = 0` — all AMBIG, all unstranded, ~84 % of the
  error mass — and its gate `v_own` is supplied almost entirely by the intron factory. It was landed
  2026-07-25, before five substrate changes. `variance_ledger.md` §2.2 says M7 absorbs the *residual*, and
  five terms have landed since; the question is whether anything is left for it.
* **The remaining ablation arms**: `nofuse × noP1e`, `nofuse × noP1d` (the fuse fix touches the peel's
  variance, so it may move the debts' verdicts), the variance-freeze reference, the duplicate λ-emission
  gate at the combine, the Schur AMBIG τ gate. Harness is built; ~15 min per round.
* **Further performance work is PAUSED by owner decision** until the full pipeline is finished.
  `PERFORMANCE.md` §6 records where the remaining time is and what the options cost.

## 5. Invariants — unchanged, still binding

`HANDOFF_15` §5 in full, plus two added this session:

* **`_relay` and `_transport` are twins — do not merge** (now also: the relay is scalar-native throughout,
  which makes a merge worse still). Likewise `residual_level` / `residual_level_scalar`,
  `peel_continue_share` / `_scalar`, `graft_frame_logvar` / `_scalar`, `_peel_share` / `_peel_share_scalar`
  — **a change to one must be mirrored in the other**, and each pair's equality is pinned bit-for-bit by a
  test in `test_enrichment_frame.py`.
* **`simplex_logodds._SOLVE_BLOCK_BYTES` is a cache knob, not a model parameter.** Every node solves
  independently and every reduction is within a row, so any value is bit-identical. Do not "tune" it for
  accuracy; there is none to gain.

## 6. Tooling added this session (all in `scratchpad/`, untracked)

| script | what it is for |
|---|---|
| `perf_scale.py` | profile/time calibration at toy **or** genome scale off the cached payloads |
| `perf_identity.py` | **the fast identity gate** — full per-node solved state, 8 conditions × both refits, ~8 s. Use this as the inner loop; `verify_clean.sh` is the ~22 min landing gate |
| `perf_shapes.py` | report the real (m, K) / (m, K, K_t) shapes each solver sees |
| `perf_twin_check.py` | bitwise battery for the scalar twins (471,824 cases) |
| `ablate_campaign.py` / `ablate_report.py` | the ablation harness: arms, the 2×2 interaction cell, held-fixed z2 |
| `p1de_firing.py` | the P1d/P1e firing diagnostic on real cfRNA |
| `verify_clean.sh` | the landing gate — now asserts both arms carry 32 rows and takes `BASE=` |
