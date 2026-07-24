# Session handoff — the MESSAGE variance model: laws SETTLED, M5 landed, single-λ combine is NEXT

**Read `docs/calibration/ROADMAP.md` first, then this. This is the LIVE handoff** — it supersedes
`SESSION_2026_07_24_HANDOFF_3.md` (which began the message-variance task). Date: 2026-07-24. Do NOT read
`docs/calibration/archive/`.

---

## 1. Status in one line

The message-variance **LAWS are derived, MC-validated, and independently verified** (M1–M5); the pure layer is
landed + unit-tested; **M5 `σ²_transfer` is wired and is a large win on the error-mass arm** (aggregate refit=1
0.1234→**0.0945**, unstr-capON 0.387→**0.163**) **but regresses the stranded arm** (0.026→0.077). Dissection
proves the regression is the **two-message per-component combine** (the M6 rank-1 double-count), not the
variance laws. **The next task is the committed fix: the single-λ combine on a λ-space relay** (§6).

## 2. What was done this session (all on branch `calib-ambig-init-wip`, UNCOMMITTED working tree)

1. **Derived the message-variance model** (`docs/calibration/message_variance_derivation.md`) on the `τ_λ`
   foundation; MC-validated every law in `scripts/debug/message_variance_mc.py` (≤4% in-regime, the arbiter).
2. **Independent verification workflow** (`scratchpad/wf_message_variance_verify.js`, run `wf_c952640d`): 3
   independent re-derivations + 2 adversaries + adjudicator, all MC. **GREENLIT**; laws SETTLED; the M6
   two-message double-count proven **rank-1, exactly 2× (up to ~7× with deep spliced)**; **single-λ combine
   committed**; sequencing validated (land core → then single-λ).
3. **Pure layer** (`enrichment_frame.py` — the M1–M5 functions `transport_seed_logvar`, `graft_rna_logvar`,
   `peel_rna_logvar`, `transfer_logvar`, `message_precision`) + closed-form/MC unit tests
   (`tests/calibration/test_enrichment_frame.py`, 61 pass). Handoff §2's "rename to `message_arithmetic.py`"
   is deferred until the combine math also lives there.
4. **Wired M5 `σ²_transfer`** into `bp_solver._unified_solve` (`logvar_tot` per node via `composition_logvar`;
   the `s2t` in `_relay`/`_transport` now `0` on the graft, `logvar_tot[dst]+logvar_tot[src]` elsewhere).
   Retires the density-uniformity NPMLE proxy. Tests green (380 calib+native pass; a `0·inf` warning fixed).

## 3. The A/B (this is the current solver state — `pass0_oracle_bench.py`, `OMP_NUM_THREADS=1`)

| condition | refit=0 base→M5 | refit=1 base→M5 | verdict |
|---|---|---|---|
| aggregate | 0.1267 → **0.1099** | 0.1234 → **0.0945** | ✅ big win |
| unstranded ss_0.50 | 0.1981 → **0.1328** | 0.1902 → **0.1068** | ✅ |
| **unstr capON** | 0.2543 → **0.1761** | 0.3874 → **0.1626** | ✅ (the landmine, −0.225) |
| **stranded ss_0.99** | 0.0262 → **0.0777** | 0.0293 → **0.0773** | ⛔ REGRESSION (gate-fail) |

Baselines shifted vs HANDOFF_3 (0.1491/0.0889) because the suite gained hard `verystrong`/`gdna1`/`gdna5`
scenarios; **gate on in-run A/B deltas, not the absolute.** Worst error mass: high-mass unstranded capON
(gdna300/gdna100) + zero-gDNA false positives (`gdna_none`, mwae ~0.48).

## 4. The stranded regression — DIAGNOSED (do not re-litigate)

* **NOT σ²_transfer placement.** A "matched-exempt" variant (σ²_transfer=0 on all matched reframes, the
  honest §3 placement) made stranded **worse** (0.0863 vs 0.0777). Reverted to graft-only exemption.
* **It is the two-message combine.** `pass0_node_dissect.py` on the worst stranded scenario: strand-alone solve
  is GOOD (boundary mwae 0.05), the imputation MESSAGES corrupt it (→0.20), concentrated on **AMBIG boundaries**
  (`tau0_lam=0` — the strand can't pin f_g for a 2-DOF node at any κ, so messages are the only info, and the
  two-message combine delivers wrong ones: a precision-0 RNA− message + a confident wrong gDNA/RNA+ pair drive
  an oracle-0.255 node to 0.995). Both refit arms regress, so the hyperprior does not rescue it.

## 5. Exact code state

* **Uncommitted** working tree on `calib-ambig-init-wip` (HEAD `438e9f0d`). Changed:
  `src/rigel/calibration/enrichment_frame.py` (M1–M5 laws), `src/rigel/calibration/bp_solver.py` (M5 wiring +
  import), `tests/calibration/test_enrichment_frame.py` (M-law tests), `docs/calibration/*` (derivation doc,
  this handoff). MC arbiter: `scripts/debug/message_variance_mc.py`; workflow: `scratchpad/wf_message_variance_verify.js`.
* Gates: `pytest tests/calibration tests/native` = **380 pass, 3 xfail, 1 xpass**. `tests/test_golden_output.py`
  will now DIFFER (M5 changes output) — **regenerate goldens LAST**, after the combine lands, not now.
* `bp_solver` still computes `mu_proj/var_proj` upstream (fed only to the `_capture` diagnostic now) — vestigial
  for the relay; safe to delete when convenient.

## 6. ▶ THE NEXT TASK — the single-λ combine on a λ-space relay (the committed fix)

**Why.** The composition is ONE DOF, `λ = logit f_g`. The two-message combine hands ψ a gDNA message (on
`log f_g`) AND an RNA message (on `log f_R`) that are **rank-1** (corr −1) — ψ double-counts → 2–7× over-confident,
and it delivers structurally-wrong messages to AMBIG boundaries (§4). Fix: emit **one λ-constraint** per source.

**The key insight (this session).** **λ is ENRICHMENT-INVARIANT** — under a matched reframe both ρ_g, ρ_R scale
by the same `r`, so `r` cancels in `λ = log(ρ_g/ρ_R)+log(E_g/E_r)`; a frame change only shifts λ by the known
eff-length ratio `log(E_g^dst E_r^src/E_r^dst E_g^src)`. So **carry a λ-belief `(λ_mean, Var(λ))` per node
through the relay** (+ a separate θ tilt belief for AMBIG), NOT three per-component densities. This is *simpler*
(no enrichment residual — the current relay laments `Σ_c f_c=75.6` at introns) AND removes the double-count.

**Steps.**
1. **λ-relay:** replace the per-component (density, precision) relay with a `(λ, Var(λ), θ, Var(θ))` relay. Each
   hop shifts λ by the eff-length ratio (enrichment-free), grafts/peels the mature (the peel is a λ-domain
   DIFFERENCE — reuse `peel_rna_logvar`'s u-weighting; the graft folds the spliced into R-total via `w_μ`), and
   adds `σ²_transfer` (M5, `transfer_logvar`) only where load-bearing (peel/partial), 0 on the matched reframe.
2. **Joint λ-precision:** each neighbour's λ-constraint uses `Var(log k)` = the T2 form
   `w_μ²(1/n+1/n_s) + (1−w_μ f_g)²·Var(λ_src)` (MC-validated in `scripts/debug/message_precision_mc.py`) — the
   λ-relay makes `Var(λ_src)` and `w_μ` available, so the joint precision is computable with NO cross-component
   covariance to track.
3. **ψ interface:** feed ψ ONE Gaussian on `λ` (`−½·λ_prec·(λ−λ_msg)²`, a direct grid-variable constraint —
   simpler than the current `log f_c` messages) + the θ message for AMBIG. Two regimes: single-λ where both
   components live; the ÷M gDNA density mode (`message_precision`) at structural anchors (`ρ_R≡0`, k singular).
4. **`struct_lock` HARD OVERRIDE + the 1/n fusion/transport split** (handoff §7, still owed): the fusion weight
   is `1/Var(λ)` (composition only), the transport seed `1/(Var(λ)+1/n)`; a `struct_lock` node ADOPTS its own
   belief (never an ∞ weight → nan). Add an interior-anchor no-nan test.
5. **A/B every step** (refit=0 AND refit=1, unstranded-capON called out); PASS iff no per-condition regression —
   specifically stranded must return to ≤ ~0.03. Then regenerate goldens.

## 7. INVARIANTS — preserve (handoff §7 unchanged) + the corrected census (HANDOFF_3 §8)

`N` enters only as power (`τ_λ` or `1/n`), never a composition vote; variances in log-odds `Var(λ)`/`Var(log f_c)`,
NEVER on simplex fractions; ONE Schur scalar `τ_λ`; fusion weight `1/Var(λ)`, transport `⊕1/n`; `struct_lock` =
hard override. Do NOT delete `rna_*_frac_var`, `strand_likelihood.py`, `lam_var`, `var_gdna` (all LIVE/referenced).

## 8. Methodology + tools

Derive → MC-validate → implement → per-condition A/B → loop; **no magic numbers**; counts are Poisson (ω=0);
`OMP_NUM_THREADS=1` for determinism; cached suite `~/Downloads/rigel_runs/ambig_dense_10mb`. Tools:
`pass0_oracle_bench.py` (the A/B, `--arm`, `P0_REFIT`, `--report`), `pass0_node_dissect.py` (ψ channel ablation),
`message_variance_mc.py` + `message_precision_mc.py` (the MC arbiters), `wf_message_variance_verify.js` (the
verification workflow).

## 9. ▶ KICKOFF PROMPT (copy-paste)

> We are developing Rigel calibration — pass-0. The message-variance LAWS are settled, MC-validated, and
> independently verified (M1–M5); the pure layer is landed; **M5 `σ²_transfer` is wired and wins big on the
> unstranded/capture arm but regresses the stranded arm** — diagnosed (dissection) as the two-message
> per-component combine's rank-1 double-count corrupting AMBIG boundaries, NOT the variance laws. Read
> `docs/calibration/ROADMAP.md`, then `docs/calibration/SESSION_2026_07_24_HANDOFF_4.md` and
> `docs/calibration/message_variance_derivation.md`. Do not read `docs/calibration/archive/`.
>
> **The task is the single-λ combine on a λ-space relay** (handoff §6): carry a λ-belief `(λ_mean, Var(λ))` per
> node (λ is enrichment-invariant — `r` cancels), emit ONE λ-constraint per source with the joint `Var(log k)`
> (T2 form), feed ψ one Gaussian on λ (+ a θ message for AMBIG), two-regime with the ÷M density mode at anchors,
> and land the `struct_lock` hard override + the 1/n fusion/transport split. Validate with the per-condition A/B
> (refit=0 AND refit=1, unstranded-capON called out) — PASS iff stranded returns to ≤~0.03 with no other
> regression; regenerate goldens last. No magic numbers. Counts are Poisson.
