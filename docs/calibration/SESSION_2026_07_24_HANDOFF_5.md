> ## ⛔ SUPERSEDED — THIS IS NOT THE LIVE HANDOFF
> The live handoff is **`SESSION_2026_07_26_HANDOFF_14.md`**; the entry point is **`ROADMAP.md`**.
> This file is kept as HISTORY — what was tried, measured and refuted. Its numbers describe the
> code as it was on its own date and many are now superseded. **Do not act on it without checking
> the live handoff first.** Its own DO-NOT-RE-RUN findings (if it has any) remain HERE and still
> stand — HANDOFF_14 §3 indexes which files carry them.

# Session handoff — the composition-mismatch (DerSimonian–Laird) cliff-precision term, pass-0 cleanup, Phase 2

*(this file's own original header — superseded)* — it supersedes
`SESSION_2026_07_24_HANDOFF_4.md` (still useful for the arc + the audit/design details; read it second). Do NOT
read `docs/calibration/archive/`. Date: 2026-07-24.

---

## 0. TL;DR — what to do next (one paragraph)

The message-variance model is derived, MC-validated, independently verified, and IMPLEMENTED: the **single-λ
combine on a three-stream relay** (double-count fixed) and the **cross-cliff precision term** (a message crossing
a large enrichment cliff must lose precision) are both landed and committed. The shipped cliff term is
`σ²_cliff = (log r)²`, which recovers the stranded regression and is a net win at refit=0, **but it over-damps
the extreme-enrichment (verystrong) arm** because prior-free it cannot tell a pure-enrichment cliff (composition
preserved — should cost nothing) from a composition-mismatch cliff (should cost). **THE NEXT TASK: replace the
`(log r)²` proxy with the honest, MC-derived DerSimonian–Laird composition-mismatch term `b̂²` (§4)** — it damps
ONLY the actual source-vs-destination composition-share mismatch, measured against the node's own self-solve, so
it fixes the stranded arm WITHOUT over-damping enrichment. It is prior-free and magic-number-free. §4 has the
exact formula + the exact code edits. Then the small pass-0 cleanup (§6) and the Phase-2 hyperprior (§7).

## 1. The owner's governing principle (do not violate)

**Pass-0 must be WEAK and correctable.** Its job is not to match the oracle — it is to be a weak, honest
approximation that trains the gDNA hyperprior, so the refit can correct it. **An over-confident message that PINS
a node wrong is worse than a weak one that is slightly off**, because high precision is nearly irreversible in
the refit (measured: the over-confident solver pins the median node at `Var(f_g)≈0`; the fixed solver drops the
median precision ~4 orders of magnitude into a correctable range AND halves the error). Every design choice below
is toward honest weakness: **precision is discrete counts over a discrete length; enrichment scales the MODE, it
must NEVER scale the PRECISION.** When in doubt, prefer the weaker (more damped, under-confident) option — that is
the count-zero-info safe direction.

## 2. What is DONE + committed (branch `calib-ambig-init-wip`)

Commits this arc (newest first):
* `813f18ba` — cross-cliff precision `σ²_cliff=(log r)²` (recovers the stranded regression).
* `e41c4e58` — three-stream relay + single-λ combine (the M6 double-count fix).
* `d78da565` — the ψ single-λ + θ message interface (`simplex_logodds`).
* `44f1ecc6` — the message-variance model: derive+verify the M1–M5 laws, pure layer, M5 σ²_transfer.

Settled facts (do NOT relitigate):
* **The message-variance laws M1–M5** (`docs/calibration/message_variance_derivation.md`): derived, MC-validated
  (`scripts/debug/message_variance_mc.py`), independently verified (workflow `wf_c952640d`).
* **The M6 double-count** (two per-component ψ messages from one rank-1 source are 2–7× over-confident) — FIXED by
  the single-λ combine.
* **The stranded regression root cause** (audit workflow `wf_fe3119fb`): the `σ²_transfer=0` graft/matched
  exemption was wrongly applied to the PRECISION, so a message crossing a 407× cliff arrived at full confidence
  and steamrolled the weak-but-correct own belief. Anchor: node 1909 (a mostly-gDNA exon, oracle 0.985, own
  self-solve 0.953 at weak τ=1.6) collapsed to 0.51.
* **The design workflow (`wf_c565d23a`, 3 derivations + 2 critiques + adjudicator, MC-validated)** produced BOTH
  the shipped `(log r)²` (derivation #1, the conservative prior-free proxy) AND the honest DL mismatch term
  (derivation #3, §4 below). Adjudicator shipped `(log r)²` as the minimal prior-free fix and flagged DL as the
  successor — but derivation #3 shows DL is implementable prior-free NOW; that is the next task.

Gates: `pytest tests/calibration tests/native` = **383 pass, 2 xfail, 2 xpass** (the AMBIG-recovery test now
xpasses — a good sign); `ruff check src/` clean. **Goldens NOT regenerated** (deferred until the solver is
final). Working tree clean except untracked `scratchpad/` (diagnostic scripts + the workflow scripts — keep).

## 3. Current A/B (the state to beat) — `pass0_oracle_bench.py`, `OMP_NUM_THREADS=1`

| condition | refit=0 base→cliff | refit=1 base→cliff | verdict |
|---|---|---|---|
| aggregate | 0.1267 → **0.1092** | 0.1234 → 0.1283 | refit=0 net win |
| stranded ss_0.99 | 0.0262 → 0.0391 (was 0.079) | 0.0293 → 0.0374 | mostly recovered |
| unstr capON | 0.2543 → 0.2186 | 0.3874 → 0.3420 | ✅ |
| **verystrong** (all `ss_0.50`, low-gDNA) | 0.2985 → 0.2779 | 0.1514 → **0.2356** | ⬆ over-damped by `(log r)²` |

`base` = pre-fix two-message (arms `base_r0`/`base_r1` in `/tmp/pass0_oracle_bench.tsv`); `cliff` = current HEAD.
**PASS target for the next task:** stranded ss_0.99 back to ≤ ~0.03 AND unstranded-capON no worse than ~0.163
(refit=1) AND verystrong no worse than the M5/three-stream level (refit=0 ~0.186, refit=1 ~0.129 — the `lam`
arm) — i.e. the DL term should keep the stranded recovery while un-doing the verystrong over-damping. **The A/B
is the arbiter, NOT the anchor** (the anchor is both a big cliff AND a mismatch, so it cannot discriminate the
law choice). Run refit=0 AND refit=1, call out unstranded-capON and verystrong explicitly.

## 4. ▶ THE NEXT TASK — the DerSimonian–Laird composition-mismatch term `b̂²` (replaces `(log r)²`)

### 4.1 Why (the derivation, machine-precision validated — derivation #3)
The per-component ÷M message mode error is EXACTLY `mo_c^msg − log f_c^dst,true = log(s_c^src/s_c^dst,true)`,
the source-vs-destination **composition-SHARE mismatch** (`s_c ≡ ρ_c/ρ_tot`; the reframe `r`, `M_dst`, and the
capture efficiency ALL cancel to machine precision). So the honest delivered variance is

```
    σ²_c,delivered  =  v_src,c  +  b_c²
      v_src,c = the source sampling variance (Var(log f_c)+1/n, or 1/n_spl for a measured junction) — COUNT/length origin
      b_c²    = the composition-mismatch BIAS² = (log s_c^src − log s_c^dst,true)²  — the third source (population gDNA prior)
```

`Var(log r)` and `(log r)²` are the WRONG object: `r`'s own scale/measurement noise is irrelevant to the share
mismatch. `(log r)²` conservatively over-attributes the whole cliff to mismatch (why it over-damps enrichment);
the honest term is `b_c²` alone.

### 4.2 The estimator (prior-free, magic-number-free)
`b_c` needs the destination truth, which we lack — but the destination has its OWN independent estimate: the
`node_init` self-solve. Treat the message and the own belief as two independent estimators of the same
composition (a two-study random-effects meta-analysis) and use the **DerSimonian–Laird between-source estimator**:

```
    b̂_c²  =  max( 0 ,  G_c²  −  v_src,c  −  v_own,c )        G_c = mo_c^msg − mo_c^own   (the observed share gap)
    p_effective,c = 1 / ( v_src,c + b̂_c² )                  (the deflated precision fed onward / to ψ)

      v_own,c  from the τ_λ FOUNDATION:  v_own,g = (1−f_g^own)²/τ_own,  v_own,R = (f_g^own)²/τ_own,  = +∞ when τ_own=0
```

**The three regimes fall out automatically (no gate, no constant):**
* **confident own belief, message conflicts (the anchor, node 1909):** large `G`, finite `v_own` → large `b̂²` →
  the wrong message is killed → the node keeps its correct own belief. **Fixes the stranded regression.**
* **`struct_lock` anchor (intergenic, composition-certain, `v_own→0`):** full `G²` damping → an intergenic node
  cannot be knocked off gDNA by a message. Correct (also subsumes the §6 struct_lock override for the message side).
* **AMBIG / unstranded (`τ_own = 0 ⇒ v_own = ∞`):** `b̂² = max(0, G² − ∞) = 0` → NO mismatch damping → the message
  propagates at its honest `v_src` (which still includes the M5 `Var(log r)` — see 4.4). **The place where messages
  are the only information is exactly where the mismatch term disables itself → the M5 unstranded/capON win + the
  verystrong non-over-damping are preserved.** This is why DL fixes verystrong that `(log r)²` broke.

MC-validated (derivation #3, `cliff_mc_3`): DL recovers the true `b²` to 0.0–0.4% for real mismatches (`b≥1.5`);
at `b≈0` it is positively biased (~0.17) — the SAFE (over-damping) direction, and harmless because an
over-damped message at `b≈0` agrees with the own belief, so the fused mode is unchanged. Matched big cliff
(`r≈425`, `b_g≈0`): full precision, fused `f_g`=0.967 (truth 0.97) — the cliff ALONE does not damp, only mismatch.

### 4.3 The EXACT code edits (`bp_solver.py`, `_unified_solve`)
**Step 1 — restore the M5 honest `Var(log r)` (undo the `(log r)²` proxy).** `(log r)²` was committed as the
proxy; DL replaces it. Revert the three `s2t` sites to the M5 `Var(log r)` (which is honest for the *scale*
sampling and correctly 0 on the graft — keep `logvar_tot`, it is already computed at ~line 321):
* `_relay` (~line 447): `s2t = 0.0 if _gr else (logvar_tot[i] + logvar_tot[s])` (the pre-`813f18ba` form).
* `_transport` (~line 531): `s2t = transfer_logvar(logvar_tot, logvar_tot[src], graft)` (re-add the
  `transfer_logvar` import). Do NOT re-add the `_matched`/`s2t_comp` exemption (that experiment regressed stranded).
* Keep the `RIGEL_S2T_OFF` diagnostic (zero both `s2t` and the DL term for isolation).

**Step 2 — add the DL mismatch deflation in `_transport`**, immediately AFTER `_pin_v` produces the message
densities `tg, tp, tn` and BEFORE the `return`. `og/op/on` (own densities from `_ni`), `E_g/E_r/M`, `_ni.f_g`,
`_ni.tau_lam` are all in scope. Per component `c ∈ {g(E_g), +(E_r), −(E_r)}` compute (vectorized):
```
    G_c   = log( t_c / o_c )                              # = mo_c^msg − mo_c^own; E_c/M cancel. Guard t_c,o_c > _EPS.
    v_own,g = (1−fg_own)²/τ_own ; v_own,R = fg_own²/τ_own # fg_own=_ni.f_g, τ_own=_ni.tau_lam; +inf where τ_own≤_EPS
    v_msg,c = 1/tp?_c  (the transported precision from _dv, i.e. 1/p_msg,c)   # per stream: tpg/tpp/tpn, tmg/tmp/tmn, ttau
    excess_c = maximum(0, G_c² − v_msg,c − v_own,c)
    deflate every precision stream for component c:  p_new = 1/(1/p_old + excess_c)   (== _damp(p_old, excess_c))
```
Apply `excess_c` to the mode-fusion (`tpg/tpp/tpn`), the measurement (`tmg/tmp/tmn`), AND the composition `ttau`
(the anchor recovers only when the composition arm is damped too — ablation-confirmed). Do the SAME in `_relay`
(scalar) if the per-hop accumulation needs it — **but first try DL in `_transport` (the combine) ONLY** and A/B;
the relay may be fine on the honest `Var(log r)` alone. Mirror twin: keep `_relay` and `_transport` consistent.

**Step 3 — anchor `G_c` on the `_ni` SELF-SOLVE densities `og/op/on`, NEVER the running `f_cur`.** `_ni` is the
node's INDEPENDENT evidence (the DL "second study" must be ⟂ the message), and using the stable self-solve avoids
feedback in the 2-iteration `_RHO_ITERS` loop. (`og/op/on` are exactly `_ni.rho_g/rho_pos/rho_neg`.)

**Step 4 — the composition arm is ONE DOF.** `b̂²` is derived per-component, but the composition is a single λ.
The pragmatic first cut (derivation #3) applies it per-component and lets ψ combine; if the stranded arm is
under/over-damped, the principled refinement is to compute ONE `G_λ = (mo_g−mo_R) − (log f_g^own − log f_R^own)`
for the composition τ-stream and per-component gaps only for the independent measurement streams. Decide by A/B.

### 4.4 MC-VALIDATE FIRST, then implement, then A/B (the standing methodology)
1. Extend `scripts/debug/message_variance_mc.py` (or reuse the workflow's `scratchpad/cliff_mc_3.py` if present) to
   pin: (a) DL recovers the true `b²` for real mismatches; (b) `b̂²=0` when the message agrees with the own belief
   (matched big cliff); (c) `τ_own=0 ⇒ b̂²=0` (AMBIG no-damp). Get it passing BEFORE solver code.
2. Implement §4.3. Unit-check node 1909 recovers to ~0.95 (`scratchpad/dump_node.py`).
3. Full A/B (refit=0 AND refit=1). PASS per §3. If a regression appears, **DISSECT** (worst scenario → worst
   nodes → trace the message provenance to root cause; `scratchpad/dump_node.py` + `RIGEL_S2T_OFF`) — do NOT
   assume a theory flaw; the theory is validated. Iterate. Regenerate goldens LAST.

### 4.5 The disagreement to be aware of (resolve empirically)
Derivation #3 says DL is shippable PRIOR-FREE (AMBIG no-damp is CORRECT — it preserves M5). The adjudicator
worried DL leaves AMBIG unprotected from a "wrong loud neighbour." Keeping the M5 `Var(log r)` (§4.3 step 1)
gives AMBIG the honest scale damping even at `b̂²=0`, which mitigates this. The A/B (verystrong + unstranded-capON)
decides; if AMBIG is dragged wrong, the fix is the hyperprior `v_own` (§7), not re-strengthening pass-0.

## 5. FALLBACK — if DL does not beat `(log r)²`
`(log r)²` is committed (`813f18ba`) and is a net win at refit=0. If DL under-performs, KEEP `(log r)²` and take
the verystrong recovery to Phase 2 (§7). Do NOT ship a re-strengthened pass-0 (violates §1).

## 6. Small pass-0 cleanup (do after §4 lands, low risk)
* **`struct_lock` HARD OVERRIDE in `_fuse`** (handoff-4 §7, still owed): a composition-certain node ADOPTS its own
  belief; never an `∞` fusion weight → nan. NOTE: with the three-stream, `tau_own=0` for intergenic and
  `mg_own=n` (finite), so there is likely NO live `∞`-weight today — verify with an interior-anchor no-nan test
  before adding machinery. The DL term (§4.2, `v_own→0` regime) already makes an intergenic node immovable by a
  MESSAGE; the override is about its own FUSION. Add the no-nan regression test regardless.
* **Retire the dead σ²_transfer NPMLE.** The `enrichment_prior` / `DensityNPMLE` was the OLD σ²_transfer source
  (`calibrate.py` fits it; `node_sweep` took `enrichment_prior`/`transfer_variance`). σ²_transfer is now `(log r)²`
  → DL, computed from `r`, NOT the NPMLE. Confirm `enrichment_prior` is unused in `node_sweep` and remove the
  dead plumbing (`calibrate.py` ~lines 331/374–388/443, the `enrichment_prior`/`transfer_variance` params). Keep
  `composition_logvar`/`logvar_tot` only if §4.3-step-1 keeps the M5 `Var(log r)`; else remove.
* **The `1/n` fusion/transport split (handoff-4 §7) is effectively DONE** by the three-stream: the composition
  τ-stream is seeded `τ_λ` (composition only, NO `1/n` — the count cancels in the ratio, a refinement of the
  handoff's blanket `⊕1/n`), and the measurement stream carries the counts. No further action; just confirm.

## 7. Phase 2 — the gDNA hyperprior refit (the ROADMAP's next workstream)
The refit machinery EXISTS (`calibrate.py`: `for it in range(calib_refit_iters): gdna_hyperprior =
_fit_gdna_hyperprior(...); belief = _sweep(gdna_hyperprior, ...)`; the hyperprior enters ψ's gDNA arm as
`gdna_prior=`). Known landmine (ROADMAP §4): the refit currently regresses unstranded-capON. Phase-2 goals:
1. **Feed the hyperprior into DL's `v_own` on AMBIG nodes.** Prior-free, AMBIG has `τ_own=0 ⇒ v_own=∞ ⇒` no
   mismatch protection. The trained hyperprior supplies a finite `v_own` there → DL protects AMBIG nodes from
   wrong loud neighbours (tightens §4.5) → recovers any residual AMBIG error + verystrong.
2. **Fix the refit unstranded-capON regression** (dissect with the same loop).
3. Then the re-solve → a ship candidate. Compare vs the legacy-with-factory baseline (mwae 0.0949) — but on the
   CURRENT suite (which gained hard verystrong/gdna1/gdna5 scenarios), so gate on in-run A/B deltas.

## 8. INVARIANTS + gotchas (must preserve)
* **Honest precision** (§1): precision = discrete counts over discrete length; enrichment scales the MODE, never
  the PRECISION; pass-0 weak + correctable; prefer under-confident when unsure.
* `N` enters only as power (`τ_λ` Fisher, or `1/n` sampling), never a composition vote.
* Variances in log-odds `Var(λ)`/`Var(log f_c)`, NEVER on simplex fractions.
* The composition is ONE DOF (single-λ message); the tilt θ is a SEPARATE DOF (its own message, AMBIG only). Do
  NOT send two per-component composition messages (the M6 rank-1 double-count).
* Counts are Poisson (ω=0); the synthetic suite is Poisson by construction — no overdispersion term.
* `OMP_NUM_THREADS=1` for A/B determinism. Do NOT delete `rna_*_frac_var`, `strand_likelihood.py`, `lam_var`,
  `var_gdna` (all LIVE/referenced — handoff-4 §8).
* Cached suite: `~/Downloads/rigel_runs/ambig_dense_10mb` (`_selfsolve_cache`, ~1 s/scenario).

## 9. Key files + tools
* Solver: `bp_solver.py` — `node_sweep` → `_unified_solve` (`_relay`, `_transport`, `_fuse`, `_pin_v`, the three
  streams: full `pg/pp/pn`→mode, measurement `mg/mp/mn`→gdna/rna_imp, composition `τ`→λ-message; the `s2t`/DL
  sites; `_ni` = `build_node_init`; the `_capture` diagnostic dict publishes `mo_*`, `cm_*`, `c_tau`).
* ψ: `simplex_logodds.py` — `_solve_nodes_logodds_all` (the `lam_imp`/`theta_imp` single-λ interface).
* Pure laws: `enrichment_frame.py` (`composition_logvar`, `transfer_logvar`, the M1–M5 variance laws).
* `node_init.py` (`build_node_init`, `tau_lam`, `rho_*`), `calibrate.py` (the refit loop).
* Docs: `message_variance_derivation.md`, `SESSION_2026_07_24_HANDOFF_4.md` (the full arc + the audit/design
  findings), this file.
* Tools (`scripts/debug/` + `scratchpad/`): `pass0_oracle_bench.py` (the A/B, `--arm`, `P0_REFIT`, `--report`),
  `pass0_node_dissect.py` (ψ ablation), `message_variance_mc.py` (the MC arbiter), `scratchpad/dump_node.py`
  (per-node message dump), `scratchpad/pin_audit.py` (the weakness/pinning audit), `scratchpad/trace_stranded.py`
  (σ²t ON/OFF node ranking), env `RIGEL_S2T_OFF=1` (disable the cliff term). The workflow scripts that produced
  the derivations: `scratchpad/wf_cliff_precision_design.js`, `wf_precision_audit.js`.

## 10. ▶ KICKOFF PROMPT (copy-paste)

> We are developing Rigel calibration — pass-0. The message-variance model, the single-λ combine (double-count
> fix), and a cross-cliff precision term are all implemented and committed. The shipped cliff term `(log r)²`
> recovers the stranded regression but OVER-damps the extreme-enrichment (verystrong) arm, because prior-free it
> cannot separate a pure-enrichment cliff (composition preserved — should cost no precision) from a
> composition-mismatch cliff (should). Read `docs/calibration/ROADMAP.md`, then
> `docs/calibration/SESSION_2026_07_24_HANDOFF_5.md` (START HERE — it has the exact plan), then
> `docs/calibration/SESSION_2026_07_24_HANDOFF_4.md` and `message_variance_derivation.md` for the arc. Do NOT read
> `docs/calibration/archive/`.
>
> **The task is to replace `(log r)²` with the honest DerSimonian–Laird composition-mismatch term `b̂²`**
> (HANDOFF_5 §4): the exact cross-cliff error is the per-component composition-SHARE mismatch `b_c`, estimated
> prior-free against the node's own `node_init` self-solve as `b̂_c² = max(0, G_c² − v_src,c − v_own,c)`,
> `G_c = log(t_c/o_c)`, `v_own` from `τ_λ`. It damps ONLY genuine composition mismatch (fixes stranded, node 1909
> → ~0.95) and disables itself where `τ_own=0` (AMBIG/unstranded → preserves the M5 win + un-does the verystrong
> over-damping). Restore the M5 `Var(log r)` for the scale-sampling term and add `b̂²` in `_transport` (the exact
> code edits are in §4.3), anchoring the gap on the `_ni` self-solve. **MC-validate the DL law FIRST**
> (`message_variance_mc.py`), then implement, then the full A/B (refit=0 AND refit=1, unstranded-capON AND
> verystrong called out) — PASS iff stranded ss_0.99 ≤ ~0.03 with no regression on unstranded-capON or verystrong.
> If a regression appears, DISSECT (worst scenario → worst nodes → trace message provenance; `scratchpad/
> dump_node.py`, `RIGEL_S2T_OFF`) — the theory is validated, find the root cause. Then the pass-0 cleanup (§6) and
> Phase 2 (§7). GOVERNING PRINCIPLE: pass-0 must be WEAK and correctable — an over-confident message that pins a
> node wrong is worse than a weak one that is slightly off; prefer the under-confident option when unsure. No
> magic numbers. Counts are Poisson. `OMP_NUM_THREADS=1`. Consider derivation/audit workflows as we did this
> session (they worked well). Regenerate goldens LAST.
