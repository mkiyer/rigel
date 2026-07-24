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

### ✅ DONE this session — the ψ λ/θ interface (Step A)
`simplex_logodds` now accepts a SINGLE Gaussian on the grid variable `λ` directly (`lam_imp_mode/prec`, both
the 1-D and 2-D paths) + a SEPARATE Gaussian on the tilt `θ` (`theta_imp_mode/prec`, AMBIG only). Unit-tested +
pinned (`tests/calibration/test_lambda_message.py`), backward-compatible (defaults None). `_local_solve` in
`bp_solver` threads them (`lam_imp`/`theta_imp`). This is the substrate the relay feeds.

### ⚠ FINDING — a combine-only single-λ CANNOT work; it must be a THREE-STREAM relay
A v1 that collapsed the existing density relay's fused densities into one `λ_msg = mo_g − mo_R` at the combine
FAILED (`test_tau_gag_fix_spliced_junction_emits_when_unstranded`): on an unstranded spliced junction the
spliced measurement stopped moving f_g, because `Var(λ) = v_g + v_R = ∞` when there is no gDNA info
(`v_g = 1/cpg = ∞`) — even though the RNA (spliced) side is strongly measured. ROOT CAUSE (traced, not
theorised): the density relay ENTANGLES two kinds of evidence that must be handled OPPOSITELY —
* **composition** (from strand/factory-solved neighbours) is RANK-1 → must be ONE λ-message (the M6 fix);
* **independent measurements** (the spliced RNA count, the anchor gDNA count) are NOT rank-1 → must be FUSED as
  separate constraints (an RNA-only spliced measurement constrains f_g via f_R with NO gDNA info needed).

Differencing `mo_g − mo_R` demands both components be known, so it drops the RNA-only channel. The fix separates
**three streams**, all fused independently at ψ:

**Steps (the three-stream relay).**
1. **Composition-τ stream:** seed `τ_own = _ni.tau_lam` (the Schur composition precision; 0 for anchors +
   unstranded non-factory nodes). Transport enrichment-free (λ shifts by the eff-length ratio; τ damped by
   `σ²_transfer` on peel/partial, unchanged on matched — reuse `peel_rna_logvar` u-weighting + `transfer_logvar`).
   Fuse additively (τ adds). → ONE **λ-message**: mode `mo_g − mo_R` (from the density relay, which stays for the
   MODE), precision `τ_fused`.
2. **Anchor-gDNA-measurement stream:** seed `mg_own = where(struct_lock, prec_g, 0)` (the ÷M anchor gDNA count).
   Transport (density reframe + `σ²_transfer`), fuse additively. → **gdna_imp** (mode `mo_g`, precision the
   measurement stream). NOT the composition — a density lower-bound.
3. **Spliced-RNA-measurement stream:** the graft's spliced precision (`SP/SN`, already isolated in the relay),
   fused additively. → **rna_imp** (mode `mo_p/mo_n`, precision the spliced measurement). The θ tilt for AMBIG
   comes from `(c_p − c_n)/(c_p + c_n)`.

   On UNSTRANDED data `τ≈0` ⇒ the λ-message is ~off and this reduces to the measurement channels — **preserving
   the M5 unstranded win**. On STRANDED data the composition rides the ONE λ-message (no double-count) while the
   measurements stay independent — **removing the stranded regression**. Keep the density relay's full-precision
   (`pg/pp/pn`) for the MODE fusion (`_fuse`); track `mg/mp/mn` (measurement) + `τ` (composition) as ADDITIONAL
   accumulators in the SAME relay loop (no extra passes).
4. **`struct_lock` HARD OVERRIDE + the fusion/transport split** (handoff §7, still owed): the composition fusion
   weight is `τ_λ` (no `1/n` — the count cancels in the ratio; a REFINEMENT of the handoff's blanket `⊕1/n`);
   a `struct_lock` node ADOPTS its own belief (never an ∞ weight → nan). Add an interior-anchor no-nan test.
5. **A/B every step** (refit=0 AND refit=1, unstranded-capON called out); PASS iff no per-condition regression —
   stranded must return to ≤ ~0.03. If a regression appears, **DISSECT** (worst scenario → worst nodes → trace
   the message propagation to root cause; do NOT assume a theory flaw — the theory is sound). Then regenerate goldens.

### ✅ DONE — the three-stream relay + single-λ combine (Steps B/C)
Implemented (`bp_solver._unified_solve`): the density relay now carries, alongside the mode densities, a
**composition-τ** accumulator (seed `_ni.tau_lam`) → ONE **λ-message** (`lam_imp`, mode `mo_g − mo_R`, precision
the fused τ), and **measurement** accumulators `mg/mp/mn` (seed `struct_lock·prec_g` for the anchor gDNA; the
spliced adds at the graft) → the independent **gdna_imp/rna_imp**; the θ tilt (AMBIG) from `(c_p−c_n)/(c_p+c_n)`.
Tests green (383); the double-count is fixed (single-λ). A diagnostic env hook `RIGEL_S2T_OFF` zeros
`σ²_transfer`.

### ⛔ STILL FAILING — the stranded regression is NOT the double-count; it is a REFRAME/PROPAGATION defect
A/B: three-stream ≈ m5 (aggregate refit=1 0.0945→0.0998, stranded 0.0777→0.0778) — **the single-λ fix did not
move the stranded regression**, so it was never the double-count. **DISSECTED** (`gdna_gdna300_ss_0.99_capON`,
node 1909, an exon oracle f_g=0.985, mostly gDNA): its OWN message-free solve `fg_loc=0.953` is CORRECT but at
WEAK precision (`τ=1.6` — a mostly-gDNA region has little RNA, so the strand tilt barely constrains it); it then
receives a **strong, WRONG RNA+ measurement** (`cm_p=26.45`, mode `f_pos=0.718`, i.e. "72 % RNA+" for a ~1.5 %-RNA
node) that overrides the weak-but-correct own belief → f_g collapses to 0.51. **Root cause:** the **reframe `r`
AMPLIFIES a tiny high-RNA-fraction neighbour into a dominant RNA message at a high-mass gDNA node** (1909 mass
66,544 next to boundary 1910 mass 14, 64 % RNA ⇒ `r = ρ_tot(1909)/ρ_tot(1910)` enormous — the "r up to 10³ into
exons" pathology the code comments already lament). σ²_transfer damps the *precision* but the *mode* is already
wrong and the precision stays high enough to win. This defect is **orthogonal to the single-λ work** — it lives
in the density relay's reframe + the measurement propagation, and is present in m5 too (which is why it improves
the unstranded/capture arm — where messages SHOULD propagate — but corrupts high-gDNA stranded nodes where a
tiny mis-composed neighbour must NOT dominate).

**▶ NEXT lead:** cap/gate the reframe so a tiny neighbour's RNA fraction cannot become a high-precision message
at a high-mass, composition-mismatched node (mass-ratio or composition-mismatch gate on the reframe or on the
measurement propagation) — this is what unblocks the stranded gate. Dissect first (`RIGEL_S2T_OFF` + a per-node
dump; `scratchpad/dump_node.py` is the template) before changing the reframe.

### ✅ ROOT-CAUSED by the precision-audit workflow (`wf_fe3119fb`, 4 agents, machine-precision reproduction)
**The `σ²_transfer = 0` exemption on matched/graft reframes — correct for the MODE (enrichment cancels) — is
wrongly applied to the PRECISION.** A message crossing a large enrichment cliff arrives at FULL precision because
(a) the graft is σ²_transfer-exempt and (b) `Var(log r) ≈ 1/n` is negligible on well-counted nodes. At node 1909
the reframe is `r ≈ 407` (ρ_tot(exon)=34 vs ρ_tot(boundary)=0.08, dominated by the exon's gDNA), so the flanking
boundary's real 47-count spliced measurement crosses the 407× cliff undamped (precision ~17–26) at a
stretched-wrong mode (f_pos=0.718) → f_g→0.51. **Channel ablation: BOTH the RNA-measurement stream (`cm_p`) AND
the composition λ-stream (`c_tau`) are corrupted identically — only dropping BOTH recovers node 1909 to 0.91.**
Secondary: the measurement stream `mp` additively accumulates FOREIGN junctions' counts (`bp_solver.py:478`,
combine `:566`). Seeds are count-honest (`SP` mass ≤ integer count); the violation is entirely in the cross-cliff
TRANSPORT. **This IS the owner's insight** — a big enrichment cliff must COST precision; today it costs nothing.

**Constraint on the fix (verifier caution):** must be composition/mass-mismatch-aware, NOT a blanket own-count
cap — on unstranded AMBIG nodes `τ_own=0` and messages are the only information (the M5 win).

### ✅ DERIVED + IMPLEMENTED — the cross-cliff precision law `σ²_cliff = (log r)²` (`wf_c565d23a`)
The design workflow (3 derivations + 2 critiques + adjudicator, MC-validated) derived: the message's EXACT error
is `log(a_c^src/a_c^dst)` (the composition-SHARE mismatch, `a_c = ρ_c/ρ_tot`; all scale/enrichment/`M_dst`
cancel to machine precision). The imputation premise IS "shares match"; the honest precision prices how suspect
that is across the cliff. The law: **`τ_msg = 1/(1/p_src + (log r)²)`** — the maximally-honest scale-free
quadratic (→0 at `r=1`, grows with the cliff), grounded in the mass-ratio (counts) + eff-lengths, coefficient
exactly 1, NO magic number. It **replaces** the wrong `Var(log r)≈1/n` and **removes the graft/matched exemption
for the precision** (kept for the mode). Applied to EVERY stream. IMPLEMENTED (`bp_solver.py` `_relay`/
`_transport`); tests green (383, +1 xpass on the AMBIG-recovery test). Anchor node 1909: **0.510 → 0.910**
(oracle 0.985).

**A/B (`pass0_oracle_bench.py`):** the stranded regression is LARGELY RECOVERED and refit=0 is a NET WIN, but
`(log r)²` OVER-damps the biggest ENRICHMENT cliffs (verystrong) at refit=1:

| condition | refit=0 base→cliff | refit=1 base→cliff |
|---|---|---|
| aggregate | 0.1267 → **0.1092** ✅ | 0.1234 → 0.1283 ⬆ |
| stranded ss_0.99 | 0.0262 → 0.0391 (was 0.079) ✅mostly | 0.0293 → 0.0374 ✅mostly |
| unstr capON | 0.2543 → 0.2186 ✅ | 0.3874 → 0.3420 ✅ |
| **verystrong** | 0.2985 → 0.2779 ✅ | 0.1514 → **0.2356** ⬆ over-damped |

**The characterized limitation (the design workflow predicted it):** prior-free, `(log r)²` cannot distinguish
a pure-ENRICHMENT cliff (composition preserved — should cost NO precision, the M5 capture win) from a
COMPOSITION-mismatch cliff (should cost). `log r` is dominated by the mass ratio (enrichment), so `(log r)²`
over-attributes enrichment to composition drift on the extreme-capture (verystrong) arm. The **clean successor
is the composition-MISMATCH term** (workflow's "D3": damp by the actual `gap²` between the transported mode and a
reference), which is un-shippable prior-free (zero protection on AMBIG `τ_own=0`) but becomes the right object
**once the gDNA hyperprior (Phase 2) supplies the reference** `v_own` on AMBIG nodes. So: `(log r)²` is the
honest prior-free fix (recovers stranded, net win refit=0); the enrichment-vs-mismatch split is Phase 2.

**▶ NEXT decision:** (a) accept `(log r)²` + move to the gDNA-hyperprior (Phase 2), where the `gap²` mismatch
term replaces the `(log r)²` proxy and removes the verystrong over-damping; or (b) a prior-free `gap²` using the
self-solve `fg_loc`/`τ_λ` as the reference (damp only where the message disagrees with a confident own belief) —
untested, may recover verystrong without the hyperprior. RIGEL_S2T_OFF zeros the cliff term for A/B isolation.

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
