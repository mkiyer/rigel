# The mature-crossing gate — the RNA message's missing structure gate

**Status:** ⚠️ **DISMANTLED 2026-07-16 (commit `5e54fdc5`)** — the gate landed, then was removed in favour of a
pristine, gate-free relay ([`message_model_derivation.md`](message_model_derivation.md) §5). This document is
retained as the diagnosis + provenance of *why* the leak exists (§1 is still the authoritative statement of the
defect the nascent factory must counter); its **prescription (the gate itself) is retired**. The gate helped
7/7 cached conditions (ALL mwae 0.169 → 0.193, concentrated in introns), so the dismantle regresses the
mature-heavy suite until the honest `σ²_transfer` + nascent factory (`RNA − mature`) land. The `mrna_active_*`
mask it introduced stays computed in the statics (the nascent factory consumes it). Written 2026-07-16,
branch `calib-ambig-init-wip`.
**Supersedes:** [`splice_junction_absorption_fix.md`](archive/splice_junction_absorption_fix.md) (item 1) — that
document's **diagnosis is retained and credited**; its **§3.2 prescription and §5 step 1 are retired** (§1.2).
**Successor:** [`dof_pie_model_fix.md`](dof_pie_model_fix.md) (item 2) — next, but its stated justification
must be re-derived (§6.3).
**Reproduce:** `scripts/debug/_ablate_replay.py` (→ promote, Phase 1), `selfsolve_diag.py`,
`intron_message_trace.py`. Substrate: `~/Downloads/rigel_runs/ambig_dense_10mb/_selfsolve_cache` (7 conditions,
gdna300, **with oracle truth pools**).

> **The one-line statement.** The RNA message's emission gate is `free_s` — a **nascent**-continuity predicate —
> but its payload `fbp[lsrc]·sm` is the source's **total** unspliced RNA, which in an exon is ~95% **mature**.
> The gate and the payload describe different molecules. Gate the payload on the **mature**-crossing predicate
> that already exists and has never been wired up.

---

## 1. What is established

### 1.1 The defect (item 1's diagnosis — survives every attempt to break it)

Exon↔intron seams carry **exactly zero** crossing mature — 0 of 1,146, reproduced digit-for-digit on all 7
cached conditions — and the RNA message fires on **100%** of them, carrying the exon's full unspliced density.
Mature either splices out (913 junction seams) or the transcript is simply unspliced there (233). Either way it
never enters the boundary's unspliced channel. The message asserts it does.

The arithmetic, traced at B740 and confirmed independently by three parties:

```
exon 739 : 975 gDNA | 29,969 mature | 0 nascent  ->  rho_unspl = 17.794   (fg_loc 0.0296, correct)
B740     :  67 gDNA |      0 mature | 0 nascent  ->  rho_unspl =  0.670   (fp_loc 0.1071, correct)
intron741:  227 gDNA |      0 mature | 0 nascent  ->  rho_unspl =  0.563   (fg_loc 0.9937, correct)
```

Every node self-solves correctly. Then exon 739 tells B740 *"your RNA density is 18.3"* — into a node whose
**entire** unspliced density is 0.670. `f_pos` saturates at 0.984 and B740 relays that into intron 741 at
`pr = 1565`. **The source's belief is honest** (`fbp[739] = 0.9693` vs oracle 0.9685). Nothing upstream is
corrupt. The **channel** is structurally invalid.

The count-zero-information reading: cross-node imputation is licensed only where the imputed quantity is the
*same* quantity on both sides. Across an exon→intron seam it is not — and the 2-component `{RNA, gDNA}` pie
cannot separate mature from nascent to fix it. Even a **perfect** nascent estimator would send
`0.0475 × 18 ≈ 0.86` into a boundary of total capacity 0.67. **The honest message is silence.**

### 1.2 Why item 1's §3.2 must not ship — measured, not argued

| finding | evidence |
|---|---|
| the message-content edit is an **exact algebraic no-op** | `rho_mature_est = SPd/ESPd ≡ rho_mat_dst` ⇒ `rho_new ≡ rho_cur`. Four independent replays: `max|Δmode| = 4.263e-14` / `1.421e-14`; **0 of 579** floored modes differ — the shipped `max(rho, 1.0/erd)` floor ([bp_solver.py:452](../../src/rigel/calibration/bp_solver.py#L452)) strictly dominates the proposed `max(·,0)` clamp (42 edges do go negative; the floor absorbs every one) |
| it **does not deliver its own prediction** | three independent end-to-end implementations: intron mwae `0.1573 → 0.1540 / 0.1517 / 0.1515` vs §4's predicted **~0.0098** — ~26× short; **regresses** `ss0.50_none_capOFF` |
| its acceptance criterion is **unachievable by construction** | §5 step 2 requires "the RNA mode on the 10 traced introns must fall to ≤0". `Δmode = 0.00e+00` |
| its only real change (`pr` from the nascent count) is a **regression** | drops `n_mat` from the precision on all 679 B→exon edges (median 84% haircut, 23.5% fully silenced) while still carrying `n_mat/esp` in `rho` — degrading the motif-stranded spliced count, the **only strand-free evidence in the system** ([bp_solver.py:429-433](../../src/rigel/calibration/bp_solver.py#L429) flags it as priority #3) |
| the estimator it rests on is **noise** | handed an *oracle-perfect* source belief: corr(est, true nascent) = **0.311**, additive bias **+837 = 4.40× the true signal**, corr(est, bias) = **0.9949**. SNR at the oracle's real nascent level = **0.13** (nascent is 4.75% of an exon's unspliced RNA). On the `nrna_none` zero control it still calls >50% of the exon's RNA "nascent" on 15.2% of faces |
| §3.1's "the 233 non-junction seams are already correct" is **false** | the exon sends `fbp·sm` regardless of a spliced count. Oracle split: silencing junction-EI alone → 0.0890; **non-junction-EI alone → 0.1456** (baseline 0.1573). §3.2 can never touch them (`SPd = 0` ⇒ the subtraction is identically zero). On **7 of the 10 worst-damaged introns** the feeding seams have `SPd = 0` |
| §5 step 1's taxonomy is **strand-blind** | it builds from `coarse_type_array`, where *exon wins over intron across both strands*. **181 of 909** coarse-EXON regions (19.9%) carry an exon bit on one strand and an intron bit on the other, on the suite defined by opposite-strand overlap. Its verification criterion is also unmeetable — the 913/233 junction split is condition-dependent (705/441 under capture); only the totals **1,146 / 430 / 121** are structural |

Two of item 1's headline numbers **do not reproduce**. `rho_spliced/rho_mature = p25 0.985 / p50 1.185 /
p75 1.506 (n=614)`: measured **n=599, p25 0.751 / p50 0.962 / p75 1.024**, geomean 0.703 — invariant across
*nine* density-definition variants. **Provenance unknown.** It is load-bearing three times (§2's "the frames
are correct", §3.1's "the mature is estimable at the seam", §5's "do NOT tune `esp`"). *The frames are in fact
correct — but on brute-force enumeration against `_accumulator_reference.py`, not on that ratio.* And §1's
`0.0098` disagrees with its own cited tool: `selfsolve_diag` reports **0.006535** (§4.4 resolves this).

**The 20× is real but is not a ranking statistic.** It holds on **1 of 7** conditions (19.71 vs
1.20/1.59/1.65/1.86/4.49/8.95), because `nrna_none × ss0.99` makes true intron `f_g = 1.0000` — any RNA message
into an intron is damage *by construction*. **Under capture the attribution inverts**: gDNA does the damage
(exon 0.0219 → 0.1320, ×6.04), and exons are 86.7% of mass. *Not* cherry-picked: on the honest increment metric
the doc's cell ranks 4 of 7.

### 1.3 The predicate already exists

`NodeStatics.mrna_active_{pos,neg}` = `exon_s(left) & exon_s(right)`, built at
[node_geometry.py:372-375](../../src/rigel/calibration/node_geometry.py#L372), gathered onto the chain at
[node_geometry.py:459-460](../../src/rigel/calibration/node_geometry.py#L459), documented at
[node_geometry.py:359](../../src/rigel/calibration/node_geometry.py#L359) as *"the tighter **mature**-crossing
gate"*, and derived from [signature.py:201](../../src/rigel/calibration/signature.py#L201):

```python
def mrna_active_strands(signature):
    return (sig & BIT_EXON_POS) != 0, (sig & BIT_EXON_NEG) != 0
```

**Zero production consumers** (`grep -rn mrna_active src/` → definition sites only; referenced solely by
`tests/calibration/test_node_prior_classifier.py`). The fix wires up dead code; it adds no state and no constant.

**It matches the stated rule exactly — 0 of 256 signature pairs disagree**, and 0 subsumption violations
(`mrna_active_s ⇒ free_s`, already asserted at
[test_node_prior_classifier.py:77](../../tests/calibration/test_node_prior_classifier.py#L77)):

> *mature passes iff the **same-stranded** exon bit is set on **both** sides. Intron bits also being set does
> not block. `EX+EX- | EX+EX-` passes on **both** strands.*

The per-strand property is the whole point on this suite. For `{EX+, IN−} | {IN+, IN−}`: the `+` strand is
exon→intron ⇒ `mrna_pos = True & False = False` ⇒ blocked; the `−` strand is intron→intron ⇒
`mrna_neg = False & False = False` ⇒ the source carries no mature, so its payload is genuinely pure nascent and
**flows freely**. A strand-blind predicate blocks both.

---

## 2. The design

> **Nascent originates in introns and at intron↔exon seams. It must propagate intron→exon from BOTH sides.
> An exon must never manufacture nascent INTO an intron.**

That is a statement about the **source's composition**, not about the edge — which is why the gate is
**asymmetric**:

```python
send_s = mrna_active_s[dst] or not mrna_active_s[src]
```

* **mature-active source → mature-incapable destination** ⇒ **blocked**. The exon's payload is
  mature-contaminated and the destination cannot host mature. These are exactly the 1,146 seams.
* **mature-active source → mature-active destination** ⇒ sent. Contiguous exon body; mature genuinely crosses
  unspliced (44,774 units measured).
* **non-mature-active source** (an intron, an intron-seam — *where nascent actually lives*) ⇒ **always sent**,
  in every direction, including into exons.

Traced on the real chain (reproduced digit-for-digit by two independent refuters): `739→740` **blocked**;
`740→741`, `741→740`, `740→739` all **allowed**. Nascent reaches exons from both sides; the exon's mature
never reaches an intron.

**Silence emerges — it is not a second Boolean.** On a gated edge `n_nasc = 0`, and `n_mat = SPs[lsrc] = 0`
(measured: **0 of 5,543** gated edges carry a live `n_mat`), so `n_src = 0` ⇒ `pr = 0/(0·vb+1) = 0` **exactly**
— [bp_solver.py:350](../../src/rigel/calibration/bp_solver.py#L350) is deliberately written in that form
("*so `n_src=0` gives `pr=0` exactly … zero density is not a measurement*"). **No new constant is introduced.**

### 2.1 The honest part: this is not the metric-optimal arm

Three arms were run head-to-head from one pinned snapshot. **D loses the benchmark:**

| arm | predicate | intron (primary) | verdict |
|---|---|---|---|
| **A** baseline | — | 0.1573 | control |
| **B** and-symmetric | `mrp[src] and mrp[dst]` | **0.0093** | beats D on intron **7/7**, on ALL in all 4 capOFF cells; loses the 3 capON ALL cells |
| **E** equality-symmetric | `mrp[dst] == mrp[src]` | ~0.0456 | **ties** D on intron 7/7, **beats** D on ALL in **7/7** — D never wins |
| **D** asymmetric ← **ship** | `mrp[dst] or not mrp[src]` | **0.0456** | improves intron **7/7**, ALL **6/7**; **beats B on all three capture-ON cells** |

**D ships because B and E are structurally wrong, not because it wins.**
* **B** deletes every RNA message touching an intron — including intron↔intron nascent relay. Measured:
  **533,708 / 166,068 / 148,446** units of true nascent relay per `nrna_present` condition (and **exactly 0**
  on all 4 `nrna_none`, as it must be) that B destroys and D preserves.
* **E** blocks whenever `mrp` *changes* across the edge — so on the `boundary→exon` hop
  (`mrp[B]=False, mrp[exon]=True`) it **blocks intron→exon nascent**, the exact direction the model requires.

D's arithmetic cost is real and should be stated: D captures **~74%** of B's intron gain
(0.1573 → 0.0456 vs 0.0093). D and B differ on **29,272 of 37,479** directed edges; D deliberately lets
seam/intron nodes relay, and their *local* beliefs still carry some exon mature. That residual is the price of
keeping the nascent channel alive.

> **⚠ "Asymmetric strictly dominates symmetric" was never measured and is refuted.** Prior work compared its own
> asymmetric numbers against *another document's reported* symmetric numbers. Do not repeat that claim.

### 2.2 D's safety margin over B — the false-negative direction

`mrna_active` is a **biased approximation**, not an exact invariant (§4.3). Its errors point in opposite
directions for the two arms:

* **False negatives** (`mrna_active = False` where mature truly crosses): D reads the predicate **at the
  source**, so an FN at the source means D **sends anyway**. Harmless under D; only **B** would silence them.
* **False positives** (`mrna_active = True` where no single transcript spans the seam — the union-over-
  transcripts construction): D's actual risk direction, an unwarranted block. Empirically dominated — D still
  improves intron 7/7. **Watch the `−` strand** (all 36 known FPs are neg-strand).

---

## 3. The exact diff

### 3a. Bind the mask — [bp_solver.py:204](../../src/rigel/calibration/bp_solver.py#L204)

```python
# BEFORE
    fp, fn = statics.free_pos, statics.free_neg
# AFTER
    fp, fn = statics.free_pos, statics.free_neg
    mrp, mrn = statics.mrna_active_pos, statics.mrna_active_neg  # mature-crossing gate (see `_scan`)
```

### 3b. The gate, `+` strand — [bp_solver.py:440-442](../../src/rigel/calibration/bp_solver.py#L440)

```python
# BEFORE
                n_nasc = (
                    fbp[lsrc] * sm
                )  # source total unspliced RNA count (nascent + exon-body mature)
# AFTER
                # MATURE-CROSSING GATE (`mrna_active_s` = exon_s(left) & exon_s(right)). The unspliced payload
                # is NASCENT. Nascent originates in introns and at intron<->exon seams and may propagate
                # intron->exon from BOTH sides; an exon must never manufacture nascent INTO an intron. So a
                # mature-active source sends only into a mature-active destination; a NON-mature-active source
                # (an intron / intron-seam -- where nascent genuinely lives) always sends. ASYMMETRIC BY DESIGN:
                # the symmetric AND over-blocks the intron->exon direction the model requires. Intron bits never
                # block -- EX+EX-|EX+EX- passes on BOTH strands. `n_mat` (the spliced MEASUREMENT) is untouched.
                n_nasc = (fbp[lsrc] * sm) if (mrp[i] or not mrp[lsrc]) else 0.0
```

### 3c. The `emit_n` mirror, `−` strand — [bp_solver.py:466](../../src/rigel/calibration/bp_solver.py#L466)

```python
# BEFORE
                n_nasc = fbn[lsrc] * sm
# AFTER
                n_nasc = (fbn[lsrc] * sm) if (mrn[i] or not mrn[lsrc]) else 0.0
```

**The mirror is not optional** — the gate is per-strand (`mrp` gates `+`, `mrn` gates `−`). An earlier
"2.79× mass inflation" error came precisely from summing both strands' truth pools onto every per-strand edge.

### 3d. Do **not** touch `emit_p` / `emit_n` (lines 400-401) — measured, exact

Static enumeration over all 7 conditions (every term is static, so this is exact, not sampled):

| | count |
|---|---|
| emitting edges | 37,485 |
| gated edges (D) | 5,543 |
| **gated edges carrying a live `n_mat`** | **0** |
| non-finite `mo` on gated edges | 0 |

Gating `emit_p` would be behaviourally equivalent **but would also silence the `n_mat` measurement path**.
**Gate the payload, not the emission.**

### 3e. The `pr = 0` combine hazard — benign, proven structurally

On a gated edge `rho = −rho_mat_dst ≤ 0` ⇒ `max(rho, 1.0/erd) = 1.0/erd` ⇒ `mo = log(1/md)`, and
`md ≥ 1e-9` ⇒ `mo ≤ log(1e9) ≈ 20.7` — **always finite**. So in `_comb`, `ap·am = 0·mo = 0`: never `0·inf`,
never NaN. **No guard is required for safety.**

But `pr = 0` is *not* inert in the running belief ([bp_solver.py:456-460](../../src/rigel/calibration/bp_solver.py#L456)):
`fbp[i] = exp(log(fp_loc[i]))` round-trips (~1 ulp drift), and `vbp[i] = 1/pp_loc[i]` **clobbers an exact
`vp_loc = 0` to `1e-9`** (measured: `vp_loc == 0` on **10,555 / 23,779** node-conditions). It propagates via
`pr = n_src/(n_src·vbp[lsrc]+1)` — at `n_src ≈ 1e4`, a 0.001% perturbation. **Measured aggregate effect: 0.0.**

```python
                pr = n_src / (n_src * vbp[lsrc] + 1.0)
                if pr > 0.0:  # zero precision is not a message — leave the running belief untouched
                    amp[i], app[i] = mo, pr
                    pt = pp_loc[i] + pr
                    fbp[i] = math.exp((pp_loc[i] * lfp_loc[i] + pr * mo) / pt)
                    vbp[i] = 1.0 / pt
```

**This is hygiene, not a defect fix** — ship it as its own bit-identical phase so the gate's delta is purely the
gate. RNA branches only; the gDNA branch carries its own bit-identity risk — leave it, note it.

### 3f. ⚠ Do **not** delete the mature-absorption term

The gate blocks the edge the absorption fires on, which invites "the gate replaces absorption, delete it."
**Measured: no.**

| absorption-firing edges (`SPd[i] > _EPS`) | 5,172 |
|---|---|
| …gated (subsumed) | 4,851 (93.8%) |
| …**still fire** under D | **321** |
| mature mass on survivors | **92,472 units** |

The survivors are exon → exon↔exon-junction edges (`mrp[dst]=True` ⇒ not gated; `SPd>0` ⇒ absorption fires).
**The term stays.** See §5 D3 for the question this opens.

---

## 4. Phases

Each phase is independently verifiable with an explicit falsifier. **Every phase runs against a pinned
snapshot, never the live tree** (§5 D5).

> **PROGRESS (2026-07-16):** Phase 0 ✓ (pinned `d73260db`), Phase 1 ✓ (`scripts/debug/ablate_replay.py` +
> `scripts/debug/gate_edge_census.py`; census reproduces `37485/5543/0/5172/4851/321` exactly; ablation
> reproduces the §1 table + fidelity `0.000e+00` on 7/7), Phase 2 ✓ (the `pr > 0.0` guard, both RNA branches;
> A/B on 7/7 conditions **bit-identical**, `max|new−old| = 0.000e+00`), Phase 3 ✓ (three tests in
> `tests/calibration/test_bp_solver.py`: `test_exon_does_not_manufacture_nascent_into_intron` **xfail(strict)**
> — genuinely FAILS today under `--runxfail` at the exact 1926.02 damage precision, flips to pass in Phase 4;
> `test_intron_relays_nascent_into_exon_both_directions` and `test_mrna_active_matches_same_strand_exon_rule`
> pass now and pin the asymmetry + the 256-pair predicate). **Phase 4 ✓ THE GATE IS IN** (`bp_solver.py` §3a-c:
> `mrp/mrn` bound + `n_nasc` gated in both RNA branches; test 1 flipped to pass, xfail removed,
> `test_mature_absorption_lowers_nascent_message_into_junction` retired). Verified: **intergenic `f_g` = 0.0000
> on all 7** (no manufactured gDNA — the falsifier), the gated intron column reproduces the predicted **D**
> digit-for-digit on all 7 (`0.1573→0.0456`, `0.0957→0.0401`, `0.0715→0.0252`, `0.3322→0.2610`, `0.4520→0.3129`,
> `0.2598→0.1766`, `0.3450→0.2290`), fidelity `0.000e+00` on 7/7, full suite **21 failed (all pre-existing
> golden debt) / 1077 passed** — zero new non-golden failures. Gated `bp_solver.py` hash =
> `c24282644fe2b62f389ddb43c8543614`. **Next: D4** (strand-aware `_spliced_faces` routing — the user's directive,
> the immediate follow-up), then **Phase 6** (production basis, fitted prior, the 13M-leak cell) → **Phase 7**
> (one golden regen). **Phase 5** (the 4-arm A/B) is documentation-only now that D1 is confirmed. Note: the
> `node_densities`/`NodeDensities` removal considered during cleanup was **cancelled** — it is a tested reference
> (2 unit tests + a use + `pass_trace.py`), not dead.

### Phase 0 — Pin the baseline
`bp_solver.py` mutated **mid-session for three separate investigators** (`95341e7d` → `d73260db` → `6b904364`).
Identical code measured **0.0433 ungated vs 0.0138 gated** on `ss0.99_none_capOFF` — a **3× baseline shift**;
two workstreams unknowingly measured different trees. **Verified at time of writing:**
`md5(bp_solver.py) = d73260dbe019fbc8e6abd836a425365b`, `grep -c mrna_active = 0` ⇒ un-gated baseline, 213
calibration tests green. Record `md5` before *and* after every measurement; discard any run whose hashes differ.

### Phase 1 — Instrument first (zero `src/` change)
Promote the census → `scripts/debug/gate_edge_census.py`; promote `_ablate_replay.py` →
`scripts/debug/ablate_replay.py`, fixing its verbatim dead expression at line 104
(`False if False else True` → `True`). Delete the other ~13 uncommitted `_adv_*.py` scratch scripts.
* **VERIFY:** census reproduces `emit=37485, gated=5543, gated_live_nmat=0, absorb=5172 / gated=4851 /
  survives=321`. `ablate_replay` prints `max|both − shipped| = 0.000e+00` on 7/7 — the tripwire.
* **FALSIFIER:** `gated_live_nmat > 0` ⇒ §3d is wrong and §3e's `pr=0` proof collapses. **Stop and re-plan.**

### Phase 2 — The `pr > 0.0` guard, alone
* **BIT-IDENTICAL: required.** `max|guard − baseline| = 0.0`, every class × all 7 conditions.
* **FALSIFIER:** *any* non-zero delta. **Prior work never tested this** — it measured gate+guard vs gate-only,
  never guard-alone vs baseline. A non-zero delta means a behaviour change masquerading as hygiene: **drop the
  guard**, proceed with §3a-c alone. The gate does not need it.

### Phase 3 — A test that FAILS on today's code
Reuse the existing `_mature_exon_chain` fixture in `tests/calibration/test_bp_solver.py`:
1. **`test_exon_does_not_manufacture_nascent_into_intron`** — backward scan, `+RNA` into `B1` (src = exon `R1`,
   `mrp=True`; dst = intron-side junction, `mrp=False`). Assert precision **exactly 0**. **Must FAIL today**
   (today asserts `> 0.0`, live at [test_bp_solver.py:743](../../tests/calibration/test_bp_solver.py#L743)).
2. **`test_intron_relays_nascent_into_exon_both_directions`** — assert `intron→boundary` and `boundary→exon`
   precisions stay `> 0.0`. **Passes today and must keep passing.** This pins the asymmetry and is the sole
   codified reason B and E are excluded.
3. **`test_mrna_active_matches_same_strand_exon_rule`** — the 256-pair table (pure, cheap).
* **FALSIFIER:** test 1 passes today ⇒ the mechanism is not what we think. Test 2 fails ⇒ the fixture is wrong.

### Phase 4 — The change (§3a-c)
* **NOT bit-identical.** This is the change.
* **VERIFY:** `intergenic` stays **0.0000** on all 7 (measured under A, B *and* D — no gDNA manufactured).
  Phase-3 tests 1+2 pass. Exactly one pre-existing test flips (§4.2).
* **FALSIFIER:** intergenic ≠ 0 ⇒ the gate manufactures gDNA. **Revert.**

### Phase 5 — The 4-arm A/B, head-to-head on ONE snapshot
Arms A / B / D / E (§2.1). **Expected D intron, 7/7** (reproduced digit-for-digit by three parties):
`0.1573→0.0456, 0.0957→0.0401, 0.0715→0.0252, 0.3322→0.2610, 0.4520→0.3129, 0.2598→0.1766, 0.3450→0.2290`.
ALL improves 6/7; `ss0.50_none_capON` regresses **+0.0018** (inside documented ~1% suite noise).
* **FALSIFIER:** D's intron column does not reproduce ⇒ contaminated snapshot (check md5) or mis-transcribed
  predicate.
* **This phase does not select the arm.** It exists to **document honestly that D loses to E and B** on this
  metric. See §5 D1.

### Phase 6 — Production basis ✅ RAN (2026-07-16, bug-finding lens — not an accuracy verdict)
Production `calibrate` (fitted prior) vs oracle over the 24-condition `_calib_cache`, gate ON vs OFF
(un-gated arm = `mrna_active_*` forced all-True). ⚠ **Harness lesson:** `rigel.calibration.calibrate` (attribute)
resolves to the FUNCTION (re-exported in `__init__`), so `import … as C; C.build_node_statics = …` patches nothing
and yields a spurious ZERO-difference; patch `sys.modules['rigel.calibration.calibrate']`. With the correct patch
the gate moves 754 regions (max Δf_g 0.92). Findings:
* **The gate reduces the gDNA→RNA leak by −7.0%** (16.86M → 15.68M fragments), mean `mwae_fg` 0.175 → 0.168,
  **12/24 improved, 8/24 worse**.
* **capture-OFF: clear win** (leak −100k…−315k per cell). **capture-ON + unstranded: slight regression** — the
  known capture-enrichment regime (the 13M-leak cell), which the message layer cannot fix without the prior work.
* **The regression that matters (a REAL, predicted effect, not a bug):** on `nrna_present × unstranded` the gate
  pushes introns toward gDNA (they lose the exon-sourced RNA message and, at κ≈½, cannot deconvolve their own
  nascent). This is `boundary_spliced_channel_design.md` §4.4's predicted muting and points at the NEXT work
  (the nascent-support channels: D4 routing + priority #3 spliced measurement), not a gate defect.
* **Node-level (prior-free, 7 oracle conditions, `scripts/debug/node_error_report.py`):** intron mwae
  0.2035 → 0.1263 pooled; the worst RESIDUALS are capture-enriched exons (true gDNA, called RNA — the known
  leak, gate-neutral); the gate-WORSENED nodes are the nrna_present introns above + a few short exons (the
  Issue-#5 `eff_rna` interaction, §6.1).
* **Not an accuracy verdict** — mid-implementation of the reference/prior; this is for finding bugs. A real
  production-accuracy pass waits until the reference + prior work lands.

### Phase 7 — One golden regen
Exactly one `pytest tests/ --update-golden`, **after** the gate lands. Not before, not per-phase.

### 4.2 Blast radius — measured

* **Unit tests: exactly ONE new failure, zero fixed** (21→22 failed, 1075→1074 passed):
  `test_bp_solver.py::test_mature_absorption_lowers_nascent_message_into_junction`. Line 743 asserts
  `app[b1] > 0.0` on the backward exon→junction `+RNA` edge — **precisely the edge D blocks**. It characterises
  the mechanism being removed: **retire it**, replaced by Phase-3 test 1 (same fixture, same edge, inverted
  assertion). Its companions (`test_mature_measurement_recovers_exon_rna`,
  `test_mature_measurement_disagreement_silenced`) test the `n_mat` measurement path, untouched by §3d — they
  must keep passing.
* **Goldens: the branch already owes a regen.** All **21** golden tests fail *before any change*
  (`test_golden_output.py` has exactly 21 tests ⇒ 100%). The gate adds **zero** new golden failures but **moves
  values** (3 of the first 10 max-abs-diffs move). **Movement is MIXED — one moves distinctly away
  (4.317 → 8.804).** Do not carry a "the gate improves the goldens" narrative; judge the regen on its own diff.
* **Flag: none. Land unconditionally.** `CalibrationConfig` has exactly 10 fields, **all numeric, zero
  booleans**; `grep ': bool' src/rigel/config.py` yields only BAM-scan semantics and an output toggle — **no
  solver-variant flag precedent**. The productionizing pass deliberately *removed* `RIGEL_MSG_MODE` /
  `RIGEL_EM_WARMSTART`.

### 4.3 Do **not** write an exact-invariant test

The review's *"TRUTH ⊆ `mrna_active_s`, zero false negatives (1,317/1,317)"* **failed to reproduce on two
independent bases**: **200 (pos) / 306 (neg)** violations, with **97.9% / 99.6%** of boundary mature mass on
`~mrna_active` boundaries; `TRUTH_n` sums to **867, never 1,317**. **`mrna_active` is a biased approximation.**
An exact assert would encode a falsehood.

**But the measurement itself is contested and cannot currently settle it** (§6.5): the oracle's truth pools are
**read-strand** keyed while `mrna_active_s` is **transcript-strand** keyed (pool split 22,309/22,456 at
κ=0.5005 vs 13,901/30,873 at κ=0.0099). Anyone comparing them reads a convention artifact as a violation.
**Resolving this needs a transcript-strand-keyed truth pool.** Until then: no invariant test; rely on §2.2's
error-direction argument, which is robust to the ambiguity because D reads the predicate at the *source*.

### 4.4 The metric basis — `0.0098`, not `0.006535`

Fully diagnosed: they are **different estimators**. `selfsolve_diag --stage self` runs `init_beliefs` alone;
`_ablate_replay`'s msgfree is phase-D `_solve_nodes_logodds_all` with message precisions forced to 0.
`selfsolve_diag`'s `0.006535` is **flattered by a truth-coincidence** — on `ss0.99_none_capOFF` intron
`true_fg = 1.000000` and AMBIG nodes lock at `{0,0,1}`, scoring free; on `nrna_present`, where intron
`true_fg = 0.684312`, the same metric degrades **8.42×** to `0.055010`. That is skill-free.
**Quote `0.0098` (`_ablate_replay`) and rename the column `mwae_init`.** Note "floor" is loose — B reaches
`0.0093`, *below* it; it is a reference level, not a bound.

---

## 5. User decisions — ALL RESOLVED (2026-07-16)

**D1 — Ship the structurally-correct arm (D). ✅ CONFIRMED.** The user: *"we must allow intron→exon nascent RNA
and nascent RNA relay. Our plan is absolutely sound. The metrics being tested are not actually relevant."* D
ships on the physics, not the benchmark. (For the record: E ties D on intron 7/7 and beats it on ALL 7/7, and B
reaches intron 0.0093 — but E blocks intron→exon nascent and B deletes 533,708 units/condition of real intron
relay, so both are wrong regardless of the metric.) Phase 5's 4-arm A/B is therefore **documentation only, not a
decision gate** — the arm is chosen.

**D2 — Issue #5's `eff_rna ≥ 1` gate: ✅ ACCEPTED on the structural argument.** *(This is a SEPARATE fix from
the gate, ships after — see §6.1.)* The user: *"An effective length is a discrete value… E is the discrete
count of admissible start positions. There must be at least ONE start position."* So the `≥ 1` is the
**quantity's own unit**, not a tuned threshold — a fragment needs ≥1 admissible start position to exist. **The
benchmark objection ("threshold 100.0 scores better") is explicitly overruled** as benchmark-chasing and is
struck from the design; do not cite it. *(Terminology for readers: "Issue #5" = §6.1 — 191 short regions whose
RNA effective length is a fraction of one start position yet carry gDNA reads that the strand tilt misreads as
RNA, fabricating a ρ_rna message. "The threshold" gates the RNA imputation on that region's `eff_rna`.)*

**D3 — Absorption always fires on splice-junction mass. ✅ CONFIRMED, and the gate already respects it.** The
user: *"Whenever there is splice junction mass at a boundary, the splice junction MUST absorb mature RNA. There
are many boundaries where the left region is (exon AND intron) → (exon) — mature contiguous on both sides PLUS a
splice junction. We must account for contiguous mature AND splice-junction-absorbed mature."* The gate touches
**only `n_nasc`** and leaves the absorption term (`− SPd/E_spl_dst`) and the measurement (`n_mat`) untouched, so
absorption keeps firing on every boundary with spliced mass — including into a mature-active exon (the 321
surviving edges). Verified for the `(exon+ AND intron+) → exon+` case: the gate ALLOWS the contiguous-mature
edges both ways (`send_p = True`) while the splice-junction mass rides the separate absorption term. **No change
to the gate; the two mechanisms coexist by construction.**

**D4 — Mature/splice routing MUST be strand-aware. ✅ DONE (2026-07-16).** Fixed in `node_geometry._spliced_faces`:
`any_exon_l/r` (strand-blind) → per-strand `exon_pos_l/r` / `exon_neg_l/r`, so a `TS_POS` junction routes mature
only to `BIT_EXON_POS` flanks and `TS_NEG` only to `BIT_EXON_NEG`. Test-first:
`test_spliced_routing_is_strand_aware_at_ambig_seams` (the user's 4-transcript locus — the −junction at 10800 and
+junction at 11000) FAILED on the strand-blind code (routed 60/80 units onto the opposite-strand exon flank),
passes after. **Metric-neutral on the 24-cond suite** (leak Δ = +5 fragments — only 93/6694 = 1.4% of faces
mis-routed, and the spliced channel is heavily down-weighted today) but a genuine correctness fix that bites once
the spliced measurement carries real weight (priority #3). *Original overview of the defect:* The
user: *"mature RNA is stranded, and splice junctions are stranded. A splice junction can only absorb mature RNA
on its strand. Mature RNA can only travel contiguously on one strand. The routing must respect strand."* The
**gate mask is already per-strand** (`mrna_active_strands` = `sig & BIT_EXON_POS` / `& BIT_EXON_NEG`), so the
gate is strand-correct. The defect is in `node_geometry._spliced_faces`
([node_geometry.py:163-164](../../src/rigel/calibration/node_geometry.py#L163)): `any_exon_l/r = (sig &
(BIT_EXON_POS|BIT_EXON_NEG)) != 0` is **strand-blind**, so a `TS_POS` junction routes its mature mass to a flank
that has *only* a `BIT_EXON_NEG` bit (93/6694 = 1.4% of spliced faces). Fix: `any_exon` becomes per-strand —
route a `TS_POS` junction's mass only to a flank carrying `BIT_EXON_POS`. This feeds the `n_mat` measurement and
the `SPd` absorption (D3's mechanism), so it makes D3 strand-correct. **Schedule: the immediate follow-up after
the gate**, its own diff + its own A/B (keep attribution clean).

**D5 — Concurrency. ✅ HANDLED (to the extent possible).** No active writer in this session; `bp_solver.py` hash
verified stable across reads at each phase (pinned per phase). A file lock against *other* sessions is not
possible from here — the mitigation is the per-phase md5 pin + `ablate_replay`'s `0.000e+00` fidelity tripwire,
which detects any drift. **Ensure no other Claude session is open on this repo during the remaining phases.**

*No new constant is introduced by the gate itself* — it is a Boolean over signature bits that already exist.
**The magic-number directive is clean for this change.**

---

## 6. Open issues — resolved, and still open

### 6.1 Issue #5 — short regions fabricating an RNA-density message *(the `_EPS`-floor hypothesis is refuted; the real bug is elsewhere)*

**Plain statement (D2's clarification).** An effective length is the expected count of admissible fragment
*start positions* in a region — a discrete quantity: base pairs, whole start sites. A region shorter than the
fragment length has an RNA effective length that is a *fraction of one* start position (`0 < eff_rna < 0.25`).
Such a region cannot host an RNA fragment, so its true RNA mass is zero — **but it still catches gDNA reads**,
and the strand tilt of those gDNA reads gets misread as a sliver of RNA, which divided by the fractional
`eff_rna` explodes into a huge fabricated RNA-density message to its neighbours. **The fix gates the RNA
imputation on `eff_rna ≥ 1`** — a region must admit at least ONE RNA start position to send an RNA density.
Per D2 this is the quantity's own unit, not a tuned threshold; the benchmark objection is overruled.

> *"When region effective length is less than the fragment length, fragments cannot be contained… We should get
> rid of the `_EPS` floor, right?"*

**No — and it does not matter either way.** Removing the floor is a **bit-identical no-op**
(`delta = 0.000e+00`, 0 NaN, 7/7, in two independent tree states) because `node_sweep` **re-floors with its own
identical `_EPS`**: `eg = EGs[lsrc] if EGs[lsrc] > _EPS else _EPS`
([bp_solver.py:405](../../src/rigel/calibration/bp_solver.py#L405)/[:439](../../src/rigel/calibration/bp_solver.py#L439)),
and both constants are exactly `1.0e-9`. **The production message path never sees the geometry floor.** The
floor *is* load-bearing for `node_densities()` (211 NaN without it) — which has **zero production consumers**.
It is a division guard for a dead function.

The premise holds: the invariant `eff == 0 ⇒ mass == 0` is verified **210/210, zero violations**. And the
frame-mismatch hypothesis is **refuted**: both `accumulator.cpp:165-172` (`if (all_same)`) and
`_accumulator_reference.py:211-215` (`if len(regions_touched)==1`) credit `region.contained` **only** on the
single-region path — an encompassing fragment goes to **boundaries**, never to `region.contained`. Numerator and
denominator count the same event, and `rho_gdna` is flat (~0.6) across every length bucket from L=100 to L=3Mb,
which it could not be if the frames disagreed.

**The real bug is RNA-channel-only and oracle-verified.** **191 region nodes** with `L` < the RNA fragment
length have `0 < eff_rna < 0.25` (p50 **0.051**) yet carry **~5,100 fragments**. The oracle says their true
`f_gdna = 1.0000` (mean *and* min) and true RNA mass = **0 EXACTLY** — but the strand solve assigns
`f_rna` p50 0.04–0.08 from the strand tilt of **gDNA** reads, so `rho_rna = f_rna·M/eff_rna` explodes to
p50 **21.8** / max **387** versus 0.13–0.29 for normal regions. A ~100–170× fabricated RNA density injected
into neighbours.

**Orthogonal to the gate — measured, one snapshot, 4 arms:** base 0.217606 → eff-only **0.210065** →
gate-only **0.203866** → both **0.200704**. Each improves 7/7 alone; they compose; disjoint node sets.
**"#5 must precede the gate" is refuted.** **Ship as its own fix**, after the gate. The threshold basis is now
settled (D2: `eff_rna ≥ 1`, the ≥1-start-position unit — no magic number); the `100.0-scores-better` benchmark
argument is struck.

### 6.2 Issue #7 (`spliced_side_eff_length`) — **real, confirmed, broader than reported, and inert**

The geometry **is** wrong — two from-scratch brute forces against `_accumulator_reference.py` reproduce the
headline exactly (`R=100/ℓ=200`: shipped **25.000** vs truth **50.000**; `R=50/ℓ=200`: **6.250** vs **25.000**),
up to **199.6×** on the real suite. The crux resolves **against the shipped code**: the half-triangle
`E[min(ℓ,R)²/2ℓ]` is correct **only when the mature transcript terminates at the flank region's far edge**;
whenever it *continues*, the accumulator's interior-slice ½ rule makes the face's expected mass exactly
`E[min(ℓ,R)]/2` = `boundary_side_eff_length`. And "the exon continues past the region" **misses the dominant
case**: if the flank's far edge is *itself* a junction (a short internal exon), the fragment spills across via a
second splice and the face still gets `R/(2ℓ)`. **Scope is ~3× the original claim: 277/1053 junction faces
(26.3%, not 8.2%), 15.7% of junction spliced mass (not 3.6%), mean 2.29×.**

**But it is functionally inert today**: an A/B through the real prior-free `node_sweep` moves mwae by
**≤0.0001** on 7/7 (≤0.0005 ungated), while deleting the *entire* spliced channel moves it by 0.008–0.014.
**"The gate raises this bug's weight, fix it first" is refuted by direct measurement in both tree states.**

> ⚠ Two tempting justifications are **dead**. The "the precision cap absorbs a 200× error" mechanism is
> **refuted** (×5.627 on the fix's own faces → ≤0.0007; the *same* magnitude on all boundary faces → **+0.0507**
> — ~70× apart at identical cap and log compression: **the faces are inert, not absorbing**). The "latent
> landmine when priority #3 un-gags the channel" prediction rested entirely on that refuted mechanism. And the
> fix **does not repair even its own targets** (±0.0007 restricted).

**Call: NOT BLOCKING. The honest case is correctness + docstring hygiene, not accuracy.** ✅ **Docstring FIXED
(2026-07-16):** [effective_length.py:86](../../src/rigel/calibration/effective_length.py#L86) no longer calls the
2×–199× error "cosmetic"; it now states the half-triangle is correct only in the terminal-exon world, is 2–199×
low under continuation (~26 % of junction faces), is inert today (A/B ≤1e-4), and becomes load-bearing once
priority #3 un-gags the spliced channel. **The geometry fix itself is DEFERRED to priority #3** (it only matters
when that channel carries weight).

### 6.3 Item 2 (`dof_pie_model_fix.md`) — its justification must be re-derived

Item 2's stated case is *"it bounds item 1's residual"* — **conditional on §3.2 shipping**. Under §3.2 the
residual is 0.1517 and that case is sound. Under the gate the capOFF residual is largely gone and the capON
remainder is the **gDNA** channel, not composition. **Item 2 must make its own case on the post-gate residual.**
It very likely still earns it — the pie incoherence is real, measured, and independent (71.7% of solvable nodes
carry `fbg+fbp+fbn ≠ 1`, max 93.35; B740 relays a 52.70 "composition"). There is also a live, never-run
disagreement to settle afterwards: item 2 §3 predicts it alone will *not* fix the bulk; the pie-bound argument
predicts intron ≲ 0.02. Cheap and decisive — but run it **after** the gate, unconfounded.

Two cheap instruments from item 2's territory cost no behaviour and can land any time: `assert n_src <= sm` on
every edge (B740 violates it 52×), and the pie probe `|fbg+fbp+fbn − 1| < 1e-9` as `xfail(strict=True)`.

### 6.4 The count term — derived, unclaimed, ship separately

`pr = n_src/(n_src·vb+1)` with `n_src = fbp[lsrc]·sm` charges **Poisson sampling** precision to a
**deconvolved sub-count** — double-booking belief as evidence. Since `M ~ Poisson(ρ_tot·E)` and
`ρ_c = f_c·ρ_tot`, the correct term is `Var(log ρ_c) = Var(log f_c) + 1/M_observed` — i.e. **`1/sm`, not
`1/(fbp·sm)`**. That is why B740 reaches `pr = 27,208` on **67 fragments** and outguns local evidence 72,091×
(`pp_loc[740] = 0.3774`). One line, derived, no constant, independently A/B-able. **Do not** fold `n_mat` into a
shared `1/n` on B→exon edges — that is the Fenton–Wilkinson split, priority #3.

### 6.5 What is NOT established

* **The provenance of `1.185 / n=614`.** Nine density-definition variants do not reach it.
* **`mrna_active`'s true false-negative rate.** The oracle **cannot** settle it as keyed (§4.3) — read-strand
  truth vs transcript-strand predicate. **Needs a transcript-strand-keyed truth pool.** Until then the
  200/306 violations are *not* proof the predicate is wrong, and the 1,317/1,317 invariant is *not* proof it is
  right. **Neither number may be quoted.**
* **What causes the capON regression.** Exons are 86.7% of capON mass with the gDNA channel dominant
  (gdna_only 0.6234 vs msgfree 0.3154) — it may not be the RNA gate at all.
* **The 36 `mrna_active` false positives** (all neg-strand, 41,062 units). `spliced_n_{pos,neg}` **cannot**
  close them: 0/36 have dst-face spliced on the emitting strand, because
  [node_geometry.py:174](../../src/rigel/calibration/node_geometry.py#L174) gates the counts on the identical
  `js == strand_val` motif test that already zeroes the mass. Closing them needs per-transcript exon-span
  evidence no plumbed array carries.
* **The D-vs-B intron gap's mechanism** ("carried by seam/intron local belief") — code-reading inference, no
  experiment. *Experiment:* ablate the local-belief contribution on non-(T,T) edges, replay phase-D, 7 conditions.
* **A suspected nascent mirror at [node_geometry.py:123](../../src/rigel/calibration/node_geometry.py#L123)**
  (`eff_rna` too high at TSS/TES) — untested hypothesis, downgraded by its own author. Do not assert the sign.
  *Experiment (cheap):* count faces where `free_s` is True at a TSS/TES flank.
* **Real data.** The 0-of-1,146 invariant is a **theorem over one GTF**, not seven measurements — the
  partition is identical across conditions *by construction*, and both truth and predicate are read off the
  same annotation. **Untested against intron retention, unannotated exons, and annotation error** — exactly
  where the predicate is wrong and the gate deletes a **true** message, silently, without graceful degradation.
  **No synthetic suite retires this.** LBX0190 cfRNA (no nascent ⇒ the premise holds and silence is a genuine
  win) is the maximally safe first real substrate; MO_3005 the honest one.

---

## 6.6 Session-close dispositions (2026-07-16)

The mature-crossing gate initiative is **COMPLETE and landed** (Phases 0–4 + D4; uncommitted on
`calib-ambig-init-wip`). Every remaining item has an explicit home:

| item | disposition |
|---|---|
| **The gate** (Phases 0–4) | ✅ DONE — landed, tested, measured. |
| **D4 strand-aware routing** | ✅ DONE — `node_geometry._spliced_faces`, test-first. |
| **Phase 2 `pr>0.0` guard** | ✅ DONE — bit-identical. |
| **Issue #7 docstring** (§6.2) | ✅ DONE — `effective_length.py:86` corrected. |
| **Issue #5 `eff_rna ≥ 1`** (§6.1) | ⏸ DEFERRED to the **calibration-roadmap session**. Basis approved (D2: the ≥1-start-position unit, no magic number); ready to implement (gate the RNA imputation on the source's `eff_rna ≥ 1`). Belongs with the eff-length / prior work. |
| **Issue #7 geometry fix** (§6.2) | ⏸ DEFERRED to **priority #3** (only matters once the spliced channel is un-gagged). |
| **Count term `1/M_observed`** (§6.4) | ⏸ DEFERRED to **`dof_pie_model_fix.md`** — it is precisely the precision surgery that document rethinks (`n_src = fbp·sm` is the pie-inflation `pr` inflates); do not touch it before the pie model. |
| **Item 2 / the pie model** (§6.3) | ➡ **NEXT SESSION, TOP PRIORITY** — `dof_pie_model_fix.md`; re-derive its justification on the post-gate residual. |
| **Nascent-sourcing regression** | ➡ **[`nascent_rna_sourcing_regression.md`](nascent_rna_sourcing_regression.md)** — a dedicated future session (density-discrepancy nascent factory). Accepted for now: gDNA siphons some nascent on unstranded AMBIG introns. |
| **Phase 5 (4-arm A/B)** | ⛔ SKIPPED — documentation-only; D1 chose the arm on the physics, not the metric. The digit-exact D reproduction (§4 Phase 4) is the record. |
| **Phase 7 (golden regen)** | ⏸ DEFERRED to a stable checkpoint — the branch already owed a regen (21 pre-existing golden failures), and more calibration work is coming; regen ONCE when the roadmap work settles, not now. |
| **Production accuracy verdict** | ⏸ NOT YET — mid-implementation of the reference/prior; Phase 6 ran only as a bug-finding lens. |

---

## 7. Established vs assumed

| | |
|---|---|
| exon↔intron seams carry 0 crossing mature; the message fires on 100% | **measured** (oracle, 1,146 seams, 7/7 conditions) |
| item 1 §3.2's content change is an exact algebraic no-op | **measured** (`max|Δmode| = 4.26e-14`; 0/579 floored modes differ; 4 replays) |
| §3.2 delivers 0.1517 vs its predicted 0.0098 | **measured** (3 independent end-to-end implementations) |
| §3.2's nascent estimator has SNR 0.13 at the oracle's nascent level | **measured** (oracle-perfect source belief) |
| `mrna_active_strands` == the stated same-strand-exon rule | **proven exhaustively** (0/256 pairs disagree) |
| D: intron improves 7/7 (0.1573→0.0456), ALL 6/7 | **measured** (reproduced digit-for-digit by 3 parties) |
| 0 of 5,543 gated edges carry a live `n_mat` ⇒ `pr = 0` exactly | **measured** (exact static enumeration, not sampled) |
| `_comb` cannot be contaminated by a `pr=0` message | **proven structurally** (`mo ≤ 20.7`) + measured (0 non-finite) |
| B destroys 533,708 units/condition of true nascent relay; D preserves it | **measured** (nrna_present; exactly 0 on nrna_none) |
| E blocks intron→exon nascent | **derived** from the predicate; consistent with the traced chain |
| D is NOT the metric optimum (E beats it on ALL 7/7) | **measured** (head-to-head, one snapshot) |
| removing the `_EPS` floor is a bit-identical no-op | **measured** (0.000e+00, 7/7, two tree states) |
| the 191 short-region RNA nodes have true RNA mass = 0 exactly | **measured** (oracle) |
| `spliced_side_eff_length` is 2–199× low, and inert (≤0.0001) | **measured** (2 from-scratch brute forces + A/B) |
| exactly 1 new test failure; 21 golden failures pre-existing | **measured** (baseline 21/1075 → 22/1074) |
| the 1,317/1,317 TRUTH ⊆ `mrna_active` invariant | **NOT ESTABLISHED** — contested, unresolvable as keyed (§4.3) |
| the gate's behaviour on real annotation error | **assumed** — no synthetic suite can test it (§6.5) |
