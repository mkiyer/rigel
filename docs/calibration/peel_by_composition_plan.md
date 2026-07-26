# The composition peel — implementation plan

**Status: PLAN, not yet implemented. Owner-directed 2026-07-25.** Read `ROADMAP.md`, then
`SESSION_2026_07_25_HANDOFF_10.md` (§9–§12 are the evidence this plan rests on), then this.

The peel's **subtraction** is the last absolute-density operation left from the pre-composition solver. It is
the only place in the message path where a scale error does not cancel, and it is what contaminates introns.
This plan replaces it with a **share** that carries its own uncertainty, so that "intronic unspliced fragments
are gDNA until proven otherwise" is **reproduced by the uncertainty** rather than asserted by a rule.

---

## 0. Objective, and what "done" means

**Objective.** The RNA-continue message across a seam must be a *scaling* of the exon's RNA by a share `w`
estimated at the boundary, with `Var(log w)` propagated into the message precision — so that a boundary whose
RNA is consistent with zero emits a claim that is both near-zero **and weak**, and a boundary with real RNA
emits a proportionate one.

**Pre-registered acceptance gates.** All five must hold to land:

| # | gate | current HEAD |
|---|---|---|
| G1 | per-condition A/B at **refit=0**: no worse than HEAD on aggregate, and ≥ as many better-than-`m8` as HEAD | 0.0884, 16 better / 11 worse / 5 flat vs `m8` |
| G2 | per-condition A/B at **refit=1**: same | 0.0678, 23 better / 4 worse / 5 flat vs `m8` |
| G3 | delivered ρ_R at `exon\|intron` boundaries in `nrna_none` (oracle = **exactly 0**) ≤ HEAD | 0.2030 |
| G4 | boundary `errQ1conf` (share of a class's error on the most-confident quartile) must not regress | single 6.6 %, AMBIG 2.2 % |
| G5 | the `√2·σ_own` DL pin-safety inequality preserved; `message_variance_mc.py` 0 failures | holds |

**Rollback rule, pre-agreed.** If G1 or G2 fails after the ablations in §6 have been run, revert to HEAD and
record the measurement. Three wirings have already lost (HANDOFF_10 §10.2); the prior on success is not high,
so **every step below is specified to produce a measurement that stands on its own even if the step is
reverted.** That is the point of the sequencing in §7.

---

## 1. The quantities, and exactly where each comes from

Per **(node i, face d ∈ {0,1}, strand s ∈ {+,−})**. All already exist; nothing new is measured.

| symbol | meaning | source |
|---|---|---|
| `ρ_μ` | spliced (leaving) density at this face | `spl_p_f[d][i]` / `spl_n_f[d][i]` = `SP[d]/E_spl[d]` |
| `n_spl` | spliced COUNT at this face | `geometry.spliced_n_pos_{left,right}` / `..._neg_...` (already in `_uni_static`) |
| `v_μ` | `Var(log ρ_μ)` | `1/n_spl` — measured, composition-certain (M1) |
| `ρ_ν` | unspliced RNA density **at the boundary** | §2 precedence |
| `v_ν` | `Var(log ρ_ν)` | §2 precedence |
| `A` | the exon's RNA density reframed to the boundary | `(rp[s] or rn[s]) · r` — unchanged |
| `v_A` | `Var(log A)` | the relay's existing damped precision ⊕ `σ²_transfer` — unchanged |

⚠ **`v_μ` must use the spliced COUNT, never the spliced MASS.** The accumulator deposits fragments
*fractionally* across nodes, so at a junction face the median count is **33** against a median mass of **11**
(audited). Using the mass would over-state `v_μ` ~3× and make every peel message spuriously weak. The count
fields (`geometry.spliced_n_{pos,neg}_{left,right}`) already exist and are already published in
`_uni_static`; the mass fields (`spliced_{pos,neg}_{left,right}`) are what forms the DENSITY `ρ_μ`. Both are
needed, for different purposes — do not conflate them.

## 2. The LEVEL — a constant-free precedence (the "first step" of the two-step)

`ρ_ν` and `v_ν` are taken from the first source that has composition evidence. Every test is a **presence
test**, no thresholds:

1. **The boundary's own evidence** — `τ_own(i) > 0` (strand tilt, or a density deconvolution if one applies
   here). Then `ρ_ν = op[i]` / `on[i]` and `v_ν = transport_seed_logvar(v_log_fr(i), n_unspl(i))`.
2. **Else the node ACROSS the seam** — the neighbour opposite the exon, if it has `τ_own > 0` (in practice the
   factory-solved intron). Then `ρ_ν = op[far]·r_{far→i}` and
   `v_ν = transport_seed_logvar(v_log_fr(far), n(far)) + σ²_transfer(far→i)`.
3. **Else no evidence** — `v_ν = +∞`. See §4.3: the message carries **zero precision**, which is the honest
   "I don't know" and reproduces the safe default without asserting it.

This is exactly the owner's *"the intron message can set the gDNA vs unspliced RNA level. Or the strand. Then
the composition scale works."* Note (1) is checked before (2) so a boundary that knows its own composition is
never overridden by a neighbour's — the same "a direct measurement outranks an imputation" ordering that the
intron study established.

⚠ **Measured expectation, so we are not surprised.** Source (2) is the estimator whose *mode* was scored in
HANDOFF_10 §10.1: |Δw| = 0.192 / **0.294** / 0.126 / 0.052 / 0.100 against a self-solve's
0.110 / **0.628** / 0.156 / 0.070 / 0.160. It is better but still over-claims on the unstranded arm. **The
variance is what must carry that** — see §3.

## 3. The LAW — mode and precision together (M10, already derived and MC-validated)

```
    w        = ρ_ν / (ρ_ν + ρ_μ)                                  enrichment-free (e cancels identically)
    w_μ      = 1 − w                                              the spliced share
    t_p      = A · w                                              a SCALING, not a difference
    Var(log t_p) = v_A + w_μ²·(v_ν + v_μ)                          convex weights ≤ 1
```

Both pure laws are **landed and unit-tested**: `enrichment_frame.peel_continue_share` and
`peel_share_logvar`. MC: `message_variance_mc.py` M10 rel 0.20–0.96 %, M10b exact to 1e-12.

**Why the variance term is the whole point.** As `ρ_ν → 0`: `w → 0` *and* `w_μ → 1`, so
`Var(log w) → v_ν + v_μ`, which is large exactly when `ρ_ν` is poorly known. A boundary reporting
`ρ_ν = 0.0006 ± 0.0012` (HANDOFF_10 §12.3 — measured, and *conservative*: z2 = 0.13–0.55) therefore emits a
near-zero claim at near-zero precision. **That is "gDNA until proven otherwise" derived, not asserted** — and
it is precisely what the three failed wirings threw away by using a ratio of modes.

## 4. Degenerate cases and structural interactions — enumerate ALL of them

Each is a required unit test (§6.1).

**4.1 `ρ_μ = 0` (no spliced flux at this face).** `w = 1`, `w_μ = 0` ⇒ `Var(log w) = 0`. Nothing splices away,
nothing is peeled, no variance added. The message passes through untouched. Already the documented limit of
`peel_continue_share`.

**4.2 `ρ_ν = 0` with `ρ_μ > 0`.** `w = 0` ⇒ `t_p = 0`. Three downstream operators read the *density*, not the
precision, and each must be checked:
* `_pin_v`: `sp = where(pp_ > 0, p, op)` — reads the precision first, so a zero-precision component correctly
  falls back to the node's own density. **Safe as written.**
* **the λ-emission gate** `ttau = where((tg > _EPS) & ((tp + tn) > _EPS), ttau, 0)` — reads the DENSITY. A
  legitimately-zero RNA density would silence the composition message. Its stated intent is "a source that
  carries only ONE component has no composition claim", i.e. a test of what was *supplied*; a component at
  zero precision is not supplied, and a component at zero density but live precision **is**. ⚠ **This gate
  must move to a precision test** — but as its own commit and its own A/B (§7 step 5), never as a rider.
* `mismatch_gap(t, o)` — `contradicted = live_t ^ live_o`, so `t_p = 0` against a live own belief marks the
  message contradicted and DL kills it. That is correct *if* the zero is a claim; it is wrong if the zero is
  "no opinion". Resolved by the same precision test.

**4.3 `v_ν = +∞` (no evidence anywhere — §2 case 3).** Then `Var(log t_p) = ∞` ⇒ precision 0. The mode is
undefined; set `w = 0` so `t_p = 0`, consistent with §4.2, and rely on the precision to make it inert. **Must
not produce a nan**: compute with the `∞`-masked-before-the-product idiom already used in
`node_init.own_precision`.

**4.4 AMBIG nodes.** `w` is per-strand; a strand that is structurally dead already has `ρ_ν = 0` and is zeroed
by the existing `fp_a`/`fn_a` gates *after* the peel. Order is unchanged.

**4.5 Mature contiguous (`exon|exon`, retained intron).** Little or nothing splices at that face ⇒ `ρ_μ` small
⇒ `w → 1` ⇒ near-no peel. The correct behaviour falls out of the measurement; **no structural gate is
needed**, which is what retires HEAD's `mrna_active` rule (§5).

**4.6 Both faces carry spliced.** Measured at 4.6–5.2 % of junction boundaries (HANDOFF_10 §9). `w` is
per-face by construction; each direction uses its own face. No special case.

**4.7 `n_spl = 0` while `ρ_μ > 0`.** Cannot occur (the density is built from the count) but guard: `v_μ = ∞`
⇒ zero precision rather than a division by zero.

**4.8 The relay's running belief.** `_relay` computes `w` from the node's *init* (`op`/`on`), not its fused
value, because when node `i` is processed its own entry has not been updated yet. This is deliberate and must
be stated: the share is a **pre-pass** quantity. Do **not** attempt to use the running belief — it would make
the share depend on the traversal direction and break the forward/backward symmetry.

## 5. What is REMOVED

* **The subtraction** `tp = max(tp − spl_*_f[df], 0)` — both sites (`_relay`, `_transport`). This also retires
  the **zero-truncation defect** (HANDOFF_7 §8.6: a fully-consumed peel emitted `ρ_ν = 0`, "no RNA continues
  past here", at a *live* precision) for free, since `w ∈ [0,1]`.
* **HEAD's `mrna_active` gate** — its justification is void (the owner's correction: mature and nascent are
  the same species; only spliced vs unspliced is observable), and §4.5 shows the share subsumes its effect.
  Remove it **only in step 4** (§7), so its contribution is measured separately.
* **NOT removed:** `peel_rna_logvar` (M3). It is unused by the solver but is the record of what the difference
  costs, and M10's convexity unit test asserts against it directly.

## 6. Validation ladder

### 6.1 Closed-form unit tests (new, in `tests/calibration/test_enrichment_frame.py`)
Four already exist for the laws. Add, for the **wiring**:
1. `ρ_μ = 0` ⇒ message identical to no-peel, and `Var(log w) = 0` (§4.1).
2. `ρ_ν = 0` ⇒ `t_p = 0` and precision 0 — and **no nan** (§4.2/4.3).
3. `v_ν = ∞` ⇒ precision 0, mode finite (§4.3).
4. enrichment invariance end-to-end: scaling `(ρ_ν, ρ_μ)` by a common `e` leaves `t_p` and `Var(log t_p)`
   unchanged to 1e-14 — the property the whole change exists for.
5. the M3→M10 contrast at a fixed operating point: `Var(log t_p)` strictly below the subtraction's, for
   `u ∈ {2, 3, 10}`.

### 6.2 MC — already done, re-run as a gate
`message_variance_mc.py` M10 (rel 0.20–0.96 %) and M10b (exact to 1e-12). **0 failures over M1–M10** is a
landing gate. No new MC is required; the law is unchanged by the wiring.

### 6.3 Replay fidelity
`scratchpad/ss_cap_3_replay.py` / `ambig_2_mechanism.py` must still reproduce the shipped ψ **bit-exactly**
(`max|Δ| = 0`) after the change. If they do not, a channel is being published inconsistently with what ψ is
handed — fix that before reading any A/B.

### 6.4 Per-condition A/B
`P0_REFIT=0|1 pass0_oracle_bench.py --arm …`, all 32 conditions, `OMP_NUM_THREADS=1`. Report better/worse/flat
against **both** `m8` (the pre-boundary-work baseline) and HEAD. Call out explicitly: `stranded ss_0.99`,
`unstranded × capON`, `verystrong`, `gdna_none`.

### 6.5 Certification matrix — by boundary class
`scratchpad/b1_boundary.py`, per flanking pair × junction present/absent. Every cell must be reported, not
just the aggregate:

| class | why it is a distinct case |
|---|---|
| `exon\|intron`, junction | the main peel case; 39–61 % of boundary error |
| `exon\|intron`, no junction | the 14–17 % whose junction was never observed |
| `exon\|exon`, junction | mature contiguous **and** splicing — the owner's TA/TB geometry; `w → 1` expected |
| `exon\|exon`, no junction | contiguous, no peel; must be untouched |
| retained intron | `coarse_type` collapses it (HANDOFF_7 §8.3) — verify it is not silently peeled |
| both enrichment orderings | boundary < exon (72–97 %) and boundary > exon (3–28 %) |

### 6.6 The regime split that decides success
Report `w` against the truth separately for **saturated** and **unsaturated** boundaries, using the
constant-free detector of HANDOFF_10 §11.4 (is the factory log-factor still rising at the top of the λ grid).
Expected, and the plan fails if not observed:
* saturated / `nrna_none`: `w → 0` with large `Var(log w)` ⇒ delivered RNA → 0 (gate G3);
* unsaturated / `nrna_present`: `w` tracks the true 0.17–0.41 rather than HEAD's 0.

⚠ **Re-measure `w_true` on junction-bearing seams only.** HANDOFF_10 §11.4 records that the earlier figure was
contaminated by boundaries with no spliced mass, where the share is undefined and the probe defaulted it to 1.

### 6.7 Standing gates
`pytest tests/ -q` (1236 pass / 2 xfail / 2 xpass), `ruff check src/ tests/ scripts/`, goldens regenerated
**last**.

## 7. Sequencing — each step independently measurable

Each step is its own commit with its own A/B, so a failure localizes and every measurement survives a revert.

| step | change | what it isolates | expectation |
|---|---|---|---|
| 1 | plumb `(ρ_ν, v_ν)` per §2 and compute `w`, `Var(log w)` — **computed, not used** | that the level precedence and the shares are well-formed; no behaviour change | bit-identical A/B |
| 2 | replace the subtraction with `t_p = A·w`, keep HEAD's gate, **precision unchanged** | the MODE effect of the share alone | ~neutral; this is close to the arm that scored 0.0893/0.0690 |
| 3 | add `w_μ²(v_ν + v_μ)` to the message log-variance | **the variance effect — the actual thesis** | the step that must pay; G3 should move |
| 4 | remove HEAD's `mrna_active` gate | whether the share subsumes the structural rule | should be ~neutral if §4.5 holds |
| 5 | λ-emission gate: density test → precision test (§4.2) | a pre-existing defect the plan surfaced | its own A/B, may be landed or reverted independently |

**Do not collapse steps 2 and 3.** The whole claim of this plan is that the mode alone is not enough and the
variance is what makes it work; measuring them together would leave that untested — and step 2 in isolation is
expected to look like the failed arm, which is the point.

## 8. Risk register

| risk | how it shows | response |
|---|---|---|
| The §2(2) intron level still over-claims on the unstranded arm | G3 fails; `unstranded × capON` regresses | this is the known |Δw| = 0.294; if the variance does not absorb it, the plan's thesis is **refuted** — record and stop |
| `w = 0` trips the λ gate / DL contradiction mask | composition messages vanish at seams; `exon\|exon` `c_tau` → 0 (the exact signature seen when the junction gate over-reached) | step 5; check `c_tau` per boundary class in §6.5 |
| The share becomes over-confident where `v_ν` is optimistic | `errQ1conf` at boundaries regresses (G4) | G4 is a hard gate, not advisory |
| Loss of the M3 `u`-weighting removes a damping that was compensating elsewhere | stranded arm regresses | ablate by restoring `peel_rna_logvar` on the same operating point |
| Traversal-order dependence | forward and backward passes disagree | §4.8 — the share is a pre-pass quantity; assert `_relay` and `_transport` compute the identical `w` |

## 9. Explicitly out of scope

* The **exon-face reframe error** (0.4–1.3 nats). It is irreducible (HANDOFF_10 §9.3 — unchanged even with the
  ORACLE `f_g`) and is precisely what the share makes harmless rather than fixes.
* The **ψ vertex** / `_JEFFREYS_REF`. HANDOFF_10 §11–§12 close this: the reference is required for properness,
  the median is already the best readout, the grid is not the problem, and the precision is already honest.
* The **hyperprior**. Nothing here needs it, and it is not the current task.
