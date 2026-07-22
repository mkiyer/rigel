# Message-Mode Implementation Plan — cautious, stepwise, verify-in-real-code

**Status (2026-07-21):** the phased plan to land the **validated mode arithmetic**
([`message_propagation_arithmetic.md`](message_propagation_arithmetic.md)) into `bp_solver._scan`, one
independently-verifiable step at a time. **Goal:** finish the pass-0 solve with the *correct message-mode
arithmetic*. **Precision is the NEXT thread** (owner directive) — this plan touches mode VALUES only, never `pr_c`.

Companions: the derivation ([`message_propagation_arithmetic.md`](message_propagation_arithmetic.md)), the
reframe ([`message_system_derivation.md`](message_system_derivation.md) §6A), the refactor survey
([`message_system_refactor_plan.md`](message_system_refactor_plan.md) — this supersedes its mode-extraction step
A2 with the gate-equality change).

---

## 0. The change in one sentence, and the honest impact expectation

**The whole arithmetic fix is a single predicate flip** at [bp_solver.py:593](../../src/rigel/calibration/bp_solver.py):

```python
   use_shift = not is_exon_node[lsrc] and not is_exon_node[i]     # NOW: exon-membership proxy
   # →
   use_shift = (fp[lsrc] == fp[i]) and (fn[lsrc] == fn[i])        # TARGET: gate-equality (§7 of the derivation)
```

Everything else the derivation validated is **already in the code** and stays untouched:
- the per-component **densities** `rho_g / rho_pos / rho_neg` (the frame-invariant currency), including the
  **mature reconciliation** `+ SPs/esp − absorb_p` ([bp_solver.py:606-610](../../src/rigel/calibration/bp_solver.py));
- the **shift** itself (`Mg = rho_g·egd`, normalize by the imputed total `_den`, `mo = log(Mg/_den)`,
  [bp_solver.py:631-641](../../src/rigel/calibration/bp_solver.py)) — bit-for-bit the derivation's §3/§4a formula.

**The graph is bipartite B-R-B-…-R-B** ([`node_chain.py`](../../src/rigel/calibration/node_chain.py)) and
`is_exon_node` marks **REGION nodes only** ([bp_solver.py:415](../../src/rigel/calibration/bp_solver.py)). So a
region never touches another region — **there is no `exon↔exon` edge**; an exon region's neighbours are always
boundaries. Every edge is `region ↔ boundary`.

**The flip has exactly TWO effects** (both correct per the derivation; both must be measured):

| # | edge class | gate | today | after `use_shift = gates_equal` | note |
|---|---|---|---|---|---|
| — | intron/intergenic ↔ boundary, **gates equal** | equal | shift | **shift** (unchanged) | already correct (§4a) — the vast majority |
| **B** | **boundary ↔ exon** (inner junction) | equal | density mode | **shift** ← the win | §4b: cliff-free under binary capture, but a **real ~2× cliff under proportional capture** where density ÷`md` *under-calls* f_g and the shift cancels it |
| **A** | intron/intergenic ↔ boundary, **gates differ** (strand transition) | **unequal** | shift | **density mode** ← a change | the shift *over-claims* across an unequal gate (§5); density is the safe placeholder. **Stage 3 measured this class EMPTY** (strand transitions live at exons, not non-exon edges) — so in practice it is a no-op |
| — | intergenic-seam ↔ exon (TSS/TES) | **unequal** | density mode | **stays density mode** | already out of the shift; the gDNA-only **lower anchor** is a ψ/one-sided change ⇒ **deferred to the precision thread** |

Effect **B** (the exon under-call receiving side) is the win; effect **A** (removing an *incorrect* shift on
gate-unequal non-exon edges) is a smaller, correctness-improving side-effect that Stage 3 will size and Stage 6
measures as its own delta. The rollout is therefore **additive** — enable the shift for progressively larger
safe subsets, ending at the pure `gates_equal` form — so each effect is isolated and measured, never flipped
blind.

**Honest impact expectation.** Under *binary* capture the boundary↔exon edge is cliff-free ⇒ density ≈ shift ⇒
the flip is nearly inert. Under *proportional* capture (the realistic model — the boundary crossing is ~0.5×
enriched vs the exon interior, `message_propagation_arithmetic.md` §6a) the density mode **under-calls the exon
f_g by the capture-cliff factor** and the shift cancels it ⇒ the flip may correct a real slice of the pass-0
under-call. **We do not know the magnitude until we measure it.** The *big* remaining levers — the seam lower
anchor and the exon-as-source robustness — are **precision-thread** work; this plan lays the correct, verified
foundation they build on.

> **Stage 1 measured a second-order gate on the impact: emission.** A mode only matters where a message is
> *emitted*. On unstranded+capture the `emit_g` gate suppresses almost all `boundary→exon` gDNA messages (τ=0),
> so the mode flip has near-zero surface *there* — the correct division of labour is: **this thread makes every
> emitted message correctly framed; the precision thread opens the emission gate (τ) so the right messages flow.**
> Neither alone fixes unstranded+capture; together they do. The mode flip's *measurable* surface is on stranded
> and no-capture-unstranded edges (where messages already emit).

---

## 1. The mode/precision entanglement (why we go slow)

The density mode was never only a mode — it was also a **robustness crutch**: on an unstranded exon **source**
whose composition is degenerate (~0.5), it anchors to the dst's *observed* total instead of trusting the source,
which is a *precision* behaviour ([`message_system_derivation.md`](message_system_derivation.md) §6A). So the
mode flip has two faces:

- **exon-as-DESTINATION** (`boundary → exon`): the exon *receives* a better-framed message. **Safe**, and it is
  exactly the direction that feeds the exon under-call. ← land first.
- **exon-as-SOURCE** (`exon → boundary`): the exon *propagates* its composition; if it is degenerate and the
  precision is not yet honest, this can revive the historical shift-on-exon regression
  ([`cliff_message_derivation.md`](cliff_message_derivation.md) §9). **Precision-limited** ← isolate, measure,
  likely defer the default to the precision thread.

Every behavioural stage is therefore **flag-gated for A/B** and split by direction. The flag is a **temporary
env var** (the `RIGEL_TAU_GAG` pattern — used for the A/B then removed), **not** a `CalibrationConfig` field: the
owner's directive is no lingering knobs (`productionize_cleanup`, the δ-pin removal A0).

---

## 2. Verification assets (all stages reuse these)

- **Node-level `_capture` hook** — `node_sweep(..., _capture=dict)` ([bp_solver.py:270](../../src/rigel/calibration/bp_solver.py))
  already dumps per-edge diagnostics; `scripts/debug/ablate_replay.py` already replays the real prior-free
  `node_sweep` once per condition capturing its inputs. Extend to emit, per edge: `(rho_g, rho_pos, rho_neg,
  E_c^dst, branch, mode_c, f_g_after)`.
- **The MC harness** `scripts/debug/chain_mode_mc.py` — the derived-formula oracle (`_fg_from_shift`,
  `_fg_from_densities`) each unit test checks against.
- **The oracle** `scripts/debug/oracle.py` + the toy driver `scripts/debug/toy_prod.py` — per-node true f_g on a
  controlled toy (owner directive: develop on toys).
- **Benchmarks** — the cached `ambig_dense_10mb` (32 scenarios) mass-weighted f_g error + the **`gdna_none`
  false-positive guard** (the phantom must not return). `scripts/debug/pass0_bench.py` A/B via the env flag.
- **Golden byte-identity** — the gate on every refactor stage (`pytest tests/` + `tests/golden/`).

---

## 3. The stages (each independently verifiable; STOP-and-review after each)

### Stage 1 — Node-level verification harness ✅ DONE *(byte-identical; `bp_solver` `_capture` hook + `scripts/debug/mode_verify.py`)*
Added a read-only per-edge record to `_scan`'s inert `_capture` hook (`_edge_modes`: src/dst densities, dst
eff-lengths, branch, modes) and `mode_verify.py` runs the real `calibrate()` (via `toy_prod.run(_debug=…)`),
recomputing the derived shift/density mode per edge and comparing to the code + the per-region oracle f_g.
**Gate met:** 21 golden + 251 calibration tests byte-identical; lint clean. **Results:**
- **CODE == DERIVED, max|Δmode_g| = 0.00e+00** across all scenarios — the code computes `message_propagation_arithmetic.md`
  §3/§4 bit-for-bit (§8's "confirmed correct" is now verified in real code, not asserted).
- **Density ≈ shift on B-safe (`boundary→exon`) without capture** (0.112 vs 0.118; shift marginally closer to the
  oracle 0.193) — consistent with cliff-free binary capture; the ~2× proportional-capture gap needs the real
  capture model to show.
- **On seam (gate-unequal) edges the SHIFT over-claims to f_g=1.0** ("pure gDNA") vs density 0.787 (oracle 0.193)
  — a live confirmation of §5 and of the plan's decision to keep seams OFF the shift (they need the lower anchor).
- **⚠ Emission-gating interaction (expectation-setter).** On the hardest scenario (**unstranded + capture**) only
  2 of ~20 edges **emit a gDNA message at all** (both seams); the `boundary→exon` edges emit nothing (`emit_g =
  sm>0 AND lam_ev`; unstranded ⇒ τ=0 ⇒ `lam_ev=False`). So **the Stage-4/5 mode flip has almost no surface on
  unstranded+capture** — that under-call is an *emission/precision* problem (the τ-gate), not a mode problem. The
  mode flip acts where messages ARE emitted (no-capture unstranded had 16 emitting edges incl. B-safe). This is
  recorded in §0's impact expectation.

### Stage 2 — Extract the mode helpers into pure, tested functions ✅ DONE *(byte-identical; golden)*
Extracted the two modes as module-level pure helpers — `_mode_shift(mass_c, den, comp_fl)` and
`_mode_density(rho_c, eff_c, md)` — and replaced all four inline call sites (gDNA / +RNA / −RNA / RNA-total) with
`_mode_shift(…) if use_shift else _mode_density(…)`. Whole expressions moved verbatim (a function wrap does not
reorder floats). **Gate met:** 279 calibration+golden+sim tests byte-identical, lint clean, sweep time in noise,
and `mode_verify` still reports `CODE==DERIVED = 0.00e+00` (the extraction changed no arithmetic). Two named
modes now match the derivation's §4a (shift) vs §8 (density) one-for-one — the Stage-4 change is a predicate swap
*outside* this tested unit. **Unit tests** (`test_bp_solver.py`): `_mode_shift` == the MC composition fraction
`Mg/(Mg+Mr)`; `_mode_density` == `log(ρE/md)`; and (review #2) both **finite under ρ_c ≤ 0** via the derived
one-fragment floor — **no added epsilon**.

> **Domain safety is already present — preserve it verbatim, add nothing.** The mature subtraction
> `rho_pos = … − absorb_p` can go ≤ 0, but every `log` is already guarded by a **derived one-fragment floor**:
> the shift path floors `Mp = max(rho_pos, 0.0)·erd` then `log(max(Mp/_den, comp_fl))` with `comp_fl = 1/md`; the
> density path is `log(max(rho_pos, 1/erd) / (md/erd))`. These are the resolution floor (a node can never be more
> certain than its one-fragment opportunity), **not** an arbitrary epsilon. Extraction moves them unchanged. A
> unit test pins `_message_mode` finite (no NaN/−inf) for `rho_pos < 0`. **Do not introduce a new `eps=1e-12`** —
> it would break byte-identity and add a magic number the owner's directive forbids; the derived floors already
> cover the domain.

### Stage 3 — Diagnose the full flip set ✅ DONE *(zero source change — post-hoc from the Stage-1 capture)*
The Stage-1 `_capture` record already carries `fp/fn/exon/use_shift` per edge, so `gates_equal` and the
classification are computed **post-hoc** in `mode_verify.flip_set_report` — no source edit at all (trivially
byte-identical), and `gates_equal` gets introduced in the code only where it is *used* (Stage 4), no dead
variable. Ran 5 scenarios (1- and 2-gene, ± strand, ± capture). **Findings:**
- **Self-consistency OK on every scenario** — the 3 changing classes {B-safe, B-src, A} are *exactly*
  `{use_shift != gates_equal}` (per-edge cross-check, 0 mismatches). The flip is fully characterized.
- **Blast radius (aggregate):** B-safe **14 edges (12 emit)**, B-src **14 (6 emit)**, **A = 0 (0)**, clean 76 /
  seam 28 (both unchanged). So the flip has real surface, concentrated in B-safe (Stage 4), the safe direction.
- **⚠ Class A is empirically EMPTY** — 0 of 76 non-exon edges change, even with two opposite-strand genes. This
  is *structurally expected*: a strand transition happens at a gene's first/last **exon** (a seam edge), never on
  a non-exon↔non-exon edge (an intron is always exon-flanked; a gene body has one strand). So the reviewer's
  effect-A concern is **negligible in practice** and Stage 6's "class-A delta" is a no-op — the additive B-safe +
  B-src rollout already *reaches* the pure `gates_equal` form. *(Verify on a real/overlapping-gene scenario
  before relying on it — antisense/overlapping loci could in principle produce a class-A edge.)*
- **Refines the Stage-1 note:** on unstranded+capture the *gDNA* message is emission-gated, but B-safe edges
  still emit **RNA** messages (the RNA-total λ-factor influences f_g), so the flip is **not fully inert** there —
  Stage 4's A/B measures the actual f_g effect.

The flip is rolled out **additively** (§0): each behavioural stage turns the shift ON for one more edge class via
a temporary `RIGEL_GATE_SHIFT` env flag (levels: `off` → `dst` → `src` → `all`), so a stage's delta is exactly
one class. The final `all` level equals the pure `use_shift = gates_equal`.

### Stage 4 — Enable the shift for **class B-safe** (`boundary → exon`) — ⛔ MEASURED, REGRESSES, NOT LANDED
`use_shift = old OR (gates_equal AND src is boundary AND dst is exon-region)`, behind `RIGEL_GATE_SHIFT=dst`
([bp_solver.py](../../src/rigel/calibration/bp_solver.py)). Gates 1–2 passed but **Gate 3 FAILED** — the flag
stays **OFF** (default), so nothing regressive shipped (byte-identical). **Result (`pass0_bench` on
`ambig_dense_10mb`, cached):**
- **Gate 1 ✅ byte-identical** flag-off (21 golden + 254 calibration).
- **Gate 2 ✅ (toy) toward-oracle** aggregate −0.0097, no toy scenario away — *but the toy only had over-call
  exons*, which hid the real effect.
- **Gate 3 ⛔ (full suite):** aggregate mwae **0.2033 → 0.2060 (+0.0027 > +0.002)**, and **7 scenarios regress
  > 0.02** (worst **+0.098**, `gdna300_ss0.50_capture_on`). The shift **uniformly lowers f_g at captured exons**:
  it *helps* low-gDNA (phantom −0.03) but *hurts* the high-gDNA **enriched-exon under-call** (+0.05…+0.10) — the
  exact cases we target. `gdna_none` did **not** regress (phantom guard held; it actually improved).
- **Diagnosis:** capture-driven (capture-off deltas ±0.003), consistent with the **§6a #3 mature differential
  cliff** under proportional capture (the shift cancels the uniform cliff, not the mature one). And on the
  enriched under-call cases the gDNA message is **emission-gated** (unstranded τ=0), so both modes under-call and
  the mode is only a second-order (here net-negative) effect. **⟹ the exon-edge shift is blocked on (a) the
  deferred mature capture-scale correction and (b) the emission/precision thread; it cannot land in isolation.**
  Stages 5–6 (which build on B-safe) are moot until then. The flag/code is retained as the A/B harness pending
  that work (revert if the mode flip is abandoned).

  *(original gate definitions, for when the blockers clear:)* **Gates (all four):**
1. **Byte-identical** with the flag `off`.
2. **Node-level** (Stage-1 harness): on the flipped edges, `code_mode == derived_shift` and f_g moves **toward**
   the oracle (not away). *This is the PRIMARY gate* (the benchmark is precision-limited).
3. **Benchmark A/B** (`ambig_dense_10mb` + `gdna_none` guard), with concrete acceptance:
   - aggregate mass-weighted f_g error is **non-inferior** (Δ ≤ +0.002, i.e. within MC noise) and ideally improves;
   - **no single scenario regresses by > 0.02** in mass-weighted f_g (a larger move is a review trigger, not an
     auto-block — check it is not a real bug);
   - **`gdna_none` false-gDNA does not increase** vs the flag-`off` baseline (it is ~10 % today, not zero — the
     guard is *non-regression*, not "exactly zero").
4. **Isolate the mechanism**: confirm movement is the capture-cliff cancellation (compare binary-≈-inert vs
   proportional-capture toy scenarios), not coincidence.

> **Short-exon caveat (Deferral #5 interaction).** The shift ratios source/dst eff-lengths, and `region_eff` (+1
> discrete) is not yet frame-consistent with `boundary/spliced` (continuous) — a 2–3 % artifact on **short exons
> (L ≲ fl_mean, ~<150 bp)** (`message_propagation_arithmetic.md` §6a #2). Expect small short-exon wiggles in Gate 3
> that are the *discretization* artifact, not the mode change; do not chase them, and re-check after Deferral #5.

### Stage 5 — Enable **class B-src** (`exon → boundary`) *(behavioural; flag=`src`; measured — the risky direction)*
Add the edges where the **source is an exon region** (dst is necessarily a boundary — bipartite). **Same four
gates**, plus the specific watch: a **degenerate unstranded exon source** must not confidently corrupt the
boundary (the §9 historical shift-on-exon regression). Expect this to be **precision-limited** — if it regresses,
that is the *diagnosis*, not a failure: it confirms the exon-as-source robustness belongs to the precision thread.
**Decision:** keep the Stage-4 win; leave Stage-5 behind the flag with a documented "revisit after precision."

*(Classes B-safe + B-src exhaust all exon-touching edges — an exon region only ever touches boundaries.)*

### Stage 6 — Switch to pure `gates_equal`, land the default, remove the flag *(behavioural; flag=`all`)*
The final level `all` sets `use_shift = gates_equal`. Its only delta beyond Stages 4–5 is **class A** — which
**Stage 3 measured empty** (a no-op on non-overlapping loci); this stage's job is therefore to (a) **confirm
class A is empty on a real/overlapping-gene scenario** (if any class-A edges appear, measure that delta on its
own — it should be neutral-to-positive: an over-claiming shift → the safe density placeholder), then (b) flip the
**default** to whatever Stages 4–5 justified (likely: B-safe ON; B-src pending precision), **delete the env
flag** (no lingering knob), regenerate goldens for the intended output change, and update
[`message_system_derivation.md`](message_system_derivation.md) §6A / [`message_system_refactor_plan.md`](message_system_refactor_plan.md)
to record what landed.

---

## 4. Explicitly DEFERRED to the precision thread (not this plan)

1. **The unequal-gate lower anchor** (intergenic-seam → exon, the TSS/TES crush). The **mode value** is already
   correct — it is the density gDNA arm `log(ρ_g·E_g/md)`. What is wrong is the **ψ application**: it enters as a
   two-sided point estimate (crush) instead of a **one-sided** lower bound `−½·pr·max(0, mode−λ)²`. One-sidedness
   + precision ⇒ precision thread ([`message_system_derivation.md`](message_system_derivation.md) §6, open #2).
   **This is the high-impact fix for the exon under-call** — sequenced deliberately after the arithmetic
   foundation is verified.
2. **Exon-as-source robustness** (Stage 5's precision-limited part) — the honest precision that makes a degenerate
   source self-silence (`τ→0 ⇒ pr→0`, §6A), retiring the density mode's crutch role.
3. **The spliced-measurement vs nascent-prediction precision merge** (the owner's open question) and the C4
   conditioning down-weight when mature ≫ nascent (`message_propagation_arithmetic.md` §4b).
4. **The mature capture-scale correction** under proportional capture (`message_propagation_arithmetic.md` §6a
   boundary #3) — the differential-cliff bias on the mature add/subtract; capture-model-specific.
5. **The eff-length discretization-consistency fix** (region `+1` vs continuous spliced/boundary,
   `message_propagation_arithmetic.md` §6a boundary #2) — an independent short-flank correctness item; touch
   `effective_length.py` carefully (it is load-bearing everywhere) and gate on golden.

---

## 5. Execution order & the one rule

**S1 → S2 → S3 → S4 → (review) → S5 → (review) → S6.** S1–S3 are byte-identical and zero-risk; do them back to
back. The behavioural stages walk the single flag `RIGEL_GATE_SHIFT` up its levels `off → dst → src → all`, one
edge class per stage, so every delta is isolated. S4 (`dst`, the under-call receiving side) is the
highest-value/lowest-risk change — **stop and review the A/B before S5.** S5 (`src`) is expected to expose the
precision boundary; S6 (`all`) adds class A and lands the default.

**The one rule:** every stage is gated on golden byte-identity (refactors) or the `gdna_none` guard +
node-level-toward-oracle (behavioural). No stage lands on the benchmark number alone — the **node-level
correctness vs the derivation** is the primary gate, because the benchmark is precision-limited until the next
thread.
