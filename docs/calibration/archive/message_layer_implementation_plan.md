# Pass-0 message layer — CONSOLIDATION implementation plan

**Status:** plan. **Date:** 2026-07-22. **Branch:** `calib-ambig-init-wip` (do **not** push to main).
**Derivations this implements:** `gdna_factory_and_seam_derivation.md` (vehicles, factory, §9 two-message
closure) and `message_arithmetic_reconciliation.md` Part E (node inventory, PEEL vs GRAFT).

---

## 0. The one mechanism

> Every node has an **active component set**. Every edge has a **per-component continuity gate**. A message
> carries **per-component densities, each with its own precision**, where *zero precision means "no
> information", never "this component is zero"*. A node **pools densities per component** across its flanks,
> then **normalizes once** against its own observed total. It solves **lazily** — only when the messages
> jointly cover its active set.

The variety of real topologies — exon↔exon seams, an intron of one transcript inside an exon of another,
opposite-strand overlaps, TSS/TES — is **data in the component set and the gates**, not branches in the solver.

### 0.1 Two facts that make this cheap

**(a) The mature-continuity gate already exists and is already correct.**
`node_geometry._boundary_strand_stats` computes

```python
mrp_l, mrn_l = mrna_active_strands(sig_l)
mrp_r, mrn_r = mrna_active_strands(sig_r)
mrna_pos = mrp_l & mrp_r      # +strand mature (contiguous exon) crosses the seam
mrna_neg = mrn_l & mrn_r
```

and it is already carried on `NodeStatics` as `mrna_active_pos` / `mrna_active_neg`. It gives the right answer
on **every** topology the owner named:

| seam | left / right | `mrna_active_s` | mature crosses? |
|---|---|---|---|
| splice junction | exon / intron | False | ❌ correct |
| **exon ↔ exon** (overlapping transcripts) | exon / exon | **True** | ✅ **correct — it does** |
| TSS / TES | intergenic / exon | False | ❌ correct |
| opposite-strand overlap | per strand | per strand | ✅ handled independently |

**`grep` shows it has no consumer in `src/`.** It is computed, documented, and unused. Step 1 is therefore
mostly *naming and wiring*, not new derivation — and it fixes the exon↔exon case for free.

**(b) The chain strictly interleaves B-R-B-R-B**, so **every edge has exactly one boundary endpoint**. The
edge's gates are the *boundary endpoint's own* arrays. No per-edge storage is needed.

### 0.2 The component model

| component | active on node | continuous across edge |
|---|---|---|
| gDNA | always | always (genomic) |
| nascent ± | `statics.free_pos/neg` | `free_s[a] & free_s[b]` |
| mature ± | `statics.mrna_active_pos/neg` | `mrna_active_s[the boundary endpoint]` |

The solver's simplex stays **3-component** `(f_pos, f_neg, f_g)` over the *unspliced* pool. Mature/nascent is a
**sub-split** that is only separable where a spliced channel exists (a splice-junction boundary). The
sub-components exist in the model to decide **vehicle validity**, not to add solver dimensions.

### 0.3 The vehicle rule, stated once

Transport a component only if the source **has** it, it is **continuous**, and the destination **has** it.

```
transported ⊆ dst.active
    transported == dst.active   ⇒  SHIFT   (normalize by the transported total; enrichment cancels)
    transported ⊊  dst.active   ⇒  DENSITY (normalize by the dst's OBSERVED total md; a partial claim)
src.active ⊋ transported        ⇒  PEEL the surplus out of the source composition first  (the c_b term)
```

This single rule reproduces every case we derived by hand:

| edge | dst active | transported | vehicle |
|---|---|---|---|
| intron → splice boundary | g, ν | g, ν | SHIFT |
| exon → splice boundary | g, ν | g, ν (**+ peel** the exon's μ) | SHIFT + `c_b` |
| **splice boundary → exon** | g, ν, **μ** | g, ν | **DENSITY** (μ can't cross) |
| **exon → exon seam boundary** | g, ν, μ | g, ν, μ | **SHIFT** (μ *does* cross) |
| intergenic → TSS seam | g | g | SHIFT (trivial, 1→1) |
| TSS seam → exon | g, ν, μ | g | **DENSITY**, gDNA only |

---

## 1. STEP 1 — the component/continuity model (no behaviour change)

**Goal:** land the model + a diagnostic that shows where the derived vehicle predicate disagrees with today's
`use_shift_g`. **Byte-identical output.** This produces the evidence Step 2 needs.

### 1.1 New: `src/rigel/calibration/component_model.py`

```python
"""Per-node ACTIVE component set + per-edge CONTINUITY gates — the single source of truth for
WHICH message vehicle is valid on an edge (docs/calibration/message_layer_implementation_plan.md §0).

The solver's simplex stays 3-component (f_pos, f_neg, f_g) over the UNSPLICED pool. Mature-vs-nascent is a
SUB-SPLIT tracked here only to decide vehicle validity: a composition SHIFT is valid iff the source can supply
every component the destination possesses, and mature is exactly the component that a splice junction cannot
carry. Nothing here is fitted; it is pure annotation structure.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .node_chain import REGION


@dataclass(frozen=True)
class ComponentModel:
    """Per-node component activity (bool[n_nodes]). gDNA is active on every node by construction."""

    active_nas_pos: np.ndarray
    active_nas_neg: np.ndarray
    active_mat_pos: np.ndarray
    active_mat_neg: np.ndarray
    is_boundary: np.ndarray

    @property
    def n_nodes(self) -> int:
        return int(self.active_nas_pos.shape[0])


def build_component_model(chain, statics) -> ComponentModel:
    """Name the arrays `NodeStatics` already carries. `mrna_active_*` is, for a boundary, the AND of its two
    flanks' exon bits — i.e. ALREADY the mature-crossing gate (`node_geometry._boundary_strand_stats`)."""
    return ComponentModel(
        active_nas_pos=np.asarray(statics.free_pos, dtype=bool),
        active_nas_neg=np.asarray(statics.free_neg, dtype=bool),
        active_mat_pos=np.asarray(statics.mrna_active_pos, dtype=bool),
        active_mat_neg=np.asarray(statics.mrna_active_neg, dtype=bool),
        is_boundary=(np.asarray(chain.kind) != REGION),
    )


def edge_vehicle(cm: ComponentModel, src: int, dst: int) -> tuple[bool, bool, bool]:
    """Return ``(shift_valid, peel_pos, peel_neg)`` for the edge ``src → dst``.

    * ``shift_valid`` — every component the DST possesses is transportable from the SRC. Otherwise the
      message is partial and must ride the DENSITY vehicle (normalized by the dst's own observed total).
    * ``peel_s``      — the SRC carries mature on strand ``s`` that the DST does not, so the source
      composition must be de-matured before transport (the ``c_b`` term).

    Exactly one endpoint of a chain edge is a boundary (B-R-B-R-B interleave), and that boundary's
    ``mrna_active_s`` IS the edge's mature-continuity gate.
    """
    b = src if cm.is_boundary[src] else dst
    mat_cross = (bool(cm.active_mat_pos[b]), bool(cm.active_mat_neg[b]))
    nas_cross = (
        bool(cm.active_nas_pos[src]) and bool(cm.active_nas_pos[dst]),
        bool(cm.active_nas_neg[src]) and bool(cm.active_nas_neg[dst]),
    )
    dst_nas = (bool(cm.active_nas_pos[dst]), bool(cm.active_nas_neg[dst]))
    dst_mat = (bool(cm.active_mat_pos[dst]), bool(cm.active_mat_neg[dst]))
    src_mat = (bool(cm.active_mat_pos[src]), bool(cm.active_mat_neg[src]))

    shift_valid = True
    for s in (0, 1):
        if dst_nas[s] and not nas_cross[s]:
            shift_valid = False          # dst has nascent the src cannot supply
        if dst_mat[s] and not mat_cross[s]:
            shift_valid = False          # dst has MATURE the src cannot supply  ← the splice-junction case
    peel = tuple(src_mat[s] and not dst_mat[s] for s in (0, 1))
    return shift_valid, peel[0], peel[1]
```

### 1.2 Wire the diagnostic into `bp_solver._scan`

Build once, before the scans (next to `is_exon_node`):

```python
comp_model = build_component_model(chain, statics)
```

and inside the edge loop, immediately after the existing `use_shift_g` assignment:

```python
if _capture is not None:                      # inert in production (_capture is None)
    _sv, _pp, _pn = edge_vehicle(comp_model, lsrc, i)
    _capture.setdefault("_vehicle", []).append(
        {"src": int(lsrc), "dst": int(i), "shift_derived": bool(_sv),
         "shift_code": bool(use_shift_g), "peel_pos": bool(_pp), "peel_neg": bool(_pn),
         "agree": bool(_sv) == bool(use_shift_g)}
    )
```

### 1.3 Tests (`tests/calibration/test_component_model.py`)

Pin the topology table of §0.1(a) directly — these are the cases that motivated the whole model:

```python
def test_mature_crosses_exon_exon_seam_but_not_a_splice_junction(): ...
def test_tss_seam_carries_gdna_only(): ...
def test_opposite_strand_overlap_gates_per_strand(): ...
def test_shift_invalid_into_an_exon_across_a_splice_junction(): ...
def test_shift_valid_into_an_intron(): ...
```

### 1.4 Gates

Output **byte-identical** (`_capture` is `None` in production). Suites green. `ruff check`.
**Deliverable:** the disagreement census — which edges the derived rule would flip, and how many.

---

## 2. STEP 2 — flip the vehicle predicate

Replace `use_shift_g` with `edge_vehicle(...)`, and route the `c_b` peel off `peel_pos/peel_neg` instead of
`_is_bnd[i] and _ex_s`.

```python
shift_ok, peel_p, peel_n = edge_vehicle(comp_model, lsrc, i)
use_shift_g = shift_ok
_c_mat = mature_dilution[i] if (peel_p or peel_n) else 0.0
_v_mat = _var_mat[i] if (peel_p or peel_n) else 0.0
```

The `−c_b` restore stays deleted (Part E.3), and the mature **graft** (`+SPs/_esp` in `rho_pos/rho_neg`) is
removed on boundary→exon edges — under the DENSITY vehicle `md_dst` already contains the exon's mature, so
grafting it double-counts.

**Expected:** boundary→exon moves onto the density vehicle. On the closed-form fixture that takes the emitted
`f_g` from **1.000000 → 0.320762** (truth 0.320762).
**Gates:** `gdna_none` delta vs HEAD; the 4 currently-passing τ-gag/vacuous tests; toy oracle direction.

---

## 3. STEP 3 — pool densities, normalize ONCE  ⟵ the load-bearing change

Today `_comb` averages two **normalized compositions**:

```python
def _comb(am_a, ap_a, am_b, ap_b):
    pc = ap_a + ap_b
    return np.where(pc > _EPS, (ap_a * am_a + ap_b * am_b) / np.maximum(pc, _EPS), 0.0), pc
```

That is what lets a one-sided flank contribute the claim "f_g = 1". Replace with: each scan emits **per-component
imputed MASS + precision**; combine the *masses* (precision-weighted, in log space, which is what the existing
arithmetic already does correctly); normalize **once** at the end.

```python
# `_scan` returns, per component, (log_mass_mode, precision) instead of (log_fraction_mode, precision).
#   log M_c^dst = log(rho_c^src * E_c^dst)      — ABSOLUTE; no normalizer inside the scan.
# A component a flank knows nothing about arrives with precision 0 and contributes NOTHING to the pool,
# instead of contributing a zero to its share.
lm_g, lp_g = _comb(a[0], a[1], b[0], b[1])      # pooled log gDNA mass
lm_p, lp_p = _comb(a[2], a[3], b[2], b[3])
lm_n, lp_n = _comb(a[4], a[5], b[4], b[5])

# normalize ONCE, against the destination's OWN observed total where any component is unpooled
tot = np.where(lp_p > _EPS, np.exp(lm_p), 0.0) + np.where(lp_n > _EPS, np.exp(lm_n), 0.0) \
    + np.where(lp_g > _EPS, np.exp(lm_g), 0.0)
covered = _pool_covers_active(comp_model, lp_g, lp_p, lp_n)     # §4
den = np.where(covered, tot, mass_global)      # SHIFT-like when covered, DENSITY-like when not
mode_g = np.log(np.maximum(np.exp(lm_g), _EPS) / np.maximum(den, _EPS))
```

**This subsumes Step 2's vehicle choice**: "normalize by the pooled total when the pool covers the active set,
by the observed total otherwise" *is* the shift/density rule, now applied to the **union of flanks** rather than
to each message alone. That is exactly the owner's point — a single message need not supply everything; the
combined messages must.

**Gates:** this is the largest behavioural change in the plan. `gdna_none` delta, full calibration suite, the
grounded toy per-node direction, and a per-condition benchmark A/B on **stranded and unstranded** arms.

---

## 4. STEP 4 — the lazy solve gate

```python
def _pool_covers_active(cm, lp_g, lp_p, lp_n) -> np.ndarray:
    """A node may solve iff its pooled messages carry real precision on every component it actually has."""
    ok = lp_g > _EPS
    ok &= ~cm.active_nas_pos | (lp_p > _EPS)
    ok &= ~cm.active_nas_neg | (lp_n > _EPS)
    return ok
```

A node whose pool does **not** cover its active set stays at **honest precision 0** and keeps its init, to be
resolved later by the gDNA hyperprior. This is *not* the refuted degrees-of-freedom gate — that one kept the
`f_g = 1` all-gDNA init, an all-gDNA lock. Coverage is a statement about **which components were informed**,
not about counts.

**Requires** the hyperprior handoff contract: a precision-0 node must be marked, not silently defaulted.

---

## 5. STEP 5 — the coverage gradient (the remaining under-call)

After 1–4 one unknown remains: `π(E)/π(B)`, the probe-coverage ratio between an exon and its flanking
boundary. Measured **2.1–3.2× under capture** and **RNA-invariant** (<6 % across nascent 50→200), so it is
learnable from pure-gDNA structure — intergenic regions and seams — with no solve.

**Blocked on** resolving the §9.3 estimator anomaly first: the seam's density reads *below* the intergenic at
capture-off, which pure geometry does not predict. Do not build on the ratio until that is explained.

---

## 6. STEP 6 — retire the patches

Each becomes dead or absorbed once 1–5 land: the `−c_b` restore (**already removed**), the mature **graft**
(Step 2), `rho_g_cross` (**already removed**), and σ²_transfer's cliff term `(μ_dst − μ_src)²` — which is an
uncorrected **mode bias** paid as variance, and moves *into the mode* only where the ratio is identified
(pure↔pure edges). **Do not delete it before Step 5** — it is currently the only thing keeping a density-based
claim safe across an enrichment change.

---

## 7. STEP 7 — the provenance merge

`pr += S` is replicate-measurement algebra on an additive decomposition. Replace with the share-weighted
delta method (`prediction_measurement_merge_derivation.md`):

```python
def _merge_parts(parts, s2t):
    """parts = [(rho_k, var_log_k), ...] — a density that is a SUM of provenance-distinct parts.
    Var(log Σρ) = Σ (ρ_k/Σρ)²·Var(log ρ_k), then the shared transfer variance."""
    rho = math.fsum(r for r, _ in parts if r > 0.0)
    if rho <= 0.0:
        return 0.0, 0.0
    v = math.fsum((r / rho) ** 2 * vv for r, vv in parts if r > 0.0)
    return rho, (0.0 if not math.isfinite(v) else 1.0 / (v + s2t))
```

with `measured spliced → v = 1/n_spliced` (**no composition term**) and `imputed unspliced →
v = 1/n_unspliced + Var(log f_c^src)`. `_pred_precision` is the single-part case of this, so it is subsumed
rather than kept alongside.

---

## 8. ASSUMPTIONS REGISTER — write these down, test what is testable

These are load-bearing and currently unstated in code. Each needs a name, an owner, and ideally a QC signal.

| # | assumption | used by | breaks when | detectable? |
|---|---|---|---|---|
| A1 | gDNA density is genomically **uniform** (`ρ_g ≡ ρ_bg`) | the factory; every DENSITY-vehicle message | copy-number change (tumour), mappability variation | yes — variance of intergenic ρ_g across the genome; already partly in `measure_background` |
| A2 | capture coverage factorizes `e_c(x) = π(x)·γ_c` (positional × component constant) | SHIFT cancellation; factory absolute scale | probe efficiency varies by target between gDNA and RNA | partly — compare γ across panels |
| A3 | intergenic is **un-enriched** (`π = 1`) **and RNA-free** | the factory's absolute scale; the sink | readthrough, unannotated transcription, eRNA | yes — intergenic spliced counts > 0 is a red flag |
| A4 | **mature cannot cross a splice junction** unspliced | the `c_b` peel; the DENSITY vehicle into an exon | never (definitional) — **but only for a true splice junction**, not for an exon↔exon annotation seam | n/a — encoded by `mrna_active_s` in Step 1 |
| A5 | the **annotation signature** correctly names each node's active components | all of Step 1 | unannotated isoforms, wrong strand, missing exons | partly — spliced motif strand vs signature disagreement |
| A6 | nascent density is **continuous** across a junction | the SHIFT on intron↔boundary | co-transcriptional splicing gradients, 5′/3′ bias | yes — that is exactly what Step 5's residual measures |

**A4 is the one that changed this session:** the claim "a boundary crossing is mature-free" is true for a
splice junction and **false** for an exon↔exon seam. Step 1 encodes the distinction so it is no longer an
assumption at all.

---

## 9. Order, gates, and the standing rules

1 → 2 → 3 → 4 → (5 blocked on the estimator) → 6 → 7. Steps 1 and 4 are safe; 2 and 3 are behavioural.

**Every behavioural step is gated by:**
* `OMP_NUM_THREADS=1 python scripts/debug/gdna_none_guard.py` — **as a delta against HEAD** (an absolute number
  is unreadable; this is what caught the +49 % regression). Reference: HEAD = 3,766,743.
* `pytest tests/calibration/ tests/native/ -q` — 306 pass, 2 xpass.
* `scripts/debug/msg_audit.py` — per-node toward/away on the grounded toy.
* Steps 2–3 additionally: per-condition benchmark A/B on **stranded and unstranded** arms.
* Goldens are regenerated **LAST**, after the benchmark A/B — never casually.

**Standing rules:** no magic numbers (pause and discuss before any constant); develop on the grounded toy with
**injected** population priors (a toy cannot fit κ / the NPMLE / ρ_bg); owner drives commits and sequencing.

---

## 10. STEP 1 — LANDED (2026-07-22). The disagreement census.

`component_model.py` + 7 topology tests + the inert per-edge census in `_scan`. Production output
**byte-identical** (`_capture` is `None`); suites **313 pass, 2 xpass**; `ruff check` clean.

**Census on the grounded full-transcript toy (capture off + on): 56 edges, the derived rule agrees with
today's `use_shift_g` on 36 (64 %).** All 20 disagreements fall into exactly three families, and each is a
case the derivation predicted:

| family | edges | now | derived | meaning |
|---|---|---|---|---|
| **`BND → exon`** | 6 | SHIFT | **DENSITY** | the theorem: no unspliced crossing carries the exon's mature. On the closed-form fixture this is the 1.000000 → 0.320762 correction. |
| **`intgc → BND`** | 5 | DENSITY | **SHIFT** | the pure↔pure factory edge — both endpoints pure gDNA, so the shift is exact (and trivially 1→1). This is the owner's "use composition invariance between the intergenic region and the intergenic-exon boundary". |
| **`exon → BND` peel** | 2 | no peel | **peel** | a TSS/TES-side boundary: the source exon holds mature the boundary cannot, so it must be peeled. |

### 10.1 A gap the model exposed — the unmeasurable peel

The third family is new information. The peel term is `c_b = log1p(S_B/D_B)`, which needs the destination's
**spliced** count `S_B`. A TSS/TES seam is **not a splice junction**, so it has no spliced channel and
`S_B = 0` ⇒ `c_b = 0` ⇒ **the peel silently does not happen**, even though the component sets say it must.

Observed consequence (grounded toy, capture on): `exon(3) → BND(2)` — a TSS seam — emits `exp(mode_g) = 1.0000`
with no peel, while the splice-junction edge `exon(3) → BND(4)` emits 1.1067 with `c_b = +0.1014`. Both are
wrong, in opposite ways: one under-peels to exactly 1, the other over-peels past 1 (a "fraction" > 1, because
the exon's RNA is suppressed from the normalizer so there is nothing left to peel *out of*).

**So a peel needs a mature-share estimate that does not depend on a junction existing.** Options, none yet
derived: use the *source* exon's own mature-vs-nascent split (not currently identified), or treat a
component-surplus source as partial in the other direction and fall back to the DENSITY vehicle — which needs
no peel at all, since `md_dst` already excludes what does not cross. **The second is likely right and would
delete `c_b` entirely**; it needs its own derivation before Step 2 flips the predicate.

### 10.2 What Step 2 must therefore cover

1. `BND → exon` → DENSITY vehicle (+ remove the mature graft, which double-counts under that vehicle).
2. `intgc → BND` → SHIFT, and no σ²_transfer on a pure↔pure edge.
3. **Resolve the peel** per §10.1 before flipping — do not flip 1 and 2 while leaving `c_b` keyed on
   `S_B > 0`, or the TSS-side exon→boundary edges keep emitting an unpeeled 1.0.

---

## 11. The PEEL derivation — and its deletion

### 11.1 Why a peel exists at all

`c_b = log1p(S_B/D_B)` is exact. With `D_B = ρ_g + ρ_ν` (the boundary's mature-free unspliced crossing) and
`S_B = ρ_μ` (its measured spliced density):

```
    f_g^B / f_g^E  =  (ρ_g+ρ_ν+ρ_μ)/(ρ_g+ρ_ν)  =  (D_B+S_B)/D_B
```

But note **where the split comes from**: the *destination's* spliced measurement. It has to, because the
**source cannot decompose its own RNA** — an exon's solve returns a lumped `f_±` covering ν and μ together
(the simplex is 3-component over the unspliced pool; nascent-vs-mature is not one of its axes). The peel is
therefore not a correction the source applies; it is a split the source **borrows**. That is exactly why it
evaporates at a TSS/TES seam: no junction ⇒ no spliced channel ⇒ `S_B = 0` ⇒ nothing is borrowed and nothing
is peeled.

### 11.2 The third validity condition

This generalizes into the condition that was missing from §0.3:

> **A source may use the SHIFT only if it can EXPRESS the transported subset in its own simplex.**

Transporting `{g, ν}` out of an exon means splitting `f_±` into ν and μ — which the source cannot do. So the
shift is unavailable whenever the source holds a component that does **not** cross, regardless of whether a
convenient measurement happens to be nearby.

Combined with §0.3's destination-side condition, the whole rule collapses to a **single set equality**:

> ### SHIFT iff `src.active == dst.active` (per component, per strand). Otherwise DENSITY.

*Destination has something the source lacks* ⇒ the shift would assert it absent. *Source has something that
does not cross* ⇒ the shift cannot subtract it. Either way: DENSITY. **There is no peel.** Continuity is
already folded in, because a boundary's active set is itself continuity-gated (`free_s` = AND of flanks,
`mrna_active_s` = AND of flank exon bits).

### 11.3 It reproduces every case — including the two that were broken

| edge | src active | dst active | vehicle | note |
|---|---|---|---|---|
| intron → splice junction | g, ν | g, ν | **SHIFT** | equal |
| exon → splice junction | g, ν, μ | g, ν | **DENSITY** | was SHIFT + `c_b` |
| splice junction → exon | g, ν | g, ν, μ | **DENSITY** | the mature theorem |
| exon ↔ exon seam | g, ν, μ | g, ν, μ | **SHIFT** | mature crosses — equal |
| intergenic → TSS seam | g | g | **SHIFT** | the pure↔pure factory edge |
| TSS seam → exon | g | g, ν, μ | **DENSITY** | |
| **exon → TSS seam** | g, ν, μ | g | **DENSITY** | **fixes the §10.1 gap** — no `S_B` needed |
| intron → intergenic-side seam | g, ν | g | **DENSITY** | source ν cannot cross |

### 11.4 Verified on the closed-form fixture — DENSITY is exact in BOTH directions

`_mature_exon_chain`, no capture, computed against the real eff-length functions:

| direction | truth | **DENSITY** | SHIFT |
|---|---|---|---|
| boundary → exon | 0.320762 | **0.320762** ✅ | 1.000000 ❌ |
| **exon → boundary** | 1.000000 | **1.000000** ✅ | 0.428571 ❌ (unpeeled) |

The density vehicle transports gDNA alone and lets the destination's **own observed total** account for
everything that did not cross — which is precisely what a peel was hand-computing. `c_b`, `absorb_p`, and the
`S_B` dependence all become dead code.

### 11.5 ⚠️ The sequencing consequence — Step 2 and Step 5 must land TOGETHER

Set-equality moves **every exon-adjacent edge** onto the density vehicle, and the density vehicle carries the
enrichment bias `log(π(dst)/π(src))`. Measured on the grounded toy, `boundary → exon`, capture ON:

| | value | oracle |
|---|---|---|
| SHIFT (today) | 0.502 | 0.754 |
| **DENSITY (derived)** | **0.308** | 0.754 |

So under capture the structurally-correct vehicle is *further* from the oracle than the structurally-wrong one,
because today's shift is wrong in the compensating direction. Without capture, density is exact and shift is
badly wrong (§11.4). **Therefore:**

* Step 2 alone will **improve** non-capture conditions and **regress** capture-on conditions.
* Step 2 must ship **with** Step 5's coverage-gradient correction, or behind a switch A/B'd per condition.
* This is a **good** outcome for the architecture: it concentrates all residual error into **one identified
  quantity** — the coverage gradient — instead of spreading it across two compensating mechanisms.

**Revised order: 1 → 3 → 4 → 5 → 2 → 6 → 7**, i.e. pool-and-normalize-once and the gradient land *before* the
vehicle flip, so the flip is measured against a solver that can actually carry a density claim correctly.

### 11.6 Census under the set-equality rule (LANDED, still inert)

`edge_vehicle` simplified to the single set-equality test; the peel return is gone. 313 tests pass, production
**byte-identical**. Re-censused on the grounded toy: **56 edges, 24 agree with today's code (43 %)** — down
from 64 %, because the rule now also flips the `exon → boundary` direction.

| family | now | derived | note |
|---|---|---|---|
| `exon → BND` (both splice-junction and TSS-side) | SHIFT (+`c_b` on 4 of 6) | **DENSITY** | the peel deletion — including the 2 TSS-side edges where `c_b` never fired |
| `BND → exon` | SHIFT | **DENSITY** | the mature theorem |
| `intgc → BND` | DENSITY | **SHIFT** | the pure↔pure factory edge |

The flip is now **symmetric**: every exon-adjacent edge rides the density vehicle in both directions, and the
only shift edges left in this topology are `intron ↔ splice junction` and `intergenic ↔ seam` — precisely the
component-set-equal ones. `c_b` / `absorb_p` / the `S_B` dependence become dead code at Step 2.

---

## 12. REFRAME — there is one RNA. The question is ROUTING, not species.

**Owner, 2026-07-22.** Two statements that reorganize §11:

> 1. *From the source node's perspective, message propagation is trivial: send the densities and the
>    precisions. The source does not have to worry about anything else.*
> 2. *There is no separation between "mature" and "nascent" RNA — that is an artificial split. It is just RNA.
>    Some of it crosses a splice junction and departs; some continues contiguously. They are not different
>    molecules. We model RNA as two strands only because transcription is strand-specific.*

### 12.1 The model, restated

**Components are `{gDNA, RNA₊, RNA₋}`. Full stop.** There is no mature component and no nascent component.

What differs between nodes is not *what the RNA is* but **where it goes**:

* at a **splice junction**, RNA arriving at the seam either **continues contiguously** along the genome (it
  lands in the unspliced crossing) or **departs** (it lands in the spliced channel and reappears in the next
  exon, skipping the intron);
* **everywhere else** — contiguous sequence, an exon↔exon seam, mid-region — RNA simply continues;
* at a **transcript end** (TSS/TES) or a strand discontinuity, RNA is structurally **not present** beyond.

**A splice junction is the only place we observe the routing decision**, because the spliced channel is the
measurement of "how much departed". That is the entire content of what we were calling mature vs nascent.

### 12.2 The accumulator already encodes the routing

The spliced channel is deposited **one-sided, on the exon flank** (`node_geometry`: `spliced_*_left` /
`spliced_*_right`, and `spliced_* ≡ 0` for every REGION node). So "where the RNA went" is already in the data:

* **toward the exon flank** — the departed RNA *reappears there*, so the exon receives contiguous **plus**
  departed RNA. The outgoing RNA density is `(unspliced RNA + spliced)/E`.
* **toward the intron flank** — the departed RNA *skipped over it*. The intron receives the contiguous part
  only.

That asymmetry is physical, not a convention, and it is why the boundary's spliced mass must be added on one
side and not the other.

### 12.3 The transmission form

This gives the uniform edge object. A message is `(ρ_c, precision_c)` per component — sending is trivial, per
(1) — and the **edge** carries a per-component transmission `τ_c ∈ [0,1]`:

```
    M_c^dst = ρ_c^src · E_c^dst · τ_c(edge)

    τ_g      = 1                        gDNA is genomic; it always continues
    τ_r,s    = 0                        strand s not present on both sides (transcript end / discontinuity)
    τ_r,s    = 1                        contiguous, no junction  (includes an exon↔exon seam)
    τ_r,s    = measured                 at a splice junction — the spliced/unspliced split, DIRECTIONAL:
                                        toward the exon the departed RNA is received; toward the intron it is not
```

**Enrichment still cancels**, because `τ` is a pure ratio and carries no `π`. That is the important
consequence: the shift survives at splice junctions.

### 12.4 This REVISES §11 — and removes the Step 2 / Step 5 coupling

§11 concluded that every exon-adjacent edge must move to the DENSITY vehicle, because "the destination holds a
component (mature) the source cannot supply". **Under the reframe that premise is wrong: mature is not a
component, it is a flux route.** The source sends one RNA density; the edge attenuates it by `τ`. Component
sets are equal, so the SHIFT applies — and with it, enrichment invariance.

Consequences:

* §11's set-equality rule **over-fires**. It treats `mrna_active_s` as a component, which the reframe removes.
* The predicted "Step 2 regresses under capture" (density's 0.308 vs shift's 0.502 vs oracle 0.754)
  **does not apply** to junction edges, because they keep the enrichment-invariant vehicle.
* **The revised order reverts to 1 → 2 → 3 → 4 → 5 → 6 → 7.** Step 2 is now "install `τ`", not "flip
  everything to density".
* **E2 (removing `−c_b` on boundary→exon) is CONFIRMED correct** under the reframe, for a better reason than
  it was made: the graft `+SPs/_esp` *is* the routing rule — "send contiguous + departed toward the exon" —
  so the additional `−c_b` frame conversion would count the departed RNA twice.
* `+c_b` on exon→boundary survives as the τ<1 attenuation in that direction: the boundary's unspliced crossing
  receives only the contiguous share.

What genuinely remains from §11 is the **structural** part: `τ_r,s = 0` at a transcript end or strand
discontinuity, which is the TSS-seam case, and there the density vehicle is still the right answer because the
source carries no RNA at all.

### 12.5 The honest limit: routing is only observable where the junction is measured

`τ` at a splice junction is estimated from the spliced channel. A **probe-depleted junction reports `S_B ≈ 0`**
and the routing looks like "everything continued" when in truth it all departed. This is exactly the
closed-form fixture (`spliced=False`, all-mature exon, gDNA-only crossings): the truth is `τ = 0`, the
measurement says `τ = 1`, and the shift transports RNA that never arrives — emitting `f_g = 0.4286` where the
truth is 1.0.

So the failure is **not** a mis-specified model; it is a **missing measurement**, and it is the same
sensitivity already recorded as trap #3 (centred probes deplete junction reads ~300×). It should be handled by
precision — `τ`'s counting variance from `S_B` — not by switching vehicles.

### 12.6 Naming teardown (follow-up, once the model settles)

`mrna_active_s` is not "mature is active"; for a boundary it means **"this seam is not a splice junction for
strand s, so no RNA departs here"**, and for a region it means "this node is exonic on s". `free_s` is "RNA on
strand s is present on both sides". Renaming to routing language (`rna_present_s`, `no_splice_loss_s`, `τ`)
should follow the model, not lead it — deferred until Step 2 lands.
