# Session 2026-08-20 (session 2) — the dissection loop, executed on the worst IN-SCOPE scenario

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** When a finding here is settled,
    MOVE it (ROADMAP §0 for numbers, TRAPS for lessons, DESIGN/EQUATIONS for rulings/derivations)
    and delete it here in the same edit.

Substrate for every number below: **the 16-condition ladder**, certified `slot_truth.npz`
(COMPOSITION+FIELD), gDNA absolute error in FRAGMENTS. No `src/` change was made this session;
every measurement is an override arm or a read. The owner drives commits — nothing was committed.

## 0. Baseline (step ⓪)

* `preflight.py`: all ✔ — toolchain, both references, panel 16/16 certified, test chromosome 8/8,
  66 instruments import, 12/12 self-tests.
* Suite: **0 failed / 3,604 passed / 9 xfail / 3,613 collected.** The `+5` over the recorded
  3,608 is ACCOUNTED, not adjusted: `+4` = `scripts/design/exon_solvability.py` is now TRACKED
  (committed in `f295a313`; the handoff's "UNTRACKED" note is stale) and carries the standard four
  meta-cases (confirmed `--collect-only`); `+1` = `test_the_premise_fit_is_not_dominated_by_the_
  sparsest_hops` added to `tests/calibration/test_currency_policy.py` by the weighted-premise fix.
* Panel re-record (`relay_pool_ab.py --arms off on currency`, all 16): **identical to
  `ladder_weighted_premise.tsv` in every cell** — the tree has not moved since the last session;
  artifact `<ladder>/benchmark/ladder_dissect_2026-08-20.tsv`.

## 1. Steps ①–② — the worst in-scope scenario

`gdna_g98_ss_0.99_nrna_mid_capture_on` is the worst scenario inside the three IN-SCOPE strata
under **every** arm (Silent 247,388 / Relay 435,846 / currency 404,197); second is
`g50_ss_0.99_capture_on` (222,761 / 390,688 / 325,185).

## 2. Step ③ — the stage walk (`calibration_walk.py`, certified truth)

| rung | ALL live | Δ |
|---|---|---|
| A init | 212,562 | |
| B strand | 855,572 | |
| C local (refit0) | 284,615 | |
| D C+messages | 507,970 | **+223,354** |
| E C+refits(3) | 247,388 | −37,227 |
| F SHIPPED (refits+msgs) | 435,846 | +188,458 |

* **The message-free local solve is the frontier**: E (the Silent number) inherits C's error;
  messages ADD error for both policies on this condition (+223 k / +188 k) — the standing
  in-scope finding, reproduced at the stage level.
* E per class: `B exon|exon` 97,558 · `B exon|intron` 70,955 · `R exon` 69,275 · `R intron` 9,600 ·
  gene edge and intergenic **exactly 0**.

## 3. Step ③ — the object dissection, and the mechanism (not an anecdote)

Per-object rung trace (cells C/E/F on one cached payload, scored against `slot_truth`; the
override/scorer recipe is §7). Under E the error is **DIFFUSE** — top 1,000 slots hold 27 %,
top 10,000 hold 81.5 % — so this is a systematic bias, and individual objects are noise.

**The decomposition, four in-scope conditions** (populations per slot, partition asserted):
P1 = `true_f_g ≥ 0.999`; P2 = refit drove past the local answer to the 0-vertex
(`errE > errC ∧ f_gE ≤ 0.2 ∧ true ≥ 0.5`); P3 = AMBIG regions not in P1/P2; P4 = rest.

| condition | E total | P1 | P2 | P3 | P4 |
|---|---|---|---|---|---|
| g98 ss0.99 ON | 247,388 | 64.5 % | 1.2 % | 4.7 % | 29.5 % |
| g50 ss0.99 ON | 222,761 | 26.1 % | 4.0 % | 8.7 % | 61.2 % |
| g98 ss0.99 OFF | 117,533 | 60.6 % | 0.1 % | 4.0 % | 35.3 % |
| g98 ss0.50 OFF | 154,942 | 61.8 % | 0.1 % | 3.4 % | 34.8 % |

* **~90 % of the Silent frontier is P1+P4: a diffuse under-call of `f_g` at exons and
  exon-adjacent boundaries, a few percent per slot, over tens of thousands of slots.** The
  spectacular-looking top rows (refits driving evidence-free `B exon|exon` to the 0-vertex) are
  P2 = **1–4 %** — an anecdote the quantification demoted. P2 is still real and named: worst at
  g50 ss0.99 ON (141 slots, refits added +5,030), essentially absent off capture.
* The split by own-composition-evidence (the required discipline): at g98 ss0.99 ON,
  WITH = 157,495 / WITHOUT = 89,893 — most of the frontier sits ON evidence-bearing slots.
* Structurally pure-gDNA classes read exactly 0 everywhere; introns are the best-solved
  deconvolved class (0.3–0.4× the n^(−1/2) yardstick).

## 4. Step ③→④ — root cause CONFIRMED by ceiling, and the blockers re-priced on THIS tree

**Arm**: ψ's reference mean set to the condition's own truth `f_lib` (strength 1), at cell E,
via `vertex_ceiling._install_ref_exponent` (fired counters checked, belief movement checked);
base re-recorded in the same process. All 9 in-scope contaminated conditions:

| condition | E base | arm ratio | P1 ratio | P4 ratio |
|---|---|---|---|---|
| g05 ss0.50 OFF | 51,102 | 1.046× | 1.494× | 0.965× |
| g05 ss0.99 OFF | 47,771 | 1.035× | 1.699× | 0.912× |
| g05 ss0.99 ON | 94,529 | **2.650× ⛔** | 3.217× | 2.499× |
| g50 ss0.50 OFF | 151,536 | 0.981× | 0.881× | 1.008× |
| g50 ss0.99 OFF | 126,512 | 0.976× | 0.873× | 1.003× |
| g50 ss0.99 ON | 222,761 | 0.813× | 0.583× | 0.894× |
| g98 ss0.50 OFF | 154,942 | **1.001×** | **0.029×** | **2.574× ⛔** |
| g98 ss0.99 OFF | 117,533 | 0.517× | 0.047× | 1.240× |
| g98 ss0.99 ON | 247,388 | **0.322×** | **0.029×** | 0.855× |

* **Root cause confirmed**: on the worst in-scope condition the reference mean alone explains
  ~68 % of the Silent frontier, and the near-vertex population melts to 0.029×.
* **The g05 blocker STANDS on the current tree and is worse under capture** (2.650× at
  g05 ss0.99 ON vs the recorded 1.43×) — a library-wide mean is refuted, again, now at rung E.
* **g98 ss0.50 OFF is a textbook cancelling pair inside a 1.001× total**: P1 0.029× against
  P4 2.574× on the same condition. A pooled number cannot see this; the P1/P4 split can.
* ⇒ **the reference must be PER-OBJECT** — one number cannot sit where every slot's truth sits.

## 5. Step ④ — the per-object candidates PRICED, and the single-mode family REFUTED

**The unpriced arm was `ref_loc=local`** (an exon takes its density from its own flanking
`exon|intron` boundaries, so capture enrichment cancels — the owner's own proposal in
`message_notes.md` §1). Priced at cell E beside a re-recorded `pooled`, one variant install per
process (stacked installs silently let the inner variant win), fired+moved asserted:

| condition (stranded × ON) | base | `local` | `pooled` | B exon\|exon damage (local / pooled) |
|---|---|---|---|---|
| g98 ss0.99 | 247,388 | **7.176×** | 8.004× | 1,116,480 / 1,116,082 |
| g50 ss0.99 | 222,761 | **5.285×** | 5.924× | 594,984 / 595,219 |
| g05 ss0.99 | 94,529 | **5.685×** | 6.267× | 211,447 / 211,486 |

* **`local` does not rescue the stratum, and the reason is a mechanism**: the two variants'
  `B exon|exon` damage is identical to 0.04 % because the flank-density change never reaches
  those slots — deep exon-corridor boundaries have no `exon|intron` flank and fall back to the
  OFF-TARGET density, which under capture is ~113× below the enriched slot it prices. The
  reference's location then lands near the 0-vertex, and **a near-vertex location is a strong
  claim in nats** (DESIGN §6b: a location IS its strength) — at g05 it overwhelms even
  evidence-bearing slots (5.9× WITH-evidence).
* So the in-gene ANCHOR under-read (2.6–3.6×) is NOT the dominant capture-ON blocker of the
  per-object family at rung E on this tree — **the dominant failure is pricing an ENRICHED
  in-gene slot with a DEPLETED background density**, i.e. not knowing the slot's enrichment
  state. `local` fixes the anchor for slots that have the right flank and it moved almost
  nothing.

**The surviving shape is exactly ROADMAP §1 rank 3's spike-and-slab**: under capture the
per-object reference must carry the enriched/depleted MIXTURE, not a point location computed
from either mode. The open quantity is the enrichment state of the deep `exon|exon` corridor —
slots with no `exon|intron` flank and no pure-gDNA object in reach (the hop census's
"NEITHER currency serves `R exon ← B gene edge[term]` under capture" is the same fact from the
message side). This session adds the measured refutation of every single-mode alternative at
the Silent rung on the current tree.

## 6. Step ⑤ — EXECUTED (owner approved): the test chromosome grid is now g00/g05/g50/g98

`configs/test_reference.yaml` now carries `rates [0.0, 0.052632, 1.0, 49.0]` (the ladder's own
semantics); the mechanism's substrate is a **gDNA level, not a transcript structure**, and each
new rung is justified by a measured transition (P1 melts 0.029× at g98 while the same arm reads
2.65× at g05 — opposite ends, so interpolation from g50 sees neither). The loop ran end to end:
**16 scenarios simulated (18 s), both caches built, all 16 certified COMPOSITION+FIELD, and the
three-arm benchmark recorded** — `benchmark_testref_16scen_2026-08-20.tsv` beside the earlier
artifacts. ⚠ `panel.py` derives the WRONG index path for this config; pass
`--index ~/Downloads/rigel_runs/test_reference/idx` explicitly or simulate/cache refuse.

What the new rows show, and honestly do not:

* The **relay pathology reproduces sharply**: `g98 ss0.99 ON` Silent 2,048 / Relay **34,949**
  (17×) / currency 2,366 — the ladder's worst-condition sign structure, now in the 1-min loop.
* The **P1 (near-vertex reference) mechanism does NOT dominate on the toy**: at the same row
  P1 is 4.7 % of cell E's 2,048 (92 slots, well-solved at 0.10× the yardstick) against 60–65 %
  on the ladder. The ladder's frontier is tens of thousands of moderately-covered slots — a
  SCALE property a 28-transcript chromosome cannot mimic.
  `TRAPS: a-toy-and-a-panel-can-disagree-in-rank`, met again: several extended-grid rows still
  rank the policies differently than the ladder (`g05 ss0.50 OFF` toy prefers the relay; the
  ladder prefers Silent). The toy is a mechanism check, never a ranking.

## 6b. The three-policy EVIDENCE SPLIT — the imputation is still too strong beside a measurement

New measurement (shipped refits, certified truth, error in fragments; net damage = Σ(policy −
Silent) over slots where the policy is worse, minus repairs):

| condition | policy | net damage WITH own evidence | net damage WITHOUT |
|---|---|---|---|
| g98 ss0.99 ON | relay | **+123,614** | +64,844 |
| g98 ss0.99 ON | currency | **+96,705** | +60,103 |
| g50 ss0.99 OFF | relay | **+23,467** | +16,608 |
| g50 ss0.99 OFF | currency | **+58,253** | +2,472 |

Both message policies' in-scope damage is CONCENTRATED on destinations that HAVE their own
composition evidence — the imputation tramples measurements (`DESIGN.md` §0c.0c's "a regression
on strand-specific data is a PRECISION defect until proven otherwise", now with per-policy
numbers). The currency policy at `g50 ss0.99 OFF` is the sharpest case: 96 % of its net damage
is on evidence-bearing slots. ⭐ This prices ROADMAP §1 rank 1a's cheap decisive question — a
policy that speaks ONLY where the destination has no own evidence would remove ~+97 k/+124 k
(g98 ON) and ~+58 k/+23 k (g50 OFF) of net damage while keeping the deferred-stratum and
zero-control wins, which live on evidence-free populations. ⚠ The WITHOUT-evidence damage at
g98 ON (+60–65 k for BOTH policies) is a VALUE problem, not a precision one: the sources'
beliefs carry the reference bias, so the message layer transports it — fix the reference first,
then re-contrast the policies (the recorded "under a corrected reference every policy gets all
four zero controls right, Silent included" is the same fact).

## 7. How to reproduce (the recipe, no scratchpad dependency)

* Rungs: cell C = `calibrate(refits=0, messages=False)`; E = shipped refits, muted; F = shipped.
  Read `_debug["capture"]`'s `fg_init` / `fg_loc` / `f_g`; score
  `Σ|f_g − true_f_g|·mass` against the oracle cache's `slot_truth.npz` (refuse an uncertified
  table). The evidence split uses `region_init.has_own_composition_evidence(cap["_tau0_lam"])`.
* Override arms: `vertex_ceiling._install_ref_exponent(a, b)` (truth-mean: `a = f_lib`,
  `b = 1 − f_lib`) and `vertex_ceiling._install_reference_location(variant)` — ONE install per
  process, base first, `_FIRED` deltas and belief movement asserted
  (`TRAPS: an-ablation-that-never-ran`).

## 8. Open, carried forward

* The spike-and-slab reference (ranks 2+3) is now the single fix target, with this session's
  refutations as design constraints. Its gates (already written in rank 3): reduce EXACTLY to
  the shipped form at capture-OFF; both zero controls; a shuffle control.
* P2 (refits → wrong vertex at evidence-free `B exon|exon`, ≤4 %, capture-ON only) — a named
  minor mechanism; do not chase before the reference.
* `exon_solvability.py` remains stamped UNDER REVIEW (repair or delete; do not quote) — not
  addressed this session.
* `relay_pool_ab.py`'s docstring still promises a `--table pipeline` that does not exist.
