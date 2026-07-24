# The §6B Solve-Gate — the DOF criterion is VALID for regions; my first evaluation was flawed

**Status (2026-07-21): ACTIVE — the concept holds for regions; the boundary rule + the metric + the prior-source
need fixing.** An earlier draft of this note declared the solve-gate "refuted"; **that conclusion was wrong** — it
rested on a flawed metric and a gate bug. The corrected picture is below.

---

## 1. The idea (§6B)

Pass-0 today solves every **structurally** solvable node (`solvable = (fp|fn) & mass>0`). §6B proposed a stricter
**DOF** gate: a node is solvable iff every **free axis** is *identified* (has ≥1 nonzero-precision source among
{strand tilt, messages, prior}). An **un**identified node should not be solved — its forced solve uses the
reference/geometry as if it were evidence, producing a value uncorrelated with truth that then pollutes neighbours
and the hyperprior fit.

## 2. Why the first evaluation was WRONG (owner's catch)

- **The metric counted arbitrary defaults.** `pass0_bench` mwae is `Σ mass·|f_g − oracle|` over *all* nodes,
  including withheld ones — but a withheld node's `f_g` is the **arbitrary init default** (`f_g=1`), not a
  prediction. Penalising a withheld node for its default measures luck, not error. *"It only matters when we
  actually solve and precision > 0."*
- **The refit test was broken two ways.** (a) The hyperprior is fit on **solved single-strand regions only**
  (`_fit_gdna_hyperprior`: `single | gonly`, explicitly **no AMBIG, no boundary** — non-circular), so the
  deferred AMBIG nodes do **not** pollute the fit. (b) The refit re-runs `node_sweep`, and my DOF gate **omitted
  the prior as an identification source**, so it re-skipped the deferred nodes and the hyperprior **never got to
  resolve them**. The "refit=1 regression" was that bug + the metric artifact, not a real result.

## 3. The correct test — CORRELATION, not error (does the forced solve mean anything?)

For each node type, the correlation of the **forced** solve `f_g` with the oracle, split by the DOF verdict
(pooled over the 20 unstranded `ambig_dense_10mb` scenarios, `scripts/debug/dof_nodetype.py`):

| node type | DOF | mass% | **corr(f_g, oracle)** | reading |
|---|---|---|---|---|
| single-strand region | SOLVABLE | 64.2 % | **0.633** | meaningfully solved ✓ |
| single-strand region | **UNSOLVABLE** | 1.8 % | **−0.037** | **coin-flip ⇒ correctly withheld** ✓ |
| AMBIG region | SOLVABLE | 12.8 % | **0.689** | meaningfully solved ✓ |
| AMBIG region | **UNSOLVABLE** | 1.9 % | **−0.066** | **coin-flip ⇒ correctly withheld** ✓ |
| boundary | SOLVABLE | 4.1 % | 0.129 | coin-flip — **wrongly solved** ✗ |
| boundary | UNSOLVABLE | 1.6 % | 0.675 | meaningful — **wrongly withheld** ✗ |

**Conclusion:** for **regions (77 % of the mass)** the DOF criterion is **correct** — the nodes it withholds are
genuinely uncorrelated with truth (their forced solve is a coin-flip), and the ones it solves are well-correlated.
So the solve-gate concept is **valid for regions**; withholding them is right, and their default `f_g` must simply
not be counted as error.

## 4. The real bug — the boundary rule is INVERTED

For **boundaries** the region-centric criterion (`λ`/`θ` identified from `tau0_lam` + `prec_g/p/n`, node class
from `fp+fn`) is *inverted*: it calls coin-flip boundaries "solvable" and meaningful boundaries "unsolvable."
Boundaries are not regions — they own one-sided spliced mass and their identification is different. **The boundary
solvability rule must be re-derived** (its own node-type analysis) before the gate can ship.

## 5. The corrected plan

1. **Metric:** evaluate the solve-gate on **correlation / precision-weighted** error over *solved* nodes — never
   the mass-weighted error that counts withheld defaults.
2. **Prior as a source:** the DOF gate must count the fitted hyperprior as an identification source, so the refit
   *resolves* the deferred nodes (not re-skips them). Only then is the refit=1 test meaningful.
3. **Boundary rule:** re-derive boundary solvability (node-type-by-node-type), fixing the inversion in §4.
4. **Then** re-run the ON/OFF comparison (correlation-based, refit-aware) to confirm the region win and the fixed
   boundary behaviour, and to quantify the confident-wrong reduction.

*(No solve-gate code currently ships — the flag-gated prototype was reverted; the criterion is validated offline
via `dof_nodetype.py`. Implementation resumes with the boundary rule + the correct metric.)*
