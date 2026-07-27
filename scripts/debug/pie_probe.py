"""The relay pie-coherence probe (S1 of docs/calibration/archive/dof_pie_model_fix.md) — measures f_g+f_pos+f_neg per node.

The calibration relay (`bp_solver._scan`) maintains the per-node running belief as THREE INDEPENDENT
log-fraction Gaussians and combines each component alone, so the relayed "pie" need not be a composition. This
probe RECONSTRUCTS the running-belief pie exactly from the sweep's own `_capture` hook (NO production patch — the
reconstruction is the solver's own combine formula, cross-validated by `intron_message_trace.py`) and reports,
per node class and condition:
  * the pie sum f_g+f_pos+f_neg distribution (born/local, forward relay, backward relay);
  * the count of nodes with a component fraction > 1 (= n_src > M, the precision-inflation source);
  * the worst node's trace.

It reproduces `docs/calibration/archive/dof_pie_model_fix.md` Sec 1 and re-derives the case on the POST-GATE residual (the gate is in the
shipped tree). Theory + fix: `docs/calibration/archive/dof_pie_relay_derivation.md`.

    OMP_NUM_THREADS=1 python scripts/debug/pie_probe.py [--suite DIR] [--conditions a,b]
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from ablate_replay import build
from selfsolve_diag import _scan_and_truth, _true_fg

from rigel.calibration.node_chain import REGION
from rigel.calibration.signature import coarse_type_array
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-9
_RTYPE = {0: "intergenic", 1: "intron", 2: "exon"}


def _running_pie(cap, scan_idx):
    """The per-node running-belief pie (fbg,fbp,fbn) the sweep actually RELAYS — read directly from the
    solver's `_relay_pie` capture (forward scan = idx 0, backward = idx 1). Under the coherent (λ,θ) relay this
    is `f_g=σ(μ_λ)`, `f_pos/f_neg` from `(f_r, θ)` — coherent by construction. (The old message-reconstruction
    is stale under the coordinate relay: the running belief is no longer a log-fraction combine of the message.)"""
    fbg, fbp, fbn = cap["_relay_pie"][scan_idx]
    return np.asarray(fbg, float), np.asarray(fbp, float), np.asarray(fbn, float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", default=None)
    a = ap.parse_args()
    suite = Path(a.suite)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    cache = suite / "_selfsolve_cache"
    conds = a.conditions.split(",") if a.conditions else sorted(p.stem for p in cache.glob("*.pkl"))

    rows = []
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, Path("/tmp/rigel_selfsolve"), cache)
        ra, sub, bsub, chain, st, geom, kappa, belief, cap, cc = build(inp, index, cfg)
        kind = np.asarray(chain.kind)
        idx = np.asarray(chain.ref_idx, np.int64)
        is_r = kind == REGION
        rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
        ri = np.clip(idx, 0, rtype.shape[0] - 1)
        cls = np.where(is_r, np.array([_RTYPE.get(int(t), "?") for t in rtype])[ri], "boundary")
        sol = np.asarray(cap["solvable"], bool)
        tf_b, m_b = _true_fg(inp["boundary_pools"])
        tf_r, m_r = _true_fg(inp["region_pools"])
        bi = np.clip(idx, 0, tf_b.shape[0] - 1)
        true_fg = np.where(is_r, tf_r[ri], tf_b[bi])
        mass = np.where(is_r, m_r[ri], m_b[bi])

        born = (np.asarray(cap["fg_loc"], float) + np.asarray(cap["fp_loc"], float)
                + np.asarray(cap["fn_loc"], float))
        fg_f, fp_f, fn_f = _running_pie(cap, 0)
        fg_b, fp_b, fn_b = _running_pie(cap, 1)
        relay_f = fg_f + fp_f + fn_f
        relay_b = fg_b + fp_b + fn_b
        pie = np.where(np.abs(relay_f - 1) >= np.abs(relay_b - 1), relay_f, relay_b)

        for label, arr in (("born", born[sol]), ("relay_fwd", relay_f[sol]),
                           ("relay_worst", pie[sol])):
            q = np.percentile(arr, [1, 50, 90, 99])
            rows.append(dict(condition=c, kappa=round(kappa, 4), which=label, n=int(sol.sum()),
                             p1=round(q[0], 3), p50=round(q[1], 3), p90=round(q[2], 3),
                             p99=round(q[3], 3), MAX=round(float(arr.max()), 1),
                             incoh=round(float(np.mean(arr > 1.01)), 3)))
        n_infl = int(np.sum((np.maximum.reduce([fg_f, fp_f, fn_f]) > 1.0) & sol))
        wi = int(np.argmax(np.where(sol, np.abs(pie - 1), -1)))
        print(f"[{c}] kappa={kappa:.4f} solvable={int(sol.sum())}  n_src>M (a component frac>1): {n_infl}  "
              f"worst={wi}({cls[wi]}) pie={pie[wi]:.1f} fbg={fg_f[wi]:.2f} fbp={fp_f[wi]:.2f} "
              f"fbn={fn_f[wi]:.2f} true_fg={true_fg[wi]:.3f} mass={mass[wi]:.0f}")

    df = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print("\n" + "=" * 110)
    print(df.to_string(index=False))
    df.to_csv("/tmp/pie_probe.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
