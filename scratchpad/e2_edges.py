"""ENRICHMENT-RATIO INTERROGATION part 2 — the owner's concern stated at the EDGE, not the node.

"A sparse node in the denominator makes the enrichment ratio unstable" is a claim about EDGES: a node with
1 count carries almost no mass itself, so a per-node mass census understates its influence — what matters is
how much DESTINATION mass is scaled by a ratio whose source is sparse.

Also sizes the unused FOURTH information source: a per-node MEAN FRAGMENT LENGTH would give a
method-of-moments f_g (E[l] = f_g*mu_g + (1-f_g)*mu_r), with lambda-axis precision n*I_FL*(f_g(1-f_g))^2.
Compared here against the composition precision the solver actually has (tau_own).

    OMP_NUM_THREADS=1 python scratchpad/e2_edges.py [cond ...]
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

for cond in CONDS:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    pg = np.asarray(inp["gdna_fl_pmf"], float)
    pr = np.asarray(inp["rna_fl_pmf"], float)
    L = max(pg.size, pr.size)
    pg, pr = np.pad(pg, (0, L - pg.size)), np.pad(pr, (0, L - pr.size))
    pg, pr = pg / max(pg.sum(), _EPS), pr / max(pr.sum(), _EPS)
    mix = 0.5 * pg + 0.5 * pr
    ok_l = mix > 1e-12
    I_FL = float(np.sum((pg[ok_l] - pr[ok_l]) ** 2 / mix[ok_l]))

    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    us = cap["_uni_static"]
    kind = np.asarray(chain.kind)
    isr = kind == REGION
    n = kind.shape[0]
    M, lv = us["M"], us["logvar_tot"]
    n_node = np.where(isr, us["n_unspl_l"], us["n_unspl_l"] + us["n_unspl_r"])
    li, ri = us["left"], us["right"]
    tau_own = us["tau_own"]
    fg = np.asarray(cap["fg_loc"], float)
    rt, _ = _node_region_type(chain, ra)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])

    print(f"\n{'=' * 108}\n{cond[5:]}\n{'=' * 108}")

    # ── the EDGE view of the sparse-denominator concern ──
    print("  [A] EDGES, by the SOURCE node's count — how much DESTINATION mass is scaled by a sparse ratio")
    print(f"      {'source n':<12}{'edges':>8}{'dst mass':>14}{'share':>8}{'sd(log r)':>11}{'src class':>28}")
    src_n, dst_m, sd_r, src_c = [], [], [], []
    for nbr in (li, ri):
        v = nbr >= 0
        s = np.clip(nbr, 0, n - 1)
        src_n.append(n_node[s][v]); dst_m.append(M[v])
        sd_r.append(np.sqrt(np.maximum(lv[v] + lv[s][v], 0.0)))
        src_c.append(cls[s][v])
    src_n = np.concatenate(src_n); dst_m = np.concatenate(dst_m)
    sd_r = np.concatenate(sd_r); src_c = np.concatenate(src_c)
    tot = dst_m.sum()
    for lo, hi, lab in ((0, 1, "n <= 1"), (2, 5, "n = 2-5"), (6, 20, "n = 6-20"),
                        (21, 100, "n = 21-100"), (101, 10**12, "n > 100")):
        m = (src_n >= lo) & (src_n <= hi)
        if not m.any():
            continue
        cc = ", ".join(f"{c} {np.mean(src_c[m] == c):.0%}" for c in ("exon", "intron", "boundary")
                       if np.mean(src_c[m] == c) > 0.05)
        print(f"      {lab:<12}{int(m.sum()):>8}{dst_m[m].sum():>14,.0f}{dst_m[m].sum() / tot:>8.1%}"
              f"{np.median(sd_r[m]):>11.3f}{cc:>28}")

    # ── the unused FOURTH source: a per-node mean fragment length ──
    print(f"\n  [B] THE UNUSED FL SIGNAL (per-fragment Fisher I_FL = {I_FL:.3f} about f_g)")
    print(f"      {'class':<12}{'nodes':>7}{'med n':>8}{'tau_own now':>13}{'tau_FL would be':>17}"
          f"{'ratio':>10}{'tau_own=0':>11}")
    fgc = np.clip(fg, 1e-3, 1 - 1e-3)
    tau_fl = n_node * I_FL * (fgc * (1.0 - fgc)) ** 2
    live = (M > _EPS) & (n_node > 0)
    for c in ("exon", "intron", "boundary"):
        m = live & (cls == c)
        if m.sum() < 3:
            continue
        t0 = tau_own[m]
        print(f"      {c:<12}{int(m.sum()):>7}{np.median(n_node[m]):>8.0f}"
              f"{np.median(t0):>13.4f}{np.median(tau_fl[m]):>17.3f}"
              f"{np.median(tau_fl[m]) / max(np.median(t0), 1e-6):>10.1f}{np.mean(t0 <= _EPS):>11.0%}")
