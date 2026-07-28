"""SUBSTRATE ACCURACY DUMP — the hyperprior's own training set, per arm.

`pass0_error_table.py` writes /tmp/pass0_state.npz but (a) filters on ``solvable`` (so intergenic is gone)
and (b) records no single/gonly split — so it cannot answer "how accurate is the population the hyperprior
is actually FITTED ON".  This dumps exactly the fields `calibrate._fit_gdna_hyperprior` selects on:

    isr    = chain.kind == REGION
    single = free_pos ^ free_neg          gonly = ~free_pos & ~free_neg      ambig = free_pos & free_neg
    live   = (eff_global > 1e-9) & (mass_global > 1e-12)
    SEL(EM path) = live & isr & (single|gonly) & ~intergenic   (the `& ~intergenic` fires iff background)

plus, per node, mass_global / eff_global / f_g / var_gdna / oracle f_g and the oracle gDNA mass, so both
    f_g accuracy   |f_g - oracle|
    density accuracy |log(rho_hat / rho_oracle)|,  rho = f_g*M/E_g
can be computed on that exact set.  It ALSO fits the hyperprior itself, twice per condition — once from the
arm's belief, once from the ORACLE f_g on the same nodes — and stores both grids so the achievable prior can
be compared against the perfect one.

    OMP_NUM_THREADS=1 python scratchpad/subacc_dump.py --arm noown2 --out /tmp/subacc_noown2.npz
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9

ap = argparse.ArgumentParser()
ap.add_argument("--arm", required=True)
ap.add_argument("--refit", type=int, default=0)
ap.add_argument("--out", required=True)
a = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

KEYS = ("cond", "mass", "eff", "fg", "var", "fo", "gor", "isr", "single", "gonly", "ambig",
        "live", "solvable", "rt", "intergenic")
C: dict[str, list] = {k: [] for k in KEYS}
PRI: dict[str, np.ndarray] = {}

for i, cond in enumerate(conds):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=a.refit), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    belief = dbg["belief"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)

    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    fg = np.asarray(cap["f_g"], float)
    var = np.asarray(cap["var_g"], float)
    isr = np.asarray(chain.kind) == REGION
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    rt, _ = _node_region_type(chain, ra)
    sig = np.asarray(ra.signature)
    ridx = np.asarray(chain.ref_idx, np.int64)
    intergenic = isr & (ridx < sig.shape[0]) & (sig[np.clip(ridx, 0, sig.shape[0] - 1)] == 0)
    live = (eff > 1.0e-9) & (mass > 1.0e-12)

    C["cond"] += [cond] * mass.shape[0]
    C["mass"] += mass.tolist()
    C["eff"] += eff.tolist()
    C["fg"] += fg.tolist()
    C["var"] += var.tolist()
    C["fo"] += fo.tolist()
    C["gor"] += (Gp + Gn).tolist()
    C["isr"] += isr.tolist()
    C["single"] += (fp ^ fn).tolist()
    C["gonly"] += (~fp & ~fn).tolist()
    C["ambig"] += (fp & fn).tolist()
    C["live"] += live.tolist()
    C["solvable"] += np.asarray(cap["solvable"], bool).tolist()
    C["rt"] += np.where(isr, rt, -1).astype(np.int64).tolist()
    C["intergenic"] += intergenic.tolist()

    # ── the hyperprior itself: fitted from the arm's belief vs from the ORACLE on the same nodes ──────────
    bg = dbg["calibration_priors"].background
    bw = float(cfg.calibration.npmle_bandwidth)

    def _fit(b):
        return calmod._fit_gdna_hyperprior(
            chain, b, st, ra, mass, eff, background=bg, bandwidth=bw,
            additive=bool(cfg.calibration.gdna_prior_additive),
        )

    fo_safe = np.where(np.isfinite(fo), fo, fg)
    orc_same = dataclasses.replace(belief, f_g=fo_safe)               # oracle mode, arm's own widths
    orc_sharp = dataclasses.replace(belief, f_g=fo_safe, var_gdna=np.zeros_like(var))  # perfect prior
    bel_sharp = dataclasses.replace(belief, var_gdna=np.zeros_like(var))  # arm's mode, NO declared width
    for tag, b in (("belief", belief), ("orc", orc_same), ("orcsharp", orc_sharp),
                   ("beliefsharp", bel_sharp)):
        p = _fit(b)
        if p is None:
            continue
        PRI[f"{cond}|{tag}|log_rho"] = p.log_rho
        PRI[f"{cond}|{tag}|logP"] = p.logP
        PRI[f"{cond}|{tag}|w"] = p.weights
        PRI[f"{cond}|{tag}|ncell"] = np.array([p.n_cells], float)
    print(f"  [{i + 1:>2}/{len(conds)}] {cond}", flush=True)

d = {k: np.asarray(v) for k, v in C.items()}
d["_arm"] = np.array([a.arm])
np.savez_compressed(a.out, **d, **PRI)
print(f"wrote {a.out}  nodes={d['mass'].shape[0]:,}")
