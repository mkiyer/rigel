"""Is the RNA-free-source λ-message handled CORRECTLY, or only where the destination happens to have evidence?

THE GEOMETRY THE OWNER ASKED ABOUT.  intergenic | seam | EXON.
  * the seam (TSS/TES) is a G1 node: RNA cannot cross a gene boundary, so it carries NO RNA — structurally.
  * the exon carries gDNA + RNA.
So the message arriving at the exon has RNA density exactly 0 while the exon's own RNA is > 0.

WHAT THE COMBINE DOES WITH THAT.  The composition (λ) message's mode is built as
    lam_msg = mo_g - mo_R ,   mo_R = log(max(cR*E_r/M, _EPS))
so when the fused message RNA density cR is exactly 0, mo_R hits the _EPS FLOOR and lam_msg becomes
"f_g = 1" — a maximally confident assertion of pure gDNA that is a NUMERICAL ARTIFACT of the floor, not a
belief anybody holds. Its precision c_tau, meanwhile, is real: it was accumulated along the relay from nodes
that never said any such thing. That is a MODE/PRECISION mismatch, and it is what dragged node 1167
(self-solve 0.643 = oracle 0.645) to f_g = 1.000.

THE QUESTION THIS SCRIPT ANSWERS. The DL term kills that message — but only via `contradicted & isfinite(v_own)`,
i.e. only where the destination HAS composition evidence (tau_own > 0). On unstranded data tau_own = 0 for
every exon and boundary (the kappa=1/2 deadband), so the artifact should survive there untouched. Measure:
how many nodes receive a floor-artifact lambda message, how much mass they carry, and what the split is
between KILLED (v_own finite) and SURVIVING (v_own = inf).

    OMP_NUM_THREADS=1 python scratchpad/seam_lambda_audit.py
"""

from __future__ import annotations

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
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}

CONDS = [
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",   # STRANDED   (tau_own > 0 on most nodes)
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",   # UNSTRANDED (tau_own = 0 except introns)
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    uni, us = cap["_uni"][-1], cap["_uni_static"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fg, mass = np.asarray(cap["f_g"]), np.asarray(cap["mass_global"])
    tau_own = np.asarray(us["tau_own"])
    own_R = np.asarray(us["op"]) + np.asarray(us["on"])
    msg_R = np.asarray(uni["cp"]) + np.asarray(uni["cn"])
    c_tau, lam = np.asarray(uni["c_tau"]), np.asarray(uni["lam_msg"])
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    ok = np.isfinite(fo) & (mass > _EPS)

    # the floor artifact: the message says "no RNA at all" while the node's own solve says it HAS RNA,
    # and the lambda stream still carries precision.
    artifact = (msg_R <= _EPS) & (own_R > _EPS) & (c_tau > _EPS)
    killed = artifact & (tau_own > _EPS)        # DL kills it: v_own finite
    survives = artifact & (tau_own <= _EPS)     # DL is inert: v_own = inf

    print(f"\n══════ {cond}")
    print(f"  nodes with mass+oracle: {int(ok.sum()):,}    mwae={np.average(np.abs(fg-fo)[ok], weights=mass[ok]):.4f}")
    print(f"  {'':<34}{'nodes':>8}{'mass':>14}{'err-mass':>12}{'mwae':>9}{'mean f_g implied':>18}")
    for lab, m in (("λ-msg is a FLOOR ARTIFACT", artifact), ("   → KILLED by DL (τ_own>0)", killed),
                   ("   → SURVIVES (τ_own=0)", survives)):
        m = m & ok
        if not m.any():
            print(f"  {lab:<34}{0:>8}")
            continue
        implied = 1.0 / (1.0 + np.exp(-lam[m]))
        print(f"  {lab:<34}{int(m.sum()):>8}{mass[m].sum():>14,.0f}"
              f"{(np.abs(fg-fo)*mass)[m].sum():>12,.0f}"
              f"{np.average(np.abs(fg-fo)[m], weights=mass[m]):>9.4f}{np.mean(implied):>18.3f}")
    s = survives & ok
    if s.any():
        print(f"  SURVIVING set: mean oracle f_g={np.average(fo[s], weights=mass[s]):.3f}  "
              f"mean solved f_g={np.average(fg[s], weights=mass[s]):.3f}  "
              f"(the artifact asserts f_g→1, so an over-estimate is its signature)")
        by = {}
        for i in np.flatnonzero(s):
            k = CLS[int(rt[i])] if kind[i] == REGION else "boundary"
            by[k] = by.get(k, 0) + mass[i]
        print("  SURVIVING mass by class: " + ", ".join(f"{k}={v:,.0f}" for k, v in sorted(by.items())))
