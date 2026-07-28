"""P-2 RESIDUAL, PART B — per-node state dump, so the residual can be diffed between two TREES.

The P-2 residual is not a property of one solve; it is the DIFFERENCE between the pinned tree and the
unpinned one.  Part A falsified the off-simplex hypothesis by looking at one tree only, so look at both:
dump every live node's solved f_g, its oracle, its incoming message packet and its own belief, then diff.

    OMP_NUM_THREADS=1 python scratchpad/p2r_b_dump.py --tag p2  --refit 1     # current tree
    <revert the pin>
    OMP_NUM_THREADS=1 python scratchpad/p2r_b_dump.py --tag pin --refit 1     # pre-P-2 tree

Writes scratchpad/p2r_dump_<tag>_r<refit>.pkl.  Conditions default to the regressing stratum plus controls.
"""
from __future__ import annotations

import argparse
import dataclasses
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9

# the 4 regressing conditions named in `pin_derivation.md` §9, their two capture-OFF stratum-mates,
# and controls: the two conditions P-2 was diagnosed on, one capture-ON, one stranded.
DEFAULT = [
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",     # +0.0181 (the worst)
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",     # +0.0078
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",  # +0.0043
    "gdna_gdna5_ss_0.50_nrna_present_capture_off",    # +0.0035
    "gdna_gdna100_ss_0.50_nrna_present_capture_off",
    "gdna_gdna1_ss_0.50_nrna_present_capture_off",
    "gdna_none_ss_0.50_nrna_none_capture_off",        # control: P-2's big win
    "gdna_gdna100_ss_0.50_nrna_none_capture_on",      # control: capture ON
    "gdna_gdna100_ss_0.99_nrna_none_capture_off",     # control: stranded capture OFF
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--refit", type=int, default=1)
    ap.add_argument("--conds", nargs="*", default=DEFAULT)
    a = ap.parse_args()

    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    out = {}
    for cond in a.conds:
        inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"),
                              SUITE / "_selfsolve_cache")
        dbg: dict = {}
        calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                  np.asarray(inp["rna_fl_pmf"]),
                  dataclasses.replace(cfg.calibration, calib_refit_iters=a.refit), _debug=dbg)
        cap, chain = dbg["capture"], dbg["chain"]
        Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
        G, R = Gp + Gn, Rp + Rn
        rt, _ = _node_region_type(chain, ra)
        S = cap["_uni_static"]
        U = cap["_uni"][-1]
        st = dbg["statics"]
        out[cond] = dict(
            fo=np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan),
            fg=np.asarray(cap["f_g"], float),
            mass=np.asarray(cap["mass_global"], float),
            eff=np.asarray(cap["eff_global"], float),
            cls=np.where(np.asarray(chain.kind) != 0, 3, rt),
            ok=np.asarray(cap["solvable_mask"], bool),
            mode_g=np.asarray(cap["mode_g"], float), prec_g=np.asarray(cap["prec_g"], float),
            mode_p=np.asarray(cap["mode_p"], float), prec_p=np.asarray(cap["prec_p"], float),
            mode_n=np.asarray(cap["mode_n"], float), prec_n=np.asarray(cap["prec_n"], float),
            og=np.asarray(S["og"], float), op=np.asarray(S["op"], float), on=np.asarray(S["on"], float),
            rho0=np.asarray(S["rho_node0"], float),
            fg_init=np.asarray(cap["fg_init"], float),
            left=np.asarray(chain.left, np.int64), right=np.asarray(chain.right, np.int64),
            G=G, R=R,
            # the PER-MESSAGE mode-fusion precisions — the pin's own "supplied" test acts on these
            apg=np.asarray(U["apg"], float), app=np.asarray(U["app"], float),
            bpg=np.asarray(U["bpg"], float), bpp=np.asarray(U["bpp"], float),
            cpg=np.asarray(U["pg"], float), cpp=np.asarray(U["pp"], float),
            cpn=np.asarray(U["pn"], float),
            cg=np.asarray(U["cg"], float), cp=np.asarray(U["cp"], float), cn=np.asarray(U["cn"], float),
            c_tau=np.asarray(U["c_tau"], float), lam_msg=np.asarray(U["lam_msg"], float),
            free_pos=np.asarray(st.free_pos, bool), free_neg=np.asarray(st.free_neg, bool),
            E_g=np.asarray(S["E_g"], float), E_r=np.asarray(S["E_r"], float),
        )
        e = np.abs(out[cond]["fg"] - out[cond]["fo"])
        m = out[cond]["ok"] & np.isfinite(out[cond]["fo"]) & (out[cond]["mass"] > _EPS)
        print(f"  {cond[5:]:50s} mwae={np.average(e[m], weights=out[cond]['mass'][m]):.4f}")
    p = Path(f"scratchpad/p2r_dump_{a.tag}_r{a.refit}.pkl")
    p.write_bytes(pickle.dumps(out))
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
