"""DISSECTION 9 — HOW BIG IS THE PIN'S PARTIAL-CLAIM BRANCH? (sizing the derivation's target)

`_pin_v` substitutes the DESTINATION's own density for any component the message does not supply. That
branch is where the destination's belief re-enters as an incoming message. Before proposing a derivation we
need to know how much of the solver's message traffic actually goes through it — on every regime, not just
the zero-gDNA one where it was found.

Reported per condition, over live destinations:
  * the share of message-carrying nodes whose incoming claim is PARTIAL (>=1 component unsupplied);
  * the share of destination MASS behind those;
  * the pin's distortion factor there, 1/(1+f_g_own) when only RNA is supplied — i.e. the fraction of the
    node's own mass budget reserved for gDNA THE MESSAGE NEVER CLAIMED.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_d9_pin_scope.py
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9


def main():
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    conds = sorted(d.name for d in SUITE.iterdir()
                   if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))
    print("=== the pin's PARTIAL-CLAIM branch: how much of the message traffic goes through it ===\n")
    print(f"{'condition':46s} {'msg nodes':>10s} {'PARTIAL':>9s} {'of mass':>9s} "
          f"{'RNA-only':>9s} {'reserved for phantom gDNA':>26s}")
    tot = np.zeros(4)
    for cond in conds:
        inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
        dbg: dict = {}
        calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                  np.asarray(inp["rna_fl_pmf"]),
                  dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
        cap = dbg["capture"]
        mass = np.asarray(cap["mass_global"])
        live = mass > _EPS
        pg = np.asarray(cap["prec_g"], float)
        pp = np.asarray(cap["prec_p"], float)
        pn = np.asarray(cap["prec_n"], float)
        st = dbg["statics"]
        fp_b, fn_b = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
        # a component is "expected" only where the strand is structurally admissible
        want = np.stack([np.ones_like(fp_b), fp_b, fn_b], -1)
        have = np.stack([pg > 0, pp > 0, pn > 0], -1)
        any_msg = live & have.any(-1)
        partial = any_msg & (want & ~have).any(-1)
        rna_only = any_msg & (pg <= 0) & ((pp > 0) | (pn > 0))
        S = cap["_uni_static"]
        og, r0 = np.asarray(S["og"], float), np.asarray(S["rho_node0"], float)
        fgo = np.where(r0 > _EPS, og / np.maximum(r0, _EPS), 0.0)
        res = 1.0 - 1.0 / (1.0 + fgo)
        w = mass[rna_only]
        rr = float(np.average(res[rna_only], weights=w)) if rna_only.any() else np.nan
        print(f"{cond[5:]:46s} {int(any_msg.sum()):10d} "
              f"{100 * partial.sum() / max(int(any_msg.sum()), 1):8.1f}% "
              f"{100 * mass[partial].sum() / max(mass[any_msg].sum(), 1):8.1f}% "
              f"{int(rna_only.sum()):9d} {100 * rr:25.1f}%")
        tot += np.array([int(any_msg.sum()), int(partial.sum()), mass[partial].sum(),
                         mass[any_msg].sum()])
    print(f"\n  SUITE: {tot[1]:.0f} of {tot[0]:.0f} message-carrying nodes take the partial branch "
          f"({100 * tot[1] / max(tot[0], 1):.1f} %), carrying {100 * tot[2] / max(tot[3], 1):.1f} % "
          f"of the destination mass.")
    print("  'reserved for phantom gDNA' = 1 − 1/(1+f_g_own): the share of the node's own mass budget the")
    print("  pin hands to a gDNA component the incoming message never claimed.")


if __name__ == "__main__":
    main()
