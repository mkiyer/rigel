"""DISSECTION 8 — THE BUG, ISOLATED AND CONFIRMED CAUSALLY: `_pin_v` feeds the node its own guess back.

THE ARITHMETIC (verified to 2.1e-16 on the median FP exon, exactly on 252/433).

`bp_solver._pin_v` rescales an incoming claim so that `sum_c rho_c * E_c = M`, and for a component the
message does NOT supply it substitutes the destination's OWN density into the budget:

    sg = where(pg > 0, g, og)          # <-- og is the DESTINATION's own gDNA density
    sp = where(pp > 0, p, op)
    sn = where(pn > 0, n, on)
    k  = M / (sg*E_g + (sp+sn)*E_r)
    return g*k, p*k, n*k

On the false-positive exons the message supplies **RNA only** (measured: gDNA-message precision is 0 on
433/433 of them). With E_g == E_r == E and an incoming RNA density equal to the node's own total density
rho_tot = M/E, the pin therefore returns

    n_pinned = n * M / (og*E + n*E) = rho_tot^2 / (rho_tot + og) = rho_tot * 1/(1 + f_g_own)

so the delivered RNA fraction is exactly **1/(1 + f_g_own)** — a function of the DESTINATION'S OWN
strand-only self-solve, with no information from the source in it at all. Measured claims
{0.6484, 0.6623, 0.6711, 0.6860} against own f_g {0.4577, 0.4902, 0.5098, 0.5423}: the four discrete
strand-only values on unstranded data.

**On a zero-gDNA library this is self-fulfilling.** The strand channel has no information, so the node's own
self-solve sits at the uninformative f_g ~ 0.5; the pin then RESERVES half the mass budget for that
imaginary gDNA; the RNA message is shrunk by 1/(1+0.5); and the solve reads back ~34 % gDNA — which is the
number it started from. The docstring's intent is sound for the opposite case ("a message carrying gDNA only
still delivers f_g < 1"), but in the RNA-only direction it injects the destination's own prior belief into
the message it is supposed to be RECEIVING.

THE CONTROL that makes it causal rather than correlational is the STRANDED arm. There the strand likelihood
resolves the node's own composition, so f_g_own collapses to 0.013 — and the pin's reservation collapses with
it, from 33.6 % to 1.2 %, taking the false-positive rate from 29.3 % to 1.4 %. **The pin's reservation IS the
false-positive rate.** Nothing else about those libraries differs.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_d8_pin_bug.py
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = [
    "gdna_none_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.50_nrna_none_capture_on",
    "gdna_none_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
]


def run(cond):
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    try:
        calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                  np.asarray(inp["rna_fl_pmf"]),
                  dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    finally:
        pass
    return dbg


def main():
    print("=== THE PIN'S SUBSTITUTION, on the nodes it fires on ===\n")
    print(f"{'condition':44s} {'mass':>12s} {'true gDNA':>11s} {'called gDNA':>12s} {'FP %':>7s} "
          f"{'msgs RNA-only':>14s}")
    for cond in CONDS:
        dbg = run(cond)
        cap = dbg["capture"]
        chain = dbg["chain"]
        mass = np.asarray(cap["mass_global"])
        fg = np.asarray(cap["f_g"])
        live = mass > _EPS
        pg = np.asarray(cap["prec_g"], float)
        pp = np.asarray(cap["prec_p"], float)
        pn = np.asarray(cap["prec_n"], float)
        rna_only = live & (pg <= 0) & ((pp > 0) | (pn > 0))
        S = cap["_uni_static"]
        og = np.asarray(S["og"], float)
        r0 = np.asarray(S["rho_node0"], float)
        # the factor the pin applies when only RNA is supplied: 1/(1+f_g_own)
        fac = 1.0 / (1.0 + np.where(r0 > _EPS, og / np.maximum(r0, _EPS), 0.0))
        isr = np.asarray(chain.kind) == REGION
        _ = isr
        called = float((fg * mass)[live].sum())
        print(f"{cond[5:]:44s} {mass[live].sum():12,.0f} {'0' if 'none_' in cond[:10] else '-':>11s} "
              f"{called:12,.0f} {100 * called / max(mass[live].sum(), 1):7.1f} "
              f"{int(rna_only.sum()):9d} nodes")
        if rna_only.any():
            w = mass[rna_only]
            print(f"{'':44s}   -> on those, the pin shrinks the RNA claim by "
                  f"{np.average(fac[rna_only], weights=w):.4f} "
                  f"(= 1/(1+f_g_own), f_g_own={np.average((og / np.maximum(r0, _EPS))[rna_only], weights=w):.4f}), "
                  f"reserving {100 * (1 - np.average(fac[rna_only], weights=w)):.1f} % of the budget for gDNA "
                  f"the message never claimed")


if __name__ == "__main__":
    main()
