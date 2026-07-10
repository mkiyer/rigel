"""Phase-3 end-to-end gate: run the PRODUCTION calibrate (enrichment-aware ê wired into node_sweep) on the
flagship and measure whether the AMBIG-exon contained gDNA recovers vs the pre-enrichment −107,509 under-call.
Also checks the recovery is not a runaway over-call (proxy for "f_g_var doesn't collapse / no over-pinning").

    OMP_NUM_THREADS=1 python scripts/debug/enrichment_phase3_check.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
TC = {0: "intergenic", 1: "intron", 2: "exon"}


def run(cond):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    cal = calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                    gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    g_pr = np.asarray(cal.mass_gdna_contained, float)
    g_or = np.asarray(CalibrationSubstrate.from_payload(blob["payload_gdna"], ra).contained.mass_unspliced, float)
    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])
    print(f"\n=== {cond} ===")
    print(f"  total contained gDNA: prod {g_pr.sum():,.0f}  oracle {g_or.sum():,.0f}  net {g_pr.sum()-g_or.sum():+,.0f}")
    print(f"  {'class':>16} {'n':>4} {'prod':>11} {'oracle':>11} {'net':>11}")
    rows = {}
    for s, t in [(3, 2), (2, 2), (1, 2), (1, 1), (2, 1), (0, 0)]:
        m = (scl == s) & (tcl == t)
        if not m.any():
            continue
        net = float(g_pr[m].sum() - g_or[m].sum())
        rows[(s, t)] = net
        print(f"  {SC[s]+'.'+TC[t]:>16} {int(m.sum()):>4} {g_pr[m].sum():>11,.0f} {g_or[m].sum():>11,.0f} {net:>+11,.0f}")
    return rows


def main():
    print("PHASE 3 end-to-end gate — production calibrate WITH enrichment-aware ê")
    conds = sys.argv[1:] or ["gdna_gdna300_ss_0.99_nrna_none_capture_on"]
    print("\n=== GATE (AMBIG.exon prior net err; pre-enrichment baseline ss0.99 was −107,509) ===")
    for c in conds:
        amb = run(c).get((3, 2), float("nan"))
        ok = abs(amb) < 40_000
        print(f"  {c}: AMBIG.exon net {amb:+,.0f}  {'RECOVERED' if ok else 'CHECK'}")


if __name__ == "__main__":
    main()
