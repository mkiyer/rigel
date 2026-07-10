"""Audit assemble_priors' boundary/region tabulation — are first-class boundary (seam) nodes counted right?

The gDNA component's effective length is an IPR over per-region CONTAINED nodes + per-boundary POOLED SEAM nodes
(priors._gdna_region_node_arrays). Each node enters at density m_n/S_n; the IPR contracts toward the high-density
nodes. If the SEAM (boundary) nodes carry a LOWER density than the contained exon nodes (because the crossing
spans into the depleted intron under capture), they DILUTE the IPR → eff-len too large → gDNA competes at too-low
a density → leak. This audits: (1) mass conservation, (2) the region-vs-seam mass/support split, (3) the
per-node DENSITY of contained exon nodes vs their flanking seam nodes (the dilution), vs the ORACLE.

    OMP_NUM_THREADS=1 python scripts/debug/fb_prior_tabulation_audit.py [cond]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.priors import _gdna_region_node_arrays  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_from_signature  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

DEFAULT = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
TC = {0: "intgc", 1: "intron", 2: "exon"}
_EPS = 1e-9


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    cal = calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                    gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())

    contained = np.asarray(cal.mass_gdna_contained, float)
    side_l = np.asarray(cal.mass_gdna_left, float)
    side_r = np.asarray(cal.mass_gdna_right, float)
    reg_eff = np.asarray(cal.gdna_region_eff_len, float)
    gdna_region, participation, support_len = _gdna_region_node_arrays(cal, ra)
    seam_mass = gdna_region - contained             # pooled seam mass keyed to region r
    seam_support = support_len - reg_eff            # seam support keyed to region r
    sig = np.asarray(ra.signature)
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])

    print(f"=== {cond} : assemble_priors boundary/region tabulation audit ===")
    # (1) mass conservation
    internal_sides = float((side_r[:-1] + side_l[1:])[np.asarray(ra.ref_id)[:-1] == np.asarray(ra.ref_id)[1:]].sum())
    print("\n(1) MASS CONSERVATION:")
    print(f"    Σ contained                  = {contained.sum():>12,.0f}")
    print(f"    Σ boundary sides (all)       = {side_l.sum() + side_r.sum():>12,.0f}")
    print(f"    Σ internal seam (pooled)     = {internal_sides:>12,.0f}")
    print(f"    Σ gdna_region (cont+seam)    = {gdna_region.sum():>12,.0f}   "
          f"(= Σcont + Σinternal? {abs(gdna_region.sum() - contained.sum() - internal_sides) < 1.0})")

    # (2) region-vs-seam split of mass & support, by region type
    print("\n(2) MASS / SUPPORT split (contained 'region' node vs pooled 'seam' node), by region type:")
    print(f"    {'type':>7} {'nreg':>5} {'cont_mass':>11} {'seam_mass':>11} {'cont_supp':>11} {'seam_supp':>11}")
    for tc in (2, 1, 0):
        m = tcl == tc
        if not m.any():
            continue
        print(f"    {TC[tc]:>7} {int(m.sum()):>5} {contained[m].sum():>11,.0f} {seam_mass[m].sum():>11,.0f} "
              f"{reg_eff[m].sum():>11,.0f} {seam_support[m].sum():>11,.0f}")

    # (3) per-node DENSITY: contained exon density vs its flanking seam density (the dilution check)
    print("\n(3) NODE DENSITY m/S (the IPR weights). EXON regions — contained vs the seam keyed to them:")
    ex = (tcl == 2) & (gdna_region > 1.0)
    cont_dens = np.where(reg_eff > _EPS, contained / np.maximum(reg_eff, _EPS), 0.0)
    seam_dens = np.where(seam_support > _EPS, seam_mass / np.maximum(seam_support, _EPS), 0.0)
    exm = ex & (reg_eff > _EPS) & (seam_support > _EPS)
    cd, sd = cont_dens[exm], seam_dens[exm]
    print(f"    exon CONTAINED density: median {np.median(cd):.2f}  (mass-wt mean {np.average(cd, weights=contained[exm]+_EPS):.2f})")
    print(f"    exon SEAM     density: median {np.median(sd):.2f}  (mass-wt mean {np.average(sd, weights=seam_mass[exm]+_EPS):.2f})")
    ratio = np.median(sd) / max(np.median(cd), _EPS)
    print(f"    → seam/contained density ratio (median) = {ratio:.2f}   "
          f"({'SEAMS DILUTE the IPR → eff-len too large → gDNA under-competes' if ratio < 0.7 else 'comparable'})")

    # (4) oracle: true contained gDNA density (payload_gdna) vs the solver's contained density
    from rigel.calibration.substrate import CalibrationSubstrate
    g_or = np.asarray(CalibrationSubstrate.from_payload(blob["payload_gdna"], ra).contained.mass_unspliced, float)
    or_dens = np.where(reg_eff > _EPS, g_or / np.maximum(reg_eff, _EPS), 0.0)
    print("\n(4) vs ORACLE contained gDNA density (exon regions):")
    print(f"    solver contained density median {np.median(cont_dens[exm]):.2f}   "
          f"oracle median {np.median(or_dens[exm]):.2f}")


if __name__ == "__main__":
    main()
