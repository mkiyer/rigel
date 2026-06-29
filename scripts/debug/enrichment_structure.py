"""Ground the hybrid-capture derivation: measure the gDNA enrichment field E(node) = ρ_g / floor
across node CLASSES, where floor = the depleted off-target gDNA density (mean ρ_g over intergenic +
intron regions). Tests: (1) is the floor well-defined + tight? (2) is E bimodal under capture (E≈1
off-target seeds vs E≫1 on-target exons) and ~1 off-capture? (3) do single-strand exons SPAN the
enriched range (→ can train ê)? (4) are AMBIG exons at the enriched level (→ need the enrichment prior,
can't impute from depleted seams)?

    OMP_NUM_THREADS=1 python scripts/debug/enrichment_structure.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration.effective_length import region_eff_length  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402

_EPS = 1e-9
CONDS = [
    ("FLAGSHIP cap-ON ss99", "gdna_gdna300_ss_0.99_nrna_none_capture_on"),
    ("cap-OFF   ss99",        "gdna_gdna300_ss_0.99_nrna_none_capture_off"),
    ("UNSTRANDED cap-ON ss50","gdna_gdna300_ss_0.50_nrna_none_capture_on"),
]


def main():
    for label, cond in CONDS:
        index, blob = build_or_load_cache(cond, False)
        ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
        sub_g = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
        g_or = np.asarray(sub_g.contained.mass_unspliced, float)
        reg_eff_g = region_eff_length(np.asarray(ra.region_size_bp, float), blob["gdna_pmf"])
        rho_g = np.where(reg_eff_g > _EPS, g_or / np.maximum(reg_eff_g, _EPS), 0.0)  # TRUE gDNA density

        sig = np.asarray(ra.signature)
        scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
        tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])  # 0=intergenic,1=intron,2=exon
        sz = np.asarray(ra.region_size_bp, float)
        big = sz > 50.0  # ignore tiny regions (unstable density)

        # FLOOR = mean ρ_g over intergenic + intron regions (the depleted off-target level).
        floor_mask = ((tcl == 0) | (tcl == 1)) & big & (reg_eff_g > _EPS)
        floor = float(np.sum(g_or[floor_mask]) / max(np.sum(reg_eff_g[floor_mask]), _EPS))
        E = rho_g / max(floor, _EPS)

        classes = [
            ("intergenic", (tcl == 0) & big),
            ("intron    ", (tcl == 1) & big),
            ("exon SS   ", (tcl == 2) & ((scl == 1) | (scl == 2)) & big),
            ("exon AMBIG", (tcl == 2) & (scl == 3) & big),
        ]
        print(f"\n===== {label} =====   gDNA floor (intergenic+intron) = {floor:.4f} /bp")
        print(f"  {'class':<11} {'n':>4} {'ρ_g med':>8} {'E med':>7} {'E p10':>7} {'E p90':>7} {'E max':>7}")
        for nm, m in classes:
            m = m & (g_or > 0)
            if not m.any():
                print(f"  {nm:<11}    0")
                continue
            e = E[m]
            print(f"  {nm:<11} {int(m.sum()):>4} {np.median(rho_g[m]):>8.3f} {np.median(e):>7.1f} "
                  f"{np.percentile(e,10):>7.1f} {np.percentile(e,90):>7.1f} {e.max():>7.1f}")


if __name__ == "__main__":
    main()
