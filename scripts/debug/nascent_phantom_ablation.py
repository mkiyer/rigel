"""Dive: WHY does the nested symmetric gDNA refit inject MORE phantom gDNA in zero-gDNA + nascent?

Zero-gDNA conditions: ALL assigned gDNA is phantom (truth=0). Nascent (nrna=rnd) regressed +40% under the
nested loop. Hypothesis: refitting the gDNA var~mean on ALL converged nodes (symmetric) learns a gDNA spread
that includes nascent-RNA-misread-as-gDNA (introns carry unspliced nascent), so the gDNA messages get stronger
and propagate phantom. Ablation isolates the cause:

  modes:
    both_K1 / both_K4     : refit BOTH var~mean (current); 1 vs 4 outer iterations (does phantom GROW w/ refits?)
    rna_only_K4           : refit RNA only, gDNA FROZEN at the seed bootstrap (the firewall) — does it fix it?

Reports total phantom gDNA (= assigned, truth 0) broken down by region type × strand class.

    OMP_NUM_THREADS=1 python scripts/debug/nascent_phantom_ablation.py [cond]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.config import CalibrationConfig  # noqa: E402

DEFAULT = "gdna_none_ss_0.50_nrna_rnd_capture_off"
SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
TC = {0: "intergenic", 1: "intron", 2: "exon"}


def _freeze_gdna():
    """Monkeypatch so the gDNA var~mean stays at the seed bootstrap (no symmetric refit). Returns a restore fn."""
    stash = {}
    orig_seed = bp._gdna_seed_estimate
    orig_fit = bp.fit_gdna_varmean

    def seed_patch(*a, **k):
        rho, gvm, vm = orig_seed(*a, **k)
        stash["gvm"] = gvm
        return rho, gvm, vm

    bp._gdna_seed_estimate = seed_patch
    bp.fit_gdna_varmean = lambda *a, **k: stash["gvm"]

    def restore():
        bp._gdna_seed_estimate = orig_seed
        bp.fit_gdna_varmean = orig_fit

    return restore


def run(cal_kwargs, max_outer, freeze_gdna):
    cfg = dataclasses.replace(CalibrationConfig(), sweep_max_outer=max_outer)
    restore = _freeze_gdna() if freeze_gdna else None
    try:
        cal = calibrate(config=cfg, **cal_kwargs)
    finally:
        if restore:
            restore()
    return cal


def report(cal, ra, tag):
    g = np.asarray(cal.mass_gdna_contained, float)  # phantom (truth gDNA=0)
    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])
    print(f"\n  --- {tag} ---   total phantom gDNA(contained) = {g.sum():,.0f}   "
          f"gdna_density_global={cal.gdna_density_global:.4g}")
    print(f"      {'class':>6} {'type':>11} {'nreg':>5} {'phantom_g':>12}")
    for tc in (1, 2, 0):  # intron first (where nascent lives), then exon, intergenic
        for sc in (3, 1, 2, 0):
            m = (tcl == tc) & (scl == sc)
            if m.any() and g[m].sum() > 1.0:
                print(f"      {SC[sc]:>6} {TC[tc]:>11} {int(m.sum()):>5} {g[m].sum():>12,.0f}")


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    kw = dict(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"])
    print(f"=== {cond} : nascent zero-gDNA phantom ablation (all assigned gDNA is phantom) ===")
    report(run(kw, 1, False), ra, "both refit, max_outer=1")
    report(run(kw, 4, False), ra, "both refit, max_outer=4 (CURRENT)")
    report(run(kw, 4, True), ra, "RNA refit only, gDNA FROZEN at seed (max_outer=4)")


if __name__ == "__main__":
    main()
