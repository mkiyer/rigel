#!/usr/bin/env python
"""Phase 5a — ρ₀ / exposure stability across seeds and depth (the redesign's vindication).

The old calibration collapsed in this exact regime (nrna_dc g20, sparse nascent): the
``ρ₀ → decode → ρ₀`` EM loop drove RNA recovery on a 20%→100% swing as ρ₀ ran away, and
it "dissolved at scale". The acyclic redesign has NO such loop — derive() is a single
feed-forward pass — so it cannot collapse by construction. This harness shows the
single-pass estimate is ALSO empirically stable: low coefficient-of-variation (CV) of ρ₀
and decoded-RNA% across seeds, holding flat as depth scales (not climbing toward 100%).

The ABSOLUTE RNA% sits below truth by the accepted nascent-confound bias (plan §0.1);
the claim here is STABILITY (low CV, flat vs depth), not accuracy.

Usage:  python scripts/debug/proto_stability.py
"""

from __future__ import annotations

import tempfile

import numpy as np

from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.derive import derive
from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.joint_decode import decode_regions, decode_sides
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario
from rigel.splice import SpliceType

_GD = GDNAConfig(abundance=20, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
_T1 = [(2000, 2500), (3000, 3500), (4000, 4500), (5000, 5500),
       (6000, 6500), (7000, 7500), (8000, 8500), (9500, 10000)]
_TC = [(12000, 12500), (13000, 13500), (14000, 14500), (15000, 15500),
       (16000, 16500), (17000, 17500), (18000, 18500), (19000, 19500)]


def run_one(seed, depth, *, nrna=70):
    tmp = tempfile.mkdtemp(prefix="s5_")
    sc = Scenario("nrna_dc", genome_length=20000, seed=seed, work_dir=tmp)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": _T1, "abundance": 100}])
    sc.add_gene("gc", "-", [{"t_id": "tc", "exons": _TC, "abundance": 0}])
    cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                        read_length=100, strand_specificity=0.65, seed=seed)
    res = sc.build_oracle(n_fragments=depth, sim_config=cfg, gdna_config=_GD, nrna_abundance=nrna)

    _, sm, flm, _, payload = scan_and_buffer(str(res.bam_path), res.index, BamScanConfig(sj_strand_tag="auto"))
    sm.finalize()
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(payload, ra)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    region_eff = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    bside = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)
    nd = node_gdna_density(sub, ra, region_eff, boundary_eff_length(fl.gdna_pmf))
    kappa = fit_strand_balance(sm).kappa_rna
    reg = decode_regions(sub, ra, nd, region_eff, kappa_rna=kappa)
    left, right = decode_sides(sub, ra, nd, bside, kappa_rna=kappa)
    de = derive(reg, left, right, region_eff, bside, ra.ref_id)
    g = float(reg.gdna_mass.sum() + left.gdna_mass.sum() + right.gdna_mass.sum())
    r = float(reg.rna_mass.sum() + left.rna_mass.sum() + right.rna_mass.sum())
    return de.rho_0, r / (g + r) if (g + r) > 0 else 0.0


def _cv(x):
    x = np.asarray(x)
    return float(x.std() / x.mean()) if x.mean() != 0 else 0.0


def main():
    seeds = [1, 7, 42, 123, 2024]
    print("\n=== Phase 5a: stability across seeds & depth (nrna_dc g20 n70) ===")
    print("  (old calibration: RNA% swung 20%→100% here and dissolved at scale)\n")
    print(f"  {'depth':>7} | {'rho_0 mean':>11} {'CV':>6} | {'RNA% mean':>10} {'CV':>6}")
    print("  " + "-" * 52)
    for depth in [2000, 8000, 20000]:
        rhos, rnas = [], []
        for s in seeds:
            rho, rf = run_one(s, depth)
            rhos.append(rho)
            rnas.append(rf)
        print(f"  {depth:>7} | {np.mean(rhos):>11.5f} {_cv(rhos):>6.1%} | "
              f"{np.mean(rnas):>9.1%} {_cv(rnas):>6.1%}")
    print("\n  → stable iff CV is small AND RNA% is flat across depth (no climb to 100%).")


if __name__ == "__main__":
    main()
