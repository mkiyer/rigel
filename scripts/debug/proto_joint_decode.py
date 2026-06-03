#!/usr/bin/env python
"""Phase 3 validation — the joint decode (count prior × Beta-Binomial strand likelihood).

Checks:
  (1) nrna_dc toy: total RNA / gDNA recovery (regions + boundaries), mass conservation,
      no NaNs. gDNA expected somewhat high (honest nascent bias, plan §0.1).
  (2) zero-gDNA scenario: π_g → 0, gDNA ≈ 0, RNA ≈ total.
  (3) zero-RNA scenario: π_g → 1, gDNA ≈ total.
  (4) stability: no NaN/inf anywhere.

Usage:  python scripts/debug/proto_joint_decode.py
"""

from __future__ import annotations

import tempfile

import numpy as np

from rigel.calibration.density_model import node_gdna_density
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
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark
from rigel.splice import SpliceType

_CFG = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                     read_length=100, strand_specificity=0.65, seed=42)
_GD = GDNAConfig(abundance=20, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
_T1 = [(2000, 2500), (3000, 3500), (4000, 4500), (5000, 5500),
       (6000, 6500), (7000, 7500), (8000, 8500), (9500, 10000)]
_TC = [(12000, 12500), (13000, 13500), (14000, 14500), (15000, 15500),
       (16000, 16500), (17000, 17500), (18000, 18500), (19000, 19500)]


def nrna_dc(tmp, *, gdna, nrna, t1_ab=100, tc_ab=0):
    sc = Scenario("nrna_dc", genome_length=20000, seed=42, work_dir=tmp)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": _T1, "abundance": t1_ab}])
    sc.add_gene("gc", "-", [{"t_id": "tc", "exons": _TC, "abundance": tc_ab}])
    gd = None if gdna == 0 else _GD
    return sc.build_oracle(n_fragments=2000, sim_config=_CFG, gdna_config=gd, nrna_abundance=nrna)


def setup(res):
    _, sm, flm, _, payload = scan_and_buffer(str(res.bam_path), res.index, BamScanConfig(sj_strand_tag="auto"))
    sm.finalize()
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(payload, ra)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    region_eff = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    mu_fl = boundary_eff_length(fl.gdna_pmf)
    bside = boundary_side_eff_length(fl.gdna_pmf, ra.region_size_bp)  # per-side density length
    nd = node_gdna_density(sub, ra, region_eff, mu_fl)
    kappa = fit_strand_balance(sm).kappa_rna
    return sub, ra, nd, region_eff, bside, kappa


def decode(sub, ra, nd, region_eff, bside, kappa):
    reg = decode_regions(sub, ra, nd, region_eff, kappa_rna=kappa)
    left, right = decode_sides(sub, ra, nd, bside, kappa_rna=kappa)
    g = float(reg.gdna_mass.sum() + left.gdna_mass.sum() + right.gdna_mass.sum())
    r = float(reg.rna_mass.sum() + left.rna_mass.sum() + right.rna_mass.sum())
    finite = all(np.all(np.isfinite(d.gdna_mass)) for d in (reg, left, right))
    print(f"    [regions g={reg.gdna_mass.sum():6.0f} r={reg.rna_mass.sum():6.0f} | "
          f"sides(L+R) g={left.gdna_mass.sum() + right.gdna_mass.sum():6.0f} "
          f"r={left.rna_mass.sum() + right.rna_mass.sum():6.0f}]")
    return g, r, reg, (left, right), finite


def main():
    print("\n=== (1) nrna_dc toy (g20 n70) — recovery + mass conservation ===")
    res = nrna_dc(tempfile.mkdtemp(prefix="jd1_"), gdna=20, nrna=70)
    sub, ra, nd, region_eff, mu_fl, kappa = setup(res)
    g, r, reg, bnd, finite = decode(sub, ra, nd, region_eff, mu_fl, kappa)
    pr = run_pipeline(res.bam_path, res.index, config=PipelineConfig(em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="auto")))
    b = run_benchmark(res, pr, scenario_name="nrna_dc")
    rna_exp = b.total_expected + b.n_nrna_expected
    total_decoded = g + r
    total_mass = float((sub.contained.mass_unspliced + sub.contained.mass_spliced).sum()
                       + sum((s.mass_unspliced + s.mass_spliced).sum() for s in (sub.left, sub.right)))
    print(f"  kappa_rna={kappa:.3f}  finite={finite}")
    print(f"  gDNA decoded = {g:7.0f}   expected = {b.n_gdna_expected:.0f}   (high = honest nascent bias)")
    print(f"  RNA  decoded = {r:7.0f}   expected = {rna_exp:.0f}   ({r / rna_exp:.0%})")
    print(f"  mass conservation: decoded {total_decoded:.0f} vs accumulator {total_mass:.0f} "
          f"(contained + left-side + right-side, each node once)")

    print("\n=== (2) ZERO gDNA (gene + nascent RNA only) — expect π_g→0 ===")
    res = nrna_dc(tempfile.mkdtemp(prefix="jd2_"), gdna=0, nrna=70)
    g, r, reg, bnd, finite = decode(*setup(res))
    print(f"  gDNA decoded = {g:.1f}   RNA decoded = {r:.0f}   finite={finite}   "
          f"mean region π_g = {reg.pi_g[reg.gdna_mass + reg.rna_mass > 0].mean():.3f} (→0 good)")

    print("\n=== (3) ZERO RNA (gDNA only, all genes silent) — expect π_g→1 ===")
    res = nrna_dc(tempfile.mkdtemp(prefix="jd3_"), gdna=20, nrna=0, t1_ab=0, tc_ab=0)
    sub3, ra3, nd3, re3, mu3, ka3 = setup(res)
    g, r, reg, bnd, finite = decode(sub3, ra3, nd3, re3, mu3, ka3)
    occ = reg.gdna_mass + reg.rna_mass > 0
    print(f"  gDNA decoded = {g:.0f}   RNA decoded = {r:.1f}   finite={finite}   "
          f"mean region π_g = {reg.pi_g[occ].mean():.3f} (→1 good)")


if __name__ == "__main__":
    main()
