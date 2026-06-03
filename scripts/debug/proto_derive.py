#!/usr/bin/env python
"""Phase 4 validation — derive ρ₀ + per-node exposure ω (the pure local/global ratio).

  (1) nrna_dc: ρ₀ bounded/sane; ω finite; report min/max/mean.
  (2) zero-gDNA: ρ₀→0, ω→1 (neutral) everywhere.
  (3) zero-RNA: ρ₀ high; ω≈1 (uniform gDNA).
  (4) capture: ω>1 at probe-target regions, <1 off-target (exposure recovers enrichment).

Usage:  python scripts/debug/proto_derive.py
"""

from __future__ import annotations

import tempfile
from pathlib import Path

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
from rigel.sim import CaptureConfig, GDNAConfig, ReadSimConfig, Scenario
from rigel.splice import SpliceType

_CFG = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                     read_length=100, strand_specificity=0.65, seed=42)
_GD = GDNAConfig(abundance=20, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
_T1 = [(2000, 2500), (3000, 3500), (4000, 4500), (5000, 5500),
       (6000, 6500), (7000, 7500), (8000, 8500), (9500, 10000)]
_TC = [(12000, 12500), (13000, 13500), (14000, 14500), (15000, 15500),
       (16000, 16500), (17000, 17500), (18000, 18500), (19000, 19500)]
_EX8 = [(2000, 2500), (3000, 3500), (4000, 4500), (5000, 5500),
        (6000, 6500), (7000, 7500), (8000, 8500), (9500, 10000)]


def setup_and_derive(res):
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
    return ra, derive(reg, left, right, region_eff, bside, ra.ref_id)


def _stats(name, de, mask=None):
    w = de.region_omega if mask is None else de.region_omega[mask]
    finite = np.all(np.isfinite(de.region_omega)) and np.all(np.isfinite(de.left_omega)) and np.all(np.isfinite(de.right_omega))
    print(f"  {name}: rho_0={de.rho_0:.5g}  finite={finite}  region ω: "
          f"min={w.min():.3f} mean={w.mean():.3f} max={w.max():.3f}")


def nrna_dc(tmp, *, gdna, nrna, t1_ab=100, tc_ab=0):
    sc = Scenario("nrna_dc", genome_length=20000, seed=42, work_dir=tmp)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": _T1, "abundance": t1_ab}])
    sc.add_gene("gc", "-", [{"t_id": "tc", "exons": _TC, "abundance": tc_ab}])
    return sc.build_oracle(n_fragments=2000, sim_config=_CFG,
                           gdna_config=(None if gdna == 0 else _GD), nrna_abundance=nrna)


def main():
    print("\n=== Phase 4: derive ρ₀ + exposure ω ===")
    _, de = setup_and_derive(nrna_dc(tempfile.mkdtemp(prefix="d1_"), gdna=20, nrna=70))
    _stats("(1) nrna_dc (g20 n70)   ", de)
    _, de = setup_and_derive(nrna_dc(tempfile.mkdtemp(prefix="d2_"), gdna=0, nrna=70))
    _stats("(2) zero-gDNA           ", de)
    _, de = setup_and_derive(nrna_dc(tempfile.mkdtemp(prefix="d3_"), gdna=20, nrna=0, t1_ab=0, tc_ab=0))
    _stats("(3) zero-RNA (silent)   ", de)

    # (4) capture: 6 genes, probe-capture the first 3 → focal gDNA enrichment there.
    tmp = tempfile.mkdtemp(prefix="d4_")
    sc = Scenario("cap", genome_length=120000, seed=42, work_dir=tmp)
    captured = {0, 1, 2}
    for i in range(6):
        o = i * 20000
        sc.add_gene(f"g{i}", "+", [{"t_id": f"t{i}", "exons": [(o + a, o + b) for a, b in _EX8], "abundance": 50}])
    probe = Path(tmp) / "probes.tsv"
    with open(probe, "w") as fh:
        fh.write("transcript_id\tstart\tend\n")
        for i in sorted(captured):
            fh.write(f"t{i}\t0\t4000\n")
    cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                        read_length=100, strand_specificity=0.9, seed=42)
    res = sc.build_oracle(n_fragments=12000, sim_config=cfg,
                          gdna_config=GDNAConfig(abundance=40, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000),
                          nrna_abundance=0,
                          capture_config=CaptureConfig(probes=str(probe), binding_per_base=30.0))
    ra, de = setup_and_derive(res)
    cap_mask = np.array([(int(s) // 20000) in captured for s in ra.start])
    print(f"\n  (4) capture: rho_0={de.rho_0:.4g}")
    _stats("      captured regions  ", de, cap_mask)
    _stats("      off-target regions", de, ~cap_mask)
    print("  → captured ω should be ≫1 (enrichment), off-target ω <1 (depletion).")


if __name__ == "__main__":
    main()
