#!/usr/bin/env python
"""Phase 1 validation — the density model (count clue).

(1) Sparse toy (nrna_dc): ρ₀ from the density model must be BOUNDED (no 5e-5 collapse) and
    Σ gDNA ≈ expected — computed with zero global-ρ₀ feedback.
(2) Overlapping locus (single-exon T1+ over multi-exon T2−, the user's example): the sweep
    must fill the AMBIG interior inward from the only decodable nodes (the locus edges).

Usage:  python scripts/debug/proto_density_model.py
"""

from __future__ import annotations

import tempfile

import numpy as np

from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import coarse_type_array
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark
from rigel.splice import SpliceType

_TS = {0: "NONE", 1: "POS", -1: "NEG", 2: "AMBIG"}
_RT = {0: "intergenic", 1: "intron", 2: "exon"}
_GD = GDNAConfig(abundance=20, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
_CFG = ReadSimConfig(
    frag_mean=200,
    frag_std=30,
    frag_min=80,
    frag_max=450,
    read_length=100,
    strand_specificity=0.65,
    seed=42,
)


def setup(res):
    scan = BamScanConfig(sj_strand_tag="auto")
    _, sm, flm, _, payload = scan_and_buffer(str(res.bam_path), res.index, scan)
    sm.finalize()
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(payload, ra)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=flm.max_size,
    )
    region_eff = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
    mu_fl = boundary_eff_length(fl.gdna_pmf)
    return sub, ra, region_eff, mu_fl


def toy():
    print("\n=== (1) SPARSE TOY (nrna_dc g20 n70 ss65) — no-collapse + recovery ===")
    tmp = tempfile.mkdtemp(prefix="dm_")
    sc = Scenario("nrna_dc", genome_length=20000, seed=42, work_dir=tmp)
    sc.add_gene(
        "g1",
        "+",
        [
            {
                "t_id": "t1",
                "exons": [
                    (2000, 2500),
                    (3000, 3500),
                    (4000, 4500),
                    (5000, 5500),
                    (6000, 6500),
                    (7000, 7500),
                    (8000, 8500),
                    (9500, 10000),
                ],
                "abundance": 100,
            }
        ],
    )
    sc.add_gene(
        "gc",
        "-",
        [
            {
                "t_id": "tc",
                "exons": [
                    (12000, 12500),
                    (13000, 13500),
                    (14000, 14500),
                    (15000, 15500),
                    (16000, 16500),
                    (17000, 17500),
                    (18000, 18500),
                    (19000, 19500),
                ],
                "abundance": 0,
            }
        ],
    )
    res = sc.build_oracle(n_fragments=2000, sim_config=_CFG, gdna_config=_GD, nrna_abundance=70)
    sub, ra, region_eff, mu_fl = setup(res)
    nd = node_gdna_density(sub, ra, region_eff, mu_fl)
    rho0_dm = float(nd.gdna_mass.sum()) / float(region_eff.sum())
    pr = run_pipeline(
        res.bam_path,
        res.index,
        config=PipelineConfig(em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="auto")),
    )
    b = run_benchmark(res, pr, scenario_name="nrna_dc")
    print(
        f"  decodable: {nd.n_region_dec} regions, {nd.n_boundary_dec} boundaries (of {ra.n_regions})"
    )
    print(f"  rho0 (density model, acyclic) = {rho0_dm:.5g}   (true ≈ 0.03)")
    print(f"  rho0 (production EM)          = {pr.calibration.rho_0:.5g}   <- the COLLAPSE")
    print(
        f"  Σ gDNA (density model) = {nd.gdna_mass.sum():.0f}   expected = {b.n_gdna_expected:.0f}"
    )


def overlapping():
    print(
        "\n=== (2) OVERLAPPING LOCUS (T1+ single-exon over T2− multi-exon) — sweep fills interior ==="
    )
    tmp = tempfile.mkdtemp(prefix="ov_")
    sc = Scenario("ov", genome_length=13000, seed=7, work_dir=tmp)
    sc.add_gene("T1", "+", [{"t_id": "t1", "exons": [(1000, 11000)], "abundance": 60}])
    sc.add_gene(
        "T2",
        "-",
        [{"t_id": "t2", "exons": [(2000, 3000), (5000, 6000), (9000, 10000)], "abundance": 40}],
    )
    res = sc.build_oracle(n_fragments=3000, sim_config=_CFG, gdna_config=_GD, nrna_abundance=30)
    sub, ra, region_eff, mu_fl = setup(res)
    nd = node_gdna_density(sub, ra, region_eff, mu_fl)
    rtype = coarse_type_array(ra.signature)
    print(
        f"  decodable: {nd.n_region_dec} regions, {nd.n_boundary_dec} boundaries (of {ra.n_regions})"
    )
    print(
        f"  {'idx':>3} {'ts':>5} {'rtype':>10} {'start':>6} {'end':>6} {'regDec':>6} {'density':>9} {'gDNA':>7}"
    )
    for i in range(ra.n_regions):
        print(
            f"  {i:>3} {_TS.get(int(ra.ts_class[i]), '?'):>5} {_RT.get(int(rtype[i]), '?'):>10} "
            f"{int(ra.start[i]):>6} {int(ra.end[i]):>6} {str(bool(nd.region_decodable[i])):>6} "
            f"{nd.density[i]:>9.4g} {nd.gdna_mass[i]:>7.1f}"
        )
    n_nan = int(np.isnan(nd.density).sum())
    print(
        f"  every region has a density? {'YES' if n_nan == 0 else f'NO ({n_nan} NaN)'} "
        f"(the AMBIG interior should be filled from the locus-edge boundaries)"
    )


def main():
    toy()
    overlapping()


if __name__ == "__main__":
    main()
