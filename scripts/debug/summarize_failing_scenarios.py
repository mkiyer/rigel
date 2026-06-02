#!/usr/bin/env python
"""Consolidating diagnostic for the 6 PR-8 failing scenario combos.

For each (topology, gdna, nrna, ss) it runs the full pipeline and prints the
calibration health (rho_0, exposure_dispersion, converged, kappa_rna, rho_r_bb)
alongside the headline routing (RNA observed vs expected, gDNA, nRNA). This
separates the two failure classes: a rho_0 runaway (rho_0 huge + not converged)
vs a mis-tuned/knife-edge strand channel.
"""

from __future__ import annotations

import tempfile

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark


def _anti_intron(work):
    sc = Scenario("anti_intron", genome_length=8000, seed=42, work_dir=work)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(1000, 1200), (2000, 2200), (5000, 5600)], "abundance": 100}])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(3000, 3800)], "abundance": 0}])
    sc.add_gene("g_ctrl", "+", [{"t_id": "t_ctrl", "exons": [(7000, 7300)], "abundance": 0}])
    return sc


def _nrna_dc(work):
    sc = Scenario("nrna_dc", genome_length=20000, seed=42, work_dir=work)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(2000, 2500), (3000, 3500), (4000, 4500),
        (5000, 5500), (6000, 6500), (7000, 7500), (8000, 8500), (9500, 10000)], "abundance": 100}])
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(12000, 12500), (13000, 13500),
        (14000, 14500), (15000, 15500), (16000, 16500), (17000, 17500), (18000, 18500),
        (19000, 19500)], "abundance": 0}])
    return sc


def _overlap(work):
    sc = Scenario("overlap_anti", genome_length=8000, seed=42, work_dir=work)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(200, 500), (1000, 1300)], "abundance": 100}])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(300, 600), (1100, 1400)], "abundance": 25}])
    sc.add_gene("g_ctrl", "+", [{"t_id": "t_ctrl", "exons": [(5500, 5800)], "abundance": 0}])
    return sc


# (label, builder, n_frag, gdna, nrna, ss)
CASES = [
    ("anti_intron  nrna30 ss100", _anti_intron, 2000, 0, 30, 1.0),
    ("anti_intron  nrna70 ss100", _anti_intron, 2000, 0, 70, 1.0),
    ("nrna_dc      g20 n0  ss90", _nrna_dc, 2000, 20, 0, 0.9),
    ("nrna_dc      g20 n70 ss65", _nrna_dc, 2000, 20, 70, 0.65),
    ("nrna_dc      g20 n70 ss90", _nrna_dc, 2000, 20, 70, 0.9),
    ("overlap_anti fc4   ss100", _overlap, 1000, 0, 0, 1.0),
]


def main() -> None:
    print(f"{'case':28s} {'rho_0':>10s} {'disp':>6s} {'conv':>4s} {'kappa':>7s} {'rho_r':>7s}"
          f" | {'RNAobs':>7s} {'RNAexp':>7s} {'gDNA':>6s} {'nRNA':>6s}")
    print("-" * 110)
    for label, builder, nfrag, gdna, nrna, ss in CASES:
        work = tempfile.mkdtemp(prefix="sum_")
        sc = builder(work)
        cfg = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                            read_length=100, strand_specificity=ss, seed=42)
        gd = None if gdna == 0 else GDNAConfig(abundance=gdna, frag_mean=350, frag_std=100,
                                               frag_min=100, frag_max=1000)
        res = sc.build_oracle(n_fragments=nfrag, sim_config=cfg, gdna_config=gd, nrna_abundance=nrna)
        pr = run_pipeline(res.bam_path, res.index,
                          config=PipelineConfig(em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="auto")))
        b = run_benchmark(res, pr, scenario_name=label)
        c = pr.calibration
        rna_exp = b.total_expected + b.n_nrna_expected
        print(f"{label:28s} {c.rho_0:10.4g} {c.exposure_dispersion:6.2f} "
              f"{str(c.converged):>4s} {c.kappa_rna:7.4g} {c.rho_r_bb:7.4g} | "
              f"{b.total_rna_observed:7.0f} {rna_exp:7.0f} {b.n_gdna_pipeline:6.0f} {b.n_nrna_pipeline:6.0f}")
        sc.cleanup()


if __name__ == "__main__":
    main()
