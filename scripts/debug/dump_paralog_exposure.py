"""Dump per-region exposure detail for the paralog scenario at gdna=20."""
from __future__ import annotations

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "tests" / "scenarios_aligned"))
from conftest import (  # noqa: E402
    GDNAConfig,
    PIPELINE_SEED,
    ReadSimConfig,
    SIM_SEED,
    run_pipeline,
)

from rigel.config import BamScanConfig, EMConfig, PipelineConfig  # noqa: E402
from rigel.sim import Scenario  # noqa: E402

logging.basicConfig(level=logging.WARNING)
work_dir = Path("/tmp/paralog_diag_dump")
if work_dir.exists():
    import shutil
    shutil.rmtree(work_dir)
work_dir.mkdir()

sc = Scenario("paralog_dump", genome_length=12000, seed=SIM_SEED, work_dir=work_dir)
sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1000)], "abundance": 100}])
sc.add_gene("g2", "+", [{"t_id": "t2", "exons": [(5000, 5500)], "abundance": 100}])
sc.genome.edit(5000, sc.genome[500:1000])
sc.add_gene("g_helper", "+", [{"t_id": "t_helper", "exons": [(8000, 8300), (8700, 9000)], "abundance": 50}])
sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0}])

sim = ReadSimConfig(frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
                   read_length=100, strand_specificity=1.0, seed=SIM_SEED)
gdna_cfg = GDNAConfig(abundance=20, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)
res = sc.build(n_fragments=3000, sim_config=sim, gdna_config=gdna_cfg, nrna_abundance=0)
cfg = PipelineConfig(em=EMConfig(seed=PIPELINE_SEED),
                    scan=BamScanConfig(sj_strand_tag="ts", include_multimap=True))
pr = run_pipeline(res.bam_path, res.index, config=cfg)

cal = pr.calibration
rc = cal.region_calibration
exp = rc.region_exposure
mass = rc.region_unspliced_mass

import numpy as np
print(f"tau2={exp.tau2:.4g}  tau2_hat={exp.tau2_hat:.4g}  method={exp.tau2_method}  pool_size={exp.tau2_pool_size}  rho0={exp.rho0:.4g}")
print(f"{'idx':>4} {'bp':>6} {'gdna':>8} {'rna':>8} {'p_unx':>6} {'raw_r':>8} {'log_raw':>8} {'v_obs':>9} {'shrink':>8} {'omega':>8}")
p_unx_arr = np.asarray(rc.p_unexpressed, dtype=np.float64)
for i in range(len(exp.omega)):
    print(f"{i:>4d} {float(mass.region_size_bp[i]):>6.0f} {float(mass.gdna_mass[i]):>8.2f} {float(mass.rna_mass[i]):>8.2f} {float(p_unx_arr[i]):>6.3f} {float(exp.raw_ratio[i]):>8.3g} {float(exp.log_raw_ratio[i]):>8.3f} {float(exp.v_obs[i]):>9.3g} {float(exp.shrink_weight[i]):>8.3g} {float(exp.omega[i]):>8.3g}")
