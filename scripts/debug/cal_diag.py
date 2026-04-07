"""Quick diagnostic: what does the old calibration produce for g0_n30_s65?"""
import logging
logging.basicConfig(level=logging.INFO, format="%(name)s: %(message)s")

from rigel.sim.scenario import Scenario, SimConfig
from rigel.pipeline import run_pipeline
from rigel.config import PipelineConfig, EMConfig, BamScanConfig
import tempfile
import pathlib

tmp = pathlib.Path(tempfile.mkdtemp())
sc = Scenario("debug", genome_length=20000, seed=12345, work_dir=tmp / "debug")
sc.add_gene("g1", "+", [
    {"t_id": "t1", "exons": [(2000, 4000), (8000, 10000)], "abundance": 100},
])
sc.add_gene("g_ctrl", "-", [
    {"t_id": "t_ctrl", "exons": [(14000, 16000), (18000, 19000)], "abundance": 0},
])
result = sc.build_oracle(
    n_fragments=2000,
    sim_config=SimConfig(gdna_abundance=0, nrna_abundance=30, strand_specificity=0.65),
)
pr = run_pipeline(
    result.bam_path, result.index,
    config=PipelineConfig(em=EMConfig(seed=42), scan=BamScanConfig(sj_strand_tag="ts")),
)
cal = pr.calibration
print(f"mixing_proportion: {cal.mixing_proportion:.4f}")
print(f"region_posteriors: {cal.region_posteriors}")
print(f"region_n_total: {cal.region_n_total}")
sc.cleanup()
