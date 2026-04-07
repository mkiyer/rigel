"""Trace ESS normalization during pipeline execution."""
import logging
import sys
import tempfile
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.WARNING)

from rigel.frag_length_model import FragmentLengthModel

_orig = FragmentLengthModel.finalize


def _trace(self, prior_counts=None, prior_ess=None):
    if prior_counts is not None:
        n = self.max_size + 1
        pc = np.zeros(n, dtype=np.float64)
        m = min(n, len(prior_counts))
        pc[:m] = prior_counts[:m]
        raw = float(pc.sum())
        if prior_ess is not None and raw > 0:
            eff = min(prior_ess, raw)
        else:
            eff = raw
        print(
            f"TRACE: cat_weight={self._total_weight:.0f} "
            f"raw_prior={raw:.0f} ess_param={prior_ess} "
            f"eff_ess={eff:.0f}",
            flush=True,
        )
    else:
        print(f"TRACE: cat_weight={self._total_weight:.0f} no prior", flush=True)
    _orig(self, prior_counts=prior_counts, prior_ess=prior_ess)


FragmentLengthModel.finalize = _trace

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

tmpdir = Path(tempfile.mkdtemp())
sc = Scenario("t", genome_length=20000, seed=42, work_dir=tmpdir / "t")
sc.add_gene(
    "g1",
    "+",
    [
        {
            "t_id": "t1",
            "exons": [(2000, 4000), (8000, 10000)],
            "abundance": 100,
        },
    ],
)
sc.add_gene(
    "g_ctrl",
    "-",
    [
        {
            "t_id": "t_ctrl",
            "exons": [(14000, 16000), (18000, 19000)],
            "abundance": 0,
        },
    ],
)

sim_cfg = SimConfig(
    strand_specificity=0.65,
    seed=42,
    read_length=75,
    frag_mean=200,
    frag_std=30,
)
result = sc.build_oracle(
    n_fragments=2000,
    sim_config=sim_cfg,
    gdna_config=GDNAConfig(),
    nrna_abundance=70,
)
config = PipelineConfig(
    em=EMConfig(seed=12345),
    scan=BamScanConfig(sj_strand_tag="auto"),
)
print(f"fl_prior_ess={config.calibration.fl_prior_ess}", flush=True)
pr = run_pipeline(result.bam_path, result.index, config=config)
print("DONE", flush=True)
