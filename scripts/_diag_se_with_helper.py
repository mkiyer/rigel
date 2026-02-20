#!/usr/bin/env python
"""Verify single-exon wipeout WITH trained strand model (helper gene)."""
import sys, tempfile
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.sim.scenario import Scenario
from hulkrna.sim.reads import SimConfig
from hulkrna.pipeline import run_pipeline
from hulkrna.sim.benchmark import run_benchmark

SIM_SEED = 42

def trace_with_helper(ss_value, n_frags=500):
    print(f"\n{'='*72}")
    print(f"SINGLE-EXON + HELPER (trained strand model), ss={ss_value}")
    print(f"{'='*72}\n")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = Scenario("se_helper", genome_length=8000, seed=SIM_SEED,
                       work_dir=tmp / "se_helper")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])
        # Multi-exon helper for strand training
        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(2500, 3000), (3500, 4000)], "abundance": 50},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(6500, 6800)], "abundance": 0},
        ])

        sim_config = SimConfig(
            read_length=100, frag_mean=200, frag_std=30,
            strand_specificity=ss_value, seed=SIM_SEED,
            frag_min=80, frag_max=450,
        )

        result = sc.build(n_fragments=n_frags, sim_config=sim_config)
        pr = run_pipeline(str(result.bam_path), result.index,
                          sj_strand_tag="ts", seed=SIM_SEED)
        bench = run_benchmark(result, pr, scenario_name=f"se_helper_ss{ss_value}")
        counter = pr.counter

        print(f"Strand model: n_obs={pr.strand_models.exonic_spliced.n_observations}, "
              f"SS={pr.strand_models.exonic_spliced.strand_specificity:.4f}")

        sense_g1 = float(counter.gene_sense_all[0])
        anti_g1 = float(counter.gene_antisense_all[0])
        print(f"g1 strand: sense={sense_g1:.0f}, anti={anti_g1:.0f}")

        for ta in bench.transcripts:
            if ta.t_id != "t_helper":
                print(f"  {ta.t_id}: expected={ta.expected:.0f}, "
                      f"observed={ta.observed:.1f}, diff={ta.abs_diff:.1f}")
        print(f"  gDNA: expected=0, pipeline={bench.n_gdna_pipeline:.1f}")
        print()
        sc.cleanup()

for ss in [0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65]:
    trace_with_helper(ss)
