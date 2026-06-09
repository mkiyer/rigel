#!/usr/bin/env python3
"""TEMP DEBUG — build a stress scenario with a RUN of consecutive non-observable regions.

An antisense overlap (+ exon over a − exon) creates the region run  ex+ | ex+|ex- | ex-  whose
MIDDLE region (the AMBIG ex+|ex- overlap) has BOTH boundaries non-observable (each shares an exon
bit with a neighbour) — so it is reachable only by sweeping inward from the run's two observable
edges. A multi-exon "seed" gene supplies spliced reads so the calibrator can identify the strand.

Writes an oracle scenario to /tmp/rigel_sim/antisense_run; inspect with the walk/autopsy tools:

  python scripts/debug/gene37_sweep_walk.py --index /tmp/rigel_sim/antisense_run/index \
    --bam /tmp/rigel_sim/antisense_run/antisense_run_oracle.bam \
    --ref antisense_run --start 4500 --end 8500
"""
from __future__ import annotations

from rigel.sim import Scenario
from rigel.sim.reads import GDNAConfig, ReadSimConfig


def main():
    sc = Scenario(
        "antisense_run", genome_length=20000, seed=7,
        work_dir="/tmp/rigel_sim/antisense_run",
    )
    # multi-exon seed gene → spliced reads (strand identifiability)
    sc.add_gene("seed", "+", [{"t_id": "seed.1",
                               "exons": [(500, 800), (1200, 1500), (2000, 2300)],
                               "abundance": 60}])
    # antisense overlap: + exon [5000,7000) over − exon [6000,8000)
    sc.add_gene("A", "+", [{"t_id": "A.1", "exons": [(5000, 7000)], "abundance": 80}])
    sc.add_gene("B", "-", [{"t_id": "B.1", "exons": [(6000, 8000)], "abundance": 80}])
    res = sc.build_oracle(
        n_fragments=300000,
        sim_config=ReadSimConfig(strand_specificity=0.99, frag_mean=250, frag_std=50, seed=7),
        gdna_config=GDNAConfig(abundance=4000.0, frag_mean=350, frag_std=100),
    )
    print("index:", res.index_dir if hasattr(res, "index_dir") else res.index)
    print("bam:  ", res.bam_path)
    print("gtf:  ", res.gtf_path)


if __name__ == "__main__":
    main()
