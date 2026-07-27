#!/usr/bin/env python3
"""Build an antisense-overlap STRESS scenario for the AMBIG (opposite-strand) calibration problem.

The standard suite has only single-strand genes, so it never exercises regions where + and −
transcripts overlap (gDNA + sense-RNA + antisense-RNA) — the genuinely-unidentifiable-by-strand
case common in the human transcriptome. This builds a richer oracle scenario:

  * THREE multi-exon single-strand genes (+,−,+) — spliced reads to train the strand model, and
    intron/intergenic regions so the global density + overdispersions calibrate well;
  * TWO opposite-strand overlap pairs — a partial overlap and a nested (contained) overlap — each
    creating an AMBIG region (ex+ | ex−), with both single-strand flanks count-observable;
  * gDNA contamination (and optional hybrid capture) on top.

Parametrized by gDNA abundance, strand specificity, and capture so the AMBIG handling can be swept
across the SS × capture spectrum (the regime where the current categorical count-vs-strand fix is
expected to floor).

  python scripts/debug/build_antisense_stress_scenario.py            # gDNA on, ss0.99, no capture
"""
from __future__ import annotations

import argparse

from rigel.sim import Scenario
from rigel.sim.reads import GDNAConfig, ReadSimConfig

GENOME = 80000


def build(gdna_abundance: float, strand_specificity: float, work_dir: str, seed: int = 11):
    sc = Scenario("antisense_stress", genome_length=GENOME, seed=seed, work_dir=work_dir)

    # --- single-strand multi-exon genes: strand training + observable deconvolvable regions ---
    sc.add_gene("S1", "+", [{"t_id": "S1.1",
                             "exons": [(5000, 5800), (7000, 7800), (9000, 9800), (11000, 12000)],
                             "abundance": 60}])
    sc.add_gene("S2", "-", [{"t_id": "S2.1",
                             "exons": [(20000, 21000), (23000, 23800), (26000, 27000)],
                             "abundance": 50}])
    sc.add_gene("S3", "+", [{"t_id": "S3.1",
                             "exons": [(35000, 36000), (38000, 38800), (41000, 42000)],
                             "abundance": 45}])

    # --- antisense overlap pair A: partial overlap [52000,54000) = AMBIG ---
    sc.add_gene("A_plus", "+", [{"t_id": "A_plus.1", "exons": [(50000, 54000)], "abundance": 70}])
    sc.add_gene("A_minus", "-", [{"t_id": "A_minus.1", "exons": [(52000, 56000)], "abundance": 70}])

    # --- antisense overlap pair B: nested overlap [62000,64000) = AMBIG, flanked by + exon ---
    sc.add_gene("B_plus", "+", [{"t_id": "B_plus.1", "exons": [(60000, 66000)], "abundance": 65}])
    sc.add_gene("B_minus", "-", [{"t_id": "B_minus.1", "exons": [(62000, 64000)], "abundance": 65}])

    res = sc.build_oracle(
        n_fragments=400000,
        sim_config=ReadSimConfig(
            strand_specificity=strand_specificity, frag_mean=250, frag_std=50, seed=seed
        ),
        gdna_config=GDNAConfig(abundance=gdna_abundance, frag_mean=350, frag_std=100),
    )
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gdna", type=float, default=2000.0)
    ap.add_argument("--ss", type=float, default=0.99)
    ap.add_argument("--work-dir", default="/tmp/rigel_sim/antisense_stress")
    args = ap.parse_args()
    res = build(args.gdna, args.ss, args.work_dir)
    print("index:", res.index_dir)
    print("bam:  ", res.bam_path)
    print("gtf:  ", res.gtf_path)
    print("AMBIG overlap regions: A [52000,54000)  B [62000,64000)")


if __name__ == "__main__":
    main()
