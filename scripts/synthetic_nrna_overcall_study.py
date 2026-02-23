#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path

import pandas as pd

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark


def make_scenario(template: str, name: str, seed: int) -> Scenario:
    scenario = Scenario(name, genome_length=9000, seed=seed)
    if template == "single_exon_plus_helper":
        scenario.add_gene("g1", "+", [{"t_id": "t1", "exons": [(600, 1700)], "abundance": 100}])
        scenario.add_gene(
            "g_helper",
            "+",
            [{"t_id": "t_helper", "exons": [(2600, 3100), (3600, 4200)], "abundance": 60}],
        )
        scenario.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(7000, 7400)], "abundance": 0}])
    elif template == "single_spliced_plus_helper":
        scenario.add_gene(
            "g1",
            "+",
            [{"t_id": "t1", "exons": [(600, 1050), (1350, 1750)], "abundance": 100}],
        )
        scenario.add_gene(
            "g_helper",
            "+",
            [{"t_id": "t_helper", "exons": [(2600, 3000), (3300, 3650), (3900, 4300)], "abundance": 40}],
        )
        scenario.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(7000, 7400)], "abundance": 0}])
    else:
        raise ValueError(f"Unknown template: {template}")
    return scenario


def run_study(outdir: Path, n_fragments: int) -> Path:
    if shutil.which("minimap2") is None or shutil.which("samtools") is None:
        raise RuntimeError("minimap2 and/or samtools not available in PATH")

    templates = ("single_exon_plus_helper", "single_spliced_plus_helper")
    nrna_grid = (0, 10, 30, 50)
    gdna_grid = (0.0, 20.0, 50.0)
    strand_grid = (1.0, 0.8)

    rows: list[dict] = []
    run_idx = 0

    for template in templates:
        for nrna_abundance in nrna_grid:
            for gdna_abundance in gdna_grid:
                for strand_specificity in strand_grid:
                    run_idx += 1
                    seed = 1000 + run_idx
                    scenario_name = (
                        f"{template}_nrna{nrna_abundance}_gdna{gdna_abundance}_ss{strand_specificity}"
                    )
                    sim_config = SimConfig(
                        frag_mean=200,
                        frag_std=30,
                        frag_min=80,
                        frag_max=450,
                        read_length=100,
                        strand_specificity=strand_specificity,
                        seed=seed,
                    )
                    gdna_cfg = (
                        None
                        if gdna_abundance <= 0
                        else GDNAConfig(
                            abundance=gdna_abundance,
                            frag_mean=350,
                            frag_std=100,
                            frag_min=100,
                            frag_max=1000,
                        )
                    )

                    with make_scenario(template, scenario_name, seed=seed) as scenario:
                        sim = scenario.build(
                            n_fragments=n_fragments,
                            sim_config=sim_config,
                            gdna_config=gdna_cfg,
                            nrna_abundance=nrna_abundance,
                        )
                        for include_multimap in (False, True):
                            pipe = run_pipeline(
                                sim.bam_path,
                                sim.index,
                                sj_strand_tag="ts",
                                seed=seed,
                                include_multimap=include_multimap,
                            )
                            bench = run_benchmark(sim, pipe, scenario_name=scenario_name)

                            mature_truth = float(bench.total_expected)
                            mature_obs = float(bench.total_observed)
                            nrna_truth = float(bench.n_nrna_expected)
                            nrna_obs = float(bench.n_nrna_pipeline)
                            gdna_truth = float(bench.n_gdna_expected)
                            gdna_obs = float(bench.n_gdna_pipeline)

                            rows.append(
                                {
                                    "template": template,
                                    "scenario_name": scenario_name,
                                    "include_multimap": include_multimap,
                                    "nrna_abundance": nrna_abundance,
                                    "gdna_abundance": gdna_abundance,
                                    "strand_specificity": strand_specificity,
                                    "alignment_rate": float(bench.alignment_rate),
                                    "n_fragments_pipeline": int(bench.n_fragments),
                                    "mature_truth": mature_truth,
                                    "mature_obs": mature_obs,
                                    "mature_signed_error": mature_obs - mature_truth,
                                    "mature_abs_error": abs(mature_obs - mature_truth),
                                    "mature_ratio": (
                                        mature_obs / mature_truth if mature_truth > 0 else math.nan
                                    ),
                                    "nrna_truth": nrna_truth,
                                    "nrna_obs": nrna_obs,
                                    "nrna_signed_error": nrna_obs - nrna_truth,
                                    "nrna_abs_error": abs(nrna_obs - nrna_truth),
                                    "nrna_ratio": nrna_obs / nrna_truth if nrna_truth > 0 else math.nan,
                                    "gdna_truth": gdna_truth,
                                    "gdna_obs": gdna_obs,
                                    "gdna_signed_error": gdna_obs - gdna_truth,
                                    "gdna_abs_error": abs(gdna_obs - gdna_truth),
                                    "gdna_ratio": gdna_obs / gdna_truth if gdna_truth > 0 else math.nan,
                                    "n_intergenic": int(bench.stats_dict.get("n_intergenic", 0)),
                                    "n_with_unannotated_sj": int(
                                        bench.stats_dict.get("n_with_unannotated_sj", 0)
                                    ),
                                    "n_multi_gene": int(bench.stats_dict.get("n_multi_gene", 0)),
                                    "n_multimapper_groups": int(
                                        bench.stats_dict.get("n_multimapper_groups", 0)
                                    ),
                                    "n_multimapper_alignments": int(
                                        bench.stats_dict.get("n_multimapper_alignments", 0)
                                    ),
                                }
                            )
                    print(f"done {scenario_name}")

    df = pd.DataFrame(rows)
    out_csv = outdir / "synthetic_nrna_overcall_study.csv"
    df.to_csv(out_csv, index=False)
    return out_csv


def main() -> int:
    parser = argparse.ArgumentParser(description="Synthetic nRNA overcall study")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--n-fragments", type=int, default=300)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_csv = run_study(args.outdir, args.n_fragments)
    print(f"WROTE {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
