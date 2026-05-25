"""Debug fused gDNA evidence and prior allocation for selected oracle scenarios."""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import numpy as np

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario, run_benchmark

SIM_SEED = 42
PIPELINE_SEED = 42


def _sim_config(strand_specificity: float) -> ReadSimConfig:
    return ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=strand_specificity,
        seed=SIM_SEED,
    )


def _gdna_config(abundance: float) -> GDNAConfig | None:
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=100,
        frag_max=1000,
    )


def _build_nrna_double(work_dir: Path) -> Scenario:
    sc = Scenario("debug_nrna_double", genome_length=20000, seed=SIM_SEED, work_dir=work_dir)
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
            },
        ],
    )
    sc.add_gene(
        "g_ctrl",
        "-",
        [
            {
                "t_id": "t_ctrl",
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
            },
        ],
    )
    return sc


def _build_wide_intron(work_dir: Path) -> Scenario:
    sc = Scenario("debug_wide_intron", genome_length=6000, seed=SIM_SEED, work_dir=work_dir)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(1000, 2000), (3000, 4000)], "abundance": 100}])
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(5000, 5300)], "abundance": 0}])
    return sc


def _build_unspliced_single(work_dir: Path) -> Scenario:
    sc = Scenario("debug_unspliced_single", genome_length=5000, seed=SIM_SEED, work_dir=work_dir)
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1500)], "abundance": 100}])
    return sc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("scenario", choices=["nrna_double", "wide_intron", "unspliced_single"])
    parser.add_argument("--n-fragments", type=int, default=2000)
    parser.add_argument("--gdna", type=float, default=0.0)
    parser.add_argument("--nrna", type=float, default=30.0)
    parser.add_argument("--ss", type=float, default=0.9)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="rigel_debug_") as tmp:
        work_dir = Path(tmp) / args.scenario
        if args.scenario == "nrna_double":
            scenario = _build_nrna_double(work_dir)
        elif args.scenario == "wide_intron":
            scenario = _build_wide_intron(work_dir)
        else:
            scenario = _build_unspliced_single(work_dir)
        result = scenario.build_oracle(
            n_fragments=args.n_fragments,
            sim_config=_sim_config(args.ss),
            gdna_config=_gdna_config(args.gdna),
            nrna_abundance=args.nrna,
        )
        pr = run_pipeline(
            result.bam_path,
            result.index,
            config=PipelineConfig(
                em=EMConfig(seed=PIPELINE_SEED),
                scan=BamScanConfig(sj_strand_tag="auto"),
            ),
        )
        bench = run_benchmark(result, pr, scenario_name=args.scenario)
        print(bench.summary())
        print("ground_truth", result.ground_truth_auto())
        print("stats", pr.stats.to_dict())
        print("calibration", pr.calibration.to_summary_dict())
        fused = pr.calibration.fused_region_gdna
        print("fused_mean_sum", float(np.sum(fused.mean_count)))
        print("fused_upper_sum", float(np.sum(fused.upper_count)))
        top = np.argsort(np.asarray(fused.upper_count))[-12:][::-1]
        print("top_fused_regions")
        region_df = result.index.region_df.reset_index(drop=True)
        for idx in top:
            if fused.upper_count[idx] <= 0 and fused.mean_count[idx] <= 0:
                continue
            row = region_df.iloc[int(idx)]
            print(
                int(idx),
                str(row.get("ref_name", "")),
                int(row.start),
                int(row.end),
                "sig",
                int(row.signature),
                "mean",
                float(fused.mean_count[idx]),
                "upper",
                float(fused.upper_count[idx]),
                "obs",
                float(fused.observed_compatible_count[idx]),
                "flags",
                int(fused.flags[idx]),
            )
        prior = pr.calibration.prior_table
        print("prior_summary", prior.to_summary_dict())
        print("prior_expected", prior.gdna_expected_count.tolist())
        print("prior_upper", prior.gdna_upper_count.tolist())
        print("enable_gdna", prior.enable_gdna.tolist())
        print("loci")
        print(pr.estimator.get_loci_df().to_string(index=False))
        print("transcripts")
        print(pr.estimator.get_counts_df(result.index).to_string(index=False))


if __name__ == "__main__":
    main()
