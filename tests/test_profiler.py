"""Focused tests for the profiling helper scripts."""

import json
import importlib.util
import sys
from pathlib import Path

from rigel.sim import Scenario, SimConfig


def _load_profiler_module():
    profiler_path = Path(__file__).resolve().parents[1] / "scripts" / "profiling" / "profiler.py"
    spec = importlib.util.spec_from_file_location("rigel_profile_script", profiler_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_profiler = _load_profiler_module()
ProfileConfig = _profiler.ProfileConfig
RigelParams = _profiler.RigelParams
run_profile = _profiler.run_profile


def test_staged_profile_writes_memory_shape_metrics(tmp_path):
    """A staged profile emits PR05 scan/CSR/partition memory metrics."""
    sc = Scenario(
        "profile_metrics",
        genome_length=5000,
        seed=42,
        work_dir=tmp_path / "profile_metrics",
    )
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
    sc.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50},
        ],
    )
    sim_config = SimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=42,
    )
    try:
        result = sc.build_oracle(n_fragments=200, sim_config=sim_config)
        outdir = tmp_path / "profile_out"
        cfg = ProfileConfig(
            bam=str(result.bam_path),
            index=str(result.index_dir),
            outdir=str(outdir),
            stages=True,
            memory_sample_interval_ms=1000,
            verbose=False,
            rigel_configs=[
                RigelParams(
                    name="metrics",
                    params={
                        "threads": 1,
                        "scan_bgzf_threads": 0,
                        "scan_buffer_size": 0.001,
                    },
                )
            ],
        )
        profiles = run_profile(cfg)
    finally:
        sc.cleanup()

    profile = profiles[0]
    assert profile.scan_config["scan_buffer_size_bytes"] == int(0.001 * 1024**3)
    assert profile.buffer_summary["chunks_finalized"] >= 1
    assert profile.buffer_summary["memory_bytes_peak"] > 0
    assert profile.scoring_csr["n_units"] >= 0
    assert "log_liks" in profile.scoring_csr["candidate_bytes"]
    assert profile.scoring_csr["total_bytes"] >= 0
    assert profile.partition_bytes_total >= 0

    data = json.loads((outdir / "profile_summary.json").read_text())
    json_profile = data["profiles"][0]
    assert "buffer_summary" in json_profile
    assert "scoring_csr" in json_profile
    assert "partition_bytes_total" in json_profile
    assert json_profile["scoring_csr"]["candidate_bytes"]["log_liks"] >= 0