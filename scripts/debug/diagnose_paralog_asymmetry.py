"""Diagnose Phase 5 paralog t1/t2 asymmetry regression.

Reproduces tests/scenarios_aligned/test_paralogs.py::test_gdna_sweep[gdna_20]
and [gdna_100], dumps calibration + locus + EM diagnostics, then re-runs with
multiple pipeline seeds to test whether asymmetry is stochastic or systematic.
"""
from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "tests" / "scenarios_aligned"))
from conftest import (  # noqa: E402  - test-fixture import
    GDNAConfig,
    N_FRAGMENTS,
    PIPELINE_SEED,
    ReadSimConfig,
    SIM_SEED,
    run_pipeline,
)

from rigel.config import BamScanConfig, EMConfig, PipelineConfig  # noqa: E402
from rigel.sim import Scenario, run_benchmark  # noqa: E402

logging.basicConfig(level=logging.WARNING, format="%(message)s")
log = logging.getLogger("paralog-diag")
log.setLevel(logging.INFO)


def make_scenario(work_dir: Path) -> Scenario:
    sc = Scenario(
        "paralogs_diag",
        genome_length=12000,
        seed=SIM_SEED,
        work_dir=work_dir,
    )
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(500, 1000)], "abundance": 100}])
    sc.add_gene("g2", "+", [{"t_id": "t2", "exons": [(5000, 5500)], "abundance": 100}])
    sc.genome.edit(5000, sc.genome[500:1000])
    sc.add_gene("g_helper", "+", [
        {"t_id": "t_helper", "exons": [(8000, 8300), (8700, 9000)], "abundance": 50},
    ])
    sc.add_gene("g_ctrl", "-", [{"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0}])
    return sc


def run_one(work_dir: Path, *, gdna: int, pipeline_seed: int = PIPELINE_SEED) -> dict:
    sc = make_scenario(work_dir)
    sim = ReadSimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=1.0, seed=SIM_SEED,
    )
    gdna_cfg = GDNAConfig(
        abundance=gdna, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000,
    ) if gdna > 0 else None
    res = sc.build(n_fragments=3000, sim_config=sim, gdna_config=gdna_cfg, nrna_abundance=0)
    cfg = PipelineConfig(
        em=EMConfig(seed=pipeline_seed),
        scan=BamScanConfig(sj_strand_tag="ts", include_multimap=True),
    )
    pr = run_pipeline(res.bam_path, res.index, config=cfg)
    bench = run_benchmark(res, pr, scenario_name=f"paralog_g{gdna}_s{pipeline_seed}")
    by_id = {t.t_id: t for t in bench.transcripts}
    out = {
        "gdna": gdna,
        "seed": pipeline_seed,
        "t1": float(by_id["t1"].observed),
        "t2": float(by_id["t2"].observed),
        "t_helper": float(by_id["t_helper"].observed),
        "t_ctrl": float(by_id["t_ctrl"].observed),
    }
    # Calibration: rho0, omega, etc.
    cal = getattr(pr, "calibration", None) or getattr(pr, "calibration_result", None)
    if cal is not None:
        try:
            summ = cal.to_summary_dict()
            out["calibration"] = {
                "rho0_mean": summ.get("region_calibration", {})
                    .get("background_density", {}).get("rho0_mean"),
                "n_regions_in_pool": summ.get("region_calibration", {})
                    .get("background_density", {}).get("n_regions_in_pool"),
                "info_histogram": summ.get("region_calibration", {})
                    .get("background_density", {}).get("info_histogram")
                    or summ.get("region_calibration", {})
                    .get("background_density", {}).get("method_histogram"),
                "fit_status": summ.get("region_calibration", {})
                    .get("background_density", {}).get("fit_status"),
                "tau2_method": summ.get("region_calibration", {})
                    .get("region_exposure", {}).get("tau2_method"),
            }
        except Exception as exc:  # pragma: no cover - diagnostic
            out["calibration_error"] = repr(exc)
    return out


def main() -> None:
    work_root = Path("/tmp/paralog_diag")
    work_root.mkdir(exist_ok=True)
    rows = []
    for gdna in (0, 20, 100):
        for seed in (42, 1, 7, 123, 999):
            run_dir = work_root / f"g{gdna}_s{seed}"
            if run_dir.exists():
                import shutil as _sh
                _sh.rmtree(run_dir)
            run_dir.mkdir()
            try:
                r = run_one(run_dir, gdna=gdna, pipeline_seed=seed)
                t1, t2 = r["t1"], r["t2"]
                tot = t1 + t2
                ratio = t1 / t2 if t2 > 0 else float("inf")
                asym = abs(t1 - t2) / tot if tot > 0 else 0.0
                cal = r.get("calibration", {})
                log.info(
                    "gdna=%3d seed=%3d  t1=%6.1f t2=%6.1f  ratio=%.3f  asym=%.3f"
                    "  rho0=%.4g  pool=%s  fit=%s  tau2=%s",
                    gdna, seed, t1, t2, ratio, asym,
                    cal.get("rho0_mean", 0.0),
                    cal.get("n_regions_in_pool"),
                    cal.get("fit_status"),
                    cal.get("tau2_method"),
                )
                rows.append({
                    "gdna": gdna, "seed": seed, "t1": t1, "t2": t2,
                    "ratio": ratio, "asym": asym,
                    "rho0": cal.get("rho0_mean"),
                    "pool": cal.get("n_regions_in_pool"),
                    "fit": cal.get("fit_status"),
                    "tau2": cal.get("tau2_method"),
                })
            except Exception as exc:
                log.exception("gdna=%d seed=%d FAILED: %s", gdna, seed, exc)
    df = pd.DataFrame(rows)
    out_csv = Path("/tmp/paralog_diag_results.csv")
    df.to_csv(out_csv, index=False)
    log.info("Wrote %s", out_csv)
    log.info("\nSummary (mean asym per gdna level):")
    log.info("\n%s", df.groupby("gdna")[["asym", "ratio"]].describe())


if __name__ == "__main__":
    main()
