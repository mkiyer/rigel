"""Dev-velocity calibration cache: scan a scenario BAM ONCE, cache the calibration inputs, re-run fast.

Re-scanning a ~500 MB BAM on every calibration iteration dominates dev wall time; calibration itself is
sub-second. This caches the exact inputs ``calibrate()`` consumes (the accumulator payload + the strand model
+ the gDNA/RNA FL pmfs, built precisely as ``pipeline.py`` does) so the inner loop is scan-free. See
``CALIBRATION_PLAN_v6.md`` §13 (the perf plan: prototype-then-optimize, thousands of runs).

    from calib_cache import cached_inputs
    inp = cached_inputs(bam, index_dir, "/tmp/cache.pkl")   # scans once, then loads
    result = calibrate(payload=inp["payload"], region_arrays=inp["region_arrays"],
                       strand_model=inp["strand_model"], gdna_fl_pmf=inp["gdna_fl_pmf"],
                       rna_fl_pmf=inp["rna_fl_pmf"], config=CalibrationConfig())

Run as a script to scan+cache a scenario and print the calibration result (with the Phase-A RNA var~mean
diagnostic at DEBUG):

    python scripts/debug/calib_cache.py <scenario_dir> [index_dir] [--bam sim_oracle.bam] [--debug]
"""

from __future__ import annotations

import dataclasses
import pickle
import sys
from pathlib import Path

from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType


def scan_inputs(bam: str, index_dir: str) -> dict:
    """Scan the BAM once; return the dict of calibration inputs (the exact pipeline.py sourcing)."""
    idx = TranscriptIndex.load(index_dir)
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, strand_model, flm, _buf, payload = scan_and_buffer(bam, idx, scan)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=flm.max_size,
    )
    region_arrays = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    return {
        "payload": payload,
        "region_arrays": region_arrays,
        "strand_model": strand_model,
        "gdna_fl_pmf": fl.gdna_pmf,
        "rna_fl_pmf": fl.rna_pmf,
    }


def cached_inputs(bam: str, index_dir: str, cache_path: str) -> dict:
    """Return cached calibration inputs, scanning + caching on first call."""
    p = Path(cache_path)
    if p.exists():
        with open(p, "rb") as f:
            return pickle.load(f)
    inp = scan_inputs(bam, index_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "wb") as f:
        pickle.dump(inp, f)
    return inp


def _main(argv: list[str]) -> None:
    import logging
    import time

    from rigel.calibration.calibrate import calibrate
    from rigel.config import CalibrationConfig

    debug = "--debug" in argv
    argv = [a for a in argv if a != "--debug"]
    bam_name = "sim_oracle.bam"
    if "--bam" in argv:
        i = argv.index("--bam")
        bam_name = argv[i + 1]
        argv = argv[:i] + argv[i + 2 :]
    scenario = argv[1]
    index_dir = argv[2] if len(argv) > 2 else str(Path(scenario).parent.parent / "refs" / "rigel_index")
    bam = f"{scenario}/{bam_name}"
    cache = f"/tmp/calib_cache_{Path(scenario).name}_{bam_name}.pkl"

    logging.basicConfig(level=logging.DEBUG if debug else logging.INFO, format="%(message)s")
    t0 = time.time()
    inp = cached_inputs(bam, index_dir, cache)
    t_scan = time.time() - t0
    print(f"[cache] inputs ready in {t_scan:.1f}s (cache={cache}, exists-reuse={t_scan < 5})")

    t0 = time.time()
    res = calibrate(
        payload=inp["payload"],
        region_arrays=inp["region_arrays"],
        strand_model=inp["strand_model"],
        gdna_fl_pmf=inp["gdna_fl_pmf"],
        rna_fl_pmf=inp["rna_fl_pmf"],
        config=CalibrationConfig(),
    )
    print(
        f"[calibrate] {time.time() - t0:.2f}s  R={res.n_regions}  "
        f"gdna_density_global={res.gdna_density_global:.4g}  rna_sense_frac={res.rna_sense_frac:.3f}  "
        f"gDNA_mass={res.mass_gdna_contained.sum():.0f}  RNA_mass={res.mass_rna_contained.sum():.0f}"
    )


if __name__ == "__main__":
    _main(sys.argv)
