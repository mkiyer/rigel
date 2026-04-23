#!/usr/bin/env python
"""Diagnostic: dump calibration region_counts + fl_table + summary pieces for one BAM.

Writes to <outdir>/{region_counts.feather, fl_table.feather}.

Purpose: enable offline re-estimation of λ_G from arbitrary region
subsets (intergenic-only, deep-intronic, etc.) without re-scanning.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from rigel.config import BamScanConfig, PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    index = TranscriptIndex.load(Path(args.index))
    scan_cfg = BamScanConfig(n_scan_threads=args.threads)
    _, strand_models, fl_models, _, region_counts, fl_table = scan_and_buffer(
        str(args.bam), index, scan_cfg,
    )
    if region_counts is None or fl_table is None:
        raise SystemExit("No region evidence collected")

    region_counts.to_feather(outdir / "region_counts.feather")
    fl_table.to_feather(outdir / "fl_table.feather")
    meta = {
        "mean_frag_len": float(fl_models.global_model.mean),
        "n_fl_obs": int(len(fl_table)),
    }
    import json
    (outdir / "scan_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"wrote {outdir}/region_counts.feather ({len(region_counts):,} rows)")
    print(f"wrote {outdir}/fl_table.feather ({len(fl_table):,} rows)")
    print(f"meta: {meta}")


if __name__ == "__main__":
    main()
