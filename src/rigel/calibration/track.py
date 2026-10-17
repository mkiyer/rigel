"""Genome-wide gDNA track from the calibration result.

Calibration solves each region for a gDNA/RNA split; joined to the region
coordinates this yields a per-region gDNA level across the genome — a QC track
that traces contamination and, on capture libraries, on-target enrichment.

This is pure persistence: every value is already in ``CalibrationResult`` (the
per-region deconvolved counts + the density-correct effective length) and
``RegionArrays`` (the coordinates). No new calibration computation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import CONSTANTS


def build_gdna_track(calibration, region_arrays, ref_names) -> pd.DataFrame:
    """Per-region gDNA track: ``(ref, start, end, count_gdna, count_rna, gdna_density, fraction_gdna)``.

    Rows are in genomic order (``RegionArrays`` is sorted by ``(ref_id, start)``),
    which is the same order the ``CalibrationResult`` per-REGION arrays are aligned
    to. ``count_gdna`` and ``count_rna`` are the region's contained fragments
    calibration deconvolves into gDNA and RNA; ``gdna_density`` is the density-correct
    contained gDNA count per bp (``count_gdna`` / ``gdna_region_eff_len``);
    ``fraction_gdna`` is gDNA's share of the region's contained count.

    Contained only, deliberately: a contiguous boundary is a 0-bp boundary, so it has no genomic
    extent to occupy a row of a track, and attributing its crossing count to a flank region would
    report a density at a position where those fragments were never contained.
    """
    gdna = np.asarray(calibration.count_gdna_region, dtype=np.float64)
    rna = np.asarray(calibration.count_rna_region, dtype=np.float64)
    efflen = np.asarray(calibration.gdna_region_eff_len, dtype=np.float64)

    floor = CONSTANTS.calibration.track_floor
    density = np.where(efflen > floor, gdna / np.maximum(efflen, floor), 0.0)
    total = gdna + rna
    frac = np.where(total > floor, gdna / np.maximum(total, floor), 0.0)

    categories = [str(x) for x in ref_names]
    ref = pd.Categorical.from_codes(np.asarray(region_arrays.ref_id, dtype=np.int64), categories)

    return pd.DataFrame(
        {
            "ref": ref,
            "start": np.asarray(region_arrays.start, dtype=np.int64),
            "end": np.asarray(region_arrays.end, dtype=np.int64),
            "count_gdna": gdna,
            "count_rna": rna,
            "gdna_density": density.astype(np.float64),
            "fraction_gdna": frac.astype(np.float64),
        }
    )


def write_bedgraph(
    track: pd.DataFrame,
    path,
    column: str = "gdna_density",
    *,
    track_name: str = "rigel_gdna_density",
) -> None:
    """Write a bedGraph of one track column (default gDNA density) for genome browsers.

    bedGraph is ``chrom  start  end  value`` (tab-separated). A ``track`` header
    line names it for IGV/UCSC. Rows follow the track's genomic order.
    """
    with open(path, "w") as fh:
        fh.write(
            f'track type=bedGraph name="{track_name}" description="Rigel per-region gDNA {column}"\n'
        )
        refs = track["ref"].astype(str).to_numpy()
        starts = track["start"].to_numpy()
        ends = track["end"].to_numpy()
        vals = track[column].to_numpy()
        for r, s, e, v in zip(refs, starts, ends, vals):
            fh.write(f"{r}\t{int(s)}\t{int(e)}\t{v:.6g}\n")
