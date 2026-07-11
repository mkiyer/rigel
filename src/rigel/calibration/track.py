"""Genome-wide gDNA track from the calibration result.

Calibration solves each region node for a gDNA/RNA split; joined to the region
coordinates this yields a per-region gDNA level across the genome — a QC track
that traces contamination and, on capture libraries, on-target enrichment.

This is pure persistence: every value is already in ``CalibrationResult`` (the
per-region deconvolved masses + the density-correct effective length) and
``RegionArrays`` (the coordinates). No new calibration computation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

_EPS = 1e-9


def build_gdna_track(calibration, region_arrays, ref_names) -> pd.DataFrame:
    """Per-region gDNA track: ``(ref, start, end, gdna_mass, rna_mass, gdna_density, gdna_frac)``.

    Rows are in genomic order (``RegionArrays`` is sorted by ``(ref_id, start)``),
    which is the same order the ``CalibrationResult`` per-region arrays are aligned
    to. ``gdna_density`` is the density-correct contained gDNA mass per bp
    (mass / ``gdna_region_eff_len``); ``gdna_frac`` is gDNA's share of the region's
    contained mass.
    """
    gdna = np.asarray(calibration.mass_gdna_contained, dtype=np.float64)
    rna = np.asarray(calibration.mass_rna_contained, dtype=np.float64)
    efflen = np.asarray(calibration.gdna_region_eff_len, dtype=np.float64)

    density = np.where(efflen > _EPS, gdna / np.maximum(efflen, _EPS), 0.0)
    total = gdna + rna
    frac = np.where(total > _EPS, gdna / np.maximum(total, _EPS), 0.0)

    categories = [str(x) for x in ref_names]
    ref = pd.Categorical.from_codes(np.asarray(region_arrays.ref_id, dtype=np.int64), categories)

    return pd.DataFrame({
        "ref": ref,
        "start": np.asarray(region_arrays.start, dtype=np.int64),
        "end": np.asarray(region_arrays.end, dtype=np.int64),
        "gdna_mass": gdna,
        "rna_mass": rna,
        "gdna_density": density.astype(np.float64),
        "gdna_frac": frac.astype(np.float64),
    })


def write_bedgraph(track: pd.DataFrame, path, column: str = "gdna_density", *, track_name: str = "rigel_gdna_density") -> None:
    """Write a bedGraph of one track column (default gDNA density) for genome browsers.

    bedGraph is ``chrom  start  end  value`` (tab-separated). A ``track`` header
    line names it for IGV/UCSC. Rows follow the track's genomic order.
    """
    with open(path, "w") as fh:
        fh.write(f'track type=bedGraph name="{track_name}" description="Rigel per-region gDNA {column}"\n')
        refs = track["ref"].astype(str).to_numpy()
        starts = track["start"].to_numpy()
        ends = track["end"].to_numpy()
        vals = track[column].to_numpy()
        for r, s, e, v in zip(refs, starts, ends, vals):
            fh.write(f"{r}\t{int(s)}\t{int(e)}\t{v:.6g}\n")
