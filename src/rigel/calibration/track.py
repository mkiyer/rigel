"""The genome-wide per-region gDNA track, a QC view of the calibration result (layer 7).

Joined to the region coordinates, calibration's per-region split traces contamination and, on capture libraries,
on-target enrichment. Every value is already in ``CalibrationResult`` and ``RegionArrays``; nothing is computed anew.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import CONSTANTS


def build_gdna_track(calibration, region_arrays, ref_names) -> pd.DataFrame:
    """Per-region gDNA track: ``(ref, start, end, count_gdna, count_rna, gdna_density, fraction_gdna)``.

    Rows are in genomic order, the order of the result's per-region arrays. ``count_gdna`` and ``count_rna`` are the
    region's contained fragments as deconvolved; ``gdna_density`` is ``count_gdna / gdna_region_eff_len``;
    ``fraction_gdna`` is gDNA's share of the contained count. Boundaries are left out: a 0-bp boundary has no extent
    to occupy a row, and assigning its crossings to a flank would report a density where those fragments were never
    contained.
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


def write_bedgraph(track: pd.DataFrame, path) -> None:
    """Write the track's gDNA density as a bedGraph for genome browsers.

    bedGraph is ``chrom  start  end  value`` (tab-separated), under a ``track`` header line that names it for
    IGV/UCSC. Rows follow the track's genomic order.
    """
    with open(path, "w") as fh:
        fh.write(
            'track type=bedGraph name="rigel_gdna_density" description="Rigel per-region gDNA gdna_density"\n'
        )
        refs = track["ref"].astype(str).to_numpy()
        starts = track["start"].to_numpy()
        ends = track["end"].to_numpy()
        vals = track["gdna_density"].to_numpy()
        for r, s, e, v in zip(refs, starts, ends, vals):
            fh.write(f"{r}\t{int(s)}\t{int(e)}\t{v:.6g}\n")
