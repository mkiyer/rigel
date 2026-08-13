"""Genome-wide gDNA track from the calibration result.

Calibration solves each region for a gDNA/RNA split; joined to the region
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
    which is the same order the ``CalibrationResult`` per-REGION arrays are aligned
    to. ``gdna_density`` is the density-correct contained gDNA mass per bp
    (mass / ``gdna_region_eff_len``); ``gdna_frac`` is gDNA's share of the region's
    contained mass.

    ⚠ **Contained only, deliberately.** A contiguous boundary is a 0-bp boundary: it has no genomic extent to
    occupy a row of a track, and attributing its crossing mass to a flank region would report a density
    at a position where that mass was never contained.
    """
    gdna = np.asarray(calibration.mass_gdna_region, dtype=np.float64)
    rna = np.asarray(calibration.mass_rna_region, dtype=np.float64)
    efflen = np.asarray(calibration.gdna_region_eff_len, dtype=np.float64)

    density = np.where(efflen > _EPS, gdna / np.maximum(efflen, _EPS), 0.0)
    total = gdna + rna
    frac = np.where(total > _EPS, gdna / np.maximum(total, _EPS), 0.0)

    categories = [str(x) for x in ref_names]
    ref = pd.Categorical.from_codes(np.asarray(region_arrays.ref_id, dtype=np.int64), categories)

    return pd.DataFrame(
        {
            "ref": ref,
            "start": np.asarray(region_arrays.start, dtype=np.int64),
            "end": np.asarray(region_arrays.end, dtype=np.int64),
            "gdna_mass": gdna,
            "rna_mass": rna,
            "gdna_density": density.astype(np.float64),
            "gdna_frac": frac.astype(np.float64),
        }
    )


#: separation (nats) above which capture is flagged "enriched" — a soft display
#: descriptor; the numbers are always reported (no hard verdict).
_CAPTURE_ENRICHED_NATS = 1.0
_KDE_N_GRID = 320


def capture_summary(track: pd.DataFrame | None, *, with_curve: bool = False) -> dict | None:
    """Mass-weighted capture-enrichment summary from the per-region gDNA track.

    On hybrid-capture RNA-seq the on-target regions are a small minority of regions
    but carry the captured gDNA *mass* at high density; the many off-target
    (expressed) genes dominate by count. So an equal-weight density KDE is
    unimodal, while a gDNA-mass-weighted KDE develops a high-density on-target
    shoulder. We report:

    * ``background_mode`` — the dominant density peak by region **count** (the
      typical region);
    * ``enriched_mode`` — the highest-density peak of the mass-weighted KDE (the
      on-target shoulder); their gap is the **peak-to-peak** fold-change;
    * ``mass_frac_ontarget`` — the fraction of gDNA **mass** sitting between the
      two peaks' midpoint and above (how much material is actually on-target).

    The peak-to-peak fold and the on-target mass fraction answer different
    questions (how enriched vs how much) — both are surfaced; no pass/fail
    verdict. Returns ``None`` if the track is empty / uninformative / SciPy is
    unavailable. Set ``with_curve`` to also return the plottable KDE curves.
    """
    if track is None or len(track) == 0:
        return None
    try:
        from scipy.stats import gaussian_kde
    except ImportError:  # pragma: no cover
        return None

    dens = np.asarray(track["gdna_density"], dtype=np.float64)
    gmass = np.asarray(track["gdna_mass"], dtype=np.float64)
    keep = (dens > 0) & (gmass > 0)
    if int(keep.sum()) < 20:
        return None
    log_rho = np.log(dens[keep])
    w = gmass[keep]
    if np.ptp(log_rho) < 1e-6:
        return None

    grid = np.linspace(float(log_rho.min()), float(log_rho.max()), _KDE_N_GRID)
    try:
        kde_count = gaussian_kde(log_rho)
        y_count = kde_count(grid)
        y_mass = gaussian_kde(log_rho, weights=w)(grid)
    except Exception:  # pragma: no cover
        return None

    background_mode = float(grid[int(np.argmax(y_count))])
    count_median = float(np.median(log_rho))
    # Enriched mode = highest-x local maximum of the mass-weighted KDE that sits
    # clearly above the median density; falls back to no enrichment. NOTE the
    # enriched mode (not the depleted mode) is the robust target — GC / mappability
    # depress the depleted mode artificially, and small panels + KDE bandwidth can
    # manufacture small modes, so we also report the bandwidth + on-target mass so
    # the modes can be judged rather than trusted blindly.
    peak = float(y_mass.max())
    hi_modes = [
        float(grid[i])
        for i in range(1, _KDE_N_GRID - 1)
        if y_mass[i] > y_mass[i - 1]
        and y_mass[i] >= y_mass[i + 1]
        and y_mass[i] > 0.02 * peak
        and grid[i] > count_median + _CAPTURE_ENRICHED_NATS
    ]
    enriched = bool(hi_modes)
    enriched_mode = max(hi_modes) if enriched else count_median

    # Two fold references (peak-to-peak uses the fragile depleted mode; vs-median
    # uses the robust median). On-target mass fraction uses the median→enriched
    # midpoint as the valley threshold (median is the robust low reference).
    sep_peak = enriched_mode - background_mode
    sep_median = enriched_mode - count_median
    midpoint = (count_median + enriched_mode) / 2.0
    mass_frac_ontarget = float(w[log_rho >= midpoint].sum() / w.sum()) if enriched else 0.0

    out = {
        "n_regions": int(keep.sum()),
        "enriched": enriched,
        "count_median_log_rho": round(count_median, 4),
        "background_mode_log_rho": round(background_mode, 4),
        "enriched_mode_log_rho": round(enriched_mode, 4),
        "fold_peak_to_peak": round(float(np.exp(sep_peak)), 2),
        "fold_vs_median": round(float(np.exp(sep_median)), 2),
        "separation_peak_nats": round(sep_peak, 4),
        "separation_median_nats": round(sep_median, 4),
        "mass_frac_ontarget": round(mass_frac_ontarget, 4),
        "kde_bandwidth_factor": round(float(kde_count.factor), 4),
    }
    if with_curve:
        cmax = float(y_count.max()) or 1.0
        mmax = peak or 1.0
        out["curve"] = [
            {
                "log_rho": round(float(grid[i]), 4),
                "by_count": round(float(y_count[i] / cmax), 5),
                "by_mass": round(float(y_mass[i] / mmax), 5),
            }
            for i in range(_KDE_N_GRID)
        ]
    return out


def write_bedgraph(
    track: pd.DataFrame,
    path,
    column: str = "gdna_density",
    *,
    track_name: str = "rigel_gdna_density",
) -> None:
    """Write a bedGraph of one track column (default gDNA density) for genome browsers.

    bedGraph is ``chrom  start  end  value`` (tab-separated). A ``track`` header
    boundary names it for IGV/UCSC. Rows follow the track's genomic order.
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
