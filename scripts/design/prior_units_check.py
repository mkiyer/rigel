"""P1 gate T3 — is the EM prior a FRAGMENT COUNT? And the in-process A/B against the raw sum.

       Unit gates: `tests/calibration/test_prior_units.py`

⭐ **THE CHECK.** ``assemble_priors`` hands the EM two additive pseudocounts that are added directly to
its own fragment counts (``G = n_gdna + a_g``). So ``Σ_loci (a_g + a_r)`` must be commensurate with the
number of fragments those loci actually contain — and ``region_start_count``, the accumulator's one
non-tautological invariant (one increment per accepted fragment, at the region holding its first covered
base), is exactly that number with no model in it.

⭐ **AND THE A/B IS EXACT AND IN-PROCESS.** Both arms are computed from ONE ``CalibrationResult``, so the
comparison has no simulation noise, no re-scan and no second code path:

    arm A (pre-P1)  = Σ_locus share · mass_c            the raw object-incidence sum
    arm B (P1)      = (Σ share·mass_c / Σ share·A_c) · span_bp

⚠ **Score against the payload's own fragment counts, never against the S5.f baseline table.** The two
defects partly cancel today, so the sign of the error flips on the gdna100 arm — the headline moves
*away* from 0.50 while moving *toward* truth (recorded before the
run). Comparing to the old number would score the cancellation, not the correction.

    python scripts/design/prior_units_check.py --index IDX --cache-root CACHE_DIR
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.priors import _component_region_arrays, _project_regions_to_loci  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402


def _arms(calibration, region_arrays, multi_loci):
    """Both arms from ONE result: (arm_A_gdna, arm_A_rna, arm_B_gdna, arm_B_rna), per locus."""
    g_mass, g_sup, _, _ = _component_region_arrays(
        calibration,
        region_arrays,
        calibration.mass_gdna_region,
        calibration.mass_gdna_edge,
        calibration.gdna_region_eff_len,
        calibration.gdna_edge_eff_len,
    )
    rna_edge = np.maximum(
        np.asarray(calibration.mass_rna_edge, np.float64)
        - np.asarray(calibration.mass_rna_spliced_edge, np.float64),
        0.0,
    )
    r_mass, r_sup, _, _ = _component_region_arrays(
        calibration,
        region_arrays,
        calibration.mass_rna_region,
        rna_edge,
        calibration.rna_region_eff_len,
        calibration.rna_edge_eff_len,
    )
    proj = _project_regions_to_loci(
        region_arrays,
        multi_loci,
        len(multi_loci),
        {
            "g_mass": g_mass,
            "g_sup": g_sup,
            "r_mass": r_mass,
            "r_sup": r_sup,
            "span_bp": np.asarray(region_arrays.region_size_bp, np.float64),
        },
    )

    def dens(m, s):
        return np.divide(m, s, out=np.zeros_like(m), where=s > 0.0) * proj["span_bp"]

    return (
        proj["g_mass"],
        proj["r_mass"],
        dens(proj["g_mass"], proj["g_sup"]),
        dens(proj["r_mass"], proj["r_sup"]),
    )


def check_one(index: TranscriptIndex, cache_dir: Path) -> dict:
    from rigel.calibration.region_arrays import RegionArrays

    cache = read_scan_cache(cache_dir, index)
    inputs = calibration_inputs(cache, index)
    cal = calibrate(**inputs, config=CalibrationConfig())
    ra: RegionArrays = inputs["region_arrays"]

    # ⭐ The locus set, WITHOUT running the EM: `build_multi_loci` needs the scored-fragment EM data, so
    # for a units check we use the index's own locus decomposition over the same region axis. What matters
    # is that both arms see the SAME loci — the comparison is internal.
    multi_loci = _index_multi_loci(index, ra)
    a_g, a_r, b_g, b_r = _arms(cal, ra, multi_loci)

    frags = float(np.asarray(cache.payload.region_start_count, np.float64).sum())
    # the unspliced population the prior arbitrates: contained + unspliced crossings, deconvolved
    unspliced_mass = float(
        np.asarray(cal.mass_gdna_region).sum()
        + np.asarray(cal.mass_rna_region).sum()
        + np.asarray(cal.mass_gdna_edge).sum()
        + np.asarray(cal.mass_rna_edge).sum()
        - np.asarray(cal.mass_rna_spliced_edge).sum()
    )
    return {
        "condition": cache_dir.name,
        "fragments": frags,
        "unspliced_incidences": unspliced_mass,
        "A_total": float(a_g.sum() + a_r.sum()),
        "B_total": float(b_g.sum() + b_r.sum()),
        "A_fgdna": float(a_g.sum() / max(a_g.sum() + a_r.sum(), 1e-9)),
        "B_fgdna": float(b_g.sum() / max(b_g.sum() + b_r.sum(), 1e-9)),
    }


def _index_multi_loci(index: TranscriptIndex, region_arrays):
    """One MultiLocus per contiguous run of non-intergenic regions — a model-free locus decomposition.

    ⚠ Deliberately NOT `locus.build_multi_loci`: that needs the scored-fragment EM data, and this check
    is about UNITS, not about the EM's locus grouping. Both arms see the same loci, so the comparison is
    unaffected by which decomposition is used.
    """
    from rigel.calibration.signature import (
        BIT_EXON_NEG,
        BIT_EXON_POS,
        BIT_INTRON_NEG,
        BIT_INTRON_POS,
    )
    from rigel.locus import Locus, MultiLocus

    bits = BIT_EXON_POS | BIT_EXON_NEG | BIT_INTRON_POS | BIT_INTRON_NEG
    genic = (np.asarray(region_arrays.signature).astype(np.int64) & bits) != 0
    ref = np.asarray(region_arrays.ref_id)
    start = np.asarray(region_arrays.start)
    end = np.asarray(region_arrays.end)

    out, i, n = [], 0, genic.shape[0]
    while i < n:
        if not genic[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and genic[j + 1] and ref[j + 1] == ref[i]:
            j += 1
        lid = len(out)
        out.append(
            MultiLocus(
                multi_locus_id=lid,
                transcript_indices=np.array([], dtype=np.int32),
                unit_indices=np.array([], dtype=np.int32),
                gdna_span=int(end[j] - start[i]),
                loci=(
                    Locus(
                        ref=str(int(ref[i])),
                        ref_id=int(ref[i]),
                        start=int(start[i]),
                        end=int(end[j]),
                    ),
                ),
            )
        )
        i = j + 1
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--cache-root", type=Path, required=True)
    ap.add_argument("--conditions", nargs="*", default=None)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(args.index))
    dirs = sorted(d for d in args.cache_root.iterdir() if (d / "payload.npz").exists())
    if args.conditions:
        dirs = [d for d in dirs if d.name in set(args.conditions)]
    rows = [check_one(index, d) for d in dirs]

    print(
        f"\n{'condition':46s} {'fragments':>11s} {'A: raw sum':>11s} {'B: P1':>11s} "
        f"{'A/frag':>7s} {'B/frag':>7s}   {'A f_gdna':>9s} {'B f_gdna':>9s}"
    )
    print("-" * 46 + " " + "-" * 74)
    for r in rows:
        print(
            f"{r['condition']:46s} {r['fragments']:11,.0f} {r['A_total']:11,.0f} "
            f"{r['B_total']:11,.0f} {r['A_total'] / r['fragments']:7.3f} "
            f"{r['B_total'] / r['fragments']:7.3f}   {r['A_fgdna']:9.4f} {r['B_fgdna']:9.4f}"
        )
    print(
        "\nT3: 'A/frag' and 'B/frag' are the prior's total against the accepted-fragment count. "
        "\n    The prior arbitrates only the UNSPLICED fragments, so a ratio below 1 is expected; "
        "\n    what must NOT happen is a ratio above 1 (a prior stronger than the whole library)."
    )


if __name__ == "__main__":
    main()
