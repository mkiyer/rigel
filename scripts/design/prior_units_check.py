"""⭐⭐ T3 — IS THE EM PRIOR A FRAGMENT COUNT? The units gate, off a plain scan cache, no oracle.

       Unit gates: `tests/calibration/test_prior_units.py`, whose docstring names THIS file as T3's home:
       *"The end-to-end conservation check against ``region_start_count`` is T3, and it lives in
       `scripts/design/prior_units_check.py`."*

⭐ **THE CHECK.** ``assemble_priors`` hands the EM two additive pseudocounts that are added directly to
its own fragment counts (``G = n_gdna + a_g``, `em_solver.cpp:apply_grouped_prior_update`). So
``Σ_loci (a_g + a_r)`` must be commensurate with the number of fragments those loci actually contain —
and ``region_start_count``, the accumulator's one non-tautological invariant (one increment per accepted
fragment, at the region holding its first covered base), is exactly that number with no model in it.

⛔ **THE VERDICT IS A DIRECTION, NOT A TARGET.** The prior arbitrates only the UNSPLICED population and
only inside a locus, so ``prior/frag`` below 1 is EXPECTED and is not a score. What must never happen is
``prior/frag`` above 1 — a prior stronger than the whole library, which is the failure this exists to
catch. ``prior/unspl`` is the tighter reading: the same total against the deconvolved unspliced fragment
count over ALL objects, so the shortfall it shows is intergenic mass falling outside every locus.

⛔⛔ **WHAT THIS FILE USED TO BE, AND WHY THAT HALF IS GONE RATHER THAN REPAIRED (2026-08-17).** It ran an
in-process A/B of two prior formulations — arm A the raw object-incidence sum, arm B
``(Σ share·mass_c / Σ share·A_c) · span_bp``, a DENSITY re-integrated over the genomic span — and it
imported ``_component_region_arrays`` from `rigel.calibration.priors` to build them. That symbol was
deleted with the formulation: the shipped prior is now a **conserved fragment count**,
``Σ_regions share·mass_c + Σ_boundaries share·mass_c·q``, with no density and no span re-integration
anywhere in it. ⛔ So both arms scored quantities nothing computes, and widening the import to the
surviving ``_region_locus_shares`` / ``_boundary_locus_shares`` would have rebuilt a dead formula out of
live parts — `TRAPS: a-gate-that-reconstructs`. ⭐ **T3 needed neither arm.** It is asked of the SHIPPED
``assemble_priors`` output directly, which is also the only quantity the EM reads.

⚠ **Score against the payload's own fragment counts, never against a remembered baseline table**
(`TRAPS: re-record-the-baseline`). Every number below is derived from the cache in front of it.

    python scripts/design/prior_units_check.py --index IDX --cache-root CACHE_DIR
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.priors import assemble_priors  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402


def _unspliced_fragments(cal) -> float:
    """The deconvolved UNSPLICED fragment count over EVERY object — the population the prior arbitrates.

    ⛔ The crossing axis is converted by the accumulator's own ``q = boundary_mass_per_crossing`` exactly
    as ``assemble_priors`` converts it, because a crossing deposits ``+1`` per boundary crossed and is
    therefore NOT a fragment count until it is multiplied by ``q``. Summing the raw crossing mass here
    would inflate the denominator and make ``prior/unspl`` read low for a reason that is arithmetic.

    ⚠ ``mass_rna_spliced_boundary`` is subtracted for the reason ``assemble_priors`` subtracts it: a
    spliced fragment has no gDNA candidate in the EM, so it is not part of the split the prior arbitrates.
    """
    q = np.asarray(cal.boundary_mass_per_crossing, np.float64)
    crossing = np.asarray(cal.mass_gdna_boundary, np.float64) + np.maximum(
        np.asarray(cal.mass_rna_boundary, np.float64)
        - np.asarray(cal.mass_rna_spliced_boundary, np.float64),
        0.0,
    )
    contained = np.asarray(cal.mass_gdna_region, np.float64) + np.asarray(
        cal.mass_rna_region, np.float64
    )
    return float(contained.sum() + (crossing * q).sum())


def check_one(index: TranscriptIndex, cache_dir: Path) -> dict:
    from rigel.calibration.region_arrays import RegionArrays

    cache = read_scan_cache(cache_dir, index)
    inputs = calibration_inputs(cache, index)
    cal = calibrate(**inputs, config=CalibrationConfig())
    ra: RegionArrays = inputs["region_arrays"]

    # ⭐ The locus set, WITHOUT running the EM: `build_multi_loci` needs the scored-fragment EM data, so
    # for a units check we use the index's own locus decomposition over the same region axis. T3 is a
    # question about MAGNITUDE, and the shipped assembler is fed the same objects either way.
    multi_loci = _index_multi_loci(index, ra)
    priors = assemble_priors(cal, ra, multi_loci)
    a_g = np.asarray(priors.gdna_prior_count, np.float64)
    a_r = np.asarray(priors.rna_prior_count, np.float64)

    frags = float(np.asarray(cache.payload.region_start_count, np.float64).sum())
    total = float(a_g.sum() + a_r.sum())
    unspliced = _unspliced_fragments(cal)
    return {
        "condition": cache_dir.name,
        "n_loci": len(multi_loci),
        "fragments": frags,
        "unspliced": unspliced,
        "prior_total": total,
        "prior_per_frag": total / max(frags, 1e-9),
        "prior_per_unspliced": total / max(unspliced, 1e-9),
        "f_gdna": float(a_g.sum() / max(total, 1e-9)),
    }


def _index_multi_loci(index: TranscriptIndex, region_arrays):
    """One MultiLocus per contiguous run of non-intergenic regions — a model-free locus decomposition.

    ⚠ Deliberately NOT `locus.build_multi_loci`: that needs the scored-fragment EM data, and this check
    is about UNITS, not about the EM's locus grouping. The assembler sees the same objects either way, so
    the total it produces is unaffected by which decomposition names the blocks.
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--cache-root", type=Path, required=True)
    ap.add_argument("--conditions", nargs="*", default=None)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(args.index))
    dirs = sorted(d for d in args.cache_root.iterdir() if (d / "payload.npz").exists())
    if args.conditions:
        dirs = [d for d in dirs if d.name in set(args.conditions)]
    if not dirs:
        raise SystemExit(
            f"⛔ no scan cache under {args.cache_root} (a condition directory holding `payload.npz`). "
            f"Build one with `scripts/design/build_scan_cache.py` or `scripts/sim/panel.py cache`."
        )
    rows = [check_one(index, d) for d in dirs]

    print(
        f"\n{'condition':46s} {'loci':>7s} {'fragments':>12s} {'unspliced':>12s} {'prior':>12s} "
        f"{'prior/frag':>10s} {'prior/unspl':>11s} {'f_gdna':>8s}"
    )
    print("-" * 46 + " " + "-" * 78)
    over = []
    for r in rows:
        flag = ""
        if r["prior_per_frag"] > 1.0:
            flag = "  ⛔ T3 FAILS"
            over.append(r["condition"])
        print(
            f"{r['condition']:46s} {r['n_loci']:7,d} {r['fragments']:12,.0f} "
            f"{r['unspliced']:12,.0f} {r['prior_total']:12,.0f} "
            f"{r['prior_per_frag']:10.3f} {r['prior_per_unspliced']:11.3f} {r['f_gdna']:8.4f}{flag}"
        )
    print(
        "\nT3: `prior/frag` is Σ(a_g + a_r) against the accepted-fragment count. The prior arbitrates only"
        "\n    the UNSPLICED fragments AND only inside a locus, so a ratio below 1 is expected and is not a"
        "\n    score; what must NOT happen is a ratio above 1 — a prior stronger than the whole library."
        "\n    `prior/unspl` is the same total against the deconvolved unspliced count over ALL objects, so"
        "\n    its shortfall is the intergenic mass that falls outside every locus."
    )
    if over:
        print(f"\n⛔⛔ T3 FAILS on {len(over)} condition(s): {', '.join(over)}")
        return 1
    print(f"\n✅ T3 holds on all {len(rows)} condition(s): the prior is never stronger than the library.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
