"""Is the EM prior in fragment units? The units gate, run off a plain scan cache with no oracle; it is
the home of the end-to-end conservation check (T3) that `tests/calibration/test_prior_units.py`
names. ``assemble_priors`` hands the EM two additive pseudocounts that are added directly to its own
fragment counts, so ``sum(a_g + a_r)`` over the loci must be commensurate with the number of
fragments those loci contain — and ``region_start_count``, the accumulator's one model-free
invariant (one increment per accepted fragment, at the region holding its first covered base), is
that number. Per cached condition it calibrates under the shipped config, assembles the shipped
priors over a model-free locus decomposition of the index (contiguous runs of genic regions; the EM's
own locus grouping needs scored fragments and the assembler sees the same objects either way), and
reports ``prior/frag`` and ``prior/unspl``. The verdict is a direction, not a target: the prior
arbitrates only the unspliced population and only inside a locus, so a ratio below 1 is expected
and is not a score; a ratio above 1 — a prior stronger than the whole library — is the failure this
exists to catch. ``prior/unspl`` is the tighter reading, the same total against the deconvolved
unspliced count over all objects, so its shortfall is the intergenic mass outside every locus. Every
number is derived from the cache in front of it, never from a remembered table.

Usage::

    python scripts/design/prior_units_check.py --index IDX --cache-root CACHE_DIR
    python scripts/design/prior_units_check.py --index IDX --cache-root CACHE_DIR --conditions NAME ...
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

    The crossing axis is converted by the accumulator's own ``q = boundary_mass_per_crossing`` exactly
    as ``assemble_priors`` converts it, because a crossing deposits ``+1`` per boundary crossed and is
    therefore NOT a fragment count until it is multiplied by ``q``. Summing the raw crossing mass here
    would inflate the denominator and make ``prior/unspl`` read low for a reason that is arithmetic.

    ``count_rna_spliced_boundary`` is subtracted for the reason ``assemble_priors`` subtracts it: a
    spliced fragment has no gDNA candidate in the EM, so it is not part of the split the prior arbitrates.
    """
    q = np.asarray(cal.boundary_mass_per_crossing, np.float64)
    crossing = np.asarray(cal.count_gdna_boundary, np.float64) + np.maximum(
        np.asarray(cal.count_rna_boundary, np.float64)
        - np.asarray(cal.count_rna_spliced_boundary, np.float64),
        0.0,
    )
    contained = np.asarray(cal.count_gdna_region, np.float64) + np.asarray(
        cal.count_rna_region, np.float64
    )
    return float(contained.sum() + (crossing * q).sum())


def check_one(index: TranscriptIndex, cache_dir: Path) -> dict:
    from rigel.calibration.region_arrays import RegionArrays

    cache = read_scan_cache(cache_dir, index)
    inputs = calibration_inputs(cache, index)
    cal = calibrate(**inputs, config=CalibrationConfig())
    ra: RegionArrays = inputs["region_arrays"]

    # the locus set WITHOUT running the EM: `build_multi_loci` needs the scored-fragment EM data, so a
    # units check uses the index's own locus decomposition over the same region axis. This is a
    # question about MAGNITUDE, and the shipped assembler is fed the same objects either way.
    multi_loci = _index_multi_loci(index, ra)
    priors = assemble_priors(cal, ra, multi_loci)
    a_g = np.asarray(priors.gdna_prior_count, np.float64)
    a_r = np.asarray(priors.rna_prior_count, np.float64)

    frags = float(np.asarray(inputs["payload"].region_start_count, np.float64).sum())
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

    Deliberately NOT `locus.build_multi_loci`: that needs the scored-fragment EM data, and this check
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
