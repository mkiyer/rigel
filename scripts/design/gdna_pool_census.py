"""The FOUR gDNA fragment-length pools, each against its OWN opportunity and against truth.

    Stage A of `docs/SUCCESS.md`: is the accumulator's gDNA length model accurate?

⭐ **The four pools, and why four.** The accumulator bins every deposited fragment into at most one pool
that is **pure by construction** (`tests/native/_accumulator_reference.py`, `FragmentPool`):

| pool | rule | dominant when |
|---|---|---|
| `DNA_INTERGENIC` | contained in exactly one **intergenic** region | capture OFF |
| `DNA_INTRONIC` | contained in exactly one **intronic** region | capture OFF |
| `DNA_INTRON_EXON` | crosses exactly one line, flanks {intron, exon} | capture ON |
| `DNA_INTERGENIC_EXON` | crosses exactly one line, flanks {intergenic, exon} | capture ON |

Off capture the library is spread over the genome and most gDNA sits wholly inside a large intergenic or
intronic region. Under capture the surviving gDNA sits beside a probe, and a fragment beside a probe
*reaches* the exon boundary — so it stops being contained and becomes crossing. **Both regimes need
covering, which is why the model is fitted from all four.**

⛔⛔ **THE POOLS MUST NOT BE POOLED RAW, AND THIS SCRIPT IS WHERE THAT IS VISIBLE.** A contained pool's
opportunity falls with length — `(ell - w + 1)+` — and a crossing pool's *rises* — roughly `(w - 1)+`.
Summing four differently-tilted histograms and applying one divisor is the defect that read a gDNA mean
of 146.05 where the pure contained pool said 88.0. Each pool is divided by **its own** opportunity.

⭐ **AND THE DIVISOR IS A PROBABILITY, NOT A COUNT.** `count(w)/A(w)` recovers the distribution lengths
were *drawn* from; every consumer needs the one the library *realises*, so the divisor is
`pi(w) = A(w) / T(w)` with `T(w) = sum_refs (L_ref - w + 1)+` the total admissible starts genome-wide.
On whole chromosomes `T(w)` is flat to ~1e-5 and the two forms coincide; on a short reference they do not,
and the probability form is the correct one either way.

⭐ **THE COMBINATION IS DERIVED, NOT CHOSEN.** Summing counts and summing opportunities,
`f(w) ~ [sum_p count_p(w)] * T(w) / [sum_p A_p(w)]`, is algebraically the **opportunity-weighted average**
of the four de-tilted pools — and under Poisson counts `Var(count_p) ∝ A_p`, so weights proportional to
`A_p` are exactly inverse-variance. There is no tunable weight.

    python scripts/design/gdna_pool_census.py [--index DIR] [--pilot DIR] [--suite DIR]
        [--conditions ...] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_PILOT = _RUNS / "suite" / "pilot" / "scan_cache"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: Coarse region types, as `signature.coarse_type_array` emits them.
TYPE_INTERGENIC, TYPE_INTRON, TYPE_EXON = 0, 1, 2

POOL_NAMES = ("DNA_INTERGENIC", "DNA_INTRONIC", "DNA_INTRON_EXON", "DNA_INTERGENIC_EXON", "RNA_SPLICED")
GDNA_POOLS = (0, 1, 2, 3)


# ── the opportunity functions ───────────────────────────────────────────────────────────────────


def _ramp(values: np.ndarray, max_w: int) -> np.ndarray:
    """``sum_i (w - 1 - values[i])+`` for every ``w`` in ``[0, max_w]``, in O(max_w + n).

    Values above ``max_w - 2`` can never contribute (they need ``x <= w - 2``), so they are clipped
    into a bin the cumulative sums never read.
    """
    clipped = np.clip(values, 0, max_w + 1).astype(np.int64)
    g = np.bincount(clipped, minlength=max_w + 2).astype(np.float64)
    index = np.arange(len(g), dtype=np.float64)
    g1, g2 = np.cumsum(g), np.cumsum(g * index)
    w = np.arange(max_w + 1, dtype=np.int64)
    out = np.zeros(max_w + 1, dtype=np.float64)
    ok = w >= 2
    m = w[ok] - 2
    out[ok] = (w[ok] - 1) * g1[m] - g2[m]
    return out


def contained_opportunity(region_lengths: np.ndarray, max_w: int) -> np.ndarray:
    """``sum_n (ell_n - w + 1)+`` — the admissible starts for a fragment contained in one of these regions.

    ⭐ Exact, and O(max_w + n) rather than O(max_w * n): regions longer than ``max_w`` always contribute,
    so they reduce to ``sum(ell) - (w - 1) * count``; the rest come from a length histogram's reverse
    cumulative sums.
    """
    lengths = np.asarray(region_lengths, dtype=np.int64)
    w = np.arange(max_w + 1, dtype=np.float64)

    big = lengths > max_w
    out = float(lengths[big].sum()) - (w - 1.0) * float(big.sum())

    small = lengths[~big]
    if small.size:
        h = np.bincount(small, minlength=max_w + 2).astype(np.float64)
        index = np.arange(len(h), dtype=np.float64)
        # Reverse cumulative sums: s1[l] = sum_{k >= l} h[k], s2[l] = sum_{k >= l} k * h[k].
        s1 = np.cumsum(h[::-1])[::-1]
        s2 = np.cumsum((h * index)[::-1])[::-1]
        idx = np.arange(max_w + 1)
        out = out + (s2[idx] - (w - 1.0) * s1[idx])
    return np.maximum(out, 0.0)


def crossing_opportunity(left: np.ndarray, right: np.ndarray, max_w: int) -> np.ndarray:
    """``sum_lines`` admissible starts that cross THIS line and no other.

    For a line at ``p`` with flanking region lengths ``a`` (left) and ``b`` (right), a fragment
    ``[s, s+w)`` crosses ``p`` for ``w - 1`` starts; of those it also crosses the previous line for
    ``(w-1-a)+`` and the next for ``(w-1-b)+``, and both for ``(w-1-a-b)+``. So

        A(w) = (w-1)+ - (w-1-a)+ - (w-1-b)+ + (w-1-a-b)+

    ⭐ The two nearest lines are the only ones that need excluding: a fragment is an interval containing
    ``p``, so if it reaches any line left of ``p-a`` it must also cross ``p-a``. The inclusion-exclusion
    over the two neighbours is therefore **exact**, not an approximation.

    ⚠ The reference ends need no special case. The partition cuts at 0 and at ``L_ref``, so the outermost
    region's length *is* the distance to the wall and the same subtraction removes the impossible starts.
    """
    w = np.arange(max_w + 1, dtype=np.float64)
    base = float(len(left)) * np.maximum(w - 1.0, 0.0)
    out = base - _ramp(left, max_w) - _ramp(right, max_w) + _ramp(left + right, max_w)
    return np.maximum(out, 0.0)


def total_opportunity(ref_lengths: np.ndarray, max_w: int) -> np.ndarray:
    """``T(w) = sum_refs (L_ref - w + 1)+`` — every admissible gDNA start in the reference."""
    return contained_opportunity(ref_lengths, max_w)


def pool_opportunities(index, max_w: int) -> dict[str, np.ndarray]:
    """The four gDNA pools' opportunities plus ``T``, all derived from the region partition alone."""
    from rigel.calibration.splice_graph import build_region_partition_arrays

    cuts, ref_cut_offsets, region_types = build_region_partition_arrays(index)
    cuts = np.asarray(cuts, dtype=np.int64)
    offsets = np.asarray(ref_cut_offsets, dtype=np.int64)
    types = np.asarray(region_types, dtype=np.int64)

    region_lengths: list[np.ndarray] = []
    line_left: list[np.ndarray] = []
    line_right: list[np.ndarray] = []
    line_types: list[np.ndarray] = []
    ref_lengths: list[int] = []
    region_base = 0
    for r in range(len(offsets) - 1):
        lo, hi = int(offsets[r]), int(offsets[r + 1])
        if hi - lo < 2:
            continue
        ref_cuts = cuts[lo:hi]
        lengths = np.diff(ref_cuts)
        n_regions = len(lengths)
        ref_lengths.append(int(ref_cuts[-1]))
        region_lengths.append(lengths)
        # Interior lines only: line i sits between region i-1 and region i, for i in [1, n_regions).
        if n_regions >= 2:
            line_left.append(lengths[:-1])
            line_right.append(lengths[1:])
            pair = np.sort(
                np.stack([types[region_base : region_base + n_regions - 1], types[region_base + 1 : region_base + n_regions]]),
                axis=0,
            )
            line_types.append(pair)
        region_base += n_regions

    all_lengths = np.concatenate(region_lengths)
    all_types = types[: len(all_lengths)]
    left = np.concatenate(line_left) if line_left else np.zeros(0, np.int64)
    right = np.concatenate(line_right) if line_right else np.zeros(0, np.int64)
    pairs = np.concatenate(line_types, axis=1) if line_types else np.zeros((2, 0), np.int64)

    intron_exon = (pairs[0] == TYPE_INTRON) & (pairs[1] == TYPE_EXON)
    intergenic_exon = (pairs[0] == TYPE_INTERGENIC) & (pairs[1] == TYPE_EXON)

    return {
        "DNA_INTERGENIC": contained_opportunity(all_lengths[all_types == TYPE_INTERGENIC], max_w),
        "DNA_INTRONIC": contained_opportunity(all_lengths[all_types == TYPE_INTRON], max_w),
        "DNA_INTRON_EXON": crossing_opportunity(left[intron_exon], right[intron_exon], max_w),
        "DNA_INTERGENIC_EXON": crossing_opportunity(
            left[intergenic_exon], right[intergenic_exon], max_w
        ),
        "T": total_opportunity(np.asarray(ref_lengths, dtype=np.int64), max_w),
        "_census": {
            "regions": int(len(all_lengths)),
            "intergenic_regions": int((all_types == TYPE_INTERGENIC).sum()),
            "intronic_regions": int((all_types == TYPE_INTRON).sum()),
            "intron_exon_lines": int(intron_exon.sum()),
            "intergenic_exon_lines": int(intergenic_exon.sum()),
        },
    }


# ── moments ─────────────────────────────────────────────────────────────────────────────────────


def moments(counts: np.ndarray) -> tuple[float, float, float]:
    """``(n, mean, sd)`` of a 1-bp histogram. ``n`` is the raw count even when weights are fractional."""
    n = float(counts.sum())
    if n <= 0:
        return 0.0, float("nan"), float("nan")
    w = np.arange(len(counts), dtype=np.float64)
    mean = float((counts * w).sum() / n)
    var = float((counts * (w - mean) ** 2).sum() / n)
    return n, mean, var**0.5 if var > 0 else 0.0


def truth_moments(condition_dir: Path, kind: str, max_w: int) -> tuple[float, float, float]:
    path = condition_dir / "truth_fragment_lengths.tsv"
    counts = np.zeros(max_w + 1, dtype=np.float64)
    with open(path) as handle:
        next(handle)
        for line in handle:
            row_kind, length, count, _frac = line.rstrip("\n").split("\t")
            if row_kind == kind and 0 <= int(length) <= max_w:
                counts[int(length)] += float(count)
    return moments(counts)


def detilt(counts: np.ndarray, opportunity: np.ndarray, total: np.ndarray) -> np.ndarray:
    """``count(w) / pi(w)`` with ``pi = A / T`` — zero where the pool had no opportunity at all."""
    pi = np.where(opportunity > 0, opportunity / np.maximum(total, 1e-300), 0.0)
    return np.where(pi > 0, counts / np.maximum(pi, 1e-300), 0.0)


def pct(value: float, truth: float) -> str:
    if not np.isfinite(value) or not np.isfinite(truth) or truth == 0:
        return "     — "
    return f"{100 * (value / truth - 1):+7.2f}%"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilot", type=Path, default=DEFAULT_PILOT)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--suite", type=Path, default=None)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--no-drain", action="store_true", help="measure pass 1 only, before the side buffer drains")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    suite = args.suite or args.pilot.parent
    from rigel.index import TranscriptIndex
    from rigel.pipeline import _drain_side_buffer
    from rigel.scan_cache import read_scan_cache

    index = TranscriptIndex.load(str(args.index))
    names = args.conditions or sorted(p.name for p in args.pilot.iterdir() if p.is_dir())

    opportunity: dict[str, np.ndarray] | None = None
    rows = []
    for name in names:
        cache = read_scan_cache(args.pilot / name, index)
        payload = cache.payload
        if not args.no_drain:
            payload = _drain_side_buffer(payload, index, cache.strand_model, seed=args.seed)
        pools = np.asarray(payload.pool_lengths, dtype=np.float64)
        max_w = pools.shape[1] - 1
        if opportunity is None:
            opportunity = pool_opportunities(index, max_w)
            census = opportunity.pop("_census")
            print("partition census (annotation-derived, condition-independent)")
            for key, value in census.items():
                print(f"  {key:<24} {value:>12,}")
            print()

        truth_gdna = truth_moments(suite / name, "gdna", max_w)
        truth_rna = truth_moments(suite / name, "rna", max_w)

        print(f"═══ {name} ═══")
        print(
            f"  {'pool':<22} {'n':>12} {'raw mean':>9} {'raw sd':>8} {'vs truth':>9} "
            f"{'⭐ de-tilted':>12} {'sd':>8} {'vs truth':>9}"
        )
        entry: dict = {"condition": name, "truth_gdna": truth_gdna, "truth_rna": truth_rna, "pools": {}}
        summed_counts = np.zeros(max_w + 1, dtype=np.float64)
        summed_opportunity = np.zeros(max_w + 1, dtype=np.float64)
        for p in GDNA_POOLS:
            pool_name = POOL_NAMES[p]
            counts = pools[p]
            n, mean, sd = moments(counts)
            fitted = detilt(counts, opportunity[pool_name], opportunity["T"])
            _fn, fmean, fsd = moments(fitted)
            entry["pools"][pool_name] = {
                "n": n, "raw_mean": mean, "raw_sd": sd, "detilted_mean": fmean, "detilted_sd": fsd,
            }
            print(
                f"  {pool_name:<22} {n:>12,.0f} {mean:>9.2f} {sd:>8.2f} {pct(mean, truth_gdna[1])} "
                f"{fmean:>12.2f} {fsd:>8.2f} {pct(fmean, truth_gdna[1])}"
            )
            summed_counts += counts
            summed_opportunity += opportunity[pool_name]

        # ⛔ The wrong combination, kept visible on purpose: pool the four histograms and treat them as
        # if they shared one opportunity. This is the defect the four-pool model exists to avoid.
        _rn, raw_mean, raw_sd = moments(summed_counts)
        combined = detilt(summed_counts, summed_opportunity, opportunity["T"])
        _cn, comb_mean, comb_sd = moments(combined)
        # The pair the shipped model uses today, for reference.
        pure = pools[0] + pools[1]
        _pn, pure_mean, pure_sd = moments(pure)
        pure_fit = detilt(pure, opportunity["DNA_INTERGENIC"] + opportunity["DNA_INTRONIC"], opportunity["T"])
        _pfn, pure_fmean, pure_fsd = moments(pure_fit)

        print(f"  {'-' * 100}")
        print(
            f"  {'SHIPPED (2 contained)':<22} {_pn:>12,.0f} {pure_mean:>9.2f} {pure_sd:>8.2f} "
            f"{pct(pure_mean, truth_gdna[1])} {pure_fmean:>12.2f} {pure_fsd:>8.2f} "
            f"{pct(pure_fmean, truth_gdna[1])}"
        )
        print(
            f"  {'⛔ 4 POOLED RAW':<22} {_rn:>12,.0f} {raw_mean:>9.2f} {raw_sd:>8.2f} "
            f"{pct(raw_mean, truth_gdna[1])} {'—':>12} {'—':>8} {'—':>9}"
        )
        print(
            f"  {'⭐ 4, EACH DE-TILTED':<22} {_rn:>12,.0f} {'—':>9} {'—':>8} {'—':>9} "
            f"{comb_mean:>12.2f} {comb_sd:>8.2f} {pct(comb_mean, truth_gdna[1])}"
        )
        print(f"  truth (gdna)          n={truth_gdna[0]:>10,.0f} mean={truth_gdna[1]:.2f} sd={truth_gdna[2]:.2f}")
        print()
        entry.update(
            shipped={"n": _pn, "raw_mean": pure_mean, "detilted_mean": pure_fmean, "detilted_sd": pure_fsd},
            pooled_raw={"n": _rn, "mean": raw_mean, "sd": raw_sd},
            four_detilted={"mean": comb_mean, "sd": comb_sd},
        )
        rows.append(entry)

    if args.json:
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True, default=float))
        print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
