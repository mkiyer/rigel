"""⛔⛔ RETIRED QUESTION, KEPT ARITHMETIC — the region-length sieve, with its own scope stamp.

⛔⛔⛔ **READ THIS BEFORE READING A NUMBER OFF THIS FILE.** The question it asks — *does the CONTAINED
count carry gDNA-vs-RNA composition information, via the two components' fragment lengths?* — **is the
fragment-length CALIBRATION COMPOSITION CHANNEL, which the owner RETIRED until after 0.8.0** (ruling,
2026-08-14). It is not a candidate, it may not be proposed and it may not be ranked. This file is kept
because it is arithmetic rather than a proposal, and because deleting it is a decision the owner owns:
it sits in `tests/test_scripts_index.py`'s ``UNDOCUMENTED_DEBT`` and is **owed a promote-or-delete
call**. ⚠ Three unrelated things also called "length" are NOT affected by that ruling and must not be
touched in its name: layer 2's ``fl`` / ``effective_length`` / ``capture_eff_length``, the second pass's
``length_likelihood``, and the fl PMFs priced by `length_ceiling.py`.

---

A molecule of length L fits inside a region of length ell in (ell - L + 1) positions, or not at all if
L > ell. So the expected contained count for component c is  rho_c * A_c(ell),  with

    A_c(ell) = E_c[ max(0, ell - L + 1) ]      <- the component's own effective length

If the boundary flux already gives the TOTAL density rho = rho_g + rho_r with no FL model (the Sum 1/L
result), then the contained count adds a SECOND equation, and the two are solvable for the split
exactly when A_g(ell) != A_r(ell).  Discriminability is the ratio A_g/A_r: 1.0 = no information,
0.0 = perfect (only RNA can fit).

⛔⛔ **THE FROZEN MEASUREMENT IS GONE AND WAS NOT REPLACED BY A GUESS (2026-08-17).** This file used to
close by PRINTING *"measured earlier this project: 8.2 % of regions < 10 bp, 56.7 % < 200 bp, median
region length 151 bp"* as though it were live. It was undated, hand-carried, and — as the word "v8" in
its own heading said — a property of **one annotation**, which is `TRAPS: re-record-the-baseline`
committed inside an instrument, where it is hardest to see. ⛔ It is not reproducible here: no human v8
index is on disk in this tree, and the partition that is (the suite carve) is a different one, so
substituting its numbers would have been a fabrication wearing the old sentence's clothes.
⭐ **The replacement is a derivation, not a number**: pass ``--index INDEX_DIR`` and the three statistics
are re-derived from that index's own ``regions.feather``. Without it the section says it is unmeasured
and names `index_census.py`, which re-derives this whole class of number.

    python scripts/design/length_sieve.py
    python scripts/design/length_sieve.py --index ~/Downloads/rigel_runs/suite/rigel_index
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

L = np.arange(1, 1201)


def eff(mu, sd, ell):
    """A_c(ell) = E[max(0, ell - L + 1)] under a clipped-normal fragment length."""
    w = np.exp(-0.5 * ((L - mu) / sd) ** 2)
    w /= w.sum()
    return float((w * np.maximum(ell - L + 1, 0)).sum())


def show(name, mug, sdg, mur, sdr):
    print(f"\n{name}:  gDNA L ~ {mug}+-{sdg}   RNA L ~ {mur}+-{sdr}")
    print(f"  {'region ell':>9} {'A_gdna':>10} {'A_rna':>10} {'A_g/A_r':>9}   discriminability")
    for ell in (25, 50, 100, 150, 167, 200, 250, 300, 500, 1000, 5000):
        ag, ar = eff(mug, sdg, ell), eff(mur, sdr, ell)
        ratio = ag / ar if ar > 1e-12 else float("nan")
        if not np.isfinite(ratio):
            bar = "both components excluded"
        elif ratio < 0.05:
            bar = "**** near-pure RNA sieve"
        elif ratio < 0.5:
            bar = "***  strong"
        elif ratio < 0.85:
            bar = "**   usable"
        elif ratio < 0.97:
            bar = "*    weak"
        else:
            bar = "     none"
        print(f"  {ell:>9} {ag:>10.3f} {ar:>10.3f} {ratio:>9.4f}   {bar}")


def informative_zone(index_dir: Path | None) -> None:
    """How much of an index's region partition sits in the informative zone — DERIVED, never quoted.

    ⛔ Reads ``regions.feather`` off the index directory, exactly as `index_census.py` does, so the two
    cannot drift into two homes for one number. No index, no claim.
    """
    print("\n--- how much of a region partition sits in the informative zone? ---")
    if index_dir is None:
        print("  ⚠ UNMEASURED: no `--index` given, and this file no longer carries a frozen number.")
        print("    Re-derive with `--index INDEX_DIR`, or run `scripts/design/index_census.py INDEX_DIR`")
        print("    for the full census. ⛔ Whatever it says is a property of ONE annotation.")
        return
    import pandas as pd

    feather = index_dir / "regions.feather"
    if not feather.exists():
        raise SystemExit(f"⛔ no region table at {feather} — is {index_dir} a rigel index?")
    length = pd.read_feather(feather)["length"].to_numpy(np.int64)
    n = int(length.size)
    print(f"  index {index_dir}")
    print(f"  regions                       {n:>12,}")
    for cut in (10, 200):
        share = float((length < cut).mean()) * 100.0
        print(f"  regions < {cut:>4} bp                {share:>11.1f} %")
    print(f"  median region length (bp)     {int(np.median(length)):>12,}")
    print(f"  mean region length (bp)       {length.mean():>12,.1f}")
    shorter = float((length < 200).mean())
    verdict = "the MAJORITY" if shorter > 0.5 else "a MINORITY"
    print(f"  -> {verdict} of regions are shorter than one ~200 bp fragment.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--index",
        type=Path,
        default=None,
        help="a rigel index directory; without it the informative-zone section reports UNMEASURED",
    )
    args = ap.parse_args()

    print("=" * 92)
    print("⛔ RETIRED QUESTION — the fragment-length CALIBRATION COMPOSITION channel is deferred past")
    print("   0.8.0 (owner, 2026-08-14). Nothing below may be proposed or ranked as a candidate.")
    print("=" * 92)

    # cfRNA: gDNA is nucleosome-protected (~167 bp mono-nucleosome, tight); RNA differs by protocol.
    show("cfRNA-like, well separated", 167, 15, 90, 35)
    show("cfRNA-like, modest separation", 167, 30, 200, 60)
    show("worst case: identical FL", 180, 50, 180, 50)

    informative_zone(args.index)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
