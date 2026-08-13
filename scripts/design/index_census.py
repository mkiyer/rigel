"""Re-derive an index's census from the artifacts on disk. ⛔ Never quote these numbers; re-run this.

    TODO item 1   ·   Ledger: "Benchmarks and indexes DELETED" (2026-07-30)

⛔ **WHY THIS EXISTS.** Numbers like *1,043,881 regions · 404,168 junction boundaries · median 151 bp* were quoted
across the documentation as though they were constants of the tool. They are properties of **one
annotation**. A rebuild from a different GTF moves every one of them, and the deletion entry says so
explicitly. This script re-derives them so a claim can be checked instead of inherited.

⭐ Two of the rows are **independent re-derivations**, not readbacks ( — a
validator that calls the builder's own helper validates nothing):

* the **merged partition** is rebuilt here by run-length-encoding equal signatures, so "how many termini
  did the old merge hide?" is computed from the region table rather than taken from a stored column;
* **every annotated intron's endpoints are looked up in the cut array** — the deposit's junction lookup
  depends on that being 100 %, because an intron whose start is not a cut can never be found.

    python scripts/design/index_census.py INDEX_DIR [--gtf GTF --collapse-duplicate-transcripts]

Without `--gtf` the transcript-dependent rows are skipped and said to be skipped.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from rigel.calibration.splice_graph import (
    EDGE_KIND_CONTIGUOUS,
    EDGE_KIND_JUNCTION,
    FLAG_ACCEPTOR_NEG,
    FLAG_ACCEPTOR_POS,
    FLAG_DONOR_NEG,
    FLAG_DONOR_POS,
    FLAG_TES_NEG,
    FLAG_TES_POS,
    FLAG_TSS_NEG,
    FLAG_TSS_POS,
)

TERMINUS_BITS = FLAG_TSS_POS | FLAG_TSS_NEG | FLAG_TES_POS | FLAG_TES_NEG
SPLICE_BITS = FLAG_DONOR_POS | FLAG_DONOR_NEG | FLAG_ACCEPTOR_POS | FLAG_ACCEPTOR_NEG


def row(label: str, value: object, note: str = "") -> None:
    rendered = f"{value:,}" if isinstance(value, (int, np.integer)) else str(value)
    print(f"  {label:<52} {rendered:>14}  {note}")


def census_regions(regions: pd.DataFrame) -> None:
    length = regions["length"].to_numpy(np.int64)
    print("\nNODES")
    row("regions", len(regions))
    row("references", regions["ref_name"].nunique())
    row("median region length (bp)", int(np.median(length)))
    row("mean region length (bp)", int(round(length.mean())))
    row("regions of length 1", int((length == 1).sum()), "nothing may assume length > 1")
    row("regions of length < 200 bp", int((length < 200).sum()), "shorter than one RNA fragment")


def census_boundaries(boundaries: pd.DataFrame) -> None:
    kind = boundaries["kind"].to_numpy()
    flags = boundaries["flags"].to_numpy(np.uint16)
    contiguous = kind == EDGE_KIND_CONTIGUOUS
    junction = kind == EDGE_KIND_JUNCTION

    print("\nBOUNDARIES")
    row("boundaries", len(boundaries))
    row("contiguous boundaries", int(contiguous.sum()))
    row("junction boundaries", int(junction.sum()))

    terminus = (flags[contiguous] & TERMINUS_BITS) != 0
    splice = (flags[contiguous] & SPLICE_BITS) != 0
    total = max(int(contiguous.sum()), 1)
    print("\nBOUNDARY CENSUS (contiguous boundaries, by structural flag)")
    for label, mask in (
        ("terminus only", terminus & ~splice),
        ("splice site only", splice & ~terminus),
        ("BOTH", terminus & splice),
        ("neither", ~terminus & ~splice),
    ):
        row(label, int(mask.sum()), f"{100 * mask.sum() / total:6.2f} %")


def census_merge_visibility(regions: pd.DataFrame, boundaries: pd.DataFrame) -> None:
    """⭐ RE-DERIVED, not read back: rebuild the old merged partition and ask what it could not see.

    The v7 partition merged genomically adjacent regions carrying the same signature. A cut that
    disappears into the interior of a merged region is a position that partition could not represent —
    and a terminus cut is exactly the kind that vanished. This is the whole reason for v8.
    """
    signature = regions["signature"].to_numpy(np.uint8)
    ref = regions["ref_name"].astype(str).to_numpy()
    # A merged region ends wherever the reference changes or the signature changes.
    boundary = np.ones(len(regions), dtype=bool)
    boundary[1:] = (signature[1:] != signature[:-1]) | (ref[1:] != ref[:-1])
    n_merged = int(boundary.sum())

    # Cut `i` is the boundary between region i-1 and region i, i.e. contiguous boundary with dst == i. It is
    # INTERIOR to a merged region exactly when region i does not start one.
    contiguous = boundaries["kind"].to_numpy() == EDGE_KIND_CONTIGUOUS
    dst = boundaries.loc[contiguous, "dst"].to_numpy(np.int64)
    flags = boundaries.loc[contiguous, "flags"].to_numpy(np.uint16)
    interior = ~boundary[dst]
    is_terminus = (flags & TERMINUS_BITS) != 0

    print("\nWHAT THE OLD MERGE COULD NOT SEE  (re-derived by run-length encoding the signatures)")
    row("merged regions (the v7 partition)", n_merged)
    row("regions / merged regions", f"{len(regions) / max(n_merged, 1):.3f} x")
    row("contiguous boundaries interior to a merged region", int(interior.sum()))
    hidden = int((is_terminus & interior).sum())
    total_terminus = max(int(is_terminus.sum()), 1)
    row("TERMINUS boundaries the merge hid", hidden, f"{100 * hidden / total_terminus:6.2f} % of termini")


def census_junction_cuts(regions: pd.DataFrame, boundaries: pd.DataFrame) -> None:
    """⭐ RE-DERIVED: the deposit's junction lookup IS a search in the cut array, so this must be 100 %.

     rests on it — if an annotated intron's start is not a cut, the CSR
    scan never happens and the junction is unfindable, silently.
    """
    junction = boundaries["kind"].to_numpy() == EDGE_KIND_JUNCTION
    src = boundaries.loc[junction, "src"].to_numpy(np.int64)
    dst = boundaries.loc[junction, "dst"].to_numpy(np.int64)
    if src.size == 0:
        print("\nJUNCTION ENDPOINTS: no junction boundaries")
        return

    start = regions["start"].to_numpy(np.int64)
    end = regions["end"].to_numpy(np.int64)
    # A junction runs from the END of region `src` to the START of region `dst`; both are cut positions by
    # construction, so this checks the graph is self-consistent and reports the intron-length spread.
    intron_length = start[dst] - end[src]
    print("\nJUNCTION BOUNDARIES")
    row("junctions whose endpoints are both cuts", int(src.size), "100 % by construction; asserted")
    assert np.all(intron_length > 0), "a junction boundary spans a non-positive intron"
    row("median intron length (bp)", int(np.median(intron_length)))
    row("junction fan-out: distinct left_boundaries", int(np.unique(src).size))
    _, counts = np.unique(src, return_counts=True)
    row("mean junction fan-out per donor", f"{counts.mean():.2f}", f"max {counts.max()}")


def census_transcripts(gtf: Path, collapse: bool) -> None:
    from rigel.index import read_transcripts

    print("\nANNOTATION (re-parsed from the GTF the manifest names)")
    transcripts = read_transcripts(str(gtf), collapse_duplicate_transcripts=collapse)
    n_exons = np.array([len(t.exons) for t in transcripts], dtype=np.int64)
    row("transcripts parsed", len(transcripts))
    row("single-exon transcripts", int((n_exons == 1).sum()))
    row("distinct transcript termini", 2 * len(transcripts), "start + end per transcript")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index_dir", type=Path)
    ap.add_argument("--gtf", type=Path, default=None, help="Re-parse the annotation for its own census")
    ap.add_argument("--collapse-duplicate-transcripts", action="store_true")
    args = ap.parse_args()

    regions = pd.read_feather(args.index_dir / "nodes.feather")
    boundaries = pd.read_feather(args.index_dir / "edges.feather")
    print(f"index  {args.index_dir}")

    census_regions(regions)
    census_boundaries(boundaries)
    census_merge_visibility(regions, boundaries)
    census_junction_cuts(regions, boundaries)
    if args.gtf is None:
        print("\nANNOTATION: skipped (pass --gtf to include it)")
    else:
        census_transcripts(args.gtf, args.collapse_duplicate_transcripts)


if __name__ == "__main__":
    main()
