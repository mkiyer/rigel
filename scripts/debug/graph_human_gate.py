"""HUMAN-SCALE GATE for the splice graph — I1-I12 (incl. I3b and the I11 walk), the census, §8 budgets.

    python scripts/debug/graph_human_gate.py [INDEX_DIR]

Reconstructs Transcript objects from the index feathers (exact — the feathers are derived from the
GTF) so the builder can be exercised on the real annotation without re-parsing GENCODE.

⚠ Its P2'/P3 half is GONE with the v7 partition (accumulator plan W1b-clean). The independent check
on the signature that P3 used to provide is now **I3b**, inside `validate_graph`: each node's
signature is recomputed from its midpoint by direct interval containment. It needs no v7, and it runs
below as part of the gate.
"""

from __future__ import annotations

import resource
import sys
import time

import numpy as np
import pandas as pd

from rigel.calibration.splice_graph import EDGE_KIND_JUNCTION, build_splice_graph, validate_graph
from rigel.index import load_ref_lengths
from rigel.transcript import Interval, Transcript

D = sys.argv[1] if len(sys.argv) > 1 else "/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index_v7"


def transcripts_from_index(d):
    tx = pd.read_feather(f"{d}/transcripts.feather")
    iv = pd.read_feather(f"{d}/intervals.feather")
    ex = iv[(iv["interval_type"] == 0) & (iv["t_index"] >= 0)].sort_values(
        ["t_index", "start"], kind="stable"
    )
    by_t: dict[int, list] = {}
    for t, s, e in zip(
        ex["t_index"].to_numpy(np.int64),
        ex["start"].to_numpy(np.int64),
        ex["end"].to_numpy(np.int64),
    ):
        by_t.setdefault(int(t), []).append(Interval(int(s), int(e)))
    return [
        Transcript(
            ref=str(row.ref),
            strand=row.strand,
            exons=by_t.get(t, []),
            t_id=row.t_id,
            t_index=t,
            is_nrna=bool(row.is_nrna),
            is_synthetic=bool(row.is_synthetic),
        )
        for t, row in enumerate(tx.itertuples(index=False))
    ]


t0 = time.time()
txs = transcripts_from_index(D)
reflen = load_ref_lengths(f"{D}/ref_lengths.feather")
print(f"reconstructed {len(txs):,} transcripts in {time.time() - t0:.1f}s")

t0 = time.time()
nodes, edges = build_splice_graph(txs, reflen)
build_s = time.time() - t0
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
print(f"\n⭐ BUILD: {build_s:.1f}s, peak RSS {rss:.2f} GB   (§8 budget: <10s, <1GB)")

t0 = time.time()
validate_graph(nodes, edges, reflen, transcripts=txs)
print(f"⭐ validate_graph I1-I13 (incl. I3b, I13 + the I11 walk): {time.time() - t0:.1f}s   (§8: <5s)")

j = edges["kind"].to_numpy() == EDGE_KIND_JUNCTION
L = nodes["length"].to_numpy()
print("\n§3.4 CENSUS")
print(f"  nodes            {len(nodes):,}")
print(
    f"  median / p25/p75 {int(np.median(L))} / {int(np.percentile(L, 25))} / "
    f"{int(np.percentile(L, 75))} bp"
)
print(f"  1 bp nodes       {int((L == 1).sum()):,}")
print(f"  <10bp / <200bp   {100 * (L < 10).mean():.1f}% / {100 * (L < 200).mean():.1f}%")
print(f"  contiguous edges {int((~j).sum()):,}")
print(f"  junction edges   {int(j.sum()):,}")
