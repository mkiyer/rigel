"""The block plumbing of the locus solve — how one block of the chain is cut out of it and how its
results are put back.

       Gate: ``tests/calibration/test_sweep_backbone.py`` (the block tests)

The chain is solved a LOCUS BLOCK at a time (`sweep.solve_chain`, `region_chain.locus_blocks`). Three
helpers, and nothing about what a solve or a message is: `view_fields` is every per-slot array a policy
may read, for the whole chain or one block; `block_slice` restricts a chain, its statics, its geometry or
a belief to one block's slots; `gather` re-assembles the blocks' diagnostic captures into the chain's.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from .messages import Received
from .region_chain import REGION, RegionChain
from .simplex_logodds import CompositionPriors

__all__ = ["block_slice", "gather", "view_fields"]


def view_fields(chain, statics, geometry, structure) -> dict:
    """Every per-slot array a policy may read, under the two belief-free headings of
    :class:`~.messages.ChainView` — for the whole chain (the library) or one block (the context)."""
    return dict(
        # observations
        eff_gdna=np.asarray(geometry.eff_gdna, np.float64),
        eff_rna=np.asarray(geometry.eff_rna, np.float64),
        sj_count=np.asarray(geometry.sj_count, np.float64),  # [n, 2] by TRANSCRIPT strand
        sj_count_lo=np.asarray(geometry.sj_count_lo, np.float64),
        sj_count_hi=np.asarray(geometry.sj_count_hi, np.float64),
        route_rate_lo=np.asarray(geometry.route_rate_lo, np.float64),
        route_rate_hi=np.asarray(geometry.route_rate_hi, np.float64),
        unspliced_count=np.asarray(geometry.unspliced_count, np.float64),  # [n, 2] by GENOME strand
        spliced_count=np.asarray(geometry.spliced_count, np.float64),
        # geometry / structure
        left=np.asarray(chain.left, np.int64),
        right=np.asarray(chain.right, np.int64),
        is_boundary=np.asarray(chain.kind) != REGION,
        is_exon_region=np.asarray(structure.is_exon_region, bool),
        free_pos=np.asarray(statics.free_pos, bool),
        free_neg=np.asarray(statics.free_neg, bool),
        exon_pos=np.asarray(structure.exon_pos, bool),
        exon_neg=np.asarray(structure.exon_neg, bool),
        boundary_flags=statics.boundary_flags,
    )


def block_slice(obj, sl: slice):
    """``obj`` — a chain, its statics, its geometry or a belief — restricted to the slots ``sl``: every
    per-slot array sliced (a view, no copy), ``n_slots`` updated, and a chain's links re-based so that
    a neighbour outside the block is no neighbour (``-1``). The block's first slot is a terminal or a
    reference start and its last read slot a terminal, so a link that leaves the block always leaves
    it at a node that receives nothing, and what that node holds is not read by any consumer."""
    fields = {f.name: getattr(obj, f.name) for f in dataclasses.fields(obj)}
    arrays = {k: v for k, v in fields.items() if isinstance(v, np.ndarray) and v.ndim >= 1}
    n = next(iter(arrays.values())).shape[0]
    new = {k: v[sl] for k, v in arrays.items() if v.shape[0] == n}
    if "n_slots" in fields:
        new["n_slots"] = len(range(*sl.indices(n)))
    if isinstance(obj, RegionChain):
        m = new["kind"].shape[0]
        for k in ("left", "right"):
            a = new[k] - sl.start
            new[k] = np.where((a >= 0) & (a < m), a, -1)
    return dataclasses.replace(obj, **new)


def gather(diagnostics: list, n: int) -> dict:
    """The blocks' diagnostic captures as the chain's: per-slot arrays concatenated over each block's
    OWNED slots, ψ's prior arms likewise, the delivered rows zero-filled where a block delivered
    nothing (``None`` only when no block delivered), the two tables and the cube rows re-keyed to the
    chain; the scalars and grids, identical in every block, taken once."""
    if not diagnostics:
        return {}
    out: dict = {}
    first = diagnostics[0][1]
    for key in first:
        vals = [(b, c[key]) for b, c in diagnostics]
        owned = [(b.stop - b.start, v) for b, v in vals]
        if key in ("from_left", "from_right"):
            out[key] = Received.concat([v.take(k) for k, v in owned])
        elif key == "cube_rows":
            cube: dict = {}
            for b, v in vals:
                for slot, row in (v or {}).items():
                    if int(slot) < b.stop - b.start:
                        cube[int(slot) + b.start] = row
            out[key] = cube or None
        elif isinstance(first[key], CompositionPriors):
            out[key] = CompositionPriors(
                *(
                    None
                    if getattr(first[key], m) is None
                    else np.concatenate([getattr(v, m)[:k] for k, v in owned], axis=0)
                    for m in ("gdna", "rna")
                )
            )
        elif key == "lam_rows":
            if all(v is None for _k, v in owned):
                out[key] = None
            else:
                K = next(np.asarray(v).shape[1] for _k, v in owned if v is not None)
                out[key] = np.concatenate(
                    [np.zeros((k, K)) if v is None else np.asarray(v)[:k] for k, v in owned], axis=0
                )
        elif isinstance(first[key], np.ndarray) and first[key].ndim >= 1:
            out[key] = np.concatenate([np.asarray(v)[:k] for k, v in owned], axis=0)
        else:
            out[key] = first[key]
    return out
