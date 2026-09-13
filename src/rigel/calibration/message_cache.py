"""The message layer's output, shared across the sweeps that would recompute it identically.

       Gate: ``tests/calibration/test_sweep_backbone.py`` (the cache tests)

`sweep._solve_block` asks the cache for a block's delivered messages before it runs the policy, and
stores them after; `calibrate` builds one cache per run and hands it to every refit sweep. What the layer
reads, and so what the key digests, is the block's :class:`~.messages.BlockContext`, the policy's library
and the policy itself — never the prior.
"""

from __future__ import annotations

import dataclasses
import hashlib

import numpy as np

from .messages import BlockContext, PsiMessage

__all__ = ["MessageCache"]


class MessageCache:
    """The message layer's output, shared across the sweeps that would recompute it identically.

    The message layer — the policy's ``prepare``, the two passes and its ``solve`` — reads only the
    context (observations, geometry, the factory rows, the incoming belief's ``belief_fg`` and the
    own-composition bits ``has_own_composition``), the library and the grid; never the prior.
    ``calibrate`` resets the belief before every sweep, so for the refit sweeps that share a grid every one of those inputs is
    identical and so is every message delivered (measured on the human chain: sweeps 1–3 deliver the
    same ψ rows and cube rows to the bit, and every node hears the same thing). Only the own-claims
    solve and the final ψ see the prior. So a refit sweep's blocks are served from here and pay only
    the two ψ solves.

    CONTENT-KEYED, so it is safe by construction rather than by trust: an entry's key is a digest of
    every input the layer reads, and a changed belief, row, count, library or grid misses. An entry
    holds the delivered rows sparsely (only the non-zero rows), the cube rows as the solve reads them,
    the block's ``held_composition`` and its assertion counts (as a plain dict; the backbone rebuilds its
    `AssertionCounts` from it). Measured on the 18.6M-fragment library: the refit sweeps served in 38 s
    instead of 176 s each, the run 0.65 of its wall, 4.1 GB held with float64 cube rows. Diagnostics never read from it: a captured sweep runs the
    whole layer.
    """

    def __init__(self):
        self._entries: dict = {}
        self.hits = 0
        self.misses = 0

    @staticmethod
    def key(ctx: BlockContext, library, policy) -> bytes:
        """The digest of everything the message layer reads for one block: every field of the block's
        context — iterated from the dataclass, so a field added to the context cannot be left out of the
        digest — then the policy's library and the policy's name and strand model."""
        h = hashlib.blake2b(digest_size=16)
        for f in dataclasses.fields(ctx):
            part = getattr(ctx, f.name)
            if isinstance(part, np.ndarray):
                a = np.ascontiguousarray(part)
                h.update(f"{f.name}{a.dtype.str}{a.shape}".encode())
                h.update(a)
            else:
                h.update(f"{f.name}={part!r}".encode())
        h.update(
            repr(
                (library, getattr(policy, "name", None), getattr(policy, "_strand", None))
            ).encode()
        )
        return h.digest()

    def get(self, key: bytes):
        entry = self._entries.get(key)
        if entry is None:
            self.misses += 1
            return None
        self.hits += 1
        return entry

    def put(self, key: bytes, msg: PsiMessage, held_composition, counts: dict) -> None:
        rows = msg.lam_rows
        if rows is None:
            sparse = None
        else:
            rows = np.asarray(rows)
            idx = np.flatnonzero(np.any(rows != 0.0, axis=1))
            sparse = (rows.shape, idx, rows[idx].copy())
        cube = None if msg.cube_rows is None else {int(k): v for k, v in msg.cube_rows.items()}
        self._entries[key] = (
            sparse,
            cube,
            np.array(held_composition, bool),
            dict(counts),
        )

    @staticmethod
    def message(entry) -> PsiMessage:
        """The stored entry as the message the solve receives — the rows dense again, zeros exact."""
        sparse, cube, _held, _counts = entry
        if sparse is None:
            rows = None
        else:
            shape, idx, kept = sparse
            rows = np.zeros(shape)
            rows[idx] = kept
        return PsiMessage(lam_rows=rows, cube_rows=None if cube is None else dict(cube))

    @property
    def nbytes(self) -> int:
        total = 0
        for sparse, cube, held, _counts in self._entries.values():
            if sparse is not None:
                total += sparse[1].nbytes + sparse[2].nbytes
            if cube:
                total += sum(
                    p.nbytes
                    for v in cube.values()
                    for p in (v.profile_pos, v.profile_neg, v.u)
                    if p is not None
                )
            total += held.nbytes
        return total
