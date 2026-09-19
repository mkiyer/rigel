"""The message layer's output, shared across the sweeps that would recompute it identically.

       Gate: ``tests/calibration/test_sweep_backbone.py`` (the cache tests)

`sweep.solve_chain` asks the cache, per block, for the block's delivery before the kernel runs, hands the
kernel what it holds, and stores what the kernel delivered for the blocks it ran; `calibrate` builds one
cache per run and hands it to every refit sweep. What the layer reads, and so what the key digests, is
the block's slice of the chain's arrays (:class:`~.messages.ChainView`), the incoming belief, the
factory's digest, the library and the policy — never the prior.
"""

from __future__ import annotations

import dataclasses
import hashlib

import numpy as np

from .messages import ChainView

__all__ = ["MessageCache"]


class MessageCache:
    """The message layer's output, shared across the sweeps that would recompute it identically.

    The message layer — the policy's builders, the two passes and its solve, run per block inside the
    kernel — reads only the block's slice of the chain's observations and geometry, the incoming belief
    (the variance freeze of every node's own strand claim), the factory's rows, the grid, the library and
    the policy; never the prior. ``calibrate`` resets the belief before every sweep, so for the refit
    sweeps that share a grid every one of those inputs is identical and so is every message delivered
    (measured on the human chain: sweeps 1–3 deliver the same ψ rows and cube rows to the bit, and every
    node hears the same thing). Only the own-claims solve and the final ψ see the prior. So a refit
    sweep's blocks are served from here and pay only the two ψ solves.

    CONTENT-KEYED, so it is safe by construction rather than by trust: an entry's key is a digest of
    every input the layer reads, and a changed belief, row, count, library or grid misses. The factory
    rows enter by the digest of their inputs (`calibrate.FactoryRows.digest`), which they are a pure
    function of, rather than of their ``(n, K)`` bytes — the same key at 1/K of the hashing; a node's
    own-evidence bit, which the layer also reads, is a function of what the key already holds (a count
    above zero, the protocol's liveness, the single-strand bits, the factory row's curvature) and needs no
    field of its own. An entry is the block's delivery exactly as the kernel returns it and takes it back
    — whether rows were delivered, the delivered rows with their slots (sparse by construction: the solve
    knows which rows it wrote), the cube table or ``None``, the owned slots' held-composition bits.
    Diagnostics never read from it: a captured sweep runs the whole layer.
    """

    def __init__(self):
        self._entries: dict = {}
        self.hits = 0
        self.misses = 0

    @staticmethod
    def key(
        view: ChainView, block, belief_fg, factory_digest: bytes | None, library, policy
    ) -> bytes:
        """The digest of everything the message layer reads for one block: the block's slice of every
        array field of the chain view — iterated from the dataclass, so a field added to the view cannot
        be left out of the digest — and its scalars, the block's slice of the incoming belief, the
        factory's digest for the block (or its absence), then the policy's library, name and strand
        model."""
        h = hashlib.blake2b(digest_size=16)
        sl = slice(int(block.start), int(block.end))

        def array(name, a):
            a = np.ascontiguousarray(a)
            h.update(f"{name}{a.dtype.str}{a.shape}".encode())
            h.update(a)

        for f in dataclasses.fields(view):
            part = getattr(view, f.name)
            if isinstance(part, np.ndarray):
                array(f.name, part[sl])
            else:
                h.update(f"{f.name}={part!r}".encode())
        array("belief_fg", np.asarray(belief_fg, np.float64)[sl])
        h.update(b"factory:" + (b"none" if factory_digest is None else factory_digest))
        h.update(repr((library, policy.name, policy.strand)).encode())
        return h.digest()

    def served(self, key: bytes):
        """The block's delivery under ``key`` — as the kernel takes it — or ``None``; counted as a hit or a
        miss."""
        entry = self._entries.get(key)
        if entry is None:
            self.misses += 1
            return None
        self.hits += 1
        return entry

    def put(self, key: bytes, delivery) -> None:
        """The kernel's delivery for a block it ran the layer on: ``(rows_delivered, slot, rows, cube or
        None, held)``."""
        self._entries[key] = delivery

    @property
    def nbytes(self) -> int:
        total = 0
        for _live, slot, rows, cube, held in self._entries.values():
            total += slot.nbytes + rows.nbytes + held.nbytes
            if cube is not None:
                total += sum(a.nbytes for a in cube)
        return total
