"""Canonical splice-site motif convention for the simulator — single source of truth.

Two call sites inject splice dinucleotides into a synthetic genome with different substrates and
error policies (``synthetic_genome.inject_splice_sites`` edits a numpy byte array and skips
too-short introns; ``annotation.GeneBuilder._inject_splice_motif`` edits a ``MutableGenome`` and
raises) — but they must agree on *which* dinucleotides define a donor/acceptor. That shared
knowledge lives here so the convention can never drift between them.
"""

from __future__ import annotations

from ..types import Strand


def splice_donor_acceptor(strand: Strand | str) -> tuple[str, str]:
    """Return the genomic-coordinate ``(donor, acceptor)`` dinucleotides for an intron.

    * ``+`` strand (canonical GT-AG): ``GT`` at the intron 5′ end, ``AG`` at the 3′ end.
    * ``-`` strand: ``CT`` / ``AC`` — the reverse complement of GT-AG, as written on the ``+``
      genomic strand.

    Accepts a :class:`Strand` enum or a ``"+"`` / ``"-"`` string (the two call sites differ).
    """
    if strand == Strand.POS or strand == "+":
        return "GT", "AG"
    if strand == Strand.NEG or strand == "-":
        return "CT", "AC"
    raise ValueError(f"cannot inject splice motifs for strand {strand!r}")
