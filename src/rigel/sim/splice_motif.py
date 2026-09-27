"""Canonical splice-site motif injection for the simulator.

:func:`splice_donor_acceptor` gives the dinucleotide convention and :func:`place_intron_motif` places it;
the caller is ``annotation.GeneBuilder._inject_splice_motif``.
"""

from __future__ import annotations

from collections.abc import Callable

from ..types import Strand


def splice_donor_acceptor(strand: Strand) -> tuple[str, str]:
    """Return the genomic-coordinate ``(donor, acceptor)`` dinucleotides for an intron.

    * ``+`` strand (canonical GT-AG): ``GT`` at the intron 5′ end, ``AG`` at the 3′ end.
    * ``-`` strand: ``CT`` / ``AC`` — the reverse complement of GT-AG, as written on the ``+``
      genomic strand.
    """
    if strand == Strand.POS:
        return "GT", "AG"
    if strand == Strand.NEG:
        return "CT", "AC"
    raise ValueError(f"cannot inject splice motifs for strand {strand!r}")


def place_intron_motif(
    edit: Callable[[int, str], None],
    intron_start: int,
    intron_end: int,
    strand: Strand,
) -> None:
    """Place the canonical ``(donor, acceptor)`` dinucleotides for ONE intron via ``edit(pos, bases)``.

    ``edit(pos, bases)`` overwrites ``len(bases)`` bases at 0-based ``pos`` in the caller's substrate.
    The donor lands at ``intron_start`` and the acceptor at ``intron_end - 2`` (genomic-coordinate
    convention, :func:`splice_donor_acceptor`). Raises on an intron shorter than 4 bp (no room for both
    dinucleotides).
    """
    if intron_end - intron_start < 4:
        raise ValueError(
            f"Intron ({intron_start},{intron_end}) too short for splice motifs "
            f"(length {intron_end - intron_start} < 4)"
        )
    donor, acceptor = splice_donor_acceptor(strand)
    edit(intron_start, donor)
    edit(intron_end - 2, acceptor)
