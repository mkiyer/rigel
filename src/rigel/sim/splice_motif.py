"""Canonical splice-site motifs for the simulator: :func:`splice_donor_acceptor` gives the dinucleotides and
:func:`place_intron_motif` writes them into a genome."""

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
    """Place one intron's ``(donor, acceptor)`` dinucleotides through ``edit(pos, bases)``, which overwrites
    ``len(bases)`` bases at 0-based ``pos``: the donor at ``intron_start``, the acceptor at ``intron_end - 2``.
    Raises on an intron shorter than 4 bp.
    """
    if intron_end - intron_start < 4:
        raise ValueError(
            f"Intron ({intron_start},{intron_end}) too short for splice motifs "
            f"(length {intron_end - intron_start} < 4)"
        )
    donor, acceptor = splice_donor_acceptor(strand)
    edit(intron_start, donor)
    edit(intron_end - 2, acceptor)
