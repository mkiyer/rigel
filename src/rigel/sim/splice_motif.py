"""Canonical splice-site motif injection for the simulator — single source of truth.

Two call sites inject splice dinucleotides into a synthetic genome with different substrates and
error policies (``synthetic_genome.inject_splice_sites`` edits a numpy byte array and *skips*
too-short introns; ``annotation.GeneBuilder._inject_splice_motif`` edits a ``MutableGenome`` and
*raises*). Both the dinucleotide *convention* (:func:`splice_donor_acceptor`) and the *placement*
(:func:`place_intron_motif`) live here so neither can drift between them; the call sites differ
only in the edit callback (their substrate) and the ``on_short`` policy.
"""

from __future__ import annotations

from collections.abc import Callable

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


def place_intron_motif(
    edit: Callable[[int, str], None],
    intron_start: int,
    intron_end: int,
    strand: Strand | str,
    *,
    on_short: str = "raise",
) -> None:
    """Place the canonical ``(donor, acceptor)`` dinucleotides for ONE intron via ``edit(pos, bases)``.

    ``edit(pos, bases)`` overwrites ``len(bases)`` bases at 0-based ``pos`` in the caller's substrate
    (a ``MutableGenome.edit``, or a numpy-array setter). The donor lands at ``intron_start`` and the
    acceptor at ``intron_end - 2`` (genomic-coordinate convention, :func:`splice_donor_acceptor`).

    ``on_short`` selects the policy for an intron shorter than 4 bp (no room for both dinucleotides):
    ``"raise"`` (the ``GeneBuilder`` contract — a construction error) or ``"skip"`` (the
    ``inject_splice_sites`` contract — leave it unmodified). Generated references never emit such
    introns, so the two policies agree in practice.
    """
    if intron_end - intron_start < 4:
        if on_short == "skip":
            return
        raise ValueError(
            f"Intron ({intron_start},{intron_end}) too short for splice motifs "
            f"(length {intron_end - intron_start} < 4)"
        )
    donor, acceptor = splice_donor_acceptor(strand)
    edit(intron_start, donor)
    edit(intron_end - 2, acceptor)
