"""rigel.calibration._diagnostics — typed breakdown of ``n_observed``.

Replaces the legacy ``dict[str, int]`` "annotation gap" diagnostic.
The eight named counters partition every observed unique-mapper
fragment by the conceptual region category it landed in (with no
reference to mask integers in the public surface).

``Diagnostics.total() == payload.n_observed`` is an accountability
sentinel checked by tests.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .scan_payload import (
    MASK_EXON,
    MASK_INTERGENIC,
    MASK_INTRON,
    CalibrationScanPayload,
)


__all__ = ["Diagnostics"]


@dataclass(frozen=True, slots=True)
class Diagnostics:
    """Named decomposition of ``payload.n_observed`` (sums to ``n_observed``)."""

    n_unannotated:            int   # mask 0b000 — no region overlap (decoy/contig)
    n_exon_only:              int   # mask 0b001 — healthy exonic (mostly mRNA)
    n_intron_only:            int   # mask 0b010 — gDNA pool
    n_exon_intron:            int   # mask 0b011 — gDNA pool (exon-intron edge)
    n_intergenic_only:        int   # mask 0b100 — gDNA pool
    n_exon_intergenic:        int   # mask 0b101 — transcription readthrough
    n_intron_intergenic:      int   # mask 0b110 — annotation gap
    n_exon_intron_intergenic: int   # mask 0b111 — annotation gap

    def total(self) -> int:
        return (
            self.n_unannotated
            + self.n_exon_only
            + self.n_intron_only
            + self.n_exon_intron
            + self.n_intergenic_only
            + self.n_exon_intergenic
            + self.n_intron_intergenic
            + self.n_exon_intron_intergenic
        )

    def to_summary_dict(self) -> dict[str, int]:
        return {k: int(v) for k, v in asdict(self).items()}

    @classmethod
    def from_payload(cls, payload: CalibrationScanPayload) -> "Diagnostics":
        gc = payload.global_counts
        return cls(
            n_unannotated            =int(gc[0]),
            n_exon_only              =int(gc[MASK_EXON]),
            n_intron_only            =int(gc[MASK_INTRON]),
            n_exon_intron            =int(gc[MASK_EXON | MASK_INTRON]),
            n_intergenic_only        =int(gc[MASK_INTERGENIC]),
            n_exon_intergenic        =int(gc[MASK_EXON | MASK_INTERGENIC]),
            n_intron_intergenic      =int(gc[MASK_INTRON | MASK_INTERGENIC]),
            n_exon_intron_intergenic =int(gc[MASK_EXON | MASK_INTRON | MASK_INTERGENIC]),
        )
