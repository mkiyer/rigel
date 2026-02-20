"""
hulkrna.sim.annotation — Gene/transcript builder with splice-motif injection.

Provides a ``GeneBuilder`` that accumulates gene and transcript
definitions, injects canonical splice-site motifs (GT-AG for +strand,
CT-AC for −strand) into the ``MutableGenome``, and writes
GENCODE-style GTF annotation files.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

from ..gtf import GTF
from ..transcript import Transcript
from ..types import Interval, Strand

from .genome import MutableGenome

logger = logging.getLogger(__name__)

__all__ = ["GeneBuilder"]


@dataclass
class _TranscriptSpec:
    """Internal specification for a transcript before it becomes a Transcript."""

    t_id: str
    g_id: str
    g_name: str
    g_type: str
    strand: Strand
    ref: str
    exons: list[Interval]
    abundance: float = 100.0
    nrna_abundance: float = 0.0


class GeneBuilder:
    """Accumulates gene/transcript annotations and edits the genome.

    For each intron (gap between consecutive exons), the builder injects
    canonical splice-site motifs into the ``MutableGenome``:

    - **Positive strand** (FR convention):
      ``genome[intron_start:intron_start+2] = 'GT'`` (5′ donor)
      ``genome[intron_end-2:intron_end] = 'AG'`` (3′ acceptor)

    - **Negative strand** (reverse complement on genomic + strand):
      ``genome[intron_start:intron_start+2] = 'CT'`` (rev-comp acceptor)
      ``genome[intron_end-2:intron_end] = 'AC'`` (rev-comp donor)

    Parameters
    ----------
    genome : MutableGenome
        The genome to edit with splice motifs.
    ref_name : str or None
        Chromosome name. Defaults to ``genome.name``.
    """

    def __init__(self, genome: MutableGenome, ref_name: str | None = None):
        self.genome = genome
        self.ref_name = ref_name or genome.name
        self._specs: list[_TranscriptSpec] = []

    # -- Gene addition --------------------------------------------------------

    def add_gene(
        self,
        gene_id: str,
        strand: str | Strand,
        transcripts: list[dict],
        *,
        gene_name: str | None = None,
        gene_type: str = "protein_coding",
    ) -> None:
        """Add a gene with one or more transcript isoforms.

        Parameters
        ----------
        gene_id : str
            Unique gene identifier (e.g. ``"g1"``).
        strand : str or Strand
            ``"+"`` / ``"-"`` or a ``Strand`` enum value.
        transcripts : list[dict]
            Each dict has keys:

            - ``"t_id"`` (str): transcript ID
            - ``"exons"`` (list[tuple[int,int]]): 0-based half-open exon coords
            - ``"abundance"`` (float, optional): mature mRNA molecular
              abundance, default 100.  Proportional to molecule count
              after length normalisation (TPM-like).
            - ``"nrna_abundance"`` (float, optional): nascent RNA molecular
              abundance, default 0.  **Independent** of ``"abundance"``;
              the two values are NOT coupled by subtraction.  Read
              fractions depend on ``abundance × spliced_eff_len`` vs.
              ``nrna_abundance × genomic_span_eff_len``.
        gene_name : str or None
            Gene symbol. Defaults to *gene_id*.
        gene_type : str
            Biotype (default ``"protein_coding"``).
        """
        if isinstance(strand, str):
            strand = Strand.from_str(strand)
        if gene_name is None:
            gene_name = gene_id

        for tdef in transcripts:
            t_id = tdef["t_id"]
            raw_exons = tdef["exons"]
            abundance = tdef.get("abundance", 100.0)
            nrna_abundance = tdef.get("nrna_abundance", 0.0)

            # Build sorted Interval list
            exons = sorted(Interval(s, e) for s, e in raw_exons)
            self._validate_exons(t_id, exons)

            # Inject splice-site motifs for each intron
            for i in range(1, len(exons)):
                intron_start = exons[i - 1].end
                intron_end = exons[i].start
                self._inject_splice_motif(intron_start, intron_end, strand)

            self._specs.append(
                _TranscriptSpec(
                    t_id=t_id,
                    g_id=gene_id,
                    g_name=gene_name,
                    g_type=gene_type,
                    strand=strand,
                    ref=self.ref_name,
                    exons=exons,
                    abundance=abundance,
                    nrna_abundance=nrna_abundance,
                )
            )

    # -- Validation -----------------------------------------------------------

    def _validate_exons(self, t_id: str, exons: list[Interval]) -> None:
        """Check exons are non-overlapping and within genome bounds."""
        genome_len = len(self.genome)
        for i, e in enumerate(exons):
            if e.start < 0 or e.end > genome_len:
                raise ValueError(
                    f"Transcript {t_id} exon {i} ({e.start},{e.end}) "
                    f"outside genome bounds [0, {genome_len})"
                )
            if e.start >= e.end:
                raise ValueError(
                    f"Transcript {t_id} exon {i} has start >= end: "
                    f"({e.start},{e.end})"
                )
        for i in range(1, len(exons)):
            if exons[i].start < exons[i - 1].end:
                raise ValueError(
                    f"Transcript {t_id} exons {i-1} and {i} overlap: "
                    f"({exons[i-1].start},{exons[i-1].end}) vs "
                    f"({exons[i].start},{exons[i].end})"
                )

    # -- Splice motif injection -----------------------------------------------

    def _inject_splice_motif(
        self, intron_start: int, intron_end: int, strand: Strand
    ) -> None:
        """Edit the genome to place canonical splice-site dinucleotides.

        For a positive-strand intron (GT-AG):
            genome[intron_start:intron_start+2] = 'GT'  (5' donor)
            genome[intron_end-2:intron_end]     = 'AG'  (3' acceptor)

        For a negative-strand intron (CT-AC on genomic + strand):
            genome[intron_start:intron_start+2] = 'CT'
            genome[intron_end-2:intron_end]     = 'AC'
        """
        intron_len = intron_end - intron_start
        if intron_len < 4:
            raise ValueError(
                f"Intron ({intron_start},{intron_end}) too short for "
                f"splice motifs (length {intron_len} < 4)"
            )
        if strand == Strand.POS:
            self.genome.edit(intron_start, "GT")
            self.genome.edit(intron_end - 2, "AG")
        elif strand == Strand.NEG:
            self.genome.edit(intron_start, "CT")
            self.genome.edit(intron_end - 2, "AC")
        else:
            raise ValueError(f"Cannot inject splice motifs for strand {strand}")

    # -- Output ---------------------------------------------------------------

    def get_transcripts(self) -> list[Transcript]:
        """Build ``Transcript`` objects from accumulated specs.

        Assigns sequential ``t_index`` and unique ``g_index`` per gene.
        Returns transcripts sorted by ``(ref, start, end, strand)``.
        """
        g_id_to_index: dict[str, int] = {}
        transcripts: list[Transcript] = []

        for spec in self._specs:
            if spec.g_id not in g_id_to_index:
                g_id_to_index[spec.g_id] = len(g_id_to_index)

            t = Transcript(
                ref=spec.ref,
                strand=spec.strand,
                exons=list(spec.exons),
                t_id=spec.t_id,
                g_id=spec.g_id,
                g_name=spec.g_name,
                g_type=spec.g_type,
                g_index=g_id_to_index[spec.g_id],
                is_basic=True,
                abundance=spec.abundance,
                nrna_abundance=spec.nrna_abundance,
            )
            t.compute_length()
            transcripts.append(t)

        # Sort by genomic position
        transcripts.sort(key=lambda t: (t.ref, t.start, t.end, t.strand))

        # Assign sequential t_index after sorting
        for ti, t in enumerate(transcripts):
            t.t_index = ti

        return transcripts

    def write_gtf(self, output_dir: Path) -> Path:
        """Write all transcripts as a GENCODE-style exon GTF file.

        Parameters
        ----------
        output_dir : Path
            Directory to write ``annotations.gtf`` into.

        Returns
        -------
        Path
            Path to the GTF file.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        gtf_path = output_dir / "annotations.gtf"

        transcripts = self.get_transcripts()

        with open(gtf_path, "w") as f:
            for t in transcripts:
                strand_str = Strand.to_str(t.strand)
                for exon in t.exons:
                    gtf_obj = GTF(
                        seqname=t.ref,
                        source="hulkrna_sim",
                        feature="exon",
                        start=exon.start,  # GTF.__str__ adds +1
                        end=exon.end,
                        score=t.abundance,
                        strand=strand_str,
                        phase=".",
                        attrs={
                            "gene_id": t.g_id,
                            "transcript_id": t.t_id,
                            "gene_name": t.g_name,
                            "gene_type": t.g_type,
                        },
                        tags={"basic"},
                    )
                    f.write(str(gtf_obj) + "\n")

        logger.info(
            f"Wrote {len(transcripts)} transcripts "
            f"({sum(len(t.exons) for t in transcripts)} exon lines) → {gtf_path}"
        )
        return gtf_path
