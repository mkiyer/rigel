"""``GeneBuilder``: accumulates gene and transcript definitions, writes their canonical splice-site motifs into a
``MutableGenome``, and writes a GENCODE-style GTF."""

import logging
from dataclasses import dataclass
from pathlib import Path

from ..gtf import GTFRecord
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


#: Every key a transcript entry of ``add_gene`` reads, and the ones it must have.
TRANSCRIPT_KEYS = frozenset({"t_id", "exons", "abundance", "nrna_abundance"})
_TRANSCRIPT_REQUIRED = frozenset({"t_id", "exons"})


def check_transcript_keys(gene_id: str, tdef: dict) -> None:
    """Refuse a transcript entry with a key ``add_gene`` does not read (a misspelt ``abundance`` would simulate the
    default) or without one it needs."""
    unknown = sorted(set(tdef) - TRANSCRIPT_KEYS)
    missing = sorted(_TRANSCRIPT_REQUIRED - set(tdef))
    if unknown or missing:
        raise ValueError(
            f"gene {gene_id!r}, transcript {tdef.get('t_id')!r}: "
            + (f"unknown key(s) {unknown}" if unknown else f"missing {missing}")
            + f"; a transcript takes {sorted(TRANSCRIPT_KEYS)}"
        )


def write_splice_motif(
    genome: MutableGenome, intron_start: int, intron_end: int, strand: Strand
) -> None:
    """Write one intron's canonical splice-site dinucleotides into *genome*, as read on the plus strand: ``GT`` at
    the intron's 5' end and ``AG`` at its 3' end on ``+``, their reverse complements ``CT`` / ``AC`` on ``-``. Raises
    on an intron shorter than 4 bp, or a strand that is neither."""
    if intron_end - intron_start < 4:
        raise ValueError(
            f"Intron ({intron_start},{intron_end}) too short for splice motifs "
            f"(length {intron_end - intron_start} < 4)"
        )
    if strand == Strand.POS:
        low_motif, high_motif = "GT", "AG"
    elif strand == Strand.NEG:
        low_motif, high_motif = "CT", "AC"
    else:
        raise ValueError(f"cannot inject splice motifs for strand {strand!r}")
    genome.edit(intron_start, low_motif)
    genome.edit(intron_end - 2, high_motif)


class GeneBuilder:
    """Accumulates gene and transcript annotations, writing each intron's splice-site motifs into the genome
    (GT…AG on the plus strand, CT…AC on the minus strand, as read on the plus strand).

    Parameters
    ----------
    genome : MutableGenome
        The genome to edit with splice motifs.
    ref_name : str or None
        Reference name. Defaults to ``genome.name``.
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
            - ``"abundance"`` (float, optional): the transcript's molecular
              abundance, default 100 (TPM-like).
            - ``"nrna_abundance"`` (float, optional): its nascent molecular
              abundance, default 0, independent of ``"abundance"``. Read
              fractions follow ``abundance × spliced_eff_len`` against
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
            check_transcript_keys(gene_id, tdef)
            t_id = tdef["t_id"]
            raw_exons = tdef["exons"]
            abundance = tdef.get("abundance", _TranscriptSpec.abundance)
            nrna_abundance = tdef.get("nrna_abundance", _TranscriptSpec.nrna_abundance)

            exons = sorted(Interval(s, e) for s, e in raw_exons)
            self._validate_exons(t_id, exons)

            for left, right in zip(exons, exons[1:]):
                write_splice_motif(self.genome, left.end, right.start, strand)

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
                    f"Transcript {t_id} exon {i} has start >= end: ({e.start},{e.end})"
                )
        for i in range(1, len(exons)):
            if exons[i].start < exons[i - 1].end:
                raise ValueError(
                    f"Transcript {t_id} exons {i - 1} and {i} overlap: "
                    f"({exons[i - 1].start},{exons[i - 1].end}) vs "
                    f"({exons[i].start},{exons[i].end})"
                )

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

        transcripts.sort(key=lambda t: (t.ref, t.start, t.end, t.strand))

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
                    gtf_obj = GTFRecord(
                        seqname=t.ref,
                        source="rigel_sim",
                        feature="exon",
                        start=exon.start,  # GTFRecord.__str__ adds +1
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
            f"({sum(len(t.exons) for t in transcripts)} exons) → {gtf_path}"
        )
        return gtf_path
