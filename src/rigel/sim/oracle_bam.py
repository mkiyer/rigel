"""
rigel.sim.oracle_bam — Oracle BAM fragment simulator.

Generates a name-sorted BAM file with perfectly aligned paired-end
reads, bypassing the FASTQ → alignment → BAM pipeline.  This isolates
estimation errors from alignment errors in benchmarking.

Coordinate projection
---------------------
The underlying ``ReadSimulator`` generates fragments in three spaces:

1. **Mature RNA** — positions on the *spliced* transcript (exon-only
   concatenation).  These must be projected back through the exon
   coordinate map to produce split CIGAR alignments with ``N``
   (reference-skip) operations for introns.

2. **Nascent RNA** — positions on the *pre-mRNA* (full genomic span,
   including introns).  Since the pre-mRNA is a contiguous genomic
   region, these produce simple ``M`` (match) alignments.

3. **gDNA** — positions already in genomic coordinates.  Simple ``M``
   CIGAR, with strand chosen randomly.

BAM conventions
~~~~~~~~~~~~~~~
- Name-sorted (``SO:queryname``) — the format ``run_pipeline`` expects.
- Proper paired-end flags: ``0x1``, ``0x2``, ``0x40``/``0x80``,
  ``0x10``/``0x20`` (reverse strand).
- CIGAR: ``M`` for contiguous alignment, ``N`` for intron skips.
- ``XS`` tag set to the donor–acceptor strand for spliced reads
  (minimap2/STAR convention) so rigel's splice-strand logic works.
- ``NH`` tag = 1 (unique mapping, since we know the true origin).

Read-name format
~~~~~~~~~~~~~~~~
Identical to ``reads.py``: ``{t_id}:{frag_start}-{frag_end}:{strand}:{idx}``.
This preserves compatibility with ``parse_truth_from_fastq``.
"""

import logging
from pathlib import Path

import numpy as np
import pysam

from ..transcript import Transcript
from ..types import Strand
from .bam import (
    BASE_R1_FLAG,
    BASE_R2_FLAG,
    FLAG_MATE_REVERSE,
    FLAG_REVERSE,
    blocks_to_cigar,
    make_aligned_segment,
    premrna_to_genomic_interval,
    take_from_left,
    take_from_right,
    transcript_to_genomic_blocks,
)
from .capture import CaptureConfig
from .genome import MutableGenome
from .reads import GDNAConfig, ReadSimulator, ReadSimConfig

logger = logging.getLogger(__name__)

__all__ = ["OracleBamSimulator"]

# ---------------------------------------------------------------------------
# OracleBamSimulator
# ---------------------------------------------------------------------------


class OracleBamSimulator:
    """Generate a name-sorted BAM with perfect alignments.

    Reuses the fragment-selection logic of :class:`ReadSimulator`
    (pool splitting, abundance weighting, fragment-length sampling)
    but writes aligned BAM records instead of FASTQ, with CIGAR
    strings derived from the known transcript structure.

    Parameters
    ----------
    genome : MutableGenome
        The reference genome (or region) from which transcripts
        are annotated.
    transcripts : list[Transcript]
        Transcript annotations with exon coordinates and abundance.
        Exon coordinates are relative to the genome/region.
    config : ReadSimConfig or None
        Simulation parameters (fragment size, read length, etc.).
    gdna_config : GDNAConfig or None
        gDNA contamination parameters. ``None`` disables gDNA.
    ref_name : str
        Reference sequence name for the BAM header (e.g. ``"chr1"``
        or a region label like ``"HBB_locus"``).
    """

    def __init__(
        self,
        genome: MutableGenome,
        transcripts: list[Transcript],
        config: ReadSimConfig | None = None,
        gdna_config: GDNAConfig | None = None,
        capture_config: CaptureConfig | None = None,
        ref_name: str = "ref",
    ):
        self.genome = genome
        self.transcripts = transcripts
        self.config = config or ReadSimConfig()
        self.gdna_config = gdna_config
        self.ref_name = ref_name

        # Internal ReadSimulator for fragment selection logic
        self._sim = ReadSimulator(
            genome,
            transcripts,
            config=self.config,
            gdna_config=gdna_config,
            capture_config=capture_config,
        )

        # Pre-compute transcript lengths and pre-mRNA lengths
        self._t_lengths = self._sim._t_lengths
        self._premrna_lengths = self._sim._premrna_lengths

    def write_bam(
        self,
        output_path: Path,
        n_fragments: int,
        *,
        name_sorted: bool = True,
        pool_split: tuple[int, int, int] | None = None,
    ) -> Path:
        """Simulate fragments and write a BAM file.

        Parameters
        ----------
        output_path : Path
            Path for the output BAM file.
        n_fragments : int
            Total number of fragments to simulate.
        name_sorted : bool
            If True (default), produce a name-sorted BAM suitable for
            ``run_pipeline``.  If False, produce a coordinate-sorted BAM.
        pool_split : tuple[int, int, int] or None
            Explicit ``(n_mrna, n_nrna, n_gdna)`` split.  When provided,
            the abundance-based pool split is bypassed and *n_fragments*
            is ignored (the sum of the split is used instead).

        Returns
        -------
        Path
            The path to the written BAM file.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # BAM header
        genome_len = len(self.genome)
        header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "queryname" if name_sorted else "coordinate"},
            "SQ": [{"SN": self.ref_name, "LN": genome_len}],
            "PG": [{"ID": "rigel_oracle_bam", "PN": "rigel_oracle_bam",
                     "VN": "1.0", "CL": "simulated"}],
        })

        ref_id = 0  # single reference
        records: list[pysam.AlignedSegment] = []

        if pool_split is not None:
            n_mrna, n_nrna, n_gdna = pool_split
        else:
            n_mrna, n_nrna, n_gdna = self._sim._compute_pool_split(n_fragments)
        rng = self._sim._rng
        ntranscripts = len(self.transcripts)
        cfg = self.config
        read_len = cfg.read_length

        # ----- Mature RNA fragments -----
        if n_mrna > 0:
            records.extend(
                self._generate_mrna_records(
                    header, ref_id, n_mrna, rng, ntranscripts, read_len,
                )
            )

        # ----- Nascent RNA fragments -----
        if n_nrna > 0:
            records.extend(
                self._generate_nrna_records(
                    header, ref_id, n_nrna, rng, ntranscripts, read_len,
                )
            )

        # ----- gDNA fragments -----
        if n_gdna > 0:
            records.extend(
                self._generate_gdna_records(
                    header, ref_id, n_gdna, rng, read_len,
                )
            )

        # Sort records
        if name_sorted:
            records.sort(key=lambda r: r.query_name)
        else:
            records.sort(key=lambda r: (r.reference_start, r.query_name))

        # Write BAM
        with pysam.AlignmentFile(str(output_path), "wb", header=header) as outf:
            for rec in records:
                outf.write(rec)

        # Index if coordinate-sorted
        if not name_sorted:
            pysam.index(str(output_path))

        n_pairs = len(records) // 2
        logger.info(
            "Wrote %d read pairs to %s (mRNA=%d, nRNA=%d, gDNA=%d)",
            n_pairs, output_path, n_mrna, n_nrna, n_gdna,
        )
        return output_path

    # -- mRNA fragment generation -------------------------------------------

    def _generate_mrna_records(
        self,
        header: pysam.AlignmentHeader,
        ref_id: int,
        n_mrna: int,
        rng: np.random.Generator,
        ntranscripts: int,
        read_len: int,
    ) -> list[pysam.AlignedSegment]:
        """Generate BAM records for mature RNA fragments."""
        records: list[pysam.AlignedSegment] = []
        sim = self._sim
        cfg = self.config

        frag_lengths = sim._sample_frag_lengths(n_mrna)
        unique_lengths, length_counts = np.unique(frag_lengths, return_counts=True)

        for frag_len, frag_count in zip(unique_lengths, length_counts):
            probs = sim._compute_probs(int(frag_len))
            t_indices = rng.choice(ntranscripts, size=int(frag_count), p=probs)
            unique_t, t_counts = np.unique(t_indices, return_counts=True)

            for t_idx, t_count in zip(unique_t, t_counts):
                t = self.transcripts[int(t_idx)]
                t_len = self._t_lengths[int(t_idx)]
                rl = min(read_len, int(frag_len))
                eff_len = t_len - int(frag_len) + 1
                if eff_len <= 0:
                    continue

                frag_starts = sim.capture.sample_starts(
                    "mrna", int(t_idx), t_len, int(frag_len), int(t_count), rng,
                )

                # Strand-flip mask
                ss = cfg.strand_specificity
                if ss < 1.0:
                    flip_mask = rng.random(int(t_count)) >= ss
                else:
                    flip_mask = None

                strand_char = "r" if t.strand == Strand.NEG else "f"

                for i in range(int(t_count)):
                    frag_start = int(frag_starts[i])
                    frag_end = frag_start + int(frag_len)

                    # Project to genomic coordinates
                    blocks = transcript_to_genomic_blocks(frag_start, frag_end, t)
                    if not blocks:
                        continue

                    # Build read sequences from the genome
                    # The fragment in genome orientation (ascending coords)
                    frag_seq_parts = [
                        self.genome[bs:be] for bs, be in blocks
                    ]
                    frag_genomic_seq = "".join(frag_seq_parts)

                    # Determine whether this fragment is flipped (antisense)
                    flipped = flip_mask is not None and flip_mask[i]

                    rname = f"{t.t_id}:{frag_start}-{frag_end}:{strand_char}:{i}"

                    recs = self._make_paired_records(
                        header=header,
                        ref_id=ref_id,
                        rname=rname,
                        blocks=blocks,
                        frag_genomic_seq=frag_genomic_seq,
                        read_len=rl,
                        transcript=t,
                        flipped=flipped,
                        is_spliced=len(blocks) > 1,
                    )
                    records.extend(recs)

        return records

    # -- nRNA fragment generation -------------------------------------------

    def _generate_nrna_records(
        self,
        header: pysam.AlignmentHeader,
        ref_id: int,
        n_nrna: int,
        rng: np.random.Generator,
        ntranscripts: int,
        read_len: int,
    ) -> list[pysam.AlignedSegment]:
        """Generate BAM records for nascent RNA fragments."""
        records: list[pysam.AlignedSegment] = []
        sim = self._sim
        cfg = self.config

        frag_lengths = sim._sample_frag_lengths(n_nrna)
        unique_lengths, length_counts = np.unique(frag_lengths, return_counts=True)

        for frag_len, frag_count in zip(unique_lengths, length_counts):
            probs = sim._compute_nrna_probs(int(frag_len))
            t_indices = rng.choice(ntranscripts, size=int(frag_count), p=probs)
            unique_t, t_counts = np.unique(t_indices, return_counts=True)

            for t_idx, t_count in zip(unique_t, t_counts):
                t = self.transcripts[int(t_idx)]
                premrna_len = self._premrna_lengths[int(t_idx)]
                rl = min(read_len, int(frag_len))
                eff_len = premrna_len - int(frag_len) + 1
                if eff_len <= 0:
                    continue

                frag_starts = sim.capture.sample_starts(
                    "nrna", int(t_idx), premrna_len, int(frag_len), int(t_count), rng,
                )

                ss = cfg.strand_specificity
                if ss < 1.0:
                    flip_mask = rng.random(int(t_count)) >= ss
                else:
                    flip_mask = None

                strand_char = "r" if t.strand == Strand.NEG else "f"

                for i in range(int(t_count)):
                    frag_start = int(frag_starts[i])
                    frag_end = frag_start + int(frag_len)

                    gstart, gend = premrna_to_genomic_interval(
                        frag_start, frag_end, t,
                    )
                    blocks = [(gstart, gend)]
                    frag_genomic_seq = self.genome[gstart:gend]

                    flipped = flip_mask is not None and flip_mask[i]
                    rname = f"nrna_{t.t_id}:{frag_start}-{frag_end}:{strand_char}:{i}"

                    recs = self._make_paired_records(
                        header=header,
                        ref_id=ref_id,
                        rname=rname,
                        blocks=blocks,
                        frag_genomic_seq=frag_genomic_seq,
                        read_len=rl,
                        transcript=t,
                        flipped=flipped,
                        is_spliced=False,
                    )
                    records.extend(recs)

        return records

    # -- gDNA fragment generation -------------------------------------------

    def _generate_gdna_records(
        self,
        header: pysam.AlignmentHeader,
        ref_id: int,
        n_gdna: int,
        rng: np.random.Generator,
        read_len: int,
    ) -> list[pysam.AlignedSegment]:
        """Generate BAM records for gDNA fragments."""
        records: list[pysam.AlignedSegment] = []
        sim = self._sim
        genome_len = len(self.genome)

        gdna_frag_lengths = sim._sample_gdna_frag_lengths(n_gdna)
        unique_lengths, length_counts = np.unique(
            gdna_frag_lengths, return_counts=True,
        )

        global_idx = 0
        for frag_len, frag_count in zip(unique_lengths, length_counts):
            eff_len = genome_len - int(frag_len) + 1
            if eff_len <= 0:
                continue

            rl = min(read_len, int(frag_len))
            frag_starts = sim.capture.sample_starts(
                "gdna", self.genome.name, genome_len, int(frag_len), int(frag_count), rng,
            )
            strands = rng.integers(0, 2, size=int(frag_count))  # 0=+, 1=-

            for i in range(int(frag_count)):
                gstart = int(frag_starts[i])
                gend = gstart + int(frag_len)
                is_reverse_strand = bool(strands[i])
                strand_char = "r" if is_reverse_strand else "f"

                frag_genomic_seq = self.genome[gstart:gend]

                rname = f"gdna:{gstart}-{gend}:{strand_char}:{global_idx}"

                # For gDNA, create a synthetic "transcript" context for
                # the paired-record builder.  gDNA strand info is in
                # strand_char — we use it to set read orientations.
                recs = self._make_gdna_paired_records(
                    header=header,
                    ref_id=ref_id,
                    rname=rname,
                    gstart=gstart,
                    gend=gend,
                    frag_genomic_seq=frag_genomic_seq,
                    read_len=rl,
                    is_reverse_strand=is_reverse_strand,
                )
                records.extend(recs)
                global_idx += 1

        return records

    # -- Paired record builders -----------------------------------------------

    def _make_paired_records(
        self,
        *,
        header: pysam.AlignmentHeader,
        ref_id: int,
        rname: str,
        blocks: list[tuple[int, int]],
        frag_genomic_seq: str,
        read_len: int,
        transcript: Transcript,
        flipped: bool,
        is_spliced: bool,
    ) -> list[pysam.AlignedSegment]:
        """Build R1 + R2 BAM records for an RNA fragment.

        R1-antisense library convention:
        - POS-strand transcript: R1 maps to − (reverse), R2 maps to + (forward)
        - NEG-strand transcript: R1 maps to + (forward), R2 maps to − (reverse)

        When ``flipped=True`` (imperfect strandedness), orientations swap.
        """
        neg_tx = transcript.strand == Strand.NEG

        # Determine read orientations:
        # Normal FR: R1 is from the 3' end (mapped reverse), R2 from 5' (forward)
        # For POS-strand tx: R1 reverse, R2 forward (in genomic coords)
        # For NEG-strand tx: R1 forward, R2 reverse (in genomic coords)
        if not flipped:
            if not neg_tx:
                r1_is_reverse = True   # POS tx: R1 maps to −
                r2_is_reverse = False  # POS tx: R2 maps to +
            else:
                r1_is_reverse = False  # NEG tx: R1 maps to +
                r2_is_reverse = True   # NEG tx: R2 maps to −
        else:
            # Flipped = swap R1↔R2 orientation
            if not neg_tx:
                r1_is_reverse = False
                r2_is_reverse = True
            else:
                r1_is_reverse = True
                r2_is_reverse = False

        # R2 aligns to the 5' genomic blocks, R1 to the 3' genomic blocks
        # In genomic coordinates (ascending), for POS-strand tx:
        #   R2 covers the leftmost bases (5' end), R1 covers rightmost (3' end)
        # For NEG-strand tx:
        #   R2 covers rightmost bases (5' end on −strand), R1 covers leftmost (3' end on −strand)

        # Compute R1 and R2 genomic blocks
        r1_blocks, r2_blocks = self._split_fragment_blocks(
            blocks, read_len, neg_tx,
        )
        if flipped:
            r1_blocks, r2_blocks = r2_blocks, r1_blocks

        # Build sequences
        # BAM SEQ is stored in forward-reference orientation: it must match
        # the reference at the aligned position (left-to-right along the
        # reference).  The 0x10 flag tells downstream tools that the
        # original sequencer read was reverse-complemented to produce this
        # stored SEQ.  Since oracle reads have zero errors, SEQ == reference.
        r1_seq = "".join(self.genome[bs:be] for bs, be in r1_blocks)
        r2_seq = "".join(self.genome[bs:be] for bs, be in r2_blocks)

        # CIGAR
        r1_cigar = blocks_to_cigar(r1_blocks)
        r2_cigar = blocks_to_cigar(r2_blocks)

        # Positions
        r1_pos = r1_blocks[0][0]
        r2_pos = r2_blocks[0][0]

        # Template length: signed distance from leftmost to rightmost
        frag_left = blocks[0][0]
        frag_right = blocks[-1][1]
        tlen = frag_right - frag_left

        # XS tag for spliced reads (splice-junction strand)
        tags: list[tuple] = [("NH", 1)]
        if is_spliced:
            xs_strand = "+" if transcript.strand == Strand.POS else "-"
            tags.append(("XS", xs_strand, "A"))

        # Build R1 flags
        r1_flag = BASE_R1_FLAG
        if r1_is_reverse:
            r1_flag |= FLAG_REVERSE
        if r2_is_reverse:
            r1_flag |= FLAG_MATE_REVERSE

        # Build R2 flags
        r2_flag = BASE_R2_FLAG
        if r2_is_reverse:
            r2_flag |= FLAG_REVERSE
        if r1_is_reverse:
            r2_flag |= FLAG_MATE_REVERSE

        r1_rec = make_aligned_segment(
            header=header,
            query_name=rname,
            query_sequence=r1_seq,
            flag=r1_flag,
            reference_id=ref_id,
            reference_start=r1_pos,
            cigar=r1_cigar,
            mate_reference_id=ref_id,
            mate_reference_start=r2_pos,
            template_length=tlen if r1_pos <= r2_pos else -tlen,
            tags=tags,
        )

        # R2 gets XS tag too if spliced
        r2_tags: list[tuple] = [("NH", 1)]
        if is_spliced and len(r2_blocks) > 1:
            xs_strand = "+" if transcript.strand == Strand.POS else "-"
            r2_tags.append(("XS", xs_strand, "A"))

        r2_rec = make_aligned_segment(
            header=header,
            query_name=rname,
            query_sequence=r2_seq,
            flag=r2_flag,
            reference_id=ref_id,
            reference_start=r2_pos,
            cigar=r2_cigar,
            mate_reference_id=ref_id,
            mate_reference_start=r1_pos,
            template_length=tlen if r2_pos <= r1_pos else -tlen,
            tags=r2_tags,
        )

        return [r1_rec, r2_rec]

    def _make_gdna_paired_records(
        self,
        *,
        header: pysam.AlignmentHeader,
        ref_id: int,
        rname: str,
        gstart: int,
        gend: int,
        frag_genomic_seq: str,
        read_len: int,
        is_reverse_strand: bool,
    ) -> list[pysam.AlignedSegment]:
        """Build R1 + R2 BAM records for a gDNA fragment.

        gDNA fragments are unspliced. The same read orientation convention applies:
        if the fragment is on the + strand, R1 is reverse-complement of the
        3' end, R2 is sense from the 5' end (and vice versa for − strand).
        """
        gend - gstart

        if not is_reverse_strand:
            # Fragment on + strand:
            # R2 = 5' end (left, forward), R1 = 3' end (right, reverse)
            r2_start = gstart
            r2_end = min(gstart + read_len, gend)
            r1_start = max(gend - read_len, gstart)
            r1_end = gend
            r1_is_reverse = True
            r2_is_reverse = False
        else:
            # Fragment on − strand:
            # R2 = 5' end (right, reverse), R1 = 3' end (left, forward)
            r1_start = gstart
            r1_end = min(gstart + read_len, gend)
            r2_start = max(gend - read_len, gstart)
            r2_end = gend
            r1_is_reverse = False
            r2_is_reverse = True

        # BAM SEQ = reference sequence for perfect oracle alignments.
        r1_seq = self.genome[r1_start:r1_end]
        r2_seq = self.genome[r2_start:r2_end]

        r1_cigar = [(pysam.CMATCH, r1_end - r1_start)]
        r2_cigar = [(pysam.CMATCH, r2_end - r2_start)]

        tlen = gend - gstart

        # R1 flags
        r1_flag = BASE_R1_FLAG
        if r1_is_reverse:
            r1_flag |= FLAG_REVERSE
        if r2_is_reverse:
            r1_flag |= FLAG_MATE_REVERSE

        # R2 flags
        r2_flag = BASE_R2_FLAG
        if r2_is_reverse:
            r2_flag |= FLAG_REVERSE
        if r1_is_reverse:
            r2_flag |= FLAG_MATE_REVERSE

        tags: list[tuple] = [("NH", 1)]

        r1_rec = make_aligned_segment(
            header=header,
            query_name=rname,
            query_sequence=r1_seq,
            flag=r1_flag,
            reference_id=ref_id,
            reference_start=r1_start,
            cigar=r1_cigar,
            mate_reference_id=ref_id,
            mate_reference_start=r2_start,
            template_length=tlen if r1_start <= r2_start else -tlen,
            tags=tags,
        )

        r2_rec = make_aligned_segment(
            header=header,
            query_name=rname,
            query_sequence=r2_seq,
            flag=r2_flag,
            reference_id=ref_id,
            reference_start=r2_start,
            cigar=r2_cigar,
            mate_reference_id=ref_id,
            mate_reference_start=r1_start,
            template_length=tlen if r2_start <= r1_start else -tlen,
            tags=tags,
        )

        return [r1_rec, r2_rec]

    @staticmethod
    def _split_fragment_blocks(
        blocks: list[tuple[int, int]],
        read_len: int,
        neg_strand_tx: bool,
    ) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
        """Split fragment genomic blocks into R1 and R2 sub-blocks.

        For a POS-strand transcript (R1-antisense convention):
        - R2 covers the leftmost `read_len` bases (5' end)
        - R1 covers the rightmost `read_len` bases (3' end)

        For a NEG-strand transcript:
        - R2 covers the rightmost `read_len` bases (5' end on − strand)
        - R1 covers the leftmost `read_len` bases (3' end on − strand)

        Returns (r1_blocks, r2_blocks) — both in ascending genomic order.
        """
        if neg_strand_tx:
            # NEG strand: 5' end is rightmost, 3' end is leftmost
            r1_blocks = take_from_left(blocks, read_len)
            r2_blocks = take_from_right(blocks, read_len)
        else:
            # POS strand: 5' end is leftmost, 3' end is rightmost
            r2_blocks = take_from_left(blocks, read_len)
            r1_blocks = take_from_right(blocks, read_len)

        return r1_blocks, r2_blocks

