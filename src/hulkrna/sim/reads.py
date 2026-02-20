"""
hulkrna.sim.reads — Paired-end RNA-seq read simulator with gDNA support.

Generates paired-end FASTQ reads from transcript sequences with
configurable fragment-size distribution, read length, and error rate.
Optionally simulates genomic DNA (gDNA) contamination from the full
genome with an independent fragment-size profile.

Read names encode ground-truth origin for post-hoc validation.

Read-name format (RNA)::

    {t_id}:{frag_start}-{frag_end}:{strand_char}:{index}/1

Read-name format (gDNA)::

    gdna:{genomic_start}-{genomic_end}:{strand_char}:{index}/1

Library-prep convention (FR orientation)
----------------------------------------
Reads simulate a standard Illumina stranded (dUTP) library:

- **Positive-strand transcript**: R1 is reverse-complement of 3′ end
  of the fragment, R2 is the sense sequence from the 5′ end.
  When aligned, R1 maps to − strand and R2 maps to + strand.
  ``Fragment.from_reads`` flips R2's strand → both report + → correct.

- **Negative-strand transcript**: R1 maps to + strand, R2 maps to − strand.
  After R2 flip → both report − → correct.

- **gDNA**: sampled from both strands equally (unstranded).
  FR library convention still applies for the orientation of reads.

Strand specificity
------------------
The ``strand_specificity`` parameter in ``SimConfig`` controls the
fraction of RNA fragments that preserve correct FR orientation.
At 1.0, all reads are perfectly stranded.  At 0.5, reads are
effectively unstranded (50/50 chance of correct orientation).
Imperfect strandedness is implemented by swapping R1↔R2 with
probability ``1 − strand_specificity`` per fragment.
"""

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..types import Strand
from ..transcript import Transcript
from .genome import MutableGenome, reverse_complement

logger = logging.getLogger(__name__)

__all__ = ["SimConfig", "GDNAConfig", "ReadSimulator"]

# DNA bases for error injection
_DNA_BASES = np.array(list("ACGT"), dtype="<U1")


@dataclass
class SimConfig:
    """Configuration for the read simulator.

    Attributes
    ----------
    frag_mean : float
        Mean fragment length (normal distribution).
    frag_std : float
        Standard deviation of fragment length.
    frag_min : int
        Minimum allowed fragment length.
    frag_max : int
        Maximum allowed fragment length.
    read_length : int
        Read length (bases). Both R1 and R2 have the same length.
    error_rate : float
        Per-base substitution error rate (0.0 = no errors).
    strand_specificity : float
        Probability that an RNA fragment preserves its correct
        FR orientation.  1.0 = perfectly stranded (dUTP library),
        0.5 = completely unstranded (no strand information).
        Intermediate values simulate imperfect dUTP incorporation.
        Implemented by swapping R1↔R2 with probability
        ``1 − strand_specificity`` per fragment.
    seed : int
        Random seed for reproducibility.
    """

    frag_mean: float = 250.0
    frag_std: float = 50.0
    frag_min: int = 50
    frag_max: int = 1000
    read_length: int = 150
    error_rate: float = 0.0
    strand_specificity: float = 1.0
    seed: int = 42


@dataclass
class GDNAConfig:
    """Configuration for genomic DNA contamination simulation.

    gDNA fragments are sampled uniformly from the full genome
    (both strands equally) with an independent fragment-size
    distribution.  The ``abundance`` parameter uses the same
    relative scale as transcript abundances — the fraction of
    total fragments that are gDNA is determined by the ratio
    of gDNA abundance × effective_length to total weight.

    Attributes
    ----------
    abundance : float
        Relative abundance of gDNA (same scale as transcript
        abundances).  For example, if a transcript has
        abundance=100 and gDNA has abundance=10, roughly
        ``10 × genome_eff_len / (100 × t_eff_len + 10 × genome_eff_len)``
        of fragments will be gDNA.
    frag_mean : float
        Mean gDNA fragment length (typically larger than RNA).
    frag_std : float
        Standard deviation of gDNA fragment length.
    frag_min : int
        Minimum gDNA fragment length.
    frag_max : int
        Maximum gDNA fragment length.
    """

    abundance: float = 10.0
    frag_mean: float = 350.0
    frag_std: float = 100.0
    frag_min: int = 100
    frag_max: int = 1000


class ReadSimulator:
    """Paired-end RNA-seq read simulator with optional gDNA contamination.

    Generates reads from transcript sequences extracted from the genome
    by concatenating exon coordinates. Fragment lengths are drawn from
    a truncated normal distribution. Transcript selection is weighted
    by ``abundance × effective_length``.

    When ``gdna_config`` is provided, gDNA fragments are also simulated.
    gDNA is sampled uniformly from the full genome (both strands equally)
    with its own fragment-size distribution. The total ``n_fragments``
    is split between RNA and gDNA proportionally based on their
    abundance-weighted effective lengths.

    Parameters
    ----------
    genome : MutableGenome
        The genome from which transcript sequences are extracted.
    transcripts : list[Transcript]
        Transcript annotations with exon coordinates and abundance.
    config : SimConfig or None
        Simulation parameters. ``None`` uses defaults.
    gdna_config : GDNAConfig or None
        gDNA contamination config. ``None`` disables gDNA simulation.
    """

    # Oversampling ratio to reduce rejection-sampling iterations
    _OVERSAMPLE_RATIO = 1.5
    _OVERSAMPLE_EXTRA = 10

    def __init__(
        self,
        genome: MutableGenome,
        transcripts: list[Transcript],
        config: SimConfig | None = None,
        gdna_config: GDNAConfig | None = None,
    ):
        self.genome = genome
        self.config = config or SimConfig()
        self.gdna_config = gdna_config
        self._rng = np.random.default_rng(self.config.seed)

        # Extract transcript sequences and compute lengths
        self.transcripts = transcripts
        self._t_seqs: list[str] = []
        self._t_lengths: list[int] = []
        max_len = 0

        for t in transcripts:
            seq = self._extract_transcript_seq(t)
            self._t_seqs.append(seq)
            self._t_lengths.append(len(seq))
            max_len = max(max_len, len(seq))

        # Extract pre-mRNA sequences (unspliced, full genomic span)
        self._premrna_seqs: list[str] = []
        self._premrna_lengths: list[int] = []
        for t in transcripts:
            premrna_seq = self._extract_premrna_seq(t)
            self._premrna_seqs.append(premrna_seq)
            self._premrna_lengths.append(len(premrna_seq))

        # Clamp frag_max to longest transcript
        self.config.frag_max = min(self.config.frag_max, max_len)

    def _extract_transcript_seq(self, t: Transcript) -> str:
        """Concatenate exonic sequences from the genome.

        For negative-strand transcripts, the concatenated exon sequence
        is reverse-complemented so that ``t_seq[0]`` is the 5′ end of
        the mRNA.
        """
        exon_seqs = [self.genome[e.start:e.end] for e in t.exons]
        seq = "".join(exon_seqs)
        if t.strand == Strand.NEG:
            seq = reverse_complement(seq)
        return seq

    def _extract_premrna_seq(self, t: Transcript) -> str:
        """Extract the unspliced pre-mRNA sequence (full genomic span).

        For multi-exon transcripts, this includes introns.  For
        single-exon transcripts, the pre-mRNA is identical to the
        spliced mRNA.  Negative-strand sequences are reverse-
        complemented.
        """
        seq = self.genome[t.start:t.end]
        if t.strand == Strand.NEG:
            seq = reverse_complement(seq)
        return seq

    # -- Fragment-length sampling ---------------------------------------------

    def _sample_frag_lengths(self, n: int) -> np.ndarray:
        """Sample *n* fragment lengths from a truncated normal."""
        cfg = self.config
        return self._sample_frag_lengths_trunc_normal(
            n, cfg.frag_mean, cfg.frag_std, cfg.frag_min, cfg.frag_max,
        )

    def _sample_gdna_frag_lengths(self, n: int) -> np.ndarray:
        """Sample *n* gDNA fragment lengths from a truncated normal."""
        gc = self.gdna_config
        assert gc is not None
        frag_max = min(gc.frag_max, len(self.genome))
        return self._sample_frag_lengths_trunc_normal(
            n, gc.frag_mean, gc.frag_std, gc.frag_min, frag_max,
        )

    def _sample_frag_lengths_trunc_normal(
        self,
        n: int,
        mean: float,
        std: float,
        frag_min: int,
        frag_max: int,
    ) -> np.ndarray:
        """Sample *n* fragment lengths from a truncated normal distribution."""
        rng = self._rng
        result = np.empty(n, dtype=int)
        filled = 0

        while filled < n:
            needed = n - filled
            sample_size = int(needed * self._OVERSAMPLE_RATIO) + self._OVERSAMPLE_EXTRA
            raw = rng.normal(mean, std, sample_size).astype(int)
            valid = raw[(raw >= frag_min) & (raw <= frag_max)]
            nkeep = min(len(valid), needed)
            result[filled : filled + nkeep] = valid[:nkeep]
            filled += nkeep

        return result

    # -- Transcript probability -----------------------------------------------

    def _compute_probs(self, frag_length: int) -> np.ndarray:
        """Compute mature RNA transcript selection probabilities.

        Weight = max(0, spliced_length - frag_length + 1) × abundance

        ``abundance`` is the mature mRNA molecular abundance and is
        independent of ``nrna_abundance``.  These are separate pools:
        the probability of drawing a mature RNA fragment is proportional
        to ``abundance × spliced_effective_length``, regardless of how
        much nascent RNA is present.
        """
        weights = np.array(
            [
                max(0, tlen - frag_length + 1) * (t.abundance or 0)
                for tlen, t in zip(self._t_lengths, self.transcripts)
            ]
        )
        total = weights.sum()
        if total == 0:
            # All transcripts shorter than frag_length — uniform fallback
            return np.ones(len(self.transcripts)) / len(self.transcripts)
        return weights / total

    def _compute_nrna_probs(self, frag_length: int) -> np.ndarray:
        """Compute nascent RNA transcript selection probabilities.

        Weight = max(0, premrna_length - frag_length + 1) × nrna_abundance

        ``nrna_abundance`` is the nascent RNA molecular abundance and is
        independent of ``abundance`` (the mature mRNA molecular abundance).
        The pre-mRNA effective length (full genomic span including introns)
        is typically much larger than the spliced effective length, so a
        small ``nrna_abundance`` can still produce a substantial fraction
        of reads.
        """
        weights = np.array(
            [
                max(0, plen - frag_length + 1) * t.nrna_abundance
                for plen, t in zip(self._premrna_lengths, self.transcripts)
            ]
        )
        total = weights.sum()
        if total == 0:
            return np.ones(len(self.transcripts)) / len(self.transcripts)
        return weights / total

    # -- Read generation ------------------------------------------------------

    def _introduce_errors(self, seq: str) -> str:
        """Introduce random substitution errors."""
        if self.config.error_rate <= 0:
            return seq
        arr = np.array(list(seq), dtype="<U1")
        mask = self._rng.random(len(seq)) < self.config.error_rate
        if not np.any(mask):
            return seq
        arr[mask] = _DNA_BASES[self._rng.integers(4, size=int(mask.sum()))]
        return "".join(arr)

    def _gen_reads_from_transcript(
        self, t_idx: int, frag_len: int, count: int
    ):
        """Generate read pairs from a transcript.

        Yields (r1_name, r1_seq, r1_qual, r2_name, r2_seq, r2_qual).
        """
        t = self.transcripts[t_idx]
        t_seq = self._t_seqs[t_idx]
        t_len = self._t_lengths[t_idx]
        read_len = min(self.config.read_length, frag_len)
        rng = self._rng

        eff_len = t_len - frag_len + 1
        if eff_len <= 0:
            return

        # Sample fragment start positions (on the mRNA/transcript)
        frag_starts = rng.integers(0, eff_len, size=count)

        # Pre-compute strand-flip mask for imperfect strandedness.
        # When strand_specificity < 1.0, some fragments get R1↔R2
        # swapped, making them appear antisense to the pipeline.
        ss = self.config.strand_specificity
        if ss < 1.0:
            flip_mask = rng.random(count) >= ss
        else:
            flip_mask = None

        strand_char = "r" if t.strand == Strand.NEG else "f"

        for i in range(count):
            frag_start = frag_starts[i]
            frag_end = frag_start + frag_len
            frag_seq = t_seq[frag_start:frag_end]

            # FR library: R1 is reverse-complement of 3′ end,
            # R2 is sense sequence from 5′ end.
            r1_seq = reverse_complement(frag_seq[-read_len:])
            r2_seq = frag_seq[:read_len]

            # Simulate imperfect strandedness: swap R1↔R2 so the
            # pipeline sees the fragment as coming from the opposite
            # strand.  This models incomplete dUTP incorporation.
            if flip_mask is not None and flip_mask[i]:
                r1_seq, r2_seq = r2_seq, r1_seq

            # Apply errors
            r1_seq = self._introduce_errors(r1_seq)
            r2_seq = self._introduce_errors(r2_seq)

            # Quality scores (constant high quality for simulation)
            quals = "I" * read_len

            # Read name encodes ground truth
            rname = f"{t.t_id}:{frag_start}-{frag_end}:{strand_char}:{i}"

            yield (
                f"{rname}/1", r1_seq, quals,
                f"{rname}/2", r2_seq, quals,
            )

    def _gen_reads_from_premrna(
        self, t_idx: int, frag_len: int, count: int
    ):
        """Generate nascent RNA read pairs from unspliced pre-mRNA.

        Reads are sampled from the full genomic span of the transcript
        (including introns).  FR library convention and strandedness
        behaviour are identical to mature RNA reads.

        Yields (r1_name, r1_seq, r1_qual, r2_name, r2_seq, r2_qual).
        """
        t = self.transcripts[t_idx]
        premrna_seq = self._premrna_seqs[t_idx]
        premrna_len = self._premrna_lengths[t_idx]
        read_len = min(self.config.read_length, frag_len)
        rng = self._rng

        eff_len = premrna_len - frag_len + 1
        if eff_len <= 0:
            return

        frag_starts = rng.integers(0, eff_len, size=count)

        ss = self.config.strand_specificity
        if ss < 1.0:
            flip_mask = rng.random(count) >= ss
        else:
            flip_mask = None

        strand_char = "r" if t.strand == Strand.NEG else "f"

        for i in range(count):
            frag_start = frag_starts[i]
            frag_end = frag_start + frag_len
            frag_seq = premrna_seq[frag_start:frag_end]

            r1_seq = reverse_complement(frag_seq[-read_len:])
            r2_seq = frag_seq[:read_len]

            if flip_mask is not None and flip_mask[i]:
                r1_seq, r2_seq = r2_seq, r1_seq

            r1_seq = self._introduce_errors(r1_seq)
            r2_seq = self._introduce_errors(r2_seq)

            quals = "I" * read_len
            rname = f"nrna_{t.t_id}:{frag_start}-{frag_end}:{strand_char}:{i}"

            yield (
                f"{rname}/1", r1_seq, quals,
                f"{rname}/2", r2_seq, quals,
            )

    def _gen_reads_from_genome(self, frag_len: int, count: int):
        """Generate gDNA read pairs from the full genome.

        Fragments are sampled uniformly from the entire genome,
        with strand chosen uniformly (+/−).  FR library convention
        is applied as for RNA reads.

        Yields (r1_name, r1_seq, r1_qual, r2_name, r2_seq, r2_qual).
        """
        genome_len = len(self.genome)
        eff_len = genome_len - frag_len + 1
        if eff_len <= 0:
            return

        read_len = min(self.config.read_length, frag_len)
        rng = self._rng
        quals = "I" * read_len

        frag_starts = rng.integers(0, eff_len, size=count)
        strands = rng.integers(0, 2, size=count)  # 0 = +, 1 = −

        for i in range(count):
            start = int(frag_starts[i])
            end = start + frag_len
            frag_seq = self.genome[start:end]

            if strands[i] == 1:
                frag_seq = reverse_complement(frag_seq)
                strand_char = "r"
            else:
                strand_char = "f"

            # FR library: R1 = revcomp of 3' end, R2 = sense from 5' end
            r1_seq = reverse_complement(frag_seq[-read_len:])
            r2_seq = frag_seq[:read_len]

            r1_seq = self._introduce_errors(r1_seq)
            r2_seq = self._introduce_errors(r2_seq)

            rname = f"gdna:{start}-{end}:{strand_char}:{i}"

            yield (
                f"{rname}/1", r1_seq, quals,
                f"{rname}/2", r2_seq, quals,
            )

    # -- Public interface -----------------------------------------------------

    def _compute_pool_split(
        self, n_fragments: int,
    ) -> tuple[int, int, int]:
        """Compute the 3-way mature RNA / nascent RNA / gDNA split.

        Each pool's weight is its molecular abundance multiplied by its
        effective length — the number of fragment-start positions that
        yield a fragment of the mean length:

        - **Mature RNA**: ``abundance × max(spliced_len − mean_frag + 1, 0)``
        - **Nascent RNA**: ``nrna_abundance × max(premrna_span − mean_frag + 1, 0)``
        - **gDNA**: ``gdna_abundance × max(genome_len − gdna_mean_frag + 1, 0)``

        ``abundance`` and ``nrna_abundance`` are independent molecular
        abundances — they are **not** derived from one another by
        subtraction.  Setting ``nrna_abundance = k × abundance`` means
        there are *k* times as many nascent RNA molecules as mature mRNA
        molecules; the actual read fractions additionally depend on the
        ratio of effective lengths.

        Returns ``(n_mrna, n_nrna, n_gdna)``.
        """
        mean_frag = int(self.config.frag_mean)

        # Mature RNA weight: abundance × spliced effective length
        mrna_weight = sum(
            (t.abundance or 0) * max(0, tlen - mean_frag + 1)
            for t, tlen in zip(self.transcripts, self._t_lengths)
        )

        # Nascent RNA weight: nrna_abundance × pre-mRNA effective length
        nrna_weight = sum(
            t.nrna_abundance * max(0, plen - mean_frag + 1)
            for t, plen in zip(self.transcripts, self._premrna_lengths)
        )

        # gDNA weight: gdna_abundance × genome effective length
        if self.gdna_config is not None:
            gc = self.gdna_config
            gdna_mean_frag = int(gc.frag_mean)
            genome_eff_len = max(0, len(self.genome) - gdna_mean_frag + 1)
            gdna_weight = gc.abundance * genome_eff_len
        else:
            gdna_weight = 0.0

        total_weight = mrna_weight + nrna_weight + gdna_weight
        if total_weight <= 0:
            return n_fragments, 0, 0

        n_nrna = int(round(n_fragments * nrna_weight / total_weight))
        n_gdna = int(round(n_fragments * gdna_weight / total_weight))
        n_mrna = max(0, n_fragments - n_nrna - n_gdna)
        return n_mrna, n_nrna, n_gdna

    def simulate(self, n_fragments: int):
        """Generate *n_fragments* paired-end read pairs.

        Fragments are split three ways — mature RNA, nascent RNA,
        and gDNA — based on abundance-weighted effective lengths.
        Mature RNA reads are emitted first, then nascent RNA, then gDNA.

        Yields ``(r1_name, r1_seq, r1_qual, r2_name, r2_seq, r2_qual)``.
        """
        n_mrna, n_nrna, n_gdna = self._compute_pool_split(n_fragments)

        rng = self._rng
        ntranscripts = len(self.transcripts)

        # --- Mature RNA reads ---
        if n_mrna > 0:
            frag_lengths = self._sample_frag_lengths(n_mrna)
            unique_lengths, length_counts = np.unique(
                frag_lengths, return_counts=True,
            )

            for frag_len, frag_count in zip(unique_lengths, length_counts):
                probs = self._compute_probs(int(frag_len))
                t_indices = rng.choice(
                    ntranscripts, size=int(frag_count), p=probs,
                )
                unique_t, t_counts = np.unique(t_indices, return_counts=True)

                for t_idx, t_count in zip(unique_t, t_counts):
                    yield from self._gen_reads_from_transcript(
                        int(t_idx), int(frag_len), int(t_count),
                    )

        # --- Nascent RNA reads ---
        if n_nrna > 0:
            frag_lengths = self._sample_frag_lengths(n_nrna)
            unique_lengths, length_counts = np.unique(
                frag_lengths, return_counts=True,
            )

            for frag_len, frag_count in zip(unique_lengths, length_counts):
                probs = self._compute_nrna_probs(int(frag_len))
                t_indices = rng.choice(
                    ntranscripts, size=int(frag_count), p=probs,
                )
                unique_t, t_counts = np.unique(t_indices, return_counts=True)

                for t_idx, t_count in zip(unique_t, t_counts):
                    yield from self._gen_reads_from_premrna(
                        int(t_idx), int(frag_len), int(t_count),
                    )

        # --- gDNA reads ---
        if n_gdna > 0:
            gdna_frag_lengths = self._sample_gdna_frag_lengths(n_gdna)
            unique_lengths, length_counts = np.unique(
                gdna_frag_lengths, return_counts=True,
            )
            for frag_len, frag_count in zip(unique_lengths, length_counts):
                yield from self._gen_reads_from_genome(
                    int(frag_len), int(frag_count),
                )

    def write_fastq(
        self,
        output_dir: Path,
        n_fragments: int,
        *,
        prefix: str = "sim",
    ) -> tuple[Path, Path]:
        """Simulate reads and write to paired FASTQ files.

        Parameters
        ----------
        output_dir : Path
            Directory for output files.
        n_fragments : int
            Number of fragments to simulate.
        prefix : str
            Filename prefix (default ``"sim"``).

        Returns
        -------
        tuple[Path, Path]
            Paths to R1 and R2 FASTQ files.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        r1_path = output_dir / f"{prefix}_R1.fq"
        r2_path = output_dir / f"{prefix}_R2.fq"

        n_written = 0
        with open(r1_path, "w") as r1_fh, open(r2_path, "w") as r2_fh:
            for r1_name, r1_seq, r1_qual, r2_name, r2_seq, r2_qual in self.simulate(
                n_fragments
            ):
                r1_fh.write(f"@{r1_name}\n{r1_seq}\n+\n{r1_qual}\n")
                r2_fh.write(f"@{r2_name}\n{r2_seq}\n+\n{r2_qual}\n")
                n_written += 1

        logger.info(f"Wrote {n_written} read pairs → {r1_path}, {r2_path}")

        # Log the pool split
        n_mrna, n_nrna, n_gdna = self._compute_pool_split(n_fragments)
        parts = [f"mRNA: {n_mrna}"]
        if n_nrna > 0:
            parts.append(f"nRNA: {n_nrna} ({n_nrna / n_fragments:.1%})")
        if n_gdna > 0:
            parts.append(f"gDNA: {n_gdna} ({n_gdna / n_fragments:.1%})")
        if n_nrna > 0 or n_gdna > 0:
            logger.info(f"  {', '.join(parts)}")

        return r1_path, r2_path
