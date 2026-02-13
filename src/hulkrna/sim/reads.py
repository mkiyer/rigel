"""
hulkrna.sim.reads — Paired-end RNA-seq read simulator.

Generates paired-end FASTQ reads from transcript sequences with
configurable fragment-size distribution, read length, and error rate.
Read names encode ground-truth origin for post-hoc validation.

Read-name format::

    {t_id}:{frag_start}-{frag_end}:{strand_char}:{index}/1
    {t_id}:{frag_start}-{frag_end}:{strand_char}:{index}/2

Library-prep convention (FR orientation)
----------------------------------------
Reads simulate a standard Illumina stranded (dUTP) library:

- **Positive-strand transcript**: R1 is reverse-complement of 3′ end
  of the fragment, R2 is the sense sequence from the 5′ end.
  When aligned, R1 maps to − strand and R2 maps to + strand.
  ``Fragment.from_reads`` flips R2's strand → both report + → correct.

- **Negative-strand transcript**: R1 maps to + strand, R2 maps to − strand.
  After R2 flip → both report − → correct.
"""

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..types import Strand
from ..transcript import Transcript
from .genome import MutableGenome, reverse_complement

logger = logging.getLogger(__name__)

__all__ = ["SimConfig", "ReadSimulator"]

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
    seed : int
        Random seed for reproducibility.
    """

    frag_mean: float = 250.0
    frag_std: float = 50.0
    frag_min: int = 50
    frag_max: int = 1000
    read_length: int = 150
    error_rate: float = 0.0
    seed: int = 42


class ReadSimulator:
    """Paired-end RNA-seq read simulator.

    Generates reads from transcript sequences extracted from the genome
    by concatenating exon coordinates. Fragment lengths are drawn from
    a truncated normal distribution. Transcript selection is weighted
    by ``abundance × effective_length``.

    Parameters
    ----------
    genome : MutableGenome
        The genome from which transcript sequences are extracted.
    transcripts : list[Transcript]
        Transcript annotations with exon coordinates and abundance.
    config : SimConfig or None
        Simulation parameters. ``None`` uses defaults.
    """

    # Oversampling ratio to reduce rejection-sampling iterations
    _OVERSAMPLE_RATIO = 1.5
    _OVERSAMPLE_EXTRA = 10

    def __init__(
        self,
        genome: MutableGenome,
        transcripts: list[Transcript],
        config: SimConfig | None = None,
    ):
        self.genome = genome
        self.config = config or SimConfig()
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

    # -- Fragment-length sampling ---------------------------------------------

    def _sample_frag_lengths(self, n: int) -> np.ndarray:
        """Sample *n* fragment lengths from a truncated normal."""
        cfg = self.config
        rng = self._rng
        result = np.empty(n, dtype=int)
        filled = 0

        while filled < n:
            needed = n - filled
            sample_size = int(needed * self._OVERSAMPLE_RATIO) + self._OVERSAMPLE_EXTRA
            raw = rng.normal(cfg.frag_mean, cfg.frag_std, sample_size).astype(int)
            valid = raw[(raw >= cfg.frag_min) & (raw <= cfg.frag_max)]
            nkeep = min(len(valid), needed)
            result[filled : filled + nkeep] = valid[:nkeep]
            filled += nkeep

        return result

    # -- Transcript probability -----------------------------------------------

    def _compute_probs(self, frag_length: int) -> np.ndarray:
        """Compute transcript selection probabilities for a given fragment length.

        Weight = max(0, transcript_length - frag_length + 1) × abundance
        """
        weights = np.array(
            [
                max(0, tlen - frag_length + 1) * t.abundance
                for tlen, t in zip(self._t_lengths, self.transcripts)
            ]
        )
        total = weights.sum()
        if total == 0:
            # All transcripts shorter than frag_length — uniform fallback
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

        strand_char = "r" if t.strand == Strand.NEG else "f"

        for i in range(count):
            frag_start = frag_starts[i]
            frag_end = frag_start + frag_len
            frag_seq = t_seq[frag_start:frag_end]

            # FR library: R1 is reverse-complement of 3′ end,
            # R2 is sense sequence from 5′ end.
            r1_seq = reverse_complement(frag_seq[-read_len:])
            r2_seq = frag_seq[:read_len]

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

    # -- Public interface -----------------------------------------------------

    def simulate(self, n_fragments: int):
        """Generate *n_fragments* paired-end read pairs.

        Yields ``(r1_name, r1_seq, r1_qual, r2_name, r2_seq, r2_qual)``.
        """
        rng = self._rng
        ntranscripts = len(self.transcripts)

        # Sample fragment lengths
        frag_lengths = self._sample_frag_lengths(n_fragments)

        # Group by unique fragment length for batched transcript selection
        unique_lengths, length_counts = np.unique(frag_lengths, return_counts=True)

        for frag_len, frag_count in zip(unique_lengths, length_counts):
            probs = self._compute_probs(int(frag_len))

            # Sample transcripts for this fragment length
            t_indices = rng.choice(ntranscripts, size=int(frag_count), p=probs)
            unique_t, t_counts = np.unique(t_indices, return_counts=True)

            for t_idx, t_count in zip(unique_t, t_counts):
                yield from self._gen_reads_from_transcript(
                    int(t_idx), int(frag_len), int(t_count)
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
        return r1_path, r2_path
