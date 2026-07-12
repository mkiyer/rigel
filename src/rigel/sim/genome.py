"""
rigel.sim.genome — Mutable genome sequence for simulation.

Generates a random DNA sequence of configurable length using a seeded
RNG, stores it as a mutable list for O(1) edits (e.g., injecting
splice-site motifs), and writes FASTA + ``samtools faidx`` index.
"""

import textwrap
from pathlib import Path

import numpy as np
import pysam

__all__ = ["MutableGenome", "reverse_complement", "random_dna_array", "write_fasta_file"]

# DNA complement table for reverse-complement
_DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Integer encoding
_DNA_BASES = np.array(list("ACGT"), dtype="<U1")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(_DNA_COMPLEMENT)[::-1]


# ── Shared genome-building mechanics ────────────────────────────────────────
# The single implementations behind BOTH synthetic-reference paths — the class path
# (``MutableGenome`` + ``GeneBuilder``, used by ``Scenario``) and the function path
# (``synthetic_genome``, used by the whole-genome suite + reference-gen scripts). The two
# paths keep distinct OUTPUT contracts (filenames + GTF format) for two ecosystems, but
# share these low-level mechanics so a fix reaches both.


def random_dna_array(length: int, rng: np.random.Generator) -> np.ndarray:
    """Random DNA as a numpy array of single-char strings (``ACGT``), drawn from *rng*."""
    return _DNA_BASES[rng.integers(0, 4, size=length)]


def write_fasta_file(seq: str, ref_name: str, fasta_path: Path) -> Path:
    """Write *seq* as a single-record FASTA (80-column wrap) at *fasta_path* + ``samtools faidx``.

    Parameters
    ----------
    seq : str
        Genome sequence.
    ref_name : str
        FASTA header / reference name.
    fasta_path : Path
        Destination ``.fa`` path (parent dir created if needed). The two paths differ only in
        the filename they choose (the suite writes ``genome.fa``; a ``MutableGenome`` writes
        ``{name}.fa``), so the caller passes the full path.
    """
    fasta_path = Path(fasta_path)
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with open(fasta_path, "w") as f:
        f.write(f">{ref_name}\n")
        f.write(textwrap.fill(seq, width=80) + "\n")
    pysam.faidx(str(fasta_path))
    return fasta_path


class MutableGenome:
    """Random DNA genome with O(1) positional editing.

    Parameters
    ----------
    length : int
        Number of bases in the genome.
    seed : int
        Seed for the random number generator (reproducibility).
    name : str
        Reference name for the FASTA header.

    Examples
    --------
    >>> g = MutableGenome(1000, seed=42, name="chr1")
    >>> g.edit(100, "GT")       # inject donor motif
    >>> g.edit(298, "AG")       # inject acceptor motif
    >>> path = g.write_fasta(Path("/tmp/test"))
    """

    def __init__(self, length: int, *, seed: int = 42, name: str = "chr1"):
        self.name = name
        self._rng = np.random.default_rng(seed)
        # Generate random DNA as a mutable list (shared generator — see random_dna_array).
        self._seq: list[str] = list(random_dna_array(length, self._rng))

    # -- Properties -----------------------------------------------------------

    @property
    def seq(self) -> str:
        """Full genome sequence as an immutable string."""
        return "".join(self._seq)

    def __len__(self) -> int:
        return len(self._seq)

    def __getitem__(self, key) -> str:
        """Return a substring. Supports int and slice indexing."""
        if isinstance(key, slice):
            return "".join(self._seq[key])
        return self._seq[key]

    # -- Editing --------------------------------------------------------------

    def edit(self, pos: int, bases: str) -> None:
        """Overwrite bases at *pos* with the given string.

        Parameters
        ----------
        pos : int
            0-based start position.
        bases : str
            Replacement bases (e.g., ``"GT"`` for a donor motif).

        Raises
        ------
        IndexError
            If the edit extends beyond the genome boundary.
        """
        end = pos + len(bases)
        if end > len(self._seq) or pos < 0:
            raise IndexError(f"Edit at {pos}:{end} exceeds genome length {len(self._seq)}")
        for i, b in enumerate(bases):
            self._seq[pos + i] = b

    # -- I/O ------------------------------------------------------------------

    def write_fasta(self, output_dir: Path) -> Path:
        """Write genome as FASTA and create ``samtools faidx`` index.

        Parameters
        ----------
        output_dir : Path
            Directory to write ``{name}.fa`` into (created if needed).

        Returns
        -------
        Path
            Path to the written FASTA file.
        """
        return write_fasta_file(self.seq, self.name, Path(output_dir) / f"{self.name}.fa")
