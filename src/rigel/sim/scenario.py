"""
rigel.sim.scenario — End-to-end simulation scenario orchestrator.

A ``Scenario`` bundles genome generation, gene annotation, read
simulation, alignment (minimap2), BAM collation (samtools), and
index building (TranscriptIndex) into a single reproducible workflow.

All artifacts are written to a working directory (temp by default)
and can be cleaned up via ``cleanup()`` or the context-manager
protocol.

Usage
-----
>>> with Scenario("test1", genome_length=5000, seed=42) as sc:
...     sc.add_gene("g1", "+", [
...         {"t_id": "t1", "exons": [(100, 300), (500, 700)], "abundance": 100},
...     ])
...     result = sc.build(n_fragments=500)
...     # result.bam_path, result.index, result.transcripts are ready
"""

import logging
import shutil
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pysam

from ..index import TranscriptIndex
from ..transcript import Transcript
from ..types import Strand

from .annotation import GeneBuilder
from .genome import MutableGenome
from .oracle_bam import OracleBamSimulator
from .reads import GDNAConfig, ReadSimulator, SimConfig

logger = logging.getLogger(__name__)

__all__ = ["Scenario", "ScenarioResult"]


@dataclass
class ScenarioResult:
    """Output artifacts from a completed simulation scenario.

    Attributes
    ----------
    fasta_path : Path
        Genome FASTA file (with .fai index).
    gtf_path : Path
        Gene annotation GTF file.
    bam_path : Path
        Name-sorted BAM file (aligned or oracle).
    index_dir : Path
        TranscriptIndex output directory (Feather files).
    index : TranscriptIndex
        Loaded TranscriptIndex ready for quantification.
    transcripts : list[Transcript]
        Transcript objects with assigned indices.
    genome : MutableGenome
        The genome (post splice-motif edits).
    fastq_r1 : Path or None
        R1 FASTQ file (None in oracle mode).
    fastq_r2 : Path or None
        R2 FASTQ file (None in oracle mode).
    n_simulated : int
        Number of fragments passed to the read simulator.
    gdna_config : GDNAConfig or None
        gDNA configuration used (None if no gDNA simulated).
    is_oracle : bool
        True when BAM was produced by the oracle simulator
        (no FASTQ/alignment step).
    """

    fasta_path: Path
    gtf_path: Path
    bam_path: Path
    index_dir: Path
    index: TranscriptIndex
    transcripts: list[Transcript]
    genome: MutableGenome
    fastq_r1: Path | None = None
    fastq_r2: Path | None = None
    n_simulated: int = 0
    gdna_config: GDNAConfig | None = None
    is_oracle: bool = False

    def ground_truth_counts(self) -> dict[str, int]:
        """Parse BAM read names to extract ground-truth fragment counts.

        The read simulator encodes the source transcript in each read
        name as ``{t_id}:{frag_start}-{frag_end}:{strand}:{idx}/1``.
        gDNA reads use ``gdna:...`` prefix and are excluded.

        Returns
        -------
        dict[str, int]
            Mapping of ``t_id → count`` for each transcript that
            contributed at least one simulated fragment.
        """
        seen: set[str] = set()
        with pysam.AlignmentFile(str(self.bam_path),
                                 "rb", check_sq=False) as bam:
            for read in bam:
                # Strip /1 or /2 suffix to get fragment identifier
                qname = read.query_name
                if qname not in seen:
                    seen.add(qname)

        counts: Counter[str] = Counter()
        for frag_id in seen:
            # Format: t_id:start-end:strand:idx  (or gdna:..., nrna_*:...)
            t_id = frag_id.split(":")[0]
            if t_id != "gdna" and not t_id.startswith("nrna_"):
                counts[t_id] += 1
        return dict(counts)

    def ground_truth_gdna_count(self) -> int:
        """Count gDNA fragments from FASTQ read names.

        Returns
        -------
        int
            Number of fragments with ``gdna:`` prefix.
        """
        count = 0
        with open(self.fastq_r1) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    if line.startswith("@gdna:"):
                        count += 1
        return count

    def ground_truth_nrna_count(self) -> int:
        """Count nascent RNA fragments from FASTQ read names.

        Returns
        -------
        int
            Number of fragments with ``nrna_`` prefix.
        """
        count = 0
        with open(self.fastq_r1) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    if line.startswith("@nrna_"):
                        count += 1
        return count

    def ground_truth_from_fastq(self) -> dict[str, int]:
        """Parse FASTQ read names to get ground-truth fragment counts.

        Unlike :meth:`ground_truth_counts`, this counts *all* simulated
        fragments regardless of alignment success, providing the true
        number of fragments the simulator produced per transcript.
        gDNA reads (``gdna:`` prefix) are excluded.

        Returns
        -------
        dict[str, int]
            Mapping of ``t_id → count`` for each transcript.
        """
        counts: Counter[str] = Counter()
        with open(self.fastq_r1) as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:  # header line
                    # Format: @t_id:start-end:strand:idx/1
                    qname = line[1:].strip()  # drop @
                    if qname.endswith("/1"):
                        qname = qname[:-2]
                    t_id = qname.split(":")[0]
                    if t_id != "gdna" and not t_id.startswith("nrna_"):
                        counts[t_id] += 1
        return dict(counts)

    # -- BAM-based ground truth (oracle mode) ---------------------------------

    def _parse_bam_qnames(self) -> list[str]:
        """Return deduplicated query names from the BAM (one per fragment)."""
        seen: set[str] = set()
        with pysam.AlignmentFile(str(self.bam_path), "rb",
                                 check_sq=False) as bam:
            for read in bam:
                seen.add(read.query_name)
        return list(seen)

    def ground_truth_from_bam(self) -> dict[str, int]:
        """Parse BAM read names for ground-truth mRNA fragment counts.

        Read-name format: ``{t_id}:{start}-{end}:{strand}:{idx}``

        Returns
        -------
        dict[str, int]
            Mapping of ``t_id → count`` for each transcript.
        """
        counts: Counter[str] = Counter()
        for qname in self._parse_bam_qnames():
            t_id = qname.split(":")[0]
            if t_id != "gdna" and not t_id.startswith("nrna_"):
                counts[t_id] += 1
        return dict(counts)

    def ground_truth_gdna_count_from_bam(self) -> int:
        """Count gDNA fragments from BAM read names."""
        return sum(
            1 for qn in self._parse_bam_qnames() if qn.startswith("gdna:")
        )

    def ground_truth_nrna_count_from_bam(self) -> int:
        """Count nascent RNA fragments from BAM read names."""
        return sum(
            1 for qn in self._parse_bam_qnames()
            if qn.split(":")[0].startswith("nrna_")
        )

    def ground_truth_auto(self) -> dict[str, int]:
        """Return ground-truth mRNA counts from FASTQ or BAM.

        Uses FASTQ-based ground truth if FASTQ files exist,
        otherwise falls back to BAM read names (oracle mode).
        """
        if self.fastq_r1 is not None and self.fastq_r1.exists():
            return self.ground_truth_from_fastq()
        return self.ground_truth_from_bam()

    def ground_truth_gdna_auto(self) -> int:
        """Return ground-truth gDNA count from FASTQ or BAM."""
        if self.fastq_r1 is not None and self.fastq_r1.exists():
            return self.ground_truth_gdna_count()
        return self.ground_truth_gdna_count_from_bam()

    def ground_truth_nrna_auto(self) -> int:
        """Return ground-truth nRNA count from FASTQ or BAM."""
        if self.fastq_r1 is not None and self.fastq_r1.exists():
            return self.ground_truth_nrna_count()
        return self.ground_truth_nrna_count_from_bam()


class Scenario:
    """End-to-end simulation scenario orchestrator.

    Parameters
    ----------
    name : str
        Scenario name (used for file naming).
    genome_length : int
        Length of the random genome (default 5000).
    seed : int
        Random seed for genome and read generation.
    work_dir : Path or None
        Working directory for output files.  If None, a temporary
        directory is created and cleaned up on ``cleanup()``.
    ref_name : str or None
        Reference name. Defaults to *name*.
    """

    def __init__(
        self,
        name: str,
        *,
        genome_length: int = 5000,
        seed: int = 42,
        work_dir: Path | None = None,
        ref_name: str | None = None,
        gdna_config: GDNAConfig | None = None,
    ):
        self.name = name
        self.seed = seed
        self.gdna_config = gdna_config
        self._owns_workdir = work_dir is None

        if work_dir is None:
            self.work_dir = Path(tempfile.mkdtemp(prefix=f"rigel_sim_{name}_"))
        else:
            self.work_dir = Path(work_dir)
            self.work_dir.mkdir(parents=True, exist_ok=True)

        self.ref_name = ref_name or name
        self.genome = MutableGenome(
            genome_length, seed=seed, name=self.ref_name
        )
        self.annotation = GeneBuilder(self.genome, ref_name=self.ref_name)

    # -- Gene addition (delegates to GeneBuilder) ----------------------------

    def add_gene(
        self,
        gene_id: str,
        strand: str | Strand,
        transcripts: list[dict],
        **kwargs,
    ) -> None:
        """Add a gene to the scenario. See ``GeneBuilder.add_gene``."""
        self.annotation.add_gene(gene_id, strand, transcripts, **kwargs)

    # -- Build ----------------------------------------------------------------

    def build(
        self,
        n_fragments: int = 1000,
        sim_config: SimConfig | None = None,
        gdna_config: GDNAConfig | None = None,
        nrna_abundance: float = 0.0,
    ) -> ScenarioResult:
        """Execute the full simulation pipeline.

        1. Write genome FASTA + index.
        2. Write gene annotations (GTF).
        3. Simulate paired-end reads (FASTQ).
        4. Align with minimap2 (``splice:sr``).
        5. Collate BAM with samtools.
        6. Build TranscriptIndex.

        Parameters
        ----------
        n_fragments : int
            Number of fragments to simulate.
        sim_config : SimConfig or None
            Read simulation config. None uses defaults with same seed.
        gdna_config : GDNAConfig or None
            gDNA contamination config.  If None, falls back to the
            ``gdna_config`` set on the ``Scenario`` constructor.
        nrna_abundance : float
            Nascent RNA molecular abundance applied to all expressed
            transcripts.  This is an **independent** abundance value —
            ``abundance`` (mature mRNA) is NOT modified.  The read pool
            split is proportional to abundance × effective_length for
            each pool.

            When set to 0.0 (default), per-transcript ``nrna_abundance``
            values from the transcript dict (set at ``add_gene`` time)
            are preserved unchanged.  When > 0, overrides all expressed
            transcripts' ``nrna_abundance`` with this value.

            Example: ``nrna_abundance=30`` with ``abundance=100`` gives
            reads drawn proportionally to ``100 × spliced_eff_len``
            (mature) vs. ``30 × genomic_span_eff_len`` (nascent RNA).

        Returns
        -------
        ScenarioResult
        """
        wdir = self.work_dir

        # Resolve gDNA config: build() param overrides constructor
        effective_gdna = gdna_config if gdna_config is not None else self.gdna_config

        # 1. Genome FASTA
        logger.info(f"[{self.name}] Writing genome FASTA...")
        fasta_path = self.genome.write_fasta(wdir)

        # 2. GTF annotation
        logger.info(f"[{self.name}] Writing GTF annotations...")
        gtf_path = self.annotation.write_gtf(wdir)
        transcripts = self.annotation.get_transcripts()

        # Apply nascent RNA abundance to all expressed transcripts.
        # nrna_abundance is set independently of abundance (mature mRNA).
        # When nrna_abundance > 0, it overrides per-transcript values.
        # When nrna_abundance == 0 (default), per-transcript values from
        # add_gene() dicts are preserved.
        if nrna_abundance > 0:
            for t in transcripts:
                if t.abundance and t.abundance > 0:
                    t.nrna_abundance = nrna_abundance

        # 3. Simulate reads
        if sim_config is None:
            sim_config = SimConfig(seed=self.seed)
        # Clamp frag_max to genome/transcript length
        sim = ReadSimulator(
            self.genome, transcripts,
            config=sim_config,
            gdna_config=effective_gdna,
        )
        logger.info(f"[{self.name}] Simulating {n_fragments} fragments...")
        r1_path, r2_path = sim.write_fastq(wdir, n_fragments)

        # 4. Align with minimap2 → SAM → name-sorted BAM
        logger.info(f"[{self.name}] Aligning with minimap2...")
        bam_path = self._align(fasta_path, r1_path, r2_path, gtf_path=gtf_path)

        # 5. Build TranscriptIndex
        logger.info(f"[{self.name}] Building TranscriptIndex...")
        index_dir = wdir / "index"
        TranscriptIndex.build(fasta_path, gtf_path, index_dir, write_tsv=False)
        index = TranscriptIndex.load(index_dir)

        logger.info(f"[{self.name}] Scenario build complete.")
        return ScenarioResult(
            fasta_path=fasta_path,
            gtf_path=gtf_path,
            bam_path=bam_path,
            index_dir=index_dir,
            index=index,
            transcripts=transcripts,
            genome=self.genome,
            fastq_r1=r1_path,
            fastq_r2=r2_path,
            n_simulated=n_fragments,
            gdna_config=effective_gdna,
        )

    def build_oracle(
        self,
        n_fragments: int = 1000,
        sim_config: SimConfig | None = None,
        gdna_config: GDNAConfig | None = None,
        nrna_abundance: float = 0.0,
        n_rna_fragments: int | None = None,
        gdna_fraction: float | None = None,
    ) -> ScenarioResult:
        """Execute the simulation pipeline using oracle (perfect) alignments.

        Bypasses FASTQ generation and external alignment (minimap2/samtools).
        Instead uses :class:`OracleBamSimulator` to produce a name-sorted
        BAM with CIGAR strings derived from known transcript structure.

        Pipeline:
        1. Write genome FASTA + index.
        2. Write gene annotations (GTF).
        3. Generate oracle BAM (perfect alignments).
        4. Build TranscriptIndex.

        Parameters
        ----------
        n_fragments : int
            Number of fragments to simulate (mode a: fixed total).
            Ignored when *n_rna_fragments* is set.
        sim_config : SimConfig or None
            Read simulation config. None uses defaults with same seed.
        gdna_config : GDNAConfig or None
            gDNA contamination config.  If None, falls back to the
            ``gdna_config`` set on the ``Scenario`` constructor.
        nrna_abundance : float
            Nascent RNA molecular abundance applied to all expressed
            transcripts (same semantics as ``build()``).
        n_rna_fragments : int or None
            Fixed RNA fragment count (mode b).  When set, exactly this
            many RNA fragments (mRNA + nRNA) are generated and gDNA
            is added on top as ``n_rna_fragments × gdna_fraction``.
        gdna_fraction : float or None
            Fraction of RNA fragments to add as gDNA contamination
            (only used when *n_rna_fragments* is set).

        Returns
        -------
        ScenarioResult
            Result with ``is_oracle=True`` and ``fastq_r1/r2=None``.
        """
        wdir = self.work_dir
        effective_gdna = gdna_config if gdna_config is not None else self.gdna_config

        # 1. Genome FASTA
        fasta_path = self.genome.write_fasta(wdir)

        # 2. GTF annotation
        gtf_path = self.annotation.write_gtf(wdir)
        transcripts = self.annotation.get_transcripts()

        # Apply nascent RNA abundance
        if nrna_abundance > 0:
            for t in transcripts:
                if t.abundance and t.abundance > 0:
                    t.nrna_abundance = nrna_abundance

        # 3. Oracle BAM (perfect alignments, no FASTQ or alignment step)
        if sim_config is None:
            sim_config = SimConfig(seed=self.seed)
        oracle = OracleBamSimulator(
            self.genome, transcripts,
            config=sim_config,
            gdna_config=effective_gdna,
            ref_name=self.ref_name,
        )
        bam_path = wdir / f"{self.name}_oracle.bam"

        # Fixed-RNA mode (b): compute explicit pool split
        pool_split = None
        if n_rna_fragments is not None:
            n_rna = n_rna_fragments
            n_gdna = (int(round(n_rna * gdna_fraction))
                      if gdna_fraction and effective_gdna else 0)
            n_mrna, n_nrna = oracle._sim.compute_rna_split(n_rna)
            pool_split = (n_mrna, n_nrna, n_gdna)

        oracle.write_bam(bam_path, n_fragments, name_sorted=True,
                         pool_split=pool_split)

        # 4. Build TranscriptIndex
        index_dir = wdir / "index"
        TranscriptIndex.build(fasta_path, gtf_path, index_dir, write_tsv=False)
        index = TranscriptIndex.load(index_dir)

        return ScenarioResult(
            fasta_path=fasta_path,
            gtf_path=gtf_path,
            bam_path=bam_path,
            index_dir=index_dir,
            index=index,
            transcripts=transcripts,
            genome=self.genome,
            fastq_r1=None,
            fastq_r2=None,
            n_simulated=sum(pool_split) if pool_split else n_fragments,
            gdna_config=effective_gdna,
            is_oracle=True,
        )

    def _align(self, fasta_path: Path, r1_path: Path, r2_path: Path,
               gtf_path: Path | None = None) -> Path:
        """Align reads with minimap2 and produce a name-sorted BAM.

        Pipeline::

            minimap2 -ax splice:sr --secondary=yes [-j ann.bed] \\
              ref.fa r1.fq r2.fq \\
              | samtools sort -n -o output.bam -

        The ``splice:sr`` preset handles short RNA-seq reads with
        splice-aware alignment.  ``--secondary=yes`` retains
        multimappers so that NH tags are populated.

        When *gtf_path* is provided, a BED12 annotation file is
        generated and passed via ``-j`` so that minimap2 can use
        known splice junctions to improve alignment accuracy.
        """
        from ..index import write_bed12, read_transcripts

        bam_path = self.work_dir / f"{self.name}.bam"

        # minimap2 → samtools sort -n (name-sorted)
        minimap2_cmd = [
            "minimap2",
            "-ax", "splice:sr",
            "--secondary=yes",
        ]

        # Add BED12 annotation if GTF is available
        if gtf_path is not None:
            bed_path = self.work_dir / f"{self.name}_annotation.bed"
            transcripts = read_transcripts(gtf_path)
            write_bed12(transcripts, bed_path)
            minimap2_cmd.extend(["-j", str(bed_path)])
            logger.info(f"Using BED12 annotation for minimap2: {bed_path}")

        minimap2_cmd.extend([
            str(fasta_path),
            str(r1_path),
            str(r2_path),
        ])

        sort_cmd = [
            "samtools", "sort",
            "-n",  # name-sort for the BAM parser
            "-o", str(bam_path),
            "-",   # read from stdin
        ]

        logger.debug(f"Running: {' '.join(minimap2_cmd)} | {' '.join(sort_cmd)}")

        # Pipe minimap2 output to samtools sort
        p1 = subprocess.Popen(
            minimap2_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p2 = subprocess.Popen(
            sort_cmd,
            stdin=p1.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Allow p1 to receive SIGPIPE if p2 exits early
        p1.stdout.close()

        _, p2_stderr = p2.communicate()

        # Drain p1 stderr to prevent deadlock (pipe buffer full)
        p1_stderr = p1.stderr.read()
        p1.wait()

        if p2.returncode != 0:
            raise RuntimeError(
                f"samtools sort failed (rc={p2.returncode}): "
                f"{p2_stderr.decode()}"
            )
        if p1.returncode != 0:
            logger.warning(
                "minimap2 exited with rc=%d: %s",
                p1.returncode, p1_stderr.decode(errors='replace'),
            )

        logger.info(f"Wrote name-sorted BAM → {bam_path}")
        return bam_path

    # -- Cleanup --------------------------------------------------------------

    def cleanup(self) -> None:
        """Remove the working directory if it was auto-created."""
        if self._owns_workdir and self.work_dir.exists():
            shutil.rmtree(self.work_dir, ignore_errors=True)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.cleanup()
