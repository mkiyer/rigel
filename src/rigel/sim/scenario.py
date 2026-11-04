"""``Scenario``: a small reproducible library end to end — random genome, annotation, index, reads, and either a
minimap2 alignment or the oracle BAM. Artifacts go to a working directory (temporary by default) that ``cleanup()``
or the context manager removes.

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
from dataclasses import dataclass
from pathlib import Path

import numpy as np


from ..index import TranscriptIndex
from ..transcript import Transcript
from ..types import Strand

from .annotation import GeneBuilder
from .capture import CaptureConfig
from .genome import MutableGenome
from .reads import GDNAConfig, ReadSimConfig
from .wgs_config import GDNASimConfig, SimulationParams
from .wgs_engine import WholeGenomeSimulator
from .whole_genome import assign_nrna_to_entities, load_transcripts_from_index
from .truth import (
    count_mrna_by_transcript_from_bam,
    count_mrna_by_transcript_from_fastq,
    count_origins_from_bam,
    count_origins_from_fastq,
)

logger = logging.getLogger(__name__)

__all__ = ["Scenario", "ScenarioResult"]


def _to_sim_params(sim_config: ReadSimConfig, n_rna: int) -> SimulationParams:
    """Translate a ``ReadSimConfig`` into the engine's ``SimulationParams``, single-worker. Strand specificity goes
    to ``WholeGenomeSimulator`` separately."""
    return SimulationParams(
        n_rna_fragments=n_rna,
        sim_seed=sim_config.seed,
        frag_mean=sim_config.frag_mean,
        frag_std=sim_config.frag_std,
        frag_min=sim_config.frag_min,
        frag_max=sim_config.frag_max,
        read_length=sim_config.read_length,
        error_rate=sim_config.error_rate,
        n_workers=1,
    )


def _gdna_overdispersion(gdna_config: GDNAConfig) -> float:
    """The gDNA strand overdispersion ``GDNAConfig`` sets; ``0`` when it is unset."""
    od = gdna_config.gdna_strand_overdispersion
    if od is not None:
        return float(od)
    return 0.0


def _to_gdna_sim(gdna_config: GDNAConfig | None) -> GDNASimConfig:
    """Translate a ``GDNAConfig`` into the engine's ``GDNASimConfig``: fragment lengths and strand overdispersion. The
    gDNA count comes from the pool split."""
    if gdna_config is None:
        return GDNASimConfig()
    return GDNASimConfig(
        frag_mean=gdna_config.frag_mean,
        frag_std=gdna_config.frag_std,
        frag_min=gdna_config.frag_min,
        frag_max=gdna_config.frag_max,
        strand_overdispersion=_gdna_overdispersion(gdna_config),
    )


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
        Name-sorted when aligned, query-grouped (``GO:query``) when oracle.
    index_dir : Path
        TranscriptIndex directory.
    index : TranscriptIndex
        The loaded index.
    transcripts : list[Transcript]
        The annotated transcripts, with index rows as ``t_index``.
    genome : MutableGenome
        The genome, splice motifs written.
    fastq_r1 : Path or None
        R1 FASTQ file (None in oracle mode).
    fastq_r2 : Path or None
        R2 FASTQ file (None in oracle mode).
    n_simulated : int
        Number of fragments passed to the read simulator.
    gdna_config : GDNAConfig or None
        gDNA configuration used (None if no gDNA simulated).
    is_oracle : bool
        True when the BAM is the simulator's oracle BAM, not an alignment.
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

    def ground_truth_gdna_count(self) -> int:
        """Count gDNA fragments from FASTQ read names."""
        return count_origins_from_fastq(self.fastq_r1)["gdna"]

    def ground_truth_nrna_count(self) -> int:
        """Count nascent RNA fragments from FASTQ read names."""
        return count_origins_from_fastq(self.fastq_r1)["nrna"]

    def ground_truth_from_fastq(self) -> dict[str, int]:
        """``t_id → count`` of the transcript-origin fragments the simulator wrote, from FASTQ read names, aligned or
        not; gDNA and nascent fragments are excluded."""
        return dict(count_mrna_by_transcript_from_fastq(self.fastq_r1))

    # -- BAM-based ground truth (oracle mode) ---------------------------------

    def ground_truth_from_bam(self) -> dict[str, int]:
        """``t_id → count`` of the transcript-origin fragments in the BAM, from its read names."""
        return dict(count_mrna_by_transcript_from_bam(self.bam_path))

    def ground_truth_gdna_count_from_bam(self) -> int:
        """Count gDNA fragments from BAM read names."""
        return count_origins_from_bam(self.bam_path)["gdna"]

    def ground_truth_nrna_count_from_bam(self) -> int:
        """Count nascent RNA fragments from BAM read names."""
        return count_origins_from_bam(self.bam_path)["nrna"]

    def ground_truth_auto(self) -> dict[str, int]:
        """Ground-truth transcript counts from the FASTQ when it exists, else from the BAM."""
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
        Output directory; ``None`` creates a temporary one that ``cleanup()`` removes.
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
        capture_config: CaptureConfig | None = None,
    ):
        self.name = name
        self.seed = seed
        self.gdna_config = gdna_config
        self.capture_config = capture_config
        self._owns_workdir = work_dir is None

        if work_dir is None:
            self.work_dir = Path(tempfile.mkdtemp(prefix=f"rigel_sim_{name}_"))
        else:
            self.work_dir = Path(work_dir)
            self.work_dir.mkdir(parents=True, exist_ok=True)

        self.ref_name = ref_name or name
        self.genome = MutableGenome(genome_length, seed=seed, name=self.ref_name)
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

    def _index_and_transcriptome(
        self, fasta_path: Path, gtf_path: Path, nrna_abundance: float
    ) -> tuple[Path, TranscriptIndex, list[Transcript], list[Transcript]]:
        """Build the index from the written GTF and return the transcriptome `rigel quant` will read, with the
        annotation's abundances carried onto it by ``t_id`` and each multi-exon transcript's nascent molecules
        pooled onto its entity. ``nrna_abundance`` is :meth:`build`'s.

        Returns ``(index_dir, index, all_rows, annotated_rows)``.
        """
        index_dir = self.work_dir / "index"
        TranscriptIndex.build(fasta_path, gtf_path, index_dir, write_tsv=False)
        rows = load_transcripts_from_index(index_dir)
        by_id = {t.t_id: t for t in self.annotation.get_transcripts()}
        per_contributor = np.zeros(len(rows))
        for i, row in enumerate(rows):
            src = by_id.get(row.t_id)
            if src is None:
                continue  # a synthetic entity: no mature molecules of its own
            row.abundance = float(src.abundance or 0.0)
            if len(row.exons) > 1 and row.abundance > 0:
                own = float(src.nrna_abundance or 0.0)
                per_contributor[i] = nrna_abundance if nrna_abundance > 0 else own
        assign_nrna_to_entities(rows, per_contributor)
        index = TranscriptIndex.load(index_dir)
        annotated = [t for t in rows if not t.is_synthetic]
        return index_dir, index, rows, annotated

    def build(
        self,
        n_fragments: int = 1000,
        sim_config: ReadSimConfig | None = None,
        gdna_config: GDNAConfig | None = None,
        capture_config: CaptureConfig | None = None,
        nrna_abundance: float = 0.0,
    ) -> ScenarioResult:
        """Write the genome and GTF, build the index, simulate paired-end FASTQ, and align it with minimap2 into a
        name-sorted BAM.

        Parameters
        ----------
        n_fragments : int
            Number of fragments to simulate.
        sim_config : ReadSimConfig or None
            Read simulation config; ``None`` uses the defaults at the scenario's seed.
        gdna_config : GDNAConfig or None
            gDNA contamination config; ``None`` uses the constructor's.
        nrna_abundance : float
            When > 0, the nascent molecular abundance of every expressed multi-exon transcript,
            independent of its ``abundance``; 0.0 (default) keeps the values from ``add_gene``.
            Reads are drawn in proportion to ``abundance × spliced_eff_len`` against
            ``nrna_abundance × genomic_span_eff_len``.

        Returns
        -------
        ScenarioResult
        """
        wdir = self.work_dir

        effective_gdna = gdna_config if gdna_config is not None else self.gdna_config
        effective_capture = capture_config if capture_config is not None else self.capture_config

        logger.info(f"[{self.name}] Writing genome FASTA...")
        fasta_path = self.genome.write_fasta(wdir)

        logger.info(f"[{self.name}] Writing GTF annotations...")
        gtf_path = self.annotation.write_gtf(wdir)
        logger.info(f"[{self.name}] Building TranscriptIndex...")
        index_dir, index, rows, transcripts = self._index_and_transcriptome(
            fasta_path, gtf_path, nrna_abundance
        )

        if sim_config is None:
            sim_config = ReadSimConfig(seed=self.seed)
        sim = WholeGenomeSimulator(
            fasta_path,
            rows,
            _to_sim_params(sim_config, n_fragments),
            _to_gdna_sim(effective_gdna),
            genomic_refs=[self.ref_name],
            strand_specificity=sim_config.strand_specificity,
            r1_sense=sim_config.r1_sense,
            capture_config=effective_capture,
        )
        gdna_abundance = effective_gdna.abundance if effective_gdna else 0.0
        n_rna, n_gdna = sim.rna_gdna_split(n_fragments, gdna_abundance)
        logger.info(f"[{self.name}] Simulating {n_fragments} fragments...")
        r1_path, r2_path, _ = sim.simulate_and_write(
            wdir,
            n_rna,
            n_gdna,
            oracle_bam=False,
            prefix=self.name,
        )
        sim.close()

        logger.info(f"[{self.name}] Aligning with minimap2...")
        bam_path = self._align(fasta_path, r1_path, r2_path, gtf_path=gtf_path)

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
        sim_config: ReadSimConfig | None = None,
        gdna_config: GDNAConfig | None = None,
        capture_config: CaptureConfig | None = None,
        nrna_abundance: float = 0.0,
        n_rna_fragments: int | None = None,
        gdna_fraction: float | None = None,
    ) -> ScenarioResult:
        """Like :meth:`build`, but with no aligner: the engine writes the oracle BAM directly, query-grouped, its
        CIGARs from the known template structure. The result carries no FASTQ.

        Parameters
        ----------
        n_fragments : int
            Fixed total, split between RNA and gDNA by abundance; ignored when *n_rna_fragments* is set.
        sim_config, gdna_config, nrna_abundance
            As for :meth:`build`.
        n_rna_fragments : int or None
            Fixed RNA fragment count, with ``n_rna_fragments × gdna_fraction`` gDNA fragments added.
        gdna_fraction : float or None
            gDNA fragments per RNA fragment, read only with *n_rna_fragments*.

        Returns
        -------
        ScenarioResult
            Result with ``is_oracle=True`` and ``fastq_r1/r2=None``.
        """
        wdir = self.work_dir
        effective_gdna = gdna_config if gdna_config is not None else self.gdna_config
        effective_capture = capture_config if capture_config is not None else self.capture_config

        fasta_path = self.genome.write_fasta(wdir)
        gtf_path = self.annotation.write_gtf(wdir)
        index_dir, index, rows, transcripts = self._index_and_transcriptome(
            fasta_path, gtf_path, nrna_abundance
        )

        if sim_config is None:
            sim_config = ReadSimConfig(seed=self.seed)
        gdna_abundance = effective_gdna.abundance if effective_gdna else 0.0
        if n_rna_fragments is not None:
            n_rna = n_rna_fragments
            n_gdna = int(round(n_rna * gdna_fraction)) if gdna_fraction and effective_gdna else 0
            sim = WholeGenomeSimulator(
                fasta_path,
                rows,
                _to_sim_params(sim_config, n_rna),
                _to_gdna_sim(effective_gdna),
                genomic_refs=[self.ref_name],
                strand_specificity=sim_config.strand_specificity,
                r1_sense=sim_config.r1_sense,
                capture_config=effective_capture,
            )
        else:
            sim = WholeGenomeSimulator(
                fasta_path,
                rows,
                _to_sim_params(sim_config, n_fragments),
                _to_gdna_sim(effective_gdna),
                genomic_refs=[self.ref_name],
                strand_specificity=sim_config.strand_specificity,
                r1_sense=sim_config.r1_sense,
                capture_config=effective_capture,
            )
            n_rna, n_gdna = sim.rna_gdna_split(n_fragments, gdna_abundance)
        _, _, bam_path = sim.simulate_and_write(
            wdir,
            n_rna,
            n_gdna,
            oracle_bam=True,
            prefix=self.name,
        )
        sim.close()
        n_simulated_total = n_rna + n_gdna

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
            n_simulated=n_simulated_total,
            gdna_config=effective_gdna,
            is_oracle=True,
        )

    def _align(
        self, fasta_path: Path, r1_path: Path, r2_path: Path, gtf_path: Path | None = None
    ) -> Path:
        """Align with ``minimap2 -ax splice:sr --secondary=yes [-j ann.bed] | samtools sort -n``, keeping secondary
        alignments, and return the name-sorted BAM. With *gtf_path*, the annotation's junctions go to minimap2 as
        BED12.
        """
        from ..index import write_bed12, read_transcripts

        bam_path = self.work_dir / f"{self.name}.bam"

        minimap2_cmd = [
            "minimap2",
            "-ax",
            "splice:sr",
            "--secondary=yes",
        ]

        if gtf_path is not None:
            bed_path = self.work_dir / f"{self.name}_annotation.bed"
            transcripts = read_transcripts(gtf_path)
            write_bed12(transcripts, bed_path)
            minimap2_cmd.extend(["-j", str(bed_path)])
            logger.info(f"Using BED12 annotation for minimap2: {bed_path}")

        minimap2_cmd.extend(
            [
                str(fasta_path),
                str(r1_path),
                str(r2_path),
            ]
        )

        sort_cmd = [
            "samtools",
            "sort",
            "-n",  # name-sort for the BAM parser
            "-o",
            str(bam_path),
            "-",  # read from stdin
        ]

        logger.debug(f"Running: {' '.join(minimap2_cmd)} | {' '.join(sort_cmd)}")

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
            raise RuntimeError(f"samtools sort failed (rc={p2.returncode}): {p2_stderr.decode()}")
        if p1.returncode != 0:
            logger.warning(
                "minimap2 exited with rc=%d: %s",
                p1.returncode,
                p1_stderr.decode(errors="replace"),
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
