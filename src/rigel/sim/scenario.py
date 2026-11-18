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

from ..config import CONSTANTS
from ..index import TranscriptIndex
from ..transcript import Transcript
from ..types import Strand

from .annotation import GeneBuilder
from .genome import MutableGenome
from .reads import GDNAConfig, ReadSimConfig
from .transcriptome import assign_nrna_to_entities, load_transcripts_from_index
from .wgs_config import GDNASimConfig, SimulationParams
from .wgs_engine import WholeGenomeSimulator

logger = logging.getLogger(__name__)

__all__ = ["Scenario", "ScenarioResult"]


def _to_sim_params(sim_config: ReadSimConfig) -> SimulationParams:
    """Translate a ``ReadSimConfig`` into the engine's ``SimulationParams``. Strand specificity and the protocol
    direction go to ``WholeGenomeSimulator`` separately, and the fragment counts to ``simulate_and_write``."""
    return SimulationParams(
        sim_seed=sim_config.seed,
        frag_mean=sim_config.frag_mean,
        frag_std=sim_config.frag_std,
        frag_min=sim_config.frag_min,
        frag_max=sim_config.frag_max,
        read_length=sim_config.read_length,
        error_rate=sim_config.error_rate,
    )


def to_gdna_sim(gdna_config: GDNAConfig | None) -> GDNASimConfig:
    """Translate a ``GDNAConfig`` into the engine's ``GDNASimConfig``: fragment lengths and strand overdispersion. The
    gDNA count comes from the pool split."""
    if gdna_config is None:
        return GDNASimConfig()
    return GDNASimConfig(
        frag_mean=gdna_config.frag_mean,
        frag_std=gdna_config.frag_std,
        frag_min=gdna_config.frag_min,
        frag_max=gdna_config.frag_max,
        strand_overdispersion=gdna_config.strand_overdispersion,
    )


@dataclass
class ScenarioResult:
    """What a scenario build wrote.

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
    n_simulated : int
        Number of fragments passed to the read simulator.
    fastq_r1, fastq_r2 : Path or None
        The simulated reads :meth:`Scenario.build` aligned; ``None`` from :meth:`Scenario.build_oracle`.
    """

    fasta_path: Path
    gtf_path: Path
    bam_path: Path
    index_dir: Path
    index: TranscriptIndex
    transcripts: list[Transcript]
    n_simulated: int
    fastq_r1: Path | None = None
    fastq_r2: Path | None = None


class Scenario:
    """End-to-end simulation scenario orchestrator.

    Parameters
    ----------
    name : str
        Scenario name (used for file naming).
    genome_length : int
        Length of the random genome (bp).
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
        genome_length: int = CONSTANTS.simulator.scenario_genome_length,
        seed: int = CONSTANTS.simulator.scenario_seed,
        work_dir: Path | None = None,
        ref_name: str | None = None,
    ):
        self.name = name
        self.seed = seed
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

    def build(
        self,
        n_fragments: int = CONSTANTS.simulator.scenario_fragments,
        sim_config: ReadSimConfig | None = None,
        gdna_config: GDNAConfig | None = None,
        nrna_abundance: float = 0.0,
    ) -> ScenarioResult:
        """Write the genome and GTF, build the index, simulate paired-end FASTQ, and align it with minimap2 into a
        name-sorted BAM.

        Parameters
        ----------
        n_fragments : int
            Number of fragments to simulate, split between RNA and gDNA by abundance.
        sim_config : ReadSimConfig or None
            Read simulation config; ``None`` uses the defaults at the scenario's seed.
        gdna_config : GDNAConfig or None
            gDNA contamination config; ``None`` simulates no gDNA.
        nrna_abundance : float
            When > 0, the nascent molecular abundance of every expressed multi-exon transcript,
            independent of its ``abundance``; 0.0 (default) keeps the values from ``add_gene``.
            Reads are drawn in proportion to ``abundance × spliced_eff_len`` against
            ``nrna_abundance × genomic_span_eff_len``.

        Returns
        -------
        ScenarioResult
        """
        fasta_path, gtf_path, index_dir, index, rows = self._reference(nrna_abundance)
        sim = self._simulator(fasta_path, rows, sim_config, gdna_config)
        n_rna, n_gdna = sim.rna_gdna_split(n_fragments, _gdna_abundance(gdna_config))
        logger.info(f"[{self.name}] Simulating {n_fragments} fragments...")
        r1_path, r2_path, _ = self._write_reads(sim, n_rna, n_gdna, oracle_bam=False)

        logger.info(f"[{self.name}] Aligning with minimap2...")
        bam_path = self._align(fasta_path, r1_path, r2_path, gtf_path=gtf_path)

        logger.info(f"[{self.name}] Scenario build complete.")
        return ScenarioResult(
            fasta_path=fasta_path,
            gtf_path=gtf_path,
            bam_path=bam_path,
            index_dir=index_dir,
            index=index,
            transcripts=_annotated(rows),
            n_simulated=n_fragments,
            fastq_r1=r1_path,
            fastq_r2=r2_path,
        )

    def build_oracle(
        self,
        n_fragments: int = CONSTANTS.simulator.scenario_fragments,
        sim_config: ReadSimConfig | None = None,
        gdna_config: GDNAConfig | None = None,
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
            gDNA fragments per RNA fragment, read only with *n_rna_fragments* and a *gdna_config*, and refused
            without them.

        Returns
        -------
        ScenarioResult
            Result with ``fastq_r1/r2=None``.
        """
        _check_gdna_fraction(gdna_fraction, n_rna_fragments, gdna_config)
        fasta_path, gtf_path, index_dir, index, rows = self._reference(nrna_abundance)
        sim = self._simulator(fasta_path, rows, sim_config, gdna_config)
        if n_rna_fragments is not None:
            n_rna = n_rna_fragments
            n_gdna = int(round(n_rna * (gdna_fraction or 0.0)))
        else:
            n_rna, n_gdna = sim.rna_gdna_split(n_fragments, _gdna_abundance(gdna_config))
        _, _, bam_path = self._write_reads(sim, n_rna, n_gdna, oracle_bam=True)
        return ScenarioResult(
            fasta_path=fasta_path,
            gtf_path=gtf_path,
            bam_path=bam_path,
            index_dir=index_dir,
            index=index,
            transcripts=_annotated(rows),
            n_simulated=n_rna + n_gdna,
        )

    # -- Build steps ------------------------------------------------------------

    def _reference(
        self, nrna_abundance: float
    ) -> tuple[Path, Path, Path, TranscriptIndex, list[Transcript]]:
        """Write the genome and the GTF, build the index, and return the transcriptome `rigel quant` will read, with
        the annotation's abundances carried onto it by ``t_id`` and each multi-exon transcript's nascent molecules
        pooled onto its entity. ``nrna_abundance`` is :meth:`build`'s.

        Returns ``(fasta_path, gtf_path, index_dir, index, rows)``, ``rows`` every index row.
        """
        logger.info(f"[{self.name}] Writing genome FASTA...")
        fasta_path = self.genome.write_fasta(self.work_dir)
        logger.info(f"[{self.name}] Writing GTF annotations...")
        gtf_path = self.annotation.write_gtf(self.work_dir)
        logger.info(f"[{self.name}] Building TranscriptIndex...")
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
            if len(row.exons) > 1:
                # a transcript's own nascent abundance is independent of its mature one; the build's override names
                # the nascent abundance of every expressed multi-exon transcript
                own = float(src.nrna_abundance or 0.0)
                per_contributor[i] = (
                    (nrna_abundance if row.abundance > 0 else 0.0) if nrna_abundance > 0 else own
                )
        assign_nrna_to_entities(rows, per_contributor)
        return fasta_path, gtf_path, index_dir, TranscriptIndex.load(index_dir), rows

    def _simulator(
        self,
        fasta_path: Path,
        rows: list[Transcript],
        sim_config: ReadSimConfig | None,
        gdna_config: GDNAConfig | None,
    ) -> WholeGenomeSimulator:
        """The engine over every index row, gDNA drawn from this scenario's one reference, single-worker."""
        if sim_config is None:
            sim_config = ReadSimConfig(seed=self.seed)
        return WholeGenomeSimulator(
            fasta_path,
            rows,
            _to_sim_params(sim_config),
            to_gdna_sim(gdna_config),
            genomic_refs=[self.ref_name],
            strand_specificity=sim_config.strand_specificity,
            r1_sense=sim_config.r1_sense,
        )

    def _write_reads(
        self, sim: WholeGenomeSimulator, n_rna: int, n_gdna: int, *, oracle_bam: bool
    ) -> tuple[Path, Path, Path | None]:
        try:
            return sim.simulate_and_write(
                self.work_dir, n_rna, n_gdna, oracle_bam=oracle_bam, prefix=self.name
            )
        finally:
            sim.close()

    def _align(self, fasta_path: Path, r1_path: Path, r2_path: Path, gtf_path: Path) -> Path:
        """Align with ``minimap2 -ax splice:sr --secondary=yes -j ann.bed | samtools sort -n``, keeping secondary
        alignments and giving minimap2 the annotation's junctions as BED12, and return the name-sorted BAM.
        """
        from ..index import read_transcripts, write_bed12

        bam_path = self.work_dir / f"{self.name}.bam"
        bed_path = self.work_dir / f"{self.name}_annotation.bed"
        write_bed12(read_transcripts(gtf_path), bed_path)
        logger.info(f"Using BED12 annotation for minimap2: {bed_path}")

        minimap2 = ["minimap2", "-ax", "splice:sr", "--secondary=yes", "-j", str(bed_path)]
        minimap2 += [str(fasta_path), str(r1_path), str(r2_path)]
        # name-sorted for the BAM parser, read from stdin
        sort = ["samtools", "sort", "-n", "-o", str(bam_path), "-"]
        logger.debug(f"Running: {' '.join(minimap2)} | {' '.join(sort)}")
        _minimap2_into_sort(minimap2, sort)
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


def _check_gdna_fraction(
    gdna_fraction: float | None, n_rna_fragments: int | None, gdna_config: GDNAConfig | None
) -> None:
    """Refuse a ``gdna_fraction`` :meth:`Scenario.build_oracle` would ignore, which would simulate no gDNA from a
    call that asked for some."""
    if gdna_fraction is None:
        return
    if n_rna_fragments is None:
        raise ValueError(
            f"gdna_fraction={gdna_fraction} is read only on top of a fixed RNA depth: pass n_rna_fragments "
            "(with a GDNAConfig) instead of n_fragments, or drop gdna_fraction and set the gDNA share of the "
            "fixed total through gdna_config's abundance"
        )
    if gdna_fraction and gdna_config is None:
        raise ValueError(
            f"gdna_fraction={gdna_fraction} needs a gdna_config to draw its gDNA from: pass "
            "gdna_config=GDNAConfig()"
        )


def _gdna_abundance(gdna_config: GDNAConfig | None) -> float:
    """The gDNA abundance the pool split weighs, 0 with no gDNA config."""
    return gdna_config.abundance if gdna_config is not None else 0.0


def _annotated(rows: list[Transcript]) -> list[Transcript]:
    return [t for t in rows if not t.is_synthetic]


def _minimap2_into_sort(minimap2: list[str], sort: list[str]) -> None:
    """Run ``minimap2 | samtools sort``. A failed sort raises; a failed minimap2 only warns."""
    p1 = subprocess.Popen(minimap2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(sort, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
            "minimap2 exited with rc=%d: %s", p1.returncode, p1_stderr.decode(errors="replace")
        )
