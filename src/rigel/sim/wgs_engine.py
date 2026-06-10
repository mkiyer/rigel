"""Whole-genome read-simulation engine (the single compute engine).

``WholeGenomeSimulator`` — vectorized (numpy), parallel (fork-based sharding), FASTA-backed
(pysam) generation of paired-end FASTQ + an optional name-sorted oracle BAM, with mRNA / nRNA /
gDNA pools, strand-specificity, hybrid-capture weighting, and gDNA strand overdispersion. Config
lives in :mod:`wgs_config`; the suite frontend (config parsing, abundance, run_simulation, CLI) is
:mod:`whole_genome`; the multi-condition loop is :mod:`orchestrator`.
"""

from __future__ import annotations

import gzip
import logging
import multiprocessing
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam

try:
    import pgzip
except ImportError:
    pgzip = None  # type: ignore[assignment]

from rigel.transcript import Transcript
from rigel.types import Strand
from rigel.sim.bam import (
    BASE_R1_FLAG,
    BASE_R2_FLAG,
    FLAG_MATE_REVERSE,
    FLAG_REVERSE,
    blocks_to_cigar,
    make_aligned_segment,
    premrna_to_genomic_interval,
    transcript_to_genomic_blocks,
)
from rigel.sim.capture import CaptureConfig, CaptureSampler
from rigel.sim.genome import reverse_complement
from rigel.sim.sampling import truncated_normal_frag_lengths
from rigel.sim.wgs_config import GDNASimConfig, SimulationParams

logger = logging.getLogger(__name__)


# Byte-level complement table for vectorized reverse-complement
_BYTE_COMPLEMENT = np.zeros(256, dtype=np.uint8)
for _c, _rc in zip(b"ACGTNacgtn", b"TGCANtgcan"):
    _BYTE_COMPLEMENT[_c] = _rc

_FASTQ_BUFFER_SIZE = 100_000


def _seq_to_bytes(seq: str) -> np.ndarray:
    """Convert a DNA string to a numpy uint8 array (writable copy)."""
    return np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()


def _batch_extract_reads(
    seq_bytes: np.ndarray,
    frag_starts: np.ndarray,
    frag_len: int,
    read_len: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized extraction of R1 and R2 read sequences.

    Given a source sequence as a byte array, extract all fragments at once
    using stride tricks, then slice R2 (first read_len bases) and R1
    (reverse complement of last read_len bases).

    Returns (r2_seqs, r1_seqs) as 2D uint8 arrays of shape (count, read_len).
    """
    offsets = np.arange(read_len, dtype=np.int64)

    # R2: first read_len bases of each fragment
    r2_indices = frag_starts[:, np.newaxis] + offsets[np.newaxis, :]
    r2_seqs = seq_bytes[r2_indices]

    # R1: reverse complement of last read_len bases
    r1_start = frag_starts + (frag_len - read_len)
    r1_indices = r1_start[:, np.newaxis] + offsets[np.newaxis, :]
    r1_seqs = _BYTE_COMPLEMENT[seq_bytes[r1_indices][:, ::-1]]

    return r2_seqs, r1_seqs


def _bam_seq_from_fastq_bytes(seq: np.ndarray, is_reverse: bool) -> str:
    """Return the SAM/BAM stored SEQ for a FASTQ-oriented read sequence."""
    if is_reverse:
        seq = _BYTE_COMPLEMENT[seq[::-1]]
    return seq.tobytes().decode("ascii")


def _open_gzip_write(path: Path, threads: int = 4):
    """Open a gzip file for text writing, using pgzip if available."""
    if pgzip is not None:
        return pgzip.open(path, "wt", thread=threads, compresslevel=3)
    return gzip.open(path, "wt", compresslevel=3)


class _FastqBuffer:
    """Buffered FASTQ writer that accumulates records before flushing.

    Reduces the number of gzip write calls from ~8 per fragment to
    ~8 per flush (every ``buf_size`` fragments).
    """

    __slots__ = ("_fh", "_parts", "_size", "_count")

    def __init__(self, fh, buf_size: int = _FASTQ_BUFFER_SIZE):
        self._fh = fh
        self._parts: list[str] = []
        self._size = buf_size
        self._count = 0

    def write_record(self, qname: str, seq: str, quals: str) -> None:
        self._parts.append(f"@{qname}\n{seq}\n+\n{quals}\n")
        self._count += 1
        if self._count >= self._size:
            self.flush()

    def write_records_batch(
        self,
        qnames: list[str],
        seqs: np.ndarray,
        quals: str,
        suffix: str,
    ) -> None:
        """Write a batch of FASTQ records from numpy byte arrays.

        Parameters
        ----------
        qnames : list of query names (without /1 or /2 suffix)
        seqs : 2D uint8 array (count x read_len)
        quals : quality string (same for all reads)
        suffix : "/1" or "/2"
        """
        parts = self._parts
        for i in range(len(qnames)):
            seq_str = seqs[i].tobytes().decode("ascii")
            parts.append(f"@{qnames[i]}{suffix}\n{seq_str}\n+\n{quals}\n")
        self._count += len(qnames)
        if self._count >= self._size:
            self.flush()

    def flush(self) -> None:
        if self._parts:
            self._fh.write("".join(self._parts))
            self._parts.clear()
            self._count = 0

    def close(self) -> None:
        self.flush()


# ═══════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════


class GenomeCache:
    """Caches full chromosome sequences as numpy byte arrays.

    Eliminates per-fragment ``pysam.FastaFile.fetch()`` overhead which
    dominates gDNA simulation time.  Chromosomes are loaded lazily on
    first access and stored as uppercase uint8 arrays.
    """

    __slots__ = ("_fasta", "_cache")

    def __init__(self, fasta: pysam.FastaFile):
        self._fasta = fasta
        self._cache: dict[str, np.ndarray] = {}

    def get(self, ref: str) -> np.ndarray:
        arr = self._cache.get(ref)
        if arr is None:
            seq = self._fasta.fetch(ref).upper()
            arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
            self._cache[ref] = arr
        return arr

    def preload(self, refs: list[str]) -> None:
        for ref in refs:
            self.get(ref)

    def clear(self) -> None:
        self._cache.clear()


# ═══════════════════════════════════════════════════════════════════
# Whole-genome read simulator
# ═══════════════════════════════════════════════════════════════════


class WholeGenomeSimulator:
    """Whole-genome paired-end RNA-seq read simulator.

    Generates gzipped FASTQ and collated oracle BAM in a single pass.

    Performance strategy
    --------------------
    1. Pre-cache all chromosome sequences as numpy byte arrays at init.
    2. Pre-extract all mRNA and pre-mRNA sequences as byte arrays.
    3. Sample fragment lengths in bulk, then sample transcript
       assignments per unique fragment length.
    4. Vectorized read extraction via numpy fancy indexing.
    5. Multi-threaded gzip compression via pgzip (if available).
    6. Unified RNA read writer handles both mRNA and nRNA.

    Architecture
    ------------
    1. Pre-compute abundance weight vectors (numpy arrays).
    2. Sample from fragment-length distribution, count fragments per
       unique length.
    3. For each fragment length, sample mRNA/nRNA transcripts
       proportional to abundance × effective_length.  Accumulate
       per-transcript fragment counts.
    4. Iterate transcripts: pull sequence ONCE, generate all
       fragments, write FASTQ.gz + BAM simultaneously.
    """

    def __init__(
        self,
        fasta_path: str | Path,
        transcripts: list[Transcript],
        sim_params: SimulationParams,
        gdna_config: GDNASimConfig,
        *,
        strand_specificity: float = 1.0,
        seed: int | None = None,
        capture_config: CaptureConfig | None = None,
    ):
        self.fasta = pysam.FastaFile(str(fasta_path))
        self.transcripts = transcripts
        self.sim_params = sim_params
        self.gdna_config = gdna_config
        self.strand_specificity = strand_specificity
        eff_seed = seed if seed is not None else sim_params.sim_seed
        self._rng = np.random.default_rng(eff_seed)
        # Nascent RNA draws on a DEDICATED, independent stream (entropy key "NRNA"=0x4E524E41) so
        # toggling nascent on/off never perturbs the mature-RNA or gDNA streams: a nascent-off run is
        # bit-identical to its nascent-on twin minus the nascent layer (head-to-head benchmarking).
        self._nrna_rng = np.random.default_rng(
            np.random.SeedSequence([int(eff_seed), 0x4E524E41])
        )

        N = len(transcripts)
        self._genome = GenomeCache(self.fasta)

        # Pre-extract spliced mRNA sequences as byte arrays
        logger.info("Pre-extracting %d mRNA sequences...", N)
        self._t_seq_bytes: list[np.ndarray] = []
        self._t_lengths = np.empty(N, dtype=np.int64)
        for i, t in enumerate(transcripts):
            seq = self._extract_mrna_seq(t)
            self._t_seq_bytes.append(_seq_to_bytes(seq))
            self._t_lengths[i] = len(seq)

        # Pre-mRNA lengths and sequences
        self._premrna_lengths = np.array(
            [t.end - t.start for t in transcripts], dtype=np.int64,
        )

        # Pre-extract pre-mRNA sequences for nRNA-eligible transcripts
        logger.info("Pre-extracting pre-mRNA sequences for nRNA...")
        self._premrna_seq_bytes: list[np.ndarray | None] = [None] * N
        for i, t in enumerate(transcripts):
            if t.nrna_abundance > 0 and len(t.exons) > 1:
                self._premrna_seq_bytes[i] = self._fetch_premrna_bytes(t)

        # Strand chars and abundance arrays
        self._t_strand_chars: list[str] = [
            "r" if t.strand == Strand.NEG else "f" for t in transcripts
        ]
        self._mrna_abund = np.array([t.abundance or 0.0 for t in transcripts])
        self._nrna_abund = np.array([t.nrna_abundance for t in transcripts])

        # Quality string cache
        self._quals_cache: dict[int, str] = {}

        # Fast-path flag for error introduction
        self._has_errors = sim_params.error_rate > 0

        # Clamp frag_max to longest transcript
        max_mrna = int(self._t_lengths.max()) if N > 0 else 0
        self._frag_max = min(sim_params.frag_max, max_mrna) if max_mrna > 0 else sim_params.frag_max

        # gDNA setup: annotated references weighted by length
        annotated_refs = {t.ref for t in transcripts}
        self._gdna_refs: list[str] = []
        self._gdna_ref_lengths: list[int] = []
        for ref in self.fasta.references:
            if ref in annotated_refs:
                self._gdna_refs.append(ref)
                self._gdna_ref_lengths.append(self.fasta.get_reference_length(ref))

        # Pre-load gDNA chromosomes into cache (major perf win)
        if self._gdna_refs:
            logger.info("Pre-loading %d chromosomes for gDNA...", len(self._gdna_refs))
            self._genome.preload(self._gdna_refs)

        # gDNA strand overdispersion: per-ref exon-derived regions each with a shared +strand
        # rate ~ Beta(a, a), so a region's gDNA strand split is Beta-Binomial(½, overdispersion)
        # — matching the calibrator's seed-region partition. 0 ⇒ uniform 50/50 (no regions built).
        self._gdna_strand_regions: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        gdna_od = float(getattr(self.gdna_config, "strand_overdispersion", 0.0) or 0.0)
        if gdna_od > 0.0 and self._gdna_refs:
            self._init_gdna_strand_regions(gdna_od)

        # BAM header and reference mapping (built once, reused)
        self._ref_names = list(self.fasta.references)
        self._ref_lengths = [self.fasta.get_reference_length(r) for r in self._ref_names]
        self._ref_name_to_id: dict[str, int] = {
            name: i for i, name in enumerate(self._ref_names)
        }
        self.capture = CaptureSampler.from_config(
            capture_config,
            transcripts,
            dict(zip(self._ref_names, self._ref_lengths)),
        )
        self._bam_header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "unsorted", "GO": "query"},
            "SQ": [
                {"SN": name, "LN": length}
                for name, length in zip(self._ref_names, self._ref_lengths)
            ],
            "PG": [{
                "ID": "rigel_sim",
                "PN": "rigel_sim",
                "VN": "1.0",
                "CL": "simulated",
            }],
        })

        # Pre-built shared tag lists for BAM records
        self._nh1_tags = [("NH", 1)]

        logger.info(
            "Simulator ready: %d transcripts, %d gDNA refs, max mRNA len=%d",
            N, len(self._gdna_refs), max_mrna,
        )

    # -- Sequence extraction ------------------------------------------------

    def _extract_mrna_seq(self, t: Transcript) -> str:
        """Extract spliced mRNA sequence (5'-to-3' oriented)."""
        exon_seqs = []
        for e in t.exons:
            exon_seqs.append(self.fasta.fetch(t.ref, e.start, e.end).upper())
        seq = "".join(exon_seqs)
        if t.strand == Strand.NEG:
            seq = reverse_complement(seq)
        return seq

    def _fetch_premrna_bytes(self, t: Transcript) -> np.ndarray:
        """Fetch unspliced pre-mRNA sequence as uint8 array."""
        seq = self.fasta.fetch(t.ref, t.start, t.end).upper()
        if t.strand == Strand.NEG:
            seq = reverse_complement(seq)
        return _seq_to_bytes(seq)

    def _get_premrna_bytes(self, t_idx: int) -> np.ndarray:
        """Get pre-mRNA bytes, fetching on demand if not pre-cached."""
        arr = self._premrna_seq_bytes[t_idx]
        if arr is None:
            arr = self._fetch_premrna_bytes(self.transcripts[t_idx])
            self._premrna_seq_bytes[t_idx] = arr
        return arr

    def _get_quals(self, read_len: int) -> str:
        """Return cached quality string of the given length."""
        q = self._quals_cache.get(read_len)
        if q is None:
            q = "I" * read_len
            self._quals_cache[read_len] = q
        return q

    # -- Fragment length sampling -------------------------------------------

    def _sample_rna_frag_lengths(self, n: int) -> np.ndarray:
        p = self.sim_params
        return truncated_normal_frag_lengths(
            self._rng, n, p.frag_mean, p.frag_std, p.frag_min, self._frag_max
        )

    def _init_gdna_strand_regions(self, overdispersion: float) -> None:
        """Build per-gDNA-ref exon-derived regions, each with a shared +strand rate ~ Beta(a, a).

        ``Beta(a, a)`` has intra-class correlation ``1/(2a + 1)``, so ``a = ½(1 − od)/od`` gives a
        region-to-region strand overdispersion of ``od``. All gDNA fragments in a region share its
        rate, so the per-region sense/antisense count is Beta-Binomial(½, od) — the distribution
        the calibrator's gDNA strand model fits.
        """
        alpha = 0.5 * (1.0 - overdispersion) / overdispersion
        exons_by_ref: dict[str, set[int]] = defaultdict(set)
        for t in self.transcripts:
            if t.ref is None:
                continue
            for e in t.exons:
                exons_by_ref[t.ref].add(int(e.start))
                exons_by_ref[t.ref].add(int(e.end))
        for ref, length in zip(self._gdna_refs, self._gdna_ref_lengths):
            bounds = {b for b in exons_by_ref.get(ref, set()) if 0 < b < length}
            boundaries = np.array([0, *sorted(bounds), int(length)], dtype=np.int64)
            n_regions = max(len(boundaries) - 1, 1)
            p_plus = self._rng.beta(alpha, alpha, size=n_regions)
            self._gdna_strand_regions[ref] = (boundaries, p_plus)

    def _sample_gdna_frag_lengths(self, n: int) -> np.ndarray:
        g = self.gdna_config
        return truncated_normal_frag_lengths(
            self._rng, n, g.frag_mean, g.frag_std, g.frag_min, g.frag_max
        )

    # -- Error introduction -------------------------------------------------

    def _introduce_errors_batch(self, seqs: np.ndarray) -> np.ndarray:
        """Apply random substitution errors to a batch of sequences in-place.

        Parameters
        ----------
        seqs : 2D uint8 array (count x read_len)
        """
        if not self._has_errors:
            return seqs
        count, read_len = seqs.shape
        mask = self._rng.random((count, read_len)) < self.sim_params.error_rate
        n_errors = int(mask.sum())
        if n_errors > 0:
            bases = np.array([65, 67, 71, 84], dtype=np.uint8)  # A, C, G, T
            seqs[mask] = bases[self._rng.integers(4, size=n_errors)]
        return seqs

    # -- Fragment count accumulation ----------------------------------------

    def _accumulate_pool(
        self,
        n_frags: int,
        abundances: np.ndarray,
        lengths: np.ndarray,
        *,
        space: str,
    ) -> dict[int, dict[int, int]]:
        """Sample fragment counts from a single abundance/length pool.

        Returns dict[t_idx, dict[frag_len, count]].
        """
        if n_frags <= 0:
            return {}

        rng = self._rng
        frag_lengths = self._sample_rna_frag_lengths(n_frags)
        unique_lengths, length_counts = np.unique(frag_lengths, return_counts=True)

        counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))

        for fl, fc in zip(unique_lengths, length_counts):
            fl, fc = int(fl), int(fc)
            eff = self.capture.partition_array(space, range(len(abundances)), lengths, fl)
            weights = abundances * eff
            total_w = weights.sum()
            if total_w <= 0:
                continue
            probs = weights / total_w
            indices = rng.choice(len(abundances), size=fc, p=probs)
            unique_idx, idx_counts = np.unique(indices, return_counts=True)
            for idx, cnt in zip(unique_idx, idx_counts):
                counts[int(idx)][fl] += int(cnt)

        return dict(counts)

    def _accumulate_rna_counts(
        self,
        n_rna: int,
        *,
        n_mrna: int | None = None,
        n_nrna: int | None = None,
    ) -> tuple[dict[int, dict[int, int]], dict[int, dict[int, int]]]:
        """Sample fragment lengths and transcript assignments.

        If *n_mrna* and *n_nrna* are given, samples each pool
        independently (fragment-count control).  Otherwise falls back
        to the combined-pool approach where *n_rna* fragments are
        drawn jointly from mRNA + nRNA weighted by abundance × length.

        Returns (mrna_counts, nrna_counts) where each is
        dict[t_idx, dict[frag_len, count]].
        """
        # Separate-pool mode: precise fragment-level control
        if n_mrna is not None and n_nrna is not None:
            mrna_counts = self._accumulate_pool(
                n_mrna, self._mrna_abund, self._t_lengths, space="mrna",
            )
            # Nascent on the dedicated stream; restore self._rng so the gDNA pool (next) and the
            # mature pool (above) are untouched by whether/how much nascent was drawn.
            _main_rng = self._rng
            self._rng = self._nrna_rng
            nrna_counts = self._accumulate_pool(
                n_nrna, self._nrna_abund, self._premrna_lengths, space="nrna",
            )
            self._rng = _main_rng
            return mrna_counts, nrna_counts

        # Combined-pool mode (original behaviour)
        if n_rna <= 0:
            return {}, {}

        N = len(self.transcripts)
        rng = self._rng

        frag_lengths = self._sample_rna_frag_lengths(n_rna)
        unique_lengths, length_counts = np.unique(frag_lengths, return_counts=True)

        mrna_counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))
        nrna_counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))

        for fl, fc in zip(unique_lengths, length_counts):
            fl, fc = int(fl), int(fc)

            mrna_eff = self.capture.partition_array("mrna", range(N), self._t_lengths, fl)
            nrna_eff = self.capture.partition_array("nrna", range(N), self._premrna_lengths, fl)
            weights = np.concatenate([
                self._mrna_abund * mrna_eff,
                self._nrna_abund * nrna_eff,
            ])
            total_w = weights.sum()
            if total_w <= 0:
                continue

            probs = weights / total_w
            indices = rng.choice(2 * N, size=fc, p=probs)
            unique_idx, idx_counts = np.unique(indices, return_counts=True)

            for idx, cnt in zip(unique_idx, idx_counts):
                idx, cnt = int(idx), int(cnt)
                if idx < N:
                    mrna_counts[idx][fl] += cnt
                else:
                    nrna_counts[idx - N][fl] += cnt

        return dict(mrna_counts), dict(nrna_counts)

    # -- Unified RNA read writer --------------------------------------------

    def _write_rna_reads(
        self,
        t_idx: int,
        len_counts: dict[int, int],
        r1_buf: _FastqBuffer,
        r2_buf: _FastqBuffer,
        bam_fh: pysam.AlignmentFile | None,
        *,
        is_nrna: bool,
    ) -> int:
        """Generate and write all mRNA or nRNA reads for one transcript.

        When ``is_nrna=False``, uses the pre-extracted spliced mRNA sequence.
        When ``is_nrna=True``, uses the pre-mRNA (unspliced) sequence.
        """
        t = self.transcripts[t_idx]
        rng = self._rng
        ss = self.strand_specificity
        strand_char = self._t_strand_chars[t_idx]
        t_id = t.t_id

        if is_nrna:
            seq_bytes = self._get_premrna_bytes(t_idx)
            seq_len = int(self._premrna_lengths[t_idx])
            qname_prefix = f"nrna_{t_id}"
        else:
            seq_bytes = self._t_seq_bytes[t_idx]
            seq_len = int(self._t_lengths[t_idx])
            qname_prefix = t_id

        ref_id = self._ref_name_to_id.get(t.ref) if bam_fh else None
        n_written = 0

        for frag_len in sorted(len_counts):
            count = len_counts[frag_len]
            read_len = min(self.sim_params.read_length, frag_len)
            eff_len = seq_len - frag_len + 1
            if eff_len <= 0:
                continue

            capture_space = "nrna" if is_nrna else "mrna"
            frag_starts = self.capture.sample_starts(
                capture_space, t_idx, seq_len, frag_len, count, rng,
            )
            flip_mask = rng.random(count) >= ss if ss < 1.0 else None
            quals = self._get_quals(read_len)

            # Vectorized read extraction
            r2_seqs, r1_seqs = _batch_extract_reads(seq_bytes, frag_starts, frag_len, read_len)

            # Strand flipping
            if flip_mask is not None:
                flip_idx = np.where(flip_mask)[0]
                if len(flip_idx) > 0:
                    r1_seqs[flip_idx], r2_seqs[flip_idx] = (
                        r2_seqs[flip_idx].copy(),
                        r1_seqs[flip_idx].copy(),
                    )

            # Batch error introduction
            if self._has_errors:
                self._introduce_errors_batch(r1_seqs)
                self._introduce_errors_batch(r2_seqs)

            # Build query names
            base_n = n_written
            qnames = [
                f"{qname_prefix}:{int(frag_starts[i])}-{int(frag_starts[i]) + frag_len}:{strand_char}:{base_n + i}"
                for i in range(count)
            ]

            # Batch FASTQ writes
            r1_buf.write_records_batch(qnames, r1_seqs, quals, "/1")
            r2_buf.write_records_batch(qnames, r2_seqs, quals, "/2")

            # BAM records (pysam API requires per-record construction)
            if bam_fh is not None and ref_id is not None:
                if is_nrna:
                    self._write_nrna_bam_batch(
                        bam_fh, qnames, r1_seqs, r2_seqs, t,
                        frag_starts, frag_len, read_len, flip_mask, ref_id,
                    )
                else:
                    self._write_mrna_bam_batch(
                        bam_fh, qnames, r1_seqs, r2_seqs, t,
                        frag_starts, frag_len, read_len, flip_mask, ref_id,
                    )

            n_written += count

        return n_written

    # -- BAM batch writers --------------------------------------------------

    def _write_mrna_bam_batch(
        self,
        bam_fh: pysam.AlignmentFile,
        qnames: list[str],
        r1_seqs: np.ndarray,
        r2_seqs: np.ndarray,
        t: Transcript,
        frag_starts: np.ndarray,
        frag_len: int,
        read_len: int,
        flip_mask: np.ndarray | None,
        ref_id: int,
    ) -> None:
        for i in range(len(qnames)):
            frag_start = int(frag_starts[i])
            frag_end = frag_start + frag_len
            flipped = flip_mask is not None and flip_mask[i]

            r2_t_end = min(frag_start + read_len, frag_end)
            r2_blocks = transcript_to_genomic_blocks(frag_start, r2_t_end, t)
            r1_t_start = max(frag_end - read_len, frag_start)
            r1_blocks = transcript_to_genomic_blocks(r1_t_start, frag_end, t)
            if not r2_blocks or not r1_blocks:
                continue

            if t.strand == Strand.POS:
                r2_is_rev, r1_is_rev = False, True
            else:
                r2_is_rev, r1_is_rev = True, False
            if flipped:
                r2_is_rev, r1_is_rev = not r2_is_rev, not r1_is_rev
                r1_blocks, r2_blocks = r2_blocks, r1_blocks

            r2_cigar = blocks_to_cigar(r2_blocks)
            r1_cigar = blocks_to_cigar(r1_blocks)
            r2_start = r2_blocks[0][0]
            r1_start = r1_blocks[0][0]

            leftmost = min(r1_blocks[0][0], r2_blocks[0][0])
            rightmost = max(r1_blocks[-1][1], r2_blocks[-1][1])
            tlen = rightmost - leftmost

            tags = self._nh1_tags
            if len(r2_blocks) > 1 or len(r1_blocks) > 1:
                xs_strand = "-" if t.strand == Strand.NEG else "+"
                tags = [("NH", 1), ("XS", xs_strand, "A")]  # type 'A' (char), per SAM convention

            r1_flag = BASE_R1_FLAG
            r2_flag = BASE_R2_FLAG
            if r1_is_rev:
                r1_flag |= FLAG_REVERSE
            if r2_is_rev:
                r1_flag |= FLAG_MATE_REVERSE
                r2_flag |= FLAG_REVERSE
            if r1_is_rev:
                r2_flag |= FLAG_MATE_REVERSE

            r1_tlen = tlen if r1_start <= r2_start else -tlen
            r2_tlen = -r1_tlen

            r1_seq_str = _bam_seq_from_fastq_bytes(r1_seqs[i], r1_is_rev)
            r2_seq_str = _bam_seq_from_fastq_bytes(r2_seqs[i], r2_is_rev)

            bam_fh.write(make_aligned_segment(
                self._bam_header, qnames[i], r1_seq_str, r1_flag, ref_id,
                r1_start, r1_cigar, ref_id, r2_start, r1_tlen, tags=tags,
            ))
            bam_fh.write(make_aligned_segment(
                self._bam_header, qnames[i], r2_seq_str, r2_flag, ref_id,
                r2_start, r2_cigar, ref_id, r1_start, r2_tlen, tags=tags,
            ))

    def _write_nrna_bam_batch(
        self,
        bam_fh: pysam.AlignmentFile,
        qnames: list[str],
        r1_seqs: np.ndarray,
        r2_seqs: np.ndarray,
        t: Transcript,
        frag_starts: np.ndarray,
        frag_len: int,
        read_len: int,
        flip_mask: np.ndarray | None,
        ref_id: int,
    ) -> None:
        for i in range(len(qnames)):
            frag_start = int(frag_starts[i])
            frag_end = frag_start + frag_len
            flipped = flip_mask is not None and flip_mask[i]

            g_start, g_end = premrna_to_genomic_interval(frag_start, frag_end, t)
            if t.strand == Strand.POS:
                r2_g_start = g_start
                r2_g_end = min(g_start + read_len, g_end)
                r1_g_start = max(g_end - read_len, g_start)
                r1_g_end = g_end
                r2_is_rev, r1_is_rev = False, True
            else:
                r1_g_start = g_start
                r1_g_end = min(g_start + read_len, g_end)
                r2_g_start = max(g_end - read_len, g_start)
                r2_g_end = g_end
                r2_is_rev, r1_is_rev = True, False
            if flipped:
                r2_is_rev, r1_is_rev = not r2_is_rev, not r1_is_rev
                r1_g_start, r2_g_start = r2_g_start, r1_g_start
                r1_g_end, r2_g_end = r2_g_end, r1_g_end

            tlen = g_end - g_start

            r1_flag = BASE_R1_FLAG
            r2_flag = BASE_R2_FLAG
            if r1_is_rev:
                r1_flag |= FLAG_REVERSE
            if r2_is_rev:
                r1_flag |= FLAG_MATE_REVERSE
                r2_flag |= FLAG_REVERSE
            if r1_is_rev:
                r2_flag |= FLAG_MATE_REVERSE

            r1_tlen = tlen if r1_g_start <= r2_g_start else -tlen
            r2_tlen = -r1_tlen

            r1_read_len = r1_g_end - r1_g_start
            r2_read_len = r2_g_end - r2_g_start

            r1_seq_str = _bam_seq_from_fastq_bytes(r1_seqs[i], r1_is_rev)
            r2_seq_str = _bam_seq_from_fastq_bytes(r2_seqs[i], r2_is_rev)

            bam_fh.write(make_aligned_segment(
                self._bam_header, qnames[i], r1_seq_str, r1_flag, ref_id,
                r1_g_start, [(pysam.CMATCH, r1_read_len)],
                ref_id, r2_g_start, r1_tlen, tags=self._nh1_tags,
            ))
            bam_fh.write(make_aligned_segment(
                self._bam_header, qnames[i], r2_seq_str, r2_flag, ref_id,
                r2_g_start, [(pysam.CMATCH, r2_read_len)],
                ref_id, r1_g_start, r2_tlen, tags=self._nh1_tags,
            ))

    # -- gDNA reads (vectorized per chromosome) -----------------------------

    def _accumulate_gdna_counts(self, n_gdna: int) -> dict[tuple[int, int], int]:
        """Sample gDNA fragment lengths + chromosomes.

        Returns dict[(ref_idx, frag_len)] = count.
        """
        if n_gdna == 0 or not self._gdna_refs:
            return {}

        rng = self._rng
        frag_lengths = self._sample_gdna_frag_lengths(n_gdna)
        unique_lengths, length_counts = np.unique(frag_lengths, return_counts=True)
        counts: dict[tuple[int, int], int] = {}
        for fl, fc in zip(unique_lengths, length_counts):
            fl, fc = int(fl), int(fc)
            chrom_eff = np.array(
                [
                    self.capture.partition("gdna", ref, ref_len, fl)
                    for ref, ref_len in zip(self._gdna_refs, self._gdna_ref_lengths)
                ],
                dtype=np.float64,
            )
            total_eff = chrom_eff.sum()
            if total_eff <= 0:
                continue
            chrom_probs = chrom_eff / total_eff
            chrom_indices = rng.choice(len(self._gdna_refs), size=fc, p=chrom_probs)
            unique_chroms, chrom_counts = np.unique(chrom_indices, return_counts=True)
            for ci, cc in zip(unique_chroms, chrom_counts):
                counts[(int(ci), fl)] = int(cc)
        return counts

    def _write_gdna_chunk(
        self,
        ref_idx: int,
        fl: int,
        count: int,
        r1_buf: _FastqBuffer,
        r2_buf: _FastqBuffer,
        bam_fh: pysam.AlignmentFile | None,
        n_offset: int = 0,
    ) -> int:
        """Generate ``count`` gDNA fragments at chromosome ref_idx with frag length fl.

        Returns number of fragments written.
        """
        if count <= 0:
            return 0
        rng = self._rng
        ref = self._gdna_refs[ref_idx]
        chrom_len = self._gdna_ref_lengths[ref_idx]
        ref_id = self._ref_name_to_id.get(ref)
        chrom_bytes = self._genome.get(ref)
        eff_len = chrom_len - fl + 1
        if eff_len <= 0:
            return 0

        read_len = min(self.sim_params.read_length, fl)
        quals = self._get_quals(read_len)

        starts = self.capture.sample_starts("gdna", ref, chrom_len, fl, count, rng)
        # Strand: uniform 50/50, or — with overdispersion — by the fragment's exon-derived region's
        # shared Beta(a, a) +strand rate (0 = +/forward, 1 = −/reverse, matching the qname encoding).
        strand_regions = self._gdna_strand_regions.get(ref)
        if strand_regions is None:
            chrom_strands = rng.integers(0, 2, size=count)
        else:
            boundaries, p_plus = strand_regions
            reg_idx = np.clip(
                np.searchsorted(boundaries, starts, side="right") - 1, 0, len(p_plus) - 1
            )
            chrom_strands = (rng.random(count) >= p_plus[reg_idx]).astype(int)

        r2_seqs, r1_seqs = _batch_extract_reads(chrom_bytes, starts, fl, read_len)

        # Negative strand: swap R1/R2 (fragment = revcomp of genomic)
        neg_mask = chrom_strands.astype(bool)
        if neg_mask.any():
            neg_idx = np.where(neg_mask)[0]
            r1_seqs[neg_idx], r2_seqs[neg_idx] = (
                r2_seqs[neg_idx].copy(),
                r1_seqs[neg_idx].copy(),
            )

        if self._has_errors:
            self._introduce_errors_batch(r1_seqs)
            self._introduce_errors_batch(r2_seqs)

        qnames = [
            f"gdna:{ref}:{int(starts[j])}-{int(starts[j]) + fl}:"
            f"{'r' if chrom_strands[j] else 'f'}:{n_offset + j}"
            for j in range(count)
        ]

        r1_buf.write_records_batch(qnames, r1_seqs, quals, "/1")
        r2_buf.write_records_batch(qnames, r2_seqs, quals, "/2")

        if bam_fh is not None and ref_id is not None:
            for j in range(count):
                start = int(starts[j])
                end = start + fl
                is_neg = bool(chrom_strands[j])

                if is_neg:
                    r1_is_rev = False
                    r2_is_rev = True
                    r1_start_pos = start
                    r2_start_pos = end - read_len
                else:
                    r1_is_rev = True
                    r2_is_rev = False
                    r1_start_pos = end - read_len
                    r2_start_pos = start

                r1_flag = BASE_R1_FLAG
                r2_flag = BASE_R2_FLAG
                if r1_is_rev:
                    r1_flag |= FLAG_REVERSE
                if r2_is_rev:
                    r1_flag |= FLAG_MATE_REVERSE
                    r2_flag |= FLAG_REVERSE
                if r1_is_rev:
                    r2_flag |= FLAG_MATE_REVERSE

                tlen = end - start
                r1_tlen = tlen if r1_start_pos <= r2_start_pos else -tlen
                r2_tlen = -r1_tlen

                r1_seq_str = _bam_seq_from_fastq_bytes(r1_seqs[j], r1_is_rev)
                r2_seq_str = _bam_seq_from_fastq_bytes(r2_seqs[j], r2_is_rev)

                bam_fh.write(make_aligned_segment(
                    self._bam_header, qnames[j], r1_seq_str,
                    r1_flag, ref_id, r1_start_pos,
                    [(pysam.CMATCH, read_len)],
                    ref_id, r2_start_pos, r1_tlen,
                    tags=self._nh1_tags,
                ))
                bam_fh.write(make_aligned_segment(
                    self._bam_header, qnames[j], r2_seq_str,
                    r2_flag, ref_id, r2_start_pos,
                    [(pysam.CMATCH, read_len)],
                    ref_id, r1_start_pos, r2_tlen,
                    tags=self._nh1_tags,
                ))
        return count

    def _write_gdna_from_counts(
        self,
        gdna_counts: dict[tuple[int, int], int],
        r1_buf: _FastqBuffer,
        r2_buf: _FastqBuffer,
        bam_fh: pysam.AlignmentFile | None,
    ) -> int:
        """Write gDNA reads from a pre-sampled (ref_idx, fl) -> count map."""
        n_written = 0
        for (ref_idx, fl), count in gdna_counts.items():
            n_written += self._write_gdna_chunk(
                ref_idx, fl, count, r1_buf, r2_buf, bam_fh, n_offset=n_written,
            )
        return n_written

    # -- Abundance-weighted pool splits (single-condition / Scenario use) ----

    def _rna_eff_weights(self) -> tuple[float, float]:
        """``(mrna_weight, nrna_weight)`` = abundance × capture-effective-length at the mean RNA
        fragment length — the basis for splitting an RNA pool into mature vs nascent."""
        mean_frag = int(self.sim_params.frag_mean)
        n = len(self.transcripts)
        mrna_eff = self.capture.partition_array("mrna", range(n), self._t_lengths, mean_frag)
        nrna_eff = self.capture.partition_array("nrna", range(n), self._premrna_lengths, mean_frag)
        return (
            float(np.sum(self._mrna_abund * mrna_eff)),
            float(np.sum(self._nrna_abund * nrna_eff)),
        )

    def rna_split(self, n_rna: int) -> tuple[int, int]:
        """Split ``n_rna`` RNA fragments into ``(n_mrna, n_nrna)`` by abundance × effective length."""
        mrna_w, nrna_w = self._rna_eff_weights()
        total = mrna_w + nrna_w
        if total <= 0:
            return 0, 0
        n_nrna = int(round(n_rna * nrna_w / total))
        return max(0, n_rna - n_nrna), n_nrna

    def pool_split(self, n_total: int, gdna_abundance: float) -> tuple[int, int, int]:
        """Split ``n_total`` fragments 3-way into ``(n_mrna, n_nrna, n_gdna)`` by abundance ×
        effective length, with gDNA weighted by ``gdna_abundance × genome effective length``
        (summed over the annotated references). The abundance-weighted equivalent of the old
        ``ReadSimulator._compute_pool_split``."""
        mrna_w, nrna_w = self._rna_eff_weights()
        gdna_mean_frag = int(self.gdna_config.frag_mean)
        genome_eff = sum(
            self.capture.partition("gdna", ref, length, gdna_mean_frag)
            for ref, length in zip(self._gdna_refs, self._gdna_ref_lengths)
        )
        gdna_w = float(gdna_abundance) * genome_eff
        total = mrna_w + nrna_w + gdna_w
        if total <= 0:
            return 0, 0, 0
        n_nrna = int(round(n_total * nrna_w / total))
        n_gdna = int(round(n_total * gdna_w / total))
        return max(0, n_total - n_nrna - n_gdna), n_nrna, n_gdna

    # -- Main entry point ---------------------------------------------------

    def simulate_and_write(
        self,
        output_dir: Path,
        n_rna: int,
        n_gdna: int = 0,
        *,
        n_mrna: int | None = None,
        n_nrna: int | None = None,
        oracle_bam: bool = True,
        prefix: str = "sim",
        n_workers: int = 1,
    ) -> tuple[Path, Path, Path | None]:
        """Single-pass simulation: accumulate counts, generate, write.

        When *n_mrna* and *n_nrna* are provided, the mRNA and nRNA
        fragment pools are sampled independently (precise fragment-count
        control).  *n_rna* is ignored in this mode but still used for
        logging.  Otherwise *n_rna* fragments are drawn from the
        combined mRNA + nRNA pool (original behaviour).

        When *n_workers* > 1 the per-transcript and per-(chrom, frag-len)
        gDNA work is sharded across worker processes (fork-based) and the
        resulting FASTQ.gz / BAM shard files are concatenated in the
        parent.

        Returns (r1_path, r2_path, bam_path | None).
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        r1_path = output_dir / f"{prefix}_R1.fq.gz"
        r2_path = output_dir / f"{prefix}_R2.fq.gz"
        bam_path = output_dir / f"{prefix}_oracle.bam" if oracle_bam else None

        # Phase 1: accumulate per-transcript / per-chunk counts (main process).
        # Uses self._rng so determinism follows the configured seed.
        logger.info("Accumulating RNA fragment counts...")
        mrna_counts, nrna_counts = self._accumulate_rna_counts(
            n_rna, n_mrna=n_mrna, n_nrna=n_nrna,
        )
        gdna_counts = self._accumulate_gdna_counts(n_gdna)

        total_mrna = sum(sum(d.values()) for d in mrna_counts.values())
        total_nrna = sum(sum(d.values()) for d in nrna_counts.values())
        total_gdna = sum(gdna_counts.values())
        logger.info(
            "Fragment counts: %d mRNA (%d txs), %d nRNA (%d txs), %d gDNA chunks",
            total_mrna, len(mrna_counts), total_nrna, len(nrna_counts), len(gdna_counts),
        )

        if n_workers <= 1:
            self._write_all_serial(
                r1_path, r2_path, bam_path, mrna_counts, nrna_counts, gdna_counts,
                oracle_bam,
            )
        else:
            self._write_all_parallel(
                output_dir, prefix, r1_path, r2_path, bam_path,
                mrna_counts, nrna_counts, gdna_counts,
                oracle_bam, n_workers,
            )

        n_written = total_mrna + total_nrna + total_gdna
        logger.info(
            "Wrote %d read pairs -> %s (RNA=%d, gDNA=%d, oracle=%s, workers=%d)",
            n_written, output_dir, n_rna, n_gdna, oracle_bam, n_workers,
        )
        return r1_path, r2_path, bam_path

    # -- Serial / parallel write paths --------------------------------------

    def _write_all_serial(
        self,
        r1_path: Path,
        r2_path: Path,
        bam_path: Path | None,
        mrna_counts: dict[int, dict[int, int]],
        nrna_counts: dict[int, dict[int, int]],
        gdna_counts: dict[tuple[int, int], int],
        oracle_bam: bool,
    ) -> None:
        bam_fh: pysam.AlignmentFile | None = None
        with (
            _open_gzip_write(r1_path) as r1_fh,
            _open_gzip_write(r2_path) as r2_fh,
        ):
            r1_buf = _FastqBuffer(r1_fh)
            r2_buf = _FastqBuffer(r2_fh)
            if oracle_bam:
                bam_fh = pysam.AlignmentFile(
                    str(bam_path), "wb", header=self._bam_header,
                )
            try:
                for t_idx in sorted(mrna_counts):
                    self._write_rna_reads(
                        t_idx, mrna_counts[t_idx], r1_buf, r2_buf, bam_fh,
                        is_nrna=False,
                    )
                # Nascent reads on the dedicated stream; restore so gDNA writing (next) is unaffected.
                _main_rng = self._rng
                self._rng = self._nrna_rng
                for t_idx in sorted(nrna_counts):
                    self._write_rna_reads(
                        t_idx, nrna_counts[t_idx], r1_buf, r2_buf, bam_fh,
                        is_nrna=True,
                    )
                self._rng = _main_rng
                self._write_gdna_from_counts(gdna_counts, r1_buf, r2_buf, bam_fh)
                r1_buf.close()
                r2_buf.close()
            finally:
                if bam_fh is not None:
                    bam_fh.close()

    def _write_all_parallel(
        self,
        output_dir: Path,
        prefix: str,
        r1_path: Path,
        r2_path: Path,
        bam_path: Path | None,
        mrna_counts: dict[int, dict[int, int]],
        nrna_counts: dict[int, dict[int, int]],
        gdna_counts: dict[tuple[int, int], int],
        oracle_bam: bool,
        n_workers: int,
    ) -> None:
        # Build per-shard task lists via greedy LPT bin-packing on fragment counts.
        mrna_items = [(t, dict(d)) for t, d in mrna_counts.items()]
        nrna_items = [(t, dict(d)) for t, d in nrna_counts.items()]
        gdna_items = list(gdna_counts.items())  # [((ref_idx, fl), count), ...]

        mrna_shards = _shard_by_count(
            mrna_items, n_workers, weight=lambda x: sum(x[1].values()),
        )
        nrna_shards = _shard_by_count(
            nrna_items, n_workers, weight=lambda x: sum(x[1].values()),
        )
        gdna_shards = _shard_by_count(
            gdna_items, n_workers, weight=lambda x: x[1],
        )

        shard_dir = output_dir / f".{prefix}_shards"
        shard_dir.mkdir(parents=True, exist_ok=True)

        base_seed = int(self._rng.integers(0, 2**31 - 1))

        tasks = []
        for k in range(n_workers):
            tasks.append((
                k,
                mrna_shards[k],
                nrna_shards[k],
                gdna_shards[k],
                str(shard_dir),
                prefix,
                oracle_bam,
                base_seed + k,
            ))

        # Bind self into a module global so fork'd children inherit it without
        # pickling.  Passing self via Pool.map would fail (pysam handles).
        global _WORKER_SIM
        _WORKER_SIM = self

        ctx = multiprocessing.get_context("fork")
        results: list[tuple[str, str, str | None]]
        try:
            with ctx.Pool(n_workers) as pool:
                results = pool.map(_shard_task, tasks)
        finally:
            _WORKER_SIM = None

        shard_r1 = [Path(r[0]) for r in results]
        shard_r2 = [Path(r[1]) for r in results]
        shard_bams = [Path(r[2]) for r in results if r[2] is not None]

        logger.info("Concatenating %d FASTQ shards...", len(shard_r1))
        _concat_files_binary(shard_r1, r1_path)
        _concat_files_binary(shard_r2, r2_path)

        if oracle_bam and shard_bams:
            logger.info("Concatenating %d BAM shards...", len(shard_bams))
            pysam.cat("-o", str(bam_path), *[str(p) for p in shard_bams])

        # Cleanup shard files
        for p in shard_r1 + shard_r2 + shard_bams:
            try:
                p.unlink()
            except OSError:
                pass
        try:
            shard_dir.rmdir()
        except OSError:
            pass

    def _run_shard(
        self,
        shard_id: int,
        mrna_items: list[tuple[int, dict[int, int]]],
        nrna_items: list[tuple[int, dict[int, int]]],
        gdna_items: list[tuple[tuple[int, int], int]],
        shard_dir: str,
        prefix: str,
        oracle_bam: bool,
        seed: int,
    ) -> tuple[str, str, str | None]:
        """Worker entry: write one shard's FASTQ + BAM files."""
        # Independent RNG per shard; nascent on its own dedicated stream (head-to-head: the mature +
        # gDNA shard reads are unaffected by whether the nascent shard is written).
        self._rng = np.random.default_rng(seed)
        self._nrna_rng = np.random.default_rng(np.random.SeedSequence([int(seed), 0x4E524E41]))

        shard_path = Path(shard_dir)
        r1_path = shard_path / f"{prefix}.shard{shard_id:03d}.R1.fq.gz"
        r2_path = shard_path / f"{prefix}.shard{shard_id:03d}.R2.fq.gz"
        bam_path = shard_path / f"{prefix}.shard{shard_id:03d}.bam" if oracle_bam else None

        bam_fh: pysam.AlignmentFile | None = None
        with (
            _open_gzip_write(r1_path, threads=1) as r1_fh,
            _open_gzip_write(r2_path, threads=1) as r2_fh,
        ):
            r1_buf = _FastqBuffer(r1_fh)
            r2_buf = _FastqBuffer(r2_fh)
            if oracle_bam:
                bam_fh = pysam.AlignmentFile(
                    str(bam_path), "wb", header=self._bam_header,
                )
            try:
                for t_idx, len_counts in mrna_items:
                    self._write_rna_reads(
                        t_idx, len_counts, r1_buf, r2_buf, bam_fh, is_nrna=False,
                    )
                _main_rng = self._rng
                self._rng = self._nrna_rng  # nascent shard on the dedicated stream
                for t_idx, len_counts in nrna_items:
                    self._write_rna_reads(
                        t_idx, len_counts, r1_buf, r2_buf, bam_fh, is_nrna=True,
                    )
                self._rng = _main_rng
                n_offset = 0
                for (ref_idx, fl), count in gdna_items:
                    n_offset += self._write_gdna_chunk(
                        ref_idx, fl, count, r1_buf, r2_buf, bam_fh, n_offset=n_offset,
                    )
                r1_buf.close()
                r2_buf.close()
            finally:
                if bam_fh is not None:
                    bam_fh.close()

        return (str(r1_path), str(r2_path), str(bam_path) if bam_path else None)

    def close(self) -> None:
        """Release resources."""
        self._genome.clear()
        self.fasta.close()


# ═══════════════════════════════════════════════════════════════════
# Parallel-shard helpers (module-level for fork-based multiprocessing)
# ═══════════════════════════════════════════════════════════════════


# Set in the parent process before Pool creation; inherited by fork()
# children so the pre-extracted sequence cache is shared via copy-on-write
# without pickling pysam handles.
_WORKER_SIM: "WholeGenomeSimulator | None" = None


def _shard_task(args):
    """Pool worker entry. Dispatches to the inherited simulator's _run_shard."""
    if _WORKER_SIM is None:
        raise RuntimeError("Worker simulator not initialized (fork inheritance failed)")
    return _WORKER_SIM._run_shard(*args)


def _shard_by_count(items, n_shards, *, weight):
    """Greedy LPT bin-packing: distribute items across shards by descending weight."""
    if n_shards <= 1:
        return [list(items)]
    sorted_items = sorted(items, key=lambda x: -weight(x))
    shards: list[list] = [[] for _ in range(n_shards)]
    loads = [0] * n_shards
    for item in sorted_items:
        idx = min(range(n_shards), key=lambda i: loads[i])
        shards[idx].append(item)
        loads[idx] += weight(item)
    return shards


def _concat_files_binary(srcs: list[Path], dst: Path) -> None:
    """Concatenate files byte-wise. Valid for gzip streams (per RFC 1952)."""
    with open(dst, "wb") as out:
        for s in srcs:
            with open(s, "rb") as f:
                shutil.copyfileobj(f, out, length=4 * 1024 * 1024)
