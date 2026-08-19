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
from typing import Callable, Sequence

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


def _mate_flags(r1_is_rev: bool, r2_is_rev: bool) -> tuple[int, int]:
    """Build the (r1_flag, r2_flag) BAM mate-orientation flags from the mate strands."""
    r1_flag = BASE_R1_FLAG
    r2_flag = BASE_R2_FLAG
    if r1_is_rev:
        r1_flag |= FLAG_REVERSE
    if r2_is_rev:
        r1_flag |= FLAG_MATE_REVERSE
        r2_flag |= FLAG_REVERSE
    if r1_is_rev:
        r2_flag |= FLAG_MATE_REVERSE
    return r1_flag, r2_flag


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

    Pipeline
    --------
    1. Pre-cache all chromosome sequences as numpy byte arrays at init,
       and pre-extract all mRNA and pre-mRNA sequences as byte arrays.
    2. Pre-compute abundance weight vectors (numpy arrays).
    3. Sample from the fragment-length distribution in bulk, counting
       fragments per unique length.
    4. For each fragment length, ONE multinomial over every RNA row — the mature row of each
       transcript and the nascent row of each nRNA ENTITY — with probability ∝ abundance ×
       capture-aware effective length. Accumulate per-row fragment counts.
    5. Iterate rows: pull sequence ONCE, extract all fragments via numpy fancy indexing, and write
       FASTQ.gz + BAM simultaneously (one RNA writer; a nascent row differs from a mature one ONLY by
       its read-name origin tag, ``nrna_<entity id>``). gzip compression is multi-threaded via pgzip
       when available.

    ⭐ **Nascent RNA is a TRANSCRIPT, not a parallel pre-mRNA space** (owner, 2026-08-19). The
    transcript list comes from the rigel index and already holds the nascent entities — single-exon
    transcripts over each clustered TSS/TES span — with the nascent molecules pooled onto them
    (`whole_genome.assign_nrna_to_entities`). Their template is their own (single-exon) sequence, the
    same one the mature writer uses, and their capture intervals are whatever probes overlap the span
    — any transcript's, any strand — exactly as for gDNA. ⛔ The previous per-transcript pre-mRNA
    space scoped a nascent molecule's capture to its OWN transcript's probes, which under-enriched
    nascent spanning another transcript's probed exon by up to 6x against the gDNA there
    (`hop_currency.py`, 2026-08-19).
    """

    def __init__(
        self,
        fasta_path: str | Path,
        transcripts: list[Transcript],
        sim_params: SimulationParams,
        gdna_config: GDNASimConfig,
        *,
        genomic_refs: Sequence[str],
        strand_specificity: float = 1.0,
        r1_sense: bool = False,
        seed: int | None = None,
        capture_config: CaptureConfig | None = None,
        capture_sampler: "CaptureSampler | None" = None,
    ):
        self.fasta = pysam.FastaFile(str(fasta_path))
        self.transcripts = transcripts
        self.sim_params = sim_params
        self.gdna_config = gdna_config
        self.strand_specificity = strand_specificity
        #: Protocol DIRECTION: False = R1-antisense (dUTP), True = R1-sense (KAPA). Independent of
        #: ``strand_specificity``, which is the FIDELITY about whichever direction is targeted.
        self.r1_sense = bool(r1_sense)
        eff_seed = seed if seed is not None else sim_params.sim_seed
        self._rng = np.random.default_rng(eff_seed)

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

        # Strand chars and abundance arrays
        self._t_strand_chars: list[str] = [
            "r" if t.strand == Strand.NEG else "f" for t in transcripts
        ]
        self._mrna_abund = np.array([t.abundance or 0.0 for t in transcripts])
        # ⛔ nascent lives on single-exon ENTITY rows only: a multi-exon row carrying nascent would be
        # sampled on a SPLICED template and written as contiguous genomic reads — wrong both ways.
        self._nrna_abund = np.array([t.nrna_abundance for t in transcripts])
        bad = [t.t_id for t in transcripts if t.nrna_abundance > 0 and len(t.exons) != 1]
        if bad:
            raise ValueError(
                f"{len(bad)} multi-exon transcripts carry nrna_abundance > 0 (e.g. {bad[:3]}); "
                "nascent RNA must be pooled onto single-exon nRNA entities "
                "(whole_genome.assign_nrna_to_entities)"
            )

        # Quality string cache
        self._quals_cache: dict[int, str] = {}

        # Fast-path flag for error introduction
        self._has_errors = sim_params.error_rate > 0

        # Clamp frag_max to longest transcript
        max_mrna = int(self._t_lengths.max()) if N > 0 else 0
        self._frag_max = min(sim_params.frag_max, max_mrna) if max_mrna > 0 else sim_params.frag_max

        # ⭐ gDNA setup: the GENOMIC references, stated by the caller, weighted by length.
        # ⛔ This used to be `annotated_refs = {t.ref for t in transcripts}` — "has an annotation"
        # standing in for "is genomic", and it is not one. An RNA-only spike-in reference carries a
        # transcript, so it qualified, and the panel filled with gDNA molecules on templates where no
        # genomic DNA exists. Every reference is genomic or RNA-only and the classification is an
        # INPUT: `tests/test_sim_genomic_refs.py` pins that, with a fixture whose RNA-only reference
        # is annotated precisely so the old proxy cannot pass it.
        unknown = sorted(set(genomic_refs) - set(self.fasta.references))
        if unknown:
            raise ValueError(
                f"genomic_refs not in the reference FASTA: {unknown}. "
                "A mis-named genomic reference must fail loudly, not silently emit no gDNA."
            )
        wanted = set(genomic_refs)
        #: The references gDNA is drawn from, in FASTA order — so reference ids stay a subsequence
        #: of the production ordering rather than a reshuffle.
        self.genomic_refs: tuple[str, ...] = tuple(
            ref for ref in self.fasta.references if ref in wanted
        )
        self._gdna_refs: list[str] = list(self.genomic_refs)
        self._gdna_ref_lengths: list[int] = [
            self.fasta.get_reference_length(ref) for ref in self._gdna_refs
        ]

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
        self._ref_name_to_id: dict[str, int] = {name: i for i, name in enumerate(self._ref_names)}
        # ⭐ A prebuilt sampler is REUSED ACROSS CONDITIONS: its probe layout and its partition memo
        # depend only on the panel and the templates, never on abundances (`CaptureSampler`), and the
        # partition is 260 s per condition otherwise.
        self.capture = capture_sampler or CaptureSampler.from_config(
            capture_config,
            transcripts,
            dict(zip(self._ref_names, self._ref_lengths)),
        )
        self._bam_header = pysam.AlignmentHeader.from_dict(
            {
                "HD": {"VN": "1.6", "SO": "unsorted", "GO": "query"},
                "SQ": [
                    {"SN": name, "LN": length}
                    for name, length in zip(self._ref_names, self._ref_lengths)
                ],
                "PG": [
                    {
                        "ID": "rigel_sim",
                        "PN": "rigel_sim",
                        "VN": "1.0",
                        "CL": "simulated",
                    }
                ],
            }
        )

        # Pre-built shared tag lists for BAM records
        self._nh1_tags = [("NH", 1)]

        logger.info(
            "Simulator ready: %d transcripts, %d gDNA refs, max mRNA len=%d",
            N,
            len(self._gdna_refs),
            max_mrna,
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

    # -- Post-capture length marginal ---------------------------------------

    def _post_capture_length_allocation(
        self,
        drawn_lengths: np.ndarray,
        weights_at: Callable[[int], np.ndarray],
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, list[np.ndarray], np.ndarray]:
        """Re-allocate a pre-capture length draw in proportion to capture-weighted opportunity.

        ⭐⭐ **THE ONE PLACE `f_post(w) = f_pre(w) . total_eff(w) / Z` IS WRITTEN**, shared by the
        mature, nascent and gDNA pools so the three cannot drift apart.

        ⛔ **The defect this replaces.** Each pool drew its length marginal FIRST, capture-blind, then
        computed the capture-aware per-template opportunity ``weights_at(w)`` and threw its TOTAL away
        by normalising within each length. Capture could then move only *where* a fragment landed,
        never *whether* its length survived, and the simulated post-capture fragment-length
        distribution came out byte-identical to the pre-capture one on every condition. Hybrid capture
        hybridises probes to sequence: a short fragment presents less sequence, binds worse, and is
        captured less efficiently — so capture SELECTS FOR LENGTH, and that selection lives entirely in
        the term that was being divided out.

        ⭐ **No new constant.** ``total_eff(w)`` is already computed by machinery that already exists:
        ``CaptureSampler.fragment_weight`` is ``off_target_weight + binding_per_base * overlap``, and
        ``overlap`` rises with fragment length until it saturates at the probe length. That is the
        physics, already parameterised by ``CaptureConfig``.

        ⚠ **It reweights off capture too, and that is also a correction.** With no probes,
        ``total_eff(w) = off_target_weight * sum_k a_k (L_k - w + 1)+`` — the ordinary effective
        length. A library CAN'T yield more fragments of a length than its templates have placements
        for it, and the old code let it. On a whole chromosome the term is flat to 1 part in 10^5, so
        the gDNA pool barely moves off capture; on transcripts of a few kb it is a real tilt.

        ⚠ **Two-stage, not analytic, on purpose.** ``f_pre`` stays defined by exactly one thing — the
        shared sampler in :mod:`sampling` — and its empirical draw is reweighted, rather than a second
        analytic copy of the truncated normal being written here for the marginal to disagree with
        later. The cost is a ratio estimator's O(1/n) bias, which at the pilot's 5 M fragments over
        ~800 lengths is a relative 1e-7, and a sqrt(2) inflation of the realised mean's Monte-Carlo
        sd: 0.063 bp instead of 0.045 bp.

        Returns ``(widths, weights_per_width, counts_per_width)`` — the ascending distinct lengths,
        each one's unnormalised per-template weight vector, and how many fragments it now carries.
        The total is preserved exactly: ``counts.sum() == len(drawn_lengths)``.
        """
        widths, drawn_counts = np.unique(drawn_lengths, return_counts=True)
        weights = [weights_at(int(width)) for width in widths]
        totals = np.array([float(vector.sum()) for vector in weights])
        posterior = drawn_counts * totals
        total = float(posterior.sum())
        if total <= 0.0:
            return widths, weights, np.zeros(len(widths), dtype=np.int64)
        counts = rng.multinomial(int(drawn_counts.sum()), posterior / total)
        logger.debug(
            "post-capture length marginal: mean %.2f -> %.2f over %d lengths",
            float(np.average(widths, weights=drawn_counts)),
            float(np.average(widths, weights=counts)) if counts.sum() else float("nan"),
            len(widths),
        )
        return widths, weights, counts

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

        counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))

        # ⭐ Only transcripts with nonzero abundance can carry a fragment: `weights` multiplies the
        # capture-aware effective length by the abundance, so a zero-abundance row contributes zero
        # whatever its effective length is. With `frac_expressed: 0.5` that is HALF the annotation, and
        # the effective length is the expensive term — `partition_array` was 96.6 % of a capture-on run.
        # ⚠ `eff` is still built at FULL length with zeros in the dead rows, so `weights`, `probs` and
        # therefore the `rng.choice` draw are bit-identical to computing every row. This is a speed
        # change, not a behaviour change, and `tests/test_sim_capture_partition.py` pins that.
        live = np.flatnonzero(abundances > 0)
        live_keys = live.tolist()
        live_lengths = lengths[live]

        def weights_at(width: int) -> np.ndarray:
            eff = np.zeros(len(abundances), dtype=np.float64)
            if live.size:
                eff[live] = self.capture.partition_array(space, live_keys, live_lengths, width)
            return abundances * eff

        # ⭐ The length marginal is the pre-capture draw reweighted by capture-weighted opportunity;
        # `weights_at(w)` is then the conditional over templates AT that length. See
        # `_post_capture_length_allocation` for why the total must not be normalised away.
        widths, weights_per_width, counts_per_width = self._post_capture_length_allocation(
            frag_lengths, weights_at, rng
        )

        for width, weights, fc in zip(widths, weights_per_width, counts_per_width):
            fl, fc = int(width), int(fc)
            if fc <= 0:
                continue
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
        self, n_rna: int
    ) -> tuple[dict[int, dict[int, int]], dict[int, dict[int, int]]]:
        """ONE multinomial over every RNA row: the mature row of each transcript (``abundance``) and
        the nascent row of each nRNA entity (``nrna_abundance``), each weighted by its abundance ×
        capture-aware effective length on its OWN template. The nascent fragment share is therefore a
        consequence of molecules and lengths, never an imposed count.

        Returns ``(mrna_counts, nrna_counts)``, each ``dict[t_idx, dict[frag_len, count]]`` — two
        dicts only because the writer tags the two rows' read names differently.
        """
        if n_rna <= 0:
            return {}, {}

        N = len(self.transcripts)
        rng = self._rng
        frag_lengths = self._sample_rna_frag_lengths(n_rna)

        mrna_counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))
        nrna_counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))

        # ⭐ Only rows with nonzero abundance can carry a fragment, and the effective length is the
        # expensive term (`partition_array` was 96.6 % of a capture-on run), so it is evaluated on the
        # LIVE rows only; `eff` is still full length with zeros in the dead rows, so the draw is
        # bit-identical to computing every row (`tests/test_sim_capture_partition.py`).
        live = np.flatnonzero((self._mrna_abund > 0) | (self._nrna_abund > 0))
        live_keys = live.tolist()
        live_lengths = self._t_lengths[live]

        def weights_at(width: int) -> np.ndarray:
            eff = np.zeros(N, dtype=np.float64)
            if live.size:
                eff[live] = self.capture.partition_array("mrna", live_keys, live_lengths, width)
            return np.concatenate([self._mrna_abund * eff, self._nrna_abund * eff])

        widths, weights_per_width, counts_per_width = self._post_capture_length_allocation(
            frag_lengths, weights_at, rng
        )

        for width, weights, fc in zip(widths, weights_per_width, counts_per_width):
            fl, fc = int(width), int(fc)
            if fc <= 0:
                continue
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
        """Generate and write all reads of one RNA row — a transcript's mature row, or an nRNA
        entity's nascent row (``is_nrna``), which differs ONLY by its read-name origin tag: the
        template is the transcript's own sequence in both cases (an entity is single-exon, so its
        sequence IS its genomic span) and the BAM blocks come from the same exon structure.
        """
        t = self.transcripts[t_idx]
        rng = self._rng
        ss = self.strand_specificity
        strand_char = self._t_strand_chars[t_idx]
        t_id = t.t_id
        if is_nrna and len(t.exons) != 1:
            raise ValueError(f"nascent row on a multi-exon transcript {t_id}")

        seq_bytes = self._t_seq_bytes[t_idx]
        seq_len = int(self._t_lengths[t_idx])
        qname_prefix = f"nrna_{t_id}" if is_nrna else t_id

        ref_id = self._ref_name_to_id.get(t.ref) if bam_fh else None
        n_written = 0

        for frag_len in sorted(len_counts):
            count = len_counts[frag_len]
            read_len = min(self.sim_params.read_length, frag_len)
            eff_len = seq_len - frag_len + 1
            if eff_len <= 0:
                continue

            frag_starts = self.capture.sample_starts(
                "mrna",
                t_idx,
                seq_len,
                frag_len,
                count,
                rng,
            )
            # The base emission is R1-ANTISENSE; a flip makes that fragment R1-sense. ``ss`` is the
            # fidelity about the TARGETED direction, so an R1-sense protocol flips the fragments the
            # R1-antisense protocol would have kept — the exact per-fragment mirror on one RNG stream.
            # ⚠ The r1_sense=False path must consume the RNG exactly as before, or every existing
            # condition's realized library moves.
            if ss < 1.0:
                flip_mask = rng.random(count) >= ss
                if self.r1_sense:
                    flip_mask = ~flip_mask
            elif self.r1_sense:
                flip_mask = np.ones(count, dtype=bool)
            else:
                flip_mask = None
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
                self._write_mrna_bam_batch(
                    bam_fh,
                    qnames,
                    r1_seqs,
                    r2_seqs,
                    t,
                    frag_starts,
                    frag_len,
                    read_len,
                    flip_mask,
                    ref_id,
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

            r1_flag, r2_flag = _mate_flags(r1_is_rev, r2_is_rev)

            r1_tlen = tlen if r1_start <= r2_start else -tlen
            r2_tlen = -r1_tlen

            r1_seq_str = _bam_seq_from_fastq_bytes(r1_seqs[i], r1_is_rev)
            r2_seq_str = _bam_seq_from_fastq_bytes(r2_seqs[i], r2_is_rev)

            bam_fh.write(
                make_aligned_segment(
                    self._bam_header,
                    qnames[i],
                    r1_seq_str,
                    r1_flag,
                    ref_id,
                    r1_start,
                    r1_cigar,
                    ref_id,
                    r2_start,
                    r1_tlen,
                    tags=tags,
                )
            )
            bam_fh.write(
                make_aligned_segment(
                    self._bam_header,
                    qnames[i],
                    r2_seq_str,
                    r2_flag,
                    ref_id,
                    r2_start,
                    r2_cigar,
                    ref_id,
                    r1_start,
                    r2_tlen,
                    tags=tags,
                )
            )

    def _accumulate_gdna_counts(self, n_gdna: int) -> dict[tuple[int, int], int]:
        """Sample gDNA fragment lengths + chromosomes.

        Returns dict[(ref_idx, frag_len)] = count.
        """
        if n_gdna == 0:
            return {}
        if not self._gdna_refs:
            # ⛔ Writing zero of a requested five million fragments is trap 20's failure mode: the
            # run completes, the truth files agree with themselves, and the deficit is invisible.
            raise ValueError(
                f"{n_gdna} gDNA fragments requested but no genomic reference was declared "
                "(genomic_refs is empty)"
            )

        rng = self._rng
        frag_lengths = self._sample_gdna_frag_lengths(n_gdna)
        counts: dict[tuple[int, int], int] = {}
        # ⭐ One batched call per fragment length, not one scalar call per (chromosome, length).
        # Profiled: this comprehension was 49,662 `partition` calls driving 6,476,550
        # `_local_overlap_weights` calls and 110.6 s, because every call re-integrates the capture
        # landscape of a WHOLE chromosome. Batched, the 93 references share one pass.
        ref_lengths = np.asarray(self._gdna_ref_lengths, dtype=np.int64)
        gdna_refs = list(self._gdna_refs)

        def weights_at(width: int) -> np.ndarray:
            return self.capture.partition_array("gdna", gdna_refs, ref_lengths, width)

        # ⭐ The per-chromosome conditional was always right; the length MARGINAL is what was thrown
        # away. Off capture this term is flat to 1 part in 10^5 on a whole chromosome — which is why
        # correcting it moves the gDNA pool essentially not at all off capture, and a great deal
        # under it. See `_post_capture_length_allocation`.
        widths, weights_per_width, counts_per_width = self._post_capture_length_allocation(
            frag_lengths, weights_at, rng
        )

        for width, chrom_eff, fc in zip(widths, weights_per_width, counts_per_width):
            fl, fc = int(width), int(fc)
            if fc <= 0:
                continue
            total_eff = chrom_eff.sum()
            if total_eff <= 0:
                continue
            chrom_probs = chrom_eff / total_eff
            chrom_indices = rng.choice(len(gdna_refs), size=fc, p=chrom_probs)
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

                r1_flag, r2_flag = _mate_flags(r1_is_rev, r2_is_rev)

                tlen = end - start
                r1_tlen = tlen if r1_start_pos <= r2_start_pos else -tlen
                r2_tlen = -r1_tlen

                r1_seq_str = _bam_seq_from_fastq_bytes(r1_seqs[j], r1_is_rev)
                r2_seq_str = _bam_seq_from_fastq_bytes(r2_seqs[j], r2_is_rev)

                bam_fh.write(
                    make_aligned_segment(
                        self._bam_header,
                        qnames[j],
                        r1_seq_str,
                        r1_flag,
                        ref_id,
                        r1_start_pos,
                        [(pysam.CMATCH, read_len)],
                        ref_id,
                        r2_start_pos,
                        r1_tlen,
                        tags=self._nh1_tags,
                    )
                )
                bam_fh.write(
                    make_aligned_segment(
                        self._bam_header,
                        qnames[j],
                        r2_seq_str,
                        r2_flag,
                        ref_id,
                        r2_start_pos,
                        [(pysam.CMATCH, read_len)],
                        ref_id,
                        r1_start_pos,
                        r2_tlen,
                        tags=self._nh1_tags,
                    )
                )
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
                ref_idx,
                fl,
                count,
                r1_buf,
                r2_buf,
                bam_fh,
                n_offset=n_written,
            )
        return n_written

    # -- Abundance-weighted pool splits (single-condition / Scenario use) ----

    def _rna_weight(self) -> float:
        """Σ over every RNA row of abundance × capture-effective length at the mean RNA fragment
        length — mature rows and nascent entity rows on the same templates."""
        mean_frag = int(self.sim_params.frag_mean)
        n = len(self.transcripts)
        eff = self.capture.partition_array("mrna", range(n), self._t_lengths, mean_frag)
        return float(np.sum((self._mrna_abund + self._nrna_abund) * eff))

    def rna_gdna_split(self, n_total: int, gdna_abundance: float) -> tuple[int, int]:
        """Split ``n_total`` fragments into ``(n_rna, n_gdna)`` by abundance × effective length, with
        gDNA weighted by ``gdna_abundance × genome effective length`` (summed over the genomic
        references). The mature/nascent split INSIDE ``n_rna`` is not imposed — it is what the one
        multinomial over the RNA rows realises."""
        rna_w = self._rna_weight()
        gdna_mean_frag = int(self.gdna_config.frag_mean)
        genome_eff = sum(
            self.capture.partition("gdna", ref, length, gdna_mean_frag)
            for ref, length in zip(self._gdna_refs, self._gdna_ref_lengths)
        )
        gdna_w = float(gdna_abundance) * genome_eff
        total = rna_w + gdna_w
        if total <= 0:
            return 0, 0
        n_gdna = int(round(n_total * gdna_w / total))
        return max(0, n_total - n_gdna), n_gdna

    # -- Main entry point ---------------------------------------------------

    def simulate_and_write(
        self,
        output_dir: Path,
        n_rna: int,
        n_gdna: int = 0,
        *,
        oracle_bam: bool = True,
        prefix: str = "sim",
        n_workers: int = 1,
    ) -> tuple[Path, Path, Path | None]:
        """Single-pass simulation: accumulate counts, generate, write.

        *n_rna* fragments are drawn from ONE multinomial over every RNA row (mature rows and nascent
        entity rows, :meth:`_accumulate_rna_counts`); *n_gdna* from the genomic references.

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
        mrna_counts, nrna_counts = self._accumulate_rna_counts(n_rna)
        gdna_counts = self._accumulate_gdna_counts(n_gdna)

        total_mrna = sum(sum(d.values()) for d in mrna_counts.values())
        total_nrna = sum(sum(d.values()) for d in nrna_counts.values())
        total_gdna = sum(gdna_counts.values())
        logger.info(
            "Fragment counts: %d mRNA (%d txs), %d nRNA (%d txs), %d gDNA chunks",
            total_mrna,
            len(mrna_counts),
            total_nrna,
            len(nrna_counts),
            len(gdna_counts),
        )

        if n_workers <= 1:
            self._write_all_serial(
                r1_path,
                r2_path,
                bam_path,
                mrna_counts,
                nrna_counts,
                gdna_counts,
                oracle_bam,
            )
        else:
            self._write_all_parallel(
                output_dir,
                prefix,
                r1_path,
                r2_path,
                bam_path,
                mrna_counts,
                nrna_counts,
                gdna_counts,
                oracle_bam,
                n_workers,
            )

        n_written = total_mrna + total_nrna + total_gdna
        logger.info(
            "Wrote %d read pairs -> %s (RNA=%d, gDNA=%d, oracle=%s, workers=%d)",
            n_written,
            output_dir,
            n_rna,
            n_gdna,
            oracle_bam,
            n_workers,
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
                    str(bam_path),
                    "wb",
                    header=self._bam_header,
                )
            try:
                for t_idx in sorted(mrna_counts):
                    self._write_rna_reads(
                        t_idx,
                        mrna_counts[t_idx],
                        r1_buf,
                        r2_buf,
                        bam_fh,
                        is_nrna=False,
                    )
                for t_idx in sorted(nrna_counts):
                    self._write_rna_reads(
                        t_idx,
                        nrna_counts[t_idx],
                        r1_buf,
                        r2_buf,
                        bam_fh,
                        is_nrna=True,
                    )
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
            mrna_items,
            n_workers,
            weight=lambda x: sum(x[1].values()),
        )
        nrna_shards = _shard_by_count(
            nrna_items,
            n_workers,
            weight=lambda x: sum(x[1].values()),
        )
        gdna_shards = _shard_by_count(
            gdna_items,
            n_workers,
            weight=lambda x: x[1],
        )

        shard_dir = output_dir / f".{prefix}_shards"
        shard_dir.mkdir(parents=True, exist_ok=True)

        base_seed = int(self._rng.integers(0, 2**31 - 1))

        tasks = []
        for k in range(n_workers):
            tasks.append(
                (
                    k,
                    mrna_shards[k],
                    nrna_shards[k],
                    gdna_shards[k],
                    str(shard_dir),
                    prefix,
                    oracle_bam,
                    base_seed + k,
                )
            )

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
        # Independent RNG per shard.
        self._rng = np.random.default_rng(seed)

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
                    str(bam_path),
                    "wb",
                    header=self._bam_header,
                )
            try:
                for t_idx, len_counts in mrna_items:
                    self._write_rna_reads(
                        t_idx,
                        len_counts,
                        r1_buf,
                        r2_buf,
                        bam_fh,
                        is_nrna=False,
                    )
                for t_idx, len_counts in nrna_items:
                    self._write_rna_reads(
                        t_idx,
                        len_counts,
                        r1_buf,
                        r2_buf,
                        bam_fh,
                        is_nrna=True,
                    )
                n_offset = 0
                for (ref_idx, fl), count in gdna_items:
                    n_offset += self._write_gdna_chunk(
                        ref_idx,
                        fl,
                        count,
                        r1_buf,
                        r2_buf,
                        bam_fh,
                        n_offset=n_offset,
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
