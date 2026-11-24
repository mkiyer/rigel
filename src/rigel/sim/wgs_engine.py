"""The read-simulation engine behind every simulated library.

``WholeGenomeSimulator`` turns a genome FASTA and a transcript list with abundances assigned into paired-end FASTQ
and, optionally, a query-grouped oracle BAM whose read names carry each fragment's origin (:mod:`.read_name`).
``simulate_and_write`` first draws every fragment count, RNA from one multinomial over the transcript and nascent
rows and gDNA from one over the genomic references, at each fragment length; then it writes the reads, serially or
in forked shards. The run is a function of the seed: each shard's seed comes from the same stream, so a given
worker count reproduces exactly.
"""

from __future__ import annotations

import contextlib
import gzip
import logging
import multiprocessing
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Callable, NamedTuple, Sequence

import numpy as np
import pysam

try:
    import pgzip
except ImportError:
    pgzip = None  # type: ignore[assignment]

from rigel import __version__ as rigel_version
from rigel.config import CONSTANTS
from rigel.transcript import Transcript
from rigel.types import Strand
from rigel.sim.bam import (
    BASE_QUALITY,
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
from rigel.sim.read_name import STRAND_FORWARD, STRAND_REVERSE, format_gdna_name, format_rna_name
from rigel.sim.sampling import truncated_normal_frag_lengths
from rigel.sim.wgs_config import GDNASimConfig, SimulationParams

logger = logging.getLogger(__name__)


_BYTE_COMPLEMENT = np.zeros(256, dtype=np.uint8)
for _c, _rc in zip(b"ACGTNacgtn", b"TGCANtgcan"):
    _BYTE_COMPLEMENT[_c] = _rc

#: The four bases in byte order, so ``searchsorted`` reads a base's index off its byte.
_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


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
    """Extract every fragment's reads at once: R2 is the first ``read_len`` bases, R1 the reverse complement of the
    last ``read_len``. Returns ``(r2_seqs, r1_seqs)``, uint8 arrays of shape ``(count, read_len)``.
    """
    offsets = np.arange(read_len, dtype=np.int64)

    r2_indices = frag_starts[:, np.newaxis] + offsets[np.newaxis, :]
    r2_seqs = seq_bytes[r2_indices]

    r1_start = frag_starts + (frag_len - read_len)
    r1_indices = r1_start[:, np.newaxis] + offsets[np.newaxis, :]
    r1_seqs = _BYTE_COMPLEMENT[seq_bytes[r1_indices][:, ::-1]]

    return r2_seqs, r1_seqs


def _bam_seq_from_fastq_bytes(seq: np.ndarray, is_reverse: bool) -> str:
    """Return the SAM/BAM stored SEQ for a FASTQ-oriented read sequence."""
    if is_reverse:
        seq = _BYTE_COMPLEMENT[seq[::-1]]
    return seq.tobytes().decode("ascii")


def _open_gzip_write(path: Path, threads: int):
    """Open a gzip file for text writing at the simulator's gzip level: with ``threads`` pgzip
    compressors when pgzip is installed, the standard library's single-threaded gzip otherwise."""
    level = CONSTANTS.resources.fastq_gzip_level
    if pgzip is not None:
        return pgzip.open(path, "wt", thread=threads, compresslevel=level)
    return gzip.open(path, "wt", compresslevel=level)


def _bam_header(ref_names: list[str], ref_lengths: list[int]) -> pysam.AlignmentHeader:
    """The oracle BAM's header: query-grouped, every reference of the FASTA."""
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "unsorted", "GO": "query"},
            "SQ": [{"SN": name, "LN": length} for name, length in zip(ref_names, ref_lengths)],
            "PG": [{"ID": "rigel_sim", "PN": "rigel_sim", "VN": rigel_version, "CL": "simulated"}],
        }
    )


class _FastqPair:
    """One writer's two FASTQ files, buffered: ``fastq_buffer_records`` records are joined into one gzip write per
    file."""

    __slots__ = ("_handles", "_parts", "_count")

    def __init__(self, r1_fh, r2_fh):
        self._handles = (r1_fh, r2_fh)
        self._parts: tuple[list[str], list[str]] = ([], [])
        self._count = 0

    def write(
        self, qnames: list[str], r1_seqs: np.ndarray, r2_seqs: np.ndarray, quals: str
    ) -> None:
        """One record per name in each file: the sequences are ``count x read_len`` uint8, ``quals`` every read's
        quality string."""
        for parts, seqs, suffix in zip(self._parts, (r1_seqs, r2_seqs), ("/1", "/2")):
            for qname, seq in zip(qnames, seqs):
                parts.append(f"@{qname}{suffix}\n{seq.tobytes().decode('ascii')}\n+\n{quals}\n")
        self._count += len(qnames)
        if self._count >= CONSTANTS.resources.fastq_buffer_records:
            self.flush()

    def flush(self) -> None:
        for fh, parts in zip(self._handles, self._parts):
            if parts:
                fh.write("".join(parts))
                parts.clear()
        self._count = 0


class _Mate(NamedTuple):
    """One read of a pair as its BAM record needs it: the sequence as the FASTQ holds it, its substitution count,
    its leftmost position, its CIGAR and whether it maps to the reverse strand."""

    seq: np.ndarray
    nm: int
    start: int
    cigar: list[tuple[int, int]]
    reverse: bool


# ═══════════════════════════════════════════════════════════════════
# Genome sequence cache
# ═══════════════════════════════════════════════════════════════════


class GenomeCache:
    """Whole chromosomes as uppercase uint8 arrays, loaded on first access, so gDNA fragments need no per-fragment
    ``fetch``."""

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

    def preload(self, refs: Sequence[str]) -> None:
        for ref in refs:
            self.get(ref)

    def clear(self) -> None:
        self._cache.clear()


# ═══════════════════════════════════════════════════════════════════
# Whole-genome read simulator
# ═══════════════════════════════════════════════════════════════════


class WholeGenomeSimulator:
    """Whole-genome paired-end RNA-seq read simulator, writing gzipped FASTQ and a query-grouped oracle BAM.

    Every transcript's spliced sequence is extracted once at init and the genomic references are preloaded. Per
    fragment length, one multinomial over every RNA row draws counts in proportion to abundance × capture-aware
    effective length; each row's fragments are then extracted at once and written to FASTQ and BAM together.

    Nascent RNA is a transcript here: the index's nascent entities are single-exon rows over each TSS/TES span,
    carrying the pooled nascent molecules (`transcriptome.assign_nrna_to_entities`). An entity's template is its own
    sequence and its capture intervals are whatever probes overlap its span, as for gDNA; a nascent row differs
    from a transcript row only by its read-name tag.
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
        capture_config: CaptureConfig | None = None,
        capture_sampler: "CaptureSampler | None" = None,
    ):
        self.fasta = pysam.FastaFile(str(fasta_path))
        self.transcripts = transcripts
        self.sim_params = sim_params
        self.gdna_config = gdna_config
        self.strand_specificity = strand_specificity
        #: Protocol direction: False is R1-antisense (dUTP), True R1-sense. ``strand_specificity`` is the fidelity
        #: about whichever direction is targeted.
        self.r1_sense = bool(r1_sense)
        self._rng = np.random.default_rng(sim_params.sim_seed)
        self._genome = GenomeCache(self.fasta)
        self._load_templates()

        self._quals_cache: dict[int, str] = {}
        self._has_errors = sim_params.error_rate > 0
        # frag_max is clamped to the longest transcript
        max_mrna = int(self._t_lengths.max()) if len(transcripts) > 0 else 0
        self._frag_max = min(sim_params.frag_max, max_mrna) if max_mrna > 0 else sim_params.frag_max

        self._load_genomic_refs(genomic_refs)

        ref_names = list(self.fasta.references)
        ref_lengths = [self.fasta.get_reference_length(r) for r in ref_names]
        self._ref_name_to_id: dict[str, int] = {name: i for i, name in enumerate(ref_names)}
        # A prebuilt sampler is shared across conditions: its probe layout, its probes' strengths and its partition
        # memo depend only on the panel, the reference and the templates.
        self.capture = capture_sampler or CaptureSampler.from_config(
            capture_config, transcripts, dict(zip(ref_names, ref_lengths)), fasta=self.fasta
        )
        # Bases the reference does not know (`N`) are never drawn: each template is masked at its unknown runs.
        self.capture.mask_unknown_sequence(self.fasta, self.genomic_refs)
        self._bam_header = _bam_header(ref_names, ref_lengths)
        self._nh1_tags = [("NH", 1)]

        logger.info(
            "Simulator ready: %d transcripts, %d gDNA refs, max sequence length=%d",
            len(transcripts),
            len(self.genomic_refs),
            max_mrna,
        )

    # -- Setup ----------------------------------------------------------------

    def _load_templates(self) -> None:
        """Every RNA row's template (its spliced sequence; an entity's is its span), length, strand letter and
        abundances."""
        transcripts = self.transcripts
        logger.info("Pre-extracting %d transcript and nascent-span sequences...", len(transcripts))
        self._t_seq_bytes: list[np.ndarray] = []
        self._t_lengths = np.empty(len(transcripts), dtype=np.int64)
        for i, t in enumerate(transcripts):
            seq = self._spliced_sequence(t)
            self._t_seq_bytes.append(_seq_to_bytes(seq))
            self._t_lengths[i] = len(seq)

        self._t_strand_chars: list[str] = [
            STRAND_REVERSE if t.strand == Strand.NEG else STRAND_FORWARD for t in transcripts
        ]
        self._mrna_abund = np.array([t.abundance or 0.0 for t in transcripts])
        # Nascent must live on single-exon entity rows: on a multi-exon row it would be drawn from a spliced
        # template.
        self._nrna_abund = np.array([t.nrna_abundance for t in transcripts])
        bad = [t.t_id for t in transcripts if t.nrna_abundance > 0 and len(t.exons) != 1]
        if bad:
            raise ValueError(
                f"{len(bad)} multi-exon transcripts carry nrna_abundance > 0 (e.g. {bad[:3]}); "
                "nascent RNA must be pooled onto single-exon nRNA entities "
                "(transcriptome.assign_nrna_to_entities)"
            )

    def _load_genomic_refs(self, genomic_refs: Sequence[str]) -> None:
        """The references gDNA is drawn from, preloaded, and their strand regions.

        They are the caller's, never inferred from the annotation: an RNA-only spike-in reference carries a
        transcript but no genomic DNA. A mis-named one raises.
        """
        unknown = sorted(set(genomic_refs) - set(self.fasta.references))
        if unknown:
            raise ValueError(
                f"genomic_refs not in the reference FASTA: {unknown}. "
                "A mis-named genomic reference must fail loudly, not silently emit no gDNA."
            )
        wanted = set(genomic_refs)
        #: The references gDNA is drawn from, in FASTA order.
        self.genomic_refs: tuple[str, ...] = tuple(
            ref for ref in self.fasta.references if ref in wanted
        )
        self._gdna_ref_lengths: list[int] = [
            self.fasta.get_reference_length(ref) for ref in self.genomic_refs
        ]
        if self.genomic_refs:
            logger.info("Pre-loading %d chromosomes for gDNA...", len(self.genomic_refs))
            self._genome.preload(self.genomic_refs)

        # gDNA strand overdispersion: exon-bounded regions, each with a shared +strand rate ~ Beta(a, a)
        # (`_init_gdna_strand_regions`); at 0 no regions are built and every split is 50/50.
        self._gdna_strand_regions: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        gdna_od = float(self.gdna_config.strand_overdispersion)
        if gdna_od > 0.0 and self.genomic_refs:
            self._init_gdna_strand_regions(gdna_od)

    def _spliced_sequence(self, t: Transcript) -> str:
        """A row's spliced sequence, 5' to 3'."""
        seq = "".join(self.fasta.fetch(t.ref, e.start, e.end).upper() for e in t.exons)
        if t.strand == Strand.NEG:
            seq = reverse_complement(seq)
        return seq

    def _get_quals(self, read_len: int) -> str:
        """Return cached quality string of the given length."""
        q = self._quals_cache.get(read_len)
        if q is None:
            q = BASE_QUALITY * read_len
            self._quals_cache[read_len] = q
        return q

    # -- Fragment length sampling -------------------------------------------

    def _sample_rna_frag_lengths(self, n: int) -> np.ndarray:
        p = self.sim_params
        return truncated_normal_frag_lengths(
            self._rng, n, p.frag_mean, p.frag_std, p.frag_min, self._frag_max
        )

    def _init_gdna_strand_regions(self, overdispersion: float) -> None:
        """Split each gDNA reference at the exon boundaries into regions, each with a shared +strand rate
        ~ Beta(a, a). Beta(a, a) has intra-class correlation ``1/(2a + 1)``, so ``a = ½(1 − od)/od`` makes the
        per-region strand count Beta-Binomial(½, od), the law the calibrator's gDNA strand model fits.
        """
        alpha = 0.5 * (1.0 - overdispersion) / overdispersion
        exons_by_ref: dict[str, set[int]] = defaultdict(set)
        for t in self.transcripts:
            if t.ref is None:
                continue
            for e in t.exons:
                exons_by_ref[t.ref].add(int(e.start))
                exons_by_ref[t.ref].add(int(e.end))
        for ref, length in zip(self.genomic_refs, self._gdna_ref_lengths):
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
        """Substitute each base of a ``count x read_len`` uint8 batch in place, with probability ``error_rate``, by one
        of the three other bases, so the realised rate is the configured one. Returns each read's substitution
        count, its NM in the oracle BAM. With no errors configured it draws nothing from the RNG."""
        count, read_len = seqs.shape
        if not self._has_errors:
            return np.zeros(count, dtype=np.int64)
        mask = self._rng.random((count, read_len)) < self.sim_params.error_rate
        n_errors = int(mask.sum())
        if n_errors > 0:
            own = np.searchsorted(_BASES, seqs[mask]) % len(_BASES)
            shift = self._rng.integers(1, len(_BASES), size=n_errors)
            seqs[mask] = _BASES[(own + shift) % len(_BASES)]
        return mask.sum(axis=1)

    # -- Post-capture length marginal ---------------------------------------

    def _post_capture_length_allocation(
        self,
        drawn_lengths: np.ndarray,
        weights_at: Callable[[int], np.ndarray],
        rng: np.random.Generator,
        pool: str,
    ) -> tuple[np.ndarray, list[np.ndarray], np.ndarray]:
        """Re-allocate a pre-capture length draw in proportion to capture-weighted opportunity:
        ``f_post(w) = f_pre(w) · total_eff(w) / Z``, the one place this is written, shared by every pool.

        The total of ``weights_at(w)`` must not be normalised away: capture binds a short fragment worse, so it
        selects for length, and that selection lives in this total. With no probes ``total_eff(w)`` is the ordinary
        effective length ``off_target_weight * Σ_k a_k (L_k − w + 1)+``, nearly flat on a whole chromosome and a
        real tilt on transcripts. ``f_pre`` stays the empirical draw of :mod:`sampling`, so there is one definition
        of the pre-capture law.

        Returns ``(widths, weights_per_width, counts_per_width)``: the ascending distinct lengths, each one's
        unnormalised per-template weights, and its new fragment count; ``counts.sum() == len(drawn_lengths)``.
        Raises when no template of the ``pool`` named can hold any drawn length: writing none of the requested
        fragments must raise, or the truth would agree with itself and hide it.
        """
        widths, drawn_counts = np.unique(drawn_lengths, return_counts=True)
        weights = [weights_at(int(width)) for width in widths]
        totals = np.array([float(vector.sum()) for vector in weights])
        posterior = drawn_counts * totals
        total = float(posterior.sum())
        if total <= 0.0:
            raise ValueError(
                f"{len(drawn_lengths)} {pool} fragments requested but no {pool} template can hold one: every "
                f"template's weight is 0 at every drawn length ({int(widths.min())}–{int(widths.max())} bp), as when "
                "no row has an abundance or every template is shorter than the fragments"
            )
        counts = rng.multinomial(int(drawn_counts.sum()), posterior / total)
        logger.debug(
            "post-capture length marginal: mean %.2f -> %.2f over %d lengths",
            float(np.average(widths, weights=drawn_counts)),
            float(np.average(widths, weights=counts)) if counts.sum() else float("nan"),
            len(widths),
        )
        return widths, weights, counts

    # -- Fragment count accumulation ----------------------------------------

    def _accumulate_rna_counts(
        self, n_rna: int
    ) -> tuple[dict[int, dict[int, int]], dict[int, dict[int, int]]]:
        """Draw ``n_rna`` fragments from one multinomial over every RNA row, each transcript's ``abundance`` and each
        entity's ``nrna_abundance``, weighted by capture-aware effective length on its own template.

        Returns ``(mrna_counts, nrna_counts)``, each ``dict[t_idx, dict[frag_len, count]]``, apart only because the
        writer names the two rows' reads differently.
        """
        if n_rna <= 0:
            return {}, {}

        N = len(self.transcripts)
        rng = self._rng
        frag_lengths = self._sample_rna_frag_lengths(n_rna)

        mrna_counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))
        nrna_counts: dict[int, dict[int, int]] = defaultdict(lambda: defaultdict(int))

        # The effective length is evaluated on the rows with abundance only; `eff` keeps zeros in the rest, so the
        # draw is bit-identical to computing every row.
        live = np.flatnonzero((self._mrna_abund > 0) | (self._nrna_abund > 0))
        live_keys = live.tolist()
        live_lengths = self._t_lengths[live]

        def weights_at(width: int) -> np.ndarray:
            eff = np.zeros(N, dtype=np.float64)
            if live.size:
                eff[live] = self.capture.partition_array("mrna", live_keys, live_lengths, width)
            return np.concatenate([self._mrna_abund * eff, self._nrna_abund * eff])

        widths, weights_per_width, counts_per_width = self._post_capture_length_allocation(
            frag_lengths, weights_at, rng, "RNA"
        )

        for width, weights, fc in zip(widths, weights_per_width, counts_per_width):
            fl, fc = int(width), int(fc)
            if fc <= 0:
                continue
            # a drawn width has weight: the allocation gives a weightless width no fragments
            probs = weights / weights.sum()
            indices = rng.choice(2 * N, size=fc, p=probs)
            unique_idx, idx_counts = np.unique(indices, return_counts=True)
            for idx, cnt in zip(unique_idx, idx_counts):
                idx, cnt = int(idx), int(cnt)
                if idx < N:
                    mrna_counts[idx][fl] += cnt
                else:
                    nrna_counts[idx - N][fl] += cnt

        return dict(mrna_counts), dict(nrna_counts)

    def _accumulate_gdna_counts(self, n_gdna: int) -> dict[tuple[int, int], int]:
        """Draw ``n_gdna`` fragments over lengths and genomic references: ``{(ref_idx, frag_len): count}``."""
        if n_gdna == 0:
            return {}
        if not self.genomic_refs:
            # Writing none of the requested fragments must raise: the truth would agree with itself and hide it.
            raise ValueError(
                f"{n_gdna} gDNA fragments requested but no genomic reference was declared "
                "(genomic_refs is empty)"
            )

        rng = self._rng
        frag_lengths = self._sample_gdna_frag_lengths(n_gdna)
        counts: dict[tuple[int, int], int] = {}
        # One batched call per fragment length covers every reference in one pass.
        ref_lengths = np.asarray(self._gdna_ref_lengths, dtype=np.int64)
        gdna_refs = list(self.genomic_refs)

        def weights_at(width: int) -> np.ndarray:
            return self.capture.partition_array("gdna", gdna_refs, ref_lengths, width)

        # Reallocate the length marginal (`_post_capture_length_allocation`); the per-reference draw follows.
        widths, weights_per_width, counts_per_width = self._post_capture_length_allocation(
            frag_lengths, weights_at, rng, "gDNA"
        )

        for width, chrom_eff, fc in zip(widths, weights_per_width, counts_per_width):
            fl, fc = int(width), int(fc)
            if fc <= 0:
                continue
            # a drawn width has weight: the allocation gives a weightless width no fragments
            chrom_probs = chrom_eff / chrom_eff.sum()
            chrom_indices = rng.choice(len(gdna_refs), size=fc, p=chrom_probs)
            unique_chroms, chrom_counts = np.unique(chrom_indices, return_counts=True)
            for ci, cc in zip(unique_chroms, chrom_counts):
                counts[(int(ci), fl)] = int(cc)
        return counts

    # -- Reads ------------------------------------------------------------------

    def _read_pairs(
        self,
        seq_bytes: np.ndarray,
        starts: np.ndarray,
        frag_len: int,
        read_len: int,
        swap: np.ndarray | None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Each fragment's two reads off its template, the mates exchanged where ``swap`` is set, then sequencing
        errors drawn into R1's batch and then R2's. Returns ``(r1_seqs, r2_seqs, r1_nm, r2_nm)``."""
        r2_seqs, r1_seqs = _batch_extract_reads(seq_bytes, starts, frag_len, read_len)
        if swap is not None:
            swapped = np.flatnonzero(swap)
            if len(swapped) > 0:
                r1_seqs[swapped], r2_seqs[swapped] = (
                    r2_seqs[swapped].copy(),
                    r1_seqs[swapped].copy(),
                )
        r1_nm = self._introduce_errors_batch(r1_seqs)
        r2_nm = self._introduce_errors_batch(r2_seqs)
        return r1_seqs, r2_seqs, r1_nm, r2_nm

    def _write_pair(
        self,
        bam_fh: pysam.AlignmentFile,
        qname: str,
        ref_id: int,
        tags: list,
        tlen: int,
        r1: _Mate,
        r2: _Mate,
    ) -> None:
        """One fragment's two oracle records, R1 then R2, each naming the other as its mate; ``tlen`` is signed from
        the leftmost mate's side."""
        r1_flag, r2_flag = _mate_flags(r1.reverse, r2.reverse)
        r1_tlen = tlen if r1.start <= r2.start else -tlen
        for mate, flag, other, mate_tlen in (
            (r1, r1_flag, r2, r1_tlen),
            (r2, r2_flag, r1, -r1_tlen),
        ):
            bam_fh.write(
                make_aligned_segment(
                    self._bam_header,
                    qname,
                    _bam_seq_from_fastq_bytes(mate.seq, mate.reverse),
                    flag,
                    ref_id,
                    mate.start,
                    mate.cigar,
                    ref_id,
                    other.start,
                    mate_tlen,
                    tags=[*tags, ("NM", int(mate.nm))],
                )
            )

    # -- RNA rows -------------------------------------------------------------

    def _rna_flips(self, count: int) -> np.ndarray | None:
        """Which RNA fragments are emitted R1-sense, ``None`` for none of them.

        The base emission is R1-antisense; a flip makes a fragment R1-sense. An R1-sense protocol flips the
        fragments the R1-antisense one keeps, the exact mirror on one RNG stream. The r1_sense=False path must keep
        consuming the RNG in this pattern, or every seeded library changes.
        """
        if self.strand_specificity < 1.0:
            flips = self._rng.random(count) >= self.strand_specificity
            return ~flips if self.r1_sense else flips
        if self.r1_sense:
            return np.ones(count, dtype=bool)
        return None

    def _write_rna_reads(
        self,
        t_idx: int,
        len_counts: dict[int, int],
        fastq: _FastqPair,
        bam_fh: pysam.AlignmentFile | None,
        *,
        is_nrna: bool,
    ) -> int:
        """Write all reads of one RNA row: a transcript's row, or an entity's nascent row (``is_nrna``), which differs
        only by its read-name tag. The template is the row's own sequence; an entity's is its genomic span.
        """
        t = self.transcripts[t_idx]
        if is_nrna and len(t.exons) != 1:
            raise ValueError(f"nascent row on a multi-exon transcript {t.t_id}")
        seq_bytes = self._t_seq_bytes[t_idx]
        seq_len = int(self._t_lengths[t_idx])
        strand_char = self._t_strand_chars[t_idx]
        ref_id = self._ref_name_to_id.get(t.ref) if bam_fh is not None else None
        n_written = 0

        for frag_len in sorted(len_counts):
            count = len_counts[frag_len]
            if seq_len - frag_len + 1 <= 0:
                continue
            read_len = min(self.sim_params.read_length, frag_len)
            starts = self.capture.sample_starts("mrna", t_idx, seq_len, frag_len, count, self._rng)
            flips = self._rna_flips(count)
            reads = self._read_pairs(seq_bytes, starts, frag_len, read_len, flips)
            qnames = [
                format_rna_name(
                    t.t_id, start, start + frag_len, strand_char, n_written + i, nascent=is_nrna
                )
                for i, start in enumerate(starts.tolist())
            ]
            fastq.write(qnames, reads[0], reads[1], self._get_quals(read_len))
            if bam_fh is not None and ref_id is not None:
                self._write_rna_bam(
                    bam_fh, qnames, reads, t, starts, frag_len, read_len, flips, ref_id
                )
            n_written += count
        return n_written

    def _write_rna_bam(
        self,
        bam_fh: pysam.AlignmentFile,
        qnames: list[str],
        reads: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        t: Transcript,
        starts: np.ndarray,
        frag_len: int,
        read_len: int,
        flips: np.ndarray | None,
        ref_id: int,
    ) -> None:
        """The oracle records of one batch of RNA fragments, each read placed by the template's exon structure."""
        r1_seqs, r2_seqs, r1_nm, r2_nm = reads
        for i in range(len(qnames)):
            frag_start = int(starts[i])
            frag_end = frag_start + frag_len
            r2_blocks = transcript_to_genomic_blocks(
                frag_start, min(frag_start + read_len, frag_end), t
            )
            r1_blocks = transcript_to_genomic_blocks(
                max(frag_end - read_len, frag_start), frag_end, t
            )
            if not r2_blocks or not r1_blocks:
                continue

            r2_is_rev, r1_is_rev = (False, True) if t.strand == Strand.POS else (True, False)
            if flips is not None and flips[i]:
                r2_is_rev, r1_is_rev = not r2_is_rev, not r1_is_rev
                r1_blocks, r2_blocks = r2_blocks, r1_blocks

            tlen = max(r1_blocks[-1][1], r2_blocks[-1][1]) - min(r1_blocks[0][0], r2_blocks[0][0])
            tags = self._nh1_tags
            if len(r2_blocks) > 1 or len(r1_blocks) > 1:
                xs_strand = "-" if t.strand == Strand.NEG else "+"
                tags = [("NH", 1), ("XS", xs_strand, "A")]  # type 'A' (char), per SAM convention
            self._write_pair(
                bam_fh,
                qnames[i],
                ref_id,
                tags,
                tlen,
                _Mate(r1_seqs[i], r1_nm[i], r1_blocks[0][0], blocks_to_cigar(r1_blocks), r1_is_rev),
                _Mate(r2_seqs[i], r2_nm[i], r2_blocks[0][0], blocks_to_cigar(r2_blocks), r2_is_rev),
            )

    # -- gDNA -----------------------------------------------------------------

    def _gdna_reverse(self, ref: str, starts: np.ndarray, count: int) -> np.ndarray:
        """Which gDNA fragments come from the reverse strand: 50/50, or with overdispersion at their region's shared
        +strand rate."""
        strand_regions = self._gdna_strand_regions.get(ref)
        if strand_regions is None:
            return self._rng.integers(0, 2, size=count).astype(bool)
        boundaries, p_plus = strand_regions
        reg_idx = np.clip(np.searchsorted(boundaries, starts, side="right") - 1, 0, len(p_plus) - 1)
        return self._rng.random(count) >= p_plus[reg_idx]

    def _write_gdna_chunk(
        self,
        ref_idx: int,
        fl: int,
        count: int,
        fastq: _FastqPair,
        bam_fh: pysam.AlignmentFile | None,
        n_offset: int,
    ) -> int:
        """Write ``count`` gDNA fragments of length ``fl`` on reference ``ref_idx``, numbering their reads on from
        ``n_offset``; return how many were written."""
        if count <= 0:
            return 0
        ref = self.genomic_refs[ref_idx]
        chrom_len = self._gdna_ref_lengths[ref_idx]
        if chrom_len - fl + 1 <= 0:
            return 0
        read_len = min(self.sim_params.read_length, fl)
        starts = self.capture.sample_starts("gdna", ref, chrom_len, fl, count, self._rng)
        reverse = self._gdna_reverse(ref, starts, count)
        # a reverse-strand fragment is the reverse complement, so its mates swap
        reads = self._read_pairs(self._genome.get(ref), starts, fl, read_len, reverse)
        qnames = [
            format_gdna_name(
                ref, start, start + fl, STRAND_REVERSE if rev else STRAND_FORWARD, n_offset + j
            )
            for j, (start, rev) in enumerate(zip(starts.tolist(), reverse.tolist()))
        ]
        fastq.write(qnames, reads[0], reads[1], self._get_quals(read_len))
        ref_id = self._ref_name_to_id.get(ref)
        if bam_fh is not None and ref_id is not None:
            self._write_gdna_bam(bam_fh, qnames, reads, starts, fl, read_len, reverse, ref_id)
        return count

    def _write_gdna_bam(
        self,
        bam_fh: pysam.AlignmentFile,
        qnames: list[str],
        reads: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        starts: np.ndarray,
        fl: int,
        read_len: int,
        reverse: np.ndarray,
        ref_id: int,
    ) -> None:
        """The oracle records of one batch of gDNA fragments: each read one unspliced block at a fragment end."""
        r1_seqs, r2_seqs, r1_nm, r2_nm = reads
        cigar = [(pysam.CMATCH, read_len)]
        for j in range(len(qnames)):
            start = int(starts[j])
            end = start + fl
            if reverse[j]:
                r1 = _Mate(r1_seqs[j], r1_nm[j], start, cigar, False)
                r2 = _Mate(r2_seqs[j], r2_nm[j], end - read_len, cigar, True)
            else:
                r1 = _Mate(r1_seqs[j], r1_nm[j], end - read_len, cigar, True)
                r2 = _Mate(r2_seqs[j], r2_nm[j], start, cigar, False)
            self._write_pair(bam_fh, qnames[j], ref_id, self._nh1_tags, end - start, r1, r2)

    # -- Abundance-weighted pool splits (single-condition / Scenario use) ----

    def _rna_weight(self) -> float:
        """Σ over every RNA row of abundance × capture-aware effective length at the mean RNA fragment length."""
        mean_frag = int(self.sim_params.frag_mean)
        n = len(self.transcripts)
        eff = self.capture.partition_array("mrna", range(n), self._t_lengths, mean_frag)
        return float(np.sum((self._mrna_abund + self._nrna_abund) * eff))

    def rna_gdna_split(self, n_total: int, gdna_abundance: float) -> tuple[int, int]:
        """Split ``n_total`` fragments into ``(n_rna, n_gdna)`` by abundance × effective length, gDNA weighing
        ``gdna_abundance`` × the genomic references' effective length. The nascent share of ``n_rna`` is not
        imposed."""
        rna_w = self._rna_weight()
        gdna_mean_frag = int(self.gdna_config.frag_mean)
        genome_eff = sum(
            self.capture.partition("gdna", ref, length, gdna_mean_frag)
            for ref, length in zip(self.genomic_refs, self._gdna_ref_lengths)
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
        """Draw *n_rna* RNA fragments (:meth:`_accumulate_rna_counts`) and *n_gdna* gDNA fragments, then write them.
        With *n_workers* > 1 the writing is sharded across forked workers and the shards concatenated. Returns
        ``(r1_path, r2_path, bam_path or None)``.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        r1_path = output_dir / f"{prefix}_R1.fq.gz"
        r2_path = output_dir / f"{prefix}_R2.fq.gz"
        bam_path = output_dir / f"{prefix}_oracle.bam" if oracle_bam else None

        # Counts are drawn in the parent from self._rng, so they follow the configured seed.
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
            self._write_files(
                r1_path,
                r2_path,
                bam_path,
                [(t_idx, mrna_counts[t_idx]) for t_idx in sorted(mrna_counts)],
                [(t_idx, nrna_counts[t_idx]) for t_idx in sorted(nrna_counts)],
                list(gdna_counts.items()),
                threads=CONSTANTS.resources.fastq_gzip_threads,
            )
        else:
            self._write_sharded(
                output_dir,
                prefix,
                r1_path,
                r2_path,
                bam_path,
                mrna_counts,
                nrna_counts,
                gdna_counts,
                n_workers,
            )

        logger.info(
            "Wrote %d read pairs -> %s (RNA=%d, gDNA=%d, oracle=%s, workers=%d)",
            total_mrna + total_nrna + total_gdna,
            output_dir,
            n_rna,
            n_gdna,
            oracle_bam,
            n_workers,
        )
        return r1_path, r2_path, bam_path

    # -- Serial / parallel write paths --------------------------------------

    def _write_files(
        self,
        r1_path: Path,
        r2_path: Path,
        bam_path: Path | None,
        mrna_items: list[tuple[int, dict[int, int]]],
        nrna_items: list[tuple[int, dict[int, int]]],
        gdna_items: list[tuple[tuple[int, int], int]],
        *,
        threads: int,
    ) -> None:
        """Write drawn fragments, in the order given, to one FASTQ pair and, with *bam_path*, one oracle BAM: the
        transcript rows, then the nascent rows, then the gDNA chunks, whose read numbers run on across chunks."""
        with (
            _open_gzip_write(r1_path, threads) as r1_fh,
            _open_gzip_write(r2_path, threads) as r2_fh,
        ):
            fastq = _FastqPair(r1_fh, r2_fh)
            bam_fh = None
            if bam_path is not None:
                bam_fh = pysam.AlignmentFile(str(bam_path), "wb", header=self._bam_header)
            try:
                for t_idx, len_counts in mrna_items:
                    self._write_rna_reads(t_idx, len_counts, fastq, bam_fh, is_nrna=False)
                for t_idx, len_counts in nrna_items:
                    self._write_rna_reads(t_idx, len_counts, fastq, bam_fh, is_nrna=True)
                n_offset = 0
                for (ref_idx, fl), count in gdna_items:
                    n_offset += self._write_gdna_chunk(ref_idx, fl, count, fastq, bam_fh, n_offset)
                fastq.flush()
            finally:
                if bam_fh is not None:
                    bam_fh.close()

    def _write_sharded(
        self,
        output_dir: Path,
        prefix: str,
        r1_path: Path,
        r2_path: Path,
        bam_path: Path | None,
        mrna_counts: dict[int, dict[int, int]],
        nrna_counts: dict[int, dict[int, int]],
        gdna_counts: dict[tuple[int, int], int],
        n_workers: int,
    ) -> None:
        """Write in *n_workers* forked shards, the rows and chunks spread by fragment count (greedy LPT) and each
        shard seeded from this simulator's stream, then concatenate the shards in order."""

        def by_row(item):
            return sum(item[1].values())

        mrna_shards = _shard_by_count(
            [(t, dict(d)) for t, d in mrna_counts.items()], n_workers, weight=by_row
        )
        nrna_shards = _shard_by_count(
            [(t, dict(d)) for t, d in nrna_counts.items()], n_workers, weight=by_row
        )
        gdna_shards = _shard_by_count(
            list(gdna_counts.items()), n_workers, weight=lambda item: item[1]
        )

        shard_dir = output_dir / f".{prefix}_shards"
        shard_dir.mkdir(parents=True, exist_ok=True)
        base_seed = int(self._rng.integers(0, 2**31 - 1))
        oracle_bam = bam_path is not None
        tasks = [
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
            for k in range(n_workers)
        ]

        # Forked children inherit self through a module global: it cannot be pickled (pysam handles).
        global _WORKER_SIM
        _WORKER_SIM = self
        ctx = multiprocessing.get_context("fork")
        try:
            with ctx.Pool(n_workers) as pool:
                results: list[tuple[str, str, str | None]] = pool.map(_shard_task, tasks)
        finally:
            _WORKER_SIM = None

        _concat_shards(results, r1_path, r2_path, bam_path)
        with contextlib.suppress(OSError):
            shard_dir.rmdir()

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
        self._rng = np.random.default_rng(seed)
        shard = Path(shard_dir) / f"{prefix}.shard{shard_id:03d}"
        r1_path = Path(f"{shard}.R1.fq.gz")
        r2_path = Path(f"{shard}.R2.fq.gz")
        bam_path = Path(f"{shard}.bam") if oracle_bam else None
        # one compressor thread per file: the shards are the parallelism, each a process of its own
        self._write_files(r1_path, r2_path, bam_path, mrna_items, nrna_items, gdna_items, threads=1)
        return (str(r1_path), str(r2_path), str(bam_path) if bam_path else None)

    def close(self) -> None:
        """Release the sequence cache and the FASTA handle."""
        self._genome.clear()
        self.fasta.close()


# ═══════════════════════════════════════════════════════════════════
# Parallel-shard helpers (module-level for fork-based multiprocessing)
# ═══════════════════════════════════════════════════════════════════


# Set in the parent before the Pool forks, so children share the sequence cache copy-on-write.
_WORKER_SIM: "WholeGenomeSimulator | None" = None


def _shard_task(args):
    """Pool worker entry. Dispatches to the inherited simulator's _run_shard."""
    if _WORKER_SIM is None:
        raise RuntimeError("Worker simulator not initialized (fork inheritance failed)")
    return _WORKER_SIM._run_shard(*args)


def _shard_by_count(items, n_shards, *, weight):
    """Greedy LPT bin-packing: distribute items across shards by descending weight."""
    sorted_items = sorted(items, key=lambda x: -weight(x))
    shards: list[list] = [[] for _ in range(n_shards)]
    loads = [0] * n_shards
    for item in sorted_items:
        idx = min(range(n_shards), key=lambda i: loads[i])
        shards[idx].append(item)
        loads[idx] += weight(item)
    return shards


def _concat_shards(
    results: list[tuple[str, str, str | None]], r1_path: Path, r2_path: Path, bam_path: Path | None
) -> None:
    """Concatenate the shards' files, in shard order, into the run's, then delete them."""
    shard_r1 = [Path(r[0]) for r in results]
    shard_r2 = [Path(r[1]) for r in results]
    shard_bams = [Path(r[2]) for r in results if r[2] is not None]

    logger.info("Concatenating %d FASTQ shards...", len(shard_r1))
    _concat_files_binary(shard_r1, r1_path)
    _concat_files_binary(shard_r2, r2_path)
    if bam_path is not None and shard_bams:
        logger.info("Concatenating %d BAM shards...", len(shard_bams))
        pysam.cat("-o", str(bam_path), *[str(p) for p in shard_bams])

    for p in shard_r1 + shard_r2 + shard_bams:
        with contextlib.suppress(OSError):
            p.unlink()


def _concat_files_binary(srcs: list[Path], dst: Path) -> None:
    """Concatenate files byte-wise. Valid for gzip streams (per RFC 1952)."""
    with open(dst, "wb") as out:
        for s in srcs:
            with open(s, "rb") as f:
                shutil.copyfileobj(f, out, length=CONSTANTS.resources.file_copy_bytes)
