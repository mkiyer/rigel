"""Every setting and every fixed number Rigel's code reads, in one place.

Two trees of frozen dataclasses, validated on construction:

* :class:`PipelineConfig` — what a run may set (the CLI, a YAML file, an instrument's ``--set``): ``em``, ``scan``,
  ``scoring`` and ``calibration``. :class:`IndexConfig` is the same for ``rigel index``.
* :data:`CONSTANTS` — what the code fixes: numerical methods, quality-control thresholds and resource budgets, by
  component. Each value is documented with why it has that value; the code reads it by name and never restates it.

Not here: identifiers and formats (enum codes, bit flags, schema versions, file names), which live with the format
they define; mathematical definitions (gDNA's strand probability ½, Jeffreys' ½), named where they are used; and the
constants the native kernels share with Python, defined in C++ and exported by :mod:`rigel.native`.
"""

from __future__ import annotations

import math
import operator
import os
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    import numpy as np


def _check_integer(where: str, value, least: int) -> None:
    """Refuse on construction an integer setting the run would otherwise refuse only after the scan, or read as another
    value: one that is not an integer (YAML reads ``1e3`` as a string), or one below ``least``."""
    try:
        operator.index(value)
    except TypeError:
        raise TypeError(f"{where} must be an integer; got {value!r}.") from None
    if value < least:
        raise ValueError(f"{where} must be >= {least}; got {value}.")


# ======================================================================
# EM algorithm configuration
# ======================================================================


@dataclass(frozen=True)
class EMConfig:
    """Configuration for the EM algorithm and posterior assignment.

    Parameters
    ----------
    seed : int
        Seed of the ``sample`` assignment's draw (default 0); another seed is another draw from the same
        posterior. Fixed, never a clock, so the draw repeats from run to run.
    mode : {"vbem", "map"}
        Algorithm variant (default ``"vbem"``).
    iterations : int
        Maximum EM iterations (default 1000).
    convergence_delta : float
        Convergence threshold for theta updates (default 1e-6).
    assignment_mode : {"sample", "fractional"}
        Post-EM fragment assignment (default ``"sample"``). ``"sample"`` gives every fragment to one
        component, count first: each component's fractional count is rounded within its EM locus, then
        each fragment is drawn from its own posterior, re-weighted toward the components still short of
        their rounded count, so whole counts land on the rounded fractional counts. ``"fractional"``
        keeps each fragment's posterior weights.
    """

    seed: int = 0
    mode: Literal["vbem", "map"] = "vbem"
    iterations: int = 1000
    convergence_delta: float = 1e-6
    assignment_mode: str = "sample"
    n_threads: int = 0
    """Number of threads for parallel locus EM.

    ``0`` (default) → use all available cores.
    ``1`` → sequential.
    Any positive value → cap at that many threads (the EM's own thread pool).
    """

    def __post_init__(self):
        if self.mode not in ("map", "vbem"):
            raise ValueError(f"Unknown EM mode: {self.mode!r}")
        if self.assignment_mode not in ("fractional", "sample"):
            raise ValueError(f"Unknown assignment mode: {self.assignment_mode!r}")
        # `None` would seed the draw from OS entropy, and one input would stop giving one answer; the draw's generator
        # refuses a negative seed
        _check_integer("EMConfig.seed", self.seed, 0)
        # the EM always takes one accelerated step, so a budget below one is never what runs
        _check_integer("EMConfig.iterations", self.iterations, 1)
        # 0 is every core; a negative count would be read as 0
        _check_integer("EMConfig.n_threads", self.n_threads, 0)
        # 0 runs the whole budget, as no change falls below it; NaN or a negative threshold would too, silently, and
        # infinity would stop at the first step
        delta = self.convergence_delta
        if (
            isinstance(delta, bool)
            or not isinstance(delta, (int, float))
            or not 0.0 <= delta < math.inf
        ):
            raise ValueError(
                f"EMConfig.convergence_delta must be a finite number >= 0; got {delta!r}."
            )


# ======================================================================
# Fragment scoring configuration
# ======================================================================


@dataclass(frozen=True)
class FragmentScoringConfig:
    """Configuration for fragment scoring penalties.

    A hit's sequence is priced by the library itself, never by a setting: its mismatches and the read bases it clips,
    each against an aligned base that matches, at the mismatch rate the scan reads on the unique unspliced hits
    (``rigel.scoring.sequence_log_penalties``).

    Parameters
    ----------
    overhang_alpha : float
        The probability kept per base of overhang, in [0, 1]: the default 0.1 cuts each base's likelihood 10×, 1 is
        no penalty, and 0 is a hard gate that excludes any candidate with overhang. Its log is the scorer's
        per-base log-penalty (:attr:`overhang_log_penalty`, −inf at the gate).
    pruning_min_posterior : float
        A candidate whose likelihood falls below this fraction of the best at its place is pruned, in [0, 1): at 1 a
        place would keep only its tied-best candidates. A multimapper's places are pruned each on its own, so every
        place stays in the EM.
    """

    overhang_alpha: float = 0.1
    pruning_min_posterior: float = 1e-4

    def __post_init__(self):
        # NaN fails both comparisons; a bool is an int to Python but never an alpha or a threshold
        value = self.overhang_alpha
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not 0.0 <= value <= 1.0
        ):
            raise ValueError(f"--overhang-alpha must be a number in [0, 1], got {value!r}")
        value = self.pruning_min_posterior
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0.0 <= value < 1.0:
            raise ValueError(f"--pruning-min-posterior must be a number in [0, 1), got {value!r}")

    @property
    def overhang_log_penalty(self) -> float:
        """``log(overhang_alpha)``: the per-base log-penalty, −inf at the hard gate."""
        return alpha_log_penalty(self.overhang_alpha)


def alpha_log_penalty(alpha: float) -> float:
    """An alpha's log-penalty: ``log(alpha)``, or −inf at 0, the hard gate."""
    return math.log(alpha) if alpha > 0.0 else -math.inf


# ======================================================================
# BAM scanning and buffering configuration
# ======================================================================


@dataclass(frozen=True)
class BamScanConfig:
    """Configuration for the BAM scanning and buffering stage.

    Parameters
    ----------
    skip_duplicates : bool
        Discard reads marked as duplicates (default True).
    include_multimap : bool
        Include read pairs with more than one usable alignment (default True); off, each is excluded and counted
        (``summary.json``'s ``read_pairs.excluded_by_option``).
    max_frag_length : int
        Maximum fragment length for histogram models (default 1000).
    sj_strand_tag : str or tuple of str
        BAM tag(s) for splice-junction strand (default ``"auto"``). ``"auto"`` reads XS and minimap2's ts on every
        spliced record: agreeing tags give the strand, disagreeing ones an ambiguous motif. A named tag, or two in
        order, reads only those, the first present winning; ``"none"`` reads no tag.
    log_every : int
        Log scoring progress (debug level) every N buffered fragments,
        checked once per buffer chunk (default 1M).
    total_threads : int
        Total thread budget available to the scan stage (default 0 = all cores).
        Rigel reserves decompression threads from this budget and uses the
        remainder for scan workers.
    bgzf_threads : int or None
        BGZF decompression threads within the scan thread budget. ``None``, the
        default, derives them from the budget (see
        :meth:`BamScanConfig.resolved_scan_threads`); an integer overrides.
    fragments_per_chunk : int
        Buffered fragments per chunk (default 1M).
    read_name_batch_size : int
        Read-name groups per native scanner input queue item (default 512).
    buffer_size_bytes : int
        Memory for the buffered fragments before the oldest chunks spill to disk (default 2 GiB). It bounds the
        buffer only: beside it each scan worker fills one chunk and the scanner's queue holds
        ``OUTPUT_CHUNKS_PER_WORKER`` finished ones per worker (the native ``constants.h``), and the spill writer up
        to ``CONSTANTS.resources.spill_queue_chunks``.
    spill_dir : Path, str, or None
        Directory for spilled buffer chunks (default None).
    """

    skip_duplicates: bool = True
    include_multimap: bool = True
    max_frag_length: int = 1000
    sj_strand_tag: str | tuple[str, ...] = "auto"
    log_every: int = 1_000_000
    total_threads: int = 0
    bgzf_threads: int | None = None
    fragments_per_chunk: int = 1_000_000
    read_name_batch_size: int = 512
    buffer_size_bytes: int = 2 * 1024**3
    spill_dir: Path | str | None = None
    """Scan buffer spill directory (default ``None`` = system temp dir)."""

    intron_overhang_tolerance: int = 3
    """Intron-overhang tolerance ``K`` (bp): how many bases an aligned block's end may lie
    inside an annotated intron, as when an aligner runs a read a base or two across a splice
    site. It acts at the edge of a gap the reads did not sequence (an implicit splice) and at
    a fragment's outer end.

    The pair rule applies it along a witness transcript, when two mates are longer than
    ``max_frag_length`` along the genome: the witness's intron in an unaligned stretch
    between them is cut whole when it lies within the stretch widened by ``K``, and the
    witness offers no path when one of its introns holds more than ``K`` aligned bases. The
    resolver applies the same two tests when it enumerates a fragment's gap hypotheses:
    a candidate's intron in an unsequenced gap is cut whole within the gap widened by ``K``,
    and a candidate one of whose introns holds more than ``K`` aligned bases joins no
    hypothesis. So ``K`` decides which pairs are one molecule, the candidates' lengths and
    calibration's tally, and the scan-cache key includes it. It never matches a junction:
    junctions are matched exactly.
    """

    def __post_init__(self) -> None:
        if self.total_threads < 0:
            raise ValueError(f"BamScanConfig.total_threads must be >= 0; got {self.total_threads}.")
        if self.bgzf_threads is not None and self.bgzf_threads < 0:
            raise ValueError(
                f"BamScanConfig.bgzf_threads must be >= 0 or None; got {self.bgzf_threads}."
            )
        if self.fragments_per_chunk < 1:
            raise ValueError(
                f"BamScanConfig.fragments_per_chunk must be >= 1; got {self.fragments_per_chunk}."
            )
        if self.read_name_batch_size < 1:
            raise ValueError(
                f"BamScanConfig.read_name_batch_size must be >= 1; got {self.read_name_batch_size}."
            )
        if self.buffer_size_bytes < 0:
            raise ValueError(
                f"BamScanConfig.buffer_size_bytes must be >= 0; got {self.buffer_size_bytes}."
            )
        if self.intron_overhang_tolerance < 0:
            raise ValueError(
                f"BamScanConfig.intron_overhang_tolerance must be >= 0; "
                f"got {self.intron_overhang_tolerance}."
            )
        # the router's progress log divides by it
        _check_integer("BamScanConfig.log_every", self.log_every, 1)
        # the scanner's tag modes; anything else is refused here rather than reinterpreted by the scanner
        tags = self.sj_strand_tag
        spec = ",".join(tags) if isinstance(tags, (list, tuple)) else tags
        if spec not in ("auto", "none", "", "XS", "ts", "XS,ts", "ts,XS"):
            raise ValueError(
                "BamScanConfig.sj_strand_tag must be auto, none, XS, ts, or XS and ts in either order; "
                f"got {tags!r}."
            )

    def resolved_total_threads(self) -> int:
        """Return the concrete total scan thread budget."""
        if self.total_threads == 0:
            return os.cpu_count() or 1
        return self.total_threads

    def resolved_scan_threads(self) -> tuple[int, int]:
        """Return ``(scan_worker_threads, bgzf_threads)`` within the budget.

        The split is a ratio, not a count: one decompression thread keeps about
        ``CONSTANTS.resources.scan_workers_per_bgzf_thread`` scan workers fed on a name-sorted BAM, so a
        decompression thread beyond that ratio is a worker taken away. ``bgzf_threads``
        (``--scan-bgzf-threads``) overrides it.
        """
        total = self.resolved_total_threads()
        want = (
            total // CONSTANTS.resources.scan_workers_per_bgzf_thread
            if self.bgzf_threads is None
            else self.bgzf_threads
        )
        bgzf = min(want, max(total - 1, 0))
        scan_workers = max(1, total - bgzf)
        return scan_workers, bgzf


# ======================================================================
# Calibration configuration
# ======================================================================


@dataclass(frozen=True)
class CalibrationConfig:
    """Configuration for the calibrator (:func:`rigel.calibration.calibrate.calibrate`).

    The calibrator is a belief-propagation sweep over the region-boundary chain, one forward and one backward pass
    per solve; each region's gDNA share is solved on a log-odds lattice (``sweep_logodds_step``,
    ``sweep_logodds_window``).
    """

    #: The λ lattice's step, in nats of log-odds: ψ's grid is ``K = round(2L/step) + 1`` points on ``[−L, L]``, so
    #: widening the bracket never coarsens it. One lattice serves every consumer (the read-out, the message rows, the
    #: AMBIG cube's λ axis, the composition prior), so nothing is regridded. At 0.2 the ½-quantile read-out puts a
    #: slot's composition within 1.25 % of its mass (quantised by at most ``n·f(1−f)·step/4`` fragments); the AMBIG
    #: cube costs ``K × (K_t + 2)`` per slot, and 0.2 is the coarsest step that loses no accuracy.
    sweep_logodds_step: float = 0.2

    #: Log-odds grid floor ``L``: ``f_g ∈ [σ(−L), σ(L)]``, the range the Beta(½,½) reference needs to stay proper. A
    #: pass solves on ``max(L, the fitted prior's demand)`` (``DensityLandscape.required_logodds_window``), because a
    #: bracket narrower than the landscape's log-dynamic range would leave ψ no coordinate for where its prior points.
    sweep_logodds_window: float = 10.0

    #: The sweep's working set: the chain, cut at every intergenic region, is solved in locus blocks of at most this
    #: many slots (``calibration.region_chain.locus_blocks``; ``None`` is one block). Memory only: the answer is the
    #: same at every value, the arena is ``16 · slots · K · 8 B · threads``, the wall time is flat in it, and 1000 is
    #: the smallest size whose memory saving still exceeds run-to-run noise.
    sweep_block_slots: int | None = 1000

    #: The message policy: what one neighbour tells another. ``"transfer"`` (the composition transfer,
    #: :class:`~rigel.calibration.messages.transfer.TransferPolicy`) is the default; ``"silent"``
    #: (:class:`~rigel.calibration.messages.silent.SilentPolicy`) sends nothing and is the floor every policy is
    #: judged against. Messages serve the slots whose own solve has no composition channel: unstranded data and the
    #: both-stranded (AMBIG) slots. An unknown name raises. A change of default must be checked with the instruments,
    #: not only the suite, because the tests install their policy themselves.
    message_policy: str = "transfer"

    #: Calibration refit iterations, the prior bootstrap: each re-fits the population gDNA landscape
    #: (:class:`~rigel.calibration.landscape.DensityLandscape`) on the current solve, resets the belief and re-solves,
    #: so only the landscape carries between iterations. ``0`` runs the prior-free pass alone; each iteration costs one
    #: landscape fit and one sweep. Why 3: the third refit reaches nearly all of the bootstrap's gain, each later one
    #: adding 2–3× less (a reading not repeated on the current landscape).
    calib_refit_iters: int = 3

    #: ψ's thread budget (``0``, the default, is every core); ``--threads`` sets it beside the scan's and the EM's. The
    #: solve is bit-identical at every count.
    n_threads: int = 0

    #: The intron abundance model's departure prior ε: the probability that a region's gDNA is not at the intergenic
    #: background (a probe's reach, a copy-number change, a mappability gap, unannotated transcription), under which any
    #: split of its count is as likely as any other (``gdna_density.BackgroundLaw``). It caps how far the background
    #: can argue against an intron's composition, ``log(1 + ((1 − ε)/ε)·(C + 1)·NB_max)`` nats for ``C`` fragments, so
    #: an intron's own decisive strand evidence outweighs a background it contradicts; and it keeps a region the law
    #: cannot explain from setting the law. 0.01 is odds of 99 to 1 that a region is at the background. The intergenic
    #: pool is always fitted at this value. A stranded library measures its introns' own ε from their strands
    #: (``calibration.intron_departure``), and this value is that estimate's prior mode; an unstranded library's introns
    #: read at this value, where calibration against truth sits within 1 % of its best for every in-scope stratum from
    #: 1e-4 to 0.05. Raised toward 1 it flattens the abundance row, which turns the model off where the intergenic
    #: regions are not a gDNA reference.
    intron_departure_prior: float = 0.01

    def __post_init__(self) -> None:
        if self.n_threads < 0:
            raise ValueError(f"CalibrationConfig.n_threads must be >= 0; got {self.n_threads}.")
        if not (0.0 < float(self.intron_departure_prior) < 1.0):
            raise ValueError(
                "CalibrationConfig.intron_departure_prior must be a probability strictly between 0 and 1; "
                f"got {self.intron_departure_prior}."
            )
        if self.calib_refit_iters < 0:
            raise ValueError(
                f"CalibrationConfig.calib_refit_iters must be >= 0; got {self.calib_refit_iters}."
            )
        if self.sweep_block_slots is not None and self.sweep_block_slots < 1:
            raise ValueError(
                f"CalibrationConfig.sweep_block_slots must be >= 1 or None; got {self.sweep_block_slots}."
            )
        if not (float(self.sweep_logodds_window) > 0.0):
            raise ValueError(
                "CalibrationConfig.sweep_logodds_window (L) must be > 0; "
                f"got {self.sweep_logodds_window}."
            )
        if not (0.0 < float(self.sweep_logodds_step) <= 2.0 * float(self.sweep_logodds_window)):
            raise ValueError(
                "CalibrationConfig.sweep_logodds_step must be > 0 and no wider than the window "
                f"(2L = {2.0 * float(self.sweep_logodds_window)}); got {self.sweep_logodds_step}."
            )


# ======================================================================
# Index configuration (``rigel index``)
# ======================================================================


@dataclass(frozen=True)
class IndexConfig:
    """What ``rigel index`` may set; the CLI's defaults and :meth:`rigel.index.TranscriptIndex.build`'s are these, and
    ``build`` validates its arguments as one.

    The index is built once and read by every run, so these are fixed for every quantification against it.
    """

    #: Compression of the index's Feather tables: ``"lz4"``, ``"zstd"`` or ``"uncompressed"``.
    feather_compression: Literal["lz4", "zstd", "uncompressed"] = "lz4"

    #: Write a human-readable TSV beside every Feather table.
    write_tsv: bool = True

    #: ``"strict"`` fails on a malformed GTF record; ``"warn-skip"`` logs it and skips it.
    gtf_parse_mode: Literal["strict", "warn-skip"] = "strict"

    #: Transcripts with identical exon coordinates are unidentifiable in quantification. ``False`` fails with a
    #: report of each group; ``True`` keeps the lexicographically smallest transcript ID per group, loss-free.
    collapse_duplicate_transcripts: bool = False

    #: Multi-exon transcripts' start and end sites within this many bp cluster into one synthetic span,
    #: which models their nascent RNA (``index.create_nrna_transcripts``).
    nascent_tolerance: int = 20

    def __post_init__(self) -> None:
        if self.feather_compression not in ("lz4", "zstd", "uncompressed"):
            raise ValueError(f"Unknown Feather compression: {self.feather_compression!r}")
        if self.gtf_parse_mode not in ("strict", "warn-skip"):
            raise ValueError(f"Unknown GTF parse mode: {self.gtf_parse_mode!r}")
        if self.nascent_tolerance < 0:
            raise ValueError(
                f"IndexConfig.nascent_tolerance must be >= 0; got {self.nascent_tolerance}."
            )


# ======================================================================
# Top-level pipeline configuration
# ======================================================================


@dataclass(frozen=True)
class PipelineConfig:
    """Top-level pipeline configuration composing all sub-configs.

    Pass a single ``PipelineConfig`` to ``run_pipeline`` instead of one
    keyword argument per setting.
    """

    em: EMConfig = field(default_factory=EMConfig)
    scan: BamScanConfig = field(default_factory=BamScanConfig)
    scoring: FragmentScoringConfig = field(default_factory=FragmentScoringConfig)
    calibration: CalibrationConfig = field(default_factory=CalibrationConfig)
    #: Seed of the second pass's draw of one (placement, hypothesis) per held fragment. Separate from ``em.seed``, so
    #: a change of EM seed never re-draws the held fragments and moves the tally the EM reads.
    second_pass_seed: int = 0
    annotated_bam_path: str | Path | None = None
    emit_locus_stats: bool = False

    def __post_init__(self) -> None:
        # the second pass's generator refuses a negative seed
        _check_integer("PipelineConfig.second_pass_seed", self.second_pass_seed, 0)

    def to_dict(self) -> dict:
        """JSON-serializable dict of all configuration fields."""
        from dataclasses import asdict

        d = asdict(self)
        # Convert Path objects to strings
        if d.get("annotated_bam_path") is not None:
            d["annotated_bam_path"] = str(d["annotated_bam_path"])
        if d.get("scan", {}).get("spill_dir") is not None:
            d["scan"]["spill_dir"] = str(d["scan"]["spill_dir"])
        return d


# ======================================================================
# CONSTANTS — what the code fixes (not user-configurable)
# ======================================================================
#
# Every number an algorithm, a quality check or a resource budget depends on, each documented once with why it has
# its value. The code reads them by name (``CONSTANTS.landscape.grid_points``) and never restates one. They are not
# run settings: a change to any of them is a change to the method. An experiment varies one on a copy
# (``dataclasses.replace``), never by editing a module.


@dataclass(frozen=True)
class JunctionFitConstants:
    """How the genuine-junction strand fit finds its maximum (:func:`rigel.junction_fit.genuine_sense_fraction`).

    The fit is an exact maximum-likelihood problem: the genuine junctions' wrong-strand rate κ, and three class
    shares (genuine, splice artifact, reversed), over the per-junction strand table. Nothing here is a model
    choice: every field sets how the maximum is found, and at these values the fit reaches the maximum an
    independent search finds.
    """

    # --- The search over κ: a scan locates the maximum, Brent's method refines it ---

    #: κ is searched on ``(floor, ½ − floor)``. Within the floor of 0 or of ½, κ × (the RNA reads) moves by less
    #: than one read for any library under 10¹² spliced reads, so the answer cannot change.
    kappa_floor: float = 1e-12

    #: The scan's points, uniform in ``u = log(κ/(½ − κ))`` across that range, so the scan is as fine near ½ as
    #: near 0: at 129 points u steps by 0.42 (κ/(½ − κ) by ×1.5). A weakly stranded library's maximum rises out
    #: of a flat all-artifact plateau just below ½, and a coarser scan can step over it. Each point costs one share
    #: solve.
    kappa_scan_points: int = 129

    #: Brent's method stops once u is known to this absolute tolerance.
    kappa_tolerance: float = 1e-12

    # --- The class shares at one κ: a concave maximisation on the 3-simplex ---

    #: A share at 0 is optimal when its KKT multiplier is at most the junction count. This is the relative
    #: rounding allowed in that comparison, for float64 sums over many junction groups.
    kkt_slack: float = 1e-12

    #: On an edge of the simplex (one share at 0), the maximum is root-found on ``t ∈ [margin, 1 − margin]``,
    #: where the slope is finite, to this tolerance.
    edge_margin: float = 1e-15

    #: An interior maximum is found by a log-barrier ascent. The barrier's weight starts at this share of the
    #: junction count, ...
    barrier_start: float = 1e-2

    #: ... is multiplied by this between stages, ...
    barrier_shrink: float = 0.1

    #: ... and ends at this absolute weight, where the barrier's own optimality gap (3 × the weight, for three
    #: shares) is below the log-likelihood's rounding.
    barrier_floor: float = 1e-11

    #: Newton's method within a stage takes at most this many steps (a termination bound) ...
    newton_steps: int = 100

    #: ... and stops when a step's predicted gain falls below this share of the objective.
    newton_tolerance: float = 1e-13

    #: A Newton step goes at most this fraction of the way to the simplex's boundary, so every share stays
    #: positive for the barrier (the fraction-to-boundary rule).
    boundary_fraction: float = 0.99

    #: A step that lowers the barrier objective is shortened by this factor ...
    line_search_shrink: float = 0.5

    #: ... until it is shorter than this, which ends the stage.
    line_search_floor: float = 1e-16

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            kappa_floor=0.0 < self.kappa_floor < 0.25,
            kappa_scan_points=self.kappa_scan_points >= 3,
            kappa_tolerance=self.kappa_tolerance > 0.0,
            kkt_slack=0.0 <= self.kkt_slack < 1.0,
            edge_margin=0.0 < self.edge_margin < 0.5,
            barrier_floor=self.barrier_floor > 0.0,
            barrier_start=self.barrier_start >= self.barrier_floor,
            barrier_shrink=0.0 < self.barrier_shrink < 1.0,
            newton_steps=self.newton_steps >= 1,
            newton_tolerance=self.newton_tolerance > 0.0,
            boundary_fraction=0.0 < self.boundary_fraction < 1.0,
            line_search_shrink=0.0 < self.line_search_shrink < 1.0,
            line_search_floor=0.0 < self.line_search_floor < 1.0,
        )


@dataclass(frozen=True)
class LandscapeConstants:
    """The population gDNA-density landscape (:mod:`rigel.calibration.landscape`): its grid and kernel grouping, and
    the modelling constants its fit was selected with.

    ``knn_scale`` and ``reliability_sd_decades`` are selected against gDNA-shaped data and validated on that component
    only; a second caller inherits them and must say so in its results.
    """

    # --- Computational budgets: discretization, not modelling (finer is slower and more exact) ---

    #: Points on the log-rate grid. Same role as the solver's λ lattice: finer is strictly more faithful and
    #: strictly slower.
    grid_points: int = 260

    #: Kernels are grouped into this many equal-count width bins and each bin convolved once, instead of one
    #: convolution per region. Pure speed: the cost goes from O(n·K²) to O(bins·K²), which is what makes the fit
    #: affordable at genome scale.
    width_bins: int = 12

    # --- Modelling constants ---

    #: Population-resolution scale for ``landscape.knn_widths``, selected by shape against a reference that is
    #: itself validated against ground truth: at 0.5 the fit renders the enriched mode at the width the truth has;
    #: below it the landscape combs, above it the two modes merge and the enriched mass collapses. It is selected by
    #: shape, never by EMD, which is monotone in smoothing at every reference and cannot discriminate.
    knn_scale: float = 0.5

    #: The reliability weight's reference spread, in decades of rate (its variance is ``(this · ln 10)²``). A tuned
    #: constant: it caps how far a confident region can be down-weighted. It is not the kernel's resolution floor,
    #: which is the grid step (about 0.025 decades) and would make the weight more aggressive and the census worse.
    reliability_sd_decades: float = 0.15

    #: The location floor on ``RegionBelief.var_gdna = Var(log f_g)`` (at fixed mass, ``Var(log count)``). The
    #: estimator's resolution wall is one fragment, and a Poisson count ``c`` has ``Var(log c) = 1/c``, so a slot wider
    #: than 1 is below the wall: it has no location and does not train the prior. The identity's value, not tuned.
    located_var: float = 1.0

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            grid_points=self.grid_points >= 2,
            width_bins=self.width_bins >= 1,
            knn_scale=self.knn_scale > 0.0,
            reliability_sd_decades=self.reliability_sd_decades > 0.0,
            located_var=self.located_var > 0.0,
        )


@dataclass(frozen=True)
class CalibrationConstants:
    """The calibration stage's remaining fixed numbers (:mod:`rigel.calibration`): one population rule, one
    bracket, and the floors below which a quantity is read as zero."""

    #: Fewest training regions the gDNA-density landscape is fitted from; below it the population is not a
    #: population, and the refit is skipped (``calibrate._fit_gdna_hyperprior``).
    min_training_regions: int = 5

    #: The pooled gDNA rate's root bracket, as a multiple of the pooled rate (``gdna_density``). The one-sided
    #: root is always below the pooled rate (contamination only inflates it), so any multiple above 1 brackets
    #: it; this is headroom, and the fit's ``bracket_ok`` reports if it ever failed.
    bracket_headroom: float = 10.0

    #: A slot whose effective length exceeds this has opportunity: it can be expressed, anchor the landscape,
    #: and belong to the landscape's grid domain (``calibrate._fit_gdna_hyperprior``).
    opportunity_floor: float = 1e-9

    #: A slot whose unspliced mass exceeds this sequenced something; at or below it, an intergenic or intronic
    #: region with opportunity is the landscape's zero-count anchor.
    mass_floor: float = 1e-12

    #: Rounding a stored fraction or capture efficiency may carry above 1 before ``CalibrationResult`` refuses it.
    fraction_tolerance: float = 1e-9

    #: The gDNA track (``calibration.track``) reports 0 density, or 0 fraction, where its divisor — an effective
    #: length, or a total — is at or below this.
    track_floor: float = 1e-9

    #: The capture efficiency's Poisson log-likelihood reads ``log(max(λ, this))``, so an empty expectation scores
    #: a finite floor rather than −∞ (``capture_efficiency``).
    rate_log_floor: float = 1e-300

    #: Levels on the geometric grid the junction cap is searched on, from the saturated level at the capture weight to
    #: the largest uncapped junction price (``capture_eff_length.fit_junction_cap``). The likelihood is flat where no
    #: junction is at the cap and concave only between the junctions' own bounds, so it is maximised by search and
    #: the first maximum taken. A discretization, not a model choice: at 600 neighbouring levels are under 1 % apart
    #: over any range up to 380-fold.
    junction_cap_grid_points: int = 600

    #: The Jeffreys prior on a Poisson rate, ``ρ^(−½)``, adds this to the intergenic pool's weighted count when the
    #: rate's own posterior width is read as a Gamma, ``Gamma(Σw·n + ½, ΣE)`` (``gdna_density.background_law``): the
    #: reference exponent ψ's Jeffreys arms carry, and the shape that keeps an empty pool's rate finite and its width
    #: proper.
    background_rate_shape: float = 0.5

    #: The background's over-dispersion search (``gdna_density.background_law``) stops when its golden-section bracket
    #: is this share of the starting one, the pool's moment estimate. A discretization, not a model choice: a millionth
    #: of the spread changes no intron's row measurably.
    background_fit_tolerance: float = 1e-6

    #: The introns' departure prior (``intron_departure``) takes each intron's background expectation at this many
    #: midpoint quantiles of the background's count. A discretization, not a model choice: at 64 the estimate is within
    #: 0.4 % of the exact sum over every count, where its own spread on a 2,000-intron library is about 6 %. The locks'
    #: departure (``lock_departure``) reads a departed lock's ``Var(log f)`` at as many quantiles of its share.
    departure_quadrature_nodes: int = 64

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            min_training_regions=self.min_training_regions >= 1,
            bracket_headroom=self.bracket_headroom > 1.0,
            opportunity_floor=self.opportunity_floor >= 0.0,
            mass_floor=self.mass_floor >= 0.0,
            fraction_tolerance=0.0 <= self.fraction_tolerance < 1.0,
            track_floor=self.track_floor > 0.0,
            rate_log_floor=self.rate_log_floor > 0.0,
            junction_cap_grid_points=self.junction_cap_grid_points >= 2,
            background_rate_shape=self.background_rate_shape > 0.0,
            background_fit_tolerance=0.0 < self.background_fit_tolerance < 1.0,
            departure_quadrature_nodes=self.departure_quadrature_nodes >= 1,
        )


@dataclass(frozen=True)
class FragmentLengthConstants:
    """The fragment-length laws (``calibration.fl``)."""

    #: Dirichlet pseudo-count of the empirical-Bayes shrink of gDNA's pool law toward the all-fragment law, with no
    #: cliff: a pool total far above it gives the empirical law, far below it the anchor, and 0 the anchor exactly.
    #: Chosen, not derived. RNA's law does not use it: it is read from its own counts
    #: (``calibration.fl.own_count_law``).
    pool_prior_ess: float = 1000.0

    #: The own-count law smooths its deviation from the all-fragment law by a random walk whose step variance (per
    #: base, on the log scale) is the evidence maximum between these two bounds. The lower is a deviation flat across
    #: any pool's range (a standard deviation of 2e-4 over 500 bp: the anchor's own shape), the upper no smoothing
    #: at all (10 per base), so the two limits are bracketed and the maximum lies between them.
    law_step_min: float = 1e-10
    law_step_max: float = 1e2

    #: The evidence search stops when its bracket on the log step is this narrow: the evidence is flat at its maximum,
    #: and a thousandth in the log step moves no law measurably.
    law_step_tolerance: float = 1e-3

    #: Newton's method on the deviation takes at most this many steps (a termination bound; the objective is convex
    #: with a tridiagonal Hessian and converges in tens) ...
    law_newton_steps: int = 200

    #: ... and stops when a step changes the penalised objective by less than this share of it.
    law_newton_tolerance: float = 1e-12

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            pool_prior_ess=self.pool_prior_ess >= 0.0,
            law_step_min=self.law_step_min > 0.0,
            law_step_max=self.law_step_max > self.law_step_min,
            law_step_tolerance=self.law_step_tolerance > 0.0,
            law_newton_steps=self.law_newton_steps >= 1,
            law_newton_tolerance=0.0 < self.law_newton_tolerance < 1.0,
        )


@dataclass(frozen=True)
class ScoringConstants:
    """Fragment scoring's numerical floors (``rigel.scoring``)."""

    #: The strand probabilities are floored here before their log, so a perfectly stranded library scores a
    #: wrong-strand fragment at log(1e-10) ≈ −23, finite, rather than −∞.
    strand_probability_floor: float = 1e-10

    #: The pruning threshold is floored here before ``−log``, so a threshold of 0 keeps every candidate within a
    #: finite margin (690 nats) rather than reading ``−log 0``.
    pruning_floor: float = 1e-300

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            strand_probability_floor=0.0 < self.strand_probability_floor < 0.5,
            pruning_floor=0.0 < self.pruning_floor < 1.0,
        )


@dataclass(frozen=True)
class QcConstants:
    """Quality-control thresholds: what is warned about, what the QC summary reports, and how the QC report
    (``rigel report``) displays it."""

    #: Fewer spliced strand observations than this draws a warning that the strand estimate may be noisy
    #: (``StrandModels``); zero draws a stronger one.
    strand_min_observations: int = 20

    #: Confidence of the upper credible limit on the wrong-strand rate the strand trainer's log line reports.
    strand_ci_confidence: float = 0.99

    #: An error or warning that lists offending items (duplicate transcript groups, …) shows at most this many.
    report_examples: int = 5

    # --- The QC report's display (``rigel report``): display choices, not measured optima. They set how a headline
    # tile is coloured, how much of a table the page embeds, how far a cell is rounded and how a chart bins and caps
    # what it draws; none reaches a number Rigel computes ---

    #: The mapping tile is good at this share of the BAM's records mapped or more, ...
    mapping_rate_good: float = 0.75

    #: ... a warning at this share or more, and critical below it.
    mapping_rate_warning: float = 0.5

    #: The gDNA tile, which reads the summary's ``fraction_gdna``, is good below this fraction, ...
    gdna_fraction_good: float = 0.1

    #: ... a warning below this one, and critical at it or above.
    gdna_fraction_warning: float = 0.3

    #: The usable-fragments tile is good at this genic share of the fragments or more, and a warning below it.
    genic_fraction_good: float = 0.6

    #: The gDNA-density panel counts the calibration track's regions whose gDNA fraction exceeds this, and its
    #: label states it.
    region_gdna_fraction: float = 0.5

    #: The gene table embeds at most this many expressed genes, the highest TPM first; every gene stays in
    #: ``gene_quant``. The cap keeps the page a sane size.
    gene_table_rows: int = 20_000

    #: Decimals a gene's TPM is rounded to.
    tpm_decimals: int = 2

    #: Decimals a count is rounded to: a gene's, and a reference's gDNA count.
    count_decimals: int = 1

    #: Decimals a reference's gDNA fraction is rounded to.
    fraction_decimals: int = 4

    #: Decimals a fragment-length category's overflow is rounded to, as a percentage.
    overflow_percent_decimals: int = 2

    #: The fragment-length overlay bins lengths this many bp wide before it normalises each pool's density.
    fl_overlay_bin_bp: int = 5

    #: The genome track facets at most this many references, those holding the most gDNA first, and its caption
    #: states the cap; the reference table lists every reference.
    genome_track_refs: int = 24

    #: The genome track draws each reference as this many equal-width bins of its span.
    genome_track_bins: int = 160

    #: Decimals a genome-track bin's position, in Mb, is rounded to.
    position_mb_decimals: int = 4

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            strand_min_observations=self.strand_min_observations >= 0,
            strand_ci_confidence=0.0 < self.strand_ci_confidence < 1.0,
            report_examples=self.report_examples >= 1,
            mapping_rate_good=0.0 < self.mapping_rate_good <= 1.0,
            mapping_rate_warning=0.0 <= self.mapping_rate_warning <= self.mapping_rate_good,
            gdna_fraction_good=0.0 < self.gdna_fraction_good <= 1.0,
            gdna_fraction_warning=self.gdna_fraction_good <= self.gdna_fraction_warning <= 1.0,
            genic_fraction_good=0.0 < self.genic_fraction_good <= 1.0,
            region_gdna_fraction=0.0 <= self.region_gdna_fraction < 1.0,
            gene_table_rows=self.gene_table_rows >= 1,
            tpm_decimals=self.tpm_decimals >= 0,
            count_decimals=self.count_decimals >= 0,
            fraction_decimals=self.fraction_decimals >= 0,
            overflow_percent_decimals=self.overflow_percent_decimals >= 0,
            fl_overlay_bin_bp=self.fl_overlay_bin_bp >= 1,
            genome_track_refs=self.genome_track_refs >= 1,
            genome_track_bins=self.genome_track_bins >= 1,
            position_mb_decimals=self.position_mb_decimals >= 0,
        )

    def examples(self, names) -> str:
        """*names* as an error or warning lists them: the first :attr:`report_examples`, and how many more."""
        names = list(names)
        k = self.report_examples
        return ", ".join(names[:k]) + (f" and {len(names) - k:,} more" if len(names) > k else "")


@dataclass(frozen=True)
class SimulatorConstants:
    """The simulator's fixed numbers (``rigel.sim``). They change what is drawn, so a change re-simulates."""

    #: Truncated fragment lengths are drawn by rejection: each pass draws ``floor(needed × ratio) + extra``
    #: candidates so the truncation usually fills the request in one pass.
    oversample_ratio: float = 1.5
    oversample_extra: int = 10

    #: A ``rigel sim`` scenario's defaults where neither its YAML nor the command line sets them (``sim.Scenario``):
    #: a genome length (bp) and a fragment count small enough that a scenario builds in seconds, and a seed, whose value
    #: is arbitrary: any fixed one makes the scenario reproducible.
    scenario_genome_length: int = 5000
    scenario_seed: int = 42
    scenario_fragments: int = 1000

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            oversample_ratio=self.oversample_ratio >= 1.0,
            oversample_extra=self.oversample_extra >= 1,
            scenario_genome_length=self.scenario_genome_length >= 1,
            scenario_seed=self.scenario_seed >= 0,
            scenario_fragments=self.scenario_fragments >= 1,
        )


@dataclass(frozen=True)
class ResourceConstants:
    """Memory, I/O and logging budgets: each sets speed, memory, file size or log volume, never a result. The
    row-tiled fits (``solve_block_bytes``, ``capture_tile_cell_bytes``) and the alignable window give the same bits at
    every tile height and window."""

    #: Scan workers one BGZF decompression thread keeps fed on a name-sorted BAM: the budget is split by this ratio
    #: (``BamScanConfig.resolved_scan_threads``), so a 4-thread run spends none on decompression and a 16-thread run
    #: spends two.
    scan_workers_per_bgzf_thread: int = 8

    #: Chunks the scan buffer's spill writer holds queued at once. Each stays in memory until it is written, so the
    #: buffer's memory can exceed ``BamScanConfig.buffer_size_bytes`` by up to this many chunks, and the scan waits
    #: when it would queue one more. A resource choice, not a measured optimum.
    spill_queue_chunks: int = 2

    #: Cache-tiling target for the row-tiled fits (the landscape's kernels, the capture efficiency), as a working-set
    #: size; ``simplex_logodds.block_rows`` turns it into rows. It orders no addition: the landscape adds every row
    #: into its group's sum in row order, whatever the tile.
    solve_block_bytes: int = 1 << 20

    #: The bytes per grid cell the capture efficiency's row tile is sized for (``block_rows``'s itemsize): two float64
    #: arrays' worth, while a tile holds four at once (the expectation, its log-likelihood, the alternative's log
    #: weight and its weight) and a temporary fifth at its peak, so a tile peaks at about two and a half times
    #: ``solve_block_bytes``. A resource choice, not a measured optimum. It sets the tile height and no bit: every
    #: sum runs along one object's own row.
    capture_tile_cell_bytes: int = 16

    #: Cuts priced at once when the transcripts' shares are built (``capture_eff_length.object_shares``): each cut's
    #: shares are its own, so the block sets the size of the temporaries (about a dozen float64 arrays of this many
    #: rows, a few tens of MB, where the whole genome's ten million cuts at once took 1.7 GB) and no bit. A resource
    #: choice, timed on the whole genome: smaller blocks pay numpy's call overhead, larger ones the memory.
    share_block_rows: int = 1 << 17

    #: Read size when the index digests its input files.
    digest_chunk_bytes: int = 1 << 20

    #: Bases of an alignable store the index reads per call when it records a reference's unobservable runs: one
    #: read-length bin's counts cost 16 bytes a base, so a window is 160 MB. It sets memory only, since a run is
    #: stitched across window edges.
    alignable_window_bp: int = 10_000_000

    #: The per-fragment annotation table is sized to the buffered fragment count plus this padding, and never
    #: below the minimum capacity; it doubles when full.
    annotation_table_padding: int = 1024
    annotation_table_min_capacity: int = 4096

    #: GTF parsing logs progress every this many exon features.
    gtf_log_interval: int = 100_000

    #: The simulator's FASTQ writer buffers this many records before a write.
    fastq_buffer_records: int = 100_000

    #: The simulator's FASTQ writer compresses with this many pgzip threads when pgzip is installed; a sharded
    #: run writes each shard with one, every shard being a process of its own.
    fastq_gzip_threads: int = 4

    #: The gzip level the simulator writes its FASTQ at. It and the threads above set the speed and the compressed
    #: bytes, never a read: resource choices, not measured optima.
    fastq_gzip_level: int = 3

    #: Chunk size when the simulator concatenates gzip streams.
    file_copy_bytes: int = 4 * 1024 * 1024

    def __post_init__(self) -> None:
        _check_ranges(self, **{f.name: getattr(self, f.name) >= 1 for f in fields(self)})


@dataclass(frozen=True)
class Constants:
    """Every fixed number Rigel's code reads, by component. The one instance is :data:`CONSTANTS`."""

    junction_fit: JunctionFitConstants = field(default_factory=JunctionFitConstants)
    landscape: LandscapeConstants = field(default_factory=LandscapeConstants)
    calibration: CalibrationConstants = field(default_factory=CalibrationConstants)
    fragment_length: FragmentLengthConstants = field(default_factory=FragmentLengthConstants)
    scoring: ScoringConstants = field(default_factory=ScoringConstants)
    qc: QcConstants = field(default_factory=QcConstants)
    simulator: SimulatorConstants = field(default_factory=SimulatorConstants)
    resources: ResourceConstants = field(default_factory=ResourceConstants)


def _check_ranges(section, **ok: bool) -> None:
    """Refuse a section whose values are out of range, naming every offending field."""
    bad = [name for name, valid in ok.items() if not valid]
    if bad:
        raise ValueError(f"{type(section).__name__}: out of range: {', '.join(bad)}")


#: The constants the code reads.
CONSTANTS = Constants()


# ======================================================================
# Pre-computed transcript geometry (not user-configurable)
# ======================================================================


@dataclass
class TranscriptGeometry:
    """The per-transcript effective lengths the EM reads, computed once at the start of
    ``quant_from_buffer`` from ``TranscriptIndex`` and the RNA
    :class:`~rigel.frag_length_model.FragmentLengthModel`. Not user-configurable.

    That model is built by ``FragmentLengthModel.from_pmf`` from ``FLModels.rna_pmf``, which is derived
    from the accumulator payload alone. The effective lengths here and the calibration divisors read the
    same pmf, so a change to it reaches every transcript in the EM, not only calibration.

    Parameters
    ----------
    effective_lengths : np.ndarray
        float64[n_transcripts] — effective transcript lengths (the output lengths, TPM's).
    effective_lengths_em : np.ndarray
        float64[n_transcripts] — the EM's effective lengths, capture-contracted; equal to
        ``effective_lengths`` off capture.
    """

    effective_lengths: np.ndarray
    effective_lengths_em: np.ndarray
