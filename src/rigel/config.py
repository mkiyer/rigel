"""rigel.config — every setting and every fixed number Rigel's code reads, in one place.

Two trees, both frozen dataclasses, validated on construction:

* :class:`PipelineConfig` — what a run may set (the CLI, a YAML file, an instrument's ``--set``): ``em``,
  ``scan``, ``scoring`` and ``calibration``. :class:`IndexConfig` is the same for ``rigel index``.
* :data:`CONSTANTS` — what the code fixes: the numerical methods, the quality-control thresholds and the resource
  budgets, by component. Not user-configurable; each value is documented with why it has that value, and the code
  reads it by name and never restates it.

What is deliberately NOT here: identifiers and formats (enum codes, bit flags, schema versions, file names), which
live with the format they define; mathematical definitions (gDNA's strand probability ½, Jeffreys' ½), named where
they are used; and the constants the native kernels share with Python, defined once in C++ and exported by
:mod:`rigel.native`.
"""

from __future__ import annotations

import math
import operator
import os
from dataclasses import dataclass, field, fields, replace
from pathlib import Path
from typing import TYPE_CHECKING, ClassVar, Literal

if TYPE_CHECKING:
    import numpy as np


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
        Post-EM fragment assignment (default ``"sample"``). ``"sample"`` gives every fragment to ONE
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
    warm_start: Literal["coverage", "prior", "uniform"] = "coverage"
    """⭐ What the EM's initial ``theta`` is derived FROM.

    ``"coverage"`` (the default and the shipped behaviour) seeds it with each component's unambiguous
    total plus a coverage-weighted share of the ambiguous fragments, and then projects that seed through
    the calibration prior.

    ``"uniform"`` seeds every component equally, so the seed asserts NOTHING — the EM's landing point
    is then a property of the likelihood rather than of where it was put down. ⭐ It is the control that
    separates *"the shipped seed steers the solver into a bad basin"* from *"the objective has one"*.

    ``"prior"`` zeroes the seed, so ``theta`` starts proportional to the prior alone. ⛔ It exists
    because the seed and the prior are two different methods and the projection MULTIPLIES them — a
    coverage-weighted share scaled by a per-transcript allocation derived some other way is neither
    method's answer. ⚠ Only meaningful together with a per-transcript weight: under the shipped
    evidence-proportional rule an all-zero seed leaves the RNA pool at zero and hands the locus to gDNA,
    because ``out[i]`` is proportional to ``raw[i]``."""

    def __post_init__(self):
        if self.mode not in ("map", "vbem"):
            raise ValueError(f"Unknown EM mode: {self.mode!r}")
        if self.assignment_mode not in ("fractional", "sample"):
            raise ValueError(f"Unknown assignment mode: {self.assignment_mode!r}")
        # `Literal` is a type-checker annotation and not a runtime constraint, so an unrecognised
        # value would otherwise fall through to the shipped path — a config field that reads as applied
        # and is not.
        if self.warm_start not in ("coverage", "prior", "uniform"):
            raise ValueError(f"Unknown warm start: {self.warm_start!r}")
        # `None` would seed the draw from the OS's entropy, and one input would stop returning one answer.
        try:
            operator.index(self.seed)
        except TypeError:
            raise TypeError(f"EMConfig.seed must be an integer; got {self.seed!r}.") from None


# ======================================================================
# Fragment scoring configuration
# ======================================================================


@dataclass(frozen=True)
class FragmentScoringConfig:
    """Configuration for fragment scoring penalties.

    All penalties are in log-space.

    Parameters
    ----------
    overhang_log_penalty : float
        Log-penalty per base of overhang.  Default ``log(0.1) ≈ −2.303``
        (i.e. each overhang base cuts probability by 10×).
    mismatch_log_penalty : float
        Log-penalty per NM mismatch.  Default ``log(0.1) ≈ −2.303``
        (i.e. each mismatch cuts probability by 10×).
    """

    #: The CLI's ``--overhang-alpha`` and ``--mismatch-alpha`` default: the probability kept per base of overhang,
    #: and per NM mismatch. The log-penalties below are their logs; the CLI keeps the alpha itself so its default
    #: prints without a log/exp round trip.
    overhang_alpha: ClassVar[float] = 0.1
    mismatch_alpha: ClassVar[float] = 0.1

    overhang_log_penalty: float = math.log(overhang_alpha)
    mismatch_log_penalty: float = math.log(mismatch_alpha)
    pruning_min_posterior: float = 1e-4


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
        Include multimapping reads (default True).
    max_frag_length : int
        Maximum fragment length for histogram models (default 1000).
    sj_strand_tag : str or tuple of str
        BAM tag(s) for splice-junction strand (default ``"auto"``).
    log_every : int
        Log scoring progress (debug level) every N buffered fragments,
        checked once per buffer chunk (default 1M).
    total_threads : int
        Total thread budget available to the scan stage (default 0 = all cores).
        Rigel reserves decompression threads from this budget and uses the
        remainder for scan workers.
    bgzf_threads : int or None
        BGZF decompression threads within the scan thread budget. ``None``, the
        default, DERIVES them from the budget (see
        :meth:`BamScanConfig.resolved_scan_threads`); an integer overrides.
    fragments_per_chunk : int
        Buffered fragments per chunk (default 1M).
    read_name_batch_size : int
        Read-name groups per native scanner input queue item (default 512).
    buffer_size_bytes : int
        Max scan-buffer memory before disk spill (default 2 GiB).
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

    splicing_anchor_tolerance: int = 3
    """Resolver-side splicing-anchor tolerance ``K`` (bp).

    Used only by implicit-splice resolution: a paired-end genomic gap
    can be treated as containing an annotated intron when the intron is
    supported with this many bp of one-sided slack.

    The calibration accumulator does not read this value: compartment,
    splice and strand are recorded directly into the per-region
    channels. It tunes the resolver alone.
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
        if self.splicing_anchor_tolerance < 0:
            raise ValueError(
                f"BamScanConfig.splicing_anchor_tolerance must be >= 0; "
                f"got {self.splicing_anchor_tolerance}."
            )

    def resolved_total_threads(self) -> int:
        """Return the concrete total scan thread budget."""
        if self.total_threads == 0:
            return os.cpu_count() or 1
        return self.total_threads

    def resolved_scan_threads(self) -> tuple[int, int]:
        """Return ``(scan_worker_threads, bgzf_threads)`` within the budget.

        The split is DERIVED from the budget, because the two sides do not scale alike: one
        decompression thread keeps about eight scan workers fed on a name-sorted BAM, so every
        thread given to decompression beyond that ratio is a worker taken away. Measured on the
        18.6M-fragment library, scan seconds by budget and decompression threads — 4: (0) 32.0, (1)
        41.2; 8: (1) 25.8, (2) 26.2, (0) 28.3, (4) 34.5; 16: (2) 17.6, (1) 19.6, (4) 20.1 — so
        ``total // 8`` names the best cell at every budget measured and the ratio, not the count, is
        the thing being set. ``bgzf_threads`` overrides it, and `--scan-bgzf-threads` is the flag.
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
    """Configuration for the calibrator (:func:`rigel.calibration.calibrate`).

    The calibrator is the belief-propagation sweep over the region-boundary chain — a single
    forward-backward pass per solve_chain call, each hop priced by both witnesses' counting plus the
    pair's disagreement beyond it (``hop_price``, `native/transfer_rows.h`); ``sweep_logodds_step`` sets the
    per-region log-odds lattice. See :func:`rigel.calibration.calibrate.calibrate`.
    """

    #: The strand overdispersion is not config: it is 0 for both components, binomial by policy.

    #: The λ lattice's STEP, in nats of log-odds. ψ's grid is ``λ ∈ [−L, L]`` at this spacing, ``K =
    #: round(2L/step) + 1`` points at whatever bracket ``L`` a pass solves on (the floor below, or the
    #: landscape prior's derived demand), so widening the bracket never coarsens the lattice. ONE lattice
    #: serves every consumer — the read-out, the message rows, the AMBIG cube's λ axis, the intron
    #: factory's rows and the composition prior (a second, finer single-strand grid with a linear
    #: regrid between the two measured worse than one grid on every in-scope stratum of both panels,
    #: and the regrid was the reason).
    #:
    #: Dimensionless, so one value serves every depth and genome, and what it guarantees is readable: the
    #: ½-quantile read-out is exact to 1 % of a step once a slot's posterior is wider than the step, and
    #: quantised by at most ``n·f(1−f)·step/4`` fragments below it — at 0.2 a slot's composition is within
    #: 1.25 % of its mass at worst, 0.6 % on average. The ladder (16 conditions, one lattice, ratio to the
    #: retired 60/256 pair on the three in-scope strata, measured when the lattice was unified): 0.69 (30
    #: points) 1.10–1.22×; 0.34 (60) 1.02–1.04×; 0.20 (101) 0.993 / 0.998 / 1.000; 0.146 (138) 0.990;
    #: 0.10 (201) 0.987. The cost grows with ``K`` through the AMBIG cube (``K × (K_t + 2)`` per slot);
    #: 0.2 is the coarsest step that loses nothing.
    sweep_logodds_step: float = 0.2

    #: Log-odds grid FLOOR ``L``: ``λ ∈ [−L, L]`` ⇒ ``f_g ∈ [σ(−L), σ(L)]``. This is the range the
    #: Beta(½,½) reference needs to stay proper, and it is a FLOOR, not the value.
    #:
    #: The bracket actually solved on is ``max(this, the fitted prior's own demand)``, and the second
    #: term is DERIVED rather than chosen:
    #: `calibration.landscape.DensityLandscape.required_logodds_window`. ψ evaluates the landscape at
    #: ``log ρ = log f + log M − log E`` and can only offer ``f ∈ [σ(−L), σ(L)]``, so a bracket narrower
    #: than the landscape's own log-dynamic range leaves ψ with no coordinate for what its own prior
    #: points at — acutely on a near-zero-gDNA library, where the floor sits orders of magnitude above
    #: the density the prior favours.
    sweep_logodds_window: float = 10.0

    #: The sweep's WORKING SET: the chain is solved one LOCUS BLOCK at a time — the chain cut at every
    #: intergenic region (where message passing ends) and the pieces merged up to this many slots per
    #: block (`calibration.region_chain.locus_blocks`). A PERFORMANCE tunable and nothing else: the
    #: answer is the same for every value (ψ's read-out is chunk-exact, gated). ``None`` solves the whole
    #: chain as one block.
    #:
    #: It SIZES THE KERNEL'S ARENA, which is what makes it a memory knob rather than a batching one: each
    #: thread holds sixteen ``(slots, grid)`` float64 tables for the block it is solving, so the arena is
    #: ``16 · slots · K · 8 B · threads`` — 1.19 GB at 5,000 slots, the refits' grid of 233 and eight
    #: threads. Measured on the deep library at eight threads, two interleaved rounds, the sweep's seconds
    #: and calibrate's peak by block size: 5,000 — 23.65 / 23.79 s, 8,431 / 8,271 MB; 2,000 — 23.63 /
    #: 23.31 s, 7,687 / 7,558 MB; 1,000 — 23.22 / 23.26 s, 7,703 / 7,368 MB; 500 — 23.54 / 23.21 s,
    #: 7,414 / 7,645 MB. The wall is FLAT and the memory is not, so the per-block overhead this used to
    #: amortise is not measurable beside the arena it pays for; the default is the smallest size whose
    #: memory gain is still larger than the run-to-run noise.
    sweep_block_slots: int | None = 1000

    #: The message policy — what one neighbour tells another, on the two-phase backbone (prepare →
    #: propagate → solve). ``"transfer"`` (:class:`~rigel.calibration.messages.transfer.TransferPolicy`,
    #: the composition transfer) is the shipped default; ``"silent"``
    #: (:class:`~rigel.calibration.messages.silent.SilentPolicy`) sends nothing, so ψ carries each slot's
    #: OWN evidence alone (its two strand counts, its spliced count, the fitted gDNA prior and the intron
    #: factory) — the measured floor every policy is judged against. Messages exist for the slots whose own
    #: solve has no composition channel: unstranded data and the both-stranded (AMBIG) slots, where the
    #: strand likelihood is flat and the local answer is a default rather than a measurement; on stranded
    #: data a sighted exon's own solve is excellent and a message can mostly only disturb it. The standing
    #: is re-derived by ``scripts/design/policy_benchmark.py``, the two halves read apart and never pooled,
    #: and priced on the composition metric by ``calibration_vs_oracle.py --set calibration.message_policy=…``.
    #: An unknown name RAISES: an arm that silently runs a policy other than the one it names is a
    #: benchmark that cannot be trusted. ⛔ Flipping this default is a config default flip — the trigger
    #: that has left instruments dead while the suite stayed green, because the TEST readers install the
    #: policy themselves. Run the instruments, not just the suite.
    message_policy: str = "transfer"

    #: Calibration refit iterations — the prior BOOTSTRAP. Each iteration re-fits the population gDNA
    #: landscape (:class:`~rigel.calibration.landscape.DensityLandscape`) on the current solved gDNA
    #: densities and belief widths, then fully RESETS the belief and re-solves with it. So nothing but
    #: the fitted landscape carries between iterations, and the prior sharpens only where the data has
    #: earned it. ``0`` ⇒ the prior-free pass-0 alone.
    #:
    #: Cost is linear — one landscape fit plus one full sweep each — so lower it if
    #: calibration wall-clock matters more than the last few percent of its accuracy.
    calib_refit_iters: int = 3

    #: ψ's thread budget: the slots of a block are solved by a pool of this many threads (``0``, the
    #: default, is every core — the locus EM's reading of the same number), and the CLI's ``--threads``
    #: sets it beside the scan's and the EM's. A resource budget, not a tunable of the answer: the solve is
    #: bit-identical at every count (`simplex_logodds._solve_regions_logodds_all`).
    n_threads: int = 0

    def __post_init__(self) -> None:
        if self.n_threads < 0:
            raise ValueError(f"CalibrationConfig.n_threads must be >= 0; got {self.n_threads}.")
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
    """What ``rigel index`` may set; the CLI's defaults and :meth:`rigel.index.TranscriptIndex.build`'s are these.

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

    #: Transcript start and end sites within this many bp cluster into one synthetic nascent-RNA span
    #: (``index.create_nrna_transcripts``).
    nrna_merge_tolerance: int = 20

    #: Fewest unique fragments per (chromosome, intron, read length) for a junction to enter the
    #: splice-artifact blacklist from the alignable store. 2 matches the alignable tool's own threshold; 1 admits
    #: singleton artifacts, higher keeps only the most reproducible.
    splice_blacklist_min_count: int = 2

    def __post_init__(self) -> None:
        if self.feather_compression not in ("lz4", "zstd", "uncompressed"):
            raise ValueError(f"Unknown Feather compression: {self.feather_compression!r}")
        if self.gtf_parse_mode not in ("strict", "warn-skip"):
            raise ValueError(f"Unknown GTF parse mode: {self.gtf_parse_mode!r}")
        if self.nrna_merge_tolerance < 0:
            raise ValueError(
                f"IndexConfig.nrna_merge_tolerance must be >= 0; got {self.nrna_merge_tolerance}."
            )
        if self.splice_blacklist_min_count < 1:
            raise ValueError(
                f"IndexConfig.splice_blacklist_min_count must be >= 1; got {self.splice_blacklist_min_count}."
            )


# ======================================================================
# Top-level pipeline configuration
# ======================================================================


@dataclass(frozen=True)
class PipelineConfig:
    """Top-level pipeline configuration composing all sub-configs.

    Pass a single ``PipelineConfig`` to ``run_pipeline`` instead of
    25+ individual keyword arguments.
    """

    em: EMConfig = field(default_factory=EMConfig)
    scan: BamScanConfig = field(default_factory=BamScanConfig)
    scoring: FragmentScoringConfig = field(default_factory=FragmentScoringConfig)
    calibration: CalibrationConfig = field(default_factory=CalibrationConfig)
    # The second pass's multinomial draw. Pass 1 holds every
    #: fragment whose unsequenced gap has more than one surviving explanation; the drain picks one
    #: hypothesis each and re-deposits, and this seeds that draw.
    #:
    #: Deliberately NOT ``em.seed``. They are two independent RNG consumers, and sharing one field
    #: would mean changing the EM's seed silently re-drew every held fragment — so an EM A/B would move
    #: the tally it was being run against.
    second_pass_seed: int = 0
    annotated_bam_path: str | Path | None = None
    emit_locus_stats: bool = False

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
# run settings: a change to any of them is a change to the method, so it re-runs the component's gates and is priced
# on the instruments before it lands. An experiment varies one on a copy (:meth:`Constants.replaced`), never by
# editing a module.


@dataclass(frozen=True)
class JunctionFitConstants:
    """How the genuine-junction strand fit finds its maximum (:func:`rigel.junction_fit.genuine_sense_fraction`).

    The fit is an exact maximum-likelihood problem: the genuine junctions' wrong-strand rate κ, and three class
    shares (genuine, splice artifact, reversed), over the per-junction strand table. Nothing here is a model
    choice: every field sets how the maximum is FOUND, and at these values the fit reaches it on every table
    checked against an independent search (``tests/test_strand_model.py``'s ``TestGenuineKappa``).
    """

    # --- The search over κ: a scan locates the maximum, Brent's method refines it ---

    #: κ is searched on ``(floor, ½ − floor)``. Within the floor of 0 or of ½, κ × (the RNA reads) moves by less
    #: than one read for any library under 10¹² spliced reads, so the answer cannot change.
    kappa_floor: float = 1e-12

    #: The scan's points, uniform in ``u = log(κ/(½ − κ))`` across that range, so the scan is as fine near ½ as
    #: near 0: at 129 points u steps by 0.42 (κ/(½ − κ) by ×1.5). A weakly stranded library's maximum rises out
    #: of a flat all-artifact plateau just below ½, and a coarser scan can step over it — 65 points spaced on
    #: log κ did (``test_a_weakly_stranded_library_is_fit_exactly``). Each point costs one share solve.
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

    ⚠ ``knn_scale`` and ``reliability_sd_decades`` were selected against gDNA-shaped data. They are not
    component-specific by construction, but they have only ever been validated on one component; a second caller
    inherits them and must say so in its results rather than discover it later.
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
    #: below it the landscape combs, above it the two modes merge and the enriched mass collapses. EMD does not
    #: discriminate here — it is monotone in smoothing at every reference — so do not re-select on it.
    knn_scale: float = 0.5

    #: The reliability weight's reference spread, in decades of rate (its variance is ``(this · ln 10)²``).
    #: ⚠ A tuning constant in disguise. It reads as "the kernel resolution floor", but the actual rendering
    #: resolution is the grid step (~0.025 decades), and substituting that makes the weight more aggressive and
    #: the census worse. What it really does is cap how far a confident region can be down-weighted. Changing it
    #: is its own measured experiment; it must not ride along with anything else.
    reliability_sd_decades: float = 0.15

    #: THE LOCATION FLOOR, in the variable every solve reports. The estimator's resolution wall is one fragment
    #: (``max(count, 1)`` centres a kernel; ``count < 1`` is location-free and E-step placed), and a Poisson count
    #: ``c`` has ``Var(log c) = 1/c``, so "below one fragment" is "the log-count is uncertain by more than one
    #: nat²". ``RegionBelief.var_gdna`` is ``Var(log f_g)`` — at fixed mass, ``Var(log count)`` — whatever produced
    #: the solve, so a slot wider than this has no location by the same floor the count rule applies, and does
    #: not train the prior. Not a tuned constant: the identity's value at the wall
    #: (``tests/calibration/test_landscape_training_population.py``).
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
    #: root is always BELOW the pooled rate (contamination only inflates it), so any multiple above 1 brackets
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
        )


@dataclass(frozen=True)
class FragmentLengthConstants:
    """The fragment-length laws (``calibration.fl``, ``frag_length_model``)."""

    #: Dirichlet pseudo-count of the smooth empirical-Bayes shrink of each pool's length law toward the global
    #: law. Not a cliff: a pool total far above it gives the empirical law, far below it the global anchor, and
    #: 0 the anchor exactly. ⚠ Not derived; weighting the pools by their own precision would replace it.
    pool_prior_ess: float = 1000.0

    #: The boundary gDNA law and its mean are refined jointly: at most this many passes (one refresh of the mean
    #: from the boundary law, measured stable) ...
    gdna_mean_passes: int = 2

    #: ... stopping early once the mean moves by less than this many bp.
    gdna_mean_tolerance_bp: float = 0.25

    #: One pseudo-observation spread across the whole length support keeps an unseen length finite without
    #: pulling a short library toward the middle of the histogram (``FragmentLengthModel``).
    unseen_smoothing_ess: float = 1.0

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            pool_prior_ess=self.pool_prior_ess >= 0.0,
            gdna_mean_passes=self.gdna_mean_passes >= 1,
            gdna_mean_tolerance_bp=self.gdna_mean_tolerance_bp > 0.0,
            unseen_smoothing_ess=self.unseen_smoothing_ess > 0.0,
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
    """Quality-control thresholds: what is warned about, and what the QC summary reports."""

    #: Fewer spliced strand observations than this draws a warning that the strand estimate may be noisy
    #: (``StrandModels``); zero draws a stronger one.
    strand_min_observations: int = 20

    #: Confidence of the upper credible limit on the wrong-strand rate the strand trainer's log line reports.
    strand_ci_confidence: float = 0.99

    #: An error or warning that lists offending items (duplicate transcript groups, …) shows at most this many.
    report_examples: int = 5

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            strand_min_observations=self.strand_min_observations >= 0,
            strand_ci_confidence=0.0 < self.strand_ci_confidence < 1.0,
            report_examples=self.report_examples >= 1,
        )


@dataclass(frozen=True)
class SimulatorConstants:
    """The simulator's fixed numbers (``rigel.sim``). They change what is drawn, so a change re-simulates."""

    #: Truncated fragment lengths are drawn by rejection: each pass draws ``ceil(needed × ratio) + extra``
    #: candidates so the truncation usually fills the request in one pass.
    oversample_ratio: float = 1.5
    oversample_extra: int = 10

    #: The capture-probe designer's default bait length in bp (``sim.capture.design``) — a standard
    #: hybrid-capture bait.
    probe_length_bp: int = 120

    def __post_init__(self) -> None:
        _check_ranges(
            self,
            oversample_ratio=self.oversample_ratio >= 1.0,
            oversample_extra=self.oversample_extra >= 1,
            probe_length_bp=self.probe_length_bp >= 1,
        )


@dataclass(frozen=True)
class ResourceConstants:
    """Memory, I/O and logging budgets. None reaches the arithmetic: every value gives the same answer, and only
    speed, memory or log volume changes."""

    #: Scan workers one BGZF decompression thread keeps fed, measured on the deep library
    #: (``BamScanConfig.resolved_scan_threads`` carries the table): the budget is split by this ratio, so a
    #: 4-thread run spends none on decompression and a 16-thread run spends two.
    scan_workers_per_bgzf_thread: int = 8

    #: Cache-tiling target for the row-tiled fits (the landscape's kernels, the capture efficiency), as a
    #: working-set size; ``simplex_logodds._block_rows`` turns it into rows. Every reduction those fits make is
    #: within a row, so the block size cannot reach the arithmetic.
    solve_block_bytes: int = 1 << 20

    #: Read size when the index digests its input files.
    digest_chunk_bytes: int = 1 << 20

    #: The per-fragment annotation table is sized to the buffered fragment count plus this padding, and never
    #: below the minimum capacity; it doubles when full, by at least the minimum growth.
    annotation_table_padding: int = 1024
    annotation_table_min_capacity: int = 4096
    annotation_table_min_growth: int = 1024

    #: GTF parsing logs progress every this many exon features.
    gtf_log_interval: int = 100_000

    #: The simulator's FASTQ writer buffers this many records before a write.
    fastq_buffer_records: int = 100_000

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

    def replaced(self, dotted: str, value) -> "Constants":
        """A copy with one value changed, named ``"section.field"`` — how an experiment or a test varies one
        constant without editing a module. The section re-validates."""
        section_name, _, field_name = dotted.partition(".")
        section = getattr(self, section_name)
        return replace(self, **{section_name: replace(section, **{field_name: value})})


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
    SAME pmf, so a change to it reaches every transcript in the EM, not only calibration.

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
