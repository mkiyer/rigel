"""rigel.config — Pipeline configuration dataclasses.

Single source of truth for all tunable parameters.  Frozen dataclasses
ensure immutability after construction.  Compose the sub-configs into
``PipelineConfig`` for clean function signatures.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal

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
    seed : int or None
        Random seed for reproducibility.
    mode : {"vbem", "map"}
        Algorithm variant (default ``"vbem"``).
    iterations : int
        Maximum EM iterations (default 1000).
    convergence_delta : float
        Convergence threshold for theta updates (default 1e-6).
    assignment_mode : str
        Post-EM fragment assignment mode: ``"fractional"`` (traditional
        EM posterior weights), ``"map"`` (assign to highest-posterior
        component), or ``"sample"`` (draw from posterior distribution).
        Default ``"sample"``.
    assignment_min_posterior : float
        Minimum posterior for a component to be eligible for discrete
        assignment (``map``/``sample`` modes only).  Components below
        this threshold are zeroed before assignment.  Default 0.01.
    """

    seed: int | None = None
    mode: Literal["vbem", "map"] = "vbem"
    iterations: int = 1000
    convergence_delta: float = 1e-6
    assignment_mode: str = "sample"
    assignment_min_posterior: float = 0.01
    n_threads: int = 0
    """Number of threads for parallel locus EM.

    ``0`` (default) → use all available cores (``omp_get_max_threads()``).
    ``1`` → sequential (no OpenMP overhead).
    Any positive value → cap at that many threads.
    Ignored when the C++ extension was built without OpenMP.
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

    gdna_em_llr_bias: float = 0.0
    """Global gDNA false-positive-aversion dial — a pure log-odds (LLR) bias added
    to the gDNA component's per-fragment weight in the locus EM (``0.0`` = neutral,
    the default). A **positive** value favors gDNA at every unspliced fragment by
    the odds factor ``exp(gdna_em_llr_bias)``: it trades the FP-deleterious gDNA→RNA
    *leak* for the FP-safe RNA→gDNA *siphon* (decreased RNA sensitivity). Use it to
    say "only call a fragment RNA when it is sufficiently more likely RNA than
    gDNA." This reaches the *EM* assignment directly (distinct from the calibration
    deconvolution). Units: nats of log-odds (e.g.
    ``log(9) ≈ 2.20`` requires ~9:1 RNA evidence)."""

    def __post_init__(self):
        if self.mode not in ("map", "vbem"):
            raise ValueError(f"Unknown EM mode: {self.mode!r}")
        if self.assignment_mode not in ("fractional", "map", "sample"):
            raise ValueError(f"Unknown assignment mode: {self.assignment_mode!r}")
        # `Literal` is a type-checker annotation and not a runtime constraint, so an unrecognised
        # value would otherwise fall through to the shipped path — a config field that reads as applied
        # and is not.
        if self.warm_start not in ("coverage", "prior", "uniform"):
            raise ValueError(f"Unknown warm start: {self.warm_start!r}")
        if not math.isfinite(float(self.gdna_em_llr_bias)):
            raise ValueError(
                f"EMConfig.gdna_em_llr_bias must be finite; got {self.gdna_em_llr_bias}."
            )


# ======================================================================
# Fragment scoring configuration
# ======================================================================


@dataclass(frozen=True)
class FragmentScoringConfig:
    """Configuration for fragment scoring penalties.

    All penalties are in log-space.  ``None`` for ``gdna_splice_penalties``
    means use the module defaults from ``scoring.py``.

    Parameters
    ----------
    overhang_log_penalty : float
        Log-penalty per base of overhang.  Default ``log(0.1) ≈ −2.303``
        (i.e. each overhang base region_bounds probability by 10×).
    mismatch_log_penalty : float
        Log-penalty per NM mismatch.  Default ``log(0.1) ≈ −2.303``
        (i.e. each mismatch region_bounds probability by 10×).
    gdna_splice_penalties : dict or None
        Per-SpliceType gDNA penalties (int keys → float values).
    """

    overhang_log_penalty: float = math.log(0.1)
    mismatch_log_penalty: float = math.log(0.1)
    gdna_splice_penalties: dict[int, float] | None = None
    pruning_min_posterior: float = 1e-4


# ======================================================================
# BAM scanning and buffering configuration
# ======================================================================


#: Scan workers one BGZF decompression thread keeps fed, measured on the deep library
#: (`BamScanConfig.resolved_scan_threads` carries the table): the budget is split by this ratio, so a
#: 4-thread run spends none on decompression and a 16-thread run spends two.
_SCAN_WORKERS_PER_BGZF_THREAD = 8


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
        Log progress every N read-name groups (default 1M).
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
        decompression thread keeps about eight scan workers fed on a coordinate-sorted BAM, so every
        thread given to decompression beyond that ratio is a worker taken away. Measured on the
        18.6M-fragment library, scan seconds by budget and decompression threads — 4: (0) 32.0, (1)
        41.2; 8: (1) 25.8, (2) 26.2, (0) 28.3, (4) 34.5; 16: (2) 17.6, (1) 19.6, (4) 20.1 — so
        ``total // 8`` names the best cell at every budget measured and the ratio, not the count, is
        the thing being set. ``bgzf_threads`` overrides it, and `--scan-bgzf-threads` is the flag.
        """
        total = self.resolved_total_threads()
        want = (
            total // _SCAN_WORKERS_PER_BGZF_THREAD
            if self.bgzf_threads is None
            else self.bgzf_threads
        )
        bgzf = min(want, max(total - 1, 0))
        scan_workers = max(1, total - bgzf)
        return scan_workers, bgzf


# ======================================================================
# Top-level pipeline configuration
# ======================================================================


@dataclass(frozen=True)
class CalibrationConfig:
    """Configuration for the calibrator (:func:`rigel.calibration.calibrate`).

    The calibrator is the belief-propagation sweep over the region-boundary chain — a single
    forward-backward pass per solve_chain call, each hop priced by both witnesses' counting plus the
    pair's disagreement beyond it (``hop_price``, `native/transfer_rows.h`); ``sweep_logodds_step`` sets the
    per-region log-odds lattice. See :func:`rigel.calibration.calibrate.calibrate`.
    """

    #: The strand-overdispersion CEILING is not config: it lives as the single asserted constant
    #: ``_CEIL_ALPHA_BETA`` in :mod:`rigel.calibration.gdna_strand`, next to the estimator it
    #: parameterises. Nor is a shrinkage target — neither component shrinks toward a constant, they
    #: shrink toward EACH OTHER by their own measured informations
    #: (``gdna_strand.reconcile_overdispersions``).

    #: The λ lattice's STEP, in nats of log-odds. ψ's grid is ``λ ∈ [−L, L]`` at this spacing, ``K =
    #: round(2L/step) + 1`` points at whatever bracket ``L`` a pass solves on (the floor below, or the
    #: landscape prior's derived demand), so widening the bracket never coarsens the lattice. ONE lattice
    #: serves every consumer — the read-out, the message rows, the AMBIG cube's λ axis, the intron
    #: factory's rows and the composition prior (ruled 2026-09-13: a second, finer single-strand grid
    #: with a linear regrid between the two measured worse than one grid on every in-scope stratum of
    #: both panels, and the regrid was the reason).
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
    #: answer is the same for every value (ψ's read-out is chunk-exact, gated), so it trades the
    #: per-sweep memory — a block's ``(slots, K)`` arrays instead of the whole chain's — against the
    #: per-block overhead. ``None`` solves the whole chain as one block. The default sits on the flat part
    #: of the sweep's peak-allocation and wall-time curves, measured over block sizes on the human chain
    #: (2.09M slots, about 420 blocks at 5,000); every size is bit-identical.
    sweep_block_slots: int | None = 5000

    #: Which (counts, exposure) pair the pooled gDNA background estimators take.
    #: ``"contained"`` (the default) pools the CONTAINED count over the gDNA contained effective
    #: length — unbiased, since ``E[count] = rho·E_contained``, but the fragment-length pmf enters the
    #: DIVISOR. ``"measured_total"`` pools the START/END banks over the region's own LENGTH
    #: (`calibration.total_abundance.region_counts_and_exposure`) — ``E[S] = rho·ell`` for EVERY
    #: fragment length, so no pmf enters at all, and double-walled regions are excluded as honestly not
    #: model-free. Both are pooled as a ratio of SUMS and both feed the SAME conjugate
    #: ``Gamma(Σcounts + ½, Σexposure)``, so this swaps the pair and not the estimator.
    #:
    #: The two agree off capture. Under capture the gDNA pmf is itself capture-distorted, so the
    #: contained divisor is mis-estimated and that pair under-reads the true gDNA rate several-fold,
    #: while a pmf-free exposure is immune. The cost of the pmf-free form is on pools that carry
    #: nascent RNA: both forms over-read there, and the START form takes slightly more of it, because a
    #: fragment starting in an intron and reaching into an exon books a START there. On a clean pool
    #: that term is absent.
    #: ``"measured_total"`` REFUSES to run unless ``calibrate`` is given ``mature_walls`` and
    #: ``boundary_reach`` — a background rate that silently changed estimator is worse than either.
    background_abundance: str = "contained"

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
        if self.background_abundance not in ("contained", "measured_total"):
            raise ValueError(
                "CalibrationConfig.background_abundance must be 'contained' or 'measured_total'; "
                f"got {self.background_abundance!r}."
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
