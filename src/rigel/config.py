"""rigel.config — Pipeline configuration dataclasses.

Single source of truth for all tunable parameters.  Frozen dataclasses
ensure immutability after construction.  Compose the sub-configs into
``PipelineConfig`` for clean function signatures.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

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
    mode : str
        Algorithm variant: ``"map"`` (default) or ``"vbem"``.
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
    mode: str = "vbem"
    iterations: int = 1000
    convergence_delta: float = 1e-6
    assignment_mode: str = "sample"
    assignment_min_posterior: float = 0.01
    n_threads: int = 0
    """Number of threads for parallel locus EM (Phase 2B).

    ``0`` (default) → use all available cores (``omp_get_max_threads()``).
    ``1`` → sequential (no OpenMP overhead).
    Any positive value → cap at that many threads.
    Ignored when the C++ extension was built without OpenMP.
    """

    def __post_init__(self):
        if self.mode not in ("map", "vbem"):
            raise ValueError(f"Unknown EM mode: {self.mode!r}")
        if self.assignment_mode not in ("fractional", "map", "sample"):
            raise ValueError(f"Unknown assignment mode: {self.assignment_mode!r}")


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
        (i.e. each overhang base cuts probability by 10×).
    mismatch_log_penalty : float
        Log-penalty per NM mismatch.  Default ``log(0.1) ≈ −2.303``
        (i.e. each mismatch cuts probability by 10×).
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
    chunk_size : int
        Fragments per buffer chunk (default 1M).
    max_memory_bytes : int
        Max memory before disk spill (default 4 GiB).
    spill_dir : Path, str, or None
        Directory for spilled buffer chunks (default None).
    """

    skip_duplicates: bool = True
    include_multimap: bool = True
    max_frag_length: int = 1000
    sj_strand_tag: str | tuple[str, ...] = "auto"
    log_every: int = 1_000_000
    chunk_size: int = 1_000_000
    max_memory_bytes: int = 4 * 1024**3
    spill_dir: Path | str | None = None
    n_scan_threads: int = 0
    n_decomp_threads: int = 4
    """Number of htslib BGZF decompression threads (default 4).

    Controls how many threads htslib uses for BAM decompression,
    independent of the worker threads used for fragment resolution.
    Set to 0 to disable multi-threaded decompression.
    """

    splicing_anchor_tolerance: int = 3
    """Minimum bp clearance K required on each side of an exon-intron
    boundary for a fragment to (a) contribute its EXON|INTRON bit to
    the per-fragment ``obs_mask`` and (b) count as a boundary-crossing
    event in ``u_left`` / ``u_right``.

    Internally enforced as ``q(K) = max(K, 1)`` so that ``K = 0``
    reproduces the pre-2026.05 strict-crossing semantics bit-for-bit.
    The matched ``B_cross(K)`` denominator in the global gDNA density
    estimator uses the same ``q(K)``; the numerator/denominator must
    agree or per-locus :math:`\\eta_g` is biased.

    Default 3 bp removes the great majority of single-bp alignment
    artefacts (soft-clip drift, indel slippage near GT-AG splice
    motifs) at the cost of ~1% of true short-overhang exposure for
    typical 350 bp gDNA fragments. Set to 0 to disable the tolerance
    and reproduce pre-2026.05 calibration outputs exactly.
    """

    def __post_init__(self) -> None:
        if self.splicing_anchor_tolerance < 0:
            raise ValueError(
                f"BamScanConfig.splicing_anchor_tolerance must be >= 0; "
                f"got {self.splicing_anchor_tolerance}."
            )


# ======================================================================
# Top-level pipeline configuration
# ======================================================================


@dataclass(frozen=True)
class CalibrationConfig:
    """Configuration for the v6 gDNA calibration orchestrator
    (:func:`rigel.calibration.calibrate`)."""

    #: Empirical-Bayes evidence strength for the FL-Dirichlet shrinkage
    #: in :func:`rigel.calibration.calibrate`.  The ``rna`` and ``gdna``
    #: per-pool FL distributions are shrunk toward the global FL with
    #: this many pseudo-observations.  Default matches
    #: :data:`rigel.calibration.fl.POOL_EB_PRIOR_ESS`.
    prior_ess: float = 1000.0

    #: Per-component nRNA-suppression weight (0 disables nRNA components
    #: in the prior; 1 treats nRNA on equal footing with mRNA).
    nrna_weight: float = 0.0

    #: Minimum SPLICED-annotated count (``rna``) and gDNA count required
    #: for the pool's per-FL distribution to be flagged ``"good"`` /
    #: ``"weak"`` respectively.  Below ``weak_threshold`` a pool is
    #: flagged ``"unusable"`` and downstream code falls back on the
    #: global FL.  Default matches
    #: :data:`rigel.calibration.fl.POOL_QUALITY_GOOD_THRESHOLD` /
    #: :data:`rigel.calibration.fl.POOL_QUALITY_WEAK_THRESHOLD`.
    pool_quality_good: int = 5_000
    pool_quality_weak: int = 200


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
    """Pre-computed transcript/gene geometry for the EM solver.

    Computed once from ``TranscriptIndex`` + ``FragmentLengthModels`` at the
    start of ``quant_from_buffer``.  Not user-configurable — these are
    derived from the reference and trained models.

    Parameters
    ----------
    effective_lengths : np.ndarray
        float64[n_transcripts] — effective transcript lengths.
    exonic_lengths : np.ndarray
        float64[n_transcripts] — spliced exonic lengths.
    t_to_g : np.ndarray
        int32[n_transcripts] — transcript-to-gene mapping.
    transcript_spans : np.ndarray
        float64[n_transcripts] — genomic transcript spans.
    """

    effective_lengths: np.ndarray
    exonic_lengths: np.ndarray
    t_to_g: np.ndarray
    transcript_spans: np.ndarray
