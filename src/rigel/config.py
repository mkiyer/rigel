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
    prior_alpha : float
        Flat Dirichlet pseudocount per eligible component (default 0.01).
    prior_gamma : float
        OVR prior scale factor (default 1.0).  0.0 disables OVR.
    mode : str
        Algorithm variant: ``"map"`` (default) or ``"vbem"``.
    iterations : int
        Maximum EM iterations (default 1000).
    convergence_delta : float
        Convergence threshold for theta updates (default 1e-6).
    prune_threshold : float or None
        Post-EM pruning evidence-ratio threshold.  Components with
        zero unambig evidence and an evidence ratio (data / alpha)
        below this value are zeroed out and a single EM iteration
        redistributes the freed mass.
        Default 0.1.  Set to ``None`` or a negative value to disable.
    confidence_threshold : float
        Posterior threshold for high-confidence assignment (default 0.95).
    """

    seed: int | None = None
    prior_alpha: float = 0.01
    prior_gamma: float = 1.0
    mode: str = "map"
    iterations: int = 1000
    convergence_delta: float = 1e-6
    prune_threshold: float | None = 0.1
    confidence_threshold: float = 0.95
    n_threads: int = 0
    """Number of threads for parallel locus EM (Phase 2B).

    ``0`` (default) → use all available cores (``omp_get_max_threads()``).
    ``1`` → sequential (no OpenMP overhead).
    Any positive value → cap at that many threads.
    Ignored when the C++ extension was built without OpenMP.
    """
    tss_window: int = 200
    """Fuzzy TSS grouping window (bp).

    Transcripts whose 5' ends lie within this distance are grouped
    together for the nRNA fraction prior hierarchy.  Set to 0
    for exact-coordinate matching.  Default 200 bp.
    """
    nrna_frac_kappa_global: float | None = None
    """Shrinkage pseudo-count pulling locus-strand nRNA fraction toward the global
    prior (nrna_frac = 0.5).  ``None`` (default) → auto-estimate from the data
    via Method of Moments.  Set a positive value to override.
    """
    nrna_frac_kappa_locus: float | None = None
    """Shrinkage pseudo-count pulling TSS-group nRNA fraction toward the (shrunk)
    locus-strand estimate.  ``None`` → auto-estimate.
    """
    nrna_frac_kappa_tss: float | None = None
    """Shrinkage pseudo-count pulling transcript nRNA fraction toward the (shrunk)
    TSS-group estimate.  Also sets the effective sample size (κ) of the
    final Beta prior passed to the EM.
    ``None`` → auto-estimate.
    """

    # -- MoM estimation advanced knobs --
    nrna_frac_mom_min_evidence_global: float = 50.0
    """Minimum fragment evidence for a locus-strand group to be
    included in the global Method-of-Moments κ estimate."""
    nrna_frac_mom_min_evidence_locus: float = 30.0
    """Minimum fragment evidence for a TSS group to be included
    in the locus-level MoM κ estimate."""
    nrna_frac_mom_min_evidence_tss: float = 20.0
    """Minimum fragment evidence for a transcript to be included
    in the TSS-level MoM κ estimate."""
    nrna_frac_kappa_min: float = 2.0
    """Lower clamp for MoM-estimated κ values."""
    nrna_frac_kappa_max: float = 200.0
    """Upper clamp for MoM-estimated κ values."""
    nrna_frac_kappa_fallback: float = 5.0
    """Fallback κ when too few features pass the evidence filter."""
    nrna_frac_kappa_min_obs: int = 20
    """Minimum number of features required for MoM κ estimation;
    fewer triggers the fallback."""

    # -- gDNA EB shrinkage knobs --
    gdna_kappa_chrom: float | None = None
    """Shrinkage pseudo-count pulling chromosome gDNA rate toward the
    global estimate.  ``None`` (default) → auto-estimate via MoM."""
    gdna_kappa_locus: float | None = None
    """Shrinkage pseudo-count pulling locus gDNA rate toward the
    (shrunk) chromosome estimate.  ``None`` → auto-estimate."""
    gdna_mom_min_evidence_chrom: float = 50.0
    """Minimum fragment evidence for a chromosome to contribute to
    the gDNA MoM κ_chrom estimate."""
    gdna_mom_min_evidence_locus: float = 30.0
    """Minimum fragment evidence for a locus to contribute to
    the gDNA MoM κ_locus estimate."""

    def __post_init__(self):
        if self.mode not in ("map", "vbem"):
            raise ValueError(f"Unknown EM mode: {self.mode!r}")


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
        Log-penalty per base of overhang.  Default ``log(0.01) ≈ −4.605``
        (i.e. each overhang base cuts probability by 100×).
    mismatch_log_penalty : float
        Log-penalty per NM mismatch.  Default ``log(0.1) ≈ −2.303``
        (i.e. each mismatch cuts probability by 10×).
    gdna_splice_penalties : dict or None
        Per-SpliceType gDNA penalties (int keys → float values).
    """

    overhang_log_penalty: float = math.log(0.01)
    mismatch_log_penalty: float = math.log(0.1)
    gdna_splice_penalties: dict[int, float] | None = None


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
    strand_prior_kappa : float
        Strand model prior pseudocount κ.  Beta prior is ``Beta(κ/2, κ/2)``
        which shrinks toward 0.5 (max entropy).  Default 2.0 → uniform
        ``Beta(1, 1)``.
    log_every : int
        Log progress every N read-name groups (default 1M).
    chunk_size : int
        Fragments per buffer chunk (default 1M).
    max_memory_bytes : int
        Max memory before disk spill (default 2 GiB).
    spill_dir : Path, str, or None
        Directory for spilled buffer chunks (default None).
    """

    skip_duplicates: bool = True
    include_multimap: bool = True
    max_frag_length: int = 1000
    sj_strand_tag: str | tuple[str, ...] = "auto"
    strand_prior_kappa: float = 2.0
    log_every: int = 1_000_000
    chunk_size: int = 1_000_000
    max_memory_bytes: int = 2 * 1024**3
    spill_dir: Path | str | None = None
    n_scan_threads: int = 0


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
    annotated_bam_path: str | Path | None = None


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
    mean_frag : float
        Mean fragment length from the trained model.
    transcript_spans : np.ndarray
        float64[n_transcripts] — genomic transcript spans.
    """

    effective_lengths: np.ndarray
    exonic_lengths: np.ndarray
    t_to_g: np.ndarray
    mean_frag: float
    transcript_spans: np.ndarray
