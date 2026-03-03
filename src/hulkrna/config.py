"""hulkrna.config — Pipeline configuration dataclasses.

Single source of truth for all tunable parameters.  Frozen dataclasses
ensure immutability after construction.  Compose the sub-configs into
``PipelineConfig`` for clean function signatures.
"""

from __future__ import annotations

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
        Post-EM pruning evidence-ratio threshold (None = disabled).
    confidence_threshold : float
        Posterior threshold for high-confidence assignment (default 0.95).
    """

    seed: int | None = None
    prior_alpha: float = 0.01
    prior_gamma: float = 1.0
    mode: str = "map"
    iterations: int = 1000
    convergence_delta: float = 1e-6
    prune_threshold: float | None = None
    confidence_threshold: float = 0.95

    def __post_init__(self):
        if self.mode not in ("map", "vbem"):
            raise ValueError(f"Unknown EM mode: {self.mode!r}")


# ======================================================================
# Fragment scoring configuration
# ======================================================================


@dataclass(frozen=True)
class FragmentScoringConfig:
    """Configuration for fragment scoring penalties.

    All penalties are in log-space.  ``None`` means use the module
    defaults from ``scoring.py``.

    Parameters
    ----------
    overhang_log_penalty : float or None
        Log-penalty per base of overhang.
    mismatch_log_penalty : float or None
        Log-penalty per NM mismatch.
    gdna_splice_penalties : dict or None
        Per-SpliceType gDNA penalties (int keys → float values).
    """

    overhang_log_penalty: float | None = None
    mismatch_log_penalty: float | None = None
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
        Include multimapping reads (default False).
    max_frag_length : int
        Maximum fragment length for histogram models (default 1000).
    sj_strand_tag : str or tuple of str
        BAM tag(s) for splice-junction strand (default ``"auto"``).
    min_spliced_observations : int
        Minimum spliced observations for strand model (default 10).
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
    include_multimap: bool = False
    max_frag_length: int = 1000
    sj_strand_tag: str | tuple[str, ...] = "auto"
    min_spliced_observations: int = 10
    log_every: int = 1_000_000
    chunk_size: int = 1_000_000
    max_memory_bytes: int = 2 * 1024**3
    spill_dir: Path | str | None = None


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
