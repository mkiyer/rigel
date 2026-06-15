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
    """Number of threads for parallel locus EM (Phase 2B).

    ``0`` (default) → use all available cores (``omp_get_max_threads()``).
    ``1`` → sequential (no OpenMP overhead).
    Any positive value → cap at that many threads.
    Ignored when the C++ extension was built without OpenMP.
    """
    gdna_em_llr_bias: float = 0.0
    """Global gDNA false-positive-aversion dial — a pure log-odds (LLR) bias added
    to the gDNA component's per-fragment weight in the locus EM (``0.0`` = neutral,
    the default). A **positive** value favors gDNA at every unspliced fragment by
    the odds factor ``exp(gdna_em_llr_bias)``: it trades the FP-deleterious gDNA→RNA
    *leak* for the FP-safe RNA→gDNA *siphon* (decreased RNA sensitivity). Use it to
    say "only call a fragment RNA when it is sufficiently more likely RNA than
    gDNA." Where the calibration ``gdna_deconv_quantile`` tilts the *strand/count
    deconvolution* (an uncertainty-aware FP-rate quantile), this reaches the *EM*
    assignment directly (the two are decoupled). Units: nats of log-odds (e.g.
    ``log(9) ≈ 2.20`` requires ~9:1 RNA evidence)."""

    def __post_init__(self):
        if self.mode not in ("map", "vbem"):
            raise ValueError(f"Unknown EM mode: {self.mode!r}")
        if self.assignment_mode not in ("fractional", "map", "sample"):
            raise ValueError(f"Unknown assignment mode: {self.assignment_mode!r}")
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
    total_threads : int
        Total thread budget available to the scan stage (default 0 = all cores).
        Rigel reserves ``bgzf_threads`` from this budget for BAM decompression
        and uses the remainder for scan workers.
    bgzf_threads : int
        Requested BGZF decompression threads within the scan thread budget
        (default 4).
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
    bgzf_threads: int = 4
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

    Under the fractional cutover the calibration accumulator no longer
    interprets this value \u2014 compartment / splice / strand are
    recorded directly into 12 fractional channels per region. The
    field is preserved here so the resolver can still be tuned
    independently. Default 3 bp matches the pre-cutover behaviour.
    """

    def __post_init__(self) -> None:
        if self.total_threads < 0:
            raise ValueError(f"BamScanConfig.total_threads must be >= 0; got {self.total_threads}.")
        if self.bgzf_threads < 0:
            raise ValueError(f"BamScanConfig.bgzf_threads must be >= 0; got {self.bgzf_threads}.")
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
        """Return ``(scan_worker_threads, bgzf_threads)`` within the budget."""
        total = self.resolved_total_threads()
        bgzf = min(self.bgzf_threads, max(total - 1, 0))
        scan_workers = max(1, total - bgzf)
        return scan_workers, bgzf


# ======================================================================
# Top-level pipeline configuration
# ======================================================================


@dataclass(frozen=True)
class CalibrationConfig:
    """Configuration for the acyclic calibrator (:func:`rigel.calibration.calibrate`).

    The calibrator is a single feed-forward pass (no EM loop, no convergence test). The per-node
    deconvolution is a **strand-vs-count handoff** (see ``strand_deconv``): a strand-observable node
    in a strand-identifiable library takes the strand module's Beta-Binomial posterior median, every
    other node takes the count module's gDNA fraction. The blended point estimate is then read at the
    FP-rate quantile ``gdna_deconv_quantile`` (default ½ ⇒ no shift). The old EM-loop knobs are gone.
    """

    #: **gDNA-deconvolution FP-rate quantile** — the false-positive-aversion dial for the per-node
    #: gDNA/RNA call (Phase 2). Each node's blended gDNA fraction is reported at this quantile of its
    #: posterior: ``gdna_frac ← clip(center + Φ⁻¹(q)·σ)``, where ``center`` is the point estimate
    #: (``w·g_strand + (1−w)·g_count``) and ``σ`` is the combined per-node std (strand BB ⊗ count).
    #:   ``0.5`` = neutral (``Φ⁻¹=0`` ⇒ no shift; the unbiased point estimate — the default, a
    #:     bit-identical no-op);
    #:   ``> 0.5`` = FP-averse — shift toward gDNA, trading the gDNA→RNA leak for the RNA→gDNA siphon.
    #:     The shift is **uncertainty-aware**: a wide-posterior (ambiguous) node moves a lot, a
    #:     confident node barely moves — so it bites where the evidence is genuinely uncertain (unlike
    #:     a uniform log-odds tilt, which moves even a confident-RNA node). ``< 0.5`` leans toward RNA
    #:     (higher sensitivity), symmetric.
    #: It only **widens** (never sharpens) — the count σ is anti-calibrated under capture, so it is
    #: used to inflate the quantile, never to correct the blend (which stays bias-robust at
    #: ``w=(2κ−1)²``). Decoupled from the EM ``gdna_em_llr_bias`` (per-fragment EM component vs
    #: per-node fraction). See ``docs/calibration/phase2_design.md``.
    gdna_deconv_quantile: float = 0.5

    #: Grid resolution of the decode posterior over the gDNA fraction on ``[0, 1]``.
    #: Advanced/technical — 200 is ample for a smooth 1-D posterior.
    n_grid: int = 200

    #: **gDNA strand-overdispersion prior** (advanced). The gDNA per-region sense rate is
    #: ``Beta(a, a)``; this is that symmetric shape ``a`` (= α = β). The fitted overdispersion is
    #: shrunk toward ``od₀ = 1/(2·a + 1)`` for sparse/low-signal libraries (the conservative,
    #: FP-safe "floor"). ``a = 3`` ⇒ ``od₀ ≈ 0.143`` (very overdispersed — good for sparse data);
    #: ``a = 2`` ⇒ ``od = 0.2`` is the **most overdispersion allowed** (the fit is capped there),
    #: so values below 2 are rejected; larger ``a`` ⇒ less overdispersion.
    gdna_strand_prior_alpha_beta: float = 3.0

    #: Strength of the overdispersion prior, in effective seed-node units (advanced). The prior is
    #: worth this many seed nodes; libraries with far more informative seed nodes follow the fit,
    #: sparse ones shrink toward the prior. Replaces the old hard min-seed-node / significance
    #: gates with continuous shrinkage.
    gdna_strand_prior_weight: float = 30.0

    #: **RNA strand-overdispersion prior** (advanced). The exact twin of
    #: ``gdna_strand_prior_alpha_beta`` for the *RNA* strand Beta-Binomial (fitted from boundary-side
    #: spliced counts). Kept at the **same default as gDNA** so that under sparse data both
    #: components collapse to the same distribution — that symmetry is what keeps an unstranded node
    #: uninformative (a gDNA-only overdispersion biases the deconvolution toward RNA). Same ``a ≥ 2``
    #: rule (``Beta(2,2)``, od=0.2, is the most overdispersion allowed).
    rna_strand_prior_alpha_beta: float = 3.0

    #: Strength of the RNA overdispersion prior, in effective seed-node units (advanced). Twin of
    #: ``gdna_strand_prior_weight``; same default.
    rna_strand_prior_weight: float = 30.0

    #: **Strand-information half-trust scale** ``I₀`` for the strand→count blend (advanced). The single
    #: blend weight is ``w = I/(I+I₀)`` with ``I = N·(2κ−1)²`` (the per-node strand information): the
    #: strand cleans a node's count by ``w·g_strand + (1−w)·1`` (``strand_deconv.cleaned_gdna_count``).
    #: ``I₀`` is "the strand information at which the strand is half-trusted" (~1 effective discriminating
    #: fragment). ``w→1`` for a confident strand (high κ, decent depth) — it then stops fighting a clearly
    #: good strand model; ``w→0`` at κ≈½ or thin — a no-op, the count floor. Validate on the suites.
    gdna_strand_info_scale: float = 10.0

    #: **Propagation lattice resolution** ``K`` for the 2-simplex grid sum-product
    #: (``simplex_sweep.deconv_regions_sweep``). Separate from ``n_grid`` (the fine 1-D strand-posterior
    #: grid used for the boundary sides): the 2-simplex has ``(K+1)(K+2)/2`` points and the propagation
    #: edge is ``K²`` per transition, so the cost grows ~quartically — ``n_grid=200`` is far too expensive
    #: here. ``K=60`` matches per-node accuracy at the tractable propagation cost (``K=20`` over-calls /
    #: under-resolves the zero-DNA case).
    sweep_n_grid: int = 60

    #: **Iterative-bootstrap pass count** for the propagation path (``CALIBRATION_PLAN_v2`` §2/§5). Each
    #: pass re-fits the gDNA ``var~mean`` + the global density ``ρ_global`` on the previous pass's gDNA
    #: estimate (Pass 0 = the all-gDNA init: every unspliced fragment is gDNA), then re-solves; the loop
    #: stops early once the per-node ``f_g`` stabilizes. A few passes suffice (strand/seed-anchored, not
    #: open EM).
    sweep_max_passes: int = 4

    def __post_init__(self) -> None:
        if not (0.0 < float(self.gdna_deconv_quantile) < 1.0):
            raise ValueError(
                f"CalibrationConfig.gdna_deconv_quantile must be in (0, 1); "
                f"got {self.gdna_deconv_quantile}."
            )
        if self.n_grid < 2:
            raise ValueError(f"CalibrationConfig.n_grid must be >= 2; got {self.n_grid}.")
        if self.sweep_n_grid < 2:
            raise ValueError(
                f"CalibrationConfig.sweep_n_grid must be >= 2; got {self.sweep_n_grid}."
            )
        if self.sweep_max_passes < 1:
            raise ValueError(
                f"CalibrationConfig.sweep_max_passes must be >= 1; got {self.sweep_max_passes}."
            )
        if self.gdna_strand_prior_alpha_beta < 2.0:
            raise ValueError(
                "CalibrationConfig.gdna_strand_prior_alpha_beta must be >= 2.0 "
                "(Beta(2,2) is the most overdispersion allowed, od=0.2); "
                f"got {self.gdna_strand_prior_alpha_beta}."
            )
        if self.gdna_strand_prior_weight < 0.0:
            raise ValueError(
                "CalibrationConfig.gdna_strand_prior_weight must be >= 0; "
                f"got {self.gdna_strand_prior_weight}."
            )
        if self.rna_strand_prior_alpha_beta < 2.0:
            raise ValueError(
                "CalibrationConfig.rna_strand_prior_alpha_beta must be >= 2.0 "
                "(Beta(2,2) is the most overdispersion allowed, od=0.2); "
                f"got {self.rna_strand_prior_alpha_beta}."
            )
        if self.rna_strand_prior_weight < 0.0:
            raise ValueError(
                "CalibrationConfig.rna_strand_prior_weight must be >= 0; "
                f"got {self.rna_strand_prior_weight}."
            )
        if self.gdna_strand_info_scale <= 0.0:
            raise ValueError(
                "CalibrationConfig.gdna_strand_info_scale must be > 0; "
                f"got {self.gdna_strand_info_scale}."
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
    effective_lengths_em : np.ndarray, optional
        float64[n_transcripts] — EM-only effective transcript lengths. When
        omitted, EM uses ``effective_lengths``.
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
    effective_lengths_em: np.ndarray | None = None
