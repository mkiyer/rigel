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
    """Configuration for the calibrator (:func:`rigel.calibration.calibrate`).

    The calibrator is the **belief-propagation sweep** over the region↔boundary chain — a single
    forward-backward pass per node_sweep call, with the belief-free Poisson disagreement-variance message
    precision (``σ²_msg = σ²_imp + 1/n_src``); ``sweep_n_grid`` sizes the per-node log-odds solve grid. See
    :func:`rigel.calibration.calibrate.calibrate`.
    """

    #: **gDNA strand-overdispersion prior** (advanced). The gDNA per-region sense rate is
    #: ``Beta(a, a)``; this is that symmetric shape ``a`` (= α = β). The fitted overdispersion is
    #: shrunk toward ``od₀ = 1/(2·a + 1)`` for sparse/low-signal libraries. The prior should be the
    #: **near-binomial null**: gDNA strand is intrinsically random (≈50/50 per region), so the expected
    #: overdispersion is small (only secondary mappability/GC/PCR effects); the MoM measures any real
    #: excess from the data. ``a = 14`` ⇒ ``od₀ ≈ 0.034``. (The old ``a = 3`` ⇒ ``od₀ ≈ 0.143`` was an
    #: over-conservative "FP-safe floor" that pulled EVERY node toward 0.5 — under-calling clean-gDNA
    #: nodes and over-calling clean-RNA nodes — because an inflated od widens the gDNA Beta-Binomial and
    #: erases its specificity at its own mean ½. Confirmed at the unit level: a pure-gDNA 50/50 node
    #: solves to f_g≈1 at od=0 but ≈0.7 at od=0.1.) ``a = 2`` ⇒ ``od = 0.2`` is the **most overdispersion
    #: allowed** (the fit is capped there), so values below 2 are rejected; larger ``a`` ⇒ less.
    gdna_strand_prior_alpha_beta: float = 14.0

    #: Strength of the overdispersion prior, in effective seed-node units (advanced). The prior is
    #: worth this many seed nodes; libraries with far more informative seed nodes follow the fit,
    #: sparse ones shrink toward the prior. Replaces the old hard min-seed-node / significance
    #: gates with continuous shrinkage.
    gdna_strand_prior_weight: float = 30.0

    #: **RNA strand-overdispersion prior** (advanced). The exact twin of
    #: ``gdna_strand_prior_alpha_beta`` for the *RNA* strand Beta-Binomial (fitted from boundary-side
    #: spliced counts). Kept at the **same default as gDNA** so that under sparse data both
    #: components collapse to the same distribution — that symmetry is what keeps an unstranded node
    #: uninformative (a gDNA-only overdispersion biases the deconvolution toward RNA). RNA-spliced
    #: strand is motif-deterministic ⇒ also near-binomial, so the same near-binomial null applies. Same
    #: ``a ≥ 2`` rule (``Beta(2,2)``, od=0.2, is the most overdispersion allowed).
    rna_strand_prior_alpha_beta: float = 14.0

    #: Strength of the RNA overdispersion prior, in effective seed-node units (advanced). Twin of
    #: ``gdna_strand_prior_weight``; same default.
    rna_strand_prior_weight: float = 30.0

    #: **Sweep grid resolution** ``K`` for the per-node log-density log-odds solve over ``λ = logit(f_g)``
    #: (``simplex_logodds``, driven by ``bp_solver.node_sweep``; single-strand nodes are exact 1-D, AMBIG
    #: nodes marginalize the RNA tilt ``τ``). ``K=60`` matches per-node accuracy at a tractable cost
    #: (``K=20`` over-calls / under-resolves the zero-DNA case).
    sweep_n_grid: int = 60

    #: **Single-strand λ-grid resolution** ``K_ss`` (Fix 3). Single-strand nodes solve a cheap 1-D λ grid
    #: (``O(m·K)``), so a fine grid is affordable there and de-quantizes the ``f_g`` readout (the coarse
    #: ``K=60`` snapped f_g to Δf_g≈0.085 steps — the dominant post-Fix-1 error on high-mass exons). Decoupled
    #: from ``sweep_n_grid`` because the AMBIG 2-D ``(λ,τ)`` cube is ``O(m·K·K_t)`` and a fine grid there is a
    #: genome-scale memory risk. Paired with the parabolic sub-grid-mode readout (``simplex_logodds``), which
    #: recovers a 4×-finer-grid accuracy at any ``K``; ``256`` is where the pair saturates (512 adds <1%).
    sweep_n_grid_single_strand: int = 256

    #: **Log-odds grid window** ``L``: ``λ ∈ [−L, L]`` ⇒ ``f_g ∈ [σ(−L), σ(L)]`` (``L=10`` ⇒
    #: ``[4.5e-5, 1−4.5e-5]``, bracketing the vertex mass Phase 0 measured).
    sweep_logodds_window: float = 10.0

    #: **Inner tilt-grid resolution** ``K_t`` for AMBIG nodes' RNA tilt ``τ`` (the 2-D ``(λ,τ)`` solve).
    #: ``None`` ⇒ reuse ``sweep_n_grid``.
    sweep_n_tilt: int | None = None

    #: **Relay fold quadrature** (the coherent ``(λ,θ)`` EP relay — ``bp_solver._fold_lambda``,
    #: ``docs/calibration/dof_pie_relay_implementation_plan.md`` §6.1). The relay folds each incoming message
    #: onto the running ``λ`` belief by a two-stage self-correcting grid moment-match: a COARSE grid locates the
    #: posterior peak + curvature, a FINE grid re-centered on it (``refine`` iterations) resolves the moments.
    #: These are **numerical-resolution knobs, NOT model dials** — the same class as ``sweep_n_grid``. The
    #: prototype (``scripts/debug/dof_pie_relay_check.py`` C8) tracks a 4001-point reference to 2.5e-3 at these
    #: defaults; they are conservative, trade only against compute, and are to be hardened on genome-scale data.
    fold_coarse_k: int = 33
    fold_fine_k: int = 33
    fold_sigma_coverage: float = 6.0
    fold_refine_iters: int = 3

    #: **NPMLE bandwidth** ``h`` (decades) for the Fixed-Kernel Poisson-lognormal Mixture NPMLE
    #: (``calibration.npmle.DensityNPMLE``) — shared by both its uses (the enrichment fit for σ²_transfer and
    #: the gDNA-hyperprior refit). ``P(log ρ)`` is a mixture of fixed-width Gaussian kernels on a log-rate grid,
    #: fit by plain EM — deterministic, no spline. The fixed ``h`` is the KDE-style bandwidth: it forbids any
    #: peak sharper than ``h`` (smooth, never a bed-of-nails) and the projected prior is extremely weak
    #: (``n_eff≈0.15`` pseudo-obs), so strand + messages dominate. Grid size / EM iters are perf-only knobs
    #: left at the fitter's defaults. Design: ``docs/calibration/npmle_struggles.md``.
    npmle_bandwidth: float = 0.15

    #: **Aggregate DNA-background floor** (`calibration.background_reference`,
    #: ``docs/calibration/background_reference_derivation.md``). Measure the genome-wide DNA background as a
    #: pooled scalar ``(log ρ_bg, σ_bg)`` and apply it as a ONE-SIDED log-floor in the gDNA prior — data-driven
    #: crush protection that is DORMANT for a DNA-free / fully-depleted library (never manufactures gDNA) and
    #: never a scale/denominator. ``False`` ⇒ no floor (the pre-background behaviour).
    background_floor: bool = True

    #: Background region set: intergenic-only (``False``, sim-safe — the sim's unrealistically abundant nascent
    #: contaminates introns) vs intergenic + intron (``True``, the real-data path — reclaims the introns' huge
    #: aggregate span for resolution; real nascent is sparse). Explore ``True`` on real data.
    background_include_introns: bool = False

    #: Optional upper MAD fence (in log-density) that drops nascent-contaminated intron outliers from the
    #: background pool before aggregation; ``None`` ⇒ no trim. Only meaningful with ``background_include_introns``.
    background_robust_trim_mad: float | None = None

    #: **gDNA intron factory** (`docs/calibration/gdna_intron_factory_design.md`). ``True`` ⇒ peel confident gDNA
    #: from INTRON nodes against the intergenic background BEFORE the pass-0 solve: a per-intron
    #: ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` λ-factor (introns are off-target ⇒ ρ_bg is their TRUE gDNA
    #: density, a two-sided estimate; peels gDNA, not RNA — strand-free). Resolves the unstranded-intron gDNA the
    #: prior-free pass-0 currently leaves at ~½ (fixes both the zero-gDNA false-positive and the gDNA under-call),
    #: and seeds the Phase-2 hyperprior fit with clean intron gDNA. ``False`` (default) ⇒ byte-identical to the
    #: pre-factory pass-0.
    intron_factory: bool = False

    #: **Calibration refit iterations.** Pass-0 bootstraps from TOTAL density: an extremely weak gDNA-rate
    #: prior + the belief-free projection ``σ²_transfer`` ⇒ deliberately GENTLE messages that nudge rather than
    #: ruin the solve. Each refit re-estimates the population model P(ρ) (`DensityNPMLE`) on the *solved* gDNA
    #: counts + belief width, so the prior sharpens where the data has earned it; message precision needs no
    #: refit (it tracks the running belief). ``0`` ⇒ the pass-0 bootstrap alone (no refit).
    calib_refit_iters: int = 1

    #: **gDNA composition prior representation (Role B).** ``False`` (default) = the EM Poisson-mixture NPMLE.
    #: ``True`` = the **additive, occupancy-weighted, fixed-bandwidth KDE** on the deconvolved gDNA + a weak
    #: 1-pseudo-observation floor (``gdna_kde_restore_plan.md``): every node deposits its own kernel so a
    #: capture-enriched minority is never competed away, and the intergenic flood is a single weak floor rather
    #: than an ``n_regions`` tower. Affects ONLY the gDNA hyperprior fit (Role B) — NOT the enrichment NPMLE
    #: (Role A / σ²_transfer) and NOT the solve's gDNA messages. Experimental; being A/B-gated before default-on.
    gdna_prior_additive: bool = False

    def __post_init__(self) -> None:
        if self.calib_refit_iters < 0:
            raise ValueError(
                f"CalibrationConfig.calib_refit_iters must be >= 0; got {self.calib_refit_iters}."
            )
        if not (0.0 < float(self.npmle_bandwidth) < 5.0):
            raise ValueError(
                "CalibrationConfig.npmle_bandwidth must be in (0, 5) decades; "
                f"got {self.npmle_bandwidth}."
            )
        if self.background_robust_trim_mad is not None and float(self.background_robust_trim_mad) <= 0.0:
            raise ValueError(
                "CalibrationConfig.background_robust_trim_mad must be > 0 (or None); "
                f"got {self.background_robust_trim_mad}."
            )
        if self.sweep_n_grid < 2:
            raise ValueError(
                f"CalibrationConfig.sweep_n_grid must be >= 2; got {self.sweep_n_grid}."
            )
        if self.sweep_n_grid_single_strand < 2:
            raise ValueError(
                "CalibrationConfig.sweep_n_grid_single_strand must be >= 2; "
                f"got {self.sweep_n_grid_single_strand}."
            )
        if not (float(self.sweep_logodds_window) > 0.0):
            raise ValueError(
                "CalibrationConfig.sweep_logodds_window (L) must be > 0; "
                f"got {self.sweep_logodds_window}."
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
