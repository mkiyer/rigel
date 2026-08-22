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
        # ⚠ `Literal` is a type-checker annotation and not a runtime constraint, so an unrecognised
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
    forward-backward pass per solve_chain call, with the belief-free Poisson disagreement-variance message
    precision (``σ²_msg = σ²_imp + 1/n_src``); ``sweep_n_grid`` sizes the per-region log-odds solve grid. See
    :func:`rigel.calibration.calibrate.calibrate`.
    """

    #: ⭐ The strand-overdispersion PRIOR and CEILING are no longer config: they live as the two
    #: asserted constants ``_PRIOR_ALPHA_BETA`` / ``_CEIL_ALPHA_BETA`` in
    #: :mod:`rigel.calibration.gdna_strand`, next to the estimator they parameterise, and the
    #: shrinkage weight is DERIVED from them rather than asserted. Four fields
    #: (``{gdna,rna}_strand_prior_{alpha_beta,weight}``) collapsed to two constants + one
    #: structural zero.

    #: **Sweep grid resolution** ``K`` for the per-region log-density log-odds solve over ``λ = logit(f_g)``
    #: (``simplex_logodds``, driven by ``sweep.solve_chain``; single-strand regions are exact 1-D, AMBIG
    #: regions marginalize the RNA tilt ``τ``). ``K=60`` matches per-region accuracy at a tractable cost
    #: (``K=20`` over-calls / under-resolves the zero-DNA case).
    sweep_n_grid: int = 60

    #: **Single-strand λ-grid resolution** ``K_ss`` (Fix 3). Single-strand regions solve a cheap 1-D λ grid
    #: (``O(m·K)``), so a fine grid is affordable there and de-quantizes the ``f_g`` readout (the coarse
    #: ``K=60`` snapped f_g to Δf_g≈0.085 steps — the dominant post-Fix-1 error on high-mass exons). Decoupled
    #: from ``sweep_n_grid`` because the AMBIG 2-D ``(λ,τ)`` cube is ``O(m·K·K_t)`` and a fine grid there is a
    #: genome-scale memory risk. Paired with the parabolic sub-grid-mode readout (``simplex_logodds``), which
    #: recovers a 4×-finer-grid accuracy at any ``K``; ``256`` is where the pair saturates (512 adds <1%).
    sweep_n_grid_single_strand: int = 256

    #: **Log-odds grid FLOOR** ``L``: ``λ ∈ [−L, L]`` ⇒ ``f_g ∈ [σ(−L), σ(L)]``. This is the range the
    #: Beta(½,½) REFERENCE needs to stay proper — under it ~0.9 % of the reference's mass lies outside
    #: ``L = 10`` — and it is a FLOOR, not the value.
    #:
    #: ⭐⭐⭐ **THE BRACKET ACTUALLY SOLVED ON IS ``max(this, the fitted prior's own demand)``**, and the
    #: second term is DERIVED rather than chosen: `calibration.landscape.DensityLandscape.required_logodds_window`.
    #: ψ evaluates the landscape at ``log ρ = log f + log M − log E`` and can only offer
    #: ``f ∈ [σ(−L), σ(L)]``; at ``L = 10`` that floor is ``4.540e-05``, which on a zero-gDNA library sits
    #: **370×** above the median density the prior points at, so ψ had no coordinate for what its own
    #: prior was telling it. The demand is the landscape's log-dynamic range — 21.0–23.0 nats on the human
    #: index — so the shipped 10 was less than half of it.
    #:
    #: ⭐⭐ **MEASURED THROUGH `calibrate` ON ALL 16 LADDER CONDITIONS**, Σ|Δ| in object-incidence
    #: fragments against the origin-split oracle, derived ÷ fixed, per stratum::
    #:
    #:     unstranded x capture OFF   0.4339      g00 0.3449 / g05 0.9917 / g50 0.9840 / g98 0.9766
    #:     stranded   x capture OFF   0.9075      g00 0.3242 / g05 0.9914 / g50 0.9820 / g98 0.9732
    #:     stranded   x capture ON    0.9674      g00 0.3951 / g05 0.9972 / g50 1.0077 / g98 0.9739
    #:     unstranded x capture ON    0.9847      (the DEFERRED stratum — reported, not targeted)
    #:
    #: ⭐ All four strata improve; 11 of the 12 IN-SCOPE conditions improve and one regresses
    #: (`g50 ss_0.99 capture_on`, 1.0077). ⭐⭐ `g05` is the tell — it regressed **1.43×** under every
    #: library-wide reference mean ever tried, the blocker that killed that whole family, and here it
    #: improves on both strand settings.
    #: ⛔ The resolution-only control — the lattice doubled at the OLD bracket — moves the OTHER way
    #: (1.0147× / 1.0046× / 1.0036×), so the effect is the BRACKET and not the lattice.
    sweep_logodds_window: float = 10.0

    #: **ψ's composition reference takes its MEAN from the ANNOTATION** — ``m → 0.75`` wherever no annotated
    #: MATURE transcript is continuous across the position (``¬mrna_active``), and the shipped neutral ½
    #: everywhere else (`calibration.simplex_logodds.structural_reference_location`).
    #: ⛔ This line read ``m → σ(L)`` until 2026-08-17, which is the LATTICE-derived form the function it
    #: names records as REPLACED and measured WORST OF ALL (9.31 nats; refute-err 2.0247 against no-prior's
    #: 0.3946). The shipped value is the one derived two paragraphs down, and the cap is all that is left
    #: of ``σ(L)``: ``m = min((a+1)/(a+b+1), σ(L))``, which needs ``L < 1.0986`` to bind.
    #:
    #: ⭐ ``a·log f_g + b·log(1−f_g)`` on the λ grid IS ``Beta(a, b)``, so the reference has always had a
    #: MEAN ``a/(a+b)`` and nobody chose it: Jeffreys' ½ **asserts the library is half gDNA**. This turns
    #: that assertion into a per-slot statement the annotation can actually make.
    #:
    #: ⭐⭐ **THE STRENGTH IS ONE PSEUDO-OBSERVATION AND INTRODUCES NO CONSTANT.** ``strength = logit(m)``,
    #: so a location IS its strength in nats; one pseudo-observation of gDNA takes ``Beta(a,b)`` to
    #: ``Beta(a+1,b)``, mean ``(a+1)/(a+b+1)`` = **0.75** at ``a = b = _JEFFREYS_REF`` — a 3:1 claim,
    #: overturned by ~1.5 stranded fragments, and refuted unstranded by the intron-vs-intergenic density
    #: factor (`density_lambda_factor`, ``τ_fac = 161.4`` at an intron).
    #:
    #: ⭐ **ON, measured 2026-08-16 through `calibrate` on all 16 ladder conditions** against a `base`
    #: re-recorded in the same session — final Σ|Δ| in fragments per stratum **0.930 / 0.908 / 0.925** on
    #: the three IN-SCOPE strata, 0.998 on the deferred one, and the ``g00`` ZERO-gDNA control
    #: **BIT-IDENTICAL on all 8 rows**. Better on all 12 contaminated conditions, worse on none.
    #:
    #: ⛔ ``False`` ⇒ ``location=None`` ⇒ the term is not written at all and every path is BIT-IDENTICAL to
    #: the one before it existed.
    structural_reference: bool = True

    #: **Inner tilt-grid resolution** ``K_t`` for AMBIG regions' RNA tilt ``τ`` (the 2-D ``(λ,τ)`` solve).
    #: ``None`` ⇒ reuse ``sweep_n_grid``.
    sweep_n_tilt: int | None = None

    #: **WHICH (counts, exposure) PAIR THE POOLED gDNA BACKGROUND ESTIMATORS TAKE.**
    #: ``"contained"`` (the default, and bit-identical to the tree before this field existed) pools the
    #: CONTAINED count over the gDNA contained effective length — unbiased, since
    #: ``E[count] = rho·E_contained``, but the fragment-length pmf enters the DIVISOR.
    #: ``"measured_total"`` pools the START/END banks over the region's own LENGTH
    #: (`calibration.total_abundance.region_counts_and_exposure`) — ``E[S] = rho·ell`` for EVERY
    #: fragment length, so no pmf enters at all, and double-walled regions are excluded as honestly not
    #: model-free. Both are pooled as a ratio of SUMS and both feed the SAME conjugate
    #: ``Gamma(Σcounts + ½, Σexposure)``, so this swaps the pair and not the estimator.
    #:
    #: ⭐⭐ **MEASURED on the 16-condition ladder against the `gdna` origin partition's own start rate,
    #: per stratum (`total_abundance_audit.py` arm ⓕ). The verdict is that it is a TIE off capture and a
    #: 1.8–4.3× repair under it:**
    #:
    #: ===============================  ==================  ==============
    #: pool                             shipped / truth     new / truth
    #: ===============================  ==================  ==============
    #: intergenic, capture-OFF          1.0000 / 1.0001     1.0000
    #: intergenic, capture-ON           **0.5627–0.5644**   **1.0818–1.0830**
    #: +introns, capture-OFF            1.1821–1.1823       1.2055–1.2058
    #: +introns, capture-ON             **0.2344–0.2348**   **1.0649–1.0654**
    #: ===============================  ==================  ==============
    #:
    #: ⭐ The capture-ON repair is the point: the gDNA pmf is itself capture-distorted, so the divisor is
    #: mis-estimated and the shipped pair under-reads the true gDNA rate by 1.8× (intergenic) to 4.3×
    #: (+introns), while a pmf-free exposure is immune. ⛔ The +introns rows show the ONE cost, and it is
    #: not the estimator's: those pools carry nascent RNA, both forms over-read, and the START form takes
    #: ~2 pp MORE of it because a fragment starting in an intron and reaching into an exon books a START
    #: there. Where the pool is clean (intergenic) that term is absent.
    #: ⛔ ``"measured_total"`` REFUSES to run unless ``calibrate`` is given ``mature_walls`` and
    #: ``boundary_reach`` — a background rate that silently changed estimator is worse than either.
    background_abundance: str = "contained"

    #: **FIT THE `AbundanceLandscape` AT CALIBRATION INIT** — the pre-pass-0 TOTAL-density field and
    #: its mode census (`calibration.abundance_landscape`). Its inputs are counts and lengths only
    #: (the wall-exact measured totals), so there is NO circularity with the solve; under capture it
    #: is bimodal and supplies `rho_0` (the depleted mode), the span `R` (the mode ratio) and a
    #: per-region enrichment responsibility `w_i` — the three quantities the measured pass-0
    #: reference consumes.
    #: ⭐⭐ **``True`` IS THE DEFAULT SINCE 2026-08-21, when the enrichment NPMLE was retired**: this
    #: object is now the SOLE source of the QC report's gDNA-density panel
    #: (`CalibrationDiagnostics.from_abundance_landscape`), and it is strictly more informative than
    #: what it replaced — basins from the mode census rather than the two tallest peaks of a curve, a
    #: depleted mode picked by an INDEPENDENT anchor measurement, a real training count, and a real
    #: rug (the npmle carried no training points, so the report's rug was always empty).
    #: ⛔ **Nothing in the SOLVE reads it**, so enabling it moves no solved number — the retirement's
    #: gate was a bit-identical deliverable.
    #: ⚠ Without the wall inputs (`mature_walls`, `boundary_reach`) the fit is SKIPPED with a WARNING
    #: and the object is ``None``, so the report simply omits the panel. ⛔ That is deliberately NOT
    #: the same policy as ``background_abundance`` above, which still REFUSES: that pair feeds ψ, so a
    #: missing input there would silently change a number the solve consumes, whereas this one is read
    #: only by the report and the debug bundle.
    abundance_landscape: bool = True

    # **gDNA intron factory**. ``True`` ⇒ deconvolve confident gDNA
    #: from INTRON regions against the intergenic background BEFORE the pass-0 solve: a per-intron
    #: ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)`` λ-factor (introns are off-target ⇒ ρ_bg is their TRUE gDNA
    #: density, a two-sided estimate; deconvolves gDNA, not RNA — strand-free). Resolves the unstranded-intron gDNA the
    #: prior-free pass-0 currently leaves at ~½ (fixes both the zero-gDNA false-positive and the gDNA under-call),
    #: and seeds the Phase-2 hyperprior fit with clean intron gDNA. ``False`` ⇒ byte-identical to the
    #: pre-factory pass-0.
    #:
    #: **DEFAULT ON since 2026-07-23**, once the factor's precision was registered as composition evidence
    #: (``messages.relay`` — ``I_factory``). Before that the factory shifted an intron's own
    #: mode but carried no ``τ``, so the intron had no standing to EMIT and the correction died one hop out
    #: (measured: intron belief +93 %, neighbour ``prec_g`` bit-identical). With the evidence channel wired,
    #: pass-0 vs oracle over the 32-scenario ambig_dense_10mb suite (⚠ DELETED — see
    #: `rigel.sim.capture.sampler`; these numbers stand as recorded and are not reproducible as written):
    #: mwae 0.1361 → 0.0949, corr 0.688 → 0.736,
    #: 20 scenarios better / 1 worse / 11 flat; intron mwae 0.1781 → 0.0117 (its share of suite error 17.0 % →
    #: 1.6 %); every stranded scenario better or flat (R4 clean).
    intron_factory: bool = True

    #: ⭐⭐⭐ **MESSAGE PROPAGATION — the belief-propagation relay between neighbouring slots. DEFAULT ON**
    #: (owner, 2026-08-18), after ~11 days muted. ``False`` installs ``messages.silent.SilentPolicy``, ψ
    #: carries each slot's OWN evidence alone — its two strand counts, its spliced count, the derived
    #: reference, the fitted gDNA prior and the intron factory. ``True`` installs
    #: ``messages.relay.RelayPolicy``, every operator behind its own named switch.
    #:
    #: ⭐ **A MEASUREMENT PUT IT OFF, not a preference.** Measured 2026-08-07 on the 36-condition ladder —
    #: ⚠ RETIRED, rebuilt at 16 conditions on 2026-08-13, so the numbers below stand as recorded and are
    #: not reproducible as written. Muting the message layer is a net IMPROVEMENT on THREE OF THE FOUR
    #: STRATA::
    #:
    #:     stranded   x capture ON    -58.3 %   16/16 conditions better
    #:     stranded   x capture OFF   -43.7 %   16/16 better
    #:     unstranded x capture OFF   -32.1 %   14/16 better
    #:     unstranded x capture ON   +154.8 %    0/16 better
    #:
    #: ⭐ **THE ``/16`` IS SCORED ROWS, NOT CONDITIONS, AND BOTH LABELS ARE CORRECT** (settled 2026-08-17
    #: by reading ``scripts/design/arm_score.py``, not by inference). Its rows are keyed
    #: ``(condition, axis)`` and every per-stratum line predicates on ``not is_g00(k)`` — the zero control
    #: gets its own line. So a 36-condition ladder is 9 rungs x 2 ss x 2 capture, of which one rung is the
    #: control: **8 scored conditions x 2 axes = 16 rows per stratum**, exactly. ⚠ On the 16-condition
    #: ladder the same arithmetic gives **6**, so a future re-price reads ``n/6`` and that is not a
    #: regression in coverage — it is a smaller panel.
    #:
    #: ⛔⛔ **AND THE PRICE IS CONCENTRATED, LARGE, AND ON THE ZERO CONTROL.** The panel total is +99.9 %
    #: worse because that one stratum carries 73 % of the error — and end-to-end it is worse than the
    #: panel suggests. On golden scenarios whose TRUTH IS ZERO gDNA, the false-positive gDNA mass goes::
    #:
    #:     antisense_contained   0.029  ->  89.93     (~3,100x)
    #:     antisense_overlap     0.005  ->   9.58     (~1,900x)
    #:     single_exon_clean     0.008  ->   0.104
    #:     multi_exon_spliced    0.001  ->   0.031
    #:
    #: ⭐ Both AMBIG loci, which is exactly where ``κ = ½`` makes the strand λ-term identically 0 so the
    #: slot has NO own composition evidence and a message is the only source there is.
    #:
    #: ⚠ **So this default is a STUDY CONFIGURATION**, and turning it back on is one word — every
    #: operator inside ``RelayPolicy`` remains individually switchable.
    #:
    #: ⛔ **The planned way out has been MEASURED AND REFUSED.** This comment used to say the exit was to
    #: give that slot its own θ-independent evidence via a fragment-length composition channel. That
    #: channel was built, priced on the drained arm and DELETED (2026-08-10): its answer is not a function
    #: of the fragment-length gap at all — closing the gap to 1e-9 bp leaves it reporting 0.59-0.72 on
    #: libraries whose truth is 0.00-0.57 — because a Gaussian log-likelihood is asymptotically LINEAR in
    #: the composition, so its argmax is a SIGN saturated at a grid endpoint. See `TRAPS.md`.
    #: ⭐⭐⭐ **WHY IT IS BACK ON (owner, 2026-08-18).** The mute was always a STUDY configuration and the
    #: numbers above are stale in two independent ways — the 36-condition ladder they were measured on was
    #: RETIRED on 2026-08-13, and many unrelated improvements have landed since. ⭐ What re-opened it: on
    #: the scoped pass-0 exon population the relay is the ONLY thing that solves a slot with no own
    #: evidence. On unstranded capture-OFF those exons have ``fg_strand == fg_loc`` exactly and a mean
    #: |error| of **0.500** — ψ's uninformative ½, i.e. no own evidence at all — and the relay takes them
    #: to **0.000087**. Whole-chain at the ``g00`` control it is **3,777,038 → 70,246 fragments** (0.019×).
    #: ⛔ **IT IS NOT UNIFORMLY BETTER AND THAT IS THE DEBUGGING TARGET, NOT A REASON TO RE-MUTE**: on
    #: STRANDED CONTAMINATED data it is worse whole-chain — ``g98 ss0.99 capture_off`` 156,767 → 255,263
    #: (1.628×), ``g98 ss0.99 capture_on`` 281,487 → 728,103 (2.587×). The relay's win is concentrated
    #: where the local solve is BLIND (unstranded, and every zero-gDNA control); its loss is where the
    #: strand channel already had the answer.
    #: ⚠ **Flipping this default is a CONFIG DEFAULT FLIP** — the trigger that once left six instruments
    #: dead while the suite stayed green, because the TEST readers install the policy themselves. Run the
    #: instruments, not just the suite (`TRAPS: a-green-suite-hid-five-dead-instruments`).
    message_propagation: bool = True

    #: **Which policy `message_propagation = True` installs** — `"relay"` (the shipped
    #: :class:`~rigel.calibration.messages.relay.RelayPolicy`, every operator behind its own switch) or
    #: `"currency"` (the Stage-3 rebuild, :class:`~rigel.calibration.messages.currency.CurrencyPolicy`,
    #: UNDER DEVELOPMENT — rung by rung on the test chromosome, `docs` carry the plan). Consulted only
    #: when `message_propagation` is True; the shipped default is byte-identical to the tree before the
    #: field existed. ⚠ The A/B between the two policies is THIS one value, which is what keeps the
    #: rebuild from touching the relay (`RelayPolicy` is not modified by the rebuild, so nothing that
    #: works can break).
    message_policy: str = "relay"

    #: **Calibration refit iterations — the prior BOOTSTRAP.** Each iteration re-fits the population gDNA
    #: landscape (:class:`~rigel.calibration.landscape.DensityLandscape`) on the *current* solved gDNA
    #: densities + belief widths, then **fully resets the belief** and re-solves with it. So nothing but the
    #: fitted landscape carries between iterations, and the prior sharpens only where the data has earned it.
    #: ``0`` ⇒ the prior-free pass-0 alone.
    #:
    #: **Default 3, measured (2026-07-28).** The bootstrap converges geometrically — suite mass-weighted
    #: mwae over the 32-condition battery (⚠ RETIRED; the panel is now the 16-condition ladder, so these
    #: stand as recorded and are not reproducible as written) goes 0.0788 → 0.0525 → 0.0486 → **0.0475** → 0.0471 → 0.0468, with
    #: successive increments shrinking 2–3× each step, and it is **monotone on every stratum including the
    #: zero-gDNA false-positive guard** (0.0667 → 0.0109), so extra iterations never trade specificity for
    #: accuracy. Iteration 3 captures **96 %** of the total available gain; past it the increments are below
    #: anything worth acting on. Cost is linear — one landscape fit plus one full sweep each, measured
    #: 46.8 s (1 iter) → 96.0 s (3 iters) on a 118 k-region real cfRNA sample. Lower it if calibration
    #: wall-clock matters more than the last ~10 % of its accuracy.
    calib_refit_iters: int = 3

    #: **gDNA hyperprior STRENGTH** — a temperature on ψ's fitted composition arm
    #: (``calibration.landscape.DensityLandscape``). ``1.0`` is exact Bayes. Below 1 tempers a prior that
    #: is, after all, fitted from *biased* pass-0 output, which is robustness rather than a fudge: it is what
    #: lets real data overcome a wrong prior, and it is the intended control for the one measured failure
    #: direction — on zero-gDNA and capture-OFF libraries the landscape places 0.2–2.4 % of its mass in the
    #: enriched region where the truth has ~0.01–1 %. Affects ONLY the fitted hyperprior refit, never
    #: the pre-solve total-density landscape (which votes on nothing) and never the solve's gDNA messages.
    gdna_prior_strength: float = 1.0

    def __post_init__(self) -> None:
        if self.calib_refit_iters < 0:
            raise ValueError(
                f"CalibrationConfig.calib_refit_iters must be >= 0; got {self.calib_refit_iters}."
            )
        if self.background_abundance not in ("contained", "measured_total"):
            raise ValueError(
                "CalibrationConfig.background_abundance must be 'contained' or 'measured_total'; "
                f"got {self.background_abundance!r}."
            )
        if float(self.gdna_prior_strength) < 0.0:
            raise ValueError(
                "CalibrationConfig.gdna_prior_strength must be >= 0 (0 disables the prior term); "
                f"got {self.gdna_prior_strength}."
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
    # ⭐ The second pass's multinomial draw. Pass 1 holds every
    #: fragment whose unsequenced gap has more than one surviving explanation; the drain picks one
    #: hypothesis each and re-deposits, and this seeds that draw.
    #:
    #: ⚠ **Deliberately NOT ``em.seed``.** They are two independent RNG consumers, and sharing one field
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
    """Pre-computed transcript/gene geometry for the EM solver.

     Computed once from ``TranscriptIndex`` + the RNA :class:`~rigel.frag_length_model.FragmentLengthModel`
     at the start of ``quant_from_buffer``.  Not user-configurable — these are
     derived from the reference and the fitted models.

     ⚠ That model is built by ``FragmentLengthModel.from_pmf`` from
     ``FLModels.rna_pmf``, which since TRAPS: pure-and-length-censored is derived from the accumulator payload alone
    the effective lengths here and the calibration divisors
     read the SAME pmf, so a change to it reaches every transcript in the EM, not only calibration.

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
