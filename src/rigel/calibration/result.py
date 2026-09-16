"""CalibrationResult — the calibrator's output schema.

    Gate: ``tests/calibration/test_result_schema.py``

THREE AXES, ONE PER ACCUMULATOR OBJECT KIND. The calibrator deconvolves a library on the splice
graph, and the graph has three kinds of object, so the result has three axes::

    regions            N            deconvolved contained count + its geometric support
    contiguous boundaries E = N − refs deconvolved crossing count + its geometric support
    sj boundaries   J            the jumping flux -- certified RNA, never deconvolved

A contiguous boundary is a 0-bp boundary with ONE set of numbers, so there is no ``left``/``right``
pair and no ½ anywhere: the count arrays are ``count_{gdna,rna}_boundary`` and the per-boundary divisor
is ``crossing_eff_length``, carried here as ``gdna_boundary_eff_len``. ⛔ Do not reintroduce a
per-face split that a consumer then pools back — that sum-then-halve pattern is exactly what hides a
factor of 2.

``count_rna_spliced_boundary`` has no region twin, structurally: ``region_contained`` is credited only
when the fragment used no sj, so a region's contained population cannot hold a spliced molecule.

``__post_init__`` enforces the intrinsic invariants (per-axis shape, dtype, finiteness, sign); mass
conservation against the raw fragment counts is checked by the calibrator / tests, since it needs the
substrate the result does not carry.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import CalibrationConfig


def _check_axis_array(arr: np.ndarray, name: str, n: int) -> None:
    """Validate one per-object array: shape against ITS OWN axis, dtype, finite, non-negative.

    The dtype gate admits exact integers as well as float64. The accumulator's primary per-object
    observable is an integer **count** (no floats anywhere in the data model), so an integer array is a
    *better* input here than a float one — ``count_rna_sj`` is the flux verbatim and arrives
    integral. Everything else the gate checks is unchanged, and ``np.isfinite`` is well defined on
    integers.

    The shape check is the load-bearing one. ``E = N − n_refs`` differs from ``N`` by only a few
    hundred genome-wide, so an array keyed to the wrong axis is a *plausible* length and nothing
    downstream would fault on it — it would just silently read the wrong object's number.
    """
    if not isinstance(arr, np.ndarray):
        raise ValueError(f"CalibrationResult.{name} must be a numpy array.")
    if arr.dtype != np.float64 and arr.dtype.kind not in ("i", "u"):
        raise ValueError(
            f"CalibrationResult.{name} must be float64 or an integer count; got {arr.dtype}."
        )
    if arr.shape != (n,):
        raise ValueError(f"CalibrationResult.{name} has shape {arr.shape}; expected ({n},).")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"CalibrationResult.{name} contains non-finite values.")
    if np.any(arr < 0.0):
        raise ValueError(f"CalibrationResult.{name} must be non-negative.")


def _check_unit_interval(value: float, name: str, *, open_upper: bool = False) -> None:
    v = float(value)
    upper_ok = v < 1.0 if open_upper else v <= 1.0
    if not (0.0 <= v and upper_ok):
        bound = "[0, 1)" if open_upper else "[0, 1]"
        raise ValueError(f"CalibrationResult.{name} must be in {bound}; got {v}.")


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Deconvolved gDNA / RNA mass on the region and contiguous-boundary axes, the sj flux, the two
    gDNA geometric supports, and the library scalars."""

    # --- the deconvolved MIXTURE, per region (float64[n_regions]) ---
    #: ``chain_region_deconv``: the region's contained unspliced count split by the converged belief.
    count_gdna_region: np.ndarray
    count_rna_region: np.ndarray

    # --- the deconvolved MIXTURE, per contiguous boundary (float64[n_boundaries]) ---
    #: ``chain_boundary_deconv``: the boundary's unspliced crossing count split by the converged belief.
    #: ``count_rna_boundary`` is spliced-INCLUSIVE — a boundary's certified-RNA crossings are RNA whatever the
    #: unspliced mixture resolves to, since gDNA cannot be spliced — so per-boundary conservation
    #: ``count_gdna_boundary + count_rna_boundary == unspliced + spliced`` holds.
    count_gdna_boundary: np.ndarray
    count_rna_boundary: np.ndarray

    #: float64[n_boundaries] — the ``boundary_spliced`` part of ``count_rna_boundary``: molecules that crossed this
    #: boundary CONTIGUOUSLY having spliced somewhere else. Carried so ``assemble_priors`` can **withhold**
    #: it from ``rna_prior_count``: a spliced fragment has no gDNA candidate in the EM (gDNA does not
    #: splice), so it is guaranteed-RNA and assigned directly — counting it in the prior would double
    #: it and inflate the RNA side of the gDNA-vs-RNA *unspliced* split, which is the only thing the
    #: prior arbitrates. ``count_rna_boundary`` itself stays spliced-inclusive so conservation is preserved.
    count_rna_spliced_boundary: np.ndarray

    #: float64[n_boundaries] — THE INCIDENCE→FRAGMENT CONVERSION, per boundary. ``mass / count`` off the
    #: accumulator's conserved-mass bank: the mean fragment-mass ONE crossing at this boundary carries.
    #:
    #: ⛔ It is GEOMETRY, not a deconvolved mass, and the distinction is load-bearing. Every array
    #: above is a gDNA/RNA split that a perfect deconvolution would change; this one is a property of
    #: the partition and the fragment-length distribution and would be identical under any split. That
    #: is why it is NOT in ``prior_vs_oracle.OVERRIDE_FIELDS``: an oracle that overrode it would be
    #: answering a different question.
    #:
    #: ``assemble_priors`` multiplies each component's per-boundary mass by it, because the accumulator
    #: deposits ``+1`` on EVERY boundary a fragment crosses — ``max(K, 1)`` of them — so a sum over boundaries is
    #: an object-incidence count and the EM adds a FRAGMENT count. It is 1.0 where both flanking regions
    #: exceed every fragment length, and falls toward the region spacing where they do not.
    boundary_mass_per_crossing: np.ndarray

    # --- the JUMPING population, per sj boundary (float64[n_sj]) ---
    #: ⛔ AN INCIDENCE COUNT DESPITE THE NAME — it is not a mass, and the gap is large. This is
    #: ``sj_count`` summed over the genome-strand columns and nothing else, and a fragment deposits
    #: ``+1`` on EVERY sj it uses, so it runs about 2× the conserved mass and the over-count grows with
    #: how many sj a fragment spans. A consumer that wants the mass wants :attr:`sj_conserved_mass`,
    #: which is this array converted; reading this one instead is wrong in proportion to how SPLICED
    #: the object is, and for a per-transcript weight that is exactly the axis the answer varies over.
    #:
    #: Never deconvolved: a sj boundary is pure RNA by construction, so there is nothing to split. It
    #: is the third population at a boundary, and it is routinely orders of magnitude larger than
    #: ``count_rna_spliced_boundary`` at the same place: at a donor boundary the sj flux is the gene's
    #: whole spliced output while the spliced crossing is the handful of molecules that read through
    #: without splicing.
    #: ``assemble_priors`` does NOT consume it, and that is deliberate: sj fragments are certified RNA
    #: in exactly the sense ``count_rna_spliced_boundary`` is withheld for, so feeding them to
    #: ``rna_prior_count`` would load the RNA side of a split that arbitrates only unspliced fragments.
    #: It is exported for QC and reporting — the calibration's output should not be silent about the
    #: population that dominates a donor boundary.
    count_rna_sj: np.ndarray

    #: float64[n_boundaries] — the ``boundary_spliced`` twin of ``boundary_mass_per_crossing``, and
    #: float64[n_sj] — the sj one, ``sj_mass / sj_count``. Same kind of quantity as
    #: ``boundary_mass_per_crossing`` in every respect: GEOMETRY, identical under any split, and therefore
    #: not in ``prior_vs_oracle.OVERRIDE_FIELDS``.
    #:
    #: ⛔ The sj one is what makes a conserved LIBRARY fragment count computable at all: a spliced
    #: fragment whose every block lies inside one region crosses no boundary and is not contained, so
    #: without the sj bank it deposits nowhere conserved — a large minority of RNA fragments, and no
    #: gDNA ones, since gDNA cannot splice.
    boundary_spliced_mass_per_crossing: np.ndarray
    sj_mass_per_crossing: np.ndarray

    # --- the gDNA geometric supports: expected admissible START POSITIONS, per component ---
    #: float64[n_regions] — ``effective_length.contained_eff_length`` on the gDNA pmf,
    #: ``E_f[(region_len − w + 1)+]``. Under uniform genomic gDNA at density ρ the expected contained
    #: mass is EXACTLY ``ρ · gdna_region_eff_len`` (a fragment must FIT to be contained), so dividing by
    #: it recovers ρ — the bedrock "factor 1 under uniform gDNA" invariant. This, NOT the genomic region
    #: length, is the density-correct divisor: the raw length ignores the fit-inside constraint and so
    #: understates a short region's density, manufacturing a spurious contraction in an unenriched library.
    gdna_region_eff_len: np.ndarray
    #: float64[n_boundaries] — ``effective_length.crossing_eff_length`` on the gDNA pmf. Uniform across
    #: boundaries, and that is physics rather than a placeholder: gDNA's template is the chromosome,
    #: so it takes ``UNBOUNDED_REACH`` on both sides at every boundary and the divisor collapses to
    #: ``mu_g − 1``. It stays a per-boundary array because that is the axis its consumers index it on.
    gdna_boundary_eff_len: np.ndarray

    # --- the RNA geometric supports: the SAME two frames, on the RNA pmf ---
    #: float64[n_regions] / float64[n_boundaries] — ``contained_eff_length`` and ``crossing_eff_length`` on the
    #: RNA pmf: the RNA population's OWN opportunity at each object.
    #:
    #: ⛔ NO CONSUMER IN ``src/`` TODAY. The prior is a conserved FRAGMENT COUNT and divides by nothing
    #: on the mass path (``tests/calibration/test_prior_units.py``), so nothing reads an RNA density
    #: divisor.
    #:
    #: They are ORPHANED, not dead, and deleting them would be the wrong repair. Their gDNA twins are
    #: load-bearing in four modules (``capture_eff_length``, ``track``, ``derive``, ``priors``); the
    #: RNA pair is the same quantity for the other population, and any prior that reasons about RNA
    #: DENSITY — per transcript or per object — needs exactly this divisor.
    #:
    #: Like their gDNA twins they are PROJECTED off ``RegionGeometry.eff_rna`` by ``_project_eff``,
    #: never recomputed, so they are byte-identically the opportunity the SOLVER used. That is the
    #: property that makes them safe to reason with.
    rna_region_eff_len: np.ndarray
    rna_boundary_eff_len: np.ndarray

    # --- THE THREE-WAY COMPOSITION — the simplex ψ actually solves, per object ---
    #: float64[n_regions] / float64[n_boundaries] — the solved ``(f_g, f_pos, f_neg)`` at each object: the
    #: gDNA share and the two RNA STRAND shares of that object's unspliced population.
    #:
    #: This is AXIOM 0's ``T(slot)``, published. The solve is over
    #: ``{gDNA} ∪ {RNA+ if free_pos} ∪ {RNA− if free_neg}`` at every slot, so the answer is three
    #: numbers; ``count_gdna_*`` and ``count_rna_*`` are that answer with the two RNA strands summed.
    #:
    #: The three close by construction: ψ solves ``(f_g, w_pos)`` and `simplex_logodds._compose` maps
    #: them onto the simplex, while a slot the solve does not reach keeps its signature-binary belief.
    #: They are published as solved, never renormalised.
    gdna_frac_region: np.ndarray
    rna_pos_frac_region: np.ndarray
    rna_neg_frac_region: np.ndarray
    gdna_frac_boundary: np.ndarray
    rna_pos_frac_boundary: np.ndarray
    rna_neg_frac_boundary: np.ndarray

    # --- library scalars ---
    gdna_density_global: float  # >= 0, global gDNA density (mass/bp); 0 in a zero-gDNA library
    #: The fully-captured gDNA density the ruler and the locus gDNA effective length contract against:
    #: the located enriched mode of the last refit's fitted landscape (`abundance_landscape.located_enriched_mode`),
    #: or ``None`` — no enriched gDNA mode, which is every capture-OFF and every gDNA-free library, and
    #: then nothing contracts. A positive finite density when present.
    gdna_reference_density: float | None
    #: The regime behind the reference: how many LOCATED kernels (`landscape.DensityLandscape.located`)
    #: the enriched mode rests on — the located population's own resolution is ``√n``, so this is the
    #: number a reader compares against; ``0`` exactly when the reference is ``None``.
    gdna_reference_members: int
    #: float64[n_regions] — each piece's capture efficiency ``E[min(ρ/ρ_ref, 1)]``, the posterior mean
    #: of its clipped gDNA density against the reference under the fitted landscape, from its own
    #: contained count and the crossings at every boundary within a fragment's reach
    #: (`capture_efficiency.capture_efficiencies`); exactly 1 everywhere when the reference is ``None``.
    #: The ruler and the locus prior read this and re-derive nothing (`capture_eff_length`, `priors`).
    gdna_capture_efficiency_region: np.ndarray
    #: float64[n_boundaries] — each boundary's own capture efficiency, the posterior mean of its clipped
    #: gDNA density from its crossing count on its crossing support; exactly 1 everywhere when the
    #: reference is ``None``. The locus prior reads it beside its count (`priors`): the count is the
    #: calibration's masses on regions and boundaries, and the length is those objects at their
    #: efficiencies. The transcript ruler never reads it — no boundary object enters a length over bases.
    gdna_capture_efficiency_boundary: np.ndarray
    rna_sense_frac: float  # in [0, 1], RNA sense fraction used by the strand clue
    gdna_strand_overdispersion: float  # in [0, 1), fitted gDNA strand Beta-Binomial dispersion
    rna_strand_overdispersion: float  # in [0, 1), fitted RNA strand Beta-Binomial dispersion

    # --- provenance: the three axis lengths, independent of each other ---
    n_regions: int
    n_boundaries: int
    n_sj: int
    config: CalibrationConfig

    def __post_init__(self) -> None:
        for axis in ("n_regions", "n_boundaries", "n_sj"):
            if int(getattr(self, axis)) < 0:
                raise ValueError(
                    f"CalibrationResult.{axis} must be >= 0; got {getattr(self, axis)}."
                )

        for name in (
            "count_gdna_region",
            "count_rna_region",
            "gdna_region_eff_len",
            "rna_region_eff_len",
            "gdna_frac_region",
            "rna_pos_frac_region",
            "rna_neg_frac_region",
        ):
            _check_axis_array(getattr(self, name), name, self.n_regions)
        for name in (
            "count_gdna_boundary",
            "count_rna_boundary",
            "count_rna_spliced_boundary",
            "boundary_mass_per_crossing",
            "boundary_spliced_mass_per_crossing",
            "gdna_boundary_eff_len",
            "rna_boundary_eff_len",
            "gdna_frac_boundary",
            "rna_pos_frac_boundary",
            "rna_neg_frac_boundary",
        ):
            _check_axis_array(getattr(self, name), name, self.n_boundaries)

        # Each component is a FRACTION, so it is bounded by 1. Closure holds by construction
        # (`simplex_logodds._compose`) and is not re-asserted here.
        for name in (
            "gdna_frac_region",
            "rna_pos_frac_region",
            "rna_neg_frac_region",
            "gdna_frac_boundary",
            "rna_pos_frac_boundary",
            "rna_neg_frac_boundary",
        ):
            arr = np.asarray(getattr(self, name), dtype=np.float64)
            if np.any(arr > 1.0 + 1e-9):
                i = int(np.flatnonzero(arr > 1.0 + 1e-9)[0])
                raise ValueError(
                    f"CalibrationResult.{name} is a fraction and must not exceed 1; index {i} is "
                    f"{float(arr[i])!r}."
                )
        for name in ("count_rna_sj", "sj_mass_per_crossing"):
            _check_axis_array(getattr(self, name), name, self.n_sj)
        _check_axis_array(
            self.gdna_capture_efficiency_region, "gdna_capture_efficiency_region", self.n_regions
        )
        c = np.asarray(self.gdna_capture_efficiency_region, dtype=np.float64)
        if np.any(c < 0.0) or np.any(c > 1.0 + 1e-9):
            raise ValueError("CalibrationResult.gdna_capture_efficiency_region must lie in [0, 1].")
        _check_axis_array(
            self.gdna_capture_efficiency_boundary,
            "gdna_capture_efficiency_boundary",
            self.n_boundaries,
        )
        cb = np.asarray(self.gdna_capture_efficiency_boundary, dtype=np.float64)
        if np.any(cb < 0.0) or np.any(cb > 1.0 + 1e-9):
            raise ValueError(
                "CalibrationResult.gdna_capture_efficiency_boundary must lie in [0, 1]."
            )
        if self.gdna_reference_density is None and (np.any(c != 1.0) or np.any(cb != 1.0)):
            raise ValueError(
                "CalibrationResult.gdna_capture_efficiency_region and _boundary must be exactly 1 "
                "everywhere when there is no reference: nothing is depleted relative to anything."
            )

        if self.gdna_reference_density is not None and not (
            np.isfinite(self.gdna_reference_density) and self.gdna_reference_density > 0.0
        ):
            raise ValueError(
                "CalibrationResult.gdna_reference_density must be None or finite and > 0; "
                f"got {self.gdna_reference_density}."
            )
        if (int(self.gdna_reference_members) > 0) != (self.gdna_reference_density is not None):
            raise ValueError(
                "CalibrationResult.gdna_reference_members must be > 0 exactly when a reference is "
                f"present; got {self.gdna_reference_members} with reference {self.gdna_reference_density}."
            )
        if not np.isfinite(self.gdna_density_global) or self.gdna_density_global < 0.0:
            raise ValueError(
                "CalibrationResult.gdna_density_global must be finite and >= 0; "
                f"got {self.gdna_density_global}."
            )
        _check_unit_interval(self.rna_sense_frac, "rna_sense_frac")
        _check_unit_interval(
            self.gdna_strand_overdispersion, "gdna_strand_overdispersion", open_upper=True
        )
        _check_unit_interval(
            self.rna_strand_overdispersion, "rna_strand_overdispersion", open_upper=True
        )

    # ── THE CONSERVED SJ MASS — the third axis in FRAGMENT units, in ONE place ──

    @property
    def sj_conserved_mass(self) -> np.ndarray:
        """float64[n_sj] — the CONSERVED fragment mass at each sj. Sums to one per spliced fragment
        across the sj it used, where :attr:`count_rna_sj` is ``+1`` on each of them.

        This is the accumulator's ``sj_mass`` bank recovered exactly: ``sj_mass_per_crossing`` is
        ``sj_mass / sj_count`` and ``count_rna_sj`` is ``sj_count``, so the product is ``sj_mass``
        identically, up to float division-then-multiplication.

        ⛔ At a sj nothing crossed the answer is 0, and the multiplication is what produces it — do
        not branch. ``mass_per_crossing`` deliberately keeps the ``1.0`` identity where the count is
        zero (there is no mass at an unobserved boundary to rescale), so a ``where(count > 0, …)``
        fallback to the factor would publish phantom mass on the large zero-count share of this axis,
        the one axis that is certified RNA by construction.

        A PROPERTY, never a stored field, and that is forced rather than preferred: ``count_rna_sj``
        is in ``prior_vs_oracle.OVERRIDE_FIELDS``, an oracle arm swaps it with
        ``dataclasses.replace``, and a cached array would survive that swap and silently describe the
        array it replaced — the same reason :attr:`library_rna_fragments` is derived.
        """
        return np.asarray(self.count_rna_sj, dtype=np.float64) * np.asarray(
            self.sj_mass_per_crossing, dtype=np.float64
        )

    # ── THE LIBRARY FRAGMENT COUNT — the deliverable, in FRAGMENT units, derived in ONE place ──

    @property
    def library_gdna_fragments(self) -> float:
        """gDNA fragments in the library — a CONSERVED count, not an object-incidence sum.

        gDNA appears on TWO axes only, because it cannot splice, and containment is exclusive — so a
        region term is already a fragment count and only the crossing term needs converting.
        """
        return float(
            np.asarray(self.count_gdna_region, dtype=np.float64).sum()
            + (
                np.asarray(self.count_gdna_boundary, dtype=np.float64)
                * np.asarray(self.boundary_mass_per_crossing, dtype=np.float64)
            ).sum()
        )

    @property
    def library_rna_fragments(self) -> float:
        """RNA fragments in the library — all THREE axes, each converted by its OWN population's ``q``.

        ⛔ Never sum the raw incidences instead. One fragment books ``max(K,1)`` boundary crossings AND
        one incidence per sj it uses, so an incidence sum over-counts, and the gDNA:RNA ratio it
        implies is biased wherever the two components' inflations differ — which they do, because gDNA
        cannot splice. These conserved counts reproduce the true origin split exactly.

        The spliced crossings and the sj flux are the SAME fragments split across two banks by
        the deposit rule — ``boundary_spliced_mass`` holds the share of a spliced fragment's bases in blocks
        that crossed a boundary and ``sj_mass`` the share in blocks that crossed none — and the two sum to
        exactly one per fragment. Adding both is conservation, not double counting.

        A PROPERTY, never a stored field. ``prior_vs_oracle`` swaps the deconvolved arrays for truth with
        ``dataclasses.replace``; a cached scalar would survive that swap and silently describe the old
        arrays. Deriving it means the
        oracle arm's count is the oracle's by construction.
        """
        unspliced_boundary = np.maximum(
            np.asarray(self.count_rna_boundary, dtype=np.float64)
            - np.asarray(self.count_rna_spliced_boundary, dtype=np.float64),
            0.0,
        )
        return float(
            np.asarray(self.count_rna_region, dtype=np.float64).sum()
            + (
                unspliced_boundary * np.asarray(self.boundary_mass_per_crossing, dtype=np.float64)
            ).sum()
            + (
                np.asarray(self.count_rna_spliced_boundary, dtype=np.float64)
                * np.asarray(self.boundary_spliced_mass_per_crossing, dtype=np.float64)
            ).sum()
            # ONE home for the sj conversion: spelling the product out here as well lets a caller
            # reading the property and a caller reading this sum disagree.
            + self.sj_conserved_mass.sum()
        )


__all__ = ["CalibrationResult"]
