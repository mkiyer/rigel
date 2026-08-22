"""CalibrationResult — the calibrator's output schema.

    Gate: ``tests/calibration/test_result_schema.py``

⭐ **THREE AXES, ONE PER ACCUMULATOR OBJECT KIND.** The calibrator deconvolves a library on the splice
graph, and the graph has three kinds of object, so the result has three axes::

    regions            N            deconvolved contained mass + its geometric support
    contiguous boundaries E = N − refs deconvolved crossing mass + its geometric support
    sj boundaries   J            the jumping flux -- certified RNA, never deconvolved

⛔ **THE ``left``/``right`` PAIR IS GONE, AND SO IS THE ½ IT CARRIED.** The predecessor carried SIX
per-region mass arrays — ``mass_{gdna,rna}_{contained,left,right}`` — because a boundary had two sides
sitting in differently-sized flanks. ``priors.assemble_priors`` then pooled two of them straight back
as ``mass_gdna_right[r] + mass_gdna_left[r+1]``, and ``capture_eff_length._pooled_boundary_arrays`` did the
identical thing. That split-then-re-pool was a no-op with a history: the same sum-then-halve pattern hid
an exact factor of 2 for months. A contiguous boundary is a 0-bp boundary with
ONE set of numbers, so the pair collapses to ``mass_{gdna,rna}_boundary`` and the pooling **disappears
rather than being re-derived** (owner ruling).

⛔ **``gdna_boundary_len`` HAS NO SUCCESSOR.** It was ``boundary_side_eff_length = E[min(ℓ,L)]/2``, a
per-FACE divisor whose ½ existed only because the face's mass was half a crossing. S5.c deleted it. The
per-boundary divisor is ``crossing_eff_length``, one number, no ½ — carried here as ``gdna_boundary_eff_len``.
Any comment claiming ``gdna_boundary_len`` "IS the halved per-side density length" describes a quantity
that no longer exists.

⚠ **``mass_rna_spliced_boundary`` has no region twin, structurally.** ``region_contained`` is credited only when
the fragment used no sj, so a region's contained population cannot hold a spliced molecule.

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

    ⚠ The dtype gate admits exact integers as well as float64. The accumulator's primary per-object
    observable is an integer **count** (no floats anywhere in the data model), so an integer array is a
    *better* input here than a float one — ``count_rna_sj`` is the flux verbatim and arrives
    integral. Everything else the gate checks is unchanged, and ``np.isfinite`` is well defined on
    integers.

    ⚠ The shape check is the load-bearing one. ``E = N − n_refs`` differs from ``N`` by only a few
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
    mass_gdna_region: np.ndarray
    mass_rna_region: np.ndarray

    # --- the deconvolved MIXTURE, per contiguous boundary (float64[n_boundaries]) ---
    #: ``chain_boundary_deconv``: the boundary's unspliced crossing count split by the converged belief.
    #: ``mass_rna_boundary`` is spliced-INCLUSIVE — an boundary's certified-RNA crossings are RNA whatever the
    #: unspliced mixture resolves to, since gDNA cannot be spliced — so per-boundary conservation
    #: ``mass_gdna_boundary + mass_rna_boundary == unspliced + spliced`` holds.
    mass_gdna_boundary: np.ndarray
    mass_rna_boundary: np.ndarray

    #: float64[n_boundaries] — the ``boundary_spliced`` part of ``mass_rna_boundary``: molecules that crossed this
    #: boundary CONTIGUOUSLY having spliced somewhere else. Carried so ``assemble_priors`` can **withhold**
    #: it from ``rna_prior_count``: a spliced fragment has no gDNA candidate in the EM (gDNA does not
    #: splice), so it is guaranteed-RNA and assigned directly — counting it in the prior would double
    #: it and inflate the RNA side of the gDNA-vs-RNA *unspliced* split, which is the only thing the
    #: prior arbitrates. ``mass_rna_boundary`` itself stays spliced-inclusive so conservation is preserved.
    mass_rna_spliced_boundary: np.ndarray

    #: float64[n_boundaries] — ⭐⭐ **THE INCIDENCE→FRAGMENT CONVERSION, per boundary.** ``mass / count`` off the
    #: accumulator's conserved-mass bank: the mean fragment-mass ONE crossing at this boundary carries.
    #:
    #: ⛔ **It is GEOMETRY, not a deconvolved mass, and the distinction is load-bearing.** Every array
    #: above is a gDNA/RNA split that a perfect deconvolution would change; this one is a property of
    #: the partition and the fragment-length distribution and would be identical under any split. That
    #: is why it is NOT in ``prior_vs_oracle.OVERRIDE_FIELDS``: an oracle that overrode it would be
    #: answering a different question.
    #:
    #: ⭐ ``assemble_priors`` multiplies each component's per-boundary mass by it, because the accumulator
    #: deposits ``+1`` on EVERY boundary a fragment crosses — ``max(K, 1)`` of them — so a sum over boundaries is
    #: an object-incidence count and the EM adds a FRAGMENT count. It is 1.0 where both flanking regions
    #: exceed every fragment length, and falls toward the region spacing where they do not.
    boundary_mass_per_crossing: np.ndarray

    # --- the JUMPING population, per sj boundary (float64[n_sj]) ---
    #: ⛔⛔ **AN INCIDENCE COUNT DESPITE THE NAME — IT IS NOT A MASS, AND THE GAP IS LARGE.** This is
    #: ``sj_count`` summed over the genome-strand columns and nothing else, and a fragment deposits
    #: ``+1`` on EVERY sj it uses. Measured on ``g00 ss0.99 capture_off``: **2.0719 incidences per
    #: unit of conserved mass** (5,668,526 against 2,735,958.8), and the over-count grows with how many
    #: sj a fragment spans. ⭐ **A consumer that wants the mass wants
    #: :attr:`sj_conserved_mass`**, which is this array converted; reading this one instead is
    #: wrong in proportion to how SPLICED the object is, and for a per-transcript weight that is exactly
    #: the axis the answer varies over. ⚠ The two remaining ``mass_*_boundary`` incidences are due the
    #: same INCIDENCE-vs-MASS rename; that is its own commit.
    #:
    #: ⭐ **Never deconvolved: a sj boundary is pure mature RNA by construction**, so there is nothing
    #: to split. It is the third population at a boundary, and it is routinely two orders of magnitude
    #: larger than ``mass_rna_spliced_boundary`` at the same place: at a donor boundary the sj flux is the
    #: gene's whole mature output while the spliced crossing is the handful of molecules that read
    #: through without splicing.
    #: ⚠ **``assemble_priors`` does NOT consume it, and that is deliberate.** SpliceJunction fragments are
    #: certified RNA in exactly the sense ``mass_rna_spliced_boundary`` is withheld for, so feeding them to
    #: ``rna_prior_count`` would load the RNA side of a split that arbitrates only unspliced fragments.
    #: It is exported for QC and reporting — the calibration's output should not be silent about the
    #: population that dominates a donor boundary (owner ruling, 2026-07-30).
    count_rna_sj: np.ndarray

    #: float64[n_boundaries] — the ``boundary_spliced`` twin of ``boundary_mass_per_crossing``, and
    #: float64[n_sj] — the sj one, ``sj_mass / sj_count``. Same kind of quantity as
    #: ``boundary_mass_per_crossing`` in every respect: GEOMETRY, identical under any split, and therefore
    #: not in ``prior_vs_oracle.OVERRIDE_FIELDS``.
    #:
    #: ⛔ **The sj one did not exist until ``sj_mass`` did, and that is why a conserved LIBRARY
    #: fragment count was not computable.** A spliced fragment whose every block lies inside one region
    #: crosses no boundary and is not contained, so it deposited on no conserved bank — **1,222,375 of
    #: 4,830,713 RNA fragments (25.3 %)** on ladder g50 capture_off, against 0 of 4,997,761 gDNA
    #: fragments, since gDNA cannot splice.
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
    #: float64[n_boundaries] — ``effective_length.crossing_eff_length`` on the gDNA pmf. ⚠ **Uniform across
    #: boundaries today, and that is physics rather than a placeholder**: gDNA's template is the chromosome,
    #: so it takes ``UNBOUNDED_REACH`` on both sides at every boundary and the divisor collapses to
    #: ``mu_g − 1``. It stays a per-boundary array because that is the axis its consumers index it on.
    gdna_boundary_eff_len: np.ndarray

    # --- the RNA geometric supports: the SAME two frames, on the RNA pmf ---
    #: float64[n_regions] / float64[n_boundaries] — ``contained_eff_length`` and ``crossing_eff_length`` on the
    #: RNA pmf: the RNA population's OWN opportunity at each object.
    #:
    #: ⛔⛔ **NO CONSUMER IN ``src/`` TODAY, AND THE REASON IS A DESIGN CHANGE, NOT AN OVERSIGHT.** This
    #: docstring used to open "they exist because a mass is not a fragment count" and state that
    #: ``assemble_priors`` divides each component's mass by its own opportunity to get a density. **That
    #: consumer was removed.** The prior was rewritten as a conserved FRAGMENT COUNT at ``d045d820``
    #: (2026-08-01) and ``5591cc01`` (2026-08-08); it now divides by nothing on the mass path —
    #: ``tests/calibration/test_prior_units.py``: *"The prior no longer divides by anything."* The
    #: quantified g:r tilt the old text cited (gDNA 1.031 incidences per fragment against RNA ~1.17) was
    #: an argument for a density prior that no longer ships.
    #:
    #: ⚠ **They are ORPHANED, not dead, and deleting them would be the wrong repair.** Their gDNA twins
    #: are load-bearing in four modules (``capture_eff_length``, ``track``, ``derive``, ``priors``); the
    #: RNA pair is the same quantity for the other population, and any future prior that reasons about
    #: RNA DENSITY — per transcript or per object — needs exactly this divisor.
    #:
    #: ⭐ Like their gDNA twins they are PROJECTED off ``RegionGeometry.eff_rna`` by ``_project_eff``,
    #: never recomputed, so they are byte-identically the opportunity the SOLVER used. That is the
    #: property that makes them safe to reason with.
    rna_region_eff_len: np.ndarray
    rna_boundary_eff_len: np.ndarray

    # --- ⭐⭐ THE THREE-WAY COMPOSITION — the simplex ψ actually solves, per object ---
    #: float64[n_regions] / float64[n_boundaries] — the solved ``(f_g, f_pos, f_neg)`` at each object: the
    #: gDNA share and the two RNA STRAND shares of that object's unspliced population.
    #:
    #: ⭐⭐ **THIS IS AXIOM 0's ``T(slot)``, PUBLISHED.** The solve is over
    #: ``{gDNA} ∪ {RNA+ if free_pos} ∪ {RNA− if free_neg}`` at every slot, so the answer has always been
    #: three numbers; ``mass_gdna_*`` and ``mass_rna_*`` are that answer with the two RNA strands summed.
    #: A consumer asking "which STRAND's RNA is here" had to re-derive it or go without.
    #:
    #: ⛔ **The crossing axis published ``0`` for both RNA strands until 2026-08-12** — ``chain_boundary_deconv``
    #: projected ``f_g`` and emitted ``np.zeros(n)`` for ``f_pos``/``f_neg``, so the boundary composition
    #: summed to ``f_g`` alone rather than to 1. Nothing consumed it, which is why it survived.
    #:
    #: ⛔⛔ **THE THREE DO NOT SUM TO 1 ON ABOUT A QUARTER OF EVERY AXIS, AND THAT IS A MEASURED DEFECT
    #: IN ψ RATHER THAN A PROPERTY OF THIS PROJECTION.** ``RegionDeconv`` asserts
    #: *"posterior means; f_pos+f_neg+gdna_frac = 1"*. Measured on ``g00 ss0.99 capture_off``
    #: (2026-08-12), with every object addressed by a chain slot:
    #:
    #: ====== ============== ================================== =========
    #: axis   sums to 1      sums into (0, 1)                   sums > 1
    #: ====== ============== ================================== =========
    #: REGION   74.72 %        **25.25 %** — median 0.978, p5 0.869   12
    #: BOUNDARY   77.24 %        **22.71 %** — median 0.979, p5 0.850   16
    #: ====== ============== ================================== =========
    #:
    #: The mechanism is visible at ``sweep.py``'s write-back: the three posterior means are
    #: ``np.clip(·, 0, 1)``-ed **INDEPENDENTLY**, and an unsolvable slot keeps an init instead. Clipping
    #: a simplex component-wise does not preserve the sum. ⚠ By linearity of expectation three posterior
    #: means over one lattice *should* close, so the deficit is not inherent to taking means.
    #:
    #: ⛔ **Published AS SOLVED and deliberately NOT renormalised.** Dividing by the sum would make the
    #: arrays look like a composition while hiding how far ψ's answer is from being one, and a consumer
    #: cannot then tell a solved object from a 15 %-short one. The schema gate therefore bounds each
    #: component to ``[0, 1]`` — which is true — and does not assert closure, which is not.
    gdna_frac_region: np.ndarray
    rna_pos_frac_region: np.ndarray
    rna_neg_frac_region: np.ndarray
    gdna_frac_boundary: np.ndarray
    rna_pos_frac_boundary: np.ndarray
    rna_neg_frac_boundary: np.ndarray

    # --- library scalars ---
    gdna_density_global: float  # >= 0, global gDNA density (mass/bp); 0 in a zero-gDNA library
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
            "mass_gdna_region",
            "mass_rna_region",
            "gdna_region_eff_len",
            "rna_region_eff_len",
            "gdna_frac_region",
            "rna_pos_frac_region",
            "rna_neg_frac_region",
        ):
            _check_axis_array(getattr(self, name), name, self.n_regions)
        for name in (
            "mass_gdna_boundary",
            "mass_rna_boundary",
            "mass_rna_spliced_boundary",
            "boundary_mass_per_crossing",
            "boundary_spliced_mass_per_crossing",
            "gdna_boundary_eff_len",
            "rna_boundary_eff_len",
            "gdna_frac_boundary",
            "rna_pos_frac_boundary",
            "rna_neg_frac_boundary",
        ):
            _check_axis_array(getattr(self, name), name, self.n_boundaries)

        # ⭐ Each component is a FRACTION, so it is bounded by 1 — the one thing that is true of all
        # three today. ⛔ Closure (`f_g + f_pos + f_neg == 1`) is deliberately NOT asserted: it fails on
        # ~25 % of both axes, for the reason the field docstring records, and a gate that fails on the
        # shipped configuration is a gate nobody can keep green. Assert it the day ψ's write-back stops
        # clipping the three independently.
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

    # ── ⭐⭐ THE CONSERVED SJ MASS — the third axis in FRAGMENT units, in ONE place ──

    @property
    def sj_conserved_mass(self) -> np.ndarray:
        """float64[n_sj] — the CONSERVED fragment mass at each sj. **Sums to one per
        spliced fragment across the sj it used**, where :attr:`count_rna_sj` is ``+1`` on
        each of them.

        ⭐⭐ **THIS IS THE ACCUMULATOR'S ``sj_mass`` BANK, RECOVERED EXACTLY.**
        ``sj_mass_per_crossing`` is ``sj_mass / sj_count`` and ``count_rna_sj`` is
        ``sj_count``, so the product is ``sj_mass`` identically. Measured on ``g00 ss0.99 capture_off``
        (13,482 sj): **12,758 elements bit-identical** to the bank and a worst element of
        **9.1e-13**, which is float division-then-multiplication and not a modelling gap.

        ⭐ **At a sj nothing crossed it is 0, and that needs the multiplication rather than a
        branch.** ``mass_per_crossing`` deliberately keeps the ``1.0`` identity where the count is zero
        — there is no mass at an unobserved boundary to rescale — and multiplying by the zero incidence is
        what turns it back into the 0 that is correct here. **4,636 of those 13,482 sj are
        zero-count**, so a ``where(count > 0, …)`` fallback to the factor would publish phantom mass on
        a third of the axis, on the one axis that is certified RNA by construction.

        ⛔ **It exists because the derivation was reachable only by knowing to perform it**, while a
        field named ``count_rna_sj`` sat next to it looking like the answer and reading **2.0719×**
        high (5,668,526 incidences against 2,735,958.8 units of mass). Every other axis publishes its
        own conversion; this one made the caller do it.

        ⚠ **A PROPERTY, never a stored field, and that is forced rather than preferred.**
        ``count_rna_sj`` is in ``prior_vs_oracle.OVERRIDE_FIELDS``: an oracle arm swaps it with
        ``dataclasses.replace``, and a cached array would survive that swap and silently describe the
        array it replaced — ``TRAPS: a-hash-that-misses-its-artifact`` in dataclass form, the same
        reason :attr:`library_rna_fragments` is derived.
        """
        return np.asarray(self.count_rna_sj, dtype=np.float64) * np.asarray(
            self.sj_mass_per_crossing, dtype=np.float64
        )

    # ── ⭐⭐⭐ THE LIBRARY FRAGMENT COUNT — the deliverable, in FRAGMENT units, derived in ONE place ──

    @property
    def library_gdna_fragments(self) -> float:
        """gDNA fragments in the library — a CONSERVED count, not an object-incidence sum.

        gDNA appears on TWO axes only, because it cannot splice, and containment is exclusive — so a
        region term is already a fragment count and only the crossing term needs converting. That is why
        this side was exact all along while the RNA side was 25.3 % short.
        """
        return float(
            np.asarray(self.mass_gdna_region, dtype=np.float64).sum()
            + (
                np.asarray(self.mass_gdna_boundary, dtype=np.float64)
                * np.asarray(self.boundary_mass_per_crossing, dtype=np.float64)
            ).sum()
        )

    @property
    def library_rna_fragments(self) -> float:
        """RNA fragments in the library — all THREE axes, each converted by its OWN population's ``q``.

        ⛔⛔ **THE TREE ONCE HELD THREE DISAGREEING ANSWERS TO THIS**, and every one summed INCIDENCES:
        ``pipeline.py`` added region + boundary + sj raw, ``calibration_truth_ab.py`` added region + boundary
        and dropped the sj, and only ``assemble_priors`` converted anything. One fragment books
        ``max(K,1)`` boundary crossings AND one incidence per sj it uses, so the sum over-counts and
        the ratio is biased wherever the two components' inflations differ — which they do. Measured on
        the origin-split oracle, ladder g50 capture_off: the incidence ratio reads ``f_gdna`` **0.3851**
        against a truth of **0.5085**, while these conserved counts reproduce **0.5085** exactly, each
        origin landing at **1.000x** its deposited fragments with **0** unaccounted.

        ⭐ The spliced crossings and the sj flux are the SAME fragments split across two banks by
        the deposit rule — ``boundary_spliced_mass`` holds the share of a spliced fragment's bases in blocks
        that crossed a boundary and ``sj_mass`` the share in blocks that crossed none — and the two sum to
        exactly one per fragment. Adding both is conservation, not double counting.

        ⚠ **A PROPERTY, never a stored field.** ``prior_vs_oracle`` swaps the mass arrays for truth with
        ``dataclasses.replace``; a cached scalar would survive that swap and silently describe the old
        arrays (``TRAPS: a-hash-that-misses-its-artifact``, in dataclass form). Deriving it means the
        oracle arm's count is the oracle's by construction.
        """
        unspliced_boundary = np.maximum(
            np.asarray(self.mass_rna_boundary, dtype=np.float64)
            - np.asarray(self.mass_rna_spliced_boundary, dtype=np.float64),
            0.0,
        )
        return float(
            np.asarray(self.mass_rna_region, dtype=np.float64).sum()
            + (
                unspliced_boundary * np.asarray(self.boundary_mass_per_crossing, dtype=np.float64)
            ).sum()
            + (
                np.asarray(self.mass_rna_spliced_boundary, dtype=np.float64)
                * np.asarray(self.boundary_spliced_mass_per_crossing, dtype=np.float64)
            ).sum()
            # ⭐ ONE home for the sj conversion — this term used to spell the product out a second
            # time, so a caller reading the property and a caller reading this sum could disagree.
            + self.sj_conserved_mass.sum()
        )


__all__ = ["CalibrationResult"]
