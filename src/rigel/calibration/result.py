"""CalibrationResult — the calibrator's output schema.

    Gate: ``tests/calibration/test_result_schema.py``

⭐ **THREE AXES, ONE PER ACCUMULATOR OBJECT KIND.** The calibrator deconvolves a library on the splice
graph, and the graph has three kinds of object, so the result has three axes::

    nodes            N            deconvolved contained mass + its geometric support
    contiguous edges E = N − refs deconvolved crossing mass + its geometric support
    junction edges   J            the jumping flux -- certified RNA, never deconvolved

⛔ **THE ``left``/``right`` PAIR IS GONE, AND SO IS THE ½ IT CARRIED.** The predecessor carried SIX
per-region mass arrays — ``mass_{gdna,rna}_{contained,left,right}`` — because a boundary had two sides
sitting in differently-sized flanks. ``priors.assemble_priors`` then pooled two of them straight back
as ``mass_gdna_right[r] + mass_gdna_left[r+1]``, and ``capture_eff_length._pooled_seam_arrays`` did the
identical thing. That split-then-re-pool was a no-op with a history: the same sum-then-halve pattern hid
an exact factor of 2 for months. A contiguous edge is a 0-bp line with
ONE set of numbers, so the pair collapses to ``mass_{gdna,rna}_edge`` and the pooling **disappears
rather than being re-derived** (owner ruling).

⛔ **``gdna_boundary_len`` HAS NO SUCCESSOR.** It was ``boundary_side_eff_length = E[min(ℓ,L)]/2``, a
per-FACE divisor whose ½ existed only because the face's mass was half a crossing. S5.c deleted it. The
per-edge divisor is ``crossing_eff_length``, one number, no ½ — carried here as ``gdna_edge_eff_len``.
Any comment claiming ``gdna_boundary_len`` "IS the halved per-side density length" describes a quantity
that no longer exists.

⚠ **``mass_rna_spliced_edge`` has no node twin, structurally.** ``node_contained`` is credited only when
the fragment used no junction, so a node's contained population cannot hold a spliced molecule.

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
    *better* input here than a float one — ``mass_rna_junction`` is the flux verbatim and arrives
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
    """Deconvolved gDNA / RNA mass on the node and contiguous-edge axes, the junction flux, the two
    gDNA geometric supports, and the library scalars."""

    # --- the deconvolved MIXTURE, per node (float64[n_nodes]) ---
    #: ``chain_node_deconv``: the node's contained unspliced count split by the converged belief.
    mass_gdna_node: np.ndarray
    mass_rna_node: np.ndarray

    # --- the deconvolved MIXTURE, per contiguous edge (float64[n_edges]) ---
    #: ``chain_edge_deconv``: the line's unspliced crossing count split by the converged belief.
    #: ``mass_rna_edge`` is spliced-INCLUSIVE — an edge's certified-RNA crossings are RNA whatever the
    #: unspliced mixture resolves to, since gDNA cannot be spliced — so per-edge conservation
    #: ``mass_gdna_edge + mass_rna_edge == unspliced + spliced`` holds.
    mass_gdna_edge: np.ndarray
    mass_rna_edge: np.ndarray

    #: float64[n_edges] — the ``edge_spliced`` part of ``mass_rna_edge``: molecules that crossed this
    #: line CONTIGUOUSLY having spliced somewhere else. Carried so ``assemble_priors`` can **withhold**
    #: it from ``rna_prior_count``: a spliced fragment has no gDNA candidate in the EM (gDNA does not
    #: splice), so it is guaranteed-RNA and assigned directly — counting it in the prior would double
    #: it and inflate the RNA side of the gDNA-vs-RNA *unspliced* split, which is the only thing the
    #: prior arbitrates. ``mass_rna_edge`` itself stays spliced-inclusive so conservation is preserved.
    mass_rna_spliced_edge: np.ndarray

    #: float64[n_edges] — ⭐⭐ **THE INCIDENCE→FRAGMENT CONVERSION, per line.** ``mass / count`` off the
    #: accumulator's conserved-mass bank: the mean fragment-mass ONE crossing at this line carries.
    #:
    #: ⛔ **It is GEOMETRY, not a deconvolved mass, and the distinction is load-bearing.** Every array
    #: above is a gDNA/RNA split that a perfect deconvolution would change; this one is a property of
    #: the partition and the fragment-length distribution and would be identical under any split. That
    #: is why it is NOT in ``prior_vs_oracle.OVERRIDE_FIELDS``: an oracle that overrode it would be
    #: answering a different question.
    #:
    #: ⭐ ``assemble_priors`` multiplies each component's per-line mass by it, because the accumulator
    #: deposits ``+1`` on EVERY line a fragment crosses — ``max(K, 1)`` of them — so a sum over lines is
    #: an object-incidence count and the EM adds a FRAGMENT count. It is 1.0 where both flanking nodes
    #: exceed every fragment length, and falls toward the node spacing where they do not.
    edge_mass_per_crossing: np.ndarray

    # --- the JUMPING population, per junction edge (float64[n_junctions]) ---
    #: ⭐ **Never deconvolved: a junction edge is pure mature RNA by construction**, so this is
    #: ``sj_count`` summed over the genome-strand columns and nothing else. It is the third population
    #: at a line, and it is routinely two orders of magnitude larger than ``mass_rna_spliced_edge`` at
    #: the same place: at a donor seam the junction flux is the gene's whole mature output while the
    #: spliced crossing is the handful of molecules that read through without splicing.
    #: ⚠ **``assemble_priors`` does NOT consume it, and that is deliberate.** Junction fragments are
    #: certified RNA in exactly the sense ``mass_rna_spliced_edge`` is withheld for, so feeding them to
    #: ``rna_prior_count`` would load the RNA side of a split that arbitrates only unspliced fragments.
    #: It is exported for QC and reporting — the calibration's output should not be silent about the
    #: population that dominates a donor seam (owner ruling, 2026-07-30).
    mass_rna_junction: np.ndarray

    #: float64[n_edges] — the ``edge_spliced`` twin of ``edge_mass_per_crossing``, and
    #: float64[n_junctions] — the junction one, ``sj_mass / sj_count``. Same kind of quantity as
    #: ``edge_mass_per_crossing`` in every respect: GEOMETRY, identical under any split, and therefore
    #: not in ``prior_vs_oracle.OVERRIDE_FIELDS``.
    #:
    #: ⛔ **The junction one did not exist until ``sj_mass`` did, and that is why a conserved LIBRARY
    #: fragment count was not computable.** A spliced fragment whose every block lies inside one node
    #: crosses no line and is not contained, so it deposited on no conserved bank — **1,222,375 of
    #: 4,830,713 RNA fragments (25.3 %)** on ladder g50 capture_off, against 0 of 4,997,761 gDNA
    #: fragments, since gDNA cannot splice.
    edge_spliced_mass_per_crossing: np.ndarray
    junction_mass_per_crossing: np.ndarray

    # --- the gDNA geometric supports: expected admissible START POSITIONS, per component ---
    #: float64[n_nodes] — ``effective_length.contained_eff_length`` on the gDNA pmf,
    #: ``E_f[(node_len − w + 1)+]``. Under uniform genomic gDNA at density ρ the expected contained
    #: mass is EXACTLY ``ρ · gdna_node_eff_len`` (a fragment must FIT to be contained), so dividing by
    #: it recovers ρ — the bedrock "factor 1 under uniform gDNA" invariant. This, NOT the genomic node
    #: length, is the density-correct divisor: the raw length ignores the fit-inside constraint and so
    #: understates a short node's density, manufacturing a spurious contraction in an unenriched library.
    gdna_node_eff_len: np.ndarray
    #: float64[n_edges] — ``effective_length.crossing_eff_length`` on the gDNA pmf. ⚠ **Uniform across
    #: edges today, and that is physics rather than a placeholder**: gDNA's template is the chromosome,
    #: so it takes ``UNBOUNDED_REACH`` on both sides at every line and the divisor collapses to
    #: ``mu_g − 1``. It stays a per-edge array because that is the axis its consumers index it on.
    gdna_edge_eff_len: np.ndarray

    # --- the RNA geometric supports: the SAME two frames, on the RNA pmf ---
    #: float64[n_nodes] / float64[n_edges] — ``contained_eff_length`` and ``crossing_eff_length`` on the
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
    #: ⭐ Like their gDNA twins they are PROJECTED off ``NodeGeometry.eff_rna`` by ``_project_eff``,
    #: never recomputed, so they are byte-identically the opportunity the SOLVER used. That is the
    #: property that makes them safe to reason with.
    rna_node_eff_len: np.ndarray
    rna_edge_eff_len: np.ndarray

    # --- library scalars ---
    gdna_density_global: float  # >= 0, global gDNA density (mass/bp); 0 in a zero-gDNA library
    rna_sense_frac: float  # in [0, 1], RNA sense fraction used by the strand clue
    gdna_strand_overdispersion: float  # in [0, 1), fitted gDNA strand Beta-Binomial dispersion
    rna_strand_overdispersion: float  # in [0, 1), fitted RNA strand Beta-Binomial dispersion

    # --- provenance: the three axis lengths, independent of each other ---
    n_nodes: int
    n_edges: int
    n_junctions: int
    config: CalibrationConfig

    def __post_init__(self) -> None:
        for axis in ("n_nodes", "n_edges", "n_junctions"):
            if int(getattr(self, axis)) < 0:
                raise ValueError(
                    f"CalibrationResult.{axis} must be >= 0; got {getattr(self, axis)}."
                )

        for name in ("mass_gdna_node", "mass_rna_node", "gdna_node_eff_len", "rna_node_eff_len"):
            _check_axis_array(getattr(self, name), name, self.n_nodes)
        for name in (
            "mass_gdna_edge",
            "mass_rna_edge",
            "mass_rna_spliced_edge",
            "edge_mass_per_crossing",
            "edge_spliced_mass_per_crossing",
            "gdna_edge_eff_len",
            "rna_edge_eff_len",
        ):
            _check_axis_array(getattr(self, name), name, self.n_edges)
        for name in ("mass_rna_junction", "junction_mass_per_crossing"):
            _check_axis_array(getattr(self, name), name, self.n_junctions)

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

    # ── ⭐⭐⭐ THE LIBRARY FRAGMENT COUNT — the deliverable, in FRAGMENT units, derived in ONE place ──

    @property
    def library_gdna_fragments(self) -> float:
        """gDNA fragments in the library — a CONSERVED count, not an object-incidence sum.

        gDNA appears on TWO axes only, because it cannot splice, and containment is exclusive — so a
        node term is already a fragment count and only the crossing term needs converting. That is why
        this side was exact all along while the RNA side was 25.3 % short.
        """
        return float(
            np.asarray(self.mass_gdna_node, dtype=np.float64).sum()
            + (
                np.asarray(self.mass_gdna_edge, dtype=np.float64)
                * np.asarray(self.edge_mass_per_crossing, dtype=np.float64)
            ).sum()
        )

    @property
    def library_rna_fragments(self) -> float:
        """RNA fragments in the library — all THREE axes, each converted by its OWN population's ``q``.

        ⛔⛔ **THE TREE ONCE HELD THREE DISAGREEING ANSWERS TO THIS**, and every one summed INCIDENCES:
        ``pipeline.py`` added node + edge + junction raw, ``calibration_truth_ab.py`` added node + edge
        and dropped the junction, and only ``assemble_priors`` converted anything. One fragment books
        ``max(K,1)`` line crossings AND one incidence per junction it uses, so the sum over-counts and
        the ratio is biased wherever the two components' inflations differ — which they do. Measured on
        the origin-split oracle, ladder g50 capture_off: the incidence ratio reads ``f_gdna`` **0.3851**
        against a truth of **0.5085**, while these conserved counts reproduce **0.5085** exactly, each
        origin landing at **1.000x** its deposited fragments with **0** unaccounted.

        ⭐ The spliced crossings and the junction flux are the SAME fragments split across two banks by
        the deposit rule — ``edge_spliced_mass`` holds the share of a spliced fragment's bases in blocks
        that crossed a line and ``sj_mass`` the share in blocks that crossed none — and the two sum to
        exactly one per fragment. Adding both is conservation, not double counting.

        ⚠ **A PROPERTY, never a stored field.** ``prior_vs_oracle`` swaps the mass arrays for truth with
        ``dataclasses.replace``; a cached scalar would survive that swap and silently describe the old
        arrays (``TRAPS: a-hash-that-misses-its-artifact``, in dataclass form). Deriving it means the
        oracle arm's count is the oracle's by construction.
        """
        unspliced_edge = np.maximum(
            np.asarray(self.mass_rna_edge, dtype=np.float64)
            - np.asarray(self.mass_rna_spliced_edge, dtype=np.float64),
            0.0,
        )
        return float(
            np.asarray(self.mass_rna_node, dtype=np.float64).sum()
            + (unspliced_edge * np.asarray(self.edge_mass_per_crossing, dtype=np.float64)).sum()
            + (
                np.asarray(self.mass_rna_spliced_edge, dtype=np.float64)
                * np.asarray(self.edge_spliced_mass_per_crossing, dtype=np.float64)
            ).sum()
            + (
                np.asarray(self.mass_rna_junction, dtype=np.float64)
                * np.asarray(self.junction_mass_per_crossing, dtype=np.float64)
            ).sum()
        )


__all__ = ["CalibrationResult"]
