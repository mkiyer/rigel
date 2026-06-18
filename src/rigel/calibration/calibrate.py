"""calibrate() — the fractional-accumulator calibrator (per-node simplex solve, iterated all-gDNA bootstrap).

The contained-region split is solved per node on the **2-simplex** ``(sense-RNA / antisense-RNA / gDNA)``
(calibration models RNA-vs-gDNA only; the per-locus EM separates nascent from mature downstream), iterated
over an **all-gDNA bootstrap**. Per the count-zero-information architecture (CALIBRATION_ARCHITECTURE.md §0),
a node's gDNA/RNA composition is determined by exactly two terms here — a **strand likelihood** (the
Beta-Binomial tilt of the per-strand counts, the only INTRINSIC per-node signal; the count enters ONLY as its
overdispersed Fisher information) and a **global gDNA prior** (the foundation, ``ρ_global`` at the MAD-spread
precision ``τ_global``, governing nodes the strand leaves undetermined)::

    substrate
      -> strand balance: rna_sense_frac (κ)
      -> node_gdna_density (RAW) -> fit gdna/rna strand Beta-Binomial overdispersions (seed)
      -> strand_deconvolve -> cleaned_gdna_count: clean the BOUNDARY crossings (nascent removal)
      -> node_gdna_density (RAW contained + CLEANED boundaries): the boundary anchors' count input
      -> deconv_sides ONCE: the fixed boundary gDNA anchors (the boundary-flux transport into priors)
      -> for each pass (signature-binary all-gDNA init):
           rho_global re-fit on the PREVIOUS pass's gDNA estimate (the iterating foundation)
           -> deconv_regions_sweep: the per-node strand + global simplex solve
           -> converge on per-node f_g
      -> gdna_density_global (the library-average density QC scalar)

Pass 0 is the **signature-binary all-gDNA init**: every node starts at ``f_gdna=1`` carrying its full unspliced
MASS (count-free; CALIBRATION_ARCHITECTURE §3) ⇒ ``ρ_global`` = the observable-node total mass density, a
deliberate over-estimate the iteration drives down as the strand removes RNA. Cross-node imputation (the
reliability-weighted neighbour propagation) is deferred to Step 2; until then a strand-weak node falls to the
global foundation (``ρ_global ≈ 0`` in a pure-RNA library ⇒ no phantom gDNA). The boundary SIDES are
deconvolved ONCE (``deconv_sides``) as the fixed gDNA anchors whose mass feeds the per-locus prior via the
boundary-flux transport (``priors``), post-loop. The library-average gDNA density (a QC scalar) is **derived**
from the aggregate deconvolved mass. The per-region gDNA length contraction under capture is the IPR of the
deconvolved gDNA mass over the per-region effective supports (``gdna_region_eff_len`` + ``gdna_boundary_len``),
applied later in ``priors.assemble_priors``.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

from .bp_solver import (
    build_node_geometry,
    build_node_statics,
    chain_boundary_side_deconv,
    chain_region_deconv,
    init_beliefs,
    node_sweep,
)
from .density_model import node_gdna_density
from .derive import gdna_density_global
from .errors import CalibrationStrandError
from .effective_length import (
    boundary_eff_length,
    boundary_side_eff_length,
    region_eff_length,
)
from .gdna_strand import (
    fit_gdna_strand_from_substrate,
    fit_rna_strand_from_substrate,
    overdispersion_for_beta,
)
from .node_chain import build_node_chain
from .result import CalibrationResult
from .strand_balance import fit_strand_balance
from .substrate import BoundarySubstrate, CalibrationSubstrate

if TYPE_CHECKING:
    from ..config import CalibrationConfig
    from ..scan_payload import AccumulatorPayload
    from ..strand_model import StrandModels
    from .region_arrays import RegionArrays

logger = logging.getLogger(__name__)


def calibrate(
    payload: "AccumulatorPayload",
    region_arrays: "RegionArrays",
    strand_model: "StrandModels",
    gdna_fl_pmf: "np.ndarray",
    rna_fl_pmf: "np.ndarray",
    config: "CalibrationConfig",
) -> CalibrationResult:
    """Deconvolve the library into gDNA / RNA per node, then derive gdna_density_global.

    Iterative all-gDNA bootstrap (``config.sweep_max_passes`` passes, re-fitting ``ρ_global`` + the
    gDNA var~mean each pass and converging on per-node ``f_g``); see the module docstring for the data
    flow. ``gdna_density_global``
    may be ``0`` (a zero-gDNA library) and a node's deconvolved gDNA mass may be ``0`` (a pure-RNA
    node); both are valid, graceful outputs — not failures.
    """
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)

    # gDNA fragment-length effective lengths (PR 4c geometry): the region-contained
    # length, the per-side boundary density length, and the region-free crossing mean.
    region_eff_len = region_eff_length(region_arrays.region_size_bp, gdna_fl_pmf)
    boundary_eff_len = boundary_side_eff_length(gdna_fl_pmf, region_arrays.region_size_bp)
    fl_mean = boundary_eff_length(gdna_fl_pmf)
    # ``rna_fl_pmf`` feeds the RNA-side effective lengths for the cross-node RNA imputation, which is
    # deferred to Step 2 (CALIBRATION_ARCHITECTURE §8); it is unused by the Step-1 strand+global solve.

    # RNA strand balance: rna_sense_frac (κ) = posterior-mean spliced sense fraction. The strand
    # channel's discriminability w=(2κ−1)² (set inside the deconv) is the smooth strand→count
    # deference weight — there is no hard identifiability gate (an unstranded library has κ≈½ ⇒ w≈0 ⇒
    # count governs, regardless of depth).
    balance = fit_strand_balance(strand_model)
    if balance.fallback_used:
        # No spliced reads at all — not a usable RNA-seq library. Fail loudly (a real RNA-seq library
        # always carries spliced reads); see CalibrationStrandError.
        raise CalibrationStrandError(
            "the library has zero spliced unique-mapper observations; this does not look like an "
            "RNA-seq library. A real RNA-seq library always carries spliced reads."
        )
    rna_sense_frac = float(balance.rna_sense_frac)

    # Count clue on RAW counts (the count module, pre-cleaning): per-region gDNA density by LOCAL
    # boundary-anchored imputation. Needed here only to fit the gDNA strand overdispersion (its seed
    # identification is pre-cleaning) — the cleaning below depends on that overdispersion, so the raw
    # pass must come first. This is NOT the region answer (the cleaned pass below is).
    node_density_raw = node_gdna_density(substrate, region_arrays, region_eff_len, fl_mean)

    # Strand-module parameters — the two Beta-Binomial overdispersions. gDNA (mean ½) fitted from the
    # count-observable seed regions/sides using the raw count-clue gDNA weight (breaks the circularity:
    # the seed weight is the strand MEAN ½, not the dispersion). RNA (mean κ) fitted from boundary-side
    # spliced counts. Both shrunk toward the SAME default prior, so under sparse data they collapse to
    # one distribution and an unstranded node (κ=½) is uninformative. See docs/em_strand/03+05.
    gdna_strand = fit_gdna_strand_from_substrate(
        substrate,
        region_arrays,
        node_density_raw,
        boundary_eff_len,
        rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(config.gdna_strand_prior_alpha_beta),
        prior_weight=config.gdna_strand_prior_weight,
    )
    gdna_strand_overdispersion = gdna_strand.gdna_strand_overdispersion
    rna_strand = fit_rna_strand_from_substrate(
        substrate,
        rna_sense_frac=rna_sense_frac,
        prior_overdispersion=overdispersion_for_beta(config.rna_strand_prior_alpha_beta),
        prior_weight=config.rna_strand_prior_weight,
    )
    rna_strand_overdispersion = rna_strand.rna_strand_overdispersion

    # THE SOLVE — the bipartite belief-propagation sweep (plan v3 / bp_solver): build the unified
    # region↔boundary chain + its per-node geometry / statics, the signature-binary init, then the directional
    # Gauss-Seidel sweep (gDNA + per-strand RNA identity-density messages, per-pass frozen-snapshot var~mean
    # reliability, the global gDNA prior). The region nodes give the per-region gDNA fraction; the boundary
    # nodes give the per-side boundary flux feeding the per-locus prior (chain_boundary_side_deconv) — the
    # first-class boundaries that retired the legacy deconv_sides / count-cleaning / I₀ path.
    boundary_substrate = BoundarySubstrate.from_payload(payload)
    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geometry = build_node_geometry(
        chain, substrate, boundary_substrate, region_arrays, gdna_fl_pmf, rna_fl_pmf
    )
    statics = build_node_statics(chain, substrate, boundary_substrate, region_arrays)
    belief = init_beliefs(
        chain, substrate, boundary_substrate, region_arrays,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_grid=config.sweep_n_grid, statics=statics,
    )
    belief, sweep_deltas = node_sweep(
        chain, statics, geometry, belief, region_arrays,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_grid=config.sweep_n_grid, max_passes=config.sweep_max_passes,
        convergence_delta=config.sweep_convergence_delta,
    )
    regions = chain_region_deconv(chain, belief, substrate)
    left, right = chain_boundary_side_deconv(chain, belief, substrate)
    logger.debug("calibration sweep: %d passes, max-Δf_g per pass=%s", len(sweep_deltas), sweep_deltas)

    # Derive gdna_density_global (the library-average density QC scalar).
    density_global = gdna_density_global(
        regions, left, right, region_eff_len, boundary_eff_len, region_arrays.ref_id
    )

    # Spliced RNA mass per region (Σ over the 3 nodes). The deconv adds the full per-node
    # spliced mass to rna_mass (rna = (1−g)·M_unspliced + M_spliced), so this is exactly the
    # spliced component of mass_rna_*. assemble_priors withholds it from rna_prior_count — a spliced
    # fragment is guaranteed-RNA in the EM (no gDNA candidate), so it must not load the RNA side of
    # the gDNA-vs-RNA *unspliced* split (the prior's a_g+a_r should equal the unspliced competing
    # mass). mass_rna_* stays spliced-inclusive so node conservation gdna+rna = total holds.
    mass_rna_spliced = (
        np.asarray(substrate.contained.mass_spliced, dtype=np.float64)
        + np.asarray(substrate.left.mass_spliced, dtype=np.float64)
        + np.asarray(substrate.right.mass_spliced, dtype=np.float64)
    )

    result = CalibrationResult(
        mass_gdna_contained=regions.gdna_mass,
        mass_rna_contained=regions.rna_mass,
        mass_gdna_left=left.gdna_mass,
        mass_rna_left=left.rna_mass,
        mass_gdna_right=right.gdna_mass,
        mass_rna_right=right.rna_mass,
        mass_rna_spliced=mass_rna_spliced,
        gdna_boundary_len=boundary_eff_len,
        gdna_region_eff_len=np.asarray(region_eff_len, dtype=np.float64),
        gdna_density_global=density_global,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        n_regions=region_arrays.n_regions,
        config=config,
    )
    # Diagnostic: the boundary-spliced sense fraction should agree with the StrandModel κ (the
    # deconv mean). A large gap flags a strand-model / accumulator mismatch (we do NOT refit the
    # mean from boundary spliced — κ stays the StrandModel posterior; this is QC only).
    spl_sense = float(
        np.sum(substrate.left.n_spliced_sense) + np.sum(substrate.right.n_spliced_sense)
    )
    spl_total = spl_sense + float(
        np.sum(substrate.left.n_spliced_antisense) + np.sum(substrate.right.n_spliced_antisense)
    )
    boundary_sense_frac = spl_sense / spl_total if spl_total > 0.0 else float("nan")
    logger.debug(
        "calibration: R=%d gdna_density_global=%.4g rna_sense_frac=%.3f "
        "gdna_strand_overdispersion=%.4g (%d seed nodes, %d frags%s) "
        "rna_strand_overdispersion=%.4g (%d seed sides, %d frags%s) "
        "[boundary-spliced sense_frac=%.3f vs κ=%.3f]",
        result.n_regions,
        result.gdna_density_global,
        rna_sense_frac,
        gdna_strand_overdispersion,
        gdna_strand.n_seed_nodes,
        gdna_strand.n_seed_fragments,
        ", FALLBACK" if gdna_strand.fallback_used else "",
        rna_strand_overdispersion,
        rna_strand.n_seed_nodes,
        rna_strand.n_seed_fragments,
        ", FALLBACK" if rna_strand.fallback_used else "",
        boundary_sense_frac,
        rna_sense_frac,
    )
    return result


__all__ = ["calibrate"]
