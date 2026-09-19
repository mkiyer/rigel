"""rigel.native — Public interface to Rigel's C++ native extensions.

All C++ functionality used by Python code is imported through this module.
This is the single boundary between Python orchestration and C++ hot paths.

Modules
-------
_bam_impl     : BAM scanning, annotation writing, SJ tag detection (htslib)
_resolve_impl : Fragment overlap resolution against the reference index
_scoring_impl : Per-fragment likelihood scoring (strand, coverage, splice)
_em_impl      : Locus-level EM solver, connected components, effective-length normalization
_solve_impl   : The calibration solve — the sweep's block pipeline in one call (`solve_blocks`), ψ (`psi_solve`
               and its pieces for the gates), the composition transfer's builders, pass and solve on tables the
               call allocates (for the gates), and `rows`, the row constructors, flag predicates and absorbed
               constructions bound for the gates
_cgranges_impl: Interval overlap queries (vendored cgranges)
"""

# -- BAM scanning ----------------------------------------------------------
from ._bam_impl import BamScanner
from ._bam_impl import BamAnnotationWriter
from ._bam_impl import detect_sj_strand_tag

# -- The accumulator --------------------------------------------------------
# The native class directly: there is no Python row-view façade in front of it, so a caller reads the
# banks the accumulator actually exports.
from ._bam_impl import Accumulator

# -- Fragment resolution ----------------------------------------------------
from ._resolve_impl import FragmentResolver
from ._resolve_impl import FragmentAccumulator
from ._resolve_impl import ResolvedFragment

# -- Scoring ----------------------------------------------------------------
from ._scoring_impl import NativeFragmentScorer
from ._scoring_impl import StreamingScorer

# -- EM solver --------------------------------------------------------------
from ._em_impl import batch_locus_em_partitioned
from ._em_impl import connected_components
from ._em_impl import build_partition_offsets
from ._em_impl import scatter_candidates_f32
from ._em_impl import scatter_candidates_f64
from ._em_impl import scatter_candidates_i32
from ._em_impl import scatter_candidates_u8
from ._em_impl import scatter_units_f32
from ._em_impl import scatter_units_f64
from ._em_impl import scatter_units_i32
from ._em_impl import scatter_units_i64
from ._em_impl import scatter_units_u8

# -- The calibration solve: the sweep's block pipeline, ψ, the transfer kernels for the gates -----
from ._solve_impl import solve_blocks
from ._solve_impl import psi_solve
from ._solve_impl import psi_cube as psi_cube_native
from ._solve_impl import posterior_median as psi_posterior_median
from ._solve_impl import compose as psi_compose
from ._solve_impl import gdna_arm as psi_gdna_arm

# the transfer's three kernels on tables the call allocates, the row constructors, the flag predicates and the
# constructions the block pipeline absorbed — bound for the gates (nothing in src/ reads them)
from ._solve_impl import transfer_prepare
from ._solve_impl import transfer_pass
from ._solve_impl import transfer_solve
from ._solve_impl import rows as transfer_rows

# -- Interval overlap -------------------------------------------------------
from ._cgranges_impl import cgranges

__all__ = [
    # BAM
    "BamScanner",
    "BamAnnotationWriter",
    "detect_sj_strand_tag",
    # The accumulator
    "Accumulator",
    # Resolution
    "FragmentResolver",
    "FragmentAccumulator",
    "ResolvedFragment",
    # Scoring
    "NativeFragmentScorer",
    "StreamingScorer",
    # EM
    "batch_locus_em_partitioned",
    "connected_components",
    "build_partition_offsets",
    "scatter_candidates_f32",
    "scatter_candidates_f64",
    "scatter_candidates_i32",
    "scatter_candidates_u8",
    "scatter_units_f32",
    "scatter_units_f64",
    "scatter_units_i32",
    "scatter_units_i64",
    "scatter_units_u8",
    # The solve
    "solve_blocks",
    "transfer_prepare",
    "transfer_pass",
    "transfer_solve",
    "transfer_rows",
    # ψ
    "psi_solve",
    "psi_cube_native",
    "psi_posterior_median",
    "psi_compose",
    "psi_gdna_arm",
    # Intervals
    "cgranges",
]
