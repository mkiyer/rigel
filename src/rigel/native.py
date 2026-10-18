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
               call allocates (for the gates), and `rows`: the row constructors, flag predicates and absorbed
               constructions bound for the gates, with `read_column` and the shared constants production reads
_cgranges_impl: Interval overlap queries (vendored cgranges)
"""

# -- BAM scanning ----------------------------------------------------------
from ._bam_impl import BamScanner
from ._bam_impl import BamAnnotationWriter

# -- The accumulator --------------------------------------------------------
# The native class directly: there is no Python row-view façade in front of it, so a caller reads the
# banks the accumulator actually exports.
from ._bam_impl import Accumulator, AccumulatorSet

# -- Fragment resolution ----------------------------------------------------
from ._resolve_impl import FragmentResolver
from ._resolve_impl import FragmentAccumulator
from ._resolve_impl import ResolvedFragment

# -- Scoring ----------------------------------------------------------------
from ._scoring_impl import NativeFragmentScorer
from ._scoring_impl import StreamingScorer

# the fragment classes the scorer routes on, which rigel.buffer reads
from ._scoring_impl import FRAG_UNAMBIG
from ._scoring_impl import FRAG_AMBIG_SAME_STRAND
from ._scoring_impl import FRAG_AMBIG_OPP_STRAND
from ._scoring_impl import FRAG_MULTIMAPPER
from ._scoring_impl import FRAG_CHIMERIC

# -- EM solver --------------------------------------------------------------
from ._em_impl import batch_locus_em_partitioned
from ._em_impl import connected_components
from ._em_impl import build_partition_offsets
from ._em_impl import scatter_candidates_f32
from ._em_impl import scatter_candidates_i32
from ._em_impl import scatter_candidates_u8
from ._em_impl import scatter_units_f32
from ._em_impl import scatter_units_i64
from ._em_impl import scatter_units_u8

# the EM's assignment and warm-start codes, which rigel.estimator reads; the EM refuses any other
from ._em_impl import ASSIGN_FRACTIONAL
from ._em_impl import ASSIGN_SAMPLE
from ._em_impl import WARM_START_COVERAGE
from ._em_impl import WARM_START_PRIOR
from ._em_impl import WARM_START_UNIFORM

# -- The calibration solve: the sweep's block pipeline, ψ, the transfer kernels for the gates -----
from ._solve_impl import solve_blocks
from ._solve_impl import psi_solve
from ._solve_impl import psi_cube as psi_cube_native
from ._solve_impl import posterior_median as psi_posterior_median
from ._solve_impl import compose as psi_compose
from ._solve_impl import gdna_arm as psi_gdna_arm

# the transfer's three kernels on tables the call allocates, bound for the gates
from ._solve_impl import transfer_prepare
from ._solve_impl import transfer_pass
from ._solve_impl import transfer_solve

# the row constructors, the flag predicates and the constructions the block pipeline absorbed,
# for the gates; and `read_column` and the shared constants, which production reads
from ._solve_impl import rows as transfer_rows

# -- Interval overlap -------------------------------------------------------
from ._cgranges_impl import cgranges

__all__ = [
    # BAM
    "BamScanner",
    "BamAnnotationWriter",
    # The accumulator
    "Accumulator",
    "AccumulatorSet",
    # Resolution
    "FragmentResolver",
    "FragmentAccumulator",
    "ResolvedFragment",
    # Scoring
    "NativeFragmentScorer",
    "StreamingScorer",
    "FRAG_UNAMBIG",
    "FRAG_AMBIG_SAME_STRAND",
    "FRAG_AMBIG_OPP_STRAND",
    "FRAG_MULTIMAPPER",
    "FRAG_CHIMERIC",
    # EM
    "batch_locus_em_partitioned",
    "connected_components",
    "build_partition_offsets",
    "scatter_candidates_f32",
    "scatter_candidates_i32",
    "scatter_candidates_u8",
    "scatter_units_f32",
    "scatter_units_i64",
    "scatter_units_u8",
    "ASSIGN_FRACTIONAL",
    "ASSIGN_SAMPLE",
    "WARM_START_COVERAGE",
    "WARM_START_PRIOR",
    "WARM_START_UNIFORM",
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
