"""rigel.native — Public interface to Rigel's C++ native extensions.

All C++ functionality used by Python code is imported through this module.
This is the single seam between Python orchestration and C++ hot paths.

Modules
-------
_bam_impl     : BAM scanning, annotation writing, SJ tag detection (htslib)
_resolve_impl : Fragment overlap resolution against the reference index
_scoring_impl : Per-fragment likelihood scoring (strand, coverage, splice)
_em_impl      : Locus-level EM solver, connected components, bias correction
_cgranges_impl: Interval overlap queries (vendored cgranges)
"""

# -- BAM scanning ----------------------------------------------------------
from ._bam_impl import BamScanner
from ._bam_impl import BamAnnotationWriter
from ._bam_impl import detect_sj_strand_tag

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
from ._em_impl import EM_PRIOR_EPSILON
from ._em_impl import build_partition_offsets
from ._em_impl import scatter_candidates_f64
from ._em_impl import scatter_candidates_i32
from ._em_impl import scatter_candidates_u8
from ._em_impl import scatter_units_f64
from ._em_impl import scatter_units_i32
from ._em_impl import scatter_units_u8

# -- Interval overlap -------------------------------------------------------
from ._cgranges_impl import cgranges

__all__ = [
    # BAM
    "BamScanner",
    "BamAnnotationWriter",
    "detect_sj_strand_tag",
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
    "EM_PRIOR_EPSILON",
    "build_partition_offsets",
    "scatter_candidates_f64",
    "scatter_candidates_i32",
    "scatter_candidates_u8",
    "scatter_units_f64",
    "scatter_units_i32",
    "scatter_units_u8",
    # Intervals
    "cgranges",
]
