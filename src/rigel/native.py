"""rigel.native — the one boundary between the Python package and its C++ extensions.

Every C++ symbol the package uses is imported through this module:

_bam_impl     : BAM scanning, the calibration accumulator, annotation writing (htslib)
_resolve_impl : fragment overlap resolution against the reference index
_scoring_impl : per-fragment likelihood scoring (strand, coverage, splice)
_em_impl      : the locus-level EM, the partition scatter that builds its loci's arrays, connected components
_solve_impl   : the calibration solve: the sweep's block pipeline in one call (``solve_blocks``), ψ (``psi_solve``),
                the cut shares of the effective lengths (``side_shares``), and ``rows``: ``read_column``,
                ``ARM_EPS`` and ``T_NATS``, the kernel's values its Python readers share

Every extension carries a ``testing`` submodule, its seams bound for the tests; production never imports it.
"""

# -- BAM scanning ----------------------------------------------------------
from ._bam_impl import BamScanner
from ._bam_impl import BamAnnotationWriter

# a record the scan cannot read as paired-end input (a ValueError)
from ._bam_impl import UnpairedRecordError

# -- The accumulator: the native class itself, so a caller reads the banks it exports -----------------
from ._bam_impl import Accumulator

# the held bank as the accumulator's batch reads take it
from ._bam_impl import HeldBank

# the kinds of object a placement's mass lands on, which the place weight reads
from ._bam_impl import PLACE_REGION, PLACE_BOUNDARY, PLACE_SJ

# -- Fragment resolution ----------------------------------------------------
from ._resolve_impl import FragmentResolver

# -- Scoring ----------------------------------------------------------------
from ._scoring_impl import NativeFragmentScorer
from ._scoring_impl import StreamingScorer

# the fragment classes the scorer routes on, which rigel.buffer reads
from ._scoring_impl import FRAG_UNAMBIG
from ._scoring_impl import FRAG_AMBIG_SAME_STRAND
from ._scoring_impl import FRAG_AMBIG_OPP_STRAND
from ._scoring_impl import FRAG_MULTIMAPPER

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

# the EM's assignment codes, which rigel.estimator reads; the EM refuses any other
from ._em_impl import ASSIGN_FRACTIONAL
from ._em_impl import ASSIGN_SAMPLE

# the winner code a unit the intergenic term wins carries, beside -2 for the locus's own gDNA
from ._em_impl import EM_WINNER_GDNA, EM_WINNER_INTERGENIC

# the columns of the EM's per-unit class masses, by name, in order
from ._em_impl import EM_CLASSES

# -- The calibration solve: the sweep's block pipeline and ψ --------------------------------------
from ._solve_impl import solve_blocks
from ._solve_impl import psi_solve

# one side's conserved cut share per cut (`effective_length._side_share`)
from ._solve_impl import side_shares

# the policy codes the block pipeline switches on; it refuses any other
from ._solve_impl import POLICY_SILENT
from ._solve_impl import POLICY_TRANSFER

# `read_column` and `ARM_EPS`, which the transfer policy and the landscape read
from ._solve_impl import rows as transfer_rows

__all__ = [
    # BAM
    "BamScanner",
    "BamAnnotationWriter",
    "UnpairedRecordError",
    # The accumulator
    "Accumulator",
    "HeldBank",
    "PLACE_REGION",
    "PLACE_BOUNDARY",
    "PLACE_SJ",
    # Resolution
    "FragmentResolver",
    # Scoring
    "NativeFragmentScorer",
    "StreamingScorer",
    "FRAG_UNAMBIG",
    "FRAG_AMBIG_SAME_STRAND",
    "FRAG_AMBIG_OPP_STRAND",
    "FRAG_MULTIMAPPER",
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
    "EM_WINNER_GDNA",
    "EM_WINNER_INTERGENIC",
    "EM_CLASSES",
    # The solve
    "solve_blocks",
    "POLICY_SILENT",
    "POLICY_TRANSFER",
    "transfer_rows",
    # ψ
    "psi_solve",
    "side_shares",
]
