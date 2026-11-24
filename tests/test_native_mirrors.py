"""The C++ mirrors of Python's codes (``native/constants.h``, ``native/transfer_kernel.h``) against their Python
definitions. A mirror that drifted would mislabel silently: one strand read as another, a ZF tag stamped with another
outcome, a terminus read as a junction. The splice types' parity is in ``test_splice.py``."""

from rigel import annotate
from rigel._resolve_impl import testing as native_resolve
from rigel._solve_impl import testing as native_solve
from rigel.calibration import splice_graph
from rigel.types import IntervalType, Strand


def test_the_native_strand_codes_are_rigel_types_strand():
    for strand in Strand:
        assert getattr(native_resolve, f"STRAND_{strand.name}") == int(strand), strand.name


def test_the_native_interval_types_are_rigel_types_interval_type():
    assert native_resolve.ITYPE_EXON == int(IntervalType.EXON)
    assert native_resolve.ITYPE_TRANSCRIPT == int(IntervalType.TRANSCRIPT)


def test_the_zf_codes_the_native_stamp_sites_write_are_rigel_annotates():
    python = {
        "AF_RESOLVED_BIT": annotate.AF_RESOLVED_BIT,
        "AF_GDNA_BIT": annotate.AF_GDNA_BIT,
        "AF_INTERGENIC_BIT": annotate.AF_INTERGENIC_BIT,
        "AF_MULTIMAPPER_DROP_BIT": annotate.AF_MULTIMAPPER_DROP_BIT,
        "AF_UNRESOLVED": annotate.AF_UNRESOLVED,
        "AF_GDNA_INTERGENIC": annotate.AF_GDNA_INTERGENIC,
        "AF_MULTIMAPPER_DROP": annotate.AF_MULTIMAPPER_DROP,
    }
    for name, value in python.items():
        assert getattr(native_resolve, name) == value, name


def test_the_boundary_flags_the_transfer_kernel_reads_are_the_splice_graphs():
    python = {
        "TSS_POS": splice_graph.FLAG_TSS_POS,
        "TSS_NEG": splice_graph.FLAG_TSS_NEG,
        "TES_POS": splice_graph.FLAG_TES_POS,
        "TES_NEG": splice_graph.FLAG_TES_NEG,
        "DON_POS": splice_graph.FLAG_DONOR_POS,
        "DON_NEG": splice_graph.FLAG_DONOR_NEG,
        "ACC_POS": splice_graph.FLAG_ACCEPTOR_POS,
        "ACC_NEG": splice_graph.FLAG_ACCEPTOR_NEG,
    }
    assert dict(native_solve.BOUNDARY_FLAGS) == {name: int(value) for name, value in python.items()}
