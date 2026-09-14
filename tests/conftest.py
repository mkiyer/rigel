"""The suite-wide pytest fixtures, and the one custom option.

One fixture built on the minimal GTF in `_index_builder`: a session-scoped `TranscriptIndex`
(`mini_index`), which lets a structural test run with no real genomic data.
`--update-golden` is registered here because a pytest option has to be declared in a conftest, and
`tests/test_golden_output.py` reads it.
"""

import pytest

from _index_builder import MINI_GTF, build_test_index

# ---------------------------------------------------------------------------
# Global pytest options
# ---------------------------------------------------------------------------


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Regenerate golden output files instead of comparing.",
    )


# ---------------------------------------------------------------------------
# Fixtures over the minimal GTF (GENCODE-style, 1-based inclusive coordinates)
# ---------------------------------------------------------------------------
@pytest.fixture(scope="session")
def mini_index(tmp_path_factory):
    """Standard MINI_GTF TranscriptIndex: t0(3-exon), t1(2-exon), t2(neg strand).

    GTF layout (0-based half-open after parse):
      g1 (+): t0 exons (99,200),(299,400),(499,600)
              t1 exons (99,200),(499,600)
      g2 (-): t2 exons (999,1100),(1199,1300)

    Transcript indices: t0=0, t1=1, t2=2
    Gene indices:       g1=0, g2=1
    """
    return build_test_index(tmp_path_factory, MINI_GTF, name="mini_idx")
