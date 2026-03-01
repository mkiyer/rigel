"""
Shared test fixtures and helpers for hulkrna.

Provides small synthetic GTF data, minimal FASTA files, and helper
functions so that tests can run without real genomic data.
"""

import textwrap
from pathlib import Path

import numpy as np
import pytest

from hulkrna.categories import SpliceStrandCol
from hulkrna.estimator import (
    AbundanceEstimator,
    EM_PRIOR_EPSILON,
    Locus,
    LocusEMInput,
)

# ---------------------------------------------------------------------------
# Minimal GTF content (GENCODE-style, 1-based inclusive coordinates)
# ---------------------------------------------------------------------------
# Two genes on chr1, one on + and one on - strand.
# Gene A (g1): two transcripts (t1 has 3 exons, t2 has 2 exons)
# Gene B (g2): one transcript (t3 has 2 exons), negative strand
#
# Coordinates (1-based inclusive as in GTF):
#   t1: exons at 100-200, 300-400, 500-600
#   t2: exons at 100-200, 500-600
#   t3: exons at 1000-1100, 1200-1300 (neg strand)

MINI_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t500\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "GeneA"; gene_type "protein_coding"; tag "basic";
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "GeneA"; gene_type "protein_coding";
    chr1\ttest\texon\t500\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t2"; gene_name "GeneA"; gene_type "protein_coding";
    chr1\ttest\texon\t1000\t1100\t.\t-\t.\tgene_id "g2"; transcript_id "t3"; gene_name "GeneB"; gene_type "lncRNA";
    chr1\ttest\texon\t1200\t1300\t.\t-\t.\tgene_id "g2"; transcript_id "t3"; gene_name "GeneB"; gene_type "lncRNA";
""")

# Expected 0-based half-open coordinates after GTF parsing:
#   t1: exons (99,200), (299,400), (499,600)
#   t2: exons (99,200), (499,600)
#   t3: exons (999,1100), (1199,1300)


@pytest.fixture
def mini_gtf_file(tmp_path: Path) -> Path:
    """Write the minimal GTF to a temp file and return its path."""
    gtf_path = tmp_path / "test.gtf"
    gtf_path.write_text(MINI_GTF)
    return gtf_path


@pytest.fixture
def mini_fasta_file(tmp_path: Path) -> Path:
    """Create a minimal FASTA + .fai so pysam can read reference lengths.

    chr1: 2000 bp, chr2: 500 bp (chr2 has no genes — tests intergenic).
    """
    fasta_path = tmp_path / "genome.fa"
    # Write a simple FASTA with 80-char line width
    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        seq1 = "N" * 2000
        for i in range(0, len(seq1), 80):
            f.write(seq1[i : i + 80] + "\n")
        f.write(">chr2\n")
        seq2 = "N" * 500
        for i in range(0, len(seq2), 80):
            f.write(seq2[i : i + 80] + "\n")

    # Create the .fai index
    import pysam
    pysam.faidx(str(fasta_path))

    return fasta_path


# ---------------------------------------------------------------------------
# Shared EM helpers (used by test_estimator and test_gdna)
# ---------------------------------------------------------------------------

# Default column for UNSPLICED_SENSE
_UNSPLICED_SENSE = int(SpliceStrandCol.UNSPLICED_SENSE)


def _make_locus_em_data(
    t_indices_per_unit,
    log_liks_per_unit=None,
    count_cols_per_unit=None,
    num_transcripts=None,
    rc=None,
    nrna_init=None,
    gdna_init=0.0,
    include_nrna=False,
    include_gdna=False,
    nrna_log_lik=-2.0,
    gdna_log_lik=0.0,
):
    """Build a LocusEMInput for unit tests (single-locus, identity mapping).

    Parameters
    ----------
    t_indices_per_unit : list[list[int]]
        LOCAL mRNA transcript indices per unit.
    rc : AbundanceEstimator or None
        If provided, unique_totals for mRNA are extracted from
        rc.unique_counts[i].sum().
    include_nrna : bool
        Add nRNA candidates for each mRNA candidate.
    include_gdna : bool
        Add a gDNA candidate for each unit.
    """
    if num_transcripts is None:
        all_t = [t for unit in t_indices_per_unit for t in unit]
        num_transcripts = (max(all_t) + 1) if all_t else 1

    n_t = num_transcripts
    n_components = 2 * n_t + 1
    gdna_idx = 2 * n_t

    offsets = [0]
    flat_t = []
    flat_lk = []
    flat_cc = []
    locus_t_list = []
    locus_cc_list = []

    for u, t_list in enumerate(t_indices_per_unit):
        for j, t_idx in enumerate(t_list):
            flat_t.append(t_idx)
            flat_lk.append(
                log_liks_per_unit[u][j] if log_liks_per_unit else 0.0
            )
            flat_cc.append(
                count_cols_per_unit[u][j]
                if count_cols_per_unit
                else _UNSPLICED_SENSE
            )

        if include_nrna:
            for t_idx in t_list:
                flat_t.append(n_t + t_idx)
                flat_lk.append(nrna_log_lik)
                flat_cc.append(_UNSPLICED_SENSE)

        if include_gdna:
            flat_t.append(gdna_idx)
            flat_lk.append(gdna_log_lik)
            flat_cc.append(_UNSPLICED_SENSE)

        offsets.append(len(flat_t))

        locus_t_list.append(t_list[0] if t_list else -1)
        cc = (
            count_cols_per_unit[u][0]
            if (count_cols_per_unit and t_list)
            else _UNSPLICED_SENSE
        )
        locus_cc_list.append(cc)

    n_units = len(t_indices_per_unit)

    ut = np.zeros(n_components, dtype=np.float64)
    if rc is not None:
        for i in range(min(n_t, rc.num_transcripts)):
            ut[i] = rc.unique_counts[i].sum()

    if nrna_init is None:
        nrna_arr = np.zeros(n_t, dtype=np.float64)
    else:
        nrna_arr = np.asarray(nrna_init, dtype=np.float64)
    for i in range(n_t):
        ut[n_t + i] = nrna_arr[i]
    ut[gdna_idx] = gdna_init

    eff_len = np.ones(n_components, dtype=np.float64)
    prior = np.full(n_components, EM_PRIOR_EPSILON, dtype=np.float64)

    locus = Locus(
        locus_id=0,
        transcript_indices=np.arange(n_t, dtype=np.int32),
        gene_indices=np.array([0], dtype=np.int32),
        unit_indices=np.arange(n_units, dtype=np.int32),
    )

    n_local_candidates = len(flat_t)
    # Build component-length array of length 1 so that with tx_starts=0
    # and tx_ends=1, frag_len=1 and eff_len = max(1-1+1,1) = 1 → no
    # correction.  Phase H: pass int64 array instead of BiasProfile list.
    bias_profiles = np.ones(n_components, dtype=np.int64)
    return LocusEMInput(
        locus=locus,
        offsets=np.array(offsets, dtype=np.int64),
        t_indices=np.array(flat_t, dtype=np.int32),
        log_liks=np.array(flat_lk, dtype=np.float64),
        count_cols=np.array(flat_cc, dtype=np.uint8),
        coverage_weights=np.ones(n_local_candidates, dtype=np.float64),
        tx_starts=np.zeros(n_local_candidates, dtype=np.int32),
        tx_ends=np.ones(n_local_candidates, dtype=np.int32),
        locus_t_indices=np.array(locus_t_list, dtype=np.int32),
        locus_count_cols=np.array(locus_cc_list, dtype=np.uint8),
        n_transcripts=n_t,
        n_components=n_components,
        local_to_global_t=np.arange(n_t, dtype=np.int32),
        unique_totals=ut,
        nrna_init=nrna_arr,
        gdna_init=gdna_init,
        effective_lengths=eff_len,
        prior=prior,
        bias_profiles=bias_profiles,
    )


def _run_and_assign(rc, locus_em, *, em_iterations=10):
    """Convenience: run locus EM then assign. Returns (theta, gdna_count)."""
    theta, _alpha = rc.run_locus_em(locus_em, em_iterations=em_iterations)
    gdna_count = rc.assign_locus_ambiguous(
        locus_em, theta,
    )
    return theta, gdna_count
