"""
hulkrna.bayes_count - Bayesian read-to-gene assignment for RNA-seq.

Implements a two-pass Bayesian approach to assigning paired-end reads
to transcripts/genes:

Pass 1 (Learning):
    - Learn the strand protocol (FR/RF) and strand specificity p_sense
      from spliced reads via StrandModel
    - Collect initial raw counts per gene using uniform assignment
      (each read split equally among candidate genes)

Pass 2 (Bayesian Assignment):
    - For each fragment, compute posterior probability of originating
      from each candidate gene/transcript using:
        P(gene_j | fragment) ∝ P(fragment | gene_j) × P(gene_j)
    where:
        P(fragment | gene_j) = P(strand | gene_j) × P(overlap | gene_j)
        P(gene_j) = prior proportional to gene abundance from Pass 1

    The strand likelihood P(strand | gene_j) uses the learned p_sense:
        - If fragment strand is concordant with gene_j: p_sense
        - If fragment strand is discordant with gene_j: 1 - p_sense

    This allows proper handling of overlapping antisense genes where
    traditional tools would discard the read or assign it arbitrarily.

Theory
------
For a fragment f overlapping candidate genes G = {g1, g2, ..., gk}:

    P(gi | f) = P(f | gi) × P(gi) / Σ_j P(f | gj) × P(gj)

where:
    P(f | gi) = P(strand_f | strand_gi, p_sense) 
                × P(splice_junctions_f | gi)
                × P(exon_overlap_f | gi)

    P(gi) ∝ abundance(gi)  [from Pass 1 or uniform]

The abundance prior allows the model to correctly assign reads in
regions where a highly-expressed sense gene overlaps a lowly-expressed
antisense gene, even with imperfect strand specificity.
"""

import logging
from typing import List, Tuple, Optional, Dict, FrozenSet
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .core import Strand, IntervalType
from .core import CountCategory, CountStrand, CountType
from .strand_model import StrandModel
from .index import TranscriptIndex, merge_sets_with_relaxation
from .fragment import Fragment


logger = logging.getLogger(__name__)


# Small constant to prevent log(0) and division by zero
_EPS = 1e-300

# Minimum fractional assignment to count (prevents noise accumulation)
_MIN_FRAC = 1e-6


@dataclass
class BayesCounterStats:
    """Tracking statistics for the Bayesian counting process."""
    n_fragments: int = 0
    n_duplicate: int = 0
    n_chimeric: int = 0
    n_ambiguous_strand: int = 0
    n_exonic: int = 0
    n_intronic: int = 0
    n_intergenic: int = 0
    n_no_feature: int = 0
    n_spliced_no_feature: int = 0
    # Bayesian-specific
    n_unique_gene: int = 0       # fragment maps to exactly 1 gene
    n_multi_gene: int = 0        # fragment maps to >1 gene
    n_opposite_strand_genes: int = 0  # fragment overlaps genes on both strands

    def to_dict(self) -> dict:
        return {k: v for k, v in self.__dict__.items()}


class BayesCounter:
    """
    Bayesian read counter that assigns fractional counts to genes/transcripts
    based on posterior probabilities incorporating strand specificity and
    relative gene abundance.

    Parameters
    ----------
    txindex : TranscriptIndex
        The transcript/gene index with interval trees and splice junction map.
    strand_model : StrandModel
        The learned strand model from Pass 1.
    abundance_prior : np.ndarray or None
        Gene abundance prior (from Pass 1). If None, uses uniform prior.
    abundance_pseudocount : float
        Pseudocount added to all gene abundances to prevent zero priors.
    """

    def __init__(
        self,
        txindex: TranscriptIndex,
        strand_model: StrandModel,
        abundance_prior: Optional[np.ndarray] = None,
        abundance_pseudocount: float = 1.0,
    ):
        self.txindex = txindex
        self.strand_model = strand_model
        self.stats = BayesCounterStats()

        # Count arrays: float32 for fractional assignment
        self.t_counts = np.zeros(
            (txindex.num_transcripts, len(CountType)), dtype=np.float64
        )
        self.g_counts = np.zeros(
            (txindex.num_genes, len(CountType)), dtype=np.float64
        )

        # Gene abundance prior (will be updated after Pass 1 if provided)
        if abundance_prior is not None:
            self._g_abundance = abundance_prior.copy() + abundance_pseudocount
        else:
            # Uniform prior
            self._g_abundance = np.ones(txindex.num_genes, dtype=np.float64)

        # Normalize abundance to probabilities
        self._g_log_prior = np.log(self._g_abundance / self._g_abundance.sum() + _EPS)

    def set_abundance_prior(self, g_counts: np.ndarray, pseudocount: float = 1.0):
        """
        Update gene abundance prior (e.g., after Pass 1).

        Parameters
        ----------
        g_counts : np.ndarray
            Raw gene count vector (sum across all CountType columns).
        pseudocount : float
            Added to each gene to prevent zero priors.
        """
        self._g_abundance = g_counts + pseudocount
        self._g_log_prior = np.log(self._g_abundance / self._g_abundance.sum() + _EPS)
        logger.info(f'Updated gene abundance prior: '
                    f'min={self._g_abundance.min():.1f}, '
                    f'max={self._g_abundance.max():.1f}, '
                    f'median={np.median(self._g_abundance):.1f}')

    def _compute_strand_log_likelihood(
        self, fragment: Fragment, g_inds: List[int]
    ) -> np.ndarray:
        """
        Compute log P(strand | gene) for each candidate gene.

        For each candidate gene, determines whether the fragment's 
        alignment strand is concordant or discordant with the gene's 
        annotated strand, then returns the appropriate probability
        from the strand model.

        Parameters
        ----------
        fragment : Fragment
            The aligned fragment.
        g_inds : list of int
            Candidate gene indices.

        Returns
        -------
        np.ndarray
            Log-likelihood for each candidate gene (same length as g_inds).
        """
        log_liks = np.zeros(len(g_inds), dtype=np.float64)

        for i, g_idx in enumerate(g_inds):
            gene_strand = Strand(self.txindex.g_to_strand_arr[g_idx])

            # Use the fragment's combined strand vs gene strand
            p = self.strand_model.strand_likelihood(
                fragment.combined_strand, gene_strand
            )
            log_liks[i] = np.log(p + _EPS)

        return log_liks

    def _compute_splice_log_likelihood(
        self, fragment: Fragment, ref: str, 
        g_inds: List[int], t_inds: List[int]
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute splice junction match log-likelihood for candidate
        genes and transcripts.

        Transcripts/genes with annotated splice junctions matching the
        fragment get a boost; those without get a penalty.

        Returns
        -------
        t_log_liks, g_log_liks : np.ndarray, np.ndarray
        """
        if not fragment.is_spliced:
            # No splice junctions → uninformative, return zeros (log(1)=0)
            return (
                np.zeros(len(t_inds), dtype=np.float64),
                np.zeros(len(g_inds), dtype=np.float64),
            )

        # Search annotated splice junctions
        sj_t_inds, sj_g_inds = self.txindex.search_splice_junctions(
            ref, fragment.sj_tuples
        )

        # For transcripts: bonus if in sj_t_inds, small penalty otherwise
        sj_t_set = set(sj_t_inds)
        t_log_liks = np.array([
            np.log(0.95) if t in sj_t_set else np.log(0.05)
            for t in t_inds
        ], dtype=np.float64) if t_inds else np.zeros(0, dtype=np.float64)

        # For genes: bonus if in sj_g_inds
        sj_g_set = set(sj_g_inds)
        g_log_liks = np.array([
            np.log(0.95) if g in sj_g_set else np.log(0.05)
            for g in g_inds
        ], dtype=np.float64) if g_inds else np.zeros(0, dtype=np.float64)

        return t_log_liks, g_log_liks

    def _bayesian_assign(
        self,
        fragment: Fragment,
        ref: str,
        t_inds: List[int],
        g_inds: List[int],
        count_category: CountCategory,
    ):
        """
        Perform Bayesian assignment of the fragment to candidate genes/transcripts.

        Computes posterior probabilities and distributes fractional counts.
        """
        # --- Gene-level assignment ---
        if len(g_inds) == 1:
            # Unique gene → assign fully, just determine sense/antisense
            self.stats.n_unique_gene += 1
            g_idx = g_inds[0]
            gene_strand = Strand(self.txindex.g_to_strand_arr[g_idx])
            strand_type = self._determine_strand_type(
                fragment.combined_strand, gene_strand
            )
            ct = CountType(len(CountStrand) * count_category + strand_type).value
            self.g_counts[g_idx, ct] += 1.0

        elif len(g_inds) > 1:
            self.stats.n_multi_gene += 1

            # Check if genes are on opposite strands
            gene_strands = [
                Strand(self.txindex.g_to_strand_arr[g]) for g in g_inds
            ]
            has_opposite = len(set(gene_strands)) > 1
            if has_opposite:
                self.stats.n_opposite_strand_genes += 1

            # Compute log-posteriors
            g_log_prior = np.array(
                [self._g_log_prior[g] for g in g_inds], dtype=np.float64
            )
            g_strand_ll = self._compute_strand_log_likelihood(fragment, g_inds)
            _, g_sj_ll = self._compute_splice_log_likelihood(
                fragment, ref, g_inds, t_inds
            )

            g_log_posterior = g_log_prior + g_strand_ll + g_sj_ll
            # Normalize (log-sum-exp trick)
            g_log_posterior -= np.max(g_log_posterior)
            g_posterior = np.exp(g_log_posterior)
            g_posterior /= g_posterior.sum()

            # Assign fractional counts
            for i, g_idx in enumerate(g_inds):
                frac = g_posterior[i]
                if frac < _MIN_FRAC:
                    continue
                gene_strand = Strand(self.txindex.g_to_strand_arr[g_idx])
                strand_type = self._determine_strand_type(
                    fragment.combined_strand, gene_strand
                )
                ct = CountType(len(CountStrand) * count_category + strand_type).value
                self.g_counts[g_idx, ct] += frac

        # --- Transcript-level assignment ---
        if len(t_inds) == 1:
            t_idx = t_inds[0]
            tx_strand = Strand(self.txindex.t_to_strand_arr[t_idx])
            strand_type = self._determine_strand_type(
                fragment.combined_strand, tx_strand
            )
            ct = CountType(len(CountStrand) * count_category + strand_type).value
            self.t_counts[t_idx, ct] += 1.0

        elif len(t_inds) > 1:
            # Compute transcript posteriors
            # Map each transcript to its gene for the abundance prior
            t_g_map = self.txindex.t_to_g_arr[t_inds]
            t_log_prior = np.array(
                [self._g_log_prior[g] for g in t_g_map], dtype=np.float64
            )
            t_strand_ll = np.array([
                np.log(
                    self.strand_model.strand_likelihood(
                        fragment.combined_strand,
                        Strand(self.txindex.t_to_strand_arr[t])
                    ) + _EPS
                )
                for t in t_inds
            ], dtype=np.float64)

            t_sj_ll, _ = self._compute_splice_log_likelihood(
                fragment, ref, g_inds, t_inds
            )

            t_log_posterior = t_log_prior + t_strand_ll + t_sj_ll
            t_log_posterior -= np.max(t_log_posterior)
            t_posterior = np.exp(t_log_posterior)
            t_posterior /= t_posterior.sum()

            for i, t_idx in enumerate(t_inds):
                frac = t_posterior[i]
                if frac < _MIN_FRAC:
                    continue
                tx_strand = Strand(self.txindex.t_to_strand_arr[t_idx])
                strand_type = self._determine_strand_type(
                    fragment.combined_strand, tx_strand
                )
                ct = CountType(len(CountStrand) * count_category + strand_type).value
                self.t_counts[t_idx, ct] += frac

    def _determine_strand_type(
        self, fragment_strand: Strand, ref_strand: Strand
    ) -> CountStrand:
        """Determine sense/antisense/ambiguous for counting."""
        if fragment_strand == Strand.NONE or ref_strand == Strand.NONE:
            return CountStrand.AMBIGUOUS
        if fragment_strand == Strand.AMBIGUOUS:
            return CountStrand.AMBIGUOUS
        if fragment_strand == ref_strand:
            return CountStrand.SENSE
        return CountStrand.ANTISENSE

    def count_fragment(
        self, fragment: Fragment
    ):
        """
        Process a single fragment and assign counts.

        Parameters
        ----------
        fragment : Fragment
            The consolidated fragment from paired-end reads.
        """
        self.stats.n_fragments += 1

        # Skip chimeric fragments
        if fragment.is_chimeric:
            self.stats.n_chimeric += 1
            return

        # Skip ambiguous strand (both strands present)
        if fragment.is_ambiguous_strand:
            self.stats.n_ambiguous_strand += 1
            return

        # Get the single reference
        ref = next(iter(fragment.refs))

        # Search for overlapping exons
        t_inds, g_inds = self.txindex.search_intervals(
            ref, fragment.exon_intervals, IntervalType.EXON
        )

        if t_inds or g_inds:
            self.stats.n_exonic += 1
            t_inds = list(t_inds)
            g_inds = list(g_inds)

            # Determine count category based on splicing
            if fragment.is_spliced:
                # Check for annotated splice junctions
                sj_t_inds, sj_g_inds = self.txindex.search_splice_junctions(
                    ref, fragment.sj_tuples
                )
                if sj_t_inds or sj_g_inds:
                    # Merge exon and SJ matches
                    t_inds_merged, g_inds_merged = merge_sets_with_relaxation(
                        [frozenset(t_inds), sj_t_inds],
                        [frozenset(g_inds), sj_g_inds],
                    )
                    t_inds = list(t_inds_merged)
                    g_inds = list(g_inds_merged)
                    count_cat = CountCategory.SPLICED_ANNOT
                else:
                    count_cat = CountCategory.SPLICED_UNANNOT
            else:
                count_cat = CountCategory.UNSPLICED

            if t_inds or g_inds:
                self._bayesian_assign(
                    fragment, ref, t_inds, g_inds, count_cat
                )
            else:
                self.stats.n_spliced_no_feature += 1
        else:
            # No exon overlap → check introns
            t_inds, g_inds = self.txindex.search_intervals(
                ref, fragment.exon_intervals, IntervalType.INTRON
            )
            if t_inds or g_inds:
                self.stats.n_intronic += 1
                self._bayesian_assign(
                    fragment, ref, list(t_inds), list(g_inds),
                    CountCategory.INTRON,
                )
            else:
                # Check intergenic
                if self.txindex.search_intergenic(ref, fragment.exon_intervals):
                    self.stats.n_intergenic += 1
                else:
                    self.stats.n_no_feature += 1

    def get_gene_abundance_vector(self) -> np.ndarray:
        """
        Return total gene counts summed across all CountType columns.
        Used as the abundance prior for Pass 2.
        """
        return self.g_counts.sum(axis=1)

    def get_t_counts_df(self) -> pd.DataFrame:
        return pd.DataFrame(self.t_counts, columns=list(CountType.columns()))

    def get_g_counts_df(self) -> pd.DataFrame:
        return pd.DataFrame(self.g_counts, columns=list(CountType.columns()))

    def summary(self) -> dict:
        d = self.stats.to_dict()
        g_sum = pd.Series(
            self.g_counts.sum(axis=0).round(),
            index=list(CountType.columns()),
        ).to_dict()
        t_sum = pd.Series(
            self.t_counts.sum(axis=0).round(),
            index=list(CountType.columns()),
        ).to_dict()
        d['g_counts'] = g_sum
        d['t_counts'] = t_sum
        return d

    def summary_by_type(self) -> pd.DataFrame:
        g_sum = pd.DataFrame(
            [self.g_counts.sum(axis=0)],
            columns=list(CountType.columns()),
            index=['g_counts'],
        )
        t_sum = pd.DataFrame(
            [self.t_counts.sum(axis=0)],
            columns=list(CountType.columns()),
            index=['t_counts'],
        )
        return pd.concat([g_sum, t_sum], axis=0)
