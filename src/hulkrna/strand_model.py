"""
hulkrna.strand_model - Bayesian strand protocol inference

Learns the probability that read1/read2 is sense vs antisense from spliced
reads whose splice-junction XS tag provides ground-truth gene strand.

The model produces a posterior Beta distribution over the strand-specificity
parameter `p_sense` (probability that a read aligning in the "expected"
orientation is truly sense).

Usage
-----
First pass (learning):
    model = StrandModel()
    for r1, r2 in parse_bam_file(bamfh, stats):
        model.observe(r1, r2)
    model.finalize()

Second pass (querying):
    # P(fragment originated from gene on strand s)
    p = model.p_sense_given_alignment(read_strand, gene_strand)

Mathematical background
-----------------------
We define:
    p_sense = P(read is on the sense strand of the gene it came from)

For a well-prepared strand-specific library (e.g. dUTP/RF), p_sense ≈ 0.95.
For unstranded data, p_sense ≈ 0.50.

We place a Beta(α, β) prior on p_sense and update it with observations from
spliced reads where the XS tag tells us the true gene strand:

    observation "concordant": read orientation matches the expected sense
    orientation for the detected protocol → increment α
    
    observation "discordant": read orientation is opposite → increment β

After the first pass:
    E[p_sense] = α / (α + β)

During counting (second pass), for a fragment overlapping genes on opposite
strands, the strand likelihood is:
    P(read | gene_strand=sense)  = p_sense
    P(read | gene_strand=antisense) = 1 - p_sense

Combined with gene abundance priors, this feeds into the Bayesian assignment.
"""

import logging
import collections
from enum import Enum
from typing import Optional, Tuple
from dataclasses import dataclass, field

import numpy as np

from .core import Strand


logger = logging.getLogger(__name__)


class ProtocolType(Enum):
    """Detected strand protocol orientation."""
    UNSTRANDED = 'unstranded'
    FR = 'fr'   # read1 sense, read2 antisense (e.g. KAPA Stranded)
    RF = 'rf'   # read1 antisense, read2 sense (e.g. dUTP, Illumina TruSeq Stranded)


@dataclass
class StrandModelResult:
    """Result of strand protocol inference with full posterior."""
    protocol: ProtocolType
    # Beta posterior parameters
    alpha: float
    beta: float
    # Derived quantities
    n_concordant: int = 0
    n_discordant: int = 0

    @property
    def n_total(self) -> int:
        return self.n_concordant + self.n_discordant

    @property
    def p_sense(self) -> float:
        """Expected probability that a read in the 'expected' orientation
        is truly sense (posterior mean of Beta distribution)."""
        return self.alpha / (self.alpha + self.beta)

    @property
    def p_antisense(self) -> float:
        """Complement: probability of strand leakage."""
        return 1.0 - self.p_sense

    @property
    def strand_specificity(self) -> float:
        """How strand-specific is the library? 
        1.0 = perfectly strand-specific, 0.5 = unstranded."""
        return max(self.p_sense, self.p_antisense)

    def posterior_variance(self) -> float:
        """Variance of the Beta posterior."""
        a, b = self.alpha, self.beta
        return (a * b) / ((a + b) ** 2 * (a + b + 1))

    def posterior_95ci(self) -> Tuple[float, float]:
        """95% credible interval for p_sense using Beta quantiles."""
        from scipy.stats import beta as beta_dist
        return (
            beta_dist.ppf(0.025, self.alpha, self.beta),
            beta_dist.ppf(0.975, self.alpha, self.beta),
        )

    def to_dict(self) -> dict:
        d = {
            'protocol': self.protocol.value,
            'p_sense': round(self.p_sense, 6),
            'p_antisense': round(self.p_antisense, 6),
            'strand_specificity': round(self.strand_specificity, 6),
            'alpha': round(self.alpha, 4),
            'beta': round(self.beta, 4),
            'n_concordant': self.n_concordant,
            'n_discordant': self.n_discordant,
            'n_total': self.n_total,
            'posterior_variance': round(self.posterior_variance(), 8),
        }
        return d

    def p_strand_given_alignment(self, read_genomic_strand: Strand, 
                                  gene_strand: Strand) -> float:
        """
        Probability that a read truly originated from `gene_strand`,
        given its observed genomic alignment strand and the learned 
        library protocol.

        For RF protocol: 
          - read1 on NEG strand → sense to POS-strand gene → p_sense
          - read1 on NEG strand → antisense to NEG-strand gene → 1-p_sense
        
        This handles the mapping through the protocol orientation.

        Parameters
        ----------
        read_genomic_strand : Strand
            The genomic strand the read aligned to (POS or NEG).
        gene_strand : Strand
            The annotated strand of the candidate gene.

        Returns
        -------
        float
            P(read came from gene on gene_strand | observed alignment strand)
        """
        if gene_strand == Strand.NONE or read_genomic_strand == Strand.NONE:
            # Unstranded gene or unstranded read → strand uninformative
            return 0.5

        if self.protocol == ProtocolType.UNSTRANDED:
            return 0.5

        # Determine if the alignment is "concordant" with this gene strand
        # under the detected protocol
        concordant = self._is_concordant_with_gene(
            read_genomic_strand, gene_strand
        )
        return self.p_sense if concordant else self.p_antisense

    def _is_concordant_with_gene(self, read_genomic_strand: Strand, 
                                  gene_strand: Strand) -> bool:
        """
        Check if the read's genomic alignment strand is concordant
        (sense) with the given gene strand under the detected protocol.

        For RF: read1 antisense to gene → concordant means 
                read_genomic_strand != gene_strand for read1.
        For FR: read1 sense to gene → concordant means 
                read_genomic_strand == gene_strand for read1.

        Since we collapse read1/read2 into a single genomic strand
        during fragment merging (using the protocol to determine the
        inferred transcript strand), concordance simply means the 
        inferred strand matches the gene strand.
        """
        # After protocol-aware strand inference, the fragment's strand
        # represents the inferred transcript strand. Concordance = match.
        return read_genomic_strand == gene_strand


class StrandModel:
    """
    First-pass strand protocol learner.

    Examines spliced reads (those with XS tags from STAR) to learn:
    1. Which protocol orientation (FR vs RF) is in use
    2. The strand-specificity probability p_sense as a Beta posterior

    Only spliced reads are used because the XS tag provides ground-truth
    for the gene strand at splice junctions.
    """

    # Weakly informative Beta(1,1) = Uniform prior
    PRIOR_ALPHA = 1.0
    PRIOR_BETA = 1.0

    def __init__(self):
        # Raw counts: (r1_matches_sj, r2_matches_sj) → count
        # r1_matches_sj: True if read1's genomic strand == splice junction strand
        # r2_matches_sj: True if read2's genomic strand == splice junction strand
        self._orientation_counts = {
            (True, True): 0,    # FF: both reads same strand as SJ
            (True, False): 0,   # FR: read1 same, read2 opposite
            (False, True): 0,   # RF: read1 opposite, read2 same
            (False, False): 0,  # RR: both reads opposite strand from SJ
        }
        self._n_fragments = 0
        self._n_spliced = 0
        self._finalized = False
        self._result: Optional[StrandModelResult] = None

    def observe(self, r1, r2):
        """
        Observe a paired-end fragment during the first pass.
        
        Only fragments where at least one read has the XS tag 
        (splice junction with unambiguous strand) are informative.

        Parameters
        ----------
        r1 : pysam.AlignedSegment or None
            Read 1 of the pair.
        r2 : pysam.AlignedSegment or None
            Read 2 of the pair.
        """
        self._n_fragments += 1

        # Determine splice junction strand from XS tag
        sj_strand = Strand.NONE
        for read in (r1, r2):
            if read is not None and read.has_tag('XS'):
                sj_strand |= Strand.from_str(read.get_tag('XS'))

        # Only use reads with unambiguous splice junction strand
        # (NONE < sj_strand < AMBIGUOUS means exactly POS or NEG)
        if not (Strand.NONE < sj_strand < Strand.AMBIGUOUS):
            return

        self._n_spliced += 1

        # Compare each read's genomic alignment strand to the SJ strand
        r1_strand = Strand.from_is_reverse(r1.is_reverse) if r1 else Strand.NONE
        r2_strand = Strand.from_is_reverse(r2.is_reverse) if r2 else Strand.NONE

        key = (r1_strand == sj_strand, r2_strand == sj_strand)
        self._orientation_counts[key] += 1

    def finalize(self, min_spliced_reads: int = 100) -> StrandModelResult:
        """
        Finalize the model after the first pass.

        Determines the protocol type (FR/RF/unstranded) and computes
        the Beta posterior for strand specificity.

        Parameters
        ----------
        min_spliced_reads : int
            Minimum number of spliced reads required for inference.
            Below this threshold, defaults to unstranded.

        Returns
        -------
        StrandModelResult
        """
        total_spliced = sum(self._orientation_counts.values())
        logger.info(f'Strand model: {self._n_fragments} total fragments, '
                     f'{self._n_spliced} spliced, '
                     f'{total_spliced} with unambiguous SJ strand')
        logger.info(f'Orientation counts: {self._orientation_counts}')

        if total_spliced < min_spliced_reads:
            logger.warning(
                f'Only {total_spliced} spliced reads (< {min_spliced_reads}). '
                f'Defaulting to unstranded.'
            )
            self._result = StrandModelResult(
                protocol=ProtocolType.UNSTRANDED,
                alpha=self.PRIOR_ALPHA,
                beta=self.PRIOR_BETA,
            )
            self._finalized = True
            return self._result

        # Determine dominant protocol by comparing FR vs RF counts
        # FR: read1 sense (matches SJ), read2 antisense
        n_fr = self._orientation_counts[(True, False)]
        # RF: read1 antisense, read2 sense (matches SJ)
        n_rf = self._orientation_counts[(False, True)]
        # FF and RR are unusual but tracked
        n_ff = self._orientation_counts[(True, True)]
        n_rr = self._orientation_counts[(False, False)]

        logger.info(f'Protocol counts: FR={n_fr}, RF={n_rf}, FF={n_ff}, RR={n_rr}')

        # The dominant protocol has the most counts
        if n_fr >= n_rf:
            protocol = ProtocolType.FR
            n_concordant = n_fr
            n_discordant = n_rf
        else:
            protocol = ProtocolType.RF
            n_concordant = n_rf
            n_discordant = n_fr

        # Check if there's enough evidence for strandedness
        # If concordant fraction < 0.55, treat as unstranded
        frac_concordant = n_concordant / (n_concordant + n_discordant) if (n_concordant + n_discordant) > 0 else 0.5
        if frac_concordant < 0.55:
            logger.warning(
                f'Concordant fraction {frac_concordant:.3f} too low. '
                f'Treating as unstranded.'
            )
            protocol = ProtocolType.UNSTRANDED

        # Build Beta posterior
        # Prior: Beta(1, 1) = Uniform
        # Likelihood: Binomial with n_concordant successes, n_discordant failures
        alpha = self.PRIOR_ALPHA + n_concordant
        beta = self.PRIOR_BETA + n_discordant

        self._result = StrandModelResult(
            protocol=protocol,
            alpha=alpha,
            beta=beta,
            n_concordant=n_concordant,
            n_discordant=n_discordant,
        )

        logger.info(
            f'Strand model result: protocol={protocol.value}, '
            f'p_sense={self._result.p_sense:.4f}, '
            f'specificity={self._result.strand_specificity:.4f}, '
            f'n_concordant={n_concordant}, n_discordant={n_discordant}'
        )

        self._finalized = True
        return self._result

    @property
    def result(self) -> StrandModelResult:
        if not self._finalized:
            raise RuntimeError('StrandModel not finalized. Call finalize() first.')
        return self._result

    @property
    def orientation_counts(self) -> dict:
        return dict(self._orientation_counts)
