"""
hulkrna.strand_model — Bayesian strand models learned from RNA-seq data.

Learns the strand distribution of a paired-end RNA-seq library by observing
how fragment alignment strands relate to annotated splice junction (SJ)
strands.  The SJ strand (from the STAR ``XS`` tag) provides ground truth
for the gene strand at each junction.

After the R2 strand flip in ``Fragment.from_reads()``, the exon alignment
strand effectively represents read 1's genomic orientation.  Comparing
this to the SJ strand tells us whether read 1 aligns in the *same*
direction as the gene (sense/FR) or the *opposite* direction (antisense/RF).

The model stores a single Beta posterior over:

    p_r1_sense = P(exon_strand == gene_strand)

and provides:

* ``strand_likelihood(exon_strand, gene_strand)`` for Bayesian counting
* ``to_dict()`` / ``write_json()`` for human-readable output including
  derived protocol flags for downstream tools.
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

from .types import Strand

logger = logging.getLogger(__name__)


@dataclass
class StrandModel:
    """Bayesian strand model learned from high-quality spliced reads.

    Accumulates a 2×2 contingency table of alignment strand
    (``exon_strand``) × SJ reference strand (``sj_strand``) from
    qualified spliced fragments.  All probabilities are derived
    from these raw counts plus a Beta prior.

    Qualification criteria (applied by the caller, not this class)
    ---------------------------------------------------------------
    Only fragments meeting *all* of these should be passed to
    :meth:`observe`:

    1. Has annotated splice junction match(es).
    2. SJ merge resolves to a unique gene.
    3. Exon strand is unambiguous (POS or NEG).
    4. SJ strand is unambiguous (POS or NEG).

    Beta posterior
    --------------
    ``alpha = prior_alpha + n_same``
    ``beta  = prior_beta  + n_opposite``

    where *n_same* counts fragments where exon and SJ strands agree
    (evidence for read 1 sense / FR) and *n_opposite* counts
    disagreements (evidence for read 1 antisense / RF).
    """

    # --- 2×2 raw counts ---
    pos_pos: int = 0    # exon POS, SJ POS
    pos_neg: int = 0    # exon POS, SJ NEG
    neg_pos: int = 0    # exon NEG, SJ POS
    neg_neg: int = 0    # exon NEG, SJ NEG

    # --- Beta prior (uniform by default) ---
    prior_alpha: float = 1.0
    prior_beta: float = 1.0

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------

    def observe(self, exon_strand: Strand, sj_strand: Strand) -> None:
        """Record one strand observation.

        Parameters
        ----------
        exon_strand : Strand
            Combined alignment strand of the fragment's exon blocks
            (POS or NEG after R2 flip; ≈ read 1 alignment strand).
        sj_strand : Strand
            Annotated SJ reference strand from the XS tag (POS or NEG).
        """
        if exon_strand == Strand.POS:
            if sj_strand == Strand.POS:
                self.pos_pos += 1
            else:
                self.pos_neg += 1
        else:
            if sj_strand == Strand.POS:
                self.neg_pos += 1
            else:
                self.neg_neg += 1

    # ------------------------------------------------------------------
    # Derived counts
    # ------------------------------------------------------------------

    @property
    def n_same(self) -> int:
        """Fragments where exon strand == SJ strand (read 1 sense)."""
        return self.pos_pos + self.neg_neg

    @property
    def n_opposite(self) -> int:
        """Fragments where exon strand != SJ strand (read 1 antisense)."""
        return self.pos_neg + self.neg_pos

    @property
    def n_observations(self) -> int:
        """Total qualified observations."""
        return self.n_same + self.n_opposite

    # ------------------------------------------------------------------
    # Posterior
    # ------------------------------------------------------------------

    @property
    def alpha(self) -> float:
        """Posterior alpha (same-direction evidence + prior)."""
        return self.prior_alpha + self.n_same

    @property
    def beta(self) -> float:
        """Posterior beta (opposite-direction evidence + prior)."""
        return self.prior_beta + self.n_opposite

    @property
    def p_r1_sense(self) -> float:
        """Posterior mean P(read 1 aligns in gene-sense direction).

        High (≈ 0.95) for FR libraries (e.g. KAPA Stranded).
        Low  (≈ 0.05) for RF libraries (e.g. dUTP / Illumina TruSeq).
        Near 0.50 for unstranded libraries.
        """
        return self.alpha / (self.alpha + self.beta)

    @property
    def p_r1_antisense(self) -> float:
        """Complement: P(read 1 aligns opposite to gene strand)."""
        return 1.0 - self.p_r1_sense

    @property
    def strand_specificity(self) -> float:
        """How strand-specific is the library?

        1.0 = perfect, 0.5 = unstranded.
        Equals ``max(p_r1_sense, p_r1_antisense)``.
        """
        return max(self.p_r1_sense, self.p_r1_antisense)

    @property
    def read1_sense(self) -> bool:
        """True if read 1 is predominantly sense (FR protocol)."""
        return self.p_r1_sense >= 0.5

    def posterior_variance(self) -> float:
        """Variance of the Beta posterior."""
        a, b = self.alpha, self.beta
        return (a * b) / ((a + b) ** 2 * (a + b + 1))

    def posterior_95ci(self) -> tuple[float, float]:
        """95% credible interval for p_r1_sense."""
        from scipy.stats import beta as beta_dist
        return (
            beta_dist.ppf(0.025, self.alpha, self.beta),
            beta_dist.ppf(0.975, self.alpha, self.beta),
        )

    # ------------------------------------------------------------------
    # Strand likelihood for Bayesian counting
    # ------------------------------------------------------------------

    def strand_likelihood(
        self, exon_strand: Strand, gene_strand: Strand,
    ) -> float:
        """Strand likelihood: P(exon_strand | fragment from gene_strand).

        Used during Bayesian counting to weight candidate genes.

        * If ``exon_strand == gene_strand`` → ``p_r1_sense``
        * If ``exon_strand != gene_strand`` → ``1 - p_r1_sense``
        * If either is NONE or AMBIGUOUS → 0.5 (uninformative)

        Parameters
        ----------
        exon_strand : Strand
            Combined alignment strand of the fragment.
        gene_strand : Strand
            Annotated strand of the candidate gene.

        Returns
        -------
        float
            Probability in [0, 1].
        """
        if (
            exon_strand not in (Strand.POS, Strand.NEG)
            or gene_strand not in (Strand.POS, Strand.NEG)
        ):
            return 0.5
        if exon_strand == gene_strand:
            return self.p_r1_sense
        return self.p_r1_antisense

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """JSON/YAML-serializable summary of the strand model.

        Includes raw counts, posterior parameters, derived
        probabilities, and protocol flags for downstream tools.
        All values are native Python types (no numpy scalars).
        """
        return {
            "observations": {
                "total": int(self.n_observations),
                "pos_pos": int(self.pos_pos),
                "pos_neg": int(self.pos_neg),
                "neg_pos": int(self.neg_pos),
                "neg_neg": int(self.neg_neg),
                "n_same": int(self.n_same),
                "n_opposite": int(self.n_opposite),
            },
            "posterior": {
                "alpha": float(round(self.alpha, 4)),
                "beta": float(round(self.beta, 4)),
                "prior_alpha": float(self.prior_alpha),
                "prior_beta": float(self.prior_beta),
                "variance": float(round(self.posterior_variance(), 10)),
            },
            "probabilities": {
                "p_r1_sense": float(round(self.p_r1_sense, 6)),
                "p_r1_antisense": float(round(self.p_r1_antisense, 6)),
                "strand_specificity": float(round(self.strand_specificity, 6)),
            },
            "protocol": {
                "read1_sense": bool(self.read1_sense),
            },
        }

    def write_json(self, path: Path | str) -> None:
        """Write the strand model to a JSON file.

        Parameters
        ----------
        path : Path or str
            Output JSON file path.
        """
        path = Path(path)
        d = self.to_dict()

        # Add 95% CI if enough observations (requires scipy)
        if self.n_observations >= 10:
            try:
                lo, hi = self.posterior_95ci()
                d["posterior"]["ci_95"] = [
                    float(round(lo, 6)), float(round(hi, 6)),
                ]
            except ImportError:
                pass

        with open(path, "w") as fh:
            json.dump(
                {"strand_model": d},
                fh,
                indent=2,
            )

        logger.info(f"Wrote strand model to {path}")

    def log_summary(self) -> None:
        """Log a human-readable summary of the learned strand model."""
        total = self.n_observations
        logger.info(
            f"Strand model: {total:,} observations "
            f"(n_same={self.n_same:,}, n_opposite={self.n_opposite:,})"
        )
        logger.info(
            f"  2×2 table: POS→POS={self.pos_pos:,}, POS→NEG={self.pos_neg:,}, "
            f"NEG→POS={self.neg_pos:,}, NEG→NEG={self.neg_neg:,}"
        )
        logger.info(
            f"  p_r1_sense={self.p_r1_sense:.4f}, "
            f"strand_specificity={self.strand_specificity:.4f}, "
            f"read1_sense={self.read1_sense}"
        )


# ======================================================================
# StrandModels — multi-region container
# ======================================================================


@dataclass
class StrandModels:
    """Container for region-specific strand models.

    Holds four independently-trained :class:`StrandModel` instances
    organized by genomic region and purity:

    * **exonic_spliced** *(pure RNA)* — trained from SPLICED_ANNOT
      fragments with unique gene assignment and unambiguous exon/SJ
      strands.  Annotated splice junctions prove RNA origin, making
      this the gold-standard measure of library strand specificity.
      This is the primary model used for transcript scoring in EM.
    * **exonic** *(RNA + gDNA mixture)* — trained from ALL exonic
      fragments (SPLICED_ANNOT + SPLICED_UNANNOT + UNSPLICED) with
      unique gene assignment.  Gene annotation strand is truth.
      Dilution toward 50% relative to ``exonic_spliced`` indicates
      gDNA contamination in exonic regions.
    * **intronic** *(nascent RNA + gDNA mixture)* — trained from
      INTRON fragments with unique gene assignment.  Gene annotation
      strand is truth.  Dilution toward 50% relative to
      ``exonic_spliced`` indicates gDNA fraction in intronic regions.
    * **intergenic** *(~100% gDNA)* — trained from intergenic
      fragments.  Uses POS as a synthetic reference strand so
      ``n_same`` counts POS-aligned and ``n_opposite`` counts
      NEG-aligned.  Expected ~50/50 since gDNA is unstranded.

    The ``strand_likelihood``, ``p_r1_sense``, etc. convenience
    methods delegate to ``exonic_spliced`` (the pure RNA model).
    """

    exonic_spliced: StrandModel = field(default_factory=StrandModel)
    exonic: StrandModel = field(default_factory=StrandModel)
    intronic: StrandModel = field(default_factory=StrandModel)
    intergenic: StrandModel = field(default_factory=StrandModel)

    # ------------------------------------------------------------------
    # Category-aware model selection
    # ------------------------------------------------------------------

    def model_for_category(self, count_cat: int) -> StrandModel:
        """Select the strand model matching a fragment's category.

        * **SPLICED_ANNOT** → ``exonic_spliced`` (pure RNA, gold standard)
        * **SPLICED_UNANNOT / UNSPLICED** → ``exonic`` (RNA + gDNA mixture)
        * **INTRON** → ``intronic`` (nascent RNA + gDNA mixture)

        The category-specific model reflects the expected gDNA dilution
        in each genomic region: exonic_spliced has the highest strand
        specificity (annotated SJs prove RNA), exonic is slightly
        diluted, and intronic is more diluted.  Using the matching
        model ensures that strand-based likelihoods accurately reflect
        each region's RNA purity.

        Parameters
        ----------
        count_cat : int
            ``CountCategory`` value (0–3).

        Returns
        -------
        StrandModel
        """
        from .categories import CountCategory

        cat = CountCategory(count_cat)
        if cat == CountCategory.SPLICED_ANNOT:
            return self.exonic_spliced
        if cat == CountCategory.INTRON:
            return self.intronic
        # UNSPLICED and SPLICED_UNANNOT → exonic (all)
        return self.exonic

    # ------------------------------------------------------------------
    # Delegation to exonic_spliced model (pure RNA gold standard)
    # ------------------------------------------------------------------

    def strand_likelihood(
        self, exon_strand: Strand, gene_strand: Strand,
    ) -> float:
        """Delegate to the exonic_spliced strand model."""
        return self.exonic_spliced.strand_likelihood(exon_strand, gene_strand)

    @property
    def p_r1_sense(self) -> float:
        """Exonic spliced model's P(read 1 is sense)."""
        return self.exonic_spliced.p_r1_sense

    @property
    def strand_specificity(self) -> float:
        """Exonic spliced model's strand specificity."""
        return self.exonic_spliced.strand_specificity

    @property
    def read1_sense(self) -> bool:
        """Exonic spliced model's read1_sense flag."""
        return self.exonic_spliced.read1_sense

    @property
    def n_observations(self) -> int:
        """Exonic spliced model's observation count."""
        return self.exonic_spliced.n_observations

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """JSON-serializable summary of all four strand models."""
        return {
            "exonic_spliced": self.exonic_spliced.to_dict(),
            "exonic": self.exonic.to_dict(),
            "intronic": self.intronic.to_dict(),
            "intergenic": self.intergenic.to_dict(),
        }

    def write_json(self, path: Path | str) -> None:
        """Write all strand models to a JSON file.

        Parameters
        ----------
        path : Path or str
            Output JSON file path.
        """
        path = Path(path)
        d = self.to_dict()

        # Add 95% CI for exonic_spliced model if enough observations
        if self.exonic_spliced.n_observations >= 10:
            try:
                lo, hi = self.exonic_spliced.posterior_95ci()
                d["exonic_spliced"]["posterior"]["ci_95"] = [
                    float(round(lo, 6)), float(round(hi, 6)),
                ]
            except ImportError:
                pass

        with open(path, "w") as fh:
            json.dump({"strand_models": d}, fh, indent=2)

        logger.info(f"Wrote strand models to {path}")

    def log_summary(self) -> None:
        """Log a human-readable summary of all strand models."""
        logger.info("Strand models:")
        logger.info(f"  [exonic_spliced]  {self.exonic_spliced.n_observations:,} obs, "
                     f"p_r1_sense={self.exonic_spliced.p_r1_sense:.4f}, "
                     f"specificity={self.exonic_spliced.strand_specificity:.4f}")
        logger.info(f"  [exonic]  {self.exonic.n_observations:,} obs, "
                     f"p_r1_sense={self.exonic.p_r1_sense:.4f}, "
                     f"specificity={self.exonic.strand_specificity:.4f}")
        logger.info(f"  [intronic]  {self.intronic.n_observations:,} obs, "
                     f"p_r1_sense={self.intronic.p_r1_sense:.4f}, "
                     f"specificity={self.intronic.strand_specificity:.4f}")
        logger.info(f"  [intergenic]  {self.intergenic.n_observations:,} obs, "
                     f"p_r1_sense={self.intergenic.p_r1_sense:.4f}, "
                     f"specificity={self.intergenic.strand_specificity:.4f}")
