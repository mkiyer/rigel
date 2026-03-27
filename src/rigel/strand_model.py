"""
rigel.strand_model — Bayesian strand models learned from RNA-seq data.

Learns the strand distribution of a paired-end RNA-seq library by observing
how fragment alignment strands relate to annotated splice junction (SJ)
strands.  The SJ strand (from the STAR ``XS`` tag) provides ground truth
for the gene strand at each junction.

After the R2 strand flip in the BAM scanner, the exon alignment
strand effectively represents read 1's genomic orientation.  Comparing
this to the SJ strand tells us whether read 1 aligns in the *same*
direction as the gene (R1-sense) or the *opposite* direction (R1-antisense).

The model stores a single Beta posterior over:

    p_r1_sense = P(exon_strand == gene_strand)

and provides:

* ``strand_likelihood(exon_strand, gene_strand)`` for Bayesian quantification
* ``to_dict()`` / ``write_json()`` for human-readable output including
  derived protocol flags for downstream tools.
"""

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .types import Strand

logger = logging.getLogger(__name__)

#: Minimum observations needed to report a 95% credible interval.
#: Below this threshold, the posterior is too diffuse to be useful.
_MIN_CI_OBSERVATIONS: int = 10


@dataclass
class StrandModel:
    """Strand model learned from high-quality spliced reads.

    Accumulates a 2×2 contingency table of alignment strand
    (``exon_strand``) × SJ reference strand (``sj_strand``) from
    qualified spliced fragments.  Probabilities are pure MLE
    from these counts, with a safe fallback to 0.5 when no
    observations are available.

    Qualification criteria (applied by the caller, not this class)
    ---------------------------------------------------------------
    Only fragments meeting *all* of these should be passed to
    :meth:`observe`:

    1. Has annotated splice junction match(es).
    2. SJ merge resolves to a unique gene.
    3. Exon strand is unambiguous (POS or NEG).
    4. SJ strand is unambiguous (POS or NEG).
    """

    # --- 2×2 raw counts ---
    pos_pos: int = 0  # exon POS, SJ POS
    pos_neg: int = 0  # exon POS, SJ NEG
    neg_pos: int = 0  # exon NEG, SJ POS
    neg_neg: int = 0  # exon NEG, SJ NEG

    # --- Cached probabilities (set by finalize()) ---
    _cached_p_sense: float = field(default=0.0, repr=False, compare=False)
    _cached_p_antisense: float = field(default=0.0, repr=False, compare=False)
    _finalized: bool = field(default=False, repr=False, compare=False)

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

    def observe_batch(
        self,
        exon_strands: "np.ndarray",
        sj_strands: "np.ndarray",
    ) -> None:
        """Record a batch of strand observations (vectorized).

        Parameters
        ----------
        exon_strands : np.ndarray
            Integer array of exon strand values (1=POS, 2=NEG).
        sj_strands : np.ndarray
            Integer array of SJ strand values (1=POS, 2=NEG).
        """
        import numpy as _np

        exon = _np.asarray(exon_strands)
        sj = _np.asarray(sj_strands)
        e_pos = exon == 1  # Strand.POS
        e_neg = ~e_pos
        s_pos = sj == 1  # Strand.POS
        s_neg = ~s_pos
        self.pos_pos += int(_np.count_nonzero(e_pos & s_pos))
        self.pos_neg += int(_np.count_nonzero(e_pos & s_neg))
        self.neg_pos += int(_np.count_nonzero(e_neg & s_pos))
        self.neg_neg += int(_np.count_nonzero(e_neg & s_neg))

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
    def p_r1_sense(self) -> float:
        """MLE P(read 1 aligns in gene-sense direction).

        High (≈ 0.95) for R1-sense libraries (e.g. KAPA Stranded).
        Low  (≈ 0.05) for R1-antisense libraries (e.g. Illumina TruSeq dUTP).
        Near 0.50 for weakly-stranded libraries.
        Returns 0.5 (uninformative) when no observations are available.
        """
        n = self.n_observations
        if n == 0:
            return 0.5
        return self.n_same / n

    @property
    def p_r1_antisense(self) -> float:
        """Complement: P(read 1 aligns opposite to gene strand)."""
        return 1.0 - self.p_r1_sense

    @property
    def strand_specificity(self) -> float:
        """How strand-specific is the library?

        1.0 = perfect, 0.5 = no strand information.
        Equals ``max(p_r1_sense, p_r1_antisense)``.
        """
        return max(self.p_r1_sense, self.p_r1_antisense)

    @property
    def read1_sense(self) -> bool:
        """True if read 1 is predominantly sense (R1-sense protocol)."""
        return self.p_r1_sense >= 0.5

    def posterior_variance(self) -> float:
        """Variance of p_r1_sense using binomial variance."""
        n = self.n_observations
        if n == 0:
            return 0.25  # max variance when unknown
        p = self.n_same / n
        return (p * (1.0 - p)) / n

    def posterior_95ci(self) -> tuple[float, float]:
        """95% confidence interval for p_r1_sense (Wald interval)."""
        import math

        n = self.n_observations
        if n == 0:
            return (0.0, 1.0)
        p = self.n_same / n
        se = math.sqrt(p * (1.0 - p) / n)
        z = 1.959964
        lo = max(0.0, p - z * se)
        hi = min(1.0, p + z * se)
        return (lo, hi)

    # ------------------------------------------------------------------
    # Finalization (call after training, before scoring)
    # ------------------------------------------------------------------

    def finalize(self) -> None:
        """Cache derived probabilities for fast scoring.

        Must be called after all ``observe()`` calls are complete and
        before any ``strand_likelihood()`` calls during EM scoring.
        Uses pure MLE; falls back to 0.5 with zero observations.
        """
        n = self.n_observations
        if n == 0:
            self._cached_p_sense = 0.5
        else:
            self._cached_p_sense = (self.pos_pos + self.neg_neg) / n
        self._cached_p_antisense = 1.0 - self._cached_p_sense
        self._finalized = True

    def strand_likelihood(self, exon_strand: int, gene_strand: int) -> float:
        """Strand likelihood: P(exon_strand | fragment from gene_strand).

        Used during Bayesian quantification to weight candidate genes.
        Accepts int strand values (``Strand`` is an IntEnum so can be
        passed directly).

        * If ``exon_strand == gene_strand`` → ``p_r1_sense``
        * If ``exon_strand != gene_strand`` → ``1 - p_r1_sense``
        * If either is NONE/AMBIGUOUS (not 1 or 2) → 0.5 (uninformative)

        Must be called after :meth:`finalize`.
        """
        # 1=POS, 2=NEG are the only informative values
        if exon_strand != 1 and exon_strand != 2:
            return 0.5
        if gene_strand != 1 and gene_strand != 2:
            return 0.5
        if exon_strand == gene_strand:
            return self._cached_p_sense
        return self._cached_p_antisense

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
            "estimate": {
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

        # Add 95% CI if enough observations
        if self.n_observations >= _MIN_CI_OBSERVATIONS:
            lo, hi = self.posterior_95ci()
            d["estimate"]["ci_95"] = [
                float(round(lo, 6)),
                float(round(hi, 6)),
            ]

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
# StrandModels — single-model container with diagnostic sub-models
# ======================================================================

#: Minimum spliced observations to consider the strand model
#: well-supported.  Below this threshold a warning is emitted.
_MIN_STRAND_OBS_WARNING: int = 20


@dataclass
class StrandModels:
    """Container for the single RNA strand model plus diagnostic sub-models.

    The primary strand model (``exonic_spliced``) is trained from
    SPLICED_ANNOT fragments with unique gene assignment and unambiguous
    exon/SJ strands.  Annotated splice junctions prove RNA origin,
    making this an uncontaminated measure of library strand specificity.
    Probabilities are pure MLE from observed counts.

    Two additional sub-models are retained **for diagnostics only** and
    are never used for scoring:

    * **exonic** — trained from ALL exonic fragments (RNA + gDNA
      mixture).  Comparing its specificity to ``exonic_spliced``
      reveals gDNA contamination in exonic regions.
    * **intergenic** — trained from intergenic fragments (~100% gDNA).
      Expected ~50/50 since gDNA has no strand bias.

    gDNA is scored with a fixed strand probability of **0.5**
    (no strand bias), not learned from intergenic data.
    """

    exonic_spliced: StrandModel = field(default_factory=StrandModel)

    # Diagnostic sub-models (not used for scoring)
    exonic: StrandModel = field(default_factory=StrandModel)
    intergenic: StrandModel = field(default_factory=StrandModel)

    # ------------------------------------------------------------------
    # Finalization (call after training, before scoring)
    # ------------------------------------------------------------------

    def finalize(self) -> None:
        """Cache derived probabilities for fast scoring.

        Finalizes all sub-models.  Uses pure MLE; falls back to
        0.5 (uninformative) when no spliced observations exist.
        """
        self.exonic_spliced.finalize()

        n_obs = self.exonic_spliced.n_observations
        if n_obs == 0:
            logger.warning(
                "No spliced strand observations — strand model is "
                "prior-only (p_r1_sense=0.5, no strand information). Is this "
                "stranded RNA-seq data?"
            )
        elif n_obs < _MIN_STRAND_OBS_WARNING:
            logger.warning(
                "Only %d spliced strand observations (< %d); "
                "strand estimates may be noisy (SS=%.4f)",
                n_obs,
                _MIN_STRAND_OBS_WARNING,
                self.exonic_spliced.strand_specificity,
            )

        # Diagnostic sub-models (finalize for reporting, never for scoring)
        self.exonic.finalize()
        self.intergenic.finalize()

    # ------------------------------------------------------------------
    # Delegation to the RNA strand model
    # ------------------------------------------------------------------

    @property
    def p_r1_sense(self) -> float:
        """RNA model's P(read 1 is sense)."""
        return self.exonic_spliced.p_r1_sense

    @property
    def strand_specificity(self) -> float:
        """RNA model's strand specificity."""
        return self.exonic_spliced.strand_specificity

    @property
    def read1_sense(self) -> bool:
        """True if R1-sense protocol (p_r1_sense ≥ 0.5)."""
        return self.exonic_spliced.read1_sense

    @property
    def n_observations(self) -> int:
        """RNA model's observation count."""
        return self.exonic_spliced.n_observations

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """JSON-serializable summary of strand models.

        The primary ``exonic_spliced`` model is used for scoring.
        ``exonic`` and ``intergenic`` are diagnostic only.
        """
        return {
            "exonic_spliced": self.exonic_spliced.to_dict(),
            "diagnostics": {
                "exonic": self.exonic.to_dict(),
                "intergenic": self.intergenic.to_dict(),
            },
        }

    def write_json(self, path: Path | str) -> None:
        """Write strand models to a JSON file.

        Parameters
        ----------
        path : Path or str
            Output JSON file path.
        """
        path = Path(path)
        d = self.to_dict()

        # Add 95% CI for exonic_spliced model if enough observations
        if self.exonic_spliced.n_observations >= _MIN_CI_OBSERVATIONS:
            lo, hi = self.exonic_spliced.posterior_95ci()
            d["exonic_spliced"]["estimate"]["ci_95"] = [
                float(round(lo, 6)),
                float(round(hi, 6)),
            ]

        with open(path, "w") as fh:
            json.dump({"strand_models": d}, fh, indent=2)

        logger.info(f"Wrote strand models to {path}")

    def log_summary(self) -> None:
        """Log a human-readable summary of the trained strand models."""
        logger.info("Strand models:")
        logger.info(
            f"  [exonic_spliced] (RNA — used for scoring)  "
            f"{self.exonic_spliced.n_observations:,} obs, "
            f"p_r1_sense={self.exonic_spliced.p_r1_sense:.4f}, "
            f"specificity={self.exonic_spliced.strand_specificity:.4f}"
        )
        logger.info(
            f"  [exonic] (diagnostic)  "
            f"{self.exonic.n_observations:,} obs, "
            f"p_r1_sense={self.exonic.p_r1_sense:.4f}, "
            f"specificity={self.exonic.strand_specificity:.4f}"
        )
        logger.info(
            f"  [intergenic] (diagnostic)  "
            f"{self.intergenic.n_observations:,} obs, "
            f"p_r1_sense={self.intergenic.p_r1_sense:.4f}, "
            f"specificity={self.intergenic.strand_specificity:.4f}"
        )
