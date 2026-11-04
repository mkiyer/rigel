"""rigel.strand_model — the library strand model, learned from annotated splice junctions.

A spliced fragment's genomic splice motif gives its true strand independently of library prep (the aligner reports it
in ``XS`` / ``ts``), so comparing motif strand with alignment orientation over annotated spliced fragments measures the
library's strand specificity, once splice artifacts (misaligned gDNA, which reads ½ in a stranded library) and
reversed junctions are separated out (:func:`rigel.junction_fit.genuine_sense_fraction`). The scanner flips R2, so the
alignment strand is read 1's genomic orientation and the estimand is ``p_r1_sense = P(align_strand == reference
strand)``.

The primary record is per junction (:class:`SJStrandTable`), which keeps how the sense split varies across junctions.
The spliced model's 2×2 (:class:`StrandModel`) is that table's marginal and is built from it, never accumulated
separately. Models are immutable: built once from the scanner's observations, then read. ``summary.json``'s
``strand_model`` block is built in :mod:`rigel.cli`.
"""

import logging
from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

from .config import CONSTANTS
from .junction_fit import genuine_sense_fraction
from .types import Strand

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class SJStrandTable:
    """Per-sj sense / antisense counts, the primary strand record.

    Six parallel arrays, one row per junction ``(reference, start, end, motif strand)``, sorted in C++ so the contents
    never depend on thread scheduling. Sense means the alignment orientation agrees with the motif strand. Each
    strand-qualified fragment adds one observation, to its leftmost annotated sj: ``sj_strand`` is one value per
    fragment, so crediting every sj it spans would count one observation K times.

    Read by the junction fit (:mod:`rigel.junction_fit`, through :attr:`StrandModel.genuine_n_same`) and by
    ``summary.json``'s strand block.
    """

    ref_id: np.ndarray  # int32[n_sj]
    start: np.ndarray  # int64[n_sj]
    end: np.ndarray  # int64[n_sj]
    motif_strand: np.ndarray  # int8[n_sj] — Strand.POS / Strand.NEG
    n_sense: np.ndarray  # int64[n_sj] — aligner orientation agrees with the motif
    n_antisense: np.ndarray  # int64[n_sj]

    @classmethod
    def empty(cls) -> "SJStrandTable":
        """A table with no sj (an unspliced or unscanned library)."""
        return cls(
            ref_id=np.empty(0, dtype=np.int32),
            start=np.empty(0, dtype=np.int64),
            end=np.empty(0, dtype=np.int64),
            motif_strand=np.empty(0, dtype=np.int8),
            n_sense=np.empty(0, dtype=np.int64),
            n_antisense=np.empty(0, dtype=np.int64),
        )

    @classmethod
    def from_arrays(cls, d: dict) -> "SJStrandTable":
        """Build from the C++ ``strand_observations`` dict (``sj_*`` keys).

        Counts arrive as ``uint64`` and are narrowed to ``int64``, since ``uint64`` silently promotes to float in mixed
        numpy arithmetic.
        """
        return cls(
            ref_id=np.asarray(d["sj_ref_id"], dtype=np.int32),
            start=np.asarray(d["sj_start"], dtype=np.int64),
            end=np.asarray(d["sj_end"], dtype=np.int64),
            motif_strand=np.asarray(d["sj_motif_strand"], dtype=np.int8),
            n_sense=np.asarray(d["sj_n_sense"], dtype=np.int64),
            n_antisense=np.asarray(d["sj_n_antisense"], dtype=np.int64),
        )

    @property
    def n_sj(self) -> int:
        """Distinct splice junctions observed."""
        return int(self.n_sense.size)

    @property
    def depth(self) -> np.ndarray:
        """Per-sj qualified fragment count ``n_j = sense_j + antisense_j``."""
        return self.n_sense + self.n_antisense

    @property
    def n_observations(self) -> int:
        """Total qualified fragments — one per fragment, so this equals the 2×2's total."""
        return int(self.depth.sum())

    def contingency(self) -> tuple[int, int, int, int]:
        """The 2×2 ``(pos_pos, pos_neg, neg_pos, neg_neg)`` this table marginalizes to.

        With ``sense ≡ (align == motif)``: over motif-POS sj sense is ``pos_pos`` and antisense ``neg_pos``; over
        motif-NEG sj sense is ``neg_neg`` and antisense ``pos_neg``. The identity is exact: one qualification rule, one
        observation per fragment.
        """
        pos = self.motif_strand == int(Strand.POS)
        neg = self.motif_strand == int(Strand.NEG)
        return (
            int(self.n_sense[pos].sum()),  # pos_pos
            int(self.n_antisense[neg].sum()),  # pos_neg
            int(self.n_antisense[pos].sum()),  # neg_pos
            int(self.n_sense[neg].sum()),  # neg_neg
        )

    def depth_quantiles(self, qs: tuple[float, ...] = (0.5, 0.9, 0.99)) -> list[int]:
        """Per-sj depth quantiles — "how deep are the sj that carry the fit"."""
        if self.n_sj == 0:
            return [0] * len(qs)
        return [int(v) for v in np.quantile(self.depth, qs)]

    def to_dict(self) -> dict:
        """JSON-serializable QC summary: how much sj evidence this library carries."""
        depth = self.depth
        q50, q90, q99 = self.depth_quantiles()
        return {
            "n_sj": self.n_sj,
            "n_observations": self.n_observations,
            "depth_median": q50,
            "depth_p90": q90,
            "depth_p99": q99,
            "depth_max": int(depth.max()) if self.n_sj else 0,
        }


@dataclass(frozen=True)
class StrandModel:
    """A 2×2 contingency table of alignment strand × reference strand.

    Probabilities are the MLE from the counts, 0.5 when there are none. ``sj_table`` is present on the spliced model
    only, where the 2×2 is its marginal (:meth:`from_sj_table`); the all-exonic diagnostic model has no table. The
    spliced model's qualification is applied in the C++ scanner (``ResolvedFragment::is_strand_qualified``).
    """

    # --- 2×2 raw counts ---
    pos_pos: int = 0  # exon POS, SJ POS
    pos_neg: int = 0  # exon POS, SJ NEG
    neg_pos: int = 0  # exon NEG, SJ POS
    neg_neg: int = 0  # exon NEG, SJ NEG

    #: The per-sj table this 2×2 marginalizes; ``None`` on the all-exonic diagnostic model.
    sj_table: SJStrandTable | None = None

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_sj_table(cls, table: SJStrandTable) -> "StrandModel":
        """Build from the per-sj table, whose marginal the 2×2 is."""
        pos_pos, pos_neg, neg_pos, neg_neg = table.contingency()
        return cls(
            pos_pos=pos_pos,
            pos_neg=pos_neg,
            neg_pos=neg_pos,
            neg_neg=neg_neg,
            sj_table=table,
        )

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
    def n_minor(self) -> int:
        """Minor-orientation observations relative to the learned protocol."""
        return min(self.n_same, self.n_opposite)

    @property
    def n_observations(self) -> int:
        """Total qualified observations."""
        return self.n_same + self.n_opposite

    # ------------------------------------------------------------------
    # Posterior
    # ------------------------------------------------------------------

    @cached_property
    def genuine_n_same(self) -> float:
        """The sense count of the genuine junctions: ``n_same`` without the splice artifacts and reversed junctions,
        scaled so that ``(genuine_n_same + 1) / (n_obs + 2)`` is :func:`genuine_sense_fraction`.

        ``n_same`` itself where there is no per-sj table, no observation, or no live strand channel (judged on the pooled
        counts by ``calibration.strand_balance.strand_discriminability``). Removing junctions at ½ or on the opposite side
        only moves κ away from ½, so this never switches a live channel off."""
        n = self.n_observations
        if self.sj_table is None or n == 0:
            return float(self.n_same)
        from .calibration.strand_balance import strand_discriminability

        if strand_discriminability((self.n_same + 1.0) / (n + 2.0), n) <= 0.0:
            return float(self.n_same)
        kappa = genuine_sense_fraction(self.sj_table.n_sense, self.sj_table.depth)
        return min(max(kappa * (n + 2.0) - 1.0, 0.0), float(n))

    @property
    def p_r1_sense(self) -> float:
        """P(read 1 aligns in gene-sense direction), from the genuine junctions (:attr:`genuine_n_same`).

        High (≈ 0.95) for R1-sense libraries (e.g. KAPA Stranded).
        Low  (≈ 0.05) for R1-antisense libraries (e.g. Illumina TruSeq dUTP).
        Near 0.50 for weakly-stranded libraries.
        Returns 0.5 (uninformative) when no observations are available.
        """
        n = self.n_observations
        if n == 0:
            return 0.5
        return self.genuine_n_same / n

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
        p = self.p_r1_sense
        return (p * (1.0 - p)) / n

    def posterior_95ci(self) -> tuple[float, float]:
        """95% Wald interval for p_r1_sense; QC only (``summary.json``'s ``strand_model.ci_95``)."""
        import math

        from scipy.special import ndtri

        n = self.n_observations
        if n == 0:
            return (0.0, 1.0)
        p = self.p_r1_sense
        se = math.sqrt(p * (1.0 - p) / n)
        z = float(ndtri(0.975))  # the two-sided 95 % normal quantile
        lo = max(0.0, p - z * se)
        hi = min(1.0, p + z * se)
        return (lo, hi)

    def strand_specificity_ci_epsilon(self) -> float:
        """Upper credible limit on the pooled 2×2's minor-orientation rate ``n_minor / n``.

        The ``CONSTANTS.qc.strand_ci_confidence`` quantile of the Beta(n_minor + 1, n − n_minor + 1) posterior, clamped
        to [0, 0.5]; 0.5 when there are no observations. ``n_minor`` comes from the pooled table, artifacts included, while
        ``strand_specificity`` comes from the genuine junctions, so this bounds ``1 − strand_specificity`` only where
        that fit falls back to the pooled counts. QC only: the strand log line in :func:`rigel.pipeline.run_pipeline`.
        """
        from scipy.special import betaincinv

        n = self.n_observations
        if n == 0:
            return 0.5
        # n_minor argues against perfect specificity, whichever side is sense.
        alpha = self.n_minor + 1.0
        beta = (n - self.n_minor) + 1.0
        # UCL on the pooled minor-orientation rate r = n_minor / n at the given confidence.
        r_ucl = float(betaincinv(alpha, beta, CONSTANTS.qc.strand_ci_confidence))
        # Clamp to [0, 0.5]: n_minor = min(n_same, n_opposite), so the minor rate is at most 0.5.
        return max(0.0, min(0.5, r_ucl))


# ======================================================================
# StrandModels — single-model container with one diagnostic sub-model
# ======================================================================


@dataclass(frozen=True)
class StrandModels:
    """The RNA strand model plus one diagnostic sub-model.

    ``exonic_spliced`` is trained from uniquely mapped, annotated spliced fragments whose candidate transcripts share
    one strand and whose alignment and junction strands are each single. An annotated junction proves RNA, so this
    measures strand specificity free of gDNA. Its 2×2 is the marginal of the :class:`SJStrandTable` it carries.

    ``exonic`` is diagnostic only, never used for scoring: every unique-mapper, unambiguous-strand fragment that
    resolves to a transcript, spliced or not. Its unspliced fragments include unstranded gDNA, so its specificity sits
    nearer ½, and the gap to ``exonic_spliced`` (``contamination_gap`` in the CLI summary) measures genic
    contamination; intergenic fragments never enter it.

    A fragment's gDNA likelihood uses a fixed strand probability of one half, not this container.
    """

    exonic_spliced: StrandModel = field(default_factory=StrandModel)

    # diagnostic only
    exonic: StrandModel = field(default_factory=StrandModel)

    @classmethod
    def from_scan(cls, strand_dict: dict) -> "StrandModels":
        """Build both sub-models from the C++ scanner's ``strand_observations`` dict.

        The spliced model comes from the per-sj table; the scanner counts the diagnostic's 2×2 cell by cell. Warns
        here, once, when the evidence is thin.
        """
        models = cls(
            exonic_spliced=StrandModel.from_sj_table(SJStrandTable.from_arrays(strand_dict)),
            exonic=StrandModel(
                pos_pos=int(strand_dict["exonic_pos_pos"]),
                pos_neg=int(strand_dict["exonic_pos_neg"]),
                neg_pos=int(strand_dict["exonic_neg_pos"]),
                neg_neg=int(strand_dict["exonic_neg_neg"]),
            ),
        )
        models._warn_if_underpowered()
        return models

    def _warn_if_underpowered(self) -> None:
        """Warn when the spliced population is too thin to identify the strand protocol."""
        n_obs = self.exonic_spliced.n_observations
        if n_obs == 0:
            logger.warning(
                "No spliced strand observations — strand model is "
                "prior-only (p_r1_sense=0.5, no strand information). Is this "
                "stranded RNA-seq data?"
            )
        elif n_obs < CONSTANTS.qc.strand_min_observations:
            logger.warning(
                "Only %d spliced strand observations (< %d); "
                "strand estimates may be noisy (SS=%.4f)",
                n_obs,
                CONSTANTS.qc.strand_min_observations,
                self.exonic_spliced.strand_specificity,
            )

    # ------------------------------------------------------------------
    # Delegation to the RNA strand model
    # ------------------------------------------------------------------

    @property
    def sj_table(self) -> SJStrandTable:
        """The RNA model's per-sj table (empty when the library was never scanned)."""
        return self.exonic_spliced.sj_table or SJStrandTable.empty()

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

    def strand_specificity_ci_epsilon(self) -> float:
        """The RNA model's :meth:`StrandModel.strand_specificity_ci_epsilon`."""
        return self.exonic_spliced.strand_specificity_ci_epsilon()

    def log_summary(self) -> None:
        """Log a human-readable summary of the trained strand models."""
        table = self.sj_table
        logger.info("Strand models:")
        logger.info(
            f"  [exonic_spliced] (RNA — used for scoring)  "
            f"{self.exonic_spliced.n_observations:,} obs, "
            f"p_r1_sense={self.exonic_spliced.p_r1_sense:.4f}, "
            f"specificity={self.exonic_spliced.strand_specificity:.4f}"
        )
        logger.info(f"    sj={table.n_sj:,} (depth median={table.depth_quantiles((0.5,))[0]:,})")
        logger.info(
            f"  [exonic] (diagnostic)  "
            f"{self.exonic.n_observations:,} obs, "
            f"p_r1_sense={self.exonic.p_r1_sense:.4f}, "
            f"specificity={self.exonic.strand_specificity:.4f}"
        )
