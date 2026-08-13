"""
rigel.strand_model — the library strand model, learned from annotated splice junctions.

Learns the strand distribution of a paired-end RNA-seq library by observing how fragment
alignment strands relate to annotated splice junction (SJ) strands.  A spliced fragment's
genomic GT/AG motif gives its true strand *independently of library prep* (STAR reports it in
the ``XS``/``ts`` tag), so comparing motif strand to aligner orientation over **annotated**
spliced fragments measures library-prep strand efficiency with no gDNA contamination — the
qualification is enforced in C++ by ``ResolvedFragment::get_is_strand_qualified()``.

After the R2 strand flip in the BAM scanner, the exon alignment strand effectively represents
read 1's genomic orientation, so the model's estimand is

    p_r1_sense = P(align_strand == reference strand)

⭐ **Two nested views of ONE population, with one source of truth.**

* :class:`SJStrandTable` — sense / antisense counts **per sj**, keyed on
  ``(ref, start, end, motif strand)``.  This is the primary record.
* :class:`StrandModel` — the 2×2 contingency table of ``align_strand × reference strand``.
  For the spliced model it is **exactly the table's marginal** and is built from it
  (:meth:`StrandModel.from_sj_table`), never accumulated separately.

The refinement exists because a **dispersion across sj** cannot be recovered from the
2×2: the RNA strand Beta-Binomial's mean (κ) and its overdispersion must be estimated from the
same population, and the overdispersion previously came from the accumulator's boundary spliced
channels, which also pool unannotated and implicit splices.  See


Models are **immutable**: built once from the scanner's arrays, then read.  There is no
observe/finalize lifecycle and therefore no way to score against a half-trained model.

⚠ **This module does not serialize itself.**  ``summary.json``'s ``strand_model`` block is
hand-built in :mod:`rigel.cli` from a handful of properties (plus
:meth:`SJStrandTable.to_dict`); a parallel ``to_dict``/``write_json`` pair lived here for a long
time with no caller and was deleted 2026-07-28.
"""

import logging
from dataclasses import dataclass, field

import numpy as np

from .types import Strand

logger = logging.getLogger(__name__)

#: Minimum spliced observations to consider the strand model well-supported.
#: Below this threshold a warning is emitted at construction.
_MIN_STRAND_OBS_WARNING: int = 20


@dataclass(frozen=True, slots=True)
class SJStrandTable:
    """Per-sj sense / antisense counts — the primary strand record.

    Six parallel arrays, one row per splice junction, sorted by
    ``(ref_id, start, end, motif_strand)`` in C++ so the contents never depend on thread
    scheduling or hash order.  A sj is uniquely specified by
    ``(reference, start, end, genomic splice-motif strand)``.

    **sense** means the aligner's fragment orientation agrees with the motif strand;
    **antisense** means it does not.  Each strand-qualified fragment contributes exactly ONE
    observation, to its leftmost annotated sj — ``sj_strand`` is read from the BAM
    ``XS``/``ts`` tag and is one value per *fragment*, so all sj a fragment spans share
    a single sense bit and crediting them all would repeat one observation K times, inflating
    the very dispersion this table exists to measure honestly.

    Consumers
    ---------
    * the **mean** κ — via the derived 2×2's ``n_same / n_observations``
      (:attr:`StrandModel.p_r1_sense`, then ``calibration.strand_balance.fit_strand_balance``);
    * the **dispersion** — the Beta-Binomial spread of ``(n_sense_j | depth_j)`` across
      sj at mean κ (``calibration.gdna_strand.fit_rna_strand_from_sj_table``).
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

        Counts arrive as ``uint64`` and are narrowed to ``int64``: they are counts of
        fragments, so ``int64`` cannot overflow before the fragment count itself does, and
        ``uint64`` silently promotes to float in mixed numpy arithmetic downstream.
        """
        if "sj_n_sense" not in d:
            return cls.empty()
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

        Writing ``sense ≡ (align == motif)``: over motif-POS sj sense is ``pos_pos``
        and antisense is ``neg_pos``; over motif-NEG sj sense is ``neg_neg`` and
        antisense is ``pos_neg``.  This identity is the correctness argument for the whole
        refinement and holds exactly — same qualification branch, one observation per fragment.
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
        """SpliceJunction-depth quantiles — "how deep are the sj that carry the fit"."""
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
            # "How many sj are deep enough to see the minority strand" is a
            # first-class question about a library: at κ ≈ 0.002 a sj needs
            # hundreds of reads before one disagreeing read is even expected.
            "n_sj_depth_ge_100": int(np.count_nonzero(depth >= 100)),
            "n_sj_depth_ge_1000": int(np.count_nonzero(depth >= 1000)),
        }


@dataclass(frozen=True)
class StrandModel:
    """A 2×2 contingency table of alignment strand × reference strand.

    Probabilities are pure MLE from the counts, with a safe fallback to 0.5 when there are no
    observations.  Immutable: build it with :meth:`from_labels` or :meth:`from_sj_table` and
    read it.

    ``sj_table`` is present on the **spliced** model only, where the 2×2 is exactly its
    marginal (:meth:`from_sj_table`) — the four counters are never maintained independently of
    it.  The all-exonic diagnostic model has no sj and therefore no table.

    Qualification (applied in C++ by ``get_is_strand_qualified()``, not here): annotated splice
    sj, unique mapper, unambiguous exon strand, unambiguous SJ strand, non-chimeric.

    ⚠ Deliberately NOT ``slots=True``: development caches under ``_selfsolve_cache`` /
    ``_calib_cache`` hold pickled instances, and a slotted class cannot restore a
    ``__dict__``-based pickle state.  The field names are load-bearing for the same reason.
    """

    # --- 2×2 raw counts ---
    pos_pos: int = 0  # exon POS, SJ POS
    pos_neg: int = 0  # exon POS, SJ NEG
    neg_pos: int = 0  # exon NEG, SJ POS
    neg_neg: int = 0  # exon NEG, SJ NEG

    #: The per-sj refinement this 2×2 marginalizes (spliced model only; ``None``
    #: on the all-exonic diagnostic model, which has no sj identity).
    sj_table: SJStrandTable | None = None

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    @classmethod
    def from_labels(cls, align_strands, sj_strands) -> "StrandModel":
        """Build the 2×2 from parallel per-fragment strand-label arrays (1=POS, 2=NEG)."""
        exon = np.asarray(align_strands)
        sj = np.asarray(sj_strands)
        e_pos = exon == int(Strand.POS)
        s_pos = sj == int(Strand.POS)
        return cls(
            pos_pos=int(np.count_nonzero(e_pos & s_pos)),
            pos_neg=int(np.count_nonzero(e_pos & ~s_pos)),
            neg_pos=int(np.count_nonzero(~e_pos & s_pos)),
            neg_neg=int(np.count_nonzero(~e_pos & ~s_pos)),
        )

    @classmethod
    def from_sj_table(cls, table: SJStrandTable) -> "StrandModel":
        """Build from the per-sj table — the 2×2 is its marginal (one source of truth)."""
        pos_pos, pos_neg, neg_pos, neg_neg = table.contingency()
        return cls(
            pos_pos=pos_pos,
            pos_neg=pos_neg,
            neg_pos=neg_pos,
            neg_neg=neg_neg,
            sj_table=table,
        )

    def contingency_matches_table(self) -> bool:
        """⭐ The invariant, made executable: the 2×2 IS the sj table's marginal.

        Trivially ``True`` when there is no table (the all-exonic diagnostic model). Nothing in the
        production path can violate it — :meth:`from_sj_table` is the only way the pair is built —
        so this exists to be asserted in tests and diagnostics, not defended against at runtime.
        """
        if self.sj_table is None:
            return True
        return self.sj_table.contingency() == (
            self.pos_pos,
            self.pos_neg,
            self.neg_pos,
            self.neg_neg,
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

    def strand_specificity_ci_epsilon(self, confidence: float = 0.99) -> float:
        """Measurement-uncertainty floor on (1 − strand_specificity).

        Returns ``ε_CI = 1 − UCL(ss)`` where ``UCL`` is the one-sided
        upper credible limit on ``strand_specificity = max(p, 1-p)``
        under a Beta(k + 1, (n − k) + 1) posterior on the minor-orientation
        rate (Jeffreys-ish / Laplace prior).

        ⚠ **QC only.** It was written for a gDNA calibration mixture that no longer exists; its
        one remaining consumer is the ``[CAL] Strand trainer`` log line in
        :func:`rigel.pipeline.run_pipeline`. Kept because "how well is the protocol pinned" is a
        fair thing to report, but do not describe it as feeding the deconvolution — it does not.

        - ``n = 0``: returns 0.5 (maximally uncertain → caps LLR completely).
        - ``n_minor = 0`` (degenerate): closed-form exact upper limit
          ``UCL = 1 − (1−conf)^(1/(n+1))`` from Beta(1, n+1).
        - general case: scipy ``betaincinv`` evaluation.
        """
        from scipy.special import betaincinv

        n = self.n_observations
        if n == 0:
            return 0.5
        # ``n_minor`` is the observation count that argues *against* perfect strand
        # specificity.  Whichever side the trainer called "sense" is irrelevant — the
        # stray-minority rate is what sizes the uncertainty.
        alpha = self.n_minor + 1.0
        beta = (n - self.n_minor) + 1.0
        # UCL on minor-orientation rate r = 1 − ss at the given confidence.
        r_ucl = float(betaincinv(alpha, beta, confidence))
        # Clamp to (0, 0.5): ss = max(p, 1-p) ≥ 0.5 so ε_CI ≤ 0.5.
        return max(0.0, min(0.5, r_ucl))


# ======================================================================
# StrandModels — single-model container with diagnostic sub-models
# ======================================================================


@dataclass(frozen=True)
class StrandModels:
    """Container for the single RNA strand model plus diagnostic sub-models.

    The primary strand model (``exonic_spliced``) is trained from
    SPLICED_ANNOT fragments with unique gene assignment and unambiguous
    exon/SJ strands.  Annotated splice junctions prove RNA origin,
    making this an uncontaminated measure of library strand specificity.
    Probabilities are pure MLE from observed counts, and its 2×2 is the marginal of the
    per-sj :class:`SJStrandTable` it carries.

    One additional sub-model is retained **for diagnostics only** and
    is never used for scoring:

    * **exonic** — trained from every unique-mapper, non-chimeric, unambiguous-strand fragment
      that RESOLVES TO A TRANSCRIPT, spliced or not. Comparing its specificity to
      ``exonic_spliced`` reveals gDNA contamination (``contamination_gap`` in the CLI summary):
      unspliced genic fragments include gDNA and nascent RNA, which are unstranded relative to
      the transcript, so the mixed estimate is dragged toward ½.
      ⚠ It is **not** "all exonic fragments (RNA + gDNA mixture)" as it was long documented —
      intergenic fragments have no transcript and never enter it. The gap it measures is
      **genic** contamination only.

    gDNA is scored with a fixed strand probability of **0.5**
    (no strand bias), not learned from intergenic data.

    ⚠ Not ``slots=True`` — see the note on :class:`StrandModel`.
    """

    exonic_spliced: StrandModel = field(default_factory=StrandModel)

    # Diagnostic sub-models (not used for scoring)
    exonic: StrandModel = field(default_factory=StrandModel)

    @classmethod
    def from_scan(cls, strand_dict: dict) -> "StrandModels":
        """Build both sub-models from the C++ scanner's ``strand_observations`` dict.

        The spliced model comes from the per-sj table (its 2×2 is the marginal); the
        all-exonic diagnostic has no sj identity and comes from its label arrays.
        Emits the low-evidence warnings once, here, where the counts first exist.
        """
        models = cls(
            exonic_spliced=StrandModel.from_sj_table(SJStrandTable.from_arrays(strand_dict)),
            exonic=StrandModel.from_labels(
                strand_dict.get("exonic_obs", []), strand_dict.get("exonic_truth", [])
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
        elif n_obs < _MIN_STRAND_OBS_WARNING:
            logger.warning(
                "Only %d spliced strand observations (< %d); "
                "strand estimates may be noisy (SS=%.4f)",
                n_obs,
                _MIN_STRAND_OBS_WARNING,
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

    def strand_specificity_ci_epsilon(self, confidence: float = 0.99) -> float:
        """Delegate: ε_CI from the primary (exonic_spliced) strand model."""
        return self.exonic_spliced.strand_specificity_ci_epsilon(confidence)

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
        logger.info(
            f"    sj={table.n_sj:,} "
            f"(depth median={table.depth_quantiles((0.5,))[0]:,}, "
            f"≥100={int(np.count_nonzero(table.depth >= 100)):,}, "
            f"≥1000={int(np.count_nonzero(table.depth >= 1000)):,})"
        )
        logger.info(
            f"  [exonic] (diagnostic)  "
            f"{self.exonic.n_observations:,} obs, "
            f"p_r1_sense={self.exonic.p_r1_sense:.4f}, "
            f"specificity={self.exonic.strand_specificity:.4f}"
        )
