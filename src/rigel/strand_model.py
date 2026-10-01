"""
rigel.strand_model — the library strand model, learned from annotated splice junctions.

Learns the strand distribution of a paired-end RNA-seq library by observing how fragment
alignment strands relate to annotated splice junction (SJ) strands.  A spliced fragment's
genomic GT/AG motif gives its true strand *independently of library prep* (STAR reports it in
the ``XS``/``ts`` tag), so comparing motif strand to aligner orientation over **annotated**
spliced fragments measures library-prep strand efficiency — once the junctions that are not genuine
RNA are separated out: a splice artifact (misaligned gDNA) reads ½ in a stranded library and a
reversed junction reads the opposite side (:func:`genuine_sense_fraction`).

After the R2 strand flip in the BAM scanner, the exon alignment strand effectively represents
read 1's genomic orientation, so the model's estimand is

    p_r1_sense = P(align_strand == reference strand)

Two nested views of ONE population, with one source of truth.

* :class:`SJStrandTable` — sense / antisense counts **per sj**, keyed on
  ``(ref, start, end, motif strand)``.  This is the primary record.
* :class:`StrandModel` — the 2×2 contingency table of ``align_strand × reference strand``.
  For the spliced model it is **exactly the table's marginal** and is built from it
  (:meth:`StrandModel.from_sj_table`), never accumulated separately.

The per-sj refinement keeps what the 2×2 loses: how the sense split varies ACROSS sj. The
accumulator's boundary spliced channels are not that population — they also pool unannotated and
implicit splices.

Models are immutable: built once from the scanner's arrays, then read.  There is no
observe/finalize lifecycle and therefore no way to score against a half-trained model.

This module does not serialize itself.  ``summary.json``'s ``strand_model`` block is hand-built
in :mod:`rigel.cli` from a handful of properties, plus :meth:`SJStrandTable.to_dict`.
"""

import logging
from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

from .types import Strand

logger = logging.getLogger(__name__)


def _class_shares(f: np.ndarray, c: np.ndarray) -> tuple[float, np.ndarray]:
    """``max_w Σ_g c_g·log(Σ_k w_k f_kg)`` over the 3-simplex, and the maximising shares. The objective is concave in
    ``w``, so the KKT conditions are necessary and sufficient: at the maximum ``Σ_g c_g f_kg/(w·f)_g`` equals ``Σc`` for
    every share in use and is at most ``Σc`` for every share at 0. A vertex, then an edge, that meets them is the
    maximum; only an interior maximum needs the log-barrier ascent, whose barrier keeps a share the maximum needs at 10⁻⁶.
    ``f`` is (3, G), each column scaled by its own constant (which only shifts the objective). A candidate on which some
    junction has no probability scores −∞, so the arithmetic warnings it raises are silenced, never acted on."""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        return _class_shares_unguarded(f, c)


#: The log-barrier's weight starts at this share of the junction count and falls tenfold per stage to an absolute
#: floor, where the barrier's own gap (3·μ for three shares) is below the log-likelihood's rounding — numeric.
_BARRIER_START, _BARRIER_FLOOR = 1e-2, 1e-11


def _class_shares_unguarded(f: np.ndarray, c: np.ndarray) -> tuple[float, np.ndarray]:
    from scipy.optimize import brentq

    def value(w):
        return float(c @ np.log(w @ f))

    total = float(c.sum())
    slack = 1e-12 * total  # rounding in the multipliers' sums, numeric only

    def satisfies(w, unused):
        mult = f @ (c / (w @ f))
        return bool(np.all(np.isfinite(mult[unused])) and np.all(mult[unused] <= total + slack))

    candidates = []
    for k in range(3):
        vertex = np.eye(3)[k]
        if satisfies(vertex, [j for j in range(3) if j != k]):
            return value(vertex), vertex
        candidates.append(vertex)
    for a_, b_, third in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):  # an edge t·e_a + (1 − t)·e_b
        dd = f[a_] - f[b_]

        def edge_slope(t, dd=dd, a_=a_, b_=b_):
            return float(c @ (dd / (t * f[a_] + (1.0 - t) * f[b_])))

        if edge_slope(1e-15) > 0.0 > edge_slope(1.0 - 1e-15):
            t = brentq(edge_slope, 1e-15, 1.0 - 1e-15, xtol=1e-15)
            e = np.zeros(3)
            e[a_], e[b_] = t, 1.0 - t
            if satisfies(e, [third]):
                return value(e), e
            candidates.append(e)
    d = f[:2] - f[2]  # the interior: free coordinates (w0, w1), w2 = 1 − w0 − w1
    w = np.full(3, 1.0 / 3.0)
    mu = _BARRIER_START * total
    while True:

        def phi(v, mu=mu):
            return value(v) + mu * float(np.log(v).sum())

        for _ in range(100):  # Newton on the barrier objective, kept strictly inside; numeric bound
            m = w @ f
            q = d / m
            g = q @ c + mu * (1.0 / w[:2] - 1.0 / w[2])
            h = -(q * c) @ q.T - mu * (np.diag(1.0 / w[:2] ** 2) + 1.0 / w[2] ** 2)
            step = -np.linalg.solve(h, g)
            dw = np.array([step[0], step[1], -step[0] - step[1]])
            shrink = dw < 0.0
            t = min(1.0, 0.99 * float(np.min(-w[shrink] / dw[shrink]))) if shrink.any() else 1.0
            now = phi(w)
            while t > 1e-16 and phi(w + t * dw) < now:
                t *= 0.5
            if t <= 1e-16:
                break
            w = w + t * dw
            if abs(float(g @ step)) * t < 1e-13 * max(1.0, abs(now)):
                break
        if mu <= _BARRIER_FLOOR:
            break
        mu = max(mu * 0.1, _BARRIER_FLOOR)
    candidates.append(w)
    best = max(candidates, key=value)
    return value(best), best


#: κ's search range, (0, ½) less the floor at each end: within it of 0 or of ½, κ̂ × (any library's RNA reads) moves
#: the fit by ≪ 1 read — numeric. The scan is uniform in log(κ/(½ − κ)), so it resolves the profile near ½ as finely as
#: near 0; it locates the basin before Brent's search — numeric.
_KAPPA_FLOOR, _SCAN_NODES = 1e-12, 129


def _genuine_mixture(
    k: np.ndarray, n: np.ndarray, c: np.ndarray
) -> tuple[float, np.ndarray, float, float, float]:
    """The three-class binomial mixture's maximum likelihood. Junction group g (``c_g`` junctions with ``k_g`` of ``n_g``
    reads on the minority side) is genuine RNA (rate κ), a splice artifact — misaligned gDNA (½) — or reversed (1 − κ).
    The class shares are solved exactly at each κ (:func:`_class_shares`); κ maximises that profile, by a scan over log(κ/(½ − κ))
    and Brent's method in its best cell — bounded work, where EM crawls without end as κ nears ½ (the classes merge) and
    can stop at a saddle. Returns ``(κ, shares, minority, reads, log-likelihood)``, ``minority`` and ``reads`` being
    the RNA classes' expected counts at the maximum."""
    from scipy.optimize import minimize_scalar

    half = n * np.log(0.5)

    def profile(u):
        kappa = 0.5 / (1.0 + float(np.exp(-u)))  # u = log(κ/(½ − κ))
        lp = np.stack(
            [
                k * np.log(kappa) + (n - k) * np.log1p(-kappa),
                half,
                k * np.log1p(-kappa) + (n - k) * np.log(kappa),
            ]
        )
        top = lp.max(axis=0)
        f = np.exp(lp - top)
        value, shares = _class_shares(f, c)
        return value + float(c @ top), shares, f

    end = float(np.log((0.5 - _KAPPA_FLOOR) / _KAPPA_FLOOR))
    grid = np.linspace(-end, end, _SCAN_NODES)
    scanned = [profile(u)[0] for u in grid]
    i = int(np.argmax(scanned))
    found = minimize_scalar(
        lambda u: -profile(u)[0],
        bounds=(grid[max(i - 1, 0)], grid[min(i + 1, grid.size - 1)]),
        method="bounded",
        options={"xatol": 1e-12},
    )
    u = float(found.x) if -found.fun >= scanned[i] else float(grid[i])
    loglik, shares, f = profile(u)
    with np.errstate(divide="ignore", invalid="ignore"):
        resp = shares[:, None] * f / (shares @ f)
    minority = float(c @ (resp[0] * k + resp[2] * (n - k)))
    reads = float(c @ ((resp[0] + resp[2]) * n))
    return 0.5 / (1.0 + float(np.exp(-u))), shares, minority, reads, loglik


def genuine_sense_fraction(n_sense, n_total) -> float:
    """The genuine junctions' sense fraction κ: the shipped ``Beta(1, 1)`` posterior mean, ``(minority + 1)/(reads +
    2)``, over the reads the three-class mixture credits to RNA at its maximum (:func:`_genuine_mixture`). Where every
    junction is genuine that is exactly ``(n_same + 1)/(n_obs + 2)``; it is never 0. The table is oriented so the
    genuine class's minority rate lies below ½. Gates: ``tests/test_strand_model.py``'s ``TestGenuineKappa``."""
    s = np.asarray(n_sense, dtype=np.float64)
    n = np.asarray(n_total, dtype=np.float64)
    keep = n > 0
    s, n = s[keep], n[keep]
    flip = s.sum() > 0.5 * n.sum()
    k = n - s if flip else s
    kn, cnt = np.unique(np.stack([k, n], axis=1), axis=0, return_counts=True)
    _, _, minority, reads, _ = _genuine_mixture(kn[:, 0], kn[:, 1], cnt.astype(np.float64))
    kappa = (minority + 1.0) / (reads + 2.0)
    return 1.0 - kappa if flip else kappa


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
      (:attr:`StrandModel.p_r1_sense`, then ``calibration.strand_balance.fit_strand_balance``).
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

    The spliced model's qualification (applied in the C++ scanner, not here): annotated splice
    junction, unique mapper, unambiguous exon strand, unambiguous SJ strand, non-chimeric.
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
        """The invariant, made executable: the 2×2 IS the sj table's marginal.

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

    @cached_property
    def genuine_n_same(self) -> float:
        """The sense count of the genuine junctions: ``n_same`` with the splice artifacts and reversed junctions taken
        out, scaled so that ``(genuine_n_same + 1)/(n_obs + 2)`` is :func:`genuine_sense_fraction`. ``n_same`` itself
        where there is no per-sj table, no observation, or no live strand channel — the gate reads the pooled counts
        (``calibration.region_init.strand_discriminability``), and removing junctions at ½ or at the opposite side
        only moves κ away from ½, so this can never switch a live channel off."""
        n = self.n_observations
        if self.sj_table is None or n == 0:
            return float(self.n_same)
        from .calibration.region_init import strand_discriminability

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
        """95% confidence interval for p_r1_sense (Wald interval)."""
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

    def strand_specificity_ci_epsilon(self, confidence: float = 0.99) -> float:
        """Upper credible limit on the minor-orientation rate ``1 − strand_specificity``.

        The ``confidence`` quantile of the Beta(n_minor + 1, n − n_minor + 1) posterior, clamped
        to [0, 0.5]; 0.5 when there are no observations.

        QC only. Its one consumer is the ``[CAL] Strand trainer`` log line in
        :func:`rigel.pipeline.run_pipeline`.
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
        # Clamp to [0, 0.5]: ss = max(p, 1-p) ≥ 0.5 so ε_CI ≤ 0.5.
        return max(0.0, min(0.5, r_ucl))


# ======================================================================
# StrandModels — single-model container with one diagnostic sub-model
# ======================================================================


@dataclass(frozen=True)
class StrandModels:
    """Container for the single RNA strand model plus one diagnostic sub-model.

    The primary strand model (``exonic_spliced``) is trained from uniquely mapped, non-chimeric
    SPLICED_ANNOT fragments whose candidate transcripts share one strand and whose alignment and
    SJ strands are each a single strand.  Annotated splice junctions prove RNA origin,
    making this an uncontaminated measure of library strand specificity.
    Probabilities are pure MLE from observed counts, and its 2×2 is the marginal of the
    per-sj :class:`SJStrandTable` it carries.

    One additional sub-model is retained **for diagnostics only** and
    is never used for scoring:

    * **exonic** — trained from every unique-mapper, non-chimeric, unambiguous-strand fragment
      that RESOLVES TO A TRANSCRIPT, spliced or not. Comparing its specificity to
      ``exonic_spliced`` reveals gDNA contamination (``contamination_gap`` in the CLI summary):
      unspliced genic fragments include gDNA, which is unstranded, so the mixed estimate is
      dragged toward ½.
      It is not "all exonic fragments" — intergenic fragments have no transcript and never enter
      it, so the gap it measures is GENIC contamination only.

    A fragment's gDNA likelihood uses a fixed strand probability of one half, not this container.
    """

    exonic_spliced: StrandModel = field(default_factory=StrandModel)

    # Diagnostic sub-model (not used for scoring)
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
            exonic=StrandModel.from_labels(strand_dict["exonic_obs"], strand_dict["exonic_truth"]),
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
