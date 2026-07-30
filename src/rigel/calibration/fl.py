"""rigel.calibration.fl — gDNA / RNA fragment-length distributions (PR 4c).

Produces the library-wide FL distributions the calibrator's **effective lengths**
need (boundary ``fl_mean``, region ``E_f[max(0, L−ℓ)]``): the **gDNA FL** (from
gDNA-dominated regions/boundaries) and the **RNA FL** (spliced fragments), each
**smoothly empirical-Bayes-shrunk** toward the global FL.

This is **not** a per-fragment FL likelihood — that channel is deliberately
excluded from calibration. BOTH FLs drive per-node effective lengths in the sweep
(``bp_solver.build_node_geometry``): gDNA eff-lengths use the gDNA FL, RNA (nascent
unspliced + spliced) eff-lengths use the RNA FL.

⭐ **Every pool is PURE BY CONSTRUCTION** (`docs/ACCUMULATOR_DESIGN.md` §8), and that purity is what
removes the circularity: each model is fitted only from a population known to be one component, so
nothing is ever estimated from the fragments it will later explain.

* **gDNA** = the two *contained* pools, intergenic + intronic. ~99 % gDNA on real data.
* **RNA** = the annotated-junction pool, splice OBSERVED. gDNA cannot be spliced.
  ⚠ ``sj_implicit`` fragments are already excluded by the accumulator — a splice that was never
  sequenced is a product of the very model this pool is used to fit.
* **splash** = the two *_EXON crossing pools, exposed for QC and **deliberately not folded into the
  gDNA model**. They are the only ON-TARGET gDNA population, and on-target gDNA runs ~42 bp shorter
  than off-target (§8.2), so they land between the two pure means. ⭐ The shipped model summed four
  differently-tilted pools and read a gDNA mean of **146.05** where the pure intergenic pool says
  **88.0** — biased long by ~40 %, by pooling exactly these.

The pool axis itself lives in :mod:`rigel.scan_payload`, with the schema, because it is the
accumulator's own enum and a private copy here is how three files come to disagree about which row is
which.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..scan_payload import (
    POOL_DNA_INTERGENIC,
    POOL_DNA_INTERGENIC_EXON,
    POOL_DNA_INTRONIC,
    POOL_DNA_INTRON_EXON,
    POOL_RNA_SPLICED,
)

if TYPE_CHECKING:
    from ..frag_length_model import FragmentLengthModel
    from ..scan_payload import AccumulatorPayload

__all__ = [
    "FLModels",
    "POOL_EB_PRIOR_ESS",
    "build_fl_models",
    "gdna_fl_mass",
    "rna_fl_mass",
    "splash_fl_mass",
]

#: The pure gDNA pools: contained in an intergenic or intronic node. ⚠ CONTAINED ONLY — see the module
#: docstring for why the two crossing "splash" pools stay out.
_GDNA_POOLS = (POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC)

#: The pure RNA pool: an OBSERVED splice across an annotated junction.
_RNA_POOLS = (POOL_RNA_SPLICED,)

#: On-target gDNA, reported rather than fitted.
_SPLASH_POOLS = (POOL_DNA_INTRON_EXON, POOL_DNA_INTERGENIC_EXON)

#: Dirichlet pseudo-count for the smooth EB shrink toward the global FL. Not a
#: cliff: ``pool_total ≫ prior_ess`` → empirical; ``≪`` → the global anchor; ``= 0``
#: → global. Revisit the default on real-data pool sizes (PR04c decision 5).
POOL_EB_PRIOR_ESS: float = 1000.0


@dataclass(frozen=True, slots=True)
class FLModels:
    """Library-wide FL distributions (float64[max_size + 1]) + their pool totals.

    Two views of each distribution are carried:

    * ``*_pmf`` — the EB-smoothed pmf used for **scoring / calibration** (RNA and
      gDNA are shrunk toward ``global_pmf``; see :func:`build_fl_models`).
    * ``*_counts`` — the **raw, unsmoothed** histograms (aligned to
      ``max_size + 1`` bins) the pmfs were built from. These are the honest
      empirical distributions surfaced as QC (:meth:`rna_model` / :meth:`gdna_model`
      / :meth:`global_model`); the EB smoothing is a <1% perturbation at real
      library scale and stays internal to scoring.

    ``gdna_counts`` is the pure contained-pool histogram (:func:`gdna_fl_mass`) — intergenic + intronic,
    so gDNA *plus* whatever nascent RNA sits in an intron, not a deconvolved pure-gDNA distribution.
    ``rna_counts`` is the annotated-junction histogram (:func:`rna_fl_mass`).
    """

    global_pmf: np.ndarray  # unconditional anchor (no prior)
    rna_pmf: np.ndarray  # spliced, EB-shrunk toward global
    gdna_pmf: np.ndarray  # gDNA pool, EB-shrunk toward global
    global_counts: np.ndarray  # raw scanner histogram (all fragments)
    rna_counts: np.ndarray  # raw spliced-annotated histogram
    gdna_counts: np.ndarray  # raw structural pool (intergenic + intronic)
    n_global: float
    n_rna: float
    n_gdna: float
    max_size: int

    def _empirical(self, counts: np.ndarray) -> "FragmentLengthModel":
        """Wrap a raw count vector as an unfinalized ``FragmentLengthModel``.

        Unfinalized on purpose: summary stats (mean/std/median/mode) and
        ``to_dict()`` then read the *raw* histogram, matching how the scanner's
        global + per-category models are reported (no scoring-side pseudo-count
        smoothing enters the QC numbers).
        """
        from ..frag_length_model import FragmentLengthModel

        return FragmentLengthModel(max_size=self.max_size, counts=counts.copy())

    def global_model(self) -> "FragmentLengthModel":
        """Empirical global FL distribution as a ``FragmentLengthModel`` (QC)."""
        return self._empirical(self.global_counts)

    def rna_model(self) -> "FragmentLengthModel":
        """Empirical RNA (spliced-annotated) FL distribution (QC)."""
        return self._empirical(self.rna_counts)

    def gdna_model(self) -> "FragmentLengthModel":
        """Empirical gDNA structural-pool FL distribution (QC).

        The pool is the intergenic + intronic CONTAINED mass (:func:`gdna_fl_mass`): a gDNA-dominated
        proxy that also includes intronic/nascent RNA, not a deconvolved pure-gDNA distribution.
        """
        return self._empirical(self.gdna_counts)


def _pool_sum(payload: "AccumulatorPayload", pools) -> np.ndarray:
    """Sum the named rows of ``payload.pool_lengths`` into one float64 histogram over ``L``.

    ⭐ Binned at ``L``, the molecule length — not at the covered length. Binning at covered length
    collapses the gDNA histogram to a spike at twice the read length, so every long gDNA fragment scores
    as RNA (`CARRY_FORWARD.md` §3 trap 8). The accumulator's ``L`` already includes the mate gap and
    excludes cut introns, so it is the molecule length for both components under one rule.
    """
    pool_lengths = payload.pool_lengths
    if pool_lengths is None:
        return np.zeros(0, dtype=np.float64)
    return np.asarray(pool_lengths, dtype=np.float64)[list(pools)].sum(axis=0)


def gdna_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The pure gDNA length histogram: fragments CONTAINED in an intergenic or intronic node."""
    return _pool_sum(payload, _GDNA_POOLS)


def rna_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The pure RNA length histogram: fragments that used an annotated junction, splice OBSERVED."""
    return _pool_sum(payload, _RNA_POOLS)


def splash_fl_mass(payload: "AccumulatorPayload") -> np.ndarray:
    """The ON-TARGET gDNA length histogram — the two crossing pools, for QC.

    ⚠ Never fold this into :func:`gdna_fl_mass`. Under capture the intergenic pool is *depleted, not
    impure*, so composition stays clean while coverage does not, and on-target gDNA fragments run ~42 bp
    shorter. A model fitted off-target is therefore mis-centred for exactly the fragments that leak —
    and having this as a named pool makes that comparison an output instead of an assumption.
    """
    return _pool_sum(payload, _SPLASH_POOLS)


def _aligned(counts: np.ndarray, max_size: int) -> np.ndarray:
    """Align a raw FL count vector to ``max_size + 1`` bins (overflow → last bin)."""
    out = np.zeros(max_size + 1, dtype=np.float64)
    c = np.asarray(counts, dtype=np.float64)
    n = min(out.size, c.size)
    out[:n] = c[:n]
    if c.size > out.size:
        out[max_size] += float(np.sum(c[out.size :], dtype=np.float64))
    return out


def _normalized(v: np.ndarray) -> np.ndarray:
    s = float(v.sum())
    return v / s if s > 0.0 else np.full(v.size, 1.0 / v.size, dtype=np.float64)


def _smooth_eb(aligned: np.ndarray, global_pmf: np.ndarray, prior_ess: float):
    """Smooth EB pmf: ``(counts + prior_ess·global_pmf) / (total + prior_ess)``.

    ``aligned`` is an already ``max_size``-aligned count vector. Continuous in the
    pool total — no quality threshold (PR04c decision 5). Returns ``(pmf, pool_total)``.
    """
    total = float(aligned.sum())
    denom = total + prior_ess
    pmf = (aligned + prior_ess * global_pmf) / denom if denom > 0.0 else global_pmf.copy()
    return _normalized(pmf), total


def build_fl_models(
    *,
    global_counts: np.ndarray,
    rna_counts: np.ndarray,
    gdna_counts: np.ndarray,
    max_size: int,
    prior_ess: float = POOL_EB_PRIOR_ESS,
) -> FLModels:
    """Build the global / RNA / gDNA FL pmfs under one smooth-EB policy.

    ``global_counts`` is the scanner's unconditional histogram (the EB anchor),
    ``rna_counts`` the spliced (SPLICED-ANNOT) histogram, ``gdna_counts`` the
    aggregated gDNA pool (:func:`gdna_fl_mass`). RNA and gDNA are EB-shrunk toward
    the global FL with a single Dirichlet ``prior_ess``.
    """
    global_aligned = _aligned(global_counts, max_size)
    rna_aligned = _aligned(rna_counts, max_size)
    gdna_aligned = _aligned(gdna_counts, max_size)
    global_pmf = _normalized(global_aligned)
    rna_pmf, n_rna = _smooth_eb(rna_aligned, global_pmf, prior_ess)
    gdna_pmf, n_gdna = _smooth_eb(gdna_aligned, global_pmf, prior_ess)
    return FLModels(
        global_pmf=global_pmf,
        rna_pmf=rna_pmf,
        gdna_pmf=gdna_pmf,
        global_counts=global_aligned,
        rna_counts=rna_aligned,
        gdna_counts=gdna_aligned,
        n_global=float(global_aligned.sum()),
        n_rna=n_rna,
        n_gdna=n_gdna,
        max_size=int(max_size),
    )
