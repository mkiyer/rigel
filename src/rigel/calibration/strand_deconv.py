"""Per-node gDNA/RNA deconvolution by a strand-vs-count handoff (decoupled architecture).

Each **node** (every region and every boundary side, kept separate) is deconvolved into gDNA / RNA
mass by **routing** it to one of two independent estimators — never a product of the two:

    route_strand = strand_identifiable (global)  AND  strand_observable (per-node)
      ├─ TRUE  → STRAND module:  Beta-Binomial posterior over the gDNA fraction, weak Beta(½,½) prior
      └─ FALSE → COUNT  module:  gDNA fraction = count_gdna_frac (density·eff_len / mass, raw)

Why a handoff, not a joint product: the strand estimator is **unbiased** (gDNA is symmetric at ½,
RNA is at κ) but noisy at low N/SS; the count estimator is **biased** under hybrid capture (exon
density imputed from depleted off-target neighbours). Mixing a biased estimate into an unbiased one
re-introduces bias, so we use strand wherever it has signal (strand-observable nodes in a
strand-identifiable library) and count only where strand structurally has none (AMBIG /
no-defined-sense nodes, and every node when the library is unstranded). The two act on disjoint node
sets and cannot conflict. See ``docs/calibration/decoupled_calibration_design.md``; the retired joint
product is archived in ``docs/calibration/archive/joint_deconvolution.md``.

The **strand module** is the Beta-Binomial posterior (the exact, clip-free robust strand
deconvolution): ``log_post(g) = log Beta(g; ½, ½) + strand_loglik(sense, antisense | g, κ, od)`` on a
grid, reported as the posterior median (then optionally log-odds-shifted by ``gdna_strand_llr_bias``,
the FP-aversion knob; default 0). As κ→½ the likelihood flattens → posterior → prior (uninformative),
but the global identifiability gate normally routes unstranded libraries to count, so the no-op is
structural. Amounts use the conserved **fractional mass** ``M``: gDNA = ``g·M``,
RNA = ``(1−g)·M`` + the deterministic spliced mass. Regions and boundaries are combined into loci
only **after** calibration (``assemble_priors``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signature import TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from .strand_likelihood import strand_loglik

_GRID_EPS = 1.0e-4
_STRAND_PRIOR = 0.5  # Jeffreys Beta(½,½): the standard uninformative prior for the strand posterior


@dataclass(frozen=True, slots=True)
class NodeDeconv:
    """Per-node deconvolution result (regions or boundary sides; kept separate)."""

    gdna_mass: np.ndarray  # float64[K]
    rna_mass: np.ndarray  # float64[K]  (= (1−gdna_frac)·M_unspliced + spliced mass)
    gdna_frac: np.ndarray  # float64[K] — reported gDNA fraction of the UNSPLICED mass
    gdna_frac_var: np.ndarray  # float64[K] — posterior variance (0 for count-routed nodes)


def _deconv_per_node(
    mass_unspl,
    mass_spliced,
    sense,
    antisense,
    strand_observable,
    count_gdna_frac,
    *,
    strand_identifiable,
    rna_sense_frac,
    gdna_strand_overdispersion,
    rna_strand_overdispersion,
    gdna_strand_llr_bias,
    n_grid,
) -> NodeDeconv:
    """Route each node to the strand posterior or the count fraction, then split its mass.

    ``strand_identifiable`` is the global gate (is κ usable?); ``strand_observable[i]`` the per-node
    gate (is the node's sense defined?). A node routes to STRAND iff both hold and it carries flux;
    otherwise it takes the COUNT module's ``count_gdna_frac[i]`` (a raw density ratio in ``[0,1]``).
    """
    k = mass_unspl.shape[0]
    grid = np.linspace(_GRID_EPS, 1.0 - _GRID_EPS, n_grid)
    log_grid = np.log(grid)
    log_1mgrid = np.log1p(-grid)
    # Weak Jeffreys prior over the gDNA fraction (strand module); the strand likelihood dominates it
    # wherever we route to strand (κ identifiable ⇒ non-flat likelihood).
    log_prior = (_STRAND_PRIOR - 1.0) * log_grid + (_STRAND_PRIOR - 1.0) * log_1mgrid
    gdna = np.zeros(k)
    rna = np.zeros(k)
    gdna_frac = np.zeros(k)
    gdna_frac_var = np.zeros(k)
    for i in range(k):
        m = float(mass_unspl[i])
        if m <= 0.0:
            rna[i] = float(mass_spliced[i])  # only deterministic spliced RNA, if any
            continue
        if strand_identifiable and strand_observable[i] and (sense[i] + antisense[i]) > 0:
            # STRAND module: Beta-Binomial posterior, weak prior, no count term.
            log_post = log_prior + strand_loglik(
                grid,
                sense[i],
                antisense[i],
                rna_sense_frac,
                gdna_strand_overdispersion=gdna_strand_overdispersion,
                rna_strand_overdispersion=rna_strand_overdispersion,
            )
            w = np.exp(log_post - log_post.max())
            w /= w.sum()
            mean = float(np.dot(grid, w))
            gdna_frac_var[i] = float(np.dot((grid - mean) ** 2, w))
            frac = float(np.clip(np.interp(0.5, np.cumsum(w), grid), 0.0, 1.0))  # posterior median
            # gDNA FP-aversion LLR bias: shift the reported fraction's log-odds by λ nats (the
            # calibration twin of the EM gdna_em_llr_bias; λ=0 ⇒ neutral, the default).
            if gdna_strand_llr_bias != 0.0 and 0.0 < frac < 1.0:
                logit_frac = np.log(frac) - np.log1p(-frac)
                frac = float(1.0 / (1.0 + np.exp(-(gdna_strand_llr_bias + logit_frac))))
        else:
            # COUNT module: the count clue's gDNA fraction directly (point estimate, no posterior).
            frac = min(max(float(count_gdna_frac[i]), 0.0), 1.0)
        gdna_frac[i] = frac
        gdna[i] = frac * m
        rna[i] = (1.0 - frac) * m + float(mass_spliced[i])
    return NodeDeconv(
        gdna_mass=gdna, rna_mass=rna, gdna_frac=gdna_frac, gdna_frac_var=gdna_frac_var
    )


def deconv_regions(
    substrate,
    region_arrays,
    node_density,
    *,
    strand_identifiable,
    rna_sense_frac,
    gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0,
    gdna_strand_llr_bias=0.0,
    n_grid=200,
) -> NodeDeconv:
    """Deconvolve each region's contained mass (a node) into gDNA / RNA by the strand/count handoff.

    A region is strand-observable iff its transcript strand is defined (``TS_POS`` / ``TS_NEG``);
    ``TS_NONE`` (intergenic) and ``TS_AMBIG`` route to count. The count fraction is the precomputed
    ``node_density.count_gdna_frac``.
    """
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    pos = c.n_unspliced_pos.astype(np.float64)
    neg = c.n_unspliced_neg.astype(np.float64)
    sense = np.where(ts == TS_NEG, neg, pos)  # orient to transcript sense
    antisense = (pos + neg) - sense
    strand_observable = (ts == TS_POS) | (ts == TS_NEG)
    return _deconv_per_node(
        c.mass_unspliced,
        c.mass_spliced,
        sense,
        antisense,
        strand_observable,
        node_density.count_gdna_frac,
        strand_identifiable=strand_identifiable,
        rna_sense_frac=rna_sense_frac,
        gdna_strand_overdispersion=gdna_strand_overdispersion,
        rna_strand_overdispersion=rna_strand_overdispersion,
        gdna_strand_llr_bias=gdna_strand_llr_bias,
        n_grid=n_grid,
    )


@dataclass(frozen=True, slots=True)
class _SideQuantities:
    """Per-region boundary-side quantities shared by the deconvolution and the strand-fit seeds."""

    sense: np.ndarray
    antisense: np.ndarray
    n_side: np.ndarray  # pos + neg (count evidence)
    mass: np.ndarray
    mass_spliced: np.ndarray
    count_gdna_frac: np.ndarray  # count module's gDNA fraction for this side (raw density ratio)
    strand_observable: np.ndarray
    count_observable: np.ndarray


def _left_right_neighbors(ts, ref_id, boundary_count_observable):
    """Per-region neighbour strand class + same-ref + boundary count-observability for both sides.

    Region ``r``'s LEFT side is the right side of boundary ``(r−1, r)`` (neighbour ``r−1``); its
    RIGHT side is the left side of boundary ``(r, r+1)`` (neighbour ``r+1``). ``boundary_count_observable[r]``
    describes boundary ``(r, r+1)``. Returns
    ``(left_same, ts_prev, left_observable, right_same, ts_next, right_observable)``.
    """
    r = ts.shape[0]
    left_same = np.zeros(r, dtype=bool)
    ts_prev = np.zeros(r, dtype=ts.dtype)
    left_observable = np.zeros(r, dtype=bool)
    right_same = np.zeros(r, dtype=bool)
    ts_next = np.zeros(r, dtype=ts.dtype)
    if r > 1:
        left_same[1:] = ref_id[1:] == ref_id[:-1]
        ts_prev[1:] = ts[:-1]
        left_observable[1:] = boundary_count_observable[:-1]
        right_same[:-1] = ref_id[:-1] == ref_id[1:]
        ts_next[:-1] = ts[1:]
    return left_same, ts_prev, left_observable, right_same, ts_next, boundary_count_observable


def _compute_side(
    view, same, ts_self, ts_other, side_count_observable, eff, region_density
) -> _SideQuantities:
    """Per-side sense split, count-prior density, and strand-observability.

    A side is **strand-observable** iff its boundary's two regions define a single, consistent
    transcript sense (so 'sense' is defined). An **intergenic (TS_NONE) region is a strand wildcard**
    — it carries no transcript, so it is compatible with either strand: a gene-edge boundary
    (stranded exon ↔ intergenic) is therefore oriented by the gene side, ``{POS,NONE}→POS``,
    ``{NEG,NONE}→NEG`` (same as ``{POS,POS}`` / ``{NEG,NEG}``). Only a genuine conflict — opposite
    strands ``{POS,NEG}`` or a ``TS_AMBIG`` flank (overlapping +/− transcripts) — leaves the sense
    undefined (CALIBRATION_TODO Issue #3; treating gene-edge sides as strand-blind leaked ~half their
    pure-gDNA crossing mass). The **count** fraction is the raw crossing density ratio: a
    count-observable side (no shared exon) reads ``count_gdna_frac → 1`` from its own crossing mass;
    otherwise it borrows the swept region density. (No strand cleaning — the count module is raw.)
    """
    mass = view.mass_unspliced
    pos = view.n_unspliced_pos.astype(np.float64)
    neg = view.n_unspliced_neg.astype(np.float64)
    # Strand sense with TS_NONE (intergenic) as a wildcard: POS-sense if each flank is POS-or-NONE
    # and at least one is POS (mirror for NEG); a TS_AMBIG flank or an opposite-strand pairing leaves
    # it undefined. Reduces to {POS,POS}/{NEG,NEG} when neither flank is intergenic.
    either_ambig = (ts_self == TS_AMBIG) | (ts_other == TS_AMBIG)
    self_pos_or_none = (ts_self == TS_POS) | (ts_self == TS_NONE)
    other_pos_or_none = (ts_other == TS_POS) | (ts_other == TS_NONE)
    self_neg_or_none = (ts_self == TS_NEG) | (ts_self == TS_NONE)
    other_neg_or_none = (ts_other == TS_NEG) | (ts_other == TS_NONE)
    cons_pos = (
        same & ~either_ambig & self_pos_or_none & other_pos_or_none
        & ((ts_self == TS_POS) | (ts_other == TS_POS))
    )
    cons_neg = (
        same & ~either_ambig & self_neg_or_none & other_neg_or_none
        & ((ts_self == TS_NEG) | (ts_other == TS_NEG))
    )
    strand_observable = cons_pos | cons_neg
    sense = np.where(cons_neg, neg, pos)
    n_side = pos + neg
    antisense = n_side - sense
    with np.errstate(divide="ignore", invalid="ignore"):
        # Raw crossing density (no strand cleaning); count-observable sides use their own crossing
        # density, others borrow the swept region density.
        own = np.where((mass > 0.0) & (eff > 0.0), mass / np.maximum(eff, 1e-12), 0.0)
        density = np.where(side_count_observable, own, region_density)
        count_gdna_frac = np.clip(
            np.where(mass > 0.0, density * eff / np.maximum(mass, 1e-12), 0.0), 0.0, 1.0
        )
    return _SideQuantities(
        sense=sense,
        antisense=antisense,
        n_side=n_side,
        mass=mass,
        mass_spliced=view.mass_spliced,
        count_gdna_frac=count_gdna_frac,
        strand_observable=strand_observable,
        count_observable=np.asarray(side_count_observable, dtype=bool),
    )


def deconv_sides(
    substrate,
    region_arrays,
    node_density,
    boundary_side_eff_len,
    *,
    strand_identifiable,
    rna_sense_frac,
    gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0,
    gdna_strand_llr_bias=0.0,
    n_grid=200,
) -> tuple[NodeDeconv, NodeDeconv]:
    """Deconvolve each boundary **side** as an independent node by the strand/count handoff.

    Region ``r`` owns the **right** side of its left boundary (``substrate.left[r]``) and the
    **left** side of its right boundary (``substrate.right[r]``); both lie inside region ``r`` and so
    use ``boundary_side_eff_len[r]``. Returns ``(left, right)`` per-region ``NodeDeconv`` results
    (zero where a side carries no mass — e.g. reference edges).
    """
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    eff = np.asarray(boundary_side_eff_len, dtype=np.float64)
    region_density = node_density.density
    l_same, ts_prev, l_obs, r_same, ts_next, r_obs = _left_right_neighbors(
        ts, ref_id, node_density.boundary_count_observable
    )

    def _deconv(view, same, ts_other, side_obs):
        sq = _compute_side(view, same, ts, ts_other, side_obs, eff, region_density)
        return _deconv_per_node(
            sq.mass,
            sq.mass_spliced,
            sq.sense,
            sq.antisense,
            sq.strand_observable,
            sq.count_gdna_frac,
            strand_identifiable=strand_identifiable,
            rna_sense_frac=rna_sense_frac,
            gdna_strand_overdispersion=gdna_strand_overdispersion,
            rna_strand_overdispersion=rna_strand_overdispersion,
            gdna_strand_llr_bias=gdna_strand_llr_bias,
            n_grid=n_grid,
        )

    left = _deconv(substrate.left, l_same, ts_prev, l_obs)
    right = _deconv(substrate.right, r_same, ts_next, r_obs)
    return left, right


def boundary_side_seeds(substrate, region_arrays, node_density, boundary_side_eff_len):
    """``(sense, total, gdna_weight)`` seed arrays from count- & strand-observable boundary sides.

    The exon–intron / exon–intergenic seam seeds for the gDNA strand-overdispersion fit
    (:mod:`gdna_strand`), complementing the contained-region seeds (needed under hybrid capture,
    which depletes off-target intergenic/intronic gDNA). The weight is the count-clue gDNA fraction
    ``count_gdna_frac`` (≈ 1 for a count-observable side, gDNA by signature) — the same quantity the
    contained-region seeds use (``node_density.count_gdna_frac``).
    """
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    eff = np.asarray(boundary_side_eff_len, dtype=np.float64)
    region_density = node_density.density
    l_same, ts_prev, l_obs, r_same, ts_next, r_obs = _left_right_neighbors(
        ts, ref_id, node_density.boundary_count_observable
    )
    sides = (
        _compute_side(substrate.left, l_same, ts, ts_prev, l_obs, eff, region_density),
        _compute_side(substrate.right, r_same, ts, ts_next, r_obs, eff, region_density),
    )
    senses, totals, weights = [], [], []
    for sq in sides:
        seed = sq.count_observable & sq.strand_observable & (sq.n_side > 0.0)
        senses.append(sq.sense[seed])
        totals.append(sq.n_side[seed])
        weights.append(sq.count_gdna_frac[seed])
    return np.concatenate(senses), np.concatenate(totals), np.concatenate(weights)


__all__ = ["NodeDeconv", "deconv_regions", "deconv_sides", "boundary_side_seeds"]
