"""rigel.calibration.total_abundance — the MEASURED TOTAL per slot, and the wall mask behind it.

⭐⭐ **The one question this module answers: how many fragments per base sit at this object, with NO
composition model anywhere in the answer?** It is not a gDNA density and not an RNA density — it is
the pooled total, which is the only density a payload can state before anything is deconvolved. What
consumes it is a landscape fitted before pass-0, so its inputs must be counts and lengths and nothing
that was solved.

⛔ **It is NEVER ``mass / effective_length``.** That divisor is a function of the composition being
solved for, so the same 100 counts in 500 bp read 0.25 as pure gDNA and 0.33 as pure RNA. Every term
below is a count over an opportunity that is the SAME for every component:

* **BOUNDARY** — the shipped reciprocal-opportunity banks. The crossing opportunity is ``w − 1`` and
  the deposit ``1/(w − 1)``, so ``E[Σ] = ρ·P(w ≥ 2) = ρ`` for any real library: exact, model-free,
  every fragment length. The sj banks complete the face (a spliced fragment crosses its sj, and
  without that arm a boundary's total is not a total), and the certified SPLICED population enters as
  an INCIDENCE at its own divisor ``mu_r − 1`` — a spliced fragment cannot be gDNA (AXIOM 0), so
  using the RNA pmf there is a certainty, not a composition assumption.
* **REGION** — the START/END banks over the region's own length. A fragment's first covered base
  falls in a region at rate ``ρ·ℓ`` for EVERY fragment length, which is what makes it a TOTAL where
  the contained bank is only a density SHAPE (``ρ·P(w ≤ ℓ)``, 11.6× at a 98 bp exon —
  ``EQUATIONS.md`` §2). The two banks are blind at OPPOSITE template ends, so the consumer
  side-selects: use the side whose wall does not bind, average where both are exact.

⭐ **The wall rule, derived and not tuned.** ``A_start(w | d) = min(ℓ, (d + ℓ − w + 1)₊)``, so the
start form is exact iff the template continues at least ``w_max − 1`` bases past the region's
genomic-HIGH bound; the end form mirrors at the LOW bound. ``w_max`` is READ from the support end of
``deposited_lengths`` — the longest fragment the library actually deposited — never chosen. The
distance is taken at the COMPONENT MINIMUM over the populations the annotation admits at that slot
(``AXIOM 0``'s ``T(slot)``): gDNA's template is the chromosome, a nascent molecule's is its genomic
span (the contiguous boundary reach), a mature molecule's is its SPLICED length
(:class:`~rigel.calibration.splice_graph.MatureWallDistances`). Where BOTH sides bind the slot is
honestly **not model-free** and says so — the total reads NaN rather than a number no consumer can
trust.

⛔ **Nothing here decides anything.** It measures; the landscape fit and the prior are separate rungs
with separate gates. Its own falsification is ``tests/calibration/test_total_abundance.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .effective_length import UNBOUNDED_REACH, crossing_eff_length
from .region_arrays import boundary_region_indices
from .region_chain import BOUNDARY, REGION
from .signature import mrna_active_strands, nrna_active_strands


__all__ = [
    "RegionWallMask",
    "TotalAbundance",
    "build_region_wall_mask",
    "build_total_abundance",
    "region_counts_and_exposure",
    "w_max_from_deposited_lengths",
]


def w_max_from_deposited_lengths(deposited_lengths) -> int:
    """The longest fragment length the library DEPOSITED — the support end of the histogram.

    ⛔ Read, never chosen: the wall rule's bar is ``w_max − 1``, and a quantile would let a
    conservative-looking choice mark a binding wall exact. ``payload.max_length`` is the scanner's
    LIMIT, not the observed maximum, and is the wrong number here.
    """
    hist = np.asarray(deposited_lengths)
    nz = np.nonzero(hist)[0]
    if nz.size == 0:
        raise ValueError(
            "deposited_lengths is empty — no fragment was deposited, so there is no w_max to read "
            "and no wall rule to apply."
        )
    return int(nz.max())


@dataclass(frozen=True, slots=True)
class RegionWallMask:
    """Per-REGION wall verdicts on the accumulator's region axis.

    ``d_low``/``d_high`` are the COMPONENT-MINIMUM template distances past the region's genomic-low
    and genomic-high bounds, in bases. ``end_exact`` is ``d_low ≥ w_max − 1`` and ``start_exact`` is
    ``d_high ≥ w_max − 1``: the START bank is blind at the DOWNSTREAM (genomic-high) end and the END
    bank at the UPSTREAM (genomic-low) one, which is why the pair closes.

    ⚠ Genomic, never transcriptional: a start at ``p`` needs ``w − 1`` bases of template to its
    genomic RIGHT whatever strand the transcript is on, so the minus strand does not flip the sides
    (only which per-strand distance column is read).
    """

    n_regions: int
    d_low: np.ndarray  # float64 (n_regions,)
    d_high: np.ndarray  # float64 (n_regions,)
    start_exact: np.ndarray  # bool (n_regions,)
    end_exact: np.ndarray  # bool (n_regions,)
    w_max: int

    @property
    def double_walled(self) -> np.ndarray:
        """Neither side exact — no deposit rule is model-free at this region."""
        return ~(self.start_exact | self.end_exact)


@dataclass(frozen=True, slots=True)
class TotalAbundance:
    """The per-slot measured total, counts/bp, on the chain axis.

    ``total`` is NaN exactly where the slot is not model-free (a double-walled REGION); ``model_free``
    is the mask a fit must select on, and ``start_used``/``end_used`` record which REGION side (or
    both) the value came from, so a consumer can never mistake an averaged slot for a side-selected
    one.
    """

    n_slots: int
    total: np.ndarray  # float64 (n_slots,)
    model_free: np.ndarray  # bool (n_slots,)
    start_used: np.ndarray  # bool (n_slots,)
    end_used: np.ndarray  # bool (n_slots,)
    w_max: int


def build_region_wall_mask(
    region_arrays,
    mature: "object",
    boundary_reach_lo: np.ndarray,
    boundary_reach_hi: np.ndarray,
    *,
    w_max: int,
) -> RegionWallMask:
    """The wall mask, at the COMPONENT MINIMUM over the populations the annotation admits.

    ``mature`` is a :class:`~rigel.calibration.splice_graph.MatureWallDistances`;
    ``boundary_reach_lo``/``_hi`` are the contiguous-boundary reach arrays ``float64[E, 2]``
    (:func:`~rigel.calibration.splice_graph.build_contiguous_boundary_reach_arrays`) — the NASCENT
    arm, genomic by construction.

    The minimum runs over, per region:

    * **gDNA** — always admitted (it is genomically continuous), template the contig: the distance is
      the region's own offset from the reference ends. This is the only arm at a structurally
      pure-gDNA slot, and it binds only at the outermost regions of each reference.
    * **RNA on strand ``s``** — admitted iff the annotation admits it there (``free_s``). Mature where
      an exon covers the region (the SPLICED distance), nascent where one does not (the genomic reach
      at the region's own wall). An admitted strand with a SHORT template is what makes an interior
      terminal exon inexact, and the mature arm must be able to bind BELOW a long nascent reach —
      that is the measured mature-differential wall.

    ⚠ The reach is keyed per BOUNDARY: a region's low wall is the boundary to its LEFT and its high
    wall the boundary to its RIGHT, and a reference's outermost regions have no such boundary. There
    the RNA arm has no reach to read, so only the gDNA arm bounds it — correct, because a molecule
    of any kind is bounded by the reference there too.
    """
    start = np.asarray(region_arrays.start, dtype=np.int64)
    end = np.asarray(region_arrays.end, dtype=np.int64)
    ref_off = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    sig = np.asarray(region_arrays.signature, dtype=np.int64)
    n = start.shape[0]

    # ── gDNA: distance to the reference ends. The partition tiles a reference, so the first region's
    # start and the last region's end ARE the reference bounds.
    ref_lo = np.repeat(start[ref_off[:-1]], np.diff(ref_off)) if n else start
    ref_hi = np.repeat(end[ref_off[1:] - 1], np.diff(ref_off)) if n else end
    d_low = (start - ref_lo).astype(np.float64)
    d_high = (ref_hi - end).astype(np.float64)

    # ── RNA, per admitted strand. `free_s` is AXIOM 0's T(slot) at a REGION: the strands the
    # annotation admits RNA on at all.
    free_pos, free_neg = nrna_active_strands(sig)
    mrna_pos, mrna_neg = mrna_active_strands(sig)
    lo_reach, hi_reach = _region_reach(region_arrays, boundary_reach_lo, boundary_reach_hi)
    m_low = np.asarray(mature.d_low, dtype=np.float64)
    m_high = np.asarray(mature.d_high, dtype=np.float64)
    m_cov = np.asarray(mature.covered, dtype=bool)

    for col, free, mrna in ((0, free_pos, mrna_pos), (1, free_neg, mrna_neg)):
        # mature where an exon of that strand covers the region; nascent (genomic reach) elsewhere.
        # ⭐ The `free_s ∧ mrna_active_s` licence is REDUNDANT with `covered` — an exon covering a
        # region puts that strand's exon bit in its signature — so it is ASSERTED rather than carried
        # as an ungated term. Where they disagree the annotation and the distances describe different
        # partitions and every verdict on that slot would be arbitrary.
        use_mature = m_cov[:, col]
        if np.any(use_mature & ~(free & mrna)):
            raise ValueError(
                "a region is marked mature-covered but carries no matching exon bit in its "
                "signature — the annotation and the wall distances disagree."
            )
        use_nascent = free & ~use_mature
        cand_low = np.where(use_mature, m_low[:, col], np.inf)
        cand_high = np.where(use_mature, m_high[:, col], np.inf)
        cand_low = np.where(use_nascent, lo_reach[:, col], cand_low)
        cand_high = np.where(use_nascent, hi_reach[:, col], cand_high)
        d_low = np.minimum(d_low, cand_low)
        d_high = np.minimum(d_high, cand_high)

    bar = float(w_max - 1)
    return RegionWallMask(
        n_regions=n,
        d_low=d_low,
        d_high=d_high,
        start_exact=d_high >= bar,
        end_exact=d_low >= bar,
        w_max=int(w_max),
    )


def _region_reach(region_arrays, boundary_reach_lo, boundary_reach_hi):
    """The contiguous-boundary reaches gathered onto the REGION axis — ``(low_wall, high_wall)``,
    each ``float64[R, 2]``, ``inf`` where the region has no boundary on that side.

    ⚠ A reach of 0 is an ANSWER (no template of that strand crosses there); ``inf`` is the absence of
    a boundary object, which is a different thing and must not read as "wall binds".
    """
    lo_e = np.asarray(boundary_reach_lo, dtype=np.float64)
    hi_e = np.asarray(boundary_reach_hi, dtype=np.float64)
    n = np.asarray(region_arrays.start).shape[0]
    low = np.full((n, 2), np.inf)
    high = np.full((n, 2), np.inf)
    lo_region, hi_region = boundary_region_indices(np.asarray(region_arrays.ref_id))
    if lo_region.size:
        # boundary e lies between region e_lo and e_lo + 1: it is the HIGH wall of the lower region
        # and the LOW wall of the upper one. reach_hi is the template above the position, reach_lo
        # below it — so the high wall reads reach_hi and the low wall reads reach_lo.
        high[lo_region] = hi_e
        low[hi_region] = lo_e
    return low, high


def region_counts_and_exposure(
    substrate, region_arrays, wall_mask: RegionWallMask
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The side-selected ``(counts, exposure, model_free)`` triple on the REGION axis.

    ⭐⭐ **This is the entry point a POOLED estimator wants, and it is deliberately not a density.** A
    pooled rate is `Σcounts / Σexposure` — a ratio of sums, never a mean of ratios
    (`TRAPS: a-mean-of-ratios-inherits-the-partition`) — and a conjugate posterior wants the same pair
    (`Gamma(Σcounts + ½, Σexposure)`), so handing a consumer a per-region density would force it to
    re-multiply and lose the pooling. Every consumer that today pools `count / E_contained` can take
    this pair instead and keep its own estimator unchanged.

    ⭐ **Why the pair is better-specified than the one it replaces**: `E[count] = ρ·E_contained` is
    unbiased but the fragment-length pmf enters the DIVISOR, while `E[S] = ρ·ℓ` for every fragment
    length — so the exposure here is a GEOMETRY (the region's own length) rather than a model output.
    Where both walls clear, the two sides are two counts of one rate at one exposure, so the exposure
    doubles with the counts and the pooled ratio is the precision-weighted combination.

    ``model_free`` is False on a double-walled region, where both counts and exposure are returned as
    0 so a caller that pools without masking gets no contribution rather than a wrong one — but a
    caller that CARES (a per-region consumer) must read the mask, and one that pools is safe either
    way.
    """
    s_bank = np.asarray(substrate.region_start_count, dtype=np.float64).sum(axis=1)
    e_bank = np.asarray(substrate.region_end_count, dtype=np.float64).sum(axis=1)
    _check_ledger(s_bank, e_bank)
    ell = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    s_ok = np.asarray(wall_mask.start_exact, dtype=bool)
    e_ok = np.asarray(wall_mask.end_exact, dtype=bool)
    counts = np.where(s_ok, s_bank, 0.0) + np.where(e_ok, e_bank, 0.0)
    exposure = ell * (s_ok.astype(np.float64) + e_ok.astype(np.float64))
    model_free = (s_ok | e_ok) & (ell > 0.0)
    return (
        np.where(model_free, counts, 0.0),
        np.where(model_free, exposure, 0.0),
        model_free,
    )


def _check_ledger(s_bank: np.ndarray, e_bank: np.ndarray) -> None:
    """⛔ The START/END ledger, CHECKED and never assumed: each bank takes one increment per deposited
    fragment, so their totals must agree. A payload where they do not is corrupted input, and
    averaging two different populations there would be a silent wrong answer."""
    if not np.isclose(s_bank.sum(), e_bank.sum(), rtol=0.0, atol=0.5):
        raise ValueError(
            f"the START/END ledger does not close: Σstart = {s_bank.sum():.0f} but "
            f"Σend = {e_bank.sum():.0f}. Each bank takes one increment per deposited fragment, so "
            "this payload is corrupted and its two sides describe different populations."
        )


def build_total_abundance(
    chain,
    substrate,
    region_arrays,
    geometry,
    wall_mask: RegionWallMask,
    rna_fl_pmf: np.ndarray,
) -> TotalAbundance:
    """Assemble the measured total per slot. No solver, no belief, no composition anywhere.

    ⛔ The START/END ledger is CHECKED, not assumed: the two banks each carry one increment per
    deposited fragment, so their totals must agree. A payload where they do not is corrupted input,
    and averaging two different populations there would be a silent wrong answer.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    n_slots = int(chain.n_slots)
    is_region = kind == REGION
    is_boundary = kind == BOUNDARY

    # ⭐ ONE derivation, shared with every pooled consumer: the side selection lives in
    # `region_counts_and_exposure` and this function only turns its pair into a per-slot rate. A second
    # copy of the selection is how a consumer and this total would drift apart.
    counts, exposure, r_free = region_counts_and_exposure(substrate, region_arrays, wall_mask)
    s_ok_r = np.asarray(wall_mask.start_exact, dtype=bool)
    e_ok_r = np.asarray(wall_mask.end_exact, dtype=bool)

    total = np.full(n_slots, np.nan, dtype=np.float64)
    start_used = np.zeros(n_slots, dtype=bool)
    end_used = np.zeros(n_slots, dtype=bool)

    r_idx = obj[is_region]
    with np.errstate(invalid="ignore", divide="ignore"):
        r_rate = np.where(
            r_free[r_idx], counts[r_idx] / np.where(exposure[r_idx] > 0.0, exposure[r_idx], 1.0), np.nan
        )
    total[is_region] = r_rate
    start_used[is_region] = s_ok_r[r_idx] & r_free[r_idx]
    end_used[is_region] = e_ok_r[r_idx] & r_free[r_idx]

    # ── BOUNDARY: the exact crossing banks + the sj faces + certified spliced at mu_r − 1.
    inv = np.asarray(geometry.inv_abundance, dtype=np.float64)
    sj_lo = np.asarray(geometry.inv_sj_lo, dtype=np.float64).sum(axis=1)
    sj_hi = np.asarray(geometry.inv_sj_hi, dtype=np.float64).sum(axis=1)
    spliced = np.asarray(geometry.spliced_count, dtype=np.float64).sum(axis=1)
    eff_spliced = float(crossing_eff_length(rna_fl_pmf, UNBOUNDED_REACH, UNBOUNDED_REACH))
    if not eff_spliced > 0.0:
        raise ValueError(
            "the certified-spliced divisor mu_r - 1 is not positive — the RNA fragment-length pmf "
            "cannot price an incidence."
        )
    total[is_boundary] = (
        inv[is_boundary]
        + sj_lo[is_boundary]
        + sj_hi[is_boundary]
        + spliced[is_boundary] / eff_spliced
    )

    return TotalAbundance(
        n_slots=n_slots,
        total=total,
        model_free=np.isfinite(total),
        start_used=start_used,
        end_used=end_used,
        w_max=int(wall_mask.w_max),
    )
