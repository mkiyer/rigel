#!/usr/bin/env python
"""⛔⛔ **THE SOFT-MIN PER-TRANSCRIPT RNA PRIOR WEIGHT — BUILT, PRICED AND REFUSED (2026-08-13).**

⛔ **This module is a RECORD, not a proposal. Read `docs/ISSUES.md` (`refused-soft-min-path-weighting`) before reusing any of it.**
16 arms on `g00 ss0.99 capture_off` and 3 on the blind stratum `g50 ss0.50 capture_on`, base re-recorded
in the same session: **worse than base at transcript level everywhere (1.26–1.60×)**, and the one good
number — gene error **0.395–0.527×** at `g00` — collapses to **1.006–1.041×** on the blind stratum.
The mechanism is `TRAPS: an-upper-bound-is-not-an-estimate`: the bound is real, but **3,644 of 4,839
silent transcripts share an object with an expressed one** and inherit its bound.

⭐⭐ **KEPT BECAUSE FOUR THINGS IN IT ARE SETTLED AND THE NEXT CANDIDATE SHOULD START FROM THEM** — the
density identity below, the measured multiplier ruling, the opportunity weighting, and the monotone dial
which says the owner's theorem is NOT what failed. ⛔ What must change is the SUPPORT, not the mean:
0.0 % of expressed transcripts are ever zeroed and the entire error is the silent ones that are not.

    python scripts/design/transcript_weights.py        # its own falsification; no data needed

`assemble_priors` produces ONE unspliced-RNA pseudocount budget per locus. Stage 5 proved that splitting
it with the TRUE relative abundances removes 84 % of transcript misassignment, and that a flipped split
is 26 % worse — so the machinery converts good weights into good answers, and the open question is only
where good weights come from.

⭐⭐ **THE OWNER'S THEOREM.** A transcript's density is bounded by the MINIMUM density along its path:
mass at any object is SHARED with every other transcript covering it, so no transcript can be denser
than the thinnest place it must pass through. A naive `min` is destroyed by sampling noise at short
objects, so the shipped form is a SOFT min.

⛔⛔ **AND THE SOFT MIN IS NOT AN UNWEIGHTED MEAN, WHICH IS WHERE THIS DEVIATES FROM THE SPEC AS
WRITTEN.** `TRAPS: a-mean-of-ratios-inherits-the-partition` was measured *while designing this prior*:
a region boundary appears wherever a signature changes — an antisense feature on the other strand splits a
uniform stretch in two — so **how many objects a transcript's path contains is an artefact of the
annotation, not biology.** Subdividing one 20 kb intron from 1 region to 4 to 10, mass and opportunity held
constant, moved a shadow span's share of its locus's prior from **9.4 % → 18.6 % → 32.4 %**. An
unweighted harmonic, geometric or arithmetic mean over path objects inherits that artefact exactly.

⭐⭐ **THE REPAIR KEEPS THE THEOREM: weight each object by its OWN OPPORTUNITY.** Splitting an object
into pieces **of the same density** splits its opportunity with it, so a weighted power mean is
unchanged for every ``p`` while an unweighted one is not — and that split is precisely the artefact,
because a signature change can region_bound a uniform stretch into 990 bp and 10 bp and an unweighted mean then
gives the sliver equal say. :func:`repartition_invariance` measures it.

⛔ **BE PRECISE ABOUT WHAT IS AND IS NOT INVARIANT, BECAUSE THE OVER-CLAIM IS TEMPTING.** Only ``p = 1``
is invariant to an *arbitrary* re-partition; at ``p < 1`` a finer partition that reveals REAL variation
does move the answer — necessarily, since responding to variation along the path is the entire content of
the theorem. What the opportunity weighting removes is the artefactual part, and that part is large:
measured below, a 10 bp sliver at low density reads **0.5455** on an unweighted harmonic mean against
**2.7523** weighted — a **5.05×** gap — and subdividing that sliver ten ways drags it on to 0.3267
(**8.4×**) while the weighted answer does not move at all.

The family is one dial::

    M_p = ( Σ_o (O_o / ΣO) · d_o^p )^(1/p)

    p =  1   arithmetic   ≡ Σmass / Σopportunity, the pooled form the trap ENDORSES — and the form
                            with no deconvolution at all: a silent transcript inherits its loud
                            neighbour's mass
    p =  0   geometric
    p = -1   harmonic     the owner's preference — the closest well-behaved approach to the min
    p = -inf min          the theorem verbatim, and maximally noise-sensitive

So `p` is not a tuned constant: it is the axis from "no deconvolution" to "the theorem", and the four
rungs are the owner's own candidate list. ⛔ Pick it by measurement, not by taste.

⭐ **ZEROS ARE OMITTED, NOT RESCUED.** A zero over a large opportunity is a strong statement, but under
hybrid capture it may equally mean *no probe here* — the same ambiguity the effective-length shrinkage
already refuses to resolve. So a zero object leaves the mean. ⛔⛔ **A transcript whose path is
ENTIRELY zero gets weight 0, and that is the load-bearing case, not an boundary case:** 4,579 of 8,750
annotated transcripts are silent at `g00`, and the refused first attempt's per-object `+½` revived
exactly that half and took false-positive mass **18.6 M → 41.6 M** (`ISSUES: refused-transcript-weights`). ⛔ There is no
floor here and there must not be one — `g00`'s lesson transposed is that doubt resolves to ABSENT.

⭐ **THE WEIGHT IS A COUNT, NOT A DENSITY**, because the budget being split is a fragment count::

    w_t = softmin_density(t) × unspliced_opportunity(t)

⭐ And the second factor is the owner's issue (2) — *transcript GEOMETRY governs how many unspliced
fragments it can support*. An unspliced fragment cannot cross an intron, so a transcript's unspliced
opportunity is `Σ_exons contained_eff_length(exon_len)` on the RNA pmf: many tiny exons ⇒ nearly all its
fragments are spliced ⇒ it can claim little of an UNSPLICED budget however abundant it is. ⛔ It is NOT
`Σ_regions rna_region_eff_len`, which would forbid crossing an interior region boundary that a real
fragment crosses freely.

⛔ **SPLICE-SJ STEPS CARRY NO WEIGHT, BY THE OWNER'S OWN KEY REALIZATION.** The prior is UNSPLICED
fragments only — a spliced fragment has no gDNA candidate in the EM and is assigned directly — so a
sj, which is spliced by construction, has zero unspliced mass and must not enter the mean. It stays
in the path because the path has other consumers.

⚠ **What this module is NOT.** It does not score end to end and must never learn to: `quant_accuracy.py`
owns that, and two scorers is how a baseline and a ceiling drift apart. Standalone it reports how well a
candidate weight vector agrees with the per-transcript truth, which is a property of the VECTOR.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]

from rigel.calibration.effective_length import contained_eff_length  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    STEP_BOUNDARY,
    STEP_REGION,
    build_transcript_path,
)
from rigel.types import IntervalType  # noqa: E402

#: The four rungs of the power-mean dial, named. ⭐ `arithmetic` is the trap-endorsed pooled form and is
#: therefore the CONTROL: it does the least deconvolution of the four, so a soft-min that cannot beat it
#: is not earning its noise.
MODES = ("arithmetic", "geometric", "harmonic", "min")

#: What the soft-min density is multiplied by to become a count. ⭐⭐ **The density's units settle this,
#: and the derivation is short enough to state here.** Transcript ``t`` with ``A_t`` fragments over an
#: effective length ``E_t`` puts ``A_t·ceff(r)/E_t`` contained fragments in region ``r``, so::
#:
#:     density(r) = mass_rna_region[r] / ceff(r) = Σ_{t ∋ r} A_t / E_t
#:
#: ⭐ **The region's own opportunity CANCELS EXACTLY** — which is what makes densities at objects of
#: wildly different sizes commensurate, and what the theorem's ``min`` is a min OF. So the soft min
#: estimates ``A_t / E_t`` and the multiplier chooses which count comes back out:
#:
#: ``total`` — ``E_t = contained_eff_length(Σ exon_len)``, the transcript's own effective length ⇒ the
#: TOTAL abundance ``A_t``. ⭐ The measured target: ``oracle_alloc`` (total truth) beats
#: ``oracle_alloc_unspliced`` **0.163× against 0.202×** at transcript level.
#: ``full`` — ``Σ_exons contained_eff_length(exon_len)`` ⇒ the UNSPLICED count. An unspliced fragment
#: cannot cross an intron, so each exon contributes separately; this is the owner's issue (2), the
#: geometry that makes a many-tiny-exon transcript unable to claim an unspliced budget. ⚠ It targets
#: the quantity the pseudocounts LITERALLY are, and the oracle arms say that is the worse share —
#: an unresolved tension, priced rather than argued.
#: ``seen`` — only the opportunity of the objects that were NONZERO. The sparsity reading: a transcript
#: with 1 of 40 objects populated is down-weighted 40×, with no new parameter.
#: ⛔ ``full`` and ``seen`` disagree most under hybrid capture, which is why both exist.
OPPORTUNITIES = ("total", "full", "seen")

#: At what granularity the soft min is taken. ⭐⭐ ``gene`` exists because the FIRST measurement of this
#: family said so, on all 12 transcript-granularity arms at once: **worse than `base` at transcript level
#: (1.32–1.60×) and much better at gene level (0.395–0.527×)**, without exception. The theorem gives an
#: UPPER BOUND, and isoforms of one gene share nearly all their objects, so they all inherit the same
#: bound and the weight then separates them by GEOMETRY — which is not abundance. A confidently wrong
#: within-gene split is worse than none.
#: ``gene`` therefore takes the soft min over the gene's objects (deduplicated) and splits it across the
#: gene's isoforms by effective length alone — the explicit *no information within the gene* statement,
#: which is what "equal molar abundance" means. ⛔ It is a RETREAT and is labelled as one:
#: `ISSUES: refused-soft-min-path-weighting` records that the gene axis is too coarse for **62.3 %** of mass, whose within-gene
#: split IS recoverable. It buys the gene-level win only if the isoform-level harm is real.
GRANULARITIES = ("transcript", "gene")


def region_unspliced_density(calibration) -> np.ndarray:
    """float64[n_regions] — deconvolved RNA fragments per admissible start position.

    ⭐ ``mass_rna_region`` is already a FRAGMENT count: containment is exclusive, so a contained fragment
    deposits on exactly one region and needs no incidence→fragment conversion. It is also unspliced by
    construction — ``region_contained`` is credited only when the fragment used no sj.
    """
    mass = np.asarray(calibration.mass_rna_region, dtype=np.float64)
    opp = np.asarray(calibration.rna_region_eff_len, dtype=np.float64)
    out = np.zeros_like(mass)
    np.divide(mass, opp, out=out, where=opp > 0.0)
    return out


def boundary_unspliced_density(calibration) -> np.ndarray:
    """float64[n_boundaries] — the same quantity on the crossing axis.

    ⛔ **Three corrections, and dropping any one changes the units.** ``mass_rna_boundary`` is
    spliced-INCLUSIVE, so the spliced part is subtracted first; the remainder is an object-INCIDENCE
    count, so it is converted to fragments by the unspliced pool's own ``boundary_mass_per_crossing``; only
    then is it a numerator commensurate with the region axis's.
    """
    unspliced = np.maximum(
        np.asarray(calibration.mass_rna_boundary, dtype=np.float64)
        - np.asarray(calibration.mass_rna_spliced_boundary, dtype=np.float64),
        0.0,
    )
    frags = unspliced * np.asarray(calibration.boundary_mass_per_crossing, dtype=np.float64)
    opp = np.asarray(calibration.rna_boundary_eff_len, dtype=np.float64)
    out = np.zeros_like(frags)
    np.divide(frags, opp, out=out, where=opp > 0.0)
    return out


def transcript_opportunities(index, rna_fl_pmf: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(unspliced, total)``, each float64[n_transcripts], on the RNA pmf.

    ⭐ **They differ by exactly one thing: whether the exon lengths are summed BEFORE or AFTER
    ``contained_eff_length``.** Summing after forbids crossing an intron, which is what an unspliced
    fragment cannot do; summing before is the ordinary transcript effective length, over the spliced
    molecule a fragment actually samples.

    ⛔ Neither is ``Σ_regions rna_region_eff_len``, and that would be the easy mistake — it would also
    forbid crossing an INTERIOR region boundary, which a real fragment crosses freely, so it understates
    every transcript whose exons the partition happened to region_bound.
    """
    iv = pd.read_feather(os.path.join(index.index_dir, "intervals.feather"))
    ex = iv[(iv["interval_type"] == int(IntervalType.EXON)) & (iv["t_index"] >= 0)]
    t_idx = ex["t_index"].to_numpy(np.int64)
    lens = (ex["end"].to_numpy(np.int64) - ex["start"].to_numpy(np.int64)).astype(np.float64)
    n_t = int(index.num_transcripts)
    unspliced = np.bincount(t_idx, weights=contained_eff_length(lens, rna_fl_pmf), minlength=n_t)
    spliced_len = np.bincount(t_idx, weights=lens, minlength=n_t)
    return unspliced.astype(np.float64), contained_eff_length(spliced_len, rna_fl_pmf)


def _power_mean(dens: np.ndarray, opp: np.ndarray, seg: np.ndarray, n_t: int, mode: str):
    """The opportunity-weighted power mean of ``dens`` within each segment, plus the seen opportunity.

    ``seg[s]`` is the transcript owning step ``s``. Only steps with ``dens > 0`` and ``opp > 0`` are
    passed in, so "omit the zeros" has already happened and a transcript with no surviving step falls
    out as exactly 0 rather than as a fallback.
    """
    tot = np.bincount(seg, weights=opp, minlength=n_t)
    seen = tot.copy()
    live = tot > 0.0
    out = np.zeros(n_t, dtype=np.float64)
    if mode == "arithmetic":
        num = np.bincount(seg, weights=opp * dens, minlength=n_t)
        np.divide(num, tot, out=out, where=live)
    elif mode == "geometric":
        num = np.bincount(seg, weights=opp * np.log(dens), minlength=n_t)
        np.divide(num, tot, out=out, where=live)
        out = np.where(live, np.exp(out, where=live, out=np.zeros_like(out)), 0.0)
    elif mode == "harmonic":
        num = np.bincount(seg, weights=opp / dens, minlength=n_t)
        np.divide(tot, num, out=out, where=num > 0.0)
    elif mode == "min":
        out = np.full(n_t, np.inf, dtype=np.float64)
        np.minimum.at(out, seg, dens)
        out = np.where(np.isfinite(out), out, 0.0)
    else:  # pragma: no cover - guarded by argparse choices
        raise ValueError(f"unknown mode {mode!r}")
    return out, seen


def build_weights(
    calibration,
    region_arrays,
    index,
    rna_fl_pmf: np.ndarray,
    *,
    mode: str = "harmonic",
    opportunity: str = "full",
    granularity: str = "transcript",
    path=None,
) -> np.ndarray:
    """float64[n_transcripts] — the per-transcript RNA prior weight. **The whole deliverable.**

    ⛔ Only RATIOS within a locus are read by the EM, so this is deliberately un-normalised.
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    if opportunity not in OPPORTUNITIES:
        raise ValueError(f"opportunity must be one of {OPPORTUNITIES}, got {opportunity!r}")
    if granularity not in GRANULARITIES:
        raise ValueError(f"granularity must be one of {GRANULARITIES}, got {granularity!r}")

    path = build_transcript_path(index, region_arrays) if path is None else path
    n_t = path.n_transcripts

    d_region = region_unspliced_density(calibration)
    d_boundary = boundary_unspliced_density(calibration)
    o_region = np.asarray(calibration.rna_region_eff_len, dtype=np.float64)
    o_boundary = np.asarray(calibration.rna_boundary_eff_len, dtype=np.float64)

    kind = np.asarray(path.kind)
    obj = np.asarray(path.obj_id, dtype=np.int64)
    counts = np.diff(np.asarray(path.offsets, dtype=np.int64))
    seg_all = np.repeat(np.arange(n_t, dtype=np.int64), counts)

    dens = np.zeros(kind.shape[0], dtype=np.float64)
    opp = np.zeros(kind.shape[0], dtype=np.float64)
    is_r = kind == STEP_REGION
    is_b = kind == STEP_BOUNDARY
    dens[is_r], opp[is_r] = d_region[obj[is_r]], o_region[obj[is_r]]
    dens[is_b], opp[is_b] = d_boundary[obj[is_b]], o_boundary[obj[is_b]]
    # ⛔ STEP_SPLICE_SJ keeps dens = opp = 0 and is dropped by the mask below — the prior is
    # UNSPLICED fragments, and a sj has none by construction.

    keep = (dens > 0.0) & (opp > 0.0)
    seg = seg_all
    if granularity == "gene":
        # ⭐ Re-segment onto the GENE axis, and DEDUPLICATE: an object shared by three isoforms is one
        # measurement, not three. Counting it per-isoform would re-import a partition artefact of the
        # same family the opportunity weighting exists to remove — this time on the annotation's isoform
        # count rather than on its boundary count.
        gene = index.t_df.sort_values("t_index")["g_index"].to_numpy(np.int64)
        seg = gene[seg_all]
        key = np.stack([seg[keep], kind[keep].astype(np.int64), obj[keep]])
        _, first = np.unique(key, axis=1, return_index=True)
        uniq = np.zeros(int(keep.sum()), dtype=bool)
        uniq[first] = True
        k2 = keep.copy()
        k2[keep] = uniq
        keep = k2
        n_seg = int(gene.max()) + 1
        soft_min_g, seen_g = _power_mean(dens[keep], opp[keep], seg[keep], n_seg, mode)
        soft_min, seen_opp = soft_min_g[gene], seen_g[gene]
    else:
        soft_min, seen_opp = _power_mean(dens[keep], opp[keep], seg[keep], n_t, mode)

    if opportunity == "seen":
        scale = seen_opp
    else:
        unspliced, total = transcript_opportunities(index, rna_fl_pmf)
        scale = total if opportunity == "total" else unspliced
    return soft_min * scale


def repartition_invariance(mode: str, *, weighted: bool = True) -> tuple[float, float, float]:
    """⭐⭐ **THE PROPERTY `TRAPS: a-mean-of-ratios-inherits-the-partition` SAYS TO TEST DIRECTLY.**

    ⛔ **The construction has to be one an estimator can FAIL, and the obvious one is not.** Splitting a
    uniform stretch into EQUAL pieces of equal density leaves every mean — weighted or not — returning
    the same number, so that version tests nothing (`TRAPS: could-the-arm-have-fired`; a fixture is an
    arm). What separates them is an UNEVEN region_bound, which is what the annotation actually does: a signature
    change from an antisense feature slices a 1,000 bp stretch into 990 bp and a 10 bp sliver.

    So: 990 bp at density 3.0 beside a 10 bp sliver at density 0.3, and the sliver is then subdivided
    into 2 and into 10 pieces **of the same density** — a re-partition that changes no data at all.
    Returns the estimate at each of the three partitions.

    ⚠ NECESSARY, not sufficient: a constant estimator passes too. Accuracy is settled by the panel.
    """
    out = []
    for k in (1, 2, 10):
        dens = np.concatenate([[3.0], np.full(k, 0.3)])
        opp = np.concatenate([[990.0], np.full(k, 10.0 / k)])
        if not weighted:
            opp = np.ones_like(opp)
        seg = np.zeros(dens.size, dtype=np.int64)
        m, _ = _power_mean(dens, opp, seg, 1, mode)
        out.append(float(m[0]))
    return tuple(out)  # type: ignore[return-value]


def weight_vs_truth(w: np.ndarray, truth: pd.DataFrame, index) -> dict:
    """How well does one candidate VECTOR agree with the per-transcript truth?

    ⛔⛔ **SCORED OVER THE ANNOTATED TRANSCRIPTS ONLY — 8,750 of the 15,669 on the EM's axis.** The other
    6,919 are SYNTHETIC nascent entities, and the solver already refuses them prior mass whatever this
    vector says (`em_solver.cpp` skips ``is_synthetic`` in both the normaliser and the distribution;
    owner ruling `2fa21540`, *absent until proven otherwise*). Pooling them in read
    ``silent_kept_zero`` **0.690** on the first stage-6 arm when the figure over the population the
    solver can actually act on is the only one that means anything.

    ⚠ ``build_weights`` deliberately does NOT zero them itself. The weight answers "how many unspliced
    fragments would this entity emit", which is a real number for a shadow span; that it may not RECEIVE
    prior is a solver ruling, and putting one ruling in two homes is how the two drift apart.

    ⛔ The two failure directions are reported SEPARATELY and must not be pooled. Inventing weight on a
    silent transcript is the graveyard's `g00` failure — every one of the eleven refused mechanisms
    lifted an evidence-free object off zero — while zeroing an expressed one is unrecoverable, since a
    zero weight is EXACTLY absorbing (measured 0.0000 in all four mode × warm-start combinations).
    """
    from scipy.stats import pearsonr, spearmanr

    t_index = dict(zip(index.t_df["t_id"].to_numpy(), index.t_df["t_index"].to_numpy(), strict=True))
    tru = np.zeros(int(index.num_transcripts), dtype=np.float64)
    for tid, n in zip(truth["transcript_id"], truth["n_unspliced"], strict=True):
        i = t_index.get(str(tid))
        if i is not None:
            tru[int(i)] += float(n)

    annotated = ~index.t_df.sort_values("t_index")["is_synthetic"].to_numpy(bool)
    tru, w = tru[annotated], np.asarray(w, dtype=np.float64)[annotated]
    expressed = tru > 0
    return {
        "spearman": float(spearmanr(tru, w).statistic),
        "pearson": float(pearsonr(np.log2(tru + 1.0), np.log2(w + 1.0)).statistic),
        "n_annotated": int(annotated.sum()),
        "n_synthetic_excluded": int((~annotated).sum()),
        "n_expressed": int(expressed.sum()),
        "n_silent": int((~expressed).sum()),
        "silent_kept_zero": float((w[~expressed] == 0.0).mean()) if (~expressed).any() else float("nan"),
        "expressed_lost": float((w[expressed] == 0.0).mean()) if expressed.any() else float("nan"),
    }


def main() -> int:
    """⭐ Standalone, this instrument runs its own FALSIFICATION and nothing else.

    ⚠ It deliberately has no data mode. Scoring a weighting function end to end belongs to
    `quant_accuracy.py`, which already owns the arms, the truth join and the noise floor — a second
    scorer here is how a baseline and a ceiling drift apart.
    """
    ap = argparse.ArgumentParser(description="the re-partition falsification for the weight family")
    ap.parse_args()

    print("  ⭐ RE-PARTITION INVARIANCE — a 10 bp sliver subdivided 1 / 2 / 10 ways, data unchanged")
    print(f"    {'mode':12s} {'OPPORTUNITY-WEIGHTED':>34s}   {'UNWEIGHTED (the artefact)':>34s}")
    ok_all = True
    for mode in MODES:
        w = repartition_invariance(mode, weighted=True)
        u = repartition_invariance(mode, weighted=False)
        ok = max(w) - min(w) < 1e-9
        ok_all &= ok
        print(f"    {'✔' if ok else '✘'} {mode:10s} "
              f"{w[0]:10.4f} {w[1]:10.4f} {w[2]:10.4f}   "
              f"{u[0]:10.4f} {u[1]:10.4f} {u[2]:10.4f}   drift {max(u) / max(min(u), 1e-12):5.2f}x")
    print(f"    {'ALL WEIGHTED FORMS INVARIANT' if ok_all else '⛔ A WEIGHTED FORM DRIFTED'}")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
