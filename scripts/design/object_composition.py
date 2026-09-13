"""Must ψ's Beta reference be one library-wide number, or can each object supply its own?

ψ's composition reference ``a·log f_g + b·log(1−f_g)`` is ``Beta(a, b)`` in ``f_g``, so ``a`` and
``b`` are pseudo-counts with a mean ``a/(a+b)``, and the shipped ½ asserts every object is half gDNA.
This instrument measures what that mean should be instead, with no solver, no EM and nothing patched in
`src/`: the target (the fragment-weighted ``f_lib`` against the object-weighted mean of the per-object
truth, the one with a prior's denominator); whether it is reachable from the payload alone (no
deconvolved array is read here, which is what makes an estimate non-circular); and the per-object
prior mean ``m_i = rho_g·E_g / (rho_g·E_g + rho_r·E_r)`` built from two pooled gDNA densities (the
off-target anchors and the in-gene ``exon|intron`` anchor, whose ratio is the capture enrichment) and
a shrunk sj-flux RNA density gated to where mature RNA can be. Each arm is scored as the fragments the
prior misplaces if believed outright, ``Σ|m_i − f_g,i|·M_i``, against the shipped ½ and per stratum,
with the truth from the origin-split oracle cache (sum-to-full gated). Every slot takes exactly one
of seven strata (:func:`strata`, asserted), the boundary axis split on whether mature RNA can cross.
Read the zero-control rows first (the gDNA-side estimator must read 0.0000 there) and the capture-ON
rows second. `PURE_GDNA_STRATA` includes ``R intron``, which the panel's sparse nascent RNA
contaminates; the shipped background pools intergenic only, so that is this instrument's anchor.

Also a library: `calibration_oracle.py`, `calibration_walk.py` and `total_abundance_audit.py`
import `strata`, `slot_counts`, `_scope`, `_SELECTIONS`, `PVO` and
the two defaults; its own tables are not the yardstick for a mechanism, `vertex_ceiling.py` is.

Usage::

    python scripts/design/object_composition.py                                # the whole ladder
    python scripts/design/object_composition.py --conditions <name> --oracle-cache <dir>
    python scripts/design/object_composition.py --self-test                    # no I/O
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


from _shared import sibling  # noqa: E402


PVO = sibling("prior_vs_oracle.py")

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import (  # noqa: E402
    build_region_geometry,
    build_region_statics,
    g1_locked,
)
from rigel.calibration.signature import (  # noqa: E402
    BIT_EXON_NEG,
    BIT_EXON_POS,
    RegionType,
    coarse_type_array,
)
from rigel.calibration.splice_graph import (  # noqa: E402
    build_boundary_flags_array,
    build_sj_geometry_arrays,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from calibration._oracle import ORIGINS, OracleTruth  # noqa: E402

DEFAULT_SUITE = PVO.DEFAULT_SUITE
DEFAULT_INDEX = PVO.DEFAULT_INDEX

_EPS = 1.0e-12

#: Either strand's exon bit. The union over strands is deliberate: a boundary is "exonic" on a flank
#: if any transcript has an exon there, and the per-strand question is answered by ``mrna_active``,
#: which :func:`strata` cross-checks against.
_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG

#: The seven populations, in report order: mutually exclusive and exhaustive over the chain, asserted
#: in :func:`strata` rather than promised here, because a stratum table whose rows overlap
#: double-counts the object weight it exists to compute and nothing downstream could tell.
#: The boundary axis is split by whether mature RNA can cross, not by whether a sj attaches: a pool
#: keyed on "a sj attaches here" lumps in ``exon|exon`` boundaries, where an alternative splice site
#: sits inside a contiguous exonic stretch and mature RNA crosses freely.
STRATA = (
    "R intergenic",
    "R intron",
    "R exon",
    "B exon|intron",
    "B exon|exon",
    "B intron|intron",
    "B gene edge",
)

#: The strata whose composition is knowable with no deconvolution, and they are not equally free:
#: ``R intergenic`` is covered by no transcript at all, so ``f_g = 1`` holds at any nascent level;
#: ``R intron`` needs unspliced nascent RNA to be sparse, and the panel's sparse nascent RNA
#: contaminates it by a factor that scales with the RNA:gDNA ratio (worst where gDNA is scarce). The
#: shipped background (``density_deconv.fit_intron_background``) pools intergenic only, so a number
#: this instrument reports off this pool is the instrument's and not shipped behaviour; whether to drop
#: ``R intron`` is a measurement-design decision, so it is recorded here rather than silently changed.
#: ``est_bg_intergenic_only`` is printed beside ``est_bg`` for that reason.
PURE_GDNA_STRATA = ("R intergenic", "R intron")

#: The gDNA anchor that sits inside genes, and therefore on-target under hybrid capture, which the
#: intergenic anchor structurally cannot be.
ONTARGET_GDNA_STRATUM = "B exon|intron"


def strata(chain, statics, geometry, region_arrays) -> dict:
    """Per-slot stratum label and the masks the derivation needs, all from the annotation alone.

    The partition is asserted: every slot must take exactly one label, because the object-weighted
    mean this file computes is a mean over objects, so a slot counted twice or not at all silently
    reweights the target.

    The boundary classification is by exon-ness of the two flanks, because a boundary's
    unspliced-crossing population is the molecules that crossed it contiguously, which mature RNA can
    do only where the template is contiguous exon on both sides:

    ================  =================================================================================
    ``B exon|intron``  exactly one flank exonic, both flanks inside a gene: mature RNA cannot cross, it
                       splices, so near-pure gDNA under sparse unspliced nascent; in-gene, so the
                       on-target gDNA anchor under capture
    ``B exon|exon``    both flanks exonic, an alternative splice site inside a contiguous exonic
                       stretch; mature RNA crosses it freely, so not an anchor
    ``B intron|intron`` neither flank exonic, both inside a gene (adjacent introns of different
                       signature); off-target, nascent-only RNA
    ``B gene edge``    at least one flank intergenic, a TSS/TES interface; the ``g1_locked`` boundary
                       class, structurally pure gDNA on both strands
    ================  =================================================================================

    ``R intergenic`` is defined by the signature and cross-checked against ``g1_locked``, the predicate
    the solver itself pins on; ``B exon|intron`` is cross-checked against the solver's own
    ``mrna_active``. If either pair separates, this file's purity claim is about a different population
    than the one the solver reasons over, and it raises.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    is_region, is_boundary = kind == REGION, kind == BOUNDARY
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    n_regions = sig.shape[0]
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    slot_type = np.where(is_region, rtype[np.clip(obj, 0, max(n_regions - 1, 0))], -1)
    locked = g1_locked(np.asarray(statics.free_pos, bool), np.asarray(statics.free_neg, bool))

    # the two flanks' signatures, through the chain's own adjacency (a BOUNDARY always has a REGION on
    # both sides, so the ``-1`` branch has no cases — `build_region_statics` makes the same argument)
    slot_sig = np.where(is_region, sig[np.clip(obj, 0, max(n_regions - 1, 0))] if n_regions else 0, 0)
    left = np.clip(np.asarray(chain.left), 0, max(int(chain.n_slots) - 1, 0))
    right = np.clip(np.asarray(chain.right), 0, max(int(chain.n_slots) - 1, 0))
    sig_l = np.where(is_boundary, slot_sig[left], 0)
    sig_r = np.where(is_boundary, slot_sig[right], 0)
    exon_l = (sig_l & _EXON_BITS) != 0
    exon_r = (sig_r & _EXON_BITS) != 0
    gene_both = is_boundary & (sig_l != 0) & (sig_r != 0)

    label = np.full(int(chain.n_slots), "", dtype=object)
    label[is_region & (slot_type == int(RegionType.INTERGENIC))] = "R intergenic"
    label[is_region & (slot_type == int(RegionType.INTRON))] = "R intron"
    label[is_region & (slot_type == int(RegionType.EXON))] = "R exon"
    label[is_boundary & ~gene_both] = "B gene edge"
    label[gene_both & exon_l & exon_r] = "B exon|exon"
    label[gene_both & (exon_l ^ exon_r)] = "B exon|intron"
    label[gene_both & ~exon_l & ~exon_r] = "B intron|intron"

    counted = sum(int(np.sum(label == s)) for s in STRATA)
    if counted != int(chain.n_slots):
        raise AssertionError(
            f"the stratum labels do not partition the chain: {counted:,} labelled of "
            f"{int(chain.n_slots):,} slots. Every slot must take exactly one label, or the "
            "object-weighted mean is a mean over the wrong denominator."
        )
    if not np.array_equal(label == "R intergenic", is_region & locked):
        raise AssertionError(
            f"`intergenic & REGION` ({int(np.sum(label == 'R intergenic')):,}) and `g1_locked & REGION` "
            f"({int(np.sum(is_region & locked)):,}) have SEPARATED on this index. They are the same "
            "population by construction — no transcript covers an intergenic region, so neither RNA "
            "strand is admissible — and this file's structural `f_g = 1` claim rests on that."
        )
    # `mrna_active_s` is the solver's own "contiguous exon on both flanks" gate, i.e. "mature RNA of
    # strand s may cross here". At an `exon|intron` boundary one flank carries no exon bit at all, so
    # it must be False on both strands; if that ever fails, the signature semantics moved and the
    # near-pure-gDNA claim is about a different population.
    mature_can_cross = np.asarray(statics.mrna_active_pos, bool) | np.asarray(
        statics.mrna_active_neg, bool
    )
    bad = (label == ONTARGET_GDNA_STRATUM) & mature_can_cross
    if bad.any():
        raise AssertionError(
            f"{int(bad.sum()):,} `{ONTARGET_GDNA_STRATUM}` boundaries report `mrna_active`, i.e. the "
            "solver thinks mature RNA may cross them contiguously. The anchor's whole claim is that it "
            "cannot, so this classification and the solver's disagree."
        )
    return {
        "label": label,
        "locked": locked,
        # on-target means "touches an exon": annotation-derived, no threshold and no capture
        # detection. It is the axis hybrid capture enriches along, so it is the axis the gDNA density
        # is allowed to differ across.
        "on_target": np.where(
            is_region,
            (np.where(is_region, sig[np.clip(obj, 0, max(n_regions - 1, 0))], 0) & _EXON_BITS) != 0,
            exon_l | exon_r,
        ),
        "has_sj": np.asarray(geometry.eff_sj, np.float64).sum(1) > 0.0,
    }


def slot_counts(payload, region_arrays, chain) -> np.ndarray:
    """One payload's unspliced/contained count per slot, the mixture ψ deconvolves and nothing else:
    ``region_contained`` at a REGION, ``boundary_unspliced`` at a BOUNDARY, exactly the populations
    :attr:`RegionGeometry.unspliced_count` carries, so a truth built from the origin partitions and an
    estimate built from the full payload are on one basis.
    """
    sub = CalibrationSubstrate.from_payload(payload, region_arrays)
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    out = np.zeros(int(chain.n_slots), np.float64)
    r, b = kind == REGION, kind == BOUNDARY
    out[r] = np.asarray(sub.region_contained.count, np.float64).sum(1)[obj[r]]
    out[b] = np.asarray(sub.boundary_unspliced.count, np.float64).sum(1)[obj[b]]
    return out


def pooled_density(mass: np.ndarray, eff: np.ndarray, select: np.ndarray) -> float:
    """``sum(mass) / sum(eff)`` over a selected population, the ratio of sums and never the mean of
    ratios (TRAPS: a-mean-of-ratios-inherits-the-partition). ``0.0`` when the population has no
    opportunity, which is the honest answer and not a floored division."""
    e = float(np.sum(np.asarray(eff, np.float64)[select]))
    if e <= 0.0:
        return 0.0
    return float(np.sum(np.asarray(mass, np.float64)[select]) / e)


def neighbour_sj_density(chain, geometry) -> np.ndarray:
    """Per-slot certified-RNA density from the sj flux at the slot and its two chain neighbours,
    pooled as ``sum(count) / sum(E)``; returns ``(density, has_flux)``.

    The flux is a density on the spliced template, so it is on the same footing as an exon REGION's
    contained RNA opportunity and not on that of a BOUNDARY's unspliced-crossing opportunity; table ③
    measures both. The neighbour set is the chain's own adjacency: a slot's neighbours are ``i-1`` and
    ``i+1`` and a reference terminal links to ``-1``.
    """
    sj_count = np.asarray(geometry.sj_count, np.float64).sum(1)
    eff_sj = np.asarray(geometry.eff_sj, np.float64).sum(1)
    num, den = sj_count.copy(), eff_sj.copy()
    for side in (np.asarray(chain.left, np.int64), np.asarray(chain.right, np.int64)):
        ok = side >= 0
        num[ok] += sj_count[side[ok]]
        den[ok] += eff_sj[side[ok]]
    return np.where(den > 0.0, num / np.maximum(den, _EPS), 0.0), den > 0.0


def object_weighted_mean(f_true: np.ndarray, live: np.ndarray) -> float:
    """The target: the mean per-object composition over objects that have any mass at all.

    Objects with no mass are excluded because an empty object has no true composition to average
    (``0/0``), and folding a fabricated value in would move the target by the share of the genome that
    happens to be empty at this gDNA level, which at the zero control is most of it.
    """
    return float(np.mean(np.asarray(f_true, np.float64)[live])) if live.any() else float("nan")


def measure_condition(index, region_arrays, sj, boundary_flags, suite: Path, oracle_cache: Path,
                      condition: str) -> dict:
    """One condition: the target, the stratum census, and every estimator. One JSON-able row."""
    start = time.perf_counter()
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    lift: dict = {}
    kw = calibration_inputs(cache, index, lift_out=lift)
    # the drained frame, the one calibration reads
    payload = kw["payload"]
    chain = build_region_chain(
        payload.ref_region_offsets, payload.ref_boundary_offsets
    )
    statics = build_region_statics(chain, region_arrays, boundary_flags)
    geometry = build_region_geometry(
        chain,
        CalibrationSubstrate.from_payload(payload, region_arrays),
        region_arrays,
        sj,
        kw["gdna_fl_pmf"],
        kw["rna_fl_pmf"],
        None,
    )
    cls = strata(chain, statics, geometry, region_arrays)
    label, locked, has_sj, on_target = (
        cls["label"], cls["locked"], cls["has_sj"], cls["on_target"]
    )

    # Through ``OracleTruth`` rather than the raw payloads, so sum-to-full runs as a hard gate on every
    # condition: a cached partition that does not reconstruct the scan calibration read is a silently
    # wrong truth source. The full payload is the scan cache's, never the oracle cache's ``_main``, so
    # sum-to-full validates the truth partition against the exact array object the estimator reads,
    # and this instrument never reads a directory ``pass0_vs_oracle.measure_condition`` writes.
    root = Path(oracle_cache) / condition
    parts = {k: read_scan_cache(root / k, index).payload for k in ORIGINS}
    truth_oracle = OracleTruth.from_cached_parts(payload, parts, lift)
    parts = truth_oracle.parts  # the drained partitions — the same frame as `payload`

    n_g = slot_counts(parts["gdna"], region_arrays, chain)
    n_r = slot_counts(parts["mrna"], region_arrays, chain) + slot_counts(
        parts["nrna"], region_arrays, chain
    )
    mass = n_g + n_r
    live = mass > 0.0
    f_true = np.where(live, n_g / np.maximum(mass, _EPS), 0.0)

    eff_g = np.asarray(geometry.eff_gdna, np.float64)
    eff_r = np.asarray(geometry.eff_rna, np.float64)

    # The estimator reads the full payload's own per-slot totals, never the origin partitions: the
    # non-circularity claim made structural.
    est_mass = np.zeros(int(chain.n_slots), np.float64)
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    full_sub = CalibrationSubstrate.from_payload(payload, region_arrays)
    est_mass[kind == REGION] = np.asarray(full_sub.region_contained.count, np.float64).sum(1)[
        obj[kind == REGION]
    ]
    est_mass[kind == BOUNDARY] = np.asarray(full_sub.boundary_unspliced.count, np.float64).sum(1)[
        obj[kind == BOUNDARY]
    ]
    if not np.allclose(est_mass, mass, rtol=0.0, atol=1e-6):
        raise AssertionError(
            "the full payload's per-slot totals and the summed origin partitions disagree. "
            "`OracleTruth.from_parts` has just validated sum-to-full on every bank, so the two "
            "projections onto the chain are reading different populations."
        )

    def gdna_side(select) -> tuple[float, float]:
        rho = pooled_density(est_mass, eff_g, select)
        est = np.clip(rho * eff_g / np.maximum(est_mass, _EPS), 0.0, 1.0)
        return rho, object_weighted_mean(est, live)

    rho_r, has_flux = neighbour_sj_density(chain, geometry)
    est_sj = np.clip(1.0 - rho_r * eff_r / np.maximum(est_mass, _EPS), 0.0, 1.0)

    anchors = np.isin(label, list(PURE_GDNA_STRATA))
    rho_bg, mean_bg = gdna_side(anchors)
    _, mean_bg_ig = gdna_side(label == "R intergenic")

    # ── the per-object prior mean `m_i` ──
    #
    #   m_i = rho_g,i * E_g,i / ( rho_g,i * E_g,i + rho_r,i * E_r,i )
    #
    # A library-wide scalar is the special case where both densities are constant and the geometry is
    # ignored; every term but the two densities is per-object and exactly known, so even a class-pooled
    # density yields a per-object prior mean, which is the claim under test. Two gDNA densities, split
    # on an annotation-derived axis (does this object touch an exon?) rather than a detected capture
    # flag: hybrid capture enriches along exactly that axis, so their ratio is the enrichment factor,
    # measured rather than declared, with no threshold and no switch.
    on_anchor = label == ONTARGET_GDNA_STRATUM
    rho_g_off = pooled_density(est_mass, eff_g, anchors)
    rho_g_on = pooled_density(est_mass, eff_g, on_anchor)
    rho_g_per_object = np.where(on_target, rho_g_on, rho_g_off)

    # the TRUE per-class densities, so the FORM can be priced apart from the ESTIMATOR
    def true_rho_g(select) -> float:
        e = float(np.sum(eff_g[select]))
        return float(np.sum(n_g[select]) / e) if e > 0.0 else 0.0

    def true_rho_r(select) -> float:
        e = float(np.sum(eff_r[select]))
        return float(np.sum(n_r[select]) / e) if e > 0.0 else 0.0

    true_g_on, true_g_off = true_rho_g(live & on_target), true_rho_g(live & ~on_target)
    true_r_on, true_r_off = true_rho_r(live & on_target), true_rho_r(live & ~on_target)

    def prior_mean(rho_g, rho_r) -> np.ndarray:
        g, r = rho_g * eff_g, rho_r * eff_r
        tot = g + r
        return np.where(tot > 0.0, g / np.maximum(tot, _EPS), 0.5)

    # The mature-RNA gate, the same predicate the boundary strata use. ``rho_r`` from the sj flux is a
    # density of mature molecules on the spliced template; handing it to a slot mature RNA cannot
    # occupy (an intron REGION, an ``exon|intron`` boundary, a gene edge) subtracts RNA that is not
    # there and calls gDNA RNA. ``mrna_active_s`` is the shipped predicate for exactly this on both
    # axes: a REGION's own exon bit, a BOUNDARY's contiguous exon on both flanks. Where mature cannot
    # be, what remains is nascent, which this design assumes sparse.
    mature_here = np.asarray(statics.mrna_active_pos, bool) | np.asarray(
        statics.mrna_active_neg, bool
    )
    # the library-wide certified-RNA density on the spliced template: Sum(flux) / Sum(E_sj), the ratio of
    # sums and never a mean of ratios (TRAPS: a-mean-of-ratios-inherits-the-partition)
    _sjc = np.asarray(geometry.sj_count, np.float64).sum(1)
    _esj = np.asarray(geometry.eff_sj, np.float64).sum(1)
    rho_r_pooled = float(_sjc.sum() / max(float(_esj.sum()), _EPS))

    # The certified spliced crossing is a subtraction, not a bound. ``boundary_spliced`` is a separate
    # bank from ``boundary_unspliced`` (the same molecules split by whether they used a sj elsewhere),
    # so the contiguous RNA crossing a boundary is ``unspliced_RNA + S`` and
    #
    #     unspliced_RNA = rho_r * E_r - S      =>      f_g = 1 - (rho_r * E_r - S) / M
    #
    # ``f_g <= 1 - S/M`` is false (S is not inside M), so S is not in the arm ladder and table ⑦ keeps
    # the identity under measurement instead.
    spliced = np.zeros(int(chain.n_slots), np.float64)
    spliced[kind == BOUNDARY] = np.asarray(full_sub.boundary_spliced.count, np.float64).sum(1)[
        obj[kind == BOUNDARY]
    ]

    # The local RNA density, three ways. The per-object sj flux is at the right level and noisy per
    # object, and setting ``rho_r = 0`` wherever no sj is in reach sends an RNA-rich exon to
    # ``m_i = 1``, a coverage artefact rather than a measurement. The fallback is the population rate,
    # and the blend is one pseudo-observation of it, the same "one pseudo-object of ignorance"
    # convention ``fit_landscape`` uses, so no constant is introduced.
    sj_num = np.asarray(geometry.sj_count, np.float64).sum(1)
    sj_den = np.asarray(geometry.eff_sj, np.float64).sum(1)
    num, den = sj_num.copy(), sj_den.copy()
    for side in (np.asarray(chain.left, np.int64), np.asarray(chain.right, np.int64)):
        ok = side >= 0
        num[ok] += sj_num[side[ok]]
        den[ok] += sj_den[side[ok]]
    #: the mean sj opportunity carried by one sj-bearing boundary: the weight of one pseudo-observation
    e_one = float(sj_den[sj_den > 0.0].mean()) if np.any(sj_den > 0.0) else 0.0
    rho_r_fallback = np.where(has_flux, rho_r, rho_r_pooled)
    rho_r_shrunk = (num + rho_r_pooled * e_one) / np.maximum(den + e_one, _EPS)

    def gated(x):
        """``rho_r`` is a mature density on the spliced template: zero where mature RNA cannot be."""
        return np.where(mature_here, x, 0.0)

    #: the arm ladder, each an `m_i` and each scored the same way. `shipped` is the constant ½ ψ
    #: carries, the baseline every ratio is against, and `TRUTH` is the class-pooled ceiling.
    #: `structural` lets the reference speak only where the annotation determines the answer (mature
    #: RNA cannot be here) and leaves ψ's ½ everywhere else.
    m_pooled = prior_mean(rho_g_per_object, gated(rho_r_pooled))
    m_arms = {
        "shipped": np.full(int(chain.n_slots), 0.5),
        "TRUTH": prior_mean(
            np.where(on_target, true_g_on, true_g_off), np.where(on_target, true_r_on, true_r_off)
        ),
        "structural": np.where(mature_here, 0.5, m_pooled),
        "pooled": m_pooled,
        "flux+fallbk": prior_mean(rho_g_per_object, gated(rho_r_fallback)),
        "flux+shrunk": prior_mean(rho_g_per_object, gated(rho_r_shrunk)),
        # the deconvolution: `f_g = rho_g * E_g / M`, RNA as the residual and never predicted. It is
        # the peak of the shipped `density_lambda_factor` ("peaked at f_g = rho_bg/rho_obs", with
        # `rho_obs = M/E_g`), so the arm measures the location the shipped factor already carries at
        # the one stratum it is switched on for. It needs no RNA density: gDNA is near-uniform and
        # predictable pre-solve, RNA is whatever mass the gDNA deconvolve leaves behind.
        "deconvolve": np.clip(rho_g_per_object * eff_g / np.maximum(est_mass, _EPS), 0.0, 1.0),
        # the hybrid: pin where the annotation determines the answer, deconvolve where it does not.
        # The structural strata are exactly ``~mature_here``, and the deconvolve's use of the observed
        # ``M`` makes it strictly worse there (a downward Poisson fluctuation reads as RNA).
        "pin+deconvolve": np.where(
            mature_here,
            np.clip(rho_g_per_object * eff_g / np.maximum(est_mass, _EPS), 0.0, 1.0),
            1.0,
        ),
        "deconvolve, 1 rho_g": np.clip(
            rho_g_off * eff_g / np.maximum(est_mass, _EPS), 0.0, 1.0
        ),
    }
    row_m = {}
    for name, m in m_arms.items():
        d = np.abs(m - f_true)
        row_m[name] = {
            "abs_err_frags": float(np.sum(d[live] * mass[live])),
            "mwae": float(np.sum(d[live] * mass[live]) / max(float(mass[live].sum()), _EPS)),
            "objw": object_weighted_mean(m, live),
            "per_stratum_frags": {
                s: float(np.sum(d[live & (label == s)] * mass[live & (label == s)]))
                for s in STRATA
            },
        }

    # ── the runaway bound, with the likelihood removed ──
    #
    # A post-solve re-estimate of the on-target gDNA density from solved exons is one scalar measured
    # on exons and applied to exons, structurally the positive-feedback loop that makes a library-wide
    # ``f_lib`` rule inadmissible: a higher rho_g raises every exonic ``m_i``, which raises the gDNA
    # mass attributed to exons, which raises rho_g. Iterating that map with each object's own
    # likelihood removed is a strict upper bound on the feedback (the prior believed outright, with no
    # data pulling against it): if this converges the real refit loop converges a fortiori, which is
    # what makes a solver-free answer admissible here at all.
    exonic = live & on_target
    rho_r_here = np.where(mature_here, rho_r_pooled, 0.0)

    def _runaway(start: float, iters: int = 24) -> list[float]:
        rho, traj = float(start), []
        denom = float(np.sum(eff_g[exonic]))
        for _ in range(iters):
            g = rho * eff_g
            tot = g + rho_r_here * eff_r
            m = np.where(tot > 0.0, g / np.maximum(tot, _EPS), 0.5)
            rho = float(np.sum((m * est_mass)[exonic]) / denom) if denom > 0.0 else 0.0
            traj.append(rho)
        return traj

    true_on = true_rho_g(exonic)
    #: four starts spanning three decades around the measured anchor: if they meet, the map has one
    #: attracting fixed point and the starting value does not matter.
    runaway = {
        k: _runaway(v)
        for k, v in (
            ("anchor", rho_g_on),
            ("10x low", rho_g_on / 10.0),
            ("10x high", rho_g_on * 10.0),
            ("truth", true_on),
        )
    }

    # ── the contiguous-RNA identity, kept under measurement ──
    # ``rho_r * E_r = unspliced_RNA + S`` at a BOUNDARY. Scored as: does the neighbouring sj flux
    # recover the true total contiguous RNA density ``(n_r + S)/E_r``, the quantity a corrected
    # estimator would need?
    contig_true = np.where(eff_r > 0.0, (n_r + spliced) / np.maximum(eff_r, _EPS), 0.0)
    has_contig = live & mature_here & has_flux & (contig_true > 0.0)

    row = {
        "condition": condition,
        "stratum": list(PVO.stratum(condition)),
        "seconds": time.perf_counter() - start,
        "n_slots": int(chain.n_slots),
        "n_live": int(live.sum()),
        "mass_total": float(mass.sum()),
        # ── the TARGET, three ways ──
        "target_objw": object_weighted_mean(f_true, live),
        # the same mean over the objects ψ actually solves a composition on: `g1_locked` slots are
        # pinned, so the reference never moves them and their share of the target is inert.
        "target_objw_unlocked": object_weighted_mean(f_true, live & ~locked),
        "f_lib": float(n_g.sum() / max(mass.sum(), _EPS)),
        # ── the ESTIMATORS ──
        "est_bg": mean_bg,
        "est_bg_intergenic_only": mean_bg_ig,
        "est_sj": object_weighted_mean(est_sj, live),
        "rho_gdna_est": rho_bg,
        "rho_gdna_true": float(n_g.sum() / max(np.sum(eff_g[live]), _EPS)),
        "sj_coverage_objects": float(np.sum(live & has_flux) / max(int(live.sum()), 1)),
        "sj_coverage_mass": float(np.sum(mass[live & has_flux]) / max(float(mass[live].sum()), _EPS)),
        # ── the two gDNA densities and the ENRICHMENT they imply ──
        "rho_g_off": rho_g_off,
        "rho_g_on": rho_g_on,
        "rho_g_off_true": true_g_off,
        "rho_g_on_true": true_g_on,
        "enrichment_est": rho_g_on / rho_g_off if rho_g_off > 0.0 else float("nan"),
        "enrichment_true": true_g_on / true_g_off if true_g_off > 0.0 else float("nan"),
        # ── the per-object prior mean, the point of the file ──
        "m_arms": row_m,
        # ── stage 0c: the runaway bound ──
        "runaway_fixed_points": {k: v[-1] for k, v in runaway.items()},
        "runaway_true_on": true_on,
        # ── the contiguous-RNA identity: rho_r * E_r = unspliced_RNA + S ──
        "contig_objects": int(has_contig.sum()),
        "contig_mass_share": float(
            np.sum(mass[has_contig]) / max(float(mass[live].sum()), _EPS)
        ),
        "contig_agg": (
            float(np.sum((rho_r * eff_r)[has_contig]) / max(np.sum((n_r + spliced)[has_contig]), _EPS))
            if has_contig.any() else float("nan")
        ),
        "spliced_share_of_contig": (
            float(np.sum(spliced[has_contig]) / max(float(np.sum((n_r + spliced)[has_contig])), _EPS))
            if has_contig.any() else float("nan")
        ),
        "strata": {},
    }

    for s in STRATA:
        m = label == s
        ml = m & live
        row["strata"][s] = {
            "objects": int(m.sum()),
            "live": int(ml.sum()),
            "mass": float(mass[m].sum()),
            "object_share": float(int(ml.sum()) / max(int(live.sum()), 1)),
            "mass_share": float(mass[m].sum() / max(float(mass.sum()), _EPS)),
            "mean_fg": float(np.mean(f_true[ml])) if ml.any() else float("nan"),
            "median_fg": float(np.median(f_true[ml])) if ml.any() else float("nan"),
            "p10_fg": float(np.percentile(f_true[ml], 10)) if ml.any() else float("nan"),
        }

    # ── ③ the sj flux as an RNA density, on two bases (the spliced template, the crossing opportunity) ──
    def flux_vs_truth(select) -> dict:
        if not select.any():
            return {"n": 0, "aggregate": float("nan"), "median": float("nan")}
        predicted = float(np.sum((rho_r * eff_r)[select]))
        actual = float(np.sum(n_r[select]))
        rho_true_r = np.where(eff_r > 0.0, n_r / np.maximum(eff_r, _EPS), 0.0)
        ok = select & (rho_true_r > 0.0)
        ratio = rho_r[ok] / np.maximum(rho_true_r[ok], _EPS)
        w = mass[ok]
        order = np.argsort(ratio)
        cum = np.cumsum(w[order]) / max(float(w.sum()), _EPS)
        return {
            "n": int(select.sum()),
            "aggregate": predicted / actual if actual > 0.0 else float("inf"),
            "median": float(ratio[order][np.searchsorted(cum, 0.50)]) if ok.any() else float("nan"),
        }

    row["flux_at_sj_boundary"] = flux_vs_truth(live & has_sj & (eff_r > 0.0))
    row["flux_at_exon_region"] = flux_vs_truth(live & (label == "R exon") & has_flux & (eff_r > 0.0))
    return row


# ── reporting ────────────────────────────────────────────────────────────────────────────────────

#: The 0.8.0 scope, stamped on every row rather than left to the reader: three strata are the
#: development target and unstranded x capture-ON is deferred but reported. Same table as
#: ``calibration_vs_oracle.py``'s; a reader ranking on a stratum that is not a target inverts the order.
_SCOPE = {
    ("stranded", "capture OFF"): "IN SCOPE",
    ("stranded", "capture ON"): "IN SCOPE",
    ("unstranded", "capture OFF"): "IN SCOPE",
    ("unstranded", "capture ON"): "DEFERRED",
}


def _scope(condition: str) -> str:
    return "CONTROL" if PVO.is_zero_gdna(condition) else _SCOPE[PVO.stratum(condition)]


#: Every selection the per-stratum tables print, in order: one list, so a stratum cannot appear on some
#: tables and not others. The zero control is its own row and is never folded into a stratum: its truth
#: is exactly 0, so every gDNA fragment there is a false positive with nothing to cancel it.
_SELECTIONS = (
    *(
        (
            f"{' x '.join(st)}  [{_SCOPE[st]}]",
            (lambda c, st=st: PVO.stratum(c) == tuple(st) and not PVO.is_zero_gdna(c)),
        )
        for st in _SCOPE
    ),
    ("⛔ g00 ZERO-gDNA control (all strata)", PVO.is_zero_gdna),
)


def _f(x, width=8, places=4) -> str:
    return f"{'—':>{width}}" if x is None or not np.isfinite(x) else f"{x:>{width}.{places}f}"


def report(rows: list[dict]) -> None:
    rows = sorted(rows, key=lambda r: r["condition"])

    print(f"\n{'=' * 124}")
    print("① THE TARGET — the OBJECT-weighted mean composition against the FRAGMENT-weighted f_lib")
    print(f"{'=' * 124}")
    print(
        f"{'condition':<40} {'scope':<9} {'objw':>8} {'objw*':>8} {'f_lib':>8} "
        f"{'objw−f_lib':>11} {'n_live':>8} {'a,b for vertex_ceiling':>26}"
    )
    print("-" * 124)
    for r in rows:
        a = r["target_objw"]
        print(
            f"{r['condition']:<40} {_scope(r['condition']):<9} {_f(a)} "
            f"{_f(r['target_objw_unlocked'])} {_f(r['f_lib'])} "
            f"{_f(a - r['f_lib'], 11)} {r['n_live']:>8,} "
            f"{f'ref_c={a:.4f},{1 - a:.4f}':>26}"
        )
    print("⚠ `objw*` excludes `g1_locked` slots — the ones psi pins, where the reference never moves.")

    print(f"\n{'=' * 124}")
    print("② THE STRATUM CENSUS — objects, mass and TRUE composition per population")
    print(f"{'=' * 124}")
    for r in rows:
        print(f"\n{r['condition']}   [{_scope(r['condition'])}]   objw = {r['target_objw']:.4f}")
        print(
            f"   {'stratum':<14} {'objects':>9} {'live':>9} {'mass':>13} {'obj share':>10} "
            f"{'mass share':>11} {'mean fg':>9} {'median':>8} {'p10':>8}"
        )
        for s in STRATA:
            d = r["strata"][s]
            print(
                f"   {s:<14} {d['objects']:>9,} {d['live']:>9,} {d['mass']:>13,.0f} "
                f"{d['object_share']:>10.3f} {d['mass_share']:>11.3f} {_f(d['mean_fg'], 9)} "
                f"{_f(d['median_fg'])} {_f(d['p10_fg'])}"
            )

    print(f"\n{'=' * 124}")
    print("③ THE sj FLUX AS AN RNA DENSITY — at the BOUNDARY it attaches to, and at the EXON beside it")
    print(f"{'=' * 124}")
    print(
        f"{'condition':<40} {'scope':<9} | {'B sj n':>8} {'B sj agg':>9} {'B sj med':>9} "
        f"| {'exon n':>8} {'exon agg':>9} {'exon med':>9}"
    )
    print("-" * 124)
    for r in rows:
        b, e = r["flux_at_sj_boundary"], r["flux_at_exon_region"]
        print(
            f"{r['condition']:<40} {_scope(r['condition']):<9} | {b['n']:>8,} "
            f"{_f(b['aggregate'], 9, 3)} {_f(b['median'], 9, 3)} | {e['n']:>8,} "
            f"{_f(e['aggregate'], 9, 3)} {_f(e['median'], 9, 3)}"
        )
    print(
        "⭐ `agg` is Σ(rho_sj·E_r) ÷ Σ(TRUE RNA count) over the population — 1.000 is a perfect\n"
        "   recovery. ⛔ At the BOUNDARY the RNA divisor is the UNSPLICED-CROSSING opportunity, built\n"
        "   with unbounded reach, and the flux is a density on the SPLICED template: two templates, and\n"
        "   the product overstates by the factor in that column."
    )

    print(f"\n{'=' * 124}")
    print("④ THE ESTIMATORS — prior-free, scored against the TARGET in column ①")
    print(f"{'=' * 124}")
    print(
        f"{'condition':<40} {'scope':<9} {'objw':>8} | {'bg':>8} {'Δ':>8} | {'bg_ig':>8} "
        f"| {'sj':>8} {'Δ':>8} | {'rho_est':>9} {'rho_true':>9}"
    )
    print("-" * 124)
    for r in rows:
        t = r["target_objw"]
        print(
            f"{r['condition']:<40} {_scope(r['condition']):<9} {_f(t)} | "
            f"{_f(r['est_bg'])} {_f(r['est_bg'] - t)} | {_f(r['est_bg_intergenic_only'])} | "
            f"{_f(r['est_sj'])} {_f(r['est_sj'] - t)} | "
            f"{_f(r['rho_gdna_est'], 9, 6)} {_f(r['rho_gdna_true'], 9, 6)}"
        )
    print(
        "⛔ Read the CONTROL rows first: `bg` must read 0.0000 there and `sj` must not.\n"
        "⛔ Read the capture-ON rows second: `rho_est` against `rho_true` is where `bg` fails, because\n"
        "   the probes deplete the off-target anchors the pooled density is measured on."
    )

    print(f"\n{'=' * 124}")
    print("⑤ THE TWO gDNA DENSITIES — off-target anchors vs the in-gene `exon|intron` anchor")
    print(f"{'=' * 124}")
    print(
        f"{'condition':<40} {'scope':<9} | {'rho_off':>10} {'true':>10} {'x':>6} | "
        f"{'rho_on':>10} {'true':>10} {'x':>6} | {'enrich':>8} {'true':>8}"
    )
    print("-" * 124)
    for r in rows:
        ro, rot = r["rho_g_off"], r["rho_g_off_true"]
        rn, rnt = r["rho_g_on"], r["rho_g_on_true"]
        print(
            f"{r['condition']:<40} {_scope(r['condition']):<9} | {_f(ro, 10, 6)} {_f(rot, 10, 6)} "
            f"{_f(ro / rot if rot > 0 else np.nan, 6, 2)} | {_f(rn, 10, 6)} {_f(rnt, 10, 6)} "
            f"{_f(rn / rnt if rnt > 0 else np.nan, 6, 2)} | "
            f"{_f(r['enrichment_est'], 8, 2)} {_f(r['enrichment_true'], 8, 2)}"
        )
    print(
        "⭐ `enrich` is measured, not detected: it is ~1 without probes and large with them, so hybrid\n"
        "   capture needs no flag and no threshold. ⛔ The `x` columns are estimate ÷ truth — the\n"
        "   in-gene anchor is the one that must survive capture."
    )

    print(f"\n{'=' * 124}")
    print("⑥ ⭐⭐⭐ THE PER-OBJECT PRIOR MEAN `m_i` AGAINST TRUTH — Σ|m−f_g| in FRAGMENTS, ratio to ½")
    print(f"{'=' * 124}")
    arms = list(rows[0]["m_arms"])
    print(f"{'condition':<40} {'scope':<9} {'shipped Σ|Δ|':>13} " + " ".join(f"{a:>21}" for a in arms[1:]))
    print("-" * 124)
    for r in rows:
        base = r["m_arms"][arms[0]]["abs_err_frags"]
        cells = []
        for a in arms[1:]:
            v = r["m_arms"][a]["abs_err_frags"]
            cells.append(f"{v / base if base > 0 else float('nan'):>21.3f}")
        print(f"{r['condition']:<40} {_scope(r['condition']):<9} {base:>13,.0f} " + " ".join(cells))
    print("\n   per stratum, summed over the conditions of each scope (ratio to the shipped ½):")
    for sel_name, pred in _SELECTIONS:
        sel = [r for r in rows if pred(r["condition"])]
        if not sel:
            continue
        print(f"\n   {sel_name}")
        print(f"      {'stratum':<16} {'shipped Σ|Δ|':>13} " + " ".join(f"{a:>21}" for a in arms[1:]))
        for s in STRATA:
            b = sum(r["m_arms"][arms[0]]["per_stratum_frags"][s] for r in sel)
            cells = []
            for a in arms[1:]:
                v = sum(r["m_arms"][a]["per_stratum_frags"][s] for r in sel)
                cells.append(f"{v / b if b > 0 else float('nan'):>21.3f}")
            print(f"      {s:<16} {b:>13,.0f} " + " ".join(cells))
        b = sum(r["m_arms"][arms[0]]["abs_err_frags"] for r in sel)
        cells = [
            f"{sum(r['m_arms'][a]['abs_err_frags'] for r in sel) / b if b > 0 else float('nan'):>21.3f}"
            for a in arms[1:]
        ]
        print(f"      {'ALL':<16} {b:>13,.0f} " + " ".join(cells))
    print(
        "\n⛔ This is the PRIOR scored on its own, with no solver: `Σ|m_i − f_g,i|·M_i` is how many\n"
        "   fragments the prior misplaces if believed outright. `class-pooled TRUTH` prices the FORM\n"
        "   (is a per-class density plus per-object geometry enough?); `prior-free` prices the whole\n"
        "   thing; `prior-free, 1 rho_g` withholds the in-gene anchor, so it prices the capture split."
    )

    print(f"\n{'=' * 124}")
    print("⑦ THE CONTIGUOUS-RNA IDENTITY — rho_r·E_r = unspliced_RNA + S, where mature RNA can be")
    print(f"{'=' * 124}")
    print(
        f"{'condition':<40} {'scope':<9} {'objects':>9} {'mass share':>11} {'flux/true':>10} "
        f"{'S share':>9}"
    )
    print("-" * 124)
    for r in rows:
        print(
            f"{r['condition']:<40} {_scope(r['condition']):<9} {r['contig_objects']:>9,} "
            f"{r['contig_mass_share']:>11.3f} {_f(r['contig_agg'], 10, 3)} "
            f"{_f(r['spliced_share_of_contig'], 9, 3)}"
        )
    print(
        "⛔ `boundary_spliced` is a SEPARATE bank from `boundary_unspliced`, not a subset of it — the\n"
        "   same molecules split by whether they used a sj ELSEWHERE. So S SUBTRACTS from the estimated\n"
        "   RNA crossing rather than bounding f_g, and `S share` says how much of the contiguous RNA it\n"
        "   accounts for. ⚠ A first draft wrote `f_g <= 1 - S/M` and the truth violated it by 302."
    )


    print(f"\n{'=' * 124}")
    print("⑧ ⭐⭐⭐ THE RUNAWAY BOUND — the post-solve update iterated with the LIKELIHOOD REMOVED")
    print(f"{'=' * 124}")
    keys = list(rows[0]["runaway_fixed_points"])
    print(
        f"{'condition':<40} {'scope':<9} {'true rho_on':>12} "
        + " ".join(f"{k:>12}" for k in keys)
        + f" {'spread':>9}"
    )
    print("-" * 124)
    for r in rows:
        fp = r["runaway_fixed_points"]
        vals = [fp[k] for k in keys]
        finite = [v for v in vals if np.isfinite(v) and v > 0]
        spread = (
            max(finite) / min(finite)
            if len(finite) == len(vals) and min(finite) > 0
            else float("nan")
        )
        print(
            f"{r['condition']:<40} {_scope(r['condition']):<9} {_f(r['runaway_true_on'], 12, 6)} "
            + " ".join(_f(v, 12, 6) for v in vals)
            + f" {_f(spread, 9, 3)}"
        )
    print(
        "⭐ Four starts spanning three decades. `spread` = max/min of the fixed points reached: **1.000\n"
        "   means ONE attracting fixed point and the starting value does not matter**, which is exactly\n"
        "   the property a runaway lacks. ⛔ This is a BOUND, not the loop: the object's own likelihood\n"
        "   is removed, so the real refit loop has strictly LESS feedback than this. A pass here is\n"
        "   necessary and not sufficient — the full test needs the solve arm (stage 2)."
    )


# ── the falsification ────────────────────────────────────────────────────────────────────────────


class _Chain:
    def __init__(self, kind, obj_idx, left, right):
        self.kind, self.obj_idx = np.asarray(kind), np.asarray(obj_idx)
        self.left, self.right = np.asarray(left), np.asarray(right)
        self.n_slots = int(self.kind.shape[0])


class _Geom:
    def __init__(self, sj_count, eff_sj):
        self.sj_count, self.eff_sj = np.asarray(sj_count), np.asarray(eff_sj)


def _toy_chain(n_regions: int = 4):
    """A REGION/BOUNDARY chain of ``2*n_regions - 1`` slots with the real alternation and adjacency."""
    n = 2 * n_regions - 1
    kind = np.array([REGION if i % 2 == 0 else BOUNDARY for i in range(n)], np.int8)
    obj = np.array([i // 2 for i in range(n)], np.int64)
    left = np.array([i - 1 for i in range(n)], np.int64)
    right = np.array([i + 1 if i + 1 < n else -1 for i in range(n)], np.int64)
    return _Chain(kind, obj, left, right)


def self_test() -> int:
    """Every comparator perturbed, with no I/O: each block asserts the honest answer and that a
    deliberate corruption changes it."""
    passed = failed = 0

    def check(name: str, ok: bool):
        nonlocal passed, failed
        print(f"   {'✅' if ok else '⛔'} {name}")
        passed, failed = passed + bool(ok), failed + (not ok)

    print("\n── pooled_density: the ratio of sums, never the mean of ratios ──")
    mass = np.array([10.0, 1.0, 0.0])
    eff = np.array([100.0, 1.0, 5.0])
    sel = np.array([True, True, False])
    check("Σmass/Σeff on the selected pair", abs(pooled_density(mass, eff, sel) - 11.0 / 101.0) < 1e-12)
    check(
        "a mean OF RATIOS would read differently — the two are not the same number",
        abs(np.mean(mass[sel] / eff[sel]) - 11.0 / 101.0) > 0.1,
    )
    check("no opportunity ⇒ 0.0, never a floored division", pooled_density(mass, eff, ~sel) == 0.0)
    check(
        "excluding the dense member LOWERS it",
        pooled_density(mass, eff, np.array([False, True, False])) > pooled_density(mass, eff, sel),
    )

    print("\n── object_weighted_mean: objects with no mass are excluded ──")
    f = np.array([1.0, 1.0, 0.0, 0.0])
    live = np.array([True, True, False, False])
    check("the mean is over LIVE objects only", object_weighted_mean(f, live) == 1.0)
    check(
        "counting the empty ones would halve it — the exclusion is load-bearing",
        object_weighted_mean(f, np.ones(4, bool)) == 0.5,
    )
    check("no live object ⇒ NaN, never a fabricated 0", not np.isfinite(object_weighted_mean(f, np.zeros(4, bool))))

    print("\n── object_weighted_mean vs the FRAGMENT-weighted mean: they differ, and that is the finding ──")
    #: one huge RNA object and three tiny pure-gDNA ones: the panel's own shape in miniature
    n_g = np.array([1.0, 1.0, 1.0, 0.0])
    n_r = np.array([0.0, 0.0, 0.0, 997.0])
    m = n_g + n_r
    fg = n_g / m
    check("object-weighted = 0.75", abs(object_weighted_mean(fg, m > 0) - 0.75) < 1e-12)
    check("fragment-weighted = 0.003", abs(n_g.sum() / m.sum() - 0.003) < 1e-12)
    check(
        "and the gap is 250x — a prior applies PER OBJECT, so the denominators are not interchangeable",
        object_weighted_mean(fg, m > 0) / (n_g.sum() / m.sum()) > 100.0,
    )

    print("\n── neighbour_sj_density: the chain's own adjacency, pooled over slot + both neighbours ──")
    chain = _toy_chain(4)
    sj_count = np.zeros((chain.n_slots, 2))
    eff_sj = np.zeros((chain.n_slots, 2))
    sj_count[1, 0], eff_sj[1, 0] = 30.0, 10.0  # one sj on the first BOUNDARY
    rho, has = neighbour_sj_density(chain, _Geom(sj_count, eff_sj))
    check("the sj boundary reads its own flux density", abs(rho[1] - 3.0) < 1e-12)
    check("both chain neighbours inherit it", abs(rho[0] - 3.0) < 1e-12 and abs(rho[2] - 3.0) < 1e-12)
    check("a slot two hops away sees nothing", rho[3] == 0.0 and not has[3])
    check("coverage is exactly the reachable set", list(np.flatnonzero(has)) == [0, 1, 2])
    sj_count[1, 0] = 60.0
    rho2, _ = neighbour_sj_density(chain, _Geom(sj_count, eff_sj))
    check("doubling the flux doubles the density — the estimator is not inert", abs(rho2[1] - 6.0) < 1e-12)
    #: two sj on one boundary are two estimates of one rate, so the pooled statement is Σcount/ΣE
    sj_count[1, 1], eff_sj[1, 1] = 0.0, 10.0
    rho3, _ = neighbour_sj_density(chain, _Geom(sj_count, eff_sj))
    check("a second, silent sj on the same boundary HALVES it (Σcount/ΣE)", abs(rho3[1] - 3.0) < 1e-12)

    print("\n── strata: the boundary axis splits on whether MATURE RNA can cross ──")

    class _Statics:
        def __init__(self, fp, fn, mp=None, mn=None):
            self.free_pos, self.free_neg = fp, fn
            self.mrna_active_pos = np.zeros_like(fp) if mp is None else mp
            self.mrna_active_neg = np.zeros_like(fn) if mn is None else mn

    class _RA:
        def __init__(self, sig):
            self.signature = np.asarray(sig, np.uint8)

    from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS

    #: REGIONs: intergenic | intron | exon | exon, so the boundaries are
    #: intergenic|intron, intron|exon, exon|exon.
    sig = np.array([0, BIT_INTRON_POS, BIT_EXON_POS, BIT_EXON_POS], np.uint8)
    fp = np.array([False, False, True, True, True, True, True], bool)
    fn = np.zeros(7, bool)
    #: mature can cross only the exon|exon boundary (slot 5) and sit in the two exon REGIONs
    mp = np.array([False, False, False, False, True, True, True], bool)
    geom = _Geom(np.zeros((7, 2)), np.zeros((7, 2)))
    out = strata(chain, _Statics(fp, fn, mp, fn.copy()), geom, _RA(sig))
    lab = out["label"]
    check(
        "each REGION takes its signature's stratum",
        list(lab[[0, 2, 4, 6]]) == ["R intergenic", "R intron", "R exon", "R exon"],
    )
    check(
        "⭐ the boundary axis splits exon|intron from exon|exon and the gene edge",
        list(lab[[1, 3, 5]]) == ["B gene edge", "B exon|intron", "B exon|exon"],
    )
    check("every slot is labelled exactly once", sum(int(np.sum(lab == s)) for s in STRATA) == 7)
    check(
        "`on_target` is exon-touching on both axes",
        list(out["on_target"].astype(bool)) == [False, False, False, True, True, True, True],
    )
    #: the g1_locked cross-check must fire when an intergenic region is not structurally locked
    try:
        strata(chain, _Statics(fp | np.eye(7, dtype=bool)[0], fn, mp, fn.copy()), geom, _RA(sig))
        check("an intergenic REGION that is NOT g1_locked raises", False)
    except AssertionError:
        check("an intergenic REGION that is NOT g1_locked raises", True)
    #: the curation's own gate: if the solver says mature RNA may cross an `exon|intron` boundary,
    #: this file's near-pure-gDNA claim is about a different population, so it must raise.
    try:
        strata(chain, _Statics(fp, fn, mp | np.eye(7, dtype=bool)[3], fn.copy()), geom, _RA(sig))
        check("⭐ an `exon|intron` boundary the solver thinks mature can cross raises", False)
    except AssertionError:
        check("⭐ an `exon|intron` boundary the solver thinks mature can cross raises", True)
    #: and the partition gate must fire on a chain the labels do not cover
    try:
        strata(chain, _Statics(fp, fn, mp, fn.copy()), geom, _RA(np.array([0, 0, 0], np.uint8)))
        check("a chain whose objects the labels do not cover raises", False)
    except (AssertionError, IndexError):
        check("a chain whose objects the labels do not cover raises", True)

    print("\n── the gDNA-side estimator: exact at a zero control, wrong under a depleted anchor ──")
    eff_g = np.array([100.0, 100.0, 100.0, 100.0])
    est_mass = np.array([0.0, 0.0, 50.0, 50.0])  # anchors empty, exons full: the `g00` shape
    anchors = np.array([True, True, False, False])
    rho = pooled_density(est_mass, eff_g, anchors)
    est = np.clip(rho * eff_g / np.maximum(est_mass, _EPS), 0.0, 1.0)
    check("empty anchors ⇒ rho = 0 ⇒ every object reads f_g = 0", rho == 0.0 and est[2] == 0.0)
    #: the capture shape: the anchors are depleted 30x relative to the on-target field
    est_mass = np.array([1.0, 1.0, 30.0, 30.0])
    rho = pooled_density(est_mass, eff_g, anchors)
    est = np.clip(rho * eff_g / np.maximum(est_mass, _EPS), 0.0, 1.0)
    check(
        "a 30x-depleted anchor under-calls the on-target objects by 30x",
        abs(est[2] - 1.0 / 30.0) < 1e-9,
    )

    print(f"\n{'=' * 60}\n{passed}/{passed + failed} gates pass\n{'=' * 60}")
    return 0 if failed == 0 else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--oracle-cache", type=Path, default=None)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    suite = Path(args.suite)
    oracle_cache = args.oracle_cache or suite / "oracle_cache"
    index = TranscriptIndex.load(str(args.index))
    region_arrays = RegionArrays.from_index(index)
    sj = build_sj_geometry_arrays(index)
    boundary_flags = build_boundary_flags_array(index)
    names = args.conditions or sorted(p.name for p in Path(oracle_cache).iterdir() if p.is_dir())

    rows = []
    for name in names:
        row = measure_condition(index, region_arrays, sj, boundary_flags, suite, oracle_cache, name)
        rows.append(row)
        print(f"  {name} {row['seconds']:.0f}s", flush=True)
    report(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
