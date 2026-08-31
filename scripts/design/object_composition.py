"""IS THE TARGET `f_lib` OR THE OBJECT-WEIGHTED MEAN COMPOSITION — AND CAN IT BE ESTIMATED WITHOUT THE
SOLVE?

⭐⭐⭐ **This is the prior-free half of the Beta(a,b) reference question, and it runs with NO SOLVER, NO
EM and NO `src/` change.** psi's composition reference ``a*log f_g + b*log(1-f_g)`` is exactly
``Beta(a, b)`` in ``f_g``, so ``a`` and ``b`` are PSEUDO-COUNTS with a strength ``a+b`` and a MEAN
``a/(a+b)``. The shipped ``_JEFFREYS_REF = 0.5`` fixes the mean at 1/2 — an assertion that the library
is half gDNA. This file measures the two things that decide what the mean should be instead:

============  =====================================================================================
**TARGET**    what quantity the reference's mean wants: the FRAGMENT-weighted library composition
              ``f_lib``, or the OBJECT-weighted mean of the per-object truth. A prior applies once
              per object, so the second is the one with the matching denominator
**ESTIMATOR**  whether that quantity is reachable from the payload alone — no deconvolved array is
              read anywhere in this file, which is what makes an estimate here non-circular
============  =====================================================================================

⛔⛔ **THE HEADLINE, AND IT IS TWO NUMBERS THAT DISAGREE BY AN ORDER OF MAGNITUDE.** At ``g05
capture-OFF`` the fragment-weighted composition is **0.067** and the object-weighted one is **0.649**;
at ``g50`` they are **0.576** and **0.791**. gDNA is spread thinly over the whole genome while RNA is
concentrated in a few exons, so the many small intergenic and intronic objects are nearly pure gDNA
while the fragment mass sits in RNA-rich exons. ⭐ The direction is the one an independent sweep of
``vertex_ceiling.py`` already implied — that the optimal pseudo-count sits ABOVE ``f_lib`` — and at
``g50`` the object-weighted value lands on that sweep's optimum. ⛔ **Neither number ranks the work on
its own: `vertex_ceiling.py` does, and this file exists to tell it what value to install.**

⭐⭐⭐ **AND THE ANSWER TO BOTH IS THAT THE QUESTION WAS WRONG: THE REFERENCE DOES NOT HAVE TO BE
LIBRARY-WIDE AT ALL.** psi solves one object at a time, and the gDNA arm's fitted term is ALREADY per
slot — ``(n_slots, K)``. The reference is the only scalar left in psi. Written per object as

    m_i  =  rho_g,i * E_g,i  /  ( rho_g,i * E_g,i  +  rho_r,i * E_r,i )

every term but the two DENSITIES is per-object and exactly known, so even POOLED densities give a
per-object prior mean. Scored as ``sum |m_i - f_g,i| * M_i`` — how many fragments the prior misplaces if
believed outright — against the shipped constant 1/2, per stratum, in FRAGMENTS (table 6):

===========================  ====================  ====================  ===============
stratum                      shipped constant 1/2  ``m_i``, prior-free   class-pooled ref
===========================  ====================  ====================  ===============
stranded x capture OFF       1.000                 **0.040**             0.081
stranded x capture ON        1.000                 **0.326**             0.346
unstranded x capture OFF     1.000                 **0.040**             0.081
unstranded x capture ON      1.000                 0.326                 0.346  (DEFERRED)
``g00`` ZERO-gDNA control    1.000                 **0.000**             0.000
===========================  ====================  ====================  ===============

⭐⭐ **25x better than the shipped constant at capture-OFF, 3x at capture-ON, and EXACTLY ZERO at both
zero controls.** ⛔ It also beats the class-pooled reference arm on every stratum, so that arm is NOT a
ceiling: it hands an on-target RNA density to objects mature RNA cannot occupy, which the ``m_i`` form
refuses to do by construction. ⚠ ``g05``, which regressed 1.43x under every library-wide mean tried,
reads **0.009**.

⭐⭐⭐ **WHAT CARRIES THE PER-OBJECT VARIATION IS THE OPPORTUNITY GEOMETRY, NOT A PER-OBJECT DENSITY —
AND AT EXONS THE PRIOR IS FRANKLY A CLASS CONSTANT.** Measured at ``g50`` capture-OFF: within ``R exon``
``m_i`` has sd **0.0021** against a true ``f_g`` sd of **0.4441**, and within ``B exon|exon`` its sd is
**exactly 0**. Replacing it by its own class mean changes nothing (168,551 vs 164,074 fragments). ⛔ **Do
not read the win at the exonic strata as per-object resolution**; it is a better CONSTANT, and the
population that needs per-object resolution — exons, ``exon|exon``, AMBIG — is precisely the one the gDNA
LANDSCAPE exists to serve. The two terms partition the object universe rather than competing on it.

**THE LOCAL RNA DENSITY, and the answer is to shrink it hard (stage 0a).** Three ways of turning the sj
flux into ``rho_r``, per stratum at capture-OFF: raw local flux with a pooled fallback **0.121**, one
pooled scalar **0.042**, and one pseudo-observation of the pooled rate blended with the local flux
**0.040**. ⭐ The local flux carries a little signal and only survives heavy shrinkage; ⛔ and the naive
version that sets ``rho_r = 0`` where no sj is in reach reads **0.088-0.130 at the zero control**, where
the shrunk form reads **0.000**.

**THE SCOPE, and the answer is EVERYWHERE (stage 0b).** Letting the reference speak only where the
annotation determines the answer — the four strata where mature RNA cannot be — and leaving psi's 1/2
elsewhere reads **0.437 / 0.827**, i.e. 10x worse than letting it speak everywhere. At pass-0 there is no
landscape, so the exonic constant is strictly better than 1/2 there.

**THE POST-SOLVE UPDATE DOES NOT RUN AWAY (stage 0c, table 8).** The one update this design proposes —
re-estimating the on-target gDNA density from solved exons — is a single scalar measured on exons and
applied to exons, structurally the positive-feedback loop that makes a library-wide ``f_lib`` rule
inadmissible. Iterated with each object's OWN LIKELIHOOD REMOVED, which is a strict UPPER BOUND on the
feedback, four starts spanning three decades converge to the SAME fixed point to six decimals on all
twelve contaminated conditions and to exactly 0.000 at the control. ⭐ The per-object geometry damps it:
raising ``rho_g`` cannot raise the share at an object whose ``rho_r * E_r`` is large. ⛔ A pass here is
NECESSARY AND NOT SUFFICIENT, and the fixed point's VALUE is not the loop's value — the likelihood is
absent from it by construction.

⭐⭐ **THE STRATA, AND WHY THE BOUNDARY AXIS SPLITS ON MATURE-RNA CROSSABILITY.** Every slot lands in
exactly one of seven populations and the partition is asserted (:func:`strata`)::

    R intergenic     no transcript covers it        f_g = 1 STRUCTURALLY - no nascent assumption
    R intron         transcript, but no exon        f_g = 1 IF unspliced nascent RNA is sparse
    R exon           the mixture                    no structural claim
    B exon|intron    one flank exonic, in-gene      mature CANNOT cross - it splices. The ON-TARGET
                                                    gDNA anchor, which intergenic can never be
    B exon|exon      both flanks exonic             an alternative splice site inside an exonic
                                                    stretch. Mature crosses FREELY - not an anchor
    B intron|intron  neither flank exonic, in-gene  off-target, nascent-only
    B gene edge      a flank is intergenic          a TSS/TES interface; the ``g1_locked`` class

⛔⛔ **THAT CURATION IS THE CORRECTION THAT MADE THIS WORK** (owner, 2026-08-15). An earlier version used
"a sj attaches here" as one stratum and measured its true ``f_g`` at **0.0000 over 955,428 fragments** at
the zero control, because that pool lumps ``exon|exon`` in. ⭐ The same predicate then does a second job:
``rho_r`` from the sj flux is a MATURE density on the spliced template, so handing it to a slot mature
cannot occupy subtracts RNA that is not there. Gating it by the solver's own ``mrna_active`` takes
``R intron`` from 0.498 to **0.000** and ``B exon|intron`` from 0.485 to **0.000**.

⭐⭐ **HYBRID CAPTURE NEEDS NO DETECTION STEP** (table 5). Measure the gDNA density at BOTH anchors — the
off-target intergenic+intron REGIONs and the in-gene ``exon|intron`` boundaries — and their RATIO is the
enrichment: **0.98 without probes and 113-114 with them**, a 116x separation, no threshold and no flag.
The off-target anchor is accurate on all 16 (**0.91-1.00x** of truth). ⛔ The in-gene anchor is a clean
DETECTOR and not yet a calibrated level: it under-reads on-target gDNA by **2.6-3.6x** under capture,
because an ``exon|intron`` boundary sits at the EDGE of the probe footprint and its crossing opportunity
reaches into unprobed intron.

⛔ **``boundary_spliced`` IS A SEPARATE BANK FROM ``boundary_unspliced``, NOT A SUBSET OF IT** — the same
molecules split by whether they used a sj ELSEWHERE. So the contiguous RNA at a boundary is
``unspliced_RNA + S``, and S SUBTRACTS from the estimated RNA crossing rather than bounding ``f_g``.
⚠ A first draft of this file wrote ``f_g <= 1 - S/M`` as an assumption-free bound; it is false, and the
truth violated it by **302**. Table 7 keeps the identity under measurement instead.

⚠ **THE ONE ASSUMPTION THIS PANEL CANNOT PRICE.** ``nrna = 0`` on all sixteen conditions, so ``R intron``
and ``B exon|intron`` reading ``f_g = 1.0000`` is the panel's own build restated. ⭐ ``R intergenic`` is
the arm that does not depend on it — no transcript covers those regions at any nascent level — and
``bg_ig`` is printed beside ``bg`` for that reason. ⚠ On a STRANDED library the nascent half of an
``exon|intron`` boundary is deconvolvable by the strand channel anyway, so the assumption is load-bearing
on unstranded data only.

Gates: ``tests/test_object_composition_self_test.py``. ``--self-test`` perturbs every comparator with no
I/O, which is the other half of the discipline (TRAPS: perturb-every-gate).

Usage::

    python scripts/design/object_composition.py                    # the whole ladder
    python scripts/design/object_composition.py --conditions NAME
    python scripts/design/object_composition.py --self-test        # no I/O
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


PVO = _sibling("prior_vs_oracle.py")

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

#: Either strand's exon bit. ⭐ The union over strands is deliberate: a boundary is "exonic" on a flank
#: if ANY transcript has an exon there, and the per-STRAND question is answered by ``mrna_active``,
#: which :func:`strata` cross-checks against.
_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG

#: ⭐⭐ The SEVEN populations, in report order. ⛔ MUTUALLY EXCLUSIVE and EXHAUSTIVE over the chain,
#: asserted in :func:`strata` rather than promised here — a stratum table whose rows overlap
#: double-counts the object weight it exists to compute, and nothing downstream could tell.
#:
#: ⛔⛔ **THE BOUNDARY AXIS IS SPLIT BY WHETHER MATURE RNA CAN CROSS IT, NOT BY WHETHER A sj ATTACHES,
#: AND THAT DISTINCTION IS THE WHOLE POINT** (owner, 2026-08-15). An earlier version of this file used
#: "has a splice junction attached" as one stratum and measured its true ``f_g`` at **0.0000 over
#: 955,428 fragments** at the zero control — because that pool lumps ``exon|exon`` boundaries in, where
#: an alternative splice site sits inside a contiguous exonic stretch and mature RNA crosses freely.
STRATA = (
    "R intergenic",
    "R intron",
    "R exon",
    "B exon|intron",
    "B exon|exon",
    "B intron|intron",
    "B gene edge",
)

#: ⭐ The strata whose composition is knowable with NO deconvolution, and they are NOT equally free:
#: ``R intergenic`` is covered by no transcript at all, so ``f_g = 1`` holds at any nascent level;
#: ``R intron`` and ``B exon|intron`` need unspliced nascent RNA to be sparse. ⚠ On a STRANDED library
#: the nascent half of an ``B exon|intron`` boundary is deconvolvable by the strand channel anyway, so
#: the assumption is load-bearing on unstranded data only.
#:
#: ⛔⛔⛔ **MEASURED 2026-08-22: THIS POOL IS CONTAMINATED ON THE CURRENT PANEL, AND THE FACTOR IS NOT A
#: CONSTANT.** The panel used to hold ``nrna = 0``, which made ``R intron`` exactly pure and this pool
#: exactly right; it now carries SPARSE nascent RNA (``DESIGN.md`` §0b) and the pool's rate over its own
#: true gDNA rate reads **4.49x at g05, 1.18x at g50, 1.004x at g98 capture-OFF** — it scales with the
#: RNA:gDNA ratio, so it is worst exactly where gDNA is scarce and the anchor matters most.
#: ⭐ ``R intergenic`` ALONE is inflated **exactly 1.0000x** on every condition, and so is ``B gene
#: edge``. ⛔ So any number this instrument reports off ``PURE_GDNA_STRATA`` — including the recorded
#: exon-reference result (see `ISSUES: measured-prior-rung-4`'s anchor-contamination constraint) — was produced with an anchor up to 4.49x high and must be
#: re-derived before it is quoted. ⭐ The SHIPPED path is clean: ``density_deconv.fit_intron_background``
#: pools intergenic only (``include_introns=False`` at both call sites), so this is an INSTRUMENT defect
#: rather than shipped behaviour. Whether to drop ``R intron`` from this pool is a MEASUREMENT-DESIGN
#: decision (it is also the on-target-vs-off-target question), so it is recorded here rather than
#: silently changed.
PURE_GDNA_STRATA = ("R intergenic", "R intron")

#: ⭐⭐ The curated gDNA anchor that sits INSIDE genes, and therefore ON-TARGET under hybrid capture —
#: the one thing the intergenic anchor structurally cannot be.
ONTARGET_GDNA_STRATUM = "B exon|intron"


def strata(chain, statics, geometry, region_arrays) -> dict:
    """Per-slot stratum label and the masks the derivation needs, all from the ANNOTATION alone.

    ⛔ **The partition is asserted, not documented.** Every slot must take exactly one label: the
    object-weighted mean this file computes is a mean over objects, so a slot counted twice or not at
    all silently reweights the target it is trying to measure.

    **The boundary classification, and the reason it is by EXON-NESS of the two flanks.** A BOUNDARY's
    unspliced-crossing population is the molecules that crossed it CONTIGUOUSLY. Mature RNA can do that
    only where the template is contiguous exon on both sides, so:

    ================  =================================================================================
    ``B exon|intron``  exactly one flank exonic, both flanks inside a gene ⇒ **mature RNA CANNOT cross,
                       it splices** ⇒ near-pure gDNA under sparse unspliced nascent. ⭐ In-gene, so it
                       is the ON-TARGET gDNA anchor under capture
    ``B exon|exon``    both flanks exonic — an alternative splice site inside a contiguous exonic
                       stretch. ⛔ Mature RNA crosses it freely; **not an anchor and never was**
    ``B intron|intron`` neither flank exonic, both inside a gene (adjacent introns of different
                       signature). Off-target, nascent-only RNA
    ``B gene edge``    at least one flank intergenic — a TSS/TES interface. This is the ``g1_locked``
                       boundary class, structurally pure gDNA on both strands
    ================  =================================================================================

    ⭐ **``R intergenic`` is defined by the SIGNATURE and cross-checked against ``g1_locked``**, which is
    the predicate the solver itself pins on. Likewise ``B exon|intron`` is cross-checked against the
    solver's own ``mrna_active`` — if the two ever separate, this file's purity claim is about a
    different population than the one the solver reasons over, and it raises.
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
    # ⛔ THE CURATION, MADE FALSIFIABLE. `mrna_active_s` is the solver's own "contiguous exon on BOTH
    #    flanks" gate — exactly "mature RNA of strand s may cross here". At an `exon|intron` boundary one
    #    flank carries no exon bit at all, so it must be False on both strands. If that ever fails, the
    #    signature semantics moved and the near-pure-gDNA claim is about a different population.
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
        # ⭐ ON-TARGET means "touches an exon" — annotation-derived, no threshold and no capture
        #   detection. It is the axis hybrid capture enriches along, so it is the axis the gDNA density
        #   is allowed to differ across.
        "on_target": np.where(
            is_region,
            (np.where(is_region, sig[np.clip(obj, 0, max(n_regions - 1, 0))], 0) & _EXON_BITS) != 0,
            exon_l | exon_r,
        ),
        "has_sj": np.asarray(geometry.eff_sj, np.float64).sum(1) > 0.0,
    }


def slot_counts(payload, region_arrays, chain) -> np.ndarray:
    """One payload's unspliced/contained count per slot — **the mixture psi deconvolves, and nothing
    else.** ``region_contained`` at a REGION, ``boundary_unspliced`` at a BOUNDARY, exactly the
    populations :attr:`RegionGeometry.unspliced_count` carries, so a truth built from the origin
    partitions and an estimate built from the full payload are on one basis.
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
    """``sum(mass) / sum(eff)`` over a selected population — the ratio of sums, never the mean of
    ratios (TRAPS: a-mean-of-ratios-inherits-the-partition). ``0.0`` when the population has no
    opportunity, which is the honest answer and not a floored division."""
    e = float(np.sum(np.asarray(eff, np.float64)[select]))
    if e <= 0.0:
        return 0.0
    return float(np.sum(np.asarray(mass, np.float64)[select]) / e)


def neighbour_sj_density(chain, geometry) -> np.ndarray:
    """Per-slot certified-RNA density from the sj flux **at the slot and its two chain neighbours**,
    pooled as ``sum(count) / sum(E)``.

    ⭐ **The flux is a density on the SPLICED template**, so it is on the same footing as an exon
    REGION's contained RNA opportunity and NOT on the same footing as a BOUNDARY's unspliced-crossing
    opportunity — table ③ measures both and the difference is a factor of 7-10.

    ⚠ The neighbour set is the chain's own adjacency: genomic order IS slot order and the chain
    alternates REGION/BOUNDARY, so a slot's neighbours are ``i-1`` and ``i+1`` and a reference terminal
    links to ``-1``.
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
    """The TARGET: the mean per-object composition over objects that have any mass at all.

    ⛔ **Objects with no mass are excluded and that is not a convenience.** An empty object has no true
    composition to average — ``0/0`` — and folding a fabricated value in would move the target by the
    share of the genome that happens to be empty at this gDNA level, which at the zero control is most
    of it.
    """
    return float(np.mean(np.asarray(f_true, np.float64)[live])) if live.any() else float("nan")


def measure_condition(index, region_arrays, sj, boundary_flags, suite: Path, oracle_cache: Path,
                      condition: str) -> dict:
    """One condition: the target, the stratum census, and every estimator. One JSON-able row."""
    start = time.perf_counter()
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    kw = calibration_inputs(cache, index)
    chain = build_region_chain(
        cache.payload.ref_region_offsets, cache.payload.ref_boundary_offsets
    )
    statics = build_region_statics(chain, region_arrays, boundary_flags)
    geometry = build_region_geometry(
        chain,
        CalibrationSubstrate.from_payload(cache.payload, region_arrays),
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

    # ⭐ Through ``OracleTruth.from_parts`` rather than the raw payloads, so sum-to-full runs as a HARD
    #   gate on every condition — a cached partition that does not reconstruct the scan calibration read
    #   is a silently wrong truth source, and it would be invisible in every number below.
    # ⭐⭐ **The full payload is the SCAN CACHE'S, never the oracle cache's ``_main``, and that is a
    #   deliberate difference from ``calibration_vs_oracle.py``.** The two are the same scan — measured
    #   bit-identical on every integer bank and 3.5e-14 apart on the six float ones, which is the
    #   recorded re-association floor (`rescan_panels.py`) — so taking the one this file ALREADY loaded
    #   makes sum-to-full validate the truth partition against the exact array object the estimator
    #   reads, rather than against a second copy of it. ⚠ It also removes a read of a directory that
    #   ``pass0_vs_oracle.measure_condition`` WRITES whenever an arm runs with ``--oracle-cache``, so
    #   this instrument cannot race a concurrent arm.
    root = Path(oracle_cache) / condition
    parts = {k: read_scan_cache(root / k, index).payload for k in ORIGINS}
    OracleTruth.from_parts(cache.payload, parts)

    n_g = slot_counts(parts["gdna"], region_arrays, chain)
    n_r = slot_counts(parts["mrna"], region_arrays, chain) + slot_counts(
        parts["nrna"], region_arrays, chain
    )
    mass = n_g + n_r
    live = mass > 0.0
    f_true = np.where(live, n_g / np.maximum(mass, _EPS), 0.0)

    eff_g = np.asarray(geometry.eff_gdna, np.float64)
    eff_r = np.asarray(geometry.eff_rna, np.float64)

    # ⛔ The estimator reads the FULL payload's own per-slot totals, never the origin partitions. That
    #   is the non-circularity claim made structural: `truth` below is the only name bound to `parts`.
    est_mass = np.zeros(int(chain.n_slots), np.float64)
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    full_sub = CalibrationSubstrate.from_payload(cache.payload, region_arrays)
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

    # ── ⭐⭐⭐ THE PER-OBJECT PRIOR MEAN `m_i`, WHICH IS WHAT THE DERIVATION ACTUALLY NEEDS ──
    #
    #   m_i = rho_g,i * E_g,i / ( rho_g,i * E_g,i + rho_r,i * E_r,i )
    #
    # ⛔ **A LIBRARY-WIDE SCALAR IS THE SPECIAL CASE WHERE BOTH DENSITIES ARE CONSTANT AND THE GEOMETRY
    #   IS IGNORED.** Every term but the two densities is per-object and exactly known, so even a
    #   CLASS-POOLED density yields a per-object prior mean — which is the claim under test.
    # ⭐ Two gDNA densities, not one, split on an ANNOTATION-derived axis (does this object touch an
    #   exon?) rather than on a detected capture flag. Hybrid capture enriches along exactly that axis,
    #   so their RATIO is the enrichment factor: ~1 without probes and large with them, measured rather
    #   than declared (TRAPS: no-magic-numbers — there is no threshold and no switch here).
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

    # ⭐⭐⭐ **THE MATURE-RNA GATE, AND IT IS THE SAME PREDICATE THE OWNER'S BOUNDARY CURATION USES.**
    #   ``rho_r`` from the sj flux is a density of MATURE molecules on the spliced template. Handing it
    #   to a slot mature RNA cannot occupy — an intron REGION, an ``exon|intron`` boundary, a gene edge —
    #   subtracts RNA that is not there and calls gDNA RNA. ``mrna_active_s`` is the shipped predicate
    #   for exactly this on BOTH axes: a REGION's own exon bit, a BOUNDARY's contiguous exon on both
    #   flanks. ⚠ Where mature cannot be, what remains is NASCENT, which this design assumes sparse and
    #   this panel cannot price (``nrna = 0``).
    mature_here = np.asarray(statics.mrna_active_pos, bool) | np.asarray(
        statics.mrna_active_neg, bool
    )
    # the library-wide certified-RNA density on the spliced template: Sum(flux) / Sum(E_sj), the ratio of
    # sums and never a mean of ratios (TRAPS: a-mean-of-ratios-inherits-the-partition)
    _sjc = np.asarray(geometry.sj_count, np.float64).sum(1)
    _esj = np.asarray(geometry.eff_sj, np.float64).sum(1)
    rho_r_pooled = float(_sjc.sum() / max(float(_esj.sum()), _EPS))

    # ⛔⛔ **AND THE CERTIFIED SPLICED CROSSING IS A SUBTRACTION, NOT A BOUND.** ``boundary_spliced`` is a
    #   SEPARATE bank from ``boundary_unspliced`` — the same molecules split by whether they used a sj
    #   ELSEWHERE — so the contiguous RNA crossing a boundary is ``unspliced_RNA + S`` and
    #
    #       unspliced_RNA = rho_r * E_r - S      =>      f_g = 1 - (rho_r * E_r - S) / M
    #
    #   ⚠ A first draft of this file wrote ``f_g <= 1 - S/M`` as an assumption-free bound. It is simply
    #   FALSE — S is not inside M — and the truth violated it by up to 302 at ``g50``. ⚠ Subtracting S
    #   inside ``m_i`` was then measured a net LOSS (0.553 -> 0.570, over-correcting at ``exon|exon``),
    #   so it is NOT in the arm ladder; table (7) keeps the identity under measurement instead.
    spliced = np.zeros(int(chain.n_slots), np.float64)
    spliced[kind == BOUNDARY] = np.asarray(full_sub.boundary_spliced.count, np.float64).sum(1)[
        obj[kind == BOUNDARY]
    ]

    # ⭐⭐ THE LOCAL RNA DENSITY, THREE WAYS — stage 0a. The per-object sj flux is at the right LEVEL
    #   and far too NOISY per object (measured: 0.515 against 0.123 for a single pooled density). But
    #   the first version set ``rho_r = 0`` wherever no sj was in reach, which sends an RNA-rich exon to
    #   ``m_i = 1``; that is a coverage artefact, not a measurement. The honest fallback is the
    #   POPULATION rate, and the honest blend is one pseudo-observation of it — the same
    #   "one pseudo-object of ignorance" convention ``fit_landscape`` uses, so no constant is introduced.
    sj_num = np.asarray(geometry.sj_count, np.float64).sum(1)
    sj_den = np.asarray(geometry.eff_sj, np.float64).sum(1)
    num, den = sj_num.copy(), sj_den.copy()
    for side in (np.asarray(chain.left, np.int64), np.asarray(chain.right, np.int64)):
        ok = side >= 0
        num[ok] += sj_num[side[ok]]
        den[ok] += sj_den[side[ok]]
    #: the mean sj opportunity carried by ONE sj-bearing boundary — the weight of one pseudo-observation
    e_one = float(sj_den[sj_den > 0.0].mean()) if np.any(sj_den > 0.0) else 0.0
    rho_r_fallback = np.where(has_flux, rho_r, rho_r_pooled)
    rho_r_shrunk = (num + rho_r_pooled * e_one) / np.maximum(den + e_one, _EPS)

    def gated(x):
        """``rho_r`` is a MATURE density on the spliced template — zero where mature RNA cannot be."""
        return np.where(mature_here, x, 0.0)

    #: ⭐ the ARM LADDER, each an `m_i` and each scored the same way. `shipped` is the constant ½ ψ
    #: carries today — the baseline every ratio is against — and `TRUTH` is the class-pooled ceiling.
    #: ⛔ `structural` is stage 0b's scope question: let the reference speak ONLY where the annotation
    #: determines the answer (mature RNA cannot be here) and leave psi's ½ everywhere else, which is
    #: the population the gDNA LANDSCAPE exists to serve.
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
        # ⭐⭐⭐ THE DECONVOLUTION: `f_g = rho_g * E_g / M`, RNA as the RESIDUAL and never predicted.
        #   ⛔ This is not a new idea — it is the PEAK of the shipped `density_lambda_factor`, whose own
        #   docstring reads "peaked at f_g = rho_bg/rho_obs". With `rho_obs = M/E_g` that is exactly this
        #   expression, so the arm measures the LOCATION the shipped NegBinom factor already carries at
        #   the one stratum it is switched on for.
        #   ⭐ It needs NO RNA density, which is the whole asymmetry: gDNA is near-uniform and
        #   predictable pre-solve; RNA spans six decades with no genomic autocorrelation and is never
        #   predicted — it is whatever mass the gDNA deconvolve leaves behind.
        "deconvolve": np.clip(rho_g_per_object * eff_g / np.maximum(est_mass, _EPS), 0.0, 1.0),
        # ⭐⭐⭐ THE HYBRID: pin where the annotation DETERMINES the answer, deconvolve where it does not.
        #   The structural strata are exactly ``~mature_here`` — mature RNA cannot be there — and the
        #   deconvolve's use of the observed ``M`` makes it strictly worse on them (a downward Poisson
        #   fluctuation reads as RNA). Everywhere else the annotation is silent and the deconvolution is the only
        #   honest statement.
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

    # ── ⭐⭐⭐ STAGE 0c — THE RUNAWAY BOUND, WITH THE LIKELIHOOD REMOVED ──
    #
    # The only post-solve update this design proposes is narrow: re-estimate the ON-TARGET gDNA density
    # from solved exons, because the pre-solve ``exon|intron`` anchor under-reads it 2.6-3.6x under
    # capture. ⛔ That is ONE SCALAR, measured on exons and applied to exons — structurally the same
    # positive-feedback loop that makes a library-wide ``f_lib`` rule inadmissible: a higher rho_g raises
    # every exonic ``m_i``, which raises the gDNA mass attributed to exons, which raises rho_g.
    #
    # ⭐⭐ **Iterating that map with each object's OWN LIKELIHOOD REMOVED is a strict UPPER BOUND on the
    # feedback** — it is the prior believed outright, with no data pulling against it. If this converges,
    # the real refit loop converges a fortiori; if it runs away, the hazard is real and the full test is
    # owed. That is what makes a solver-free answer admissible here at all.
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
    #: ⭐ four starts spanning three decades around the measured anchor: if they meet, the map has ONE
    #: attracting fixed point and the starting value does not matter (TRAPS: perturb-every-gate).
    runaway = {
        k: _runaway(v)
        for k, v in (
            ("anchor", rho_g_on),
            ("10x low", rho_g_on / 10.0),
            ("10x high", rho_g_on * 10.0),
            ("truth", true_on),
        )
    }

    # ── THE CONTIGUOUS-RNA IDENTITY, kept under measurement ──
    # ``rho_r * E_r = unspliced_RNA + S`` at a BOUNDARY. Scored as: does the neighbouring sj flux
    # recover the TRUE total contiguous RNA density ``(n_r + S)/E_r``? ⭐ That is the quantity the
    # corrected estimator needs, and it is better posed than the ``n_r/E_r`` this file asked for first.
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
        # ⚠ the same mean over the objects psi actually SOLVES a composition on: `g1_locked` slots are
        #   pinned, so the reference never moves them and their share of the target is inert.
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

    # ── ③ the sj flux as an RNA density, on TWO bases that differ by a factor of 7-10 ──
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

#: ⭐⭐ The 0.8.0 scope, stamped on every row rather than left to the reader — three strata are the
#: development target and unstranded x capture-ON is DEFERRED-but-REPORTED. Same table as
#: ``calibration_vs_oracle.py``'s; a reader ranking on a stratum that is not a target inverts the order.
_SCOPE = {
    ("stranded", "capture OFF"): "IN SCOPE",
    ("stranded", "capture ON"): "IN SCOPE",
    ("unstranded", "capture OFF"): "IN SCOPE",
    ("unstranded", "capture ON"): "DEFERRED",
}


def _scope(condition: str) -> str:
    return "CONTROL" if PVO.is_zero_gdna(condition) else _SCOPE[PVO.stratum(condition)]


#: Every selection the per-stratum tables print, in order — one list, so a stratum cannot appear on some
#: tables and not others. ⛔ The zero control is its own row and is never folded into a stratum: its truth
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
    """⛔ **Every comparator perturbed, with no I/O.** A gate that cannot fail is not a gate — each
    block below asserts the honest answer AND that a deliberate corruption changes it."""
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
    #: one huge RNA object and three tiny pure-gDNA ones — the panel's own shape in miniature
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
    #: ⭐ two sj on one boundary are two estimates of ONE rate, so the pooled statement is Σcount/ΣE
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

    #: REGIONs: intergenic | intron | exon | exon — so the boundaries are
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
    #: ⛔ the g1_locked cross-check must FIRE when an intergenic region is not structurally locked
    try:
        strata(chain, _Statics(fp | np.eye(7, dtype=bool)[0], fn, mp, fn.copy()), geom, _RA(sig))
        check("an intergenic REGION that is NOT g1_locked raises", False)
    except AssertionError:
        check("an intergenic REGION that is NOT g1_locked raises", True)
    #: ⛔⛔ AND THE CURATION'S OWN GATE: if the solver says mature RNA may cross an `exon|intron`
    #:    boundary, this file's near-pure-gDNA claim is about a different population. Must raise.
    try:
        strata(chain, _Statics(fp, fn, mp | np.eye(7, dtype=bool)[3], fn.copy()), geom, _RA(sig))
        check("⭐ an `exon|intron` boundary the solver thinks mature can cross raises", False)
    except AssertionError:
        check("⭐ an `exon|intron` boundary the solver thinks mature can cross raises", True)
    #: ⛔ and the partition gate must FIRE on a chain the labels do not cover
    try:
        strata(chain, _Statics(fp, fn, mp, fn.copy()), geom, _RA(np.array([0, 0, 0], np.uint8)))
        check("a chain whose objects the labels do not cover raises", False)
    except (AssertionError, IndexError):
        check("a chain whose objects the labels do not cover raises", True)

    print("\n── the gDNA-side estimator: exact at a zero control, wrong under a depleted anchor ──")
    eff_g = np.array([100.0, 100.0, 100.0, 100.0])
    est_mass = np.array([0.0, 0.0, 50.0, 50.0])  # anchors EMPTY, exons full — the `g00` shape
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
