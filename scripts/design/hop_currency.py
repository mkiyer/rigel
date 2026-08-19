"""THE HOP CURRENCY MAP — which currency, LEVEL or COMPOSITION, each hop type carries, MEASURED.

⭐⭐⭐ **STAGE 2 OF THE MESSAGE-PROPAGATION REBUILD: a MEASUREMENT, no solver runs, no `src/` change.**
A hop between two adjacent chain slots carries either a gDNA **LEVEL** (``rho_g``, fragments per
placement — invariant to a POPULATION change, destroyed by an ENRICHMENT change) or a gDNA
**COMPOSITION** (``f_g``, the gDNA share of the crossing population — the reverse). The two invariances
are complementary, so at every hop at least one survives; this instrument lets the certified truth say
WHICH, per hop type, before a line of policy code is written. Three questions, on every condition:

**① THE HOP CENSUS** — every ordered adjacent pair ``dst <- src`` classified by the two slots' OBJECT
CLASSES (`object_composition.strata`'s seven labels — ``R intergenic`` … ``B gene edge``, partition-
asserted — which that module calls "strata"; here they are CLASSES, and "stratum" is reserved for the
panel's strand x capture axis) with hop count and DESTINATION MASS per type; plus the DEPTH distribution — hops from each ``mrna_active`` slot to the nearest slot no mature
transcript crosses, mass-weighted — which says whether a chain is needed at all.

**② THE CURRENCY ORACLE** — for every hop, the source's TRUE value is transported BOTH ways and scored
at the destination in FRAGMENTS, ``sum |f_hat - true_f_g| * M``::

    LEVEL:        f_hat(dst) = clip( rho_g_true(src) * E_g(dst) / M(dst), 0, 1 )
    COMPOSITION:  f_hat(dst) = f_g_true(src population that ENTERS dst)

The source is perfect, so what remains is the RULE's error, and the smaller column IS that hop type's
currency. ⭐ **The population of a message is DIRECTION-DEPENDENT** — BOUNDARY -> REGION carries what
physically crosses INTO the region: the unspliced crossings PLUS the contiguous spliced crossings PLUS
the sj flux whose body lies on the destination's side (``sj_count_lo`` enters the LEFT neighbour,
``sj_count_hi`` the RIGHT); REGION -> BOUNDARY carries the region's contained population, which is
unspliced by construction. ⛔ Using a boundary's raw unspliced composition into an exon is the
prototype defect this file replaces.

⭐ **Every error is printed beside its NOISE FLOOR** — the same statistic under the rule's OWN null
(equal density, or equal composition) with the source and destination re-drawn as Poisson / binomial
counts. A hop type whose error sits at its floor has a CORRECT rule and finite counts; one far above
it has a wrong rule. Without the floor, ~1e2 fragments per hop reads as a several-percent "residual".

**③ THE RESIDUAL CENSUS** — where ``min(LEVEL, COMPOSITION)`` is still far above its floor, neither
currency works and a third MECHANISM is needed; those hop types are NAMED before any is proposed.

⛔ Scored on the certified per-slot truth (`calibration_oracle.py`'s ``slot_truth.npz``; REFUSED if
absent) plus the origin-split spliced/sj banks, with the ``gdna-cannot-splice`` gate asserted rather
than assumed. Reported per condition and never pooled across panel strata; ``g00`` rows are stamped — both
currencies are exact there by construction, so the zero control is uninformative about currency.

⭐⭐⭐ **THE VERDICT — 16 conditions, 2026-08-19, the rebuilt panel, 64/64 gates green.** Three findings
before the table, because the table is only legible with them:

* ⛔⛔ **THE ARMS DISAGREE, WHICH IS WHY THE MAP HAD TO BE PER ARM** (owner ruling). At
  `R exon <- B exon|exon[term]` off capture the gDNA arm wants a **LEVEL** (0.0 % excess over its
  floor) while both RNA arms want a **COMPOSITION** (5.2 % / 4.1 %) — the pooled first pass called that
  hop "LEVEL" outright and hid it. Same shape at `B exon|exon[term] <- R exon` under capture.
* ⛔⛔ **THE SEVEN OBJECT CLASSES DO NOT NAME THE HOP TYPES. A BOUNDARY'S STRUCTURAL FLAGS DO.** Split by
  the splice graph's flags, the composition error into exons sits on the **TERMINUS** boundaries (a
  TSS/TES lying inside another transcript's intron) while the **SPLICE-SITE** ones are near-exact. A
  `B exon|intron` is two hop types (15,480 `[sj]` vs 3,867 `[term]`), and so is a `B exon|exon`
  (4,838 / 7,850). The key is ``object class x {sj, term, sj+term}``.
* ⛔ **EVERY ERROR NEEDS ITS NOISE FLOOR, AND THE FLOOR MUST REPRODUCE THE SCORER'S OWN SELECTION.** A
  rule that is exactly right still scores the sampling noise of a ~1e2-fragment source; and scoring only
  non-empty destinations zero-truncates an object expecting ~1.2 crossings, which reads as a 10 % "rule
  error" that is not there.

**The measured map — modal currency across the contaminated conditions, median excess over floor:**

======================================  =======  ==============================================
hop (dst <- src)                        arms     currency
======================================  =======  ==============================================
exon <- SPLICE-SITE boundary            all 3    **COMPOSITION** (the SPLICE-IN population).
 (`exon|exon[sj]`, `exon|intron[sj]`)            Under capture 0.4-1.1 % on every arm; at floor off it
exon <- TERMINUS boundary               gDNA     **LEVEL** (0.0-0.7 %) ...
 (`exon|exon[term]`)                    RNA      ... but **COMPOSITION** for both RNA arms (4-8 %)
boundary <- EXON (SPLICE OUT)           gDNA     **LEVEL**, at its floor off capture, 0.8-0.9 % on
 (`exon|exon[term]`, `exon|exon[sj]`)   RNA      **COMPOSITION** under capture (1.4-2.6 %)
`exon|intron`/gene edge <- its own      all 3    **LEVEL** at its floor off capture; the intron's own
 intron / intergenic                             composition where the populations match
⛔ exon <- `B gene edge[term]`,         all 3    **NEITHER** — ~10 % excess on every arm, the residual
 UNDER CAPTURE                                   census, and `ROADMAP.md` §1 rank 3's spike-and-slab
======================================  =======  ==============================================

⭐ **Depth**: capture-OFF puts 23-29 % of imputed mass >= 9 hops from any measured gDNA and depth 1
reaches only ~31 %; capture-ON inverts the measured/imputed ratio with ~40 % at depth 1. The chain stays.

Usage::

    python scripts/design/hop_currency.py                      # all 32 conditions + summary
    python scripts/design/hop_currency.py --condition NAME     # one condition, full detail
    python scripts/design/hop_currency.py --out hops.tsv       # every (condition x hop type) row
    python scripts/design/hop_currency.py --self-test          # perturb every comparator, no I/O
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
from scipy.stats import poisson  # noqa: E402


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


OC = _sibling("object_composition.py")
PVO = OC.PVO

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import build_region_geometry, build_region_statics  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    build_boundary_flags_array,
    build_sj_geometry_arrays,
    is_splice_site,
    is_terminus,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402
from rigel.types import Strand  # noqa: E402

DEFAULT_SUITE = OC.DEFAULT_SUITE
DEFAULT_INDEX = OC.DEFAULT_INDEX
LEVEL, COMP = "LEVEL", "COMP"
#: ⭐⭐⭐ **THE THREE ARMS, NATIVELY — AXIOM 0, and the owner ruled the map must carry them** (2026-08-19).
#: Order is the component axis of every per-slot array here. ⛔ A map that pools the two RNA arms has
#: not measured what the policy carries: the message is three DENSITIES (`EQUATIONS.md` §3.5e).
COMPONENTS = ("gDNA", "RNA+", "RNA-")
#: which origin partition supplies each component's TRUE counts. ⭐ `rna_pos`/`rna_neg` are the RNA reads
#: split by TRANSCRIPT strand (`calibration/_oracle.RNA_STRAND_ORIGINS`), so a component's spliced and sj
#: banks are restricted to the right strand BY THE PARTITION — no bank-column convention is involved,
#: which matters because those banks use two different strand keyings.
COMPONENT_ORIGIN = ("gdna", "rna_pos", "rna_neg")
#: Monte-Carlo replicates for the noise floor. The floor is a DIAGNOSTIC column, and its own sampling
#: error at R replicates over ~1e4 hops is ~1/sqrt(R * 1e4) relative — far below anything read off
#: it; the replicate spread is printed beside it so that claim is checkable rather than asserted.
_FLOOR_REPS = 8
#: depth buckets for the census: exact to 4, then 5-8, then 9+ (the pre-fix finding was "~28 % >= 9").
_DEPTH_EDGES = (1, 2, 3, 4, 8)


# ── ① the census: hops and depth ─────────────────────────────────────────────────────────────────


def enumerate_hops(left: np.ndarray, right: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Every ORDERED adjacent pair ``(dst, src)``, each exactly once, from the chain's own adjacency.

    ⭐ A hop is a message from ``src`` into ``dst``; both directions of every adjacency are hops, so
    ``n_hops == 2 * n_adjacencies``. A reference terminal links to ``-1`` and is not a hop.
    """
    left = np.asarray(left, np.int64)
    right = np.asarray(right, np.int64)
    i = np.arange(left.shape[0], dtype=np.int64)
    fwd = right >= 0  # src i -> dst right[i]
    bwd = left >= 0  # src i -> dst left[i]
    dst = np.concatenate([right[fwd], left[bwd]])
    src = np.concatenate([i[fwd], i[bwd]])
    if dst.shape[0] != 2 * int(fwd.sum()):
        raise AssertionError("adjacency is not symmetric: left/right disagree on the number of pairs")
    return dst, src


def hop_type(label: np.ndarray, dst: np.ndarray, src: np.ndarray) -> np.ndarray:
    """``"<dst class> <- <src class>"`` per hop."""
    lab = np.asarray(label).astype(str)
    return np.char.add(np.char.add(lab[dst], " <- "), lab[src])


def boundary_structure(kind: np.ndarray, flags: np.ndarray) -> np.ndarray:
    """Per slot, the BOUNDARY's structural class from the splice graph's flags, either strand:
    ``"sj"`` (a donor/acceptor only), ``"term"`` (a TSS/TES only), ``"sj+term"`` (both), ``""`` on a
    REGION. ⛔ A boundary carrying neither is a graph defect and raises.

    ⭐ **This split is what the seven object classes cannot see and the currency oracle needs** (measured
    2026-08-19): a ``B exon|intron`` is a splice site — RNA ENTERS the exon across it, by the sj —
    or a transcript's TSS/TES lying inside another transcript's intron — RNA ORIGINATES there and
    crosses nothing. The first is a composition-preserving hop; the second is a POPULATION change.
    """
    kind = np.asarray(kind)
    flags = np.asarray(flags, np.uint16)
    term = is_terminus(flags, Strand.POS) | is_terminus(flags, Strand.NEG)
    spl = is_splice_site(flags, Strand.POS) | is_splice_site(flags, Strand.NEG)
    is_b = kind == BOUNDARY
    if np.any(is_b & ~term & ~spl):
        raise AssertionError("a BOUNDARY carries neither a terminus nor a splice-site flag")
    out = np.full(kind.shape[0], "", dtype=object)
    out[is_b & spl & ~term] = "sj"
    out[is_b & term & ~spl] = "term"
    out[is_b & term & spl] = "sj+term"
    return out.astype(str)


def refine_labels(label: np.ndarray, structure: np.ndarray) -> np.ndarray:
    """``"B exon|intron"`` -> ``"B exon|intron[sj]"`` etc.; REGION labels unchanged."""
    lab = np.asarray(label).astype(str)
    st = np.asarray(structure).astype(str)
    return np.where(st != "", np.char.add(np.char.add(np.char.add(lab, "["), st), "]"), lab)


def depth_to_measured(left: np.ndarray, right: np.ndarray, ref: np.ndarray, ma: np.ndarray) -> np.ndarray:
    """Per slot, the number of hops to the nearest ``~ma`` slot in the same reference (0 on a ``~ma``
    slot itself, ``inf`` when the whole reference is ``ma``). Computed in both directions and the
    minimum taken — the chain is a line, so the nearest measured gDNA is on one side or the other.
    """
    n = int(np.asarray(left).shape[0])
    ma = np.asarray(ma, bool)
    ref = np.asarray(ref)
    idx = np.flatnonzero(~ma)
    depth = np.full(n, np.inf)
    depth[~ma] = 0.0
    if idx.size == 0:
        return depth
    pos = np.arange(n)
    # nearest measured slot at or below, and at or above, by position; a ref is a contiguous block
    # of slot ids, so "same ref" is a check, not a search
    k = np.searchsorted(idx, pos, side="right") - 1
    lo_ok = k >= 0
    lo = np.where(lo_ok, idx[np.clip(k, 0, idx.size - 1)], -1)
    lo_ok &= (lo >= 0) & (ref[np.clip(lo, 0, n - 1)] == ref)
    k2 = np.searchsorted(idx, pos, side="left")
    hi_ok = k2 < idx.size
    hi = np.where(hi_ok, idx[np.clip(k2, 0, idx.size - 1)], -1)
    hi_ok &= (hi >= 0) & (ref[np.clip(hi, 0, n - 1)] == ref)
    d_lo = np.where(lo_ok, pos - lo, np.inf)
    d_hi = np.where(hi_ok, hi - pos, np.inf)
    d = np.minimum(d_lo, d_hi)
    depth[ma] = d[ma]
    return depth


def depth_census(depth: np.ndarray, mass: np.ndarray, ma: np.ndarray) -> dict:
    """Mass-weighted depth distribution over the ``ma`` (imputed) slots, bucketed by ``_DEPTH_EDGES``."""
    mass = np.asarray(mass, np.float64)
    ma = np.asarray(ma, bool)
    imputed = float(mass[ma].sum())
    measured = float(mass[~ma].sum())
    d = depth[ma]
    m = mass[ma]
    buckets = {}
    lo = 1
    for hi in _DEPTH_EDGES:
        sel = (d >= lo) & (d <= hi)
        buckets[f"{lo}" if lo == hi else f"{lo}-{hi}"] = float(m[sel].sum())
        lo = hi + 1
    buckets[f">={lo}"] = float(m[(d >= lo) & np.isfinite(d)].sum())
    buckets["unreachable"] = float(m[~np.isfinite(d)].sum())
    return {"measured_mass": measured, "imputed_mass": imputed, "buckets": buckets}


# ── ② the currency oracle ────────────────────────────────────────────────────────────────────────


def source_population(kind, left, right, dst, src, M, spliced, sj_lo, sj_hi):
    """The COMPOSITION source population per hop, DIRECTION-DEPENDENT — what physically enters ``dst``.

    ========================  =====================================================================
    src is a BOUNDARY         unspliced crossings + contiguous spliced crossings + the sj flux whose
    (dst a REGION)            body lies on ``dst``'s side: ``sj_lo`` if ``dst`` is the boundary's
                              LEFT neighbour, ``sj_hi`` if its RIGHT (`build_region_geometry`
                              attaches a sj's genomic-LOW end to the boundary right of its source
                              region, so ``_lo`` flux has its body in the LEFT flank)
    src is a REGION           its contained population — unspliced by construction, a contained
    (dst a BOUNDARY)          fragment used no sj
    ========================  =====================================================================

    ⭐ Every input carries a trailing COMPONENT axis (:data:`COMPONENTS`), so the population is
    returned per arm: ``pop[hop, c]``. gDNA's spliced and sj banks are identically zero (gated), so its
    entering population is its unspliced crossings whatever the direction — which is the arithmetic
    statement of "gDNA cannot splice".
    """
    kind = np.asarray(kind)
    src_is_b = (kind[src] == BOUNDARY)[:, None]
    dst_is_left = (np.asarray(left)[src] == dst)[:, None]
    side = np.where(dst_is_left, sj_lo[src], sj_hi[src])
    return np.where(src_is_b, M[src] + spliced[src] + side, M[src])


def transport(dst, src, M, E, pop) -> dict:
    """Both currencies transported from a PERFECT source and scored at the destination, PER ARM.

    ``M[slot, c]`` is the component's true unspliced count, ``E[slot, c]`` its opportunity
    (``eff_gdna`` for gDNA, ``eff_rna`` for both RNA arms), ``pop[hop, c]`` the population entering the
    destination. Per hop and per component::

        LEVEL:        count_hat = min( M_c(src)/E_c(src) * E_c(dst),  M_total(dst) )
        COMPOSITION:  count_hat = pop_c(src)/pop_total(src) * M_total(dst)

    and the error is ``|count_hat − M_c(dst)|`` in FRAGMENTS. ⛔ A LEVEL claim is CLIPPED at the
    destination's observed TOTAL: a component cannot hold more fragments than the object has, so the
    excess is certainly wrong and is reported separately as ``clipped`` rather than inflating the score
    without bound.
    """
    tot_src = pop.sum(axis=1)
    M_tot_dst = M[dst].sum(axis=1)
    E_s, E_d = E[src], E[dst]
    level_ok = E_s > 0.0
    comp_ok = (tot_src > 0.0)[:, None]
    rho = np.where(level_ok, M[src] / np.where(level_ok, E_s, 1.0), 0.0)
    claim = rho * E_d
    clipped = np.maximum(claim - M_tot_dst[:, None], 0.0)
    count_level = np.minimum(claim, M_tot_dst[:, None])
    f_comp = np.where(comp_ok, pop / np.where(comp_ok, tot_src[:, None], 1.0), 0.0)
    count_comp = f_comp * M_tot_dst[:, None]
    # ⛔ a hop is scored only where the destination HOLDS mass and the source can supply BOTH claims;
    #    scoring one currency where the other is undefined would compare them on different populations
    scorable = (M_tot_dst > 0.0) & level_ok.all(axis=1) & (tot_src > 0.0)
    return {"err_level": np.abs(count_level - M[dst]), "err_comp": np.abs(count_comp - M[dst]),
            "scorable": scorable, "clipped": clipped,
            "count_level": count_level, "count_comp": count_comp}


def noise_floor(dst, src, M, E, pop, rng, reps: int = _FLOOR_REPS, truncate: bool = True) -> dict:
    """The per-hop, per-ARM error a PERFECT rule would show from sampling alone, under each currency's
    own null.

    LEVEL null: one density per component, ``rho_c = (M_c(src) + M_c(dst)) / (E_c(src) + E_c(dst))``;
    redraw both counts Poisson and rescore. COMPOSITION null: one share per component,
    ``f_c = (pop_c + M_c(dst)) / (pop_tot + M_tot(dst))``; redraw both multinomially and rescore.

    ⛔ **The null must reproduce the SELECTION, or the floor is biased low.** A hop is scored only if
    its destination holds mass, so an object whose ONLY population is one component — every pure-gDNA
    boundary at ``g05``, expecting ~1 crossing — is a ZERO-TRUNCATED Poisson (mean 1.6 where the rate
    is 1.2), and an unconditional redraw reads that truncation as a 10 % rule error (measured
    2026-08-19). Such destinations are redrawn conditional on ``>= 1`` by the inverse CDF.
    """
    n_hops, n_comp = pop.shape
    M_dst, M_src = M[dst], M[src]
    M_tot_dst = M_dst.sum(axis=1)
    E_s, E_d = E[src], E[dst]
    lv = E_s > 0.0
    rho = np.where(lv, (M_src + M_dst) / np.where(lv, E_s + E_d, 1.0), 0.0)
    lam_d = rho * E_d
    tot_src = pop.sum(axis=1)
    denom = tot_src + M_tot_dst
    f = np.where((denom > 0)[:, None], (pop + M_dst) / np.where((denom > 0)[:, None], denom[:, None], 1.0), 0.0)
    pop_i = np.rint(tot_src).astype(np.int64)
    M_i = np.rint(M_tot_dst).astype(np.int64)
    # a destination whose whole population is ONE component was selected by that component being
    # non-zero, so its redraw is zero-truncated (the other components contribute no floor of their own)
    lone = (M_dst > 0).sum(axis=1) == 1
    fl = np.zeros((n_hops, n_comp))
    fc = np.zeros((n_hops, n_comp))
    reps_l, reps_c = [], []
    for _ in range(reps):
        ns = rng.poisson(rho * E_s)
        nd = rng.poisson(lam_d).astype(np.float64)
        if truncate and lone.any():
            for c in range(n_comp):
                sel = lone & (M_dst[:, c] > 0)
                if not sel.any():
                    continue
                lt = lam_d[sel, c]
                u = np.exp(-lt) + rng.random(lt.shape[0]) * (1.0 - np.exp(-lt))
                nd[sel, c] = np.maximum(poisson.ppf(u, lt), 1.0)
        Md_tot = nd.sum(axis=1)
        claim = np.where(lv, ns / np.where(lv, E_s, 1.0), 0.0) * E_d
        el = np.abs(np.minimum(claim, Md_tot[:, None]) - nd)
        bs = rng.binomial(np.repeat(pop_i[:, None], n_comp, axis=1), f)
        bd = rng.binomial(np.repeat(M_i[:, None], n_comp, axis=1), f)
        share = np.where((bs.sum(axis=1) > 0)[:, None],
                         bs / np.maximum(bs.sum(axis=1), 1)[:, None], 0.0)
        ec = np.abs(share * M_tot_dst[:, None] - bd)
        fl += el
        fc += ec
        reps_l.append(el)
        reps_c.append(ec)
    return {"floor_level": fl / reps, "floor_comp": fc / reps,
            "reps_level": np.stack(reps_l), "reps_comp": np.stack(reps_c)}


def per_type(htype: np.ndarray, scorable: np.ndarray, M_dst_tot: np.ndarray, tr: dict,
             fl: dict | None) -> list[dict]:
    """Aggregate per ``(hop type, COMPONENT)``, in FRAGMENTS.

    ⭐⭐ **A CURRENCY PER ARM, NOT PER HOP TYPE** (owner, 2026-08-19). A hop may perfectly well carry a
    gDNA LEVEL and an RNA COMPOSITION — the two have complementary invariances (`DESIGN.md` §0c.0) and
    nothing makes the three arms agree — so collapsing them to one winner per hop type would hide
    exactly the structure the policy's table needs.

    ``excess`` is the currency's error minus its own noise floor, as a share of the scored destination
    mass, so it is comparable across arms of very different size.
    """
    rows = []
    for t in sorted(set(htype.tolist())):
        sel = htype == t
        sc = sel & scorable
        for ci, comp in enumerate(COMPONENTS):
            row = {
                "hop_type": t, "component": comp,
                "n_hops": int(sel.sum()), "dst_mass": float(M_dst_tot[sel].sum()),
                "n_scored": int(sc.sum()), "mass_scored": float(M_dst_tot[sc].sum()),
                "err_level": float(tr["err_level"][sc, ci].sum()),
                "err_comp": float(tr["err_comp"][sc, ci].sum()),
                "clipped": float(tr["clipped"][sc, ci].sum()),
                "floor_level": float(fl["floor_level"][sc, ci].sum()) if fl else float("nan"),
                "floor_comp": float(fl["floor_comp"][sc, ci].sum()) if fl else float("nan"),
                "floor_level_sd": float(fl["reps_level"][:, sc, ci].sum(1).std()) if fl else float("nan"),
                "floor_comp_sd": float(fl["reps_comp"][:, sc, ci].sum(1).std()) if fl else float("nan"),
            }
            el, ec = row["err_level"], row["err_comp"]
            if row["mass_scored"] <= 0.0:
                row["winner"] = "—"
                row["excess"] = row["excess_level"] = row["excess_comp"] = float("nan")
                rows.append(row)
                continue
            ms = row["mass_scored"]
            row["excess_level"] = max(el - (row["floor_level"] if fl else 0.0), 0.0) / ms
            row["excess_comp"] = max(ec - (row["floor_comp"] if fl else 0.0), 0.0) / ms
            xl, xc = row["excess_level"], row["excess_comp"]
            # the currency is the one whose error is explained by sampling alone. Each floor replicate
            # is one realisation of the total under the null, so "within 2 sd of the floor" is the
            # consistency test; when BOTH pass, either serves and the cell says so rather than picking
            # on noise; when neither does, the smaller excess wins and the residual census shows it
            if fl is None:  # no floor (``--no-floor``, or a noise-free toy): exact zero is the only pass
                ok_l, ok_c = el == 0.0, ec == 0.0
            else:
                ok_l = el <= row["floor_level"] + 2.0 * row["floor_level_sd"]
                ok_c = ec <= row["floor_comp"] + 2.0 * row["floor_comp_sd"]
            if ok_l and ok_c:
                row["winner"] = "either"
            elif ok_l or ok_c:
                row["winner"] = LEVEL if ok_l else COMP
            else:
                row["winner"] = LEVEL if xl < xc else (COMP if xc < xl else "tie")
            row["excess"] = min(xl, xc)
            rows.append(row)
    return rows


# ── the condition loader, and the gates on its inputs ─────────────────────────────────────────────


def gate_gdna_cannot_splice(spliced_g, sj_lo_g, sj_hi_g) -> dict:
    """**gdna-cannot-splice** — the gDNA origin's spliced and sj banks must be EXACTLY zero. If not, the
    SPLICE-IN population's gDNA part is not the boundary's unspliced gDNA and the oracle is wrong."""
    bad = float(np.sum(spliced_g) + np.sum(sj_lo_g) + np.sum(sj_hi_g))
    return {"gate": "gate gdna-cannot-splice", "ok": bool(bad == 0.0), "mass": bad}


def gate_components_close(M, n_rna_pos, n_rna_neg, atol: float = 1e-6) -> dict:
    """**components-close** — the three per-component partitions must reproduce the certified table's
    own RNA strand split at every slot. Two independent readings of the same reads (this instrument's
    geometries vs `calibration_oracle.py`'s), so a disagreement is a partition bug."""
    gap = np.abs(M[:, 1] - n_rna_pos) + np.abs(M[:, 2] - n_rna_neg)
    worst = float(gap.max()) if gap.size else 0.0
    return {"gate": "gate components-close", "ok": bool(worst <= atol), "worst_gap": worst}


def gate_hops_stay_in_reference(ref, dst, src) -> dict:
    """**hops-stay-in-reference** — the chain never links across a reference."""
    bad = int(np.sum(np.asarray(ref)[dst] != np.asarray(ref)[src]))
    return {"gate": "gate hops-stay-in-reference", "ok": bad == 0, "n_bad": bad}


def gate_truth_agrees(M, n_g, count_truth, n_g_truth) -> dict:
    """**truth-agrees** — the counts this instrument re-derives from the origin geometries equal the
    certified table's, slot for slot. One truth source, or the two tables describe different scans."""
    bad = int(np.sum(M != count_truth) + np.sum(n_g != n_g_truth))
    return {"gate": "gate truth-agrees", "ok": bad == 0, "n_bad": bad}


def load_condition(index, region_arrays, sj, bflags, suite: Path, condition: str) -> dict:
    """Everything one condition needs, as per-slot arrays. REFUSES a condition without a certified
    ``slot_truth.npz`` — the oracle it reads must itself have been gated."""
    truth_path = Path(suite) / "oracle_cache" / condition / "slot_truth.npz"
    if not truth_path.is_file():
        raise FileNotFoundError(
            f"{truth_path} is missing — run `calibration_oracle.py --condition {condition}` first; "
            "this instrument scores against CERTIFIED truth only."
        )
    truth = np.load(truth_path, allow_pickle=True)
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    kw = calibration_inputs(cache, index)
    chain = build_region_chain(cache.payload.ref_region_offsets, cache.payload.ref_boundary_offsets)
    statics = build_region_statics(chain, region_arrays, bflags)

    def geom_of(payload):
        return build_region_geometry(chain, CalibrationSubstrate.from_payload(payload, region_arrays),
                                     region_arrays, sj, kw["gdna_fl_pmf"], kw["rna_fl_pmf"], None)

    geom = geom_of(cache.payload)
    root = Path(suite) / "oracle_cache" / condition
    try:
        geoms = {k: geom_of(read_scan_cache(root / k, index).payload) for k in COMPONENT_ORIGIN}
    except Exception as exc:  # noqa: BLE001
        raise FileNotFoundError(
            f"{root} lacks the per-COMPONENT partitions {list(COMPONENT_ORIGIN)} "
            f"({type(exc).__name__}). `rna_pos`/`rna_neg` are the RNA reads split by TRANSCRIPT strand, "
            "built by `panel.py cache`; without them the three arms cannot be scored and a pooled-RNA "
            "map is not what the policy carries."
        ) from exc
    nrna_geom = geom_of(read_scan_cache(root / "nrna", index).payload)
    label = OC.strata(chain, statics, geom, region_arrays)["label"].astype(str)
    if not np.array_equal(label, truth["stratum"].astype(str)):
        raise AssertionError("the object classes re-derived here differ from the certified table's")

    def s1(a):
        return np.asarray(a, np.float64).sum(1)

    def per_component(attr):
        """``[n_slots, 3]`` — one column per :data:`COMPONENTS`, from that component's own partition."""
        return np.stack([s1(getattr(geoms[k], attr)) for k in COMPONENT_ORIGIN], axis=1)

    E_g = np.asarray(geom.eff_gdna, np.float64)
    E_r = np.asarray(geom.eff_rna, np.float64)
    d = {
        "condition": condition, "n_slots": int(chain.n_slots),
        "kind": np.asarray(chain.kind), "left": np.asarray(chain.left, np.int64),
        "right": np.asarray(chain.right, np.int64), "ref": np.asarray(truth["ref"]),
        "label": label,
        "structure": boundary_structure(np.asarray(chain.kind), np.asarray(statics.boundary_flags)),
        "ma": np.asarray(statics.mrna_active_pos, bool) | np.asarray(statics.mrna_active_neg, bool),
        # ⭐ every count carries the COMPONENT axis; the opportunity is `eff_gdna` for gDNA and
        #   `eff_rna` for BOTH RNA arms (the same geometry, a different strand)
        "M": per_component("unspliced_count"),
        "E": np.stack([E_g, E_r, E_r], axis=1),
        "spliced": per_component("spliced_count"),
        "sj_lo": per_component("sj_count_lo"),
        "sj_hi": per_component("sj_count_hi"),
        "field_certified": bool(truth["field_certified"]),
        "n_nrna_total": float(s1(nrna_geom.unspliced_count).sum()),
    }
    d["gates"] = [
        gate_truth_agrees(d["M"].sum(axis=1), d["M"][:, 0],
                          np.asarray(truth["count"], np.float64),
                          np.asarray(truth["n_gdna"], np.float64)),
        gate_gdna_cannot_splice(d["spliced"][:, 0], d["sj_lo"][:, 0], d["sj_hi"][:, 0]),
        gate_components_close(d["M"], np.asarray(truth["n_rna_pos"], np.float64),
                              np.asarray(truth["n_rna_neg"], np.float64)),
    ]
    return d


def measure(d: dict, rng) -> dict:
    """①②③ on one condition's arrays. Pure: no I/O, so the self-test feeds it synthetic slots."""
    dst, src = enumerate_hops(d["left"], d["right"])
    gates = list(d.get("gates", [])) + [gate_hops_stay_in_reference(d["ref"], dst, src)]
    htype = hop_type(refine_labels(d["label"], d["structure"]), dst, src)
    depth = depth_to_measured(d["left"], d["right"], d["ref"], d["ma"])
    dep = depth_census(depth, d["M"].sum(axis=1), d["ma"])
    pop = source_population(d["kind"], d["left"], d["right"], dst, src, d["M"], d["spliced"],
                            d["sj_lo"], d["sj_hi"])
    tr = transport(dst, src, d["M"], d["E"], pop)
    fl = noise_floor(dst, src, d["M"], d["E"], pop, rng) if rng is not None else None
    rows = per_type(htype, tr["scorable"], d["M"].sum(axis=1)[dst], tr, fl)
    return {"dst": dst, "src": src, "htype": htype, "depth": dep, "rows": rows, "gates": gates,
            "transport": tr, "floor": fl, "pop": pop}


# ── reporting ────────────────────────────────────────────────────────────────────────────────────


def _axes(cond: str) -> dict:
    ss, cap = PVO.stratum(cond)
    g = cond.split("_")[1]
    return {"gdna": g, "strand": ss, "capture": cap, "nrna": "mid" if "nrna_mid" in cond else "none",
            "scope": OC._scope(cond)}


def _k(x, width: int = 11) -> str:
    return f"{'—':>{width}}" if x is None or not np.isfinite(x) else f"{x:>{width},.0f}"


def _pct(x) -> str:
    return "     —" if x is None or not np.isfinite(x) else f"{100 * x:5.1f}%"


def report_condition(d: dict, m: dict, verbose: bool) -> None:
    ax = _axes(d["condition"])
    print(f"\n== {d['condition']}   [{ax['scope']}]   field-certified {d['field_certified']}   "
          f"nascent truth {d['n_nrna_total']:,.0f} fragments")
    for g in m["gates"]:
        print(f"   {'✔' if g['ok'] else '⛔'} {g['gate']}" + ("" if g["ok"] else f"   {g}"))
    dep = m["depth"]
    tot = dep["measured_mass"] + dep["imputed_mass"]
    print(f"   depth: gDNA-measuring mass {dep['measured_mass']:>13,.0f} ({100 * dep['measured_mass'] / max(tot, 1):.1f} %)"
          f"   imputed mass {dep['imputed_mass']:>13,.0f} ({100 * dep['imputed_mass'] / max(tot, 1):.1f} %)")
    im = max(dep["imputed_mass"], 1.0)
    print("   imputed mass by hops to the nearest measured slot: "
          + "  ".join(f"{k}: {100 * v / im:.1f}%" for k, v in dep["buckets"].items()))
    zero = PVO.is_zero_gdna(d["condition"])
    print(f"   {'dst <- src':<34} {'arm':<5} {'hops':>6} {'dst mass':>11} "
          f"{'LEVEL err':>11} {'floor':>9}  {'COMP err':>11} {'floor':>9}  winner  excess  clipped")
    for r in sorted(m["rows"], key=lambda r: (-r["dst_mass"], r["hop_type"], COMPONENTS.index(r["component"]))):
        if not verbose and r["n_scored"] == 0:
            continue
        stamp = "  (g00: gDNA exact by construction)" if zero and r["component"] == "gDNA" else ""
        label = r["hop_type"] if r["component"] == COMPONENTS[0] else ""
        print(f"   {label:<34} {r['component']:<5} {r['n_hops']:>6} {r['dst_mass']:>11,.0f} "
              f"{_k(r['err_level'])} {_k(r['floor_level'], 9)}  {_k(r['err_comp'])} {_k(r['floor_comp'], 9)}  "
              f"{r['winner']:<6} {_pct(r['excess'])}  {r['clipped']:>9,.0f}{stamp}")


def report_summary(results: list[tuple[dict, dict]]) -> None:
    """One table per (strand x capture x nascent): hop types down, gDNA levels across, each cell the two
    excesses and the winner. Never pooled across panel strata or across gDNA levels."""
    by = defaultdict(dict)
    types = set()
    for d, m in results:
        ax = _axes(d["condition"])
        for r in m["rows"]:
            if r["n_scored"] == 0:
                continue
            by[(ax["strand"], ax["capture"], ax["nrna"])].setdefault(
                (r["hop_type"], r["component"]), {})[ax["gdna"]] = r
            types.add(r["hop_type"])
    levels = ("g05", "g50", "g98")
    print("\n\n" + "=" * 140)
    print("SUMMARY — the currency per (hop type x ARM). Cell: LEVEL excess / COMP excess (each error minus")
    print("its own noise floor, as % of scored destination mass) -> the currency; 'either' = both at floor.")
    print("=" * 140)
    for key in sorted(by):
        strand, cap, nrna = key
        scope = OC._SCOPE[(strand, cap)]
        print(f"\n-- {strand} x {cap} x nrna_{nrna}   [{scope}]")
        print(f"   {'dst <- src':<34} {'arm':<5}" + "".join(f" | {g:^24}" for g in levels))
        tab = by[key]
        order = sorted(tab, key=lambda t: (-max(r["dst_mass"] for r in tab[t].values()),
                                           t[0], COMPONENTS.index(t[1])))
        for t in order:
            cells = []
            for g in levels:
                r = tab[t].get(g)
                if r is None:
                    cells.append(f" | {'—':^24}")
                    continue
                cells.append(f" | {_pct(r['excess_level'])} {_pct(r['excess_comp'])} {r['winner']:>6}")
            print(f"   {t[0] if t[1] == COMPONENTS[0] else '':<34} {t[1]:<5}" + "".join(cells))


def report_residuals(results: list[tuple[dict, dict]]) -> None:
    """③ the hop types where the WINNING currency's excess over its floor is still large, ranked, on
    the contaminated conditions — named, not explained."""
    rows = []
    for d, m in results:
        if PVO.is_zero_gdna(d["condition"]):
            continue
        ax = _axes(d["condition"])
        for r in m["rows"]:
            if r["n_scored"] and np.isfinite(r["excess"]):
                rows.append((r["excess"], r["excess"] * r["mass_scored"], ax, r))
    rows.sort(key=lambda x: -x[1])
    print("\n\n" + "=" * 140)
    print("③ RESIDUAL CENSUS — hop type x condition ranked by the WINNER's excess over its floor, in fragments")
    print("   (a row here is a hop where NEITHER currency is enough: a third MECHANISM, to be named before proposed)")
    print("=" * 140)
    print(f"   {'dst <- src':<34} {'arm':<5} {'condition':<26} {'winner':>6} {'excess':>8} {'frags':>11}")
    for ex, frag, ax, r in rows[:40]:
        cond = f"{ax['gdna']} {ax['strand'][:3]} {ax['capture'][-3:]:<3}"
        print(f"   {r['hop_type']:<34} {r['component']:<5} {cond:<26} {r['winner']:>6} {_pct(ex):>8} {frag:>11,.0f}")


def write_tsv(path: Path, results: list[tuple[dict, dict]]) -> None:
    cols = ["condition", "gdna", "strand", "capture", "nrna", "scope", "hop_type", "component", "n_hops", "dst_mass",
            "n_scored", "mass_scored", "err_level", "floor_level", "floor_level_sd", "err_comp",
            "floor_comp", "floor_comp_sd", "clipped", "winner", "excess_level", "excess_comp", "excess"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for d, m in results:
            ax = _axes(d["condition"])
            for r in m["rows"]:
                vals = {**ax, **r, "condition": d["condition"]}
                fh.write("\t".join(str(vals[c]) for c in cols) + "\n")
    print(f"\n→ {path}")


# ── self-test: every comparator perturbed, no I/O ────────────────────────────────────────────────


def _toy(n_g_exon_frac: float = 0.3, enrich_exon: float = 1.0, sj_lo_at_5: float = 0.0,
         sj_hi_at_5: float = 0.0, spliced_g_at_5: float = 0.0, two_refs: bool = True,
         n_refs: int | None = None, rho: float = 1.0, neg_share: float = 0.0) -> dict:
    """A hand-built chain: ``intergenic | edge | exon | e|i | intron | e|i | exon | edge | intergenic``
    per reference, with EXACT (noise-free) counts so the oracle's errors are structural, not sampled.
    gDNA: density 1.0 everywhere, times ``enrich_exon`` on exon REGIONs; RNA: only in exons, at a share
    that makes the exon's true f_g == ``n_g_exon_frac`` when ``enrich_exon == 1``.
    """
    labels = ["R intergenic", "B gene edge", "R exon", "B exon|intron", "R intron",
              "B exon|intron", "R exon", "B gene edge", "R intergenic"]
    kinds = [REGION if lab.startswith("R") else BOUNDARY for lab in labels]
    E = np.array([1000.0, 100.0, 200.0, 100.0, 500.0, 100.0, 200.0, 100.0, 1000.0])
    if n_refs is None:
        n_refs = 2 if two_refs else 1
    n = len(labels) * n_refs
    label = np.array(labels * n_refs)
    kind = np.array(kinds * n_refs, np.int8)
    ref = np.repeat(np.arange(n_refs), len(labels))
    left = np.arange(n) - 1
    right = np.arange(n) + 1
    left[ref != np.roll(ref, 1)] = -1
    left[0] = -1
    right[ref != np.roll(ref, -1)] = -1
    right[-1] = -1
    E_g = np.tile(E, n_refs)
    n_g = rho * E_g
    is_exon = label == "R exon"
    n_g[is_exon] *= enrich_exon
    rna = np.zeros(n)
    rna[is_exon] = E_g[is_exon] * (1.0 - n_g_exon_frac) / n_g_exon_frac  # f_g = n_g/(n_g+rna)
    ma = np.isin(label, ["R exon", "B exon|exon"])
    structure = np.array(["", "term", "", "sj", "", "sj", "", "term", ""] * n_refs)
    # ⭐ the COMPONENT axis: [gDNA, RNA+, RNA−]. ``neg_share`` moves that share of the RNA onto the
    #   minus arm, which is what exercises the two RNA arms against each other.
    M = np.stack([n_g, rna * (1.0 - neg_share), rna * neg_share], axis=1)
    E_arr = np.stack([E_g, E_g, E_g], axis=1)
    z = np.zeros((n, 3))
    d = {"condition": "toy", "n_slots": n, "kind": kind, "left": left, "right": right, "ref": ref,
         "label": label, "structure": structure, "ma": ma, "M": M, "E": E_arr,
         "spliced": z.copy(), "sj_lo": z.copy(), "sj_hi": z.copy(),
         "field_certified": True, "n_nrna_total": 0.0}
    d["sj_lo"][5, 1] = sj_lo_at_5
    d["sj_hi"][5, 1] = sj_hi_at_5
    d["spliced"][5, 0] = spliced_g_at_5  # gDNA cannot splice — the gate's perturbation handle
    d["gates"] = [gate_gdna_cannot_splice(d["spliced"][:, 0], d["sj_lo"][:, 0], d["sj_hi"][:, 0])]
    return d


def self_test() -> int:
    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"   ⛔ {name}")

    def rows_of(m, comp: str = "gDNA"):
        """One component's rows, keyed by hop type — ``gDNA`` by default, so every check written for
        the single-arm map keeps testing exactly what it tested."""
        return {r["hop_type"]: r for r in m["rows"] if r["component"] == comp}

    # ① the census partitions the hops: both directions of every adjacency, once each
    d = _toy()
    m = measure(d, None)
    n_adj = int((d["right"] >= 0).sum())
    check("every ordered adjacent pair is one hop, n_hops == 2 * adjacencies",
          m["dst"].shape[0] == 2 * n_adj == 2 * (d["n_slots"] - 2))
    check("hops stay inside a reference", all(g["ok"] for g in m["gates"]))
    check("hop types sum to the hop count",
          sum(r["n_hops"] for r in rows_of(m).values()) == m["dst"].shape[0])
    one_ref = measure(_toy(two_refs=False), None)
    check("a second reference doubles the hops and adds none across the gap",
          m["dst"].shape[0] == 2 * one_ref["dst"].shape[0])

    # ① depth: exons are 1 hop from their flanking boundaries (which are measured here); the two exons
    #   are the only ma slots, so depth mass is all at depth 1
    check("depth census: every imputed slot is at depth 1 on this chain",
          m["depth"]["buckets"]["1"] == m["depth"]["imputed_mass"] > 0)
    # a chain whose middle is a long exonic stretch: exon, e|e, exon, e|e, exon ... depth grows
    dd = _toy()
    dd["label"] = dd["label"].copy()
    dd["ma"] = dd["ma"].copy()
    dd["label"][3:6] = ["B exon|exon", "R exon", "B exon|exon"]
    dd["ma"][3:6] = True
    md = measure(dd, None)
    b = md["depth"]["buckets"]
    check("depth census: a 5-slot imputed stretch reads depths 1,2,3,2,1 (mass-weighted)",
          b["1"] == dd["M"][2].sum() + dd["M"][6].sum() + dd["M"][11].sum() + dd["M"][15].sum()
          and b["2"] == dd["M"][3].sum() + dd["M"][5].sum() and b["3"] == dd["M"][4].sum())
    # depth unreachable when the whole reference is imputed
    du = _toy(two_refs=False)
    du["ma"] = np.ones(du["n_slots"], bool)
    check("depth census: a wholly imputed reference is UNREACHABLE, not depth 0",
          measure(du, None)["depth"]["buckets"]["unreachable"] == du["M"].sum())

    # ② LEVEL is EXACT on a uniform field (no noise), at EVERY hop type; COMP is exact only where the
    #   population does not change
    r = rows_of(m)
    check("LEVEL error is 0 at every hop on a uniform noise-free field",
          all(row["err_level"] == 0.0 for row in r.values()))
    check("COMP error is 0 intron -> e|i (same population, pure gDNA both)", r["B exon|intron[sj] <- R intron"]["err_comp"] == 0.0)
    check("COMP error is > 0 e|i -> exon (raw unspliced crossing is pure gDNA, exon holds RNA)",
          r["R exon <- B exon|intron[sj]"]["err_comp"] > 0.0)
    check("COMP error e|i -> exon equals the exon's RNA mass (the whole population change)",
          abs(r["R exon <- B exon|intron[sj]"]["err_comp"] - 2 * 2 * d["M"][2, 1:].sum()) < 1e-9)
    check("winner named LEVEL where COMP fails and LEVEL is exact", r["R exon <- B exon|intron[sj]"]["winner"] == LEVEL)

    # ② an ENRICHMENT edge breaks LEVEL at exactly the hops that cross it, and nowhere else
    de = _toy(enrich_exon=4.0)
    re_ = rows_of(measure(de, None))
    check("enrichment: LEVEL fails into the exon", re_["R exon <- B exon|intron[sj]"]["err_level"] > 0.0)
    # out of the enriched exon the claim is 4x the boundary's whole mass; clipped, it says "all gDNA" —
    # which is TRUE of a pure-gDNA boundary (err 0, clipped > 0), and FALSE the moment the boundary
    # holds RNA: then every RNA fragment there is called gDNA
    check("enrichment: LEVEL out of the exon into a pure-gDNA boundary is clipped to a TRUE claim",
          re_["B exon|intron[sj] <- R exon"]["err_level"] == 0.0)
    dr = _toy(enrich_exon=4.0)
    dr["M"][5, 1] += 50.0  # 50 RNA crossings at boundary 5
    rr = rows_of(measure(dr, None))
    check("enrichment: the same over-claim into a boundary holding RNA calls all of it gDNA",
          rr["B exon|intron[sj] <- R exon"]["err_level"] == 50.0)
    check("enrichment: LEVEL still exact intron -> e|i (no probe edge there)",
          re_["B exon|intron[sj] <- R intron"]["err_level"] == 0.0)
    check("enrichment: a LEVEL over-claim out of the exon is CLIPPED and counted",
          re_["B exon|intron[sj] <- R exon"]["clipped"] > 0.0 and re_["B exon|intron[sj] <- R intron"]["clipped"] == 0.0)

    # ② the SPLICE-IN population: sj_lo at boundary 5 enters its LEFT neighbour (the intron, slot 4),
    #   sj_hi its RIGHT neighbour (the exon, slot 6) — and the side rule is falsifiable
    base = rows_of(measure(_toy(), None))["R exon <- B exon|intron[sj]"]["err_comp"]
    hi = rows_of(measure(_toy(sj_hi_at_5=1000.0), None))["R exon <- B exon|intron[sj]"]["err_comp"]
    lo = rows_of(measure(_toy(sj_lo_at_5=1000.0), None))["R exon <- B exon|intron[sj]"]["err_comp"]
    check("sj flux on the exon's side (hi at slot 5 -> slot 6) lowers COMP's error into the exon", hi < base)
    check("sj flux on the OTHER side (lo at slot 5 -> slot 4) leaves the exon hop untouched", lo == base)
    lo_int = rows_of(measure(_toy(sj_lo_at_5=1000.0), None))["R intron <- B exon|intron[sj]"]["err_comp"]
    base_int = rows_of(measure(_toy(), None))["R intron <- B exon|intron[sj]"]["err_comp"]
    check("…and the same lo flux DOES change the hop into the intron it enters", lo_int != base_int)
    # the exact SPLICE-IN arithmetic: with sj_hi = rna(exon) * E_cross/E_contained the composition
    # entering equals the exon's own, so the error vanishes
    dx = _toy()
    rna_exon = dx["M"][6, 1]
    dx["sj_hi"][5, 1] = rna_exon * dx["E"][5, 1] / dx["E"][6, 1]
    dx["sj_hi"][14, 1] = rna_exon * dx["E"][14, 1] / dx["E"][15, 1]
    mx = measure(dx, None)
    sel = (mx["htype"] == "R exon <- B exon|intron[sj]") & (mx["dst"] == 6)
    check("SPLICE-IN at the matching density reproduces the exon's composition exactly",
          np.allclose(mx["transport"]["err_comp"][sel], 0.0))
    check("gate gdna-cannot-splice fires on 1 spliced gDNA fragment",
          not measure(_toy(spliced_g_at_5=1.0), None)["gates"][0]["ok"])

    # ① the structural split: a boundary's flags name its class, a REGION's stays empty, and a boundary
    #   with neither bit raises
    from rigel.calibration.splice_graph import FLAG_ACCEPTOR_NEG, FLAG_DONOR_POS, FLAG_TES_NEG, FLAG_TSS_POS
    kk = np.array([REGION, BOUNDARY, BOUNDARY, BOUNDARY, REGION], np.int8)
    ff = np.array([0, FLAG_TSS_POS, FLAG_DONOR_POS, FLAG_TES_NEG | FLAG_ACCEPTOR_NEG, 0], np.uint16)
    check("boundary_structure reads term / sj / sj+term from the flags, '' on regions",
          boundary_structure(kk, ff).tolist() == ["", "term", "sj", "sj+term", ""])
    check("refine_labels suffixes boundaries only",
          refine_labels(np.array(["R exon", "B gene edge"]), np.array(["", "term"])).tolist() == ["R exon", "B gene edge[term]"])
    try:
        boundary_structure(kk, np.zeros(5, np.uint16))
        check("a flagless boundary raises", False)
    except AssertionError:
        check("a flagless boundary raises", True)
    # and the split is what separates the hop types: relabel boundary 5 as a terminus and the exon hop
    # from it moves class, carrying its error with it
    dt = _toy()
    dt["structure"] = dt["structure"].copy()
    dt["structure"][5] = "term"
    rt = rows_of(measure(dt, None))
    check("relabelling one boundary as a terminus moves its hops into the [term] class",
          "R exon <- B exon|intron[term]" in rt and rt["R exon <- B exon|intron[term]"]["n_hops"] == 1
          and rt["R exon <- B exon|intron[sj]"]["n_hops"] == r["R exon <- B exon|intron[sj]"]["n_hops"] - 1)

    # ② scorability: an empty destination contributes no score; a source with E_g = 0 supplies no LEVEL
    dz = _toy()
    dz["M"][2, :] = 0.0
    mz = measure(dz, None)
    t = "R exon <- B gene edge[term]"
    check("an empty destination is excluded from the scored hops and the scored mass",
          rows_of(mz)[t]["n_scored"] == rows_of(m)[t]["n_scored"] - 1
          and rows_of(mz)[t]["mass_scored"] == rows_of(m)[t]["mass_scored"] - d["M"][2].sum())
    dn = _toy()
    dn["E"][4, :] = 0.0
    mn = measure(dn, None)
    check("a source with no gDNA opportunity supplies no LEVEL and the hop is not scored",
          rows_of(mn)["B exon|intron[sj] <- R intron"]["n_scored"] == rows_of(m)["B exon|intron[sj] <- R intron"]["n_scored"] - 2)

    # the noise floor: positive where counts are, zero where the null density is zero, and DISTINCT
    # from the error (a noise-free exact field has error 0 and floor > 0)
    rng = np.random.default_rng(0)
    mf = measure(_toy(), rng)
    rf = rows_of(mf)
    check("floor > 0 on a populated hop while the noise-free error is exactly 0",
          rf["B exon|intron[sj] <- R intron"]["floor_level"] > 0.0
          and rf["B exon|intron[sj] <- R intron"]["err_level"] == 0.0)
    dz0 = _toy()
    dz0["M"][:, 0] = 0.0  # no gDNA anywhere — the `g00` shape
    rz = rows_of(measure(dz0, np.random.default_rng(0)))
    check("the gDNA arm's floor is 0 where there is no gDNA at all (the g00 shape)",
          all(row["floor_level"] == 0.0 and row["floor_comp"] == 0.0 for row in rz.values()))
    # and the floor is what a Poisson-noisy PERFECT rule actually scores: draw one noisy realisation
    # of the exact toy and its LEVEL error must sit within the floor's own spread
    # ⚠ scored on the gDNA ARM only, and that is not a convenience: the toy's gDNA field is uniform by
    #   construction, so a LEVEL is the perfect rule for it — while its RNA lives ONLY in exons, so a
    #   LEVEL is structurally wrong for the RNA arms there and would not be "a perfect rule" at all
    dp = _toy()
    g = np.random.default_rng(1)
    dp["M"][:, 0] = g.poisson(dp["M"][:, 0]).astype(float)
    mp = measure(dp, np.random.default_rng(2))
    gd = [r_ for r_ in mp["rows"] if r_["component"] == "gDNA"]
    tot_err = sum(r_["err_level"] for r_ in gd)
    tot_floor = sum(r_["floor_level"] for r_ in gd)
    check("a noisy perfect rule scores within 2x of its floor (gDNA arm, totals)",
          0.5 * tot_floor < tot_err < 2.0 * tot_floor)
    # ⭐ and the RNA arms are where a LEVEL is structurally wrong on this toy — the check that the
    #   per-arm split is doing something rather than copying the gDNA answer
    rna_err = sum(r_["err_level"] for r_ in mp["rows"] if r_["component"] == "RNA+")
    check("the RNA+ arm's LEVEL error is far above the gDNA arm's — the arms are scored separately",
          rna_err > 5.0 * tot_err)

    # ⛔ the floor must reproduce the SELECTION: at ~1 expected crossing per boundary, scoring only
    #   non-empty destinations truncates the Poisson; the truncated null matches a perfect rule's
    #   realized error and the unconditional one reads the truncation as a rule error
    dq = _toy(n_refs=400, rho=0.012)  # boundary expectation 1.2, like g05
    g = np.random.default_rng(3)
    dq["M"][:, 0] = g.poisson(dq["M"][:, 0]).astype(float)
    dq_dst, dq_src = enumerate_hops(dq["left"], dq["right"])
    hq = hop_type(refine_labels(dq["label"], dq["structure"]), dq_dst, dq_src)
    pq = source_population(dq["kind"], dq["left"], dq["right"], dq_dst, dq_src, dq["M"],
                           dq["spliced"], dq["sj_lo"], dq["sj_hi"])
    tq = transport(dq_dst, dq_src, dq["M"], dq["E"], pq)
    into_edge = (hq == "B gene edge[term] <- R intergenic") & tq["scorable"]
    f_tr = noise_floor(dq_dst, dq_src, dq["M"], dq["E"], pq, np.random.default_rng(4), truncate=True)
    f_un = noise_floor(dq_dst, dq_src, dq["M"], dq["E"], pq, np.random.default_rng(4), truncate=False)
    e_real = tq["err_level"][into_edge, 0].sum()
    e_tr = f_tr["floor_level"][into_edge, 0].sum()
    e_un = f_un["floor_level"][into_edge, 0].sum()
    check("truncated LEVEL floor sits within 15 % of a perfect rule's realized error at ~1 crossing/boundary",
          abs(e_tr - e_real) < 0.15 * e_real)
    check("…and the UNconditional floor reads that truncation as a >20 % 'rule error'",
          e_un < 0.8 * e_real)

    # ⭐⭐ THE PER-ARM SPLIT ITSELF: with RNA on BOTH strands the two RNA arms must carry their own
    #   mass and their own errors, and moving mass between them must move only those two rows
    d2 = _toy(neg_share=0.25)
    m2 = measure(d2, None)
    pos = {r["hop_type"]: r for r in m2["rows"] if r["component"] == "RNA+"}
    neg = {r["hop_type"]: r for r in m2["rows"] if r["component"] == "RNA-"}
    t2 = "R exon <- B exon|intron[sj]"
    check("both RNA arms carry mass when the toy is two-stranded",
          pos[t2]["err_comp"] > 0.0 and neg[t2]["err_comp"] > 0.0)
    check("the minus arm is a QUARTER of the plus arm's error, as its share says",
          abs(neg[t2]["err_comp"] / pos[t2]["err_comp"] - 0.25 / 0.75) < 1e-9)
    one = {r["hop_type"]: r for r in measure(_toy(), None)["rows"] if r["component"] == "RNA-"}
    check("a single-stranded toy leaves the minus arm empty, reported not invented",
          one[t2]["err_comp"] == 0.0 and one[t2]["dst_mass"] > 0.0)
    gd2 = {r["hop_type"]: r for r in m2["rows"] if r["component"] == "gDNA"}
    gd1 = {r["hop_type"]: r for r in measure(_toy(), None)["rows"] if r["component"] == "gDNA"}
    check("moving RNA between the strands does not touch the gDNA arm",
          abs(gd2[t2]["err_comp"] - gd1[t2]["err_comp"]) < 1e-9)

    # ⛔ every hop type reports exactly the three arms, in order — a missing arm would read as a hop
    #   type that simply has no RNA rather than as a scoring hole
    per_hop = {}
    for r in m2["rows"]:
        per_hop.setdefault(r["hop_type"], []).append(r["component"])
    check("every hop type carries all three arms in COMPONENTS order",
          all(v == list(COMPONENTS) for v in per_hop.values()))

    # the winner and excess are derived from the same numbers the table prints
    # the 'either' verdict: both within 2 sd of their floors on the noise-free exact toy? No — there the
    # error is 0 and the floor > 0, so both are consistent with the null and the cell says 'either'
    check("a noise-free uniform field reads 'either' where both currencies are exact",
          rf["B exon|intron[sj] <- R intron"]["winner"] == "either")
    check("…and LEVEL where COMP's error is structural", rf["R exon <- B gene edge[term]"]["winner"] == LEVEL)
    check("excess is the winner's error over its floor as a share of scored mass",
          all(abs(row["excess"] - min(row["excess_level"], row["excess_comp"])) < 1e-12
              for row in rf.values() if np.isfinite(row["excess"])))

    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None, help="one condition (default: every scan_cache entry)")
    ap.add_argument("--out", type=Path, default=None, help="write every (condition x hop type) row as TSV")
    ap.add_argument("--verbose", action="store_true", help="print hop types with no scorable hop too")
    ap.add_argument("--no-floor", action="store_true", help="skip the Monte-Carlo noise floor")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    index = TranscriptIndex.load(args.index)
    region_arrays = RegionArrays.from_index(index)
    sj = build_sj_geometry_arrays(index)
    bflags = build_boundary_flags_array(index)
    conds = [args.condition] if args.condition else sorted(p.name for p in (args.suite / "scan_cache").iterdir())
    results = []
    bad = 0
    for c in conds:
        d = load_condition(index, region_arrays, sj, bflags, args.suite, c)
        m = measure(d, None if args.no_floor else np.random.default_rng(0))
        report_condition(d, m, args.verbose)
        bad += sum(not g["ok"] for g in m["gates"])
        results.append((d, m))
    if len(results) > 1:
        report_summary(results)
        report_residuals(results)
    if args.out:
        write_tsv(args.out, results)
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
