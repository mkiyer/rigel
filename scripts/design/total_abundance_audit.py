"""IS THE MEASURED TOTAL A TRUE TOTAL? — the composition-free per-slot abundance against certified
truth, on every condition, with the wall mask checked against an independent enumeration.

⭐⭐⭐ **THE VALIDATION RUNG. A MEASUREMENT: no solver runs, no EM, no BAM re-scan, and nothing in
``src/`` is patched.** ``rigel.calibration.total_abundance`` claims that a slot's TOTAL fragment
density is measurable with no composition model in it — the START/END region banks over ``ell``,
side-selected by the wall rule, plus the exact reciprocal-opportunity banks and the certified spliced
incidence at a BOUNDARY. This instrument is what decides whether that claim survives contact with the
certified per-slot truth, before any consumer reads the number.

Five arms, and they answer different questions on purpose — read ⓔ first.

**ⓔ THE START/END AGREEMENT ARM — the decisive test of the WALL RULE, and the only one that is
FIELD-FREE.** ``S/ell`` and ``E/ell`` estimate the SAME per-region rate with the SAME opportunity
``ell`` for every fragment length, reading whatever field capture produced at that very slot. So no
uniformity assumption and no length distribution enter, and the rule's content is a DIRECTION::

    both sides exact         ->  the two must AGREE            ( Sum E / Sum S = 1 )
    start exact, end BOUND   ->  E is depressed                ( < 1 )
    end exact, start BOUND   ->  S is depressed                ( > 1 )

⛔ A mask assigned at random reads ~1 on all three; a mask with the sides BACKWARDS inverts the last
two. ⭐ **Measured on the ladder under capture — the hardest case — 1.0060 / 0.8018 / 1.2461.**

**ⓐ THE TRUNCATION LAW — the START bank is a TOTAL where the CONTAINED bank is a SHAPE.** Per
component and per region the two are estimators of the same density differing only by support::

    E[S_c / ell]      =  rho_c                     <- the total
    E[Cinv_c]         =  rho_c * P_c(w <= ell)     <- the shipped bank, truncated

so their ratio must reproduce ``P_c(w <= ell)``, which each origin partition's own deposit histogram
states exactly. ⛔⛔ **IT ASSUMES LOCAL FIELD UNIFORMITY AND CAPTURE VIOLATES THAT, WHICH IS A FINDING
RATHER THAN A FAILURE**: the two banks weight the region's positions differently (a start anywhere in
the region, versus a fragment wholly inside it), so under a probe-shaped field they are different
weighted averages of it — legitimately, and by up to 30x on intron/intergenic regions. Read it at
capture-OFF, where it is the claim's confirmation (**1.00025 / 1.00045** on the ladder over 14.4 M
fragments), and read ⓔ under capture. Aggregated as a RATIO OF SUMS inside ``ell`` bins
(``TRAPS: a-mean-of-ratios-inherits-the-partition``).

**ⓑ THE UNIFORM-FIELD ARM — the absolute total, in FRAGMENTS, and STAMPED VACUOUS under capture.**
Against the per-reference per-component realized rate ``rho_c = Sigma S_c / Sigma ell``, whose sum over
components is the slot's true total density where the field is uniform. Error is reported as
``|measured - true| * ell`` — misplaced fragments, not a ratio. ⛔ Under capture the gDNA field is
deliberately non-uniform, so this arm is not a measurement there and says so, exactly as
``calibration_oracle.py`` stamps its own field gate.

**ⓒ THE WALL MASK against an INDEPENDENT ENUMERATION.** The shipped mask reads its mature distances
from a vectorised kernel; this arm recomputes them with a per-transcript Python loop over
``index.get_exon_intervals`` — a second implementation sharing no helper with the first — and demands
EXACT agreement. It then reports the exposure (what fraction of START mass sits on a wall-bound or
double-walled region) so the wall is a number rather than a worry.

**ⓓ THE BOUNDARY ARM.** A crossing's opportunity is ``w - 1`` and the deposit ``1/(w - 1)``, so
``E[Sigma] = rho`` for any real library — no support factor and no wall. Scored against the same
uniform-field truth, and its residual is the control for ⓑ: the boundary axis has no truncation and no
wall, so an error there is a field error and not an estimator error.

⛔ Every arm is scored per stratum with the 0.8.0 scope stamped, ``g00`` on its own row, and the
DEFERRED stratum reported rather than dropped. Gates refuse the run rather than warning: the ledger
must close on every payload, the partitions must sum to the whole, and the mask's two implementations
must agree.

Usage::

    python scripts/design/total_abundance_audit.py                    # every cached condition
    python scripts/design/total_abundance_audit.py --condition NAME   # one condition, full detail
    python scripts/design/total_abundance_audit.py --suite ... --index ...   # the other panel
    python scripts/design/total_abundance_audit.py --out rows.tsv     # one row per condition x arm
    python scripts/design/total_abundance_audit.py --self-test        # perturbed, no I/O
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
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


OC = _sibling("object_composition.py")
PVO = OC.PVO

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import (  # noqa: E402
    build_region_geometry,
    build_region_statics,
)
from rigel.calibration.splice_graph import (  # noqa: E402
    MatureWallDistances,
    build_boundary_flags_array,
    build_contiguous_boundary_reach_arrays,
    build_mature_wall_distances,
    build_sj_geometry_arrays,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.calibration.total_abundance import (  # noqa: E402
    build_region_wall_mask,
    build_total_abundance,
    w_max_from_deposited_lengths,
)
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402
from rigel.types import Strand  # noqa: E402

DEFAULT_SUITE = OC.DEFAULT_SUITE
DEFAULT_INDEX = OC.DEFAULT_INDEX

#: the three origin partitions whose START banks must sum to the whole. ⭐ gdna/mrna/nrna, not the
#: transcript-strand split: the truncation law is per POPULATION and the fragment-length pmf that
#: truncates it is that population's own.
ORIGINS = ("gdna", "mrna", "nrna")
#: ell bins for the truncation law, in bp. ⚠ Bins rather than per-region ratios because the law is a
#: ratio of expectations: a per-region ratio of two small counts has no usable distribution.
ELL_BINS = (0, 50, 100, 200, 400, 800, 1600, 1 << 62)


# ---------------------------------------------------------------------------
# ⓒ the independent enumeration — a per-transcript Python loop, sharing no helper with the kernel
# ---------------------------------------------------------------------------


def enumerate_mature_walls(index, region_arrays) -> MatureWallDistances:
    """The mature wall distances, recomputed the slow obvious way: for every transcript, walk its own
    exons, accumulate spliced offsets by hand, and write MAX per (region, strand).

    ⛔ Deliberately NOT vectorised and deliberately not sharing a line with
    :func:`~rigel.calibration.splice_graph.mature_wall_distances_kernel` — an enumeration that reuses
    the implementation's own machinery gates nothing. This is the RANK-3 collapse stated directly:
    ``d_high(r) = MAX`` downstream spliced template over the templates covering ``r``.
    """
    starts = np.asarray(region_arrays.start, np.int64)
    ends = np.asarray(region_arrays.end, np.int64)
    ref_id = np.asarray(region_arrays.ref_id, np.int64)
    n = starts.shape[0]
    d_low = np.zeros((n, 2), np.float64)
    d_high = np.zeros((n, 2), np.float64)
    covered = np.zeros((n, 2), bool)

    name_to_id = index.ref_name_to_id
    tdf = index.t_df
    if tdf is None:
        return MatureWallDistances(d_low=d_low, d_high=d_high, covered=covered)

    # a plain per-reference list of (start, end, row); no CSR arithmetic, so an offset bug cannot be
    # shared with the kernel. ⚠ Sorted and bisected rather than scanned: at panel scale a full scan
    # per exon is ~1e8 Python steps. The bisect locates the first candidate; the walk below is an
    # explicit loop, which is what keeps this a second implementation rather than a second call.
    import bisect

    rows_of_ref: dict[int, list[tuple[int, int, int]]] = {}
    for i in range(n):
        rows_of_ref.setdefault(int(ref_id[i]), []).append((int(starts[i]), int(ends[i]), i))
    for v in rows_of_ref.values():
        v.sort()
    starts_of_ref = {k: [x[0] for x in v] for k, v in rows_of_ref.items()}

    # ⛔ SYNTHETIC spans are excluded, and this filter is the POPULATION rule rather than a detail:
    # a synthetic entity is a manufactured nascent template and a nascent molecule extends
    # GENOMICALLY, which is the contiguous reach's arm. ⚠ `get_exon_intervals` DOES return a
    # synthetic span's own interval, so omitting this reads 9 nascent spans as mature templates —
    # measured, on the test chromosome, as 57 differing distances against the shipped kernel.
    synthetic = (
        tdf["is_synthetic"].to_numpy(dtype=bool)
        if "is_synthetic" in tdf.columns
        else np.zeros(len(tdf), bool)
    )
    for t_index, ref_name, strand, is_syn in zip(
        tdf["t_index"].to_numpy(np.int64),
        tdf["ref"].to_numpy(),
        tdf["strand"].to_numpy(np.int64),
        synthetic,
        strict=True,
    ):
        if is_syn:
            continue
        exons = index.get_exon_intervals(int(t_index))
        if exons is None or len(exons) == 0:
            continue
        rid = name_to_id.get(str(ref_name))
        if rid is None:
            continue
        col = 0 if int(strand) == int(Strand.POS) else 1
        blocks = sorted((int(a), int(b)) for a, b in exons if int(b) > int(a))
        total = sum(b - a for a, b in blocks)
        seen = 0
        rows = rows_of_ref.get(int(rid), ())
        keys = starts_of_ref.get(int(rid), ())
        for a, b in blocks:
            j = bisect.bisect_left(keys, a)
            while j < len(rows):
                lo, hi, i = rows[j]
                if lo >= b:
                    break
                j += 1
                if lo < a or hi > b:
                    continue  # this region is not inside this exon
                below = seen + (lo - a)
                above = total - seen - (hi - a)
                covered[i, col] = True
                d_low[i, col] = max(d_low[i, col], float(below))
                d_high[i, col] = max(d_high[i, col], float(above))
            seen += b - a
    return MatureWallDistances(d_low=d_low, d_high=d_high, covered=covered)


# ---------------------------------------------------------------------------
# the per-condition measurement
# ---------------------------------------------------------------------------


def _pmf_cdf(deposited_lengths: np.ndarray, ell: np.ndarray) -> np.ndarray:
    """``P(w <= ell)`` from a deposit histogram, evaluated at each region length."""
    h = np.asarray(deposited_lengths, np.float64)
    tot = h.sum()
    if tot <= 0.0:
        return np.zeros(ell.shape[0], np.float64)
    cdf = np.cumsum(h) / tot
    idx = np.clip(np.asarray(ell, np.int64), 0, cdf.shape[0] - 1)
    return cdf[idx]


def _unspliced_hist(payload) -> np.ndarray:
    """This partition's deposit histogram with the SPLICED pool removed — the support the contained
    bank actually has. ⛔ Refuses a negative residue rather than clipping it: the spliced pool is a
    subset of the deposits by construction, and a negative would mean the two disagree."""
    from rigel.scan_payload import POOL_RNA_SPLICED

    dep = np.asarray(payload.deposited_lengths, np.float64)
    spl = np.asarray(payload.pool_lengths, np.float64)[POOL_RNA_SPLICED]
    out = dep - spl[: dep.shape[0]]
    if np.any(out < -1e-9):
        raise ValueError(
            "the spliced pool exceeds the deposit histogram at some length — the two describe "
            "different populations and the truncation law cannot be formed."
        )
    return np.maximum(out, 0.0)


def _spliced_share(payload) -> float:
    from rigel.scan_payload import POOL_RNA_SPLICED

    dep = float(np.asarray(payload.deposited_lengths, np.float64).sum())
    if dep <= 0.0:
        return float("nan")
    return float(np.asarray(payload.pool_lengths, np.float64)[POOL_RNA_SPLICED].sum() / dep)


def _ratio_of_sums(num: np.ndarray, den: np.ndarray, select: np.ndarray) -> float:
    d = float(np.sum(np.asarray(den, np.float64)[select]))
    if d <= 0.0:
        return float("nan")
    return float(np.sum(np.asarray(num, np.float64)[select]) / d)


def measure(payload, parts, index_parts, condition: str) -> dict:
    """One condition. ``parts`` maps origin -> that partition's payload; ``index_parts`` carries the
    index-derived arrays (built once per index by the caller). PURE apart from the arithmetic, so the
    self-test feeds it synthetic payloads."""
    ra = index_parts["region_arrays"]
    chain = index_parts["chain"]
    geometry = index_parts["geometry"]
    substrate = index_parts["substrate"]
    mature = index_parts["mature"]
    reach_lo, reach_hi = index_parts["reach"]
    label = index_parts["label"]
    rna_pmf = index_parts["rna_pmf"]

    gates: list[dict] = []
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    is_region = kind == REGION
    is_boundary = kind == BOUNDARY
    ell = np.asarray(ra.region_size_bp, np.float64)
    ref_id = np.asarray(ra.ref_id, np.int64)

    S = {k: np.asarray(p.region_start_count, np.float64).sum(1) for k, p in parts.items()}
    E = {k: np.asarray(p.region_end_count, np.float64).sum(1) for k, p in parts.items()}
    C = {
        k: np.asarray(p.region_contained_inv_opportunity_sum, np.float64) for k, p in parts.items()
    }
    S_all = np.asarray(payload.region_start_count, np.float64).sum(1)
    E_all = np.asarray(payload.region_end_count, np.float64).sum(1)

    # ── gate: the ledger closes on every payload, twice over
    for name, p in (("full", payload), *parts.items()):
        s = float(np.asarray(p.region_start_count, np.float64).sum())
        e = float(np.asarray(p.region_end_count, np.float64).sum())
        d = float(p.qc.deposited)
        gates.append(
            {
                "gate": f"gate ledger-closes[{name}]",
                "ok": abs(s - d) < 0.5 and abs(e - d) < 0.5,
                "detail": f"S={s:.0f} E={e:.0f} deposited={d:.0f}",
            }
        )

    # ── gate: the partitions sum to the whole, on BOTH new banks
    s_sum = sum(S.values())
    e_sum = sum(E.values())
    gates.append(
        {
            "gate": "gate partitions-sum-to-full[start]",
            "ok": bool(np.array_equal(s_sum, S_all)),
            "detail": f"max|delta| = {float(np.max(np.abs(s_sum - S_all))):.0f}",
        }
    )
    gates.append(
        {
            "gate": "gate partitions-sum-to-full[end]",
            "ok": bool(np.array_equal(e_sum, E_all)),
            "detail": f"max|delta| = {float(np.max(np.abs(e_sum - E_all))):.0f}",
        }
    )

    # ── the mask and the measured total
    w_max = w_max_from_deposited_lengths(payload.deposited_lengths)
    mask = build_region_wall_mask(ra, mature, reach_lo, reach_hi, w_max=w_max)
    total = build_total_abundance(chain, substrate, ra, geometry, mask, rna_pmf)

    # ── ⓔ THE START/END AGREEMENT ARM — the decisive WALL test, and the only arm that is FIELD-FREE
    #
    # ⭐⭐⭐ **THIS IS THE ARM THAT VALIDATES THE WALL RULE, AND IT WORKS UNDER CAPTURE.** S and E are
    # two estimators of the SAME per-region rate, both with opportunity `ell` for every fragment
    # length, both reading whatever field capture produced at that slot — so no uniformity assumption
    # enters and no length distribution is needed. The rule's content is exactly:
    #
    #   both sides exact  ->  S/ell and E/ell must AGREE          (ratio 1.000)
    #   start exact, end bound  ->  E is depressed, so E/S < 1
    #   end exact, start bound  ->  S is depressed, so E/S > 1
    #
    # ⛔ A mask that classified at random would give ratio ~1 in all three rows; a mask that has the
    # sides BACKWARDS inverts rows 2 and 3. So the arm has a direction, not just a magnitude.
    se: list[dict] = []
    s_ok_all = np.asarray(mask.start_exact, bool)
    e_ok_all = np.asarray(mask.end_exact, bool)
    live_r = (ell > 0) & ((S_all + E_all) > 0)
    for name, sel in (
        ("BOTH exact", s_ok_all & e_ok_all & live_r),
        ("start exact, end BOUND", s_ok_all & ~e_ok_all & live_r),
        ("end exact, start BOUND", e_ok_all & ~s_ok_all & live_r),
        ("DOUBLE-walled", ~s_ok_all & ~e_ok_all & live_r),
    ):
        if not sel.any():
            continue
        se.append(
            {
                "wall": name,
                "n_regions": int(sel.sum()),
                "start_mass": float(S_all[sel].sum()),
                "end_over_start": _ratio_of_sums(E_all, S_all, sel),
            }
        )

    # ── ⓐ the truncation law, per component, SPLIT BY WALL CLASS — UNIVERSAL
    #
    # ⭐⭐ The split is the decisive part, not a refinement. On a region whose START side is wall-EXACT
    # the law must read 1.000 for EVERY component; where the wall binds, S is depressed and the ratio
    # rises — and it rises by COMPONENT, because gDNA's template is the chromosome (no wall at all)
    # while a mature template ends at its own TES. A pooled ratio mixes the two and reads neither.
    trunc: list[dict] = []
    for k in parts:
        # ⛔⛔ **THE SUPPORT IS THE UNSPLICED LENGTH DISTRIBUTION, AND THE REASON IS A CONFOUND WORTH
        # NAMING**: a SPLICED fragment books a START (its first covered base lands in some region)
        # but can NEVER be CONTAINED — both endpoints of an annotated intron are region bounds, so
        # no spliced fragment touches the contained axis at all (`DESIGN.md` §3.1). Using the full
        # histogram therefore prices C's support with a population C cannot hold. ⚠ Removing it from
        # the DIVISOR is exact; the residue in the NUMERATOR `S` is not removable from the banks —
        # `region_start_count` has no spliced/unspliced split — so the RNA arms keep a confound
        # proportional to their spliced share, which is printed beside them. ⭐ gDNA CANNOT SPLICE
        # (AXIOM 0, and the oracle asserts it), so for gDNA the law is EXACTLY formable and the gDNA
        # rows are the decisive ones.
        unspl = _unspliced_hist(parts[k])
        pk = _pmf_cdf(unspl, ell.astype(np.int64))
        pred = (S[k] / np.where(ell > 0, ell, 1.0)) * pk
        spliced_share = _spliced_share(parts[k])
        for wall_name, wall_sel in (
            ("start-EXACT", np.asarray(mask.start_exact, bool)),
            ("start-BOUND", ~np.asarray(mask.start_exact, bool)),
        ):
            for i in range(len(ELL_BINS) - 1):
                lo, hi = ELL_BINS[i], ELL_BINS[i + 1]
                sel = (ell >= lo) & (ell < hi) & (S[k] > 0.0) & (pk > 0.0) & wall_sel
                if not sel.any():
                    continue
                trunc.append(
                    {
                        "component": k,
                        "wall": wall_name,
                        "ell_bin": f"[{lo},{hi if hi < (1 << 61) else 'inf'})",
                        "n_regions": int(sel.sum()),
                        "start_mass": float(S[k][sel].sum()),
                        "ratio": _ratio_of_sums(C[k], pred, sel),
                        "spliced_share": spliced_share,
                        "decisive": k == "gdna",
                    }
                )

    # ── THE REALISED FRAGMENT-LENGTH GAP, per condition, from the partitions' own deposit histograms.
    # ⭐⭐ On the fl-gap SIDE panel every conclusion hangs off this number, and "configured" is not
    # "realised": the sampler truncates to [frag_min, frag_max], and an mRNA fragment must additionally
    # fit inside its transcript while gDNA need not. So it is MEASURED here and printed on every row
    # (`gdna_ladder.yaml`'s header makes the same point about its equal arms).
    def _mu(hist):
        h = np.asarray(hist, np.float64)
        tot = h.sum()
        return float(np.arange(h.shape[0]) @ h / tot) if tot > 0 else float("nan")

    rna_hist = np.asarray(parts["mrna"].deposited_lengths, np.float64) + np.asarray(
        parts["nrna"].deposited_lengths, np.float64
    )
    fl = {
        "mu_gdna": _mu(parts["gdna"].deposited_lengths),
        "mu_rna": _mu(rna_hist),
        "mu_mrna": _mu(parts["mrna"].deposited_lengths),
        "mu_nrna": _mu(parts["nrna"].deposited_lengths),
    }
    fl["gap"] = fl["mu_rna"] - fl["mu_gdna"]

    # ── ⓕ THE CONSUMER POOLS — what would actually CHANGE if each consumer swapped estimators
    #
    # ⭐⭐⭐ **THIS IS THE ARM THAT PRICES A SWAP, AND IT NEEDS NO `src/` CHANGE.** Every consumer below
    # forms a pooled gDNA rate today as `Σcount / ΣE_contained` — unbiased, but the fragment-length pmf
    # enters the DIVISOR. The candidate replacement is `Σcounts / Σexposure` from
    # `total_abundance.region_counts_and_exposure`, where the exposure is the region's own length and no
    # pmf enters at all. ⛔ Both are scored against the SAME truth: the gDNA-only start rate over the
    # same population, taken from the `gdna` origin partition — which is untruncated and pmf-free, and
    # is therefore not the target of either estimator's own machinery.
    # ⚠ Each population's predicate is IMPORTED from the module that owns it, never restated here: a
    # second copy of a selection is how a consumer and its own audit drift apart.
    from rigel.calibration.signature import (
        BIT_EXON_NEG,
        BIT_EXON_POS,
        BIT_INTRON_NEG,
        BIT_INTRON_POS,
    )
    from rigel.calibration.density_model import count_observable_masks
    from rigel.calibration.signature import RegionType, coarse_type_array
    from rigel.calibration.total_abundance import region_counts_and_exposure

    sig_i = np.asarray(ra.signature, np.int64)
    _EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG
    _GENE_BITS = BIT_INTRON_POS | BIT_INTRON_NEG | _EXON_BITS
    del (
        _GENE_BITS
    )  # the intergenic pool is `coarse == INTERGENIC` below; kept for the exon mask only
    coarse = coarse_type_array(sig_i)
    # ⭐ The shipped divisor, obtained by PROJECTION exactly as `calibrate` does it, never recomputed
    # here: `contained_eff_length` was already applied at REGION slots by `build_region_geometry`, and a
    # second call would put two implementations of one quantity in the tree (the reason `_project_eff`
    # exists at all).
    eff_g_region = np.zeros(ell.shape[0], np.float64)
    eff_g_region[obj[is_region]] = np.asarray(geometry.eff_gdna, np.float64)[is_region]
    contained = np.asarray(substrate.region_contained.count, np.float64).sum(axis=1)
    new_counts, new_exposure, new_free = region_counts_and_exposure(substrate, ra, mask)
    obs_region, _obs_boundary = count_observable_masks(
        np.asarray(ra.signature), np.asarray(ra.ref_id)
    )

    pools = (
        (
            "fit_intron_background (→ psi, intron lambda-factor)",
            (coarse == int(RegionType.INTERGENIC)) & (eff_g_region > 0.0),
        ),
        (
            "density_model own branch (count-observable)",
            np.asarray(obs_region, bool) & (eff_g_region > 0.0),
        ),
    )
    # ⭐⭐ **THE EXON ROW IS WHERE THE DEFECT ACTUALLY LIVES, and it is scored on the TOTAL rather than
    # on a gDNA rate** — an exon's total is not a gDNA density (RNA is there by definition), so the
    # pure-gDNA truth above does not apply. What IS measurable is the shipped RECIPROCAL bank's
    # under-read of the total: `inv_abundance` at a REGION reads `rho*P(w<=ell)` and the exon median is
    # ~98 bp, so this row prints the truncation factor on the population that carries it.
    # ⛔ `new/truth` on this row is CIRCULAR BY CONSTRUCTION (both are S over ell, differing only by
    # the wall restriction) and is printed only to show the wall's cost — it validates nothing. The
    # non-circular number is `shipped/truth`.
    exon_sel = (coarse == int(RegionType.EXON)) & (ell > 0)
    exon_row = None
    if exon_sel.any():
        e_truth = _ratio_of_sums(S_all, ell, exon_sel)
        e_ship = _ratio_of_sums(
            np.asarray(substrate.region_contained.inv_opportunity_sum, np.float64),
            np.ones_like(ell),
            exon_sel,
        )
        sel_new_e = exon_sel & new_free
        exon_row = {
            "pool": "R exon — the TRUNCATION defect, on the TOTAL",
            "n_regions": int(exon_sel.sum()),
            "n_regions_model_free": int(sel_new_e.sum()),
            "coverage": float(ell[sel_new_e].sum() / ell[exon_sel].sum()),
            "truth_rho_total": e_truth,
            "shipped_reciprocal_over_truth": e_ship / e_truth if e_truth > 0 else float("nan"),
            "new_over_truth_CIRCULAR": (
                _ratio_of_sums(new_counts, new_exposure, sel_new_e) / e_truth
                if e_truth > 0
                else float("nan")
            ),
        }

    consumer: list[dict] = []
    for name, sel in pools:
        if not sel.any():
            continue
        truth = _ratio_of_sums(S["gdna"], ell, sel)
        shipped = _ratio_of_sums(contained, eff_g_region, sel)
        sel_new = sel & new_free
        new = _ratio_of_sums(new_counts, new_exposure, sel_new)
        consumer.append(
            {
                "pool": name,
                "n_regions": int(sel.sum()),
                "n_regions_model_free": int(sel_new.sum()),
                "coverage": float(
                    ell[sel_new].sum() / ell[sel].sum() if ell[sel].sum() > 0 else np.nan
                ),
                "truth_rho_gdna": truth,
                "shipped": shipped,
                "new_total": new,
                "shipped_over_truth": shipped / truth if truth > 0 else float("nan"),
                "new_over_truth": new / truth if truth > 0 else float("nan"),
            }
        )

    # ── ⓐ′ the same law with a PER-REGION-CLASS support — the capture diagnostic
    #
    # ⭐⭐ **WHY THIS ARM EXISTS: ⓐ's support is a GLOBAL histogram, and capture makes the fragment
    # length distribution vary SPATIALLY.** The contained bank is truncated by the LOCAL distribution,
    # so a global CDF prices it wrongly wherever capture has selected lengths — which is a comparator
    # limitation, not necessarily an estimator defect, and the two are separated here rather than
    # argued about. The accumulator already keeps the gDNA deposit histogram PER REGION CLASS
    # (`pool_lengths`: intergenic-contained and intron-contained are the two whose population is
    # exactly a contained gDNA fragment), so restricting to those two classes and using each one's
    # OWN pool is the same law with the right support and no new model.
    from rigel.scan_payload import POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC

    g_pools = np.asarray(parts["gdna"].pool_lengths, np.float64)
    per_class: list[dict] = []
    for cls_name, cls_id, pool in (
        ("R intergenic", int(RegionType.INTERGENIC), POOL_DNA_INTERGENIC),
        ("R intron", int(RegionType.INTRON), POOL_DNA_INTRONIC),
    ):
        pk_local = _pmf_cdf(g_pools[pool], ell.astype(np.int64))
        pred_local = (S["gdna"] / np.where(ell > 0, ell, 1.0)) * pk_local
        for wall_name, wall_sel in (
            ("start-EXACT", np.asarray(mask.start_exact, bool)),
            ("start-BOUND", ~np.asarray(mask.start_exact, bool)),
        ):
            sel = (coarse == cls_id) & (S["gdna"] > 0.0) & (pk_local > 0.0) & wall_sel
            if not sel.any():
                continue
            per_class.append(
                {
                    "region_class": cls_name,
                    "wall": wall_name,
                    "n_regions": int(sel.sum()),
                    "start_mass": float(S["gdna"][sel].sum()),
                    "ratio_local": _ratio_of_sums(C["gdna"], pred_local, sel),
                    "ratio_global": _ratio_of_sums(
                        C["gdna"],
                        (S["gdna"] / np.where(ell > 0, ell, 1.0))
                        * _pmf_cdf(_unspliced_hist(parts["gdna"]), ell.astype(np.int64)),
                        sel,
                    ),
                }
            )

    # ── ⓑ the uniform-field arm — gDNA ONLY, the absolute rate, split by wall class
    #
    # ⛔⛔ **gDNA IS THE ONLY COMPONENT WITH A UNIFORM-FIELD TRUTH, and pooling the three was the first
    # draft's error.** RNA density varies per transcript BY DESIGN (that is what abundance means), so
    # a per-reference pooled RNA rate is not that region's truth and scoring against it measures the
    # panel's expression profile, not the estimator. gDNA off capture is uniform — that is exactly
    # what `calibration_oracle.py`'s gdna-field-uniformity gate certifies — so the gDNA start rate
    # has a truth that does NOT come from the region's own bank, and the comparison is not circular.
    n_refs = int(ref_id.max()) + 1 if ref_id.size else 0
    num = np.bincount(ref_id, weights=S["gdna"], minlength=n_refs)
    den = np.bincount(ref_id, weights=ell, minlength=n_refs)
    rho_g_true = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)[ref_id]

    ok_ell = ell > 0
    g_meas = S["gdna"] / np.where(ok_ell, ell, 1.0)
    field = {}
    for wall_name, wall_sel in (
        ("start-EXACT", np.asarray(mask.start_exact, bool)),
        ("start-BOUND", ~np.asarray(mask.start_exact, bool)),
    ):
        sel = wall_sel & ok_ell & (rho_g_true > 0)
        w = ell[sel]
        err = float(np.sum(np.abs(g_meas[sel] - rho_g_true[sel]) * w))
        mass = float(np.sum(rho_g_true[sel] * w))
        # the SHIPPED contained bank, scored identically — the comparator this replaces
        c_err = float(np.sum(np.abs(C["gdna"][sel] - rho_g_true[sel]) * w))
        field[wall_name] = {
            "n_regions": int(sel.sum()),
            "mass": mass,
            "err": err,
            "rel": err / mass if mass > 0 else float("nan"),
            "shipped_err": c_err,
            "shipped_rel": c_err / mass if mass > 0 else float("nan"),
            "bias": float(np.sum((g_meas[sel] - rho_g_true[sel]) * w)),
        }

    live = total.model_free[is_region] & (ell[obj[is_region]] > 0)
    field_err = field["start-EXACT"]["err"]
    field_mass = field["start-EXACT"]["mass"]
    shipped_err = field["start-EXACT"]["shipped_err"]

    # ── ⓓ the boundary arm — no truncation, no wall
    b_total = total.total[is_boundary]
    b_ok = np.isfinite(b_total)

    # ── ⓒ the wall exposure
    exposure = {
        "n_regions": int(mask.n_regions),
        "w_max": int(w_max),
        "start_exact_frac_mass": _mass_frac(S_all, mask.start_exact),
        "end_exact_frac_mass": _mass_frac(S_all, mask.end_exact),
        "double_walled_frac_mass": _mass_frac(S_all, mask.double_walled),
        "both_exact_frac_mass": _mass_frac(S_all, mask.start_exact & mask.end_exact),
    }

    return {
        "condition": condition,
        "stratum": PVO.stratum(condition),
        "scope": OC._scope(condition),
        "zero_gdna": PVO.is_zero_gdna(condition),
        "gates": gates,
        "truncation": trunc,
        "per_class": per_class,
        "consumer": consumer,
        "fl": fl,
        "exon_row": exon_row,
        "start_end": se,
        "field": field,
        "field_err_fragments": field_err,
        "field_mass_fragments": field_mass,
        "field_rel": field_err / field_mass if field_mass > 0 else float("nan"),
        "shipped_err_fragments": shipped_err,
        "shipped_rel": shipped_err / field_mass if field_mass > 0 else float("nan"),
        "n_region_live": int(live.sum()),
        "n_boundary_finite": int(b_ok.sum()),
        "boundary_total_sum": float(np.nansum(b_total)),
        "exposure": exposure,
        "label": label,
    }


def _mass_frac(mass: np.ndarray, sel: np.ndarray) -> float:
    tot = float(np.asarray(mass, np.float64).sum())
    if tot <= 0.0:
        return float("nan")
    return float(np.asarray(mass, np.float64)[np.asarray(sel, bool)].sum() / tot)


def load_condition(suite: Path, condition: str, index, cached: dict) -> tuple[object, dict, dict]:
    """The full payload, the three origin payloads, and the index-derived parts (built once)."""
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    root = Path(suite) / "oracle_cache" / condition
    parts = {}
    for k in ORIGINS:
        d = root / k
        if not d.is_dir():
            raise FileNotFoundError(
                f"{d} is missing — run `panel.py cache` (or rescan_panels.py --conditions "
                f"{condition}) first; this instrument scores against the origin partitions only."
            )
        parts[k] = read_scan_cache(d, index).payload

    lift: dict = {}
    kw = calibration_inputs(cache, index, lift_out=lift)
    # ⭐ the DRAINED frame (the 2026-08-31 frame ruling): the audited total and the origin partitions
    # it is scored against are drained consistently (`lift_drain_parts` replays the whole's choices).
    payload = kw["payload"]
    from calibration._oracle import lift_drain_parts

    parts = dict(zip(ORIGINS, lift_drain_parts(lift, [parts[k] for k in ORIGINS])[0]))
    ra = cached["region_arrays"]
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    substrate = CalibrationSubstrate.from_payload(payload, ra)
    geometry = build_region_geometry(
        chain, substrate, ra, cached["sj"], kw["gdna_fl_pmf"], kw["rna_fl_pmf"], None
    )
    statics = build_region_statics(chain, ra, cached["bflags"])
    label = OC.strata(chain, statics, geometry, ra)["label"].astype(str)
    return (
        payload,
        parts,
        {
            "region_arrays": ra,
            "chain": chain,
            "substrate": substrate,
            "geometry": geometry,
            "statics": statics,
            "mature": cached["mature"],
            "reach": cached["reach"],
            "label": label,
            "rna_pmf": kw["rna_fl_pmf"],
        },
    )


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------


def report(rows: list[dict]) -> int:
    bad = 0
    for d in rows:
        ss, cap = d["stratum"]
        print(f"\n── {d['condition']}  [{ss} x {cap}]  {d['scope']}")
        for g in d["gates"]:
            flag = "✔" if g["ok"] else "⛔"
            if not g["ok"]:
                bad += 1
            print(f"   {flag} {g['gate']}: {g['detail']}")
        ex = d["exposure"]
        print(
            f"   wall: w_max {ex['w_max']}  start-exact {ex['start_exact_frac_mass']:.4f}  "
            f"end-exact {ex['end_exact_frac_mass']:.4f}  both {ex['both_exact_frac_mass']:.4f}  "
            f"⛔ double-walled {ex['double_walled_frac_mass']:.4f}  (START-mass weighted)"
        )
        vac = cap == "capture ON"
        stamp = (
            "  ⛔ VACUOUS here (capture makes the gDNA field non-uniform BY DESIGN)" if vac else ""
        )
        zc = "  ⛔ g00: no gDNA at all, so this arm has no population" if d["zero_gdna"] else ""
        print(f"   ⓑ uniform-field, gDNA ONLY — |S/ell − rho_g|·ell, in fragments{stamp}{zc}")
        for w, f in d["field"].items():
            print(
                f"      {w:>11}  n={f['n_regions']:>6}  mass={f['mass']:>13,.0f}  "
                f"measured {f['err']:>12,.0f} ({f['rel']:.4f})  signed {f['bias']:>+13,.0f}  |  "
                f"shipped contained {f['shipped_err']:>12,.0f} ({f['shipped_rel']:.4f})"
            )
        f = d["fl"]
        print(
            f"   fl REALISED: gDNA {f['mu_gdna']:.2f}  RNA {f['mu_rna']:.2f} "
            f"(mature {f['mu_mrna']:.2f} / nascent {f['mu_nrna']:.2f})  GAP {f['gap']:+.2f} bp"
        )
        print(
            "   ⓕ CONSUMER POOLS — pooled gDNA rate: SHIPPED (count/E_contained, pmf in the divisor) "
            "vs NEW (S or E over ell, pmf-free), both against the gdna partition's own start rate"
        )
        for c in d["consumer"]:
            print(
                f"      {c['pool']:<48} n={c['n_regions']:>6} (model-free {c['n_regions_model_free']:>6}, "
                f"cov {c['coverage']:.3f})"
            )
            print(
                f"         truth {c['truth_rho_gdna']:.6g}   shipped/truth {c['shipped_over_truth']:.4f}"
                f"   new/truth {c['new_over_truth']:.4f}"
            )
        if d["exon_row"]:
            x = d["exon_row"]
            print(
                f"      {x['pool']:<48} n={x['n_regions']:>6} (model-free {x['n_regions_model_free']:>6}, "
                f"cov {x['coverage']:.3f})"
            )
            print(
                f"         truth(total) {x['truth_rho_total']:.6g}   ⛔ shipped RECIPROCAL bank / truth "
                f"{x['shipped_reciprocal_over_truth']:.4f}  <- the truncation, on the population that "
                f"carries it   (new/truth {x['new_over_truth_CIRCULAR']:.4f} is CIRCULAR)"
            )
        print("   ⓔ START/END agreement — FIELD-FREE, capture-valid: Σ E / Σ S by wall class")
        for x in d["start_end"]:
            print(
                f"      {x['wall']:>24}  n={x['n_regions']:>6}  mass={x['start_mass']:>13,.0f}  "
                f"E/S = {x['end_over_start']:.4f}"
            )
        print(
            "   ⓐ′ the SAME law with each region class's OWN gDNA length pool — separates a"
            " comparator limitation under capture from an estimator defect"
        )
        for pc in d["per_class"]:
            print(
                f"      {pc['region_class']:>13} {pc['wall']:>11}  n={pc['n_regions']:>6}  "
                f"mass={pc['start_mass']:>13,.0f}  local-support {pc['ratio_local']:.4f}  "
                f"vs global-support {pc['ratio_global']:.4f}"
            )
        print(
            "   ⓐ truncation law  ratio = Σ Cinv / Σ (S/ell)·P_unspliced(w≤ell)   [1.000 is the "
            "claim; ⭐ = DECISIVE (gDNA cannot splice, so the law is exactly formable)]"
        )
        for t in d["truncation"]:
            star = "⭐" if t["decisive"] else "  "
            print(
                f"    {star} {t['component']:>5} {t['wall']:>11}  ell {t['ell_bin']:>12}  "
                f"n={t['n_regions']:>5}  mass={t['start_mass']:>12,.0f}  ratio={t['ratio']:.4f}"
                f"  (spliced share of this partition {t['spliced_share']:.3f})"
            )
    return bad


def summarise(rows: list[dict]) -> None:
    print("\n" + "=" * 100)
    print(
        "SUMMARY — per stratum, ratio of SUMS, the zero control on its own row, DEFERRED reported"
    )
    print("=" * 100)
    for label, pred in OC._SELECTIONS if hasattr(OC, "_SELECTIONS") else []:
        if label is None:
            print("   " + "-" * 90)
            continue
        sel = [d for d in rows if pred(d["condition"])]
        if not sel:
            continue
        num = sum(d["field_err_fragments"] for d in sel)
        den = sum(d["field_mass_fragments"] for d in sel)
        snum = sum(d["shipped_err_fragments"] for d in sel)
        dw = np.mean([d["exposure"]["double_walled_frac_mass"] for d in sel])
        print(
            f"   {label:<42} n={len(sel):>2}  measured {num / den if den else float('nan'):.4f}  "
            f"shipped {snum / den if den else float('nan'):.4f}  double-walled {dw:.4f}"
        )
    print(
        "\n⛔ ⓑ is VACUOUS on every capture-ON row (the field is non-uniform by design); read ⓐ there."
    )


# ---------------------------------------------------------------------------
# self-test — perturbed, no I/O
# ---------------------------------------------------------------------------


def self_test() -> int:
    import dataclasses

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests" / "calibration"))
    from _synthetic import delta_pmf, make_synthetic_payload, make_synthetic_sj

    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"⛔ {name}")

    base, ra = make_synthetic_payload()
    chain = build_region_chain(base.ref_region_offsets, base.ref_boundary_offsets)
    ell = np.asarray(ra.region_size_bp, np.float64)
    n_r = ell.shape[0]
    z = np.zeros((n_r, 2))
    mature = MatureWallDistances(d_low=z.copy(), d_high=z.copy(), covered=np.zeros((n_r, 2), bool))
    #: ⭐ GENEROUS reach on purpose: with a zero reach every RNA-admitting region is double-walled and
    #: the exposure arm has nothing to report. Here only the CONTIG walls bind, which is the shape the
    #: gDNA-only case has on a real reference.
    reach = (np.full((2, 2), 5000.0), np.full((2, 2), 5000.0))

    def condition_of(start, end, *, contained=None, pmf_at=60):
        """A consistent synthetic condition: banks stated, the contained bank derived to satisfy the
        truncation law exactly at ``P(w <= ell) == 1``, and the three partitions split so they sum to
        the whole and each is internally consistent."""
        s2 = np.asarray(start, np.uint32)
        e2 = np.asarray(end, np.uint32)
        dep = int(s2.sum())
        cinv = (
            s2.astype(np.float64).sum(1) / ell
            if contained is None
            else np.asarray(contained, float)
        )
        p = dataclasses.replace(
            base,
            region_start_count=s2,
            region_end_count=e2,
            region_contained_inv_opportunity_sum=cinv,
            deposited_lengths=_hist({pmf_at: dep}, 201),
            qc=dataclasses.replace(base.qc, deposited=dep),
        )
        return p, _split_three(p, pmf_at=pmf_at)

    def parts_for(payload):
        substrate = CalibrationSubstrate.from_payload(payload, ra)
        geometry = build_region_geometry(
            chain, substrate, ra, make_synthetic_sj(), delta_pmf(50, 200), delta_pmf(80, 200)
        )
        return {
            "region_arrays": ra,
            "chain": chain,
            "substrate": substrate,
            "geometry": geometry,
            "mature": mature,
            "reach": reach,
            "label": np.array(["R exon"] * int(chain.n_slots)),
            "rna_pmf": delta_pmf(80, 200),
        }

    COND = "gdna_g50_ss_0.99_nrna_mid_capture_off"
    #: a UNIFORM field: every region the same rate, so the arm's truth and its measurement coincide.
    uniform = [[6, 6], [6, 6], [6, 6]]
    p_u, parts_u = condition_of(uniform, uniform)
    ip_u = parts_for(p_u)
    d = measure(p_u, parts_u, ip_u, COND)
    check("every gate passes on a consistent synthetic condition", all(g["ok"] for g in d["gates"]))
    check("the stratum parses", d["stratum"] == ("stranded", "capture OFF"))
    check("the scope is stamped", d["scope"] == "IN SCOPE")

    # ⛔ a partition that does not sum to the whole must FAIL the gate, not be absorbed
    broken = dict(parts_u)
    broken["gdna"] = dataclasses.replace(
        parts_u["gdna"],
        region_start_count=np.asarray(parts_u["gdna"].region_start_count) + np.uint32(1),
    )
    d2 = measure(p_u, broken, ip_u, COND)
    check(
        "a partition that oversums FAILS sum-to-full",
        not next(g for g in d2["gates"] if "sum-to-full[start]" in g["gate"])["ok"],
    )

    # ⛔ a payload whose ledger does not close must fail its own gate
    d3 = measure(
        dataclasses.replace(
            p_u, region_end_count=np.array([[6, 6], [6, 6], [6, 7]], dtype=np.uint32)
        ),
        parts_u,
        ip_u,
        COND,
    )
    check(
        "a broken ledger FAILS ledger-closes[full]",
        not next(g for g in d3["gates"] if "ledger-closes[full]" in g["gate"])["ok"],
    )

    # ⭐ the truncation law: a DELTA pmf at 60 with ell = 100 gives P(w <= ell) = 1, so the contained
    # bank and S/ell must agree exactly — the law with its factor equal to one.
    ratios = [t["ratio"] for t in d["truncation"] if np.isfinite(t["ratio"])]
    check(
        "the truncation law reads 1.000 when P(w<=ell) == 1",
        bool(ratios) and max(abs(r - 1.0) for r in ratios) < 1e-9,
    )

    # ⛔ and it must MOVE when the contained bank is scaled — a law that cannot fail measures nothing
    s_u = np.asarray(uniform, np.float64).sum(1)
    p_h, parts_h = condition_of(uniform, uniform, contained=(s_u / ell) * 0.5)
    d5 = measure(p_h, parts_h, parts_for(p_h), COND)
    r5 = [t["ratio"] for t in d5["truncation"] if np.isfinite(t["ratio"])]
    check(
        "halving the contained bank halves the ratio",
        bool(r5) and max(abs(r - 0.5) for r in r5) < 1e-9,
    )

    # ⭐ the truncation law must also READ a real truncation: a pmf at 150 against ell = 100 gives
    # P(w <= ell) = 0, so the contained bank is 0 and the ratio must be 0 — the 11.6×-at-98 bp defect
    # in its extreme form.
    p_t, parts_t = condition_of(uniform, uniform, contained=np.zeros(n_r), pmf_at=150)
    d_t = measure(p_t, parts_t, parts_for(p_t), COND)
    check(
        "a pmf longer than every region reads ratio 0 (P(w<=ell) == 0)",
        all(abs(t["ratio"]) < 1e-12 for t in d_t["truncation"] if np.isfinite(t["ratio"]))
        or not d_t["truncation"],
    )

    # ⭐ the uniform-field arm: uniform by construction, so the error must be exactly 0
    check(
        "the uniform-field error is ~0 on a uniform synthetic field",
        d["field_err_fragments"] < 1e-9,
    )

    # ⛔ and it must be nonzero when the field is made lumpy
    lumpy = [[15, 15], [3, 0], [3, 0]]
    p_l, parts_l = condition_of(lumpy, lumpy)
    d6 = measure(p_l, parts_l, parts_for(p_l), COND)
    check("a lumpy field gives a nonzero uniform-field error", d6["field_err_fragments"] > 1.0)

    # ⓒ the exposure: only the contig walls bind, so the outer regions are one-sided and the middle
    # one is exact on both sides — strictly between 0 and 1, and nothing double-walled.
    check(
        "the wall exposure reports the outer regions as one-sided",
        0.0 < d["exposure"]["both_exact_frac_mass"] < 1.0,
    )
    check(
        "nothing is double-walled with a generous reach on a 300 bp reference",
        d["exposure"]["double_walled_frac_mass"] == 0.0,
    )
    # ⛔ ... and a longer w_max must move the verdict, by an amount stated in advance. At w_max = 200
    # the bar is 199: the INTERIOR region has 100 bases either side and goes double-walled, while the
    # two outer ones still clear on their one open side (200 >= 199). So exactly its third of the
    # START mass, and no more — an exposure that cannot move measures nothing.
    p_w, parts_w = condition_of(uniform, uniform, pmf_at=200)
    d7 = measure(p_w, parts_w, parts_for(p_w), COND)
    check(
        "a longer w_max double-walls exactly the interior region",
        abs(d7["exposure"]["double_walled_frac_mass"] - 1.0 / 3.0) < 1e-12,
    )

    # ⛔ the independent enumeration must agree with the kernel on a hand-built annotation
    check("the enumeration agrees with the kernel", _enumeration_agrees())

    # ⛔ ... and must DISAGREE when the kernel is fed a different annotation — a comparator that
    # cannot fail is not a comparator
    check("the enumeration comparator can fail", _enumeration_can_fail())

    print(f"self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def _hist(pairs: dict, size: int) -> np.ndarray:
    h = np.zeros(size, np.uint32)
    for k, v in pairs.items():
        h[k] = v
    return h


def _split_three(payload, *, pmf_at: int = 60):
    """Split a payload's banks three ways so the parts sum to the whole exactly, each part INTERNALLY
    consistent: the contained bank is split in PROPORTION to that part's own start bank.

    ⚠ A test fixture, not a truth: it exists so the sum-to-full gate has something consistent to pass
    and a deliberate over-sum has something to break. ⛔ The proportional split is the point — an
    equal-thirds contained bank against an unequal start split makes the truncation law read a
    fixture artefact rather than the estimator (measured: it read 0.92/1.00/1.08 and the first draft
    of this self-test failed on its own fixture).
    """
    import dataclasses

    out = {}
    banks = ("region_start_count", "region_end_count")
    s_full = np.asarray(payload.region_start_count, np.float64).sum(1)
    cinv_full = np.asarray(payload.region_contained_inv_opportunity_sum, np.float64)
    for i, k in enumerate(ORIGINS):
        kw = {}
        for b in banks:
            a = np.asarray(getattr(payload, b), np.int64)
            q, r = a // 3, a % 3
            kw[b] = (q + (r > i).astype(np.int64)).astype(np.uint32)
        share = np.divide(
            np.asarray(kw["region_start_count"], np.float64).sum(1),
            s_full,
            out=np.zeros_like(cinv_full),
            where=s_full > 0,
        )
        dep = int(np.asarray(kw["region_start_count"]).sum())
        kw["region_contained_inv_opportunity_sum"] = cinv_full * share
        kw["deposited_lengths"] = _hist({pmf_at: max(dep, 1)}, 201)
        kw["qc"] = dataclasses.replace(payload.qc, deposited=dep)
        out[k] = dataclasses.replace(payload, **kw)
    return out


def _hand_annotation():
    """A 5-region reference with three transcripts, and the wall answers by hand.

    bounds 0/100/200/300/400/500; T0 POS exons [0,200)+[300,400); T1 POS exon [0,100);
    T2 NEG exon [300,500).
    """
    import pandas as pd

    from rigel.calibration.region_arrays import RegionArrays

    bounds = np.array([0, 100, 200, 300, 400, 500], np.int64)
    frame = pd.DataFrame(
        {
            "region_id": np.arange(5, dtype=np.int64),
            "ref_name": pd.array(["c"] * 5, dtype="string"),
            "start": bounds[:-1],
            "end": bounds[1:],
            "length": bounds[1:] - bounds[:-1],
            "signature": np.zeros(5, np.uint8),
        }
    )
    ra = RegionArrays.from_frame(frame, {"c": 0})
    exons = {0: [(0, 200), (300, 400)], 1: [(0, 100)], 2: [(300, 500)]}
    strand = {0: int(Strand.POS), 1: int(Strand.POS), 2: int(Strand.NEG)}
    return ra, exons, strand


class _FakeIndex:
    def __init__(self, exons, strand):
        import pandas as pd

        self._exons = exons
        self.ref_name_to_id = {"c": 0}
        self.num_transcripts = len(exons)
        self.t_df = pd.DataFrame(
            {
                "t_index": np.array(sorted(exons), np.int64),
                "ref": ["c"] * len(exons),
                "strand": np.array([strand[t] for t in sorted(exons)], np.int64),
            }
        )

    def get_exon_intervals(self, t):
        e = self._exons.get(int(t))
        return None if e is None else np.asarray(e, np.int64)


def _enumeration_agrees() -> bool:
    from rigel.calibration.splice_graph import mature_wall_distances_kernel

    ra, exons, strand = _hand_annotation()
    slow = enumerate_mature_walls(_FakeIndex(exons, strand), ra)
    t, r, a, b = [], [], [], []
    for tid, blocks in exons.items():
        for s, e in blocks:
            t.append(tid)
            r.append(0)
            a.append(s)
            b.append(e)
    sof = np.zeros(max(exons) + 1, np.int64)
    for tid, s in strand.items():
        sof[tid] = s
    fast = mature_wall_distances_kernel(
        np.array(t, np.int64),
        np.array(r, np.int64),
        np.array(a, np.int64),
        np.array(b, np.int64),
        sof,
        ra,
    )
    return (
        np.array_equal(slow.d_low, fast.d_low)
        and np.array_equal(slow.d_high, fast.d_high)
        and np.array_equal(slow.covered, fast.covered)
    )


def _enumeration_can_fail() -> bool:
    from rigel.calibration.splice_graph import mature_wall_distances_kernel

    ra, exons, strand = _hand_annotation()
    slow = enumerate_mature_walls(_FakeIndex(exons, strand), ra)
    # drop T0's second exon from the kernel's input: the collapse must change
    t, r, a, b = [], [], [], []
    for tid, blocks in exons.items():
        for s, e in blocks:
            if tid == 0 and s == 300:
                continue
            t.append(tid)
            r.append(0)
            a.append(s)
            b.append(e)
    sof = np.zeros(max(exons) + 1, np.int64)
    for tid, s in strand.items():
        sof[tid] = s
    fast = mature_wall_distances_kernel(
        np.array(t, np.int64),
        np.array(r, np.int64),
        np.array(a, np.int64),
        np.array(b, np.int64),
        sof,
        ra,
    )
    return not np.array_equal(slow.d_high, fast.d_high)


# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None, help="default: every scan_cache entry")
    ap.add_argument("--out", type=Path, default=None, help="TSV, one row per condition")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args(argv)

    if args.self_test:
        return self_test()

    index = TranscriptIndex.load(args.index)
    ra = RegionArrays.from_index(index)
    cached = {
        "region_arrays": ra,
        "sj": build_sj_geometry_arrays(index),
        "bflags": build_boundary_flags_array(index),
        "reach": build_contiguous_boundary_reach_arrays(index),
        "mature": build_mature_wall_distances(index, ra),
    }
    print("⭐ checking the wall mask's mature distances against the independent enumeration …")
    slow = enumerate_mature_walls(index, ra)
    fast = cached["mature"]
    agree = (
        np.array_equal(slow.d_low, fast.d_low)
        and np.array_equal(slow.d_high, fast.d_high)
        and np.array_equal(slow.covered, fast.covered)
    )
    if not agree:
        n = int(np.sum(slow.d_high != fast.d_high) + np.sum(slow.d_low != fast.d_low))
        print(
            f"⛔ gate mask-vs-enumeration FAILED: {n} distance(s) differ between the vectorised "
            "kernel and the per-transcript enumeration. The wall mask is not validated."
        )
        return 1
    print(
        f"   ✔ gate mask-vs-enumeration: {int(fast.covered.any(axis=1).sum()):,} mature-covered "
        f"regions, every distance identical in two implementations"
    )

    conds = (
        [args.condition]
        if args.condition
        else sorted(p.name for p in (args.suite / "scan_cache").iterdir() if p.is_dir())
    )
    rows = []
    for c in conds:
        payload, parts, ip = load_condition(args.suite, c, index, cached)
        rows.append(measure(payload, parts, ip, c))
    bad = report(rows)
    summarise(rows)

    if args.out:
        # ⭐ ONE ROW PER (condition x consumer pool) — the shape arm ⓕ's question has. The per-condition
        # scalars (the wall exposure, the realised fl gap) are repeated on every row of that condition so
        # the file needs no join, and the fl columns are here because on the fl-gap SIDE panel every
        # conclusion hangs off the REALISED gap rather than the configured one.
        cols = [
            "condition",
            "strand",
            "capture",
            "scope",
            "mu_gdna",
            "mu_rna",
            "fl_gap",
            "w_max",
            "double_walled_frac_mass",
            "both_exact_frac_mass",
            "pool",
            "n_regions",
            "n_regions_model_free",
            "coverage",
            "truth",
            "shipped_over_truth",
            "new_over_truth",
        ]
        with open(args.out, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for d in rows:
                ss, cap = d["stratum"]
                f = d["fl"]
                ex = d["exposure"]
                head = [
                    d["condition"],
                    ss,
                    cap,
                    d["scope"],
                    f"{f['mu_gdna']:.4f}",
                    f"{f['mu_rna']:.4f}",
                    f"{f['gap']:.4f}",
                    ex["w_max"],
                    ex["double_walled_frac_mass"],
                    ex["both_exact_frac_mass"],
                ]
                for c in d["consumer"]:
                    fh.write(
                        "\t".join(
                            str(x)
                            for x in head
                            + [
                                c["pool"],
                                c["n_regions"],
                                c["n_regions_model_free"],
                                c["coverage"],
                                c["truth_rho_gdna"],
                                c["shipped_over_truth"],
                                c["new_over_truth"],
                            ]
                        )
                        + "\n"
                    )
                # ⚠ the exon row is scored on the TOTAL, not on a gDNA rate, so its `truth` column is a
                # different quantity and its `new_over_truth` is circular by construction — both are
                # labelled in the pool name rather than silently sharing the numeric columns.
                x = d["exon_row"]
                if x:
                    fh.write(
                        "\t".join(
                            str(v)
                            for v in head
                            + [
                                x["pool"],
                                x["n_regions"],
                                x["n_regions_model_free"],
                                x["coverage"],
                                x["truth_rho_total"],
                                x["shipped_reciprocal_over_truth"],
                                x["new_over_truth_CIRCULAR"],
                            ]
                        )
                        + "\n"
                    )
        print(f"\nwrote {args.out}")
    if bad:
        print(f"\n⛔ {bad} gate failure(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
