#!/usr/bin/env python
"""⭐⭐⭐ IS A ZERO-COUNT ANCHOR'S DENSITY CLAIM TRUE OF ITS NEIGHBOURHOOD? — no solver runs.

⛔ **THE QUESTION.** A structurally pure-gDNA REGION (intergenic, ``struct_lock``) with **zero** counts now
forms and transmits the claim ``rho_g = 0`` at precision ``1/trigamma(1/2) = 0.2026`` (TRAPS: a-zero-count-is-a-measurement/TRAPS: a-ratio-cannot-carry-zero).
Off capture that is the strongest true statement in a gDNA-free library. **Under capture an empty
intergenic region may mean "no probe here" rather than "no gDNA here"** — and the relay's currency is a
density RATIO, so a source claiming zero multiplies every destination it reaches to zero.

This measures, from the origin-split oracle alone, whether that claim is true of the neighbourhood the
message reaches. Everything below is an identity on truth; nothing is fitted and no solve happens.

::

    rho_g_true[slot] = gdna_count[slot] / eff_gdna[slot]

⭐ ``eff_gdna`` is built by the SHIPPED `build_region_geometry` from the gDNA partition's OWN fitted length
model — i.e. the TRUE post-capture gDNA pmf — so the divisor is the solver's own function evaluated at the
exact pmf `SUCCESS.md` calls the ceiling. A census that wrote its own opportunity would be measuring a
different program.

⛔⛔ **THE POPULATION IS `g1_locked` ∧ REGION, AND IT IS *NOT* ``strand_evidence``'s ``struct_lock``.** This
docstring and the comment at the mask both said "exactly ``struct_lock``" until 2026-08-11 and the two
masks differ by **15–23×**:

============================================  ==========  =======================================
mask                                          size        what it is
============================================  ==========  =======================================
``anchor`` — what this instrument measures         1,312   ``g1_locked ∧ REGION``: RNA is inadmissible
                                                           on both strands, so the composition is
                                                           structurally certain
``struct_lock`` — what the solver actually    14,875 to    ``~((free_pos|free_neg) ∧ n_region>0) ∧
locks (`region_init.py`)                            30,423   REGION`` — which is ``g1_locked`` **∪ every
                                                           zero-count REGION**, empty exons and
                                                           introns included
============================================  ==========  =======================================

⭐ ``anchor`` is the mask the docstring of ``struct_lock`` *describes* ("scoped to true intergenic REGION
regions") and the one the standing xfail wants it rescoped to. The shipped code grants certainty to
**13,563–29,111 further slots** whose ``f_g`` is a default belief and whose evidence is nothing at all.

⛔ **So this instrument's 346× verdict is a statement about 1,312 anchors, not about the population the
solver locks.** It is not thereby wrong — the anchors are what the question is about — but its SCOPE was
overstated wherever it was quoted, and every run now prints both sizes so the reading cannot drift again.

⚠ Measured on four ladder conditions 2026-08-11: the ``∧ intergenic`` term is **redundant** —
``g1_locked ∧ REGION`` and ``g1_locked ∧ intergenic ∧ REGION`` are the same 1,312 slots on all four. It is
kept and ASSERTED rather than deleted, so an index that ever separates them says so instead of silently
widening this instrument's population.

**The three columns that answer it**, per condition, over that anchor population:

* ``empty%``     — share of anchors holding zero gDNA counts. These are the ones that claim ``rho_g = 0``.
* ``rho_nb``     — the mass-weighted TRUE gDNA density at the chain neighbours of the EMPTY anchors.
                   ⭐ If this is ~0 the claim is true and the anchor is the 39 % win. If it is far from 0
                   the anchor is asserting a zero into a gDNA-rich neighbourhood.
* ``lie x``      — ``rho_nb / rho_anchor_nonempty``: how many-fold the zero understates what a NON-empty
                   anchor in the same library would have said. A pure depletion artefact shows up here as
                   the probe-retention factor and nowhere else.

⚠ **The neighbour set is the chain's own adjacency** (`chain.left`/`chain.right`), not a re-derived
traversal: genomic order IS slot order and the chain strictly alternates REGION/BOUNDARY, so an anchor's
neighbours are ``i-1`` and ``i+1`` and a reference terminal links to ``-1``.

⛔ **This instrument cannot say the anchor caused the regression** — it says only whether the anchor's
claim is true. A claim being false is necessary for the hypothesis and not sufficient; the panel arm is.

Usage::

    python scripts/design/anchor_opportunity_census.py                 # the whole ladder
    python scripts/design/anchor_opportunity_census.py --conditions gdna_g75_ss_0.99_nrna_none_capture_on
"""

from __future__ import annotations

import argparse
import dataclasses
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import (  # noqa: E402
    build_region_geometry,
    build_region_statics,
    g1_locked,
)
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _native_detect_sj_tag  # noqa: E402
from rigel.scan_cache import read_scan_cache  # noqa: E402

DEFAULT_SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
DEFAULT_INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"
ORIGINS = ("gdna", "mrna", "nrna")
INTERGENIC = 0
_EPS = 1.0e-12


def _payloads(cache_root: Path, suite: Path, tag: str, index, cfg):
    """The three origin partitions, through the SHIPPED ``read_scan_cache`` — which refuses a payload
    whose ``graph_hash`` / ``reach_digest`` / ``payload_schema_digest`` does not describe this index, so
    a stale cache is refused loudly rather than feeding a truth number."""
    bam = str(suite / tag / "sim_oracle.bam")
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    return {k: read_scan_cache(cache_root / tag / k, index, scan).payload for k in ORIGINS}


def measure(parts, index, ra, sj, boundary_flags) -> dict:
    """Per-slot TRUE gDNA density, and the anchor population's own view of it."""
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.sj_opportunity import crossing_probability_from_index

    gp = parts["gdna"]
    chain = build_region_chain(gp.ref_region_offsets, gp.ref_boundary_offsets)
    statics = build_region_statics(chain, ra, boundary_flags)

    # ⭐ the divisor comes from the gDNA partition's OWN length model — the true post-capture gDNA pmf.
    size = int(gp.max_length)
    fl = build_fl_models(
        gp,
        sj_opportunity=crossing_probability_from_index(index, size),
        gdna_opportunity=gdna_opportunity_from_index(index, size),
    )
    geom = build_region_geometry(
        chain,
        CalibrationSubstrate.from_payload(gp, ra),
        ra,
        sj,
        fl.gdna_pmf,
        fl.rna_pmf,
        None,
    )
    eff_g = np.asarray(geom.eff_gdna, np.float64)

    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)

    def _counts(payload) -> np.ndarray:
        """One origin's unspliced count on every slot — the same populations ``eff_gdna`` divides."""
        sub = CalibrationSubstrate.from_payload(payload, ra)
        out = np.zeros(chain.n_slots, np.float64)
        out[kind == REGION] = np.asarray(sub.region_contained.count, np.float64).sum(1)[
            obj[kind == REGION]
        ]
        out[kind == BOUNDARY] = np.asarray(sub.boundary_unspliced.count, np.float64).sum(1)[
            obj[kind == BOUNDARY]
        ]
        return out

    n_g = _counts(parts["gdna"])
    n_r = _counts(parts["mrna"]) + _counts(parts["nrna"])
    live = eff_g > 0.0
    rho_g = np.where(live, n_g / np.where(live, eff_g, 1.0), 0.0)

    # ⛔ THE ANCHOR POPULATION IS ``g1_locked ∧ REGION`` — NOT ``strand_evidence``'s ``struct_lock``, which
    # is 15-23x larger. See the module docstring's table; both claimed "exactly struct_lock" until
    # 2026-08-11 and neither mask matched it.
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    is_region = kind == REGION
    ntype = np.where(is_region, rtype[np.clip(obj, 0, rtype.shape[0] - 1)], -1)
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    locked = g1_locked(fp, fn)
    anchor = is_region & (ntype == INTERGENIC) & locked

    # ⭐ the ``∧ intergenic`` term is measured REDUNDANT (2026-08-11, four conditions). Asserted rather
    # than deleted: if an index ever admits a G1 REGION that is not intergenic, this instrument's
    # population would silently narrow and its verdict would change meaning with no tell.
    if not np.array_equal(anchor, is_region & locked):
        raise AssertionError(
            f"`g1_locked & REGION` ({int((is_region & locked).sum()):,}) and `g1_locked & intergenic & "
            f"REGION` ({int(anchor.sum()):,}) have SEPARATED on this index. They were the same 1,312 "
            f"slots when this instrument's population was defined; re-read the docstring's table "
            f"before trusting any number below."
        )

    # the mask the SOLVER actually locks, so the scope of every verdict below is printed beside it
    # rather than asserted in prose (`region_init.py`: `locked = ~solvable`, `solvable = (fp|fn) & n>0`).
    mass = n_g + n_r
    struct_lock = (~((fp | fn) & (mass > 0.0))) & is_region

    left, right = np.asarray(chain.left, np.int64), np.asarray(chain.right, np.int64)
    empty = anchor & (n_g <= 0.0)
    nb = np.zeros(chain.n_slots, bool)
    for side in (left, right):
        idx = side[empty]
        nb[idx[idx >= 0]] = True
    nb &= ~anchor  # a neighbour that is itself an anchor says nothing new

    mw = lambda sel: (  # noqa: E731
        float(np.sum(mass[sel] * rho_g[sel]) / np.sum(mass[sel])) if np.sum(mass[sel]) > 0 else 0.0
    )
    full = anchor & (n_g > 0.0)
    return dict(
        n_anchor=int(anchor.sum()),
        # ⭐ the two population sizes, side by side. `n_lock_extra` is the count of slots the solver
        # treats as composition-CERTAIN and this instrument never looked at.
        n_lock=int(struct_lock.sum()),
        n_lock_extra=int((struct_lock & ~anchor).sum()),
        n_empty=int(empty.sum()),
        eff_empty=float(eff_g[empty].sum()),
        eff_anchor=float(eff_g[anchor].sum()),
        rho_full=mw(full),
        rho_nb=mw(nb),
        n_nb=int(nb.sum()),
        lib_gdna=float(n_g.sum()),
        lib_rna=float(n_r.sum()),
        # ⭐ the library-wide gDNA density, the honest reference for "what should an anchor say?"
        rho_lib=float(n_g.sum() / max(eff_g.sum(), _EPS)),
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default=str(DEFAULT_SUITE))
    ap.add_argument("--index", default=str(DEFAULT_INDEX))
    ap.add_argument("--cache", default=None, help="default: <suite>/oracle_cache")
    ap.add_argument("--conditions", nargs="*", default=None)
    args = ap.parse_args()

    suite = Path(args.suite)
    cache = Path(args.cache) if args.cache else suite / "oracle_cache"
    index = TranscriptIndex.load(args.index)
    cfg = PipelineConfig()
    ra = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    from rigel.calibration.splice_graph import (
        build_boundary_flags_array,
        build_sj_geometry_arrays,
    )

    sj = build_sj_geometry_arrays(index)
    boundary_flags = build_boundary_flags_array(index)
    conds = args.conditions or sorted(p.name for p in cache.iterdir() if p.is_dir())

    # ⭐ `anchor` is this instrument's population; `lock` is what the SOLVER locks. They are not the same
    # mask and the ratio is the scope of every verdict on the row (see the module docstring).
    print(
        f"{'condition':<42} {'anchor':>7} {'lock':>7} {'lock/anch':>9} {'empty':>7} {'empty%':>7} "
        f"{'Mb empty':>9} {'rho_full':>9} {'rho_nb':>9} {'rho_lib':>9} {'lie x':>7}"
    )
    print("-" * 136)
    rows = []
    for tag in conds:
        try:
            m = measure(_payloads(cache, suite, tag, index, cfg), index, ra, sj, boundary_flags)
        except Exception as e:  # noqa: BLE001
            print(f"{tag:<42} SKIP  {type(e).__name__}: {e}")
            continue
        pct = 100.0 * m["n_empty"] / max(m["n_anchor"], 1)
        lie = m["rho_nb"] / m["rho_full"] if m["rho_full"] > 0 else float("inf")
        print(
            f"{tag:<42} {m['n_anchor']:>7,} {m['n_lock']:>7,} "
            f"{m['n_lock'] / max(m['n_anchor'], 1):>8.1f}x {m['n_empty']:>7,} {pct:>6.1f}% "
            f"{m['eff_empty'] / 1e6:>9.2f} {m['rho_full']:>9.5f} {m['rho_nb']:>9.5f} "
            f"{m['rho_lib']:>9.5f} {lie:>7.1f}"
        )
        rows.append((tag, pct, m))

    if not rows:
        print("\nnothing measured — is the oracle cache present?")
        return 1

    print("\n" + "=" * 136)
    print(
        "⛔ SCOPE: `anchor` is `g1_locked & REGION`; `lock` is the solver's own `struct_lock`, which also\n"
        "   contains every ZERO-COUNT region (empty exons and introns). Every verdict on a row is a\n"
        "   statement about the `anchor` column only — `lock/anch` is the factor by which the solver\n"
        "   locks more than this instrument measured.\n"
        "⭐ rho_full = mass-weighted TRUE gDNA density at anchors that DO hold counts;\n"
        "   rho_nb   = the same at the chain neighbours of the EMPTY anchors — what the zero is asserted"
        " ONTO;\n"
        "   lie x    = rho_nb / rho_full. ⛔ A large value means the zero-count anchor's claim is false"
        " of\n"
        "              its own neighbourhood by that factor."
    )
    for axis in ("capture_off", "capture_on", "ss_0.50", "ss_0.99"):
        sub = [r for r in rows if axis in r[0]]
        if sub:
            print(
                f"   {axis:<12} mean empty% {np.mean([r[1] for r in sub]):>5.1f}   "
                f"mean lie x {np.mean([min(r[2]['rho_nb'] / r[2]['rho_full'], 1e6) if r[2]['rho_full'] > 0 else 0.0 for r in sub]):>8.1f}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
