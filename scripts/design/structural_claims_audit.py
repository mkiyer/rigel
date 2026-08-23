"""IS EVERY SLOT THE STAGE-0 SUBSTRATE ADMITS TRULY WHAT IT CLAIMS? — a confusion matrix against
certified slot truth, no solver.

The stage-0 substrate (`rigel.calibration.structural_claims`) admits slots on STRUCTURE alone, and
every class carries a claim the certified truth can test directly, in FRAGMENTS:

    intergenic           n_nrna + n_mrna == 0   no RNA strand is admissible at all
    ss_intron_region     n_mrna == 0            contiguous mature RNA does not fit inside
    ss_intron_boundary   n_mrna == 0            an unspliced crossing has no mature term
    solvable_exon        n_mrna == 0            **at the LICENSING FLANK(S)** — the claim is the
                                                flank's, never the exon's own mass

Per condition it rebuilds the chain and statics from the index and scan cache (exactly as
`calibration_oracle.py` does), derives the claims, reads ``slot_truth.npz`` beside the oracle cache
(⛔ REFUSED if absent — run `calibration_oracle.py` first; a merely plausible truth is how a predicate
bug survives), and reports per class: slots admitted, live slots, claimed fragments, violating slots,
violating fragments, and the worst offenders. A single violating fragment falsifies the predicate for
that class — that is the point of deriving the substrate from structure.

⛔ JUDGE ONLY THE CLAIMED SLOTS (owner): the un-admitted remainder appears once, as COVERAGE context,
and is never scored. ⚠ Nascent fragments inside an ss_intron slot are NOT a violation — they are the
population stage 2 exists to deconvolve; only a MATURE fragment there is.

Exit status: nonzero iff any claim is violated anywhere. ``--self-test`` runs the checker against
synthetic truth with injected violations and must catch every one of them (and must NOT score a
violation parked on an unclaimed slot).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import RegionStatics, build_region_statics  # noqa: E402
from rigel.calibration.splice_graph import build_boundary_flags_array  # noqa: E402
from rigel.calibration.structural_claims import build_structural_claims  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import read_scan_cache  # noqa: E402

#: (class name, forbidden-population keys) — the claim each class makes, in truth-table columns.
CLAIMS = (
    ("intergenic", ("n_nrna", "n_mrna")),
    ("ss_intron_region", ("n_mrna",)),
    ("ss_intron_boundary", ("n_mrna",)),
    ("solvable_exon", ("n_mrna",)),
)


def _licensing_flanks(chain, claims) -> np.ndarray:
    """The UNIQUE slot ids of every flank that licenses a solvable exon — where that class's claim
    actually lives. Unique, because two exons may be licensed through one boundary and a fragment
    must not be counted twice."""
    left = np.asarray(chain.left, np.int64)[np.asarray(claims.exon_flank_left, bool)]
    right = np.asarray(chain.right, np.int64)[np.asarray(claims.exon_flank_right, bool)]
    return np.unique(np.concatenate([left, right]))


def confusion(chain, claims, truth: dict) -> list[dict]:
    """The per-class rows: each class's admitted slots scored against ITS OWN claim, plus one
    unscored coverage row. ``truth`` needs ``count``/``n_nrna``/``n_mrna`` slot-keyed arrays."""
    count = np.asarray(truth["count"], np.float64)
    rows = []
    for name, forbidden in CLAIMS:
        label = name
        if name == "solvable_exon":
            idx = _licensing_flanks(chain, claims)
            sel = np.zeros(int(chain.n_slots), bool)
            sel[idx] = True
            n_exons = int(np.asarray(claims.solvable_exon, bool).sum())
            label = f"solvable_exon ({n_exons} exons, at their flanks)"
        else:
            sel = np.asarray(getattr(claims, name), bool)
        bad_mass = np.zeros(int(chain.n_slots), np.float64)
        for key in forbidden:
            bad_mass += np.asarray(truth[key], np.float64)
        bad_mass = np.where(sel, bad_mass, 0.0)
        violating = bad_mass > 0
        order = np.argsort(bad_mass)[::-1]
        worst = [(int(s), float(bad_mass[s])) for s in order[:3] if bad_mass[s] > 0]
        rows.append(
            {
                "class": name,
                "label": label,
                "n_slots": int(sel.sum()),
                "n_live": int((sel & (count > 0)).sum()),
                "claimed_fragments": float(count[sel].sum()),
                "violating_slots": int(violating.sum()),
                "violating_fragments": float(bad_mass.sum()),
                "worst": worst,
            }
        )
    claimed = np.asarray(claims.claimed, bool)
    lib = float(count.sum())
    rows.append(
        {
            "class": "(coverage — context only, never scored)",
            "n_slots": int(claimed.sum()),
            "n_live": int((claimed & (count > 0)).sum()),
            "claimed_fragments": float(count[claimed].sum()),
            "violating_slots": None,
            "violating_fragments": None,
            "worst": [],
            "library_fragments": lib,
        }
    )
    return rows


def audit(index, region_arrays, suite: Path, condition: str) -> list[dict]:
    """One condition: rebuild chain + statics, derive the claims, score against certified truth."""
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    chain = build_region_chain(
        cache.payload.ref_region_offsets, cache.payload.ref_boundary_offsets
    )
    statics = build_region_statics(chain, region_arrays, build_boundary_flags_array(index))
    claims = build_structural_claims(chain, statics)

    npz = Path(suite) / "oracle_cache" / condition / "slot_truth.npz"
    if not npz.exists():
        raise FileNotFoundError(
            f"{npz} is missing — the certified truth is a precondition, not an option. "
            f"Build it: python scripts/design/calibration_oracle.py --suite {suite} "
            f"--condition {condition}"
        )
    truth = dict(np.load(npz))
    if truth["kind"].shape[0] != chain.n_slots or not np.array_equal(
        truth["kind"], np.asarray(chain.kind)
    ):
        raise ValueError(
            f"{npz} is not slot-aligned with the chain rebuilt from this index + scan cache "
            f"({truth['kind'].shape[0]} truth slots vs {chain.n_slots}). The oracle table predates a "
            f"rebuild — regenerate it with calibration_oracle.py before trusting anything here."
        )
    return confusion(chain, claims, truth)


def report(condition: str, rows: list[dict]) -> int:
    print(f"\n== {condition}")
    print(
        f"   {'class':<42} {'slots':>8} {'live':>8} {'fragments':>14} {'bad slots':>10} {'bad frags':>12}"
    )
    n_bad = 0
    for r in rows:
        vs = "—" if r["violating_slots"] is None else str(r["violating_slots"])
        vf = "—" if r["violating_fragments"] is None else f"{r['violating_fragments']:.0f}"
        print(
            f"   {r.get('label', r['class']):<42} {r['n_slots']:>8} {r['n_live']:>8} "
            f"{r['claimed_fragments']:>14.0f} {vs:>10} {vf:>12}"
        )
        if r["violating_fragments"]:
            n_bad += 1
            for slot, mass in r["worst"]:
                print(f"      ⛔ worst offender: slot {slot} carries {mass:.0f} forbidden fragments")
        if "library_fragments" in r and r["library_fragments"] > 0:
            share = r["claimed_fragments"] / r["library_fragments"]
            print(
                f"      coverage: {share:.1%} of the library's {r['library_fragments']:.0f} "
                f"unspliced/contained fragments sit on claimed slots"
            )
    return n_bad


# ── self-test: falsify the CHECKER — injected violations must be caught, and only on claimed slots ──


def _synthetic() -> tuple:
    """One reference, five regions: intergenic · exon · intron · exon · intergenic (a two-exon
    transcript), built directly as chain + statics with a DONOR flag at the exon|intron boundary."""
    from rigel.calibration.splice_graph import FLAG_DONOR_POS, FLAG_TSS_POS

    chain = build_region_chain(np.array([0, 5]), np.array([0, 4]))
    is_region = np.asarray(chain.kind) == REGION
    # slots 0..8: N(ig) B(edge) N(exon) B(donor) N(intron) B(acceptor) N(exon) B(edge) N(ig).
    # the transcript spans its own intron, so continuity holds at slots 2..6; the gene edges do not.
    fp = np.array([0, 0, 1, 1, 1, 1, 1, 0, 0], bool)
    fn = np.zeros(9, bool)
    mp = np.array([0, 0, 1, 0, 0, 0, 1, 0, 0], bool)
    mn = np.zeros(9, bool)
    bflags = np.zeros(9, np.uint16)
    bflags[1] = FLAG_TSS_POS  # gene edge — a terminus, licenses nothing
    bflags[3] = FLAG_DONOR_POS  # exon|intron — licenses the exon at slot 2
    bflags[5] = FLAG_DONOR_POS  # intron|exon — licenses the exon at slot 6
    bflags[7] = FLAG_TSS_POS
    statics = RegionStatics(
        n_slots=9,
        free_pos=fp,
        free_neg=fn,
        mrna_active_pos=mp,
        mrna_active_neg=mn,
        boundary_flags=np.where(~is_region, bflags, 0).astype(np.uint16),
    )
    claims = build_structural_claims(chain, statics)
    truth = {
        "count": np.full(9, 10.0),
        "n_nrna": np.zeros(9),
        "n_mrna": np.zeros(9),
    }
    truth["n_mrna"][2] = truth["n_mrna"][6] = 10.0  # the exons' own mature mass — NOT a violation
    truth["count"][0] = truth["count"][8] = 3.0
    return chain, claims, truth


def self_test() -> int:
    ok = 0

    def check(name: str, cond: bool):
        nonlocal ok
        print(f"   {'✔' if cond else '✘'} {name}")
        if not cond:
            raise SystemExit(f"self-test FAILED at: {name}")
        ok += 1

    chain, claims, truth = _synthetic()
    rows = {r["class"]: r for r in confusion(chain, claims, truth)}
    check(
        "the designed substrate is admitted (intergenic 4, intron 1+2, both exons)",
        rows["intergenic"]["n_slots"] == 4
        and rows["ss_intron_region"]["n_slots"] == 1
        and rows["ss_intron_boundary"]["n_slots"] == 2
        and rows["solvable_exon"]["n_slots"] == 2,
    )
    check(
        "a clean truth scores ZERO violations on every class",
        all(not rows[n]["violating_fragments"] for n, _ in CLAIMS),
    )
    check(
        "mature mass on the exon's OWN slot is not scored — the claim is the flank's",
        rows["solvable_exon"]["violating_fragments"] == 0.0,
    )

    t = {k: v.copy() for k, v in truth.items()}
    t["n_mrna"][4] = 7.0  # mature parked on the INTRON region slot
    r = {x["class"]: x for x in confusion(chain, claims, t)}
    check(
        "an injected mature fragment inside the intron REGION is caught, exactly",
        r["ss_intron_region"]["violating_fragments"] == 7.0
        and r["ss_intron_region"]["worst"] == [(4, 7.0)],
    )

    t = {k: v.copy() for k, v in truth.items()}
    t["n_mrna"][3] = 2.0  # mature crossing at the licensing DONOR flank
    r = {x["class"]: x for x in confusion(chain, claims, t)}
    check(
        "an injected mature crossing at a licensing flank is caught on BOTH its classes",
        r["solvable_exon"]["violating_fragments"] == 2.0
        and r["ss_intron_boundary"]["violating_fragments"] == 2.0,
    )

    t = {k: v.copy() for k, v in truth.items()}
    t["n_nrna"][0] = 5.0  # nascent parked on intergenic
    r = {x["class"]: x for x in confusion(chain, claims, t)}
    check(
        "injected nascent on an intergenic slot is caught (that class forbids ALL RNA)",
        r["intergenic"]["violating_fragments"] == 5.0,
    )

    t = {k: v.copy() for k, v in truth.items()}
    t["n_nrna"][4] = 9.0  # nascent inside the intron — stage 2's population, NOT a violation
    r = {x["class"]: x for x in confusion(chain, claims, t)}
    check(
        "nascent inside an ss intron is NOT a violation — it is what stage 2 deconvolves",
        r["ss_intron_region"]["violating_fragments"] == 0.0,
    )

    t = {k: v.copy() for k, v in truth.items()}
    t["n_mrna"][2] = 999.0  # a huge violation parked on an UNCLAIMED slot (the exon's own mass)
    r = {x["class"]: x for x in confusion(chain, claims, t)}
    check(
        "a violation parked on an unclaimed slot is scored NOWHERE (judge only what is claimed)",
        all(not r[n]["violating_fragments"] for n, _ in CLAIMS),
    )

    print(f"self-test: {ok}/{ok}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    index = TranscriptIndex.load(args.index)
    region_arrays = RegionArrays.from_index(index)
    conds = (
        [args.condition]
        if args.condition
        else sorted(p.name for p in (args.suite / "scan_cache").iterdir())
    )
    bad = 0
    for c in conds:
        bad += report(c, audit(index, region_arrays, args.suite, c))
    if bad:
        print(f"\n⛔ {bad} class×condition claims VIOLATED — the stage-0 predicate is falsified there.")
    else:
        print("\n⭐ every claim holds on every condition — the substrate is what it says it is.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
