"""HOW WELL DOES PASS-0 SOLVE THE SLOTS IT CLAIMS, PER POLICY? — silent / relay / fanout at the
stage-0 substrate's claimed populations, judged on NOTHING else (the owner's rule).

Two claimed populations, each scored as misplaced gDNA fragments ``Σ|est_gdna − true_gdna|`` against
certified `slot_truth`, per condition and per policy, DELIVER/REFUTE split and never pooled
(`TRAPS: a-refutability-test-needs-the-refuting-channel-in-the-fixture`):

    B   ``ss_intron_boundary``  (stage 3's destinations), on the boundary axis
    E   ``solvable_exon``       (stage 4's destinations), on the region axis

DELIVER rows are the slots whose certified truth is EXACTLY pure gDNA (the claim is true and must be
delivered); REFUTE rows carry real RNA (the claim must be overturned by evidence). ⛔ A pooled number
over the whole library cannot judge pass-0: outside its claimed slots the fan-out is silence BY
DESIGN, so a whole-library comparison charges it for slots it deliberately leaves unsolved — that
comparison lives in `ladder_arm_ab.py --arm policy_fanout`, as context.

Each policy runs the FULL pipeline in-process off the cached scan and oracle
(`pass0_vs_oracle.measure_condition`), so the three columns differ by exactly the message policy.
``--self-test`` falsifies the SCORER: injected error at a claimed slot must be caught in the right
population and split, and error parked on an unclaimed slot must be scored NOWHERE.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import build_region_statics  # noqa: E402
from rigel.calibration.splice_graph import build_boundary_flags_array  # noqa: E402
from rigel.calibration.structural_claims import build_structural_claims  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import read_scan_cache  # noqa: E402

#: the three policies, as configs — the SHIPPED config with exactly one axis varied.
POLICIES = {
    "silent": lambda: CalibrationConfig(message_propagation=False),
    "relay": lambda: CalibrationConfig(),
    "fanout": lambda: CalibrationConfig(message_policy="fanout"),
}


def claimed_masks(chain, claims, truth: dict) -> dict:
    """The two claimed populations on their OWN axes, plus the DELIVER/REFUTE split from certified
    truth (pure ⇔ every RNA column exactly zero at the slot)."""
    out = {}
    for tag, kind_val, slot_mask in (
        ("B", BOUNDARY, np.asarray(claims.ss_intron_boundary, bool)),
        ("E", REGION, np.asarray(claims.solvable_exon, bool)),
    ):
        is_k = np.asarray(truth["kind"]) == kind_val
        obj = np.asarray(truth["obj"], np.int64)[is_k]
        n_axis = int(obj.max()) + 1 if obj.size else 0
        true_g = np.zeros(n_axis)
        true_g[obj] = np.asarray(truth["n_gdna"], np.float64)[is_k]
        rna = np.zeros(n_axis)
        rna[obj] = (
            np.asarray(truth["n_nrna"], np.float64) + np.asarray(truth["n_mrna"], np.float64)
        )[is_k]
        mask = np.zeros(n_axis, bool)
        slot_is_k = np.asarray(chain.kind) == kind_val
        mask[np.asarray(chain.obj_idx, np.int64)[slot_is_k]] = slot_mask[slot_is_k]
        out[tag] = {
            "mask": mask,
            "true_gdna": true_g,
            "deliver": mask & (rna == 0.0),
            "refute": mask & (rna > 0.0),
        }
    return out


def score(est_gdna_by_axis: dict, masks: dict) -> dict:
    """``Σ|est − true|`` per population and split — the estimate arrays come per axis
    (``mass_gdna_boundary`` for B, ``mass_gdna_region`` for E)."""
    rows = {}
    for tag, m in masks.items():
        err = np.abs(np.asarray(est_gdna_by_axis[tag], np.float64) - m["true_gdna"])
        rows[tag] = {
            "claimed": float(err[m["mask"]].sum()),
            "deliver": float(err[m["deliver"]].sum()),
            "refute": float(err[m["refute"]].sum()),
        }
    return rows


def audit(index, region_arrays, suite: Path, condition: str, work_dir: Path, policies) -> dict:
    P0 = _sibling("pass0_vs_oracle.py")
    t = dict(np.load(Path(suite) / "oracle_cache" / condition / "slot_truth.npz"))
    cache = read_scan_cache(Path(suite) / "oracle_cache" / condition / "_main", index)
    chain = build_region_chain(cache.payload.ref_region_offsets, cache.payload.ref_boundary_offsets)
    statics = build_region_statics(chain, region_arrays, build_boundary_flags_array(index))
    masks = claimed_masks(chain, build_structural_claims(chain, statics), t)
    out = {}
    for name in policies:
        m = P0.measure_condition(
            bam=str(Path(suite) / condition / "sim_oracle.bam"),
            index=index,
            pipeline_config=PipelineConfig(),
            calibration_config=POLICIES[name](),
            work_dir=work_dir,
            tag=condition,
            truth_pmfs=None,
            oracle_cache=Path(suite) / "oracle_cache",
        )
        res = m.arms["final"] if "final" in m.arms else m.arms[sorted(m.arms)[0]]
        out[name] = score({"B": res.mass_gdna_boundary, "E": res.mass_gdna_region}, masks)
    return out


def report(condition: str, rows: dict, policies) -> None:
    print(f"\n== {condition}")
    print(f"   {'pop':<4} {'split':<9}" + "".join(f"{p:>12}" for p in policies))
    for tag in ("B", "E"):
        for split in ("claimed", "deliver", "refute"):
            print(
                f"   {tag if split == 'claimed' else '':<4} {split:<9}"
                + "".join(f"{rows[p][tag][split]:>12.0f}" for p in policies)
            )


def self_test() -> int:
    ok = 0

    def check(name, cond):
        nonlocal ok
        print(f"   {'✔' if cond else '✘'} {name}")
        if not cond:
            raise SystemExit(f"self-test FAILED at: {name}")
        ok += 1

    # a synthetic 9-slot chain: ig B exon B intron B exon B ig — one claimed boundary pair, and the
    # claims derived by the REAL builder so the masks cannot drift from the substrate's definition.
    from rigel.calibration.region_geometry import RegionStatics
    from rigel.calibration.splice_graph import FLAG_DONOR_POS, FLAG_TSS_POS

    chain = build_region_chain(np.array([0, 5]), np.array([0, 4]))
    fp = np.array([0, 0, 1, 1, 1, 1, 1, 0, 0], bool)
    mp = np.array([0, 0, 1, 0, 0, 0, 1, 0, 0], bool)
    bflags = np.zeros(9, np.uint16)
    bflags[1] = bflags[7] = FLAG_TSS_POS
    bflags[3] = bflags[5] = FLAG_DONOR_POS
    statics = RegionStatics(
        n_slots=9,
        free_pos=fp,
        free_neg=np.zeros(9, bool),
        mrna_active_pos=mp,
        mrna_active_neg=np.zeros(9, bool),
        boundary_flags=np.where(np.asarray(chain.kind) == BOUNDARY, bflags, 0).astype(np.uint16),
    )
    claims = build_structural_claims(chain, statics)
    truth = {
        "kind": np.asarray(chain.kind),
        "obj": np.asarray(chain.obj_idx),
        "n_gdna": np.array([9, 2, 1, 8, 6, 8, 1, 2, 9], float),
        "n_nrna": np.array([0, 0, 0, 0, 3, 0, 0, 0, 0], float),
        "n_mrna": np.array([0, 0, 5, 0, 0, 0, 5, 0, 0], float),
    }
    masks = claimed_masks(chain, claims, truth)
    check(
        "the claimed populations are the substrate's (2 boundaries, 2 exons)",
        masks["B"]["mask"].sum() == 2 and masks["E"]["mask"].sum() == 2,
    )
    check(
        "DELIVER and REFUTE split by certified purity and never overlap",
        masks["B"]["deliver"].sum() == 2
        and masks["E"]["refute"].sum() == 2
        and not (masks["B"]["deliver"] & masks["B"]["refute"]).any(),
    )
    # perfect estimates score zero; an injected error is caught in the right population and split
    truth_b = np.zeros(4)
    truth_b[np.asarray(chain.obj_idx)[np.asarray(chain.kind) == BOUNDARY]] = truth["n_gdna"][
        np.asarray(chain.kind) == BOUNDARY
    ]
    perfect = {"B": truth_b, "E": truth["n_gdna"][np.asarray(chain.kind) == REGION]}
    check(
        "perfect estimates score exactly zero everywhere",
        all(v == 0.0 for row in score(perfect, masks).values() for v in row.values()),
    )
    bad = {k: v.copy().astype(float) for k, v in perfect.items()}
    bad["B"][1] += 7.0  # boundary obj 1 = the claimed donor (slot 3): DELIVER side
    bad["E"][2] += 5.0  # region obj 2 = the intron — NOT a claimed exon: must score nowhere
    got = score(bad, masks)
    check(
        "an injected boundary error lands in B/deliver, exactly",
        got["B"]["claimed"] == 7.0 and got["B"]["deliver"] == 7.0 and got["B"]["refute"] == 0.0,
    )
    check(
        "error parked on an UNCLAIMED slot is scored nowhere (judge only what is claimed)",
        got["E"]["claimed"] == 0.0,
    )
    bad2 = {k: v.copy().astype(float) for k, v in perfect.items()}
    bad2["E"][2 + 0] = bad2["E"][2] + 0.0  # no-op guard
    bad2["E"][np.flatnonzero(masks["E"]["refute"])[0]] += 3.0
    got2 = score(bad2, masks)
    check(
        "an injected exon error lands in E/refute, exactly",
        got2["E"]["refute"] == 3.0 and got2["E"]["deliver"] == 0.0,
    )
    print(f"self-test: {ok}/{ok}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None)
    ap.add_argument("--policies", nargs="*", default=list(POLICIES), choices=list(POLICIES))
    ap.add_argument(
        "--work-dir",
        type=Path,
        default=Path.home() / ".cache" / "rigel_pass0_claimed",
    )
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
    for c in conds:
        report(
            c,
            audit(index, region_arrays, args.suite, c, args.work_dir, args.policies),
            args.policies,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
