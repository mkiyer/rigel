#!/usr/bin/env python3
"""How does a prototype message policy score, per gene type and per slot, against certified truth?
This is the harness a message mechanism is developed on before it touches `src/`. A prototype is a
Python class with the backbone's ``Policy`` shape (``prepare(ctx)`` returning an object with
``propagate`` / ``solve``), named in an ``ARMS = {"my_arm": PolicyClass, ...}`` table in the module
``--module`` points at; the harness installs it in place of ``calibrate``'s ``TransferPolicy`` for
the ``transfer`` arm and scores it beside ``silent`` and the shipped ``transfer`` on a cached,
certified condition, with the shipped class restored afterwards. The error is |gDNA estimate -
truth| in fragments against `slot_truth.npz`; no EM runs and nothing is re-scanned. Three views,
never pooled with each other: the whole-library number per condition with its region/boundary
split and a per-gene-type table (a gene's type is the token after ``gB<k>_`` in its test-chromosome
``gene_id``, so a moved number names its structure); ``--by-class``, the same error summed per node
class (certified stratum, a boundary's terminus and junction flags, an exon's reach: licensed /
edge / walled), which judges a message at its destinations where the whole-library number carries
the refit prior's response; and ``dissect``, every slot of one gene type with its truth beside every
arm. A prototype that subclasses the shipped policy inherits every mechanism the parent later
gains, so compare `src` against `src` across a landing. `backbone_parity.py` imports `load_arms`.

Usage::

    python scripts/design/policy_prototype.py --panel test --arms transfer my_arm \\
        --module path/to/proto.py --conditions gdna_g50_ss_0.99_nrna_file_capture_on
    python scripts/design/policy_prototype.py --panel test --arms transfer my_arm --module p.py --all
    python scripts/design/policy_prototype.py --panel test --arms transfer my_arm --module p.py --all --by-class
    python scripts/design/policy_prototype.py dissect --panel test --arms transfer my_arm --module p.py \\
        --condition gdna_g50_ss_0.99_nrna_file_capture_on --type capnasc
    python scripts/design/policy_prototype.py --self-test
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import numpy as np  # noqa: E402

from rigel.calibration.messages.transfer import TransferPolicy  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    build_boundary_flags_array,
    build_sj_geometry_arrays,
)
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

RUNS = Path.home() / "Downloads" / "rigel_runs"
PANELS = {
    "test": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios"),
    "test_sparse": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios_probes_sparse"),
    "test_junction": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios_probes_junction"),
    "ladder": (RUNS / "suite" / "rigel_index", RUNS / "suite" / "ladder"),
}
#: the MODULE (the package exports a function of the same name — patching that would be a silent no-op)
CALMOD = importlib.import_module("rigel.calibration.calibrate")


# ── gene types from the index ─────────────────────────────────────────────────────────────────────


def gene_type_token(gene_id: str) -> str:
    """The structure token of a test-chromosome gene: ``gB3_capnasc`` -> ``capnasc``; anything else
    (a real annotation's gene ids) -> ``-``."""
    if gene_id.startswith("gB") and "_" in gene_id:
        head, _, tail = gene_id.partition("_")
        if head[2:].isdigit() and tail:
            return tail
    return "-"


def slot_types(index, ra, chain, kind, obj) -> np.ndarray:
    """Per slot, the gene type of the annotated gene whose span covers it (a REGION by containment, a
    BOUNDARY by its left region), ``-`` where none does. Derived from the index's gene table, never
    from coordinates."""
    starts = np.asarray(ra.start, np.int64)
    ends = np.asarray(ra.end, np.int64)
    ref_id = np.asarray(ra.ref_id)
    region_type = np.full(starts.shape[0], "-", dtype=object)
    genes = index.g_df
    genes = genes[~genes["is_synthetic"].astype(bool)] if "is_synthetic" in genes.columns else genes
    for g in genes.itertuples(index=False):
        rid = index.ref_name_to_id.get(g.ref)
        if rid is None:
            continue
        inside = (ref_id == rid) & (starts >= int(g.start)) & (ends <= int(g.end)) & (region_type == "-")
        region_type[inside] = gene_type_token(str(g.g_id))
    left = np.asarray(chain.left, np.int64)
    out = np.full(kind.shape[0], "-", dtype=object)
    is_r = kind == REGION
    out[is_r] = region_type[np.clip(obj[is_r], 0, starts.shape[0] - 1)]
    for b in np.flatnonzero(kind == BOUNDARY):
        if left[b] >= 0:
            out[b] = out[left[b]]
    return out.astype(str)


# ── scoring ──────────────────────────────────────────────────────────────────────────────────────


def load_panel(panel):
    idx_dir, suite = PANELS[panel]
    index = TranscriptIndex.load(str(idx_dir))
    ra = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    return index, ra, build_sj_geometry_arrays(index), build_boundary_flags_array(index), suite


def condition_setup(index, ra, sj, bflags, suite, cond):
    cache_dir = suite / "oracle_cache" / cond
    cache = read_scan_cache(cache_dir / "_main", index)
    kw = calibration_inputs(cache, index)
    payload = kw["payload"]
    truth = np.load(cache_dir / "slot_truth.npz", allow_pickle=True)
    kind = np.asarray(truth["kind"])
    obj = np.asarray(truth["obj"], np.int64)
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    types = slot_types(index, ra, chain, kind, obj)
    kwargs = dict(
        region_arrays=ra,
        strand_model=kw["strand_model"],
        gdna_fl_pmf=kw["gdna_fl_pmf"],
        rna_fl_pmf=kw["rna_fl_pmf"],
        sj=sj,
        boundary_flags=bflags,
    )
    return dict(
        payload=payload,
        kwargs=kwargs,
        kind=kind,
        obj=obj,
        strata=np.asarray(truth["stratum"]).astype(str),
        truth_gdna=np.asarray(truth["n_gdna"], float),
        count=np.asarray(truth["count"], float),
        types=types,
        chain=chain,
    )


def run_arm(arm, arms, c):
    """One arm's per-slot gDNA estimate. ``silent`` and ``transfer`` are the shipped policies; any other
    name is looked up in ``arms`` and installed in place of ``TransferPolicy``."""
    if arm == "silent":
        cfg = CalibrationConfig(message_propagation=False)
        CALMOD.TransferPolicy = TransferPolicy
    else:
        cfg = CalibrationConfig(message_propagation=True, message_policy="transfer")
        CALMOD.TransferPolicy = TransferPolicy if arm == "transfer" else arms[arm]
    try:
        res = CALMOD.calibrate(payload=c["payload"], config=cfg, **c["kwargs"])
    finally:
        CALMOD.TransferPolicy = TransferPolicy
    kind, obj = c["kind"], c["obj"]
    est = np.zeros(kind.shape[0])
    est[kind == REGION] = np.asarray(res.mass_gdna_region, float)[obj[kind == REGION]]
    est[kind == BOUNDARY] = np.asarray(res.mass_gdna_boundary, float)[obj[kind == BOUNDARY]]
    return est


def slot_classes(c, bflags):
    """One label per slot: the certified stratum (``R exon``, ``B exon|exon`` …), a boundary's terminus /
    junction flags, and an exon's reach — ``licensed`` (an intron|exon face without a terminus),
    ``edge`` (an intergenic|exon edge and no licensed face), or ``walled`` (neither)."""
    from rigel.calibration.messages.transfer_rows import SJ_FLAGS, TERMINUS

    kind, obj, chain, strata = c["kind"], c["obj"], c["chain"], c["strata"]
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    flags = np.asarray(bflags, np.uint16)
    labels = np.array(list(strata), dtype=object)
    is_bnd = kind == BOUNDARY
    for b in np.flatnonzero(is_bnd):
        f = int(flags[obj[b]])
        labels[b] = strata[b] + (" [term]" if f & TERMINUS else "") + (" [sj]" if f & SJ_FLAGS else "")
    for e in np.flatnonzero(strata == "R exon"):
        reach = "walled"
        for b in (left[e], right[e]):
            if b < 0 or not is_bnd[b]:
                continue
            o = right[b] if left[b] == e else left[b]
            if o < 0:
                continue
            if strata[o] == "R intron" and not (int(flags[obj[b]]) & TERMINUS):
                reach = "licensed"
                break
            if strata[o] == "R intergenic":
                reach = "edge"
        labels[e] = f"R exon ({reach})"
    return labels


def score(panel, arm_names, arms, conds, by_class=False):
    index, ra, sj, bflags, suite = load_panel(panel)
    if conds == ["all"]:
        conds = sorted(d.name for d in (suite / "oracle_cache").iterdir() if (d / "slot_truth.npz").exists())
    for cond in conds:
        c = condition_setup(index, ra, sj, bflags, suite, cond)
        line, per, axis = {}, {}, {}
        for arm in ["silent", *arm_names]:
            t0 = time.perf_counter()
            err = np.abs(run_arm(arm, arms, c) - c["truth_gdna"])
            line[arm] = (float(err.sum()), time.perf_counter() - t0)
            per[arm] = {t: float(err[c["types"] == t].sum()) for t in np.unique(c["types"])}
            axis[arm] = (float(err[c["kind"] == REGION].sum()), float(err[c["kind"] == BOUNDARY].sum()))
        print(f"\n== {cond}")
        print("  whole-library |err|: " + "  ".join(f"{a} {v[0]:,.0f} ({v[1]:.1f}s)" for a, v in line.items()))
        print("  region/boundary:     " + "  ".join(f"{a} {axis[a][0]:,.0f}/{axis[a][1]:,.0f}" for a in line))
        print(f"  {'type':<14}" + "".join(f"{a:>11}" for a in line))
        for t in sorted(per["silent"]):
            if t != "-":
                print(f"  {t:<14}" + "".join(f"{per[a][t]:>11,.0f}" for a in line))
        if by_class:
            labels = slot_classes(c, bflags)
            errs = {arm: np.abs(run_arm(arm, arms, c) - c["truth_gdna"]) for arm in ["silent", *arm_names]}
            mass = c["count"]
            print(f"  {'node class':<36}{'slots':>7}{'mass':>11}" + "".join(f"{a:>12}" for a in errs))
            for lab in sorted(set(labels), key=lambda x: -float(mass[labels == x].sum())):
                m = labels == lab
                print(f"  {lab:<36}{int(m.sum()):>7}{float(mass[m].sum()):>11,.0f}" + "".join(f"{float(errs[a][m].sum()):>12,.0f}" for a in errs))


def dissect(panel, arm_names, arms, cond, gtype):
    index, ra, sj, bflags, suite = load_panel(panel)
    c = condition_setup(index, ra, sj, bflags, suite, cond)
    ests = {arm: run_arm(arm, arms, c) for arm in ["silent", *arm_names]}
    kind, cnt, tg, types = c["kind"], c["count"], c["truth_gdna"], c["types"]
    starts = np.asarray(ra.start, np.int64)
    ends = np.asarray(ra.end, np.int64)
    obj = c["obj"]
    print(f"\n== DISSECT {cond}  type {gtype}  (R region / B boundary; f_g = gDNA share; |err| per arm)")
    print(f"  {'slot':>6}{'pos':>10} {'k':<3}{'count':>7}{'true_fg':>8}" + "".join(f"{a:>9}" for a in ests) + "   |err|")
    for i in np.flatnonzero(types == gtype):
        pos = int(starts[obj[i]]) if kind[i] == REGION else int(ends[obj[i]])
        fg = tg[i] / cnt[i] if cnt[i] > 0 else float("nan")
        row = f"  {i:>6}{pos:>10} {'R' if kind[i] == REGION else 'B':<3}{cnt[i]:>7.0f}{fg:>8.3f}"
        row += "".join(f"{(ests[a][i] / cnt[i] if cnt[i] > 0 else float('nan')):>9.3f}" for a in ests)
        row += "   " + " ".join(f"{abs(ests[a][i] - tg[i]):>5.0f}" for a in ests)
        print(row)


def load_arms(module_path):
    if module_path is None:
        return {}
    spec = importlib.util.spec_from_file_location("policy_prototype_arms", module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return dict(module.ARMS)


def self_test() -> int:
    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"   ⛔ {name}")

    check("a test-chromosome gene id yields its type token", gene_type_token("gB3_capnasc") == "capnasc")
    check("an isoform id keeps the whole token", gene_type_token("gB1_capaltstart") == "capaltstart")
    check("a real annotation's gene id yields no type", gene_type_token("ENSG00000123") == "-")
    check("a malformed block id yields no type", gene_type_token("gBx_clean") == "-")
    # the arm installation is scoped: after run_arm the shipped policy is back in place
    check("the shipped policy is the module's default", CALMOD.TransferPolicy is TransferPolicy)
    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", nargs="?", default="score", choices=("score", "dissect"))
    ap.add_argument("--panel", default="test", choices=sorted(PANELS))
    ap.add_argument("--arms", nargs="+", default=["transfer"])
    ap.add_argument("--module", type=Path, default=None, help="a .py defining ARMS = {name: PolicyClass}")
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--condition", default=None)
    ap.add_argument("--type", default=None)
    ap.add_argument("--by-class", action="store_true", help="also sum each arm's error per node class")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()
    arms = load_arms(args.module)
    if args.mode == "dissect":
        if not (args.condition and args.type):
            raise SystemExit("dissect needs --condition and --type")
        dissect(args.panel, args.arms, arms, args.condition, args.type)
        return 0
    conds = ["all"] if args.all else (args.conditions or [])
    if not conds:
        raise SystemExit("name --conditions or pass --all")
    score(args.panel, args.arms, arms, conds, by_class=args.by_class)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
