#!/usr/bin/env python
"""Is this rename stage numerically a no-op? Every stage proven bit-identical to a frozen reference.

A rename, refactor or speed-up must change no number, and this makes that claim falsifiable. `--freeze` runs
the whole pipeline once on one condition and reduces it to digests; `--check` runs it again after a stage and
compares. It compares content, never names, because the names are what changes. Two views:
- the content multiset of the calibration payload and of the calibration result: every value either holds,
  recursing into the dataclasses they nest (the payload's deferred bank, gap census, scan QC and drain QC), each
  array as `dtype|shape|sha256` and each scalar as `type|value` at full precision (the capture weight, the
  saturation, the bulk density, the sense fraction, every count). Field names are discarded, so a pure rename is
  invisible and a moved value is not. A field of any other type stops the run: a skipped value would be certified
  unseen;
- the tool's output: the transcript, gene, loci and nascent tables, each sorted by its id (its first column) at
  full precision.
The pair is needed: the multiset cannot see two values swapping names, and the table digests can. Every source of
run-to-run variation is pinned (one scan thread, one BGZF thread, `OMP_NUM_THREADS=1`, a fixed EM seed), so a
`--check` is a statement about the code, not about the shipped thread configuration. The reference is frozen,
never rolling: every stage compares to the same capture, so one stage's defect cannot become the next stage's
truth, and `--freeze` refuses to overwrite it. A reference frozen by an older digest does not hold this record's
fields and never certifies it: re-freeze from the tree the stage starts on.

Usage::

    python scripts/design/rename_identity.py --self-test          # perturb the comparator; no I/O
    python scripts/design/rename_identity.py --freeze             # capture the reference, once
    python scripts/design/rename_identity.py --check --stage s1   # after every stage
    python scripts/design/rename_identity.py --check --suite DIR --index DIR --condition NAME
    python scripts/design/rename_identity.py --freeze --bam real.bam --index DIR --reference REF.json
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import os
import sys
from collections import Counter
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np  # noqa: E402


from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import run_pipeline  # noqa: E402

from _shared import DEFAULT_INDEX, DEFAULT_SUITE, ORACLE_BAM, RUNS  # noqa: E402

DEFAULT_CONDITION = "gdna_g00_ss_0.99_nrna_mid_capture_off"
REFERENCE = RUNS / "arms" / "rename_identity_reference.json"

#: every source of run-to-run variation, pinned; removing any one makes "bit-identical" unprovable.
#: ``total_threads=1`` is the load-bearing one: float addition is not associative across worker threads.
DETERMINISM = dict(total_threads=1, bgzf_threads=1)
EM_SEED = 20260807

#: the tool's output tables, every one the estimator writes; a record that lacks one gates less than it claims
TABLES = ("transcript", "gene", "loci", "nascent")
#: the record's fields; a record missing any one (a reference frozen by an older digest) certifies nothing
RECORD_KEYS = ("condition", "tables", "payload_multiset", "calibration_multiset")


def _sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]


def table_digest(table) -> str:
    """One output table as one digest: sorted by its first column (the table's id) and rendered at full precision
    without the header, so it states the numbers and their column order, not the row order, the float formatting or
    a column's name. A renamed column is invisible; a moved value or a reordered column is not."""
    text = table.sort_values(table.columns[0], kind="stable").to_csv(index=False, header=False, float_format="%.17g")
    return hashlib.sha256(text.encode()).hexdigest()


def content_multiset(obj) -> list[str]:
    """Every value a dataclass holds, recursing into the dataclasses it nests, sorted and with the field names
    discarded: an array as ``dtype|shape|sha256``, a scalar as ``type|value`` at full precision, an absent optional
    as ``None``.

    Discarding the name is what survives a rename. Sorting makes the result a multiset, so a reordering of the
    dataclass fields is invisible while a changed value is not. A value of any other type raises: skipping it
    would certify a value the digest never saw.
    """
    out: list[str] = []
    _digest_fields(obj, type(obj).__name__, out)
    return sorted(out)


def _digest_fields(obj, where: str, out: list[str]) -> None:
    for f in dataclasses.fields(obj):
        _digest_value(getattr(obj, f.name), f"{where}.{f.name}", out)


def _digest_value(v, where: str, out: list[str]) -> None:
    if isinstance(v, np.ndarray):
        if v.dtype.hasobject:
            raise TypeError(f"⛔ {where} is an object array: its bytes are pointers, not values, so it has no digest")
        out.append(f"{v.dtype}|{tuple(v.shape)}|{_sha(v)}")
    elif dataclasses.is_dataclass(v) and not isinstance(v, type):
        _digest_fields(v, where, out)
    elif v is None:
        out.append("None")
    elif isinstance(v, (bool, np.bool_)):
        out.append(f"bool|{bool(v)}")
    elif isinstance(v, (int, np.integer)):
        out.append(f"int|{int(v)}")
    elif isinstance(v, (float, np.floating)):
        out.append(f"float|{float(v)!r}")
    elif isinstance(v, str):
        out.append(f"str|{hashlib.sha256(v.encode()).hexdigest()[:16]}")
    else:
        raise TypeError(
            f"⛔ {where} is a {type(v).__name__}, which the digest cannot read; teach `_digest_value` its content "
            "rather than skip it, because a skipped value is certified unseen"
        )


def capture(suite: Path, index_dir: Path, condition: str, bam: Path | None = None) -> dict:
    """Run the whole pipeline once, deterministically, and reduce it to comparable digests.

    ``bam`` names any BAM directly (a real library); without it the BAM is the panel condition's oracle.
    """
    bam = str(bam) if bam is not None else str(suite / condition / ORACLE_BAM)
    index = TranscriptIndex.load(str(index_dir))
    base = PipelineConfig()
    cfg = dataclasses.replace(
        base,
        scan=dataclasses.replace(base.scan, **DETERMINISM),
        em=dataclasses.replace(base.em, seed=EM_SEED),
    )
    # the payload is not on `PipelineResult`, so it is captured at the `calibrate` boundary; the pipeline
    # imports `calibrate` function-locally, so patching the module attribute is picked up at call time
    import rigel.calibration as CAL

    original = CAL.calibrate
    box: dict = {}

    def _wrapper(*a, **kw):
        box["payload"] = kw["payload"]
        return original(*a, **kw)

    CAL.calibrate = _wrapper
    try:
        result = run_pipeline(bam, index, cfg)
    finally:
        CAL.calibrate = original
    if "payload" not in box:
        raise RuntimeError("⛔ calibrate was never called — the capture boundary moved, and a capture "
                           "missing the payload would silently gate less than it claims")

    est = result.estimator
    tables = {
        "transcript": est.get_counts_df(index),
        "gene": est.get_gene_counts_df(index),
        "loci": est.get_loci_df(index),
        "nascent": est.get_nascent_quant_df(index),
    }
    return reduce_run(condition, tables, box["payload"], result.calibration)


def reduce_run(condition: str, tables: dict, payload, calibration) -> dict:
    """One run's outputs reduced to the comparable record: every output table's digest and row count, and the
    content multisets of the payload and the calibration result."""
    if set(tables) != set(TABLES):
        raise ValueError(f"⛔ the record digests the tables {TABLES}; this run offered {sorted(tables)}")
    return {
        "condition": condition,
        "tables": {name: {"rows": int(len(tables[name])), "sha": table_digest(tables[name])} for name in TABLES},
        "payload_multiset": content_multiset(payload),
        "calibration_multiset": content_multiset(calibration),
    }


def compare(ref: dict, now: dict) -> list[str]:
    """Every difference, named. An empty list means the stage is a proven numeric no-op."""
    missing = sorted({k for k in RECORD_KEYS if k not in ref} | {k for k in RECORD_KEYS if k not in now})
    if missing:
        return [f"⛔ the two records do not hold the same fields (missing {missing}): a reference frozen by an older "
                "digest cannot certify this run; re-freeze one from the tree the stage starts on"]
    bad = []
    if ref["condition"] != now["condition"]:
        bad.append(f"condition: {ref['condition']!r} -> {now['condition']!r}")
    for name in sorted(set(ref["tables"]) | set(now["tables"])):
        a, b = ref["tables"].get(name), now["tables"].get(name)
        if a is None or b is None:
            bad.append(f"⛔ the {name} table: present in one capture and not the other")
            continue
        if a["rows"] != b["rows"]:
            bad.append(f"the {name} table's rows: {a['rows']:,} -> {b['rows']:,}")
        if a["sha"] != b["sha"]:
            bad.append(f"⛔ THE TOOL'S OUTPUT MOVED: the {name} table — this is not a rename")
    for name in ("payload_multiset", "calibration_multiset"):
        a, b = Counter(ref[name]), Counter(now[name])
        if a != b:
            only_ref, only_now = sorted((a - b).elements()), sorted((b - a).elements())
            bad.append(
                f"⛔ {name}: {len(only_ref)} value(s) gone, {len(only_now)} appeared "
                f"(a pure rename moves NEITHER)"
            )
            for x in only_ref[:4]:
                bad.append(f"      only in the reference: {x}")
            for x in only_now[:4]:
                bad.append(f"      only in this run:      {x}")
    return bad


def _run_perturbations() -> list[tuple[str, bool]]:
    """A synthetic run shaped like the pipeline's (four output tables, a payload with nested banks and scalars, a
    calibration with library scalars), reduced, then perturbed one value at a time. Each perturbation must move the
    record; the unperturbed copy must not."""
    import pandas as pd

    Bank = dataclasses.make_dataclass("Bank", [("offsets", np.ndarray), ("values", np.ndarray)])
    Census = dataclasses.make_dataclass("Census", [("resolved", int), ("deferred", int)])
    Drain = dataclasses.make_dataclass("Drain", [("offered", int), ("census_before", Census)])
    Payload = dataclasses.make_dataclass("Payload", [
        ("counts", np.ndarray), ("deferred", Bank), ("gap_resolution", Census), ("drain", object),
        ("max_length", int), ("graph_hash", object),
    ])
    Calibration = dataclasses.make_dataclass("Calibration", [
        ("count_gdna_region", np.ndarray), ("gdna_bulk_density", object), ("gdna_capture_weight", float),
        ("gdna_capture_saturation", float), ("rna_sense_frac", float), ("n_regions", int),
    ])

    def run():
        tables = {
            "transcript": pd.DataFrame({"transcript_id": ["t2", "t1"], "locus_id": [0, 0], "count": [1.0, 2.0]}),
            "gene": pd.DataFrame({"gene_id": ["g1"], "locus_id": [0], "count": [3.0]}),
            "loci": pd.DataFrame({"locus_id": [0], "count_gdna": [4.0], "has_gdna": [True]}),
            "nascent": pd.DataFrame({"nascent_id": ["n1"], "locus_id": [0], "count": [5.0]}),
        }
        payload = Payload(
            counts=np.arange(4, dtype=np.uint32),
            deferred=Bank(offsets=np.array([0, 2], dtype=np.int64), values=np.array([7.0, 8.0])),
            gap_resolution=Census(resolved=3, deferred=1),
            drain=Drain(offered=5, census_before=Census(resolved=2, deferred=2)),
            max_length=1000, graph_hash="abc",
        )
        calibration = Calibration(
            count_gdna_region=np.array([0.5, 1.5]), gdna_bulk_density=0.01, gdna_capture_weight=0.75,
            gdna_capture_saturation=4.0, rna_sense_frac=0.99, n_regions=2,
        )
        return tables, payload, calibration

    base = reduce_run("c", *run())
    tiny = 2.0**-40

    def moved(plant) -> bool:
        tables, payload, calibration = run()
        plant(tables, payload, calibration)
        return compare(base, reduce_run("c", tables, payload, calibration)) != []

    def nested(p, path, value):
        *head, last = path
        obj = p
        for name in head:
            obj = getattr(obj, name)
        setattr(obj, last, value)

    out = [("a whole run, reduced twice => silent", compare(base, reduce_run("c", *run())) == [])]
    out.append(("a planted nested-bank value (payload.deferred) => caught",
                moved(lambda t, p, c: p.deferred.values.__setitem__(0, 7.0 + tiny))))
    out.append(("a planted nested scalar (payload.drain.census_before) => caught",
                moved(lambda t, p, c: nested(p, ["drain", "census_before", "resolved"], 99))))
    out.append(("a planted top-level payload scalar => caught",
                moved(lambda t, p, c: setattr(p, "max_length", 1001))))
    out.append(("a nested bank appearing where there was None => caught",
                moved(lambda t, p, c: setattr(p, "drain", None))))
    for field in ("gdna_capture_weight", "gdna_capture_saturation", "gdna_bulk_density", "rna_sense_frac"):
        out.append((f"a planted calibration scalar ({field}) => caught",
                    moved(lambda t, p, c, f=field: setattr(c, f, getattr(c, f) + tiny))))
    # the transcript table is the control: the old digest saw it, so this plant proves the plant mechanism works
    for table, column in (("transcript", "count"), ("gene", "count"), ("loci", "count_gdna"), ("nascent", "count")):
        out.append((f"a planted {table}-table cell => caught",
                    moved(lambda t, p, c, k=table, col=column: t.update({k: t[k].assign(**{col: t[k][col] + tiny})}))))
    out.append(("a flipped loci-table flag => caught",
                moved(lambda t, p, c: t.update(loci=t["loci"].assign(has_gdna=[False])))))
    out.append(("a transcript-table row order alone => silent",
                not moved(lambda t, p, c: t.__setitem__("transcript", t["transcript"].iloc[::-1]))))

    # a field the digest cannot read must stop the run, never be skipped: a skip is a value certified unseen
    tables, payload, calibration = run()
    payload.graph_hash = {"a": 1}
    try:
        reduce_run("c", tables, payload, calibration)
        refused = False
    except TypeError:
        refused = True
    out.append(("a field the digest cannot read => refused", refused))
    return out


def self_test() -> int:
    """The comparator's own falsification. No pipeline, no I/O, pure perturbation."""
    ref = {
        "condition": "c",
        "tables": {name: {"rows": 10, "sha": f"{name}-abc"} for name in TABLES},
        "payload_multiset": sorted(["int64|(3,)|aaa", "float64|(4,)|bbb", "int|7"]),
        "calibration_multiset": sorted(["float64|(2,)|ccc", "float|0.5"]),
    }
    checks = []

    checks.append(("identical => silent", compare(ref, dict(ref)) == []))

    # a renamed field must be invisible; that is the whole design: the multiset holds no names
    A = dataclasses.make_dataclass("A", [("x", np.ndarray), ("y", np.ndarray)])
    B = dataclasses.make_dataclass("B", [("renamed_x", np.ndarray), ("renamed_y", np.ndarray)])
    u, v = np.arange(3), np.ones(4)
    checks.append(("a pure rename => silent", content_multiset(A(u, v)) == content_multiset(B(u, v))))

    # an output table's digest: a renamed column is invisible, a moved value is not
    import pandas as pd

    table = pd.DataFrame({"transcript_id": ["t2", "t1"], "count": [1.0, 2.0]})
    checks.append(("a renamed transcript-table column => silent",
                   table_digest(table) == table_digest(table.rename(columns={"count": "count_annotated"}))))
    checks.append(("a moved transcript-table value => caught",
                   table_digest(table) != table_digest(table.assign(count=[1.0, 2.0 + 2**-40]))))

    # a changed value must be caught, however small; the hash has no tolerance
    moved = dict(ref, payload_multiset=sorted(["int64|(3,)|aaa", "float64|(4,)|bbZ", "int|7"]))
    checks.append(("one array's CONTENT moved => caught", any("payload_multiset" in x
                                                              for x in compare(ref, moved))))

    # a changed shape is a different array
    shaped = dict(ref, payload_multiset=sorted(["int64|(3,)|aaa", "float64|(5,)|bbb", "int|7"]))
    checks.append(("a SHAPE change => caught", compare(ref, shaped) != []))

    # a value repeated one time more or less is a changed multiset, though its set is the same
    doubled = dict(ref, calibration_multiset=sorted(["float64|(2,)|ccc", "float|0.5", "float|0.5"]))
    checks.append(("a value's MULTIPLICITY changed => caught",
                   any("1 appeared" in x for x in compare(ref, doubled))))

    # the tool's output is the second, independent view, table by table
    for name in TABLES:
        tables = dict(ref["tables"], **{name: {"rows": 10, "sha": "zzz"}})
        checks.append((f"the {name} table's digest moved => caught",
                       any("OUTPUT MOVED" in x and name in x for x in compare(ref, dict(ref, tables=tables)))))
    dropped_table = {k: t for k, t in ref["tables"].items() if k != "gene"}
    checks.append(("a table DISAPPEARING => caught", compare(ref, dict(ref, tables=dropped_table)) != []))

    # a dropped array must not read as agreement (an empty comparison scoring as a pass)
    dropped = dict(ref, payload_multiset=["int64|(3,)|aaa", "int|7"])
    checks.append(("an array DISAPPEARING => caught", compare(ref, dropped) != []))

    # an empty capture is the sharpest version of the same lie
    checks.append(("an EMPTY capture => caught", compare(ref, dict(ref, payload_multiset=[])) != []))

    # a reference frozen by an older digest has fewer fields; it must refuse, never pass on what it shares
    old = {"condition": "c", "quant_sha": "abc", "quant_rows": 10,
           "payload_multiset": ref["payload_multiset"], "calibration_multiset": ref["calibration_multiset"]}
    checks.append(("a reference of an older digest => refused", compare(old, ref) != []))

    # the known hole, pinned: the multiset cannot see two arrays swapping names (it has no names); the
    # table digests cover that
    checks.append(("a NAME SWAP is invisible to the multiset (documented hole)",
                   content_multiset(A(u, v)) == content_multiset(A(v, u))))

    # a whole run through the reduction `capture` uses: every value the record certifies must move it
    for name, ok in _run_perturbations():
        checks.append((name, ok))

    width = max(len(n) for n, _ in checks)
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}s}")
    n = sum(ok for _, ok in checks)
    print(f"\n{n}/{len(checks)} comparator perturbations fire")
    return 0 if n == len(checks) else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--freeze", action="store_true", help="capture the reference (do this ONCE)")
    ap.add_argument("--check", action="store_true", help="compare the current tree to the reference")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=DEFAULT_CONDITION)
    ap.add_argument("--bam", type=Path, default=None,
                    help="a BAM to capture directly (a real library); the label is its library directory")
    ap.add_argument("--reference", type=Path, default=REFERENCE)
    ap.add_argument("--stage", default="", help="label recorded with a --check, for the log")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if not (args.freeze or args.check):
        raise SystemExit("one of --freeze / --check / --self-test is required")

    condition = args.condition if args.bam is None else args.bam.resolve().parent.parent.name
    now = capture(args.suite, args.index, condition, bam=args.bam)

    if args.freeze:
        if args.reference.is_file():
            raise SystemExit(
                f"⛔ {args.reference} already exists. The reference is FROZEN on purpose — every stage "
                "compares to the SAME capture, because the stages compound and a rolling baseline "
                "would let one stage's defect become the accepted truth for every later stage. Delete "
                "it deliberately if the capture's schema genuinely moved."
            )
        args.reference.parent.mkdir(parents=True, exist_ok=True)
        args.reference.write_text(json.dumps(now, indent=2))
        print(f"  ⭐ frozen -> {args.reference}")
        for name, t in now["tables"].items():
            print(f"     {name:<10s} {t['sha'][:16]}…  ({t['rows']:,} rows)")
        for k in ("payload_multiset", "calibration_multiset"):
            print(f"     {k}: {len(now[k])} values")
        return 0

    if not args.reference.is_file():
        raise SystemExit(f"⛔ no reference at {args.reference} — run --freeze first")
    ref = json.loads(args.reference.read_text())
    bad = compare(ref, now)
    label = f" [{args.stage}]" if args.stage else ""
    if bad:
        print(f"  ⛔ NOT BIT-IDENTICAL{label} — this stage changed a number, so it is not a rename")
        for b in bad:
            print(f"     {b}")
        return 1
    print(f"  ✅ BIT-IDENTICAL{label} — every output table and every payload and calibration value unchanged")
    return 0


if __name__ == "__main__":
    sys.exit(main())
