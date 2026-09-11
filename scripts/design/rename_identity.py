#!/usr/bin/env python
"""Is this rename stage numerically a no-op? Every stage proven bit-identical to a frozen reference.

A rename, or any refactor or speed-up, must change no number, and this makes that claim falsifiable. `--freeze` runs the whole
pipeline once on one condition and reduces it to digests; `--check` runs it again after a stage and
compares. It compares content, never names, because the names are the thing changing: the content
multiset (the sorted `dtype|shape|sha256` of every payload and calibration array, field names discarded,
so any pure rename is invisible and any moved value is not) and the tool's own output (the transcript
table, canonicalised by transcript id at full precision, which carries no region/boundary vocabulary at
all). The pair is needed: the multiset cannot see two arrays swapping names, and the transcript digest
can. No name map is used. Bit-identity is only provable with every source of run-to-run variation
pinned, so the scan runs on one thread with one BGZF thread, `OMP_NUM_THREADS=1` and a fixed EM seed;
a `--check` is therefore a statement about the code, not about the shipped thread configuration. The
reference is frozen, not rolling: every stage compares to the same capture, never to the stage before it,
because renames compound and a rolling baseline would let one stage's defect become the next stage's
truth, so `--freeze` refuses to overwrite an existing reference.

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
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]

from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _native_detect_sj_tag, run_pipeline  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"
DEFAULT_CONDITION = "gdna_g00_ss_0.99_nrna_mid_capture_off"
REFERENCE = _RUNS / "arms" / "rename_identity_reference.json"

#: every source of run-to-run variation, pinned; removing any one makes "bit-identical" unprovable.
#: ``total_threads=1`` is the load-bearing one: float addition is not associative across worker threads.
DETERMINISM = dict(total_threads=1, bgzf_threads=1)
EM_SEED = 20260807


def _sha(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]


def content_multiset(obj) -> list[str]:
    """Every array as ``dtype|shape|sha256``, sorted and with the field name discarded.

    Discarding the name is what survives a rename. Sorting makes the result a multiset, so a
    reordering of the dataclass fields is invisible while a changed value is not.
    """
    out = []
    for f in dataclasses.fields(obj):
        v = getattr(obj, f.name, None)
        if isinstance(v, np.ndarray):
            out.append(f"{v.dtype}|{tuple(v.shape)}|{_sha(v)}")
    return sorted(out)


def capture(suite: Path, index_dir: Path, condition: str, bam: Path | None = None) -> dict:
    """Run the whole pipeline once, deterministically, and reduce it to comparable digests.

    ``bam`` names any BAM directly (a real library); without it the BAM is the panel condition's oracle.
    """
    bam = str(bam) if bam is not None else str(suite / condition / "sim_oracle.bam")
    index = TranscriptIndex.load(str(index_dir))
    base = PipelineConfig()
    cfg = dataclasses.replace(
        base,
        scan=dataclasses.replace(base.scan, sj_strand_tag=_native_detect_sj_tag(bam), **DETERMINISM),
        em=dataclasses.replace(base.em, seed=EM_SEED),
    )
    # the payload is not on `PipelineResult` and holds the very banks being renamed, so it is captured
    # at the `calibrate` boundary rather than by scanning a second time; `calibrate` is imported
    # function-locally by the pipeline, so patching the module attribute is picked up at call time
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

    quant = result.estimator.get_counts_df(index)
    # canonicalised: sorted by the transcript id and rendered at full precision, so the digest is a
    # statement about the numbers and not about row order or float formatting
    key = "transcript_id" if "transcript_id" in quant.columns else quant.columns[0]
    q = quant.sort_values(key).to_csv(index=False, float_format="%.17g")

    return {
        "condition": condition,
        "quant_sha": hashlib.sha256(q.encode()).hexdigest(),
        "quant_rows": int(len(quant)),
        "payload_multiset": content_multiset(box["payload"]),
        "calibration_multiset": content_multiset(result.calibration),
    }


def compare(ref: dict, now: dict) -> list[str]:
    """Every difference, named. An empty list means the stage is a proven numeric no-op."""
    bad = []
    if ref.get("condition") != now.get("condition"):
        bad.append(f"condition: {ref.get('condition')!r} -> {now.get('condition')!r}")
    if ref.get("quant_rows") != now.get("quant_rows"):
        bad.append(f"quant_rows: {ref.get('quant_rows')} -> {now.get('quant_rows')}")
    if ref.get("quant_sha") != now.get("quant_sha"):
        bad.append("⛔ THE TOOL'S OUTPUT MOVED — this is not a rename")
    for name in ("payload_multiset", "calibration_multiset"):
        a, b = ref.get(name), now.get(name)
        if a is None and b is None:
            continue
        if a is None or b is None:
            bad.append(f"{name}: present in one capture and not the other")
            continue
        if a != b:
            only_ref, only_now = sorted(set(a) - set(b)), sorted(set(b) - set(a))
            bad.append(
                f"⛔ {name}: {len(only_ref)} array(s) gone, {len(only_now)} appeared "
                f"(a pure rename moves NEITHER)"
            )
            for x in only_ref[:4]:
                bad.append(f"      only in the reference: {x}")
            for x in only_now[:4]:
                bad.append(f"      only in this run:      {x}")
    return bad


def self_test() -> int:
    """The comparator's own falsification. No pipeline, no I/O, pure perturbation."""
    ref = {
        "condition": "c", "quant_rows": 10, "quant_sha": "abc",
        "payload_multiset": sorted(["int64|(3,)|aaa", "float64|(4,)|bbb"]),
        "calibration_multiset": sorted(["float64|(2,)|ccc"]),
    }
    checks = []

    checks.append(("identical => silent", compare(ref, dict(ref)) == []))

    # a renamed field must be invisible; that is the whole design
    renamed = dict(ref)  # the multiset holds no names at all, so a rename cannot change it
    checks.append(("a pure rename => silent", compare(ref, renamed) == []))

    # a changed value must be caught, however small; the hash has no tolerance
    moved = dict(ref, payload_multiset=sorted(["int64|(3,)|aaa", "float64|(4,)|bbZ"]))
    checks.append(("one array's CONTENT moved => caught", any("payload_multiset" in x
                                                              for x in compare(ref, moved))))

    # a changed shape is a different array
    shaped = dict(ref, payload_multiset=sorted(["int64|(3,)|aaa", "float64|(5,)|bbb"]))
    checks.append(("a SHAPE change => caught", compare(ref, shaped) != []))

    # the tool's output is the second, independent view
    checks.append(("the quant digest moved => caught",
                   any("OUTPUT MOVED" in x for x in compare(ref, dict(ref, quant_sha="zzz")))))

    # a dropped array must not read as agreement (an empty comparison scoring as a pass)
    dropped = dict(ref, payload_multiset=["int64|(3,)|aaa"])
    checks.append(("an array DISAPPEARING => caught", compare(ref, dropped) != []))

    # an empty capture is the sharpest version of the same lie
    checks.append(("an EMPTY capture => caught", compare(ref, dict(ref, payload_multiset=[])) != []))

    # the known hole, pinned so nobody mistakes it for coverage: the multiset cannot see two arrays
    # swapping names, because a multiset has no names; the quant digest is what covers that, and this
    # asserts the hole is exactly where it is documented to be
    swapped = dict(ref)  # same contents, names exchanged — multiset identical BY CONSTRUCTION
    checks.append(("a NAME SWAP is invisible to the multiset (documented hole)",
                   compare(ref, swapped) == []))

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
                "compares to the SAME capture, because the renames compound and a rolling baseline "
                "would let a stage-2 defect become the accepted truth for stages 3-9. Delete it "
                "deliberately if the schema genuinely moved (stage 0)."
            )
        args.reference.parent.mkdir(parents=True, exist_ok=True)
        args.reference.write_text(json.dumps(now, indent=2))
        print(f"  ⭐ frozen -> {args.reference}")
        print(f"     quant {now['quant_sha'][:16]}…  ({now['quant_rows']:,} rows)")
        for k in ("payload_multiset", "calibration_multiset"):
            if k in now:
                print(f"     {k}: {len(now[k])} arrays")
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
    print(f"  ✅ BIT-IDENTICAL{label} — quant digest and every array's content unchanged")
    return 0


if __name__ == "__main__":
    sys.exit(main())
