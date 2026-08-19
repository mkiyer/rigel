#!/usr/bin/env python3
"""PREFLIGHT — can this session actually run and REGENERATE everything? One command, one verdict.

⭐⭐⭐ **WHAT THIS IS FOR.** Owner, 2026-08-19: *"I want the next session to verify that they have access
to all of the necessary tools and can rerun and regenerate everything."* A session that starts by
assuming its environment is fine discovers otherwise three hours in, usually as a confusing wrong number
rather than an error. This checks — before any of that — that the toolchain, the reference data and the
instruments are present and working, and for anything MISSING it prints the exact command that rebuilds
it.

⛔ **IT CHANGES NOTHING AND MEASURES NOTHING.** Every check is a read or an import. A ✘ is never fatal
here: much of this data is DERIVED and the point is to say which command derives it.

==========================  =========================================================================
**toolchain**               the `rigel` conda env is active, the compiled native extension imports,
                            and the `rigel` CLI is on PATH
**reference data**          the panel's genome/GTF/index/probes, and the method-development test
                            reference — each with the command that regenerates it
**panel**                   the 16 simulated conditions, their scan caches, and the oracle caches
                            with all FIVE partitions (the three ORIGINS plus `rna_pos`/`rna_neg`)
**certification**           every condition's `slot_truth.npz` and its stamp
**instruments**             every `scripts/design/` file imports, and every one that has a
                            `--self-test` passes it
==========================  =========================================================================

⚠ **The instrument sweep is the slow half (~2 min)**; `--fast` skips it and checks everything else.

Usage::

    python scripts/design/preflight.py            # everything
    python scripts/design/preflight.py --fast     # skip the instrument self-test sweep
    python scripts/design/preflight.py --self-test
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = Path(__file__).resolve().parents[2]
RUNS = Path.home() / "Downloads" / "rigel_runs"
SUITE = RUNS / "suite"
LADDER = SUITE / "ladder"
TESTREF = RUNS / "test_reference"
#: the five partitions an oracle-cached condition must carry — three ORIGINS plus the per-STRAND RNA
#: pair the three-arm map needs (`calibration/_oracle.RNA_STRAND_ORIGINS`)
ORACLE_PARTS = ("gdna", "mrna", "nrna", "rna_pos", "rna_neg")


class Report:
    """Every check as ``(ok, what, detail, fix)`` — and a MISSING thing always names its command."""

    def __init__(self) -> None:
        self.rows: list[tuple[bool, str, str, str]] = []

    def add(self, ok: bool, what: str, detail: str = "", fix: str = "") -> bool:
        self.rows.append((bool(ok), what, detail, fix))
        return bool(ok)

    def print(self, title: str) -> None:
        print(f"\n══ {title}")
        for ok, what, detail, fix in self.rows:
            print(f"   {'✔' if ok else '✘'} {what:<58} {detail}")
            if not ok and fix:
                print(f"       ↳ {fix}")

    @property
    def n_bad(self) -> int:
        return sum(1 for ok, *_ in self.rows if not ok)


def check_toolchain(rep: Report) -> None:
    rep.add(os.environ.get("CONDA_DEFAULT_ENV") == "rigel"
            or "rigel" in sys.prefix,
            "conda env `rigel` is active", sys.prefix,
            'source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel')
    try:
        import rigel  # noqa: F401
        from rigel import native  # noqa: F401
        ok, detail = True, "rigel + native extension import"
    except Exception as exc:  # noqa: BLE001
        ok, detail = False, f"{type(exc).__name__}: {exc}"
    rep.add(ok, "the compiled native extension imports", detail,
            'pip install --no-build-isolation -e ".[dev]"')
    cli = shutil.which("rigel")
    rep.add(bool(cli), "the `rigel` CLI is on PATH", cli or "not found",
            'pip install --no-build-isolation -e ".[dev]"')


def check_reference(rep: Report) -> None:
    rep.add((SUITE / "reference" / "genome.fa").is_file(), "panel reference genome + GTF",
            str(SUITE / "reference"),
            "scripts/sim/build_suite_reference.py --fasta <source.fa> --gtf <source.gtf> "
            "--refs chr21 chr22 --ercc -o <reference>")
    rep.add((SUITE / "rigel_index").is_dir(), "panel rigel index", str(SUITE / "rigel_index"),
            "python scripts/sim/panel.py build --config scripts/sim/configs/gdna_ladder.yaml")
    rep.add((SUITE / "reference" / "capture_panel.tsv").is_file(), "panel capture probes", "",
            "python scripts/sim/panel.py build --config scripts/sim/configs/gdna_ladder.yaml")
    # the METHOD-DEVELOPMENT reference — hand-edited sources in the repo, everything else derived
    gtf = REPO / "scripts" / "sim" / "test_reference" / "test_chr.gtf"
    rep.add(gtf.is_file(), "test chromosome GTF (hand-edited, in the repo)", str(gtf))
    rep.add((TESTREF / "test_chr.fa").is_file(), "test chromosome FASTA (DERIVED)", str(TESTREF),
            "python scripts/sim/build_test_reference.py")
    rep.add((TESTREF / "idx").is_dir(), "test chromosome rigel index (DERIVED)", "",
            "rigel index --fasta <T>/test_chr.fa --gtf <T>/test_chr.gtf --no-mappability "
            "--no-tsv -o <T>/idx")


def _conditions(root: Path) -> list[str]:
    d = root / "scan_cache"
    return sorted(p.name for p in d.iterdir()) if d.is_dir() else []


def check_panel(rep: Report, root: Path, label: str, rebuild: str) -> None:
    conds = _conditions(root)
    rep.add(bool(conds), f"{label}: scan caches", f"{len(conds)} conditions", rebuild)
    if not conds:
        return
    oracle = root / "oracle_cache"
    complete, partial = [], []
    for c in conds:
        have = {p.name for p in (oracle / c).iterdir()} if (oracle / c).is_dir() else set()
        (complete if set(ORACLE_PARTS) <= have else partial).append(c)
    rep.add(not partial, f"{label}: oracle caches, all five partitions",
            f"{len(complete)}/{len(conds)} complete"
            + (f"; missing on {partial[0]}…" if partial else ""),
            rebuild.replace("simulate", "cache"))
    stamped = [c for c in conds if (oracle / c / "slot_truth.npz").is_file()]
    rep.add(len(stamped) == len(conds), f"{label}: certified slot_truth", f"{len(stamped)}/{len(conds)}",
            "python scripts/design/calibration_oracle.py")


def check_instruments(rep: Report, fast: bool) -> None:
    import importlib.util

    design = REPO / "scripts" / "design"
    files = sorted(p for p in design.glob("*.py") if p.name != "__init__.py")
    broken = []
    for path in files:
        name = f"_preflight_{path.stem}"
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        sys.path.insert(0, str(path.parent))
        try:
            spec.loader.exec_module(module)
        except SystemExit:
            pass
        except Exception as exc:  # noqa: BLE001
            broken.append(f"{path.name}: {type(exc).__name__}")
        finally:
            sys.path.remove(str(path.parent))
            sys.modules.pop(name, None)
    rep.add(not broken, "every scripts/design/ instrument imports",
            f"{len(files)} files" + (f"; {broken[:2]}" if broken else ""),
            "an instrument that cannot import cannot be run — repair or delete it")
    if fast:
        rep.add(True, "instrument --self-test sweep", "SKIPPED (--fast)")
        return
    have = [p for p in files if '"--self-test"' in p.read_text()]
    failed = []
    for path in have:
        r = subprocess.run([sys.executable, str(path), "--self-test"], capture_output=True,
                           text=True, cwd=str(REPO))
        if r.returncode != 0:
            failed.append(path.name)
    rep.add(not failed, "every instrument's --self-test passes",
            f"{len(have) - len(failed)}/{len(have)}" + (f"; {failed}" if failed else ""),
            "run the named instrument's --self-test and read its output")


def self_test() -> int:
    """⛔ The reporter perturbed, with no I/O — a preflight that cannot report a failure is decoration."""
    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"   ⛔ {name}")

    r = Report()
    check("a passing check counts as good", r.add(True, "x") is True and r.n_bad == 0)
    check("a failing check is counted", r.add(False, "y", "d", "fix me") is False and r.n_bad == 1)
    check("a second failure accumulates", r.add(False, "z") is False and r.n_bad == 2)
    check("truthiness is normalised to bool", r.add([], "empty") is False)

    import io
    import contextlib

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        r.print("T")
    out = buf.getvalue()
    check("a failure prints its fix", "fix me" in out)
    check("a pass does not print a fix line", out.count("↳") == 1)
    check("both marks are rendered", "✔" in out and "✘" in out)
    check("the five oracle partitions are named", set(ORACLE_PARTS) >= {"rna_pos", "rna_neg"})
    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fast", action="store_true", help="skip the instrument --self-test sweep")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    tool, ref, panel, testref, inst = Report(), Report(), Report(), Report(), Report()
    check_toolchain(tool)
    check_reference(ref)
    check_panel(panel, LADDER, "panel (16 conditions)",
                "python scripts/sim/panel.py simulate --config scripts/sim/configs/gdna_ladder.yaml --jobs 16")
    check_panel(testref, TESTREF / "scenarios", "test chromosome (8 scenarios)",
                "python scripts/sim/simulate_reads.py --config scripts/sim/configs/test_reference.yaml -j 8")
    check_instruments(inst, args.fast)

    tool.print("TOOLCHAIN")
    ref.print("REFERENCE DATA")
    panel.print("THE PANEL — the benchmark")
    testref.print("THE TEST CHROMOSOME — where the policy is developed")
    inst.print("INSTRUMENTS")

    bad = sum(r.n_bad for r in (tool, ref, panel, testref, inst))
    print()
    if bad:
        print(f"⛔ {bad} check(s) failed. Each one above names the command that regenerates it.")
        print("   ⚠ A missing DERIVED artifact is not damage — it is a command you have not run yet.")
    else:
        print("⭐ everything present and working: the toolchain, both references, both panels "
              "(five oracle partitions each), and every instrument.")
    print("\n⚠ This checks PRESENCE and IMPORTABILITY. It does not run the suite — do that too:")
    print("   python -m pytest tests/ -q   (baseline: 0 failed / 3,554 passed / 9 xfailed)")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
