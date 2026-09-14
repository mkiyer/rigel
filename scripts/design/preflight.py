#!/usr/bin/env python3
"""Can this session run and regenerate everything? One command, one verdict, before anything else.

A session that assumes its environment is fine discovers otherwise hours in, usually as a confusing
wrong number rather than an error. This checks first that the toolchain (the `rigel` conda env, the
compiled native extension, the `rigel` CLI), both references (the panel's genome/GTF/index/probes and
the test chromosome's YAML with its renders), both panels (scan caches, oracle caches carrying every part
`panel.py` requires, the certified `slot_truth.npz`) and every `scripts/design/` instrument (each one imports)
are present and working. It changes nothing and measures nothing: every check is a read or an import,
and a failed check prints the exact command that regenerates the missing artifact, because most of this
data is derived and a missing derived artifact is a command not yet run rather than damage. An import
break is the check that actually rots after a `src/` deletion or a rename, so the default path is the
fast one; the `--self-test` sweep over every instrument is opt-in via `--full`, runs each one in a
subprocess, and is what to run after a deposit-rule change, a default flip, or before a commit that
touches many instruments.

Usage::

    python scripts/design/preflight.py            # the fast path: toolchain, data, imports
    python scripts/design/preflight.py --full     # + every instrument's --self-test, in parallel
    python scripts/design/preflight.py --self-test
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
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


class Report:
    """Every check as ``(ok, what, detail, fix)``; a missing thing always names its command."""

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
    # the method-development reference: one hand-edited YAML in the repo, the GTFs / abundances /
    # probe panels rendered beside it, everything else derived
    spec = REPO / "scripts" / "sim" / "test_reference" / "test_chr.yaml"
    gtf = spec.parent / "test_chr.gtf"
    rep.add(spec.is_file(), "test chromosome YAML (the ONE hand-edited file)", str(spec))
    rep.add(test_chromosome_renders_in_sync(spec), "test chromosome renders match the YAML", str(gtf),
            "python scripts/sim/build_test_reference.py")
    rep.add((TESTREF / "test_chr.fa").is_file(), "test chromosome FASTA (DERIVED)", str(TESTREF),
            "python scripts/sim/build_test_reference.py")
    rep.add((TESTREF / "idx").is_dir(), "test chromosome rigel index (DERIVED)", "",
            "rigel index --fasta <T>/test_chr.fa --gtf <T>/test_chr.gtf --no-mappability "
            "--no-tsv -o <T>/idx")


def _sim_script(name: str):
    """A ``scripts/sim/`` script loaded by path, so a check here reads the rule from the script that
    builds what it checks instead of restating it."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / "sim" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_chromosome_renders_in_sync(spec: Path) -> bool:
    """Do the checked-in renders (GTFs, abundances, probe panels) match what the YAML renders now?
    In-process and read-only — the builder's own `--check`, without a subprocess."""
    if not spec.is_file():
        return False
    module = _sim_script("build_test_reference")
    loaded = module.load_spec(spec)
    return not module.check_spec(loaded) and not module.check_renders(loaded, spec.parent)


def _conditions(root: Path) -> list[str]:
    d = root / "scan_cache"
    return sorted(p.name for p in d.iterdir()) if d.is_dir() else []


def check_panel(rep: Report, root: Path, label: str, rebuild: str) -> None:
    conds = _conditions(root)
    rep.add(bool(conds), f"{label}: scan caches", f"{len(conds)} conditions", rebuild)
    if not conds:
        return
    oracle = root / "oracle_cache"
    parts = _sim_script("panel").ORACLE_PARTS
    complete, partial = [], []
    for c in conds:
        done = all((oracle / c / part / "payload.npz").is_file() for part in parts)
        (complete if done else partial).append(c)
    rep.add(not partial, f"{label}: oracle caches, every part",
            f"{len(complete)}/{len(conds)} complete"
            + (f"; missing on {partial[0]}…" if partial else ""),
            rebuild.replace("simulate", "cache"))
    stamped = [c for c in conds if (oracle / c / "slot_truth.npz").is_file()]
    rep.add(len(stamped) == len(conds), f"{label}: certified slot_truth", f"{len(stamped)}/{len(conds)}",
            "python scripts/design/calibration_oracle.py")


def check_instruments(rep: Report, full: bool) -> None:
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
    if not full:
        rep.add(True, "instrument --self-test sweep",
                "SKIPPED (default; --full runs it)")
        return
    # the outer pool is deliberately narrow because several instruments parallelise themselves (by
    # arm or by condition), so a wide outer pool oversubscribes the machine and starves the one that
    # dominates the wall time; three lets the quick instruments overlap while leaving cores for it
    have = [p for p in files if '"--self-test"' in p.read_text()]
    workers = max(1, min(3, (os.cpu_count() or 2) - 1))

    def _run(path):
        r = subprocess.run([sys.executable, str(path), "--self-test"], capture_output=True,
                           text=True, cwd=str(REPO))
        return path.name, r.returncode

    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_run, have))
    failed = sorted(name for name, code in results if code != 0)
    rep.add(not failed, "every instrument's --self-test passes",
            f"{len(have) - len(failed)}/{len(have)}" + (f"; {failed}" if failed else ""),
            "run the named instrument's --self-test and read its output")


def self_test() -> int:
    """The reporter perturbed, with no I/O: a preflight that cannot report a failure is decoration."""
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
    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--full", action="store_true",
                    help="also run every instrument's --self-test, in parallel")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    tool, ref, panel, testref, inst = Report(), Report(), Report(), Report(), Report()
    check_toolchain(tool)
    check_reference(ref)
    check_panel(panel, LADDER, "panel (16 conditions)",
                "python scripts/sim/panel.py simulate --config scripts/sim/configs/gdna_ladder.yaml --jobs 16")
    check_panel(testref, TESTREF / "scenarios", "test chromosome",
                "python scripts/sim/panel.py simulate --config scripts/sim/configs/test_reference.yaml")
    check_instruments(inst, args.full)

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
        print("⭐ everything present and working: the toolchain, the references, both panels and every "
              "instrument.")
    print("\n⚠ This checks PRESENCE and IMPORTABILITY. It does not run the suite — do that too:")
    print("   python -m pytest tests/ -q   (ANY failure is a regression; the baseline count lives in CLAUDE.md)")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
