#!/usr/bin/env python3
"""Can this session run and regenerate everything? One command, one verdict, before anything else.

A session that assumes its environment is fine discovers otherwise hours in, usually as a confusing
wrong number rather than an error. This checks first that the toolchain (the `rigel` conda env, the
compiled native extension, the `rigel` CLI), both references (the panel's genome/GTF/index/probes and
the test chromosome's YAML with its renders), both panels (scan caches and oracle caches carrying every
part `panel.py` requires, each one VALID: its key checked against this build and the panel's index by the
cache module's own check, so a stale cache never reads as present; the certified `slot_truth.npz`) and
every `scripts/design/` instrument (each one imports) are present and working. It changes nothing and
measures nothing: every check is a read or an import, and a failed check prints the exact command that
regenerates the missing or stale artifact, because most of this data is derived and a missing derived
artifact is a command not yet run rather than damage. An import break is the check that actually rots
after a `src/` deletion or a rename, so the default path is the fast one; the `--self-test` sweep over
every instrument is opt-in via `--full`, runs each one in a subprocess, and is what to run after a
deposit-rule change, a default flip, or before a commit that touches many instruments.

Usage::

    python scripts/design/preflight.py            # the fast path: toolchain, data and cache keys, imports
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

from _shared import RUNS  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
SUITE = RUNS / "suite"
SUITE_INDEX = SUITE / "rigel_index"
LADDER = SUITE / "ladder"
TESTREF = RUNS / "test_reference"
TESTREF_INDEX = TESTREF / "idx"


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
    rep.add(SUITE_INDEX.is_dir(), "panel rigel index", str(SUITE_INDEX),
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
    rep.add(TESTREF_INDEX.is_dir(), "test chromosome rigel index (DERIVED)", "",
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
    return sorted(p.name for p in d.iterdir() if p.is_dir()) if d.is_dir() else []


def _refusal(check, *args) -> str:
    """Why ``check(*args)`` refuses a cache, as the first sentence of its refusal, or "" when it accepts.

    A refusal is what `build_scan_cache.py` re-scans and the instruments refuse to read: a key that does
    not describe this build and index, a manifest without a key it needs, a file that is not there.
    Anything else escapes, because it is not a stale cache."""
    from rigel.scan_cache import ScanCacheKeyError

    try:
        check(*args)
    except ScanCacheKeyError as exc:
        return str(exc).split(". ")[0]
    except (FileNotFoundError, KeyError) as exc:
        return f"{type(exc).__name__}: {exc}"
    return ""


def _cache_refusal(cache_dir: Path, index) -> str:
    """"" when the cache at *cache_dir* is one this build reads against *index*, else why it is not.

    The cache module's own check, from the manifest alone: `scan_cache.check_scan_cache`, the one
    `build_scan_cache.py` decides a skip by, under the scan settings every instrument reads a cache with. The BAM is
    compared at its recorded path when it is still there, and its bytes are read only when its size matches and its
    modification time does not."""
    from rigel.config import BamScanConfig
    from rigel.scan_cache import check_scan_cache

    if not (cache_dir / "payload.npz").is_file():
        return "missing"
    return _refusal(check_scan_cache, cache_dir, index, None, BamScanConfig())


def _first_part_refusal(cond_dir: Path, parts, index) -> tuple[str, str]:
    """``(name, why)`` for one condition's oracle cache: its first part that is not valid, or the
    condition with "" when every part is."""
    for part in parts:
        why = _cache_refusal(cond_dir / part, index)
        if why:
            return f"{cond_dir.name}/{part}", why
    return cond_dir.name, ""


def cache_verdict(refusals: list[tuple[str, str]]) -> tuple[bool, str]:
    """``(ok, detail)`` over caches given as ``(name, why)``, ``why`` "" for a valid one: how many are
    valid, and the first stale one named with its reason. No cache at all is not a pass."""
    stale = [(name, why) for name, why in refusals if why]
    detail = f"{len(refusals) - len(stale)}/{len(refusals)} valid"
    if stale:
        detail += f"; {stale[0][0]}: {stale[0][1]}"
    return bool(refusals) and not stale, detail


def _load_index(index_dir: Path):
    """``(index, "")``, or ``(None, why)`` when the index is not there or does not load."""
    from rigel.index import IndexLoadError, TranscriptIndex

    if not index_dir.is_dir():
        return None, f"no index at {index_dir}"
    try:
        return TranscriptIndex.load(str(index_dir)), ""
    except (IndexLoadError, ValueError) as exc:
        return None, f"the index does not load: {str(exc).split('. ')[0]}"


def check_panel(rep: Report, root: Path, index_dir: Path, label: str, rebuild: str) -> None:
    """A panel's caches, each VALID: keyed to this build and the panel's index, as the instruments that read it
    require, so a stale cache never reads as present."""
    conds = _conditions(root)
    if not conds:
        rep.add(False, f"{label}: scan caches", "0 conditions", rebuild)
        return
    recache = rebuild.replace("simulate", "cache")
    oracle = root / "oracle_cache"
    index, why = _load_index(index_dir)
    if index is None:
        for what in ("scan caches", "oracle caches, every part"):
            rep.add(False, f"{label}: {what}", f"{len(conds)} conditions, keys unchecked: {why}",
                    "the index must load first; REFERENCE DATA names the command that builds it")
    else:
        ok, detail = cache_verdict([(c, _cache_refusal(root / "scan_cache" / c, index)) for c in conds])
        rep.add(ok, f"{label}: scan caches", detail, recache)
        parts = _sim_script("panel").ORACLE_PARTS
        ok, detail = cache_verdict([_first_part_refusal(oracle / c, parts, index) for c in conds])
        rep.add(ok, f"{label}: oracle caches, every part", detail, recache)
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
    """The reporter and the cache verdict perturbed, with no I/O: a preflight that cannot report a
    failure is decoration."""
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

    # the cache rows: a stale cache is counted and the first one named with the reason the cache module
    # gave; an error that is not a refusal is not turned into a verdict
    from rigel.scan_cache import ScanCacheKeyError

    check("every cache valid passes", cache_verdict([("a", ""), ("b", "")]) == (True, "2/2 valid"))
    check("a stale cache fails the row, counted, the FIRST named with its reason",
          cache_verdict([("a", ""), ("b", "moved"), ("c", "gone")]) == (False, "1/3 valid; b: moved"))
    check("no cache at all is not a pass", cache_verdict([])[0] is False)

    def refuse():
        raise ScanCacheKeyError("cache reach digest a != index reach digest b. The reaches moved. Re-scan.")

    def broken():
        raise ValueError("not a stale cache")

    check("a refusal reads as its first sentence",
          _refusal(refuse) == "cache reach digest a != index reach digest b")
    check("an accepted cache has no refusal", _refusal(lambda: None) == "")
    check("a manifest without a key is a refusal",
          _refusal(lambda: {}["reach_digest"]) == "KeyError: 'reach_digest'")
    try:
        _refusal(broken)
        escaped = False
    except ValueError:
        escaped = True
    check("an error that is not a refusal escapes", escaped)
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
    check_panel(panel, LADDER, SUITE_INDEX, "panel (16 conditions)",
                "python scripts/sim/panel.py simulate --config scripts/sim/configs/gdna_ladder.yaml --jobs 16")
    check_panel(testref, TESTREF / "scenarios", TESTREF_INDEX, "test chromosome",
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
        print("   ⚠ A missing or stale DERIVED artifact is not damage — it is a command you have not run yet.")
    else:
        print("⭐ everything present and working: the toolchain, the references, both panels and every "
              "instrument.")
    print("\n⚠ This checks PRESENCE, the caches' KEYS and IMPORTABILITY. It does not run the suite — do "
          "that too:")
    print("   python -m pytest tests/ -q   (ANY failure is a regression; the baseline count lives in CLAUDE.md)")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
