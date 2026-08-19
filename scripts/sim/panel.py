#!/usr/bin/env python
"""⭐⭐⭐ THE SIMULATION + BENCHMARKING WORKFLOW, AS ONE COMMAND PER STAGE.

    panel.py status    what exists and what is missing — run this FIRST
    panel.py build     reference + index + capture probes
    panel.py simulate  the conditions
    panel.py cache     the scan cache AND the oracle cache
    panel.py score     run rigel end to end, score against per-fragment truth
    panel.py report    the tables

⛔⛔ **WHY THIS EXISTS: THE RECIPE STOPPED HALFWAY AND THE MISSING HALF WAS THE POINT.** `TESTING.md` §2
documented five manual shell steps ending at "cache the scans". **Running the tool and scoring it against
truth — the entire purpose — was in no recipe anywhere**, and had to be reassembled from a 56-row
instrument table. Worse, the ORACLE cache (the origin-split truth every scoring instrument reads) had no
step at all: it was a *side effect* of running `pass0_vs_oracle.py --oracle-cache`, so a reader who
followed the documented steps ended up with a panel that every scorer refused.

⭐ **THIS SCRIPT ADDS NO MEASUREMENT CODE.** Every stage shells out to the instrument that already owns
it — the simulator engine in `rigel.sim`, `build_scan_cache.py`, `pass0_vs_oracle.py`,
`quant_accuracy.py`. Duplicating a scorer is how a baseline and a ceiling drift apart
(`TRAPS: score-the-consumers-own-count`), so the value here is *sequencing and prerequisites*, not new
arithmetic. Anything this prints about a number, the underlying instrument printed first.

⭐ **EVERY PATH DERIVES FROM THE CONFIG**, which already carries `genome`, `gtf` and `outdir`. The index
is the one exception — it is a sibling of the reference by convention — so it is derived and
`--index` overrides it.

⚠ **`status` is the command to run when you do not know where you are.** Each stage is expensive and
resumable, and the failure mode this replaces is discovering after 20 minutes that step 3 never ran.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import yaml  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
DESIGN = REPO / "scripts" / "design"
SIM = REPO / "scripts" / "sim"

#: the three origin partitions the oracle cache holds, plus the undrained full payload.
ORACLE_PARTS = ("gdna", "mrna", "nrna", "_main")


class Panel:
    """Every path this workflow touches, derived from one config file."""

    def __init__(self, config: Path, index: Path | None = None):
        self.config = Path(config).resolve()
        if not self.config.is_file():
            raise SystemExit(f"⛔ no such config: {self.config}")
        cfg = yaml.safe_load(self.config.read_text()) or {}
        for key in ("genome", "gtf", "outdir"):
            if key not in cfg:
                raise SystemExit(
                    f"⛔ {self.config.name} has no `{key}:`. This workflow derives every path from the "
                    f"config; a panel that does not say where it lives cannot be driven by it."
                )
        self.cfg = cfg
        self.genome = Path(cfg["genome"])
        self.gtf = Path(cfg["gtf"])
        self.dir = Path(cfg["outdir"])
        self.reference = self.gtf.parent
        # ⚠ DERIVED, not configured: the index is a sibling of the reference directory. `--index`
        # overrides it rather than this guessing twice.
        self.index = Path(index) if index else self.reference.parent / "rigel_index"
        self.probes = self.reference / "capture_panel.tsv"
        self.scan_cache = self.dir / "scan_cache"
        self.oracle_cache = self.dir / "oracle_cache"
        self.arms = Path(os.environ.get("RIGEL_ARMS", self.dir / "arms"))

    @property
    def conditions(self) -> list[str]:
        """The conditions actually ON DISK — a simulated condition is one with an oracle BAM."""
        if not self.dir.is_dir():
            return []
        return sorted(p.name for p in self.dir.iterdir() if (p / "sim_oracle.bam").is_file())

    @property
    def planned(self) -> int:
        """How many conditions the config asks for, so `simulate` can report partial progress."""
        return len(self.cfg.get("conditions", []) or []) or 0


def run(cmd: list[str], *, what: str) -> None:
    """One stage. ⛔ A non-zero exit STOPS the workflow — a stage that failed must not look like a
    stage that was skipped, which is how a partial panel gets scored as a complete one."""
    print(f"\n\033[1m── {what}\033[0m\n   $ {' '.join(str(c) for c in cmd)}\n", flush=True)
    rc = subprocess.run([str(c) for c in cmd], cwd=str(REPO)).returncode
    if rc != 0:
        raise SystemExit(f"⛔ {what} FAILED (rc={rc}). Nothing after it has run.")


def need(ok: bool, what: str, fix: str) -> None:
    """A prerequisite, named, with the command that satisfies it."""
    if not ok:
        raise SystemExit(f"⛔ missing: {what}\n   fix: {fix}")


# ── the stages ───────────────────────────────────────────────────────────────────────────────────


def cmd_status(p: Panel, args) -> int:
    conds = p.conditions
    n_scan = len(list(p.scan_cache.glob("*/payload.npz"))) if p.scan_cache.is_dir() else 0
    # ⭐ an oracle condition is complete only when ALL FOUR parts are present; counting directories
    # would call a half-written condition done (`TRAPS: could-the-arm-have-fired`, storage form).
    n_oracle = 0
    if p.oracle_cache.is_dir():
        n_oracle = sum(
            all((p.oracle_cache / c / part / "payload.npz").is_file() for part in ORACLE_PARTS)
            for c in conds
        )
    arms = sorted(p.arms.glob("*.jsonl")) if p.arms.is_dir() else []

    def mark(ok):
        return "\033[32m✔\033[0m" if ok else "\033[31m✘\033[0m"

    print(f"\n  PANEL  {p.config.name}   ->  {p.dir}")
    print(f"  {'':<4}{'stage':<12} {'state':<34} path")
    print("  " + "-" * 96)
    rows = [
        ("build", p.genome.is_file() and p.gtf.is_file(), "reference", p.reference),
        ("build", p.index.is_dir(), "index", p.index),
        ("build", p.probes.is_file(), "capture probes", p.probes),
        ("simulate", bool(conds), f"{len(conds)}/{p.planned or '?'} conditions", p.dir),
        ("cache", n_scan >= len(conds) and n_scan > 0, f"scan cache {n_scan}/{len(conds)}",
         p.scan_cache),
        ("cache", n_oracle >= len(conds) and n_oracle > 0, f"oracle cache {n_oracle}/{len(conds)}",
         p.oracle_cache),
        ("score", bool(arms), f"{len(arms)} arm file(s)", p.arms),
    ]
    for stage, ok, what, path in rows:
        print(f"  {mark(ok)}  {stage:<12} {what:<34} {path}")
    if arms:
        print()
        for a in arms:
            print(f"       {a.name:<40} {a.stat().st_size / 1024:>8.0f} KB")
    nxt = next((s for s, ok, _w, _p in rows if not ok), None)
    print(f"\n  ⭐ next: `panel.py {nxt}`\n" if nxt else "\n  ⭐ every stage complete\n")
    return 0


def cmd_build(p: Panel, args) -> int:
    """Reference, index, capture probes. ⚠ The reference build is NOT driven from the panel config —
    it needs the source genome/GTF the panel was carved from, which the panel config does not name."""
    need(p.genome.is_file() and p.gtf.is_file(),
         f"the reference ({p.genome}, {p.gtf})",
         f"python {SIM / 'build_suite_reference.py'} --fasta <source.fa> --gtf <source.gtf> "
         f"--refs chr21 chr22 --ercc -o {p.reference}")
    if not p.index.is_dir() or args.force:
        run(["rigel", "index", "--fasta", p.genome, "--gtf", p.gtf,
             "--collapse-duplicate-transcripts", "--no-mappability", "--no-tsv", "-o", p.index],
            what=f"index -> {p.index}")
    else:
        print(f"   ✔ index exists: {p.index}  (--force to rebuild)")
    if not p.probes.is_file() or args.force:
        run([sys.executable, SIM / "design_suite_probes.py", "--gtf", p.gtf,
             "-o", p.probes, "--capture-fraction", "0.5"],
            what=f"capture probes -> {p.probes}")
    else:
        print(f"   ✔ probes exist: {p.probes}  (--force to rebuild)")
    return 0


def cmd_simulate(p: Panel, args) -> int:
    need(p.index.is_dir(), f"the index ({p.index})", "panel.py build")
    run([sys.executable, SIM / "simulate_reads.py", "--config", p.config, "-j", str(args.jobs)],
        what=f"simulate -> {p.dir}")
    return 0


def cmd_cache(p: Panel, args) -> int:
    """⛔ BOTH caches, and the oracle one is the reason this stage exists as a named step.

    The scan cache makes calibration re-runnable without rescanning. The ORACLE cache is the
    origin-split truth — `gdna` / `mrna` / `nrna` partitions plus the undrained `_main` payload — and
    every truth-scoring instrument refuses to run without it. It is populated as a side effect of
    `pass0_vs_oracle.py`, which is why it was invisible in the documented recipe."""
    need(bool(p.conditions), f"simulated conditions in {p.dir}", "panel.py simulate")
    conds = args.conditions or p.conditions
    # ⛔ `--force` MUST reach the scan cache, and this used to be the one stage it did not. A scan cache
    #    is keyed by the index, the graph and the scan CONFIG — never by the accumulator's output — so a
    #    change to what the scanner DEPOSITS leaves the key identical and the stale cache is silently
    #    reused. Measured 2026-08-19: after the chimera repair every condition reported `cached / skip`
    #    and the rebuild was a no-op (`TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule`
    #    is the change that exposed it).
    run([sys.executable, DESIGN / "build_scan_cache.py", "--index", p.index, "--suite", p.dir,
         "--out", p.scan_cache, *(["--force"] if args.force else []),
         *(["--conditions", *conds] if args.conditions else [])],
        what=f"scan cache -> {p.scan_cache}")
    # ⭐ `--jobs` reaches the oracle stage too. Building one condition's cache is a BAM split plus four
    #    scans, ~95 % of this stage's wall clock, and it saturates exactly ONE core at ~2 GB of real
    #    memory — so the stage is core-bound and shards cleanly. Measured 2026-08-19: serial, the stage
    #    used 1 of 16 cores.
    run([sys.executable, DESIGN / "pass0_vs_oracle.py", "--suite", p.dir, "--index", p.index,
         "--oracle-cache", p.oracle_cache, "--jobs", str(args.jobs),
         *(["--conditions", *conds] if args.conditions else [])],
        what=f"oracle cache (origin-split truth) -> {p.oracle_cache}")
    return 0


def cmd_score(p: Panel, args) -> int:
    need(p.oracle_cache.is_dir(), f"the oracle cache ({p.oracle_cache})", "panel.py cache")
    p.arms.mkdir(parents=True, exist_ok=True)
    for arm in args.arms:
        run([sys.executable, DESIGN / "quant_accuracy.py", "--arm", arm, "--jobs", str(args.jobs),
             "--suite", p.dir, "--index", p.index, "--oracle-cache", p.oracle_cache,
             "--out", p.arms / f"qa_{p.dir.name}_{arm}.jsonl",
             *(["--conditions", *args.conditions] if args.conditions else [])],
            what=f"score --arm {arm}")
    return 0


def cmd_report(p: Panel, args) -> int:
    files = [p.arms / f"qa_{p.dir.name}_{a}.jsonl" for a in args.arms]
    missing = [f for f in files if not f.is_file()]
    need(not missing, f"arm file(s) {[m.name for m in missing]}",
         f"panel.py score --arms {' '.join(args.arms)}")
    run([sys.executable, DESIGN / "quant_accuracy.py", "--report", *files], what="report")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="stages: status | build | simulate | cache | score | report",
    )
    ap.add_argument("stage", choices=("status", "build", "simulate", "cache", "score", "report"))
    ap.add_argument("--config", type=Path, required=True, help="a panel YAML in sim/configs/")
    ap.add_argument("--index", type=Path, default=None, help="override the derived index path")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--conditions", nargs="*", default=None, help="a subset, for a quick loop")
    ap.add_argument("--arms", nargs="+", default=["base"],
                    help="quant_accuracy arms: base, oracle, noop, base_reseed, oracle_gdna, …")
    ap.add_argument("--force", action="store_true", help="rebuild a stage that already exists")
    args = ap.parse_args()
    p = Panel(args.config, args.index)
    return {
        "status": cmd_status, "build": cmd_build, "simulate": cmd_simulate,
        "cache": cmd_cache, "score": cmd_score, "report": cmd_report,
    }[args.stage](p, args)


if __name__ == "__main__":
    raise SystemExit(main())
