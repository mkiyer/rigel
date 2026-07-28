"""Survey 5: run every scripts/debug/*.py with --help (short timeout) and classify failures."""

import pathlib
import subprocess
import sys

ROOT = pathlib.Path("/Users/mkiyer/proj/rigel")
scripts = sorted(ROOT.glob("scripts/debug/*.py"))
bad = []
for p in scripts:
    try:
        r = subprocess.run(
            [sys.executable, str(p), "--help"],
            capture_output=True, text=True, timeout=90, cwd=str(ROOT),
        )
    except subprocess.TimeoutExpired:
        bad.append((p.name, "TIMEOUT(90s) on --help"))
        continue
    err = r.stderr
    # A script that does not use argparse will just run; that is fine unless it raises an
    # import/name/attribute error, which is what we care about.
    for marker in ("ImportError", "ModuleNotFoundError", "NameError", "AttributeError",
                   "SyntaxError", "TypeError", "ValueError", "KeyError", "IndexError",
                   "FileNotFoundError", "RuntimeError", "AssertionError"):
        if f"{marker}:" in err:
            last = [ln for ln in err.strip().splitlines() if ln.startswith(marker)]
            bad.append((p.name, (last[-1] if last else marker)[:170]))
            break
    else:
        if r.returncode not in (0, 1, 2) and err.strip():
            bad.append((p.name, f"rc={r.returncode} " + err.strip().splitlines()[-1][:150]))

for n, e in bad:
    print(f"{n}\t{e}")
print(f"\n{len(bad)} / {len(scripts)} failed", file=sys.stderr)
