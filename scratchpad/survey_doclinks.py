"""Survey 6: every .md filename mentioned inside a LIVE calibration doc — live / archive-only / missing."""

import pathlib
import re
from collections import defaultdict

ROOT = pathlib.Path("/Users/mkiyer/proj/rigel")
D = ROOT / "docs/calibration"

LIVE = [p for p in sorted(D.glob("*.md"))]
allmd = {p.name: p for p in ROOT.glob("docs/**/*.md")}
live_names = {p.name for p in ROOT.glob("docs/**/*.md") if "archive" not in p.parts}
arch_names = {p.name for p in ROOT.glob("docs/**/*.md") if "archive" in p.parts}

pat = re.compile(r"[A-Za-z0-9_\-]+\.md")
bad = defaultdict(lambda: defaultdict(list))
for p in LIVE:
    txt = p.read_text()
    for m in pat.finditer(txt):
        n = m.group(0)
        if n == p.name:
            continue
        ln = txt[: m.start()].count("\n") + 1
        if n in live_names:
            continue
        st = "ARCHIVE-ONLY" if n in arch_names else "MISSING"
        bad[p.name][f"{st}|{n}"].append(ln)

for f in sorted(bad):
    print(f"### {f}")
    for k in sorted(bad[f]):
        st, n = k.split("|")
        lns = bad[f][k]
        print(f"    {st:13s} {n:55s} lines {lns[:8]}{'...' if len(lns) > 8 else ''}")
