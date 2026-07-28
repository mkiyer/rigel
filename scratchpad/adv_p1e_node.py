"""Verify the section-4.3 worked-node SOLVED values (2651 / 2441 / 1405) with and without RIGEL_P1E."""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import p1e_lib as L  # noqa: E402

inp, dbg = L.solve(L.CONDS["gdna300_ss0.50_present_capOFF"])
nf = L.node_frame(inp, dbg)
tag = "P1E=" + (os.environ.get("RIGEL_P1E") or "OFF")
print(f"  {tag:<12}{'node':>7}{'oracle':>10}{'solved':>10}{'sd':>10}{'z':>9}")
for j in (2651, 2441, 1405):
    sd = float(np.sqrt(max(nf["var_g"][j], 0.0)))
    z = abs(nf["fg"][j] - nf["fo"][j]) / max(sd, 1e-12)
    print(f"  {'':<12}{j:>7}{nf['fo'][j]:>10.4f}{nf['fg'][j]:>10.4f}{sd:>10.4f}{z:>9.2f}")
