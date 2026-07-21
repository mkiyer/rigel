"""Cross-suite validation of a landscape recipe (generalization test). Evaluates gdna_landscape_recipe.py on
each cached suite (ambig, quick, ...) — the real test of whether the derived recipe transfers off the suite it
was tuned on. Usage: OMP_NUM_THREADS=1 python scripts/debug/gdna_validate_recipe.py [suite ...]"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
import gdna_explore_lib as L  # noqa: E402

exec(Path("/Users/mkiyer/proj/rigel/scripts/debug/gdna_landscape_recipe.py").read_text())  # defines recipe(s)

suites = sys.argv[1:] or ["ambig", "quick"]


def r_v25(s):
    mk = L.masks(s)
    sel = L.vpercentile(s, mk["base"] & (mk["single"] | mk["gonly"]), 25) | mk["struct_zero"]
    return L.fit_poisson(s["g_hat"][sel], s["eff"][sel])


for suite in suites:
    try:
        scen = L.load_scenarios(suite)
    except FileNotFoundError:
        print(f"{suite:8s}: cache not built yet")
        continue
    r = L.evaluate(recipe, scen)  # noqa: F821  (recipe defined via exec above)
    b = L.evaluate(r_v25, scen)
    print(f"{suite:8s} ({len(scen)} scen)  RECIPE mean_emd={r['mean_emd']:.3f} max_emd={r['max_emd']:.3f} "
          f"mean_enr={r['mean_enr']:.3f} max_enr={r['max_enr']:.3f}   | v25 mean_emd={b['mean_emd']:.3f} "
          f"max_emd={b['max_emd']:.3f}")
    worst = sorted(r["per_group_emd"].items(), key=lambda kv: -kv[1])[:4]
    print("          worst groups: " + "  ".join(f"{'/'.join(g)}={v:.2f}" for g, v in worst))
