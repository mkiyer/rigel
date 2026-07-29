"""THE ABLATION CAMPAIGN — is the solver over-engineered, and which terms interact?

One term off at a time, plus the 2x2 factorial cell that answers the interaction question the owner asked:
if removing A and removing B each cost X and Y, does removing BOTH cost X+Y (independent), less
(REDUNDANT — they were double-paying for the same failure), or more (SYNERGISTIC — each needs the other)?

⚠ EVERY ARM HERE CHANGES RESULTS BY CONSTRUCTION. This is a modelling audit, not a bit-identity refactor;
`verify_clean.sh` must fail on all of them.

Two mechanisms, because not every term is a function:
  * MONKEYPATCH — for terms that are module-level calls (`graft_premise_logvar`, `composition_logvar`).
  * SOURCE      — for terms written inline (P1e, the level fuse). Run with `--source-arm NAME` AFTER
                  applying that arm's patch to src/ by hand; the driver only labels the rows.

The per-condition ORACLE is the expensive part (it reads sim_oracle.bam), so it is cached to disk and every
arm re-uses it — the whole campaign then costs about one bench run.

    python scratchpad/ablate_campaign.py --arms base,noP1d,noM5c
    python scratchpad/ablate_campaign.py --source-arm noP1e --arms noP1e,noP1e_noP1d
"""
from __future__ import annotations

import argparse
import dataclasses
import importlib
import os
import sys
from contextlib import contextmanager
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

import rigel.calibration.bp_solver as bps  # noqa: E402
from rigel.calibration.node_chain import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
OUT = Path("/tmp/rigel_ablate")
ORACLE = OUT / "_oracle_cache"
_EPS = 1e-9


# ── the monkeypatch arms ───────────────────────────────────────────────────────────────────────────────
@contextmanager
def _nothing():
    yield


@contextmanager
def _no_p1d():
    """P1d `ω_graft` OFF — the pooled premise log-variance goes to 0, so the graft's RNA claim is no longer
    damped by the seam-pair disagreement. `_seam_pair` still runs; only the fitted scalar is neutralised."""
    orig = bps.graft_premise_logvar
    bps.graft_premise_logvar = lambda *a, **k: (orig(*a, **k)[0], 0.0)
    try:
        yield
    finally:
        bps.graft_premise_logvar = orig


@contextmanager
def _no_m5_comp():
    """M5's COMPOSITION half OFF — `Var(log ρ_tot)` keeps only the counting term `1/n`. Verified already
    identically 0 at every AMBIG node (f_g_ref = 1 ⇒ f_g(1−f_g) = 0), so this only bites single-strand."""
    orig = bps.composition_logvar
    bps.composition_logvar = lambda f_g, E_g, E_r, var_fg, n: orig(
        f_g, E_g, E_r, np.zeros_like(np.asarray(var_fg, np.float64)), n
    )
    try:
        yield
    finally:
        bps.composition_logvar = orig


@contextmanager
def _no_p1d_src():
    with _no_p1d():
        yield


ARMS = {
    "base": _nothing,
    "noP1d": _no_p1d,
    "noM5c": _no_m5_comp,
    # source-patched arms: the patch is applied to src/ by hand; these only add the monkeypatch on top
    "noP1e": _nothing,
    "noP1e_noP1d": _no_p1d_src,
    "nofuse": _nothing,
}

ap = argparse.ArgumentParser()
ap.add_argument("--arms", default="base")
ap.add_argument("--source-arm", default=None, help="label only — assert the source patch is applied")
a = ap.parse_args()
arms = [x for x in a.arms.split(",") if x]
assert all(x in ARMS for x in arms), (arms, list(ARMS))

OUT.mkdir(exist_ok=True)
ORACLE.mkdir(exist_ok=True)
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
rows: dict[tuple[str, int], dict[str, list]] = {(m, r): {} for m in arms for r in (0, 1)}

for ci, cond in enumerate(conds):
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    ocache = ORACLE / f"{cond}.npz"
    chain = None
    for arm in arms:
        for refit in (0, 1):
            dbg: dict = {}
            cc = dataclasses.replace(cfg.calibration, calib_refit_iters=refit)
            with ARMS[arm]():
                calmod.calibrate(
                    inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                    np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
                )
            chain, cap = dbg["chain"], dbg["capture"]
            if not ocache.exists():
                Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
                np.savez_compressed(ocache, G=Gp + Gn, R=Rp + Rn)
            o = np.load(ocache)
            G, R = o["G"], o["R"]
            fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
            rt = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
            ri = np.clip(np.asarray(chain.ref_idx, np.int64), 0, rt.shape[0] - 1)
            is_reg = np.asarray(chain.kind) == REGION
            cls = np.where(~is_reg, "boundary", np.where(rt[ri] == 2, "exon",
                                                        np.where(rt[ri] == 1, "intron", "other")))
            fp, fn = np.asarray(cap["free_pos"], bool), np.asarray(cap["free_neg"], bool)
            d = rows[(arm, refit)]
            for k, v in (("fg", np.asarray(cap["f_g"])), ("fo", fo),
                         ("mass", np.asarray(cap["mass_global"])), ("var", np.asarray(cap["var_g"])),
                         ("cls", cls), ("amb", fp & fn), ("reg", is_reg), ("ss", fp ^ fn),
                         ("solv", np.asarray(cap["solvable"], bool))):
                d.setdefault(k, []).append(v)
    print(f"  [{ci + 1:>2}/{len(conds)}] {cond}", flush=True)

for (arm, refit), d in rows.items():
    p = OUT / f"{arm}_r{refit}.npz"
    np.savez_compressed(p, **{k: np.concatenate(v) for k, v in d.items()})
    print(f"wrote {p}")
