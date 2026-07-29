"""Pass-0 vs oracle error diagnostic — WHERE is the error, and is it ATTACKABLE by precision?

Single-arm production pass-0 over all 32 cached ambig_dense_10mb scenarios. Three questions:
  (1) Rank scenarios by mass-weighted |Δf_g| vs oracle, split by node type (boundary / intron / exon).
  (2) Decompose the error budget: CONFIDENTLY-WRONG (var<1 & |Δf_g|>0.3 — a precision over-credit we could
      attack) vs UNDER-DETERMINED (var≥1 — honest high-variance, an identifiability floor) vs small-error.
  (3) The over-crediting test: are the confidently-wrong nodes disproportionately SHORT (small eff-length) or
      HIGH-SPLICED (over-credited spliced measurement)? If not, precision is not the lever.
Run in the rigel env, OMP_NUM_THREADS=1.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth, _true_fg  # noqa: E402
from rigel.calibration.bp_solver import REGION, node_sweep, node_global_geometry  # noqa: E402
from rigel.calibration.calibrate import _build_intron_prior  # noqa: E402
from rigel.calibration.effective_length import region_eff_length  # noqa: E402
from rigel.calibration.node_chain import build_node_chain  # noqa: E402
from rigel.calibration.node_geometry import build_node_geometry, build_node_statics, init_beliefs  # noqa: E402
from rigel.calibration.npmle import DensityNPMLE  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.calibration.strand_balance import fit_strand_balance  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb"); cache = suite / "_selfsolve_cache"
index = TranscriptIndex.load(str(suite / "rigel_index")); cfg = PipelineConfig(); cc = cfg.calibration
ra = RegionArrays.from_index(index)
conds = sorted(p.stem for p in cache.glob("*.pkl"))


def solve(inp):
    """Production pass-0. Returns per-node dict: f_g, tf(true), mass, var, ntype(0=bnd,1=intron,2=exon,-1),
    eff(gDNA support eff-length), spl(spliced mass)."""
    pl = inp["payload"]; sub = CalibrationSubstrate.from_payload(pl, ra); bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    st = build_node_statics(chain, sub, bsub, ra)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    bal = fit_strand_balance(inp["strand_model"]); kappa = float(bal.rna_sense_frac)
    inter = coarse_type_array(np.asarray(ra.signature)) == 0
    gp = float(np.sum(np.asarray(sub.contained.n_unspliced_pos, float)[inter]))
    gn = float(np.sum(np.asarray(sub.contained.n_unspliced_neg, float)[inter]))
    mass_g, eff_g = [np.asarray(x, float) for x in node_global_geometry(chain, geom)]
    ep = DensityNPMLE.fit(mass_g, eff_g, bandwidth=cc.npmle_bandwidth)
    reff = region_eff_length(ra.region_size_bp, inp["gdna_fl_pmf"])
    ip = _build_intron_prior(chain, sub, ra, reff, cc)
    be = init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa, n_grid=cc.sweep_n_grid,
                      n_grid_ss=cc.sweep_n_grid_single_strand, logodds_window=cc.sweep_logodds_window, statics=st)
    fb = node_sweep(chain, st, geom, be, ra, rna_sense_frac=kappa,
                    rna_strand_overdispersion=0.0011, gdna_strand_overdispersion=0.0012,
                    n_gdna_obs=gp + gn, n_rna_obs=float(bal.n_observations), n_grid=cc.sweep_n_grid,
                    n_grid_ss=cc.sweep_n_grid_single_strand, logodds_window=cc.sweep_logodds_window,
                    n_tilt=cc.sweep_n_tilt, gdna_prior=None, intron_prior=ip)
    idx = np.asarray(chain.ref_idx, np.int64); is_r = np.asarray(chain.kind) == REGION
    tf_r, m_r = _true_fg(inp["region_pools"]); tf_b, m_b = _true_fg(inp["boundary_pools"])

    def ci(a):
        return np.clip(idx, 0, a.shape[0] - 1)

    tf = np.where(is_r, tf_r[ci(tf_r)], tf_b[ci(tf_b)])
    mass = np.where(is_r, m_r[ci(m_r)], m_b[ci(m_b)])
    rt = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    ntype = np.where(is_r, rt[ci(rt)], -1).astype(np.int64)  # region: 0=intergenic 1=intron 2=exon; bnd=-1
    spl = (np.asarray(geom.spliced_pos_left) + np.asarray(geom.spliced_pos_right)
           + np.asarray(geom.spliced_neg_left) + np.asarray(geom.spliced_neg_right))
    return dict(f=np.asarray(fb.f_g, float), tf=tf, mass=mass, var=np.asarray(fb.var_gdna, float),
                is_r=is_r, ntype=ntype, eff=eff_g, spl=spl,
                ok=np.isfinite(tf) & (mass > _EPS) & np.isfinite(np.asarray(fb.var_gdna, float)))


def werr(d, sel):
    m = d["mass"][sel]
    return float((m * np.abs(d["f"][sel] - d["tf"][sel])).sum() / max(m.sum(), _EPS)) if m.sum() > 0 else np.nan


# ---- gather all scenarios ----
rows = []
allnodes = {k: [] for k in ("f", "tf", "mass", "var", "eff", "spl", "ntype", "is_bnd")}
for c in conds:
    inp = _scan_and_truth(suite, c, index, cfg, Path("/tmp/rigel_selfsolve"), cache)
    d = solve(inp)
    ok = d["ok"]
    bnd = ok & ~d["is_r"]; intron = ok & d["is_r"] & (d["ntype"] == 1); exon = ok & d["is_r"] & (d["ntype"] == 2)
    err_all = werr(d, ok)
    rows.append((c, err_all, werr(d, bnd), werr(d, intron), werr(d, exon),
                 d["mass"][ok & (d["var"] < 1.0) & (np.abs(d["f"] - d["tf"]) > 0.3)].sum()
                 / max(d["mass"][ok].sum(), _EPS) * 100))
    for k, key in (("f", "f"), ("tf", "tf"), ("mass", "mass"), ("var", "var"), ("eff", "eff"),
                   ("spl", "spl"), ("ntype", "ntype")):
        allnodes[k].append(d[key][ok])
    allnodes["is_bnd"].append((~d["is_r"])[ok])

A = {k: np.concatenate(v) for k, v in allnodes.items()}

# ---- PART 1: ranked scenario table ----
print("=" * 92)
print("PART 1 — pass-0 |Δf_g| vs oracle, ranked WORST-first (mass-weighted). CW% = confidently-wrong mass.")
print(f"{'scenario':44s} {'ALL':>6s} {'bnd':>6s} {'intr':>6s} {'exon':>6s} {'CW%':>6s}")
for c, ea, eb, ei, ee, cw in sorted(rows, key=lambda r: -r[1]):
    lbl = c.replace("gdna_", "").replace("_nrna", " ").replace("_capture", " ")
    print(f"{lbl:44s} {ea:6.3f} {eb:6.3f} {ei:6.3f} {ee:6.3f} {cw:6.1f}")

# ---- PART 2: error-budget decomposition (whole suite) ----
print("\n" + "=" * 92)
print("PART 2 — error BUDGET: of the total mass-weighted error, how much is ATTACKABLE (confident-wrong) vs")
print("         UNDER-DETERMINED (honest high-variance = identifiability floor)?")
err = np.abs(A["f"] - A["tf"]); m = A["mass"]; tot_err_mass = float((m * err).sum())
big = err > 0.3
conf_wrong = big & (A["var"] < 1.0)          # confident AND wrong  -> a precision over-credit (attackable)
under_det = big & (A["var"] >= 1.0)          # wrong but honest high-variance -> identifiability
small = ~big
for name, sel in (("confident-wrong (var<1, |Δ|>0.3)", conf_wrong),
                  ("under-determined (var≥1, |Δ|>0.3)", under_det),
                  ("small error (|Δ|≤0.3)", small)):
    em = float((m[sel] * err[sel]).sum())
    print(f"  {name:38s}: {em / max(tot_err_mass, _EPS) * 100:5.1f}% of total error mass   "
          f"({m[sel].sum() / m.sum() * 100:5.1f}% of node mass)")

# ---- PART 3: the over-crediting test ----
print("\n" + "=" * 92)
print("PART 3 — over-crediting test: are CONFIDENT-WRONG nodes SHORT / HIGH-SPLICED vs the population?")
print("         If they are NOT, short-node/spliced over-credit is not the driver ⇒ precision is not the lever.")


def q(x, w):  # mass-weighted median + IQR
    if w.sum() <= 0:
        return (np.nan, np.nan, np.nan)
    o = np.argsort(x); xs, ws = x[o], w[o]; cw = np.cumsum(ws) / ws.sum()
    return tuple(float(np.interp(p, cw, xs)) for p in (0.25, 0.5, 0.75))


for label, arr in (("eff-length (gDNA support)", A["eff"]), ("spliced mass", A["spl"]),
                   ("total mass", A["mass"])):
    cwq = q(arr[conf_wrong], m[conf_wrong]); allq = q(arr, m)
    print(f"  {label:26s}  confident-wrong q25/med/q75 = {cwq[0]:9.1f}/{cwq[1]:9.1f}/{cwq[2]:9.1f}"
          f"   | ALL = {allq[0]:9.1f}/{allq[1]:9.1f}/{allq[2]:9.1f}")
print(f"\n  confident-wrong nodes: {int(conf_wrong.sum())} of {conf_wrong.size} "
      f"({m[conf_wrong].sum()/m.sum()*100:.2f}% of mass). fraction that are BOUNDARIES: "
      f"{m[conf_wrong & A['is_bnd']].sum()/max(m[conf_wrong].sum(),_EPS)*100:.0f}%  "
      f"(boundaries are {m[A['is_bnd']].sum()/m.sum()*100:.0f}% of all mass)")
print(f"  fraction of confident-wrong mass with spliced>0: "
      f"{m[conf_wrong & (A['spl']>_EPS)].sum()/max(m[conf_wrong].sum(),_EPS)*100:.0f}%  "
      f"(vs {m[A['spl']>_EPS].sum()/m.sum()*100:.0f}% of all mass has spliced>0)")
