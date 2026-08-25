#!/usr/bin/env python
"""⭐⭐⭐ WHERE DOES THE ERROR COME FROM? — one scenario, every object, every stage, and an ATTRIBUTION.

`toy_dissect.py` prints the solver's three-rung ladder. This goes further: it walks the error from the
per-object initialisation, through each relay hop, into the combine, and finally ABLATES each psi channel
one at a time to say which one is responsible.

**The five sections, in the order the solver computes them:**

1. **EVERY REGION AND EVERY BOUNDARY**, with truth in FRAGMENTS beside the answer, ranked by error MASS — not
   by error rate, because a small rate on a big object is what the deliverable actually pays for.
2. **THE FOUR INIT SOURCES**, per object: is there any own evidence at all (structural certainty, the
   intron factory, the strand tilt, the length channel), and what did the message-free self-solve make of
   it? ⭐ If the error is already here it is an initialisation defect and no message caused it.
3. **THE RELAY, HOP BY HOP** in both directions: the reframe ``r``, whether the composition licence was
   granted, the gDNA density each hop delivers and at what precision, and the running level. ⭐ This is
   the "how does it percolate" step — the level is visible arriving.
4. **THE COMBINE AND psi** at the worst object: every channel's mode and precision as psi receives them.
5. ⭐⭐ **THE ABLATION**: re-run psi with one channel removed at a time. ⛔ It first replays with EVERY
   channel and requires the shipped ``f_g`` back to ~1e-6 — a replay that cannot reproduce the answer is
   not measuring the answer. `region_sweep` captures ``global_lp`` / ``solve_grid`` / ``fg_init`` precisely
   so this replay is possible without re-solving the chain.

⚠ **It ablates, it does not fix.** Removing a channel tells you what that channel contributed; it is not
a proposal, and a channel whose removal helps here may be load-bearing elsewhere (TRAPS: honesty-metrics-reward-ignorance — read
an accuracy number, not an honesty number, when judging an ablation).
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

DESIGN = Path(__file__).resolve().parent
_s = importlib.util.spec_from_file_location("toy_ceiling", DESIGN / "toy_ceiling.py")
TC = importlib.util.module_from_spec(_s)
sys.modules["toy_ceiling"] = TC
_s.loader.exec_module(TC)
TH = TC.TH

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.calibration.simplex_logodds import _solve_regions_logodds_all  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import index_derived_inputs  # noqa: E402

SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"
TYPES = {0: "intergenic", 1: "intron", 2: "exon"}


def _labels(chain, ra):
    kind, obj = np.asarray(chain.kind), np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    starts, sizes = np.asarray(ra.start, np.int64), np.asarray(ra.region_size_bp, np.int64)
    out = []
    for s in range(int(chain.n_slots)):
        i = int(obj[s])
        if kind[s] == REGION:
            out.append(f"{TYPES[int(rtype[i])]}[{starts[i]:,},{starts[i] + sizes[i]:,})")
        else:
            hi, lo = s + 1, s - 1
            b = int(rtype[obj[hi]]) if hi < int(chain.n_slots) and kind[hi] == REGION else -1
            a = int(rtype[obj[lo]]) if lo >= 0 and kind[lo] == REGION else -1
            pair = "|".join(TYPES.get(x, "?") for x in sorted((a, b)) if x >= 0)
            out.append(f"{pair}@{starts[obj[hi]]:,}" if b >= 0 else f"boundary#{i}")
    return out


def run(cond, *, n_rna, genome_length, nrna, work_dir, messages=True):
    index = TranscriptIndex.load(str(INDEX))
    cfg = TH.with_messages(
        dataclasses.replace(CalibrationConfig(), calib_refit_iters=0), messages
    )
    donor = TH.harvest(SUITE / cond, index, config=cfg)
    spec = TH.SPECS["spliced_exons"]
    if genome_length:
        spec = dataclasses.replace(spec, genome_length=genome_length)
    spec = dataclasses.replace(spec, n_rna_fragments=n_rna, nrna_abundance=nrna,
                               name=f"tr_{cond[5:12]}_{n_rna}")
    sub = TC.simulate(spec, donor, work_dir)
    dbg: dict = {}
    calibrate(payload=sub.payload, strand_model=sub.strand_model, gdna_fl_pmf=donor.gdna_fl_pmf,
              rna_fl_pmf=donor.rna_fl_pmf, config=cfg, injected_priors=donor.priors, _debug=dbg,
              **index_derived_inputs(sub.index))
    return donor, spec, sub, dbg, cfg


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--n-rna", type=int, required=True)
    ap.add_argument("--genome-length", type=int, default=0)
    ap.add_argument("--nrna", type=float, default=0.0)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_trace"))
    # ⭐ DEFAULTS TO `on`, AND THAT IS THE HONEST DEFAULT HERE: sections 3, 4, 5 and 6 ARE the relay
    # (hop by hop, what ψ receives, the channel ablation, and withholding the composition licence in
    # both twins). Under the shipped mute there are no hops to walk and no channels to ablate, so this
    # instrument RUNS SECTION 1 AND THEN REFUSES rather than printing zeros for five sections.
    TH.add_messages_flag(ap, default=True)
    args = ap.parse_args()
    messages = TH.messages_on(args)

    donor, spec, sub, dbg, cfg = run(args.condition, n_rna=args.n_rna,
                                    genome_length=args.genome_length, nrna=args.nrna,
                                    work_dir=args.work_dir, messages=messages)
    cap, chain = dbg["capture"], dbg["chain"]
    # ⛔ `_uni` is written only by `messages/relay.py`. Sections 2–6 all read it (or a head-only
    # `_uni_static` key such as `og` / `pg_own` / `struct_lock`), so the refusal is deferred to just
    # after section 1 — the one section that is policy-independent.
    st = TH.relay_static(cap)
    uni = TH.relay_channels(cap)
    ra = sub.region_arrays
    lab = _labels(chain, ra)
    kind, obj = np.asarray(chain.kind), np.asarray(chain.obj_idx, np.int64)
    n_slots = int(chain.n_slots)

    ov = sub.truth.override_masses(ra)
    T = {"region": (np.asarray(ov["mass_gdna_region"], float), np.asarray(ov["mass_rna_region"], float)),
         "boundary": (np.asarray(ov["mass_gdna_boundary"], float), np.asarray(ov["mass_rna_boundary"], float))}
    # the shipped per-object answer, in the same currency as the truth
    res_g = {"region": np.asarray(cap["f_g"], float), "boundary": np.asarray(cap["f_g"], float)}

    print("=" * 132)
    print(f"⭐⭐⭐ ERROR TRACE — {args.condition}")
    print(f"   toy {spec.genome_length:,} bp · {args.n_rna:,} RNA fragments · "
          f"{sub.n_gdna_target:,} gDNA fragments · nascent {args.nrna:g}")
    print(f"   kappa {donor.priors.rna_sense_frac:.6f} · strand specificity "
          f"{donor.strand_specificity:g} · capture {'ON' if donor.capture_on else 'off'} · "
          f"pass-0 (prior-free)")
    print(TH.messages_stamp(messages))
    print("=" * 132)

    # ── 1. every object ──────────────────────────────────────────────────────────────────────────
    print("\n── 1. EVERY REGION AND EVERY BOUNDARY, ranked by error MASS (fragments of gDNA mis-assigned) ──")
    print(f"\n   {'slot':<28} {'n':>6} {'gDNA':>6} {'RNA':>6} {'true f_g':>9} {'f_g':>8} "
          f"{'Δf_g':>8} {'err frags':>10} {'share':>7}")
    print("   " + "-" * 118)
    recs = []
    for s in range(n_slots):
        ax = "region" if kind[s] == REGION else "boundary"
        i = int(obj[s])
        g, r = T[ax][0][i], T[ax][1][i]
        tot = g + r
        n = float(np.asarray(cap["count"], float).sum(axis=1)[s])
        fg = float(res_g[ax][s])
        truth = g / tot if tot > 0 else float("nan")
        err = abs(fg * tot - g) if tot > 0 else 0.0
        recs.append({"s": s, "lab": lab[s], "n": n, "g": g, "r": r, "truth": truth, "fg": fg,
                     "err": err, "mass": tot})
    tot_err = sum(x["err"] for x in recs) or 1.0
    for x in sorted(recs, key=lambda z: -z["err"]):
        d = x["fg"] - x["truth"]
        print(f"   {x['lab']:<28} {x['n']:>6.0f} {x['g']:>6.0f} {x['r']:>6.0f} "
              f"{x['truth']:>9.4f} {x['fg']:>8.4f} {d:>+8.4f} {x['err']:>10.1f} "
              f"{x['err'] / tot_err:>6.1%}")
    worst = max(recs, key=lambda z: z["err"])
    print(f"\n   ⭐ HIGHEST ERROR MASS: {worst['lab']}  "
          f"({worst['err']:.1f} of {tot_err:.1f} fragments, {worst['err'] / tot_err:.1%})")

    # ⛔⛔ EVERYTHING BELOW IS A MESSAGE-LAYER MEASUREMENT, so refuse rather than print zeros for five
    # sections (the `ladder_arm_ab.py` contract: an arm the policy cannot express is REFUSED up front).
    # ⚠ Section 1 above is policy-independent and has already printed, which is why the refusal is here
    # and not at argparse time.
    TH.require_relay(cap, what="sections 2–6 of the error trace (init sources, the relay hop by hop, "
                               "what ψ receives, the channel ablation, the licence arm)")

    # ── 2. the four init sources ─────────────────────────────────────────────────────────────────
    print("\n── 2. THE FOUR INIT SOURCES — does the object have ANY own evidence? ───────────────────")
    print("   struct = structural certainty (pure-gDNA object) · tau_own = the λ-axis evidence it earned")
    print("   ⭐ fg_strand is the strand likelihood ALONE; fg_loc adds the intron factory and ψ's reference")
    print("   ⭐⭐ `prec_+` is the SUPPLY conjunct of the composition licence: a source may may_share_composition a")
    print("      composition only if it supplied BOTH components at nonzero precision. A source whose")
    print("      RNA precision is 0 cannot may_share_composition, and its gDNA LEVEL then crosses UNSCALED.")
    print(f"\n   {'slot':<28} {'u+':>6} {'u−':>6} {'struct':>7} {'tau_own':>10} {'fg_strand':>10} "
          f"{'fg_loc':>8} {'own rho_g':>10} {'prec_g':>9} {'prec_+':>9} {'|Δ loc|':>8}")
    print("   " + "-" * 118)
    cnt = np.asarray(cap["count"], float)
    for x in recs:
        s = x["s"]
        dl = abs(cap["fg_loc"][s] - x["truth"]) if np.isfinite(x["truth"]) else float("nan")
        print(f"   {x['lab']:<28} {cnt[s, 0]:>6.0f} {cnt[s, 1]:>6.0f} "
              f"{str(bool(st['struct_lock'][s])):>7} {st['tau_own'][s]:>10.3g} "
              f"{cap['fg_strand'][s]:>10.4f} {cap['fg_loc'][s]:>8.4f} {st['og'][s]:>10.4g} "
              f"{st['pg_own'][s]:>9.3g} {st['pp_own'][s]:>9.3g} {dl:>8.4f}")

    # ── 3. the relay, hop by hop ─────────────────────────────────────────────────────────────────
    print("\n── 3. THE RELAY, HOP BY HOP — how the gDNA level travels, and where it is rescaled ─────")
    print("   `may_share_composition` = the composition licence. When it is FALSE the level crosses UNSCALED and the")
    print("   mass rescale is off; when TRUE the source's SHARE is imputed onto the destination's total.")
    for name, pin in zip(("forward L→R", "backward R→L"), cap["_rescale"]):
        src = np.asarray(pin["src"], np.int64)
        val = np.asarray(pin["valid"], bool)
        print(f"\n   {name}")
        print(f"      {'dst slot':<28} {'← src':<28} {'r':>8} {'may_share_composition':>6} {'tg delivered':>13} "
              f"{'tpg':>9} {'SPLICE IN J':>8} {'spl_prec':>9}")
        for s in range(n_slots):
            if not val[s]:
                continue
            print(f"      {lab[s]:<28} {lab[int(src[s])]:<28} {float(pin['r'][s]):>8.4f} "
                  f"{str(bool(pin['may_share_composition'][s])):>6} {float(pin['tg'][s]):>13.5g} "
                  f"{float(pin['tpg'][s]):>9.3g} "
                  f"{float(pin['spl_p'][s]) + float(pin['spl_n'][s]):>8.0f} "
                  f"{float(pin['spl_prec'][s]):>9.3g}")
    print("\n   the running level that ARRIVES at each slot (own rho_g → fwd → bwd):")
    print(f"      {'slot':<28} {'own rho_g':>11} {'fwd_g':>11} {'bwd_g':>11} {'fused cg':>11} "
          f"{'TRUE rho_g':>11}")
    for x in recs:
        s = x["s"]
        Eg = st["E_g"][s]
        true_rho = x["g"] / Eg if Eg > 0 else float("nan")
        print(f"      {x['lab']:<28} {st['og'][s]:>11.5g} {st['fwd_g'][s]:>11.5g} "
              f"{st['bwd_g'][s]:>11.5g} {uni['cg'][s]:>11.5g} {true_rho:>11.5g}")

    # ── 4. the combine, at the worst object ──────────────────────────────────────────────────────
    ws = worst["s"]
    print(f"\n── 4. WHAT ψ RECEIVES AT {worst['lab']} ───────────────────────────────────────────────")
    print(f"      truth f_g {worst['truth']:.4f}   shipped f_g {worst['fg']:.4f}   "
          f"self-solve f_g {cap['fg_loc'][ws]:.4f}")
    for nm, mo, pr in (("gDNA measurement (mo_g, cm_g)", uni["mo_g"][ws], uni["cm_g"][ws]),
                       ("RNA + measurement (mo_p, cm_p)", uni["mo_p"][ws], uni["cm_p"][ws]),
                       ("RNA − measurement (mo_n, cm_n)", uni["mo_n"][ws], uni["cm_n"][ws]),
                       ("composition λ (lam_msg, c_tau)", uni["lam_msg"][ws], uni["c_tau"][ws])):
        print(f"      {nm:<34} mode {float(mo):>10.4f}   precision {float(pr):>10.4g}"
              + ("   ⛔ DEAD" if float(pr) <= 0 else ""))

    # ── 5. the ablation ──────────────────────────────────────────────────────────────────────────
    print("\n── 5. ⭐⭐ ABLATION — re-run ψ at every slot with ONE channel removed ──────────────────")
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    n_slot = cnt.sum(axis=1)
    base_kw = dict(
        u_pos=cnt[:, 0], u_neg=cnt[:, 1], allow_pos=fp, allow_neg=fn, mass_unspl=n_slot,
        mass_spliced=np.asarray(cap["spliced"], float),
        kappa=float(donor.priors.rna_sense_frac),
        od_g=float(donor.priors.gdna_strand_overdispersion or 0.0),
        od_r=float(donor.priors.rna_strand_overdispersion or 0.0),
        n_grid=int(cfg.sweep_n_grid), L=float(cfg.sweep_logodds_window),
        n_tilt=cfg.sweep_n_tilt, n_grid_ss=cfg.sweep_n_grid_single_strand,
        # ⛔⛔ A SECOND, INDEPENDENT DEATH IN THIS FILE, and it is NOT the `_uni` one: this arm called
        # `_solve_regions_logodds_all(global_logprior=…, length_loglik=None)` and NEITHER parameter
        # exists any more. `global_logprior` became `priors` — one `CompositionPriors` bundling the gDNA
        # arm and the RNA arm (the reference LOCATION was deleted 2026-08-24), which is exactly what
        # `sweep.py` publishes as `global_lp`. `length_loglik` was PURGED with the fragment-length composition channel
        # (`b7ed7a0b`) and has no successor, so it is simply gone rather than renamed.
        # ⚠ `TRAPS: a-green-suite-hid-five-dead-instruments`: the import gate cannot see a stale kwarg
        # inside a function body.
        priors=cap["global_lp"], lam_logprior=cap["intron_prior"],
        fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"],
    )
    ch = dict(gdna_imp_mode=uni["mo_g"], gdna_imp_prec=uni["cm_g"],
              rna_imp_mode=(uni["mo_p"], uni["mo_n"]), rna_imp_prec=(uni["cm_p"], uni["cm_n"]),
              lam_imp_mode=uni["lam_msg"], lam_imp_prec=uni["c_tau"],
              theta_imp_mode=None, theta_imp_prec=None)

    def solve(**over):
        kw = dict(base_kw)
        kw.update(ch)
        kw.update(over)
        return np.asarray(_solve_regions_logodds_all(**kw).gdna_frac, float)

    solvable = np.asarray(cap["solvable"], bool)
    init = np.asarray(cap["fg_init"], float)

    def as_shipped(f):
        """⛔ `region_sweep` writes back ONLY solvable slots; a G1 object keeps its pinned init of 1.0.
        Scoring psi's raw output at those slots reads a structurally-certain object as 0 % gDNA and
        swamps the table — the first run of this file did exactly that and reported 218 error
        fragments where the gene carries 55."""
        return np.where(solvable, f, init)

    full = as_shipped(solve())
    resid = float(np.max(np.abs(full[solvable] - np.asarray(cap["f_g"], float)[solvable])))
    print(f"      ⛔ FIDELITY CHECK: replaying with every channel reproduces the shipped f_g to "
          f"{resid:.2e}")
    if resid > 1e-6:
        print("      ⛔⛔ THE REPLAY DOES NOT REPRODUCE THE ANSWER — the ablations below are not")
        print("          measuring the shipped solve. Do not read them.")
    OFF_G = dict(gdna_imp_mode=None, gdna_imp_prec=None)
    OFF_R = dict(rna_imp_mode=None, rna_imp_prec=None)
    OFF_L = dict(lam_imp_mode=None, lam_imp_prec=None)
    arms = [
        ("full (the shipped solve)", {}),
        ("− the gDNA measurement", OFF_G),
        ("− the RNA measurement", OFF_R),
        ("− the composition λ message", OFF_L),
        ("− the intron factory", dict(lam_logprior=None)),
        ("⭐ ONLY the gDNA measurement", {**OFF_R, **OFF_L}),
        ("⭐ ONLY the RNA measurement", {**OFF_G, **OFF_L}),
        ("⭐ ONLY the λ message", {**OFF_G, **OFF_R}),
        ("− ALL messages (= the self-solve)", {**OFF_G, **OFF_R, **OFF_L}),
    ]
    keep = [x for x in recs if np.isfinite(x["truth"]) and x["mass"] > 0]
    print(f"\n      {'arm':<36}" + "".join(f"{x['lab'][:13]:>15}" for x in keep))
    print(f"      {'truth':<36}" + "".join(f"{x['truth']:>15.4f}" for x in keep))
    print("      " + "-" * (36 + 15 * len(keep)))
    for nm, over in arms:
        f = as_shipped(solve(**over))
        print(f"      {nm:<36}" + "".join(f"{f[x['s']]:>15.4f}" for x in keep))
    print(f"\n      {'ERROR MASS over the whole gene (fragments of gDNA mis-assigned)':<52}")
    for nm, over in arms:
        f = as_shipped(solve(**over))
        e = sum(abs(f[x["s"]] * x["mass"] - x["g"]) for x in keep)
        print(f"      {nm:<52}{e:>10.1f}")
    print("\n      ⚠ Ablating a psi CHANNEL cannot test the reframe that BUILT the channel's value. The")
    print("        delivered level is set in the relay, upstream of anything psi can be asked to ignore.")

    # ── 6. the RELAY-level test the channel ablations cannot do ──────────────────────────────────
    print("\n── 6. ⭐⭐⭐ THE RELAY-LEVEL ARM — withhold the COMPOSITION LICENCE everywhere ────────────")
    print("   `_may_share = pop & (the source supplied both components)`. When it is TRUE the reframe")
    print("   multiplies the gDNA LEVEL by `r = rho_tot(dst)/rho_tot(src)`; when FALSE the level crosses")
    print("   UNSCALED. Forcing the POPULATION conjunct to fail everywhere makes `_may_share` false on every")
    print("   hop, in BOTH twins at once, and is the only way to ask whether the licensed reframe is")
    print("   what inflates the delivered level. ⚠ A DIAGNOSTIC, not a proposal: the licence is")
    print("   load-bearing elsewhere (`EQUATIONS.md` §3.5b/§3.5c).")
    # ⛔⛔ A THIRD, INDEPENDENT DEATH IN THIS FILE — the same species as the two above and again invisible
    # to an import gate. This patched `rigel.calibration.sweep.terminus_flank_gain`, a name that is NOT
    # THERE: the licence's only consumer is `messages/relay.py`, which binds the function at import time
    # (`from ..region_geometry import … terminus_flank_gain`). ⭐ Patch the CONSUMER's binding, because
    # patching `region_geometry` would be too late for a name already bound.
    from rigel.calibration.messages import relay as BP

    orig = BP.terminus_flank_gain

    def all_gain(boundary_flags):
        rgain, lgain = orig(boundary_flags)
        return (np.ones_like(np.asarray(rgain, bool)),
                np.ones_like(np.asarray(lgain, bool)))

    BP.terminus_flank_gain = all_gain
    try:
        dbg2: dict = {}
        calibrate(payload=sub.payload, strand_model=sub.strand_model,
                  gdna_fl_pmf=donor.gdna_fl_pmf, rna_fl_pmf=donor.rna_fl_pmf, config=cfg,
                  injected_priors=donor.priors, _debug=dbg2, **index_derived_inputs(sub.index))
    finally:
        BP.terminus_flank_gain = orig
    cap2 = dbg2["capture"]
    uni2 = TH.require_relay(cap2, what="the licence-withheld twin arm")
    print(f"\n      {'slot':<28} {'TRUE rho_g':>11} {'cg SHIPPED':>11} {'cg no-licence':>14} "
          f"{'f_g SHIPPED':>12} {'f_g no-lic':>11} {'truth':>8}")
    for x in recs:
        s = x["s"]
        Eg = st["E_g"][s]
        print(f"      {x['lab']:<28} {(x['g'] / Eg if Eg > 0 else float('nan')):>11.5g} "
              f"{uni['cg'][s]:>11.5g} {uni2['cg'][s]:>14.5g} {x['fg']:>12.4f} "
              f"{float(cap2['f_g'][s]):>11.4f} {x['truth']:>8.4f}")
    e1 = sum(abs(x["fg"] * x["mass"] - x["g"]) for x in keep)
    e2 = sum(abs(float(cap2["f_g"][x["s"]]) * x["mass"] - x["g"]) for x in keep)
    e0 = sum(abs(as_shipped(solve(**{**OFF_G, **OFF_R, **OFF_L}))[x["s"]] * x["mass"] - x["g"])
             for x in keep)
    print(f"\n      error mass  SHIPPED {e1:.1f}   ·   licence WITHHELD {e2:.1f}   ·   "
          f"no messages at all {e0:.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
