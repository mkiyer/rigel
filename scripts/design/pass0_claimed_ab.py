"""HOW WELL DOES PASS-0 SOLVE THE SLOTS IT CLAIMS, PER POLICY? — silent / relay / fanout at the
stage-0 substrate's claimed populations, judged on NOTHING else (the owner's rule).

Two claimed populations, each scored as misplaced gDNA fragments ``Σ|est_gdna − true_gdna|`` against
certified `slot_truth`, per condition and per policy, DELIVER/REFUTE split and never pooled
(`TRAPS: a-refutability-test-needs-the-refuting-channel-in-the-fixture`):

    B   ``ss_intron_boundary``  (stage 3's destinations), on the boundary axis
    E   ``solvable_exon``       (stage 4's destinations), on the region axis

DELIVER rows are the slots whose certified truth is EXACTLY pure gDNA (the claim is true and must be
delivered); REFUTE rows carry real RNA (the claim must be overturned by evidence). ⛔ A pooled number
over the whole library cannot judge pass-0: outside its claimed slots the fan-out is silence BY
DESIGN, so a whole-library comparison charges it for slots it deliberately leaves unsolved — that
comparison lives in `ladder_arm_ab.py --arm policy_fanout`, as context.

Each policy runs the FULL pipeline in-process off the cached scan and oracle
(`pass0_vs_oracle.measure_condition`), so the three columns differ by exactly the message policy.
``--self-test`` falsifies the SCORER: injected error at a claimed slot must be caught in the right
population and split, and error parked on an unclaimed slot must be scored NOWHERE.

⭐⭐ ``--dissect`` — WHY DOES THE FAN-OUT LAND WHERE IT DOES AT THE CLAIMED EXONS? The per-slot
survey that produced the 2026-08-23 residual verdict, so that verdict is re-derivable rather than
quoted. Four readings off ONE cached payload, all on the exon population:

    STAGE     silent / fanout / fanout-with-the-exon-λ-muted / fanout-with-the-one-sided-bound-muted,
              which says WHICH stage of the fan-out moved the slots at all
    LEDGER    per slot: the exon's own solve (``f_g``, ``tau_lam``), the REAL delivered message
              (``lam_mode``, ``lam_prec``), the posterior, and certified truth — BIAS and PRECISION
              separately, because a claim can be wrong, over-confident, or both
    CALIBRATION  the delivered message's ``z = (lam_mode − lam_true)·sqrt(lam_prec)`` against
              certified truth: a message whose precision is honest has ``median |z| ≈ 0.67``
    DRIFT     the flank's CERTIFIED composition pushed through :func:`~rigel.calibration.messages.fanout.flank_to_exon_lambda`
              and scored against the exon's certified composition — the frame error ALONE, with the
              gDNA route and the merged RNA route reported separately

⛔ The message is RECORDED from the policy's own ``deliver``, never re-implemented, and an IDENTITY
GATE refuses the survey unless the recorded call belongs to the sweep the debug capture describes.
⛔ ``--dissect`` reads ``calib_refit_iters = 0`` as well as the shipped count, because a policy that
changes pass-0 changes the FITTED prior — a local-versus-global attribution needs both.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np

_EPS = 1.0e-12
_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_SUITE = _RUNS / "suite" / "ladder"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import build_region_statics  # noqa: E402
from rigel.calibration.splice_graph import build_boundary_flags_array  # noqa: E402
from rigel.calibration.structural_claims import build_structural_claims  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import read_scan_cache  # noqa: E402

#: the three policies, as configs — the SHIPPED config with exactly one axis varied.
POLICIES = {
    "silent": lambda: CalibrationConfig(message_propagation=False),
    "relay": lambda: CalibrationConfig(),
    "fanout": lambda: CalibrationConfig(message_policy="fanout"),
}


def claimed_masks(chain, claims, truth: dict) -> dict:
    """The two claimed populations on their OWN axes, plus the DELIVER/REFUTE split from certified
    truth (pure ⇔ every RNA column exactly zero at the slot)."""
    out = {}
    for tag, kind_val, slot_mask in (
        ("B", BOUNDARY, np.asarray(claims.ss_intron_boundary, bool)),
        ("E", REGION, np.asarray(claims.solvable_exon, bool)),
    ):
        is_k = np.asarray(truth["kind"]) == kind_val
        obj = np.asarray(truth["obj"], np.int64)[is_k]
        n_axis = int(obj.max()) + 1 if obj.size else 0
        true_g = np.zeros(n_axis)
        true_g[obj] = np.asarray(truth["n_gdna"], np.float64)[is_k]
        rna = np.zeros(n_axis)
        rna[obj] = (
            np.asarray(truth["n_nrna"], np.float64) + np.asarray(truth["n_mrna"], np.float64)
        )[is_k]
        mask = np.zeros(n_axis, bool)
        slot_is_k = np.asarray(chain.kind) == kind_val
        mask[np.asarray(chain.obj_idx, np.int64)[slot_is_k]] = slot_mask[slot_is_k]
        out[tag] = {
            "mask": mask,
            "true_gdna": true_g,
            "deliver": mask & (rna == 0.0),
            "refute": mask & (rna > 0.0),
        }
    return out


def score(est_gdna_by_axis: dict, masks: dict) -> dict:
    """``Σ|est − true|`` per population and split — the estimate arrays come per axis
    (``mass_gdna_boundary`` for B, ``mass_gdna_region`` for E)."""
    rows = {}
    for tag, m in masks.items():
        err = np.abs(np.asarray(est_gdna_by_axis[tag], np.float64) - m["true_gdna"])
        rows[tag] = {
            "claimed": float(err[m["mask"]].sum()),
            "deliver": float(err[m["deliver"]].sum()),
            "refute": float(err[m["refute"]].sum()),
        }
    return rows


def audit(index, region_arrays, suite: Path, condition: str, work_dir: Path, policies) -> dict:
    P0 = _sibling("pass0_vs_oracle.py")
    t = dict(np.load(Path(suite) / "oracle_cache" / condition / "slot_truth.npz"))
    cache = read_scan_cache(Path(suite) / "oracle_cache" / condition / "_main", index)
    chain = build_region_chain(cache.payload.ref_region_offsets, cache.payload.ref_boundary_offsets)
    statics = build_region_statics(chain, region_arrays, build_boundary_flags_array(index))
    masks = claimed_masks(chain, build_structural_claims(chain, statics), t)
    out = {}
    for name in policies:
        m = P0.measure_condition(
            bam=str(Path(suite) / condition / "sim_oracle.bam"),
            index=index,
            pipeline_config=PipelineConfig(),
            calibration_config=POLICIES[name](),
            work_dir=work_dir,
            tag=condition,
            truth_pmfs=None,
            oracle_cache=Path(suite) / "oracle_cache",
        )
        res = m.arms["final"] if "final" in m.arms else m.arms[sorted(m.arms)[0]]
        out[name] = score({"B": res.mass_gdna_boundary, "E": res.mass_gdna_region}, masks)
    return out


def transport_composition(f_src, U, S, eff_sj, eff_gx, eff_rx, eff_gc, eff_rc):
    """The transfer's VALUE alone — :func:`flank_to_exon_lambda`'s route-merge ∘ reframe with the
    source composition supplied rather than solved, as a composition in the exon's frame.

    Two structural identities make it falsifiable with no data: with ``S = 0`` and matching frames it
    is the IDENTITY (a pure reparameterisation returns ``f_src`` exactly), and adding spliced flux
    moves the answer monotonically toward RNA. Used by the DRIFT reading, where ``f_src`` is the
    flank's CERTIFIED composition, so whatever error survives is the FRAME's and not the source's.
    """
    f = np.clip(np.asarray(f_src, np.float64), 0.0, 1.0)
    U = np.asarray(U, np.float64)
    rho_g = f * U / np.maximum(np.asarray(eff_gx, np.float64), _EPS)
    rho_nu = (1.0 - f) * U / np.maximum(np.asarray(eff_rx, np.float64), _EPS)
    S = np.asarray(S, np.float64)
    e_sj = np.asarray(eff_sj, np.float64)
    rho_mu = np.where((S > 0.0) & (e_sj > 0.0), S / np.maximum(e_sj, _EPS), 0.0)
    g = rho_g * np.asarray(eff_gc, np.float64)
    r = (rho_nu + rho_mu) * np.asarray(eff_rc, np.float64)
    return np.where(g + r > 0.0, g / np.maximum(g + r, _EPS), 0.0)


def z_calibration(mode, prec, true_value, live):
    """``median |z|`` of a delivered claim against certified truth, ``z = (mode − true)·sqrt(prec)``.

    ⭐ The reference is **0.6745** — the median ``|z|`` of an honestly-calibrated Gaussian — so the
    variance is over-claimed by ``(median|z| / 0.6745)²``. ⛔ This is the only reading that judges a
    message's PRECISION against anything; ``lam_prec`` compared with the own side's says which voice
    wins the fuse, never whether either is entitled to be heard.
    """
    m = np.asarray(live, bool) & (np.asarray(prec, np.float64) > 0.0)
    if not m.any():
        return float("nan"), 0
    z = (np.asarray(mode, np.float64)[m] - np.asarray(true_value, np.float64)[m]) * np.sqrt(
        np.asarray(prec, np.float64)[m]
    )
    return float(np.median(np.abs(z))), int(m.sum())


class _DeliverRecorder:
    """Records the fan-out's OWN ``deliver`` output and every ``flank_to_exon_lambda`` argument.

    ⛔ Wrapping rather than re-implementing is the point: a survey built from a second copy of the
    transfer would agree with itself and with nothing else."""

    def __init__(self):
        self.calls: list = []
        self._pending: list = []

    def __enter__(self):
        from rigel.calibration.messages import fanout as FO

        self._FO = FO
        self._real = (FO.FanOutPolicy.prepare, FO.flank_to_exon_lambda)
        rec = self

        def transfer(*a, **k):
            out = self._real[1](*a, **k)
            rec._pending.append((a, out))
            return out

        def prepare(policy, ctx):
            relay = rec._real[0](policy, ctx)
            real_deliver = relay.deliver

            def deliver(left, right):
                rec._pending = []
                msg = real_deliver(left, right)
                rec.calls.append(
                    {"ctx": ctx, "msg": msg, "left": left, "right": right,
                     "transfer": list(rec._pending)}
                )
                return msg

            relay.deliver = deliver
            return relay

        FO.FanOutPolicy.prepare = prepare
        FO.flank_to_exon_lambda = transfer
        return self

    def __exit__(self, *exc):
        self._FO.FanOutPolicy.prepare, self._FO.flank_to_exon_lambda = self._real
        return False


def _mute_exon_lambda(is_region: np.ndarray, one_sided_too: bool):
    """A policy wrapper that drops one delivered channel at REGION slots — the STAGE attribution.

    ⛔ It patches the delivered MESSAGE, never the solver, so an arm differs from the shipped fan-out
    by exactly the channel named and nothing else."""
    from rigel.calibration.messages import PsiMessage
    from rigel.calibration.messages import fanout as FO

    real = FO.FanOutPolicy.prepare

    def prepare(policy, ctx):
        relay = real(policy, ctx)
        real_deliver = relay.deliver

        def deliver(left, right):
            m = real_deliver(left, right)
            if one_sided_too:
                return PsiMessage(lam_mode=m.lam_mode, lam_prec=m.lam_prec)
            return PsiMessage(
                lam_mode=np.where(is_region, 0.0, m.lam_mode),
                lam_prec=np.where(is_region, 0.0, m.lam_prec),
            )

        relay.deliver = deliver
        return relay

    FO.FanOutPolicy.prepare = prepare
    return real


def dissect(index, region_arrays, suite: Path, condition: str, sweep: bool = False) -> None:
    """The four readings of ``--dissect``, off ONE cached payload. See the module docstring."""
    from rigel.calibration.calibrate import calibrate
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.messages import fanout as FO
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import build_sj_geometry_arrays

    truth = dict(np.load(Path(suite) / "oracle_cache" / condition / "slot_truth.npz"))
    cache = read_scan_cache(Path(suite) / "oracle_cache" / condition / "_main", index)
    payload, strand_model = cache.payload, cache.strand_model
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    statics = build_region_statics(chain, region_arrays, build_boundary_flags_array(index))
    claims = build_structural_claims(chain, statics)
    is_region = np.asarray(chain.kind) == REGION
    obj = np.asarray(chain.obj_idx, np.int64)
    max_size = int(payload.max_length)
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, max_size),
        gdna_opportunity=gdna_opportunity_from_index(index, max_size),
    )
    kwargs = dict(
        region_arrays=region_arrays,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        sj=build_sj_geometry_arrays(index),
        boundary_flags=build_boundary_flags_array(index),
    )
    cnt = np.asarray(truth["count"], np.float64)
    ng = np.asarray(truth["n_gdna"], np.float64)
    tf = np.asarray(truth["true_f_g"], np.float64)
    exon = np.asarray(claims.solvable_exon, bool) & is_region

    def slot_fg(res):
        fg = np.zeros(int(chain.n_slots))
        fg[is_region] = np.asarray(res.gdna_frac_region, np.float64)[obj[is_region]]
        fg[~is_region] = np.asarray(res.gdna_frac_boundary, np.float64)[obj[~is_region]]
        return fg

    print(f"\n== {condition}   DISSECTION at the {int(exon.sum())} claimed exons")

    # ── ① STAGE — which part of the fan-out moved these slots at all ──────────────────────────────
    print(f"\n   ① STAGE   {'arm':<20}{'refits=0':>12}{'refits=shipped':>16}")
    shipped_refits = CalibrationConfig().calib_refit_iters
    stage: dict = {}
    for name in ("silent", "fanout", "fanout_no_exon_lam", "fanout_no_one_sided"):
        row = []
        for refits in (0, shipped_refits):
            if name == "silent":
                cfg = CalibrationConfig(message_propagation=False, calib_refit_iters=refits)
                restore = None
            else:
                cfg = CalibrationConfig(message_policy="fanout", calib_refit_iters=refits)
                restore = (
                    _mute_exon_lambda(is_region, name == "fanout_no_one_sided")
                    if name != "fanout"
                    else None
                )
            try:
                res = calibrate(payload=payload, config=cfg, **kwargs)
            finally:
                if restore is not None:
                    FO.FanOutPolicy.prepare = restore
            fg = slot_fg(res)
            # ⛔ the BASIS gate: the per-slot reconstruction must be the scored array itself
            proj = np.asarray(res.mass_gdna_region, np.float64)
            recon = np.zeros_like(proj)
            recon[obj[is_region]] = (fg * cnt)[is_region]
            if not np.allclose(recon, proj, rtol=1e-9, atol=1e-6):
                raise SystemExit("basis gate FAILED: f_g·count is not mass_gdna_region")
            row.append(float(np.abs(fg * cnt - ng)[exon].sum()))
            stage[(name, refits)] = fg
        print(f"             {name:<20}{row[0]:>12.0f}{row[1]:>16.0f}")
    print(
        "             ⭐ 'fanout_no_exon_lam' equal to 'silent' ⇒ stage 3 is inert at exons;\n"
        "                a gap between refits=0 and shipped ⇒ the FITTED PRIOR is the competitor."
    )

    # ── ② LEDGER + ③ CALIBRATION + ④ DRIFT, off one recorded fan-out sweep ───────────────────────
    with _DeliverRecorder() as rec:
        res = calibrate(
            payload=payload,
            config=CalibrationConfig(message_policy="fanout", calib_refit_iters=0),
            _debug=(dbg := {}),
            **kwargs,
        )
    call = rec.calls[-1]
    ctx, msg = call["ctx"], call["msg"]
    if not np.allclose(np.asarray(ctx.own.f_g), np.asarray(dbg["capture"]["fg_loc"])):
        raise SystemExit("identity gate FAILED: the recorded deliver is not the captured sweep")
    print("\n   identity gate ✔ the recorded deliver IS the captured sweep")

    own_fg = np.asarray(ctx.own.f_g, np.float64)
    tau_own = np.asarray(ctx.own.tau_lam, np.float64)
    lam_m = np.asarray(msg.lam_mode, np.float64)
    lam_p = np.asarray(msg.lam_prec, np.float64)
    fg_fan = slot_fg(res)
    fg_sil = stage[("silent", 0)]
    err_s = np.abs(fg_sil * cnt - ng)
    err_f = np.abs(fg_fan * cnt - ng)
    d = err_f - err_s
    heard = exon & (lam_p > 0.0)
    print(
        f"\n   ② LEDGER  {int(heard.sum())} of {int(exon.sum())} claimed exons hear a lambda claim; "
        f"they carry {d[heard].sum():+.0f} of the {d[exon].sum():+.0f} error change"
    )
    w = heard & (cnt > 0.0)
    tot = max(float(cnt[w].sum()), _EPS)
    f_msg = 1.0 / (1.0 + np.exp(-lam_m))

    def wmae(x):
        return float(np.sum(cnt[w] * np.abs(x[w] - tf[w])) / tot)

    def wsig(x):
        return float(np.sum(cnt[w] * (x[w] - tf[w])) / tot)

    print(f"             mass-weighted |f_g − true|  own {wmae(own_fg):.4f}   message {wmae(f_msg):.4f}")
    print(f"             mass-weighted SIGNED        own {wsig(own_fg):+.4f}   message {wsig(f_msg):+.4f}")
    fin = w & (tau_own > 0.0)
    ratio = lam_p[fin] / np.maximum(tau_own[fin], _EPS)
    print(
        f"             lam_prec / tau_lam  median {np.median(ratio) if ratio.size else float('nan'):.3g}"
        f"   >1 at {float((ratio > 1).mean()) if ratio.size else float('nan'):.1%} of slots"
        f"   (own tau EXACTLY 0 at {int((w & (tau_own <= 0)).sum())} — kappa=1/2 kills the strand term)"
    )

    lam_true = np.log(np.clip(tf, _EPS, 1 - _EPS)) - np.log1p(-np.clip(tf, _EPS, 1 - _EPS))
    interior = heard & (tf > 0.0) & (tf < 1.0)
    mz, nz = z_calibration(lam_m, lam_p, lam_true, interior)
    print(
        f"\n   ③ CALIBRATION  median |z| {mz:.3f} over {nz} interior-truth exons  ⇒  the delivered "
        f"variance is over-claimed {(mz / 0.6745) ** 2:.1f}x  (an honest claim reads 0.674)"
    )

    print(
        "\n   ④ DRIFT   the flank's composition pushed through the transfer, SOLVED vs CERTIFIED.\n"
        "             ⭐ If the CERTIFIED source delivers a WORSE value than the solved one, the\n"
        "                source's own error is CANCELLING the frame's, and the frame is the fault."
    )
    lam_b_of = {}
    for _side, (_a, _r) in zip(("L", "R"), call["transfer"]):
        lam_b_of[_side] = np.asarray(_a[0], np.float64)
    g = ctx.geometry
    eg = np.asarray(ctx.eff_gdna, np.float64)
    er = np.asarray(ctx.eff_rna, np.float64)
    U = np.asarray(ctx.n_slot, np.float64)
    col = np.where(np.asarray(ctx.free_pos, bool), 0, 1)
    for side, flank, nb, S_all, E_all in (
        ("L", claims.exon_flank_left, call["left"], g.sj_count_hi, g.eff_sj_hi),
        ("R", claims.exon_flank_right, call["right"], g.sj_count_lo, g.eff_sj_lo),
    ):
        src = np.asarray(nb.src, np.int64)
        # ⛔ restricted to the exons that actually HEARD a claim, so ④ describes the same population
        #   ② and ③ do — an unheard licensed flank is not part of the verdict.
        act = (
            heard
            & np.asarray(flank, bool)
            & np.asarray(nb.valid, bool)
            & (tf > 0.0)
            & (cnt > 0.0)
        )
        if not act.any():
            continue
        S = np.asarray(S_all, np.float64)[src, col]
        E_sj = np.asarray(E_all, np.float64)[src, col]
        got = transport_composition(tf[src], U[src], S, E_sj, eg[src], er[src], eg, er)
        f_solved = 1.0 / (1.0 + np.exp(-np.clip(lam_b_of[side], -ctx.logodds_window,
                                                ctx.logodds_window)))
        got_s = transport_composition(f_solved, U[src], S, E_sj, eg[src], er[src], eg, er)
        ww = cnt[act]
        rho_g_x = tf[src] * U[src] / np.maximum(eg[src], _EPS)
        rho_g_true = ng / np.maximum(eg, _EPS)
        ok = act & (rho_g_x > 0) & (rho_g_true > 0)
        lg = np.log(rho_g_x[ok]) - np.log(rho_g_true[ok])
        print(
            f"      side {side}: {int(act.sum())} RNA-carrying claims — delivered f_g − exon truth:"
            f"   SOLVED source {float(np.sum(ww * (got_s[act] - tf[act])) / ww.sum()):+.4f}"
            f"   CERTIFIED source {float(np.sum(ww * (got[act] - tf[act])) / ww.sum()):+.4f}"
        )
        print(
            f"         the gDNA route arrives at log ratio "
            f"{float(np.sum(cnt[ok] * lg) / max(cnt[ok].sum(), _EPS)):+.3f} "
            f"(median {np.median(lg):+.3f}) against the exon's own certified gDNA density"
        )
    print(
        "\n   ⛔ Read ② and ③ together: ② says whether the claim is RIGHT, ③ whether it is entitled\n"
        "      to the confidence it carries, and ④ whether the FRAME rather than the source is at fault."
    )
    if sweep:
        _sweep_readings(
            payload, kwargs, chain, claims, is_region, obj, cnt, ng, exon, suite, condition,
            index, region_arrays, float(np.abs(fg_sil * cnt - ng)[exon].sum()), heard,
        )


def _sweep_readings(
    payload, kwargs, chain, claims, is_region, obj, cnt, ng, exon, suite, condition,
    index, region_arrays, silent_floor, heard,
):
    """⑤ and ④b — what each candidate DIRECTION would have to overcome. ⛔ Both price a hypothetical
    rather than describing the tool, so they are opt-in (``--sweep``) and cost ~15 `calibrate` runs.

    ⑤ sweeps the delivered λ precision and, separately, a uniform shift of the delivered λ mode. If
    the fault were PRECISION, some scaling beats silence; if the claim is BIASED, the best scaling is
    the one that silences it. ④b rebuilds the transfer's frame constant
    ``K = log(E_rx·E_gc) − log(E_gx·E_rc)`` from the simulator's OWN post-capture fl pmfs and delivers
    the difference as a pure mode shift — ``K`` enters λ with coefficient exactly 1 and carries no
    variance term, so its error is a bias with no confession. ⛔ A DIAGNOSTIC: a truth pmf does not
    exist at run time; this prices the constant's error, it does not propose reading truth.
    """
    from rigel.calibration.calibrate import calibrate
    from rigel.calibration.messages import PsiMessage
    from rigel.calibration.messages import fanout as FO
    from rigel.calibration.region_geometry import build_region_geometry
    from rigel.calibration.splice_graph import build_sj_geometry_arrays
    from rigel.calibration.substrate import CalibrationSubstrate

    real = FO.FanOutPolicy.prepare

    def run(scale: float, shift, ):
        def prepare(policy, ctx):
            relay = real(policy, ctx)
            rd = relay.deliver
            L = float(ctx.logodds_window)

            def deliver(left, right):
                m = rd(left, right)
                lp = np.asarray(m.lam_prec, np.float64)
                lm = np.asarray(m.lam_mode, np.float64)
                b = shift if np.ndim(shift) else float(shift)
                return PsiMessage(
                    lam_mode=np.where(is_region & (lp > 0.0), np.clip(lm + b, -L, L), lm),
                    lam_prec=np.where(is_region, lp * scale, lp),
                    rna_mode=m.rna_mode,
                    rna_prec=m.rna_prec,
                    rna_one_sided=m.rna_one_sided,
                )

            relay.deliver = deliver
            return relay

        FO.FanOutPolicy.prepare = prepare
        try:
            res = calibrate(
                payload=payload,
                config=CalibrationConfig(message_policy="fanout", calib_refit_iters=0),
                **kwargs,
            )
        finally:
            FO.FanOutPolicy.prepare = real
        fg = np.zeros(int(chain.n_slots))
        fg[is_region] = np.asarray(res.gdna_frac_region, np.float64)[obj[is_region]]
        return float(np.abs(fg * cnt - ng)[exon].sum())

    print(f"\n   ⑤ SWEEP   claimed-exon error; the silent floor is {silent_floor:.0f}")
    print("      lam_prec x  " + "".join(f"{c:>10g}" for c in (1.0, 0.3, 0.1, 0.03, 0.01, 0.0)))
    print(
        "                  "
        + "".join(f"{run(c, 0.0):>10.0f}" for c in (1.0, 0.3, 0.1, 0.03, 0.01, 0.0))
    )
    print("      lam_mode +  " + "".join(f"{b:>10g}" for b in (0.0, 0.5, 0.8, 1.0, 1.5)))
    print(
        "                  "
        + "".join(f"{run(1.0, b):>10.0f}" for b in (0.0, 0.5, 0.8, 1.0, 1.5))
    )
    print(
        "      ⭐ if some scaling BEATS the floor the fault is precision; if the best scaling is the\n"
        "         one that silences the claim, the claim itself is wrong (a-variance-cannot-fix-a-bias)."
    )

    # ── ④b the frame constant, rebuilt from the simulator's own pmfs ──
    P0 = _sibling("pass0_vs_oracle.py")
    cond_dir = Path(suite) / condition
    max_size = int(payload.max_length)
    sub = CalibrationSubstrate.from_payload(payload, region_arrays)
    sj_arrays = build_sj_geometry_arrays(index)

    def frames(gp, rp):
        g = build_region_geometry(chain, sub, region_arrays, sj_arrays, gp, rp)
        return np.asarray(g.eff_gdna, np.float64), np.asarray(g.eff_rna, np.float64)

    eg_f, er_f = frames(kwargs["gdna_fl_pmf"], kwargs["rna_fl_pmf"])
    g_true = P0.truth_length_pmf(cond_dir, "gdna", max_size)
    r_true = P0.truth_length_pmf(cond_dir, "rna", max_size)
    eg_t, er_t = frames(g_true, r_true)
    left_idx = np.clip(np.asarray(chain.left, np.int64), 0, int(chain.n_slots) - 1)
    right_idx = np.clip(np.asarray(chain.right, np.int64), 0, int(chain.n_slots) - 1)

    def K(eg, er, src):
        return (
            np.log(np.maximum(er[src], _EPS))
            + np.log(np.maximum(eg, _EPS))
            - np.log(np.maximum(eg[src], _EPS))
            - np.log(np.maximum(er, _EPS))
        )

    dK = np.zeros(int(chain.n_slots))
    for flank, src in (
        (np.asarray(claims.exon_flank_left, bool), left_idx),
        (np.asarray(claims.exon_flank_right, bool), right_idx),
    ):
        use = exon & flank & (dK == 0.0)
        dK = np.where(use, K(eg_t, er_t, src) - K(eg_f, er_f, src), dK)
    # ⛔ reported over the exons the shift actually LANDS on (a claim is delivered there), not over
    #   every licensed one — the two differ, and the wider set understates the correction.
    m = heard & (dK != 0.0)
    got = run(1.0, dK)
    print(
        f"\n   ④b FRAME CONSTANT  fitted fl pmf means gDNA "
        f"{float((np.arange(len(kwargs['gdna_fl_pmf'])) * kwargs['gdna_fl_pmf']).sum()):.2f} / RNA "
        f"{float((np.arange(len(kwargs['rna_fl_pmf'])) * kwargs['rna_fl_pmf']).sum()):.2f}   "
        f"truth {float((np.arange(len(g_true)) * g_true).sum()):.2f} / "
        f"{float((np.arange(len(r_true)) * r_true).sum()):.2f}"
    )
    print(
        f"      median dK over {int(m.sum())} DELIVERING exons {np.median(dK[m]):+.4f} nats  ⇒  "
        f"claimed-exon error {run(1.0, 0.0):.0f} → {got:.0f} with K rebuilt from truth pmfs"
    )
    print(
        "      ⛔ K enters lambda with coefficient exactly 1 and carries NO variance term, so its\n"
        "         error is a bias that never confesses. This prices it; it does not read truth at run time."
    )


def report(condition: str, rows: dict, policies) -> None:
    print(f"\n== {condition}")
    print(f"   {'pop':<4} {'split':<9}" + "".join(f"{p:>12}" for p in policies))
    for tag in ("B", "E"):
        for split in ("claimed", "deliver", "refute"):
            print(
                f"   {tag if split == 'claimed' else '':<4} {split:<9}"
                + "".join(f"{rows[p][tag][split]:>12.0f}" for p in policies)
            )


def self_test() -> int:
    ok = 0

    def check(name, cond):
        nonlocal ok
        print(f"   {'✔' if cond else '✘'} {name}")
        if not cond:
            raise SystemExit(f"self-test FAILED at: {name}")
        ok += 1

    # a synthetic 9-slot chain: ig B exon B intron B exon B ig — one claimed boundary pair, and the
    # claims derived by the REAL builder so the masks cannot drift from the substrate's definition.
    from rigel.calibration.region_geometry import RegionStatics
    from rigel.calibration.splice_graph import FLAG_DONOR_POS, FLAG_TSS_POS

    chain = build_region_chain(np.array([0, 5]), np.array([0, 4]))
    fp = np.array([0, 0, 1, 1, 1, 1, 1, 0, 0], bool)
    mp = np.array([0, 0, 1, 0, 0, 0, 1, 0, 0], bool)
    bflags = np.zeros(9, np.uint16)
    bflags[1] = bflags[7] = FLAG_TSS_POS
    bflags[3] = bflags[5] = FLAG_DONOR_POS
    statics = RegionStatics(
        n_slots=9,
        free_pos=fp,
        free_neg=np.zeros(9, bool),
        mrna_active_pos=mp,
        mrna_active_neg=np.zeros(9, bool),
        boundary_flags=np.where(np.asarray(chain.kind) == BOUNDARY, bflags, 0).astype(np.uint16),
    )
    claims = build_structural_claims(chain, statics)
    truth = {
        "kind": np.asarray(chain.kind),
        "obj": np.asarray(chain.obj_idx),
        "n_gdna": np.array([9, 2, 1, 8, 6, 8, 1, 2, 9], float),
        "n_nrna": np.array([0, 0, 0, 0, 3, 0, 0, 0, 0], float),
        "n_mrna": np.array([0, 0, 5, 0, 0, 0, 5, 0, 0], float),
    }
    masks = claimed_masks(chain, claims, truth)
    check(
        "the claimed populations are the substrate's (2 boundaries, 2 exons)",
        masks["B"]["mask"].sum() == 2 and masks["E"]["mask"].sum() == 2,
    )
    check(
        "DELIVER and REFUTE split by certified purity and never overlap",
        masks["B"]["deliver"].sum() == 2
        and masks["E"]["refute"].sum() == 2
        and not (masks["B"]["deliver"] & masks["B"]["refute"]).any(),
    )
    # perfect estimates score zero; an injected error is caught in the right population and split
    truth_b = np.zeros(4)
    truth_b[np.asarray(chain.obj_idx)[np.asarray(chain.kind) == BOUNDARY]] = truth["n_gdna"][
        np.asarray(chain.kind) == BOUNDARY
    ]
    perfect = {"B": truth_b, "E": truth["n_gdna"][np.asarray(chain.kind) == REGION]}
    check(
        "perfect estimates score exactly zero everywhere",
        all(v == 0.0 for row in score(perfect, masks).values() for v in row.values()),
    )
    bad = {k: v.copy().astype(float) for k, v in perfect.items()}
    bad["B"][1] += 7.0  # boundary obj 1 = the claimed donor (slot 3): DELIVER side
    bad["E"][2] += 5.0  # region obj 2 = the intron — NOT a claimed exon: must score nowhere
    got = score(bad, masks)
    check(
        "an injected boundary error lands in B/deliver, exactly",
        got["B"]["claimed"] == 7.0 and got["B"]["deliver"] == 7.0 and got["B"]["refute"] == 0.0,
    )
    check(
        "error parked on an UNCLAIMED slot is scored nowhere (judge only what is claimed)",
        got["E"]["claimed"] == 0.0,
    )
    bad2 = {k: v.copy().astype(float) for k, v in perfect.items()}
    bad2["E"][2 + 0] = bad2["E"][2] + 0.0  # no-op guard
    bad2["E"][np.flatnonzero(masks["E"]["refute"])[0]] += 3.0
    got2 = score(bad2, masks)
    check(
        "an injected exon error lands in E/refute, exactly",
        got2["E"]["refute"] == 3.0 and got2["E"]["deliver"] == 0.0,
    )

    # ── --dissect's two pure readings, each falsified against something it cannot fake ──
    one = np.ones(1)
    ident = transport_composition(
        np.array([0.3]), np.array([40.0]), np.zeros(1), 100.0 * one, 200.0 * one, 200.0 * one,
        200.0 * one, 200.0 * one,
    )
    check(
        "the transfer with no spliced flux and matching frames is the IDENTITY on composition",
        np.isclose(ident[0], 0.3, atol=1e-12),
    )
    scaled = transport_composition(
        np.array([0.3]), np.array([40.0]), np.zeros(1), 100.0 * one, 200.0 * one, 200.0 * one,
        7.0 * 200.0 * one, 7.0 * 200.0 * one,
    )
    check(
        "…and it is CAPTURE-INVARIANT under one uniform pull: scaling both destination frames "
        "together leaves the composition exactly unchanged",
        np.isclose(scaled[0], ident[0], atol=1e-12),
    )
    lop = transport_composition(
        np.array([0.3]), np.array([40.0]), np.zeros(1), 100.0 * one, 200.0 * one, 200.0 * one,
        200.0 * one, 2.0 * 200.0 * one,
    )
    check(
        "⛔ …but NOT under a pull that differs between the two components — the reframe cancels a "
        "COMMON factor only, which is the premise the drift reading tests",
        lop[0] < ident[0] - 0.1,
    )
    spliced = transport_composition(
        np.array([0.3]), np.array([40.0]), np.array([50.0]), 100.0 * one, 200.0 * one, 200.0 * one,
        200.0 * one, 200.0 * one,
    )
    # ⭐ pinned to a HAND-DERIVED value, not merely to a direction: with every frame equal to E the
    #   delivered composition is ``f·U / (U + S·E/E_sj)``, so ``f=0.3, U=40, S=50, E=200, E_sj=100``
    #   gives ``12/140 = 3/35`` exactly. ⛔ A direction-only check cannot see ``U`` at all — it
    #   cancels between the two unspliced routes — and the perturbation sweep proved that hole.
    check(
        "adding spliced flux moves the composition toward RNA, by exactly f·U/(U + S·E/E_sj)",
        spliced[0] < ident[0] and np.isclose(spliced[0], 3.0 / 35.0, atol=1e-12),
    )
    # ⭐ the SOURCE frames must be separable from each other, or a fixture with ``gx == rx`` cannot
    #   tell a dropped source frame from a kept one — the perturbation sweep found exactly that hole.
    #   Hand-derived: with ``S = 0`` and ``gc == rc`` the odds transform by exactly ``rx/gx``, so
    #   ``f = 0.3`` (odds 3/7) at ``rx = 2·gx`` must deliver odds 6/7, i.e. ``f = 6/13``.
    lopsided_src = transport_composition(
        np.array([0.3]), np.array([40.0]), np.zeros(1), 100.0 * one, 200.0 * one, 400.0 * one,
        200.0 * one, 200.0 * one,
    )
    check(
        "the SOURCE frames enter as their own ratio: rx = 2·gx transforms the odds by exactly 2",
        np.isclose(lopsided_src[0], 6.0 / 13.0, atol=1e-12),
    )
    rng = np.random.default_rng(3)
    # ⛔ a NON-ZERO, non-constant truth, or a reading that ignores its truth argument scores
    #   identically to one that uses it — the second hole the sweep found.
    tru = rng.normal(1.5, 2.0, 20000)
    honest = tru + rng.normal(0.0, 1.0 / np.sqrt(4.0), 20000)
    mz, n = z_calibration(honest, np.full(20000, 4.0), tru, np.ones(20000, bool))
    check(
        "an HONESTLY-calibrated claim reads median |z| = 0.6745 (the reading's own reference)",
        n == 20000 and abs(mz - 0.6745) < 0.02,
    )
    mz4, _ = z_calibration(tru + 2.0 * (honest - tru), np.full(20000, 4.0), tru, np.ones(20000, bool))
    check(
        "a claim whose variance is 4x over-claimed reads 2x the reference — the reading has the "
        "resolution the verdict rests on",
        abs((mz4 / 0.6745) ** 2 - 4.0) < 0.4,
    )
    mzb, _ = z_calibration(
        honest + 3.0, np.full(20000, 4.0), tru, np.ones(20000, bool)
    )
    check(
        "⛔ and a pure BIAS at honest precision also inflates |z| — so ③ alone cannot separate the "
        "two, which is why ② is read beside it",
        mzb > 3.0 * 0.6745,
    )
    check(
        "no live claim ⇒ nan and a zero count, never a silent 0.0",
        np.isnan(z_calibration(honest, np.zeros(20000), tru, np.ones(20000, bool))[0]),
    )
    print(f"self-test: {ok}/{ok}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None)
    ap.add_argument("--policies", nargs="*", default=list(POLICIES), choices=list(POLICIES))
    ap.add_argument(
        "--work-dir",
        type=Path,
        default=Path.home() / ".cache" / "rigel_pass0_claimed",
    )
    ap.add_argument(
        "--dissect",
        action="store_true",
        help="the per-slot survey behind the claimed-exon verdict (see the module docstring)",
    )
    ap.add_argument(
        "--sweep",
        action="store_true",
        help="with --dissect: add readings ⑤ and ④b, which price the candidate DIRECTIONS (~15 solves)",
    )
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    index = TranscriptIndex.load(args.index)
    region_arrays = RegionArrays.from_index(index)
    conds = (
        [args.condition]
        if args.condition
        else sorted(p.name for p in (args.suite / "scan_cache").iterdir())
    )
    for c in conds:
        if args.dissect:
            dissect(index, region_arrays, args.suite, c, sweep=args.sweep)
        else:
            report(
                c,
                audit(index, region_arrays, args.suite, c, args.work_dir, args.policies),
                args.policies,
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
