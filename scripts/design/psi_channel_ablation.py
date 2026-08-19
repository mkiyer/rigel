#!/usr/bin/env python
"""⭐⭐⭐ WHICH ψ CHANNEL IS DOING THE WORK? — the ablation, at the ψ boundary, on the REAL ladder.

Step 3 of the debug loop, one level below `worst_objects.py`. That instrument says WHICH OBJECTS carry
the error; this says WHICH CHANNEL put it there. It records every argument of the final ψ combine, then
re-solves with one imputed channel nulled at a time — so an attribution is a re-solve of the real call
rather than a second implementation of it.

⚠ **This line said "on the SHIPPED solver" until 2026-08-17, and the paragraph two below now refutes
it**: the arms exist only under the relay, so the instrument defaults to a configuration the tool does
NOT ship. A one-line summary that contradicts its own body is the failure `CLAUDE.md`'s own preamble
records — *the first sentence won* — so the claim is deleted here rather than qualified further down.

⛔⛔ **THIS IS A MESSAGE-LAYER INSTRUMENT AND IT DEFAULTS TO ``--messages on``, WHICH IS *NOT* THE
SHIPPED CONFIGURATION.** All four arms null a ``*_imp_*`` channel, and those keyword arguments are
written by ``messages/relay.py`` alone — under the shipped ``CalibrationConfig.message_propagation =
False`` (``SilentPolicy``) the solver is never handed one, so every arm is structurally inert. ⚠ Until
2026-08-17 the file had no way to express its own policy: it ran the shipped config, recorded **0**
combines carrying a non-``None`` ``lam_imp_prec``, and died on ``IndexError: list index out of range``
selecting the last of an empty list. ⭐ The contract is ``ladder_arm_ab.py``'s — the message setting is
PART OF THE ARM, it is STAMPED into the output, and a policy that cannot express the arm is REFUSED UP
FRONT rather than reported as a null result (`TRAPS: an-ablation-that-never-ran`).

⛔⛔ **ψ's ``gdna_frac`` IS THE FRACTION OF THE *UNSPLICED* BANK, SO THE SCORER'S DENOMINATOR IS THE
UNSPLICED COUNT — NEVER THE OBJECT'S TOTAL.** On the BOUNDARY axis a truth total is
``boundary_unspliced + boundary_spliced``, and a spliced crossing is certified RNA that ψ never
deconvolves. Scoring ``f_g · (unspliced + spliced)`` therefore inflated every BOUNDARY row: measured
2026-08-17 at ``g05 ss0.50 capture_on``, **2,760,340** unspliced crossings against **2,182,619**
certified-spliced ones, so the old denominator over-read the BOUNDARY axis **1.79×** while the REGION
axis — which carries no spliced molecule at all — read **1.00×**. ⭐ Both factors are printed on every
run, and the partition is ASSERTED against the truth arm's own per-object total, so the basis is a gate
rather than a claim.

⛔ **TRAPS: byte-identity-gate — the ``as-is`` arm must reproduce the run BIT-IDENTICALLY, and reproducing it means reproducing
the WRITE-BACK too.** `solve_chain` keeps the incoming belief wherever ``solvable`` is False, so comparing
raw ψ output against the stored belief differs by up to 1.0 at every unsolvable slot. The first version of
this script reported ``max |Δ| = 1.0`` and its numbers were on a different basis than the panel's; the
`wb()` helper is that fix and the identity boundary is the gate.

⭐⭐ **AND THE GATE IS NOW THE *SELECTOR*, WHICH IS WHAT MAKES ``--arm pass0`` MEAN ANYTHING.**
``pass0_vs_oracle.measure_condition`` runs the ``pass0`` arm and then the ``final`` one into ONE
recorder, so the old ``combine[-1]`` was **always the final arm's** — an ``--arm pass0`` ablation
attributed the final arm's channels to pass-0, and said so only in a printed ``False``. The combine that
produced an arm's answer is exactly the one whose re-solve reproduces that answer bit-for-bit, so the
search runs backwards and takes the first that closes. ⭐ Measured 2026-08-17 at ``g05 ss0.50
capture_on``: ``final`` closes on combine **5 of 5** and ``pass0`` on combine **2 of 5**, with the three
later ones failing — the mis-attribution was real, not hypothetical.
⛔ **TRAPS: could-the-arm-have-fired — every ablation prints how many LIVE slots it COULD have moved**, so "no effect" is separable
from "never fired".

⭐ **Read the per-slot table, not only the totals.** TRAPS: all-small-singly-large-jointly: when every single ablation is small
and the joint one is large, the channels are all built from one upstream quantity and ablating consumers
tells you nothing. Re-measured 2026-08-17 on the rebuilt ladder at ``g05 ss0.50 capture_on``,
``--arm final --messages on``: **−9.5 / −15.3 / +0.1 / −0.2 %** singly (LEVEL / certified RNA /
composition / tilt) against **−32.2 %** jointly.

⛔ **THE PERCENTAGES THIS PARAGRAPH USED TO QUOTE ARE UNREPRODUCIBLE AND ARE NOT REPLACED BY THE ONES
ABOVE.** They read *"−2.3 / −6.8 / −4.9 / −0.1 % singly against −60.7 % jointly"* at ``g01 ss0.50
capture_on``, and every leg of that citation is gone: ``g01`` was deleted with the 36-condition ladder on
2026-08-13 (the rebuilt rungs are ``g00/g05/g50/g98``), and the scorer that produced them used the
inflated denominator above. What would reproduce them is a panel carrying a ``g01`` rung; nothing on
disk does.

Usage::

    python scripts/design/psi_channel_ablation.py --condition NAME [--arm pass0|final]
    python scripts/design/psi_channel_ablation.py --messages off      # REFUSED, with the reason
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_DESIGN = Path("/Users/mkiyer/proj/rigel/scripts/design")


def _sibling(name: str):
    key = name[:-3]
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, _DESIGN / name)
    m = importlib.util.module_from_spec(spec)
    sys.modules[key] = m
    spec.loader.exec_module(m)
    return m


P0 = _sibling("pass0_vs_oracle.py")
#: ⭐ The relay contract has ONE home and this is it — `toy_harness.py`'s `--messages` block, which
#: `reframe_walk.py` / `toy_trace_error.py` / `zero_controls.py` / `toy_dissect.py` / `toy_ceiling.py`
#: / `certified_rna_audit.py` already share. ⛔ Re-declaring the flag, the shipped-default lookup and
#: the stamp here would be a second home for one rule, which is how two instruments come to disagree
#: about which configuration produced their numbers. ⚠ This file runs on the REAL ladder and uses none
#: of the toy machinery; the import is for the contract alone.
TH = _sibling("toy_harness.py")

import rigel.calibration.region_init as NI  # noqa: E402
import rigel.calibration.simplex_logodds as SL  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"
CACHE = SUITE / "oracle_cache"
TYPE_NAME = {0: "intergenic", 1: "intron", 2: "exon"}

CALLS: list[dict] = []
_orig = SL._solve_regions_logodds_all


def _rec(*a, **kw):
    CALLS.append({"args": a, "kw": dict(kw)})
    return _orig(*a, **kw)


def main() -> int:
    ap = argparse.ArgumentParser()
    # ⚠ `g01` until 2026-08-13; the rebuilt ladder's low-gDNA rung is `g05` (the one that satisfies
    # `suite_resolves.py`'s requirement (f), 0 < rate <= 0.10).
    ap.add_argument("--condition", default="gdna_g05_ss_0.50_nrna_none_capture_on")
    ap.add_argument("--suite", type=Path, default=SUITE)
    ap.add_argument("--index", type=Path, default=INDEX)
    ap.add_argument("--oracle-cache", type=Path, default=CACHE)
    ap.add_argument("--arm", default="final", choices=("pass0", "final"))
    # ⭐ DEFAULTS TO `on`, AND THAT IS THE ONLY HONEST DEFAULT HERE: every arm below nulls a channel the
    # relay is the sole writer of, so under the shipped mute this instrument has nothing to ablate.
    TH.add_messages_flag(ap, default=True)
    args = ap.parse_args()
    messages = TH.messages_on(args)

    # ⛔ REFUSED UP FRONT, before the scan — `ladder_arm_ab._require`'s contract. An arm the policy
    # cannot express must not be run and reported as "no effect" (`TRAPS: an-ablation-that-never-ran`).
    if not messages:
        raise SystemExit(
            "⛔ every arm of psi_channel_ablation nulls a `*_imp_*` channel, and `messages/relay.py` is "
            "the ONLY writer of those\n"
            "   keyword arguments. With message_propagation=False (SilentPolicy) the solver is handed "
            "none of them, so all four\n"
            "   arms are structurally inert and the ablation would report a null result it never "
            "measured.\n"
            "   Re-run with `--messages on`. ⚠ That is NOT the shipped configuration "
            f"(CalibrationConfig.message_propagation = {TH.MESSAGES_SHIPPED})\n"
            "   and the stamp on the output will say so."
        )

    import rigel.calibration.sweep as BP
    patched = []
    for mod in (SL, NI, BP):
        if hasattr(mod, "_solve_regions_logodds_all"):
            setattr(mod, "_solve_regions_logodds_all", _rec)
            patched.append(mod.__name__)
    print(f"   [TRAPS: an-ablation-that-never-ran] recording ψ in {patched}")

    index = TranscriptIndex.load(str(args.index))
    # ⚠ `measure_condition` derives its pass-0 config as `replace(calibration_config, refit=0)`, so the
    # policy threads through BOTH arms from this one field — no monkeypatching.
    m = P0.measure_condition(
        bam=str(args.suite / args.condition / "sim_oracle.bam"), index=index,
        pipeline_config=PipelineConfig(),
        calibration_config=TH.with_messages(CalibrationConfig(), messages),
        work_dir=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "eta_dissect",
        tag=args.condition, truth_pmfs=None, oracle_cache=args.oracle_cache,
    )
    print(f"\n⭐⭐ HEAD channel ablation — {args.condition}  arm={args.arm}")
    print(TH.messages_stamp(messages))
    print(f"   library f_gdna  truth {m.library_f_gdna['T']:.4f}   HEAD {m.library_f_gdna[args.arm]:.4f}")

    dbg = m.debug_final if args.arm == "final" else m.debug_pass0
    chain, cap = dbg["chain"], dbg["capture"]
    ra = dbg["region_arrays"]
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    n_slots = kind.shape[0]
    n_regions, n_boundaries = int(m.payload.n_regions), int(m.payload.n_boundaries)
    region_slot = np.full(n_regions, -1, np.int64)
    boundary_slot = np.full(n_boundaries, -1, np.int64)
    region_slot[obj[kind == REGION]] = np.flatnonzero(kind == REGION)
    boundary_slot[obj[kind != REGION]] = np.flatnonzero(kind != REGION)

    def field(res, which):
        out = np.zeros(n_slots)
        for ax, sl in (("region", region_slot), ("boundary", boundary_slot)):
            ok = sl >= 0
            out[sl[ok]] = np.asarray(getattr(res, f"mass_{which}_{ax}"), np.float64)[ok]
        return out

    tg, tr = field(m.truth, "gdna"), field(m.truth, "rna")
    pg = field(m.arms[args.arm], "gdna")
    total = tg + tr

    # ── ⛔⛔ THE SCORER'S BASIS, AND IT IS A GATE RATHER THAN A CLAIM ─────────────────────────────
    # ψ answers `gdna_frac` — the gDNA share of the UNSPLICED bank. The truth arm's per-object total is
    # `region_contained` on the region axis and `boundary_unspliced + boundary_spliced` on the boundary
    # one (`pass0_vs_oracle.check_same_basis` states it in those words), and a certified-spliced
    # crossing is RNA that gDNA cannot produce and ψ never arbitrates. So the denominator is the
    # unspliced count. ⛔ If these two banks stop summing to the truth's own total, every number below
    # is on a basis nothing else in the tree shares — so it raises rather than prints.
    unspl = np.asarray(cap["count"], np.float64).sum(axis=1)
    spl = np.asarray(cap["spliced"], np.float64)
    off = np.abs(unspl + spl - total)
    if float(off.max()) > 1e-6:
        w = int(np.argmax(off))
        raise SystemExit(
            "⛔ BASIS GATE: the solver's own (unspliced + certified-spliced) count is not the truth "
            "arm's per-object total.\n"
            f"   worst slot {w} ({'REGION' if kind[w] == REGION else 'BOUNDARY'}): "
            f"{unspl[w]:.6g} + {spl[w]:.6g} = {unspl[w] + spl[w]:.6g} against truth {total[w]:.6g}.\n"
            "   Nothing below can be scored until that is explained."
        )
    live = unspl > 0.0
    err = np.where(live, pg - tg, 0.0)
    is_b = kind != REGION
    print("   ⛔ SCORER BASIS — ψ's `gdna_frac` is a fraction of the UNSPLICED bank, so the UNSPLICED "
          "count is the denominator")
    for nm, sel in (("REGION", ~is_b), ("BOUNDARY", is_b)):
        u, s = float(unspl[sel].sum()), float(spl[sel].sum())
        print(f"      {nm:<9} unspliced {u:>14,.0f}   certified spliced {s:>14,.0f}   "
              f"⇒ the object TOTAL would over-read it {(u + s) / max(u, 1.0):.2f}×")

    combine = [c for c in CALLS if c["kw"].get("lam_imp_prec") is not None]
    if not combine:
        raise SystemExit(
            f"⛔ NOT ONE of the {len(CALLS):,} recorded `_solve_regions_logodds_all` calls carries a "
            f"non-None `lam_imp_prec`, so the\n"
            f"   relay published no imputed channel and there is nothing to ablate "
            f"(messages={'on' if messages else 'off'}).\n"
            "   ⚠ The ARTIFACT is the witness here, not the flag: if this fires with `--messages on` "
            "then the policy stopped\n"
            "   handing ψ its channels, and the ablation must NOT be reported as a null result "
            "(`TRAPS: an-ablation-that-never-ran`)."
        )
    print("   capture keys:", sorted(k for k in cap if isinstance(cap[k], np.ndarray)))
    fg = np.asarray(cap["f_g"], np.float64)
    # ⛔ the WRITE-BACK, reproduced: `solve_chain` keeps the incoming belief wherever `solvable` is False,
    # so comparing raw psi output against the stored belief differs by up to 1.0 at every unsolvable slot.
    solvable = (np.asarray(cap["free_pos"], bool) | np.asarray(cap["free_neg"], bool)) & live
    def wb(x):
        return np.where(solvable, np.clip(x, 0.0, 1.0), fg)

    def resolve(C, drop=()):
        kk = dict(C["kw"])
        for n in drop:
            kk[f"{n}_imp_mode"] = None
            kk[f"{n}_imp_prec"] = None
        return np.asarray(_orig(*C["args"], **kk).gdna_frac, np.float64)

    # ── ⛔⛔ WHICH RECORDED COMBINE IS *THIS* ARM'S? — THE IDENTITY GATE IS THE SELECTOR ──────────
    # `measure_condition` runs `pass0` and then `final` into ONE `CALLS` list, so `combine[-1]` is
    # always the FINAL arm's. Taking it unconditionally made `--arm pass0` ablate the final arm's
    # channels and report them as pass-0's — a silent mis-attribution whose only witness was a printed
    # `False`. ⭐ The combine that produced an arm's answer is exactly the one whose re-solve reproduces
    # that answer BIT-FOR-BIT, so the gate this file already owed is also the way to FIND it: search
    # backwards and take the first that closes. ⛔ No tolerance — `TRAPS: byte-identity-gate`.
    C, tried = None, []
    for cand in reversed(combine):
        d = np.abs(wb(resolve(cand)) - fg)
        tried.append(float(d.max()))
        if float(d.max()) == 0.0:
            C = cand
            break
    if C is None:
        raise SystemExit(
            f"⛔ TRAPS: byte-identity-gate — NONE of the {len(combine)} recorded combines re-solves to "
            f"the `{args.arm}` arm's own answer.\n"
            f"   best max |Δ| = {min(tried):.3e} over {len(tried)} tried, on {n_slots:,} slots. Every "
            f"ablation below would be measured\n"
            f"   against a belief the run never held, so nothing is printed. ⚠ Either the recorder no "
            f"longer sees the combine that\n"
            f"   produced this arm, or the write-back `wb()` reproduces is no longer what `solve_chain` "
            f"does."
        )
    # ⚠ The SOLVABLE COUNT is printed with the verdict because `wb()` passes `fg` straight through
    # wherever `solvable` is False: if that mask were empty, EVERY candidate would close and the
    # selector would silently degrade back to `combine[-1]` while still printing a pass. The count is
    # what separates "the gate had teeth" from "the gate could not fail" — the shape of
    # `TRAPS: an-ablation-that-never-ran`, applied to a gate rather than to an arm.
    print(f"   ⭐ TRAPS: byte-identity-gate — combine {len(combine) - len(tried) + 1} of "
          f"{len(combine)} re-solves to the `{args.arm}` arm's own answer BIT-FOR-BIT "
          f"(max |Δ| = 0.0e+00; the {len(tried) - 1} later one(s) did not)")
    print(f"      over {int(solvable.sum()):,} SOLVABLE of {n_slots:,} slots — the ones `wb()` does "
          f"NOT pass through, so the gate could have failed")

    def solve(drop=()):
        return resolve(C, drop)

    def sigma(x):
        return float(np.abs(x[live] * unspl[live] - tg[live]).sum())

    def could(drop):
        n = 0
        for nm in drop:
            p = C["kw"].get(f"{nm}_imp_prec")
            if p is None:
                continue
            for q in (p if isinstance(p, tuple) else (p,)):
                n = max(n, int(((np.asarray(q, np.float64) > 0) & live).sum()))
        return n

    print("\n   ⭐⭐ WHICH CHANNEL IS DOING THE WORK IN HEAD")
    print(f"      {'arm':<32} {'Σ|err| frag':>14} {'vs as-is':>10} {'could move':>11}")
    ref = sigma(wb(solve()))
    for name, drop in (("as-is (HEAD)", ()), ("− gdna_imp (the LEVEL)", ("gdna",)),
                       ("− rna_imp (certified RNA)", ("rna",)), ("− lam_imp (composition)", ("lam",)),
                       ("− theta_imp (tilt)", ("theta",)),
                       ("− ALL messages", ("gdna", "rna", "lam", "theta"))):
        s = sigma(wb(solve(drop)))
        print(f"      {name:<32} {s:>14,.0f} {(s - ref) / max(ref, 1) * 100:>+9.1f}% {could(drop):>11,}")

    # per-channel reach, by slot type
    rt = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    st = np.where(kind == REGION, rt[np.clip(obj, 0, rt.shape[0] - 1)], -1)
    print("\n   ⭐ WHO RECEIVES EACH CHANNEL IN HEAD (nonzero precision), by slot type")
    print(f"      {'channel':<26} {'BOUNDARY':>9} {'intergenic':>11} {'intron':>9} {'exon':>9}")
    chans = {
        "gdna_imp (LEVEL)": np.asarray(C["kw"]["gdna_imp_prec"], np.float64),
        "rna_imp (certified)": np.asarray(C["kw"]["rna_imp_prec"][0], np.float64)
        + np.asarray(C["kw"]["rna_imp_prec"][1], np.float64),
        "lam_imp (composition)": np.asarray(C["kw"]["lam_imp_prec"], np.float64),
        "theta_imp (tilt)": np.asarray(C["kw"]["theta_imp_prec"], np.float64),
    }
    for nm, p in chans.items():
        row = [int(((st == c) & (p > 0)).sum()) for c in (-1, 0, 1, 2)]
        print(f"      {nm:<26} {row[0]:>9,} {row[1]:>11,} {row[2]:>9,} {row[3]:>9,}")
    tot = [int((st == c).sum()) for c in (-1, 0, 1, 2)]
    print(f"      {'(slots of each type)':<26} {tot[0]:>9,} {tot[1]:>11,} {tot[2]:>9,} {tot[3]:>9,}")

    # the top slots, with each of HEAD's channels visible
    order = np.argsort(-np.abs(err))[:12]
    x = {n: wb(solve((n,))) for n in ("gdna", "rna", "lam", "theta")}
    fgl = np.asarray(cap["fg_loc"], np.float64)
    print("\n   ⭐ TOP 12 SLOTS — HEAD's own channels, and what removing each does")
    print(f"   {'#':>3} {'slot':>7} {'type':<11} {'M':>10} {'trueFg':>7} {'own':>7} {'HEAD':>7} "
          f"{'pLevel':>9} {'pRna':>9} {'pLam':>9} {'pThe':>9} "
          f"{'noLvl':>7} {'noRna':>7} {'noLam':>7}")
    print("   " + "-" * 140)
    for r, s in enumerate(order, 1):
        s = int(s)
        # ⚠ `M` and `trueFg` are the UNSPLICED bank, the population ψ's answer is a fraction OF.
        print(f"   {r:>3} {s:>7,} {TYPE_NAME.get(int(st[s]), 'BOUNDARY'):<11} {unspl[s]:>10,.0f} "
              f"{tg[s] / max(unspl[s], 1.0):>7.4f} {fgl[s]:>7.4f} {fg[s]:>7.4f} "
              f"{chans['gdna_imp (LEVEL)'][s]:>9.3g} {chans['rna_imp (certified)'][s]:>9.3g} "
              f"{chans['lam_imp (composition)'][s]:>9.3g} {chans['theta_imp (tilt)'][s]:>9.3g} "
              f"{x['gdna'][s]:>7.4f} {x['rna'][s]:>7.4f} {x['lam'][s]:>7.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
