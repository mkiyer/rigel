"""What is the certified per-object truth? Run this before debugging calibration against anything.

For every REGION and BOUNDARY slot of the chain this instrument derives the accumulator's own tally
(``count``, the mixture calibration deconvolves), the realized ``n_gdna`` / ``n_nrna`` / ``n_mrna``
from the origin-split oracle caches (each partition scanned by the same accumulator), the same RNA
split again by transcript strand (``n_rna_pos`` / ``n_rna_neg``), and ``true_f_g``, all in the drained
frame production calibrates. It refuses to write the table unless its named gates pass: sum-to-full
(the partitions reconstruct every bank of the full scan), partition-projects-exactly (``n_gdna + n_nrna
+ n_mrna == count`` at every slot), gdna-field-uniformity (at capture-OFF every slot's expected gDNA
count is ``rho_ref x E_g``, scored as Poisson z per slot and a ratio-z per object class, the signature
of a counting or divisor bug), exact-zeros (a bank that must be empty is identically zero, and an empty
bank is vacuous rather than a pass), nascent-in-annotation (nascent RNA only where a transcript span
admits RNA) and rna-strands-close (the strand split and the mature/nascent split are two partitions of
the same reads). Two certification levels: COMPOSITION (``true_f_g``, no opportunity model anywhere in
it) and FIELD (gdna-field-uniformity as well, so a density ``n/E`` may be trusted). A merely plausible
oracle is how a calibration bug and a truth bug survive each other. The analytic per-object RNA
expectation is not certified here; the RNA truth is realized. No solver runs. Writes ``slot_truth.npz``
beside each oracle cache, which every slot-scored instrument reads.

``--build`` first BUILDS the origin-split caches every truth instrument reads — the oracle BAM split by
read-name origin, each partition scanned by the production scanner, the two transcript-strand
partitions beside them, and the whole scan copied in as ``_main`` — keyed by the shipped scan-cache
loader so a stale cache is refused rather than reused, in parallel over conditions with ``--jobs``
(one condition saturates one core at ~2 GB), then certifies. ``panel.py cache`` runs it.

Usage::

    python scripts/design/calibration_oracle.py --condition <name>     # certify one condition
    python scripts/design/calibration_oracle.py                        # certify the whole ladder
    python scripts/design/calibration_oracle.py --build --jobs 8       # build every cache, then certify
    python scripts/design/calibration_oracle.py --condition <name> --out truth.npz
    python scripts/design/calibration_oracle.py --self-test            # perturb every gate, no I/O
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


from _shared import DEFAULT_INDEX, DEFAULT_SUITE, sibling  # noqa: E402



from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import (  # noqa: E402
    build_region_geometry,
    build_region_statics,
    g1_locked,
)
from rigel.calibration.signature import (  # noqa: E402
    BIT_EXON_NEG,
    BIT_EXON_POS,
    RegionType,
    coarse_type_array,
)
from rigel.calibration.splice_graph import build_boundary_flags_array, build_sj_geometry_arrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from calibration._oracle import ORIGINS, RNA_STRAND_ORIGINS, OracleTruth, lift_drain_parts  # noqa: E402

_EPS = 1.0e-12
#: |z| above which a slot is flagged under gdna-field-uniformity. 4 sigma two-sided is ~6e-5 expected
#: false flags per slot; the gate is on the flag rate and the class ratio, not on any single slot.
_Z_FLAG = 4.0
#: Poisson slots below this expectation are pooled into their class rather than z-scored singly:
#: the normal approximation is not honest there, and a per-slot z on a tiny lambda flags noise.
_MIN_LAMBDA = 5.0


# ── the gates, arrays in / verdict out, so the self-test can feed synthetic inputs ────────────────


def gate_partition(count: np.ndarray, n_g, n_n, n_m, atol: float = 1e-6) -> dict:
    """partition-projects-exactly: the origin partition must reproduce the accumulator count at every slot."""
    gap = np.abs((n_g + n_n + n_m) - count)
    worst = float(gap.max()) if gap.size else 0.0
    return {"gate": "gate partition-projects-exactly", "ok": bool(worst <= atol), "worst_gap": worst,
            "n_bad": int((gap > atol).sum())}


def gate_uniformity(n_g: np.ndarray, eff_g: np.ndarray, classes: np.ndarray,
                    select: np.ndarray) -> dict:
    """gdna-field-uniformity: Poisson uniformity of the realized gDNA field against ``rho x E_g``, per class.

    The class ratio-z is the counting-bug detector: a deposit-rule or divisor bug is per kind
    (regions vs boundaries vs a signature class), so it shows as one class sitting off the shared
    rate while per-slot flags stay unremarkable. Vacuous when the selected field is empty.
    """
    sel = select & (eff_g > 0.0)
    tot_n, tot_e = float(n_g[sel].sum()), float(eff_g[sel].sum())
    if tot_n <= 0.0:
        return {"gate": "gate gdna-field-uniformity", "ok": True, "vacuous": True}
    rho = tot_n / tot_e
    lam = rho * eff_g
    zable = sel & (lam >= _MIN_LAMBDA)
    z = (n_g[zable] - lam[zable]) / np.sqrt(lam[zable])
    flag_rate = float(np.mean(np.abs(z) > _Z_FLAG)) if z.size else 0.0
    rows = []
    ok = True
    for cls in sorted(set(classes[sel].tolist())):
        m = sel & (classes == cls)
        cn, ce = float(n_g[m].sum()), float(rho * eff_g[m].sum())
        cz = (cn - ce) / np.sqrt(max(ce, _EPS))
        rows.append({"class": cls, "n": cn, "expected": ce, "ratio": cn / max(ce, _EPS), "z": cz})
        # class totals are large, so a real bias is hundreds of sigma; 6 allows the slight
        # non-independence of the shared rho-hat without admitting any real bias.
        if abs(cz) > 6.0:
            ok = False
    return {"gate": "gate gdna-field-uniformity", "ok": ok and flag_rate < 1e-3, "vacuous": False,
            "rho": rho, "slot_flag_rate": flag_rate, "n_z_scored": int(z.size), "classes": rows}


def gate_zero(n: np.ndarray, select: np.ndarray, what: str) -> dict:
    """exact-zeros: a bank that must be empty, checked as exactly zero. Vacuous if nothing is selected."""
    if not select.any():
        return {"gate": f"gate exact-zeros:{what}", "ok": True, "vacuous": True}
    bad = float(n[select].sum())
    return {"gate": f"gate exact-zeros:{what}", "ok": bool(bad == 0.0), "vacuous": False, "mass": bad}


def gate_rna_strands_close(n_pos: np.ndarray, n_neg: np.ndarray, n_m: np.ndarray,
                           n_n: np.ndarray, atol: float = 1e-6) -> dict:
    """rna-strands-close: the transcript-strand split and the mature/nascent split are two partitions
    of the same RNA reads, so they must agree at every slot."""
    gap = np.abs((n_pos + n_neg) - (n_m + n_n))
    worst = float(gap.max()) if gap.size else 0.0
    return {"gate": "gate rna-strands-close", "ok": bool(worst <= atol), "worst_gap": worst,
            "n_bad": int((gap > atol).sum())}


def gate_nascent_scope(n_n: np.ndarray, rna_admissible: np.ndarray) -> dict:
    """nascent-in-annotation: nascent RNA may exist only where the annotation admits any RNA."""
    bad = float(n_n[~rna_admissible].sum())
    return {"gate": "gate nascent-in-annotation", "ok": bool(bad == 0.0), "out_of_scope_mass": bad}


# ── the per-slot strata and counts the truth table is keyed by, from the annotation and the payload ──

#: Either strand's exon bit: a boundary flank is "exonic" if any transcript has an exon there.
_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG

#: The seven slot populations, in report order: mutually exclusive and exhaustive over the chain, asserted
#: in :func:`stratum_labels`. The boundary axis is split by whether mature RNA can cross, not by whether a
#: sj attaches: an ``exon|exon`` boundary sits inside a contiguous exonic stretch that mature RNA crosses.
STRATA = (
    "R intergenic",
    "R intron",
    "R exon",
    "B exon|intron",
    "B exon|exon",
    "B intron|intron",
    "B gene edge",
)

#: The gDNA anchor that sits inside genes, and therefore on-target under hybrid capture.
ONTARGET_GDNA_STRATUM = "B exon|intron"


def stratum_labels(chain, statics, region_arrays) -> np.ndarray:
    """Per-slot stratum label, from the annotation alone, asserted to partition the chain.

    The boundary classification is by exon-ness of the two flanks, because a boundary's
    unspliced-crossing population is the molecules that crossed it contiguously, which mature RNA can
    do only where the template is contiguous exon on both sides:

    ================  =================================================================================
    ``B exon|intron``  exactly one flank exonic, both flanks inside a gene: mature RNA cannot cross, it
                       splices, so near-pure gDNA under sparse unspliced nascent; in-gene, so the
                       on-target gDNA anchor under capture
    ``B exon|exon``    both flanks exonic, an alternative splice site inside a contiguous exonic
                       stretch; mature RNA crosses it freely, so not an anchor
    ``B intron|intron`` neither flank exonic, both inside a gene (adjacent introns of different
                       signature); off-target, nascent-only RNA
    ``B gene edge``    at least one flank intergenic, a TSS/TES interface; the ``g1_locked`` boundary
                       class, structurally pure gDNA on both strands
    ================  =================================================================================

    ``R intergenic`` is cross-checked against ``g1_locked``, the predicate the solver pins on, and
    ``B exon|intron`` against the solver's own ``mrna_active``; if either pair separates, the labels
    describe a different population than the one the solver reasons over, and this raises.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    is_region, is_boundary = kind == REGION, kind == BOUNDARY
    sig = np.asarray(region_arrays.signature).astype(np.int64)
    n_regions = sig.shape[0]
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    slot_type = np.where(is_region, rtype[np.clip(obj, 0, max(n_regions - 1, 0))], -1)
    locked = g1_locked(np.asarray(statics.free_pos, bool), np.asarray(statics.free_neg, bool))

    # the two flanks' signatures, through the chain's own adjacency (a BOUNDARY always has a REGION on
    # both sides, so the ``-1`` branch has no cases — `build_region_statics` makes the same argument)
    slot_sig = np.where(is_region, sig[np.clip(obj, 0, max(n_regions - 1, 0))] if n_regions else 0, 0)
    left = np.clip(np.asarray(chain.left), 0, max(int(chain.n_slots) - 1, 0))
    right = np.clip(np.asarray(chain.right), 0, max(int(chain.n_slots) - 1, 0))
    sig_l = np.where(is_boundary, slot_sig[left], 0)
    sig_r = np.where(is_boundary, slot_sig[right], 0)
    exon_l = (sig_l & _EXON_BITS) != 0
    exon_r = (sig_r & _EXON_BITS) != 0
    gene_both = is_boundary & (sig_l != 0) & (sig_r != 0)

    label = np.full(int(chain.n_slots), "", dtype=object)
    label[is_region & (slot_type == int(RegionType.INTERGENIC))] = "R intergenic"
    label[is_region & (slot_type == int(RegionType.INTRON))] = "R intron"
    label[is_region & (slot_type == int(RegionType.EXON))] = "R exon"
    label[is_boundary & ~gene_both] = "B gene edge"
    label[gene_both & exon_l & exon_r] = "B exon|exon"
    label[gene_both & (exon_l ^ exon_r)] = "B exon|intron"
    label[gene_both & ~exon_l & ~exon_r] = "B intron|intron"

    counted = sum(int(np.sum(label == s)) for s in STRATA)
    if counted != int(chain.n_slots):
        raise AssertionError(
            f"the stratum labels do not partition the chain: {counted:,} labelled of "
            f"{int(chain.n_slots):,} slots. Every slot must take exactly one label."
        )
    if not np.array_equal(label == "R intergenic", is_region & locked):
        raise AssertionError(
            f"`intergenic & REGION` ({int(np.sum(label == 'R intergenic')):,}) and `g1_locked & REGION` "
            f"({int(np.sum(is_region & locked)):,}) have SEPARATED on this index. They are the same "
            "population by construction — no transcript covers an intergenic region, so neither RNA "
            "strand is admissible."
        )
    # `mrna_active_s` is the solver's own "contiguous exon on both flanks" gate, i.e. "mature RNA of
    # strand s may cross here". At an `exon|intron` boundary one flank carries no exon bit at all, so
    # it must be False on both strands.
    mature_can_cross = np.asarray(statics.mrna_active_pos, bool) | np.asarray(
        statics.mrna_active_neg, bool
    )
    bad = (label == ONTARGET_GDNA_STRATUM) & mature_can_cross
    if bad.any():
        raise AssertionError(
            f"{int(bad.sum()):,} `{ONTARGET_GDNA_STRATUM}` boundaries report `mrna_active`, i.e. the "
            "solver thinks mature RNA may cross them contiguously. The anchor's whole claim is that it "
            "cannot, so this classification and the solver's disagree."
        )
    return label


def slot_counts(payload, region_arrays, chain) -> np.ndarray:
    """One payload's unspliced/contained count per slot, the mixture ψ deconvolves and nothing else:
    ``region_contained`` at a REGION, ``boundary_unspliced`` at a BOUNDARY, exactly the populations
    :attr:`RegionGeometry.unspliced_count` carries, so a truth built from the origin partitions and an
    estimate built from the full payload are on one basis.
    """
    sub = CalibrationSubstrate.from_payload(payload, region_arrays)
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    out = np.zeros(int(chain.n_slots), np.float64)
    r, b = kind == REGION, kind == BOUNDARY
    out[r] = np.asarray(sub.region_contained.count, np.float64).sum(1)[obj[r]]
    out[b] = np.asarray(sub.boundary_unspliced.count, np.float64).sum(1)[obj[b]]
    return out


# ── the derivation ────────────────────────────────────────────────────────────────────────────────


def derive(index, region_arrays, suite: Path, condition: str) -> tuple[dict, list[dict]]:
    """One condition's certified per-slot truth table plus its gate verdicts."""
    sj = build_sj_geometry_arrays(index)
    bflags = build_boundary_flags_array(index)
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    lift: dict = {}
    kw = calibration_inputs(cache, index, lift_out=lift)
    # the drained frame: the truth certified here must describe the tally production calibrates, so
    # the whole and every partition below are drained consistently.
    payload = kw["payload"]
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    statics = build_region_statics(chain, region_arrays, bflags)
    geom = build_region_geometry(
        chain, CalibrationSubstrate.from_payload(payload, region_arrays),
        region_arrays, sj, kw["gdna_fl_pmf"], kw["rna_fl_pmf"],
    )
    label = stratum_labels(chain, statics, region_arrays)

    # sum-to-full is a hard gate inside from_parts, asserted on the drained frame, which makes it
    # the lift's own end-to-end identity check; an exception is the verdict. The parts are loaded
    # pass-one and drained by replaying the whole's choices (`from_cached_parts`).
    root = Path(suite) / "oracle_cache" / condition
    parts = {k: read_scan_cache(root / k, index).payload for k in ORIGINS}
    truth = OracleTruth.from_cached_parts(payload, parts, lift)

    count = slot_counts(payload, region_arrays, chain)
    n_g = slot_counts(truth.parts["gdna"], region_arrays, chain)
    n_n = slot_counts(truth.parts["nrna"], region_arrays, chain)
    n_m = slot_counts(truth.parts["mrna"], region_arrays, chain)
    # the same RNA reads again, keyed by transcript strand: the per-component truth. Refused rather
    # than skipped if absent: an instrument that silently drops to two arms would measure a different thing.
    try:
        strand_parts = {k: read_scan_cache(root / k, index).payload for k in RNA_STRAND_ORIGINS}
    except Exception as exc:  # noqa: BLE001
        raise FileNotFoundError(
            f"{root} has no {list(RNA_STRAND_ORIGINS)} partitions ({type(exc).__name__}). They are "
            "built by `calibration_oracle.py --build` (via `panel.py cache`) alongside the three ORIGINS ones; "
            "without them there is no per-strand RNA truth and the three-arm map cannot be scored."
        ) from exc
    # the second exact partitioning of the same whole, drained with the same choice queue, gdna
    # first in both lists so the shared member takes the identical choice slice and the two
    # partitionings' RNA remainders stay consistent (`lift_drain_parts`' docstring). `drain` is pure,
    # so the undrained `parts["gdna"]` loaded above is reusable here.
    strand_drained, n_amb_strand = lift_drain_parts(
        lift, [parts["gdna"], strand_parts["rna_pos"], strand_parts["rna_neg"]]
    )
    n_rp = slot_counts(strand_drained[1], region_arrays, chain)
    n_rn = slot_counts(strand_drained[2], region_arrays, chain)

    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    is_region = kind == REGION
    left = np.clip(np.asarray(chain.left, np.int64), 0, int(chain.n_slots) - 1)
    # a slot's reference: a REGION's own; a BOUNDARY sits between two regions of ONE reference,
    # so its left flank's region names it.
    ref = np.where(is_region, region_arrays.ref_id[np.clip(obj, 0, len(region_arrays.ref_id) - 1)], -1)
    ref = np.where(is_region, ref, ref[left])
    eff_g = np.asarray(geom.eff_gdna, np.float64)
    rna_ok = np.asarray(statics.free_pos, bool) | np.asarray(statics.free_neg, bool)

    capture_on = condition.endswith("_capture_on")
    gdna_refs = np.array(sorted({int(r) for r in ref[n_g > 0]}), np.int64)

    verdicts = [gate_partition(count, n_g, n_n, n_m), gate_rna_strands_close(n_rp, n_rn, n_m, n_n)]
    # the drained frame's own report, informational and never a pass/fail: the leak is production's
    # behaviour (`ISSUES: drain-contaminates-certified-rna`) and the ambiguity bounds the lift.
    verdicts.append({
        "gate": "drained-frame report", "ok": True, "vacuous": False,
        "gdna_spliced_leak": truth.gdna_spliced_leak,
        "n_ambiguous_origins": int(truth.n_ambiguous),
        "n_ambiguous_strands": int(n_amb_strand),
    })
    if capture_on:
        verdicts.append({"gate": "gate gdna-field-uniformity", "ok": True, "vacuous": True,
                         "note": "capture-ON: the field is deliberately non-uniform; gate not applicable"})
    else:
        for r in gdna_refs:
            v = gate_uniformity(n_g, eff_g, label, ref == r)
            v["ref"] = int(r)
            verdicts.append(v)
        if gdna_refs.size == 0:
            verdicts.append({"gate": "gate gdna-field-uniformity", "ok": True, "vacuous": True,
                             "note": "no gDNA anywhere (a zero condition)"})
    # refs that carry no gDNA at all must be exactly zero (ERCC backbone, and every ref at g00)
    verdicts.append(gate_zero(n_g, ~np.isin(ref, gdna_refs), "gdna outside genomic refs"))
    verdicts.append(gate_nascent_scope(n_n, rna_ok))

    mass = n_g + n_n + n_m
    table = {
        "condition": condition, "kind": kind, "obj": obj, "stratum": label.astype(str),
        "ref": ref, "eff_g": eff_g, "eff_r": np.asarray(geom.eff_rna, np.float64),
        "count": count, "n_gdna": n_g, "n_nrna": n_n, "n_mrna": n_m,
        "n_rna_pos": n_rp, "n_rna_neg": n_rn,
        "true_f_g": np.where(mass > 0, n_g / np.maximum(mass, _EPS), 0.0),
        "live": mass > 0,
    }
    return table, verdicts


def report(condition: str, verdicts: list[dict]) -> tuple[bool, bool]:
    """Print the verdicts and return ``(composition_ok, field_ok)``, two certification levels because
    the gates certify two different things.

    COMPOSITION-certified (sum-to-full, partition-projects-exactly, exact-zeros, rna-strands-close,
    nascent-in-annotation): ``true_f_g`` is sound, realized counts against realized counts through one
    accumulator, with no opportunity model anywhere in it. FIELD-certified (gdna-field-uniformity as
    well): the deposit geometry also matches the opportunity model, so a density ``n/E`` may be
    trusted too. The table is stamped with both verdicts rather than one.
    """
    comp_ok, field_ok = True, True
    print(f"\n== {condition}")
    for v in verdicts:
        ok = v["ok"]
        if v["gate"].startswith("gate gdna-field"):
            field_ok &= ok
        else:
            comp_ok &= ok
        mark = "✔" if ok else "⛔"
        extra = ""
        if v.get("vacuous"):
            extra = "  (VACUOUS — the selected field is empty; this is not a pass)"
        elif v["gate"].startswith("gate gdna-field") and not v.get("vacuous"):
            extra = (f"  rho={v['rho']:.6g}  slot flags {v['slot_flag_rate']:.2e} "
                     f"over {v['n_z_scored']:,} z-scored slots")
        ref_tag = f"  ref {v['ref']}" if "ref" in v else ""
        print(f"   {mark} {v['gate']:<34}{ref_tag}{extra}")
        if v["gate"].startswith("gate gdna-field") and not v.get("vacuous"):
            for c in v["classes"]:
                flag = "" if abs(c["z"]) <= 6.0 else "   ⛔ CLASS BIAS"
                print(f"        {c['class']:<18} n {c['n']:>13,.0f}  expected {c['expected']:>13,.0f}"
                      f"  ratio {c['ratio']:.4f}  z {c['z']:+8.2f}{flag}")
    return comp_ok, field_ok


# ── self-test: every gate shown to fire on the defect it exists for ───────────────────────────────


def self_test() -> int:
    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"   ⛔ {name}")

    rng = np.random.default_rng(7)
    n = 4000
    eff = rng.uniform(50, 5000, n)
    cls = np.array(["R intron", "R exon", "B exon|exon", "B exon|intron"])[np.arange(n) % 4]
    rho = 0.05
    n_g = rng.poisson(rho * eff).astype(float)
    sel = np.ones(n, bool)

    v = gate_uniformity(n_g, eff, cls, sel)
    check("uniformity passes on a genuinely uniform field", v["ok"] and not v["vacuous"])

    # ① the counting-bug shape: one class deposited at half weight
    bad = n_g.copy()
    bad[cls == "B exon|exon"] *= 0.5
    check("a half-weight class fires the class ratio-z",
          not gate_uniformity(bad, eff, cls, sel)["ok"])
    # ② a wrong divisor on one class (eff doubled) fires the same way
    bad_e = eff.copy()
    bad_e[cls == "R intron"] *= 2.0
    check("a doubled divisor on one class fires", not gate_uniformity(n_g, bad_e, cls, sel)["ok"])
    # ③ vacuous is not a pass
    check("an empty field reports VACUOUS", gate_uniformity(np.zeros(n), eff, cls, sel)["vacuous"])

    n_n = rng.poisson(2.0, n).astype(float)
    n_m = rng.poisson(10.0, n).astype(float)
    count = n_g + n_n + n_m
    check("partition gate passes when exact", gate_partition(count, n_g, n_n, n_m)["ok"])
    check("one lost fragment fires the partition gate",
          not gate_partition(count, n_g - (np.arange(n) == 17), n_n, n_m)["ok"])

    admissible = np.arange(n) % 3 != 0
    n_n2 = np.where(admissible, n_n, 0.0)
    check("rna-strand closure passes when the two partitions agree",
          gate_rna_strands_close(n_m * 0.4, n_m * 0.6 + n_n, n_m, n_n)["ok"])
    check("one fragment on the wrong strand partition fires the closure gate",
          not gate_rna_strands_close(n_m * 0.4, n_m * 0.6 + n_n - (np.arange(n) == 5), n_m, n_n)["ok"])

    check("nascent-scope passes when confined", gate_nascent_scope(n_n2, admissible)["ok"])
    leak = n_n2.copy()
    leak[0] = 1.0  # slot 0 is inadmissible
    check("one intergenic nascent fragment fires", not gate_nascent_scope(leak, admissible)["ok"])

    z = np.zeros(n)
    check("zero gate passes on exact zeros", gate_zero(z, ~admissible, "t")["ok"])
    z[3] = 1e-9
    check("1e-9 of forbidden mass fires the zero gate", not gate_zero(z, ~admissible, "t")["ok"])
    check("zero gate on empty selection is VACUOUS", gate_zero(z, np.zeros(n, bool), "t")["vacuous"])

    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def build_one(index, suite: Path, condition: str, work_dir: Path) -> None:
    """One condition's oracle cache under ``<suite>/oracle_cache/<condition>``: the three origin
    partitions and the two transcript-strand ones (`_oracle_arms.load_or_build_oracle`, which
    re-runs sum-to-full on a cache hit and rebuilds on a miss), and ``_main`` — the whole scan — copied
    from the scan cache when absent. The drained frame: the partitions are lifted by replaying the
    whole's drain, exactly as `derive` reads them."""
    OA = sibling("_oracle_arms.py")
    root = Path(suite) / "oracle_cache" / condition
    scan_dir = Path(suite) / "scan_cache" / condition
    if not (scan_dir / "payload.npz").is_file():
        raise FileNotFoundError(f"no scan cache at {scan_dir} — run build_scan_cache.py first")
    cache = read_scan_cache(scan_dir, index)
    lift: dict = {}
    kw = calibration_inputs(cache, index, lift_out=lift)
    bam = str(Path(suite) / condition / "sim_oracle.bam")
    OA.load_or_build_oracle(
        bam, index, PipelineConfig(), Path(work_dir) / f"w_{condition}", condition, kw["payload"],
        Path(suite) / "oracle_cache", lift,
    )
    if not (root / "_main" / "payload.npz").is_file():
        shutil.copytree(scan_dir, root / "_main", dirs_exist_ok=True)


def build(index, suite: Path, conds: list[str], jobs: int, work_dir: Path, args_index: Path) -> None:
    """Every condition's cache, in parallel over conditions when ``jobs > 1`` — each worker is this
    script on one condition with ``--skip-certify``; a worker that fails is rebuilt serially here, so
    the numbers do not depend on ``jobs``."""
    todo = list(conds)
    if jobs > 1 and len(todo) > 1:
        import concurrent.futures as cf
        import subprocess

        base = [sys.executable, str(Path(__file__).resolve()), "--suite", str(suite), "--index",
                str(args_index), "--work-dir", str(work_dir), "--build", "--skip-certify", "--condition"]

        def one(c):
            return c, subprocess.run(base + [c], capture_output=True, text=True).returncode

        n = max(1, min(int(jobs), len(todo)))
        print(f"  building {len(todo)} oracle cache(s), {n} worker(s) …", flush=True)
        with cf.ThreadPoolExecutor(max_workers=n) as ex:
            for c, rc in ex.map(one, todo):
                print(f"    {'✔' if rc == 0 else '⚠ rebuilding serially'} {c}", flush=True)
    for c in todo:
        build_one(index, suite, c, work_dir)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", nargs="*", default=None, help="one or more conditions (default: every one with a scan cache)")
    ap.add_argument("--out", type=Path, default=None,
                    help="write the certified per-slot table as .npz (default: <suite>/oracle_cache/<condition>/slot_truth.npz)")
    ap.add_argument("--build", action="store_true", help="build the origin-split caches before certifying")
    ap.add_argument("--jobs", type=int, default=1, help="--build's worker processes, one condition each")
    ap.add_argument("--work-dir", type=Path, default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_oracle_build",
                    help="--build's scratch for the split BAMs")
    ap.add_argument("--skip-certify", action="store_true", help="--build only (the worker half of --jobs)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    index = TranscriptIndex.load(args.index)
    region_arrays = RegionArrays.from_index(index)
    conds = args.condition or sorted(
        p.name for p in (args.suite / "scan_cache").iterdir()
    )
    if args.build:
        build(index, args.suite, conds, args.jobs, args.work_dir, args.index)
        if args.skip_certify:
            return 0
    bad = 0
    for c in conds:
        try:
            table, verdicts = derive(index, region_arrays, args.suite, c)
        except Exception as exc:  # noqa: BLE001 — **sum-to-full** failures surface here and must be a verdict
            print(f"\n== {c}\n   ⛔ gate sum-to-full / cache validation FAILED: {type(exc).__name__}: {exc}")
            bad += 1
            continue
        comp_ok, field_ok = report(c, verdicts)
        if comp_ok:
            out = args.out or (args.suite / "oracle_cache" / c / "slot_truth.npz")
            table["field_certified"] = bool(field_ok)
            np.savez_compressed(out, **{k: v for k, v in table.items() if k != "condition"})
            level = "COMPOSITION + FIELD" if field_ok else "COMPOSITION only (⛔ field gate failed — densities not certified)"
            print(f"   → table written [{level}]: {out}")
        if not comp_ok:
            bad += 1
            print("   ⛔ NOT CERTIFIED — no table written. Fix the gate, not the gate's bound.")
        elif not field_ok:
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
