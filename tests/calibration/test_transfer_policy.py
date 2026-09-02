"""Gates for THE TRANSFER POLICY (`calibration.messages.transfer`) — the first rung of the
ground-up message rebuild (owner ruling, 2026-09-01): the intron -> intron|exon boundary
COMPOSITION TRANSFER, delivered as the intron's own factory row, strictly one hop, exon-side
message off.

The rung's anchors, in the fail-first order they were built:

* an evidence-free transfer is byte-identical to `SilentPolicy` through the real backbone —
  the rung-0 identity every new rung must keep;
* a live transfer delivers `lam_rows` at exactly the intron|exon PAIR boundaries and nowhere
  else, each row its intron source's row VERBATIM (the sender publishes its claim unchanged —
  the hop cost was measured at zero beyond the row's own ``alpha_eff`` width, so no blur
  ships; the module docstring carries the measurement), and the policy relays NOTHING (one hop
  is structural);
* `message_policy = "transfer"` actually installs the policy, and the unknown-name refusal
  stays intact.

⚠ Falsification note (2026-09-01, watched): dropping the exon-flank requirement from the pair
predicate is NOT catchable on this toy — every intron-flanking boundary here also has an exon
flank — so the perturbation that must fire instead is the SOURCE-SIDE flip (exon as source),
which the pair gate catches.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

import rigel.calibration.sweep as SW
from rigel.calibration.messages import Policy
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import REGION


def _mp():
    spec = importlib.util.spec_from_file_location(
        "tmp_for_transfer_policy", Path(__file__).parent / "test_message_policy.py"
    )
    m = importlib.util.module_from_spec(spec)
    sys.modules["tmp_for_transfer_policy"] = m
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def sweep_inputs(tmp_path_factory):
    """The message-policy gate file's captured `solve_chain` inputs, reused verbatim so every
    policy gate in this package runs on byte-identical inputs."""
    return _mp().sweep_inputs.__wrapped__(tmp_path_factory)


def _run(si, policy, capture=None):
    kw = dict(si["kw"])
    if capture is not None:
        kw["_capture"] = capture
    out = SW.solve_chain(*si["args"], **kw, policy=policy)
    return {
        f: np.asarray(getattr(out, f))
        for f in ("f_g", "f_pos", "f_neg", "var_gdna", "var_pos", "var_neg")
    }


def _expected_pairs(si):
    """The intron|exon pairs derived INDEPENDENTLY of the policy (its falsification power):
    a BOUNDARY whose one flank is an exon REGION and whose other flank is an intron REGION
    that admits RNA."""
    from rigel.calibration.signature import coarse_type_array

    chain, statics, _geometry, _belief, region_arrays = si["args"]
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    is_exon = is_reg & (rtype[np.clip(obj, 0, rtype.shape[0] - 1)] == 2)
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    is_intron = is_reg & ~is_exon & (fp | fn)
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    pairs = []
    for i in np.flatnonzero(~is_reg):
        lo, hi = left[i], right[i]
        if lo < 0 or hi < 0:
            continue
        if is_exon[lo] and is_intron[hi]:
            pairs.append((int(i), int(hi)))
        elif is_exon[hi] and is_intron[lo]:
            pairs.append((int(i), int(lo)))
    return pairs, int(chain.n_slots)


def _live_provider(si, n_grid, window):
    """A synthetic factory-row provider: a distinct non-flat row at every intron REGION slot."""
    from rigel.calibration.simplex_logodds import _logodds_grid

    pairs, n_slots = _expected_pairs(si)
    lam, _ = _logodds_grid(n_grid, window)
    rows = np.zeros((n_slots, lam.shape[0]))
    for _b, j in pairs:
        rows[j] = -0.05 * (lam - (0.1 * (j % 7) - 0.3)) ** 2  # non-flat, slot-distinct
    return lambda g, w: rows if (int(g), float(w)) == (int(n_grid), float(window)) else None


def test_the_transfer_policy_satisfies_the_backbone_protocol():
    from rigel.calibration.messages.transfer import TransferPolicy

    assert isinstance(TransferPolicy(lambda g, w: None), Policy)


def test_an_evidence_free_transfer_is_byte_identical_to_silence(sweep_inputs):
    """THE RUNG-0 IDENTITY: with no factory rows to transfer, the policy must reproduce
    `SilentPolicy` byte-for-byte through the real backbone — and it must do so by delivering a
    SILENT message, never zero-filled channel arrays (the measured 1-ULP non-identity of a
    present-but-zero channel)."""
    from rigel.calibration.messages.transfer import TransferPolicy

    a = _run(sweep_inputs, SilentPolicy())
    b = _run(sweep_inputs, TransferPolicy(lambda g, w: None))
    for f in a:
        np.testing.assert_array_equal(a[f], b[f], err_msg=f)


def test_the_transfer_delivers_only_at_pair_boundaries_and_relays_nothing(sweep_inputs):
    """One hop, one destination class: the delivered `lam_rows` are non-flat at exactly the
    independently-derived intron|exon pair boundaries, zero everywhere else, each row equal to
    its OWN intron source's row VERBATIM (the sender's claim arrives unchanged) — and the
    policy's scan is None, so nothing can travel a second hop by construction."""
    from rigel.calibration.messages.transfer import TransferPolicy

    n_grid = int(sweep_inputs["kw"]["n_grid"])
    window = float(sweep_inputs["kw"]["logodds_window"])
    provider = _live_provider(sweep_inputs, n_grid, window)
    pairs, n_slots = _expected_pairs(sweep_inputs)
    assert pairs, "the toy must carry at least one intron|exon pair or this gate proves nothing"

    pol = TransferPolicy(provider)
    prepared = pol.prepare(_ctx_of(sweep_inputs))
    assert prepared.scan(backward=False) is None and prepared.scan(backward=True) is None
    import rigel.calibration.messages as M

    msg = prepared.deliver(_dummy_nb(n_slots), _dummy_nb(n_slots))
    assert isinstance(msg, M.PsiMessage) and msg.lam_rows is not None
    rows = np.asarray(msg.lam_rows)
    src = provider(n_grid, window)
    expected = {b for b, _j in pairs}
    for i in range(n_slots):
        if i in expected:
            j = dict(pairs)[i]
            np.testing.assert_array_equal(rows[i], src[j] - src[j].max(), err_msg=f"slot {i}")
        else:
            assert not rows[i].any(), f"slot {i} received a transfer it is not licensed for"


def _ctx_of(si):
    """Rebuild the StepContext exactly as the backbone would, by running a silent sweep with a
    capture and reading nothing — instead we call solve_chain's own construction path via a spy
    on the policy prepare."""
    grabbed = []

    class _Spy:
        name = "ctx-spy"

        def prepare(self, ctx):
            grabbed.append(ctx)
            return SilentPolicy().prepare(ctx)

    _run(si, _Spy())
    assert grabbed, "the spy never fired"
    return grabbed[0]


def _dummy_nb(n_slots):
    from rigel.calibration.messages import NeighbourState

    idx = np.zeros(n_slots, np.int64)
    return NeighbourState(state=(), valid=np.zeros(n_slots, bool), src=idx)


def test_the_policy_name_installs_the_transfer_policy(sweep_inputs):
    """`message_policy = "transfer"` must install `TransferPolicy` (an unreadable knob is worse
    than no knob), and the unknown-name refusal must survive the new branch."""
    import dataclasses as _dc

    from rigel.calibration.messages.transfer import TransferPolicy
    from rigel.config import CalibrationConfig

    calibrate_mod = sys.modules["rigel.calibration.calibrate"]
    seen: list = []
    orig = SW.solve_chain

    def spy(*a, **kw):
        seen.append(kw.get("policy"))
        return orig(*a, **kw)

    calibrate_mod.solve_chain = spy
    try:
        calibrate_mod.calibrate(
            payload=sweep_inputs["payload"],
            config=_dc.replace(
                CalibrationConfig(),
                message_propagation=True,
                message_policy="transfer",
                rna_anchor=False,
            ),
            **sweep_inputs["calibrate_kw"],
        )
    finally:
        calibrate_mod.solve_chain = orig
    assert seen and all(isinstance(p, TransferPolicy) for p in seen)
    with pytest.raises(ValueError, match="unknown message_policy"):
        calibrate_mod.calibrate(
            payload=sweep_inputs["payload"],
            config=_dc.replace(
                CalibrationConfig(), message_propagation=True, message_policy="no-such-policy"
            ),
            **sweep_inputs["calibrate_kw"],
        )
