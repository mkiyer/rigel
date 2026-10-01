"""``rigel.config.CONSTANTS`` — the fixed numbers' one home.

Every section refuses an out-of-range value by name: :data:`BROKEN` breaks every field of every section once, and
must name every field, so a constant added without a range check fails here. An experiment varies one constant on a
copy (``Constants.replaced``), never by editing the shipped instance, which is frozen.
"""

import ast
import dataclasses
import re
from pathlib import Path

import pytest

from rigel.config import CONSTANTS, Constants

#: One out-of-range value per field, per section.
BROKEN = {
    "junction_fit": {
        "kappa_floor": 0.25,
        "kappa_scan_points": 2,
        "kappa_tolerance": 0.0,
        "kkt_slack": 1.0,
        "edge_margin": 0.5,
        "barrier_floor": 0.0,
        "barrier_start": CONSTANTS.junction_fit.barrier_floor / 2.0,
        "barrier_shrink": 1.0,
        "newton_steps": 0,
        "newton_tolerance": 0.0,
        "boundary_fraction": 1.0,
        "line_search_shrink": 1.0,
        "line_search_floor": 0.0,
    },
    "landscape": {
        "grid_points": 1,
        "width_bins": 0,
        "knn_scale": 0.0,
        "reliability_sd_decades": 0.0,
        "located_var": 0.0,
    },
    "calibration": {
        "min_training_regions": 0,
        "bracket_headroom": 1.0,
        "opportunity_floor": -1.0,
        "mass_floor": -1.0,
        "fraction_tolerance": 1.0,
        "track_floor": 0.0,
        "rate_log_floor": 0.0,
    },
    "fragment_length": {
        "pool_prior_ess": -1.0,
        "gdna_mean_passes": 0,
        "gdna_mean_tolerance_bp": 0.0,
        "unseen_smoothing_ess": 0.0,
    },
    "scoring": {"strand_probability_floor": 0.5, "pruning_floor": 0.0},
    "qc": {"strand_min_observations": -1, "strand_ci_confidence": 1.0, "report_examples": 0},
    "simulator": {"oversample_ratio": 0.5, "oversample_extra": 0, "probe_length_bp": 0},
    "resources": {f.name: 0 for f in dataclasses.fields(CONSTANTS.resources)},
}


def test_every_constant_has_a_broken_value():
    assert set(BROKEN) == {f.name for f in dataclasses.fields(Constants)}
    for section, broken in BROKEN.items():
        assert set(broken) == {f.name for f in dataclasses.fields(getattr(CONSTANTS, section))}, (
            section
        )


@pytest.mark.parametrize("section", sorted(BROKEN))
def test_every_constant_refuses_an_out_of_range_value(section):
    for name, value in BROKEN[section].items():
        with pytest.raises(ValueError, match=name):
            CONSTANTS.replaced(f"{section}.{name}", value)


def test_replaced_changes_one_value_on_a_copy():
    shipped = CONSTANTS.landscape.grid_points
    changed = CONSTANTS.replaced("landscape.grid_points", shipped + 1)
    assert changed.landscape.grid_points == shipped + 1
    assert CONSTANTS.landscape.grid_points == shipped
    assert changed.landscape.knn_scale == CONSTANTS.landscape.knn_scale
    assert changed.junction_fit is CONSTANTS.junction_fit


def test_the_shipped_constants_cannot_be_reassigned():
    with pytest.raises(dataclasses.FrozenInstanceError):
        CONSTANTS.landscape.grid_points = 1
    with pytest.raises(dataclasses.FrozenInstanceError):
        CONSTANTS.landscape = CONSTANTS.landscape


# ── The gate: a fixed number lives in ``rigel.config``, never as a bare module constant ─────────────────────────

_SRC = Path(__file__).resolve().parents[1] / "src" / "rigel"

#: Modules whose numeric module constants are identifiers or formats — codes, flags, schema versions, output
#: precision — defined with the format they name, not settings.
FORMAT_MODULES = {
    "annotate.py",  # the ZF per-fragment bitfield
    "buffer.py",  # fragment class codes
    "cli.py",  # summary.json's schema version and precision
    "index.py",  # the index format version
    "report/substrate.py",  # the summary schema it reads
    "scan_payload.py",  # pool codes and axis sizes
    "splice.py",  # splice-type codes
    "types.py",
    "calibration/region_chain.py",  # node kinds
    "calibration/signature.py",  # signature bits
    "calibration/splice_graph.py",  # boundary flags and edge kinds
    "sim/bam.py",  # SAM flags
}

#: Named definitions: mathematics, or the model's own statement — not a choice anyone tunes.
DEFINITIONS = {
    ("second_pass.py", "P_ORIENTATION_GIVEN_GDNA"),  # gDNA has no strand: ½
    ("calibration/density_deconv.py", "_JEFFREYS_SHAPE"),  # the Jeffreys prior's shape, ½
    ("calibration/landscape.py", "_LN10"),  # ln 10, decades to nats
    (
        "calibration/effective_length.py",
        "UNBOUNDED_REACH",
    ),  # the finite stand-in for an unbounded reach
    (
        "sim/wgs_engine.py",
        "_BYTE_COMPLEMENT",
    ),  # a byte-indexed lookup table: one entry per byte value
}

_CONSTANT_NAME = re.compile(r"_?[A-Z][A-Z0-9_]*")


def _written_down(value: ast.expr) -> bool:
    """True when every value in the expression is a numeric literal, through any calls and operators: a number
    written down rather than derived from something named (``np`` and ``math`` count as functions, not values)."""
    nodes = list(ast.walk(value))
    numbers = [n for n in nodes if isinstance(n, ast.Constant) and type(n.value) in (int, float)]
    other_literals = [
        n for n in nodes if isinstance(n, ast.Constant) and type(n.value) not in (int, float)
    ]
    called = {n.func.id for n in nodes if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    names = {n.id for n in nodes if isinstance(n, ast.Name)}
    return bool(numbers) and not other_literals and names <= called | {"np", "math"}


def _module_constants(tree: ast.Module):
    """``(name, value)`` for each module-level single-name assignment."""
    for node in tree.body:
        if isinstance(node, ast.AnnAssign) and node.value is not None:
            target, value = node.target, node.value
        elif isinstance(node, ast.Assign) and len(node.targets) == 1:
            target, value = node.targets[0], node.value
        else:
            continue
        if isinstance(target, ast.Name):
            yield target.id, value


def test_no_bare_numeric_constant_outside_config():
    offenders = []
    for path in sorted(_SRC.rglob("*.py")):
        rel = path.relative_to(_SRC).as_posix()
        if rel == "config.py" or rel in FORMAT_MODULES:
            continue
        for name, value in _module_constants(ast.parse(path.read_text())):
            if (
                _CONSTANT_NAME.fullmatch(name)
                and _written_down(value)
                and (rel, name) not in DEFINITIONS
            ):
                offenders.append(f"{rel}: {name} = {ast.unparse(value)}")
    assert not offenders, (
        "a fixed number belongs in rigel.config.CONSTANTS, documented:\n" + "\n".join(offenders)
    )


def test_the_gate_tells_a_written_number_from_a_derived_one():
    """The gate's own perturbation: a planted constant of each written form is caught; a derived one is not."""
    written = ("1e-12", "0.25", "(0.15 * 2.0) ** 2", "1 << 20", "np.log(10.0)", "-3")
    derived = ("float(_rows.EPS)", "-np.log(np.finfo(np.float64).eps)", "2 * _Y", '("a", 1)')
    for text in written:
        assert _written_down(ast.parse(text, mode="eval").body), text
    for text in derived:
        assert not _written_down(ast.parse(text, mode="eval").body), text
