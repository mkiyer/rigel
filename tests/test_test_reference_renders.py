"""THE TEST CHROMOSOME HAS ONE SOURCE — `scripts/sim/test_reference/test_chr.yaml` — AND THE
CHECKED-IN RENDERS MUST MATCH IT (owner ruling, 2026-09-02).

`test_chr.gtf`, `test_shadow.gtf`, `test_abundances.tsv` and the three probe panels are DERIVED by
`scripts/sim/build_test_reference.py` from the YAML and versioned beside it so a reader sees the
annotation without running anything. A hand edit to a render, or a YAML edit without re-running the
builder, leaves the benchmark describing two different chromosomes — this gate refuses both.

Falsified by writing: a hand-edited render is caught (the builder's own self-test perturbs it), and the
strand-balance rule below is watched firing on an all-plus chromosome.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SIM = REPO / "scripts" / "sim"
SPEC = SIM / "test_reference" / "test_chr.yaml"


@pytest.fixture(scope="module")
def builder():
    spec = importlib.util.spec_from_file_location(
        "build_test_reference", SIM / "build_test_reference.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_yaml_is_well_formed(builder):
    spec = builder.load_spec(SPEC)
    assert builder.check_spec(spec) == []


def test_every_checked_in_render_matches_the_yaml(builder):
    spec = builder.load_spec(SPEC)
    drift = builder.check_renders(spec, SPEC.parent)
    assert drift == [], (
        f"{drift} differ from what {SPEC.name} renders — run "
        "`python scripts/sim/build_test_reference.py` (never edit a render by hand)"
    )


def test_both_strands_are_equally_represented(builder):
    """⭐ Owner, 2026-09-02: every gene carries an explicit strand and the chromosome keeps equal
    representation, so a sign error in a strand-dependent rule cannot hide on a one-strand substrate."""
    spec = builder.load_spec(SPEC)
    pos = sum(1 for g in spec.genes if g.strand == "+")
    neg = len(spec.genes) - pos
    assert pos > 0 and neg > 0 and abs(pos - neg) <= 1, (pos, neg)


def test_every_type_sits_on_both_strands(builder):
    """A block type that only ever appears on one strand is the one-strand substrate in disguise.
    A both-stranded locus is TWO genes on opposite strands sharing one type (2026-09-07): their ids carry
    a role suffix, ``_H`` the host and ``_A`` the antisense, which is not part of the type."""
    spec = builder.load_spec(SPEC)
    by_type: dict[str, set[str]] = {}
    for g in spec.genes:
        parts = g.gene_id.split("_", 1)
        if len(parts) == 2 and parts[0].startswith("gB"):
            token = parts[1]
            if token.endswith(("_H", "_A")):  # a both-stranded locus's role suffix is not its type
                token = token[:-2]
            by_type.setdefault(token, set()).add(g.strand)
    assert by_type, "no block genes found"
    one_sided = sorted(t for t, s in by_type.items() if s != {"+", "-"})
    assert one_sided == [], one_sided


def test_the_balance_rule_fires_on_a_one_strand_chromosome(builder):
    g = builder.Gene("g", "+", False, [])
    assert builder.check_strand_balance([g, g, g]) != []
    assert builder.check_strand_balance([g, builder.Gene("h", "-", False, [])]) == []
