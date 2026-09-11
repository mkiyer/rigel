"""The sandbox boundary: `docs/dev/` may exist freely and nothing outside it may depend on it.

Working docs are encouraged, and a note's length, staleness or wrongness is harmless. What is not
harmless is a dev doc quietly becoming the state — which has happened, with a new session pointed at a
provisional note while the permanent docs went stale beside it. The mechanism of that failure is a
citation: a note nobody cites costs nothing and can be deleted at any time, but the moment a permanent
doc, a docstring or a test points into the sandbox the note is load-bearing. That single property, plus
the shape of the allowlist that carves out the files whose job is to describe the sandbox, is all this
file holds.
"""

from __future__ import annotations

import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEV = ROOT / "docs" / "dev"

#: anything that points into the sandbox — `docs/dev/x.md`, `dev/x.md`, or a bare sibling filename cited
#: from outside. The first two are what a real citation looks like.
_CITES_DEV = re.compile(r"docs/dev/|(?<![\w/])dev/[a-zA-Z0-9_.-]+\.md")

#: files that are ALLOWED to name the sandbox — the ones whose job is to describe it.
_MAY_NAME_IT = frozenset(
    {
        "CLAUDE.md",
        "docs/dev/README.md",
        "tests/test_docs_boundary.py",
    }
)

#: The allowlist is EXCLUDED from the parametrisation rather than skipped inside it. A generated case that
#: skips is a case that did not run, and reads in the tally as indistinguishable from one that could not
#: run for an environmental reason. The exemption stays explicit and is gated by
#: :func:`test_every_ALLOWLISTED_file_actually_names_the_sandbox`, so it cannot quietly become a blanket.
SEARCHED = [
    p
    for p in sorted(
        p
        for base in ("docs", "src", "tests", "scripts")
        for p in (ROOT / base).rglob("*")
        if p.suffix in (".py", ".md", ".h", ".cpp") and p.is_file() and DEV not in p.parents
    )
    + [ROOT / "CLAUDE.md"]
    # The exclusion is applied to the WHOLE list, CLAUDE.md included: filtering only the globbed part
    # leaves the one hand-appended file under the gate, and it is allowlisted, so it fails.
    if str(p.relative_to(ROOT)) not in _MAY_NAME_IT
]


@pytest.mark.parametrize("path", SEARCHED, ids=lambda p: str(p.relative_to(ROOT)))
def test_nothing_outside_the_sandbox_cites_into_it(path: pathlib.Path):
    """A citation is what turns a working note into a dependency. Everything else about a dev doc —
    length, staleness, being wrong — is harmless and is explicitly permitted."""
    rel = str(path.relative_to(ROOT))
    hits = sorted(set(_CITES_DEV.findall(path.read_text(errors="ignore"))))
    assert not hits, (
        f"{rel} cites the sandbox ({hits}). Nothing outside `docs/dev/` may depend on it. If the finding "
        f"has settled, MOVE it to its permanent home — an issue or refusal to ISSUES.md, a lesson to TRAPS.md, a "
        f"ruling to DESIGN.md, a derivation to EQUATIONS.md — and delete it from the dev doc in the same "
        f"edit. Copying is what creates two homes; moving does not."
    )


def test_every_ALLOWLISTED_file_actually_names_the_sandbox():
    """An exemption nobody re-checks becomes a blanket. Each allowlisted file is exempt because its job
    is to describe the sandbox, so each one must actually do that; a file that does not need the
    exemption belongs back under the gate above. Stronger than skipping the allowlisted cases would be,
    because a skip cannot fail and an entry for a deleted or repurposed file would sit there forever.
    """
    for rel in sorted(_MAY_NAME_IT):
        path = ROOT / rel
        assert path.is_file(), (
            f"{rel} is allowlisted to cite the sandbox but does not exist — the exemption is stale"
        )
        hits = sorted(set(_CITES_DEV.findall(path.read_text(errors="ignore"))))
        assert hits, (
            f"{rel} is exempt from the sandbox-citation rule because its job is to DESCRIBE the sandbox, "
            f"and it does not name it at all. Remove it from _MAY_NAME_IT and let the gate cover it."
        )


def test_the_sandbox_exists_and_says_what_it_is_for():
    """An empty unexplained directory is an invitation to guess at the rules."""
    assert DEV.is_dir(), "docs/dev/ is the sanctioned sandbox and should exist"
    readme = DEV / "README.md"
    assert readme.is_file(), "docs/dev/README.md must state the two rules"
    txt = readme.read_text()
    assert "MOVE it out" in txt and "may cite" in txt


def test_the_permanent_set_is_still_nine_files():
    """The sandbox is additive, and must not become a tenth permanent doc by accretion — if a working
    note has become something everyone reads, that is the signal to promote it, not to bless it. Adding
    to the permanent set is an owner decision, which is why the list is written out here."""
    top = sorted(p.name for p in (ROOT / "docs").glob("*.md"))
    assert top == [
        "DESIGN.md",
        "EQUATIONS.md",
        "ISSUES.md",
        "MANUAL.md",
        "PUBLISHING.md",
        "ROADMAP.md",
        "SUCCESS.md",
        "TESTING.md",
        "TRAPS.md",
    ], f"the permanent set changed: {top}. Adding a tenth is an owner decision, not a side effect."


def test_this_gate_is_not_vacuous():
    """`TRAPS: could-the-arm-have-fired` — the pattern must actually match a citation."""
    assert _CITES_DEV.findall("see `docs/dev/message_notes.md` for the sketch")
    assert _CITES_DEV.findall("dev/foo.md")
    assert not _CITES_DEV.findall("docs/DESIGN.md and src/rigel/dev_utils.py")
