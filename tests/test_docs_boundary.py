"""⭐⭐ THE SANDBOX BOUNDARY — `docs/dev/` may exist freely; nothing may DEPEND on it.

    Owner ruling, 2026-08-07: *"our docs system and rules is too strict. we need to be allowed to generate
    working docs. this is how i communicate with collaborators. we need a sandbox for docs as we are
    developing new features. they are not considered reliable."*

⭐ So working docs are **encouraged**, and the previous rule ("if you open a working doc, that is a bug")
is withdrawn. What is enforced instead is the one property whose loss actually caused harm.

⛔⛔ **A DEV DOC MUST NEVER BECOME THE STATE.** Two working docs once reached **1,181 boundaries between them —
larger than DESIGN + ROADMAP + SUCCESS combined** — and a new session was pointed at one as "THE STATE".
The failure was never that they existed. It was that nothing was ever MOVED out of them, so the provisional
copy became authoritative while the permanent docs went stale beside it.

⭐ **The mechanism of that failure is a CITATION.** A note nobody cites can be as wrong and as long as it
likes; it costs nothing and can be deleted at any time. The moment a permanent doc, a docstring or a test
points *into* the sandbox, the note has become load-bearing and cannot be deleted without breaking
something. That is the only thing this file forbids.
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

#: ⭐ The allowlist is EXCLUDED from the parametrisation rather than skipped inside it. A generated case
#: that skips is a case that DID NOT RUN, and two of them sat in the suite tally reading as "2 skipped" —
#: indistinguishable from a test that could not run for an environmental reason. The exemption is still
#: explicit, and it is now gated by :func:`test_every_ALLOWLISTED_file_actually_names_the_sandbox`, so it
#: cannot quietly become a blanket.
SEARCHED = [
    p
    for p in sorted(
        p
        for base in ("docs", "src", "tests", "scripts")
        for p in (ROOT / base).rglob("*")
        if p.suffix in (".py", ".md", ".h", ".cpp") and p.is_file() and DEV not in p.parents
    )
    + [ROOT / "CLAUDE.md"]
    # ⛔ The exclusion is applied to the WHOLE list, CLAUDE.md included. Filtering only the globbed part
    # left the one file that is appended by hand still in — and it is allowlisted, so it failed.
    if str(p.relative_to(ROOT)) not in _MAY_NAME_IT
]


@pytest.mark.parametrize("path", SEARCHED, ids=lambda p: str(p.relative_to(ROOT)))
def test_nothing_outside_the_sandbox_cites_into_it(path: pathlib.Path):
    """⛔ A citation is what turns a working note into a dependency. Everything else about a dev doc —
    length, staleness, being wrong — is harmless and is explicitly permitted."""
    rel = str(path.relative_to(ROOT))
    hits = sorted(set(_CITES_DEV.findall(path.read_text(errors="ignore"))))
    assert not hits, (
        f"{rel} cites the sandbox ({hits}). Nothing outside `docs/dev/` may depend on it. If the finding "
        f"has settled, MOVE it to its permanent home — a number to ROADMAP.md, a lesson to TRAPS.md, a "
        f"ruling to DESIGN.md, a derivation to EQUATIONS.md — and delete it from the dev doc in the same "
        f"edit. Copying is what creates two homes; moving does not."
    )


def test_every_ALLOWLISTED_file_actually_names_the_sandbox():
    """⛔ An exemption nobody re-checks becomes a blanket. Each allowlisted file is exempt because its JOB
    is to describe the sandbox — so each one must actually do that. A file that does not need the
    exemption should be back under the gate above.

    ⭐ This replaces two ``pytest.skip`` calls, and it is strictly stronger than they were: a skip cannot
    fail, so an allowlist entry for a deleted or repurposed file would have sat there indefinitely reading
    as "2 skipped".
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
    """⚠ An empty unexplained directory is an invitation to guess at the rules."""
    assert DEV.is_dir(), "docs/dev/ is the sanctioned sandbox and should exist"
    readme = DEV / "README.md"
    assert readme.is_file(), "docs/dev/README.md must state the two rules"
    txt = readme.read_text()
    assert "MOVE it out" in txt and "may cite" in txt


def test_the_permanent_set_is_still_eight_files():
    """⭐ The sandbox is additive. It must not become a ninth permanent doc by accretion — if something in
    `docs/dev/` has become something everyone reads, that is the signal to promote it, not to bless it."""
    top = sorted(p.name for p in (ROOT / "docs").glob("*.md"))
    assert top == [
        "DESIGN.md",
        "EQUATIONS.md",
        "MANUAL.md",
        "PUBLISHING.md",
        "ROADMAP.md",
        "SUCCESS.md",
        "TESTING.md",
        "TRAPS.md",
    ], f"the permanent set changed: {top}. Adding a ninth is an owner decision, not a side effect."


def test_this_gate_is_not_vacuous():
    """`TRAPS: could-the-arm-have-fired` — the pattern must actually match a citation."""
    assert _CITES_DEV.findall("see `docs/dev/message_notes.md` for the sketch")
    assert _CITES_DEV.findall("dev/foo.md")
    assert not _CITES_DEV.findall("docs/DESIGN.md and src/rigel/dev_utils.py")
