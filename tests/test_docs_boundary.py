"""The docs' citation boundaries: nothing outside `docs/dev/` may depend on the sandbox, nothing under
`src/` may cite a doc, and every rule has one name.

A working note nobody cites costs nothing and can be deleted at any time; once a permanent doc, a
docstring or a test cites it, it is load-bearing. So the sandbox gate holds one property, no citation
into `docs/dev/` from outside, plus an allowlist of the files whose job is to describe the sandbox.
Each gate is one case that lists every offender, so adding a file never moves the collected count.
"""

from __future__ import annotations

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEV = ROOT / "docs" / "dev"

#: anything that points into the sandbox — `docs/dev/x.md`, `dev/x.md`, or a bare sibling filename cited
#: from outside. The first two are what a real citation looks like.
_CITES_DEV = re.compile(r"docs/dev/|(?<![\w/])dev/[a-zA-Z0-9_.-]+\.md")

#: files that are allowed to name the sandbox: the ones whose job is to describe it.
_MAY_NAME_IT = frozenset(
    {
        "CLAUDE.md",
        "docs/dev/README.md",
        "tests/test_docs_boundary.py",
    }
)

#: every text file the gates read, the sandbox included
TEXT = sorted(
    p
    for base in ("docs", "src", "tests", "scripts")
    for p in (ROOT / base).rglob("*")
    if p.suffix in (".py", ".md", ".h", ".cpp") and p.is_file()
)

#: the sandbox gate's reach. The allowlist is excluded here and held to its job by
#: :func:`test_every_allowlisted_file_actually_names_the_sandbox`, so it cannot quietly become a blanket.
SEARCHED = [
    p for p in TEXT if DEV not in p.parents and str(p.relative_to(ROOT)) not in _MAY_NAME_IT
]


def test_nothing_outside_the_sandbox_cites_into_it():
    """A citation is what turns a working note into a dependency; anything else about a dev doc is
    permitted."""
    hits = {
        str(p.relative_to(ROOT)): found
        for p in SEARCHED
        if (found := sorted(set(_CITES_DEV.findall(p.read_text(errors="ignore")))))
    }
    assert not hits, (
        f"these cite the sandbox: {hits}. Nothing outside `docs/dev/` may depend on it. If the finding "
        f"has settled, MOVE it to its permanent home — an issue or refusal to ISSUES.md, a lesson to TRAPS.md, a "
        f"ruling to DESIGN.md, a derivation to EQUATIONS.md — and delete it from the dev doc in the same "
        f"edit. Copying is what creates two homes; moving does not."
    )


def test_every_allowlisted_file_actually_names_the_sandbox():
    """Each allowlisted file is exempt because its job is to describe the sandbox, so each one must
    actually name it; a file that does not belongs back under the gate above. Unlike a skip, this fails
    on an entry for a deleted or repurposed file.
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


def test_the_permanent_set_is_nine_files():
    """The permanent docs are exactly the nine written out here, so the sandbox cannot add a tenth by
    accretion: a working note everyone reads is promoted, not blessed."""
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
    """The pattern must actually match a citation, or the gate above checks nothing."""
    assert _CITES_DEV.findall("see `docs/dev/message_notes.md` for the sketch")
    assert _CITES_DEV.findall("dev/foo.md")
    assert not _CITES_DEV.findall("docs/DESIGN.md and src/rigel/dev_utils.py")


DOCS = ROOT / "docs"

#: every text source under `src/`
SOURCE = sorted(
    p
    for p in (ROOT / "src").rglob("*")
    if p.suffix in (".py", ".h", ".cpp", ".c", ".js", ".css") and p.is_file()
)

#: a markdown file's name or path, a `§` with the word before it, a path into a `docs/` tree, or an
#: `ISSUES: name` / `TRAPS: name` citation, whether or not the file or the name exists
_CITES_A_DOC = re.compile(
    r"[\w./-]*\w\.md\b"
    r"|(?:\w+\s*)?§[\w.]*"
    r"|\bdocs/[\w./-]*"
    r"|\b(?:ISSUES|TRAPS):\s*`?[a-z0-9]+(?:-[a-z0-9]+)+"
)

#: a whole kebab token: not part of a longer token or an identifier
_KEBAB = re.compile(r"(?<![\w-])[a-z0-9]+(?:-[a-z0-9]+)+(?![\w-])")

#: an ISSUES entry is a `### name` heading, a TRAPS rule a column-0 `**name.`
_ISSUE_HEADING = re.compile(r"^### ([a-z0-9-]+)\s*$", re.M)
_TRAP_HEADING = re.compile(r"^\*\*([a-z0-9-]+)\.", re.M)


def _entry_names() -> set[str]:
    return set(_ISSUE_HEADING.findall((DOCS / "ISSUES.md").read_text())) | set(
        _TRAP_HEADING.findall((DOCS / "TRAPS.md").read_text())
    )


def _doc_citations(txt: str, names: set[str]) -> list[str]:
    return sorted(set(_CITES_A_DOC.findall(txt)) | {t for t in _KEBAB.findall(txt) if t in names})


def test_the_source_cites_no_doc():
    """A doc's file name or path, a `§`, an `ISSUES:` / `TRAPS:` citation and a bare entry name are all
    citations."""
    names = _entry_names()
    hits = {
        str(p.relative_to(ROOT)): found
        for p in SOURCE
        if (found := _doc_citations(p.read_text(errors="ignore"), names))
    }
    assert not hits, (
        f"the source cites the docs: {hits}. Drop the pointer and keep the content (a hit that is not a "
        f"citation is a name colliding with the code's own text)."
    )


def test_the_doc_citation_gate_is_not_vacuous():
    """Every citation form matches, ordinary text does not, and every entry heading parses as a name,
    or the gate above checks nothing."""
    issues = (DOCS / "ISSUES.md").read_text()
    traps = (DOCS / "TRAPS.md").read_text()
    assert 0 < len(_ISSUE_HEADING.findall(issues)) == len(re.findall(r"^### ", issues, re.M))
    assert 0 < len(_TRAP_HEADING.findall(traps)) == len(re.findall(r"^\*\*", traps, re.M))

    names = _entry_names()
    name = min(names)
    assert _doc_citations("derived in docs/EQUATIONS.md", names) == ["docs/EQUATIONS.md"]
    assert _doc_citations("see README.md", names) == ["README.md"]  # the sandbox's own
    assert _doc_citations("see docs/ARCHITECTURE.md", names) == ["docs/ARCHITECTURE.md"]  # absent
    assert _doc_citations("(`variance_ledger.md`)", names) == ["variance_ledger.md"]  # absent
    assert _doc_citations("(`splice_graph`, DESIGN §2); METHODS §3", names) == [
        "DESIGN §2",
        "METHODS §3",
    ]
    assert _doc_citations("# not continuous across this edge (§12.9)", names) == ["§12.9"]
    assert _doc_citations("// see design §9.1", names) == ["design §9.1"]
    assert _doc_citations("# See docs/em_strand/03+05 + gdna_strand.py.", names) == [
        "docs/em_strand/03"
    ]
    assert _doc_citations("(`TRAPS: a-deleted-rule`)", names) == ["TRAPS: a-deleted-rule"]
    assert _doc_citations(f"ISSUES: {name}", names) == sorted([f"ISSUES: {name}", name])
    assert _doc_citations(f"// as {name} says", names) == [name]
    assert not _doc_citations(
        f"BY DESIGN; KNOWN ISSUES: none; --{name}; {name}-x; {name}_x; {name.replace('-', '_')}",
        names,
    )


def test_the_trap_names_are_unique():
    """Two rules under one name would be the collision naming them exists to end."""
    heads = _TRAP_HEADING.findall((DOCS / "TRAPS.md").read_text())
    dupes = sorted({h for h in heads if heads.count(h) > 1})
    assert not dupes, f"duplicate rule names: {dupes}"
