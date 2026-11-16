"""The docs' citation boundaries: nothing outside `docs/dev/` may depend on the sandbox, nothing under
`src/` may cite a doc or an issue, nothing tracked cites the retired issue log, and every rule has one name.

A working note nobody cites costs nothing and can be deleted at any time; once a permanent doc, a
docstring or a test cites it, it is load-bearing. So the sandbox gate holds one property, no citation
into `docs/dev/` from outside, plus an allowlist of the files whose job is to describe the sandbox.
Each gate is one case that lists every offender, so adding a file never moves the collected count.
"""

from __future__ import annotations

import pathlib
import re
import subprocess

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
        f"has settled, MOVE it to its permanent home — an open problem to a GitHub issue, a lesson to TRAPS.md, a "
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


def test_the_permanent_set_is_eight_files():
    """The permanent docs are exactly the eight written out here, so the sandbox cannot add a ninth by
    accretion: a working note everyone reads is promoted, not blessed."""
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

#: a markdown file's name or path, a `§` with the word before it, a path into a `docs/` tree, a
#: `TRAPS: name` citation, whether or not the file or the name exists, or an `issue #N`
_CITES_A_DOC = re.compile(
    r"[\w./-]*\w\.md\b"
    r"|(?:\w+\s*)?§[\w.]*"
    r"|\bdocs/[\w./-]*"
    r"|\bTRAPS:\s*`?[a-z0-9]+(?:-[a-z0-9]+)+"
    r"|\b[Ii]ssues?\s*#\d+"
)

#: a whole kebab token: not part of a longer token or an identifier
_KEBAB = re.compile(r"(?<![\w-])[a-z0-9]+(?:-[a-z0-9]+)+(?![\w-])")

#: a TRAPS rule is a column-0 `**name.`
_TRAP_HEADING = re.compile(r"^\*\*([a-z0-9-]+)\.", re.M)


def _entry_names() -> set[str]:
    return set(_TRAP_HEADING.findall((DOCS / "TRAPS.md").read_text()))


def _doc_citations(txt: str, names: set[str]) -> list[str]:
    return sorted(set(_CITES_A_DOC.findall(txt)) | {t for t in _KEBAB.findall(txt) if t in names})


def test_the_source_cites_no_doc():
    """A doc's file name or path, a `§`, a `TRAPS:` citation, a bare rule name and an `issue #N` are all
    citations."""
    names = _entry_names()
    hits = {
        str(p.relative_to(ROOT)): found
        for p in SOURCE
        if (found := _doc_citations(p.read_text(errors="ignore"), names))
    }
    assert not hits, (
        f"the source cites the docs or an issue: {hits}. Drop the pointer and keep the content (a hit that is not a "
        f"citation is a name colliding with the code's own text)."
    )


def test_the_doc_citation_gate_is_not_vacuous():
    """Every citation form matches, ordinary text does not, and every rule heading parses as a name,
    or the gate above checks nothing."""
    traps = (DOCS / "TRAPS.md").read_text()
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
    assert _doc_citations(f"TRAPS: {name}", names) == sorted([f"TRAPS: {name}", name])
    assert _doc_citations(f"// as {name} says", names) == [name]
    assert _doc_citations("# the fix for issue #12", names) == ["issue #12"]
    assert _doc_citations("// see issues #3 and #4", names) == ["issues #3"]
    assert not _doc_citations(
        f"BY DESIGN; KNOWN ISSUES: none; #include <vector>; --{name}; {name}-x; {name}_x; "
        f"{name.replace('-', '_')}",
        names,
    )


def test_the_trap_names_are_unique():
    """Two rules under one name would be the collision naming them exists to end."""
    heads = _TRAP_HEADING.findall((DOCS / "TRAPS.md").read_text())
    dupes = sorted({h for h in heads if heads.count(h) > 1})
    assert not dupes, f"duplicate rule names: {dupes}"


#: a citation of the retired issue log: `ISSUES: name`, also across a comment's line break
_CITES_THE_ISSUE_LOG = re.compile(r"\bISSUES:[\s#*/]*`?[a-z0-9]+(?:-[a-z0-9]+)*")


def test_no_issue_log_citations():
    """The issue log moved to GitHub Issues on 2026-11-16: cite `issue #N`, not a backlog name."""
    tracked = subprocess.run(
        ["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.split("\0")
    files = [ROOT / rel for rel in tracked if rel and (ROOT / rel).is_file()]
    assert files, "git ls-files listed nothing: this gate reads the tracked files"
    found = []
    for path in files:
        if path == pathlib.Path(__file__).resolve():
            continue
        text = path.read_text(errors="ignore")
        for match in _CITES_THE_ISSUE_LOG.finditer(text):
            found.append(f"{path.relative_to(ROOT)}:{text.count(chr(10), 0, match.start()) + 1}")
    assert not found, "cites a backlog name; cite `issue #N`:\n" + "\n".join(found)
    # the pattern itself, or the gate checks nothing
    assert _CITES_THE_ISSUE_LOG.search("(`ISSUES: the-release-builds`)")
    assert _CITES_THE_ISSUE_LOG.search("ISSUES: x")
    assert _CITES_THE_ISSUE_LOG.search("# see ISSUES:\n# the-release-builds")
    assert not _CITES_THE_ISSUE_LOG.search("issue #11; docs/ROADMAP.md; MISSUES: none")
