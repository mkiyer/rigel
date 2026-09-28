"""The docs' citation boundaries: nothing outside `docs/dev/` may depend on the sandbox, nothing under
`src/` may cite a doc, and a rule is cited by its name, never by a number.

Working docs are encouraged, and a note's length, staleness or wrongness is harmless. What is not
harmless is a dev doc quietly becoming the state — which has happened, with a new session pointed at a
provisional note while the permanent docs went stale beside it. The mechanism of that failure is a
citation: a note nobody cites costs nothing and can be deleted at any time, but the moment a permanent
doc, a docstring or a test points into the sandbox the note is load-bearing. That single property, plus
the shape of the allowlist that carves out the files whose job is to describe the sandbox, is all the
sandbox half holds.

`TRAPS.md`'s rules were once headed `A16`, `D4j`, `C0b`. The numbers were ambiguous as well as opaque —
one token meant a rule and also a success criterion, a moment variable or a structurally pure-gDNA
object — so the label gate scans every text file in the tree and exempts each real collision by exact
file, never in bulk.

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

#: files that are ALLOWED to name the sandbox — the ones whose job is to describe it.
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
#: :func:`test_every_ALLOWLISTED_file_actually_names_the_sandbox`, so it cannot quietly become a blanket.
SEARCHED = [
    p for p in TEXT if DEV not in p.parents and str(p.relative_to(ROOT)) not in _MAY_NAME_IT
]


def test_nothing_outside_the_sandbox_cites_into_it():
    """A citation is what turns a working note into a dependency. Everything else about a dev doc —
    length, staleness, being wrong — is harmless and is explicitly permitted."""
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
    """`TRAPS: could-the-arm-have-fired` — every citation form matches, ordinary text does not, and
    every entry heading parses as a name."""
    issues = (DOCS / "ISSUES.md").read_text()
    traps = (DOCS / "TRAPS.md").read_text()
    assert 0 < len(_ISSUE_HEADING.findall(issues)) == len(re.findall(r"^### ", issues, re.M))
    assert 0 < len(_TRAP_HEADING.findall(traps)) == len(re.findall(r"^\*\*", traps, re.M))

    names = _entry_names()
    name = min(names)
    assert _doc_citations("derived in docs/EQUATIONS.md", names) == ["docs/EQUATIONS.md"]
    assert _doc_citations("see README.md", names) == ["README.md"]  # the sandbox's own
    assert _doc_citations("see docs/ARCHITECTURE.md", names) == ["docs/ARCHITECTURE.md"]  # deleted
    assert _doc_citations("(`variance_ledger.md`)", names) == ["variance_ledger.md"]  # deleted
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


#: a rule label: a family letter and a number, optionally a sub-letter — ``A16``, ``D4j``, ``C0b``.
LABEL = re.compile(r"(?<![\w./-])([A-G]\d{1,2}[a-z]?)(?![\w-])")

#: A chain diagram, where ``E1`` is a slot id and not a rule: ``N0 E0 N1 E1 … E(k-2) N(k-1)``. Rewriting
#: one of these into prose is the same collision met from the other direction, so a line carrying other
#: slot ids is read as a diagram.
DIAGRAM = re.compile(r"\bN0\b|\bE0\b|\bB0\b|\bR0\b|\bN1\b|\bE\(k|\bR\d\b|\(N\d\|")

#: Where a bare letter+digit token is not a rule citation. Each entry is a collision that was actually
#: met, and the list is deliberately short and specific — a broad exemption would make this gate vacuous.
ALLOWED: dict[str, tuple[str, ...]] = {
    # the signature belief CLASSES on the chain. `G1` here is "a structurally pure-gDNA object", a
    # vocabulary term, and has nothing to do with the process rule that once carried the same name.
    "G1": ("*",),
    "G2": ("*",),
    "G3": ("*",),
    # SUCCESS.md numbers its own Stage-A acceptance criteria A1/A2/A3 (FIDELITY / BIAS / SUFFICIENCY) —
    # a different kind of thing from a trap, and scoped to the one file that defines them.
    # A1/A2/B1/B2 are ALSO fixture ids and a gate series in the test files below, and omitting those
    # scopes is what corrupts them: a bulk rename meets `"A1"` in a GTF gene id and in `── GATE A1:`,
    # which are data rather than citations, and rewriting the data passes the gate while destroying the
    # fixture. The scoped exemption is the remedy this file's own assertion message names.
    "A1": (
        "docs/SUCCESS.md",
        "docs/TRAPS.md",
        "tests/test_second_pass_scoring.py",
        "tests/test_sim_genomic_refs.py",
    ),
    "A2": (
        "docs/SUCCESS.md",
        "tests/test_second_pass_scoring.py",
        "tests/test_sim_genomic_refs.py",
    ),
    "A3": ("docs/SUCCESS.md",),
    "B1": ("tests/test_second_pass_scoring.py", "tests/test_sim_genomic_refs.py"),
    "B2": ("tests/test_second_pass_scoring.py",),
    # gene ids in the implicit-splice GTF fixture, two rungs of `test_vertex_reference.py`'s own G1…G6
    # case series, and the method-development test chromosome's gene ids — all beside their definitions,
    # and all naming genes or cases rather than rules. Scoped to the files that carry them, never blanket.
    # Each key is declared ONCE: a duplicate key in a dict literal is not an error, it is a lost edit,
    # and the later literal silently wins so the new scopes have no effect at all.
    "G4": (
        "tests/test_implicit_splice.py",
        "tests/calibration/test_vertex_reference.py",
        "docs/TESTING.md",
    ),
    "G5": (
        "tests/test_implicit_splice.py",
        "tests/calibration/test_vertex_reference.py",
        "docs/TESTING.md",
    ),
    # moment variables in the opportunity-tilted length quadrature
    "C2": (
        "src/rigel/calibration/effective_length.py",
        "src/rigel/native/fast_exp.h",
        "docs/TRAPS.md",
        "tests/calibration/test_certified_rna_licence.py",
    ),  # quoted inside the migration lesson
    # A test file's own case ids, scoped to the file that defines them. A local id a reader meets beside
    # its definition is legible; the same id cited from another document is not, and those citations
    # belong in prose rather than here.
    # A case-id series is scoped WHOLE. `C2b` and `G1b` are rungs of C1…C5 and G1…G6, each mirrored by
    # the `def test_<id>_…` names that carry them, so exempting only the sub-lettered rungs leaves a bulk
    # rename free to rewrite their plain-numbered siblings and break the series in half.
    "C2b": ("tests/calibration/test_certified_rna_licence.py",),
    "G1b": ("tests/calibration/test_vertex_reference.py",),
    "G14": ("src/rigel/calibration/splice_graph.py", "tests/calibration/test_splice_graph.py"),
    "G17": (
        "src/rigel/calibration/splice_graph.py",
        "tests/calibration/test_splice_graph.py",
    ),
    "G18": ("src/rigel/calibration/splice_graph.py", "tests/calibration/test_splice_graph.py"),
    # The C++ has its own letter+digit tokens and none of them is a rule, which is why the scan covers
    # `.h` and `.cpp` as well as `.py` and `.md`: a Python test can assert on a header's text, so a
    # rename that skips the header leaves the two describing different things.
    "C1": (
        "src/rigel/calibration/effective_length.py",
        "docs/TRAPS.md",
        "src/rigel/native/fast_exp.h",
        "tests/calibration/test_certified_rna_licence.py",
    ),  # Taylor coefficients 1/n!
    "C3": ("src/rigel/native/fast_exp.h", "tests/calibration/test_certified_rna_licence.py"),
    "C4": ("src/rigel/native/fast_exp.h", "tests/calibration/test_certified_rna_licence.py"),
    "C5": ("src/rigel/native/fast_exp.h", "tests/calibration/test_certified_rna_licence.py"),
    "C6": ("src/rigel/native/fast_exp.h",),
    "C7": ("src/rigel/native/fast_exp.h",),
    "C8": ("src/rigel/native/fast_exp.h",),
    "C9": ("src/rigel/native/fast_exp.h",),
    "C10": ("src/rigel/native/fast_exp.h",),
    "C11": ("src/rigel/native/fast_exp.h",),
    "F32": ("src/rigel/native/em_solver.cpp",),
    "F64": ("src/rigel/native/em_solver.cpp",),
    # TRAPS.md's own header explains what was renamed and must name the old labels to do it — one
    # paragraph, in the canonical home.
    "A16": ("docs/TRAPS.md",),
    "D4j": ("docs/TRAPS.md",),
    "C0b": ("docs/TRAPS.md",),
    # a test file's own case ids, beside their definitions
    "G6": (
        "tests/test_implicit_splice.py",
        "tests/test_scanner_accumulator_integration.py",
        "tests/calibration/test_vertex_reference.py",
        "docs/TESTING.md",
    ),
    # more of the test chromosome's gene ids — genes, not rules, exactly as `G4`/`G5`/`G6` above. Scoped
    # to the doc that tables them; the GTF itself is data and is not scanned. Declared once each — see
    # the duplicate-key note on `G4`/`G5` above before touching these.
    "G9": ("docs/TESTING.md",),
    "G10": ("docs/TESTING.md",),
    "G11": ("docs/TESTING.md",),
    "G12": ("docs/TESTING.md",),
    "G13": ("docs/TESTING.md",),
    # this file, which must name the banned labels in order to ban them
    "*": ("tests/test_docs_boundary.py",),
}


def _permitted(label: str, rel: str) -> bool:
    if rel in ALLOWED.get("*", ()):
        return True
    scopes = ALLOWED.get(label)
    return bool(scopes) and ("*" in scopes or rel in scopes)


def _numbered_labels(txt: str, rel: str) -> list[str]:
    bad = set()
    for m in LABEL.finditer(txt):
        if _permitted(m.group(1), rel):
            continue
        line = txt[txt.rfind("\n", 0, m.start()) + 1 : txt.find("\n", m.end())]
        if not DIAGRAM.search(line):  # a slot id in a chain diagram is not a rule citation
            bad.add(m.group(1))
    return sorted(bad)


def test_no_numbered_rule_labels():
    """Cite a rule by its name. `TRAPS.md` carries the full list, and every one reads as English."""
    hits = {}
    for p in TEXT:
        rel = str(p.relative_to(ROOT))
        if bad := _numbered_labels(p.read_text(errors="ignore"), rel):
            hits[rel] = bad
    assert not hits, (
        f"numbered rule labels: {hits}. Rules have NAMES — see TRAPS.md, e.g. "
        f"`TRAPS: off-grid-message-mode` rather than A16. If one of these is not a rule citation at all "
        f"(a signature class, a success criterion, a variable), add the exact scope to ALLOWED here with "
        f"the reason — never a blanket exemption."
    )


def test_PERTURBATION_a_reintroduced_label_is_caught():
    """PERTURBATION: the gate is worth nothing until it is shown to fire, so a fresh docstring citing
    `A16` must be caught, and a chain diagram's slot id must not be."""
    assert _numbered_labels('"""A docstring citing A16 the old way."""\n', "regression.py") == [
        "A16"
    ]
    assert not _numbered_labels("N0 E0 N1 E1 N2\n", "regression.py")


def test_the_trap_names_are_unique():
    """Two rules under one name would be the collision naming them exists to end."""
    heads = _TRAP_HEADING.findall((DOCS / "TRAPS.md").read_text())
    dupes = sorted({h for h in heads if heads.count(h) > 1})
    assert not dupes, f"duplicate rule names: {dupes}"


def test_the_allowlist_is_scoped_not_blanket():
    """`TRAPS: could-the-arm-have-fired` applied to this file: an allowlist that exempted everything
    would pass every test above and forbid nothing. Only the three signature belief classes may be
    repo-wide, and each of those is a vocabulary term rather than a rule."""
    wide = {k for k, v in ALLOWED.items() if "*" in v and k != "*"}
    assert wide == {"G1", "G2", "G3"}, (
        f"repo-wide exemptions are {sorted(wide)}; only the signature belief classes qualify, because "
        f"they are a vocabulary term that predates the labels and appears in hundreds of places."
    )
