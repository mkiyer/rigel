"""⛔⛔⛔ NUMBERED RULE LABELS ARE BANNED, AND THIS IS WHAT MAKES THAT TRUE.

    Owner ruling, 2026-08-07: *"the roadmap and consolidate documents are complete jargon. 'd4j', a16. it
    is unintelligible to humans. Outlaw jargon references, so step2j, d4j, a15 disappear and we have
    actual issue names."*

`TRAPS.md`'s rules used to be `A16`, `D4j`, `C0b`. They are now **named**: `off-grid-message-mode`,
`a-cancelling-defect-pair`, `frame-free-is-not-assumption-free`. A citation says what it means without a
lookup, and it is still one greppable string with one home.

⚠ **The numbers were not merely opaque — they were AMBIGUOUS, and that is the part a style preference
would not have caught.** Measured over the whole tree before the rename:

============  ==========================================  ===============================================
label         meant, in `TRAPS.md`                        …and ALSO meant
============  ==========================================  ===============================================
``A1``        a validator that checks itself              the **FIDELITY** criterion in `SUCCESS.md`
``C1``        a purity filter is a length filter          a **moment variable** in `effective_length.py`
``G1``        **no magic numbers**                        a **structurally pure-gDNA object** — 201 uses
============  ==========================================  ===============================================

⭐ A name cannot collide that way: nothing else in this repo is called `no-magic-numbers`.

⛔ **Why a test and not a convention.** The rename rewrote **980 citations across 114 files**. A convention
would let the next one back in one docstring at a time, and the reason the labels became unreadable in the
first place is that nothing ever said no.
"""

from __future__ import annotations

import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]

#: a rule label: a family letter and a number, optionally a sub-letter — ``A16``, ``D4j``, ``C0b``.
LABEL = re.compile(r"(?<![\w./-])([A-G]\d{1,2}[a-z]?)(?![\w-])")

#: ⛔ A CHAIN DIAGRAM, where ``E1`` is a SLOT ID and not a rule. The first pass of the rename rewrote 20 of
#: these into citations — ``N0 E0 N1 E1 … E(k-2) N(k-1)`` became prose — which is the same collision the
#: rename exists to end, met in the other direction. A boundary carrying other slot ids is a diagram.
DIAGRAM = re.compile(r"\bN0\b|\bE0\b|\bB0\b|\bR0\b|\bN1\b|\bE\(k|\bR\d\b|\(N\d\|")

#: ⭐ Where a bare letter+digit token is NOT a rule citation. Each entry is a MEASURED collision, and the
#: list is deliberately short and specific — a broad exemption would make this gate vacuous.
ALLOWED: dict[str, tuple[str, ...]] = {
    # the signature belief CLASSES on the chain. `G1` here is "a structurally pure-gDNA object", the
    # DESIGN.md vocabulary term, and has nothing to do with the process rule of the same old name.
    "G1": ("*",),
    "G2": ("*",),
    "G3": ("*",),
    # SUCCESS.md numbers its own Stage-A acceptance criteria A1/A2/A3 (FIDELITY / BIAS / SUFFICIENCY).
    # ⚠ They are a different KIND of thing from a trap, and they are scoped to one file.
    # ⛔⛔ A1/A2/B1/B2 ARE ALSO FIXTURE IDS AND A GATE SERIES, AND OMITTING THOSE SCOPES IS WHAT CORRUPTED
    # THEM. The rename met `"A1"` in a GTF gene id and in `── GATE A1:` — data, not citations — and
    # resolved the collision by REWRITING THE DATA to `"TRAPS: self-checking-validator"`. 127 gene ids
    # across 15 test files and two script labels went that way, and the suite stayed green because a
    # gene name is arbitrary and the substitution was injective. The scoped exemption is the remedy this
    # file's own assertion message names; destroying the fixture is not.
    # ⚠ `verify_toy_substrate.py`'s series is S1/S2/S3/A1/A2 — the S rungs survived only because `S`
    # falls outside LABEL's `[A-G]` class, which is what proves A1/A2 there are the same KIND of label.
    "A1": (
        "docs/SUCCESS.md",
        "docs/TRAPS.md",
        "tests/test_second_pass_scoring.py",
        "tests/test_sim_genomic_refs.py",
        "scripts/design/verify_toy_substrate.py",
    ),
    "A2": (
        "docs/SUCCESS.md",
        "tests/test_second_pass_scoring.py",
        "tests/test_sim_genomic_refs.py",
        "scripts/design/verify_toy_substrate.py",
    ),
    "A3": ("docs/SUCCESS.md",),
    "B1": ("tests/test_second_pass_scoring.py", "tests/test_sim_genomic_refs.py"),
    "B2": ("tests/test_second_pass_scoring.py",),
    # gene ids in the implicit-splice GTF fixture, and two rungs of `test_vertex_reference.py`'s own
    # G1…G6 case series — both beside their definitions
    # ⭐ AND the METHOD-DEVELOPMENT TEST CHROMOSOME's gene ids (owner's convention, 2026-08-19:
    #   *"Genes with capital G"*), which name genes rather than rules exactly as `G1`/`G2`/`G3` above
    #   name signature classes. ⛔ Scoped to the files that carry them; never a blanket exemption.
    #   ⚠ These keys are declared ONCE — `G4`/`G5`/`G6` appeared twice in this dict when the test
    #   chromosome was added and the LATER literal silently won, so the new scopes had no effect and
    #   the gate stayed red. A duplicate key in a dict literal is not an error; it is a lost edit.
    "G4": ("tests/test_implicit_splice.py", "tests/calibration/test_vertex_reference.py",
           "docs/TESTING.md"),
    "G5": ("tests/test_implicit_splice.py", "tests/calibration/test_vertex_reference.py",
           "docs/TESTING.md"),
    # moment variables in the opportunity-tilted length quadrature
    "C2": (
        "src/rigel/calibration/effective_length.py",
        "src/rigel/native/fast_exp.h",
        "docs/TRAPS.md",
        "tests/calibration/test_certified_rna_licence.py",
    ),  # quoted inside the migration lesson
    # ⚠ A TEST FILE'S OWN CASE IDS, scoped to the file that defines them. A local id a reader meets beside
    # its definition is legible; the SAME id cited from another document is not, and those citations were
    # rewritten into prose rather than exempted here.
    # ⛔⛔ AND SCOPING ONLY THE SUB-LETTERED SURVIVORS IS WHAT DESTROYED THEIR PLAIN-NUMBERED SIBLINGS.
    # `C2b` and `G1b` are two rungs of two series — C1…C5 and G1…G6 — each mirrored in this repo by the
    # `def test_<id>_…` names that still carry them. The bulk rename rewrote the plain-numbered rungs of
    # both into trap names and left the sub-lettered ones alone, so the surviving exemption below was
    # evidence that the whole series belonged here. Restored 2026-08-17 with the series scoped whole.
    "C2b": ("tests/calibration/test_certified_rna_licence.py",),
    "G1b": ("tests/calibration/test_vertex_reference.py",),
    "D1": ("tests/calibration/test_solvability_denominator.py",),
    "D2": ("tests/calibration/test_solvability_denominator.py",),
    "D3": ("tests/calibration/test_solvability_denominator.py",),
    "D4": ("tests/calibration/test_solvability_denominator.py",),
    "G14": ("src/rigel/calibration/splice_graph.py", "tests/calibration/test_splice_graph.py"),
    "G17": (
        "src/rigel/calibration/splice_graph.py",
        "tests/calibration/test_splice_graph.py",
        "tests/conftest.py",
    ),
    "G18": ("src/rigel/calibration/splice_graph.py", "tests/calibration/test_splice_graph.py"),
    # ⛔ THE C++ HAS ITS OWN LETTER+DIGIT TOKENS AND NONE OF THEM IS A RULE. Missing them is how the
    # migration was found INCOMPLETE: a Python test asserting a deletion record said "DELETED by C2" and
    # the header it reads still did, because the first pass covered only .py and .md.
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
    # ⭐ TRAPS.md's own header explains WHAT WAS RENAMED and must name the old labels to do it — that is
    # the migration note a reader needs, and it is one paragraph in the canonical home.
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
    # the test chromosome's STAGE 3 gene ids (2026-08-19: the gap-closing structures, `TESTING.md`
    # §0a) — genes, not rules, exactly as `G4`/`G5`/`G6` above. ⛔ Scoped to the doc that tables them;
    # the GTF itself is data and is not scanned. ⚠ Declared ONCE each — see the duplicate-key lesson
    # on `G4`/`G5` above before touching these.
    "G9": ("docs/TESTING.md",),
    "G10": ("docs/TESTING.md",),
    "G11": ("docs/TESTING.md",),
    "G12": ("docs/TESTING.md",),
    "G13": ("docs/TESTING.md",),
    # this file, which must name the banned labels in order to ban them
    "*": ("tests/test_no_jargon_labels.py",),
}

FILES = sorted(
    p
    for base in ("docs", "src", "tests", "scripts")
    for p in (ROOT / base).rglob("*")
    if p.suffix in (".py", ".md", ".h", ".cpp") and p.is_file()
)


def _permitted(label: str, rel: str) -> bool:
    if rel in ALLOWED.get("*", ()):
        return True
    scopes = ALLOWED.get(label)
    return bool(scopes) and ("*" in scopes or rel in scopes)


@pytest.mark.parametrize("path", FILES, ids=lambda p: str(p.relative_to(ROOT)))
def test_no_numbered_rule_labels(path: pathlib.Path):
    """⛔ Cite a rule by its NAME. `TRAPS.md` carries the full list, and every one reads as English."""
    rel = str(path.relative_to(ROOT))
    txt = path.read_text(errors="ignore")
    bad = set()
    for m in LABEL.finditer(txt):
        if _permitted(m.group(1), rel):
            continue
        boundary = txt[txt.rfind("\n", 0, m.start()) + 1 : txt.find("\n", m.end())]
        if DIAGRAM.search(boundary):
            continue  # a slot id in a chain diagram, not a rule citation
        bad.add(m.group(1))
    bad = sorted(bad)
    assert not bad, (
        f"{rel} cites numbered rule labels {bad}. Rules have NAMES — see TRAPS.md, e.g. "
        f"`TRAPS: off-grid-message-mode` rather than A16. If one of these is not a rule citation at all "
        f"(a signature class, a success criterion, a variable), add the exact scope to ALLOWED here with "
        f"the reason — never a blanket exemption."
    )


def test_every_traps_rule_has_a_name_not_a_number():
    """The canonical home. Each entry's heading is now its name, and the name is the identifier."""
    txt = (ROOT / "docs" / "TRAPS.md").read_text()
    heads = re.findall(r"^\*\*([a-z0-9-]+)\.", txt, re.M)
    assert len(heads) >= 90, f"only {len(heads)} rule headings found — the file shape changed"
    numbered = re.findall(r"^\*\*([A-G]\d{1,2}[a-z]?)\.", txt, re.M)
    assert not numbered, f"rules still headed by a number: {numbered}"
    for h in heads:
        assert re.fullmatch(r"[a-z0-9-]+", h.strip()), (
            f"rule name {h.strip()!r} is not lower-case-kebab — the names are identifiers as well as "
            f"prose, and one greppable form is what lets a citation have exactly one home."
        )


def test_the_names_are_unique():
    """Two rules under one name would be the collision this rename exists to end."""
    heads = [
        h.strip()
        for h in re.findall(r"^\*\*([a-z0-9-]+)\.", (ROOT / "docs" / "TRAPS.md").read_text(), re.M)
    ]
    dupes = sorted({h for h in heads if heads.count(h) > 1})
    assert not dupes, f"duplicate rule names: {dupes}"


def test_PERTURBATION_a_reintroduced_label_is_caught(tmp_path):
    """⛔ The gate is worth nothing until it is shown to fire. A new docstring citing `A16` must fail."""
    p = tmp_path / "regression.py"
    p.write_text('"""A docstring citing A16 the old way."""\n')
    bad = sorted(
        {
            m.group(1)
            for m in LABEL.finditer(p.read_text())
            if not _permitted(m.group(1), "regression.py")
        }
    )
    assert bad == ["A16"]


def test_the_allowlist_is_scoped_not_blanket():
    """⚠ `TRAPS: could-the-arm-have-fired` applied to this file: an allowlist that exempted everything would
    pass every test above and forbid nothing. Only the three signature classes may be repo-wide, and each
    of those is a DESIGN.md vocabulary term rather than a rule."""
    wide = {k for k, v in ALLOWED.items() if "*" in v and k != "*"}
    assert wide == {"G1", "G2", "G3"}, (
        f"repo-wide exemptions are {sorted(wide)}; only the signature belief classes qualify, because "
        f"they are a vocabulary term that predates the labels and appears in hundreds of places."
    )
