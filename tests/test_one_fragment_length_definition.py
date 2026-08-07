"""⭐ TRAPS: pure-and-length-censored — there is ONE definition of fragment length in this tree, and one place that measures it.

 found **three** live definitions of "fragment length", two of them
summed into a single array called ``global_model`` and used as the empirical-Bayes anchor for pools
measured by the third. This module is the standing gate that they do not come back.

⛔ **A partial delete that still compiles is the failure mode**, which is why the gate is a search over
the source rather than a behavioural assertion: the scanner could keep recording a second histogram
forever without any test noticing, exactly as it did.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"

#: Names that existed only to serve the scanner's own fragment-length histogram. Each one is a
#: distinct piece of the machinery, listed separately so a failure names what survived.
DELETED_BY_C2 = {
    "FragmentLengthModels": (
        "the scanner's plural container (global + per-SpliceType raw histograms). "
        "⚠ FragmentLengthModel (SINGULAR) stays — it is the scoring / effective-length model built "
        "by from_pmf, and deleting it would delete the scorer."
    ),
    "FragLenObservations": "the C++ per-worker observation arrays",
    "frag_length_observations": "the scan result's fragment-length observation block",
    "_replay_fraglen_observations": "the Python-side replay of those arrays",
    "fraglen_obs": "the per-worker observation buffer",
    "n_frag_length_unambiguous": "site B's accepted-observation counter",
    "n_frag_length_ambiguous": "site B's discarded-observation counter",
    "get_unique_frag_length_mrna": "definition B itself — the transcript-space, unanimity-gated length",
}


#: A tombstone — a comment that records the deletion and why — is not a live reference, and this
#: codebase keeps them deliberately. ⚠ The exemption is a LITERAL, exact marker rather than "any
#: comment": a stale comment that merely *mentions* a deleted symbol is precisely how the next reader
#: concludes it still exists, and three such were found and fixed when this gate first ran.
TOMBSTONE = "DELETED by TRAPS: pure-and-length-censored"


def _sources() -> list[Path]:
    return [p for p in SRC.rglob("*") if p.suffix in {".py", ".cpp", ".h"} and p.is_file()]


@pytest.mark.parametrize("name", sorted(DELETED_BY_C2))
def test_the_scanners_fragment_length_machinery_is_GONE(name: str) -> None:
    """No occurrence of *name* survives anywhere under ``src/`` — comments included.

    ⚠ Comments are deliberately in scope. A dangling reference to a deleted subsystem is how the next
    reader concludes it still exists.
    """
    why = DELETED_BY_C2[name]
    pattern = re.compile(rf"\b{re.escape(name)}\b")
    hits = [
        f"  {path.relative_to(SRC)}:{i}: {line.strip()}"
        for path in _sources()
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if pattern.search(line) and TOMBSTONE not in line
    ]
    assert not hits, (
        f"`{name}` survives TRAPS: pure-and-length-censored in {len(hits)} place(s) — {why}\n"
        f"(a deletion record is exempt only on a line containing {TOMBSTONE!r})\n"
        + "\n".join(hits[:20])
    )


def test_the_singular_model_SURVIVES() -> None:
    """⛔ The converse, and it is not a formality.

    ``FragmentLengthModel`` and ``FragmentLengthModels`` differ by one character. The singular is the
    scorer — ``pipeline`` builds the RNA and gDNA scoring models and every transcript's effective
    length from it via ``from_pmf``. A search-and-delete on the plural that caught the singular would
    remove the fragment-length term from scoring altogether, and the suite would report it as a
    numerical drift rather than as a deletion.
    """
    from rigel.frag_length_model import FragmentLengthModel

    assert hasattr(FragmentLengthModel, "from_pmf")
    assert hasattr(FragmentLengthModel, "compute_all_transcript_eff_lens")


def test_the_deleted_names_are_gone_from_the_BUILT_extension_too() -> None:
    """The C++ is compiled, so a stale source name could still be shipping.

    ⚠ Greps over ``src/`` cannot see a stale build. This asserts the scan result itself no longer
    carries the observation block — the thing the Python side used to replay.
    """
    from rigel.native import BamScanner  # noqa: F401  — import is the check that the build is current

    proc = subprocess.run(
        ["git", "grep", "-n", "frag_length_observations", "--", "src/"],
        capture_output=True,
        text=True,
        cwd=SRC.parent,
    )
    assert proc.stdout.strip() == "", f"still referenced:\n{proc.stdout}"
