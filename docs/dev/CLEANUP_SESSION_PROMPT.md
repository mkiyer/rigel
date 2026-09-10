# THE CLEANUP SESSION — the kickoff prompt (prepared 2026-09-10; paste the block below into a new session)

⚠ A dev doc: the prompt itself, kept so it can be edited before it is used. Nothing cites it.

---

Start from `main` after the 2026-09-10 work is committed (the landscape prior's two landings, the
instrument re-runs, the doc updates). Read `CLAUDE.md` first, then `docs/dev/NEXT_SESSION.md`. Run
`python scripts/design/preflight.py` and the suite before touching anything; the baseline line in
`CLAUDE.md` is the count you must reproduce.

THIS SESSION IS A CODE REVIEW, A CODE CLEANUP AND A DOCUMENTATION CLEANUP. No mechanism changes, no
new numbers, no experiments. The calibration phase is functional and measured; the job now is to make
the pipeline solid, the code beautiful, and the documentation readable by a human. The owner's rulings:
no backwards compatibility, no legacy, no "kept for comparison"; stale code and stale prose are deleted,
not commented out; the changelog is git; converge and delete.

THE SIZE OF WHAT YOU ARE CLEANING (re-derive; these are the 2026-09-10 counts). Python source 38,300
lines in `src/rigel/` (the calibration package 13,900 lines in 39 modules across seven layers; the eight
largest files: `index.py` 1,795, `calibration/splice_graph.py` 1,635, `cli.py` 1,556, `sim/wgs_engine.py`
1,388, `sim/whole_genome.py` 1,304, `pipeline.py` 1,111, `estimator.py` 1,056, `sim/locus_sweep.py` 1,020)
plus 12,300 lines of C++ in `src/rigel/native/`. Tests 54,200 lines in 166 files. Scripts 31,600 lines,
69 files in `scripts/design/` alone. Docs 10,300 lines across nine permanent files (`DESIGN.md` 2,431,
`EQUATIONS.md` 1,766, `TRAPS.md` 1,607, `TESTING.md` 1,216, `ISSUES.md` ~990, `MANUAL.md` 1,054), plus
`CLAUDE.md` at 457. The owner has read `TRAPS.md`, `TESTING.md` and `EQUATIONS.md` and finds them
unreadable: a huge accumulation, overwhelming, written for the sessions that produced them rather than
for a reader.

THE ORDER, AND THE GATE AT THE END OF EACH PHASE. Every phase ends with: the suite at zero failures with
the collected count re-derived from `CLAUDE.md`'s table (never adjusted), `preflight.py --full`, `ruff
check src/ tests/ scripts/` and `ruff format src/ tests/` (never `scripts/`), and — for any phase that
touches `src/` — `rename_identity.py --check` against a reference frozen at the start of the session
(`--freeze` on the ladder's identity pair named in `NEXT_SESSION.md`; array content and the transcript
table, bit for bit). The owner commits at the end of each phase; you do not commit. Present each phase's
plan as a short list of what will be deleted, what will be rewritten and what will be kept, and wait for
the owner's go before executing it.

PHASE 0 — THE INVENTORY (no edits). Produce one census, as a table the owner can read in five minutes:
per permanent doc, its stated purpose, its length, how many of its sections are rulings or derivations
that code or tests depend on versus history, measurements and narrative; per `TRAPS.md` rule, whether it
is cited anywhere (`grep -o "TRAPS: [a-z0-9-]*"` across `src/ scripts/ tests/ docs/ CLAUDE.md` finds ~130
distinct names cited) and whether it still changes what a session does; per `scripts/design/` file,
whether its question is still open, answered-and-recorded (a refusal whose verdict lives in `ISSUES.md`),
or dead (the eight in `tests/test_scripts_index.py`'s `UNDOCUMENTED_DEBT`, and any whose consumer is
gone); per `src/` module, `module_census.py`'s flags (upward imports, docstrings naming a sibling with no
import, dead public surface) plus the docstrings that are histories rather than contracts (dates,
measurements, "this said X until 2026-08-17"); per test file, what it gates and whether that thing still
exists. End with a recommended cut list and a length budget per doc.

PHASE 1 — THE DOCUMENTATION. Each permanent doc gets a one-paragraph purpose at its top and is cut to
that purpose: `DESIGN.md` = the rulings and the design as built, each ruling with its one measurement and
its date, nothing narrated; `EQUATIONS.md` = the derivations the code depends on, each named from the
code or `DESIGN.md`, nothing else; `TRAPS.md` = rules that change what the next session does, one
paragraph each, cited by name — an uncited rule that no longer changes behaviour is deleted, a cited one is
rewritten to a paragraph a human can read; `TESTING.md` = a manual (how to build each panel and reference,
how to run each gate, what the suite can and cannot judge), not a history of the panels; `SUCCESS.md` =
how performance is judged, Stage A and Stage B, with the instruments; `ISSUES.md` = open entries cut to
their substance (the question, the priority, the instrument), and the CLOSED / REFUSED record compacted to
name, verdict and the killing number (the append-only rule keeps the numbers, not the prose);
`MANUAL.md` and `README.md` = a user's documents, checked against the CLI as it is. `CLAUDE.md` stays
directions, not history, and shrinks with the rest. Move, never copy: when a passage's substance survives
elsewhere, delete it here in the same edit. Git is the archive; nothing is moved to an "archive" folder.

PHASE 2 — THE SCRIPTS. Retire every `scripts/design/` instrument whose question is answered and whose
verdict is recorded in `ISSUES.md` or `DESIGN.md` (the instrument's docstring is not a permanent record;
the doc is), and every dead one. Each survivor keeps a docstring of one purpose paragraph and one usage
block. `CLAUDE.md`'s instrument table and `tests/test_scripts_index.py` shrink with them; the `--self-test`
sweep in `preflight.py --full` must still pass.

PHASE 3 — THE SOURCE, one calibration layer at a time from `rigel/calibration/_layers.py`, then the
pipeline, the scanner's Python side, the CLI, and the simulator. The review's questions per module: is
its docstring a contract (what it computes, what it guarantees, what it refuses) rather than a diary; is
every comment about the code as it is; is there a function no caller reaches; is there a parameter,
branch, flag or config field that exists for a mechanism that was deleted or refused; are two things
called by one name or one thing by two; is the module the right size for one concept (layer 4 is five
modules for one concept; `splice_graph.py` and `index.py` are the largest and the oldest). Measurements
and dates in `src/` comments move to `DESIGN.md` or are deleted; the source cites tests and the executable
references, never a doc. Every refactor is proved a numeric no-op by the identity check before the next
one starts; a refactor that is not a no-op is a mechanism change and is out of scope for this session.

PHASE 4 — THE TESTS. Delete gates for mechanisms that no longer exist, merge duplicated fixtures, and
make every remaining test file's docstring say what it gates in one paragraph. Keep the falsification
discipline: a gate that is kept must still fire when its mechanism is broken (spot-check by perturbation
where the file's docstring claims it).

RULES THAT STAND THROUGHOUT: work in the `rigel` conda env with `OMP_NUM_THREADS=1`; never bare `pytest`;
re-derive every collected count from `CLAUDE.md`'s table; the two xfails are executable records of
defects and are not touched; goldens are regenerated only after a column-by-column diff shows a
numeric no-op (a cleanup that moves a golden is a mechanism change); the deferred stratum stays reported;
the owner drives commits; if a deletion needs a judgement about whether something is still a ruling, ask,
in plain words, one question at a time. Finish by rewriting `docs/dev/NEXT_SESSION.md` as the state and
updating `CLAUDE.md`'s baseline line.
