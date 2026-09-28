---
name: issues-report
description: Build and publish the "Rigel Open Issues" page, generated from docs/ISSUES.md (every OPEN entry grouped now / next / later / parked, each collapsible, with a filter; the CLOSED / REFUSED names at the bottom) and docs/ROADMAP.md's next-steps order, to one Artifact that updates in place. Use after any change to docs/ISSUES.md or to ROADMAP.md's next steps, or when the owner asks for the open issues page.
---

# The open issues page

The page is GENERATED from the issue log and never edited by hand (owner, 2026-09-28). To change what
it says, edit `docs/ISSUES.md` or `docs/ROADMAP.md`, then rebuild and republish. The script runs nothing
and measures nothing: it reads the two docs, and `git` for the stamp.

## When

After every change to `docs/ISSUES.md`, and after any change to `ROADMAP.md`'s `## Next` section. A page
older than the log shows a stale state and looks current.

## 1. Build

```bash
python .claude/skills/issues-report/build_report.py --html <scratchpad>/rigel_open_issues.html
```

`<scratchpad>` is the session's scratchpad directory, where an artifact's page file belongs.

Standard library only; any Python 3.10+ runs it. It prints the counts per priority and any warnings.

It REFUSES to write the page, and names each problem, when:
- an OPEN entry has no `priority:` line directly under its heading (a blank line between counts as none),
  or its word is not exactly one of `now` / `next` / `later` / `parked` followed by ` — `, ` · ` or the
  closing backtick;
- a name is not kebab-case, or the same name heads two entries (OPEN and CLOSED count together);
- a `## ` heading in `ISSUES.md` is neither OPEN nor CLOSED, since every entry after it would drop off
  the page;
- a line inside a section starts `###` (or more) but is not `### name`, since its entry would fold into
  the one above;
- `ROADMAP.md` has no `## Next` section with steps, or that section cites `ISSUES: name` for a name no
  entry holds;
- an entry would appear on the page other than exactly once.

It WARNS, and still writes the page, on a citation `ISSUES: name` in `ISSUES.md` to a name no entry
holds (shown on the page with a wavy underline), an OPEN entry out of tier order, and a priority line
without `kind:`.
Fix a warning in the log, not in the page.

## 2. Publish, to the SAME artifact

**`https://claude.ai/artifact/Ah9RoZaM1sHAHoTXWRseFC`**

```
Artifact(action="publish", url="https://claude.ai/artifact/Ah9RoZaM1sHAHoTXWRseFC",
         file_path="<scratchpad>/rigel_open_issues.html")
```

If this conversation has not published or read it yet, read it first (`action="read"`). Publishing
without `url` creates a second page, and the owner is left with two links and no way to tell which is
current.

## What the page reads

| from | what |
|---|---|
| `## OPEN` in `ISSUES.md` | each `### name`, then the `` `priority: WORD — reason · kind: K · date` `` line, then the body, rendered as markdown. `ISSUES: name` citations become in-page links |
| `## CLOSED / REFUSED` | each name, a verdict chip (the leading CLOSED / REFUSED / RULED word) and the first sentence of its body |
| `## Next` in `ROADMAP.md` | the order: each column-0 **bold** line is a group heading, each top-level `1.` or `- ` item under it a step (its bold title and the entries it cites), and a group written as one paragraph is one step holding its citations. Steps are numbered in sequence across the groups; an entry carries a `step N` chip for every step that cites it (`step 4 · 5`) |

The markdown converter is the script's own (paragraphs, lists, bold, italics, code, links, tables,
fences). Two rules follow the log's own conventions rather than plain markdown: a line starting `(a)`,
a circled number, `Instrument:`, or a short capitalised label ending in a colon (`BATCHED:`) starts its own
paragraph; and an unindented line after a list item ends the list.

## Files

| | |
|---|---|
| `build_report.py` | parses the two docs, runs the checks above, fills the template |
| `report_template.html` | the page: tokens for both themes, the filter, expand/collapse, citation links, copy-citation. `__CHIPS__`, `__MAIN__` and `__STAMP__` are the placeholders |

The visual language is `ladder-report`'s (IBM Plex, the same neutrals). The priority colours are fixed by
tier, not by rank: now red, next amber, later blue, parked grey.
