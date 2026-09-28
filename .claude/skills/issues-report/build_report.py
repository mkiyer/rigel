#!/usr/bin/env python
"""Build the Rigel Open Issues page from docs/ISSUES.md and docs/ROADMAP.md.

The page is a rendering of the issue log and never a second home for it: it is regenerated after every
ISSUES.md change and never edited by hand. It runs nothing and measures nothing; it reads the two docs
and `git` for the stamp. Standard library only: the entry bodies go through the small markdown
converter below (paragraphs, lists, bold, italics, code, links, tables, fences).

It fails loudly on what would make the page lie about the log: an OPEN entry without a vocabulary
priority (now / next / later / parked) directly under its heading, a name that is not kebab-case, a name
used twice, a `## ` heading other than OPEN or CLOSED (every entry after it would drop off the page), a
heading line inside a section that is not `### name` (its entry would fold into the one above), a
ROADMAP.md with no next steps, and a citation in ROADMAP.md's next steps to a name no entry holds. Other
dangling `ISSUES:` citations and an OPEN section out of tier order are printed as warnings.
"""

from __future__ import annotations

import argparse
import html
import re
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SKILL = Path(__file__).resolve().parent
TIERS = ("now", "next", "later", "parked")
PLACEHOLDERS = ("__CHIPS__", "__MAIN__", "__STAMP__")

KEBAB = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
NAME = r"[a-z0-9]+(?:-[a-z0-9]+)*"
DATE = re.compile(r"\d{4}-\d{2}(?:-\d{2})?")

# ── the markdown converter ───────────────────────────────────────────────────────────────────────

BULLET = re.compile(r"^( {0,3})([-*+]) +(.*)$")
ORDERED = re.compile(r"^( {0,3})(\d{1,4})[.)] +(.*)$")
FENCE = re.compile(r"^ {0,3}(```|~~~)")
HEADING = re.compile(r"^ {0,3}(#{1,6}) +(.*?)[ #]*$")
RULE = re.compile(r"^ {0,3}([-*_])(?: *\1){2,} *$")
TABLE_SEP = re.compile(r"^ *\|? *:?-{3,}:? *(?:\| *:?-{3,}:? *)*\|? *$")
#: a line that starts its own paragraph inside a run of text: the log's `(a)` enumerators, its circled
#: numbers, and the `Instrument:` line that closes an entry. Plain markdown would fold them into the
#: paragraph above; the page keeps them apart because the log is written that way.
OWN_LINE = re.compile(r"^(?:\([a-z]\) |[①-⑳] |Instruments?: )")
INSTRUMENT = re.compile(r"^(Instruments?):\s*")
#: a short line that opens with a capitalised word and ends in a colon (`BATCHED:`, `LATENT (no current
#: path reaches it):`) is a label for what follows, so it stands alone
LABEL = re.compile(r"^[A-Z]{3,}\b(?![.'’])[^\n]{0,56}:$")

CODE_SPAN = re.compile(r"(`+)(.+?)\1")
CODE_CITE = re.compile(rf"^ISSUES: ({NAME})\b")
PLAIN_CITE = re.compile(rf"\bISSUES: ({NAME})")
BOLD = re.compile(r"\*\*(.+?)\*\*")
ITALIC = re.compile(r"(?<![\w*])\*(?![\s*])([^*]+?)(?<!\s)\*(?![\w*])")
LINK = re.compile(r"\[([^\]]+)\]\(([^)\s]+)\)")
SLOT = re.compile(r"\x00(\d+)\x00")


def esc(s: str) -> str:
    return html.escape(s, quote=True)


def indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


class Markdown:
    """Block and inline rendering, with `ISSUES: name` turned into in-page links.

    `known` maps every entry name to its section; a citation to a name not in it is kept as text,
    marked, and collected in `dangling` so the build can report it.
    """

    def __init__(self, known: dict[str, str]):
        self.known = known
        self.dangling: dict[str, int] = {}

    # inline
    def xref(self, name: str, label: str, code: bool) -> str:
        inner = f"<code>{label}</code>" if code else label
        if name in self.known:
            return f'<a class="xref {self.known[name]}" href="#{name}">{inner}</a>'
        self.dangling[name] = self.dangling.get(name, 0) + 1
        return f'<span class="dangling" title="no entry by this name">{inner}</span>'

    def inline(self, text: str) -> str:
        slots: list[str] = []

        def keep(fragment: str) -> str:
            slots.append(fragment)
            return f"\x00{len(slots) - 1}\x00"

        def code(m: re.Match) -> str:
            body = m[2]
            if len(body) > 2 and body[0] == body[-1] == " ":
                body = body[1:-1]
            c = CODE_CITE.match(body)
            if c:
                return keep(self.xref(c[1], esc(body), code=True))
            return keep(f"<code>{esc(body)}</code>")

        def link(m: re.Match) -> str:
            href = m[2]
            if href.startswith(("http://", "https://")):
                return keep(f'<a href="{href}" target="_blank" rel="noopener">{m[1]}</a>')
            return m[1]

        out = esc(CODE_SPAN.sub(code, text))
        out = LINK.sub(link, out)
        out = PLAIN_CITE.sub(lambda m: keep(self.xref(m[1], m[0], code=False)), out)
        out = BOLD.sub(r"<strong>\1</strong>", out)
        out = ITALIC.sub(r"<em>\1</em>", out)
        while SLOT.search(out):
            out = SLOT.sub(lambda m: slots[int(m[1])], out)
        return out

    # blocks
    def is_table(self, lines: list[str], i: int) -> bool:
        return "|" in lines[i] and i + 1 < len(lines) and bool(TABLE_SEP.match(lines[i + 1]))

    def interrupts(self, lines: list[str], i: int) -> bool:
        """Does line `i` end the paragraph running above it?"""
        line = lines[i]
        if FENCE.match(line) or HEADING.match(line) or RULE.match(line) or BULLET.match(line):
            return True
        if OWN_LINE.match(line) or LABEL.match(line.strip()) or self.is_table(lines, i):
            return True
        m = ORDERED.match(line)
        # an ordered item breaks a paragraph only at 1 or after a colon, so a wrapped line that
        # happens to start "429. The ladder's…" stays inside its sentence
        return bool(m) and (m[2] == "1" or lines[i - 1].rstrip().endswith(":"))

    def table(self, rows: list[str]) -> str:
        def cells(row: str) -> list[str]:
            row = row.strip()
            row = row[1:] if row.startswith("|") else row
            row = row[:-1] if row.endswith("|") and not row.endswith("\\|") else row
            return [c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", row)]

        head = "".join(f"<th>{self.inline(c)}</th>" for c in cells(rows[0]))
        body = "".join(
            "<tr>" + "".join(f"<td>{self.inline(c)}</td>" for c in cells(r)) + "</tr>"
            for r in rows[2:]
        )
        return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'

    def list_block(self, lines: list[str], i: int) -> tuple[str, int]:
        pattern = ORDERED if ORDERED.match(lines[i]) else BULLET
        base = len(pattern.match(lines[i])[1])
        items: list[list[str]] = []
        start = None
        loose = False
        n = len(lines)
        while i < n:
            m = pattern.match(lines[i])
            if not m or len(m[1]) != base:
                break
            if pattern is ORDERED and start is None:
                start = int(m[2])
            offset = len(lines[i]) - len(m[3])
            item = [m[3]]
            i += 1
            while i < n:
                if not lines[i].strip():
                    j = i
                    while j < n and not lines[j].strip():
                        j += 1
                    if j < n and indent(lines[j]) >= base + 2 and not pattern.match(lines[j]):
                        item.extend([""] * (j - i))
                        loose = True
                        i = j
                        continue
                    break
                # an unindented line ends the list: the log writes its closing sentence flush left
                if indent(lines[i]) >= base + 2:
                    item.append(lines[i][min(indent(lines[i]), offset) :])
                    i += 1
                    continue
                break
            items.append(item)
            j = i
            while j < n and not lines[j].strip():
                j += 1
            if j > i and j < n:
                m2 = pattern.match(lines[j])
                if m2 and len(m2[1]) == base:
                    loose = True
                    i = j
                    continue
                break
        lis = []
        for item in items:
            inner = self.blocks(item)
            if not loose and inner.startswith("<p>") and inner.count("<p") == 1:
                inner = inner[3:].replace("</p>", "", 1)
            lis.append(f"<li>{inner}</li>")
        if pattern is ORDERED:
            attr = f' start="{start}"' if start not in (None, 1) else ""
            return f"<ol{attr}>{''.join(lis)}</ol>", i
        return f"<ul>{''.join(lis)}</ul>", i

    def blocks(self, lines: list[str]) -> str:
        out: list[str] = []
        i, n = 0, len(lines)
        while i < n:
            line = lines[i]
            if not line.strip():
                i += 1
                continue
            if FENCE.match(line):
                j = i + 1
                while j < n and not FENCE.match(lines[j]):
                    j += 1
                out.append(
                    f'<pre class="scroll"><code>{esc(chr(10).join(lines[i + 1 : j]))}</code></pre>'
                )
                i = j + 1
                continue
            if RULE.match(line):
                out.append("<hr>")
                i += 1
                continue
            m = HEADING.match(line)
            if m:
                out.append(f'<h4 class="md-h">{self.inline(m[2])}</h4>')
                i += 1
                continue
            if self.is_table(lines, i):
                j = i + 2
                while j < n and "|" in lines[j] and lines[j].strip():
                    j += 1
                out.append(self.table(lines[i:j]))
                i = j
                continue
            if BULLET.match(line) or ORDERED.match(line):
                fragment, i = self.list_block(lines, i)
                out.append(fragment)
                continue
            para = [line.strip()]
            i += 1
            if LABEL.match(para[0]):
                out.append(f'<p class="label">{self.inline(para[0])}</p>')
                continue
            while i < n and lines[i].strip() and not self.interrupts(lines, i):
                para.append(lines[i].strip())
                i += 1
            text = " ".join(para)
            m = INSTRUMENT.match(text)
            if m:
                rest = self.inline(text[m.end() :])
                out.append(f'<p class="instrument"><span class="lbl">{m[1]}</span> {rest}</p>')
            else:
                out.append(f"<p>{self.inline(text)}</p>")
        return "\n".join(out)


# ── the log ──────────────────────────────────────────────────────────────────────────────────────


def fail(problems: list[str]) -> None:
    raise SystemExit(
        "⛔ the issue log cannot be rendered:\n" + "\n".join(f"  - {p}" for p in problems)
    )


def trim_tail(body: list[str]) -> list[str]:
    while body and (not body[-1].strip() or RULE.match(body[-1])):
        body.pop()
    return body


def trim(body: list[str]) -> list[str]:
    body = trim_tail(body)
    while body and not body[0].strip():
        body.pop(0)
    return body


def parse_issues(path: Path) -> tuple[list[dict], list[dict]]:
    """The OPEN and CLOSED / REFUSED entries, in file order."""
    entries: list[dict] = []
    problems: list[str] = []
    section = None
    seen_sections = set()
    current = None
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if line.startswith("## "):
            title = line[3:].strip()
            section = (
                "open"
                if title.startswith("OPEN")
                else ("closed" if title.startswith("CLOSED") else None)
            )
            if section is None:
                problems.append(
                    f"line {number}: '## {title}' is neither OPEN nor CLOSED, so every entry after it "
                    "would drop off the page"
                )
            seen_sections.add(section)
            current = None
            continue
        m = re.match(r"^### (\S.*?)\s*$", line)
        if m and section:
            current = {"name": m[1], "section": section, "line": number, "body": []}
            entries.append(current)
            continue
        if section and re.match(r"^#{3,}", line):
            problems.append(
                f"line {number}: '{line.strip()[:40]}' is not an entry heading (`### name`), so it "
                "would fold into the entry above"
            )
        if current is not None:
            current["body"].append(line)

    problems += [
        f"no '## {s.upper()}' section in {path.name}"
        for s in ("open", "closed")
        if s not in seen_sections
    ]
    counts: dict[str, list[int]] = {}
    for e in entries:
        # an OPEN entry's priority line must sit directly under its heading, so only its tail is trimmed
        e["body"] = trim_tail(e["body"]) if e["section"] == "open" else trim(e["body"])
        counts.setdefault(e["name"], []).append(e["line"])
        if not KEBAB.match(e["name"]):
            problems.append(f"line {e['line']}: '### {e['name']}' is not a kebab-case name")
    problems += [
        f"'{name}' is used {len(at)} times (lines {', '.join(map(str, at))})"
        for name, at in counts.items()
        if len(at) > 1
    ]

    open_entries = [e for e in entries if e["section"] == "open"]
    for e in open_entries:
        problem = read_priority(e)
        if problem:
            problems.append(f"line {e['line']}: '{e['name']}' {problem}")
    if problems:
        fail(problems)
    return open_entries, [e for e in entries if e["section"] == "closed"]


def read_priority(e: dict) -> str | None:
    """Split the entry's `priority: WORD — reason · kind: K · date` line off its body."""
    body = e["body"]
    if not body or not body[0].startswith("`priority:"):
        return "has no `priority:` line directly under its heading"
    k = 1
    first = body[0].rstrip()
    while not first.endswith("`") and k < len(body) and body[k].strip():
        first += " " + body[k].strip()
        k += 1
    if not first.endswith("`"):
        return "has a priority line with no closing backtick"
    m = re.match(r"^priority: (\S+?)(?= — | · |$)(.*)$", first[1:-1])
    if not m or m[1] not in TIERS:
        word = m[1] if m else first[:30]
        return f"has priority '{word}', which is not one of {' / '.join(TIERS)}"
    parts = m[2].split(" · ")
    at = next((i for i, p in enumerate(parts) if p.startswith("kind: ")), None)
    reason = " · ".join(parts if at is None else parts[:at])
    e["tier"] = m[1]
    e["reason"] = re.sub(r"^\s*[—–-]\s*", "", reason).strip()
    e["kind"] = parts[at][len("kind: ") :].strip() if at is not None else ""
    e["dated"] = " · ".join(parts[at + 1 :]).strip() if at is not None else ""
    e["body"] = trim(body[k:])
    return None


SENTENCE_END = re.compile(r"(?<=[.!?])\s+(?=\S)")
ABBREVIATIONS = ("e.g.", "i.e.", "vs.", "cf.", "etc.")


def first_sentence(lines: list[str], cap: int = 300) -> str:
    """The first sentence of an entry's first paragraph: a CLOSED entry's verdict line."""
    para: list[str] = []
    for line in lines:
        if not line.strip() or (para and (BULLET.match(line) or line.startswith("|"))):
            break
        para.append(line.strip())
    text = " ".join(para)
    for m in SENTENCE_END.finditer(text):
        head = text[: m.start()]
        if len(head) >= 25 and not head.endswith(ABBREVIATIONS) and head.count("`") % 2 == 0:
            text = head
            break
    if len(text) > cap:
        text = text[:cap].rsplit(" ", 1)[0]
        if text.count("`") % 2:
            text = text[: text.rfind("`")].rstrip()
        text += " …"
    return text


GROUP = re.compile(r"^\*\*(.+?)\*\*")
ITEM = re.compile(r"^(?:\d+\.|[-*+]) +(.*)$")


def cites_in(text: str, names: set[str]) -> list[str]:
    """The entries a step names, in order: every `ISSUES: name`, and a bare `name` in code that an entry holds."""
    cites: list[str] = []
    for c in re.finditer(rf"ISSUES: ({NAME})|`({NAME})`", text):
        name = c[1] or c[2]
        if (c[1] or name in names) and name not in cites:
            cites.append(name)
    return cites


def parse_roadmap(path: Path, names: set[str]) -> tuple[str, list[dict], list[str]]:
    """ROADMAP.md's next-steps section: its heading, its groups of steps, its lines.

    A line opening in bold at column 0 starts a group (`**The local track.**`); each top-level item under
    it (`1. …` or `- …`) is a step, and a group with no items is one step holding its paragraph's
    citations. The steps are numbered in sequence across the groups, so every part of the order shows.
    """
    lines = path.read_text().splitlines()
    start = next((i for i, line in enumerate(lines) if re.match(r"^## Next\b", line)), None)
    if start is None:
        fail([f"{path.name} has no '## Next' section, so the page has no order to show"])
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")), len(lines))
    section = lines[start + 1 : end]
    dangling = sorted(
        {m[1] for m in re.finditer(rf"ISSUES: ({NAME})", "\n".join(section)) if m[1] not in names}
    )
    if dangling:
        fail(
            [f"{path.name}'s next steps cite `ISSUES: {n}`, which no entry holds" for n in dangling]
        )

    groups: list[dict] = []
    i, n = 0, len(section)
    while i < n:
        line = section[i]
        g = GROUP.match(line)
        if g:
            para = [line.strip()]
            i += 1
            while i < n and section[i].strip() and not ITEM.match(section[i]):
                para.append(section[i].strip())
                i += 1
            groups.append({"title": g[1].strip().rstrip("."), "para": " ".join(para), "steps": []})
            continue
        m = ITEM.match(line)
        if not m:
            i += 1
            continue
        if not groups:
            groups.append({"title": "", "para": "", "steps": []})
        item = [m[1]]
        i += 1
        while i < n and section[i].startswith(" ") and section[i].strip():
            item.append(section[i].strip())
            i += 1
        text = " ".join(item)
        bold = re.match(r"^\*\*(.+?)\*\*", text)
        groups[-1]["steps"].append(
            {
                "title": bold[1] if bold else first_sentence([text], 120),
                "cites": cites_in(text, names),
            }
        )
    k = 0
    for g in groups:
        if not g["steps"]:
            # a group written as one paragraph is one step: its entries, under the group's heading
            g["steps"].append({"title": "", "cites": cites_in(g["para"], names)})
        for step in g["steps"]:
            k += 1
            step["n"] = k
    if not k:
        fail([f"{path.name}'s '{lines[start][3:].strip()}' section has no steps"])
    return lines[start][3:].strip(), groups, section


def git_stamp(paths: list[Path]) -> tuple[str, list[str]]:
    """The HEAD commit, and which of `paths` differ from it in the working tree."""
    try:
        sha = subprocess.run(
            ["git", "-C", str(REPO), "rev-parse", "--short=8", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        rel = [str(p.relative_to(REPO)) for p in paths if p.is_relative_to(REPO)]
        status = subprocess.run(
            ["git", "-C", str(REPO), "status", "--porcelain", "--", *rel],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return "unknown", []
    return sha, [Path(line[3:]).name for line in status.splitlines() if line.strip()]


# ── the page ─────────────────────────────────────────────────────────────────────────────────────


def render_issue(e: dict, md: Markdown, step_of: dict[str, list[int]]) -> str:
    date = DATE.search(e["dated"])
    meta = []
    if e["kind"]:
        meta.append(f'<span class="chip">{esc(e["kind"])}</span>')
    if e["name"] in step_of:
        at = " · ".join(str(n) for n in step_of[e["name"]])
        meta.append(f'<a class="chip step" href="#order">step {at}</a>')
    if date:
        meta.append(f'<span class="s-date num">{date[0]}</span>')
    dated = f'<span class="dated">{md.inline(e["dated"])}</span>' if e["dated"] else ""
    cite = f"ISSUES: {e['name']}"
    return (
        f'<details class="issue" id="{e["name"]}" data-tier="{e["tier"]}">'
        f'<summary><span class="chev" aria-hidden="true"></span>'
        f'<span class="s-main"><span class="s-name">{e["name"]}</span>'
        f'<span class="s-why">{md.inline(e["reason"])}</span></span>'
        f'<span class="s-meta">{"".join(meta)}</span></summary>'
        f'<div class="entry"><div class="md">{md.blocks(e["body"])}</div>'
        f'<div class="entry-foot"><code class="cite">{cite}</code>'
        f'<button type="button" class="btn copy" data-cite="{cite}">Copy citation</button>'
        f'{dated}<span class="src num">docs/ISSUES.md:{e["line"]}</span></div></div></details>'
    )


def render_closed(e: dict, md: Markdown) -> str:
    sentence = first_sentence(e["body"])
    verdict = re.match(r"^([A-Z]{4,})\b", sentence) or re.search(
        r"\b(closed|refused|refuted|ruled|retired|moved)\b", sentence, re.I
    )
    date = DATE.search(sentence)
    chips = f'<span class="chip">{verdict[1].upper()}</span>' if verdict else ""
    chips += f'<span class="s-date num">{date[0]}</span>' if date else ""
    return (
        f'<li class="closed" id="{e["name"]}"><div class="c-head"><span class="c-name">'
        f'{e["name"]}</span>{chips}</div><p class="c-why">{md.inline(sentence)}</p></li>'
    )


def build(issues: Path, roadmap: Path) -> tuple[str, dict]:
    open_entries, closed = parse_issues(issues)
    known = {e["name"]: "open" for e in open_entries} | {e["name"]: "closed" for e in closed}
    md = Markdown(known)
    heading, step_groups, section = parse_roadmap(roadmap, set(known))
    steps = [s for g in step_groups for s in g["steps"]]
    # every step that names an entry, so an entry the order touches twice shows both
    step_of: dict[str, list[int]] = {}
    for s in steps:
        for name in s["cites"]:
            at = step_of.setdefault(name, [])
            if s["n"] not in at:
                at.append(s["n"])

    warnings = []
    ranks = [TIERS.index(e["tier"]) for e in open_entries]
    for a, b, e in zip(ranks, ranks[1:], open_entries[1:], strict=False):
        if b < a:
            warnings.append(
                f"'{e['name']}' ({e['tier']}) sits below a {TIERS[a]} entry: the OPEN "
                "section is out of tier order"
            )
    warnings += [
        f"'{e['name']}' has no `kind:` in its priority line" for e in open_entries if not e["kind"]
    ]

    by_tier = {t: [e for e in open_entries if e["tier"] == t] for t in TIERS}
    tiles = "".join(
        f'<a class="tile" href="#tier-{t}" data-tier="{t}"><span class="tile-l">'
        f'<span class="dot" aria-hidden="true"></span>{t.capitalize()}</span>'
        f'<span class="tile-v num">{len(by_tier[t])}</span></a>'
        for t in TIERS
    ) + (
        f'<a class="tile muted" href="#closed"><span class="tile-l">Closed / refused</span>'
        f'<span class="tile-v num">{len(closed)}</span></a>'
    )

    def chip_for(name: str) -> str:
        where = next((e["tier"] for e in open_entries if e["name"] == name), "closed")
        return (
            f'<a class="xchip" href="#{name}" data-tier="{where}">'
            f'<span class="dot" aria-hidden="true"></span>{name}</a>'
        )

    def step_item(s: dict) -> str:
        return (
            f'<li value="{s["n"]}"><span class="step-n num">{s["n"]}</span><div class="step-b">'
            f'<div class="step-t">{md.inline(s["title"])}</div>'
            f'<div class="step-c">{"".join(chip_for(n) for n in s["cites"] if n in known)}</div>'
            f"</div></li>"
        )

    step_items = "".join(
        (f'<h3 class="step-group">{md.inline(g["title"])}</h3>' if g["title"] else "")
        + f'<ol class="steps" start="{g["steps"][0]["n"]}">'
        + "".join(step_item(s) for s in g["steps"])
        + "</ol>"
        for g in step_groups
    )
    order = (
        f'<section class="order" id="order"><div class="sec-h"><h2>The order</h2>'
        f'<span class="sec-note">from ROADMAP.md, “{esc(heading)}”</span></div>'
        f"{step_items}"
        f'<details class="fulltext"><summary>The section as written</summary>'
        f'<div class="md">{md.blocks(section)}</div></details></section>'
    )
    groups = "".join(
        f'<section class="tier" id="tier-{t}" data-tier="{t}"><div class="sec-h">'
        f'<span class="dot" aria-hidden="true"></span><h2>{t.capitalize()}</h2>'
        f'<span class="count num" data-n="{len(by_tier[t])}">{len(by_tier[t])}</span></div>'
        f'<div class="issues">{"".join(render_issue(e, md, step_of) for e in by_tier[t])}</div>'
        f"</section>"
        for t in TIERS
        if by_tier[t]
    )
    closed_html = (
        f'<section id="closed"><details class="closed-wrap"><summary>'
        f'<span class="chev" aria-hidden="true"></span><h2>Closed / refused</h2>'
        f'<span class="count num" data-n="{len(closed)}">{len(closed)}</span>'
        f'<span class="sec-note">measured and turned down: read before proposing a mechanism'
        f"</span></summary>"
        f'<ul class="closed-items">{"".join(render_closed(e, md) for e in closed)}</ul>'
        f"</details></section>"
    )
    empty = (
        '<p class="empty" id="empty" hidden>Nothing matches the filter. Clear it to see '
        "every entry.</p>"
    )
    main = (
        f'<nav class="tiles" aria-label="Entries by priority">{tiles}</nav>{order}{empty}{groups}'
    )
    main += closed_html

    warnings += [
        f"`ISSUES: {name}` is cited {k}× but no entry has that name"
        for name, k in sorted(md.dangling.items())
    ]
    sha, dirty = git_stamp([issues, roadmap])
    built = time.strftime("%Y-%m-%d")
    edits = f" + uncommitted edits to {', '.join(dirty)}" if dirty else ""
    chips = "".join(
        f'<span class="chip">{c}</span>'
        for c in (
            f"commit {sha}{edits}",
            f"built {built}",
            f"{len(open_entries)} open",
            f"{len(closed)} closed",
        )
    )
    stamp = (
        f"Built {built} from <code>docs/ISSUES.md</code> and <code>docs/ROADMAP.md</code> at "
        f"commit <code>{sha}</code>{esc(edits)} by "
        "<code>.claude/skills/issues-report/build_report.py</code>."
    )
    tpl = (SKILL / "report_template.html").read_text()
    missing = [p for p in PLACEHOLDERS if tpl.count(p) != 1]
    if missing:
        fail([f"report_template.html must hold {p} exactly once" for p in missing])
    page = tpl.replace("__CHIPS__", chips).replace("__STAMP__", stamp).replace("__MAIN__", main)

    # every entry is on the page exactly once, as its own block
    placed = {n: page.count(f'id="{n}"') for n in known}
    lost = [f"'{n}' is on the page {k} times as a block" for n, k in placed.items() if k != 1]
    if lost:
        fail(lost)
    summary = {
        "tiers": {t: len(by_tier[t]) for t in TIERS},
        "closed": len(closed),
        "steps": len(steps),
        "warnings": warnings,
    }
    return page, summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--issues", type=Path, default=REPO / "docs" / "ISSUES.md")
    ap.add_argument("--roadmap", type=Path, default=REPO / "docs" / "ROADMAP.md")
    ap.add_argument("--html", type=Path, required=True, help="where to write the page")
    args = ap.parse_args()

    page, s = build(args.issues.resolve(), args.roadmap.resolve())
    args.html.parent.mkdir(parents=True, exist_ok=True)
    args.html.write_text(page)
    for w in s["warnings"]:
        print(f"  ⚠ {w}", file=sys.stderr)
    counts = " / ".join(f"{t} {k}" for t, k in s["tiers"].items())
    print(
        f"  ⭐ page -> {args.html}  ({args.html.stat().st_size:,} bytes; {counts}; "
        f"{s['closed']} closed; {s['steps']} roadmap steps)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
