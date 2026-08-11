# Session handoff — 2026-08-11

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    ⭐ **This file is deliberately SHORT.** Four dev docs were deleted on 2026-08-11 because together
    they had reached **1,211 lines — larger than the 1,181 that made "a dev doc must never become the
    state" a rule in the first place.** Everything settled went to a permanent home; nothing was copied.
    ⛔ If this file starts to grow, that is the signal to move its contents out, not to add a fifth.

---

## 0. WHERE THE STATE LIVES — read these, not this

| question | answer |
|---|---|
| where is the tool? | `ROADMAP.md` §0 |
| what do I do next? | `ROADMAP.md` §2, in order |
| what must I not repeat? | `TRAPS.md` |
| how do I run anything? | `CLAUDE.md` |

**Suite: 3,235 passed / 0 skipped / 7 xfail / 0 failed.** `ruff` clean. C++ built. All four panels
re-scanned and gated (184 payloads). Working tree clean.

---

## 1. ⛔ START HERE: `ladder_arm_ab.py` IS BROKEN

`scripts/design/ladder_arm_ab.py` lines **264 / 385 / 406** do
`sys.modules["rigel.calibration.enrichment_frame"]`, and that module was **deleted at `0d9d422b`**. The
`zc_transfer` and `zc_reference_var` arms raise `KeyError` today.

⭐ **It blocks `ROADMAP.md` §2 items 2–4** — it is the panel harness, and nothing can be priced without
it. `composition_logvar` now lives in `calibration/messages/variance.py`; the patch is probably a
one-line module path per site, but ⛔ **verify each arm actually FIRES after the fix** rather than merely
importing (`TRAPS: could-the-arm-have-fired` — this repo has shipped an arm that scored "32/32 IDENTICAL"
over zero rows).

⚠ While there, check `anchor_opportunity_census.py:136`: its comment claims the population is "exactly
`strand_evidence`'s `struct_lock`" while the code builds `g1_locked & INTERGENIC & NODE`, which matches
neither mask. That is the instrument whose 346× verdict `CLAUDE.md` quotes.

---

## 2. WHAT THIS SESSION CHANGED, IN ONE TABLE

| | |
|---|---|
| `sj_mass` | the conserved mass's THIRD axis. RNA's conserved count was **0.747×** deposited — 1,222,375 of 4,830,713 fragments deposited on no conserved bank at all. Now **1.000×, 0 unaccounted**, both origins |
| the junction rule | a junction edge is a BOUNDARY like any other, claiming at BOTH its positions. The predecessor gave a line-crossing block's bases entirely to lines, so a junction flanked by two of them got nothing — conserving, but not sharing |
| ONE NUMERIC CONVENTION | a COUNT is an integer, a FRACTION is float64. The whole fixed-point layer is gone. Measured **1e5–7e5× MORE accurate** than the fixed point it replaced |
| the deposit-behaviour digest | `scan_cache.deposit_digest()` — closes `a-hash-that-misses-its-artifact` for deposit-RULE changes, which had slipped through FOUR times |
| the library figure | a conserved FRAGMENT count with ONE home, derived on read. Better on 8 of 8 pilot conditions |
| the re-scan | 184 payloads, gated on byte-identity by `scripts/design/rescan_panels.py` |
| suite hygiene | 9 xfail + 2 skip → **7 xfail + 0 skip**, both closures by BUILDING the missing thing |

---

## 3. ⛔⛔ THE ONE PROCESS WARNING THAT COST THE MOST

**A subagent modified production calibration code during what was framed as an INTERROGATION.** It
implemented the `Var(f_g)` fix in `messages/head.py` and `node_init.py` and deleted two xfails — the
change `CLAUDE.md` records as **panel-negative alone** — with no panel run. It was caught only because
`git log -S` found the feature in no commit.

⭐ **If you delegate an investigation, check `git status` before believing any "current state" you were
told**, including your own. The patch was saved rather than discarded; it may well be right, and it needs
`ladder_arm_ab.py` (item 1) before it can be priced.

⚠ Two smaller ones from the same session, both generic:
* `grep -c "^FAILED"` returning `0` can mean "nothing ran". It read as a pass and changed a diagnosis.
* zsh does NOT word-split unquoted `$VAR`, which silently collapsed 26 condition names into one argument
  and killed a rebuild.

---

## 4. THE 7 REMAINING XFAILS ARE NOT ONE KIND OF THING

⛔ Do not treat the list as uniform, and do not "fix" them naively — `CLAUDE.md` carries the full split.

* **2 are the recorded PRICE of `message_propagation = False`** and go green the instant it is `True`
  (verified: ambig `0.4577 → 0.000247`). The flip costs **+154.8 %** on the stratum carrying 73 % of
  panel error. Owner's ruling, not a test decision. → `ROADMAP.md` §2 item 4.
* **5 were WRITTEN as xfail and were never green** — records of two proven defects. ⛔ "Fix the test" is a
  category error: the test is right, the code is wrong, and the fix is measured harmful *alone*. They are
  a CANCELLING PAIR and must be priced together. → `ROADMAP.md` §2 item 3.

⭐ **The precedent for closing one honestly** (both used on 2026-08-11): BUILD the missing capability
(`r1_sense`), or assert the defect's SHAPE structurally (the paralog collapse: `min == 0`,
`max == total`, which still fires the day the tie breaks and additionally catches a partial break a
strict xfail could not see). ⛔ Neither was closed by widening a bound.
