# NEXT SESSION — the MESSAGE LAYER, developed on the test chromosome

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand and how to start, not what is
    settled — rulings are `DESIGN.md`, the ranked list is `ROADMAP.md`, the open problems are
    `ISSUES.md`. MOVE anything that settles into those and DELETE this file.

## The thread, and why it is next

**The message layer is the largest measured in-scope defect, and it ships ON.** Two independent
routes arrived at it on 2026-09-01:

* the calibration walk — at the worst in-scope condition the SHIPPED tool is **456,838** fragments of
  error against **298,597** with messages off. Messages ADD **+158,241 (+53 %)**; at `g98 ss.99 OFF`
  they add **+81,839 (+65 %)**. Re-derive with `calibration_walk.py --condition <name>` (rung E vs F);
* the vertex ceiling — handing every reachable vertex-truth object the ORACLE answer and re-solving
  the whole chain nets **≈0** (boundary −20,208, region +18,887). **A perfect local answer does not
  survive propagation**, which is why no local repair — coordinate, prior, or atom — can pay off until
  this is fixed. The arcsine coordinate thread was run to a verdict and refused on the same day;
  `ISSUES.md`'s CLOSED / REFUSED carries it, and ⛔ do not re-derive it.

⭐ The relay is not simply bad — its value is real and CONCENTRATED where the local solve is blind,
and its harm is concentrated where the destination already has good own evidence. That is the whole
shape of the problem, and it is exactly the charter's failure mode: **"minimal harm on stranded" is
currently not met.**

## ⭐⭐⭐ THE PARADIGM — develop on the test chromosome, decide on the ladder

The substrate is `docs/TESTING.md` §0a — **read it first**; it owns the description, the commands and
the two bars, and nothing is duplicated here. In one line: the ANCHORED TWIN BLOCK is one shape
repeated over five gene types × five abundance blocks, and the five types are the message layer's own
controls — `clean` (a pure-gDNA message SOURCE), `nasc` (a contaminated one), `cap` (the enrichment
cliff between source and destination), `capnasc` (both), `silent` (every object pure gDNA, so any
claimed RNA is a false positive).

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                          # FIRST — ~2 s, can this session run?

python scripts/design/policy_benchmark.py --panel test      # THE LOOP — 30 conditions, seconds
python scripts/design/policy_benchmark.py --panel ladder    # THE JUDGEMENT — 16 conditions, minutes
```

⛔ **THE OWNER'S RULING ON HOW TO USE IT: the test chromosome is where the policy is DEVELOPED, and a
result is only believed once the ladder confirms it** (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`
— a small substrate's objects mostly lack own evidence, so a message is nearly free there and a policy
looks better than it is). ⭐ **When a condition indicts a structure, ADD THAT STRUCTURE TO THE TEST
CHROMOSOME** — that is the loop, and the substrate is designed to grow one structure at a time
(`TESTING.md` §0a has the edit-then-rebuild sequence; the three hand-edited files are the only inputs).

⛔ **NEVER POOL THE TWO HALVES.** Unstranded rows are where a policy must WIN; stranded rows are where
it must do as little HARM as possible. They are judged against DIFFERENT bars.

## What the baseline says about where to aim

The full numbers are the dated snapshot beside this file. The three that set the direction:

1. **On the test chromosome the split is stark.** Relay beats silence by **0.01–0.32×** on unstranded
   × capture-ON — a large, real win exactly where the local solve is blind — and LOSES on stranded ×
   capture-ON by **2.1–4.96×** (and 58× on one near-zero-error row, where the ratio is not the thing
   to read). Off capture it is roughly neutral, degrading to **1.25–1.67×** harm as gDNA rises.
2. **On the ladder the same shape**: relay wins the blind rows and all four zero controls, and hurts
   the solvable set on 9/12 contaminated conditions.
3. ⛔ **The declared precision is NOT EARNED on 11/12 rows** (`solvability_audit.py`). That is the
   suspected mechanism: a message that claims more precision than it has will overwhelm good own
   evidence, which is precisely the measured stranded harm.

## The first two steps, in order

**① Measure the split — no new code.** Score silent vs relay **split by whether the destination had
its own composition evidence**. That is the unfinished half of `ISSUES: message-value-for-blind-slots`
and it converts "relay hurts 9/12" into "relay hurts THIS destination class by THIS much", which is
what says whether an emission gate alone recovers the harm. ⚠ Its no-solver companion is already
measured and recorded in that issue: an ORACLE neighbour message has POSITIVE skill only in a narrow
regime (mid-contamination and `g98` off capture) and is actively MISLEADING at low gDNA — so a policy
that speaks everywhere cannot win, however well engineered.

**② The precision ledger, per node type** — the owner's own question in `message_notes.md`, and the
right one: *when the strand model solves an exon, what precision does that exon actually get?* Three
numbers per node type: what `τ_λ` the strand channel EARNS (`region_init.build_region_init`), what the
relay DECLARES after `hop_logvar` damping (`messages/relay.py`), and their ratio. The owner's node
order — exon regions → `intron|exon` boundaries → intron regions — is right, and exons first is right
because that is where the harm concentrates.

**Then one mechanism at a time**, whichever the ledger indicts: the emission gate (speak only into
blind destinations) or the precision correction (declare what you earned). ⛔ Not both —
`CLAUDE.md`'s DERIVE → DESIGN → PLAN → PROTOTYPE → A/B, and a change that cannot be A/B'd alone
cannot be judged alone.

⛔ **What NOT to do: flip the default to silent.** It is tempting at 1.65× worst harm, but relay wins
all four zero controls and the blind rows, so a blanket flip trades those away. The destination gate
is the shape that keeps both.

## Standing cautions

* ⛔ **Read `HONEST_PRECISION.md` beside this file before proposing any mechanism** — it is the record
  of the 2026-08-27 campaign that was torn down, so a refuted experiment is not re-run.
* ⛔ `MessagePolicy` is byte-identical to `SilentPolicy` at rung 0 and confirmed so on all 46
  conditions. That identity is the foundation spec's gate — keep it as the falsification for every
  new rung (`TRAPS: an-ablation-that-never-ran`).
* ⚠ The certified-flux stream rides inside the relay (`RelaySwitches.certified_flux`, final solve
  only) and `config.rna_anchor` is live iff propagation is on AND the policy is relay — so a
  silent-vs-relay delta bundles the anchor. Separate it with
  `ladder_arm_ab.py --arm anchor_off` or `backbone_parity.py --arm-b no_certified_flux`.
* ⚠ Nothing measured before 2026-08-31 is comparable to current numbers (the drained-frame ruling,
  `DESIGN.md` §4.3).
