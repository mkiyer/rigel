# THE ISOFORM BLOCK — the substrate for holes ① and ② of the `transfer` policy (proposal, 2026-09-02)

    ⚠ A DEV DOC — a proposal for the owner, with the ladder evidence that shaped it. Nothing here
    is settled. When the block is authored, the design moves into the GTF header (the paradigm:
    the GTF IS the benchmark) and this file shrinks to the rung's working notes.

## 1. What the ladder says the holes ARE (census, no solver, `g50 ss.50 OFF` chain)

The audit's "12,811 exon|exon slots" split by the splice graph's flags, and the split changes the
design: **the dominant class is the internal TERMINUS, not the alternative splice site.**

| exon\|exon boundary class | slots | crossing mass |
|---|---|---|
| **TERMINUS, same strand** (a TSS/TES of one isoform inside another's exon) | **7,604** | **950,516** |
| splice site, same strand (an alternative 5'/3' site; LICENSED by the rung-2 predicate) | 4,838 | 457,797 |
| terminus, strand change | 246 | 26,529 |
| sj+term | 123 | 8,442 |

Nearly every terminus boundary carries ONE flag (TSS+ 2,089 · TES+ 1,937 · TSS− 1,850 · TES− 1,726),
so the terminating transcripts' DIRECTION is unambiguous at ~99 % of them.

The exons, by their two faces (24,018 exons on the chain):

| exon class | exons | mass | reached today? |
|---|---|---|---|
| both faces licensed intron (the twin-block shape) | 4,062 | 76,619 | rung 2 |
| both faces TERMINUS exon\|exon (walled) | **3,917** | **345,965** | ⛔ no |
| terminus exon\|exon + licensed intron | 2,813 | 271,458 | rung 2 (one witness) |
| licensed exon\|exon + terminus exon\|exon | 2,564 | 262,823 | ⛔ no |
| licensed exon\|exon + licensed intron | 2,117 | 58,347 | rung 2 (one witness) |
| both faces licensed exon\|exon | 1,804 | 190,035 | ⛔ no — multi-hop only |
| edge + intron / edge + edge (the mono shape) | 974 / 396 | 132,131 / 68,771 | rung 3 |

Reach by a chain of LICENSED faces from the nearest intron factory, walking through exons: depth 1 =
11,418 exons (658,790), depth 2 = 1,582 (136,005), depth 3+ = 759 (83,388), and **UNREACHABLE by any
licensed chain: 10,259 exons, 1,051,310 mass** — walled by terminus faces (11,630 of their faces are
`exon|exon[term]`, 1,478 are `exon|intron[term]`). So hole ① is two holes: a small licensed multi-hop
one (the sj class) and a large terminus one; hole ② is the same terminus mechanism at an intron flank.
The strand-change classes are small (246 + 888 + 15 walled exons) and stay a later rung.

## 2. THE DIRECTED LICENCE — certified (stage 0)

Measured on ladder truth and re-measured per type on the built block: the OUTSIDE flank of a terminus
boundary shares the crossing's composition to within counting noise, the INSIDE flank does not. The
table, the derivation and the prototype results MOVED to `COMPOSITION_TRANSFER_STAGE01.md` (RUNG 4).

## 3. THE BLOCK — AUTHORED, then TRIMMED to `altstart` (owner rulings 2026-09-02)

⚠ The four-structure build below ran once and was then trimmed to ONE structure (`altstart`, probed and
unprobed) per the owner's incremental ruling; the others are queued in `MESSAGE_RUNGS.md`. The layout,
strand rule and verification stand for each structure as it is re-added.

Decided and built, in `scripts/sim/test_reference/test_chr.yaml` (the YAML header is the design
record): **8 types × 5 blocks = 40 genes, 80 transcripts** — four structures, each probed and unprobed
(`altss` · `altstart` · `nest` · `instart` · `capaltss` · `capaltstart` · `capnest` · `capinstart`),
host `U` by the mature ladder, `T` by the anti-ladder, no nascent. Depth held by raising
`n_total_fragments` 200 k → 480 k in all seven panel configs. Genome 2.9 Mb, blank contig kept at 1.5 Mb.

**Strands (owner refinement):** not "all +" and not "one type on −" — every gene on the chromosome
gets an explicit strand by `(type index + block) mod 2`, so every type sits on both strands across its
five blocks and the chromosome is balanced (42 + / 43 −; twin and mono blocks re-stranded too). On a −
gene the isoform geometry is mirrored within the 17 kb span so each type keeps its biological meaning
(a − `altstart` still carries an alternative TSS, flagged TSS− at the mirrored position). Both-stranded
(overlapping, opposite-strand) loci are a later step and remain absent. Verified on the built index: the
intended boundary classes appear on both strands (`altss`: ACC/DON exon|exon pairs; `altstart`: one
TSS± exon|exon; `nest`: TES±+TSS± walls; `instart`: TSS± intron|exon + ACC± intron|exon).

Offsets from the gene start `s` on a + gene (U's exons `[s, s+1000) [s+8000, s+9000) [s+16000, s+17000)`):

| type | T's exons | creates | indicts |
|---|---|---|---|
| `altss` | `[s, s+1000) [s+7500, s+9500) [s+16000, s+17000)` | two licensed exon\|exon (ACC of U at s+8000, DON of U at s+9000); the core `[s+8000, s+9000)` has ONLY exon\|exon faces | licensed multi-hop (1,804 + 2,117 ladder exons) |
| `altstart` | `[s+8500, s+9000) [s+16000, s+17000)` | one exon\|exon TERMINUS (TSS at s+8500); outside `[s+8000, s+8500)` U-only, inside U+T | terminus + licensed intron (2,813) |
| `nest` | `[s+8250, s+8750)` | TSS and TES inside U's exon: the middle walled by termini both sides | the walled class (3,917 exons, 346 k) |
| `instart` | `[s+4000, s+5000) [s+8000, s+9000) [s+16000, s+17000)` | intron\|exon TERMINUS (TSS at s+4000) + intron\|exon sj (DON at s+5000) | hole ② (2,979 slots) |

The rebuild recipe lives in `docs/TESTING.md` §0a now (moved). The old derived set is at
`~/Downloads/rigel_runs/test_reference_STALE_45tx_2026-09-02/`.

## 4. State (end of 2026-09-02)

Substrate built and certified on all seven panels; the rung-4 prototype ran on the benign panel
(record in the thread doc); the ladder run is the pending decision. The stage-0 instruments live in the
session scratchpad (`exon_exon_census.py`, `terminus_pair_gap.py`, `directed_reach_census.py`,
`rung4_proto.py` under `/private/tmp/claude-503/-Users-mkiyer-proj-rigel/<session>/scratchpad/`) —
promotion into `scripts/design/` with `--self-test` is owed if the rung proceeds.
