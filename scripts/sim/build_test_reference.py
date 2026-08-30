#!/usr/bin/env python3
"""BUILD THE METHOD-DEVELOPMENT TEST REFERENCE — one 100 kb chromosome that grows one transcript at a time.

⭐⭐⭐ **THIS IS THE BENCHMARK THE MESSAGE-PROPAGATION POLICY IS DEVELOPED AGAINST** (owner, 2026-08-19):
*"over the course of this method development, we will add transcripts to a 'test chromosome' … the test
reference chromosome, test transcript GTF, and test abundances comprise a critical method development
benchmark."* It starts at **ZERO transcripts** and every addition is a deliberate step: the policy must
solve the new structure AND everything already there.

**THE THREE PIECES, AND WHICH ONE IS EDITED BY HAND**

===========================  =========================================================================
`test_chr.gtf`               ⭐ **HAND-EDITED.** One `exon` line per exon. This is where a transcript
                             is added
`test_abundances.tsv`        ⭐ **HAND-EDITED.** `transcript_id / mrna_abundance / nrna_abundance` —
                             relative molecular sampling weights (owner: *"abundances … are relative
                             molecular abundance levels"*)
`test_chr.fa` + the index    ⛔ **DERIVED — never hand-edited.** This script writes them
===========================  =========================================================================

⛔ **WHY THE FASTA IS DERIVED RATHER THAN VERSIONED.** A spliced transcript needs a GT..AG at every
intron or the aligner and the simulator disagree with the annotation. Deriving the sequence FROM the GTF
makes that impossible to get wrong: the chromosome is regenerated from a fixed seed and every intron the
GTF declares gets its motif injected. A versioned FASTA would have to be re-edited by hand on every
addition, and the first missed edit is a silent, confusing failure.

⭐ **The nascent RNA needs no declaration.** `rigel index` creates a single-exon nascent ENTITY spanning
each multi-exon transcript (`index.create_nrna_transcripts`), which is exactly the owner's *"just add a
single-exon transcript spanning the multi-exon transcripts"*. Give an entity abundance by naming its
CONTRIBUTOR's `nrna_abundance`; the simulator pools it onto the entity.

Usage::

    python scripts/sim/build_test_reference.py                 # build into the default runs dir
    python scripts/sim/build_test_reference.py --out DIR       # elsewhere
    python scripts/sim/build_test_reference.py --self-test     # no I/O
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from rigel.sim.genome import MutableGenome  # noqa: E402
from rigel.transcript import Transcript  # noqa: E402
from rigel.types import Strand  # noqa: E402

#: ⭐ 1 Mb. It started at the owner's 100 kb — *"a scratchpad to add transcripts"* — and was raised the
#: same day (*"in hindsight we are going to need more than 100kb"*) to leave room for the loci this
#: benchmark will grow: overlapping genes, long introns, alternate TSS/TES, nested isoforms. ⚠ Still
#: tiny to simulate, and the gDNA depth is what sets the density, not the length.
GENOME_LENGTH = 1_000_000
REF_NAME = "test_chr"
#: the chromosome's random sequence is fixed by this seed, so the same GTF always gives the same FASTA
GENOME_SEED = 20260819
#: ⭐⭐ THE BLANK CHROMOSOME (owner design, 2026-08-29): a second contig with NO annotation in `test_chr.gtf`,
#: so `rigel index` knows nothing on it — while the SIMULATOR draws unannotated "shadow" transcripts from
#: `test_shadow.gtf` there. Its own seed; its own length. gDNA is drawn on it too (`gdna.genomic_refs`).
BLANK_REF_NAME = "test_blank"
BLANK_GENOME_LENGTH = 1_000_000
BLANK_GENOME_SEED = 20260829
REF_LENGTHS = {REF_NAME: GENOME_LENGTH, BLANK_REF_NAME: BLANK_GENOME_LENGTH}

HERE = Path(__file__).resolve().parent / "test_reference"
DEFAULT_GTF = HERE / "test_chr.gtf"
DEFAULT_SHADOW_GTF = HERE / "test_shadow.gtf"
DEFAULT_ABUNDANCES = HERE / "test_abundances.tsv"
DEFAULT_OUT = Path.home() / "Downloads" / "rigel_runs" / "test_reference"


def read_test_gtf(path: Path) -> list[Transcript]:
    """The hand-edited GTF, as transcripts. An empty file is legal and is where this benchmark starts."""
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist")
    body = [ln for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    if not body:
        return []
    return Transcript.read_gtf(str(path), parse_mode="strict")


def read_abundances(path: Path) -> dict[str, tuple[float, float]]:
    """``transcript_id -> (mrna_abundance, nrna_abundance)``. Header-only is legal (zero transcripts)."""
    rows: dict[str, tuple[float, float]] = {}
    for i, line in enumerate(path.read_text().splitlines()):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if i == 0 and fields[0] == "transcript_id":
            continue
        if len(fields) < 3:
            raise ValueError(f"{path}:{i + 1}: expected transcript_id / mrna / nrna, got {line!r}")
        rows[fields[0]] = (float(fields[1]), float(fields[2]))
    return rows


def intron_spans(transcripts: list[Transcript]) -> list[tuple[int, int, Strand]]:
    """Every intron the GTF declares, as ``(start, end, strand)`` — what needs a splice motif."""
    out: list[tuple[int, int, Strand]] = []
    for t in transcripts:
        exons = sorted(t.exons, key=lambda e: e.start)
        for left, right in zip(exons, exons[1:]):
            if right.start > left.end:
                out.append((left.end, right.start, t.strand))
    return out


def inject_motifs(genome: MutableGenome, introns: list[tuple[int, int, Strand]]) -> int:
    """GT..AG at every intron, reverse-complemented on the minus strand (CT..AC in genomic terms).

    ⛔ The motif is written in GENOMIC orientation, which is what an aligner reads. Getting this
    backwards on ``−`` is invisible in the FASTA and shows up only as unaligned spliced reads.
    """
    for start, end, strand in introns:
        donor, acceptor = ("GT", "AG") if strand != Strand.NEG else ("CT", "AC")
        genome.edit(start, donor)
        genome.edit(end - 2, acceptor)
    return len(introns)


def check_transcripts(transcripts: list[Transcript], abundances: dict, shadows: list[Transcript] | None = None) -> list[str]:
    """Everything that would make a transcript unusable, reported together rather than one per run.

    ``transcripts`` must live on ``test_chr`` (the annotated chromosome, the one the index is built from);
    ``shadows`` must live on ``test_blank`` (the BLANK chromosome — unannotated by construction); ids are
    unique across BOTH files (a shadow the index would know is not a shadow); every transcript of either
    file has an abundance row and every row names a transcript of one of them."""
    problems: list[str] = []
    seen: dict[str, Transcript] = {}
    shadows = shadows or []
    for t, expected_ref, length in [(t, REF_NAME, GENOME_LENGTH) for t in transcripts] + [
        (t, BLANK_REF_NAME, BLANK_GENOME_LENGTH) for t in shadows
    ]:
        if t.ref != expected_ref:
            problems.append(f"{t.t_id}: reference {t.ref!r}, expected {expected_ref!r}")
        if t.t_id in seen:
            problems.append(f"{t.t_id}: declared twice")
        seen[t.t_id] = t
        exons = sorted(t.exons, key=lambda e: e.start)
        for e in exons:
            if e.start < 0 or e.end > length:
                problems.append(f"{t.t_id}: exon [{e.start}, {e.end}) outside [0, {length})")
            if e.end <= e.start:
                problems.append(f"{t.t_id}: empty exon [{e.start}, {e.end})")
        for left, right in zip(exons, exons[1:]):
            if right.start < left.end:
                problems.append(f"{t.t_id}: exons overlap at {left.end}/{right.start}")
            elif right.start - left.end < 4:
                # ⛔ a GT..AG needs 4 bases; a shorter gap cannot carry a motif
                problems.append(f"{t.t_id}: intron [{left.end}, {right.start}) is under 4 bp")
        if t.t_id not in abundances:
            problems.append(f"{t.t_id}: no row in the abundances file")
    for t_id in abundances:
        if t_id not in seen:
            problems.append(f"{t_id}: has an abundance but is not in either GTF")
    return problems


def _write_two_record_fasta(genomes: list[MutableGenome], path: Path) -> Path:
    """One FASTA, two records, 80-column wrap, + ``samtools faidx``. ⚠ The file keeps the name
    ``test_chr.fa`` so every panel config's ``genome:`` path stays valid; it now carries BOTH contigs."""
    import pysam

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for g in genomes:
            f.write(f">{g.name}\n")
            seq = g.seq
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")
    for stale in (path.with_suffix(path.suffix + ".fai"),):
        if stale.exists():
            stale.unlink()
    pysam.faidx(str(path))
    return path


def build(gtf: Path, abundances_path: Path, out: Path, shadow_gtf: Path | None = None) -> dict:
    """Write the derived FASTA (both contigs) and report what the benchmark currently holds."""
    transcripts = read_test_gtf(gtf)
    shadows = read_test_gtf(shadow_gtf) if shadow_gtf is not None and shadow_gtf.is_file() else []
    abundances = read_abundances(abundances_path)
    problems = check_transcripts(transcripts, abundances, shadows)
    if problems:
        raise ValueError(
            f"{len(problems)} problem(s) in the test reference:\n  " + "\n  ".join(problems)
        )
    genome = MutableGenome(GENOME_LENGTH, seed=GENOME_SEED, name=REF_NAME)
    introns = intron_spans(transcripts)
    inject_motifs(genome, introns)
    blank = MutableGenome(BLANK_GENOME_LENGTH, seed=BLANK_GENOME_SEED, name=BLANK_REF_NAME)
    shadow_introns = intron_spans(shadows)
    inject_motifs(blank, shadow_introns)
    out.mkdir(parents=True, exist_ok=True)
    fasta = _write_two_record_fasta([genome, blank], out / f"{REF_NAME}.fa")
    # the GTFs and abundances the tools read are COPIES of the hand-edited ones, so the built reference is
    # self-contained. ⛔ `test_shadow.gtf` is copied for the SIMULATOR only — the index never reads it.
    (out / "test_chr.gtf").write_text(gtf.read_text())
    (out / "test_abundances.tsv").write_text(abundances_path.read_text())
    if shadow_gtf is not None and shadow_gtf.is_file():
        (out / "test_shadow.gtf").write_text(shadow_gtf.read_text())
    multi = sum(1 for t in transcripts if len(t.exons) > 1)
    return {"fasta": fasta, "n_transcripts": len(transcripts), "n_multi_exon": multi,
            "n_introns": len(introns), "length": GENOME_LENGTH,
            "n_nascent_entities_expected": multi,
            "n_shadow_transcripts": len(shadows), "n_shadow_introns": len(shadow_introns),
            "blank_ref": BLANK_REF_NAME, "blank_length": BLANK_GENOME_LENGTH}


def self_test() -> int:
    """⛔ Every check perturbed, no I/O against the real reference."""
    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"   ⛔ {name}")

    from rigel.types import Interval

    def tx(t_id, exons, strand=Strand.POS, ref=REF_NAME):
        t = Transcript(ref=ref, strand=strand, exons=[Interval(a, b) for a, b in exons], t_id=t_id,
                       g_id="g_" + t_id)
        t.length = t.compute_length()
        return t

    # ── the genome is reproducible, and the seed is what makes it so
    a = MutableGenome(2000, seed=GENOME_SEED, name=REF_NAME).seq
    b = MutableGenome(2000, seed=GENOME_SEED, name=REF_NAME).seq
    check("the chromosome is reproducible from its seed", a == b)
    check("a different seed gives a different chromosome",
          MutableGenome(2000, seed=GENOME_SEED + 1, name=REF_NAME).seq != a)

    # ── motifs: written, in genomic orientation, and only where the GTF says
    g = MutableGenome(2000, seed=1, name=REF_NAME)
    n = inject_motifs(g, intron_spans([tx("T", [(100, 200), (400, 500)])]))
    check("one intron gives one motif pair", n == 1)
    check("a + intron gets GT..AG", g[200:202] == "GT" and g[398:400] == "AG")
    gm = MutableGenome(2000, seed=1, name=REF_NAME)
    inject_motifs(gm, intron_spans([tx("T", [(100, 200), (400, 500)], strand=Strand.NEG)]))
    check("a − intron gets CT..AC in genomic orientation",
          gm[200:202] == "CT" and gm[398:400] == "AC")
    check("a single-exon transcript declares no intron", intron_spans([tx("S", [(10, 900)])]) == [])

    # ── the checks fire on what would break a build, and pass on what would not
    good = [tx("T1", [(1000, 2000), (5000, 6000)])]
    ab = {"T1": (100.0, 25.0)}
    check("a well-formed transcript with an abundance passes", check_transcripts(good, ab) == [])
    check("an exon past the end of the chromosome is caught",
          any("outside" in p for p in check_transcripts([tx("X", [(0, GENOME_LENGTH + 1)])], {"X": (1.0, 0.0)})))
    check("overlapping exons are caught",
          any("overlap" in p for p in check_transcripts([tx("X", [(100, 500), (400, 900)])], {"X": (1.0, 0.0)})))
    check("an intron too short to carry a motif is caught",
          any("under 4 bp" in p for p in check_transcripts([tx("X", [(100, 200), (202, 400)])], {"X": (1.0, 0.0)})))
    check("a transcript with no abundance row is caught",
          any("no row" in p for p in check_transcripts(good, {})))
    check("an abundance for a transcript that is not in either GTF is caught",
          any("not in either GTF" in p for p in check_transcripts(good, {**ab, "GHOST": (1.0, 0.0)})))
    check("a duplicate transcript id is caught",
          any("twice" in p for p in check_transcripts(good + good, ab)))
    check("the wrong reference name is caught",
          any("expected" in p for p in check_transcripts([tx("X", [(10, 90)], ref="chr1")], {"X": (1.0, 0.0)})))

    # ── the BLANK chromosome and its shadows (2026-08-29)
    sh = [tx("shadow_A", [(1000, 2000), (5000, 6000)], ref=BLANK_REF_NAME)]
    both = {"T1": (100.0, 25.0), "shadow_A": (5.0, 0.0)}
    check("a shadow on the blank chromosome with an abundance row passes", check_transcripts(good, both, sh) == [])
    check("a shadow placed on the ANNOTATED chromosome is caught",
          any("expected 'test_blank'" in p for p in check_transcripts(good, both, [tx("shadow_A", [(1000, 2000), (5000, 6000)])])))
    check("an annotated transcript placed on the blank chromosome is caught",
          any("expected 'test_chr'" in p for p in check_transcripts([tx("T1", [(1000, 2000), (5000, 6000)], ref=BLANK_REF_NAME)], both, sh)))
    check("a shadow id the annotation already declares is caught (it would not be a shadow)",
          any("twice" in p for p in check_transcripts(good, both, [tx("T1", [(1000, 2000), (5000, 6000)], ref=BLANK_REF_NAME)])))
    check("a shadow with no abundance row is caught", any("no row" in p for p in check_transcripts(good, ab, sh)))
    gb = MutableGenome(8000, seed=BLANK_GENOME_SEED, name=BLANK_REF_NAME)
    inject_motifs(gb, intron_spans(sh))
    check("a shadow intron gets its motif on the BLANK genome", gb[2000:2002] == "GT" and gb[4998:5000] == "AG")
    check("the two contigs are distinct sequences",
          MutableGenome(2000, seed=GENOME_SEED, name=REF_NAME).seq != MutableGenome(2000, seed=BLANK_GENOME_SEED, name=BLANK_REF_NAME).seq)

    # ── zero transcripts is the STARTING state and must be legal
    check("zero transcripts is legal", check_transcripts([], {}) == [])
    check("zero transcripts declares no intron", intron_spans([]) == [])

    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    ap.add_argument("--abundances", type=Path, default=DEFAULT_ABUNDANCES)
    ap.add_argument("--shadow-gtf", type=Path, default=DEFAULT_SHADOW_GTF,
                    help="unannotated SHADOW transcripts on the blank chromosome (simulator-only)")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    info = build(args.gtf, args.abundances, args.out, args.shadow_gtf)
    print(f"\n⭐ THE TEST REFERENCE — {REF_NAME}, {info['length']:,} bp")
    print(f"   fasta        {info['fasta']}")
    print(f"   gtf          {args.gtf}")
    print(f"   abundances   {args.abundances}")
    print(f"   transcripts  {info['n_transcripts']}  ({info['n_multi_exon']} multi-exon)")
    print(f"   introns      {info['n_introns']} (each with a GT..AG injected)")
    print(f"   ⭐ nascent entities `rigel index` will create: {info['n_nascent_entities_expected']}"
          " — one single-exon transcript spanning each multi-exon one")
    print(f"   ⭐⭐ BLANK chromosome {info['blank_ref']}, {info['blank_length']:,} bp — NO annotation; "
          f"{info['n_shadow_transcripts']} SHADOW transcripts ({info['n_shadow_introns']} shadow introns with motifs) "
          f"simulated from {args.shadow_gtf.name}, never given to the index")
    if info["n_transcripts"] == 0:
        print("\n   ⭐ ZERO TRANSCRIPTS — the benchmark's starting state. Add one to "
              f"{args.gtf.name} and a row to {args.abundances.name}, then re-run this.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
