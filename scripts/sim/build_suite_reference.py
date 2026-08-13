"""Carve a benchmark-suite reference out of the real human genome, optionally appending a synthetic one.

    TODO item 2 · objective 2

⭐ **WHY A REAL BACKBONE.** The suite this replaces was a generated mini-genome, and it could not judge
what it was used to judge: zero fragment-length variance, Poisson by
construction, and a fine region set row-for-row identical to its merged one. Real human genes give
calibration a real fragment-length distribution, a real strand model, and — the point of the v8 partition
— alternative TSS/TES that fall strictly inside exons. The owner's plan is one chromosome as the training
substrate, with a synthetic stress chromosome piggybacked on top.

⚠ **The ERCC controls are kept deliberately, and not as filler.**: a
single-reference synthetic index hid a reference-id-space mismatch that silently dropped **476,719 of
476,732** real fragments inside `deposit()` while every golden test passed. 92 tiny references cost
~83 kb and make the ref-id space non-trivial, which is the configuration that would have caught it.

    python scripts/sim/build_suite_reference.py --fasta G.fa.bgz --gtf A.gtf -o OUT \\
        --refs chr22 --ercc [--append-fasta S.fa --append-gtf S.gtf]

Writes `OUT/genome.fa` (+ `.fai`) and `OUT/genes.gtf`, and prints the census of what it kept.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pysam

#: FASTA line width. A formatting choice with no effect on any downstream sequence.
_FASTA_BOUNDARY_WIDTH = 60


def selected_references(fai_path: Path, refs: list[str], ercc: bool) -> list[str]:
    """Resolve the requested reference names against the source `.fai`, preserving its order.

    Source order matters: reference ids are assigned by it, and keeping it means the suite's ref-id
    space is a subsequence of the production one rather than a reshuffle.
    """
    available = [line.split("\t")[0] for line in fai_path.read_text().splitlines() if line]
    wanted = set(refs)
    if ercc:
        wanted |= {name for name in available if name.startswith("ERCC-")}

    missing = wanted - set(available)
    if missing:
        raise SystemExit(f"references not in {fai_path}: {sorted(missing)}")
    return [name for name in available if name in wanted]


def write_fasta(source: Path, names: list[str], out_path: Path, append: Path | None) -> None:
    fasta = pysam.FastaFile(str(source))
    with open(out_path, "w") as fh:
        for name in names:
            sequence = fasta.fetch(name)
            fh.write(f">{name}\n")
            for i in range(0, len(sequence), _FASTA_BOUNDARY_WIDTH):
                fh.write(sequence[i : i + _FASTA_BOUNDARY_WIDTH] + "\n")
        if append is not None:
            fh.write(append.read_text())
    pysam.faidx(str(out_path))


def write_gtf(source: Path, names: list[str], out_path: Path, append: Path | None) -> int:
    keep = set(names)
    written = 0
    with open(source) as src, open(out_path, "w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            if line.split("\t", 1)[0] in keep:
                dst.write(line)
                written += 1
        if append is not None:
            for line in append.read_text().splitlines(keepends=True):
                if not line.startswith("#"):
                    written += 1
                dst.write(line)
    return written


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fasta", type=Path, required=True, help="Source genome FASTA (faidx'd)")
    ap.add_argument("--gtf", type=Path, required=True, help="Source annotation GTF")
    ap.add_argument("-o", "--outdir", type=Path, required=True)
    ap.add_argument("--refs", nargs="*", default=[], help="Reference names to keep, e.g. chr22")
    ap.add_argument("--ercc", action="store_true", help="Also keep every ERCC-* spike-in reference")
    ap.add_argument("--append-fasta", type=Path, default=None, help="Synthetic FASTA to append")
    ap.add_argument("--append-gtf", type=Path, default=None, help="Synthetic GTF to append")
    args = ap.parse_args()

    if (args.append_fasta is None) != (args.append_gtf is None):
        raise SystemExit("--append-fasta and --append-gtf must be given together")

    args.outdir.mkdir(parents=True, exist_ok=True)
    names = selected_references(Path(f"{args.fasta}.fai"), args.refs, args.ercc)

    out_fasta = args.outdir / "genome.fa"
    out_gtf = args.outdir / "genes.gtf"
    write_fasta(args.fasta, names, out_fasta, args.append_fasta)
    n_gtf_lines = write_gtf(args.gtf, names, out_gtf, args.append_gtf)

    fai = [line.split("\t") for line in Path(f"{out_fasta}.fai").read_text().splitlines() if line]
    total_bp = sum(int(row[1]) for row in fai)
    print(f"references  {len(fai):,}  ({len(names):,} carved" + (
        f" + {len(fai) - len(names):,} appended)" if len(fai) > len(names) else ")"))
    print(f"total bp    {total_bp:,}")
    print(f"gtf lines   {n_gtf_lines:,}")
    print(f"written     {out_fasta}\n            {out_gtf}")


if __name__ == "__main__":
    main()
