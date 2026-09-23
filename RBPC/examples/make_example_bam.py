#!/usr/bin/env python3
"""Synthesize a small, deterministic BAM for demoing / testing the peak caller.

Layout on a single 20 kb contig ``chr_test`` (+ strand):
  * a broad, sparsely-covered background region (~2,000-8,000);
  * TWO adjacent MOUNTAINS near 5,000 and 5,300 -- read starts concentrated
    (Gaussian) so each rises and falls; a prominence-based caller must return
    them as TWO peaks, not one;
  * a weak bump near 12,000 (low prominence -> ignored);
  * sparse genome-wide noise;
  * one spliced read across an intron to exercise splice-aware coverage.
No external reference is needed -- reads carry a dummy sequence.
"""
from __future__ import annotations

import argparse
import random

import pysam

CHROM = "chr_test"
CHROM_LEN = 20_000
READ_LEN = 30


def _add(reads, pos, n, strand="+", cigar=None):
    for _ in range(n):
        reads.append((pos, strand, cigar))


def _mountain(reads, center, n, sd):
    for _ in range(n):
        s = int(random.gauss(center, sd))
        _add(reads, max(0, min(CHROM_LEN - READ_LEN, s)), 1)


def build(path: str, seed: int = 0) -> None:
    random.seed(seed)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}

    reads: list[tuple[int, str, list | None]] = []
    for pos in range(2000, 8000, 40):          # sparse background
        _add(reads, pos, 2)
    _mountain(reads, 5000, 400, 45)            # mountain 1
    _mountain(reads, 5300, 300, 40)            # mountain 2 (adjacent -> must split)
    for pos in range(11800, 12200, 40):        # weak bump (low prominence)
        _add(reads, pos, 3)
    for _ in range(300):                       # noise
        _add(reads, random.randint(0, CHROM_LEN - READ_LEN), 1)
    _add(reads, 9000, 5, cigar=[(0, 10), (3, 2000), (0, 20)])   # spliced

    reads.sort(key=lambda r: r[0])
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for i, (pos, strand, cigar) in enumerate(reads):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 60
            a.flag = 16 if strand == "-" else 0
            if cigar is None:
                a.cigartuples = [(0, READ_LEN)]
                a.query_sequence = "A" * READ_LEN
            else:
                a.cigartuples = cigar
                a.query_sequence = "A" * sum(l for op, l in cigar if op in (0, 1, 4))
            a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
            bam.write(a)

    pysam.index(path)
    print(f"wrote {path} ({len(reads)} reads) + index")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="example.bam")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    build(args.output, args.seed)
