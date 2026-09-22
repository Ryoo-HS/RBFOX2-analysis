"""BAM -> per-strand, per-position coverage (splice-aware)."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pysam

FLAG_UNMAPPED = 0x4
FLAG_SECONDARY = 0x100
FLAG_QCFAIL = 0x200
FLAG_DUP = 0x400
FLAG_SUPPLEMENTARY = 0x800


@dataclass
class CoverageResult:
    chroms: dict
    chrom_lengths: dict
    total_reads: int
    total_aligned_bp: int
    strands: tuple
    meta: dict = field(default_factory=dict)


def _feature_strand(is_reverse, library):
    if library == "unstranded":
        return "."
    genomic = "-" if is_reverse else "+"
    if library == "reverse":
        return "+" if genomic == "-" else "-"
    return genomic


def build_coverage(bam_path, *, library="forward", min_mapq=0, keep_dup=False,
                   min_read_length=0, chroms=None):
    """Whole-chromosome coverage arrays (used by tests / small regions).

    For large / whole-genome BAMs the pipeline uses the memory-efficient
    ``chrom_islands`` instead, which allocates only island-sized arrays.
    """
    if library not in ("forward", "reverse", "unstranded"):
        raise ValueError(f"unknown library: {library}")
    strands = (".",) if library == "unstranded" else ("+", "-")
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref_lengths = dict(zip(bam.references, bam.lengths))
    wanted = chroms if chroms is not None else list(bam.references)
    drop_flags = FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_QCFAIL
    if not keep_dup:
        drop_flags |= FLAG_DUP
    cov = {}
    total_reads = 0
    total_aligned_bp = 0
    for chrom in wanted:
        length = ref_lengths[chrom]
        arrays = {s: np.zeros(length, dtype=np.int32) for s in strands}
        for read in bam.fetch(chrom):
            if read.flag & drop_flags:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if min_read_length and (read.query_length or 0) < min_read_length:
                continue
            st = _feature_strand(read.is_reverse, library)
            arr = arrays[st]
            for b_start, b_end in read.get_blocks():
                arr[b_start:b_end] += 1
                total_aligned_bp += (b_end - b_start)
            total_reads += 1
        if any(a.any() for a in arrays.values()):
            cov[chrom] = arrays
    bam.close()
    return CoverageResult(chroms=cov, chrom_lengths=ref_lengths, total_reads=total_reads,
                          total_aligned_bp=total_aligned_bp, strands=strands,
                          meta={"library": library})


def _drop_flags(keep_dup: bool) -> int:
    d = FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_QCFAIL
    if not keep_dup:
        d |= FLAG_DUP
    return d


def _passes(read, drop_flags, min_mapq, min_read_length):
    if read.flag & drop_flags:
        return False
    if read.mapping_quality < min_mapq:
        return False
    if min_read_length and (read.query_length or 0) < min_read_length:
        return False
    return True


def count_library(bam_path, *, min_mapq=0, keep_dup=False, min_read_length=0, chroms=None):
    """One light pass (no arrays): total passing reads and aligned bp."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    drop = _drop_flags(keep_dup)
    wanted = chroms if chroms is not None else list(bam.references)
    n = 0
    abp = 0
    for chrom in wanted:
        for r in bam.fetch(chrom):
            if not _passes(r, drop, min_mapq, min_read_length):
                continue
            n += 1
            abp += (r.query_alignment_length or 0)
    bam.close()
    return n, abp


def _grouped_blocks(bam_path, chrom, *, library="forward", min_mapq=0, keep_dup=False,
                    min_read_length=0, coverage_gap=200):
    """Group raw read blocks into (strand, start, end, block_list) islands.
    Shared by chrom_islands (builds arrays eagerly) and chrom_island_blocks
    (keeps the raw blocks, e.g. for a permutation null model)."""
    if library not in ("forward", "reverse", "unstranded"):
        raise ValueError(f"unknown library: {library}")
    bam = pysam.AlignmentFile(bam_path, "rb")
    drop = _drop_flags(keep_dup)
    strands = (".",) if library == "unstranded" else ("+", "-")
    blocks = {s: [] for s in strands}
    n = 0
    abp = 0
    for r in bam.fetch(chrom):
        if not _passes(r, drop, min_mapq, min_read_length):
            continue
        st = _feature_strand(r.is_reverse, library)
        for b in r.get_blocks():
            blocks[st].append(b)
            abp += b[1] - b[0]
        n += 1
    bam.close()
    islands = []
    for st, bl in blocks.items():
        if not bl:
            continue
        bl.sort()
        cs, ce, group = bl[0][0], bl[0][1], [bl[0]]
        groups = []
        for bs, be in bl[1:]:
            if bs - ce <= coverage_gap:
                ce = max(ce, be)
                group.append((bs, be))
            else:
                groups.append((cs, ce, group))
                cs, ce, group = bs, be, [(bs, be)]
        groups.append((cs, ce, group))
        for s, e, gb in groups:
            islands.append((st, s, e, gb))
    return islands, n, abp


def chrom_island_blocks(bam_path, chrom, *, library="forward", min_mapq=0, keep_dup=False,
                        min_read_length=0, coverage_gap=200):
    """Like chrom_islands, but returns the raw (start, end) read blocks per
    island instead of a built array -- e.g. for a read-block permutation null
    model (see calibrate.auto_min_prominence).
    Returns (islands, n_reads, aligned_bp) where each island is
    (strand, start, end, block_list)."""
    return _grouped_blocks(bam_path, chrom, library=library, min_mapq=min_mapq,
                           keep_dup=keep_dup, min_read_length=min_read_length,
                           coverage_gap=coverage_gap)


def blocks_to_array(start, end, blocks):
    """Per-bp depth array for one island, built from its raw read blocks."""
    arr = np.zeros(end - start, dtype=np.int32)
    for bs, be in blocks:
        arr[bs - start:be - start] += 1
    return arr


def chrom_islands(bam_path, chrom, *, library="forward", min_mapq=0, keep_dup=False,
                  min_read_length=0, coverage_gap=200):
    """Memory-efficient per-chromosome coverage as (strand, start, end, array)
    for each covered island only -- never a whole-chromosome array.
    ``coverage_gap`` merges raw read blocks within this many bp into one array
    chunk -- a memory/chunking concern, independent of the stage-1 region gap
    (see ``island.find_islands``'s ``region_gap``).
    Returns (islands, n_reads, aligned_bp)."""
    grouped, n, abp = _grouped_blocks(
        bam_path, chrom, library=library, min_mapq=min_mapq, keep_dup=keep_dup,
        min_read_length=min_read_length, coverage_gap=coverage_gap)
    islands = [(st, s, e, blocks_to_array(s, e, gb)) for st, s, e, gb in grouped]
    return islands, n, abp


def stream_bedgraph(bam_path, out_base, *, library="forward", min_mapq=0, keep_dup=False,
                    min_read_length=0, coverage_gap=200, factor=1.0, chroms=None):
    """Write per-strand bedGraph tracks by streaming per-chromosome islands (low memory)."""
    from .normalize import to_rpm
    from .output import write_bedgraph
    bam = pysam.AlignmentFile(bam_path, "rb")
    wanted = chroms if chroms is not None else list(bam.references)
    bam.close()
    strands = (".",) if library == "unstranded" else ("+", "-")
    tag = {"+": "plus", "-": "minus", ".": "unstranded"}
    handles = {st: open(f"{out_base}.{tag[st]}.bdg", "w") for st in strands}
    for chrom in wanted:
        isls, _, _ = chrom_islands(bam_path, chrom, library=library, min_mapq=min_mapq,
                                   keep_dup=keep_dup, min_read_length=min_read_length,
                                   coverage_gap=coverage_gap)
        for st, s, e, arr in isls:
            write_bedgraph(chrom, to_rpm(arr, factor), handles[st], strand=st, offset=s)
    paths = [h.name for h in handles.values()]
    for h in handles.values():
        h.close()
    return paths
