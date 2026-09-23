#!/usr/bin/env python3
"""Diagnose --region-min-depth: compare a background-depth statistic
(median / upper-trimmed mean of *covered* bases) against the existing
calibrate.auto_region_min_depth (75th percentile of per-island MAX depth),
per BAM, per strand.

Read-only diagnostic -- does not touch any pipeline code. `background_stats`
below is the single definition of "bg" this diagnostic uses, meant to be
reused as-is by the future --region-min-fold-bg implementation (stage 2)
so the two paths can never compute two different things called "bg".

Usage:
    python scripts/diag_bg.py --mode subset   # chr16+chr6 subsets in ~/test_subset
    python scripts/diag_bg.py --mode full     # whole-genome, original 4 BAMs
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pysam

from rbpc import calibrate
from rbpc.calibrate import background_stats  # noqa: F401 -- canonical "bg" impl,
                                                    # shared with auto_summit_floor
                                                    # (--summit-margins, stage 2)
from rbpc.coverage import _drop_flags, _feature_strand, chrom_islands
from rbpc.normalize import scaling_factor

TRIM_PCT = 2.0  # exclude the top 2% of covered-base depths before averaging

# --- preset BAM locations -----------------------------------------------
# library="forward" for all four: confirmed on ESRP1/RBFOX2/QKI via motif
# orientation + gene-strand concordance (project memory, 2026-08-02). YBX1's
# forward assumption is NOT independently verified the same way -- flagged
# in the report. This does not bias region-min-depth itself (pooled across
# strands, see auto_region_min_depth), only which physical reads land in the
# "+" vs "-" bucket of the per-strand rows below.
SUBSET_BAMS = {
    "ESRP1": os.path.expanduser("~/test_subset/ESRP1_chr16_chr6.bam"),
    "RBFOX2": os.path.expanduser("~/test_subset/RBFOX2_chr16_chr6.bam"),
    "YBX1": os.path.expanduser("~/test_subset/YBX1_chr16_chr6.bam"),
    "QKI": os.path.expanduser("~/test_subset/QKI_chr16_chr6.bam"),
}
SUBSET_CHROMS = ["chr16", "chr6"]
# Subset BAMs only contain reads on chr16/chr6, but their headers still list
# all 22 autosomes. auto_region_min_depth's default sample_n=4 would sample
# the 4 *largest* autosomes by header length (chr1-4), which are empty here
# and would silently fall back to 3.0. Passing sample_n=22 (all autosomes)
# guarantees chr6 (rank 6 by size) and chr16 (rank 16) are among the chroms
# it scans; the other 20 empty chroms just contribute nothing. This is a
# subset-only workaround -- the "full" mode below uses the real default.
SUBSET_AUTO_SAMPLE_N = 22

FULL_BAMS = {
    "ESRP1": os.path.expanduser(
        "~/ESRP1/ESRP1_hg19_hairpin_miscRNA_rRNA_tRNA_snoRNA_snRNA.sorted.bam"),
    "RBFOX2": os.path.expanduser(
        "~/RBFOX2/RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA_hg19.sorted.bam"),
    "YBX1": os.path.expanduser(
        "~/YBX1/YBX1_hg19_hairpin_miscRNA_rRNA_tRNA_snoRNA_snRNA_hg19.bam"),
    "QKI": os.path.expanduser(
        "~/QKI/QKI_hg19_hairpin_miscRNA_rRNA_tRNA_snRNA_snoRNA.sorted.bam"),
}

LIBRARY = {"ESRP1": "forward", "RBFOX2": "forward", "YBX1": "forward", "QKI": "forward"}


def covered_depths(bam_path, chroms, strand, *, library, min_mapq=0, keep_dup=False,
                    min_read_length=0, coverage_gap=200) -> np.ndarray:
    """Concatenated raw per-bp depth values (depth >= 1 only) for `strand`,
    pooled across `chroms`, built from the pipeline's own chrom_islands (same
    function pipeline.run() uses -- no reimplementation)."""
    parts = []
    for c in chroms:
        isls, _n, _abp = chrom_islands(
            bam_path, c, library=library, min_mapq=min_mapq, keep_dup=keep_dup,
            min_read_length=min_read_length, coverage_gap=coverage_gap)
        for st, _s, _e, arr in isls:
            if st == strand:
                covered = arr[arr > 0]
                if covered.size:
                    parts.append(covered)
    return np.concatenate(parts) if parts else np.array([], dtype=np.int64)


def count_reads(bam_path, chroms, *, library=None, strand=None, min_mapq=0, keep_dup=False,
                min_read_length=0) -> int:
    """Total passing reads over `chroms`. If `strand` is given, restrict to
    reads assigned that strand under `library` (same _feature_strand the
    pipeline uses)."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    drop = _drop_flags(keep_dup)
    n = 0
    for c in chroms:
        for r in bam.fetch(c):
            if r.flag & drop or r.mapping_quality < min_mapq:
                continue
            if min_read_length and (r.query_length or 0) < min_read_length:
                continue
            if strand is not None and _feature_strand(r.is_reverse, library) != strand:
                continue
            n += 1
    bam.close()
    return n


def diagnose(bam_path, gene, chroms, *, library="forward", auto_sample_n=None,
            min_mapq=0, keep_dup=False, min_read_length=0, coverage_gap=200,
            region_min_depth_pct=calibrate.DEFAULT_PCT,
            trim_pcts=(TRIM_PCT,)) -> list[dict]:
    """Per-BAM: one pooled (both-strand) auto_region_min_depth call (that
    function never splits by strand -- it pools all islands, see
    calibrate.py), reported identically on both strand rows below; plus
    per-strand total reads / covered-base bg stats.

    `trim_pcts`: bg_trim is computed at every percentage in this list from
    the SAME covered-depth array (one pass), so a trim-sensitivity sweep
    costs nothing extra beyond the first value. The primary bg_trim/fold_trim
    columns use trim_pcts[0] (kept at 2% by default -- unchanged behavior)."""
    total_reads = count_reads(bam_path, chroms, min_mapq=min_mapq, keep_dup=keep_dup,
                              min_read_length=min_read_length)
    factor = scaling_factor(total_reads, 1e6, "rpm") if total_reads else 0.0

    auto_kwargs = dict(library=library, min_mapq=min_mapq, keep_dup=keep_dup,
                       min_read_length=min_read_length, coverage_gap=coverage_gap,
                       pct=region_min_depth_pct)
    if auto_sample_n is not None:
        auto_kwargs["sample_n"] = auto_sample_n
    auto_value, sampled = calibrate.auto_region_min_depth(bam_path, **auto_kwargs)

    rows = []
    for strand in ("+", "-"):
        depths = covered_depths(bam_path, chroms, strand, library=library, min_mapq=min_mapq,
                                keep_dup=keep_dup, min_read_length=min_read_length,
                                coverage_gap=coverage_gap)
        reads_this_strand = count_reads(bam_path, chroms, library=library, strand=strand,
                                        min_mapq=min_mapq, keep_dup=keep_dup,
                                        min_read_length=min_read_length)
        trims_by_pct = {p: background_stats(depths, trim_pct=p)["trim"] for p in trim_pcts}
        bg = background_stats(depths, trim_pct=trim_pcts[0])
        rows.append({
            "gene": gene, "strand": strand,
            "total_reads": reads_this_strand,
            "covered_bases": bg["n"],
            "bg_median": bg["median"], "bg_trim": bg["trim"],
            "auto_value": auto_value,
            "fold_median": (auto_value / bg["median"]) if bg["median"] > 0 else float("inf"),
            "fold_trim": (auto_value / bg["trim"]) if bg["trim"] > 0 else float("inf"),
            "rpm_factor": factor,
            "bg_median_rpm": bg["median"] * factor,
            "bg_trim_rpm": bg["trim"] * factor,
            "auto_value_rpm": auto_value * factor,
            "sampled_chroms_for_auto": ",".join(sampled),
            "bg_trim_by_pct": trims_by_pct,
        })
    return rows


def print_trim_sensitivity(rows: list[dict], trim_pcts) -> None:
    """Table of bg_trim at each trim% plus each row's rank (1 = highest
    bg_trim) at that trim%, and whether ranks shift across trim% choices."""
    print("\n--- trim sensitivity (bg_trim at each trim%, and its rank) ---")
    ranks_by_pct = {}
    for p in trim_pcts:
        ordered = sorted(range(len(rows)), key=lambda i: -rows[i]["bg_trim_by_pct"][p])
        ranks_by_pct[p] = {idx: rank + 1 for rank, idx in enumerate(ordered)}

    label_w = max(len(f"{r['gene']} {r['strand']}") for r in rows)
    val_hdr = "  ".join(f"trim{p:g}%".rjust(10) for p in trim_pcts)
    rank_hdr = "  ".join(f"rank{p:g}%".rjust(8) for p in trim_pcts)
    print(f"{'gene strand'.ljust(label_w)}  {val_hdr}  {rank_hdr}")
    shifts = []
    for i, r in enumerate(rows):
        label = f"{r['gene']} {r['strand']}".ljust(label_w)
        vals = "  ".join(f"{r['bg_trim_by_pct'][p]:.4g}".rjust(10) for p in trim_pcts)
        rks = [ranks_by_pct[p][i] for p in trim_pcts]
        rank_str = "  ".join(str(rk).rjust(8) for rk in rks)
        print(f"{label}  {vals}  {rank_str}")
        shifts.append(max(rks) - min(rks))

    max_shift = max(shifts) if shifts else 0
    n_changed = sum(1 for s in shifts if s > 0)
    print(f"\nrank stability across trim% {list(trim_pcts)}: "
          f"{n_changed}/{len(rows)} rows change rank, max shift = {max_shift} position(s) "
          f"({'ranking is stable' if max_shift <= 1 else 'ranking shifts notably'})")


COLUMNS = ["gene", "strand", "total_reads", "covered_bases", "bg_median", "bg_trim",
          "auto_value", "fold_median", "fold_trim", "rpm_factor", "bg_median_rpm",
          "bg_trim_rpm", "auto_value_rpm", "sampled_chroms_for_auto"]


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def print_table(rows: list[dict]) -> None:
    widths = {c: max(len(c), *(len(_fmt(r[c])) for r in rows)) for c in COLUMNS}
    header = "  ".join(c.ljust(widths[c]) for c in COLUMNS)
    print(header)
    print("-" * len(header))
    for r in rows:
        print("  ".join(_fmt(r[c]).ljust(widths[c]) for c in COLUMNS))


def write_tsv(rows: list[dict], path: str) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r[c]) for c in COLUMNS) + "\n")


def print_judgment(rows: list[dict]) -> None:
    medians = [r["bg_median"] for r in rows]
    flattened = all(m <= 1.0 for m in medians)
    print(f"\n1. bg_median all == 1 across {len(rows)} rows? "
          f"{'YES -- flattened to 1, uninformative as bg' if flattened else 'NO'} "
          f"(values: {[round(m, 3) for m in medians]})")

    by_gene = {}
    for r in rows:
        by_gene.setdefault(r["gene"], []).append(r)
    ratios = []
    for gene, grs in by_gene.items():
        total = sum(r["total_reads"] for r in grs)
        trim = sum(r["bg_trim"] * r["covered_bases"] for r in grs) / max(
            sum(r["covered_bases"] for r in grs), 1)
        ratios.append((gene, total, trim, trim / total if total else float("nan")))
    spread = max(r[3] for r in ratios) / min(r[3] for r in ratios) if ratios else float("nan")
    print(f"2. bg_trim roughly proportional to total_reads across genes? "
          f"ratio bg_trim/total_reads spread (max/min) = {spread:.3g} "
          f"({'roughly proportional' if spread < 3 else 'NOT proportional'})")
    for gene, total, trim, ratio in ratios:
        print(f"   {gene}: total_reads={total} bg_trim(pooled)={trim:.4g} "
              f"ratio={ratio:.4g}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["subset", "full"], default="subset")
    ap.add_argument("--out-tsv", default=None)
    ap.add_argument("--gene", action="append", default=None,
                    help="restrict to these genes (default: all 4)")
    ap.add_argument("--trim-pcts", default=str(TRIM_PCT),
                    help="comma-separated trim percentages for bg_trim sensitivity "
                         "(e.g. 2,10,25); first value is used for the primary "
                         "bg_trim/fold_trim columns (default: %(default)s)")
    args = ap.parse_args(argv)
    trim_pcts = [float(x) for x in args.trim_pcts.split(",")]

    if args.mode == "subset":
        bams, chroms, sample_n = SUBSET_BAMS, SUBSET_CHROMS, SUBSET_AUTO_SAMPLE_N
        out_tsv = args.out_tsv or os.path.expanduser("~/test_subset/diag_bg.tsv")
    else:
        bams, chroms, sample_n = FULL_BAMS, None, None  # None chroms = whole genome
        out_tsv = args.out_tsv or os.path.expanduser("~/test_subset/diag_bg_full.tsv")

    genes = args.gene or list(bams.keys())
    all_rows = []
    for gene in genes:
        bam_path = bams[gene]
        if not os.path.exists(bam_path):
            print(f"error: BAM not found for {gene}: {bam_path}", file=sys.stderr)
            return 1
        gene_chroms = chroms
        if gene_chroms is None:
            bam = pysam.AlignmentFile(bam_path, "rb")
            gene_chroms = list(bam.references)
            bam.close()
        rows = diagnose(bam_path, gene, gene_chroms, library=LIBRARY[gene],
                        auto_sample_n=sample_n, trim_pcts=trim_pcts)
        all_rows.extend(rows)
        print(f"[{gene}] done", file=sys.stderr)

    print_table(all_rows)
    write_tsv(all_rows, out_tsv)
    print(f"\nwrote {out_tsv}")
    print_judgment(all_rows)
    if len(trim_pcts) > 1:
        print_trim_sensitivity(all_rows, trim_pcts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
