"""Auto-calibrate stage-1 thresholds (``region-min-depth``, ``min-region-width``)
from the BAM's own read population, instead of requiring hand-picked absolute
numbers.

Both calibrations are deliberately *descriptive* -- they measure a property the
library actually has (its island depth population, its read length
distribution) and turn it into a threshold. Neither is a significance test:
this caller has no p-value/q-value by design (see README / SPEC_REVIEW.md), and
that applies to the auto-calibration as much as to the detection itself.

(The stage-2 side of the same problem -- humps built from a few heavily
amplified molecules -- is handled by the read positional diversity filter in
``complexity.py``, whose threshold is likewise derived, not tuned.)

Two approaches were tested on real data before landing on this one:

1. A single fixed number (or a value derived from a *different* library's
   tuning) does not transfer across BAMs of different sequencing depth.
2. A per-island *self*-referential percentile (compare an island's positions
   against that same island's own depth distribution) was tested and
   rejected: it always keeps roughly the same fraction of any island
   regardless of whether the whole island is real signal or pure noise. On
   real data, a 1-read/22bp noise island passed 100% of itself under its own
   median -- self-comparison cannot ever discard an island wholesale.

What works: compare each island's own scale (its max depth) against the
population of *all other* islands. On real data (ESRP1, chr1), ~66% of
coverage islands are a single overlapping read (max depth 1) -- unambiguous
background -- while real peak-bearing islands have max depth in the hundreds
to thousands. A modest percentile of this population cleanly separates them
without needing gene annotation or a statistical background model.
"""
from __future__ import annotations

import logging
import math

import numpy as np
import pysam

from .coverage import (
    FLAG_DUP,
    FLAG_QCFAIL,
    FLAG_SECONDARY,
    FLAG_SUPPLEMENTARY,
    FLAG_UNMAPPED,
    chrom_island_blocks,
    chrom_islands,
)

log = logging.getLogger("peakcaller")

DEFAULT_PCT = 75.0     # modest by design: on real data this already excludes
                        # single/double-read noise (~75% of islands) without
                        # cutting into real-but-low-expression loci
DEFAULT_SAMPLE_N = 4
DEFAULT_WIDTH_FACTOR = 2.0      # region must fit this many typical reads side by side
DEFAULT_READ_SAMPLE_N = 50_000  # reads sampled per chromosome for the length distribution
_READ_LEN_PCT = 90.0   # upper-tail read length used as "one read's width" (see below)
DEFAULT_TARGET_FDR = 0.05
DEFAULT_N_SHUFFLES = 3
MIN_REAL_COUNT = 30  # ignore FDR estimates backed by fewer real summits than this
_FIND_REGIONS_FLOOR = 5.0  # low floor used only to enumerate candidate regions for sampling
DEFAULT_SUMMIT_FLOOR_TRIM_PCT = 25.0  # trim% used by auto_summit_floor (--summit-margins)


def _sample_chroms(bam_path: str, n: int) -> list[str]:
    bam = pysam.AlignmentFile(bam_path, "rb")
    refs = list(zip(bam.references, bam.lengths))
    bam.close()
    autosomes = [r for r in refs if r[0].replace("chr", "").isdigit()]
    autosomes.sort(key=lambda r: -r[1])
    chosen = autosomes[:n] if autosomes else refs[:n]
    return [r[0] for r in chosen]


def auto_region_min_depth(
    bam_path: str,
    *,
    library: str = "forward",
    min_mapq: int = 0,
    keep_dup: bool = False,
    min_read_length: int = 0,
    coverage_gap: int = 200,
    pct: float = DEFAULT_PCT,
    sample_n: int = DEFAULT_SAMPLE_N,
) -> tuple[float, list[str]]:
    """Return ``(region_min_depth, sampled_chroms)``.

    ``region_min_depth`` is the ``pct``-th percentile of per-coverage-island
    max depth, measured across a sample of the largest autosomes (or the
    first ``sample_n`` references if none look like autosomes).
    """
    chroms = _sample_chroms(bam_path, sample_n)
    max_depths: list[float] = []
    for c in chroms:
        isls, _n, _abp = chrom_islands(
            bam_path, c, library=library, min_mapq=min_mapq, keep_dup=keep_dup,
            min_read_length=min_read_length, coverage_gap=coverage_gap)
        max_depths.extend(float(arr.max()) for _st, _s, _e, arr in isls)
    if not max_depths:
        return 3.0, chroms  # nothing sampled -- fall back to the old static default
    value = max(float(np.percentile(max_depths, pct)), 1.0)
    return value, chroms


# ---------------------------------------------------------------------------
# min_region_width: derived from the library's own read length distribution
# ---------------------------------------------------------------------------
#
# Why width (and not more depth) is the right stage-1 filter:
#
# A stack of reads that all start at the *same* position -- PCR/optical
# duplicates and near-duplicates, the dominant artifact in CLIP libraries --
# produces a hump purely from read-edge geometry: depth rises at the common
# start, plateaus, and tapers off as the individually shorter reads end. That
# hump is indistinguishable from a real peak by any *relative* measure, because
# it sits on zero baseline: its prominence equals its full height, so
# `prominence >= prominence_frac * region_scale` is satisfied trivially
# (region_scale is that same stack's own depth). No amount of stage-2 tuning
# can reject it -- but its WIDTH gives it away.
#
# The geometric fact: a region covered only by reads sharing one start position
# cannot be wider than the longest read in the library. Measured on real data
# (ESRP1 chr1, region_min_depth=8, 19,751 stage-1 regions): the fraction of
# regions consisting of a single distinct read start is 68% below 30bp, still
# 45-48% at 60-75bp, and drops to exactly 0% at 75bp -- the library's longest
# read is 70bp. The cutoff is not a tuned number; it is where single-start
# stacks become geometrically impossible.
#
# Formula: min_region_width = round(factor * p90(read length)), factor = 2.
#   * p90 rather than max/p99: robust to a handful of outlier-long alignments,
#     which a max would let dictate the threshold for the whole library.
#   * factor 2: the region must be able to hold two typical reads end to end
#     without overlapping, i.e. it must contain genuine positional spread, not
#     one clump. This also lands comfortably past the longest read (ESRP1:
#     2 x 44 = 88bp vs a 70bp longest read), which is exactly where the
#     measurement above shows single-start regions disappear.
#
# Rejected alternatives:
#   * a fixed bp constant (e.g. 50): read length is a per-library property
#     (trimming/adapter settings); the same constant is too strict for a 100bp
#     library and too lax for a 20bp one.
#   * a percentile of the observed *region width* distribution: self-
#     referential in the same way the rejected per-island depth percentile was
#     -- it keeps a fixed fraction of regions no matter how duplicated the
#     library is, and it has no elbow to aim at (the width distribution decays
#     smoothly).
#   * max read length exactly: correct in principle, but one stray long
#     alignment moves the threshold for the entire run.
#
# What this filter does and does not do (measured on ESRP1): it removes ~78% of
# stage-1 regions -- essentially all single-start stacks -- but only ~6% of
# called peaks (whole genome, default settings: 78,636 -> 73,678), because most
# peaks live in regions far wider than any read. It is therefore a
# duplication/artifact filter and a stage-2 work pruner, NOT a peak-count dial:
# overall stringency still comes from --region-min-depth(-pct)
# (ESRP1 p75/p90/p99 -> 73,678 / 49,885 / 20,315 peaks). Raising the factor does
# not change that (factor 2/3/4 -> 73,678 / 71,981 / 70,039), which is why 2 --
# the point where single-start stacks become impossible -- is the default rather
# than something larger.


def _sample_read_lengths(bam_path, chroms, *, min_mapq, keep_dup, min_read_length,
                         reads_per_chrom) -> np.ndarray:
    """Aligned read lengths sampled from the head of each chromosome.

    Read length is a library-preparation property (trimming), not a positional
    one, so sampling the first N reads per chromosome is enough -- verified on
    real data: the first 30k reads of ESRP1 chr1 give a median of 27bp vs 25bp
    for the whole chromosome (847k reads), a 2bp difference on a threshold of
    ~88bp.
    """
    drop = FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_QCFAIL
    if not keep_dup:
        drop |= FLAG_DUP
    lengths: list[int] = []
    bam = pysam.AlignmentFile(bam_path, "rb")
    for c in chroms:
        n = 0
        for r in bam.fetch(c):
            if r.flag & drop or r.mapping_quality < min_mapq:
                continue
            if min_read_length and (r.query_length or 0) < min_read_length:
                continue
            lengths.append(r.query_alignment_length or 0)
            n += 1
            if n >= reads_per_chrom:
                break
    bam.close()
    return np.array([x for x in lengths if x > 0])


def auto_min_region_width(
    bam_path: str,
    *,
    min_mapq: int = 0,
    keep_dup: bool = False,
    min_read_length: int = 0,
    factor: float = DEFAULT_WIDTH_FACTOR,
    sample_n: int = DEFAULT_SAMPLE_N,
    reads_per_chrom: int = DEFAULT_READ_SAMPLE_N,
) -> tuple[int, dict]:
    """Return ``(min_region_width, diagnostics)``.

    The width floor is ``factor`` x the ``_READ_LEN_PCT``-th percentile of this
    BAM's own aligned read lengths -- "a region must be at least this many
    typical reads wide". See the module notes above for why this, and not a
    depth percentile, is what stage 1 should filter on.
    """
    chroms = _sample_chroms(bam_path, sample_n)
    lengths = _sample_read_lengths(bam_path, chroms, min_mapq=min_mapq, keep_dup=keep_dup,
                                   min_read_length=min_read_length,
                                   reads_per_chrom=reads_per_chrom)
    diagnostics = {"sampled_chroms": chroms, "n_reads": int(lengths.size), "factor": factor}
    if lengths.size == 0:
        return 0, diagnostics  # nothing to measure -- leave the filter off
    read_len = float(np.percentile(lengths, _READ_LEN_PCT))
    diagnostics["read_length_pct"] = _READ_LEN_PCT
    diagnostics["read_length"] = read_len
    diagnostics["median_read_length"] = float(np.median(lengths))
    return int(round(factor * read_len)), diagnostics


# ---------------------------------------------------------------------------
# min_prominence: permutation-based empirical FDR (OPT-IN ONLY -- not default)
# ---------------------------------------------------------------------------
#
# Kept because it is implemented and validated, but deliberately NOT the
# default path any more. Rationale (2026-07-29): even though the FDR here is
# empirical (a read-block permutation null) rather than a parametric
# Poisson/NB model, "accept peaks whose FDR <= 5%" is still a statistical
# significance test, which contradicts this project's stated design -- no
# p-value, no q-value, no significance testing. The default pipeline now uses
# a low static min_prominence and lets the read-length-derived
# --min-region-width (above) reject the baseline~=0 stacks that an absolute
# prominence floor used to be needed for. Reach this code by passing
# min_prominence=None (CLI: --min-prominence-auto-fdr).
#
# Three approaches were tried for auto-calibrating min_prominence, in order:
#
# 1. A percentile of the real prominence distribution (as done for
#    region_min_depth above). Rejected: unlike per-island max depth, the real
#    prominence distribution is smooth and heavy-tailed with no natural
#    "noise vs signal" separation at any percentile -- there is no principled
#    percentile to pick, and different BAMs need wildly different ones
#    (RBFOX2 needed ~p99, YBX1 needed ~p98 for a comparable peak count) with
#    no way to know which in advance.
# 2. Knee/elbow detection on the sorted prominence curve. Rejected: the same
#    heavy tail makes the "knee" degenerate to the single most extreme point.
# 3. Permutation null (this one). Build a "null" signal per coverage island
#    by keeping every read block's LENGTH but reassigning its position
#    uniformly at random within the island's own span, then run the exact
#    same prominence detection on it. This measures "what prominence would
#    appear here by pure chance". Compare the real prominence distribution
#    against this null to get an empirical false discovery rate (FDR) at each
#    candidate threshold, and pick the smallest threshold with FDR <=
#    target_fdr for it and every stricter threshold above it.
#
#    A naive first attempt shuffled the raw per-bp depth VALUES instead of
#    read blocks -- rejected, because it also destroys the local smoothness
#    that real reads (which span many bp) naturally create, making the null
#    *noisier* than the real data (FDR > 1 everywhere). Shuffling read BLOCKS
#    (not values) preserves that structural smoothness and only destroys
#    genuine positional clustering, which is what we actually want to test.


def _build_array(start: int, end: int, blocks) -> np.ndarray:
    arr = np.zeros(end - start, dtype=np.int32)
    for bs, be in blocks:
        arr[bs - start:be - start] += 1
    return arr


def _shuffle_blocks(start: int, end: int, blocks, rng: np.random.Generator):
    span = end - start
    shuffled = []
    for bs, be in blocks:
        length = be - bs
        if length >= span:
            shuffled.append((start, end))
            continue
        new_start = start + int(rng.integers(0, span - length + 1))
        shuffled.append((new_start, new_start + length))
    return shuffled


def _sample_prominences(bam_path, chroms, *, library, min_mapq, keep_dup, min_read_length,
                        coverage_gap, region_gap, shuffle, rng) -> np.ndarray:
    from . import core
    from .island import find_islands

    proms: list[float] = []
    for c in chroms:
        isls, _n, _abp = chrom_island_blocks(
            bam_path, c, library=library, min_mapq=min_mapq, keep_dup=keep_dup,
            min_read_length=min_read_length, coverage_gap=coverage_gap)
        for _st, s, e, blocks in isls:
            if shuffle:
                blocks = _shuffle_blocks(s, e, blocks, rng)
            sig = _build_array(s, e, blocks).astype(np.float64)
            regions = find_islands(sig, region_gap=region_gap, min_signal=_FIND_REGIONS_FLOOR)
            if not regions:
                continue
            cores = core.call_cores(
                sig, regions, min_prominence=0.0, prominence_frac=0.0, min_steepness=0.0,
                min_summit_reads=0.0, min_distance=1, min_peak_width=1, rel_height=0.3,
                max_peaks_per_region=0)
            proms.extend(c["prominence"] for c in cores)
    return np.array(proms)


def auto_min_prominence(
    bam_path: str,
    *,
    library: str = "forward",
    min_mapq: int = 0,
    keep_dup: bool = False,
    min_read_length: int = 0,
    coverage_gap: int = 200,
    region_gap: int = 200,
    target_fdr: float = DEFAULT_TARGET_FDR,
    n_shuffles: int = DEFAULT_N_SHUFFLES,
    sample_n: int = DEFAULT_SAMPLE_N,
    seed: int = 0,
) -> tuple[float, dict]:
    """Return ``(min_prominence, diagnostics)`` chosen via a permutation-based
    empirical FDR (see module notes above for why percentile/knee were
    rejected first).
    """
    chroms = _sample_chroms(bam_path, sample_n)
    rng = np.random.default_rng(seed)
    real = _sample_prominences(bam_path, chroms, library=library, min_mapq=min_mapq,
                               keep_dup=keep_dup, min_read_length=min_read_length,
                               coverage_gap=coverage_gap, region_gap=region_gap,
                               shuffle=False, rng=rng)
    nulls = [
        _sample_prominences(bam_path, chroms, library=library, min_mapq=min_mapq,
                            keep_dup=keep_dup, min_read_length=min_read_length,
                            coverage_gap=coverage_gap, region_gap=region_gap,
                            shuffle=True, rng=rng)
        for _ in range(n_shuffles)
    ]
    diagnostics = {
        "sampled_chroms": chroms, "n_real": int(real.size),
        "n_null_mean": float(np.mean([n.size for n in nulls])) if nulls else 0.0,
        "target_fdr": target_fdr, "n_shuffles": n_shuffles,
    }
    if real.size == 0:
        return 5.0, diagnostics

    real_sorted = np.sort(real)
    null_sorted = [np.sort(n) for n in nulls]

    def count_at_least(sorted_arr: np.ndarray, x: float) -> int:
        return len(sorted_arr) - int(np.searchsorted(sorted_arr, x, side="left"))

    # Only trust FDR estimates backed by enough real summits -- at very high
    # thresholds n_real gets tiny and a single stray null summit (pure chance
    # in only 3 shuffles) swings the ratio wildly (checked on real data: FDR
    # jumped from ~0.05 to 1.0 once n_real dropped under ~40). Restricting the
    # search to a well-supported region avoids picking a threshold 10x higher
    # than warranted just because of one noisy point deep in a sparse tail.
    candidates = np.unique(np.round(np.geomspace(1.0, float(real.max()), num=200), 3))
    n_reals = np.array([count_at_least(real_sorted, x) for x in candidates])
    reliable = n_reals >= MIN_REAL_COUNT
    if not reliable.any():
        diagnostics["note"] = "too few real summits for a reliable FDR estimate"
        return max(float(candidates[-1]), 5.0), diagnostics
    candidates = candidates[reliable]

    fdrs = np.empty(candidates.size)
    for i, x in enumerate(candidates):
        n_real = count_at_least(real_sorted, x)
        n_null = sum(count_at_least(ns, x) for ns in null_sorted) / max(n_shuffles, 1)
        fdrs[i] = n_null / n_real

    ok = fdrs <= target_fdr
    diagnostics["target_fdr_achieved"] = bool(ok.any())
    chosen = float(candidates[ok][0]) if ok.any() else float(candidates[-1])
    diagnostics["fdr_at_chosen"] = float(fdrs[candidates == chosen][0])
    diagnostics["n_real_at_chosen"] = int(count_at_least(real_sorted, chosen))
    if not ok.any():
        diagnostics["note"] = (
            f"target FDR {target_fdr:.3g} not reachable within the reliable range "
            f"(n_real >= {MIN_REAL_COUNT}); used the strictest reliable threshold instead "
            f"(achieved FDR {diagnostics['fdr_at_chosen']:.3g}) -- this gene's signal may be "
            f"noisier than others, or --fdr-n-shuffles may need to be higher for a more "
            f"stable estimate")
    return max(chosen, 5.0), diagnostics


# ---------------------------------------------------------------------------
# auto_summit_floor: descriptive background floor for --summit-margins
# ---------------------------------------------------------------------------
#
# Same "descriptive, not a significance test" stance as auto_region_min_depth
# above: this measures the BAM's own covered-base depth population and turns
# it into a raw-read floor, rather than picking a fixed constant or testing
# significance. The canonical "background depth" computation
# (``background_stats``) is shared with scripts/diag_bg.py's stage-1
# region-min-depth diagnostic -- one definition of "bg" for both.


def background_stats(depths: np.ndarray, trim_pct: float = 2.0) -> dict:
    """Median and upper-trimmed mean (drop the top ``trim_pct``% before
    averaging) of a covered-base (depth >= 1) population."""
    if depths.size == 0:
        return {"median": 0.0, "trim": 0.0, "n": 0}
    median = float(np.median(depths))
    cutoff = np.percentile(depths, 100.0 - trim_pct)
    kept = depths[depths <= cutoff]
    trim_mean = float(kept.mean()) if kept.size else median
    return {"median": median, "trim": trim_mean, "n": int(depths.size)}


def auto_summit_floor(
    bam_path: str,
    *,
    library: str = "forward",
    min_mapq: int = 0,
    keep_dup: bool = False,
    min_read_length: int = 0,
    coverage_gap: int = 200,
    trim_pct: float = DEFAULT_SUMMIT_FLOOR_TRIM_PCT,
    sample_n: int = DEFAULT_SAMPLE_N,
) -> tuple[float, dict]:
    """Return ``(floor, diagnostics)``.

    ``floor = ceil(bg_trim)``, where ``bg_trim`` is the ``trim_pct``-th
    upper-trimmed mean (default 25%) of this BAM's own covered-base depth,
    pooled across both strands and a sample of the largest autosomes (same
    sampling as ``auto_region_min_depth``). Raw depth units -- the caller
    (``pipeline.run``) scales this by ``factor`` before comparing against
    normalized signal, the same way it already does for
    ``region_min_depth``/``min_prominence``.
    """
    chroms = _sample_chroms(bam_path, sample_n)
    parts = []
    for c in chroms:
        isls, _n, _abp = chrom_islands(
            bam_path, c, library=library, min_mapq=min_mapq, keep_dup=keep_dup,
            min_read_length=min_read_length, coverage_gap=coverage_gap)
        for _st, _s, _e, arr in isls:
            covered = arr[arr > 0]
            if covered.size:
                parts.append(covered)
    depths = np.concatenate(parts) if parts else np.array([], dtype=np.int64)
    stats = background_stats(depths, trim_pct=trim_pct)
    floor = float(math.ceil(stats["trim"])) if stats["n"] else 1.0
    diagnostics = {
        "sampled_chroms": chroms, "trim_pct": trim_pct, "bg_trim": stats["trim"],
        "bg_median": stats["median"], "n_covered_bases": stats["n"],
    }
    return floor, diagnostics


def tier_for_signal(signal: float, thresholds: list[float]) -> int:
    """Highest 1-indexed tier (T1, T2, ...) whose ``thresholds[i]`` ``signal``
    clears. ``thresholds`` must be ascending; index 0 is assumed always
    satisfied (candidates are generated at that floor in the first place)."""
    return max(i + 1 for i, thr in enumerate(thresholds) if signal >= thr)
