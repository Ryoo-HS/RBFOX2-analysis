"""End-to-end pipeline (memory-efficient, per-island coverage)."""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

from . import calibrate, complexity, context, core, island, normalize
from .coverage import blocks_to_array, chrom_island_blocks, count_library
from .peaks import Peak

log = logging.getLogger("peakcaller")


@dataclass
class Config:
    # normalization
    normalize_method: str = "none"     # none = raw counts (default); rpm = per-million
    scale_to: float = 1e6
    # peak detection (stage 2): topographic prominence
    min_prominence: float | None = 5.0  # absolute prominence floor (signal units). Deliberately
                                        # low: selectivity is meant to come from the relative
                                        # prominence_frac x region_scale test plus stage-1's
                                        # read-length-derived min_region_width, not from this
                                        # floor. None = opt in to the permutation-FDR
                                        # calibration (calibrate.auto_min_prominence) -- kept
                                        # working but no longer the default, since a
                                        # significance test contradicts the project's design.
    target_fdr: float = calibrate.DEFAULT_TARGET_FDR  # used only when min_prominence is None
    fdr_n_shuffles: int = calibrate.DEFAULT_N_SHUFFLES  # used only when min_prominence is None
    prominence_frac: float = 0.0       # region-adaptive: require prom >= frac * region_scale
    region_scale_pct: float = 90.0     # region scale = this percentile of covered depth
    min_steepness: float = 0.0         # min prominence / width (rise per bp)
    min_summit_reads: float = 1.0      # absolute height floor at the summit
    min_distance: int = 15             # min bp between two summits
    min_peak_width: int = 5            # min peak width (bp)
    rel_height: float = 0.3            # boundary width point (0=tip, 1=base); real-data-tuned
    max_peaks_per_region: int = 0      # 0 = unlimited; N = keep top-N per region/gene
    min_complexity: float | None = None  # min read positional diversity of a called peak:
                                        # distinct read start positions (= independent
                                        # molecules) per bp of peak width. None = use
                                        # min_steepness, which width_complexity is a (tight)
                                        # upper bound on once duplicates are collapsed --
                                        # see complexity.py. 0 = off; with the default
                                        # min_steepness=0 the filter is therefore off unless
                                        # the user asks for steepness at all.
    flank: int = 10                    # widen each called peak by this many bp on both sides
                                        # (e.g. for motif search near, not just at, the summit)
    context_window: int = context.DEFAULT_CONTEXT_WINDOW  # bp each side used as the local
                                        # background for the depth_ratio *scoring* component
                                        # (see context.py). Never gates a peak; 0 = skip the
                                        # measurement entirely.
    min_corrected_steepness: float = 0.0  # OFF by default (unlike min_complexity, this is not
                                        # tied to --min-steepness automatically): reject a peak
                                        # if (peak_depth - bg_depth) / (width/2) is below this,
                                        # i.e. min_steepness applied to depth measured against
                                        # the same fixed-width local background as depth_ratio
                                        # instead of the peak's own self-referential region_scale.
                                        # Explicit opt-in because, unlike width/complexity, this
                                        # was tested on real data and only helped 2 of 3
                                        # libraries (ESRP1 inverted) -- see context.py notes and
                                        # project memory. Requires --context-window > 0.
    # region (stage 1)
    region_min_depth: float | None = 3.0   # a region needs depth >= this (enrichment gate).
                                        # None = auto-compute from this BAM's own coverage-island
                                        # depth population (see calibrate.auto_region_min_depth).
    region_min_depth_pct: float = calibrate.DEFAULT_PCT  # percentile used only when auto (None)
    region_gap: int = 200              # stage-1: gap (bp) below region_min_depth merged into one region
    min_region_reads: float = 0.0
    min_region_width: int | None = 0   # drop regions narrower than this many bp (0 = off).
                                        # None = auto-compute from this BAM's own read length
                                        # distribution (see calibrate.auto_min_region_width) --
                                        # this is the primary stage-1 artifact filter: a region
                                        # no wider than one read is a single (near-)duplicate
                                        # read stack, not positional evidence of binding.
    min_region_width_factor: float = calibrate.DEFAULT_WIDTH_FACTOR  # used only when auto (None)
    # coverage building (stage 0, memory chunking -- not a biological parameter)
    coverage_gap: int = 200            # merge raw read blocks within this many bp into one array chunk
    # BAM handling
    library: str = "forward"
    min_mapq: int = 0
    keep_dup: bool = False
    min_read_length: int = 0
    meta: dict = field(default_factory=dict)


@dataclass
class Result:
    peaks: list
    factor: float
    total_reads: int = 0
    strands: tuple = ("+", "-")


def run(bam_path: str, cfg: Config, *, chroms: list[str] | None = None) -> Result:
    """Call peaks with bounded memory: coverage is built per chromosome as small
    per-island arrays (never a whole-chromosome array)."""
    import pysam

    bam = pysam.AlignmentFile(bam_path, "rb")
    all_chroms = list(bam.references)
    chrom_lengths = dict(zip(bam.references, bam.lengths))
    bam.close()
    if chroms is not None:
        missing = [c for c in chroms if c not in all_chroms]
        if missing:
            raise ValueError(f"chromosome(s) not found in BAM: {', '.join(missing)}")
    wanted = list(dict.fromkeys(chroms)) if chroms is not None else all_chroms
    strands = (".",) if cfg.library == "unstranded" else ("+", "-")

    if cfg.coverage_gap < cfg.region_gap:
        log.warning(
            "--coverage-gap (%d) < --region-gap (%d): coverage array chunks are built "
            "before stage-1 region merging, so a smaller --coverage-gap silently caps "
            "how far apart two regions can ever be merged, regardless of --region-gap. "
            "Recommend --coverage-gap >= --region-gap.",
            cfg.coverage_gap, cfg.region_gap)

    if cfg.min_corrected_steepness > 0 and not cfg.context_window:
        log.warning(
            "--min-corrected-steepness=%.4g has no effect with --context-window=0: it needs "
            "the local-background measurement that flag disables.",
            cfg.min_corrected_steepness)

    if cfg.normalize_method == "rpm" or cfg.min_region_reads:
        total_reads, total_aligned_bp = count_library(
            bam_path, min_mapq=cfg.min_mapq, keep_dup=cfg.keep_dup,
            min_read_length=cfg.min_read_length, chroms=wanted)
    else:
        total_reads, total_aligned_bp = 0, 0
    factor = normalize.scaling_factor(total_reads, cfg.scale_to, cfg.normalize_method)
    mrl = (total_aligned_bp / total_reads) if total_reads else 1.0
    log.info("scaling_factor=%.6g", factor)

    if cfg.region_min_depth is None:
        region_min_depth, sampled = calibrate.auto_region_min_depth(
            bam_path, library=cfg.library, min_mapq=cfg.min_mapq, keep_dup=cfg.keep_dup,
            min_read_length=cfg.min_read_length, coverage_gap=cfg.coverage_gap,
            pct=cfg.region_min_depth_pct)
        region_min_depth *= factor  # raw-depth-based -> same units as `signal` below
        log.info("auto-computed region-min-depth=%.4g (region-min-depth-pct=%.1f, "
                  "sampled chroms=%s)", region_min_depth, cfg.region_min_depth_pct, sampled)
    else:
        region_min_depth = cfg.region_min_depth

    if cfg.min_region_width is None:
        # bp, not signal units -- deliberately NOT multiplied by `factor`
        # (normalization rescales depth, never genomic width).
        min_region_width, width_diag = calibrate.auto_min_region_width(
            bam_path, min_mapq=cfg.min_mapq, keep_dup=cfg.keep_dup,
            min_read_length=cfg.min_read_length, factor=cfg.min_region_width_factor)
        log.info("auto-computed min-region-width=%d bp (%.4g x p%.0f read length %.4g, "
                 "median read length %.4g, %d reads sampled from %s)",
                 min_region_width, cfg.min_region_width_factor,
                 width_diag.get("read_length_pct", 0.0), width_diag.get("read_length", 0.0),
                 width_diag.get("median_read_length", 0.0), width_diag.get("n_reads", 0),
                 width_diag.get("sampled_chroms"))
    else:
        min_region_width = cfg.min_region_width

    if cfg.min_prominence is None:
        min_prominence, fdr_diag = calibrate.auto_min_prominence(
            bam_path, library=cfg.library, min_mapq=cfg.min_mapq, keep_dup=cfg.keep_dup,
            min_read_length=cfg.min_read_length, coverage_gap=cfg.coverage_gap,
            region_gap=cfg.region_gap, target_fdr=cfg.target_fdr, n_shuffles=cfg.fdr_n_shuffles)
        min_prominence *= factor  # raw-depth-based -> same units as `signal` below
        log_fn = log.info if fdr_diag.get("target_fdr_achieved", True) else log.warning
        log_fn("auto-computed min-prominence=%.4g (target-fdr=%.3g, achieved-fdr=%.3g, "
               "n-shuffles=%d, sampled chroms=%s)", min_prominence, cfg.target_fdr,
               fdr_diag.get("fdr_at_chosen", float("nan")), cfg.fdr_n_shuffles,
               fdr_diag.get("sampled_chroms"))
        if fdr_diag.get("note"):
            log.warning(fdr_diag["note"])
    else:
        min_prominence = cfg.min_prominence

    min_complexity = complexity.resolve_min_complexity(
        cfg.min_complexity, cfg.min_steepness, factor)
    if min_complexity > 0:
        log.info("min-complexity=%.4g distinct read starts per bp%s", min_complexity,
                 "" if cfg.min_complexity is not None else " (= --min-steepness)")

    if cfg.min_corrected_steepness > 0:
        log.info("min-corrected-steepness=%.4g (peak_depth - bg_depth) / (width/2), raw "
                 "depth units", cfg.min_corrected_steepness)

    peaks: list[Peak] = []
    seen_reads = 0
    n_low_complexity = 0
    n_weak_isolated = 0
    for chrom in sorted(wanted):
        # blocks (not just the summed array) so stage-2 output can be filtered on
        # how many *distinct* read start positions -- i.e. independent molecules --
        # produced each peak (see complexity.py). Same memory profile: chrom_islands
        # groups the very same blocks internally before summing them.
        isls, n, _abp = chrom_island_blocks(
            bam_path, chrom, library=cfg.library, min_mapq=cfg.min_mapq,
            keep_dup=cfg.keep_dup, min_read_length=cfg.min_read_length,
            coverage_gap=cfg.coverage_gap)
        seen_reads += n
        chrom_peaks: list[Peak] = []
        blocks_by_strand: dict[str, list] = {}
        for st, cs, ce, blocks in isls:
            if cfg.context_window:
                blocks_by_strand.setdefault(st, []).extend(blocks)
            arr = blocks_to_array(cs, ce, blocks)
            signal = normalize.to_rpm(arr, factor)
            sub = island.find_islands(signal, region_gap=cfg.region_gap,
                                      min_signal=region_min_depth)
            if cfg.min_region_reads or min_region_width:
                qualified = []
                for s0, e0 in sub:
                    approx_reads = float(arr[s0:e0].sum()) / mrl if mrl else 0.0
                    if approx_reads < cfg.min_region_reads:
                        continue
                    if (e0 - s0) < min_region_width:
                        continue
                    qualified.append((s0, e0))
                sub = qualified
            if not sub:
                continue
            cores = core.call_cores(
                signal, sub,
                min_prominence=min_prominence, prominence_frac=cfg.prominence_frac,
                region_scale_pct=cfg.region_scale_pct, min_steepness=cfg.min_steepness,
                min_summit_reads=cfg.min_summit_reads, min_distance=cfg.min_distance,
                min_peak_width=cfg.min_peak_width, rel_height=cfg.rel_height,
                max_peaks_per_region=cfg.max_peaks_per_region)
            # Post-filter, so core.py's prominence math is untouched. Note this
            # runs *after* max_peaks_per_region has already taken the top N, so a
            # region can end up with fewer than N peaks if some of its strongest
            # were duplicate-driven -- which is the intended reading of "N best".
            if min_complexity > 0 and cores:
                n_before = len(cores)
                cores = complexity.filter_cores(cores, cs, blocks, min_complexity)
                n_low_complexity += n_before - len(cores)
            for c in cores:
                chrom_peaks.append(Peak(
                    chrom=chrom, start=cs + c["start"], end=cs + c["end"], strand=st,
                    summit=cs + c["summit"], signal=c["signal"], baseline=c["baseline"],
                    fold=c["fold"], region_scale=c.get("region_scale", 0.0)))

        # --flank first, then depth_ratio, so the local-background measurement
        # describes the interval that is actually written out. Both are done per
        # chromosome, while this chromosome's read blocks are still in memory --
        # depth_ratio needs a fixed +-context_window view that the per-island
        # arrays (chunked by --coverage-gap) cannot supply, but it does not need
        # a second pass over the BAM.
        if cfg.flank:
            chrom_len = chrom_lengths.get(chrom)
            for p in chrom_peaks:
                p.start = max(0, p.start - cfg.flank)
                p.end = (p.end + cfg.flank if chrom_len is None
                         else min(chrom_len, p.end + cfg.flank))
        if cfg.context_window:
            context.annotate(chrom_peaks, blocks_by_strand, window=cfg.context_window)
            if cfg.min_corrected_steepness > 0:
                n_before = len(chrom_peaks)
                chrom_peaks = [
                    p for p in chrom_peaks
                    if (p.peak_depth - p.bg_depth) / max((p.end - p.start) / 2.0, 1.0)
                    >= cfg.min_corrected_steepness
                ]
                n_weak_isolated += n_before - len(chrom_peaks)
        peaks.extend(chrom_peaks)
        del blocks_by_strand, isls

    if n_low_complexity:
        log.info("dropped %d peak(s) with fewer than %.4g distinct read starts per bp "
                 "(duplicate-driven humps; --min-complexity)",
                 n_low_complexity, min_complexity)
    if n_weak_isolated:
        log.info("dropped %d peak(s) below --min-corrected-steepness=%.4g "
                 "(weak against their local background)",
                 n_weak_isolated, cfg.min_corrected_steepness)
    peaks.sort(key=lambda p: (p.chrom, p.start, p.end, p.strand))
    for i, p in enumerate(peaks, 1):
        p.name = f"peak_{i}"
    log.info("called %d peaks", len(peaks))
    return Result(peaks=peaks, factor=factor,
                  total_reads=total_reads or seen_reads, strands=strands)
