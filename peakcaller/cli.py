"""Command-line interface."""
from __future__ import annotations

import argparse
import logging
import os
import sys

from . import __version__
from .output import write_bed6, write_bedgraph, write_narrowpeak
from .pipeline import Config, run


def build_parser():
    p = argparse.ArgumentParser(
        prog="peakcaller",
        description="Rule-based (prominence) peak caller for CLIP-seq / RBP data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--bam", required=True, help="input BAM (coordinate sorted + indexed)")
    p.add_argument("-o", "--output", required=True, help="output peak file")
    p.add_argument("--format", choices=["narrowPeak", "bed"], default="narrowPeak")
    g = p.add_argument_group("peak detection (prominence)")
    g.add_argument("--min-prominence", type=float, default=Config().min_prominence,
                   help="absolute floor: rise above flanking valleys (signal units); "
                        "kept deliberately low -- selectivity comes from --prominence-frac "
                        "(relative to the region's own scale) and from stage-1's "
                        "--min-region-width, not from this floor")
    g.add_argument("--prominence-frac", type=float, default=0.0,
                   help="region-adaptive: prominence >= frac x region scale (e.g. 0.2)")
    g.add_argument("--region-scale-pct", type=float, default=90.0,
                   help="percentile of region covered depth used as its scale")
    g.add_argument("--min-steepness", type=float, default=0.0,
                   help="min prominence/width; drops gentle bumps (e.g. 0.5)")
    g.add_argument("--min-summit-reads", type=float, default=1.0, help="height floor at summit")
    g.add_argument("--summit-margins", type=str, default=None,
                   help="comma-separated ascending non-negative-int margins (raw reads, "
                        "max 3, e.g. 18,48,98 -- roughly 20/50/100 at a typical floor=2, "
                        "spacing chosen from observed summit distributions, not "
                        "significance thresholds) above the auto background floor "
                        "(ceil of the trimmed-mean covered-base depth at trim=25%% of "
                        "this BAM -- see calibrate.auto_summit_floor). Candidates are "
                        "generated once at the T1 threshold; each called peak is tagged "
                        "with the highest tier (T1..Tn) its summit height clears and the "
                        "name column gets a _T<n> suffix. AND'd with --min-summit-reads: "
                        "effective T1 = max(--min-summit-reads, floor + margins[0]). "
                        "Default off (no tiers; existing single-threshold behavior)")
    g.add_argument("--split-tiers", action="store_true",
                   help="also write one output file per tier (<output>.T1.<ext> etc.) "
                        "in addition to the combined file; requires --summit-margins")
    g.add_argument("--context-window", type=int, default=1000,
                   help="bp each side used as local background for the depth_ratio "
                        "SCORING component: peak mean depth / mean depth of the same-strand "
                        "window with all called peaks masked out. Never filters a peak, only "
                        "ranks it. Unlike --prominence-frac's region scale this window is a "
                        "fixed width, so an isolated hump cannot become its own baseline. "
                        "0 = skip the measurement")
    g.add_argument("--min-complexity", type=float, default=None,
                   help="min read positional diversity of a called peak: distinct read "
                        "start positions (= independent molecules, PCR duplicates collapse "
                        "to one) per bp of peak width; rejects tall humps made of few "
                        "heavily amplified molecules. 0 = off; default: --min-steepness, "
                        "which this quantity is an upper bound on once duplicates are "
                        "collapsed (so a peak below it provably cannot stay steep enough)")
    g.add_argument("--min-corrected-steepness", type=float,
                   default=Config().min_corrected_steepness,
                   help="reject a peak if (peak_depth - bg_depth) / (width/2) is below this, "
                        "where bg_depth is the same fixed-window local background used by "
                        "--context-window's depth_ratio score. Unlike --min-complexity this is "
                        "OFF by default (0) and not tied to --min-steepness: real-data testing "
                        "found it helps for some libraries but inverts on others, so it is an "
                        "explicit opt-in, not a default gate. Raw depth units (unaffected by "
                        "--normalize-method). Needs --context-window > 0.")
    g.add_argument("--min-distance", type=int, default=15, help="min bp between summits")
    g.add_argument("--rel-height", type=float, default=0.3,
                   help="boundary width point (0=tip,1=base); lower=tighter")
    g.add_argument("--max-peaks-per-region", type=int, default=0,
                   help="0=unlimited; N=top-N per region/gene")
    g.add_argument("--min-peak-width", type=int, default=5, help="min peak width (bp)")
    g.add_argument("--flank", type=int, default=10,
                   help="widen each called peak by N bp on both sides in the output "
                        "(e.g. for motif search near, not just at, the summit); "
                        "does not affect detection")
    g = p.add_argument_group("region (stage 1)")
    g.add_argument("--region-min-depth", type=float, default=None,
                   help="a region needs depth >= this (enrichment gate); "
                        "default: auto-computed from this BAM's own coverage "
                        "(see --region-min-depth-pct)")
    g.add_argument("--region-min-depth-pct", type=float, default=Config().region_min_depth_pct,
                   help="percentile of per-coverage-island max depth used to "
                        "auto-compute --region-min-depth; only used when "
                        "--region-min-depth is not given")
    g.add_argument("--region-gap", type=int, default=200,
                   help="max gap (bp, depth < region-min-depth) merged into one region")
    g.add_argument("--min-region-reads", type=float, default=0.0)
    g.add_argument("--min-region-width", type=int, default=None,
                   help="drop regions narrower than N bp (0 = off); default: auto-computed "
                        "from this BAM's own read lengths -- a region no wider than one read "
                        "is a single (near-)duplicate read stack, not evidence of binding "
                        "(see --min-region-width-factor)")
    g.add_argument("--min-region-width-factor", type=float,
                   default=Config().min_region_width_factor,
                   help="auto --min-region-width = this many read lengths (90th percentile "
                        "of this BAM's aligned read lengths); only used when "
                        "--min-region-width is not given")
    g = p.add_argument_group(
        "opt-in: empirical-FDR calibration of --min-prominence (not the default; this "
        "caller otherwise computes no statistical significance at all)")
    g.add_argument("--min-prominence-auto-fdr", action="store_true",
                   help="ignore --min-prominence and instead pick the smallest prominence "
                        "floor whose empirical FDR against a read-block-shuffled null is "
                        "<= --target-fdr")
    g.add_argument("--target-fdr", type=float, default=Config().target_fdr,
                   help="target empirical FDR; only used with --min-prominence-auto-fdr")
    g.add_argument("--fdr-n-shuffles", type=int, default=Config().fdr_n_shuffles,
                   help="number of shuffled nulls averaged for the FDR estimate; "
                        "only used with --min-prominence-auto-fdr")
    g = p.add_argument_group("coverage building (memory chunking, not biological)")
    g.add_argument("--coverage-gap", type=int, default=200,
                   help="merge raw read blocks within this many bp into one coverage "
                        "chunk; should be >= --region-gap")
    g = p.add_argument_group("normalization")
    g.add_argument("--normalize-method", choices=["none", "rpm"], default="none",
                   help="none=raw counts (default); rpm=per-million (MACS SPMR)")
    g.add_argument("--scale-to", type=float, default=1e6)
    g = p.add_argument_group("BAM handling")
    g.add_argument("--library-type", choices=["forward", "reverse", "unstranded"],
                   default="forward", dest="library")
    g.add_argument("--min-mapq", type=int, default=0)
    g.add_argument("--min-read-length", type=int, default=0)
    g.add_argument("--keep-dup", action="store_true")
    g = p.add_argument_group("extra outputs")
    g.add_argument("--bdg", action="store_true", help="also write per-strand bedGraph")
    p.add_argument("--chrom", action="append", default=None)
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def _config_from_args(a):
    return Config(
        normalize_method=a.normalize_method, scale_to=a.scale_to,
        min_prominence=None if a.min_prominence_auto_fdr else a.min_prominence,
        target_fdr=a.target_fdr,
        fdr_n_shuffles=a.fdr_n_shuffles, prominence_frac=a.prominence_frac,
        region_scale_pct=a.region_scale_pct, min_steepness=a.min_steepness,
        min_summit_reads=a.min_summit_reads, summit_margins=a.summit_margins,
        split_tiers=a.split_tiers, min_distance=a.min_distance,
        min_peak_width=a.min_peak_width, flank=a.flank, min_complexity=a.min_complexity,
        context_window=a.context_window, min_corrected_steepness=a.min_corrected_steepness,
        rel_height=a.rel_height, max_peaks_per_region=a.max_peaks_per_region,
        region_min_depth=a.region_min_depth, region_min_depth_pct=a.region_min_depth_pct,
        region_gap=a.region_gap,
        min_region_reads=a.min_region_reads, min_region_width=a.min_region_width,
        min_region_width_factor=a.min_region_width_factor, coverage_gap=a.coverage_gap,
        library=a.library, min_mapq=a.min_mapq, keep_dup=a.keep_dup,
        min_read_length=a.min_read_length)


def _write_bedgraphs(bam_path, args, result):
    from .coverage import stream_bedgraph
    base = os.path.splitext(args.output)[0]
    paths = stream_bedgraph(bam_path, base, library=args.library, min_mapq=args.min_mapq,
                            keep_dup=args.keep_dup, min_read_length=args.min_read_length,
                            coverage_gap=args.coverage_gap, factor=result.factor,
                            chroms=args.chrom)
    for p in paths:
        logging.getLogger("peakcaller").info("wrote %s", p)


def _parse_summit_margins(s: str, parser):
    try:
        vals = [int(x) for x in s.split(",")]
    except ValueError:
        parser.error(f"--summit-margins must be comma-separated integers (got {s!r})")
    if not 1 <= len(vals) <= 3:
        parser.error(f"--summit-margins must have 1-3 values (got {len(vals)}: {vals})")
    if any(v < 0 for v in vals):
        parser.error(f"--summit-margins values must be >= 0 (got {vals})")
    if any(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)):
        parser.error(f"--summit-margins must be strictly ascending (got {vals})")
    return tuple(vals)


def _validate(parser, a):
    nn = {"prominence-frac": a.prominence_frac,
          "min-steepness": a.min_steepness, "min-summit-reads": a.min_summit_reads,
          "min-distance": a.min_distance, "min-peak-width": a.min_peak_width, "flank": a.flank,
          "region-gap": a.region_gap, "coverage-gap": a.coverage_gap,
          "min-region-reads": a.min_region_reads,
          "min-region-width-factor": a.min_region_width_factor,
          "min-mapq": a.min_mapq, "min-read-length": a.min_read_length,
          "max-peaks-per-region": a.max_peaks_per_region, "scale-to": a.scale_to,
          "min-corrected-steepness": a.min_corrected_steepness}
    if a.region_min_depth is not None:
        nn["region-min-depth"] = a.region_min_depth
    if a.min_region_width is not None:
        nn["min-region-width"] = a.min_region_width
    if a.min_prominence is not None:
        nn["min-prominence"] = a.min_prominence
    if a.min_complexity is not None:
        nn["min-complexity"] = a.min_complexity
    for name, val in nn.items():
        if val < 0:
            parser.error(f"--{name} must be >= 0 (got {val})")
    if not 0.0 <= a.rel_height <= 1.0:
        parser.error(f"--rel-height must be between 0 and 1 (got {a.rel_height})")
    if not 0.0 <= a.region_scale_pct <= 100.0:
        parser.error(f"--region-scale-pct must be between 0 and 100 (got {a.region_scale_pct})")
    if not 0.0 <= a.region_min_depth_pct <= 100.0:
        parser.error(f"--region-min-depth-pct must be between 0 and 100 "
                     f"(got {a.region_min_depth_pct})")
    if not 0.0 < a.target_fdr <= 1.0:
        parser.error(f"--target-fdr must be between 0 (exclusive) and 1 (got {a.target_fdr})")
    if a.fdr_n_shuffles < 1:
        parser.error(f"--fdr-n-shuffles must be >= 1 (got {a.fdr_n_shuffles})")
    if a.summit_margins is not None:
        a.summit_margins = _parse_summit_margins(a.summit_margins, parser)
    if a.split_tiers and a.summit_margins is None:
        parser.error("--split-tiers requires --summit-margins")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    _validate(parser, args)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="[%(levelname)s] %(message)s")
    if not os.path.exists(args.bam):
        print(f"error: BAM not found: {args.bam}", file=sys.stderr)
        return 1
    try:
        result = run(args.bam, _config_from_args(args), chroms=args.chrom)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    writer = write_narrowpeak if args.format == "narrowPeak" else write_bed6
    header_lines = None
    if result.summit_tier_diag:
        d = result.summit_tier_diag
        header_lines = [
            f"summit-margins floor={d['floor_raw']:.4g} raw "
            f"(bg_trim@trim={d['trim_pct']:.0f}%={d['bg_trim_raw']:.4g} raw)",
            f"summit-margins margins={d['margins']} raw_thresholds={d['raw_thresholds']}",
            f"summit-margins effective_thresholds(signal units)="
            f"{[round(t, 4) for t in d['effective_thresholds']]}",
        ]
    with open(args.output, "w") as fh:
        n = writer(result.peaks, fh, header_lines=header_lines)
    print(f"{n} peaks written to {args.output}")
    if args.split_tiers and result.summit_tier_diag:
        base, ext = os.path.splitext(args.output)
        n_tiers = len(result.summit_tier_diag["effective_thresholds"])
        for t in range(1, n_tiers + 1):
            tier_peaks = [p for p in result.peaks if p.tier == t]
            tier_path = f"{base}.T{t}{ext}"
            with open(tier_path, "w") as fh:
                writer(tier_peaks, fh, header_lines=header_lines)
            print(f"{len(tier_peaks)} T{t} peaks written to {tier_path}")
    if args.bdg:
        _write_bedgraphs(args.bam, args, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
