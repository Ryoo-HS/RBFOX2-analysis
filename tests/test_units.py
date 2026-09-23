import numpy as np
import pytest

from rbpc import core, island, normalize
from rbpc.cli import build_parser
from rbpc.pipeline import Config


def test_flank_default_is_10bp():
    assert Config().flank == 10
    args = build_parser().parse_args(["--bam", "x.bam", "-o", "out.narrowPeak"])
    assert args.flank == 10


def test_cli_defaults_auto_stage1_but_static_min_prominence():
    args = build_parser().parse_args(["--bam", "x.bam", "-o", "out.narrowPeak"])
    # stage 1 auto-calibrates from the BAM itself (None = auto)...
    assert args.region_min_depth is None
    assert args.min_region_width is None
    assert args.min_region_width_factor == 2.0
    # ...while min-prominence is a low static floor again: the permutation-FDR
    # calibration is opt-in only (a significance test is against this project's design).
    assert args.min_prominence == 5.0
    assert args.min_prominence_auto_fdr is False


def test_cli_min_prominence_auto_fdr_is_opt_in():
    from rbpc.cli import _config_from_args
    p = build_parser()
    assert _config_from_args(p.parse_args(
        ["--bam", "x.bam", "-o", "o"])).min_prominence == 5.0
    assert _config_from_args(p.parse_args(
        ["--bam", "x.bam", "-o", "o", "--min-prominence-auto-fdr"])).min_prominence is None


def test_scaling_factor_rpm():
    assert normalize.scaling_factor(2_000_000, 1e6, "rpm") == 0.5
    assert normalize.scaling_factor(0, 1e6, "rpm") == 0.0
    assert normalize.scaling_factor(123, 1e6, "none") == 1.0


def test_find_islands_merges_within_gap():
    sig = np.zeros(100)
    sig[10:20] = 1
    sig[25:30] = 1
    sig[60:70] = 1
    assert island.find_islands(sig, region_gap=10) == [(10, 30), (60, 70)]
    assert island.find_islands(sig, region_gap=2) == [(10, 20), (25, 30), (60, 70)]


def test_find_islands_min_signal_is_inclusive():
    # depth exactly equal to min_signal must count as covered ("depth >= this").
    sig = np.zeros(20)
    sig[5:15] = 3.0
    assert island.find_islands(sig, region_gap=200, min_signal=3.0) == [(5, 15)]
    assert island.find_islands(sig, region_gap=200, min_signal=3.1) == []


def test_find_islands_zero_gap_not_covered_even_at_min_signal_zero():
    # min_signal=0 ("no depth floor") must still treat true zero-depth as a gap,
    # not silently merge everything into one island.
    sig = np.zeros(300)
    sig[10:20] = 1
    sig[150:160] = 1
    assert island.find_islands(sig, region_gap=50, min_signal=0.0) == [(10, 20), (150, 160)]


def test_prominence_keeps_high_bump_drops_low():
    x = np.arange(1000)
    sig = np.maximum(0, 30 - np.abs(x - 300) * 0.8)   # prominence ~30
    sig = sig + np.maximum(0, 5 - np.abs(x - 700) * 0.5)  # prominence ~5
    cores = core.call_cores(sig, [(0, 1000)], min_prominence=20, min_summit_reads=1,
                            min_distance=15, min_peak_width=3)
    assert len(cores) == 1
    assert abs(cores[0]["summit"] - 300) < 20


def test_prominence_splits_two_mountains():
    x = np.arange(1000)
    sig = (np.maximum(0, 40 - np.abs(x - 300) * 0.8)
           + np.maximum(0, 40 - np.abs(x - 500) * 0.8))
    cores = core.call_cores(sig, [(0, 1000)], min_prominence=20, min_summit_reads=1,
                            min_distance=15, min_peak_width=3)
    summits = sorted(c["summit"] for c in cores)
    assert len(cores) == 2
    assert abs(summits[0] - 300) < 20 and abs(summits[1] - 500) < 20


def test_width_complexity_counts_distinct_starts_per_bp():
    from rbpc.complexity import block_bounds, width_complexity
    # 6 reads, but only 2 distinct start positions -> 2 molecules, amplified.
    blocks = [(100, 130)] * 4 + [(110, 140)] * 2
    starts, ends = block_bounds(blocks)
    cx, nd = width_complexity(starts, ends, 100, 120)
    assert nd == 2
    assert cx == 2 / 20
    # a block that ends before the peak starts must not be counted
    starts, ends = block_bounds([(0, 50), (100, 130)])
    cx, nd = width_complexity(starts, ends, 100, 110)
    assert nd == 1


def test_width_complexity_empty_and_non_overlapping():
    from rbpc.complexity import block_bounds, width_complexity
    assert width_complexity(*block_bounds([]), 0, 10) == (0.0, 0)
    assert width_complexity(*block_bounds([(500, 530)]), 0, 10) == (0.0, 0)


def test_min_complexity_defaults_to_min_steepness():
    from rbpc.complexity import resolve_min_complexity
    # None = "the same bar, applied to molecules instead of reads"
    assert resolve_min_complexity(None, 0.5) == 0.5
    # ...so with the default min_steepness=0 the filter stays off (backward compatible)
    assert resolve_min_complexity(None, 0.0) == 0.0
    # an explicit value wins, and 0 turns it off
    assert resolve_min_complexity(1.5, 0.5) == 1.5
    assert resolve_min_complexity(0.0, 0.5) == 0.0


def test_min_complexity_is_converted_to_molecule_units_under_rpm():
    from rbpc.complexity import resolve_min_complexity
    # min_steepness is in signal units per bp; under rpm one molecule is worth
    # `factor` signal units, so the molecules-per-bp floor is steepness/factor.
    assert resolve_min_complexity(None, 0.5, factor=0.5) == 1.0
    assert resolve_min_complexity(None, 0.5, factor=1.0) == 0.5


def test_cli_min_complexity_defaults_to_auto():
    args = build_parser().parse_args(["--bam", "x.bam", "-o", "out.narrowPeak"])
    assert args.min_complexity is None      # None = tied to --min-steepness
    args = build_parser().parse_args(
        ["--bam", "x.bam", "-o", "out.narrowPeak", "--min-complexity", "0.8"])
    assert args.min_complexity == 0.8


# --- local background contrast (context.py / depth_ratio) --------------------

def _idx(blocks):
    from rbpc.context import build_block_index
    return build_block_index(blocks)


def test_window_depth_matches_naive_block_sum():
    from rbpc.context import window_depth
    blocks = [(10, 20), (12, 18), (30, 40)]
    starts, ends, mx = _idx(blocks)
    arr = window_depth(starts, ends, mx, 0, 50)
    naive = np.zeros(50, dtype=np.int32)
    for s, e in blocks:
        naive[s:e] += 1
    assert np.array_equal(arr, naive)


def test_depth_ratio_is_peak_mean_over_background_mean():
    from rbpc.context import depth_ratio
    # peak [100,110) covered 10 deep; background is a flat 1-deep carpet.
    blocks = [(100, 110)] * 10 + [(0, 300)]
    starts, ends, mx = _idx(blocks)
    empty = (np.array([100], dtype=np.int64), np.array([110], dtype=np.int64))
    r, pk, bg = depth_ratio(starts, ends, mx, 100, 110, empty, window=100)
    assert pk == pytest.approx(11.0)      # 10 stacked + the carpet
    assert bg == pytest.approx(1.0)       # carpet only, peak masked out
    assert r == pytest.approx(11.0)


def test_depth_ratio_masks_other_called_peaks_out_of_the_background():
    """The CDH2 case: a neighbouring real peak must not inflate the background
    and drag its own cluster's contrast down."""
    from rbpc.context import depth_ratio
    blocks = [(100, 110)] * 10 + [(200, 210)] * 10 + [(0, 300)]
    starts, ends, mx = _idx(blocks)
    alone = (np.array([100], dtype=np.int64), np.array([110], dtype=np.int64))
    both = (np.array([100, 200], dtype=np.int64),
            np.array([110, 210], dtype=np.int64))
    r_unmasked, _, bg_unmasked = depth_ratio(starts, ends, mx, 100, 110, alone, window=100)
    r_masked, _, bg_masked = depth_ratio(starts, ends, mx, 100, 110, both, window=100)
    assert bg_masked < bg_unmasked        # neighbour removed from denominator
    assert r_masked > r_unmasked


def test_depth_ratio_is_scale_free_under_uniform_depth_change():
    """A ratio of depths, so sequencing depth / normalization cancels -- the
    property an absolute magnitude floor lacks."""
    from rbpc.context import depth_ratio
    base = [(100, 110)] * 5 + [(0, 300)]
    starts, ends, mx = _idx(base)
    deep_starts, deep_ends, deep_mx = _idx(base * 7)
    bounds = (np.array([100], dtype=np.int64), np.array([110], dtype=np.int64))
    r1, _, _ = depth_ratio(starts, ends, mx, 100, 110, bounds, window=100)
    r2, _, _ = depth_ratio(deep_starts, deep_ends, deep_mx, 100, 110, bounds, window=100)
    assert r1 == pytest.approx(r2)


def test_depth_ratio_reports_inf_when_no_background_remains():
    from rbpc.context import depth_ratio
    starts, ends, mx = _idx([(100, 110)] * 3)
    bounds = (np.array([100], dtype=np.int64), np.array([110], dtype=np.int64))
    r, _, bg = depth_ratio(starts, ends, mx, 100, 110, bounds, window=50)
    assert bg == 0.0 and r == float("inf")


def test_score_weights_sum_to_one_and_include_depth_ratio():
    import inspect

    from rbpc import output
    src = inspect.getsource(output._scores)
    assert "0.40 * p_comp + 0.15 * r_comp + 0.25 * s_comp + 0.20 * c_comp" in src


def test_score_rises_with_depth_ratio_all_else_equal():
    from rbpc.output import _scores
    from rbpc.peaks import Peak

    def mk(dr):
        return Peak(chrom="c", start=0, end=20, strand="+", summit=10, signal=100.0,
                    baseline=1.0, fold=100.0, region_scale=50.0, depth_ratio=dr)
    lo, hi = _scores([mk(1.0), mk(50.0)])
    assert hi > lo


def test_infinite_depth_ratio_does_not_poison_the_reference():
    from rbpc.output import _scores
    from rbpc.peaks import Peak

    def mk(dr):
        return Peak(chrom="c", start=0, end=20, strand="+", summit=10, signal=100.0,
                    baseline=1.0, fold=100.0, region_scale=50.0, depth_ratio=dr)
    scores = _scores([mk(1.0), mk(4.0), mk(float("inf"))])
    assert all(np.isfinite(s) and s >= 0 for s in scores)
    assert scores[2] >= scores[1] >= scores[0]


def test_cli_context_window_defaults_to_1000():
    args = build_parser().parse_args(["--bam", "x.bam", "-o", "out.narrowPeak"])
    assert args.context_window == 1000
    args = build_parser().parse_args(
        ["--bam", "x.bam", "-o", "out.narrowPeak", "--context-window", "0"])
    assert args.context_window == 0


# --- --summit-margins (stage-2 tiered summit floor) --------------------------

def test_background_stats_median_and_upper_trimmed_mean():
    from rbpc.calibrate import background_stats
    depths = np.array([1, 1, 1, 2, 2, 3, 100], dtype=np.int64)
    stats = background_stats(depths, trim_pct=(1 / 7) * 100)  # drop just the 100
    assert stats["n"] == 7
    assert stats["median"] == 2.0
    assert stats["trim"] == pytest.approx(np.mean([1, 1, 1, 2, 2, 3]))


def test_background_stats_empty():
    from rbpc.calibrate import background_stats
    stats = background_stats(np.array([], dtype=np.int64))
    assert stats == {"median": 0.0, "trim": 0.0, "n": 0}


def test_tier_for_signal_picks_highest_cleared_tier():
    from rbpc.calibrate import tier_for_signal
    thresholds = [5.0, 10.0, 20.0]
    assert tier_for_signal(5.0, thresholds) == 1
    assert tier_for_signal(9.9, thresholds) == 1
    assert tier_for_signal(10.0, thresholds) == 2
    assert tier_for_signal(19.9, thresholds) == 2
    assert tier_for_signal(20.0, thresholds) == 3
    assert tier_for_signal(1000.0, thresholds) == 3


def test_cli_summit_margins_default_is_off():
    args = build_parser().parse_args(["--bam", "x.bam", "-o", "out.narrowPeak"])
    assert args.summit_margins is None
    assert args.split_tiers is False


def test_cli_summit_margins_parses_ascending_ints():
    from rbpc.cli import _validate
    parser = build_parser()
    args = parser.parse_args(
        ["--bam", "x.bam", "-o", "out.narrowPeak", "--summit-margins", "3,8,18"])
    _validate(parser, args)
    assert args.summit_margins == (3, 8, 18)


def test_cli_summit_margins_single_value_ok():
    from rbpc.cli import _validate
    parser = build_parser()
    args = parser.parse_args(
        ["--bam", "x.bam", "-o", "out.narrowPeak", "--summit-margins", "5"])
    _validate(parser, args)
    assert args.summit_margins == (5,)


@pytest.mark.parametrize("bad", ["8,3", "3,3,8", "-1,3,8", "1,2,3,4", "a,b,c", ""])
def test_cli_summit_margins_rejects_invalid(bad):
    # argparse itself rejects a leading "-1,..." at parse time (looks like an
    # unknown option), everything else is rejected by _validate -- either way
    # it must be a SystemExit, so both calls are wrapped together.
    from rbpc.cli import _validate
    parser = build_parser()
    with pytest.raises(SystemExit):
        args = parser.parse_args(
            ["--bam", "x.bam", "-o", "out.narrowPeak", "--summit-margins", bad])
        _validate(parser, args)


def test_cli_split_tiers_requires_summit_margins():
    from rbpc.cli import _validate
    parser = build_parser()
    args = parser.parse_args(
        ["--bam", "x.bam", "-o", "out.narrowPeak", "--split-tiers"])
    with pytest.raises(SystemExit):
        _validate(parser, args)
