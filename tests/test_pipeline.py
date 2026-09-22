import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))
import make_example_bam  # noqa: E402

from peakcaller import coverage  # noqa: E402
from peakcaller.pipeline import Config, run  # noqa: E402


@pytest.fixture(scope="module")
def bam(tmp_path_factory):
    path = str(tmp_path_factory.mktemp("data") / "example.bam")
    make_example_bam.build(path, seed=1)
    return path


def test_splice_aware_coverage(bam):
    cov = coverage.build_coverage(bam, library="forward")
    plus = cov.chroms["chr_test"]["+"]
    assert plus[9005] >= 5
    assert plus[11015] >= 5
    assert plus[10000] < plus[9005]        # intron not filled


def test_prominence_detects_mountains_not_weak_bump(bam):
    cfg = Config(normalize_method="none", min_prominence=20, min_summit_reads=5,
                 region_gap=200, min_distance=15, min_peak_width=5)
    res = run(bam, cfg)
    summits = np.array([p.summit for p in res.peaks if p.strand == "+"])
    assert np.any(np.abs(summits - 5000) < 120)        # mountain 1
    assert np.any(np.abs(summits - 5300) < 120)        # mountain 2 (split out)
    assert not np.any(np.abs(summits - 12000) < 250)   # weak bump rejected


def test_flank_widens_symmetrically_and_clamps(bam):
    base = run(bam, Config(normalize_method="none", min_prominence=20, min_summit_reads=5,
                            flank=0))
    flanked = run(bam, Config(normalize_method="none", min_prominence=20, min_summit_reads=5,
                               flank=20))
    base_by_summit = {p.summit: p for p in base.peaks}
    for fp in flanked.peaks:
        bp = base_by_summit[fp.summit]
        assert fp.start == max(0, bp.start - 20)
        assert fp.end == bp.end + 20  # test peaks are far from chrom end, so no clamping
        assert fp.summit_offset == fp.summit - fp.start


def test_region_min_depth_default_is_still_explicit_not_auto():
    # Config() must stay backward-compatible for library users who don't opt in.
    assert Config().region_min_depth == 3.0


def test_region_min_depth_none_triggers_auto_computation(bam):
    cfg = Config(normalize_method="none", min_prominence=20, min_summit_reads=5,
                 region_min_depth=None)
    res = run(bam, cfg)
    # auto mode must still find the same two mountains as the explicit default.
    summits = np.array([p.summit for p in res.peaks if p.strand == "+"])
    assert np.any(np.abs(summits - 5000) < 120)
    assert np.any(np.abs(summits - 5300) < 120)


def test_auto_region_min_depth_rejects_self_referential_noise(bam):
    from peakcaller import calibrate
    value, chroms = calibrate.auto_region_min_depth(bam, library="forward")
    assert value >= 1.0
    assert chroms


def test_min_region_width_default_is_still_off_not_auto():
    # Config() must stay backward-compatible for library users who don't opt in;
    # the CLI is what defaults to auto (see test_units).
    assert Config().min_region_width == 0


def test_min_region_width_none_triggers_auto_computation(bam):
    cfg = Config(normalize_method="none", min_prominence=20, min_summit_reads=5,
                 min_region_width=None)
    res = run(bam, cfg)
    # auto mode must still find the same two mountains as the explicit default:
    # they sit in regions far wider than one read, so the width floor can't touch them.
    summits = np.array([p.summit for p in res.peaks if p.strand == "+"])
    assert np.any(np.abs(summits - 5000) < 120)
    assert np.any(np.abs(summits - 5300) < 120)


def test_auto_min_region_width_is_two_read_lengths(bam):
    from peakcaller import calibrate
    width, diag = calibrate.auto_min_region_width(bam)
    assert diag["n_reads"] > 0
    assert diag["sampled_chroms"]
    # the whole point: the floor is derived from this BAM's own read lengths,
    # so a region must be wider than any single read stack can ever be.
    assert width == round(diag["factor"] * diag["read_length"])
    assert width > diag["read_length"]


def test_auto_min_region_width_scales_with_factor(bam):
    from peakcaller import calibrate
    w2, _ = calibrate.auto_min_region_width(bam, factor=2.0)
    w4, _ = calibrate.auto_min_region_width(bam, factor=4.0)
    assert w4 == 2 * w2


@pytest.fixture(scope="module")
def stack_bam(tmp_path_factory):
    """A near-duplicate read stack (the artifact --min-region-width targets)
    plus one genuinely wide mountain, so the filter's selectivity is testable."""
    import random

    import pysam
    random.seed(7)
    path = str(tmp_path_factory.mktemp("stack") / "stack.bam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr_test", "LN": 20_000}]}
    starts = [int(random.gauss(3000, 3)) for _ in range(60)]          # one clump, ~46bp wide
    starts += [int(random.gauss(10_000, 60)) for _ in range(400)]     # real mountain, ~500bp wide
    with pysam.AlignmentFile(path, "wb", header=header) as bam_out:
        for i, pos in enumerate(sorted(starts)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 60
            a.flag = 0
            a.cigartuples = [(0, 30)]
            a.query_sequence = "A" * 30
            a.query_qualities = pysam.qualitystring_to_array("I" * 30)
            bam_out.write(a)
    pysam.index(path)
    return path


def test_min_region_width_drops_near_duplicate_stack(stack_bam):
    from peakcaller import calibrate
    width, diag = calibrate.auto_min_region_width(stack_bam)
    assert width == 60  # 2 x 30bp reads

    def summits(min_region_width):
        res = run(stack_bam, Config(normalize_method="none", min_prominence=1,
                                    min_summit_reads=2, min_peak_width=1, min_distance=5,
                                    region_min_depth=3.0, min_region_width=min_region_width))
        return [p.summit for p in res.peaks]

    off = summits(0)
    auto = summits(None)          # None = auto -> 60bp
    # with the filter off the one-clump artifact is called; with it on, only the
    # real mountain survives (its region is far wider than any single read).
    assert any(abs(s - 3000) < 60 for s in off)
    assert not any(abs(s - 3000) < 60 for s in auto)
    assert any(abs(s - 10_000) < 150 for s in auto)


def test_min_prominence_default_is_still_explicit_not_auto():
    # Config() must stay backward-compatible for library users who don't opt in.
    assert Config().min_prominence == 5.0


def test_min_prominence_none_triggers_auto_computation(bam):
    cfg = Config(normalize_method="none", min_summit_reads=5, min_prominence=None,
                 fdr_n_shuffles=1)
    res = run(bam, cfg)
    # auto mode must still find the same two mountains as the explicit default.
    summits = np.array([p.summit for p in res.peaks if p.strand == "+"])
    assert np.any(np.abs(summits - 5000) < 120)
    assert np.any(np.abs(summits - 5300) < 120)


def test_auto_min_prominence_rejects_noisier_null_than_real(bam):
    from peakcaller import calibrate
    value, diag = calibrate.auto_min_prominence(bam, library="forward", n_shuffles=1)
    assert value >= 5.0
    assert diag["sampled_chroms"]


def test_narrowpeak_columns(bam, tmp_path):
    from peakcaller.output import write_narrowpeak
    res = run(bam, Config(normalize_method="none", min_prominence=20, min_summit_reads=5))
    out = tmp_path / "peaks.narrowPeak"
    with open(out, "w") as fh:
        write_narrowpeak(res.peaks, fh)
    for line in [l for l in out.read_text().splitlines() if l.strip()]:
        cols = line.split("\t")
        assert len(cols) == 10
        assert int(cols[1]) < int(cols[2])
        assert cols[5] in ("+", "-", ".")


@pytest.fixture(scope="module")
def dup_stack_bam(tmp_path_factory):
    """The artifact --min-complexity targets, which --min-region-width cannot see:
    a *wide* region (so it survives the stage-1 width floor) that is nevertheless
    built from only three heavily amplified molecules, next to a genuinely diverse
    mountain of the same shape."""
    import random

    import pysam
    random.seed(11)
    path = str(tmp_path_factory.mktemp("dup") / "dup_stack.bam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr_test", "LN": 20_000}]}
    # 3 distinct start positions, 100 PCR copies each, spread over 130bp:
    # the region is far wider than one read, so stage 1 keeps it, and the middle
    # stack is a clean, steep hump -- but it is 3 molecules, not 300.
    starts = [3000] * 20 + [3050] * 100 + [3100] * 20
    # a real mountain: 400 molecules, each start position its own
    starts += [int(random.gauss(10_000, 60)) for _ in range(400)]
    with pysam.AlignmentFile(path, "wb", header=header) as bam_out:
        for i, pos in enumerate(sorted(starts)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 60
            a.flag = 0
            a.cigartuples = [(0, 30)]
            a.query_sequence = "A" * 30
            a.query_qualities = pysam.qualitystring_to_array("I" * 30)
            bam_out.write(a)
    pysam.index(path)
    return path


def _dup_stack_summits(bam_path, **over):
    cfg = dict(normalize_method="none", min_prominence=1, min_summit_reads=2,
               min_peak_width=1, min_distance=5, min_steepness=0.5,
               region_min_depth=3.0, min_region_width=None)
    cfg.update(over)
    return [p.summit for p in run(bam_path, Config(**cfg)).peaks]


def test_min_complexity_drops_amplified_stack_that_width_filter_cannot_see(dup_stack_bam):
    off = _dup_stack_summits(dup_stack_bam, min_complexity=0.0)
    on = _dup_stack_summits(dup_stack_bam)          # None = auto -> min_steepness (0.5)
    # the amplified stack is tall, steep and in a wide region, so every other
    # gate lets it through...
    assert any(2990 <= s <= 3140 for s in off)
    # ...and only positional diversity rejects it, while the real mountain --
    # same shape, but built from hundreds of distinct molecules -- survives.
    assert not any(2990 <= s <= 3140 for s in on)
    assert any(abs(s - 10_000) < 150 for s in on)


def test_min_complexity_is_off_when_no_steepness_is_asked_for(dup_stack_bam):
    # backward compatibility: Config() defaults (min_steepness=0) leave it off.
    assert Config().min_complexity is None
    assert any(2990 <= s <= 3140
               for s in _dup_stack_summits(dup_stack_bam, min_steepness=0.0))


def test_min_complexity_explicit_value_overrides_the_steepness_tie(dup_stack_bam):
    # a floor above the real mountain's own diversity drops everything;
    # this is the knob's stringency direction, and it is independent of steepness.
    assert _dup_stack_summits(dup_stack_bam, min_complexity=100.0) == []


# --- local background contrast (depth_ratio) --------------------------------

@pytest.fixture(scope="module")
def context_bam(tmp_path_factory):
    """Two identically shaped humps: one isolated on empty sequence, one sitting
    on a broad carpet of background coverage. The stage-2 region test cannot
    tell them apart -- the isolated one becomes its own region, so its
    region_scale is its own height -- but the local-background ratio can."""
    import random

    import pysam
    random.seed(7)
    path = str(tmp_path_factory.mktemp("ctx") / "context.bam")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr_test", "LN": 30_000}]}
    starts = [int(random.gauss(5_000, 15)) for _ in range(300)]     # isolated hump
    starts += [int(random.gauss(20_000, 15)) for _ in range(300)]   # same hump...
    starts += list(range(17_000, 23_000, 8))                        # ...on a carpet
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for i, pos in enumerate(sorted(starts)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = max(0, pos)
            a.mapping_quality = 60
            a.flag = 0
            a.cigartuples = [(0, 30)]
            a.query_sequence = "A" * 30
            a.query_qualities = pysam.qualitystring_to_array("I" * 30)
            out.write(a)
    pysam.index(path)
    return path


def _ctx_cfg(**over):
    cfg = dict(normalize_method="none", min_prominence=5, min_summit_reads=5,
               min_peak_width=1, min_distance=20, min_steepness=0.0,
               prominence_frac=0.2, region_min_depth=3.0, min_region_width=0,
               min_complexity=0.0, library="forward")
    cfg.update(over)
    return Config(**cfg)


def _nearest(peaks, pos):
    return min(peaks, key=lambda p: abs(p.summit - pos))


def test_depth_ratio_separates_isolated_hump_from_one_on_a_carpet(context_bam):
    peaks = run(context_bam, _ctx_cfg()).peaks
    isolated = _nearest(peaks, 5_000)
    on_carpet = _nearest(peaks, 20_000)
    # The point of the metric: same hump, different neighbourhood.
    assert isolated.depth_ratio > on_carpet.depth_ratio
    # ...and it is a real difference in the data, not an artifact of the hump:
    # the two humps are the same size.
    assert isolated.peak_depth == pytest.approx(on_carpet.peak_depth, rel=0.2)
    assert isolated.bg_depth < on_carpet.bg_depth
    # The existing region-relative test cannot make this distinction, because
    # the isolated hump *is* its own region: region_scale is measured from the
    # hump itself, so prominence >= prominence_frac * region_scale reduces to
    # comparing the hump to a fraction of itself and is satisfied trivially.
    iso_prom = isolated.signal - isolated.baseline
    assert isolated.region_scale >= 0.5 * iso_prom
    assert iso_prom >= 0.2 * isolated.region_scale


def test_context_window_does_not_change_which_peaks_are_called(context_bam):
    """depth_ratio is a scoring component only: it must never gate a peak."""
    with_ctx = run(context_bam, _ctx_cfg()).peaks
    without = run(context_bam, _ctx_cfg(context_window=0)).peaks
    assert [(p.chrom, p.start, p.end, p.summit) for p in with_ctx] == \
           [(p.chrom, p.start, p.end, p.summit) for p in without]
    assert [p.signal - p.baseline for p in with_ctx] == \
           [p.signal - p.baseline for p in without]


def test_context_window_zero_leaves_depth_ratio_neutral(context_bam):
    peaks = run(context_bam, _ctx_cfg(context_window=0)).peaks
    assert peaks and all(p.depth_ratio == 1.0 for p in peaks)


def test_depth_ratio_is_measured_on_the_flanked_interval(context_bam):
    """--flank widens the written interval, so the measurement must follow it."""
    a = _nearest(run(context_bam, _ctx_cfg(flank=0)).peaks, 5_000)
    b = _nearest(run(context_bam, _ctx_cfg(flank=200)).peaks, 5_000)
    assert b.end - b.start > a.end - a.start
    assert b.depth_ratio != a.depth_ratio
