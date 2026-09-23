"""Output formatters: narrowPeak (ENCODE), BED6, and bedGraph."""
from __future__ import annotations

from typing import Iterable, List, TextIO, Tuple

import numpy as np

from .peaks import Peak

# Percentile (within this run's called peaks) used as each component's
# reference scale -- see _scores() for why this replaced fixed constants.
_SCORE_REF_PCT = 95.0


def _raw_components(peak: Peak) -> Tuple[float, float, float, float]:
    prom = max(peak.signal - peak.baseline, 0.0)
    width = max(peak.end - peak.start, 1)
    rel = prom / peak.region_scale if peak.region_scale > 0 else 1.0
    steep = prom / width
    ctx = float(getattr(peak, "depth_ratio", 1.0))
    return prom, rel, steep, ctx


def _scores(peaks: List[Peak]) -> List[int]:
    """Composite 0-1000 score: prominence + region significance + steepness +
    local background contrast, weighted 0.40 / 0.15 / 0.25 / 0.20.

    Each component is scaled by its own magnitude (log for prominence, linear
    for the three ratios) against a *reference* value -- the _SCORE_REF_PCT-th
    percentile of that component across all peaks called in this run, not a
    fixed constant. Fixed constants (previously 3000 / 1.0 / 3.0, tuned on a
    small subset) saturate almost every peak to 1000 on real full-genome data
    where typical steepness/region-relative values run far higher; deriving
    the reference from the run's own peak population keeps the score a
    magnitude-based (not rank-based) but well-spread 0-1000 range regardless
    of the library's absolute depth scale.

    Where the fourth component's weight comes from
    ----------------------------------------------
    ``depth_ratio`` (see context.py) and ``rel`` are two attempts at the *same*
    question -- how far does this peak stand above its surroundings -- but
    ``rel`` answers it self-referentially, because ``region_scale`` is measured
    inside a region whose boundary the hump itself drew, which is exactly the
    "baseline is approximately 0" trap. ``depth_ratio`` measures it against a
    fixed-width window with all called peaks masked out, so it cannot be gamed
    that way. The new weight is therefore taken from ``rel`` (0.35 -> 0.15)
    rather than from ``prominence`` or ``steepness``, whose weights came out of
    real-data tuning and are left untouched.

    The 0.20 share is deliberately not larger. Against motif ground truth
    ``depth_ratio`` outranks every other component measured, including
    prominence, and it is the only one that still discriminates among weak
    peaks -- but that evidence comes from a single proxy label on three
    libraries, and the confirmed CDH2 cluster sits only mid-distribution on it.
    A minority share improves ranking without letting one proxy-validated
    metric dominate a score whose other weights were set on inspected data.
    """
    n = len(peaks)
    if n == 0:
        return []
    prom = np.empty(n)
    rel = np.empty(n)
    steep = np.empty(n)
    ctx = np.empty(n)
    for i, p in enumerate(peaks):
        prom[i], rel[i], steep[i], ctx[i] = _raw_components(p)

    # A peak with no usable local background (whole window masked, or no reads
    # outside it) has an infinite ratio; treat it as the strongest observed
    # rather than letting it poison the percentile reference.
    finite = np.isfinite(ctx)
    ctx = np.where(finite, ctx, float(ctx[finite].max()) if finite.any() else 1.0)

    def ref(x: np.ndarray) -> float:
        return max(float(np.percentile(x, _SCORE_REF_PCT)), 1e-9)

    # No 1000 ceiling: a peak sitting exactly at the _SCORE_REF_PCT reference
    # on all four axes scores ~1000, but a peak well beyond that reference
    # (the strongest of the strong) is allowed to score higher, rather than
    # collapsing the whole top tail onto a single flat ceiling value.
    p_comp = np.log1p(prom) / np.log1p(ref(prom))
    r_comp = rel / ref(rel)
    s_comp = steep / ref(steep)
    c_comp = ctx / ref(ctx)
    combined = 0.40 * p_comp + 0.15 * r_comp + 0.25 * s_comp + 0.20 * c_comp
    scores = np.maximum(0, np.round(combined * 1000)).astype(int)
    return scores.tolist()


def _write_header(handle: TextIO, header_lines) -> None:
    if not header_lines:
        return
    for line in header_lines:
        handle.write(f"# {line}\n")


def write_narrowpeak(peaks: Iterable[Peak], handle: TextIO, *, header_lines=None) -> int:
    peaks = list(peaks)
    _write_header(handle, header_lines)
    scores = _scores(peaks)
    for p, score in zip(peaks, scores):
        handle.write(
            f"{p.chrom}\t{p.start}\t{p.end}\t{p.name}\t{score}\t{p.strand}\t"
            f"{(p.signal - p.baseline):.5f}\t-1\t-1\t{p.summit_offset}\n"
        )
    return len(peaks)


def write_bed6(peaks: Iterable[Peak], handle: TextIO, *, header_lines=None) -> int:
    peaks = list(peaks)
    _write_header(handle, header_lines)
    scores = _scores(peaks)
    for p, score in zip(peaks, scores):
        handle.write(f"{p.chrom}\t{p.start}\t{p.end}\t{p.name}\t{score}\t{p.strand}\n")
    return len(peaks)


def write_bedgraph(chrom, signal, handle, *, strand=".", offset=0):
    n = signal.size
    if n == 0:
        return
    change = np.flatnonzero(np.diff(signal)) + 1
    starts = np.concatenate(([0], change))
    ends = np.concatenate((change, [n]))
    vals = signal[starts]
    nz = vals != 0
    for s, e, v in zip(starts[nz], ends[nz], vals[nz]):
        handle.write(f"{chrom}\t{int(s) + offset}\t{int(e) + offset}\t{float(v):.5f}\n")
