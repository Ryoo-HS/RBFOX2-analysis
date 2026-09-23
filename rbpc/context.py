"""Local background contrast (``depth_ratio``) -- a scoring component.

For each called peak, compare its mean depth to the mean depth of the
*non-peak* sequence in a fixed-width window around it, on the same strand:

    depth_ratio = mean_depth(peak) / mean_depth(window minus all called peaks)

``window`` is ``[start - context_window, end + context_window)`` (default
1000bp each side). Every position that belongs to *any* called peak on that
strand is removed from the denominator, so a cluster of genuine adjacent sites
does not suppress its own members.

Why this is not another region-relative test
--------------------------------------------
``core.call_cores`` already has a region-relative test, ``prominence >=
prominence_frac * region_scale``, where ``region_scale`` is a percentile of the
covered depth *inside the peak's own stage-1 region*. That test is
**self-referential** for an isolated hump: the region boundary is drawn by the
data, so a hump standing alone becomes its own region, ``region_scale`` is
computed from the hump itself, and the test is satisfied no matter what. This
is the project's long-standing "baseline is approximately 0" trap.

``depth_ratio`` closes that trap structurally rather than by tuning: the
denominator is a **fixed physical window**, not a data-drawn region, and the
peak's own positions (and those of every other called peak) are excluded from
it. The denominator therefore *cannot* be the hump itself.

It is also a pure ratio of depths, so it is invariant to sequencing depth and
to ``--normalize-method`` -- the normalization factor cancels. That is the
property an absolute magnitude floor lacks, and the reason a magnitude floor
was rejected twice: raw depth counts PCR copies, not molecules, and does not
transfer between libraries.

What the measurement showed (3 eCLIP libraries, motif ground truth)
-------------------------------------------------------------------
Validated against an independent external label -- presence of each protein's
canonical motif within 25bp of the summit (RBFOX2 ``TGCATG``, QKI ``ACTAAC``,
ESRP1 ``GGTGGT``/``GGTGAT``/``GATGGT``), scored as AUROC minus the AUROC of the
same metric against shifted control windows, which subtracts off any sequence
or expression-level confounding:

  | metric                  | RBFOX2 | QKI    | ESRP1  |
  |-------------------------|--------|--------|--------|
  | depth_ratio             | +0.278 | +0.245 | +0.085 |
  | prominence (magnitude)  | +0.094 | +0.089 | -0.028 |

Restricted to the bottom 40% of peaks by prominence -- the weak isolated humps
that the relative test cannot reject and that ``min_complexity`` passes -- the
gap widens into a qualitative difference:

  | metric                  | RBFOX2 | QKI    | ESRP1  |
  |-------------------------|--------|--------|--------|
  | depth_ratio             | +0.213 | +0.217 | +0.108 |
  | prominence              | -0.002 | +0.033 | -0.022 |
  | distinct molecules      | -0.007 | -0.062 | +0.120 |

So among weak peaks, magnitude carries essentially no information and
``depth_ratio`` still does. Motif rate rises monotonically across its
quartiles (RBFOX2 3.3% -> 16.5%, against a 1.4% background rate).

Why it is a score component and not a filter
---------------------------------------------
No threshold on this metric is defensible yet. Its distribution decays
smoothly with no elbow -- the same result the depth, region-width,
``width_complexity`` and steepness distributions all gave -- so any cutoff
would be a stringency dial, not a measurement. More decisively, the project's
one externally confirmed anchor argues against a hard gate: the CDH2 cluster
(ESRP1, chr18:25.53Mb, 8 real peaks) sits at percentiles 11, 38, 45, 50, 57,
61, 66 and 67 of ESRP1's own ``depth_ratio`` distribution, because CDH2 is a
highly expressed locus with a correspondingly high local background. A cut at
the median would delete half of a cluster known to be real. Ranking is
therefore the honest use of this signal; gating is not.

Rejected alternatives
---------------------
  * Unmasked ratio (peak reads / all flanking reads). Clustered real sites
    suppress each other: each CDH2 peak's background is inflated by its 8
    neighbours, dropping the whole cluster's ratio from 2.8-12.3 to 2.1-8.7.
    Masking called peaks out of the denominator is what makes the measurement
    fair to clusters.
  * Widening ``region_gap``/``coverage_gap`` so the existing ``region_scale``
    sees more context. This does not fix anything: the boundary would still be
    drawn by the data, so it stays self-referential, and it would also change
    which peaks are detected at all.
  * Flanking coverage as a positive indicator ("a real site sits in an
    expressed transcript"). Measured, and the data says the opposite -- peaks
    in the most densely covered neighbourhoods have the *lowest* motif rate
    (delta AUROC -0.15 / -0.17 / -0.11). Splice-junction support in the window
    gave AUROC exactly 0.500, i.e. no information at all.
  * Strand-aware 5'-end diversity, and 5'-end concentration as a crosslink
    truncation signature. Both measured and rejected: the 5'-aware variant of
    ``width_complexity`` scored *worse* than the existing genomic-left version
    (+0.011 vs +0.023), and 5'-end concentration was flat (+0.003) with its
    sign inverted from the hypothesis -- concentrated 5' ends indicate PCR
    duplication, not crosslinking.
"""
from __future__ import annotations

import numpy as np

DEFAULT_CONTEXT_WINDOW = 1000


def build_block_index(blocks) -> tuple[np.ndarray, np.ndarray, int]:
    """Sort one strand's read blocks into ``(starts, ends, max_len)``."""
    bl = list(blocks)
    if not bl:
        return (np.empty(0, np.int64), np.empty(0, np.int64), 0)
    arr = np.asarray(bl, dtype=np.int64)
    order = np.argsort(arr[:, 0], kind="stable")
    starts = arr[order, 0].copy()
    ends = arr[order, 1].copy()
    max_len = int((ends - starts).max()) if starts.size else 0
    return starts, ends, max_len


def window_depth(starts: np.ndarray, ends: np.ndarray, max_len: int,
                 lo: int, hi: int) -> np.ndarray:
    """Per-bp raw depth over ``[lo, hi)`` from sorted read blocks."""
    n = hi - lo
    if n <= 0 or starts.size == 0:
        return np.zeros(max(n, 0), dtype=np.int32)
    i0 = int(np.searchsorted(starts, lo - max_len, side="left"))
    i1 = int(np.searchsorted(starts, hi, side="left"))
    if i1 <= i0:
        return np.zeros(n, dtype=np.int32)
    # Vectorized +1/-1 sweep rather than a per-block slice loop: this runs once
    # per called peak over every read block in the window, so it is the hot path.
    a = np.clip(starts[i0:i1] - lo, 0, n)
    b = np.clip(ends[i0:i1] - lo, 0, n)
    diff = np.zeros(n + 1, dtype=np.int32)
    np.add.at(diff, a, 1)
    np.add.at(diff, b, -1)
    return np.cumsum(diff[:n], dtype=np.int32)


def depth_ratio(starts: np.ndarray, ends: np.ndarray, max_len: int,
                peak_start: int, peak_end: int,
                peak_bounds: tuple[np.ndarray, np.ndarray],
                window: int = DEFAULT_CONTEXT_WINDOW,
                ) -> tuple[float, float, float]:
    """Return ``(depth_ratio, peak_mean_depth, background_mean_depth)``.

    ``peak_bounds`` is ``(starts, ends)`` of every called peak on this strand
    and chromosome, sorted by start; all of them are masked out of the
    background. Returns ``inf`` for the ratio when no usable background
    remains (the caller replaces it with the run's maximum finite value).
    """
    lo = max(0, int(peak_start) - window)
    hi = int(peak_end) + window
    arr = window_depth(starts, ends, max_len, lo, hi)
    if arr.size == 0:
        return float("inf"), 0.0, 0.0
    ps, pe = peak_bounds
    usable = np.ones(arr.size, dtype=bool)
    if ps.size:
        # Peaks are narrow relative to the window; bound the scan by the widest
        # one so no overlapping peak is missed regardless of how many there are.
        span = int((pe - ps).max())
        j0 = int(np.searchsorted(ps, lo - span, side="left"))
        j1 = int(np.searchsorted(ps, hi, side="left"))
        if j1 > j0:
            a = np.clip(ps[j0:j1] - lo, 0, arr.size)
            b = np.clip(pe[j0:j1] - lo, 0, arr.size)
            m = np.zeros(arr.size + 1, dtype=np.int32)
            np.add.at(m, a, 1)
            np.add.at(m, b, -1)
            usable &= np.cumsum(m[:arr.size]) == 0
    a0 = max(0, int(peak_start) - lo)
    b0 = min(arr.size, int(peak_end) - lo)
    usable[a0:b0] = False
    pk_depth = float(arr[a0:b0].mean()) if b0 > a0 else 0.0
    if not usable.any():
        return float("inf"), pk_depth, 0.0
    bg_depth = float(arr[usable].mean())
    if bg_depth <= 0:
        return float("inf"), pk_depth, 0.0
    return pk_depth / bg_depth, pk_depth, bg_depth


def annotate(peaks, blocks_by_strand, window: int = DEFAULT_CONTEXT_WINDOW) -> None:
    """Attach ``depth_ratio`` to every peak of one chromosome, in place.

    ``blocks_by_strand`` maps strand -> list of raw ``(start, end)`` read
    blocks for this chromosome. Peaks must already carry their final
    coordinates (i.e. ``--flank`` applied), so the measurement matches the
    intervals that are written out.
    """
    if not peaks:
        return
    by_strand: dict[str, list] = {}
    for p in peaks:
        by_strand.setdefault(p.strand, []).append(p)
    for strand, plist in by_strand.items():
        starts, ends, max_len = build_block_index(blocks_by_strand.get(strand, ()))
        order = sorted(range(len(plist)), key=lambda i: plist[i].start)
        ps = np.array([plist[i].start for i in order], dtype=np.int64)
        pe = np.array([plist[i].end for i in order], dtype=np.int64)
        for p in plist:
            r, pk, bg = depth_ratio(starts, ends, max_len, p.start, p.end,
                                    (ps, pe), window=window)
            p.depth_ratio = r
            p.peak_depth = pk
            p.bg_depth = bg
