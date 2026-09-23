"""Stage-2 post-filter: read positional diversity (``width_complexity``).

For each called peak, count how many *distinct genomic start positions* the
reads overlapping it have, and divide by the peak's width:

    width_complexity = n_distinct_read_starts / (peak_end - peak_start)

For single-end data the alignment start position is the standard duplicate
signature (it is what ``samtools markdup`` and MACS's ``--keep-dup`` collapse
on), so ``n_distinct_read_starts`` is the number of *independent molecules*
that produced the peak, and ``width_complexity`` is "independent molecules per
bp of peak".

Why this and not a taller prominence floor
------------------------------------------
An isolated hump sits on a zero local baseline, so its prominence equals its
own full height and it satisfies ``prominence >= prominence_frac *
region_scale`` trivially -- ``region_scale`` is that same hump. The relative
test is scale-free and therefore cannot reject it, and the stage-1
``min_region_width`` filter only catches stacks narrower than a read, not a
handful of amplified stacks spread over a few hundred bp. Raising an absolute
``min_prominence`` was tried and rejected twice: a magnitude floor does not
transfer across libraries of different sequencing depth, because depth counts
PCR copies, not molecules.

Positional diversity is a different *kind* of measurement: it asks how many
independent molecules contributed, not how tall the pile is. Measured on real
data (ESRP1, whole genome, default settings), the two answers actively
disagree, in both directions:

  * dropped: ``chr4:113,569,643`` prominence **1631** from **7** distinct
    starts -- 1631 reads / 7 molecules = 233 PCR copies each.
  * kept: peaks of prominence **5** built from 5-11 distinct molecules.

A magnitude floor keeps the first and drops the second; this filter does the
opposite, which is the intended behaviour.

Where the threshold comes from (not a tuned number)
---------------------------------------------------
``width_complexity`` is an **upper bound on the peak's own steepness measured
on deduplicated coverage**. Proof: build the coverage that would result from
keeping one read per distinct start ("one read per molecule"). Every position
inside the peak is then covered by at most ``n_distinct_starts`` reads, so the
deduplicated depth -- and therefore the deduplicated prominence, which is depth
minus a non-negative valley -- is at most ``n_distinct_starts``. Dividing by
the same width gives

    steepness_deduplicated <= n_distinct_starts / width = width_complexity

So ``width_complexity < min_steepness`` proves the peak **cannot** meet the
caller's own ``--min-steepness`` bar once PCR duplicates are collapsed: it is
steep only because some molecules were sequenced many times. That is why the
default threshold is ``--min-steepness`` itself rather than a new constant --
the same bar the user already set, applied to molecules instead of reads.

Verified on real data (ESRP1 chr18, 1,106 peaks): the bound holds for 100% of
peaks, and it is *tight* -- the median ratio of measured deduplicated steepness
to ``width_complexity`` is 1.00, because for these narrow peaks essentially all
distinct-start reads cover the summit and the flanking valleys are ~0. So the
filter is not merely conservative in principle; in practice it agrees with the
deduplicated measurement it bounds.

Rejected alternatives
---------------------
  * A percentile of the observed ``width_complexity`` distribution. Rejected
    for the same reason the depth and region-width percentiles were: the
    distribution decays smoothly with no elbow (ESRP1 whole genome: 83% / 57% /
    38% / 20% / 5.5% of peaks survive at 0.1 / 0.2 / 0.3 / 0.5 / 1.0), so any
    percentile is a stringency dial dressed up as a measurement.
  * A "geometric impossibility" anchor like the one behind
    ``min_region_width`` -- e.g. "all read starts fall within one read length,
    so this could be one stack". Tested and rejected: because peaks are narrow
    (median width 20bp), that fraction decays smoothly and never reaches 0
    (10% even at ``width_complexity`` > 3), so there is no clean cutoff.
  * Duplication *rate* (reads / distinct starts). Rejected on measurement: a
    genuinely strong site is also highly amplified (the CDH2 cluster runs 16
    reads per molecule), so this flags real signal as strongly as artifacts.
  * Actually rebuilding deduplicated coverage and re-running detection on it.
    This is the exact version of the test that ``width_complexity`` bounds, and
    it works, but it amounts to running a second, different caller inside the
    first (with its own region gate, so its own new thresholds), and it changes
    what ``core.call_cores`` sees. The bound is a necessary condition computed
    from the peaks the caller actually called, and -- per the measurement above
    -- is tight, so it buys the same selectivity without a second detector.

What this filter does and does not do (measured, whole genome)
--------------------------------------------------------------
It is a *duplication* filter, so how much it removes depends on how duplicated
the library is -- it is not a universal peak-count dial:

  | library | median copies/molecule | peaks | at min-steepness 0.5 |
  |---------|------------------------|-------|----------------------|
  | ESRP1   | 6.5                    |  73,678 | 15,042  (-80%)     |
  | RBFOX2  | 8.0                    | 415,711 | 48,187  (-88%)     |
  | QKI     | 11.0                   |  78,523 |  6,889  (-91%)     |
  | YBX1    | 2.4                    | 383,478 | 288,951 (-25%)     |

YBX1 is barely touched because its library genuinely is not duplicated (median
2.4 copies per molecule, median 14 distinct molecules per peak, vs 2-4 distinct
molecules for the others). That is the filter reporting a property of the data,
not failing: where there is little amplification there are few duplicate-driven
humps to remove, and that library's peak count stays high for reasons this
filter is not designed to address.

Sensitivity check against the project's known anchor: the CDH2 cluster (ESRP1,
chr18:25.53Mb, previously confirmed real) keeps all 8 of its peaks with
prominences 64-467 unchanged -- their width_complexity runs 0.889-3.222, clear
of the 0.5 bar, while the weak peaks immediately around them sit at 0.079-0.5.
"""
from __future__ import annotations

from typing import Iterable

import numpy as np


def resolve_min_complexity(min_complexity: float | None, min_steepness: float,
                           factor: float = 1.0) -> float:
    """Return the ``width_complexity`` floor to apply (0 = filter off).

    ``None`` means "use ``min_steepness``", the threshold justified in the
    module notes above: ``width_complexity`` bounds the peak's steepness on
    deduplicated coverage, so requiring the two to be equal asks the peak to
    still clear the steepness bar the user already set when every molecule is
    counted once.

    ``factor`` is the normalization scaling factor. ``min_steepness`` is in
    signal units per bp, while ``width_complexity`` is in *molecules* per bp;
    under ``--normalize-method rpm`` one molecule is worth ``factor`` signal
    units, so the comparable molecule-per-bp floor is ``min_steepness /
    factor``. With the default ``none`` normalization ``factor`` is 1.0 and
    this is a no-op.
    """
    if min_complexity is not None:
        return max(float(min_complexity), 0.0)
    if factor <= 0:
        return max(float(min_steepness), 0.0)
    return max(float(min_steepness) / factor, 0.0)


def block_bounds(blocks: Iterable[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray]:
    """Split an island's ``[(start, end), ...]`` read blocks into two arrays."""
    bl = list(blocks)
    if not bl:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    arr = np.asarray(bl, dtype=np.int64)
    return arr[:, 0], arr[:, 1]


def width_complexity(starts: np.ndarray, ends: np.ndarray, peak_start: int,
                     peak_end: int) -> tuple[float, int]:
    """Return ``(width_complexity, n_distinct_starts)`` for one peak.

    ``starts``/``ends`` are the read blocks of the island containing the peak,
    in absolute genomic coordinates; only blocks overlapping
    ``[peak_start, peak_end)`` are counted.
    """
    width = max(int(peak_end) - int(peak_start), 1)
    if starts.size == 0:
        return 0.0, 0
    overlapping = (ends > peak_start) & (starts < peak_end)
    if not overlapping.any():
        return 0.0, 0
    n_distinct = int(np.unique(starts[overlapping]).size)
    return n_distinct / width, n_distinct


def filter_cores(cores: list[dict], island_start: int, blocks, min_complexity: float,
                 ) -> list[dict]:
    """Drop cores whose read positional diversity is below ``min_complexity``.

    ``cores`` are :func:`core.call_cores` results with island-relative
    ``start``/``end``; ``island_start`` converts them to the absolute
    coordinates the blocks are in. Each surviving core gains ``complexity`` and
    ``n_distinct_starts`` keys. ``min_complexity <= 0`` disables the filter
    (the annotations are then not computed, so this stays free when off).
    """
    if min_complexity <= 0 or not cores:
        return cores
    starts, ends = block_bounds(blocks)
    kept = []
    for c in cores:
        cx, nd = width_complexity(starts, ends, island_start + c["start"],
                                  island_start + c["end"])
        if cx < min_complexity:
            continue
        c["complexity"] = cx
        c["n_distinct_starts"] = nd
        kept.append(c)
    return kept
