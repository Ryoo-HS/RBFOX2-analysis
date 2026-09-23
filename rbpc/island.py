"""Stage 1: island detection."""
from __future__ import annotations

import numpy as np


def find_islands(signal, *, region_gap=200, min_signal=0.0, min_region_len=1):
    # depth 0 is never "covered", even when min_signal == 0 (no threshold requested) --
    # otherwise every true zero-coverage gap would count as covered and swallow region_gap.
    covered = (signal > 0) & (signal >= min_signal)
    if not covered.any():
        return []
    idx = np.flatnonzero(covered)
    breaks = np.flatnonzero(np.diff(idx) > 1)
    run_starts = np.concatenate(([idx[0]], idx[breaks + 1]))
    run_ends = np.concatenate((idx[breaks] + 1, [idx[-1] + 1]))
    islands = [[int(run_starts[0]), int(run_ends[0])]]
    for s, e in zip(run_starts[1:], run_ends[1:]):
        if s - islands[-1][1] <= region_gap:
            islands[-1][1] = int(e)
        else:
            islands.append([int(s), int(e)])
    return [(s, e) for s, e in islands if (e - s) >= min_region_len]
