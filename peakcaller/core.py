"""Stage 2: prominent local maxima (mountains) via topographic prominence.

A peak is a local maximum whose topographic prominence -- how far it rises above
the higher of the two valleys that flank it -- is at least ``min_prominence``.
Prominence naturally splits a broad multi-hump region into one peak per hump.
Implementation uses :func:`scipy.signal.find_peaks`.
"""
from __future__ import annotations

import logging
import warnings

import numpy as np
from scipy.signal import find_peaks

_EPS = 1e-9
log = logging.getLogger("peakcaller")


def call_cores(
    signal: np.ndarray,
    islands: list[tuple[int, int]],
    *,
    min_prominence: float = 5.0,
    prominence_frac: float = 0.0,
    region_scale_pct: float = 90.0,
    min_steepness: float = 0.0,
    min_summit_reads: float = 1.0,
    min_distance: int = 15,
    min_peak_width: int = 5,
    rel_height: float = 0.3,
    max_peaks_per_region: int = 0,
    **_compat,
) -> list[dict]:
    """Detect prominent local-maxima peaks inside each region.

    min_prominence : absolute floor on how far a summit rises above its flanking
        valleys (signal units).
    prominence_frac : region-adaptive requirement -- prominence >= prominence_frac
        * region_scale (a percentile of the region's covered depth). Drops ripples
        and shoulder bumps on a tall gene; 0 disables it.
    region_scale_pct : percentile of the region's covered depth used as region scale.
    min_steepness : minimum prominence / width (rise per bp); 0 disables it.
    min_summit_reads : absolute height floor at the summit.
    min_distance : minimum bp between two summits.
    min_peak_width : minimum peak width (bp).
    rel_height : where the boundary is measured (0 = tip, 1 = base); lower = tighter.
    max_peaks_per_region : if > 0, keep only the strongest N peaks per region.
    """
    cores: list[dict] = []
    for isl_start, isl_end in islands:
        sig = signal[isl_start:isl_end].astype(np.float64)
        if sig.size < 3 or sig.max() < min_summit_reads:
            continue
        _covered = sig[sig > 0]
        region_scale = (float(np.percentile(_covered, region_scale_pct))
                        if _covered.size else 0.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            pk, props = find_peaks(
                sig, prominence=min_prominence, height=min_summit_reads,
                distance=max(1, min_distance), width=max(1, min_peak_width),
                rel_height=rel_height,
            )
        if any("width" in str(w.message) for w in caught):
            log.warning(
                "island %d-%d: candidate peak(s) had near-zero measured width at "
                "--rel-height=%.3g and may have been dropped by --min-peak-width=%d "
                "(a low --rel-height measures width right at the summit tip, which "
                "conflicts with a positive --min-peak-width; raise --rel-height or "
                "lower --min-peak-width)",
                isl_start, isl_end, rel_height, min_peak_width)
        if len(pk) == 0:
            continue
        lefts = props["left_ips"]
        rights = props["right_ips"]
        proms = props["prominences"]
        island_cores = []
        for i, p in enumerate(pk):
            summit_val = float(sig[p])
            if float(proms[i]) < max(min_prominence, prominence_frac * region_scale):
                continue
            width = max(float(rights[i]) - float(lefts[i]), 1.0)
            if float(proms[i]) / width < min_steepness:
                continue
            start_i = isl_start + int(np.floor(lefts[i]))
            end_i = isl_start + int(np.ceil(rights[i]))
            summit = isl_start + int(p)
            base = max(summit_val - float(proms[i]), _EPS)
            island_cores.append(
                {
                    "start": start_i,
                    "end": max(end_i, start_i + 1),
                    "summit": summit,
                    "signal": summit_val,
                    "baseline": base,
                    "fold": summit_val / base,
                    "prominence": float(proms[i]),
                    "region_scale": region_scale,
                }
            )
        if max_peaks_per_region and len(island_cores) > max_peaks_per_region:
            island_cores.sort(key=lambda d: d["prominence"], reverse=True)
            island_cores = island_cores[:max_peaks_per_region]
        cores.extend(island_cores)
    return cores
