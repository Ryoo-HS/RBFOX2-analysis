"""Core data structures."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Peak:
    chrom: str
    start: int
    end: int
    strand: str
    summit: int
    signal: float
    baseline: float
    fold: float
    region_scale: float = 0.0
    # Local background contrast (see context.py). inf = no usable background;
    # output._scores() replaces that with the run's maximum finite value.
    depth_ratio: float = 1.0
    peak_depth: float = 0.0
    bg_depth: float = 0.0
    name: str = "peak"
    # Highest --summit-margins tier (1 = T1, 2 = T2, ...) this peak's summit
    # height clears. None = --summit-margins not requested (unchanged default).
    tier: int | None = None

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def summit_offset(self) -> int:
        return self.summit - self.start
