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

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def summit_offset(self) -> int:
        return self.summit - self.start
