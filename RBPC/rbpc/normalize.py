"""Depth normalization (MACS-style SPMR / RPM)."""
from __future__ import annotations

import numpy as np


def scaling_factor(total_reads, scale_to=1e6, method="rpm"):
    if method == "none":
        return 1.0
    if method != "rpm":
        raise ValueError(f"unknown normalize method: {method}")
    if total_reads <= 0:
        return 0.0
    return scale_to / float(total_reads)


def to_rpm(depth, factor):
    return depth.astype(np.float32) * np.float32(factor)
