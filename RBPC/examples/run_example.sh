#!/usr/bin/env bash
# Build the demo BAM (two adjacent mountains + a weak bump + a spliced read)
# and call peaks on it.
set -euo pipefail
cd "$(dirname "$0")"

python make_example_bam.py -o example.bam

rbpc \
    --bam example.bam \
    --normalize-method none \
    --min-prominence 20 \
    --min-summit-reads 30 \
    --min-distance 15 \
    --region-min-depth 3 \
    --rel-height 0.3 \
    --library-type forward \
    --bdg \
    -v \
    -o example.narrowPeak

echo "---- called peaks (expect the two mountains ~5000 and ~5300, split) ----"
cat example.narrowPeak
