# prominence-peakcaller

A rule-based peak caller for CLIP-seq / RBP-binding data. It follows the way a
researcher reads a coverage track in IGV, as a **two-stage** decision:

1. **Region (stage 1)** — first find the stretches that are genuinely covered
   ("this region has signal"), the way existing callers delimit enriched
   regions.
2. **Peaks (stage 2)** — then, *inside* each region, pick the spots that rise
   sharply above their **local surroundings** and come back down — the real
   mountains. "Prominence relative to the surroundings" is the core idea.

No p-value / q-value, no training, no annotation. Detection is built on
`scipy.signal.find_peaks` (topographic prominence), which naturally splits a
broad multi-hump region into one peak per hump.

## Install

```bash
pip install .
# or for development:
pip install -e ".[test]"
```

Requires Python ≥ 3.9, `pysam`, `numpy`, `scipy`. Linux. The BAM must be
coordinate-sorted and indexed (`samtools sort`, `samtools index`).

## Usage

```bash
peakcaller \
    --bam input.bam \
    --min-prominence 20 \
    --min-summit-reads 30 \
    --prominence-frac 0.2 \
    --min-steepness 0.5 \
    --rel-height 0.3 \
    --library-type forward \
    -o output.narrowPeak
```

`peakcaller -h` prints the full help.

## How it works

```
BAM ─▶ coverage (per-strand, splice-aware, pysam)
    ─▶ [normalize: none = raw counts (default) | rpm = per-million]
    ─▶ STAGE 1  regions = covered stretches with depth >= region-min-depth,
                          gaps <= region-gap merged, then regions narrower than
                          min-region-width (one read's width) dropped as
                          duplicate read stacks
    ─▶ STAGE 2  inside each region, scipy.find_peaks by topographic prominence:
                 a summit is a peak if it rises above its flanking valleys by
                   >= max(min-prominence, prominence-frac x region-scale)   AND
                 is steep enough (>= min-steepness) AND tall enough (min-summit-reads).
                 boundary width is set by rel-height; keep top max-peaks-per-region.
                 then drop peaks whose reads come from too few distinct start
                 positions (min-complexity) -- amplified, not independent, evidence.
    ─▶ SCORING  each surviving peak is also measured against a fixed +-context-window
                 local background with all called peaks masked out (depth_ratio);
                 this only ranks peaks, it never removes any.
    ─▶ narrowPeak / BED6  (+ optional bedGraph)
```

**Prominence** is how far a summit rises above the higher of the two valleys
that flank it — a *local* measure, so the same threshold behaves consistently
across regions of different depth. **`prominence-frac`** makes the bar
region-adaptive: on a tall gene a peak must rise a real fraction of the region's
scale (a high percentile of its covered depth), so small ripples and
shoulder-bumps on a dominant peak are dropped, while peaks in low-coverage
regions still pass via the absolute `min-prominence` floor.

## Options

| option | default | meaning |
|---|---|---|
| `--bam` | — | input BAM (sorted + indexed) |
| `-o, --output` | — | output peak file |
| `--format` | `narrowPeak` | `narrowPeak` (ENCODE BED6+4) or `bed` (BED6) |
| **peak detection** | | |
| `--min-prominence` | `5` | absolute floor: rise above flanking valleys (signal units); deliberately low — selectivity is meant to come from `--prominence-frac` and stage 1 |
| `--prominence-frac` | `0` | region-adaptive: also require prominence ≥ frac × region scale (e.g. 0.2) |
| `--region-scale-pct` | `90` | percentile of the region's covered depth used as its scale |
| `--min-steepness` | `0` | min prominence / width (rise per bp); drops gentle broad bumps (e.g. 0.5) |
| `--min-summit-reads` | `1` | absolute height floor at the summit |
| `--min-complexity` | *= `--min-steepness`* | min distinct read start positions (independent molecules) per bp of peak width; rejects tall humps built from a few heavily amplified molecules. `0` = off, see below |
| `--context-window` | `1000` | bp each side used as the local background for the `depth_ratio` **scoring** component (peak depth ÷ background depth, all called peaks masked out). Never filters a peak. `0` = skip, see Output |
| `--min-distance` | `15` | minimum bp between two summits |
| `--min-peak-width` | `5` | minimum peak width (bp) |
| `--flank` | `10` | widen each called peak by N bp on both sides in the output (e.g. for motif search near, not just at, the summit); does not affect detection |
| `--rel-height` | `0.3` | boundary width point (0 = tip, 1 = base); lower = tighter peak |
| `--max-peaks-per-region` | `0` | 0 = unlimited; N = keep only the strongest N per region/gene |
| **region (stage 1)** | | |
| `--region-min-depth` | *auto* | a region needs depth ≥ this (enrichment gate; prunes stage-2 work). Auto = a percentile of this BAM's own per-island max depth |
| `--region-min-depth-pct` | `75` | percentile used for the auto `--region-min-depth`; raise it for fewer, stronger regions |
| `--region-gap` | `200` | max gap (bp, depth < region-min-depth) merged into one region |
| `--min-region-reads` | `0` | drop regions with fewer than ~N reads (0 = off) |
| `--min-region-width` | *auto* | drop regions narrower than N bp (0 = off). Auto = `factor` × this BAM's read length — the main artifact filter, see below |
| `--min-region-width-factor` | `2` | auto `--min-region-width` = this many read lengths (90th pct of aligned read length) |
| **coverage building** (memory chunking, not biological) | | |
| `--coverage-gap` | `200` | merge raw read blocks within this many bp into one coverage chunk; keep ≥ `--region-gap` |
| **normalization** | | |
| `--normalize-method` | `none` | `none` = raw counts (thresholds mean "reads"); `rpm` = per-million |
| `--scale-to` | `1e6` | reads to scale to under rpm |
| **BAM handling** | | |
| `--library-type` | `forward` | `forward`, `reverse` (typical single-end eCLIP), `unstranded` |
| `--min-mapq` | `0` | minimum mapping quality |
| `--min-read-length` | `0` | drop reads shorter than this (bp) |
| `--keep-dup` | off | keep duplicate-flagged reads (default: drop) |
| `--bdg` | off | also write per-strand bedGraph tracks |
| `--chrom` | all | restrict to a chromosome (repeatable) |

> **Strand — verify this per library, do not assume.** Single-end eCLIP is
> *often* reverse-stranded, but preprocessing can already have flipped the
> reads. This only relabels column 6 (peak coordinates, prominence and score
> are byte-identical either way), but getting it wrong inverts every downstream
> motif or gene assignment.
>
> Two cheap checks that don't need IGV, both of which decided it for the three
> BAMs in `~/test` (all three turned out to be **`forward`**, because the
> upstream pipeline had already reverse-complemented the reads — see
> `~/test/README.md`):
> 1. **Gene-strand agreement** — in loci with no antisense gene, count reads
>    aligning on the annotated gene's strand. 86–88% agreement means forward;
>    a similar majority the other way means reverse.
> 2. **Motif orientation** — scan enriched k-mers near summits on both
>    orientations. The one that recovers the protein's canonical motif
>    (e.g. RBFOX2 `TGCATG`, QKI `ACTAAC`) is the correct strand.

> **`--coverage-gap` vs `--region-gap`.** These look similar but serve different
> stages: `--coverage-gap` merges raw read blocks into per-chromosome array chunks
> (stage 0, purely a memory/chunking concern — never a whole-chromosome array);
> `--region-gap` merges depth-below-threshold gaps into one *region* (stage 1, a
> biologically meaningful choice). Because stage 0 runs first, a `--coverage-gap`
> smaller than `--region-gap` silently caps how far apart two regions can ever be
> merged. Keep `--coverage-gap >= --region-gap` (the default keeps both at 200).

> **Stage-1 auto-calibration.** Both stage-1 gates default to a value measured
> from the BAM itself, so a new library needs no hand-tuning:
> `--region-min-depth` is a percentile (`--region-min-depth-pct`, default 75) of
> the per-coverage-island **max depth** population — on real data ~66–75% of
> coverage islands are a single overlapping read, and this discards them;
> `--min-region-width` is `--min-region-width-factor` × the 90th percentile of
> this BAM's **aligned read length**. Neither is a statistical test — they are
> descriptions of what the library contains. Pass an explicit number to either
> option to override (`--min-region-width 0` turns the width filter off).
>
> **Why width is the primary artifact filter.** A pile of reads that all start
> at the *same* position — PCR/optical near-duplicates, the dominant CLIP
> artifact — makes a hump purely from read-edge geometry, and that hump sits on
> a zero baseline, so its prominence equals its full height and it passes any
> *relative* test trivially (`region_scale` is its own depth). Width is what
> gives it away: a region covered only by reads sharing one start position can
> never be wider than the longest read in the library. Measured on real data
> (ESRP1, chr1): 68% of sub-30bp regions are a single distinct read start, and
> that fraction hits exactly 0% at 75bp — the library's longest read is 70bp.
> The default (2 × p90 read length ≈ 88bp there) sits just past that point.
> In practice this drops ~78% of stage-1 regions but only ~6% of called peaks:
> it is an artifact/compute filter, **not** a peak-count dial. Use
> `--region-min-depth-pct` for overall stringency.
>
> **`--min-complexity` — how many *molecules*, not how many reads.** The width
> filter above only catches stacks narrower than a read. A handful of amplified
> stacks a few hundred bp apart merge into a wide region, pass stage 1, and then
> pass stage 2 for the same reason: each hump sits on a zero baseline, so its
> prominence *is* its height and every relative test is satisfied trivially.
> Making the prominence floor taller does not fix this — depth counts PCR copies,
> so a magnitude floor does not transfer between libraries of different depth.
>
> So instead of asking how tall a peak is, this filter asks **how many
> independent molecules built it**. For single-end reads the alignment start
> position is the standard duplicate signature (what `samtools markdup`
> collapses on), so
>
> ```
> width_complexity = (distinct read start positions overlapping the peak) / peak width
> ```
>
> is "independent molecules per bp". The threshold is **not** a new tuned
> number: `width_complexity` is an *upper bound on the peak's own steepness
> measured on deduplicated coverage*. Keep one read per distinct start and every
> position inside the peak is covered at most `n_distinct` times, so the
> deduplicated prominence — and hence `prominence / width` — cannot exceed
> `n_distinct / width`. A peak below `--min-steepness` on this measure therefore
> **provably** fails the steepness bar you already set, once duplicates are
> collapsed; that is why the default is `--min-steepness` itself. (Measured on
> ESRP1 chr18 the bound holds for 100% of peaks and is tight: the median ratio of
> actual deduplicated steepness to `width_complexity` is 1.00.) With the default
> `--min-steepness 0` the filter is off, so nothing changes unless you ask for
> steepness at all. Under `--normalize-method rpm` the floor is converted to
> molecule units automatically.
>
> It reverses a magnitude floor's decisions in both directions, which is the
> point: on ESRP1 it **drops** a prominence-**1631** peak built from **7**
> distinct molecules (233 PCR copies each) and **keeps** prominence-5 peaks built
> from 5–11 independent ones. All 8 peaks of the known CDH2 cluster
> (chr18:25.53Mb) survive, at `width_complexity` 0.89–3.22.
>
> **It is a duplication filter, so its effect depends on the library, not on a
> dial setting** (whole genome, `--min-steepness 0.5`): ESRP1 73,678 → 15,042
> (−80%), RBFOX2 415,711 → 48,187 (−88%), QKI 78,523 → 6,889 (−91%) — but YBX1
> only 383,478 → 288,951 (−25%), because that library genuinely is not
> duplicated (median 2.4 reads per molecule, versus 6.5–11 for the other
> three). That is the metric reporting a
> property of the data rather than imposing a target count.
>
> **`--min-complexity`: count molecules, not reads.** Width (above) catches a
> stack narrower than one read, but not a *wide* region built from a handful of
> heavily amplified molecules — and such a hump sits on a zero baseline too, so
> its prominence equals its full height and it passes every *relative* test
> trivially. The fix is not a taller prominence floor (a magnitude floor does
> not transfer across libraries, because depth counts PCR copies, not
> molecules). It is to ask **how many independent molecules made this peak**:
> for single-end data the alignment start position is the standard duplicate
> signature, so
>
> ```
> width_complexity = distinct read start positions / peak width   (molecules per bp)
> ```
>
> **The threshold is derived, not tuned.** `width_complexity` is an upper bound
> on the peak's own steepness measured on *deduplicated* coverage: keep one read
> per distinct start and every position in the peak is covered at most
> `n_distinct` times, so `steepness_dedup <= n_distinct / width`. A peak below
> `--min-steepness` on this scale therefore **provably cannot** clear the
> steepness bar you already set once PCR duplicates are collapsed — which is why
> the default is `--min-steepness` itself rather than a new constant. Measured on
> ESRP1 chr18 the bound holds for 100% of peaks and is *tight* (median ratio of
> real deduplicated steepness to `width_complexity` = 1.00). With the default
> `--min-steepness 0` the filter is off; pass a number to set it independently,
> or `0` to disable.
>
> It filters by *kind* of evidence, not amount — on ESRP1 it drops a
> prominence-**1631** hump made of **7** molecules (233 PCR copies each) while
> keeping prominence-**5** peaks made of 5–11 independent molecules; a magnitude
> floor does exactly the reverse. Being a duplication filter, how much it removes
> depends on how duplicated the library is (real data, `--min-steepness 0.5`):
> ESRP1 73,678 → 15,042, RBFOX2 415,711 → 48,187, QKI 78,523 → 6,889, but YBX1
> 383,478 → 288,951 because that library genuinely is not duplicated (median 2.4
> copies per molecule vs 6.5–11 for the others). Like `--min-region-width`, it is
> an artifact filter, not a peak-count dial.

> **No significance testing — including in the calibration.** There is an
> opt-in `--min-prominence-auto-fdr` that picks a prominence floor by empirical
> FDR against a read-block-shuffled null. It works, but it is deliberately not
> the default: "accept peaks with FDR ≤ 5%" is a significance test, and this
> caller is designed without one (columns 8/9 stay `-1`).

> **`--flank`.** Widens each already-called peak symmetrically after detection
> (e.g. `--flank 20` turns a 30 bp peak into a 70 bp one, centered the same) —
> useful for motif search, since the motif may sit near the summit rather than
> exactly inside the tight `--rel-height`-bounded core. Purely a post-hoc output
> transform: it does not change which peaks get called, only the reported
> `start`/`end` (clamped to `[0, chrom length)`). Neighboring flanked peaks are
> not merged even if they end up overlapping.

## Output — narrowPeak columns

Standard ENCODE narrowPeak (BED6+4):

| col | field | value |
|---|---|---|
| 5 | score | composite, ~0–1000 (no hard ceiling, see below) |
| 7 | signalValue | the peak's **prominence** |
| 8, 9 | pValue, qValue | `-1` (no statistical test, by design) |
| 10 | peak | summit offset from the peak start |

The **score** combines four normalized factors of peak quality, since no
p/q-value is computed: prominence (log-scaled), region significance
(prominence ÷ region scale), steepness (prominence ÷ width), and local
background contrast (`depth_ratio`, see below), weighted
0.40 / 0.15 / 0.25 / 0.20. Each factor is scaled against the 95th percentile of that
same factor across **all peaks called in this run** (not a fixed constant) —
a peak sitting exactly at that reference on all four axes scores ~1000; a
peak well beyond it (the strongest of the strong) scores higher still, rather
than every strong peak collapsing onto a flat ceiling the way fixed constants
did on real (non-subset) data, where nearly every called peak used to hit
1000. There is deliberately **no hard cap at 1000** for this reason.
**Score is therefore only comparable *within* one run/output file** — the
same peak can get a different score in a different run (e.g. a different
`--chrom` selection or threshold set changes the peak population it's scaled
against).

> **`depth_ratio` — contrast against a *fixed* window, not a self-drawn one.**
> The score's fourth factor is the peak's mean depth divided by the mean depth
> of the same-strand sequence within `--context-window` bp (default 1000) on
> either side, with **every called peak masked out of that background** (so a
> cluster of genuine neighbouring sites does not suppress its own members).
>
> It exists because the region-relative factor above is *self-referential* for
> an isolated hump: `region_scale` is measured inside a region whose boundary
> the hump itself drew, so a hump standing alone becomes its own baseline and
> the test `prominence >= prominence-frac x region-scale` is satisfied no
> matter what. A fixed-width window cannot be gamed that way. Being a ratio of
> depths, it is also invariant to sequencing depth and to `--normalize-method`
> — the property an absolute magnitude floor lacks.
>
> Measured against an independent label (presence of each protein's canonical
> motif near the summit, corrected for local sequence composition) on three
> eCLIP libraries, it out-ranked every other component including prominence,
> and it was the **only** one that still discriminated among *weak* peaks,
> where prominence carries essentially no information. Adding it improved the
> score's own ranking on all three libraries. Full numbers and the rejected
> alternatives are in `peakcaller/context.py`.
>
> It is **only a scoring component — it never filters a peak.** Its
> distribution has no elbow, and the project's one externally confirmed
> anchor (the CDH2 cluster) sits mid-distribution on it, so a hard cutoff
> would delete peaks known to be real. `--context-window 0` skips the
> measurement (costs ~7% runtime).

## Tuning — important

The absolute thresholds (`--min-prominence`, `--min-summit-reads`, `--region-min-depth`)
are in **signal units**, so their meaning depends on how the run is normalized:

- `--normalize-method none` (default): units are **raw reads** — `--min-prominence 20`
  means "rises 20 reads above its surroundings". Interpretable per sample.
- `--normalize-method rpm`: units are **RPM**, which depend on library size.

Because coverage depth scales with sequencing depth, **these absolute thresholds
must be tuned per dataset**. A worked example on a ~50k-read multi-gene subset
that matched a researcher's IGV calls:

```
--min-prominence 20  --min-summit-reads 30  --prominence-frac 0.2
--min-steepness 0.5  --rel-height 0.3  --min-distance 20  --region-min-depth 3
```

Moving to a full BAM (millions of reads), raise the *absolute* knobs
(`--min-prominence`, `--min-summit-reads`) proportionally to the larger depths; the
*relative* knobs (`--prominence-frac`, `--min-steepness`, `--rel-height`,
`--region-scale-pct`) are ratios and carry over unchanged.

### Worked example: scaling by library size

The easiest recipe is to scale the two absolute knobs by the **depth ratio**
between libraries. Pick one clear peak you trust and read its coverage depth in
IGV for each library:

| | reads in BAM | depth at a typical peak | `--min-prominence` | `--min-summit-reads` |
|---|---|---|---|---|
| tuned subset | ~50k | ~40 | 20 | 30 |
| full library | ~20M | ~1600 (≈ 40×) | **20 × 40 = 800** | **30 × 40 = 1200** |

So for the ~20M-read full BAM you would use roughly:

```bash
peakcaller --bam full.bam \
    --min-prominence 800 --min-summit-reads 1200 \
    --prominence-frac 0.2 --min-steepness 0.5 --rel-height 0.3 \
    --region-min-depth 120 --min-distance 20 \
    --library-type forward -o full.narrowPeak
```

(`--region-min-depth` scales the same way, 3 × 40 ≈ 120.) The exact multiplier
is the ratio of typical peak depths, not necessarily the read-count ratio, so
eyeball one or two peaks in IGV and adjust. Alternatively, run with
`--normalize-method rpm` and tune the thresholds once in RPM units — then the
same numbers transfer across libraries of different depth.

Quick guide: fewer/cleaner peaks → raise `--min-prominence`, `--min-summit-reads`,
`--min-steepness`, or `--prominence-frac`. Fewer *duplication-driven* peaks
specifically (independent of how tall they are) → raise `--min-complexity`. Tighter boundaries → lower
`--rel-height`. One peak per gene → `--max-peaks-per-region 1`.

## Example

```bash
cd examples && bash run_example.sh
```

Builds a small synthetic BAM (two adjacent mountains that must be split, a weak
bump that must be ignored, a spliced read) and calls peaks on it.

## Tests

```bash
pytest
```

## Memory / scaling

Coverage is built **per chromosome as small per-island arrays** (only the
covered stretches are allocated, never a whole-chromosome array), so peak memory
is bounded by the largest covered island — typically a single gene (kilobases).
Large / whole-genome BAMs therefore run without large allocations; RAM does not
grow with BAM size. (`--chrom` restricts processing; under `rpm` the scaling
factor is then computed from the fetched chromosomes only.)

## License

MIT.
