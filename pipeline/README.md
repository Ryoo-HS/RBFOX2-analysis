# CLIP-seq Peak Calling & Annotation Pipeline

RBFOX2, ESRP1, QKI (그리고 추후 PUM1) CLIP-seq 데이터에 공통으로 적용하는 범용 파이프라인입니다. RBP별 실행 파라미터와 결과 해석은 각 RBP 폴더(`../RBFOX2/`, `../ESRP1/`, `../QKI/`)의 README를 참고하세요.

## Pipeline Steps

### 1. Peak calling & noise filtering (`01_peak_calling.R`)
BAM → strand별 coverage → peak 콜링 후, 아래 순서로 노이즈 제거:
`width/max artifact` → `5S/7S/tRNA/rRNA repeat` → `lncRNA` (선택) → `snoRNA/miRNA` → `NUMT` (선택)

- **width/max artifact**: 폭이 좁은데 signal이 비정상적으로 높은 peak을 제거합니다. mapping artifact 또는 PCR 과증폭(PCR duplication/amplification bias)으로 특정 위치에 read가 비정상적으로 쌓인 경우로 판정합니다.
- **`repeat_hg19`/`sno_mi`는 항상 적용되는 기본 필터, `lnc_gr`/`numt_gr`(NUMT: 핵 내 미토콘드리아 유래 서열)는 RBP·분석 목적에 따라 선택적으로 켜고 끄는 필터**입니다. 해당 인자를 넘기지 않으면 그 단계는 자동으로 스킵됩니다.
- 이번 분석에서는 `lnc_gr`, `numt_gr` 모두 미투입 상태로 실행함 (필요 시 추후 투입 가능)

### 2. Signal threshold filtering (`02_threshold_filter.R`)
RBP별 max signal threshold로 최종 peak set을 확정합니다.

| RBP | Threshold | 최종 peak 수 |
|---|---|---|
| RBFOX2 | max > 65 | 23,652 |
| ESRP1 | max > 30 | 20,514 |
| QKI | max > 35 | 22,447 |

**방법론 한계**: 이 threshold는 통계적 유의성(p-value/FDR) 기준이 아니라, RBP 간 raw peak 개수 편차(QKI 116k / RBFOX2 352k / ESRP1 54.5k)를 downstream 비교(co-binding 등)에서 공정하게 다루기 위해 각 RBP의 최종 peak 수를 비슷한 규모(~2만 개)로 맞춘 **경험적(empirical) 컷오프**입니다. Saturation curve, IDR 등 원칙적 기준은 적용하지 않았습니다.

### 3. Annotation & GO/KEGG (`03_annotation_go_kegg.R`)
ChIPseeker annotation, region(5'UTR/exon/intron/3'UTR)별 GO(BP/CC/MF)·KEGG enrichment, ORFik metagene plot을 생성합니다.

## Requirements

```
GenomicAlignments, GenomicRanges, chipseq,
TxDb.Hsapiens.UCSC.hg19.knownGene, GenomicFeatures,
ChIPseeker, ORFik, data.table, clusterProfiler, org.Hs.eg.db
```

## Usage

```r
source("scripts/01_peak_calling.R")
source("scripts/02_threshold_filter.R")
source("scripts/03_annotation_go_kegg.R")
```
