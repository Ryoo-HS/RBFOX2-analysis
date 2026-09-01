# RBP-Mediated Alternative Splicing Regulation in EMT

## Background

Epithelial-mesenchymal transition (EMT)의 스플라이싱 조절 네트워크를 RNA-binding protein(RBP) 관점에서 분석하는 프로젝트입니다.

- **RBFOX2**는 결합 위치에 따라 스플라이싱을 활성화하거나 억제하는 position-dependent splicing 조절자로 알려져 있습니다 (Yeo et al., 2009).
- **ESRP1**은 epithelial splicing regulator로서 RBFOX2와 길항적으로 작용하며 EMT 과정에서 상반된 스플라이싱 프로그램을 형성한다는 프레임워크가 제안되어 왔습니다 (Shapiro et al., 2011; Yang et al., 2016; Fici et al., 2017).
- 이 외에도 **QKI**, (예정: **PUM1**) 등 EMT 관련 RBP들의 결합 패턴과 상호작용을 CLIP-seq 기반으로 정리합니다.

CLIP-seq, RNA-seq 그리고 다중 RBP 간 co-binding 분석을 정리했습니다.

## Goals

1. hESC-CLIP-seq data기반 peak calling 및 분석 파이프라인 구축
2. RNA-seq data를 이용한 splicing pattern, expression change 분석
3. RBFOX2–ESRP1–QKI(–PUM1) 간 co-binding 패턴 분석 및 EMT 연관 타겟 유전자 도출

## Repository Structure

```
RBFOX2/
├── RNA-seq/       RBFOX2 KD vs control RNA-seq 검증 (sashimi plot)
└── CLIP/          RBFOX2 CLIP-seq 피크콜링
ESRP1/
└── CLIP/          ESRP1 CLIP-seq 피크콜링
QKI/
└── CLIP/          QKI CLIP-seq 피크콜링
PUM1/              (예정)
cobinding_analysis/  RBP 간 co-binding 통계 분석 (regioneR permutation 등)
docs/              전체 워크플로우 다이어그램
environment.yml    재현을 위한 conda 환경 정의
```

각 분석 폴더에는 별도 README가 있으며, 배경/방법/핵심 결과/재현 방법을 설명합니다.

## Methods Overview

- **Alignment / Noise filtering**: bowtie1 기반 6단계 순차 필터링(hairpin → miscRNA → rRNA → snRNA → snoRNA → tRNA) 후 hg19 genome mapping
- **Peak calling**: Coverage-based, MACS, CLIPper 세 가지 방법 비교 및 strand imputation
- **Annotation**: ChIPseeker, GENCODE v19 GTF 기반
- **Co-binding statistics**: regioneR permutation test
- **Splicing / expression validation**: rMATS, DESeq2, IGV sashimi plot

## Reproducibility

`environment.yml`에 명시된 conda 환경을 사용합니다. 원본 탐색 과정의 코드는 포함하지 않으며, 검증이 끝난 최종 스크립트만 정리해 업로드합니다.

## Status

진행 중 — 분석이 진행될 때마다 각 폴더와 이 README를 업데이트합니다.
