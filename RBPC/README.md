# RBPC

CLIP-seq 데이터에서 RNA-binding protein의 결합 부위(peak)를 찾는 도구.

통계적 유의성 검정(p-value / q-value) 대신 **주변부 대비 신호가 얼마나 솟아 있는지(prominence)** 로 peak을 판정한다.
IGV에서 coverage track을 눈으로 보고 판단하는 방식을 코드로 옮긴 것에 가깝다.

> **개발 중단 (2026-09).** 학부 연구를 마무리하면서 개발을 멈췄다.
> 배포용이 아니라 작업 기록을 남기는 목적의 저장소다.
> 아래 "한계"에 적힌 검증들을 하지 않았으므로 결과를 그대로 연구에 쓰기는 어렵다.

---

## 동작 방식

두 단계로 나뉜다.

**Stage 1 — 후보 구간 찾기**
read가 일정 깊이 이상 깔린 구간(region)을 찾는다. 절대 depth 기준.

**Stage 2 — peak 판정**
각 구간 안에서 주변보다 뚜렷하게 솟은 지점을 peak으로 부른다. 모양(prominence) 기준.

p-value / q-value는 계산하지 않는다. narrowPeak 출력의 8, 9번 컬럼은 `-1`로 고정된다.

---

## 자동으로 정해지는 값

근거를 세울 수 있었던 것만 자동 계산한다.

| 값 | 계산 방식 | 근거 |
|--|--|--|
| `min-region-width` | `2 × p90(read length)` | 이 폭보다 넓으면 PCR duplication 스택으로 설명할 수 없음 |
| `min-complexity` | `min-steepness`와 동일 | dedup 시 steepness의 상한이라는 수학적 관계 |
| summit floor | `ceil(bg_trim25)` | BAM의 노이즈 바닥 (covered base depth에서 상위 25% 제외한 평균) |

나머지(`prominence-frac`, `min-steepness`, `rel-height`)는 사용자가 정하는 값이다.
자동 결정을 세 번 시도했지만 depth 분포에 뚜렷한 경계(elbow)가 없어 근거를 만들 수 없었다.

---

## `--summit-margins`: 결과를 3단계로 나누기

threshold 하나를 근거 있게 정하는 게 불가능해서, **여러 엄격도의 결과를 한 번에 내놓는** 방식을 택했다.

```bash
rbpc --bam sample.bam -o out.narrowPeak \
  --library-type forward --summit-margins 18,48,98 --split-tiers
```

1. BAM에서 노이즈 바닥(floor)을 자동으로 구한다
2. 실제 기준 = floor + margin (항상 더하므로 노이즈 바닥 아래로는 못 내려간다)
3. 각 peak에 통과한 최고 단계를 붙인다 → 이름이 `peak_123_T2` 형태가 된다
4. `--split-tiers`를 주면 단계별 파일도 따로 나온다

예시값 `18,48,98`은 floor가 2일 때 20 / 50 / 100 reads가 된다.
이 값은 관측된 summit 높이 분포를 보고 **세 구간이 고르게 나뉘도록** 잡은 것이지,
"20 reads부터 진짜 결합"이라는 근거가 있는 건 아니다. 엄격도 라벨로 보면 된다.

옵션을 지정하지 않으면 기존 동작(단일 기준) 그대로다.

### 다른 방식들을 기각한 이유

| 시도한 방식 | 결과 |
|--|--|
| percentile로 자동 결정 | 분포가 매끄러워 경계가 없음. 값의 근거를 못 세움 |
| 노이즈 바닥의 k배 | trim 비율을 바꾸면 라이브러리 간 순위가 뒤집힘 → 노이즈 값에 신호가 섞여 있었음 |
| permutation 기반 FDR | 유의성 검정이라 이 도구의 방향과 맞지 않음 (코드는 `--min-prominence-auto-fdr` 뒤에 남겨둠) |

진단해보니 4개 라이브러리의 노이즈 바닥이 covered base 당 1.4~1.8 reads로 거의 같았다.
라이브러리별 보정이 필요한 문제가 아니라, 노이즈 바닥 위에 어디에 선을 그을지가 문제였다.

---

## 옵션

### 주요 옵션

아래 실행 결과에서 실제로 지정한 것들이다.

| 옵션 | 설명 | 기본값 |
|--|--|--|
| `--bam` | 입력 BAM (정렬 + 인덱스 필요) | 필수 |
| `-o` | 출력 파일 | 필수 |
| `--library-type` | forward / reverse / unstranded. **가정하지 말고 확인할 것** | forward |
| `--normalize-method` | none / rpm | none |
| `--prominence-frac` | region 자체 규모 대비 prominence 비율 | 0.0 (off) |
| `--min-steepness` | 기울기 하한 (prominence / width) | 0.0 (off) |
| `--rel-height` | peak 경계를 잡는 높이 지점 (0=꼭대기, 1=바닥) | 0.3 |
| `--min-distance` | summit 간 최소 거리 (bp) | 15 |
| `--summit-margins` | 3단계 summit 기준. 예: `18,48,98` | 없음(off) |
| `--split-tiers` | 단계별 파일 추가 출력 | off |

`--prominence-frac`과 `--min-steepness`는 기본값이 off다. 아래 실행 결과처럼 쓰려면 직접 지정해야 한다.

<details>
<summary><b>전체 옵션</b> (대부분 기본값으로 동작)</summary>

**Stage 1 — 구간 탐색**

| 옵션 | 설명 | 기본값 |
|--|--|--|
| `--region-min-depth` | 구간으로 인정할 깊이 기준 | auto |
| `--region-min-depth-pct` | auto 계산에 쓰는 percentile | 75.0 |
| `--min-region-width` | 구간 최소 길이 (bp) | auto |
| `--min-region-width-factor` | auto 계산 시 read length의 몇 배로 할지 | 2.0 |
| `--region-gap` | 이 간격 이내는 한 구간으로 병합 (bp) | 200 |
| `--min-region-reads` | 구간 최소 read 수 | 없음 |
| `--coverage-gap` | coverage 청크 병합 간격 (bp, 메모리용) | 200 |

**Stage 2 — peak 판정**

| 옵션 | 설명 | 기본값 |
|--|--|--|
| `--min-prominence` | prominence 절대 하한 | 5.0 |
| `--region-scale-pct` | region 규모를 정하는 percentile | 90.0 |
| `--min-summit-reads` | summit 높이 하한 | 1.0 |
| `--min-complexity` | read 시작 위치 다양성 하한 (PCR 증폭 덩어리 제거) | `--min-steepness` |
| `--min-corrected-steepness` | background를 뺀 기울기 하한 | 0.0 (off) |
| `--min-peak-width` | peak 최소 너비 (bp) | 5 |
| `--max-peaks-per-region` | region 당 최대 peak 수 (0=무제한) | 0 |
| `--context-window` | depth_ratio 계산용 local background 범위 (bp) | 1000 |

**입출력**

| 옵션 | 설명 | 기본값 |
|--|--|--|
| `--format` | narrowPeak / bed | narrowPeak |
| `--flank` | 출력 시 peak을 양쪽으로 넓힘 (bp, 검출에는 영향 없음) | 10 |
| `--bdg` | bedGraph 추가 출력 | off |
| `--chrom` | 특정 염색체만 처리 | 전체 |
| `--min-mapq` | MAPQ 하한 | 0 |
| `--min-read-length` | read 길이 하한 | 0 |
| `--keep-dup` | duplicate read 유지 | off |
| `--scale-to` | rpm 환산 기준값 | 1000000 |

**opt-in (기본 비활성)**

`--min-prominence-auto-fdr`, `--target-fdr`, `--fdr-n-shuffles` — permutation 기반 empirical FDR로 `--min-prominence`를 자동 결정한다. 유의성 검정이라 이 도구의 방향과 맞지 않아 기본값에서 제외했고, 코드만 남겨뒀다.

</details>

`--summit-margins`와 `--min-summit-reads`는 AND 조건이다. 실제 T1 기준 = 둘 중 큰 값.

배경 depth를 미리 확인하려면 `scripts/diag_bg.py`를 쓰면 된다.

strand 검증 방법, coverage-gap과 region-gap의 차이, min-complexity 기본값의 유도 근거 등
상세한 설계 메모는 [docs/design-notes.md](docs/design-notes.md)에 있다.

### Score

```
score = 0.40 × 강도 + 0.15 × region 상대값 + 0.25 × 가파름 + 0.20 × depth_ratio
```

`depth_ratio`(±1kb 주변 대비)는 필터가 아니라 점수 성분으로만 쓴다.

---

## 실행 결과 (hg19, whole genome)

```
--library-type forward --normalize-method none --prominence-frac 0.2
--min-steepness 0.5 --rel-height 0.3 --min-distance 20 --summit-margins 18,48,98
```
(`region-min-depth`, `min-region-width`는 지정하지 않아 auto)

| RBP | region-min-depth | min-region-width | floor | 전체 peak | T1 | T2 | T3 | 시간 |
|--|--|--|--|--|--|--|--|--|
| ESRP1 | 2 | 88 bp | 2 | 12,763 | 5,189 | 3,873 | 3,701 | 1m51s |
| RBFOX2 | 2 | 84 bp | 2 | 43,612 | 14,318 | 13,749 | 15,545 | 6m24s |
| YBX1 | 3 | 100 bp | 2 | 233,631 | 115,433 | 65,115 | 53,083 | 7m07s |
| QKI | 2 | 94 bp | 2 | 6,136 | 2,038 | 1,821 | 2,277 | 2m09s |

### Summit 높이 분포 (raw reads)

`samtools depth`로 summit 좌표를 재조회한 근사값이라 내부 계산과 미세한 차이가 있을 수 있다.

| RBP | p10 | p25 | p50 | p75 | p90 |
|--|--|--|--|--|--|
| ESRP1 | 25 | 35 | 60 | 110 | 201 |
| RBFOX2 | 28 | 42 | 72 | 132 | 256 |
| YBX1 | 24 | 32 | 50 | 92 | 180 |
| QKI | 29 | 42 | 73 | 143 | 309 |

### 알려진 타깃 확인

EMT 관련 유전자 GJA1, CDH1에서 peak이 남아 있는지만 확인했다. 놓치지 않았다는 확인이고, 단계 구분이 맞는지를 검증한 건 아니다.

| RBP | GJA1 | CDH1 |
|--|--|--|
| ESRP1 | 12 (전부 T3) | 13 (T1×1, T2×4, T3×8) |
| RBFOX2 | 12 (전부 T3) | 8 (T1×1, T2×1, T3×6) |
| YBX1 | 31 (전부 T3) | 52 (T1×1, T3×51) |
| QKI | 3 (T1×1, T2×2) | 1 (T2×1) |

QKI만 두 유전자 모두 T3가 없다. 결합 강도 차이일 수 있으나 확인하지 않았다.

---

## 한계 / 하지 못한 것

1. **IGV 대조를 하지 않았다.** 단계 경계 부근의 peak이 실제로 결합 부위처럼 보이는지 눈으로 확인하지 않았다. 남은 작업 중 가장 중요한 것.
2. **생물학적 검증이 부족하다.** motif 농축은 일부만 확인했고, conservation 분석은 하지 않았다. replicate가 없어 재현성 확인도 불가능했다.
3. **기존 도구와 비교하지 않았다.** MACS3, CLIPper, PureCLIP과의 벤치마크 미수행.
4. **RBP 간 peak 수를 직접 비교할 수 없다.** 절대 depth 기준이라 깊게 시퀀싱된 라이브러리가 더 많이 통과한다. 비교하려면 downsampling이 필요하다. 위 표의 YBX1(233,631개)이 다른 RBP보다 크게 많은데, 실제 결합 특성인지 아티팩트인지 판단하지 않았다.
5. **사용자가 정해야 하는 값이 남아 있다.** `prominence-frac`, `min-steepness`, `rel-height`는 자동화에 실패해 사용자 입력으로 남겼다.
6. **depth_ratio 기반 기준을 못 만들었다.** background를 뺀 뒤 기울기를 재적용하는 방안을 검토했지만 ESRP1에서 신호와 역상관이 나와 채택하지 않았다. ESRP1은 모티프(GU-rich)가 약해 기준으로 삼기 어렵다.

---

## 설치 / 사용

```bash
pip install -e .

rbpc --bam sample.sorted.bam -o out.narrowPeak --library-type forward
```

테스트: `pytest`

---

## 참고

- `--summit-margins`의 floor는 가장 큰 상염색체 4개를 샘플링해 계산한다. 특정 염색체만 담긴 subset BAM에서는 샘플이 비어 fallback 값이 쓰인다. whole-genome BAM에서는 문제없다.
- `--normalize-method rpm`을 쓰면 raw 단위로 계산된 floor가 신호 단위로 자동 환산된다.
