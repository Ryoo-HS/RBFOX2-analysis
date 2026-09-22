# Peak Caller 스펙 검토 및 보충

> **최종 설계 (실데이터 반복 튜닝 후 확정)** — 아래 원본 검토 이후, 형의 실데이터(ESRP1/CDH2 등)와 IGV 라벨로 여러 번 튜닝하면서 핵심 알고리즘이 **topographic prominence 기반**으로 확정되었다. 최종 구조는 README 참조. 요약:
>
> 1. **Stage 1 — region**: depth ≥ `region-min-depth`인 구간(gap ≤ `region-gap` 병합)만 골라 Stage 2로 넘김 (enriched 구역 + 계산 효율).
> 2. **Stage 2 — 봉우리**: `scipy.signal.find_peaks`로 각 국소 최대점의 **prominence(양옆 계곡 대비 솟은 높이)**를 재서, `max(min-prominence, prominence-frac × region-스케일)`을 넘고 충분히 가파른(`min-steepness`)·높은(`min-summit-reads`) 것만 채택. 경계 폭은 `rel-height`, gene당 개수는 `max-peaks-per-region`.
> 3. **출력**: narrowPeak/BED, 종합 score(강도+region상대+가파름 가중합).
>
> "prominence = 주변부 대비"가 핵심이고, `prominence-frac`이 region 스케일 대비 상대 기준을 담당한다(큰 유전자의 잔봉우리·어깨봉우리 제거). fold/plateau/baseline-window 방식은 이 과정에서 prominence 방식으로 대체되었다.
>
> ---
>
> **아래는 초기 프로토타입 검토 기록** (역사적 참고용 — 일부는 위 최종 설계로 대체됨).

원본 `peak_caller_spec.md`를 검토하고, 프로토타입 구현 과정에서 **불명확하거나 빠져 있던 지점**을 구체화한 기록. 원본의 설계 철학(규칙 기반 local prominence + 절대 floor, island→core 2단계, RPM 정규화, 통계모델 배제)은 그대로 유지했고, MACS3는 개념(SPMR/RPM, local background window)만 참고했다.

## 1. 원본이 잘 잡아둔 부분 (그대로 유지)

- **두 조건 AND 설계**(3.1): fold_threshold(상대) + min_reads(절대)를 AND로 묶어, "baseline≈0에서 1→3 read가 3배로 잡히는 함정"과 "depth 높은 유전자의 완만한 굴곡이 다 잡히는 문제"를 동시에 방지 — 이 핵심 논리는 구현의 중심축.
- **island→core 2단계**(4): 넓은 island 전체가 아니라 내부의 강한 core만 출력.
- **RPM 정규화**(3.2), **fold 만족 연속구간 = peak 경계**(3.3), **q-value 배제**(3.4).
- **narrowPeak/BED6 출력, pysam+numpy 최소 의존성, pip 배포**(8).

## 2. 불명확해서 구체화한 지점

| 원본에서 모호했던 점 | 구현에서 확정한 내용 |
|---|---|
| "local baseline = baseline-window 내 최소값 또는 양끝 계곡"(5-4) — min? 계곡? 어느 범위? | island의 **양끝 flank를 포함한 국소 window**(스펙 "양끝 계곡")에서 **저백분위(기본 10th, `--baseline-method percentile`)** 또는 **min**으로 계산. 실데이터 검증에서 CLIP peak가 날카로운 단일 stack이라 "island 내부"만 보면 계곡이 없어 fold≈1로 놓치는 문제가 드러나 flank 포함 방식으로 수정. 속도를 위해 전체 염색체가 아니라 island+padding 구간에서만 계산. |
| baseline가 0에 가까울 때 fold 계산이 발산(6×10¹¹ 같은 값) | fold 분모를 **min_reads로 clamp**. 즉 저배경 구역에선 `fold_threshold × min_reads`를 넘겨야 호출됨 — 원본 3.1이 걱정한 "baseline≈0 함정"을 수치적으로도 방어. |
| island "gap 기준 분리"(4,5)에서 gap 임계값 | `--island-gap`(기본 50bp)로 노출. |
| core의 최소 길이/미세 dip 처리 | `--min-length`(peak 폭 최소, 기본 10bp), `--max-gap`(core 내부 sub-threshold dip 연결, 기본 10bp) 추가. `--min-length`를 read 길이보다 크게 주면 단일 read-폭 stack을 아티팩트로 걸러낼 수 있음. |
| 실험별 read 길이 QC | `--min-read-length`(bp) 추가 — adapter trimming 후 너무 짧은 read 제거. |
| "strand별 분리"(5-1)에서 라이브러리 방향성 | `--library-type {forward,reverse,unstranded}`. **eCLIP single-end는 대개 reverse-stranded** — 이걸 틀리면 peak가 반대 strand에 찍힘. |

## 3. 원본에 없어서 보충한 항목 (정확도/실용성)

- **Splice-aware coverage.** RNA 기반 CLIP에서 spliced read가 intron을 채우면 안 됨. `read.get_blocks()`로 정렬 블록만 depth에 반영(intron/deletion 제외). 원본 spec에는 명시가 없었으나 RBP CLIP엔 필수.
- **기본 BAM QC 필터.** unmapped/secondary/supplementary/QC-fail 제외, PCR **duplicate 기본 제거**(`--keep-dup`로 유지 가능), `--min-mapq`. CLIP에서 PCR 중복은 흔한 인공물이라 기본 제거가 안전.
- **narrowPeak 컬럼 정의 확정.** signalValue = summit fold, pValue/qValue = -1(통계 미계산, 의도적), peak = summit offset, score = fold의 0–1000 표시용 매핑.
- **bedGraph 출력**(`--bdg`). strand별 RPM 트랙 — 형이 IGV에서 결과를 눈으로 대조하는 워크플로에 직접 도움.
- **결정론적 정렬/이름**(peak_1..N), `--chrom`으로 특정 염색체만.

## 4. MACS3에서 "참고만" 한 것 vs "안 가져온" 것

참고: SPMR(=RPM) 정규화 개념, "국소 배경을 window로 본다"는 아이디어(단, MACS의 mean-λ가 아니라 **저백분위 valley**).

안 가져옴(원본 철학 유지): p/q-value·Poisson/NB 모델, `--slocal/--llocal` 이중 λ, `--extsize/--shift`·read 확장, broad/gapped 모드, subcommand 구조. 즉 옵션 표면을 그대로 이식하지 않고 원본 CLI(§7)에 맞춰 최소화.

## 5. 여전히 스코프 밖 (원본 6절 유지)

Motif/conservation/replicate 재현성 등 생물학적 검증, 통계적 유의성 계산은 이번에도 제외. 추가로 지금 단계에서 **넣지 않은 실무 항목**(추후 후보): blacklist 영역 제외 옵션, 멀티프로세싱(염색체 병렬), paired-end 조각 단위 처리, 실데이터 기반 기본값 튜닝.

## 6. 검증

합성 BAM(넓은 island + 5kb 강한 core + 12kb 약한 decoy bump + 노이즈 + spliced read)으로 end-to-end 확인:
- splice-aware coverage(intron 안 채움), 강한 core만 호출·약한 bump 배제, narrowPeak 컬럼 유효성 — pytest 7개 통과.
- 적정 min_reads에서 CLI가 5kb core 하나만 fold≈18.7로 정확히 호출.

다음 단계는 원본 §9-6대로 **형의 실데이터(ESRP1/RBFOX2/YBX1 BAM)로 파라미터 튜닝 + IGV 대조**.
