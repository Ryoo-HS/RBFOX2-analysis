# RBFOX2 CLIP-seq Read Mapping Pipeline (hg19)

ESRP1 매핑 파이프라인(교수님 제공, bowtie 1.0.1 기준)과 동일한 알고리즘/파라미터를 사용하여
RBFOX2 CLIP-seq raw fastq를 noise filtering 후 hg19 genome에 매핑하는 과정을 기록한다.

## Environment

- bowtie 1.3.0 (conda env `bowtie_env`)
- samtools 1.23
- hg19 bowtie1 index: `~/index/bowtie1/hg19/hg19.*.ebwt`
- Noise reference indexes: `~/index/bowtie1/{rRNA,tRNA,snRNA,snoRNA,miscRNA,hairpin_human}_index`
  - rRNA: NCBI RefSeq rRNA (U13369.1 + 5S/5.8S, 5 sequences)
  - tRNA: GtRNAdb hg19 genomic tRNA (419 sequences)
  - snRNA / snoRNA / miscRNA: GENCODE v19 annotation, `gene_type` 필터링 후 `gffread`로 fasta 추출
  - hairpin: miRBase `hairpin.fa`, human(hsa) 서열만 추출

## Input

```
RBFOX2_trimmed_rc_unique_hg19.fastq   (48,510,106 reads)
```

## Step 1–6: Noise filtering (rRNA / snRNA / snoRNA / tRNA / miscRNA / hairpin)

순차적으로 각 reference에 매핑된 read를 제거(`--un`으로 unmapped만 다음 단계로 전달)한다.
Mismatch 모델은 `-v 1`(end-to-end, mismatch 1개 허용), circular CLIP 라이브러리 특성을 반영하여
5′ 1bp / 3′ 5bp 트리밍(`-5 1 -3 5`)을 적용했다.

```bash
# 1. hairpin (miRBase precursor)
bowtie -v 1 -5 1 -3 5 -x ~/index/bowtie1/hairpin_human_index \
    --un RBFOX2_hg19_hairpin.fastq \
    RBFOX2_trimmed_rc_unique_hg19.fastq > /dev/null

# 2. miscRNA (GENCODE v19 misc_RNA)
bowtie -v 1 -5 1 -3 5 -x ~/index/bowtie1/miscRNA_index \
    --un RBFOX2_hg19_hairpin_miscRNA.fastq \
    RBFOX2_hg19_hairpin.fastq > /dev/null

# 3. rRNA
bowtie -v 1 -5 1 -3 5 -x ~/index/bowtie1/rRNA_index \
    --un RBFOX2_hg19_hairpin_miscRNA_rRNA.fastq \
    RBFOX2_hg19_hairpin_miscRNA.fastq > /dev/null

# 4. snRNA (GENCODE v19)
bowtie -v 1 -5 1 -3 5 -x ~/index/bowtie1/snRNA_index \
    --un RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA.fastq \
    RBFOX2_hg19_hairpin_miscRNA_rRNA.fastq > /dev/null

# 5. snoRNA (GENCODE v19)
bowtie -v 1 -5 1 -3 5 -x ~/index/bowtie1/snoRNA_index \
    --un RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA.fastq \
    RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA.fastq > /dev/null

# 6. tRNA (GtRNAdb)
bowtie -v 1 -5 1 -3 5 -x ~/index/bowtie1/tRNA_index \
    --un RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA.fastq \
    RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA.fastq > /dev/null
```

### Filtering results

| Step | Reference | Reads in | Removed (matched) | Remaining |
|---|---|---|---|---|
| 1 | hairpin | 48,510,106 | 14,725 (0.03%) | 48,495,381 |
| 2 | miscRNA | 48,495,381 | 629,733 (1.30%) | 47,865,648 |
| 3 | rRNA | 47,865,648 | 905,380 (1.89%) | 46,960,268 |
| 4 | snRNA | 46,960,268 | 182,904 (0.39%) | 46,777,364 |
| 5 | snoRNA | 46,777,364 | 96,642 (0.21%) | 46,680,722 |
| 6 | tRNA | 46,680,722 | 26,270 (0.06%) | **46,654,452** |

전체 noise 제거 비율: 1,855,654 / 48,510,106 = **3.83%**

최종 noise-filtered fastq: `RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA.fastq`

## Step 7: Genome mapping (hg19)

ESRP1과 동일한 mismatch 모델(`-n 1 -l 71`, quality-aware seed-based) 및 repeat 처리 옵션
(`-m 5 -a --best --strata`)을 사용했다.

```bash
bowtie -n 1 -m 5 -l 71 -a --best --strata -5 1 -3 5 -S -p 16 \
    -x ~/index/bowtie1/hg19/hg19 \
    RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA.fastq \
    RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA_hg19.sam \
    --un RBFOX2_hg19_genome_n1_unmapped.fastq
```

### Mapping results

```
reads processed: 46,654,452
reads with at least one alignment: 35,481,395 (76.05%)
reads that failed to align: 11,173,057 (23.95%)
reads with alignments suppressed due to -m: 3,788,045 (8.12%)
Reported alignments: 39,807,523
```

> `-a --best --strata` 사용으로 best stratum 내에서 동등한 alignment가 여러 개인 read는
> 여러 줄로 출력될 수 있어, "Reported alignments"(39,807,523)는
> "reads with at least one alignment"(35,481,395)보다 많다.

## Step 8: SAM → BAM → sort → index

```bash
samtools view -bS RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA_hg19.sam \
    > RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA_hg19.bam

samtools sort -@ 6 \
    -o RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA_hg19.sorted.bam \
    RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA_hg19.bam

samtools index RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA_hg19.sorted.bam
```

### Final BAM stats (`samtools flagstat`)

```
50,980,580 + 0 in total (QC-passed reads + QC-failed reads)
50,980,580 + 0 primary
39,807,523 + 0 mapped (78.08%)
39,807,523 + 0 primary mapped (78.08%)
```

> Total record 수(50,980,580)는 mapped alignment lines(39,807,523, 일부 read는
> tie로 다중 라인) + unmapped record lines(11,173,057)의 합과 일치한다.

## Final output

```
RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA_hg19.sorted.bam
RBFOX2_hg19_hairpin_miscRNA_rRNA_snRNA_snoRNA_tRNA_hg19.sorted.bam.bai
```

## Notes

- Noise filtering reference 중 `ATLAS tRNA`(ESRP1 원본 표기)는 정확한 출처를 확인하지 못해
  GtRNAdb genomic tRNA로 대체했다.
- snRNA / snoRNA / miscRNA는 GENCODE v19 `gene_type` 필터링으로 새로 추출한 reference이며,
  ESRP1 측이 사용한 reference와 100% 동일하다는 보장은 없다.
- Genome mapping은 `-v 1`(end-to-end)과 `-n 1 -l 71`(seed-based, quality-aware) 두 가지로
  비교 테스트했으며, 본 데이터셋에서는 두 결과(read 수, mapping rate)가 완전히 동일하게 나왔다.
  ESRP1의 알고리즘과 정확히 일치시키기 위해 `-n 1 -l 71`(ver2) 결과를 최종으로 채택했다.
