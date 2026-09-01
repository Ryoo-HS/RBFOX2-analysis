# ============================================================
# CLIP-seq Peak Calling & Cleaning Functions  (v2.1)
# RBFOX2 / QKI / PTBP1 / ESRP1 공용 파이프라인
#
# 사용법:
#   source("scripts/01_peak_calling.R")
#   summary <- process_all_rbps(bam_paths, repeat_hg19, lnc_gr, sno_mi, numt_gr,
#                                save = TRUE, save_dir = "data")
#
# 필요 패키지: GenomicAlignments, GenomicRanges, IRanges, chipseq
# 필요 객체 (세션에 미리 로드되어 있어야 함):
#   - repeat_hg19 : RepeatMasker GRanges (hg19), $gene_id 컬럼에 repeat class
#   - lnc_gr      : GENCODE v19 lncRNA 계열 GRanges (선택, NULL이면 스킵)
#                   (lincRNA/antisense/processed_transcript/
#                    sense_intronic/sense_overlapping/3prime_overlapping_ncRNA)
#   - sno_mi      : snoRNA/miRNA GRanges
#   - numt_gr     : NUMT(핵 내 미토콘드리아 유래 서열) GRanges (선택, NULL이면 스킵)
#
# 필터 설계 원칙:
#   repeat_hg19 / sno_mi는 항상 적용되는 기본 필터이고, lnc_gr / numt_gr는
#   RBP·분석 목적에 따라 선택적으로 켜고 끄는 필터임 (해당 인자를 NULL로
#   두면 그 단계는 자동으로 스킵되고 로그에 스킵 사실이 남음). 즉 모든
#   RBP·모든 실행에 동일한 필터 세트를 강제하지 않고, 그때그때 목적에
#   맞는 필터만 골라 적용하도록 설계함.
#
# v2.1 변경사항:
#   - noise_filters 리스트 통합 버전(v3)에서 이름 누락 시 필터가 통째로
#     스킵되는 문제가 있어 repeat_hg19 / lnc_gr / sno_mi 개별 인자 방식으로
#     되돌림 (v2 필터링 로직 그대로)
#   - process_all_rbps에 save(T/F) 인자는 유지.
#       save = TRUE  : save_dir에 저장 (save_dir 안 주면 현재 working directory)
#       save = FALSE : 디스크에 안 쓰고 결과 GRanges를 list로 반환
#
# v2 변경사항:
#   - sno_mi(snoRNA/miRNA) 필터링 단계 추가
#   - 함수 내부에서 다 쓴 무거운 중간 객체(reads, coverage 등) 즉시 rm()+gc()
#   - process_all_rbps가 lapply 대신 for-loop + 즉시 저장 방식으로 변경
#     (중간에 메모리 에러 나도 그 전까지 처리된 RBP는 디스크에 남음)
# ============================================================

library(GenomicAlignments)
library(GenomicRanges)
library(chipseq)   # peakSummary() 함수 제공

#' BAM 파일에서 peak을 콜링하고 노이즈(repeat, lncRNA, snoRNA/miRNA, artifact)를 제거
#'
#' @param bam_path    BAM 파일 경로 (문자열)
#' @param repeat_hg19 RepeatMasker GRanges 객체
#' @param lnc_gr      lncRNA GRanges 객체
#' @param sno_mi      snoRNA/miRNA GRanges 객체
#' @param numt_gr     NUMT(핵 내 미토콘드리아 유래 서열) GRanges 객체. NULL이면 스킵
#' @param lower_thresh coverage slice cutoff (기본값 20)
#' @param min_width   width/max artifact 필터의 width 기준 (기본값 26, read length 기준)
#' @param max_signal_for_short_peak  width < min_width인 peak을 살릴 수 있는 max 상한 (기본값 200)
#'
#' @return 정제된 peak GRanges 객체
call_and_clean_peaks <- function(bam_path, repeat_hg19, lnc_gr = NULL, sno_mi, numt_gr = NULL,
                                  lower_thresh = 20,
                                  min_width = 26,
                                  max_signal_for_short_peak = 200) {

  cat("=== BAM 읽는 중:", bam_path, "===\n")
  reads <- readGAlignments(bam_path)

  # 1. 표준 염색체만 (scaffold/alt contig 제거)
  reads <- reads[grepl("chr[1234567890XY]", seqnames(reads))]
  reads <- reads[!grepl("random", seqnames(reads))]
  reads <- reads[!grepl("hap", seqnames(reads))]
  reads <- reads[!grepl("Un", seqnames(reads))]

  # 2. strand별 coverage + peak
  cov_pos <- coverage(reads[strand(reads) == "+"])
  cov_neg <- coverage(reads[strand(reads) == "-"])

  # reads 다 썼으니 바로 정리 (제일 무거운 객체)
  rm(reads)
  gc()

  peak_pos <- peakSummary(slice(cov_pos, lower = lower_thresh))
  peak_neg <- peakSummary(slice(cov_neg, lower = lower_thresh))

  # coverage 객체도 다 썼으니 정리
  rm(cov_pos, cov_neg)
  gc()

  strand(peak_pos) <- "+"
  strand(peak_neg) <- "-"

  peaks_gr <- c(peak_pos, peak_neg)
  peaks_gr <- peaks_gr[order(peaks_gr$max, decreasing = TRUE)]

  rm(peak_pos, peak_neg)
  gc()

  n0 <- length(peaks_gr)
  cat("초기 peak:", n0, "\n")

  # 3. width/max artifact 필터
  #    (폭이 좁은데 signal이 비정상적으로 높은 peak = mapping artifact 또는
  #     PCR 과증폭(PCR duplication/amplification bias)으로 특정 위치에
  #     read가 비정상적으로 쌓인 것으로 의심)
  peaks_gr <- peaks_gr[width(peaks_gr) >= min_width |
                        peaks_gr$max <= max_signal_for_short_peak]
  cat("width/max 필터 후:", length(peaks_gr), "\n")

  # 4. repeat 제거 (5S, 7S, tRNA, rRNA)
  for (pat in c("5S", "7S", "tRNA", "rRNA")) {
    counts <- countOverlaps(peaks_gr, repeat_hg19[grepl(pat, repeat_hg19$gene_id)])
    peaks_gr <- peaks_gr[counts == 0]
    cat(pat, "제거 후:", length(peaks_gr), "\n")
  }

  # 5. lncRNA 제거 (lnc_gr 안 넘기면 스킵)
  if (!is.null(lnc_gr)) {
    counts <- countOverlaps(peaks_gr, lnc_gr)
    peaks_gr <- peaks_gr[counts == 0]
    cat("lncRNA 제거 후:", length(peaks_gr), "\n")
  } else {
    cat("lnc_gr 없음 -> lncRNA 필터링 스킵\n")
  }

  # 6. snoRNA/miRNA 제거
  counts <- countOverlaps(peaks_gr, sno_mi)
  peaks_gr <- peaks_gr[counts == 0]
  cat("snoRNA/miRNA 제거 후:", length(peaks_gr), "\n")

  # 7. NUMT(핵 내 미토콘드리아 유래 서열) 제거 (numt_gr 안 넘기면 스킵)
  if (!is.null(numt_gr)) {
    counts <- countOverlaps(peaks_gr, numt_gr)
    peaks_gr <- peaks_gr[counts == 0]
    cat("NUMT 제거 후:", length(peaks_gr), "\n")
  } else {
    cat("numt_gr 없음 -> NUMT 필터링 스킵\n")
  }

  cat("=== 최종:", n0, "->", length(peaks_gr), "===\n\n")

  rm(counts)
  gc()

  return(peaks_gr)
}


#' 여러 RBP를 순차 처리 + (선택적) 즉시 저장하는 래퍼
#'
#' lapply 대신 for-loop을 써서, 한 RBP 처리가 끝나면 바로 다음으로 넘어감.
#' save = TRUE면 .rds로 즉시 저장 후 메모리에서 비움 (RBP 하나 처리 중
#' 에러가 나도 그 전까지 결과는 디스크에 남음). save = FALSE면 저장 없이
#' 전부 메모리에 들고 있다가 list로 반환함.
#'
#' @param bam_paths   named list, 이름 = RBP명, 값 = BAM 경로
#' @param repeat_hg19 RepeatMasker GRanges
#' @param lnc_gr      lncRNA GRanges
#' @param sno_mi      snoRNA/miRNA GRanges
#' @param numt_gr     NUMT GRanges. NULL이면 스킵
#' @param save        TRUE/FALSE. TRUE면 디스크에 저장, FALSE면 저장 없이
#'                    peak GRanges들을 list로 반환 (기본값 TRUE)
#' @param save_dir    save = TRUE일 때 저장할 디렉토리. NULL이면 현재
#'                    working directory(getwd())에 저장 (기본값 NULL)
#' @param ...         call_and_clean_peaks에 전달할 추가 인자 (lower_thresh 등)
#'
#' @return save = TRUE  : RBP별 최종 peak 개수 요약 (list, invisible)
#'         save = FALSE : RBP별 peak GRanges 객체 (list, invisible)
process_all_rbps <- function(bam_paths, repeat_hg19, lnc_gr = NULL, sno_mi, numt_gr = NULL,
                              save = TRUE, save_dir = NULL, ...) {

  if (isTRUE(save)) {
    if (is.null(save_dir)) save_dir <- getwd()
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    cat("저장 경로:", normalizePath(save_dir), "\n")
  } else {
    cat("save = FALSE -> 디스크 저장 없이 결과를 메모리에 보관함\n")
  }

  results_summary <- list()
  peaks_list <- list()

  for (name in names(bam_paths)) {
    cat("\n##########", name, "시작 ##########\n")

    peaks_gr <- call_and_clean_peaks(bam_paths[[name]], repeat_hg19, lnc_gr, sno_mi, numt_gr, ...)

    results_summary[[name]] <- length(peaks_gr)

    if (isTRUE(save)) {
      saveRDS(peaks_gr, file.path(save_dir, paste0(name, "_peaks.rds")))
      rm(peaks_gr)
      gc()
    } else {
      peaks_list[[name]] <- peaks_gr
    }
  }

  cat("\n=== 전체 완료 ===\n")
  print(unlist(results_summary))

  if (isTRUE(save)) {
    invisible(results_summary)
  } else {
    invisible(peaks_list)
  }
}


# ============================================================
# 사용 예시 (직접 실행하지 않고 참고용으로 남겨둠)
# 경로는 전부 프로젝트 루트 기준 상대경로 예시임 - 본인 환경에 맞게 수정
# ============================================================
#
# repeat_hg19 <- readRDS("data/annotation/repeat_hg19.rds")
# lnc_gr      <- readRDS("data/annotation/lnc_gr.rds")       # 선택 필터 - 이번 분석에서는 미투입 (스킵)
# sno_mi      <- readRDS("data/annotation/sno_mi.rds")
# numt_gr     <- readRDS("data/annotation/numt_gr.rds")      # 선택 필터 - 이번 분석에서는 미투입 (스킵)
#
# bam_paths <- list(
#   RBFOX2 = "data/bam/RBFOX2.sorted.bam",
#   QKI    = "data/bam/QKI.sorted.bam",
#   PTBP1  = "data/bam/PTBP1.sorted.bam",
#   ESRP1  = "data/bam/ESRP1.sorted.bam"
# )
#
# process_all_rbps(bam_paths, repeat_hg19, lnc_gr, sno_mi, numt_gr,
#                   save = TRUE, save_dir = "data/peaks")
#
# # save = FALSE로 쓰면 저장 없이 바로 list로 받을 수 있음
# all_peaks <- process_all_rbps(bam_paths, repeat_hg19, lnc_gr, sno_mi, numt_gr, save = FALSE)
# QKI_peaks_gr <- all_peaks$QKI
#
# # 저장한 경우엔 나중에 개별 .rds로 불러오기
# RBFOX2_peaks_gr <- readRDS("data/peaks/RBFOX2_peaks.rds")
# QKI_peaks_gr    <- readRDS("data/peaks/QKI_peaks.rds")
# PTBP1_peaks_gr  <- readRDS("data/peaks/PTBP1_peaks.rds")
# ESRP1_peaks_gr  <- readRDS("data/peaks/ESRP1_peaks.rds")
