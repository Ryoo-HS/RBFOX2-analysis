########################################################################
# CLIP-seq peak 자동화 파이프라인: ChIPseeker annotation + region별 GO/KEGG + ORFik metagene plot
#
# 입력: peak_list (named list of GRanges) 하나뿐. 예)
#   peak_list <- list(RBFOX2 = RBFOX2_peaks_gr, ESRP1 = ESRP1_peaks_gr, QKI = QKI_peaks_gr)
#
# 필요한 annotation 객체(txdb, exon_range, intron_range, three_range, five_range,
# leaders, cds, trailers, universe_entrez)는 세션에 없으면
# run_annotation_GO_KEGG() 실행 시 setup_annotation_env()가 자동으로 만들어서
# .GlobalEnv에 올려줌 (아래 SETUP 섹션 참고). 이미 있는 건 다시 안 만들고 스킵함.
#
# 생성 객체 (RBP 이름별, .GlobalEnv에 assign):
#   <RBP>_anno  : ChIPseeker annotatePeak 결과
#   <RBP>_go    : list(five_utr=, exon=, intron=, three_utr=)
#                 각각 list(BP=, CC=, MF=, BP_universe=) enrichGO 결과
#                 BP/CC/MF = universe 없이 (기본), BP_universe = universe 적용한 BP (비교용)
#   <RBP>_kegg  : list(five_utr=, exon=, intron=, three_utr=) 각각 enrichKEGG 결과 (universe 없음)
#
# GO/KEGG는 plot 없이 객체만 생성. ORFik만 title 달아서 즉시 plot 출력.
########################################################################

## ---------------------------------------------------------------------
## 라이브러리 (파일 source 시 바로 로드됨)
## ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(GenomicFeatures)
  library(ChIPseeker)
  library(ORFik)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

## ---------------------------------------------------------------------
## SETUP: annotation 객체 자동 생성
##
## run_annotation_GO_KEGG()가 필요한 객체를 세션에서 못 찾으면 이 함수를
## 자동으로 호출함. 이미 있는 객체는 다시 안 만들고 스킵하므로, 필요하면
## 직접 미리 한 번 불러도 되고(setup_annotation_env()), 그냥 냅둬도
## run_annotation_GO_KEGG() 첫 호출 때 알아서 채워짐.
##
## universe_entrez는 "TPM >= 1 발현 유전자의 Entrez ID 벡터(unique 처리 끝난
## 최종 벡터)"를 저장해둔 .rds 파일에서 읽어옴. 기본 경로가 안 맞으면
## setup_annotation_env(universe_rds = "다른/경로.rds")로 직접 넘기면 됨.
##
## 기본값은 프로젝트 루트 기준 상대경로 예시임 - 본인 환경에 맞게 수정
## ---------------------------------------------------------------------
setup_annotation_env <- function(
    universe_rds  = "data/annotation/universe_gene.rds",
    min_five_utr  = 100,
    min_cds       = 100,
    min_three_utr = 100,
    force = FALSE) {

  cat("=== annotation 환경 세팅 시작 ===\n")

  make_if_missing <- function(name, build_fun) {
    if (force || !exists(name, envir = .GlobalEnv, inherits = FALSE)) {
      assign(name, build_fun(), envir = .GlobalEnv)
      cat("  -", name, ": 생성 완료\n")
    } else {
      cat("  -", name, ": 이미 있음 (스킵)\n")
    }
  }

  make_if_missing("txdb", function() TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- get("txdb", envir = .GlobalEnv)

  make_if_missing("exon_range",   function() exonsBy(txdb, by = "gene"))
  make_if_missing("intron_range", function() intronsByTranscript(txdb))
  make_if_missing("three_range",  function() threeUTRsByTranscript(txdb))
  make_if_missing("five_range",   function() fiveUTRsByTranscript(txdb))

  if (force || !exists("txNames", envir = .GlobalEnv, inherits = FALSE)) {
    txNames <- filterTranscripts(txdb, min_five_utr, min_cds, min_three_utr)
    assign("txNames", txNames, envir = .GlobalEnv)
    cat("  - txNames : 생성 완료 (", length(txNames), "개 transcript)\n")
  } else {
    txNames <- get("txNames", envir = .GlobalEnv)
    cat("  - txNames : 이미 있음 (스킵)\n")
  }

  region_objs_ready <- all(sapply(c("leaders", "cds", "trailers"),
                                   exists, envir = .GlobalEnv, inherits = FALSE))
  if (force || !region_objs_ready) {
    loadRegions(txdb, parts = c("leaders", "cds", "trailers"),
                names.keep = txNames, envir = .GlobalEnv)
    cat("  - leaders / cds / trailers : 생성 완료\n")
  } else {
    cat("  - leaders / cds / trailers : 이미 있음 (스킵)\n")
  }

  make_if_missing("universe_entrez", function() {
    if (!file.exists(universe_rds)) {
      stop("universe_entrez용 rds 파일을 못 찾음: ", universe_rds,
           "\n  -- setup_annotation_env(universe_rds = \"실제경로.rds\")로 다시 호출하세요.")
    }
    readRDS(universe_rds)
  })

  cat("=== annotation 환경 세팅 끝 ===\n\n")

  invisible(NULL)
}

## ---------------------------------------------------------------------
## 내부 함수: transcript 단위 region(GRangesList by TXID)에서 peak overlap → Entrez ID
## ---------------------------------------------------------------------
.get_entrez_from_tx_region <- function(region_gr, peaks_gr, txdb) {
  counts <- countOverlaps(region_gr, peaks_gr)
  counts <- counts[counts > 0]
  if (length(counts) == 0) return(character(0))
  gene_id <- AnnotationDbi::select(txdb, keys = names(counts),
                                    keytype = "TXID", columns = c("TXID", "GENEID"))
  unique(na.omit(gene_id$GENEID))
}

## ---------------------------------------------------------------------
## 내부 함수: exon_range(by gene)는 names()가 이미 Entrez ID라 바로 뽑음
## ---------------------------------------------------------------------
.get_entrez_from_exon_region <- function(region_gr, peaks_gr) {
  counts <- countOverlaps(region_gr, peaks_gr)
  counts <- counts[counts > 0]
  unique(names(counts))
}

## ---------------------------------------------------------------------
## 내부 함수: 한 유전자 set에 대해 BP/CC/MF enrichGO 세 개 돌리기 (기본 = universe 없음)
## 추가로 universe 적용한 BP 하나를 BP_universe로 같이 넣어서 비교 가능하게 함
## ---------------------------------------------------------------------
.run_go_bp_cc_mf <- function(entrez_ids, universe_entrez_ids) {
  if (length(entrez_ids) == 0) return(list(BP = NULL, CC = NULL, MF = NULL, BP_universe = NULL))
  run_one <- function(ont, universe = NULL) {
    tryCatch(
      enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
               ont = ont, universe = universe,
               pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.1,
               minGSSize = 10, maxGSSize = 500, readable = TRUE, pool = FALSE),
      error = function(e) NULL
    )
  }
  list(
    BP          = run_one("BP"),                                    # 기본: universe 없음
    CC          = run_one("CC"),                                    # 기본: universe 없음
    MF          = run_one("MF"),                                    # 기본: universe 없음
    BP_universe = run_one("BP", universe = universe_entrez_ids)      # 비교용: universe 적용
  )
}

## ---------------------------------------------------------------------
## 내부 함수: 한 유전자 set에 대해 KEGG (기본 = universe 없음)
## ---------------------------------------------------------------------
.run_kegg <- function(entrez_ids, universe_entrez_ids) {
  if (length(entrez_ids) == 0) return(NULL)
  tryCatch(
    enrichKEGG(gene = entrez_ids, organism = "hsa",
               pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.1),
    error = function(e) NULL
  )
}

## ---------------------------------------------------------------------
## 메인 함수
## 입력은 peak_list 하나뿐. txdb, exon_range, intron_range, three_range, five_range,
## leaders, cds, trailers, universe_entrez가 세션(전역 환경)에 없으면
## setup_annotation_env()를 자동으로 호출해서 만든 다음 진행함.
## ---------------------------------------------------------------------
run_annotation_GO_KEGG <- function(peak_list, tssRegion = c(-3000, 3000)) {

  ## 필요한 전역 객체가 세션에 다 있는지 확인, 없으면 자동으로 만듦
  required_objs <- c("txdb", "exon_range", "intron_range", "three_range", "five_range",
                      "leaders", "cds", "trailers", "universe_entrez")
  missing_objs <- required_objs[!sapply(required_objs, exists, envir = .GlobalEnv)]
  if (length(missing_objs) > 0) {
    cat("세션에 없는 객체 감지:", paste(missing_objs, collapse = ", "), "-> 자동 생성 중...\n")
    setup_annotation_env()

    still_missing <- required_objs[!sapply(required_objs, exists, envir = .GlobalEnv)]
    if (length(still_missing) > 0) {
      stop("자동 생성 후에도 없는 객체: ", paste(still_missing, collapse = ", "),
           " -- setup_annotation_env()를 직접 호출해서 원인을 확인하세요.")
    }
  }

  region_list <- list(
    five_utr  = five_range,
    exon      = exon_range,   # 얘만 gene-level (names = Entrez), 나머지는 transcript-level (names = TXID)
    intron    = intron_range,
    three_utr = three_range
  )

  for (rbp in names(peak_list)) {

    cat("=====", rbp, "=====\n")
    peaks_gr <- peak_list[[rbp]]

    ## 1) ChIPseeker annotation -----------------------------------------
    anno <- annotatePeak(peaks_gr, TxDb = txdb, tssRegion = tssRegion, verbose = FALSE)
    assign(paste0(rbp, "_anno"), anno, envir = .GlobalEnv)
    cat("  - annotation 완료\n")

    ## 2) region별 GO / KEGG ---------------------------------------------
    go_list   <- list()
    kegg_list <- list()

    for (region_name in names(region_list)) {
      region_gr <- region_list[[region_name]]

      entrez_ids <- if (region_name == "exon") {
        .get_entrez_from_exon_region(region_gr, peaks_gr)
      } else {
        .get_entrez_from_tx_region(region_gr, peaks_gr, txdb)
      }

      go_list[[region_name]]   <- .run_go_bp_cc_mf(entrez_ids, universe_entrez)
      kegg_list[[region_name]] <- .run_kegg(entrez_ids, universe_entrez)

      cat("  -", region_name, ": gene 수 =", length(entrez_ids), "\n")
    }

    assign(paste0(rbp, "_go"),   go_list,   envir = .GlobalEnv)
    assign(paste0(rbp, "_kegg"), kegg_list, envir = .GlobalEnv)

    ## 3) ORFik metagene plot (title 포함, 바로 출력) ---------------------
    leaderCov  <- metaWindow(peaks_gr, leaders,  scoring = NULL, feature = "leaders")
    cdsCov     <- metaWindow(peaks_gr, cds,      scoring = NULL, feature = "cds")
    trailerCov <- metaWindow(peaks_gr, trailers, scoring = NULL, feature = "trailers")

    df <- data.table::rbindlist(list(leaderCov, cdsCov, trailerCov))
    df[, `:=`(fraction = rbp)]

    print(windowCoveragePlot(df, scoring = "zscore", title = paste0(rbp, "_ORFik_metagene")))
    cat("  - ORFik plot 출력 완료\n")
  }

  invisible(NULL)
}

## ---------------------------------------------------------------------
## 실행 예시
## ---------------------------------------------------------------------
# peak_list <- list(RBFOX2 = RBFOX2_peaks_gr, ESRP1 = ESRP1_peaks_gr, QKI = QKI_peaks_gr)
#
# # 세션에 txdb/exon_range/... 등이 하나도 없어도 바로 실행 가능함.
# # 첫 호출 때 setup_annotation_env()가 자동으로 돌면서 필요한 객체를 다 만듦
# # (exonsBy/intronsByTranscript 등은 데이터가 커서 시간이 좀 걸릴 수 있음).
# run_annotation_GO_KEGG(peak_list)
#
# # universe_entrez rds 경로가 기본값과 다르면 미리 한 번 이렇게 불러두면 됨
# # setup_annotation_env(universe_rds = "다른/경로/universe_gene.rds")
# # run_annotation_GO_KEGG(peak_list)
#
# # 결과 예: RBFOX2_anno, RBFOX2_go$exon$BP, RBFOX2_go$three_utr$MF, RBFOX2_kegg$intron 등
