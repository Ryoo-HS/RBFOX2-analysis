# ============================================================
# Step 2: Signal Threshold Filtering (empirical cutoff)
#
# 목적:
#   01_peak_calling.R에서 노이즈 필터링까지 끝난 peak set은
#   RBP마다 raw peak 개수 편차가 매우 큼
#     - QKI    : 116,003 peaks
#     - RBFOX2 : 352,435 peaks
#     - ESRP1  :  54,532 peaks
#   이 상태로 downstream 비교(co-binding 등)를 하면 peak set
#   크기 차이가 결과를 왜곡할 수 있어, 각 RBP의 max signal
#   threshold를 조정해 최종 peak 개수를 비슷한 규모(~2만 개)로
#   맞춤.
#
# 방법론 한계 (README에도 동일하게 명시):
#   이 threshold는 p-value/FDR 등 통계적 유의성 기준이 아니라
#   "RBP 간 peak set 크기를 표준화"하기 위한 경험적(empirical)
#   컷오프임. 값 자체(65 / 30 / 35)는 RBP별로 peak 개수가
#   ~2만 개 선에 오도록 반복 확인하며 정함.
# ============================================================

# --- RBFOX2 -------------------------------------------------
# max > 65  ->  23,652 peaks
RBFOX2_peaks_gr <- RBFOX2_peaks_gr[RBFOX2_peaks_gr$max > 65]

# --- ESRP1 ----------------------------------------------------
# max > 30  ->  20,514 peaks
ESRP1_peaks_gr <- ESRP1_peaks_gr[ESRP1_peaks_gr$max > 30]

# --- QKI --------------------------------------------------------
# max > 35  ->  22,447 peaks
# 참고: QKI는 lncRNA(lnc_gr) 필터링을 01단계에서 적용하지 않음
#       (QKI 전용 lncRNA annotation 미비로 전체 필터링에서 제외)
QKI_peaks_gr <- QKI_peaks_gr[QKI_peaks_gr$max > 35]

cat("RBFOX2:", length(RBFOX2_peaks_gr), "\n")
cat("ESRP1 :", length(ESRP1_peaks_gr), "\n")
cat("QKI   :", length(QKI_peaks_gr), "\n")

# 이후 03_annotation_go_kegg.R에서 사용할 peak_list
peak_list <- list(RBFOX2 = RBFOX2_peaks_gr, ESRP1 = ESRP1_peaks_gr, QKI = QKI_peaks_gr)
