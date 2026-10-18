#!/usr/bin/env Rscript
# 11: Publication-grade "Top 50 differential LR pairs" heatmap.
#     Direct source: lr_intensity_matrix (26 samples x 57 LR pairs) from the original
#     differential analysis; rows ranked by |limma t| (top 50); row z-score clipped +/-2;
#     annotations: sample Group + per-pair LogFC and Significance.
suppressMessages({library(pheatmap)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
figd <- file.path(pkg, "Figures")
e <- new.env()
load(file.path(base, "cell_communication_analysis_results/lr_intensity_matrix.RData"), envir = e)
M <- e$lr_intensity_matrix            # samples x pairs
ann_col_orig <- e$annotation_col      # samples x Group
cat("M dims:", nrow(M), "x", ncol(M), "\n")
cat("pair cols head:", paste(head(colnames(M), 6), collapse = " | "), "\n")
cat("Group table:\n"); print(table(ann_col_orig[[1]]))

diff <- read.csv(file.path(base, "cellchat_matrices/differential_lr_analysis/tables/differential_lr_pairs.csv"),
                 stringsAsFactors = FALSE)
# pair-level stats: aggregate t/logFC/significant by lr_pair_id (any source-target)
pstats <- aggregate(cbind(t = diff$t, logFC = diff$logFC) ~ lr_pair_id, diff, FUN = function(x) x[which.max(abs(x))])
sig <- aggregate(significant ~ lr_pair_id, diff, function(x) names(sort(table(x), decreasing = TRUE))[1])
pstats <- merge(pstats, sig, by = "lr_pair_id")
pstats <- pstats[order(-abs(pstats$t)), ]
pairs_top <- intersect(pstats$lr_pair_id, rownames(M))
pstats <- pstats[match(pairs_top, pstats$lr_pair_id), ]
pairs_top50 <- head(pstats$lr_pair_id, 50)
pstats <- pstats[match(pairs_top50, pstats$lr_pair_id), ]

X <- M[pairs_top50, , drop = FALSE]  # pairs x samples
X[!is.finite(X)] <- 0
rz <- t(scale(t(X)))
rz[!is.finite(rz)] <- 0
rz[rz > 2] <- 2; rz[rz < -2] <- -2

# corrected grouping (fixes original mislabels: T1Dx & PRE-T are malignant, not normal)
normal_set <- c("CD45 BM", "GSM3872442_PBMMC_1", "GSM3872443_PBMMC_2", "GSM3872444_PBMMC_3",
                "GSM4664009_0064", "GSM4664010_3958", "GSM4664011_5286", "GSM4664012_7903")
is_normal <- vapply(colnames(rz), function(s) s %in% normal_set, logical(1))
grp_s <- ifelse(is_normal, "Normal", "Malignant")
ord <- order(factor(grp_s, levels = c("Malignant", "Normal")))
rz <- rz[, ord, drop = FALSE]
grp_s <- grp_s[ord]
ann_col <- data.frame(Group = factor(grp_s, levels = c("Malignant", "Normal")), row.names = colnames(rz))

ann_row <- data.frame(LogFC = pstats$logFC,
                      Significance = factor(ifelse(pstats$significant %in% c("Down", "Up"), pstats$significant, "NS"),
                                            levels = c("Down", "NS", "Up")),
                      row.names = rownames(rz))
ann_cols <- list(Group = c("Malignant" = "#e34948", "Normal" = "#2a78d6"),
                 Significance = c("Down" = "#2a78d6", "NS" = "#c3c2b7", "Up" = "#e34948"))
breaks <- c(seq(-2, 0, length.out = 25), seq(0.05, 2, length.out = 25))
pal <- colorRampPalette(c("#104281", "#5598e7", "#fcfcfb", "#eb6834", "#8c2d1a"))(49)

lab_rows <- sub("_CD74_CD44", "–(CD74+CD44)", rownames(rz))
lab_rows <- sub("_CD74_CXCR4", "–(CD74+CXCR4)", lab_rows)
lab_rows <- sub("_PTPRC", "–CD45", lab_rows)

print(table(grp_s, useNA = "always"))
gaps_val <- sum(grp_s == "Malignant", na.rm = TRUE)
if (!is.finite(gaps_val) || gaps_val < 1 || gaps_val >= ncol(rz)) gaps_val <- NULL
cat("gaps_val:", gaps_val, " ncol:", ncol(rz), "
")
cairo_pdf(file.path(figd, "SuppFig28_Top50_Differential_LR_Heatmap.pdf"), width = 10.5, height = 8.5, family = "Arial")
pheatmap(rz, cluster_rows = FALSE, cluster_cols = FALSE,
         color = pal, breaks = breaks,
         annotation_col = ann_col, annotation_row = ann_row,
         annotation_colors = ann_cols,
         labels_row = lab_rows, labels_col = gsub("GSM[0-9]+_", "", colnames(rz)),
         fontsize = 6.5, fontsize_row = 6, fontsize_col = 6, fontfamily = "Arial",
         angle_col = 45, border_color = "white",
         cellwidth = 9.5, cellheight = 7.5,
         legend_breaks = c(-2, -1, 0, 1, 2), legend_labels = c("≤-2", "-1", "0", "1", "≥2"),
         gaps_col = gaps_val,
         main = "Top 50 differential LR pairs (row-scaled pair probability)")
dev.off()
system2("pdftoppm", c("-png", "-r", "300", "-singlefile",
                      file.path(figd, "SuppFig28_Top50_Differential_LR_Heatmap.pdf"),
                      file.path(figd, "SuppFig28_Top50_Differential_LR_Heatmap")))
cat("SuppFig28 written; rows:", nrow(rz), "cols:", ncol(rz), "\n")
