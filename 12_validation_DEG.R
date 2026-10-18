#!/usr/bin/env Rscript
# 12: Validation-plan items D, E, G (pure computation):
#   D: CDKN2A/CDKN2B expression across CNV bins (chr9p21 context)
#   E: CNV burden persists within G1 (non-cycling) cells -> not a pure S-phase artifact
#   G: correlation of MIF/LGALS9/CD74/HAVCR2 expression with InferCNVcount per cell
suppressMessages({library(Seurat)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
qa   <- file.path(pkg, "Analysis_QA")
log <- file(file.path(qa, "12_validation_DEG_log.txt"), open = "wt")
sink(log, type = "output"); sink(log, type = "message")

obj <- readRDS(file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/marker-v2/combined_TALL_22679_prediction_label_seurat.rds"))
md <- obj@meta.data
genes <- c("CDKN2A", "CDKN2B", "MIF", "LGALS9", "CD74", "HAVCR2", "TYMS", "DHFR", "UHRF1", "MKI67")
genes <- genes[genes %in% rownames(obj)]
cts <- GetAssayData(obj, layer = "counts")[genes, , drop = FALSE]
lib <- colSums(GetAssayData(obj, layer = "counts"))
norm <- t(t(cts) / lib * 1e4)
md$bin <- cut(md$InferCNVcount, breaks = c(0, 0.25, 0.75, 1.5, Inf),
              labels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"),
              include.lowest = TRUE, right = FALSE)

# ---- D: CDKN2A/B across bins ----
D <- do.call(rbind, lapply(levels(md$bin), function(b) {
  cc <- rownames(md)[md$bin == b]
  m <- norm[, cc, drop = FALSE]
  r <- data.frame(bin = b, n = length(cc))
  for (g in genes) { r[[paste0("mean_", g)]] <- mean(m[g, ]); r[[paste0("pct_", g)]] <- 100 * mean(m[g, ] > 0) }
  r
}))
write.csv(D, file.path(qa, "40_gene_expression_by_bin.csv"), row.names = FALSE)
show_cols <- c("bin", "n", paste0("pct_", intersect(c("CDKN2A", "CDKN2B", "MIF", "LGALS9", "TYMS"), genes)))
Dshow <- D[show_cols]; Dshow[, -(1:2)] <- round(Dshow[, -(1:2)], 2); print(Dshow)
cat("genes found:", paste(genes, collapse = ","), "
")

# ---- E: burden within G1 cells ----
G1 <- md[md$Phase == "G1", ]
E <- do.call(rbind, lapply(levels(G1$bin), function(b) {
  sub <- G1[G1$bin == b, ]
  data.frame(bin = b, n_G1 = nrow(sub), mean_CNVcount_in_G1 = mean(sub$InferCNVcount))
}))
write.csv(E, file.path(qa, "41_CNVcount_within_G1_by_bin.csv"), row.names = FALSE)
print(E)
cat("G1-only cells:", nrow(G1), "of", nrow(md), "\n")

# ---- G: per-cell correlation CNVcount vs gene expression ----
G <- data.frame(gene = character(), spearman_rho = numeric(), p = numeric())
for (g in intersect(c("MIF", "LGALS9", "CD74", "HAVCR2", "TYMS", "DHFR", "UHRF1", "MKI67"), genes)) {
  ct <- suppressWarnings(cor.test(md$InferCNVcount, norm[g, ], method = "spearman"))
  G <- rbind(G, data.frame(gene = g, spearman_rho = unname(ct$estimate), p = ct$p.value))
}
write.csv(G, file.path(qa, "42_CNVcount_gene_correlation.csv"), row.names = FALSE)
print(G)

sink(type = "message"); sink(type = "output")
cat("DONE\n")
