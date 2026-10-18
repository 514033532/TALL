#!/usr/bin/env Rscript
# Lineage check: T-cell vs NK markers across prediction labels (Table 10)
suppressMessages(library(Seurat))
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
obj <- readRDS(file.path(base,"cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/marker-v2/combined_TALL_22679_prediction_label_seurat.rds"))
md <- obj@meta.data
genes <- c("CD3D","CD3E","TRAC","CD4","CD8A","NCAM1","FCGR3A","KLRD1","KLRB1","NKG7","GNLY","PRF1","GZMA",
           "MIF","CD74","CD44","CXCR4","LGALS9","HAVCR2","PTPRC","TYMS","MKI67","UHRF1","SOX4")
genes <- genes[genes %in% rownames(obj)]
cts <- GetAssayData(obj, layer="counts")[genes, , drop=FALSE]
lib <- colSums(GetAssayData(obj, layer="counts"))
norm <- t(t(cts)/lib*1e4)
agg <- do.call(rbind, lapply(sort(unique(md$prediction_label)), function(L){
  cc <- rownames(md)[md$prediction_label==L]
  m <- norm[, cc, drop=FALSE]
  data.frame(label=L, gene=genes, mean_expr=as.numeric(rowMeans(m)), pct=as.numeric(rowMeans(m>0)*100))
}))
wide <- reshape(agg, idvar="gene", timevar="label", direction="wide")
out <- file.path(base,"TALL_manuscript_package/Analysis_QA/30_lineage_markers_by_label.csv")
write.csv(wide, out, row.names=FALSE)
cat("written:", out, "\n")
