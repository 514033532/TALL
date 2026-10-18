#!/usr/bin/env Rscript
# Probe: verify CellChat works with harmonized bins after data fix; test 2-group behavior
suppressMessages({library(Seurat); library(CellChat)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
obj <- readRDS(file.path(base,"cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/marker-v2/combined_TALL_22679_prediction_label_seurat.rds"))
obj <- JoinLayers(obj)
if (!"data" %in% Layers(obj[["RNA"]])) obj <- NormalizeData(obj)
obj$GSM <- regmatches(colnames(obj), regexpr("GSM[0-9]+", colnames(obj)))
obj$bin <- cut(obj$InferCNVcount, breaks = c(0, 0.25, 0.75, 1.5, Inf),
               labels = c("B1", "B2", "B3", "B4"), include.lowest = TRUE, right = FALSE)
data.all <- GetAssayData(obj, layer = "data")
cat("data dim:", paste(dim(data.all), collapse = "x"), "\n")

run_cc <- function(g) {
  cells <- colnames(obj)[obj$GSM == g]
  labs <- droplevels(factor(obj$bin[cells]))
  tab <- table(labs); keep <- names(tab)[tab >= 10]
  cells <- cells[labs %in% keep]; labs <- droplevels(labs[labs %in% keep])
  cat("==", g, "groups:", paste(levels(labs), collapse = "/"), "cells:", length(cells), "\n")
  meta <- data.frame(labels = labs, row.names = cells)
  cc <- createCellChat(object = data.all[, cells, drop = FALSE], meta = meta)
  cc <- setIdent(cc, ident.use = "labels")
  cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  r <- tryCatch({cc <- computeCommunProb(cc, raw.use = TRUE); "computeOK"},
                error = function(e) paste("computeERR:", conditionMessage(e)))
  cat("  ->", r, "\n")
  if (r == "computeOK") {
    sc <- tryCatch(subsetCommunication(cc), error = function(e) paste("subERR:", conditionMessage(e)))
    if (is.character(sc)) { cat("  subset:", sc, "\n"); return(NULL) }
    cat("  net interactions:", nrow(sc), " cols:", paste(colnames(sc), collapse = ","), "\n")
    pw <- aggregate(prob ~ pathway_name + source + target, sc, sum)
    print(head(pw[order(-pw$prob), ], 8))
    return(pw)
  }
  NULL
}
pw3 <- run_cc("GSM7091996")   # 3 bins
pw2 <- run_cc("GSM3872441")   # 2 bins
cat("PROBE DONE\n")
