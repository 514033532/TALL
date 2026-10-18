#!/usr/bin/env Rscript
# 05c: Harmonized CellChat (identity = CNV burden bins) — fixed pipeline.
# Aggregates ALL interactions (thresh=NULL) to pathway strength (netP-like), p-value independent.
suppressMessages({library(Seurat); library(CellChat)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
out  <- file.path(pkg, "Analysis_QA")
dir.create(file.path(out, "cellchat_bins"), showWarnings = FALSE, recursive = TRUE)
log <- file(file.path(out, "05c_cellchat_log.txt"), open = "wt")
sink(log, type = "output"); sink(log, type = "message")

t0 <- Sys.time()
obj <- readRDS(file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/marker-v2/combined_TALL_22679_prediction_label_seurat.rds"))
obj <- JoinLayers(obj)
if (!"data" %in% Layers(obj[["RNA"]])) obj <- NormalizeData(obj)
obj$GSM <- regmatches(colnames(obj), regexpr("GSM[0-9]+", colnames(obj)))
obj$bin <- cut(obj$InferCNVcount, breaks = c(0, 0.25, 0.75, 1.5, Inf),
               labels = c("B1_0-0.25", "B2_0.25-0.75", "B3_0.75-1.5", "B4_1.5+"),
               include.lowest = TRUE, right = FALSE)
data.all <- GetAssayData(obj, layer = "data")
cat("loaded in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min; data dim:",
    paste(dim(data.all), collapse = "x"), "\n")

gsms <- sort(unique(obj$GSM))
path_str <- list(); net_long <- list(); sm_list <- list()
for (g in gsms) {
  cells <- colnames(obj)[obj$GSM == g]
  labs <- droplevels(factor(obj$bin[cells]))
  tab <- table(labs); keep <- names(tab)[tab >= 10]
  cells <- cells[labs %in% keep]; labs <- droplevels(labs[labs %in% keep])
  cat("==", g, "groups:", paste(levels(labs), collapse = "/"), "cells:", length(cells), "\n")
  if (length(cells) < 50 || length(unique(labs)) < 2) { cat("  skip\n"); next }
  ok <- tryCatch({
    meta <- data.frame(labels = labs, row.names = cells)
    cc <- createCellChat(object = data.all[, cells, drop = FALSE], meta = meta)
    cc <- setIdent(cc, ident.use = "labels")
    cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, raw.use = TRUE)
    TRUE
  }, error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); FALSE })
  if (!ok) next
  sc <- tryCatch(subsetCommunication(cc, slot.name = "net", thresh = NULL), error = function(e) NULL)
  if (is.null(sc) || nrow(sc) == 0) { cat("  no interactions\n"); next }
  saveRDS(cc, file.path(out, "cellchat_bins", paste0(g, "_cellchat_bins.rds")))
  # pathway-level strength (netP-like, all pairs)
  pw <- aggregate(prob ~ pathway_name + source + target, sc, sum)
  pw$sample <- g
  net_long[[g]] <- pw
  tot <- aggregate(prob ~ pathway_name, pw, sum)
  tot$sample <- g
  path_str[[g]] <- tot
  sm_list[[g]] <- data.frame(sample = g, total_strength = sum(pw$prob),
                             n_pathways = length(unique(sc$pathway_name)),
                             n_pairs = nrow(sc))
  cat("  done; pairs:", nrow(sc), "pathways:", length(unique(sc$pathway_name)),
      "top:", paste(head(tot$pathway_name[order(-tot$prob)], 3), collapse = "/"), "\n")
}
if (length(net_long) > 0) {
  net_df <- do.call(rbind, net_long)
  write.csv(net_df, file.path(out, "33_harmonized_cellchat_netP_long.csv"), row.names = FALSE)
  ps_df <- do.call(rbind, path_str)
  write.csv(ps_df, file.path(out, "34_harmonized_pathway_strength.csv"), row.names = FALSE)
  sm_df <- do.call(rbind, sm_list)
  write.csv(sm_df, file.path(out, "35_harmonized_sample_metrics.csv"), row.names = FALSE)
  print(sm_df)
  gp <- aggregate(prob ~ source + pathway_name, net_df, sum)
  gp <- gp[order(-gp$prob), ]
  write.csv(gp, file.path(out, "36_harmonized_sender_pathway.csv"), row.names = FALSE)
  print(head(gp, 25))
}
sink(type = "message"); sink(type = "output")
cat("DONE\n")
