#!/usr/bin/env Rscript
# 05: (A) Harmonized CellChat re-run (identity = CNV burden bins) per GSM sample
#     (B) MRD+/- re-analysis of EOI samples (MRD status from Bhasin et al. 2023 Table S1)
suppressMessages({library(Seurat); library(CellChat)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
out  <- file.path(pkg, "Analysis_QA")
dir.create(file.path(out, "cellchat_bins"), showWarnings = FALSE, recursive = TRUE)
log <- file(file.path(out, "05_cellchat_MRD_log.txt"), open = "wt")
sink(log, type = "output"); sink(log, type = "message")

t0 <- Sys.time()
obj <- readRDS(file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/marker-v2/combined_TALL_22679_prediction_label_seurat.rds"))
cat("loaded", dim(obj)[1], "x", dim(obj)[2], "in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

# GSM-level identity from cell names
obj$GSM <- sub("^.*(GSM[0-9]+).*$", "\\1", colnames(obj))
cat("GSMs:\n"); print(table(obj$GSM))

# CNV burden bins (harmonized identity)
obj$bin <- cut(obj$InferCNVcount, breaks = c(0, 0.25, 0.75, 1.5, Inf),
               labels = c("B1_0-0.25", "B2_0.25-0.75", "B3_0.75-1.5", "B4_1.5+"),
               include.lowest = TRUE, right = FALSE)

# ensure normalized data layer
if (!"data" %in% Layers(obj[["RNA"]])) {
  cat("normalizing data...\n"); obj <- NormalizeData(obj)
}
data.all <- GetAssayData(obj, layer = "data")

## ============ Part A: harmonized CellChat per GSM ============
gsms <- sort(unique(obj$GSM))
netP_long <- list()
for (g in gsms) {
  cells <- colnames(obj)[obj$GSM == g]
  labs <- obj$bin[cells]
  tab <- table(labs)
  keep <- names(tab)[tab >= 10]
  cells <- cells[labs %in% keep]
  labs <- droplevels(factor(labs[cells]))
  cat("== ", g, ": cells", length(cells), "bins", paste(keep, collapse = "/"), "\n")
  if (length(cells) < 50 || length(unique(labs)) < 2) { cat("  skip\n"); next }
  tryCatch({
    data.input <- data.all[, cells, drop = FALSE]
    meta <- data.frame(labels = labs, row.names = cells)
    cc <- createCellChat(object = data.input, meta = meta)
    cc <- setIdent(cc, ident.use = "labels")
    cc@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, raw.use = TRUE)
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    # save long-format netP
    np <- subsetCommunication(cc, slot.name = "netP")
    np$sample <- g
    netP_long[[g]] <- np[, c("sample", "source", "target", "pathway", "prob")]
    saveRDS(cc, file.path(out, "cellchat_bins", paste0(g, "_cellchat_bins.rds")))
    cat("  done; interactions:", nrow(subsetCommunication(cc)), "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}
if (length(netP_long) > 0) {
  netP_df <- do.call(rbind, netP_long)
  write.csv(netP_df, file.path(out, "33_harmonized_cellchat_netP_long.csv"), row.names = FALSE)
  # pathway strength per sample
  pw <- aggregate(prob ~ sample + pathway, netP_df, sum)
  write.csv(pw, file.path(out, "34_harmonized_pathway_strength.csv"), row.names = FALSE)
  sm <- aggregate(prob ~ sample, netP_df, sum); names(sm)[2] <- "total_strength"
  sm$n_pathways <- as.integer(table(netP_df$sample)[sm$sample])
  write.csv(sm, file.path(out, "35_harmonized_sample_metrics.csv"), row.names = FALSE)
  print(sm)
  # sender x pathway by bin
  gp <- aggregate(prob ~ source + pathway, netP_df, sum)
  gp <- gp[order(-gp$prob), ]
  write.csv(gp, file.path(out, "36_harmonized_sender_pathway.csv"), row.names = FALSE)
  print(head(gp, 25))
}

## ============ Part B: MRD +/- re-analysis (EOI samples) ============
mrd_map <- c("GSM7092002" = "MRD+ (5.4%)", "GSM7092004" = "MRD- (0%)",
             "GSM7092006" = "MRD+ (8%)",  "GSM7092008" = "MRD+ (5%)",
             "GSM7092010" = "MRD- (0%)")
eoi <- obj[, obj$GSM %in% names(mrd_map)]
md_e <- eoi@meta.data
md_e$MRD <- unname(mrd_map[md_e$GSM])
md_e$MRD_grp <- ifelse(grepl("MRD\\+", md_e$MRD), "MRD_pos", "MRD_neg")
per_s <- do.call(rbind, lapply(sort(unique(md_e$GSM)), function(g) {
  sub <- md_e[md_e$GSM == g, ]
  data.frame(GSM = g, MRD = unique(sub$MRD), MRD_grp = unique(sub$MRD_grp),
             n_cells = nrow(sub),
             pct_CNV_Pos = 100 * mean(sub$InferCNV == "Positive"),
             mean_CNVcount = mean(sub$InferCNVcount),
             pct_SG2M = 100 * mean(sub$Phase %in% c("S", "G2M")),
             n_ML_mag = sum(sub$prediction_label == "ML_maglinant"),
             n_ML_norm = sum(sub$prediction_label == "ML_normal"))
}))
write.csv(per_s, file.path(out, "31_MRD_EOI_per_sample.csv"), row.names = FALSE)
print(per_s)

# gene expression per GSM (EOI)
genes <- c("PDCD1", "HAVCR2", "MIF", "LGALS9", "CD74", "CD44", "CXCR4", "SOX4", "TYMS", "MKI67", "UHRF1", "DHFR", "CCL5", "PRF1")
genes <- genes[genes %in% rownames(eoi)]
cts <- GetAssayData(eoi, layer = "counts")[genes, , drop = FALSE]
lib <- colSums(GetAssayData(eoi, layer = "counts"))
norm <- t(t(cts) / lib * 1e4)
gexp <- do.call(rbind, lapply(sort(unique(md_e$GSM)), function(g) {
  cc <- rownames(md_e)[md_e$GSM == g]
  m <- norm[, cc, drop = FALSE]
  r <- data.frame(GSM = g, MRD = unique(md_e$MRD[md_e$GSM == g]), MRD_grp = unique(md_e$MRD_grp[md_e$GSM == g]))
  for (gg in genes) {
    r[[paste0("mean_", gg)]] <- mean(m[gg, ])
    r[[paste0("pct_", gg)]] <- 100 * mean(m[gg, ] > 0)
  }
  r
}))
write.csv(gexp, file.path(out, "32_MRD_EOI_gene_expression.csv"), row.names = FALSE)
print(gexp[, 1:12])

# group comparison MRD+ vs MRD-
gcmp <- do.call(rbind, lapply(c("pct_CNV_Pos", "mean_CNVcount", "pct_SG2M", "n_ML_mag", "n_ML_norm"), function(v) {
  a <- per_s[[v]][per_s$MRD_grp == "MRD_pos"]; b <- per_s[[v]][per_s$MRD_grp == "MRD_neg"]
  data.frame(metric = v, mean_MRDpos = mean(a), mean_MRDneg = mean(b),
             MRDpos_values = paste(a, collapse = "; "), MRDneg_values = paste(b, collapse = "; "),
             wilcox_p = if (length(a) >= 2 && length(b) >= 2) tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA) else NA)
}))
write.csv(gcmp, file.path(out, "37_MRD_group_comparison.csv"), row.names = FALSE)
print(gcmp)

sink(type = "message"); sink(type = "output")
cat("DONE\n")
