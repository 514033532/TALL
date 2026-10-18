#!/usr/bin/env Rscript
# Authoritative summary + QA for T-ALL manuscript package
# Loads the master Seurat objects and CellChat matrices; writes clean tables.
suppressMessages({library(Seurat)})

base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
out  <- file.path(pkg, "Analysis_QA")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

log <- file(file.path(out, "01_R_summary_log.txt"), open = "wt")
sink(log, type = "output"); sink(log, type = "message")

cat("========== 1. Master object: malignant_tcells_updated_with_cnv.rds ==========\n")
t0 <- Sys.time()
obj <- readRDS(file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/malignant_tcells_updated_with_cnv.rds"))
cat("dim:", dim(obj)[1], "genes x", dim(obj)[2], "cells | load time (min):",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
md <- obj@meta.data
cat("--- metadata columns ---\n"); print(colnames(md))

sample_col <- if ("sample_id" %in% colnames(md)) "sample_id" else if ("orig.ident" %in% colnames(md)) "orig.ident" else NULL
cat("sample column:", sample_col, "\n")

write_summary <- function(x, file) write.csv(x, file.path(out, file), row.names = FALSE)

# per-sample summary
if (!is.null(sample_col) && "InferCNV" %in% colnames(md)) {
  tab_pos  <- table(md[[sample_col]], md$InferCNV)
  if (!"Positive" %in% colnames(tab_pos)) tab_pos <- cbind(tab_pos, Positive = 0)
  if (!"Negative" %in% colnames(tab_pos)) tab_pos <- cbind(tab_pos, Negative = 0)
  ss <- data.frame(sample = rownames(tab_pos),
                   n_cells   = rowSums(tab_pos),
                   n_CNV_Positive = tab_pos[, "Positive"],
                   n_CNV_Negative = tab_pos[, "Negative"],
                   pct_Positive = round(100 * tab_pos[, "Positive"] / rowSums(tab_pos), 2))
  if ("InferCNVcount" %in% colnames(md)) {
    ss$mean_InferCNVcount <- round(tapply(md$InferCNVcount, md[[sample_col]], mean, na.rm = TRUE)[rownames(ss)], 4)
    ss$median_InferCNVcount <- round(tapply(md$InferCNVcount, md[[sample_col]], median, na.rm = TRUE)[rownames(ss)], 4)
  }
  ss <- ss[order(-ss$pct_Positive), ]
  write_summary(ss, "02_per_sample_CNV_summary.csv")
  print(ss)
}

# global bin / status / phase
if ("InferCNV_bin" %in% colnames(md)) {
  write_summary(as.data.frame(table(md$InferCNV_bin, dnn = "Bin")), "03_bin_counts.csv")
  print(table(md$InferCNV_bin))
}
if ("Phase" %in% colnames(md)) {
  write_summary(as.data.frame(table(md$Phase, dnn = "Phase")), "04_phase_counts.csv")
  print(table(md$Phase))
  if ("InferCNV_bin" %in% colnames(md)) {
    ct <- table(md$Phase, md$InferCNV_bin)
    write_summary(as.data.frame.matrix(ct), "05_phase_x_bin_crosstab.csv")
    ct_pct <- prop.table(ct, 2)
    write_summary(round(as.data.frame.matrix(ct_pct) * 100, 2), "06_phase_x_bin_pct_bybin.csv")
    print(ct)
    if (min(dim(ct)) > 1) {
      cs <- chisq.test(ct, simulate.p.value = TRUE, B = 1e5)
      cat("chi2 phase x bin: p =", cs$p.value, "\n")
    }
  }
  if ("InferCNV" %in% colnames(md)) {
    ct2 <- table(md$Phase, md$InferCNV)
    print(ct2)
    cs2 <- chisq.test(ct2, simulate.p.value = TRUE, B = 1e5)
    cat("chi2 phase x CNV status: p =", cs2$p.value, "\n")
  }
}

# S/G2M fraction per bin (proliferation index)
if (all(c("Phase", "InferCNV_bin") %in% colnames(md))) {
  bins <- levels(factor(md$InferCNV_bin))
  prol <- t(sapply(bins, function(b) {
    sub <- md[md$InferCNV_bin == b, ]
    c(n = nrow(sub),
      pct_S = 100 * mean(sub$Phase == "S", na.rm = TRUE),
      pct_G2M = 100 * mean(sub$Phase == "G2M", na.rm = TRUE),
      pct_G1 = 100 * mean(sub$Phase == "G1", na.rm = TRUE))
  }))
  write_summary(cbind(Bin = bins, round(prol, 2)), "07_proliferation_by_bin.csv")
  print(prol)
}

# export full per-cell annotation (compressed)
keep <- intersect(c(sample_col, "InferCNV", "InferCNVcount", "InferCNV_bin", "Phase",
                    "seurat_clusters", "prediction_label", "sample_type", "cell_type"), colnames(md))
write.csv(cbind(cell_id = rownames(md), md[, keep, drop = FALSE]),
          gzfile(file.path(out, "08_cell_annotation_export.csv.gz")), row.names = FALSE)
cat("exported", nrow(md), "cells x", length(keep) + 1, "columns\n")

cat("\n========== 2. Prediction-label object (22,679 cells) ==========\n")
t1 <- Sys.time()
obj2 <- readRDS(file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/marker-v2/combined_TALL_22679_prediction_label_seurat.rds"))
cat("dim:", dim(obj2)[1], "x", dim(obj2)[2], "| load (min):",
    round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "\n")
md2 <- obj2@meta.data
cat("columns:"); print(colnames(md2))
lab_col <- intersect(c("prediction_label", "prediction", "Prediction_Label"), colnames(md2))[1]
if (!is.na(lab_col)) {
  write_summary(as.data.frame(table(md2[[lab_col]], dnn = "label")), "09_prediction_label_counts.csv")
  print(table(md2[[lab_col]]))
  s2_col <- if ("sample_id" %in% colnames(md2)) "sample_id" else "orig.ident"
  if (s2_col %in% colnames(md2)) {
    tab2 <- table(md2[[s2_col]], md2[[lab_col]])
    write_summary(as.data.frame.matrix(tab2), "10_prediction_label_by_sample.csv")
    print(tab2)
  }
}
rm(obj2); gc()

cat("\n========== 3. CellChat RData structure ==========\n")
cc_env <- new.env()
load(file.path(base, "cellchat_matrices/all_cellchat_matrices_summary.RData"), envir = cc_env)
cat("objects:", paste(ls(cc_env), collapse = ", "), "\n")
for (nm in ls(cc_env)) {
  o <- get(nm, envir = cc_env)
  cat("--", nm, "class:", paste(class(o), collapse = "/"), "dim:", paste(dim(o), collapse = "x"), "\n")
}

cat("\n========== 4. Pathway-level CellChat differential (supplement) ==========\n")
# netP matrices: cell-type-pair x pathway. Aggregate to pathway-level strength per sample.
read_netP <- function(f) {
  m <- read.csv(f, check.names = FALSE, row.names = 1)
  if (ncol(m) == 0) m <- read.csv(f, check.names = FALSE)
  m
}
dirs <- list(malignant = file.path(base, "cellchat_matrices/malignant"),
             normal    = file.path(base, "cellchat_matrices/NORMAL"))
sum_list <- list()
for (g in names(dirs)) {
  fs <- list.files(dirs[[g]], pattern = "netP_matrix\\.csv$", full.names = TRUE)
  for (f in fs) {
    s <- sub("_netP_matrix\\.csv", "", basename(f))
    m <- read_netP(f)
    tot <- colSums(as.matrix(m), na.rm = TRUE)          # per-pathway total strength
    sum_list[[paste0(g, "|", s)]] <- data.frame(group = g, sample = s,
                                                 pathway = names(tot), strength = as.numeric(tot))
  }
}
sum_df <- do.call(rbind, sum_list)
write_summary(sum_df, "11_pathway_strength_by_sample.csv")
pathways <- sort(unique(sum_df$pathway))
res <- do.call(rbind, lapply(pathways, function(p) {
  a <- sum_df$strength[sum_df$pathway == p & sum_df$group == "malignant"]
  b <- sum_df$strength[sum_df$pathway == p & sum_df$group == "normal"]
  if (length(a) < 2 || length(b) < 2) return(NULL)
  w <- wilcox.test(a, b)
  data.frame(pathway = p, mean_malignant = mean(a), mean_normal = mean(b),
             log2FC = log2((mean(a) + 1e-6) / (mean(b) + 1e-6)),
             wilcox_p = w$p.value, n_mal = length(a), n_norm = length(b))
}))
if (!is.null(res)) {
  res$p_adj <- p.adjust(res$wilcox_p, "BH")
  res <- res[order(res$p_adj), ]
  write_summary(res, "12_pathway_differential_wilcoxon.csv")
  print(head(res, 30))
}

sink(type = "message"); sink(type = "output")
cat("DONE\n")
