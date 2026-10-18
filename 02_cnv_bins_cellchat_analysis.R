#!/usr/bin/env Rscript
# Supplement 2: CNV-bin x cell-cycle, subcluster landscape, CellChat network analysis
# All CSV-based (fast). Outputs to Analysis_QA.
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
out  <- file.path(pkg, "Analysis_QA")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
log <- file(file.path(out, "02_analysis_log.txt"), open = "wt")
sink(log, type = "output"); sink(log, type = "message")

cat("========== 1. CNV bins x cell cycle (from annotation export) ==========\n")
md <- read.csv(gzfile(file.path(out, "08_cell_annotation_export.csv.gz")), stringsAsFactors = FALSE)
cat("cells:", nrow(md), "columns:", paste(colnames(md), collapse = ","), "\n")
md$InferCNV_bin <- cut(md$InferCNVcount, breaks = c(0, 0.25, 0.75, 1.5, Inf),
                      labels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"),
                      include.lowest = TRUE, right = FALSE)
bin_tab <- table(md$InferCNV_bin)
write.csv(as.data.frame(bin_tab), file.path(out, "13_bin_counts_recomputed.csv"), row.names = FALSE)
print(bin_tab)

# phase x bin
phxbin <- table(md$Phase, md$InferCNV_bin)
write.csv(as.data.frame.matrix(phxbin), file.path(out, "14_phase_x_bin_counts.csv"))
phxbin_pct <- round(prop.table(phxbin, 2) * 100, 2)
write.csv(as.data.frame.matrix(phxbin_pct), file.path(out, "15_phase_x_bin_pct.csv"))
print(phxbin); print(phxbin_pct)
cs <- chisq.test(phxbin, simulate.p.value = TRUE, B = 1e5)
cat("chi2(phase x bin) p =", cs$p.value, "\n")

prol <- t(sapply(levels(md$InferCNV_bin), function(b) {
  sub <- md[md$InferCNV_bin == b, ]
  c(n = nrow(sub), pct_S = 100 * mean(sub$Phase == "S"), pct_G2M = 100 * mean(sub$Phase == "G2M"),
    pct_G1 = 100 * mean(sub$Phase == "G1"))
}))
write.csv(cbind(Bin = rownames(prol), round(prol, 2)), file.path(out, "16_proliferation_by_bin.csv"), row.names = FALSE)
print(prol)

# phase x CNV status
phxst <- table(md$Phase, md$InferCNV)
print(phxst)
cs2 <- chisq.test(phxst, simulate.p.value = TRUE, B = 1e5)
cat("chi2(phase x CNV status) p =", cs2$p.value, "\n")

# S+G2M fraction in CNV_Positive vs Negative
for (st in c("Positive", "Negative")) {
  sub <- md[md$InferCNV == st, ]
  cat(st, "n =", nrow(sub), "S+G2M% =", round(100 * mean(sub$Phase %in% c("S", "G2M")), 2), "\n")
}

cat("\n========== 2. infercnv subclusters (dC) ==========\n")
meta_csv <- read.csv(file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/combined_cnv_metadata.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
names(meta_csv)[1] <- "cell_id"
cat("meta_csv cells:", nrow(meta_csv), "\n")
if ("infercnv_subcluster" %in% colnames(meta_csv)) {
  sc_tab <- sort(table(meta_csv$infercnv_subcluster), decreasing = TRUE)
  write.csv(as.data.frame(sc_tab), file.path(out, "17_subcluster_counts.csv"), row.names = FALSE)
  print(head(sc_tab, 30))
  # subcluster x sample
  scxsm <- table(meta_csv$infercnv_subcluster, meta_csv$sample_id)
  write.csv(as.data.frame.matrix(scxsm), file.path(out, "18_subcluster_x_sample.csv"))
  # subcluster CNV burden
  agg <- aggregate(InferCNVcount ~ infercnv_subcluster, meta_csv, mean)
  names(agg)[2] <- "mean_InferCNVcount"
  agg$n <- as.integer(table(meta_csv$infercnv_subcluster)[agg$infercnv_subcluster])
  agg <- agg[order(-agg$n), ]
  write.csv(agg, file.path(out, "19_subcluster_burden.csv"), row.names = FALSE)
  print(head(agg, 30))
}

cat("\n========== 3. CellChat: pathway strength per sample (netP long format) ==========\n")
# netP CSVs are mislabeled: actual order = source, target, pathway, probability
read_netP <- function(f) {
  m <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  names(m) <- c("source", "target", "pathway", "prob")
  m
}
samples <- data.frame(
  sample = c("GSM3872440_PRE-T_1","GSM3872441_PRE-T_2","GSM7091996_T1Dx","GSM7091997_T2Dx",
             "GSM7091998_T3Dx","GSM7091999_T4Dx","GSM7092000_T5Dx","GSM7092001_T6Dx",
             "GSM7092002_T6EOI","GSM7092003_T7Dx","GSM7092004_T7EOI","GSM7092005_T8Dx",
             "GSM7092006_T8EOI","GSM7092007_T9Dx","GSM7092008_T9EOI","GSM7092009_T10Dx",
             "GSM7092010_T10EOI","GSM7092011_T11Rel"),
  phase = c("PreT","PreT","Dx","Dx","Dx","Dx","Dx","Dx","EOI","Dx","EOI","Dx","EOI","Dx","EOI","Dx","EOI","Rel"),
  stringsAsFactors = FALSE)
maldir <- file.path(base, "cellchat_matrices/malignant")
path_str <- list()
for (i in seq_len(nrow(samples))) {
  s <- samples$sample[i]
  m <- read_netP(file.path(maldir, paste0(s, "_netP_matrix.csv")))
  tot <- aggregate(prob ~ pathway, m, sum)
  path_str[[s]] <- data.frame(sample = s, phase = samples$phase[i], tot)
}
ps <- do.call(rbind, path_str)
write.csv(ps, file.path(out, "20_pathway_strength_malignant.csv"), row.names = FALSE)
# sample-level metrics
sm <- aggregate(prob ~ sample + phase, ps, sum)
names(sm)[3] <- "total_strength"
n_path <- aggregate(pathway ~ sample + phase, ps, length)
sm$n_pathways <- n_path$pathway
write.csv(sm, file.path(out, "21_sample_communication_metrics.csv"), row.names = FALSE)
print(sm[order(sm$phase), ])
# phase comparison Dx vs EOI per pathway (pathways present in >=3 samples per phase)
pw <- as.data.frame(xtabs(prob ~ pathway + phase, data = ps))
wide <- reshape(pw, idvar = "pathway", timevar = "phase", direction = "wide")
wide[is.na(wide)] <- 0
names(wide) <- sub("Freq\\.", "", names(wide))
wide$diff_EOI_Dx <- wide$EOI - wide$Dx
wide <- wide[order(-wide$Dx), ]
write.csv(wide, file.path(out, "22_pathway_phase_comparison.csv"), row.names = FALSE)
print(head(wide, 25))
# per-pathway wilcox Dx vs EOI (sample-level values)
res_pw <- do.call(rbind, lapply(sort(unique(ps$pathway)), function(p) {
  a <- ps$prob[ps$pathway == p & ps$phase == "Dx"]
  b <- ps$prob[ps$pathway == p & ps$phase == "EOI"]
  if (length(a) < 2 || length(b) < 2) return(NULL)
  w <- wilcox.test(a, b)
  data.frame(pathway = p, mean_Dx = mean(a), mean_EOI = mean(b),
             n_Dx = length(a), n_EOI = length(b), wilcox_p = w$p.value)
}))
res_pw$p_adj <- p.adjust(res_pw$wilcox_p, "BH")
res_pw <- res_pw[order(res_pw$p_adj), ]
write.csv(res_pw, file.path(out, "23_pathway_Dx_vs_EOI_wilcoxon.csv"), row.names = FALSE)
print(res_pw)

cat("\n========== 4. Sender-group contribution (pooled, malignant) ==========\n")
all_netP <- do.call(rbind, lapply(seq_len(nrow(samples)), function(i) {
  m <- read_netP(file.path(maldir, paste0(samples$sample[i], "_netP_matrix.csv")))
  m$sample <- samples$sample[i]; m$phase <- samples$phase[i]; m
}))
sender <- aggregate(prob ~ source, all_netP, sum)
sender <- sender[order(-sender$prob), ]
write.csv(sender, file.path(out, "24_sender_group_total.csv"), row.names = FALSE)
print(head(sender, 20))
# group x pathway matrix (pooled)
gp <- as.data.frame(xtabs(prob ~ source + pathway, data = all_netP))
gp_wide <- reshape(gp, idvar = "source", timevar = "pathway", direction = "wide")
gp_wide[is.na(gp_wide)] <- 0
names(gp_wide) <- sub("Freq\\.", "", names(gp_wide))
write.csv(gp_wide, file.path(out, "25_sender_x_pathway.csv"), row.names = FALSE)
# top sender-pathway combos
gp <- gp[order(-gp$Freq), ]
write.csv(gp, file.path(out, "26_sender_pathway_top.csv"), row.names = FALSE)
print(head(gp, 30))

cat("\n========== 5. Significant LR pairs within T-ALL (lr_pairs, pooled) ==========\n")
lrs <- do.call(rbind, lapply(seq_len(nrow(samples)), function(i) {
  m <- read.csv(file.path(maldir, paste0(samples$sample[i], "_lr_pairs.csv")), stringsAsFactors = FALSE)
  m$sample <- samples$sample[i]; m$phase <- samples$phase[i]; m
}))
lrs_sig <- lrs[lrs$pval < 0.05, ]
cat("total LR rows:", nrow(lrs), "significant (p<0.05):", nrow(lrs_sig), "\n")
agg_lr <- aggregate(prob ~ interaction_name_2 + pathway_name + ligand + receptor + source + target,
                    lrs_sig, FUN = function(x) c(n = length(x), mean_prob = mean(x)))
agg_lr <- do.call(data.frame, agg_lr)
names(agg_lr)[7:8] <- c("n_samples", "mean_prob")
agg_lr <- agg_lr[order(-agg_lr$n_samples, -agg_lr$mean_prob), ]
write.csv(agg_lr, file.path(out, "27_top_significant_LR_pairs.csv"), row.names = FALSE)
print(head(agg_lr, 40))

cat("\n========== 6. Normal BM reference (netP pathway strength) ==========\n")
nordir <- file.path(base, "cellchat_matrices/NORMAL")
nfiles <- list.files(nordir, pattern = "netP_matrix\\.csv$", full.names = TRUE)
norm_ps <- do.call(rbind, lapply(nfiles, function(f) {
  m <- read_netP(f)
  tot <- aggregate(prob ~ pathway, m, sum)
  tot$sample <- sub("_netP_matrix\\.csv", "", basename(f))
  tot
}))
write.csv(norm_ps, file.path(out, "28_pathway_strength_normal.csv"), row.names = FALSE)
print(head(aggregate(prob ~ pathway, norm_ps, mean)[order(-aggregate(prob ~ pathway, norm_ps, mean)$prob), ], 15))

sink(type = "message"); sink(type = "output")
cat("DONE\n")
