#!/usr/bin/env Rscript
# 09: Revised main figures per user feedback (Rev2)
#   Fig3: A + C only (B dropped; S+G2M annotated on A), widths 1:1.4, wrapped GO terms, full axis
#   Fig4: A model colors in legend; B density plot; C only 1020 reclassified cells, sorted desc
#   Fig6: longitudinal only (pathway x phase 2x3; CNV proportion by sample)
#   Fig7: NEW independent MRD figure (3 panels)
suppressMessages({library(Seurat); library(ggplot2)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
qa   <- file.path(pkg, "Analysis_QA"); figd <- file.path(pkg, "Figures")

pub_theme <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "Arial") %+replace%
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "#444444", fill = NA, linewidth = 0.4),
          axis.title = element_text(face = "bold", size = base_size, colour = "black"),
          axis.text  = element_text(size = base_size - 0.5, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0.5),
          legend.title = element_text(face = "bold", size = base_size),
          legend.text  = element_text(size = base_size - 0.5),
          plot.tag    = element_text(face = "bold", size = base_size + 3),
          plot.tag.position = "topleft",
          plot.margin = margin(3, 3, 3, 3, "mm"))
}
c_blue <- "#2a78d6"; c_orange <- "#eb6834"; c_aqua <- "#1baf7a"; c_violet <- "#4a3aa7"; c_red <- "#e34948"
phase_cols <- c("G1" = c_blue, "S" = c_orange, "G2M" = c_aqua)
wrap_txt <- function(s, w = 38) vapply(s, function(x) paste(strwrap(x, w), collapse = "\n"), "")

compose <- function(plots, nrow, ncol, file, w = 6.5, h = 3, widths = NULL) {
  cairo_pdf(file, width = w, height = h, family = "Arial")
  grid::grid.newpage()
  if (is.null(widths)) widths <- rep(1/ncol, ncol)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol,
    widths = grid::unit(widths, "npc"), heights = grid::unit(rep(1/nrow, nrow), "npc"))))
  for (i in seq_along(plots)) {
    rr <- ceiling(i / ncol); cc <- (i - 1) %% ncol + 1
    print(plots[[i]], vp = grid::viewport(layout.pos.row = rr, layout.pos.col = cc))
  }
  grid::popViewport(); dev.off()
  png_file <- sub("\\.pdf$", ".png", file)
  system2("pdftoppm", c("-png", "-r", "300", "-singlefile", file, sub("\\.png$", "", png_file)))
  cat("wrote", file, "\n")
}

# ---------- data ----------
md <- read.csv(gzfile(file.path(qa, "08_cell_annotation_export.csv.gz")), stringsAsFactors = FALSE)
md$bin <- factor(cut(md$InferCNVcount, breaks = c(0, 0.25, 0.75, 1.5, Inf),
                     labels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"),
                     include.lowest = TRUE, right = FALSE),
                 levels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"))
md$Phase <- factor(md$Phase, levels = c("G1", "S", "G2M"))

# ================= Figure 3 (A + C) =================
phx <- as.data.frame(prop.table(table(md$bin, md$Phase), 1) * 100)
names(phx) <- c("bin", "phase", "pct")
sg2m <- aggregate(pct ~ bin, phx[phx$phase %in% c("S", "G2M"), ], sum)
p3A <- ggplot(phx, aes(bin, pct, fill = phase)) + geom_col(width = 0.62, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = phase_cols, name = "Phase") +
  geom_text(data = sg2m, aes(bin, 100 + 2.5, label = sprintf("%.1f%%", pct)),
            inherit.aes = FALSE, size = 2.6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 110), expand = expansion(mult = c(0, 0))) +
  labs(tag = "A", x = "InferCNV count bin", y = "Cell proportion (%)",
       title = "Cell-cycle composition per CNV bin (numbers = S+G2M)") +
  pub_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7))

src <- file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final")
gof <- c("0-0.25" = "GO_BP_Enrichment_0-0.25_Bin.csv", "0.25-0.75" = "GO_BP_Enrichment_0.25-0.75_Bin.csv",
         "0.75-1.5" = "GO_BP_Enrichment_0.75-1.5_Bin.csv")
go_df <- do.call(rbind, lapply(names(gof), function(b) {
  d <- read.csv(file.path(src, gof[[b]]), check.names = FALSE)
  d <- d[order(d$p.adjust), ][1:6, ]
  data.frame(bin = b, term = wrap_txt(sub(" signaling pathway", " signaling", sub(" regulation of", " reg. of", d$Description)), 34),
             nlogP = -log10(d$p.adjust), Count = d$Count)
}))
go_df$bin <- factor(go_df$bin, levels = names(gof))
go_df$term <- factor(go_df$term, levels = rev(unique(go_df$term)))
p3C <- ggplot(go_df, aes(nlogP, term, size = Count)) + geom_point(colour = c_blue, alpha = 0.85) +
  facet_wrap(~ bin, ncol = 3, scales = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(tag = "C", x = "-log10(adjusted P)", y = "GO biological process",
       title = "Enrichment per CNV bin", size = "Gene count") +
  pub_theme(base_size = 6.5) + theme(strip.text = element_text(face = "bold", size = 7),
                                     axis.text.y = element_text(size = 5.5))
compose(list(p3A, p3C), 1, 2, file.path(figd, "Figure3_CNV_CellCycle_Programs.pdf"), w = 9.5, h = 3.6,
        widths = c(1, 1.4))

# ================= Figure 4 =================
perf <- read.csv(file.path(src, "expression_matrices/模型详细性能表.csv"), check.names = FALSE, fileEncoding = "UTF-8-BOM")
models <- c("SVM", "LR", "RF")
perf_long <- data.frame(model = factor(rep(models, 4), levels = models),
                        metric = factor(rep(c("Accuracy", "F1 score", "AUC", "LOOPO accuracy"), each = 3),
                                        levels = c("Accuracy", "F1 score", "AUC", "LOOPO accuracy")),
                        value = c(0.9946, 0.9943, 0.9910, 0.9960, 0.9958, 0.9933, 0.9992, 0.9988, 0.9975,
                                  NA, 0.9863, 0.9916))
p4A <- ggplot(perf_long, aes(metric, value, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.4f", value)), position = position_dodge(width = 0.7),
            vjust = -0.4, size = 2.2) +
  scale_fill_manual(values = c("SVM" = c_blue, "LR" = c_orange, "RF" = c_aqua), name = "Model") +
  scale_y_continuous(limits = c(0.95, 1.01), breaks = c(0.95, 0.98, 1.0)) +
  labs(tag = "A", x = "Metric", y = "Score", title = "Classifier performance (5-fold CV)") +
  pub_theme(base_size = 7) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, size = 7))

pred <- read.csv(file.path(src, "expression_matrices/negative_data_predictions.csv"), stringsAsFactors = FALSE)
p4B <- ggplot(pred, aes(malignant_probability, fill = prediction_label)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6,
                 colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = c("ML_maglinant" = c_orange, "ML_normal" = c_aqua), name = "Predicted label") +
  annotate("text", x = 0.15, y = Inf, vjust = 1.6, size = 2.4, hjust = 0,
           label = sprintf("ML_maglinant: %d (75.7%%)\nML_normal: %d (24.3%%)",
                           sum(pred$prediction_label == "ML_maglinant"),
                           sum(pred$prediction_label == "ML_normal"))) +
  labs(tag = "B", x = "Predicted malignant probability", y = "Number of cells",
       title = "Test set (1,020 CNV-negative cells)") +
  pub_theme(base_size = 7) + theme(legend.position = "bottom")

lab <- read.csv(file.path(qa, "10_prediction_label_by_sample.csv"), row.names = 1, check.names = FALSE)
lab$sample <- as.character(rownames(lab))
lab <- lab[lab$sample != "NA" & !is.na(lab$sample), ]
lab$tot <- lab$ML_maglinant + lab$ML_normal
lab <- lab[lab$tot > 0, ][order(-lab$tot), ]
lab$sample <- factor(lab$sample, levels = lab$sample)
ll <- data.frame(sample = rep(lab$sample, 2),
                 label = factor(rep(c("ML_maglinant", "ML_normal"), each = nrow(lab)),
                                levels = c("ML_maglinant", "ML_normal")),
                 n = c(lab$ML_maglinant, lab$ML_normal))
p4C <- ggplot(ll, aes(sample, n, fill = label)) + geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("ML_maglinant" = c_orange, "ML_normal" = c_aqua), name = "Predicted label") +
  labs(tag = "C", x = "Sample", y = "Number of cells (CNV-negative)",
       title = "772 CNV-silent cells re-classified malignant") +
  pub_theme(base_size = 7) + theme(legend.position = "bottom")
compose(list(p4A, p4B, p4C), 1, 3, file.path(figd, "Figure4_ML_Classification.pdf"), w = 10.5, h = 3.6)

# ================= Figure 6 (longitudinal only) =================
pw_orig <- read.csv(file.path(qa, "20_pathway_strength_malignant.csv"), stringsAsFactors = FALSE)
ph_o <- c("GSM3872440_PRE-T_1" = "PreT", "GSM3872441_PRE-T_2" = "PreT", "GSM7091996_T1Dx" = "Dx",
          "GSM7091997_T2Dx" = "Dx", "GSM7091998_T3Dx" = "Dx", "GSM7091999_T4Dx" = "Dx", "GSM7092000_T5Dx" = "Dx",
          "GSM7092001_T6Dx" = "Dx", "GSM7092003_T7Dx" = "Dx", "GSM7092005_T8Dx" = "Dx", "GSM7092007_T9Dx" = "Dx",
          "GSM7092009_T10Dx" = "Dx", "GSM7092002_T6EOI" = "EOI", "GSM7092004_T7EOI" = "EOI", "GSM7092006_T8EOI" = "EOI",
          "GSM7092008_T9EOI" = "EOI", "GSM7092010_T10EOI" = "EOI", "GSM7092011_T11Rel" = "Rel")
pw_orig$phase <- factor(unname(ph_o[pw_orig$sample]), levels = c("PreT", "Dx", "EOI", "Rel"))
top6 <- names(sort(tapply(pw_orig$prob, pw_orig$pathway, sum), decreasing = TRUE))[1:6]
pw2 <- pw_orig[pw_orig$pathway %in% top6, ]; pw2$pathway <- factor(pw2$pathway, levels = top6)
summ <- aggregate(prob ~ pathway + phase, pw2, function(x) c(m = mean(x), se = sd(x) / sqrt(length(x))))
summ <- do.call(data.frame, summ); names(summ)[3:4] <- c("m", "se")
p6A <- ggplot(summ, aes(phase, m)) + geom_col(fill = c_blue, width = 0.6, colour = "white", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(0, m - se), ymax = m + se), width = 0.2, linewidth = 0.3) +
  facet_wrap(~ pathway, scales = "free_y", nrow = 2) +
  labs(tag = "A", x = "Disease phase", y = "Mean pathway strength",
       title = "Communication by phase (subcluster design, mean ± SE)") +
  pub_theme(base_size = 7) + theme(axis.text.x = element_text(size = 6.5))
prop_s <- as.data.frame(prop.table(table(md$sample_id, md$InferCNV), 1) * 100)
names(prop_s) <- c("sample", "status", "pct")
sord <- c(sort(unique(md$sample_id[grepl("Dx", md$sample_id) & !grepl("EOI|Rel", md$sample_id)])),
          sort(unique(md$sample_id[grepl("EOI", md$sample_id)])),
          unique(md$sample_id[grepl("Rel", md$sample_id)]), "PRE-T")
prop_s$sample <- factor(prop_s$sample, levels = sord)
p6B <- ggplot(prop_s, aes(sample, pct, fill = status)) + geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Positive" = c_red, "Negative" = c_blue), name = "inferCNV status") +
  labs(tag = "B", x = "Sample (Dx → EOI → Rel)", y = "Proportion of cells (%)",
       title = "CNV-positive proportion per sample") +
  pub_theme(base_size = 7) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
                                   legend.position = "bottom")
compose(list(p6A, p6B), 1, 2, file.path(figd, "Figure6_Longitudinal_Dynamics.pdf"), w = 10.5, h = 3.4)

# ================= Figure 7 (independent MRD figure) =================
per_s <- read.csv(file.path(qa, "31_MRD_EOI_per_sample.csv"), stringsAsFactors = FALSE)
per_s$MRD_grp <- factor(ifelse(grepl("MRD\\+", per_s$MRD), "MRD+", "MRD-"), levels = c("MRD+", "MRD-"))
per_s$GSM <- factor(per_s$GSM, levels = per_s$GSM[order(per_s$GSM)])
per_s$lab <- sprintf("%s\nMRD %s", gsub("GSM70920(0[2-9]|10)", "T\\1-", per_s$GSM), sub(" \\(.*", "", per_s$MRD))
ps_long <- data.frame(GSM = rep(per_s$GSM, 2), MRD = rep(per_s$MRD_grp, 2),
                      lab = rep(per_s$lab, 2),
                      metric = factor(rep(c("pct_CNV_Pos", "pct_SG2M"), each = nrow(per_s)),
                                      levels = c("pct_CNV_Pos", "pct_SG2M"),
                                      labels = c("CNV-positive cells (%)", "S+G2M cells (%)")),
                      value = c(per_s$pct_CNV_Pos, per_s$pct_SG2M))
p7A <- ggplot(ps_long, aes(GSM, value, fill = MRD)) + geom_col(width = 0.6, colour = "white", linewidth = 0.4) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("MRD+" = c_red, "MRD-" = c_blue), name = "EOI MRD status") +
  labs(tag = "A", x = "EOI sample", y = "Value (%)",
       title = "CNV burden and proliferation of EOI residual cells") +
  pub_theme(base_size = 7) + theme(axis.text.x = element_text(size = 6.5), legend.position = "bottom")

gexp <- read.csv(file.path(qa, "32_MRD_EOI_gene_expression.csv"), stringsAsFactors = FALSE)
gexp$MRD_grp <- factor(ifelse(grepl("MRD\\+", gexp$MRD), "MRD+", "MRD-"), levels = c("MRD+", "MRD-"))
g4 <- c("mean_PDCD1", "mean_HAVCR2", "mean_LGALS9", "mean_CD74", "mean_MIF", "mean_SOX4")
g4n <- c("PD-1", "TIM-3", "LGALS9", "CD74", "MIF", "SOX4")
gg <- do.call(rbind, lapply(seq_along(g4), function(i) {
  data.frame(MRD = gexp$MRD_grp, gene = factor(g4n[i], levels = g4n), value = gexp[[g4[i]]])
}))
p7B <- ggplot(gg, aes(gene, value, fill = MRD)) + geom_boxplot(width = 0.55, outlier.size = 0.5, linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.75), size = 1, colour = "black") +
  scale_fill_manual(values = c("MRD+" = c_red, "MRD-" = c_blue), name = "EOI MRD status") +
  labs(tag = "B", x = "Gene", y = "Mean expression (CP10K)",
       title = "MRD-associated communication molecules (EOI)") +
  pub_theme(base_size = 7) + theme(legend.position = "bottom")

eo_samples <- c("GSM7092002_T6EOI", "GSM7092004_T7EOI", "GSM7092006_T8EOI",
               "GSM7092008_T9EOI", "GSM7092010_T10EOI")
mrd_pos <- c("GSM7092002_T6EOI", "GSM7092006_T8EOI", "GSM7092008_T9EOI")
comm <- pw_orig[pw_orig$sample %in% eo_samples, ]
comm$MRD_grp <- factor(ifelse(comm$sample %in% mrd_pos, "MRD+", "MRD-"), levels = c("MRD+", "MRD-"))
comm$sample2 <- factor(sub("GSM70920(0[2-9]|10)_", "", comm$sample),
                       levels = c("T6EOI", "T7EOI", "T8EOI", "T9EOI", "T10EOI"))
cs_tot <- aggregate(prob ~ sample2 + MRD_grp, comm, sum); names(cs_tot)[3] <- "total"
cm_path <- comm[comm$pathway %in% c("MIF", "GALECTIN"), ]
cm2 <- aggregate(prob ~ sample2 + MRD_grp + pathway, cm_path, sum)
p7C <- ggplot(cm2, aes(sample2, prob, fill = pathway)) + geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("MIF" = c_blue, "GALECTIN" = c_orange), name = "Pathway") +
  labs(tag = "C", x = "EOI sample", y = "Pathway strength",
       title = "MIF/GALECTIN signaling at EOI by MRD") +
  pub_theme(base_size = 7) + theme(axis.text.x = element_text(size = 6.5), legend.position = "bottom") +
  annotate("text", x = 5.4, y = max(cm2$prob) * 1.02, vjust = 1, hjust = 1, size = 2.2,
           label = paste("MRD+:", paste(sort(unique(cs_tot$sample2[cs_tot$MRD_grp == "MRD+"])), collapse = ", "),
                         "  MRD-:", paste(sort(unique(cs_tot$sample2[cs_tot$MRD_grp == "MRD-"])), collapse = ", ")))
compose(list(p7A, p7B, p7C), 1, 3, file.path(figd, "Figure7_MRD_Reanalysis.pdf"), w = 10.5, h = 3.6)
cat("ALL REV2 FIGURES DONE\n")
