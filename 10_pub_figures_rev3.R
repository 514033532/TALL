#!/usr/bin/env Rscript
# 10: Figure Rev3 — fixes every item of the user's per-figure SCI-format critique.
#   - ggarrange-based composition (robust rendering, bold top-left A/B/C tags)
#   - rotated+shortened x labels; value labels with gaps; legend placement fixes
#   - UMAP shuffle+alpha; histogram re-binned; heatmap white borders
#   - Fig3 = A+C (B dropped); Fig4 LOOPO panel; Fig5BC full render; Fig7 fixes
suppressMessages({library(Seurat); library(ggplot2); library(ggpubr); library(pheatmap)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
qa   <- file.path(pkg, "Analysis_QA"); figd <- file.path(pkg, "Figures")

pub_theme <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "Arial") %+replace%
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "#444444", fill = NA, linewidth = 0.4),
          axis.title = element_text(face = "bold", size = base_size, colour = "black",
                                    margin = margin(t = 4, r = 4)),
          axis.title.x = element_text(face = "bold", size = base_size, colour = "black",
                                      angle = 0, margin = margin(t = 6)),
          axis.title.y = element_text(face = "bold", size = base_size, colour = "black",
                                      angle = 90, margin = margin(r = 6)),
          axis.text  = element_text(size = base_size - 0.5, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0.5,
                                    margin = margin(b = 8)),
          legend.title = element_text(face = "bold", size = base_size),
          legend.text  = element_text(size = base_size - 0.5),
          legend.margin = margin(t = 2, b = 2),
          plot.margin = margin(6, 6, 6, 6, "mm"))
}
c_blue <- "#2a78d6"; c_orange <- "#eb6834"; c_aqua <- "#1baf7a"; c_violet <- "#4a3aa7"; c_red <- "#e34948"
phase_cols <- c("G1" = c_blue, "S" = c_orange, "G2M" = c_aqua)
short_nm <- function(x) gsub("_", " ", sub("PRE-T", "PRE-T", x))

comp <- function(plots, labels, file, w, h, ncol = length(plots), nrow = 1,
                 widths = NULL, common.legend = FALSE, legend = NULL) {
  args <- list(plotlist = plots, labels = labels, ncol = ncol, nrow = nrow,
               common.legend = common.legend, legend = legend,
               font.label = list(size = 11, face = "bold", family = "Arial"),
               label.x = 0.02, label.y = 0.97, hjust = 0, vjust = 1, align = "hv")
  if (!is.null(widths)) args$widths <- widths
  arr <- tryCatch(do.call(ggarrange, args),
    error = function(e) { cat("COMP FAIL:", file, "->", conditionMessage(e), "\n"); NULL })
  if (is.null(arr)) return(invisible(NULL))
  ggsave(file, arr, device = cairo_pdf, width = w, height = h)
  png_file <- sub("\\.pdf$", ".png", file)
  system2("pdftoppm", c("-png", "-r", "300", "-singlefile", file, sub("\\.png$", "", png_file)))
  cat("wrote", file, "\n")
}
wrap_txt <- function(s, w = 30) vapply(s, function(x) paste(strwrap(x, w), collapse = "\n"), "")

# ---------- data ----------
md <- read.csv(gzfile(file.path(qa, "08_cell_annotation_export.csv.gz")), stringsAsFactors = FALSE)
md$bin <- factor(cut(md$InferCNVcount, breaks = c(0, 0.25, 0.75, 1.5, Inf),
                     labels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"),
                     include.lowest = TRUE, right = FALSE),
                 levels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"))
md$Phase <- factor(md$Phase, levels = c("G1", "S", "G2M"))
md$sample_id <- factor(md$sample_id, levels = sort(unique(md$sample_id)))
obj <- readRDS(file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/marker-v2/combined_TALL_22679_prediction_label_seurat.rds"))
obj$sample_id <- factor(obj$sample_id, levels = levels(md$sample_id))

# ================= Figure 1 =================
sm_cnt <- as.data.frame(table(md$sample_id)); names(sm_cnt) <- c("sample", "n")
p1A <- ggplot(sm_cnt, aes(sample, n)) + geom_col(fill = c_blue, width = 0.6, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n), vjust = -0.5, size = 2) +
  scale_x_discrete(labels = short_nm) +
  labs(x = "Sample", y = "Number of cells", title = "Cells per sample (n = 22,679)") +
  pub_theme(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
                       plot.margin = margin(4, 4, 8, 4, "mm"))
st <- ifelse(grepl("PRE-T", md$sample_id), "Pre-T", ifelse(grepl("EOI", md$sample_id), "EOI",
        ifelse(grepl("Rel", md$sample_id), "Relapse", "Diagnosis")))
st <- factor(st, levels = c("Diagnosis", "EOI", "Relapse", "Pre-T"))
p1B <- ggplot(data.frame(st), aes(st)) + geom_bar(fill = c_orange, width = 0.6, colour = "white", linewidth = 0.4) +
  labs(x = "Sample type", y = "Number of cells", title = "Cells per sample type") +
  pub_theme(7) + theme(axis.text.x = element_text(angle = 0, size = 7))
p1C <- DimPlot(obj, reduction = "umap", group.by = "sample_id", pt.size = 0.25,
               shuffle = TRUE, alpha = 0.6, label = FALSE) +
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(12),
                     labels = short_nm, name = "Sample") +
  labs(x = "UMAP_1", y = "UMAP_2", title = "UMAP by sample") +
  pub_theme(7) + theme(legend.text = element_text(size = 5.5), legend.key.size = unit(0.25, "cm"),
                       legend.position = "right")
comp(list(p1A, p1B, p1C), c("A", "B", "C"), file.path(figd, "Figure1_Cohort_Overview.pdf"), w = 10.5, h = 3.4)

# ================= Figure 2 =================
p2A <- ggplot(md, aes(InferCNV)) + geom_bar(aes(fill = InferCNV), width = 0.6, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = c_red, "Negative" = c_blue), name = "inferCNV status") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 2.2) +
  labs(x = "inferCNV status", y = "Number of cells", title = "CNV status") +
  pub_theme(7) + theme(legend.position = "none")
p2B <- ggplot(md, aes(InferCNVcount)) + geom_histogram(bins = 30, fill = c_blue, colour = "white", linewidth = 0.2) +
  labs(x = "InferCNV count (CNVs per cell)", y = "Number of cells", title = "CNV burden distribution") +
  pub_theme(7)
p2C <- DimPlot(obj, reduction = "umap", group.by = "InferCNV", pt.size = 0.3,
               shuffle = TRUE, alpha = 0.5) +
  scale_color_manual(values = c("Positive" = c_red, "Negative" = c_blue), name = "inferCNV status") +
  labs(x = "UMAP_1", y = "UMAP_2", title = "UMAP by CNV status") +
  pub_theme(7) + theme(legend.position = "none")
p2D <- FeaturePlot(obj, features = "InferCNVcount", reduction = "umap", pt.size = 0.3,
                   order = FALSE) +
  scale_color_gradientn(colours = c("#cde2fb", "#5598e7", "#104281")) +
  labs(x = "UMAP_1", y = "UMAP_2", title = "UMAP by CNV burden", colour = "InferCNV\ncount") +
  pub_theme(7) + theme(legend.position = "right", legend.key.height = unit(0.3, "cm"))
prop_s <- as.data.frame(prop.table(table(md$sample_id, md$InferCNV), 1) * 100)
names(prop_s) <- c("sample", "status", "pct")
p2E <- ggplot(prop_s, aes(sample, pct, fill = status)) + geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Positive" = c_red, "Negative" = c_blue), name = "inferCNV status") +
  scale_x_discrete(labels = short_nm) +
  labs(x = "Sample", y = "Proportion of cells (%)", title = "CNV-positive proportion") +
  pub_theme(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
                       legend.position = "bottom", plot.margin = margin(4, 4, 8, 4, "mm"))
comp(list(p2A, p2B, p2C, p2D, p2E), c("A", "B", "C", "D", "E"),
     file.path(figd, "Figure2_CNV_Landscape.pdf"), w = 11, h = 7, ncol = 3, nrow = 2)

# ================= Figure 3 (A + C) =================
phx <- as.data.frame(prop.table(table(md$bin, md$Phase), 1) * 100)
names(phx) <- c("bin", "phase", "pct")
sg2m <- aggregate(pct ~ bin, phx[phx$phase %in% c("S", "G2M"), ], sum)
p3A <- ggplot(phx, aes(bin, pct, fill = phase)) + geom_col(width = 0.62, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = phase_cols, name = "Phase") +
  geom_text(data = sg2m, aes(bin, 100 + 3, label = sprintf("%.1f%%", pct)),
            inherit.aes = FALSE, size = 2.6, fontface = "bold") +
  scale_y_continuous(limits = c(0, 112), expand = expansion(mult = c(0, 0))) +
  labs(x = "InferCNV count bin", y = "Cell proportion (%)",
       title = "Cell-cycle composition (numbers = S+G2M)") +
  pub_theme(7.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
                         plot.margin = margin(4, 4, 8, 4, "mm"))
src <- file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final")
gof <- c("0-0.25" = "GO_BP_Enrichment_0-0.25_Bin.csv", "0.25-0.75" = "GO_BP_Enrichment_0.25-0.75_Bin.csv",
         "0.75-1.5" = "GO_BP_Enrichment_0.75-1.5_Bin.csv")
go_df <- do.call(rbind, lapply(names(gof), function(b) {
  d <- read.csv(file.path(src, gof[[b]]), check.names = FALSE)
  d <- d[order(d$p.adjust), ][1:6, ]
  data.frame(bin = b, term = wrap_txt(sub(" signaling pathway", " signaling",
              sub(" regulation of", " reg. of", d$Description)), 28),
             nlogP = -log10(d$p.adjust), Count = d$Count)
}))
go_df$bin <- factor(go_df$bin, levels = names(gof))
go_df$term <- factor(go_df$term, levels = rev(unique(go_df$term)))
p3C <- ggplot(go_df, aes(nlogP, term, size = Count)) + geom_point(colour = c_blue, alpha = 0.85) +
  facet_wrap(~ bin, ncol = 3, scales = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.1))) +
  labs(x = "-log10(adjusted P)", y = "GO biological process",
       title = "Enrichment per CNV bin", size = "Gene\ncount") +
  pub_theme(7) + theme(strip.text = element_text(face = "bold", size = 7.5),
                         axis.text.y = element_text(size = 6))
comp(list(p3A, p3C), c("A", "B"), file.path(figd, "Figure3_CNV_CellCycle_Programs.pdf"),
     w = 11.5, h = 4.6, ncol = 2, widths = c(1, 2))

# ================= Figure 4 =================
# unified model colors: SVM blue, LR violet, RF orange
mcols <- c("SVM" = c_blue, "LR" = c_violet, "RF" = c_orange)
models <- c("SVM", "LR", "RF")
perf_long <- data.frame(model = factor(rep(models, 4), levels = models),
                        metric = factor(rep(c("Accuracy", "F1", "AUC", "LOOPO"), each = 3),
                                        levels = c("Accuracy", "F1", "AUC", "LOOPO")),
                        value = c(0.9946, 0.9943, 0.9910, 0.9960, 0.9958, 0.9933, 0.9992, 0.9988, 0.9975,
                                  NA, 0.9863, 0.9916))
perf_long <- perf_long[!is.na(perf_long$value), ]
perf_long$bar_y <- perf_long$value - 0.94   # baseline-shifted so bars are fully visible in zoom
p4A <- ggplot(perf_long, aes(metric, bar_y, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.4f", value)), position = position_dodge(width = 0.75),
            vjust = -0.5, size = 1.7, angle = 0, colour = "black") +
  scale_fill_manual(values = mcols, name = "Model") +
  scale_y_continuous(limits = c(0, 0.078), breaks = c(0, 0.04, 0.06),
                     labels = c("0.94", "0.98", "1.00")) +
  labs(x = "Metric", y = "Score (axis starts at 0.94)",
       title = "Classifier performance: 5-fold CV and leave-one-patient-out (LOOPO)") +
  pub_theme(7.5) + theme(legend.position = "bottom", axis.text.x = element_text(size = 7))
pred <- read.csv(file.path(src, "expression_matrices/negative_data_predictions.csv"), stringsAsFactors = FALSE)
p4B <- ggplot(pred, aes(malignant_probability, fill = prediction_label)) +
  geom_histogram(binwidth = 0.05, position = "identity", alpha = 0.6, colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3) +
  scale_y_continuous(trans = scales::log1p_trans()) +
  scale_fill_manual(values = c("ML_maglinant" = c_blue, "ML_normal" = c_red),
                    labels = c("Malignant", "Normal"), name = "Predicted label") +
  labs(x = "Predicted malignant probability", y = "log1p(number of cells)",
       title = "Test set (n = 1,020): 772 malignant / 248 normal\ndashed line = decision threshold (0.5)") +
  pub_theme(7.5) + theme(legend.position = "bottom")
lab <- read.csv(file.path(qa, "10_prediction_label_by_sample.csv"), row.names = 1, check.names = FALSE)
lab$sample <- as.character(rownames(lab))
lab <- lab[!is.na(lab$sample) & lab$sample != "" & lab$sample != "NA", ]
lab$tot <- lab$ML_maglinant + lab$ML_normal
lab <- lab[!is.na(lab$tot) & lab$tot > 0, ]
time_order <- c("PRE-T", "T1Dx", "T2Dx", "T3Dx", "T4Dx", "T6Dx_T6EOI", "T7Dx_T7EOI",
                "T8Dx_T8EOI", "T9Dx_T9EOI", "T10Dx_T10EOI", "T11Rel")
lab$sample <- factor(lab$sample, levels = intersect(time_order, lab$sample))
ll <- data.frame(sample = rep(lab$sample, 2),
                 label = factor(rep(c("ML_maglinant", "ML_normal"), each = nrow(lab)),
                                levels = c("ML_maglinant", "ML_normal")),
                 n = c(lab$ML_maglinant, lab$ML_normal))
p4C <- ggplot(ll, aes(sample, n, fill = label)) + geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  coord_flip() +
  scale_x_discrete(labels = short_nm) +
  scale_fill_manual(values = c("ML_maglinant" = c_blue, "ML_normal" = c_red),
                    labels = c("Malignant", "Normal"), name = "Predicted label") +
  geom_text(aes(label = ifelse(n > 0, n, "")), position = position_stack(vjust = 0.5),
            size = 2, colour = "white", fontface = "bold") +
  labs(x = "Sample (time order)", y = "Number of cells (CNV-negative)",
       title = "Re-classified CNV-silent cells") +
  pub_theme(7.5) + theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5))
comp(list(p4A, p4B, p4C), c("A", "B", "C"), file.path(figd, "Figure4_ML_Classification.pdf"),
     w = 11.2, h = 4.6, ncol = 3, widths = c(1.05, 0.9, 0.95))

# ================= Figure 5A (heatmap) =================
pw <- read.csv(file.path(qa, "34_harmonized_pathway_strength.csv"), stringsAsFactors = FALSE)
phase_map <- c("GSM3872440" = "PreT", "GSM3872441" = "PreT", "GSM7091996" = "Dx", "GSM7091997" = "Dx",
               "GSM7091998" = "Dx", "GSM7091999" = "Dx", "GSM7092000" = "Dx", "GSM7092001" = "Dx",
               "GSM7092003" = "Dx", "GSM7092005" = "Dx", "GSM7092007" = "Dx", "GSM7092009" = "Dx",
               "GSM7092002" = "EOI", "GSM7092004" = "EOI", "GSM7092006" = "EOI", "GSM7092008" = "EOI",
               "GSM7092010" = "EOI", "GSM7092011" = "Rel")
samp_ord <- c("GSM3872440","GSM3872441","GSM7091996","GSM7091997","GSM7091998","GSM7091999","GSM7092000",
              "GSM7092001","GSM7092003","GSM7092005","GSM7092007","GSM7092009","GSM7092002","GSM7092004",
              "GSM7092006","GSM7092008","GSM7092010","GSM7092011")
short5 <- c("PRE-T_1","PRE-T_2","T1Dx","T2Dx","T3Dx","T4Dx","T5Dx","T6Dx","T7Dx","T8Dx","T9Dx","T10Dx",
            "T6EOI","T7EOI","T8EOI","T9EOI","T10EOI","T11Rel"); names(short5) <- samp_ord
avail <- intersect(samp_ord, unique(pw$sample))
wide <- reshape(pw[pw$sample %in% avail, c("sample", "pathway_name", "prob")],
                idvar = "pathway_name", timevar = "sample", direction = "wide")
names(wide) <- sub("prob\\.", "", names(wide)); rownames(wide) <- wide$pathway_name; wide$pathway_name <- NULL
wide <- wide[, avail]; wide[is.na(wide)] <- 0
wide <- wide[rowSums(wide > 0) >= 1, , drop = FALSE]
logm <- log10(as.matrix(wide) + 1e-3)
phases <- data.frame(Phase = unname(phase_map[colnames(logm)]), row.names = colnames(logm))
gaps <- which(!duplicated(phases$Phase[order(match(phases$Phase, c("PreT", "Dx", "EOI", "Rel")))])) - 1
gaps <- gaps[gaps > 0]
ht <- pheatmap(logm, cluster_rows = FALSE, cluster_cols = FALSE,
               color = c("#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7","#3987e5","#2a78d6",
                         "#256abf","#1c5cab","#184f95","#104281"),
               annotation_col = phases, labels_col = short5[colnames(logm)],
               annotation_colors = list(Phase = c("PreT" = c_violet, "Dx" = c_blue, "EOI" = c_orange, "Rel" = c_red)),
               fontsize = 8, fontsize_col = 7, angle_col = 45, fontfamily = "Arial",
               breaks = seq(-3, 2, length.out = 13), border_color = "white", cellwidth = 13, cellheight = 12,
               legend = TRUE, annotation_legend = TRUE,
               gaps_col = if (length(gaps) > 0) gaps else NULL, silent = TRUE)
cairo_pdf(file.path(figd, "Figure5A_harmonized_heatmap.pdf"), width = 7, height = 4.8, family = "Arial")
grid::grid.newpage()
grid::grid.draw(ht$gtable)
grid::grid.text("A", x = unit(0.012, "npc"), y = unit(0.985, "npc"), just = c("left", "top"),
                gp = grid::gpar(fontsize = 11, fontface = "bold", fontfamily = "Arial"))
dev.off()
system2("pdftoppm", c("-png", "-r", "300", "-singlefile",
                      file.path(figd, "Figure5A_harmonized_heatmap.pdf"), file.path(figd, "Figure5A_harmonized_heatmap")))
cat("wrote Figure5A\n")

# ================= Figure 5BC =================
gp <- read.csv(file.path(qa, "36_harmonized_sender_pathway.csv"), stringsAsFactors = FALSE)
snd <- aggregate(prob ~ source, gp, sum); snd <- snd[order(snd$prob), ]
snd$source <- factor(snd$source, levels = snd$source)
p5B <- ggplot(snd, aes(source, prob)) + geom_col(fill = c_blue, width = 0.6, colour = "white", linewidth = 0.4) +
  coord_flip() +
  geom_text(aes(label = sprintf("%.1f", prob)), hjust = -0.15, size = 2.4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Sender (CNV burden bin)", y = "Total outgoing strength",
       title = "Signal output by CNV bin") + pub_theme(7)
top_pw <- names(sort(tapply(gp$prob, gp$pathway_name, sum), decreasing = TRUE))[1:6]
gpt <- gp[gp$pathway_name %in% top_pw, ]
gpt$pathway_name <- factor(gpt$pathway_name, levels = top_pw)
p5C <- ggplot(gpt, aes(source, prob, fill = pathway_name)) +
  geom_col(position = "stack", width = 0.6, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(c_blue, c_orange, c_aqua, c_violet, c_red, "#e87ba4"), name = "Pathway") +
  labs(x = "Sender (CNV burden bin)", y = "Pathway strength",
       title = "Top pathways by sender bin") +
  pub_theme(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.5),
                       legend.position = "bottom")
comp(list(p5B, p5C), c("B", "C"), file.path(figd, "Figure5BC_sender_contribution.pdf"), w = 8, h = 3.4)

# ================= Figure 6 =================
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
p6A <- ggplot(summ, aes(phase, m)) + geom_col(fill = c_blue, width = 0.62, colour = "white", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(0, m - se), ymax = m + se), width = 0.22, linewidth = 0.35) +
  facet_wrap(~ pathway, scales = "free_y", nrow = 2) +
  labs(x = "Disease phase", y = "Mean pathway strength",
       title = "Communication by phase (mean ± SE; y-scales differ per panel)") +
  pub_theme(7) + theme(axis.text.x = element_text(size = 6.5))
prop_s <- as.data.frame(prop.table(table(md$sample_id, md$InferCNV), 1) * 100)
names(prop_s) <- c("sample", "status", "pct")
sord <- c(sort(unique(as.character(md$sample_id[grepl("Dx", md$sample_id) & !grepl("EOI|Rel", md$sample_id)]))),
          sort(unique(as.character(md$sample_id[grepl("EOI", md$sample_id)]))),
          "T11Rel", "PRE-T")
prop_s$sample <- factor(prop_s$sample, levels = sord)
p6B <- ggplot(prop_s, aes(sample, pct, fill = status)) + geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Positive" = c_red, "Negative" = c_blue), name = "inferCNV status") +
  scale_x_discrete(labels = short_nm) +
  labs(x = "Sample (Dx → EOI → Rel)", y = "Proportion of cells (%)",
       title = "CNV-positive proportion per sample") +
  pub_theme(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
                       legend.position = "bottom", plot.margin = margin(4, 4, 8, 4, "mm"))
comp(list(p6A, p6B), c("A", "B"), file.path(figd, "Figure6_Longitudinal_Dynamics.pdf"), w = 10.5, h = 3.8,
     ncol = 2, widths = c(1.25, 1))

# ================= Figure 7 =================
per_s <- read.csv(file.path(qa, "31_MRD_EOI_per_sample.csv"), stringsAsFactors = FALSE)
per_s$MRD_grp <- factor(ifelse(grepl("MRD\\+", per_s$MRD), "MRD+", "MRD-"), levels = c("MRD+", "MRD-"))
per_s$sname <- factor(gsub("GSM70920(0[2-9]|10)", "T", per_s$GSM), levels = c("T6", "T7", "T8", "T9", "T10"))
per_s$mrdlab <- sprintf("MRD %s%%", sub("^MRD. \\(([0-9.]+)%\\)$", "\\1", per_s$MRD))
ps_long <- data.frame(sname = rep(per_s$sname, 2), MRD = rep(per_s$MRD_grp, 2),
                      mrdlab = rep(per_s$mrdlab, 2),
                      metric = factor(rep(c("pct_CNV_Pos", "pct_SG2M"), each = nrow(per_s)),
                                      levels = c("pct_CNV_Pos", "pct_SG2M"),
                                      labels = c("CNV-positive cells (%)", "S+G2M cells (%)")),
                      value = c(per_s$pct_CNV_Pos, per_s$pct_SG2M))
p7A <- ggplot(ps_long, aes(sname, value, fill = MRD)) + geom_col(width = 0.6, colour = "white", linewidth = 0.4) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  geom_text(aes(label = sprintf("%.1f", value)), vjust = -0.5, size = 2) +
  scale_fill_manual(values = c("MRD+" = c_red, "MRD-" = c_blue), name = "EOI MRD status") +
  labs(x = "EOI sample", y = "Value (%)",
       title = "Residual-cell CNV and proliferation (labels = EOI MRD %)") +
  pub_theme(7) + theme(legend.position = "none", axis.text.x = element_text(size = 7))
gexp <- read.csv(file.path(qa, "32_MRD_EOI_gene_expression.csv"), stringsAsFactors = FALSE)
gexp$MRD_grp <- factor(ifelse(grepl("MRD\\+", gexp$MRD), "MRD+", "MRD-"), levels = c("MRD+", "MRD-"))
g4 <- c("mean_PDCD1", "mean_HAVCR2", "mean_LGALS9", "mean_CD74", "mean_MIF", "mean_SOX4")
g4n <- c("PD-1", "TIM-3", "LGALS9", "CD74", "MIF", "SOX4")
gg <- do.call(rbind, lapply(seq_along(g4), function(i) {
  data.frame(MRD = gexp$MRD_grp, gene = factor(g4n[i], levels = g4n), value = gexp[[g4[i]]])
}))
p7B <- ggplot(gg, aes(gene, value, fill = MRD)) + geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
              size = 0.7, colour = "black", alpha = 0.7) +
  scale_fill_manual(values = c("MRD+" = c_red, "MRD-" = c_blue), name = "EOI MRD status") +
  labs(x = "Gene", y = "Mean expression (CP10K)",
       title = "MRD-associated communication molecules") +
  pub_theme(7) + theme(legend.position = "bottom")
eo_samples <- c("GSM7092002_T6EOI", "GSM7092004_T7EOI", "GSM7092006_T8EOI",
                "GSM7092008_T9EOI", "GSM7092010_T10EOI")
mrd_pos <- c("GSM7092002_T6EOI", "GSM7092006_T8EOI", "GSM7092008_T9EOI")
comm <- pw_orig[pw_orig$sample %in% eo_samples, ]
comm$MRD_grp <- factor(ifelse(comm$sample %in% mrd_pos, "MRD+", "MRD-"), levels = c("MRD+", "MRD-"))
comm$sname <- factor(sub("GSM70920(0[2-9]|10)_", "", comm$sample),
                     levels = c("T6EOI", "T7EOI", "T8EOI", "T9EOI", "T10EOI"))
cm_path <- comm[comm$pathway %in% c("MIF", "GALECTIN"), ]
cm2 <- aggregate(prob ~ sname + pathway, cm_path, sum)
p7C <- ggplot(cm2, aes(sname, prob, fill = pathway)) + geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", prob)), position = position_stack(vjust = 0.5), size = 2) +
  scale_fill_manual(values = c("MIF" = c_blue, "GALECTIN" = c_orange), name = "Pathway") +
  labs(x = "EOI sample", y = "Pathway strength", title = "MIF/GALECTIN signaling at EOI") +
  pub_theme(7) + theme(legend.position = "bottom", axis.text.x = element_text(size = 7))
comp(list(p7A, p7B, p7C), c("A", "B", "C"), file.path(figd, "Figure7_MRD_Reanalysis.pdf"), w = 10.5, h = 3.8)
cat("ALL REV3 FIGURES DONE\n")
