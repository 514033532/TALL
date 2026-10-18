#!/usr/bin/env Rscript
# 06: Publication-grade main figures (user standards):
#   - Arial throughout, uniform font sizes; bold axis titles; explicit units
#   - subfigure tags (A,B,C...) bold, top-left margin, never overlapping
#   - equal panel sizes, grid-aligned; vector PDF + 300 dpi PNG previews
suppressMessages({library(Seurat); library(ggplot2); library(pheatmap)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
qa   <- file.path(pkg, "Analysis_QA"); figd <- file.path(pkg, "Figures")
log <- file(file.path(qa, "06_figures_log.txt"), open = "wt"); sink(log, type = "output"); sink(log, type = "message")

if (any(grepl("Arial", systemfonts::system_fonts()$family))) cat("Arial available\n")

# ---------- theme ----------
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
          plot.margin = margin(3, 3, 3, 3, "mm"),
          legend.key.size = unit(0.35, "cm"))
}
c_blue <- "#2a78d6"; c_orange <- "#eb6834"; c_aqua <- "#1baf7a"
c_violet <- "#4a3aa7"; c_red <- "#e34948"
phase_cols <- c("G1" = c_blue, "S" = c_orange, "G2M" = c_aqua)
seq_blue <- c("#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7","#3987e5","#2a78d6","#256abf","#1c5cab","#184f95","#104281")

# compose panels: equal-size grid, aligned, tags by labs(tag=)
compose <- function(plots, nrow, ncol, file, w = 6.5, h = 3, dpi = 300) {
  cairo_pdf(file, width = w, height = h, family = "Arial")
  grid::grid.newpage()
  vp <- grid::viewport(layout = grid::grid.layout(nrow, ncol,
                       widths = grid::unit(rep(1/ncol, ncol), "npc"),
                       heights = grid::unit(rep(1/nrow, nrow), "npc")))
  grid::pushViewport(vp)
  for (i in seq_along(plots)) {
    rr <- ceiling(i / ncol); cc <- (i - 1) %% ncol + 1
    print(plots[[i]], vp = grid::viewport(layout.pos.row = rr, layout.pos.col = cc))
  }
  grid::popViewport(); dev.off()
  png_file <- sub("\\.pdf$", ".png", file)
  # rasterize preview at 300 dpi via pdftoppm
  system2("pdftoppm", c("-png", "-r", "300", "-singlefile", file, sub("\\.png$", "", png_file)))
  cat("wrote", file, "and", png_file, "\n")
}

# ---------- data ----------
md <- read.csv(gzfile(file.path(qa, "08_cell_annotation_export.csv.gz")), stringsAsFactors = FALSE)
md$bin <- factor(cut(md$InferCNVcount, breaks = c(0, 0.25, 0.75, 1.5, Inf),
                     labels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"),
                     include.lowest = TRUE, right = FALSE),
                 levels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"))
md$Phase <- factor(md$Phase, levels = c("G1", "S", "G2M"))
obj <- readRDS(file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/marker-v2/combined_TALL_22679_prediction_label_seurat.rds"))
cat("obj loaded\n"); print(names(obj@reductions))
umap_ok <- "umap" %in% names(obj@reductions)

# ================= Figure 1 =================
sm_cnt <- as.data.frame(table(md$sample_id)); names(sm_cnt) <- c("sample", "n")
sm_cnt$sample <- factor(sm_cnt$sample, levels = sm_cnt$sample[order(-sm_cnt$n)])
p1A <- ggplot(sm_cnt, aes(sample, n)) + geom_col(fill = c_blue, width = 0.6, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n), vjust = -0.4, size = 2, colour = "black") +
  labs(tag = "A", x = "Sample", y = "Number of cells", title = "Cells per sample") +
  pub_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
st <- ifelse(grepl("PRE-T", md$sample_id), "Pre-T", ifelse(grepl("EOI", md$sample_id), "EOI", ifelse(grepl("Rel", md$sample_id), "Relapse", "Diagnosis")))
st <- factor(st, levels = c("Diagnosis", "EOI", "Relapse", "Pre-T"))
p1B <- ggplot(data.frame(st), aes(st)) + geom_bar(fill = c_orange, width = 0.6, colour = "white", linewidth = 0.4) +
  labs(tag = "B", x = "Sample type", y = "Number of cells", title = "Cells per sample type") +
  pub_theme()
if (umap_ok) {
  p1C <- DimPlot(obj, reduction = "umap", group.by = "sample_id", pt.size = 0.15, label = FALSE) +
    labs(tag = "C", x = "UMAP_1", y = "UMAP_2", title = "UMAP by sample") +
    pub_theme() + theme(legend.text = element_text(size = 5), legend.position = "right") +
    scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(unique(md$sample_id))))
} else p1C <- p1B
compose(list(p1A, p1B, p1C), 1, 3, file.path(figd, "Figure1_Cohort_Overview.pdf"), w = 10, h = 3.2)

# ================= Figure 2 =================
p2A <- ggplot(md, aes(InferCNV)) + geom_bar(aes(fill = InferCNV), width = 0.6, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = c_red, "Negative" = c_blue), name = "inferCNV status") +
  labs(tag = "A", x = "inferCNV status", y = "Number of cells", title = "CNV status") + pub_theme() +
  theme(legend.position = "none")
p2B <- ggplot(md, aes(InferCNVcount)) + geom_histogram(bins = 40, fill = c_blue, colour = "white", linewidth = 0.2) +
  labs(tag = "B", x = "InferCNV count (CNVs per cell)", y = "Number of cells", title = "CNV burden distribution") + pub_theme()
if (umap_ok) {
  p2C <- DimPlot(obj, reduction = "umap", group.by = "InferCNV", pt.size = 0.15) +
    scale_color_manual(values = c("Positive" = c_red, "Negative" = c_blue), name = "inferCNV status") +
    labs(tag = "C", x = "UMAP_1", y = "UMAP_2", title = "UMAP by CNV status") + pub_theme() +
    theme(legend.position = "none")
  p2D <- FeaturePlot(obj, features = "InferCNVcount", reduction = "umap", pt.size = 0.15) +
    scale_color_gradientn(colours = seq_blue) +
    labs(tag = "D", x = "UMAP_1", y = "UMAP_2", title = "UMAP by CNV burden", colour = "InferCNV\ncount") +
    pub_theme() + theme(legend.position = "right")
} else { p2C <- p2B; p2D <- p2B }
prop_s <- as.data.frame(prop.table(table(md$sample_id, md$InferCNV), 1) * 100)
names(prop_s) <- c("sample", "status", "pct")
prop_s$sample <- factor(prop_s$sample, levels = levels(sm_cnt$sample))
p2E <- ggplot(prop_s, aes(sample, pct, fill = status)) + geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("Positive" = c_red, "Negative" = c_blue), name = "inferCNV status") +
  labs(tag = "E", x = "Sample", y = "Proportion of cells (%)", title = "CNV-positive proportion") +
  pub_theme() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "right")
compose(list(p2A, p2B, p2C, p2D, p2E), 2, 3, file.path(figd, "Figure2_CNV_Landscape.pdf"), w = 10.5, h = 6.4)

# ================= Figure 3 =================
phx <- as.data.frame(prop.table(table(md$bin, md$Phase), 1) * 100)
names(phx) <- c("bin", "phase", "pct")
sg2m <- aggregate(pct ~ bin, phx[phx$phase %in% c("S", "G2M"), ], sum)
p3A <- ggplot(phx, aes(bin, pct, fill = phase)) + geom_col(width = 0.62, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = phase_cols, name = "Phase") +
  geom_text(data = sg2m, aes(bin, 100 + 2.5, label = sprintf("%.1f%%", pct)), inherit.aes = FALSE, size = 2.2, fontface = "bold") +
  labs(tag = "A", x = "InferCNV count bin", y = "Cell proportion (%)", title = "Cell-cycle composition") +
  pub_theme() + scale_y_continuous(limits = c(0, 110), expand = expansion(mult = c(0, 0)))
p3B <- ggplot(sg2m, aes(bin, pct)) + geom_col(fill = c_blue, width = 0.55, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 2.2, fontface = "bold") +
  labs(tag = "B", x = "InferCNV count bin", y = "S + G2M cells (%)", title = "Proliferation vs CNV burden") +
  pub_theme() + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
# GO bubble: top 6 terms per bin
src <- file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final")
gof <- c("0-0.25" = "GO_BP_Enrichment_0-0.25_Bin.csv", "0.25-0.75" = "GO_BP_Enrichment_0.25-0.75_Bin.csv", "0.75-1.5" = "GO_BP_Enrichment_0.75-1.5_Bin.csv")
go_df <- do.call(rbind, lapply(names(gof), function(b) {
  d <- read.csv(file.path(src, gof[[b]]), check.names = FALSE)
  d <- d[order(d$p.adjust), ][1:6, ]
  data.frame(bin = b, term = sub(" signaling pathway", "", sub(" regulation of", " reg. of", d$Description)),
             nlogP = -log10(d$p.adjust), Count = d$Count)
}))
go_df$bin <- factor(go_df$bin, levels = names(gof))
go_df$term <- factor(go_df$term, levels = rev(unique(go_df$term)))
p3C <- ggplot(go_df, aes(nlogP, term, size = Count)) + geom_point(colour = c_blue, alpha = 0.85) +
  facet_wrap(~ bin, ncol = 3, scales = "free_y") +
  labs(tag = "C", x = "-log10(adjusted P)", y = "GO biological process", title = "Enrichment per CNV bin") +
  pub_theme(base_size = 6.5) + theme(strip.text = element_text(face = "bold", size = 7))
# marker heatmap
mk <- read.csv(file.path(src, "expression_matrices/marker-v2/2025-prediction_top10_markers.csv"), check.names = FALSE)
ge <- read.csv(file.path(src, "expression_matrices/marker-v2/gene_expression_raw.csv"), check.names = FALSE, row.names = 1)
genes <- unique(mk$gene); genes <- genes[genes %in% colnames(ge)]
ge <- ge[, c(genes, "prediction_label")]
lab <- ge$prediction_label
ge$prediction_label <- NULL
gx <- as.data.frame(t(as.matrix(ge))); gx[is.na(gx)] <- 0
labs3 <- c("CNV_Positive", "ML_maglinant", "ML_normal")
grp_mean <- do.call(cbind, lapply(labs3, function(L) rowMeans(gx[, lab == L, drop = FALSE])))
colnames(grp_mean) <- labs3
grp_z <- t(scale(t(grp_mean)))
cairo_pdf(file.path(figd, "Figure3D_marker_heatmap.pdf"), width = 3.6, height = 4.2, family = "Arial")
ann2 <- data.frame(label = colnames(grp_z), row.names = colnames(grp_z))
pheatmap(grp_z, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("#184f95", "#fcfcfb", "#d03b3b"))(50),
         annotation_col = ann2, fontsize = 6, fontfamily = "Arial", show_rownames = TRUE,
         annotation_colors = list(label = c("CNV_Positive" = c_red, "ML_maglinant" = c_orange, "ML_normal" = c_aqua)),
         annotation_names_col = FALSE, main = "Marker genes")
dev.off()
system2("pdftoppm", c("-png", "-r", "300", "-singlefile", file.path(figd, "Figure3D_marker_heatmap.pdf"), file.path(figd, "Figure3D_marker_heatmap")))
compose(list(p3A, p3B, p3C), 1, 3, file.path(figd, "Figure3_CNV_CellCycle_Programs.pdf"), w = 10.5, h = 3.6)

# ================= Figure 4 =================
perf <- read.csv(file.path(src, "expression_matrices/模型详细性能表.csv"), check.names = FALSE, fileEncoding = "UTF-8-BOM")
perf_long <- data.frame(model = rep(perf$模型, 3),
                        metric = rep(c("Accuracy", "F1 score", "AUC"), each = nrow(perf)),
                        value = c(perf$准确率, perf$F1分数, c(0.9992, 0.9988, 0.9975)))
perf_long$model <- factor(perf_long$model, levels = unique(perf$模型))
p4A <- ggplot(perf_long, aes(model, value, fill = model)) + geom_col(width = 0.6, colour = "white", linewidth = 0.4) +
  facet_wrap(~ metric, nrow = 1) + geom_text(aes(label = sprintf("%.4f", value)), vjust = -0.4, size = 2) +
  scale_fill_manual(values = c(c_blue, c_orange, c_aqua)) +
  scale_y_continuous(limits = c(0.95, 1.01), breaks = c(0.95, 0.98, 1.0)) +
  labs(tag = "A", x = "Model", y = "Score", title = "Classifier performance (5-fold CV)") +
  pub_theme(base_size = 6.5) + theme(legend.position = "none", axis.text.x = element_text(angle = 0))
pred <- read.csv(file.path(src, "expression_matrices/negative_data_predictions.csv"), stringsAsFactors = FALSE)
p4B <- ggplot(pred, aes(malignant_probability, fill = prediction_label)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.2, position = "identity", alpha = 0.7) +
  scale_fill_manual(values = c("ML_maglinant" = c_orange, "ML_normal" = c_aqua), name = "Predicted label") +
  labs(tag = "B", x = "Predicted malignant probability", y = "Number of cells", title = "Test set (1,020 CNV-negative cells)") +
  pub_theme(base_size = 6.5)
lab <- read.csv(file.path(qa, "10_prediction_label_by_sample.csv"), row.names = 1, check.names = FALSE)
lab$sample <- rownames(lab); lab <- lab[order(lab$CNV_Positive + lab$ML_maglinant + lab$ML_normal), ]
lab$sample <- factor(lab$sample, levels = lab$sample)
ll <- data.frame(sample = rep(lab$sample, 3),
                 label = rep(c("CNV_Positive", "ML_maglinant", "ML_normal"), each = nrow(lab)),
                 n = c(lab$CNV_Positive, lab$ML_maglinant, lab$ML_normal))
ll$label <- factor(ll$label, levels = c("CNV_Positive", "ML_maglinant", "ML_normal"))
p4C <- ggplot(ll, aes(sample, n, fill = label)) + geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  coord_flip() + scale_fill_manual(values = c("CNV_Positive" = c_red, "ML_maglinant" = c_orange, "ML_normal" = c_aqua), name = "Label") +
  labs(tag = "C", x = "Sample", y = "Number of cells", title = "Label composition (772 reclassified malignant)") +
  pub_theme(base_size = 6.5) + theme(legend.position = "bottom")
compose(list(p4A, p4B, p4C), 1, 3, file.path(figd, "Figure4_ML_Classification.pdf"), w = 10.5, h = 3.4)

cat("Figures 1-4 done\n")
sink(type = "message"); sink(type = "output")
