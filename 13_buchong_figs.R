#!/usr/bin/env Rscript
# 13: Regenerate the 3 keepable supplementary figures (SuppFig18/19/20) to SCI standards:
#     Arial, English, bold axis titles, unified color scales, label colors
#     CNV_Positive=red, ML_maglinant=orange, ML_normal=blue (user-specified)
suppressMessages({library(Seurat); library(ggplot2); library(ggpubr)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
figd <- file.path(pkg, "Figures", "buchong")
c_blue <- "#2a78d6"; c_orange <- "#eb6834"; c_aqua <- "#1baf7a"; c_red <- "#e34948"
lab_cols <- c("CNV_Positive" = c_red, "ML_maglinant" = c_orange, "ML_normal" = c_blue)

pub_theme <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "Arial") %+replace%
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "#444444", fill = NA, linewidth = 0.4),
          axis.title = element_text(face = "bold", size = base_size, colour = "black"),
          axis.title.y = element_text(angle = 90, margin = margin(r = 6)),
          axis.title.x = element_text(angle = 0, margin = margin(t = 6)),
          axis.text  = element_text(size = base_size - 0.5, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0.5, margin = margin(b = 8)),
          legend.title = element_text(face = "bold", size = base_size),
          legend.text  = element_text(size = base_size - 0.5),
          strip.text = element_text(face = "bold", size = base_size - 0.5),
          plot.margin = margin(6, 6, 6, 6, "mm"))
}
save_png <- function(pdf_file) {
  system2("pdftoppm", c("-png", "-r", "300", "-singlefile", pdf_file, sub("\\.pdf$", "", pdf_file)))
}

obj <- readRDS(file.path(base, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/marker-v2/combined_TALL_22679_prediction_label_seurat.rds"))
obj$prediction_label <- factor(obj$prediction_label, levels = c("CNV_Positive", "ML_maglinant", "ML_normal"))

genes_core <- c("TYMS", "MKI67", "UHRF1", "MIF", "LGALS9", "CD74", "SOX4", "CCL5")

# ============ SuppFig18: Feature plots on UMAP (unified scale 0-4) ============
fps <- lapply(genes_core, function(g) {
  FeaturePlot(obj, features = g, reduction = "umap", pt.size = 0.25, order = TRUE,
              min.cutoff = 0, max.cutoff = 4) +
    scale_color_gradientn(colours = c("#cde2fb", "#5598e7", "#104281"), limits = c(0, 4)) +
    labs(title = g, x = "UMAP_1", y = "UMAP_2", colour = "Expression\n(CP10K)") +
    pub_theme(7)
})
arr18 <- ggarrange(plotlist = fps, labels = LETTERS[1:8], ncol = 4, nrow = 2,
                   font.label = list(size = 10, face = "bold", family = "Arial"),
                   label.x = 0.02, label.y = 0.97, hjust = 0, vjust = 1,
                   common.legend = TRUE, legend = "bottom")
f18 <- file.path(figd, "SuppFig18_Feature_Plot_UMAP.pdf")
ggsave(f18, arr18, device = cairo_pdf, width = 11, height = 6); save_png(f18)
cat("wrote SuppFig18\n")

# ============ SuppFig19: Dot plot (uniform scale) ============
dp <- DotPlot(obj, features = genes_core, group.by = "prediction_label",
              cols = c("#cde2fb", "#104281"), col.min = 0, col.max = 4, dot.scale = 7) +
  scale_color_gradientn(colours = c("#cde2fb", "#5598e7", "#104281"), limits = c(0, 4),
                        name = "Average\nexpression") +
  labs(x = "Gene", y = "Cell label", title = "Marker expression across cellular states") +
  pub_theme(8) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
                       plot.margin = margin(6, 6, 10, 6, "mm"))
f19 <- file.path(figd, "SuppFig19_Marker_Dot_Plot.pdf")
ggsave(f19, dp, device = cairo_pdf, width = 6.5, height = 3.6); save_png(f19)
cat("wrote SuppFig19\n")

# ============ SuppFig20: Violin plots (9 genes x 3 labels, log1p) ============
genes_v <- c("TYMS", "MKI67", "UHRF1", "MIF", "LGALS9", "CD74", "SOX4", "PRF1", "S100A8")
vdf <- do.call(rbind, lapply(genes_v, function(g) {
  e <- GetAssayData(obj, layer = "data")[g, ]
  data.frame(gene = factor(g, levels = genes_v), label = obj$prediction_label,
             expr = as.numeric(log1p(e)))
}))
vp <- ggplot(vdf, aes(label, expr, fill = label)) +
  geom_violin(scale = "width", linewidth = 0.3, colour = "white") +
  geom_boxplot(width = 0.12, outlier.size = 0.3, linewidth = 0.25, colour = "black", fill = "white", alpha = 0.8) +
  facet_wrap(~ gene, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = lab_cols, name = "Cell label") +
  scale_x_discrete(labels = c("CNV+\n(malignant)", "ML\nmalignant", "ML\nnormal")) +
  labs(x = "Cell label", y = "log1p(expression)", title = "Marker gene expression by cellular state") +
  pub_theme(7) + theme(legend.position = "bottom", axis.text.x = element_text(size = 6))
f20 <- file.path(figd, "SuppFig20_Marker_Violin.pdf")
ggsave(f20, vp, device = cairo_pdf, width = 8, height = 6.5); save_png(f20)
cat("wrote SuppFig20\n")
cat("ALL BUCHONG FIGS DONE\n")
