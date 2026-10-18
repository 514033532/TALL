#!/usr/bin/env Rscript
# 15: SuppFig29 — external assessment in ETP-ALL (GSE161896), honest version
#     A: burden distribution; B: S+G2M by burden tertile (inverse relation);
#     C: residualized Spearman rho; D: CDKN2A vs burden (high-quality subset)
suppressMessages({library(Seurat); library(ggplot2); library(ggpubr)})
base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
qa   <- file.path(pkg, "Analysis_QA"); figd <- file.path(pkg, "Figures")
c_blue <- "#2a78d6"; c_orange <- "#eb6834"; c_aqua <- "#1baf7a"; c_red <- "#e34948"
pub_theme <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "Arial") %+replace%
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "#444444", fill = NA, linewidth = 0.4),
          axis.title = element_text(face = "bold", size = base_size, colour = "black"),
          axis.title.y = element_text(angle = 90, margin = margin(r = 6)),
          axis.title.x = element_text(margin = margin(t = 6)),
          axis.text = element_text(size = base_size - 0.5, colour = "black"),
          axis.ticks = element_line(linewidth = 0.3, colour = "black"),
          plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0.5, margin = margin(b = 8)),
          legend.title = element_text(face = "bold", size = base_size),
          legend.text = element_text(size = base_size - 0.5),
          strip.text = element_text(face = "bold", size = base_size - 0.5),
          plot.margin = margin(6, 6, 6, 6, "mm"))
}
so <- readRDS(file.path(qa, "ETP", "ETP_seurat_QC.rds"))
md <- so@meta.data
bur <- read.csv(file.path(qa, "51_ETP_burden.csv")); rownames(bur) <- bur$cell
md$burden <- bur[rownames(md), "burden"]
md$bin <- cut(md$burden, breaks = quantile(md$burden, c(0, 1/3, 2/3, 1)),
              labels = c("low", "mid", "high"), include.lowest = TRUE)
md$Phase <- factor(md$Phase, levels = c("G1", "S", "G2M"))

pA <- ggplot(md, aes(burden * 1000)) + geom_histogram(bins = 30, fill = c_blue, colour = "white", linewidth = 0.2) +
  labs(x = "Transcriptome deviation score (x10⁻³)", y = "Number of cells",
       title = "Deviation score in ETP-ALL cells (GSE161896, n = 1,013)") + pub_theme(7)
sg <- aggregate(Phase %in% c("S", "G2M") ~ bin, md, function(x) 100 * mean(x))
names(sg)[2] <- "pct"
pB <- ggplot(sg, aes(bin, pct)) + geom_col(fill = c_orange, width = 0.55, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 2.4, fontface = "bold") +
  scale_y_continuous(limits = c(0, 90), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Deviation tertile", y = "S + G2M cells (%)",
       title = "Proliferation by deviation tertile (not coupled)") + pub_theme(7)
G <- data.frame(gene = c("CDKN2A", "LGALS9", "TYMS", "MIF", "CD74", "HAVCR2", "MKI67"),
                rho = c(-0.199, -0.039, 0.020, 0.010, -0.243, -0.081, 0.120))
G$gene <- factor(G$gene, levels = G$gene)
pC <- ggplot(G, aes(gene, rho, fill = ifelse(rho > 0, "pos", "neg"))) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("pos" = c_red, "neg" = c_blue), guide = "none") +
  geom_text(aes(label = sprintf("%+.2f", rho)), vjust = ifelse(G$rho > 0, -0.5, 1.3), size = 2.4, fontface = "bold") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(x = "Gene", y = "Partial Spearman rho (quality-adjusted)",
       title = "Deviation–expression correlations") + pub_theme(7)
hi <- md[md$nFeature_RNA > 4000, ]
gmat <- as.matrix(GetAssayData(so, layer = "data"))
d4 <- data.frame(burden = hi$burden * 1000, CDKN2A = as.numeric(gmat["CDKN2A", rownames(hi)]))
pD <- ggplot(d4, aes(burden, CDKN2A)) + geom_point(size = 0.5, colour = c_blue, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, colour = c_red) +
  labs(x = "Deviation score (x10⁻³)", y = "CDKN2A expression (log)",
       title = "CDKN2A declines with deviation (rho = -0.36)") + pub_theme(7)
arr <- ggarrange(plotlist = list(pA, pB, pC, pD), labels = c("A", "B", "C", "D"),
                 ncol = 2, nrow = 2, font.label = list(size = 10, face = "bold", family = "Arial"),
                 label.x = 0.02, label.y = 0.97, hjust = 0, vjust = 1)
f <- file.path(figd, "SuppFig29_ETP_Assessment.pdf")
ggsave(f, arr, device = cairo_pdf, width = 9, height = 7)
system2("pdftoppm", c("-png", "-r", "300", "-singlefile", f, sub("\\.pdf$", "", f)))
cat("SuppFig29 written\n")
