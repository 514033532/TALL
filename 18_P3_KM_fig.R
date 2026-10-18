#!/usr/bin/env Rscript
# 18: SuppFig30 — TARGET-ALL-P3 T-lineage KM curves (DHFR, CD74, SOX4, MIF)
suppressMessages({library(survival); library(survminer); library(ggplot2); library(ggpubr)})
qa <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results/TALL_manuscript_package/Analysis_QA"
figd <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results/TALL_manuscript_package/Figures"
kmd <- read.csv(file.path(qa, "61_P3_TALL_KM_data.csv"), check.names = FALSE)
theme_km <- theme_bw(base_size = 8, base_family = "Arial") %+replace%
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "#444444", fill = NA, linewidth = 0.4),
        axis.title = element_text(face = "bold", size = 8, colour = "black"),
        axis.text = element_text(size = 7.5, colour = "black"),
        plot.title = element_text(face = "bold", size = 9, hjust = 0.5, margin = margin(b = 8)),
        legend.position = c(0.82, 0.85), legend.text = element_text(size = 7),
        legend.title = element_blank(), legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(6, 6, 6, 6, "mm"))
mk <- function(g) {
  kmd$grp <- factor(ifelse(kmd[[g]] > median(kmd[[g]]), "High", "Low"), levels = c("High", "Low"))
  fit <- survfit(Surv(time, event) ~ grp, data = kmd)
  pv <- surv_pvalue(fit, data = kmd)$pval
  ggsurvplot(fit, data = kmd, palette = c("#e34948", "#2a78d6"), size = 0.5,
             censor.size = 2, legend.labs = c("High", "Low"),
             xlab = "Time (days)", ylab = "Overall survival",
             title = sprintf("%s (log-rank P = %.3f)", g, pv),
             ggtheme = theme_km, risk.table = FALSE)$plot
}
p1 <- mk("DHFR"); p2 <- mk("CD74"); p3 <- mk("SOX4"); p4 <- mk("MIF")
arr <- ggarrange(plotlist = list(p1, p2, p3, p4), labels = c("A", "B", "C", "D"),
                 ncol = 2, nrow = 2, font.label = list(size = 10, face = "bold", family = "Arial"),
                 label.x = 0.02, label.y = 0.97, hjust = 0, vjust = 1)
f <- file.path(figd, "SuppFig30_TARGET_P3_KM.pdf")
ggsave(f, arr, device = cairo_pdf, width = 8, height = 7)
system2("pdftoppm", c("-png", "-r", "300", "-singlefile", f, sub("\\.pdf$", "", f)))
cat("SuppFig30 written\n")
