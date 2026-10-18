#!/usr/bin/env Rscript
# 07: Figures 5-6 from harmonized CellChat + MRD reanalysis (user figure standards)
suppressMessages({library(ggplot2); library(pheatmap)})
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
seq_blue <- c("#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7","#3987e5","#2a78d6","#256abf","#1c5cab","#184f95","#104281")

compose <- function(plots, nrow, ncol, file, w = 6.5, h = 3) {
  cairo_pdf(file, width = w, height = h, family = "Arial")
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol,
    widths = grid::unit(rep(1/ncol, ncol), "npc"), heights = grid::unit(rep(1/nrow, nrow), "npc"))))
  for (i in seq_along(plots)) {
    rr <- ceiling(i / ncol); cc <- (i - 1) %% ncol + 1
    print(plots[[i]], vp = grid::viewport(layout.pos.row = rr, layout.pos.col = cc))
  }
  grid::popViewport(); dev.off()
  png_file <- sub("\\.pdf$", ".png", file)
  system2("pdftoppm", c("-png", "-r", "300", "-singlefile", file, sub("\\.png$", "", png_file)))
  cat("wrote", file, "\n")
}

# ---------- harmonized CellChat data ----------
pw <- read.csv(file.path(qa, "34_harmonized_pathway_strength.csv"), stringsAsFactors = FALSE)
phase_map <- c("GSM3872440" = "PreT", "GSM3872441" = "PreT", "GSM7091996" = "Dx", "GSM7091997" = "Dx",
               "GSM7091998" = "Dx", "GSM7091999" = "Dx", "GSM7092000" = "Dx", "GSM7092001" = "Dx",
               "GSM7092003" = "Dx", "GSM7092005" = "Dx", "GSM7092007" = "Dx", "GSM7092009" = "Dx",
               "GSM7092002" = "EOI", "GSM7092004" = "EOI", "GSM7092006" = "EOI", "GSM7092008" = "EOI",
               "GSM7092010" = "EOI", "GSM7092011" = "Rel")
pw$phase <- factor(unname(phase_map[pw$sample]), levels = c("PreT", "Dx", "EOI", "Rel"))
samp_ord <- c("GSM3872440","GSM3872441","GSM7091996","GSM7091997","GSM7091998","GSM7091999","GSM7092000",
              "GSM7092001","GSM7092003","GSM7092005","GSM7092007","GSM7092009","GSM7092002","GSM7092004",
              "GSM7092006","GSM7092008","GSM7092010","GSM7092011")
short <- c("PRE-T_1","PRE-T_2","T1Dx","T2Dx","T3Dx","T4Dx","T5Dx","T6Dx","T7Dx","T8Dx","T9Dx","T10Dx",
           "T6EOI","T7EOI","T8EOI","T9EOI","T10EOI","T11Rel"); names(short) <- samp_ord

# Fig 5A: pathway x sample heatmap (harmonized)
avail <- intersect(samp_ord, unique(pw$sample))
wide <- reshape(pw[pw$sample %in% avail, c("sample", "pathway_name", "prob")], idvar = "pathway_name", timevar = "sample", direction = "wide")
names(wide) <- sub("prob\\.", "", names(wide)); rownames(wide) <- wide$pathway_name; wide$pathway_name <- NULL
wide <- wide[, avail]; wide[is.na(wide)] <- 0
wide <- wide[rowSums(wide > 0) >= 1, , drop = FALSE]
logm <- log10(as.matrix(wide) + 1e-3)
phases <- data.frame(Phase = unname(phase_map[colnames(logm)]), row.names = colnames(logm))
gaps <- which(!duplicated(phases$Phase[order(match(phases$Phase, c("PreT","Dx","EOI","Rel")))])) - 1
gaps <- gaps[gaps > 0]
cairo_pdf(file.path(figd, "Figure5A_harmonized_heatmap.pdf"), width = 7, height = 4.6, family = "Arial")
pheatmap(logm, cluster_rows = FALSE, cluster_cols = FALSE, color = seq_blue,
         annotation_col = phases, labels_col = short[colnames(logm)],
         annotation_colors = list(Phase = c("PreT" = c_violet, "Dx" = c_blue, "EOI" = c_orange, "Rel" = c_red)),
         fontsize = 7, fontsize_col = 6.5, angle_col = 45, fontfamily = "Arial",
         breaks = seq(-3, 2, length.out = 13), border_color = "#fcfcfb",
         main = "log10(1+pathway strength), harmonized CNV-bin identities", gaps_col = if (length(gaps) > 0) gaps else NULL)
dev.off()
system2("pdftoppm", c("-png", "-r", "300", "-singlefile", file.path(figd, "Figure5A_harmonized_heatmap.pdf"), file.path(figd, "Figure5A_harmonized_heatmap")))

# Fig 5B/C: sender (bin) contribution + sender x pathway
gp <- read.csv(file.path(qa, "36_harmonized_sender_pathway.csv"), stringsAsFactors = FALSE)
snd <- aggregate(prob ~ source, gp, sum); snd <- snd[order(snd$prob), ]
snd$source <- factor(snd$source, levels = snd$source)
p5B <- ggplot(snd, aes(source, prob)) + geom_col(fill = c_blue, width = 0.6, colour = "white", linewidth = 0.4) +
  coord_flip() + geom_text(aes(label = sprintf("%.1f", prob)), hjust = -0.15, size = 2.4, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(tag = "B", x = "Sender (CNV burden bin)", y = "Total outgoing strength", title = "Signal output by CNV bin") +
  pub_theme()
top_pw <- names(sort(tapply(gp$prob, gp$pathway_name, sum), decreasing = TRUE))[1:6]
gpt <- gp[gp$pathway_name %in% top_pw, ]
gpt$pathway_name <- factor(gpt$pathway_name, levels = top_pw)
p5C <- ggplot(gpt, aes(source, prob, fill = pathway_name)) + geom_col(position = "stack", width = 0.6, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(c_blue, c_orange, c_aqua, c_violet, c_red, "#e87ba4"), name = "Pathway") +
  labs(tag = "C", x = "Sender (CNV burden bin)", y = "Pathway strength", title = "Top pathways by sender bin") +
  pub_theme(base_size = 6.5) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
compose(list(p5B, p5C), 1, 2, file.path(figd, "Figure5BC_sender_contribution.pdf"), w = 7.5, h = 3)

# ---------- Figure 6: longitudinal + MRD ----------
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
  facet_wrap(~ pathway, scales = "free_y", nrow = 1) +
  labs(tag = "A", x = "Disease phase", y = "Mean pathway strength", title = "Communication by phase (subcluster design)") +
  pub_theme(base_size = 6.5) + theme(axis.text.x = element_text(size = 6))

per_s <- read.csv(file.path(qa, "31_MRD_EOI_per_sample.csv"), stringsAsFactors = FALSE)
per_s$MRD_grp <- factor(ifelse(grepl("MRD\\+", per_s$MRD), "MRD+", "MRD-"), levels = c("MRD+", "MRD-"))
per_s$GSM <- factor(per_s$GSM, levels = per_s$GSM[order(per_s$MRD_grp, per_s$GSM)])
ps_long <- data.frame(GSM = rep(per_s$GSM, 2), MRD = rep(per_s$MRD_grp, 2),
                      metric = rep(c("pct_CNV_Pos", "pct_SG2M"), each = nrow(per_s)),
                      value = c(per_s$pct_CNV_Pos, per_s$pct_SG2M))
ps_long$metric <- factor(ps_long$metric, levels = c("pct_CNV_Pos", "pct_SG2M"),
                         labels = c("CNV-positive cells (%)", "S+G2M cells (%)"))
p6B <- ggplot(ps_long, aes(GSM, value, fill = MRD)) + geom_col(width = 0.6, colour = "white", linewidth = 0.4) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("MRD+" = c_red, "MRD-" = c_blue), name = "EOI MRD") +
  labs(tag = "B", x = "EOI sample", y = "Value (%)", title = "CNV burden and proliferation by MRD") +
  pub_theme(base_size = 6.5) + theme(axis.text.x = element_text(angle = 0, size = 6))

gexp <- read.csv(file.path(qa, "32_MRD_EOI_gene_expression.csv"), stringsAsFactors = FALSE)
gexp$MRD_grp <- factor(ifelse(grepl("MRD\\+", gexp$MRD), "MRD+", "MRD-"), levels = c("MRD+", "MRD-"))
g4 <- c("mean_PDCD1", "mean_HAVCR2", "mean_MIF", "mean_LGALS9", "mean_CD74", "mean_SOX4")
g4n <- c("PD-1", "TIM-3", "MIF", "LGALS9", "CD74", "SOX4")
gg <- do.call(rbind, lapply(seq_along(g4), function(i) {
  data.frame(GSM = gexp$GSM, MRD = gexp$MRD_grp, gene = g4n[i], value = gexp[[g4[i]]])
}))
gg$gene <- factor(gg$gene, levels = g4n)
p6C <- ggplot(gg, aes(gene, value, fill = MRD)) + geom_boxplot(width = 0.5, outlier.size = 0.5, linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.75), size = 0.8, colour = "black") +
  scale_fill_manual(values = c("MRD+" = c_red, "MRD-" = c_blue), name = "EOI MRD") +
  labs(tag = "C", x = "Gene", y = "Mean expression (CP10K)", title = "Key genes by MRD (EOI)") +
  pub_theme(base_size = 6.5)
compose(list(p6A, p6B, p6C), 1, 3, file.path(figd, "Figure6_Longitudinal_MRD.pdf"), w = 10.5, h = 3.3)
cat("Figures 5-6 done\n")
