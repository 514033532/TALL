#!/usr/bin/env Rscript
# Generate supplementary figures for T-ALL manuscript (dataviz-compliant)
suppressMessages({library(ggplot2); library(pheatmap)})

base <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
pkg  <- file.path(base, "TALL_manuscript_package")
qa   <- file.path(pkg, "Analysis_QA")
figd <- file.path(pkg, "Analysis_QA", "figs")
dir.create(figd, showWarnings = FALSE, recursive = TRUE)

# ---- palette (validated dataviz reference) ----
c_blue   <- "#2a78d6"; c_orange <- "#eb6834"; c_aqua <- "#1baf7a"
c_mag    <- "#e87ba4"; c_red  <- "#e34948"; c_violet <- "#4a3aa7"
ink1 <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"
gridc <- "#e1e0d9"; axisc <- "#c3c2b7"; surf <- "#fcfcfb"
seq_blue <- c("#cde2fb","#b7d3f6","#9ec5f4","#86b6ef","#6da7ec","#5598e7",
              "#3987e5","#2a78d6","#256abf","#1c5cab","#184f95","#104281")

paper_theme <- theme_bw(base_size = 10) +
  theme(panel.background = element_rect(fill = surf), plot.background = element_rect(fill = surf),
        panel.grid.major = element_line(color = gridc, linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = axisc, linewidth = 0.5),
        axis.line = element_blank(), axis.ticks = element_line(color = axisc, linewidth = 0.4),
        text = element_text(color = ink1),
        axis.text = element_text(color = ink2), axis.title = element_text(color = ink1, face = "bold"),
        legend.text = element_text(color = ink2), legend.title = element_text(color = ink1, face = "bold"),
        plot.title = element_text(color = ink1, face = "bold", size = 11),
        strip.background = element_rect(fill = "#f0efec", color = axisc),
        strip.text = element_text(color = ink1, size = 9))

savefig <- function(p, f, w, h) ggsave(file.path(figd, f), p, width = w, height = h, units = "in", dpi = 300)

# ============ data ============
md <- read.csv(gzfile(file.path(qa, "08_cell_annotation_export.csv.gz")), stringsAsFactors = FALSE)
md$bin <- factor(cut(md$InferCNVcount, breaks = c(0, 0.25, 0.75, 1.5, Inf),
                     labels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"),
                     include.lowest = TRUE, right = FALSE),
                 levels = c("0-0.25", "0.25-0.75", "0.75-1.5", "1.5+"))
md$Phase <- factor(md$Phase, levels = c("G1", "S", "G2M"))
phase_cols <- c("G1" = c_blue, "S" = c_orange, "G2M" = c_aqua)

# ---- Supp Fig 1A: phase composition by CNV bin (100% stacked) ----
phx <- as.data.frame(prop.table(table(md$bin, md$Phase), 1) * 100)
names(phx) <- c("bin", "phase", "pct")
prol <- aggregate(pct ~ bin, phx[phx$phase %in% c("S", "G2M"), ], sum)
pA <- ggplot(phx, aes(x = bin, y = pct, fill = phase)) +
  geom_col(width = 0.62, color = surf, linewidth = 0.6) +
  geom_text(data = prol, aes(x = bin, y = 100 + 3, label = sprintf("%.1f%%", pct)),
            inherit.aes = FALSE, size = 3, color = ink2, fontface = "bold") +
  scale_fill_manual(values = phase_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, 108)) +
  labs(x = "InferCNV count bin", y = "Cell proportion (%)",
       fill = "Cell-cycle phase",
       title = "Cell-cycle composition across CNV burden bins (n = 22,679 cells)") +
  paper_theme
savefig(pA, "SuppFig1A_phase_by_CNVbin.pdf", 5.2, 3.4)

# ---- Supp Fig 1B: S+G2M gradient + CNV status ----
sg <- aggregate(Phase %in% c("S", "G2M") ~ bin, md, function(x) 100 * mean(x))
names(sg)[2] <- "pct"
pB <- ggplot(sg, aes(x = bin, y = pct)) +
  geom_col(fill = c_blue, width = 0.55, color = surf, linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.6, size = 3, color = ink1, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "InferCNV count bin", y = "S + G2M cells (%)",
       title = "Proliferating-cell fraction rises with CNV burden\n(chi-square p < 1e-5)") +
  paper_theme
savefig(pB, "SuppFig1B_SG2M_gradient.pdf", 4.2, 3.4)

# ---- Supp Fig 2: pathway x sample heatmap ----
ps <- read.csv(file.path(qa, "20_pathway_strength_malignant.csv"), stringsAsFactors = FALSE)
ph_order <- c("PreT", "Dx", "EOI", "Rel")
samp_order <- c("GSM3872440_PRE-T_1","GSM3872441_PRE-T_2","GSM7091996_T1Dx","GSM7091997_T2Dx",
                "GSM7091998_T3Dx","GSM7091999_T4Dx","GSM7092000_T5Dx","GSM7092001_T6Dx",
                "GSM7092003_T7Dx","GSM7092005_T8Dx","GSM7092007_T9Dx","GSM7092009_T10Dx",
                "GSM7092002_T6EOI","GSM7092004_T7EOI","GSM7092006_T8EOI","GSM7092008_T9EOI",
                "GSM7092010_T10EOI","GSM7092011_T11Rel")
short <- c("PRE-T_1","PRE-T_2","T1Dx","T2Dx","T3Dx","T4Dx","T5Dx","T6Dx","T7Dx","T8Dx","T9Dx","T10Dx",
           "T6EOI","T7EOI","T8EOI","T9EOI","T10EOI","T11Rel")
names(short) <- samp_order
wide <- reshape(ps[, c("sample", "pathway", "prob")], idvar = "pathway", timevar = "sample", direction = "wide")
names(wide) <- sub("prob\\.", "", names(wide))
rownames(wide) <- wide$pathway; wide$pathway <- NULL
wide <- wide[, samp_order]
wide[is.na(wide)] <- 0
# keep pathways with any activity; order by mean Dx
mDx <- rowMeans(wide[, grepl("Dx", colnames(wide)) | colnames(wide) %in% c("GSM7091996_T1Dx","GSM7091997_T2Dx","GSM7091998_T3Dx","GSM7091999_T4Dx","GSM7092000_T5Dx","GSM7092001_T6Dx","GSM7092003_T7Dx","GSM7092005_T8Dx","GSM7092007_T9Dx","GSM7092009_T10Dx")], na.rm = TRUE)
wide <- wide[order(-mDx), ]
keep <- rownames(wide)[rowSums(wide > 0) >= 2]
wide <- wide[keep, ]
logm <- log10(as.matrix(wide) + 1e-3)
phases <- data.frame(Phase = c(rep("PreT", 2), rep("Dx", 10), rep("EOI", 5), "Rel"))
rownames(phases) <- colnames(logm)
ann_cols <- list(Phase = c("PreT" = c_violet, "Dx" = c_blue, "EOI" = c_orange, "Rel" = c_red))
pdf(file.path(figd, "SuppFig2_pathway_heatmap.pdf"), width = 7.5, height = 5.2)
pheatmap(logm, cluster_rows = FALSE, cluster_cols = FALSE,
         color = seq_blue, annotation_col = phases, annotation_colors = ann_cols,
         labels_col = short[colnames(logm)],
         fontsize = 8, fontsize_col = 7, angle_col = 45,
         breaks = seq(-3, 2, length.out = 13),
         main = "log10(1+pathway strength) across samples",
         border_color = surf, gaps_col = c(2, 12, 17))
dev.off()

# ---- Supp Fig 3: longitudinal pathway dynamics (small multiples) ----
top8 <- c("MIF", "GALECTIN", "RESISTIN", "ANNEXIN", "MK", "VISFATIN", "CXCL", "TGFb")
ps2 <- ps[ps$pathway %in% top8, ]
ps2$pathway <- factor(ps2$pathway, levels = top8)
ps2$phase <- factor(ps2$phase, levels = ph_order)
summ <- aggregate(prob ~ pathway + phase, ps2, FUN = function(x) c(m = mean(x), se = sd(x) / sqrt(length(x))))
summ <- do.call(data.frame, summ); names(summ)[3:4] <- c("m", "se")
pC <- ggplot(summ, aes(x = phase, y = m)) +
  geom_col(fill = c_blue, width = 0.6, color = surf, linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(0, m - se), ymax = m + se), width = 0.25, color = ink1, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", m)), vjust = -0.7, size = 2.2, color = ink2) +
  facet_wrap(~ pathway, scales = "free_y", nrow = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Disease phase", y = "Mean pathway strength",
       title = "Communication pathway strength across disease phases") +
  paper_theme + theme(axis.text.x = element_text(angle = 0, size = 8))
savefig(pC, "SuppFig3_longitudinal_pathways.pdf", 7.5, 4.2)

# ---- Supp Fig 4: sender-group contribution ----
snd <- read.csv(file.path(qa, "24_sender_group_total.csv"), stringsAsFactors = FALSE)
snd <- snd[order(snd$prob), ]
snd$source <- factor(snd$source, levels = snd$source)
pD <- ggplot(snd, aes(x = source, y = prob)) +
  geom_col(fill = c_blue, width = 0.6, color = surf, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", prob)), hjust = -0.2, size = 2.6, color = ink1) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Sender subpopulation", y = "Total outgoing communication strength",
       title = "CNV_Positive and subclones dC2/dC0 dominate signal output") +
  paper_theme
savefig(pD, "SuppFig4_sender_groups.pdf", 5, 3.8)

# ---- Supp Fig 5: label composition per sample ----
lab <- read.csv(file.path(qa, "10_prediction_label_by_sample.csv"), stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)
lab$sample <- rownames(lab)
lab <- lab[order(lab$CNV_Positive + lab$ML_maglinant + lab$ML_normal), ]
lab$sample <- factor(lab$sample, levels = lab$sample)
ll <- data.frame(sample = rep(lab$sample, 3),
                 label = rep(c("CNV_Positive", "ML_maglinant", "ML_normal"), each = nrow(lab)),
                 n = c(lab$CNV_Positive, lab$ML_maglinant, lab$ML_normal))
ll$label <- factor(ll$label, levels = c("CNV_Positive", "ML_maglinant", "ML_normal"))
pE <- ggplot(ll, aes(x = sample, y = n, fill = label)) +
  geom_col(width = 0.68, color = surf, linewidth = 0.5) +
  scale_fill_manual(values = c("CNV_Positive" = c_blue, "ML_maglinant" = c_orange, "ML_normal" = c_aqua)) +
  coord_flip() +
  labs(x = "Sample", y = "Number of cells", fill = "Label",
       title = "Cell labeling: 21,659 CNV_Positive; 772 CNV-negative cells reclassified malignant by ML (24.3% normal)") +
  paper_theme + theme(legend.position = "top")
savefig(pE, "SuppFig5_label_composition.pdf", 5.5, 4.5)

cat("Figures written to", figd, "\n")
