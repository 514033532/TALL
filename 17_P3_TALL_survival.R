#!/usr/bin/env Rscript
# 17: TARGET-ALL-P3 (downloaded subset, 135 samples): expression lineage + survival join
suppressMessages({library(Seurat); library(survival)})
B <- "/media/zhh2/00f14a6a-b492-4d39-bd84-bf715b467ad8/home/zhh/10genome/data-ALL/gdc-atcc-all/TARGET-P3"
qa <- "/media/zhh2/Newsmy/TALL研究汇总/analysis_results/TALL_manuscript_package/Analysis_QA"
log <- file(file.path(qa, "17_P3_log.txt"), open = "wt")
sink(log, type = "output"); sink(log, type = "message")
e <- new.env(); load(file.path(B, "pbmc-targetp3-fpkm.Rdata"), envir = e)
pbmc <- e$pbmc
cts0 <- GetAssayData(pbmc, layer = "counts")
colnames(cts0) <- colnames(pbmc)
pbmc <- CreateSeuratObject(counts = cts0)
cat("pbmc:", dim(pbmc), "\n")
# FPKM data slot
mat <- GetAssayData(pbmc, layer = "data")
if (nrow(mat) == 0) mat <- GetAssayData(pbmc, layer = "counts")
genes <- rownames(pbmc)
t_mk <- intersect(c("CD3D","CD3E","CD3G","CD2","CD7","TRAC","TRBC1","TRBC2","LCK","LAT","TCF7","LEF1","CD5","CD1A"), genes)
b_mk <- intersect(c("CD19","CD79A","PAX5","MS4A1","VPREB1","EBF1","CD22","CD24","IGLL1","CD79B"), genes)
ts <- colMeans(log1p(as.matrix(mat[t_mk, , drop = FALSE])))
bs <- colMeans(log1p(as.matrix(mat[b_mk, , drop = FALSE])))
lin <- ifelse(ts - bs > 0.8 & ts > 0.5, "T",
       ifelse(bs - ts > 0.8 & bs > 0.5, "B", "Ambiguous"))
print(table(lin))
# map File Name -> Sample ID -> survival
ss <- read.delim(file.path(B, "gdc_sample_sheet.2020-11-11.tsv"), check.names = FALSE)
sm2smp <- setNames(ss$`Sample ID`, ss$`File ID`)
smp <- sm2smp[colnames(pbmc)]
cat("mapped to Sample ID:", sum(!is.na(smp)), "\n")
surv <- read.delim(file.path(B, "TARGET-ALL-P3.survival.tsv"), check.names = FALSE)
surv$pat <- sub("-0[0-9][A-Z].*$", "", surv$sample)   # case id
sv <- data.frame(sample = surv$sample, OS = surv$OS, OS.time = surv$OS.time, stringsAsFactors = FALSE)
df <- data.frame(cell = colnames(pbmc), sample = smp, lineage = lin, tscore = ts, bscore = bs, stringsAsFactors = FALSE)
df <- merge(df, sv, by = "sample", all.x = TRUE)
cat("有生存数据:", sum(!is.na(df$OS.time)), "| T 谱系且有生存:", sum(df$lineage == "T" & !is.na(df$OS.time)), "\n")
print(df[df$lineage == "T", c("sample", "tscore", "OS", "OS.time")])
write.csv(df, file.path(qa, "59_P3_lineage_survival.csv"), row.names = FALSE)

## --- gene survival on T-lineage subset ---
tt <- df[df$lineage == "T" & !is.na(df$OS.time) & df$OS.time > 0, ]
cat("T-lineage survival subset n =", nrow(tt), "| events =", sum(tt$OS), "\n")
genes_ck <- intersect(c("MIF","LGALS9","TYMS","DHFR","CDKN2A","SOX4","CD74","HAVCR2"), genes)
cat("DEBUG mat dim:", paste(dim(mat), collapse="x"), "| genes_ck:", length(genes_ck), "| tt cells:", length(tt$cell), "\n")
cat("DEBUG tt$cell head:", paste(head(tt$cell, 3), collapse=" | "), "\n")
cat("DEBUG cell in colnames(mat):", sum(tt$cell %in% colnames(mat)), "| gene in rownames:", sum(genes_ck %in% rownames(mat)), "\n")
gex <- as.matrix(mat[genes_ck, tt$cell])
res <- do.call(rbind, lapply(genes_ck, function(g) {
  hi <- gex[g, ] > median(gex[g, ])
  dd <- data.frame(time = tt$OS.time, event = tt$OS, hi = hi)
  lr <- survdiff(Surv(time, event) ~ hi, data = dd)
  p <- 1 - pchisq(lr$chisq, 1)
  cx <- coxph(Surv(time, event) ~ hi, data = dd)
  data.frame(gene = g, n = nrow(dd), events = sum(dd$event),
             HR = exp(coef(cx))[1], HR_low = exp(confint(cx))[1,1], HR_high = exp(confint(cx))[1,2],
             logrank_p = p)
}))
write.csv(res, file.path(qa, "60_P3_TALL_survival_analysis.csv"), row.names = FALSE)
print(res)
kmd <- data.frame(sample = tt$sample, time = tt$OS.time, event = tt$OS,
                  t(gex), check.names = FALSE)
write.csv(kmd, file.path(qa, "61_P3_TALL_KM_data.csv"), row.names = FALSE)
sink(type = "message"); sink(type = "output")
cat("DONE\n")
