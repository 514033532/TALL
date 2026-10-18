#!/usr/bin/env Rscript
# 16: TARGET/MP2PRT T-ALL identification by expression lineage + survival validation of
#     MIF/LGALS9/TYMS/DHFR/CDKN2A/SOX4/CD74/HAVCR2 (OS, log-rank + Cox)
suppressMessages({library(Seurat); library(survival); library(survminer); library(org.Hs.eg.db)})
base <- "/media/zhh2/Newsmy/TALL研究汇总"
pkg  <- file.path(base, "analysis_results/TALL_manuscript_package")
qa   <- file.path(pkg, "Analysis_QA")
B    <- "/media/zhh2/00f14a6a-b492-4d39-bd84-bf715b467ad8/home/zhh/10genome/data-ALL"
log <- file(file.path(qa, "16_TARGET_log.txt"), open = "wt")
sink(log, type = "output"); sink(log, type = "message")
t0 <- Sys.time()

## 1. Load expression + clinical
e <- new.env(); load(file.path(B, "gdc-atcc-all/MP2PRT/MP2PRT-202407.Rdata"), envir = e)
pbmc <- e$pbmc
gids <- sub("\\.\\d+$", "", rownames(pbmc))
sym <- AnnotationDbi::select(org.Hs.eg.db, keys = gids, columns = "SYMBOL", keytype = "ENSEMBL")
sym <- sym[!is.na(sym$SYMBOL) & !duplicated(sym$ENSEMBL), ]
keep <- gids %in% sym$ENSEMBL
pbmc <- pbmc[keep, ]
rownames(pbmc) <- sym$SYMBOL[match(gids[keep], sym$ENSEMBL)]
pbmc <- pbmc[!duplicated(rownames(pbmc)), ]
cat("pbmc:", dim(pbmc), "\n")
clin <- read.delim(file.path(B, "gdc-atcc-all/MP2PRT/clinical.project-mp2prt-all.2024-07-24/clinical.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
ss <- read.delim(file.path(B, "gdc-atcc-all/MP2PRT/gdc_sample_sheet.2024-08-02-gene-count.tsv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
cat("clinical:", nrow(clin), "| sample sheet:", nrow(ss), "\n")

# map sample -> case
sm2case <- setNames(ss$`Case ID`, ss$`File Name`)
cells <- colnames(pbmc)
case_of <- sm2case[cells]
cat("matched cases:", sum(!is.na(case_of)), "\n")

## 2. Lineage calling by expression
genes <- rownames(pbmc)
t_mk <- intersect(c("CD3D","CD3E","CD3G","CD2","CD7","TRAC","TRBC1","TRBC2","LCK","LAT","TCF7","LEF1","CD5","CD1A"), genes)
b_mk <- intersect(c("CD19","CD79A","PAX5","MS4A1","VPREB1","EBF1","CD22","CD24","IGLL1","CD79B"), genes)
cat("T markers:", length(t_mk), "| B markers:", length(b_mk), "\n")
# counts or normalized? use counts -> log1p CPM
cts <- GetAssayData(pbmc, layer = "counts")
lib <- colSums(cts); lib[lib == 0] <- 1
tscore <- rowMeans(log1p(t(cts[t_mk, , drop = FALSE]) / lib * 1e6))
bscore <- rowMeans(log1p(t(cts[b_mk, , drop = FALSE]) / lib * 1e6))
lineage <- ifelse(tscore - bscore > 1.2 & tscore > 1.2, "T",
           ifelse(bscore - tscore > 1.2 & bscore > 1.2, "B", "Ambiguous"))
cat("expression lineage:\n"); print(table(lineage))

# clinical morphology join
clin$morphology[is.na(clin$morphology)] <- ""
mor <- setNames(clin$morphology, clin$case_submitter_id)
mor[is.na(mor)] <- ""
case_mor <- mor[case_of]
case_mor[is.na(case_mor)] <- ""
is_BALL_byMorph <- grepl("9836|9811", case_mor)
tab <- table(morph_BALL = is_BALL_byMorph, expr_lineage = lineage)
print(tab)
cat("morphology-B-ALL 中被表达判定为 B 的比例:",
    round(100 * sum(is_BALL_byMorph & lineage == "B") / sum(is_BALL_byMorph), 1), "%\n")

## 3. T-ALL set: expression-T lineage (regardless of morphology)
tall_cells <- cells[lineage == "T"]
cat("T-ALL 样本数:", length(tall_cells), "\n")
# dedupe per case (prefer BM samples)
tall_df <- data.frame(cell = tall_cells, case = case_of[tall_cells],
                      tscore = tscore[tall_cells], bscore = bscore[tall_cells], stringsAsFactors = FALSE)
tall_df$stype <- ss$`Sample Type`[match(tall_df$cell, ss$`File Name`)]
stype_rank <- c("Blood Derived Cancer - Bone Marrow" = 1, "Blood Derived Cancer - Peripheral Blood" = 2)
tall_df$rank <- ifelse(is.na(stype_rank[tall_df$stype]), 3, stype_rank[tall_df$stype])
tall_df <- tall_df[order(tall_df$case, tall_df$rank, -tall_df$tscore), ]
tall_df <- tall_df[!duplicated(tall_df$case), ]
cat("去重后 T-ALL 病例数:", nrow(tall_df), "\n")

## 4. Survival data
surv <- data.frame(case = clin$case_submitter_id,
                   os_time = ifelse(clin$vital_status == "Dead", clin$days_to_death, clin$days_to_last_follow_up),
                   os_event = as.integer(clin$vital_status == "Dead"),
                   age = clin$age_at_diagnosis, gender = clin$gender, stringsAsFactors = FALSE)
surv$os_time <- as.numeric(surv$os_time)
surv <- surv[!is.na(surv$os_time) & surv$os_time > 0, ]
tall_df <- merge(tall_df, surv, by = "case", all.x = TRUE)
tall_df <- tall_df[!is.na(tall_df$os_time), ]
cat("有生存数据的 T-ALL:", nrow(tall_df), "| 事件数:", sum(tall_df$os_event), "\n")

## 5. Gene survival analysis (median split, log-rank + Cox)
expr <- GetAssayData(pbmc, layer = "counts")
genes_ck <- c("MIF","LGALS9","TYMS","DHFR","CDKN2A","SOX4","CD74","HAVCR2")
genes_ck <- genes_ck[genes_ck %in% rownames(expr)]
gex <- t(as.matrix(expr[genes_ck, tall_df$cell]))
gex <- log1p(t(t(gex) / colSums(expr)[tall_df$cell] * 1e6))
colnames(gex) <- genes_ck
res <- do.call(rbind, lapply(genes_ck, function(g) {
  hi <- gex[, g] > median(gex[, g])
  df <- data.frame(time = tall_df$os_time, event = tall_df$os_event, hi = hi,
                   age = tall_df$age, gender = tall_df$gender)
  fit <- survfit(Surv(time, event) ~ hi, data = df)
  lr <- survdiff(Surv(time, event) ~ hi, data = df)
  p <- 1 - pchisq(lr$chisq, 1)
  cx <- coxph(Surv(time, event) ~ hi, data = df)
  cx_adj <- tryCatch(coxph(Surv(time, event) ~ hi + age + gender, data = df), error = function(e) NULL)
  data.frame(gene = g, n = nrow(df), events = sum(df$event),
             HR = exp(coef(cx))[1], HR_low = exp(confint(cx))[1, 1], HR_high = exp(confint(cx))[1, 2],
             logrank_p = p,
             HR_adj = if (!is.null(cx_adj)) exp(coef(cx_adj))[1] else NA,
             adj_p = if (!is.null(cx_adj)) summary(cx_adj)$coefficients[1, 5] else NA)
}))
write.csv(res, file.path(qa, "56_TALL_survival_analysis.csv"), row.names = FALSE)
print(res)

## 6. T-ALL vs B-ALL expression comparison
ball_df <- data.frame(cell = cells[lineage == "B"], case = case_of[lineage == "B"], stringsAsFactors = FALSE)
ball_df <- ball_df[!duplicated(ball_df$case), ]
gexB <- t(as.matrix(expr[genes_ck, ball_df$cell]))
gexB <- log1p(t(t(gexB) / colSums(expr)[ball_df$cell] * 1e6))
cmp <- do.call(rbind, lapply(genes_ck, function(g) {
  w <- wilcox.test(gex[, g], gexB[, g])
  data.frame(gene = g, mean_T = mean(gex[, g]), mean_B = mean(gexB[, g]),
             fold = mean(gex[, g]) - mean(gexB[, g]), wilcox_p = w$p.value)
}))
write.csv(cmp, file.path(qa, "57_TALL_vs_BALL_gene_expression.csv"), row.names = FALSE)
print(cmp)

## 7. KM data for figure
write.csv(data.frame(case = tall_df$case, cell = tall_df$cell, time = tall_df$os_time,
                     event = tall_df$os_event, gex, check.names = FALSE),
          file.path(qa, "58_TALL_survival_KM_data.csv"), row.names = FALSE)

sink(type = "message"); sink(type = "output")
cat("DONE in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
