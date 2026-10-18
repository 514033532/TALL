#!/usr/bin/env Rscript
# 14: External validation in ETP-ALL single cells (GSE161896, SMART-seq2, hg19)
#     Pipeline: QC -> ENSG->SYMBOL -> cell-cycle scoring -> inferCNV (refs = normal T cells)
#     -> CNV bins -> LGALS9/MIF/TYMS/CDKN2A gradients -> correlations.
#     Plus quick pass over GSE168594 (4x 10x h5, human/mouse mixed PDX).
suppressMessages({library(Seurat); library(infercnv); library(data.table); library(org.Hs.eg.db)})
base <- "/media/zhh2/Newsmy/TALL研究汇总"
pkg  <- file.path(base, "analysis_results/TALL_manuscript_package")
qa   <- file.path(pkg, "Analysis_QA")
dir.create(file.path(qa, "ETP"), showWarnings = FALSE, recursive = TRUE)
log <- file(file.path(qa, "14_ETP_validation_log.txt"), open = "wt")
sink(log, type = "output"); sink(log, type = "message")

t0 <- Sys.time()
## ============ 1. Read htseq files (sparse) ============
fs <- list.files(file.path(base, "GSE161896_RAW"), pattern = "\\.htseq-count\\.txt\\.gz$", full.names = TRUE)
cat("files:", length(fs), "\n")
cell_ids <- sub("\\.htseq-count\\.txt\\.gz$", "", basename(fs))
sample   <- sub("^GSM[0-9]+_", "", cell_ids)
sample   <- sub("_SingleCell_.*", "", sample)
sample   <- sub("_S[0-9]+$", "", sample)
grp      <- ifelse(grepl("^ETP", sample), "Patient_ETP", "PDX_TALL")
cat("groups:\n"); print(table(sample))

genes <- character(); rows_list <- list()
# first pass: collect all ENSG ids present
allg <- character()
for (f in fs[seq(1, length(fs), by = 50)]) {
  d <- tryCatch(fread(f, header = FALSE, sep = "\t"), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !"V2" %in% names(d)) next
  allg <- union(allg, d$V1[!grepl("^__", d$V1)])
}
cat("union genes:", length(allg), "\n")
gi <- setNames(seq_along(allg), allg)

mat <- Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(),
                            dims = c(length(allg), length(fs)),
                            dimnames = list(allg, cell_ids))
for (k in seq_along(fs)) {
  d <- tryCatch(fread(fs[k], header = FALSE, sep = "\t"), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !"V2" %in% names(d)) { cat("SKIP bad file:", basename(fs[k]), "\n"); next }
  d <- d[!grepl("^__", d$V1)][V2 > 0]
  if (nrow(d) > 0) mat[gi[d$V1], k] <- d$V2
  if (k %% 200 == 0) cat("read", k, "files\n")
}
cat("matrix built:", dim(mat), "in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

## ============ 2. ENSG -> SYMBOL ============
ensg <- rownames(mat)
sym  <- AnnotationDbi::select(org.Hs.eg.db, keys = ensg, columns = "SYMBOL", keytype = "ENSEMBL")
sym  <- sym[!is.na(sym$SYMBOL) & !duplicated(sym$ENSEMBL), ]
mapped <- mat[sym$ENSEMBL, ]
rownames(mapped) <- sym$SYMBOL[match(rownames(mapped), sym$ENSEMBL)]
# collapse duplicated symbols
dup <- unique(rownames(mapped)[duplicated(rownames(mapped))])
if (length(dup) > 0) {
  keep <- !(rownames(mapped) %in% dup)
  agg <- as(rowsum(as(mapped[rownames(mapped) %in% dup, , drop = FALSE], "dgCMatrix"), rownames(mapped)[rownames(mapped) %in% dup]), "dgCMatrix")
  mapped <- rbind(mapped[keep, ], agg)
}
cat("symbol matrix:", dim(mapped), "\n")

## ============ 3. QC + Seurat ============
ngene <- Matrix::colSums(mapped > 0)
ncount <- Matrix::colSums(mapped)
cat("ngene quantiles:", paste(round(quantile(ngene, c(0.1, 0.25, 0.5, 0.9))), collapse = " "),
    "| ncount quantiles:", paste(round(quantile(ncount, c(0.1, 0.25, 0.5, 0.9))), collapse = " "), "\n")
qc <- ngene >= 1000 & ncount >= 20000
cat("cells passing QC:", sum(qc), "of", length(qc), "\n")
so <- CreateSeuratObject(counts = mapped[, qc], project = "ETP")
so$sample <- sample[qc]; so$group <- grp[qc]
so <- NormalizeData(so, scale.factor = 1e6)
so <- CellCycleScoring(so, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
saveRDS(so, file.path(qa, "ETP", "ETP_seurat_QC.rds"))
cat("QC object:", dim(so), "| Phase table:\n"); print(table(so$Phase))
write.csv(so@meta.data[, c("sample", "group", "Phase", "S.Score", "G2M.Score", "nCount_RNA", "nFeature_RNA")],
          file.path(qa, "50_ETP_cell_metadata.csv"))

## ============ 4. inferCNV (refs = 3000 normal T cells from our cohort) ============
ref_df <- fread(file.path(base, "analysis_results/cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices/t_cells_subset_normal_genes_expression.csv"),
                colClasses = list(character = 1, numeric = 2:10555))
ref_genes <- ref_df[[1]]
ref_m <- as.matrix(ref_df[, -1]); rownames(ref_m) <- ref_genes; rm(ref_df); gc()
storage.mode(ref_m) <- "double"
ref_m[!is.finite(ref_m)] <- 0
set.seed(42); ref_cells <- sample(colnames(ref_m), 3000)
ref_m <- ref_m[, ref_cells]
gene_pos_all <- fread(file.path(base, "analysis_results/cell_cycle_analysis/gencode_v19_gene_pos.txt"), header = FALSE)
common <- intersect(intersect(rownames(mapped), rownames(ref_m)), gene_pos_all$V1)
cat("common genes with refs:", length(common), "\n")
cnts <- as.matrix(mapped[common, colnames(so)])
cnts <- cbind(cnts, ref_m[common, ref_cells, drop = FALSE])
ann <- data.frame(cell = c(rep("Malignant", ncol(so)), rep("Normal", length(ref_cells))),
                  row.names = colnames(cnts))
cat("cnts:", dim(cnts), "\n")
gene_order_path <- file.path(base, "analysis_results/cell_cycle_analysis/gencode_v19_gene_pos.txt")
obj <- CreateInfercnvObject(raw_counts_matrix = cnts, annotations_file = ann,
                            gene_order_file = gene_order_path, ref_group_names = "Normal")
obj <- infercnv::run(obj, num_threads = 24, cutoff = 0.1,
                     out_dir = file.path(qa, "ETP", "infercnv_out"),
                     cluster_by_groups = TRUE, plot_steps = FALSE, HMM = TRUE, denoise = TRUE)
saveRDS(obj, file.path(qa, "ETP", "ETP_infercnv_obj.rds"))
cnv_vals <- obj@expr.data
per_chr_score <- apply(cnv_vals, 2, function(x) mean((x - 1)^2))
counts_chr <- per_chr_score * 0
# approximate per-cell count of CNV-bearing chromosomes via SD across genes per chr is heavy;
# use overall deviation-based burden instead (documented):
burden <- per_chr_score
so$InferCNV_burden <- as.numeric(burden[colnames(so)])
write.csv(data.frame(cell = colnames(so), burden = so$InferCNV_burden), file.path(qa, "51_ETP_burden.csv"), row.names = FALSE)
cat("burden quantiles:", paste(round(quantile(so$InferCNV_burden, c(0.25, 0.5, 0.75, 0.9)), 4), collapse = " "), "\n")

## ============ 5. Gradients & correlations ============
so$bin <- cut(so$InferCNV_burden, breaks = quantile(so$InferCNV_burden, c(0, 0.33, 0.66, 1)),
              labels = c("low", "mid", "high"), include.lowest = TRUE)
tb <- table(so$Phase, so$bin)
write.csv(as.data.frame.matrix(tb), file.path(qa, "52_ETP_phase_by_bin.csv"))
print(tb)
genes_ck <- c("LGALS9", "MIF", "TYMS", "MKI67", "CDKN2A", "HAVCR2", "CD74")
gexp <- GetAssayData(so, layer = "data")
G <- do.call(rbind, lapply(intersect(genes_ck, rownames(gexp)), function(g) {
  ct <- suppressWarnings(cor.test(so$InferCNV_burden, as.numeric(gexp[g, ]), method = "spearman"))
  data.frame(gene = g, rho = unname(ct$estimate), p = ct$p.value,
             mean_low = mean(gexp[g, so$bin == "low"]), mean_mid = mean(gexp[g, so$bin == "mid"]),
             mean_high = mean(gexp[g, so$bin == "high"]))
}))
write.csv(G, file.path(qa, "53_ETP_burden_gene_correlation.csv"), row.names = FALSE)
print(G)

## ============ 6. GSE168594 quick pass (human cells) ============
h5s <- list.files(file.path(base, "GSE168594_RAW"), pattern = "\\.h5$", full.names = TRUE)
pd <- do.call(rbind, lapply(h5s, function(f) {
  m <- Read10X_h5(f)
  hum <- grepl("^hg19", rownames(m)); mou <- grepl("^mm10", rownames(m))
  hum_umi <- Matrix::colSums(m[hum, , drop = FALSE]); tot_umi <- Matrix::colSums(m)
  is_hum <- hum_umi / tot_umi > 0.5
  hum_cells <- m[hum, is_hum]; rownames(hum_cells) <- sub("^hg19_", "", rownames(hum_cells))
  sm <- colSums(hum_cells)
  data.frame(sample = sub("_filtered_feature_bc_matrix.h5", "", basename(f)),
             n_cells_total = ncol(m), n_human = sum(is_hum),
             mean_LGALS9 = if ("LGALS9" %in% rownames(hum_cells)) mean(hum_cells["LGALS9", ] / sm * 1e4) else NA,
             mean_MIF = if ("MIF" %in% rownames(hum_cells)) mean(hum_cells["MIF", ] / sm * 1e4) else NA,
             mean_HAVCR2 = if ("HAVCR2" %in% rownames(hum_cells)) mean(hum_cells["HAVCR2", ] / sm * 1e4) else NA)
}))
write.csv(pd, file.path(qa, "54_PDX_h5_human_summary.csv"), row.names = FALSE)
print(pd)

sink(type = "message"); sink(type = "output")
cat("DONE in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
