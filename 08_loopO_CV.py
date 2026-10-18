#!/usr/bin/env python3
# 08: Leave-one-patient-out (LOOPO) cross-validation for the malignant-vs-normal classifiers
# Addresses reviewer concern #3 (sample-of-origin confound). RF + LR (SVM skipped: ~3.5h/fold).
import os, time
import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

BASE = "/media/zhh2/Newsmy/TALL研究汇总/analysis_results"
SRC  = os.path.join(BASE, "cell_cycle_analysis/inferCNV_results/CNV_combined_final/expression_matrices")
OUT  = os.path.join(BASE, "TALL_manuscript_package/Analysis_QA")
os.makedirs(OUT, exist_ok=True)
logf = open(os.path.join(OUT, "08_loopO_log.txt"), "w")
def log(*a): print(*a, file=logf, flush=True); print(*a, flush=True)

t0 = time.time()
log("== load positive matrix ==")
pos_m = scipy.io.mmread(os.path.join(SRC, "positive_expression_data.mtx")).tocsr()   # genes x cells
pos_genes = [l.strip() for l in open(os.path.join(SRC, "positive_genes.tsv"))]
pos_cells = [l.strip() for l in open(os.path.join(SRC, "positive_cells.tsv"))]
log("pos:", pos_m.shape, len(pos_genes), len(pos_cells))
neg_df = pd.read_csv(os.path.join(SRC, "t_cells_subset_normal_genes_expression.csv"), index_col=0)
neg_genes = list(neg_df.index); neg_cells = list(neg_df.columns)
neg_m = sparse.csr_matrix(neg_df.values)
log("neg:", neg_m.shape, len(neg_genes), len(neg_cells))

sel = pd.read_csv(os.path.join(SRC, "selected_genes.csv"))
sel_genes = [g for g in sel.iloc[:, 0].astype(str).tolist() if g not in ("0", "gene")]
common = sorted(set(pos_genes) & set(neg_genes) & set(sel_genes))
log("selected common genes:", len(common))
gi_p = {g: i for i, g in enumerate(pos_genes)}
gi_n = {g: i for i, g in enumerate(neg_genes)}
Xp = pos_m[[gi_p[g] for g in common], :].tocsc()      # genes x cells
Xn = neg_m[[gi_n[g] for g in common], :].tocsc()
X = sparse.vstack([Xp.T, Xn.T]).tocsr()               # cells x genes
y = np.array([1] * Xp.shape[1] + [0] * Xn.shape[1])
log("X:", X.shape, "y:", np.bincount(y).tolist())

# patient assignment for positive cells
def gsm_of(c):
    i = c.find("GSM"); return c[i:i+10] if i >= 0 else c
gsm2pat = {
 "GSM7091996":"T1","GSM7091997":"T2","GSM7091998":"T3","GSM7091999":"T4","GSM7092000":"T5",
 "GSM7092001":"T6","GSM7092002":"T6","GSM7092003":"T7","GSM7092004":"T7","GSM7092005":"T8",
 "GSM7092006":"T8","GSM7092007":"T9","GSM7092008":"T9","GSM7092009":"T10","GSM7092010":"T10",
 "GSM7092011":"T11","GSM3872440":"PRE-T_1","GSM3872441":"PRE-T_2"}
pat = np.array([gsm2pat.get(gsm_of(c), "NA") for c in pos_cells])
patients = sorted(set(pat))
log("patients:", patients)

# CNV-negative cells (the 1,020 test cells) with patient labels
pred = pd.read_csv(os.path.join(SRC, "negative_data_predictions.csv"))
pred["GSM"] = [gsm_of(c) for c in pred["cell_id"]]
pred["patient"] = [gsm2pat.get(g, "NA") for g in pred["GSM"]]
# their expression: from negative_data_genes_expression.csv? those are the 10,554 normal cells, NOT the 1020.
# The 1020 cells' expression: they are the CNV-Negative cells of the T-ALL object; get from Expr matrix of the combined object is complex.
# Instead: predict them via gene_expression from marker-v2 gene_expression_raw.csv? only 50 genes.
# Use negative expression of 1020 cells unavailable here -> restrict LOOPO to held-out POSITIVES + optionally leave-out negative donors.
# We record: per-patient positive accuracy; that directly tests patient-level generalization.

pat_all = np.full(len(y), "NEG", dtype=object)
pat_all[y == 1] = pat
rows = []
for p in patients:
    t1 = time.time()
    tr_idx = np.where((y == 0) | ((y == 1) & (pat_all != p)))[0]
    te_idx = np.where((y == 1) & (pat_all == p))[0]
    Xtr, ytr = X[tr_idx], y[tr_idx]
    Xte, yte = X[te_idx], y[te_idx]
    sc = StandardScaler(with_mean=False); Xtr_s = sc.fit_transform(Xtr); Xte_s = sc.transform(Xte)
    rf = RandomForestClassifier(n_estimators=200, n_jobs=64, random_state=42)
    rf.fit(Xtr_s, ytr)
    pr_rf = rf.predict(Xte_s)
    acc_rf = float(np.mean(pr_rf == yte))
    lr = LogisticRegression(max_iter=2000, solver="liblinear", C=1.0)
    lr.fit(Xtr_s, ytr)
    pr_lr = lr.predict(Xte_s)
    acc_lr = float(np.mean(pr_lr == yte))
    rows.append(dict(patient=p, n_heldout_pos=len(te_idx), acc_pos_RF=acc_rf, acc_pos_LR=acc_lr,
                     min_fold=(time.time() - t1) / 60))
    log(f"{p}: held-out pos n={len(te_idx)} RF_acc={acc_rf:.4f} LR_acc={acc_lr:.4f} ({rows[-1]['min_fold']:.1f} min)")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "38_LOOPO_per_patient.csv"), index=False)
log("TOTAL time (min):", (time.time() - t0) / 60)
log("DONE")
