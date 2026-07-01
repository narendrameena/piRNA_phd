#!/usr/bin/env Rscript
# ACCURATE expression-coverage decomposition (per timepoint), using the SAME edgeR method as edger16.R
# (filterByExpr -> calcNormFactors -> tagwise estimateDisp -> glmQLFit -> per-strain glmQLFTest; NO shortcuts).
# For each strain X computes the EXPRESSION COVERAGE = (reads in X of the set) / (total piRNA reads in X), and counts:
#   DA_only        : edgeR FDR<0.05 & logFC>0  (strain X vs mean of other 15) — strain-ENRICHED, not specific
#   presence_loose2: present >=2/3 reps in X AND every other strain present <2/3 reps  (= current pipeline's presence/absence)
#   presence_strict2: present >=2/3 reps in X AND 0 reads in ANY rep of ANY other strain (strict absence)
#   presence_strict1: present >=1 rep    in X AND 0 reads in ANY rep of ANY other strain (the 1-rep vs 2-rep test)
#   intersection   : DA_only & presence_loose2  (= Fig_strain_specific_DA16 final set; cross-checked vs saved csv)
# total piRNA expression in X = sum of reads in X's 3 libraries over ALL sequences (full matrix). Output: <tp>.coverage_full.csv
suppressMessages({library(edgeR); library(data.table)})
tp  <- commandArgs(TRUE)[1]
DIR <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
pre <- file.path(DIR, tp)
cnt  <- as.matrix(fread(cmd=paste("zcat", paste0(pre,".counts.tsv.gz"))))
samp <- fread(paste0(pre,".samples.tsv"))
seqs <- fread(cmd=paste("zcat", paste0(pre,".seqs.txt.gz")), header=FALSE)$V1
stopifnot(nrow(cnt)==length(seqs), ncol(cnt)==nrow(samp))
strain <- factor(samp$strain); strains <- levels(strain)
y <- DGEList(counts=cnt, group=strain, lib.size=samp$libsize_window)
keep <- filterByExpr(y, group=strain)
cat(sprintf("[%s] %d -> %d after filterByExpr\n", tp, nrow(y), sum(keep))); flush.console()
y <- y[keep,, keep.lib.sizes=TRUE]; seqs_k <- seqs[keep]; cnt_k <- cnt[keep,,drop=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0+strain); colnames(design) <- levels(strain)
y <- estimateDisp(y, design); fit <- glmQLFit(y, design)              # the ~2 h step (tagwise dispersion)
cat(sprintf("[%s] dispersion+fit done\n", tp)); flush.console()
det <- cnt_k >= 1
sidx_k   <- split(seq_len(ncol(cnt_k)), strain)
sidx_all <- split(seq_len(ncol(cnt)),  strain)
nrep <- sapply(strains, function(s) rowSums(det[, sidx_k[[s]], drop=FALSE])); p1 <- nrep>=1; p2 <- nrep>=2
spec <- fread(cmd=paste("zcat", paste0(pre,".strain_specific_DA.csv.gz")))   # saved intersection set, for consistency check
rows <- list()
for (X in strains){
  o <- setdiff(strains, X)
  con <- setNames(rep(0, length(strains)), strains); con[X] <- 1; con[o] <- -1/length(o)
  qlf <- glmQLFTest(fit, contrast=con[colnames(design)]); tt <- topTags(qlf, n=Inf, sort.by="none")$table
  da  <- tt$FDR < 0.05 & tt$logFC > 0
  po_loose2  <- p2[,X] & rowSums(p2[,o,drop=FALSE])==0
  po_strict2 <- p2[,X] & rowSums(p1[,o,drop=FALSE])==0
  po_strict1 <- p1[,X] & rowSums(p1[,o,drop=FALSE])==0
  inter      <- da & po_loose2
  Xexpr <- rowSums(cnt_k[, sidx_k[[X]], drop=FALSE]); total_all <- sum(cnt[, sidx_all[[X]]])
  cov <- function(m) round(100*sum(Xexpr[m])/total_all, 4)
  rows[[X]] <- data.frame(strain=X, timepoint=tp, total_reads=total_all,
    n_DA=sum(da),                 cov_DA=cov(da),
    n_presence_loose2=sum(po_loose2),   cov_presence_loose2=cov(po_loose2),
    n_presence_strict2=sum(po_strict2), cov_presence_strict2=cov(po_strict2),
    n_presence_strict1=sum(po_strict1), cov_presence_strict1=cov(po_strict1),
    n_intersection=sum(inter),    cov_intersection=cov(inter),
    n_intersection_saved=sum(seqs_k %in% spec[strain==X]$sequence))
}
res <- do.call(rbind, rows); fwrite(res, paste0(pre,".coverage_full.csv"))
cat(sprintf("[%s] DONE. intersection recomputed==saved? %s\n", tp, all(res$n_intersection==res$n_intersection_saved)))
print(res[,c("strain","cov_DA","cov_presence_loose2","cov_presence_strict2","cov_presence_strict1","cov_intersection")])
