#!/usr/bin/env Rscript
# Scale-16: strain-specific piRNA via edgeR (per timepoint, up to 16 strains x 3 reps = 48 samples).
# Identical data-driven logic to the 3-strain pilot edger_da.R (strain-generic; loops over all strains):
#   1. filterByExpr(group=strain): low-count floor from library sizes + min group size (no magic #)
#   2. QL F-test, strain X vs mean(others 15); FDR(BH)<0.05 & logFC>0 = significantly X-enriched
#   3. intersect with presence/absence on FILTERED counts (>=2/3 reps in X, absent in ALL other strains)
#      => strain-SPECIFIC among the 16 strains, not merely enriched.
# lib.size = 24-32 nt (piRNA-window) total per sample. Only DIR differs from the pilot (edger16/).
suppressMessages({library(edgeR); library(data.table)})
tp  <- commandArgs(TRUE)[1]
DIR <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
pre <- file.path(DIR, tp)
cnt  <- as.matrix(fread(cmd=paste("zcat", paste0(pre,".counts.tsv.gz"))))
samp <- fread(paste0(pre,".samples.tsv"))
seqs <- fread(cmd=paste("zcat", paste0(pre,".seqs.txt.gz")), header=FALSE)$V1
stopifnot(nrow(cnt)==length(seqs), ncol(cnt)==nrow(samp))
strain <- factor(samp$strain)

y <- DGEList(counts=cnt, group=strain, lib.size=samp$libsize_window)
keep <- filterByExpr(y, group=strain)
cat(sprintf("[%s] features %d -> %d after filterByExpr (data-driven floor)\n", tp, nrow(y), sum(keep)))
y <- y[keep,, keep.lib.sizes=TRUE]; seqs_k <- seqs[keep]; cnt_k <- cnt[keep,,drop=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0+strain); colnames(design) <- levels(strain)
y <- estimateDisp(y, design); fit <- glmQLFit(y, design)

# presence/absence on filtered counts: detected (>=1 read) in >=2 of a strain's reps
det <- cnt_k >= 1
sidx <- split(seq_len(ncol(cnt_k)), strain)
present <- sapply(levels(strain), function(s) rowSums(det[, sidx[[s]], drop=FALSE]) >= 2)

strains <- levels(strain); res <- list()
for (X in strains){
  others <- setdiff(strains, X)
  con <- setNames(rep(0,length(strains)), strains); con[X] <- 1; con[others] <- -1/length(others)
  qlf <- glmQLFTest(fit, contrast=con[colnames(design)])
  tt  <- topTags(qlf, n=Inf, sort.by="none")$table
  da_sig <- tt$FDR < 0.05 & tt$logFC > 0
  absent_others <- rowSums(present[, others, drop=FALSE]) == 0
  specific <- da_sig & present[, X] & absent_others
  out <- data.frame(sequence=seqs_k[specific], strain=X, timepoint=tp,
                    length=nchar(seqs_k[specific]),
                    logFC=round(tt$logFC[specific],3), FDR=signif(tt$FDR[specific],3))
  cat(sprintf("  %s: DA-enriched=%d | +present/absent => strain-specific=%d\n", X, sum(da_sig), nrow(out)))
  res[[X]] <- out
}
res <- do.call(rbind, res)
fwrite(res, paste0(pre,".strain_specific_DA.csv.gz"))
cat(sprintf("[%s] wrote %s.strain_specific_DA.csv.gz (%d strain-specific piRNAs across %d strains)\n",
            tp, pre, nrow(res), length(strains)))
