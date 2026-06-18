#!/usr/bin/env Rscript
# DATA-DRIVEN PROBE for the presence/absence threshold question (does the definition change the call?).
# Per timepoint, on filterByExpr-kept counts, count strain-specific piRNAs per strain under 3 definitions:
#   current        : present(X) >=2/3 reps  AND  every other strain present <2/3 reps   (the pipeline; "absent" allows up to 1/3 reps elsewhere)
#   strict_abs     : present(X) >=2/3 reps  AND  ZERO reads in ANY rep of ANY other strain  (biologically strict absence)
#   rep1_strictabs : present(X) >=1 rep     AND  ZERO reads in ANY rep of ANY other strain  ("at least one rep" present + truly absent)
suppressMessages({library(edgeR); library(data.table)})
tp  <- commandArgs(TRUE)[1]
DIR <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
pre <- file.path(DIR, tp)
cnt  <- as.matrix(fread(cmd=paste("zcat", paste0(pre,".counts.tsv.gz"))))
samp <- fread(paste0(pre,".samples.tsv")); strain <- factor(samp$strain); strains <- levels(strain)
y    <- DGEList(counts=cnt, group=strain, lib.size=samp$libsize_window)
keep <- filterByExpr(y, group=strain); cnt_k <- cnt[keep,,drop=FALSE]; det <- cnt_k >= 1
sidx <- split(seq_len(ncol(cnt_k)), strain)
nrep <- sapply(strains, function(s) rowSums(det[, sidx[[s]], drop=FALSE]))   # reps detecting (0..3) per strain
p2 <- nrep >= 2; p1 <- nrep >= 1
res <- data.frame(strain=strains, timepoint=tp,
  current        = sapply(strains, function(X){ o<-setdiff(strains,X); sum(p2[,X] & rowSums(p2[,o,drop=FALSE])==0) }),
  strict_abs     = sapply(strains, function(X){ o<-setdiff(strains,X); sum(p2[,X] & rowSums(p1[,o,drop=FALSE])==0) }),
  rep1_strictabs = sapply(strains, function(X){ o<-setdiff(strains,X); sum(p1[,X] & rowSums(p1[,o,drop=FALSE])==0) }))
fwrite(res, paste0(pre,".threshold_probe.csv"))
cat(sprintf("[%s] kept=%d\n", tp, sum(keep))); print(res)
cat(sprintf("[%s] TOTALS  current=%d  strict_abs=%d (%.0f%% of current)  rep1_strictabs=%d\n",
            tp, sum(res$current), sum(res$strict_abs), 100*sum(res$strict_abs)/sum(res$current), sum(res$rep1_strictabs)))
