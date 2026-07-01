#!/usr/bin/env Rscript
# Presence/absence-ONLY strain-specific counts (the component of Fig_strain_specific_DA16 WITHOUT the edgeR DA
# filter), pipeline-consistent: same filterByExpr floor + same presence rule as edger16.R, but skip glmQLFit.
# A piRNA is presence/absence-only strain-specific for X if present (>=1 read in >=2/3 reps) in X AND in 0 of the
# other 15 strains. Run per timepoint. Output: edger16/<tp>.presence_only_counts.csv
suppressMessages({library(edgeR); library(data.table)})
tp  <- commandArgs(TRUE)[1]
DIR <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
pre <- file.path(DIR, tp)
cnt  <- as.matrix(fread(cmd=paste("zcat", paste0(pre,".counts.tsv.gz"))))
samp <- fread(paste0(pre,".samples.tsv"))
stopifnot(ncol(cnt)==nrow(samp))
strain <- factor(samp$strain)
y    <- DGEList(counts=cnt, group=strain, lib.size=samp$libsize_window)
keep <- filterByExpr(y, group=strain)             # identical data-driven floor to edger16.R
cnt_k <- cnt[keep,,drop=FALSE]
det  <- cnt_k >= 1
sidx <- split(seq_len(ncol(cnt_k)), strain)
present <- sapply(levels(strain), function(s) rowSums(det[, sidx[[s]], drop=FALSE]) >= 2)   # present in >=2/3 reps
strains <- levels(strain)
po <- sapply(strains, function(X){ others <- setdiff(strains, X)
            sum(present[,X] & rowSums(present[,others,drop=FALSE])==0) })                   # present in X, absent in all others
res <- data.frame(strain=strains, timepoint=tp, presence_only=as.integer(po))
fwrite(res, paste0(pre,".presence_only_counts.csv"))
cat(sprintf("[%s] features kept=%d ; presence/absence-only strain-specific per strain:\n", tp, sum(keep)))
print(res)
