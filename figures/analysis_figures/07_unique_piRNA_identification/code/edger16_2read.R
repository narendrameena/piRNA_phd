#!/usr/bin/env Rscript
# edger16_2read.R — ADOPTED >=2-read absence definition (BioMNI-signed-off 2026-06-17, 3/3 agents YES).
# Identical to edger16.R EXCEPT the "absent in other strains" rule:
#   OLD (loose):  absent in Y = Y has <2/3 replicates detected (a single replicate, any read count, still = "absent")
#   NEW (>=2read): absent in Y = Y has <2 READS total (one read anywhere = noise/contamination, >=2 reads = present)
# Focal presence unchanged (present in X = detected in >=2 of X's 3 reps). Writes to *_2read.csv.gz (NON-destructive:
# the loose *.strain_specific_DA.csv.gz is left intact). Same edgeR method (filterByExpr -> calcNormFactors ->
# tagwise estimateDisp -> glmQLFit -> per-strain QL F-test, FDR<0.05 & logFC>0).
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
cat(sprintf("[%s] features %d -> %d after filterByExpr\n", tp, nrow(y), sum(keep))); flush.console()
y <- y[keep,, keep.lib.sizes=TRUE]; seqs_k <- seqs[keep]; cnt_k <- cnt[keep,,drop=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~0+strain); colnames(design) <- levels(strain)
y <- estimateDisp(y, design); fit <- glmQLFit(y, design)
cat(sprintf("[%s] dispersion+fit done\n", tp)); flush.console()
det  <- cnt_k >= 1
sidx <- split(seq_len(ncol(cnt_k)), strain)
present  <- sapply(levels(strain), function(s) rowSums(det[, sidx[[s]], drop=FALSE]) >= 2)        # focal presence: >=2/3 reps
totreads <- sapply(levels(strain), function(s) rowSums(cnt_k[, sidx[[s]], drop=FALSE]))
present2 <- totreads >= 2                                                                          # present-elsewhere: >=2 reads total
strains <- levels(strain); res <- list()
for (X in strains){
  others <- setdiff(strains, X)
  con <- setNames(rep(0,length(strains)), strains); con[X] <- 1; con[others] <- -1/length(others)
  qlf <- glmQLFTest(fit, contrast=con[colnames(design)]); tt <- topTags(qlf, n=Inf, sort.by="none")$table
  da_sig <- tt$FDR < 0.05 & tt$logFC > 0
  absent_others <- rowSums(present2[, others, drop=FALSE]) == 0                                    # >=2-read absence
  specific <- da_sig & present[, X] & absent_others
  spq <- seqs_k[specific]
  out <- data.frame(sequence=spq, strain=rep(X,length(spq)), timepoint=rep(tp,length(spq)),
                    length=nchar(spq), logFC=round(tt$logFC[specific],3), FDR=signif(tt$FDR[specific],3))
  cat(sprintf("  %s: DA-enriched=%d | + >=2-read presence/absence => strain-specific=%d\n", X, sum(da_sig), nrow(out)))
  res[[X]] <- out
}
res <- do.call(rbind, res)
fwrite(res, paste0(pre,".strain_specific_DA_2read.csv.gz"))
cat(sprintf("[%s] wrote %s.strain_specific_DA_2read.csv.gz (%d strain-specific piRNAs)\n", tp, pre, nrow(res)))
