#!/usr/bin/env Rscript
# deseq_filter_order_test.R — EFFECT of WHERE the 27/30 nt stage-peak length filter is applied vs DESeq2.
#   Order A (BEFORE = production, user 2026-06-23): raw -> length-filter(tp) -> filterByExpr -> DESeq2
#       one-vs-rest -> +>=2-read presence/absence -> candidates_A
#   Order B (AFTER): raw -> filterByExpr(full 25-32) -> DESeq2 one-vs-rest -> +>=2-read -> SUBSET to tp-lengths
#       -> candidates_B
# The two differ ONLY in size factors, dispersion trend, and the BH denominator (computed on the 27/30 subset
# vs the full 25-32 universe). >=2-read presence/absence is per-sequence (order-independent). Compare A vs B.
# tp-lengths: 16.5dpc=27 ; 12.5dpp=27&30 ; 20.5dpp=30.
# Usage: deseq_filter_order_test.R <tp> [orderB_maxf(0=full, >0 = smoke cap)]
suppressMessages({library(DESeq2); library(edgeR); library(data.table)})
a <- commandArgs(TRUE); tp <- a[1]; BCAP <- if (length(a)>=2 && !is.na(a[2])) as.integer(a[2]) else 0
KEEP <- list("16.5dpc"=c(27),"12.5dpp"=c(27,30),"20.5dpp"=c(30))[[tp]]
DIR <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
OUT <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/deseq16_lenfilt"
dir.create(OUT, showWarnings=FALSE); pre <- file.path(DIR, tp)
C <- as.matrix(fread(cmd=paste("zcat", paste0(pre,".counts.tsv.gz"))))
samp <- fread(paste0(pre,".samples.tsv")); seqs <- fread(cmd=paste("zcat", paste0(pre,".seqs.txt.gz")), header=FALSE)$V1
stopifnot(nrow(C)==length(seqs)); storage.mode(C) <- "integer"
strain <- factor(samp$strain); strains <- levels(strain); n <- length(strains)
call_specific <- function(cnt, sq){            # DESeq2 one-vs-rest + >=2-read presence/absence
  rownames(cnt) <- sq; colnames(cnt) <- samp$sample
  dds <- DESeq(DESeqDataSetFromMatrix(cnt, data.frame(strain=strain, row.names=samp$sample), ~0+strain), quiet=TRUE)
  rn <- resultsNames(dds); sidx <- split(seq_len(ncol(cnt)), strain); det <- cnt>=1
  present <- sapply(strains, function(s) rowSums(det[,sidx[[s]],drop=FALSE])>=2)          # focal >=2/3 reps
  pres2   <- sapply(strains, function(s) rowSums(cnt[,sidx[[s]],drop=FALSE]))>=2           # strain present = >=2 reads
  out <- list()
  for (X in strains){
    con <- setNames(rep(-1/(n-1), length(rn)), rn); con[paste0("strain",X)] <- 1
    r <- results(dds, contrast=con)
    da <- !is.na(r$padj) & r$padj<0.05 & r$log2FoldChange>0
    ab <- rowSums(pres2[,setdiff(strains,X),drop=FALSE])==0
    k <- which(da & present[,X] & ab)
    out[[X]] <- data.frame(sequence=sq[k], strain=rep(X,length(k)), length=nchar(sq[k]),
                           log2FC=round(r$log2FoldChange[k],3), padj=signif(r$padj[k],3))
  }
  do.call(rbind, out)
}
# ---- Order A: length filter BEFORE DESeq2 ----
lf <- nchar(seqs) %in% KEEP; Ca <- C[lf,,drop=FALSE]; sa <- seqs[lf]
keepA <- filterByExpr(DGEList(Ca, group=strain, lib.size=samp$libsize_window), group=strain)
cat(sprintf("[A %s] length-filter %s nt: %d -> filterByExpr %d\n", tp, paste(KEEP,collapse="&"), sum(lf), sum(keepA))); flush.console()
candA <- call_specific(Ca[keepA,,drop=FALSE], sa[keepA])
cat(sprintf("[A %s] candidates_A = %d\n", tp, nrow(candA))); flush.console()
# ---- Order B: DESeq2 on FULL 25-32, length filter AFTER ----
keepB <- filterByExpr(DGEList(C, group=strain, lib.size=samp$libsize_window), group=strain)
idxB <- which(keepB)
if (BCAP>0 && length(idxB)>BCAP){ idxB <- sort(union(idxB[nchar(seqs[idxB]) %in% KEEP][1:min(BCAP,sum(nchar(seqs[idxB])%in%KEEP))], sample(idxB, BCAP)))
  cat(sprintf("[B %s] SMOKE cap -> %d features\n", tp, length(idxB))) }
cat(sprintf("[B %s] full filterByExpr %d (DESeq2 on full set)\n", tp, length(idxB))); flush.console()
candB_all <- call_specific(C[idxB,,drop=FALSE], seqs[idxB])
candB <- candB_all[nchar(candB_all$sequence) %in% KEEP,]
# ---- compare ----
idA <- paste(candA$strain, candA$sequence); idB <- paste(candB$strain, candB$sequence)
sh <- length(intersect(idA,idB)); un <- length(union(idA,idB))
cat(sprintf("\n== FILTER-ORDER EFFECT (%s, lengths %s) ==\n", tp, paste(KEEP,collapse="&")))
cat(sprintf("Order A (filter BEFORE): %d stage-peak candidates\n", length(idA)))
cat(sprintf("Order B (filter AFTER):  %d stage-peak candidates  (from %d full-set strain-specific)\n", length(idB), nrow(candB_all)))
cat(sprintf("shared=%d  A-only=%d  B-only=%d  Jaccard=%.3f\n", sh, length(setdiff(idA,idB)), length(setdiff(idB,idA)), sh/un))
fwrite(candA, file.path(OUT, paste0(tp,".candidates_orderA_before.csv.gz")))
fwrite(candB, file.path(OUT, paste0(tp,".candidates_orderB_after.csv.gz")))
perstrain <- data.frame(strain=strains,
  A=as.integer(table(factor(candA$strain, strains))), B=as.integer(table(factor(candB$strain, strains))))
fwrite(perstrain, file.path(OUT, paste0(tp,".filterorder_perstrain.csv")))
saveRDS(list(tp=tp, A=length(idA), B=length(idB), shared=sh, jaccard=sh/un, nB_all=nrow(candB_all),
             keepA=sum(keepA), keepB=length(idxB)), file.path(OUT, paste0(tp,".filterorder.rds")))
cat("done\n")
