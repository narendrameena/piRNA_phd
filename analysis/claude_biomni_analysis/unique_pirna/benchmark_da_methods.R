#!/usr/bin/env Rscript
# benchmark_da_methods.R — DATA-BASED benchmark: edgeR vs DESeq2 for strain-specific piRNA DA.
# Design: one-vs-rest per strain (16 strains x 3 reps). Benchmark design BioMNI 3/3-verified 2026-06-23
# (Soneson & Delorenzi 2013; Schurch 2016; Love 2014; Robinson 2010). Metrics:
#   (1) CONCORDANCE of real-data calls (overlap/Jaccard per strain + overall).
#   (2) PERMUTATION NULL false-positive control: shuffle strain labels (true structure destroyed) ->
#       count "strain-specific" calls; FEWER = better FDR control. Same permutations fed to both methods.
#   (3) P-VALUE CALIBRATION under null: fraction of tests with p<0.05 (target 0.05) + KS vs uniform.
# Both methods use the SAME filterByExpr testable universe; each its own standard normalization
# (edgeR TMM on libsize_window; DESeq2 median-of-ratios). Sig = (FDR/padj<0.05 & logFC>0), one-vs-rest.
# Usage: benchmark_da_methods.R <tp> <n_perm> <seed> <max_features(0=all)>
suppressMessages({library(edgeR); library(DESeq2); library(data.table)})
a <- commandArgs(TRUE); tp <- a[1]; NPERM <- as.integer(a[2])
SEED <- if (length(a)>=3 && !is.na(a[3])) as.integer(a[3]) else 1
MAXF <- if (length(a)>=4 && !is.na(a[4])) as.integer(a[4]) else 0
DIR <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
OUT <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/bench_da"
dir.create(OUT, showWarnings=FALSE); pre <- file.path(DIR, tp)
cnt  <- as.matrix(fread(cmd=paste("zcat", paste0(pre,".counts.tsv.gz"))))
samp <- fread(paste0(pre,".samples.tsv"))
seqs <- fread(cmd=paste("zcat", paste0(pre,".seqs.txt.gz")), header=FALSE)$V1
stopifnot(nrow(cnt)==length(seqs), ncol(cnt)==nrow(samp))
strain <- factor(samp$strain); strains <- levels(strain)
y0 <- DGEList(counts=cnt, group=strain, lib.size=samp$libsize_window)
keep <- filterByExpr(y0, group=strain)
cat(sprintf("[%s] %d features -> %d after filterByExpr\n", tp, nrow(cnt), sum(keep))); flush.console()
cnt <- cnt[keep,,drop=FALSE]; seqs <- seqs[keep]
if (MAXF>0 && nrow(cnt)>MAXF){ s <- sort(sample(nrow(cnt), MAXF)); cnt <- cnt[s,,drop=FALSE]; seqs <- seqs[s]
  cat(sprintf("  [smoke] subsampled to %d features\n", MAXF)) }
storage.mode(cnt) <- "integer"; rownames(cnt) <- seqs; colnames(cnt) <- samp$sample

run_edger <- function(cntm, fac){
  y <- DGEList(counts=cntm, group=fac); y <- calcNormFactors(y)
  des <- model.matrix(~0+fac); colnames(des) <- levels(fac)
  y <- estimateDisp(y, des); fit <- glmQLFit(y, des); lv <- levels(fac)
  pv <- c(); sig <- list()
  for (X in lv){ con <- setNames(rep(-1/(length(lv)-1), length(lv)), lv); con[X] <- 1
    tt <- topTags(glmQLFTest(fit, contrast=con[colnames(des)]), n=Inf, sort.by="none")$table
    pv <- c(pv, tt$PValue); sig[[X]] <- rownames(cntm)[tt$FDR < 0.05 & tt$logFC > 0] }
  list(pv=pv, sig=sig)
}
run_deseq <- function(cntm, fac){
  dds <- DESeqDataSetFromMatrix(cntm, data.frame(strain=fac, row.names=colnames(cntm)), ~0+strain)
  dds <- DESeq(dds, quiet=TRUE); rn <- resultsNames(dds); lv <- levels(fac)
  pv <- c(); sig <- list()
  for (X in lv){ con <- setNames(rep(-1/(length(lv)-1), length(rn)), rn); con[paste0("strain",X)] <- 1
    res <- results(dds, contrast=con)
    pv <- c(pv, res$pvalue); sig[[X]] <- rownames(cntm)[which(res$padj < 0.05 & res$log2FoldChange > 0)] }
  list(pv=pv, sig=sig)
}
cat("== REAL-data fits ==\n"); flush.console()
re <- run_edger(cnt, strain); rd <- run_deseq(cnt, strain)
conc <- do.call(rbind, lapply(strains, function(X){ e<-re$sig[[X]]; d<-rd$sig[[X]]
  i<-length(intersect(e,d)); u<-length(union(e,d))
  data.frame(strain=X, edgeR=length(e), DESeq2=length(d), shared=i, jaccard=ifelse(u>0,i/u,NA)) }))
print(conc)
allE <- unique(unlist(re$sig)); allD <- unique(unlist(rd$sig))
ov <- length(intersect(allE,allD)); un <- length(union(allE,allD))
cat(sprintf("OVERALL real: edgeR=%d DESeq2=%d shared=%d jaccard=%.3f (edgeR-only=%d DESeq2-only=%d)\n",
            length(allE), length(allD), ov, ov/un, length(setdiff(allE,allD)), length(setdiff(allD,allE)))); flush.console()
cat(sprintf("== PERMUTATION NULL (%d perms, seed %d) ==\n", NPERM, SEED)); set.seed(SEED)
perm <- do.call(rbind, lapply(1:NPERM, function(p){
  fac <- factor(sample(as.character(strain)))
  pe <- run_edger(cnt, fac); pd <- run_deseq(cnt, fac)
  r <- data.frame(perm=p, edgeR_FP=sum(lengths(pe$sig)), DESeq2_FP=sum(lengths(pd$sig)),
                  edgeR_p05=mean(pe$pv<0.05, na.rm=TRUE), DESeq2_p05=mean(pd$pv<0.05, na.rm=TRUE),
                  edgeR_KS=ks.test(pe$pv,"punif")$statistic, DESeq2_KS=ks.test(pd$pv,"punif")$statistic)
  cat(sprintf("perm %d: FP edgeR=%d DESeq2=%d | p<.05 edgeR=%.4f DESeq2=%.4f | KS edgeR=%.3f DESeq2=%.3f\n",
              p, r$edgeR_FP, r$DESeq2_FP, r$edgeR_p05, r$DESeq2_p05, r$edgeR_KS, r$DESeq2_KS)); flush.console()
  if (p==1) saveRDS(list(edgeR=pe$pv, DESeq2=pd$pv), file.path(OUT, paste0(tp,".nullpv.rds")))
  r }))
fwrite(conc, file.path(OUT, paste0(tp,".concordance.csv")))
fwrite(perm, file.path(OUT, paste0(tp,".permnull.csv")))
cat("\n== SUMMARY (",tp,") ==\n")
cat(sprintf("real overall calls: edgeR=%d  DESeq2=%d  (jaccard %.3f)\n", length(allE), length(allD), ov/un))
cat(sprintf("NULL false-positives (mean): edgeR=%.0f  DESeq2=%.0f   <- fewer is better\n", mean(perm$edgeR_FP), mean(perm$DESeq2_FP)))
cat(sprintf("NULL p<0.05 (mean, target .05): edgeR=%.4f  DESeq2=%.4f\n", mean(perm$edgeR_p05), mean(perm$DESeq2_p05)))
cat(sprintf("NULL KS-to-uniform (mean, lower=better): edgeR=%.3f  DESeq2=%.3f\n", mean(perm$edgeR_KS), mean(perm$DESeq2_KS)))
saveRDS(list(conc=conc, perm=perm, allE=allE, allD=allD), file.path(OUT, paste0(tp,".bench.rds")))
cat("done\n")
