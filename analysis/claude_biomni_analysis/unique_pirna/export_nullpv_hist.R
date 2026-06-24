#!/usr/bin/env Rscript
# Export binned null p-value histograms (from benchmark_da_methods.R *.nullpv.rds) for the calibration panel.
# Under a correct null the p-values are ~uniform -> flat histogram; an anti-conservative method spikes near 0.
suppressMessages(library(data.table))
OUT <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/bench_da"
for (tp in c("16.5dpc","12.5dpp","20.5dpp")){
  f <- file.path(OUT, paste0(tp,".nullpv.rds")); if (!file.exists(f)) next
  pv <- readRDS(f); br <- seq(0,1,0.05)
  he <- hist(pv$edgeR[!is.na(pv$edgeR)], breaks=br, plot=FALSE)$counts
  hd <- hist(pv$DESeq2[!is.na(pv$DESeq2)], breaks=br, plot=FALSE)$counts
  dt <- data.table(bin_lo=head(br,-1), bin_hi=tail(br,-1),
                   edgeR_frac=he/sum(he), DESeq2_frac=hd/sum(hd))
  fwrite(dt, file.path(OUT, paste0(tp,".nullpv_hist.csv")))
  cat(sprintf("%s: edgeR n=%d DESeq2 n=%d (frac in [0,0.05]: edgeR=%.3f DESeq2=%.3f; uniform=0.05)\n",
              tp, sum(he), sum(hd), he[1]/sum(he), hd[1]/sum(hd)))
}
