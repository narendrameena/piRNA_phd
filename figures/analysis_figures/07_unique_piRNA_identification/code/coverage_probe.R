#!/usr/bin/env Rscript
# How much of TOTAL piRNA EXPRESSION (read mass) is covered by each decomposition set, per strain x timepoint.
# total piRNA expression in strain X = sum of reads in X's 3 libraries over ALL piRNA sequences (full matrix).
# coverage(set) = sum of reads (in X) of the set's piRNAs / total. Sets (all on filterByExpr-kept counts):
#   kept_all        : every kept piRNA (shows how much expression survives the low-count floor)
#   presence_only   : present in X, absent elsewhere — under (>=2 reps, loose absence=<2/3 elsewhere) = current pipeline;
#                     (>=2 reps, strict absence=0 reads anywhere); (>=1 rep, strict absence)  <-- the 1-rep vs 2-rep test
#   intersection    : the saved strain-specific set (DA & presence/absence) = Fig_strain_specific_DA16
# DA-only coverage is NOT here: the per-piRNA edgeR DA flags were not saved (only counts) -> needs a full glmQLFit re-run.
suppressMessages({library(edgeR); library(data.table)})
tp  <- commandArgs(TRUE)[1]
DIR <- "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
pre <- file.path(DIR, tp)
cnt  <- as.matrix(fread(cmd=paste("zcat", paste0(pre,".counts.tsv.gz"))))
samp <- fread(paste0(pre,".samples.tsv")); seqs <- fread(cmd=paste("zcat", paste0(pre,".seqs.txt.gz")), header=FALSE)$V1
strain <- factor(samp$strain); strains <- levels(strain)
y <- DGEList(counts=cnt, group=strain, lib.size=samp$libsize_window)
keep <- filterByExpr(y, group=strain); cnt_k <- cnt[keep,,drop=FALSE]; seqs_k <- seqs[keep]; det <- cnt_k >= 1
sidx_all <- split(seq_len(ncol(cnt)), strain)             # columns per strain in the FULL matrix
sidx_k   <- split(seq_len(ncol(cnt_k)), strain)           # same columns, kept-row matrix
nrep <- sapply(strains, function(s) rowSums(det[, sidx_k[[s]], drop=FALSE])); p1 <- nrep>=1; p2 <- nrep>=2
totreads <- sapply(strains, function(s) rowSums(cnt_k[, sidx_k[[s]], drop=FALSE])); present2reads <- totreads >= 2   # "present elsewhere" = >=2 reads total (1 read = noise guard)
spec <- fread(cmd=paste("zcat", paste0(pre,".strain_specific_DA.csv.gz")))         # intersection (DA ∩ loose absence) = current pipeline
sp2f <- paste0(pre,".strain_specific_DA_2read.csv.gz")
spec2 <- if (file.exists(sp2f)) fread(cmd=paste("zcat", sp2f)) else data.table::data.table(sequence=character(), strain=character())  # DA ∩ >=2-read absence = ADOPTED intersection
rows <- list()
for (X in strains){
  Xcols <- sidx_all[[X]]; o <- setdiff(strains, X)
  Xexpr <- rowSums(cnt_k[, sidx_k[[X]], drop=FALSE])       # per kept-piRNA expression in X (summed over X's reps)
  total_all <- sum(cnt[, Xcols])                           # TOTAL piRNA expression in X (all sequences)
  cov <- function(mask) 100*sum(Xexpr[mask])/total_all
  po_loose2  <- p2[,X] & rowSums(p2[,o,drop=FALSE])==0
  po_strict2 <- p2[,X] & rowSums(p1[,o,drop=FALSE])==0
  po_strict1 <- p1[,X] & rowSums(p1[,o,drop=FALSE])==0
  po_2read2  <- p2[,X] & rowSums(present2reads[,o,drop=FALSE])==0   # focal >=2/3 reps, every other strain <2 reads total (intermediate absence)
  po_2read1  <- p1[,X] & rowSums(present2reads[,o,drop=FALSE])==0   # focal >=1 rep
  is_spec    <- seqs_k %in% spec[strain==X]$sequence
  is_spec2   <- seqs_k %in% spec2[strain==X]$sequence                    # DA ∩ >=2-read (ADOPTED final set)
  is_specS   <- is_spec2 & (rowSums(p1[,o,drop=FALSE])==0)               # DA ∩ strict = >=2-read set ∩ 0 reads elsewhere
  rows[[X]] <- data.frame(strain=X, timepoint=tp, total_reads=total_all,
    cov_kept_all      = round(cov(rep(TRUE,nrow(cnt_k))),4),
    cov_presence_loose2  = round(cov(po_loose2),4),  n_presence_loose2  = sum(po_loose2),
    cov_presence_2read2  = round(cov(po_2read2),4),  n_presence_2read2  = sum(po_2read2),
    cov_presence_2read1  = round(cov(po_2read1),4),  n_presence_2read1  = sum(po_2read1),
    cov_presence_strict2 = round(cov(po_strict2),4), n_presence_strict2 = sum(po_strict2),
    cov_presence_strict1 = round(cov(po_strict1),4), n_presence_strict1 = sum(po_strict1),
    cov_intersection       = round(cov(is_spec),4),  n_intersection       = sum(is_spec),
    cov_intersection_2read = round(cov(is_spec2),4), n_intersection_2read = sum(is_spec2),
    cov_intersection_strict= round(cov(is_specS),4), n_intersection_strict= sum(is_specS))
}
res <- do.call(rbind, rows); fwrite(res, paste0(pre,".coverage_probe.csv"))
# per-REPLICATE coverage (for error bars): each library's strain-private read fraction (cnt_k==cnt columns, only rows filtered)
libr <- list()
for (X in strains){
  o <- setdiff(strains, X)
  po_loose2  <- p2[,X] & rowSums(p2[,o,drop=FALSE])==0
  po_2read2  <- p2[,X] & rowSums(present2reads[,o,drop=FALSE])==0
  po_strict2 <- p2[,X] & rowSums(p1[,o,drop=FALSE])==0
  is_spec    <- seqs_k %in% spec[strain==X]$sequence
  is_spec2   <- seqs_k %in% spec2[strain==X]$sequence
  is_specS   <- is_spec2 & (rowSums(p1[,o,drop=FALSE])==0)
  for (i in seq_along(sidx_all[[X]])){
    c <- sidx_all[[X]][i]; tot <- sum(cnt[, c]); if (tot==0) next
    libr[[length(libr)+1]] <- data.frame(strain=X, timepoint=tp, rep=i,
      cov_loose =100*sum(cnt_k[po_loose2,  c])/tot,
      cov_2read =100*sum(cnt_k[po_2read2,  c])/tot,
      cov_strict=100*sum(cnt_k[po_strict2, c])/tot,
      cov_intersection=100*sum(cnt_k[is_spec, c])/tot,
      cov_intersection_2read=100*sum(cnt_k[is_spec2, c])/tot,
      cov_intersection_strict=100*sum(cnt_k[is_specS, c])/tot, total_reads=tot)
  }
}
fwrite(do.call(rbind, libr), paste0(pre,".coverage_perlib.csv"))
cat(sprintf("[%s] kept=%d\n", tp, sum(keep)))
print(res[,c("strain","cov_kept_all","cov_presence_loose2","cov_presence_strict2","cov_presence_strict1","cov_intersection")])
