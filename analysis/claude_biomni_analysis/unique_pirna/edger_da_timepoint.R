#!/usr/bin/env Rscript
# TIMEPOINT-specific piRNA via edgeR (per strain X, 9 samples = 3 timepoints x 3 reps). Same data-driven
# floor (filterByExpr) + QL F-test (timepoint T vs mean(others), FDR<0.05 & logFC>0) intersected with
# presence/absence on filtered counts (present >=2/3 reps at T, absent <2/3 at each other stage). Within
# one strain/genome => no cross-genome/SNP step needed (a piRNA absent at other stages IS stage-specific).
suppressMessages({library(edgeR); library(data.table)})
X <- commandArgs(TRUE)[1]
DIR<-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger_tp"
pre<-file.path(DIR,X)
cnt <-as.matrix(fread(cmd=paste("zcat",paste0(pre,".counts.tsv.gz"))))
samp<-fread(paste0(pre,".samples.tsv"))
seqs<-fread(cmd=paste("zcat",paste0(pre,".seqs.txt.gz")),header=FALSE)$V1
tp<-factor(samp$timepoint,levels=c("16.5dpc","12.5dpp","20.5dpp")); colnames(cnt)<-samp$sample
y<-DGEList(cnt,group=tp,lib.size=samp$libsize_window); keep<-filterByExpr(y,group=tp)
cat(sprintf("[%s] features %d -> %d after filterByExpr\n",X,nrow(y),sum(keep)))
y<-y[keep,,keep.lib.sizes=TRUE]; seqs_k<-seqs[keep]; cnt_k<-cnt[keep,,drop=FALSE]
y<-calcNormFactors(y); design<-model.matrix(~0+tp); colnames(design)<-levels(tp)
y<-estimateDisp(y,design); fit<-glmQLFit(y,design)
det<-cnt_k>=1; tidx<-split(seq_len(ncol(cnt_k)),tp)
present<-sapply(levels(tp),function(t) rowSums(det[,tidx[[t]],drop=FALSE])>=2)
tps<-levels(tp); res<-list()
for(T in tps){
  others<-setdiff(tps,T); con<-setNames(rep(0,length(tps)),tps); con[T]<-1; con[others]<- -1/length(others)
  qlf<-glmQLFTest(fit,contrast=con[colnames(design)]); tt<-topTags(qlf,n=Inf,sort.by="none")$table
  da<-tt$FDR<0.05 & tt$logFC>0; absent_others<-rowSums(present[,others,drop=FALSE])==0
  spec<-da & present[,T] & absent_others
  res[[T]]<-data.frame(sequence=seqs_k[spec],strain=X,timepoint=T,length=nchar(seqs_k[spec]),
                       logFC=round(tt$logFC[spec],3),FDR=signif(tt$FDR[spec],3))
  cat(sprintf("  %s %s: DA-enriched=%d -> timepoint-specific=%d\n",X,T,sum(da),nrow(res[[T]])))
}
res<-do.call(rbind,res); fwrite(res,paste0(pre,".timepoint_specific_DA.csv.gz"))
cat(sprintf("[%s] wrote %d timepoint-specific piRNAs -> %s.timepoint_specific_DA.csv.gz\n",X,nrow(res),pre))
