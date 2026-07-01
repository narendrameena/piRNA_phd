#!/usr/bin/env Rscript
# Combined all-timepoints-pooled PCA, computed SEPARATELY for each feature-set class (companion bottom row to the
# 9 per-timepoint panels in Fig_pca_unique16). Pool E16.5+P12.5+P20.5 on shared sequences, DESeq2 size-factor
# normalise (log2), then PCA for: all_full (every shared feature), all_expressed (top-500 most variable),
# unique (genuinely-unique = klass5 startswith "unique", ≥2-read adopted classification). Gram-matrix PCA
# (X Xᵀ, 144×144) avoids the huge rotation. Output: pca16/combined_byclass.csv (+ feature_set column).
suppressMessages({library(DESeq2); library(data.table)})
DIR<-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
U  <-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
OUT<-file.path(U,"pca16"); tps<-c("16.5dpc","12.5dpp","20.5dpp"); labs<-c("E16.5","P12.5","P20.5")
mats<-list(); samps<-list()
for(i in seq_along(tps)){
  tp<-tps[i]
  cnt<-as.matrix(fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".counts.tsv.gz")))))
  samp<-fread(file.path(DIR,paste0(tp,".samples.tsv")))
  seqs<-fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".seqs.txt.gz"))),header=FALSE)$V1
  rownames(cnt)<-seqs; colnames(cnt)<-paste0(samp$sample,"|",labs[i])
  samp$tp<-labs[i]; samp$col<-colnames(cnt); mats[[i]]<-cnt; samps[[i]]<-samp
}
common<-Reduce(intersect, lapply(mats, rownames)); cat(length(common),"shared sequences\n")
M<-do.call(cbind, lapply(mats, function(m) m[common,,drop=FALSE])); SAMP<-rbindlist(samps, fill=TRUE)
cd<-DataFrame(strain=factor(SAMP$strain), tp=factor(SAMP$tp), row.names=SAMP$col)
dds<-DESeqDataSetFromMatrix(M, cd, ~1); dds<-estimateSizeFactors(dds)
logn<-log2(counts(dds,normalized=TRUE)+1)
# genuinely-unique sequences (ADOPTED >=2-read classification)
uq<-fread(cmd=paste("zcat",file.path(U,"unique16/final_classified_clean_2read.csv.gz")))
uqseq<-unique(uq[grepl("^unique",klass5)]$sequence)
nC<-ncol(logn); v<-(rowSums(logn^2)-rowSums(logn)^2/nC)/(nC-1); top500<-rownames(logn)[order(v,decreasing=TRUE)][1:min(500,nrow(logn))]
fsets<-list(all_full=rownames(logn), all_expressed=top500, unique=intersect(common,uqseq))
pca<-function(sub){
  m<-scale(t(logn[sub,,drop=FALSE]),center=TRUE,scale=FALSE)
  G<-tcrossprod(m); e<-eigen(G,symmetric=TRUE); ve<-round(100*e$values/sum(e$values),1)
  sc<-e$vectors %*% diag(sqrt(pmax(e$values,0))); list(PC1=sc[,1],PC2=sc[,2],v1=ve[1],v2=ve[2],n=length(sub))
}
res<-list()
for(fs in names(fsets)){ sub<-fsets[[fs]]; if(length(sub)<3){cat("skip",fs,length(sub),"\n");next}
  p<-pca(sub); res[[fs]]<-data.frame(sample=SAMP$sample,strain=SAMP$strain,tp=SAMP$tp,PC1=p$PC1,PC2=p$PC2,
    pc1_var=p$v1,pc2_var=p$v2,n_features=p$n,feature_set=fs)
  cat(sprintf("%s: %d features, PC1 %.1f%% PC2 %.1f%%\n",fs,p$n,p$v1,p$v2)) }
fwrite(rbindlist(res), file.path(OUT,"combined_byclass.csv")); cat("wrote combined_byclass.csv\n")
