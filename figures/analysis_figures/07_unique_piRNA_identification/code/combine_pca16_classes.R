#!/usr/bin/env Rscript
# PCA of EACH of the 5 classification classes (klass5, ≥2-read adopted) SEPARATELY — per-timepoint and combined.
# For class k: subset the DESeq2 size-factor-normalised (log2) expression to the sequences labelled k, Gram-matrix
# PCA (X Xᵀ). Per-timepoint = that tp's libraries (coloured by strain in the figure); combined = pooled 144 libs on
# shared sequences (coloured by timepoint). Output: pca16/classes_pca.csv (columns + klass5 + view).
suppressMessages({library(DESeq2); library(data.table)})
DIR<-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
U  <-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
OUT<-file.path(U,"pca16"); tps<-c("16.5dpc","12.5dpp","20.5dpp"); labs<-c("E16.5","P12.5","P20.5")
KL5<-c("expressed elsewhere (exact)","SNP-variant (1-3mm)","low-quality: no mm0 own-genome locus","unique: conserved-but-silent","unique: strain-private locus")
fc<-fread(cmd=paste("zcat",file.path(U,"unique16/final_classified_clean_2read.csv.gz")))
gpca<-function(sl){ m<-scale(t(sl),center=TRUE,scale=FALSE); G<-tcrossprod(m); e<-eigen(G,symmetric=TRUE)
  ve<-round(100*e$values/sum(e$values),1); sc<-e$vectors%*%diag(sqrt(pmax(e$values,0))); list(PC1=sc[,1],PC2=sc[,2],v1=ve[1],v2=ve[2]) }
res<-list(); matl<-list(); sampl<-list()
for(i in seq_along(tps)){ tp<-tps[i]
  cnt<-as.matrix(fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".counts.tsv.gz")))))
  samp<-fread(file.path(DIR,paste0(tp,".samples.tsv"))); seqs<-fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".seqs.txt.gz"))),header=FALSE)$V1
  rownames(cnt)<-seqs; colnames(cnt)<-paste0(samp$sample,"|",labs[i])
  dds<-DESeqDataSetFromMatrix(cnt, DataFrame(strain=factor(samp$strain),row.names=colnames(cnt)), ~1); dds<-estimateSizeFactors(dds)
  logn<-log2(counts(dds,normalized=TRUE)+1)
  for(k in KL5){ ks<-intersect(unique(fc[timepoint==tp & klass5==k]$sequence), rownames(logn)); if(length(ks)<3) next
    p<-gpca(logn[ks,,drop=FALSE]); res[[length(res)+1]]<-data.frame(sample=samp$sample,strain=samp$strain,tp=labs[i],
      PC1=p$PC1,PC2=p$PC2,pc1_var=p$v1,pc2_var=p$v2,n_features=length(ks),klass5=k,view="timepoint") }
  matl[[i]]<-cnt; samp$tp<-labs[i]; sampl[[i]]<-samp; cat("tp",tp,"per-class PCA done\n")
}
common<-Reduce(intersect, lapply(matl, rownames)); M<-do.call(cbind, lapply(matl, function(m) m[common,,drop=FALSE])); SAMP<-rbindlist(sampl,fill=TRUE)
ddsC<-DESeqDataSetFromMatrix(M, DataFrame(strain=factor(SAMP$strain),tp=factor(SAMP$tp),row.names=colnames(M)), ~1); ddsC<-estimateSizeFactors(ddsC)
lognC<-log2(counts(ddsC,normalized=TRUE)+1)
for(k in KL5){ ks<-intersect(unique(fc[klass5==k]$sequence), common); if(length(ks)<3) next
  p<-gpca(lognC[ks,,drop=FALSE]); res[[length(res)+1]]<-data.frame(sample=SAMP$sample,strain=SAMP$strain,tp=SAMP$tp,
    PC1=p$PC1,PC2=p$PC2,pc1_var=p$v1,pc2_var=p$v2,n_features=length(ks),klass5=k,view="combined")
  cat(sprintf("combined %s: %d features, PC1 %.1f%%\n",k,length(ks),p$v1)) }
fwrite(rbindlist(res,fill=TRUE), file.path(OUT,"classes_pca.csv")); cat("wrote classes_pca.csv (",length(res),"PCAs )\n")
