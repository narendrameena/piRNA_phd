#!/usr/bin/env Rscript
# Combined all-timepoints expressed-piRNA PCA (companion panel to pca_unique16.R): pool E16.5+P12.5+P20.5 samples
# on shared sequences, DESeq2 size-factor normalise, top-500 variable, prcomp. Shows strain grouping across
# development in one space. Output = pca16/combined.pca.csv.
suppressMessages({library(DESeq2); library(data.table)})
DIR<-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
OUT<-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pca16"
tps<-c("16.5dpc","12.5dpp","20.5dpp"); labs<-c("E16.5","P12.5","P20.5")
mats<-list(); samps<-list()
for(i in seq_along(tps)){
  tp<-tps[i]
  cnt<-as.matrix(fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".counts.tsv.gz")))))
  samp<-fread(file.path(DIR,paste0(tp,".samples.tsv")))
  seqs<-fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".seqs.txt.gz"))),header=FALSE)$V1
  rownames(cnt)<-seqs; colnames(cnt)<-paste0(samp$sample,"|",labs[i])
  samp$tp<-labs[i]; samp$col<-colnames(cnt); mats[[i]]<-cnt; samps[[i]]<-samp
}
common<-Reduce(intersect, lapply(mats, rownames))
cat(length(common),"shared sequences across the 3 timepoints\n")
M<-do.call(cbind, lapply(mats, function(m) m[common,,drop=FALSE]))
SAMP<-rbindlist(samps, fill=TRUE)
cd<-DataFrame(strain=factor(SAMP$strain), tp=factor(SAMP$tp), row.names=SAMP$col)
dds<-DESeqDataSetFromMatrix(M, cd, ~1); dds<-estimateSizeFactors(dds)
logn<-log2(counts(dds,normalized=TRUE)+1)
# ALL shared expressed piRNAs (no top-500 filter); Gram-matrix (X Xᵀ, 144×144) PCA avoids the huge feature-space rotation matrix
m<-scale(t(logn),center=TRUE,scale=FALSE)               # samples × features, feature-centred
G<-tcrossprod(m); e<-eigen(G,symmetric=TRUE); ve<-round(100*e$values/sum(e$values),1)
sc<-e$vectors %*% diag(sqrt(pmax(e$values,0)))          # PC scores (n×n): PC_k = sc[,k]
out<-data.frame(sample=SAMP$sample,strain=SAMP$strain,tp=SAMP$tp,PC1=sc[,1],PC2=sc[,2],pc1_var=ve[1],pc2_var=ve[2],n_features=nrow(logn))
fwrite(out, file.path(OUT,"combined.pca.csv"))
cat("combined PCA done:",nrow(logn),"features,",ncol(M),"samples; PC1",ve[1],"% PC2",ve[2],"%\n")
