#!/usr/bin/env Rscript
# 16-strain per-timepoint PCA of piRNA-sequence expression (thesis Fig 5.21 / METHODS §7 method:
# DESeq2 1.42.1 size-factor normalisation -> prcomp on log2(norm+1)). Two feature sets: all_expressed
# (top-500 most-variable) and genuinely-unique (final_classified klass ^unique at this timepoint, union
# over the 16 strains; normalised with the same full-library size factors). Input = edger16/ (16-strain
# count matrices). Output = pca16/{tp}.pca.csv.
suppressMessages({library(edgeR); library(DESeq2); library(data.table)})
a<-commandArgs(TRUE); tp<-a[1]; lab<-a[2]
DIR<-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
FC <-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/unique16/final_classified_clean.csv.gz"  # CANONICAL 5-class (klass5): mm0-clean strain-private
OUT<-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pca16"
dir.create(OUT,showWarnings=FALSE)
cnt <-as.matrix(fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".counts.tsv.gz")))))
samp<-fread(file.path(DIR,paste0(tp,".samples.tsv")))
seqs<-fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".seqs.txt.gz"))),header=FALSE)$V1
colnames(cnt)<-samp$sample; strain<-factor(samp$strain)
y<-DGEList(cnt,group=strain,lib.size=samp$libsize_window); keep<-filterByExpr(y,group=strain)
cnt<-cnt[keep,,drop=FALSE]; seqs<-seqs[keep]
cat(sprintf("[%s] %d features after filterByExpr ; %d samples\n",lab,nrow(cnt),ncol(cnt)))
cd<-DataFrame(strain=strain,row.names=samp$sample)
dds<-DESeqDataSetFromMatrix(cnt,cd,~strain); dds<-estimateSizeFactors(dds)
logn<-log2(counts(dds,normalized=TRUE)+1); rownames(logn)<-seqs
fc<-fread(cmd=paste("zcat",FC)); uni<-unique(fc[grepl("^unique",klass5) & timepoint==tp]$sequence)  # genuinely-unique = CBS + mm0-clean strain-private (klass5; excludes SNP-variant AND low-quality)
uni<-intersect(uni,seqs)
do_pca<-function(M,tag){
  if(tag=="all_expressed"){v<-rowSums((M-rowMeans(M))^2); M<-M[order(-v)[1:min(500,nrow(M))],,drop=FALSE]}
  p<-prcomp(t(M),scale.=FALSE); ve<-round(100*p$sdev^2/sum(p$sdev^2),1)
  data.frame(sample=samp$sample,strain=samp$strain,rep=samp$rep,PC1=p$x[,1],PC2=p$x[,2],
             pc1_var=ve[1],pc2_var=ve[2],feature_set=tag,timepoint=lab,n_features=nrow(M))
}
res<-rbind(do_pca(logn,"all_expressed"),do_pca(logn,"all_full"),do_pca(logn[uni,,drop=FALSE],"unique"))   # all_full = ALL expressed piRNAs (no top-500 filter)
fwrite(res,file.path(OUT,paste0(tp,".pca.csv")))
cat(sprintf("[%s] unique features=%d ; wrote %s.pca.csv\n",lab,length(uni),tp))
