#!/usr/bin/env Rscript
# Per-timepoint PCA of piRNA-sequence expression across the 3 pilot strains (thesis Fig 5.21 method:
# DESeq2 size-factor normalisation -> PCA). TWO feature sets, side by side:
#   all_expressed = every 24-32 nt piRNA passing filterByExpr (top-500 most-variable, DESeq2 plotPCA
#                   convention) -> reproduces Fig 5.21 with our data.
#   unique        = genuinely-unique (Step-4) piRNAs specific at this timepoint (union over strains),
#                   ALL of them. Normalised with the SAME full-library size factors (valid; the unique
#                   set alone would violate DESeq2's most-features-not-DE assumption).
# DESeq2 1.42.1 (= thesis v1.42). filterByExpr = the data-driven floor used in the DA step.
suppressMessages({library(edgeR); library(DESeq2); library(data.table)})
a<-commandArgs(TRUE); tp<-a[1]; lab<-a[2]
DIR<-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger"
S4 <-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4"
OUT<-"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pca"
dir.create(OUT,showWarnings=FALSE)
cnt <-as.matrix(fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".counts.tsv.gz")))))
samp<-fread(file.path(DIR,paste0(tp,".samples.tsv")))
seqs<-fread(cmd=paste("zcat",file.path(DIR,paste0(tp,".seqs.txt.gz"))),header=FALSE)$V1
colnames(cnt)<-samp$sample; strain<-factor(samp$strain)
y<-DGEList(cnt,group=strain,lib.size=samp$libsize_window); keep<-filterByExpr(y,group=strain)
cnt<-cnt[keep,,drop=FALSE]; seqs<-seqs[keep]
cat(sprintf("[%s] %d features after filterByExpr\n",lab,nrow(cnt)))
cd<-DataFrame(strain=strain,row.names=samp$sample)
dds<-DESeqDataSetFromMatrix(cnt,cd,~strain); dds<-estimateSizeFactors(dds)
logn<-log2(counts(dds,normalized=TRUE)+1); rownames(logn)<-seqs
# genuinely-unique seqs specific at this timepoint (union over strains)
uni<-character(0)
for(X in c("C57BL_6NJ","CAST_EiJ","SPRET_EiJ")){
  d<-fread(cmd=paste("zcat",file.path(S4,paste0(X,".step4_classified.csv.gz"))))
  d<-d[grepl("^unique",klass) & grepl(lab,timepoints,fixed=TRUE)]
  uni<-union(uni,d$sequence)
}
uni<-intersect(uni,seqs)
do_pca<-function(M,tag){
  if(tag=="all_expressed"){v<-rowSums((M-rowMeans(M))^2); M<-M[order(-v)[1:min(500,nrow(M))],,drop=FALSE]}
  p<-prcomp(t(M),scale.=FALSE); ve<-round(100*p$sdev^2/sum(p$sdev^2),1)
  data.frame(sample=samp$sample,strain=samp$strain,rep=samp$rep,PC1=p$x[,1],PC2=p$x[,2],
             pc1_var=ve[1],pc2_var=ve[2],feature_set=tag,timepoint=lab,n_features=nrow(M))
}
res<-rbind(do_pca(logn,"all_expressed"),do_pca(logn[uni,,drop=FALSE],"unique"))
fwrite(res,file.path(OUT,paste0(tp,".pca.csv")))
cat(sprintf("[%s] unique features=%d ; wrote %s.pca.csv\n",lab,length(uni),tp))
