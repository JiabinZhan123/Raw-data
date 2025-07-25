##############GEO
rm(list = ls())
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
library(magrittr)
library(dplyr)
library(GEOquery)
library(data.table)
library(tidyr)
mg_count2TPMs=function(exp,gff3=NULL,gene_len=NULL){
  cut=c("__alignment_not_unique","__not_aligned","__too_low_aQual","__ambiguous","__no_feature")
  exp=exp[which(!row.names(exp)%in%cut),]
  if(is.null(gff3)&!is.null(gene_len)){
    exp1=exp[row.names(exp)%in%gene_len[,1],]
    anno=gene_len[match(row.names(exp1),gene_len[,1]),]
    lens=anno[,2]
    exp1=exp1/lens
    ct=apply(exp1,2,sum)
    exp1=t(t(exp1)/ct)
    exp1=exp1*1e6
    return(exp1)
  }else if(is.null(gff3)){
    anno=readMatrix(paste0(MG_Grobal_baseFolder,'/source/gencode.v22.ensg.genelen.tab'))
    row.names(anno)=gsub('\\..*','',row.names(anno))
    row.names(exp)=gsub('\\..*','',row.names(exp))
    exp1=exp[row.names(exp)%in%anno$V4,]
    anno=anno[match(row.names(exp1),anno$V4),]
    lens=anno[,4]
    exp1=exp1/lens
    ct=apply(exp1, 2, sum)
    exp1=t(t(exp1)/ct)
    exp1=exp1*1e6
    return(exp1)
  }else{
    library(GenomicFeatures)
    txdb <- makeTxDbFromGFF(gff3,format="gff3")
    exons_gene <- exonsBy(txdb, by = "gene")
    exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
    exons_gene_lens=unlist(exons_gene_lens)
    exp=exp[row.names(exp)%in%names(exons_gene_lens),]
    
    lens=exons_gene_lens[match(row.names(exp),names(exons_gene_lens))]
    exp1=exp/lens
    ct=apply(exp1, 2, sum)
    exp1=t(t(exp1)/ct)
    exp1=exp1*1e6
    return(exp1)
  }
}
##GSE75011 ############
#GSE75011=getGEO('GSE75011')
#GSE75011=GSE75011[[1]]
#saveRDS(GSE75011,file = "00_origin_datas/GEO/GSE75011.rds")
GSE75011=readRDS('00_origin_datas/GEO/GSE75011.rds')
GSE75011.pheno=pData(GSE75011)
head(GSE75011.pheno)
GSE75011.pheno=data.frame(Samples=GSE75011.pheno$title  ,
                          tissue=GSE75011.pheno$`disease group:ch1` )
rownames(GSE75011.pheno)=GSE75011.pheno$Samples 
head(GSE75011.pheno)
table(GSE75011.pheno$tissue)
GSE75011.pheno=GSE75011.pheno[which(GSE75011.pheno$tissue%in%c("Allergic rhinitis","Healthy non allergic")),]
GSE75011.pheno$tissue=as.vector(GSE75011.pheno$tissue)

GSE75011.pheno$group="Control"
GSE75011.pheno$group[which(GSE75011.pheno$tissue=="Allergic rhinitis")]="Case"
table(GSE75011.pheno$group)

# 
GSE75011.exp <- data.table::fread("00_origin_datas/GEO/GSE75011_Raw_counts.tsv/GSE75011_Raw_counts.tsv",data.table = F,header=T)
rownames(GSE75011.exp)=GSE75011.exp$V1
GSE75011.exp=GSE75011.exp[,-1]
GSE75011.exp[1:5,1:5]
dim(GSE75011.exp)
GSE75011.exp=GSE75011.exp[,as.vector(GSE75011.pheno$Samples)]
identical(colnames(GSE75011.exp),as.vector(GSE75011.pheno$Samples))
GSE75011.exp[1:5,1:5]
range(GSE75011.exp)
dim(GSE75011.exp)
exp_tpm <- mg_count2TPMs(GSE75011.exp)
exp_tpm[1:5,1:5]
range(exp_tpm)
exp_tpm=log2(exp_tpm+1)
GSE75011.exp=exp_tpm
save(GSE75011.exp,file='00_origin_datas/GEO/GSE75011.exp.RData')
save(GSE75011.pheno,file='00_origin_datas/GEO/GSE75011.pheno.RData')







###########GSE50223########
##GSE50223 ############
#GSE50223=getGEOExpData('GSE50223')
#save(GSE50223,file='00_origin_datas/GEO/GSE50223.RData')
load('00_origin_datas/GEO/GSE50223.RData')
GSE50223.pheno=GSE50223$Sample
head(GSE50223.pheno)
GSE50223.pheno=data.frame(Samples=GSE50223.pheno$Acc ,
                          tissue=GSE50223.pheno$diagnosis)
rownames(GSE50223.pheno)=GSE50223.pheno$Samples
head(GSE50223.pheno)
table(GSE50223.pheno$tissue)
GSE50223.pheno$tissue=as.vector(GSE50223.pheno$tissue)
GSE50223.pheno$group="Control"
GSE50223.pheno$group[which(GSE50223.pheno$tissue=="P")]="Case"
table(GSE50223.pheno$group)

# 
GSE50223.exp=GSE50223$Exp$GPL6884_48803_Data_col1
GSE50223.exp[1:5,1:5]
rownames(GSE50223.exp)
dim(GSE50223.exp)
GSE50223.exp=GSE50223.exp[,as.vector(GSE50223.pheno$Samples)]

identical(rownames(GSE50223.pheno),colnames(GSE50223.exp))

probe2ID=GSE50223$Anno$GPL6884[,c(1,13)]
GSE50223.exp=exp_probe2symbol_v2(datExpr = GSE50223.exp,anno = probe2ID,method = 'mean')
identical(colnames(GSE50223.exp),as.vector(GSE50223.pheno$Samples))
GSE50223.exp[1:5,1:5]
range(GSE50223.exp)
dim(GSE50223.exp)
boxplot(GSE50223.exp)
save(GSE50223.exp,file='00_origin_datas/GEO/GSE50223.exp.RData')
save(GSE50223.pheno,file='00_origin_datas/GEO/GSE50223.pheno.RData')




