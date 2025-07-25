rm(list = ls())

library(tidyverse)
library(cowplot)
library(clusterProfiler)
library(enrichplot)
library(plyr)
library(ggrepel)        
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(limma)
library(dplyr)
#"ODC1" "H1FX"
load('00_origin_datas/GEO/GSE75011.exp.RData')
load('00_origin_datas/GEO/GSE75011.pheno.RData')
identical(rownames(GSE75011.pheno), colnames(GSE75011.exp))
table(GSE75011.pheno$group)
geo_exp=as.data.frame(GSE75011.exp[,which(GSE75011.pheno$group=="Case")])
a<-geo_exp
a<-as.data.frame(t(a))
a<-a[order(a$H1FX, decreasing = T),]
a<-as.data.frame(t(a))
ncol(a)/2
list <- c(rep("high", 13), rep("low",12)) %>% factor(., levels = c("high", "low"), ordered = F)
list <- model.matrix(~factor(list)+0)
colnames(list) <- c("high", "low")
df.fit <- lmFit(a, list)
df.matrix <- makeContrasts(high - low, levels = list)
fit <- contrasts.fit(df.fit, df.matrix)
fit <- eBayes(fit)
tempOutput <- topTable(fit,coef=1,n = Inf, adjust = "BH")
nrDEG = na.omit(tempOutput)
diffsig <- nrDEG  
#write.csv(diffsig, "all.limmaOut.csv")

gsym.fc<-diffsig
gsym.fc$SYMBOL<-rownames(gsym.fc)
gsym.id <- bitr(rownames(diffsig), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)

gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]
id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID
id.fc[1:10]

kk=gseKEGG(
  id.fc,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = T,
  by = "fgsea")
saveRDS(kk,file = "04_GSEA/tcga.gsea.KEGG.RDS")
dim(kk)
head(kk)
kk.res=kk@result

table(kk.res$pvalue<0.05 & kk.res$NES<0)
table(kk.res$pvalue<0.05 & kk.res$NES>0)
kk.gsym <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T),]
sortkk[,1:2]
write.csv(sortkk,"04_GSEA/gsea_pathway_H1FX.csv", quote = F, row.names = F)

geneSetID <- c("hsa04370","hsa04662","hsa04917", "hsa04380", "hsa04926",
               "hsa04910","hsa04010","hsa03040")
suppressMessages(library(ggsci))

fig_kegg=enrichplot::gseaplot2(kk,geneSetID, pvalue_table = F,title ='KEGG enrichment')

##########################hallmark#########
pathway<-read.gmt("00_origin_datas/h.all.v2024.1.Hs.entrez.gmt")
y <- GSEA(id.fc,TERM2GENE =pathway,pvalueCutoff =0.05,minGSSize = 2)
hh.res=y@result

table(hh.res$pvalue<0.05 & hh.res$NES<0)
table(hh.res$pvalue<0.05 & hh.res$NES>0)

hh.gsym <- setReadable(y, 'org.Hs.eg.db', 'ENTREZID')

sorthh <- hh.gsym[order(hh.gsym$enrichmentScore, decreasing = T),]
write.csv(sorthh,"04_GSEA/gsea_hallmark_H1FX.csv", quote = F, row.names = F)
sorthh$ID
fig_hall=enrichplot::gseaplot2(y,sorthh$ID, pvalue_table = F,title ='Hallmark enrichment')

pdf('04_GSEA/fig_kegg.pdf', width = 9, height = 6)
fig_kegg
dev.off()

pdf('04_GSEA/fig_hall.pdf', width = 9, height = 6)
fig_hall
dev.off()
save.image("project_GSEA.Rdata")
