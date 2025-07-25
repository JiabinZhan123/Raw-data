rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(devtools))
suppressMessages(library(customLayout))
suppressMessages(library(stringr))
suppressMessages(library(ConsensusClusterPlus))
suppressMessages(library(tidydr))
suppressMessages(library(openxlsx))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(fgsea))
suppressMessages(library(corrplot))
suppressMessages(library(colorspace))
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(maftools))
suppressMessages(library(vegan))
suppressMessages(library(forcats))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(rstatix))
suppressMessages(library(ggstatsplot))
suppressMessages(library(ggcor))
suppressMessages(library(ggstance))
suppressMessages(library(tidyverse))
suppressMessages(library(GOplot))
suppressMessages(library(caret))
suppressMessages(library(writexl))
suppressMessages(library(rcartocolor))
suppressMessages(library(ggcorrplot))
suppressMessages(library(psych))
suppressMessages(library(clusterProfiler))
suppressMessages(library(dplyr))
suppressMessages(library(cols4all))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(scales))
suppressMessages(library(oncoPredict))
suppressMessages(library(gghalves))
suppressMessages(library(cowplot))
suppressMessages(library(IOBR))
suppressMessages(library(estimate))
suppressMessages(library(UpSetR))
suppressMessages(library(ggbiplot))
suppressMessages(library(ggsci))
suppressMessages(library(WGCNA))
suppressMessages(library(circlize))
suppressMessages(library(rJava))
suppressMessages(library(xlsxjars))
suppressMessages(library(xlsx))
suppressMessages(library(glmnet))
suppressMessages(library(tidyr))
suppressMessages(library(pROC))
suppressMessages(library(ROCR))
####GSE75011
load('00_origin_datas/GEO/GSE75011.exp.RData')
load('00_origin_datas/GEO/GSE75011.pheno.RData')


######## WGCNA ssGSEA########

exp <- as.data.frame(GSE75011.exp)
cli <- GSE75011.pheno
identical(as.vector(cli$Samples), colnames(exp))
table(cli$group)
genesets=read.table("00_origin_datas/277_PMID_37398672.txt",sep="\t",header=F,as.is=T,quote="\"",fill=T,check.names = F,stringsAsFactors = F)

geneset <- as.vector(genesets$V1)
geneset_list <- list()
geneset_list[['ssGSEA score']] <- as.vector(unique(geneset))
ssgsea_score <- ssGSEAScore_by_muti_group_genes(gene.exp = exp,
                                                genelist = geneset_list)
GSE_cli <- cbind(cli, t(ssgsea_score))
head(GSE_cli)
colnames(GSE_cli)=c("Samples" ,"tissue ", "Group" ,"ssGSEA score")
GSE_cli$Group=factor(GSE_cli$Group,levels = c("Control","Case" ))
fig1a <- mg_violin(GSE_cli[, c("Group", "ssGSEA score")]
                   ,melt = T
                   ,xlab = ''
                   ,legend.pos = 'tl'
                   ,ylab = 'ssGSEA score')
fig1a
dev.off()
ggsave(plot = fig1a,
       filename = '01_WGCNA/Fig1a.pdf',
       width = 5, height = 5)
write.table(t(ssgsea_score),file="01_WGCNA/ssgsea_score.txt", sep="\t", quote=F, row.names = T)

#############WGCNA#################
ssgsea_score <- read.table("01_WGCNA/ssgsea_score.txt",header = T,check.names = F,fill=T,sep = "\t")

allowWGCNAThreads(nThreads = 36)
enableWGCNAThreads(nThreads = 36)
tcga_exp=exp

tcga.mads=apply(tcga_exp, 1, mad)
tpm_T2=tcga_exp[which(tcga.mads>quantile(tcga.mads, probs=seq(0, 1, 0.25))[2]),]
dim(tpm_T2)
tpm_T2=t(tpm_T2)
dim(tpm_T2)
range(tpm_T2)
pdf('01_WGCNA/FigABC.pdf',width = 8,height = 8)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()

tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.25,
                                 minModuleSize=80)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))

pdf('01_WGCNA/FigD.pdf',height = 6,width = 6)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

zx=data.frame(tpm_T2.module$Modules)
write.table(zx,'01_WGCNA/tcga.wgcna.module.genes.txt',sep = "\t", row.names = TRUE,quote = FALSE,col.names = T)

fig2e=mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                       ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                       ,point_sizes = 2
                       ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                       ,xlab = 'Number of Genes',legend.pos = NULL)
fig2e
dev.off()
ggsave('01_WGCNA/Fig2e.pdf',fig2e,height = 6,width = 6)

##

# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('01_WGCNA/Fig2f.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()


identical(rownames(tpm_T2),rownames(ssgsea_score))
spms <- ssgsea_score
MEs_col<-tpm_T2.module$MEs
dim(MEs_col)
identical(rownames(spms),rownames(MEs_col))
modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'p')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])

textMatrix = paste(signif(modTraitCor, 2), "\n", " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)


pdf('01_WGCNA/Figf-2.pdf',width =5,height = 7)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = rownames(modTraitCor), 
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))

geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))


module = "pink"
#module = "purple"

pheno = "ssGSEA score"
modNames<-colnames(geneModuleMembership)

module_column = match(module, modNames)
pheno_column = match(pheno, colnames(spms))

gene_M=as.data.frame(tpm_T2.module$Modules)
identical(rownames(geneModuleMembership),rownames(gene_M))
moduleGenes <-  gene_M[,2]== module

#par(mar = c(1, 3, 3, 3))
pdf('01_WGCNA/GS_MM.pdf',width = 6,height = 7)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitSignificance[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance"),cex=1,pch = 16,
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'pink')
abline(v=0.5,h=0.3,col="red")
dev.off()

zx1=as.data.frame(geneModuleMembership[moduleGenes, module_column])
rownames(zx1)=rownames(geneModuleMembership)[moduleGenes]

zx2=as.data.frame(geneTraitSignificance[moduleGenes, 1])
rownames(zx2)=rownames(geneTraitSignificance)[moduleGenes]
zx3=intersect(rownames(zx1)[abs(zx1[,1])>0.5],
              rownames(zx2)[abs(zx2[,1])>0.3])

length(zx3)
write.table(zx3, file="01_WGCNA/WGCNA_genes.txt",quote = FALSE, row.names=F, col.names=F)

library(rJava)
library(xlsxjars)
library(xlsx)

write.table(geneModuleMembership, file="01_WGCNA/geneModuleMembership.txt",quote = FALSE, row.names=T)
write.table(MMPvalue, file="01_WGCNA/MMPvalue.txt",quote = FALSE,   row.names=T)
write.table(geneTraitSignificance, file="01_WGCNA/geneTraitSignificance.txt",quote = FALSE, row.names=T)
write.table(GSPvalue, file="01_WGCNA/GSPvalue.txt",quote = FALSE,  row.names=T)






##############DEG##########
load('00_origin_datas/GEO/GSE75011.exp.RData')
load('00_origin_datas/GEO/GSE75011.pheno.RData')
identical(rownames(GSE75011.pheno), colnames(GSE75011.exp))

table(GSE75011.pheno$group)

geo.limma=mg_limma_DEG(exp = GSE75011.exp,
                       group = GSE75011.pheno$group,
                       ulab = 'Case',dlab = 'Control')
geo.limma$Summary

geo.degs=geo.limma$DEG[which(geo.limma$DEG$P.Value<0.05 & abs(geo.limma$DEG$logFC) > log2(1.5)),]
head(geo.degs)
geo.degs$group=ifelse(geo.degs$logFC>0,'up','down')
geo.degs$gene=rownames(geo.degs)
table(geo.degs$group)
write.table(geo.degs,file="02_DEGs//diffSig_mRNA_log2FC1.5_P.Vale0.05.txt",sep="\t",quote=F)

######################DEGsplot

fig2a=my_volcano(geo.limma,p_cutoff = 0.05,fc_cutoff = log2(1.5),col = c("#FDB462","#8BACD1",'grey'))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        text = element_text(color = "black",family = 'Times',size = 14),
        axis.text = element_text(color = "black",family = 'Times',size = 14),
        legend.position = 'top')
fig2a
ggsave('02_DEGs/fig2a.pdf',fig2a,height = 6,width = 6.5)

dev.off()


enrichment=mg_clusterProfiler(as.vector(rownames(geo.degs)))
write.table(enrichment$Enrich_tab,file = '02_DEGs/enrichment.txt',sep = '\t',quote = F,row.names = T,col.names = T)

kegg_dot=enrichplot::dotplot(enrichment$KEGG)+scale_y_discrete(labels=function(y)str_wrap(y,width = 25))+ scale_x_continuous(limits = c(0.05, 0.18))

ggsave('02_DEGs/kegg_dot.pdf',kegg_dot,height = 6,width = 8)


hub=intersect(rownames(geo.degs),zx3)
length(hub)
write.table(hub,file="02_DEGs//hub.txt",sep="\t",quote=F,row.names = F,col.names = F)


pdf(file="02_DEGs/Venn.pdf",width=7,height=6)
mg_venn_plot(list("DEGs"=rownames(geo.degs)
                  , "WGCNA_genes" = zx3))
dev.off()



############################Lasso##########################

library(glmnet)
library(tidyr)
library(pROC)

exp=GSE75011.exp
cli=GSE75011.pheno
identical(as.vector(cli$Samples), colnames(exp))
hub

OADE=exp[which(rownames(exp)%in%hub),]
table(cli$group)
identical(as.vector(cli$Samples), colnames(OADE))

data_ML <- as.data.frame(t(OADE))
group.v <- cli$group  %>% factor(.,levels = c("Control","Case"),ordered = F)
design.v <- model.matrix(~factor(group.v))   
data_ML<-cbind(design.v[,2],data_ML)
colnames(data_ML)[1]<- 'group'
data_ML <- as.matrix(data_ML)

set.seed(123456)

fit_stage.v <-cv.glmnet(as.matrix(data_ML[,-1]) ,as.factor(data_ML[,'group']),
                        family="binomial",type.measure="default",nfolds = 10, gamma = 1)
fit.v <- fit_stage.v$glmnet.fit
pdf(file = "03_Lasso/lasso.Binomial.Deviance.pdf",height = 8,width = 9)
par(mgp = c(2.5,1,0),mar=c(4,5,3,3))
plot(fit_stage.v,xlab='Log Lambda',cex.lab = 1)+
  text(x = log(fit_stage.v$lambda.min),y = 3,
       paste('Lambda.min\n',round(fit_stage.v$lambda.min,4)),cex=0.8,adj=0.1)+
  text(x = log(fit_stage.v$lambda.1se),y = 1.5,
       paste('Lambda.lse\n',round(fit_stage.v$lambda.1se,4)),cex=1)
dev.off()

pdf(file = "03_Lasso/lasso.voefficients.venalty.pdf",height = 8,width = 10)
par(mgp = c(4,1,0),mai=c(2,2,1,1))
plot(fit.v, xvar="lambda",cex.lab = 1)+
  abline(v = c(log(fit_stage.v$lambda.min), log(fit_stage.v$lambda.1se)),lty=2)+
  text(x = log(fit_stage.v$lambda.min),y = 6,
       paste('Lambda.min\n',round(fit_stage.v$lambda.min,4)),cex=1,adj=0.1)+
  text(x = log(fit_stage.v$lambda.1se),y = -6,
       paste('Lambda.lse\n',round(fit_stage.v$lambda.1se,4)),cex=1,adj=0.1)
dev.off()

lam = fit_stage.v$lambda.min

coefficients<-coef(fit_stage.v,s=lam)
Active.Index<-coefficients@i
Active.name.train <- colnames(data_ML[,-1])[Active.Index[-1]]
write.csv(Active.name.train,"03_Lasso/lasso.csv")


##SVM

library(e1071)
library(kernlab)
library(caret)

set.seed(123456)
Profile=rfe(x=(data_ML[,-1]),
            y=(as.numeric( data_ML[,'group'])),
            sizes =c(2,4,6,8, seq(10,50,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")
x = Profile$results$Variables
y = Profile$results$RMSE
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
pdf(file="03_Lasso/SVM-RFE.pdf", width=8, height=8)
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
points(wmin.x, wmin.y, col="darkblue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=1, col=1)
dev.off()

featureGenes=Profile$optVariables
write.table(file="03_Lasso/SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)



sig=intersect(featureGenes,Active.name.train)

###############
write.table(file="03_Lasso/signaturegenes.txt", sig, sep="\t", quote=F, row.names=F, col.names=F)

pdf(file="03_Lasso/Venn_ML.pdf",width=7,height=6)
mg_venn_plot(list("Lasso"=Active.name.train
                  , "SVM_RFE" = featureGenes 
))
dev.off()

####################
identical(as.vector(cli$Samples), colnames(exp))
#colnames(cli)=c("Samples","group")

sigexp=exp[which(rownames(exp)%in%sig),]
sigexp=t(sigexp)
identical(rownames(sigexp), rownames(cli))

sigexp <- cbind(cli,sigexp)
library(reshape2)
library(gghalves)
sigexp1 <- melt(sigexp,id.vars=c("Samples","group"),
                measure.vars = c(sig),variable.name = "Signaturegene",value.name = "Expression")

sig_plot=ggplot(sigexp1,aes(x=Signaturegene,y=Expression,fill=group))+geom_boxplot()+
  xlab("Signaturegenes")+ylab("Expression")+theme(axis.ticks = element_blank())+
  ggsci::scale_fill_nejm()+theme_bw()+ggsci::scale_fill_nejm (name="group")+
  ggpubr::stat_compare_means(method = 't.test',label= "p.signif")

sig_plot
dev.off()
ggsave(plot = sig_plot, filename = '03_Lasso/sigexp.pdf',width = 6, height = 5)


####################
##################glm

temp_data=exp[which(rownames(exp)%in%sig),]

temp_data <- as.data.frame(t(temp_data))
temp_data$"Samples"=rownames(temp_data)
temp_data=merge(temp_data,cli,by="Samples")
rownames(temp_data)=temp_data$Samples
#data=OADE[,c(ncol(OADE),2:(ncol(OADE)-1))]
data=temp_data[,c("group",sig)]
data <- na.omit(data)
data$group=factor(data$group,levels = c("Control","Case"))

##########################ROC


rocobj <- roc(group ~ ., data = data)
rocobj_legend <- c()
for (l in 1:length(rocobj)) {
  auc <- auc(rocobj[[l]])
  legend <- paste0(colnames(data)[l+1],"_AUC: ",round(auc,3))
  rocobj_legend <- c(rocobj_legend,legend)
}
rocobj_legend

cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
pdf(file = paste0("03_Lasso/roc_train", ".pdf"), height = 8, width = 8, onefile = FALSE)
ggroc(rocobj,size = 1.1,legacy.axes=T) + 
  labs(x = "False positiverate", y = "True positiverate",)+
  scale_color_manual(values = cbPalette[1:length(sig)],breaks = colnames(data)[-1],labels = rocobj_legend)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed",size = 1)+
  theme(legend.title=element_blank(), legend.position = c(0.7,0.25), legend.text = element_text(size = 18),text = element_text(size = 20),axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"), axis.title = element_text(size = 20, hjust = 0.5, colour = "black"), axis.text = element_text(size = 20, color = "black"), panel.background = element_rect(fill = NA),panel.border = element_rect(fill = NA, size = 1.5))
dev.off()


#############GSE50223#############


load('00_origin_datas/GEO/GSE50223.exp.RData')
load('00_origin_datas/GEO/GSE50223.pheno.RData')
cepu <-GSE50223.exp
cecli <- GSE50223.pheno[,c(1,3)]
table(cecli$group)
identical(colnames(cepu), rownames(cecli))

cepu <- cepu[,as.vector(cecli$Samples)]
range(cepu)
identical(as.vector(cecli$Samples), colnames(cepu))
sig

cepu <- as.data.frame(t(cepu))
which(colnames(cepu)%in%sig)
cepu <- cepu[,which(colnames(cepu)%in%sig)]
cepu$"Samples"=rownames(cepu)
cepu=merge(cepu,cecli,by="Samples")
rownames(cepu)=cepu$geo_accession
sigcepu1 <- melt(cepu,id.vars=c("Samples","group"),
                 measure.vars = c(sig),variable.name = "Signaturegene",value.name = "Expression")

sigcepu_plot=ggplot(sigcepu1,aes(x=Signaturegene,y=Expression,fill=group))+geom_boxplot()+
  xlab("Signaturegenes")+ylab("Expression")+theme(axis.ticks = element_blank())+
  ggsci::scale_fill_nejm()+theme_bw() +ggsci::scale_fill_nejm(name="group")+
  ggpubr::stat_compare_means(method = 't.test',label= "p.signif")
ggsave(plot = sigcepu_plot,
       filename = '03_Lasso/sigcepu_plot_GSE20347.pdf',
       width = 6, height = 5)


data=cepu[,c(ncol(cepu),2:(ncol(cepu)-1))]
data$group=factor(data$group,levels = c("Control","Case"))


rocobj <- roc(group ~ ., data = data)
rocobj_legend <- c()
for (l in 1:length(rocobj)) {
  auc <- auc(rocobj[[l]])
  legend <- paste0(colnames(data)[l+1],"_AUC: ",round(auc,3))
  rocobj_legend <- c(rocobj_legend,legend)
}
rocobj_legend

cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
pdf(file = paste0("03_Lasso/roc_train_GSEtest", ".pdf"), height = 8, width = 8, onefile = FALSE)
ggroc(rocobj,size = 1.1,legacy.axes=T) + 
  labs(x = "False positiverate", y = "True positiverate",)+
  scale_color_manual(values = cbPalette[1:length(sig)],breaks = colnames(data)[-1],labels = rocobj_legend)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed",size = 1)+
  theme(legend.title=element_blank(),
        legend.position = c(0.7,0.25),
        legend.text = element_text(size = 18, margin = margin (b = 12, unit = "pt")),
        text = element_text(size = 20),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title = element_text(size = 20, hjust = 0.5, colour = "black"),
        axis.text = element_text(size = 20, color = "black"),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, size = 1.5))
#scale_color_discrete(breaks = colnames(data)[-1],labels = rocobj_legend)
dev.off()




# #"ODC1" "H1FX"
# load('00_origin_datas/GEO/GSE75011.exp.RData')
# load('00_origin_datas/GEO/GSE75011.pheno.RData')
# identical(rownames(GSE75011.pheno), colnames(GSE75011.exp))
# 
# table(GSE75011.pheno$group)
# geo_exp=as.data.frame(GSE75011.exp[,which(GSE75011.pheno$group=="Case")])
# geo_cli=GSE75011.pheno[colnames(geo_exp),]
# site=which(rownames(geo_exp)=="H1FX")
# med=median(as.numeric(geo_exp[site,]))
# identical(rownames(geo_cli), colnames(geo_exp))
# which(geo_exp[site,]>=med)
# which(geo_exp[site,]<med)
# geo_cli$"Type"="High"
# geo_cli$"Type"[which(geo_exp[site,]<med)]="Low"
# table(geo_cli$"Type")
# tcga.geneList=getGeneFC(gene.exp=geo_exp[,as.vector(geo_cli$Samples)], group=geo_cli$Type,ulab='High',dlab='Low')
# tcga.geneList[1:10]
# 
# saveRDS(tcga.geneList,file = "04_GSEA/tcga.geneList.RDS")
#############imm##########
####immune_ #############
load('00_origin_datas/GEO/GSE75011.exp.RData')
load('00_origin_datas/GEO/GSE75011.pheno.RData')
identical(rownames(GSE75011.pheno), colnames(GSE75011.exp))


cli=GSE75011.pheno
exp=GSE75011.exp
colnames(cli)[1]="Samples"
head(cli)
cli$Samples=as.vector(cli$Samples)
table(cli$"group")
identical(colnames(exp),as.vector(cli$Samples))

geo.immu.ssgsea=immu_ssgsea(exp = exp,isTCGA=F)
head(geo.immu.ssgsea)
saveRDS(geo.immu.ssgsea,"05_imm/geo.immu.ssgsea.rds")
geo.immu.ssgsea <- readRDS("05_imm/geo.immu.ssgsea.rds")

head(geo.immu.ssgsea)
tme.df2=geo.immu.ssgsea[cli$Samples,]
tme.df2=as.data.frame(tme.df2)
tme.df2$type=cli$group
tme.df2=melt(tme.df2)
head(tme.df2)
pdf('05_imm/Fig5a.pdf',height = 5,width = 9)
ggplot(tme.df2,aes(x=variable,y=value,fill=type))+
  geom_boxplot()+stat_compare_means(aes(group=type), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =c("#CE3D32","#5050FF"))+
  xlab('')+ylab('Score')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   axis.text.x = element_text(color = "black", size = 12,angle = 90,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
dev.off()
cli_anno=cli[order(cli$group),'group',drop=F]
cli$Samples=as.vector(cli$Samples)
tme.df1=data.frame(t(exp[sig,cli$Samples[cli$group=='Case']]),
                   geo.immu.ssgsea[cli$Samples[cli$group=='Case'],])

colnames(tme.df1)=gsub('[.]',' ',colnames(tme.df1))
cor_res <- Hmisc::rcorr(as.matrix(tme.df1),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 0
library(corrplot)

pdf('05_imm/Fig5b.pdf',height = 6,width = 10.5,onefile = F)
corrplot(as.matrix((cor_res$r[sig,-c(1:length(sig))])),
         p.mat = as.matrix((cor_res$P[sig,-c(1:length(sig))])),
         mar = c(0,0,1,1),
         col=colorRampPalette(c('#226ED1', 'white','#E01010'))(10),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[2],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[1],
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()


key.genes_immunescore_cor_list <- list()
colnames(geo.immu.ssgsea)
geo.immu.ssgsea=as.data.frame(geo.immu.ssgsea)

list_imm=c("Activated B cell","Activated CD4 T cell","Central memory CD4 T cell" )
for (ge in list_imm) {
  print(ge)
  tmp <- cor_point(x = as.numeric(geo.immu.ssgsea[cli$Samples[cli$group=='Case'], ge]),
                   y = as.numeric(exp["H1FX", cli$Samples[cli$group=='Case']]),
                   method='Pearson',
                   #top_col=mycolors[1],
                   #right_col=mycolors[2],
                   ylab=paste0("H1FX", ' expression'),
                   xlab=paste0(ge, ' Score'))
  key.genes_immunescore_cor_list[[ge]] <- tmp
}
pdf('05_imm/key.genes_cor_plot.pdf', width = 12, height = 4)
cowplot::plot_grid(plotlist = key.genes_immunescore_cor_list,
                   ncol = length(list_imm),nrow = 1)
dev.off()
library('enrichR')
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyverse) 
symbol=sig[2]

dbs <- listEnrichrDbs()
dbs<- c("DSigDB")
enrichr<- enrichr(symbol, dbs)
enrichr_DSigDB=enrichr$DSigDB
#enrichr_DSigDB=enrichr_DSigDB[which(enrichr_DSigDB$P.value<0.05),]

write.table(enrichr_DSigDB,file = '06_drug/DSigDB_result.txt',sep = '\t',quote = F,row.names = F,col.names = T)

