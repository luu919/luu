#DGEs#
rt <- read.table("GLTP_exp.txt", sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
rt_transposed <- t(rt)
write.table(rt_transposed, file = "GLTP_exp_transposed.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
set.seed(123456)
kmeans_result <- kmeans(rt_transposed, centers = 2, nstart = 25)
cluster_data <- data.frame(Sample = rownames(rt_transposed), Cluster = factor(kmeans_result$cluster))
write.table(cluster_data, file = "cluster_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
BangerBox_data <- read.table("BangerBox_log2.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
transposed_cluster <- read.table("transposed_cluster.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
Merged_Data <- rbind(transposed_cluster, BangerBox_data)
write.table(Merged_Data, "fenzu_all_exp_Data.txt", sep = "\t", quote = FALSE, row.names = TRUE)
BangerBox_data_high <- read.table("Exp_allgene_gltphigh.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
BangerBox_data_low <- read.table("Exp_allgene_gltplow.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
BangerBox_data_high[, -1] <- sapply(BangerBox_data_high[, -1], function(x) as.numeric(as.character(x)))
BangerBox_data_low[, -1] <- sapply(BangerBox_data_low[, -1], function(x) as.numeric(as.character(x)))
t_test_results_all_gene <- data.frame(
  Gene = rownames(BangerBox_data_high), 
  p_value = rep(NA, nrow(BangerBox_data_high)), 
  logFC = rep(NA, nrow(BangerBox_data_high))
)
for (i in 1:nrow(BangerBox_data_high)) {
  low_expr <- as.numeric(as.character(BangerBox_data_low[i, -1]))
  high_expr <- as.numeric(as.character(BangerBox_data_high[i, -1]))
  
  if (any(is.na(low_expr)) || any(is.na(high_expr))) next
  
  t_test_result <- t.test(high_expr, low_expr, na.rm = TRUE)
  t_test_results_all_gene$p_value[i] <- t_test_result$p.value
  
  mean_low <- mean(low_expr, na.rm = TRUE)
  mean_high <- mean(high_expr, na.rm = TRUE)
  
  if (!is.na(mean_low) && !is.na(mean_high) && mean_low != 0 && mean_high != 0) {
    t_test_results_all_gene$logFC[i] <- log2(mean_low/mean_high)
  }
}
t_test_results_all_gene$adj_p_value <- p.adjust(t_test_results_all_gene$p_value, method = "fdr")
write.table(t_test_results_all_gene, file = "t_test_results_all.txt", sep = "\t", row.names = TRUE, quote = FALSE)
t_test_results_all_gene <- read.table("t_test_results_all.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
t_test_results_saixuan <- t_test_results_all_gene[t_test_results_all_gene$adj_p_value < 0.05 & abs(t_test_results_all_gene$logFC) > 0.58496, ]
write.table(t_test_results_saixuan, file = "t_test_results_saixuan_1.2.txt", sep = "\t", row.names = TRUE, quote = FALSE)





#Enrichment analysis#
t_test_results <- read.table("t_test_results_all.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
t_test_results$neg_log10_adj_p_val <- -log10(t_test_results$adj_p_value)
t_test_results <- t_test_results[t_test_results$neg_log10_adj_p_val <= 50, ]
t_test_results$color <- ifelse(t_test_results$logFC > 0.2630344058337 & t_test_results$adj_p_value < 0.05, "Up",                              ifelse(t_test_results$logFC < -0.2630344058337 & t_test_results$adj_p_value < 0.05, "Down", "Not sig"))
write.table(t_test_results, file = "t_test_results_6ge.txt", sep = "\t", row.names = FALSE, quote = FALSE)
library(ggplot2)
ggplot(t_test_results, aes(x = logFC, y = -log10(adj_p_value), color = color)) +
  geom_point(alpha = 0.6, size = 1.2) + 
  scale_color_manual(values = c("Up" = "#F98787", "Down" = "#70DEF4", "Not sig" = "#C0C0C0")) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + 
  theme(panel.grid.major = element_line(color = "grey90", size = 0.5), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.line = element_line(color = "black"), 
        legend.position = "right") + 
  labs(x = "Log2(Fold Change)", y = "-Log10(FDR)")
# GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
data = read.table("t_test_results_all.txt", header = TRUE, sep = "\t")
colnames(data)[1] = "SYMBOL"
gene = data$SYMBOL
gene = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
gene = dplyr::distinct(gene, SYMBOL, .keep_all = TRUE)
data_all <- data %>% inner_join(gene, by = "SYMBOL")
data_all_sort <- data_all %>% arrange(desc(logFC))
geneList = data_all_sort$logFC
names(geneList) <- data_all_sort$ENTREZID
gsea <- gseKEGG(geneList, organism = "hsa", pvalueCutoff = 0.05)
gsea <- setReadable(gsea, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
dotplot(gsea)
ridgeplot(gsea, label_format = 100)
custom_colors <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C")
gseaplot2(gsea, 1:5, pvalue_table = T, color = custom_colors)
# GO
rt <- read.table("t_test_results_saixuan_12.txt", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
genes = rownames(rt)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
out <- cbind(rt, entrezID = entrezIDs)
write.table(out, file = "GO-id.txt", sep = "\t", quote = F, row.names = T)
rt <- read.table("GO-id.txt", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
gene = rt$entrezID
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "all", readable = T)
write.table(kk, file = "GO.txt", sep = "\t", quote = F, row.names = F)
pdf(file = "GO-bubble13_8.pdf", width = 13, height = 10)
dotplot(kk, showCategory = 10, label_format = 100, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = 'free')
dev.off()
# KEGG
rt <- read.table("t_test_results_saixuan_12.txt", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
genes = rownames(rt)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
out = cbind(rt, entrezID = entrezIDs)
write.table(out, file = "KEGG-id.txt", sep = "\t", quote = F, row.names = T)
rt <- read.table("KEGG-id.txt", sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
gene = rt$entrezID
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.table(kk, file = "KEGG.txt", sep = "\t", quote = F, row.names = F)
ego <- read.table("KEGG.txt", sep = "\t", check.names = F, header = T)
go = data.frame(Category = "ALL", ID = ego$ID, Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)
id.fc = rt
genelist <- data.frame(ID = id.fc$entrezID, logFC = id.fc$logFC, gene = rownames(id.fc))
library(GOplot)
circ <- circle_dat(go, genelist)
termNum = 5
geneNum = nrow(genelist)
chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
genelist = na.omit(genelist)
rownames(genelist) = genelist$ID
sameid = intersect(rownames(chord), rownames(genelist))
chord = chord[sameid,]
genelist = genelist[sameid,]
rownames(chord) = genelist$gene
pdf(file = "KEGG_circ.pdf", width = 12, height = 11)
GOChord(chord, space = 0.001, gene.order = 'logFC', gene.space = 0.15, gene.size = 1, border.size = 0.1, process.label = 9)
dev.off()







#Characterisation of the immune landscape
setwd("TCGA-HNSC")
setwd("CIBERSORT")
library(e1071)
library(preprocessCore)
library(tidyverse)
source("CIBERSORT.R")
sig_matrix <- "LM22.txt"
mixture_file = 'BangerBox_log2.txt'
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
res_cibersort <- res_cibersort[,1:22]
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]
ciber.res <- as.data.frame(ciber.res)
write.table(ciber.res,"ciber.res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7)
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)), border = NA, names.arg = rep("",nrow(ciber.res)), 
        yaxt = "n", ylab = "Relative percentage", col = mycol)
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, par("usr")[4], legend = colnames(ciber.res), xpd = T,
       fill = mycol, cex = 0.6, border = NA, y.intersp = 1, x.intersp = 0.2, bty = "n")

library(tidyr)
library(reshape2)
library(ggpubr)
a <- read.csv("CIBERSORT_Results.csv", row.names=1, header=T)
a$X <- rownames(a)
fen <- read.csv("group.csv")
a <- merge(fen, a, "X")
mydata1 <- melt(a, id.vars=c("X","Group"), variable.name="immunecell", value.name="tpm")
colnames(mydata1) <- c("Sample", "Groups", "immunecell","tpm")

pvalues <- sapply(unique(mydata1$immunecell), function(x) {
  res <- t.test(tpm ~ Groups, data = subset(mydata1, immunecell == x))
  res$p.value
})
pv <- data.frame(gene = unique(mydata1$immunecell), pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001,0.001,0.01,0.05,1), labels=c('****','***','**','*','ns'))

p.box <- ggplot(mydata1, aes(x=immunecell, y=tpm, color=Groups, fill=Groups)) +
  geom_boxplot(alpha = .7) + theme_classic() + 
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_text(aes(x=gene, y=max(mydata1$tpm)*1.1, label=sigcode), data=pv, inherit.aes=F) +
  ylab("Immunecell expression")
p.box


data <- read.csv("gltp_exp.csv", row.names=1, check.names=FALSE)
transposed_data <- t(data)
colnames(transposed_data) <- colnames(data)
write.csv(transposed_data, "transposed_gltp_exp.csv", row.names=TRUE)
a <- read.csv("transposed_gltp_exp.csv", row.names=1, header=T)
a$X <- rownames(a)
fen <- read.csv("group.csv")
a <- merge(fen, a, "X")
mydata1 <- melt(a, id.vars=c("X","Groups"), variable.name="gene", value.name="tpm")
colnames(mydata1) <- c("Sample", "Groups", "gene","tpm")

pvalues <- sapply(unique(mydata1$gene), function(x) {
  res <- t.test(tpm ~ Groups, data = subset(mydata1, gene == x))
  res$p.value
})
pv <- data.frame(gene = unique(mydata1$gene), pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001,0.001,0.01,0.05,1), labels=c('****','***','**','*','ns'))

p.violin <- ggplot(mydata1, aes(x=gene, y=tpm, color=Groups, fill=Groups)) +
  geom_violin(alpha = .5) + theme_classic() + 
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_text(aes(x=gene, y=max(mydata1$tpm)*1.1, label=sigcode), data=pv, inherit.aes=F) +
  ylab("GLTP expression")
p.violin


a <- read.csv("transposed_t_test_FDR_HLA.csv", row.names=1, header=T)
a$X <- rownames(a)
fen <- read.csv("group.csv")
a <- merge(fen, a, "X")
mydata1 <- melt(a, id.vars=c("X","Group"), variable.name="gene", value.name="tpm")
colnames(mydata1) <- c("Sample", "Groups", "gene","tpm")

pvalues <- sapply(unique(mydata1$gene), function(x) {
  res <- t.test(tpm ~ Groups, data = subset(mydata1, gene == x))
  res$p.value
})
pv <- data.frame(gene = unique(mydata1$gene), pvalue = pvalues)
pv$sigcode <- cut(pv$pvalue, c(0,0.0001,0.001,0.01,0.05,1), labels=c('****','***','**','*','ns'))

mydata1$gene <- gsub("\\.", "-", mydata1$gene)
pv$gene <- gsub("\\.", "-", pv$gene)

p.box <- ggplot(mydata1, aes(x=gene, y=tpm, color=Groups, fill=Groups)) +
  geom_boxplot(alpha = .5) + theme_classic() + 
  scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_text(aes(x=gene, y=max(mydata1$tpm)*1.1, label=sigcode), data=pv, inherit.aes=F) +
  ylab("Expression of HLA related genes")
p.box

#Prognostic model
library(pheatmap)
library(survival)
library(tidyverse)
library(glmnet)
library(caret)
library(limma)
library(survminer)
library(forestplot)

pFilter=0.01
exp=read.table("BangerBox_log2.txt",row.names=1,header=T,sep="\t",check.names=F)
exp_tumor= exp%>%as.data.frame()%>%dplyr::select(str_which(colnames(.), "-01|-06"))
hdiff=read.table("t_test_results_saixuan_1.2.txt",header=T,sep="\t",row.names=1)
diffexp=exp_tumor[rownames(hdiff),]
diffexp=t(diffexp)
rownames(diffexp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(diffexp))
clinical=read.table("clinical_1.txt",header=T,sep="\t",row.names=1)
samesample=intersect(rownames(diffexp),rownames(clinical))
diffexp=as.data.frame(diffexp[samesample,])
clinical=as.data.frame(clinical[samesample,])
rt=data.frame(clinical[,1:2],diffexp)
rt$futime=rt$futime/365
write.table(rt,file="clinical and exp.txt",sep="\t",row.names=T,quote=F)

rt=read.table("clinical and exp.txt",row.names=1,header=T,sep="\t",check.names=F)

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
  cox<-coxph(Surv(futime,fustat)~rt[,i],data=rt)
  coxSummary=summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
  }
}

write.table(outTab,file="UniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp_data=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp_data,file="UniSigExp.txt",sep="\t",row.names=F,quote=F)
rownames(outTab)=outTab$id

k<-read.csv("预后数据.csv",row.names=1)
x<-read.csv("gene.csv",row.names=1)
x<-as.matrix(x)
cvfit=cv.glmnet(x,Surv(k$OS.time,k$OS),nfold=10,family="cox")
plot(cvfit)
fit<-glmnet(x,Surv(k$OS.time,k$OS),family="cox")

coef.min=coef(cvfit,s="lambda.min")
active.min=which(coef.min!=0)
geneids<-colnames(x)[active.min]
index.min=coef.min[active.min]
combine<-cbind(geneids,index.min)
write.csv(combine,"gene_index.csv")

x_df<-as.data.frame(x)
m<-cbind(k,x_df)
res<-coxph(Surv(k$OS.time,k$OS)~A2MP1+AC007014.2+AC007570.1+AC010336.2+AC018781.1+AC020634.2+AC079753.1+AC079921.1+
             AC090587.2+AC090617.4+AC093367.1+AC102953.1+AC079921.1+AC090587.2+AC090617.4+AC093367.1+AC102953.1+AC104083.1+
             AC125807.1+AC135050.2+AC144831.1+AL034417.3+AL360181.4+AL365475.1+AL451007.1+AL589745.1
           +ALOX15P1+AP005432.2+BUB1B.PAK6+C2CD4D+CAMK2N1+CD177P1+CDC20P1+CETN4P+CORO2B+DNAJC12+EPHX3+FMN1+FOXB1+
             FUT4+GP5+GRAP+GRHL3.AS1+HNRNPA1P25+HOXB9+IFIT1B+IL10+IMPDH1P5+LINC00659+LINC01963+LINC02345+ME1+METTL7B+
             MIR4520.1+NBPF13P+PARP4P2+PCOLCE2+PTCHD4+RAB11FIP1+RNF225+RPSAP4+SEMA5A+SLC13A2+SLC23A1+SLC2A3+SLC5A9+SNORA14A+SPA17P1+
             +STARD9+TPRG1+TRAV36DV7+TRIML2+UPK2+USP6NL+VAT1L,data=m)

summary_res<-summary(res)
coefficients<-summary_res$coefficients
forest_data<-data.frame(
  Variable=rownames(coefficients),
  HR=coefficients[,"exp(coef)"],
  Lower_CI=coefficients[,"lower .95"],
  Upper_CI=coefficients[,"upper .95"],
  P_value=coefficients[,"Pr(>|z|)"]
)

ping<-as.data.frame(-2.554*m$AC010336.2+2.489*m$AC020634.2+7.433*m$AC079753.1-2.581*m$AC079921.1-7.021*m$AL034417.3
                    -1.084*m$ALOX15P1-18.51*m$BUB1B.PAK6+1.581*m$FMN1-1.143*m$FUT4-9.768*m$HNRNPA1P25+0.522*m$HOXB9+0.549*m$LINC02345
                    +0.311*m$PCOLCE2+0.264*m$SLC2A3-4.898*m$SLC5A9+0.837*m$SNORA14A+0.427*m$UPK2+1.496*m$VAT1L)
rownames(ping)<-rownames(m)
colnames(ping)<-"score"

f<-cbind(k,ping)
f<-f[order(f$score),]
list<-c(rep("low",189),rep("high",188))%>%factor(.,levels=c("low","high"),ordered=F)
list<-model.matrix(~factor(list)+0)
colnames(list)<-c("low","high")
list<-as.data.frame(list)
f<-cbind(f,list$high)
colnames(f)[4]<-"fen"
fit<-survfit(Surv(OS.time,OS)~fen,data=f)
ggsurvplot(fit,data=f,conf.int=F,risk.table=TRUE,surv.median.line="hv",pval=TRUE)
ggsurvplot(fit,data=f,conf.int=TRUE,pval=TRUE,add.all=F)

write.table(diffexp,file="diffexp.txt",sep="\t",quote=F,row.names=T)
write.table(f,file="fen.txt",sep="\t",quote=F,row.names=T)
write.table(k,file="os.txt",sep="\t",quote=F,row.names=T)
write.table(m,file="os_exp.txt",sep="\t",quote=F,row.names=T)
write.table(ping,file="score.txt",sep="\t",quote=F,row.names=T)

rt=data.frame(uniSigExp[,1:2],uniSigExp[,lassoGene])
cox<-coxph(Surv(futime,fustat)~.,data=rt)
cox=step(cox,direction="both")
riskScore=predict(cox,type="risk",newdata=rt)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,1:2],riskScore,risk)),cbind(rt[,1:2],riskScore,risk)),file="risk.txt",sep="\t",quote=F,row.names=F)

coxSummary=summary(cox)
outTab=cbind(HR=coxSummary$conf.int[,"exp(coef)"],
             HR.95L=coxSummary$conf.int[,"lower .95"],
             HR.95H=coxSummary$conf.int[,"upper .95"],
             pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
risk=data.frame(rt[,1:2],rt[,rownames(outTab)],riskScore,risk)
rt=as.data.frame(outTab)
gene<-rownames(rt)
hr<-sprintf("%.3f",as.numeric(rt$"HR"))
hrLow<-sprintf("%.3f",as.numeric(rt$"HR.95L"))
hrHigh<-sprintf("%.3f",as.numeric(rt$"HR.95H"))
Hazard.ratio<-paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal<-ifelse(as.numeric(rt$pvalue)<0.001,"<0.001",sprintf("%.3f",as.numeric(rt$pvalue)))

pdf(file="mulforest2.pdf",width=6,height=6)
n<-nrow(rt)
nRow<-n+1
ylim<-c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))
xlim=c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim=c((min(as.numeric(hrLow))-0.1),(max(as.numeric(hrHigh))+0.1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor=ifelse(as.numeric(hr)>1,'red','green')
points(as.numeric(hr),n:1,pch=15,col=boxcolor,cex=1.3)
axis(1)
dev.off()



surv <- read.table("surival_data.txt", sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
surv$group <- ifelse(surv$RiskScore > median(surv$RiskScore), "High", "Low")
surv$group <- factor(surv$group, levels = c("Low", "High"))
table(surv$group)

library(survival)
library(survminer)

fitd <- survdiff(Surv(OS.time, OS) ~ group, data = surv, na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

fit <- survfit(Surv(OS.time, OS) ~ group, data = surv)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ", round(pValue, 3))))

ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           palette = "aaas",
           legend.labs = c("Low", "High"),
           size = 1,
           xlim = c(0, 20),
           break.time.by = 1,
           legend.title = "RiskScore",
           surv.median.line = "hv",
           ylab = "Survival probability (%)",
           xlab = "Time (Years)",
           ncensor.plot = TRUE,
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)



library(ggplot2)
library(survival)
library(survminer)

coxSummary <- summary(cox)
outTab <- cbind(HR = coxSummary$conf.int[, "exp(coef)"],
                HR.95L = coxSummary$conf.int[, "lower .95"],
                HR.95H = coxSummary$conf.int[, "upper .95"],
                pvalue = coxSummary$coefficients[, "Pr(>|z|)"])

riskScore <- predict(cox, type = "risk", newdata = surv.dat)
surv.dat$RiskScore <- riskScore
surv.dat$RiskGroup <- ifelse(surv.dat$RiskScore > median(surv.dat$RiskScore), "High", "Low")

write.table(surv.dat, "surv_dat_with_risk_score.txt", row.names = TRUE, quote = FALSE, sep = "\t")

fit <- survfit(Surv(OS.time, OS) ~ RiskGroup, data = surv.dat)
ggsurvplot(fit, data = surv.dat, pval = TRUE, risk.table = TRUE)

rt <- as.data.frame(outTab)
gene <- rownames(rt)
hr <- sprintf("%.3f", as.numeric(rt$HR))
hrLow <- sprintf("%.3f", as.numeric(rt$HR.95L))
hrHigh <- sprintf("%.3f", as.numeric(rt$HR.95H))
Hazard.ratio <- paste0(hr, " (", hrLow, "-", hrHigh, ")")
pVal <- ifelse(as.numeric(rt$pvalue) < 0.001, "<0.001", sprintf("%.3f", as.numeric(rt$pvalue)))

write.table(rt, "多因素森林图数据.txt", row.names = TRUE, quote = FALSE, sep = "\t")

pdf(file = "mulforest2.pdf", width = 6, height = 6)
n <- nrow(rt)
nRow <- n + 1
ylim <- c(1, nRow)
layout(matrix(c(1, 2), nc = 2), width = c(3, 2))

par(mar = c(4, 2.5, 2, 1))
plot(1, xlim = c(0, 3), ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
text.cex <- 0.8
text(0, n:1, gene, adj = 0, cex = text.cex)
text(1.5 - 0.5 * 0.2, n:1, pVal, adj = 1, cex = text.cex)
text(1.5 - 0.5 * 0.2, n + 1, 'pvalue', cex = text.cex, font = 2, adj = 1)
text(3, n:1, Hazard.ratio, adj = 1, cex = text.cex)
text(3, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)

par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
xlim <- c(min(as.numeric(hrLow)) - 0.1, max(as.numeric(hrHigh)) + 0.1)
plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, ylab = "", xaxs = "i", xlab = "Hazard ratio")
arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
abline(v = 1, col = "black", lty = 2, lwd = 2)
boxcolor <- ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.3)
axis(1)
dev.off()



```r
risk <- read.csv("患者分布可视化数据.csv", row.names = 1, check.names = FALSE)
rt <- risk
rt <- rt[order(rt$riskScore),]

riskClass <- rt[, "risk"]
lowLength <- length(riskClass[riskClass == "low"])
highLength <- length(riskClass[riskClass == "high"])
line <- rt[, "riskScore"]
line[line > 10] <- 10

pdf(file = "riskScore11.pdf", width = 10, height = 3.5)
plot(line, type = "p", pch = 20,
     xlab = "Patients (increasing risk socre)", ylab = "Risk score",
     col = c(rep("green", lowLength), rep("red", highLength)))
abline(h = median(rt$riskScore), v = lowLength, lty = 2)
legend("topleft", c("High risk", "Low risk"), bty = "n", pch = 19, col = c("red", "green"), cex = 1.2)
dev.off()

color <- as.vector(rt$fustat)
color[color == 1] <- "red"
color[color == 0] <- "green"
pdf(file = "survstat.pdf", width = 10, height = 3.5)
plot(rt$futime, pch = 19,
     xlab = "Patients (increasing risk socre)", ylab = "Survival time (years)",
     col = color)
legend("topleft", c("Dead", "Alive"), bty = "n", pch = 19, col = c("red", "green"), cex = 1.2)
abline(v = lowLength, lty = 2)
dev.off()

rt1 <- rt[c(3:(ncol(rt) - 2))]
rt1 <- log2(rt1 + 1)
rt1 <- t(rt1)
annotation <- data.frame(type = rt[, ncol(rt)])
rownames(annotation) <- rownames(rt)

library(pheatmap)
pdf(file = "heatmap.pdf", width = 10, height = 10)
pheatmap(rt1,
         annotation = annotation,
         cluster_cols = FALSE,
         fontsize_row = 11,
         show_colnames = F,
         fontsize_col = 3,
         scale = "row",
         color = colorRampPalette(c("green", "black", "red"))(100))
dev.off()


library(pROC)
library(ggplot2)
library(survival)
library(survminer)
library(timeROC)

risk <- read.table("risk.txt", header = T, sep = "\t", check.names = F, row.names = 1)
risk <- risk[, c("time", "state", "riskScore")]

cli <- read.table("clinical.txt", header = T, sep = "\t", check.names = F, row.names = 1)

samSample <- intersect(row.names(risk), row.names(cli))
risk1 <- risk[samSample, , drop = F]
cli <- cli[samSample, , drop = F]
rt <- cbind(risk1, cli)

bioCol <- c("DarkOrchid", "Orange2", "MediumSeaGreen", "NavyBlue", "#8B668B", "#FF4500", "#135612", "#561214")

roc1 <- roc(rt$state ~ rt$riskScore)
pdf(file = "ROC.riskscore.pdf", width = 5, height = 5)
plot(roc1, print.auc = TRUE, col = bioCol, legacy.axes = T)
dev.off()

ROC_rt <- timeROC(T = risk$time, delta = risk$state,
                  marker = risk$riskScore, cause = 1,
                  weighting = 'aalen',
                  times = c(1, 3, 5), ROC = TRUE)
pdf(file = "ROC.all.pdf", width = 5, height = 5)
plot(ROC_rt, time = 1, col = bioCol[1], title = FALSE, lwd = 4)
plot(ROC_rt, time = 3, col = bioCol[2], add = TRUE, title = FALSE, lwd = 4)
plot(ROC_rt, time = 5, col = bioCol[3], add = TRUE, title = FALSE, lwd = 4)
legend('bottomright',
       c(paste0('AUC at 1 years: ', sprintf("%.03f", ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ', sprintf("%.03f", ROC_rt$AUC[3]))),
       col = bioCol[1:3], lwd = 4, bty = 'n', title = "All set")
dev.off()

predictTime <- 5
aucText <- c()
pdf(file = "cliROC.all.pdf", width = 5.5, height = 5.5)

ROC_rt <- timeROC(T = risk$time, delta = risk$state,
                  marker = risk$riskScore, cause = 1,
                  weighting = 'aalen',
                  times = c(predictTime), ROC = TRUE)
plot(ROC_rt, time = predictTime, col = bioCol[1], title = FALSE, lwd = 4)
aucText <- c(paste0("Risk", ", AUC=", sprintf("%.3f", ROC_rt$AUC[2])))
abline(0, 1)

for (i in 4:ncol(rt)) {
  ROC_rt <- timeROC(T = rt$time, delta = rt$state,
                    marker = rt[, i], cause = 1,
                    weighting = 'aalen',
                    times = c(predictTime), ROC = TRUE)
  plot(ROC_rt, time = predictTime, col = bioCol[i - 2], title = FALSE, lwd = 4, add = TRUE)
  aucText <- c(aucText, paste0(colnames(rt)[i], ", AUC=", sprintf("%.3f", ROC_rt$AUC[2])))
}

legend("bottomright", aucText, lwd = 4, bty = "n", col = bioCol[1:(ncol(rt) - 1)], title = "All set")
dev.off()


library(limma)
library(ggplot2)
library(scatterplot3d)
library(pheatmap)

rt <- read.table("risk.txt", header = T, sep = "\t", check.names = F, row.names = 1)
risk <- as.vector(rt$risk)
data <- rt[, 3:(ncol(rt) - 2)]

data.pca <- prcomp(data, scale. = TRUE)
pcaPredict <- predict(data.pca)

PCA <- data.frame(PC1 = pcaPredict[, 1], PC2 = pcaPredict[, 2], risk = risk)
PCA.mean <- aggregate(PCA[, 1:2], list(risk = PCA$risk), mean)

pdf(file = "PCA.2d.pdf", width = 5.5, height = 4.75)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk, shape = risk)) +
  scale_colour_manual(name = "risk", values = c("DarkOrchid", "Orange2")) +
  theme_bw() +
  labs(title = "PCA") +
  theme(plot.margin = unit(rep(1.5, 4), 'lines'), plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = PCA.mean$PC1, y = PCA.mean$PC2, label = PCA.mean$risk, cex = 7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

color <- ifelse(risk == "high", "DarkOrchid", "Orange2")
pdf(file = "PCA.3d.pdf", width = 7, height = 7)
par(oma = c(1, 1, 2.5, 1))
s3d <- scatterplot3d(pcaPredict[, 1:3], pch = 16, color = color, angle = 60)
legend("top", legend = c("High risk", "Low risk"), pch = 16, inset = -0.2, 
       box.col = "white", xpd = TRUE, horiz = TRUE, col = c("DarkOrchid", "Orange2"))
dev.off()

rt <- rt[order(rt$riskScore), ]
riskClass <- rt[, "risk"]
lowLength <- length(riskClass[riskClass == "low"])
highLength <- length(riskClass[riskClass == "high"])
lowMax <- max(rt$riskScore[riskClass == "low"])
line <- rt[, "riskScore"]
line[line > 10] <- 10

pdf(file = "riskScore.pdf", width = 7, height = 4)
plot(line, type = "p", pch = 20,
     xlab = "Patients (increasing risk socre)",
     ylab = "Risk score",
     col = c(rep("MediumSeaGreen", lowLength), rep("Firebrick3", highLength)))
abline(h = lowMax, v = lowLength, lty = 2)
legend("topleft", c("High risk", "Low risk"), bty = "n", pch = 19, 
       col = c("Firebrick3", "MediumSeaGreen"), cex = 1.2)
dev.off()

color <- as.vector(rt$state)
color[color == 1] = "Firebrick3"
color[color == 0] = "MediumSeaGreen"
pdf(file = "survStat.pdf", width = 7, height = 4)
plot(rt$time, pch = 19,
     xlab = "Patients (increasing risk socre)",
     ylab = "Survival time (years)",
     col = color)
legend("topleft", c("Dead", "Alive"), bty = "n", pch = 19, 
       col = c("Firebrick3", "MediumSeaGreen"), cex = 1.2)
abline(v = lowLength, lty = 2)
dev.off()

ann_colors <- list()
bioCol <- c("MediumSeaGreen", "Firebrick3")
names(bioCol) <- c("low", "high")
ann_colors[["Risk"]] <- bioCol

rt1 <- rt[c(3:(ncol(rt) - 2))]
rt1 <- t(rt1)
annotation <- data.frame(Risk = rt[, ncol(rt)])
rownames(annotation) <- rownames(rt)

pdf(file = "heatmap.pdf", width = 6, height = 5)
pheatmap(rt1,
         annotation = annotation,
         annotation_colors = ann_colors,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         show_colnames = F,
         scale = "row",
         color = colorRampPalette(c(rep("MediumSeaGreen", 3.5), "white", rep("Firebrick3", 3.5)))(50),
         fontsize_col = 3,
         fontsize = 7,
         fontsize_row = 8)
dev.off()