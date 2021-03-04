#Si quieres sacar anotación de arrays Illumina
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)

annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#Código sin tener en cuenta universo de genes.

testgenes<-unique(as.character(dfleca$gene))
testgenes<-testgenes[!testgenes==""]
eg = bitr(testgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
pdf("gene_enrichment_dfleca.pdf")
kk <- enrichKEGG(gene= eg$ENTREZID,organism= 'hsa',pvalueCutoff = 0.05)
dotplot(kk,font.size = 8,title = "KEGG")
egomf <- enrichGO(gene= eg$ENTREZID,OrgDb= org.Hs.eg.db,ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable= TRUE)
dotplot(egomf,font.size = 8,title = "GO-MF")
egobp <- enrichGO(gene= eg$ENTREZID,OrgDb= org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable= TRUE)
dotplot(egobp,font.size = 8,title = "GO-BP")
egocc <- enrichGO(gene= eg$ENTREZID,OrgDb= org.Hs.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable= TRUE)
dotplot(egocc,font.size = 8,title = "GO-CC")
react<-enrichPathway(eg$ENTREZID)
dotplot(react,font.size = 8,title = "REACTOME")
gmtfile <- "/illumina/databases/references/gsea_db/c5.bp.v6.2.entrez.gmt"
c5 <- read.gmt(gmtfile)
egmt <- enricher(eg$ENTREZID, TERM2GENE=c5)
dotplot(egmt,font.size = 8,title = "MSigDB GO-BP")
write.table(as.data.frame(egmt),"msigdb_bp.csv",col.names=T,row.names=F,sep="\t",quote=F)
gmtfile <- "/illumina/databases/references/gsea_db/c5.mf.v6.2.entrez.gmt"
c5 <- read.gmt(gmtfile)
egmt <- enricher(eg$ENTREZID, TERM2GENE=c5)
dotplot(egmt,font.size = 8,title = "MSigDB GO-MF")
write.table(as.data.frame(egmt),"msigdb_mf.csv",col.names=T,row.names=F,sep="\t",quote=F)
gmtfile <- "/illumina/databases/references/gsea_db/c5.cc.v6.2.entrez.gmt"
c5 <- read.gmt(gmtfile)
egmt <- enricher(eg$ENTREZID, TERM2GENE=c5)
dotplot(egmt,font.size = 8,title = "MSigDB GO-CC")
write.table(as.data.frame(egmt),"msigdb_cc.csv",col.names=T,row.names=F,sep="\t",quote=F)
gmtfile <- "/illumina/databases/references/gsea_db/c2.cp.kegg.v6.2.entrez.gmt"
c5 <- read.gmt(gmtfile)
egmt <- enricher(eg$ENTREZID, TERM2GENE=c5)
dotplot(egmt,font.size = 8,title="MSigDB KEGG")
write.table(as.data.frame(egmt),"msigdb_kegg.csv",col.names=T,row.names=F,sep="\t",quote=F)
dev.off()

write.table(as.data.frame(kk),"kegg.csv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(as.data.frame(egomf),"GO_MF.csv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(as.data.frame(egobp),"GO_BP.csv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(as.data.frame(egocc),"GO_CC.csv",col.names=T,row.names=F,sep="\t",quote=F)
write.table(as.data.frame(react),"reactome.csv",col.names=T,row.names=F,sep="\t",quote=F)






#Código teniendo en cuenta universo de genes array:
  #FUNCTIONAL GENE ENRICHMENT
  library(clusterProfiler)
genelist<-fDataEPICunnested$UCSC_RefGene_Name
genelist<-unique(genelist[!genelist==""])
genelisteg<-bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
testgenes<-fDataEPICunnested[fDataEPICunnested$IlmnID %in% diffbetas_mia_all_subset$Row.names,"UCSC_RefGene_Name"]
testgenes<-unique(testgenes[!testgenes== ""])
testgenes_promoter<-unique(fDataEPICunnested[fDataEPICunnested$IlmnID %in% diffbetas_mia_all_subset$Row.names & (fDataEPICunnested$UCSC_RefGene_Group=="TSS1500" | fDataEPICunnested$UCSC_RefGene_Group=="TSS200" | fDataEPICunnested$UCSC_RefGene_Group=="5'UTR" | fDataEPICunnested$UCSC_RefGene_Group=="1stExon"),"UCSC_RefGene_Name"])

eg = bitr(testgenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

kk <- enrichKEGG(gene= eg$ENTREZID,organism= 'hsa',pvalueCutoff = 0.05)
pdf("kegg.pdf")
dotplot(kk)
dev.off()
write.table(as.data.frame(kk),"kegg.csv",col.names=T,row.names=F,sep="\t",quote=F)

egomf <- enrichGO(gene= eg$ENTREZID,universe= genelisteg$ENTREZID,OrgDb= org.Hs.eg.db,ont= "MF",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable= TRUE)
pdf("GO_MF.pdf")
dotplot(egomf)
dev.off()
write.table(as.data.frame(egomf),"GO_MF.csv",col.names=T,row.names=F,sep="\t",quote=F)

egobp <- enrichGO(gene= eg$ENTREZID,universe= genelisteg$ENTREZID,OrgDb= org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable= TRUE)
pdf("GO_BP.pdf")
dotplot(egobp)
dev.off()
write.table(as.data.frame(egobp),"GO_BP.csv",col.names=T,row.names=F,sep="\t",quote=F)

egocc <- enrichGO(gene= eg$ENTREZID,universe= genelisteg$ENTREZID,OrgDb= org.Hs.eg.db,ont= "CC",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable= TRUE)
pdf("GO_CC.pdf")
dotplot(egocc)
dev.off()
write.table(as.data.frame(egocc),"GO_CC.csv",col.names=T,row.names=F,sep="\t",quote=F)

react<-enrichPathway(eg$ENTREZID)
pdf("reactome.pdf")
dotplot(react)
dev.off()
write.table(as.data.frame(react),"reactome.csv",col.names=T,row.names=F,sep="\t",quote=F)

browseKEGG(kk,"hsa04010")
egmt <- enricher(eg$ENTREZID, TERM2GENE=c5)
mkk <- enrichMKEGG(gene =eg$ENTREZID ,organism = 'hsa')

#PLOT GSEA DEFINITIVO
gmtfile <- system.file("extdata", "c5.bp.v6.2.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
egmt <- enricher(eg$ENTREZID, TERM2GENE=c5)
png("kegg.png",width = 16, height = 8, units = 'in', res = 300)
barplot(egmt,drop=T,showCategory=12)
dev.off()