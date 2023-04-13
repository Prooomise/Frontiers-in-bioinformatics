library(DESeq2)
library(ggplot2)
library(gplots)
library(ggrepel)
library(org.Mm.eg.db)
library(clusterProfiler)
library(fgsea)
library(dplyr)
library(msigdbr)

#build suitable expression matrix
countData5824 <- read.table("SRR23405824countResult.txt", sep = "\t", header = TRUE, row.names=1)
countData5829 <- read.table("SRR23405829countResult.txt", sep = "\t", header = TRUE, row.names=1)
count5824 <- countData5824[c('SRR23405824.bam')]
count5829 <- countData5829[c('SRR23405829.bam')]
rm(countData5824)
rm(countData5829)
colnames(count5824) <- c("CK")
colnames(count5829) <- c("Treatment")

Data <- cbind(count5824, count5824, count5829, count5829)
colnames(Data) <- c("CK_1", "CK_2", "Treatment_1", "Treatment_2") 
condition <- factor(c("CK", "CK", "Treatment", "Treatment"))
coldata <- data.frame(row.names=colnames(Data), condition)
rm(count5824)
rm(count5829)

#build dds matrix
dds <- DESeq2::DESeqDataSetFromMatrix(countData=Data, 
                                      colData=coldata, design= ~condition)

#filter low quality data
dds <- dds[rowSums(BiocGenerics::counts(dds)) > 100, ]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds) 
dispersions(dds) <- mcols(dds)$dispGeneEst 
dds <- nbinomWaldTest(dds)
resdata <- data.frame(results(dds, lfcThreshold=1, alpha=0.05)) 
resdata <- tibble::rownames_to_column(resdata, "ENSEMBL")
#write.csv(resdata,file= "DESeq2_logFC.csv")

#draw volcano plot
resdata$threshold <- as.factor(ifelse(resdata$padj < 0.001 & abs(resdata$log2FoldChange) >= 2 ,
                              ifelse(resdata$log2FoldChange >= 2 ,
                                     'Up','Down'),'None')) 
gene.df <- bitr(resdata$ENSEMBL, "ENSEMBL", "SYMBOL", OrgDb = org.Mm.eg.db)
resdata <- right_join(resdata, gene.df, by="ENSEMBL", multiple = "all")
rm(gene.df)
resdata$label <- ""
resdata[resdata$padj < 0.001 & abs(resdata$log2FoldChange) >= 5,]$label <- resdata[resdata$padj < 0.001 & abs(resdata$log2FoldChange) >= 5,]$SYMBOL
resdata[resdata$log2FoldChange >= 5,]$label <- resdata[resdata$log2FoldChange >= 5,]$SYMBOL

ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
  ylab(bquote("-log"[10]~"(qvalue)")) +
  xlab(bquote("log"[2]~"(Fold Change)"))+
  geom_point()+
  geom_text_repel(
    aes(label = label),
    size = 3,
    segment.color = "black", show.legend = FALSE )+
  ylim(0,50) + xlim(-12,12) +
  scale_color_manual(values=c("#2F5688", "#BBBBBB", "#CC0000"))+
  theme_minimal()+
  geom_hline(yintercept = -log10(0.001), linetype = "dashed")+
  geom_vline(xintercept = c(-2, 2), linetype = "dashed")+
  theme(legend.position = c(0.9, 0.8))+
  labs(colour = "")
#ggsave("volcano.png", dpi = 600, width = 8, height = 6)

subset(resdata,pvalue < 0.001) -> diff
subset(diff,log2FoldChange < -2) -> down
subset(diff,log2FoldChange > 2) -> up

down_geo <- clusterProfiler::enrichGO(gene = down$ENSEMBL, keyType = "ENSEMBL",
                                      OrgDb = BiocGenerics::get("org.Mm.eg.db"), 
                                      ont = "BP",
                                      pAdjustMethod = "BH", qvalueCutoff = 0.05)

enrichplot::dotplot(down_geo, showCategory = 10)
#ggsave("Down_GO.png", dpi = 600, width = 8, height = 6)

up_geo <- clusterProfiler::enrichGO(gene = up$ENSEMBL, keyType = "ENSEMBL",
                                    OrgDb = BiocGenerics::get("org.Mm.eg.db"), 
                                    ont = "BP",
                                    pAdjustMethod = "BH", qvalueCutoff = 0.05)

enrichplot::dotplot(up_geo, showCategory =10)
#ggsave("Up_GO.png", dpi = 600, width = 8, height = 6)
rm(up, down, diff)

df_id <- bitr(resdata$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
df_all<-merge(resdata, df_id, by="SYMBOL", all=F)
df_all_sort <- df_all[order(df_all$log2FoldChange, decreasing = T),]
gene_fc <- df_all_sort$log2FoldChange
names(gene_fc) = df_id$SYMBOL
rm(df_id, df_all, df_all_sort)
geneSet <- msigdbr(species = "Mus musculus", category = "M8") #!!!
geneSet <- split(geneSet$gene_symbol, geneSet$gs_name)
KEGG <- fgsea(geneSet, gene_fc)

Up <- KEGG[ES > 0][head(order(pval), n=10), pathway]
Down <- KEGG[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(Up, rev(Down))
plotGseaTable(topPathways, gene_fc, KEGG, gseaParam = 0.5)
plotGseaTable()
