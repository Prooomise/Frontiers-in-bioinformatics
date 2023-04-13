# Read in the data and packages
library(Seurat)
library(harmony)
library(MAST)
library(ggplot2)
library(ggrepel)
library(org.Mm.eg.db)
library(clusterProfiler)
scRNA1 <- readRDS("./SRR23407930.rds")
scRNA2 <- readRDS("./SRR23407934.rds")

# Construct the dataset
combined <- merge(x = scRNA1, y = c(scRNA2))
Group <- c(rep("Treatment", dim(scRNA1)[2]), rep("CK", dim(scRNA2)[2]))
combined <- AddMetaData(combined, Group, col.name = "split")
DefaultAssay(combined) <- "RNA"

# Normalize the data
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
s.genes <- capwords(tolower(cc.genes$s.genes))
g2m.genes <- capwords(tolower(cc.genes$g2m.genes))
combined <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
combined <- RunHarmony(combined, dims = 1:8, "split",verbose = FALSE)
combined <- RunUMAP(combined,reduction = "harmony", dims = 1:30, verbose = FALSE)
DimPlot(combined, reduction = "umap", group.by = "split")
DimPlot(combined, reduction = "umap", group.by = "Phase")

##ggsave("./UMAP.png", height = 6, width = 8, dpi = 600)
rm(scRNA1, scRNA2)
##saveRDS(combined, "./combined.rds")

# Clustering the cells and find markers
combined <- FindNeighbors(combined, dims = 1:20, reduction = "harmony")
combined <- FindClusters(combined, resolution = 0.1)
DimPlot(combined, reduction = "umap", label=TRUE)
diff_gene <- FindMarkers(combined, ident.1 = "Treatment", ident.2 = "CK", logfc.threshold = 0,min.pct = 0, 
                         test.use="MAST", group.by="split")

# Volcano plot
degdf <- subset(diff_gene, p_val_adj<0.05 & abs(avg_log2FC)>0.15)
degdf$sign <- ifelse(degdf$p_val_adj < 0.0005 & abs(degdf$avg_log2FC) > 1,
                     rownames(degdf), NA)
P.Value_t = 0.01
degdf$change = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < -0.5,"down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 0.5,"up",
                             "stable"))
ggplot(degdf, aes(avg_log2FC,  -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=3, aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.3) +
  geom_vline(xintercept = 0.5,lty=2,col="black",lwd=0.3)+
  geom_vline(xintercept = -0.5,lty=2,col="black",lwd=0.3)+
  geom_label_repel(aes(label=sign), fontface="bold", color="grey50", box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), segment.colour = "grey50")+
  theme_bw()
#ggsave("./scRNA-volcano.png", width = 8, height = 6, dpi = 600)

# GO analysis
sig_gene <- subset(diff_gene, p_val_adj<0.05 & abs(avg_log2FC)>0.15)
gene_up=rownames(sig_gene[sig_gene$avg_log2FC > 0.15,])
gene_down=rownames(sig_gene[sig_gene$avg_log2FC < -0.15,])

gene_up=as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                   keys = gene_up,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))
gene_down=as.character(na.omit(AnnotationDbi::select(org.Mm.eg.db,
                                                     keys = gene_down,
                                                     columns = 'ENTREZID',
                                                     keytype = 'SYMBOL')[,2]))
##GO-up
geneset <- gene_up
deg_type <- "up"
ego_CC <- enrichGO(gene = geneset, OrgDb=org.Mm.eg.db, keyType = "ENTREZID",
                   ont = "CC", pAdjustMethod = "BH", minGSSize = 1
                   , pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
ego_BP <- enrichGO(gene = geneset, OrgDb=org.Mm.eg.db, keyType = "ENTREZID",
                   ont = "BP", pAdjustMethod = "BH", minGSSize = 1,
                   pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
ego_MF <- enrichGO(gene = geneset, OrgDb=org.Mm.eg.db, keyType = "ENTREZID",
                   ont = "MF", pAdjustMethod = "BH", minGSSize = 1,
                   pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

display_number = c(10,10,10)
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID), Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:10],collapse = " ") 
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}
go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = paste0("The Most Enriched GO Terms ",deg_type))+
  theme_bw()
## ggsave("./GO-up.png", width = 12, height = 6, dpi = 600)

##GO-down
geneset <- gene_down
deg_type <- "down"
ego_CC <- enrichGO(gene = geneset, OrgDb=org.Mm.eg.db, keyType = "ENTREZID",
                   ont = "CC", pAdjustMethod = "BH", minGSSize = 1,
                   pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
ego_BP <- enrichGO(gene = geneset, OrgDb=org.Mm.eg.db, keyType = "ENTREZID",
                   ont = "BP", pAdjustMethod = "BH", minGSSize = 1,
                   pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
ego_MF <- enrichGO(gene = geneset, OrgDb=org.Mm.eg.db, keyType = "ENTREZID",
                   ont = "MF", pAdjustMethod = "BH", minGSSize = 1,
                   pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)
display_number = c(10,10,10)
ego_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]

go_enrich_df <- data.frame(
  ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID), Description=c(ego_result_BP$Description,ego_result_CC$Description,ego_result_MF$Description),
  GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
  type=factor(c(rep("biological process", display_number[1]), 
                rep("cellular component", display_number[2]),
                rep("molecular function", display_number[3])), 
              levels=c("biological process", "cellular component","molecular function" )))

for(i in 1:nrow(go_enrich_df)){
  description_splite=strsplit(go_enrich_df$Description[i],split = " ")
  description_collapse=paste(description_splite[[1]][1:10],collapse = " ") 
  go_enrich_df$Description[i]=description_collapse
  go_enrich_df$Description=gsub(pattern = "NA","",go_enrich_df$Description)
}

go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + 
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = COLS) + 
  coord_flip() + 
  xlab("GO term") + 
  ylab("Gene_Number") + 
  labs(title = paste0("The Most Enriched GO Terms ",deg_type))+
  theme_bw()
## ggsave("./GO-down.png", width = 12, height = 6, dpi = 600)

# KEGG -> under construction, need to be continued
kk <- enrichKEGG(gene = geneset,keyType = "kegg",organism= "mmn", 
                 qvalueCutoff = 0.05, pvalueCutoff=0.01)
hh <- as.data.frame(kk)
rownames(hh) <- 1:nrow(hh)
hh$order=factor(rev(as.integer(rownames(hh))),labels = rev(hh$Description))
ggplot(hh,aes(y=order,x=Count,fill=p.adjust))+
  geom_point(stat = "identity",width=0.7)+
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = paste0("KEGG Pathways Enrichment ",deg_type),
       x = "Gene numbers", 
       y = "Pathways")+
  theme(axis.title.x = element_text(face = "bold",size = 8),
        axis.title.y = element_text(face = "bold",size = 8),
        legend.title = element_text(face = "bold",size = 8),
        axis.text.y = element_text(size=3))
