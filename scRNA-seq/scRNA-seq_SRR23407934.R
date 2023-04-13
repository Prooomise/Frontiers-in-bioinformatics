# Read in the packages and data
library(dplyr)
library(Seurat)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
rawdata <- Read10X((data.dir = "SRR23407934"), gene.column = 2)
scObject <- CreateSeuratObject(rawdata)

# Pre-processing - Transform gene name
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
s.genes <- capwords(tolower(cc.genes$s.genes))
g2m.genes <- capwords(tolower(cc.genes$g2m.genes))
scObject <- CellCycleScoring(scObject, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)


# Pre-processing - wash the data
library(patchwork)
scObject[["percent.mt"]] <- PercentageFeatureSet(scObject, pattern = "^mt-")
VlnPlot(scObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.2)
## ggsave("./pre-process-violin-before.png", height = 4, width = 6, dpi = 600)
scObject <- subset(scObject, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 20000)
plot1 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
## ggsave("./pre-process-violin-after.png", height = 4, width = 6, dpi = 600)
CombinePlots(plots = list(plot1, plot2))
## ggsave("./pre-process-point.png", height = 4, width = 8, dpi = 600)
rm(plot1, plot2)

# Pre-processing 
scObject <- NormalizeData(scObject, normalization.method = "RC", scale.factor = 10000)
scObject <- FindVariableFeatures(scObject, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
Top5 <- head(VariableFeatures(scObject), 5)
p1 <- VariableFeaturePlot(scObject)
LabelPoints(plot = p1, points = Top5, repel = TRUE)  
## ggsave("./VariableGene.png", height = 6, width = 9, dpi = 600)
rm(p1)

scObject <- ScaleData(scObject, features = rownames(scObject))
scObject <- SCTransform(scObject, verbose = TRUE, vars.to.regress = c("percent.mt","nCount_RNA"), variable.features.n = 3000)
scObject <- RunPCA(scObject, verbose = TRUE, features = VariableFeatures(scObject))
ElbowPlot(scObject, ndims = 50)

scObject <- RunUMAP(scObject, dims = 1:8, n.neighbors = 30L, min.dist = 0.01, verbose = FALSE)
DimPlot(scObject, group.by = "orig.ident", reduction = "umap") +NoLegend()
## ggsave("./umap1.png", width = 8, height = 6, dpi = 600)
DimPlot(scObject, group.by = "Phase", reduction = "umap")
## saveRDS(scObject, "SRR23407934.rds")

# Clustering
scObject <- FindNeighbors(scObject, dims = 1:8, verbose = TRUE)
scObject <- FindClusters(scObject, verbose = TRUE,resolution = 0.1)
DimPlot(scObject, group.by = "seurat_clusters",reduction = "umap")
## ggsave("./Cluster.png", width = 8, height = 6, dpi = 600)

#Violin plots for phase marker (Mki67) and neuronal precursor marker (Nestin)
VlnPlot(scObject, features = "Mki67",log = TRUE, pt.size = 0, assay = "RNA") + 
  theme(plot.title = element_text(hjust = 0.5))
VlnPlot(scObject, features = "Nes", log = TRUE, pt.size = 0, assay = "RNA") + 
  theme(plot.title = element_text(hjust = 0.5))
#Violin plots for phase marker (Mki67) and neuronal precursor marker (Nestin) grouped by Phase with boxplots
VlnPlot(scObject, features = "Mki67",log = TRUE,pt.size = 0,group.by = "Phase",assay = "RNA") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_boxplot(width=0.075,fill="grey")
VlnPlot(scObject, features = "Nes",log = TRUE,pt.size = 0,group.by = "Phase",assay = "RNA") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_boxplot(width=0.075,fill="grey")

# Find marker
S.marker <- FindMarkers(scObject, assay = "RNA", group.by = "Phase", ident.1 = "S", min.diff.pct = 0.2, verbose = FALSE)
G1.marker <- FindMarkers(scObject, assay = "RNA", group.by = "Phase", ident.1 = "G1", min.diff.pct = 0.2, verbose = FALSE)
G2M.marker <- FindMarkers(scObject, assay = "RNA", group.by = "Phase", ident.1 = "G2M", min.diff.pct = 0.2, verbose = FALSE)
markers <- c("Celf5","Id2","Hmga2","Irx3","Scn1b","Hells", "Pola1", "Mki67")
DotPlot(scObject, features = markers, group.by = "Phase")
DotPlot(scObject, features = markers, group.by = "seurat_clusters")


