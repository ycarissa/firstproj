library(Seurat)
library(ggplot2)
library(glmGamPoi)
library(dplyr)

#### Single Cell RNAseq Analysis based on the GEO database dataset GSE282630 ####
## Comparison of WT to NOTCH2 gain-of-function mice ##


### Directories
work_dir <- "C:/Users/caris/Desktop/Projects/GSE292807"
save_dir <- paste0(work_dir, "/Seuratv5_Analysis")
dir.create(save_dir)


### Create Seurat Object
counts_het <- Read10X(data.dir = paste0(work_dir, "/Het_NOTCH"))
counts_wt  <- Read10X(data.dir = paste0(work_dir, "/WT_NOTCH"))

Het.s <- CreateSeuratObject(counts = counts_het, project = "Het", min.cells = 3, min.features = 200)
WT.s  <- CreateSeuratObject(counts = counts_wt,  project = "WT", min.cells = 3, min.features = 200)


### QC filtering
qc_dir <- paste0(save_dir, "/QC_filtering")
dir.create(qc_dir)

Het.s[["percent.mt"]] <- PercentageFeatureSet(Het.s, pattern = "^mt-")
WT.s[["percent.mt"]] <- PercentageFeatureSet(WT.s, pattern = "^mt-")

#Violin plots for filtering for Seurat Object Het.s
p1 <- VlnPlot(Het.s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + labs(x = "Identity", y = "Expression Level")
p2 <- FeatureScatter(Het.s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5) 
p3 <- FeatureScatter(Het.s, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5) 
p4 <- p2 + p3

ggsave("QC_violin_Het.png", plot = p1, path = qc_dir, width = 10, height = 8, dpi = 600)
ggsave("QC_scatter_Het.png", plot = p4, path = qc_dir, width = 10, height = 8, dpi = 600)


#Violin plots for filtering for Seurat Object WT.s
p1 <- VlnPlot(WT.s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + labs(x = "Identity", y = "Expression Level")
p2 <- FeatureScatter(WT.s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5) 
p3 <- FeatureScatter(WT.s, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5) 
p4 <- p2 + p3

ggsave("QC_violin_WT.png", plot = p1, path = qc_dir, width = 10, height = 8, dpi = 600)
ggsave("QC_scatter_WT.png", plot = p4, path = qc_dir, width = 10, height = 8, dpi = 600)

### Subset based on Cells with 15000 > nFeature > 250, >10% percent mt and 15000 > nCount.
Het_subset.s <- subset(Het.s, subset = nFeature_RNA > 250 & nFeature_RNA < 10000 & nCount_RNA < 15000 & percent.mt < 10)
WT_subset.s <- subset(WT.s, subset = nFeature_RNA > 250 & nFeature_RNA < 10000 & nCount_RNA < 15000 & percent.mt < 10)

### Integration
combined.s <- merge(x = Het_subset.s, y = WT_subset.s, project = "combined", add.cell.ids = c("Het", "WT"))
Idents(combined.s) <- "orig.ident"

combined.s <- SCTransform(combined.s, vars.to.regress = "percent.mt", verbose = FALSE)
combined.s <- RunPCA(combined.s, assay = "SCT", verbose = FALSE)
ElbowPlot(combined.s) ### 1:10 was chosen


combined.s <- IntegrateLayers(object = combined.s, method = CCAIntegration, assay = "SCT", normalization.method = "SCT", orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

combined.s <- RunUMAP(combined.s, reduction = "integrated.cca", dims = 1:10)
combined.s <- FindNeighbors(combined.s, reduction = "integrated.cca", dims = 1:10)
combined.s <- FindClusters(combined.s, resolution = 0.5)
DimPlot(combined.s, reduction = "umap")

###Saving plots
u1 <- DimPlot(combined.s, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
  labs(title = "UMAP by Clusters")
u2 <- DimPlot(combined.s, reduction = "umap", label = TRUE, group.by = "orig.ident") +
  labs(title = "UMAP by Identity")

ggsave(filename = file.path(save_dir, "UMAP_by_Clusters.png"), plot = u1, width = 8, height = 6, dpi=600)
ggsave(filename = file.path(save_dir, "UMAP_by_Identity.png"), plot = u2, width = 8, height = 6, dpi=600)

