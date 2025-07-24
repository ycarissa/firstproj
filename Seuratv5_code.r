library(Seurat)
library (ggplot2)

#### Single Cell RNAseq Analysis based on the GEO database dataset GSE292807 ####
## Comparison of WT to NOTCH2 gain-of-function mice ##


### Directories
work_dir <- "C:/Users/caris/Desktop/Projects/GSE292807"
save_dir <- paste0(work_dir, "/Seuratv5_Analysis")
dir.create(save_dir)


### Create Seurat Object
counts_het <- Read10X(data.dir = paste0(work_dir, "/Het_NOTCH"))
counts_wt  <- Read10X(data.dir = paste0(work_dir, "/WT_NOTCH"))

Het.s <- CreateSeuratObject(counts = counts_het, min.cells = 3, min.features = 200)
WT.s  <- CreateSeuratObject(counts = counts_wt,  min.cells = 3, min.features = 200)


### QC filtering
qc_dir <- paste0(save_dir, "/QC_filtering")
dir.create(qc_dir)

Het.s[["percent.mt"]] <- PercentageFeatureSet(Het.s, pattern = "^MT-")
WT.s[["percent.mt"]] <- PercentageFeatureSet(WT.s, pattern = "^MT-")

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

### Subset based on Cells with >250 nFeature, >10% percent mt and >10,000 genes/cell based on original analysis
Het_subset.s <- subset(Het.s, subset = nFeature_RNA > 250 & nFeature_RNA < 10000 & nCount_RNA < 10000 & percent.mt < 5)
WT_subset.s <- subset(WT.s, subset = nFeature_RNA > 250 & nFeature_RNA < 10000 & nCount_RNA < 10000 & percent.mt < 5)




