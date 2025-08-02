library(Seurat)
library(ggplot2)
library(glmGamPoi)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(AnnotationHub)
library(clusterProfiler)
library(org.Mm.eg.db)

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

ggsave(filename = file.path(qc_dir,"QC_violin_Het.png"), plot = p1, width = 10, height = 8, dpi = 600)
ggsave(filename = file.path(qc_dir,"QC_scatter_Het.png"), plot = p4, width = 10, height = 8, dpi = 600)


#Violin plots for filtering for Seurat Object WT.s
p1 <- VlnPlot(WT.s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) + labs(x = "Identity", y = "Expression Level")
p2 <- FeatureScatter(WT.s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5) 
p3 <- FeatureScatter(WT.s, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5) 
p4 <- p2 + p3

ggsave(filename = file.path(qc_dir,"QC_violin_WT.png"), plot = p1, width = 10, height = 8, dpi = 600)
ggsave(filename = file.path(qc_dir,"QC_scatter_WT.png"), plot = p4, width = 10, height = 8, dpi = 600)

### Subset based on Cells with 15000 > nFeature > 250, >10% percent mt and 15000 > nCount.
Het_subset.s <- subset(Het.s, subset = nFeature_RNA > 250 & nFeature_RNA < 10000 & nCount_RNA < 15000 & percent.mt < 10)
WT_subset.s <- subset(WT.s, subset = nFeature_RNA > 250 & nFeature_RNA < 10000 & nCount_RNA < 15000 & percent.mt < 10)

### Integration
combined.s <- merge(x = Het_subset.s, y = WT_subset.s, project = "combined", add.cell.ids = c("Het", "WT"))
Idents(combined.s) <- "orig.ident"

combined.s <- SCTransform(combined.s, vars.to.regress = "percent.mt", assay = "RNA")

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

### Cell Type Labeling
Idents(combined.s) <- "seurat_clusters"
combined.s <- PrepSCTFindMarkers(combined.s)

marker_dir <- paste0(save_dir, "/Marker_Genes")
dir.create(marker_dir)

clusters <- levels(combined.s)

### Loop to 1. find marker genes per cluster and 2. make vlnplots for the top 5 genes in each
for (cluster in clusters) {
  markers <- FindMarkers(combined.s, 
                                ident.1 = cluster,
                                only.pos = TRUE, 
                                min.pct = 0.25,            
                                logfc.threshold = 0.25,
                                assay = "SCT")
  
  write.csv(markers, file = file.path(marker_dir, "/cluster_", cluster, "_markers.csv"), row.names = TRUE)
  
  ###Violin plots of the top 5 genes of each cluster
  top5_genes <- head(rownames(markers), 5)
  vln_plots <- lapply(top5_genes, function(gene) {
    VlnPlot(combined.s, features = gene, pt.size = 0) +
      labs(title = gene, x = "Identity", y = "Expression Level") + theme(legend.position = "none")
  })
  
  all_5_plots <- wrap_plots(vln_plots, ncol = 1)
  ggsave(filename = file.path(marker_dir, "cluster_", cluster, "_VlnPlot.png"), plot = all_5_plots, width = 6, height = 10, dpi = 600)
}

### Cell type labeling was performed using PanglaoDB and literature as reference. ###

### Assigning cell types to seurat clusters
Idents(combined.s) <- "seurat_clusters"

cell_type_naming <- c(
  "0" = "Mature Chondrocytes",
  "1" = "Fibroblasts",
  "2" = "Limb Mesenchyme",
  "3" = "Fibroblasts", 
  "4" = "Early Chondrocytes",
  "5" = "Early Pregenitor Chondrocytes",
  "6" = "Myofibroblasts",
  "7" = "Fibrochondrocytes",
  "8" = "Osteogenic Chondrocytes",
  "9" = "Mitotic Chondrocytes",
  "10" = "Macrophages",
  "11" = "Hypertrophic Chondrocytes")

combined.s <- RenameIdents(combined.s, cell_type_naming)
combined.s$cell_type <- Idents(combined.s)

### Remaking plots from the original paper
#1 UMAP of cell types
u3 <- DimPlot(combined.s, reduction = "umap", label = TRUE, group.by = "cell_type") + ggtitle("UMAP with Labeled Cell Types")
ggsave(filename = file.path(save_dir, "UMAP_by_Cell_Type.png"), plot = u3, width = 8, height = 6, dpi=600)

#2 Bar graph displaying number of cells per cell type
cell_counts <- combined.s@meta.data %>%
  group_by(orig.ident, cell_type) %>%
  summarise(n = n(), .groups = "drop")

b1 <- ggplot(cell_counts, aes(x = orig.ident, y = n, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(x = "Condition", y = "Number of Cells", fill = "Cell Type") +
  theme_minimal()

ggsave(filename = file.path(save_dir, "stacked_celltype_barplot.png"), plot = b1, width = 6, height = 5, dpi = 600)

### Differential Expression Analysis  ###
Idents(combined.s) <- "cell_type"
cell_type <- levels(combined.s) 
de_dir <- paste0(save_dir, "/DGE_Analysis")
dir.create(de_dir)

combined.s <- PrepSCTFindMarkers(combined.s, assay = "SCT")
for (cell in cell_type) {
  celltype_subset <- subset(combined.s, idents = cell)
  
  de_results <- FindMarkers(celltype_subset,
                            ident.1 = "WT",
                            ident.2 = "Het",
                            group.by = "orig.ident",  
                            assay = "SCT",
                            logfc.threshold = 0.25,
                            min.pct = 0.1,
                            only.pos = FALSE,
                            recorrect_umi = FALSE)

  current_cluster <- paste0("DE_", gsub(" ", "_", cell), "_WT_vs_Het")
  assign(current_cluster, de_results)
  
  if (nrow(de_results) > 0) {
    write.csv(de_results, file = file.path(de_dir, paste0("DE_", gsub(" ", "_", cell), "_WT_vs_Het.csv")))
  }
}

### For showcase purposes, I have decided to focus on the DE of Fibroblasts.
fibro_only <- subset(combined.s, subset = cell_type == "Fibroblasts")

fibro_top_genes <- DE_Fibroblasts_WT_vs_Het[abs(DE_Fibroblasts_WT_vs_Het$avg_log2FC) > 1 
                                            & DE_Fibroblasts_WT_vs_Het$p_val_adj < 0.05,]

# Dotplot
d1 <- DotPlot(fibro_only, features = rownames(fibro_top_genes), group.by = "orig.ident") +
  labs(title = "DEGs Fibroblast only")
d2 <- DotPlot(combined.s, features = rownames(fibro_top_genes), group.by = "orig.ident") +
  labs(title = "DEGs Combined")

ggsave(filename = file.path(de_dir, "dotplot_DEGs_fibroblast_only.png"), plot = d1, width = 8, height = 6, dpi=600)
ggsave(filename = file.path(de_dir, "dotplot_DEGs_combined.png"), plot = d2, width = 8, height = 6, dpi=600)

### Pathway Enrichment
pe_dir <- paste0(de_dir, "/Pathway_enrichment")
dir.create(pe_dir)

fibro_sig_genes <- DE_Fibroblasts_WT_vs_Het[abs(DE_Fibroblasts_WT_vs_Het$avg_log2FC) > 0.25 
                                            & DE_Fibroblasts_WT_vs_Het$p_val_adj < 0.05,]


fibro_sig_genes$entrez_id <- mapIds(org.Mm.eg.db,
                     keys = rownames(fibro_sig_genes),
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")


## Running GO Enrichment on all 3 gene lists (Total, up, down) and ont (BP, MF, CC) based on gene_list (complete list of significant DEGs)
#Assumes that no list comes out to 0 genes
run_go <- function(gene_list) {
  
  onts <- c("BP", "CC", "MF")
  
  for (ont in onts) {
    
    ont_dir <- paste0(pe_dir, "/", ont)
  
    if (!dir.exists(ont_dir)) {
     dir.create(ont_dir, recursive = TRUE)
   }
  
    #Total
    ego <- enrichGO(gene  = gene_list$entrez_id,
           OrgDb        = org.Mm.eg.db,
           keyType      = "ENTREZID",
           ont          = ont,
           pAdjustMethod= "BH",
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2,
           readable     = TRUE)
  
    if (!is.null(ego) && nrow(ego@result) > 0) {
  
    d3 <- dotplot(ego, showCategory=10) + ggtitle(paste("Total GO Enrichment -", ont))
    ggsave(filename = file.path(ont_dir, paste0("Total_Fibroblast_GO_", ont, "_Enrichment_Dotplot.png")), plot = d3, width = 8, height = 6, dpi=600)
    write.csv(as.data.frame(ego), file = file.path(ont_dir, paste0("Total_Fibroblast_GO_", ont, "_Enrichment.csv")))
    }
  
    genes_up <- gene_list$entrez_id[gene_list$avg_log2FC > 0]
    genes_down <- gene_list$entrez_id[gene_list$avg_log2FC < 0]
  
    #Up 
    up <- enrichGO(gene  = genes_up,
                    OrgDb        = org.Mm.eg.db,
                    keyType      = "ENTREZID",
                    ont          = ont,
                    pAdjustMethod= "BH",
                   pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                  readable     = TRUE)
    
    if (!is.null(up) && nrow(up@result) > 0) {
    d3 <- dotplot(up, showCategory=10) + ggtitle(paste("Upregulated GO Enrichment -", ont))
    ggsave(filename = file.path(ont_dir, paste0("Upregulated_Fibroblast_GO_", ont, "_Enrichment_Dotplot.png")), plot = d3, width = 8, height = 6, dpi=600)
    write.csv(as.data.frame(up), file = file.path(ont_dir, paste0("Upregulated_Fibroblast_GO_", ont, "_Enrichment.csv")))
    }
  
    #Down
    down <- enrichGO(gene  = genes_down,
                 OrgDb        = org.Mm.eg.db,
                 keyType      = "ENTREZID",
                 ont          = ont,
                 pAdjustMethod= "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 readable     = TRUE)
  
    if (!is.null(down) && nrow(down@result) > 0) {
    d3 <- dotplot(down, showCategory=10) + ggtitle(paste("Downregulated GO Enrichment -", ont))
    ggsave(filename = file.path(ont_dir, paste0("Downregulated_Fibroblast_GO_", ont, "_Enrichment_Dotplot.png")), plot = d3, width = 8, height = 6, dpi=600)
    write.csv(as.data.frame(down), file = file.path(ont_dir, paste0("Downregulated_Fibroblast_GO_", ont, "_Enrichment.csv")))
    }
  }
}

##################################

run_go(fibro_sig_genes)


run_kegg <- function(gene_list) {
  
  ont_dir <- file.path(pe_dir, "KEGG")
  
  if (!dir.exists(ont_dir)) {
    dir.create(ont_dir, recursive = TRUE)
  }
  
  #Total
  kegg <- enrichKEGG(gene         = gene_list$entrez_id,
                           organism     = 'mmu',            
                           keyType      = "kegg",           
                           pAdjustMethod= "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.2)
  
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    d3 <- dotplot(kegg, showCategory = 10) + ggtitle("Total KEGG Enrichment")
    ggsave(filename = file.path(ont_dir, "Total_Fibroblast_KEGG_Enrichment_Dotplot.png"),
           plot = d3, width = 8, height = 6, dpi = 600)
    write.csv(as.data.frame(kegg),
              file = file.path(ont_dir, "Total_Fibroblast_KEGG_Enrichment.csv"))
  }
  
  genes_up <- gene_list$entrez_id[gene_list$avg_log2FC > 0]
  genes_down <- gene_list$entrez_id[gene_list$avg_log2FC < 0]
  
  #Up
  kegg_up <- enrichKEGG(gene         = genes_up,
                          organism     = 'mmu',
                          keyType      = "kegg",
                          pAdjustMethod= "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
  if (!is.null(kegg_up) && nrow(kegg_up@result) > 0) {
    d3 <- dotplot(kegg_up, showCategory = 10) + ggtitle("Upregulated KEGG Enrichment")
    ggsave(filename = file.path(ont_dir, "Upregulated_Fibroblast_KEGG_Enrichment_Dotplot.png"),
            plot = d3, width = 8, height = 6, dpi = 600)
    write.csv(as.data.frame(kegg_up),
                file = file.path(ont_dir, "Upregulated_Fibroblast_KEGG_Enrichment.csv"))
    }
  
  
  #Down
    kegg_down <- enrichKEGG(gene         = genes_down,
                            organism     = 'mmu',
                            keyType      = "kegg",
                            pAdjustMethod= "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.2)
    if (!is.null(kegg_down) && nrow(kegg_down@result) > 0) {
      d3 <- dotplot(kegg_down, showCategory = 10) + ggtitle("Downregulated KEGG Enrichment")
      ggsave(filename = file.path(ont_dir, "Downregulated_Fibroblast_KEGG_Enrichment_Dotplot.png"),
             plot = d3, width = 8, height = 6, dpi = 600)
      write.csv(as.data.frame(kegg_down),
                file = file.path(ont_dir, "Downregulated_Fibroblast_KEGG_Enrichment.csv"))
    }
  
}

##################################

run_kegg(fibro_sig_genes)
