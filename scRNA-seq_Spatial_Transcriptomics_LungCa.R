# scRNA-seq Spatial Transcriptomics Analysis
# 10X Genomics Dataset: Visium HD Spatial Gene Expression Libraries, Post-Xenium, Human Lung Cancer (FFPE)
# For this experiment, FFPE-preserved tissue was purchased from Avaden Biosciences (Lung Cancer, Invasive Acinar Adenocarcinoma, IB, T2a N0 MX). 
# M. GAGE
# 08/2025



#----SETUP THE SEURAT OBJECT----

# Load in the appropriate packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(arrow)

# Load in the Dataset. Here, I am loading the 8 um resolution data
object <- Load10X_Spatial(data.dir = "C:/Users/micha/OneDrive/Documents/Bioinformatics Masters/Applied Project Course/Lung_Cancer_Spatial", 
                          bin.size = c(8))

# Viewing the assay and the meta data (barcodes, ident, counts, and features)
Assays(object)
View(object@meta.data)



#----VISUALLY INSPECTING THE QUALITY OF DATA----

# Violin plots of counts and features
vln.plot <- VlnPlot(object, features = c("nFeature_Spatial.008um","nCount_Spatial.008um"), pt.size = 0, raster = FALSE)

# Counts feature plot. There are several regions of low counts due in part to low cellular density
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um", alpha = c(0.3, 1 ))

# Plots together
vln.plot | count.plot



#----NORMALIZE THE DATA----

# Standard log-normalization for spatial data
object <- NormalizeData(object)



#----VISUALIZE GENE EXPRESSION----

geneex <- SpatialFeaturePlot(object, features = "SCGB1A1", pt.size.factor = 2, alpha = c(0.4, 1)) + ggtitle("SCGB1A1 Expression (8um)")
geneex



#----UNSUPERVISED CLUSTERING----

# Note that data is already normalized
DefaultAssay(object) <- "Spatial.008um"
object <- FindVariableFeatures(object)
object <- ScaleData(object)

# We select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "Uniform",
  sketched.assay = "sketch"
)

# Switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# Perform clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:25)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 3)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:25)


# Create PCA on the full 8µm assay using features learned from the sketch

DefaultAssay(object) <- "Spatial.008um"
object <- RunPCA(
  object,
  assay = "Spatial.008um",
  features = VariableFeatures(object[["sketch"]]),
  reduction.name = "full.pca.sketch"
)

# Now we can project the cluster labels, and dimensional reductions (PCA and UMAP) that we learned 
# from the 50,000 sketched cells - to the entire dataset, using the ProjectData function.

object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:25,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# Visualize the clustering results for sketched cells and full dataset

DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# Switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2

# Visualizing the unsupervised clusters based on spatial location
SpatialDimPlot(object, label = T, repel = T, label.size = 4)

# Plotting the spatial location of one or more selected clusters. 
Idents(object) <- "seurat_cluster.projected"
cells <- CellsByIdentities(object, idents = c(7))
p <- SpatialDimPlot(object,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T,
                    alpha = c(0.3, 1 )
) + NoLegend()
p

# Finding and visualizing the top gene expression markers for each cluster

# Crete downsampled object to make visualization either
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

# Scaling the data subset and visualizing in a large heatmap
object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p



#----IDENTIFYING SPATIALLY-DEFINED TISSUE DOMAINS----

# BANKSY performs multiple tasks, but we find it particularly valuable for identifying and segmenting tissue domains. 
# When performing clustering, BANKSY augments a spot’s expression pattern with both the mean and the gradient of 
# gene expression levels in a spot’s broader neighborhood. In other words, Banksy generates clusters based on cell 
# context, not just their transcriptomes (like Seurat alone). 

library(SeuratWrappers)
library(Banksy)

# The RunBanksy function creates a new BANKSY assay, which can be used for dimensional reduction and clustering
object <- RunBanksy(object,
                    lambda = 0.8, verbose = TRUE,
                    assay = "Spatial.008um", slot = "data", features = "variable",
                    k_geom = 50
)

DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)

Idents(object) <- "banksy_cluster"
p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4)
p

# As with unsupervised clustering, we can highlight the spatial location of each tissue domain individually:
banksy_cells <- CellsByIdentities(object)

p <- SpatialDimPlot(
  object,
  cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")][0:30],
  cols.highlight = c("#FFFF00", "grey50"),
  facet.highlight = TRUE,
  combine = TRUE
) + NoLegend()
p




