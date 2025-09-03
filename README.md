# 🧬 Spatial Transcriptomics of Lung Cancer Tissue

## 📖 Project Overview  
This project applies **spatial transcriptomics** analysis to lung cancer tissue, integrating histology with spatially resolved gene expression. Using **10x Genomics Visium HD**, **Seurat**, and **BANKSY spatial clustering**, we explore the spatial organization of gene expression, immune infiltration, and tumor microenvironment heterogeneity.  

The workflow demonstrates how computational clustering, spatial feature plots, and marker-based annotation can be combined to uncover tissue architecture and molecular programs in lung cancer.  

---

## 🧪 Methods & Workflow  
1. **Data Preprocessing**  
   - Input: Visium HD lung cancer dataset  
   - Filtering, normalization, and QC with Seurat  

2. **Dimensionality Reduction & Clustering**  
   - UMAP sketch clustering (50k cells) → projected full dataset  
   - BANKSY spatially-aware clustering for tissue segmentation  

3. **Marker Gene Analysis**  
   - **CD3D** → T-cell infiltration  
   - **SCGB1A1** → airway epithelial/club cells  
   - Heatmap of top marker genes  

4. **Visualization**  
   - SpatialDimPlots of clusters  
   - Spatial feature maps (CD3D, SCGB1A1)  
   - Violin plots of QC metrics  

---

## 📂 Repository Contents  

### 🔬 Analysis Notebooks & Scripts  
- `Lung_Cancer_RMarkdown.Rmd` — Main workflow (RMarkdown)  
- `scRNA-seq_Spatial_Transcriptomics_LungCa.R` — Supporting R script  

### 📊 Figures & Results  
- `Plot_BANKSY_Spatial_Separate_Clusters.jpeg` – BANKSY clusters visualized separately  
- `Plot_CD3D_Expression.jpeg` – CD3D expression distribution  
- `Plot_ClusteringSketch&Full.jpeg` – UMAP sketch vs full clustering  
- `Plot_HeatmapTop5.jpeg` – Top 5 marker genes heatmap  
- `Plot_SpatialCluster7.jpeg` – Spatial overlay for cluster 7  
- `Plot_SpatialClusters.jpeg` – Tissue-wide spatial clustering  
- `Plot_SpatialDim_BANKSY_Clusters.jpeg` – BANKSY clustering overlay  
- `Plot_SpatialFeatureSCGB1A1.jpeg` – SCGB1A1 expression distribution  
- `Plot_Violin&SpatialFeature.jpeg` – Violin plots & feature intensity  

---

## 🔧 Tools & Dependencies  
- **Language**: R (≥4.3)  
- **Packages**:  
  - [Seurat](https://satijalab.org/seurat/)  
  - [Banksy](https://github.com/)  
  - [ggplot2](https://ggplot2.tidyverse.org/)  
  - [dplyr](https://dplyr.tidyverse.org/)  
