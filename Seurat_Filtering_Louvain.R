library(dplyr)
library(Seurat)
library(patchwork)
library(ggrepel)


setwd("~/Desktop/Data for R/filtered_feature_bc_matrix")

########## Control #####

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "Control/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "control")
pbmc

pbmc <- RenameCells(object = pbmc,new.names = paste0("ct_", Cells(x = pbmc)))

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)

# VlnPlot(pbmc, features = c("nFeature_RNA"), ncol = 3, pt.size = 0.01, y.max = 200)
dim(pbmc) # 32285  4352
median(pbmc$nFeature_RNA) #374
min(pbmc$nFeature_RNA) #22

pbmc <- subset(pbmc, subset = nFeature_RNA > 50 & nFeature_RNA < 1500 & percent.mt < 5)

dim(pbmc) # 32285  4113
# Visualize QC metrics as a violin plot after filtering
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

pbmc <- SCTransform(pbmc, vars.to.regress = c("nCount_RNA","percent.mt"), verbose = FALSE, variable.features.n = 1000)

top8 <- head(VariableFeatures(pbmc), 8)
plot1 <- VariableFeaturePlot(pbmc, selection.method = "sct", assay = "SCT")
plot2 <- LabelPoints(plot = plot1, points = top8)
plot1
plot2

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca")
FeaturePlot(pbmc, features = "nFeature_RNA", reduction = "pca")
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.4) #0.4

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:15)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", label = TRUE, group.by = "SCT_snn_res.0.4")
FeaturePlot(pbmc, features = "Ins1")

new.cluster.ids.control <- c("Beta-cells", "pDc", "Delta-cells","Alpha", "T-cells", "MAC", "cDC", "PP cells", "Ductal cells", "double positive", "double positive")
names (new.cluster.ids.control) <- levels(pbmc)
control <- RenameIdents(pbmc, new.cluster.ids.control)
DimPlot(control, reduction = "umap", label = TRUE, pt.size = 0.5)+NoLegend()

features <- c("Ins1", "Ins2", "Gcg", "Ppy", "Sst","Cd3e", "Ptprf", "Ptpn2", "Ptpn22")
VlnPlot(control, features = features) + NoLegend()


pbmc$louvain <- paste("ct", pbmc$seurat_clusters, sep = "_")

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)%>% write.csv(.,file = "Data/AllClusters/Control_cluster_DE_genes.csv")

top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

saveRDS(pbmc, "Data/AllClusters/Control.rds")
control <- pbmc
write.csv(control$louvain, "Data/AllClusters/Control_leiden.csv")

########## Special diet #####
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "SpecialDiet/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "diet")
pbmc
rm(pbmc.data)

pbmc <- RenameCells(object = pbmc,new.names = paste0("diet_", Cells(x = pbmc)))

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

dim(pbmc) # 32285  5022
median(pbmc$nFeature_RNA) #361
min(pbmc$nFeature_RNA) #31

pbmc <- subset(pbmc, subset = nFeature_RNA > 50 & nFeature_RNA < 1500 & percent.mt < 5)

# Visualize QC metrics as a violin plot after filtering
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

dim(pbmc) # 32285  4850
pbmc <- SCTransform(pbmc, vars.to.regress = c("nCount_RNA","percent.mt"), verbose = FALSE, variable.features.n = 1000)

top8 <- head(VariableFeatures(pbmc), 8)
plot1 <- VariableFeaturePlot(pbmc, selection.method = "sct", assay = "SCT")
plot2 <- LabelPoints(plot = plot1, points = top8)
plot2

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.4)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap",  label = TRUE)
FeaturePlot(pbmc, features = "Ins1")

new.cluster.ids.diets <- c("pDC", "B-cells", "Beta-cells","T-cells", "Mac", "Delta-cells", "Alpha-cells", "Ductal cells", "PP cells", "Acinar")
names (new.cluster.ids.diets) <- levels(pbmc)
diet <- RenameIdents(pbmc, new.cluster.ids.diets)
DimPlot(diet, reduction = "umap", label = TRUE, pt.size = 0.5)+NoLegend()

features <- c("Ins1", "Ins2", "Gcg", "Ppy", "Sst","Cd3e", "Ptprf", "Ptpn2", "Ptpn22")
VlnPlot(diet, features = features) + NoLegend()

pbmc$louvain <- paste("diet", pbmc$seurat_clusters, sep = "_")


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)%>% write.csv(.,file = "Data/AllClusters/Diet_cluster_DE_genes_res1.3.csv")

saveRDS(pbmc, "Data/AllClusters/Diet.rds")
diet <- pbmc
write.csv(diet$louvain, "Data/AllClusters/Diet_leiden.csv")


####### merge ######

# control
pbmc.data <- Read10X(data.dir = "Control/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "control")
pbmc <- RenameCells(object = pbmc,new.names = paste0("ct_", Cells(x = pbmc)))

#transfer louvain info
metadata <- read.csv('Data/Control_leiden.csv', row.names = 1)
head(metadata)

cells <- rownames(metadata)
cell_types <- metadata$x
pbmc <- subset(pbmc, cells = cells)
pbmc <- AddMetaData(object = pbmc, metadata = cell_types, col.name = 'louvain')
control <- pbmc


## diet 
pbmc.data <- Read10X(data.dir = "SpecialDiet/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "diet")
pbmc <- RenameCells(object = pbmc,new.names = paste0("diet_", Cells(x = pbmc)))

#transfer louvain info
metadata <- read.csv('Data/Diet_leiden.csv', row.names = 1)
head(metadata)

cells <- rownames(metadata)
cell_types <- metadata$x
pbmc <- subset(pbmc, cells = cells)
pbmc <- AddMetaData(object = pbmc, metadata = cell_types, col.name = 'louvain')
diet <- pbmc
rm(pbmc)

# merge data
pbmc.combined <- sp::merge(control, y = diet, project = "Esteban")
pbmc.combined

pbmc.combined[["percent.mt"]] <- PercentageFeatureSet(pbmc.combined, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

pbmc.combined <- SCTransform(pbmc.combined, vars.to.regress = c("nCount_RNA","percent.mt"), verbose = FALSE, variable.features.n = 1000)

pbmc.combined <- RunPCA(pbmc.combined, features = VariableFeatures(object = pbmc.combined))
ElbowPlot(pbmc.combined)

dims.use = 1:15
pbmc.combined <- FindNeighbors(pbmc.combined, dims = dims.use)
pbmc.combined <- FindClusters(pbmc.combined, resolution = 0.4)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc.combined <- RunUMAP(pbmc.combined, dims = dims.use)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc.combined, reduction = "umap", group.by = "orig.ident")

DimPlot(pbmc.combined, reduction = "umap", group.by = "louvain")
FeaturePlot(pbmc.combined, features = "Ins1")
DimPlot(pbmc.combined, reduction = "umap")

new.cluster.ids.combined <- c("Beta-cells", "pDc", "Delta-cells","Alpha", "T-cells", "MAC", "cDC", "PP cells", "Ductal cells", "double positive", "double positive", "Extra", "extra")
names (new.cluster.ids.combined) <- levels(pbmc.combined)
combined <- RenameIdents(pbmc.combined, new.cluster.ids.combined)
DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5)+NoLegend()

features <- c("Ins1", "Ins2", "Gcg", "Ppy", "Sst","Cd3e", "Ptprf", "Ptpn2", "Ptpn22")
VlnPlot(combined, features = features) + NoLegend()

Idents(pbmc.combined) <- pbmc.combined$orig.ident

table(pbmc.combined$seurat_clusters)
table(pbmc.combined$seurat_clusters, pbmc.combined$orig.ident)

pbmc.markers <- FindAllMarkers(pbmc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)%>% write.csv(.,file = "Data/AllClusters/Combined_cluster_DE_genes.csv")

diet.markers <- FindMarkers(pbmc.combined, ident.1 = "control",ident.2 = "diet", logfc.threshold = 0.25, test.use = "roc")
head(diet.markers, n= 5)

write.csv(diet.markers, file = "diet_vs_control.csv")

ct.markers <- FindMarkers(pbmc.combined, ident.1 = "diet",ident.2 = "control", logfc.threshold = 0.25, test.use = "roc")
write.csv(ct.markers, file = "Data/AllClusters/control_vs_diet.csv")

saveRDS(pbmc.combined, "Data/AllClusters/Control_Diet_merged.rds")

features <- c("ins1", "ins2")
VlnPlot(pbmc, features = features, ncol = 2)
