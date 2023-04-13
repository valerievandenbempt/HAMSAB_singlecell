set.seed(1732)
library(viridis)
library(dplyr)
library(Seurat)
library(scales)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(sctransform)
library(scDblFinder)
library(SingleCellExperiment)
library(cluster)
library(factoextra)
library(intrinsicDimension)
library(xlsx)
library(tibble)
colmaps <- c("#a6cee3",
             "#1f78b4", 
             "#b2df8a",
             "#33a02c",
             "#fb9a99",
             "#e31a1c",
             "#fdbf6f",
             "#ff7f00",
             "#cab2d6",
             "#6a3d9a",
             "#ffd92f",
             "#b15928",
             "#737373")
colmaps2 <- c("#ff7f0e","#2ca02c","#d62728","#9467bd")
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
pal <- c("#f7fcfd","#e0ecf4","#bfd3e6", "#9ebcda","#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b")

setwd("~/Elif/Collaborations/Esteban/scRNA/Run1/from_Valerie_final")
setwd("~/Desktop/Data for R/filtered_feature_bc_matrix")

# 1) Loading data (4352 cells with 32285 captured genes)

data.raw <- Read10X(data.dir = "./SpecialDiet/")

data <- CreateSeuratObject(counts = data.raw, names.field = 2, names.delim = "-")

rm(data.raw)

#Compute the % of MT and Ribosomal genes for each cell
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^Rp[sl][[:digit:]]")

# 2) QC plots before filtering

VlnPlot(data, features = c("nFeature_RNA","nCount_RNA", "percent.mt", "percent.ribo" ),
        ncol = 4, pt.size = 0.001)

FeatureScatter(data, feature1 = "percent.mt", feature2 = "percent.ribo")

data

# The percentage of MT expressed genes is not correlated with the percentage of Ribosomal expressed genes,
# so filtering by MT will not especially affect the other
# 
# 3) QC filtering
# 
# 3.1 low RNA counts 
VlnPlot(data, features = "nFeature_RNA", y.max = 500)
VlnPlot(data, features = "percent.mt", y.max = 30)

data1 <- subset(data, subset = nFeature_RNA > 100 & percent.mt < 10)
VlnPlot(data1, features = c("nFeature_RNA","nCount_RNA", "percent.mt", "percent.ribo" ),
        ncol = 4, pt.size = 0.001)

data1
# 4301 cells and 32285 genes

# 3.2 Filtering doublets/multiplets with scDblFinder package (42 potential doublets found)

#Isolate Seurat raw count matrix
counts <- as.matrix(data1@assays$RNA@data)
#Convert into sce object
sce <- SingleCellExperiment(assays = list(counts = counts))
#Run scDblFinder package on it
sce <- scDblFinder(sce)

#184 DOUBLETS (4.3% of cells) FOUND


#Transpose the matrix to work easily on it
counts_t <- t(counts)
#Then turn it into dataframe
counts_t <- as.data.frame(counts_t)

#Convert as matrix the scDblFinder results (singlet >< doublet)
matrix <- as.data.frame(colData(sce)$scDblFinder.class)

#Add a state variable for each cell (singlet >< doublet)
counts_t$state <- matrix$`colData(sce)$scDblFinder.class`
#Filter to only keep singlets
counts_t <- filter(counts_t, state == 'singlet')
#Remove the state column
counts_t$state <- NULL

#Create a new df of the count matrix for further analysis
counts_try <- t(counts_t)
counts_try <- as.data.frame(counts_try)

# 3.3 Re-Organize the data for some plots
#Create a df which sums the count for each gene and the number of cells who express it
genes <- tibble(
  gene = rownames(counts_try),
  count = Matrix::rowSums(counts_try),
  cells = Matrix::rowSums(counts_try != 0)
)

head(genes, n = 10)

###Try as best as possible to make a plot as in the scanpy tutorial, preproccessing part chunk 8

#Order genes by their count
genes_top <- genes[order(-genes$count),]
#Select the top 20
top20genes <- genes_top[1:20,]
#Compute the total count
total_count <- sum(genes$count) 

#Add a variable which computes the percentage of total count for the top 20 genes
top20genes$per_tot_count <- with(top20genes, count/total_count * 100)

#Necessary for the incoming plot
top20genes$gene <- factor(top20genes$gene, levels = top20genes$gene)

#Plot the result with wonderful colors
ggplot(data=top20genes, aes(x=per_tot_count, y=gene, fill = gene)) + geom_bar(stat="identity") + theme(legend.position = "none")     

# 3.4 Filter genes which are expressed in less than 2 cells (15910 genes kept over 32285)

#List of kept genes
genes_to_keep <- genes %>% dplyr::filter(cells >= 2) %>% pull(gene)
#Length of the list
length(genes_to_keep)

#Remove filtered genes in the final matrix
liste <- as.vector(genes_to_keep)
counts <- counts[rownames(counts_try) %in% liste ,]
# 3.5 Restructure the count matrix and convert it into Seurat object for final QC filtered violin plots

#Restructure the matrix
counts <- as.data.frame(counts_try)

#Compute nCount_rna for incoming violin plots
nCount_rna <- as.numeric(colSums(counts_try))

#Convert the matrix into Seurat

data_QC <- CreateSeuratObject(counts = counts)
data_QC <- RenameCells(object = data_QC,new.names = paste0("diet_", Cells(x = data_QC)))

data_QC
#Replace the nCount_rna inh the same "directory" than non- previously converted Seurat objects
data_QC[["nCount_RNA"]] <- nCount_rna

#Compute the % of MT and Ribosomal genes for each cell
data_QC[["percent.mt"]] <- PercentageFeatureSet(data_QC, pattern = "^mt-")
data_QC[["percent.ribo"]] <- PercentageFeatureSet(data_QC, pattern = "^Rp[sl][[:digit:]]")
VlnPlot(data_QC, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"),
        ncol = 4) + theme(legend.position = "none") 

data_QC


# 4) Normalization using SCTransform package + plots
# 
# It returns the Pearson residuals from regularized negative binomial regression, 
# where cellular sequencing depth is utilized as a covariate in a generalized linear model. 
# It aims to remove the influence of technical characteristics from downstream analyses 
# while preserving biological heterogeneity.

#SCT normalisation with top 1000 variable features
data_QC <- SCTransform(data_QC, variable.features.n = 1000, verbose = F)
#Variable feature plot
VariableFeaturePlot(data_QC, selection.method = "sct")

# #isolate the residual variance from the Seurat object
# residual_variance <- as.data.frame(data_QC@assays[["SCT"]]@meta.features[["sct.residual_variance"]])
# colnames(residual_variance) <- "residual_variance"
# 
# #Plot the ECDF of the residual variance
# ggplot(residual_variance, aes(residual_variance)) + stat_ecdf(geom = "point")

# 5) Dimensional reductions

# I tested the dimensionality of the dataset using intrinsicDimension package. 
# The result was 9 (I had manually set it to 11 because it was giving good clustering score instead of 
#                   15, 20, 25, 50,… for which there had important outliers in the clustering) and 
# by selecting the first 9 PCs, it gave a silhouette score (see later) of 0.63, which is really good. 
# We finally find 11 clusters at a resolution of 0.5 (actually the one who gives the best silouette score)

#Run PCA with variable features
data_QC <- RunPCA(data_QC, assay = "SCT",features = VariableFeatures(object = data_QC))

#Test the dimensionality
intrinsicDimension::maxLikGlobalDimEst(data_QC@reductions$pca@cell.embeddings,
                                       k = 10)

#Elbow plot
ElbowPlot(data_QC, ndims = 50)

dims.use = 1:11
#Umap
data_QC <- RunUMAP(data_QC, dims = dims.use, verbose = F)
#FindNeighbors
data_QC <- FindNeighbors(data_QC, dims = dims.use, verbose = F)
#Findclusters
data_QC <- FindClusters(data_QC, resolution = 0.5, verbose = F)

#Final Umap
DimPlot(data_QC, label = TRUE)


# 6) Clustering quality using silhouette plot
# 
# Some metrics such as silhouette coefficient can be used in order to evaluate the goodness of clustering 
# algorithm results. It processes as follows:
#   
#   -- Calculation of the average dissimilarity ai between each observation i and 
# all points belonging to i
# 
#   -- Calculation of the average dissimilarity d(i,C) for all clusters where i does not belong to.
# 
#   -- Define the silhouette width (Si)
# 
# The silhouette width varies between -1 and 1, for which a negative score means that observations may 
# probably be in the wrong cluster, while positive score means that observations are well clustered. 


#Compute distance matrix to UMAP coordinates
distance_matrix <- dist(Embeddings(data_QC[['umap']])[, 1:2])

#Isolate cluster nameq
clusters <- data_QC@active.ident

#Compute silhouette score for each cluster
silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
data_QC@meta.data$silhouette_score <- silhouette[,3]

#Plot
fviz_silhouette(silhouette)

saveRDS(data_QC, "Diet_doubletremoved.rds")

#### merge data ####

control <- readRDS("Control_doubletremoved.rds")
diet <- readRDS("Diet_doubletremoved.rds")
control$orig.ident <- rep(c("HAMS"))
diet$orig.ident <- rep(c("HAMSAB"))
combined2 <- merge(control, y = diet)
combined2

combined2 <- SCTransform(combined2, verbose = FALSE, variable.features.n = 1000)

combined2 <- RunPCA(combined2, features = VariableFeatures(object = combined2))


intrinsicDimension::maxLikGlobalDimEst(combined2@reductions$pca@cell.embeddings,
                                       k = 10)
#11 PC

#Elbow plot
ElbowPlot(combined2, ndims = 50)

dims.use = 1:16
#Umap
combined2 <- RunUMAP(combined2, dims = dims.use, verbose = F)
#FindNeighbors
combined2 <- FindNeighbors(combined2, dims = dims.use, verbose = F)
#Findclusters
combined2 <- FindClusters(combined2, resolution = 0.6, verbose = F)



#Final Umap
DimPlot(combined2, label = TRUE)



saveRDS(combined2, "Combined_doubletremoved.rds")

#### annotations ####

combined <- readRDS("Combined_doubletremoved_finetuned.rds")
combined2 <- readRDS("Combined_doubletremoved.rds")

#Compute distance matrix to UMAP coordinates
distance_matrix <- dist(Embeddings(combined[['umap']])[, 1:2])

#Isolate cluster nameq
clusters <- combined@active.ident

#Compute silhouette score for each cluster
silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
combined@meta.data$silhouette_score <- silhouette[,3]

#Plot
fviz_silhouette(silhouette)
#average silhouette width: 0.44
png(filename="UMAP.png", height = 2800, width = 4000,res=300, bg="transparent")
DimPlot(combined2, label = TRUE, repel = FALSE, pt.size = 1.5, label.size = 10) +NoLegend() &  theme(title = element_text(size = 30),
                                                      axis.text.x=element_text(hjust=0.5 , size=20),
                                                      axis.title = element_text(size=20),
                                                      axis.text.y = element_text(hjust = 0.5, size = 20),
                                                      axis.line = element_line(size=1.5), 
                                                      legend.text = element_text(size = 20), 
                                                      axis.ticks = element_line(size = 3),
                                                      axis.ticks.length = unit(7, "pt"),
                                                      legend.background = element_rect(fill = "transparent"),
                                                      panel.background = element_rect(fill = "transparent"),
                                                      plot.background = element_rect(fill = "transparent"))
dev.off()

png(filename="Trbc2.png", height = 2800, width = 4000, res = 300)
FeaturePlot(combined2, "Cd79a", pt.size = 1.5) + theme(title = element_text(size = 30),
                                                     axis.text.x=element_text(hjust=0.5 , size=20),
                                                     axis.title = element_text(size=20),
                                                     axis.text.y = element_text(hjust = 0.5, size = 20),
                                                     axis.line = element_line(size=1.5),
                                                    legend.text = element_text(size = 15), 
                                                    axis.ticks = element_line(size = 3),
                                                    axis.ticks.length = unit(7, "pt"))
dev.off()

combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

metadata <- read.csv("Combined_metadata.csv")

filter(metadata, celltypes == "Acinar cells") %>% select(X) %>% head(1)

plot <- DimPlot(combined, label = TRUE,group.by = "SCT_snn_res.0.6", repel = TRUE) 
LabelPoints(plot = plot, points = "ct_AAAGTCCTCTCAGGCG-1", repel = TRUE)

FeaturePlot(combined, c("Actb"))

FeaturePlot(combined, c("Gcg", "Ins1","Sst","Ppy")) + plot

#4 b cells
#0 beta
#8 beta
#6 alpha
#5 delta
#2 t
#1 cdc
#7 cdc
#3 mf
#9 acinar
#11 pdc
#10 pp
#14 endo
#13 b cells
#12 alpha/delta
#15 alpha/delta/beta


new.cluster.ids <- c("Beta cells","cDC","T cells","MF","B cells","Delta cells","Alpha cells","cDC","Beta cells",
                     "Acinar cells","PP cells","pDC","Alpha/Delta","B cells","Endothelial","Alpha/Delta/Beta")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)

Idents(combined) -> combined$celltypes

colmaps <- c("#a6cee3",
             "#1f78b4", 
             "#b2df8a",
             "#33a02c",
             "#fb9a99",
             "#e31a1c",
             "#fdbf6f",
             "#ff7f00",
             "#cab2d6",
             "#6a3d9a",
             "#ffd92f",
             "#b15928",
             "#737373")

DimPlot(combined2, group.by = "celltypes",split.by = "orig.ident", label=T, repel = T, cols = colmaps)
DimPlot(combined2, group.by = "celltypes",label = T, repel = T) 

combined$UMAP_1 <- combined@reductions[["umap"]]@cell.embeddings[,1]
combined$UMAP_2 <- combined@reductions[["umap"]]@cell.embeddings[,2]


saveRDS(combined, "Combined_doubletremoved.rds")

write.csv(combined@meta.data, "Combined_metadata_doubletremoved.csv")

plot <- FeaturePlot(combined, features= "Cd8a")
plot
#### fine tuning ####
combined <- readRDS('Combined_doubletremoved.rds')
DimPlot(combined, group.by = "celltypes", label=T, repel = T, cols = colmaps)

VlnPlot(combined, features = c('nCount_RNA', 'nFeature_RNA'), group.by = "celltypes")
VlnPlot(combined, features = c('nCount_RNA', 'nFeature_RNA'), group.by = "orig.ident")


combined2 <- subset(combined, subset = nFeature_RNA < 1500 & nCount_RNA <2500)
combined2
# 8716 cells left from 8867 cells x 47427 features
DimPlot(combined2, group.by = "celltypes", label=T, repel = T, cols = colmaps)

combined2 <- SCTransform(combined2, verbose = FALSE, variable.features.n = 1000)
combined2 <- RunPCA(combined2, features = VariableFeatures(object = combined2))
intrinsicDimension::maxLikGlobalDimEst(combined2@reductions$pca@cell.embeddings,
                                       k = 10)
# 11 PC

#Elbow plot
ElbowPlot(combined2, ndims = 50)

dims.use = 1:16
#Umap
combined2 <- RunUMAP(combined2, dims = dims.use, verbose = F)
#FindNeighbors
combined2 <- FindNeighbors(combined2, dims = dims.use, verbose = F)
#Findclusters
combined2 <- FindClusters(combined2, resolution = 0.6, verbose = F)

DimPlot(combined2, label = T)
DimPlot(combined2,group.by = 'celltypes', label = T)

new.cluster.ids <- c('cDC','Beta cells','B cells','Delta cells', 'MF', 'Alpha cells',
                     'T cells', 'T cells','cDC','Beta cells','Acinar cells','PP cells',
                     'pDC','Alpha/Delta','B cells','Endothelial cells')
names(new.cluster.ids) <- levels(combined2)
combined2 <- RenameIdents(combined2, new.cluster.ids)

Idents(combined2) -> combined2$celltypes
DimPlot(combined,group.by = 'celltypes', label = T)
colmaps <- c("#1f78b4",
             "#a6cee3",
             "#fb9a99",
             "#e31a1c",
             "#33a02c",
             "#fdbf6f",
             "#b2df8a",
             "#ff7f00",
             "#cab2d6",
             "#6a3d9a",
             "#ffd92f",
             "#b15928")

DimPlot(combined2, group.by = "celltypes",split.by = "orig.ident", label=T, repel = T, cols = colmaps)
DimPlot(combined2, group.by = "celltypes",label = T, repel = T, cols = colmaps) 

combined2$UMAP_1 <- combined2@reductions[["umap"]]@cell.embeddings[,1]
combined2$UMAP_2 <- combined2@reductions[["umap"]]@cell.embeddings[,2]


VlnPlot(combined2, features = c('nCount_RNA', 'nFeature_RNA'), group.by = "celltypes", cols = colmaps)
VlnPlot(combined2, features = c('nCount_RNA', 'nFeature_RNA'), group.by = "orig.ident")


saveRDS(combined2, "combined_doubletremoved_finetuned.rds")

write.csv(combined2@meta.data, "combined_metadata_doubletremoved_finetuned.csv")

combined <- readRDS('combined_doubletremoved_finetuned.rds')

#### T cells ####
tcell <- readRDS("tcells.rds")
tcells <- subset(tcell)
DimPlot(tcells)
# plot <- DimPlot(tcells)
# select.cells <- CellSelector(plot = plot)
# tcells <- subset(tcells, cells = select.cells)
tcells <- NormalizeData(tcells)
tcells <- FindVariableFeatures(tcells, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tcells)
tcells <- ScaleData(tcells, features = all.genes)

png("FoxP3_Violin.png",  width = 1500, height = 1500, res = 300)
cd4 <- subset(tcells, ident="1")
VlnPlot(cd4, features = "Foxp3", group.by = "orig.ident")
dev.off()

DoHeatmap(subset(cd4, downsample = 200), features = "Foxp3")

tcells <- RunPCA(tcells, features = VariableFeatures(object = tcells))
# intrinsicDimension::maxLikGlobalDimEst(tcells@reductions$pca@cell.embeddings,
#                                      k = 10)
# #21 PC
# 
# #Elbow plot
# ElbowPlot(tcells, ndims = 50)

dims.use = 1:21
#Umap
tcells <- RunUMAP(tcells, dims = dims.use, verbose = F)
#FindNeighbors
tcells <- FindNeighbors(tcells, dims = dims.use, verbose = F)
 #Findclusters
tcells <- FindClusters(tcells, resolution = 0.4, verbose = F)

png("foxp3.png", width = 1500, height = 1500, res = 300)
plot <- FeaturePlot(tcells, features= "Foxp3", group.by = "orig.ident")
plot
dev.off()

#Final Umap
# DimPlot(tcells, label = TRUE)

# 
# #Compute distance matrix to UMAP coordinates
# distance_matrix <- dist(Embeddings(tcells[['umap']])[, 1:2])
# 
# #Isolate cluster nameq
# clusters <- tcells@active.ident
# 
# #Compute silhouette score for each cluster
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# tcells@meta.data$silhouette_score <- silhouette[,3]
# 
# #Plot
# fviz_silhouette(silhouette)
# # 0.27
png("UMAP_tcells_Diet.png", width = 2700, height = 2700, res = 250, bg="transparent")
p2 <- DimPlot(tcells, group.by = "orig.ident", pt.size = 6.5, cols = c("#B19CD9", "#bfe6ff")) + labs(title = "Diet") &  theme(title = element_text(size=55), 
                                                                                                                             legend.position = c(0.65,1.07),
                                                                                                                             axis.text.x=element_text(hjust=1, size=55),
                                                                                                                             axis.title = element_text(size=55),
                                                                                                                             axis.ticks = element_line(size = 3),
                                                                                                                             axis.ticks.length = unit(7, "pt"),
                                                                                                                             axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                             legend.text=element_text(size=40),
                                                                                                                             legend.title=element_text(size=50),
                                                                                                                             legend.key = element_rect(size = 100),
                                                                                                                             axis.line = element_line(size=1.5), 
                                                                                                                             legend.background = element_rect(fill = "transparent"),
                                                                                                                             panel.background = element_rect(fill = "transparent"),
                                                                                                                             plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
p2
dev.off()
png("UMAP_tcells_SeuratClusters.png", width = 2700, height = 2700, res = 250, bg="transparent")
p3 <- DimPlot(tcells, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) + labs(title = "Subclusters") &  theme(title = element_text(size=55), 
                                                                                                                         legend.position = c(0.88, 1.0),
                                                                                                                         axis.text.x=element_text(hjust=1, size=55),
                                                                                                                         axis.title = element_text(size=55),
                                                                                                                         axis.text.y = element_text(hjust=0.5, size = 55),
                                                                                                                         axis.ticks = element_line(size = 3),
                                                                                                                         axis.ticks.length = unit(10, "pt"),
                                                                                                                         legend.text=element_text(size=40),
                                                                                                                         legend.title=element_text(size=55),
                                                                                                                         legend.key.height= unit(2, "cm"),
                                                                                                                         axis.line = element_line(size=1.5), 
                                                                                                                         legend.background = element_rect(fill = "transparent"),
                                                                                                                         panel.background = element_rect(fill = "transparent"),
                                                                                                                         plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )

p3
dev.off()



p1 <- DimPlot(tcells, group.by = "celltypes", pt.size = 3, cols = "#b2df8a") + labs(title = "T cells")
p2 <- DimPlot(tcells, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) + labs(title = "Diet") 
p3 <- DimPlot(tcells, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) + labs(title = "Subclusters")

plot_grid(p1,p2,p3, ncol = 3)


png("figures/combined/tcells.png", width = 2800, height = 1000, res = 200)

p1 <- DimPlot(tcells, group.by = "celltypes", pt.size = 3, cols = "#b2df8a") + labs(title = "T cells")
p2 <- DimPlot(tcells, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) + labs(title = "Diet") 
p3 <- DimPlot(tcells, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) + labs(title = "Subclusters")

plot_grid(p1,p2,p3, ncol = 3)

dev.off()


Idents(tcells) <- tcells$seurat_clusters
t.markers <- FindAllMarkers(tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(subset(tcells, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  scale_fill_gradientn(colours = rev(mapal))  + labs(title = 'DE genes between subclusters')

# 1. Open jpeg file
png("figures/combined/tcells_heatmap.png", width = 3000, height = 1500, res = 300)
# 2. Create the plot
DoHeatmap(subset(tcells, downsample = 200), features = top10$gene, group.colors = colmaps2) + 
  scale_fill_gradientn(colours = rev(mapal)) + labs(title = 'DE genes between subclusters')
# 3. Close the file
dev.off()

write.csv(t.markers, 'DEG/tcells_DEG_subclusters.csv')

t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

features <- c("Top2a","Hmgb2","Sept11","Zfp367",
              "Maf","Ikzf2","Malt1","Zfp36l1","Ptprc","H2-Eb1","Cd74","Gcg")


VlnPlot(tcells, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
VlnPlot(tcells, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident", ncol = 5)+ 
  theme(legend.position = 'right')
p1 <- DotPlot(tcells, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between subclusters')
p2 <- DotPlot(tcells, features = top10$gene, split.by = "orig.ident", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)
DoHeatmap(subset(tcells, downsample = 200), features = top10$gene, group.colors = colmaps2) + scale_fill_gradientn(colours = rev(mapal))


Idents(tcells) <- tcells$orig.ident
t.markers <- FindAllMarkers(tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(t.markers, 'DEG/tcells_DEG_diet.csv')
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

DoHeatmap(subset(tcells, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

VlnPlot(tcells, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
p1 <- DotPlot(tcells, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between diet')
# VlnPlot(tcells, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident")+ theme(legend.position = 'right')
p2 <- DotPlot(tcells, features = top10$gene, split.by = "seurat_clusters", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)

png("figures/combined/tcells_heatmap2.png", width = 2000, height = 1000, res = 200)

DoHeatmap(subset(tcells, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

dev.off()

saveRDS(tcells, "tcells.rds")
#### B cells ####
bcell <- readRDS("bcells.RDS")
bcells <- subset(bcell)
# DimPlot(bcells)
# plot <- DimPlot(bcells)
# select.cells <- CellSelector(plot = plot)
# bcells <- subset(bcells, cells = select.cells)

bcells <- NormalizeData(bcells)
bcells <- FindVariableFeatures(bcells, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(bcells)
bcells <- ScaleData(bcells, features = all.genes)

bcells <- RunPCA(bcells, features = VariableFeatures(object = bcells))
intrinsicDimension::maxLikGlobalDimEst(bcells@reductions$pca@cell.embeddings,
                                       k = 10)
#19 PC

#Elbow plot
ElbowPlot(bcells, ndims = 50)

dims.use = 1:19
#Umap
bcells <- RunUMAP(bcells, dims = dims.use, verbose = F)
#FindNeighbors
bcells <- FindNeighbors(bcells, dims = dims.use, verbose = F)
#Findclusters
bcells <- FindClusters(bcells, resolution = 0.4, verbose = F)

#Final Umap
DimPlot(bcells, label = TRUE)

png("UMAP_bcells_Diet.png", width = 2700, height = 2700, res = 250, bg="transparent")
p2 <- DimPlot(bcells, group.by = "orig.ident", pt.size = 6.5, cols = c("#B19CD9", "#bfe6ff")) + labs(title = "Diet") &  theme(title = element_text(size=55), 
                                                                                                                              legend.position = c(0.65,1.07),
                                                                                                                              axis.text.x=element_text(hjust=1, size=55),
                                                                                                                              axis.title = element_text(size=55),
                                                                                                                              axis.ticks = element_line(size = 3),
                                                                                                                              axis.ticks.length = unit(7, "pt"),
                                                                                                                              axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                              legend.text=element_text(size=40),
                                                                                                                              legend.title=element_text(size=50),
                                                                                                                              legend.key = element_rect(size = 100),
                                                                                                                              axis.line = element_line(size=1.5), 
                                                                                                                              legend.background = element_rect(fill = "transparent"),
                                                                                                                              panel.background = element_rect(fill = "transparent"),
                                                                                                                              plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
p2
dev.off()
png("UMAP_bcells_SeuratClusters.png", width = 2700, height = 2700, res = 250, bg= "transparent")
p3 <- DimPlot(bcells, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) + labs(title = "Subclusters") & theme(title = element_text(size=55), 
                                                                                                                          legend.position = c(0.88, 1.0),
                                                                                                                          axis.text.x=element_text(hjust=1, size=55),
                                                                                                                          axis.title = element_text(size=55),
                                                                                                                          axis.text.y = element_text(hjust=0.5, size = 55),
                                                                                                                          axis.ticks = element_line(size = 3),
                                                                                                                          axis.ticks.length = unit(10, "pt"),
                                                                                                                          legend.text=element_text(size=40),
                                                                                                                          legend.title=element_text(size=55),
                                                                                                                          legend.key.height= unit(2, "cm"),
                                                                                                                          axis.line = element_line(size=1.5),
                                                                                                                          legend.background = element_rect(fill = "transparent"),
                                                                                                                          panel.background = element_rect(fill = "transparent"),
                                                                                                                          plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )

p3
dev.off()


# #Compute distance matrix to UMAP coordinates
# distance_matrix <- dist(Embeddings(bcells[['umap']])[, 1:2])
# 
# #Isolate cluster nameq
# clusters <- bcells@active.ident
# 
# #Compute silhouette score for each cluster
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# bcells@meta.data$silhouette_score <- silhouette[,3]
# 
# #Plot
# fviz_silhouette(silhouette)
# # 0.57

p1 <- DimPlot(bcells, group.by = "celltypes", pt.size = 3, cols = "#fb9a99") + labs(title = "B cells")
p2 <- DimPlot(bcells, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) + labs(title = "Diet") 
p3 <- DimPlot(bcells, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) + labs(title = "Subclusters")

plot_grid(p1,p2,p3, ncol = 3)

png("figures/combined/bcells.png", width = 2800, height = 1000, res = 200)

p1 <- DimPlot(bcells, group.by = "celltypes", pt.size = 3, cols = "#fb9a99") + labs(title = "B cells")
p2 <- DimPlot(bcells, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) + labs(title = "Diet") 
p3 <- DimPlot(bcells, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) + labs(title = "Subclusters")

plot_grid(p1,p2,p3, ncol = 3)

dev.off()

Idents(bcells) <- bcells$seurat_clusters
t.markers <- FindAllMarkers(bcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

write.csv(t.markers, 'DEG/bcells_DEG_subclusters.csv')

DoHeatmap(subset(bcells, downsample = 200), features = top10$gene, group.colors = colmaps2) + 
  scale_fill_gradientn(colours = rev(mapal)) + labs(title = 'DE genes between subclusters')
VlnPlot(bcells, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
VlnPlot(bcells, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident", ncol = 5)+ 
  theme(legend.position = 'right')
p1 <- DotPlot(bcells, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between subclusters')
p2 <- DotPlot(bcells, features = top10$gene, split.by = "orig.ident", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)



Idents(bcells) <- bcells$orig.ident
t.markers <- FindAllMarkers(bcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(t.markers, 'DEG/bcells_DEG_diet.csv')
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

DoHeatmap(subset(bcells, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

VlnPlot(bcells, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
p1 <- DotPlot(bcells, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between diet')
# VlnPlot(bcells, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident")+ theme(legend.position = 'right')
p2 <- DotPlot(bcells, features = top10$gene, split.by = "seurat_clusters", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)

png("figures/combined/bcells_heatmap2.png", width = 2000, height = 1000, res = 200)

DoHeatmap(subset(bcells, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

dev.off()

saveRDS(bcells, 'bcells.rds')
#### cDC  ####

cdc <- subset(combined, ident = "cDC")
# DimPlot(cdc)
# plot <- DimPlot(cdc)
# select.cells <- CellSelector(plot = plot)
# cdc <- subset(cdc, cells = select.cells)

cdc <- NormalizeData(cdc)
cdc <- FindVariableFeatures(cdc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cdc)
cdc <- ScaleData(cdc, features = all.genes)

cdc <- RunPCA(cdc, features = VariableFeatures(object = cdc))
# intrinsicDimension::maxLikGlobalDimEst(cdc@reductions$pca@cell.embeddings,
#                                        k = 10)
# #20 PC

# #Elbow plot
# ElbowPlot(cdc, ndims = 50)

dims.use = 1:16
#Umap
cdc <- RunUMAP(cdc, dims = dims.use, verbose = F)
#FindNeighbors
cdc <- FindNeighbors(cdc, dims = dims.use, verbose = F)
#Findclusters
cdc <- FindClusters(cdc, resolution = 0.1, verbose = F)

# #Final Umap
# DimPlot(cdc, label = TRUE)

png("UMAP_cdc_Diet.png", width = 2700, height = 2700, res = 250)
p2 <- DimPlot(cdc, group.by = "orig.ident", pt.size = 6.5, cols = c("#B19CD9", "#bfe6ff")) + labs(title = "Diet") +  theme(title = element_text(size=55), 
                                                                                                                              legend.position = c(0.65,1.07),
                                                                                                                              axis.text.x=element_text(hjust=1, size=55),
                                                                                                                              axis.title = element_text(size=55),
                                                                                                                              axis.ticks = element_line(size = 3),
                                                                                                                              axis.ticks.length = unit(7, "pt"),
                                                                                                                              axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                              legend.text=element_text(size=40),
                                                                                                                              legend.title=element_text(size=50),
                                                                                                                              legend.key = element_rect(size = 100),
                                                                                                                              axis.line = element_line(size=1.5)) +
  guides(color = guide_legend(override.aes = list(size = 10) ) )
p2
dev.off()
png("UMAP_cdc_SeuratClusters.png", width = 2700, height = 2700, res = 250)
p3 <- DimPlot(cdc, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) + labs(title = "Subclusters")+  theme(title = element_text(size=55), 
                                                                                                                          legend.position = c(0.88, 1.0),
                                                                                                                          axis.text.x=element_text(hjust=1, size=55),
                                                                                                                          axis.title = element_text(size=55),
                                                                                                                          axis.text.y = element_text(hjust=0.5, size = 55),
                                                                                                                          axis.ticks = element_line(size = 3),
                                                                                                                          axis.ticks.length = unit(10, "pt"),
                                                                                                                          legend.text=element_text(size=40),
                                                                                                                          legend.title=element_text(size=55),
                                                                                                                          legend.key.height= unit(2, "cm"),
                                                                                                                          axis.line = element_line(size=1.5)) +
  guides(color = guide_legend(override.aes = list(size = 10) ) )

p3
dev.off()


# #Compute distance matrix to UMAP coordinates
# distance_matrix <- dist(Embeddings(cdc[['umap']])[, 1:2])
# 
# #Isolate cluster nameq
# clusters <- cdc@active.ident
# 
# #Compute silhouette score for each cluster
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# cdc@meta.data$silhouette_score <- silhouette[,3]
# 
# #Plot
# fviz_silhouette(silhouette)
# # 0.31

DimPlot(cdc, group.by = "celltypes", pt.size = 3, cols = c("#1f78b4")) +labs(title = "cDC") + NoLegend() +
  DimPlot(cdc, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) +
  DimPlot(cdc, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) 

png("figures/combined/cdc.png", width = 2800, height = 1000, res = 200)

p1 <- DimPlot(cdc, group.by = "celltypes", pt.size = 3, cols = "#1f78b4") + labs(title = "cDC")
p2 <- DimPlot(cdc, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) + labs(title = "Diet") 
p3 <- DimPlot(cdc, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) + labs(title = "Subclusters")

plot_grid(p1,p2,p3, ncol = 3)

dev.off()

Idents(cdc) <- cdc$seurat_clusters
t.markers <- FindAllMarkers(cdc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

write.csv(t.markers, 'DEG/cdc_DEG_subclusters.csv')

DoHeatmap(subset(cdc, downsample = 200), features = top10$gene, group.colors = colmaps2) + 
  scale_fill_gradientn(colours = rev(mapal)) +labs(title = 'DE genes between subclusters')

VlnPlot(cdc, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
VlnPlot(cdc, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident", ncol = 5)+ 
  theme(legend.position = 'right')
p1 <- DotPlot(cdc, features = top10$gene[-27], cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between subclusters')
p2 <- DotPlot(cdc, features = top10$gene[-27], split.by = "orig.ident", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)


Idents(cdc) <- cdc$orig.ident
t.markers <- FindAllMarkers(cdc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(t.markers, 'DEG/cdc_DEG_diet.csv')
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

DoHeatmap(subset(cdc, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

VlnPlot(cdc, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
p1 <- DotPlot(cdc, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between diet')
# VlnPlot(cdc, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident")+ theme(legend.position = 'right')
p2 <- DotPlot(cdc, features = top10$gene, split.by = "seurat_clusters", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)

png("figures/combined/cdc_heatmap2.png", width = 2000, height = 1000, res = 200)

DoHeatmap(subset(cdc, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

dev.off()

write.csv(t.markers, 'DEG/cdc_DEG_subclusters.csv')

saveRDS(cdc, "cdc.rds")

#### MF ####

mf <- subset(combined, ident = "MF")
# DimPlot(mf)
# plot <- DimPlot(mf)
# select.cells <- CellSelector(plot = plot)
# mf <- subset(mf, cells = select.cells)

mf <- NormalizeData(mf)
mf <- FindVariableFeatures(mf, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mf)
mf <- ScaleData(mf, features = all.genes)

mf <- RunPCA(mf, features = VariableFeatures(object = mf))
# intrinsicDimension::maxLikGlobalDimEst(mf@reductions$pca@cell.embeddings,
#                                        k = 10)
#21 PC

# #Elbow plot
# ElbowPlot(mf, ndims = 50)

dims.use = 1:21
#Umap
mf <- RunUMAP(mf, dims = dims.use, verbose = F)
#FindNeighbors
mf <- FindNeighbors(mf, dims = dims.use, verbose = F)
#Findclusters
mf <- FindClusters(mf, resolution = 0.4, verbose = F)
# #Final Umap
# DimPlot(mf, label = TRUE)
# 

png("UMAP_mf_Diet.png", width = 2700, height = 2700, res = 250)
p2 <- DimPlot(mf, group.by = "orig.ident", pt.size = 6.5, cols = c("#B19CD9", "#bfe6ff")) + labs(title = "Diet") +  theme(title = element_text(size=55), 
                                                                                                                              legend.position = c(0.65,1.07),
                                                                                                                              axis.text.x=element_text(hjust=1, size=55),
                                                                                                                              axis.title = element_text(size=55),
                                                                                                                              axis.ticks = element_line(size = 3),
                                                                                                                              axis.ticks.length = unit(7, "pt"),
                                                                                                                              axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                              legend.text=element_text(size=40),
                                                                                                                              legend.title=element_text(size=50),
                                                                                                                              legend.key = element_rect(size = 100),
                                                                                                                              axis.line = element_line(size=1.5)) +
  guides(color = guide_legend(override.aes = list(size = 10) ) )
p2
dev.off()
png("UMAP_mf_SeuratClusters.png", width = 2700, height = 2700, res = 250)
p3 <- DimPlot(mf, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) + labs(title = "Subclusters")+  theme(title = element_text(size=55), 
                                                                                                                          legend.position = c(0.88, 1.0),
                                                                                                                          axis.text.x=element_text(hjust=1, size=55),
                                                                                                                          axis.title = element_text(size=55),
                                                                                                                          axis.text.y = element_text(hjust=0.5, size = 55),
                                                                                                                          axis.ticks = element_line(size = 3),
                                                                                                                          axis.ticks.length = unit(10, "pt"),
                                                                                                                          legend.text=element_text(size=40),
                                                                                                                          legend.title=element_text(size=55),
                                                                                                                          legend.key.height= unit(2, "cm"),
                                                                                                                          axis.line = element_line(size=1.5)) +
  guides(color = guide_legend(override.aes = list(size = 10) ) )

p3
dev.off() 

# #Compute distance matrix to UMAP coordinates
# distance_matrix <- dist(Embeddings(mf[['umap']])[, 1:2])
# 
# #Isolate cluster nameq
# clusters <- mf@active.ident
# 
# #Compute silhouette score for each cluster
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# mf@meta.data$silhouette_score <- silhouette[,3]
# 
# #Plot
# fviz_silhouette(silhouette)
# # 0.21

DimPlot(mf, group.by = "celltypes", pt.size = 3, cols = c("#33a02c")) +labs(title = "MF") + NoLegend() +
  DimPlot(mf, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) +
  DimPlot(mf, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) 

png("figures/combined/mf.png", width = 2800, height = 1000, res = 200)

p1 <- DimPlot(mf, group.by = "celltypes", pt.size = 3, cols = "#33a02c") + labs(title = "MF")
p2 <- DimPlot(mf, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) + labs(title = "Diet") 
p3 <- DimPlot(mf, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) + labs(title = "Subclusters")

plot_grid(p1,p2,p3, ncol = 3)

dev.off()

t.markers <- FindAllMarkers(mf, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(t.markers, 'DEG/mf_DEG_subclusters.csv')

DoHeatmap(subset(mf, downsample = 200), features = top10$gene, group.colors = colmaps2) + 
  scale_fill_gradientn(colours = rev(mapal)) +labs(title = 'DE genes between subclusters')

VlnPlot(mf, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
VlnPlot(mf, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident", ncol = 5)+ 
  theme(legend.position = 'right')
p1 <- DotPlot(mf, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between subclusters')
p2 <- DotPlot(mf, features = top10$gene, split.by = "orig.ident", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)


Idents(mf) <- mf$orig.ident
t.markers <- FindAllMarkers(mf, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(t.markers, 'DEG/mf_DEG_diet.csv')
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

DoHeatmap(subset(mf, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

VlnPlot(mf, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
p1 <- DotPlot(mf, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between diet')
# VlnPlot(mf, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident")+ theme(legend.position = 'right')
p2 <- DotPlot(mf, features = top10$gene, split.by = "seurat_clusters", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)

png("figures/combined/mf_heatmap2.png", width = 2000, height = 1000, res = 200)

DoHeatmap(subset(mf, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

dev.off()


saveRDS(mf, "mf.rds")

#### pdc ####

pdc <- subset(combined, ident = "pDC")
# DimPlot(pdc)
# plot <- DimPlot(pdc)
# select.cells <- CellSelector(plot = plot)
# pdc <- subset(pdc, cells = select.cells)

pdc <- NormalizeData(pdc)
pdc <- FindVariableFeatures(pdc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pdc)
pdc <- ScaleData(pdc, features = all.genes)
pdc <- RunPCA(pdc, features = VariableFeatures(object = pdc))
intrinsicDimension::maxLikGlobalDimEst(pdc@reductions$pca@cell.embeddings,
                                       k = 10)
#16 PC

# #Elbow plot
# ElbowPlot(pdc, ndims = 50)

dims.use = 1:16
#Umap
pdc <- RunUMAP(pdc, dims = dims.use, verbose = F)
#FindNeighbors
pdc <- FindNeighbors(pdc, dims = dims.use, verbose = F)
#Findclusters
pdc <- FindClusters(pdc, resolution = 0.8, verbose = F)

# # #Final Umap
# DimPlot(pdc, label = TRUE)

png("UMAP_pdc_Diet.png", width = 2700, height = 2700, res = 250)
p2 <- DimPlot(pdc, group.by = "orig.ident", pt.size = 6.5, cols = c("#B19CD9", "#bfe6ff")) + labs(title = "Diet") +  theme(title = element_text(size=55), 
                                                                                                                              legend.position = c(0.65,1.07),
                                                                                                                              axis.text.x=element_text(hjust=1, size=55),
                                                                                                                              axis.title = element_text(size=55),
                                                                                                                              axis.ticks = element_line(size = 3),
                                                                                                                              axis.ticks.length = unit(7, "pt"),
                                                                                                                              axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                              legend.text=element_text(size=40),
                                                                                                                              legend.title=element_text(size=50),
                                                                                                                              legend.key = element_rect(size = 100),
                                                                                                                              axis.line = element_line(size=1.5)) +
  guides(color = guide_legend(override.aes = list(size = 10) ) )
p2
dev.off()
png("UMAP_pdc_SeuratClusters.png", width = 2700, height = 2700, res = 250)
p3 <- DimPlot(pdc, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) + labs(title = "Subclusters")+  theme(title = element_text(size=55), 
                                                                                                                          legend.position = c(0.88, 1.0),
                                                                                                                          axis.text.x=element_text(hjust=1, size=55),
                                                                                                                          axis.title = element_text(size=55),
                                                                                                                          axis.text.y = element_text(hjust=0.5, size = 55),
                                                                                                                          axis.ticks = element_line(size = 3),
                                                                                                                          axis.ticks.length = unit(10, "pt"),
                                                                                                                          legend.text=element_text(size=40),
                                                                                                                          legend.title=element_text(size=55),
                                                                                                                          legend.key.height= unit(2, "cm"),
                                                                                                                          axis.line = element_line(size=1.5)) +
  guides(color = guide_legend(override.aes = list(size = 10) ) )

p3
dev.off()

# # 
# #Compute distance matrix to UMAP coordinates
# distance_matrix <- dist(Embeddings(pdc[['umap']])[, 1:2])
# 
# #Isolate cluster nameq
# clusters <- pdc@active.ident
# 
# #Compute silhouette score for each cluster
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# pdc@meta.data$silhouette_score <- silhouette[,3]
# 
# #Plot
# fviz_silhouette(silhouette)
# # 0.63

DimPlot(pdc, group.by = "celltypes", pt.size = 8, cols = c("#6a3d9a")) +labs(title = "pDC") + NoLegend() +
  DimPlot(pdc, group.by = "orig.ident", pt.size = 8, cols = c("#1f77b4","#d62728")) +
  DimPlot(pdc, group.by = "seurat_clusters", pt.size = 8, cols = colmaps2) 

png("figures/combined/pdc.png", width = 2800, height = 1000, res = 200)

p1 <- DimPlot(pdc, group.by = "celltypes", pt.size = 3, cols = "#6a3d9a") + labs(title = "pdc")
p2 <- DimPlot(pdc, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) + labs(title = "Diet") 
p3 <- DimPlot(pdc, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) + labs(title = "Subclusters")

plot_grid(p1,p2,p3, ncol = 3)

dev.off()

t.markers <- FindAllMarkers(pdc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(t.markers, 'DEG/pdc_DEG_subclusters.csv')


DoHeatmap(subset(pdc, downsample = 200), features = top10$gene, group.colors = colmaps2) + 
  scale_fill_gradientn(colours = rev(mapal)) +labs(title = 'DE genes between subclusters')

VlnPlot(pdc, features = top10$gene, cols = colmaps2, ncol = 4) + theme(legend.position = 'right')
VlnPlot(pdc, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident", ncol = 4)+ 
  theme(legend.position = 'right')
p1 <- DotPlot(pdc, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between subclusters')
p2 <- DotPlot(pdc, features = top10$gene, split.by = "orig.ident", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)


Idents(pdc) <- pdc$orig.ident
t.markers <- FindAllMarkers(pdc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(t.markers, 'DEG/pdc_DEG_diet.csv')
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

DoHeatmap(subset(pdc, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

VlnPlot(pdc, features = top10$gene, cols = colmaps2, ncol = 7) + theme(legend.position = 'right')
p1 <- DotPlot(pdc, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between diet')
# VlnPlot(pdc, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident")+ theme(legend.position = 'right')
p2 <- DotPlot(pdc, features = top10$gene, split.by = "seurat_clusters", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)

png("figures/combined/pdc_heatmap2.png", width = 2000, height = 1000, res = 200)

DoHeatmap(subset(pdc, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

dev.off()


saveRDS(pdc, "pdc.rds")

#### Alpha cells ####
Alph <- readRDS("Alpha.rds")
Alpha <- subset(Alph)

Alpha <- NormalizeData(Alpha)
Alpha <- FindVariableFeatures(Alpha, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Alpha)
Alpha <- ScaleData(Alpha, features = all.genes)
Alpha <- RunPCA(Alpha, features = VariableFeatures(object = Alpha))
# intrinsicDimension::maxLikGlobalDimEst(Alpha@reductions$pca@cell.embeddings,
#                                        k = 10)
#10 PC

# #Elbow plot
# ElbowPlot(Alpha, ndims = 50)

dims.use = 1:19
#Umap
Alpha <- RunUMAP(Alpha, dims = dims.use, verbose = F)
#FindNeighbors
Alpha <- FindNeighbors(Alpha, dims = dims.use, verbose = F)
#Findclusters
Alpha <- FindClusters(Alpha, resolution = 0.4, verbose = F)

# # #Final Umap
# DimPlot(Alpha, label = TRUE)
# # 
# # 
# #Compute distance matrix to UMAP coordinates
# distance_matrix <- dist(Embeddings(Alpha[['umap']])[, 1:2])
# 
# #Isolate cluster nameq
# clusters <- Alpha@active.ident
# 
# #Compute silhouette score for each cluster
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# Alpha@meta.data$silhouette_score <- silhouette[,3]
# 
# #Plot
# fviz_silhouette(silhouette)
# # 0.63

DimPlot(Alpha, group.by = "celltypes", pt.size = 6.5, cols = c("#fdbf6f")) +labs(title = "Alpha cells") + NoLegend() +
  DimPlot(Alpha, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) +
  DimPlot(Alpha, group.by = "orig.ident", pt.size = 6.5, cols = c("#1f77b4","#d62728")) 
  

png("Alpha.png", width = 2800, height = 1000, res = 500)

#UMAPS for JDRF - Alpha
p1 <- DimPlot(Alpha, group.by = "celltypes", pt.size = 10, cols = "#fdbf6f") + labs(title = "Alpha cells")  +  theme(title = element_text(size = 30),
                                                                                                                                  axis.text.x=element_text(hjust=0.5 , size=20),
                                                                                                                                  axis.title = element_text(size=20),
                                                                                                                                  axis.text.y = element_text(hjust = 0.5, size = 20),
                                                                                                                                  axis.line = element_line(size=1.5))
png("UMAP_Alpha_Diet.png", width = 2700, height = 2700, res = 250, bg="transparent")
p2 <- DimPlot(Alpha, group.by = "orig.ident", pt.size = 6.5, cols = c("#B19CD9", "#bfe6ff")) + labs(title = "Diet") &  theme(title = element_text(size=55), 
                                                                                                                             legend.position = c(0.52, 0.9),
                                                                                                                            axis.text.x=element_text(hjust=1, size=55),
                                                                                                                            axis.title = element_text(size=55),
                                                                                                                            axis.ticks = element_line(size = 3),
                                                                                                                            axis.ticks.length = unit(7, "pt"),
                                                                                                                            axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                            legend.text=element_text(size=55),
                                                                                                                            legend.title=element_text(size=50),
                                                                                                                          legend.key = element_rect(size = 100),
                                                                                                                            axis.line = element_line(size=1.5))
                                                                                                                          & guides(color = guide_legend(override.aes = list(size = 10) ) )
p2
dev.off()
png("UMAP_Alpha_SeuratClusters.png", width = 2700, height = 2700, res = 250, bg="transparent")
p3 <- DimPlot(Alpha, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) + labs(title = "Subclusters") &  theme(title = element_text(size=55), 
                                                                                                                         legend.position = c(0.8, 0.9),
                                                                                                                         axis.text.x=element_text(hjust=1, size=55),
                                                                                                                         axis.title = element_text(size=55),
                                                                                                                         axis.text.y = element_text(hjust=0.5, size = 55),
                                                                                                                         axis.ticks = element_line(size = 3),
                                                                                                                         axis.ticks.length = unit(10, "pt"),
                                                                                                                         legend.text=element_text(size=55),
                                                                                                                         legend.title=element_text(size=55),
                                                                                                                         legend.key.height= unit(2, "cm"),
                                                                                                                         axis.line = element_line(size=1.5),
                                                                                                                         legend.background = element_rect(fill = "transparent"),
                                                                                                                         panel.background = element_rect(fill = "transparent"),
                                                                                                                         plot.background = element_rect(fill = "transparent")) &
                                                                                                                          guides(color = guide_legend(override.aes = list(size = 10) ) )

p3
dev.off()

plot_grid(p3,p2)

Idents(Alpha) <- Alpha$seurat_clusters
t.markers <- FindAllMarkers(Alpha, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

t.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
write.csv(t.markers, 'DEG/Alpha_DEG_subclusters.csv')

DoHeatmap(subset(Alpha, downsample = 200), features = top10$gene, group.colors = colmaps2) + 
  scale_fill_gradientn(colours = rev(mapal)) +labs(title = 'DE genes between subclusters')

#Violinpots for JDRF - Alpha
Features_AlphaCells_HAMSAB <- c("Gcg","Gpx3", "B2m","Cpe", "Pdcd4","Btg2","Ins2", "Fos" )
Dotplot <- c("Gcg", "Gnas", "Gm20721", "Malat1","Ins1", "Ins2", "Iapp", "Mafa" )

VlnPlot(Alpha, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
VlnPlot(Alpha, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident", ncol = 5)+ 
  theme(legend.position = 'right')

Idents(Alpha) <- Alpha$orig.ident
png("VlnPlot_Alpha_New.png", height = 1800, width = 6500, res = 250)
VlnPlot(Alpha, features = Features_AlphaCells_HAMSAB, cols = c("#B19CD9", "#bfe6ff"), ncol=8, y.max = 8) + theme(legend.position = "right") & 
  theme(title = element_text(vjust=2), text = element_text(size = 30), axis.ticks.y = element_line(size = 3), axis.ticks.length.y = unit(7, "pt"),
        axis.line = element_line(size=1.5), axis.text = element_text(size = 30), axis.title = element_text(size = 35))
dev.off()
png("VlnPlot_Alpha_1.png", height = 1700, width = 3500, res = 250)
VlnPlot(Alpha, features = Features_AlphaCells_HAMSAB_1, cols = c("hotpink4", "coral1"), ncol = 5, split.by = "orig.ident", idents = 1) + theme(legend.position = 'right')& 
  theme(text = element_text(size = 30), axis.line = element_line(size=1.5), axis.text = element_text(size = 30), axis.title= element_text(size=30), legend.text = element_text(size = 30))
dev.off()

png("DEG_dotplot_Alpha.png", height = 2000, width = 3500, res = 250, bg="transparent")
p1 <- DotPlot(Alpha, features = top5$gene, cols = c("#1f77b4","#d62728"), dot.scale = 30)+ RotatedAxis() + labs(title = NULL, x=NULL, y = "Subclusters") + theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 30), 
                                                       axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 30), legend.title = element_text( size = 30), legend.text = element_text(size=20), axis.line = element_line(size=1.5))
p1
dev.off()
p2 <- DotPlot(Alpha, features = top10$gene, split.by = "orig.ident", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)


Idents(Alpha) <- Alpha$orig.ident
t.markers <- FindAllMarkers(Alpha, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(t.markers, 'DEG/Alpha_DEG_diet.csv')
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

DoHeatmap(subset(Alpha, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

#Marker Violinplot - Alpha/Gcg
png("Violinplot_Gcg.png", height = 2000, width = 2000, res = 300)
VlnPlot(Alpha, features = "Gcg", cols = c("darkslateblue", "cadetblue2")) + NoLegend() & theme(text = element_text(size = 30), axis.line = element_line(size=1.5), axis.text = element_text(size = 30), 
                                                                                                      legend.text = element_text(size = 30))

dev.off()

p1 <- DotPlot(Alpha, features = top5$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between diet')
# VlnPlot(Alpha, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident")+ theme(legend.position = 'right')
p2 <- DotPlot(Alpha, features = top5$gene, split.by = "seurat_clusters", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)

png("figures/Alpha_heatmap2.png", width = 2000, height = 1000, res = 200, bg="transparent")

DoHeatmap(subset(Alpha, downsample = 200), features = top5$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

dev.off()

p1 <-DimPlot(Alpha, group.by = "celltypes", pt.size = 2, cols = c("#fdbf6f")) +labs(title = "Alpha cells") + NoLegend() 
p2 <-DimPlot(Alpha, group.by = "seurat_clusters", pt.size = 2, cols = colmaps2)
p3<- FeaturePlot(Alpha, features = c("Gcg"),pt.size = 1)& 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))
p4 <- VlnPlot(Alpha, features = "Gcg", cols = colmaps2)
plot_grid(p1,p2,p3,p4, ncol = 2)


saveRDS(Alpha, "Alpha.rds")



#### Delta cells ####
delt <- readRDS("Delta.rds")
Delta <- subset(delt)

Delta <- NormalizeData(Delta)
Delta <- FindVariableFeatures(Delta, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Delta)
Delta <- ScaleData(Delta, features = all.genes)
Delta <- RunPCA(Delta, features = VariableFeatures(object = Delta))
# intrinsicDimension::maxLikGlobalDimEst(Delta@reductions$pca@cell.embeddings,
#                                        k = 10)
#10 PC

# #Elbow plot
# ElbowPlot(Delta, ndims = 50)

dims.use = 1:19
#Umap
Delta <- RunUMAP(Delta, dims = dims.use, verbose = F)
#FindNeighbors
Delta <- FindNeighbors(Delta, dims = dims.use, verbose = F)
#Findclusters
Delta <- FindClusters(Delta, resolution = 0.3, verbose = F)

# # #Final Umap
# DimPlot(Delta, label = TRUE)
# # 
# # 
# #Compute distance matrix to UMAP coordinates
# distance_matrix <- dist(Embeddings(Delta[['umap']])[, 1:2])
# 
# #Isolate cluster nameq
# clusters <- Delta@active.ident
# 
# #Compute silhouette score for each cluster
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# Delta@meta.data$silhouette_score <- silhouette[,3]
# 
# #Plot
# fviz_silhouette(silhouette)
# # 0.63

DimPlot(Delta, group.by = "celltypes", pt.size = 8, cols = c("#e31a1c")) +labs(title = "Delta cells") + NoLegend() +
  DimPlot(Delta, group.by = "orig.ident", pt.size = 8, cols = c("#1f77b4","#d62728")) +
  DimPlot(Delta, group.by = "seurat_clusters", pt.size = 8, cols = colmaps2) )

png("figures/combined/Delta.png", width = 2800, height = 1000, res = 200)

#UMAPs for JDRF - Delta
p1 <- DimPlot(Delta, group.by = "celltypes", pt.size = 3, cols = "#e31a1c") + labs(title = "Delta cells")
png("UMAP_Delta_Diet.png", width = 2700, height = 2700, res = 250, bg="transparent")
p2 <- DimPlot(Delta, group.by = "orig.ident", pt.size = 6.5, cols = c("#B19CD9", "#bfe6ff")) + labs(title = "Diet") & theme(title = element_text(size=55),
                                                                                                                             legend.position = c(0.6,1.06),
                                                                                                                          axis.text.x=element_text(hjust=1, size=55),
                                                                                                                          axis.title = element_text(size=55),
                                                                                                                          axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                          axis.ticks = element_line(size = 3),
                                                                                                                          axis.ticks.length = unit(7, "pt"),
                                                                                                                          legend.text=element_text(size=50),
                                                                                                                          legend.title=element_text(size=55),
                                                                                                                          legend.key = element_rect(size = 100),
                                                                                                                          axis.line = element_line(size=1.5),legend.background = element_rect(fill = "transparent"),
panel.background = element_rect(fill = "transparent"),
plot.background = element_rect(fill = "transparent")) & guides(color = guide_legend(override.aes = list(size = 10) ) )
p2
dev.off()
png("UMAP_Delta_SeuratCluster.png", width = 2700, height = 2700, res = 250, bg="transparent")
p3 <- DimPlot(Delta, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) + labs(title = "Subclusters") &  theme(title = element_text(size=55),
                                                                                                                         legend.position = c(0.85,1.03),
                                                                                                                       axis.text.x=element_text(hjust=1, size=55),
                                                                                                                       axis.title = element_text(size=55),
                                                                                                                       axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                       axis.ticks = element_line(size = 3),
                                                                                                                       axis.ticks.length = unit(7, "pt"),
                                                                                                                       legend.text=element_text(size=55),
                                                                                                                       legend.title=element_text(size=55),
                                                                                                                       legend.key = element_rect(size = 100),
                                                                                                                       axis.line = element_line(size=1.5),legend.background = element_rect(fill = "transparent"),
                                                                                                                       panel.background = element_rect(fill = "transparent"),
                                                                                                                       plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
p3
dev.off()
plot_grid(p3,p2)

Idents(Delta) <- Delta$seurat_clusters
t.markers <- FindAllMarkers(Delta, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top5
write.csv(t.markers, 'DEG/Delta_DEG_subclusters.csv')

DoHeatmap(subset(Delta, downsample = 200), features = top10$gene, group.colors = colmaps2) + 
  scale_fill_gradientn(colours = rev(mapal)) +labs(title = 'DE genes between subclusters')

VlnPlot(Delta, features = top10$gene, cols = colmaps2, ncol = 5) + theme(legend.position = 'right')
VlnPlot(Delta, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident", ncol = 5)+ 
  theme(legend.position = 'right')
#p1 = Dotplot needed for JDRF - Delta
png("DEG_dotplot_Delta.png", height = 2000, width = 3500, res = 250, bg="transparent")
p1 <- DotPlot(Delta, features = top5$gene, cols = c("#1f77b4","#d62728"), dot.scale = 28)+ RotatedAxis()+
  labs(title = NULL, x=NULL, y="Subclusters") & theme(title = element_text(hjust = 0.5, size=40), axis.text.x = element_text(size = 40), axis.title.x = element_text(size = 40), 
                                                     axis.text.y = element_text(size = 40), axis.title.y = element_text(size = 40), legend.title = element_text( size = 30), legend.text = element_text(size=20), axis.line = element_line(size=1.5),legend.background = element_rect(fill = "transparent"),
                                                     panel.background = element_rect(fill = "transparent"),
                                                     plot.background = element_rect(fill = "transparent"))
p1
dev.off()
p2 <- DotPlot(Delta, features = top10$gene, split.by = "orig.ident", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)

##ViolinPlots for JDRF - Delta
Idents(Delta) <- Delta$orig.ident
Features_DeltaCells_HAMSAB <- c("Sst", "Cd74","Insm1", "Sod1", "Gbp7", "Dnajb1", "Fos","Ins2")
png("VlnPlot_Delta.png", height = 1800, width = 6500, res = 250)
VlnPlot(Delta, features = Features_DeltaCells_HAMSAB, cols = c("#B19CD9", "#bfe6ff"), ncol = 8, y.max = 8) + theme(legend.position = "right") & 
  theme(text = element_text(size = 30), axis.ticks.y = element_line(size = 3), axis.ticks.length.y = unit(7, "pt"),axis.line = element_line(size=1.5), axis.text = element_text(size = 30), axis.title = element_text(size = 35)) 
dev.off()

Idents(Delta) <- Delta$orig.ident
t.markers <- FindAllMarkers(Delta, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(t.markers, 'DEG/Delta_DEG_diet.csv')
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

DoHeatmap(subset(Delta, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

#Violin plot for marker - Delta/sst
png("Violinplot_Sst.png", height = 2000, width = 2000, res = 300)
VlnPlot(Delta, features = "Sst", cols = c("darkslateblue", "cadetblue2")) + NoLegend() + theme(text = element_text(size = 30), axis.line = element_line(size=1.5), axis.text = element_text(size = 30), 
                                                                                                                   legend.text = element_text(size = 30))
dev.off()

p1 <- DotPlot(Delta, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between diet')
p2 <- DotPlot(Delta, features = top10$gene, split.by = "seurat_clusters", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)

png("figures/combined/Delta_heatmap2.png", width = 2000, height = 1000, res = 200)

DoHeatmap(subset(Delta, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

dev.off()


p1<-FeaturePlot(Delta, features = c("Sst"),pt.size = 1)& 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))
p2<-VlnPlot(Delta, features = c("Sst"), cols = colmaps2)
plot_grid(p1, p2, ncol=2)

DimPlot(Delta, group.by = "celltypes", pt.size = 5, cols = c("#e31a1c")) +labs(title = "Delta cells") + NoLegend() +
  DimPlot(Delta, group.by = "seurat_clusters", pt.size = 5, cols = colmaps2) 

p1 <-DimPlot(Delta, group.by = "celltypes", pt.size = 2, cols = c("#e31a1c")) +labs(title = "Delta cells") + NoLegend() 
p2 <-DimPlot(Delta, group.by = "seurat_clusters", pt.size = 2, cols = colmaps2)
p3<- FeaturePlot(Delta, features = c("Sst"),pt.size = 1)& 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))
p4 <- VlnPlot(Delta, features = "Sst", cols = colmaps2)
plot_grid(p1,p2,p3,p4, ncol = 2)

saveRDS(Delta, "Delta.rds")


#### Beta cells ####
bet <- readRDS("beta.rds")
Beta <- subset(bet)

Beta <- NormalizeData(Beta)
Beta <- FindVariableFeatures(Beta, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Beta)
Beta <- ScaleData(Beta, features = all.genes)
Beta <- RunPCA(Beta, features = VariableFeatures(object = Beta))
# intrinsicDimension::maxLikGlobalDimEst(Beta@reductions$pca@cell.embeddings,
#                                        k = 10)
#10 PC

# #Elbow plot
# ElbowPlot(Beta, ndims = 50)

dims.use = 1:20
#Umap
Beta <- RunUMAP(Beta, dims = dims.use, verbose = F)
#FindNeighbors
Beta <- FindNeighbors(Beta, dims = dims.use, verbose = F)
#Findclusters
Beta <- FindClusters(Beta, resolution = 0.2, verbose = F)

# # #Final Umap
# DimPlot(Beta, label = TRUE)
# # 
# # 
# #Compute distance matrix to UMAP coordinates
# distance_matrix <- dist(Embeddings(Beta[['umap']])[, 1:2])
# 
# #Isolate cluster nameq
# clusters <- Beta@active.ident
# 
# #Compute silhouette score for each cluster
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# Beta@meta.data$silhouette_score <- silhouette[,3]
# 
# #Plot
# fviz_silhouette(silhouette)
# # 0.63

DimPlot(Beta, group.by = "celltypes", pt.size = 8, cols = c("#a6cee3")) +labs(title = "Beta cells") + NoLegend() +
  DimPlot(Beta, group.by = "orig.ident", pt.size = 8, cols = c("#1f77b4","#d62728")) +
  DimPlot(Beta, group.by = "seurat_clusters", pt.size = 8, cols = colmaps2) 

png("figures/combined/Beta.png", width = 3500, height = 2700, res = 200)

p1 <- DimPlot(Beta, group.by = "celltypes", pt.size = 3, cols = "#a6cee3") + labs(title = "Beta cells")

png("UMAP_beta_Diet.png", width = 2700, height = 2700, res = 250, bg = "transparent")
DimPlot(Beta, group.by = "orig.ident", pt.size = 6.5, cols = c("#B19CD9", "#bfe6ff")) + labs(title = "Diet") &  theme(title = element_text(size=55),
                                                                                                                            legend.position = c(0.60,1.06),
                                                                                                                           axis.text.x=element_text(hjust=1, size=55),
                                                                                                                           axis.title = element_text(size=55),
                                                                                                                           axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                           axis.ticks = element_line(size = 3),
                                                                                                                           axis.ticks.length = unit(7, "pt"),
                                                                                                                           legend.text=element_text(size=50),
                                                                                                                           legend.title=element_text(size=50),
                                                                                                                           legend.key = element_rect(size = 100),
                                                                                                                           axis.line = element_line(size=1.5),
                                                                                                                      legend.background = element_rect(fill = "transparent"),
                                                                                                                      panel.background = element_rect(fill = "transparent"),
                                                                                                                      plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
dev.off()
png("UMAP_beta_seuratclusters.png",width = 2700, height = 2700, res = 250, bg = "transparent")
DimPlot(Beta, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) + labs(title = "Subclusters") & theme(title = element_text(size=55), 
                                                                                                                        legend.position = c(0.8, 1),
                                                                                                                        axis.text.x=element_text(hjust=1, size=55),
                                                                                                                        axis.title = element_text(size=55),
                                                                                                                        axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                        axis.ticks = element_line(size = 3),
                                                                                                                        axis.ticks.length = unit(7, "pt"),
                                                                                                                        legend.text=element_text(size=55),
                                                                                                                        legend.title=element_text(size=55),
                                                                                                                        legend.key = element_rect(size = 100),
                                                                                                                        axis.line = element_line(size=1.5),
                                                                                                                        legend.background = element_rect(fill = "transparent"),
                                                                                                                        panel.background = element_rect(fill = "transparent"),
                                                                                                                        plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
dev.off()

plot_grid(p3,p2, ncol = 2)

Idents(Beta) <- Beta$seurat_clusters
t.markers <- FindAllMarkers(Beta, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

t.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
write.csv(t.markers, 'DEG/Beta_DEG_subclusters.csv')

DoHeatmap(subset(Beta, downsample = 200), features = top10$gene, group.colors = colmaps2) + 
  scale_fill_gradientn(colours = rev(mapal)) +labs(title = 'DE genes between subclusters')

Features_Betacell_HAMSAB <- c("Ins1", "Ins2","Iapp", "G6pc2", "Ssr3", "Fos","Impact", "Dnajb1")

Idents(Beta) <- Beta$orig.ident
png("VlnPlot_Beta.png", height = 1800, width = 6500, res = 250, bg="transparent")
VlnPlot(Beta, features = Features_Betacell_HAMSAB, cols = c("#B19CD9", "#bfe6ff"), ncol = 8, y.max = 8) + theme(legend.position = 'right') & 
  theme(text = element_text(size = 30), axis.ticks.y = element_line(size = 3), axis.ticks.length.y = unit(7, "pt"),
        axis.line = element_line(size=1.5), axis.text = element_text(size = 30), axis.title = element_text(size = 35),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent")) 

dev.off()
png("VlnPlot_Beta_cluster1.png", height = 1700, width = 3500, res = 250)
VlnPlot(Beta, features = Features_Betacell_cluster1, cols = c("darkslateblue", "cadetblue2"), split.by = "orig.ident", ncol = 5, idents = 1) + theme(legend.position = 'right') & 
  theme(text = element_text(size = 30), axis.line.y = element_line(size=1.5), axis.text.y = element_text(size = 30), axis.title = element_text(size = 35))
dev.off()
VlnPlot(Beta, features = top10$gene, cols = c("#1f77b4","#d62728"), split.by = "orig.ident", ncol = 5)+ 
  theme(legend.position = 'right')

png("DEG_dotplot_Beta.png", height = 2000, width = 3500, res = 250, bg="transparent")
p1 <- DotPlot(Beta, features = top5$gene, cols = c("#1f77b4","#d62728"), dot.scale = 30, scale.max = 100)+ RotatedAxis()+
  labs(title = NULL, x=NULL, y="Subclusters")+ theme(axis.text.x = element_text(size = 30), axis.title.x = element_text(size = 30), 
                                                                               axis.text.y = element_text(size = 30), axis.title.y = element_text(size = 30), legend.title = element_text( size = 30), legend.text = element_text(size=20), axis.line = element_line(size=1.5))
p1
dev.off()

p2 <- DotPlot(Beta, features = top10$gene, split.by = "orig.ident", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)


Idents(Beta) <- Beta$orig.ident
t.markers <- FindAllMarkers(Beta, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(t.markers, 'DEG/Beta_DEG_diet.csv')
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
t.markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = avg_log2FC) -> top4

DoHeatmap(subset(Beta, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

png("Violinplot_Ins.png", height = 2000, width = 2000, res = 300)
VlnPlot(Beta, features = c("Ins1", "Ins2"), cols = c("darkslateblue", "cadetblue2")) + NoLegend() & theme(text = element_text(size = 30), axis.line = element_line(size=1.5), axis.text = element_text(size = 30), 
                                                                                                legend.text = element_text(size = 30))
dev.off()

p1 <- DotPlot(Beta, features = top10$gene, cols = "RdBu", dot.scale = 15)+ RotatedAxis()+
  labs(title = 'DE genes between diet')
p2 <- DotPlot(Beta, features = top10$gene, split.by = "seurat_clusters", cols = "RdBu", dot.scale = 15)+ 
  RotatedAxis()
plot_grid(p1,p2, ncol = 1)

png("figures/combined/Beta_heatmap2.png", width = 2000, height = 1000, res = 200)

DoHeatmap(subset(Beta, downsample = 200), features = top10$gene, group.colors = colmaps2)+ 
  labs(title = 'DE genes between diet') + scale_fill_gradientn(colours = rev(mapal))

dev.off()

FeaturePlot(Beta, features = c("Ins1","Ins2"),pt.size = 1, ncol = 1)& 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))
VlnPlot(Beta, features = c("Ins1","Ins2"), cols = colmaps2)

DimPlot(Beta, group.by = "celltypes", pt.size = 1, cols = c("#a6cee3")) +labs(title = "Beta cells") + NoLegend() +
  DimPlot(Beta, group.by = "seurat_clusters", pt.size = 1, cols = colmaps2) 

p1 <-DimPlot(Beta, group.by = "celltypes", pt.size = 2, cols = c("#a6cee3")) +labs(title = "Beta cells") + NoLegend() 
p2 <-DimPlot(Beta, group.by = "seurat_clusters", pt.size = 2, cols = colmaps2)
p3<- FeaturePlot(Beta, features = c("Ins1",'Ins2'),pt.size = 1)& 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))
p4 <- VlnPlot(Beta, features = c("Ins1",'Ins2'), cols = colmaps2)
plot_grid(p3,p4, ncol = 1)

saveRDS(Beta, "Beta.rds")

#### mix ####
mix <- readRDS("Alpha_Delta.rds")
mix <- subset(rds, ident = c("Alpha/Delta"))
DimPlot(mix, cols = c("#ffd92f"), pt.size = 3)
FeaturePlot(mix, features = c("Sst","Ins1","Gcg"),pt.size = 2, ncol = 3)& 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu")) 

mix <- NormalizeData(mix)
mix <- FindVariableFeatures(mix, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mix)
mix <- ScaleData(mix, features = all.genes)
mix <- RunPCA(mix, features = VariableFeatures(object = mix))
intrinsicDimension::maxLikGlobalDimEst(mix@reductions$pca@cell.embeddings,
                                       k = 10)
#14 PC

# #Elbow plot
# ElbowPlot(mix, ndims = 50)

dims.use = 1:14
#Umap
mix <- RunUMAP(mix, dims = dims.use, verbose = F)
#FindNeighbors
mix <- FindNeighbors(mix, dims = dims.use, verbose = F)
#Findclusters
mix <- FindClusters(mix, resolution = 0.4, verbose = F)

# # #Final Umap
# DimPlot(mix, label = TRUE)
# # 
# # 
# #Compute distance matrix to UMAP coordinates
# distance_matrix <- dist(Embeddings(mix[['umap']])[, 1:2])
# 
# #Isolate cluster nameq
# clusters <- mix@active.ident
# 
# #Compute silhouette score for each cluster
# silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
# mix@meta.data$silhouette_score <- silhouette[,3]
# 
# #Plot
# fviz_silhouette(silhouette)
# # 0.63

png("UMAP_Alpha_Delta_seuratclusters.png",width = 2700, height = 2700, res = 250)
p3 <- DimPlot(mix, group.by = "seurat_clusters", pt.size = 6.5, cols = colmaps2) + labs(title = "Subclusters")+  theme(title = element_text(size=55), 
                                                                                                                        legend.position = c(0.9, 1),
                                                                                                                        axis.text.x=element_text(hjust=1, size=55),
                                                                                                                        axis.title = element_text(size=55),
                                                                                                                        axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                        axis.ticks = element_line(size = 3),
                                                                                                                        axis.ticks.length = unit(7, "pt"),
                                                                                                                        legend.text=element_text(size=55),
                                                                                                                        legend.title=element_text(size=55),
                                                                                                                        legend.key = element_rect(size = 100),
                                                                                                                        axis.line = element_line(size=1.5)) +
  guides(color = guide_legend(override.aes = list(size = 10) ) )
p3
dev.off()

png("UMAP_AlphaDelta_Diet.png", width = 2700, height = 2700, res = 250)
p2 <- DimPlot(mix, group.by = "orig.ident", pt.size = 6.5, cols = c("#B19CD9", "#bfe6ff")) + labs(title = "Diet") +  theme(title = element_text(size=55),
                                                                                                                            legend.position = c(0.60,1.06),
                                                                                                                            axis.text.x=element_text(hjust=1, size=55),
                                                                                                                            axis.title = element_text(size=55),
                                                                                                                            axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                                            axis.ticks = element_line(size = 3),
                                                                                                                            axis.ticks.length = unit(7, "pt"),
                                                                                                                            legend.text=element_text(size=50),
                                                                                                                            legend.title=element_text(size=50),
                                                                                                                            legend.key = element_rect(size = 100),
                                                                                                                            axis.line = element_line(size=1.5)) +
  guides(color = guide_legend(override.aes = list(size = 10) ) )
p2
dev.off()

DimPlot(mix, group.by = "celltypes", pt.size = 8, cols = c("#ffd92f")) +labs(title = "Alpha/Delta") + NoLegend() +
  DimPlot(mix, group.by = "orig.ident", pt.size = 8, cols = c("#1f77b4","#d62728")) +
  DimPlot(mix, group.by = "seurat_clusters", pt.size = 8, cols = colmaps2) 

png("figures/combined/alpha_delta.png", width = 2800, height = 1000, res = 200)

p1 <- DimPlot(mix, group.by = "celltypes", pt.size = 3, cols = "#ffd92f") + labs(title = "Alpha/Delta")
p2 <- DimPlot(mix, group.by = "orig.ident", pt.size = 3, cols = c("#1f77b4","#d62728")) + labs(title = "Diet") 
p3 <- DimPlot(mix, group.by = "seurat_clusters", pt.size = 3, cols = colmaps2) + labs(title = "Subclusters")

plot_grid(p1,p2,p3, ncol = 3)

dev.off()

t.markers <- FindAllMarkers(mix, min.pct = 0, logfc.threshold = 0.25)


FeaturePlot(mix, features = c("Sst","Ins1","Gcg"),pt.size = 6, ncol = 3, coord.fixed = T, order = T)& 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu")) 

FeaturePlot(combined, features = c("nCount_RNA","nFeature_RNA"),pt.size = 1, ncol = 3)& 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu")) 

plot<-DimPlot(combined, cols = colmaps, label = T, repel = T)
plot2 <-FeaturePlot(combined, features = c("nCount_RNA","nFeature_RNA"),pt.size = 1, ncol = 2)& 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu")) 
plot_grid(plot, plot2, ncol = 1)

VlnPlot(combined, features = c("nCount_RNA","nFeature_RNA"))





#### cell cluster percentage - 100% is clusters####
Alpha$cell_cluster <- paste(Alpha$seurat_clusters, Alpha$celltypes,sep = '_')
alpha_table1 <- prop.table(table(Alpha$cell_cluster, Alpha$orig.ident), margin = 2)
Delta$cell_cluster <- paste(Delta$celltypes, Delta$seurat_clusters, sep = '_')
Delta_table1 <- prop.table(table(Delta$cell_cluster, Delta$orig.ident), margin = 2)
Beta$cell_cluster <- paste(Beta$celltypes, Beta$seurat_clusters, sep = '_')
Beta_table1 <- prop.table(table(Beta$cell_cluster, Beta$orig.ident), margin = 2)
tcells$cell_cluster <- paste(tcells$seurat_clusters, tcells$celltypes,sep = '_')
tcells_table <- prop.table(table(tcells$cell_cluster,tcells$orig.ident), margin = 2)

cell_table <- rbind(alpha_table1, Delta_table1, Beta_table1)*100
cell_table
write.csv(cell_table, "Percentages_clusters.csv")

tcells$cell_cluster <- paste( tcells$seurat_clusters, tcells$celltypes,sep = '_')
tcells_table1 <- prop.table(table( tcells$cell_cluster, tcells$orig.ident), margin = 2)
bcells$cell_cluster <- paste(bcells$celltypes, bcells$seurat_clusters, sep = '_')
bcells_table1 <- prop.table(table(bcells$cell_cluster, bcells$orig.ident), margin = 2)
cdc$cell_cluster <- paste(cdc$celltypes, cdc$seurat_clusters, sep = '_')
cdc_table1 <- prop.table(table(cdc$cell_cluster, cdc$orig.ident), margin = 2)
mf$cell_cluster <- paste(mf$seurat_clusters, mf$celltypes,sep = '_')
mf_table <- prop.table(table(mf$cell_cluster, mf$orig.ident), margin = 2)
pdc$cell_cluster <- paste(pdc$celltypes, pdc$seurat_clusters, sep = '_')
pdc_table1 <- prop.table(table(pdc$cell_cluster, pdc$orig.ident), margin = 2)

cell_table <- rbind(tcells_table1, bcells_table1, cdc_table1, mf_table, pdc_table1)*100
cell_table
write.csv(cell_table, "Percentages_Immune_clusters.csv")

#### cell cluster percentage - 100% is diet####
Alpha$cell_cluster <- paste(Alpha$celltypes, Alpha$seurat_clusters, sep = '_')
alpha_table1 <- prop.table(table(Alpha$orig.ident, Alpha$cell_cluster), margin = 2)
Delta$cell_cluster <- paste(Delta$celltypes, Delta$seurat_clusters, sep = '_')
Delta_table1 <- prop.table(table(Delta$orig.ident, Delta$cell_cluster), margin = 2)
Beta$cell_cluster <- paste(Beta$celltypes, Beta$seurat_clusters, sep = '_')
Beta_table1 <- prop.table(table(Beta$orig.ident, Beta$cell_cluster), margin = 2)

cell_table <- cbind(alpha_table1,Delta_table1, Beta_table1)*100
cell_table
write.csv(cell_table, "Percentages_Endocrine_clusters_diet.csv")

tcells$cell_cluster <- paste(tcells$celltypes, tcells$seurat_clusters, sep = '_')
t_table1 <- prop.table(table(tcells$orig.ident, tcells$cell_cluster), margin = 2)
bcells$cell_cluster <- paste(bcells$celltypes, bcells$seurat_clusters, sep='_')
b_table1 <- prop.table(table(bcells$orig.ident, bcells$cell_cluster), margin = 2)
cdc$cell_cluster <- paste(cdc$celltypes, cdc$seurat_clusters, sep = '_')
cdc_table1 <- prop.table(table(cdc$orig.ident, cdc$cell_cluster), margin = 2)
mf$cell_cluster <- paste(mf$celltypes, mf$seurat_clusters, sep = '_')
mf_table1 <- prop.table(table(mf$orig.ident, mf$cell_cluster), margin = 2)
pdc$cell_cluster <- paste(pdc$celltypes, pdc$seurat_clusters, sep = '_')
pdc_table1 <- prop.table(table(pdc$orig.ident, pdc$cell_cluster), margin = 2)

table <- cbind(t_table1, b_table1, cdc_table1, mf_table1, pdc_table1)*100
table
write.csv(table, "Percentages_Immune_clusters_diet.csv")

#### cell cluster percentage - all cells ####
combined$cell_cluster <- paste(combined$celltypes, sep = '_')
combined_table1 <- prop.table(table(combined$orig.ident, combined$cell_cluster), margin = 2)
combined_table2 <- combined_table1*100
write.csv(combined_table2, "Percentages_clusters_allcells.csv")

#### cell cluster percentage - per diet- all cells ####
combined$cell_cluster <- paste(combined$celltypes, sep = '_')
combined_table1 <- prop.table(table( combined$cell_cluster, combined$orig.ident), margin = 2)
combined_table2 <- combined_table1*100
write.csv(combined_table2, "Percentages_clusters_allcells_diet.csv")


