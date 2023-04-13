set.seed(1732)

setwd("~/Desktop/Data for R/filtered_feature_bc_matrix")

library(ggplot2)
library(Seurat)

#all cells ####
#finding the cell cycle effect###
Seuratobj <- readRDS("combined_doubletremoved_finetuned.rds")

load("cycle.rda")

rownames(Seuratobj@assays$RNA@data) <- toupper(rownames(Seuratobj@assays$RNA@data))

Seuratobj<- CellCycleScoring(Seuratobj, 
                                 g2m.features=g2m_genes, 
                                 s.features=s_genes)

head(Seuratobj@meta.data)

Seuratobj <- FindVariableFeatures(Seuratobj)

Seuratobj <- ScaleData(Seuratobj)

Seuratobj <- RunPCA(Seuratobj, 
                    features = VariableFeatures(Seuratobj), 
                    nfeatures.print = 10)

DimPlot(object = Seuratobj, reduction = "pca", group.by = "Phase", split.by = "Phase")


#Imaging Cell cycle effect###
library(ggpubr)
Seuratobj <- RunTSNE(Seuratobj, dims = 1:20)

png(filename="UMAP_CellCycle_All_order_2.png", height = 2800, width = 4000,res=300,bg="transparent")
DimPlot(object = Seuratobj, group.by = "Phase",pt.size = 1, label.size = 10, order = c("S", "G2M", "G1")) &  theme(title = element_text(size = 30),
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

Seuratobj <- FindNeighbors(Seuratobj, dims = 1:20)

Seuratobj <- FindClusters(Seuratobj, resolution = 0.6)

p2 <- DimPlot(Seuratobj, label = TRUE) + NoLegend()
ggarrange(p2)

Seuratobj$cell_cluster <- paste(Seuratobj$celltypes, sep = '_')
Seuratobj_table1 <- prop.table(table(Seuratobj$Phase, Seuratobj$cell_cluster), margin = 2)
Seuratobj_table2 <- Seuratobj_table1*100
write.csv(Seuratobj_table2, "Percentages_clusters_allcells_cellcycle.csv")


#T-Cells####
#finding the cell cycle effect###
Seuratobj_T <- readRDS("tcells.rds")

load("cycle.rda")

rownames(Seuratobj_T@assays$RNA@data) <- toupper(rownames(Seuratobj_T@assays$RNA@data))

Seuratobj_T<- CellCycleScoring(Seuratobj_T, 
                             g2m.features=g2m_genes, 
                             s.features=s_genes)

head(Seuratobj_T@meta.data)

Seuratobj_T <- FindVariableFeatures(Seuratobj_T)

Seuratobj_T <- ScaleData(Seuratobj_T)

Seuratobj_T <- RunPCA(Seuratobj_T, 
                    features = VariableFeatures(Seuratobj_T), 
                    nfeatures.print = 10)

DimPlot(object = Seuratobj_T, reduction = "pca", group.by = "Phase", split.by = "orig.ident")


#Imaging Cell cycle effect####
Seuratobj_T <- RunTSNE(Seuratobj_T, dims = 1:21)

png("UMAP_Tcells_cellcylce.png", width = 2700, height = 2700, res = 250, bg="transparent")
DimPlot(object = Seuratobj_T, group.by = "Phase", pt.size = 6) & NoLegend() &  theme(axis.text.x=element_text(hjust=1, size=55),
                                                                axis.title = element_text(size=55),
                                                                axis.ticks = element_line(size = 3),
                                                                axis.ticks.length = unit(7, "pt"),
                                                                axis.text.y = element_text(hjust=0.5,size = 55),
                                                                axis.line = element_line(size=1.5), 
                                                                legend.background = element_rect(fill = "transparent"),
                                                                panel.background = element_rect(fill = "transparent"),
                                                                plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
dev.off()

Seuratobj_T <- FindNeighbors(Seuratobj_T, dims = 1:21)

Seuratobj_T <- FindClusters(Seuratobj_T, resolution = 0.4)

p2 <- DimPlot(Seuratobj_T, label = TRUE) + NoLegend()
ggarrange(p1,p2)

#Percentage####
Seuratobj_T$cell_cluster <- paste(Seuratobj_T$celltypes, sep = '_')
Seuratobj.T_table1 <- prop.table(table(Seuratobj_T$Phase, Seuratobj_T$orig.ident), margin = 2)
Seuratobj.T_table2 <- Seuratobj.T_table1*100
write.csv(Seuratobj.T_table2, "Percentages_clusters_Tcells_cellcycle.csv")

#B-Cells####
#finding the cell cycle effect###
Seuratobj_B <- readRDS("bcells.rds")

load("cycle.rda")

rownames(Seuratobj_B@assays$RNA@data) <- toupper(rownames(Seuratobj_B@assays$RNA@data))

Seuratobj_B<- CellCycleScoring(Seuratobj_B, 
                                  g2m.features=g2m_genes, 
                                  s.features=s_genes)

head(Seuratobj_B@meta.data)

Seuratobj_B <- FindVariableFeatures(Seuratobj_B)

Seuratobj_B <- ScaleData(Seuratobj_B)

Seuratobj_B <- RunPCA(Seuratobj_B, 
                         features = VariableFeatures(Seuratobj_B), 
                         nfeatures.print = 10)

DimPlot(object = Seuratobj_B, reduction = "pca", group.by = "Phase", split.by = "orig.ident")


#Imaging Cell cycle effect####
Seuratobj_B <- RunTSNE(Seuratobj_B, dims = 1:19)

png("UMAP_Bcells_cellcylce.png", width = 2700, height = 2700, res = 250, bg="transparent")
DimPlot(object = Seuratobj_B, group.by = "Phase", pt.size = 6) & NoLegend() &  theme(axis.text.x=element_text(hjust=1, size=55),
                                                                                         axis.title = element_text(size=55),
                                                                                         axis.ticks = element_line(size = 3),
                                                                                         axis.ticks.length = unit(7, "pt"),
                                                                                         axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                         axis.line = element_line(size=1.5),
                                                                        legend.background = element_rect(fill = "transparent"),
                                                                        panel.background = element_rect(fill = "transparent"),
                                                                        plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
dev.off()

Seuratobj_B <- FindNeighbors(Seuratobj_B, dims = 1:19)

Seuratobj_B <- FindClusters(Seuratobj_B, resolution = 0.4)

p2 <- DimPlot(Seuratobj_B, label = TRUE) + NoLegend()
ggarrange(p1,p2)

#Percentage####
Seuratobj_B$cell_cluster <- paste(Seuratobj_B$celltypes, sep = '_')
Seuratobj.B_table1 <- prop.table(table(Seuratobj_B$Phase, Seuratobj_B$orig.ident), margin = 2)
Seuratobj.B_table2 <- Seuratobj.B_table1*100
write.csv(Seuratobj.B_table2, "Percentages_clusters_Bcells_cellcycle.csv")

# Macrophages ####
#finding the cell cycle effect###
Seuratobj_mf <- readRDS("mf.rds")

load("cycle.rda")

rownames(Seuratobj_mf@assays$RNA@data) <- toupper(rownames(Seuratobj_mf@assays$RNA@data))

Seuratobj_mf<- CellCycleScoring(Seuratobj_mf, 
                               g2m.features=g2m_genes, 
                               s.features=s_genes)

head(Seuratobj_mf@meta.data)

Seuratobj_mf <- FindVariableFeatures(Seuratobj_mf)

Seuratobj_mf <- ScaleData(Seuratobj_mf)

Seuratobj_mf <- RunPCA(Seuratobj_mf, 
                      features = VariableFeatures(Seuratobj_mf), 
                      nfeatures.print = 10)

DimPlot(object = Seuratobj_mf, reduction = "pca", group.by = "Phase", split.by = "orig.ident")


#Imaging Cell cycle effect####
Seuratobj_mf <- RunTSNE(Seuratobj_mf, dims = 1:19)

png("UMAP_mfcells_cellcylce.png", width = 2700, height = 2700, res = 250, bg="transparent")
DimPlot(object = Seuratobj_mf, group.by = "Phase", pt.size = 6) & NoLegend() &  theme(axis.text.x=element_text(hjust=1, size=55),
                                                                                                    axis.title = element_text(size=55),
                                                                                                    axis.ticks = element_line(size = 3),
                                                                                                    axis.ticks.length = unit(7, "pt"),
                                                                                                    axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                    axis.line = element_line(size=1.5),
                                                                                                    legend.background = element_rect(fill = "transparent"),
                                                                                                    panel.background = element_rect(fill = "transparent"),
                                                                                                    plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
dev.off()

Seuratobj_mf <- FindNeighbors(Seuratobj_mf, dims = 1:19)

Seuratobj_mf <- FindClusters(Seuratobj_mf, resolution = 0.4)

p2 <- DimPlot(Seuratobj_mf, label = TRUE) + NoLegend()
ggarrange(p1,p2)

#Percentage####
Seuratobj_mf$cell_cluster <- paste(Seuratobj_mf$celltypes, sep = '_')
Seuratobj.mf_table1 <- prop.table(table(Seuratobj_mf$Phase, Seuratobj_mf$orig.ident), margin = 2)
Seuratobj.mf_table2 <- Seuratobj.mf_table1*100
write.csv(Seuratobj.mf_table2, "Percentages_clusters_mf_cellcycle.csv")

# CD4+ ####
#finding the cell cycle effect###
tcells <- readRDS("tcells.rds")
Seuratobj_cd4 <- subset(tcells, ident="1")

load("cycle.rda")

rownames(Seuratobj_cd4@assays$RNA@data) <- toupper(rownames(Seuratobj_cd4@assays$RNA@data))

Seuratobj_cd4<- CellCycleScoring(Seuratobj_cd4, 
                                g2m.features=g2m_genes, 
                                s.features=s_genes)

head(Seuratobj_cd4@meta.data)

Seuratobj_cd4 <- FindVariableFeatures(Seuratobj_cd4)

Seuratobj_cd4 <- ScaleData(Seuratobj_cd4)

Seuratobj_cd4 <- RunPCA(Seuratobj_cd4, 
                       features = VariableFeatures(Seuratobj_cd4), 
                       nfeatures.print = 10)

DimPlot(object = Seuratobj_cd4, reduction = "pca", group.by = "Phase", split.by = "orig.ident")


#Imaging Cell cycle effect####
Seuratobj_cd4 <- RunUMAP(Seuratobj_cd4, dims = 1:5)

png("UMAP_cd4cells_cellcylce.png", width = 2700, height = 2700, res = 250, bg="transparent")
DimPlot(object = Seuratobj_cd4, group.by = "Phase", pt.size = 6) & NoLegend() &  theme( axis.text.x=element_text(hjust=1, size=55),
                                                                                                     axis.title = element_text(size=55),
                                                                                                     axis.ticks = element_line(size = 3),
                                                                                                     axis.ticks.length = unit(7, "pt"),
                                                                                                     axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                                     axis.line = element_line(size=1.5),
                                                                                                     legend.background = element_rect(fill = "transparent"),
                                                                                                     panel.background = element_rect(fill = "transparent"),
                                                                                                     plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
dev.off()

Seuratobj_cd4 <- FindNeighbors(Seuratobj_cd4, dims = 1:6)

Seuratobj_cd4 <- FindClusters(Seuratobj_cd4, resolution = 0.4)

p2 <- DimPlot(Seuratobj_cd4, label = TRUE) + NoLegend()
ggarrange(p1,p2)

#Percentage####
Seuratobj_cd4$cell_cluster <- paste(Seuratobj_cd4$celltypes, sep = '_')
Seuratobj.cd4_table1 <- prop.table(table(Seuratobj_cd4$Phase, Seuratobj_cd4$orig.ident), margin = 2)
Seuratobj.cd4_table2 <- Seuratobj.cd4_table1*100
write.csv(Seuratobj.cd4_table2, "Percentages_clusters_cd4_cellcycle.csv")

# CD8+ ####
#finding the cell cycle effect###
Seuratobj_tcells <- readRDS("tcells.rds")
Seuratobj_cd8 <- subset(tcells, ident="0")

load("cycle.rda")

rownames(Seuratobj_cd8@assays$RNA@data) <- toupper(rownames(Seuratobj_cd8@assays$RNA@data))

Seuratobj_cd8<- CellCycleScoring(Seuratobj_cd8, 
                                 g2m.features=g2m_genes, 
                                 s.features=s_genes)

head(Seuratobj_cd8@meta.data)

Seuratobj_cd8 <- FindVariableFeatures(Seuratobj_cd8)

Seuratobj_cd8 <- ScaleData(Seuratobj_cd8)

Seuratobj_cd8 <- RunPCA(Seuratobj_cd8, 
                        features = VariableFeatures(Seuratobj_cd8), 
                        nfeatures.print = 10)

DimPlot(object = Seuratobj_cd8, reduction = "pca", group.by = "Phase", split.by = "orig.ident")


#Imaging Cell cycle effect####
Seuratobj_cd8 <- RunUMAP(Seuratobj_cd8, dims = 1:5)

png("UMAP_cd8cells_cellcylce.png", width = 2700, height = 2700, res = 250, bg="transparent")
DimPlot(object = Seuratobj_cd8, group.by = "Phase", pt.size = 6) & NoLegend() &  theme( axis.text.x=element_text(hjust=1, size=55),
                                                                                                      axis.title = element_text(size=55),
                                                                                                      axis.ticks = element_line(size = 3),
                                                                                                      axis.ticks.length = unit(7, "pt"),
                                                                                                      axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                        axis.line = element_line(size=1.5),
                                                                                                      legend.background = element_rect(fill = "transparent"),
                                                                                                      panel.background = element_rect(fill = "transparent"),
                                                                                                      plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
dev.off()

Seuratobj_cd8 <- FindNeighbors(Seuratobj_cd8, dims = 1:6)

Seuratobj_cd8 <- FindClusters(Seuratobj_cd8, resolution = 0.4)

p2 <- DimPlot(Seuratobj_cd8, label = TRUE) + NoLegend()
ggarrange(p1,p2)

#Percentage####
Seuratobj_cd8$cell_cluster <- paste(Seuratobj_cd8$celltypes, sep = '_')
Seuratobj.cd8_table1 <- prop.table(table(Seuratobj_cd8$Phase, Seuratobj_cd8$orig.ident), margin = 2)
Seuratobj.cd8_table2 <- Seuratobj.cd8_table1*100
write.csv(Seuratobj.cd8_table2, "Percentages_clusters_cd8_cellcycle.csv")

# cDC ####
#finding the cell cycle effect###
Seuratobj_cdc <- readRDS("cdc.rds")
Seuratobj_cdc <- subset(cdc, ident="0")

load("cycle.rda")

rownames(Seuratobj_cdc@assays$RNA@data) <- toupper(rownames(Seuratobj_cdc@assays$RNA@data))

Seuratobj_cdc<- CellCycleScoring(Seuratobj_cdc, 
                                 g2m.features=g2m_genes, 
                                 s.features=s_genes)

head(Seuratobj_cdc@meta.data)

Seuratobj_cdc <- FindVariableFeatures(Seuratobj_cdc)

Seuratobj_cdc <- ScaleData(Seuratobj_cdc)

Seuratobj_cdc <- RunPCA(Seuratobj_cdc, 
                        features = VariableFeatures(Seuratobj_cdc), 
                        nfeatures.print = 10)

DimPlot(object = Seuratobj_cdc, reduction = "pca", group.by = "Phase", split.by = "orig.ident")


#Imaging Cell cycle effect####
Seuratobj_cdc <- RunTSNE(Seuratobj_cdc, dims = 1:5)

png("UMAP_cdc_cellcylce.png", width = 2700, height = 2700, res = 250, bg="transparent")
DimPlot(object = Seuratobj_cdc, group.by = "Phase", pt.size = 6) & NoLegend() &  theme( axis.text.x=element_text(hjust=1, size=55),
                                                                                        axis.title = element_text(size=55),
                                                                                        axis.ticks = element_line(size = 3),
                                                                                        axis.ticks.length = unit(7, "pt"),
                                                                                        axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                        axis.line = element_line(size=1.5),
                                                                                        legend.background = element_rect(fill = "transparent"),
                                                                                        panel.background = element_rect(fill = "transparent"),
                                                                                        plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
dev.off()

Seuratobj_cdc <- FindNeighbors(Seuratobj_cdc, dims = 1:6)

Seuratobj_cdc <- FindClusters(Seuratobj_cdc, resolution = 0.4)

p2 <- DimPlot(Seuratobj_cdc, label = TRUE) + NoLegend()
ggarrange(p1,p2)

#Percentage####
Seuratobj_cdc$cell_cluster <- paste(Seuratobj_cdc$celltypes, sep = '_')
Seuratobj.cdc_table1 <- prop.table(table(Seuratobj_cdc$Phase, Seuratobj_cdc$orig.ident), margin = 2)
Seuratobj.cdc_table2 <- Seuratobj.cdc_table1*100
write.csv(Seuratobj.cdc_table2, "Percentages_clusters_cdc_cellcycle.csv")


# PP ####
#finding the cell cycle effect###
Seuratobj_PP <- readRDS("PP.rds")

load("cycle.rda")

rownames(Seuratobj_PP@assays$RNA@data) <- toupper(rownames(Seuratobj_PP@assays$RNA@data))

Seuratobj_PP<- CellCycleScoring(Seuratobj_PP, 
                                 g2m.features=g2m_genes, 
                                 s.features=s_genes)

head(Seuratobj_PP@meta.data)

Seuratobj_PP <- FindVariableFeatures(Seuratobj_PP)

Seuratobj_PP <- ScaleData(Seuratobj_PP)

Seuratobj_PP <- RunPCA(Seuratobj_PP, 
                        features = VariableFeatures(Seuratobj_PP), 
                        nfeatures.print = 10)

DimPlot(object = Seuratobj_PP, reduction = "pca", group.by = "Phase", split.by = "orig.ident")


#Imaging Cell cycle effect####
Seuratobj_PP <- RunTSNE(Seuratobj_PP, dims = 1:30)

png("UMAP_PP_cellcylce.png", width = 2700, height = 2700, res = 250, bg="transparent")
DimPlot(object = Seuratobj_PP, group.by = "Phase", pt.size = 6) & NoLegend() &  theme( axis.text.x=element_text(hjust=1, size=55),
                                                                                        axis.title = element_text(size=55),
                                                                                        axis.ticks = element_line(size = 3),
                                                                                        axis.ticks.length = unit(7, "pt"),
                                                                                        axis.text.y = element_text(hjust=0.5,size = 55),
                                                                                        axis.line = element_line(size=1.5),
                                                                                        legend.background = element_rect(fill = "transparent"),
                                                                                        panel.background = element_rect(fill = "transparent"),
                                                                                        plot.background = element_rect(fill = "transparent")) &
  guides(color = guide_legend(override.aes = list(size = 10) ) )
dev.off()

Seuratobj_PP <- FindNeighbors(Seuratobj_PP, dims = 1:30)

Seuratobj_PP <- FindClusters(Seuratobj_PP, resolution = 0.4)

p2 <- DimPlot(Seuratobj_PP, label = TRUE) + NoLegend()
ggarrange(p1,p2)

#Percentage####
Seuratobj_PP$cell_cluster <- paste(Seuratobj_PP$celltypes, sep = '_')
Seuratobj.PP_table1 <- prop.table(table(Seuratobj_PP$Phase, Seuratobj_PP$orig.ident), margin = 2)
Seuratobj.PP_table2 <- Seuratobj.PP_table1*100
write.csv(Seuratobj.PP_table2, "Percentages_clusters_PP_cellcycle.csv")


