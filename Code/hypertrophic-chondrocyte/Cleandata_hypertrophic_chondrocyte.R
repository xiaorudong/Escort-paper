# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8893718/
# figure 3f

rm(list = ls())

library(SeuratWrappers)
library(monocle3)
library(Seurat)
library(ggplot2)

setwd("path/to/your_file")

mydata <- readRDS("DATA/bone/GSM5726705_Final.rds")
# startcs <- colnames(mydata[, mydata@meta.data[["seurat_clusters"]]=="7"])
# save(startcs, file="RDATA/example_realData_bone/orig_startcs.RData")
DimPlot(mydata, reduction = "umap", label = T)

cds <- as.cell_data_set(mydata)
fData(cds)$gene_short_name <- rownames(fData(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- mydata@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- mydata@reductions$umap@cell.embeddings

cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")

cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds, reduction_method = "UMAP")
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

save.image("RDATA/res_bone_SSPC.RData")



VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

subdata <- subset(mydata, seurat_clusters!=2)
subdata <- RenameIdents(object = subdata, "0" = "Cluster 3", "1" = "Cluster 5",
                        "3" = "Cluster 6", "4" = "Cluster 8",
                        "5" = "Cluster 2", "6" = "Cluster 4", "7" = "Cluster 1")

old_cls <- droplevels(subdata@active.ident)
DimPlot(subdata, reduction = "umap", label = T)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
subdata <- NormalizeData(subdata)
subdata <- CellCycleScoring(subdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
subdata = SCTransform(subdata, vars.to.regress= c("S.Score", "G2M.Score"))

subdata <- RunPCA(subdata, features = head(VariableFeatures(subdata), 2000))
subdata <- FindNeighbors(subdata, dims = 1:35)
subdata <- FindClusters(subdata, resolution = 0.5)
subdata <- RunUMAP(subdata, dims = 1:35)

subdata@meta.data[["old_cls"]] <- old_cls
DimPlot(subdata, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(subdata, reduction = "umap", label = TRUE, repel = TRUE, group.by = "old_cls")
table(subdata@active.ident, subdata@meta.data[["old_cls"]])

subdata@active.ident <- subdata@meta.data[["old_cls"]]
list.cluster <- subdata@active.ident
names(list.cluster) <- rownames(subdata@meta.data)
startcs <- colnames(subdata[, subdata@active.ident=="Cluster 1"])

save(subdata, file="RDATA/Manuscript/res_bone_SSPC_afterQC.RData")
save(list.cluster, file="RDATA/Manuscript/orig_bone_SSPC_cls_afterQC.RData")
save(startcs, file="RDATA/Manuscript/orig_bone_SSPC_startcs_afterQC.RData")

cds <- as.cell_data_set(subdata)
fData(cds)$gene_short_name <- rownames(fData(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- subdata@reductions$umap@cell.embeddings

cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                 group_label_size = 5) + theme(legend.position = "right")

cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds, reduction_method = "UMAP")
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)




# extract trajectory info
library(igraph)
library(RColorBrewer)
umap_dim <- reducedDims(cds)@listData[["UMAP"]]
curve_igraph <- as_long_data_frame(cds@principal_graph@listData[["UMAP"]])
mst <- t(cds@principal_graph_aux@listData[["UMAP"]][["dp_mst"]])
curve_igraph_sub <- curve_igraph[,4:5]
multi_segments <- as.data.frame(matrix(NA, nrow = nrow(curve_igraph_sub), ncol=4))
colnames(multi_segments) <- c("x0", "y0", "x1", "y1")

multi_segments$x0 <- mst[curve_igraph_sub[,1],1]
multi_segments$y0 <- mst[curve_igraph_sub[,1],2]
multi_segments$x1 <- mst[curve_igraph_sub[,2],1]
multi_segments$y1 <- mst[curve_igraph_sub[,2],2]

rawpseorig <- pseudotime(cds, reduction_method = "UMAP")
fitLineorig <- as.data.frame(multi_segments)

save(rawpseorigl, fitLineorig, file = "RDATA/Manuscript/orig_bone_SSPC_m3_afterQC.RData")



