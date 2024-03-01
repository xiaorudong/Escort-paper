rm(list = ls())

library(scran)
library(Seurat)
library(dplyr)
library(RColorBrewer)

# https://github.com/dynverse/dynbenchmark
# https://zenodo.org/record/1443566#.YYMvcC-B19c
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")

source("CODE/MethodScripts/Funtions/FunctionsFORStep0.R")


load("RDATA/step0_clean_realGold_pancreatic.RData")
load("RDATA/step0_dr_realGold_pancreatic.RData")

#######################################
######### Trajectory Analysis #########
#######################################

# select genes
library(scuttle)
library(scran)
gene.var <- modelGeneVar(norm_counts)
genes.HVGs_1 <- getTopHVGs(gene.var, prop=0.2, var.threshold=NULL)
genes.HVGs_2 <- getTopHVGs(gene.var, prop=0.4, var.threshold=NULL)
genes.HVGs_3 <- getTopHVGs(gene.var, prop=1, var.threshold=NULL)


# Take using Monocle3 to do trajectory analysis as an example
library(monocle3)
library(igraph)
library(scales)

# for subset 1
cds <- new_cell_data_set(rawcounts[genes.HVGs_1,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- dimred7
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds, root_cells = names(which.min(scdata_time)))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# extract trajectory info
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

pt <- pseudotime(cds, reduction_method = "UMAP")
pse <- pt/max(pt)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]
plot(umap_dim, pch=20, col = alpha(plotcol, 0.5))
fitLine <- segments(x0 = multi_segments$x0,                   
                    y0 = multi_segments$y0,
                    x1 = multi_segments$x1,
                    y1 = multi_segments$y1, lwd = 2)

ti_out7 <- cds
dimred7 <- as.data.frame(umap_dim)
rawpse7 <- pt
pse7 <- pse
fitLine7 <- as.data.frame(multi_segments)

# for subset 2
cds <- new_cell_data_set(rawcounts[genes.HVGs_2,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- dimred8
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds, root_cells = names(which.min(scdata_time)))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# extract trajectory info
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

pt <- pseudotime(cds, reduction_method = "UMAP")
pse <- pt/max(pt)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]
plot(umap_dim, pch=20, col = alpha(plotcol, 0.5))
fitLine <- segments(x0 = multi_segments$x0,                   
                    y0 = multi_segments$y0,
                    x1 = multi_segments$x1,
                    y1 = multi_segments$y1, lwd = 2)

ti_out8 <- cds
dimred8 <- as.data.frame(umap_dim)
rawpse8 <- pt
pse8 <- pse
fitLine8 <- as.data.frame(multi_segments)

# for subset 3
cds <- new_cell_data_set(rawcounts[genes.HVGs_3,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- dimred9
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds, root_cells = names(which.min(scdata_time)))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# extract trajectory info
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

pt <- pseudotime(cds, reduction_method = "UMAP")
pse <- pt/max(pt)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]
plot(umap_dim, pch=20, col = alpha(plotcol, 0.5))
fitLine <- segments(x0 = multi_segments$x0,                   
                    y0 = multi_segments$y0,
                    x1 = multi_segments$x1,
                    y1 = multi_segments$y1, lwd = 2)

ti_out9 <- cds
dimred9 <- as.data.frame(umap_dim)
rawpse9 <- pt
pse9 <- pse
fitLine9 <- as.data.frame(multi_segments)

####
#mds:
# for subset 1
cds <- new_cell_data_set(rawcounts[genes.HVGs_1,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- dimred1
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds, root_cells = names(which.min(scdata_time)))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# extract trajectory info
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

pt <- pseudotime(cds, reduction_method = "UMAP")
pse <- pt/max(pt)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]
plot(umap_dim, pch=20, col = alpha(plotcol, 0.5))
fitLine <- segments(x0 = multi_segments$x0,                   
                    y0 = multi_segments$y0,
                    x1 = multi_segments$x1,
                    y1 = multi_segments$y1, lwd = 2)

ti_out1 <- cds
dimred1 <- as.data.frame(umap_dim)
rawpse1 <- pt
pse1 <- pse
fitLine1 <- as.data.frame(multi_segments)

# for subset 2
cds <- new_cell_data_set(rawcounts[genes.HVGs_2,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- dimred2
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds, root_cells = names(which.min(scdata_time)))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# extract trajectory info
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

pt <- pseudotime(cds, reduction_method = "UMAP")
pse <- pt/max(pt)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]
plot(umap_dim, pch=20, col = alpha(plotcol, 0.5))
fitLine <- segments(x0 = multi_segments$x0,                   
                    y0 = multi_segments$y0,
                    x1 = multi_segments$x1,
                    y1 = multi_segments$y1, lwd = 2)

ti_out2 <- cds
dimred2 <- as.data.frame(umap_dim)
rawpse2 <- pt
pse2 <- pse
fitLine2 <- as.data.frame(multi_segments)

# for subset 3
cds <- new_cell_data_set(rawcounts[genes.HVGs_3,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- dimred3
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds, root_cells = names(which.min(scdata_time)))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# extract trajectory info
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

pt <- pseudotime(cds, reduction_method = "UMAP")
pse <- pt/max(pt)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]
plot(umap_dim, pch=20, col = alpha(plotcol, 0.5))
fitLine <- segments(x0 = multi_segments$x0,                   
                    y0 = multi_segments$y0,
                    x1 = multi_segments$x1,
                    y1 = multi_segments$y1, lwd = 2)

ti_out3 <- cds
dimred3 <- as.data.frame(umap_dim)
rawpse3 <- pt
pse3 <- pse
fitLine3 <- as.data.frame(multi_segments)


####
#tsne:
# for subset 1
cds <- new_cell_data_set(rawcounts[genes.HVGs_1,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- dimred4
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds, root_cells = names(which.min(scdata_time)))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# extract trajectory info
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

pt <- pseudotime(cds, reduction_method = "UMAP")
pse <- pt/max(pt)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]
plot(umap_dim, pch=20, col = alpha(plotcol, 0.5))
fitLine <- segments(x0 = multi_segments$x0,                   
                    y0 = multi_segments$y0,
                    x1 = multi_segments$x1,
                    y1 = multi_segments$y1, lwd = 2)

ti_out4 <- cds
dimred4 <- as.data.frame(umap_dim)
rawpse4 <- pt
pse4 <- pse
fitLine4 <- as.data.frame(multi_segments)

# for subset 2
cds <- new_cell_data_set(rawcounts[genes.HVGs_2,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- dimred5
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds, root_cells = names(which.min(scdata_time)))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# extract trajectory info
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

pt <- pseudotime(cds, reduction_method = "UMAP")
pse <- pt/max(pt)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]
plot(umap_dim, pch=20, col = alpha(plotcol, 0.5))
fitLine <- segments(x0 = multi_segments$x0,                   
                    y0 = multi_segments$y0,
                    x1 = multi_segments$x1,
                    y1 = multi_segments$y1, lwd = 2)

ti_out5 <- cds
dimred5 <- as.data.frame(umap_dim)
rawpse5 <- pt
pse5 <- pse
fitLine5 <- as.data.frame(multi_segments)

# for subset 3
cds <- new_cell_data_set(rawcounts[genes.HVGs_3,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- dimred6
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds, root_cells = names(which.min(scdata_time)))
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
# extract trajectory info
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

pt <- pseudotime(cds, reduction_method = "UMAP")
pse <- pt/max(pt)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pt, breaks=100)]
plot(umap_dim, pch=20, col = alpha(plotcol, 0.5))
fitLine <- segments(x0 = multi_segments$x0,                   
                    y0 = multi_segments$y0,
                    x1 = multi_segments$x1,
                    y1 = multi_segments$y1, lwd = 2)

ti_out6 <- cds
dimred6 <- as.data.frame(umap_dim)
rawpse6 <- pt
pse6 <- pse
fitLine6 <- as.data.frame(multi_segments)




# generate the traj objs

M3_obj1 <- info_traj(dimred1, rawpse = rawpse1, fitLine = fitLine1)
M3_obj2 <- info_traj(dimred2, rawpse = rawpse2, fitLine = fitLine2)
M3_obj3 <- info_traj(dimred3, rawpse = rawpse3, fitLine = fitLine3)
M3_obj4 <- info_traj(dimred4, rawpse = rawpse4, fitLine = fitLine4)
M3_obj5 <- info_traj(dimred5, rawpse = rawpse5, fitLine = fitLine5)
M3_obj6 <- info_traj(dimred6, rawpse = rawpse6, fitLine = fitLine6)
M3_obj7 <- info_traj(dimred7, rawpse = rawpse7, fitLine = fitLine7)
M3_obj8 <- info_traj(dimred8, rawpse = rawpse8, fitLine = fitLine8)
M3_obj9 <- info_traj(dimred9, rawpse = rawpse9, fitLine = fitLine9)


MDS1_eval <- TEvalobj(dimred1, traj = M3_obj1)
MDS2_eval <- TEvalobj(dimred2, traj = M3_obj2)
MDS3_eval <- TEvalobj(dimred3, traj = M3_obj3)
TSNE4_eval <- TEvalobj(dimred4, traj = M3_obj4)
TSNE5_eval <- TEvalobj(dimred5, traj = M3_obj5)
TSNE6_eval <- TEvalobj(dimred6, traj = M3_obj6)
UMAP7_eval <- TEvalobj(dimred7, traj = M3_obj7)
UMAP8_eval <- TEvalobj(dimred8, traj = M3_obj8)
UMAP9_eval <- TEvalobj(dimred9, traj = M3_obj9)



save(MDS1_eval, MDS2_eval, MDS3_eval,
     TSNE4_eval, TSNE5_eval, TSNE6_eval,
     UMAP7_eval, UMAP8_eval, UMAP9_eval,
     file = "RDATA/step0_m3_realGold_pancreatic.RData")







