
library(Escort)

setwd("path/to/your_file")
load("RDATA/Manuscript/clean_bone_SSPC.RData")
load("RDATA/Manuscript/dr_bone_SSPC.RData")
load("RDATA/Manuscript/orig_bone_SSPC_cls_afterQC.RData")


############################
######### Monocle3 #########
############################

library(monocle3)
library(igraph)
library(scales)

gene.var <- quick_model_gene_var(norm_counts)
genes.HVGs_1 <- getTopHVGs(gene.var, prop=0.1, var.threshold=NULL)
genes.HVGs_2 <- getTopHVGs(gene.var, prop=0.2, var.threshold=NULL)
genes.HVGs_3 <- getTopHVGs(gene.var, prop=0.5, var.threshold=NULL)
genes.HVGs_ls <- list(genes.HVGs_1, genes.HVGs_2, genes.HVGs_3)


# for subset 1
cds <- new_cell_data_set(rawcounts[genes.HVGs_1,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[7]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse7 <- pseudotime(cds, reduction_method = "UMAP")
fitLine7 <- as.data.frame(multi_segments)

# for subset 2
cds <- new_cell_data_set(rawcounts[genes.HVGs_2,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[8]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse8 <- pseudotime(cds, reduction_method = "UMAP")
fitLine8 <- as.data.frame(multi_segments)


# for subset 3
cds <- new_cell_data_set(rawcounts[genes.HVGs_3,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[9]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse9 <- pseudotime(cds, reduction_method = "UMAP")
fitLine9 <- as.data.frame(multi_segments)



####
#mds:
# for subset 1
cds <- new_cell_data_set(rawcounts[genes.HVGs_1,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[1]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse1 <- pseudotime(cds, reduction_method = "UMAP")
fitLine1 <- as.data.frame(multi_segments)


# for subset 2
cds <- new_cell_data_set(rawcounts[genes.HVGs_2,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[2]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse2 <- pseudotime(cds, reduction_method = "UMAP")
fitLine2 <- as.data.frame(multi_segments)



# for subset 3
cds <- new_cell_data_set(rawcounts[genes.HVGs_3,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[3]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse3 <- pseudotime(cds, reduction_method = "UMAP")
fitLine3 <- as.data.frame(multi_segments)



####
#tsne:
# for subset 1
cds <- new_cell_data_set(rawcounts[genes.HVGs_1,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[4]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse4 <- pseudotime(cds, reduction_method = "UMAP")
fitLine4 <- as.data.frame(multi_segments)



# for subset 2
cds <- new_cell_data_set(rawcounts[genes.HVGs_2,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[5]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse5 <- pseudotime(cds, reduction_method = "UMAP")
fitLine5 <- as.data.frame(multi_segments)



# for subset 3
cds <- new_cell_data_set(rawcounts[genes.HVGs_3,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[6]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse6 <- pseudotime(cds, reduction_method = "UMAP")
fitLine6 <- as.data.frame(multi_segments)



##
# pca
# for subset 1
cds <- new_cell_data_set(rawcounts[genes.HVGs_1,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[10]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse10 <- pseudotime(cds, reduction_method = "UMAP")
fitLine10 <- as.data.frame(multi_segments)


# for subset 2
cds <- new_cell_data_set(rawcounts[genes.HVGs_2,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[11]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse11 <- pseudotime(cds, reduction_method = "UMAP")
fitLine11 <- as.data.frame(multi_segments)



# for subset 3
cds <- new_cell_data_set(rawcounts[genes.HVGs_3,])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
reducedDims(cds)@listData[["UMAP"]] <- embeddings[[12]]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, use_partition=F)
cds <- order_cells(cds)
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

rawpse12 <- pseudotime(cds, reduction_method = "UMAP")
fitLine12 <- as.data.frame(multi_segments)

fitLines <- list(fitLine1, fitLine2, fitLine3,
                 fitLine4, fitLine5, fitLine6,
                 fitLine7, fitLine8, fitLine9,
                 fitLine10, fitLine11, fitLine12)
rawpses <- list(rawpse1, rawpse2, rawpse3,
                rawpse4, rawpse5, rawpse6,
                rawpse7, rawpse8, rawpse9,
                rawpse10, rawpse11, rawpse12)

load("RDATA/Manuscript/orig_bone_SSPC_m3_afterQC.RData")
fitLines[[13]] <- fitLineorig
rawpses[[13]] <- rawpseorig


library(parallel)
resobjs <- mclapply(1:9, function(x) {
  prepTraj(embeddings[[x]], PT=rawpses[[x]], fitLine=fitLines[[x]])
})


save(fitLines, resobjs, file = "RDATA/Manuscript/monocle3_bone_SSPC.RData")



