
library(Escort)

setwd("path/to/your_file")
load("RDATA/Manuscript/clean_DEG20.RData")
load("RDATA/Manuscript/dr_DEG20.RData")


############################
######### Monocle3 #########
############################

library(monocle3)
library(igraph)
library(scales)

startcell <- names(which.min(scdata_time))
gene.var <- quick_model_gene_var(norm_counts)
genes.HVGs_1 <- rownames(gene.var)[1:round(nrow(gene.var)*0.2)]
genes.HVGs_2 <- rownames(gene.var)[1:round(nrow(gene.var)*0.4)]
genes.HVGs_3 <- rownames(gene.var)
genes.HVGs_ls <- list(genes.HVGs_1, genes.HVGs_2, genes.HVGs_3)


fitLines <- list()
rawpses <- list()
for(i in 1:9) {
  gene_number <- ifelse(i%%3==0, 3, i%%3)
  genes <- genes.HVGs_ls[[gene_number]]
  cds <- new_cell_data_set(rawcounts)
  cds <- preprocess_cds(cds, num_dim = 100, use_genes = genes)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  
  reducedDims(cds)@listData[["UMAP"]] <- as.matrix(embeddings[[i]])
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
  cds <- learn_graph(cds, use_partition=F)
  cds <- order_cells(cds, root_cells = startcell)
  
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
  
  rawpse <- pseudotime(cds, reduction_method = "UMAP")
  rawpses[[i]] <- rawpse
  fitLines[[i]] <- as.data.frame(multi_segments)
}

library(parallel)
resobjs <- mclapply(1:9, function(x) {
  prepTraj(embeddings[[x]], PT=rawpses[[x]], fitLine=fitLines[[x]])
})


save(fitLines, resobjs, file = "RDATA/Manuscript/monocle3_DEG20.RData")



