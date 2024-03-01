
# setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")


pickHVGs <- function(norm_counts, num) {
  if (!require("scran")) install.packages("scran")
  library(scran)
  
  gene.var <- modelGeneVar(norm_counts)
  genes.HVGs_1 <- getTopHVGs(gene.var, n=1000)
  genes.HVGs_2 <- getTopHVGs(gene.var, n=2500)
  genes.HVGs_3 <- rownames(gene.var[order(gene.var$bio, decreasing = T),])
  
  pickgenes <- lapply(num, function(x) getTopHVGs(gene.var, n=x))
  
  return(pickgenes)
}



### Generate 2D Embeddings by different methods:
# Here we provide 4 methods: MDS, UMAP, PCA, TSNE
# norm_counts: a normalized count data matrix: row:genes, column:cells
# method: choose one RD method: MDS, UMAP, PCA, TSNE
# output is coordinates of the embedding.
DR_2D <- function(norm_counts, method) {
  
  if (method=="MDS") {
    if (!require("SCORPIUS")) install.packages("SCORPIUS")
    dimred <- SCORPIUS::reduce_dimensionality(t(norm_counts), dist="spearman", ndim = 2)
  }
  
  if(method=="UMAP") {
    if (!require("umap")) install.packages("umap")
    dimred <- umap::umap(t(norm_counts))$layout
  }
  
  if(method=="PCA") {
    if (!require("stats")) install.packages("stats")
    dimred <- stats::prcomp(t(norm_counts), rank. = 2)$x
    rownames(dimred) <- rownames(t(norm_counts))
  }
  
  if(method=="TSNE") {
    if (!require("Rtsne")) install.packages("Rtsne")
    dimred <- Rtsne::Rtsne(t(norm_counts), dims = 2)$Y
    rownames(dimred) <- rownames(t(norm_counts))
  }
  
  if(method=="ICA") {
    if (!require("fastICA")) install.packages("fastICA")
    dimred <- fastICA::fastICA(t(norm_counts), 2)$S
  }
  dimred <- as.data.frame(dimred)
  return(dimred)
}


### Trajectory Evaluation Object Generation
# norm_counts: a normalized count data matrix: row:genes, column:cells. This data set didn't apply feature selection.
# dimred: A data frame. a 2D embedding from a DR method
# pse: A data frame. estimated pseudotime from TI. Each column contains the PT for one lineage.
# fitLine: the fitted curve line got from TI.
# output is a list containing:
#           + Normcounts: normalized count data matrix: row:genes, column:cells. This data set didn't apply feature selection.
#           + Embedding: A data frame. a 2D embedding from a DR method


TEvalobj <- function(dimred, traj=NULL) {
  if(is.null(traj)) {
    traj=info_traj(dimred)
  }
  rawpse <- as.data.frame(traj$PT)
  pse <- rawpse / max(rawpse, na.rm=TRUE)
  pse <- pse[match(rownames(dimred), rownames(pse)),]
  pse <- as.data.frame(pse)
  rownames(pse) <- rownames(dimred)
  
  obj <- list(Embedding=dimred, pse=pse, fitLine=traj$fitLine)
  return(obj)
}


info_traj <- function(dimred, rawpse=NULL, fitLine=NULL, K=NULL, max.nc=5) {
  
  if(is.null(fitLine) | is.null(rawpse)) {
    
    if (!require("plyr")) install.packages("plyr")
    library(plyr)
    if (!require("NbClust")) install.packages("NbClust")
    library(NbClust)
    if (!require("TSCAN")) BiocManager::install("TSCAN")
    library(TSCAN)
    
    if((nrow(dimred))/100<10) {
      num <- 10
    } else {
      num <- round_any((nrow(dimred))/100, 10, f = round)
    }
    # get clusters
    if(is.null(K)) {
      res.nbclust <- NbClust(dimred, distance = "euclidean", min.nc = 2, max.nc = max.nc, method = "complete", index ="all")
      K <- length(unique(res.nbclust$Best.partition))
      c_cl <- res.nbclust$Best.partition
      # names(c_cl) <- rownames(DRdims)
    } else {
      res.hc <- hclust(dist_mat, method = "complete")
      c_cl <- cutree(res.hc, k = K)
      # names(c_cl) <- rownames(DRdims)
    }
    
    # PT from MST
    out <- TSCAN::quickPseudotime(dimred, clusters=c_cl)
    rawpse <- rowMeans(pathStat(out$ordering), na.rm=TRUE)
    
    edges <- reportEdges(x=NULL, out$mst, combined=FALSE)
    df <- data.frame(x0=edges$start[,2], y0=edges$start[,3], 
                     x1=edges$end[,2], y1=edges$end[,3])
    
    obj <- list(PT=rawpse, fitLine=df)
  } else {
    obj <- list(PT=rawpse, fitLine=fitLine)
  }
  return(obj)
}






