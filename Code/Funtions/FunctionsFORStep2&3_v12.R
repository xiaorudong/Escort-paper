


### Similarity between cells
# obj: a trajectory evaluation object from TEvalobj() function
# Cluters: object from DCClusterscheck() fucntion.
# norm_counts: a normalized count data matrix: row:genes, column:cells
# dimred: A data frame. a 2D embedding from a DR method
# dist_mat: distance matrix of t(norm_counts) by using canberra distance.
# K: the number of clusters. This is could be got from step 1.
# output is a list containing:
#           + SimilarityLevel: a table contain if each cell's neighbors are similar as it. The higher frequency of the smaller value (e.g. 0, 1), the better results
#           + Fig: Results visualization
#           + all: nerighbors' information
#           + cls: cells' clusters

Similaritycheck <- function(norm_counts, obj, Cluters, dimred, PLOTname) {
  if (!missing(obj)) {
    dimred <- obj$Embedding
    PLOTname <- deparse(substitute(obj))
  }
  cut_avg <- Cluters$Clusters
  K <- Cluters$K

  if (!require("FNN")) install.packages("FNN")
  library(FNN)
  if (!require("scales")) install.packages("scales")
  library(scales)

  knn_index <- get.knn(dimred, k=3)$nn.index
  knn_group <- matrix(apply(knn_index, 1, function(x) cut_avg[x]), ncol = 3, byrow = T)
  knn_df <- cbind(cut_avg, knn_group)
  knn_overlap <- apply(knn_df, 1, function(x) sum(x[2:4]!=x[1]))
  # overlap_area <- sum(knn_overlap>1)
  t_knn <- table(knn_overlap) # not great results

  good_rate <- sum(t_knn[as.numeric(names(t_knn))<=1])/ncol(norm_counts)

  plotcol <- as.numeric(as.factor(cut_avg))

  plot(dimred, col = alpha(plotcol,0.7), pch=16, main=PLOTname)

  return(list(SimilarityLevel=t_knn, GoodRate=good_rate, all=knn_overlap, cls=cut_avg))
}






### Check the shape of trajectory curves
# obj: a trajectory evaluation object from TEvalobj() function
# norm_counts: a normalized count data matrix: row:genes, column:cells
# dimred: A data frame. a 2D embedding from a DR method
# K: the number of clusters. This is could be got from step 1.
# output is a list containing:
#           + SimilarityLevel: a table contain if each cell's neighbors are similar as it. The higher frequency of the smaller value (e.g. 0, 1), the better results
#           + Fig: Results visualization
#           + all: nerighbors' information
#           + cls: cells' clusters

UshapeDetector <- function(obj, dimred, PLOTname, outlierdetect='adjbox') {
  if (!require("FNN")) install.packages("FNN")
  library(FNN)
  if (!require("scales")) install.packages("scales")
  library(scales)
  if (!require("univOutl")) install.packages("univOutl")
  library(univOutl)
  if (!require("TSCAN")) BiocManager::install("TSCAN")
  library(TSCAN)
  if (!require("plyr")) install.packages("plyr")
  library(plyr)
  if (!require("RColorBrewer")) install.packages("RColorBrewer")
  library(RColorBrewer)

  
  if (!missing(obj)) {
    dimred <- obj$Embedding
    pse <- obj$pse
    fitLine <- obj$fitLine
    PLOTname <- deparse(substitute(obj))
  }
  pse <- as.data.frame(pse)
  
  nlineage <- ncol(pse)
  HR_cells_vec <- c()
  
  if((nrow(dimred)/nlineage)/100<10) {
    num <- 10
  } else {
    num <- round_any((nrow(dimred)/nlineage)/100, 10, f = round)
  }
  for (i in 1:nlineage) {
    lineage_pse <- pse[!is.na(pse[,i]),i]
    names(lineage_pse) <- rownames(pse)[!is.na(pse[,i])]
    lineage_dimred <- dimred[rownames(dimred) %in% names(lineage_pse),]
    lineage_pse <- lineage_pse[match(rownames(lineage_dimred), names(lineage_pse))]
    
    knn_obj <- get.knn(lineage_dimred, k=num)
    knn_index <- knn_obj$nn.index
    knn_dist <- knn_obj$nn.dist
    
    outdect <- boxB(x = as.vector(knn_dist), k = 1.5, method = 'adjbox')
    upperB <- outdect$fences[2]
    filter_ind <- which(knn_dist > upperB, arr.ind = TRUE)
    if(sum(knn_dist > upperB)>0) {
      for (i in 1:nrow(filter_ind)) {
        knn_dist[filter_ind[i, 1], filter_ind[i, 2]] <- NA
        knn_index[filter_ind[i, 1], filter_ind[i, 2]] <- NA
      }
    }
    
    knn_pt <- matrix(apply(knn_index, 1, function(x) lineage_pse[x]), ncol = num, byrow = T)
    knn_sd <- apply(knn_pt, 1, function(x) sd(x, na.rm = T))
    outdect <- boxB(x = knn_sd, k = 1.5, method = outlierdetect)
    upperB <- outdect$fences[2]
    if (any(knn_sd>upperB, na.rm = T)) HR <- which(knn_sd>upperB) else HR <- NULL
    # rm_cells <- as.numeric(names(table(filter_ind[,1]))[(num-table(filter_ind[,1]))<num/2])
    rm_cells <- as.numeric(names(table(filter_ind[,1])))

    HR <- HR[! HR %in% rm_cells]
    HR_cells <- rownames(lineage_dimred)[HR]
    HR_cells_vec <- c(HR_cells_vec, HR_cells)
  }
  allHR_cells <- unique(HR_cells_vec)
  allHR.no <- length(unique(HR_cells_vec))
  
  
  
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  pse$Ave <- rowMeans(pse, na.rm = T)
  plotcol <- colors[cut(pse$Ave, breaks=100)]
  plot(dimred, col = alpha(plotcol, 0.7), pch=16)
  segments(x0 = fitLine$x0,                   
           y0 = fitLine$y0,
           x1 = fitLine$x1,
           y1 = fitLine$y1, lwd = 3)
  points(dimred[allHR_cells,1], dimred[allHR_cells,2], col = "black", pch=16)

  return(list(AmbiguousCells=allHR_cells, NoACells=allHR.no))
}



### Goodness of Fit
# obj: a trajectory evaluation object from TEvalobj() function
# dimred: A data frame. a 2D embedding from a DR method
# output is a list containing:
#           + cellArea: the cell occupied area
#           + cellRate: the ratio of cell occupied area and minimum circle area
#           + NoDC: number of disconnected cells

GOFeval <- function(obj, dimred, PLOTname, alpha=NULL) {

  if (!require("alphahull")) install.packages("alphahull")
  library(alphahull)
  if (!require("elbow")) devtools::install_github("ahasverus/elbow", build_vignettes = TRUE)
  library(elbow)
  if (!require("shotGroups")) install.packages("shotGroups")
  library(shotGroups)
  if (!require("rstatix")) install.packages("rstatix")
  library(rstatix)

  if (!missing(obj)) {
    dimred <- obj$Embedding
    PLOTname <- deparse(substitute(obj))
  }

  # remove the outliers:
  outlier_md <- mahalanobis_distance(dimred)
  outliers <- rownames(outlier_md)[outlier_md$is.outlier]

  clean_dimred <- dimred
  if(length(outliers)>0) clean_dimred <- dimred[!rownames(dimred) %in% outliers,]

  minCircle <- getMinCircle(as.matrix(clean_dimred))
  circle_rad <- minCircle$rad
  circlearea <- pi*minCircle$rad^2

  if(is.null(alpha)) {

    upper_bound <- round(circle_rad+0.01, 2)

    if(circle_rad<1) {
      seq2 <- seq(0, upper_bound, by=0.05)
    } else {
      seq2 <- seq(0, upper_bound, length.out=15)
    }
    area <- c()
    for (a in seq2) {
      outline_dimred <- ahull(clean_dimred, alpha=a)
      area <- c(area, areaahull(outline_dimred))
    }
    # norm_area <- (area-min(area))/(max(area)-min(area))

    data <- data.frame(x=seq2, y=area)
    data <- data[order(data$x, decreasing = T),]
    # plot(data[,1:2])
    data$slope = c(NA, diff(data$y)/diff(data$x))
    data$slope_chg = c(NA, round(diff(data$slope),5))
    data$change = ifelse(data$slope_chg != 0, "change","")

    big_slope <- min(head(order(data$slope_chg), 2))
    alpha=data$x[big_slope]

    outline_dimred <- ahull(clean_dimred, alpha=alpha)
    area <- areaahull(outline_dimred)
    plot(outline_dimred, main=PLOTname)
    points(dimred[rownames(dimred) %in% outliers,1], dimred[rownames(dimred) %in% outliers,2])
  } else {
    outline_dimred <- ahull(clean_dimred, alpha=alpha)
    area <- areaahull(outline_dimred)
    plot(outline_dimred, main=PLOTname)
    points(dimred[rownames(dimred) %in% outliers,1], dimred[rownames(dimred) %in% outliers,2])
  }

  op_rate <- round(area/circlearea, 3)

  return(list(cellArea=area, alpha=alpha,
              alphahull_obj = outline_dimred,
              occupiedRate=op_rate))
}




source("CODE/MethodScripts/Funtions/FunctionsFORStep1_v18(dc).R")
LD_DCClusterscheck <- function(dist_mat, DRdims, cutoff=0.1,
                               max.nc=5, K=NULL, checkcells=NULL,
                               connectedCells=NULL, checksize=NULL) {

  if (!require("NbClust")) install.packages("NbClust")

  if(is.null(K)) {
    library(NbClust)
    res.nbclust <- NbClust(DRdims, distance = "euclidean", min.nc = 2, max.nc = max.nc, method = "complete", index ="all")
    K <- length(unique(res.nbclust$Best.partition))
    c_cl <- res.nbclust$Best.partition
    # names(c_cl) <- rownames(DRdims)
  } else {
    res.hc <- hclust(dist_mat, method = "complete")
    c_cl <- cutree(res.hc, k = K)
    # names(c_cl) <- rownames(DRdims)
  }
  if (any(table(c_cl)<=5)) {
    return(list(DCcheck="There are small clusters deteced. Please do clustering again",
                K=K, Clusters=c_cl))
  } else {
    BWClusters_Determination(dist_mat=dist_mat, K=K, c_cl=c_cl, cutoff=cutoff,
                             checkcells=checkcells, connectedCells=connectedCells,
                             checksize=checksize)
  }
}



score_cal <- function(scoredf, SimiRetaincutoff=0.8) {
  if (all(is.na(scoredf))) {
    final_scoredf <- scoredf
    final_scoredf$score <- rep(NA, nrow(final_scoredf))
    final_scoredf$decision <- rep("very bad(1)", nrow(final_scoredf))
    return(final_scoredf)
  }
  if(sum(scoredf$DCcheck & scoredf$SimiRetain>=SimiRetaincutoff)>0) {
    df <- scoredf[scoredf$DCcheck & scoredf$SimiRetain>=SimiRetaincutoff,]
    if (nrow(df)==1) {
      df$score <- 999
      df$ranking <- 999
      df$decision <- "Only choice"
      df <- df[,5:7]
      final_scoredf <- merge(scoredf, df, by="row.names", all=T)
    } else {
      df_norm <- as.data.frame(apply(df[,3:4], 2, function(x) -(x-max(x))/(max(x)-min(x))))
      colnames(df_norm) <- paste0("Norm_", colnames(df_norm))
      df_norm$score <- rowSums(df_norm)*df$SimiRetain
      final_scoredf <- merge(scoredf, df_norm, by="row.names", all=T)
      final_scoredf$ranking[!is.na(final_scoredf$score)] <- rank(-final_scoredf$score, na.last = NA)
      final_scoredf$decision <- ifelse(final_scoredf$score>1, "good", "bad")
      final_scoredf$decision[is.na(final_scoredf$decision)] <- "very bad"
    }
  }else {
    final_scoredf <- scoredf
    final_scoredf$score <- rep(NA, nrow(final_scoredf))
    final_scoredf$decision <- rep("very bad(1)", nrow(final_scoredf))
  }
  return(final_scoredf)
}




