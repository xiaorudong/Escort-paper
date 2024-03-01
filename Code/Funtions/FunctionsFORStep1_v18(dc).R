

### Test Homogeneous: 
# HVGs: highly variable genes
# norm_counts: a normalized count data matrix: row:genes, column:cells
# output is a correlation matrix
testHomogeneous <- function(HVGs, norm_counts, n=100, pct_cutoff=0.46, num.sim = 20000) {
  if(length(HVGs)<20) return("Please input more highly variable genes.")
  
  if (!require("factoextra")) install.packages("factoextra")
  if (!require("jmuOutlier")) install.packages("jmuOutlier")
  library(factoextra)
  pcDat <- prcomp(t(norm_counts)) 
  pca_cells <- pcDat$x
  est_pt <- sort(pca_cells[,1])
  
  library(jmuOutlier)
  time_df <- data.frame(Cell=names(est_pt), PT=est_pt)
  time_df$Cell <- factor(time_df$Cell, levels=colnames(norm_counts))
  time_df <- time_df[order(time_df$Cell),]
  comb_df <- cbind(time_df, t(norm_counts[head(HVGs, n),]))
  comb_df <- comb_df[order(comb_df$PT),]
  
  p_vec_perm <- apply(comb_df[,3:ncol(comb_df)], 2, 
                      function(x) perm.cor.test(comb_df[,2], x, method = "spearman", num.sim = num.sim)$p.value)
  padj <- p.adjust(p_vec_perm, method = "fdr")
  # res_perm_p <- mean(p_vec_perm<0.05)
  res_perm_padj <- mean(padj<0.05)
  decision <- ifelse(res_perm_padj>=pct_cutoff, 
                     "The trajectory signal is detected.",
                     "No trajectory signal is detected.")
  
  return(list(signal_pct = res_perm_padj, decision=decision))
}



### Differentiate Clusters vs Lineages
# dist_mat: a distance object that contains the distance between cells from dist() function
# norm_counts: a normalized count data matrix: row:genes, column:cells
# K: the number of clusters
# output is a list containing: 
#           + DCcheck: whether there exist null spaces between clusters.
#           + ElbowP: a elbow plot to choose the cluster number
#           + K: the number of clusters
#           + HCluster: clusters to which each cell is allocated by perform hierarchical clustering
#           + WithinDist: a distance vector containing cell-to-cell distances within clusters
#           + BetweenDist: a distance vector containing cell-to-cell distances between clusters
HD_DCClusterscheck <- function(dist_mat, rawcounts, 
                               K=NULL, clust.max=10, 
                               cutoff=0.3, checkcells=NULL, 
                               connectedCells=NULL, checksize=NULL, seed=1111) {
 
  if (!require("devtools")) install.packages("devtools")
  if (!require("cluster")) install.packages("cluster")
  if (!require("mclust")) install.packages("mclust")
  if (!require("RMTstat")) install.packages("RMTstat")
  
  library(devtools)
  
  set.seed(seed)
  if(is.null(K)) {
    # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02622-0
    if (!require("scLCA")) install_bitbucket("scLCA/single_cell_lca")
    library(scLCA)
    myclust.res <- myscLCA(rawcounts, clust.max=clust.max)
    c_cl <- myclust.res[[1]]
    names(c_cl) <- colnames(rawcounts)
    K <- max(c_cl)
  } else {
    # if (!require("SIMLR")) BiocManager::install("SIMLR")
    # library(SIMLR)
    # myclust.res <- SIMLR(X=rawcounts, c=K, normalize=T)
    # c_cl <- myclust.res$y$cluster
    # names(c_cl) <- colnames(rawcounts)
    # K <- max(c_cl)
    # https://link.springer.com/article/10.1186/s12859-021-04210-8#Sec6
    
    # library(devtools)
    # if (!require("SOUP")) devtools::install_github("lingxuez/SOUP")
    # library(SOUP)
    # log_exp <- log2(rawcounts+1)
    # res = SOUP(t(log_exp), Ks=K, type = "log")
    # c_cl_raw <- res$major.labels[[1]]
    # c_cl <- as.numeric(factor(c_cl_raw, levels = sort(unique(c_cl_raw))))
    # names(c_cl) <- names(c_cl_raw)
    # K <- max(c_cl)
    
    # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02622-0
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7444317/
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    if (!require("SingleCellExperiment")) BiocManager::install("SingleCellExperiment")
    # if (!require("SC3")) devtools::install_github("hemberg-lab/SC3")
    if (!require("sc3")) BiocManager::install("SC3")
    
    
    library(SingleCellExperiment)
    library(SC3)
    # create a SingleCellExperiment object
    sce <- SingleCellExperiment(
      assays = list(
        counts = as.matrix(rawcounts),
        logcounts = log2(rawcounts + 1)))
    
    # define feature names in feature_symbol column
    rowData(sce)$feature_symbol <- rownames(sce)
    res <- sc3(sce, ks = K, biology = TRUE)
    c_cl_raw <- colData(res)[,1]
    c_cl <- as.numeric(factor(c_cl_raw, levels = sort(unique(c_cl_raw))))
    names(c_cl) <- rownames(colData(res))
    K <- max(c_cl)
  }
  

  if(K==1) {
    DCdecision <- "We only detect one cluster. Please do it again"
    return(list(DCcheck=DCdecision, K=K, Clusters=c_cl, ifConnected=NA))
  } 
  if (K>1) {
    if (any(table(c_cl)<=10)) {
      return(list(DCcheck="There are small clusters deteced. Please do clustering again",
                  K=K, Clusters=c_cl, ifConnected=NA))
    }
    if (any(table(c_cl)>10)) {
      BWClusters_Determination(dist_mat=dist_mat, K=K, c_cl=c_cl, 
                               cutoff=cutoff, checkcells=checkcells,
                               connectedCells=connectedCells, checksize=checksize)
    }
  }
}



####

BWClusters_Determination <- function(dist_mat, K, c_cl, cutoff, checkcells, connectedCells, checksize) {
  if (!require("reshape2")) install.packages("reshape2")
  library(reshape2)
  # if (!require("rstatix")) install.packages("rstatix")
  # library(rstatix)
  
  res <- BWClusters(dist_mat, c_cl)
  JaccardIndex_ls <- res$JaccardIndex
  JaccardIndex_op_ls <- res$signJaccardIndex  
  Dists_W <- res$WithinDist
  Dists_B <- res$BetweenDist
  
  if(is.null(checkcells)) {
    # checkcells <- min(round(mean(table(c_cl)/5)), min(table(c_cl)))
    # checkcells <- round(mean(table(c_cl)/5))
    checkcells <- max(round(min(table(c_cl)/5)), 10)
  } 
  topcc <-lapply(JaccardIndex_ls, function(x) lapply(x, function(y) head(sort(y, decreasing = T),checkcells)))
  # topcc <-lapply(topcc, function(x) lapply(x, function(y) y[!is_outlier(y, coef = 1.5)]))
  
  topcc_tls <-lapply(1:length(topcc), function(x) lapply(1:length(topcc[[x]]), function(y) {
    subls <- topcc[[x]][[y]]
    df <- data.frame(value=names(subls))
    df$ToGroup <- rep(names(topcc[[x]])[y], nrow(df))
    df$InGroup <- rep(names(topcc)[x], nrow(df))
    return(df)
  }))
  
  topcc_tls <- lapply(topcc_tls, function(x) do.call(rbind, x))
  topcc_df <- as.data.frame(do.call(rbind, topcc_tls))
  
  combnt <- as.data.frame(t(combn(1:K, 2)))
  ggcomparison_ls <- list()
  connectedc_ls <- list()
  dist_matmat <- as.matrix(dist_mat)
  if(is.null(checksize)) {
    checksize <- checkcells
  } 
  for (i in 1:nrow(combnt)) {
    findg <- paste("Group", combnt[i,], sep = "_")
    needrow <- topcc_df$ToGroup %in% findg & topcc_df$InGroup %in% findg
    needc <- topcc_df[needrow,]
    ingroupc <- lapply(1:nrow(needc), function(x) {
      a <- names(head(sort(dist_matmat[needc$value[x], names(c_cl)[c_cl==as.numeric(strsplit(needc$InGroup[x], split = "_")[[1]][2])]]), checksize))
      b <- rep(as.numeric(strsplit(needc$InGroup[x], split = "_")[[1]][2]), length(a))
      names(b) <- a
      return(b)
    })
    
    betweengroupc <- lapply(1:nrow(needc), function(x) {
      a <- names(head(sort(dist_matmat[needc$value[x], names(c_cl)[c_cl==as.numeric(strsplit(as.character(needc$ToGroup[x]), split = "_")[[1]][2])]]), checksize))
      b <- rep(as.numeric(strsplit(as.character(needc$ToGroup[x]), split = "_")[[1]][2]), length(a))
      names(b) <- a
      return(b)
    })
    
    connectedc <- list()
    for (each in 1:nrow(needc)) {
      sub_ls <- c(ingroupc[[each]], betweengroupc[[each]])
      c_cl_sub <- sub_ls
      dist_mat_sub <- dist_matmat[names(c_cl_sub), names(c_cl_sub)]
      res_eachc <- BWClusters(dist_mat_sub, c_cl_sub)
      res_JaccardIndex <- sapply(res_eachc$JaccardIndex, cbind)
      sub <- sapply(res_JaccardIndex, function(x) names(x)[which(x>cutoff)])
      connectedc[[each]] <- as.vector(unlist(sub))
    }
    connectedc_ls[[i]] <- unique(unlist(connectedc))
    needc$ConnectedCellID <- unlist(lapply(connectedc, function(x) paste(x, collapse = ",")))
    ggcomparison_ls[[i]] <- needc
  }
  combnt$Connected <- sapply(connectedc_ls, function(x) {
    if(is.null(connectedCells)) {
      # connectedCells=10
      connectedCells=round(min(table(c_cl)/5))
      if(connectedCells<3) connectedCells==3
      if(connectedCells>10) connectedCells==10
    }
    des_check <- ifelse(length(x)>connectedCells, T, F)
    return(des_check)
    })
  
  check_mat <- combnt[combnt$Connected,c(1,2)]
  check1 <- length(unique(as.vector(unlist(check_mat))))==K
  check2 <- nrow(check_mat)>=K-1
  check3 <- apply(check_mat, 2, function(x) any(duplicated(x)))
  if(sum(check3)==0) {
    check <- check1 & check2
  } 
  if(sum(check3)!= 0) {
    dup_V1 <- check_mat$V1[duplicated(check_mat$V1)]
    dup_V2 <- check_mat$V2[duplicated(check_mat$V2)]
    cg_check <- unique(unlist(check_mat[check_mat$V1 %in% dup_V1 | check_mat$V2 %in% dup_V2, ]))
    dg_check <- unique(unlist(check_mat[!(check_mat$V1 %in% dup_V1 | check_mat$V2 %in% dup_V2), ]))
    
    check4 <- ifelse(sum(duplicated(c(cg_check, dg_check)))==0, length(cg_check), length(unique(c(cg_check, dg_check))))
    check4 <- check4==K
    
    check <- check1 & check2 & check4
  }
  
  DCdecision <- ifelse(check,"Congratulations! We didn't find null spaces between clusters. Go ahead doing trajectory anlaysis.", 
                       "There are null spaces between clusters. We don't recommend doing trajectory analysis on this dataset")
  
  # DC_vec <- as.vector(na.omit(DC_vec))
  # DCdecision <- ifelse(any(DC_vec), "There are null spaces between clusters. We don't recommend doing trajectory analysis on this dataset",
  #                      "Congratulations! We didn't find null spaces between clusters. Go ahead doing trajectory anlaysis.")
  return(list(DCcheck=DCdecision, ifConnected=check, Jaccardsummary=combnt, 
              clusterLocation=check_mat, K=K, Clusters=c_cl, 
              ConnectedCells=ggcomparison_ls, No_ConnectedCells=connectedc_ls))
  
}


BWClusters <- function(dist_mat, c_cl) {
  
  dist_mat2 <- as.matrix(dist_mat)
  
  Dists_W <- list()
  Dists_B <- list()
  JaccardIndex_ls <- list()
  JaccardIndex_op_ls <- list()
  
  # if (!require("univOutl")) install.packages("univOutl")
  # library(univOutl)
  for (i in unique(c_cl)) {
    sub_c <- dist_mat2[which(c_cl==i), which(c_cl==i)]
    # sub_sd <- sd(sub_c[lower.tri(sub_c)])
    Dists_W_cells <- lapply(1:nrow(sub_c), function(x) sub_c[x,-x])
    names(Dists_W_cells) <- rownames(sub_c)
    Dists_W[[paste("Group", i, sep="_")]] <- Dists_W_cells
    
    notsub_c_ls <- list()
    gp <- c(unique(c_cl))[! unique(c_cl) %in% i]
    
    for (gp1 in gp) {
      notsub_c <- dist_mat2[which(c_cl==i), which(c_cl==gp1)]
      Dists_B_cells <- lapply(1:nrow(notsub_c), function(x) notsub_c[x,])
      names(Dists_B_cells) <- rownames(notsub_c)
      notsub_c_ls[[paste("Group", gp1, sep = "_")]] <- Dists_B_cells
    }
    Dists_B[[paste("Group", i, sep="_")]] <- notsub_c_ls
    
    JaccardIndex_gp_ls <- list()
    JaccardIndex_gp_op_ls <- list()
    for (gp1 in names((Dists_B[[paste("Group", i, sep="_")]]))) {
      Dists_B_cells <- Dists_B[[paste("Group", i, sep="_")]][[gp1]]
      JaccardIndex_vec <- c()
      JaccardIndex_op_vec <- c()
      for (c in 1:length(Dists_B_cells)) {
        w_vec <- Dists_W_cells[[c]]
        b_vec <- Dists_B_cells[[c]]
        
        Jaccard <- JaccardIndex_fun(w_vec, b_vec, plot=F)
        JaccardIndex_vec <- c(JaccardIndex_vec, Jaccard$JaccardIndex)
        JaccardIndex_op_vec <- c(JaccardIndex_op_vec, Jaccard$signJaccardIndex)
      }
      names(JaccardIndex_vec) <- names(Dists_B_cells)
      names(JaccardIndex_op_vec) <- names(Dists_B_cells)
      JaccardIndex_vec <- ifelse(JaccardIndex_vec<0, 1, JaccardIndex_vec)
      JaccardIndex_gp_ls[[gp1]] <- JaccardIndex_vec
      JaccardIndex_gp_op_ls[[gp1]] <- JaccardIndex_op_vec
    }
    JaccardIndex_ls[[paste("Group", i, sep="_")]] <- JaccardIndex_gp_ls
    JaccardIndex_op_ls[[paste("Group", i, sep="_")]] <- JaccardIndex_gp_op_ls
  }
  return(list(JaccardIndex=JaccardIndex_ls, signJaccardIndex=JaccardIndex_op_ls,
              WithinDist=Dists_W, BetweenDist=Dists_B))
}






###
JaccardIndex_fun <- function(w_vec, b_vec, plot=F) {
  
  if (!require("sfsmisc")) install.packages("sfsmisc")
  library(sfsmisc)

  lower <- min(c(w_vec, b_vec)) - 0.2
  upper <- max(c(w_vec, b_vec)) + 0.2
  
  if(length(w_vec)==1 | length(b_vec)==1) {
    JaccardIndex=0
    JaccardIndex_op=0
  } else {
    # generate kernel densities
    dw <- density(w_vec, from=lower, to=upper)
    db <- density(b_vec, from=lower, to=upper)
    d <- data.frame(x=dw$x, a=dw$y, b=db$y)
    if (plot) {
      plot(dw, main="overlaid density plots")
      lines(db, col="blue")
      legend("topleft", legend=c(deparse(substitute(w_vec)), deparse(substitute(b_vec))),
             col=c("black", "blue"), lty=1:1, cex=0.5)
    }
    
    # calculate intersection densities
    d$w <- pmin(d$a, d$b)
    
    # integrate areas under curves
    total <- integrate.xy(d$x, d$a) + integrate.xy(d$x, d$b)
    intersection <- integrate.xy(d$x, d$w)
    
    # compute overlap coefficient
    JaccardIndex <- intersection/(total-intersection)
    JaccardIndex_op <- ifelse(dw$x[which.max(dw$y)]<=db$x[which.max(db$y)], JaccardIndex, -1*JaccardIndex)
  }  
  return(list(JaccardIndex=JaccardIndex, signJaccardIndex=JaccardIndex_op))
}



tailsArea_fun <- function(w_vec, b_vec, tails=0.1) {
  
  if (!require("sfsmisc")) install.packages("sfsmisc")
  library(sfsmisc)

  # integrate areas between -0.1 and 0.1
  dw <- density(w_vec, from=-1*tails, to=tails)
  db <- density(b_vec, from=-1*tails, to=tails)
  d <- data.frame(x=dw$x, a=dw$y, b=db$y)
  
  if (min(w_vec)>-1*tails & max(w_vec) <tails) area_dw <- 0 else area_dw <- 1-integrate.xy(d$x, d$a)
  if (min(b_vec)>-1*tails & max(b_vec) <tails) area_db <- 0 else area_db <- 1-integrate.xy(d$x, d$b)

  # integrate areas under curves: 
  lower <- min(c(w_vec, b_vec)) - 0.2
  upper <- max(c(w_vec, b_vec)) + 0.2
  
  dw <- density(w_vec, from=lower, to=upper)
  db <- density(b_vec, from=lower, to=upper)
  d <- data.frame(x=dw$x, a=dw$y, b=db$y)
  
  allarea_dw <- integrate.xy(d$x, d$a)
  allarea_db <- integrate.xy(d$x, d$b)
  
  return(c(area_dw/allarea_dw, area_db/allarea_db))
}

