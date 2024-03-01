rm(list = ls())

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORStep1_v18(dc).R")


load("RDATA/step0_clean_dyntoy_L3.RData")


library(scran)
library(slingshot)
library(scales)
library(Seurat)
library(dplyr)

#####################################################
######### Disconnected Clusters vs Lineages #########
#####################################################


library(parallelDist)
par(mfrow = c(2, 2))
# dist_mat <- dist(t(norm_counts), method = 'manhattan')
dist_mat <- parDist(t(norm_counts), method = "manhattan")
LvsC <- HD_DCClusterscheck(dist_mat=dist_mat, rawcounts=rawcounts)

LvsC$DCcheck
# LvsC$Jaccardsummary
LvsC$clusterLocation
K <- LvsC$K



####################################
######### Test Homogeneous #########
####################################

gene.var <- modelGeneVar(norm_counts)
HVGs <- getTopHVGs(gene.var, n=100)
cor_test <- testHomogeneous(HVGs=HVGs, norm_counts=norm_counts)
cor_test$signal_pct
cor_test$decision


######
# Visualization
par(mfrow=c(1,2))
library(umap)
plotcol <- as.factor(LvsC$Clusters)
dimred_umap <- umap::umap(t(norm_counts))$layout
plot(dimred_umap, col = alpha(plotcol,0.7), pch=16, main="UMAP")
legend("topright", legend=as.character(1:K), col=c(1:K), pch=16, cex = 0.5)

library(Rtsne)
dimred_tsne <- Rtsne::Rtsne(t(norm_counts), dims = 2)$Y
rownames(dimred_tsne) <- rownames(t(norm_counts))
plot(dimred_tsne, col = alpha(plotcol,0.7), pch=16, main="TSNE")
legend("topleft", legend=as.character(1:K), col=c(1:K), pch=16, cex = 0.5)



save(LvsC, cor_test, file = "RDATA/step1_dyntoy_L3.RData")





