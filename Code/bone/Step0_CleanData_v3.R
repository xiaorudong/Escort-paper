rm(list = ls())

library(scran)
library(Seurat)
library(dplyr)

# https://github.com/dynverse/dynbenchmark
# https://zenodo.org/record/1443566#.YYMvcC-B19c
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")

source("CODE/MethodScripts/Funtions/FunctionsFORStep0.R")
load("RDATA/res_bone_SSPC_afterQC.RData")

norm_counts <- as.matrix(subdata@assays[["RNA"]]@data)
rawcounts <- as.matrix(subdata@assays[["RNA"]]@counts)


save(norm_counts, rawcounts,
     file = "RDATA/example_realData_bone/step0_clean_afterQC.RData")


###################################
######### Select Features #########
###################################

# select genes
library(scuttle)
library(scran)
gene.var <- modelGeneVar(norm_counts)
genes.HVGs_1 <- getTopHVGs(gene.var, prop=0.1, var.threshold=NULL)
genes.HVGs_2 <- getTopHVGs(gene.var, prop=0.2, var.threshold=NULL)
genes.HVGs_3 <- getTopHVGs(gene.var, prop=0.5, var.threshold=NULL)


sub_counts1 <- norm_counts[genes.HVGs_1,]
sub_counts2 <- norm_counts[genes.HVGs_2,]
sub_counts3 <- norm_counts[genes.HVGs_3,]


#######################################
######### Dimension Reduction #########
#######################################

set.seed(111)
# get the DR
dimred1 <- DR_2D(sub_counts1, "MDS")
dimred2 <- DR_2D(sub_counts2, "MDS")
dimred3 <- DR_2D(sub_counts3, "MDS")
dimred4 <- DR_2D(sub_counts1, "TSNE")
dimred5 <- DR_2D(sub_counts2, "TSNE")
dimred6 <- DR_2D(sub_counts3, "TSNE")
dimred7 <- DR_2D(sub_counts1, "UMAP")
dimred8 <- DR_2D(sub_counts2, "UMAP")
dimred9 <- DR_2D(sub_counts3, "UMAP")
dimred10 <- DR_2D(sub_counts1, "PCA")
dimred11 <- DR_2D(sub_counts2, "PCA")
dimred12 <- DR_2D(sub_counts3, "PCA")


save(dimred1, dimred2, dimred3,
     dimred4, dimred5, dimred6,
     dimred7, dimred8, dimred9,
     dimred10, dimred11, dimred12,
     file = "RDATA/example_realData_bone/step0_dr_afterQC.RData")
