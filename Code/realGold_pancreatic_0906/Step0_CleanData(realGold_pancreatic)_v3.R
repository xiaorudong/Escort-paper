rm(list = ls())

library(scran)
library(Seurat)
library(dplyr)

# https://github.com/dynverse/dynbenchmark
# https://zenodo.org/record/1443566#.YYMvcC-B19c
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")

source("CODE/MethodScripts/Funtions/FunctionsFORStep0.R")


##############################
######### Clean Data #########
##############################

scdata <- readRDS("DATA/sc_benchmark_real/pancreatic-alpha-cell-maturation_zhang.rds")
scdata_time <- scdata[["prior_information"]][["timecourse_continuous"]]

rawcounts <- t(scdata[["counts"]])
scdata <- CreateSeuratObject(counts = rawcounts)
scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT-")
VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

scdata <- subset(scdata, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000)
selected_f <- rownames(scdata)[Matrix::rowSums(scdata) > 3]
scdata <- subset(scdata, features = selected_f)
scdata <- NormalizeData(scdata)

norm_counts <- as.matrix(scdata@assays[["RNA"]]@data)
rawcounts <- as.matrix(scdata@assays[["RNA"]]@counts)

a <- match(colnames(norm_counts), names(scdata_time))
scdata_time <- scdata_time[a]

save(norm_counts, rawcounts, scdata_time, 
     file = "RDATA/step0_clean_realGold_pancreatic.RData")


###################################
######### Select Features #########
###################################

# select genes
library(scuttle)
library(scran)
gene.var <- modelGeneVar(norm_counts)
genes.HVGs_1 <- getTopHVGs(gene.var, prop=0.2, var.threshold=NULL)
genes.HVGs_2 <- getTopHVGs(gene.var, prop=0.4, var.threshold=NULL)
genes.HVGs_3 <- getTopHVGs(gene.var, prop=1, var.threshold=NULL)


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


save(dimred1, dimred2, dimred3,
     dimred4, dimred5, dimred6,
     dimred7, dimred8, dimred9,
     file = "RDATA/step0_dr_realGold_pancreatic.RData")
