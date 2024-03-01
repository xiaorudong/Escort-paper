rm(list = ls())

library(scran)
library(dplyr)

# https://github.com/dynverse/dynbenchmark
# https://zenodo.org/record/1443566#.YYMvcC-B19c
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")

source("CODE/MethodScripts/Funtions/FunctionsFORStep0.R")


##############################
######### Clean Data #########
##############################

library(Seurat)
scdata <- readRDS("DATA/sc_benchmark_real/aging-hsc-young_kowalczyk.rds")
scdata_time <- scdata[["prior_information"]][["timecourse_continuous"]]

library(scran)
library(org.Mm.eg.db)
sce <- SingleCellExperiment(assays = list(counts = t(scdata[["counts"]])))
filter.sce <- quickPerCellQC(sce)
filter.sce <- logNormCounts(filter.sce)

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
anno <- select(org.Mm.eg.db, keys = rownames(filter.sce), keytype = "SYMBOL", column = "ENSEMBL")
ensembl <- anno$ENSEMBL[match(rownames(filter.sce), anno$SYMBOL)]
rowData(filter.sce)$ENSEMBL <- ensembl
assignments <- cyclone(filter.sce, mm.pairs, gene.names=rowData(filter.sce)$ENSEMBL)
# design <- model.matrix(~as.matrix(assignments$scores))
# library(batchelor)
# reg.nocycle <- regressBatches(filter.sce, design=design)

filter.sce <- runPCA(filter.sce)
plotPCA(filter.sce, colour_by=I(assignments$phases), point_size=3)


library(scater)
diff <- getVarianceExplained(filter.sce, DataFrame(assignments$phases))
discard <- diff > 5
summary(discard)

top.discard <- getTopHVGs(filter.sce[which(!discard),], n=1000)
filter.sce.discard <- runPCA(filter.sce, subset_row=top.discard)

plotPCA(filter.sce.discard, colour_by=I(assignments$phases), point_size=3)







rawcounts <- t(scdata[["counts"]])
scdata <- CreateSeuratObject(counts = rawcounts)
scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT-")
VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

scdata <- subset(scdata, subset = nCount_RNA < 5*10^6)
selected_f <- rownames(scdata)[Matrix::rowSums(scdata) > 3]
scdata <- subset(scdata, features = selected_f)
scdata <- NormalizeData(scdata)
scdata <- ScaleData(scdata, features = rownames(scdata))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scdata <- CellCycleScoring(scdata, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
scdata <- ScaleData(scdata, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(scdata))

correct_counts <- as.matrix(scdata@assays$RNA$scale.data)
norm_counts <- as.matrix(scdata@assays$RNA$data)
rawcounts <- as.matrix(scdata@assays$RNA$counts)

a <- match(colnames(norm_counts), names(scdata_time))
scdata_time <- scdata_time[a]

save(correct_counts, norm_counts, rawcounts, scdata_time,
     file = "RDATA/step0_clean_removeCC.RData")


###################################
######### Select Features #########
###################################

# select genes
library(scran)
gene.var <- modelGeneVar(norm_counts)
genes.HVGs_1 <- getTopHVGs(gene.var, prop=0.2, var.threshold=NULL)
genes.HVGs_2 <- getTopHVGs(gene.var, prop=0.4, var.threshold=NULL)
genes.HVGs_3 <- getTopHVGs(gene.var, prop=1, var.threshold=NULL)


sub_counts1 <- correct_counts[genes.HVGs_1,]
sub_counts2 <- correct_counts[genes.HVGs_2,]
sub_counts3 <- correct_counts[genes.HVGs_3,]


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
     file = "RDATA/step0_dr_removeCC.RData")
