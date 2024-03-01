# forRuby:
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORStep0.R")

library(scRNAseq)
library(scaffold)
library(scater)
library(dplyr)

# ref_data <- scRNAseq::ZilionisLungData(which = "human", filter = TRUE)
# ref_data <- ref_data[, stringr::str_detect(ref_data$`Major cell type`, "tNeutrophils")]
# ref_data <- ref_data[which(rowSums(counts(ref_data) > 0) >= 5), ]


# sce <- scRNAseq::LaMannoBrainData(which = "human-embryo")
# sce <- sce[rowSums(SingleCellExperiment::counts(sce) > 0) >= 3, colSums(SingleCellExperiment::counts(sce)) > 0]
# sce <- sce[Rfast::Sort(rownames(sce)), ]  # necessary to play nice with scaffold


sce <- scRNAseq::BaronPancreasData(which = "human")
sce <- sce[rowSums(SingleCellExperiment::counts(sce) > 0) >= 3, colSums(SingleCellExperiment::counts(sce)) > 0]
sce <- sce[Rfast::Sort(rownames(sce)), ]

# #
# sce <- logNormCounts(sce)
# norm_counts <- as.matrix(sce@assays@data@listData[["logcounts"]])
# rawcounts <- as.matrix(sce@assays@data@listData[["counts"]])
# 
# dimred1 <- DR_2D(norm_counts, "MDS")
# dimred2 <- DR_2D(norm_counts, "UMAP")
# dimred3 <- DR_2D(norm_counts, "TSNE")
# 
# dimred <- dimred1
# plot(dimred, col = as.factor(sce@colData@listData[["label"]]), pch=16)
# 
# dimred <- dimred2
# plot(dimred, col = as.factor(sce@colData@listData[["label"]]), pch=16)
# 
# dimred <- dimred3
# plot(dimred, col = as.factor(sce@colData@listData[["label"]]), pch=16)
# 
# 
# library(grDevices)
# library(RColorBrewer)
# dimred <- dimred1
# clust <- kmeans(dimred, centers = 1, nstart = kmeans_nstart)
# cl1 <- clust$cluster
# ti_out <- slingshot(data=dimred, clusterLabels=cl1)
# rawpse <- slingPseudotime(ti_out, na=T)
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(rawpse, breaks=100)]
# plot(dimred, col = plotcol, pch=16, asp = 1)
# 
# dimred <- dimred2
# clust <- kmeans(dimred, centers = 1, nstart = kmeans_nstart)
# cl1 <- clust$cluster
# ti_out <- slingshot(data=dimred, clusterLabels=cl1)
# rawpse <- slingPseudotime(ti_out, na=T)
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(rawpse, breaks=100)]
# plot(dimred, col = plotcol, pch=16, asp = 1)
# 
# dimred <- dimred3
# clust <- kmeans(dimred, centers = 1, nstart = kmeans_nstart)
# cl1 <- clust$cluster
# ti_out <- slingshot(data=dimred, clusterLabels=cl1)
# rawpse <- slingPseudotime(ti_out, na=T)
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(rawpse, breaks=100)]
# plot(dimred, col = plotcol, pch=16, asp = 1)


ref.dataset <- sce

perc.dyn.genes = 0.001
n.cells = 500

# set up simulation parameters
set.seed(1)
n_dyn_genes <- ceiling(perc.dyn.genes * nrow(ref.dataset))
my_knots <- matrix(runif(2 * n_dyn_genes, 0, 1), ncol = 2, nrow = n_dyn_genes)
my_theta <- matrix(rnorm(5, 5, 5), ncol = 5, nrow = n_dyn_genes)
dynamic_params <- list(propGenes = perc.dyn.genes,
                       degree = 2,
                       knots = my_knots,
                       theta = my_theta)

# simulate 10X dataset
boxplot(colSums(counts(ref.dataset)))
ref.dataset@metadata$seqdepth <- colSums(counts(ref.dataset))
ref.dataset <- ref.dataset[,colSums(counts(ref.dataset)) >3000 & colSums(counts(ref.dataset)) < 8000]
dim(ref.dataset)
scaffold_params <- estimateScaffoldParameters(sce = ref.dataset,
                                              sceUMI = TRUE,
                                              useUMI = TRUE,
                                              protocol = "droplet",
                                              numCells = n.cells,
                                              popHet = c(1, 1))
sim_data <- simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.dataset)
sim_data <- sim_data[rowSums(counts(sim_data)) > 0, ]  # only non-zero genes

# typical scran + scater pre-processing pipeline
colData(sim_data) <- colData(sim_data) %>%
  as.data.frame() %>%
  dplyr::mutate(cell_time = as.numeric(gsub("Cell_", "", rownames(.))),
                cell_time_normed = cell_time / max(cell_time)) %>%
  S4Vectors::DataFrame()


saveRDS(sim_data, file = "DATA/simScaffold_forRuby_sim_DEG0_500.rds")


