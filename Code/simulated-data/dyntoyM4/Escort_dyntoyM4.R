rm(list = ls())

library(Escort)

setwd("path/to/your_file")



##############################
######### Clean Data #########
##############################

scdata <- readRDS("DATA/sc_benchmark/synthetic_dyntoy_multifurcating_4.rds")
scdata_time <- scdata[["prior_information"]][["timecourse_continuous"]]
rawcounts <- t(scdata[["counts"]])

# remove cells with duplicated time
scdata_time <- scdata_time[!duplicated(scdata_time)]
rawcounts <- rawcounts[,colnames(rawcounts) %in% names(scdata_time)]

library(Seurat)
scdata <- CreateSeuratObject(counts = rawcounts)
scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = "^MT-")
VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

selected_f <- rownames(scdata)[Matrix::rowSums(scdata) > 3]
scdata <- subset(scdata, features = selected_f)
scdata <- NormalizeData(scdata)

norm_counts <- as.matrix(scdata@assays$RNA$data)
rawcounts <- as.matrix(scdata@assays$RNA$counts)

a <- match(colnames(norm_counts), names(scdata_time))
scdata_time <- scdata_time[a]

save(norm_counts, rawcounts, scdata_time, 
     file = "RDATA/Manuscript/clean_dyntoy_M4.RData")



##########################
######### Step 1 #########
##########################

# Disconnected Clusters vs Lineages
LvsC <- HD_DCClusterscheck(normcounts=norm_counts, rawcounts=rawcounts)
LvsC$DCcheck

# Test Homogeneous
cor_test <- step1_testHomogeneous(normcounts=norm_counts)
cor_test$decision
cor_test$signal_pct



###################################
######### Select Features #########
###################################

# select genes
gene.var <- quick_model_gene_var(norm_counts)
genes.HVGs_1 <- rownames(gene.var)[1:round(nrow(gene.var)*0.2)]
genes.HVGs_2 <- rownames(gene.var)[1:round(nrow(gene.var)*0.4)]
genes.HVGs_3 <- rownames(gene.var)

sub_counts1 <- norm_counts[genes.HVGs_1,]
sub_counts2 <- norm_counts[genes.HVGs_2,]
sub_counts3 <- norm_counts[genes.HVGs_3,]


#######################################
######### Dimension Reduction #########
#######################################

set.seed(111)
# get the DR
dimred1 <- getDR_2D(sub_counts1, "MDS")
dimred2 <- getDR_2D(sub_counts2, "MDS")
dimred3 <- getDR_2D(sub_counts3, "MDS")
dimred4 <- getDR_2D(sub_counts1, "TSNE")
dimred5 <- getDR_2D(sub_counts2, "TSNE")
dimred6 <- getDR_2D(sub_counts3, "TSNE")
dimred7 <- getDR_2D(sub_counts1, "UMAP")
dimred8 <- getDR_2D(sub_counts2, "UMAP")
dimred9 <- getDR_2D(sub_counts3, "UMAP")

embeddings <- list(dimred1, dimred2, dimred3,
                   dimred4, dimred5, dimred6,
                   dimred7, dimred8, dimred9)



##########################
######### Step 2 #########
##########################

# Check Disconnected Clusters in Embedding
library(parallel)
DRLvsCs <- mclapply(embeddings, LD_DCClusterscheck)
lapply(DRLvsCs, function(x) x$DCcheck)


# Check the similarity within cells
simi_cells <- mclapply(embeddings, function(x) {
  Similaritycheck(normcounts=norm_counts, dimred=x, Cluters=LvsC)
})
lapply(simi_cells, function(x) x$GoodRate)

# Check goodness of fit
gof_evals <- mclapply(embeddings, GOFeval)
lapply(gof_evals, function(x) x$occupiedRate)

save(embeddings, file = "RDATA/Manuscript/dr_dyntoy_M4.RData")



##########################
######### Step 3 #########
##########################

# from Slingshot:
load("RDATA/Manuscript/slingshot_dyntoyM4.RData")
ushap_ss_evals <- mclapply(resobjs, UshapeDetector)
lapply(ushap_ss_evals, function(x) x$Ambpct)


# from Monocle3:
load("RDATA/Manuscript/Monocle3_dyntoyM4.RData")
ushap_m3_evals <- mclapply(resobjs, function(x) UshapeDetector(x, outlierdetect = "resistant"))
lapply(ushap_m3_evals, function(x) x$Ambpct)



##################################
######### Scoring System #########
##################################

scoredf_ss <- data.frame(DCcheck=sapply(DRLvsCs, function(x) x$ifConnected),
                      SimiRetain=sapply(simi_cells, function(x) x$GoodRate),
                      GOF=sapply(gof_evals, function(x) x$occupiedRate), 
                      USHAPE=sapply(ushap_ss_evals, function(x) x$Ambpct))
rownames(scoredf_ss) <- c("MDS; 20% genes", "MDS; 40% genes", "MDS; all genes",
                       "TSNE; 20% genes", "TSNE; 40% genes", "TSNE; all genes",
                       "UMAP; 20% genes", "UMAP; 40% genes", "UMAP; all genes")
calcScore(scoredf_ss)



scoredf_m3 <- data.frame(DCcheck=sapply(DRLvsCs, function(x) x$ifConnected),
                      SimiRetain=sapply(simi_cells, function(x) x$GoodRate),
                      GOF=sapply(gof_evals, function(x) x$occupiedRate), 
                      USHAPE=sapply(ushap_m3_evals, function(x) x$Ambpct))
rownames(scoredf_m3) <- c("MDS; 20% genes", "MDS; 40% genes", "MDS; all genes",
                       "TSNE; 20% genes", "TSNE; 40% genes", "TSNE; all genes",
                       "UMAP; 20% genes", "UMAP; 40% genes", "UMAP; all genes")
calcScore(scoredf_m3)





