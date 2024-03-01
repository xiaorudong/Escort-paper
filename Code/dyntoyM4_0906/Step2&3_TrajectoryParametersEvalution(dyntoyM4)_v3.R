
rm(list = ls())

setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")
source("CODE/MethodScripts/Funtions/FunctionsFORStep2&3_v12.R")
# dir.create("PLOTS/Script_dyntoy_M4")


#################################
######### Load Data Set #########
#################################

load("RDATA/step0_clean_dyntoy_M4.RData")
load("RDATA/step0_ss_dyntoy_M4.RData")
load("RDATA/step1_dyntoy_M4.RData")



######################## STEP 2:########################


###########################################################
######## Check Disconnected Clusters in Embedding #########
###########################################################


library(scales)
par(mfrow=c(2,2))
DRLvsC1 <- LD_DCClusterscheck(dist_mat=dist(MDS1_eval$Embedding, method = "euclidean"), DRdims=MDS1_eval$Embedding, connectedCells = 1)
plot(MDS1_eval$Embedding, col=alpha(DRLvsC1$Clusters,0.5), pch=16, main=deparse(substitute(MDS1_eval)))
legend("bottomright", legend=as.character(1:DRLvsC1$K), col=c(1:DRLvsC1$K), pch=16, cex = 0.5)
DRLvsC1$DCcheck
DRLvsC1$clusterLocation


DRLvsC2 <- LD_DCClusterscheck(dist_mat=dist(MDS2_eval$Embedding, method = "euclidean"), DRdims=MDS2_eval$Embedding, connectedCells = 1)
plot(MDS2_eval$Embedding, col=alpha(DRLvsC2$Clusters,0.5), pch=16, main=deparse(substitute(MDS2_eval)))
legend("topright", legend=as.character(1:DRLvsC2$K), col=c(1:DRLvsC2$K), pch=16, cex = 0.5)
DRLvsC2$DCcheck
DRLvsC2$clusterLocation

DRLvsC3 <- LD_DCClusterscheck(dist_mat=dist(MDS3_eval$Embedding, method = "euclidean"), DRdims=MDS3_eval$Embedding, connectedCells = 1)
plot(MDS3_eval$Embedding, col=alpha(DRLvsC3$Clusters,0.5), pch=16, main=deparse(substitute(MDS3_eval)))
legend("topleft", legend=as.character(1:DRLvsC3$K), col=c(1:DRLvsC3$K), pch=16, cex = 0.5)
DRLvsC3$DCcheck
DRLvsC3$clusterLocation

DRLvsC4 <- LD_DCClusterscheck(dist_mat=dist(TSNE4_eval$Embedding, method = "euclidean"), DRdims=TSNE4_eval$Embedding, connectedCells = 1)
plot(TSNE4_eval$Embedding, col=alpha(DRLvsC4$Clusters,0.5), pch=16, main=deparse(substitute(TSNE4_eval)))
legend("topleft", legend=as.character(1:DRLvsC4$K), col=c(1:DRLvsC4$K), pch=16, cex = 0.5)
DRLvsC4$DCcheck
DRLvsC4$clusterLocation

DRLvsC5 <- LD_DCClusterscheck(dist_mat=dist(TSNE5_eval$Embedding, method = "euclidean"), DRdims=TSNE5_eval$Embedding, connectedCells = 1)
plot(TSNE5_eval$Embedding, col=alpha(DRLvsC5$Clusters,0.5), pch=16, main=deparse(substitute(TSNE5_eval)))
legend("bottomleft", legend=as.character(1:DRLvsC5$K), col=c(1:DRLvsC5$K), pch=16, cex = 0.5)
DRLvsC5$DCcheck
DRLvsC5$clusterLocation

DRLvsC6 <- LD_DCClusterscheck(dist_mat=dist(TSNE6_eval$Embedding, method = "euclidean"), DRdims=TSNE6_eval$Embedding, connectedCells = 1)
plot(TSNE6_eval$Embedding, col=alpha(DRLvsC6$Clusters,0.5), pch=16, main=deparse(substitute(TSNE6_eval)))
legend("bottomleft", legend=as.character(1:DRLvsC6$K), col=c(1:DRLvsC6$K), pch=16, cex = 0.5)
DRLvsC6$DCcheck
DRLvsC6$clusterLocation

DRLvsC7 <- LD_DCClusterscheck(dist_mat=dist(UMAP7_eval$Embedding, method = "euclidean"), DRdims=UMAP7_eval$Embedding, connectedCells = 1)
plot(UMAP7_eval$Embedding, col=alpha(DRLvsC7$Clusters,0.5), pch=16, main=deparse(substitute(UMAP7_eval)))
legend("topleft", legend=as.character(1:DRLvsC7$K), col=c(1:DRLvsC7$K), pch=16, cex = 0.5)
DRLvsC7$DCcheck
DRLvsC7$clusterLocation

DRLvsC8 <- LD_DCClusterscheck(dist_mat=dist(UMAP8_eval$Embedding, method = "euclidean"), DRdims=UMAP8_eval$Embedding, connectedCells = 1)
plot(UMAP8_eval$Embedding, col=alpha(DRLvsC8$Clusters,0.5), pch=16, main=deparse(substitute(UMAP8_eval)))
legend("topleft", legend=as.character(1:DRLvsC8$K), col=c(1:DRLvsC8$K), pch=16, cex = 0.5)
DRLvsC8$DCcheck
DRLvsC8$clusterLocation

DRLvsC9 <- LD_DCClusterscheck(dist_mat=dist(UMAP9_eval$Embedding, method = "euclidean"), DRdims=UMAP9_eval$Embedding, connectedCells = 1)
plot(UMAP9_eval$Embedding, col=alpha(DRLvsC9$Clusters,0.5), pch=16, main=deparse(substitute(UMAP9_eval)))
legend("topleft", legend=as.character(1:DRLvsC9$K), col=c(1:DRLvsC9$K), pch=16, cex = 0.5)
DRLvsC9$DCcheck
DRLvsC9$clusterLocation



####################################################
######## Check the similarity within cells #########
####################################################

par(mfrow=c(2,2))
simiRetain <- c()

simi_cells <- Similaritycheck(norm_counts=norm_counts, obj=MDS1_eval, Cluters=LvsC) # we use all norm counts instead of counts data after feature selection
simi_cells$SimilarityLevel
(rate_good <- sum(simi_cells$SimilarityLevel[as.numeric(names(simi_cells$SimilarityLevel))<=1])/ncol(norm_counts))
simiRetain <- c(simiRetain, rate_good)
simi_cells1 <- simi_cells

simi_cells <- Similaritycheck(norm_counts=norm_counts, obj=MDS2_eval, Cluters=LvsC) # we use all norm counts instead of counts data after feature selection
simi_cells$SimilarityLevel
(rate_good <- sum(simi_cells$SimilarityLevel[as.numeric(names(simi_cells$SimilarityLevel))<=1])/ncol(norm_counts))
simiRetain <- c(simiRetain, rate_good)
simi_cells2 <- simi_cells

simi_cells <- Similaritycheck(norm_counts=norm_counts, obj=MDS3_eval, Cluters=LvsC) # we use all norm counts instead of counts data after feature selection
simi_cells$SimilarityLevel
(rate_good <- sum(simi_cells$SimilarityLevel[as.numeric(names(simi_cells$SimilarityLevel))<=1])/ncol(norm_counts))
simiRetain <- c(simiRetain, rate_good)
simi_cells3 <- simi_cells

simi_cells <- Similaritycheck(norm_counts=norm_counts, obj=TSNE4_eval, Cluters=LvsC) # we use all norm counts instead of counts data after feature selection
simi_cells$SimilarityLevel
(rate_good <- sum(simi_cells$SimilarityLevel[as.numeric(names(simi_cells$SimilarityLevel))<=1])/ncol(norm_counts))
simiRetain <- c(simiRetain, rate_good)
simi_cells4 <- simi_cells

simi_cells <- Similaritycheck(norm_counts=norm_counts, obj=TSNE5_eval, Cluters=LvsC) # we use all norm counts instead of counts data after feature selection
simi_cells$SimilarityLevel
(rate_good <- sum(simi_cells$SimilarityLevel[as.numeric(names(simi_cells$SimilarityLevel))<=1])/ncol(norm_counts))
simiRetain <- c(simiRetain, rate_good)
simi_cells5 <- simi_cells

simi_cells <- Similaritycheck(norm_counts=norm_counts, obj=TSNE6_eval, Cluters=LvsC) # we use all norm counts instead of counts data after feature selection
simi_cells$SimilarityLevel
(rate_good <- sum(simi_cells$SimilarityLevel[as.numeric(names(simi_cells$SimilarityLevel))<=1])/ncol(norm_counts))
simiRetain <- c(simiRetain, rate_good)
simi_cells6 <- simi_cells

simi_cells <- Similaritycheck(norm_counts=norm_counts, obj=UMAP7_eval, Cluters=LvsC) # we use all norm counts instead of counts data after feature selection
simi_cells$SimilarityLevel
(rate_good <- sum(simi_cells$SimilarityLevel[as.numeric(names(simi_cells$SimilarityLevel))<=1])/ncol(norm_counts))
simiRetain <- c(simiRetain, rate_good)
simi_cells7 <- simi_cells

simi_cells <- Similaritycheck(norm_counts=norm_counts, obj=UMAP8_eval, Cluters=LvsC) # we use all norm counts instead of counts data after feature selection
simi_cells$SimilarityLevel
(rate_good <- sum(simi_cells$SimilarityLevel[as.numeric(names(simi_cells$SimilarityLevel))<=1])/ncol(norm_counts))
simiRetain <- c(simiRetain, rate_good)
simi_cells8 <- simi_cells

simi_cells <- Similaritycheck(norm_counts=norm_counts, obj=UMAP9_eval, Cluters=LvsC) # we use all norm counts instead of counts data after feature selection
simi_cells$SimilarityLevel
(rate_good <- sum(simi_cells$SimilarityLevel[as.numeric(names(simi_cells$SimilarityLevel))<=1])/ncol(norm_counts))
simiRetain <- c(simiRetain, rate_good)
simi_cells9 <- simi_cells

simiRetain
simiRetain >=0.8


####################################################################
######## Filter Methods cannot reflect cells' relationship #########
####################################################################

DRLvsC1$DCcheck
DRLvsC2$DCcheck
DRLvsC3$DCcheck
DRLvsC4$DCcheck
DRLvsC5$DCcheck
DRLvsC6$DCcheck
DRLvsC7$DCcheck
DRLvsC8$DCcheck
DRLvsC9$DCcheck 
ifConn <- c(DRLvsC1$ifConnected, DRLvsC2$ifConnected, DRLvsC3$ifConnected,
            DRLvsC4$ifConnected, DRLvsC5$ifConnected, DRLvsC6$ifConnected,
            DRLvsC7$ifConnected, DRLvsC8$ifConnected, DRLvsC9$ifConnected)

simiRetain
simiRetain >=0.8 


scoredf <- data.frame(DCcheck=ifConn, SimiRetain=simiRetain)
rownames(scoredf) <- c("MDS1_eval", "MDS2_eval", "MDS3_eval",
                       "TSNE4_eval", "TSNE5_eval", "TSNE6_eval",
                       "UMAP7_eval", "UMAP8_eval", "UMAP9_eval")  

scoredf



########################################
######## Check goodness of fit #########
########################################

pdf("PLOTS/Script_dyntoy_M4/gof_ss_evaluation.pdf", onefile = TRUE, width = 8, height = 4.5)
par(mfrow=c(1, 2))
gof_eval1 <- GOFeval(MDS1_eval)
gof_eval2 <- GOFeval(MDS2_eval)
gof_eval3 <- GOFeval(MDS3_eval)
gof_eval4 <- GOFeval(TSNE4_eval)
gof_eval5 <- GOFeval(TSNE5_eval)
gof_eval6 <- GOFeval(TSNE6_eval)
gof_eval7 <- GOFeval(UMAP7_eval)
gof_eval8 <- GOFeval(UMAP8_eval)
gof_eval9 <- GOFeval(UMAP9_eval)
dev.off()

scoredf$GOF <- c(gof_eval1$occupiedRate, gof_eval2$occupiedRate, gof_eval3$occupiedRate, 
                 gof_eval4$occupiedRate, gof_eval5$occupiedRate, gof_eval6$occupiedRate,
                 gof_eval7$occupiedRate, gof_eval8$occupiedRate, gof_eval9$occupiedRate) 




######################## STEP 3:########################


#######################################################
######## Check the shape of trajectory curves #########
#######################################################

source("CODE/MethodScripts/Funtions/updated_ushape_cal.R")
pdf("PLOTS/Script_dyntoy_M4/ushape_ss_evaluation.pdf", onefile = TRUE, width = 8, height = 4.5)
par(mfrow=c(1, 2))
ushap_eval1 <- UshapeDetector(MDS1_eval)
ushap_eval2 <- UshapeDetector(MDS2_eval)
ushap_eval3 <- UshapeDetector(MDS3_eval)
ushap_eval4 <- UshapeDetector(TSNE4_eval)
ushap_eval5 <- UshapeDetector(TSNE5_eval)
ushap_eval6 <- UshapeDetector(TSNE6_eval)
ushap_eval7 <- UshapeDetector(UMAP7_eval)
ushap_eval8 <- UshapeDetector(UMAP8_eval)
ushap_eval9 <- UshapeDetector(UMAP9_eval)
dev.off()

goodmed <-  c(ushap_eval1$NoACells, ushap_eval2$NoACells, ushap_eval3$NoACells,
              ushap_eval4$NoACells, ushap_eval5$NoACells, ushap_eval6$NoACells,
              ushap_eval7$NoACells, ushap_eval8$NoACells, ushap_eval9$NoACells) 
scoredf$USHAPE <- goodmed



###################################
######## Calculate Scores #########
###################################

final_df <- score_cal(scoredf)

save(scoredf, final_df, file = "RDATA/res_dyntoy_M4.RData")
save.image(file = "RDATA/allres_dyntoy_M4.RData")










