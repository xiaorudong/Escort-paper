rm(list = ls())

library(scran)
library(Seurat)
library(dplyr)
library(mclust)

# https://github.com/dynverse/dynbenchmark
# https://zenodo.org/record/1443566#.YYMvcC-B19c
setwd("/Volumes/bachergroup/TrajectoryEvaluation_Method/")

source("CODE/MethodScripts/Funtions/FunctionsFORStep0.R")


load("RDATA/example_realData_bone/step0_clean_afterQC.RData")
load("RDATA/example_realData_bone/step0_dr_afterQC.RData")
load("RDATA/res_bone_SSPC_afterQC.RData")



#######################################
######### Trajectory Analysis #########
#######################################

# Take using Slingshot to do trajectory analysis as an example
library(slingshot)
library(DescTools)
library(vioplot)
library(scales)
library(mclust)
set.seed(123)
# kmeans_nstart= ifelse(ncol(norm_counts)>2000, 50, 1000)

# for subset 1
dimred <- dimred1
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out1 <- ti_out
rawpse1 <- rawpse
pse1 <- pse
fitLine1 <- as.data.frame(fitLine)


# for subset 2
dimred <- dimred2
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out2 <- ti_out
rawpse2 <- rawpse
pse2 <- pse
fitLine2 <- as.data.frame(fitLine)

# for subset 3
dimred <- dimred3
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out3 <- ti_out
rawpse3 <- rawpse
pse3 <- pse
fitLine3 <- as.data.frame(fitLine)

# for subset 4
dimred <- dimred4
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out4 <- ti_out
rawpse4 <- rawpse
pse4 <- pse
fitLine4 <- as.data.frame(fitLine)

# for subset 5
dimred <- dimred5
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out5 <- ti_out
rawpse5 <- rawpse
pse5 <- pse
fitLine5 <- as.data.frame(fitLine)

# for subset 6
dimred <- dimred6
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out6 <- ti_out
rawpse6 <- rawpse
pse6 <- pse
fitLine6 <- as.data.frame(fitLine)

# for subset 7
dimred <- dimred7
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out7 <- ti_out
rawpse7 <- rawpse
pse7 <- pse
fitLine7 <- as.data.frame(fitLine)

# for subset 8
dimred <- dimred8
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out8 <- ti_out
rawpse8 <- rawpse
pse8 <- pse
fitLine8 <- as.data.frame(fitLine)

# for subset 9
dimred <- dimred9
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out9 <- ti_out
rawpse9 <- rawpse
pse9 <- pse
fitLine9 <- as.data.frame(fitLine)


# for subset 10
dimred <- dimred10
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out10 <- ti_out
rawpse10 <- rawpse
pse10 <- pse
fitLine10 <- as.data.frame(fitLine)

# for subset 11
dimred <- dimred11
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out11 <- ti_out
rawpse11 <- rawpse
pse11 <- pse
fitLine11 <- as.data.frame(fitLine)

# for subset 12
dimred <- dimred12
cl1 <- list.cluster
ti_out <- slingshot(data=dimred, clusterLabels=cl1, start.clus = "Cluster 1")
# ti_out <- slingshot(data=dimred, clusterLabels=rep(1, nrow(dimred)))
rawpse <- slingPseudotime(ti_out, na=T)
pse <- as.data.frame(rawpse / max(rawpse, na.rm=TRUE))
ls_fitLine <- lapply(slingCurves(ti_out), function(x) x$s[x$ord,])
fitLine <- do.call(rbind, lapply(ls_fitLine, function(x) {
  df_seg <- cbind(x[-nrow(x),],x[-1,])
  colnames(df_seg) <- c("x0", "y0", "x1", "y1")
  return(df_seg)
}))
ti_out12 <- ti_out
rawpse12 <- rawpse
pse12 <- pse
fitLine12 <- as.data.frame(fitLine)

#plot:
library(grDevices)
library(RColorBrewer)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
pse$Ave <- rowMeans(pse, na.rm = T)
plotcol <- colors[cut(pse$Ave, breaks=100)]
plot(dimred, col = cl1, pch=16,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
segments(as.data.frame(fitLine)$x0, as.data.frame(fitLine)$y0, 
         as.data.frame(fitLine)$x1, as.data.frame(fitLine)$y1, 
         lwd=2, col='black')


# generate the traj objs

SS_obj1 <- info_traj(dimred1, rawpse = rawpse1, fitLine = fitLine1)
SS_obj2 <- info_traj(dimred2, rawpse = rawpse2, fitLine = fitLine2)
SS_obj3 <- info_traj(dimred3, rawpse = rawpse3, fitLine = fitLine3)
SS_obj4 <- info_traj(dimred4, rawpse = rawpse4, fitLine = fitLine4)
SS_obj5 <- info_traj(dimred5, rawpse = rawpse5, fitLine = fitLine5)
SS_obj6 <- info_traj(dimred6, rawpse = rawpse6, fitLine = fitLine6)
SS_obj7 <- info_traj(dimred7, rawpse = rawpse7, fitLine = fitLine7)
SS_obj8 <- info_traj(dimred8, rawpse = rawpse8, fitLine = fitLine8)
SS_obj9 <- info_traj(dimred9, rawpse = rawpse9, fitLine = fitLine9)
SS_obj10 <- info_traj(dimred10, rawpse = rawpse10, fitLine = fitLine10)
SS_obj11 <- info_traj(dimred11, rawpse = rawpse11, fitLine = fitLine11)
SS_obj12 <- info_traj(dimred12, rawpse = rawpse12, fitLine = fitLine12)

MDS1_eval <- TEvalobj(dimred1, traj = SS_obj1)
MDS2_eval <- TEvalobj(dimred2, traj = SS_obj2)
MDS3_eval <- TEvalobj(dimred3, traj = SS_obj3)
TSNE4_eval <- TEvalobj(dimred4, traj = SS_obj4)
TSNE5_eval <- TEvalobj(dimred5, traj = SS_obj5)
TSNE6_eval <- TEvalobj(dimred6, traj = SS_obj6)
UMAP7_eval <- TEvalobj(dimred7, traj = SS_obj7)
UMAP8_eval <- TEvalobj(dimred8, traj = SS_obj8)
UMAP9_eval <- TEvalobj(dimred9, traj = SS_obj9)
PCA10_eval <- TEvalobj(dimred10, traj = SS_obj10)
PCA11_eval <- TEvalobj(dimred11, traj = SS_obj11)
PCA12_eval <- TEvalobj(dimred12, traj = SS_obj12)



save(MDS1_eval, MDS2_eval, MDS3_eval,
     TSNE4_eval, TSNE5_eval, TSNE6_eval,
     UMAP7_eval, UMAP8_eval, UMAP9_eval,
     PCA10_eval, PCA11_eval, PCA12_eval,
     file = "RDATA/example_realData_bone/step0_ss_cls_afterQC.RData")







